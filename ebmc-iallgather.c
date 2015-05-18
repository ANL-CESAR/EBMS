/*
 Copyright (c) 2015 UChicago Argonne, LLC

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include "comm.h"

int main(int argc, char **argv)
{
    long  npg;                // global number of particles
    int   r;                 // number of nodes per memory group
    int   nb;                 // num energy bands
    long  npl;                // local number of particles on each proc
    long  gsizekb;            // global size of xs data in kilobytes
    float tracking_rate;      // empirical tracking rate (particles/sec)
    Range *ranges;            // max and min energy for each band
    float **scattering_matrix;// intra-group scattering probability, or size nb x nb
    float *xsdata;     // cross section data (mimicked)
    char *buf1, *buf2; // Double buffering to receive remote xsdata
    char *xscomm; // Point to buf1 or buf2. It is the receive buf for a pending communication
    char *xswork; // Contains xs data we can use for computation. It may point to
                  // buf1, buf2, or a band in xsdata[] if the band is local
    MPI_Request req1 = MPI_REQUEST_NULL;
    MPI_Request req2 = MPI_REQUEST_NULL;
    MPI_Request *reqcomm, *reqwork; // point to req1 or req2
    int flip = 1; // choosing between xscomm/xswork, reqcomm/reqwork
    const float epsilon = 1.0e-8;

    Particle *p;              // local list of particles
    char     sm_file[128];    // path to file containing scattering matrix
    int use_file = FALSE;

    MPI_Datatype dtype_kbytes; // An MPI contiguous datatype of 1024 bytes
    long my_total_alive, global_total_alive, *n_alive;
    double ts, tf, ttot, tave, tmax;
    double ran_val, probability;

    MPI_Comm shmcomm, dsmcomm;
    int my_global_rank, nprocs;
    int my_shm_rank, my_shm_size;
    int my_dsm_rank, my_dsm_size;
    int rc;
    MPI_Aint size;
    int disp_unit;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    srand48(100);

    // Creating a big data type so that the count argument in MPI_Send/Recv/Bcast
    // won't run out of range of an integer, i.e., 2GB
    MPI_Type_contiguous(1024, MPI_CHAR, &dtype_kbytes);
    MPI_Type_commit(&dtype_kbytes);

    // Rank 0 reads user input and then broadcast it
    fflush(stdout);
    if (my_global_rank == 0) {
        process_input(argc, argv, &r, &nb, &gsizekb, &npg,
		              &tracking_rate, &use_file, sm_file);
    }

    MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nb, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gsizekb, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&npg, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tracking_rate, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    assert (nb >= r);

    npl = npg / nprocs;
    if (my_global_rank < npg % nprocs) npl++;

    // init scattering matrix on rank 0 and broadcast
    scattering_matrix = matrix(0,nb-1,0,nb-1);
    if (my_global_rank == 0){
        if (use_file == TRUE){
            read_scattering_matrix(sm_file,scattering_matrix,nb);
        } else {
            create_scattering_matrix(scattering_matrix,nb);
        }
    }
    MPI_Bcast(&scattering_matrix[0][0], nb*nb, MPI_FLOAT, 0, MPI_COMM_WORLD);

    assert ((n_alive = (long *) malloc(nb*sizeof(long))) != NULL);

    // set energy ranges associated with EACH BAND (min..max). This
    // isn't really needed in this communication kernel but is included
    // for conceptual purposes in case we choose to extend kernel
    assert( (ranges =(Range *) malloc(nb*sizeof(Range))) != NULL);
    set_energy_ranges(ranges,nb);

    // allocate local list of particles and initialize all to highest
    // energy band (ie band 0). Exact energies are not important for this
    // kernel, only band.
    assert( (p = (Particle *) malloc(npl*sizeof(Particle))) != NULL);
    init_particles(p, npl, my_global_rank);

    ////////////////////////////////////////////////////////////////////////////
    // Making various communicators out of MPI_COMM_WORLD.
    //
    // * shmcomm: it is a communicator for processes in a node. Processes in
    //     a shmcomm can share memory. There is such a communicator per node.
    // * shm_leader_comm: It is a communicator containing leaders (rank 0's)
    //     of all shmcomm's. There is only one such communicator in the world.
    // * dsmcomm: It is a communicator containing leaders (rank 0's) of
    //      r shmcomm's. There is such a communicator per r nodes.
    ////////////////////////////////////////////////////////////////////////////
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, my_global_rank, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &my_shm_rank);
    MPI_Comm_size(shmcomm, &my_shm_size);

    int is_shm_leader = (my_shm_rank == 0) ? 1 : 0;
    int color = is_shm_leader ? 0 : MPI_UNDEFINED;
    MPI_Comm shm_leader_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, my_global_rank, &shm_leader_comm);

    if (is_shm_leader) {
        int num_shms;
        int shm_idx;
        MPI_Comm_size(shm_leader_comm, &num_shms);

        if (my_global_rank == 0) printf("num_shms = %d\n", num_shms);
        assert(num_shms % r == 0); /* TODO: relax this assumption */

        MPI_Comm_rank(shm_leader_comm, &shm_idx);
        MPI_Comm_split(shm_leader_comm, shm_idx/r, shm_idx, &dsmcomm);
    }

    // Compute dsm rank and size
    int tmp[2];
    if (is_shm_leader) {
        MPI_Comm_rank(dsmcomm, &my_dsm_rank);
        MPI_Comm_size(dsmcomm, &my_dsm_size);
        tmp[0] = my_dsm_rank;
        tmp[1] = my_dsm_size;
    }

    MPI_Bcast(tmp, 2, MPI_INT, 0, shmcomm);
    my_dsm_rank = tmp[0];
    my_dsm_size = tmp[1];

    ////////////////////////////////////////////////////////////////////////////
    //
    // Every node allocates bands it owns in shared memory
    //
    ////////////////////////////////////////////////////////////////////////////

    // For simplicity, assume xsdata can be evenly divided in nb bands
    assert(gsizekb%nb == 0);
    long long bandsizekb = gsizekb/nb;

    // For simplicity, assume a band can be evenly distributed on r nodes
    assert(bandsizekb % r == 0);
    long long slicesizekb = bandsizekb/r; // slice size in kilobytes

    // SHM leader allocates xsdata in its memory, since only leaders do communication.
    // All other processes in the node only touch data in buf1 or buf2 (see below).
    if (is_shm_leader) {
        xsdata = (float*)malloc(slicesizekb*1024*nb);
        assert(xsdata);
    }

    if (is_shm_leader && my_dsm_rank == 0) {
        for (int i = 0; i < nb; i++) {
            int *p = (int *) ((char*)xsdata + slicesizekb*1024*i);
            *p = i;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    // Every node allocates the two buffers in shared memory
    //
    ////////////////////////////////////////////////////////////////////////////
    // For simplicity, assume a band can evenly distributed on cores in a shm.
    // buf1/2 are divided by processes in a shm to even out numa-effect.
    MPI_Win shmwin1, shmwin2, winwork;
    assert(bandsizekb*1024 % my_shm_size == 0);
    rc = MPI_Win_allocate_shared(bandsizekb*1024/my_shm_size, 1, MPI_INFO_NULL, shmcomm, &buf1, &shmwin1); assert(!rc);
    rc = MPI_Win_allocate_shared(bandsizekb*1024/my_shm_size, 1, MPI_INFO_NULL, shmcomm, &buf2, &shmwin2); assert(!rc);

    // Query the start address of the buf on rank 0 of shmcomm
    rc = MPI_Win_shared_query(shmwin1, 0, &size, &disp_unit, &buf1); assert(!rc);
    rc = MPI_Win_shared_query(shmwin2, 0, &size, &disp_unit, &buf2); assert(!rc);

    ////////////////////////////////////////////////////////////////////////////
    //
    // Start tracking partices
    //
    ////////////////////////////////////////////////////////////////////////////

    // need to count how many alive in each band so
    // that we don't load memory on subsequent passes
    // when there are none alive. On first pass all neutrons are in band 0
    n_alive[0] = npl;
    for (int i=1; i<nb; ++i) n_alive[i] = 0;

    global_total_alive = npg; // assign a non-zero value
    my_total_alive = npl;

    ts = MPI_Wtime();
    int trips = 0;
    const double absorption_threshold = 1.0/nb;

    // Open a passive target epoch on shared memory windows, since we will call
    // MPI_Win_sync on them to sync shm members. Per MPI-3.0 p.449,
    // "All flush and sync functions can be called only within passive target epochs."
    MPI_Win_lock_all(MPI_MODE_NOCHECK, shmwin1);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, shmwin2);

    while (global_total_alive > 0) {
        trips++;
        // Start the pipepline by communicating band 0
        flip = 1;
        if (is_shm_leader) {
            char* sendbuf = (char*)xsdata;
            xscomm = buf1;
            reqcomm = &req1;
            rc = MPI_Iallgather(sendbuf, slicesizekb, dtype_kbytes, xscomm, slicesizekb, dtype_kbytes, dsmcomm, reqcomm);
            assert(!rc);
        }

        // loop over bands starting with highest energy
        for (int k = 0; k < nb; k++) {
            // SHM mebmers align here since the SHM leader
            // will receive band k+1 and overwrite one of the two bufs. We need to
            // make sure all SHM members are done with the data in the buf.
            MPI_Barrier(shmcomm);

            // SHM leader starts communicating band k + 1 within DSM
            if (is_shm_leader && (k+1 < nb)) {
                xscomm = flip ? buf2 : buf1;
                reqcomm = flip ? &req2 : &req1;
                char* sendbuf = (char*)xsdata + slicesizekb*1024*(k+1);
                rc = MPI_Iallgather(sendbuf, slicesizekb, dtype_kbytes, xscomm, slicesizekb, dtype_kbytes, dsmcomm, reqcomm);
                assert(!rc);
            }

            //if (!my_global_rank) printf("rank 0 is working on band %d\n", k);

            // Get xswork and then reverse
            reqwork = flip ? &req1 : &req2;
            winwork = flip ? shmwin1 : shmwin2;
            xswork = flip ? buf1 : buf2;
            flip = !flip;

            // SHM members wait for their leader to complete nonblocking comm.
            if (is_shm_leader) {
                MPI_Wait(reqwork, MPI_STATUS_IGNORE);
                MPI_Win_sync(winwork); // Push stores to the window from private to public
            }
            MPI_Barrier(shmcomm);
            //    MPI_Win_sync(winwork); // Push stores to the window from private to public
            if (!is_shm_leader) MPI_Win_sync(winwork); // Pull stores to the window from public to private

            assert(*(int*)xswork == k); // Make sure we are using correct cross sections

            if (n_alive[k] != 0) { // skip if no alive particles in band k
                for (int i=0; i< npl; ++i) {
                    while (p[i].band == k && !p[i].absorbed) {
                        ran_val = drand48();
                        if (ran_val <= absorption_threshold) {
                            usleep( (1.0/tracking_rate)*1000000 ); // sleep in microseconds
                            p[i].absorbed = TRUE;
                            --n_alive[k];
                        } else {
                            int j = 0;
                            ran_val = drand48(); // ran_val is in [0.0, 1.0)
                            probability = scattering_matrix[k][j];
                            while (ran_val + epsilon > probability && j < nb) {
                                ++j;
                                probability = scattering_matrix[k][j];
                            }
                            p[i].band = j; // move particle to band to which it was scattered
                            --n_alive[k]; ++n_alive[j]; // adjust band-dependent alive counts
                        }
                    }
                    // To make progress of nonblocking comm. of band k+1
                    int flag;
                    if (is_shm_leader && i % 64 == 0) MPI_Test(reqcomm, &flag, MPI_STATUS_IGNORE);
                }
            }
        }
        my_total_alive = tot(n_alive, nb);  // total alive (on one process) across all bands after one sweep
        // TODO: Do allreduce on all processes within dsmcomm instead of MPI_COMM_WORLD, since
        // particle tracking of different dsmcomms is independant.
        MPI_Allreduce(&my_total_alive, &global_total_alive, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    MPI_Win_unlock_all(shmwin1);
    MPI_Win_unlock_all(shmwin2);

    tf = MPI_Wtime();
    ttot = tf - ts;

    DoubleIntPair in, minout, maxout;
    in.val = ttot;
    in.rank = my_global_rank;
    MPI_Reduce(&ttot, &tave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&in, &maxout, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    MPI_Reduce(&in, &minout, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

    if (my_global_rank== 0){
        printf("nprocs = %d\n", nprocs);
        printf("npg = %.1f million\n", npg/1000000.0);
        printf("xsdata = %.1fGB\n", gsizekb/1048576.0);
        printf("energy bands per DSM(i.e., nb) = %d\n", nb);
        printf("energy bands per node = %d\n", nb/r);
        printf("nodes in a DSM (i.e., r) = %d\n", r);
        printf("message size per ibcast = %.1f(MB)\n\n", bandsizekb/1024.0);

        if (use_file)
            printf("scattering matrix file = %s\n", sm_file);
        else
            printf("scattering matrix is self-created\n");

        printf("input tracking rate = %.1f particles/s\n", tracking_rate);
        printf("Go through the energy bands %d times\n", trips);
        printf("Average total tracking time(s)         = %.6f\n", tave/nprocs);
        printf("Minimal total tracking time on rank %d = %.6f\n", minout.rank, minout.val);
        printf("Maximal total tracking time on rank %d = Tebmc(s) = %.6f\n", maxout.rank, maxout.val);
        printf("Tclassic(s)                            = %.1f\n", npg/tracking_rate/nprocs);
        printf("Tebmc/Tclassic                         = %.1f\n", maxout.val/(npg/tracking_rate/nprocs));
        fflush(stdout);
    }

    free(n_alive);
    free(ranges);
    free(p);
    MPI_Comm_free(&shmcomm);
    if (is_shm_leader) {
        MPI_Comm_free(&shm_leader_comm);
        MPI_Comm_free(&dsmcomm);
        free(xsdata);
    }
    MPI_Win_free(&shmwin1);
    MPI_Win_free(&shmwin2);
    matrix_free(scattering_matrix, 0, nb-1, 0, nb-1);
    MPI_Type_free(&dtype_kbytes);
    MPI_Finalize();
    return 0;
}
