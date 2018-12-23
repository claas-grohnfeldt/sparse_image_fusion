#include "mpi_counter.h"

using namespace std;

struct mpi_counter_t *create_counter(int hostrank) {
    struct mpi_counter_t *count;

    count = (struct mpi_counter_t *)malloc(sizeof(struct mpi_counter_t));
    count->hostrank = hostrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &(count->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(count->size));
    

    if (count->rank == hostrank) {
        MPI_Alloc_mem(count->size * sizeof(int), MPI_INFO_NULL, &(count->data));
        for (int i=0; i<count->size; i++) count->data[i] = 0;
        MPI_Win_create(count->data, count->size * sizeof(int), sizeof(int),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &(count->win));
    } else {
        count->data = NULL;
        MPI_Win_create(count->data, 0, 1,
                       MPI_INFO_NULL, MPI_COMM_WORLD, &(count->win));
    }
    count -> myval = 0;

    return count;
}

struct mpi_counter_t *create_counter(int hostrank, int initialValue, MPI_Comm comm) {
    struct mpi_counter_t *count;

    count = (struct mpi_counter_t *)malloc(sizeof(struct mpi_counter_t));
    count->hostrank = hostrank;
    MPI_Comm_rank(comm, &(count->rank));
    MPI_Comm_size(comm, &(count->size));
    

    if (count->rank == hostrank) {
        MPI_Alloc_mem(count->size * sizeof(int), MPI_INFO_NULL, &(count->data));
        for (int i=0; i<count->size; i++) count->data[i] = 0;
        MPI_Win_create(count->data, count->size * sizeof(int), sizeof(int),
                       MPI_INFO_NULL, comm, &(count->win));
    } else {
        count->data = NULL;
        MPI_Win_create(count->data, 0, 1,
                       MPI_INFO_NULL, comm, &(count->win));
    }
    count -> myval = initialValue;

    return count;
}

int increment_counter(struct mpi_counter_t *count, int increment) {
  
    int *vals = (int *)malloc( count->size * sizeof(int) );
    int val;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, count->hostrank, 0, count->win);

    for (int i=0; i<count->size; i++) {

        if (i == count->rank) {
            MPI_Accumulate(&increment, 1, MPI_INT, 0, i, 1, MPI_INT, MPI_SUM,
                           count->win);
        } else {
            MPI_Get(&vals[i], 1, MPI_INT, 0, i, 1, MPI_INT, count->win);
        }
    }

    MPI_Win_unlock(0, count->win);
    
    count->myval += increment;
    
    vals[count->rank] = count->myval;
    val = 0;
    for (int i=0; i<count->size; i++)
        val += vals[i];
    free(vals);
    return val;
}

void delete_counter(struct mpi_counter_t **count) {
    if ((*count)->rank == (*count)->hostrank) {
        MPI_Free_mem((*count)->data);
    }
    MPI_Win_free(&((*count)->win));
    free((*count));
    *count = NULL;

    return;
}

void print_counter(struct mpi_counter_t *count) {
    if (count->rank == count->hostrank) {
        for (int i=0; i<count->size; i++) {
            printf("%2d ", count->data[i]);
        }
        puts("");
    }
}