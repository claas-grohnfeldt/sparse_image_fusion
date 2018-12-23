/** Solution of an MPI global counter adapted from Jonathan Dursi
 *  Source: http://stackoverflow.com/questions/4948788/creating-a-counter-that-stays-synchronized-across-mpi-processes
 **/

#ifndef MPI_COUNTER_H
#define MPI_COUNTER_H

#include "includes.h"

// counter struct
struct mpi_counter_t {
    MPI_Win win;
    int  hostrank ;
    int  myval;
    int *data;
    int rank, size;
};

struct mpi_counter_t *create_counter(int hostrank);
struct mpi_counter_t *create_counter(int hostrank, int initialValue, MPI_Comm comm);
int increment_counter(struct mpi_counter_t *count, int increment);
void delete_counter(struct mpi_counter_t **count);
void print_counter(struct mpi_counter_t *count);

#endif // MPI_COUNTER_H
