#ifndef MPISTUFF_H
#define MPISTUFF_H

#ifdef NO_MPI
# include "mpi_tommie.h"
#else
# include <mpi.h>
#endif

#ifdef MPI_OFFSET_IS_LONG_LONG_INT
# define MPI_OFFSET_DATATYPE MPI_LONG_LONG
#else
# define MPI_OFFSET_DATATYPE MPI_LONG
#endif

namespace MpiStuff {

  // these namespace members are declared "extern" so we can include
  // the namespace definition in a header file included by
  // multiple routines
  
  extern int mpi_rank;
  extern int mpi_size;
  extern MPI_Comm mpi_comm;

  // call this method if you really are running mpi...
  extern void initMpiStuff();
  extern void initMpiStuff(MPI_Comm& comm);
  extern void MPI_Pause(char * message);

};

#endif
