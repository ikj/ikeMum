/*
 * JoeUgpWithCvCompFlow.h
 *
 *  Created on: Sep 3, 2014
 *      Author: ikj
 */

// This file will be included in JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h

#ifndef JOEUGPWITHCVCOMPFLOW_H_
#define JOEUGPWITHCVCOMPFLOW_H_

// =========================
//  UTILITY FUNCTIONS
// =========================

/*
 * Method: showBarrieredMessage
 * ----------------------------
 * One CPU core shows its message (usually, error or warning message).
 * MPI_Barrier() function will be called both at the beginning and in the end.
 *
 * Note: 1. All the CPU cores must call this method (Otherwise, communication error)
 *       2. If the string is empty, NO message will be shown
 */
void showBarrieredMessage(const char message[], const int mpi_rank_to_show) {
	MPI_Barrier(mpi_comm);
	for(int i=0; i<100000; ++i) {} // Pause for a while
	if(mpi_rank == mpi_rank_to_show)
		cout<<message<<endl;
	MPI_Barrier(mpi_comm);
}
void showBarrieredMessage(string& message, const int mpi_rank_to_show) {
	MPI_Barrier(mpi_comm);
	for(int i=0; i<100000; ++i) {} // Pause for a while
	if(mpi_rank == mpi_rank_to_show)
		cout<<message<<endl;
	MPI_Barrier(mpi_comm);
}

#endif /* JOEUGPWITHCVCOMPFLOW_H_ */
