/*
    2014/10/01 Naoto Kawahara
  mpi_mm.c for CAPI

  This code from github.
  https://github.com/hiranya911/bsp-mpi/blob/master/mpi/mpi_mm.cpp#L176
  //                    Matrix Multiplication with MPI                       //
  //                         CS240A Assignment 1                             //
  //                                                                         //
  // Author: Hiranya Jayathilaka (hiranya@cs.ucsb.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include "carbon_user.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include "omp.h"

#define N 256
#define CORE 4

using namespace std;

int num_threads;
carbon_barrier_t rank_barrier;
int send_completed[CORE] ;
double* A ;//= new double[N * N];
double* B ;//= new double[N * N];

void* mpi_mm(void* threadid);
void do_worker(long rank);
void check_results(double* C, double* D);
void print(double* matrix, int n);

void CAPI_message_receive_w_bl(int sender, int receiver, char *buffer, int size, int *completed);
void do_bcast(int rank);
void receive_from_all(int rank);

void* null_thread(void* threadid);

int main(int argc, char **argv)
{
   int limit_num = N*N;
   int worker_count;
   A = new double[limit_num];
   B = new double[limit_num];

  #pragma omp parallel for
  for(int i=0; i < limit_num; i++){
    A[i] =  2 * drand48( ) - 1;
    B[i] =  2 * drand48( ) - 1;
  }    

  /*
  for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
         A[i*N + j] = 2 * drand48( ) - 1;
         B[i*N + j] = 2 * drand48( ) - 1;
       }
     }
  */
 
  
  CarbonStartSim(argc, argv);
  printf("Starting mpi_mm!\n");

  if (argc != 2)
   {
      fprintf(stderr, "Usage: ./mpi_mm <num_threads>\n");
      CarbonStopSim();
      exit(-1);
   }
 
  num_threads = atoi(argv[1]) - 1;
  worker_count = CORE;

   carbon_thread_t threads[worker_count];

   CarbonBarrierInit(&rank_barrier, worker_count);
   
   for(int i = 0; i < worker_count; i++){
     send_completed[i] = 0;
     threads[i] = CarbonSpawnThread(mpi_mm, (void *) i);
   }
   
   for(int i = worker_count; i < num_threads; i++){
     threads[i] = CarbonSpawnThread(null_thread, (void *) i);
   }
   
   for(int i = 0; i < num_threads; i++)
     CarbonJoinThread(threads[i]);
   
   printf("Ending all_to_all!\n");
   CarbonStopSim();

   delete [] A;
   delete [] B;

   return 0;
}

void* mpi_mm(void* threadid)
{
  long rank = (long)threadid;
  int process_count = CORE;
  double mflops;
  UInt64 start_time, end_time, time;
  
  CAPI_Initialize(rank);
   if (rank == 0) {
     //double* C = new double[N * N];
     double* D = new double[N * N];
     int rows, offset;
     int worker_count = process_count;// - 1;

     cout << "Worker count: " << worker_count << endl;

    /*
     start_time = CarbonGetTime();
     
     for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
	 double cij = 0.0;
	 for (int k = 0; k < N; k++) {
	   cij += A[i*N + k] * B[k*N + j];
	 }
	 C[i*N + j] = cij;
       }
     }

     UInt64 end_time = CarbonGetTime();
     UInt64 time = end_time - start_time;

     double mflops = (2e+4 * N * N * N) /(double)(time);
     cout << endl;
     cout << "Serial Execution" << endl;
     cout << "================" << endl;
     cout << "Time Elapsed: " << time << " nanoseconds" << endl;
     cout << "Mflop/s: " << mflops << endl;
    */    
 
     CarbonBarrierWait(&rank_barrier);
     
     int rows_per_worker = N / worker_count;
     int extra_rows = N % worker_count;
     offset = 0;
     rows = N /worker_count;     

     start_time = CarbonGetTime();
     
     //for (int worker = worker_count - 1 ; worker >= 0; worker--) {
     for(int worker = 0; worker < worker_count ; worker++){
       // Try to equally divide the available rows among workers
       //rows = (worker <= extra_rows) ? rows_per_worker + 1 : rows_per_worker; 

       CAPI_message_send_w(0, worker , (char*) &rows, sizeof(int));
       CAPI_message_send_w(0, worker , (char*) &offset, sizeof(int));
       CAPI_message_send_w(0, worker , (char*) &A[offset*N], sizeof(double)*rows*N);
       CAPI_message_send_w(0, worker , (char*) B, sizeof(double)*N*N);
      // cout << rank << ")Send: " << worker << " " << rows << " " << offset << " " << endl;
       send_completed[worker] = 1;
       offset += rows;
     }
 
    //cout << rank << ")Send: " << rows << " " << offset << " " << endl;
     do_worker(rank);
    
     rows = 0;
     offset = 0;
    //for(int worker = 0; worker < worker_count ; worker++){
    for (int worker = worker_count - 1 ; worker >= 0; worker--) {
       // Retrieve the results from workers

       CAPI_message_receive_w_bl(worker, 0, (char *) &rows, sizeof(int), send_completed + worker);
       CAPI_message_receive_w(worker, 0, (char *) &offset, sizeof(int));
      // cout << rank << ")Receive: " << worker << " " << rows << " " << offset << " " << endl;
       CAPI_message_receive_w(worker, 0, (char *) &D[offset*N], sizeof(double)*rows*N);
       //print(D,N);
       send_completed[worker] = 0;
     }
     
     end_time = CarbonGetTime();
     time = end_time - start_time;
   mflops = (2e-6 * N * N * N) /((double)(time)/1000000000.0);
   cout << endl;
     cout << "Parallel Execution" << endl;
     cout << "==================" << endl;
     cout << "Time Elapsed: " << time << " nanoseconds" << endl;
     cout << "Mflop/s: " << mflops << endl;

     //check_results(C, D);
     
     //delete [] C;
     delete [] D;
   } else {
       CarbonBarrierWait(&rank_barrier);
     do_worker(rank);
   }
   
   return NULL;
}

/*
 * This functon should be only called by worker nodes
 */
void do_worker(long rank) {
  int rows, offset;
  // MPI_Status status;
  double* AA = new double[N * N];
  double* BB = new double[N * N];
  double* CC = new double[N * N];

  CAPI_message_receive_w_bl(0, rank, (char *) &rows, sizeof(int),send_completed + rank);
  send_completed[rank] = 0;
  CAPI_message_receive_w(0, rank, (char *) &offset, sizeof(int));
  CAPI_message_receive_w(0, rank, (char *) AA, sizeof(double)*rows*N);
  CAPI_message_receive_w(0, rank, (char *) BB, sizeof(double)*N*N);
  //cout << rank << ")Receive: " << rows << " " << offset << " " << endl;
  
  // Multiply the set of rows assigned to me
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < N; j++) {
      double cij = 0.0;
      for (int k = 0; k < N; k++) {
	cij += AA[i*N + k] * BB[k*N + j];
      }
      CC[i*N + j] = cij;
    }
  }
  
  CAPI_message_send_w(rank, 0 , (char *) &rows, sizeof(int));
  CAPI_message_send_w(rank, 0 , (char *) &offset, sizeof(int));
  CAPI_message_send_w(rank, 0 , (char *) CC, sizeof(double)*rows*N);
  //cout << rank << ")Send: 0 " << rows << " " << offset << " " << endl;
  send_completed[rank] = 1;

  //print(C,rows);
  
  delete [] AA;
  delete [] BB;
  delete [] CC;
}

void do_bcast(int rank)
{
  int sent_message_payload = rank; 

   CAPI_endpoint_t sender_rank = (CAPI_endpoint_t) rank;
   CAPI_endpoint_t receiver_rank = (CAPI_endpoint_t) CAPI_ENDPOINT_ALL;
		
   printf("[do_bcast] rank %d is broadcasting out message; value = %d\n", rank, sent_message_payload);
   CAPI_message_send_w(sender_rank, receiver_rank, (char *) &sent_message_payload, sizeof(int));
   printf("[do_bcast] rank %d finished broadcasting out message\n", rank);
}

void receive_from_all(int rank)
{
   for (int i = 0; i < num_threads; i++)
   {
      int received_message_payload = -1;
      CAPI_endpoint_t sender_rank = (CAPI_endpoint_t) CAPI_ENDPOINT_ANY;
      CAPI_endpoint_t receiver_rank = (CAPI_endpoint_t) rank;

      printf("[receive_from_all] rank %d is about to receive message\n", rank);
      CAPI_message_receive_w(sender_rank, receiver_rank, (char*) &received_message_payload, sizeof(int));
      printf("[receive_from_all] rank %d received message; value = %d\n", rank, received_message_payload);
   }
}


void CAPI_message_receive_w_bl(int sender, int receiver, char *buffer, int size, int *completed)
{
  while(!*completed){
    usleep(2);
  }
  CAPI_message_receive_w(sender, receiver, buffer, size);
}

void check_results(double* C, double* D) {
  cout << endl << "Verifying the multiplication results..." << endl;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double diff = C[i*N + j] - D[i*N + j];
      if (diff < 0) {
	diff = diff * -1;
      }
      if (diff > 0.0001) {
	cout << "Multiplication error detected" << endl;
	cout << C[i*N+j] << ", " << D[i*N+j] << endl;
	return;
      }
    }
  }
  cout << "All good..." << endl;
  return;
}

void print(double* matrix, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << matrix[i*n + j] << " ";
    }
    cout << endl;
  }
  cout << endl << endl;
}


void* null_thread(void* threadid)
{
  return NULL;
}
