/*
  2014/10/22 Naoto Kawahara
  main.c(sort program) for CAPI

  This code from masato yoshimi OpenMP.
 */

#include <stdio.h>
#include <stdlib.h>
#include "carbon_user.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <math.h>
#include <time.h>
//#include "mpi.h"

#define N 1024*64
#define CORE 8
#define ASC 0
#define DESC 1
//#define DEBUG 1

using namespace std;

int num_threads;
carbon_barrier_t rank_barrier;
int send_completed[CORE] ;

void* mpi_mm(void* threadid);
void do_worker(long rank);
void check_results(double* C, double* D);
void print(double* matrix, int n);

void CAPI_message_receive_w_bl(int sender, int receiver, char *buffer, int size, int *completed);
void do_bcast(int rank);
void receive_from_all(int rank);
void merge_sort(float sa1[], float sa2[], float da1[], float da2[], int num);
void quick_sort(float a[], int left, int right, int flag);
float med3(float x, float y, float z);
double get_dtime(void);
void* null_thread(void* threadid);

int main(int argc, char **argv)
{
  int worker_count;

   CarbonStartSim(argc, argv);
   printf("Starting mpi_mm!\n");
 
   if (argc != 2)
   {
      fprintf(stderr, "Usage: ./mpi_sort <num_threads>\n");
      CarbonStopSim();
      exit(-1);
   }
   
   num_threads = atoi(argv[1]) -1;
   worker_count = CORE ;

   carbon_thread_t threads[worker_count];

   CarbonBarrierInit(&rank_barrier, worker_count);

   for(int i = 0; i < worker_count; i++)
   {
     send_completed[i] = 0;
     printf("Spawning thread: %d\n", i);
     threads[i] = CarbonSpawnThread(mpi_mm, (void *) i);
   }

   for(int i = worker_count; i < num_threads; i++)
   {
     printf("Spawning network thread: %d\n", i);
     threads[i] = CarbonSpawnThread(null_thread, (void *) i);
   }
   
   for(int i = 0; i < num_threads; i++)
     CarbonJoinThread(threads[i]);
   
   printf("Ending mpi_sort!\n");
   CarbonStopSim();
   return 0;
}

void* mpi_mm(void* threadid)
{
  // Vars init
  long rank = (long)threadid;
  int process_count = CORE;
  double mflops;
  UInt64 start_time, end_time, d1_time;
  int i, num = N, k, offset;
  float* array1;
  float* tmp_array;
  float* array2;
  float* array3;
  float* array_quiq;

  array1 = (float *)malloc(sizeof(float)*num);
  array_quiq =  (float *)malloc(sizeof(float)*num/process_count);
  array2 = (float *)malloc(sizeof(float)*num);

  //  array3 = (float *)malloc(sizeof(float)*num);

  int order = ASC;
  
  CAPI_Initialize(rank);
   if (rank == 0) {
     //double d0, d1, d2; //time 

     cout << "Worker count: " << process_count  << endl;
     if(N >= process_count ){
       // cout << argv[0] << " [num]" << endl;
       CarbonStopSim();
     }

     srand((unsigned int)time(NULL));    

     cout << "generate radum number."<< endl;
     for(i=0;i<num;i++){
       array1[i] = (float)rand() / (float)RAND_MAX;
       //array1[i] = 0.0;
     }
    
#ifdef DEBUG
     for(i=0;i<num;i++){
       cout << i << " : " << array1[i] << endl;
     }
#endif
    }

   CarbonBarrierWait(&rank_barrier);
   //array2 = (float *)malloc(sizeof(float)*num);

   //cout << "Start Sorting" << endl;
   if (rank == 0){
     //d0 = get_dtime();    
     start_time = CarbonGetTime();
   }

   int dnum    = process_count;
   int blkelem = num / dnum;

   // distributed send num / process_count
   // send recive
   if (rank == 0){
     for (int worker = dnum - 1; worker >= 0; worker--) {
       // Try to equally divide the available rows among works
       //rows = (worker <= extra_rows) ? rows_per_worker + 1 : rows_per_worker;
       offset = blkelem * worker;
       CAPI_message_send_w(0, worker , (char*) &array1[offset], sizeof(float)*blkelem);
       send_completed[worker] = 1;
       //offset += rows;
     }
#ifdef DEBUG
     cout << "End Sending data" << endl;
#endif
   }

   CAPI_message_receive_w_bl(0, rank, (char *) &array_quiq[0], sizeof(float)*blkelem,send_completed + rank);
   send_completed[rank] = 0;

   quick_sort(array_quiq, 0, blkelem-1, ASC);  
   

   /*
   //int offset = num * (rank - 1)
     
   // distributed work
   if(rank == 0){
     //for(k=0;k<dnum;k++){
     //}
   //CarbonBarrierWait(&rank_barrier);
   }   
   else{
     //cout << "Quick_sort" << endl;
     quick_sort(array1, 0, blkelem-1, ASC);
     //CarbonBarrierWait(&rank_barrier);
     //cout << "end quick_sort" <<endl;
     //d1_time = CarbonGetTime();
     */
     
     /*
#ifdef DEBUG       
	 for(int i=0;i<num;i++){
	   cout << "[" << i << "] " << array1[i] << endl;
	 }
#endif 
     */       
   //   }

   //   CarbonBarrierWait(&rank_barrier);
   CAPI_message_send_w(rank, 0, (char *) &array_quiq[0], sizeof(float)*blkelem);
   send_completed[rank] = 1;

    // merge data
   if (rank == 0){
     for (int worker = dnum-1; worker >= 0; worker--) {
      // Try to equally divide the available rows among workers                                                                        
       offset = blkelem * worker;
       CAPI_message_receive_w_bl(worker, 0, (char *) &array1[offset], sizeof(float)*blkelem, send_completed + worker);
        send_completed[worker] = 0;
	}
     d1_time = CarbonGetTime();
     }
   
   /*
#ifdef DEBUG       
       if(rank ==0){
	 for(int i=0;i<num;i++){
	   cout << "[" << i << "] " << array1[i] << endl;
	 }
       }
#endif        
   */
   
   //bitonic sort step2
   for(int fb=1; fb<=(int)log2(dnum); fb++) {
     for(int sb=fb-1; sb>=0; sb--) {
       //distributed work
       /*
	 for(int i = 0; i < dnum; i++) {
	 if((i>>fb)&1^(i>>sb)&1){
	 //cout << "I tell number rank (i^(1<<sb)), i"<< rank << " " << (i^(1<<sb)) << " " << i<< endl;
	 merge_sort(&array1[(i^(1<<sb))*blkelem],  &array1[i*blkelem],
	 &array2[(i^(1<<sb))*blkelem], &array2[i*blkelem],
	 blkelem);
	 }
	 }
       */
       for(int i = 0; i < dnum; i++) {
	 if((i>>fb)&1^(i>>sb)&1){
	   if ( rank  == 0 ){
	     //offset = (i^(1<<sb))*blkelem;
	     //cout << "Worker Choice (" <<  i+1 << endl;
	     CAPI_message_send_w(0, i , (char*) &array1[(i^(1<<sb))*blkelem], sizeof(float)*blkelem);
	     CAPI_message_send_w(0, i , (char*) &array1[i*blkelem], sizeof(float)*blkelem);
	     send_completed[i] = 1;
	   }
	   else if( rank == i ){  
	     CAPI_message_receive_w_bl(0, rank, (char *) &array1[(i^(1<<sb))*blkelem], sizeof(float)*blkelem,send_completed + rank);
	     CAPI_message_receive_w(0, rank, (char *) &array1[i*blkelem], sizeof(float)*blkelem);
	     send_completed[rank] = 0;
	     merge_sort(&array1[(i^(1<<sb))*blkelem],  &array1[i*blkelem],
			&array2[(i^(1<<sb))*blkelem], &array2[i*blkelem],
			blkelem);
	   }
#ifdef DEBUG       
	   if ( rank == 0){
	     cout << "End Sending data" << i << endl;
	   }
#endif
	 }
       }

       //CarbonBarrierWait(&rank_barrier);

       for(int i = 0; i < dnum; i++) {
	 if((i>>fb)&1^(i>>sb)&1){
	   if( rank == i ){
	     CAPI_message_send_w(rank, 0, (char*) &array2[(i^(1<<sb))*blkelem], sizeof(float)*blkelem);
	     CAPI_message_send_w(rank , 0, (char*) &array2[i*blkelem], sizeof(float)*blkelem);
	     send_completed[rank] = 1;
	   }
	   else if ( rank  == 0 ) {
	     //cout << "Recive Worker Choice (" <<  i+1 << endl;
	     CAPI_message_receive_w_bl(i, 0, (char *) &array2[(i^(1<<sb))*blkelem], sizeof(float)*blkelem,send_completed + i);
	     CAPI_message_receive_w(i, 0, (char *) &array2[i*blkelem], sizeof(float)*blkelem);
	     send_completed[i] = 0;
	   }
#ifdef DEBUG       
	   if ( rank == 0){
	     cout << "End Small merginge data" << endl;
	   }
#endif
	 }	 
       }
       if ( rank == 0){
	 cout << "End Merge data" << endl;
	        tmp_array = array2;
		array2 = array1;
		array1 = tmp_array;
       }       
       //CarbonBarrierWait(&rank_barrier);
     
       
       /*       
#ifdef DEBUG       
       if(rank ==0){
	 for(int i=0;i<num;i++){
	   cout << "[" << i << "] " << array1[i] << endl;
	 }
       }
#endif        
              //cout << endl;
	      */
     }
   }
   
   /*
   // merge data
   if (rank == 0){
     for (int worker = 1; worker <= dnum; worker++) {
      // Try to equally divide the available rows among workers                                                                        
       offset = dnum * (worker - 1 );
       CAPI_message_receive_w_bl(worker, 0, (char *) &array1[offset], sizeof(float)*dnum, send_completed + worker);
        send_completed[rank] = 0;
	}
     }
     else{
       CAPI_message_send_w(rank, 0, (char *) &array1, sizeof(float)*dnum);
       send_completed[rank] = 1;
     }
   */
     CarbonBarrierWait(&rank_barrier);

     
     // end result check
   if (rank == 0){
     end_time = CarbonGetTime();
     
#ifdef DEBUG
     cout << endl;
     for(i=0;i<num;i++){
       cout << i << " : " << array1[i] << endl;
     }
#endif

     for(i=1;i<num;i++){
       if(order == ASC){
	 if(array1[i-1] > array1[i]){
	   cout << "num[" << i-1 << "] : " << array1[i-1]
		<< " is larger than num[" << i
		<< "] : " << array1[i] <<endl;
	 }
       }
       else{
	 if(array1[i-1] < array1[i]){
	   cout << "num[" << i-1 << "] : " << array1[i-1]
		<< " is smaller than num[" << i
		<< "] : " << array1[i] <<endl;
	 }
       }
     }
     cout << "Sorting  : ,     " << end_time - start_time << ", [sec],  (, " << (d1_time-start_time) << ", " << (end_time-d1_time) << ", )"<< endl;
   }
   
   CarbonBarrierWait(&rank_barrier);
   
   //MPI_Finalize();
   return NULL;
}



/*
 * This functon should be only called by worker nodes
 */
void do_worker(long rank) {

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
    usleep(1);
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

double get_dtime(void){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return ((double)(tv.tv_sec)
            + (double)(tv.tv_usec) * 0.001 * 0.001);
}

float med3(float x, float y, float z){
    float ret;
    if(x<y)
        if(y<z) ret = y; else if(z<x) ret = x; else ret = z;
    else
        if(z<y) ret = y; else if(x<z) ret = x; else ret = z;
    return ret;
}

void quick_sort(float a[], int left, int right, int flag){
    int i, j;
    float tmp, pivot;
    if(left < right){
        i = left; j = right;
        pivot = med3(a[i], a[(i+j)/2], a[j]);
        while(1){
            if(flag == ASC) while(a[i] < pivot) i++;
            else            while(a[i] > pivot) i++;
            if(flag == ASC) while(pivot < a[j]) j--;
            else            while(pivot > a[j]) j--;
            if(i >= j) break;
            tmp = a[i]; a[i] = a[j]; a[j] = tmp;
            i++; j--;
        }
        quick_sort(a, left, i-1, flag);
        quick_sort(a, j+1, right, flag);
    }
}

void merge_sort(float sa1[], float sa2[], float da1[], float da2[], int num){
    int i;
    int ps, qs;
    
    ps = 0; qs = 0;

    for(i=0;i<num*2;i++){
        if(ps >= num){
            if(i < num) da1[i] = sa2[qs];
            else        da2[i-num] = sa2[qs];
            qs++;
        }
        else if(qs >= num){
            if(i < num) da1[i] = sa1[ps];
            else        da2[i-num] = sa1[ps];
            ps++;
        }
        else{
            if(sa1[ps] < sa2[qs]){
                if(i < num) da1[i] = sa1[ps];
                else        da2[i-num] = sa1[ps];
                ps++;
            }
            else{
                if(i < num) da1[i] = sa2[qs];
                else        da2[i-num] = sa2[qs];
                qs++;
            }
        }
    }
}

void* null_thread(void* threadid)
{
  return NULL;
}
