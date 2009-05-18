#include "routine_replace.h"
#include "simulator.h"
#include "thread_manager.h"
#include "log.h"
#include "carbon_user.h"
#include "thread_support_private.h"

// -------------------------------------
// Begin Memory redirection stuff
#include "config.h"
#include "simulator.h"
#include "core.h"
#include "core_manager.h"
#include "redirect_memory.h"
#include "thread_start.h"
// End Memory redirection stuff
// --------------------------------------


int CarbonPthreadCreateWrapperReplacement(CONTEXT *ctx, AFUNPTR orig_fp, void *pthread_t_p, void *pthread_attr_t_p, void *routine_p, void* arg_p)
{
   fprintf(stderr, "Create pthread called from pin.\n");
   // Get the function for the thread spawner
   PIN_LockClient();
   AFUNPTR pthread_create_function;
   IMG img = IMG_FindByAddress((ADDRINT)orig_fp);
   RTN rtn = RTN_FindByName(img, "pthread_create");
   pthread_create_function = RTN_Funptr(rtn);
   PIN_UnlockClient();

   fprintf(stderr, "pthread_create_function: %x\n", (int)pthread_create_function);

   int res;
   PIN_CallApplicationFunction(ctx,
         PIN_ThreadId(),
         CALLINGSTD_DEFAULT,
         pthread_create_function,
         PIN_PARG(int), &res,
         PIN_PARG(void*), pthread_t_p,
         PIN_PARG(void*), pthread_attr_t_p,
         PIN_PARG(void*), routine_p,
         PIN_PARG(void*), arg_p,
         PIN_PARG_END());

   return res;
}

// ---------------------------------------------------------
// Memory Redirection
//
// Routine replacement needs to work a little differently
// We can't use RTN_Replace or RTN_ReplaceSignature because
// the application executes out of a different memory space
// (and hence also uses a different stack) that Pin isn't aware of
// Our solution is to insert a call before the routine executes 
// to a wrapper function that extracts the correct arguments from 
// the user stack (including populating memory pointers with data
// from simulated memory) and then calls the "replacement functions"
// within the pintool. This wrapper also writes back results to user
// memory before returning control to the correct return address
// in user space
//
// --------------------------------------------------------------

bool replaceUserAPIFunction(RTN& rtn, string& name)
{
   AFUNPTR msg_ptr = NULL;
   PROTO proto = NULL;

   // TODO: Check that the starting stack is located below the text segment
   // thread management
   if (name == "_start") msg_ptr = AFUNPTR (replacement_start);
   else if (name == "main") msg_ptr = AFUNPTR (replacementMain);
   else if (name == "CarbonGetThreadToSpawn") msg_ptr = AFUNPTR(replacementGetThreadToSpawn);
   else if (name == "CarbonThreadStart") msg_ptr = AFUNPTR (replacementThreadStartNull);
   else if (name == "CarbonThreadExit") msg_ptr = AFUNPTR (replacementThreadExitNull);
   else if (name == "CarbonGetCoreId") msg_ptr = AFUNPTR(replacementGetCoreId);
   else if (name == "CarbonDequeueThreadSpawnReq") msg_ptr = AFUNPTR (replacementDequeueThreadSpawnRequest);

   // PIN specific stack management
   else if (name == "CarbonPthreadAttrInitOtherAttr") msg_ptr = AFUNPTR(replacementPthreadAttrInitOtherAttr);

   // Carbon API

   else if (name == "CarbonStartSim") msg_ptr = AFUNPTR(replacementStartSimNull); 
   else if (name == "CarbonStopSim") msg_ptr = AFUNPTR(replacementStopSim);
   else if (name == "CarbonSpawnThread") msg_ptr = AFUNPTR(replacementSpawnThread);
   else if (name == "CarbonJoinThread") msg_ptr = AFUNPTR(replacementJoinThread);

   // CAPI
   else if (name == "CAPI_Initialize") msg_ptr = AFUNPTR(replacement_CAPI_Initialize);
   else if (name == "CAPI_rank") msg_ptr = AFUNPTR(replacement_CAPI_rank);
   else if (name == "CAPI_message_send_w") msg_ptr = AFUNPTR(replacement_CAPI_message_send_w);
   else if (name == "CAPI_message_receive_w") msg_ptr = AFUNPTR(replacement_CAPI_message_receive_w);

   // synchronization
   else if (name == "CarbonMutexInit") msg_ptr = AFUNPTR(replacementMutexInit);
   else if (name == "CarbonMutexLock") msg_ptr = AFUNPTR(replacementMutexLock);
   else if (name == "CarbonMutexUnlock") msg_ptr = AFUNPTR(replacementMutexUnlock);
   else if (name == "CarbonCondInit") msg_ptr = AFUNPTR(replacementCondInit);
   else if (name == "CarbonCondWait") msg_ptr = AFUNPTR(replacementCondWait);
   else if (name == "CarbonCondSignal") msg_ptr = AFUNPTR(replacementCondSignal);
   else if (name == "CarbonCondBroadcast") msg_ptr = AFUNPTR(replacementCondBroadcast);
   else if (name == "CarbonBarrierInit") msg_ptr = AFUNPTR(replacementBarrierInit);
   else if (name == "CarbonBarrierWait") msg_ptr = AFUNPTR(replacementBarrierWait);

   // pthread wrappers
//   else if (name == "CarbonPthreadCreateWrapper") msg_ptr = AFUNPTR(replacementPthreadCreateWrapperReplacement);
//   else if (name.find("pthread_create") != std::string::npos) msg_ptr = AFUNPTR(replacementPthreadCreate);
//   else if (name.find("pthread_join") != std::string::npos) msg_ptr = AFUNPTR(replacementPthreadJoin);
   else if (name.find("pthread_exit") != std::string::npos) msg_ptr = AFUNPTR(replacementPthreadExitNull);

   // do replacement
   if (msg_ptr == AFUNPTR(CarbonPthreadCreateWrapperReplacement))
   {
      proto = PROTO_Allocate(PIN_PARG(int),
                             CALLINGSTD_DEFAULT,
                             name.c_str(),
                             PIN_PARG(void*),
                             PIN_PARG(void*),
                             PIN_PARG(void*),
                             PIN_PARG(void*),
                             PIN_PARG_END());
      RTN_ReplaceSignature(rtn, msg_ptr,
                           IARG_PROTOTYPE, proto,
                           IARG_CONTEXT,
                           IARG_ORIG_FUNCPTR,
                           IARG_FUNCARG_ENTRYPOINT_VALUE, 0,
                           IARG_FUNCARG_ENTRYPOINT_VALUE, 1,
                           IARG_FUNCARG_ENTRYPOINT_VALUE, 2,
                           IARG_FUNCARG_ENTRYPOINT_VALUE, 3,
                           IARG_END);
      PROTO_Free(proto);
      return true;
   }
   else if (msg_ptr != NULL)
   {
      RTN_Open (rtn);

      RTN_InsertCall (rtn, IPOINT_BEFORE,
            msg_ptr,
            IARG_CONTEXT,
            IARG_END);

      RTN_Close (rtn);

      return true;
   }
   else
   {
      return false;
   }
}

// void SimSpawnThreadSpawner(CONTEXT *ctx, AFUNPTR fp_main)
// {
//    // Get the function for the thread spawner
//    PIN_LockClient();
//    AFUNPTR thread_spawner;
//    IMG img = IMG_FindByAddress((ADDRINT)fp_main);
//    RTN rtn = RTN_FindByName(img, "CarbonSpawnThreadSpawner");
//    thread_spawner = RTN_Funptr(rtn);
//    PIN_UnlockClient();
// 
//    // Get the address of the thread spawner
//    int res;
//    PIN_CallApplicationFunction(ctx,
//          PIN_ThreadId(),
//          CALLINGSTD_DEFAULT,
//          thread_spawner,
//          PIN_PARG(int), &res,
//          PIN_PARG_END());
// 
// }

void replacement_start (CONTEXT *ctxt)
{
   // FIXME: 
   return;


   if (Sim()->getConfig()->getCurrentProcessNum() == 0)
      return;

   else
   {
      Core *core = Sim()->getCoreManager()->getCurrentCore();
      core->getNetwork()->netRecv (0, SYSTEM_INITIALIZATION_NOTIFY);
      
      int res;
      ADDRINT reg_eip = PIN_GetContextReg (ctxt, REG_INST_PTR);

      PIN_LockClient();

      AFUNPTR thread_spawner;
      IMG img = IMG_FindByAddress(reg_eip);
      RTN rtn = RTN_FindByName(img, "CarbonThreadSpawner");
      thread_spawner = RTN_Funptr(rtn);

      PIN_UnlockClient();
      
      PIN_CallApplicationFunction (ctxt,
            PIN_ThreadId(),
            CALLINGSTD_DEFAULT,
            thread_spawner,
            PIN_PARG(int), &res,
            PIN_PARG(void*), NULL,
            PIN_PARG_END());

      // Should have some ack for the core_manager and thread_manager here
      // to notify them that the core is finished and is exiting
      exit (0);
   }
}


void replacementMain (CONTEXT *ctxt)
{
   if (Sim()->getConfig()->getCurrentProcessNum() == 0)
   {
      spawnThreadSpawner(ctxt);

      Core *core = Sim()->getCoreManager()->getCurrentCore();
      UInt32 num_processes = Sim()->getConfig()->getProcessCount();
      for (UInt32 i = 1; i < num_processes; i++)
      {
         core->getNetwork()->netSend (Sim()->getConfig()->getThreadSpawnerCoreNum (i), SYSTEM_INITIALIZATION_NOTIFY, NULL, 0);
      }

      return;
   }
   else
   {
      int res;
      ADDRINT reg_eip = PIN_GetContextReg (ctxt, REG_INST_PTR);

      PIN_LockClient();

      AFUNPTR thread_spawner;
      IMG img = IMG_FindByAddress(reg_eip);
      RTN rtn = RTN_FindByName(img, "CarbonThreadSpawner");
      thread_spawner = RTN_Funptr(rtn);

      PIN_UnlockClient();
      
      PIN_CallApplicationFunction (ctxt,
            PIN_ThreadId(),
            CALLINGSTD_DEFAULT,
            thread_spawner,
            PIN_PARG(int), &res,
            PIN_PARG(void*), NULL,
            PIN_PARG_END());

      Sim()->getThreadManager()->onThreadExit();

      while (!Sim()->finished())
         sched_yield();

      Simulator::release();

      exit (0);
   }
}

void replacementGetThreadToSpawn (CONTEXT *ctxt)
{
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core != NULL);
   
   ThreadSpawnRequest *req;
   ThreadSpawnRequest req_buf;
   initialize_replacement_args (ctxt,
         IARG_PTR, &req,
         IARG_END);
   
   // Preserve REG_GAX across function call with
   // void return type
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);

   CarbonGetThreadToSpawn (&req_buf);

   core->accessMemory(Core::NONE, WRITE, (IntPtr) req, (char*) &req_buf, sizeof (ThreadSpawnRequest));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementThreadStartNull (CONTEXT *ctxt)
{
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   retFromReplacedRtn (ctxt, ret_val);
}

void replacementThreadExitNull (CONTEXT *ctxt)
{
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   retFromReplacedRtn (ctxt, ret_val);
}

void replacementGetCoreId (CONTEXT *ctxt)
{
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   CarbonGetCoreId ();

   retFromReplacedRtn (ctxt, ret_val);
}


void replacementDequeueThreadSpawnRequest (CONTEXT *ctxt)
{
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);

   ThreadSpawnRequest *thread_req;
   ThreadSpawnRequest thread_req_buf;
   
   initialize_replacement_args (ctxt,
         IARG_PTR, &thread_req,
         IARG_END);
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   CarbonDequeueThreadSpawnReq (&thread_req_buf);

   core->accessMemory (Core::NONE, WRITE, (IntPtr) thread_req, (char*) &thread_req_buf, sizeof (ThreadSpawnRequest));

   retFromReplacedRtn (ctxt, ret_val);
}

// PIN specific stack management
void replacementPthreadAttrInitOtherAttr(CONTEXT *ctxt)
{
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert(core != NULL);

   pthread_attr_t *attr;
   pthread_attr_t attr_buf;

   initialize_replacement_args(ctxt,
         IARG_PTR, &attr,
         IARG_END);

   core->accessMemory (Core::NONE, READ, (ADDRINT) attr, (char*) &attr_buf, sizeof (pthread_attr_t));

   ADDRINT ret_val = PIN_GetContextReg(ctxt, REG_GAX);

   SimPthreadAttrInitOtherAttr(&attr_buf);

   core->accessMemory (Core::NONE, WRITE, (ADDRINT) attr, (char*) &attr_buf, sizeof (pthread_attr_t));

   retFromReplacedRtn(ctxt, ret_val);
}


void replacementStartSimNull (CONTEXT *ctxt)
{
   ADDRINT ret_val = 0; 
   retFromReplacedRtn (ctxt, ret_val);
}

void replacementStopSim (CONTEXT *ctxt)
{
   terminateThreadSpawner ();

   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   retFromReplacedRtn (ctxt, ret_val);
}

void replacementSpawnThread (CONTEXT *ctxt)
{
   thread_func_t func;
   void *arg;

   initialize_replacement_args (ctxt,
         IARG_PTR, &func,
         IARG_PTR, &arg,
         IARG_END);

   LOG_PRINT("Calling SimSpawnThread");
   ADDRINT ret_val = (ADDRINT) CarbonSpawnThread (func, arg);

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementJoinThread (CONTEXT *ctxt)
{
   ADDRINT tid;

   initialize_replacement_args (ctxt,
         IARG_ADDRINT, &tid,
         IARG_END);

   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);

   CarbonJoinThread ((int) tid);

   retFromReplacedRtn (ctxt, ret_val);
}

void replacement_CAPI_Initialize (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the CAPI communication API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   int comm_id;
   initialize_replacement_args (ctxt,
         IARG_UINT32, &comm_id,
         IARG_END);

   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);

   CAPI_Initialize (comm_id);

   retFromReplacedRtn (ctxt, ret_val);
}

void replacement_CAPI_rank (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the CAPI communication API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   int *core_id;
   int core_id_buf;
   CAPI_return_t ret_val;

   initialize_replacement_args (ctxt,
         IARG_PTR, &core_id,
         IARG_END);

   ret_val = CAPI_rank (&core_id_buf);
   
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) core_id, (char*) &core_id_buf, sizeof (int));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacement_CAPI_message_send_w (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the CAPI communication API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   CAPI_endpoint_t sender;
   CAPI_endpoint_t receiver;
   char *buffer;
   int size;
   CAPI_return_t ret_val = 0;

   initialize_replacement_args (ctxt,
         IARG_ADDRINT, &sender,
         IARG_ADDRINT, &receiver,
         IARG_PTR, &buffer,
         IARG_ADDRINT, &size,
         IARG_END);

   char *buf = new char [size];
   core->accessMemory (Core::NONE, READ, (ADDRINT) buffer, buf, size);
   ret_val = CAPI_message_send_w (sender, receiver, buf, size);

   retFromReplacedRtn (ctxt, ret_val);
}

void replacement_CAPI_message_receive_w (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the CAPI communication API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   CAPI_endpoint_t sender;
   CAPI_endpoint_t receiver;
   char *buffer;
   int size;
   CAPI_return_t ret_val = 0;

   initialize_replacement_args (ctxt,
         IARG_ADDRINT, &sender,
         IARG_ADDRINT, &receiver,
         IARG_PTR, &buffer,
         IARG_ADDRINT, &size,
         IARG_END);

   char *buf = new char [size];
   ret_val = CAPI_message_receive_w (sender, receiver, buf, size);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) buffer, buf, size);

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementMutexInit (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_mutex_t *mux;
   initialize_replacement_args (ctxt,
         IARG_PTR, &mux,
         IARG_END);
   
   carbon_mutex_t mux_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));
   CarbonMutexInit (&mux_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementMutexLock (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_mutex_t *mux;
   initialize_replacement_args (ctxt,
         IARG_PTR, &mux,
         IARG_END);
   
   carbon_mutex_t mux_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));
   CarbonMutexLock (&mux_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementMutexUnlock (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_mutex_t *mux;
   initialize_replacement_args (ctxt,
         IARG_PTR, &mux,
         IARG_END);
   
   carbon_mutex_t mux_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));
   CarbonMutexUnlock (&mux_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementCondInit (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_cond_t *cond;
   initialize_replacement_args (ctxt,
         IARG_PTR, &cond,
         IARG_END);
   
   carbon_cond_t cond_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));
   CarbonCondInit (&cond_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementCondWait (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_cond_t *cond;
   carbon_mutex_t *mux;
   initialize_replacement_args (ctxt,
         IARG_PTR, &cond,
         IARG_PTR, &mux,
         IARG_END);
   
   carbon_cond_t cond_buf;
   carbon_mutex_t mux_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));
   core->accessMemory (Core::NONE, READ, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));
   CarbonCondWait (&cond_buf, &mux_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) mux, (char*) &mux_buf, sizeof (mux));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementCondSignal (CONTEXT *ctxt)
{
   // Only the user-threads (all of which are cores) call
   // the Carbon synch API functions
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   
   carbon_cond_t *cond;
   initialize_replacement_args (ctxt,
         IARG_PTR, &cond,
         IARG_END);
   
   carbon_cond_t cond_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   core->accessMemory (Core::NONE, READ, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));
   CarbonCondSignal (&cond_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementCondBroadcast (CONTEXT *ctxt)
{
   carbon_cond_t *cond;
   initialize_replacement_args (ctxt,
         IARG_PTR, &cond,
         IARG_END);
   
   carbon_cond_t cond_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   core->accessMemory (Core::NONE, READ, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));
   CarbonCondBroadcast (&cond_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) cond, (char*) &cond_buf, sizeof (cond));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementBarrierInit (CONTEXT *ctxt)
{
   carbon_barrier_t *barrier;
   UINT32 count;
   initialize_replacement_args (ctxt,
         IARG_PTR, &barrier,
         IARG_UINT32, &count,
         IARG_END);
   
   carbon_barrier_t barrier_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   core->accessMemory (Core::NONE, READ, (ADDRINT) barrier, (char*) &barrier_buf, sizeof (barrier));
   CarbonBarrierInit (&barrier_buf, count);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) barrier, (char*) &barrier_buf, sizeof (barrier));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementBarrierWait (CONTEXT *ctxt)
{
   carbon_barrier_t *barrier;
   initialize_replacement_args (ctxt,
         IARG_ADDRINT, &barrier,
         IARG_END);
   
   carbon_barrier_t barrier_buf;
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   
   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);
   core->accessMemory (Core::NONE, READ, (ADDRINT) barrier, (char*) &barrier_buf, sizeof (barrier));
   CarbonBarrierWait (&barrier_buf);
   core->accessMemory (Core::NONE, WRITE, (ADDRINT) barrier, (char*) &barrier_buf, sizeof (barrier));

   retFromReplacedRtn (ctxt, ret_val);
}

void replacementPthreadExitNull (CONTEXT *ctxt)
{
   ADDRINT ret_val = PIN_GetContextReg (ctxt, REG_GAX);
   retFromReplacedRtn (ctxt, ret_val);
}

void initialize_replacement_args (CONTEXT *ctxt, ...)
{
   va_list vl;
   va_start (vl, ctxt);
   int type;
   ADDRINT addr;
   ADDRINT ptr;
   ADDRINT buffer;
   unsigned int count = 0;
   // Core *core = Sim()->getCoreManager()->getCurrentCore();
   // assert (core);
   do
   {
      type = va_arg (vl, int);
      addr = PIN_GetContextReg (ctxt, REG_STACK_PTR) + ((count + 1) * sizeof (ADDRINT));
     
      memOp (Core::NONE, READ, addr, (char*) &buffer, sizeof (ADDRINT));
      switch (type)
      {
         case IARG_ADDRINT:
            ptr = va_arg (vl, ADDRINT);
            * ((ADDRINT*) ptr) = buffer;
            count++;
            break;

         case IARG_PTR:
            ptr = va_arg (vl, ADDRINT);
            * ((ADDRINT*) ptr) = buffer;
            count++;
            break;

         case IARG_UINT32:
            ptr = va_arg (vl, ADDRINT);
            * ((UINT32*) ptr) = (UINT32) buffer;
            break;

         case IARG_END:
            break;

         default:
            assert (false);
            break;
      }
   } while (type != IARG_END);
}

void retFromReplacedRtn (CONTEXT *ctxt, ADDRINT ret_val)
{
   ADDRINT esp = PIN_GetContextReg (ctxt, REG_STACK_PTR);
   ADDRINT next_ip = emuRet (&esp, 0, sizeof (ADDRINT));

   PIN_SetContextReg (ctxt, REG_GAX, ret_val);
   PIN_SetContextReg (ctxt, REG_STACK_PTR, esp);
   PIN_SetContextReg (ctxt, REG_INST_PTR, next_ip);

   PIN_ExecuteAt (ctxt);
}

void setupCarbonSpawnThreadSpawnerStack (CONTEXT *ctx)
{
   // FIXME: 
   // This will clearly need to change somewhat in the multi-process case
   // We can go back to our original scheme of having the "main" thread 
   // on processes other than 0 execute the thread spawner, in which case
   // this will probably just work as is

   ADDRINT esp = PIN_GetContextReg (ctx, REG_STACK_PTR);
   ADDRINT ret_ip = * (ADDRINT*) esp;

   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);

   core->accessMemory(Core::NONE, WRITE, (IntPtr) esp, (char*) &ret_ip, sizeof (ADDRINT));
}

void setupCarbonThreadSpawnerStack (CONTEXT *ctx)
{
   if (Sim()->getConfig()->getCurrentProcessNum() == 0)
      return;

   ADDRINT esp = PIN_GetContextReg (ctx, REG_STACK_PTR);
   ADDRINT ret_ip = * (ADDRINT*) esp;
   ADDRINT p = * (ADDRINT*) (esp + sizeof (ADDRINT));

   Core *core = Sim()->getCoreManager()->getCurrentCore();
   assert (core);

   core->accessMemory (Core::NONE, WRITE, (IntPtr) esp, (char*) &ret_ip, sizeof (ADDRINT));
   core->accessMemory (Core::NONE, WRITE, (IntPtr) (esp + sizeof (ADDRINT)), (char*) &p, sizeof (ADDRINT));
}


