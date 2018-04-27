// Matt Matto
// Lab B - The Bonding Lab
// 
// DESCRIPTION: Implements pthreads to determine which Hydrogen and Oxygen 
//              molecules to bond together into water molecules
//
// USAGE: ./bonding seed num_molecules max_outstanding verbosity

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <dllist.h>
#include "bonding.h"

struct global_info {
  pthread_mutex_t *lock;  // mutex lock to protecc the Dllists
  Dllist waitingHydro;    // Dllist for hydrogen atoms waiting to be bonded
  Dllist waitingOxy;      // Dllist for oxygen atoms waiting to be bonded 
};

struct thread_info {
  pthread_cond_t *cond;   // conditino variable to be signal when a thread should call Bond()
  int tid;                // thread id
  int h1;                 // thread id of the first hydrogen atom
  int h2;                 // thread id of the second hydrogen atom
  int o;                  // thread id of the oxygen atom
};

void *initialize_v(char *verbosity) {
  struct global_info *g;  // global information shared between all atom threads

  // malloc & initialize the global info
  g = (struct global_info *) malloc(sizeof(struct global_info));
  g->lock = new_mutex();
  g->waitingHydro = new_dllist();
  g->waitingOxy = new_dllist();

  return (void *) g;
}

void *hydrogen(void *arg) {
  struct bonding_arg *b;        // pointer to the passed data structure
  struct global_info *g;        // pointer to the global info
  struct thread_info *t;        // this hydrogen atom's thread information
  struct thread_info *e1, *e2;  // pointers to this threads this hydrogen will bond to 
  Dllist tmp1, tmp2;            // used to extract thread information from Dllists
  char *rv;                     // return value string

  // set b, g, and initialize t
  b = (struct bonding_arg *) arg;
  g = (struct global_info *) b->v;
  t = (struct thread_info *) malloc(sizeof(struct thread_info));
  t->cond = new_cond();
  t->tid = b->id;

  // lock the thread and check if I can make a water molecule
  pthread_mutex_lock(g->lock);
  if((dll_first(g->waitingHydro) != dll_nil(g->waitingHydro)) && (dll_first(g->waitingOxy) != dll_nil(g->waitingOxy))) {
    
    // pulls the first H and O off each list
    tmp1 = dll_first(g->waitingHydro);
    e1 = tmp1->val.v;
    tmp2 = dll_first(g->waitingOxy);
    e2 = tmp2->val.v;
    
    // set tids for each thread
    t->h1 = t->tid;
    t->h2 = e1->tid;
    t->o = e2->tid;
    e1->h1 = t->tid;
    e1->h2 = e1->tid;
    e1->o = e2->tid;
    e2->h1 = t->tid;
    e2->h2 = e1->tid;
    e2->o = e2->tid;

    // remove the threads from the lists
    dll_delete_node(tmp1);
    dll_delete_node(tmp2);

    // signal the other threads
    pthread_cond_signal(e1->cond);
    pthread_cond_signal(e2->cond);

    // unlock and Bond()
    pthread_mutex_unlock(g->lock);
    rv = Bond(t->h1, t->h2, t->o);

    // free malloc'ed stuff
    free(t->cond);
    free(t);

    return rv;
  } else {
    // else we could not make a new water molecule
    // add this thread to the waiting list, then block
    dll_append(g->waitingHydro, new_jval_v((void *) t));
    pthread_cond_wait(t->cond, g->lock);
    pthread_mutex_unlock(g->lock);

    // after getting a signal, Bond
    rv = Bond(t->h1, t->h2, t->o);

    // free malloc'ed stuff
    free(t->cond);
    free(t);

    return rv;
  }
}

// works very similar to the Hydrogen atom thread
void *oxygen(void *arg) {
  struct bonding_arg *b;        // pointer to the passed data structure
  struct global_info *g;        // pointer to the global info
  struct thread_info *t;        // this hydrogen atom's thread information
  struct thread_info *e1, *e2;  // pointers to this threads this hydrogen will bond to 
  Dllist tmp1, tmp2;            // used to extract thread information from Dllists
  char *rv;                     // return value string

  // set b, g, and initialize t
  b = (struct bonding_arg *) arg;
  g = (struct global_info *) b->v;
  t = (struct thread_info *) malloc(sizeof(struct thread_info));
  t->cond = new_cond();
  t->tid = b->id;

  // lock the thread and check if I can make a water molecule
  pthread_mutex_lock(g->lock);
  if ((dll_first(g->waitingHydro) != dll_nil(g->waitingHydro)) && (dll_next(dll_first(g->waitingHydro)) != dll_nil(g->waitingHydro))) {
    
    // pulls the first and second H off the H list
    tmp1 = dll_first(g->waitingHydro);
    e1 = tmp1->val.v;
    tmp2 = dll_next(dll_first(g->waitingHydro));
    e2 = tmp2->val.v;
    
    // set tids for each thread
    t->h1 = e1->tid;
    t->h2 = e2->tid;
    t->o = t->tid;
    e1->h1 = e1->tid;
    e1->h2 = e2->tid;
    e1->o = t->tid;
    e2->h1 = e1->tid;
    e2->h2 = e2->tid;
    e2->o = t->tid;

    // remove the threads from the lists
    dll_delete_node(tmp1);
    dll_delete_node(tmp2);

    // signal the other threads
    pthread_cond_signal(e1->cond);
    pthread_cond_signal(e2->cond);

    // unlock and Bond()
    pthread_mutex_unlock(g->lock);
    rv = Bond(t->h1, t->h2, t->o);
   
    // free malloc'ed stuff
    free(t->cond);
    free(t);
    
    return rv;
  } else {
    // else we could not make a new water molecule
    // add this thread to the waiting list, then block
    dll_append(g->waitingOxy, new_jval_v((void *) t));
    pthread_cond_wait(t->cond, g->lock);
    pthread_mutex_unlock(g->lock);
    
    // after getting a signal, Bond
    rv = Bond(t->h1, t->h2, t->o);
    
    // free malloc'ed stuff
    free(t->cond);
    free(t);

    return rv;
  }
}
