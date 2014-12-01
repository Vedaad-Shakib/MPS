/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "queue.c": Queue implementation
**
*******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************/

#include "sys.h"
#include "que.h"

/*******************************************************************************
 * "queNew": creates a new queue
 *******************************************************************************/
QueHd queNew() {
  QueHd queHd; // queue structure

  /* allocate structure */
  queHd = memNew(Que, 1);

  queHd->head = 0;
  queHd->tail = -1;
  queHd->maxSize = 2;
  queHd->body = memNew(double, queHd->maxSize);

  /* return structure */
  return (queHd);
}

/*******************************************************************************
 * "queIsEmpty": checks if the queue is empty
 *
 * Parameters:
 *            queHd: the queue to be checked
 *******************************************************************************/
bool queIsEmpty(QueHd queHd) {
  return (queHd->head > queHd->tail);
}

/*******************************************************************************
 * "queGetSize": returns the size of the queue
 *
 * Parameters:
 *            queHd: the queue with size to be calculated
 *******************************************************************************/
int queGetSize(QueHd queHd) {
  return (queHd->tail-queHd->head+1);
}

/*******************************************************************************
 * "queExpand": expands the queue by its maximum size + 100
 *
 * Parameters:
 *            queHd: the queue to be expanded
 *******************************************************************************/
void queExpand(QueHd queHd) {
  double *newBody; /* the expanded body */
  double tmpHead; /* the temporary head */
  double tmpTail; /* the temporary tail */
  int maxSize; /* new max size */

  maxSize = queHd->maxSize * 2 + 100 ;

  newBody = memNew(double, maxSize);
  queHd->maxSize = maxSize ;
  tmpHead = queHd->head;
  tmpTail = queHd->tail;
  queHd->head = 0;
  queHd->tail = -1;

  for (int i = tmpHead; i <= tmpTail; i++) {
    queHd->tail++;
    newBody[queHd->tail] = queHd->body[i];
  }

  free(queHd->body);
  queHd->body = newBody;
}

/*******************************************************************************
 * "quePush": inserts an element at the end of the queue
 *
 * Parameters:
 *            queHd: the queue to which the element is added
 *            val: the value to be added
 *******************************************************************************/
void quePush(QueHd queHd, double val) {
  if ( queHd->tail == queHd->maxSize-1 ) {
    queExpand(queHd);
  }

  queHd->tail++;
  queHd->body[queHd->tail] = val;
}

/*******************************************************************************
 * "quePop": deletes an element from the start of the queue
 *
 * Parameters:
 *            queHd: the queue to which the element is added
 *******************************************************************************/
double quePop(QueHd queHd) {
  assert(!queIsEmpty(queHd));
  queHd->head++;
  return queHd->body[queHd->head-1];
}

/*******************************************************************************
 * "quePrintQueue": prints the queue for debugging purposes
 *
 * Parameters:
 *            queHd: the queue to be printed
 *******************************************************************************/
void quePrintQueue(QueHd queHd) {
  printf("\n");
  for (int i = queHd->head; i <= queHd->tail; i++) {
    printf("%le\n", queHd->body[i]);
  }
}

/*******************************************************************************
 * "queFree": frees the memory of a queue
 *
 * Parameters:
 *            queHd: the queue to be freed
 *******************************************************************************/
void queFree(QueHd queHd) {
  free(queHd->body);
  free(queHd);
}
