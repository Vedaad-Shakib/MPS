/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "que.h": Queue implementation
**
*******************************************************************************/

/*******************************************************************************
 * Queue struct definition
 *******************************************************************************
 */
typedef struct {
  int		head;     	/* the head of the queue		*/
  int		tail;     	/* the tail of the queue		*/
  int		maxSize;  	/* the maximum size of the queue	*/
  double*	body;		/* storage for queue elements		*/
} Que;

typedef Que* QueHd;

/*******************************************************************************
 * Abstract functions
 *******************************************************************************
 */

QueHd queNew();
bool queIsEmpty(QueHd q);
int queGetSize(QueHd q);
void queExpand(QueHd q);
void quePush(QueHd q, double val);
double quePop(QueHd q);
void queFree(QueHd q);
