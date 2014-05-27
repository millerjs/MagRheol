/******************************************************************************/
/*                                                                            */
/* FILE    : cll.c                                                            */
/* DESCRIP : header file for linked lists                                     */
/* AUTHOR  : Joshua Miller                                                    */
/* PROJECT : Project 4                                                        */
/* CLASS   : Parallel Computing - Winter 2014                                 */
/*                                                                            */
/******************************************************************************/


#ifndef _CLL_H
#define _CLL_H

#include <stdlib.h>
#include <unistd.h>

#include "libjosh.h"

#define ERR_LISTFULL -2
#define ERR_LISTEMPTY -3

typedef struct ll_node_t ll_node_t;

struct ll_node_t {
    int key;
    void *data;
    ll_node_t *next;
} ;

typedef struct ll_t {
    int len;
    int maxlen;
    ll_node_t *head;
    ll_node_t *tail;    
} ll_t; 

ll_t * ll_new(int maxlen);;
int ll_push(ll_t *list, void *item, int key);
void *ll_pop(ll_t *list);
int ll_contains(ll_t *list, int key);
void ll_print(ll_t *list);
void *ll_remove(ll_t *list, int key);
ll_node_t *ll_popnode(ll_t *list);
int ll_pushnode(ll_t *list, ll_node_t *new);
int ll_len(ll_t *list);
ll_node_t *ll_node_new(void *item, int key);
void *ll_Lamport_pop(ll_t *list);

#endif
