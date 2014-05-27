/******************************************************************************/
/*                                                                            */
/* FILE    : cll.c                                                            */
/* DESCRIP : header file for linked lists                                     */
/* AUTHOR  : Joshua Miller                                                    */
/* PROJECT : Project 4                                                        */
/* CLASS   : Parallel Computing - Winter 2014                                 */
/*                                                                            */
/******************************************************************************/


#include <stdlib.h>

#include "libjosh.h"
#include "cll.h"

ll_t *ll_new(int maxlen)
{
    ll_t * new = Malloc(ll_t, 1);
    new->maxlen = maxlen;
    new->len = 0;
    new->head = NULL;
    new->tail = NULL;
    ERROR_IF(!new, ERR_NOMEM);
    return new;
}

ll_node_t *ll_node_new(void *item, int key)
{
   ll_node_t * new = Malloc(ll_node_t, 1);
   new->key = key;
   new->data = item;
   new->next = NULL;
   return new;
}

int ll_len(ll_t *list)
{
   return list->len;
}

int ll_push(ll_t *list, void *item, int key)
{
    ERROR_IF(!list, ERR_NOMEM);
    ll_node_t *new = ll_node_new(item, key);
    ll_node_t *tail = list->tail;
    list->tail = new;
    if (tail) tail->next = new;    
    else      list->head = new;
    __sync_fetch_and_add(&list->len, 1);
    return 0;
}

int ll_pushnode(ll_t *list, ll_node_t *new)
{
    ERROR_IF(!list, ERR_NOMEM);
    ll_node_t *tail = list->tail;
    list->tail = new;
    if (tail) tail->next = new;    
    else      list->head = new;
    __sync_fetch_and_add(&list->len, 1);
    return list->len;
}

void *ll_Lamport_pop(ll_t *list)
{

    if (!list) return NULL;
    
    ERROR_IF(!list, ERR_NOMEM);

    if (list->head == list->tail)
        return NULL;
    if (!list->head)
        return NULL;

    ll_node_t * head = list->head;
    list->head = list->head->next;
    void *data = head->data;
    __sync_fetch_and_sub(&list->len, 1);
    return data;
}

void *ll_pop(ll_t *list)
{
    ERROR_IF(!list, ERR_NOMEM);

    ll_node_t * head = list->head;

    if (!head) return NULL;

    list->head = list->head->next;
    void *data = head->data;
    __sync_fetch_and_sub(&list->len, 1);
    return data;
}

ll_node_t *ll_popnode(ll_t *list)
{
    ERROR_IF(!list, ERR_NOMEM);

    ll_node_t * head = list->head;

    if (!head) return NULL;

    list->head = list->head->next;

    if (!list->head) list->tail = NULL;
    __sync_fetch_and_sub(&list->len, 1);
    
    return head;
}


int ll_contains(ll_t *list, int key)
{
    ERROR_IF(!list, "passed null list");
    ll_node_t *curr = list->head;
    while (curr){
        if (curr->key == key){
            return 1;
        }
        curr = curr->next;
    }
    return 0;
}

void ll_print(ll_t *list)
{
    if (!list) return;
    ll_node_t *curr = list->head;
    while (curr){
        fprintf(stderr, "[%d]->", curr->key);
        curr = curr->next;
    }
}

void *ll_remove(ll_t *list, int key)
{
    ERROR_IF(!list, "passed null list");
    ll_node_t *curr = list->head;
    ll_node_t *last = NULL;
    while (curr){
        if (curr->key == key){
            __sync_fetch_and_sub(&list->len, 1);
            if (last) {
                last->next = curr->next;
            } else {
                list->head = curr->next;
            }
            if (list->tail == curr)
                list->tail = last;
            void *data = curr->data;
            return data;
        }
        last = curr;
        curr = curr->next;
    }
    return NULL;
}
