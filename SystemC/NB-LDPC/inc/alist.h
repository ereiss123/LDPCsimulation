/*==========================================================================================
** alist.h
**
** Date: March, 2024
**
** Authors: Eric Reiss, Chris Winstead
**          Utah State University
**
** Based on original source by Radford Neal and David MacKay
**
** Description:
**   Defines struct and functions for processing non-binary "alist" files, 
**   which represent the positions and values of non-zero elements in a 
**   sparse matrix. The alist file is typically used to specify NB-LDPC 
**   codes.
==============================================================================================*/

#ifndef ALIST_H
#define ALIST_H

#include <cstdlib>
#include <cstdio>
#include "r.h"

typedef struct {
	int N , M ;      /* size of the matrix */
	int **mlist;     /* list of integer coordinates in the m direction where the non-zero entries are */
	int **mvals;     /* list of integer values corresponding to the mlist positions                   */
	int **nlist;     /* list of integer coordinates in the n direction where the non-zero entries are */
	int **nvals;     /* list of integer values corresponding to the nlist positions                   */
	int *num_mlist;  /* weight of each row, m */
	int *num_nlist;  /* weight of each column n */
	int *l_up_to ;
	int *u_up_to ;
	int *norder ;
	int biggest_num_m ;       /* actual biggest sizes             */
	int biggest_num_n ;
	int biggest_num_m_alloc ; /* sizes used for memory allocation */
	int biggest_num_n_alloc ;
	int tot ;
	int same_length ;         /* whether all vectors in mlist and nlist have same length */
	int q;                    /* Order of the Galois field                               */
} alist_struct ;


alist_struct loadFile(const char * fileName);
void printAlist(alist_struct alist);
void freeAlist(alist_struct alist);
int fread_nbmatrix 
(
 int **b ,
 int **c ,
 int l1,
 int h1,
 int l2,
 int h2,
 FILE *fp 
 );

#endif
