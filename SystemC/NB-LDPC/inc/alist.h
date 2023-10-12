/*==========================================================================================
** alist.h
** By Chris Winstead
   Based on original code by Radford Neal and David MacKay

** Description: 
   Defines struct and functions for processing "alist" files, which
   represent the locations of 1's in a large sparse binary matrix.
   The alist file is typically used to specify LDPC codes.

** Usage:
    ldpcsim <alist_fname> <stim_fname> <iterations> <clock_cycles> <vcd_fname>
==============================================================================================*/

#ifndef ALIST_H
#define ALIST_H

//#include "string.h"


typedef struct {
	int N , M ;      /* size of the matrix */
	int **mlist;     /* list of integer coordinates in the m direction where the non-zero entries are */
	int **nlist;     /* list of integer coordinates in the n direction where the non-zero entries are */
	int *num_mlist;  /* weight of each row, m */
	int *num_nlist;  /* weight of each column n */
	int *l_up_to ;
	int *u_up_to ;
	int *norder ;
	int biggest_num_m ;       /* actual biggest sizes */
	int biggest_num_n ; 
	int biggest_num_m_alloc ; /* sizes used for memory allocation */
	int biggest_num_n_alloc ; 
	int tot ; 
	int same_length ;  /* whether all vectors in mlist and nlist have same length */
} alist_struct ;


alist_struct loadFile(const char * fileName);
void printAlist(alist_struct alist);
void freeAlist(alist_struct alist);

#endif
