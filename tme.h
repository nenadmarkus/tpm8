#pragma once

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <malloc.h>

/*
	
*/

#define _FIXED_POINT_SCALE_ 256

void compute_rcs_transformation(int* T, int r, int c, int s);
void compute_rcso_transformation(int* T, int r, int c, int s, float o);

/*
	
*/

int learn_template(int32_t template[], int maxnumtests, int useorientation, float s2p, int r, int c, int s, uint8_t pixels[], uint8_t mask[], int nrows, int ncols, int ldim, int threshold);

#define SAVE_TEMPLATE(t, file) fwrite((t), sizeof(int32_t), (t)[0]+1, (file))
#define LOAD_TEMPLATE(t, file) fread(&(t)[0], sizeof(int32_t), 1, (file)), ((t)[0]==0)?0:fread(&(t)[1], sizeof(int32_t), (t)[0], (file))

int match_template_at(int32_t template[], int threshold, int* T, int* pn1, int n0max, uint8_t pixels[], int nrows, int ncols, int ldim);

/*
	
*/

typedef struct _tnode
{
	int leaf, tag;

	int32_t* template;

	struct _tnode* subtree1;
	struct _tnode* subtree2;

} tnode;

tnode* grow_tree(int* Ts[], int32_t* templates[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int tags[], int inds[], int numinds, int32_t tcodepool[], int tcodepoolsize, int maxnumtestspernode, int perturbationstrength);

tnode* load_tree_from_file(FILE* file);
int save_tree_to_file(tnode* root, FILE* file);

int get_tree_output(tnode* root, int threshold, int n0max, int* T, uint8_t pixels[], int nrows, int ncols, int ldim, int tags[], int* n);
