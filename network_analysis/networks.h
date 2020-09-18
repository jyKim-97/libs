#ifndef FUNC_CALL

#define FUNC_CALL
#include <stdio.h>

// typedef struct _Network Network;
typedef struct _Network{
    int N;
    int **adj_list;
    int *n_edges;
    int cluster_id_max;
    int *cluster_id;
} Network;

typedef struct _Root{
    int node;
    struct _Root *front;
} Root;

// 
void InitNetwork(Network *ntk, int N);
void DelNetwork(Network *ntk);
// create ntk
void create_random_network(Network *ntk, double p);
void create_random_network_fixed_edges(Network *ntk, double p);
void create_BA_network(Network *ntk, int m);
// caclulate network properties
void getClusterID(Network *ntk);
double *calc_betweenness_centrality(Network ntk);
int *calc_distance(Network ntk, int s0);
void calc_bc_once(Network ntk, int s0, double **g_list);
// subfunctions
void connect_i2j(Network *ntk, int i, int j);
int rand_val_max(int N);
void save_network(Network *ntk, FILE *fid);
void read_network(Network *ntk, FILE *fid);

#endif
