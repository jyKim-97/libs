#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mt64.h"
#include "networks.h"

// initialize random
// init_genrand64((int)time(NULL));

void InitNetwork(Network *ntk, int N){
    ntk->N = N;
    ntk->adj_list = (int**) malloc(sizeof(int*) * N);
    ntk->n_edges = (int*) malloc(sizeof(int) * N);
    for (int i=0; i<N; i++){
        ntk->adj_list[i] = (int*) malloc(sizeof(int));
        ntk->n_edges[i] = 0;
    }
    ntk->cluster_id_max = -1;
}

void DelNetwork(Network *ntk){
    int N = ntk->N;
    free(ntk->n_edges);
    for (int i=0; i<N; i++){
        free(ntk->adj_list[i]);
    }
    free(ntk->adj_list);
}
// ------------------Building Network Section
void create_random_network(Network *ntk, double p){
    int N = ntk->N;
    double p_rand;
    // work
    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){
            p_rand = genrand64_real2();
            if (p_rand < p){ // connect two nodes; i-j
                connect_i2j(ntk, i, j);
                connect_i2j(ntk, j, i);
            }
        }
    }
}

void create_random_network_fixed_edges(Network *ntk, double p){
    int num_of_edges, num_edges_i, flag, i, j;
    int N = ntk->N;
    // work
    num_of_edges = (int) (round(N*p*N/2));
    while (num_of_edges != 0){
        i = rand_val_max(N);
        j = rand_val_max(N);
        flag = 0;
        num_edges_i = ntk->n_edges[i];
        if (i != j){ // check self-loop
            for (int k=0; k<num_edges_i; k++){ // check does the edge exist
                if (ntk->adj_list[i][k] == j){
                    flag = 1;
                    break;
                }
            }
            if (flag == 0){
                connect_i2j(ntk, i, j);
                connect_i2j(ntk, j, i);
                num_of_edges--;
            }
        }
    }
}

void create_BA_network(Network *ntk, int m){
    int sum_edges=0, n_b, n_acc;
    int N = ntk->N;
    int *num_nodes = (int*) malloc(sizeof(int) * N); // save node number
    int *pre_nodes = (int*) malloc(sizeof(int) * m);
    int sum_pre_nodes;
    // prepare 
    for (int i=0; i<N; i++){
        num_nodes[i] = 1;
    }
    // step1. connect 3 nodes
    for (int i=0; i<m+1; i++){
        for (int j=0; j<m+1; j++){
            if (i != j){
                connect_i2j(ntk, i, j);
                connect_i2j(ntk, j, i);
                sum_edges += 2;
            }
        }
    }
    // step2. connect nodes
    for (int i=m+1; i<N; i++){
        sum_pre_nodes = 0;
        for (int n=0; n<m; n++){ // connect n nodes
            n_acc = rand_val_max(sum_edges-sum_pre_nodes); // [0, N)
            pre_nodes[n] = -1;
            n_b = 1;
            for (int j=0; j<i; j++){
                if (num_nodes[j] == 1){ // if the # node == -1 (used) -> pass
                    n_b += ntk->n_edges[j];
                    if (n_acc <= n_b){
                        connect_i2j(ntk, i, j);
                        connect_i2j(ntk, j, i);
                        sum_edges += 2;
                        // prev nodes
                        sum_pre_nodes += ntk->n_edges[j];
                        pre_nodes[n] = j;
                        num_nodes[j] = -1;
                        break;
                    }
                }
            }
        }
        // undo
        for (int n=0; n<m; n++){
            num_nodes[pre_nodes[n]] = 1;
        }
    }
    free(num_nodes);
    free(pre_nodes);
}

// ------------------Network Algorithm


void calc_bc_once(Network ntk, int s0, double **g_list){
    /*
    calculate betweenness centrality from "source (s)"
    add to *
    correspond to frac{\sigma_{s?}(v)}{\sigma_{s?}
    */
    int flag=1, nr0, nr1, nr0_nxt, nr1_nxt, N=ntk.N;
    int s, t, dptr;
    int sz_root = N;
    int *num_used_nodes = (int*) malloc(sizeof(int) * N);
    Root *root_list = (Root*) malloc(sizeof(Root) * N);
    Root *root_list_prev;
    Root *up_root_adr;
    //
    for (int i=0; i<N; i++) num_used_nodes[i] = 0;
    num_used_nodes[s0] = 1;
    // initialize
    root_list[0].node = s0;
    root_list[0].front = NULL;
    nr0 = 0; nr1 = 1;
    // run
    while (flag == 1){
        flag = 0;
        // search
        nr0_nxt = nr1;
        nr1_nxt = nr1;
        for (int n=nr0; n<nr1; n++){
            s = root_list[n].node; // node number
            for (int i=0; i<ntk.n_edges[s]; i++){
                t = ntk.adj_list[s][i];
                if (num_used_nodes[t] == 0){ // didn't use
                    flag = 1;
                    root_list[nr1_nxt].node = t;
                    root_list[nr1_nxt].front = &(root_list[n]); // save the upper node address
                    nr1_nxt++;
                    if (nr1_nxt == sz_root){
                        sz_root += N;
                        // realloc을 하면서 주소 바뀔거같은데
                        // root_list = realloc(root_list, sizeof(Root) * sz_root);
                        root_list_prev = root_list;
                        root_list = realloc(root_list, sizeof(Root) * sz_root);
                        if (root_list != root_list_prev){ // exchange front address
                            dptr = root_list - root_list_prev;
                            for (int k=1; k<nr1_nxt; k++){ // except source
                                root_list[k].front += dptr;
                            }
                        }
                    }
                }
            }
        }
        // record used nodes
        for (int n=nr0_nxt; n<nr1_nxt; n++){
            t = root_list[n].node;
            num_used_nodes[t]++;
        }
        // update g_list
        for (int n=nr0_nxt; n<nr1_nxt; n++){
            t = root_list[n].node;
            up_root_adr = root_list[n].front;
            while ((*up_root_adr).front != NULL){ // top root : NULL
                (*g_list)[(*up_root_adr).node] += 1 / (double) num_used_nodes[t];
                up_root_adr = (*up_root_adr).front;
            }
        }
        // exchange
        nr0 = nr0_nxt;
        nr1 = nr1_nxt;
        // update depth
        // d++;
    }
    free(num_used_nodes);
    free(root_list);
}

int *calc_distance(Network ntk, int s0){
    // use BFS Algorithm
    int *ptr0, *ptr1, flag, *nodes_arr, s, t, d=1, N=ntk.N;
    int *d_list;
    // init
    d_list = (int*) malloc(sizeof(int) * N);
    for (int i=0; i<N; i++) d_list[i] = 0;
    d_list[s0] = -1;
    nodes_arr = (int*) malloc(sizeof(int) * N);
    ptr0 = nodes_arr;
    ptr1 = nodes_arr+1;
    *ptr0 = s0;
    // run
    while (ptr1 != ptr0){
        s = *ptr0;
        for (int i=0; i<ntk.n_edges[s]; i++){
            t = ntk.adj_list[s][i];
            if (d_list[t] == 0){
                d_list[t] = d;
                *ptr1 = t;
                ptr1++;
            }
        }
        ptr0++;
        d++;
    }
    // d_list[s0] = 0;
    free(nodes_arr);
    return d_list;
}

void getClusterID(Network *ntk){
    int flag=1, *d_list, N=ntk->N, s0=0;
    ntk->cluster_id = (int*) malloc(sizeof(int) * N);
    for (int i=0; i<N; i++) ntk->cluster_id[i] = -1;
    ntk->cluster_id_max = 0;
    while (flag == 1){
        flag = 0;
        d_list = calc_distance(*ntk, s0);
        for (int i=0; i<N; i++){
            if (d_list[i] != 0){
                ntk->cluster_id[i] = ntk->cluster_id_max;
            }
        }
        ntk->cluster_id_max++;
        // search next id
        for (int i=0; i<N; i++){
            if (ntk->cluster_id[i] == -1){
                flag = 1;
                s0 = i;
                break;
            }
        }
        free(d_list);
    }
    ntk->cluster_id_max--;
}

double *calc_betweenness_centrality(Network ntk){
    double *g_list = (double*) malloc(sizeof(double) * ntk.N);
    for (int i=0; i<ntk.N; i++) g_list[i] = 0;
    // run
    for (int s=0; s<ntk.N; s++){
        calc_bc_once(ntk, s, &g_list);
    }
    for (int s=0; s<ntk.N; s++){
        g_list[s] /= 2;
    }
    return g_list;
}
// ------------------Network visualization




// ------------------Sub functions
void connect_i2j(Network *ntk, int i, int j){
    int num_edges_i = ntk->n_edges[i];
    // work
    ntk->adj_list[i][num_edges_i] = j;
    ntk->n_edges[i]++;
    // instead of using num_edges ++, use +2
    ntk->adj_list[i] = (int*) realloc(ntk->adj_list[i], sizeof(int)*(num_edges_i+2));
}

int rand_val_max(int N){
    // return value for [0, N)
    return (int) floor(genrand64_real2() * N);
}

void save_network(Network *ntk, FILE *fid){
    int num_edges_i;
    int N = ntk->N;
    // work
    for (int i=0; i<N; i++){
        fprintf(fid, "%d", i);
        num_edges_i = ntk->n_edges[i];
        for (int j=0; j<num_edges_i; j++){
            fprintf(fid, ";%d", ntk->adj_list[i][j]);
        }
        fprintf(fid, "\n");
    }
    fprintf(fid, "EOF");
}

void read_network(Network *ntk, FILE *fid){
    int N, n, k, row=0, len, t_node;
    char *buffer = (char*) malloc(sizeof(char));
    char num[50];
    ntk->adj_list = (int**) malloc(sizeof(int*));
    ntk->adj_list[0] = (int*) malloc(sizeof(int));
    ntk->n_edges = (int*) malloc(sizeof(int));
    fscanf(fid, "%[^\n]\n", buffer);
    while (strcmp(buffer, "EOF") != 0){
        n = k = 0;
        len = 0;
        if (N > 0){
            ntk->adj_list = realloc(ntk->adj_list, sizeof(int*) * (N+1));
            ntk->adj_list[N] = (int*) malloc(sizeof(int));
            ntk->n_edges = realloc(ntk->n_edges, sizeof(int) * (N+1));
        }
        while (*(buffer+n) != '\0'){
            num[k] = *(buffer+n);
            if (num[k] == ';'){
                num[k] = '\0';
                t_node = atoi(num);
                k = 0;
                if (len > 0){
                    ntk->adj_list[row][len-1] = t_node;
                    ntk->adj_list[row] = realloc(ntk->adj_list[row], sizeof(int) * (len+1));
                }
                len++;
                n++;
                continue;
            }
            k++;
            n++;
        }
        ntk->adj_list[row][len-1] = atoi(num);
        ntk->n_edges[N] = len;
        row++;
        fscanf(fid, "%[^\n]\n", buffer);
        N++;
    }
    ntk->N = N;
}
