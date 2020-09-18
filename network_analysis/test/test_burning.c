#include <stdlib.h>
#include <stdio.h>

typedef struct _Network{
    int N;
    int **adj_list;
    int *n_edges;
    // clustering coefficient
    int max_cluster;
    int *clustering_coef;
    // max_cluster = -1;
} Network;

void InitNetwork(Network *ntk, int N);
int* BFSalgorithm(Network *ntk, int source);
void calc_clustering_coef(Network *ntk);
void connect_i2j(Network *ntk, int i, int j);

int main(){
    Network ntk;
    InitNetwork(&ntk, 10);
    for (int i=0; i<6; i++){
        for (int j=i+1; j<6; j++){
            connect_i2j(&ntk, i, j);
            connect_i2j(&ntk, j, i);
        }
    }
    printf("===================\n");
    for (int i=0; i<10; i++){
        printf("Node#%d: ", i);
        for (int j=0; j<ntk.n_edges[i]; j++){
            printf("%d,", j);
        }
        printf("\n");
    }
    printf("===================\n");
    calc_clustering_coef(&ntk);
    printf("=====cluster id=====\n");
    for (int i=0; i<ntk.N; i++){
        printf("%d,", ntk.clustering_coef[i]);
    }
    printf("\n===================\n");
    printf("Done\n");
}

void InitNetwork(Network *ntk, int N){
    ntk->N = N;
    ntk->adj_list = (int**) malloc(sizeof(int*) * N);
    ntk->n_edges = (int*) malloc(sizeof(int) * N);
    for (int i=0; i<N; i++){
        ntk->adj_list[i] = (int*) malloc(sizeof(int));
        ntk->n_edges[i] = 0;
    }
}

// ------------------Network Algorithm
void calc_clustering_coef(Network *ntk){
    int *d_list, s, id=0, flag=1;
    int N = ntk->N;
    ntk->clustering_coef = (int*) malloc(sizeof(int)*N);
    for (int i=0; i<N; i++){
        ntk->clustering_coef[i] = -1;
    }
    s = 0;
    while (flag == 1){
        flag = 0;
        // if flag == 0 -> all nodes have cluster id;
        d_list = BFSalgorithm(ntk, s);
        printf("id%d: ", id);
        for (int i=0; i<N; i++){
            printf("%2d, ", d_list[i]);
            if (d_list[i] != -1){
                ntk->clustering_coef[i] = id;
            }
        }
        printf("\n");
        // search next source
        for (s=0; s<N; s++){
            if (ntk->clustering_coef[s] == -1){
                flag = 1;
                break;
            }
        }
        id++;
    }
}


int* BFSalgorithm(Network *ntk, int source){
    // source means source node
    // target means target node
    int *d_list, *ptr_read, *ptr_write, *queue;
    int target, num_edges_i, d = 0;
    int N = ntk->N;
    // initialize
    queue = (int*) malloc(sizeof(int) * N);
    d_list = (int*) malloc(sizeof(int) * N);
    for (int i=0; i<N; i++){
        d_list[i] = -1;
    }
    // save source with distance = 0 (=d)
    queue[0] = source;
    d_list[source] = d;
    ptr_read = queue;
    ptr_write = queue+1;
    // run
    while (ptr_write > ptr_read){
        // search nearest neighbors
        source = *(ptr_read); // start: queue[0]
        d++;
        num_edges_i = ntk->n_edges[source];
        for (int i=0; i<num_edges_i; i++){
            target = ntk->adj_list[source][i];
            if (d_list[target] == -1){
                d_list[target] = d;
                (*ptr_write) = target; // save node number at queue
                ptr_write++;
            }
        }
        ptr_read++;
    }
    free(queue);
    return d_list;
}

void connect_i2j(Network *ntk, int i, int j){
    int num_edges_i = ntk->n_edges[i];
    // work
    ntk->adj_list[i][num_edges_i] = j;
    ntk->n_edges[i]++;
    // instead of using num_edges ++, use +2
    ntk->adj_list[i] = (int*) realloc(ntk->adj_list[i], sizeof(int)*(num_edges_i+2));
}
