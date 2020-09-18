#include <stdio.h>
#include <stdlib.h>
#include "mt64.h"
#include "networks.h"

void test_fn(int a, int b);

int main(){
    Network ntk;
    int N = 10000;
    int sum_km = 0; // mean degree
    
    ntk.N = N;
    InitNetwork(&(ntk), N);
    create_BA_network(&(ntk), 3);

    for (int i=0; i<N; i++){
        sum_km += ntk.n_edges[i];
    }
    printf("mean degree = %d\n", sum_km);

    printf("Done\n");
}

