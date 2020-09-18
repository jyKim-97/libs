#include <stdio.h>
#include <math.h>
#include "mt64.h"

int rand_val_max(int N);

int main(){
    int N=0, Nmax=100;
    int tmp_N;
    for (int n=0; n<100000; n++){
        tmp_N = rand_val_max(Nmax);
        if (tmp_N > N){
            N = tmp_N;
        }
    }
    printf("max: %d\n", N);

}

int rand_val_max(int N){
    // return value for [0, N)
    return (int) floor(genrand64_real2() * N);
}