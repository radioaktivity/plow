#include "algos.h"
float get_min_from_array(float * input, int n){

  float min = input[0];

    for(int i = 0; i < n; i++){
      if(input[i] < min){
        min = input[i];
      }
    }
  return min;
}
