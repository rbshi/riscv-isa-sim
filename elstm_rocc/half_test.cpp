#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "half.hpp"
using namespace half_float;

int main(){

    uint16_t aa = 13098;
    half* c = reinterpret_cast<half*>(&aa);
    std::cout << *c << std::endl;

}