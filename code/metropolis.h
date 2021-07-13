#include <iostream>
#include "su3.h"

extern uniform_real_distribution<> rdist;
extern uniform_int_distribution<int> idist;

int update(ranlux48& rnd);
void print();
