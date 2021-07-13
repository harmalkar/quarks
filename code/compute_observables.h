#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "su3.h"

#define NMAX 1000000
#define B 100

bool read_configuration();
void bootstrap(double *, double *, double *, int);
