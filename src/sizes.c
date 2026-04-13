#include <stdint.h>
#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const int NC = 3;
const int T = 64; // lattice temporal size
const int L = 32; // lattice spacial size
const int D = 4; // spacetime dimensions
const int X = 18; // SU(3) gauge theory, 18 reals
const int S_X = 12; // for short storage
const uint64_t CLENGTH   = T*L*L*L*D*X;   // this has to be 64 bits to handle larger lattices  
const int S_CLENGTH = T*L*L*L*D*S_X; // 
const uint64_t OLENGTH   = T*L*L*L*X;     // 
const int S_OLENGTH = T*L*L*L*S_X;   // 
const int LLENGTH   = T*L*L*L*D;     // 
const int LVOLUME   = T*L*L*L;       // 
const int LCUBE     = L*L*L;         // 
const int LSQUARE   = L*L;           // 
const int NUM_DATAS = 2;


