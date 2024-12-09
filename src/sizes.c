#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const int T = 64; // lattice temporal size
const int L = 32; // lattice spacial size
const int D = 4; // spacetime dimensions
const int X = 18; // SU(3) gauge theory, 18 reals
const int S_X = 12; // for short storage
const int CLENGTH = T*L*L*L*D*X; // 
const int S_CLENGTH = T*L*L*L*D*S_X; // short storage 
const int OLENGTH = T*L*L*L*X; //
const int S_OLENGTH = T*L*L*L*S_X; //
const int LLENGTH = T*L*L*L*D;
const int LVOLUME = T*L*L*L;
const int LCUBE   = L*L*L;
const int LSQUARE = L*L;

