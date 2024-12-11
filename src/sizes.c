#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const int T = 64; // lattice temporal size
const int L = 32; // lattice spacial size
const int D = 4; // spacetime dimensions
const int X = 18; // SU(3) gauge theory, 18 reals
const int S_X = 12; // for short storage
const int CLENGTH   = 150994944; // = T*L*L*L*D*X; // 
const int S_CLENGTH = 100663296; // = T*L*L*L*D*S_X; // short storage 
const int OLENGTH   = 37748736; //  = T*L*L*L*X; //
const int S_OLENGTH = 25165824; //  = T*L*L*L*S_X; //
const int LLENGTH   = 8388608; //   = T*L*L*L*D;
const int LVOLUME   = 2097152; //   = T*L*L*L;
const int LCUBE     = 32768; //     = L*L*L;
const int LSQUARE   = 1024; //      = L*L;


