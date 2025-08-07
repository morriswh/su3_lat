#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const int NC = 3;
const int T = 64; // lattice temporal size
const int L = 32; // lattice spacial size
const int D = 4; // spacetime dimensions
const int X = 18; // SU(3) gauge theory, 18 reals
const int S_X = 12; // for short storage
const int CLENGTH   = T*L*L*L*D*X;   //= 150994944; //  
const int S_CLENGTH = T*L*L*L*D*S_X; //= 100663296; // 
const int OLENGTH   = T*L*L*L*X;     //= 37748736;  // 
const int S_OLENGTH = T*L*L*L*S_X;   //= 25165824;  // 
const int LLENGTH   = T*L*L*L*D;     //= 8388608;   // 
const int LVOLUME   = T*L*L*L;       //= 2097152;   // 
const int LCUBE     = L*L*L;         //= 32768;     // 
const int LSQUARE   = L*L;           //= 1024;      // 
const int NUM_DATAS = 6;


