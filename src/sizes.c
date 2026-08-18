#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const long int NC = 3;
const long int T = 64; // lattice temporal size
const long int L = 32; // lattice spacial size
const long int D = 4; // spacetime dimensions
const long int X = 18; // SU(3) gauge theory, 18 reals
const long int S_X = 12; // for short storage
const long int CLENGTH   = T*L*L*L*D*X;   //= 150994944; //  
const long int S_CLENGTH = T*L*L*L*D*S_X; //= 100663296; // 
const long int OLENGTH   = T*L*L*L*X;     //= 37748736;  // 
const long int S_OLENGTH = T*L*L*L*S_X;   //= 25165824;  // 
const long int LLENGTH   = T*L*L*L*D;     //= 8388608;   // 
const long int LVOLUME   = T*L*L*L;       //= 2097152;   // 
const long int LCUBE     = L*L*L;         //= 32768;     // 
const long int LSQUARE   = L*L;           //= 1024;      // 
const long int NUM_DATAS = 2;


