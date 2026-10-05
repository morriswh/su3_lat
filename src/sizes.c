#include <inttypes.h>
#include <stddef.h>
#include "sizes.h"

const double beta = 0.3; // 1/g^2
const double PI = 3.141592653589793;
const int64_t NC = 3;
const int64_t T = 64; // lattice temporal size
const int64_t L = 32; // lattice spacial size
const int64_t D = 4; // spacetime dimensions
const int32_t X = 18; // SU(3) gauge theory, 18 reals
const int64_t S_X = 12; // for short storage
const int64_t CLENGTH   = T*L*L*L*D*X;   //= 150994944; //  
const int64_t S_CLENGTH = T*L*L*L*D*S_X; //= 100663296; // 
const int64_t OLENGTH   = T*L*L*L*X;     //= 37748736;  // 
const int64_t S_OLENGTH = T*L*L*L*S_X;   //= 25165824;  // 
const int64_t LLENGTH   = T*L*L*L*D;     //= 8388608;   // 
const int64_t LVOLUME   = T*L*L*L;       //= 2097152;   // 
const int64_t LCUBE     = L*L*L;         //= 32768;     // 
const int64_t LSQUARE   = L*L;           //= 1024;      // 
const size_t NUM_DATAS = 2;


