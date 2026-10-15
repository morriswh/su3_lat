#ifndef SIZES_H
#define SIZES_H

#define beta 0.3
#define PI 3.141592653589793
#define NC 3
#define T 64 // lattice temporal size
#define L 32 // lattice spacial size
#define D 4 // spacetime dimensions
#define X 18 // SU(3) gauge theory, 12 reals
#define S_X 12 // for short storage
#define CLENGTH T*L*L*L*D*X
#define S_CLENGTH T*L*L*L*D*S_X
#define OLENGTH T*L*L*L*X    
#define S_OLENGTH T*L*L*L*S_X  
#define LLENGTH T*L*L*L*D    
#define LVOLUME T*L*L*L      
#define LCUBE  L*L*L        
#define LSQUARE  L*L          
#define NUM_DATAS 2

#endif // !SIZES_H
