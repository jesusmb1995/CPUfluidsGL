#include "defines.h"

int DIM   = 512;            // Square size of solver domain
int DS    = (DIM*DIM);      // Total domain size
int CPADW = (DIM/2+1);      // Padded width for real->complex in-place FFT
int RPADW = (2*(DIM/2+1));  // Padded width for real->complex in-place FFT
int PDS   = (DIM*CPADW);    // Padded total domain size

float DT    =  0.09f;          // Delta T for interative solver
float VIS   =  0.0025f;        // Viscosity constant
float FORCE =  (5.8f*DIM);     // Force scale factor 
int FR    =  4;              // Force update radius

void updateVariables() {
	DS    = (DIM*DIM);
	CPADW = (DIM/2+1);
	RPADW = (2*(DIM/2+1));
	PDS   = (DIM*CPADW);
	FORCE =  (5.8f*DIM);
}
