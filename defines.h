/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
 
#ifndef DEFINES_H
#define DEFINES_H

extern int DIM;    // Square size of solver domain
extern int DS;     // Total domain size
extern int CPADW;  // Padded width for real->complex in-place FFT
extern int RPADW;  // Padded width for real->complex in-place FFT
extern int PDS;    // Padded total domain size

extern float DT;          // Delta T for extern interative solver
extern float VIS;         // Viscosity constant
extern float FORCE;       // Force scale factor 
extern int FR;          // Force update radius

extern int TILEX; // Tile width
extern int TILEY; // Tile height
extern int TIDSX; // Tids in X
extern int TIDSY; // Tids in Y

void updateVariables();

// Vector data type used to velocity and force fields
typedef struct cData {
	float x,y;
} cData;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#endif
