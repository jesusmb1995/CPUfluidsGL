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
 
// This is a ported CPU version of the CUDA fluidsGL demo
// originally published by NVidia.
//
// Port made by Jesús Martín Berlanga
//
// Info and documentation at rtfluids.bdevel.org
 
#ifndef __STABLEFLUIDS_CPU_CH_
#define __STABLEFLUIDS_CPU_CH_

#include "defines.h"

// This method adds constant force vectors to the velocity field
// stored in 'v' according to v(x,t+1) = v(x,t) + dt * f.
void
addForces_cpu(cData *v, int dx, int dy, int spx, int spy, float fx,
	 float fy, int r, size_t pitch, int nx, int ny);

// This method performs the velocity advection step, where we
// trace velocity vectors back in time to update each grid cell.
// That is, v(x,t+1) = v(p(x,-dt),t). Here we perform bilinear
// interpolation in the velocity space.
void
advectVelocity_cpu(cData *v, float *vx, float *vy,
	int dx, int pdx, int dy, float dt, int lb);

// This method performs velocity diffusion and forces mass conservation
// in the frequency domain. The inputs 'vx' and 'vy' are complex-valued
// arrays holding the Fourier coefficients of the velocity field in
// X and Y. Diffusion in this space takes a simple form described as:
// v(k,t) = v(k,t) / (1 + visc * dt * k^2), where visc is the viscosity,
// and k is the wavenumber. The projection step forces the Fourier
// velocity vectors to be orthogonal to the wave wave vectors for each
// wavenumber: v(k,t) = v(k,t) - ((k dot v(k,t) * k) / k^2.
void
diffuseProject_cpu(cData *vx, cData *vy, int dx, int dy, float dt,
	float visc, int lb);

// This method updates the velocity field 'v' using the two complex
// arrays from the previous step: 'vx' and 'vy'. Here we scale the
// real components by 1/(dx*dy) to account for an unnormalized FFT.
void
updateVelocity_cpu(cData *v, float *vx, float *vy,
	int dx, int pdx, int dy, int lb);

// This method updates the particles by moving particle positions
// according to the velocity field and time step. That is, for each
// particle: p(t+1) = p(t) + dt * v(p(t)).
void
advectParticles_cpu(cData *part, cData *v, int dx, int dy,
	float dt, int lb);

#endif

