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

#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>

// Bilinear interpolation (Se necesita para versión CPU)
#include "bilinear_interpolation.hpp"
#include <math.h> // Para usar ceil y floor

// FluidsGL CPU definitions
#include "fluidsGL_cpu.hpp"

#define TRACE_ADD_FORCES
#define TRACE_ADVECT_VELOCITY
#define TRACE_ADVECT_PARTICLES

//#define DUMPEABLE
#ifdef DUMPEABLE
extern FILE* dump_file;
#endif

// Note that these kernels are designed to work with arbitrary
// domain sizes, not just domains that are multiples of the tile
// size. Therefore, we have extra code that checks to make sure
// a given thread location falls within the domain boundaries in
// both X and Y. Also, the domain is covered by looping over
// multiple elements in the Y direction, while there is a one-to-one
// mapping between threads in X and the tile size in X.
// Nolan Goodnight 9/22/06

// This method adds constant force vectors to the velocity field
// stored in 'v' according to v(x,t+1) = v(x,t) + dt * f.
void addForces_cpu(cData *v, int dx, int dy, int spx, int spy, float fx, float fy, int r)
{
	while(spx < 0) {
                // Transform negative coordinates
                // into possitive ones
		spx = dx + spx;
	}
	while(spy < 0) {
		spy = dy + spy;
	}
	int nx = 2*r+1;
	int ny = 2*r+1;
	int x;
	int y;

	#if defined(TRACE_ADD_FORCES) && defined(DUMPEABLE)
        fprintf(dump_file, "\n''addForces_cpu''''''''''''''''''''\n");
	fprintf(dump_file, "nx = %d\n", nx);
	fprintf(dump_file, "ny = %d\n", ny);
	fprintf(dump_file, "spy = %d\n", spx);
	fprintf(dump_file, "spx = %d\n", spy);
	fprintf(dump_file, "dx = %d\n", dx);
	fprintf(dump_file, "dy = %d\n", dy);
	fprintf(dump_file, "fx = %f\n", fx);
	fprintf(dump_file, "fy = %f\n", fy);
	#endif

	for(x=0; x < nx; x++) {
		for(y=0; y < ny; y++)
		{
			int viy = y + spy;
			int vix = x + spx;
			while(vix > dx-1) {
				// If we have a coordinate 
				// outside canvas
				vix = vix - dx;
			}
			while(viy > dy-1) {
				viy = viy - dy;
			}
			int vi = viy*dx + vix;
			cData vterm = v[vi];
			#if defined(TRACE_ADD_FORCES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tvterm = v[%d] = (%f,%f)", vi, v[vi].x, v[vi].y);
			#endif
	
			int xmr = x-r;
			int ymr = y-r;

			// Compute force smooth as the point is far from the radious center: 1/(1 + x^4 + y^4)
			float s = 1.f / (1.f + xmr*xmr*xmr*xmr + ymr*ymr*ymr*ymr);
			vterm.x += s * fx;
			vterm.y += s * fy;
			v[vi] = vterm;

			#if defined(TRACE_ADD_FORCES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tv[%d] += (%f*%f,%f*%f)", vi, s, fx, s, fy);
			fprintf(dump_file, "\tv[%d] = (%f,%f)\n", vi, vterm.x, vterm.y);
			#endif
		}
		#if defined(TRACE_ADD_FORCES) && defined(DUMPEABLE)
		fprintf(dump_file,"\n");
		#endif
	}

	#if defined(TRACE_ADD_FORCES) && defined(DUMPEABLE)
        fprintf(dump_file, "\n'''''''''''''''''''''''''''''''''''\n");
	#endif
}

// This method performs the velocity advection step, where we
// trace velocity vectors back in time to update each grid cell.
// That is, v(x,t+1) = v(p(x,-dt),t). Here we perform bilinear
// interpolation in the velocity space.
void
advectVelocity_cpu(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy, float dt)
{
    cData vterm, ploc;
    float vxterm, vyterm;

    #if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
    fprintf(dump_file, "\n''advectVelocity_cpu''''''''''''''''''''\n");
    #endif

	int x, y;
	for(x=0; x < dx; x++) {
		for(y=0; y < dy; y++) {
			int f = y*dx + x;
			vterm = v[f];
			ploc.x = x - (dt * vterm.x * dx);
			ploc.y = y - (dt * vterm.y * dy);
			#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
			fprintf(dump_file, "vterm = v[%d] = (%f,%f)\t",f,vterm.x,vterm.y);
			fprintf(dump_file, "ploc = (%f,%f)\t",ploc.x,ploc.y);
			#endif

			// Calculate 4 points to interpolate
			int x1 = floor(ploc.x);
			int y1 = floor(ploc.y);
			// Force xy1 to be different than xy2
                        // don't use ceil as for ploc.xy == .0
                        // they will match
			int x2 = x1+1; // int x2 = ceil(ploc.x)
			int y2 = y1+1; // int y2 = ceil(ploc.y)

			// transform outside coordinates to
                        // coordinates inside our size boundaries
			int x1d = x1 % dx;
			int x2d = x2 % dx;
			int y2d = y2 % dy;
			int y1d = y1 % dy;
			// Force possitive coordinates
			if(x1d < 0) {
				x1d = dx+x1d;
			}
			if(x2d < 0) {
            			x2d = dx+x2d;
			}
            		if(y1d < 0) {
            			y1d = dy+y1d;
			}
			if(y2d < 0) {
				y2d = dy+y2d;
			}

			#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
			fprintf(dump_file, "(x1,y2,x2,y2) = (%d,%d,%d,%d)\t",x1,y1,x2,y2);
			#endif

			int q11f = y1d*dx + x1d;
			int q12f = y2d*dx + x1d;
			int q21f = y1d*dx + x2d;
			int q22f = y2d*dx + x2d;
			cData q11v = v[q11f];
			cData q12v = v[q12f];
			cData q21v = v[q21f];
			cData q22v = v[q22f];
			#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
			fprintf(dump_file, "(q11v,q12v,q21v,q22v) = ((%f,%f),(%f,%f)"
				",(%f,%f),(%f,%f))\t",q11v.x,q11v.y,q12v.x,q12v.y,
				q21v.x,q21v.y,q22v.x,q22v.y);
			#endif
                        // Individually interpolate x and y components
			vterm.x = BilinearInterpolation(
					q11v.x, q12v.x, q21v.x, q22v.x,
					x1,x2,y1,y2,
					ploc.x, ploc.y
				);
			vterm.y = BilinearInterpolation(
					q11v.y, q12v.y, q21v.y, q22v.y,
					x1,x2,y1,y2,
					ploc.x, ploc.y
				);
			vxterm = vterm.x;
			vyterm = vterm.y;
			// The next step will use complex numbers
                        // save the result in vx and vy
			int fj = y*pdx + x;
			vx[fj] = vxterm;
			vy[fj] = vyterm;
			#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
			fprintf(dump_file, "vx[%d] = %f\t",fj,vxterm);
			fprintf(dump_file, "vy[%d] = %f\t",fj,vyterm);
			fprintf(dump_file,"\n");
			#endif
		}
		#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
		fprintf(dump_file,"\n");
		#endif
	}
	#if defined(TRACE_ADVECT_VELOCITY) && defined(DUMPEABLE)
        fprintf(dump_file, "\n'''''''''''''''''''''''''''''''''''\n");
	#endif
}

// This method performs velocity diffusion and forces mass conservation
// in the frequency domain. The inputs 'vx' and 'vy' are complex-valued
// arrays holding the Fourier coefficients of the velocity field in
// X and Y. Diffusion in this space takes a simple form described as:
// v(k,t) = v(k,t) / (1 + visc * dt * k^2), where visc is the viscosity,
// and k is the wavenumber. The projection step forces the Fourier
// velocity vectors to be orthogonal to the vectors for each
// wavenumber: v(k,t) = v(k,t) - ((k dot v(k,t) * k) / k^2.
void
diffuseProject_cpu(cData *vx, cData *vy, int dx, int dy, float dt,
                 float visc)
{
	cData xterm, yterm;

	int x, y;
	for(x=0; x < dx; x++) {
		for(y=0; y < dy; y++) {
			int f = y*dx + x;
			xterm = vx[f];
			yterm = vy[f];

			// Compute the index of the wavenumber based on the
			// data order produced by a standard NN FFT.
			int iix = x;
			int iiy = (y>dy/2)?(y-(dy)):y;

			// Velocity diffusion
			// v(k,t) = v(k,t) / (1 + visc * dt * k^2)
			float kk = (float)(iix * iix + iiy * iiy); // k^2
			float diff = 1.f / (1.f + visc * dt * kk);
			xterm.x *= diff;
			xterm.y *= diff;
			yterm.x *= diff;
			yterm.y *= diff;

			// Velocity projection
			if (kk > 0.f)
			{
			    float rkk = 1.f / kk;
			    // Real portion of velocity projection
			    float rkp = (iix * xterm.x + iiy * yterm.x);
			    // Imaginary portion of velocity projection
			    float ikp = (iix * xterm.y + iiy * yterm.y);
			    xterm.x -= rkk * rkp * iix;
			    xterm.y -= rkk * ikp * iix;
			    yterm.x -= rkk * rkp * iiy;
			    yterm.y -= rkk * ikp * iiy;
			}

			vx[f] = xterm;
			vy[f] = yterm;
		}
	}
}

// This method updates the velocity field 'v' using the two complex
// arrays from the previous step: 'vx' and 'vy'. Here we scale the
// real components by 1/(dx*dy) to account for an unnormalized FFT.
void
updateVelocity_cpu(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy)
{
    float vxterm, vyterm;
    cData nvterm;

	int x, y;
	for(x=0; x < dx; x++) {
		for(y=0; y < dy; y++) {
			int fr = y*pdx + x;
			vxterm = vx[fr];
			vyterm = vy[fr];

			// Normalize the result of the inverse FFT
			float scale = 1.f / (dx * dy);
			nvterm.x = vxterm * scale;
			nvterm.y = vyterm * scale;

			int f = y*dx + x;
			v[f] = nvterm;
		}
	}
}

// This method updates the particles by moving particle positions
// according to the velocity field and time step. That is, for each
// particle: p(t+1) = p(t) + dt * v(p(t)).
void
advectParticles_cpu(cData *part, cData *v, int dx, int dy,
                  float dt)
{
	cData pterm, vterm;

	#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
       	fprintf(dump_file, "\n''advectParticles_cpu''''''''''''''''''''\n");
	fprintf(dump_file, "dt = %f\n", dt);
	fprintf(dump_file, "dx = %d\n", dx);
	fprintf(dump_file, "dy = %d\n", dy);
	#endif

	int x, y;
	for(x=0; x < dx; x++) {
		for(y=0; y < dy; y++) {
			int f = y*dx + x;
			pterm = part[f];
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tpterm = part[%d] = (%f,%f)\n", f, pterm.x, pterm.y);
			#endif

			int xvi = ((int)(pterm.x * dx));
			int yvi = ((int)(pterm.y * dy));
			int vi = yvi*dx + xvi;
			vterm = v[vi];
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tv[%d*%d+%d]=v[%d]=(%f,%f)\n", xvi,dx,yvi,vi,vterm.x,vterm.y);
			#endif

			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tpterm.x += %f", dt*vterm.x);
			#endif
			pterm.x += dt * vterm.x;
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, " -> %f\n", pterm.x);
			fprintf(dump_file, "\t%f - %d", pterm.x, (int)pterm.x);
			#endif
			pterm.x = pterm.x - (int)pterm.x;
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "= %f\n", pterm.x);
			#endif
			pterm.x += 1.f;
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tpterm.x + 1 = %f\n", pterm.x);
			#endif
			pterm.x = pterm.x - (int)pterm.x;
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tpterm.x - (int)pterm.x 1 = %f\n", pterm.x);
			#endif

			pterm.y += dt * vterm.y;
			pterm.y = pterm.y - (int)pterm.y;
			pterm.y += 1.f;
			pterm.y = pterm.y - (int)pterm.y;

			part[f] = pterm;
			#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
			fprintf(dump_file, "\tpart[%d] = (%f,%f)\n", f, part[f].x, part[f].y);
			fprintf(dump_file, "\n");
			#endif
		}
		#if defined(TRACE_ADVECT_PARTICLES) && defined(DUMPEABLE)
		fprintf(dump_file, "\n\n");
		#endif
	}
}

// These are the external function calls necessary for launching fluid simulation
extern "C"
void addForces(cData *v, int dx, int dy, int spx, int spy, float fx, float fy, int r)
{
   addForces_cpu(v, dx, dy, spx, spy, fx, fy, r);
}

extern "C"
void advectVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy, float dt)
{
   advectVelocity_cpu(v, vx, vy, dx, pdx, dy, dt);
}

// Ver: http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
extern "C"
void diffuseProject(cData *vx, cData *vy, int dx, int dy, float dt, float visc)
{
	fftwf_plan plan;

    // Forward FFT
	// fftw_plan fftwf_plan_dft_r2c_1d(int n, float *in, fftwf_complex *out, unsigned flags);
	// typedef float fftwf_complex[2]
	plan = fftwf_plan_dft_r2c_2d(DIM, DIM, (float*) vx, (float (*)[2]) vx, 0);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	plan = fftwf_plan_dft_r2c_2d(DIM, DIM, (float*) vy, (float (*)[2]) vy, 0);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);

    diffuseProject_cpu(vx, vy, dx, dy, dt, visc);

    // Inverse FFT
	plan = fftwf_plan_dft_c2r_2d(DIM, DIM, (float (*)[2]) vx, (float*) vx, 0);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	plan = fftwf_plan_dft_c2r_2d(DIM, DIM, (float (*)[2]) vy, (float*) vy, 0);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
}

extern "C"
void updateVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy)
{
    updateVelocity_cpu(v, vx, vy, dx, pdx, dy);
}

extern "C"
void advectParticles(cData* particles, cData *v, int dx, int dy, float dt)
{
    advectParticles_cpu(particles, v, dx, dy, dt);
}
