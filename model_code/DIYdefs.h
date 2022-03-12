// Header file for QGDefs.c
#ifndef _DIYDEFS_H_
#define _DIYDEFS_H_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>

#include "nsode.h"

// Max name length for files containing input parameter arrays
#define MAX_PARAMETER_FILENAME_LENGTH 256

// Time-integration method identifiers
#define METHOD_RK4 0
#define METHOD_AB3 1

// Azimuthal finite-difference method identifiers
#define METHOD_PS 0
#define METHOD_FD 1

// Handy constants
#define _PI  3.14159265358979323846
#define _2PI 6.28318530717958647692

// Efficient way to square and cube numbers
#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// Handy macros to operate on complex types
#define CMULT(r,x,y)                    \
    r[0] = x[0]*y[0] - x[1]*y[1];       \
    r[1] = x[1]*y[0] + x[0]*y[1];
#define CSUB(r,x,y)                     \
    r[0] = x[0]-y[0];                   \
    r[1] = x[1]-y[1];
#define CADD(r,x,y)                     \
    r[0] = x[0]+y[0];                   \
    r[1] = x[1]+y[1];

typedef fftw_complex complex;

// Input parameter data
typedef struct paramdata
{
    char * name;
    char * type;
    void * pvar;
    bool read;
}
paramdata;

// Stores the position of a tracer particle
typedef struct tracer
{
    real r; // Radial position
    real a; // Azimuthal position
}
tracer;

// To define input parameters
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read);
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm);

// Data I/O
void printMatrix (FILE * outfile, real ** mat, uint m, uint n);
void printVector (FILE * outfile, real * vec, uint m);
bool readMatrix4 (char * fname, real **** mat4, uint Nr, uint Nc, uint Ns, uint Nt, FILE * errstrm);
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm);
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm);
bool readVector (char * fname, real * vec, uint len, FILE * errstrm);

// Wrappers for allocating/freeing FFTW arrays
real ** rmatalloc (uint m, uint n);
complex ** cmatalloc (uint m, uint n);
void rmatfree (real ** mat);
void cmatfree (complex ** mat);
void vec2mat (real * vec, real *** pmat, uint m, uint n);

// Simple convenience methods
void normalise (real * mat, uint size, uint norm);

// Thomas algorithm
void rthomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n);
void cthomas (const real * a, const real * b, real * c, complex * d, complex * x, unsigned int n);

// van der Corput generator
real vanderCorput (const uint n);
uint bitReverse (const uint n);

// Kahan summation
real kahanSum (real * vec, int N);

#endif
