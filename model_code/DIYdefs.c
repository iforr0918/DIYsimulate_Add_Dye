/**
 * DIYdefs.c
 *
 * Andrew Stewart
 *
 * Contains a series of subroutines that are of practical use in all
 * annular simulation programs.
 *
 */

#include "DIYdefs.h"


/**
 * printMatrix
 *
 * Convenience method to write a matrix of data (mat, dimensions m by n)
 * to a data file (outfile).
 *
 */
void printMatrix (FILE * outfile, real ** mat, uint m, uint n)
{
    uint i,j;

    // Write array of data
    for (j = 0; j < n; j ++)
    {
        for (i = 0; i < m; i ++)
        {
            //fprintf(outfile,EXP_FMT,mat[i][j]);
            //fprintf(outfile," ");
            fwrite(&mat[i][j],sizeof(real),1,outfile);
        }
        //fprintf(outfile,"\r\n");
    }
    //fflush(outfile);
}


/**
 * printVector
 *
 * Convenience method to write a vector of data (vec) of length m
 * to a data file (outfile).
 *
 */
void printVector (FILE * outfile, real * vec, uint m)
{
    uint i;
      
    for (i = 0; i < m; i ++)
    {
        fwrite(&vec[i],sizeof(real),1,outfile);
    }
}


/**
 * readVector
 *
 * Reads 'len' real values into the vector 'vec' from the file
 * specified by 'fname'. Any error will result in 'false' being
 * returned, and an error message being printed to the specified
 * FILE/stream 'errstrm'.
 *
 */
bool readVector (char * fname, real * vec, uint len, FILE * errstrm)
{
    // Attempt to open the topography file
    FILE * pfile = NULL;
    int i = 0;
    bool readok = true;
    
    // Basic error checking
    if (fname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readVector'\r\n");
        }
        return false;
    }
    if (vec == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL vector supplied to 'readVector'\r\n");
        }
        return false;
    }
    
    // Attempt to open the file
    pfile = fopen(fname,"r");
    
    // Check that the file exists
    if (pfile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
        }
        return false;
    }
    
    // Read in 'len' real values from the file
    for (i = 0; i < len; i ++)
    {
        //if (fscanf(pfile,FLT_FMT,&vec[i]) == EOF)
        if (fread(&vec[i],sizeof(real),1,pfile) < 1)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",i,len,fname);
            }
            readok = false;
            break;
        }
    }
    
    // Close the file and exit
    fclose(pfile);
    return readok;
}


/**
 * readMatrix
 *
 * Reads 'm'x'n' real values into the matrix 'mat' in column-major
 * order (i.e. reads mat[1][1], mat[2][1], mat[3][1], ... , mat[2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm)
{
    // Attempt to open the topography file
    FILE * pfile = NULL;
    int i = 0;
    int j = 0;
    bool readok = true;
    
    // Basic error checking
    if (fname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix'\r\n");
        }
        return false;
    }
    if (mat == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix'\r\n");
        }
        return false;
    }
    for (i = 0; i < m; i ++)
    {
        if (mat[i] == NULL)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix'\r\n");
            }
            return false;
        }
    }
    
    // Attempt to open the file
    pfile = fopen(fname,"r");
    
    // Check that the file exists
    if (pfile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
        }
        return false;
    }
    
    // Read in 'len' real values from the file
    for (j = 0; j < n; j ++)
    {
        for (i = 0; i < m; i ++)
        {
            //if (fscanf(pfile,FLT_FMT,&mat[i][j]) == EOF)
            if (fread(&mat[i][j],sizeof(real),1,pfile) < 1)
            {
                if (errstrm != NULL)
                {
                    fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",j*m+i,m*n,fname);
                }
                readok = false;
                break;
            }
        }
    }
    
    // Close the file and exit
    fclose(pfile);
    return readok;
}


/**
 * readMatrix3
 *
 * Reads 'Nr'x'Nc'x'Ns' real values into the matrix 'mat3' in column-major
 * order (i.e. reads mat[1][1][1], mat[2][1][1], mat[3][1][1], ... , mat[1][2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm)
{
  FILE * pfile = NULL;
  int i = 0;
  int j = 0;
  int k = 0;
  bool readok = true;
  
  // Basic error checking
  if (fname == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  if (mat3 == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  for (i = 0; i < Nr; i ++)
  {
    if (mat3[i] == NULL)
    {
      if (errstrm != NULL)
      {
        fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
      }
      return false;
    }
    for (j = 0; j < Nc; j ++)
    {
      if (mat3[i][j] == NULL)
      {
        if (errstrm != NULL)
        {
          fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix3'\r\n");
        }
        return false;
      }
    }
  }
  
  // Attempt to open the file
  pfile = fopen(fname,"r");
  
  // Check that the file exists
  if (pfile == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
    }
    return false;
  }
  
  // Read data from the file
  for (k = 0; k < Ns; k ++)
  {
    for (j = 0; j < Nc; j ++)
    {
      for (i = 0; i < Nr; i ++)
      {
        //if (fscanf(pfile,"%le",&mat3[i][j][k]) == EOF)
        if (fread(&mat3[i][j][k],sizeof(real),1,pfile) < 1)
        {
          if (errstrm != NULL)
          {
            fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",k*Nr*Nc+j*Nr+i,Nr*Nc*Ns,fname);
          }
          readok = false;
          break;
        }
      }
    }
  }
  
  // Close the file and exit
  fclose(pfile);
  return readok;
}


/**
 * readMatrix4
 *
 * Reads 'Nr'x'Nc'x'Ns'x'Nt' real values into the matrix 'mat4' in column-major
 * order (i.e. reads mat[1][1][1], mat[2][1][1], mat[3][1][1], ... , mat[1][2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix4 (char * fname, real **** mat4, uint Nr, uint Nc, uint Ns, uint Nt, FILE * errstrm)
{
  FILE * pfile = NULL;
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  bool readok = true;
  
  // Basic error checking
  if (fname == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix4'\r\n");
    }
    return false;
  }
  if (mat4 == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix4'\r\n");
    }
    return false;
  }
  for (i = 0; i < Nr; i ++)
  {
    if (mat4[i] == NULL)
    {
      if (errstrm != NULL)
      {
        fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix4'\r\n");
      }
      return false;
    }
    for (j = 0; j < Nc; j ++)
    {
      if (mat4[i][j] == NULL)
      {
        if (errstrm != NULL)
        {
          fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix4'\r\n");
        }
        return false;
      }
      for (k = 0; k < Ns; k ++)
      {
        if (mat4[i][j][k] == NULL)
        {
          if (errstrm != NULL)
          {
            fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix4'\r\n");
          }
          return false;
        }
      }
    }
  }
  
  // Attempt to open the file
  pfile = fopen(fname,"r");
  
  // Check that the file exists
  if (pfile == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
    }
    return false;
  }
  
  // Read data from the file
  for (l = 0; l < Nt; l ++)
  {
    for (k = 0; k < Ns; k ++)
    {
      for (j = 0; j < Nc; j ++)
      {
        for (i = 0; i < Nr; i ++)
        {
          if (fread(&mat4[i][j][k][l],sizeof(real),1,pfile) < 1)
          {
            if (errstrm != NULL)
            {
              fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",l*Nr*Nc*Ns+k*Nr*Nc+j*Nr+i,Nr*Nc*Ns*Nt,fname);
            }
            readok = false;
            break;
          }
        }
      }
    }
  }
  
  // Close the file and exit
  fclose(pfile);
  return readok;
}



/**
 *  rmatalloc
 *
 *  Convenience method to allocate a real matrix with dimensions
 *  (m,n). Returns a pointer to an array of pointers that point to the
 *  beginning of each row in the matrix.
 *
 */
real ** rmatalloc (uint m, uint n)
{
    real * data = fftw_malloc(m*n*sizeof(real));
    real ** mat = fftw_malloc(m*sizeof(real *));
    uint i;

    if (data == NULL || mat == NULL)
    {
        return NULL;
    }

    mat[0] = data;

    for (i = 1; i < m; i ++)
    {
        mat[i] = mat[i-1] + n;
    }

    return mat;
}


/**
 *  cmatalloc
 *
 *  Convenience method to allocate a complex matrix with dimensions
 *  (m,n). Returns a pointer to an array of pointers that point to the
 *  beginning of each row in the matrix.
 *
 */
complex ** cmatalloc (uint m, uint n)
{
    complex * data = fftw_malloc(m*n*sizeof(complex));
    complex ** mat = fftw_malloc(m*sizeof(complex *));
    uint i;

    if (data == NULL || mat == NULL)
    {
        return NULL;
    }

    mat[0] = data;

    for (i=1;i<m;i++)
    {
        mat[i] = mat[i-1] + n;
    }

    return mat;
}


/**
 *  rmatfree
 *
 *  Frees memory used by a matrix allocated with rmatalloc.
 *
 */
void rmatfree (real ** mat)
{
    fftw_free(*mat);
    fftw_free(mat);
}


/**
 *  cmatfree
 *
 *  Frees memory used by a matrix allocated with cmatalloc.
 *
 */
void cmatfree (complex ** mat)
{
    fftw_free(*mat);
    fftw_free(mat);
}


/**
 *  vec2mat
 *
 * Creates a matrix from a vector of length m*n. If the matrix pointed
 * to by pmat (*pmat) is NULL, a vector of pointers is created that must
 * later be freed.
 */
void vec2mat (real * vec, real *** pmat, uint m, uint n)
{
    int i;
    
    if (*pmat == NULL)
    {
        *pmat = fftw_malloc(m*sizeof(real *));
    }
    
    (*pmat)[0] = vec;
    
    for (i = 1; i < m; i ++)
    {
        (*pmat)[i] = (*pmat)[i-1] + n;
    }
}


/**
 * normalise
 *
 * Convenience method to divide all data in the specified array (mat)
 * of size (size) by an integer normaliser (norm).
 *
 */
void normalise (real * mat, uint size, uint norm)
{
    uint i;

    for (i = 0; i < size; i ++)
    {
        mat[i] /= norm;
    }
}


/**
 * rthomas
 *
 * Performs the thomas algorithm. Here a, b and c are the sub-diagonal, diagonal
 * and super-diagonal in the tridiagonal coefficient matrix, and d is the vector
 * on the right hand side. The solution is filled into x.
 * Warning: will modify c and d!
 */
void rthomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n)
{
    real id;
    int i;

	// Modify the coefficients
	c[0] /= b[0];	// Division by zero risk
	d[0] /= b[0];	// Division by zero would imply a singular matrix
	for (i = 1; i < n; i ++)
    {
		id = 1 / (b[i] - c[i-1]*a[i]);     // Division by zero risk
		c[i] *= id;	                       // Last value calculated is redundant
		d[i] = (d[i] - d[i-1]*a[i]) * id;
	}

	// Now back substitute
	x[n - 1] = d[n - 1];
	for (i = n-2; i >= 0; i --)
    {
		x[i] = d[i] - c[i] * x[i + 1];
    }
}


/**
 * cthomas
 *
 * Performs the thomas algorithm. Here a, b and c are the sub-diagonal, diagonal
 * and super-diagonal in the tridiagonal coefficient matrix, and d is the vector
 * on the right hand side. The solution is filled into x.
 * Warning: will modify c and d!
 */
void cthomas (const real * a, const real * b, real * c, complex * d, complex * x, unsigned int n)
{
    real id;
    int i;

	// Modify the coefficients
	c[0] /= b[0];
	d[0][0] /= b[0];	// Division by zero risk - this would imply a singular matrix
	d[0][1] /= b[0];
	for (i = 1; i < n; i ++)
    {
		id = 1 / (b[i] - c[i-1]*a[i]);     // Division by zero risk
		c[i] *= id;	                       // Last value calculated is redundant
		d[i][0] = (d[i][0] - d[i-1][0]*a[i]) * id;
		d[i][1] = (d[i][1] - d[i-1][1]*a[i]) * id;
	}

	// Now back substitute
	x[n-1][0] = d[n-1][0];
	x[n-1][1] = d[n-1][1];
	for (i = n-2; i >= 0; i --)
    {
		x[i][0] = d[i][0] - c[i] * x[i+1][0];
		x[i][1] = d[i][1] - c[i] * x[i+1][1];
    }
}


/**
 * vanderCorput
 *
 * Returns the nth value in the binary van der Corput sequence.
 *
 */
real vanderCorput (const uint n)
{
    return bitReverse(n)/pow(2,32);
}


/**
 * bitReverse
 *
 * Reverses the bits in the specified unsigned int.
 *
 */
uint bitReverse (const uint n)
{
    uint bits = 32;
    uint pos = 0;
    uint out = 0;
  
    // For each bit
    for (pos = 0; pos < bits; pos ++)
    {
        // If the bit is set...
        if ((n & (1<<pos)) != 0)
        {
            // ...then set the corresponding bit in the reverse
            out |= 1 << (bits-1-pos);
        }
    }
    
    return out;
}


/**
 * setParam
 *
 * Convenience method to set the values in the paramdata struct at the specified
 * index (idx) in an array (params). The method sets the parameter name (name)
 * parameter type (type, corresponding to the format string required for fscanf),
 * the pointer to the variable that will hold the parameter value (pvar) and
 * whether the parameter should be marked as having been read (read).
 *
 */
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read)
{
    params[idx].name = name;
    params[idx].type = type;
    params[idx].pvar = pvar;
    params[idx].read = read;
}



/**
 * readParams
 *
 * Reads data into the specified array of paramdata structures
 * 'params' of length 'nparams' from the file specified by
 * 'infname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 */
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm)
{
    FILE * infile = NULL;
    bool readok = true;
    char inbuf[256];
    int i = 0;

    // Basic error checking
    if (infname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readParams'\r\n");
        }
        return false;
    }
    if (params == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL parameter information supplied to 'readParams'\r\n");
        }
        return false;
    }

    // Open the input file and check that it can be read
    infile = fopen(infname,"r");
    if (infile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not find input parameter file\r\n");
        }
        return false;
    }

    // Scan input file for parameter values
    while (fscanf(infile,"%s",inbuf) != EOF)
    {
        // Check whether the string read matches any parameter names
        for (i = 0; i < nparams; i ++)
        {
            if (strcmp(inbuf,params[i].name) == 0)
            {
                if (fscanf(infile,params[i].type,params[i].pvar) == EOF)
                {
                    if (errstrm != NULL)
                    {
                        fprintf(errstrm,"ERROR: Could not read parameter %s\r\n",inbuf);
                    }
                    readok = false;
                }
                else
                {
                    params[i].read = true;
                }
                break;
            }
        }


        // Parameter name doesn't match any of those in the list
        if (i == nparams)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Parameter not recognised: %s\r\n",inbuf);
            }
            readok = false;
        }
    }


    // Check that required parameter values have been specified
    for (i = 0; i < nparams; i ++)
    {
        if (!params[i].read)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Required parameter not specified: %s\r\n",params[i].name);
            }
            readok = false;
        }
    }

    // We don't need the input file any more
    fclose(infile);
    return readok;
}



/**
 * kahanSum
 *
 * Performs Kahan summation over an input vector 'vec' of length 'N'.
 *
 */
real kahanSum (real * vec, int N)
{
  real sum = 0.0;
  real c = 0.0;                  // A running compensation for lost low-order bits.
  int i;
  real y;
  real t;
  
  for (i = 0; i < N; i ++)
  {
    y = vec[i] - c;     // So far, so good: c is zero.
    t = sum + y;          // Alas, sum is big, y small, so low-order digits of y are lost.
    c = (t - sum) - y;    // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
    sum = t;              // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    // Next time around, the lost low part will be added to y in a fresh attempt.
  }
  
  return sum;
}

