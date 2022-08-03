/**
 * DIYsimulate.c
 *
 * Andrew Stewart
 *
 * 2D simulation of quasigeostrophic flow in a regular annulus.
 *
 * TODO maybe we should be using adaptive time stepping to make the code more robust?
 * TODO a smagorinsky viscosity might be a good call for simulations not able to resolve the molecular viscosity
 * TODO no-slip or no-stress boundary conditions?
 * TODO add temperature variable following Warneford and Dellar 2013?
 * Add Dye Tracers - Isaac Forrest
 */
#include <time.h>

#include "rk.h"
#include "ab.h"
#include "DIYdefs.h"

#define NPARAMS 20

#define RMATALLOC(mat,m,n)                                      \
    mat = rmatalloc(m,n);                                       \
    if (mat == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\n"); \
        return 0;                                               \
    }

#define CMATALLOC(mat,m,n)                                      \
    mat = cmatalloc(m,n);                                       \
    if (mat == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\n"); \
        return 0;                                               \
    }

#define CHECKALLOC(ptr,size)                                    \
    ptr = fftw_malloc(size);                                    \
    if (ptr == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\n"); \
        return 0;                                               \
    }


// DFT plans
fftw_plan forward_transform = NULL;
fftw_plan backward_transform = NULL;

// Data arrays for potential vorticity, vorticity and streamfunction,
// with spectral representations
real ** qvort_r = NULL;
real ** vort_r = NULL;
real ** psi_r = NULL;
real ** redred_r = NULL;

complex ** qvort_s = NULL;
complex ** vort_s = NULL;
complex ** psi_s = NULL;
complex ** redred_s = NULL;

// Used for calculation of advective terms
real ** vq_r = NULL;
complex ** vq_s = NULL;

real ** vc_Red_r = NULL;
complex ** vc_Red_s = NULL;

//

// Real and spectral derivatives with respect to theta
real ** da_vq_r = NULL;
real ** da_vc_Red_r = NULL;
real ** da_psi_r = NULL;
real ** dr_psi_r = NULL;
complex ** da_vq_s = NULL;
complex ** da_vc_Red_s = NULL;
complex ** da_psi_s = NULL;
complex ** dada_vort_s = NULL;
complex ** dada_redred_s = NULL;
real ** dada_redred_r = NULL;
real ** dada_vort_r = NULL;



// Used for boundary circulation evolution
real * u_bdy = NULL;
real * dt_Gamma_wrk = NULL;

// Time derivative of real vorticity - for use in time
// derivative function
real ** dt_qvort_r = NULL;
real ** dt_redred_r = NULL;

// Coefficient arrays used in the update_psi method
real ** psicoeffb = NULL;
real * psicoeffa = NULL;
real * psicoeffc = NULL;
real * psitempc = NULL;
complex * psitempd = NULL;
complex * psitempx = NULL;

// Physical parameters
real rmin = 0;    // Inner and outer radii
real rmax = 0;
real Lr = 0;
real La = 0;
real f = 0.0;     // Coriolis parameter
real H = 1.0;     // Max. water depth
real nu = 0.0;// Viscosity coefficient
real nu_rgb = 0.0;
real ** hh = NULL; // Height of the bottom topography

// To avoid repeated calculation
real _drsq = 0;
real _2dr = 0;
real _fH = 0;
real da = 0;
real _dasq = 0;
real _2da = 0;
real _dr_da = 0;
real _4drda = 0;
real _dr_dasq = 0;

// Dissipation coefficient (i.e. frictional drag)
real kap0 = 0.0;
real ** kap = NULL;

// Quantities related to determining the streamfunction on the outer boundary
real ** qvortL = NULL;
real psiL0 = 0;
real ** psiL = NULL;
real GammaL = 0;
real ** psiP = NULL;

// Radial grid points
real * rr = NULL;
real * rrsq = NULL;

// Array of wave numbers
int * kk = NULL;
int * kksq = NULL;

// Spectral Fourier smoothing filter
const int fil_m = 36;
const int fil_a = 36;
real * filter = NULL;

// Grid size
uint Nr = 0;
uint Na = 0;
uint N = 0;
uint Np = 0;

// Fraction of the annulus under consideration is (1/amult).
int amult = 1;
int amultsq = 1;

// Radial grid point spacing
real dr = 0;

// Time-stepping and azimuthal differencing methods to use
const uint method_t = METHOD_AB3;
const uint method_a = METHOD_FD;

// Output filenames
static const char OUTN_PSI[] = "PSI";
static const char OUTN_PV[] = "PV";
static const char OUTN_TRACER[] = "TRACER";

static const char OUTN_TFILE[] = "time.txt";
//Isaac
static const char OUTN_RED[] = "RED";
//End Isaac


/**
 * calcGamma
 *
 * Calculates the circulation around the inner wall of the annulus.
 *
 */
real calcGamma (real ** psi)
{
    int j = 0;
    real Gamma = 0;
    real dr_psi = 0;
    
    for (j = 0; j < Na; j ++)
    {
        u_bdy[j] = (-psi[2][j]+4*psi[1][j]-3*psi[0][j]) / _2dr;
    }
  
    Gamma = kahanSum(u_bdy,Na) * rmin * da;
    
    return Gamma;
}



/**
 * invertPsi
 *
 * Inverts the specified array of potential vorticities to obtain the
 * streamfunction. Fills in the vort_r and vort_s matrices with the
 * corresponding relative vorticities, an then calculates psi in spectral
 * space, filling it in to psi_s. If keepSpectral is true, it is left in
 * this form, otherwise it is transformed back to real space and filled into psi_r.
 *
 */
void invertPsi (real ** qvort, real psi0, bool keepSpectral)
{
    // Looping variables
    int i, j;

    // Calculate 'true' vorticity from potential vorticity
    for (i = 0; i < Nr; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            vort_r[i][j] = qvort[i][j] - _fH*hh[i][j];
        }
    }

    // Transform to find spectral vorticity vort_r->vort_s
    fftw_execute_dft_r2c(forward_transform,*vort_r,*vort_s);

    // For each wave number, we calculate psi_s by approximating the
    // radial derivatives with central differences and solving the
    // resulting tridiagonal linear system with the Thomas algorithm.
    for (j = 0; j < Na; j ++)
    {
        // First set the boundary values - psi=0 on r=rmin
        // and psi=const on r=rmax. Note the multiplication
        // by Na because the Fourier transforms are unnormalised
        // until we return to real space.
        psi_s[0][j][0] = 0;
        psi_s[0][j][1] = 0;
        psi_s[Nr-1][j][1] = 0;
        if (j == 0)
        {
            psi_s[Nr-1][j][0] = psi0*Na;
        }
        else
        {
            psi_s[Nr-1][j][0] = 0;
        }

        // Set up 'c' and 'd' arrays for the Thomas algorithm
        memcpy(psitempc,psicoeffc,Nr*sizeof(real));
        for (i = 1; i < Nr-1; i ++)
        {
            psitempd[i][0] = _drsq*vort_s[i][j][0];
            psitempd[i][1] = _drsq*vort_s[i][j][1];
        }

        // Here we account for contributions of the boundary psi-values to the
        // RHS of the tridiagonal matrix for the Thomas algorithm. We need only
        // consider wavenumber zero because psi=constant implies that only wavenumber
        // zero is nonzero on the boundary. We don't need to worry about the boundary
        // at r=rmin because we impose psi=0 there always.
        if (j == 0)
        {
            psitempd[Nr-2][0] -= (1+dr/(2*rr[Nr-2]))*psi_s[Nr-1][0][0];
        }

        // Perform the Thomas algorithm for this wave number
        // NOTE: We supply each array at starting index 1, because only (Nr-2)
        // psi-values need to be solved for
        cthomas(psicoeffa+1,psicoeffb[j]+1,psitempc+1,psitempd+1,psitempx+1,Nr-2);

        // Copy calculated values to psi_s
        // This could be done more efficiently by incorporating the Thomas
        // algorithm directly, but this is clearer and doesn't waste much time
        for (i = 1; i < Nr-1; i ++)
        {
            psi_s[i][j][0] = psitempx[i][0];
            psi_s[i][j][1] = psitempx[i][1];
        }
    }

    // If required, transform back to real space
    if (!keepSpectral)
    {
        fftw_execute_dft_c2r(backward_transform,*psi_s,*psi_r);
        normalise(*psi_r,Nr*Na,Na);
    }
}



/**
 * calcPsi
 *
 * Calculates the streamfunction at any given time, given the PV
 * everywhere and the circulation around the inner wall.
 *
 */
void calcPsi (real ** qvort, real Gamma, real ** psi)
{
    // Laplacian coefficient. The streamfunction is constructed as a
    // sum of solutions
    //    psi = psiP + eta*psiL,
    // where psiP and psiL are respectively solutions to the Poisson equation
    //    Del psiP = omega, psiP = 0 on r=min and r=rmax,
    // and the Laplace equation
    //    Del psiL = 0, psiL = 0 on r=rmin, psiL=psiL0 on r=rmax.
    real eta = 0;
    real GammaP = 0;
    int i = 0;
    int j = 0;
    
    // Calculate psiP, the solution to the Poisson
    // equation with psi=0 on both boundaries
    invertPsi(qvort,0,false);
    memcpy(*psiP,*psi_r,Nr*Na*sizeof(real));
    
    // Calculate GammaP, the inner-wall circulation of the Poisson streamfunction
    GammaP = calcGamma(psiP);
    
    // Calculate eta such that the sum GammaP+eta*GammaL constructs the correct
    // circulation around the inner wall
    eta = (Gamma-GammaP)/GammaL;
    
    // Construct the streamfunction as psi = psiP + eta*psiL
    for (i = 0; i < Nr; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            psi[i][j] = psiP[i][j] + eta*psiL[i][j];
        }
    }
}



/**
 * tderiv
 *
 * Calculates the time derivative of the potential vorticity field (qvort)
 * by evaluating the terms in the QG vorticity equation at the
 * current time (t). The calculated time derivative is put in the
 * 'dt_vort' output array. The parameter 'numvars' specifies the
 * length of the real arrays 'qvort' and 'dt_qvort'.
 *
 * This function has been formatted as a DERIVATIVE_FUNCTION, so
 * that is can be supplied to the general-purpose ODE solvers.
 * In fact, the real arrays 'qvort' and 'dt_qvort' will contain
 * real matrices of size Nr by Na. The parameters 't' and
 * 'numvars' are, in fact, unused, because the format of the
 * data passed in via 'qvort' is known.
 *
 * NOTE: We are using INPUT-DESTROYING complex-to-real transforms,
 * so complex spectral arrays MUST NOT be used in this function after
 * they have been transformed back to real space.
 *
 */
void tderiv (const real t, const real * vars, real * dt_vars, const uint numvars)
{
    // Looping variables
    int i, j, p;

    // Indexing variables for azimuthal finite-differencing
    int jp1, jm1;

    // Temporary storage variables
    real dr_vort = 0;
    real dr_qvort = 0;
    real drdr_vort = 0;
    real drdr_qvort = 0;
    real da_vq = 0;
    real dr_ruq = 0;
    real dada_vort = 0;
    real dada_qvort = 0;
    
    real da_vc_Red = 0;
    real dr_ruc_Red = 0;
    real dr_redred = 0;
    real drdr_redred = 0;
    real dada_redred = 0;
    
    // Pointer to circulation data
    real * Gamma_r = 0;
    real * dt_Gamma_r = NULL;
    
    // Pointers to tracer data
    tracer * dyeline_r = NULL;
    tracer * dt_dyeline_r = NULL;
    
    //Pointers to Red Dye data
    // Pointer to circulation data
    //real * Red_r = 0;
    //real * dt_Red_r = NULL;
    
    // Variables for tracer tendency calculation
    real rpos = 0;
    real apos = 0;
    int iprev = 0;
    int inext = 0;
    int jprev = 0;
    int jnext = 0;
    real rprev = 0;
    real rnext = 0;
    real aprev = 0;
    real anext = 0;
    real wrprev = 0;
    real wrnext = 0;
    real waprev = 0;
    real wanext = 0;
    real dr_psi_tp = 0;
    real da_psi_tp = 0;
  
    // Copy the vorticity data into the transform matrix
    memcpy(*qvort_r,vars+1,Nr*Na*sizeof(real));
  
    // Copy the red dye data into the transform matrix
    memcpy(*redred_r,vars+1+Nr*Na+ Np*sizeof(tracer),Nr*Na*sizeof(real));
    
    // Make sure time derivatives are zero before we set them
    memset(dt_vars,0,numvars*sizeof(real));
    
    // Create pointers for the various subsections of the input data
    Gamma_r = (real *) vars;
    dt_Gamma_r = (real *) dt_vars;
    dyeline_r = (tracer *) (vars + 1 + Nr*Na);
    dt_dyeline_r = (tracer *) (dt_vars + 1 + Nr*Na);
    //Red_r = (real *) (vars + 1 + Nr*Na + Np*sizeof(tracer)) ;
    //dt_Red_r = (real *) (dt_vars + 1 + Nr*Na + Np*sizeof(tracer));

    // CALCULATION OF THE STREAMFUNCTION


    calcPsi(qvort_r,*Gamma_r,psi_r);


    // CALCULATE STREAMFUNCTION THETA-DERIVATIVES
    

    if (method_a == METHOD_PS)
    {
        // Compute psi in spectral space (N.B. this destroys psi_r)
        fftw_execute_dft_r2c(forward_transform,*psi_r,*psi_s);
      
        // Calculate theta-derivatives in spectral space
        for (i = 0; i < Nr; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                da_psi_s[i][j][0] = -amult*kk[j]*psi_s[i][j][1];
                da_psi_s[i][j][1] = amult*kk[j]*psi_s[i][j][0];
            }
        }

        // Transform back to get theta-derivatives in real space
        fftw_execute_dft_c2r(backward_transform,*da_psi_s,*da_psi_r);
        normalise(*da_psi_r,Nr*Na,Na);
      
        // Now transform back to get psi_r (N.B. this destroys psi_s)
        fftw_execute_dft_c2r(backward_transform,*psi_s,*psi_r);
        normalise(*psi_r,Nr*Na,Na);
    }

    // Now we can get the finite-difference theta-derivative using psi_r
    if (method_a == METHOD_FD)
    {
        for (j = 0; j < Na; j ++)
        {
            jp1 = (j+1) % Na;
            jm1 = (j+Na-1) % Na;

            // At each gridpoint...
            for (i = 1; i < Nr-1; i ++)
            {
                 da_psi_r[i][j] = (psi_r[i][jp1]-psi_r[i][jm1]) / _2da;
            }
        }
    }
     
    // Streamfunction is constant along boundaries
    for (j = 0; j < Na; j ++)
    {
        da_psi_r[0][j] = 0;
        da_psi_r[Nr-1][j] = 0;
    }
  
    // Radial derivative of psi
    for (j = 0; j < Na; j ++)
    {
        for (i = 1; i < Nr-1; i ++)
        {
            dr_psi_r[i][j] = (psi_r[i+1][j]-psi_r[i-1][j]) / _2dr;
        }

        // Radial derivatives calculated at first order on boundaries.
        // If these are used, it's only in advecting passive tracers.
        dr_psi_r[0][j] = (psi_r[1][j]-psi_r[0][j]) / dr;
        dr_psi_r[Nr-1][j] = (psi_r[Nr-1][j]-psi_r[Nr-2][j]) / dr;
    }


    // CALCULATION OF THETA DERIVATIVES


    // Calculate q*v in real space
    for (i = 1; i < Nr-1; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            vq_r[i][j] = qvort_r[i][j]*(psi_r[i+1][j]-psi_r[i-1][j])/_2dr;
        }
    }
    for (i = 1; i < Nr-1; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            vc_Red_r[i][j] = redred_r[i][j]*(psi_r[i+1][j]-psi_r[i-1][j])/_2dr;
        }
    }
    

    // If we're using a pseudospectral algorithm, calculate other
    // theta derivatives now
    if (method_a == METHOD_PS)
    {
        // Transform qvort_r -> qvort_s and vq_r -> vq_s
        fftw_execute_dft_r2c(forward_transform,*vq_r,*vq_s);
        fftw_execute_dft_r2c(forward_transform,*qvort_r,*qvort_s);
        
        // Transform redreds_r -> redred_s and vc_Red_r -> vc_Red_s
        fftw_execute_dft_r2c(forward_transform,*vc_Red_r,*vc_Red_s);
        fftw_execute_dft_r2c(forward_transform,*redred_r,*redred_s);
        
        // Calculate theta-derivatives in spectral space:
        // da_psi_s, da_vort_s and dada_vort_s
        for (i = 0; i < Nr; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                da_vq_s[i][j][0] = -amult*kk[j]*vq_s[i][j][1];
                da_vq_s[i][j][1] = amult*kk[j]*vq_s[i][j][0];

                if (nu != 0)
                {
                    dada_vort_s[i][j][0] = -amultsq*kksq[j]*vort_s[i][j][0];
                    dada_vort_s[i][j][1] = -amultsq*kksq[j]*vort_s[i][j][1];
                }
            }
        }
        
        for (i = 0; i < Nr; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                da_vc_Red_s[i][j][0] = -amult*kk[j]*vc_Red_s[i][j][1];
                da_vc_Red_s[i][j][1] = amult*kk[j]*vc_Red_s[i][j][0];

                if (nu_rgb != 0)
                {
                    dada_redred_s[i][j][0] = -amultsq*kksq[j]*redred_s[i][j][0];
                    dada_redred_s[i][j][1] = -amultsq*kksq[j]*redred_s[i][j][1];
                }
            }
        }

        // Transform back to get da_vq_r, dada_vort_r
        fftw_execute_dft_c2r(backward_transform,*da_vq_s,*da_vq_r);
        fftw_execute_dft_c2r(backward_transform,*da_vc_Red_s,*da_vc_Red_r);
        normalise(*da_vq_r,Nr*Na,Na);
        normalise(*da_vc_Red_r,Nr*Na,Na);
        
        if (nu != 0)
        {
            fftw_execute_dft_c2r(backward_transform,*dada_vort_s,*dada_vort_r);
            normalise(*dada_vort_r,Nr*Na,Na);
        }
        if (nu_rgb != 0)
        {
            fftw_execute_dft_c2r(backward_transform,*dada_redred_s,*dada_redred_r);
            normalise(*dada_redred_r,Nr*Na,Na);
        }
    }


    // CALCULATION OF TIME DERIVATIVE
    // This section is coded for clarity rather than efficiency


    // For each radial row of gridpoints...
    for (j = 0; j < Na; j ++)
    {
        jp1 = (j+1) % Na;
        jm1 = (j+Na-1) % Na;

        // At each gridpoint...
        for (i = 1; i < Nr-1; i ++)
        {
            // Pseudospectral method
            if (method_a == METHOD_PS)
            {
                // Theta-derivatives calculated via FFT
                da_vq = da_vq_r[i][j];
                
                da_vc_Red = da_vc_Red_r[i][j];
                
                if (nu != 0)
                {
                    dada_vort = dada_vort_r[i][j];
                }
                
                if (nu_rgb != 0)
                {
                    dada_redred = dada_redred_r[i][j];
                }
            }

            // Finite-differencing method
            else if (method_a == METHOD_FD)
            {
                // Calculate theta-derivatives using second order central differencing
                da_vq = (vq_r[i][jp1]-vq_r[i][jm1]) / _2da;
                
                da_vc_Red = (vc_Red_r[i][jp1]-vc_Red_r[i][jm1]) / _2da;
                
                if (nu != 0)
                {
                    dada_vort = (vort_r[i][jp1]-2*vort_r[i][j]+vort_r[i][jm1]) / _dasq;
                }
                
                if (nu_rgb != 0)
                {
                    dada_redred = (redred_r[i][jp1]-2*redred_r[i][j]+redred_r[i][jm1]) / _dasq;
                }
            }

            // r-derivatives
            dr_ruq = (qvort_r[i-1][j]*da_psi_r[i-1][j]-qvort_r[i+1][j]*da_psi_r[i+1][j]) / _2dr;
            
            dr_ruc_Red = (redred_r[i-1][j]*da_psi_r[i-1][j]-redred_r[i+1][j]*da_psi_r[i+1][j]) / _2dr;

            if (nu != 0)
            {
                dr_vort = (vort_r[i+1][j]-vort_r[i-1][j]) / _2dr;
                drdr_vort = (vort_r[i+1][j]-2*vort_r[i][j]+vort_r[i-1][j]) / _drsq;
            }
            
            if (nu_rgb != 0)
            {
                dr_redred = (redred_r[i+1][j]-redred_r[i-1][j]) / _2dr;
                drdr_redred = (redred_r[i+1][j]-2*redred_r[i][j]+redred_r[i-1][j]) / _drsq;
            }

            // Calculate d(qvort)/dt
            dt_qvort_r[i][j] = - (dr_ruq + da_vq) / rr[i];

            if (kap0 != 0)
            {
                dt_qvort_r[i][j] -= kap[i][j]*vort_r[i][j];
            }
            if (nu != 0)
            {
                // NOTE: This calculation of the r-derivatives is q_rr + (1/r)q_r,
                // the discrete formulation of which is identical to that of (1/r)*(r*q_r)_r
                dt_qvort_r[i][j] += nu*(drdr_vort+dr_vort/rr[i]+dada_vort/rrsq[i]);                
            }
            
            // Calculate d(c_Red)/dt
            dt_redred_r[i][j] = - (dr_ruc_Red + da_vc_Red) / rr[i];

            //if (kap0 != 0)
            //{
            //    dt_qvort_r[i][j] -= kap[i][j]*vort_r[i][j];
            //}
            if (nu_rgb != 0)
            {
                //NOTE: This calculation of the r-derivatives is q_rr + (1/r)q_r,
                //the discrete formulation of which is identical to that of (1/r)*(r*q_r)_r
                dt_redred_r[i][j] += nu_rgb*(drdr_redred+dr_redred/rr[i]+dada_redred/rrsq[i]);
            }
        }
    }
    
    
    
    // Copy the calculated time derivative values into the output array
    memcpy(dt_vars+1,*dt_qvort_r,Nr*Na*sizeof(real));
    memcpy(dt_vars+ 1 + Nr*Na + Np*sizeof(tracer),*dt_redred_r,Nr*Na*sizeof(real));
    
    // BOUNDARY CIRCULATION (see Mcwilliams 1977, Stewart et al. 2014)
          
  
    // Loop over boundary points
    for (j = 0; j < Na; j ++)
    {
        // First-order approximation to linear drag term
        // dt_Gamma_wrk[j] = - kap[i][j] * (psi_r[1][j]-psi_r[0][j])/dr;
      
        // Second-order approximation to linear drag term
        dt_Gamma_wrk[j] = - kap[i][j] * (-psi_r[2][j]+4*psi_r[1][j]-3*psi_r[0][j]) / _2dr;
        
        // Second-order approximation to nu*(d^3psi/dr^3 + 1/r d^psi/dr^2)
        if (nu != 0)
        {
            // dt_Gamma_wrk[j] += (nu/(dr*dr*dr)) * (-psi_r[2][j]+2*psi_r[1][j]-psi_r[0][j]);
          dt_Gamma_wrk[j] += (nu/(dr*dr*dr)) * (-1.5*psi_r[4][j]+7*psi_r[3][j]-12*psi_r[2][j]+9*psi_r[1][j]-2.5*psi_r[0][j])
                           + (nu/(rmin*_drsq)) * (-psi_r[3][j]+4*psi_r[2][j]-5*psi_r[1][j]+2*psi_r[0][j]);
        }
    }
    *dt_Gamma_r = kahanSum(dt_Gamma_wrk,Na)*rmin*da;
  
  
    // TRACER PARTICLE EVOLUTION
  
  
    // If Np > 0 then there are tracer particles to be advected
    if (Np > 0)
    {
        // Calculate an offset so that we can treat all tracer
        // azimuthal positions as positive numbers
        for (p = 0; p < Np; p ++)
        {
            // Extract coordinates in transformed space
            rpos = dyeline_r[p].r;
            apos = dyeline_r[p].a;

            // Find theta-positions of nearest gridpoints. We allow the azimuthal
            // positions to take values outside of [0,La), but obviously the
            // corresponding indices must not lie outside [0,Na-1].
            iprev = floor(Nr*(rpos-rmin)/Lr);
            inext = ceil(Nr*(rpos-rmin)/Lr);
            jprev = floor(Na*apos/La);
            jnext = ceil(Na*apos/La);
            
            // Particle has exited the domain: this should never happen
            // if the code is working correctly, but if it does then this
            // clause avoids a segmentation fault.
            if ((inext > Nr-1) || (iprev < 0))
            {
                continue;
            }
            
            // Positions of nearest grid lines
            rprev = rr[iprev];
            rnext = rr[inext];
            aprev = jprev*La/Na;
            anext = jnext*La/Na;

            // Wrap indices so that they point to the correct positions in Rh
            jprev = jprev - (int) (Na * floor(jprev*1.0/Na));
            jnext = jnext - (int) (Na * floor(jnext*1.0/Na));

            // If the tracer lies exactly on a grid-line, then the positions of
            // the adjacent gridpoints will coincide. In this case we add a space
            // increment so that the linear interpolation below just calculates
            // exactly the value on the coincident gridpoint.
            if (rprev == rnext)
            {
                rnext += dr;
            }
            if (aprev == anext)
            {
                anext += da;
            }

            // Linear interpolation weights
            wrnext = (rpos-rprev)/dr;
            wrprev = (rnext-rpos)/dr;
            wanext = (apos-aprev)/da;
            waprev = (anext-apos)/da;

            // Interpolate to get streamfunction derivatives at the particle position
            dr_psi_tp = dr_psi_r[inext][jnext]*wrnext*wanext
                      + dr_psi_r[inext][jprev]*wrnext*waprev
                      + dr_psi_r[iprev][jnext]*wrprev*wanext
                      + dr_psi_r[iprev][jprev]*wrprev*waprev;

            da_psi_tp = da_psi_r[inext][jnext]*wrnext*wanext
                      + da_psi_r[inext][jprev]*wrnext*waprev
                      + da_psi_r[iprev][jnext]*wrprev*wanext
                      + da_psi_r[iprev][jprev]*wrprev*waprev;

            // Calculate position tendencies
            dt_dyeline_r[p].r = - da_psi_tp / rpos;
            dt_dyeline_r[p].a = dr_psi_tp / rpos;
        }
    }

}



/**
 *
 * writeModelState
 *
 * Writes model variables to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
bool writeModelState (const int t, const int n, real ** qvort, real ** psi, real ** red_c,char * outdir)
{
    char outfile[MAX_PARAMETER_FILENAME_LENGTH];
    char nstr[MAX_PARAMETER_FILENAME_LENGTH];
    FILE * outfd = NULL;
    
    // Create a string for the current iteration number
    sprintf(nstr,"%d",n);
    
    // Write potential vorticity
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_PV);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,qvort,Nr,Na);
    fclose(outfd);

    // Write streamfunction
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_PSI);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,psi,Nr,Na);
    fclose(outfd);
    
    //Isaac
    
    // Write Red Dye Concentration
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_RED);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,red_c,Nr,Na);
    fclose(outfd);
    //End Isaac
    
    return true;
    
    
    
}



/**
 *
 * writeTracerState
 *
 * Writes tracer particle positions to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
bool writeTracerState (const int t, const int n, tracer * dyeline, char * outdir)
{
    char outfile[MAX_PARAMETER_FILENAME_LENGTH];
    char nstr[MAX_PARAMETER_FILENAME_LENGTH];
    FILE * outfd = NULL;
    
    // Create a string for the current iteration number
    sprintf(nstr,"%d",n);
    
    // Write tracer positions
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_TRACER);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printVector(outfd,(real *) dyeline,2*Np);
    fclose(outfd);
    
    return true;
}



/**
 *  printUsage
 *
 *  Prints an information message describing correct usage of the
 *  program in terms of required parameters.
 *
 */
void printUsage()
{
    printf
    (
        "USAGE: DIYsimulate <infile> <outdir> \n"
        "  \n"
        "  where <infile> is the name of the input parameter file, and\n"
        "  <outdir> is the directory into which output files should be\n"
        "  written.\n"
        "  \n"
        "  The input file must specify a series of input parameters\n"
        "  in name-value form, separated by whitespace, i.e.\n"
        "  <name> <value> [<name> <value> [...]]\n"
        "  \n"
        "  All matrix parameters should be formatted in column-major\n"
        "  order.\n"
        "  \n"
        "  name             value\n"
        "  \n"
        "  Nr                Number of r-gridpoints.\n"
        "  Na                Number of theta-gridpoints.\n"
        "  Nt                Number of time-values.\n"
        "  tmax              Time at which integration should stop.\n"
        "                    (Units s.) Must be >0.\n"
        "  rmin              Lower bound of r-range (units m). Must be >0.\n"
        "  rmax              Upper bound of r-range (units m). Must be >rmin.\n"
        "  savefreq          Data storage frequency. An iteration will\n"
        "                    be saved every (savefreq) seconds.\n"
        "                    Must be >0.\n"
        "  savefreqP         Particle output frequency. Particle positions will\n"
        "                    be saved every (savefreqP) seconds. Must be >0.\n"
        "                    Optional - default is savefreq.\n"
        "  f                 Coriolis parameter (units 1/s).\n"
        "                    Optional - default value is 0. \n"
        "  H                 Max. fluid depth (units m). Must be > 0.\n"
        "  kappa             Dissipation constant (units 1/s).\n"
        "                    Optional - default value is 0.\n"
        "  nu                Dynamic viscosity (units m^2/s).\n"
        "                    Optional - default value is 0.\n"
        "  amult             If amult is set, then only a fraction 1/amult\n"
        "                    of the annulus will be computed and treated as\n"
        "                    periodic. Must be an integer > 0.\n"
        "  transportInit     Initial transport along the channel in the positive\n"
        "                    azimuthal-direction (units m^2/s). Superseded \n"
        "                    by streamFuncFile. Default is 0.\n"
        "  vortInitFile      Name of a file containing an Nr x Na matrix.\n"
        "                    Specifies initial relative vorticity at all grid\n"
        "                    points. Cannot be specified with psiInitFile.\n"
        "  psiInitFile       Name of a file containing an Nr x Na matrix.\n"
        "                    Specifies initial streamfunction at all grid points.\n"
        "                    Cannot be specified with vortInitFile.\n"
        "  bathyFile         Name of a file containing an Nr x Na matrix.\n"
        "                    Specifies the elevation of the ocean bed above\n"
        "                    the reference depth H. Optional - default is zeros.\n"
        "  tracInitFile      Name of a file containing an Np x 2 matrix.\n"
        "                    Specifies the initial positions of tracer particles.\n"
        "                    Required if Np>0.\n"
    );
}



/**
 * main
 *
 * Program entry point. Initialises all of the required memory,
 * reads in input parameters, performs iterations of the required
 * numerical method, and finall cleans up.
 *
 */
int main (int argc, char ** argv)
{
    // Time parameters
    real tmin = 0.0;        // Start time
    real tmax = 0;       // End time
    real t = 0;             // Time counter
    real dt = 0;            // Time step
    uint Nt = 0;           // Number of steps
    real dt_s = 0;          // Data save time step
    real t_next = 0;        // Next save time
    uint n_saves = 0;       // Save counter
    real wc = 0;            // Time-interpolation weights
    real wn = 0;

    // Dimensions of arrays to by FFT'd - input parameter for DFT planner
    int tdims[1];

    // Arrays of current variables and buffers for time-interpolation
    // Note: these will point to indices in vars and vars_buf
    real ** qvort = NULL;
    real ** qvort_buf = NULL;
    real ** qvort_out = NULL;
    real * Gamma = NULL;
    real * Gamma_buf = NULL;
    real * Gamma_out = NULL;
    real ** psi_out = NULL;
  
    // For storing tracer positions
    tracer * dyeline = NULL;
    tracer * dyeline_buf = NULL;
    tracer * dyeline_out = NULL;
    
    //Isaac
    //For storing dye vars
    real ** redred = NULL;
    real ** redred_buf = NULL;
    real ** redred_out = NULL;
    //Isaac End
    
    // Work arrays for reading initial conditions
    real ** vortInit = NULL;
    real ** psiInit = NULL;
    complex ** psiInit_s = NULL;
    complex ** dada_psiInit_s = NULL;
    real ** dada_psiInit_r = NULL;
    //Isaac
    //For reading initial dye cond
    real ** redredInit = NULL;
    //Isaac End
    
    // Work arrays for Runge-Kutta algorithms
    real * vars = NULL;
    real * vars_out = NULL;
    real * vars_buf = NULL;
    real * k13 = NULL;
    real * k24 = NULL;

    // Time-stepping arrays - used to store previous iterations
    // for Adams-Bashforth method
    real * dt_vars = NULL;
    real * dt_vars_1 = NULL;
    real * dt_vars_2 = NULL;

    // File descriptors for input parameter file and output data file
    char * infname = argv[1];
    FILE * infile = NULL;
    char * outdir = argv[2];
    char tfilename[MAX_PARAMETER_FILENAME_LENGTH];
    FILE * tfile = NULL;
    time_t now; // To store current time
  
    // For tracer output
    bool writeTracers = false;   // Set true if a tracer output file is specified
    real dt_p = 0;               // Tracer save time step
    real t_next_p = 0;           // Next tracer save time
    uint n_p = 0;                // Tracer save counter

    // Looping variables
    uint i = 0;
    uint j = 0;
    uint n = 0;
    uint jm1 = 0;
    uint jp1 = 0;

    // For creating initial conditions
    bool initVort = false;
    bool initPsi = false;
    real psi_offset = 0;
    real psi0Init = 0;
    real dr_psiInit = 0;
    real drdr_psiInit = 0;

    // Set this flag to perform a check for instability
    bool checkstability = false;

    // Stores data required for parsing input parameters
    paramdata params[NPARAMS];
  
    // Filename holders for input parameter arrays
    char bathyFile[MAX_PARAMETER_FILENAME_LENGTH];
    char vortInitFile[MAX_PARAMETER_FILENAME_LENGTH];
    char psiInitFile[MAX_PARAMETER_FILENAME_LENGTH];
    char tracInitFile[MAX_PARAMETER_FILENAME_LENGTH];
    //Isaac
    // Filename holders for dye input parameter arrays
    char redInitFile[MAX_PARAMETER_FILENAME_LENGTH];
    //Isaac End
    // Default file name parameters - zero-length strings
    bathyFile[0] = '\0';
    vortInitFile[0] = '\0';
    psiInitFile[0] = '\0';
    tracInitFile[0] = '\0';
    //Isaac
    redInitFile[0] = '\0';
    //Isaac End
    // Define input parameter data
    setParam(params,0,"Nr","%u",&Nr,false);
    setParam(params,1,"Na","%u",&Na,false);
    setParam(params,2,"Nt","%u",&Nt,false);
    setParam(params,3,"Np","%u",&Np,false);
    setParam(params,4,"savefreq","%le",&dt_s,true);
    setParam(params,5,"savefreqP","%le",&dt_p,true);
    setParam(params,6,"amult","%u",&amult,true);
    setParam(params,7,"tmax","%lf",&tmax,false);
    setParam(params,8,"rmin","%lf",&rmin,false);
    setParam(params,9,"rmax","%lf",&rmax,false);
    setParam(params,10,"f","%lf",&f,true);
    setParam(params,11,"H","%lf",&H,true);
    setParam(params,12,"kappa","%lf",&kap0,true);
    setParam(params,13,"nu","%lf",&nu,true);
    setParam(params,14,"transportInit","%lf",&psi0Init,true);
    setParam(params,15,"bathyFile","%s",&bathyFile,true);
    setParam(params,16,"vortInitFile","%s",&vortInitFile,true);
    setParam(params,17,"psiInitFile","%s",&psiInitFile,true);
    setParam(params,18,"tracInitFile","%s",&tracInitFile,true);
    //Isaac
    setParam(params,19,"redTracInitFile","%s",&redInitFile,true);
    //Isaac End
    
    // Check that a file name has been specified
    if (argc < 3)
    {
        fprintf(stderr,"ERROR: No input parameter file name supplied\n");
        printUsage();
        return 0;
    }

    // Calculate elapsed time and write to a file
    strcpy(tfilename,outdir);
    strcat(tfilename,"/");
    strcat(tfilename,OUTN_TFILE);
    tfile = fopen(tfilename,"w");
    if (tfile != NULL)
    {
      time(&now);
      fprintf(tfile,"Program started at %s\n", ctime(&now));
      fflush(tfile);
    }
  
  
  
    //////////////////////////////////////////
    ///// BEGIN READING INPUT PARAMETERS /////
    //////////////////////////////////////////
    
    // Attempt to read input parameter data. Errors will be printed
    // to stderr, and will result in 'false' being returned.
    if (!readParams(argv[1],params,NPARAMS,stderr))
    {
        printUsage();
        return 0;
    }

    // Check that required parameters take legal values
    if ( (Na == 0) || (Nr == 0) || (Nt == 0) || (tmax <= 0.0) ||
         (rmin <= 0) || (rmax <= rmin) || (amult == 0) || (H <= 0.0) )
    {
        fprintf(stderr,"ERROR: Invalid input parameter values\n");
        printUsage();
        return 0;
    }
  
    // Calculate important values from the parameters
    N = 1 + Nr*Na + Np*sizeof(tracer) + Nr*Na;
    //N = 1 + Nr*Na + Np*sizeof(tracer);

    dr = (rmax-rmin)/(Nr-1);
    da = _2PI / (Na*amult);
    dt = (tmax-tmin)/(Nt-1);
    t = tmin;
    t_next = tmin + dt_s;
    Lr = rmax-rmin; // Channel width in rho/phi-space
    La = _2PI/amult; // Channel angular length
    
    // To avoid repeated calculation
    _dr_da = dr/da;
    _4drda = 4*dr*da;
    _dr_dasq = _dr_da*_dr_da;
    _drsq = dr*dr;
    _2dr = 2*dr;
    _fH = f/H;
    _dasq = da*da;
    _2da = 2*da;
    amultsq = amult*amult;
  
    // If no output frequency specified, use the time step
    if (dt_s <= 0.0)
    {
        dt_s = dt;
    }
  
    // If tracer output frequency isn't specified then it is set to dt_s
    if (dt_p <= 0.0)
    {
        dt_p = dt_s;
    }
  
    // Check that initial conditions have been specified
    initVort = strlen(vortInitFile)>0;
    initPsi = strlen(psiInitFile)>0;
    if (!initVort && !initPsi)
    {
      fprintf(stderr,"ERROR: Must specify initial streamfunction or relative vorticity\n");
      printUsage();
      return 0;
    }
    
    // Check we're not trying to specify both psi and PV
    if (initVort && initPsi)
    {
        fprintf(stderr,"ERROR: Can only specify initial streamfunction OR relative vorticity\n");
        printUsage();
        return 0;
    }


  
    ////////////////////////////////////////
    ///// END READING INPUT PARAMETERS /////
    ////////////////////////////////////////

  
  
    ///////////////////////////////////
    ///// BEGIN MEMORY ALLOCATION /////
    ///////////////////////////////////


    // All-variable storage buffers
    CHECKALLOC(vars,N*sizeof(real));
    CHECKALLOC(vars_out,N*sizeof(real));
    CHECKALLOC(vars_buf,N*sizeof(real));
    CHECKALLOC(dt_vars,N*sizeof(real));
    CHECKALLOC(dt_vars_1,N*sizeof(real));
    CHECKALLOC(dt_vars_2,N*sizeof(real));
    CHECKALLOC(k13,N*sizeof(real));
    CHECKALLOC(k24,N*sizeof(real));
  
    // Pointers to data stored in 'vars' buffer
    CHECKALLOC(qvort,Nr*sizeof(real *));
    CHECKALLOC(qvort_buf,Nr*sizeof(real *));
    CHECKALLOC(qvort_out,Nr*sizeof(real *));
    CHECKALLOC(redred,Nr*sizeof(real *));
    CHECKALLOC(redred_buf,Nr*sizeof(real *));
    CHECKALLOC(redred_out,Nr*sizeof(real *));
    Gamma = vars;
    vec2mat(vars+1,&qvort,Nr,Na);
    dyeline = (tracer *) (vars+1+Nr*Na);
    vec2mat(vars+1+Nr*Na+Np*sizeof(tracer),&redred,Nr,Na);
    
    Gamma_buf = vars_buf;
    vec2mat(vars_buf+1,&qvort_buf,Nr,Na);
    dyeline_buf = (tracer *) (vars_buf+1+Nr*Na);
    vec2mat(vars_buf+1+Nr*Na+Np*sizeof(tracer),&redred_buf,Nr,Na);
    
    Gamma_out = vars_out;
    vec2mat(vars_out+1,&qvort_out,Nr,Na);
    dyeline_out = (tracer *) (vars_out+1+Nr*Na);
    vec2mat(vars_out+1+Nr*Na+Np*sizeof(tracer),&redred_out,Nr,Na);
    
    RMATALLOC(psi_out,Nr,Na);
  
    // Initial conditions
    RMATALLOC(vortInit,Nr,Na);
    RMATALLOC(psiInit,Nr,Na);
    RMATALLOC(dada_psiInit_r,Nr,Na);
    CMATALLOC(psiInit_s,Nr,Na);
    CMATALLOC(dada_psiInit_s,Nr,Na);
    
    CHECKALLOC(redredInit,Nr*sizeof(real *)); //do i need this Jordyn??
    RMATALLOC(redredInit,Nr,Na);
    // Wrapper matrices for input arrays in 'tderiv' function
    CHECKALLOC(dt_qvort_r,Nr*sizeof(real *));
    CHECKALLOC(qvort_r,Nr*sizeof(real *));
    
    // Wrapper matrices for input arrays in 'tderiv' function for the dye part
    CHECKALLOC(redred_r,Nr*sizeof(real *));
    CHECKALLOC(dt_redred_r,Nr*sizeof(real *));
    
    // Model state variables in real and spectral space
    RMATALLOC(qvort_r,Nr,Na);
    RMATALLOC(dt_qvort_r,Nr,Na);
    RMATALLOC(vort_r,Nr,Na);
    RMATALLOC(psi_r,Nr,Na);
    
    RMATALLOC(redred_r,Nr,Na);
    RMATALLOC(dt_redred_r,Nr,Na);
    
    CMATALLOC(vort_s,Nr,Na);
    CMATALLOC(qvort_s,Nr,Na);
    CMATALLOC(psi_s,Nr,Na);
    
    CMATALLOC(redred_s,Nr,Na);
    
    // Work arrays for tderiv
    RMATALLOC(vq_r,Nr,Na);
    CMATALLOC(vq_s,Nr,Na);
    RMATALLOC(vc_Red_r,Nr,Na);
    CMATALLOC(vc_Red_s,Nr,Na);
    RMATALLOC(da_vq_r,Nr,Na);
    RMATALLOC(da_vc_Red_r,Nr,Na);
    RMATALLOC(da_psi_r,Nr,Na);
    RMATALLOC(dr_psi_r,Nr,Na);
    CMATALLOC(da_vq_s,Nr,Na);
    CMATALLOC(da_vc_Red_s,Nr,Na);
    CMATALLOC(da_psi_s,Nr,Na);
    RMATALLOC(dada_vort_r,Nr,Na);
    CMATALLOC(dada_vort_s,Nr,Na);
    RMATALLOC(dada_redred_r,Nr,Na);
    CMATALLOC(dada_redred_s,Nr,Na);
    
  
    // Matrices for streamfunction decomposition
    RMATALLOC(qvortL,Nr,Na);
    RMATALLOC(psiL,Nr,Na);
    RMATALLOC(psiP,Nr,Na);
  
    // For boundary circulation evolution
    CHECKALLOC(u_bdy,Na*sizeof(real));
    CHECKALLOC(dt_Gamma_wrk,Na*sizeof(real));

    // Matrices for streamfunction inversion
    RMATALLOC(psicoeffb,Na,Nr);
    CHECKALLOC(psicoeffa,Nr*sizeof(real));
    CHECKALLOC(psicoeffc,Nr*sizeof(real));
    CHECKALLOC(psitempc,Nr*sizeof(real));
    CHECKALLOC(psitempd,Nr*sizeof(complex));
    CHECKALLOC(psitempx,Nr*sizeof(complex));

    // Grid matrices
    CHECKALLOC(rr,Nr*sizeof(real));
    CHECKALLOC(rrsq,Nr*sizeof(real));
    CHECKALLOC(kk,Na*sizeof(real));
    CHECKALLOC(kksq,Na*sizeof(real));
    RMATALLOC(hh,Nr,Na);

    RMATALLOC(kap,Nr,Na);
    CHECKALLOC(filter,Na*sizeof(real));

    // Create FFT plans
    tdims[0] = Na;
    forward_transform = fftw_plan_many_dft_r2c(1,tdims,Nr,*vort_r,NULL,1,Na,*vort_s,NULL,1,Na,FFTW_MEASURE);
    backward_transform = fftw_plan_many_dft_c2r(1,tdims,Nr,*vort_s,NULL,1,Na,*vort_r,NULL,1,Na,FFTW_MEASURE);
    if (forward_transform == NULL || backward_transform == NULL)
    {
        fprintf(stderr,"ERROR: Unable to create DFT plans\n");
        return 0;
    }
  
  
    /////////////////////////////
    ///// BEGIN WORK ARRAYS /////
    /////////////////////////////
    nu_rgb = nu/8;
    // Calculate r-values on gridpoints
    rr[0] = rmin;
    rrsq[0] = SQUARE(rmin);
    for (i = 1; i < Nr; i++)
    {
        rr[i] = rr[i-1]+dr;
        rrsq[i] = SQUARE(rr[i]);
    }

    // Calculate k-values - wave numbers in DFT
    for (j = 0; j < Na; j++)
    {
        if (j < (Na-Na/2))
        {
            kk[j] = j;
        }
        else
        {
            kk[j] = j-Na;
        }

        kksq[j] = SQUARE(kk[j]);
    }

    // Calculate the Fourier smoothing filter values (see Li and Hou 2007)
    for (j = 0; j < Na; j ++)
    {
        filter[j] = exp(-fil_a*pow(abs(kk[j])/(Na/2.0),fil_m));
    }

    // Set coefficient arrays used for calculation of the streamfunction
    // NOTE: 0th and (Nr-1)th indices are not used
    for (i = 1; i < Nr-1; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            psicoeffb[j][i] = -(2 + dr*dr*amultsq*kksq[j]/rrsq[i]);
        }

        psicoeffa[i] = 1-dr/(2*rr[i]);
        psicoeffc[i] = 1+dr/(2*rr[i]);
    }
  
    // Dissipation coefficient
    for (i = 0; i < Nr; i ++)
    {
        for (j = 0; j < Na; j ++)
        {
            kap[i][j] = kap0;
        }
    }

    ///////////////////////////
    ///// END WORK ARRAYS /////
    ///////////////////////////
  

    ////////////////////////////////////////
    ///// BEGIN READING PARAMETER DATA /////
    ////////////////////////////////////////
  
  
    // Read bottom elevation
    if (strlen(bathyFile) > 0)
    {
        if (!readMatrix(bathyFile,hh,Nr,Na,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",bathyFile);
            printUsage();
            return 0;
        }
    }
    else
    {
        // Default to no bottom topography
        for (i = 0; i < Nr; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                hh[i][j] = 0;
            }
        }
    }
    //Isaac
    //read red dye initial concentration
    if (strlen(redInitFile) > 0)
    {
        if (!readMatrix(redInitFile,redredInit,Nr,Na,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",redInitFile);
            printUsage();
            return 0;
        }
    }
    else
    {
        // Default to no dye topography
        for (i = 0; i < Nr; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                redredInit[i][j] = 0;
            }
        }
    }
    //End Isaac

    // Read relative vorticity initialization data and compute initial streamfunction
    if (initVort)
    {
        if (!readMatrix(vortInitFile,vortInit,Nr,Na,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",vortInitFile);
            printUsage();
            return 0;
        }
      
        // Ensure relative vorticity is zero on boundaries (free slip condition)
        // NOTE: remember that qvort is storing relative vorticity at this point in the code
        for (j = 0; j < Na; j ++)
        {
            if ((vortInit[0][j] != 0) || (vortInit[Nr-1][j] != 0))
            {
                fprintf(stderr,"ERROR: Relative vorticity must be zero on the boundaries\n");
                printUsage();
                return 0;
            }
        }
        
        // Construct PV
        for (i = 0; i < Nr; i++)
        {
            for (j = 0; j < Na; j++)
            {
                qvort[i][j] = vortInit[i][j] + f*hh[i][j]/H;
            }
        }
        
        // Calculate the initial streamfunction - we'll need this to compute the initial circulation
        // around the offshore wall
        invertPsi(qvort,psi0Init,false);
        memcpy(*psiInit,*psi_r,Nr*Na*sizeof(real));
      
    } // end if (initVort)
    
    // Read streamfunction initialization data and compute relative vorticity
    if (initPsi)
    {
        if (!readMatrix(psiInitFile,psiInit,Nr,Na,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",psiInitFile);
            printUsage();
            return 0;
        }
      
        // Ensure streamfunction is uniform on boundaries (necessary for our scheme)
        for (j = 0; j < Na; j ++)
        {
            jm1 = (j+Na-1) % Na;
            if ((psiInit[0][j] != psiInit[0][jm1]) || (psiInit[Nr-1][j] != psiInit[Nr-1][jm1]))
            {
                fprintf(stderr,"ERROR: Streamfunction must be uniform on the boundaries\n");
                printUsage();
                return 0;
            }
        }
       
        // Without loss of generality we can add a constant to the streamfunction
        // so that it is zero along the offshore wall (x=0)
        psi_offset = psiInit[0][0];
        if (psi_offset != 0)
        {
            for (i = 0; i < Nr; i ++)
            {
                for  (j = 0; j < Na; j ++)
                {
                    psiInit[i][j] = psiInit[i][j] - psi_offset;
                }
            }
        }
    
        // Streamfunction on the outer boundary - initial along-channel transport
        psi0Init = psiInit[Nr-1][0];
      
        // First compute second derivative in azimuthal direction
        if (method_a == METHOD_PS)
        {
            // Transform qvort_r -> qvort_s and vq_r -> vq_s
            fftw_execute_dft_r2c(forward_transform,*psiInit,*psiInit_s);

            // Calculate theta-derivatives in spectral space:
            // da_psi_s, da_vort_s and dada_vort_s
            for (i = 0; i < Nr; i ++)
            {
                for (j = 0; j < Na; j ++)
                {
                    dada_psiInit_s[i][j][0] = -amultsq*kksq[j]*psiInit_s[i][j][0];
                    dada_psiInit_s[i][j][1] = -amultsq*kksq[j]*psiInit_s[i][j][1];
                }
            }

            // Transform back to get dada_psiInit_r
            fftw_execute_dft_c2r(backward_transform,*dada_psiInit_s,*dada_psiInit_r);
            normalise(*dada_psiInit_r,Nr*Na,Na);
        
        }
        // Finite-differencing method
        else if (method_a == METHOD_FD)
        {
            for (i = 0; i < Nr; i ++)
            {
                for (j = 0; j < Na; j ++)
                {
                    jp1 = (j+1) % Na;
                    jm1 = (j+Na-1) % Na;
                    dada_psiInit_r[i][j] = (psiInit[i][jp1]-2*psiInit[i][j]+psiInit[i][jm1]) / _dasq;
                }
            }
        }
                  
        // Add second radial derivative to compute initial relative vorticity
        for (i = 1; i < Nr-1; i ++)
        {
            for (j = 0; j < Na; j ++)
            {
                dr_psiInit = (psiInit[i+1][j]-psiInit[i-1][j]) / _2dr;
                drdr_psiInit = (psiInit[i+1][j]-2*psiInit[i][j]+psiInit[i-1][j]) / _drsq;
                vortInit[i][j] = drdr_psiInit + dr_psiInit/rr[i] + dada_psiInit_r[i][j]/rrsq[i];
            }
        }
      
        // Free-slip boundary conditions
        for (j = 0; j < Na; j ++)
        {
            vortInit[0][j] = 0;
            vortInit[Nr-1][j] = 0;
        }
        
        // Add planetary PV
        for (i = 0; i < Nr; i++)
        {
            for (j = 0; j < Na; j++)
            {
                qvort[i][j] = vortInit[i][j] + f*hh[i][j]/H;
            }
        }
      
    } // end if (initPsi)
    
        
    // Calculate the initial circulation around the southern wall
    *Gamma = calcGamma(psiInit);
    
    // Initial tracer positions
    if (Np > 0)
    {
        if (!readVector(tracInitFile,(real *) dyeline,2*Np,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",tracInitFile);
            printUsage();
            return 0;
        }
    }
    
    //Initial red dye concentration
    if (strlen(redInitFile)>0)
    {
        if (!readMatrix(redInitFile,redredInit,Nr,Na,stderr))
        {
            fprintf(stderr,"ERROR: Could not read data from %s\n",redInitFile);
            printUsage();
            return 0;
        }
      
        // Ensure relative dye concentration is zero on boundaries (free slip condition)
        // NOTE: remember that qvort is storing relative vorticity at this point in the code
        for (j = 0; j < Na; j ++)
        {
            //enforce violently bc i cant figure out how to get it to work in matlab
            
            redredInit[0][j] = 0;
            redredInit[Nr-1][j] = 0;
            
        }
        for (j = 0; j < Na; j ++)
        {
            
            if ((redredInit[0][j] != 0) || (redredInit[Nr-1][j] != 0))
            {
                fprintf(stderr,"ERROR: Red Dye concentration must be zero on the boundaries\n");
                printUsage();
                return 0;
            }
        }
        
        
        
        for (i = 0; i < Nr; i++)
        {
            for (j = 0; j < Na; j++)
            {
                redred[i][j] = redredInit[i][j];
            }
        }
      
    } // end if (strlen(redInitFile)>0)
    
    
    //////////////////////////////////////
    ///// END READING PARAMETER DATA /////
    //////////////////////////////////////


  
    /////////////////////////////////////////////
    ///// BEGIN STREAMFUNCTION SOLVER SETUP /////
    /////////////////////////////////////////////
        
    // Initial streamfunction on the outer boundary. We use this as a
    // guess for the typical transport in the channel when calculating
    // the streamfunction at later times.
    psiL0 = psi0Init;
    if (psiL0 == 0)
    {
        // Must not be zero
        psiL0 = 1;
    }
    
    // Calculate psiL, the solution to Laplace's equation - N.B. the invertPsi
    // function subtracts off qvortL from the supplied PV array, so this is
    // guaranteed to compute the solution to Laplace's equation
    invertPsi(qvortL,psiL0,false);
    memcpy(*psiL,*psi_r,Nr*Na*sizeof(real));
    
    // Calculate GammaL, the circulation around the inner wall
    // in the solution to Laplace's equation
    GammaL = calcGamma(psiL);
    
    ///////////////////////////////////////////
    ///// END STREAMFUNCTION SOLVER SETUP /////
    ///////////////////////////////////////////
  
  

    // Open the tracer file
    if (Np > 0)
    {
        writeTracers = true;
    }
    
    // Write out the initial data
    if (!writeModelState(t,n_saves,qvort,psiInit,redred,outdir))//Isaac modified
    {
        fprintf(stderr,"Unable to write model state");
        printUsage();
        return 0;
    }
    
    // Write initial tracer data
    if (writeTracers)
    {
        if (!writeTracerState(t,n_saves,dyeline,outdir))
        {
            fprintf(stderr,"Unable to write tracer state");
            printUsage();
            return 0;
        }
    }
      
    // Increment recorded number of output writes
    n_saves ++;
  
    // Initial time record
    if (tfile != NULL)
    {
        fprintf(tfile,"%f ",t);
        fflush(tfile);
    }

    // Numerical time-integration loop
    for (n = 1; n < Nt; n++)
    {
        if (method_t == METHOD_RK4)
        {
            // Generate the next iteration and copy back to the 'qvort' matrix
            rk4 (   &t,
                    vars,
                    vars_out,
                    k13,
                    k24,
                    dt,
                    N,
                    &tderiv             );

            memcpy(*qvort,*qvort_out,N*sizeof(real));
            memcpy(*redred,*redred_out,N*sizeof(real));//Hey why do we only need to do this in rk4?
        }
        else if (method_t == METHOD_AB3)
        {
            // Generate the first couple of iterations
            // with lower-order methods
            if (n == 1)
            {
                ab1 (&t,
                     vars,
                     vars_out,
                     dt_vars,
                     dt,
                     N,
                     &tderiv             );
            }
            else if (n == 2)
            {
                ab2 (&t,
                     vars,
                     vars_out,
                     dt_vars,
                     dt_vars_1,
                     dt,
                     N,
                     &tderiv             );
              
            }
            // Once we have previous iterations, use third order Adams-Bashforth
            else
            {
                ab3 (&t,
                     vars,
                     vars_out,
                     dt_vars,
                     dt_vars_1,
                     dt_vars_2,
                     dt,
                     N,
                     &tderiv             );
            }
        }
        else
        {
            fprintf(stderr,"ERROR: Unknown time-integration method\n");
            break;
        }

        // If we're using a pseudospectral method, apply spectral filtering here
        if (method_a == METHOD_PS)
        {
            fftw_execute_dft_r2c(forward_transform,*qvort,*qvort_s);

            fftw_execute_dft_r2c(forward_transform,*redred,*redred_s);
            
            for (i = 0; i < Nr; i ++)
            {
                for (j = 0; j < Na; j ++)
                {
                    qvort_s[i][j][0] *= filter[j];
                    qvort_s[i][j][1] *= filter[j];
                    
                    redred_s[i][j][0] *= filter[j];
                    redred_s[i][j][1] *= filter[j];
                }
            }

            fftw_execute_dft_c2r(backward_transform,*qvort_s,*qvort);
            normalise(*qvort,Nr*Na,Na);
            
            fftw_execute_dft_c2r(backward_transform,*redred_s,*redred);
            normalise(*redred,Nr*Na,Na);
        }

        // If the time step has taken us past a save point (or multiple
        // save points), interpolate and write out the data
        while (t >= t_next)
        {
            // Write out the latest model state
            calcPsi(qvort_out,*Gamma_out,psi_out);
            
            if (!writeModelState(t,n_saves,qvort_out,psi_out,redred_out,outdir))//Isaac Modified
            {
                fprintf(stderr,"Unable to write model state");
                printUsage();
                return 0;
            }
            
            // Increment recorded number of output writes
            n_saves ++;
            
            // Log time in time file
            if (tfile != NULL)
            {
                fprintf(tfile,"%e ",t_next);
                fflush(tfile);
            }
            
            // Show progress
            fprintf(stdout,"Simulation progress: %d out of %d time steps\n",n,Nt);
            fflush(stdout);
            
            // Update the next save time
            t_next = n_saves*dt_s;
        }
        
        // If the time step has taken us past a tracer save point (or multiple
        // save points), interpolate and write out the data
        while (writeTracers && (t >= t_next_p))
        {
            wc = (t-t_next_p) / dt;
            wn = 1 - wc;
          
            // Interpolate values at t_next_p into the dyeline_buf array
            for (i = 0; i < Np; i ++)
            {
                dyeline_buf[i].r = wc*dyeline[i].r + wn*dyeline_out[i].r;
                dyeline_buf[i].a = wc*dyeline[i].a + wn*dyeline_out[i].a;
            }
            
            // Write out the interpolated data from the dyeline_buf array
            if (!writeTracerState(t,n_p,dyeline_buf,outdir))
            {
                fprintf(stderr,"Unable to write tracer state");
                printUsage();
                return 0;
            }
            
            // Increment recorded number of tracer saves
            n_p ++;
            
            // Update the next save time
            t_next_p = n_p*dt_p;
        }
      
        // Update with the new iteration data
        memcpy(vars,vars_out,N*sizeof(real));
    }

    // Print completion time
    if (tfile != NULL)
    {
        time(&now);
        fprintf(tfile,"\nProgram completed at %s\n", ctime(&now));
        fclose(tfile);
    }


    return 0;
}
