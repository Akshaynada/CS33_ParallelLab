//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name :
 * UCLA ID:
 * Email id:
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int OMP_xMax;
#define xMax OMP_xMax
int OMP_yMax;
#define yMax OMP_yMax
int OMP_zMax;
#define zMax OMP_zMax

int OMP_Index(int x, int y, int z)
{
    return ((z * yMax + y) * xMax + x);
}
#define Index(x, y, z) OMP_Index(x, y, z)

double OMP_SQR(double x)
{
    return pow(x, 2.0);
}
#define SQR(x) OMP_SQR(x)

double* OMP_conv;
double* OMP_g;

void OMP_Initialize(int xM, int yM, int zM)
{
    xMax = xM;
    yMax = yM;
    zMax = zM;
    assert(OMP_conv = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
    assert(OMP_g = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
}
void OMP_Finish()
{
    free(OMP_conv);
    free(OMP_g);
}
void OMP_GaussianBlur(double *u, double Ksigma, int stepCount)
{
    double lambda = (Ksigma * Ksigma) / (double)(2 * stepCount);
    double nu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    int x, y, z, step;
    double boundryScale = 1.0 / (1.0 - nu);
    double postScale = pow(nu / lambda, (double)(3 * stepCount));

    // My own constants
    int yMaxTxMax = xMax * yMax;

    // 1st for loop block
    for(step = 0; step < stepCount; step++)
    {
        /*
        // A - Operate on the x=0 "yz-face" with no dependencies
        // Changed from y,z to z,y for slightly improved spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                u[Index(0, y, z)] *= boundryScale;
            }
        }
        */

        // B - Operate on yz-slices where x=n depends on value of x=n-1
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                int currIndex = Index(0, y, z);
                // The Loop A above can be moved into Loop B
                // by moving the boundry operation on the "x=0" face within the yz-loop
                u[currIndex] *= boundryScale;

                for(x = 1; x < xMax; x++)
                {
                    u[currIndex + x] += u[currIndex + x - 1] * nu;
                }

                // The Loop C below can be moved into Loop B
                // by moving the following statement within the yz-loop
                u[currIndex] *= boundryScale;
            }
        }

        /*
        // C - Operate on the x=0 "yz-face" (Same behavior as Loop A)
        // Changed from y,z to z,y for slightly improved spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                u[Index(0, y, z)] *= boundryScale;
            }
        }
        */

        // D - Operate on yz-slices where x=n depends on value of x=n+1
        // Changed from x,y,z to z,y,x
        for(z = 0; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                for(x = xMax - 2; x >= 0; x--)
                {
                    int currIndex = Index(x, y, z);
                    // Changed from Index(x + 1, y, z) to currIndex + 1
                    u[currIndex] += u[currIndex + 1] * nu;
                }
            }
        }

        // E - Operate on the y=0 "xz-face"
        // Changed from x,z to z,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(x = 0; x < xMax; x++)
            {
                u[Index(x, 0, z)] *= boundryScale;
            }
        }

        // F - Operate on xz-slices where y=n depends on value of y=n-1
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = 1; y < yMax; y++)
            {
                for(x = 0; x < xMax; x++)
                {
                    int currIndex = Index(x, y, z);
                    u[currIndex] += u[currIndex - xMax] * nu;
                }
            }
        }

        // G - Operate on y = yMax - 1 (xz-face) with no dependencies
        // Changed from x,y to y,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(x = 0; x < xMax; x++)
            {
                u[Index(x, yMax - 1, z)] *= boundryScale;
            }
        }

        // H - Operate on xz-slices where y=n depends on value of y=n+1
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = yMax - 2; y >= 0; y--)
            {
                for(x = 0; x < xMax; x++)
                {
                    int currIndex = Index(x, y, z);
                    u[currIndex] += u[currIndex + xMax] * nu;
                }
            }
        }

        // I - Operate on z = 0 (xy-face) with no dependencies
        // Changed from x,y to y,x for spatial locality
        for(y = 0; y < yMax; y++)
        {
            for(x = 0; x < xMax; x++)
            {
                u[Index(x, y, 0)] *= boundryScale;
            }
        }

        // J - Operate on xy-slices where z=n depends on value of z=n-1
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 1; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                for(x = 0; x < xMax; x++)
                {
                    int currIndex = Index(x, y, z);
                    u[currIndex] = u[currIndex - yMaxTxMax] * nu;
                }
            }
        }

        // K - Operate on z = zMax - 1 (xy-face) with no dependencies
        // Changed from x,y to y,x
        for(y = 0; y < yMax; y++)
        {
            for(x = 0; x < xMax; x++)
            {
                u[Index(x, y, zMax - 1)] *= boundryScale;
            }
        }

        // L - Operate on xy-slices where z=n depends on value of z=n+1
        // Changed from x, y, z to z, y, x for improved spatial locality
        for (z = zMax - 2; z >= 0; z--)
        {
            for (y = 0; y < yMax; y++)
            {
                for (x = 0; x < xMax; x++)
                {
                    // Switched Index(x,y,z+1) to currIndex + (yMax*xMax)
                    int currIndex = Index(x, y, z);
                    u[currIndex] += u[currIndex + yMaxTxMax] * nu;
                }
            }
        }
    }

    // 2nd for loop block
    // Changed x, y, z to z, y, x to take advantage of caching
    for (z = 0; z < zMax; z++)
    {
        for (y = 0; y < yMax; y++)
        {
            for (x = 0; x <= xMax; x++)
            {
                u[Index(x, y, z)] *= postScale;
            }
        }
    }
}
void OMP_Deblur(double* u, const double* f, int maxIterations, double dt, double gamma, double sigma, double Ksigma)
{
    double epsilon = 1.0e-7;
    double sigma2 = SQR(sigma);
    int x, y, z, iteration;
    int converged = 0;
    int lastConverged = 0;
    int fullyConverged = (xMax - 1) * (yMax - 1) * (zMax - 1);
    double* conv = OMP_conv;
    double* g = OMP_g;

    // My own constants
    int xMaxTyMax = xMax * yMax;

    for(iteration = 0; iteration < maxIterations && converged != fullyConverged; iteration++)
    {
        // 1st for loop block
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 1; z < zMax - 1; z++)
        {
            for(y = 1; y < yMax - 1; y++)
            {
                for(x = 1; x < xMax - 1; x++)
                {
                    g[Index(x, y, z)] = 1.0 / sqrt(epsilon + 
						SQR(u[Index(x, y, z)] - u[Index(x + 1, y, z)]) + 
						SQR(u[Index(x, y, z)] - u[Index(x - 1, y, z)]) + 
						SQR(u[Index(x, y, z)] - u[Index(x, y + 1, z)]) + 
						SQR(u[Index(x, y, z)] - u[Index(x, y - 1, z)]) + 
						SQR(u[Index(x, y, z)] - u[Index(x, y, z + 1)]) + 
						SQR(u[Index(x, y, z)] - u[Index(x, y, z - 1)]));
                    /*
                    int currIndex = Index(x, y, z);
                    g[Index(x, y, z)] = 1.0 / sqrt(epsilon + 
                        SQR(u[currIndex] - u[currIndex + 1]) + 
                        SQR(u[currIndex] - u[currIndex - 1]) + 
                        SQR(u[currIndex] - u[currIndex + xMax]) + 
                        SQR(u[currIndex] - u[currIndex - xMax]) + 
                        SQR(u[currIndex] - u[currIndex + xMaxTyMax]) + 
                        SQR(u[currIndex] - u[currIndex - xMaxTyMax]));
                        */
                }
            }
        }
        memcpy(conv, u, sizeof(double) * xMax * yMax * zMax);
        OMP_GaussianBlur(conv, Ksigma, 3);

        // 2nd for loop block
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 0; z < zMax; z++)
        {
            for(y = 0; y < yMax; y++)
            {
                for(x = 0; x < xMax; x++)
                {
                    double r = conv[Index(x, y, z)] * f[Index(x, y, z)] / sigma2;
                    r = (r * (2.38944 + r * (0.950037 + r))) / (4.65314 + r * (2.57541 + r * (1.48937 + r)));
                    conv[Index(x, y, z)] -= f[Index(x, y, z)] * r;
                }
            }
        }
        OMP_GaussianBlur(conv, Ksigma, 3);
        converged = 0;

        // 3rd for loop block
        // Changed from x,y,z to z,y,x for spatial locality
        for(z = 1; z < zMax - 1; z++)
        {
            for(y = 1; y < yMax - 1; y++)
            {
                for(x = 1; x < xMax - 1; x++)
                {
                    double oldVal = u[Index(x, y, z)];
                    double newVal = (u[Index(x, y, z)] + dt * ( 
                        u[Index(x - 1, y, z)] * g[Index(x - 1, y, z)] + 
                        u[Index(x + 1, y, z)] * g[Index(x + 1, y, z)] + 
                        u[Index(x, y - 1, z)] * g[Index(x, y - 1, z)] + 
                        u[Index(x, y + 1, z)] * g[Index(x, y + 1, z)] + 
                        u[Index(x, y, z - 1)] * g[Index(x, y, z - 1)] + 
                        u[Index(x, y, z + 1)] * g[Index(x, y, z + 1)] - gamma * conv[Index(x, y, z)])) /
                        (1.0 + dt * (g[Index(x + 1, y, z)] + g[Index(x - 1, y, z)] + g[Index(x, y + 1, z)] + g[Index(x, y - 1, z)] + g[Index(x, y, z + 1)] + g[Index(x, y, z - 1)]));
                    if(fabs(oldVal - newVal) < epsilon)
                    {
                        converged++;
                    }
                    u[Index(x, y, z)] = newVal;
                }
            }
        }
        if(converged > lastConverged)
        {
            printf("%d pixels have converged on iteration %d\n", converged, iteration);
            lastConverged = converged;
        }
    }
}
