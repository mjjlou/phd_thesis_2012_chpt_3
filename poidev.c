#include <math.h>
#define PI 3.141592654
float poidev(float xm, long *idum)
//Returns as a floating-point number an integer value that is a random deviate drawn from a
//Poisson distribution of mean xm, using ran1(idum) as a source of uniform random deviates.
{
    float gammln(float xx);
    float ran1(long *idum);
    static float sq,alxm,g,oldm=(-1.0);
    float em, t, y;

    if (xm < 12.0) {
        if (xm != oldm) {
            oldm=xm;
            g=exp(-xm); //If xm is new, compute the exponential.
        }
        em = -1;
        t=1.0;
        do { //Instead of adding exponential deviates it is equivalent
            //to multiply uniform deviates. We never
            //actually have to take the log, merely compare
            //to the pre-computed exponential.
            ++em;
            t *= ran1(idum);
        } while (t > g);
    } else { //Use rejection method.
        if (xm != oldm) { //If xm has changed since the last call, then pre-
            oldm=xm; // compute some functions that occur below.
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
            //The function gammln is the natural log of the gamma function, as given in §6.1.
        }
        do {
            do { //y is a deviate from a Lorentzian comparison function
                y=tan(PI*ran1(idum));
                em=sq*y+xm; //em is y, shifted and scaled.
            } while (em < 0.0); //Reject if in regime of zero probability.
            em=floor(em); //The trick for integer-valued distributions.
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
            //The ratio of the desired distribution to the comparison function; we accept or
            //reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
            //that t never exceeds 1.
        } while (ran1(idum) > t);
    }
    return em;
}
