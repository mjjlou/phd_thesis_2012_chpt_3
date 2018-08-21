#include <math.h>

// A random deviate is a transformation (to a desired distribution) of
// the generated random number.
// Consider discrete events occurring in continuous time which occur
// with a known average rate. The exponential distribution expresses
// the waiting time between events. It is a function of a non-negative
// real number. It depends on the average rate of events, $\lambda$.

float expdev(long *idum) {
	float ran1(long *idum);
	float dum;

	do dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}
