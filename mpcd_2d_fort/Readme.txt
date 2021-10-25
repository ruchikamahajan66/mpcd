There are 6 subroutines to implement MPCD in 2D

1. double streaming(int step,double *pos);
	used to stream the particles with time 
2. void posinit(double *pos);
	it initializes the particles position using ran2(int *idum) subroutine.
3. void velinit(double *vel, int step);
	velocity initialization using ran2(int *idum) subroutine.
4. void celllist(double *temp_pos);
	It divided the simulation box in the number of cells
5. void collision(int step, double *vel);
	This subroutine implements the logic for elastic collision.
6. double ran2(int *idum);
	It gives us the unformly distributed random number.


The above all is for general MPCD method.
 	In the streaming subroutine I have applied the periodic boundary condition in x direction and bounce back condition in the y direction. And I have applied force to the particles in x-y direction.