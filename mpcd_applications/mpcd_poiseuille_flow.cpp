#include<stdlib.h>
#include<stdio.h>
#include <iostream>
#include<math.h>

#define nl cout<<endl;

using namespace std;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

int const nofp =8000,max_part=50,lx=20,ly=20;
double llx= static_cast<double>(lx), lly=static_cast<double>(ly);
double *pos,*vel,*temp_pos;
double kb_temp=1.0,mass=1.0,dt=0.1,mass_inv=1.0,force=0.01;
double vconst=sqrt(12.0*kb_temp/mass),nofpp = double(nofp);
int zzz=-4280145,zzzz=87383143,nofiter=5000,nofiter1=2000;
int part_no[lx*ly],cell_part[max_part][lx*ly];
double energy,momentum;
double avg_outlet_vel[lx*ly];
double cell_vel1[lx*ly],cell_vel2[lx*ly];
double avg_last_vel[ly];

double streaming(int step,double *pos);
void posinit(double *pos);
void velinit(double *vel, int step);
void celllist(double *temp_pos);
void printxy_avg(int step);

void collision(int step, double *vel);
double ran2(int *idum);
int main()
{
        cout.setf(std::ios_base::scientific);
        cout.precision(16);
	for (int i = 0; i < lx; ++i)
	{
		avg_last_vel[i]=0.0;
	}

	int step=0;
	double ran1,zz1;
        zz1=ran2(&zzz);// subroutine call
	
	pos=(double*)calloc(2*nofp,sizeof(double));
	vel=(double*)calloc(2*nofp,sizeof(double));
	temp_pos=(double*)calloc(2*nofp,sizeof(double));

	posinit(pos);//Subroutine call
	velinit(vel,step);//subroutine call

        FILE *fp=fopen("energy.dat","w");
       	for(step=1;step<=nofiter;step++){

		energy=streaming(step,pos);
		celllist(&temp_pos[0]);
		collision(step,vel);
		 if(step%100==0)
		 {
			printxy_avg(step);

		}   
		fprintf(fp,"  Step \t Energy \t Momentum \n");
                fprintf(fp,"  %d\t %lf\t %lf\n",step,energy,momentum);
               
	}
         
      // cout<<"_____________Avg Last  vel____________________"<<endl;
      // for (int i = 0; i < lx; i++)
      // {
      // 	cout<<19.5 <<" "<<i+0.5<<" "<<avg_last_vel[i]/(nofiter/100)<<endl;
      // }
	fclose(fp);

	return 0;
}

double returnAbsolute(double val)
{
	if (val<0)
	{
		val=-1*val;
	}
	return val;
}

// position
void posinit(double *pos)
{
	for(int i=0;i<nofp;i++)
	{
		pos[2*i]=ran2(&zzzz)*llx;
		pos[2*i+1]=ran2(&zzzz)*lly;
	}
}

//velocity
void velinit(double *vel,int step)
{
	double av_velx,av_vely;
       
	av_velx=0.0,av_vely=0.0,energy=0.0,momentum =0.0;

	for(int i=0;i<nofp;i++){
          vel[2*i] = vconst*(ran2(&zzzz)-0.5);
          av_velx = av_velx + vel[2*i];
          vel[2*i+1] =vconst*(ran2(&zzzz)-0.5);
          av_vely = av_vely + vel[2*i+1];
	}

	av_velx/=nofpp;
	av_vely/=nofpp;
       for(int i=0;i<nofp;i++){
          vel[2*i] = vel[2*i] - av_velx;
          vel[2*i+1] = vel[2*i+1] - av_vely;
                  
        }
         
	for(int i=0;i<nofp;i++)
        {
          energy += 0.5*mass*(vel[2*i]*vel[2*i]+vel[2*i+1]*vel[2*i+1]); 
          momentum = momentum + mass*sqrt(vel[2*i]*vel[2*i]+vel[2*i+1]*vel[2*i+1]);
	}
        energy/=(double)nofpp;
        momentum/=(double)nofpp;
        
}

// celllist
void celllist(double *temp_pos)
{
	int cell_no;
	double r1,r2;

	for(int i=0;i<lx*ly;i++) part_no[i]=0;

	r1=ran2(&zzzz)-0.5;
	r2=ran2(&zzzz)-0.5;
	for(int i=0;i<nofp;i++)
	{
		temp_pos[2*i]   = pos[2*i]-r1;
		temp_pos[2*i+1] = pos[2*i+1]-r2;

		if(temp_pos[2*i]<0.0)
		{
			temp_pos[2*i]=temp_pos[2*i]+llx;
		}
		else if (temp_pos[2*i]>llx)
		{
			temp_pos[2*i]=temp_pos[2*i]-llx;
		}
		if(temp_pos[2*i+1]<0.0)
		{
		     temp_pos[2*i+1]=temp_pos[2*i+1]+lly;
                                          
		}
		else if (temp_pos[2*i+1]>lly)
		{
			temp_pos[2*i+1]=temp_pos[2*i+1]-lly;
                                             
		}
	}

	for(int i=0;i<nofp;i++)
	{
		cell_no=int(temp_pos[2*i])+lx*int(temp_pos[2*i+1]);
		part_no[cell_no] = part_no[cell_no]+1;
		
		int j=part_no[cell_no];
		cell_part[j][cell_no] = i;
                
	}
}

//collision
void collision(int step, double *vel){

        double rx,ry,theta,ran1,rho;
        double rot11,rot12,rot21,rot22;
        double del_vx,del_vy,var,scale_fac;
        double del_vx1[nofp],del_vy1[nofp];

	for(int i=0;i<lx*ly;i++){
          cell_vel1[i]=0.0;
          cell_vel2[i]=0.0;
       }

	for(int i=0;i<lx*ly;i++){
          if(part_no[i]>1){
	   for(int j=1;j<=part_no[i];j++){
             int k=cell_part[j][i];
             cell_vel1[i]=cell_vel1[i]+vel[2*k];
             cell_vel2[i]=cell_vel2[i]+vel[2*k+1];	
           }

           cell_vel1[i]=cell_vel1[i]/double(part_no[i]);
           cell_vel2[i]=cell_vel2[i]/double(part_no[i]);
	
	       
         
            theta=2.0*M_PI*ran2(&zzzz);

           rot11 =  cos(theta);
           rot12 = -sin(theta);
           rot21 =  sin(theta);
           rot22 =  cos(theta);

           for(int j=1;j<=part_no[i];j++){
             int k=cell_part[j][i];
             del_vx=vel[2*k]-cell_vel1[i];
             del_vy=vel[2*k+1]-cell_vel2[i];

            vel[2*k]=cell_vel1[i]+rot11*del_vx+rot12*del_vy;
            vel[2*k+1]=cell_vel2[i]+rot21*del_vx+rot22*del_vy;

	   }
          }
        }
  
if(step%25==0){
          for(int i=0;i<lx*ly;i++){
            var=0.0;
            for(int j=1;j<=part_no[i];j++){
               int k=cell_part[j][i];
               del_vx1[k]=vel[2*k]-cell_vel1[i];
               del_vy1[k]=vel[2*k+1]-cell_vel2[i];
               var=var + del_vx1[k]*del_vx1[k] + del_vy1[k]*del_vy1[k];
            }
            scale_fac=sqrt(2.0*float(part_no[i]-1)*kb_temp*mass_inv/var);
            for(int j=1;j<part_no[i];j++){
              int k=cell_part[j][i];
              vel[2*k]   = cell_vel1[i] + del_vx1[k]*scale_fac;
              vel[2*k+1] = cell_vel2[i] + del_vy1[k]*scale_fac;
            }
          }
	}  
    }

// streaming
double streaming(int step,double *pos)
{
      
	for(int i=0;i<nofp;i++)
	{
		
                pos[2*i]=pos[2*i]+(vel[2*i]*dt)+(dt*dt*force)/(2.0*mass);
		pos[2*i+1]=pos[2*i+1]+(vel[2*i+1]*dt);

		vel[2*i]=vel[2*i]+(force*dt*mass_inv);

		if(pos[2*i]<0.0)
		{
			pos[2*i]=pos[2*i]+llx;
		}
		else if(pos[2*i]>llx)
		{
			pos[2*i]=pos[2*i]-llx;
		}

		if(pos[2*i+1]<0.0)
		{

			double dtafter=returnAbsolute(pos[2*i+1])/returnAbsolute(vel[2*i+1]);

			pos[2*i+1]=pos[2*i+1]-2*dtafter*vel[2*i+1];
			vel[2*i+1]=-1*vel[2*i+1];
		
		}
       
		else if(pos[2*i+1]>lly)
		{

			double dtafter=(returnAbsolute(pos[2*i+1])-lly)/returnAbsolute(vel[2*i+1]);

			pos[2*i+1]=pos[2*i+1]-2*dtafter*vel[2*i+1];
			vel[2*i+1]=-1*vel[2*i+1];
		}
                 
                  energy += 0.5*mass*(vel[2*i]*vel[2*i]+vel[2*i+1]*vel[2*i+1]);
                  momentum +=  mass*sqrt(vel[2*i]*vel[2*i]+vel[2*i+1]*vel[2*i+1]);

	     } 
                //cout<<step<<" "<<energy<< " " << momentum<<endl;
	energy = energy/(double)nofpp;
        momentum = momentum/(double)nofpp;

	return energy;
}

//RAN1 returns a unifom random deviate on the interval [0,1]
double ran2(int *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
void printxy_avg(int step)
{
	//cout << "step" <<step<< endl;
	//cout <<"-------------------"<<endl;
	for (int i = 0; i < lx; i++)
	{
		for (int j = 0; j < ly; j++)
		{
			if (i==lx-1)
			{
				cout << i+0.5<<" "<<j+0.5<< "  "<< cell_vel1[i*lx+j] << "  "<<cell_vel2[i*ly+j]<< endl;
				avg_last_vel[j]+=cell_vel1[i*lx+j];
			}
		}
	}
	//cout << "------------------------------"<<endl;
}
