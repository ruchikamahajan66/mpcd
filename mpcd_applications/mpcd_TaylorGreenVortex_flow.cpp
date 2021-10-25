#include<stdlib.h>
#include<stdio.h>
#include <iostream>
#include<math.h>

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

int const nofp =80000,max_part=40,lx=20,ly=20,lz=20;
double llx= static_cast<double>(lx), lly=static_cast<double>(ly),llz=static_cast<double>(lz);
double *pos,*vel,*temp_pos;
double kb_temp=1.0,mass=1.0,dt=0.1,force=0.00348,mass_inv=1.0,U=1.0,N=10.0,a=1.0;
double vconst=sqrt(12.0*kb_temp/mass),nofpp = double(nofp);
double llxby2=llx/2.0,llyby2=lly/2.0,llzby2=llz/2.0;
int zzz=-4280145,zzzz=87383143,nofiter=5000,nofiter1=2000;
int part_no[lx*ly*lz],cell_part[max_part][lx*ly*lz];
double energy,momentum,collviscosity,theta,fenergy;


double streaming(int step,double *pos);
void posinit(double *pos);
void velinit(double *vel, int step);
void celllist(double *temp_pos);
void collision(int step, double *vel, double *theta);
double ran2(int *idum);
double viscosity(int step, double *theta,double *collviscosity);
double velfinal(double *vel,int step,double *collviscosity);
int main()
{
        std::cout.setf(std::ios_base::scientific);
        std::cout.precision(16);
	int step=0;
	double ran1,zz1;
        zz1=ran2(&zzz);// subroutine call
         // printf("%lf\n",zz1);
	pos=(double*)calloc(3*nofp,sizeof(double));
	vel=(double*)calloc(3*nofp,sizeof(double));
	temp_pos=(double*)calloc(3*nofp,sizeof(double));

	posinit(pos);//Subroutine call
	velinit(vel,step);//subroutine call

        FILE *fp=fopen("energy.dat","w");
        FILE *np=fopen("viscosity.dat","w");
         FILE *qp=fopen("energyf.dat","w");
       	for(step=1;step<=nofiter;step++){

		energy=streaming(step,pos);
		celllist(&temp_pos[0]);//subroutine call
		collision(step,vel,&theta);//subroutine call
                 viscosity(step,&theta,&collviscosity);
                 velfinal(vel,step,&collviscosity);
                //printf("%d\n",step);
                fprintf(fp,"%d %lf %lf\n",step,energy,momentum);
                fprintf(np,"%d\t %lf\n",step,collviscosity);
              fprintf(qp,"%d\t %lf\n",step,fenergy);   
                              
	}
	fclose(fp);
        fclose(np);
        fclose(qp);
	return 0;
}

// position
void posinit(double *pos)
{
	for(int i=0;i<nofp;i++)
	{
		pos[3*i]=ran2(&zzzz)*llx;
		pos[3*i+1]=ran2(&zzzz)*lly;
		pos[3*i+2]=ran2(&zzzz)*llz;
	}
}

//velocity
void velinit(double *vel,int step)
{
	double av_velx,av_vely,av_velz,energy;
	av_velx=0.0,av_vely=0.0,av_velz=0.0,energy=0.0;

        FILE *fp=fopen("velocity.dat","w");

	for(int i=0;i<nofp;i++){
          vel[3*i] = U*sin(2*M_PI*pos[3*i])*cos(2*M_PI*pos[3*i+1])*cos(2*M_PI*pos[3*i+2]);
          //vel[3*i] = vconst*(ran2(&zzzz)-0.5);
          av_velx = av_velx + vel[3*i];
          vel[3*i+1] = -U*cos(2*M_PI*pos[3*i])*sin(2*M_PI*pos[3*i+1])*cos(2*M_PI*pos[3*i+2]);  
         //vel[3*i+1] = vconst*(ran2(&zzzz)-0.5);
          av_vely = av_vely + vel[3*i+1];
          //vel[3*i+2] = vconst*(ran2(&zzzz)-0.5);
          vel[3*i+2] = 0.0;
          av_velz = av_velz + vel[3*i+2];

         fprintf(fp," %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", pos[3*i],pos[3*i+1],pos[3*i+2],vel[3*i],vel[3*i+1],vel[3*i+2]);
      // std::cout <<  vel[3*i]<< "   " << vel[3*i+1] << "  " << vel[3*i+2] << "   " <<  std::endl;
	}
       fclose(fp);

	av_velx/=nofpp;
	av_vely/=nofpp;
	av_velz/=nofpp;
	for(int i=0;i<nofp;i++){
          vel[3*i] = vel[3*i] - av_velx;
          vel[3*i+1] = vel[3*i+1] - av_vely;
          vel[3*i+2] = vel[3*i+2] - av_velz;
          energy += 0.5*mass*(vel[3*i]*vel[3*i]+vel[3*i+1]*vel[3*i+1]+vel[3*i+2]*vel[3*i+2]);
	}
        energy/=(double)nofpp;
     // std::cout << step << "  " << energy << std::endl;
}

// celllist
void celllist(double *temp_pos)
{
	int cell_no;
	double r1,r2,r3;

	for(int i=0;i<lx*ly*lz;i++) part_no[i]=0;

	r1=ran2(&zzzz)-0.5;
	r2=ran2(&zzzz)-0.5;
	r3=ran2(&zzzz)-0.5;
	for(int i=0;i<nofp;i++)
	{
		temp_pos[3*i]   = pos[3*i]-r1;
		temp_pos[3*i+1] = pos[3*i+1]-r2;
		temp_pos[3*i+2] = pos[3*i+2]-r3;

		if(temp_pos[3*i]<0.0)
		{
			temp_pos[3*i]=temp_pos[3*i]+llx;
		}
		else if (temp_pos[3*i]>llx)
		{
			temp_pos[3*i]=temp_pos[3*i]-llx;
		}
		if(temp_pos[3*i+1]<0.0)
		{
			temp_pos[3*i+1]=temp_pos[3*i+1]+lly;
		}
		else if (temp_pos[3*i+1]>lly)
		{
			temp_pos[3*i+1]=temp_pos[3*i+1]-lly;
		}

		if(temp_pos[3*i+2]<0.0)
		{
			temp_pos[3*i+2]=temp_pos[3*i+2]+llz;
		}
		else if (temp_pos[3*i+2]>llz)
		{
			temp_pos[3*i+2]=temp_pos[3*i+2]-llz;
		}
                //printf("%lf\t %lf\t %lf\n",temp_pos[3*i],temp_pos[3*i+1], temp_pos[3*i+2]);
	}
               FILE *fp=fopen("cellparticle1.dat","w");

	for(int i=0;i<nofp;i++)
	{
		cell_no=int(temp_pos[3*i])+lx*int(temp_pos[3*i+1])+lx*ly*int(temp_pos[3*i+2]);
		part_no[cell_no] = part_no[cell_no]+1;
		int j=part_no[cell_no];
		cell_part[j][cell_no] = i;
              // std::cout << step << "   " << i << "  " << j << "   " << cell_no << "   " << cell_part[j][cell_no] << std::endl;
                fprintf(fp,"%d\n",j);
	}
                fclose(fp);
}

//collision
void collision(int step, double *vel, double *theta){

        double cell_vel1[lx*ly*lz],cell_vel2[lx*ly*lz],cell_vel3[lx*ly*lz];
        double rx,ry,rz,ran1,phi,rho;
        double rot11,rot12,rot13,rot21,rot22,rot23,rot31,rot32,rot33;
        double del_vx,del_vy,del_vz,var,scale_fac;
        double del_vx1[nofp],del_vy1[nofp],del_vz1[nofp];

	for(int i=0;i<lx*ly*lz;i++){
          cell_vel1[i]=0.0;
          cell_vel2[i]=0.0;
          cell_vel3[i]=0.0;
       }

	for(int i=0;i<lx*ly*lz;i++){
          if(part_no[i]>1){
	   for(int j=1;j<=part_no[i];j++){
             int k=cell_part[j][i];
             cell_vel1[i]=cell_vel1[i]+vel[3*k];
             cell_vel2[i]=cell_vel2[i]+vel[3*k+1];
             cell_vel3[i]=cell_vel3[i]+vel[3*k+2];	
           }
           cell_vel1[i]=cell_vel1[i]/double(part_no[i]);
           cell_vel2[i]=cell_vel2[i]/double(part_no[i]);
           cell_vel3[i]=cell_vel3[i]/double(part_no[i]);

           rho=2.0*ran2(&zzzz)-1.0;
           phi=4.0*asin(1.0)*ran2(&zzzz);
           rx=cos(phi)*sqrt(1-rho*rho);
           ry=sin(phi)*sqrt(1-rho*rho);
           rz=rho;
           *theta=2.0*asin(1.0)*130.0/180.0;
           // std::cout << rho  << "  " << theta<< "  " << phi << std::endl;

           rot11 = (1.0-cos(*theta)) *rx*rx + cos(*theta);
           rot12 = (1.0-cos(*theta)) *rx*ry - sin(*theta)*rz;
           rot13 = (1.0-cos(*theta)) *rx*rz + sin(*theta)*ry;
           rot21 = (1.0-cos(*theta)) *ry*rx + sin(*theta)*rz;
           rot22 = (1.0-cos(*theta)) *ry*ry + cos(*theta);
           rot23 = (1.0-cos(*theta)) *ry*rz - sin(*theta)*rx;
           rot31 = (1.0-cos(*theta)) *rz*rx - sin(*theta)*ry;
           rot32 = (1.0-cos(*theta)) *rz*ry + sin(*theta)*rx;
           rot33 = (1.0-cos(*theta)) *rz*rz + cos(*theta);
//std::cout << rot11<< "   " << rot12 << "  " << rot21 << "   " << rot22 << "   " << std::endl;

           for(int j=1;j<=part_no[i];j++){
             int k=cell_part[j][i];
             del_vx=vel[3*k]-cell_vel1[i];
             del_vy=vel[3*k+1]-cell_vel2[i];
             del_vz=vel[3*k+2]-cell_vel3[i];
   //  std::cout << del_vx << "   " << del_vy << "  " << std::endl;

             vel[3*k]=cell_vel1[i]+rot11*del_vx+rot12*del_vy+rot13*del_vz;
             vel[3*k+1]=cell_vel2[i]+rot21*del_vx+rot22*del_vy+rot23*del_vz;
             vel[3*k+2]=cell_vel3[i]+rot31*del_vx+rot32*del_vy+rot33*del_vz;
    //std::cout << vel[3*k] << "   " << vel[3*k+1]<< "  " << std::endl;

	   }
          }// if
        }//i

	/*if(step%25==0){
          for(int i=0;i<lx*ly*lz;i++){
            var=0.0;
            for(int j=1;j<=part_no[i];j++){
               int k=cell_part[j][i];
               del_vx1[k]=vel[3*k]-cell_vel1[i];
               del_vy1[k]=vel[3*k+1]-cell_vel2[i];
               del_vz1[k]=vel[3*k+2]-cell_vel3[i];
               var=var + del_vx1[k]*del_vx1[k] + del_vy1[k]*del_vy1[k] + del_vz1[k]*del_vz1[k];
            }
            scale_fac=sqrt(3.0*float(part_no[i]-1)*kb_temp*mass_inv/var);
            for(int j=1;j<part_no[i];j++){
              int k=cell_part[j][i];
              vel[3*k]   = cell_vel1[i] + del_vx1[k]*scale_fac;
              vel[3*k+1] = cell_vel2[i] + del_vy1[k]*scale_fac;
              vel[3*k+2] = cell_vel3[i] + del_vz1[k]*scale_fac;
            }
          }//for
	}//if*/
}
// streaming
double streaming(int step,double *pos)
{
	energy = 0.0;
        momentum = 0.0;

	for(int i=0;i<nofp;i++)
	{
		
                pos[3*i]=pos[3*i]+vel[3*i]*dt;
		pos[3*i+1]=pos[3*i+1]+vel[3*i+1]*dt;
		pos[3*i+2]=pos[3*i+2]+vel[3*i+2]*dt;

		if(pos[3*i]<0.0)
		{
			pos[3*i]=pos[3*i]+llx;
		}
		else if(pos[3*i]>llx)
		{
			pos[3*i]=pos[3*i]-llx;
		}

		if(pos[3*i+1]<0.0)
		{
			pos[3*i+1]=pos[3*i+1]+lly;
		}
		else if(pos[3*i+1]>lly)
		{
			pos[3*i+1]=pos[3*i+1]-lly;
		}

		if(pos[3*i+2]<0.0)
		{
			pos[3*i+2]=pos[3*i+2]+llz;
		}
		else if(pos[3*i+2]>llz)
		{
			pos[3*i+2]=pos[3*i+2]-llz;
		}
                energy +=0.5*mass*(vel[3*i]*vel[3*i]+vel[3*i+1]*vel[3*i+1]+ vel[3*i+2]*vel[3*i+2]);
                momentum = momentum + mass*sqrt(vel[3*i]*vel[3*i]+vel[3*i+1]*vel[3*i+1]+ vel[3*i+2]*vel[3*i+2]);

	}
	        // Per particle Energy
	        energy/=(double)nofpp;
               // std::cout << step  << "  " << momentum<< "  " << energy << std::endl;
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
double viscosity(int step, double *theta,double *collviscosity)
{
        *collviscosity = 0.0;
        double mult_fac,div_fac;
               mult_fac = N - 1 + exp(-N), div_fac = 18*a*dt;
        *collviscosity = *collviscosity + mass*(1-cos(*theta))* mult_fac/div_fac;
// printf("%lf\n", *collviscosity);
 return *collviscosity;
}
double velfinal(double *vel,int step,double *collviscosity)
{
 //double fenergy;
	fenergy=0.0;

	for(int i=0;i<nofp;i++){
          vel[3*i] = U*exp(*collviscosity*dt)*sin(2*M_PI*pos[3*i])*cos(2*M_PI*pos[3*i+1])*cos(2*M_PI*pos[3*i+2]);
          vel[3*i+1] = -U*exp(*collviscosity*dt)*cos(2*M_PI*pos[3*i])*sin(2*M_PI*pos[3*i+1])*cos(2*M_PI*pos[3*i+2]);  
          vel[3*i+2] = 0.0;
           fenergy +=0.5*mass*(vel[3*i]*vel[3*i]+vel[3*i+1]*vel[3*i+1]+ vel[3*i+2]*vel[3*i+2]);
      // std::cout << vel[3*i]<< "   " << vel[3*i+1] << "  " << vel[3*i+2] << "   " <<  std::endl;
	}
          fenergy/=(double)nofpp;   
          printf("%d\t %lf\n",step,fenergy);
           
return   fenergy;
}
