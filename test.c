#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<time.h>
#include"allvars.h"

int main(int argc,char **argv)
{
    float z,dl,dA,rad,size,theta,mu,mr,Kcor,Mr;
    float IncludedAngle,ra,dec;
    double *array,median;
    int n,i;
    char hms[20],dms[20];

    if (argc!=3)
    {
	printf("please enter the redshift and number of array.\n");
	printf("./a.out 0.06 10\n");
	exit(1);
    }

    z=atof(argv[1]);
    mu=19.04;
    mr=16.25;
    
    /*distance*/ 
    printf("*******distance*********\n");
    dl=Cosmos_Luminosity_Distance(z);	//Mpc
    printf("dl=%f\n",dl);	
    dA=Cosmos_Angular_Distance(z);	//Mpc
    printf("dA=%f\n",dA);
    rad=radians(1./3600);	//unit in radian
    size=rad*dA*1.0e3;	//kpc
    printf("size=%f\n",size);
    theta=Cosmos_Angle(size/1.0e3,z);
    printf("rad=%e\ttheta=%e\n",rad,theta);
    IncludedAngle=galaxy_included_angle(136.810196,52.062061,136.953033,52.095886);
    printf("angle=%f\n",IncludedAngle);
    
    /*absolute magnitude*/
    printf("*******absolute magnitude**********\n");
    Kcor=k_correction("r",z,"u-r",mu-mr);
    Mr=mr-5*(log10(dl*1.0e6)-1.0)-Kcor;
    printf("Mr=%f\n",Mr);
    Mr=absolute_magnitude(10.098123,'r');
    printf("Mr=%f\n",Mr);
    
    /*K correction*/
    printf("*******K correction********\n");
    printf("%f\n",k_correction("g",0.2,"g-r",1.1));

    /*deg to hms or hms to deg or deg to dms or dms to deg*/
    printf("******convert between decimal and sexagesimal*******\n");
    ra=hms_to_deg(9,7,40.05);
    dec=dms_to_deg(56,4,13.85);
    printf("ra=%f,dec=%f\n",ra,dec);
    deg_to_hms(136.916907919,hms);
    deg_to_dms(56.070515026,dms);
    printf("%s %s\n",hms,dms);
    
    /*median value*/
    printf("******median value**********\n");
    srand(time(NULL));
    n=atoi(argv[2]);
    array=(double *)malloc(sizeof(double)*n);
    for (i=0;i<n;i++)
	array[i]=rand()%10+1;
    median=median_value(array,n);
    printf("median=%lf\n",median);
    
    return 0;
}
