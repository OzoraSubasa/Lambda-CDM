#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include<stdbool.h>
#include"allvars.h"

float E_z(float z)
{
    float zplus1;
    zplus1=1+z;
    return sqrt(OmegaLambda+(1-Omega-OmegaLambda)*pow(zplus1,2.0)+Omega*pow(zplus1,3.0));
}

float Omega_z(float z)
{
    return Omega*pow(1+z,3.0)/E_z(z)/E_z(z);
}

float OmegaLambda_z(float z)
{
    return OmegaLambda/E_z(z)/E_z(z);
}

float H_z(float z)
{
    return H0*E_z(z);
}

float growth_factor(float z)
{
    float g,omegaz,omegalambdaz;
    omegaz=Omega_z(z);
    omegalambdaz=OmegaLambda_z(z);
    g=5./2.*omegaz/(pow(omegaz,4./7.)-omegalambdaz+(1+omegaz/2.)*(1+omegalambdaz/70.));
    return g/(1+z);
}

float Rho_crit()
{
    return 3*pow(H0,2.0)/8/PI/G;	//27.75505 (h^2.10^10Msun.Mpc^-3)
}

float integral(float low,float high,float (*func)(float x))
{
    int i,n=10000;
    float interval,I=0.0,temp;
    if (high<low)
    {
	temp=low;
	low=high;
	high=temp;
    }
    if (high==low)
	return I;
    temp=0.0;
    interval=(high-low)*1.0/n;
    for (i=1;i<=n-1;i++)
	temp+=(*func)(low+i*interval);
    I=interval*(0.5*(*func)(low)+temp+0.5*(*func)(high));
    return I;
}

float lookback_time(float z)
{
    float t;
    t=1.0/((1+z)*sqrt(OmegaLambda+Omega*pow(1+z,3.0)));
    return t;
}

float luminosity_distance(float z)
{
    float Dl;
    Dl=c/H0/sqrt(OmegaLambda+Omega*pow(1+z,3.0));
    return Dl;
}

float Cosmos_Time(float z)
{
    float t;
    t=integral(0,z,lookback_time);
    return t;
}

float Cosmos_Luminosity_Distance(float z)
{
    float DL;
    DL=(1+z)*integral(0,z,luminosity_distance);
    return DL;
}

float Cosmos_Angular_Distance(float z)
{
    return Cosmos_Luminosity_Distance(z)/pow(1+z,2.0);
}

float k_correction(char *FilterName,float z,char *ColorName,float Color)
{
    /*SDSS K-correction calculator in C.See http://kcor.sai.msu.ru for the 
     * reference.*/
    char buf[100],options[100];
    int i,j,Row;
    for (i=0,j=0;i<strlen(ColorName);i++)
    {
	if (ColorName[i]=='-')
	    continue;
	buf[j++]=ColorName[i];
    }
    buf[j]='\0';
    sprintf(options,"%s_%s",FilterName,buf);
//    printf("%s\n",options);
    if (strcmp(options,"g_gi")==0)
    {
	float Coff[][4]={{0,0,0,0},{1.59269,-2.97991,7.31089,-3.46913},
	{-27.5631,-9.89034,15.4693,6.53131},{161.969,-76.171,-56.1923,0},
	{-204.457,217.977,0,0},{-50.6269,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"g_gz")==0)
    {
	float Coff[][4]={{0,0,0,0},{2.37454,-4.39943,7.29383,-2.90691},
	{-28.7217,-20.7783,18.3055,5.04468},{220.097,-81.883,-55.8349,0},
	{-290.86,253.677,0,0},{-73.5316,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"g_gr")==0)
    {
	float Coff[][4]={{0,0,0,0},{-2.45204,4.10188,10.5258,-13.5889},
	{56.7969,-140.913,144.572,57.2155},{-466.949,222.789,-917.46,-78.0591},
	{2906.77,1500.8,1689.97,30.889},{-10453.7,-4419.56,-1011.01,0},
	{17568,3236.68,0,0},{-10820.7,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"i_gi")==0)
    {
	float Coff[][4]={{0,0,0,0},{-2.21853,3.94007,0.678402,-1.24751},
	{-15.7929,-19.3587,15.0137,2.27779},{118.791,-40.0709,-30.6727,0},
	{-134.571,125.799,0,0},{-55.4483,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"i_ui")==0)
    {
	float Coff[][4]={{0,0,0,0},{-3.91949,3.20431,-0.431124,-0.000912813},
	{-14.776,-6.56405,1.15975,0.0429679},{135.273,-1.30583,-1.81687,0},
	{-264.69,15.2846,0,0},{142.624,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"r_gr")==0)
    {
	float Coff[][4]={{0,0,0,0},{1.83285,-2.71446,4.97336,-3.66864},
	{-19.7595,10.5033,18.8196,6.07785},{33.6059,-120.713,-49.299,0},
	{144.371,216.453,0,0},{-295.39,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"r_ur")==0)
    {
	float Coff[][4]={{0,0,0,0},{3.03458,-1.50775,0.576228,-0.0754155},
	{-47.8362,19.0053,-3.15116,0.286009},{154.986,-35.6633,1.09562,0},
	{-188.094,28.1876,0,0},{68.9867,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"u_ur")==0)
    {
	float Coff[][4]={{0,0,0,0},{10.3686,-6.12658,2.58748,-0.299322},
	{-138.069,45.0511,-10.8074,0.95854},{540.494,-43.7644,3.84259,0},
	{-1005.28,10.9763,0,0},{710.482,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"u_ui")==0)
    {
	float Coff[][4]={{0,0,0,0},{11.0679,-6.43368,2.4874,-0.276358},
	{-134.36,36.0764,-8.06881,0.788515},{528.447,-26.7358,0.324884,0},
	{-1023.1,13.8118,0,0},{721.096,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"u_uz")==0)
    {
	float Coff[][4]={{0,0,0,0},{11.9853,-6.71644,2.31366,-0.234388},
	{-137.024,35.7475,-7.48653,0.655665},{519.365,-20.9797,0.670477,0},
	{-1028.36,2.79717,0,0},{767.552,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"z_gz")==0)
    {
	float Coff[][4]={{0,0,0,0},{0.30146,-0.623614,1.40008,-0.534053},
	{-10.9584,-4.515,2.17456,0.913877},{66.0541,4.18323,-8.42098,0},
	{-169.494,14.5628,0,0},{144.021,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"z_rz")==0)
    {
	float Coff[][4]={{0,0,0,0},{0.669031,-3.08016,9.87081,-7.07135},
	{-18.6165,8.24314,-14.2716,13.8663},{94.1113,11.2971,-11.9588,0},
	{-225.428,-17.8509,0,0},{197.505,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else if (strcmp(options,"z_uz")==0)
    {
	float Coff[][4]={{0,0,0,0},{0.623441,-0.293199,0.16293,-0.0134639},
	{-21.567,5.93194,-1.41235,0.0714143},{82.8481,-0.245694,0.849976,0},
	{-185.812,-7.9729,0,0},{168.691,0,0,0}};
	Row=sizeof(Coff)/sizeof(float)/4;
	return calc_kcor(Coff,Row,z,Color);
    }
    else
    {
	fprintf(stderr,"options has no '%s'.\n",options);
	exit(10);
    }
}

float calc_kcor(float Coff[][4],int Row,float z,float Color)
{
    int i,j;
    float Kcor=0.0;

//    printf("Row=%d\n",Row);
    for (i=0;i<Row;i++)
    {
	for (j=0;j<4;j++)
	    Kcor+=Coff[i][j]*pow(z,i)*pow(Color,j);
    }
    return Kcor;
}

float absolute_magnitude(float lum,char band)
{
    float SunV=4.83;
    float SunR=SunV-0.36,SunB=SunV+0.65,SunI=SunV-0.72;
    float SunU=SunB+0.14;

    if (band=='r')
	return -2.5*lum+SunR;
    else if (band=='u')
	return -2.5*lum+SunU;
    else if (band=='i')
	return -2.5*lum+SunI;
    else 
    {
	fprintf(stderr,"band has no '%c'.\n",band);
	exit(20);
    }
}

float radians(float deg)
{
    return deg/180.0*PI;
}

float Cosmos_Angle(float size,float z)
{
    return size/Cosmos_Angular_Distance(z);
}

float galaxy_included_angle(double raA,double decA,double raB,double decB)
{
    float angle;

    decA=radians(decA);
    raA=radians(raA);
    decB=radians(decB);
    raB=radians(raB);
    angle=acos(cos(decA)*cos(decB)*cos(raA-raB)+sin(decA)*sin(decB));
    return angle;
}

float hms_to_deg(float h1,float m,float s)
{
    float deg;

    if (h1<0.0 || h1>=24.0)
    {
	fprintf(stderr,"h must is between 0.0 and 24.0.\n");
	exit(1);
    }
    if (m<0.0 || m>=60.0)
    {
	fprintf(stderr,"m must is between 0.0 and 60.0.\n");
	exit(2);
    }
    if (s<0.0 || s>=60.0)
    {
	fprintf(stderr,"s must is between 0.0 and 60.0.\n");
	exit(3);
    }
    deg=h1*15.0+m*0.25+s*(0.25/60.0);
    return deg;
}
void deg_to_hms(float d,char *hms)
{
    int h1,m;
    float s;
    if (d<0.0 || d>=360.0)
    {
	fprintf(stderr,"d must is between 0.0 and 360.0.\n");
	exit(4);
    }
    h1=(int)(d/15.0);
    d-=h1*15.0;
    m=(int)(d/0.25);
    d-=m*0.25;
    s=d/(0.25/60.0);
    sprintf(hms,"%d:%d:%f",h1,m,s);
}
float dms_to_deg(float d,float m,float s)
{
    int sign=1;
    float deg;

    if (m<0.0 || m>60.0)
    {
	fprintf(stderr,"m must is between 0.0 and 60.0.\n");
	exit(5);
    }
    if (s<0.0 || s>=60.0)
    {
	fprintf(stderr,"s must is between 0.0 and 60.0.\n");
	exit(6);
    }
    if (d<0.0)
	sign=-1;
    deg=d+sign*(m/60.0)+sign*(s/3600.0);
    return deg;
}

void deg_to_dms(float g,char *dms)
{
    int d,m,sign=1;
    float s;
    char symbol;
    if (g<0.0)
	sign=-1;
    g=fabs(g);
    d=(int)(g)*sign;
    g-=(int)(g);
    m=(int)(g*60.0);
    g-=m/60.0;
    s=g*3600.0;
    switch (sign)
    {
	case -1:symbol='-';break;
	case 1:symbol='+';break;
    }
    sprintf(dms,"%c%d:%d:%f",symbol,d,m,s);
}

int comparefunc(const void *arg1,const void *arg2)
{
    const double *ptr1=(const double *)arg1;
    const double *ptr2=(const double *)arg2;
    double val1=*ptr1;
    double val2=*ptr2;
    if (val1<val2)
	return -1;
    if (val1==val2)
	return 0;
    return 1;
}

double median_value(const double *array,int n)
{
    double *temp;
    double median;
    int i;

    temp=(double *)malloc(sizeof(double)*n);
    assert(temp);
    for (i=0;i<n;i++)
	temp[i]=array[i];
    qsort(temp,n,sizeof(double),comparefunc);
    for (i=0;i<n;i++)
	printf("%lf ",temp[i]);
    printf("\n");
    if (n%2==0)
    {
	median=(temp[n/2-1]+temp[n/2])/2;
	free(temp);
	return median;
    }
    else
    {
	median=temp[(n+1)/2-1];
	free(temp);
	return median;
    }
}
