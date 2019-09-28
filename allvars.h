#define OmegaLambda 0.693
#define Omega 0.307
#define h 0.677
#define c 2.9979e5	//unit in (km/s)
#define H0 67.7
#define PI 3.141592654
#define G 43007.1

float Cosmos_Luminosity_Distance(float z);
float Cosmos_Time(float z);
float luminosity_distance(float z);
float lookback_time(float z);
float integral(float low,float high,float (*func)(float x));
float Rho_crit();
float growth_factor(float z);
float H_z(float z);
float OmegaLambda_z(float z);
float Omega_z(float z);
float E_z(float z);
float k_correction(char *FilterName,float z,char *ColorName,float Color);
float calc_kcor(float Coff[][4],int Row,float z,float Color);
float absolute_magnitude(float lum,char band);
float radians(float deg);
float Cosmos_Angular_Distance(float z);
float Cosmos_Angle(float size,float z);
float galaxy_included_angle(double raA,double decA,double raB,double decB);
void stats_number_in_pro_radius(double **Data,int Row,int Col,int id,int *N);
float hms_to_deg(float h1,float m,float s);
void deg_to_hms(float d,char *hms);
float dms_to_deg(float d,float m,float s);
void deg_to_dms(float g,char *dms);
int comparefunc(const void *arg1,const void *arg2);
double median_value(const double *array,int n);
