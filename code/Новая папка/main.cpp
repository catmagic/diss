#include <iostream>
#include<stdlib.h>
#include <fstream>
#include <boost/multiprecision/float128.hpp>
#define _USE_MATH_DEFINES
#include<math.h>
#include <cmath>
#include <complex>
#include<string>
#include <boost/multiprecision/complex128.hpp>
#include"quadmath.h"
using namespace std;
using namespace boost::multiprecision;
#define  EPSILON 1e-16Q
float128 convert_F_to_K(float128 T_f)
{
    return (T_f-32.0)*5.0/9.0+273.15;
}

float128 convert_psi_to_bar(float128 P_psi)
{
    return P_psi/14.504;
}

void   normalization_P(float128* P,float128 * P_crit,float128 *P_r ,int n)
{
    for(int i=0; i<n ;i++)
    {
        P_r[i]=P[i]/P_crit[i];
    }
}

void   normalization_T(float128* T,float128 * T_crit,float128 *T_r ,int n)
{
    for(int i=0; i<n ;i++)
    {
        T_r[i]=T[i]/T_crit[i];
    }
}
void Omega_a_calc(float128 Omega_a0,float128* omega,float128 * Omega_a,float128* T_r, int n)
{
    float128 temp_calc;
    for(int i=0; i<n;i++)
    {
		//temp_calc=1+(0.379642+1.48503*omega[i]-0.164423*omega[i]*omega[i]+0.016666*omega[i]*omega[i]*omega[i])*(1-sqrt(T_r[i]));
        temp_calc=1+(0.37464+1.54226*omega[i]-0.26992*omega[i]*omega[i])*(1-sqrtq((__float128)T_r[i]));
        Omega_a[i]=Omega_a0*temp_calc*temp_calc;
	//	cout<<Omega_a[i]<<'\n';
    }
}
void Omega_b_calc(float128 Omega_b0,float128* Omega_b, int n)
{
    for(int i=0; i<n;i++)
    {
        Omega_b[i]=Omega_b0;
    }
}
void calc_A(float128** A_i_j ,float128 *c,float128 &A,int n)
{
    A=0;
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n;j++)
        {
            A+=A_i_j[i][j]*c[i]*c[j];
        }
    }

}
void calc_B(float128* B_i ,float128 *c,float128 &B,int n)
{
    B=0;
    for(int i=0; i<n;i++)
    {
        B+=B_i[i]*c[i];
    }
}
void calc_A_i(float128* Omega_A,float128* P_r,float128* T_r,float128 *A_i ,int n)
{
    for(int i=0; i<n;i++)
    {
        A_i[i]=Omega_A[i]*P_r[i]/(T_r[i]*T_r[i]);
	//	cout<<A_i[i]<<'\n';
    }
}
void calc_B_i(float128* Omega_B,float128* P_r,float128* T_r,float128 *B_i ,int n)
{
    for(int i=0; i<n;i++)
    {
        B_i[i]=Omega_B[i]*P_r[i]/T_r[i];
    }
}
void calc_A_i_j(float128** Betta,float128* A_i,float128** A_i_j ,int n)
{
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n;j++)
        {
            A_i_j[i][j]=(1-Betta[i][j])*sqrtq((__float128)(A_i[i]*A_i[j]));
            //cout<<A_i_j[i][j]<<"    "; //A_i[i]<<"    "<<A_i[j]<<"     "<<Betta[i][j]<<'\n';
			
        }
		//cout<<'\n';
    }int temp;
		//	cin>>temp;
}
void calc_V(float128 *k_i,float128 *z,float128 &v,int n)
{
    float128 k_max=k_i[0],k_min=k_i[0],I_max,I_min,l,r;
    for(int i=0;i<n;i++)
    {
        if(k_max<k_i[i]){k_max=k_i[i];}
        if(k_min>k_i[i]){k_min=k_i[i];}
    }
    I_min=1/(1-k_max);//+infinity
    I_max=1/(1-k_min);//-infinity
	//cout<<I_min<<" "<<I_max<<" "<<k_min<<" "<<k_max<<'\n';
    //if(k_min>1.0)
    l=I_min;
    r=I_max;
    float128 sum,temp;

    if(k_min>=1.0)
    {
        r=1.0;

        do{
        sum=0.0;
//cin>>temp;
 //cout<<l<<"     "<<r<<"   "<<sum<<'\n';
        for(int i=0;i<n;i++)
        {
            sum+=((k_i[i]-1)*z[i])/(1+(k_i[i]-1)*r);
        }
        r=r*10;


        }while(sum>EPSILON/2.0);
       // cout<<l<<"     "<<r<<"   "<<sum<<'\n';
	v=r;
	return;
    }
if(k_max<=1.0)
    {
        l=-1.0;

        do{
        sum=0.0;
//cin>>temp;
// cout<<l<<"     "<<r<<"   "<<sum<<'\n';
        for(int i=0;i<n;i++)
        {
            sum+=((k_i[i]-1)*z[i])/(1+(k_i[i]-1)*l);
        }
        l=l*10;


        }while(sum<-EPSILON/2.0);
	v=l;
	return;
       // cout<<l<<"     "<<r<<"   "<<sum<<'\n';

    }

    v=(l+r)/2.0;
    do{

       sum=0.0;
        for(int i=0;i<n;i++)
        {
            sum+=((k_i[i]-1)*z[i])/(1+(k_i[i]-1)*v);
        }
        if(10000000*abs(sum)<EPSILON){/* cout<<"L"<<l<<" M"<<v<<" R"<<r<<" sum"<<k_min<<'\n';*/return;}
//cout<<"L"<<l<<" M"<<v<<" R"<<r<<" sum"<<k_min<<'\n';
    if(sum>0.0){l=v;v=(l+r)/2.0;}
    else {r=v;v=(l+r)/2.0;}
  // cout<<"L"<<l<<" M"<<v<<" R"<<r<<" sum"<<k_min<<'\n';
//cin>>temp;

    }while((10000*abs(sum)>EPSILON)&&(abs(r-l)>EPSILON));

}
void calc_x(float128 *k_i,float128 *z,float128*x, float128 v,int n)
{
    for(int i=0;i<n;i++)
    {
        x[i]=z[i]/(1+(k_i[i]-1)*v);
    }

}
void calc_y(float128 *k_i,float128 *z,float128*y, float128 v,int n)
{
    for(int i=0;i<n;i++)
    {
        y[i]=z[i]*k_i[i]/(1+(k_i[i]-1)*v);
    }
}
void cube_root(complex128 a,complex128 root[3])
{

 float128 angle=arg(a),radius=abs(a),rotation=1;
 if(angle<0.0){rotation=-1;}
 //std::complex<__float128> temp_complex;
 //temp_complex=polar((__float128)cbrtq((__float128)radius),(__float128)angle/3.0Q);
 float128 r=cbrtq((__float128)radius);
 float128 ang=angle/3.0Q;
 float128 ang_plus_rot=(angle+rotation*2*M_PIq)/3.0Q;
 float128 ang_mines_rot=(angle-rotation*2*M_PIq)/3.0Q;
 
 root[0]=complex128 (r*cosq((__float128)ang),r*sinq((__float128)ang));
 root[1]=complex128 (r*cosq((__float128)ang_plus_rot),r*sinq((__float128)ang_plus_rot));
 root[2]=complex128 (r*cosq((__float128)ang_mines_rot),r*sinq((__float128)ang_mines_rot));
 /*temp_complex=polar((__float128)cbrtq((__float128)radius),(__float128)(angle+rotation*2*M_PIq)/3.0Q);
 root[1]=complex128 (temp_complex.real(),temp_complex.imag());
 temp_complex=polar((__float128)cbrtq((__float128)radius),(__float128)(angle-rotation*2*M_PIq)/3.0Q);
 root[2]=complex128 (temp_complex.real(),temp_complex.imag());*/
}
void normalization_c(float128 *c,int n)
{
    float128 sum=0.0;
    for(int i=0;i<n;i++)
    {
        sum+=c[i];
    }
    for(int i=0;i<n;i++)
    {
        c[i]=c[i]/sum;
    }
}
void normalization_c(float128 *c,int n,float128 sum)
{
   // float128 sum=0.0;
    /*for(int i=0;i<n;i++)
    {
        sum+=c[i];
    }*/
    for(int i=0;i<n;i++)
    {
        c[i]=c[i]/sum;
    }
}
void cube_solver(float128 E2,float128 E1,float128 E0,float128 root[3],int &count_float128_root)
{
    float128 p,q;
    p=(3.0*E1-E2*E2)/3.0;                   //canonical form
    q=(2.0*E2*E2*E2-9.0*E2*E1+27.0*E0)/27.0;//cube equation
    float128 Q;
    Q=(p/3.0)*(p/3.0)*(p/3.0)+(q/2.0)*(q/2.0);
    complex128 temp_root1[3],temp_root2[3],a;
    if(abs(Q)<EPSILON)//Q=0
    {
        //matching roots
        if((abs(p)<EPSILON)&&(abs(q)<EPSILON))
        {
            //3 matching roots
            count_float128_root=3;
            root[0]=-E2/3.0;
            root[1]=-E2/3.0;
            root[2]=-E2/3.0;
        }
        else
        {
            //2 matching roots
            count_float128_root=3;
	    a=complex128(q/2.0Q,0.0Q);
            //a.real()=q/2.0;
            //a.imag()=0.0;
            cube_root(a,temp_root1);
            cube_root(a,temp_root2);
            root[0]=-2.0*cbrt((q/2.0))-E2/3.0;
            root[1]=-temp_root1[1].real()-temp_root1[2].real()-E2/3.0;
            root[2]=-temp_root1[1].real()-temp_root1[2].real()-E2/3.0;
        }
    }
    else
    {
        if(Q>0.0)
        {
            //1 real ,2 complex root
            count_float128_root=1;
            root[0]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
            root[1]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
            root[2]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
        }
        else
        {
            //3 real different root
            count_float128_root=3;
	    a=complex128(-q/2.0,sqrt((-Q)));
            //a.real()=-q/2.0;
            //a.imag()=sqrt((-Q));
            cube_root(a,temp_root1);
	    a=complex128(-q/2.0,-sqrt((-Q)));
            //a.real()=-q/2.0;
            //a.imag()=-sqrt((-Q));
            cube_root(a,temp_root2);
            root[0]=temp_root1[0].real()+temp_root2[0].real()-E2/3.0;
            root[1]=temp_root1[1].real()+temp_root2[1].real()-E2/3.0;
            root[2]=temp_root1[2].real()+temp_root2[2].real()-E2/3.0;
        }
    }
}
void calc_phi(float128 C,float128 * S_i,float128 *B_i,float128* Phi_i,float128 m1,float128 m2,float128 Z,float128 B,float128 A,int n)
{
    float128 M=logq((__float128)((Z+m2*B)/(Z+m1*B)));
    float128 *Ln_phi_i;
    Ln_phi_i=(float128 *)malloc(n*sizeof(float128));

    for(int i=0;i<n;i++)
    {
        Ln_phi_i[i]=logq((__float128)(C/abs(Z-B)))+M*A/((m1-m2)*B)*(2.0*S_i[i]/A-B_i[i]/B)+B_i[i]*(Z-C)/B;
        Phi_i[i]=expq((__float128)Ln_phi_i[i]);
        //cout<<Z<<"   "<< B<<'\n';
    }
	free(Ln_phi_i);
}
void coeff(float128 &E2,float128 &E1,float128 &E0,float128 m1,float128 m2,float128 A,float128 B,float128 C)
{
    E2=(m1+m2)*B-(B+C);
    E1=A+m1*m2*B*B-(m1+m2)*B*(B+C);
    E0=-A*B-m1*m2*B*B*(B+C);

}
void calc_c(float128 *c_i,float128 &c,int n)
{
    c=0.0;
    for(int i=0;i<n;i++)
    {
        c+=c_i[i];
    }

}
void calc_si(float128* si,float128** A_ij,float128 *c,int n)
{
    for(int i=0;i<n; i++)
    {
        si[i]=0.0;
    }

    for(int i=0;i<n ;i++)
    {

        for(int j=0;j<n;j++)
        {
            si[i]+=A_ij[i][j]*c[j];
        }
    }
}
float128 delta(int i,int j)
{
    return float128(i==j);
}
void practical_liq(float128* z,float128 *K,float128 V,float128** practical_x_practical_K,int n)
{
    float128 sum=0.0;
    float128 temp_i,temp_j;
    for(int i=0;i<n;i++)
    {
        sum+=((K[i]-1)*(K[i]-1)*z[i])/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
    }
    for(int i=0;i<n;i++)
    {
        temp_i=z[i]/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
         
        for(int j=0;j<n;j++)
        {temp_j=z[j]/((1+(K[j]-1)*V)*(1+(K[j]-1)*V));
            practical_x_practical_K[i][j]=temp_i*(-V*delta(i,j)-(K[i]-1)*temp_j/sum);
        }
    }
}
void practical_gas(float128* z,float128 *K,float128 V,float128** practical_y_practical_K,int n)
{
    float128 sum=0.0;
    float128 temp_i,temp_j;
    for(int i=0;i<n;i++)
    {
        sum+=((K[i]-1)*(K[i]-1)*z[i])/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
    }
    for(int i=0;i<n;i++)
    {
       temp_i=z[i]/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
        
        for(int j=0;j<n;j++)
        { temp_j=z[j]/((1+(K[j]-1)*V)*(1+(K[j]-1)*V));
            practical_y_practical_K[i][j]=temp_i*((1-V)*delta(i,j)-K[i]*(K[i]-1)*temp_j/sum);
        }
    }
}
void practical_Z(float128 Z,float128 M1,float128 M2,float128 E2,float128 E1,float128 A,float128 B,float128 C,float128 *B_i,float128 *S_i,float128 *practical_Z_practical_C,int n)
{
    for(int i=0;i<n ; i++)
    {
        float128 dE2_dc,dE1_dc,dE0_dc;
        dE2_dc=(M1+M2-1.0)*B_i[i]-1.0;
        dE1_dc=2.0*S_i[i]+2.0*M1*M2*B*B_i[i]-(M1+M2)*(B_i[i]*(B+C)+B*(B_i[i]+1.0));
        dE0_dc=-2.0*S_i[i]*B-A*B_i[i]-M1*M2*(2.0*B*B_i[i]*(B+C)+B*B*(B_i[i]+1.0));
        practical_Z_practical_C[i]=-(Z*Z*dE2_dc+Z*dE1_dc+dE0_dc)/(3.0*Z*Z+2.0*Z*E2+E1);
    }
}
void practical_ln_phi(float128 M1,float128 M2,float128 A,float128 B,float128 C,float128 Z,float128 *B_i,float128* S_i,float128 *practical_Z_practical_C,float128 **practical_c_practical_K,float128 **practical_ln_phi_practical_K,float128 **A_ij, int n)
{
    float128 W1,W2,W3,W4,W5,W6,M,W7;
    float128 **practical_ln_phi_practical_c;
    practical_ln_phi_practical_c        =   (float128**)malloc(n*sizeof(float128*));
    for(int i=0;i<n;i++)
    {
        practical_ln_phi_practical_c[i]         =   (float128*)malloc(n*sizeof(float128));}
    W1=A/((M1-M2)*B);
    M=logq((__float128)((Z+M2*B)/(Z+M1*B)));
    for(int i=0;i<n;i++)
    {
//cin>>W7;
        W5=2.0*S_i[i]/A-B_i[i]/B;

        for(int j=0;j<n;j++)
        {
            W2=1.0/C-(practical_Z_practical_C[j]-B_i[j])/(Z-B);

            W3=(practical_Z_practical_C[j]+M2*B_i[j])/(Z+M2*B)-(practical_Z_practical_C[j]+M1*B_i[j])/(Z+M1*B);
            W4=M*(2.0*S_i[j]*B-A*B_i[j])/((M1-M2)*B*B);
            W6=(B_i[i]*B_i[j])/(B*B);//cout<<W6;
            practical_ln_phi_practical_c[i][j]=W2+W1*W5*W3+W4*W5+M*W1*((2.0*A*A_ij[i][j]-4.0*S_i[i]*S_i[j])/(A*A)+W6)-(Z-C)*W6+(B_i[i]/B)*(practical_Z_practical_C[j]-1);
            //cout<< practical_ln_phi_practical_c[i][j]<<" ";
        }
        //cout<<'\n';
    } //cout<<'\n';cout<<'\n';cout<<'\n';

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            practical_ln_phi_practical_K[i][j]=0.0;

            for(int k=0;k<n;k++)
            {
                practical_ln_phi_practical_K[i][j]+=practical_ln_phi_practical_c[i][k]*practical_c_practical_K[k][j];
            }
            //cout<< practical_ln_phi_practical_K[i][j]<<" ";
           // cout<<practical_ln_phi_practical_K[i][j];
        } //cout<<'\n';cout<<'\n';cout<<'\n';
    }
    for(int i=0;i<n;i++)
    {


	    free(practical_ln_phi_practical_c[i]);


    }
free(practical_ln_phi_practical_c);

}
void A_linar_calc(float128 **A_linar_ij,float128** practical_ln_phi_practical_K_liq,float128 **practical_ln_phi_practical_K_gas,float128* K_i,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int  j=0;j<n;j++)
        {
            A_linar_ij[i][j]=delta(i,j)+K_i[j]*(-practical_ln_phi_practical_K_liq[i][j]+practical_ln_phi_practical_K_gas[i][j]);
        }
    }
}
void B_linar_calc(float128*B_linar_i, float128* phi_liq,float128* phi_gas,float128* K_i,int n)
{
    for(int i=0;i<n;i++)
    {
        B_linar_i[i]=logq((__float128)(phi_liq[i]/(K_i[i]*phi_gas[i])));
    }
}
void K_calk(float128* K_i,float128* tetta_i,int n)
{
    for(int i=0;i<n;i++)
    {
        K_i[i]*=expq((__float128)tetta_i[i]);
//cout<<exp(tetta_i[i])<<"   "<<tetta_i[i]<<'\n';
    }

}
void Gauss(float128** A,float128* B,float128* x,int n)
{

    for(int i=0;i<n;i++)
    {


            bool flag=false;
            int k=i;
            for(;(k<n)&&!flag;k++)
            {
                if(abs(A[k][i])>EPSILON){flag=true;}
            }k--;
            float128 swap_temp;
            float128 a,b;
            for(int l=i;l<n;l++)
            {
                swap_temp=A[k][l];
                A[k][l]=A[i][l];
                A[i][l]=swap_temp;
            }

            swap_temp=B[k];
            B[k]=B[i];
            B[i]=swap_temp;
            a=A[i][i];
            for(k=i+1;k<n;k++)
            {
                b=A[k][i];
                for(int l=i;l<n;l++)
                {
                    A[k][l]=a*A[k][l]-b*A[i][l];
                }
                B[k]=a*B[k]-b*B[i];
            }


    }

    for(int i=n-1;i>=0;i--)
    {

        float128 a,b;
        a=A[i][i];
        for(int j=i-1;j>=0;j--)
        {
            b=A[j][i];
            A[j][i]=A[j][i]-b*A[i][i]/a;
            B[j]=B[j]-b/a*B[i];

        }
    }
    for(int i=0;i<n;i++)
    {
        x[i]=B[i]/A[i][i];
    }


}
int main()
{
    ifstream data;
    data.open("c1-c4-c10.txt");
    ofstream data_final;
    ofstream data_final_k_i;
    float128 *p,*p_crit,*p_r,*s_liq_i;
    float128 *t,*t_crit,*t_r,*s_gas_i;
    float128 *omega,*Omega_A,*Omega_B,*A_i,*B_i;
    float128 **A_ij,**Betta,**A_linar;
    float128 A_liq=0.0,B_liq=0.0,A_gas=0.0,B_gas=0.0,C_liq,C_gas,temp,T,P,v=0.0,Omega_A0=0.457235529,Omega_B0=0.077796074;
    float128 *x,*y,*z,*c,*k_i,*B_linar,*x_solve;
    float128 m1,m2;
    float128 root_liq[3],root_gas[3];
    m1=1.0+sqrt(2.0);
    m2=1.0-sqrt(2.0);
    float128 E2_liq,E1_liq,E0_liq,z_liq,P_min,T_min,P_max,T_max,y_sum,x_sum;
    float128 E2_gas,E1_gas,E0_gas,z_gas;
    int n,count_root_gas,count_root_liq,count_step_p,count_step_t;
    int templ,count_step_interation;
    float128 *phi_liq_i,*phi_gas_i,*dz_dc_liq,*dz_dc_gas,**d_ln_phi_dK_liq,**d_ln_phi_dK_gas,**dx_dK,**dy_dK;



    data>>n>>T_min>>P_min>>T_max>>P_max>>count_step_t>>count_step_p;
	cout<<"1";
    dz_dc_liq   =   (float128*)malloc(n*sizeof(float128));
    dz_dc_gas   =   (float128*)malloc(n*sizeof(float128));
    B_linar     =   (float128*)malloc(n*sizeof(float128));
    x_solve     =   (float128*)malloc(n*sizeof(float128));
    s_liq_i     =   (float128*)malloc(n*sizeof(float128));
    s_gas_i     =   (float128*)malloc(n*sizeof(float128));
    phi_liq_i   =   (float128*)malloc(n*sizeof(float128));
    phi_gas_i   =   (float128*)malloc(n*sizeof(float128));
    p           =   (float128*)malloc(n*sizeof(float128));
    p_crit      =   (float128*)malloc(n*sizeof(float128));
    p_r         =   (float128*)malloc(n*sizeof(float128));
    t           =   (float128*)malloc(n*sizeof(float128));
    t_crit      =   (float128*)malloc(n*sizeof(float128));
    t_r         =   (float128*)malloc(n*sizeof(float128));
    omega       =   (float128*)malloc(n*sizeof(float128));
    Omega_A     =   (float128*)malloc(n*sizeof(float128));
    Omega_B     =   (float128*)malloc(n*sizeof(float128));
    A_i         =   (float128*)malloc(n*sizeof(float128));
    B_i         =   (float128*)malloc(n*sizeof(float128));
    x           =   (float128*)malloc(n*sizeof(float128));
    y           =   (float128*)malloc(n*sizeof(float128));
    z           =   (float128*)malloc(n*sizeof(float128));
    c           =   (float128*)malloc(n*sizeof(float128));
    k_i         =   (float128*)malloc(n*sizeof(float128));



    A_ij                =   (float128**)malloc(n*sizeof(float128*));
    Betta               =   (float128**)malloc(n*sizeof(float128*));
    A_linar             =   (float128**)malloc(n*sizeof(float128*));
    dy_dK               =   (float128**)malloc(n*sizeof(float128));
    dx_dK               =   (float128**)malloc(n*sizeof(float128));
    d_ln_phi_dK_liq     =   (float128**)malloc(n*sizeof(float128*));
    d_ln_phi_dK_gas     =   (float128**)malloc(n*sizeof(float128*));
float128 z_c1;
    int i=0,j=0;
data_final_k_i.open("C1-C4-C10__T=410.95.txt");
for (i = 0; i < n; i++)
{
    data >> t_crit[i] >> p_crit[i] >> omega[i] >> temp >> c[i] >> z[i];

    
 

    A_ij[i] = (float128*)malloc(n * sizeof(float128));
    Betta[i] = (float128*)malloc(n * sizeof(float128));
    A_linar[i] = (float128*)malloc(n * sizeof(float128));
    dy_dK[i] = (float128*)malloc(n * sizeof(float128));
    dx_dK[i] = (float128*)malloc(n * sizeof(float128));
    d_ln_phi_dK_liq[i] = (float128*)malloc(n * sizeof(float128*));
    d_ln_phi_dK_gas[i] = (float128*)malloc(n * sizeof(float128*));
}
for (i = 0; i < n; i++)
{
    for (j = 0; j < n; j++)
    {
        data >> Betta[i][j];
    }
}
for( z_c1=0.33;z_c1<=0.99;z_c1+=0.91)
{
for(float128 z_c2=0.33;z_c1+z_c2<=0.99;z_c2+=0.91)
{z[0]=z_c1;
            z[1]=z_c2;
            z[2]=1.-z_c1-z_c2;
            cout<<z[0]<<" "<<z[1]<<" "<<z[2]<<'\n';
    for(T=T_min;T<=T_max;T+=(T_max-T_min+1)/count_step_t)
    {cout<<"11n";
        for (i = 0; i < n; i++) 
        {
            
            k_i[i] = p_crit[i] / P_min * expq((__float128)(5.31 * (1 + omega[i]) * (1 - t_crit[i] / T)));cout<<k_i[i]<<'\n';
        }char tempinput;
        //cin>>tempinput;
		float128 dp=(P_max-P_min)/count_step_p;
        for(P=P_min;P<=P_max;P+=dp)
        {
			cout<<"\n"<<P<<'\n';
			
            count_step_interation=0;
           /* string name="z_P=";
            char temp_str[15];
            cout<<"test11\n";
            name+=to_string(double(P));
            cout<<"test12";
            quadmath_snprintf ( temp_str,15, ".10Q",46,P);
            name+=string(temp_str);
            cout<<"test12";
            quadmath_snprintf ( temp_str,15, "%+-#*.10Q",46,T);
            name+=" T=";name+=string(temp_str);
            quadmath_snprintf ( temp_str,15, "%+-#*.10Q",46,z_c1);
            name+=" C=";name+=string(temp_str);name+=".txt";*/
            //data_final.open(name);
           // cout<<"test\n";
            for (i = 0; i < n; i++) 
            {
                t[i] = T;
                p[i] = P;

            }//cout<<z[0];
			//cin>>z[0];
            
	float128 kriterij;
    	normalization_P(p,p_crit,p_r,n);
    	normalization_T(t,t_crit,t_r,n);
        
    	
	do
	{
       
    		normalization_c(z,n);
    		calc_V(k_i,z,v,n);
    		//+
           //if (v <= 0.0) { v = 0.0; }
           //if (v >= 1.0) { v = 1.0; }
            //v = 0.5;
    		//cout<<"v"<<v<<'\n';
            float128 y0 = (1 - k_i[1]) / (k_i[0] - k_i[1]);
            
    		calc_x(k_i,z,x,v,n);//+
    		calc_y(k_i,z,y,v,n);//+
            //normalization_c(x, 2);
            //normalization_c(y, 2);
          
            //z[0] = (x[0] + y[0]) / 2;
            //z[1] = (x[1] + y[1]) / 2;
    		Omega_a_calc(Omega_A0,omega,Omega_A,t_r,n);
    		Omega_b_calc(Omega_B0,Omega_B,n);
             //cout<<"test\n";
    		calc_B_i(Omega_B,p_r,t_r,B_i,n);
    		calc_A_i(Omega_A,p_r,t_r,A_i,n);
    		calc_A_i_j(Betta,A_i,A_ij,n);
    		calc_A(A_ij,x,A_liq,n);
    		calc_A(A_ij,y,A_gas,n);
    		calc_B(B_i,x,B_liq,n);
    		calc_B(B_i,y,B_gas,n);
    		calc_c(x,C_liq,n);
    		calc_c(y,C_gas,n);
    		coeff(E2_liq,E1_liq,E0_liq,m1,m2,A_liq,B_liq,C_liq);
    		coeff(E2_gas,E1_gas,E0_gas,m1,m2,A_gas,B_gas,C_gas);


    		cube_solver(E2_liq,E1_liq,E0_liq,root_liq,count_root_liq);
  		cube_solver(E2_gas,E1_gas,E0_gas,root_gas,count_root_gas);
    		z_liq=max(max(root_liq[0],root_liq[1]),root_liq[2]);
    		if((z_liq>root_liq[0])&&(root_liq[0]>B_liq)){z_liq=root_liq[0];}
    		if((z_liq>root_liq[1])&&(root_liq[1]>B_liq)){z_liq=root_liq[1];}
    		if((z_liq>root_liq[2])&&(root_liq[2]>B_liq)){z_liq=root_liq[2];}
    		//z min positive rootd?
    		z_gas=max(max(root_gas[0],root_gas[1]),root_gas[2]);
    		for(int k=0;k<n;k++)
    		{
			//z_liq-=z[k]*c[k];
			//z_gas-=z[k]*c[k];
    		}
    		//
    		calc_si(s_liq_i,A_ij,x,n);
    		calc_si(s_gas_i,A_ij,y,n);


    		calc_phi(C_liq,s_liq_i,B_i,phi_liq_i,m1,m2,z_liq,B_liq,A_liq,n);
    		calc_phi(C_gas,s_gas_i,B_i,phi_gas_i,m1,m2,z_gas,B_gas,A_gas,n);

     kriterij=0.0;
    
    for(int k=0;k<n;k++)
    {
        kriterij+=cabsq((__float128)((x[k]* phi_liq_i[k] )/(y[k]* phi_gas_i[k]) -1.0));
       // cout<<phi_gas_i[k];
       // cout << k_i[k] << "   " << phi_liq_i[k] / phi_gas_i[k] << '\n';
       // float128 qw;
       // cin >> qw;
    }//cout<<"nya";
   // cout<<"Kriterij: "<<kriterij<<'\n';
    if(10*kriterij>EPSILON)
    {
        practical_Z(z_liq,m1,m2,E2_liq,E1_liq,A_liq,B_liq,C_liq,B_i,s_liq_i,dz_dc_liq,n);
        practical_Z(z_gas,m1,m2,E2_gas,E1_gas,A_gas,B_gas,C_gas,B_i,s_gas_i,dz_dc_gas,n);
        practical_liq(z,k_i,v,dx_dK,n);
        practical_gas(z,k_i,v,dy_dK,n);//cout<<"nya";

   // cin>>templ;
        practical_ln_phi(m1,m2,A_liq,B_liq,C_liq,z_liq,B_i,s_liq_i,dz_dc_liq,dx_dK,d_ln_phi_dK_liq,A_ij,n);
          //  cin>>templ;
        practical_ln_phi(m1,m2,A_gas,B_gas,C_gas,z_gas,B_i,s_gas_i,dz_dc_gas,dy_dK,d_ln_phi_dK_gas,A_ij,n);
        A_linar_calc(A_linar,d_ln_phi_dK_liq,d_ln_phi_dK_gas,k_i,n);
        B_linar_calc(B_linar,phi_liq_i,phi_gas_i,k_i,n);
        Gauss(A_linar,B_linar,x_solve,n);
//cout<<"K_calk";
        K_calk(k_i,x_solve,n);
        
    }
  //  cout<<"k_i:\n";

    for(int w=0;w<n;w++)
    {//cout<<"K_i: "<<k_i[w]<<"  X_i:"<<x[w]<<"  Y_i:"<<y[w]<<'\n';
//data_final <<count_step_interation<<"    "<<k_i[0]<<"    "<<x[0]<<"    "<<y[0]<<"    "<<k_i[1]<<"    "<<x[1]<<"    "<<y[1]<<"    "<<v<<'\n';
count_step_interation++;
    }//cin>>templ;
   // cout<<"Kriterij: "<<kriterij<<'\n';
	}while((kriterij>EPSILON)&&(count_step_interation<1000));
    complex128 a,root[3];
    float128 root_d[3];
	a=complex128(0.0,-1);
    //a.real()=0.0;
    //a.imag()=-1.0;
    cube_root(a,root);
     cube_solver(-8.0,21.0,-18.0,root_d,i);
//	cout<<"z_liq  "<<z_liq<<"  z_gas  "<<z_gas<<'\n';
//	cout<<"A_liq  "<<A_liq<<"  A_gas  "<<A_gas<<'\n';
//	cout<<"B_liq  "<<B_liq<<"  B_gas  "<<B_gas<<'\n';
    //cout<<root_d[0]<<"   "<<'\n';
    //cout<<root_d[1]<<"   "<<'\n';
    //cout<<root_d[2]<<"   "<<'\n';
    //cout<<root[0].real()<<"   "<<root[0].imag()<<'\n';
    //cout<<root[1].real()<<"   "<<root[1].imag()<<'\n';
    //cout<<root[2].real()<<"   "<<root[2].imag()<<'\n';
    //cout<<p_r[0];

  /*  data_final<< P<<" "<<T<<" ";

    for(int w=0;w<n;w++)
    {data_final<<k_i[w]<<" ";
    }
    data_final<<'\n';*/

data_final.close();
if(k_i[0]!=k_i[0]){k_i[0]=-1;k_i[1]=-1;}
calc_c(x,x_sum,2);
calc_c(y,y_sum,2);
//z[0]=(x[0]+y[0])/2;
//z[1]=(x[1]+y[1])/2;
//cout<<z[0];
//			cin>>z[0];
//normalization_c(x,2);
//normalization_c(y,2);
float128 k_w[2];
for (i = 0; i < n; i++) 
{
            k_w[i] = p_crit[i] / P * expq((__float128)(5.31 * (1 + omega[i]) * (1 - t_crit[i] / T)));
			cout<<k_w[i]<<'\n';
}
float128 v_w=0.0;
calc_V(k_w, z, v_w, n);
float128 x_w[2];
float128 y_w[2];
calc_x(k_w, z, x_w, v_w, n);
calc_y(k_w, z, y_w, v_w, n);
data_final_k_i<<z[0]<<"    "<<z[1]<<"    "<<z[2]<<"    "<<P<<"    "<<T<<"    "<<k_i[0]<<"    "<<x[0]<<"    "<<y[0]<<"    "<<k_i[1]<<"    "<<x[1]<<"    "<<y[1]<<"    "<<k_i[2]<<"    "<<x[2]<<"    "<<y[2]<<"    "<<v<<"    "<<x_w[0]<<"    "<<y_w[0]<<"    "<<x_w[1]<<"    "<<y_w[1]<<'\n';
}
data_final_k_i<<'\n';
}cout<<"12345";}}
data_final_k_i.close();
    free(dz_dc_liq   );
    free(dz_dc_gas   );
    free(B_linar     );
    free(x_solve     );
    free(s_liq_i     );
    free(s_gas_i     );
    free(phi_liq_i   );
    free(phi_gas_i   );
    free(p           );
    free(p_crit      );
    free(p_r         );
    free(t           );
    free(t_crit      );
    free(t_r         );
    free(omega       );
    free(Omega_A     );
    free(Omega_B     );
    free(A_i         );
    free(B_i         );
    free(x           );
    free(y           );
    free(z           );
    free(c           );
    free(k_i         );


    for(i=0;i<n;i++)
    {

        free(A_ij[i]);
        free(Betta[i]);
        free(A_linar[i]);
        free(dy_dK[i]  );
        free(dx_dK[i]  );
        free(d_ln_phi_dK_liq[i] );
        free(d_ln_phi_dK_gas[i] );
    }
	free(A_ij);
        free(Betta);
        free(A_linar);
        free(dy_dK  );
        free(dx_dK  );
        free(d_ln_phi_dK_liq );
        free(d_ln_phi_dK_gas );



    return 0;
}
