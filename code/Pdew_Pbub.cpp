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
#define  EPSILON 1e-18Q
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
        if(omega[i]>0.49)
        {
            temp_calc=1+(0.379642+1.48503*omega[i]-0.164423*omega[i]*omega[i]+0.016666*omega[i]*omega[i]*omega[i])*(1-sqrtq((__float128)T_r[i]));
        }
        else
        {
            temp_calc=1+(0.37464+1.54226*omega[i]-0.26992*omega[i]*omega[i])*(1-sqrtq((__float128)T_r[i]));
        }
        
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
        cout<<k_i[i]<<" ";
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
    cout<<"Ln_phi_i"<<M<<" "<<Z<<" "<<B<<"\n";
    for(int i=0;i<n;i++)
    {
        Ln_phi_i[i]=logq((__float128)(C/abs(Z-B)))+M*A/((m1-m2)*B)*(2.0*S_i[i]/A-B_i[i]/B)+B_i[i]*(Z-C)/B;
       // cout<<Ln_phi_i[i]<<"\n";
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
void practical_liq(float128* z,float128 *K,float128** practical_x_practical_K,int n)
{

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            practical_x_practical_K[i][j]=-z[i]*delta(i,j)/K[i]/K[i];
        }
    }
}

void practical_gas(float128* z,float128 *K,float128** practical_y_practical_K,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            practical_y_practical_K[i][j]=z[i]*delta(i,j);
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
void practical_Z_dp(float128 Z,float128 M1,float128 M2,float128 E2,float128 E1,float128 A,float128 B,float128 C,float128 &practical_Z_dp,float128 p,int n)
{
    
    float128 dE2_dp,dE1_dp,dE0_dp;
    dE2_dp=(M1+M2-1.0)*B/p;
    dE1_dp=(A+2.0*M1*M2*B*B-(M1+M2)*(2*B*B+B*C))/p;
    dE0_dp=(-2*A*B-M1*M2*(3*B*B*B+2*B*B*C))/p;
    practical_Z_dp=-(Z*Z*dE2_dp+Z*dE1_dp+dE0_dp)/(3.0*Z*Z+2.0*Z*E2+E1);
}
void practical_ln_phi(float128 M1,float128 M2,float128 A,float128 B,float128 C,float128 Z,float128 *B_i,float128* S_i,float128 *practical_Z_practical_C,float128 **practical_c_practical_K,float128 **&practical_ln_phi_practical_K,float128 **A_ij, int n)
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
            practical_ln_phi_practical_K[i][j]=0.0Q;
             

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
void practical_ln_phi_dp(float128 z,float128 dz_dp,float128 M1,float128 M2,float128 A,float128 *S_i,float128 B,float128 *B_i,float128 * practical_ln_phi_practical_p,float128 p,int n)
{
    float128 W1,W5,W7,W8;
    W1=A/((M1-M2)*B);
    W7=(dz_dp-B/p)/(z-B);
    W8=(dz_dp+M2*B/p)/(z+M2*B)-(dz_dp+M1*B/p)/(z+M1*B);
    for(int i=0;i<n;i++)
    {
        W5=(2.0Q*S_i[i]/A)-B_i[i]/B;
        practical_ln_phi_practical_p[i]=W7+W1*W5*W8+B_i[i]*dz_dp/B;

    }
}
void A_linar_calc(float128 **A_linar_ij,  float128 ** practical_ln_phi_practical_K_liq,float128 **practical_ln_phi_practical_K_gas,float128* K_i,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int  j=0;j<n;j++)
        {
            A_linar_ij[i][j]=delta(i,j)+K_i[j]*(-practical_ln_phi_practical_K_liq[i][j]+practical_ln_phi_practical_K_gas[i][j]);
        }
    }
}
void A_linar_calc_dew(float128 **A_linar,float128 *d_ln_phi_dp_gas,float128 *d_ln_phi_dp_liq,float128 **d_ln_phi_dK,float128 *K_i,float128 *z_i,float128 P,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            A_linar[i][j]=delta(i,j)-K_i[j]*d_ln_phi_dK[i][j];
        }
    }
    for(int j=0;j<n;j++)
    {
        A_linar[n][j]=-z_i[j]/K_i[j];
        A_linar[j][n]=P*(d_ln_phi_dp_gas[j]-d_ln_phi_dp_liq[j]);
    }
    A_linar[n][n]=0;
}
void A_linar_calc_bub(float128 **A_linar,float128 *d_ln_phi_dp_gas,float128 *d_ln_phi_dp_liq,float128 **d_ln_phi_dK,float128 *K_i,float128 *z_i,float128 P,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            A_linar[i][j]=delta(i,j)+K_i[j]*d_ln_phi_dK[i][j];
        }
    }
    for(int j=0;j<n;j++)
    {
        A_linar[n][j]=z_i[j]*K_i[j];
        A_linar[j][n]=P*(d_ln_phi_dp_gas[j]-d_ln_phi_dp_liq[j]);
    }
    A_linar[n][n]=0;
}
void B_linar_calc_dew(float128*B_linar_i, float128* phi_liq,float128* phi_gas,float128* K_i,float128* z_i,int n)
{
    float128 sum=0.0Q;
    for(int i=0;i<n;i++)
    {
        B_linar_i[i]=logq((__float128)(phi_liq[i]/(K_i[i]*phi_gas[i])));
        sum+=z_i[i]/K_i[i];
    }
    B_linar_i[n]=1.0Q-sum;

}
void B_linar_calc_bub(float128*B_linar_i, float128* phi_liq,float128* phi_gas,float128* K_i,float128* z_i,int n)
{
    float128 sum=0.0Q;
    for(int i=0;i<n;i++)
    {
        B_linar_i[i]=logq((__float128)(phi_liq[i]/(K_i[i]*phi_gas[i])));
        sum+=z_i[i]*K_i[i];
    }
    B_linar_i[n]=1.0Q-sum;

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
    data.open("for_2_comp.txt");
    ofstream data_final;
    ofstream data_final_k_i;
    float128 *p,*p_crit,*p_r_dew,*p_r_bub,*s_liq_i_dew,*s_liq_i_bub;
    float128 *t,*t_crit,*t_r,*s_gas_i_dew,*s_gas_i_bub;
    float128 *omega,*Omega_A,*Omega_B,*A_i_dew,*A_i_bub,*B_i_dew,*B_i_bub;
    float128 **A_ij_dew,**A_ij_bub,**Betta,**A_linar_dew,**A_linar_bub;
    float128 A_liq_dew=0.0,A_liq_bub=0.0,B_liq_dew=0.0,B_liq_bub=0.0,
            A_gas_dew=0.0, A_gas_bub=0.0,B_gas_dew=0.0,B_gas_bub=0.0,
            C_liq_dew,C_liq_bub,C_gas_dew,C_gas_bub,temp,T,P,v=0.0,Omega_A0=0.457235529Q,Omega_B0=0.077796074Q;
    float128 *x_dew,*x_bub,*y_dew,*y_bub,*z,*c,*k_i_dew,*k_i_bub,*B_linar_dew,*B_linar_bub,*x_solve_dew,*x_solve_bub;
    float128 m1,m2;
    float128 root_liq_dew[3],root_gas_dew[3],root_liq_bub[3],root_gas_bub[3];
    m1=1.0+sqrt(2.0);
    m2=1.0-sqrt(2.0);
    float128 P_min,T_min,P_max,T_max,y_sum,x_sum;
    float128 E2_liq_dew,E1_liq_dew,E0_liq_dew,z_liq_dew;
    float128 E2_gas_dew,E1_gas_dew,E0_gas_dew,z_gas_dew;
    float128 E2_liq_bub,E1_liq_bub,E0_liq_bub,z_liq_bub;
    float128 E2_gas_bub,E1_gas_bub,E0_gas_bub,z_gas_bub;
    float128 *practical_ln_phi_practical_p_liq_dew,*practical_ln_phi_practical_p_liq_bub,*practical_ln_phi_practical_p_gas_dew,*practical_ln_phi_practical_p_gas_bub;
    int n,count_root_gas_dew,count_root_liq_dew,count_root_gas_bub,count_root_liq_bub,count_step_p,count_step_t;
    int templ,count_step_interation;
    float128 *phi_liq_i_dew,*phi_gas_i_dew,*phi_liq_i_bub,*phi_gas_i_bub,*dz_dc_liq_dew,*dz_dc_gas_bub,**d_ln_phi_dK_liq,**d_ln_phi_dK_gas,**dx_dK,**dy_dK;



    data>>n>>T_min>>P_min>>T_max>>P_max>>count_step_t>>count_step_p;
	cout<<"1";
	
	
    practical_ln_phi_practical_p_liq_dew   =   (float128*)malloc(n*sizeof(float128));
    practical_ln_phi_practical_p_liq_bub   =   (float128*)malloc(n*sizeof(float128));
    practical_ln_phi_practical_p_gas_dew   =   (float128*)malloc(n*sizeof(float128));
    practical_ln_phi_practical_p_gas_bub   =   (float128*)malloc(n*sizeof(float128));
    dz_dc_liq_dew   =   (float128*)malloc(n*sizeof(float128));
    dz_dc_gas_bub   =   (float128*)malloc(n*sizeof(float128));
    B_linar_dew     =   (float128*)malloc((n+1)*sizeof(float128));
    B_linar_bub     =   (float128*)malloc((n+1)*sizeof(float128));
    x_solve_dew     =   (float128*)malloc((n+1)*sizeof(float128));
    x_solve_bub     =   (float128*)malloc((n+1)*sizeof(float128));
    s_liq_i_dew     =   (float128*)malloc(n*sizeof(float128));
    s_gas_i_dew     =   (float128*)malloc(n*sizeof(float128));
    s_liq_i_bub     =   (float128*)malloc(n*sizeof(float128));
    s_gas_i_bub     =   (float128*)malloc(n*sizeof(float128));
    phi_liq_i_dew   =   (float128*)malloc(n*sizeof(float128));
    phi_gas_i_dew   =   (float128*)malloc(n*sizeof(float128));
    phi_liq_i_bub   =   (float128*)malloc(n*sizeof(float128));
    phi_gas_i_bub   =   (float128*)malloc(n*sizeof(float128));
    p               =   (float128*)malloc(n*sizeof(float128));
    p_crit          =   (float128*)malloc(n*sizeof(float128));
    p_r_dew         =   (float128*)malloc(n*sizeof(float128));
    p_r_bub         =   (float128*)malloc(n*sizeof(float128));
    t               =   (float128*)malloc(n*sizeof(float128));
    t_crit          =   (float128*)malloc(n*sizeof(float128));
    t_r             =   (float128*)malloc(n*sizeof(float128));
    omega           =   (float128*)malloc(n*sizeof(float128));
    Omega_A         =   (float128*)malloc(n*sizeof(float128));
    Omega_B         =   (float128*)malloc(n*sizeof(float128));
    A_i_dew         =   (float128*)malloc(n*sizeof(float128));
    B_i_dew         =   (float128*)malloc(n*sizeof(float128));
    A_i_bub         =   (float128*)malloc(n*sizeof(float128));
    B_i_bub         =   (float128*)malloc(n*sizeof(float128));
    x_dew           =   (float128*)malloc(n*sizeof(float128));
    x_bub           =   (float128*)malloc(n*sizeof(float128));
    y_dew           =   (float128*)malloc(n*sizeof(float128));
    y_bub           =   (float128*)malloc(n*sizeof(float128));
    z               =   (float128*)malloc(n*sizeof(float128));
    c               =   (float128*)malloc(n*sizeof(float128));
    k_i_dew         =   (float128*)malloc(n*sizeof(float128));
    k_i_bub         =   (float128*)malloc(n*sizeof(float128));



    A_ij_dew            =   (float128**)malloc(n*sizeof(float128*));
    A_ij_bub            =   (float128**)malloc(n*sizeof(float128*));
    Betta               =   (float128**)malloc(n*sizeof(float128*));
    A_linar_dew         =   (float128**)malloc((n+1)*sizeof(float128*));
    A_linar_bub         =   (float128**)malloc((n+1)*sizeof(float128*));
    dy_dK               =   (float128**)malloc(n*sizeof(float128));
    dx_dK               =   (float128**)malloc(n*sizeof(float128));
    d_ln_phi_dK_liq     =   (float128**)malloc(n*sizeof(float128*));
    d_ln_phi_dK_gas     =   (float128**)malloc(n*sizeof(float128*));

    int i=0,j=0;
    //data_final_k_i.open("ki_for_10_comp_T=410.95.txt");
    data_final.open("P_for_10_comp_T=300-400.txt");
    for (i = 0; i < n; i++)
    {
        data >> t_crit[i] >> p_crit[i] >> omega[i] >> temp >> c[i] >> z[i];

        A_ij_dew[i] = (float128*)malloc(n * sizeof(float128));
        A_ij_bub[i] = (float128*)malloc(n * sizeof(float128));
        Betta[i] = (float128*)malloc(n * sizeof(float128));
        A_linar_dew[i] = (float128*)malloc((n+1) * sizeof(float128));
        A_linar_bub[i] = (float128*)malloc((n+1) * sizeof(float128));
        dy_dK[i] = (float128*)malloc(n * sizeof(float128));
        dx_dK[i] = (float128*)malloc(n * sizeof(float128));
        d_ln_phi_dK_liq[i] = (float128*)malloc(n * sizeof(float128));
        d_ln_phi_dK_gas[i] = (float128*)malloc(n * sizeof(float128));
    }
    A_linar_dew[n] = (float128*)malloc((n+1) * sizeof(float128));
    A_linar_bub[n] = (float128*)malloc((n+1) * sizeof(float128));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            data >> Betta[i][j];
        }
    }
    normalization_c(z,n);

    float128 p_dew,p_bub;
    p_dew=p_bub=P_min;
    for(T=T_min;T<=T_max;T+=(T_max-T_min)/count_step_t)
    {
        count_step_interation=0;
        p_dew=p_bub=(p_dew+p_bub)/2.0Q;
        p_dew=p_bub=P_min;
        for (i = 0; i < n; i++)
        {
            
            k_i_dew[i] = p_crit[i] / (p_dew) * expq((__float128)(5.31 * (1 + omega[i]) * (1 - t_crit[i] / T)));//cout<<k_i_dew[i]<<'\n';
            k_i_bub[i] = p_crit[i] / (p_bub) * expq((__float128)(5.31 * (1 + omega[i]) * (1 - t_crit[i] / T)));//cout<<k_i_bub[i]<<'\n';

        }
        for (i = 0; i < n; i++)
        {
            t[i] = T;
        }
        float128 kriterij_dew,kriterij_bub;
    	normalization_T(t,t_crit,t_r,n);//+
    	
	do
	{
            for (i = 0; i < n; i++)
            {
                p[i] = p_dew;
            }
            normalization_P(p,p_crit,p_r_dew,n);//+
           
            
            for (i = 0; i < n; i++)
            {
                p[i] = p_bub;
            }
            normalization_P(p,p_crit,p_r_bub,n);//+
            char mmmmmm;
            //xi=zi/(1+(ki-1)v)
            //v=0.0,xi=zi,pbub
            //v=1.0,xi=zi/ki,pdew
            calc_x(k_i_dew,z,x_dew,1.0Q,n);//+
            calc_y(k_i_dew,z,y_dew,1.0Q,n);//+
            calc_x(k_i_bub,z,x_bub,0.0Q,n);//+
            calc_y(k_i_bub,z,y_bub,0.0Q,n);//+
            for(int i=0;i<n;i++)
            {
                cout<<x_dew[i]<<"    "<<y_dew[i]<<"   "<<k_i_dew[i]<<"   "<<x_bub[i]<<"   "<<y_bub[i]<<"   "<<k_i_bub[i]<<"   "<<z[i]<<'\n';
            }
             
            cout<<'\n';
           /* char temp_nya;
            cin>>temp_nya;*/
            Omega_a_calc(Omega_A0,omega,Omega_A,t_r,n);//+
            Omega_b_calc(Omega_B0,Omega_B,n);//*
          /*  normalization_c(x_dew,n);
            normalization_c(x_bub,n);
            normalization_c(y_dew,n);
            normalization_c(y_bub,n);*/
            
            calc_B_i(Omega_B,p_r_dew,t_r,B_i_dew,n);//+
            calc_A_i(Omega_A,p_r_dew,t_r,A_i_dew,n);//+
            
            calc_B_i(Omega_B,p_r_bub,t_r,B_i_bub,n);//+
            calc_A_i(Omega_A,p_r_bub,t_r,A_i_bub,n);//+

            calc_A_i_j(Betta,A_i_dew,A_ij_dew,n);//+
            calc_A_i_j(Betta,A_i_bub,A_ij_bub,n);//+
            

            calc_A(A_ij_dew,x_dew,A_liq_dew,n);
            calc_A(A_ij_dew,y_dew,A_gas_dew,n);
            calc_A(A_ij_bub,x_bub,A_liq_bub,n);
            calc_A(A_ij_bub,y_bub,A_gas_bub,n);

            calc_B(B_i_dew,x_dew,B_liq_dew,n);//+
            calc_B(B_i_dew,y_dew,B_gas_dew,n);//+
            calc_B(B_i_bub,x_bub,B_liq_bub,n);//+
            calc_B(B_i_bub,y_bub,B_gas_bub,n);//+
            
            calc_c(x_dew,C_liq_dew,n);//+
            calc_c(y_dew,C_gas_dew,n);//+
            calc_c(x_bub,C_liq_bub,n);//+
            calc_c(y_bub,C_gas_bub,n);//+
            cout<<"C";
            cout<<C_liq_dew<<" "<<C_gas_dew<<"\n";
            cout<<C_liq_bub<<" "<<C_gas_bub<<"\n";
            cin>>mmmmmm;


            coeff(E2_liq_dew,E1_liq_dew,E0_liq_dew,m1,m2,A_liq_dew,B_liq_dew,C_liq_dew);
            coeff(E2_gas_dew,E1_gas_dew,E0_gas_dew,m1,m2,A_gas_dew,B_gas_dew,C_gas_dew);
            coeff(E2_liq_bub,E1_liq_bub,E0_liq_bub,m1,m2,A_liq_bub,B_liq_bub,C_liq_bub);
            coeff(E2_gas_bub,E1_gas_bub,E0_gas_bub,m1,m2,A_gas_bub,B_gas_bub,C_gas_bub);


            cube_solver(E2_liq_dew,E1_liq_dew,E0_liq_dew,root_liq_dew,count_root_liq_dew);
            cube_solver(E2_gas_dew,E1_gas_dew,E0_gas_dew,root_gas_dew,count_root_gas_dew);
            cube_solver(E2_liq_bub,E1_liq_bub,E0_liq_bub,root_liq_bub,count_root_liq_bub);
            cube_solver(E2_gas_bub,E1_gas_bub,E0_gas_bub,root_gas_bub,count_root_gas_bub);
            cout<<"\n\n\n";
            
            z_liq_dew=max(max(root_liq_dew[0],root_liq_dew[1]),root_liq_dew[2]);
            if((z_liq_dew>root_liq_dew[0])&&(root_liq_dew[0]>B_liq_dew)){z_liq_dew=root_liq_dew[0];}
            if((z_liq_dew>root_liq_dew[1])&&(root_liq_dew[1]>B_liq_dew)){z_liq_dew=root_liq_dew[1];}
            if((z_liq_dew>root_liq_dew[2])&&(root_liq_dew[2]>B_liq_dew)){z_liq_dew=root_liq_dew[2];}
            z_gas_dew=max(max(root_gas_dew[0],root_gas_dew[1]),root_gas_dew[2]);

            z_liq_bub=max(max(root_liq_bub[0],root_liq_bub[1]),root_liq_bub[2]);
            if((z_liq_bub>root_liq_bub[0])&&(root_liq_bub[0]>B_liq_bub)){z_liq_bub=root_liq_bub[0];}
            if((z_liq_bub>root_liq_bub[1])&&(root_liq_bub[1]>B_liq_bub)){z_liq_bub=root_liq_bub[1];}
            if((z_liq_bub>root_liq_bub[2])&&(root_liq_bub[2]>B_liq_bub)){z_liq_bub=root_liq_bub[2];}
            z_gas_bub=max(max(root_gas_bub[0],root_gas_bub[1]),root_gas_bub[2]);

            
            cout<<root_liq_dew[0]<<"  "<<root_liq_dew[1]<<"  "<<root_liq_dew[2]<<"     "<<z_liq_dew<<"\n";
            cout<<root_gas_dew[0]<<"  "<<root_gas_dew[1]<<"  "<<root_gas_dew[2]<<"     "<<z_gas_dew<<"\n";
            cout<<root_liq_bub[0]<<"  "<<root_liq_bub[1]<<"  "<<root_liq_bub[2]<<"     "<<z_liq_bub<<"\n";
            cout<<root_gas_bub[0]<<"  "<<root_gas_bub[1]<<"  "<<root_gas_bub[2]<<"     "<<z_gas_bub<<"\n";
            calc_si(s_liq_i_dew,A_ij_dew,x_dew,n);
            calc_si(s_gas_i_dew,A_ij_dew,y_dew,n);
            calc_si(s_liq_i_bub,A_ij_bub,x_bub,n);
            calc_si(s_gas_i_bub,A_ij_bub,y_bub,n);


            calc_phi(C_liq_dew,s_liq_i_dew,B_i_dew,phi_liq_i_dew,m1,m2,z_liq_dew,B_liq_dew,A_liq_dew,n);
            calc_phi(C_gas_dew,s_gas_i_dew,B_i_dew,phi_gas_i_dew,m1,m2,z_gas_dew,B_gas_dew,A_gas_dew,n);
            calc_phi(C_liq_bub,s_liq_i_bub,B_i_bub,phi_liq_i_bub,m1,m2,z_liq_bub,B_liq_bub,A_liq_bub,n);
            calc_phi(C_gas_bub,s_gas_i_bub,B_i_bub,phi_gas_i_bub,m1,m2,z_gas_bub,B_gas_bub,A_gas_bub,n);

            kriterij_dew=0.0;
            kriterij_bub=0.0;
            
            for(int k=0;k<n;k++)
            {
                kriterij_dew+=cabsq((__float128)((x_dew[k]* phi_liq_i_dew[k] )/(y_dew[k]* phi_gas_i_dew[k]) -1.0));
                cout<<phi_liq_i_bub[k]<<"  "<<phi_gas_i_bub[k]<<"   "<<(x_dew[k]* phi_liq_i_dew[k] )/(y_dew[k]* phi_gas_i_dew[k])<<'\n';
                kriterij_bub+=cabsq((__float128)((x_bub[k]* phi_liq_i_bub[k] )/(y_bub[k]* phi_gas_i_bub[k]) -1.0));
            }
            cout<<'\n';
            if(kriterij_dew+kriterij_bub>EPSILON)
            {
                float128 dZ_dp_liq_dew,dZ_dp_gas_dew,dZ_dp_liq_bub,dZ_dp_gas_bub;

                practical_Z_dp(z_liq_dew,m1,m2,E2_liq_dew,E1_liq_dew,A_liq_dew,B_liq_dew,C_liq_dew,dZ_dp_liq_dew,p_dew,n);
                cout<<"\nz_liq_dew "<<z_liq_dew;
                cout<<"\ndZ_dp_liq_dew"<<dZ_dp_liq_dew<<'\n';
                cout<<"\nE2_liq_dew "<<E2_liq_dew;
                cout<<"\nE1_liq_dew "<<E1_liq_dew<<"\n";
                cin>>mmmmmm;
                practical_Z_dp(z_gas_dew,m1,m2,E2_gas_dew,E1_gas_dew,A_gas_dew,B_gas_dew,C_gas_dew,dZ_dp_gas_dew,p_dew,n);
                practical_Z_dp(z_liq_bub,m1,m2,E2_liq_bub,E1_liq_bub,A_liq_bub,B_liq_bub,C_liq_bub,dZ_dp_liq_bub,p_bub,n);
                practical_Z_dp(z_gas_bub,m1,m2,E2_gas_bub,E1_gas_bub,A_gas_bub,B_gas_bub,C_gas_bub,dZ_dp_gas_bub,p_bub,n);


                practical_Z(z_liq_dew,m1,m2,E2_liq_dew,E1_liq_dew,A_liq_dew,B_liq_dew,C_liq_dew,B_i_dew,s_liq_i_dew,dz_dc_liq_dew,n);
                practical_Z(z_gas_bub,m1,m2,E2_gas_bub,E1_gas_bub,A_gas_bub,B_gas_bub,C_gas_bub,B_i_bub,s_gas_i_bub,dz_dc_gas_bub,n);
                practical_liq(z,k_i_dew,dx_dK,n);
                practical_gas(z,k_i_bub,dy_dK,n);//cout<<"nya";
   // cin>>templ;
                practical_ln_phi(m1,m2,A_liq_dew,B_liq_dew,C_liq_dew,z_liq_dew,B_i_dew,s_liq_i_dew,dz_dc_liq_dew,dx_dK,d_ln_phi_dK_liq,A_ij_dew,n);
                practical_ln_phi(m1,m2,A_gas_bub,B_gas_bub,C_gas_bub,z_gas_bub,B_i_bub,s_gas_i_bub,dz_dc_gas_bub,dy_dK,d_ln_phi_dK_gas,A_ij_bub,n);
                practical_ln_phi_dp(z_liq_dew,dZ_dp_liq_dew,m1,m2,A_liq_dew,s_liq_i_dew, B_liq_dew,
                                B_i_dew,practical_ln_phi_practical_p_liq_dew,p_dew,n);
                practical_ln_phi_dp(z_gas_dew,dZ_dp_gas_dew,m1,m2,A_gas_dew,s_gas_i_dew, B_gas_dew,
                                B_i_dew,practical_ln_phi_practical_p_gas_dew,p_dew,n);
                practical_ln_phi_dp(z_liq_bub,dZ_dp_liq_bub,m1,m2,A_liq_bub,s_liq_i_bub, B_liq_bub,
                                B_i_bub,practical_ln_phi_practical_p_liq_bub,p_bub,n);
                practical_ln_phi_dp(z_gas_bub,dZ_dp_gas_bub,m1,m2,A_gas_bub,s_gas_i_bub, B_gas_bub,
                                B_i_bub,practical_ln_phi_practical_p_gas_bub,p_bub,n);
                A_linar_calc_dew(A_linar_dew,practical_ln_phi_practical_p_gas_dew,practical_ln_phi_practical_p_liq_dew,
                             d_ln_phi_dK_liq,k_i_dew,z,p_dew,n);
                B_linar_calc_dew(B_linar_dew,phi_liq_i_dew,phi_gas_i_dew,k_i_dew,z,n);
                A_linar_calc_bub(A_linar_bub,practical_ln_phi_practical_p_gas_bub,practical_ln_phi_practical_p_liq_bub,
                             d_ln_phi_dK_gas,k_i_bub,z,p_bub,n);
                B_linar_calc_bub(B_linar_bub,phi_liq_i_bub,phi_gas_i_bub,k_i_bub,z,n);
                
                cout<<"A_linar_bub\n";
                for(i=0;i<=n;i++)
                {
                    for(j=0;j<=n;j++)
                    {
                        cout<<A_linar_bub[i][j]<<"   ";
                    }cout<<"\n";
                }
                char mmm;
                cin>>mmm;
                Gauss(A_linar_dew,B_linar_dew,x_solve_dew,n+1);
                Gauss(A_linar_bub,B_linar_bub,x_solve_bub,n+1);
                
//cout<<"K_calk";
                K_calk(k_i_dew,x_solve_dew,n);
                K_calk(k_i_bub,x_solve_bub,n);
                
             p_dew=p_dew*expq((__float128)x_solve_dew[n]);
                p_bub=p_bub*expq((__float128)x_solve_bub[n]);}

                count_step_interation++;
            
            kriterij_dew=0.0;
            kriterij_bub=0.0;
    
            for(int k=0;k<n;k++)
            {
                kriterij_dew+=cabsq((__float128)((x_dew[k]* phi_liq_i_dew[k] )/(y_dew[k]* phi_gas_i_dew[k]) -1.0));
                kriterij_bub+=cabsq((__float128)((x_bub[k]* phi_liq_i_bub[k] )/(y_bub[k]* phi_gas_i_bub[k]) -1.0));
            }
            cout<<count_step_interation<<" "<<kriterij_dew<<" "<<T<<"\n\n";
            cout<<"Pdew "<<p_dew<<"  Pbub  "<<p_bub<<"\n\n";
            count_step_interation++;
        }while(/*((kriterij_dew!=kriterij_dew)||((kriterij_dew+kriterij_bub)>EPSILON))&&*/(count_step_interation<100));
        data_final<<"    "<< T<<"    "<<p_dew<<"    "<<p_bub<<'\n';


    }
data_final_k_i.close();
data_final.close();
 /*  free(dz_dc_liq   );
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
        free(d_ln_phi_dK_gas );*/



    return 0;
}
