#include <stdio.h>                  // C ���Ա�׼�������������
#include <math.h>                   // C ��ѧ������

#include "mathlib.h"                // DSP ��ѧ������
#include "dsplib.h"                 // DSP ������

/*************************�궨��*************************/

#define     PI      3.141592654
#define     Fs      100000.0		//����Ƶ��
#define     N       8			//��Ԫ����
#define     R       0.075		//Բ��뾶
#define     Beta    2*PI/N
#define     C		1500.0		//ˮ������

#define     x_l     600		//�źŽ�ȡ����
#define     b_l     100		//���˲��㳤��

#define     NF      1024		//FFT����

/*************************ȫ�ֱ���*************************/

#pragma DATA_ALIGN(CFFT, 8);
float CFFT[2*NF]={0};
#pragma DATA_ALIGN(Cw, 8);
float Cw[2*NF]={0};

double xin[N][1000]={0};			//��������
double x_d[N][2*b_l+x_l]={0};
double D[20]={0};

/*************************��������*************************/
double absd(double n);
int argmax(double *l, int n);
void bit_rev(float* x, int n);
void gen_w_r2(float* w, int n);
void FFT(double *sig_in, int n,  double *sig_afft);
double find_f(double *sig, int n);

int main(void)
{
	int i, j, t, T, n, delta=30;
	double theta, delay, sum_c=0, sum_r=0;
	double v1, v2, v3, a_max, b_max, phi, u, direction=0, f=0;
	double P[2]={0}, ks[2]={0}, W[N]={0};
	f = find_f(xin[0], 1000);
	T = 360/delta;
	for (t=0; t<T; t++)
	{
		theta = (-180 + t * delta) * PI / 180.0;
		for (i=0; i<N; i++)
		{
			delay = Fs * R * cos(theta - i * Beta)/C;
			int start_p = b_l + round(delay);
			int end_p = start_p + x_l;
			for (j=0; j<start_p; j++ )
				x_d[i][j] = 0;
			for (j=start_p; j<end_p; j++ )
				x_d[i][j] = xin[i][j-start_p];
			for (j=end_p; j<(2*b_l+x_l); j++ )
				x_d[i][j] = 0;

			double t_e = 2 * PI * i / N;
			double lambda = C / f;
			P[0] = R * cos(t_e);
			P[1] = R * sin(t_e);
			ks[0] = -2 * PI / lambda * cos(theta);
			ks[1] = -2 * PI / lambda * sin(theta);
			W[i] = cos(P[0] * ks[0] + P[1] * ks[1]) / N;
		}

		sum_r=0;
		for (i=0; i<2*b_l+x_l; i++)
		{
			sum_c=0;
			for(j=0; j<N; j++)
			{
				sum_c += x_d[j][i] * W[j];
			}
			sum_r += sum_c;
		}
		D[t] = absd(sum_r);
		sum_r = 0;
	}
	n = argmax(D, T);
	if (n == 0)
		n = n + 6;
	if (n == 11)
		n = n - 6;

	v1 = D[(n-1+T)%T];
	v2 = D[n];
	v3 = D[(n+1)%T];
	a_max = v1 - v3;
	b_max = v1 + v3 - 2 * v2;
	phi = a_max / b_max;
	u = -135 + n * delta;
	direction = u + phi/2.0*delta;
	if (direction > 90)
		direction = direction - 180;
	else if(direction < -90)
		direction = 180 - direction;
	printf("����ֵ��%f\n",direction);
	asm(" SWBP 0 ");
	return 0;
}

double absd(double n)
{
	if (n>=0)
		return n;
	else
		return -n;
}

int argmax(double *l, int n)
{
	int i, temp=0;
	double max = l[0];
	for (i=0; i<n; i++)
		if(l[i]>max)
		{
			max = l[i];
			temp = i;
		}
	return temp;
}

// ����λ��ת
void bit_rev(float* x, int n)
{
    int i,j,k;
    float rtemp,itemp;
    j=0;
    for(i=1;i<(n-1);i++)
    {
        k=n>>1;
        while(k<=j)
        {
            j-=k;
            k>>=1;
        }
        j+=k;
        if(i<j)
        {
            rtemp=x[j*2];
            x[j*2]=x[i*2];
            x[i*2]=rtemp;
            itemp=x[j*2+1];
            x[j*2+1]=x[i*2+1];
            x[i*2+1]=itemp;
        }
    }
}

// ������ת����
void gen_w_r2(float* w, int n)
{
    int i;

    float e=PI*2.0/n;

    for(i=0;i<(n>>1);i++)
    {
        w[2*i]=cossp(i*e);
        w[2*i+1]=sinsp(i*e);
    }
}

//���ٸ���Ҷ�任
void FFT(double *sig_in, int n,  double *sig_afft)
{
	/* input sig_in �ź�ʱ��
	 * input n �źų���
	 * output sig_afft �ź�Ƶ��ķ���
	 */
	unsigned int i;
	for(i=0;i<n;i++)
	{
		CFFT[2*i]=(float)sig_in[i];
		CFFT[2*i+1]=0;
	}
	for(i=2*n;i<2*NF;i++)
		CFFT[i]=0;
	gen_w_r2(Cw, NF);
	bit_rev(Cw, NF >> 1);
	DSPF_sp_cfftr2_dit(CFFT, Cw, NF);
	bit_rev(CFFT, NF);//FFT�������Ч��MATLAB��fft���
	for(i=0;i<NF;i++)
		sig_afft[i]=sqrt(CFFT[2*i] * CFFT[2*i] + CFFT[2*i+1] * CFFT[2*i+1]);
}

double find_f(double *sig, int n)
{
	double sig_afft[NF]={0}, fc;
	int m_p;
	FFT(sig, n, sig_afft);
	m_p = argmax(sig_afft, NF);
	fc = (double)Fs / (double)NF * m_p;
	return fc;
}
