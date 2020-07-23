#include <stdio.h>                  // C ���Ա�׼�������������
#include <math.h>                   // C ��ѧ������

#include "mathlib.h"                // DSP ��ѧ������
#include "dsplib.h"                 // DSP ������

/*************************�궨��*************************/

#define     PI      3.141592654
#define     Fs      200000		//����Ƶ��
#define     N       8			//��Ԫ����
#define     R       0.24		//Բ��뾶
#define     Beta    2*PI/N
#define     C		1500.0		//ˮ������

#define     x_l     600		//�źŽ�ȡ����
#define     b_l     100		//���˲��㳤��

/*************************ȫ�ֱ���*************************/

double xin[N][1000]={0};			//��������
double x_d[N][2*b_l+x_l]={0};
double W[N] = {1,1,1,1,1,1,1,1};
double D[20] = {0};

/*************************��������*************************/
void CRS(double *x, double *y, int l, int n);
double absd(double n);
int argmax(double *l, int n);

int main(void)
{
	int i, j, t, T, n, delta=30;
	double theta, delay, sum_c=0, sum_r=0;
	double v1, v2, v3, a_max, b_max, phi, u, direction=0;
	T = 360/delta;
	for (t=0; t<T; t++)
	{

		double w[N]={0};
		CRS(w, W, N, t);
		theta = (-180 + t * delta) * PI / 180.0;
		for (i=0; i<N; i++)
		{
			delay = Fs * R * cos(theta - i * Beta)/C;
			int start_p = b_l + round(delay);
			int end_p = start_p + x_l;
			for (j=0; j<start_p; j++ )
				x_d[i][j] = 0;
			for (j=start_p; j<end_p; j++ )
				x_d[i][j] = xin[i][j-start_p] * w[i];
			for (j=end_p; j<(2*b_l+x_l); j++ )
				x_d[i][j] = 0;
		}

		sum_r=0;
		for (i=0; i<2*b_l+x_l; i++)
		{
			sum_c=0;
			for(j=0; j<N; j++)
			{
				sum_c += x_d[j][i];
			}
			sum_r += absd(sum_c);
		}
		D[t] = sum_r;
		sum_r = 0;
	}
	n = argmax(D, T);
	v1 = D[n-1];
	v2 = D[n];
	v3 = D[n+1];
	a_max = v1 - v3;
	b_max = v1 + v3 - 2 * v2;
	phi = a_max / b_max;
	u = -180 + n * delta;
	direction = u + phi/2.0*delta;
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

void CRS(double *x, double *y, int l, int n)
{
	int i;
	for (i=0; i<l; i++)
	{
		x[(i+n)%l] = y[i];
	}
}
