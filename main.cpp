// Genetic Algorithm for VRP
// Written by Microsoft Visual C++
// Copyright by UTLab @ Tsinghua University
// http://orsc.edu.cn/UTLab

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "UTLab.h"
#include <time.h>
#include <string.h>
#define node  18   
#define vehicle 4 
#define M 1 
#define TYPE -1 
#define GEN 5000
#define POP_SIZE 30    
#define P_MUTATION 0.6   
#define P_CROSSOVER 0.4  
#define Num 5000 
#define alpha 0.9 
//#define N 60

static void  initialization(void);      
static void  evaluation(int gen);     
static void  selection(void);          
static void  crossover(void);           
static void  mutation(void);            
static void  objective_function(void);  
static int constraint_check(int x_temp[node+1], int y_temp[vehicle+1], double t_temp[vehicle+1]);
int rdint(int a,int b);                 
void sortminmax(int p[vehicle+1],int min,int max);      
double myu(double a, double b); 
static double Cr(int x[node+1], int y[vehicle+1], double t[vehicle+1]);  




double  OBJECTIVE[POP_SIZE+1][M+1];    
double  q[POP_SIZE+1];           
int x[POP_SIZE+1][node+1];       
int y[POP_SIZE+1][vehicle+1];    
double t[POP_SIZE+1][vehicle+1];  


double d[node+1][node+1]=    
{ 
0.00, 19.0, 17.5, 28.0, 24.0, 24.5, 31.2, 31.0, 21.0, 18.0, 21.5, 36.5, 31.5, 23.0, 28.0, 34.5, 30.0, 18.5, 24.0,
19.0, 0.00, 6.00, 11.0, 21.0, 32.0, 44.5, 48.5, 37.5, 36.0, 40.0, 55.0, 46.5, 38.5, 38.5, 40.0, 29.5, 16.5, 14.0,
17.5, 6.00, 0.00, 10.5, 15.0, 26.0, 39.5, 45.0, 33.5, 33.0, 39.0, 54.0, 48.0, 39.0, 41.5, 44.0, 34.5, 21.3, 20.0,
28.0, 11.0, 10.5, 0.00, 20.0, 34.0, 49.0, 55.5, 44.0, 44.0, 49.5, 65.0, 57.0, 48.5, 49.0, 50.0, 38.0, 26.5, 22.0,
24.0, 21.0, 15.0, 20.0, 0.00, 15.5, 31.0, 41.5, 30.0, 38.0, 43.0, 56.0, 55.5, 47.0, 52.0, 56.5, 48.5, 35.0, 35.0,
24.5, 32.0, 26.0, 34.0, 15.5, 0.00, 16.0, 28.5, 18.5, 24.0, 36.5, 46.0, 51.0, 44.0, 51.0, 58.5, 54.0, 41.0, 44.0,
31.2, 44.5, 39.5, 49.0, 31.0, 16.0, 0.00, 16.0, 13.0, 20.0, 32.5, 37.0, 48.5, 43.0, 52.5, 62.0, 60.5, 50.0, 54.5,
31.0, 48.5, 45.0, 55.5, 41.5, 28.5, 16.0, 0.00, 11.5, 13.5, 21.0, 21.0, 36.0, 32.5, 43.5, 54.0, 56.0, 48.0, 55.0,
21.0, 37.5, 33.5, 44.0, 30.0, 18.5, 13.0, 11.5, 0.00, 7.00, 20.0, 28.0, 36.0, 30.0, 40.0, 50.0, 49.0, 39.0, 45.0,
18.0, 36.0, 33.0, 44.0, 38.0, 24.0, 20.0, 13.5, 7.00, 0.00, 13.0, 23.0, 29.0, 23.0, 33.0, 43.0, 43.5, 35.0, 41.5,
21.5, 40.0, 39.0, 49.5, 43.0, 36.5, 32.5, 21.0, 20.0, 13.0, 0.00, 15.5, 16.0, 12.0, 22.5, 34.0, 38.0, 33.0, 41.0,
36.5, 55.0, 54.0, 65.0, 56.0, 46.0, 37.0, 21.0, 28.0, 23.0, 15.5, 0.00, 21.5, 23.0, 33.0, 44.5, 51.5, 48.0, 56.5,
31.5, 46.5, 48.0, 57.0, 55.5, 51.0, 48.5, 36.0, 36.0, 29.0, 16.0, 21.5, 0.00, 9.00, 13.0, 24.0, 33.0, 34.0, 43.0,
23.0, 38.5, 39.0, 48.5, 47.0, 44.0, 43.0, 32.5, 30.0, 23.0, 12.0, 23.0, 9.00, 0.00, 11.0, 22.0, 28.0, 27.0, 35.0,
28.0, 38.5, 41.5, 49.0, 52.0, 51.0, 52.5, 43.5, 40.0, 33.0, 22.5, 33.0, 13.0, 11.0, 0.00, 11.5, 20.5, 24.0, 32.0,
34.5, 40.0, 44.0, 50.0, 56.5, 58.5, 62.0, 54.0, 50.0, 43.0, 34.0, 44.5, 24.0, 22.0, 11.5, 0.00, 14.0, 23.5, 36.0,
30.0, 29.5, 34.5, 38.0, 48.5, 54.0, 60.5, 56.0, 49.0, 43.5, 38.0, 51.5, 33.0, 28.0, 20.5, 14.0, 0.00, 13.5, 17.0,
18.5, 16.5, 21.3, 26.5, 35.0, 41.0, 50.0, 48.0, 39.0, 35.0, 33.0, 48.0, 34.0, 27.0, 24.0, 23.5, 13.5, 0.00, 8.50,
24.0, 14.0, 20.0, 22.0, 35.0, 44.0, 54.5, 55.0, 45.0, 41.5, 41.0, 56.5, 43.0, 35.0, 32.0, 36.0, 17.0, 8.50, 0.00,
};        
//The distances between nodes

double tw[node+1][2]=
{9.0,			 10.0,
 9,			 15+5.0/6,
 9+2.0/6,    15.5,
 9+4.0/6,    14+4.0/6,
 9+2.0/6,    14.5,
 9,			 15+2.0/6,
 9,			 14+2.0/6,
 9.5,		 14+1.0/6,
 9,			 15.5,
 9+1.0/6,    15+5.0/6,
 9+4.0/6,    13+2.0/6,
 9+1.0/6,    14+2.0/6,
 9,			 15+2.0/6,
 9+2.0/6,    16.5,
 9+2.0/6,    16,
 9+2.0/6,    14.5,
 9,			 15+1.0/6,
 9,			 16+2.0/6,
 9,			 16+2.0/6,
};
//The time windows of nodes.  
/*
double tw[node+1][2]=
{9.0,			 10.0,
 11,			 15+5.0/6,
 11+2.0/6,    15.5,      //11+2.0/6,    15.5,为experiment2.txt结果
 9+4.0/6,    12+4.0/6,
 9+2.0/6,    12.5,
 12,			 15+2.0/6,
 11,			 14+2.0/6,
 9.5,		 14+1.0/6,
 11,			 15.5,
 12+1.0/6,    15+5.0/6,
 9+4.0/6,    13+2.0/6,
 9+1.0/6,    12+2.0/6,
 9,			 15+2.0/6,
 12+2.0/6,    16.5,
 12+2.0/6,    16,
 9+2.0/6,    12.5,
 9,			 12+1.0/6,
 13,			 16+2.0/6,
 13,			 16+2.0/6,
};
*/

//The time windows of nodes. 缩小时间窗


double tt[node+1][node+1][3]=   
{
00, 00, 00,		25, 50, 75,		05, 10, 15, 	25, 50, 75,		07, 15, 23, 	25, 50, 75, 	25, 50, 75, 	12, 25, 38, 	07, 15, 23, 	 25, 50, 75, 	 10, 20, 30, 	25, 50, 75, 	27, 55, 83, 	05, 10, 15, 	25, 50, 75, 	22, 45, 68, 	07, 15, 23, 	15, 30, 45, 	25, 50, 75, 
25, 50, 75,		00, 00, 00, 	20, 40, 60, 	05, 10, 15, 	25, 50, 75, 	17, 35, 53, 	07, 15, 23, 	20, 40, 60, 	20, 40, 60,		 07, 15, 23,	 22, 45, 68, 	05, 10, 15, 	17, 35, 53, 	20, 40, 60, 	05, 10, 15, 	05, 10, 15, 	22, 45, 68, 	20, 40, 60, 	05, 10, 15,
05, 10, 15, 	20, 40, 60, 	00, 00, 00, 	20, 40, 60, 	07, 15, 23, 	17, 35, 53, 	20, 40, 60, 	15, 30, 45, 	05, 10, 15,		 22, 45, 68,	 12, 25, 38, 	17, 35, 53, 	17, 35, 53, 	05, 10, 15, 	20, 40, 60, 	20, 40, 60, 	07, 15, 23, 	12, 25, 38, 	22, 45, 68,
25, 50, 75,		05, 10, 15, 	20, 40, 60, 	00, 00, 00, 	22, 45, 68, 	15, 30, 45, 	02, 05, 8, 		17, 35, 53, 	22, 45, 68,		 05, 10, 15,	 22, 45, 68, 	15, 30, 45, 	15, 30, 45, 	20, 40, 60, 	02, 05, 8, 		05, 10, 15, 	22, 45, 68, 	20, 40, 60, 	05, 10, 15,
07, 15, 23, 	25, 50, 75, 	07, 15, 23, 	22, 45, 68, 	00, 00, 00, 	17, 35, 53, 	22, 45, 68, 	07, 15, 23, 	10, 20, 30,		 22, 45, 68,	 07, 15, 23, 	17, 35, 53, 	17, 35, 53, 	07, 15, 23, 	22, 45, 68, 	22, 45, 68, 	10, 20, 30, 	10, 20, 30, 	25, 50, 75,
25, 50, 75, 	17, 35, 53, 	17, 35, 53, 	15, 30, 45, 	17, 35, 53, 	00, 00, 00, 	15, 30, 45, 	12, 25, 38, 	17, 35, 53,		 15, 30, 45,	 15, 30, 45, 	05, 10, 15, 	02, 05, 8, 		15, 30, 45, 	15, 30, 45, 	15, 30, 45, 	15, 30, 45, 	12, 25, 38, 	15, 30, 45,
25, 50, 75, 	07, 15, 23, 	20, 40, 60, 	02, 05, 8, 		22, 45, 68, 	15, 30, 45, 	00, 00, 00, 	17, 35, 53, 	20, 40, 60,		 05, 10, 15,	 20, 40, 60, 	15, 30, 45, 	15, 30, 45, 	17, 35, 53, 	02, 05, 8, 		05, 10, 15, 	22, 45, 68, 	17, 35, 53, 	07, 15, 23,
12, 25, 38, 	20, 40, 60, 	15, 30, 45, 	17, 35, 53, 	07, 15, 23, 	12, 25, 38, 	17, 35, 53, 	00, 00, 00, 	17, 35, 53,		 20, 40, 60,	 05, 10, 15, 	05, 10, 15, 	07, 15, 23, 	17, 35, 53, 	17, 35, 53, 	17, 35, 53, 	17, 35, 53, 	02, 05, 8, 	17, 35, 53,
07, 15, 23, 	20, 40, 60, 	05, 10, 15, 	22, 45, 68, 	10, 20, 30, 	17, 35, 53, 	20, 40, 60, 	17, 35, 53, 	00, 00, 00,		 20, 40, 60,	 12, 25, 38, 	17, 35, 53, 	17, 35, 53, 	05, 10, 15, 	17, 35, 53, 	17, 35, 53, 	10, 20, 30, 	12, 25, 38, 	20, 40, 60,
25, 50, 75, 	07, 15, 23, 	22, 45, 68, 	05, 10, 15, 	22, 45, 68, 	15, 30, 45, 	05, 10, 15, 	20, 40, 60, 	20, 40, 60,		 00, 00, 00,	 22, 45, 68, 	17, 35, 53, 	17, 35, 53, 	20, 40, 60, 	05, 10, 15, 	02, 05, 8, 		22, 45, 68, 	20, 40, 60, 	07, 15, 23,
10, 20, 30, 	22, 45, 68, 	12, 25, 38, 	22, 45, 68, 	07, 15, 23, 	15, 30, 45, 	20, 40, 60, 	05, 10, 15, 	12, 25, 38,		 22, 45, 68,	 00, 00, 00, 	15, 30, 45, 	12, 25, 38, 	20, 40, 60, 	20, 40, 60, 	20, 40, 60, 	15, 30, 45, 	05, 10, 15, 	20, 40, 60,
25, 50, 75, 	05, 10, 15, 	17, 35, 53, 	15, 30, 45, 	17, 35, 53, 	05, 10, 15, 	15, 30, 45, 	05, 10, 15, 	17, 35, 53,		 17, 35, 53,	 15, 30, 45, 	00, 00, 00, 	07, 15, 23, 	17, 35, 53, 	15, 30, 45, 	15, 30, 45, 	17, 35, 53, 	12, 25, 38, 	15, 30, 45,
27, 55, 83, 	17, 35, 53, 	17, 35, 53, 	15, 30, 45, 	17, 35, 53, 	02, 05, 8, 		15, 30, 45, 	07, 15, 23, 	17, 35, 53,		 17, 35, 53,	 12, 25, 38, 	07, 15, 23, 	00, 00, 00, 	17, 35, 53, 	15, 30, 45, 	20, 40, 60, 	07, 15, 23, 	12, 25, 38, 	20, 40, 60,
05, 10, 15, 	20, 40, 60, 	05, 10, 15, 	20, 40, 60, 	07, 15, 23, 	15, 30, 45, 	17, 35, 53, 	17, 35, 53, 	05, 10, 15,		 20, 40, 60,	 20, 40, 60, 	17, 35, 53, 	17, 35, 53, 	00, 00, 00, 	20, 40, 60, 	20, 40, 60, 	07, 15, 23, 	12, 25, 38, 	20, 40, 60,
25, 50, 75, 	05, 10, 15, 	20, 40, 60, 	02, 05, 8, 		22, 45, 68, 	15, 30, 45, 	02, 05, 8, 		17, 35, 53, 	17, 35, 53,		 05, 10, 15,	 20, 40, 60, 	15, 30, 45, 	15, 30, 45, 	20, 40, 60, 	00, 00, 00, 	02, 05, 8, 		22, 45, 68, 	17, 35, 53, 	17, 35, 53,
22, 45, 68, 	05, 10, 15, 	20, 40, 60, 	05, 10, 15, 	22, 45, 68, 	15, 30, 45, 	05, 10, 15, 	17, 35, 53, 	17, 35, 53,		 02, 05, 8,		 20, 40, 60, 	15, 30, 45, 	20, 40, 60, 	20, 40, 60, 	02, 05, 8, 		00, 00, 00, 	22, 45, 68, 	17, 35, 53, 	07, 15, 23,
07, 15, 23, 	22, 45, 68, 	07, 15, 23, 	22, 45, 68, 	10, 20, 30, 	15, 30, 45, 	22, 45, 68, 	17, 35, 53, 	10, 20, 30,		 22, 45, 68,	 15, 30, 45, 	17, 35, 53, 	07, 15, 23, 	07, 15, 23, 	22, 45, 68, 	22, 45, 68, 	00, 00, 00, 	12, 25, 38, 	07, 15, 23,
15, 30, 45, 	20, 40, 60, 	12, 25, 38, 	20, 40, 60, 	10, 20, 30, 	12, 25, 38, 	17, 35, 53, 	02, 05, 8, 		12, 25, 38,		 20, 40, 60,	 05, 10, 15, 	12, 25, 38, 	12, 25, 38, 	12, 25, 38, 	17, 35, 53, 	17, 35, 53, 	12, 25, 38, 	00, 00, 00, 	20, 40, 60,
25, 50, 75, 	05, 10, 15, 	22, 45, 68, 	05, 10, 15, 	25, 50, 75, 	15, 30, 45, 	07, 15, 23, 	17, 35, 53, 	20, 40, 60,		 07, 15, 23,	 20, 40, 60, 	15, 30, 45, 	20, 40, 60, 	20, 40, 60, 	07, 15, 23, 	07, 15, 23, 	20, 40, 60, 	20, 40, 60, 	00, 00, 00,
};   
// The travelling times between Node i and Node j, denoted by triangular fuzzy numbers. 节点i和j之间的行驶时间


double S[node+1]=     
{ 0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,0.25, 0.25, 0.25, 0.25,0.25, 0.25, 0.25, 0.25,0.25, 0.25, 0.25,};        
//The unloading times for customers.客户卸货时间（已知固定15mins）

int Q[vehicle+1]=     
{ 0, 1000, 1000, 1000, 1000,
};        
//The capacities of vehicles.车辆容量

double qq[node+1]=     
{ 0, 200, 100, 140, 160, 200, 60, 200, 135, 120, 140, 100, 200, 80, 60, 200, 90, 200, 100, 
};        


static double max(double a,double b)
{
 if(a>b)
  return a;
 else 
  return b;
}

static void objective_function(void) 
{
	 int x_temp[node+1], y_temp[vehicle+1]; 
	double g[vehicle+1]; 
	double g_temp=0;
	double temp;      
	int i, j, k;
	for(i = 1; i <=POP_SIZE; i++) 
	{
		g_temp=0;     //20160806-2
		for(j=0; j<=node; j++)
		{
			x_temp[j]=x[i][j];
		}
		for(k=0;k<=vehicle;k++)
		{
			y_temp[k]=y[i][k];
		}
		for(k=1;k<=vehicle;k++)
		{
			if(y_temp[k]==y_temp[k-1]) 
			{
				g[k]=0;
			}
			else     
			{
				temp=d[0][x_temp[y_temp[k-1]+1]]+d[x_temp[y_temp[k]]][0];
				for(j=(y_temp[k-1]+1);j<=(y_temp[k]-1);j++)
				{
					temp=temp+d[x_temp[j]][x_temp[j+1]];
				}
				g[k]=temp;
			}                                                                   
			g_temp=g_temp+g[k]; 
		}
		OBJECTIVE[i][1]=g_temp;        
		OBJECTIVE[i][0]=OBJECTIVE[i][1]; 
	}
}

static int constraint_check(int x_temp[node+1], int y_temp[vehicle+1], double t_temp[vehicle+1]) 
{
	double q_temp;
	int k, j; 
	for(k=1; k<=vehicle; k++)
	{
		q_temp=0;
		for(j=y_temp[k-1]+1; j<=y_temp[k]; j++)
		{
			q_temp=q_temp+qq[x_temp[j]];  
		}
		if((q_temp>Q[k])||(t_temp[k]<tw[0][0])||(t_temp[k]>tw[0][1])) return 0;
	}
	if(Cr(x_temp, y_temp, t_temp)<alpha)
	return 0;
	else
	return 1;
}


static void initialization(void)
{
	int x_temp[node+1], y_temp[vehicle+1];
	double t_temp[vehicle+1]; 
	int temp;
	int i,j,k;
	for(i=0; i<=POP_SIZE; i++)
	{
		mark:              
	    for(j=1;j<=node;j++)
		{
			x_temp[j]=j;
		}
		for(j=1;j<=node;j++)  
		{
			k=rdint(j,node);  
			temp=x_temp[j];
			x_temp[j]=x_temp[k];
			x_temp[k]=temp;
		}
		for(k=1;k<vehicle;k++)  
		{
			y_temp[k]=rdint(0,node);
		}
		y_temp[0]=0;
		y_temp[vehicle]=node;    
		sortminmax(y_temp,1,vehicle-1);  
		for(k=0;k<=vehicle;k++)    
			y[i][k]=y_temp[k];     
		
		for(k=1;k<=vehicle;k++)    
		
		{
			if(y_temp[k-1]==y_temp[k])
				t_temp[k]=0;
			else
				t_temp[k]=myu(tw[0][0],tw[0][1]);  
		}
		if(constraint_check(x_temp,y_temp,t_temp)==0) goto mark; 
		else
		{
			for(j=1;j<=node;j++)
			{
				x[i][j]=x_temp[j];	
			}
			for(k=0;k<=vehicle;k++)
			{
				y[i][k]=y_temp[k];
				t[i][k]=t_temp[k]; 
			}
		}	
	}
}

main()
{
	clock_t start,end;  
	start = clock();    
 
  int i, j, k;
  double a;
  srand(time(0));
  for(i=0;i<=node;i++)
	  for(j=0;j<=node;j++)
		  for(k=0;k<=2;k++)
		  tt[i][j][k]=tt[i][j][k]/60;
  q[0]=0.05; a=0.05;     
  for(i=1; i<=POP_SIZE; i++) {a=a*0.95; q[i]=q[i-1]+a;}   
  initialization();
  evaluation(0);
//FILE *fp = fopen("1.txt", "a+");//a+表示下次会接着上次的结果写
  FILE *fp = fopen("17.txt", "w");//

 if (fp==0) { printf("can't open file\n"); return 0;}
 for(i=1; i<=GEN; i++) 
 {
	 selection(); 
	 crossover();   
	 mutation();   
	 evaluation(i);
	 fprintf(fp,"\nGeneration NO.%d  ",i);
//	 printf("\nGeneration NO.%d\n", i);
	 fprintf(fp,"%3.4f\n",OBJECTIVE[0][1]);
//	 printf("%3.4f\n",OBJECTIVE[0][1]);
 }
 	 	end = clock();  
	printf("time=%f\n",(double)(end-start)/CLK_TCK); 
	fprintf(fp,"time=%f\n",(double)(end-start)/CLK_TCK); 
	  printf("x=(");
	  fprintf(fp,"x=(");
	  for(j=1; j<=node; j++) {
		  if(j<node) {printf("%d,",x[0][j]);fprintf(fp,"%d,",x[0][j]);}
		  else {printf("%d",x[0][j]);fprintf(fp,"%d",x[0][j]);}
	  }
	  printf(")\n");fprintf(fp,")\n");
	  printf("y=(");fprintf(fp,"y=(");
	  for(k=0; k<=vehicle; k++) {
		  if(k<vehicle){ printf("%d,",y[0][k]);fprintf(fp,"%d,",y[0][k]);}
		  else {printf("%d",y[0][k]);fprintf(fp,"%d",y[0][k]);}
	  }
	  printf(")\n");fprintf(fp,")\n");
	  printf("t=(");fprintf(fp,"t=(");
	  for(k=1; k<=vehicle; k++) {
		  if(k<vehicle) {printf("%3.4f,",t[0][k]);fprintf(fp,"%3.4f,",t[0][k]);}
		  else {printf("%3.4f",t[0][k]);fprintf(fp,"%3.4f",t[0][k]);}
	  }
	  if(M==1){ printf(")\nf=%3.4f\n", OBJECTIVE[0][1]);fprintf(fp,")\nf=%3.4f\n", OBJECTIVE[0][1]);}
	  else {
	      printf(")\nf=(");fprintf(fp,")\nf=(");
	      for(j=1; j<=M; j++) {
			  if(j<M) {printf("%3.4f,", OBJECTIVE[0][j]);fprintf(fp,"%3.4f,", OBJECTIVE[0][j]);}
			  else {printf("%3.4f", OBJECTIVE[0][j]); fprintf(fp,"%3.4f", OBJECTIVE[0][j]);}
		  }
          printf(")  Aggregating Value=%3.4f\n",OBJECTIVE[0][0]);
		  fprintf(fp,")  Aggregating Value=%3.4f\n",OBJECTIVE[0][0]);
	  
  }
  printf("\n");fprintf(fp,"\n");
  fclose(fp);
  return 1;  
}

static void evaluation(int gen)
{
  double a;
  int   i, j, k, label,b;
  objective_function();
 if(gen==0)
  {
	 for(k=0; k<=M; k++) OBJECTIVE[0][k]=OBJECTIVE[1][k];
//	 for(j = 1; j <= node; j++) x[0][j]=x[1][j];
//	 for(k = 0; k <vehicle; k++) 
//	{
//		y[0][k]=y[1][k];
//		t[0][k]=t[1][k]; 
//	 }  
  }
  for(i=0; i<POP_SIZE; i++)
  {
	  label=0;  a=OBJECTIVE[i][0];
	  for(j=i+1; j<=POP_SIZE; j++)
		  if((TYPE*a)<(TYPE*OBJECTIVE[j][0])) 
		  {  //选择最小的g(x,y)
			 a=OBJECTIVE[j][0];
			 label=j;
		  }
		  
		  if(label!=0) {
		  for(k=0; k<=M; k++) {  //将objective[i][k]和objective[label][k]进行交换
			 a=OBJECTIVE[i][k];
			 OBJECTIVE[i][k]=OBJECTIVE[label][k];
			 OBJECTIVE[label][k]=a;
		 }
		  for(j=1; j<=node; j++) {  //将x[i][j]和x[label][j]进行交换
			 a=x[i][j];
			 x[i][j]=x[label][j];
			 x[label][j]=a;
		 }
		  for(j=1; j<=vehicle; j++) { //将y[i][j]和y[label][j]进行交换
			 b=y[i][j];
			 y[i][j]=y[label][j];
			 y[label][j]=b;
		 }
		  for(j=1; j<=vehicle; j++) {  //将t[i][j]和t[label][j]进行交换
			 a=t[i][j];
			 t[i][j]=t[label][j];
			 t[label][j]=a;
		 }
		  	 
	  }
  }
}

static void selection()  
{
  double r;
  //double temp[POP_SIZE+1][N+1];
  int tempx[POP_SIZE+1][node+1],tempy[POP_SIZE+1][vehicle+1];
  double tempt[POP_SIZE+1][vehicle+1];     //老师没有定义数组，这是后面补上的
  int   i, j, k;
  for(i=1; i<=POP_SIZE; i++) {
	  r=myu(0, q[POP_SIZE]);  //从0到q[30]之间产生一个随机值
	  for(j=0; j<=POP_SIZE; j++) 
	  {
		  if(r<=q[j]) 
		  {
			  for(k=1; k<=node; k++) tempx[i][k]=x[j][k];
			  for(k=0; k<=vehicle; k++) tempy[i][k]=y[j][k];
			  for(k=1; k<=vehicle; k++) tempt[i][k]=t[j][k];
			  break;
		  }
	  }//tempx,y,t相当于把q[i]<r的部分删掉，重新给x,y,t赋值。
  }
  for(i=1; i<=POP_SIZE; i++)
  {
	  for(k=1; k<=node; k++) x[i][k]=tempx[i][k];
	  for(k=0; k<=vehicle; k++) y[i][k]=tempy[i][k];
	  for(k=1; k<=vehicle; k++) t[i][k]=tempt[i][k];
//	   printf("%d  ",y[i][4]);
  }
 
}

static void crossover()
{
	int i,j,jj,k,pop;
	int x_temp1[node+1],x_temp2[node+1],y_temp1[vehicle+1],y_temp2[vehicle+1];
	double t_temp1[vehicle+1],t_temp2[vehicle+1];
	double r;
	pop=POP_SIZE/2;
	for(i=1;i<=pop;i++) 
	{
		if(myu(0,1)>P_CROSSOVER) continue;
		j=(int)myu(1,POP_SIZE);
		jj=(int)myu(1,POP_SIZE);  
		r=myu(0,1);
		for(k=0;k<=vehicle;k++)
		{
			t_temp1[k]=r*t[j][k]+(1-r)*t[jj][k];
			t_temp2[k]=(1-r)*t[j][k]+r*t[jj][k];
		}
		for(k=0;k<=node;k++)
		{
			x_temp1[k]=x[j][k];  
			x_temp2[k]=x[jj][k];
		}
		for(k=0;k<=vehicle;k++)
		{
			y_temp1[k]=y[jj][k];
			y_temp2[k]=y[j][k]; 
		}

		if(constraint_check(x_temp1, y_temp1, t_temp1)==1)  
		{
			for(k=1;k<=node;k++)
			{
				x[j][k]=x_temp1[k];
			}
			for(k=1;k<=vehicle;k++)
			{
				y[j][k]=y_temp1[k];
				t[j][k]=t_temp1[k];
			}
		}
		if(constraint_check(x_temp2, y_temp2, t_temp2)==1)
		{
			for(k=1;k<=node;k++)
			{
				x[jj][k]=x_temp2[k];
			}
			for(k=1;k<=vehicle;k++)
			{
				y[jj][k]=y_temp2[k];
				t[jj][k]=t_temp2[k];
			}
		}
	}
}



static void mutation()   //对染色体进行变异处理
{	
	int i,j,k,r,a;
	int x_temp[node+1],y_temp[vehicle+1],direction[vehicle+1];
	double t_temp[vehicle+1];
	double temp;
	for(i=1;i<=POP_SIZE;i++) 
	{
		if(myu(0,1)>P_MUTATION) continue;  
		for(k=0;k<=vehicle;k++)
		{
			if(myu(0,1)<0.5) 
				direction[k]=(int)myu(-2,2);
			else
				direction[k]=0;
		}
		for(j=1;j<=vehicle;j++) 
		{
			t_temp[j]=t[i][j]+direction[j];  
		}
		r=(int)myu(1,node+1);  
		for(j=1;j<=r;j++)
		{
			k=rdint(j,r);    
			a=x[i][j];
			x[i][j]=x[i][k];
			x[i][k]=a;   
		}
		for(j=1;j<=node;j++)
			x_temp[j]=x[i][j];  
		r=(int)myu(1,vehicle+1);
		for(j=0;j<=vehicle;j++)
		{
			if(j<=r)  y_temp[j]=rdint(0,node); 
			else 
				y_temp[j]=y[i][j];
		}
		y_temp[0]=0;y_temp[4]=18;
	//	sortminmax(y_temp,1,r);  
		if(constraint_check(x_temp,y_temp,t_temp)==1){
			for(k=1;k<=node;k++) x[i][k]=x_temp[k];
			for(k=0;k<=vehicle;k++) y[i][k]=y_temp[k];
			for(k=1;k<=vehicle;k++) t[i][k]=t_temp[k];
			break;
		}
	}
}

//***************************************************************************************
// *******************************cr模拟改进——二分法***********************************
//***************************************************************************************


static double Cr(int x[node+1], int y[vehicle+1], double t[vehicle+1])//主要分析函数的可信度是否满足条件，是约束条件之一
{
	double r=1;
	int k,j;
	double pr[node+1][node+1]; //保存phi逆的值
	double f[node+1];
	double a=0.9;
	for(k=1;k<=vehicle; k++)
	{
		if(y[k]>y[k-1])
		{
			pr[0][x[y[k-1]+1]]=(2*tt[0][x[y[k-1]+1]][2]-2*tt[0][x[y[k-1]+1]][1])*a
										+2*tt[0][x[y[k-1]+1]][1]-tt[0][x[y[k-1]+1]][2];
			f[x[y[k-1]+1]]=t[k]+pr[0][x[y[k-1]+1]];
			if((f[x[y[k-1]+1]]>tw[x[y[k-1]+1]][1]))
			{
				r=0;
				goto mark0;
			}
			for(j=2;j<=(y[k]-y[k-1]);j++)
			{
				pr[x[y[k-1]+j-1]][x[y[k-1]+j]]=2*a*(tt[x[y[k-1]+j-1]][x[y[k-1]+j]][2]-tt[x[y[k-1]+j-1]][x[y[k-1]+j]][1])
										+2*tt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]-tt[x[y[k-1]+j-1]][x[y[k-1]+j]][2];
//				printf("%f\n",pr[x[y[k-1]+j-1]][x[y[k-1]+j]]);
				f[x[y[k-1]+j]]=max(f[x[y[k-1]+j-1]], tw[x[y[k-1]+j-1]][0])+S[x[y[k-1]+j-1]]+pr[x[y[k-1]+j-1]][x[y[k-1]+j]];
				if((f[x[y[k-1]+j]]>tw[x[y[k-1]+j]][1]))
				{
					r=0;
					goto mark0;
				}
			}
		}
	}
mark0:
	return r;
}
//cr simulatuon_bisection
/*
static double Cr(int x[node+1], int y[vehicle+1], double t[vehicle+1])//主要分析函数的可信度是否满足条件，是约束条件之一
{
	double mu1=1.0 ,mu2=0.0;
	double tt_temp[node+1][node+1], f[node+1];
	double ttt[node+1][node+1][3];
	double mutt_temp[node+1][node+1];  
	int signal,s1,s2,i,k,j,m,n,r;
	signal=1;
	for(k=1;k<=vehicle; k++)
	{
		signal=1;
		if(y[k]>y[k-1])
		{
			tt_temp[0][x[y[k-1]+1]]=tt[0][x[y[k-1]+1]][1];
			f[x[y[k-1]+1]]=t[k]+tt_temp[0][x[y[k-1]+1]];
			if((f[x[y[k-1]+1]]>tw[x[y[k-1]+1]][1]))
			{
				signal=0;
				goto mark0;
			}
			for(j=2;j<=(y[k]-y[k-1]);j++)
			{
				tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]=tt[x[y[k-1]+j-1]][x[y[k-1]+j]][1];
				f[x[y[k-1]+j]]=max(f[x[y[k-1]+j-1]], tw[x[y[k-1]+j-1]][0])+S[x[y[k-1]+j-1]]+tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]];
				if((f[x[y[k-1]+j]]>tw[x[y[k-1]+j]][1]))
				{
					signal=0;
					goto mark0;
				}
			}
		}
	}
mark0:
	if(signal==0)
	{
		for(m=0;m<=node;m++)
			for(n=0;n<=node;n++)
				for(r=0;r<=2;r++)
					ttt[m][n][r]=tt[m][n][r];
		mu2=1;
		mu1=0;
		s1=1;
		for(i=1;i<=13;i++)
		{
			for(k=1;k<=vehicle; k++)
			{
				if(y[k]>y[k-1])
				{
					tt_temp[0][x[y[k-1]+1]]=(ttt[0][x[y[k-1]+1]][0]+ttt[0][x[y[k-1]+1]][1])/2;
					f[x[y[k-1]+1]]=t[k]+tt_temp[0][x[y[k-1]+1]];
					mu1=triangle(tt_temp[0][x[y[k-1]+1]], tt[0][x[y[k-1]+1]][0],tt[0][x[y[k-1]+1]][1], tt[0][x[y[k-1]+1]][2]);
					if((f[x[y[k-1]+1]]>tw[x[y[k-1]+1]][1]))	
					{   
						s1=0; 
						goto mark1;
					}
					for(j=2;j<=y[k]-y[k-1];j++)
					{
						tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][0]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1])/2;
						f[x[y[k-1]+j]]=max(f[x[y[k-1]+j-1]], tw[x[y[k-1]+j-1]][0])+S[x[y[k-1]+j-1]]+tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]];
						if((f[x[y[k-1]+j]]>tw[x[y[k-1]+j]][1]))  
						{  
							s1=0;  
							goto mark1;
						}
					}
				}
			}
mark1:
			if(s1==0)
			{
				for(k=1;k<=vehicle; k++)
				{
					if(y[k]>y[k-1])
					{
						ttt[0][x[y[k-1]+1]][1]=(ttt[0][x[y[k-1]+1]][0]+ttt[0][x[y[k-1]+1]][1])/2;
						for(j=2;j<=y[k]-y[k-1]; j++)
							ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][0]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1])/2;
					}
				}
			}
			if(s1==1)
			{
				for(k=1;k<=vehicle; k++)
				{
					if(y[k]>y[k-1])
					{
						ttt[0][x[y[k-1]+1]][0]=(ttt[0][x[y[k-1]+1]][0]+ttt[0][x[y[k-1]+1]][1])/2;
						for(j=2;j<=y[k]-y[k-1]; j++)
							ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][0]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][0]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1])/2;
					}
				}
			}
		}
	}
	if(signal==1)
	{
		for(m=0;m<=node;m++)
			for(n=0;n<=node;n++)
				for(r=0;r<=2;r++)
					ttt[m][n][r]=tt[m][n][r];
		mu2=0;
		mu1=1;
		s2=1;
		for(i=1;i<=13;i++)
		{
			for(k=1;k<=vehicle; k++)
			{
				if(y[k]>y[k-1])
				{
					tt_temp[0][x[y[k-1]+1]]=(ttt[0][x[y[k-1]+1]][1]+ttt[0][x[y[k-1]+1]][2])/2;
					f[x[y[k-1]+1]]=t[k]+tt_temp[0][x[y[k-1]+1]];
					mu2=triangle(tt_temp[0][x[y[k-1]+1]], tt[0][x[y[k-1]+1]][0],tt[0][x[y[k-1]+1]][1], tt[0][x[y[k-1]+1]][2]);
					if(f[x[y[k-1]+1]]>tw[x[y[k-1]+1]][1])	
					{   
						s2=0; 
						goto mark2;
					}
					for(j=2;j<=y[k]-y[k-1];j++)
					{
						tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][2])/2;
						f[x[y[k-1]+j]]=max(f[x[y[k-1]+j-1]], tw[x[y[k-1]+j-1]][0])+S[x[y[k-1]+j-1]]+tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]];
						if((f[x[y[k-1]+j]]>tw[x[y[k-1]+j]][1]))  
						{  
							s2=0;  
							goto mark2;
						}
					}
				}
			}
mark2:
			if(s2==0)
			{
				for(k=1;k<=vehicle; k++)
				{
					if(y[k]>y[k-1])
					{
						ttt[0][x[y[k-1]+1]][2]=(ttt[0][x[y[k-1]+1]][1]+ttt[0][x[y[k-1]+1]][2])/2;
						for(j=2;j<=y[k]-y[k-1]; j++)
							ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][2]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][2])/2;
					}
				}
			}
			if(s2==1)
			{
				for(k=1;k<=vehicle; k++)
				{
					if(y[k]>y[k-1])
					{
						ttt[0][x[y[k-1]+1]][1]=(ttt[0][x[y[k-1]+1]][1]+ttt[0][x[y[k-1]+1]][2])/2;
						for(j=2;j<=y[k]-y[k-1]; j++)
							ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]=(ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][1]+ttt[x[y[k-1]+j-1]][x[y[k-1]+j]][2])/2;
					}
				}
			}
		}
	}
	return((mu1+1-mu2)/2);
}
*/
/**/
//***************************************************************************************
// *******************************  cr模拟_原始方法   ***********************************
//***************************************************************************************

//cr simulation_original
/*
static double Cr(int x[node+1], int y[vehicle+1], double t[vehicle+1])//主要分析函数的可信度是否满足条件，是约束条件之一
{
	double mu;
	double tt_temp[node+1][node+1], f[node+1];
	double mutt_temp[node+1][node+1];  //这两个里面具体的下标不知道是多少啊？？？
	int signal,i,k,j;
	double v1=0.0, v2=0.0; //v1代表POS，v2代表NEC
      	for(i=1; i<=Num; i++)  
	{
		mu=1;
		signal=1;
		for(k=1;k<=vehicle;k++)
		{
			if(y[k]>y[k-1])  //tt[i][j][k]是车辆行驶时间的模糊变量
			{
				tt_temp[0][x[y[k-1]+1]]=myu(tt[0][x[y[k-1]+1]][0], tt[0][x[y[k-1]+1]][2]); //随机产生一个介于两者之间的数，这是出发时间
				mutt_temp[0][x[y[k-1]+1]]=triangle(tt_temp[0][x[y[k-1]+1]], tt[0][x[y[k-1]+1]][0], tt[0][x[y[k-1]+1]][1], tt[0][x[y[k-1]+1]][2]); //traingle主要求的是隶属度
				if(mu>mutt_temp[0][x[y[k-1]+1]]) mu=mutt_temp[0][x[y[k-1]+1]]; //第k辆车到达第一个用户的时间点
				f[x[y[k-1]+1]]=t[k]+tt_temp[0][x[y[k-1]+1]]; //没有定义，前面先定义一下，这是到达时间
				if((f[x[y[k-1]+1]]>tw[x[y[k-1]+1]][1]))
				{
					signal=0;
					goto mark;
				}
				for(j=2;j<=(y[k]-y[k-1]);j++) //第k辆车到达第j个用户的时间点
				{
					tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]=myu(tt[x[y[k-1]+j-1]][x[y[k-1]+j]][0], tt[x[y[k-1]+j-1]][x[y[k-1]+j]][2]);
					mutt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]=triangle(tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]], tt[x[y[k-1]+j-1]][x[y[k-1]+j]][0], tt[x[y[k-1]+j-1]][x[y[k-1]+j]][1], tt[x[y[k-1]+j-1]][x[y[k-1]+j]][2]);
					if(mu>mutt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]) mu=mutt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]]; //mu得到了最小的隶属度
					f[x[y[k-1]+j]]=max(f[x[y[k-1]+j-1]], tw[x[y[k-1]+j-1]][0])+S[x[y[k-1]+j-1]]+tt_temp[x[y[k-1]+j-1]][x[y[k-1]+j]];  
					if((f[x[y[k-1]+j]]>tw[x[y[k-1]+j]][1]))
					{
						signal=0;
						goto mark;
					}
				}
			}
		}
		mark:
		if((signal==0)&&(v2<mu))   v2=mu;  //这边等于是不是表判断，是不是应该两个=？？？	
		if((signal==1)&&(v1<mu))   v1=mu;	
	}
     	if(v1>v2) v1=1; else v2=1;
   	return (v1+1-v2)/2;
}
  
*/
int rdint(int a,int b)
{
	double x,temp2;
	int temp1;
	x=myu(a,b+0.99);      
	temp2=floor(x);   
	temp1=(int) temp2;   
	return(temp1);
}

void sortminmax(int p[vehicle+1],int min,int max)  //对p进行从小到大进行排序
{
  int i,j,w;
  for(i=min; i<max; i++)
    for(j=i+1; j<=max; j++)
       if(p[i] > p[j]) 
       {
	      w = p[i];
	      p[i] = p[j];
	      p[j] = w;
       }
}


