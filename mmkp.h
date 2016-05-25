#ifndef _MMKPH
#define _MMKPH
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include"ilcplex\ilocplex.h"
#define EPSILON (1e-6)
#define FABS(x) ((x)<0? -(x):(x))
#define MAX(a,b) ((a)>(b)? (a):(b))
#define MIN(a,b) ((a)<(b)? (a):(b))

int n, l, m;
typedef struct Group_st{
	int nItem;
	double *profit;
	int **consume;
}Group;

 Group *group;
 int *Gb;

 int n_var;
 int *sum_n_item;//(1,2,3,4) -> (0,1,3,6), sum_n_item[i] indicates the sum before i (excluding)
 int *bestSolVec;
 int *fixed_one; //fixed_one[i] gives if group i has fixed an item to 1.

 double UB;
 double LB;
 double *rdCosts;
 double *primal_sol;
 
typedef struct ReducedProblem_st {
	int len;
	int *posi;
}ReducedProblem;

ReducedProblem rp;

typedef struct ValuePosi{
	double value;
	int posi;
}ValuePosi;

ValuePosi *rdCost_sort;
ValuePosi *rdCost_sort1;

int len_rdCost0;
int len_rdCost1;

int **pre_v0;
int **pre_v1;
int *len_v0;
int *len_v1;

int *sum_pre_v1;
int maxIte;

void solveLP();
void solveIP(double time_limit, int ite, double gap);
void solveMIP(double);
int biased_selection();
void subproblem(double ration, double maxRdc);
 void initMemory(){
	UB=INT_MAX;
	LB=-1;
	rdCosts=new double[n_var];
	primal_sol=new double[n_var];
	bestSolVec=new int[n];
	rp.len=0;
	rp.posi=new int[n_var];

	for(int i=0;i<n_var;i++){
		rp.posi[rp.len++]=i;
	}

	fixed_one=new int[n];
	memset(fixed_one, 0, n*sizeof(int));
	rdCost_sort=new ValuePosi[n_var];
	len_rdCost0=0;

	rdCost_sort1=new ValuePosi[n_var];
	len_rdCost1=0;
	maxIte=10000;
}
 void deleteAll(){
	 delete []primal_sol;
	 delete []rdCosts;
	 delete []bestSolVec;
	 delete []rp.posi;
	 delete []fixed_one;
	 delete []rdCost_sort;
	 delete []rdCost_sort1;
	 delete []Gb;
	 delete []sum_n_item;
	 int i,j,k;
	 for(i=0;i<n;i++){
		 delete [] group[i].profit;
		 for(j=0;j<group[i].nItem;j++){
			 delete []group[i].consume[j];
		 }
		 delete []group[i].consume;
	 }
	 delete []group;
	 for(i=0;i<maxIte;i++){
		 delete []pre_v0[i];
		 delete []pre_v1[i];
	 }
	 delete []pre_v0;
	 delete []pre_v1;
	 delete []len_v0;
	 delete []len_v1;
	 delete []sum_pre_v1;
 }
 void read_inst(char *file){
	int i,j,k;
	freopen(file, "r",stdin);
	scanf("%d%d%d", &n, &l, &m);
	group=new Group[n];
	Gb=new int[m];

	for(i=0;i<m;i++){
		scanf("%d",&Gb[i]);
	}
	n_var=0;
	sum_n_item=new int[n];
	memset(sum_n_item, 0, n*sizeof(int));
	for(i=0;i<n;i++){
		scanf("%d", &j);
		n_var+=l;
		if(i>0) sum_n_item[i]=sum_n_item[i-1]+l;

		group[i].nItem=l;
		group[i].profit=new double[l];
		group[i].consume=new int*[l];

		for(j=0;j<l;j++){
			scanf("%lf",&group[i].profit[j]);
			group[i].consume[j]=new int[m];
			for(k=0;k<m;k++){
				scanf("%d",&group[i].consume[j][k]);
			}
		}
	}
	initMemory();
}
#endif
