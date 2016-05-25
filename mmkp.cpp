
#include "mmkp.h"
#include "./cpu_time.c"

void solveMILP(double time_limit){

}
int cmp(const void *p, const void *q){
	ValuePosi *a=(ValuePosi*)p;
	ValuePosi *b=(ValuePosi*)q;
	if(a->value > b->value) return -1;
	else if( a->value==b->value) return 0;
	else return 1;
}

void solveLP(){
	IloEnv env;
	IloModel model(env);
	IloNumVarArray x(env, n_var, 0,1);
	IloRangeArray con(env);
	IloExpr obj(env);
	int i,j,k;
	for(i=0;i<n;i++){
		for(j=0;j<group[i].nItem;j++){
			obj+=group[i].profit[j]*x[sum_n_item[i]+j];
		}
	}
	model.add(IloMaximize(env,obj));
	obj.end();
	for(k=0;k<m;k++){
		IloExpr expr(env);
		for(i=0;i<n;i++){
			for(j=0;j<group[i].nItem;j++){
				expr+=group[i].consume[j][k]*x[sum_n_item[i]+j];
			}
		}
		con.add(expr<=Gb[k]);
		expr.end();
	}
	for(i=0;i<n;i++){
		IloExpr expr(env);
		for(j=0;j<group[i].nItem;j++){
			expr+=x[sum_n_item[i]+j];
		}
		con.add(expr==1);
		expr.end();
	}
	
	model.add(con);
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads,1);
	if(cplex.solve()){
		printf("LP solved with new value: %lf \n\n", cplex.getObjValue());
		UB=cplex.getObjValue();
		len_rdCost0=len_rdCost1=0;
		for(i=0;i<n_var;i++) {
			primal_sol[i]=(double)cplex.getValue(x[i]);
			rdCosts[i]=FABS((double)cplex.getReducedCost(x[i]));
			if(primal_sol[i]<EPSILON && rdCosts[i]>EPSILON){// non-basic variables at 0
				rdCost_sort[len_rdCost0].value=rdCosts[i];
				rdCost_sort[len_rdCost0++].posi=i;
			} else if(primal_sol[i]>1-EPSILON && rdCosts[i]>EPSILON){// non-basic variables at 1
				rdCost_sort1[len_rdCost1].value=rdCosts[i];
				rdCost_sort1[len_rdCost1++].posi=i;
			}
		}
		qsort(rdCost_sort, len_rdCost0,sizeof(ValuePosi), cmp);
		qsort(rdCost_sort1, len_rdCost1, sizeof(ValuePosi), cmp);
		
		len_v1=new int[maxIte];
		len_v0=new int[maxIte];

		pre_v0=new int*[maxIte];
		pre_v1=new int*[maxIte];
		for(i=0;i<maxIte;i++){
			pre_v1[i]=new int[len_rdCost1];
			pre_v0[i]=new int[len_rdCost0];
		}
		sum_pre_v1=new int[maxIte];
		memset(len_v0, 0,maxIte*sizeof(int));
		memset(len_v1, 0, maxIte*sizeof(int));
		memset(sum_pre_v1, 0, maxIte*sizeof(int));
	} else{ 
		printf("LP not solved\n");
	}
}

//on the reduced problem
void solveIP(double time_limit, int ite, double gap){
	IloEnv env;
	IloModel model(env);
	IloIntVarArray x(env, n_var, 0, 1);
	IloRangeArray con(env);

	IloExpr obj(env);
	int i,j,k;
	for(i=0;i<n;i++){
		for(j=0;j<group[i].nItem;j++){
			obj+=group[i].profit[j]*x[sum_n_item[i]+j];
		}
	}
	model.add(IloMaximize(env, obj));
	obj.end();
	
	for(k=0;k<m;k++){
		IloExpr expr(env);
		for(i=0;i<n;i++){
			for(j=0;j<group[i].nItem;j++){
				expr+=group[i].consume[j][k]*x[sum_n_item[i]+j];
			}
		}
		con.add(expr<=Gb[k]);
		expr.end();
	}

	for(i=0;i<n;i++){
		IloExpr expr(env);
		for(j=0;j<group[i].nItem;j++){
			expr+=x[sum_n_item[i]+j];
		}
		con.add(expr==1);
		expr.end();
	}

	/*
	int *tempFix=new int[n_var];
	memset(tempFix, 0, n_var*sizeof(int));
	for(i=0;i<rp.len;i++){
		tempFix[rp.posi[i]]=1;
	}
	//force those not appeared in rp to 0
	for(i=0;i<n;i++){
		IloExpr expr(env);
		for(j=0;j<group[i].nItem;j++){
			if(tempFix[sum_n_item[i]+j]==0){
				expr+=x[sum_n_item[i]+j];
			}
		}
		con.add(expr==0);
		expr.end();
	}
	delete []tempFix;
	*/
	//add constrains from the last iteration
	if(ite>0){
		for(k=0;k<ite;k++){
			IloExpr expr(env);
			for(i=0;i<len_v1[k];i++){
				expr+=x[pre_v1[k][i]];
			}
			for(i=0;i<len_v0[k];i++){
				expr-=x[pre_v0[k][i]];
			}
			con.add(expr<=sum_pre_v1[k]-1);
			expr.end();
		}
	}
	if(0&&len_rdCost0>0&&len_rdCost1>0){
		IloExpr expr(env);
		for(i=0;i<len_rdCost1;i++){
			expr+=(1-x[rdCost_sort1[i].posi])*rdCost_sort1[i].value;
		}

		for(i=0;i<len_rdCost0;i++){
			expr+=x[rdCost_sort[i].posi]*rdCost_sort[i].value;
		}
		con.add(expr<=gap);
		expr.end();
	}
	if(len_rdCost0>0){
		printf("0 rdCostSort_len:%d. \n",len_rdCost0);
		k=1;
		int start=0, end;
		int cnt=0;
		//for(i=0;i<len_rdCost0;i++) printf("%lf, ", rdCost_sort[i].value);
		while(start<len_rdCost0){
			double sum=0;
			for(i=start;i<len_rdCost0-k;i++){
				sum=0;
				for(int ii=i;ii<i+k;ii++){
					sum+=rdCost_sort[ii].value;
				}
				if(sum<gap) break;
			}
			if(i>=len_rdCost0-k && sum<gap) break;
			//printf("cut %d, start=%d, end=%d\n",k-1, start, i+k);
			end=i+k;
			IloExpr expr1(env);
			for(i=start;i<end;i++){
				expr1+=x[rdCost_sort[i].posi];
			}
			con.add(expr1<=k-1);
			expr1.end();
			for(i=start;i<end;i++){
				pre_v0[ite][cnt++]=rdCost_sort[i].posi;
			}
			start=end;
			k++;
		}
		len_v0[ite]=cnt;
	}
	if(len_rdCost1>0){
		printf("1 rdCost_sort length:%d. \n",len_rdCost1);
		k=1;
		int start=0, end;
		int cnt=0;
		sum_pre_v1[ite]=0;
		//for(i=0;i<len_rdCost1;i++) printf("%lf, ", rdCost_sort1[i].value);
		while(start<len_rdCost1){
			double sum=0;
			for(i=start;i<len_rdCost1-k;i++){
				sum=0;
				for(int ii=i;ii<i+k;ii++){
					sum+=rdCost_sort1[ii].value;
				}
				if(sum<gap) break;
			}
			if(i>=len_rdCost1-k && sum<gap) break;
			//printf("cut %d, start=%d, end=%d\n",k-1, start, i+k);
			end=i+k;
			IloExpr expr1(env);
			for(i=start;i<end;i++){
				expr1+=x[rdCost_sort1[i].posi];
			}
			int cnt1=(end-start)-(k-1);
			con.add(expr1>=cnt1);
			expr1.end();
			for(i=start;i<end;i++){
				pre_v1[ite][cnt++]=rdCost_sort1[i].posi;
			}
			sum_pre_v1[ite]=sum_pre_v1[ite]+(end-start);
			start=end;
			k++;
		}
		len_v1[ite]=cnt;
	}
	model.add(con);
	IloCplex cplex(model);
	cplex.setParam(IloCplex::TiLim, time_limit);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::NodeFileInd,2);
	//cplex.setOut(env.getNullStream());
	if(LB>0){
		printf("\t add start point\n");
		IloNumVarArray start_var(env);
		IloNumArray start_val(env);
		int k=0;
		for(i=0;i<n;i++){
			for(j=0;j<group[i].nItem;j++){
				if(bestSolVec[i]==j){
					start_val.add(1);
				} else start_val.add(0);
				start_var.add(x[k]);
				k++;
			}
		}
		cplex.addMIPStart(start_var, start_val);
		start_var.end();
		start_val.end();
	}
	if(!cplex.solve()){
		std::cout<<"not solved IP."<<std::endl;
	} else {
		std::cout<<"solution status "<<cplex.getStatus()<<std::endl;
		std::cout<<"solution value "<<cplex.getObjValue()<<std::endl;
		
		if((double)cplex.getObjValue()>LB){
			LB=(double)cplex.getObjValue();	
			IloNumArray vals(env);
			cplex.getValues(vals, x);
			printf("new best solution vector is: ");
			for(i=0;i<n;i++){
				for(j=0;j<group[i].nItem;j++){
					if(cplex.getValue(x[sum_n_item[i]+j])>1-EPSILON){
					printf("%d ", j);
					bestSolVec[i]=j;
					}
				}
			}
			printf("\n\n");
		}
	}
	env.end();
}

int main(int argc, char *argv[]){
	if(argc<3){
		printf("useage:./solve file time_limit\n");
		return 0;
	}
	srand(11);
	read_inst(argv[1]);
	int time_limit=atoi(argv[2]);
	//solveIP(10);
	solveLP();
	cpu_time();
	int cnt=0;
	double ratio=0.5;
	int i=len_rdCost1-1;
	double maxRdc=0.0;
	maxIte=n*10;
	while(cpu_time()<time_limit && maxRdc<UB-LB-EPSILON){
		double t_remain=time_limit-cpu_time();
		if(t_remain<EPSILON) break;
		if(i<0) maxRdc= UB - LB;
		else maxRdc=rdCost_sort1[i].value-EPSILON;
		printf("iteration %d, cup_time:%lf, maxRdc %lf, LB =%lf\n", cnt, cpu_time(),maxRdc, LB);
		solveIP(t_remain,cnt, maxRdc);
		cnt++;
		printf("--------------\n\n");
		i--;
		//while(i>=0&&rdCost_sort1[i].value<maxRdc+0.1){i--;}
	}
	printf("best LB is %lf, cpu_time %lf\n", LB, cpu_time());
	deleteAll();
	return 0;
}
