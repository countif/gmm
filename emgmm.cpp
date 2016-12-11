#include<iostream>
#include<stdio.h>
#include<fstream>
#define K 3
#define M 100
#define D 4
using namespace std;
double Gamma[M][K];
double Data[M][D];
double Phi[K];
int ylabel[M];
double Mean[K][D];
double Variance[K][D];
int iteration = 10;

void Init(){
	memset(Gamma,0,sizeof(Gamma));
	memset(Data,0,sizeof(Data));
	memset(Phi,0,sizeof(Phi));
}
void GetData(string path){
	ifstream in(path.c_str(),ios::in);
	for (int i = 0;i<M;++i){
		for (int d = 0;d<D;++d){
			in>>Data[i][d];
		}
		in >>ylabel[i];
	}
	in.close();
}
double Gassuian(double feature[D],int k){
	double innerproduct = 1.0;
	double tepsum = 0;
	for (int d = 0;d<D;++d){
		tepsum += (feature[d] - Mean[k][d])*(feature[d] - Mean[k][d])/Variance[k][d];
		innerproduct *= sqrt(2*3.141592*Variance[k][d]);
	}
	tepsum = exp(-0.5*tepsum);
	double prefix = 1.0/innerproduct;
	return prefix * tepsum;

}
double getrand(){
	double res = rand();
	res /= RAND_MAX;
	return res;
}
void RandomInit(){
	for (int i = 0;i<M;++i){
		double tepsum = 0;
		for (int k = 0;k<K;++k){
			Gamma[i][k] = getrand();
			tepsum += Gamma[i][k];
		}
		for (int k = 0;k<K;++k)
			Gamma[i][k] /= tepsum;
	}


	for (int i = 0;i<K;++i){
		for (int k = 0;k<D;++k){
			Mean[i][k] = getrand();
			Variance[i][k] = getrand();
		}
	}

	double phisum = 0;
	for (int k = 0;k<K;++k){
		Phi[k] = getrand();
		phisum += Phi[k];
	}
	for (int k = 0;k<K;++k){
		Phi[k] /= phisum;
	}


	

}
void estamate(){


	for (int i = 0;i<M;++i){
		double tepsum = 0;

		for (int k = 0 ;k<K;++k){
			Gamma[i][k] = Phi[k] * Gassuian(Data[i],k);
			tepsum += Gamma[i][k];
		}
		for (int k = 0;k<K;++k){
			Gamma[i][k] /= tepsum;
		}
	}

}
void mstep(){
	double tepK[K];
	for (int k = 0;k<K;++k){
		for (int i  = 0;i<M;++i){
			tepK[k] += Gamma[i][k] ;
		}

	}

	double total = 0;
	for (int i = 0;i<K;++i){
		total += tepK[i] ;
	}


	//update mean
	for (int k = 0;k<K;++k){
		double MeanTepRes[D];
		memset(MeanTepRes,0,sizeof(MeanTepRes));
		for (int i = 0;i<M;++i){
			for (int d = 0;d<D;++d){
				MeanTepRes[d] += Gamma[i][k] * Data[i][d];
			}
		}

		for (int d = 0;d<D;++d){
			Mean[k][d] = MeanTepRes[d] / tepK[k];
		}

	}

	//update Variance
	for (int k = 0;k<K;++k){
			double VarianceTepRes[D];
			memset(VarianceTepRes,0,sizeof(VarianceTepRes));

			for (int i = 0;i<M;++i){
				double tepsum = 0;
				for (int d = 0;d<D;++d){
					VarianceTepRes[d] += Gamma[i][k] * (Data[i][d] - Mean[k][d]) * (Data[i][d] - Mean[k][d]);
				}
				
			}
			for (int d = 0 ;d <D;++d)
				Variance[k][d] = VarianceTepRes[d] / tepK[k];
		
	}

	// update phi
	for (int k = 0;k<K;++k){
		Phi[k] /= total;
	}
	

}

void train(){
	for (int iter = 0 ;iter < iteration;++iter){
		estamate();
		mstep();
	}
}
int main(){
	RandomInit();
	GetData("test.txt");
	train();
	return 0;
}