#include <iostream>
#include <Rcpp.h>
#include <functional>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
using namespace Rcpp;
/* random number functions */
std::vector<double> rnorm_cpp2(int n, double mean, double stdev){
	std::normal_distribution<double> distribution(mean, stdev);
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::vector<double> out(n);
	for(int i=0; i < n; i++){
		out[i] = distribution(generator);
	}
	return out;
}
std::vector<int> rbin01_2(int n, double pr){
	std::uniform_real_distribution<double> distribution(0,1);
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::vector<int> p(n);
	for(int i=0; i < n; i++){
		double tmp = distribution(generator);
		if(tmp <= pr){
			p[i] = 1;
		} else {
			p[i] = 0;
		}
	}
	return p;
}
int rbin01_single(double pr){
	std::uniform_real_distribution<double> distribution(0,1);
	std::random_device rd;
	std::default_random_engine generator(rd());
	int p;
	double tmp = distribution(generator);
	if(tmp <= pr){
		p=1;
	}else{
		p=0;
	}
	return p;
}
int adjfunc2_cpp(std::vector<int> x, std::vector<int> xst){
	int n = x.size();
	std::vector<bool> smp(n);
	int ind1=0;
	int x1;
	for(int i=0; i < n; i++){
		if(xst[x[i]] > 0){
			smp[i] = TRUE;
			ind1++;
		}
	}
	if(ind1 == 0){
		int k = rand() % n;
		x1 = x[k];
	} else {
		std::vector<int> xvec(ind1);
		int ind2=0;
		for(int j=0; j < n; j++){
			if(xst[x[j]] > 0){
				xvec[ind2] = x[j];
				ind2++;
			}
		}
		int k = rand() % xvec.size();
		x1 = xvec[k];
	}
	return x1;
} // randomly select a neighbor with a non-zero amount

/* standard need-based open-ended network (osotua) condition */
std::vector<int> needbased_Trades(std::vector<int> st, std::vector<double> gr, std::vector<double> vc, Rcpp::List adjlist){
	int n = st.size();
	int min_st = 64;
	std::vector<bool> nlist(n);
	int yt = 0;
	for(int i=0; i < n; i++){
		st[i] += (st[i]*gr[i]) - (st[i]*vc[i]) + 0.5;
//		st[i] -= (st[i]*vc[i]) + 0.5;
		if((st[i] < min_st) && (st[i] > 1)){
			nlist[i] = TRUE;
			yt++;
		}
	}
	if(yt != 0){
	//	std::vector<int> ramt(n);
		std::vector<int> surpl(n);
		std::vector<int> iso_rec(yt);  // 'recipient' or receiver of request
		std::vector<int> iso_init(yt);  // 'initiates' as the requester
		std::vector<int> rlist(yt);  // request list
		// for re-arranging the queue to randomize risk burden per node
		std::vector<int> xord(yt);
		
		int ind = 0;
		for(int k=0; k < n; k++){
			if(nlist[k] == TRUE){
			//	ramt[k] = min_st - st[k];
				rlist[ind] = min_st - st[k];
				iso_rec[ind] = adjfunc2_cpp(adjlist[k], st);
				iso_init[ind] = k;
				xord[ind] = ind;
				ind++;
			} else {
				surpl[k] = st[k] - min_st;
			}
		}
		// randomize xord here
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(xord.begin(), xord.end(), g);
		for(int z=0; z < yt; z++){
			if(surpl[iso_rec[xord[z]]] > 0){
				if(surpl[iso_rec[xord[z]]] > rlist[xord[z]] ){
					st[iso_rec[xord[z]]] -= rlist[xord[z]];
					st[iso_init[xord[z]]] += rlist[xord[z]];
					surpl[iso_rec[xord[z]]] -= rlist[xord[z]];
				} else {
					st[iso_rec[xord[z]]] -= surpl[iso_rec[xord[z]]];
					st[iso_init[xord[z]]] += surpl[iso_rec[xord[z]]];
					surpl[iso_rec[xord[z]]] = 0;
				}	
			}
		}
	}
	return st;
}
std::vector<int> network_osotuaStandard(int dur, int n, int init_st, std::vector<double> glist, std::vector<double> vlist, Rcpp::List adjlist){	
	// initialization: dur(ation), n-players, st_vec(stock, n), g(rowth) list, v(olatility) list, adjlist (possible neighbors per n person)
	std::vector<int> st_vec(n);
	std::vector<int> ut(n);
	std::vector<int> fin(n);
	for(int q=0; q < n; q++){
		st_vec[q] = init_st;
		fin[q] = dur;
	}
	// beginning of the simuation loop
	for(int y=0; y < dur; y++){
		
		std::vector<double> gr_vec(n);
		std::vector<double> vc_vec(n);
		for(int i=0; i < n; i++){
			gr_vec[i] = glist[(y*n)+i];
			vc_vec[i] = vlist[(y*n)+i];
		}
		st_vec = needbased_Trades(st_vec, gr_vec, vc_vec, adjlist);
		for (int p = 0; p < n; p++) {
			if ((st_vec[p] < 64) && (st_vec[p] != 0)) {
				ut[p]++;
			}
			if (st_vec[p] > 600) {
				st_vec[p] = 600;
			}
		}
		for (int h = 0; h < n; h++) {
			if (ut[h] >= 2) {
				st_vec[h] = 0;
				fin[h] = (y + 1);
				ut[h] = 0;
			}
		}
		int meta_herd = 0;
		for(int u=0; u < n; u++){
			meta_herd += st_vec[u];
		}
		if(meta_herd == 0){
			break;
		}
	}
	// end of the simulation loop
	std::vector<int> out_vec(2*n);
	for(int j=0; j < n; j++){
		out_vec[j] = st_vec[j];
		out_vec[j+n] = fin[j];
	}
	return out_vec;
}

/* no exchange control/ defect (stingy) condition */
std::vector<int> needbased_Defect(std::vector<int> st, std::vector<double> gr, std::vector<double> vc, std::vector<int> dft, Rcpp::List adjlist){
	int n = st.size();
	int min_st = 64;
	std::vector<bool> nlist(n);
	int yt = 0;
	for(int i=0; i < n; i++){
		st[i] += (st[i]*gr[i]) - (st[i]*vc[i]) + 0.5;
//		st[i] -= (st[i]*vc[i]) + 0.5;
		if((st[i] < min_st) && (st[i] > 1)){
			nlist[i] = TRUE;
			yt++;
		}
	}
	
	if(yt != 0){
	//	std::vector<int> ramt(n);
		std::vector<int> surpl(n);
		std::vector<int> iso_rec(yt);  // 'recipient' or receiver of request
		std::vector<int> iso_init(yt);  // 'initiates' as the requester
		std::vector<int> rlist(yt);  // request list
		
		// for re-arranging the queue to randomize risk burden per node
		std::vector<int> xord(yt);
		
		int ind = 0;
		for(int k=0; k < n; k++){
			if(nlist[k] == TRUE){
			//	ramt[k] = min_st - st[k];
				rlist[ind] = min_st - st[k];
				iso_rec[ind] = adjfunc2_cpp(adjlist[k], st);
				iso_init[ind] = k;
				xord[ind] = ind;
				ind++;
			} else {
				surpl[k] = st[k] - min_st;
			}
		}
		
		// randomize xord here
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(xord.begin(), xord.end(), g);
		
		for(int z=0; z < yt; z++){
			if((surpl[iso_rec[xord[z]]] > 0) && (dft[iso_rec[xord[z]]] == 0)){
				if(surpl[iso_rec[xord[z]]] > rlist[xord[z]] ){
					st[iso_rec[xord[z]]] -= rlist[xord[z]];
					st[iso_init[xord[z]]] += rlist[xord[z]];
					surpl[iso_rec[xord[z]]] -= rlist[xord[z]];
				} else {
					st[iso_rec[xord[z]]] -= surpl[iso_rec[xord[z]]];
					st[iso_init[xord[z]]] += surpl[iso_rec[xord[z]]];
					surpl[iso_rec[xord[z]]] = 0;
				}
				
			}
		}
		
	}
	return st;
}
std::vector<int> network_osotuaDefect(int dur, int n, int init_st, std::vector<double> glist, std::vector<double> vlist, std::vector<int> dftlist, Rcpp::List adjlist){	
	// initialization: dur(ation), n-players, st_vec(stock, n), g(rowth) list, v(olatility) list, adjlist (possible neighbors per n person)
	std::vector<int> st_vec(n);
	std::vector<int> ut(n);
	std::vector<int> fin(n);
	for(int q=0; q < n; q++){
		st_vec[q] = init_st;
		fin[q] = dur;
	}
	// beginning of the simuation loop
	for(int y=0; y < dur; y++){
		
		std::vector<double> gr_vec(n);
		std::vector<double> vc_vec(n);
		std::vector<int> dft_vec(n);
		for(int i=0; i < n; i++){
			gr_vec[i] = glist[(y*n)+i];
			vc_vec[i] = vlist[(y*n)+i];
			dft_vec[i] = dftlist[(y*n)+i];
		}
		
		st_vec = needbased_Defect(st_vec, gr_vec, vc_vec, dft_vec, adjlist);
		
		for (int p = 0; p < n; p++) {
			if ((st_vec[p] < 64) && (st_vec[p] != 0)) {
				ut[p]++;
			}
			if (st_vec[p] > 600) {
				st_vec[p] = 600;
			}
		}
		for (int h = 0; h < n; h++) {
			if (ut[h] >= 2) {
				st_vec[h] = 0;
				fin[h] = (y + 1);
				ut[h] = 0;
			}
		}
		
		int meta_herd = 0;
		for(int u=0; u < n; u++){
			meta_herd += st_vec[u];
		}
		if(meta_herd == 0){
			break;
		}
	}
	// end of the simulation loop
	
	std::vector<int> out_vec(2*n);
	for(int j=0; j < n; j++){
		out_vec[j] = st_vec[j];
		out_vec[j+n] = fin[j];
	}
	return out_vec;
}

/*  standard debt-based open-ended network (esile) condition */
// NOTE: ** ONLY DYADS IN CURRENT FORM **
std::vector<int> debtbased_Trades(std::vector<int> st, std::vector<double> gr, std::vector<double> vc, Rcpp::List adjlist, std::vector<int> owes, std::vector<int> creditor, std::vector<double> pr_repay, std::vector<int> standing, std::vector<int> delay, int tolerated, int credit_limit){
	int n = st.size();    // *** NOTE: might need pairwise only, though ideally avoided -- USE LIST instead of vector for owes
	int min_st = 64;
	std::vector<bool> nlist(n);
	int yt = 0;
	for(int i=0; i < n; i++){
		st[i] += (st[i]*gr[i]) - (st[i]*vc[i]) + 0.5;
		if((st[i] < min_st) && (st[i] > 1)){
			nlist[i] = TRUE;
			yt++;
		}
	}
	
	// BEGINNING OF ACCOUNT-KEEPING SECTION ***
	// payback rule
	std::vector<int> rsurp(n);
	for(int i=0; i < n; i++){
		if(st[i] > min_st){
			rsurp[i] = st[i] - min_st;
		}
	}
	
/*	int ownum =0;
	for(int i=0; i < n; i++){
		if((owes[i] > 0) && (rsurp[i]>=owes[i])){
			ownum++;
		}
	}
	std::vector<int> repaym(ownum);	
	int prtmp = 1;
	std::vector<int> repaym(1); */
	
	
	for(int i=0; i < n; i++){
		if(owes[i] > 0){
			if(rsurp[i]>=owes[i]){  
				// prtmp = rand() % pr_repay[i];
				//prtmp = pr_repay[i];
				int repaym = rbin01_single(pr_repay[i]);
			//	repaym = rbin01_2(1,prtmp);
				if(repaym==1){
					st[creditor[i]] += owes[i];
					st[i] -= owes[i];
					owes[i] = 0;
					standing[i] = 1;
					delay[i] = 0;
				} else{
					delay[i]++; // didn't pay though able
				}
			} else {
				delay[i]++;  // didn't pay because unable
			}
		}	
	}
	
	// standing check
	for(int i=0; i < n; i++){
		if((delay[i] <= tolerated) && (owes[i] <= credit_limit)){
			standing[i] = 1;
		} else {
			standing[i] = 0;
		}
	}
	
	// END OF ACCOUNT-KEEPING SECTION ***
	
	// partner credit check rule && need-based request
	
	if(yt != 0){
		std::vector<int> surpl(n);
		std::vector<int> iso_rec(yt);  // 'recipient' or receiver of request
		std::vector<int> iso_init(yt);  // 'initiates' as the requester
		std::vector<int> rlist(yt);  // request list
		// for re-arranging the queue to randomize risk burden per node
		std::vector<int> xord(yt);
		
		int ind = 0;      // this is where the need-based requests are coordinated between init and recipient (same in DBT vs. NBT), order is randomized, and surplus is determined if no need exists for individuals
		for(int k=0; k < n; k++){
			if(nlist[k] == TRUE){
			//	ramt[k] = min_st - st[k];
				rlist[ind] = min_st - st[k];
				iso_rec[ind] = adjfunc2_cpp(adjlist[k], st);
				iso_init[ind] = k;
				xord[ind] = ind;
				ind++;
			} else {
				surpl[k] = st[k] - min_st;
			}
		}
		// randomize xord here
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(xord.begin(), xord.end(), g);
		
		// this is where surpluses/requests are fulfilled, conditional only on non-zero surplus (NBT) and possibly non-defector strategy (Defect); conditional added here: non-zero surplus+ (1) not exceeding tolerated delay AND (2) not exceeding (capping at?) credit size (i.e., good standing, standing == 1)
		for(int z=0; z < yt; z++){
			if(surpl[iso_rec[xord[z]]] > 0 && standing[iso_init[xord[z]]]==1){
				if(surpl[iso_rec[xord[z]]] > rlist[xord[z]]){
					if(rlist[xord[z]] > credit_limit){
						st[iso_rec[xord[z]]] -= credit_limit;
						st[iso_init[xord[z]]] += credit_limit;
						surpl[iso_rec[xord[z]]] -= credit_limit;
						owes[iso_init[xord[z]]] += credit_limit;
					} else {
						st[iso_rec[xord[z]]] -= rlist[xord[z]];
						st[iso_init[xord[z]]] += rlist[xord[z]];
						surpl[iso_rec[xord[z]]] -= rlist[xord[z]];
						owes[iso_init[xord[z]]] += rlist[xord[z]];
					}
				} else {
					if(surpl[iso_rec[xord[z]]] > credit_limit){
						st[iso_rec[xord[z]]] -= credit_limit;
						st[iso_init[xord[z]]] += credit_limit;
						surpl[iso_rec[xord[z]]] = 0;
						owes[iso_init[xord[z]]] += credit_limit;
					} else {
						st[iso_rec[xord[z]]] -= surpl[iso_rec[xord[z]]];
						st[iso_init[xord[z]]] += surpl[iso_rec[xord[z]]];
						surpl[iso_rec[xord[z]]] = 0;
						owes[iso_init[xord[z]]] += surpl[iso_rec[xord[z]]];
					}
				}	
			}
		}
		
		
	}
	
	st.insert(st.end(), delay.begin(), delay.end());
	st.insert(st.end(), owes.begin(), owes.end());
	
/*	std::vector<int> testv;
	std::vector<int>::iterator mItr(st.begin()+st.size()/2);
	for(auto it=st.begin(); it != st.end(); it++){
		if(std::distance(it,mItr) > 0){
			testv.push_back(*it);
		}
	} 
	std::vector<int> testv1(n);
	std::vector<int> testv2(n);
	std::vector<int> testv3(n);
	for(int i=0; i < n; i++){
		testv1[i] = st[i];
		testv2[i] = st[i+n];
		testv3[i] = st[i+(2*n)];
	} */
	
	
//	std::vector<int> testv(st.begin(), st.begin() + st.end() / 2);
	// append (and account for) st and delay vectors
	
	return st;
	//return st;
//	return st;
}
std::vector<int> network_Esile(int dur, int n, int init_st, std::vector<double> glist, std::vector<double> vlist, Rcpp::List adjlist, std::vector<int> creditor, std::vector<double> pr_repay, int tolerated, int credit_limit){	
	// initialization: dur(ation), n-players, st_vec(stock, n), g(rowth) list, v(olatility) list, adjlist (possible neighbors per n person)
	std::vector<int> st_vec(n);
	std::vector<int> ut(n);
	std::vector<int> fin(n);
	for(int q=0; q < n; q++){
		st_vec[q] = init_st;
		fin[q] = dur;
	}
	// npackvec, owe_vec, stand_vec, delay_vec
	std::vector<int> npack_vec(3*n);
	std::vector<int> owe_vec(n);
	std::vector<int> stand_vec(n);
	std::vector<int> delay_vec(n);
	for(int i=0; i < n; i++){
		owe_vec[i] = 0;
		stand_vec[i] =1;
		delay_vec[i] = 0;
	}
	// beginning of the simuation loop
	for(int y=0; y < dur; y++){
		std::vector<double> gr_vec(n);
		std::vector<double> vc_vec(n);
//		std::vector<int> dft_vec(n);
		for(int i=0; i < n; i++){
			gr_vec[i] = glist[(y*n)+i];
			vc_vec[i] = vlist[(y*n)+i];
//			dft_vec[i] = dftlist[(y*n)+i];
		}
		npack_vec = debtbased_Trades(st_vec, gr_vec, vc_vec, adjlist, owe_vec, creditor, pr_repay, stand_vec, delay_vec, tolerated, credit_limit);
//		std::vector<int> st_vtmp(n);
//		std::vector<int> del_vtmp(n);
//		std::vector<int> owe_vtmp(n); 		
		for(int i=0; i < n; i++){
			st_vec[i] = npack_vec[i];
			delay_vec[i] = npack_vec[i+n];
			owe_vec[i] = npack_vec[i+(2*n)];
		}		
/*		std::vector<int>::iterator mItr(npack_vec.begin()+npack_vec.size()/2);				
		for(auto it=npack_vec.begin(); it != npack_vec.end(); it++){
			if(std::distance(it,mItr) > 0){
				st_vtmp.push_back(*it);
			} else {
				del_vtmp.push_back(*it);
			}
		}  */
//		st_vec = st_vtmp;
//		delay_vec = del_vtmp;
		
// debtbased_Trades(std::vector<int> st, std::vector<double> gr, std::vector<double> vc, Rcpp::List adjlist, std::vector<int> owes, std::vector<int> creditor, std::vector<double> pr_repay, std::vector<int> standing, std::vector<int> delay, int tolerated, int credit_limit)		
		
//		st_vec = needbased_Defect(st_vec, gr_vec, vc_vec, dft_vec, adjlist);
		
		for (int p = 0; p < n; p++) {
			if ((st_vec[p] < 64) && (st_vec[p] != 0)) {
				ut[p]++;
			}
			if (st_vec[p] > 600) {
				st_vec[p] = 600;
			}
		}
		for (int h = 0; h < n; h++) {
			if (ut[h] >= 2) {
				st_vec[h] = 0;
				fin[h] = (y + 1);
				ut[h] = 0;
			}
		}
		
		int meta_herd = 0;
		for(int u=0; u < n; u++){
			meta_herd += st_vec[u];
		}
		if(meta_herd == 0){
			break;
		}
	}
	// end of the simulation loop
	
	std::vector<int> out_vec(2*n);
	for(int j=0; j < n; j++){
		out_vec[j] = st_vec[j];
		out_vec[j+n] = fin[j];
	}
	return out_vec;
}

/* need-based defect (greedy) condition */
/* feign open-ended network osotua condition ('t1' = single time step) */
std::vector<int> needbased_Greedy(std::vector<int> st, std::vector<double> gr, std::vector<double> vc, std::vector<int> dft, Rcpp::List adjlist){
	int n = st.size();
	int min_st = 64;
	std::vector<bool> nlist(n);
	int yt = 0;
	for(int i=0; i < n; i++){
		st[i] += (st[i]*gr[i]) - (st[i]*vc[i]) + 0.5;
//		st[i] -= (st[i]*vc[i]) + 0.5;
		if((st[i] < min_st) && (st[i] > 1)){
			nlist[i] = TRUE;
			yt++;
		}
		if((dft[i] == 1) && (nlist[i] == FALSE)){
			nlist[i] = TRUE;
			yt++;
		}
	}
	
	if(yt != 0){
	//	std::vector<int> ramt(n);
		std::vector<int> surpl(n);
		std::vector<int> iso_rec(yt);  // 'recipient' or receiver of request
		std::vector<int> iso_init(yt);  // 'initiates' as the requester
		std::vector<int> rlist(yt);  // request list
		
		// for re-arranging the queue to randomize risk burden per node
		std::vector<int> xord(yt);
		
		int ind = 0;
		for(int k=0; k < n; k++){
			if(nlist[k] == TRUE){
			//	ramt[k] = min_st - st[k];
				if((dft[k] == 1) && (st[k] > min_st)){
					std::default_random_engine generator;
					std::normal_distribution<double> distribution(11.73, 7.66); // taken from bootstrapped sample
					rlist[ind] = distribution(generator);
			
			//		rlist[ind] = 12;
				} else {
					rlist[ind] = min_st - st[k];
				}
				
		//		rlist[ind] = min_st - st[k];
				iso_rec[ind] = adjfunc2_cpp(adjlist[k], st);
				iso_init[ind] = k;
				xord[ind] = ind;
				ind++;
			} else {
				surpl[k] = st[k] - min_st;
			}
		}
		
		// randomize xord here
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(xord.begin(), xord.end(), g);
		
		for(int z=0; z < yt; z++){
			if((surpl[iso_rec[xord[z]]] > 0) && (dft[iso_rec[xord[z]]] == 0)){
				if(surpl[iso_rec[xord[z]]] > rlist[xord[z]] ){
					st[iso_rec[xord[z]]] -= rlist[xord[z]];
					st[iso_init[xord[z]]] += rlist[xord[z]];
					surpl[iso_rec[xord[z]]] -= rlist[xord[z]];
				} else {
					st[iso_rec[xord[z]]] -= surpl[iso_rec[xord[z]]];
					st[iso_init[xord[z]]] += surpl[iso_rec[xord[z]]];
					surpl[iso_rec[xord[z]]] = 0;
				}
				
			}
		}
		
	}
	return st;
}
std::vector<int> network_osotuaGreedy(int dur, int n, int init_st, std::vector<double> glist, std::vector<double> vlist, std::vector<int> dftlist, Rcpp::List adjlist){	
	// initialization: dur(ation), n-players, st_vec(stock, n), g(rowth) list, v(olatility) list, adjlist (possible neighbors per n person)
	std::vector<int> st_vec(n);
	std::vector<int> ut(n);
	std::vector<int> fin(n);
	for(int q=0; q < n; q++){
		st_vec[q] = init_st;
		fin[q] = dur;
	}
	// beginning of the simuation loop
	for(int y=0; y < dur; y++){
		
		std::vector<double> gr_vec(n);
		std::vector<double> vc_vec(n);
		std::vector<int> dft_vec(n);
		for(int i=0; i < n; i++){
			gr_vec[i] = glist[(y*n)+i];
			vc_vec[i] = vlist[(y*n)+i];
			dft_vec[i] = dftlist[(y*n)+i];
		}
		
		st_vec = needbased_Greedy(st_vec, gr_vec, vc_vec, dft_vec, adjlist);
		
		for (int p = 0; p < n; p++) {
			if ((st_vec[p] < 64) && (st_vec[p] != 0)) {
				ut[p]++;
			}
			if (st_vec[p] > 600) {
				st_vec[p] = 600;
			}
		}
		for (int h = 0; h < n; h++) {
			if (ut[h] >= 2) {
				st_vec[h] = 0;
				fin[h] = (y + 1);
				ut[h] = 0;
			}
		}
		
		int meta_herd = 0;
		for(int u=0; u < n; u++){
			meta_herd += st_vec[u];
		}
		if(meta_herd == 0){
			break;
		}
	}
	// end of the simulation loop
	
	std::vector<int> out_vec(2*n);
	for(int j=0; j < n; j++){
		out_vec[j] = st_vec[j];
		out_vec[j+n] = fin[j];
	}
	return out_vec;
}

// [[Rcpp::export]]
int EsileStandard(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, std::vector<int> creditor, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, std::vector<double> pr_repay, int tolerated, int credit_limit, bool newfile){
	std::ofstream outfile("EsileStandard.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("EsileStandard.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		vlint = rbin01_2(xlen, vtrate);
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);
		
		outvecf = network_Esile(dur, n, init_st, glist, vlist, adjlist, creditor, pr_repay, tolerated, credit_limit);
		
// int dur, int n, int init_st, std::vector<double> glist, std::vector<double> vlist, Rcpp::List adjlist, std::vector<int> creditor, std::vector<double> pr_repay, int tolerated, int credit_limit		
//		outvecf = network_osotuaStandard(dur, n, init_st, glist, vlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}

// [[Rcpp::export]]
int OsotuaStandard(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, bool newfile){
	std::ofstream outfile("OsotuaStandard.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("OsotuaStandard.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		vlint = rbin01_2(xlen, vtrate);
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);
		outvecf = network_osotuaStandard(dur, n, init_st, glist, vlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}

// [[Rcpp::export]]
int OsotuaStingy(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, std::vector<int> wdefect, double prdefect, bool newfile){
	
//const char filename = 'OsotuaDefect.csv';
//	if(newfile==TRUE){		
	// cin >> charname;
//	}
	
	std::ofstream outfile("OsotuaStingy.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("OsotuaStingy.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
//		std::vector<int> defectlist(xlen);
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		
		// ** setting vlist **
		vlint = rbin01_2(xlen, vtrate);     
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		
		// ** setting grlist **
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);      
		
		// ** setting defectlist **
		std::vector<int> dftmp(n);
		int ind = 0;
		for(int z=0; z < n; z++){
			if(z==wdefect[ind]){
				dftmp[z] = 1;
				ind++;
			} else {
				dftmp[z] = 0;
			}
		}
		std::vector<int> defect1(xlen);
		for(int z=0; z < (xlen/n); z++){
			for(int dk=0; dk < n; dk++){
				defect1[(z*n)+dk] = dftmp[dk];
			}
		}
		int dfsum = 0;
		for(int z=0; z < xlen; z++){
			dfsum += defect1[z];
		}
		std::vector<int> defect2(dfsum);
		defect2 = rbin01_2(dfsum, prdefect);
		std::vector<int> defectlist(xlen);
		ind=0;
		for(int z=0; z < xlen; z++){
			if(defect1[z]==1){
				defectlist[z] = defect2[ind];
				ind++;
			}
		}
		outvecf = network_osotuaDefect(dur, n, init_st, glist, vlist, defectlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}

// [[Rcpp::export]]
int OsotuaDefect2(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, std::vector<int> wdefect, double prdefect, bool newfile){
	
//const char filename = 'OsotuaDefect.csv';
//	if(newfile==TRUE){		
	// cin >> charname;
//	}
	
	std::ofstream outfile("OsotuaDefect2.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("OsotuaDefect2.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
//		std::vector<int> defectlist(xlen);
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		
		// ** setting vlist **
		vlint = rbin01_2(xlen, vtrate);     
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		
		// ** setting grlist **
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);      
		
		// ** setting defectlist **
		std::vector<int> dftmp(n);
		int ind = 0;
		for(int z=0; z < n; z++){
			if(z==wdefect[ind]){
				dftmp[z] = 1;
				ind++;
			} else {
				dftmp[z] = 0;
			}
		}
		std::vector<int> defect1(xlen);
		for(int z=0; z < (xlen/n); z++){
			for(int dk=0; dk < n; dk++){
				defect1[(z*n)+dk] = dftmp[dk];
			}
		}
		int dfsum = 0;
		for(int z=0; z < xlen; z++){
			dfsum += defect1[z];
		}
		std::vector<int> defect2(dfsum);
		defect2 = rbin01_2(dfsum, prdefect);
		std::vector<int> defectlist(xlen);
		ind=0;
		for(int z=0; z < xlen; z++){
			if(defect1[z]==1){
				defectlist[z] = defect2[ind];
				ind++;
			}
		}
		outvecf = network_osotuaDefect(dur, n, init_st, glist, vlist, defectlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}

// [[Rcpp::export]]
int OsotuaGreedy(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, std::vector<int> wdefect, double prdefect, bool newfile){
	
//const char filename = 'OsotuaDefect.csv';
//	if(newfile==TRUE){		
	// cin >> charname;
//	}
	
	std::ofstream outfile("OsotuaGreedy.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("OsotuaGreedy.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
//		std::vector<int> defectlist(xlen);
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		
		// ** setting vlist **
		vlint = rbin01_2(xlen, vtrate);     
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		
		// ** setting grlist **
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);      
		
		// ** setting defectlist **
		std::vector<int> dftmp(n);
		int ind = 0;
		for(int z=0; z < n; z++){
			if(z==wdefect[ind]){
				dftmp[z] = 1;
				ind++;
			} else {
				dftmp[z] = 0;
			}
		}
		std::vector<int> defect1(xlen);
		for(int z=0; z < (xlen/n); z++){
			for(int dk=0; dk < n; dk++){
				defect1[(z*n)+dk] = dftmp[dk];
			}
		}
		int dfsum = 0;
		for(int z=0; z < xlen; z++){
			dfsum += defect1[z];
		}
		std::vector<int> defect2(dfsum);
		defect2 = rbin01_2(dfsum, prdefect);
		std::vector<int> defectlist(xlen);
		ind=0;
		for(int z=0; z < xlen; z++){
			if(defect1[z]==1){
				defectlist[z] = defect2[ind];
				ind++;
			}
		}
		outvecf = network_osotuaGreedy(dur, n, init_st, glist, vlist, defectlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}

// [[Rcpp::export]]
int OsotuaGreedy2(int tnum, int dur, int n, int init_st, Rcpp::List adjlist, double vtrate, double vt_mean, double vt_sd, double gr_mean, double gr_sd, std::vector<int> wdefect, double prdefect, bool newfile){
	
//const char filename = 'OsotuaDefect.csv';
//	if(newfile==TRUE){		
	// cin >> charname;
//	}
	
	std::ofstream outfile("OsotuaGreedy2.csv", std::ios_base::app);	
	if(newfile == TRUE){
		std::ofstream outfile("OsotuaGreedy2.csv", std::ofstream::out);
	}
	for(int t=0; t < tnum; t++){	
		int xlen = dur*n;
		std::vector<double> vlist(xlen);   // disaster calendar
		std::vector<int> vlint(xlen);
		std::vector<double> glist(xlen);   // growth schedule
//		std::vector<int> defectlist(xlen);
		
//		std::vector<double> outvecf;
		std::vector<int> outvecf;   // results ('out') vector
		
		// tmp01 = rbin01_2(flen, vtrate);
		
		// ** setting vlist **
		vlint = rbin01_2(xlen, vtrate);     
		int sumv = 0;
		for(int z=0; z < xlen; z++){
			sumv += vlint[z];
		}
		std::vector<double> tmplist(sumv);
		tmplist = rnorm_cpp2(sumv, vt_mean, vt_sd);
		int indc = 0;
		for(int z=0; z < xlen; z++){
			if(vlint[z]==1){
				vlist[z] = tmplist[indc];
				indc++;
			}
		}		
//		vlist = covar_gen(xlen, eco_covar, vtrate, vt_mean, vt_sd);
		
		// ** setting grlist **
		glist = rnorm_cpp2(xlen, gr_mean, gr_sd);      
		
		// ** setting defectlist **
		std::vector<int> dftmp(n);
		int ind = 0;
		for(int z=0; z < n; z++){
			if(z==wdefect[ind]){
				dftmp[z] = 1;
				ind++;
			} else {
				dftmp[z] = 0;
			}
		}
		std::vector<int> defect1(xlen);
		for(int z=0; z < (xlen/n); z++){
			for(int dk=0; dk < n; dk++){
				defect1[(z*n)+dk] = dftmp[dk];
			}
		}
		int dfsum = 0;
		for(int z=0; z < xlen; z++){
			dfsum += defect1[z];
		}
		std::vector<int> defect2(dfsum);
		defect2 = rbin01_2(dfsum, prdefect);
		std::vector<int> defectlist(xlen);
		ind=0;
		for(int z=0; z < xlen; z++){
			if(defect1[z]==1){
				defectlist[z] = defect2[ind];
				ind++;
			}
		}
		outvecf = network_osotuaGreedy(dur, n, init_st, glist, vlist, defectlist, adjlist);
		
//		outvecf.push_back(eco_covar*100);		
//		outvecf[outvecf.size()] = eco_covar;
		
		std::stringstream vecst2;
		for(int m=0; m < outvecf.size(); m++){
			vecst2 << outvecf[m];
			if(m != (outvecf.size()-1)){
				vecst2 << ",";
			}
		}
		std::string vecst3 = vecst2.str();
		outfile << vecst3 << std::endl;
    }
	outfile.close();
	return 0;	
}





/*
std::vector<int> gadget (int dur, int n, std::vector<int> wdefect, double prdefect){
	int xlen =dur*n;
	std::vector<int> dftmp(n);
	int ind = 0;
	for(int z=0; z < n; z++){
		if(z==wdefect[ind]){
			dftmp[z] = 1;
			ind++;
		} else {
			dftmp[z] = 0;
		}
	}
	std::vector<int> defect1(xlen);
	for(int z=0; z < (xlen/n); z++){
		for(int dk=0; dk < n; dk++){
			defect1[(z*n)+dk] = dftmp[dk];
		}
	}
	int dfsum = 0;
	for(int z=0; z < xlen; z++){
		dfsum += defect1[z];
	}
	std::vector<int> defect2(dfsum);
	defect2 = rbin01_2(dfsum, prdefect);
	std::vector<int> defectlist(xlen);
	ind=0;
	for(int z=0; z < xlen; z++){
		if(defect1[z]==1){
			defectlist[z] = defect2[ind];
			ind++;
		}
	}
	return defectlist;
	
}
*/

