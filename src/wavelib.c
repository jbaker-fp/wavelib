/*
  Copyright (c) 2014, Rafat Hussain
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wavelib.h"
#include "wtmath.h"

wave_object wave_init(const char* wname) {
	wave_object obj = NULL;
	int retval;
	retval = 0;

	if (wname != NULL) {
		retval = filtlength(wname);
		//obj->filtlength = retval;
		//strcopy(obj->wname, wname);
	}

	obj = (wave_object)malloc(sizeof(struct wave_set) + sizeof(double)* 4 * retval);

	obj->filtlength = retval;
	obj->lpd_len = obj->hpd_len = obj->lpr_len = obj->hpr_len = obj->filtlength;
	strcpy(obj->wname, wname);
	if (wname != NULL) {
		filtcoef(wname,obj->params,obj->params+retval,obj->params+2*retval,obj->params+3*retval);
	}
	obj->lpd = &obj->params[0];
	obj->hpd = &obj->params[retval];
	obj->lpr = &obj->params[2 * retval];
	obj->hpr = &obj->params[3 * retval];
	return obj;
}

wpt_object wpt_init(wave_object wave, int siglength, int J) {
	int size, i, MaxIter, temp, nodes,elength,p2,N,lp;
	wpt_object obj = NULL;

	size = wave->filtlength;

	if (J > 100) {
		printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
		exit(-1);
	}


	MaxIter = wmaxiter(siglength, size);
	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
		exit(-1);
	}
	temp = 1;
	nodes = 0;
	for (i = 0; i < J; ++i) {
		temp *= 2;
		nodes += temp;
	}

	i = J;
	p2 = 2;
	N = siglength;
	lp = size;
	elength = 0;
	while (i > 0) {
		N = N + lp - 2;
		N = (int)ceil((double)N / 2.0);
		elength = p2 * N;
		i--;
		p2 *= 2;
	}
	//printf("elength %d", elength);

	obj = (wpt_object)malloc(sizeof(struct wpt_set) + sizeof(double)* (elength + 4 * nodes + 2 * J + 6));
	obj->outlength = siglength + 2 * (J + 1) * (size + 1);
	strcpy(obj->ext, "sym");
	strcpy(obj->entropy, "shannon");
	obj->eparam = 0.0;

	obj->wave = wave;
	obj->siglength = siglength;
	obj->J = J;
	obj->MaxIter = MaxIter;

	if (siglength % 2 == 0) {
		obj->even = 1;
	}
	else {
		obj->even = 0;
	}

	obj->cobj = NULL;
	obj->nodes = nodes;

	obj->lenlength = J + 2;
	obj->output = &obj->params[0];
	obj->costvalues = &obj->params[elength];
	obj->basisvector = &obj->params[elength + nodes + 1];
	obj->nodeindex = (int*)&obj->params[elength + 2*nodes + 2];
	obj->numnodeslevel = (int*)&obj->params[elength + 4 * nodes + 4];
	obj->coeflength = (int*)&obj->params[elength + 4 * nodes + J + 5];

	for (i = 0; i < elength + 4 * nodes + 2 * J + 6; ++i) {
		obj->params[i] = 0.0;
	}

	//wave_summary(obj->wave);

	return obj;
}

static void dwpt_per(wpt_object wt, double *inp, int N, double *cA, int len_cA, double *cD) {
	int l, l2, isodd, i, t, len_avg;

	len_avg = wt->wave->lpd_len;
	l2 = len_avg / 2;
	isodd = N % 2;

	for (i = 0; i < len_cA; ++i) {
		t = 2 * i + l2;
		cA[i] = 0.0;
		cD[i] = 0.0;
		for (l = 0; l < len_avg; ++l) {
			if ((t - l) >= l2 && (t - l) < N) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < l2 && (t - l) >= 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < 0 && isodd == 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l + N];
				cD[i] += wt->wave->hpd[l] * inp[t - l + N];
			}
			else if ((t - l) < 0 && isodd == 1) {
				if ((t - l) != -1) {
					cA[i] += wt->wave->lpd[l] * inp[t - l + N + 1];
					cD[i] += wt->wave->hpd[l] * inp[t - l + N + 1];
				}
				else {
					cA[i] += wt->wave->lpd[l] * inp[N - 1];
					cD[i] += wt->wave->hpd[l] * inp[N - 1];
				}
			}
			else if ((t - l) >= N && isodd == 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l - N];
				cD[i] += wt->wave->hpd[l] * inp[t - l - N];
			}
			else if ((t - l) >= N && isodd == 1) {
				if (t - l != N) {
					cA[i] += wt->wave->lpd[l] * inp[t - l - (N + 1)];
					cD[i] += wt->wave->hpd[l] * inp[t - l - (N + 1)];
				}
				else {
					cA[i] += wt->wave->lpd[l] * inp[N - 1];
					cD[i] += wt->wave->hpd[l] * inp[N - 1];
				}
			}

		}
	}



}

static void dwpt_sym(wpt_object wt, double *inp, int N, double *cA, int len_cA, double *cD) {
	int i, l, t, len_avg;

	len_avg = wt->wave->lpd_len;
	
	for (i = 0; i < len_cA; ++i) {
		t = i*2 + 1;
		cA[i] = 0.0;
		cD[i] = 0.0;

		for (l = 0; l <= MIN(t-N, len_avg - 1); l++) {
			cA[i] += wt->wave->lpd[l] * inp[2 * N - t - 1 + l];
			cD[i] += wt->wave->hpd[l] * inp[2 * N - t - 1 + l];
		}
		for (;l <= MIN(t, len_avg - 1); l++) {
			cA[i] += wt->wave->lpd[l] * inp[t - l];
			cD[i] += wt->wave->hpd[l] * inp[t - l];
		}
		for (; l < len_avg; l++) {
			cA[i] += wt->wave->lpd[l] * inp[-t + l - 1];
			cD[i] += wt->wave->hpd[l] * inp[-t + l - 1];
		}
	}
}
static int ipow2(int n) {
	int p,i;
	p = 1;

	for (i = 0; i < n; ++i) {
		p *= 2;
	}

	return p;
}

void dwpt(wpt_object wt, const double *inp) {
	int i, J, temp_len, iter, N, lp, p2, k, N2, Np;
	int temp, elength, temp2,size,nodes,llb,n1,j;
	double eparam,v1,v2;
	int len_cA, t, t2, it1,it2;
	double *orig,*tree;
	int *nodelength;

	temp_len = wt->siglength;
	J = wt->J;
	wt->length[J + 1] = temp_len;
	wt->outlength = 0;
	temp = 1;
	elength = 0;
	size = wt->wave->filtlength;
	nodes = wt->nodes;
	n1 = nodes + 1;
	for (i = 0; i < J; ++i) {
		temp *= 2;
		temp2 = (size - 2) * (temp - 1);
		elength += temp2;
	}
	eparam = wt->eparam;
	orig = (double*)malloc(sizeof(double)* temp_len);
	tree = (double*)malloc(sizeof(double)* (temp_len * (J + 1) + elength));
	nodelength = (int*)malloc(sizeof(int)* nodes);

	for (i = 0; i < wt->siglength; ++i) {
		orig[i] = inp[i];
	}

	for (i = 0; i < temp_len * (J + 1) + elength; ++i) {
		tree[i] = 0.0;
	}

	for (i = 0; i < nodes + 1; ++i) {
		wt->basisvector[i] = 0.0;
		wt->costvalues[i] = 0.0;
	}

	N = temp_len;
	lp = wt->wave->lpd_len;
	p2 = 1;

	//set eparam value here
	it2 = 1;
	if (!strcmp(wt->ext, "per")) {
		i = J;
		p2 = 2;
		while (i > 0) {
			N = (int)ceil((double)N / 2.0);
			wt->length[i] = N;
			wt->outlength += p2 * (wt->length[i]);
			i--;
			p2 *= 2;
		}
		wt->length[0] = wt->length[1];

		N2 = N = wt->outlength;
		p2 = 1;
		for (iter = 0; iter < J; ++iter) {
			len_cA = wt->length[J - iter];
			N2 -= 2 * p2 * len_cA;
			N = N2;
			for (k = 0; k < p2; ++k) {
				if (iter == 0) {
					dwpt_per(wt, orig, temp_len, tree + N, len_cA, tree + N + len_cA);
				}
				else {
					dwpt_per(wt, tree + Np + k * temp_len, temp_len, tree + N, len_cA, tree + N + len_cA);
				}
				wt->costvalues[it2] = costfunc(tree + N, len_cA, wt->entropy, eparam);
				it2++;
				wt->costvalues[it2] = costfunc(tree + N +len_cA, len_cA, wt->entropy, eparam);
				it2++;
				N += 2 * len_cA;
			}

			temp_len = wt->length[J - iter];
			p2 = 2 * p2;
			Np = N2;
		}
	}
	else if (!strcmp(wt->ext, "sym")) {
		//printf("\n YES %s \n", wt->ext);
		i = J;
		p2 = 2;
		while (i > 0) {
			N = N + lp - 2;
			N = (int)ceil((double)N / 2.0);
			wt->length[i] = N;
			wt->outlength += p2 * (wt->length[i]);
			i--;
			p2 *= 2;
		}
		wt->length[0] = wt->length[1];

		N2 = N = wt->outlength;
		p2 = 1;

		for (iter = 0; iter < J; ++iter) {
			len_cA = wt->length[J - iter];
			N2 -= 2 * p2 * len_cA;
			N = N2;
			for (k = 0; k < p2; ++k) {
				if (iter == 0) {
					dwpt_sym(wt, orig, temp_len, tree + N, len_cA, tree + N + len_cA);
				}
				else {
					dwpt_sym(wt, tree + Np + k * temp_len, temp_len, tree + N, len_cA, tree + N + len_cA);
				}
				N += 2 * len_cA;
			}

			temp_len = wt->length[J - iter];
			p2 = 2 * p2;
			Np = N2;
		}

	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}

	J = wt->J;
	t2 = wt->outlength - 2 * wt->length[J];
	p2 = 2;
	it1 = 0;
	for (i = 0; i < J; ++i) {
		t = t2;
		for (k = 0; k < p2; ++k) {
			nodelength[it1] = t;
			it1++;
			t += wt->length[J - i];
		}
		p2 *= 2;
		t2 = t2 - p2 * wt->length[J - i - 1];
	}


	J = wt->J;
	llb = 1;
	for (i = 0; i < J; ++i) {
		llb *= 2;
	}

	for (k = n1-llb; k < n1; k++) {
		wt->basisvector[k] = 1; 
	}

	N2 = 0;
	it1 = n1;
	it2 = 0;
	wt->nodes = 0;
	wt->numnodeslevel[0] = 0;

	// We only ever use the nodes from the last decomposition level now
	// Removed what I know I can remove Otherwise I left it
	// This could be simplified since we only need to loop 2^J times, but I have it working
	// and I really do not understand what is going on very well here
	// it is usually just checking a conditional so may not be worth changing
	if (wt->basisvector[0] == 1) {
		wt->outlength = wt->siglength;
		for (i = 0; i < wt->siglength; ++i) {
			wt->output[i] = inp[i];
		}
		wt->nodes = 1;
		wt->nodeindex[0] = 0;
		wt->nodeindex[1] = 0;
		wt->numnodeslevel[0] = 1;
	}
	else {
		for (i = J; i > 0; --i) { 
			llb = ipow2(i);
			it1 -= llb;
			wt->numnodeslevel[i] = 0;
			for (j = 0; j < llb; ++j) {
				if (wt->basisvector[it1 + j] == 1) {
					//printf("NODE %d %d %d \n", i, j, wt->length[J - i + 1]);
					wt->nodeindex[2 * wt->nodes] = i;
					wt->nodeindex[2 * wt->nodes + 1] = j;
					wt->nodes += 1;
					wt->numnodeslevel[i] += 1;
					for (k = 0; k < wt->length[J - i + 1]; ++k) {
						wt->output[it2 + k] = tree[nodelength[it1 - 1 + j] + k];// access tree
					}
					it2 += wt->length[J - i + 1];
				}
			}
		}
		wt->outlength = it2;
	}

	wt->coeflength[0] = wt->siglength;

	for (i = 1; i < J + 1; ++i) {
		wt->coeflength[i] = wt->length[J - i + 1];
	}

	free(orig);
	free(tree);
	free(nodelength);
}

void getDWPTCoeffs(wpt_object wt, int X, int Y, double *coeffs, int N) {
	int ymax, i;
	int np,citer;
	int flag;

	if (X <= 0 || X > wt->J) {
		printf("X co-ordinate must be >= 1 and <= %d", wt->J);
		exit(-1);
	}
	ymax = 1;
	for (i = 0; i < X; ++i) {
		ymax *= 2;
	}

	ymax -= 1;

	if (Y < 0 || Y > ymax) {
		printf("Y co-ordinate must be >= 0 and <= %d", ymax);
		exit(-1);
	}

	np = 0;
	citer = 0;

	for (i = wt->J; i > X; --i) {
		np += wt->numnodeslevel[i];
		citer += wt->numnodeslevel[i] * wt->coeflength[i];
	}

	i = 0;
	flag = 0;
	for (i = 0; i < wt->numnodeslevel[X]; ++i) {
		if (wt->nodeindex[2 * np + 1] == Y) {
			flag = 1;
			break;
		}
		np++;
		citer += wt->coeflength[X];
	}

	if (flag == 0) {
		printf("The Node is Not Part Of The Best Basis Tree Use wpt_summary function to list available nodes \n");
		exit(-1);
	}

	for (i = 0; i < N; ++i) {
		coeffs[i] = wt->output[citer + i];
	}

}

static void idwpt_per(wpt_object wt, double *cA, int len_cA, double *cD, double *X) {
	int len_avg, i, l, m, n, t, l2;

	len_avg = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	l2 = len_avg / 2;
	m = -2;
	n = -1;

	for (i = 0; i < len_cA + l2 - 1; ++i) {
		m += 2;
		n += 2;
		X[m] = 0.0;
		X[n] = 0.0;
		for (l = 0; l < l2; ++l) {
			t = 2 * l;
			if ((i - l) >= 0 && (i - l) < len_cA) {
				X[m] += wt->wave->lpr[t] * cA[i - l] + wt->wave->hpr[t] * cD[i - l];
				X[n] += wt->wave->lpr[t + 1] * cA[i - l] + wt->wave->hpr[t + 1] * cD[i - l];
			}
			else if ((i - l) >= len_cA && (i - l) < len_cA + len_avg - 1) {
				X[m] += wt->wave->lpr[t] * cA[i - l - len_cA] + wt->wave->hpr[t] * cD[i - l - len_cA];
				X[n] += wt->wave->lpr[t + 1] * cA[i - l - len_cA] + wt->wave->hpr[t + 1] * cD[i - l - len_cA];
			}
			else if ((i - l) < 0 && (i - l) > -l2) {
				X[m] += wt->wave->lpr[t] * cA[len_cA + i - l] + wt->wave->hpr[t] * cD[len_cA + i - l];
				X[n] += wt->wave->lpr[t + 1] * cA[len_cA + i - l] + wt->wave->hpr[t + 1] * cD[len_cA + i - l];
			}
		}
	}
}

static void idwpt_sym(wpt_object wt, double *cA, int len_cA, double *cD, double *X) {
	int len_avg, i, m, n, t, v;
	len_avg = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	m = -2;
	n = -1;

	for (v = 0; v < len_cA; ++v) {
		i = v;
		m += 2;
		n += 2;
		X[m] = 0.0;
		X[n] = 0.0;

		for (int l = MAX(i-len_cA, 0); l < MIN(i+1, len_avg / 2); ++l) {
			t = 2 * l;

			X[m] += wt->wave->lpr[t] * cA[i-l] + wt->wave->hpr[t] * cD[i-l];
			X[n] += wt->wave->lpr[t + 1] * cA[i-l] + wt->wave->hpr[t + 1] * cD[i-l];
		}
	}
}

void idwpt(wpt_object wt, double *dwtop) {
	int J, i, lf, k,p,l;
	int app_len, det_len, index, n1, llb, index2, index3, index4,indexp,xlen;
	double *X_lp, *X,  *out, *out2;
	int *prep,*ptemp;
	J = wt->J;
	app_len = wt->length[0];
	p = ipow2(J);
	lf = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	xlen = p * (app_len + 2 * lf);

	X_lp = (double*)malloc(sizeof(double)* 2 * (wt->length[J] + lf));
	X = (double*)malloc(sizeof(double)* xlen);
	out = (double*)malloc(sizeof(double)* wt->length[J]);
	out2 = (double*)malloc(sizeof(double)* wt->length[J]);
	prep = (int*)malloc(sizeof(int)* p);
	ptemp = (int*)malloc(sizeof(int)* p);
	n1 = 1;
	llb = 1;
	index2 = xlen / p;
	indexp = 0;
	if (wt->basisvector[0] == 1) {
		for (i = 0; i < wt->siglength; ++i) {
			dwtop[i] = wt->output[i];
		}

	}
	else {
		for (i = 0; i < J; ++i) {
			llb *= 2;
			n1 += llb;
		}

		for (i = 0; i < xlen; ++i) {
			X[i] = 0.0;
		}

		for (i = 0; i < llb; ++i) {
			prep[i] = (int)wt->basisvector[n1 - llb + i];
			ptemp[i] = 0;
		}

		if (!strcmp(wt->ext, "per")) {
			app_len = wt->length[0];
			det_len = wt->length[1];
			index = 0;


			for (i = 0; i < J; ++i) {
				p = ipow2(J - i - 1);
				det_len = wt->length[i + 1];
				index2 *= 2;
				index3 = 0;
				index4 = 0;
				//idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);
				n1 -= llb;
				for (l = 0; l < llb; ++l) {
					if (ptemp[l] != 2) {
						prep[l] = (int)wt->basisvector[n1 + l];
					}
					else {
						prep[l] = ptemp[l];
					}
					ptemp[l] = 0;
				}


				for (l = 0; l < p; ++l) {
					if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
						for (k = 0; k < det_len; ++k) {
							out[k] = wt->output[index + k];
							out2[k] = wt->output[index + det_len + k];
						}
						idwpt_per(wt, out, det_len, out2, X_lp);
						for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
							X[index3 + k - lf / 2 + 1] = X_lp[k];
						}
						index += 2 * det_len;
						index3 += index2;
						index4 += 2 * indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
						index4 += indexp;
						for (k = 0; k < det_len; ++k) {
							out[k] = wt->output[index + k];
							out2[k] = X[index4 + k];
						}
						idwpt_per(wt, out, det_len, out2, X_lp);
						for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
							X[index3 + k - lf / 2 + 1] = X_lp[k];
						}
						index += det_len;
						index3 += index2;
						index4 += indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
						for (k = 0; k < det_len; ++k) {
							out[k] = X[index4 + k];
							out2[k] = wt->output[index + k];
						}
						idwpt_per(wt, out, det_len, out2, X_lp);
						for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
							X[index3 + k - lf / 2 + 1] = X_lp[k];
						}
						index += det_len;
						index3 += index2;
						index4 += 2 * indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
						for (k = 0; k < det_len; ++k) {
							out[k] = X[index4 + k];
							out2[k] = X[index4 + indexp + k];
						}
						idwpt_per(wt, out, det_len, out2, X_lp);
						for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
							X[index3 + k - lf / 2 + 1] = X_lp[k];
						}
						index4 += 2 * indexp;
						index3 += index2;
						ptemp[l] = 2;
					}
					else {
						index3 += index2;
						index4 += 2 * indexp;
					}

				}


				/*
				idwt_per(wt, out, det_len, wt->output + iter, det_len, X_lp);
				for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
				out[k - lf / 2 + 1] = X_lp[k];
				}

				iter += det_len;
				det_len = wt->length[i + 2];
				*/
				llb /= 2;
				indexp = index2;
			}

			//free(X_lp);

		}
		else if (!strcmp(wt->ext, "sym")) {
			app_len = wt->length[0];
			det_len = wt->length[1];

			//X_lp = (double*)malloc(sizeof(double)* (N + 2 * lf - 1));
			index = 0;

			for (i = 0; i < J; ++i) {
				p = ipow2(J - i - 1);
				det_len = wt->length[i + 1];
				index2 *= 2;
				index3 = 0;
				index4 = 0;
				//idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);
				n1 -= llb;
				for (l = 0; l < llb; ++l) {
					if (ptemp[l] != 2) {
						prep[l] = (int)wt->basisvector[n1 + l];
					}
					else {
						prep[l] = ptemp[l];
					}
					ptemp[l] = 0;
				}


				for (l = 0; l < p; ++l) {
					if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
						for (k = 0; k < det_len; ++k) {
							out[k] = wt->output[index + k];
							out2[k] = wt->output[index + det_len + k];
						}
						idwpt_sym(wt, out, det_len, out2, X_lp);
						for (k = lf - 2; k < 2 * det_len; ++k) {
							X[index3 + k - lf + 2] = X_lp[k];
						}
						index += 2 * det_len;
						index3 += index2;
						index4 += 2 * indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
						index4 += indexp;
						for (k = 0; k < det_len; ++k) {
							out[k] = wt->output[index + k];
							out2[k] = X[index4 + k];
						}
						idwpt_sym(wt, out, det_len, out2, X_lp);
						for (k = lf - 2; k < 2 * det_len; ++k) {
							X[index3 + k - lf + 2] = X_lp[k];
						}
						index += det_len;
						index3 += index2;
						index4 += indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
						for (k = 0; k < det_len; ++k) {
							out[k] = X[index4 + k];
							out2[k] = wt->output[index + k];
						}
						idwpt_sym(wt, out, det_len, out2, X_lp);
						for (k = lf - 2; k < 2 * det_len; ++k) {
							X[index3 + k - lf + 2] = X_lp[k];
						}
						index += det_len;
						index3 += index2;
						index4 += 2 * indexp;
						ptemp[l] = 2;
					}
					else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
						for (k = 0; k < det_len; ++k) {
							out[k] = X[index4 + k];
							out2[k] = X[index4 + indexp + k];
						}
						idwpt_sym(wt, out, det_len, out2, X_lp);
						for (k = lf - 2; k < 2 * det_len; ++k) {
							X[index3 + k - lf + 2] = X_lp[k];
						}
						index4 += 2 * indexp;
						index3 += index2;
						ptemp[l] = 2;
					}
					else {
						index3 += index2;
						index4 += 2 * indexp;
					}

				}

				//idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);
				/*
				idwpt_sym(wt, out, det_len, wt->output + iter, det_len, X_lp);
				for (k = lf - 2; k < 2 * det_len; ++k) {
				out[k - lf + 2] = X_lp[k];
				}

				iter += det_len;
				det_len = wt->length[i + 2];
				*/
				llb /= 2;
				indexp = index2;
			}

			//free(X_lp);

		}
		else {
			printf("Signal extension can be either per or sym");
			exit(-1);
		}

		for (i = 0; i < wt->siglength; ++i) {
			//printf("%g ", X[i]);
			dwtop[i] = X[i];
		}

	}


	free(out);
	free(X_lp);
	free(X);
	free(out2);
	free(prep);
	free(ptemp);
}

void setDWPTExtension(wpt_object wt, const char *extension) {
	if (!strcmp(extension, "sym")) {
		strcpy(wt->ext, "sym");
	}
	else if (!strcmp(extension, "per")) {
		strcpy(wt->ext, "per");
	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}
}

void setDWPTEntropy(wpt_object wt, const char *entropy, double eparam) {
	if (!strcmp(entropy, "shannon")) {
		strcpy(wt->entropy, "shannon");
	}
	else if (!strcmp(entropy, "threshold")) {
		strcpy(wt->entropy, "threshold");
		wt ->eparam = eparam;
	}
	else if (!strcmp(entropy, "norm")) {
		strcpy(wt->entropy, "norm");
		wt->eparam = eparam;
	}
	else if (!strcmp(entropy, "logenergy") || !strcmp(entropy, "log energy") || !strcmp(entropy, "energy")) {
		strcpy(wt->entropy, "logenergy");
	}
	else {
		printf("Entropy should be one of shannon, threshold, norm or logenergy");
		exit(-1);
	}
}

void wpt_summary(wpt_object wt) {
	int i, k, p2;
	int J, it1,it2;
	J = wt->J;
	wave_summary(wt->wave);
	printf("\n");
	printf("Signal Extension : %s \n", wt->ext);
	printf("\n");
	printf("Entropy : %s \n", wt->entropy);
	printf("\n");
	printf("Number of Decomposition Levels %d \n", wt->J);
	printf("\n");
	printf("Number of Active Nodes %d \n", wt->nodes);
	printf("\n");
	printf("Length of Input Signal %d \n", wt->siglength);
	printf("\n");
	printf("Length of WT Output Vector %d \n", wt->outlength);
	printf("\n");
	printf("Wavelet Coefficients are contained in vector : %s \n", "output");
	printf("\n");
	printf("Coefficients Access \n");
	it1 = 1;
	it2 = 0;
	for (i = 0; i < J; ++i) {
		it1 += ipow2(i + 1);
	}
	for (i = J; i > 0; --i) {
		p2 = ipow2(i);
		it1 -= p2;
		for (k = 0; k < p2; ++k) {
			if (wt->basisvector[it1 + k] == 1) {
				printf("Node %d %d Access : output[%d] Length : %d \n", i, k, it2, wt->length[J - i + 1]);
				it2 += wt->length[J - i + 1];
			}
		}
	}

	printf("\n");

}

void wave_summary(wave_object obj) {
	int i,N;
	N = obj->filtlength;
	printf("\n");
	printf("Wavelet Name : %s \n",obj->wname);
	printf("\n");
	printf("Wavelet Filters \n\n");
	printf("lpd : [");
	for (i = 0; i < N-1; ++i) {
		printf("%g,", obj->lpd[i]);
	}
	printf("%g", obj->lpd[N-1]);
	printf("] \n\n");
	printf("hpd : [");
	for (i = 0; i < N-1; ++i) {
		printf("%g,", obj->hpd[i]);
	}
	printf("%g", obj->hpd[N - 1]);
	printf("] \n\n");
	printf("lpr : [");
	for (i = 0; i < N-1; ++i) {
		printf("%g,", obj->lpr[i]);
	}
	printf("%g", obj->lpr[N - 1]);
	printf("] \n\n");
	printf("hpr : [");
	for (i = 0; i < N-1; ++i) {
		printf("%g,", obj->hpr[i]);
	}
	printf("%g", obj->hpr[N - 1]);
	printf("] \n\n");
}

void wpt_free(wpt_object object) {
	free(object);
}

void wave_free(wave_object object) {
	free(object);
}
