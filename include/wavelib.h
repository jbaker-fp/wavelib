#ifndef WAVELIB_H_
#define WAVELIB_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_MSC_VER)
#pragma warning(disable : 4200)
#pragma warning(disable : 4996)
#endif

#ifndef fft_type
#define fft_type double
#endif

#ifndef cplx_type
#define cplx_type double
#endif


typedef struct cplx_t {
	cplx_type re;
	cplx_type im;
} cplx_data;

typedef struct wave_set* wave_object;

wave_object wave_init(const char* wname);

struct wave_set{
	char wname[50];
	int filtlength;// When all filters are of the same length. [Matlab uses zero-padding to make all filters of the same length]
	int lpd_len;// Default filtlength = lpd_len = lpr_len = hpd_len = hpr_len
	int hpd_len;
	int lpr_len;
	int hpr_len;
	double *lpd;
	double *hpd;
	double *lpr;
	double *hpr;
	double params[0];
};

typedef struct fft_t {
  fft_type re;
  fft_type im;
} fft_data;

typedef struct fft_set* fft_object;

fft_object fft_init(int N, int sgn);

struct fft_set{
	int N;
	int sgn;
	int factors[64];
	int lf;
	int lt;
	fft_data twiddle[1];
};

typedef struct fft_real_set* fft_real_object;

fft_real_object fft_real_init(int N, int sgn);

struct fft_real_set{
	fft_object cobj;
	fft_data twiddle2[1];
};

typedef struct conv_set* conv_object;

conv_object conv_init(int N, int L);

struct conv_set{
	fft_real_object fobj;
	fft_real_object iobj;
	int ilen1;
	int ilen2;
	int clen;
};

typedef struct wpt_set* wpt_object;

wpt_object wpt_init(wave_object wave, int siglength, int J);

struct wpt_set{
	wave_object wave;
	conv_object cobj;
	int siglength;// Length of the original signal.
	int outlength;// Length of the output DWT vector
	int lenlength;// Length of the Output Dimension Vector "length"
	int J; // Number of decomposition Levels
	int MaxIter;// Maximum Iterations J <= MaxIter
	int even;// even = 1 if signal is of even length. even = 0 otherwise
	char ext[10];// Type of Extension used - "per" or "sym"
	char entropy[20];
	double eparam;

	int N; //
	int nodes;
	int length[102];
	double *output;
	double *costvalues;
	double *basisvector;
	int *nodeindex;
	int *numnodeslevel;
	int *coeflength;
	double params[0];
};

void dwpt(wpt_object wt, const double *inp);

void idwpt(wpt_object wt, double *dwtop);

void setDWPTExtension(wpt_object wt, const char *extension);

void setDWPTEntropy(wpt_object wt, const char *entropy, double eparam);

int getDWPTNodelength(wpt_object wt, int X);

void getDWPTCoeffs(wpt_object wt, int X, int Y, double *coeffs, int N);

void wpt_summary(wpt_object wt);

void wpt_free(wpt_object object);

void wave_summary(wave_object obj);

void wave_free(wave_object object);


#ifdef __cplusplus
}
#endif


#endif /* WAVELIB_H_ */
