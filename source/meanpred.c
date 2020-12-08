//gcc -g meanpreds.c -o meanpreds
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "my_misc.c"
#define isnum(c) (((c>='0') && (c<='9'))||(c=='-')||(c=='+')||(c==' '))
//static char pr[1024];
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
//#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
//#include "my_misc.c"
//#define GSL 1
#undef GSL
#define GSL 1
//#define GSL 0
int GSLtag;
void least_gsl_error_handler( const char * reason,
			      const char * file,
			      int line,
			      int gsl_errno)
{
  //  printf("%s:%d: %s (gsl_error: %s)\n", 
  //	 file, line, reason, gsl_strerror(gsl_errno));
  printf("%s:%d: %s (gsl_error %d)\n", 
	 file, line, reason,gsl_errno);
  GSLtag=1;
}
double VarY;
void calc_Ainvb(double *M, double *A_data, double *b_data, int nx, int ndata,int nxmarume)
{
  /* b= A M  --> M= Ainv b
     b  in R^{ndata x 1}
     A  in R^{ndata x nx}   (ndata>nx)
     M  in R^{nx    x 1}
   */
  int i;
  gsl_vector *S;
  gsl_vector *work;
  gsl_vector *x;
  gsl_vector *b;
  gsl_matrix *V;
  gsl_matrix *A;
  if(ndata<=0){
    M[0]=1;
    for(i=1;i<nx;i++) M[i]=1./nx;
  }
  else{
    A =gsl_matrix_alloc(ndata,nx);
    b =gsl_vector_alloc(ndata);b->data=b_data;
    A->data=A_data;  
    //  for(i=0;i<ndata;i++){printf("%e %e %d #bi,Ai\n",A->data[i],b->data[i],i); }
    
    V = gsl_matrix_alloc (nx,nx); 
    S = gsl_vector_alloc (nx); 
    work = gsl_vector_alloc(nx); 
    x = gsl_vector_alloc(nx); 
    
    gsl_set_error_handler( &least_gsl_error_handler );GSLtag=0;
    
#if GSL == 1
    gsl_linalg_SV_decomp_jacobi (A, V, S); //higer accuracy than Golub-Reinsh
#elif GSL == 2
    gsl_linalg_SV_decomp (A, V, S, work); //Golub-Reinsh
#elif GSL == 3
    {//faster
      gsl_matrix *X;
      X = gsl_matrix_alloc (nx,nx); 
      gsl_linalg_SV_decomp_mod (A, X, V, S, work); 
      gsl_matrix_free(X);
    }
#endif
    //    for(i=nxmarume;i<nx;i++){
    
    fprintf(stdout,"Lambda");
    for(i=0;i<nx;i++){
      //      if(nx>40) break;
      //      fprintf(stdout,"S->data[%d]=%e;",i,S->data[i]);
      fprintf(stdout,"(%d)%e",i,S->data[i]);
    }
    fprintf(stdout,"\n");
    VarY=0;
    if(GSLtag==1){
      fprintf(stderr,"\nEigenValue=");for(i=0;i<nx;i++) fprintf(stderr,"%e ",S->data[i]); fprintf(stderr,"\n");
      for(i=0;i<nx;i++) M[i]=1./nx;
    }
    else{
      if(nxmarume>=0){//Principal Component Analysis
	for(i=nxmarume;i<nx;i++) S->data[i]=0;
	for(i=0;i<nxmarume;i++) VarY+=S->data[i];
	gsl_linalg_SV_solve(A,V,S,b,x);
	for(i=0;i<nx;i++) M[i]=gsl_vector_get(x,i);
	//    {//check
	//      double err=0,yhat;
	//      int j;
	//      for(i=0;i<ndata;i++){
	//	if(i>10) break;
	//	yhat=b_data[i];
	//	for(j=0;j<nx;j++){
	//	  yhat-= M[j]*A_data[i*nx+j];
	//	}
	//	err+=square(yhat);
	//	fprintf(stderr,"err=%e=%e/%d. yhat=%e=(%e",err/(i+1.),err,i,yhat,b_data[i]);
	//	for(j=0;j<nx;j++) fprintf(stderr,"-(%e)*(%e)",M[j],A_data[i*nx+j]);
	//	fprintf(stderr,")^2\n");
	//      }
	//    }
      }
      else{//Minor Component Analysis
	int ii=0;
	int im=-nxmarume;
	for(i=nx-1;i>=0;i--){
	  if(S->data[i]>1e-20){
	    ii++;
	    VarY+=S->data[i];
	    if(ii>=im) break;
	  }
	}
	int j;

	double M0=-V->data[i];
	for(j=1;j<nx;j++){
	  M[j-1]=V->data[j*nx+i]/M0;
	}
//	double M0=-V->data[i*nx];
//	for(j=1;j<nx;j++){
//	  M[j-1]=V->data[i*nx+j]/M0;
//	}

      }
    }
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(x);
    gsl_matrix_free(V);
    gsl_matrix_free(A);
    gsl_vector_free(b);
  }
}
void calc_AtinvbMCA(double *M, double *At_data, double *b_data, int nx, int ndata,int nxmarume)
{ 
  /* b= A M  --> M= Ainv b
     b  in R^{ndata x 1}
     A=At  in R^{ndata x nx} 
     M  in R^{nx    x 1}
   */
  int nx1=nx+1,i,j,jnx1;
  gsl_matrix *A  =gsl_matrix_alloc(ndata,nx1);
  for(j=0;j<ndata;j++){
    jnx1=j*nx1;
    A->data[jnx1]=b_data[j];
    for(i=0;i<nx;i++){
      A->data[jnx1+i+1]=At_data[i*ndata+j];
    }
  }
  calc_Ainvb(M, (double *)(A->data), b_data, nx1, ndata, -nxmarume);
  
  gsl_matrix_free(A);
}
void calc_Atinvb(double *M, double *At_data, double *b_data, int nx, int ndata,int nxmarume)
{ 
  /* b= A M  --> M= Ainv b
     b  in R^{ndata x 1}
     A=At  in R^{ndata x nx} 
     M  in R^{nx    x 1}
   */
  gsl_matrix *At =gsl_matrix_alloc(nx,ndata);
  gsl_matrix *A  =gsl_matrix_alloc(ndata,nx);
  At->data=At_data;
  gsl_matrix_transpose_memcpy(A, At);
  calc_Ainvb(M, (double *)(A->data), b_data, nx, ndata, nxmarume);
  gsl_matrix_free(A);
  gsl_matrix_free(At);
}
void calc_yp(double *y, double *a_data, int nx, int ndata,int nxmarume)
{//M[i]=Ainv b
  int i, l,j;
  double sumS=0,sumSj=0;
  gsl_vector *S;
  gsl_vector *work;
  gsl_vector *x;
  //  gsl_vector *b;
  gsl_matrix *V;
  gsl_matrix *A;
  A =gsl_matrix_alloc(ndata,nx);
  A->data=a_data;  
  //    b =gsl_vector_alloc(ndata);
  //    b->data=b_data;
  V = gsl_matrix_alloc (nx,nx); 
  S = gsl_vector_alloc (nx); 
  work = gsl_vector_alloc(nx); 
  x = gsl_vector_alloc(nx); 
  
  gsl_set_error_handler( &least_gsl_error_handler );GSLtag=0;
  
#if GSL == 1
  gsl_linalg_SV_decomp_jacobi (A, V, S); //higer accuracy than Golub-Reinsh
#elif GSL == 2
  gsl_linalg_SV_decomp (A, V, S, work); //Golub-Reinsh
#elif GSL == 3
  {//faster
    gsl_matrix *X;
    X = gsl_matrix_alloc (nx,nx); 
    gsl_linalg_SV_decomp_mod (A, X, V, S, work); 
    gsl_matrix_free(X);
  }
#endif
  //    for(i=nxmarume;i<nx;i++){
  sumS=sumSj=0;
  for(l=0;l<nx;l++) sumS+=gsl_vector_get(S,l);
  if(nxmarume>=0){
    for(i=0;i<ndata;i++){
      y[i]=0;
      for(l=0;l<nxmarume;l++){
	//	fprintf(stderr,"S%d=%e.\n",l,gsl_vector_get(S,l));
	for(j=0;j<nx;j++){
	  //	y[i]+=S->data[l]*A->data[i*nx+l]*V->data[j*nx+l];
	  y[i]+=gsl_vector_get(S,l)*gsl_matrix_get(A,i,l)*gsl_matrix_get(V,j,l);
	}
      }
      y[i]/=nx;
      if(i<5){
	sumSj+=gsl_vector_get(S,i);
	fprintf(stderr,"[%d]y=%e, s=%e, CumulativeProp=%f=%f/%f\n",i,y[i],gsl_vector_get(S,i),sumSj/sumS,sumSj,sumS);
      }
    }
  }
  else{
    int ll=0;
    for(i=0;i<ndata;i++){
      y[i]=0;
      for(l=nx;;l--){
	if(gsl_vector_get(S,l)<=1e-30) continue;
	if(i==0) fprintf(stderr,"S%d=%e.\n",l,gsl_vector_get(S,l));
	for(j=0;j<nx;j++){
	  y[i]+=gsl_vector_get(S,l)*gsl_matrix_get(A,i,l)*gsl_matrix_get(V,j,l);
	}
	if(++ll>=-nxmarume) break;
      }
      y[i]/=nx;
      if(i<5){
	sumSj+=gsl_vector_get(S,i);
	if(i<5) fprintf(stderr,"[%d]y=%e, s=%e, CumulativeProp=%f=%f/%f\n",i,y[i],gsl_vector_get(S,i),sumSj/sumS,sumSj,sumS);
      }
    }
  }
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_vector_free(x);
  gsl_matrix_free(V);
}

////For Maximizing the Entropy of the training data
////For Minimizing the minus of the Entropy of the training data
//Lagrange Method
int
my_df_Entropy(const gsl_vector *v, void *params, gsl_vector *df)
{
  double **dp = (double **)params;
  int nd=(int)((double)*dp[0]);
  int nl=(int)((double)*dp[1]);
  double *y=(double *)dp[2];
  double *f=(double *)dp[3];
  int sgn=(int)((double)*dp[4]);
  double *yp=(double*)malloc(sizeof(double)*nd);
  double *e2=(double*)malloc(sizeof(double)*nd);
  double *e =(double*)malloc(sizeof(double)*nd);
  double *dHdp=(double*)malloc(sizeof(double)*nd);
  double *dpde2=(double*)malloc(sizeof(double)*nd*nd);
  double E2,lognd=log(nd),pi;
  int i,j,k,l;

  E2=0;
  for(i=0;i<nd;i++){
    yp[i]=0; for(j=0;j<nl;j++) yp[i]+=square(v->data[j])*f[j*nd+i];
    e[i]=(y[i]-yp[i]);
    e2[i]=e[i]*e[i];
    E2+=e2[i];
  }
  for(i=0;i<nd;i++){
    pi=e2[i]/E2;
    dHdp[i]=(log(pi)+1.)/lognd;//for the minus Entropy
    for(k=0;k<nd;k++){
      if(i == k) dpde2[i*nd+k]=1./E2-e2[i]/E2/E2;
      else dpde2[i*nd+k]=-e2[i]/E2/E2;
    }
  }
  //  double dfnorm=0;
  //  fprintf(stdout,"df=");
  for(l=0;l<nl;l++){
    df->data[l]=0;
    for(i=0;i<nd;i++){
      for(k=0;k<nd;k++){
	df->data[l]+=dHdp[i]*dpde2[i*nd+k]*(2.*e[k]*(-2.*v->data[l]*f[l*nd+k]));
	//	df->data[l]+=dHdp[i]*dpde2[i*nd+k]*(-2.*e[k]*f[l*nd+k]);
      }
      df->data[l]+= 2.*v->data[nl]*v->data[l]*sgn;//for Penalty
      //      df->data[l]-= 2.*v->data[nl]*v->data[l];//for Penalty
      //      df->data[l] = C0*(-2.*e[i]*f[l*nd+i])/E0 +C1*df->data[l];//for H=C0*E2/E0+C1*Entropy
    }
    //    dfnorm+=square(df->data[l]);
    //    fprintf(stdout,"%.4f ",df->data[l]);
  }
  //Lagrange for sumv=1. v->data[nd]=lambda 
  double sumv=0;//for Lagrange
  for(j=0;j<nl;j++) sumv+=square(v->data[j]);//for Lagrange
  df->data[nl]=(sumv-1.)*sgn;//for Lagrange L=H+lambda*(sum c_l-1)
  //  df->data[nl]=1.-sumv;//for Lagrange L=H-lambda*(sum c_l-1)
  //  fprintf(stdout,";dfnorm=%.4e\n",dfnorm);
  free(yp);
  free(e);
  free(e2);
  return GSL_SUCCESS;
}
int
gsl_print_state (size_t iter, gsl_multiroot_fsolver * s, void *params)
{
  double **dp = (double **)params;
  int nd=(int)((double)*dp[0]);
  int nl=(int)((double)*dp[1]);
  double *y=(double *)dp[2];
  double *f=(double *)dp[3];
  double *yp=(double*)malloc(sizeof(double)*nd);
  double *e2=(double*)malloc(sizeof(double)*nd);
  double *e =(double*)malloc(sizeof(double)*nd);
  double E2,pi;
  int i,j;
  
  E2=0;
  for(i=0;i<nd;i++){
    yp[i]=0; for(j=0;j<nl;j++) yp[i]+=square(gsl_vector_get (s->x, j))*f[j*nd+i];
    e[i]=(y[i]-yp[i]);
    e2[i]=e[i]*e[i];
    E2+=e2[i];
  }
  double H=0;//
  for(i=0;i<nd;i++){
    pi=e2[i]/E2;
    H += pi*log(pi);//the minus Entropy
  }
  H/=log(nd);
  printf ("iter = %3u H=%f x0=%.3e...%.3e "
	  "dL/dx=%.3e...%3e\n",
	  (unsigned int)iter,-H,
	  square(gsl_vector_get (s->x, 0)),
	  gsl_vector_get (s->x, nl),
	  gsl_vector_get (s->f, 0),
	  gsl_vector_get (s->f, nl));
  return 0;
}

//////For Maximizing the Entropy of the training data
//////For Minimizing the minus of the Entropy of the training data
//double
//my_f_Entropy(const gsl_vector *v, void *params)
//{
//  //params
//  /*
//    e_i=y_i- sum_l v[l] f[l][i] (l=0:nx-1)
//    p_i= e_i^2/ sum_j e_j^2
//    H= sum_i p_i log(p_i);
//
//  */
//  double **dp = (double **)params;
//  int nd=(int)((double)*dp[0]);
//  int nl=(int)((double)*dp[1]);
//  double *y=(double *)dp[2];
//  double *f=(double *)dp[3];
//  double C0=(double)*dp[4];
////  double C1=(double)*dp[5];
////  double E0=(double)*dp[6];
//  double *e2=(double*)malloc(sizeof(double)*nd);
//  double E2;
//  double H,ypi,pi;
//  int i,j;
//  E2=0;
//  for(i=0;i<nd;i++){
//    ypi=0; for(j=0;j<nl;j++) ypi+=square(v->data[j])*f[j*nd+i];
//    e2[i]=square(ypi-y[i]);
//    E2+=e2[i];
//  }
//  H=0;
//  for(i=0;i<nd;i++){
//    pi=e2[i]/E2;
//    H += pi*log(pi);//the minus Entropy
//  }
//  H/=log(nd);
//
//  //  H=C0*E2/E0+C1*H;//for H=C0*E2/E0+C1*Entropy
//
//  //Lagrange for sumv=1. v->data[nd]=lambda 
//  double sumv=0;//for Lagrange
//  for(j=0;j<nl;j++) sumv+=square(v->data[j]);//for Lagrange
//  //  if(C0*(sumv-1.)<0) C0*=-1;
//  //  H+=v->data[nl]*(sumv-1.);//for Lagrange 
//  fprintf(stderr,"H=%f,sumv=%e",H,sumv);
//  H+=C0*(sumv-1.);//for Penalty method
//  /////
//  free(e2);
//  return H;
//}
//
///* The gradient of f, df = (df/dx, df/dy). */
//void
//my_df_Entropy(const gsl_vector *v, void *params,
//       gsl_vector *df)
//{
//  double **dp = (double **)params;
//  int nd=(int)((double)*dp[0]);
//  int nl=(int)((double)*dp[1]);
//  double *y=(double *)dp[2];
//  double *f=(double *)dp[3];
//  double C0=(double)*dp[4];
//  //  double C1=(double)*dp[5];
//  //  double E0=(double)*dp[6];
//  double *yp=(double*)malloc(sizeof(double)*nd);
//  double *e2=(double*)malloc(sizeof(double)*nd);
//  double *e =(double*)malloc(sizeof(double)*nd);
//  double *dHdp=(double*)malloc(sizeof(double)*nd);
//  double *dpde2=(double*)malloc(sizeof(double)*nd*nd);
//  double E2,lognd=log(nd),pi;
//  int i,j,k,l;
////  {
////    double sumv;
////    for(j=0;j<nl;j++) sumv+=square(v->data[j]);//for Lagrange
////    if(C0*(sumv-1.)<0) C0*=-1;
////  }
//  E2=0;
//  for(i=0;i<nd;i++){
//    yp[i]=0; for(j=0;j<nl;j++) yp[i]+=square(v->data[j])*f[j*nd+i];
//    e[i]=(y[i]-yp[i]);
//    e2[i]=e[i]*e[i];
//    E2+=e2[i];
//  }
//  for(i=0;i<nd;i++){
//    pi=e2[i]/E2;
//    //    dHdp[i]=(-log(pi)-1.)/lognd;
//    dHdp[i]=(log(pi)+1.)/lognd;//for the minus Entropy
//    for(k=0;k<nd;k++){
//      if(i == k) dpde2[i*nd+k]=1./E2-e2[i]/E2/E2;
//      else dpde2[i*nd+k]=-e2[i]/E2/E2;
//    }
//  }
//  double dfnorm=0;
//  //  fprintf(stdout,"df=");
//  for(l=0;l<nl;l++){
//    df->data[l]=0;
//    for(i=0;i<nd;i++){
//      for(k=0;k<nd;k++){
//	df->data[l]+=dHdp[i]*dpde2[i*nd+k]*(2.*e[k]*(-2.*v->data[l]*f[l*nd+k]));
//	//	df->data[l]+=dHdp[i]*dpde2[i*nd+k]*(-2.*e[k]*f[l*nd+k]);
//      }
//      df->data[l]+= 2.*C0*v->data[l];//for Penalty
//      //      df->data[l] = C0*(-2.*e[i]*f[l*nd+i])/E0 +C1*df->data[l];//for H=C0*E2/E0+C1*Entropy
//    }
//    dfnorm+=square(df->data[l]);
//    //    fprintf(stdout,"%.4f ",df->data[l]);
//  }
//  //Lagrange for sumv=1. v->data[nd]=lambda 
//  //  double sumv=0;//for Lagrange
//  //  for(j=0;j<nl;j++) sumv+=square(v->data[j]);//for Lagrange
//  //  df->data[nl]=sumv-1.;//for Lagrange
//
//  fprintf(stdout,";dfnorm=%.4e\n",dfnorm);
//  free(yp);
//  free(e);
//  free(e2);
//}
//
///* Compute both f and df together. */
//void
//my_fdf_Entropy(const gsl_vector *x, void *params,
//	double *f, gsl_vector *df)
//{
//  *f = my_f_Entropy(x, params);
//  my_df_Entropy(x, params, df);
//}

////
typedef struct{
  int r1;
  int r2;
  double r3;
  int nr;
  double r12;
  double *r;
  double ymin;
  double ymax;
} RESOLUTION;
#define ZERO 1e-20
int init_res(RESOLUTION *res)
{
  int i;
  char buff[32];
  if(res->r2==0) return(-1);
  if(fabs(res->ymax-res->ymin)<ZERO) {
    res->r2=0;
    return(-1);
  }
  if(res->r3<-ZERO){//negrect resolution for the meanpred
    res->r2=0;
    return(-1);
  }
  res->r12=(double)res->r1/res->r2;
  //  res->nr=(res->ymax-res->ymin)/res->r12;
  res->nr=(res->ymax-res->ymin)/res->r12+0.5;//??
  //  res->r12=(res->ymax-res->ymin)/res->nr;//for ymin and ymax more reliable than r1 and r2 for ijcnn06
  res->r=(double*)malloc(sizeof(double)*(res->nr+1));
  if(fabs(res->r3)>ZERO){
    for(i=0;i<=res->nr;i++){
      sprintf(buff,"%.10le",((int)((i*res->r12+res->ymin)/res->r3+0.5))*res->r3);//r3桁で四捨五入
      sscanf(buff,"%lf",&res->r[i]);
    }
  }
  else{
    for(i=0;i<=res->nr;i++){
      sprintf(buff,"%.10le",(double)i*res->r12+res->ymin);
      sscanf(buff,"%lf",&res->r[i]);
    }
  }
  return(0);
}
double y_res(double y,RESOLUTION *res)
{
  int ii;
  if(res->r2==0) return(y);
  //  ii=(int)((y-ymin)/r12+0.5);//4sha5nyu
  ii=floor((y-res->ymin)/res->r12+0.5);//4sha5nyu
  if(y<res->ymin) return((double)res->r[0]);
  if(y>res->ymax) return((double)res->r[res->nr]);
  return((double)res->r[ii]);
}
#define buffsize 1024
int main(int argc, char **argv)
{
  int i,j,l;
  FILE *fp;
  char buff[buffsize];
  double *yp,_yp;
  double *y,_y,*yr,_yr;
  int *n,_n;  
  int num=0,num0=0;
  //  char *fnpredict="predict+.dat";
  //  char *fnloss="loss+.dat";
  int intmode=0, j0=1,nfiles=argc-1;
  char fn_is[256],fn_pred[256],fn_pred0[256],fn_prob[256];
  double _is,_mseis,_n_msetrain,_msetrain,_n_cells2,_nCells1,_n_mse,_mse,_msetrainis;
  double *prob,*sigma2,meansigma2,probtotal;
  FILE *fp0;
  int bayesmode=0; //0 (for flat probability) or 1 (for probatility decaying according to MSE).
  double bayesopt1=0;//for adjusting probability (1 for same probatility)
  double bayesopt2=0;//0 for original L, 1 for L_derivative requiring nfiles>=3.
  double bayesopt3=0;//??
  double bayesopt4=0;//??
  //  double *ds;
  //  int *s;
  double Lvar,Lvarval,Lval,LD=0,_LD;//
  double skew,kurtosis,entropy,skew2,skewa;//skew=3rd moment, kurtosis=4th moment
  double m2,m3,m4,m5,m6,m3a,m5a;//5th moment, 6th moment
  double m3p,m3m;//m3plus,m3minus
  int nm3p,nm3m;
  double mcov;//mean covariance
  double Lvar0;
  double Ltrain,LtrainE,H0;
  double *var,*vv;
  int *cc,ccmax=0;
  double _var,_vv,_cc;
  double *kyorix2,_kyorix2,_e2;
  int LDmode=2;
  RESOLUTION res={0,0,0,};
  double *yp00;
  int nop=0;
  if(argc<2){
    fprintf(stderr,"Usage: %s [ib:<intmode>:<basesmode>[:<probfile>]] <pred&is-body#1> <pred&is-body#1> ... \n",argv[0]);
    fprintf(stderr,"<intmode>:1 for intmode, 0 for real.\n");
    fprintf(stderr,"<bayesmode>:0 simple ensemble.?\n");
    fprintf(stderr,"<bayesmode>:1 for bayes mode1 exp.\n");
    fprintf(stderr,"<bayesmode>:2 for bayes mode2 YpInvY.\n");
    fprintf(stderr,"<bayesmode>:-1 for reading from probfile.\n");
    fprintf(stderr,"res:<r1>:<r2>:<ymin>:<ymax>.\n");
    fprintf(stderr,"Old Usage: %s [<int>:<bases>] <pred&is-body#1> <pred&is-body#1> ... \n",argv[0]);
    fprintf(stderr,"par:paramfile.\n");
    exit(1);
  }
  //  fprintf(stderr,"last arg[%d]='%s'\n",argc,argv[argc-1]);
  int parfile=0;
  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"par:",4)==0){//for compatibility to old version
      parfile=1;
      break;
    }
  }
  char *_argv0;
  if(parfile==1){
    char *fn=&argv[i][4];
    FILE *fp;
    if((fp=fopen(fn,"r"))==NULL){
      fprintf(stderr,"Param file (%s) open error.\n",fn);
      return(-1);
    }
    
    fseek( fp, 0L, SEEK_END );
    int fsize=ftell(fp);
    fseek( fp, 0L, SEEK_SET );//rewind(fp);//
    //    fprintf(stderr,"fsize=%d\n",fsize);
    _argv0=(char*)malloc(sizeof(char)*fsize);
    fgets(_argv0,fsize,fp);
    fclose(fp);
    argc=1;
    char *p=_argv0;
    for(j=0;j<fsize;j++) if(p[j]==' ') argc++;
    argv=(char**)malloc(sizeof(char*)*(argc));
    for(j=1;j<argc;j++){
      for(p++;;p++) if(*p==' ') {*p=0;break;}
      for(p++;;p++) if(*p!=' ') break;
      argv[j]=p;
    }
  }
  intmode=0;
  bayesmode=0;
  j0=1;
  nfiles=argc-1;
  for(i=j0;i<argc;i++){
    if(strncmp(argv[i],"int",3)==0){//for compatibility to old version
      intmode=1;
      j0++; 
      nfiles--;
    }
    else if(strncmp(argv[i],"ib:",3)==0){
      sscanf(&(argv[i][3]),"%d:%d:%lf:%lf:%lf:%lf",&intmode,&bayesmode,&bayesopt1,&bayesopt2,&bayesopt3,&bayesopt4);
      if(bayesmode==-1) sprintf(fn_prob,"%s",&(argv[1][8]));
      j0++; 
      nfiles--;
    }
    else if(strncmp(argv[i],"ry:",3)==0){
      sscanf(&argv[i][3],"%d:%d:%lf:%lf:%lf",&res.r1,&res.r2,&res.r3,&res.ymin,&res.ymax);
      j0++;
      nfiles--;
      init_res(&res);
    }
    else if(strncmp(argv[i],"LDm:",4)==0){
      sscanf(&argv[i][4],"%d",&LDmode);
      j0++;
      nfiles--;
      init_res(&res);
    }
    else if(strncmp(argv[i],"nop:",4)==0){
      sscanf(&argv[i][4],"%d",&nop);
      j0++;
      nfiles--;
    }
  }
  ////

  int (*net_printf)(char *);
  if(nop==0) net_printf=printf1;
  else net_printf=noprintf;

  prob =(double*) malloc(sizeof(double)*nfiles);
  sigma2=(double*) malloc(sizeof(double)*nfiles);
  meansigma2=0;
  for(j=0;j<nfiles;j++) argv[j+j0][strlen(argv[j+j0])-8]=0; //argv[j+j0]=fn_bodypred.dat
  for(j=0;j<nfiles;j++){
    sprintf(fn_is,"%sis.dat",argv[j+j0]);
    if((fp=fopen(fn_is,"r"))==NULL){
      fprintf(stderr,"file(%s) could not be opened!\n",fn_is);
    }
    else{
      fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",&_is,&_mseis,&_n_msetrain,&_msetrain,&_n_cells2,&_nCells1,&_n_mse,&_mse,&_msetrainis);
      //      num=_n_msetrain; sigma2[j]=_msetrain; meansigma2+=_msetrain;//050902??
      //      num=_n_mse; sigma2[j]=_mse; meansigma2+=_mse;//050823
      num=_n_mse; sigma2[j]=_msetrain; meansigma2+=_msetrain;//050914
      num0=_n_msetrain;
      //      fprintf(stderr,"sigma2[%d]%e\n",j,sigma2[j]);
      fclose(fp);
    }
  }
  meansigma2/=nfiles;
  
  ////////////////////////////////
  if(bayesmode>=100){}//end of bayesmode>=100
  else{//bayesmode <100
    if(bayesmode==0){//same probability
      probtotal=nfiles;
      for(j=0;j<nfiles;j++)  prob[j]=1./probtotal;
      probtotal=0; for(j=0;j<nfiles;j++) probtotal+=prob[j];
    }
    else if(bayesmode==1){//probatility depending on MSE(=sigma2)
      int jj;double probj;//probj is 1/prob[j] for avoiding marumegosa.
      probtotal=0;
      for(j=0;j<nfiles;j++){
	probj=0;
	for(jj=0;jj<nfiles;jj++){
	  if(jj==j) probj+=1;
	  else probj+= exp((sigma2[j]-sigma2[jj])*bayesopt1/meansigma2);
	}
	prob[j]=1./probj;
	probtotal+=prob[j];
      }
      for(j=0;j<nfiles;j++)  prob[j]/=probtotal;
      probtotal=0; for(j=0;j<nfiles;j++) probtotal+=prob[j];
    }
    ////////////////////////////////
    int m_LD=bayesopt2,nfiles2;//bayesopt2=m for LD= Loss_D = sum_{j=1}^m ((f_{l0+j}-f_{l0})^2 +(f_{l0-j}-f_{l0})^2)
    if(bayesmode==3) m_LD=0;
    else m_LD=bayesopt2;
    if(nfiles<m_LD*2+1) m_LD=(nfiles-1)/2;// (nfiles,m_LD)=(2,0),(3,1),(4,1),(5,2),(6,2);
    nfiles2=nfiles-m_LD*2;//for LD of generalization error
    int nmean=nfiles-nfiles2+1;//nmean=2*m_LD+1
    //    double *ww =(double*) malloc(sizeof(double)*nfiles*(nmean));
    //#define LDmode 2
    //#define LDmode 4
    //	//LDmode==4 is only for ib:0:0:0:m ???
    //#if LDmode == 2
    //    double *ww =(double*) malloc(sizeof(double)*(nfiles+1)*(nmean));
    //#elif  LDmode == 4
    //    double *ww0 =(double*) malloc(sizeof(double)*(nfiles+1)*(nmean)*2);
    //    double *ww =(double*) malloc(sizeof(double)*(nfiles+1)*(nmean)*(nmean));
    double *ww =(double*) malloc(sizeof(double)*(nfiles)*(nmean)*(nmean));
    //#endif
    //    double *ww=&ww0[1];
    double *ww0=&ww[(m_LD+m_LD*nmean)*nfiles];

    {//////////for Lvar0
      //      double *yp0=(double*) malloc(sizeof(double)*num0*nfiles);
      //      double *yt0=(double*) malloc(sizeof(double)*(nfiles+1)*num0+2);
      //      double *yt=(double*) &yt0[2];//yt[0]=num0,yt[1]=nmean;

      double *yt=(double*) malloc(sizeof(double)*(nfiles+1)*num0);
      double *yp01=(double*) &yt[num0];
      double *yp =(double*) malloc(sizeof(double)*num0);
      //      nfiles2=nfiles-2*bayesopt2;//for LD of generalization error
      if(nfiles<3) {m_LD=bayesopt2=0;nfiles2=nfiles;}//training ?
      //      else {m_LD=1;nfiles2=nfiles-2*bayesopt2;}//for LD of generalization error
      //      probtotal=0;
      //      for(j=m_LD;j<nfiles2+m_LD;j++){//prediction of training data?
      for(j=0;j<nfiles;j++){//prediction of training data?
#define NEWp0
#ifdef NEWp0
	{
	  char *p1,*p;
	  sprintf(fn_pred0,"%s",argv[j+j0]);
	  for(p=fn_pred0;;p++){if(*p=='+') break;}
	  p1=++p;
	  for(;;p++){if(*p=='+') break;}
	  ++p;
	  sprintf(p1,"%spred0.dat",p);
	}
#else
	sprintf(fn_pred0,"%spred0.dat",argv[j+j0]);
#endif
	fp0=fopen(fn_pred0,"r");
	for(i=0;i<num0;i++){
	  fgets(buff,buffsize,fp0);
	  sscanf(buff,"%lf%d%lf%lf%lf%lf%lf",&_yp,&_n,&_y,&_yr,&_var,&_cc,&_vv);//sscanf(buff,"%lf%d%lf%lf",&_yp,&_n,&_y,&_yr);
	  if(res.r3<0) _yp=_yr; //neglect resolution for the meanpred
	  else _yp=y_res(_yp,&res);
	  //	  yp0[i*nfiles+j] =_yp;
	  yp01[j*num0+i] =_yp;
	  yt[i]=_y;
	}
	fclose(fp0);
      }
      if(bayesmode==2){
	//	if(1==0 && nfiles==1) ww[0]=1;
	//	else
	{
	  for(j=0;j<nmean;j++){
	    for(l=0;l<nmean;l++){
	      //LDmode = 2,4,8 neighbours
	      if(LDmode==2 && l!=m_LD) continue;//only yoko
	      if(LDmode==4 && (l-m_LD)*(j-m_LD)!=0 ) continue;//only tate&yoko
	      //calc_AtinvbMCA((double *)(&ww[j*nfiles]),(double *)(&yp01[j*num0]),(double *)(&yt[0]),(int)(nfiles2),(int)num0,(int)bayesopt1);
	      calc_Atinvb((double *)(&ww[(l*nmean+j)*nfiles]),(double *)(&yp01[j*num0]),(double *)(&yt[0]),(int)(nfiles2+l-m_LD),(int)num0,(int)bayesopt1);
	      //	  calc_Ainvb((double *)(&ww[j*nfiles]),(double *)(&yp0[j]),(double *)(&yt[j]),(int)(nfiles2),(int)num0,(int)bayesopt1);
	      //	  calc_Ainvb((double *)(&ww[j*nfiles]),(double *)(&yp01[j*num0]),(double *)(&yt[0]),(int)(nfiles2),(int)num0,(int)bayesopt1);
	      //	    calc_Ainvb((double *)(&ww[j*nfiles]),(double *)(&yp01[j*num0]),(double *)(&yt[0]),(int)(nfiles2),(int)num0,(int)bayesopt1);
	      //	    for(i=0;i<=nfiles2;i++) ww[j*nfiles+i]=1./nfiles2;
	      //	    for(i=0;i<10;i++) fprintf(stderr,"j%2d,i%2d)yt%+e yp%+e %+e\n",j,i,yt[i],yp01[(j)*num0+i],yp01[(j+1)*num0+i]);
	      if(fabs(VarY)<1e-10){
		for(i=0;i<nfiles2+l-m_LD;i++) ww[(l*nmean+j)*nfiles+i]=1./(nfiles2+l-m_LD);
	      }
	      fprintf(stderr,"w[%d,%d]:",j,l);
	      for(i=0;i<nfiles2+l-m_LD;i++) fprintf(stderr,"%+.4f ",ww[(l*nmean+j)*nfiles+i]);
	      fprintf(stderr,"\n");
	    }
	  }
	}
      }

      else if(bayesmode==3){//ib:0:3:10:0:0 ib:0:3:itermax:0:0 H=C0*E2/E0+C1*Entropy
	for(j=0;j<nmean;j++){
	  for(l=0;l<nmean;l++){
	    if(LDmode==2 && l!=m_LD) continue;//only yoko
	    if(LDmode==4 && (l-m_LD)*(j-m_LD)!=0 ) continue;//only tate&yoko
	    {
	      size_t iter = 0;
	      size_t itermax=bayesopt1;if(itermax<1e-20) itermax=100;
	      int status;
	      double nl=(double)(nfiles2+l-m_LD);
	      double nd=(double)num0;
	      int nl1=nl+1;
	      double *par[7];
	      double sgn;
	      if(bayesopt2>=0) sgn=1; else sgn=-1;
	      par[0]=(double *)&nd;
	      par[1]=&nl;
	      par[2]=&yt[0];
	      par[3]=&yp01[j*num0];
	      par[4]=&sgn;
	      const gsl_multiroot_fsolver_type *T;
	      gsl_multiroot_fsolver *s;
	      gsl_multiroot_function df = {&my_df_Entropy, nl1, &par};
	      gsl_vector *x = gsl_vector_alloc (nl1); //&ww[(l*nmean+j)*nfiles];
	      for(i=0;i<=nl;i++) gsl_vector_set (x, i, sqrt(1./nl));
	      gsl_vector_set (x, nl, 0.0);//Good!?
	      T = gsl_multiroot_fsolver_hybrids;
	      s = gsl_multiroot_fsolver_alloc (T, nl1);
	      gsl_multiroot_fsolver_set (s, &df, x);
	      gsl_print_state (iter, s, par);
	      do {
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		gsl_print_state (iter, s, par);
		/* check if solver is stuck */
		if (status) {
		  printf ("Minimum not found (ERR=%d)? See /usr/include/gsl/gsl_errno.h for details.",status);
		  if (status==GSL_EBADFUNC) {printf ("Encountered a singular point.\n");}
		  else if (status==GSL_EBADFUNC) {printf ("Noprogress.\n");}
		  else if (status==GSL_ENOPROGJ) {printf (" ERR=%d:Jacobian is not improving the solution.\n",GSL_ENOPROGJ);}
		  break;}//break at Sucess or failure ?
		//		status = gsl_multiroot_test_residual (s->f, 1e-7);
		status = gsl_multiroot_test_residual (s->f, 1./nl);
	      }  while (status == GSL_CONTINUE && iter < itermax);
	      for(i=0;i<nl;i++) ww[(l*nmean+j)*nfiles+i]=square(gsl_vector_get (s->x, i));
	      gsl_multiroot_fsolver_free (s);
	      gsl_vector_free (x);
	      double sumw=0;
	      fprintf(stderr,"w[%d,%d]:",j,l);
	      for(i=0;i<nl;i++) fprintf(stderr,"%+.4f ",ww[(l*nmean+j)*nfiles+i]);
	      for(i=0;i<nl;i++) sumw+=ww[(l*nmean+j)*nfiles+i];
	      fprintf(stderr,"sumw=%f\n",sumw);
	    }
	  }
	}
      }
      else{//for bayesmode!=2
	for(l=0;l<nmean;l++){
	  if(LDmode==2 && l!=m_LD) continue;//only yoko
	  for(j=0;j<nmean;j++){
	    if(LDmode==4 && (l-m_LD)*(j-m_LD)!=0 ) continue;//only tate&yoko
	    for(i=0;i<nfiles2+l-m_LD;i++) ww[(l*nmean+j)*nfiles+i]=1./(nfiles2+l-m_LD);
	  }
	  //	  fprintf(stderr,"w[0-%d,0-%d]=%+.4f\n",nfiles2,nfiles2+l-m_LD,1./(nfiles2+l-m_LD));
	}
      }
      LtrainE=Ltrain=Lvar0=0;
      for(i=0;i<num0;i++){
	//      yp[i]=0;for(j=m_LD;j<nfiles2+m_LD;j++) yp[i]+=(yp01[j*num0+i]*ww[m_LD*nfiles+j-m_LD]);//prediction of training data
	yp[i]=0;
	for(j=0;j<nfiles2;j++){
	  yp[i]+=(yp01[(j+m_LD)*num0+i]*ww0[j]);//prediction of training data
	}
	//	for(j=m_LD;j<nfiles2+m_LD;j++){
	//	  //	  yp[i]+=(yp01[j*num0+i]*ww[(m_LD+m_LD*nmean)*nfiles+j-m_LD]);//prediction of training data
	//	  yp[i]+=(yp01[j*num0+i]*ww0[j-m_LD]);//prediction of training data
	//	}
	yp[i]=y_res(yp[i],&res);
	for(j=m_LD;j<nfiles2+m_LD;j++){
	  Lvar0+= square(yp[i]-yp01[j*num0+i]);
	  LtrainE += square(yp01[(j)*num0+i]-yt[i]);
	}
	//	for(j=m_LD;j<nfiles2+m_LD;j++) Lvar0+= exp(square(yp[i]-yp01[j*num0+i]));//???060530
	Ltrain+= square(yt[i]-yp[i]);
      }
      //      Lvar0 /=(nfiles2*num0); //      Lvar0 = VarY;//?????
      //Lvar0 /=(nfiles2*num0-1+1e-10); //unbiased estimator
      if(nfiles2>=2) Lvar0 /=((num0)*(nfiles2)*(nfiles2-1));//estimated variance of ensemble prediction
      else Lvar0 /=(num0); //average of unbiased estimator
      Ltrain/=(num0);
      LtrainE/=(nfiles*num0);
      //      free(yp01);//??why bat result if active?
#define Ltrain4Entropy
#ifdef Ltrain4Entropy
      {
	double pn,E2sum=Ltrain*num0;//Entropy
	H0=0;
	for(i=0;i<num0;i++){
	  pn=square(yp[i]-yt[i])/E2sum;
	  H0 -= pn*log(pn);
	}
	H0/=log(num0);
	//	Ltrain=H0;
      }
#endif
      free(yt);//
      free(yp);//
    }//end of "for Lvar0"
#define BAYES2GE 1
#if BAYES2GE == 1
    {//
      int li,l;
      yp00=(double*) malloc(sizeof(double)*num*nfiles);
      //#if LDmode == 2
      //      double *ypbuff=(double*) malloc(sizeof(double)*num*nmean);
      //#elif LDmode == 4
      //      double *ypbuff=(double*) malloc(sizeof(double)*num*nmean*2);
      double *ypbuff=(double*) malloc(sizeof(double)*num*nmean*(nmean));
      //#endif
      y=(double*) malloc(sizeof(double)*num);
      yr=(double*) malloc(sizeof(double)*num);
      n=(int*) malloc(sizeof(int)*num);
      var=(double*) malloc(sizeof(double)*num);
      vv=(double*) malloc(sizeof(double)*num);
      //      cc=(int*) malloc(sizeof(int)*num);
      cc=(int*) malloc(sizeof(int)*num*nfiles);
      kyorix2=(double*)malloc(sizeof(double)*num);
      //      yp=&ypbuff[num*(int)(bayesopt2)];//??what happen when bayesopt2==0??
      yp=&ypbuff[(m_LD+m_LD*nmean)*num];//center of nmean*nmean=(2*m_LD+1)*(2*m_LD+1)
      // yp=&ypbuff[m_LD*num];//
      for(i=0;i<num;i++){yp[i]= n[i]= y[i]= yr[i]=var[i]=vv[i]=0; }
      {///???
	for(j=0;j<nfiles;j++){
	  sprintf(fn_pred,"%spred.dat",argv[j+j0]);
	  fp=fopen(fn_pred,"r");
	  for(i=0;i<num;i++){
	    fgets(buff,buffsize,fp);
	    sscanf(buff,"%lf%d%lf%lf%lf%lf%lf%lf%lf",&_yp,&_n,&_y,&_yr,&_var,&_cc,&_vv,&_e2,&_kyorix2);
	    //	    yp0[i*nfiles+j] =y_res(_yp,&res);
	    if(res.r3<0) yp00[j*num+i]=_yr; //neglect resolution for the meanpred
	    else yp00[j*num+i] =y_res(_yp,&res);
	    //	    yp00[j*num+i] =y_res(_yp,&res);
	    cc[j*num+i] =_cc;//??
	    if(ccmax<_cc) ccmax=_cc;
	    if(j>=m_LD && j<nfiles-m_LD){
	      //	      yp[i]+=_yp*prob[j];//sum for j in [1,...,nfiles-2]
	      n[i]+=_n;
	      y[i]+=_y;//sum for j in [1,...,nfiles-2]
	      yr[i]+=_yr;
	      var[i]+=(_var*_vv);
	      vv[i]+=_vv;
	      //	      cc[i]=_cc;
	      kyorix2[i]=_kyorix2;
	    }
	  }
	  //	  if(j>=m_LD && j<nfiles-m_LD) probtotal+=prob[j];
	  fclose(fp);
	}
	//////////////
	Lvar=Lvarval=Lval=0;
	skew=kurtosis=entropy=mcov=skew2=0,skewa=0;
	m3p=m3m=m2=m3a=m4=m5a=m6=0;nm3p=nm3m=0;
	for(i=0;i<num;i++){//num=1?
	  //	  yp[i]/=probtotal;//mean for j in [1,...,nfiles-2]
	  //	  yp[i]=y_res(yp[i]/probtotal,&res);
	  yp[i]=0;
	  //	  for(j=0;j<nfiles2;j++) {yp[i] +=yp0[i*nfiles+j]*ww[(m_LD)*nfiles+j];}
	  //	  for(j=0;j<nfiles2;j++) {yp[i] +=yp00[(j+m_LD)*num+i]*ww[(m_LD)*nfiles+j];}
	  //	  for(j=m_LD;j<nfiles2+m_LD;j++) yp[i]+=(yp00[j*num+i]*ww[m_LD*nfiles+j]);//prediction of center of test data
	  //	  for(j=0;j<nfiles2;j++) yp[i]+=(yp00[(j+m_LD)*num+i]*ww[(m_LD+n_LD*nmean)*nfiles+j]);//prediction of center of test data
	  for(j=0;j<nfiles2;j++) yp[i]+=(yp00[(j+m_LD)*num+i]*ww0[j]);//prediction of center of test data
	  yp[i]=y_res(yp[i],&res);
	  n[i] = (int)((double)n[i]/nfiles2+0.5);
	  y[i] = y[i]/nfiles2;//mean for j in [1,...,nfiles-2]
	  yr[i]= yr[i]/nfiles2;
	  //	  for(j=m_LD;j<nfiles-m_LD;j++){
	  //	  for(j=m_LD;j<nfiles2+m_LD;j++){
	  double skewj=0,sigma2j=0;
	  double m2j=0,m3j=0,m4j=0,m5j=0,m6j=0;
	  double m3pj=0,m3mj=0;
	  //	  int nm3pj=0,nm3mj=0;
	  double nm3pj=0,nm3mj=0;
	  for(j=0;j<nfiles2;j++){//number of bags
	    double err=yp00[(j+m_LD)*num+i]-yp[i];//inversed 091008
	    double err2_=err*err;
	    m2j+=err2_;//sigma2j+=err2;//m2j
	    skew2+=(err>0)?(err2_):(-err2_);//?
	    m3j+=(err2_*=err);//skewj+=(errn=err2*err);//m3j
	    if(err>0){
	      nm3pj++;
	      m3pj+=err2_;//err^3
	    }
	    else if(err<0){
	      nm3mj++;
	      m3mj+=err2_;//err^3
	    }
	    m4j+=(err2_*=err);//kurtosis+=(errn=errn*err);//m4j
	    m5j+=(err2_*=err);
	    m6j+=(err2_*=err);
	    //Lvar    += square(yp[i]-yp00[(j+m_LD)*num+i]);
	    //Lvar    += square(yp[i]-yp00[(j+m_LD)*num+i]);
	    //Lvar    += exp(square(yp[i]-yp00[(j+m_LD)*num+i]));//???060530
	    Lvarval += square(yp00[(j+m_LD)*num+i]-y[i]);
	    {//mean covariance
	      int l;
	      for(l=0;l<nfiles2;l++){
		if(l!=j) mcov+=(err*(yp00[(l+m_LD)*num+i]-yp[i]));
	      }
	    }
	  }//closing for(j=0;j<nfiles2;j++){
	  Lval += square(yp[i]-y[i]);
	  Lvar+=m2j;//Lvar+=sigma2j;
	  double sj=sqrt(m2j/nfiles2);// double sj=sqrt(sigma2j/nfiles);
	  sj+=ZERO;
	  nm3pj+=ZERO;
	  nm3mj+=ZERO;
	  int mnorm=0;
	  if(mnorm==0){//original
	    m2+=(m2j/nfiles2);
	    m3+=(m3j/nfiles2)/pow(sj,3.);//skew+=skewj;//m3
	    m3a+=fabs(m3j/nfiles2)/pow(sj,3.);//skewa+=fabs(skewj/nfiles2)/pow(sj,3);//m3a
	    m3p+=(m3pj/(nm3pj))/pow(sj,3.);
	    m3m+=(m3mj/(nm3mj))/pow(sj,3.);
	    nm3p+=nm3pj;
	    nm3m+=nm3mj;
	    m4 +=    (m4j/nfiles2)/pow(sj,4.);//
	    m5a+=fabs(m5j/nfiles2)/pow(sj,5.);//?5/2
	    m6 +=    (m6j/nfiles2)/pow(sj,6.);//?3=6/2
	  }
	  else if(mnorm==1){//new for normalize m3,m4,m5,m6
	    m2+=(m2j/nfiles2);
	    m3+=pow(m3j/nfiles2,1./3.)/sj;//skew+=skewj;//m3
	    m3a+=pow(fabs(m3j/nfiles2),1./3.)/sj;//skewa+=fabs(skewj/nfiles2)/pow(sj,3);//m3a
	    m4 +=pow(m4j/nfiles2,1./4.)/sj;//
	    m5a+=pow(fabs(m5j/nfiles2),1./5.)/sj;//?5/2
	    m6 +=pow(m6j/nfiles2,1./6.)/sj;//?3=6/2
	  }
	  else if(mnorm==2){//this does not work for skew and kurtosis
	    m2+=(m2j/nfiles2);
	    m3+=pow(m3j/nfiles2,1./3.);//skew+=skewj;//m3
	    m3a+=pow(fabs(m3j/nfiles2),1./3.);//skewa+=fabs(skewj/nfiles2)/pow(sj,3);//m3a
	    m4 +=pow(m4j/nfiles2,1./4.);//
	    m5a+=pow(fabs(m5j/nfiles2),1./5.);//
	    m6 +=pow(m6j/nfiles2,1./6.);//
	  }
	  //	  if(i<10) fprintf(stderr,"[%d]Lval=%9.3e,yp=%.3e,y=%.3e\n",i,Lval,yp[i],y[i]);
	  //	  if(i>1900) fprintf(stderr,"[%d]Lval=%9.3e,yp=%.3e,y=%.3e\n",i,Lval/num,yp[i],y[i]);
	}//closing for(i=0;i<num;i++){
	int num2=num*nfiles2;
	//Lvar /= (nfiles2*num);
	//??	if(nfiles2>=2) Lvar /=((num)*(nfiles2)*(nfiles2-1.));//estimated variance of ensemble prediction
	//if(nfiles2>=2) Lvar /=num2;//estimated variance of ensemble prediction
	//else Lvar/=(num);
	for(i=0;i<num;i++){
	  for(j=0;j<nfiles2;j++){
	    double err=yp[i]-yp00[(j+m_LD)*num+i];
	    double p=(err*err)/(Lvar+ZERO);
	    p+=ZERO;
	    entropy -= p*log(p);
	  }
	}
	entropy/=log(num2+ZERO);
	//
	Lvar/=num2;//?see above
	mcov/=(num2*(nfiles2-1));
	double sigma=sqrt(Lvar);
	double sigma3=sigma*sigma*sigma;
	m2/=num;
	m3/=num;//skew=(skew/num2)/sigma3;//	skewa/=num;
	m3a/=num;//	skewa/=num;
	m3p/=num;
	m3m/=num;
	m4/=num;//kurtosis=(kurtosis/num2)/sigma3/sigma;//!modified 20130418
	m5a/=num;
	m6=m6/num;
	skew=m3;skewa=m3a;kurtosis=m4;
	Lvar=m2+ZERO;//??
	//	skew2=skew2/num2/Lvar;
	//	kurtosis=(kurtosis/num2)/sigma3/sigma-3;//!!!m4-3
	/////////////
	//Lvar /=((nfiles2*num)*(nfiles2*num-1.));//A-type Uncertainty
	//	Lvar /=((nfiles2*num)*(nfiles2*num-1.));//A-type Uncertainty
	Lvarval /= (num2);
	Lval /= (num);
	//	fprintf(stderr,"#Lvar=%e,Lvarval=%e,Lval=%e,_LD:",Lvar,Lvarval,Lval);
	//	fprintf(stderr,"#Lvar=%e,Lvarval=%e,Lval=%e,num=%d,_LD:",Lvar,Lvarval,Lval,num);
	//fprintf(stderr,"%g %g %g %g %g %g #Ltst,Lvar,skew,krt,ent,log|s*k|#n%db%d\n",Lval,Lvar,skew,kurtosis,entropy,log(fabs(skew*kurtosis)),num,nfiles);
	{
	  FILE *fp=fopen("loss.dat","w");
                                                   
	  fprintf(fp,"%g %g %g %g %g %d %d %g %g %g %d %g %g %g %g %d %g %d\n#Ltst,Lvar,skw,krt,ent,n,b,cor,skwa,yp0,num,skw2,m5a,m6,m3p,nm3p,m3m,nm3m\n",
		  Lval,Lvar,m3,m4,entropy,num,nfiles,mcov/Lvar,m3a,yp[0],num,skew2,m5a,m6,m3p,nm3p,m3m,nm3m);//4=K,9=sa,13=m5a,14=m6
	  //      1    2    3  4  5       6   7      8         9   10    11  12    13  14 15  16   17  18
	  //	  fprintf(fp,"%g %g %g %g %g %d %d %g %g %g %d %g %g %g\n#Ltst,Lvar,skw,krt,ent,n,b,cor,skwa,yp0,num,skw2\n",Lval,Lvar,skew,kurtosis,entropy,num,nfiles,mcov/Lvar,skewa,yp[0],num,skew2,m5a,m6);//4=K,9=sa,13=m5a,14=m6
	  //	  fprintf(fp,"%g %g %g %g %g %d %d %g %g %g %d\n#Ltst,Lvar,skw,krt,ent,n,b,cor,skew2,yp0,num\n",Lval,Lvar,skew,kurtosis,entropy,num,nfiles,mcov/Lvar,skew2,yp[0],num);
	  //	  fprintf(fp,"%.4g %.4g %.4g %.4g %.4g %d %d\n#Ltst,Lvar,skew,krt,ent,n,b\n",Lval,Lvar,skew,kurtosis,entropy,num,nfiles);
	  fclose(fp);
//	  fprintf(stderr,"###cat loss.dat###\n");
//	  system("cat loss.dat");
//	  fprintf(stderr,"#######\n");
	}
	//	fprintf(stderr,"#Lvar%e,Lvarval%e,Ltst%e,skew%.3g,krt%g,ent%.3g,n%d,b%d,",Lvar,Lvarval,Lval,skew,kurtosis,entropy,num,nfiles);
	int jj;
	for(i=0;i<num;i++){
	  for(j=-m_LD;j<=m_LD;j++){//l1
	    for(l=-m_LD;l<=m_LD;l++){//ld
	      if(j==0 && l==0 ) continue;//already obtained above
	      if(LDmode==2 && l!=0) continue;//only yoko
	      if(LDmode==4 && (l)*(j)!=0 ) continue;//only tate&yoko
	      li=(l*nmean+j)*num+i;
	      yp[li]=0;
	      //	    for(j=0;j<nfiles2;j++) {yp[li] +=yp0[i*nfiles+j]*ww[(m_LD+l)*nfiles+j];}
	      //	    for(j=0;j<nfiles2;j++) {yp[li] +=yp00[(j+m_LD+l)*num+i]*ww[(m_LD+l)*nfiles+j];}
	      for(jj=0;jj<nfiles2+l;jj++) {
		yp[li] +=yp00[(jj+m_LD+j)*num+i]*ww0[(l*nmean+j)*nfiles+jj];
		//		fprintf(stderr,"jj+m_LD+j=%d+%d+%d=%d\n",jj,m_LD,j,jj+m_LD+j);
	      }
	      yp[li]=y_res(yp[li],&res);
	      //	    li=(l+m_LD+1)*num+i;
	      //	    yp[li]=0;
	      //	    for(j=0;j<nfiles2+l;j++) {yp[li] +=yp00[(j+m_LD)*num+i]*ww[((l*nmean+m_LD)*nfiles+j];}
	      //	    yp[li]=y_res(yp[li],&res);
	    }
	  }
	}
	if(1==0){//unnecessary????
	  LD=0;
	  for(j=-m_LD;j<=m_LD;j++){//l1
	    for(l=-m_LD;l<=m_LD;l++){
	      if(j==0 && l==0 ) continue;//already obtained above
	      if(LDmode==2 && l!=0) continue;//only yoko
	      if(LDmode==4 && (l)*(j)!=0 ) continue;//only tate&yoko
	      _LD=0;
	      for(i=0;i<num;i++) {
		_LD += square(yp[i]-yp[(l*nmean+j)*num+i]);
		//		_LD += square(yp[i]-yp[l*num+i]);
		//		_LD += square(yp[i]-yp[(l+m_LD+1)*num+i]);
	      }
	      _LD/=num;
	      LD+=_LD;
	      fprintf(stderr,"_LD%9.3e ",_LD);
	    }
	  }
	  if(LDmode==2) LD/=(m_LD*4.);//
	  else if(LDmode==4) LD/=(m_LD*8.);//
	  else LD/=(nmean*nmean-1);//
	  //	  LD/=(m_LD*8.);//LD/=(m_LD*2.);
	  fprintf(stderr,"LD%e\n",LD);
	}//end of if(1==0)
      }//end of ???
    }//end of for BAYES2GE
#endif
  }//end of bayesmode<100

  //obtain heterosedastic err only vaid for known y[i] ?
#define newvari
#ifdef newvari //20150115
  double *bias=(double*) malloc(sizeof(double)*num);
  double *err=(double*) malloc(sizeof(double)*num);
  int    *nVi=(int*) malloc(sizeof(int)*num);
  {
    int ji,ii;
    for(i=0;i<num;i++) {
      bias[i]=var[i]=nVi[i]=0;
    }
    double *varcc=(double*)malloc(sizeof(double)*ccmax);
    double *biascc=(double*)malloc(sizeof(double)*ccmax);
    int *nvarcc=(int*)malloc(sizeof(int)*ccmax);
    //    free *varcc;
    for(j=0;j<nfiles;j++){//bags
      int cci;
      for(cci=0;cci<ccmax;cci++) varcc[cci]=biascc[cci]=nvarcc[cci]=0;
      for(i=0;i<num;i++){
	err[i]=yp00[j*num+i]-yp[i];
	varcc[cc[i]]+=square(err[i]);
	biascc[cci]=err[i];
	nvarcc[cc[i]]++;
      }
      for(i=0;i<num;i++){
	var[i]+=varcc[cc[i]];
	bias[i]+=biascc[cc[i]];
	nVi[i]+=nvarcc[cc[i]];
      }
    }
    for(i=0;i<num;i++){
      bias[i]/=nVi[i];
      var[i]/=nVi[i];
    }
    net_printf("##################new variance###\n");
  }
#else //old var[i]
  double *bias=(double*) malloc(sizeof(double)*num);
  double *err=(double*) malloc(sizeof(double)*num);
  int    *nVi=(int*) malloc(sizeof(int)*num);
  int    *nVib=(int*) malloc(sizeof(int)*num);
  {
    int ji,ii;
    for(i=0;i<num;i++) {
      err[i]=yp[i]-y[i];
      bias[i]=var[i]=nViv[i]=nVib[i]=0;
    }
    for(j=0;j<nfiles;j++){//bags
      for(i=0;i<num;i++){
	bias[i]+=err[i];        nVib[i]++;
	var[i] +=square(err[i]);nVi[i]++;
	//	if(1==0){//for check Lval0=9.126e-02
	{//for obtaining the average of the noise
	  ji=j*num+i;
	  for(ii=i+1;ii<num;ii++){
	    if(cc[ji]==cc[j*num+ii]) {//error for the same cell 
	      //bias[i]+=err[ii];  nVib[i]++;
	      //bias[ii]+=err[i];  nVib[ii]++;
	      var[i] +=square(err[ii]); nVi[i]++;
	      var[ii]+=square(err[i]);  nVi[ii]++;
	    }
	  }
	}
      }
    }
    for(i=0;i<num;i++){
      bias[i]/=nVib[i];
      var[i]/=nViv[i];
    }
  }
#endif
  //////
  //  fp=fopen(fnpredict,"w");
  fp=fopen("predict+.dat","w");
  for(i=0;i<num;i++){
    if(intmode==1){
      //      fprintf(fp,"%d %d %.7e %.7e %.7e %d %.7e %.7e %.7e #Y^,t,Y,yr,var,c,v,k,b\n",(int)(yp[i]+0.5),n[i],y[i],yr[i],var[i]/vv[i],cc[i],vv[i],kyorix2[i],bias[i]);
      fprintf(fp,"%d %d %.7e %.7e %.7e %d %.7e %.7e %.7e %d #Y^,t,Y,yr,var,c,v,k,b,nvar\n",(int)(yp[i]+0.5),n[i],y[i],yr[i],var[i],cc[i],vv[i],kyorix2[i],bias[i],nVi[i]);
    }
    else{
      //      fprintf(fp,"%.7e %d %.7e %.7e %.7e %d %.7e %.7e %.7e #Y^,t,Y,yr,var,c,v,k,b\n",yp[i],n[i],y[i],yr[i],var[i]/vv[i],cc[i],vv[i],kyorix2[i],bias[i]);
      fprintf(fp,"%.7e %d %.7e %.7e %.7e %d %.7e %.7e %.7e %d #Y^,t,Y,yr,var,c,v,k,b,nvar\n",yp[i],n[i],y[i],yr[i],var[i],cc[i],vv[i],kyorix2[i],bias[i],nVi[i]);
    }
  }
  fclose(fp);
  //  fp=fopen(fnloss,"w");
  fp=fopen("loss+.dat","w");
  fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e %.7e #LD Lval Lvar Lvarval Lvar0 nfiles num Ltrain LtrainE H0\n",LD,Lval,Lvar,Lvarval,Lvar0,nfiles,num,Ltrain,LtrainE,H0);
  fclose(fp);
  //  fprintf(stderr,"Created '%s'=%s",fnpredict,fnbody(argv[1],pr));
  //  for(j=2;j<=nfiles;j++) fprintf(stderr,"+%s" ,fnbody(argv[j],pr));
  //  fprintf(stderr,"\n");
  net_printf(strx3("Created '%s' and '%s'. #LDmode=%d.\n","predict+.dat","loss+.dat",LDmode));
  return(1);
}
