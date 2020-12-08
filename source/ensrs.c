/* Resampling Methods
   Leave-Many-Out cross-validation (K=1 for HoldOut, K>1 for Monte-Carlo)
   Multifold cross-validation (K=n for Leave-one-out, K<n for K-fold (V-fold))
   Bootstrap 
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include "my_misc.c"
////////////////
///
//#define RAND ZMTRAND
//#define RAND DRAND48
//#define RAND RONDOM
//#define RAND MYROND

#include "randoms.c"
//#if RAND == ZMTRAND
//#include "share/zmtrand.c"
//#elif RAND == MYRAND
//#include "random.c"
//#endif

#define DIRECT       -1
#define LEAVEMANYOUT 0
#define MULTIFOLD    1
#define BOOTSTRAP    2
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

char *methodname1[4]={"Direct","LeaveManyOut","MultiFold","BootStrap"};

#undef isnump
#define isnump(c) ((c >= '0' && c <= '9')||c=='.'||c=='-'||c=='+'||c==' '||c=='e')
#define nFoldsmax 100000
//#define buffsize 20000
#define buffsize 10240
//#define Kmax 1200

int NC=6;       //number of comparing cells
int K=2;        //dimension
double vt=0.2;  //v_threshold
int vm=3;       //   # number of compensation vectors
int vm2=0;      //  # number of minimum vectors to be neglected
double vr=5;  // # v_ratio
double w=0.2;   // # width
double ymin=0,ymax=0,ymin1=0,ymax1=0;
double *xmin,*xmax,*xmin1,*xmax1;
int r1=0,r2=0;
int r0=0;
double r3=0;
double gamma0=0.05,entropy_thresh=0.7;
double alpha;
double BIAS=1;
int LINESIZE=10240;
void usage(int argc,char **argv)
{
  fprintf(stderr,"Ensemble & Resampling Program for CAN2.\n");
  //  fprintf(stderr,"Usage: %s <fn_given> <m>:<nFolds | fn_valid>[:<nValidData>] <<N1>:<Nens>:<NStep>|<N1>-<N2>:<NStep>> [<options>]\n",argv[0]);
  fprintf(stderr,"Usage: %s <fn_given> <method> <nCells> [<options>]\n",argv[0]);
  fprintf(stderr,
	  "<fn_given>        = file name of training data\n"
          "  <method of Validation>;\n"
	  "  0:<nFolds>:<alpha>[:<seed>]= <nFolds>-fold LeaveManyOut with  and nTrainData=<alpha>*nGivenData.\n"
	  "  1:<nFolds>                = <nFolds>-fold CrossValidation.\n"
	  "  2:<nFolds>:<alpha>[:<seed>]= <nFolds>-fold BootStrap with nTrainData=<alpha>*nGivenData.\n"
	  " -1:<fn_validdata>  = file name of validation data for DirectValidation. \n"
	  "<N1>[:<NN>[:<ND>]]  = number of cells N=N1+ND*i for i=0,1,..,NN-1 .\n"
	  "<N1>-<N2>[:<ND>]]   = number of cells N=N1+ND*i for i=0,1,... until N<=N2.\n"
	  "lc:<last command>  = last command placed in the commands for the CAN2.\n"
	  "options; \n"
	  "N:<N1>:<N2>...    = number of cells overwrite above.\n"
	  "k:<k>             = dimension of input vector[2]\n"
	  "T:<T>             = number of learning iterations [100:default]  \n"
	  "g:<gamma>         = gamma[0.05] \n"
	  "e:<entropy_thresh>= entropy_thresh[0.7]  \n"
	  "y:<ymin>:<ymax>:<ymin1>:<ymax1>    = transform [ymin,ymax] to [ymin1,ymax1] [0:0:0:0]\n"
	  "x<i>:<xmin>:<xmax>:<xmin1:<xmax1> = transform [x_i_min,x_i_min] to [x_i_min1,x_i_max1] \n"
	  "r:<r1>:<r2>:<r3>  = output resolution r1/r2  [0:0]\n"
	  "vm:<vm>           = vmin [3]\n"
	  "vm2:<vm2>         = vmin2 [0]  \n"
	  "NC:<NC>           = number of comparing units [6]\n"
	  "vt:<vt>           = v_thresh	    [0.2]\n"
	  "vr:<vr>           = v_ratio      [5]   \n"
	  "w:<w>             = width        [0.2]\n" 
	  "Tsk:<Task>[:<t1>:<t2>:<t0>] <Task>==1 for regression, ==0 for time-series [t0:t1+t0-1] for training;[t1+t0:t2+t0] for test.\n" 
	  //	  "E:<err>           = 0 for MSE, 1 for IREM, 2 for FULLIREM, 3 for IREM, 4 for FULLIREM\n" 
	  "i:<imode>         = imode=1 for 'meanpred int', imode=0  for 'meanpred'\n" 
	  "rm:<rm>           =<rm>=1 or 2 for range data\n" 
	  "tau:<tau_c>:<tau_h>:eta1 = tau_c for weighting biases, tau_h for weighting boosting iterations.\n" 
	  );
}
static char pr[buffsize];//for fnbody
static char pr1[buffsize];//for fnbody
double *v_walkeralias;
int *a_walkeralias;
void init_walkeralias(int n,double *p)
{
  int k,j,i;
  int *S=(int *)malloc(sizeof(int)*n);
  int *G=(int *)malloc(sizeof(int)*n);
  v_walkeralias=(double *)malloc(sizeof(double)*n);
  a_walkeralias=(int *)malloc(sizeof(int)*n);
  i=j=-1;
  for(k=0;k<n;k++){
    v_walkeralias[k]=n*p[k];
    if(v_walkeralias[k]>=1) G[++i]=k;//
    else S[++j]=k;//
  }
  for(;j>=0&&i>=0;){   //  for(;j>=0;){
    a_walkeralias[S[j]]=G[i];
    v_walkeralias[G[i]]-=(1.-v_walkeralias[S[j]]);
    if(v_walkeralias[G[i]]<1){
      S[j]=G[i];
      --i;//if(--i<0)	break;
    }
    else --j;
    //    if(j==0){//check
    //      fprintf(stderr,"Owarimae?G[%d]=%d,S[%d]=%d\n",i,G[i],j,S[j]);
    //    }
    //    fprintf(stderr,"G[%d]=%d,S[%d]=%d\n",i,G[i],j,S[j]);
  }
  //  if(i>=0 && j>=0){
  //    fprintf(stderr,"?? strainge !? Wakler alias (i=%d,j=%d)",i,j);
  //  }
  free(S);
  free(G);
}
int walkeralias(int n){
  //int walkeralias(int n,double *p){
//#if RAND == ZMTRAND
//  //  double V=n*NextUnif();
//#elif RAND == MYRAND
#if RAND == MYRAND
#define NextUnif() myrandom()
  //  double V=n*myrandom();
#elif RAND == DRAND48
#define NextUnif() drand48()
  //  double V=n*drand48();
  //  double V=n*lrand48()/(RAND_MAX+1.0);
#elif RAND == RANDOM
#define NextUnif() (random()/(RAND_MAX+1.0))
  //  double V=n*(random()/(RAND_MAX+1.0));
#endif
  double V=n*NextUnif();
  int k=(int)V;
  double u=V-k;
  if(u<=v_walkeralias[k]) return(k);
  else return(a_walkeralias[k]);
}
double lossboost(double c, void *params){
  double **dp=(double **)params;
  int n=(int)(*dp[0]);
  double *wbst=dp[1];
  double *err2=dp[2];
  double f=0;
  double c2=1./sqrt(c);
  int i;
  for(i=0;i<n;i++) f+=wbst[i]*c2*exp(c*err2[i]);
  //  for(i=0;i<n;i++) f+=c2*exp(c*err2[i]);//for check
  return(f);
}
////for make fntrain fntest
int nTrainData;
int nGivenData;//,datanum;
int *nTestData;
int nFolds;
int nFolds1;
int method;
int *nr;//resample
int *nc;
double *yp0,*yt;//for BAGGING==2
#define NoBAG 0
#define BAGwithVal 1
#define BAGwithoutVal 2
#define NoBoost 0
#define EmBoost 1
#define GbBoost 2
int BAGGING=NoBAG;
int Boost=NoBoost;
double meannTestData=0;
int nValData=0;
int t_boost=0;//apply boosting for t_boost>=1, Gradient-based boosting for t_boost==-2
int chkomit=0;//int chkomit=1;
#define fnsize 256
char **fntrain=NULL;//fntrain[nFoldsmax][256];
char **fntest=NULL;// char fntest[nFoldsmax][256];
double err4propagate=0;
double err4terminate=0;
int nop=0;
int setfntraintest(unsigned int seed,char *fn_givendata,char **argv,char *fupdate)
//int setfntraintest(char **fntrain,char **fntest,unsigned int seed,char *fn_givendata,char **argv)
{
  int i;
  fntrain=(char**) malloc(sizeof(char*)*(nFolds));
  for(i=0;i<nFolds;i++) fntrain[i]=(char*) malloc(sizeof(char)*fnsize);
  fntest=(char**) malloc(sizeof(char*)*(nFolds));
  for(i=0;i<nFolds;i++) fntest[i]=(char*) malloc(sizeof(char)*fnsize);

  if(method==-1){
    nTrainData=nGivenData;
    sprintf(fntrain[0],"%s",fn_givendata);
    sprintf(fntest[0],"%s",&argv[2][3]);
  }
  else{//method!=-1
    int j;
    for(j=0;j<nFolds;j++){
      sprintf(fntrain[j],"./tmp/train%03d.dat",j);
      sprintf(fntest[j],"./tmp/test%03d.dat",j);
    }
    if(*fupdate=='0') return(0);
    //    if(method==BOOTSTRAP){
    {
      //bootstrap or resampling
#if RAND == ZMTRAND
      InitMt((unsigned long)seed);
#elif RAND == MYRAND
      //  _randn=1;
      //      _randn=0;
      _randn=seed;
      //  srandom(1);
      //srandom(10);
#elif RAND == DRAND48
      srand48((long int)seed);
#elif RAND == RANDOM
      srandom(seed);
#endif
    }
    //    if(opendir("./tmp/")==NULL) system("mkdir ./tmp");
    //    system("mkdir ./tmp");
    {//////////////////////////obtain training data, and separate them for CV
      FILE *fp1;
      if((fp1=fopen(fn_givendata,"r"))==NULL){
	fprintf(stderr,"file(%s) open error\n",fn_givendata);
	return(-1);
      }
      //      int j;
      FILE *fptrain_xy;
      FILE *fptest_xy;
      //      int nn;
      for(j=0;j<nFolds;j++){
	//	sprintf(fntrain[j],"./tmp/%strain%03d.dat",fnbody_givendata,j);
	//	sprintf(fntrain[j],"./tmp/train%03d.dat",j);
	fptrain_xy=fopen(fntrain[j],"w+");
	fclose(fptrain_xy);
	//	sprintf(fntest[j],"./tmp/%stest%03d.dat",fnbody_givendata,j);
	//	sprintf(fntest[j],"./tmp/test%03d.dat",j);
	fptest_xy=fopen(fntest[j],"w+");
	fclose(fptest_xy);
      }
      {
	//	int nn,*nr,*nc,*nTestData;//resample
	//	int nn,*nr,*nTestData;//resample
	//	int nn,*nr;//resample
	int nn,n;//resample
	if(method==LEAVEMANYOUT){
	  //	  nTrainData=nGivenData-alpha;
	  nTrainData=nGivenData*alpha+0.5;
	}
	else if(method==MULTIFOLD){
	  nTrainData=nGivenData*(nFolds1-1.)/nFolds1+0.5;
	  //	  if(nTrainData>=nGivenData) nTrainData=nGivenData-1;//050818??
	  //	  else if(nTrainData*(nFolds1)/(nFolds1-1.) <nGivenData) nTrainData++;
	  if(nTrainData*(nFolds1)/(nFolds1-1.) <nGivenData) nTrainData++;
	  //	  nFolds1=nFolds;
	  //	  nFolds=nFolds1*alpha;
	}
	else if(method==BOOTSTRAP){
	  nTrainData=alpha*nGivenData+0.5;
	  //nTrainData=nGivenData;//
	}
	nr=(int*) malloc(sizeof(int)*nTrainData*nFolds);
	nc=(int*) malloc(sizeof(int)*nGivenData*nFolds);
	for(j=0;j<nFolds;j++){
	  nTestData[j]=0;
	  for(nn=0;nn<nGivenData;nn++) nc[j*nGivenData+nn]=-1;
	  for(n=0;n<nTrainData;n++) nr[j*nTrainData+n]=-1;
	  if(method==BOOTSTRAP){
	    for(n=0;n<nTrainData;n++){
	      //	      if(t_boost<0){//no boosting
	      //	      if(Boost==NoBoost){//no boosting _061117
	      if(Boost==NoBoost || t_boost==1){//no boosting 061117
#if RAND == ZMTRAND
	      nn=NextIntEx(nGivenData);
#elif RAND == MYRAND
	      nn=((int)(myrandom()*nGivenData))%nGivenData;
	      //	    nn=nGivenData*myrandom();	    nn%=nGivenData;
#elif RAND == DRAND48
	      nn=(int)((double)nGivenData*lrand48()/(RAND_MAX+1.0));
	      //	      nn=(int)((double)nGivenData*drand48());
#elif RAND == RANDOM
	      nn=(int)((double)nGivenData*random()/(RAND_MAX+1.0));
#endif
	      }
	      else{//apply boosting if t_boost>=0
		nn=walkeralias(nGivenData);//for bootstrap resampling with pbst
		//		nn=walkeralias(nGivenData,pbst);//for bootstrap resampling with pbst
	      }

	      nr[j*nTrainData+n]=nn;//resampled n-th data is nn-th data
	      nc[j*nGivenData+nn]=1;//nn-th data is in training set
	    }
	  }
	  else if(method==LEAVEMANYOUT){
	    for(n=0;n<nTrainData;){
#if RAND == ZMTRAND
	      nn=NextIntEx(nGivenData);
#elif RAND == MYRAND
	      nn=((int)(myrandom()*nGivenData))%nGivenData;
#elif RAND ==DRAND48
	      nn=(int)((double)nGivenData*lrand48()/(RAND_MAX+1.0));
	      //	      nn=(int)((double)nGivenData*drand48());
#elif RAND == RANDOM
	      nn=(int)((double)nGivenData*random()/(RAND_MAX+1.0));
#endif
	      if(nc[j*nGivenData+nn]==-1){
		nr[j*nTrainData+n]=nn;//resampled n-th data is nn-th givendata
		nc[j*nGivenData+nn]=1;//nn-th data is in training set
		n++;
	      }
	    }
	  }
	  else if(method==MULTIFOLD){//K-fold CV
	    int *nn1,*nn1c;
	    nn1=(int*)malloc(sizeof(int)*nGivenData);
	    nn1c=(int*)malloc(sizeof(int)*nGivenData);
	    if(j<nFolds1) for(n=0;n<nGivenData;n++) nn1[n]=n;//for j<nFolds1
	    else{
	      for(n=0;n<nGivenData;n++) nn1c[n]=-1;
	      for(n=0;;){
#if RAND == ZMTRAND
	      nn=NextIntEx(nGivenData);
#elif RAND == MYRAND
		nn=((int)(myrandom()*nGivenData))%nGivenData;
#elif RAND == DRAND48
	      nn=(int)((double)nGivenData*lrand48()/(RAND_MAX+1.0));
	      //	      nn=(int)((double)nGivenData*drand48());
#elif RAND == RANDOM
		nn=(int)((double)nGivenData*random()/(RAND_MAX+1.0));
#endif
		if(nn1c[nn]<0){
//		  if(n==81 || n==165 || nn==165){
//		    fprintf(stdout,"n:%d nn:%d j:%d nn1c[%d]=%d.strange????\n",n,nn,j,nn,nn1c[nn]);
//		  }
		  nn1[n]=nn;
		  nn1c[nn]=1;
		  if(++n>=nGivenData) break;
		}
	      }
	    }
	    n=0;
	    for(nn=0;nn<nGivenData;nn++){
	      //		if((nn%nFolds1)!=j%(int)alpha){
	      if((nn%nFolds1)!=j%(int)nFolds1){
//		if(nn1[nn]==165){
//		  fprintf(stdout,"nn:%d j:%d nn1:%d strange????\n",nn,j,nn1[nn]);
//		}
		if(nc[j*nGivenData+nn1[nn]]>0) fprintf(stdout,"nn:%d j:%d nn1:%d strange\n",nn,j,nn1[nn]);
		nr[j*nTrainData+n]=nn1[nn];//resampled n-th data is nn-th givendata
		nc[j*nGivenData+nn1[nn]]=1;//nn-th data is in training set
		n++;
	      }
	    }
	    //	    fprintf(stdout,"n=%d.\n",n);fflush(stdout);
	    free(nn1);
	    free(nn1c);
	  }
//	  if(j==3){
//	    for(n=0;n<nTrainData;n++) fprintf(stdout,"%d %d #n nr\n",n,nr[j*nTrainData+n]);
//	    for(n=0;n<nGivenData;n++) fprintf(stdout,"%d %d #n nc\n",n,nc[j*nGivenData+n]);
//	  }
	}
	for(n=0;n<nGivenData;n++){
	  char *line1=(char*)malloc(sizeof(char)*LINESIZE);//char line1[buffsize];
	  line1[0]=0;
	  fgets(line1,LINESIZE,fp1);
	  if(feof(fp1)) break;
	  if(line1[0]==0) break;
	  
	  for(j=0;j<nFolds;j++){
	    if(nc[j*nGivenData+n]!=1){//n-th data is not in training set
	      nTestData[j]++;
	      fptest_xy=fopen(fntest[j],"a+");
	      fprintf(fptest_xy,"%s",line1); 
	      fclose(fptest_xy);
	    }
	    else{
	      for(nn=0;nn<nTrainData;nn++){
		if(nr[j*nTrainData+nn]==n){
		  fptrain_xy=fopen(fntrain[j],"a+");
		  fprintf(fptrain_xy,"%s",line1); 
		  //		  fprintf(stderr,"[%3d,%3d]%s,nTrainData=%d",n,nn,line1,nTrainData); 
		  fclose(fptrain_xy);
		}
	      }
	    }
	  }
	}
	{ 
	  for(j=0;j<nFolds;j++) {
	    //	    fprintf(stderr,"nTestData[%d]=%d.",j,nTestData[j]);
	    meannTestData+=nTestData[j];
	  }
	  meannTestData/=nFolds;
	  if(nop==0) printf("#Mean nTestData=%g in nGivenData=%d.\n",meannTestData,nGivenData);
	}
	//	free(nr);
	//	free(nc);
	//	free(nTestData);
      }
      fclose(fp1);
    }
  }//closing else{
  //the above are for preparing training datasets; closeing of else{//method!=-1
  return(0);
}//endof setfntraintest
int main(int argc,char **argv)
{
  FILE *fp1,*fp0;//,*fp2;
  //  FILE *fptrain_xy;
  //  FILE *fptest_xy;
  //  printf("start ensrs");
  int i;
  for(i=4;i<argc;i++){
    if(strncmp(argv[i],"LINESIZE:",9)==0) sscanf(&argv[i][9],"%d",&LINESIZE);
  }
  char *line1=(char*)malloc(sizeof(char)*LINESIZE);//char line1[buffsize];
  char *line2=(char*)malloc(sizeof(char)*LINESIZE);//char line1[buffsize];
  //  char line1[buffsize],line2[buffsize];
  //  char p;
  //  int nFolds;
  int nCells,_nCells;
  int nens; //nens=num is number of ensemble
  int T=100;
  //  int i;//, j;//, k;
  //  int nV;
  //  double *y,*yp,err,err2,*err2p;
  double err;//yj,ypj,err,err2;
  //  double *x1,*x2,*y1;
  double msesum=0,msemean,msestdp=0,msestdm=0;
  int n_msesum=0,n_msestdp=0,n_msestdm=0;
  double *mse;//mse[nFoldsmax];
  double *Lhat,Lhatmean=0;
  int *n_mse;//n_mse[nFoldsmax];
  double msetrainsum=0,*msetrain,msetrainmean;//
  int n_msetrainsum=0,*n_msetrain;//
  int *n_cells2,n_cells2sum=0;
  int *is;//is[nFoldsmax];
  //#define fnsize 256
  char fn_givendata[fnsize];
  char fnbody_givendata[fnsize];
  char fn_target[fnsize];
  char fn_pred[fnsize];
  char fn_pred0[fnsize];
  char fn_is[fnsize];
  char fn_enspred[fnsize];
  char fn_enspred0[fnsize];
  char fn_ensis[fnsize];
  char fb[fnsize];
#define NEWp0
#ifdef NEWp0
  char fb0[fnsize];
#endif
  char fn_boost[fnsize];
  int n;
  FILE *fp;
  double *sse;
  int j;//,n1;//  int i1,i2,j,n1;
  //  int method;//,nV;
  int intmode=0;
  char ibmode0[20];
  char *ibmode=ibmode0;
  int NN[1000];
  //  int num;
  char *lcom=NULL;
  //  unsigned int seed=1;
  unsigned int seed=0;
  int DISP=0;
  int Tpinv=999;
  //  int nop=0;
#define TIMESERIES 0
#define REGRESSION 1
  //  int Task=1,t1=0,t2=0,t0=0;
  int Task=1,tr0=0,tr1=0,tp0=0,tp1=0;
#define MSE 0
#define IREM 1
#define FULLIREM 2
#define IREM2 3
#define FULLIREM2 4
  //  int ERR=MSE,p0=0;
  //  int ERR=MSE;
  //double gamma=0.04,entropy_thresh=0.7;
//#define cmdsize 7000
//#define cmd1size 30000
//#define cmd2size 60000
  //  char cmd[cmdsize];
  int cmdsize=5000;
  int cmd1size=30000;
  int cmd2size=30000;
  char *cmd=(char*)malloc(sizeof(char)*cmdsize);
  char *cmd1=(char*)malloc(sizeof(char)*cmd1size);
  char *cmd2=(char*)malloc(sizeof(char)*cmd2size);
  //  char cmd1[cmd1size];
  //  char cmd2[32000];//for bagging;
  //  char cmd2[50000];//safe for b=500 increase for big number of baggs;
  //  char cmd2[cmd2size];//increase for big number of baggs;
  char rdir[buffsize];
  char *rdir0="result-ensrs2ge";
  DIR *pdir;
  int LDmode=2;
  int NStep;
  int nob_thresh=0;
  char *fupdate="1";
  char *pupdate="1";
  sprintf(ibmode,"%s","ib:0:0");
  if(argc<4){
    usage(argc,argv);
    return(-1);
  }

  sprintf(fn_givendata,"%s",argv[1]);
  {
    int ncol=0;
    char *p=argv[2];
    //  sscanf(argv[2],"%d:%d:%lf",&method,&nFolds,&alpha);
    //    sscanf(argv[2],"%d:%d:%lf:%d",&method,&nFolds,&alpha,&seed);
    for(;*p!=0;p++) if(*p==':') ncol++;
    if(ncol==1) sscanf(argv[2],"%d",&method);//second element is fntest[0] 
    else if(ncol==2) sscanf(argv[2],"%d:%d:%lf",&method,&nFolds,&alpha);
    else if(ncol==3) sscanf(argv[2],"%d:%d:%lf:%d",&method,&nFolds,&alpha,&seed);
    else sscanf(argv[2],"%d",&method);//second element is fntest[0] 
  }
  if(method==-1){
    nFolds=1;
  }
  else if(nFolds<1 || (method==0 && alpha<1e-10)){
    usage(argc,argv);
    return(-1);
  }

  if(method==MULTIFOLD){
    if(alpha<1) alpha=1;//?
    nFolds1=nFolds;
    nFolds=nFolds1*alpha;
  }


  {
    char *p;
    int N1,N2,ptn=0;
    p=argv[3];// arg[3]=<N1>:<nens>:<NStep> or <N1>-<N2>:<NStep> 
    sscanf(p,"%d",&N1); N2=N1; nens=1; NStep=1;
    for(;;){
      if(*p==0) break;
      if(ptn==0){//unknown pattern
	if(*p==':') {
	  sscanf(++p,"%d",&nens);N2=N1+nens-1; 
	  ptn=1;
	}
	if(*p=='-') {
	  sscanf(++p,"%d",&N2); 
	  ptn=2;
	}
      }
      else if(ptn>=1){
	if(*p==':'){
	  sscanf(++p,"%d",&NStep);
	  if(NStep<1) NStep=1;
	  if(ptn==1){N2=N1+nens*NStep-1;}
	  break;
	}
      }
      p++;
    }
    if(N1<1){              //nCellsのエラー処理
      fprintf(stderr,"nCells must be bigger than %d.\n",0);
      N1=N2=NStep=1;
    }
    nens=0;
    for(nCells=N1;nCells<=N2;nCells+=NStep){
      NN[nens++]=nCells;
    }
  }

//  if(1==0){//orig
//    char *p;
//    int ncolon=0,N1,N2,ND;
//    for(p=argv[3];;){
//      if(*p==0) break;
//      if(*p==':') {ncolon++;break;}
//      p++;
//    }
//    if(ncolon==0){
//      sscanf(argv[3],"%d",&N1);
//      N2=N1;ND=1;
//    }
//    else if(ncolon==1){
//      sscanf(argv[3],"%d:%d",&N1,&N2);
//      ND=1;
//    }
//    else if(ncolon==2) sscanf(argv[3],"%d:%d:%d",&N1,&N2,&ND);
//
//    if(N1<1){              //nCellsのエラー処理
//      fprintf(stderr,"nCells must be bigger than %d.\n",0);
//      N1=N2=ND=1;
//    }
//    nens=0;
//    for(nCells=N1;nCells<=N2;nCells+=ND){
//      NN[nens++]=nCells;
//    }
//  }
  //  char fn_bagging[fnsize];
  char *fn_bagging;
  int ibg=0;
  double rethresh_boost=12.;
  int rangedatamode=0;
  //double *cbst,*ypt,*wbst,*wbst1,tau_c=4.,tau_h=4;
  double *cbst,*ypt,*wbst,*wbst1,tau_c=8.,tau_h=8.,eta1=2.0;
  double Lstd=0;int Lstdm=2;//  double Lstd=0;int Lstdm=0;
  char *lossall=NULL;
  int Bayes=0,BayesSeed=0,iBayes=0;
  double BayesLambdaL=0,BayesLambdaS=0; int BayesUseAllData=0, Bayesseed=1;

  for(i=4;i<argc;i++) if(strncmp(argv[i],"k:",2)==0) sscanf(&argv[i][2],"%d",&K);
  xmin=(double *)malloc(sizeof(double)*K);
  xmax=(double *)malloc(sizeof(double)*K);
  xmin1=(double *)malloc(sizeof(double)*K);
  xmax1=(double *)malloc(sizeof(double)*K);
  for(j=0;j<K;j++){xmin[j]=xmin1[j]=0;xmax[j]=xmax1[j]=0;}
  double *x0m=(double *)malloc(sizeof(double)*K);
  double *x1m=(double *)malloc(sizeof(double)*K);
  double *x0M=(double *)malloc(sizeof(double)*K);
  double *x1M=(double *)malloc(sizeof(double)*K);
  double y0m=0.,y1m=0.,y0M=0,y1M=0.;
  int ssp=0;
  for(j=0;j<K;j++){x0m[j]=x0M[j]=0;x1m[j]=x1M[j]=0;}
  //  for(i=0;i<K;i++) fprintf(stderr,"%e %e %e %e   #%d ???入力の正規化 x0min0 x0max0 →x0min x0max \n",xmin[i],xmax[i],xmin1[i],xmax1[i],i);
#ifdef CHKTIME
  mytimer_start();
#endif
  for(i=4;i<argc;i++){
    double _xm,_xM,_xm1,_xM1;
    //    double _xm,_xM;
    if(strncmp(argv[i],"vm2:",4)==0) sscanf(&argv[i][4],"%d",&vm2);
    else if(strncmp(argv[i],"chkomit:",8)==0) sscanf(&argv[i][8],"%d",&chkomit);
    else if(strncmp(argv[i],"bg:",3)==0){
      //      sprintf(fn_bagging,"%s",&argv[i][3]);
      fn_bagging=&argv[i][3];
      ibg=i;
      if(strncmp(fn_bagging,"/dev/null",9)==0) BAGGING=BAGwithoutVal;//bagging without validfile
      else if((fp=fopen(fn_bagging,"r"))!=NULL){
	BAGGING=BAGwithVal;//bagging with validfile
	nValData=0;
	for(;;){
	  fgets(line1,LINESIZE,fp);
	  if(feof(fp)) break;
	  if(line1[0]==0) break;
	  //      printf("nGivenData=%d,line=%s",nGivenData,line1);
	  nValData++;
	}
	fclose(fp);
      }
      else{
	fprintf(stderr,"There is no bagging file (%s)\n",fn_bagging);
	return(-1);
      }
    }
    else if(strncmp(argv[i],"x",1)==0){
      char *p,ncol=0;
      p=&argv[i][3];
      for(;*p!=0;p++) if(*p==':') ncol++;
//      for(ncol=0;;p++) {
//	if(*p==0) break;
//	if(*p==':') ++ncol;
//	//	if(*p==':') if(++ncol>=3) break;
//      }
      if(ncol>=5) {
	sscanf(&argv[i][1],"%d",&j);
	sscanf(&argv[i][1],"%d:%lf:%lf:%lf:%lf:%lf:%lf:%lf:%lf",&j,&xmin[j],&xmax[j],&xmin1[j],&xmax1[j],&x0m[j],&x1m[j],&x0M[j],&x1M[j]);
	//	xmin[j]=_xm;xmax[j]=_xM;xmin1[j]=_xm1;xmax1[j]=_xM1;
      }
      else if(ncol>=3) {
	if(argv[i][1]!=':'){
	  sscanf(&argv[i][1],"%d:%lf:%lf:%lf:%lf",&j,&_xm,&_xM,&_xm1,&_xM1);
	  xmin[j]=_xm;xmax[j]=_xM;xmin1[j]=_xm1;xmax1[j]=_xM1;
	}
	else{
	  sscanf(&argv[i][2],"%lf:%lf:%lf:%lf",&_xm,&_xM,&_xm1,&_xM1);
	  for(j=0;j<K;j++){xmin[j]=_xm;xmax[j]=_xM;xmin1[j]=_xm1;xmax1[j]=_xM1;}
	}
      }
      else {
	sscanf(&argv[i][1],"%d:%lf:%lf",&j,&_xm,&_xM);xmin[j]=_xm;xmax[j]=_xM;xmax1[j]=1;xmin1[j]=0;
      }
    }
    else if(strncmp(argv[i],"y:",2)==0){
      char *p,ncol=0;
      p=&argv[i][3];
      for(;*p!=0;p++) if(*p==':') ncol++;
      if(ncol>=4) {
	sscanf(&argv[i][2],"%lf:%lf:%lf:%lf:%lf:%lf:%lf:%lf",&ymin,&ymax,&ymin1,&ymax1,&y0m,&y1m,&y0M,&y1M);
      }
      if(ncol>=2) {sscanf(&argv[i][2],"%lf:%lf:%lf:%lf",&ymin,&ymax,&ymin1,&ymax1);}
      else {sscanf(&argv[i][2],"%lf:%lf",&ymin,&ymax);ymin1=0;ymax1=1;}
    }
    else if(strncmp(argv[i],"Tsk:",4)==0){
      char *p,ncol=0,nm=0;
      p=&argv[i][4];
      for(;*p!=0;p++) if(*p==':') ncol++;
      p=&argv[i][4];
      for(;*p!=0;p++) if(*p=='-') nm++;
      if(nm==2 && ncol==2) {sscanf(&argv[i][4],"%d:%d-%d:%d-%d",&Task,&tr0,&tr1,&tp0,&tp1);}
      else{
	int t0,t1,t2;
	if(ncol>=3) {sscanf(&argv[i][4],"%d:%d:%d:%d",&Task,&t1,&t2,&t0);}
	else if(ncol>=2) {sscanf(&argv[i][4],"%d:%d:%d",&Task,&t1,&t2);t0=0;}
	else {sscanf(&argv[i][2],"%d",&Task);t1=t2=t0=0;}
	tr0=t0;
	tr1=t0+t1;
	tp0=t0+t1;
	tp1=t0+t2;
      }
    }
    //    else if(strncmp(argv[i],"y:",2)==0) sscanf(&argv[i][2],"%lf:%lf:%lf:%lf",&ymin,&ymax,&ymin1,&ymax1);
    else if(strncmp(argv[i],"Lstd:",5)==0) sscanf(&argv[i][5],"%lf:%d",&Lstd,&Lstdm);
    else if(strncmp(argv[i],"lc:",3)==0) lcom=&argv[i][3];
    else if(strncmp(argv[i],"lossall:",8)==0) lossall=&argv[i][8];
    else if(strncmp(argv[i],"r0",2)==0) r0=1;
    else if(strncmp(argv[i],"r:",2)==0) sscanf(&argv[i][2],"%d:%d:%lf",&r1,&r2,&r3);
    else if(strncmp(argv[i],"NC:",3)==0) sscanf(&argv[i][3],"%d",&NC);
    else if(strncmp(argv[i],"vt:",3)==0) sscanf(&argv[i][3],"%lf",&vt);
    else if(strncmp(argv[i],"vm:",3)==0) sscanf(&argv[i][3],"%d",&vm);
    else if(strncmp(argv[i],"vr:",3)==0) sscanf(&argv[i][3],"%lf",&vr);
    else if(strncmp(argv[i],"ib:",3)==0){
      int _ibm;
      sscanf(&argv[i][3],"%d",&_ibm);
      if(_ibm==-1) ibmode=NULL;
      else sprintf(ibmode,"%s",argv[i]);
    }
    else if(strncmp(argv[i],"LDm:",4)==0) sscanf(&argv[i][4],"%d",&LDmode);
    else if(strncmp(argv[i],"sd:",3)==0) sscanf(&argv[i][3],"%u",&seed);
    //    else if(strncmp(argv[i],"bst:",4)==0){sscanf(&argv[i][4],"%d",&t_boost);}
    else if(strncmp(argv[i],"bst:",4)==0){
      int _t,_B;
      int ncols=sscanf(&argv[i][4],"%d:%d",&_t,&_B);
      if(ncols>=1) t_boost=_t;
      if(ncols>=2) Boost=_B;
//      char *p;
//      p=&argv[i][4];
//      for(j=0;;p++) {
//	if(*p==0) break;
//	if(*p==':') if(++j>=1) break;
//      }
//      //      if(j>=1) {sscanf(&argv[i][4],"%d:%lf",&t_boost,&rethresh_boost);}
//      if(j>=1) {sscanf(&argv[i][4],"%d:%d",&t_boost,&Boost);}
//      else {sscanf(&argv[i][4],"%d",&t_boost);}
//      //      if(t_boost<0) {Boost=GbBoost;t_boost*=-1;} else Boost=EmBoost;
//      //      sscanf(&argv[i][4],"%d",&t_boost);
    }
    else if(strncmp(argv[i],"t:",2)==0) sprintf(fn_target,"%s",&argv[i][2]);
    else if(strncmp(argv[i],"g:",2)==0) sscanf(&argv[i][2],"%lf",&gamma0);
    else if(strncmp(argv[i],"e:",2)==0) sscanf(&argv[i][2],"%lf",&entropy_thresh);
    else if(strncmp(argv[i],"T:",2)==0) sscanf(&argv[i][2],"%d",&T);
    else if(strncmp(argv[i],"w:",2)==0) sscanf(&argv[i][2],"%lf",&w);
    //    else if(strncmp(argv[i],"E:",2)==0) sscanf(&argv[i][2],"%d",&ERR);
    //    else if(strncmp(argv[i],"p0:",3)==0) sscanf(&argv[i][3],"%d",&p0);
    else if(strncmp(argv[i],"i:",2)==0) sscanf(&argv[i][2],"%d",&intmode);
    else if(strncmp(argv[i],"N:",2)==0){
      char *p=&argv[i][2];
      //      sprintf(fnsubmit,"%s+pred%s.dat",fnbody(fntest,pr),argv[i]);
      //      rsmo[i]=0;
      for(nens=0;;){
	sscanf(p,"%d",&NN[nens++]);
	for(;;) if(*p==':'|| *p=='-' || *p==0) break; else p++;
	if(*p==0) break;
	if(*p=='-'){
	  int _n,j;
	  sscanf(++p,"%d",&_n);
	  for(j=NN[nens-1]+1;j<=_n;j++) NN[nens++]=j;
	  for(;;) if(*p==':'|| *p=='-' || *p==0) break; else p++;
	  if(*p==0) break;
	}
	p++;
      }
    }
    else if(strncmp(argv[i],"rdm:",4)==0) sscanf(&argv[i][3],"%d",&rangedatamode);
    else if(strncmp(argv[i],"ssp:",4)==0){sscanf(&argv[i][4],"%d",&ssp);}
    else if(strncmp(argv[i],"tau:",4)==0){
      char *p,ncol=0;
      p=&argv[i][3];
      for(;*p!=0;p++) if(*p==':') ncol++;
//      for(ncol=0;;p++){
//	if(*p==0) break;
//	if(*p==':') ++ncol;
//      }
      if(ncol>=2) sscanf(&argv[i][4],"%lf:%lf:%lf",&tau_c,&tau_h,&eta1);
      else sscanf(&argv[i][4],"%lf:%lf",&tau_c,&tau_h);
    }
    else if(strncmp(argv[i],"bayes:",6)==0){
      iBayes=i;
      sscanf(&argv[i][6],"%d:%lf:%lf:%d:%d",&Bayes,&BayesLambdaL,&BayesLambdaS,&BayesUseAllData,&Bayesseed);
   } 
    else if(strncmp(argv[i],"nobt:",5)==0){sscanf(&argv[i][5],"%d",&nob_thresh);}
    else if(strncmp(argv[i],"BIAS:",5)==0) sscanf(&argv[i][5],"%lf",&BIAS);
    else if(strncmp(argv[i],"seed:",5)==0) sscanf(&argv[i][5],"%u",&seed);
    else if(strncmp(argv[i],"DISP:",5)==0) sscanf(&argv[i][5],"%d",&DISP);
    else if(strncmp(argv[i],"Tpinv:",6)==0) sscanf(&argv[i][6],"%d",&Tpinv);
    else if(strncmp(argv[i],"nop:",4)==0) sscanf(&argv[i][4],"%d",&nop);
    else if(strncmp(argv[i],"fupdate:",8)==0) fupdate=&argv[i][8];
    else if(strncmp(argv[i],"pupdate:",8)==0) pupdate=&argv[i][8];
    else if(strncmp(argv[i],"e4t:",4)==0) sscanf(&argv[i][4],"%lf",&err4terminate);
    else if(strncmp(argv[i],"e4p:",4)==0) sscanf(&argv[i][4],"%lf",&err4propagate);
  }//end of  for(i=4;i<argc;i++){
#ifdef CHKTIME
  fprintf(stderr,"laptime: reading argv:%f[sec]\n",mytimer_lap());
#endif
  int (*net_printf)(char *);
  if(nop==0) net_printf=printf1;
  else net_printf=noprintf;

  if(Bayes>0 && strncmp(fn_bagging,"/dev/null",9)!=0){
    for(i=0;i<argc;i++){
      if(i==ibg) sprintf(&cmd[strlen(cmd)],"bg:/dev/null");
      else {
	int lc=strlen(cmd),la=strlen(argv[i]);//char *oldcmd=cmd;
	if(cmdsize<lc+la){
	  cmdsize=lc+10*la;
	  cmd=(char*)realloc(cmd,cmdsize);
	  //	  if(cmd!=oldcmd) {int j; for(j=0;j<lc;j++) cmd[j]=oldcmd[j]; }
	}
	sprintf(&cmd[lc],"%s ",argv[i]);
      }
    }
    fprintf(stderr,"Recursive ensrs '%s'.\n",cmd); 
//    if(strlen(cmd)>=cmdsize){
//      fprintf(stderr,"enlarge cmdsize(%d) larger than %d.\n",cmdsize,strlen(cmd));
//      return(-1);
//    }
    system(cmd);//recursive call for ensrs ?
  }
  //  for(i=0;i<K;i++) fprintf(stderr,"%e %e %e %e   #%d ???2入力の正規化 x0min0 x0max0 →x0min x0max \n",xmin[i],xmax[i],xmin1[i],xmax1[i],i);
#ifdef CHKTIME
  fprintf(stderr,"laptime: recursive Bayes:%f[sec]\n",mytimer_lap());
#endif
  {///get num of given data
    if((fp1=fopen(fn_givendata,"r"))==NULL){
      printf("file(%s) open error\n",fn_givendata);
      return(-1);
    }
    nGivenData=0;
    for(;;){
      line1[0]=line2[0]=0;
      fgets(line1,LINESIZE,fp1);
      if(feof(fp1)) break;
      if(line1[0]==0) break;
      nGivenData++;
    }
    fclose(fp1);
  }
  if(nFolds>=nFoldsmax){            //nFoldsのエラー処理
    if(method==MULTIFOLD && nFolds>nGivenData ){            //nFoldsのエラー処理
      net_printf(strx2("nFolds must be 0<nFolds<%d<nGivenData=%d. But no ristriction",nFoldsmax,nGivenData));
      nFolds=nGivenData-1;//050819??
    }
    else nFolds=nFoldsmax;
  }
    
  net_printf(strx2("number of given data =%d, dimension of the input vector=%d\n",nGivenData,K));

  sse=(double*) malloc(sizeof(double)*(nGivenData));

  is=(int*) malloc(sizeof(int)*(nFolds));
  mse=(double*) malloc(sizeof(double)*(nFolds));
  Lhat=(double*) malloc(sizeof(double)*(nFolds));
  n_mse=(int*) malloc(sizeof(int)*(nFolds));
  msetrain=(double*) malloc(sizeof(double)*(nFolds));
  n_msetrain=(int*) malloc(sizeof(int)*(nFolds));
  n_cells2=(int*) malloc(sizeof(int)*(nFolds));
  //  int *nc; double *yp0,*yt,*yr;//for BAGGING==2
  //  int *nc; double *yp0,*yt;//for BAGGING==2
  double *pbst,*pbst0,*pbst1,*pbst2,*pbst_1;
  double *ebst,*ebst0,*ebst1,*ebst2;
  double *dEdp;
  //  double *cbst,*ypt,*wbst,*wbst1,tau_c=8.,tau_h=4;
  double *bias2,*var,*errt;
  double pthresh=0.01/nGivenData;//good
  //      double pthresh=0.1/nGivenData;//good
  //      double pthresh=0.0;
  //  if(t_boost==-1);
  //  //  else if(t_boost==-2){//Gradient-based boosting
  //  else if(t_boost==GbBoost){//Gradient-based boosting
  if(Boost==GbBoost){//Gradient-based boosting
    cbst=(double*)malloc(sizeof(double)*(t_boost+1));
    if(BAGGING==BAGwithoutVal){
      ypt=(double*)malloc(sizeof(double)*nGivenData*(t_boost+2));
      bias2=(double*)malloc(sizeof(double)*nGivenData*(t_boost+2));
      errt=(double*)malloc(sizeof(double)*nGivenData);
      var=(double*)malloc(sizeof(double)*nGivenData*(t_boost+2));
    }
    else ypt=(double*)malloc(sizeof(double)*nValData*(t_boost+2));
    wbst=(double*)malloc(sizeof(double)*(nGivenData)*2);
    wbst1=&wbst[nGivenData];
    pbst_1=(double*)malloc(sizeof(double)*(nGivenData)*3);
    pbst=&pbst_1[nGivenData];
    pbst1=&pbst_1[nGivenData*2];
    int t;
    //    int sizeofline1=sizeof(line1);

    if(BAGGING==BAGwithoutVal){//
      double _dum;
      int t_boost1=t_boost-1,t_boost2=t_boost-2;
      for(t=0;t<t_boost;t++){//t_boost=1,2,3,...
	sprintf(fn_boost,"%s/tmp/BstGb:%s+m%da%gs%dN:%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t);//
	//fn_boost involves 
	//cbst[t]
	//pbst[(t)*nGivenData+i],wbst[(t)*nGivenData+i],yp[t*nGivenData+i],pbst[(t+1)*nGivenData+i],wbst[(t+1)*nGivenData+i]
	if((fp=fopen(fn_boost,"r"))==NULL){
	  double pbst00=1./nGivenData;
	  for(i=0;i<nGivenData;i++) pbst[i]=wbst[i]=pbst00;
	  cbst[t]=0.0;//this will be overwritten.
	} 
	else {
	  fgets(line1,LINESIZE,fp);sscanf(line1,"%lf",&cbst[t]);
	  if(t==t_boost1){
	    for(i=0;i<nGivenData;i++){
	      fgets(line1,LINESIZE,fp);//	      fgets(line1,sizeofline1,fp);
	      sscanf(line1,"%lf%lf%lf%lf%lf%lf",&_dum,&_dum,&ypt[t*nGivenData+i],&pbst[i],&wbst[i],&bias2[t*nGivenData+i]);//needs pbst[t_boost] and wbst[t_boost]
	    }
	  }
	  else if(t==t_boost2){
	    for(i=0;i<nGivenData;i++){
	      fgets(line1,LINESIZE,fp);//	      fgets(line1,sizeofline1,fp);
	      sscanf(line1,"%lf%lf%lf%lf%lf%lf",&_dum,&_dum,&ypt[t*nGivenData+i],&pbst_1[i],&_dum,&bias2[t*nGivenData+i]);//needs pbst[t_boost] and wbst[t_boost]
	    }
	  }
	  else {
	    for(i=0;i<nGivenData;i++){
	      fgets(line1,LINESIZE,fp);
	      //	      sscanf(line1,"%lf",&ypt[t*nGivenData+i]);
	      sscanf(line1,"%lf%lf%lf%lf%lf%lf",&_dum,&_dum,&ypt[t*nGivenData+i],&_dum,&_dum,&bias2[t*nGivenData+i]);
	    }
	    fclose(fp);
	  }
	}
      }
    }
    else if(BAGGING==BAGwithVal){//after BAGwithoutVal for(t=1;t<=t_boost;t++)
      for(t=1;t<=t_boost;t++){//t_boost=1,2,3,...
	sprintf(fn_boost,"%s/tmp/BstGb:%s+m%da%gs%dN:%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t);//
	//fn_boost for t involves 
	//cbst[t]
	//pbst[(t)*nGivenData+i],wbst[(t)*nGivenData+i],yp[t*nGivenData+i],pbst[(t+1)*nGivenData+i],wbst[(t+1)*nGivenData+i]
	if((fp=fopen(fn_boost,"r"))==NULL){
	  double pbst00=1./nGivenData;
	  for(i=0;i<nGivenData;i++) pbst[i]=wbst[i]=pbst00;
	  cbst[t]=0.0;//will not encounter
	}
	else {
	  fgets(line1,LINESIZE,fp);sscanf(line1,"%lf",&cbst[t]);
	  if(t==t_boost){//needs only pbst[t_boost] for BAGwithVal, overwritten for overwritten for BAGwithoutVal
	    for(i=0;i<nGivenData;i++){
	      fgets(line1,LINESIZE,fp);
	      sscanf(line1,"%lf",&pbst[i]);//needs pbst[t_boost]
	    }
	  }
	}
      }
      for(t=1;t<t_boost;t++){//t_boost=1,2,3,...
	sprintf(fn_boost,"%s/tmp/BstGb:%s+%s+m%da%gs%dN:%d:%db%dpred.dat",rdir0,fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,seed,NN[0],NStep,t);//
	if((fp=fopen(fn_boost,"r"))==NULL);
	else{
	  for(i=0;i<nValData;i++) {
	    fgets(line1,sizeof(line1),fp); 
	    sscanf(line1,"%lf",&ypt[t*nValData+i]);
	  }
	  fclose(fp);
	}
      }
    }
  }
  else if(Boost==EmBoost){// if(t_boost>=0){//Entropy minimization?
    pbst=(double*) malloc(sizeof(double)*(nGivenData)*5);
    ebst=(double*) malloc(sizeof(double)*(nGivenData)*5);
    pbst0=&pbst[nGivenData];
    ebst0=&ebst[nGivenData];
    pbst1=&pbst[nGivenData*2];
    ebst1=&ebst[nGivenData*2];
    pbst2=&pbst[nGivenData*3];
    ebst2=&ebst[nGivenData*3];
    dEdp =&pbst[nGivenData*4];
    sprintf(fn_boost,"%s/tmp/BstEM:%s+m%da%gs%dN:%d-%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NN[nens-1],NStep,t_boost-1);//
    ///
    if((fp=fopen(fn_boost,"r"))==NULL){
      double pbst00=1./nGivenData;
      for(n=0;n<nGivenData;n++){
	pbst[n]=pbst0[n]=pbst1[n]=pbst2[n]=pbst00;
	ebst0[n]=ebst1[n]=ebst2[n]=0.0;
      }
    }
    //    if(fp=fopen(fn_boost,"r")!=NULL){//
    else{//usually t_boost>=2, or t_boost==1 for outlier
      //      sprintf(fn_boost,"%s/tmp/%s+m%da%.2fs%dF%dN:%d-%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,nFolds,NN[0],NN[nens-1],NStep,t_boost-1);//
      double emean0,emean1,emax=-1e20,emin=1e20;
      double esum0=0,esum1=0;
#define ENTROPY4BOOST
#ifdef ENTROPY4BOOST
      double H0,H1;//entropy
#endif
      //      double remax=-1e20;
      int nGivenData1=0;
      //      int sizeofline1=sizeof(line1);
      for(n=0;n<nGivenData;n++){
	fgets(line1,LINESIZE,fp);
	//	sscanf(line1,"%lf%lf%lf%lf%lf%lf",&pbst0[n],&ebst0[n],&pbst1[n],&ebst1[n],&pbst2[n],&ebst2[n]);
	sscanf(line1,"%lf%lf%lf%lf%lf%lf",&pbst0[n],&pbst1[n],&pbst2[n],&ebst0[n],&ebst1[n],&ebst2[n]);
	//	if(n>=200 && t_boost>0) pbst0[n]=0;//????

	if(pbst0[n]>pthresh){//??
	  esum0 += fabs(ebst0[n]);
	  esum1 += fabs(ebst1[n]);
	  if(emin>ebst0[n]) emin=ebst0[n];
	  if(emax<ebst0[n]) emax=ebst0[n];
	  nGivenData1++;
	}
      }
      fclose(fp);
      emean0=esum0/nGivenData1;
      emean1=esum1/nGivenData1;
      if(t_boost==1){//After t_boost==0 for outlier processing
	double pbst00=1./nGivenData;
	for(n=0;n<nGivenData;n++){
	  if(pbst0[n]>pthresh && ebst0[n]/emean0<rethresh_boost){//??
	    pbst[n]=pbst0[n]=pbst1[n]=pbst2[n]=pbst00;
	    ebst0[n]=ebst1[n]=ebst2[n]=0.0;
	  }
	  else{
	    fprintf(stderr,"n=%d is supposed to be outlier.\n",n);
	    pbst[n]=pbst0[n]=0.0;
	  }
	}      
      }

#ifdef ENTROPY4BOOST
      double pn;
      H0=H1=0;
      for(n=0;n<nGivenData;n++){
	if(pbst0[n]>pthresh){//??
	  pn=ebst0[n]/esum0+1e-20;H0 -= (pn)*log(pn);
	  pn=ebst1[n]/esum1+1e-20;H1 -= (pn)*log(pn);
	}
      }
      H0/=log(nGivenData1);
      H1/=log(nGivenData1);
#endif
      //      if(t_boost==2){
      //	double emean00=0;
      //	nGivenData1=0;
      //	for(n=0;n<nGivenData;n++){
      //	  if(ebst0[n]/emean0<rethresh_boost){
      //	    emean00+=ebst0[n];
      //	    nGivenData1++;
      //	  }
      //	}
      //	emean0=emean00/nGivenData1;
      //      }
#define NewdEdp
#ifdef NewdEdp
      double dEdpmax=-1e30,dEdpmin=1e30;
      double pmax=-1e30,pmin=1e30;
      //#ifdef NewdEdpOld
      {
	double dpthresh=0.000001/nGivenData;
	//	double dpthresh=0.00001/nGivenData;
	//	double dpthresh=0.0001/nGivenData;
	//	double dpthresh=0.01/nGivenData;//Cf. pthresh=0.01/nGivenData;//good
	for(n=0;n<nGivenData;n++){
	  if(pbst0[n]<=pthresh) dEdp[n]=0;
	  else if(ebst0[n]/emean0>rethresh_boost) dEdp[n]=0;
	  else{//if(pbst0[0]>pthresh) dEdp[n]=0;
	    double dp0=pbst0[n]-pbst1[n];
	    if(t_boost<=2){
	      //if(t_boost<=3){//for outlier ??
#ifdef ENTROPY4BOOST
	      dEdp[n]=-(ebst0[n]-emean0);//for insufficient learning
#else	      
	      //dEdp[n]=emean0-ebst0[n]/(1.+exp(ebst0[n]-emean0*9.));//for ???
	      dEdp[n]=-(ebst0[n]-emean0);//for insufficient learning
	      //	      dEdp[n]=(ebst0[n]-emean0);//for outliers??
	      //	      dEdp[n]/=(1.+exp(ebst0[n]-emean0*16.));//for ???
	      //	      dEdp[n]/=(1.+exp(ebst0[n]-emean0*9.));//for ???
	      //		dEdp[n]=(ebst0[n]-emean0)*log(1.+ebst0[n]/emean0/4.);//for ???
	      //	    dEdp[n]=(ebst0[n]-emean0);//for outliers 
	      //	    if(ebst0[n]<0) fprintf(stderr,"strange ebst0[%d]=%e\n",n,ebst0[n]);//because using ebst[n]-=uc[n]
#endif //ENTROPY4BOOST
	    }
	    //	  if(t_boost<=2) dEdp[n]=-(fabs(ebst0[n])-emean0);
	    //if(t_boost<=2) dEdp[n]=1./ebst0[n];
	    else{
	      if(fabs(dp0)<dpthresh){
		fprintf(stderr,"!!!Optimal? n%d dp0=%e=%e-%e<%e, e0=%e e1=%e.!!!\n",n,dp0,pbst0[n],pbst1[n],dpthresh,ebst0[n],ebst1[n]);
		dEdp[n]=0;
		//		if(dp0>0) dEdp[n]=dpthresh;else dEdp[n]=-dpthresh;
		//	      dEdp[n]=-0.001*(ebst0[n]-emean0)/(emax-emin);//for insufficient learning
		//	      dEdp[n]=0.001*(1./nGivenData-pbst0[n]);
		//	      if(ebst0[n]>ebst1[n]) pbst0[n]=pbst1[n];
	      }
	      //	  else{
	      //	    if(fabs(dp0)<1e-30){
	      //	      fprintf(stderr,"!!!!!!!!!!!!!!!!!Optimal? n%d p%e e0=%e e1=%e.!!!!!!!!\n",n,pbst0[0],ebst0[n],ebst1[n]);
	      //	      //	      dEdp[n]=0;
	      //	      dEdp[n]=0.01*(1./nGivenData-pbst0[n]);
	      //	      //	      if(ebst0[n]>ebst1[n]) pbst0[n]=pbst1[n];
	      //	    }
#define ERR2BOOST
#ifdef ERR2BOOST
	      else{
#ifdef ENTROPY4BOOST
		dEdp[n]=-(H0-H1)/dp0;//for minimizing -H0 == maximizing H0
		//dEdp[n]=(H0-H1)/dp0;//for minimizing H0
		//dEdp[n]=-(ebst0[n]-emean0);//for insufficient learning
#else //not ENTROPY4BOOST
		dEdp[n]=(emean0-emean1)/dp0;//good 060731
		//		dEdp[n]=-(ebst0[n]-emean0);//SoSo 060810
		//		dEdp[n]=(emean0-emean1-(ebst0[n]-ebst1[n])/nGivenData1)/dp0;//good for underlearning
		//		dEdp[n]=-(ebst0[n]-emean0);//for insufficient learning
		//		dEdp[n]/=(1.+exp(ebst0[n]-emean0*16.));//for ???
		//		dEdp[n]/=(1.+exp(ebst0[n]-emean0*9.));//for ???
		//dEdp[n]=(emean0-emean1-(ebst0[n]-ebst1[n])/nGivenData1)/dp0;//good for underlearning
		//		dEdp[n]=(emean0-emean1)/dp0;//better for outlier and furtherlearning?
		//	      if(dp0>0) dEdp[n]=((emean0-emean1)-(ebst0[n]-ebst1[n])/dp0/nGivenData);
		//	      else dEdp[n]=-((emean0-emean1)-(ebst0[n]-ebst1[n])/dp0/nGivenData);
#endif  //ENTROPY4BOOST
	      }
	      //	    else dEdp[n]=(emean0-emean1)/dp0;//??for insufficient learning also, this works well for neglecting outliers
	      //	    else dEdp[n]=-pbst0[n]*(ebst0[n]-emean0);//
	      //	    else dEdp[n]=-(ebst0[n]-emean0);//good
	      //	    else dEdp[n]=1./ebst0[n];//
	      //	    else dEdp[n]=-(ebst0[n]-emean0)*dp0;//
	      //	    else dEdp[n]=-(emean0-emean1)/dp0;//this works well for neglecting outliers?
	      //	    else dEdp[n]=-(ebst0[n]-ebst1[n])/dp0;
#else //not ERR2BOOST
	      else dEdp[n]=-ebst0[n]*(ebst0[n]-ebst1[n])/dp0;
	      //	    else{
	      //	      if(fabs(pbst1[n]-pbst2[n])<1e-30) {
	      //		dEdp[n]=0;
	      //		fprintf(stderr,"???????????? Optimal? n%d p%e e0=%e e1=%e.!!!!!!!!\n",n,pbst0[0],ebst0[n],ebst1[n]);
	      //		if(ebst0[n]>ebst1[n]) pbst0[n]=pbst1[n];
	      //	      }
	      //	      else dEdp[n]= square(ebst0[n]-ebst1[n])/dp0+ebst0[n]*((ebst0[n]-ebst1[n])/dp0-(ebst1[n]-ebst2[n])/(pbst1[n]-pbst2[n]));
	      //	    }
#endif//ifdefERR2BOOST
	      //	    dEdp[n]-=0.001/nGivenData;//for outlier? decrease if 
	    }
	    //	  if(dEdpmax<dEdp[n]) dEdpmax=dEdp[n];
	    //	  if(dEdpmax<fabs(dEdp[n])) dEdpmax=fabs(dEdp[n]);
	    //	  }
	    if(dEdpmax<dEdp[n]) dEdpmax=dEdp[n];
	    if(dEdpmin>dEdp[n]) dEdpmin=dEdp[n];
	    if(pmax<pbst0[n]) pmax=pbst0[n];
	    if(pmin>pbst0[n]) pmin=pbst0[n];
	    //	  if(n>=199){
	    //	    fprintf(stderr,">>>>de(%d)=%e=(%e)-(%e).dEdp=%e=%e-%e,emean0=%e,emean1=%e\n",n,ebst0[n]-ebst1[n],ebst0[n],ebst1[n],pbst0[n]-pbst1[n],pbst0[n],pbst1[n],emean0,emean1);
	  }
	}//closing of for(n=0;n<nGivenData;n++){
      }
      //      if(fabs(dEdpmax-dEdpmin)<1.0e-20) dEdpmax=dEdpmin+1;
      //      double lambda=(1./nGivenData)*0.8/(dEdpmax-dEdpmin)/(t_boost-1.);
      //      double lambda=(1./nGivenData)*1.0/(-dEdpmin)/(t_boost-1.);
      //      double lambda=pmax*1.0/t_boost/dEdpmax;
      //      double lambda=pmax*1.5/t_boost/dEdpmax;
      //      double lambda=pmax*0.25/dEdpmax;
      //      double lambda=pmax*0.9/t_boost/dEdpmax;
      //      double lambda=pmax*1.5/(1.+t_boost/1.)/dEdpmax;
      //      double lambda=pmax*0.4*(1./t_boost)/dEdpmax;
      //      double lambda=pmax/(1.+t_boost/1.)/dEdpmax;
      //      double lambda=pmax/(1.+t_boost/10.)/dEdpmax;
      //      double lambda=pmax/(2.+t_boost/2.)/dEdpmax;
      //      double lambda=pmax*0.1/dEdpmax;
      //      double ethresh=emean0*0.0001;
      //      double lambda=(0.9/nGivenData)/dEdpwidth/(t_boost-1.);
      //      double lambda=(1.0/nGivenData)/dEdpwidth/(t_boost-1.);
      //      double lambda=1.0/dEdpwidth/(t_boost-1.);
      //      double dEdpwidth=dEdpmax-dEdpmin;double lambda=1.0/dEdpwidth/(1.+(t_boost-2.)/0.2);
      //      double lambda=1.0/(emax-emin)/(1.+(t_boost-2.)/0.2);
      //      double lambdae=1.0/(emax-emin)/(1.+(t_boost-2.)/0.2);//good1
      //      double lambdad=1.0/(dEdpmax-dEdpmin)/(1.+(t_boost-2.)/0.2);//good0 for no-outlier
      double lambdae,lambdad;
      //      if(t_boost<=2) lambda=1.0;
      //#endif //ifdef NewdEdpOLD
      //      double lambda=1.0/(emax-emin)/(1.+(t_boost-2.)/0.2);
      double psum=0;
      double dEdpw=dEdpmax-dEdpmin;
      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.2))/dEdpw;//OK for all
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.08))/dEdpw;//OK for all
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/dEdpw;//OK underlearning but slow NG for rsa:2:0.7:4321:100bst
      //if(dEdpw<1.e-20) lambdad=0;else if(t_boost<=2) lambdad=1./dEdpw; else lambdad=(1./(1.+(t_boost-3.)/0.1))/dEdpw;
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else if(t_boost<=2) lambdad=1./(dEdpmax-dEdpmin); else lambdad=(1./(1.+(t_boost-3.)/0.1))/(dEdpmax-dEdpmin);//
      //if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/(dEdpmax-dEdpmin);//OK underlearning but slow
      //if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.1)/(dEdpmax-dEdpmin);//SoSo Ltst0=9.190e-02 for t:sin2e2noise0.5 v:sin1e3noise0.5test 
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.01)/(dEdpmax-dEdpmin);//SoSo Ltst0=9.152e-02 for t:sin2e2noise0.5 v:sin1e3noise0.5test
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.001)/(dEdpmax-dEdpmin);//NG Ltst0=9.362e-02 for t:sin2e2noise0.5 v:sin1e3noise0.5test 
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/(dEdpmax-dEdpmin);//OK underlearning but slow
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.2))/(dEdpmax-dEdpmin);//SoSo underlearning but slow
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.05))/(dEdpmax-dEdpmin);//Not So Good underlearning
      //if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.98/(1.+(t_boost-2.)/0.01)+0.02)/(dEdpmax-dEdpmin);//NG underlearning
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.999/(1.+(t_boost-2.)/0.01)+0.001)/(dEdpmax-dEdpmin)/nGivenData1;//??underlearning
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(0.97/(1.+(t_boost-2.)/0.01)+0.03)/(dEdpmax-dEdpmin)/nGivenData1;//??NGunderlearning
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.05))/(dEdpmax-dEdpmin);//underlearning
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0; else lambdad=(0.999/(1.+(t_boost-2.)/0.1)+0.001)/(dEdpmax-dEdpmin);//underlearning
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0; else lambdad=(0.99/(1.+(t_boost-2.)/0.1)+0.01)/(dEdpmax-dEdpmin);//underlearningNG
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0; else lambdad=(0.95/(1.+(t_boost-2.)/0.05)+0.05)/(dEdpmax-dEdpmin);//good0 for no-outlier
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0; else lambdad=(0.98/(1.+(t_boost-2.)/0.1)+0.02)/(dEdpmax-dEdpmin);//underlearning NG
      if(emax-emin<1.e-20) lambdae=0;else lambdae=(1./(1.+(t_boost-2.)/0.1))/(emax-emin);//for outlier
      //      if(emax-emin<1.e-20) lambdae=0.; else lambdae=(0.999/(1.+(t_boost-2.)/0.1)+0.001)/(emax-emin);//
      //      if(emax-emin<1.e-20) lambdae=0;else lambdae=(0.99/(1.+(t_boost-2.)/0.1)+0.01)/(emax-emin);//for outlierNG
      //      if(emax-emin<1.e-20) lambdae=0;else lambdae=(0.95/(1.+(t_boost-2.)/0.1)+0.05)/(emax-emin);//for outlierNG
      //      if(emax-emin<1.e-20) lambdae=0;else lambdae=(0.95/(1.+(t_boost-2.)/0.1)+0.05)/(emax-emin);//for outlier
      //      if(emax-emin<1.e-20) lambdae=0; else lambdae=(0.999/(1.+(t_boost-2.)/0.1)+0.001)/(emax-emin);//for outlier
      //      if(dEdpmax-dEdpmin<1.e-20) lambdad=0; else lambdad=(0.98/(1.+(t_boost-2.)/0.1)+0.02)/(dEdpmax-dEdpmin);//good0 for no-outlier
      //      if(emax-emin<1.e-20) lambdae=0; else lambdae=(0.98/(1.+(t_boost-2.)/0.1)+0.02)/(emax-emin);
      //      fprintf(stderr,"lambda=%e,%e\n",lambdad,lambdae);
      //      if(lambdad<0.01) lambdad=0.01;if(lambdae<0.01) lambdae=0.01;
      //      {//for check
      //	double lambdadEdp=1.0/(dEdpmax-dEdpmin)/(1.+(t_boost-2.)/0.2);
      //	double lambdae =1.0/(emax-emin)/(1.+(t_boost-2.)/0.2);
      //	fprintf(stderr,"****width=%e=(%e)-(%e) ewidth=%e=(%e)-(%e)",dEdpmax-dEdpmin,dEdpmax,dEdpmin,emax-emin,emax,emin);
      //	fprintf(stderr,"****lambda=%e,%e\n",lambdadEdp,lambdae);
      //      }
      //if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-3.)/0.1))/dEdpw;//SoSo-1
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-1.)/0.1))/dEdpw;//SoSo
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/dEdpw/nGivenData;//??
      //if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(2.+(0)/0.1))/dEdpw/nGivenData;//??SoSo?
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/dEdpw/nGivenData;//??
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-1.)/0.1))/dEdpw/nGivenData;//??
      //if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+t_boost/0.1))/dEdpw/nGivenData;//??
      //      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/dEdpw/nGivenData;//??
      if(dEdpw<1.e-20) lambdad=0;else lambdad=(1./(1.+(t_boost-2.)/0.1))/dEdpw;//SoSo lambdad=0 if t_boost==1
      for(n=0;n<nGivenData;n++){
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.05*ebst0[n]/emean0)/(10.+exp(ebst0[n]/emean0-9.));//
	pbst[n]=pbst0[n]-lambdad*dEdp[n]/nGivenData;
	//	pbst[n]=(pbst0[n]-lambdad*dEdp[n]);
	//	pbst[n]=(pbst0[n]-lambdad*dEdp[n])*(1.+0.05*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.05*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//good for underlearng & outlier
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]);//Not So Good
	//	pbst[n]=(1.-lambdad*dEdp[n])*(1.+0.05*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//not so good for underlearning
	//pbst[n]=(1.-lambdad*dEdp[n])*(1.+0.05*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//NG 
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]);
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.1*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.04*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.03*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])*(1.+0.02*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n])/(1.+exp(ebst0[n]-emean0*9.));//for ???
//	if(ebst0[n]/emean0<rethresh_boost) pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]);//for underlearning
//	else{
//	  //	  pbst[n]=pbst0[n]*(1.-lambdae*(ebst0[n]-emin));//sosogood for outlier
//	  fprintf(stderr,"***boost lambdae n%d,re%6.2f=%.3e/%.3e>%.3fe, p%.3edp%+.3e\n",n,ebst0[n]/emean0,ebst0[n],emean0,rethresh_boost,pbst0[n],pbst0[n]-pbst1[n]);
//	  if(t_boost>=2) pbst[n]=0.0;
//	  //if(t_boost>=3) pbst[n]=0.0;
//	}
//	pbst[n]/=(1.+exp(ebst0[n]-emean0*12.));//for ???
//	pbst[n]/=(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	if(t_boost<=2) pbst[n]=pbst0[n]*(1.+0.02*ebst0[n]/emean0)/(1.+exp(ebst0[n]-emean0*9.));//for ???
	//	double wp=1./(1.+exp(ebst0[n]-emean0*16.0));pbst[n]=(1.-wp)*wp+wp*pbst0[n]*(1.-lambdad*dEdp[n]);//(#1-1+#2) aggregated
	//	double wp=1./(1.+exp(ebst0[n]-emean0*9.0));pbst[n]=(1.-wp)*(pbst0[n]*(1.-lambdae*(ebst0[n]-emin)))+wp*pbst0[n]*(1.-lambdad*dEdp[n]);//(#1-1+#2) aggregated
	//	pbst[n]/=(1.+exp(ebst0[n]-emean0*16.0));//for ???

	//	pbst[n]=(pbst0[n]*(1.-lambdad*dEdp[n]));//#2 for undertrainig
	//	pbst[n]=(pbst0[n]*(1.-lambdae*(ebst0[n]-emin)));//#1-1 better for outliers?
	//	if(ebst0[n]/emean0<12.) pbst[n]=pbst0[n]-lambdad*dEdp[n];//for underlearning
	//	pbst[n]=pbst0[n]-lambdad*dEdp[n];//#2 for undertrainigNG
	//	if(ebst0[n]<12.*emean0) pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]); else pbst[n]=pbst0[n]*(1.-lambdae*(ebst0[n]));//NG
	//	if(ebst0[n]<12.*emean0) pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]); else pbst[n]=pbst0[n]*(1.-lambdae*(ebst0[n]-12.*emean0));//NG
//	if(n>198){
//	  if(ebst0[n]<12.*emean0) fprintf(stderr,"%d)e%.3e,p%.3e,dp%+.3e,*ldEdp=%+.3e=%.3e*%+.3e, le=%+.3e=%.3e*%+.3e\n",n,ebst0[n]/emean0,pbst0[n],pbst0[n]-pbst1[n],lambdad*dEdp[n],lambdad,dEdp[n],lambdae*(ebst0[n]-emin),lambdae,ebst0[n]-emin);
//	  else fprintf(stderr,"%d)e%.3e,p%.3edp%+.3e ldEdp=%+.3e=%.3e*%+.3e, *le=%+.3e=%.3e*%+.3e\n",n,ebst0[n]/emean0,pbst0[n],pbst0[n]-pbst1[n],lambdad*dEdp[n],lambdad,dEdp[n],lambdae*(ebst0[n]-emin),lambdae,ebst0[n]-emin);
//	}
	//	pbst[n]=0.999*pbst0[n] - lambdad*dEdp[n]*pbst0[n];//#3 for undertrainig and outliers?
	//	pbst[n]=0.9999*pbst0[n] - lambdad*dEdp[n];//#3 for undertrainig and outliers?
	//	pbst[n]=0.9999*pbst0[n] - lambdad*dEdp[n]*pbst0[n];//#3 for undertrainig and outliers? 0.9999 is same as 1.0 untilt=4
	//	pbst[n]=0.99*pbst0[n] - lambdad*dEdp[n]*pbst0[n];//#3 for undertrainig and outliers?
	//pbst[n]=(pbst0[n]*(1.-lambdae*(ebst0[n]-emean0)));//#1-2 good for outliers?
	//	if(ebst0[n]<16.*emean0) pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]); else pbst[n]=pbst0[n]*(1.-lambdae*(ebst0[n]-emin));
	//	if(ebst0[n]<16.*emean0) pbst[n]=0; else pbst[n]=pbst0[n]*(1.-lambdae*(ebst0[n]-emin));
	//	double wp=exp(-ebst0[n]/(emean0*16.));pbst[n]=(1.-wp)*(pbst0[n]*(1.-lambdae*(ebst0[n]-emin)))+wp*pbst0[n]*(1.-lambdad*dEdp[n]);//(#1-1+#2) aggregated
	//	pbst[n]=(pbst0[n]*(1.-lambdad*dEdp[n]))-0.1*pbst0[n]*lambdae*(ebst0[n]-emin);//#2-#3 for undertraining and outlier?NG
	//	double wp=exp(-ebst0[n]/(emean0*9));pbst[n]=(1.-wp)*(pbst0[n]*(1.-lambdae*(ebst0[n]-emin)))+wp*pbst0[n]*(1.-lambdad*dEdp[n]);//(#1-1+#2) aggregated
	//	if(ebst0[n]<16*emean0) pbst[n]=(pbst0[n]*(1.-lambdad*dEdp[n]));else{pbst[n]=(pbst0[n]*(1.-lambdae*(ebst0[n]-emin)));fprintf(stderr,"?? outlier? e[%d]/emean=%e\n",n,ebst0[n]/emean0);}
	//pbst[n]=(pbst0[n]*(1.+lambdae*(ebst0[n]-emin)));//??
	//	pbst[n]=(pbst0[n]*(1.-lambdae*(ebst0[n]-emin))+pbst0[n]*(1.-lambdad*dEdp[n]));//(#1-1+#2) aggregated
	//	pbst[n]=(pbst0[n]*(1.-lambdae*(ebst0[n]-emean0))+pbst0[n]*(1.-lambdad*dEdp[n]));//(#1-2+#2) aggregated
	//	if(fabs(dEdp[n])>1e-20) pbst[n]=pbst0[n]*(1.-lambdad*dEdp[n]);else{pbst[n]=pbst0[n]+0.01*(1./nGivenData-pbst0[n]);}//good0 for no-outlier
	//pbst[n]=pbst0[n]*(1.-lambda*(ebst0[n]-emean0));//good1
	//	pbst[n]=pbst0[n]*(1.+lambda*(ebst0[n]-emean0));//-good1
	//pbst[n]=pbst0[n]*(1.+lambda*(ebst0[n]-emin));//for no-outliers?
	//	pbst[n]=pbst0[n]*(1.-lambda*(ebst0[n]-emin));//good for outliers
	//	  pbst[n]=pbst0[n]-lambda*pbst0[n]*(ebst0[n]-emean0);
	//	  pbst[n]=pbst0[n]+lambda*dEdp[n];
	//	  pbst[n]=pbst0[n]*(1.+lambda*dEdp[n]);
	//	if(fabs(dEdp[n])>1e-20) pbst[n]=pbst0[n]+lambda*dEdp[n];else pbst[n]=pbst0[n]+0.01*(1./nGivenData-pbst0[n]);//soso
	//	if(t_boost==2) pbst[n]=1./ebst0[n];else pbst[n]=pbst0[n]+lambda*dEdp[n];
	//	if(ebst0[n]>ethresh) pbst[n]=1./(ebst0[n]);else pbst[n]=1./ethresh;
	//	else pbst[n]=pbst0[n]+0.1*(1./nGivendata-pbst0[n]);
	///////
	//	if(pbst[n]<=pthresh) pbst[n]=pthresh;
	if(pbst[n]<=pthresh) pbst[n]=0.0; 
	//	if(pbst[n]<=pthresh) pbst[n]=0.0; else pbst[n]=1;//Ltst0=1.330e-02 for 
	/////
	//	fprintf(stderr,"%d %e #n pbst\n",n,pbst[n]);
	psum+=pbst[n];
      }
      //      fprintf(stderr,"psum%flambda%.3e,pmax%.3epmin%.3edEdpmax%.3e.\n",psum,lambda,pmax,pmin,dEdpmax);
      fprintf(stderr,"psum%fp[%.3e,%.3e]dEdp[%.3e,%.3e]lambdad%.3e,lambdae%.3e,remax=%.3f=%.3e/%.3e\n",psum,pmin,pmax,dEdpmin,dEdpmax,lambdad,lambdae,emax/emean0,emax,emean0);
      for(n=0;n<nGivenData;n++) pbst[n]/=psum;
#else
#endif//#ifdef OldEdp
    }
    //    init_walkeralias(nGivenData,pbst);//for bootstrap resampling with pbst
  }//closing else if(t_boost>=0)
  if(Boost!=NoBoost) init_walkeralias(nGivenData,pbst);//for bootstrap resampling with pbst
#ifdef CHKTIME
  fprintf(stderr,"laptime: after init_walkeralias:%f[sec]\n",mytimer_lap());
#endif
  if((pdir=opendir("./tmp"))==NULL) system("mkdir ./tmp"); else closedir(pdir);
  sprintf(fnbody_givendata,"%s",fnbody(fn_givendata,pr));
  nTestData=(int*) malloc(sizeof(int)*nFolds);
  ////make fntrain fntest
  if(chkomit==0) setfntraintest(seed,fn_givendata,argv,fupdate);
  //  if(chkomit==0) setfntraintest(fntrain,fntest,seed,fn_givendata,argv);
  if(BAGGING==BAGwithoutVal){
    yp0=(double*)malloc(sizeof(double)*nGivenData*nFolds);//for BAGGING==2 & EmBoost
    yt=(double*)malloc(sizeof(double)*nGivenData);//for BAGGING==2
  }
  else yt=(double*)malloc(sizeof(double)*nValData);//for BAGGING==2
  {//
    int kk,ii;
    sprintf(rdir,"%s",rdir0);
    if((pdir=opendir(rdir))==NULL) {sprintf(cmd,"mkdir %s",rdir);system(cmd);} else closedir(pdir);
    sprintf(rdir,"%s/tmp",rdir0);
    if((pdir=opendir(rdir))==NULL) {sprintf(cmd,"mkdir %s",rdir);system(cmd);} else closedir(pdir);
    sprintf(rdir,"%s/net",rdir0);
    if((pdir=opendir(rdir))==NULL) {sprintf(cmd,"mkdir %s",rdir);system(cmd);} else closedir(pdir);

    if(BAGGING){//ensemble of all resampled-training folds (datasets)
      if(ibmode) sprintf(cmd2,"./meanpred %s ry:%d:%d:%.7e:%.7e:%.7e LDm:%d nop:%d",ibmode,r1,r2,r3,ymin,ymax,LDmode,nop);//r3桁で四捨五入
      else if(intmode) sprintf(cmd2,"./meanpred int");
      else sprintf(cmd2,"./meanpred");
    }
    int *bayes;
    char *fn_bayes;
    //    char *fn_bayesopt;
    //    int n_bayesopt=0;
    int bused=0;
    if(Bayes>0){
      bayes=(int*)malloc(sizeof(int)*nFolds);
      if(BAGGING==BAGwithVal){
	fn_bayes=(char*)malloc(sizeof(char)*fnsize);
	//	fn_bayesopt=(char*)malloc(sizeof(char)*fnsize);
	fn_bayes[0]=0;
	FILE *fp;
	if((fp=fopen("tmp/bayes.dat","r"))==NULL){
	  fprintf(stderr,"#Could not open 'tmp/bayes.dat'. Do Bagging with bg:/dev/null, first.\n");
	  exit(1);
	}
	char buff[80];
	for(j=0;j<nFolds;j++){
	  fgets(buff,80,fp);
	  sscanf(buff,"%d",&bayes[j]);
	}
	fclose(fp);
	sprintf(fb,"%s+%s+m%da%gb%d:%ds%dj%d",fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,t_boost,Boost,seed,bayes[0]);//
	sprintf(fn_bayes,  "%s/tmp/%sk%dN%d-%d:%denspred.dat"  ,rdir0,fb,K,NN[0],NN[nens-1],NStep);
      }
    }//closing of if(Bayes>0)
    for(j=0;j<nFolds;j++){//for all resampled-training datasets
      {//ensemble of the CAN2s with different units for each training fold
	if(ibmode) sprintf(cmd1,"./meanpred %s ry:%d:%d:%.7e:%.7e:%.7e LDm:%d nop:%d",ibmode,r1,r2,r3,ymin,ymax,LDmode,nop);//r3桁で四捨五入
	else if(intmode) sprintf(cmd1,"./meanpred int");
	else sprintf(cmd1,"./meanpred");
	msetrain[j]=n_msetrainsum=0;
      }

      if(Task==TIMESERIES){
	sprintf(fb,"%s+t%d-%d:%d-%d+s%d+",fnbody(fn_givendata,pr),tr0,tr1,tp0,tp1,seed);//20191210
	sprintf(fb0,"%s+s%d+",fnbody(fn_givendata,pr),seed);//
	// tradisional
	//	sprintf(fb,"%s+t%d-%d:%d-%d+",fnbody(fn_givendata,pr),tr0,tr1,tp0,tp1);
	//	sprintf(fb0,"%s+",fnbody(fn_givendata,pr));//
	//
	//	sprintf(fb,"%s+%d-%d+",fnbody(fn_givendata,pr),t1,t2);
	//	sprintf(fb,"%s+t%d:%d:%d+",fnbody(fn_givendata,pr),t1,t2,t0);
      }
      else if(BAGGING){
	sprintf(fb,"%s+%s+m%da%gb%d:%ds%dj%d",fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,t_boost,Boost,seed,j);//
#ifdef NEWp0
	sprintf(fb0,"%s+m%da%gb%d:%ds%dj%d",fnbody(fn_givendata,pr),method,alpha,t_boost,Boost,seed,j);//
#endif
      }
      else{//nobagging
	sprintf(fb,"%s+%s+s%dj%d",fnbody(fn_givendata,pr),fnbody(argv[2],pr1),seed,j);//argv[2]=method//argv[3]=nCells:nens
#ifdef NEWp0
	sprintf(fb0,"%s+s%dj%d",fnbody(fn_givendata,pr),seed,j);//argv[2]=method//argv[3]=nCells:nens
#endif
      }
      for(kk=0;kk<strlen(fb);kk++) if(fb[kk]=='/') fb[kk]='@';
      sprintf(fn_enspred,  "%s/tmp/%sk%dN%d-%d:%denspred.dat"  ,rdir0,fb,K,NN[0],NN[nens-1],NStep);
#ifdef NEWp0
      sprintf(fn_enspred0, "%s/tmp/%sk%dN%d-%d:%denspred0.dat" ,rdir0,fb0,K,NN[0],NN[nens-1],NStep);
#else
      sprintf(fn_enspred0, "%s/tmp/%sk%dN%d-%d:%denspred0.dat" ,rdir0,fb,K,NN[0],NN[nens-1],NStep);
#endif
      sprintf(fn_ensis,    "%s/tmp/%sk%dN%d-%d:%densis.dat"    ,rdir0,fb,K,NN[0],NN[nens-1],NStep);
      //      sprintf(fn_targetens,"%s/tmp/%starget.dat"   ,rdir0,fb);
      sprintf(fn_target,"%s/tmp/%starget.dat"   ,rdir0,fb);
      
      if(((fp=fopen(fn_enspred,"r"))!=NULL)&&((fp1=fopen(fn_ensis,"r"))!=NULL)&&((fp0=fopen(fn_enspred0,"r"))!=NULL)){
	fclose(fp);fclose(fp1);fclose(fp0);
	fprintf(stderr,"'%s' exists. Omit recalc.\n",fn_enspred);
	sprintf(cmd,"cp %s predict+.dat",fn_enspred);system(cmd);
      }////use previous fn_pred 
      else{//id1 fn_enspred 
	//insert preparing training dataset here ?
	//	if(chkomit && j==0) setfntraintest(fntrain,fntest,seed,fn_givendata,argv);
	if(chkomit!=0 && fntrain==NULL){
	  setfntraintest(seed,fn_givendata,argv,fupdate);
#ifdef CHKTIME
	  fprintf(stderr,"laptime: after setfntraintest:%f[sec]\n",mytimer_lap());
#endif
	}
	//	if(chkomit!=0 && fntrain==NULL) setfntraintest(fntrain,fntest,seed,fn_givendata,argv);
	for(ii=0;ii<nens;ii++){
	  _nCells=NN[ii];
	  //	  for(kk=0;kk<strlen(fb);kk++) if(fb[kk]=='/') fb[kk]='@';
	  sprintf(fn_pred,  "%s/tmp/%sk%dN%dpred.dat"  ,rdir0,fb,K,_nCells);
#ifdef NEWp0
	  sprintf(fn_pred0, "%s/tmp/%sk%dN%dpred0.dat" ,rdir0,fb0,K,_nCells);
#else
	  sprintf(fn_pred0, "%s/tmp/%sk%dN%dpred0.dat" ,rdir0,fb,K,_nCells);
#endif
	  sprintf(fn_is,    "%s/tmp/%sk%dN%dis.dat"    ,rdir0,fb,K,_nCells);

	  if(((fp=fopen(fn_pred,"r"))!=NULL)&&((fp1=fopen(fn_is,"r"))!=NULL)&&((fp0=fopen(fn_pred0,"r"))!=NULL)){
	    fclose(fp);fclose(fp1);fclose(fp0);
	    fprintf(stderr,"'%s' exists. Omit recalc.\n",fn_pred);
	  }////use previous fn_pred 
	  else{//making new fn_pred
	    net_printf(strx1("'%s' getting created.\n",fn_pred));
	    char fn_net[fnsize];
	    int fn_net_exist=0;
	    sprintf(fn_net,"result-ensrs2ge/net/%s+m%da%gb%d:%ds%dj%dN%dK%d.net",fnbody_givendata,method,alpha,t_boost,Boost,seed,j,_nCells,K);//
	    //	    fprintf(stderr,"@@@@@@@@kurocheck fn_net=%s@@@\n",fn_net);//only check
	    if((fp=fopen(fn_net,"r"))!=NULL){
	      fclose(fp);
	      fn_net_exist=1;
	    } 
	    fp=fopen("param.dat","w+");
	    if(Task==TIMESERIES){
	      //	      nTestData[j]=t2-t1;
	      nTestData[j]=tp1-tp0;
	      fprintf(fp,"0        #0:time-series,1:function,3:ijcnn04,4:range-data           \n");
	      fprintf(fp,"%d 0     #dimensionality k1, k2				       \n",K);
	      fprintf(fp,"%s						       	        \n",fntrain[j]);
	      if(BAGGING==NoBAG)  fprintf(fp,"%d-%d:%d-%d			        \n",tr0,tr1,tp0,tp1);
	      else{ fprintf(stderr,"Bagging for Time Series is not used now.\n");}
	      if(ymax-ymin<1e-20 && ymax1-ymin1<1e-20){
		fprintf(fp,"0 0 0 0   #normalization-of-output ymin0 ymax0 ->ymin1 ymax1	        \n"); 
	      }
	      else{
		fprintf(fp,"%e %e %e %e %e:%e:%e:%e #y-normalization[ymin0:ymax0]->[ymin1:ymax1]\n",ymin,ymax,ymin1,ymax1,y0m,y1m,y0M,y1M); 
		fprintf(fp,"%d %d %e #r1 r2 					        \n",r1,r2,r3);
	      }
	    }
	    else if(Task==REGRESSION){
	      fprintf(fp,"1          #0:time-series,1:function,3:ijcnn04,4:range-data,No.1\n");
	      fprintf(fp,"%d 0       #dimensionality k1 k2					\n",K);
	      fprintf(fp,"%s						       	        \n",fntrain[j]);
	      if(BAGGING==NoBAG)  fprintf(fp,"%s				        \n",fntest[j]);
	      else if(BAGGING==BAGwithVal) fprintf(fp,"%s 				\n",fn_bagging);//see below fn_givendata
	      else if(BAGGING==BAGwithoutVal) fprintf(fp,"%s 				\n",fn_givendata);//bagging test
	      if(ymax-ymin<1e-20 && ymax1-ymin1<1e-20){
		fprintf(fp,"%e %e %e %e %e:%e:%e:%e   #y-normalization[ymin0:ymax0]->[ymin1:ymax1]\n",ymin,ymax,ymin1,ymax1,y0m,y1m,y0M,y1M); 
	      }
	      else
		{
		  fprintf(fp,"%e %e %e %e %e:%e:%e:%e   #y-normalization[ymin0:ymax0]->[ymin1:ymax1]\n",ymin,ymax,ymin1,ymax1,y0m,y1m,y0M,y1M); 
		  for(i=0;i<K;i++) fprintf(fp,"%e %e %e %e %e:%e:%e:%e #x%d-normalization xmin0:xmax0->xmin:xmax\n",xmin[i],xmax[i],xmin1[i],xmax1[i],x0m[i],x1M[i],x0M[i],x1M[i],i);
		  fprintf(fp,"%d %d %e #r1 r2 					        \n",r1,r2,r3);
		}
	    }//if(Task==REGRESSION
	    if(fn_net_exist){
	      fprintf(stderr,"'%s' exists. Omit re-learning.\n",fn_net);
	      fprintf(fp,"nl         net-load		                \n");
	      fprintf(fp,"%s\n",fn_net);
	      if(r0!=0) fprintf(fp,"r0\n");
	    }
	    else{
	      net_printf(strx1("'%s' being created.\n",fn_net));
	      fprintf(fp,"in          init-net			                \n");
	      fprintf(fp,"%d          number-of-cells				          	\n",_nCells);
	      fprintf(fp,"%d           n_compare					\n",NC);
	      fprintf(fp,"%e %d %d     v_thresh vmin vmin2	         		\n",vt,vm,vm2);
	      fprintf(fp,"%e           v_ratio				                \n",vr);
	      fprintf(fp,"%e %e %e   width,gamma  #window width                         \n",w,gamma0,entropy_thresh);
	      fprintf(fp,"ex          execution						\n");
	      //	      fprintf(fp,"ex          実行						\n");
	      fprintf(fp,"1 %e %e   #<0:online,1:batch>, <gamma0>, <entropy_thresh>     \n",gamma0,entropy_thresh);	  
	      fprintf(fp,"%d        i_times number of learning iterations	        \n",T);
	      //      fprintf(fp,"%d        i_times 学習回数				        \n",T);
	      fprintf(fp,"%d 50 350    number-of-display, rot_x, rot_z			       \n",T);
	      //     fprintf(fp,"%d 50 350    表示回数, rot_x, rot_z			       \n",T);
	    }
	    //if(1==1){ fprintf(fp,"ssp_\n");}//ssp without display
	    if(Task==TIMESERIES){
	      if(ssp==1) fprintf(fp,"ssp_\n");//no effect 091027?
	      else fprintf(fp,"msp_%g %g\n",err4terminate,err4propagate);//20150218
	      //	      else fprintf(fp,"msp_\n");
	    }
	    else{
	      if(fn_net_exist){fprintf(fp,"ssp_\n");}//ssp without display
	    }
	    if(BAGGING==NoBAG){//??	    if(!BAGGING){
	      //	      fprintf(fp,"!cp mse.dat mse_.dat\n");//for saving
	      if(DISP==0){
		fprintf(fp,"ssp0_      #prediction of training data -> predict0.dat    \n");//for Bayesian??
	      }
	      else{
		fprintf(fp,"ssp0      #prediction of training data -> predict0.dat    \n");//for Bayesian??
	      }
	      //	      fprintf(fp,"ssp0_      #訓練データの予測 > predict0.dat    \n");//for Bayesian??
	      if(lcom) fprintf(fp,"%s            \n",lcom);//for last command
	    }
	    fprintf(fp,"qu         \n");
	    fclose(fp);
            if(nop==0){fprintf(stdout,"--- param.dat ---\n");system("cat param.dat");fprintf(stdout,"--- param.dat ---\n");}
	    if(fabs(BIAS-1.)>1e-10) {
              char _cmd[80]; sprintf(_cmd,"./can2 BIAS:%g LINESIZE:%d seed:%u Tpinv:%d nop:%d DISP:%d < param.dat",BIAS,LINESIZE,seed,Tpinv,nop,DISP);
              net_printf(strx1("#%s\n",_cmd));
	      system(_cmd);//>>predict0.dat ??
	    }
	    else{
              char _cmd[80]; sprintf(_cmd,"./can2 LINESIZE:%d seed:%u Tpinv:%d nop:%d DISP:%d< param.dat",LINESIZE,seed,Tpinv,nop,DISP);
	      system(_cmd);//--err
              net_printf(strx1("#%s\n",_cmd));
	    }
#ifdef CHKTIME
	    fprintf(stderr,"laptime: after can2:%f[sec]\n",mytimer_lap());
#endif
	    //	    if(fn_net_exist==0){
	    if(fn_net_exist==0){
	      sprintf(cmd,"cp last.net %s",fn_net);system(cmd);//if(fn_net_exist==0) sprintf(cmd,"cp bestssp.net %s",fn_net);system(cmd); 
	    }
	    if(method!=-1 || Task==TIMESERIES){
	      fp=fopen("is.dat","w+");
	      if(BAGGING==NoBAG){
		FILE *fpmse=fopen("mse0.dat","r");fscanf(fpmse,"%lf",&msetrain[j]);fclose(fpmse);
		fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.,nTrainData,msetrain[j],_nCells,_nCells,nTestData[j],0.,0.);
	      }
	      else if(BAGGING==BAGwithVal){
		fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.,nTrainData,0.,_nCells,_nCells,nValData,0.,0.);
	      }
	      else if(BAGGING==BAGwithoutVal){
		//		fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.,nTrainData,0.,_nCells,_nCells,nTrainData,0.,0.);
		fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.,nTrainData,0.,_nCells,_nCells,nGivenData,0.,0.);//070221
	      }
	      fclose(fp);
	    }
	    sprintf(cmd,"cp predict.dat %s",fn_pred);system(cmd);
	    sprintf(cmd,"cp is.dat %s",fn_is);system(cmd);

	    if(BAGGING==BAGwithVal){//for predict0.dat for bagging. training loss:  Lvar0, Lib, Lemp for the gn_givendata
	      fp=fopen("param.dat","w+");
	      fprintf(fp,"1          #0:time-series,1:function,3:ijcnn04,4:range-data,No.2\n");
	      //	      fprintf(fp,"1          #0:時系列,1:関数,3:ijcnn04,4:距離データ,No.2\n");
	      fprintf(fp,"%d 0     #dimensionality k1 k2				        \n",K);
	      //	      fprintf(fp,"%d 0     #次元 k1 k2				        \n",K);
	      //dummy
	      fprintf(fp,"%s 				\n",fn_bagging);//see below fn_givendata
	      //if(BAGGING==NoBAG)  fprintf(fp,"%s				        \n",fntest[j]);
	      //else if(BAGGING==BAGwithVal) fprintf(fp,"%s 				\n",fn_bagging);//see below fn_givendata
	      //else if(BAGGING==BAGwithoutVal) fprintf(fp,"%s 				\n",fn_givendata);//bagging test
	      fprintf(fp,"%s \n",fn_givendata);
	      //??090810	      fprintf(fp,"%s						       	        \n","/dev/null");//??
	      if(ymax-ymin<1e-20){
		fprintf(fp,"0 0 0 0   #normalization-of-output ymin0 ymax0 -> ymin1 ymax1	        \n"); 
		//		fprintf(fp,"0 0 0 0   出力の正規化 ymin0 ymax0 →ymin1 ymax1	        \n"); 
	      }
	      else{
		fprintf(fp,"%e %e %e %e    #normalization-of-output ymin0 ymax0 -> ymin1 ymax1	        \n",ymin,ymax,ymin1,ymax1); 
		//		fprintf(fp,"%e %e %e %e    出力の正規化 ymin0 ymax0 →ymin1 ymax1	        \n",ymin,ymax,ymin1,ymax1); 
		for(i=0;i<K;i++) fprintf(fp,"%e %e %e %e #input-%d-normalization xmin0 xmax0->xmin xmax\n",xmin[i],xmax[i],xmin1[i],xmax1[i],i);
		//for(i=0;i<K;i++) fprintf(fp,"%e %e %e %e #入力%dの正規化 xmin0 xmax0→xmin xmax\n",xmin[i],xmax[i],xmin1[i],xmax1[i],i);
		fprintf(fp,"%d %d %e #r1 r2 					        \n",r1,r2,r3);
	      }
	      //	      if(fabs(BIAS-1.)>1e-10) fprintf(fp,"BIAS:%e		        \n",BIAS);
	      fprintf(fp,"nl         net-load			                \n");
	      //	      fprintf(fp,"nl         ネットのロード			                \n");
	      fprintf(fp,"%s\n",fn_net);//	      fprintf(fp,"last.net\n");
	      fprintf(fp,"ssp0_      #prediction of training data -> predict0.dat    \n");//for Bayesian??
	      //	      fprintf(fp,"ssp0_      #訓練データの予測 > predict0.dat    \n");//for Bayesian??
	      if(lcom) fprintf(fp,"%s            \n",lcom);//for last command
	      fprintf(fp,"qu         \n");
	      fclose(fp);
              if(nop==0){
                printf("--- param.dat ---\n");
                system("cat param.dat");
                printf("--- param.dat ---\n");
              }
	      if(fabs(BIAS-1.)>1e-10) {
		char _cmd[40]; sprintf(_cmd,"./can2 BIAS:%g LINESIZE:%d seed:%d Tpinv:%d nop:%d DISP:%d< param.dat",BIAS,LINESIZE,seed,Tpinv,nop,DISP);system(_cmd);//>>predict0.dat ??
	      }
	      else{
		char _cmd[40]; sprintf(_cmd,"./can2 LINESIZE:%d seed:%d Tpinv:%d nop:%d DISP:%d< param.dat",LINESIZE,seed,Tpinv,nop,DISP);system(_cmd);//
		//		system("./can2 LINESIZE:%d < param.dat");//>>predict0.dat ??
	      }
	    }
#ifdef CHKTIME
	    fprintf(stderr,"laptime: after 2nd can2:%f[sec]\n",mytimer_lap());
#endif
	    if(BAGGING!=BAGwithoutVal) {sprintf(cmd,"cp predict0.dat %s",fn_pred0);system(cmd);}//???
	    else {sprintf(cmd,"cp /dev/null %s",fn_pred0);system(cmd);}
	    //	    if(1==1 || BAGGING!=BAGwithoutVal) {sprintf(cmd,"cp predict0.dat %s",fn_pred0);system(cmd);}//???
	    //	    else {sprintf(cmd,"cp /dev/null %s",fn_pred0);system(cmd);}
	  }//closing making fn_pred 
	  /////////////////////////////////////////////////////////////////////
	  {
	    int lc=strlen(cmd1),la=strlen(fn_pred);//char *oldcmd1=cmd1;
	    if(cmd1size<lc+la){
	      cmd1size=lc+10*la;
	      cmd1=realloc(cmd1,cmd1size);
	      //	      if(cmd1!=oldcmd1) {int j; for(j=0;j<lc;j++) cmd1[j]=oldcmd1[j];}
	    }
	    sprintf(&cmd1[lc]," %s",fn_pred);//combine meanpred 
	  }
	  if(BAGGING==NoBAG){
	    {//	if(ERR==MSE){
	      double _mseis,_mse,_msetrain,_msetrainis; 
	      int _nCells1;
	      fp=fopen(fn_is,"r");
	      fscanf(fp,"%d %lf %d %lf %d %d %d %lf %lf",&is[j],&_mseis,&n_msetrain[j],&_msetrain,&n_cells2[j],&_nCells1,&n_mse[j],&_mse,&_msetrainis);
	      msetrain[j]+=(_msetrain*n_msetrain[j]);//rational
	      n_msetrainsum+=n_msetrain[j];
	      fclose(fp);
	    }
	  }
	  if(ii==(int)(nens/2.)){//center???
	    sprintf(cmd,"cp %s %s",fn_is,fn_ensis);system(cmd);//??
	    sprintf(cmd,"cp %s %s",fn_pred0,fn_enspred0);system(cmd);//??
	  }
	}//closing  for(ii=0;ii<nens;ii++){
	{//ensemble,mse.dat,loss+.dat for nFolds
	  msetrain[j]/=(n_msetrainsum);//Modified at 060517
	  n_msetrain[j]=n_msetrainsum;//modified at 080707??
	  net_printf(strx1("Doing '%s'\n",cmd1));system(cmd1);//exec meanpred for ensemble
	  sprintf(cmd,"cp predict+.dat %s",fn_enspred);system(cmd);
	  sprintf(cmd,"cat %s | awk '{print $%d}' > %s",fntest[j],K+1,fn_target);system(cmd);
	  sprintf(&cmd[strlen(cmd)]," ;./pred2y_ts predict+.dat %s > /dev/null",fn_target);system(cmd);
	  fp=fopen("mse.dat","r"); fscanf(fp,"%lf",&mse[j]);fclose(fp);
	  fp=fopen("loss+.dat","r"); fscanf(fp,"%lf",&Lhat[j]);fclose(fp);
	}
      }//closing else id1 fn_enspred 
      //if((BAGGING==BAGwithoutVal||BAGGING==NoBAG) && Boost!=GbBoost){//for training(empirical) and prediction(generalization) error
      if(BAGGING==BAGwithoutVal && Boost!=GbBoost){//for training(empirical) and prediction(generalization) error
	//if(BAGGING==BAGwithoutVal){//for training(empirical) and prediction(generalization) error
	//if(BAGGING==BAGwithoutVal && Boost==EmBoost){//for training(empirical) and prediction(generalization) error
	//Boost==EmBoost && Boost == NoBoost for manipulate training and test predictions
	int _idum;
	fp=fopen(fn_enspred,"r");
	for(n=0;n<nGivenData;n++) {
	  fgets(line1,LINESIZE,fp);
	  sscanf(line1,"%lf%d%lf",&yp0[j*nGivenData+n],&_idum,&yt[n]);
	}
	fclose(fp);
      }
      int lc=strlen(cmd2),la=strlen(fn_enspred);//char *oldcmd2=cmd2;
      if(cmd2size<lc+la){
	cmd2size=lc+10*la;
	cmd2=realloc(cmd2,cmd2size);
	//	if(cmd2!=oldcmd2) {int j; for(j=0;j<lc;j++) cmd2[j]=oldcmd2[j];}
      }
	    
      if(Bayes>0 && BAGGING==BAGwithVal){
	if(bayes[j]==j){
	  sprintf(&cmd2[lc]," %s",fn_enspred);
	  sprintf(fn_bayes,"%s",fn_enspred);
	  bused++;
	}
	else{
	  if(strlen(fn_bayes)>=1) sprintf(&cmd2[lc]," %s",fn_bayes);
	  //	  else n_bayesopt++;
	}
	//if(bayes[0]==j) sprintf(fn_bayesopt,"%s",fn_enspred);
      }
      else sprintf(&cmd2[lc]," %s",fn_enspred);//combine meanpred 
    }//end of for(j=0;j<nFolds;j++) for learning CAN2
//    if(Bayes>0 && BAGGING==BAGwithVal){//compensate until nFolds
//      for(j=0;j<n_bayesopt;j++){
//	int lc=strlen(cmd2),la=strlen(fn_bayesopt);//
//	if(cmd2size<lc+la){
//	  cmd2size=lc+10*la;
//	  cmd2=realloc(cmd2,cmd2size);
//	}
//	sprintf(&cmd2[lc]," %s",fn_bayesopt);
//      }
//    }
//////////////////////////////////////
#ifdef CHKTIME
    fprintf(stderr,"laptime: before averaging:%f[sec]\n",mytimer_lap());
#endif
    //below are for averaging
    if(Boost==NoBoost){//no boosting
      if(BAGGING==BAGwithVal){
	if(nFolds>100){
	  FILE *fpm=fopen("tmp/meanpred.sh","w");
	  fprintf(fpm,"%s\n",cmd2);
	  fclose(fpm);
	  sprintf(cmd2,"meanpred par:tmp/meanpred.sh");
	}
	net_printf(strx1("Doing cmd2='%s'\n",cmd2));
	system(cmd2);//exec meanpred for ensemble
	//	fprintf(stderr,"#ensrs with b%d/%d,Bt%gBs%d\n",nFolds,bused,BayesLambdaL,Bayesseed);
	if(lossall){
	  FILE *fp=fopen("loss.dat","r");
	  char buff[2][256];
	  fgets(buff[0],256,fp);
	  fgets(buff[1],256,fp);
	  fclose(fp);
	  double Lbayesob;int nobmin;
	  if((fp=fopen("tmp/bayesL.dat","r"))!=NULL){fscanf(fp,"%lf%d",&Lbayesob,&nobmin);fclose(fp);}
	  else {Lbayesob=nobmin=-1;}
	  fp=fopen("lossall.dat","a+");
	  buff[0][strlen(buff[0])-1]=0; buff[1][strlen(buff[1])-1]=0;
	  fprintf(fp,"%s %d %g %d %d %g %d %s,bu,Bt,Bs,bs,Lob,nobmin\n",buff[0],bused,BayesLambdaL,Bayesseed,seed,Lbayesob,nobmin,buff[1]);
	  fclose(fp);
	  fprintf(stderr,"#'rm lossall.dat' if necessary.\n");
	}
	//	fprintf(stderr,"#b%dB%g,bused%d\n",nFolds,BayesLambdaL,bused);
    //	printf("#%s",fn_bagging);
    //	printf("#%d",K+1);
    //	printf("#%s",fn_target);
    printf("cat %s | awk '{print $%d}' > %s",fn_bagging,K+1,fn_target);
    //    printf("#cmd=%s.\n",cmd);//for check
	sprintf(cmd,"cat %s | awk '{print $%d}' > %s",fn_bagging,K+1,fn_target);
    system(cmd);
	sprintf(&cmd[strlen(cmd)]," ;./pred2y_ts predict+.dat %s > /dev/null",fn_target);system(cmd);
      }//end of if(BAGGING==BAGwithVal)
      else if(BAGGING==BAGwithoutVal){
	double LD=0,Ltst,Lvar,Lvarval=0,Lvar0,Lib,Lemp,H0,L0,Ltr;
	int nfiles,num;
	double *yptest=(double*)malloc(sizeof(double)*nGivenData);
	double *yptrain=(double*)malloc(sizeof(double)*nGivenData);
	double *yptotal=(double*)malloc(sizeof(double)*nGivenData);
	int *ntrain=(int*)malloc(sizeof(int)*nGivenData);
	int *ntest=(int*)malloc(sizeof(int)*nGivenData);
	int ntrainall,ntestall;
	double err2;
	double UC,ucn;//uncertainty
	double *uc=(double*)malloc(sizeof(double)*nGivenData);
	int nGivenData1=0,nGivenData0=0;
	double _ntrain=0,_ntest=0;
	double Lbayes;
	int nbayes=0;
	double Lbmean;
	Ltst=Lib=Lemp=0;
	if(Bayes>0){
#if RAND == ZMTRAND
	  InitMt((unsigned long)BayesSeed);
#elif RAND == DRAND48
          seed48((unsigned short*)&BayesSeed);
#endif
	  double *Lb=(double*)malloc(sizeof(double)*(nFolds+1));
	  double *Sb=(double*)malloc(sizeof(double)*(nFolds+1));
	  double _nob,_nib;
	  
#if RAND == ZMTRAND
	  InitMt((unsigned long)Bayesseed);
#elif RAND == MYRAND
	  _randn=Bayesseed;
	  
#elif RAND == DRAND48
	  srand48((long int)Bayesseed);
#elif RAND == RANDOM
	  srandom(Bayesseed);
#endif
	  if(Bayes==3){//use L=MSEmean, 
	    double LambdaL;
	    _nob=_nib=0;
	    for(j=0;j<nFolds;j++){
	      Lb[j]=0;
	      int nob=0,nib=0;
	      for(n=0;n<nGivenData;n++){
		if(BayesUseAllData || nc[j*nGivenData+n]!=1){//Lemp(for BayesUseAllData==1) or Lob(not in Dtrain)
		  //if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  double err=yp0[j*nGivenData+n]-yt[n];
		  Lb[j]+=(err*err);
		}
		else nib++;
	      }
	      if(nob==0){
		Lb[j]=1e20;
	      }
	      if(nob>0) Lb[j]/=nob;
	      _nob+=nob;
	      _nib+=nib;
	    }
	    _nob/=nGivenData;
	    _nib/=nGivenData;
	    Lbmean=0;//double Lbm=1e20,LbM=-1;
	    int nLbmean=0;
#define metropolisNEW //use minimum Lb[jm] for the initial net
	    //#undef metropolisNEW //use minimum Lb[jm] for the initial net
#ifdef metropolisNEW 
	    int jm=0;
#endif
	    for(j=0;j<nFolds;j++){
	      if(Lb[j]<1e20){
		Lbmean+=Lb[j];
		nLbmean++;
	      }
#ifdef metropolisNEW
	      if(Lb[j]<Lb[jm]) jm=j;
#endif
		//if(Lbm>Lb[j]) Lbm=Lb[j]; if(LbM<Lb[j] && Lb[j]<1e20) LbM=Lb[j];
	    }
	    Lbmean/=nLbmean;
	    if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./Lbmean;
	    else LambdaL=1./Lbmean;
//	    if(LbM>Lbm){
//	      if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./(LbM-Lbm);
//	      else LambdaL=1./(LbM-Lbm);
//	    }
//	    else LambdaL=0;
	    fprintf(stderr,"##%sLmbdL%g(Lbm%g)nob%.0fnib%.0f",argv[iBayes],LambdaL,Lbmean,_nob,_nib);
	    //	    fprintf(stderr,"#%sLambdaL%g=%g/2/%.2enob%.0fnib%.0f",argv[iBayes],LambdaL,BayesLambdaL,Lbmean,_nob,_nib);
	    //	    fprintf(stderr,"#%s. LambdaL=%g=%g/(2(%g-%g)). ",argv[iBayes],BayesLambdaL,LambdaL,LbM,Lbm);
	    FILE *fpb=fopen("tmp/bayes.dat","w");
	    int b1;// b1=j=0;fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]); bayes[j]=b1;
	    Lbmean=0;
	    nLbmean=0;

#ifdef metropolisNEW
	    b1=jm;
	    for(j=0;j<nFolds;j++){
	      if(Lb[j]<Lb[b1] || NextUnif()<exp((Lb[b1]-Lb[j])*LambdaL)){
		if(Lb[j]<1e20) {
		  b1=j;nbayes++;
		}
	      }
	      fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]);
	      bayes[j]=b1;
	      if(Lb[b1]<1e20){ Lbmean+=Lb[b1];nLbmean++;}
	    }
#else
	    for(j=0;j<nFolds;j++){
	      if(j==0){
		b1=j;
	      }
	      //	      else if(Lb[j]<Lb[j-1] || NextUnif()<exp((Lb[j-1]-Lb[j])*LambdaL)){
	      else if(Lb[j]<Lb[b1] || NextUnif()<exp((Lb[b1]-Lb[j])*LambdaL)){
		if(Lb[j]<1e20){
		  b1=j; nbayes++;
		}
	      }
	      fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]);
	      bayes[j]=b1;
	      if(Lb[b1]<1e20){ Lbmean+=Lb[b1];nLbmean++;}
	    }
#endif
	    fclose(fpb);
	    Lbmean/=nLbmean;
	  }//end of if(Bayes==3)
	  else if(Bayes==1){//use L=MSEmax-MSEmin, 
	    double LambdaL;
	    _nob=_nib=0;
	    for(j=0;j<nFolds;j++){
	      Lb[j]=0;
	      int nob=0,nib=0;
	      for(n=0;n<nGivenData;n++){
		if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  double err=yp0[j*nGivenData+n]-yt[n];
		  Lb[j]+=(err*err);
		}
		else nib++;
	      }
	      if(nob==0){
		Lb[j]=1e20;
	      }
	      if(nob>0) Lb[j]/=nob;
	      _nob+=nob;
	      _nib+=nib;
	    }
	    _nob/=nGivenData;
	    _nib/=nGivenData;
	    double Lbm=1e20,LbM=-1;
	    for(j=0;j<nFolds;j++){
	      if(Lbm>Lb[j]) Lbm=Lb[j];
	      if(LbM<Lb[j] && Lb[j]<1e20) LbM=Lb[j];
	    }
	    if(LbM>Lbm){
	      if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./(LbM-Lbm);
	      else LambdaL=1./(LbM-Lbm);
	    }
	    else LambdaL=0;
	    fprintf(stderr,"##%sLambdaL%g=%g/(2(%g-%g))nob%.0fnib%.0f",argv[iBayes],LambdaL,BayesLambdaL,LbM,Lbm,_nob,_nib);
	    //	    fprintf(stderr,"#%s. LambdaL=%g=%g/(2(%g-%g)). ",argv[iBayes],BayesLambdaL,LambdaL,LbM,Lbm);
	    FILE *fpb=fopen("tmp/bayes.dat","w");
	    int b1; b1=j=0;fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]); bayes[j]=j;
	    for(j=1;j<nFolds;j++){
	      if(Lb[j]<Lb[j-1] || NextUnif()<exp((Lb[j-1]-Lb[j])*LambdaL))	b1=j;
	      fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]);
	      bayes[j]=b1;
	    }
	    fclose(fpb);
	  }//end of if(Bayes==1)
	  else if(Bayes==4){//use L=MSEmean, and S=
	    //skew
	    double *_yp=(double*)malloc(sizeof(double)*nGivenData);
	    for(n=0;n<nGivenData;n++){
	      double _yob=0,_yib=0;
	      int nob=0,nib=0;
	      for(j=0;j<nFolds;j++){
		if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  _yob+=yp0[j*nGivenData+n];
		}
		else{
		  nib++;
		  _yib+=yp0[j*nGivenData+n];
		}
	      }
	      if(nob>0) _yp[n]=_yob/nob;
	      else _yp[n]=_yib/nib;//?
	    }
	    double LambdaL,LambdaS;
	    _nob=_nib=0;
	    double skew=0,sigma=0;
	    for(j=0;j<nFolds;j++){
	      Lb[j]=0;
	      int nob=0,nib=0;
	      for(n=0;n<nGivenData;n++){
		if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  double err=yp0[j*nGivenData+n]-yt[n];
		  Lb[j]+=(err*err);
		  double err1=yp0[j*nGivenData+n]-_yp[n];
		  double err2=err1*err1;
		  double err3=err2*err1;
		  sigma+=err2;
		  skew+=err3;
		}
		else nib++;
	      }
	      if(nob>0){
		Lb[j]/=nob;
		sigma=sqrt(sigma/nob);
		Sb[j]=(skew/nob)/sigma/sigma/sigma;
	      }
	      else { //if(nob==0){
		Lb[j]=Sb[j]=1e20;
	      }
	      _nob+=nob;
	      _nib+=nib;
	    }
	    free(_yp);

	    _nob/=nGivenData;
	    _nib/=nGivenData;
	    Lbmean=0;//double Lbm=1e20,LbM=-1;
	    int nLbmean=0;
#define metropolisNEW
#ifdef metropolisNEW
	    //use minimum Lb[jm] for the initial net
	    int jm=0;
#endif
	    for(j=0;j<nFolds;j++){
	      if(Lb[j]<1e20){
		Lbmean+=Lb[j];
		nLbmean++;
	      }
#ifdef metropolisNEW
	      if(Lb[jm]>Lb[j]) jm=j;
#endif
		//if(Lbm>Lb[j]) Lbm=Lb[j]; if(LbM<Lb[j] && Lb[j]<1e20) LbM=Lb[j];
	    }
	    Lbmean/=nLbmean;
	    if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./Lbmean;
	    else LambdaL=1./Lbmean;
	    if(BayesLambdaS<1000) LambdaS=BayesLambdaS;
	    else LambdaS=0;

	    fprintf(stderr,"##%sLambdaL%gnob%.0fnib%.0f",argv[iBayes],LambdaL,_nob,_nib);
	    //	    fprintf(stderr,"#%sLambdaL%g=%g/2/%.2enob%.0fnib%.0f",argv[iBayes],LambdaL,BayesLambdaL,Lbmean,_nob,_nib);
	    //	    fprintf(stderr,"#%s. LambdaL=%g=%g/(2(%g-%g)). ",argv[iBayes],BayesLambdaL,LambdaL,LbM,Lbm);
	    FILE *fpb=fopen("tmp/bayes.dat","w");
	    int b1;// b1=j=0;fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]); bayes[j]=b1;
	    Lbmean=0;
	    nLbmean=0;
#ifdef metropolisNEW
	    b1=jm;
	    for(j=0;j<nFolds;j++){
	      double p=(Lb[b1]-Lb[j])*LambdaL+(Sb[b1]-Sb[j])*LambdaS;
	      if(p>0 || NextUnif()<exp(p)){
		//if(Lb[j]<Lb[b1] || NextUnif()<exp((Lb[b1]-Lb[j])*LambdaL+(Sb[b1]-Sb[j])*LambdaS)){
		if(Lb[j]<1e20) {
		  b1=j;nbayes++;
		}
	      }
	      fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]);
	      bayes[j]=b1;
	      if(Lb[b1]<1e20){ Lbmean+=Lb[b1];nLbmean++;}
	    }
#endif
	    fclose(fpb);
	    Lbmean/=nLbmean;
	    //	    fprintf(stderr,"Lb=%g",Lbmean);
	  }//endof if(Bayes==4)
	  else if(Bayes==2){//use variance of prediction
	    double LambdaL;
	    double *_yp=(double*)malloc(sizeof(double)*nGivenData);
	
	    for(n=0;n<nGivenData;n++){
	      double _yob=0,_yib=0;
	      int nob=0,nib=0;
	      for(j=0;j<nFolds;j++){
		if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  _yob+=yp0[j*nGivenData+n];
		}
		else{
		  nib++;
		  _yib+=yp0[j*nGivenData+n];
		}
		
	      }
	      if(nob>0)	_yp[n]=_yob/nob;
	      else _yp[n]=_yib/nib;//?
	    }
	    for(j=0;j<nFolds;j++){
	      int nob=0,nib=0;
	      Lb[j]=0;
	      for(n=0;n<nGivenData;n++){
		if(nc[j*nGivenData+n]!=1){//not in Dtrain
		  nob++;
		  double err=yp0[j*nGivenData+n]-_yp[n];
		  Lb[j]+=(err*err);
		}
		else nib++;
	      }
	      if(nob>0)	Lb[j]/=nob;
	      else Lb[j]=1e20;
	      _nob+=nob;
	      _nib+=nib;
	    }
	    _nob/=nGivenData;
	    _nib/=nGivenData;
#define LambdaLOrig
#ifndef LambdaLOrig
	    double _Lb=0;int n_Lb;
	    for(j=0;j<nFolds;j++){
	      if(Lb[j]<1e20) {
		_Lb+=Lb[j];
		n_Lb++;
	      }
	    }
	    _Lb/=n_Lb;
	    if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./_Lb;
	    fprintf(stderr,"##%sLambdaL%g=%g/2/%.2enob%.0fnib%.0f",argv[iBayes],LambdaL,BayesLambdaL,_Lb,_nob,_nib);
#else
	    double Lbm=1e20,LbM=-1;
	    for(j=0;j<nFolds;j++){
	      if(Lbm>Lb[j]) Lbm=Lb[j];
	      if(LbM<Lb[j] && Lb[j]<1e20) LbM=Lb[j];
	    }
	    if(LbM>Lbm){
	      if(BayesLambdaL>=0) LambdaL=BayesLambdaL/2./(LbM-Lbm);
	      else LambdaL=1./(LbM-Lbm);
	    }
	    else LambdaL=0;
	    fprintf(stderr,"##%sLambdaL%g=%g/2/(%.2e-%.2e)nob%.0fnib%.0f",argv[iBayes],LambdaL,BayesLambdaL,LbM,Lbm,_nob,_nib);
#endif
	    FILE *fpb=fopen("tmp/bayes.dat","w");
	    int b1; b1=j=0;fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]); bayes[j]=j;
	    for(j=1;j<nFolds;j++){
	      if(Lb[j]<Lb[j-1] || NextUnif()<exp((Lb[j-1]-Lb[j])*LambdaL)) b1=j;//metropolis
	      fprintf(fpb,"%d %d %.5e %.5e\n",b1,j,Lb[j],Lb[b1]);
	      bayes[j]=b1;
	    }
	    fclose(fpb);
	  }//endof if(Bayes==2)
	  free(Lb);
	  //endof all if(Bayes==*)
	  double Lob=0,Lib=0;Lemp=0;
	  int nob_=0,nib_=0,nobmin=nFolds,nemp_=0;
	  for(n=0;n<nGivenData;n++){
	    double yptest=0,yib=0,yemp=0;
	    int nob=0,nib=0,nemp=0;
	    int jj;for(jj=0;jj<nFolds;jj++){
	      j=bayes[jj];
	      double yp0jn=yp0[j*nGivenData+n];
	      if(nc[j*nGivenData+n]!=1){//not in Dtrain
		yptest+=yp0jn;
		nob++;
	      }
	      else{
		yib+=yp0jn;
		nib++;
	      }
	      yemp+=yp0jn;
	      nemp++;
	    }
	    if(nob<nobmin) nobmin=nob;
	    //	    if(nob>=nob_thresh){//accept stable elements
	    if(nob>0){//accept stable elements
	      nob_++;
	      yptest/=(nob);
	      double err=yptest-yt[n];
	      Lob+=(err*err);
	      fprintf(stderr,"Lob%g,n%d,nob=%d,nFolds=%d.OK\n",Lob,n,nob,nFolds);
	    }
	    else{//neglect unstable elements
	      //nob_++;
	      //Lob+=1e20;//???
	      fprintf(stderr,"Lob%g,n%d.nob=%d,nFolds=%d.NG\n",Lob,n,nob,nFolds);
	    }
	    if(nib>0){
	      nib_++;
	      yib/=(nib);
	      double err=yib-yt[n];
	      Lib+=(err*err);
	    }
	    {
	      nemp_++;
	      yemp/=(nemp);
	      double err=yemp-yt[n];
	      Lemp+=(err*err);
	    }
	  }//endof for(n=0;n<nGivenData;n++)
	  Lob/=nob_;
	  Lib/=nib_;
	  Lemp/=nemp_;
	  Lbayes=Lob;
	  fprintf(stderr,"#Bys:Lob%.3eLib%.3eLemp%.3eLbB%.3eenb%dnobm%d",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin);
	  //	  fprintf(stderr,"#Bys:Lob%gLib%gLemp%gLbB%genb%dnobm%d",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin);
	  //	  fprintf(stderr,"#Bys:Lob%gLib%gLemp%gLbB%genb%dnobm%d k%dN%darg:%s",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin,K,nCells-1,argv[2]);
	  if(lossall){
	    FILE *fp=fopen("lossall.dat","a+");
	    fprintf(fp,"#Bys:Lob%.3eLib%.2eLem%.2eLbB%.2eenb%dnob%d",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin);
	    //	    fprintf(fp,"#Bys:Lob%.3eLib%.3eLemp%.3eLbB%.3eenb%dnobm%d",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin);
	    //	    fprintf(fp,"#Bys:Lob%gLib%gLemp%gLbB%genb%dnobm%d",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin);
	    //	    fprintf(fp,"#Bys:Lob%gLib%gLemp%gLbB%genb%dnobm%d k%dN%darg:%s",Lob,Lib,Lemp,Lbmean-Lbayes,nbayes,nobmin,K,nCells-1,argv[2]);
	    fclose(fp);
	  }
//	}//if(Bayes>0)
//      }//else if(BAGGING==BAGwithoutVal){
//    }//    if(Boost==NoBoost){//no boosting
	  //	  fprintf(stderr,"#Bayes:Lob(nb%d>%d)Lib%.2eLbB%.2e",Lbayes,nbayes,nobmin,Lib,Lbmean-Lbayes);
	  //calc distortion of variance
	  double *_yp=(double*)malloc(sizeof(double)*nGivenData);
	  for(n=0;n<nGivenData;n++){
	    double _yob=0,_yib=0;
	    int nob=0,nib=0;
	    int jj;for(jj=0;jj<nFolds;jj++){
	      j=bayes[jj];
	      if(nc[j*nGivenData+n]!=1){//not in Dtrain
		nob++;
		_yob+=yp0[j*nGivenData+n];
	      }
	      else{
		nib++;
		_yib+=yp0[j*nGivenData+n];
	      }
	    }
	    if(nob>0) _yp[n]=_yob/nob;
	    else _yp[n]=_yib/nib;//?
	    //	    _yp[n]=_yib/nib;//?
	  }
	  //	  double Dvp=0,Dvm=0;
	  double sigma2=0,skew=0,kurtosis=0,Merr=0,cov=0,skew2=0,skew2p=2./2.,skewa=0;//skew2 signed sigma2 
	  //	  double sigma2=0,skew=0,kurtosis=0,Merr=0,cov=0,skew2=0,skew2p=2.5/2.;//skew2 signed sigma2 
	  //	  double sigma2=0,skew=0,kurtosis=0,Merr=0,cov=0,skew2=0,skew2p=1.25;//skew2 signed sigma2 
	  //	  double Lvar;//Diversity=sigma
	  int nob=0,nib=0;
	  for(n=0;n<nGivenData;n++){
	    double skewj=0,sigma2j=0;int nobj=0;
	    int jj;for(jj=0;jj<nFolds;jj++){
	      j=bayes[jj];
	      if(nc[j*nGivenData+n]!=1){//not in Dtrain
		nob++;nobj++;
		double err=yp0[j*nGivenData+n]-_yp[n];
		double err2=err*err;//double err3=err2*err;
		Merr+=err;
		sigma2j+=err2;//sigma2+=err2;//diversity?
		skewj+=err2*err;//		skew+=err2*err;
		skew2+=(err>0)?pow(err2,skew2p):(-pow(err2,skew2p));//sigma^{1.5}
		//skew2+=(err>0)?pow(err2,0.75):(-pow(err2,0.75));//sigma^{1.5}
		//skew2+=(err>0)?(err2):(-err2);
		kurtosis+=err2*err2;
		int ll;for(ll=0;ll<nFolds;ll++){
		  int l=bayes[ll];
		  if(l!=j){
		    cov+=(err*(yp0[l*nGivenData+n]-_yp[n]));
		  }
		}
	      }
	      else nib++;
	    }//int jj;for(jj=0;jj<nFolds;jj++){
	    skew+=skewj;
	    sigma2+=sigma2j;
	    skewa+=fabs(skewj/nobj)/pow(sigma2j/nobj,1.5);
	  }//for(n=0;n<nGivenData;n++){
	  free(_yp);
	  if(nob>0){
	    //Dvp/=nob;Dvm/=nob;
	    sigma2/=nob; //	    sigma=sqrt(sigma/nob);
	    cov/=(nob*(nFolds-1));
	    double sigma1=sqrt(sigma2);
	    double sigma3=sigma1*sigma1*sigma1;
	    skew=(skew/nob)/sigma3;
	    skew2=skew2/nob/pow(sigma2,skew2p);
	    skewa/=nGivenData;
	    //	    skew2=skew2/nob/sigma2;
	    kurtosis=(kurtosis/nob)/sigma3/sigma1-3;
	    Merr/=nob;
	  }
	  //	  fprintf(stderr,"skw%+.4f;%+.4fkrt%+.0eVar%.1eCor%+g%sk%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov/sigma2,argv[iBayes],K,nCells-1,argv[2]);
	  //	  fprintf(stderr,"D%g=%g-%g.\n",Dvp-Dvm,Dvp,Dvm);
	  //	  fprintf(stderr,"D%gskw%g\n",Dvp-Dvm,skew);
	  //	  fprintf(stderr,"skw%gkrt%gMerr=%g#%s\n",skew,kurtosis,Merr,argv[iBayes]);
	  //	  fprintf(stderr,"skw%+.1ess%+.1ekrt%+.0evar%.1ecov%+.1e%sk%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov,argv[iBayes],K,nCells-1,argv[2]);
	    //	  fprintf(stderr,"skw%g_%gkrt%gvar%gcov%g%s k%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov,argv[iBayes],K,nCells-1,argv[2]);//Div=diversity
	  //	  fprintf(stderr,"skw%gkrt%gMerr=%g#%s k%dN%dA2=%s\n",skew,kurtosis,Merr,argv[iBayes],K,nCells-1,argv[2]);
	  if(lossall){
	    FILE *fp=fopen("lossall.dat","a+");
	    fprintf(fp,"Skw,Skwa,Krt,Var:%+.2e %g %g %g %sk%dN%d_%s\n",skew,skewa,kurtosis,sigma2,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%+.1eskw2%+.1eskwa%gkrt%+.0eVar%.1e%sk%dN%d_%s\n",skew,skew2,skewa,kurtosis,sigma2,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%+.1ess%+.1ekrt%+.0eVar%.1eCor%+g%sk%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov/sigma2,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%+.1ess%+.1ekrt%+.0eVar%.1eCor%+g%sk%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov/sigma2,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%+.4fss%+.4fkrt%+.1evar%.1ecov%+.1e%sk%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%+.3e_%+.3ekrt%+.3evar%.3ecov%+.3e%s k%dN%d_%s\n",skew,skew2,kurtosis,sigma2,cov,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%gkrt%gvar%gcov%g#%s k%dN%dA2=%s\n",skew,kurtosis,sigma2,cov,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%gkrt%gMerr=%g#%s k%dN%dA2=%s\n",skew,kurtosis,Merr,argv[iBayes],K,nCells-1,argv[2]);
	    //	    fprintf(fp,"skw%gkrt%gMerr=%g#%s\n",skew,kurtosis,Merr,argv[iBayes]);
	    fclose(fp);
	  }
	  FILE *fp=fopen("tmp/bayesL.dat","w");
	  fprintf(fp,"%.7e %d #Lbayes,nobmin\n",Lbayes,nobmin);
	  fclose(fp);
	}//end of if(Bayes>0)

	for(n=0;n<nGivenData;n++){
	  yptest[n]=yptrain[n]=ntrain[n]=ntest[n]=yptotal[n]=0;
	  for(j=0;j<nFolds;j++){
	    if(nc[j*nGivenData+n]!=1){//not in Dtrain
	      ntest[n]++;
	      _ntest++;
	      yptest[n]+=yp0[j*nGivenData+n];
	    }
	    else{//in Dtrain
	      ntrain[n]++;
	      _ntrain++;
	      yptrain[n]+=yp0[j*nGivenData+n];
	    }
	    yptotal[n]+=yp0[j*nGivenData+n];
	  }
	  //////
	  if(ntest[n]>0){//not in Dtrain
	    nGivenData0++;
	    yptest[n]/=(ntest[n]);//avoiding zero division?
	    err=yptest[n]-yt[n];
	    err2=err*err;
	    Ltst+= err2;//V2
	    //	    fprintf(stderr,"#n%d,err%g=%g-%g, Ltst=%g,nG=%d\n",n,err,yptest[n],yt[n],Ltst,nGivenData0);
	  }
	  else{
	    fprintf(stderr,"ntest=0 for x[%d] >> Increase the number of folds(bags) or reduce the alpha.\n",n);
	  }
	  yptotal[n]/=nFolds;
	  if(ntrain[n]>0){
	    nGivenData1++;
	    yptrain[n]/=(ntrain[n]);//avoiding zero division?
	    Lib+=square(yptrain[n]-yt[n]);
	    Lemp+=square(yptotal[n]-yt[n]);
	    //	    fprintf(stdout,"yerr[%d]=%g ntrain=%d.Lib=%g=%g/%d.\n",n,square(yptrain[n]-yt[n]),ntrain[n],Lib/nGivenData1,Lib,nGivenData1);//for check
	  }
	  else{
	    fprintf(stderr,"ntrain=0 for x[%d] >> Increase the number of folds(bags) or increase the alpha.\n",n);
	  }
	}
	{
	  _ntest/=nGivenData0;
	  Ltst/=nGivenData0;
	  _ntrain/=nGivenData1;
	  Lib/=nGivenData1;
	  Lemp/=nGivenData1;
	}

#ifdef CHECK
	{//only for check
	  fprintf(stderr,"training vectors:");
	  for(n=0;n<nGivenData;n++){
	    fprintf(stderr,"\n#%4d)%3d+%3d=",n,ntrain[n],ntest[n]);
	    for(j=0;j<nFolds;j++){
	      //	      fprintf(stderr,"%1d ",nc[j*nGivenData+n]);
	      if(nc[j*nGivenData+n]!=1) fprintf(stderr,"%1d ",0);//not in Dtrain
	      else fprintf(stderr,"%1d ",1);
	    }
	  }
	}
#endif
	{//Generalization Loss based on BIC criterion ?
	  int n;
	  if(Lstdm==0){//same as below when nGivenData==nGivenData1==nGivenData2
	    n=nGivenData; 
	    Ltst=(nGivenData1*Lib+nGivenData0*Ltst)/(nGivenData1+nGivenData0);//soso
	  }
	  else if(Lstdm==1){
	    n=nGivenData; 
	    Ltst=(Lib+Ltst)/2.; //good
	  }
	  else if(Lstdm==2){//out-of-bag error
	    n=nGivenData0;//original good
	  }
	  //good
	  //	  Ltst=Ltst+(2.*Lstd/n)*(NN[0]*K*2)*log(n)/2.;//BIC?N:17 K*2= dim(w_i)+dim(M_i)??
	  if(Lstd>0) Ltst=(Ltst+(2.*Lstd/n)*(NN[0]*K)*log(n)/2.)/2.;//BIC?N:17 K= dim(M_i) ??070220
	  
	  //	  Ltst=Ltst+(2.*Ltst/n)*(NN[0]*K)*log(n)/2.;//BIC?N:17 K= dim(M_i) ??070220
	  //L=-sum_{i=1}^n log(p(y_i|x))+ d*log(n)/2        [ p(y_i|x)= (1/sqrt(2 pi s)) exp(-(y_i-hat{y}_i)^2/(2*s))]
          // = sum_{i=1}^n (y_i-hat{y}_i)^2/(2*s) + d*log(n)/2
	  // = Ltst*n/(2*s) + (N*K*2)*log(n)/2              [ d= N*K*2; number of parameters ]
          //L*(2*s/n) = Ltst + (2*s/n)(N*K*2)*log(n)/2
	  //	  double s2=0;  Ltst=Ltst+ (2.*s2/nGivenData0)*(NN[0]*K*2)*log(nGivenData0)/2.;//BIC?N:17 good
	  //	  double s2=0.07787; Ltst=Ltst+ (2.*s2/nGivenData0)*(NN[0]*K*2)*log(nGivenData0)/2.;//BIC?N:17 good
	  //also good but fist step gets the minimum at the biggest N???
	  //	  double s2=0;  Ltst=Lib+ (2.*s2/nGivenData0)*(NN[0]*K*2)*log(nGivenData0)/2.;//BIC?N:17 good
	  //	  double s2=0.03511; Ltst=Ltst+ (2.*s2/nGivenData0)*(NN[0]*K*2)*log(nGivenData0)/2.;//BIC?N:17 good
	}
	{
	  double pn,E2sum=Lemp*nGivenData1;
	  H0=0;
	  for(n=0;n<nGivenData;n++){
	    if(ntrain[n]>0){
	      pn=square(yptotal[n]-yt[n])/E2sum+1e-20;
	      H0-= pn*log(pn);
	    }
	  }
	  H0/=log(nGivenData1);
	}
	
	Lvar0=Lvar=ntrainall=ntestall=L0=Ltr=0;
	UC=ucn=0;
	for(n=0;n<nGivenData;n++){
	  uc[n]=0;
	  for(j=0;j<nFolds;j++){
	    if(nc[j*nGivenData+n]!=1){//not in Dtrain
	      err2=square(yp0[j*nGivenData+n]-yptest[n]);
	      Lvar+=err2;
	      uc[n]+=err2;
	      L0+=square(yp0[j*nGivenData+n]-yt[n]);//original out-of-bag loss
	    }
	    else{
	      Lvar0+=square(yp0[j*nGivenData+n]-yptrain[n]);
	      // Ltr+=square(yp0[j*nGivenData+n]-yt[n]);//original training loss
	    }
	  }
	  if(ntest[n]>2){
	    uc[n]/=(ntest[n]*(ntest[n]-1));//uncertainty
	    UC+=uc[n];//uncertainty
	    ucn++;
	  }
	  ntestall+=ntest[n];
	  ntrainall+=ntrain[n];
	}

	UC/=ucn;//average uncertainty
	
	Lvar/=(ntestall);
	Lvar0/=(ntrainall);
	L0/=(ntestall);
	//	Ltr/=(ntrainall);

	nfiles=nFolds;
	num=nGivenData;
	LD=UC;
	//	LD=Ltst;
	int nn;
	for(nn=0;nn<nTrainData;nn++){//NewLtr
	  for(j=0;j<nFolds;j++){
	    n=nr[j*nTrainData+nn];
	    Ltr+=square(yp0[j*nGivenData+n]-yt[n]);//original training loss
	  }
	}
	Ltr/=(nFolds*nTrainData);
	if(Bayes>0){
	  Ltst=Lbayes;
	}
	{
	  fp=fopen("predict+oob.dat","w");
	  for(n=0;n<nGivenData;n++){
            fprintf(fp,"%.7e %d %.7e #Y^,t,Y\n",yptest[n],n,yt[n]);
	  }
          fclose(fp);
	}
	fp=fopen("./loss+.dat","w");
	//	fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib LtrE H0 L0\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0,L0);
	fprintf(fp,"%.7e %.7e %.7e %.7e %.7e %d %d %.7e %.7e %.7e %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib LtrE H0 L0\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0,L0,Ltr);
	fclose(fp);
	//	fprintf(stderr,"#method:NoBoost Ltst(n%d)%.3e Lib(n%d)%.3e %.3e #Ltst,Lib,H0\n",nGivenData0,Ltst,nGivenData1,Lib,H0);
	fprintf(stderr,"#Ltst(n%.1f*%d)%.3e Lib(n%.1f*%d)%.3e H%.3e, nTest%.1f=%.1f=exp(-%.2f)*nGiven%dmethod:%sN%s.\n",_ntest,nGivenData0,Ltst,_ntrain,nGivenData1,Lib,H0,meannTestData,exp(-alpha)*nGivenData,alpha,nGivenData,argv[2],argv[3]);
	//	fprintf(stderr,"#method:NoBoost Ltst(n%.1f*%d)%.7e Lib(n%.1f*%d)%.7e H%.7e, nTest%.1f=%.1f=(1-exp(a%.2f))*nGiven%d.\n",_ntest,nGivenData0,Ltst,_ntrain,nGivenData1,Lib,H0,meannTestData,(1.-exp(-alpha))*nGivenData,alpha,nGivenData);
	//	fprintf(stderr,"#BAGGING2#%13.7e %.3e %1.3e %.3e %d %d %.3e %.3e %.4e #Ltst Lvar Lvarval Lvar0 nfiles num Ltr LtrE ntrain/n=%.3f.\n",Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0,((double)ntrainall/(ntrainall+ntestall)));
      }//closing of if(BAGGING==BAGwithoutVal)
      else if(BAGGING==NoBAG){//?? NoBoost+NoBag??
	msesum=Lhatmean=0;
	for(n=0;n<nFolds;n++){
	  msesum+=mse[n];
	  Lhatmean+=Lhat[n];
	}
	msemean=msesum/nFolds;
	Lhatmean /=nFolds;
	for(n=0;n<nFolds;n++){
	  if((err=mse[n]-msemean)>0){
	    msestdp+=err*err;
	    n_msestdp++;
	  }
	  else {
	    msestdm+=err*err;
	    n_msestdm++;
	  }
	}
	if(n_msestdp>0) msestdp=sqrt(msestdp/n_msestdp); else msestdp=0;
	if(n_msestdm>0) msestdm=sqrt(msestdm/n_msestdm); else msestdm=0;
	int nn;
	n_msesum=msetrainsum=n_msetrainsum=n_cells2sum=0;//??20130123
	for(nn=0;nn<nFolds;nn++){
	  n_msesum+=n_mse[nn];
	  msetrainsum+=msetrain[nn];
	  n_msetrainsum+=n_msetrain[nn];
	  n_cells2sum+=n_cells2[nn];//
	  net_printf(strx7("#%3d %13.7e %13.7e %3d %3d %13.7e N%d#n Ltest Lib bestIT Nr Lhat? N\n",
			   nn,mse[nn],msetrain[nn],is[nn],n_cells2[nn],Lhat[nn],NN[0]));//modified 051016
	}
	msetrainmean=msetrainsum/nFolds;
	net_printf(strx9("#%10.5e %8.3e %8.3e %10.5e %d %d N%d %3d %10.5e #<Ltest Lstd Lstd2 Lib?(see ./loss+.dat) n nFolds N Nr Lhat>\n",
                        msemean,
                        msestdp,msestdm,
                        msetrainmean,
                        nGivenData,
                        nFolds,_nCells,(int)n_cells2sum/nFolds,
                        Lhatmean));
	fp=fopen("./tmp/rsresult.dat","w+");
	fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %d %d %3d %3d %d %13.7e #Ltest Lstd Lstd2 Lib n nFolds N Nr a*n Lhat\n",
		msemean,//msesum/n_msesum,
		msestdp,msestdm,
		msetrainmean,//msetrainsum/n_msetrainsum,
		nGivenData,
		nFolds,_nCells,
		(int)n_cells2sum/nFolds,
		(int)(nGivenData-meannTestData+0.5),
		Lhatmean);
	fclose(fp);
	//      fprintf(stderr,"method=%s, <nTrainData>=%d, <nTestData>=%d. E:%d.\n",argv[2],(int)(n_msetrainsum/nFolds),(int)(n_msesum/nFolds),ERR);
	net_printf(strx3("#method=%s, <nTrainData>=%d, <nTestData>=%d.\n",argv[2],(int)(n_msetrainsum/nFolds),(int)(n_msesum/nFolds)));
      }
    }//closing  if(Boost==NoBoost){//gradient-based boosting
    else if(Boost==GbBoost){//gradient-based boosting
      net_printf(strx1("Doing '%s'\n",cmd2));
      system(cmd2);//exec meanpred for ensemble
      sprintf(cmd,"cat %s | awk '{print $%d}' > %s",fn_bagging,K+1,fn_target);system(cmd);
      sprintf(&cmd[strlen(cmd)]," ;./pred2y_ts predict+.dat %s > /dev/null",fn_target);system(cmd);
      //fprintf(stderr,"Prediction saved in '%s', Loss in'./loss+.dat'. Bagging with method=%s. \n",fn_target,argv[2]);
      
      if(BAGGING==BAGwithVal){
	int t;
	double csum=0;
	double LD=0,Ltst,Lvar=0,Lvarval=0,Lvar0=0,Lib,Lemp,H0;
	double *yp=(double*)malloc(sizeof(double)*nValData);
	{
	  fp=fopen("predict+.dat","r");
	  int _i;
	  for(i=0;i<nValData;i++) {
	    fgets(line1,LINESIZE,fp);
	    sscanf(line1,"%lf%d%lf",&ypt[t_boost*nValData+i],&_i,&yt[i]);
	  }
	  fclose(fp);
	}
	fprintf(stderr,"ct=");for(t=1;t<=t_boost;t++) fprintf(stderr,"(%d)%f ",t,cbst[t]); fprintf(stderr,"\n");
#define NEWct
#ifdef NEWct
	{
	  double *hbst=(double*)malloc(sizeof(double)*(t_boost+1));
	  if(rangedatamode==0){
	    for(t=1;t<=t_boost;t++){
	      hbst[t]=exp(-(t-1.)/tau_h);
	      //	      hbst[t]=exp(-2.*(t-1.)/t_boost);
	      csum+=hbst[t];
	    }
	  }
	  else if(rangedatamode==1){
	    for(t=1;t<=t_boost;t++){
	      hbst[t]=exp((t-1.)/tau_h);
	      csum+=hbst[t];
	    }
	  }
	  else if(rangedatamode==2){
	    for(t=1;t<t_boost;t++){
	      hbst[t]=0;
	    }
	    csum=hbst[t_boost]=1.;
	  }
	  //	  for(t=1;t<=t_boost;t++){
	  //	    hbst[t]=exp(-(t-1.)/tau_h);//
	  //	    //	    hbst[t]=exp(-2.*(t-1.)/t_boost);
	  //	    csum+=hbst[t];
	  //	  }
	  
	  fprintf(stderr,"ht=");  for(t=1;t<=t_boost;t++) fprintf(stderr,"(%d)%f ",t,hbst[t]); fprintf(stderr,"\n");
	  for(t=1;t<=t_boost;t++) hbst[t]/=csum;
	  for(i=0;i<nValData;i++){
	    yp[i]=0;for(t=1;t<=t_boost;t++) yp[i]+=hbst[t]*ypt[t*nValData+i];
	  }
	  free(hbst);
	}
#else
	for(t=1;t<=t_boost;t++) csum+=cbst[t];
	for(i=0;i<nValData;i++){
	  yp[i]=0;
#define GBBOOSTTEST 0
	  //#define GBBOOSTTEST 1
	  //#define GBBOOSTTEST 2
#if GBBOOSTTEST==0
	  for(t=1;t<=t_boost;t++) yp[i]+=cbst[t]*ypt[t*nValData+i];
	  yp[i]/=csum;
#elif GBBOOSTTEST==1 //equivalent weight  boosting
	  for(t=1;t<=t_boost;t++) yp[i]+=ypt[t*nValData+i];
	  yp[i]/=t_boost;
#elif GBBOOSTTEST==2  //gradient boosting
	  yp[i]=ypt[t_boost*nValData+i];
#endif
	}
#endif
	double E2sum,e2;
	int ntrain;
	Ltst=E2sum=0;
	H0=0;
	for(i=0;i<nValData;i++){
	  e2 = square(yp[i]-yt[i]);
	  Ltst+=e2;
	  if(nc[i]==1){
	    E2sum+=e2;
	    ntrain++;
	  }
	}
	Ltst/=nValData;
	Lib=E2sum/ntrain;
	Lemp=Lib;
	double pi;
	H0=0;
	for(i=0;i<nValData;i++){
	  if(nc[i]==1){
	    pi=square(yp[i]-yt[i])/E2sum+1e-20;
	    H0-= pi*log(pi);
	  }
	}
	if(ntrain==1) H0=1; else H0/=log(ntrain);
	int nfiles=t_boost,num=nGivenData;
	sprintf(fn_boost,"%s/tmp/BstGb:%s+%s+m%da%gs%dN:%d:%db%dpred.dat",rdir0,fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,seed,NN[0],NStep,t_boost);//
	fp=fopen(fn_boost,"w");
	//fn_boost(for val) involves 
	//yp[t*nGivenData+i]
	for(i=0;i<nValData;i++) fprintf(fp,"%.7e %d %.7e\n",ypt[t_boost*nValData+i],i,yt[i]);
	fclose(fp);
	fprintf(stderr,">>'%s' for boosting is created\n",fn_boost);
	{
#ifdef NEWp0
	  //pred0 for BAGwithVal is already created through BAGwithoutVal?
#else
	  char fn_boost0[fnsize];
	  //	  sprintf(fn_boost0,"%s/tmp/BstGb:%s+m%da%.2fs%dN:%d:%db%dpred0.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t_boost);//
	  sprintf(fn_boost0,      "BstGb:%s+m%da%gs%dN:%d:%db%dpred0.dat",fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t_boost);//
	  sprintf(fn_boost,"%s/tmp/BstGb:%s+%s+m%da%gs%dN:%d:%db%dpred0.dat",rdir0,fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,seed,NN[0],NStep,t_boost);//
	  sprintf(cmd,"ln -s %s %s",fn_boost0,fn_boost);system(cmd);
#endif
	  //	  sprintf(cmd,"if [ ! -f %s ];then ln -s %s %s ; fi",fn_boost,fn_boost0,fn_boost);system(cmd);

	  sprintf(fn_boost,"%s/tmp/BstGb:%s+%s+m%da%gs%dN:%d:%db%dis.dat",rdir0,fnbody(fn_givendata,pr),fnbody(fn_bagging,pr1),method,alpha,seed,NN[0],NStep,t_boost);//
	  fp=fopen(fn_boost,"w");
	  fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.0,nGivenData,0.0,0,0,nValData,0.0,0.0);
	  fclose(fp);
	}

	fp=fopen("./loss+.dat","w");
	fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0);
	fprintf(stderr,"#Method:Bagwithval boost%d %.7e %.7e %.7e #Ltst,Lib,H0\n",t_boost,Ltst,Lib,H0);
	fclose(fp);
	
	fp=fopen("./predict+.dat","w");
	for(i=0;i<nValData;i++){
	  if(intmode) fprintf(fp,"%d %d %.7f #Y^ t Y\n",(int)(yp[i]+0.5),i,yt[i]);
	  else fprintf(fp,"%.7f %d %.7f # Y^ t Y\n",yp[i],i,yt[i]);
	}
	fclose(fp);
	sprintf(cmd,"cat %s | awk '{print $%d}' > %s",fn_bagging,K+1,fn_target);system(cmd);
	sprintf(&cmd[strlen(cmd)]," ;./pred2y_ts predict+.dat %s > /dev/null",fn_target);system(cmd);
	//fprintf(stderr,"Prediction saved in '%s', Loss in'./loss+.dat'. Bagging with method=%s. \n",fn_target,argv[2]);
      }
      else if(BAGGING==BAGwithoutVal){
	double LD=0,Ltst,Lvar=0,Lvarval=0,Lvar0=0,Lib,Lemp,H0;
	double *yp=(double*)malloc(sizeof(double)*nGivenData);
	double bias2meanGlobal;
	//	cbst[-1-t_boost]=0.5;??
#define BiasV1 1
#define BiasV2 2
#define BiasMode BiasV2 
#undef BiasEmp
#define BiasEmp
#if BiasMode == BiasV1
	{
	  fp=fopen("predict+.dat","r");
	  int _i;
	  for(i=0;i<nGivenData;i++) {
	    fgets(line1,LINESIZE,fp);
	    sscanf(line1,"%lf%d%lf",&ypt[t_boost*nGivenData+i],&_i,&yt[i]);
	    errt[i]=ypt[t_boost*nGivenData+i]-yt[i];
	  }
	  fclose(fp);
	}
	{//new hypothesis
	  bias2meanGlobal=0;
	  for(i=0;i<nGivenData;i++){
	    double bias2mean=0;int t;
	    for(t=1;t<=t_boost;t++) bias2mean += (ypt[t*nGivenData+i]-yt[i]);
	    bias2mean/=t_boost;
	    //for(t=t_boost;t<=t_boost;t++) bias2mean += (ypt[t*nGivenData+i]-yt[i]);//simple error
	    bias2[i]=square(bias2mean);
	    bias2meanGlobal+=bias2[i];
	  }
	  bias2meanGlobal/=nGivenData;
	}
#elif BiasMode == BiasV2
	{
	  fp=fopen("predict+.dat","r");
	  int _i;
	  double varmax=0;
	  for(i=0;i<nGivenData;i++) {
	    fgets(line1,LINESIZE,fp);
	    double _d,_bias,_var;
	    sscanf(line1,"%lf%d%lf%lf%lf%d%lf%lf%lf",&ypt[t_boost*nGivenData+i],&_i,&yt[i],&_d,&_var,&_i,&_d,&_d,&_bias);
	    bias2[t_boost*nGivenData+i]=_bias;
	    errt[i]=ypt[t_boost*nGivenData+i]-yt[i];
	    if(varmax<_var) varmax=_var;
	    var[i]=_var;
	  }
	  fclose(fp);
	  //new hypothesis
	  bias2meanGlobal=0;
	  for(i=0;i<nGivenData;i++){
	    double bias2mean=0;int t;
//	    for(t=1;t<=t_boost;t++){
//	      if(fabs(bias2[t*nGivenData+i]-(ypt[t*nGivenData+i]-yt[i]))>1e-7){
//		fprintf(stderr,"t%d,i%d)0!=%.7e=%f-(%f)(=%f-%f)\n",t,i,
//			bias2[t*nGivenData+i]-(ypt[t*nGivenData+i]-yt[i]),
//			bias2[t*nGivenData+i],
//			ypt[t*nGivenData+i]-yt[i],
//			ypt[t*nGivenData+i],yt[i]);
//	      }
//	    }
	    bias2mean=0;
#ifdef BiasEmp
	    for(t=1;t<=t_boost;t++) bias2mean += (ypt[t*nGivenData+i]-yt[i]);
	    bias2mean/=t_boost;
	    bias2[i]=square(bias2mean)+varmax-var[i];//good 
	    //bias2[i]=square(bias2mean);//soso
	    //	    fprintf(stderr,"[%d]bias2=%f=bias2(%f)+varmax(%f)-vari(%f),%f(bias)=%f(yp)-%f(y)\n",i,bias2[i],square(bias2mean),varmax,var[i],bias2mean,ypt[t_boost*nGivenData+i],yt[i]);
#else //BiasTheotical
	    for(t=1;t<=t_boost;t++) bias2mean += bias2[t*nGivenData+i];
	    bias2mean/=t_boost;
	    bias2[i]=square(bias2mean);
	    //	    fprintf(stderr,"[%d]bias2=%f=bias2(%f),varmax(%f),vari(%f),%f(bias)=%f(yp)-%f(y)\n",i,bias2[i],square(bias2mean),varmax,var[i],bias2mean,ypt[t_boost*nGivenData+i],yt[i]);
#endif
	    //	    if(bias2[i]<0) bias2[i]=0;
	    bias2meanGlobal+=bias2[i];
	  }
	  bias2meanGlobal/=nGivenData;
	}

//	if(1==0){
//	  fp=fopen("predict+.dat","r");
//	  int _i;
//	  for(i=0;i<nGivenData;i++) {
//	    fgets(line1,buffsize,fp);
//	    double _d,_bias;
//	    sscanf(line1,"%lf%d%lf%lf%lf%d%lf%lf%lf",&ypt[t_boost*nGivenData+i],&_i,&yt[i],&_d,&_d,&_i,&_d,&_d,&_bias);
//	    bias2[t_boost*nGivenData+i]=_bias;
//	    //	    sscanf(line1,"%lf%d%lf%lf%lf%d%lf%lf%lf",&ypt[t_boost*nValData+i],&_i,&yt[n],&_d,&_bias,&_i,&_d,&_d,&_d);
//	    //	    bias2[i]=_bias;
//	  }
//	  fclose(fp);
//	  //new hypothesis
//	  bias2meanGlobal=0;
//	  for(i=0;i<nGivenData;i++){
//	    double bias2mean=0;int t;
//	    for(t=1;t<=t_boost;t++) bias2mean += bias2[t*nGivenData+i];
//	    //for(t=1;t<=t_boost;t++) bias2mean += (ypt[t*nGivenData+i]-yt[i]);//original
//	    bias2mean/=t_boost;
//	    bias2[i]=square(bias2mean);
//	    //	    bias2[i]=square(bias2mean)+10.0;
//	    bias2meanGlobal+=bias2[i];
//	  }
//	  bias2meanGlobal/=nGivenData;
//	}
#endif

#undef useETA
#ifdef useETA
	double *eta=(double*)malloc(sizeof(double)*nGivenData);
	{//new hypothesis2
	  if(t_boost==1) for(i=0;i<nGivenData;i++) eta[i]=1.;
	  else{
	    double etasum=0;
	    double dpthresh=0.000001/nGivenData;
	    double *bias2_1=(double*)malloc(sizeof(double)*nGivenData);
	    for(i=0;i<nGivenData;i++){
	      double bias2mean=0;int t;
	      for(t=1;t<t_boost;t++) bias2mean += (ypt[t*nGivenData+i]-yt[i]);
	      bias2mean/=(t_boost-1);
	      bias2_1[i]=bias2mean*bias2mean;
	      double dB=square(bias2[i]-bias2_1[i]);
	      double dp=pbst[i]-pbst_1[i];
//	      if(i>=145){
//		fprintf(stderr,"check\n");
//	      }
	      if(fabs(dp)>dpthresh) eta[i]= -dB/dp;
	      else eta[i]=0;
	      if(eta[i]>0) eta[i]=1;
	      else if(eta[i]<0) eta[i]=-1;
	      //	      if(eta[i]<0) eta[i]=0;
	      etasum+=eta[i];
	    }
	    //	    etasum/=nGivenData;
	    //	    for(i=0;i<nGivenData;i++) eta[i]/=etasum;
	    free(bias2_1);
	  }
	}
#endif
	{
#define GoldenSection 1
#define LineSearchGSL 2
#define EntropySaidaika 3
	  //#define Linesearch GoldenSection
	  //#define LineSearch LineSearchGSL 
#define LineSearch EntropySaidaika
#if GBBOOSTTEST==1
	  cbst[t_boost]=1./nGivenData;
#else
	  ///////////////
	  double mse2=0;
	  for(i=0;i<nGivenData;i++) mse2+=bias2[i];
	  mse2/=nGivenData;
	  {
	    double bias2max=0;
	    for(i=0;i<nGivenData;i++) if(bias2max<bias2[i]) bias2max=bias2[i];
	    double bias2min=1e20;
	    for(i=0;i<nGivenData;i++) if(bias2min>bias2[i]) bias2min=bias2[i];
	    fprintf(stderr,"bias2 min=%e max=%e\n",bias2min,bias2max);
	  }

#if LineSearch == GoldenSection
	  double r=2/(3.+sqrt(5.));
	  double a,b,c,d,c2,d2;
	  double fc,fd;
	  double tolerance=0.00001;
	  int t;

#if GBBOOSTTEST==2 
	  b=1.;tolerance=0.00001; a=0.00001;
#else
	  //	  b=1.;tolerance=0.00001; a=0.00001;
	  b=2.*(1./(2.*mse2)); tolerance=0.001*b; a=0.001*b;
#endif
	  c=a+r*(b-a); d=b-r*(b-a);
	  //#define ReducePredictionError 
#ifdef ReducePredictionError 
	  fc=0;c2=1./sqrt(c);for(i=0;i<nGivenData;i++) if(nc[i]!=1) fc+=wbst[i]*c2*exp(c*bias2[i]);
	  fd=0;d2=1./sqrt(d);for(i=0;i<nGivenData;i++) if(nc[i]!=1) fd+=wbst[i]*d2*exp(d*bias2[i]);
	  for(t=0;t<1000;t++){
	    if(fc>fd){
	      a=c;c=d;fc=fd;d=b-r*(b-a);
	      if(d-c<= tolerance) {cbst[t_boost]=c;break;}
	      fd=0;d2=1./sqrt(d);for(i=0;i<nGivenData;i++) if(nc[i]!=1) fd+=wbst[i]*d2*exp(d*bias2[i]);
	    } else {
	      b=d;d=c;fd=fc;c=a+r*(b-a);
	      if(d-c<=tolerance) {cbst[t_boost]=d;break;}
	      fc=0;c2=1./sqrt(c);for(i=0;i<nGivenData;i++) if(nc[i]!=1) fc+=wbst[i]*c2*exp(c*bias2[i]);
	    }
	  }
#else
	  fc=0;c2=1./sqrt(c);for(i=0;i<nGivenData;i++) fc+=wbst[i]*c2*exp(c*bias2[i]);
	  fd=0;d2=1./sqrt(d);for(i=0;i<nGivenData;i++) fd+=wbst[i]*d2*exp(d*bias2[i]);
	  for(t=0;t<1000;t++){
	    if(fc>fd){
	      a=c;c=d;fc=fd;d=b-r*(b-a);
	      if(d-c<= tolerance) {cbst[t_boost]=c;break;}
	      fd=0;d2=1./sqrt(d);for(i=0;i<nGivenData;i++) fd+=wbst[i]*d2*exp(d*bias2[i]);
	    } else {
	      b=d;d=c;fd=fc;c=a+r*(b-a);
	      if(d-c<=tolerance) {cbst[t_boost]=d;break;}
	      fc=0;c2=1./sqrt(c);for(i=0;i<nGivenData;i++) fc+=wbst[i]*c2*exp(c*bias2[i]);
	    }
	  }
	  fprintf(stderr,"[%d]ct=%f,fc[%f]=%f,fd[%f]=%f,mse2=%e,1/(2*mse2)=%e",t_boost,cbst[t_boost],c,fc,d,fd,mse2,1./(2.*mse2));
#endif //closing ReducePrediction

#elif LineSearch == LineSearchGSL
	  double *par[4];
	  double _nGivenData=nGivenData;
	  par[0]=(double *)&_nGivenData;
	  par[1]=&wbst[0];
	  par[2]=&bias2[0];
	  int status;
	  int iter=0,max_iter=10000;
	  double a,b,m,tolerance;
	  //	  b=1.; a=0.00001; m=0.5;tolerance=0.001*b;
	  b=2.*(1./(2.*mse2));a=1e-7*b;m=1./(2.*mse2);tolerance=1e-5*b;
	  {//for error check
	    double fa,fm,fb;
	    for(;;){
	      fa=lossboost(a,&par);
	      fm=lossboost(m,&par);
	      fb=lossboost(b,&par);
	      fprintf(stderr,"??f(a=%.7f)=%.7f>f(m=%.7f)=%.7f<f(b=%.7f)=%.7f\n",a,fa,m,fm,b,fb);
	      if(fm>fb){m=b;b=2.*m;}
	      else break;
	    }
	  }
	  const gsl_min_fminimizer_type *T;
	  gsl_min_fminimizer *s;
	  gsl_function F;
	  F.function = &lossboost;
	  F.params = (void *)&par;
	  T = gsl_min_fminimizer_brent;
	  //T= gsl_min_fminimizer_goldensection;
	  s = gsl_min_fminimizer_alloc (T);
	  gsl_min_fminimizer_set (s, &F, m, a, b);
	  do{
	    iter++;
	    status = gsl_min_fminimizer_iterate (s);
	    
	    m = gsl_min_fminimizer_x_minimum (s);
	    a = gsl_min_fminimizer_x_lower (s);
	    b = gsl_min_fminimizer_x_upper (s);
	    
	    status 
	      = gsl_min_test_interval (a, b, tolerance, 0.0);
	    
	    if (status == GSL_SUCCESS) printf ("Converged:\n");
	    
	    printf ("[%5d] a%.7f b%.7f m%.7f b-a=%+.7f\n",
		    iter, a, b, m, b - a);
	  } while (status == GSL_CONTINUE && iter < max_iter);
	  //	  cbst[t_boost]=1.0*(2./log((double)nGivenData))/mse2/(1.+(t_boost-1)/4.);//for check! 
	  cbst[t_boost]=m;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check! good
	  //cbst[t_boost]=0.2/mse2/(1.+(t_boost-1)/4.);//for check! not so good
	  //	  cbst[t_boost]=0.4/mse2/(1.+(t_boost-1)/4.);//for check! good for sin2e2noise0.5
	  //	  cbst[t_boost]=0.3/mse2/(1.+(t_boost-1)/4.);//for check! good
	  //	  cbst[t_boost]=30.*(2./log(nGivenData))/mse2/(1.+(t_boost-1)/4.);//for check! good
	  //	  cbst[t_boost]=(2./nGivenData)/mse2;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check!
	  //	  cbst[t_boost]=m=m/10.;//for check
	  //	  cbst[t_boost]=m=m/(1.+(t_boost-1)/4.);//for check
	  //cbst[t_boost]=m=1./3./mse2/(1.+(t_boost-1)/4.);//for check
	  //	  cbst[t_boost]=m=1./2./mse2/(1.+(t_boost-1)/4.);//for check
	  //	  cbst[t_boost]=m=1./2./mse2/(1.+(t_boost-1)/4.);//for check
	  //	  cbst[t_boost]=1./bias2max;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check ct=(1)2.926910
	  //	  cbst[t_boost]=1./mse2;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check 
	  //	  cbst[t_boost]=0.5/mse2;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check good
	  //	  cbst[t_boost]=m;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/8.);//for check good
	  //	  cbst[t_boost]=m;//original
	  //	  cbst[t_boost]=m;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);if(t_boost>5) cbst[t_boost]=m=cbst[5];//for check good
	  //	  cbst[t_boost]=2.*m;if(t_boost>1) cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check good
	  //	  cbst[t_boost]=m=cbst[1]/(1.+(t_boost+3.)/4.);//???
	  //	  cbst[t_boost]=m=cbst[1]/(1.+(t_boost-1)/4.);//for check
	  //	  cbst[t_boost]=m=m/(1.+(t_boost-1)/4.);//for check soso?
	  //	  m/=(t_boost/2.);//for check soso?
	  //	  m/=(1.+t_boost/2.);//for check soso?
	  //	  m/=t_boost;//for check -> seems good but slow
	  //	  cbst[t_boost]=1.;
	  //	  if(t_boost>=2) cbst[t_boost]=cbst[1];
	  //	  cbst[t_boost]=m=t_boost/(2.*mse2);//??for check bad
	  //	  cbst[t_boost]=m=t_boost/(2.*mse2);//??for check bad
	  //	  cbst[t_boost]=m=1./(2.*mse2)/t_boost;//??for check bad
	  //	  cbst[t_boost]=1./(2.*mse2);//??for check bad
	  //	  cbst[t_boost]=m=1.;//for check
	  //	  fprintf(stderr,"[%d]ct=%f,f(%f)=%f,iter=%d,1/2mse=%f\n",t_boost,cbst[t_boost],m,lossboost(m,&par),iter,1./2./mse2);
	  fprintf(stderr,"[%d]ct=%f,f(%f)=%f,iter=%d,.5/mse2=%f\n",t_boost,cbst[t_boost],m,lossboost(m,&par),iter,0.5/mse2);
	  //	  fprintf(stderr,"log(%d)=%e\n",nGivenData,log((double)nGivenData));
	  gsl_min_fminimizer_free(s);
#elif LineSearch == EntropySaidaika
#if BiasMode == BiasV1
	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2;if(t_boost>1) cbst[t_boost]=cbst[1]/(1.+(t_boost-1)/4.);//for check! 
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/8.);//for check! 
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//good for sin2e2
	  //	cbst[t_boost]=2.0/120./log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! 
	//	cbst[t_boost]=1.0/(log((double)nGivenData)/2.)/mse2/(1.+(t_boost-1)/4.);//for check! 
	//	cbst[t_boost]=1.0/(log((double)nGivenData)/2.)/mse2/(1.+(t_boost-1)/4.);//for check! 
	//	cbst[t_boost]=1.0*(2./log((double)nGivenData))/mse2/(1.+(t_boost-1)/4.);//for check! 
#elif BiasMode == BiasV2
	  cbst[t_boost]=eta1/log((double)nGivenData)/mse2/(1.+(t_boost-1)/tau_c);//for check! 
	  //	  if(t_boost==1) cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/tau_c);//Ltst0=8.917e-02 t=20
	  //	  else cbst[t_boost]=cbst[1]/(1.+(t_boost-1)/tau_c);
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/8.);//for check! 
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2;//for NEWct check >> NG
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! 
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/2.);//for check for friedman! not so
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2;if(t_boost>1) cbst[t_boost]=cbst[1]/(1.+(t_boost-1)/4.);//for check! 
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/8.);//for check for friedman! not so
	  //	  cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/3.);//for check! Ltst0=9.472e-02
	  //cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/2.);//for check! Ltst0=9.491e-02
	  //cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! Ltst0=1.111e-01
	  //cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! Ltst0=9.494e-02??
	  //cbst[t_boost]=1.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! Ltst0=9.561e-02
	  //cbst[t_boost]=2.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/1.);//for check! Ltst0=9.620e-02
	  //	  cbst[t_boost]=3.0/log((double)nGivenData)/mse2/(1.+(t_boost-1)/4.);//for check! Ltst0=9.836e-02
#endif
#endif //closeing #if LineSearch == LineSearchGSL
#endif
	}
	{
	  //	  double wsum=0,emean=0;
	  double wsum=0;
	  //	  double cbst2=1./sqrt(cbst[t_boost]);
	  double cT=cbst[t_boost];
	  //	  double cT2=1./sqrt(cT);
	  for(i=0;i<nGivenData;i++){
#ifdef useETA
	    wbst1[i]=wbst[i]*cT2*exp(cT*bias2[i]*eta[i]);//for new hypothesis1&2
#else
	    wbst1[i]=wbst[i]*exp(cT*bias2[i]);//for wbst1 with wbst(original)
	    //	    wbst1[i]=wbst[i]*cT2*exp(cT*bias2[i]);//for wbst1 with wbst(original)
#endif
	    //	    wbst1[i]*=(1.)/(1.+exp(3.*(square(errt[i])/var[i]-9.)));//good noisecanceller;Ltst0=8.958e-02sin2e2outlier Ltst0=9.012e-02sin2e2
	    //	    wbst1[i]*=(1.)/(1.+exp(3.*(bias2[i]/bias2meanGlobal-9.)));//good noisecanceller;Ltst0=8.958e-02sin2e2outlier Ltst0=9.012e-02sin2e2
	    //	    wbst1[i]*=(1.)/(1.+exp(2.*(bias2[i]/bias2meanGlobal-9.)));//SoSo Ltst0=9.016e-02
	    //wbst1[i]*=(1.)/(1.+exp(4.*(bias2[i]/bias2meanGlobal-9.)));//SoSo Ltst0=8.958e-02sin2e2outlier Ltst0=9.080e-02sin2e2
	    //wbst1[i]*=(1.)/(1.+exp(3.*(bias2[i]/bias2meanGlobal-16.)));//NG for Ltst0=1.941e-01sin2e2outlier Good for Ltst0=9.080e-02sin2e2
	    //wbst1[i]*=(1.)/(1.+exp(3.*(bias2[i]/bias2meanGlobal-7.)));//good Ltst0=8.958e-02sin2e2outlier Ltst0=9.250e-02sin2e2
	    //	    wbst1[i]=cT2*exp(cT*bias2[i]);//for wbst1 without wbst(check)
	    //	    wbst1[i]*=(1.+0.05*bias2[i]/bias2meanGlobal)/(10.+exp(4.*(bias2[i]/bias2meanGlobal-9.)));//good
	    //	    for(i=0;i<nGivenData;i++) {
	    //	      if(bias2[i]>bias2meanGlobal*16.){
	    //		fprintf(stderr,"bias2[%d]=%e>%e*16:outlier?\n",i,bias2[i],bias2meanGlobal);
	    //		pbst1[i]=0;
	    //	      }
	    //	    }
	    wsum+=wbst1[i];
	    //	    emean+=bias2[i];
	  }
	  for(i=0;i<nGivenData;i++){
	    pbst1[i]=wbst1[i]/wsum;
//	    if(bias2[i]>bias2meanGlobal*16.||i>=199){
//	      fprintf(stderr,"normalized bias2[%d]=%f. p=%f. outlier?\n",i,bias2[i]/bias2meanGlobal,pbst1[i]);
//	    }
	  }
	}
	///
	//	emean/=nGivenData;
	//	double pmean=0;
	//	for(i=0;i<nGivenData;i++){
	//	  pbst1[i] *=(1+0.05*bias2[i]/emean)/(1.+exp(bias2[i]-emean*9.));
	//	  pmean+=pbst1[i];
	//	}
	//	for(i=0;i<nGivenData;i++) pbst1[i]/=pmean;
	{
	  int t;
	  fprintf(stderr,"ct =");  for(t=1;t<=t_boost;t++) fprintf(stderr,"(%d)%f ",t,cbst[t]); fprintf(stderr,"\n");
	  double csum=0;
#ifdef NEWct
	  {
	    double *hbst=(double*)malloc(sizeof(double)*(t_boost+1));
	    if(rangedatamode==0){
	      for(t=1;t<=t_boost;t++){
		hbst[t]=exp(-(t-1.)/tau_h);
		//	      hbst[t]=exp(-2.*(t-1.)/t_boost);
		csum+=hbst[t];
	      }
	    }
	    else if(rangedatamode==1){
	      for(t=1;t<=t_boost;t++){
		hbst[t]=exp((t-1.)/tau_h);
		csum+=hbst[t];
	      }
	    }
	    else if(rangedatamode==2){
	      for(t=1;t<t_boost;t++){
		hbst[t]=0;
	      }
	      csum=hbst[t_boost]=1.;
	    }
	    fprintf(stderr,"ht=");  for(t=1;t<=t_boost;t++) fprintf(stderr,"(%d)%f ",t,hbst[t]); fprintf(stderr,"\n");
	    for(t=1;t<=t_boost;t++) hbst[t]/=csum;
	    for(i=0;i<nGivenData;i++){
	      yp[i]=0;for(t=1;t<=t_boost;t++) yp[i]+=hbst[t]*ypt[t*nGivenData+i];
	    }
	    free(hbst);
	  }
#else
	  for(t=1;t<=t_boost;t++) csum+=cbst[t];
	  for(i=0;i<nGivenData;i++){
	    yp[i]=0;for(t=1;t<=t_boost;t++) yp[i]+=cbst[t]*ypt[t*nGivenData+i];
	    yp[i]/=csum;
	  }
#endif

//	  double csum=0;
//	  for(t=1;t<=t_boost;t++) csum+=cbst[t];
//	  for(i=0;i<nGivenData;i++){
//	    yp[i]=0;
//	    for(t=1;t<=t_boost;t++) yp[i]+=cbst[t]*ypt[t*nGivenData+i];
//	    yp[i]/=csum;
//	  }
	}
	double E2sum,e2;
	int ntrain=0;
	Ltst=E2sum=0;
	for(i=0;i<nGivenData;i++){
	  //	  e2 = square(yp[i]-yt[i]);//for entro yp
	  e2 = square(ypt[t_boost*nGivenData+i]-yt[i]);//for entro ypt
	  Ltst+=e2;
	  if(nc[i]==1){
#ifdef Bias24H0
	    E2sum+=bias2[i];
#else
	    E2sum+=e2;
#endif
	    ntrain++;
	  }
	}
	Ltst/=nGivenData;
	Lib=E2sum/ntrain;
	Lemp=Lib;
	double pi;
	H0=0;
	for(i=0;i<nGivenData;i++){
	  if(nc[i]==1){
	    //	    pi=square(yp[i]-yt[i])/E2sum;//for entro yp
#ifdef Bias4H0
	    pi=bias2[i]/E2sum+1e-20;
#else
	    pi=square(ypt[t_boost*nGivenData+i]-yt[i])/E2sum+1e-20;//for entro ypt
#endif
	    H0-= pi*log(pi);
	  }
	}
	H0/=log(ntrain);
	int nfiles=t_boost,num=nGivenData;
	sprintf(fn_boost,"%s/tmp/BstGb:%s+m%da%gs%dN:%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t_boost);//
	fp=fopen(fn_boost,"w");
	fprintf(fp,"%.7e\n",cbst[t_boost]);
	for(i=0;i<nGivenData;i++) fprintf(fp,"%.7e %.7e %.7e %.7e %.7e %.7e\n",pbst[i],wbst[i],ypt[t_boost*nGivenData+i],pbst1[i],wbst1[i],bias2[t_boost*nGivenData+i]);//yp0-->ypt
	fclose(fp);
	fprintf(stderr,">>'%s' for boosting is created\n",fn_boost);
	{
	  sprintf(fn_boost,"%s/tmp/BstGb:%s+m%da%gs%dN:%d:%db%dpred0.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t_boost);//
	  fp=fopen(fn_boost,"w");
	  //	  fprintf(fp,"%.7e\n",cbst[t_boost]);
	  for(i=0;i<nGivenData;i++) fprintf(fp,"%.7e %d %.7e\n",ypt[t_boost*nGivenData+i],i,yt[i]);//yp0-->ypt
	  fclose(fp);
	  sprintf(fn_boost,"%s/tmp/BstGb:%s+m%da%gs%dN:%d:%db%dis.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NStep,t_boost);//
	  fp=fopen(fn_boost,"w");
	  fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",0,0.0,nGivenData,0.0,0,0,nGivenData,0.0,0.0);
	  fclose(fp);
	}
	fp=fopen("./loss+.dat","w");
	fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0);
	fprintf(stderr,"#Method:BagwithoutVal boost%d %.7e %.7e %.7e #Ltst,Lib,H0\n",t_boost,Ltst,Lib,H0);
	fclose(fp);
	//	fprintf(stderr,"#Boost#%13.7e %.3e %1.3e %.3e %d %d %.3e %.3e #Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp ntrain/n=%.3f.\n",Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,((double)ntrainall/(ntrainall+ntestall)));
      }
    }//closing Boost==GbBoost
    ////
    else if(Boost==EmBoost){
      if(BAGGING==BAGwithVal){
	net_printf(strx1("Doing '%s'\n",cmd2));
	system(cmd2);//exec meanpred for ensemble
	sprintf(cmd,"cat %s | awk '{print $%d}' > %s",fn_bagging,K+1,fn_target);system(cmd);
	sprintf(&cmd[strlen(cmd)]," ;./pred2y_ts predict+.dat %s > /dev/null",fn_target);system(cmd);
	//fprintf(stderr,"Prediction saved in '%s', Loss in'./loss+.dat'. Bagging with method=%s. \n",fn_target,argv[2]);
      }
      else if(BAGGING==BAGwithoutVal){
	double LD=0,Ltst,Lvar=0,Lvarval=0,Lvar0=0,Lib,Lemp,H0;
	int nfiles,num;
	double *yptest=(double*)malloc(sizeof(double)*nGivenData);
	double *yptrain=(double*)malloc(sizeof(double)*nGivenData);
	double *yptotal=(double*)malloc(sizeof(double)*nGivenData);
	int *ntrain=(int*)malloc(sizeof(int)*nGivenData);
	int *ntest=(int*)malloc(sizeof(int)*nGivenData);
	int ntrainall,ntestall;
	double err2;
	double UC,ucn;//uncertainty
	double *uc=(double*)malloc(sizeof(double)*nGivenData);
	int nGivenData1=0,nGivenData0=0;
	//      double kc=2; //kc=coverage factor for extended uncertainty
	///////////////
	Ltst=Lib=Lemp=0;
	for(n=0;n<nGivenData;n++){
	  yptest[n]=yptrain[n]=ntrain[n]=ntest[n]=0;
	  for(j=0;j<nFolds;j++){
	    if(nc[j*nGivenData+n]!=1){//not in Dtrain
	      ntest[n]++;
	      yptest[n]+=yp0[j*nGivenData+n];
	      //	    yrtest[n]+=yr[j*nGivenData+n];
	    }
	    else{
	      ntrain[n]++;
	      yptrain[n]+=yp0[j*nGivenData+n];
	      //	    yrtrain[n]+=yr[j*nGivenData+n];
	      //	    Lemp+=square(yp0[j*nGivenData+n]-yt[n]);
	    }
	    yptotal[n]+=yp0[j*nGivenData+n];
	    //	  Lemp+=square(yp0[j*nGivenData+n]-yt[n]);
	  }
	  //////
	  if(ntest[n]>0){
	    nGivenData0++;
	    yptest[n]/=(ntest[n]);//avoiding zero division?
	    //	  yrtest[n]/=(ntest[n]);//avoiding zero division?
	    err=yptest[n]-yt[n];
	    err2=err*err;
	    Ltst+= err2;
#ifdef ERR2BOOST
	    if(t_boost>=0) ebst[n]=err2;
#else
	    if(t_boost>=0) ebst[n]=err;
#endif
	  }
	  else{
	    fprintf(stderr,"ntest=0 for x[%d] >> Increase the number of folds or reduce the alpha.\n",n);
	    if(t_boost>=0) ebst[n]=ebst0[n];//
	    //	  yptest[n]=yt[n];
	  }
	  yptotal[n]/=nFolds;
	  if(ntrain[n]>0){
	    nGivenData1++;
	    yptrain[n]/=(ntrain[n]);//avoiding zero division?
	    //	  yrtrain[n]/=(ntrain[n]);//avoiding zero division?
	    Lib+=square(yptrain[n]-yt[n]);
	    Lemp+=square(yptotal[n]-yt[n]);
	  }
	  else{
	    fprintf(stderr,"ntrain=0 for x[%d] >> Increase the number of folds or increase the alpha.\n",n);
	    //	  yptrain[n]=yt[n];
	  }
	}
	Ltst/=nGivenData0;
	Lib/=nGivenData1;
	Lemp/=nGivenData1;
	{
	  double pn,E2sum=Lemp*nGivenData1;
	  H0=0;
	  for(n=0;n<nGivenData;n++){
	    if(ntrain[n]>0){
	      pn=square(yptotal[n]-yt[n])/E2sum+1e-20;
	      H0-= pn*log(pn);
	    }
	  }
	  H0/=log(nGivenData1);
	}
	
	if(t_boost>=0){
	  //	if(1==0 && t_boost<=1){//remove data with ebst > rethresh_boost
	  //	  double psum=0;
	  //	  for(n=0;n<nGivenData;n++){
	  //	    if(ebst[n]/Ltst>rethresh_boost) pbst[n]=0.0;
	  //	    psum+=pbst[n];
	  //	  }
	  //	  for(n=0;n<nGivenData;n++) pbst[n]/=psum;
	  //	}
	  //new Ltst try!!!!!060709
	  //	double a=0.0,a1=1.-a;//good
	  double a=1.0,a1=1.-a;//
	  //double a=0.5,a1=1.-a;
	  //	double a=0.9,a1=1.-a;
	  //	double a=0.3,a1=1.-a;
	  //	double a=0.1,a1=1.-a;
	  int nGD=0;
	  Ltst=0;
	  for(n=0;n<nGivenData;n++){
	    if(pbst[n]>pthresh){
	      //	    {//entropy as error <<---trial 060810
	      //	      err=(yptotal[n]-yt[n]);//
	      //	      err2=err*err/Lemp;
	      //	      ebst[n]=-err2*log(err2);//
	      //	    }
	      {//good060810
		//	      err=(yptotal[n]-yt[n]);//
		err=a*(yptest[n]-yt[n])+a1*(yptrain[n]-yt[n]);//average??
		//	  err=((yptest[n]-yt[n])+(yptrain[n]-yt[n]))/2.;//average??
		//	  err=(yptrain[n]-yt[n]);//
		err2=err*err;
		ebst[n]=err2;
	      }
	      //	    Ltst+=err2;
	      Ltst+=square(yptest[n]-yt[n]);
	      nGD++;
	    }
	  }
	  Ltst/=nGD;
	}
	
	Lvar0=Lvar=ntrainall=ntestall=0;
	UC=ucn=0;
	for(n=0;n<nGivenData;n++){
	  uc[n]=0;
	  for(j=0;j<nFolds;j++){
	    if(nc[j*nGivenData+n]!=1){//not in Dtrain
	      err2=square(yp0[j*nGivenData+n]-yptest[n]);
	      Lvar+=err2;
	      uc[n]+=err2;
	    }
	    else{
	      Lvar0+=square(yp0[j*nGivenData+n]-yptrain[n]);
	    }
	  }
	  if(ntest[n]>2){
	    //	  uc[n]/=uc[n]/ntest[n]/(ntest[n]-1);//uncertainty
	    uc[n]/=(ntest[n]*(ntest[n]-1));//uncertainty
	    UC+=uc[n];//uncertainty
	    //	  UC+=uc[n]/ntest[n]/(ntest[n]-1);//uncertainty
	    ucn++;
	  }
	  ntestall+=ntest[n];
	  ntrainall+=ntrain[n];
	}
	UC/=ucn;//average uncertainty
	
	Lvar/=(ntestall);
	Lvar0/=(ntrainall);
	
	///
	nfiles=nFolds;
	num=nGivenData;
	LD=UC;
	if(t_boost>=0){//boosting
	  //???	for(n=0;n<nGivenData;n++) if((ebst[n]-=uc[n])<0) ebst[n]=0; //new Ltst-uncertainty=Lbias+err
	  //      double UCB;//Uncertainty for the real ensemble with nFolds
	  //      UCB=UC*(ntestall/nGivenData)/nFolds;
	  //	sprintf(fn_boost,"%s/tmp/%s+m%da%.2fs%dF%dN:%d-%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,nFolds,NN[0],NN[nens-1],NStep,t_boost);//
	  sprintf(fn_boost,"%s/tmp/BstEM:%s+m%da%gs%dN:%d-%d:%db%d.dat",rdir0,fnbody(fn_givendata,pr),method,alpha,seed,NN[0],NN[nens-1],NStep,t_boost);//
	  fp=fopen(fn_boost,"w");
	  for(n=0;n<nGivenData;n++) fprintf(fp,"%e %e %e %+e %+e %+e uc%eLtst%eUC%entrain%dntest%dn%d\n",pbst[n],pbst0[n],pbst1[n],ebst[n],ebst0[n],ebst1[n],uc[n],ebst[n],UC,ntrain[n],ntest[n],n);//
	  //	for(n=0;n<nGivenData;n++) fprintf(fp,"%e %e %e %+e %+e %+e uc%e Ltst%e UC%e\n",pbst[n],pbst0[n],pbst1[n],ebst[n],ebst0[n],ebst1[n],uc[n],ebst[n],UC);//
	  //for(n=0;n<nGivenData;n++) fprintf(fp,"%e %e %e %+e %+e %+e uc%e Ltst%e UC%e\n",pbst[n],pbst0[n],pbst1[n],ebst[n]-uc[n],ebst0[n],ebst1[n],uc[n],ebst[n],UC);//NG?
	  //for(n=0;n<nGivenData;n++) fprintf(fp,"%e %e %e %+e %+e %+e uc%e Ltst%e UC%e\n",pbst[n],pbst0[n],pbst1[n],ebst[n]-UC,ebst0[n],ebst1[n],uc[n],ebst[n],UC);
	  //for(n=0;n<nGivenData;n++) fprintf(fp,"%.7f %.7f %.7f %+e %+e %+e\n",pbst[n],pbst0[n],pbst1[n],ebst[n],ebst0[n],ebst1[n]);
	  //for(n=0;n<nGivenData;n++) fprintf(fp,"%e %e %e %+e %+e %+e\n",pbst[n],pbst0[n],pbst1[n],ebst[n],ebst0[n],ebst1[n]);
	  fclose(fp);
	  fprintf(stderr,">>'%s' for boosting is created\n",fn_boost);
	}
	//      Lvar0=Ltst-kc*UC; //kc=coverage factor for extended uncertainty
	//      Lvar=Ltst+kc*UC;//
	
	//      Lemp/=nGivenData1;
	//      Lemp=Entro;//
	
	fp=fopen("./loss+.dat","w");
	fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,H0);
	fclose(fp);
	fprintf(stderr,"#BAGGING2#%13.7e %.3e %1.3e %.3e %d %d %.3e %.3e #Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp ntrain/n=%.3f.H0=%e\n",Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,((double)ntrainall/(ntrainall+ntestall)),H0);
      }//closing of if(BAGGING==BAGwithoutVal)
    }
    ///////////////
    //    else if(BAGGING==BAGwithoutVal){
    //      double LD=0,Ltst,Lvar,Lvarval,Lvar0,Lib,Lemp;
    //      int nfiles,num;
    //      double *yptest=(double*)malloc(sizeof(double)*nGivenData);
    //      double *yptrain=(double*)malloc(sizeof(double)*nGivenData);
    //      int *ntrain=(int*)malloc(sizeof(int)*nGivenData);
    //      int *ntest=(int*)malloc(sizeof(int)*nGivenData);
    //      int ntrainall,ntestall;
    //      double err2;
    //      double UC,ucn;//uncertainty
    //      double *uc=(double*)malloc(sizeof(double)*nGivenData);
    //      //      double kc=2; //kc=coverage factor for extended uncertainty
    //      Ltst=Lib=Lemp=0;
    //
    //      for(n=0;n<nGivenData;n++){
    //	yptest[n]=yptrain[n]=ntrain[n]=ntest[n]=0;
    //	for(j=0;j<nFolds;j++){
    //	  if(nc[j*nGivenData+n]!=1){
    //	    ntest[n]++;
    //	    yptest[n]+=yp0[j*nGivenData+n];
    //	  }
    //	  else{
    //	    ntrain[n]++;
    //	    yptrain[n]+=yp0[j*nGivenData+n];
    //	    Lemp+=square(yp0[j*nGivenData+n]-yt[n]);
    //	  }
    //	}
    //	if(ntest[n]>0){
    //	  yptest[n]/=(ntest[n]);//avoiding zero division?
    //	  Ltst+=square(yptest[n]-yt[n]);
    //	  //	  if(n<10) fprintf(stderr,"%d ntest=%d, Ltst=%e\n",n,ntest[n],Ltst);//for check
    //	}
    //	else{
    //	  fprintf(stderr,"ntest=0 >> Increase the number of folds or reduce the alpha.\n");
    //	  //	  yptest[n]=yt[n];
    //	}
    //	if(ntrain[n]>0){
    //	  yptrain[n]/=(ntrain[n]);//avoiding zero division?
    //	  Lib+=square(yptrain[n]-yt[n]);
    //	}
    //	else{
    //	  fprintf(stderr,"ntrain=0 >> Increase the number of folds or increase the alpha.\n");
    //	  //	  yptrain[n]=yt[n];
    //	}
    //	//check
    //	//fprintf(stderr,"%d)Lib=%e,Ltst=%e,yptrain=%e,yptest=%e,yt=%e,ntrain=%d,ntest=%d.\n",n,Lib,Ltst,yptrain[n],yptest[n],yt[n],ntrain[n],ntest[n]);
    //      }
    //      Lib/=nGivenData;
    //      Ltst/=nGivenData;
    //      
    //      Lvar0=Lvar=ntrainall=ntestall=0;
    //      UC=ucn=0;
    //      for(n=0;n<nGivenData;n++){
    //	uc[n]=0;
    //	for(j=0;j<nFolds;j++){
    //	  if(nc[j*nGivenData+n]!=1){
    //	    err2=square(yp0[j*nGivenData+n]-yptest[n]);
    //	    Lvar+=err2;
    //	    uc[n]+=err2;
    //	  }
    //	  else{
    //	    Lvar0+=square(yp0[j*nGivenData+n]-yptrain[n]);
    //	  }
    //	}
    //	if(ntest[n]>2){
    //	  UC+=uc[n]/ntest[n]/(ntest[n]-1);
    //	  ucn++;
    //	}
    //	ntestall+=ntest[n];
    //	ntrainall+=ntrain[n];
    //      }
    //      UC/=ucn;
    //      Lvar/=(ntestall);
    //      Lvar0/=(ntrainall);
    //      ///
    //      nfiles=nFolds;
    //      num=nGivenData;
    //      LD=UC;
    //      //      Lvar0=Ltst-kc*UC; //kc=coverage factor for extended uncertainty
    //      //      Lvar=Ltst+kc*UC;//
    //      fp=fopen("./loss+.dat","w");
    //      fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %13.7e %d %d %.7e %.7e #LD Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp\n",LD,Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp);
    //      fclose(fp);
    //      fprintf(stderr,"#BAGGING2#%13.7e %.3e %1.3e %.3e %d %d %.3e %.3e #Ltst Lvar Lvarval Lvar0 nfiles num Lib Lemp ntrain/n=%.3f.\n",Ltst,Lvar,Lvarval,Lvar0,nfiles,num,Lib,Lemp,((double)ntrainall/(ntrainall+ntestall)));
    //    }//closing of if(BAGGING==2)
    else if(BAGGING==NoBAG){
      for(n=0;n<nFolds;n++){
	msesum+=mse[n];
	Lhatmean+=Lhat[n];
      }
      msemean=msesum/nFolds;
      Lhatmean /=nFolds;
      for(n=0;n<nFolds;n++){
	if((err=mse[n]-msemean)>0){
	  msestdp+=err*err;
	  n_msestdp++;
	}
	else {
	  msestdm+=err*err;
	  n_msestdm++;
	}
      }
      if(n_msestdp>0) msestdp=sqrt(msestdp/n_msestdp); else msestdp=0;
      if(n_msestdm>0) msestdm=sqrt(msestdm/n_msestdm); else msestdm=0;
      int nn;
      for(nn=0;nn<nFolds;nn++){
	n_msesum+=n_mse[nn];
	msetrainsum+=msetrain[nn];
	n_msetrainsum+=n_msetrain[nn];
	n_cells2sum+=n_cells2[nn];//
 	fprintf(stderr,"%3d %13.7e %13.7e %3d %3d %13.7e #n Ltest Lib?(see ./loss+.dat) bestIT Nr Lhat?+\n",
		nn,mse[nn],msetrain[nn],is[nn],n_cells2[nn],Lhat[nn]);//modified 051016
      }
      msetrainmean=msetrainsum/nFolds;
      fprintf(stderr,"%10.5e %8.3e %8.3e %10.5e %d %d %3d %3d %10.5e #<Ltest Lstd Lstd2 Lib?(see ./loss+.dat) n nFolds N Nr Lhat>\n",
	      msemean,//msesum/n_msesum,
	      msestdp,msestdm,
	      msetrainmean,//msetrainsum/n_msetrainsum,
	      nGivenData,
	      nFolds,nCells,(int)n_cells2sum/nFolds,
	      Lhatmean);
      fp=fopen("./tmp/rsresult.dat","w+");
      fprintf(fp,"%13.7e %13.7e %13.7e %13.7e %d %d %3d %3d %d %13.7e #Ltest Lstd Lstd2 Lib n nFolds N Nr a*n Lhat\n",
	      msemean,//msesum/n_msesum,
	      msestdp,msestdm,
	      msetrainmean,//msetrainsum/n_msetrainsum,
	      nGivenData,
	      nFolds,nCells,
	      (int)n_cells2sum/nFolds,
	      (int)(nGivenData-meannTestData+0.5),
	      Lhatmean);
      fclose(fp);
      //      fprintf(stderr,"method=%s, <nTrainData>=%d, <nTestData>=%d. E:%d.\n",argv[2],(int)(n_msetrainsum/nFolds),(int)(n_msesum/nFolds),ERR);
      fprintf(stderr,"method=%s, <nTrainData>=%d, <nTestData>=%d.\n",argv[2],(int)(n_msetrainsum/nFolds),(int)(n_msesum/nFolds));
    }
  }
  return(0);
}

