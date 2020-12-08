/*
 * am.c 
 * usage:
 * (1) AM am;
 * (2) init_AM(&am,3,2);
 * (3) calc_AM(am);
 */
#include "am.h"

#define am_lambda 1.0
#ifdef FLT
#define ZERO_AM 1.0e-5	    // 零 for float
#else
#define ZERO_AM 1.0e-10	    // 零 for double
//#define ZERO_AM 1.0e-11	    // 零 for double
//#define ZERO_AM 1.0e-5	    // 零 for double
#endif

/***********************************************************************
 * LMS Algorithm
 ***********************************************************************/
#if AM_VER == 0

void init_AM(AM *q, int nx, int ny)
{
  int i,j;
  q[0].nx=nx;
  q[0].ny=ny;
  q[0].r=0.002;
  q[0].M=(FLOAT **)malloc(ny*sizeof(FLOAT));
  q[0].y=(FLOAT *)malloc(ny*sizeof(FLOAT));
  q[0].x=(FLOAT *)malloc(nx*sizeof(FLOAT));
  for(i=0;i<ny;i++){
    q[0].M[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
    for(j=0;j<nx;j++) q[0].M[i][j]=0;
  }
  return;
}

void free_AM(AM *q)
{
  int i;

  for(i=0;i<q[0].ny;i++){
    free(q[0].M[i]);
  }
  free(q[0].M);
  free(q[0].x);
  free(q[0].y);
  return;
} 

void calc_AM(AM q)
{
  FLOAT err;
  int i,j;
  for(i=0;i<q.ny;i++){
    for(err=q.y[i],j=0;j<q.nx;j++) err -= (q.M[i][j]*q.x[j]);
    for(j=0;j<q.nx;j++){
      q.M[i][j]+=2.*q.r*err*q.x[j];
    }
  }
  return;
}

/***********************************************************************
 * RLS Algorithm
 ***********************************************************************/
#elif AM_VER == 1

void init_AMdata(AM *q)
{
  int i,j,nx,ny;
  nx=q[0].nx;
  ny=q[0].ny;
  for(i=0;i<ny;i++){
    for(j=0;j<nx;j++) q[0].M[i][j]=0;
  }
  for(i=0;i<nx;i++){
    for(j=0;j<nx;j++){
      //      if(i==j) q[0].P[i][j]=1.0e7;
      //      if(i==j) q[0].P[i][j]=1.0e5;
      if(i==j) q[0].P[i][j]=10000.0;//original
      else q[0].P[i][j]=0;
    }
  }
  return;
}

void init_AM(AM *q, int nx, int ny)
{
  int i;
  q[0].nx=nx;
  q[0].ny=ny;
  q[0].M=(FLOAT **) malloc(ny*sizeof(FLOAT));
  q[0].P=(FLOAT **) malloc(nx*sizeof(FLOAT));
  q[0].x=(FLOAT *) malloc(nx*sizeof(FLOAT));
  q[0].y=(FLOAT *) malloc(ny*sizeof(FLOAT));
  for(i=0;i<ny;i++){
    q[0].M[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
  }
  for(i=0;i<nx;i++){
    q[0].P[i]=(FLOAT *) malloc(nx*sizeof(FLOAT));
  }
  init_AMdata(q);
  return;
}

void free_AM(AM *q)
{
  int i;
  for(i=0;i<q[0].nx;i++){
    free(q[0].P[i]);
  }
  for(i=0;i<q[0].ny;i++){
    free(q[0].M[i]);
  }
  free(q[0].P);
  free(q[0].M);
  free(q[0].x);
  free(q[0].y);
  return;
} 

void calc_AM(AM q)
{
  FLOAT err;
  int i,j;
  FLOAT *p,*g,xPx1;

  p=(FLOAT *) malloc(q.nx*sizeof(FLOAT));
  g=(FLOAT *) malloc(q.nx*sizeof(FLOAT));
  
  for(i=0;i<q.nx;i++)for(g[i]=0,j=0;j<q.nx;j++) g[i]+= q.P[i][j]*q.x[j];
  for(xPx1=1.0,i=0;i<q.nx;i++) xPx1 += q.x[i]*g[i];
  for(i=0;i<q.nx;i++) p[i] = g[i]/xPx1;

  for(i=0;i<q.nx;i++){
    for(j=0;j<q.nx;j++){
      q.P[i][j] -= p[i]*g[j];
    }
  }
  for(i=0;i<q.ny;i++){
    for(err=q.y[i],j=0;j<q.nx;j++) err -= (q.M[i][j]*q.x[j]);
    for(j=0;j<q.nx;j++) q.M[i][j] += err*p[j];
  }
  /*  for(i=0;i<q.nx;i++){
    printf("[%d]g=%f,x=%f,w=%f,p=%f\n",i,g[i],q.x[i],q.M[0][i],p[i]);
    }*/
  
  free(p);
  free(g);
  return;
}

/***********************************************************************
 * 逐次的自己回帰モデル
 * Cf. Kohonen(1977)Associative Memory,Sect.3.3.5
 * 中谷和夫訳、コホネン(1977)連想記憶、3.3.5節(pp.164-165)
 ***********************************************************************/
#elif AM_VER == 2
//int dimdim=0;

void init_AMdata(AM *q)
{
  int i,j,nx,ny;

  nx=q[0].nx;
  ny=q[0].ny;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++) q[0].M[j][i]=0.0;
    for(j=0;j<nx;j++){
      q[0].Q[i][j]=0;
      if(i==j) q[0].P[i][j]=1.0;
      else q[0].P[i][j]=0;
    }
  }
  return;
}

/*
 * init_AM(&am) is a standard call of this subroutine. 
 * Here, *q is the pointer to the variable am with struct AM). 
 */
void init_AM(AM *q, int nx, int ny)
{
  int i,j;

  q[0].nx=nx;
  q[0].ny=ny;
  q[0].x=(FLOAT *)malloc(nx*sizeof(FLOAT));
  q[0].y=(FLOAT *)malloc(ny*sizeof(FLOAT));
  q[0].M=(FLOAT **)malloc(ny*sizeof(FLOAT *));
  q[0].P=(FLOAT **)malloc(nx*sizeof(FLOAT *));
  q[0].Q=(FLOAT **)malloc(nx*sizeof(FLOAT *));
  for(i=0;i<ny;i++){
    q[0].M[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
  }
  for(i=0;i<nx;i++){
    q[0].P[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
    q[0].Q[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
  }
  init_AMdata(q);
  return;
}

void free_AM(AM *q)
{
  int i;
  for(i=0;i<q[0].nx;i++){
    free(q[0].P[i]);
    free(q[0].Q[i]);
  }
  for(i=0;i<q[0].ny;i++){
    free(q[0].M[i]);
  }
  free(q[0].P);
  free(q[0].Q);
  free(q[0].M);
  free(q[0].x);
  free(q[0].y);
  return;
} 

/*
 * recursive algorithm
 */
int calc_AM(AM q)
{
  FLOAT *h, *g,*p,yhat;
  int i,j,k;
  FLOAT hth,gtx1,e;

  p=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  h=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  g=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  for(i=0;i<q.nx;i++) for(h[i]=0,j=0;j<q.nx;j++) h[i] += q.P[i][j]*q.x[j];
  for(i=0;i<q.nx;i++) for(g[i]=0,j=0;j<q.nx;j++) g[i] += q.Q[i][j]*q.x[j];
  for(hth=0,i=0;i<q.nx;i++) hth += h[i]*h[i];
  if(hth<ZERO_AM){/*線形従属のとき*/
    for(gtx1=am_lambda,i=0;i<q.nx;i++) gtx1 += g[i]*q.x[i];
    for(i=0;i<q.nx;i++) p[i] = g[i]/gtx1;/**/
    for(i=0;i<q.nx;i++){
      for(j=0;j<q.nx;j++){
	q.Q[i][j] = (q.Q[i][j] - g[i]*p[j])/am_lambda;
      }
    }
  }
  else{/*線形独立のとき*/
    /*    ++dimdim;*/
    /*    printf("\nhth=%e,",hth);*/
    for(gtx1=am_lambda,i=0;i<q.nx;i++) gtx1 += g[i]*q.x[i];
    for(i=0;i<q.nx;i++) p[i] = h[i]/hth;
    for(i=0;i<q.nx;i++){
      for(j=0;j<q.nx;j++){
	q.P[i][j] -= h[i]*p[j];
	q.Q[i][j] = (q.Q[i][j]
		     +gtx1*p[i]*p[j]
		     -g[i]*p[j]
		     -p[i]*g[j])/am_lambda;
      }
    }
  }
  
  for(j=0;j<q.ny;j++){
    yhat=0;
    for(i=0;i<q.nx;i++) yhat += (q.M[j][i]*q.x[i]);
    e=q.y[j]-yhat;
    for(i=0;i<q.nx;i++) q.M[j][i]+=e*p[i];
  }
  free(h);
  free(g);
  free(p);
  /*  return(hth<ZERO_AM);*/
  return(hth>=ZERO_AM);
}

/***********************************************************************
 * UD フィルタ
 * 飯國，適応信号処理アルゴリズム(2000),培風館,p125,or p205
 ***********************************************************************/
#elif AM_VER == 3

void init_AMdata(AM *q)
{
  int i,j,nx,ny;
  nx=q[0].nx;
  ny=q[0].ny;
  for(i=0;i<ny;i++){
    for(j=0;j<nx;j++) q[0].M[i][j]=0;
  }
  for(i=0;i<nx;i++){
    q[0].D[i]=1000.0;
    for(j=0;j<nx;j++){
      if(i==j) q[0].U[i][j]=1.0;
      else q[0].U[i][j]=0;
    }
  }
  return;
}

void init_AM(AM *q, int nx, int ny)
{
  int i,j;
  q[0].nx=nx;
  q[0].ny=ny;
  q[0].M=(FLOAT **) malloc(ny*sizeof(FLOAT));
  q[0].U=(FLOAT **) malloc(nx*sizeof(FLOAT));
  q[0].D=(FLOAT *) malloc(nx*sizeof(FLOAT));
  q[0].x=(FLOAT *) malloc(nx*sizeof(FLOAT));
  q[0].y=(FLOAT *) malloc(ny*sizeof(FLOAT));
  for(i=0;i<ny;i++){
    q[0].M[i]=(FLOAT *)malloc(nx*sizeof(FLOAT));
  }
  for(i=0;i<nx;i++){
    q[0].U[i]=(FLOAT *) malloc(nx*sizeof(FLOAT));
  }
  init_AMdata(q);
  return;
}

void free_AM(AM *q)
{
  int i;
  for(i=0;i<q[0].nx;i++){
    free(q[0].U[i]);
  }
  for(i=0;i<q[0].ny;i++){
    free(q[0].M[i]);
  }
  free(q[0].U);
  free(q[0].M);
  free(q[0].D);
  free(q[0].x);
  free(q[0].y);
  return;
}

int calc_AM(AM q)
{
  FLOAT *f,*g,*alpha,*k,*v,*p,*mu;
  int i,j;
  FLOAT gtx1,uij,e;

  f=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  g=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  alpha=(FLOAT *)malloc((q.nx+1)*sizeof(FLOAT));
  k=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  v=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  mu=(FLOAT *)malloc(q.nx*sizeof(FLOAT));
  p=(FLOAT *)malloc(q.nx*sizeof(FLOAT));

  for(i=0;i<q.nx;i++)for(f[i]=0,j=0;j<q.nx;j++) f[i] += q.U[j][i]*q.x[j];
  for(i=0;i<q.nx;i++) g[i]=q.D[i]*f[i];
  alpha[0]=am_lambda;
  for(j=0;j<q.nx;j++){
    alpha[j+1]=alpha[j]+f[j]*g[j];
    q.D[j]=alpha[j]/alpha[j+1]*q.D[j]/am_lambda;
    v[j]=g[j];
    mu[j]=-f[j]/alpha[j];
    if(j>0){/*三角行列*/
      for(i=0;i<j;i++){
	uij=q.U[i][j];
	q.U[i][j] = uij + v[i]*mu[j];
	v[i] = v[i] + uij*v[j];
      }
    }
  }
  for(i=0;i<q.ny;i++){
    for(e=q.y[i],j=0;j<q.nx;j++) e -= (q.M[i][j]*q.x[j]);
    for(j=0;j<q.nx;j++) q.M[i][j] += e*v[j]/alpha[q.nx];;/*orig*/
  }

  free(f);
  free(g);
  free(alpha);
  free(k);
  free(v);
  free(mu);
  free(p);
  return(1);
}

/***********************************************************************
 * Levinson-Durbin
 ***********************************************************************/
#elif AM_VER == 4

void init_AM(AM *q, int p, int N)
{
  int i,j,p1;
  p1=p+1;
  q[0].p=p;
  q[0].N=N;
  q[0].aa=(FLOAT **) malloc(p1*sizeof(FLOAT));
  q[0].rho=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].r=(FLOAT *) malloc(p1*sizeof(FLOAT));
  for(i=0;i<p1;i++){
    q[0].aa[i]=(FLOAT *)malloc(p1*sizeof(FLOAT));
  }
  return;
}

void free_AM(AM *q)
{
  int i;
  for(i=0;i<q[0].p+1;i++){
    free(q[0].aa[i]);
  }
  free(q[0].aa);
  free(q[0].rho);
  free(q[0].r);
  return;
}

void calc_AM(AM q)
{
  FLOAT err;
  int i,j,m,p1,n;
  FLOAT *r,*sigma2,Dm1;

  p1=q.p+1;
  /*
  auto_cor(q.x,q.N);
  r=&x_re[0];*/
  r=q.r;
  for(m=0;m<p1;m++){
    r[m]=0;n=0;
    for(j=m;j<q.N;j++){
      r[m]+=q.x[j]*q.x[j-m];
      n++;
    }
    r[m]/=n;
  }/**/

  sigma2=(FLOAT *)malloc(p1*sizeof(FLOAT));
  sigma2[0]=r[0];
  for(m=0;m<q.p;m++){
    Dm1=r[m+1];
    for(i=1;i<=m;i++){ 
      Dm1 += q.aa[m][i]*r[m+1-i];
    }
    q.rho[m+1]=q.aa[m+1][m+1]=-Dm1/sigma2[m];
    sigma2[m+1]=sigma2[m]*(1.0-q.rho[m+1]*q.rho[m+1]);
    for(i=1;i<=m;i++){ 
      q.aa[m+1][i]=q.aa[m][i]+q.rho[m+1]*q.aa[m][m+1-i];
    }
  }
  free(sigma2);
  return;
}

/***********************************************************************
 * LSL
 ***********************************************************************/
#elif AM_VER == 5

void init_AM(AM *q, int p, int N)
{
  int i,j,p1;
  p1=p+1;
  q[0].p=p;
  q[0].N=N;
  q[0].R=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].Delta=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].f=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].r=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].F=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].theta=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].R1=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].r1=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].alpha=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].beta=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].n=(FLOAT *) malloc(1*sizeof(FLOAT));
  return;
}

void free_AM(AM *q)
{
  int i;
  free(q[0].Delta);
  free(q[0].f);
  free(q[0].R);
  free(q[0].r);
  free(q[0].F);
  free(q[0].theta);
  free(q[0].R1);
  free(q[0].r1);
  free(q[0].alpha);
  free(q[0].beta);
  return;
}

void calc_AM(AM q)
{
  int i,j,m,p,p1,n,ip;
  FLOAT *R,*Delta,*f,*r,*F;
  FLOAT *theta,*R1,*r1,*alpha,*beta;

  p=q.p;
  p1=p+1;
  Delta=q.Delta;
  f=q.f;
  F=q.F;
  theta=q.theta;
  R=q.R;
  r=q.r;
  R1=q.R1;
  r1=q.r1;
  alpha=q.alpha;
  beta=q.beta;

  if(q.n[0]==0){
    F[0]=0;
    for(m=0;m<=p-1;m++)  Delta[m+1]=0;
  }
  f[0]=r[0]=q.x[0];
  F[0]=R[0]= am_lambda * F[0]+q.x[0]*q.x[0];/**/
  theta[0]=1;
  if(q.n[0]<q.p) ip=q.n[0]; else ip=q.p;
  for(m=0;m<=ip-1;m++){
    Delta[m+1]=am_lambda *Delta[m+1]+f[m]*r1[m]/theta[m];
    if(fabs(R1[m])<ZERO_AM){
      printf("R1[%d]=%e\n",m,R1[m]);
      R1[m]=ZERO_AM;
    }
    if(fabs(F[m])<ZERO_AM){
      printf("F[%d]=%e\n",m,F[m]);
      F[m]=ZERO_AM;
    }
    F[m+1]    =F[m] -Delta[m+1]*Delta[m+1]/R1[m];
    R[m+1]    =R1[m]-Delta[m+1]*Delta[m+1]/F[m];
    alpha[m+1]=-Delta[m+1]/R1[m];
    beta[m+1] =-Delta[m+1]/ F[m];
    f[m+1]    =f[m] +alpha[m+1]*r[m];
    r[m+1]    =r1[m]+beta[m+1] *f[m];
    theta[m+1]=theta[m]-r1[m]*r1[m]/R1[m];
  }
  for(m=0;m<=p;m++){
    R1[m]=R[m];/**/
    r1[m]=r[m];
  }
  q.n[0]++;
  return;
}

void calc_AM0(AM q)
{
  int i,j,m,p,p1,n,ip;
  FLOAT *R,*Delta,*f,*r,*F;
  FLOAT *theta,*R1,*r1,*alpha,*beta;

  p=q.p;
  p1=p+1;
  Delta=q.Delta;
  f=q.f;
  F=q.F;
  theta=q.theta;
  R=q.R;
  r=q.r;
  R1=q.R1;
  r1=q.r1;
  alpha=q.alpha;
  beta=q.beta;

  F[0]=0;
  for(m=0;m<=p-1;m++)  Delta[m+1]=0;
  for(n=0;n<q.N;n++){
    f[0]=r[0]=q.x[n];
    F[0]=R[0]= am_lambda * F[0]+q.x[n]*q.x[n];/**/
    theta[0]=1;
    if(n<q.p) ip=n;else ip=q.p;
    for(m=0;m<=ip-1;m++){
      Delta[m+1]=am_lambda *Delta[m+1]+f[m]*r1[m]/theta[m];
      if(fabs(R1[m])<ZERO_AM){
	printf("R1[%d]=%e\n",m,R1[m]);
	R1[m]=ZERO_AM;
	/*	if(R1[m]>=0) R1[m]=ZERO_AM;
		else R1[m]=-ZERO_AM;*/
      }
      if(fabs(F[m])<ZERO_AM){
	printf("F[%d]=%e\n",m,F[m]);
	F[m]=ZERO_AM;
	/*if(F[m]>=0) F[m]=ZERO_AM;
	  else F[m]=-ZERO_AM;*/
      }
      F[m+1]    =F[m] -Delta[m+1]*Delta[m+1]/R1[m];
      R[m+1]    =R1[m]-Delta[m+1]*Delta[m+1]/F[m];
      alpha[m+1]=-Delta[m+1]/R1[m];
      beta[m+1] =-Delta[m+1]/ F[m];
      f[m+1]    =f[m] +alpha[m+1]*r[m];
      r[m+1]    =r1[m]+beta[m+1] *f[m];
      theta[m+1]=theta[m]-r1[m]*r1[m]/R1[m];
    }
    /*for(m=0;m<=p;m++){
      if(fabs(R[m])<ZERO_AM){
	printf("R[%d]=%e\n",m,R[m]);
	R1[m]=ZERO_AM;
	if(R[m]>=0) R1[m]=1; else R1[m]=-1;
	if(R[m]>=0) R1[m]=ZERO_AM; else R1[m]=-ZERO_AM;
      }
      else R1[m]=R[m];*/
      R1[m]=R[m];/**/
      r1[m]=r[m];
      /*}*/
      /*if(F[0]<ZERO_AM) F[0]=ZERO_AM;*/
      /*printf("f=%f %f %f,ip=%d,%d\n",f[0],f[1],f[2],ip,n);*/
  }
  return;
}

/***********************************************************************
 * 正規化LSL
 ***********************************************************************/
#elif AM_VER == 6

void init_AM(AM *q, int p, int N)
{
  int i,j,p1;
  p1=p+1;
  q[0].p=p;
  q[0].N=N;
  q[0].rho=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].rho1=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].f=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].r=(FLOAT *) malloc(p1*sizeof(FLOAT));
  q[0].r1=(FLOAT *) malloc(p1*sizeof(FLOAT));
  return;
}

void free_AM(AM *q)
{
  int i;
  free(q[0].rho);
  free(q[0].rho1);
  free(q[0].f);
  free(q[0].r);
  free(q[0].r1);
  return;
}

void calc_AM(AM q)
{
  int i,j,m,p,p1,n,ip;
  FLOAT *rho1,*rho,*f,*r,*r1,F,a1,a2,b1,b2;
  p=q.p;
  p1=p+1;
  rho=q.rho;
  rho1=q.rho1;
  f=q.f;
  r=q.r;
  r1=q.r1;

  F=0;
  for(n=0;n<q.N;n++){
    F= am_lambda * F+q.x[n]*q.x[n];
    if(F<ZERO_AM) F=ZERO_AM;
    f[0]=r[0]=q.x[n]/sqrt(F);
    for(m=0;m<=p-1;m++) {
      rho[m+1]=0;
    }
    if(n<q.p) ip=n;else ip=q.p;
    for(m=0;m<=ip-1;m++){
      if(f[m]>=1) f[m]=1-ZERO_AM; else if(f[m]<=-1) f[m]=-1+ZERO_AM;
      if(r1[m]>=1) r1[m]=1-ZERO_AM; else if(r1[m]<=-1) r1[m]=-1+ZERO_AM;
      rho[m+1]=sqrt(1.0-f[m]*f[m])*sqrt(1.0-r1[m]*r1[m])*rho1[m+1]-f[m]*r1[m];
      a1=sqrt(1.0-rho[m+1]*rho[m+1]);
      a2=sqrt(1.0-r1[m]*r1[m]);
      f[m+1]=(f[m]+rho[m+1]*r1[m+1])/a1/a2;
      b1=sqrt(1.0-rho[m+1]*rho[m+1]);
      b2=sqrt(1.0-f[m]*f[m]);
      r[m+1]=(r1[m]+rho[m+1]*f[m])/b1/b2;
    }
    for(m=0;m<=p;m++){
      r1[m]=r[m];
      rho1[m]=rho[m];
    }
/*    printf("2r%d=%f,%f,%f,%f\n",n,r1[0],r1[1],r[0],r[1]);*/
      /*    printf("f=%f %f %f,ip=%d,%d\n",f[0],f[1],f[2],ip,n);*/
  }
  return;
}

/***********************************************************************
 * FTF
 ***********************************************************************/
#elif AM_VER == 7

void init_AM(AM *q, int nx, int ny)
{
  int i,j,nx1;
  q[0].nx=nx;
  q[0].ny=ny;
  nx1=nx+1;
  q[0].w=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].gp=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].g=(FLOAT *) malloc((nx+2)*sizeof(FLOAT));
  q[0].newgp=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].alpha=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].beta=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].rr=(FLOAT *) malloc(sizeof(FLOAT));
  q[0].ff=(FLOAT *) malloc(sizeof(FLOAT));
  q[0].thetap=(FLOAT *) malloc(sizeof(FLOAT));
  q[0].x=(FLOAT *) malloc(nx1*sizeof(FLOAT));
  q[0].y=(FLOAT *) malloc(ny*sizeof(FLOAT));
  q[0].n=(FLOAT *) malloc(1*sizeof(FLOAT));
  q[0].n[0]=0;
  return;
}

void free_AM(AM *q)
{
  free(q[0].w);
  return;
}

void calc_AM(AM q)
{
  int i,j,m,p,p1,n,ip,t,T;
  FLOAT *x,*y,tmp;
  FLOAT be,e;
  FLOAT oldff, theta, br, r,bf,f,fff;
  p=q.nx;
  p1=p+1;
  x=q.x;
  y=q.y;

  n=q.n[0];
  if(n==0){
    q.ff[0]=x[0]*x[0];
    tmp=x[0];
    if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.w[0]=-y[0]/tmp;
    q.thetap[0]=1;
    q.n[0]++;
  }
  else if(n<=p){
    bf=x[0];
    for(i=1;i<=n-1; i++) bf = bf+x[i]*q.alpha[i];
    f=bf*q.thetap[0];
    q.ff[0]= am_lambda*q.ff[0];
    fff=q.ff[0]+bf*f;
    tmp=fff*q.thetap[0];
    if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.thetap[0]=q.ff[0]/tmp;
    tmp=q.ff[0];
    if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.newgp[1]=bf/tmp;
    for(i=2;i<=n;i++){
      q.newgp[i]=q.gp[i-1]+q.newgp[1]*q.alpha[i-1];
    }
    for(i=1;i<=n;i++) q.gp[i]=q.newgp[i];
    be=y[0];
    for(i=0;i<=n-1;i++) be=be+x[i]*q.w[i];
    e=be*q.thetap[0];
    tmp=x[n];
    if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.alpha[n]=-bf/tmp;
    q.w[n]=-be/tmp;
    if(n==p){
      q.rr[0]=x[p]*x[p]*q.thetap[0];
      for(i=1;i<=p;i++) q.beta[i]=-x[p]*q.thetap[0]*q.gp[i];
    }
    q.n[0]++;
  }
  else if(n>=p1){
    bf=x[0];
    for(i=1;i<=p;i++) bf=bf+x[i]*q.alpha[i];
    f=bf*q.thetap[0];
    oldff=q.ff[0];
    q.ff[0]=am_lambda*oldff+bf*f;
    tmp=q.ff[0];if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    theta=am_lambda*(oldff/tmp)*q.thetap[0];
    tmp=am_lambda*oldff;if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.g[1]=bf/tmp;
    for(i=2;i<=p+1;i++){
      q.g[i]=q.gp[i-1]+q.g[1]*q.alpha[i-1];
    }
    br=am_lambda*q.g[p+1]*q.rr[0];
    tmp=1-q.g[p+1]*theta*br;if(fabs(tmp)<ZERO_AM) tmp=ZERO_AM;
    q.thetap[0]=theta/tmp;
    r=br*q.thetap[0];
    q.rr[0]=am_lambda*q.rr[0]+br*r;
    for(i=1;i<=p;i++) q.alpha[i]=q.alpha[i]-f*q.gp[i];
    for(i=1;i<=p;i++) q.gp[i]=q.g[i]-q.g[p+1]*q.beta[i];
    for(i=1;i<=p;i++) q.beta[i]=q.beta[i]-r*q.gp[i];
    be=y[0];
    for(i=0;i<=p;i++) be=be+x[i]*q.w[i];
    e=be*theta;
    for(i=0;i<=p;i++) q.w[i]=q.w[i]-e*q.g[i+1];
    q.n[0]++;
  }
  printf("n=%f,x=%e,y=%e,w=%f,%f,%f,e=%f\n",q.n[0],q.x[0],q.y[0],q.w[0],q.w[1],q.w[2],e);
  return;
}

#endif /* AM_VER */
