/*
 *
 * 競合連想ネット Competitive Associative Networks (CAN2)
 *
 *
 * <学習モード>
 * L_MODE 1 alpha_i    を一定化  区分的2次関数近似の平均二乗誤差
 * L_MODE 2 S_i/n_i を一定化  区分的2次関数近似の平均二乗誤差
 * L_MODE 3 n_i     を一定化  発火回数
 * L_MODE 4 |y-y^| < errlimit
 *
 */
#include "my_plinn.h"

////////////////////////////
#define ZERO 1.0e-10
#define L_MODE 1
//#define E_CHECK
//#define alpha1 1.2
#define alpha1 1.1
#define alphaT1 5.0
//#define alphaT1 4.0
//#define alphaT1 8.0
//#define alphaT1 9.0
//#define alphaT1 15.0
//#define alphaT1 20.0 //good for 02train
//#define alphaT1 40.0
//#define alphaT1 2.0

#include "randoms.h"

#define LA(x,y1,x1,y2,x2) ((y1)+((x)-(x1))*((y2)-(y1))/((x2)-(x1))) //Linear Approximation

//0.95-(0.95-0.8)*(500-net->n_cells)/(500-100);
FLOAT alphat_e(int t,NET *net)
{
  return(1.-net->nentropy);
}
FLOAT alphat(int t)
{
  FLOAT a;
  a=alpha1/(1.+t/alphaT1);//
  //(n_train/n_cells)の関数?
  //02train(1,1)114ssp1.2939e-02(5.1940e-02)N500k2w0.1vm2
  //02train(1,2)24ssp1.2795e-02(5.1363e-02)N500k2w0.1vm2
  //02train(1,3)18ssp1.2177e-02(4.8881e-02)N500k2w0.1vm2
  //02train(1,4) 2ssp1.3949e-02(5.5995e-02)N500k2w0.1vm2
  //02train(1,5)50ssp1.0796e-02(4.3338e-02)N500k2w0.1vm2
  //02train(1,6) 2ssp1.3197e-02(5.2975e-02)N500k2w0.1vm2
  //02train(1,7)>10ssp1.1872e-02(4.7658e-02)N500k2w0.1vm2
  //02train(1,8)42ssp1.2867e-02(5.1651e-02)N500k2w0.1vm2

  //  a=1./t-1./1000.;//best
  // 52ssp1.1347e-02(4.5550e-02)N500k2w0.1vm2
  //a=1./t-1./10000.;//
  //02train 17ssp1.3171e-02(5.2873e-02)N500k2w0.1vm2
  //  a=1./t;//best52ssp1.1347e-02(4.5550e-02)N500k2w0.1vm2
  //02train 11ssp1.3323e-02(5.3483e-02)N500k2w0.1vm2
  //  a=1*exp(-t/70.);//02train8ssp1.4051e-02(5.6404e-02)N500k2w0.1vm2
  //  a=1*exp(-t/75.);//02train95ssp1.4747e-02(5.9197e-02)N500k2w0.1vm2
  //  a=1*exp(-t/65.);//02train28ssp1.3182e-02(5.2915e-02)N500k2w0.1vm2
  //
  //  a=0.999/t;//NG
  //    a=0.99*exp(-t/70.);//good 02train86ssp1.1720e-02(4.7050e-02)N500k2w0.1vm2
  //  a=0.99*exp(-t/75.);//??02train25ssp1.3589e-02(5.4550e-02)N500k2w0.1vm2
  //a=0.99*exp(-t/65.);//??02train88ssp1.4367e-02(5.7674e-02)N500k2w0.1vm2
  //a=5e-3;//soso
  //  a=1./t;
  //  a=1./t-1./100.;
  //if(a<1e-3) a=1e-3;
  if(a<1e-5) a=1e-5;
    //if(a<1e-20) a=1e-20;
  return(a);
}
#define min(a, b) ((a) < (b) ? (a) : (b))

// 外部変数
FLOAT entropy = 0.0;
FLOAT rate_w = 1.0;


/*====================================================================*
 * init_net
 * ネットを初期化
 *====================================================================*/
NET* init_net(NET *net) {
  char name[32] = "init_net";
  int i, k;
  int n_channels=net->k;
  {
    net->printf("================================================================\n");
    net->printf("	Generate Networks\n");
    net->printf("================================================================\n");
    net->printf("	<learning mode list>\n");
    net->printf(strx1("	%d :	online mode\n", ONLINE_MODE));
    net->printf(strx1("	%d :	batch mode\n", BATCH_MODE));
    net->printf("----------------------------------------------------------------\n");
    net->printf("* Number Of Units  : ");scanf1("%d", &net->n_cells);
    net->printf("* Number Of Compare Cells?   : "); scanf1("%d", &(net->n_compare));
    net->printf("* v_thresh in [0,1], v_min, vmin2: "); scanf3("%lf%d%d",&net->v_thresh,&net->vmin,&net->vmin2);
    if(net->vmin>net->n_cells) net->vmin=net->n_cells;
    net->printf("* Value Ratio for reinit         : "); scanf1("%lf", &net->v_ratio);
    net->printf("* Window Width?                  : "); scanf1("%lf", &net->width);
    net->printf("\n");
  }

  if (net == NULL) { net->printf("memory allocation error #1.\n"); exit(-1); }
  net->k = n_channels;
  if (net->n_cells < net->n_compare) net->n_compare = net->n_cells;
  net->t_ri = 0;

  net->cell = (CELL *)malloc(sizeof(CELL)*net->n_cells);
  if (net->cell == NULL) {
    net->printf(strx1("error: memory allocation error #?. (%s)\n", name));
    exit(-1);
  }
  net->w  = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);
  net->dw = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);
  for (i=0; i<net->n_cells; i++) {
    net->cell[i].w  = net->w[i]  = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    net->cell[i].dw = net->dw[i] = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    if (net->cell[i].w == NULL) { net->printf("memory allocation error #3.\n"); exit(-1); }
    for (k=0; k<n_channels; k++) net->cell[i].w[k] = 0.0;
    net->cell[i].v = 0.0;
    net->cell[i].S = 0.0;
    net->cell[i].alpha = 0.0;
    net->cell[i].enough = 1;
    init_AM(&(net->cell[i].am), n_channels+1, 1);
  }
  net->tau_E = (FLOAT)(net->n_cells*800);
  net->r_E = (FLOAT)exp(-1./net->tau_E);
  net->init=1;
  net->nEns=1;
  return(net);
}


NET* init_nets(NET *net) {
  char name[32] = "init_net";
  int i, k,j;
  int n_channels=net->k;
  int N1,N2,dN,nEns;
  int n_compare,vmin,vmin2;
  FLOAT v_thresh,v_ratio,width;

  net->printf("================================================================\n");
  net->printf("	Generate Networks\n");
  net->printf("================================================================\n");
  net->printf("	<learning mode list>\n");
  net->printf(strx1("	%d :	online mode\n", ONLINE_MODE));
  net->printf(strx1("	%d :	batch mode\n", BATCH_MODE));
  net->printf("----------------------------------------------------------------\n");
  net->printf("* Number Of Units (N1-Nn:dN for N=N1+i*dN for ensemble)  : ");
  {
    int ptn=0;
    char buff[20],*p;
    fgets(buff,20,stdin);
    for(p=buff;;){
      if(*p==0) break;
      if(ptn==0){if(*p=='-') {ptn=1;}}
      else if(ptn==1){ if(*p==':') {ptn=2;break;}}
      p++;
    }
    if(ptn==2) {sscanf(buff,"%d-%d:%d",&N1,&N2,&dN);;}
    else if(ptn==1) {sscanf(buff,"%d-%d",&N1,&N2);dN=1;}
    else if(ptn==0) {sscanf(buff,"%d",&N1);N2=N1,dN=1;}
    if(net!=NULL) free(net);
    nEns=(N2-N1)/dN+1;
    net=(NET *)malloc(sizeof(NET)*nEns);
  }
  net->printf("* Number Of Compare Cells?   : "); scanf1("%d", &(n_compare));
  net->printf("* v_thresh in [0,1], v_min, vmin2: "); scanf3("%lf%d%d",&v_thresh,&vmin,&vmin2);
  net->printf("* Value Ratio for reinit         : "); scanf1("%lf", &v_ratio);
  net->printf("* Window Width?                  : "); scanf1("%lf", &width);
  net->printf("\n");

  for(j=0;j<nEns;j++){
    net[j].n_cells=N1+j*dN;
    net[j].k = n_channels;
    net[j].n_compare=n_compare;
    net[j].v_thresh=v_thresh;
    net[j].vmin=vmin;
    net[j].vmin2=vmin2;
    net[j].v_ratio=v_ratio;
    net[j].width=width;
    if (net[j].n_compare>net[j].n_cells) net[j].n_compare = net[j].n_cells;
    net[j].t_ri = 0;
    
    net[j].cell = (CELL *)malloc(sizeof(CELL)*net[j].n_cells);
    if (net[j].cell == NULL) {
      net[j].printf(strx1("error: memory allocation error #?. (%s)\n", name));
      exit(-1);
    }
    net[j].w  = (FLOAT **)malloc(sizeof(FLOAT *)*net[j].n_cells);
    net[j].dw = (FLOAT **)malloc(sizeof(FLOAT *)*net[j].n_cells);
    for (i=0; i<net[j].n_cells; i++) {
      net[j].cell[i].w  = net[j].w[i]  = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
      net[j].cell[i].dw = net[j].dw[i] = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
      if (net[j].cell[i].w == NULL) { net[j].printf("memory allocation error #3.\n"); exit(-1); }
      for (k=0; k<n_channels; k++) net[j].cell[i].w[k] = 0.0;
      net[j].cell[i].v = 0.0;
      net[j].cell[i].S = 0.0;
      net[j].cell[i].alpha = 0.0;
      net[j].cell[i].enough = 1;
      init_AM(&(net[j].cell[i].am), n_channels+1, 1);
    }
    net[j].tau_E = (FLOAT)(net[j].n_cells*800);
    net[j].r_E = (FLOAT)exp(-1./net[j].tau_E);
    net[j].init=1;
    net[j].nEns=nEns;
  }
  return(net);
}


/*====================================================================*
 * show_net_parms
 * ネットの各パラメータを表示する
 *====================================================================*/
void show_net_parms(NET *net) {
  //char str[64];

  net->printf("================================================================\n");
  net->printf("	Show The Network Prameters\n");
  net->printf("================================================================\n");
  net->printf(strx1("> Learning Mode                   : %d\n", net->l_mode));
  net->printf(strx1("> Iteration Times                 : %d\n", net->i_times));
  net->printf(strx1("> Display Times With Learning     : %d\n", net->d_times));
  net->printf(strx1("> Number Of Storing Vectors       : %d\n", net->n_cells));
  net->printf(strx1("> Number Of Pattern Channels      : %d\n", net->k));
  net->printf(strx1("> Number Of Compare Cells         : %d\n", net->n_compare));
  net->printf(strx1("> Value Threshold                 : %f\n", net->v_thresh));
  net->printf(strx1("> Value Ratio                     : %f\n", net->v_ratio));
  net->printf(strx1("> Window Width                    : %f\n", net->width));
  if(1==0){
    int i, j;
    net->printf("w=\n");
    
    for(i=0;i<net->n_cells;i++){
      net->printf(strx1("%d ",i));
      for(j=0;j<net->k;j++){
	net->printf(strx1("%e ",net->cell[i].w[j]));
      }
      net->printf("\n");
    }
  }
    
  return;
}

/*====================================================================*
 * sort_weights
 * x_pの中のベクトルを入力xに近い順にソートする
 *====================================================================*/
/*
  int sort_weights(FLOAT *x, NET *net, int *s, FLOAT *v,
  int nn, int n_fires, int mode, int n0) {
  FLOAT temp;
  //FLOAT vmin;
  int ivmin = 0;
  int i, j;
  
  if (mode == 0) { // sの初期化，v[],s[0]の計算
  n_fires = 0;
  for (i=0; i<net->n_cells; i++) s[i] = i;
  // 入力xとの距離を求める
  for(i=0;i<net->n_cells;i++){
  v[i] = 0.0;
  if (net->cell[i].v > 0) {
  n_fires++;
  for (j=0;j<net->k;j++) v[i] += (FLOAT)square(net->cell[i].w[j]-x[j]);
  if(v[i] < v[ivmin]) ivmin = i;
  }
  else continue;
  }
  s[ivmin] = 0;
  s[0] = ivmin;
  }
  else {
  for (i=n0; i<nn; i++) {
  for (j=i+1; j<n_fires; j++) {
  if(v[(int)s[i]] > v[(int)s[j]]) {
  temp = s[i];
  s[i]= s[j];
  s[j]= temp;
  }
  }
  }
  }
  return(n_fires);
  }
*/

/*====================================================================*
 * init_memory
 * 配列s[]，v[]を初期化する
 *
 * s[i] ：荷重ベクトルcell[i].wのインデックスを入力ベクトルxに
 *        近い順に並たものを持つ配列
 * v[i] ：荷重ベクトルcell[i].wと入力ベクトルxとの距離を持つ配列
 *====================================================================*/
inline void init_sort_weights(NET *net, int *s, FLOAT *v, FLOAT *x) {
  FLOAT v_min = 1e+12;
  int i_v_min = 0;
  int n_cells = net->n_cells, n_channels = net->k;
  int n_fires = 0;
  int i, k;
  
  // 初期化
  for (i=0; i<n_cells; i++) s[i] = i;
  n_fires = 0;
  
  for (i=0; i<n_cells; i++) {
    v[i] = 0.0;
    if (net->cell[i].v > 0) {
      //    if(1==1){
      ++n_fires;
      for (k=0; k<n_channels; k++) v[i] += square(net->cell[i].w[k]-x[k]);
      if (v[i] < v_min) {
	i_v_min = i;
	v_min = v[i];
      }
    }
    else v[i] = 1e+12;
  }
  net->n_fires = n_fires;
  // 再近隣セルのインデックス
  s[i_v_min] = 0;
  s[0] = i_v_min;

  return;
}

/*====================================================================*
 * sort_weights
 * 配列s[]を並び替える
 *
 * s[i] ：荷重ベクトルcell[i].wのインデックスを入力ベクトルxに
 *        近い順に並たものを持つ配列
 * v[i] ：荷重ベクトルcell[i].wと入力ベクトルxとの距離を持つ配列
 *====================================================================*/
inline void sort_weights(NET *net, int *s, FLOAT *v, int start, int end) {
  //  char name[32] = "sort_weights";
  int n_fires = net->n_fires;
  int i, j, temp;
  
  if (n_fires < end) return;
  // ソート
  for (i=start; i<end; i++) {
    for (j=i+1; j<n_fires; j++) {
      if (v[(int)s[i]] > v[(int)s[j]]) {
	temp = s[i];
	s[i] = s[j];
	s[j] = temp;
      }
    }
  }
  return;
}
inline FLOAT cosw0xw1x(FLOAT *w1, FLOAT *w0, FLOAT *x,int k) {
  FLOAT w01=0,w0n=0,w1n=0,d0,d1;
  int i;
  for(i=0;i<k;i++){
    w01 += ((d0=w0[i]-x[i])*(d1=w1[i]-x[i]));
    w0n += d0*d0;
    w1n += d1*d1;
  }
  if(w0n<1e-10) return(2);
//  if (w01/sqrt(w0n)/sqrt(w1n)<-0.1) {
//    printf("ct=%f=%e/%e/%e\n",w01/sqrt(w0n)/sqrt(w1n),w01,sqrt(w0n),sqrt(w1n));
//  }
  return(w01/sqrt(w0n)/sqrt(w1n));
}
  
 
/*====================================================================*
 * calc_output
 * 区分的最小自乗法
 *====================================================================*/
FLOAT calc_output_c(NET *net, FLOAT *x, 
		  FLOAT *yr,//normalized && !digitized for recursive call
		  FLOAT *y, //normalized && digitized  
		  FLOAT *Y  //!normalized && digitized for error evaluate
		  ){
  FLOAT *v;
  int *s;
  //  FLOAT y;
  FLOAT ymin,ymax;
  int n_channels=net->k,n_fires = 0.0;
  int i, j,i0=-1,i1=-1;

  s = (int *)malloc(sizeof(int)*net->n_cells);
  v = (FLOAT *)malloc(sizeof(FLOAT)*net->n_cells);
  if (s == NULL) { printf("memory allocation error #4.\n"); exit(-1); }
  if (v == NULL) { printf("memory allocation error #5.\n"); exit(-1); }
  //n_fires = sort_weights(x, net, s, v, n_fires, n_fires, 0, 0);

  {
    init_sort_weights(net, s, v, x);
    n_fires = net->n_fires;
    for (i=0; i<n_fires;) {//?? vが正のユニットを探す
      if(net->cell[s[i]].v > ZERO){
	if(i0<0){
	  i0=i;
	  if(distance2(net->cell[s[i0]].w,x,net->k)<1e-20){//w0=x
	    i1=-2;
	    break;
	  }
	}
	//	else if(cosw0xw1x(net->cell[s[i]].w,net->cell[s[i0]].w,x,net->k)<-0.5) {//w1 is oppositeside
	else if(cosw0xw1x(net->cell[s[i]].w,net->cell[s[i0]].w,x,net->k)<net->cost) {//w1 is oppositeside
	//	else if(cosw0xw1x(net->cell[s[i]].w,net->cell[s[i0]].w,x,net->k)<-0.5) {//w1 is oppositeside
	//	else if(cosw0xw1x(net->cell[s[i]].w,net->cell[s[i0]].w,x,net->k)<-0.2) {//w1 is oppositeside
	//	else if(cosw0xw1x(net->cell[s[i]].w,net->cell[s[i0]].w,x,net->k)<-0.1) {//w1 is oppositeside
	  i1=i;
	  break;//after normalized v 030628ekuro???
	}
      }
      if(++i >= n_fires) break;
      sort_weights(net, s, v, i, i+1);
    }
    if(i0<0) i0=0;
    net->c=s[i0];//for NIPS04
    for(*y = 0.0,j=0; j<=n_channels; j++) *y += net->cell[s[i0]].am.M[0][j]*x[j];

    if(i1<0) net->dc=i1;//no opposite weight vector
    else {
      double yM0x,yM1x,yM0w0,yM1w0;
      //      double e1=1e-4,e2=1e-4,dy0,dyx,t;
      double e1=1e-10,e2=1e-4,dy0,dyx,t;
      yM0x=yM1x=yM0w0=yM1w0=0;
      yM0x=*y;
      for(j=0; j<=n_channels; j++) yM1x  += net->cell[s[i1]].am.M[0][j]*x[j];
      for(j=0; j< n_channels; j++) yM0w0 += net->cell[s[i0]].am.M[0][j]*net->cell[s[i0]].w[j];
      yM0w0 += net->cell[s[i0]].am.M[0][n_channels];
      for(j=0; j< n_channels; j++) yM1w0 += net->cell[s[i1]].am.M[0][j]*net->cell[s[i0]].w[j];
      yM1w0 += net->cell[s[i1]].am.M[0][n_channels];
      dy0=yM0w0-yM1w0;
      dyx=yM0x-yM1x;
      if(fabs(dy0-dyx)<e1){
	if(fabs(dyx)<e2){//continuous
	  //	  *y=yM0x;
	  net->dc=0;
	}
	else{//discontinuous
	  //	  *y=yM0x;
	  net->dc=1;
	}
      }
      else{
	t=dy0/(dy0-dyx);
	if(t<0) {//discontinuous
	  //	  *y=yM0x;
	  //	  *y=yM1x;//??
	  net->dc=2;
	  //	  if(t<-20) net->dc=5;
	  if(fabs(dyx)>1)  net->dc=5;
	}
	else if(t>1){//do not care
	  //	  *y=yM0x;
	  //	  *y=yM1x;//??
	  net->dc=3;
	  //	  if(t>20) net->dc=6;
	  if(fabs(dyx)>1)  net->dc=6;
	}
	else {//cont
	  *y=yM1x;
	  net->c=s[i1];
	  net->dc=4;
	}
      }
    }
  }

  ymin=net->ymin; ymax=net->ymax;
  if(*y<ymin) *y=ymin;
  if(*y>ymax) *y=ymax;//check
  //resolution 
//  if(net->r3>0){//for nonlinear parameter;old version modified at 060211
//    *y=net->ywidth*pow((*y-ymin)/net->ywidth,1./net->r3);
//  }
  *yr=*y;
  if(net->r1>0){
    if(*yr<net->r[0]) *yr=net->r[0];//
    else if(*yr>net->r[net->nr]) *yr=net->r[net->nr];
  }
  if(net->r1>0){
    int ii;
    //    int ii=(int)((*y-net->ymin)/net->r12+0.5);
    ii=floor((*y-ymin)/net->r12+0.5);//4sha5nyu
    if(ii>net->nr){
      net->printf("******************* set net->nr larger, next time *********\n");
      ii=net->nr;
    }
    else if(ii<0){
      net->printf("******************* set net->ymin smaller, next time *********\n");
      ii=0;
    }
    *y=net->r[ii];
    *Y=net->R[ii];
  }
  else{
    *Y=moverange(*y,ymin,ymax,net->ymin0,net->ymax0);
  }
  free(v);
  free(s);
  return(*yr);
}

FLOAT calc_output(NET *net, FLOAT *x, 
		  FLOAT *yr,//normalized && !digitized for recursive call
		  FLOAT *y, //normalized && digitized  
		  FLOAT *Y  //!normalized && digitized for error evaluate
		  ) {
  FLOAT *v;
  int *s;
  //  FLOAT y;
  FLOAT ymin,ymax;
  int n_channels=net->k,n_fires = 0.0;
  int i, j;

  s = (int *)malloc(sizeof(int)*net->n_cells);
  v = (FLOAT *)malloc(sizeof(FLOAT)*net->n_cells);
  if (s == NULL) { printf("memory allocation error #4.\n"); exit(-1); }
  if (v == NULL) { printf("memory allocation error #5.\n"); exit(-1); }
  //n_fires = sort_weights(x, net, s, v, n_fires, n_fires, 0, 0);
  init_sort_weights(net, s, v, x);
  n_fires = net->n_fires;

  for (i=0; i<n_fires;) {//??
    //if(net->cell[s[i]].v > net->v_thresh) break;
    //if(net->cell[s[i]].v >= net->k) break;//kuro???
    //if(net->cell[s[i]].v > 0.1) break;//after normalized v 030628ekuro???
    if(net->cell[s[i]].v > ZERO) break;//after normalized v 030628ekuro???
    if(++i >= n_fires) { i = 0; break; }
    sort_weights(net, s, v, i, i+1);
  }

  *y = 0.0;

  net->c=s[i];//for NIPS04

  for(j=0; j<=n_channels; j++) {
    *y += net->cell[(int)s[i]].am.M[0][j]*x[j];
  }
  //  ymin=net->ymin-0.1*net->ywidth;  ymax=net->ymax+0.1*net->ywidth;
  //  if(net->l_mode==BATCH_MODE){ ymin=net->ymin; ymax=net->ymax;} else {ymin=-1e30;ymax=1e30;}
//  {  // for debug
////////    printf("\nM=");for(j=0; j<=n_channels; j++) printf("%e ",net->cell[(int)s[i]].am.M[0][j]);
////////    printf("\nx=");for(j=0; j<=n_channels; j++) printf("%e ",x[j]);
////////    printf("y=%e",*y);
//////	printf("\n");
//	  for(j=0; j<=n_channels; j++) printf("%e ",x[j]);printf("%e",*y);
//  }
  ymin=net->ymin; ymax=net->ymax;
  if(*y<ymin) *y=ymin;
  if(*y>ymax) *y=ymax;//check
  //resolution 
  *yr=*y;
  if(net->r1>0){
    if(*yr<net->r[0]) *yr=net->r[0];//
    else if(*yr>net->r[net->nr]) *yr=net->r[net->nr];
  }
  if(net->r1>0){
    //    int ii=(int)((*y-net->ymin)/net->r12+0.5);
    int ii=floor((*y-net->ymin)/net->r12+0.5);//4sha5nyu
    if(ii>net->nr){
      printf("******************* set net->nr larger, next time *********\n");
      ii=net->nr;
    }
    else if(ii<0){
      printf("******************* set net->ymin smaller, next time *********\n");
      ii=0;
    }
    *y=net->r[ii];
    *Y=net->R[ii];
  }
  else{
    *Y=moverange(*y,net->ymin,net->ymax,net->ymin0,net->ymax0);
  }
//    {//check
//	    FLOAT yy;int j;
//	    for(yy=0,j=0;j<=net->k;j++) yy += net->cell[0].am.M[0][j]*x[j];
//	    if(fabs(*Y-yy)>1e-4){
//	      printf("check %+e=%+e\n",*Y,yy);
//	    }
//	    *Y=yy;
//    }

  //for data precision is "int".8bit
  free(v);
  free(s);
  return(*yr);
}
//FLOAT calc_output_real(NET *net, FLOAT *x) {
//  FLOAT *v;
//  int *s;
//  FLOAT y;
//  FLOAT ymin,ymax;
//  int n_channels=net->k,n_fires = 0.0;
//  int i, j;
//
//  s = (int *)malloc(sizeof(int)*net->n_cells);
//  v = (FLOAT *)malloc(sizeof(FLOAT)*net->n_cells);
//  if (s == NULL) { printf("memory allocation error #4.\n"); exit(-1); }
//  if (v == NULL) { printf("memory allocation error #5.\n"); exit(-1); }
//  //n_fires = sort_weights(x, net, s, v, n_fires, n_fires, 0, 0);
//  init_sort_weights(net, s, v, x);
//  n_fires = net->n_fires;
//
//  for (i=0; i<n_fires;) {//??
//    //if(net->cell[s[i]].v > net->v_thresh) break;
//    //if(net->cell[s[i]].v >= net->k) break;//kuro???
//    //if(net->cell[s[i]].v > 0.1) break;//after normalized v 030628ekuro???
//    if(net->cell[s[i]].v > ZERO) break;//after normalized v 030628ekuro???
//    if(++i >= n_fires) { i = 0; break; }
//    sort_weights(net, s, v, i, i+1);
//  }
//
//  y = 0.0;
//  for(j=0; j<=n_channels; j++) {
//    y += net->cell[(int)s[i]].am.M[0][j]*x[j];
//  }
//    ymin=net->ymin-0.1*net->ywidth;  ymax=net->ymax+0.1*net->ywidth;
//  //  ymin=net->ymin;  ymax=net->ymax;
//  if(y<ymin) y=ymin;
//  if(y>ymax) y=ymax;
//
//  if(net->r1>0){
//    if(y<0) y=0; //if(y>255) y=255;
//  }
//  //for data precision is "int".8bit
//  free(v);
//  free(s);
//  return(y);
//}
//
/*====================================================================*
 * calc_value
 *
 *====================================================================*/
inline void calc_value(NET *net) {//online
  int n_cells = net->n_cells;
  //  int n_fires = net->n_fires;
  FLOAT  sigma2_hat, alpha_min, alpha_bar;
  //  FLOAT v_total, alpha_hat, alpha_var;
  int i_alpha_min, alpha_NG;
  //  AM am;
  int i;
  //  FLOAT v_max=0,v_max3;
  //  int n_v_max3;
  //  int ivmax;//for debug
  FLOAT SvMin=1e30;
  VORONOI *V=net->V;//?????
  int n_alphai=0;
#define KURO
  //kuro:debug
  for (i=0; i<n_cells; i++) {
    if(net->cell[i].v < net->v_thresh) continue;
    //    if(net->cell[i].v*V->vmax < 2*n_channels) continue;//050825
    if(SvMin>(net->cell[i].S/net->cell[i].v)){
      SvMin=(net->cell[i].S/net->cell[i].v);
      n_alphai++;
    }
  }
  SvMin /= V->vmax;//040122 corrected
  //  printf("SvMin=%e\n",SvMin);
  net->sigma2_hat = sigma2_hat=SvMin;
  // 価値alphaiの平均値を求める
  alpha_bar = 0;
  for (i=0; i<n_cells; i++) {
        net->cell[i].alpha = 
	  (net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax))/n_alphai;//040122 corrected
	if(net->cell[i].v < net->v_thresh) continue;
	//	if(net->cell[i].v*V->vmax < 2*n_channels) continue;//050825
	alpha_bar += net->cell[i].alpha;
  }

  alpha_bar /= (FLOAT) n_alphai;
  if (alpha_bar< 1e-20) alpha_NG = 1;
  else alpha_NG = 0;//
  net->alpha_bar = alpha_bar;
  net->alpha_NG = alpha_NG;
  //
  net->printf(strx3("alpha_bar%3.2e,sigma2hat%3.2env%d\n",alpha_bar,sigma2_hat,n_alphai));
  // 価値alpha_iを求める
  i_alpha_min = 0;
  //alpha_min = net->cell[0].alpha - alpha_bar;
  alpha_min = net->cell[0].alpha;
  for (i=0; i<n_cells; i++) {
    if (net->cell[i].alpha < alpha_min) {
      i_alpha_min = i;
      alpha_min = net->cell[i].alpha;
    }
  }
  net->l = i_alpha_min;
  net->alpha_min = alpha_min;

  return;
}
//int calc_alpha_051201(NET *net, FLOAT **x, FLOAT *y, int n_train){//batch
//int calc_alpha_051201(NET *net, FLOAT **x, FLOAT *y, int n_train){//batch
//  int n_cells = net->n_cells;
//  FLOAT  sigma2_hat, alpha_bar;
//  int alpha_NG;
//  FLOAT y_m;
//  int i,j,k,t;
//  int n_channels=net->k;
//  FLOAT SvMin;
//  VORONOI *V=net->V;
//  FLOAT S=0.0,NDS;
//  double n_alphai=1e-20;
//  FLOAT alpha_sum=1e-20,aa,H=0;
//  for (i=0; i<n_cells; i++) {
//    net->cell[i].S = 0.0;
//    for (j=0; j<V->i2v0[i]; j++) {
//  //for (j=0; j<V->i2v[i]; j++) {
//      t=V->ij2t[i][j];
//      y_m=0;
//      for (k=0; k<=(n_channels); k++) y_m += net->cell[i].am.M[0][k]*x[t][k];
//      net->cell[i].S += (FLOAT)square(y[t]-y_m);
//    }
//#if GSL >= 1 
//    if(net->cell[i].S>1e5) net->cell[i].S=1e5;//???for infinite err of GSL //040211
//#endif
//    S += net->cell[i].S;
//  }
//#ifdef LAPTIME
//  printf("laptime:calc_alpha1=%5.3f\n",mytimer_lap());
//#endif
//
//  net->S[(int)GlobalTime] = S;
//  NDS=S/net->S[(int)GlobalTime-1]-1.;
//
//  printf(" %+f=NDS=S/S(t-1)-1=%e/%e-1\n",
//	 NDS, S,net->S[(int)GlobalTime-1]);
//  //  printf("laptime:calc_alpha2=%5.3f\n",mytimer_lap());
//  //  if(GlobalTime<ReinitTime+2) return(0);
//  //  if(GlobalTime<ReinitTime+GlobalTimeMax/10) return(0);
//  //  if(GlobalTime<ReinitTime+10) return(0);
//  //  if(GlobalTime<ReinitTime+10) return(0);
//  //  if(GlobalTime<ReinitTime+5) return(0);
//  //  if(GlobalTime<ReinitTime+10) return(0);//省く040123
//  //  if(GlobalTime<ReinitTime+1) return(0);//省く040123
//  //  if(NDS<0.001) return(0);//下へ移動040123
//  //  if((net->S[(int)GlobalTime] < net->S[(int)GlobalTime-1])&&(GlobalTime>GlobalTimeMax/5)) return(0);
//#undef actvthresh
//#define actvthresh
//  //kuro:debug
//  SvMin=1e30;
//  for (i=0; i<n_cells; i++) {
//#ifdef actvthresh
//    if(net->cell[i].v < net->v_thresh) continue;
//#else
//    if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//    n_alphai++;
//    if(SvMin>(net->cell[i].S/net->cell[i].v)){
//      SvMin=(net->cell[i].S/net->cell[i].v);
//    }
//  }
//  //  printf("sigma2hat=%e=(%e)/(%d)\n",SvMin/V->vmax,SvMin,V->vmax);
//  SvMin /= V->vmax;//040122 corrected
//  net->sigma2_hat = sigma2_hat=SvMin;
//  // 価値alpha_iの平均値を求める
//  for (i=0; i<n_cells; i++) {
//        net->cell[i].alpha = 
//	  (net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax))/n_alphai;//040122 corrected
//#ifdef actvthresh
//	if(net->cell[i].v < net->v_thresh) continue;
//#else
//	if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//	//	if(net->cell[i].alpha<0) continue;//050922
//	if(net->cell[i].alpha<0) net->cell[i].alpha*=-1;//050922
//	alpha_sum += net->cell[i].alpha;
//	//	printf("!!!!!!!!!!!alpha[%d]=%e\n",i,net->cell[i].alpha);
//  }
//  {//calc entropy
//    for(i=0; i<n_cells; i++){
//#ifdef actvthresh
//      if(net->cell[i].v < net->v_thresh) continue;
//#else 
//      if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//      aa = net->cell[i].alpha/alpha_sum;
//      //      aa =fabs(aa);
//      //      if(aa>1) aa=1;
//      //      if(aa<=0) continue;
//      H -= aa*(FLOAT)log(aa);
//    }
//    //net->nentropy = H/(FLOAT)logl(n_alphai);
//    if((int)n_alphai>=2) net->nentropy = H/(FLOAT)log(n_alphai); else net->nentropy =1;
//  }
////  if(GlobalTime>=23){
////    printf("check \n");
////  }
////  if((long int)GlobalTime ==44){
////    printf("check here");
////  }
//  alpha_bar = alpha_sum/(FLOAT) n_alphai;
//  if (alpha_bar< 1e-20) alpha_NG = 1;
//  else alpha_NG = 0;//
//  net->alpha_bar = alpha_bar;
//  net->alpha_NG = alpha_NG;
//  //
//  printf(" sigma2hat%3.2e#alphai%3.0f,<alphai>%3.2e,NormalizedEntropy=%f\n",
//	 sigma2_hat,n_alphai,alpha_bar,net->nentropy);
//  //debug
//  //  if(net->nentropy<0.95) return(1);//1 for reinit
//  if(NDS<0.0) return(0);//  if(NDS<0.001) return(0);
//  if(net->nentropy<net->nentropy_thresh) return(1);//1 for reinit original
//  //    if(net->nentropy<1-(1-net->nentropy_thresh)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
//  //    if(net->nentropy<0.5+(net->nentropy_thresh-0.5)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
//    else return(0);
//  //  if(alpha_NG==0) return(1);  return(0);
//}
int calc_alpha(NET *net, FLOAT **x, FLOAT *y, int n_train){//batch
  int n_cells = net->n_cells;
  FLOAT  sigma2_hat, alpha_bar;
  int alpha_NG;
  FLOAT y_m;
  int i,j,k,t;
  int n_channels=net->k;
  FLOAT SvMin;
  VORONOI *V=net->V;
  FLOAT S=0.0,NDS;
  double n_alphai=1e-20;
  FLOAT alpha_sum=1e-20,aa,H=0;
  for (i=0; i<n_cells; i++) {
    net->cell[i].S = 0.0;
    for (j=0; j<V->i2v0[i]; j++) {
      //for (j=0; j<V->i2v[i]; j++) {
      t=V->ij2t[i][j];
      y_m=0;
      for (k=0; k<=(n_channels); k++) y_m += net->cell[i].am.M[0][k]*x[t][k];
      net->cell[i].S += (FLOAT)square(y[t]-y_m);
      //      printf("%d)S%e,e%e tY%e gY%e\n",i,net->cell[i].S,y[t]-y_m,y[t],y_m);//check060517
    }
#if GSL >= 1 
    if(net->cell[i].S>1e5) net->cell[i].S=1e5;//???for infinite err of GSL //040211
#endif
    S += net->cell[i].S;
    //    printf("S=%e\n",S);//check060517
  }
#ifdef LAPTIME
  printf("laptime:calc_alpha1=%5.3f\n",mytimer_lap());
#endif

  net->S[(int)GlobalTime] = S;
  NDS=S/net->S[(int)GlobalTime-1]-1.;

//20191003  printf(" %+f=NDS=S/S(t-1)-1=%e/%e-1\n",
//20191003	 NDS, S,net->S[(int)GlobalTime-1]);
  //  printf("laptime:calc_alpha2=%5.3f\n",mytimer_lap());
  //  if(GlobalTime<ReinitTime+2) return(0);
  //  if(GlobalTime<ReinitTime+GlobalTimeMax/10) return(0);
  //  if(GlobalTime<ReinitTime+10) return(0);
  //  if(GlobalTime<ReinitTime+10) return(0);
  //  if(GlobalTime<ReinitTime+5) return(0);
  //  if(GlobalTime<ReinitTime+10) return(0);//省く040123
  //  if(GlobalTime<ReinitTime+1) return(0);//省く040123
  //  if(NDS<0.001) return(0);//下へ移動040123
  //  if((net->S[(int)GlobalTime] < net->S[(int)GlobalTime-1])&&(GlobalTime>GlobalTimeMax/5)) return(0);
#undef actvthresh
#define actvthresh
  //kuro:debug
  SvMin=1e30;
  for (i=0; i<n_cells; i++) {
#ifdef actvthresh
    if(net->cell[i].v < net->v_thresh) continue;
#else
    if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
#endif
    //        if(net->cell[i].alpha<ZERO) continue;//051207-L750
    n_alphai++;
    if(SvMin>(net->cell[i].S/net->cell[i].v)){
      SvMin=(net->cell[i].S/net->cell[i].v);
    }
  }
  //  printf("sigma2hat=%e=(%e)/(%d)\n",SvMin/V->vmax,SvMin,V->vmax);
  SvMin /= V->vmax;//040122 corrected
  //   printf("S[%d]=%e\n",i,net->cell[i].S);
  net->sigma2_hat = sigma2_hat=SvMin;
  // 価値alpha_iの平均値を求める
  for (i=0; i<n_cells; i++) {
    net->cell[i].alpha=(net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax))/n_alphai;//040122 corrected
#ifdef actvthresh
	if(net->cell[i].v < net->v_thresh) continue;
#else
	if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
#endif
	//	if(net->cell[i].alpha<0) continue;//050922
	//	if(net->cell[i].alpha<0) net->cell[i].alpha*=-1;//050922
	if(net->cell[i].alpha<ZERO) continue;//051207
	alpha_sum += net->cell[i].alpha;
	//	printf("!!!!!!!!!!!alpha[%d]=%e\n",i,net->cell[i].alpha);
  }
  {//calc entropy
    for(i=0; i<n_cells; i++){
#ifdef actvthresh
      if(net->cell[i].v < net->v_thresh) continue;
#else 
      if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
#endif
      if(net->cell[i].alpha<ZERO) continue;//051207
      aa = net->cell[i].alpha/alpha_sum;
      //      aa =fabs(aa);
      //      if(aa>1) aa=1;
      //      if(aa<=0) continue;
      H -= aa*(FLOAT)log(aa);
    }
    //net->nentropy = H/(FLOAT)logl(n_alphai);
    if((int)n_alphai>=2) net->nentropy = H/(FLOAT)log(n_alphai); else net->nentropy =1;
  }
//  if(GlobalTime>=23){
//    printf("check \n");
//  }
//  if((long int)GlobalTime ==44){
//    printf("check here");
//  }
  alpha_bar = alpha_sum/(FLOAT) n_alphai;
  if (alpha_bar< 1e-20) alpha_NG = 1;
  else alpha_NG = 0;//
  net->alpha_bar = alpha_bar;
  net->alpha_NG = alpha_NG;
  //
//20191003  printf(" sigma2hat%3.2e#alphai%3.0f,<alphai>%3.2e,NormalizedEntropy=%f\n",
//20191003	 sigma2_hat,n_alphai,alpha_bar,net->nentropy);
  //debug
  //  if(net->nentropy<0.95) return(1);//1 for reinit
  double NDS_thresh=0.0;
  int reinit_flag1=NDS>=NDS_thresh;
  int reinit_flag2=net->nentropy<net->nentropy_thresh;
  int reinit_flag=reinit_flag1 && reinit_flag2;
  char *NDSmes="";
  char *RImes="False";
  if(NDS>1.0) NDSmes="(NG!)";
  if(reinit_flag) RImes="True";
  net->printf(strx10("%d #check S=%g NDS%s%g<%g:%d entropy%g>=%g:%d-> Reinit=%s\n",net->GlobalTime,S,NDSmes,NDS,NDS_thresh,NDS<NDS_thresh,net->nentropy,net->nentropy_thresh,net->nentropy>=net->nentropy_thresh,RImes));
  return(reinit_flag);
  //  if(NDS<NDS_thresh) return(0);//  if(NDS<0.001) return(0);
  //  if(net->nentropy<net->nentropy_thresh) return(1);//1 for reinit original
  //    if(net->nentropy<1-(1-net->nentropy_thresh)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
  //    if(net->nentropy<0.5+(net->nentropy_thresh-0.5)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
    //    else return(0);
  //  if(alpha_NG==0) return(1);  return(0);
}

//int calc_alpha051208(NET *net, FLOAT **x, FLOAT *y, int n_train){//batch
//  int n_cells = net->n_cells;
//  FLOAT  sigma2_hat, alpha_bar;
//  int alpha_NG;
//  FLOAT y_m;
//  int i,j,k,t;
//  int n_channels=net->k;
//  FLOAT SvMin;
//  VORONOI *V=net->V;
//  FLOAT S=0.0,NDS;
//  double n_alphai;
//  FLOAT alpha_sum,aa,H=0;
//  for (i=0; i<n_cells; i++) {
//    net->cell[i].S = 0.0;
//    for (j=0; j<V->i2v0[i]; j++) {
//  //for (j=0; j<V->i2v[i]; j++) {
//      t=V->ij2t[i][j];
//      y_m=0;
//      for (k=0; k<=(n_channels); k++) y_m += net->cell[i].am.M[0][k]*x[t][k];
//      net->cell[i].S += (FLOAT)square(y[t]-y_m);
//    }
//#if GSL >= 1 
//    if(net->cell[i].S>1e5) net->cell[i].S=1e5;//???for infinite err of GSL //040211
//#endif
//    S += net->cell[i].S;
//  }
//#ifdef LAPTIME
//  printf("laptime:calc_alpha1=%5.3f\n",mytimer_lap());
//#endif
//
//  net->S[(int)GlobalTime] = S;
//  NDS=S/net->S[(int)GlobalTime-1]-1.;
//
//  printf(" %+f=NDS=S/S(t-1)-1=%e/%e-1\n",
//	 NDS, S,net->S[(int)GlobalTime-1]);
//#undef actvthresh
//#define actvthresh
//  //kuro:debug
//  SvMin=1e30;
//  for (i=0; i<n_cells; i++) {
//#ifdef actvthresh
//    if(net->cell[i].v < net->v_thresh) continue;
//#else
//    if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//    //        if(net->cell[i].alpha<ZERO) continue;//051207-L750
//    if(SvMin>(net->cell[i].S/net->cell[i].v)){
//      SvMin=(net->cell[i].S/net->cell[i].v);
//    }
//  }
//  //  printf("sigma2hat=%e=(%e)/(%d)\n",SvMin/V->vmax,SvMin,V->vmax);
//  SvMin /= V->vmax;//040122 corrected
//  //   printf("S[%d]=%e\n",i,net->cell[i].S);
//  net->sigma2_hat = sigma2_hat=SvMin;
//  // 価値alpha_iの平均値を求める
//  {
//    n_alphai=0;
//    for (i=0; i<n_cells; i++) {
//      net->cell[i].alpha=(net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax));//040122 corrected
//#ifdef actvthresh
//      if(net->cell[i].v < net->v_thresh) continue;
//#else
//      if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//      if(net->cell[i].alpha>ZERO) n_alphai++;//051207
//    }
//  }
//  if((int)n_alphai>2) {//calc entropy
//    for (i=0; i<n_cells; i++) {
//      //net->cell[i].alpha=(net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax))/n_alphai;//040122 corrected
//      //    net->cell[i].alpha=(net->cell[i].S-sigma2_hat*(net->cell[i].v*V->vmax));//040122 corrected
//      net->cell[i].alpha/=n_alphai;//040122 corrected
//#ifdef actvthresh
//      if(net->cell[i].v < net->v_thresh) continue;
//#else
//      if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//      //	if(net->cell[i].alpha<0) continue;//050922
//      //	if(net->cell[i].alpha<0) net->cell[i].alpha*=-1;//050922
//      if(net->cell[i].alpha>ZERO) alpha_sum += net->cell[i].alpha;//051207
//      //	printf("!!!!!!!!!!!alpha[%d]=%e\n",i,net->cell[i].alpha);
//    }
//    for(i=0; i<n_cells; i++){
//#ifdef actvthresh
//      if(net->cell[i].v < net->v_thresh) continue;
//#else 
//      if(net->cell[i].v*V->vmax < n_channels+5) continue;//050825
//#endif
//      if(net->cell[i].alpha>ZERO){
//	aa = net->cell[i].alpha/alpha_sum;
//	H -= aa*(FLOAT)log(aa);
//      }
//    }
//    //net->nentropy = H/(FLOAT)logl(n_alphai);
//    net->nentropy = H/(FLOAT)log(n_alphai);
//    net->alpha_bar = alpha_bar = alpha_sum/(FLOAT) n_alphai;
//    if (alpha_bar< 1e-20) net->alpha_NG = alpha_NG = 1;
//    else  net->alpha_NG = alpha_NG = 0;//
//  }
//  else{
//    net->nentropy =1;
//    net->alpha_bar = alpha_bar=0;
//    net->alpha_NG = alpha_NG= 1;
//  }
//  //    if((int)n_alphai>=2) net->nentropy = H/(FLOAT)log(n_alphai); else net->nentropy =1;
//  //    if((int)n_alphai>=2) net->nentropy = H/(FLOAT)log(n_alphai); else net->nentropy =1;
//  //  if(GlobalTime>=23){
//  //    printf("check \n");
//  //  }
//  //  if((long int)GlobalTime ==44){
//  //    printf("check here");
//  //  }
//  //
//  printf(" sigma2hat%3.2e#alphai%3.0f,<alphai>%3.2e,NormalizedEntropy=%f\n",
//	 sigma2_hat,n_alphai,alpha_bar,net->nentropy);
//  //debug
//  //  if(net->nentropy<0.95) return(1);//1 for reinit
//  if(NDS<0.0) return(0);//  if(NDS<0.001) return(0);
//  if(net->nentropy<net->nentropy_thresh) return(1);//1 for reinit original
//  //    if(net->nentropy<1-(1-net->nentropy_thresh)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
//  //    if(net->nentropy<0.5+(net->nentropy_thresh-0.5)/(1.+GlobalTime/net->Tgamma)) return(1);//1 for reinit 050517
//    else return(0);
//  //  if(alpha_NG==0) return(1);  return(0);
//}

//inline void calc_value0(NET *net) {
//  int n_cells = net->n_cells;
//  int n_fires = net->n_fires;
//  FLOAT v_total, alpha_hat, sigma2_hat, alpha_min, alpha_bar, alpha_var;
//  int i_alpha_min, alpha_NG;
//  AM am;
//  int i;
//  FLOAT v_max=0,v_max3;
//  int n_v_max3;
//  //  int ivmax;//for debug
//  
//  // v_total（n_hat）, alpha_hat, sigma^2_hat を求める
//#define KURO
//  //kuro
//  //  for (i=0; i<n_cells; i++){
//  //    if(v_max<net->cell[i].v){
//  //	v_max=net->cell[i].v;
//  //	ivmax=i;
//  //    }
//  //  }
//  //    net->v_thresh=v_max3=v_max*2./3.;//intuition???
//  //  net->v_thresh=v_max3=v_max*1./3.;//intuition???
//  //  net->v_thresh=v_max3=v_max*0.5;//intuition???
//  //    net->v_thresh=v_max3=v_max*0.9999;//intuition???
//  {//kuro:debug
//    FLOAT SvMin=1e30;
//    for (i=0; i<n_cells; i++) {
//	if(SvMin>(net->cell[i].S/net->cell[i].v)){
//	  SvMin=(net->cell[i].S/net->cell[i].v);
//	}
//    }
//    printf("SvMin=%e\n",SvMin);
//  }
//  v_total=0;n_v_max3 = 0;
//  for (i=0; i<n_cells; i++){
//    v_total += net->cell[i].v;
//    if(net->cell[i].v>=net->v_thresh){
//	n_v_max3++;
//    }
//  }
//  
//  init_AM(&am, 2, 1);
//  am.x[0] = v_total;
//  //  am.x[0] = n_v_max3;
//  for (i=0; i<n_cells; i++) {
//    if(net->cell[i].v < net->v_thresh) continue;
//    am.x[1] = net->cell[i].v;
//    am.y[0] = net->cell[i].S;
//    calc_AM(am);
//  }
//  
//  alpha_hat = am.M[0][0];
//  sigma2_hat = am.M[0][1];
//  free_AM(&am);
//  //  if(alpha_hat<0) alpha_hat=0;
//  net->alpha_hat = alpha_hat;
//  net->sigma2_hat = sigma2_hat;
//
//  // 価値alpha_iの平均値を求める
//  alpha_bar = 0;
//  for (i=0; i<n_cells; i++) {
//	  net->cell[i].alpha = (net->cell[i].S-sigma2_hat*(net->cell[i].v))/v_total;
//    // net->cell[i].alpha = (net->cell[i].S-sigma2_hat*(net->cell[i].v))/n_v_max3;
//    if(net->cell[i].v < net->v_thresh) continue;
//    alpha_bar += net->cell[i].alpha;
//    //    printf("%d]alpha%3.2e S%3.2e v%3.2e\n",i,net->cell[i].alpha,net->cell[i].S,net->cell[i].v);
//  }
//
//  alpha_bar /= (FLOAT) n_v_max3;
//  if (alpha_bar< 1e-20) alpha_NG = 1;
//  else alpha_NG = 0;//
//  net->alpha_bar = alpha_bar;
//  net->alpha_NG = alpha_NG;
//  //
//  //  printf("alpha_bar%3.2e,alhat%3.2e,sigma2hat%3.2e,vmax%3.2e,ivmax%d,nv%d\n",alpha_bar,alpha_hat,sigma2_hat,v_max,ivmax,n_v_max3);
//  printf("alpha_bar%3.2e,alhat%3.2e,sigma2hat%3.2env%d\n",alpha_bar,alpha_hat,sigma2_hat,n_v_max3);
//  // 価値alpha_iを求める
//  i_alpha_min = 0;
//  //alpha_min = net->cell[0].alpha - alpha_bar;
//  alpha_min = net->cell[0].alpha;
//  for (i=0; i<n_cells; i++) {
//    if (net->cell[i].alpha < alpha_min) {
//	i_alpha_min = i;
//	alpha_min = net->cell[i].alpha;
//    }
//  }
//  net->l = i_alpha_min;
//  net->alpha_min = alpha_min;
//  
//  return;
//}

/*====================================================================*
 * modify_x
 *
 *====================================================================*/
void modify_x(NET *net, int m, FLOAT *x, int y) {
  FLOAT alpha;
  int n_channels = net->k;
  int k;

  /*
#ifdef OLD
  // v==0 なら net->cell[m].w = x(初期化)
  // v==inf なら net->cell[m].w は変化なし
  for(j=0; j<(n_channels+1); j++){
    net->cell[m].w[j] = alpha*net->cell[m].w[j]+(1.-alpha)*x[j];
  }
#endif
#if L_MODE == 5
  rate_w = 1.0;
#endif
  */
  alpha = rate_w/(net->cell[m].v+1.);
  
  if(alpha > 1.0) alpha = 1.0;
#if L_MODE==1
  if(alpha < -0.1) alpha = -0.1;
#else
  if(alpha < -0.5) alpha = -0.5;
#endif

  for(k=0; k<n_channels; k++){
    net->cell[m].w[k] = (1.0-alpha)*net->cell[m].w[k]+alpha*x[k];
  }

  return;
}

/*====================================================================*
 * modify_v
 *
 *====================================================================*/
void modify_v(NET *net, int m, FLOAT *x, FLOAT y) {
  FLOAT y_m;
  int n_channels = net->k;
  int i, j;

  //  if (net->cell[m].v >= net->v_thresh) {//発火回数が多いユニットの誤差でsigmaを計算するため?
    y_m = 0.0;
    for(j=0; j<=n_channels; j++) y_m += net->cell[m].am.M[0][j]*x[j];
    net->cell[m].S += (FLOAT)square(y_m-y);
    //  }
  for (i=0; i<net->n_cells; i++) {
    net->cell[i].S *= net->r_E;
    net->cell[i].v *= net->r_E;
  }
  net->cell[m].v++; // 価値を減衰
  return;
}

/*====================================================================*
 * modify_M
 *
 *====================================================================*/
void modify_M(NET *net, int m, FLOAT *x, FLOAT y) {
  int n_channels = net->k;
  int k;

  //recursive_pinv(net->cell[m].P,net->cell[m].Q,net->cell[m].M,x,y);
  //net->cell[m].am.x = x;
  for (k=0; k<=n_channels; k++) net->cell[m].am.x[k] = x[k];
  net->cell[m].am.y[0] = y;
  calc_AM(net->cell[m].am);
  //  printf("m=%d.",m);  show_weights(net);//check
  return;
}

/*====================================================================*
 * modify_cell
 *
 *====================================================================*/
void modify_cell(NET *net, int m, FLOAT *x, FLOAT y) {
  modify_x(net, m, x, y);
  modify_v(net, m, x, y);
  modify_M(net, m, x, y);
  return;
}

/*====================================================================*
 * in_window
 *
 *====================================================================*/
int in_window(FLOAT *x1, FLOAT *x2, FLOAT *x, FLOAT width, int n_channels) {
  FLOAT l2 = 0.0, l3 = 0.0, l3l2, www;
  int ret;
  int i;

  for (i=0; i<n_channels; i++) {
    l2 += (FLOAT)square(x1[i]-x2[i]);
    l3 += (x[i]-x2[i])*(x1[i]-x2[i]);
  }
  if (l2 < 0.0) ret = 0;
  else {
    www = (1.-width)/2.;
    l3l2 = l3/l2;
    ret = (int)((www < l3l2) && (l3l2 < 1.-www));
  }
  return(ret);
}

#ifdef E_CHECK
/*====================================================================*
 * calc_entropy
 *
 *====================================================================*/
void calc_entropy(FLOAT *d, int N) {
  FLOAT H, sumd, dd;
  int i;

  sumd = 0.0;
  H = 0.0 ;
  for(i=0; i<N; i++) sumd+=d[i];
  for(i=0; i<N; i++){
    dd = d[i]/sumd;
    H -= dd*(FLOAT)log(dd);
  }
  entropy = H/(FLOAT)log(N);
  return;
}

/*====================================================================*
 * check_entropy_of_alpha
 *
 *====================================================================*/
void check_entropy_of_alpha(NET *net) {
  FLOAT *alpha, alpha_bar, alpha_hat, sigma2_hat, vv;
  AM am;
  int i;

  alpha = (FLOAT *)malloc(sizeof(FLOAT)*net->n_cells);
  if (alpha == NULL) { printf("memory allocation error #6.\n"); exit(-1); }

  vv = 0;
  for (i=0; i<net->n_cells; i++) vv += net->cell[i].v;
  init_AM(&am,2,1);
  am.x[0] = vv;
  for (i=0; i<net->n_cells; i++) {
    am.x[1] = net->cell[i].v;
    am.y[0] = net->cell[i].S;
    calc_AM(am);
  }
  alpha_hat = am.M[0][0];
  sigma2_hat = am.M[0][1];
  free_AM(&am);

  alpha_bar = 0;
  for (i=0; i<net->n_cells; i++) {
    alpha[i] = (net->cell[i].S-sigma2_hat*net->cell[i].v)/vv; // alpha_iを一定
    // if(alpha[i] < 0.0) alpha[i] = 0.0;
  }
  calc_entropy(alpha, net->n_cells);
  free(alpha);
  return;
}
#endif

//#include "my_plinn_online.c"
/*====================================================================*
 * store_vector
 *
 *====================================================================*/
FLOAT TIME_FOR_REINIT=0;

int store_vector(NET *net, FLOAT *x, FLOAT y) 
{//online
  int i,j,nmax=0,modified=-1;
  //  FLOAT xmax,xmin;
  FLOAT *v;
  int *s;//,*ss;
  //  FLOAT errmin=100://,yi://,temp;
  //  int ierrmin=0,iclosest=0;
  int s0;
  AM am;
  FLOAT vv;
  FLOAT value_hat,sigma2_hat,value_bar,value_var,value_min;
  FLOAT *value;
  FLOAT err0,erri;
  int i_value_min=0;
  int alpha_NG=0;
  int n_cells =net->n_cells;
  int tau_E=net->tau_E;
  FLOAT vt=0;
  int   n_vt;
  //  FLOAT WindowWidth=net->width;/**/
  
//    if(GlobalTime>=52){
//	show_weights(net);
//	printf("check\n");
//    }
  v=(FLOAT *)malloc(sizeof(FLOAT)*n_cells);
  s=(int *)malloc(sizeof(int)*n_cells);
  //  nmax=sort_memory(x,net,s,v,nmax,nmax,0,0);
  init_sort_weights(net, s, v, x);
  nmax=net->n_fires;
  s0=(int)s[0];
  
  if(net->cell[s0].v<net->v_thresh){/*まだ十分記憶してない*/
    modify_cell(net,s0,x,y);
    modified=s0;
#ifdef DEBUG
    if(debug==2){
      printf("\nsv12(%d,%f)",modified,errmin);/*DEBUG*/
    }
#endif
  }
#if L_MODE!=-1
  /***L_MODE==1の場合，全セルをはじめに初期化しておかないと
      value[i](=alpha[i])の計算が不安定*/
  else if(nmax<n_cells){
    modify_cell(net,nmax,x,y);
    modified=nmax;
  }
#else
#endif
  else{/* if(net->cell[s0].v>=net->v_thresh){*/
    value=(FLOAT*)malloc(sizeof(FLOAT)*n_cells);
    for(vv=0,i=0;i<n_cells;i++) vv+=net->cell[i].v;

#if L_MODE==1 || L_MODE ==14
    if(1==1){//original
      init_AM(&am,2,1);/**/
      am.x[0]=vv;
      for(i=0;i<n_cells;i++){
	am.x[1]=net->cell[i].v;
	am.y[0]=net->cell[i].S;
	calc_AM(am);
      }
      value_hat=am.M[0][0];
      sigma2_hat=am.M[0][1];
      free_AM(&am);
    }
    else{//new
      vt=n_vt=0;
      for(i=0;i<n_cells;i++)if(vt<net->cell[i].v) vt=net->cell[i].v;
      sigma2_hat=1e30;
      vt/=2;
      for(i=0;i<n_cells;i++){
	if(net->cell[i].v<vt) continue;
	n_vt++;
	if(sigma2_hat>(net->cell[i].S/net->cell[i].v)) sigma2_hat=(net->cell[i].S/net->cell[i].v);
      }
    }
    //printf("sigma2hat=%e\n",sigma2_hat);
#else
    value_hat=sigma2_hat=0;
#endif
    
#if L_MODE==1 || L_MODE==14
    /*alpha_i*/
    for(value_bar=0,i=0;i<nmax;i++){
      value[i]=(net->cell[i].S-sigma2_hat*net->cell[i].v)/vv;/*alpha_i 一定*/
      if(value[i]<0) alpha_NG=1;
      value_bar+=value[i];
    }
#elif L_MODE==2
    for(value_bar=0,i=0;i<n_cells;i++){
      value[i]=net->cell[i].S/(net->cell[i].v+1);/*S_i/n_i一定*/
      value_bar+=value[i];
    }
#elif L_MODE==3 || L_MODE==4 || L_MODE==34
    for(value_bar=0,i=0;i<n_cells;i++){
      value[i]=net->cell[i].v;/*等発火(n_i)*/
      value_bar+=value[i];
    }
#endif
#ifdef CHECK
#if L_MODE==1 || L_MODE==14
    calc_entropy(value,n_cells);/**/
    //      /*check_entropy_of_alpha(net);/**/
#else
    if(dj==0){
      calc_entropy(value,n_cells);/**/
      ///*check_entropy_of_alpha(net);/**/
    }
#endif      
#endif      
    value_bar/=n_cells;//value_bar/=n_vt; //   
    value_var=0;
    i_value_min=0;
    value_min=value[0]-value_bar;
    for(i=0;i<n_cells;i++){
      //value[i]-=value_bar;
      net->cell[i].alpha=(value[i]-=value_bar);//ueno
      value_var+=value[i]*value[i];
      /*      if(value[i]<value_min&&net->cell[i].v>=net->v_thresh){*/
      if(value[i]<value_min){
	value_min=value[i];
	i_value_min=i;
      }
    }
    value_var=sqrt(value_var/n_cells);
    TIME_FOR_REINIT+=1;
#if L_MODE==1 
    if(TIME_FOR_REINIT>tau_E/4. && alpha_NG==0 && value[s0]>2.0*value_bar){/**/
       TIME_FOR_REINIT=0;
#elif L_MODE==2 ||L_MODE==3 
    if(value[s0]>3.0*value_var){/*&&net->cell[s0].v>=net->v_thresh*/
#elif L_MODE==4
	for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
    if(fabs(err0)>errlimit){/*&&net->cell[s0].v>=net->v_thresh*/
#elif L_MODE==34
    for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
    if(TIME_FOR_REINIT>tau_E/4.&&fabs(err0)>errlimit||value[s0]>2.0*value_bar){
	TIME_FOR_REINIT=0;
#elif L_MODE==14
    for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
    if(nmax<n_cells) {modify_cell(net,nmax,x,y);modified=nmax;}
    else if(fabs(err0)>errlimit||value[s0]>3.0*value_var){/**/
	///*    else if((fabs(err0)>errlimit||value[s0]>3.0*value_var)&&Entropy<0.985){/**/
#elif L_MODE==5
    if(nmax<n_cells) {modify_cell(net,nmax,x,y);modified=nmax;}
    if(1==0){
#endif 
      /*一定化すべき変数が一定でないとき*/
      if(nmax<n_cells){/* 記憶容量未満なら使ってないセルでそのまま記憶 */
	modify_cell(net,nmax,x,y);
	modified=nmax;/*??*/
      }
      else{
	/*最小価値のセルで記憶 
	  nmax>n_cells&&value[s0]>1.1*value_var&&net->cell[s0].v>=net->v_thresh*/
	CELL *c0,*c1;
	c0=&net->cell[s0];
	c1=&net->cell[i_value_min];
	c1->v=0;/*initialize i_value_min*/
	init_AMdata(&(net->cell[i_value_min].am));
	///*	if(_TIME<100000)	modify_cell(net,i_value_min,x,y);*/
	fprintf(stderr,"Reinit(%d,%e,%e)\n",i_value_min,value[s0],value_var);/*DEBUG*/
	modify_cell(net,i_value_min,x,y);
	//	printf("Reinit(%d,%3.1f)",i_value_min,value[s0]/value_var);/*DEBUG*/
	modified=i_value_min;
      }
    }
    else {/*uniform のとき vales[s0]<1.1*value_var&&net->cell[s0].v>=net->v_thresh*/
      FLOAT delta_wic;
      int si;
#if L_MODE==1 || L_MODE==2 || L_MODE==3
      for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];//check
      err0=square(err0);
#elif L_MODE==4 || L_MODE==34 || L_MODE==14 || L_MODE==5
      err0=square(err0);
#endif
      for(i=1;i<net->n_compare;i++){/*境界の調節*/
	sort_weights(net, s, v, i, i+1);//	nmax=sort_memory(x,net,s,v,i+1,nmax,1,i);
	si=s[i];
	if(net->cell[si].v<net->v_thresh) continue;
	//	WindowWidth=0.1;/**/
	//	WindowWidth=0.5;
	if(in_window(net->cell[si].w,net->cell[s0].w, x, net->width,net->k)>0){/*境界*/
	  for(delta_wic=0,j=0;j<=net->k;j++) 
	    delta_wic+=square(net->cell[s0].w[j]-net->cell[si].w[j]);
	  //	    if((delta_wic=sqrt(delta_wic))<ZERO) delta_wic=ZERO;
	  if((delta_wic=sqrt(delta_wic))<ZERO) break;/**/
	  for(erri=-y,j=0;j<=net->k;j++) erri += net->cell[s[i]].am.M[0][j]*x[j];
	  erri=square(erri);
	  ///*	  rate_w=-(err0-erri)/delta_wic/0.001;/**f_07: MSE=0.00229494 MODE3*/
	  rate_w=-(err0-erri)/delta_wic/0.01;/*f_07: MSE=0.00018038 MODE1
					      *f_07: MSE=0.00042767 MODE3****/
	  ///*	  rate_w=-(err0-erri)/delta_wic/0.1; /*f_07: MSE=0.00013 MODE1***
	  //					    *f_07: MSE=0.00050 MODE3*/
	  ///*	  rate_w=-(err0-erri)/delta_wic/1.0; /*f_07: MSE=0.00052056  MODE1*/
//	    if(GlobalTime>=53){//check
//	      printf("err=%e,%e\n",err0,erri);
//	    }

	  modify_cell(net,s0,x,y);/**/
	  modified=s0;
	  rate_w=1.0;
	  break;
	}
	else continue;
      }
    }
    if(modified<0){/*上で変更されていない時*/
      /*		  modify_x(net,s0,x,y);*/
      modify_v(net,s0,x,y);
      modify_M(net,s0,x,y);
#ifdef DEBUG
      if(debug==2){
	modified=s0;
	printf("\nsv24(%d,%f)",modified,errmin);/*DEBUG*/
      }
#endif
    }
    free(value);
  }
  free(v);
  free(s);
  return(0);
}

//int store_vector1(NET *net, FLOAT *x, FLOAT y) 
//{//online
//  int i,j,nmax=0,modified=-1;
//  //  FLOAT xmax,xmin;
//  FLOAT *v;
//  int *s;//,*ss;
//  //  FLOAT errmin=100://,yi://,temp;
//  //  int ierrmin=0,iclosest=0;
//  int s0;
//  AM am;
//  FLOAT vv;
//  FLOAT value_hat,sigma2_hat,value_bar,value_var,value_min;
//  FLOAT *value;
//  FLOAT err0,erri;
//  int i_value_min=0;
//  int alpha_NG=0;
//  int n_cells =net->n_cells;
//  int tau_E=net->tau_E;
//  FLOAT vt=0;
//  int   n_vt;
//  //  FLOAT WindowWidth=net->width;/**/
//  
////    if(GlobalTime>=52){
////	show_weights(net);
////	printf("check\n");
////    }
//  v=(FLOAT *)malloc(sizeof(FLOAT)*n_cells);
//  s=(int *)malloc(sizeof(int)*n_cells);
//  //  nmax=sort_memory(x,net,s,v,nmax,nmax,0,0);
//  init_sort_weights(net, s, v, x);
//  nmax=net->n_fires;
//  s0=(int)s[0];
//  
//  if(net->cell[s0].v<net->v_thresh){/*まだ十分記憶してない*/
//    modify_cell(net,s0,x,y);
//    modified=s0;
//#ifdef DEBUG
//    if(debug==2){
//	printf("\nsv12(%d,%f)",modified,errmin);/*DEBUG*/
//    }
//#endif
//  }
//#if L_MODE!=-1
//  /***L_MODE==1の場合，全セルをはじめに初期化しておかないと
//	value[i](=alpha[i])の計算が不安定*/
//  else if(nmax<n_cells){
//    modify_cell(net,nmax,x,y);
//    modified=nmax;
//  }
//#else
//#endif
//  else{/* if(net->cell[s0].v>=net->v_thresh){*/
//    value=(FLOAT*)malloc(sizeof(FLOAT)*n_cells);
//    for(vv=0,i=0;i<n_cells;i++) vv+=net->cell[i].v;
//
//#if L_MODE==1 || L_MODE ==14
//    if(1==0){//original
//	init_AM(&am,2,1);/**/
//	am.x[0]=vv;
//	for(i=0;i<n_cells;i++){
//	  am.x[1]=net->cell[i].v;
//	  am.y[0]=net->cell[i].S;
//	  calc_AM(am);
//	}
//	value_hat=am.M[0][0];
//	sigma2_hat=am.M[0][1];
//	free_AM(&am);
//    }
//    else{//new
//	vt=n_vt=0;
//	for(i=0;i<n_cells;i++)if(vt<net->cell[i].v) vt=net->cell[i].v;
//	sigma2_hat=1e30;
//	vt/=2;
//	for(i=0;i<n_cells;i++){
//	  if(net->cell[i].v<vt) continue;
//	  n_vt++;
//	  if(sigma2_hat>(net->cell[i].S/net->cell[i].v)) sigma2_hat=(net->cell[i].S/net->cell[i].v);
//	}
//    }
//    //printf("sigma2hat=%e\n",sigma2_hat);
//#else
//    value_hat=sigma2_hat=0;
//#endif
//    
//#if L_MODE==1 || L_MODE==14
//    /*alpha_i*/
//    for(value_bar=0,i=0;i<nmax;i++){
//	value[i]=(net->cell[i].S-sigma2_hat*net->cell[i].v)/vv;/*alpha_i 一定*/
//	if(value[i]<0) alpha_NG=1;
//	value_bar+=value[i];
//    }
//#elif L_MODE==2
//    for(value_bar=0,i=0;i<n_cells;i++){
//	value[i]=net->cell[i].S/(net->cell[i].v+1);/*S_i/n_i一定*/
//	value_bar+=value[i];
//    }
//#elif L_MODE==3 || L_MODE==4 || L_MODE==34
//    for(value_bar=0,i=0;i<n_cells;i++){
//	value[i]=net->cell[i].v;/*等発火(n_i)*/
//	value_bar+=value[i];
//    }
//#endif
//#ifdef CHECK
//#if L_MODE==1 || L_MODE==14
//    calc_entropy(value,n_cells);/**/
//    //      /*check_entropy_of_alpha(net);/**/
//#else
//    if(dj==0){
//	calc_entropy(value,n_cells);/**/
//	///*check_entropy_of_alpha(net);/**/
//    }
//#endif      
//#endif      
//    value_bar/=n_vt; // value_bar/=n_cells;//
//    if(1==0){
//	value_var=0;
//	i_value_min=0;
//	value_min=value[0]-value_bar;
//	for(i=0;i<n_cells;i++){
//	  value[i]-=value_bar;
//	  value_var+=value[i]*value[i];
//	  /*      if(value[i]<value_min&&net->cell[i].v>=net->v_thresh){*/
//	  if(value[i]<value_min){
//	    value_min=value[i];
//	    i_value_min=i;
//	  }
//	}
//	value_var=sqrt(value_var/n_cells);
//    }
//    TIME_FOR_REINIT+=1;
//#if L_MODE==1 
//    if(TIME_FOR_REINIT>tau_E/4. && alpha_NG==0 && value[s0]>2.0*value_bar){/**/
//	 TIME_FOR_REINIT=0;
//#elif L_MODE==2 ||L_MODE==3 
//    if(value[s0]>3.0*value_var){/*&&net->cell[s0].v>=net->v_thresh*/
//#elif L_MODE==4
//	  for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
//    if(fabs(err0)>errlimit){/*&&net->cell[s0].v>=net->v_thresh*/
//#elif L_MODE==34
//    for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
//    if(TIME_FOR_REINIT>tau_E/4.&&fabs(err0)>errlimit||value[s0]>2.0*value_bar){
//	  TIME_FOR_REINIT=0;
//#elif L_MODE==14
//    for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];
//    if(nmax<n_cells) {modify_cell(net,nmax,x,y);modified=nmax;}
//    else if(fabs(err0)>errlimit||value[s0]>3.0*value_var){/**/
//	  ///*    else if((fabs(err0)>errlimit||value[s0]>3.0*value_var)&&Entropy<0.985){/**/
//#elif L_MODE==5
//    if(nmax<n_cells) {modify_cell(net,nmax,x,y);modified=nmax;}
//    if(1==0){
//#endif 
//	/*一定化すべき変数が一定でないとき*/
//	if(nmax<n_cells){/* 記憶容量未満なら使ってないセルでそのまま記憶 */
//	  modify_cell(net,nmax,x,y);
//	  modified=nmax;/*??*/
//	}
//	else{
//	  /*最小価値のセルで記憶 
//	    nmax>n_cells&&value[s0]>1.1*value_var&&net->cell[s0].v>=net->v_thresh*/
//	  CELL *c0,*c1;
//	  c0=&net->cell[s0];
//	  c1=&net->cell[i_value_min];
//	  c1->v=0;/*initialize i_value_min*/
//	  init_AMdata(&(net->cell[i_value_min].am));
//	  ///*	if(_TIME<100000)	modify_cell(net,i_value_min,x,y);
//	  printf("%d)Reinit(%d,%e,%e)",GlobalTime,i_value_min,value[s0],value_var);/*DEBUG*/
//	  modify_cell(net,i_value_min,x,y);
//	  //	printf("Reinit(%d,%3.1f)",i_value_min,value[s0]/value_var);/*DEBUG*/
//	  modified=i_value_min;
//	}
//    }
//    else {/*uniform のとき vales[s0]<1.1*value_var&&net->cell[s0].v>=net->v_thresh*/
//	FLOAT delta_wic;
//	int si;
//#if L_MODE==1 || L_MODE==2 || L_MODE==3
//	for(err0=-y,j=0;j<=net->k;j++) err0 += net->cell[s[0]].am.M[0][j]*x[j];//check
//	err0=square(err0);
//#elif L_MODE==4 || L_MODE==34 || L_MODE==14 || L_MODE==5
//	err0=square(err0);
//#endif
//	for(i=1;i<net->n_compare;i++){/*境界の調節*/
//	  sort_weights(net, s, v, i, i+1);//	nmax=sort_memory(x,net,s,v,i+1,nmax,1,i);
//	  si=s[i];
//	  if(net->cell[si].v<net->v_thresh) continue;
//	  //	WindowWidth=0.1;/**/
//	  //	WindowWidth=0.5;
//	  if(in_window(net->cell[si].w,net->cell[s0].w, x, net->width,net->k)>0){/*境界*/
//	    for(delta_wic=0,j=0;j<=net->k;j++) 
//	      delta_wic+=square(net->cell[s0].w[j]-net->cell[si].w[j]);
//	    //	    if((delta_wic=sqrt(delta_wic))<ZERO) delta_wic=ZERO;
//	    if((delta_wic=sqrt(delta_wic))<ZERO) break;/**/
//	    for(erri=-y,j=0;j<=net->k;j++) erri += net->cell[s[i]].am.M[0][j]*x[j];
//	    erri=square(erri);
//	    ///*	  rate_w=-(err0-erri)/delta_wic/0.001;/**f_07: MSE=0.00229494 MODE3*/
//	    rate_w=-(err0-erri)/delta_wic/0.01;/*f_07: MSE=0.00018038 MODE1
//						*f_07: MSE=0.00042767 MODE3****/
//	    ///*	  rate_w=-(err0-erri)/delta_wic/0.1; /*f_07: MSE=0.00013 MODE1***
//	    //					    *f_07: MSE=0.00050 MODE3*/
//	    ///*	  rate_w=-(err0-erri)/delta_wic/1.0; /*f_07: MSE=0.00052056  MODE1*/
////	    if(GlobalTime>=53){//check
////	      printf("err=%e,%e\n",err0,erri);
////	    }
//
//	    modify_cell(net,s0,x,y);/**/
//	    modified=s0;
//	    rate_w=1.0;
//	    break;
//	  }
//	  else continue;
//	}
//    }
//    if(modified<0){/*上で変更されていない時*/
//	/*		  modify_x(net,s0,x,y);*/
//	modify_v(net,s0,x,y);
//	modify_M(net,s0,x,y);
//#ifdef DEBUG
//	if(debug==2){
//	  modified=s0;
//	  printf("\nsv24(%d,%f)",modified,errmin);/*DEBUG*/
//	}
//#endif
//    }
//    free(value);
//  }
//  free(v);
//  free(s);
//  return(0);
//}


/*====================================================================*
 * show_weights
 *
 *====================================================================*/
void save_wM(NET *net)
{
  int i,k,j,n;
  int nEns=net->nEns;
  FILE *fp=fopen("./tmp/M.dat","w");
  for(n=0;n<nEns;n++){
    for(k=0;k<net[n].n_cells;k++){
      for(j=0;j<=net[n].k;j++){
	fprintf(fp,"%+e ",net[n].cell[k].am.M[0][j]);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  fp=fopen("./tmp/w.dat","w");
  for(n=0;n<nEns;n++){
    for(k=0;k<net[n].n_cells;k++){
      for(j=0;j<=net[n].k;j++){
	fprintf(fp,"%+e ",net[n].cell[k].w[i]);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  net->printf("\n*** tmp/M.dat and w.dat are saved ***\n");
}
void save_M(NET *net)
{
  int i,k,j,n;
  int nEns=net->nEns;
  FILE *fp=fopen("./tmp/M.dat","w");
  for(n=0;n<nEns;n++){
    for(k=0;k<net[n].n_cells;k++){
      for(j=0;j<=net[n].k;j++){
	fprintf(fp,"%+e ",net[n].cell[k].am.M[0][j]);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  net->printf(strx2("\n***tmp/M.dat saved N%d b%d***\n",net->n_cells,nEns));
  //  printf("\n*** tmp/M.dat saved N%d b%d ***\n",net->n_cells,net->nEns);//check20191114
  //  system("wc tmp/M.dat");
}
void save_wm(NET *net)
{
  int i,k,j;
  int nEns=net->nEns;
  int n;
  FILE *fp=fopen("./tmp/wm.dat","w");
  for(n=0;n<nEns;n++){
    for(k=0;k<net[n].n_cells;k++){
      for(i=0;i<net[n].k;i++){
	fprintf(fp,"%+e ",net[n].cell[k].w[i]);
      }
      for(j=0;j<net->k;j++){
	fprintf(fp,"%+e ",net[n].cell[k].am.M[0][j]);
      }
      fprintf(fp,"%+e ",net[n].cell[k].am.M[0][net->k]);
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  net->printf("\n***  w and M are saved in tmp/wm.dat. ***\n");
//
//  FLOAT maxwi=-1.;
//  FLOAT minwi=1e3;
//  FLOAT widthwi;
//  FLOAT dM;
//
//  for(k=0;k<net->n_cells;k++){
//    for(i=0;i<net->k;i++){
//      if(net->cell[k].w[i]>maxwi) maxwi=net->cell[k].w[i];
//      if(net->cell[k].w[i]<minwi) minwi=net->cell[k].w[i];
//    }
//  }
//  maxwi+=1.;
//  widthwi=maxwi-minwi;
//  //  maxwi=1.;
//
//  for(k=0;k<net->n_cells;k++){
//    dM=0;
//    for(i=0;i<net->k;i++){
//      fprintf(fp,"%f ",(net->cell[k].w[i]-minwi)/widthwi);
//    }
//    for(j=0;j<net->k;j++){
//      fprintf(fp,"%f ",net->cell[k].am.M[0][j]*widthwi);
//      dM+=(minwi*net->cell[k].am.M[0][j]);
//    }
//    fprintf(fp,"%f ",net->cell[k].am.M[0][net->k]+dM);
//
//    fprintf(fp,"\n");
//    //    fprintf(fp,"%e\n",maxwi);
//  }
//  fclose(fp);

//  fprintf(stdout,"\n***  w and M are saved in tmp/wm.dat. ***\n");
}

void show_weights(NET *net,int nn)
{
#define NEW20181006
#ifdef NEW20181006
  int k,j,kk;
  if(net->n_cells<nn) kk=net->n_cells; else kk=nn;
  for(k=0;k<kk;k++){
    fprintf(stderr,"M[%d]=",k);
    for(j=0;j<=net->k;j++){
      fprintf(stderr,"%+10.3e ",net->cell[k].am.M[0][j]);
    }
    fprintf(stderr,"\n");
  }
#else
  int i,k,j,kk;
  printf("### x0,x1,x2,b0,b1,b2,value\n");
  if(net->n_cells<nn) kk=net->n_cells; else kk=nn;
  //  for(k=0;k<net->n_cells;k++){
  for(k=0;k<kk;k++){
    printf("\n%3d)w=",k);
    for(i=0;i<net->k;i++){
	printf("%+11.3f ",net->cell[k].w[i]);
    }
    printf("\nM=");
    //mat_dumpf(net->cell[k].M,"%+20.10f ");
    for(j=0;j<=net->k;j++){
	printf("%+10.3e ",net->cell[k].am.M[0][j]);
    }
    printf("nVi=%8.1f S=%8.4f", net->cell[k].v*net->vmax,net->cell[k].S);
    //printf("nVi=%8.1f S=%8.4f", net->cell[k].v*net->V->vmax,net->cell[k].S);
    //    printf("V=%8.4f S=%8.4f", net->cell[k].v,net->cell[k].S);
  }
  printf("\nVmax=%d.\n",net->vmax);
#endif
  //  printf("\nVmax=%d.\n",net->V->vmax);
}
//void search_planes(NET *net)
//{
//  int j,nn;
//  int u1,u2;
//  int k=net->k;//input channels;
//  int k1=k+1;
//  int k2=k+2;
//  int nEns=net[0].nEns;//number of ensemble
//  int n_cells=net[0].n_cells;
//  FLOAT *nv=(FLOAT *)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//normal (hosen) vectors
//  FLOAT *nvmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//normal (hosen) vectors
//  int *sameplanes=(int*)malloc(sizeof(int)*nEns*n_cells*n_cells);//normal (hosen) vectors
//  int *color=(int*)malloc(sizeof(int)*nEns*n_cells);//number of planes
//  int *ncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  int *sncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  FLOAT nvk1;
//  FLOAT t1,cost1,sint1;
//  int maxnp=1;
//  char fname[20];
//  int col=0;
//  fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90] and number of planes to extract[<%d]:",n_cells);
//  scanf2("%lf%d",&t1,&maxnp);
//  if(t1<1)t1=1; else if(t1>89) t1=89;
//  cost1=cos(t1*PI/180.);
//  sint1=sin(t1*PI/180.);
//  FILE *fp;
//  for(j=0;j<maxnp;j++){
//    sprintf(fname,"tmp/plane%d.dat",j);
//    fp=fopen(fname,"w");
//    fclose(fp);
//  }
//  for(nn=0;nn<nEns;nn++){//calc normal vectors
//    FLOAT *nvnn=&nv[(nn*n_cells)*k2];
//    FLOAT *nvmeannn=&nvmean[(nn*n_cells)*k2];
//    //    int *sameplanesnn=&sameplanes[(nn*n_cells)*n_cells];
//    int *colornn=&color[nn*n_cells];
//    int *ncolnn=&ncol[nn*n_cells];
//    int *sncolnn=&sncol[nn*n_cells];
//    NET *netnn=&net[nn];
//    for(u1=0;u1<n_cells;u1++){
//      FLOAT sumM2=0;
//      for(j=0;j<k;j++) sumM2+=square(netnn->cell[u1].am.M[0][j]);
//      nvk1=nvnn[u1*k2+k]=sqrt(1.0/(1.0+sumM2));
//      for(j=0;j<k;j++)
//	nvnn[(u1)*k2+j]=-netnn->cell[u1].am.M[0][j]*nvk1;
//      
//      nvnn[(u1)*k2+k1]=-netnn->cell[u1].am.M[0][k]*nvk1;//distance from the origin;
//      colornn[u1]=0;
//      ncolnn[u1]=0;
//    }
//    for(u1=0;u1<n_cells;u1++){
//      if(colornn[u1]!=0) continue;
//      col++;
//      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] = nvnn[(u1)*k2+j];
//      ncolnn[col]++;
//
//      for(u2=0;u2<n_cells;u2++){
//	if(colornn[u2]!=0) continue;
//	else{
//	  FLOAT n1n2=0.0;
//	  for(j=0;j<k1;j++) n1n2+=(nvnn[(u1)*k2+j]*nvnn[(u2)*k2+j]);
//	  if(n1n2>=cost1){
//	    FLOAT n1w12=0.0;
//	    FLOAT n2w12=0.0;
//	    FLOAT w12=0,w12j;
//	    FLOAT y1=0;
//	    FLOAT y2=0;
//	    for(j=0;j<k;j++){
//	      w12j=netnn->cell[u1].w[j]-netnn->cell[u2].w[j];
//	      w12   += w12j*w12j;
//	      n1w12 += (nvnn[(u1)*k2+j]*w12j);
//	      n2w12 += (nvnn[(u2)*k2+j]*w12j);
//	    }
//	    {
//	      for(j=0;j<k;j++){
//		y1 += (netnn->cell[u1].am.M[0][j]*netnn->cell[u1].w[j]);
//		y2 += (netnn->cell[u2].am.M[0][j]*netnn->cell[u2].w[j]);
//	      }
//	      y1 += (netnn->cell[u1].am.M[0][k]);//constant element
//	      y2 += (netnn->cell[u2].am.M[0][k]);//constant element
//
//	      w12j=y1-y2;
//	      w12 += w12j*w12j;
//	      n1w12 += (nvnn[(u1)*k2+k]*w12j);
//	      n2w12 += (nvnn[(u2)*k2+k]*w12j);
//	    }
//	    w12=sqrt(w12);
//	    n1w12=fabs(n1w12/w12);
//	    n2w12=fabs(n2w12/w12);
//
//	    if(n1w12<=sint1 && n1w12<=sint1){
//	      colornn[u2]=col;
//	      ncolnn[col]++;
//	      for(j=0;j<=k1;j++) nvmeannn[col*k2+j]+= nvnn[(u2)*k2+j];
//	      //	    if(u1<10) fprintf(stderr,"check u1=%d,u2=%d.\n",u1,u2);
//	    }
//	  }
//	}
//      }
//      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] /=(ncolnn[col]);
//      //      fprintf(stderr,"check npnn[%d]=%d\n",u1,npnn[u1]);
//    }
//    int i;
//    int colmax=col;
//    //    if(colmax<maxnp) maxnp=colmax;
//    for(i=0;i<colmax;i++) sncolnn[i]=i;
//    for(i=0;i<maxnp-1;i++){//sort
//      for(j=i+1;j<colmax;j++){
//	if(ncolnn[sncolnn[j]]>ncolnn[sncolnn[i]]){
//	  int ntmp=sncolnn[j];
//	  sncolnn[j]=sncolnn[i];
//	  sncolnn[i]=ntmp;
//	}
//      }
//    }
//    for(i=0;i<maxnp;i++){//
//      col=sncolnn[i];
//      if(col==0) continue;//no color
//      fprintf(stdout,"Plane#%d (number of units for this plane=%d;col=%d);",i,ncolnn[col],col);
//      fprintf(stdout,"\nNormal Vector (Hohsen)=");
//
//      for(j=0;j<k1;j++) fprintf(stdout,"%+.5f ",nvmeannn[col*k2+j]);
//
//      fprintf(stdout,"Distance from origin=%.5f\n",nvmeannn[col*k2+k1]);
//
//      sprintf(fname,"tmp/plane%d.dat",i);
//      fp=fopen(fname,"a+");
//
//      for(u1=0;u1<n_cells;u1++){
//	if(colornn[u1]==col){
//	  FLOAT *w=&(netnn->cell[u1].w[0]);
//	  fprintf(stderr,"u%d",u1);//check
//	  for(j=0;j<k;j++){
//	    fprintf(fp,"%e ",w[j]);
//	  }
//	  FLOAT y=0;
//	  for(j=0;j<k;j++) y += (-nvmean[col*k2+j]/nvmean[col*k2+k])*w[j];
//	  y+=(-nvmean[col*k2+k1]/nvmean[col*k2+k]);
////	  {
////	    FLOAT *x=(FLOAT*)malloc(sizeof(FLOAT)*k1),_yt,_y,_Yt;
////	    for(j=0;j<k;j++) x[j]=w[j];
////	    x[k]=1;
////	    calc_output(&net[nn], x,&_yt,&_y,&_Yt);
////	    y=_Yt;
////	  }
//
//	  fprintf(fp,"%e\n",y);
//	  if(col==1){//check
//	    fprintf(stderr,"\nw=");
//	    for(j=0;j<k;j++) fprintf(stderr,"%f ",netnn->cell[u1].w[j]);
//	    fprintf(stderr,"y=%f",y);
//
//	    fprintf(stderr,"\nm=");
//	    for(j=0;j<k;j++) fprintf(stderr,"%f ",(-nvmean[col*k2+j]/nvmean[col*k2+k]));
//	    fprintf(stderr,"%f ",(-nvmean[col*k2+k1]/nvmean[col*k2+k]));
//
//
//	    fprintf(stderr,"\nM=");
//	    for(j=0;j<=k;j++) fprintf(stderr,"%f ",netnn->cell[u1].am.M[0][j]);
//	    fprintf(stderr,"v=%f\nn=",netnn->cell[u1].v);
//
//	    for(j=0;j<=k1;j++) fprintf(stderr,"%f ",nvnn[(u1)*k2+j]);
//	    fprintf(stderr,"\n");
//
//
//	  }
//	}
//      }
//      fclose(fp);
//    }
//  }
//  //  fclose(fp);
//  free(nv);
//  free(nvmean);
//  free(sameplanes);
//  free(color);
//  free(ncol);
//  free(sncol);
//}
//void show_wy(NET *net)
//{
//  int i,k,j,kk;
//  printf("### y and x0,x1,x2,...\n");
//  for(k=0;k<net->n_cells;k++){
//    calc_output(net,net->cell[k].w,&yr,&y,&Y);
//    printf("%+e ",Y);
//    for(i=0;i<net->k;i++){
//	  printf("%+e ",net->cell[k].w[i]);
//    }
//    printf("\n");
//  }
//}

/*====================================================================*
 * change_net_parameters
 *
 *====================================================================*/
/*
change_net_parameters(NET* net) {
  sprintf(ffff,"%s%s%s%s%s","%d%d",LF,LF,LF,LF);
  printf( "v_thresh,n_compare,errlimit,ZERO_EQUIV,ZERO(%d,%d,%e,%e,%e):",
	   net->v_thresh,net->n_compare,errlimit,ZERO_EQUIV,ZERO);
  //myscanf("%d%d%lf%lf%lf",&net->v_thresh,&net->n_compare,&errlimit,&ZERO_EQUIV,&ZERO);
  //myscanf(ffff,&net->v_thresh,&net->n_compare,&errlimit,&ZERO_EQUIV,&ZERO,&net->v_ratio);
  if(net->n_compare>=net->n_cells) net->n_compare=net->n_cells;
}
*/


/*====================================================================*
 * get_w_batch
 * 荷重ベクトルを取りだす
 *====================================================================*/
void get_w_batch(NET *net, FLOAT **w) {
  int n_cells = net->n_cells, n_channels = net->k;
  int i, k;

  for (i=0; i<n_cells; i++)
    for (k=0; k<n_channels; k++)
      w[i][k] = net->cell[i].w[k];

  return;
}

/*====================================================================*
 * set_w_batch
 *
 *====================================================================*/
void set_w_batch(NET *net, FLOAT **w) {
  int n_cells = net->n_cells;
  int n_channels = net->k;
  //  int n_ivectors_thresh;
  int i, k;
  
  //  n_ivectors_thresh = (int)((n_channels+1)*2);
  
  for (i=0; i<n_cells; i++)
    for (k=0; k<n_channels; k++)
      net->cell[i].w[k] = w[i][k];
  
  return;
}

///*====================================================================*
// * calc_value_batch1
// *
// *====================================================================*/
//inline void calc_value_batch1(NET *net) {
//  int n_cells = net->n_cells;
//  int n_fires = net->n_fires;
//  FLOAT v_total, alpha_hat, sigma2_hat, alpha_min, alpha_bar, alpha_var;
//  int i_alpha_min, alpha_NG;
//  AM am;
//  int i;
//
//  // v_total（n_hat）, alpha_hat, sigma^2_hat を求める
//  v_total = 0.;
//  for (i=0; i<n_cells; i++) v_total += net->cell[i].v;
//
//  init_AM(&am, 2, 1);
//  for (i=0; i<n_cells; i++) {
//    /*
//    printf("> v[%d]: %3.2e, S[%d]: %3.2e, enough?: %d\n",
//	     i, net->cell[i].v, i, net->cell[i].S, net->cell[i].enough);
//    */
//    am.x[0] = v_total/net->cell[i].v;
//    am.x[1] = 1.0/net->cell[i].v;
//    am.y[0] = net->cell[i].S;
//    calc_AM(am);
//  }
//  alpha_hat = am.M[0][0];
//  sigma2_hat = am.M[0][1];
//  free_AM(&am);
//  net->alpha_hat = alpha_hat;
//  net->sigma2_hat = sigma2_hat;
//
//  // 価値alpha_iの平均値を求める
//  alpha_bar = 0.;
//  alpha_NG = 0;
//  for (i=0; i<n_fires; i++) {
//    net->cell[i].alpha = net->cell[i].S*net->cell[i].v/v_total-sigma2_hat/v_total;
//    alpha_bar += net->cell[i].alpha;
//
//    // BAD!! Miscalculate the Value for Each Unit
//    if (net->cell[i].alpha < 0.0) alpha_NG = 1;
//  }
//  alpha_bar /= (FLOAT)n_cells;
//  net->alpha_bar = alpha_bar;
//  net->alpha_NG = alpha_NG;
//
//  // 価値alpha_iを求める
//  i_alpha_min = 0;
//  //alpha_min = net->cell[0].alpha - alpha_bar;
//  alpha_min = net->cell[0].alpha;
//  for (i=0; i<n_cells; i++) {
//    // Caution!! Every Value is Substracted from Average Value 
//    //net->cell[i].alpha -= alpha_bar;
//    if (net->cell[i].alpha < alpha_min) {
//	i_alpha_min = i;
//	alpha_min = net->cell[i].alpha;
//    }
//  }
//  net->l = i_alpha_min;
//  net->alpha_min = alpha_min;
//
//  //printf("> v_total: %3.2e, alpha_hat: %3.2e, sigma2_hat %3.2e, alpha_bar: %3.2e\n",
//  // v_total, alpha_hat, sigma2_hat, alpha_bar);
//
//  return;
//}
//
///*====================================================================*
// * calc_value_batch2
// *
// *====================================================================*/
//inline void calc_value_batch2(NET *net) {
//  int n_cells = net->n_cells;
//  int n_fires = net->n_fires;
//  FLOAT v_total, alpha_hat, sigma2_hat, alpha_min, alpha_bar, alpha_var;
//  int i_alpha_min, alpha_NG;
//  AM am;
//  int i;
//
//  // v_total（n_hat）, alpha_hat, sigma^2_hat を求める
//  v_total = 0.;
//  for (i=0; i<n_cells; i++) v_total += net->cell[i].v;
//
//  init_AM(&am, 2, 1);
//  for (i=0; i<n_cells; i++) {
//    /*
//    printf("> v[%d]: %3.2e, S[%d]: %3.2e, enough?: %d\n",
//	     i, net->cell[i].v, i, net->cell[i].S, net->cell[i].enough);
//    */
//    am.x[0] = 1.0;
//    am.x[1] = 1.0/v_total;
//    am.y[0] = net->cell[i].S;
//    calc_AM(am);
//  }
//  alpha_hat = am.M[0][0];
//  sigma2_hat = am.M[0][1];
//  free_AM(&am);
//  net->alpha_hat = alpha_hat;
//  net->sigma2_hat = sigma2_hat;
//
//  // 価値alpha_iの平均値を求める
//  alpha_bar = 0.;
//  alpha_NG = 0;
//  for (i=0; i<n_fires; i++) {
//    net->cell[i].alpha = net->cell[i].S-sigma2_hat/v_total;
//    alpha_bar += net->cell[i].alpha;
//
//    // BAD!! Miscalculate the Value for Each Unit
//    if (net->cell[i].alpha < 0.0) alpha_NG = 1;
//  }
//  alpha_bar /= (FLOAT)n_cells;
//  net->alpha_bar = alpha_bar;
//  net->alpha_NG = alpha_NG;
//
//  // 価値alpha_iを求める
//  i_alpha_min = 0;
//  //alpha_min = net->cell[0].alpha - alpha_bar;
//  alpha_min = net->cell[0].alpha;
//  for (i=0; i<n_cells; i++) {
//    // Caution!! Every Value is Substracted from Average Value 
//    //net->cell[i].alpha -= alpha_bar;
//    if (net->cell[i].alpha < alpha_min) {
//	i_alpha_min = i;
//	alpha_min = net->cell[i].alpha;
//    }
//  }
//  net->l = i_alpha_min;
//  net->alpha_min = alpha_min;
//
//  //printf("> v_total: %3.2e, alpha_hat: %3.2e, sigma2_hat %3.2e, alpha_bar: %3.2e\n",
//  // v_total, alpha_hat, sigma2_hat, alpha_bar);
//
//  return;
//}

/*====================================================================*
 * modify_M_batch
 * 連想行列の更新
 *====================================================================*/
//#define GSL
// make clean;make "GSL=0"
// make clean;make "GSL=1"
// make clean;make "GSL=2"
//#undef GSL
#if GSL >= 1
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
int GSLtag;
void least_gsl_error_handler( const char * reason,
			      const char * file,
			      int line,
			      int gsl_errno)
{
  fprintf(stderr,"%s:%d: %s (gsl_error: %s)\n", 
	 file, line, reason, gsl_strerror(gsl_errno));
  GSLtag=1;
}
void calc_Ainvb(FLOAT *M, FLOAT *a_data, FLOAT *b_data, int nx, int ndata)
{
  int i;
  gsl_vector *S;
  gsl_vector *work;
  gsl_vector *x;
  gsl_vector *b;
  gsl_matrix *V;
  gsl_matrix *A;
  gsl_matrix_view Av= gsl_matrix_view_array(a_data, ndata, nx);
  gsl_vector_view bv = gsl_vector_view_array(b_data, ndata);;

  A =&(Av.matrix);
  b =&(bv.vector);

  V = gsl_matrix_alloc (nx,nx); 
  S = gsl_vector_alloc (nx); 
  work = gsl_vector_alloc(nx); 
  x = gsl_vector_alloc(nx); 

  gsl_set_error_handler( &least_gsl_error_handler );GSLtag=0;

  //  printf("laptime:gsl_linalg_SV_decomp=%5.3f\n",mytimer_lap());
  //#define GSL 2
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
  //    printf("laptime:gsl_linalg_SV_solve=%5.3f\n",mytimer_lap());
  gsl_linalg_SV_solve(A,V,S,b,x);
  for(i=0;i<nx;i++) M[i]=gsl_vector_get(x,i);
  if(GSLtag==1){
    //    printf("\nEigenValue=");for(i=0;i<nx;i++) printf("%e ",S[i]); printf("\n");
    //    for(i=0;i<nx;i++) M[i]=0;
  }

  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_vector_free(x);
  gsl_matrix_free(V);
}
void modify_M_batch_pinv(NET *net, FLOAT **x, FLOAT *y) {
  int n_channels = net->k;
  int i,j,k, t;
  int k1=net->k+1;
  int j1;
  int vmax;
  VORONOI *V=net->V;
  FLOAT *A;
  FLOAT *b;

  vmax=0;
  for(i=0;i<net->n_cells;i++) if(vmax<V->i2v[i]) vmax=V->i2v[i];
  A=(FLOAT *)malloc(sizeof(FLOAT)*k1*vmax);
  b=(FLOAT *)malloc(sizeof(FLOAT)*vmax);

  for(i=0;i<net->n_cells;i++){
    j1=V->i2v[i];
    if(j1<=0){
      //      printf("????????????????????modify_M_batch????????????????????\n");
      continue;
      j1=1;
      for (j=0; j<j1; j++) {
	for (k=0; k<=n_channels; k++) A[j*k1+k] = 0;
	b[j]=0;
      }
    }
    else{
      if(j1<k1) j1=k1;
      for (j=0; j<j1; j++) {
	//    for (j=0; j<V->i2v[i]; j++) {
	t=V->ij2t[i][j%V->i2v[i]];
	//      t=V->ij2t[i][j];
	for (k=0; k<=n_channels; k++) A[j*k1+k] = x[t][k];
	b[j]=y[t];
      }
    }
    //    calc_Ainvb(net->cell[i].am.M[0],A,b,k1,V->i2v[i]);
    calc_Ainvb(net->cell[i].am.M[0],A,b,k1,j1);
//  if(GlobalTime==5 && i==3){
//    //	printf("\n%d)y%ex",j,b[j]);
//    //	for (k=0; k<=n_channels; k++) printf("%e ",A[j*k1+k]);
//    printf("\nM=");
//    for (k=0; k<=n_channels; k++) printf("%e ",net->cell[i].am.M[0][k]);
//  }
    
  }
  free(A);
  free(b);
  return;
}
//void modify_M_batch_hen(NET *net, FLOAT **x, FLOAT *y) {
//  int n_channels = net->k;
//  int i,j,k, t;
//  int k1=net->k+1;
//  int j1;
//  int vmax;
//  VORONOI *V=net->V;
//  FLOAT *A;
//  FLOAT *b;
//
//  vmax=0;
//  for(i=0;i<net->n_cells;i++) if(vmax<V->i2v[i]) vmax=V->i2v[i];
//  A=(FLOAT *)malloc(sizeof(FLOAT)*k1*vmax);
//  b=(FLOAT *)malloc(sizeof(FLOAT)*vmax);
//  
//  for(i=0;i<net->n_cells;i++){
//    if(V->i2v[i]<=0) continue;
//    j1=j<V->i2v[i];
//    if(j1<k1) j1=k1;
//    for (j=0; j<j1; j++) {
//	//    for (j=0; j<V->i2v[i]; j++) {
//	t=V->ij2t[i][j%V->i2v[i]];
//	//      t=V->ij2t[i][j];
//	for (k=0; k<=n_channels; k++) A[j*k1+k] = x[t][k];
//	b[j]=y[t];
//    }
//    //    calc_Ainvb(net->cell[i].am.M[0],A,b,k1,V->i2v[i]);
//    calc_Ainvb(net->cell[i].am.M[0],A,b,k1,j1);
//  }
//  free(A);
//  free(b);
//  return;
//}
#endif
//#else //#if GSL >= 1
void modify_M_batch_RLS(NET *net, FLOAT **x, FLOAT *y) {
  int n_channels = net->k;
  int i,j,k, t;
  VORONOI *V=net->V;

  for(i=0;i<net->n_cells;i++){
    //init_AMdata(&(net->cell[i].am));//kuro040115test
    for (j=0; j<V->i2v[i]; j++) {
    //for (j=0; j<V->i2v0[i]; j++) {
      t=V->ij2t[i][j];
      for (k=0; k<=n_channels; k++) net->cell[i].am.x[k] = x[t][k];
      net->cell[i].am.y[0] = y[t];
      calc_AM(net->cell[i].am);
      //check 20190115
//	{//debug
//	  printf("\n");
//	  for (k=0; k<=n_channels; k++) printf("%e ",x[t][k]);
//	  printf("%e # %d", y[t],t);
//	}      
    }
  }
  //check M 20190115
#ifdef DEBUG
  if(debug==1){
    net->printf("finish modify_M_batch\n");
    for(i=0;i<net->n_cells;i++){
      net->printf(strx4("%d w%g M%g %g\n",i,net->cell[i].w[0],net->cell[i].am.M[0][0],net->cell[i].am.M[0][1]));
    }
  }
#endif
  return;
}
//#endif //#if GSL >= 1
//void modify_M_batch0(NET *net, int i, FLOAT **x, FLOAT *y, int n_ivectors) {
//  int n_channels = net->k;
//  int j, k;
//
//  //  init_AMdata(&(net->cell[i].am));//kuro
//  for (j=0; j<n_ivectors; j++) {
//    for (k=0; k<=n_channels; k++) net->cell[i].am.x[k] = x[j][k];
//    net->cell[i].am.y[0] = y[j];
//    calc_AM(net->cell[i].am);
////    if(i==39){//kuro:debug
////	printf("y[%d]=%e",j,y[j]);
////	for (k=0; k<=(n_channels); k++){
////	  printf("(%3.2e)*(%3.2e)",net->cell[i].am.M[0][k],x[j][k]);
////	}
////	printf("\n");
////    }
//  }
//  return;
//}
//20191025#include <gsl/gsl_poly.h>
//20191025void calc_poles_of_M(NET *net)
//20191025{
//20191025  int i,n,j,k1,k,nn;
//20191025  double *a,*z;
//20191025  //  double *zmean;
//20191025  char *fname="./tmp/poles.dat";
//20191025  char *fname2="./tmp/poles2.dat";
//20191025  char *fname3="./tmp/poles3.dat";
//20191025  FILE *fp=fopen(fname,"w");
//20191025  FILE *fp2=fopen(fname2,"w");
//20191025  FILE *fp3=fopen(fname3,"w");
//20191025  int sign;
//20191025  gsl_poly_complex_workspace *w;
//20191025  gsl_poly_complex_workspace *w2;
//20191025  k=net[0].k;
//20191025  k1=k+1;
//20191025  w= gsl_poly_complex_workspace_alloc (k1);//for poles
//20191025  w2= gsl_poly_complex_workspace_alloc (k1+1);//for LSP
//20191025  z=(double*)malloc(sizeof(double)*(k+1)*2);//for LSP
//20191025  a=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025  for(nn=0;nn<net[0].nEns;nn++){
//20191025    //poles of M
//20191025    for(n=0;n<net[nn].n_cells;n++){
//20191025      a[k]=-1.;
//20191025      for(j=0;j<k;j++) a[j]=net[nn].cell[n].am.M[0][k-1-j];
//20191025      gsl_poly_complex_solve (a, k1, w, z);
//20191025      for (i = 0; i < k; i++){
//20191025	fprintf (fp,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*i], z[2*i+1],sqrt(z[2*i]*z[2*i]+z[2*i+1]*z[2*i+1]),atan2(z[2*i+1],z[2*i]),net[nn].cell[n].v,i,k,nn);
//20191025	fprintf (fp2,"%+e %+e ",z[2*i], z[2*i+1]);
//20191025      }
//20191025      fprintf (fp2,"%+e %d %d\n",net[nn].cell[n].v,k,nn);
//20191025      //LSP:Line Spectral Pairs
//20191025      sign=1;a[k+1]=-1.;a[0]=-sign;
//20191025      for(j=0;j<k;j++) a[j+1]=(net[nn].cell[n].am.M[0][k-1-j]+sign*net[nn].cell[n].am.M[0][j]);
//20191025      //      for(i=0;i<=k+1;i++) printf("a[%d]=%e\n",i,a[i]);
//20191025      gsl_poly_complex_solve (a, k1+1, w2, z);
//20191025      for (i = 0; i <=k; i++){
//20191025	fprintf (fp3,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*i], z[2*i+1],sqrt(z[2*i]*z[2*i]+z[2*i+1]*z[2*i+1]),atan2(z[2*i+1],z[2*i]),net[nn].cell[n].v,i,k,nn);
//20191025      }
//20191025      sign=-1;a[k+1]=-1.;a[0]=-sign;
//20191025      for(j=0;j<k;j++) a[j+1]=(net[nn].cell[n].am.M[0][k-1-j]+sign*net[nn].cell[n].am.M[0][j]);
//20191025      //      for(i=0;i<=k+1;i++) printf("a[%d]=%e\n",i,a[i]);
//20191025      gsl_poly_complex_solve (a, k1+1, w2, z);
//20191025      for (i = 0; i <=k; i++){	//	printf("a[%d]=%e\n",i,a[i]);
//20191025	fprintf (fp3,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*i], z[2*i+1],sqrt(z[2*i]*z[2*i]+z[2*i+1]*z[2*i+1]),atan2(z[2*i+1],z[2*i]),net[nn].cell[n].v,i,k,nn);
//20191025      }
//20191025    }
//20191025  }
//20191025  fclose(fp3);
//20191025  fclose(fp2);
//20191025  fclose(fp);
//20191025  fprintf(stdout,"Poles of M are saved in '%s','%s','%s'.",fname,fname2,fname3);
//20191025  gsl_poly_complex_workspace_free (w);
//20191025  gsl_poly_complex_workspace_free (w2);
//20191025  free(a);
//20191025  free(z);
//20191025}
//20191025
//20191025void calc_poles_of_Mmean(NET *net)
//20191025{
//20191025  int i,n,j,nn;
//20191025  double rz2,vs;
//20191025  char *fname="./tmp/poles.dat";
//20191025  char *fname2="./tmp/poles2.dat";
//20191025  char *fname3="./tmp/poles3.dat";
//20191025  FILE *fp=fopen(fname,"w");
//20191025  FILE *fp2=fopen(fname2,"w");
//20191025  FILE *fp3=fopen(fname3,"w");
//20191025  int sign;
//20191025  int NG;
//20191025  int k=net[0].k;
//20191025  int k1=k+1;
//20191025  double *z=(double*)malloc(sizeof(double)*(2*k1)*2);//for LSP of PQ
//20191025  double *As=(double*)malloc(sizeof(double)*(k1));//for LPC
//20191025  //  gsl_poly_complex_workspace *w;//
//20191025
//20191025  {    //poles of Mmean
//20191025    double *A=(double*)malloc(sizeof(double)*(k1));//for LPC
//20191025    gsl_poly_complex_workspace *wk= gsl_poly_complex_workspace_alloc (k+1);//for poles
//20191025    A[k]=-1.;vs=0;
//20191025    for(j=0;j<=k;j++) As[j]=0;
//20191025    for(nn=0;nn<net[0].nEns;nn++){
//20191025      for(n=0;n<net[nn].n_cells;n++){
//20191025	//poles of M
//20191025	if(net[nn].cell[n].v>0){
//20191025	  for(j=0;j<k;j++) A[j]=net[nn].cell[n].am.M[0][k-1-j];
//20191025	  gsl_poly_complex_solve (A, k+1, wk, z);
//20191025	  NG=0;
//20191025	  for(j=0;j<k;j++){
//20191025	    rz2=z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1];
//20191025	    //	  if(rz2>1 || rz2<0.707) {NG=1; break;}
//20191025	    if(rz2>=1 || rz2<0.5) {NG=1; break;}
//20191025	    //	  if(rz2<0.5) {NG=1; break;}
//20191025	    //	  if(rz2>1) {NG=1; break;}
//20191025	  }
//20191025	  //	if(1==1 || NG==0){
//20191025	  if(1==1 || NG==0){
//20191025	    //	    for(i=0;i<=k;i++) As[i]+=A[i]*net[nn].cell[n].v;
//20191025	    for(i=0;i<=k;i++) As[i]+=A[i]*net[nn].cell[n].v;
//20191025	  }
//20191025	  vs+=net[nn].cell[n].v;
//20191025	}
//20191025      }
//20191025    }
//20191025    for(i=0;i<k1;i++) As[i]/=vs;
//20191025    //poles of Mmean
//20191025    gsl_poly_complex_solve (As, k+1, wk, z);
//20191025    gsl_poly_complex_workspace_free (wk);
//20191025    for(j=0;j<k;j++){
//20191025      fprintf (fp,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) i k\n", 
//20191025	       z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025    }
//20191025    {
//20191025      FILE *fp4m=fopen("./tmp/poles4m.dat","w");//original poles, poles in each line
//20191025      for(j=0;j<k;j++){
//20191025	for(j=0;j<k;j++) fprintf(fp4m,"%e %e ",z[2*j],z[2*j+1]);
//20191025	fprintf(fp4m,"%e\n",net[nn].cell[n].v);
//20191025      }
//20191025      fclose(fp4m);
//20191025    }
//20191025    free(A);
//20191025  }
//20191025  //
//20191025  {
//20191025    double *P=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025    double *Q=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025    //LSP:Line Spectral Pairs
//20191025    if(1==0){
//20191025      gsl_poly_complex_workspace *wk1= gsl_poly_complex_workspace_alloc (k1+1);//for LSP
//20191025      sign=1;P[k1]=-1.;P[0]=-sign;
//20191025      for(j=0;j<k;j++) P[j+1]=As[j]+sign*As[k-1-j];   //  for(i=0;i<=k1;i++) printf("P[%d]=%e\n",i,P[i]);
//20191025      //    w= gsl_poly_complex_workspace_alloc (k2);//for poles
//20191025      gsl_poly_complex_solve (P, k1+1, wk1, z);
//20191025      for(j=0;j<k1;j++){
//20191025	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025      }
//20191025      sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025      for(j=0;j<k;j++) Q[j+1]=As[j]+sign*As[k-1-j]; //  for(i=0;i<=k1;i++) printf("Q[%d]=%e\n",i,P[i]);
//20191025      //    w= gsl_poly_complex_workspace_alloc (k2);//for poles
//20191025      gsl_poly_complex_solve (Q, k1+1, wk1, z);
//20191025      gsl_poly_complex_workspace_free (wk1);
//20191025      for(j=0;j<k1;j++){
//20191025	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025      }
//20191025      gsl_poly_complex_workspace_free (wk1);
//20191025      free(P);
//20191025      free(Q);
//20191025    }
//20191025    
//20191025    if(1==1){//k3=2*k1+1
//20191025      double *PQ=(double*)malloc(sizeof(double)*(2*k1+1));//for LSP
//20191025      gsl_poly_complex_workspace *wk2=gsl_poly_complex_workspace_alloc (2*k1+1);//for LSP
//20191025      sign=1;P[k1]=-1.;P[0]=-sign;
//20191025      for(j=0;j<k;j++) P[j+1]=As[j]+sign*As[k-1-j];
//20191025      sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025      for(j=0;j<k;j++) Q[j+1]=As[j]+sign*As[k-1-j];
//20191025      for(j=0;j<=2*k1;j++) PQ[j]=0;
//20191025      for(j=0;j<=k1;j++){
//20191025	for(i=0;i<=k1;i++){
//20191025	  PQ[j+i] +=P[i]*Q[j];
//20191025	}
//20191025      }
//20191025      //      for(i=0;i<2*k1;i++) printf("PQ[%d]=%e\n",i,PQ[i]);
//20191025      gsl_poly_complex_solve (PQ, 2*k1+1, wk2, z);
//20191025      //      gsl_poly_complex_workspace_free (w);
//20191025      //    gsl_poly_complex_solve (Q, k1+1, w2, z);
//20191025      for(j=0;j<2*k1;j++){
//20191025	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025      }
//20191025      gsl_poly_complex_workspace_free(wk2);
//20191025      free(PQ);
//20191025    }
//20191025  }
//20191025  free(As);
//20191025  fclose(fp3);
//20191025  fclose(fp2);
//20191025  fclose(fp);
//20191025  fprintf(stdout,"Poles of M are saved in '%s','%s'.",fname,fname3);
//20191025  free(z);
//20191025}
//20191025
//##########
//20191025void poles_of_M_shrink(NET *net)
//20191025{
//20191025  int i,n,j,nn;
//20191025  FILE *fp=fopen("./tmp/poles.dat","w");
//20191025  FILE *fp0=fopen("./tmp/poles0.dat","w");//original poles, a pole in a separate line
//20191025  FILE *fp2=fopen("./tmp/poles2.dat","w");//shrinked LPC poles in each CAN2, a pole in a separate line
//20191025  FILE *fp3=fopen("./tmp/poles3.dat","w");//shrinked LSP poles, a pole in a separate line
//20191025  FILE *fp4=fopen("./tmp/poles4.dat","w");//original poles, poles in each line
//20191025  FILE *fp5=fopen("tmp/poles5.dat","w");//shrinked poles, poles in each line
//20191025  FILE *fp6=fopen("./tmp/poles6.dat","w");//shrinced LSP
//20191025  FILE *fp7=fopen("./tmp/M4.dat","w");//original LPC (Associative Matrix), LPCs in each line
//20191025  int k=net[0].k;
//20191025  int k1=k+1;
//20191025  int nEns=net[0].nEns;
//20191025  int maxk=2*k1;
//20191025  double *vn=(double*)malloc(sizeof(double)*(nEns+1));//number of data in Vi
//20191025  double *An=(double*)malloc(sizeof(double)*(k1)*(nEns+1));//Coefficients of LPC
//20191025  double *zn=(double*)malloc(sizeof(double)*maxk*(nEns+1));//for LSP of PQ
//20191025  double *theta_zn=(double*)malloc(sizeof(double)*maxk*(nEns+1));//for LPC
//20191025  double *rho_zn=(double*)malloc(sizeof(double)*maxk*(nEns+1));//for LPC
//20191025  int *i_zn=(int*)malloc(sizeof(int)*maxk*(nEns+1));//for LPC
//20191025  int *c_zn=(int*)malloc(sizeof(int)*maxk*(nEns+1));//for LPC
//20191025  double *v_zn=(double*)malloc(sizeof(double)*maxk);//for LSP of PQ
//20191025  gsl_poly_complex_workspace *wk2=gsl_poly_complex_workspace_alloc (2*k1+1);//for LSP
//20191025  gsl_poly_complex_workspace *wk= gsl_poly_complex_workspace_alloc (k+1);//for poles
//20191025  {
//20191025    //poles of M[nn]
//20191025    //each pole zn[nn]
//20191025    double *A=(double*)malloc(sizeof(double)*(k1));//for LPC
//20191025    double *z=(double*)malloc(sizeof(double)*maxk);//for LSP of PQ
//20191025    double *P=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025    double *Q=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025    double *PQ=(double*)malloc(sizeof(double)*(2*k1+1));//for LSP
//20191025    int sign;
//20191025
//20191025    A[k]=-1.;
//20191025#define SELECTMETHOD 2
//20191025#if SELECTMETHOD == 2
//20191025    double *v_backup=(double*)malloc(sizeof(double)*nEns*net[0].n_cells);//number of data in Vi
//20191025    for(nn=0;nn<nEns;nn++) for(n=0;n<net[nn].n_cells;n++) v_backup[n*nEns+nn]=net[nn].cell[n].v;
//20191025    for(nn=0;nn<nEns;nn++){
//20191025      vn[nn]=0;for(j=0;j<=k;j++) An[nn*k1+j]=0;
//20191025      for(n=0;n<net[nn].n_cells;n++){
//20191025	//poles of M
//20191025	for(j=0;j<k;j++) A[j]=net[nn].cell[n].am.M[0][k-1-j];
//20191025	gsl_poly_complex_solve (A, k+1, wk, z);
//20191025	for(j=0;j<k;j++){
//20191025	  fprintf (fp0,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		   z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),net[nn].cell[n].v,j,k,nn);
//20191025	}
//20191025	//	fprintf(stderr,"check v=%.0f=%f*%d\n",net[nn].cell[n].v*net[nn].vmax,net[nn].cell[n].v,net[nn].vmax);
//20191025	//	if(net[nn].cell[n].v*net[nn].vmax > k ){
//20191025	if(net[nn].cell[n].v>0){
//20191025	  int NG=0;
//20191025	  for(j=0;j<k;j++){
//20191025	    double rz2=z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1];
//20191025	    if(rz2>1.0 || rz2<0.35) {NG=1; break;}
//20191025	    //	    if(rz2>1.0) {NG=1; break;}
//20191025	    //	    if(rz2>1.0 || rz2<0.3) {NG=1; break;}
//20191025	    //	    if(rz2<0.3) {NG=1; break;}//OK but few features?
//20191025	  }
//20191025	  if(net[nn].n_cells !=1 && NG==1) {
//20191025	    net[nn].cell[n].v=2;//v>1 for excluding calculation of MSE. see exec_ssp_test_ensemble in sim.c 
//20191025	    //	    fprintf(stdout,"check: net[%d].cell[%d].v=2 is not used for evaluating MSE.\n",nn,n);
//20191025	  }
//20191025	  //	    if(net[nn].n_cells !=1 && NG==1) net[nn].cell[n].v=-1;//v<0 for neglecting this unit
//20191025	  //	  if(net[nn].n_cells==1 || NG==0){
//20191025	  if(1==1){
//20191025	    for(j=0;j<=k;j++) An[nn  *k1+j]+=A[j]*net[nn].cell[n].v;
//20191025	    vn[nn]  +=net[nn].cell[n].v;
//20191025	    ///poles4.dat LPC
//20191025	    for(j=0;j<k;j++) fprintf(fp4,"%e %e ",z[2*j],z[2*j+1]);
//20191025	    //	    for(j=0;j<k;j++) fprintf(fp4,"%e %e ",log(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]));
//20191025	    //	    for(j=0;j<k;j++) fprintf(fp4,"%e %e ",z[2*j],z[2*j+1]);
//20191025	    //	    for(j=0;j<k;j++) fprintf(fp4,"%e %e ",0.5*log(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),log(0.01+fabs(atan2(z[2*j+1],z[2*j]))));
//20191025	    fprintf(fp4,"%e\n",net[nn].cell[n].v);
//20191025
//20191025	    for(j=0;j<k;j++) fprintf(fp7,"%e ",net[nn].cell[n].am.M[0][j]);
//20191025	    fprintf(fp7,"%e\n",net[nn].cell[n].v);
//20191025	    //poles6.dat LSP
//20191025	    sign=1;P[k1]=-1.;P[0]=-sign;
//20191025	    for(j=0;j<k;j++) P[j+1]=A[j]+sign*A[k-1-j];
//20191025	    sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025	    for(j=0;j<k;j++) Q[j+1]=A[j]+sign*A[k-1-j];
//20191025	    for(j=0;j<=2*k1;j++) PQ[j]=0;
//20191025	    for(j=0;j<=k1;j++){
//20191025	      for(i=0;i<=k1;i++){
//20191025		PQ[j+i] +=P[i]*Q[j];
//20191025	      }
//20191025	    }
//20191025	    gsl_poly_complex_solve (PQ, 2*k1+1, wk2, zn);
//20191025	    for(j=0;j<2*k1;j++){
//20191025	      fprintf (fp6,"%e %e ", zn[2*j], zn[2*j+1]);
//20191025	    }
//20191025	    fprintf (fp6,"%e \n",net[nn].cell[n].v);
//20191025	  }
//20191025	}
//20191025      }
//20191025
//20191025      gsl_poly_complex_solve (&An[nn*k1], k+1, wk, &zn[nn*maxk]);
//20191025      net_save(net,"tmp/shrink2.net");
//20191025    }
//20191025    free(PQ);
//20191025    free(P);
//20191025    free(Q);
//20191025//#elif SELECTMETHOD == 1
//20191025//    for(nn=0;nn<nEns;nn++){
//20191025//      vn[nn]=0;//for(j=0;j<=k;j++) An[nn*k1+j]=0;
//20191025//      int nM0max=0,nM0min=0;
//20191025//      for(n=0;n<net[nn].n_cells;n++){	//save the original poles of M
//20191025//	for(j=0;j<k;j++) A[j]=net[nn].cell[n].am.M[0][k-1-j];
//20191025//	gsl_poly_complex_solve (A, k+1, wk, z);
//20191025//	for(j=0;j<k;j++){
//20191025//	  fprintf (fp0,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025//		   z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),net[nn].cell[n].v,j,k,nn);
//20191025//	}
//20191025//      }
//20191025//      for(n=0;n<net[nn].n_cells;n++){//search the 
//20191025//	if(net[nn].cell[n].v*net[nn].vmax > k ){
//20191025//	  if(fabs(net[nn].cell[n].am.M[0][k])<fabs(net[nn].cell[nM0min].am.M[0][k])) nM0min=n;
//20191025//	  if(fabs(net[nn].cell[n].am.M[0][k])>fabs(net[nn].cell[nM0max].am.M[0][k])) nM0max=n;
//20191025//	}
//20191025//      }
//20191025//      fprintf(stderr,"check nM0max=%f, v=%.0f=(%e)*(%d)\n",net[nn].cell[nM0max].am.M[0][k],net[nn].cell[nM0max].v*net[nn].vmax,net[nn].cell[nM0max].v,net[nn].vmax);
//20191025//      fprintf(stderr,"check nM0min=%f, v=%.0f=(%e)*(%d)\n",net[nn].cell[nM0min].am.M[0][k],net[nn].cell[nM0min].v*net[nn].vmax,net[nn].cell[nM0min].v,net[nn].vmax);
//20191025//      for(j=0;j<k;j++) A[j]=net[nn].cell[nM0max].am.M[0][k-1-j];//use nM0max for robustness to noise?
//20191025//      //      for(j=0;j<k;j++) A[j]=net[nn].cell[nM0min].am.M[0][k-1-j];//use nM0max for robustness to noise?
//20191025//      vn[nn]  ++;
//20191025//      for(j=0;j<=k;j++) An[nn  *k1+j]+=A[j];
//20191025//      gsl_poly_complex_solve (&An[nn*k1], k+1, wk, &zn[nn*maxk]);
//20191025////#ifdef NEWSELECTMETHOD
//20191025////      {
//20191025////	FILE *fp=fopen("tmp/M1.dat","w");
//20191025////	for(j=0;j<k-1;j++) fprintf(fp,"%e\n",net[nn].cell[nM0max].am.M[0][j]);
//20191025////	//	for(j=0;j<k-1;j++) fprintf(fp,"%e\n",net[nn].cell[nM0min].am.M[0][j]);
//20191025////	fclose(fp);
//20191025////      }
//20191025////#endif
//20191025//    }
//20191025//#eif SELECTMETHOD == 0 //original
//20191025//    //    gsl_poly_complex_workspace_free(wk2);
//20191025//    for(nn=0;nn<nEns;nn++) for(n=0;n<net[nn].n_cells;n++) net[nn].cell[n].v=v_backup[n*nEns+nn];
//20191025//    for(nn=0;nn<nEns;nn++){
//20191025//      vn[nn]=0;for(j=0;j<=k;j++) An[nn*k1+j]=0;
//20191025//      //      if(nn==3){
//20191025//      //	fprintf(stderr,"check\n");
//20191025//      //      }
//20191025//      for(n=0;n<net[nn].n_cells;n++){
//20191025//	//poles of M
//20191025//	if(net[nn].cell[n].v>0){
//20191025//	  for(j=0;j<k;j++) A[j]=net[nn].cell[n].am.M[0][k-1-j];
//20191025//	  gsl_poly_complex_solve (A, k+1, wk, z);
//20191025//	  for(j=0;j<k;j++){
//20191025//	    fprintf (fp0,"%+e %+e %+e %+e %+e %d %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025//		     z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),net[nn].cell[n].v,j,k,nn);
//20191025//	    fprintf(fp4,"%e %e ",z[2*j],z[2*j+1]);
//20191025//	  }
//20191025//	  fprintf(fp4,"\n");
//20191025//	  {
//20191025//	    int NG=0;
//20191025//	    for(j=0;j<k;j++){
//20191025//	      //	    if(fabs(net[nn].cell[n].am.M[0][k])>1e-5) {NG=1; break;}//constant term 
//20191025//	      double rz2=z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1];
//20191025//	      //	    if(rz2>=1 || rz2<0.5) {NG=1; break;  }
//20191025//	      //	  if(rz2<0.5) {NG=1; break;}
//20191025//	      if(rz2>1.01) {NG=1; break;}
//20191025//	      //	    if(fabs(z[2*j+1])<ZERO) {NG=1; break;}
//20191025//	    }
//20191025//	    //	  if(nEns==1 || NG==0){
//20191025//	    if(net[nn].n_cells==1 || NG==0){
//20191025//	      //	  if(1==1 || nEns==1 || NG==0){
//20191025//	      //	  if(nEns==1 || NG==0){
//20191025//	      //	  if(1==1 || nEns==1 || NG==0){
//20191025//	      for(j=0;j<=k;j++){
//20191025//		//	      An[nn  *k1+j]+=A[j]*net[nn].cell[n].E;
//20191025//		An[nn  *k1+j]+=A[j]*net[nn].cell[n].v;
//20191025//	      }
//20191025//	      //	    vn[nn]  +=net[nn].cell[n].E;
//20191025//	      vn[nn]  +=net[nn].cell[n].v;
//20191025//	    }
//20191025//	  }
//20191025//	}
//20191025//      }
//20191025//      //      if(fabs(An[nn*k1])<ZERO){
//20191025//      //	fprintf(stderr,"Error nn=%d\n",nn);
//20191025//      //      }
//20191025//      gsl_poly_complex_solve (&An[nn*k1], k+1, wk, &zn[nn*maxk]);
//20191025//    }
//20191025#endif
//20191025    fclose(fp0);
//20191025    //////////////////////////////
//20191025    //mean poles zn[nEns]
//20191025    vn[nEns]=0;for(j=0;j<=k;j++) An[nEns*k1+j]=0;// mean An
//20191025    for(nn=0;nn<nEns;nn++){
//20191025      for(j=0;j<k1;j++) An[nEns*k1+j] += An[nn*k1+j];
//20191025      vn[nEns]+=vn[nn];
//20191025    }
//20191025    for(nn=0;nn<=nEns;nn++){for(j=0;j<k1;j++) An[nn*k1+j]/=vn[nn];}
//20191025    //    M0mean/=vn[nEns];
//20191025    gsl_poly_complex_solve (&An[nEns*k1], k+1, wk, &zn[nEns*maxk]);//poles of mean An
//20191025    //    gsl_poly_complex_workspace_free (wk);
//20191025    
//20191025    //order w.r.t. rho and angle
//20191025    for(nn=0;nn<=nEns;nn++){
//20191025      double *znn=&zn[nn*maxk];
//20191025      double *rhon=&rho_zn[nn*maxk];
//20191025      double *thetan=&theta_zn[nn*maxk];
//20191025      int *in=&i_zn[nn*maxk];
//20191025      int itmp;
//20191025      for(j=0;j<k;j++){
//20191025	rhon[j]=sqrt(square(znn[2*j])+square(znn[2*j+1]));
//20191025	//	thetan[j]=atan2(znn[2*j+1],znn[2*j]);
//20191025	//	thetan[j]=fabs(atan2(znn[2*j+1],znn[2*j]));
//20191025	thetan[j]=(atan2(znn[2*j+1],znn[2*j]));
//20191025      }
//20191025      for(j=0;j<k;j++) in[j]=j;
//20191025      for(j=0;j<k-1;j++){
//20191025      	for(i=j+1;i<k;i++){
//20191025      	  if(thetan[in[j]]>thetan[in[i]] || //sort from smallest theta
//20191025      	     (fabs(thetan[in[j]]-thetan[in[i]])<ZERO && rhon[in[j]]> rhon[in[i]])){//sort from smallest rho
//20191025      	    itmp=in[j];
//20191025      	    in[j]=in[i];
//20191025      	    in[i]=itmp;
//20191025      	  }
//20191025      	}
//20191025      }
//20191025    }
//20191025    //nearest poles
//20191025    //    double *thetaEns=&theta_zn[nEns*k];
//20191025    double *znEns=&zn[nEns*maxk];
//20191025    double dd,ddmin;
//20191025    int ii,jj;
//20191025    for(nn=0;nn<=nEns;nn++){//search nearest poles
//20191025      double *znn=&zn[nn*maxk];
//20191025      for(ii=0;ii<k;ii++){
//20191025	c_zn[nn*maxk+ii]=-1;
//20191025	i=i_zn[nEns*maxk+ii];
//20191025	ddmin=1e20;
//20191025	for(j=0;j<k;j++){
//20191025	  dd=square(znn[2*j]-znEns[2*i])+square(znn[2*j+1]-znEns[2*i+1]);
//20191025	  if(ddmin>dd){
//20191025	    ddmin=dd;
//20191025	    c_zn[nn*maxk+ii]=j;
//20191025	  }
//20191025	}
//20191025      }
//20191025    }
//20191025    
//20191025    for(ii=0;ii<k;ii++){
//20191025      v_zn[ii]=0;
//20191025      i=i_zn[nEns*maxk+ii];
//20191025      for(nn=0;nn<nEns;nn++){//mean of poles
//20191025	double *znn=&zn[nn*maxk];
//20191025	j=c_zn[nn*maxk+ii];
//20191025	dd=square(znn[2*j]-znEns[2*i])+square(znn[2*j+1]-znEns[2*i+1]);
//20191025	v_zn[ii]+=(dd*dd);
//20191025      }
//20191025      v_zn[ii]/=nEns;//mean distance
//20191025    }
//20191025    
//20191025    for(nn=0;nn<nEns;nn++){
//20191025      for(j=0;j<k;j++){
//20191025	ii=-1;for(i=0;i<k;i++) if(j==c_zn[nn*maxk+i]) {ii=i;break;}
//20191025	jj=-1;for(i=0;i<k;i++) if(j==i_zn[nn*maxk+i]) {jj=i;break;}
//20191025	fprintf (fp,"%+e %+e %+e %+e %e %d %d %d %d #Re(z) Im(z) |z| Arg(z) v i c j nn\n", 
//20191025		 zn[nn*maxk+2*j],zn[nn*maxk+2*j+1],rho_zn[nn*maxk+j],theta_zn[nn*maxk+j],v_zn[ii],jj,ii,j,nn);
//20191025      }
//20191025    }
//20191025
//20191025//    {
//20191025//      FILE *fp4m=fopen("./tmp/poles4m.dat","w");//original poles, poles in each line
//20191025//      for(nn=0;nn<=nEns;nn++){
//20191025//	double *znn=&zn[nn*maxk];
//20191025//	for(ii=0;ii<k;ii++){
//20191025//	  j=c_zn[nn*maxk+ii];
//20191025//	  fprintf (fp4m,"%g %g ", znn[2*j], znn[2*j+1]);
//20191025//	}
//20191025//	fprintf(fp4m,"%g\n",1.0);
//20191025//      }
//20191025//      fclose(fp4m);
//20191025//    }
//20191025
//20191025    {
//20191025      for(nn=0;nn<=nEns;nn++){
//20191025	double *znn=&zn[nn*maxk];
//20191025	for(ii=0;ii<k;ii++){
//20191025	  j=c_zn[nn*maxk+ii];
//20191025	  fprintf (fp2,"%+e %+e %e %d ", znn[2*j], znn[2*j+1],sqrt(v_zn[ii]),ii);
//20191025	}
//20191025	fprintf (fp2,"%d \n",nn);
//20191025      }
//20191025    }
//20191025//    {
//20191025//      FILE *fp=fopen("tmp/poles4.dat","w");
//20191025//      for(nn=0;nn<nEns;nn++){
//20191025//	double *znn=&zn[nn*maxk];
//20191025//	for(j=0;j<k;j++){
//20191025//	  fprintf (fp,"%+e %+e ", znn[2*j], znn[2*j+1]);
//20191025//	}
//20191025//	fprintf (fp,"\n");
//20191025//      }
//20191025//      fclose(fp);
//20191025//    }
//20191025    {
//20191025      nn=nEns;
//20191025      {
//20191025	double *znn=&zn[nn*maxk];
//20191025	for(j=0;j<k;j++){
//20191025	  fprintf (fp5,"%e %e ", znn[2*j], znn[2*j+1]);
//20191025	}
//20191025	//	fprintf(fp5,"%e\n",net[nn].cell[n].v);
//20191025	fprintf(fp5,"%e\n",1.0);
//20191025      }
//20191025      fclose(fp5);
//20191025    }
//20191025
//20191025    {//remove poles with only real part
//20191025      FLOAT *M0=(FLOAT*)malloc(sizeof(double)*k1);
//20191025      FLOAT *M1=(FLOAT*)malloc(sizeof(double)*k1);
//20191025      M1[0]=1;for(i=1;i<k1;i++) M1[i]=0;
//20191025      for(j=0;j<k;j++){
//20191025	for(i=0;i<k1;i++) M0[i]=M1[i];
//20191025	if(fabs(znEns[2*j+1])<1e-20){
//20191025	  //	  for(i=0;i<k1-1;i++) M1[i+1]-=M0[i]*znEns[2*j];
//20191025	  continue; //remove poles with only real part
//20191025	}
//20191025	else if(znEns[2*j+1]>0){
//20191025	  for(i=0;i<k1-1;i++) M1[i+1]-=(M0[i]*2.0*znEns[2*j]);
//20191025	  for(i=0;i<k1-2;i++) M1[i+2]+=(M0[i]*(square(znEns[2*j])+square(znEns[2*j+1])));
//20191025//	  double re=cos(atan2(znEns[2*j+1],znEns[2*j]));//magnify to 1
//20191025//	  for(i=0;i<k1-1;i++) M1[i+1]-=(M0[i]*2.0*re);
//20191025//	  for(i=0;i<k1-2;i++) M1[i+2]+=(M0[i]);
//20191025	}
//20191025      }
//20191025
//20191025      for(j=0;j<k;j++) net[0].cell[0].am.M[0][j]=-M1[j+1];
//20191025      net[0].cell[0].am.M[0][k]=0;
//20191025      //      net[0].cell[0].am.M[0][k]=M0mean;
//20191025      //      net[0].cell[0].am.M[0][j]=0;//??
//20191025      net[0].nEns=1;
//20191025      net[0].n_cells=1;
//20191025      net_save(net,"tmp/shrink1.net");
//20191025      //#ifndef NEWSELECTMETHOD
//20191025      {
//20191025	FILE *fp=fopen("tmp/M1.dat","w");
//20191025	for(j=0;j<k;j++) fprintf(fp,"%e\n",net[0].cell[0].am.M[0][j]);
//20191025	fclose(fp);
//20191025      }
//20191025      //#endif
//20191025      {
//20191025	double *thetanEns=&theta_zn[nEns*maxk];
//20191025	int *inEns=&i_zn[nEns*maxk];
//20191025	FILE *fp=fopen("tmp/theta1.dat","w");
//20191025	for(j=0;j<k;j++) fprintf(fp,"%e\n",thetanEns[inEns[j]]);
//20191025	fclose(fp);
//20191025      }
//20191025
//20191025      fprintf(stdout,"Net is shrinked to single unit and saved in shrink1.net.\n");
//20191025      free(M0);
//20191025      free(M1);
//20191025    }
//20191025
//20191025    free(A);
//20191025    free(z);
//20191025  }
//20191025  /// LSP
//20191025  {
//20191025    //LSP:Line Spectral Pairs
//20191025    //    if(1==0){
//20191025    //      gsl_poly_complex_workspace *wk1= gsl_poly_complex_workspace_alloc (k1+1);//for LSP
//20191025    //      sign=1;P[k1]=-1.;P[0]=-sign;
//20191025    //      for(j=0;j<k;j++) P[j+1]=As[j]+sign*As[k-1-j];   //  for(i=0;i<=k1;i++) printf("P[%d]=%e\n",i,P[i]);
//20191025    //      //    w= gsl_poly_complex_workspace_alloc (k2);//for poles
//20191025    //      gsl_poly_complex_solve (P, k1+1, wk1, z);
//20191025    //      for(j=0;j<k1;j++){
//20191025    //	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025    //		 z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025    //      }
//20191025    //      sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025    //      for(j=0;j<k;j++) Q[j+1]=As[j]+sign*As[k-1-j]; //  for(i=0;i<=k1;i++) printf("Q[%d]=%e\n",i,P[i]);
//20191025    //      //    w= gsl_poly_complex_workspace_alloc (k2);//for poles
//20191025    //      gsl_poly_complex_solve (Q, k1+1, wk1, z);
//20191025    //      gsl_poly_complex_workspace_free (wk1);
//20191025    //      for(j=0;j<k1;j++){
//20191025    //	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025    //		 z[2*j], z[2*j+1],sqrt(z[2*j]*z[2*j]+z[2*j+1]*z[2*j+1]),atan2(z[2*j+1],z[2*j]),j,k);
//20191025    //      }
//20191025    //      gsl_poly_complex_workspace_free (wk1);
//20191025    //      free(P);
//20191025    //      free(Q);
//20191025    //    }
//20191025    
//20191025    {//k3=2*k1+1
//20191025      double *P=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025      double *Q=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025      double *PQ=(double*)malloc(sizeof(double)*(2*k1+1));//for LSP
//20191025      double *As=&An[nEns*k1];
//20191025      int sign;
//20191025      //      gsl_poly_complex_workspace *wk2=gsl_poly_complex_workspace_alloc (2*k1+1);//for LSP
//20191025      sign=1;P[k1]=-1.;P[0]=-sign;
//20191025      for(j=0;j<k;j++) P[j+1]=As[j]+sign*As[k-1-j];
//20191025      sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025      for(j=0;j<k;j++) Q[j+1]=As[j]+sign*As[k-1-j];
//20191025      for(j=0;j<=2*k1;j++) PQ[j]=0;
//20191025      for(j=0;j<=k1;j++){
//20191025	for(i=0;i<=k1;i++){
//20191025	  PQ[j+i] +=P[i]*Q[j];
//20191025	}
//20191025      }
//20191025      //      for(i=0;i<2*k1;i++) printf("PQ[%d]=%e\n",i,PQ[i]);
//20191025      gsl_poly_complex_solve (PQ, 2*k1+1, wk2, zn);
//20191025      //      gsl_poly_complex_workspace_free (w);
//20191025      //    gsl_poly_complex_solve (Q, k1+1, w2, z);
//20191025      for(j=0;j<2*k1;j++){
//20191025	fprintf (fp3,"%+e %+e %+e %+e %d %d #Re(z) Im(z) |z| Arg(z) v i k net\n", 
//20191025		 zn[2*j], zn[2*j+1],sqrt(zn[2*j]*zn[2*j]+zn[2*j+1]*zn[2*j+1]),atan2(zn[2*j+1],zn[2*j]),j,k);
//20191025      }
//20191025      //      gsl_poly_complex_workspace_free(wk2);
//20191025      free(PQ);
//20191025      free(P);
//20191025      free(Q);
//20191025    }
//20191025
//20191025//    {//k3=2*k1+1
//20191025//      double *P=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025//      double *Q=(double*)malloc(sizeof(double)*(k1+1));//for LSP
//20191025//      FILE *fp6=fopen("./tmp/poles6.dat","w");
//20191025//      double *PQ=(double*)malloc(sizeof(double)*(2*k1+1));//for LSP
//20191025//      int sign;
//20191025//      gsl_poly_complex_workspace *wk2=gsl_poly_complex_workspace_alloc (2*k1+1);//for LSP
//20191025//      for(nn=0;nn<nEns;nn++){
//20191025//	double *As=&An[nn*k1];
//20191025//	sign=1;P[k1]=-1.;P[0]=-sign;
//20191025//	for(j=0;j<k;j++) P[j+1]=As[j]+sign*As[k-1-j];
//20191025//	sign=-1;Q[k1]=-1.;Q[0]=-sign;
//20191025//	for(j=0;j<k;j++) Q[j+1]=As[j]+sign*As[k-1-j];
//20191025//	for(j=0;j<=2*k1;j++) PQ[j]=0;
//20191025//	for(j=0;j<=k1;j++){
//20191025//	  for(i=0;i<=k1;i++){
//20191025//	    PQ[j+i] +=P[i]*Q[j];
//20191025//	  }
//20191025//	}
//20191025//	gsl_poly_complex_solve (PQ, 2*k1+1, wk2, zn);
//20191025//	for(j=0;j<2*k1;j++){
//20191025//	  fprintf (fp6,"%e %e ", zn[2*j], zn[2*j+1]);
//20191025//	}
//20191025//	fprintf (fp6,"%e \n",net[nn].cell[n].v);
//20191025//      }
//20191025//      gsl_poly_complex_workspace_free(wk2);
//20191025//      fclose(fp6);
//20191025//      free(PQ);
//20191025//      free(P);
//20191025//      free(Q);
//20191025//    }
//20191025  }
//20191025
//20191025  {//shrink the net
//20191025    FLOAT vsum=0;
//20191025    FLOAT *M0=(FLOAT*)malloc(sizeof(FLOAT)*(k+1));
//20191025    for(j=0;j<=k;j++) M0[j]=0;
//20191025    for(nn=0;nn<nEns;nn++){
//20191025      for(n=0;n<net[nn].n_cells;n++){
//20191025	//	if(fabs(net[nn].cell[n].am.M[0][k])>1e-5) continue;
//20191025	for(j=0;j<=k;j++) M0[j]+=(net[nn].cell[n].am.M[0][j]*net[nn].cell[n].v);
//20191025	vsum+=net[nn].cell[n].v;
//20191025      }
//20191025    }
//20191025    for(j=0;j<=k;j++) net[0].cell[0].am.M[0][j]=M0[j]/vsum;
//20191025    net[0].nEns=1;
//20191025    net[0].n_cells=1;
//20191025    net_save(net,"tmp/shrink.net");
//20191025    fprintf(stdout,"Net is shrinked to single unit and saved in shrink.net.\n");
//20191025    free(M0);
//20191025  }
//20191025
//20191025  gsl_poly_complex_workspace_free (wk);
//20191025  gsl_poly_complex_workspace_free(wk2);
//20191025  free(vn);
//20191025  free(An);
//20191025  free(zn);
//20191025  free(theta_zn);
//20191025  free(rho_zn);
//20191025  free(i_zn);
//20191025  free(c_zn);
//20191025  free(v_zn);
//20191025  fclose(fp3);
//20191025  fclose(fp2);
//20191025  fclose(fp);
//20191025  fclose(fp4);
//20191025  fclose(fp6);
//20191025  fclose(fp7);
//20191025  //  fprintf(stdout,"Poles of M are saved in '%s','%s','%s'.",fname,fname2,fname3);
//20191025  fprintf(stdout,"Poles of M are saved in './tmp/poles?.dat'.");
//20191025}
/*====================================================================*
 * modify_w_batch
 * 荷重ベクトルの更新
 *====================================================================*/
void modify_w_batch(NET *net, FLOAT **x, FLOAT *y,int n_train) {
  char name[32] = "modify_w_batch";
  int *s = NULL;
  FLOAT *d2 = NULL;
  //,*dw,dwNorm;
  FLOAT err_i, err_xi, delta_w_ic;
  //  FLOAT delta_w_x;
  FLOAT y_hat_i, y_hat_xi;
  FLOAT alpha;
  FLOAT width = net->width;
  //  FLOAT beta = 1.0e-1;
  int n_cells = net->n_cells, n_channels = net->k, n_compare = net->n_compare;
  int i,j, k,t, xi,sxi;
  //  int nnn=0;
  FLOAT dwNorm=0;
  VORONOI *V=net->V;
#ifdef DEBUG
  if(debug==1){
    printf("start modify_w_batch\n");
    for(i=0;i<net->n_cells;i++){
      printf("i%d w%g M%g %g\n",i,net->cell[i].w[0],net->cell[i].am.M[0][0],net->cell[i].am.M[0][1]);
    }
  }
#endif
  s = (int *)malloc(sizeof(int)*n_cells);
  d2 = (FLOAT *)malloc(sizeof(FLOAT)*n_cells);
  if (s == NULL||d2==NULL) {
    printf("error: memory allocation error #?. (%s)\n", name);
    exit(-1);
  }
  for(i=0;i<net->n_cells;i++)for(k=0; k<n_channels; k++) net->cell[i].dw[k] =0;

  for(i=0;i<net->n_cells;i++){
    for (j=0; j<V->i2v0[i]; j++) {
      //for (j=0; j<V->i2v[i]; j++) {
      t=V->ij2t[i][j];
      init_sort_weights(net, s, d2, x[t]);
      for (xi=1; xi<n_compare; xi++) {
      //      for (xi=0; xi<n_compare; xi++) {
      //      for (xi=0; xi<net->n_cells; xi++) {
	sort_weights(net, s, d2, xi, xi+1);
	sxi=s[xi];
	if(sxi==i) continue;
	// 境界の中？
	if(in_window(net->cell[sxi].w, net->cell[i].w,x[t],width,n_channels) > 0) {
	  //	  printf("%e %e\n",x[t][0],x[t][1]);//for debug
	  y_hat_i = y_hat_xi = delta_w_ic = 0.0;
	  for (k=0; k<=n_channels; k++) {
	    y_hat_i  += net->cell[i].am.M[0][k]*x[t][k];
	    y_hat_xi += net->cell[sxi].am.M[0][k]*x[t][k];
	  }
	  for (k=0;k<n_channels; k++) 
	    // delta_w_ic += (FLOAT)square(x[t][k]-net->cell[i].w[k]);//031205
	    //	  delta_w_ic = (FLOAT)sqrt(delta_w_ic)*2.;//031205
	    delta_w_ic += (FLOAT)square(net->cell[sxi].w[k]-net->cell[i].w[k]);//orig
	  delta_w_ic = (FLOAT)sqrt(delta_w_ic);//orig
	  
	  err_xi = (y_hat_xi-y[t])/net->ywidth;
	  err_i  = (y_hat_i -y[t])/net->ywidth;
	  if (delta_w_ic <1e-20) continue;
	  alpha = ((FLOAT)square(err_xi)-(FLOAT)square(err_i))/delta_w_ic;
//031204 from 
	  for (k=0; k<n_channels; k++) {
	    net->cell[i].dw[k]   += alpha*(x[t][k]-net->cell[i].w[k]);
	    net->cell[sxi].dw[k] -= alpha*(x[t][k]-net->cell[sxi].w[k]);
	  }
	  break;
	}
      }
    }
    for(k=0;k<(n_channels);k++){
#ifdef dwNormOrig
      dwNorm+= square(net->cell[i].dw[k]);
#else
      if(dwNorm<fabs(net->cell[i].dw[k])) dwNorm=fabs(net->cell[i].dw[k]); //max?
#endif
    }
  }
#ifdef dwNormOrig
  dwNorm=sqrt(dwNorm);
#endif
  
//  {////////030709 近づき過ぎは離す
//    FLOAT dw01,dwk;
//    int s1;
//    for(i=0;i<net->n_cells;i++){
//	if(net->cell[i].v<ZERO){
//	  printf("***v(%d)=%e",i,net->cell[i].v);
//	  init_sort_weights(net,s,d2,net->cell[i].w);
//	  //      sort_weights(net, s, d2, 0, 2);
//	  sort_weights(net, s, d2, 1, 2);
//	  if(s[0]!=i){
//	    s1=s[0];
//	    printf("hen i,s0,s1=%d,%d,%d,v(i)=%e,ds0=%e,ds1=%e!!!!!!!\n",
//		   i,s[0],s[1],net->cell[i].v,d2[s[0]],d2[s[1]]);
//	  }
//	  else{
//	    s1=s[1];
//	  }
//	  if((dw01=sqrt(d2[s1]))<net->xwidth*0.005){
//	    if(dw01<ZERO){
//	      for(k=0;k<(n_channels);k++){
//		dwk=0.005*net->xwidth*(2.*(rand()/(RAND_MAX+1.0))-1.);
//		net->cell[i].dw[k] += dwk;
//		//	    net->cell[s1].dw[k] += dwk;
//	      }
//	    }
//	    else{
//	      for(k=0;k<(n_channels);k++){
//		dwk=0.1*(net->cell[i].dw[k]-net->cell[s1].dw[k]);
//		net->cell[i].dw[k] += dwk;
//		//	    net->cell[s1].dw[k] += dwk;
//	      }
//	    }
//	  }
//	}	      
//    }
//  }
  ////////
//  alpha=0.10*(1./GlobalTime-1./GlobalTimeMax);//good 
//  alpha=0.05*(1./GlobalTime-1./GlobalTimeMax)+0.01;//bad
//  alpha=0.10*(1.-GlobalTime/GlobalTimeMax)+0.001;//bad
//  alpha=0.05*(1./GlobalTime-1./GlobalTimeMax)+0.001;//good?
//  alpha=0.05*(1./GlobalTime-1./GlobalTimeMax)+1e-4;//so so
  //  alpha=0.05*(1./GlobalTime-1./GlobalTimeMax);//so so
  //  alpha=0.05*(1./GlobalTime-1./100);//so so
  //  alpha=0.01*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//
  //  alpha=0.04*(1-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);//soso//100MSE012 1.734e-02 6.19e-02 3.65e-02NMSE06.96e-02MSEtr2.44e-03N500k2w0.1vm3:0
  //  alpha=0.09*(1.-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//soso100MSE012 2.030e-02
  //  alpha=0.05*(1-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);//soso
  //  alpha=0.1*(1-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);//bad

  //alpha *= 0.33*alphat(GlobalTime)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);
  //  alpha=0.05*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//bestforN500
  //  alpha=0.08*(1.-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//pretty
  //  alpha=0.04*(1-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);//soso
  //  alpha=0.04*(1-net->nentropy)*(net->xwidth)/(dwNorm+ZERO);  printf("!!!!!dwNorm=%e\n",dwNorm);//soso
  //  alpha=0.15*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//bestforN100
  //  alpha=LA(net->n_cells,0.15,100,0.05,500)*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//bestforN100
  //  alpha=LA(net->n_cells,0.15,100,0.10,500)*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);printf("!!!!!dwNorm=%e\n",dwNorm);//bestforN100
  //  alpha=LA(net->n_cells,0.15,100,0.10,500)*alphat(GlobalTime)* (net->xwidth)/(dwNorm+ZERO);//good
  //  alpha=LA(net->n_cells,0.15,100,0.10,500)*(net->n_cells*100./n_train)/(1.+GlobalTime/alphaT1)*(net->xwidth)/(dwNorm+ZERO);//good

  //  printf("!!!!!dwNorm=%e\n",dwNorm);
  //  printf("check w \n");//check w 20190115
  if(dwNorm>0){
    alpha=net->gamma0/(1.+GlobalTime/net->Tgamma)*(net->xwidth)/(dwNorm);//good
    //    alpha=net->gamma0/(1.+GlobalTime/net->Tgamma)*(net->xwidth)/(dwNorm+ZERO);//good
    for (i=0; i<n_cells; i++){
      //    printf("[%d]=",i);
      for (k=0; k<(n_channels); k++) {
	net->cell[i].w[k] += (net->cell[i].dw[k]*alpha);
	//      net->cell[i].w[k] = net->cell[i].w[k]+ net->cell[i].dw[k]*alpha;
	//      printf("w%g += %g * %g",net->cell[i].w[k],net->cell[i].dw[k],alpha);
      }
      //    printf("\n");
    }
    net->printf(strx3("%d #modify_w_batch alpha=%.3g max_dw=%.3g\n",net->GlobalTime,alpha,dwNorm));
  }
  free(s); free(d2);
}

//void modify_w_batch2(NET *net, int i, FLOAT **x, FLOAT *y, int n_ivectors) {
//  char name[32] = "modify_w_batch";
//  int *s = NULL;
//  FLOAT *v = NULL;
//  //,*dw,dwNorm;
//  FLOAT err_i, err_xi, delta_w_ic;
//  //  FLOAT delta_w_x;
//  FLOAT y_hat_i, y_hat_xi;
//  FLOAT alpha,alphasum;
//  FLOAT width = net->width;
//  //  FLOAT beta = 1.0e-1;
//  int n_cells = net->n_cells, n_channels = net->k, n_compare = net->n_compare;
//  int j, k, xi,nnn=0,sxi;
//  
//  //  dw= (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
//  //  if (dw == NULL) {
//  //    printf("error: memory allocation error #?. (%s)\n", name);
//  //    exit(-1);
//  //  }
//  //  for (k=0; k<(n_channels); k++) dw[k]=0;
//#define COMPARE
//#ifdef COMPARE
//  s = (int *)malloc(sizeof(int)*n_cells);
//  v = (FLOAT *)malloc(sizeof(FLOAT)*n_cells);
//  if (s == NULL) {
//    printf("error: memory allocation error #?. (%s)\n", name);
//    exit(-1);
//  }
//  if (v == NULL) {
//    printf("error: memory allocation error #?. (%s)\n", name);
//    exit(-1);
//  }
//  alphasum=0;
//  for (j=0; j<n_ivectors; j++) {
//    init_sort_weights(net, s, v, x[j]);
//    for (xi=1; xi<n_compare; xi++) {
//	sxi=s[xi];
//#else
//  for (j=0; j<n_ivectors; j++) {
//    for (sxi=0; sxi<n_cells; sxi++) {
//#endif
//	if(sxi==i) continue;
//	// 境界の中？
//	if (in_window(net->cell[sxi].w, net->cell[i].w, x[j], width, n_channels) > 0) {
//	  y_hat_i = y_hat_xi = delta_w_ic = 0.0;
//	  for (k=0; k<=n_channels; k++) {
//	    y_hat_i += net->cell[i].am.M[0][k]*x[j][k];
//	    y_hat_xi += net->cell[sxi].am.M[0][k]*x[j][k];
//	  }
//	  for (k=0;k<n_channels; k++) 
//	    delta_w_ic += (FLOAT)square(net->cell[sxi].w[k]-net->cell[i].w[k]);
//	  
//	  err_xi = (y_hat_xi-y[j])/net->ywidth;
//	  err_i  = (y_hat_i -y[j])/net->ywidth;
//	  delta_w_ic = (FLOAT)sqrt(delta_w_ic);
//	  if (delta_w_ic <1e-20) continue;
//	  alpha = ((FLOAT)square(err_xi)-(FLOAT)square(err_i))/delta_w_ic;
//	  //	if(alpha<0) alpha=0;else if(alpha>=1) alpha=1;
//	  if(alpha<-0.9) alpha=-0.9;else if(alpha>=0.9) alpha=0.9;
//	  for (k=0; k<n_channels; k++) {
//	    net->cell[i].dw[k] += alpha*(x[j][k]-net->cell[i].w[k]);
//	    net->cell[sxi].dw[k] -= alpha*(x[j][k]-net->cell[sxi].w[k]);
//	  }
//	  nnn++;
//	  //	break;
//	}
//    }
//    //    printf("nnn=%d\n",nnn);
//  }
//  //  free(dw);
//  //  return;
//#ifdef BAKA
//    }
//  }
//#endif
//}
    

//void modify_w_batch0(NET *net, int i, FLOAT **x, FLOAT *y, int n_ivectors) {
//  char name[32] = "modify_w_batch";
//  int *s = NULL;
//  FLOAT *v = NULL;
//  FLOAT err_i, err_xi, delta_w_x, delta_w_ic;
//  FLOAT y_hat_i, y_hat_xi;
//  FLOAT alpha;
//  FLOAT width = net->width;
//  FLOAT beta = 1.0e-1;
//  int n_cells = net->n_cells, n_channels = net->k, n_compare = net->n_compare;
//  int j, k, xi;
//  
//  s = (int *)malloc(sizeof(int)*n_cells);
//  v = (FLOAT *)malloc(sizeof(FLOAT)*n_cells);
//  if (s == NULL) {
//    printf("error: memory allocation error #?. (%s)\n", name);
//    exit(-1);
//  }
//  if (v == NULL) {
//    printf("error: memory allocation error #?. (%s)\n", name);
//    exit(-1);
//  }
//  
//  for (j=0; j<n_ivectors; j++) {
//    init_sort_weights(net, s, v, x[j]);
//    alpha = 0.0;
//    
//    for (xi=1; xi<n_compare; xi++) {
//	  sort_weights(net, s, v, xi, xi+1);
//	  err_i = err_xi = delta_w_x = delta_w_ic = 0.0;
//	  
//	  // 境界の中？
//	  if (in_window(net->cell[s[xi]].w, net->cell[i].w, x[j], width, n_channels) > 0) {
//	    y_hat_i = y_hat_xi = 0.0;
//	    for (k=0; k<=n_channels; k++) {
//	      y_hat_i += net->cell[i].am.M[0][k]*x[j][k];
//	      y_hat_xi += net->cell[s[xi]].am.M[0][k]*x[j][k];
//	    }
//	    for (k=0; k<n_channels; k++)
//	      delta_w_ic += (FLOAT)square(net->cell[s[xi]].w[k]-net->cell[i].w[k]);
//	    err_xi = y[j]-y_hat_xi;
//	    err_i  = y[j]-y_hat_i;
//	    //delta_w_x  += (FLOAT)(x[j][k]-net->cell[xi].w[k]);
//	    delta_w_ic = (FLOAT)sqrt(delta_w_ic);
//	  }
//	  if (delta_w_ic > 0.0)
//	    alpha += ((FLOAT)square(err_xi)-(FLOAT)square(err_i))/delta_w_ic/beta;
//    }
//    for (k=0; k<n_channels; k++) {
//	  net->cell[i].w[k] = (1.0-alpha)*net->cell[i].w[k] + alpha*x[j][k];
//    }
//  }
//
//  /*
//  for (k=0; k<(n_channels+1); k++)
//    printf("w[%2d][%d] = %3.2e, ", i, k, net->cell[i].w[k]);
//  printf("\n");
//  */
//
//  free(v);
//  free(s);
//  return;
//}

/*====================================================================*
 * cond_reinit_batch
 * 再初期化条件を満たすかチェックする
 *====================================================================*/
//int cond_reinit_batch(NET *net, int m, int n) {
//  int ret = 0;
//  
//  ret = ((m != n)
//	   && (net->alpha_NG == 0)
//	   && (net->cell[m].alpha > net->v_ratio*net->alpha_bar)
//	   && (net->cell[n].v < net->v_thresh)
//	   //&& (net->cell[m].alpha > net->v_ratio*net->alpha_hat)
//	   );
//  return(ret);
//}
/*====================================================================*
 * reinit_cell_batch
 * セル（ユニットの）再初期化を行う
 *====================================================================*/
int reinit_cell_batch(NET *net, FLOAT **x,FLOAT *y,int n_train) {

  FLOAT min_length, length;
  //  int n_cells = net->n_cells;
  int n_channels = net->k;
  int xi_j = 0;
  int i,j,k,N_i;
  int rho_i, rho_N_i, temp,reinit=0;
  int n_cells=net->n_cells;
  int *rho = (int *)malloc(sizeof(int)*n_cells);
  int iii=0;//restrict reinit units

#ifdef DEBUG
  if(debug==1){
    printf("start reinit_cell_batch\n");
    for(i=0;i<net->n_cells;i++){
      printf("i%d w%g M%g %g\n",i,net->cell[i].w[0],net->cell[i].am.M[0][0],net->cell[i].am.M[0][1]);
    }
  }
#endif
  for (i=0; i<n_cells; i++) rho[i] = i;
  // Sort via alpha from biggest
  for (i=0; i<n_cells-1; i++) {
    for (j=i+1; j<n_cells; j++) {
      if (net->cell[(int)rho[i]].alpha < net->cell[(int)rho[j]].alpha) {
	temp = rho[j];
	rho[j] = rho[i];
	rho[i] = temp;
      }
    }
  }
//  printf("laptime:reinit_cell_batch1=%5.3f\n",mytimer_lap());
//  printf("!!!!!!!!!!!alpha[%d]=%e\n",1,net->cell[1].alpha);
//  printf("!!!!!!!!!!!alpha[%d]=%e\n",229,net->cell[229].alpha);
  N_i=n_cells-1;
  for (i=0; i<n_cells;) {
    rho_i = rho[i++];
    //    printf("rho[%d]=%d\n",i,rho[i]);
    if(net->cell[rho_i].alpha > net->v_ratio*net->alpha_bar){
      for(;N_i>=0;){
	if(i>=N_i) break;
	rho_N_i = rho[N_i--];
	if(net->cell[rho_N_i].v >= net->v_thresh) continue;
	//	if(net->cell[rho_N_i].v >= net->v_thresh) continue;
	//	if(net->cell[rho_N_i].v >= net->v_thresh+(1.-net->v_thresh)*0.2) continue;
	// for Debug
	//	if(++iii>=net->n_cells*0.10 && iii>=2) break;//040119??
	//	if(++iii>=net->n_cells*0.30 && iii>=2) break;//good for N500
	//	if(++iii>=net->n_cells*LA(net->n_cells,0.1,100,0.3,500) && iii>=2) break;//good for 
	if(++iii>=net->N_reinit_max && iii>=2) break;//good for 
	net->printf(strx9("%d #reinit w[%d]=[%.2g,]ga%.2g>w[%d]=[%.2g,]a%.2ga_%.2g,t%d\n",
		(int)GlobalTime,
		rho_N_i,net->cell[rho_N_i].w[0],
		net->cell[rho_N_i].alpha,
		rho_i, net->cell[rho_i].w[0],
		net->cell[rho_i].alpha,
                         net->alpha_bar,xi_j));
		//		i,
		//		net->cell[rho_N_i].w[1],
		//		net->cell[rho_i].w[1],
//	fprintf(stderr,"%d #reinit w[%d]=[%.2g,%.2g]a%.2g>w[%d]=[%.2g,%.2g]a%.2ga_%.2g,t%d\n",
//		(int)GlobalTime,
//		//		i,
//		rho_N_i,net->cell[rho_N_i].w[0],
//		net->cell[rho_N_i].w[1],
//		net->cell[rho_N_i].alpha,
//		rho_i, net->cell[rho_i].w[0],
//		net->cell[rho_i].w[1],
//		net->cell[rho_i].alpha,
//		net->alpha_bar,xi_j);
	// Initialize w_(rho_N_i)
	// Search the vector nearest to w_(rho_i)
#undef move2x
#define move2x
#ifdef move2x
	{
	  VORONOI *V=net->V;
	  int t;
	  min_length = 1e+20;
	  for (j=0; j<V->i2v0[rho_i]; j++) {
	    t=V->ij2t[rho_i][j];
	    length = distance2(net->cell[rho_i].w, x[t], n_channels);
	    if (length>1e-10 && length < min_length) {
	      min_length = length;
	      xi_j = t;
	    }
            //if(t>70) printf("t%d length=%g=|%g-%g|,min=%g %d\n",t,length,net->cell[rho_i].w[0],x[t][0],min_length,xi_j);
	  }
	  //	printf("***********xi_j=%d\n",xi_j);
	  //for (k=0; k<(n_channels); k++){//reinit to w[rho_N_i] w[rho_i] 
	  //net->cell[rho_N_i].w[k] = (net->cell[rho_i].w[k]+x[xi_j][k])/2.;
	  //}
          //          printf("#check reinit: w[%d]=%.2g -> w[%d]=%.2g x=%.2g",rho_N_i,net->cell[rho_N_i].w[0],rho_i,net->cell[rho_i].w[0],x[xi_j][0]);
	  for (k=0; k<(n_channels); k++){//reinit to w[rho_N_i] w[rho_i] 
	    net->cell[rho_N_i].w[k]= 
	      net->cell[rho_i].w[k]+1.1*(x[xi_j][k]-net->cell[rho_i].w[k]);//good theta_r=1.1
	    //net->cell[rho_i].w[k]+1.9*(x[xi_j][k]-net->cell[rho_i].w[k]);//good theta_r=1.9
            //	      net->cell[rho_i].w[k]+0.9*(x[xi_j][k]-net->cell[rho_i].w[k]);//good theta_r=0.9
	      //net->cell[rho_i].w[k]+1.61803398875*(x[xi_j][k]-net->cell[rho_i].w[k]);//good
	    //	      net->cell[rho_i].w[k]+1.5*(x[xi_j][k]-net->cell[rho_i].w[k]);//good
	    //net->cell[rho_i].w[k]+3.*(x[xi_j][k]-net->cell[rho_i].w[k]);
	    //net->cell[rho_i].w[k]+4.*(x[xi_j][k]-net->cell[rho_i].w[k]);
	    //net->cell[rho_i].w[k]+2.*(x[xi_j][k]-net->cell[rho_i].w[k]);
	  }
	  //	  net->cell[rho_N_i].w[n_channels]=1;
	}
#else
	min_length = 1e+20;
	for (j=0; j<net->n_cells; j++) {
	  if(j==rho_i) continue;
	  length = distance2(net->cell[rho_i].w, net->cell[j].w, n_channels);
	  if (length < min_length) {
	    min_length = length;
	    xi_j = j;
	  }
	}

	for (k=0; k<(n_channels); k++){//reinit to w[rho_N_i] w[rho_i] 
	  net->cell[rho_N_i].w[k] = 
	    (net->cell[rho_i].w[k]+net->cell[xi_j].w[k])/2.;
	}
#endif
#define REINITMODE 2
#if REINITMODE == 0
	for (k=0; k<=n_channels; k++) {//replace AM[rho_N_i] by AM[rho_i]
	  net->cell[rho_N_i].am.M[0][k] = net->cell[rho_i].am.M[0][k];//kuro
	  for (j=0; j<=(n_channels); j++){//kuro
	    net->cell[rho_N_i].am.P[k][j] = net->cell[rho_i].am.P[k][j];//kuro
	  }
	}
#elif REINITMODE == 1
	//	modify_M(net, rho_N_i,x[xi_j],y[xi_j]);
	init_AMdata(&(net->cell[rho_N_i].am));
//	modify_M(net, rho_N_i,x[xi_j],y[xi_j]);
	init_AMdata(&(net->cell[rho_i].am));
	net->cell[rho_N_i].v0=0;
	net->cell[rho_i].v0=0;
#elif REINITMODE == 2
	for (k=0; k<=n_channels; k++) {//replace AM[rho_N_i] by AM[rho_i]
	  net->cell[rho_N_i].am.M[0][k] = net->cell[rho_i].am.M[0][k];//kuro
	  for (j=0; j<=(n_channels); j++){//kuro
	    if(k==j){
	      net->cell[rho_N_i].am.P[k][j] = 1e4;
	      net->cell[rho_i].am.P[k][j] = 1e4;
	    }
	    else{
	      net->cell[rho_N_i].am.P[k][j] = 0;
	      net->cell[rho_i].am.P[k][j] = 0;
	    }
	  }
	}
	net->cell[rho_N_i].v0=0;
	net->cell[rho_i].v0=0;
#endif
	reinit=1;
	//	N_i--;
	break;
      }
    }
  }
#ifdef DEBUG
  if(debug==1){
    printf("finish reinit_cell_batch\n");
    for(i=0;i<net->n_cells;i++){
      printf("i%d w%g M%g %g\n",i,net->cell[i].w[0],net->cell[i].am.M[0][0],net->cell[i].am.M[0][1]);
    }
  }
#endif
  free(rho);
  ReinitTime=GlobalTime;
  return(reinit);
}
//void reinit_cell_batch2(NET *net, int rho_j, int rho_N_j, 
//FLOAT **x,FLOAT *y,int n_ivectors) {
//  FLOAT min_length, length;
//  int n_cells = net->n_cells;
//  int n_channels = net->k;
//  int xi_j = 0;
//  int i = 0, k = 0;
//
//  // for Debug
//  printf(">%3.0f)reinit cell[%3d]w%3.2f,%3.2falpha%3.2e>[%3d]w%3.2f,%3.2falpha%3.2ealpha_%3.2e,[%d]w\n",GlobalTime,
//	   rho_N_j,net->cell[rho_N_j].w[0],net->cell[rho_N_j].w[1],net->cell[rho_N_j].alpha,
//	   rho_j, net->cell[rho_j].w[0],net->cell[rho_j].w[1],net->cell[rho_j].alpha,
//	   net->alpha_bar,xi_j);
//  // Initialize w_(rho_N_j)
//
//  // Search the vector nearest to w_(rho_j)
//  min_length = 1e+20;
//  for (i=0; i<n_cells; i++) {
//    if (i == rho_j) continue;
//    length = calc_length(net->cell[rho_j].w, net->cell[i].w, n_channels);
//    if (length < min_length) {
//	min_length = length;
//	xi_j = i;
//    }
//  }
//  for (k=0; k<(n_channels); k++)
//    net->cell[rho_N_j].w[k] = (net->cell[rho_j].w[k]+net->cell[xi_j].w[k])/2;
//
//  for (k=0; k<=n_channels; k++) {//kuro
//    net->cell[rho_N_j].am.M[0][k] = net->cell[rho_j].am.M[0][k];//kuro
//    for (i=0; i<=(n_channels); i++){//kuro
//	net->cell[rho_N_j].am.P[k][i] = net->cell[rho_j].am.P[k][i];//kuro
//    }
//  }
//  ///////////////////////
//  min_length = 1e+20;
//  for (i=0; i<n_ivectors; i++) {
//    length = calc_length(net->cell[rho_N_j].w, x[i], n_channels);
//    if (length < min_length) {
//	  min_length = length;
//	  xi_j = i;
//    }
//  }
//  printf("***********xi_j=%d\n",xi_j);
//  for (k=0; k<(n_channels); k++){//reinit to w[rho_N_j] w[rho_j] 
//    //    net->cell[rho_N_j].w[k] = (2.*net->cell[rho_j].w[k]+x[xi_j][k])/3.;
//    net->cell[rho_N_j].w[k] = (net->cell[rho_N_j].w[k]+x[xi_j][k])/2.;
//  }
//  init_AMdata(&(net->cell[rho_N_j].am));
//  //modify_M(net, rho_j,x[xi_j],y[xi_j]);
//  return;
//}
//void reinit_cell_batch1(NET *net, int rho_j, int rho_N_j,FLOAT **x,FLOAT *y,int n_ivectors, FLOAT **train_x,FLOAT *train_y) {
//  FLOAT min_length, length;
//  //  int n_cells = net->n_cells;
//  int n_channels = net->k;
//  int xi_j = 0;
//  int i = 0, k = 0;
//  
//  // for Debug
//  printf(">%3.0f)reinit cell[%3d]w%3.2f,%3.2falpha%3.2e>[%3d]w%3.2f,%3.2falpha%3.2ealpha_%3.2e,[%d]w\n",GlobalTime,
//	   rho_N_j,net->cell[rho_N_j].w[0],net->cell[rho_N_j].w[1],net->cell[rho_N_j].alpha,
//	   rho_j, net->cell[rho_j].w[0],net->cell[rho_j].w[1],net->cell[rho_j].alpha,
//	   net->alpha_bar,xi_j);
//  // Initialize w_(rho_N_j)
//
//  // Search the vector nearest to w_(rho_j)
//  min_length = 1e+20;
//  for (i=0; i<n_ivectors; i++) {
//    length = calc_length(net->cell[rho_j].w, x[i], n_channels);
//    if (length < min_length) {
//	min_length = length;
//	xi_j = i;
//    }
//  }
//  printf("***********xi_j=%d\n",xi_j);
//  for (k=0; k<(n_channels); k++){//reinit to w[rho_N_j] w[rho_j] 
//    //    net->cell[rho_N_j].w[k] = (2.*net->cell[rho_j].w[k]+x[xi_j][k])/3.;
//    net->cell[rho_N_j].w[k] = (net->cell[rho_j].w[k]+x[xi_j][k])/2.;
//  }
//
//  for (k=0; k<=n_channels; k++) {//replace AM[rho_N_j] by AM[rho_j]
//    net->cell[rho_N_j].am.M[0][k] = net->cell[rho_j].am.M[0][k];//kuro
//    for (i=0; i<=(n_channels); i++){//kuro
//	  net->cell[rho_N_j].am.P[k][i] = net->cell[rho_j].am.P[k][i];//kuro
//    }
//  }
//  modify_M(net, rho_j,x[xi_j],y[xi_j]);
//  //  init_AMdata(&(net->cell[rho_N_j].am));modify_M(net, rho_N_j,x[xi_j],y[xi_j]);
//  //  init_AMdata(&(net->cell[rho_j].am));modify_M(net, rho_j,x[xi_j],y[xi_j]);
//  
//  //  for (k=0; k<=n_channels; k++) net->cell[rho_j].am.x[k] = x[xi_j][k];
//  //  net->cell[rho_j].am.y[0] = y[xi_j];
//  //  calc_AM(net->cell[rho_j].am);
//  
//  
//  return;
//}
//
//void reinit_cell_batch0(NET *net, int rho_j, int rho_N_j,FLOAT **x,FLOAT *y,int n_ivectors) {
//  FLOAT min_length, length;
//  int n_cells = net->n_cells;
//  int n_channels = net->k;
//  int xi_j = 0;
//  int i = 0, k = 0;
//  
//  // for Debug
//  printf(">%3.0f)reinit cell[%3d]w%3.2f,%3.2falpha%3.2e>[%3d]w%3.2f,%3.2falpha%3.2ealpha_%3.2e,[%d]w\n",GlobalTime,
//	   rho_N_j,net->cell[rho_N_j].w[0],net->cell[rho_N_j].w[1],net->cell[rho_N_j].alpha,
//	   rho_j, net->cell[rho_j].w[0],net->cell[rho_j].w[1],net->cell[rho_j].alpha,
//	   net->alpha_bar,xi_j);
//  // Initialize w_(rho_N_j)
//  
//  // Search the vector nearest to w_(rho_j)
//  min_length = 1e+20;
//  for (i=0; i<n_cells; i++) {
//    if (i == rho_j) continue;
//    length = calc_length(net->cell[rho_j].w, net->cell[i].w, n_channels);
//    if (length < min_length) {
//	min_length = length;
//	xi_j = i;
//    }
//  }
//  //#define ReinitRandom
//#ifdef ReinitRandom
//  {
//    FLOAT *dw,dwsum=0,aaa;
//    dw=(FLOAT *)malloc(sizeof(FLOAT)*n_channels);
//    for (k=0; k<(n_channels); k++){
//#ifdef MYRANDOM
//	dw[k]=2.*myrandom()-1.;
//#else
//	dw[k]=2.*(rand()/(RAND_MAX+1.))-1.;
//#endif//MYRANDOM
//	dwsum += dw[k]*dw[k];
//    }
//    aaa = 0.5*min_length/sqrt(dwsum);
//    // Initialize w_(rho_N_j)
//    for (k=0; k<n_channels; k++) {
//	net->cell[rho_N_j].w[k] = net->cell[rho_j].w[k]+aaa*dw[k];
//    }
//    free(dw);
//  }
//#else
//  for (k=0; k<(n_channels); k++)
//    net->cell[rho_N_j].w[k] = (net->cell[rho_j].w[k]+net->cell[xi_j].w[k])/2;
//#endif//#ifdef ReinitRandom
////#define ReinitRandom
////#ifdef ReinitRandom
////  {
////    FLOAT aaa=1e-8*(net->xmax-net->xmin)/pow((FLOAT)n_cells,1./n_channels);
////    // Initialize w_(rho_N_j)
////    for (k=0; k<(n_channels); k++) {
////	net->cell[rho_N_j].w[k] = net->cell[rho_j].w[k]
////#ifdef MYRANDOM
////	  +aaa*(2.*myrandom()-1.);
////#else
////	+aaa*(2.*(rand()/(RAND_MAX+1.))-1.);
////#endif//MYRANDOM
////    }
////  }
////#else
////  // Search the vector nearest to w_(rho_j)
////  min_length = 1e+20;
////  for (i=0; i<n_cells; i++) {
////    if (i == rho_j) continue;
////    length = calc_length(net->cell[rho_j].w, net->cell[i].w, n_channels);
////    if (length < min_length) {
////	min_length = length;
////	xi_j = i;
////    }
////  }
////  // for Debug
////  for (k=0; k<(n_channels); k++)
////    net->cell[rho_N_j].w[k] = (net->cell[rho_j].w[k]+net->cell[xi_j].w[k])/2;
////#endif//#ifdef ReinitRandom
//
//  for (k=0; k<=n_channels; k++) {//kuro
//    net->cell[rho_N_j].am.M[0][k] = net->cell[rho_j].am.M[0][k];//kuro
//    for (i=0; i<=(n_channels); i++){//kuro
//#if AM_VER != 3
//	net->cell[rho_N_j].am.P[k][i] = net->cell[rho_j].am.P[k][i];//kuro
//#endif
//    }
//  }
//  /////////////////////////////////////////
//
//  min_length = 1e+20;
//  for (i=0; i<n_ivectors; i++) {
//    length = calc_length(net->cell[rho_N_j].w, x[i], n_channels);
//    if (length < min_length) {
//	min_length = length;
//	xi_j = i;
//    }
//  }
//  printf("***********xi_j=%d\n",xi_j);
//  for (k=0; k<(n_channels); k++){//reinit to w[rho_N_j] w[rho_j] 
//    //    net->cell[rho_N_j].w[k] = (2.*net->cell[rho_j].w[k]+x[xi_j][k])/3.;
//    //    net->cell[rho_N_j].w[k] = (net->cell[rho_j].w[k]+x[xi_j][k])/2.;
//    net->cell[rho_N_j].w[k] = (net->cell[rho_N_j].w[k]+x[xi_j][k])/2.;
//  }
//
//  //  modify_M(net, rho_N_j,x[xi_j],y[xi_j]);
//
//  /////////////////////////////////////////
//
//  // for Debug
//  //  printf("> reinit cell: alpha(rho_N_j:%3.2e[%d], rho_j:%3.2e[%d]w(%3.2f,%3.2f), bar:%3.2e), x_ij:%d\n",
//  //	   net->cell[rho_N_j].alpha, rho_N_j,
//  //	   net->cell[rho_j].alpha, rho_j, net->cell[rho_j].w[1],net->cell[rho_j].w[2],
//  //	   net->alpha_bar, xi_j);
//  //  printf("w[%d]=%3.2f,%3.2f\n",rho_j,net->cell[rho_j].w[1],net->cell[rho_j].w[2]);
//  return;
//}
//
/*====================================================================*
 * store_vector_batch
 *
 *====================================================================*/
//int store_vector_batch0(NET *net, FLOAT ***x, FLOAT **y, int *n_ivectors, int i_times,int *n0_ivectors, FLOAT **train_x, FLOAT *train_y,int n_train) {
//  char name[32] = "store_vector_batch";
//  FLOAT y_m,vmax;
//  int n_cells = net->n_cells;
//  int n_channels = net->k;
//  int i, j, k;
//  //  int t;
//
//  net->c = -1;
//
//  vmax=0;
//  for (i=0; i<n_cells; i++){
//    net->cell[i].v = (FLOAT)n_ivectors[i];
//    if(net->cell[i].v>vmax) vmax=net->cell[i].v;
//  }
//  for (i=0; i<n_cells; i++) net->cell[i].v/=vmax;//kuro
//
//  // STEP 1:  Modify the Associative Matrix
//  //  for (i=0; i<n_cells; i++) modify_M_batch(net, i, x[i], y[i], n_ivectors[i]);
//  for (i=0; i<n_cells; i++) modify_M_batch(net, i, x[i], y[i], n_ivectors[i]);
//  
//  // STEP 2:  Modify the Weight Vector
//  //  for (i=0; i<n_cells; i++) modify_w_batch(net, i, x[i], y[i], n0_ivectors[i]);
//  {
//    FLOAT alpha,dwNorm;
//    dwNorm=0; 
//    for (i=0; i<n_cells; i++){
//	for(k=0;k<n_channels;k++) net->cell[i].dw[k]=0;
//
//	modify_w_batch(net, i, x[i], y[i], n_ivectors[i]);//commentout for test
//
//	for(k=0;k<(n_channels);k++) dwNorm+= square(net->cell[i].dw[k]);
//    }
//    dwNorm=sqrt(dwNorm);
//    // alpha=1./(GlobalTime*100+1);if(alpha<1e-5)alpha=1e-5;//MSE6.596191e-09;1.085559e-06
//    //  alpha=1./(GlobalTime*10);if(alpha<1e-5)alpha=1e-5;//MSE1.334734e-08;5.219017e-06
//    //  alpha=1./(GlobalTime*100);if(alpha<1e-5)alpha=1e-5;//
//    //    alpha=1.e-5+1./(GlobalTime*50.);//f02-> 100 MSE1.354559e-02
//    //    alpha=1.e-5+1./(GlobalTime*5.);//
//    //    alpha=(1.-GlobalTime/GlobalTimeMax);//bad
//    //    alpha=0.10*(1.-GlobalTime/GlobalTimeMax);//not good
//    //    alpha=0.05*(1.-GlobalTime/GlobalTimeMax);//no
//    //    alpha=0.02*(1.-GlobalTime/GlobalTimeMax);//f02MSE1.31e-02;f07MSE7.327e-05
//    alpha=0.05*(1./GlobalTime-1./GlobalTimeMax);//f07MSE6.576e-05,f02MSE1.303592e-02;
//    alpha *= (net->xwidth)/(dwNorm+ZERO)/pow((FLOAT)n_cells,1./n_channels);
//    for (i=0; i<n_cells; i++){
//	for (k=0; k<(n_channels); k++) {
//	  net->cell[i].w[k] = net->cell[i].w[k]+ net->cell[i].dw[k]*alpha;
//	}
//    }
//  }
//  
//  
//  // STEP 3:  Calculate the Total Square Error
//  net->S[i_times] = 0.0;
//
////  for(t=0;t<net->m;t++){//kuro
////    y_m=calc_output(net,train->x[k]);
////    net->cell[i].S += (FLOAT)square(y[i][j]-y_m);
//
//  //#ifdef OLDKURO
//  for (i=0; i<n_cells; i++) {
//    net->cell[i].S = 0.0;
//    for (j=0; j<n_ivectors[i]; j++) {
//	y_m=0;
//	for (k=0; k<=(n_channels); k++) y_m += net->cell[i].am.M[0][k]*x[i][j][k];
//	net->cell[i].S += (FLOAT)square(y[i][j]-y_m);
////	if(i==39){//kuro:debug
////	  printf("y_m=%e,w=%3.2e %3.2e\n",y_m,net->cell[i].w[0],net->cell[i].w[1]);
////	  for (k=0; k<=(n_channels); k++){
////	    printf("(%3.2e)*(%3.2e)",net->cell[i].am.M[0][k],x[i][j][k]);
////	  }
////	  printf("\n");
////	}
//    }
//    //net->cell[i].S /= n_ivectors[i];
//    net->S[i_times] += net->cell[i].S;
//  }
//  //#endif  
//  // STEP 4: Calculate the Asymptotic Optimality. In the Case, Initialize Unit
//  // Calculate to the Value for Every Unit
//  calc_value(net);//calc alpha
//  //calc_value_batch1(net);
//  //calc_value_batch2(net);
//  printf("S[%d]=%e\n",i_times,net->S[i_times]);
//  //  if (net->S[i_times] >= net->S[i_times-1]//
//  //      && GlobalTime>=ReinitTime+5){//
//  //    if (((net->S[i_times] >= net->S[i_times-1])||(GlobalTime<7))//
//  //        && GlobalTime>=ReinitTime+3){//>50 MSE3.610958e-05 r=2
//  //    if (((net->S[i_times] >= net->S[i_times-1])||(GlobalTime<GlobalTimeMax/2.))//
//  //        && GlobalTime>=ReinitTime+4){//>50 MSE7.095745e-05 r=5
//  //    if (((net->S[i_times] >= net->S[i_times-1])||(GlobalTime<GlobalTimeMax/2.))//
//  //        && GlobalTime>=ReinitTime+5){//>50 MSE3.831274e-05 r=5
//  //    if (((net->S[i_times] >= net->S[i_times-1])||(GlobalTime<GlobalTimeMax/4.))
//  //        && GlobalTime>=ReinitTime+GlobalTimeMax/16.){//>50 MSE3.669258e-05
//  if (((net->S[i_times] >= net->S[i_times-1])||(GlobalTime<GlobalTimeMax/5))
//	&& GlobalTime>=ReinitTime+GlobalTimeMax/20){//>  50 MSE3.528878e-05
//    int *rho = NULL;
//    int rho_j, rho_N_j, temp;
//    // ρ_j
//    rho = (int *)malloc(sizeof(int)*n_cells);
//    if (rho == NULL) {
//	printf("error: memory allocation error #?. (%s)\n", name);
//	exit(-1);
//    }
//    for (i=0; i<n_cells; i++) rho[i] = i;
//    
//    // Sort
//    for (j=0; j<n_cells-1; j++) {
//	for (k=j+1; k<n_cells; k++) {
//	  if (net->cell[(int)rho[j]].alpha < net->cell[(int)rho[k]].alpha) {//大きい順 kuro '<' not '>'
//	    temp = rho[k];
//	    rho[k] = rho[j];
//	    rho[j] = temp;
//	  }
//	}
//    }
//    
//    /*
//	// for Debug
//	for (j=0; j<n_cells; j++) {
//	printf("j:%d, alpha[%d]:%3.2e\n", j, rho[j], net->cell[rho[j]].alpha);
//	}
//    */
//    
//    // Condition of Reinit Cell
//    for (j=0; j<n_cells; j++) {
//	rho_j = rho[j]; rho_N_j = rho[n_cells-(j+1)];
//	if((net->alpha_NG==0)
//	   //	 && rho_j!= rho_N_j
//	   && (net->cell[rho_j].alpha > net->v_ratio*net->alpha_bar)
//	   //	 && (net->cell[rho_N_j].v < net->v_thresh/2.)){
//	   && (net->cell[rho_N_j].v < net->v_thresh)){
//	   //	 && (net->cell[rho_N_j].v < net->v_thresh))
//	  reinit_cell_batch(net, rho_j, rho_N_j,x[rho_j],y[rho_j],n_ivectors[rho_j],train_x,train_y,n_train);
//	  //	break;//kuro comment out kuro
//	  ReinitTime=GlobalTime;
//	}
//	//	else {
//	//	  printf("NG=%d!=0 or %3.2e<%3.2e*%3.2e or %3.2e>%3.2e\n",net->alpha_NG, 
//	//		 net->cell[rho_j].alpha,net->v_ratio,net->alpha_bar,
//	//		 net->cell[rho_N_j].v,net->v_thresh);
//	//	}
//    }
//    free(rho);
//  }
//    
//  return(0);
//}
//-->init_batch_ivectors in my_function.c

void init_Voronoi(VORONOI *V, int n_cells, int n_channels, int n_train)
{
  int i;
  V->i2v=(int *)malloc(sizeof(int)*n_cells);//
  V->i2v0=(int *)malloc(sizeof(int)*n_cells);//
  V->t2i=(int *)malloc(sizeof(int)*n_train);//
  V->ij2t =(int **)malloc(sizeof(int*)*(n_cells));
  //  V->it2d2=(FLOAT **)malloc(sizeof(FLOAT *)*(n_cells));
  for (i=0; i<n_cells; i++) {
    V->ij2t[i]  = (int *)malloc(sizeof(int)*(n_train));
  }
//  for (i=0; i<n_train; i++) {
//    V->it2d2[i] = (FLOAT *)malloc(sizeof(FLOAT)*(n_train));
//  }
}
//void free_Voronoi(VORONOI *V,int n_cells, int n_channels, int n_train)
//{
//  int i;
//  for (i=0; i<n_cells; i++) {
//    free(V->ij2t[i]);
//  }
////  for (i=0; i<n_train; i++) {
////    free(V->it2d2[i]);
////  }
//  free(V->i2v);
//  free(V->i2v0);
//  free(V->t2i);
//  free(V->ij2t);
//  //  free(V->it2d2);
//}
//void calc_Voronoi(VORONOI *V, FLOAT **x, FLOAT *y, FLOAT **w, 
//		      int n_cells, int n_channels, 
//		      int n_train, NET *net) {
void calc_Voronoi(NET *net, FLOAT **x, FLOAT *y, int n_train) {
  char name[32] = "init_batch_ivector";
  FLOAT length, min_length;
  int i_min_length;
  int i, j, t;
  //  int *m;
  FLOAT **w=net->w;
  int n_cells = net->n_cells;
  int n_channels = net->k;
  VORONOI *V=net->V;
  int ii;//040109 for remove Unit
  
  // bidata->m，bidata->sの初期化
  for (i=0; i<n_cells; i++){
    if(V->i2v[i]<0) continue; //removed unit
    V->i2v[i] = 0;//Viに属すxの個数
  }
  
//  for (t=0; t<n_train; t++) fprintf(stdout,"x[%d]=%g %g %g\n",t,x[t][0],x[t][1],x[t][2]); //20190107
//  for (i=0; i<n_cells; i++) fprintf(stdout,"w[%d]=%g %g\n",i,w[i][0],w[i][1]); //20190107

  for (t=0; t<n_train; t++) {
    min_length = 1e30;//
    i_min_length = 0;
    
    for (i=0; i<n_cells; i++) {
      //      V->it2d2[i][t]=length = distance2(x[t], w[i], n_channels);
      //      if(w[i][n_channels]<0) continue;//040110<<damedayo
      if(V->i2v[i]<0) continue;//040115removed unit
      length = distance2(x[t], w[i], n_channels);
      if (length < min_length) {
	min_length = length;
	i_min_length = i;
      }
    }
    V->t2i[t] = i_min_length;
    V->ij2t[i_min_length][V->i2v[i_min_length]]=t;
    ++(V->i2v[i_min_length]);
  }

  //040109 for remove Unit
  //vmin2=0 >14ssp1.3229e-06(1.3387e-04)N500k2w0.8vm3
  for (ii=0; ii<n_cells; ii++){
    net->n_cells2=n_cells;
    if(V->i2v[ii]<=net->vmin2){//vmin2 is for neglecting the cell
      for (j=0; j<V->i2v[ii]; j++){
	t=V->ij2t[ii][j];
	min_length = 1e30;//
	i_min_length = 0;
	for (i=0; i<n_cells; i++){
	  if(ii==i) continue;
	  length = distance2(x[t], w[i], n_channels);
	  if (length < min_length) {
	    min_length = length;
	    i_min_length = i;
	  }
	}
	V->t2i[t] = i_min_length;
	V->ij2t[i_min_length][V->i2v[i_min_length]]=t;
	++(V->i2v[i_min_length]);
      }
      //      fprintf(stderr,"(%d)%d<vmin2",ii,V->i2v[ii]); 
#ifdef DEBUGvm
      fprintf(stderr,"(%d:vm2)",ii); 
#endif
      V->i2v[ii]=-1e3;//040212 remove?
      //V->i2v[ii]=-1.;//060609
      net->n_cells2--;//040520
      //      V->i2v[ii]=0;//040110
      for(j=0;j<n_channels;j++) w[ii][j]=0;//040110
      //      w[ii][n_channels]=-1;//040110ダメ
      //      for(j=0;j<n_channels;j++) w[ii][j]=1e10;//040110
      //      w[ii][n_channels]=-1e10;//040110
      net->n_cells2--;//040524
    }
  }
  
  // Complete the Input Vectors
  // 荷重ベクトルが参照する入力ベクトルの個数は，
  // 最低でも入力次数以上ないと，連想行列Mの値が保証されない
  {
    int n1, temp;
    int n2,k,tt,n0,n10;
    FLOAT *d2 = (FLOAT *)malloc(sizeof(FLOAT)*n_train), ltemp;
    int *s = (int *)malloc(sizeof(int)*n_train);
    
    if (d2 == NULL || s == NULL) {
      printf("error: memory allocation error #2. (%s)", name);exit(-1);
    }
    V->vmax=0;
#define NEWN1
#undef NEWN1
#ifdef NEWN1
    for (i=0; i<n_cells; i++) {
      V->i2v0[i]=V->i2v[i];
      n1=n_channels+1;
      if(V->i2v[i] >= n1) {//入力次数以上ないとき
	net->cell[i].v0=1;
      }
      else if(net->cell[i].v0==0){
#ifdef emacseditor
      }}
#endif
#else
    //      n1=n_channels+1;
    //      if(SuccessiveLearn) n1=V->i2v[i];
    //      n1*=2;//???
    //      n1=V->i2v[i]+(n_channels+1.-V->i2v[i])*alphat(GlobalTime);
    //n1=(V->i2v[i]+2)*alphat(GlobalTime);
    //      n1=4;//good for mackey-glass
    n1=net->vmin;//good for mackey-glass
    //    n1=net->vmin/(1.+GlobalTime/40.)+0.5;printf("!!!!!!!!!n1=%d\n",n1);//040107
    //      n1=3;//
    for (i=0; i<n_cells; i++) {
      V->i2v0[i]=V->i2v[i];
      //      if(n1>n_channels+1) n1=n_channels+1;
      //            if(V->i2v[i]<n1 && GlobalTime<=1){//good}
      if(V->i2v[i]<=net->vmin2){//040109 for remove Unit
	;//	printf("{%d}",i); 
      }
      else if(V->i2v[i]<n1){//better
	//n1=(n_channels+1.)*alphat(GlobalTime);
	//	n1=V->i2v[i]+(n_channels+1.-V->i2v[i])*alphat(GlobalTime);
	//	n1=(3)*alphat(GlobalTime);
	//	n1=V->i2v[i];//??
	//	n1=(V->i2v[i]+3)*alphat(GlobalTime/2.);
	//	n1=(V->i2v[i]+3)*alphat(GlobalTime);
	//	printf("<%d,%d>",V->i2v[i],n1);
	//	fprintf(stderr,"[%d]%d<vmin",i,V->i2v[i]);
#ifdef DEBUGvm
	fprintf(stderr,"[%d:%d<vm]",i,V->i2v[i]);
#endif
	net->n_cells2--;//040524
	//	n1=(10)*alphat(GlobalTime);
#endif
	net->cell[i].v0=1;
	n0=V->i2v[i];
	n2=n1-n0;
	//	if(n2>=n_cells) n2=n_cells;//check 060110
	//n2=n1*3-n0;
	//n2=n1*2-n0;
	n10=n1-n0;
	tt=0;
	//入力ベクトルと荷重ベクトルの距離による入力ベクトルのソート
	//	for (t=0; t<n_train; t++) s[t]=t;//check 060110
	for (t=0; t<n_train; t++){
	  //if(V->t2i[t]!=i) d2[t] = distance2(w[i],x[t],n_channels);//距離
	  //if(V->t2i[t]==i) continue;//080612 comment out
	  d2[t] = distance2(w[i],x[t],n_channels);//距離
	  s[tt++] = t;//ソート準備
	}
	// ソート
	for (t=0; t<n2; t++) {
	  for (j=t+1; j<tt; j++) {
	    if (d2[(int)s[t]] > d2[(int)s[j]]) {
	      temp = s[j];
	      s[j] = s[t];
	      s[t] = temp;
	    }
	  }
	}
	//kuro
	if(GlobalTime>1){//kuro
	  //近似誤差による入力ベクトルのソート
	  for(j=0;j<n2;j++){
	    ltemp =y[s[j]];
	    for(k=0; k<=(n_channels); k++) 
	      ltemp -= net->cell[i].am.M[0][k]*x[s[j]][k];
	    d2[s[j]]=fabs(ltemp); //誤差
	  }
	}
	for(j=0;j<n10;j++){//Sort
	  for(k=j+1;k<n2;k++){
	    if(d2[s[j]]>d2[s[k]]){
	      temp=s[j];
	      s[j]=s[k];
	      s[k]=temp;
	    }
	  }
	}
	for (j=0; j<n10; j++) V->ij2t[i][j+n0] = s[j];
	V->i2v[i] = n1;
      }
      if(V->vmax<V->i2v[i]) V->vmax=V->i2v[i];
    }
    //    printf("V->vmax=%d\n",V->vmax);
    free(s);
    free(d2);
  }
  for (i=0; i<n_cells; i++) net->cell[i].v=(FLOAT)V->i2v[i]/V->vmax;
  //  if(net->v_thresh<(net->k+1.)/V->vmax) net->v_thresh=(net->k+1)/V->vmax;
  //  printf("v_thresh=%e\n",net->v_thresh);
  net->vmax=V->vmax;
  return;
}
/*====================================================================*
 * create_batch_wvector
 * 荷重ベクトルの初期値データを作成
 *====================================================================*/
FLOAT** create_batch_wvector(int n_cells, int n_channels) {
  char name[32] = "create_batch_wvector";
  FLOAT **w = NULL;
  int i;
  
  w = (FLOAT **)malloc(sizeof(FLOAT *)*n_cells);
  if (w == NULL) {
    printf("error: memory allocation error #2. (%s)\n", name);
    exit(-1);
  }
  for (i=0; i<n_cells; i++) {
    w[i] = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    if (w[i] == NULL) {
      printf("error: memory allocation error #3-%d. (%s)\n", i, name);
      exit(-1);
    }
  }

  return(w);
}
/*====================================================================*
 * init_batch_wvector
 * 荷重ベクトルの初期値データを初期化
 *====================================================================*/
 void init_batch_wvector(NET *net, FLOAT **x, int n_train) {
  //  FLOAT **w = NULL;
  int *t = NULL;
  int time, dupli;
  int i, j, k;
  int n_cells = net->n_cells;
  int n_channels = net->k;
  FLOAT **w=net->w;

  if(n_cells>n_train*0.632){
    printf("*** n_cells(%d) may be too big w.r.t n_train(%3d) !!?\n",n_cells,n_train);
    //    if(n_cells>n_train) {fprintf(stderr,"*** n_cells(%d) bigger than n_train(%3d) !!?\n",n_cells,n_train);exit(-1);}
    //    if(n_cells>n_train*0.80) {fprintf(stderr,"*** n_cells(%d) bigger than n_train(%3d) !!?\n",n_cells,n_train);}
  }
  //  w = create_batch_wvector(n_cells, n_channels);
  t = (int *)malloc(sizeof(int)*n_cells);

  i = 0;
  int smallntrain=(n_cells<n_train*0.632);
  while (i < n_cells) {
    dupli = 0;
    // t は 0〜(n_train-1) の間からの乱数
    int RND=RAND;
#if RAND == RANDOM
    //    double RNDM=random(); unsigned long RM=RAND_MAX; int RND=RAND;
    time = (int)(rand()/(RAND_MAX+1.0)*(n_train));//
#elif RAND == DRAND48
    time = (int) (drand48()*(n_train));
#elif RAND == MYRAND
    //    double RNDM=myrandom(); unsigned long RM=RAND_MAX; 
    time = (int)((double)myrandom()*(n_train));
#elif RAND == ZMTRAND
    time = (int)(NextUnif()*(n_train));
#endif

//#ifdef MYRANDOM
//    time = (int)(myrandom()*(n_train));
//#else
//    //    time = (int)(rand()%(n_train));//bad! see man rand.
//    time = (int)(rand()/(RAND_MAX+1.0)*(n_train));//
//#endif //MYRANDOM

    // Avoid duplication
    for (j=0; j<i; j++) {
      if (time == t[j]) {
	dupli = 1;
	break;
      }
    }
    if (dupli == 1 && smallntrain) continue;
    //    if (dupli == 1) continue;
      
    // Set Data
    //    printf("%d)\n",i);
    //    printf("net['w'][%d,:]=x[%d][0:n_channels]\n",i,time); //check c2py20190115    
    for (k=0; k<(n_channels); k++){
      //      x[time][k]=0;
      w[i][k] = x[time][k];
    }
    t[i] = time;
    ++i;
  }
  free(t);
    //{fprintf(stderr,"1GTime=%d\n",GlobalTime);int *kkkk=(int*)malloc(sizeof(int)*1000000); fprintf(stderr,"&k=%p,i=%d\n",kkkk,i);free(kkkk);}//test fedora core error
  return;
}

/*====================================================================*
 * remove_batch_wvector
 * 荷重ベクトルの初期値データを破棄
 *====================================================================*/
void remove_batch_wvector(FLOAT **w, int n_cells) {
  int i;

  for (i=0; i<n_cells; i++) free(w[i]);
  free(w);
  return;
}

void init_net_batch(NET *net, FLOAT **x, int n_train){
  int i;
  if(net->winit==0){
    init_batch_wvector(net, x, n_train);//
    net->winit=1;
  }
  if(net->Vinit==0){
    net->V=(VORONOI *)malloc(sizeof(VORONOI)*1);
    init_Voronoi(net->V,net->n_cells,net->k,n_train);
    for(i=0;i<net->n_cells;i++) net->cell[i].v0=0;
    net->Vinit=1;
  }
  //  net->nentropy_thresh=0.95-(0.95-0.85)*(500-net->n_cells)/(500-100);//adhoc??
  //  net->nentropy_thresh=0.95-(0.95-0.8)*(500-net->n_cells)/(500-100);//adhoc??
  //  net->nentropy_thresh=LA(net->n_cells,0.80,100,0.95,500);
  //  net->nentropy_thresh=LA(net->n_cells,0.80,200,0.90,500);
  if(net->nentropy_thresh<0)  net->nentropy_thresh=LA(net->n_cells,0.75,100,0.90,500);
  //  net->N_reinit_max   =LA(net->n_cells,0.10,100,0.30,500)*net->n_cells;
  net->N_reinit_max   =LA(net->n_cells,0.10,100,0.30,500)*net->n_cells;
  //  net->gamma0         =0.01; //good
  //  net->gamma0         =0.001; //good
  //  net->gamma0         =0.05; //good
  //  net->gamma0         =LA(net->n_cells,0.24,100,0.11,500); //good
  //    net->gamma0         =2.7*pow(net->n_cells,-1./net->k); //good50MSE1.559e-06N50050MSE2.174e-02N100
  //net->gamma0         =2.6*pow(net->n_cells,-1./net->k); //good50MSE1.591e-06N50050MSE2.261e-02N100
  //net->gamma0         =2.8*pow(net->n_cells,-1./net->k);     //good50MSE1.672e-06N500
  //ボロノイ領域の半径はpow(net->n_cells,-1./net->k)に比例。学習率もそれに比例させる?
  net->Tgamma         =5;
  //  net->Tgamma         =LA(net->n_cells,3.0,100,5.0,500); 
  net->printf(strx4("entropy_thresh=%f,N_reinit_max=%d,gamma0=%f,Tgamma=%f\n",
                   net->nentropy_thresh, net->N_reinit_max,net->gamma0,net->Tgamma));
}
void store_vector_batch(NET *net, FLOAT **x_train, FLOAT *y_train,int n_train, int phase){

  int ret2=0;

  if(phase==0){
    //STEP 0:
    //    calc_Voronoi(V,train->x,train->y,net->w,n_cells, n_channels, n_train,net);
    //    printf("laptime=%5.3f\n",mytimer_lap());
//    {
//	    int i;
//	    for(i=0;i<n_train;i++) printf("y[%d]=%e\n",i,y_train[i]);
//    }
    calc_Voronoi(net,x_train,y_train, n_train);
#ifdef LAPTIME
    printf("laptime:calc_Voronoi=%5.3f\n",mytimer_lap());
#endif
    // STEP 1:  Modify the Associative Matrix
#if GSL ==0
    modify_M_batch_RLS(net, x_train, y_train);
#else
    if(net->Tpinv>0 && GlobalTime>=net->Tpinv){
      modify_M_batch_pinv(net, x_train, y_train);
      net->printf(strx1("%d #modify_M_batch_pinv\n",net->GlobalTime));
    }
    else {
      modify_M_batch_RLS(net, x_train, y_train);
      net->printf(strx1("%d #modify_M_batch_RLS\n",net->GlobalTime));
    } 
#endif

#ifdef LAPTIME
    printf("laptime:modify_M_batch=%5.3f\n",mytimer_lap());
#endif
    // STEP 2:  Modify the Weight Vector
    // STEP 3:  Calculate the Total Square Error
    // STEP 4: Calculate the Asymptotic Optimality. In the Case, Initialize Unit
    // Calculate to the Value for Every Unit
    net->ret1=calc_alpha(net, x_train, y_train,n_train);
#ifdef LAPTIME
    printf("laptime:calc_alpha=%5.3f\n",mytimer_lap());
#endif
  }
  if(phase==1){
    if(net->ret1){
      ret2=reinit_cell_batch(net, x_train, y_train,n_train) ;
#ifdef LAPTIME
      printf("laptime:reinit_cell =%5.3f\n",mytimer_lap());
#endif
    }
    //    printf("check ret=%d,%d\n",net->ret1,ret2);
    if(net->ret1==0 || ret2==0){
      modify_w_batch(net, x_train, y_train, n_train);
//      int i,j;
//      for(i=0;i<net->n_cells;i++){
//        printf("check w[%d]=",i);
//        for(j=0;j<net->k;j++){
//          printf("%g ",net->cell[i].w[j]);
//        }
//        printf("\n");
//      }
#ifdef LAPTIME
      printf("laptime:modify_w_batch =%5.3f\n",mytimer_lap());
#endif
    }
  }
  
  //    if(calc_alpha(net, x_train, y_train,n_train) && 
  //	 reinit_cell_batch(net, x_train, y_train,n_train)) ;
  //    else modify_w_batch(net, x_train, y_train, n_train);
  //    printf("reinit or modify_w_batch laptime=%5.3f\n",mytimer_lap());
}
void clear_net(NET *net)
{
  int i;
  if(net->Vinit){
    for (i=0; i<net->n_cells; i++) free((net->V)->ij2t[i]);
    free((net->V)->ij2t);
    free((net->V)->t2i);
    free((net->V)->i2v0);
    free((net->V)->i2v);
    free(net->V);
    net->Vinit=0;
  }
  if(net->init){
    for (i=0; i<net->n_cells; i++) {
      free(net->w[i]);
      free(net->dw[i]);
      free_AM(&(net->cell[i].am));
    }
    free(net->w);
    free(net->dw);
    free(net->cell);
    net->init=0;
  }
}
//#define HEN
//#ifdef HEN
//int net_save(NET *net,char *fname){
//  int i,k,l;
//  char buff[256];
//  FILE *fp;
//  if(fname==NULL){
//    printf("Name of the Net File:");
//    buff[0]=0;scanf1("%s",buff);
//    if(buff[0]==0) return(1);
//    if((fp = open_file(buff, "wt"))==NULL) return(1);
//  }
//  else{
//    if((fp = open_file(fname, "wt"))==NULL) return(1);
//  }
//  fprintf(fp,"%d %d %d %d\n",net->n_cells,net->winit,GlobalTime,ReinitTime);
//  //
//  for(i=0;i<net->n_cells;i++){
//    for(k=0; k<net->k; k++) fprintf(fp,"%e ",net->cell[i].w[k]);
//    fprintf(fp,"\n");
//    for(k=0; k<=net->k; k++) fprintf(fp,"%e ",net->cell[i].am.M[0][k]);
//    fprintf(fp,"\n");
//    for(k=0; k<=net->k; k++){
//	for(l=0; l<=net->k; l++) fprintf(fp,"%e ",net->cell[i].am.P[k][l]);
//	fprintf(fp,"\n");
//    }
//  }
//  fclose(fp);
//  return(1);
//}
//int net_load(NET *net,char *fname){
//  //  FLOAT i_times, d_times, n_cells, n_compare, v_thresh, v_ratio, width;
//  //  int n_channels;
//  int i, k, l;
//  FILE *fp;
//  char buff[1024];
//  FLOAT dummy;
//
//  if(fname==NULL){
//    printf("Name of the Net File:");
//    buff[0]=0;scanf1("%s",buff);
//    if(buff[0]==0) return(0);
//    if((fp = open_file(buff, "r"))==NULL) return(1);
//  }
//  else{
//    if((fp = open_file(fname, "r"))==NULL) return(1);
//  }
//  clear_net(net);
////  fread(net,sizeof(NET),1,fp);
//  fscanf(fp,"%d %d %d %d\n",&net->n_cells,&net->winit,&GlobalTime,&ReinitTime);
//  GlobalTimeMax=GlobalTime;
//  //  ReinitTime=GlobalTime;
//  net->init=1;
//  net->Vinit=0;
//
//  //  net->k1 = (int)n_channels+1;
//  
//  net->cell = (CELL *)malloc(sizeof(CELL)*net->n_cells);
//  net->w  = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);
//  net->dw = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);
//
//  for (i=0; i<net->n_cells; i++) {
//    net->cell[i].w  = net->w[i]  = (FLOAT *)malloc(sizeof(FLOAT)*net->k);
//    net->cell[i].dw = net->dw[i] = (FLOAT *)malloc(sizeof(FLOAT)*net->k);
//    //    for (k=0; k<net->k; k++) net->cell[i].w[k] = 0.0;
//    net->cell[i].S = 0.0;
//    net->cell[i].alpha = 0.0;
//    net->cell[i].enough = 1;
//    init_AM(&(net->cell[i].am), net->k+1, 1);
//  }
//  //
//  //    for (k=0; k<net->k; k++) fprintf(fp,"%e",net->cell[i].w[k]);
//  for(i=0;i<net->n_cells;i++){
//    for(k=0; k<net->k; k++) fscanf(fp,"%lf",&(net->cell[i].w[k]));
//    for(k=0; k<=net->k; k++) fscanf(fp,"%lf",&(net->cell[i].am.M[0][k]));
//    for(k=0; k<=net->k; k++){
//	for(l=0; l<=net->k; l++){
//	  fscanf(fp,"%lf",&(net->cell[i].am.P[k][l]));
//	  //	fscanf("%lf",&(dummy));
//	  //	  if(k==l) net->cell[i].am.P[k][l]=1e4;//??040116
//	  //	  if(k==l) net->cell[i].am.P[k][l]=1e3;//??040116
//	  //	  else net->cell[i].am.P[k][l]=0;//??040116
//	}
//    }
//  }
//  fclose(fp);
//  return(1);
//}
//#endif//#ifdef HEN
//#ifdef FURUI
int net_save(NET *net,char *fname){
  int i,k,l;
  char buff[256];
  FILE *fp;
  if(fname==NULL){
    net->printf("Name of the Net File:");
    buff[0]=0;scanf1("%s",buff);
    if(buff[0]==0) return(1);
    if((fp = open_file(buff, "wt"))==NULL) return(1);
  }
  else{
    if((fp = open_file(fname, "wt"))==NULL) return(1);
  }
  //  fwrite(net,sizeof(NET),1,fp);// net->winit=1;

  //(gdb print net[0]
  //$2 = {cell = 0x85066e8, S = 0x85306e0, l_mode = 1, i_times = 10, d_times = 10, k = 2, k1 = 0, c = 0, l = 0, 
  //  n_cells = 500, n_compare = 6, n_fires = 500, t_ri = 0, alpha_min = 0, alpha_bar = 0, alpha_hat = 0, sigma2_hat = 0, 
  //  alpha_NG = 0, v_thresh = 0.5, v_ratio = 5, tau_E = 400000, r_E = 0.99999750000312504, width = 0.80000000000000004, 
  //  xmin = 5.5506830000000002e-06, xmax = 0.99999570000000004, xwidth = 0.99999014931700003, 
  //  ymin = -0.0099984819999999995, ymax = 0.23221910000000001, ywidth = 0.24221758200000001, 
  //  xmin0 = 5.5506830000000002e-06, xmax0 = 0.99999570000000004, xwidth0 = 0, ymin0 = -0.0099984819999999995, 
  //  ymax0 = 0.23221910000000001, ywidth0 = 0, w = 0x850e3f0, dw = 0x850ebc8, V = 0x8575608, r1 = 0, r2 = 0, nr = 0, 
  //  r = 0x0, r12 = 0, R = 0x0, R12 = 0, vmin = 3, vmin2 = -1, init = 1, Vinit = 1, winit = 1}
  //  //  
  //typedef struct {
  //  CELL *cell;		// セル（ユニット）構造体
  //  FLOAT *S;		// セルの平均２乗誤差和
  //  int l_mode;// 学習モード（ONLINE or BATCH MODE）
  //  int i_times;		// 学習（繰り返し）回数
  //  int d_times;		// 表示回数（学習してる時に予測結果を表示）
  //  int k;		// 入力次数
  //  int k1;		// 入力次数+1
  //  int c;		// 競合に勝ったセルのインデックス
  //  int l;		// 価値alpha_iの最小となるインデックス
  //  int n_cells;		// ネット全体のセル数
  //  int n_compare;	// 入力が来たとき比較する周囲のセル数
  //  int n_fires;		// 発火したセル数
  //  int t_ri;		// 再初期化のための学習回数
  //  FLOAT alpha_min;	// 価値alpha_iの最小値
  //  FLOAT alpha_bar;	// 価値alpha_iの平均値
  //  FLOAT alpha_hat;
  //  FLOAT sigma2_hat;
  //  int alpha_NG;		// 価値alpha_iのいずれか１つでも負数になったとき真
  //  FLOAT v_thresh;	// 連想記憶が最小限記憶すべきベクトルの数
  //  FLOAT v_ratio;	// θ_alpha（再初期化条件のしきい値）
  //  FLOAT tau_E;		// ？（再初期化で用いる定数）
  //  FLOAT r_E;		// 忘却係数
  //  FLOAT width;		// Voronoi領域の境界の幅
  //  FLOAT xmin,xmax,xwidth;
  //  FLOAT ymin,ymax,ywidth;
  //  FLOAT xmin0,xmax0,xwidth0;
  //  FLOAT ymin0,ymax0,ywidth0;
  //  FLOAT xmin1,xmax1,xwidth1;
  //  FLOAT ymin1,ymax1,ywidth1;
  //  FLOAT **w,**dw;
  //  VORONOI *V;
  //  int r1,r2,nr; //p1/p2がデータの分解能(p1=0→実数, p1=p2=1→整数型)
  //  FLOAT *r,r12;
  //  FLOAT *R,R12;
  //  int vmin;
  //  int vmin2;
  //  int init;
  //  int Vinit;
  //  int winit;
  //}NET;
  fwrite(&(net->n_cells),sizeof(int),1,fp);
  fwrite(&(net->winit),sizeof(int),1,fp);
  fwrite(&(GlobalTime),sizeof(int),1,fp);
  fwrite(&(ReinitTime),sizeof(int),1,fp);
  //from here, new data 051208
  fwrite(&(net->k),sizeof(int),1,fp);
  fwrite(&(net->n_compare),sizeof(int),1,fp);
  fwrite(&(net->v_thresh),sizeof(FLOAT),1,fp);
  fwrite(&(net->vmin),sizeof(int),1,fp);
  fwrite(&(net->vmin2),sizeof(int),1,fp);
  fwrite(&(net->vmax),sizeof(int),1,fp);
  fwrite(&(net->v_ratio),sizeof(FLOAT),1,fp);
  fwrite(&(net->width),sizeof(FLOAT),1,fp);
  fwrite(&(net->ymin),sizeof(FLOAT),1,fp);
  fwrite(&(net->ymax),sizeof(FLOAT),1,fp);
  fwrite(&(net->ywidth),sizeof(FLOAT),1,fp);
  fwrite(&(net->ymin0),sizeof(FLOAT),1,fp);
  fwrite(&(net->ymax0),sizeof(FLOAT),1,fp);
  fwrite(&(net->ywidth),sizeof(FLOAT),1,fp);
  fwrite(&(net->xmin),sizeof(FLOAT),1,fp);
  fwrite(&(net->xmax),sizeof(FLOAT),1,fp);
  fwrite(&(net->xwidth),sizeof(FLOAT),1,fp);
  fwrite(&(net->r1),sizeof(int),1,fp);
  if(net->r1>0){
    fwrite(&(net->r2),sizeof(int),1,fp);
    fwrite(&(net->nr),sizeof(int),1,fp);
    fwrite(&(net->r12),sizeof(FLOAT),1,fp);
    for(i=0;i<=net->nr;i++){
      fwrite(&(net->R[i]),sizeof(FLOAT),1,fp);
      fwrite(&(net->r[i]),sizeof(FLOAT),1,fp);
    }
  }
  //to here, new data 051208
  //
  for(i=0;i<net->n_cells;i++){
    for(k=0; k<net->k; k++) fwrite(&(net->cell[i].w[k]),sizeof(FLOAT),1,fp);
    for(k=0; k<=net->k; k++) fwrite(&(net->cell[i].am.M[0][k]),sizeof(FLOAT),1,fp);
    for(k=0; k<=net->k; k++){
      for(l=0; l<=net->k; l++){
	fwrite(&(net->cell[i].am.P[k][l]),sizeof(FLOAT),1,fp);
      }
    }
    fwrite(&(net->cell[i].v),sizeof(FLOAT),1,fp);
    fwrite(&(net->cell[i].E),sizeof(FLOAT),1,fp);
    fwrite(&(net->cell[i].S),sizeof(FLOAT),1,fp);
  }

  fclose(fp);
  return(1);
}
int net_load(NET *net,char *fname){
  //  FLOAT i_times, d_times, n_cells, n_compare, v_thresh, v_ratio, width;
  //  int n_channels;
  int i, k, l;
  FILE *fp;
  char buff[1024];
  //  FLOAT dummy;

  if(fname==NULL){
    net->printf("Name of the Net File:");
    buff[0]=0;scanf1("%s",buff);
    if(buff[0]==0) return(0);
    if((fp = open_file(buff, "r"))==NULL) return(1);
  }
  else{
    if((fp = open_file(fname, "r"))==NULL) return(1);
  }
  
  clear_net(net);
  fread(&(net->n_cells),sizeof(int),1,fp);
  fread(&(net->winit),sizeof(int),1,fp);  
  fread(&(GlobalTime),sizeof(int),1,fp);
  fread(&(ReinitTime),sizeof(int),1,fp);
  //from here, new data 051208
  fread(&(net->k),sizeof(int),1,fp);
  fread(&(net->n_compare),sizeof(int),1,fp);
  fread(&(net->v_thresh),sizeof(double),1,fp);
  fread(&(net->vmin),sizeof(int),1,fp);
  fread(&(net->vmin2),sizeof(int),1,fp);
  fread(&(net->vmax),sizeof(int),1,fp);
  fread(&(net->v_ratio),sizeof(double),1,fp);
  fread(&(net->width),sizeof(double),1,fp);
  fread(&(net->ymin),sizeof(FLOAT),1,fp);
  fread(&(net->ymax),sizeof(FLOAT),1,fp);
  fread(&(net->ywidth),sizeof(FLOAT),1,fp);
  fread(&(net->ymin0),sizeof(FLOAT),1,fp);
  fread(&(net->ymax0),sizeof(FLOAT),1,fp);
  fread(&(net->ywidth),sizeof(FLOAT),1,fp);
  fread(&(net->xmin),sizeof(FLOAT),1,fp);
  fread(&(net->xmax),sizeof(FLOAT),1,fp);
  fread(&(net->xwidth),sizeof(FLOAT),1,fp);
  fread(&(net->r1),sizeof(int),1,fp);
  if(net->r1>0){
    fread(&(net->r2),sizeof(int),1,fp);
    fread(&(net->nr),sizeof(int),1,fp);
    fread(&(net->r12),sizeof(FLOAT),1,fp);
    for(i=0;i<=net->nr;i++){
      fread(&(net->R[i]),sizeof(FLOAT),1,fp);
      fread(&(net->r[i]),sizeof(FLOAT),1,fp);
    }
  }
  //to here, new data 051208

  GlobalTimeMax=GlobalTime;
  net->init=1;
  net->Vinit=0;

  net->cell = (CELL *)malloc(sizeof(CELL)*net->n_cells);
  net->w  = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);
  net->dw = (FLOAT **)malloc(sizeof(FLOAT *)*net->n_cells);

  for (i=0; i<net->n_cells; i++) {
    net->cell[i].w  = net->w[i]  = (FLOAT *)malloc(sizeof(FLOAT)*net->k);
    net->cell[i].dw = net->dw[i] = (FLOAT *)malloc(sizeof(FLOAT)*net->k);
    net->cell[i].S = 0.0;
    net->cell[i].alpha = 0.0;
    net->cell[i].enough = 1;
    init_AM(&(net->cell[i].am), net->k+1, 1);
  }
  for(i=0;i<net->n_cells;i++){
    for(k=0; k<net->k; k++) fread(&(net->cell[i].w[k]),sizeof(FLOAT),1,fp);
    for(k=0; k<=net->k; k++) fread(&(net->cell[i].am.M[0][k]),sizeof(FLOAT),1,fp);
    for(k=0; k<=net->k; k++){
      for(l=0; l<=net->k; l++){
	fread(&(net->cell[i].am.P[k][l]),sizeof(FLOAT),1,fp);
      }
    }
    fread(&(net->cell[i].v),sizeof(FLOAT),1,fp);
    fread(&(net->cell[i].S),sizeof(FLOAT),1,fp);
    fread(&(net->cell[i].E),sizeof(FLOAT),1,fp);
  }
  fclose(fp);
  net->nEns=1;

  return(1);
}

NET *net_loads(NET *net){
  char buff[5000];
  char *p,*fname[1024];
  int j,nEns=0;

  net->printf("Name(s) of Net File(s):");
  fgets(buff,5000,stdin);
  {
    p=buff;
    for(;;){
      if(*p==0) break;
      for(;;p++) if(isgraph(*p) || *p==0) break;//search next no-space
      if(*p==0) break;
      fname[nEns++]=p;
      p++;
      for(;;p++) if(!isgraph(*p) || *p==0) break;//search next space
      *p++=0;
    }
  }
  //  fprintf(stdout,"Number net files = %d.\n",nEns);
  if(nEns<0) return(NULL);

  net=(NET *)malloc(sizeof(NET)*nEns);
  for(j=0;j<nEns;j++){
    net_load(&net[j],fname[j]);
    net[j].nEns=nEns;
  }
  return(net);
}
//#endif //
FLOAT moverange(FLOAT y1,FLOAT y1min, FLOAT y1max, FLOAT y0min, FLOAT y0max)
{
//from y1 to y0
  FLOAT div=y1max-y1min;
  if(div>-1e-20 && div<1e-20) return((y0max+y0min)/2.);
  FLOAT y0;
  y0=((y1)-(y1min))*((y0max)-(y0min))/((y1max)-(y1min))+(y0min);
  if(y0>y0max) return(y0max);  else if(y0<y0min) return(y0min); else return(y0);
  //  if(y0>y0max) y0=y0max;  else if(y0<y0min) y0=y0min;  return(y0);
}
//
/* this file ends here,,, */
