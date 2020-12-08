/*
 * $id: time-stamp: "2004/02/25 05:10:46"$ 
 * "2004/02/13" kuro
 * sim.c
 *
 */
#include "randoms.h"
#include "sim.h"
//#define GNUPLOT "gnuplot -geometry 400x300"
#define GNUPLOT "gnuplot -geometry 320x240"
#define Infty 1e10
#define Degree 57.2957795131

//#define GSL 1
//#include "my_plinn.h"
//void calc_Ainvb(FLOAT *M, FLOAT *a_data, FLOAT *b_data, int nx, int ndata);
//#include "mytimer.c"
/*
typedef struct {
  FLOAT *MSEssp;
  FLOAT *MSEmsp;
  FLOAT *MSEtrain;
  FLOAT *MSE1;
  FLOAT *MSE2;
  FLOAT *Entro;
  int   *MSEtime;
  int   *c_times;
  int   *i_times;
} MSEbank;
*/
/*====================================================================*
 * disp_yhat
 * 予測値、予測誤差を画面/ファイルに保存する
 *====================================================================*/
void disp_yhat(NET *net,DATA *givendata, DATA *result, FILE *fpgd0, char *gpl, char *dat, char *obj,char *title, char *type,int disperr,char *fb)
{//yhat,ert,y vs x1x2
  int t;
  char str1[256],str2[256];
  FILE *fpd,*fpgd,*fpgf;
//  char fb[126];
//  sprintf(fb,"%s",gpl);
//  fb[strlen(gpl)-4]=0;

  //fprintf(fpgf,"set term postscript eps enhanced color;set output \"%s.eps\"\n",fb);

  if(givendata->data_class==TIME_SERIES){
    sprintf(str1, "set style data linespoints;set grid;set xlabel \"t\"\nset title \"%s(MSE%3.2eNMSE%3.2et%d) H%g(Ey%g)\";set term postscript eps enhanced color;set output \"%s.eps\"\n",title,result->MSE,result->NMSE,GlobalTime,net->tpH,net->tpEy,fb);
    //    sprintf(str1, "set style data linespoints\nset grid\nset xlabel \"t\"\nset title \"%s(MSE%3.2eNMSE%3.2et%d)\"\n",title,result->MSE,result->NMSE,GlobalTime);
    if(disperr){
      sprintf(str2,"plot \"%s\" using 4:3 t \"err[k]\"\n",dat);
    }
    else{
      sprintf(str2,"plot \"%s\" using 4:2 w l t \"y\", \"\" using 4:1 w l t \"yp\",\"\" using 4:($1-$2) w l t\"yp-y\" \n",dat);
      //      sprintf(str2,"plot \"%s\" using 4:2 t \"y[k]\", \"%s\" using 4:1 t \"y[k]^\"\n",dat, dat);
    }

    if(dat[0]=='-') fpd=fpgd=fpgd0;
    else {
      fpgd = open_pipe(GNUPLOT, "w");
      fpd  = open_file(dat, "w");
    }
    //display
    fprintf(fpgd, str1,NULL);
    fprintf(fpgd, str2,NULL);
    //data
    fprintf(fpd, "#yhat y error time (MSE%e,NMSE%e,time=%d)\n",result->MSE,result->NMSE,GlobalTime);
    if (PREDTRAIN==3){// Prediction
      for (t=1; t<givendata->n_total; t++){
	fprintf(fpd, "%e %e %e %d \n", result->Y[t],givendata->Y[t],result->e[t],t+givendata->t0);
      }
    }
    else if (PREDTRAIN){// Prediction
      for (t=1; t<givendata->n_train; t++){
	fprintf(fpd, "%e %e %e %d \n", result->Y[t],givendata->Y[t],result->e[t],t+givendata->t0);
      }
    }
    else if (strcmp(type, SSP) == 0 || (strcmp(type, MSP) == 0)){// Prediction
      for (t=givendata->n_train; t<givendata->n_total; t++){
	if(givendata->tp0==0) fprintf(fpd, "%e %e %e %d \n", result->Y[t],givendata->Y[t],result->e[t],t+givendata->t0);
	else fprintf(fpd, "%e %e %e %d \n", result->Y[t],givendata->Y[t],result->e[t],t-givendata->n_train+givendata->tp0);//
      }
    }
    if(dat[0]=='-') fprintf(fpd, "e\n");
    else close_file(fpd);

    //display the data
    fflush(fpgd);
    //close_file(fpg);//??
  }
  else if(givendata->data_class==FUNCTION_APPROX){
    FLOAT x0;
    if(net->xwidth>100) return;
    sprintf(str1,"set style data lines\n###set yrange [0:1]\n###set xrange [0:1]\nset style data lines\nset xlabel \"x0\"\nset ylabel \"x1\"\nset title \"%s MSE%3.2eNMSE%3.2et%d\"\nset term postscript eps enhanced color;set output \"%s.eps\"\n",title,result->MSE,result->NMSE,GlobalTime,fb);
    //sprintf(str1,"set style data lines\n###set yrange [0:1]\n###set xrange [0:1]\nset style data lines\nset xlabel \"x0\"\nset ylabel \"x1\"\nset terminal x11\nset title \"%s MSE%3.2eNMSE%3.2et%d\"\n",title,result->MSE,result->NMSE,GlobalTime);
    if(disperr){
    sprintf(str2,"splot '%s' using 1:2:4 t \"err\"\n",dat);
    }
    else{
    sprintf(str2,"splot '%s' using 1:2:5 t \"y\", '%s' using 1:2:3 t \"yhat\"\n",dat,dat);
    }
    if(dat[0]=='-') fpd=fpgd=fpgd0;
    else {
      fpgd =open_pipe(GNUPLOT, "w");
      fpd =open_file(dat,"w");
    }
    //display (command execution)
    fprintf(fpgd,str1,NULL);
    fprintf(fpgd,str2,NULL);
    //data
    fprintf(fpd, "#x1 x2 yhat error y(MSE=%e,NMSE=%e,time=%d)\n",result->MSE,result->NMSE,GlobalTime);
    x0=givendata->X[givendata->n_train][0];
    for (t=givendata->n_train; t<givendata->n_total; t++) {
      //      if(fabs(givendata->X[t][0]-givendata->X[t-1][0])>1e-10) fprintf(fpd, "\n");//kuro
      //      if(fabs(givendata->X[t][0]-givendata->X[t-1][0])>1e-3) fprintf(fpd, "\n");
      if(application==RANGE_APPLI){
	//	if(fabs(givendata->X[t][0]-givendata->X[t-1][0])>net->xwidth*0.6) fprintf(fpd, "\n");//for rangedata040513
	if(fabs(givendata->X[t][0]-x0)>1e-5) fprintf(fpd, "\n");x0=givendata->X[t][0];//kuro
      }
      //      if(fabs(givendata->X[t][0]<0 && givendata->X[t-1][0])>0) fprintf(fpd, "\n");//for rangedata040513
      fprintf(fpd,"%e %e %+e %+f %e\n",
	      givendata->X[t][0],givendata->X[t][1],result->Y[t],result->e[t],givendata->Y[t]);
    }
    if(dat[0]=='-') fprintf(fpd, "e\n");
    else close_file(fpd);
    //display the data
    fflush(fpgd);
    //close_file(fpg);//??
  }
  //file
  if(gpl[0]!=0){
    fpgf = open_file(gpl, "w");
    fprintf(fpgf,str1,NULL);//
    fprintf(fpgf,str2,NULL);
    //    fprintf(fpgf,"set term postscript eps enhanced color;set output \"%s.eps\"\n",fb);
      //    fprintf(fpgf, "set terminal tgif;set output \"%s\"\n", obj);
      //    fprintf(fpgf, "pause -1 \"Hit Return to Save %s\"\n", obj);
    fprintf(fpgf,str2,NULL);
    close_file(fpgf);//
    printf("-> Do 'gnuplot %s;gv %s.eps'. Data is saved in '%s'\n",gpl,fb,dat);
  }
  if(net->DISP!=0){
    char cmd[126];
    sprintf(cmd,"gv %s.eps&",fb);
    system(cmd);
  }
}

/*====================================================================*
 * disp_mse
 * 近似・予測曲線データを画面/ファイルに保存する
 *====================================================================*/
void disp_mse(NET *net,DATA *givendata, DATA *test, FILE *fpgd0, char *gpl, char *dat, char *obj,char *title, char *type, MSEbank msebank,char *fb)
{//yhat,ert,y vs x1x2
  //  int i_times=msebank.i_times[0];
  int c_times=msebank.c_times[0];
  FLOAT *MSEtrain=msebank.MSEtrain;
  FLOAT *MSEssp=msebank.MSEssp;
  FLOAT *MSEmsp=msebank.MSEmsp;
  int *MSEtime=msebank.MSEtime;
  FLOAT *MSE1=msebank.MSE1;
  FLOAT *MSE2=msebank.MSE2;
  FLOAT *Entro=msebank.Entro;
  int t;
  char str1[256],str2[256];
  FILE *fpd,*fpgd,*fpgf;
//  char fb[126];
//  sprintf(fb,"%s",gpl);
//  fb[strlen(gpl)-4]=0;

  sprintf(str1, "set title \"MSE,NMSE vs. t\"\nset style data linespoints\nset grid\nset logscale x\nset logscale y\n#set format y \"%s\"\nset term postscript eps enhanced color;set output \"%s.eps\"\n","%6.1e",fb);
  //  sprintf(str1, "set title \"MSE,NMSE vs. t\"\nset style data linespoints\nset grid\nset logscale x\nset logscale y\n#set format y \"%s\"\nset terminal x11\n","%6.1e");
  sprintf(str2,"plot \"%s\" using 1:2 t \"MSEssp\", \"%s\" using 1:3 t \"MSEmsp\", \"%s\" using 1:4 t \"NMSEssp\", \"%s\" using 1:5 t \"NMSEmsp\", \"%s\" using 1:6 t \"MSEtrain\", \"%s\" using 1:9 t \"Entro\"\n",dat, dat, dat, dat, dat,dat);
  
  if(dat[0]=='-') fpd=fpgd=fpgd0;
  else {
    fpgd = open_pipe(GNUPLOT, "w");
    fpd  = open_file(dat, "w");
  }
  //display
  fprintf(fpgd,str1,NULL);
  fprintf(fpgd,str2,NULL);
  //data
  fprintf(fpd, "#time MSEssp, MSEmsp, NMSEssp, NMSEmsp, MSEtrain,MSE1,MSE2,Entro,1-Entro\n");
  for (t=0; t<c_times; t++) {//      for (t=0; t<=c_times; t++) {
    fprintf(fpd, "%e %e %e %e %e %e %e %e %e %e\n", 
	    (FLOAT)MSEtime[t], 
	    MSEssp[t], MSEmsp[t],
	    MSEssp[t]/givendata->VARtest, //
	    MSEmsp[t]/givendata->VARtest,
	    MSEtrain[t],
	    MSE1[t],
	    MSE2[t],
	    Entro[t],
	    1.-Entro[t]
	    );
    //	    log(net->nentropy));
    //    printf("----------------log entropy=%e,%e,%e\n",log(net->nentropy),exp(net->nentropy),net->nentropy);
  }
  //data close
  if(dat[0]=='-') fprintf(fpd, "e\n");
  else close_file(fpd);
  //display the data
  fflush(fpgd);

  //file
  if(gpl[0]!=0){
    fpgf = open_file(gpl, "w");
    fprintf(fpgf,str1,NULL);//
    fprintf(fpgf,str2,NULL);
    fprintf(fpgf, "set terminal tgif\n");
    fprintf(fpgf, "set output \"%s\"\n", obj);
    fprintf(fpgf, "pause -1 \"Hit Return to Save %s\"\n", obj);
    fprintf(fpgf,str2,NULL);
    close_file(fpgf);//
    printf("-> Do 'gnuplot %s'. Data is saved in '%s'\n",gpl,dat);
    //    printf("Do 'gnuplot %s'.\n",gpl);
  }
  if(net->DISP!=0){
    char cmd[126];
    sprintf(cmd,"gv %s.eps&",fb);
    system(cmd);
  }
}
void disp_weights(NET *net, FILE *fpgd0, char *gpl, char *dat, char *obj,char *title, char *type,char *fb)
{//yhat,ert,y vs x1x2
  int i,j;
  char str1[256],str2[256];
  FILE *fpd,*fpgd,*fpgf;
//  char fb[126];
//  sprintf(fb,"%s",gpl);
//  fb[strlen(gpl)-4]=0;

  sprintf(str1, "set style data points\nset xlabel \"w0\"\nset ylabel \"w1\"\nset title \"%s time%d\"\nset term postscript eps enhanced color;set output \"%s.eps\"\n",title,GlobalTime,fb);
  if(net->k<=2){
    sprintf(str2,"plot \"%s\" using 1:2 t \"(w1,w2)\"\n",dat);
  }
  else{
    sprintf(str2,"splot \"%s\" using 1:2:3 t \"(w1,w2,w3)\"\n", dat);
  }

  if(dat[0]=='-') fpd=fpgd=fpgd0;
  else {
    fpgd = open_pipe(GNUPLOT, "w");
    fpd  = open_file(dat, "w");
  }
  //display
  fprintf(fpgd,str1,NULL);
  fprintf(fpgd,str2,NULL);
  //data
  fprintf(fpd, "#w1,w2,w3... (time=%d)\n",GlobalTime);
  for(i=0;i<net->n_cells;i++){
    for(j=0;j<net->k;j++){
      fprintf(fpd,"%+e ",net->cell[i].w[j]);
    }
    fprintf(fpd,"\n");
  }
  if(dat[0]=='-') fprintf(fpd, "e\n");
  else close_file(fpd);
  
  //display the data
  fflush(fpgd);
    //close_file(fpg);//??
  //file
  if(gpl[0]!=0){
    fpgf = open_file(gpl, "w");
    fprintf(fpgf,str1,NULL);//
    fprintf(fpgf,str2,NULL);
//    fprintf(fpgf, "set terminal tgif\n");
//    fprintf(fpgf, "set output \"%s\"\n", obj);
//    fprintf(fpgf, "pause -1 \"Hit Return to Save %s\"\n", obj);
    //    fprintf(fpgf,str2,NULL);
    close_file(fpgf);//
    printf("-> Do 'gnuplot %s;gv %s&'. Data is saved in '%s'\n",gpl,fb,dat);
  }
  if(net->DISP!=0){
    char cmd[126];
    sprintf(cmd,"gv %s.eps&",fb);
    system(cmd);
  }
}

/*====================================================================*
 * save_plot_weights
 * 荷重ベクトルをファイルに保存する
 *====================================================================*/
void save_plot_weights(NET *net, char *dat,char *type) {
  char name[32] = "save_plot_weights";
  FILE *fp = NULL;
  //  char type[64];
  //  int n_total = givendata->n_total;
  //  int n_channels = givendata->k;
  //  int n_train = givendata->n_train;
  int i = 0;

  printf("> debug: plot error data \"%s\". (%s)\n", dat, name);
  fp = open_file(dat, "w");
  for (i=0; i<net->n_cells; i++) {
    fprintf(fp, "%e %e\n", net->cell[i].w[0],net->cell[i].w[1]);
  }
  close_file(fp);
  return;
}
/*====================================================================*
 * exec_plot_weights
 * ファイルからデータを読みとりGNUPLOTで誤差結果をプロットする(時系列)
 *====================================================================*/
void exec_plot_weights(NET *net, char *dat, char *gpl, char *obj, char *title) {
  //char name[32] = "exec_plot_weights";
  FILE *gp = NULL, *fp = NULL;
  //  FLOAT MSE = result->MSE, NMSE = result->NMSE;
  //  char buff[128];

  // Display
  gp = open_pipe(GNUPLOT, "w");
  fprintf(gp, "set title \"Weights\"\n");
  fprintf(gp, "set style data points\n");
  fprintf(gp, "set xlabel \"w0\"\n");
  fprintf(gp, "set ylabel \"w1\"\n");
  //  fprintf(gp, "set yrange [%e:%e]\n", 0.,1.);
    fprintf(gp, "###set xrange [%e:%e]\n", 0.,1.);
  fprintf(gp, "plot \"%s\" using 1:2 t \"w\"\n",dat);
  fflush(gp);

  fp = open_file(gpl, "w");
  fprintf(fp, "set title \"Weights\"\n");
  fprintf(fp, "set style data points\n");
  fprintf(fp, "set xlabel \"w0\"\n");
  fprintf(fp, "set ylabel \"w1\"\n");
  //  fprintf(fp, "set yrange [%e:%e]\n", 0.,1.);
  //  fprintf(fp, "set xrange [%e:%e]\n", 0.,1.);
  fprintf(fp, "plot \"%s\" using 1:2 t \"w\" w lp\n",dat);
  // for Output X11
  fprintf(fp, "set terminal x11\n");
  fprintf(fp, "pause -1 \"Hit Return to Save %s\"\n", obj);
  // for Output Tgif
  fprintf(fp, "set terminal tgif\n");
  fprintf(fp, "set output \"%s\"\n", obj);
  //  fprintf(fp, buff);//  
  fprintf(fp, "plot \"%s\" using 1:2 t \"w\"\n",dat);
  fprintf(fp, "pause -1 \"Hit Return\"\n");
  close_file(fp);

  return;
}

/*====================================================================*
 * exec_plot
 * データをプロットする(時系列)
 *====================================================================*/
void exec_plot(DATA *givendata, DATA *result, char *title, char *type,NET *net) {
  //  char type[64];
  char dat[64], gpl[64], obj[64], fb[64];

  // Data Type
  //  strcpy(type, result->type);

  {
    sprintf(fb, "w");
    sprintf(dat, "w.dat");
    sprintf(gpl, "w.gpl");
    sprintf(obj, "w.obj");
    disp_weights(net, 0, gpl, dat, obj,"weights", "",fb);
  }
//  {
//    // Save And Plot Approximation/Prediction Error Data
//    sprintf(fb, "e_%s",type);
//    sprintf(dat, "y_%s.dat", type);
//    sprintf(gpl, "e_%s.gpl", type);
//    sprintf(obj, "e_%s.obj", type);
//    
//    disp_yhat(net,       //NET *net,
//	      givendata, //DATA *givendata, 
//	      result,      //DATA *result, 
//	      0,       //FILE *fpgd0,
//	      gpl,      //char *gpl,
//	      dat,       //char *dat, 
//	      obj,         //char *obj,
//	      "error",    //char *title, 
//	      SSP,        //char *type SSP,MSP or PREDTRAIN
//	      1,fb         //int disperr
//	      );
//  }
  {
    // Save And Plot Approximation/Prediction Curve Data
    sprintf(fb, "y_%s",type);
    sprintf(dat, "%s.dat", fb);
    sprintf(gpl, "%s.gpl", fb);
    sprintf(obj, "%s.obj", fb);
//    sprintf(dat, "y_%s.dat", type);
//    sprintf(gpl, "y_%s.gpl", type);
//    sprintf(obj, "y_%s.obj", type);
    disp_yhat(net,       //NET *net,
	      givendata, //DATA *givendata, 
	      result,      //DATA *result, 
	      0,       //FILE *fpgd0,
	      gpl,      //char *gpl,
	      dat,       //char *dat, 
	      obj,         //char *obj,
	      title,//	   "yhat",    //char *title, 
	      SSP,        //char *type SSP,MSP or PREDTRAIN
	      0,fb         //int disperr
	      );
  }
  {
    sprintf(fb,"tmp/mse_t");//type=ssp,msp,???
    sprintf(gpl, "%s.gpl", fb);
    FILE *fp=fopen(gpl,"w");
    fprintf(fp,"set style data lines;set title \"T:%d,%d N=%d seed=%lu\"\n",GlobalTime,net->Tpinv,net->n_cells,net->seed);
    fprintf(fp,"set term postscript eps enhanced color;set output \"%s.eps\"\n",fb);
    fprintf(fp,"set logscale y;set format y \"%%.0e\";set grid\n");
    fprintf(fp,"plot \"%s.csv\" using 0:1 w l t \"MSEtr\", \"\" using 0:2 w l t \"MSE\"\n",fb);
    fclose(fp);
    char cmd[126];
    sprintf(cmd,"gnuplot %s",gpl);
    system(cmd);
    if(net->DISP!=0){
      sprintf(cmd,"gv %s.eps&",fb);
      system(cmd);
    }

  }
  return;
}

/*====================================================================*
 * exec_aprx
 * 時系列関数の近似を行う
 *====================================================================*/
//void exec_aprx(NET *net, DATA *givendata, DATA *aprx) {
//  //char name[32] = "exec_aprx";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_total = givendata->n_total;
//  int n_channels = givendata->k;
//  int n_train = givendata->n_train;
//  //int n_pred = givendata->n;
//  int t=0;
//  int k=0;
//  FLOAT yr;
//  // 入力をコピー
//  for (t=0; t<(n_total); t++) {
//    for (k=0; k<=(n_channels); k++) {
//	aprx->x[t][k] = givendata->x[t][k];
//    }
//  }
//  aprx->ymax = givendata->ymax;
//  aprx->ymin = givendata->ymin;
//
//  // 出力および誤差を計算
//  for (t=0; t<n_train; t++) {
//    aprx->y[t] = calc_output(net, givendata->x[t],&yr);
//    aprx->e[t] = (aprx->y[t]-givendata->y[t]);
//    MSE += square(aprx->e[t]);
//  }
//  MSE /= n_train;
//  NMSE = MSE/(max-min);
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//  aprx->MSE = MSE;
//  aprx->NMSE = NMSE;
//  return;
//}

/*====================================================================*
 * exec_ssp_test(時系列+関数近似も)
 * 一段予測を行う
 *====================================================================*/
//void exec_ssp_test_rt_p2r(NET *net, DATA *givendata, DATA *test) {
//  //char name[32] = "exec_ssp_test";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_test = givendata->n_test;
//  FLOAT y;
//  FLOAT MSE1=0;
//  FLOAT MSE2=0;
//  FILE *fp;
//  int i;
//  char line[256];
//  int GasosuH,GasosuV;
//  FLOAT GakakuH,GakakuV,xc,yc,zc,roll,pitch,yaw;
//  FLOAT *r,*x1,*x2,dd;
//  int *dc;
//  int ViewV=60,ViewH=330;
//  printf("cos(theta) for search neighboring weight vector(2 for y0x):");
//  fgets(line,256,stdin);
//  sscanf(line,"%lf",&(net->cost));
//
//  printf("Camera parameters (GasosuH,GasosuV,GakakuH,GakakuV):");
//  fgets(line,80,stdin);
//  sscanf(line,"%d%d%lf%lf",&GasosuH,&GasosuV,&GakakuH,&GakakuV);
//
//  printf("Camera pose (xc,yc,zc,roll,pitch(tilt),yaw(pan):");
//  fgets(line,256,stdin);
//  sscanf(line,"%lf%lf%lf%lf%lf%lf",&xc,&yc,&zc,&roll,&pitch,&yaw);
//
//  fp=open_file("rangegot.dat","w");  //write predict.dat
//  {//higher resolution
//    FLOAT xxxx[3],yyyy,YYYY;
//    FLOAT cr,sr,cp,sp,cy,sy;
//    FLOAT JV,JH,JD,yw,Nx,Ny,Nz,Xwx,Xwy,Xwz,Xcwx,Xcwy,Xcwz,t;
//    FLOAT Zjx,Zjy,Zjz;
//    FLOAT Xpx,Xpy,Xpz,Xpr;
//    FLOAT Xpx1,Xpx2;
//    int jV,jH,ii;
//    int atta;
//
//    //    r=(FLOAT *)malloc(sizeof(FLOAT)*GasosuH*GasosuV);
//    r=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x1=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x2=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    dc=(int *)calloc(GasosuH*GasosuV,sizeof(int));
//    
//    xxxx[2]=1;
//    cr=cos(roll/Degree);
//    sr=sin(roll/Degree);
//    cp=cos(pitch/Degree);
//    sp=sin(pitch/Degree);
//    cy=cos(yaw/Degree);
//    sy=sin(yaw/Degree);
//    for(jV=0;jV<GasosuV;jV++){
//	JV=jV*GakakuV/GasosuV;
//	for(jH=0;jH<GasosuH;jH++){
//	  JH=jH*GakakuH/GasosuV;
//	  atta=0;
//
////	  Zjx = cp*JV*sr + cy*sp*sr - cr*sy + JH*(cr*cy + sp*sr*sy);
////	  Zjy = cp*cy - JV*sp + cp*JH*sy; 
////	  Zjz = cp*cr*JV + cr*cy*sp + sr*sy + JH*(-(cy*sr) + cr*sp*sy);
//	  Zjx=cp JV sr + JD (cy sp sr - cr sy) + JH (cr cy + sp sr sy); 
//	  Zjy=cp cy JD - JV sp + cp JH sy;
//	  Zjz=cp cr JV + JH (-(cy sr) + cr sp sy) + JD (cr cy sp + sr sy);
//	  for(i=0;i<net->n_cells;i++){
//	    Nx=-cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree) ;//Normal( Hosen) Vector
//	    Ny= cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree) ;
//	    Nz=-sin(net->cell[i].am.M[0][1]/Degree);
//	    
//	    yw =net->cell[i].am.M[0][2]
//	      +net->cell[i].am.M[0][0]*net->cell[i].w[0]
//	      +net->cell[i].am.M[0][1]*net->cell[i].w[1];
//	    Xwx= (yw*cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree)) ;
//	    Xwy= (yw*cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree)) ;
//	    Xwz= (yw*sin(net->cell[i].am.M[0][1]/Degree));
//	    Xcwx= xc-Xwx;
//	    Xcwy= yc-Xwy;
//	    Xcwz= zc-Xwz;
//
//	    dd=Zjx*Nx + Zjy*Ny + Zjz*Nz; if(fabs(dd)<1e-20) dd=1e-20;
//	    t=(-(Nx*Xcwx) - Ny*Xcwy - Nz*Xcwz)/dd;//Solve[{nn.(t Zj + Xcw) == 0}, t] 
//
//	    Xpx= t*Zjx+xc;
//	    Xpy= t*Zjy+yc;
//	    Xpz= t*Zjz+zc;
//	    //	  dd=((d=(Xpx-Xwx))*d)+((d=(Xpy-Xwy))*d)+((d=(Xpy-Xwy))*d);
//	    Xpr= sqrt(Xpx*Xpx+Xpy*Xpy+Xpz*Xpz);
//	    Xpx1=atan2(Xpx,Xpy)*Degree;
//	    Xpx2=atan2(Xpz,sqrt(Xpx*Xpx+Xpy*Xpy))*Degree;
//	    xxxx[0]=moverange(Xpx1, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
//	    xxxx[1]=moverange(Xpx2, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
//	    //{double d;
//	    //	  dd=sqrt(((d=(net->cell[i].am.M[0][0]-xxxx[0]))*d)+((d=(net->cell[i].am.M[0][1]-xxxx[1]))*d));
//	    //	  if(dd>0.2*net->xwidth) continue;
//	    //}
//	    if(fabs(net->cost)>1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
//	    else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
//	    if(jV==0 && jH == 16){
//	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx1,Xpx2,Xpr,YYYY,i,net->c);
//	      printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx,Xpy,Xpr,YYYY,i,net->c);
//	    }
//	    if(net->c==i){ 
////	      printf("cell %3d w=(%+6.3f,%+6.3f) for (%3d,%3d) xp=(%+6.3f,%+8.3f)<-Xp=(%+6.3f,%+6.3f)\n",
////		     i,net->cell[i].w[0],net->cell[i].w[1],
////		     jH,jV,xxxx[0],xxxx[1], Xpx1, Xpx2
////		     );
//	      atta=1;
//	      if(r[jV*GasosuH+jH]>1){
//		printf("(%d,%d) new r(%d)=%7.2f old r(%d)=%7.2f\n",jV,jH, i,YYYY, dc[jV*GasosuH+jH],r[jV*GasosuH+jH]);
//	      }
//	      if(r[jV*GasosuH+jH]<=1 || YYYY<r[jV*GasosuH+jH]){
//		r[jV*GasosuH+jH]=YYYY;
//		x1[jV*GasosuH+jH]=xxxx[0];
//		x2[jV*GasosuH+jH]=xxxx[1];
//		dc[jV*GasosuH+jH]=net->dc;
//	      }
//	      //	    break; 
//	    }
//	  }
//	  if(atta==0) {
//	    printf("Hendayo!??? nakatta (%d,%d)\n",jV,jH);
//	  }
//	}
//	fprintf(stdout,"."); fflush(stdout);
//    }
//    for(jV=0;jV<GasosuV;jV++){
//	for(jH=0;jH<GasosuH;jH++){
//	  fprintf(fp,"%d %d %13.7e", 
//		  jH,
//		  jV,
//		  r[jV*GasosuH+jH]);
//	  for(ii=-2;ii<=8;ii++){
//	    //	  if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
//	    if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %d",jH);
//	    else fprintf(fp," -");
//	  }
//	  fprintf(fp," %13.7e %13.7e", x1[jV*GasosuH+jH],x2[jV*GasosuH+jH]);
//	  fprintf(fp,"\n");
//	}
//	fprintf(fp,"\n");
//    }
//  }
//  fclose(fp);
//  
//  fp=open_file("rangegot.gpl","w");
//  
//  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
//  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"set term x11\nset view 60,330\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"set nocontour\nset data style lines\n");
//  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p 0 0, \"rangegot.dat\" using 1:3:2 w p 1 0\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"splot \"rangegot.dat\" using 1:3:2 w l 1 0\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"set contour\nset view 60,350\n");
//  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//  //  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3 w l\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:2 t \"no opposite weight (open boundary)\" w p 2 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:2 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:2 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:2 t \"xi(fu renzoku?) with y0x for t<0\"          w p 1 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:2 t \"xi(   renzoku?) with y0x for t>1\"          w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 10:2 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 11:2 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p 3 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:2 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p 3 8");
//  fprintf(fp,", \"weights.dat\" u 1:2 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l 3 1");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//
//  //////////////////////////////////////
//  {
//    int Hmin,Hmax,Vmin,Vmax;
//    for(;;){
//	printf("Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin=-1 for quit:");
//	fgets(line,256,stdin);
//	sscanf(line,"%d%d%d%d%d%d",&Hmin,&Hmax,&Vmin,&Vmax,&ViewH,&ViewV);
//	if(Hmin==Hmax) break;
//	if(Hmin<0) break;
//	if(Vmin<0) Vmin=0;
//	if(Hmax>GasosuH) Hmax=GasosuH;
//	if(Vmax>GasosuV) Vmax=GasosuV;
//	ViewH = (ViewH+360)%360;
//	ViewV =(ViewV+180)%180;
//
//	fp = open_pipe(GNUPLOT, "w");
//	//      fpd  = open_file("rangedisp.dat", "w");
//	
//	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
//	fprintf(fp,"set nocontour\nset data style lines\n");
//	fprintf(fp,"splot [%d:%d][][%d:%d] \"rangegot.dat\" using 1:3:2 w l 1 0\n",Hmin,Hmax,Vmin,Vmax);
//	//      fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//	//	for(jV=0;jV<GasosuV;jV++){
//	//	  for(jH=0;jH<GasosuH;jH++){
//	//	    fprintf(fp,"%d %d %13.7e", 
//	//		    jH,
//	//		    jV,
//	//		    r[jV*GasosuH+jH]);
//	//	    for(ii=-2;ii<=8;ii++){
//	//	      if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %e",jH);
//	//	      else fprintf(fp," -");
//	//	    }
//	//	    fprintf(fp,"\n");
//	//	  }
//	//	  fprintf(fp,"\n");
//	//	}
//	//	fprintf(fp, "e\n");
//	fflush(fp);
////
////	fp=open_file("rangedisp.gpl","w");
////	
////	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
////	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
////	fprintf(fp,"set nocontour\nset data style lines\n");
////	fprintf(fp,"splot [%e:%e][][%e:%e] \"rangegot.dat\" using 1:3:2 w p 1 0\n",Hmin,Hmax,Vmin,Vmax);
////	fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
////	fclose(fp);
////
////	system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//    }
//    fclose(fp);
//  }
//  //////////////////////////////////////
//
//  MSE2=MSE1=MSE /= n_test;
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//
//  printf("=>MSEtest=%f\n",MSE);
//
//  test->MSE2 =test->MSE1 =test->MSE  = MSE;
//  test->NMSE = NMSE;
//  free(r);
//  free(x1);
//  free(x2);
//  free(dc);
//  return;
//}

#define Polar2Polar      0
#define Polar2Rect       1
#define Rect2Polar       3 
#define Rect2Rect        4  

#define Polar       0  
#define Rect        1  

FLOAT moverange1(FLOAT y1,FLOAT y1min, FLOAT y1max, FLOAT y0min, FLOAT y0max)
{//from y1 to y0
  FLOAT div=y1max-y1min;
  if(div>-1e-20 && div<1e-20){
    return((y0max+y0min)/2.);
  }
  FLOAT y0;
  y0=((y1)-(y1min))*((y0max)-(y0min))/((y1max)-(y1min))+(y0min);
  //  if(y0>y0max) y0=y0max;
  //  else if(y0<y0min) y0=y0min;
  return(y0);
}

//void search_planes_2d_old(NET *net, DATA *givendata, DATA *test,FLOAT t1,int maxnp,int maxoptit)
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
//  //  FLOAT *Mmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//normal (hosen) vectors
//  FLOAT *_Mmean=(FLOAT*)malloc(sizeof(FLOAT)*k2);//normal (hosen) vectors
//  FLOAT *wmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//center of weight vectors
//  //  FLOAT *_nvmean=(FLOAT*)malloc(sizeof(FLOAT)*k2);//normal (hosen) vectors
//  FLOAT *_wmean=(FLOAT*)malloc(sizeof(FLOAT)*k2);//center of weight vectors
//  int *sameplanes=(int*)malloc(sizeof(int)*nEns*n_cells*n_cells);//normal (hosen) vectors
//  int *color=(int*)malloc(sizeof(int)*nEns*n_cells);//number of planes
//  int *ncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  int *sncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  FLOAT cost1,sint1;
//  char fname[20];
//  int col=0;
//  int it,i;
//  FILE *fp;
//  if(maxoptit<1) maxoptit=1;
//  //  if(t1>1 && t1<90) ;
//  if(t1>1e-5 && t1<90-1e-5) ;
//  else {
//    //    fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90] and number of planes to extract[<%d]:",n_cells);
//    //    scanf2("%lf%d",&t1,&maxnp);
//    fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90]:");
//    scanf1("%lf",&t1);
//  }
//  if(t1<1)t1=1; else if(t1>89) t1=89;
//  cost1=cos(t1*PI/180.);
//  sint1=sin(t1*PI/180.);
//
//  for(nn=0;nn<nEns;nn++){//calc normal vectors
//    //    FLOAT *_nvmeannn=&nvmean[(nn*n_cells+6)*k2]; FLOAT *nvmeannn=&_nvmeannn[(nn*6)*k2];    //
//    FLOAT *nvmeannn=&nvmean[(nn*n_cells)*k2];
//    FLOAT *nvnn=&nv[(nn*n_cells)*k2];
//    //    FLOAT *Mmeannn=&Mmean[(nn*n_cells)*k2];
//    FLOAT *wmeannn=&wmean[(nn*n_cells)*k2];
//    //    int *sameplanesnn=&sameplanes[(nn*n_cells)*n_cells];
//    int *colornn=&color[nn*n_cells];
//    int *ncolnn=&ncol[nn*n_cells];
//    //    int *_sncolnn=&sncol[nn*(n_cells+6)];int *sncolnn=&_sncolnn[nn*(6)];    //    
//    int *sncolnn=&sncol[nn*n_cells];
//    NET *netnn=&net[nn];
//    for(u1=0;u1<n_cells;u1++){
//      FLOAT sumM2=0;
//      for(j=0;j<k;j++) sumM2+=square(netnn->cell[u1].am.M[0][j]);
//      FLOAT nvk=nvnn[u1*k2+k]=sqrt(1.0/(1.0+sumM2));
//      for(j=0;j<k;j++)
//	nvnn[(u1)*k2+j]=-netnn->cell[u1].am.M[0][j]*nvk;
//      nvnn[(u1)*k2+k1]=netnn->cell[u1].am.M[0][k]*nvk;//distance from the origin;
//      colornn[u1]=0;//no color for each unit
//      ncolnn[u1]=0;//number of units for a color
//    }
//
//    for(col=1;col<n_cells;col++){
//      for(u1=0;u1<n_cells;u1++) if(colornn[u1]==0 && netnn->cell[u1].v>0) break;
//      if(u1>=n_cells-1) break;
//      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] = nvnn[(u1)*k2+j];
//      for(j=0;j<k  ;j++) wmeannn[col*k2+j] = netnn->cell[u1].w[j];
//
//      for(it=0;it<maxoptit;it++){
//	for(u2=0;u2<n_cells;u2++) if(colornn[u2]==col) colornn[u2]=0;
//	ncolnn[col]=0;
//	//	for(j=0;j<=k1;j++) _nvmean[j]=_wmean[j] =0;
//	for(j=0;j<=k1;j++) _Mmean[j]=_wmean[j] =0;
//	FLOAT vsum=0.0;
//
//	FLOAT y1=0, y2=0;
//	for(i=0;i<k;i++) y1 += (-nvmeannn[col*k2+i]/nvmeannn[col*k2+k])*wmeannn[col*k2+i];
//	y1+=(nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);
//	//
//
//	for(u2=0;u2<n_cells;u2++){
//	  if(colornn[u2]==0 && netnn->cell[u2].v>0){
//	    FLOAT n1n2=0.0;
//	    for(j=0;j<=k;j++) n1n2+=(nvmeannn[col*k2+j]*nvnn[(u2)*k2+j]);
//	    
//	    if(n1n2>=cost1){//condition #1
//	      FLOAT n1w12=0.0, n2w12=0.0, w12=0.0,w12j;
//	      for(j=0;j<k;j++){
//		w12j=wmeannn[col*k2+j]-netnn->cell[u2].w[j];
//		w12   += (w12j*w12j);
//		n1w12 += (nvmeannn[col*k2+j]*w12j);
//		n2w12 += (nvnn[(u2)*k2+j]*w12j);
//	      }
//	      {//for j=k
//		y2=0;
//		for(i=0;i<k;i++) y2 += (netnn->cell[u2].am.M[0][i]*netnn->cell[u2].w[i]);
//		y2 += (netnn->cell[u2].am.M[0][k]);//constant element
//		
//		w12j=y1-y2;
//		w12 += w12j*w12j;
//		n1w12 += (nvmeannn[col*k2+k]*w12j);
//		n2w12 += (nvnn[(u2)*k2+k]*w12j);
//	      }
//	      w12=sqrt(w12);
//	      //	      n1w12=fabs(n1w12/w12);
//	      //	      n2w12=fabs(n2w12/w12);
//		
//	      //	      if(fabs(n1w12/w12)<=sint1 && fabs(n2w12/w12)<=sint1){
//	      if(fabs(n1w12/w12)<=sint1 && fabs(n2w12/w12)<=sint1){//condition #2,#3
//		FLOAT y3=0;
//		if(it>10){//no use??
//		  for(i=0;i<k;i++) y3 += (nvmeannn[col*k2+i]/nvmeannn[col*k2+k])*netnn->cell[u2].w[i];
//		  y3 += (nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);//constant element
//		  if(fabs(y3-y2)>100) continue;
//		}
//
//		colornn[u2]=col;
//		ncolnn[col]++;
//		for(j=0;j<=k;j++) _Mmean[j]+= (netnn->cell[u2].am.M[0][j]*netnn->cell[u2].v);
//		for(j=0;j<k  ;j++) _wmean[j]+= (netnn->cell[u2].w[j]*netnn->cell[u2].v);
//		vsum+=netnn->cell[u2].v;
////		if(col==10){
////		  fprintf(stderr,"check2 M%+.4f %+.4f %+.4f M%+.4f %+.4f %+.4f w%+3.0f %+3.0f w%+3.0f %+3.0f v%f u%d nw=%+.2f %+.2f y=%3.0f %3.0f\n",
////			  _Mmean[0]/vsum,_Mmean[1]/vsum,_Mmean[2]/vsum,
////			  netnn->cell[u2].am.M[0][0],netnn->cell[u2].am.M[0][1],netnn->cell[u2].am.M[0][2],
////			  _wmean[0]/vsum,_wmean[1]/vsum,
////			  netnn->cell[u2].w[0],netnn->cell[u2].w[1],
////			  netnn->cell[u2].v,u2,
////			  n1w12/w12,n2w12/w12,y2,y3);
////		}
////
////		if(col==10 && (u2==380) ){//check
////		  fprintf(stderr,"check y1=%f y2=%f n1w12=%f n2w12=%f sint1=%f\n",y1,y2,n1w12,n2w12,sint1);
////		  //		  fprintf(stderr,"[%d] %f %f %d\n",u2,netnn->cell[u2].w[0],netnn->cell[u2].w[1]);
////		}
//	      }
//	    }
//	  }
//	}
//	if(ncolnn[col]<=0){
//	  ncolnn[col]=1;
//	  colornn[u1]=col;
//	  break;
//	}
//	else{
//	  //	  for(j=0;j<=k1;j++) nvmeannn[col*k2+j]= _nvmean[j];
//	  //	  for(j=0;j<=k  ;j++) _Mmean[j] /= ncolnn[col];
//	  for(j=0;j<=k  ;j++) _Mmean[j] /= vsum;
//	  FLOAT sumM2=0;
//	  for(j=0;j<k;j++) sumM2+=square(_Mmean[j]);
//	  FLOAT nvk=nvmeannn[col*k2+k]=sqrt(1.0/(1.0+sumM2));
//	  for(j=0;j<k;j++) nvmeannn[col*k2+j]=-_Mmean[j]*nvk;
//	  nvmeannn[col*k2+k1]=+_Mmean[k]*nvk;//distance from the origin;
//	  //	  for(j=0;j<k  ;j++) wmeannn[col*k2+j] =_wmean[j]/ncolnn[col]; 
//	  for(j=0;j<k  ;j++) wmeannn[col*k2+j] =_wmean[j]/vsum; 
//	}
//      }
//    }
//    int colmax=col;
//    if(maxnp>colmax) maxnp=colmax;
//    //    maxnp=colmax;//for check
//    fp=fopen("tmp/wmc.dat","w");
//    for(u1=0;u1<n_cells;u1++){
//      if(netnn->cell[u1].v<=0) continue;
//      for(j=0;j<k;j++) fprintf(fp,"%+e ",netnn->cell[u1].w[j]);
//      for(j=0;j<=k;j++) fprintf(fp,"%+e ",netnn->cell[u1].am.M[0][j]);
//      //      double y=netnn->cell[u1].am.M[0][k];
//      //      for(j=0;j<k;j++) y+=netnn->cell[u1].am.M[0][j]*netnn->cell[u1].w[j];
//      //      fprintf(fp,"%e ",y);
//      fprintf(fp,"%d\n",colornn[u1]);
//    }
//    fclose(fp);
//
//    fp=fopen("tmp/wmy.dat","w");
//    for(u1=0;u1<n_cells;u1++){
//      if(netnn->cell[u1].v<=0) continue;
//      for(j=0;j<k;j++) fprintf(fp,"%e ",netnn->cell[u1].w[j]);
//      //      for(j=0;j<=k;j++) fprintf(fp,"%e ",netnn->cell[u1].am.M[0][j]);
//      double y=netnn->cell[u1].am.M[0][k];
//      for(j=0;j<k;j++) y+=netnn->cell[u1].am.M[0][j]*netnn->cell[u1].w[j];
//      fprintf(fp,"%e ",y);
//      fprintf(fp,"%d\n",colornn[u1]);
//    }
//    fclose(fp);
//
//    //sorting below
//    for(i=0;i<colmax;i++) sncolnn[i]=i;//prepare
//    for(i=0;i<colmax-1;i++){//sort
//      for(j=i+1;j<colmax;j++){
//	if(ncolnn[sncolnn[j]]>ncolnn[sncolnn[i]]){
//	  int ntmp=sncolnn[j];
//	  sncolnn[j]=sncolnn[i];
//	  sncolnn[i]=ntmp;
//	}
//      }
//    } //sorting above
//
//    {
//      FILE *fpp=fopen("tmp/planes.dat","w");
//      col=colmax;
//      fprintf(fpp,"1 0 0 %f %d\n",net->xmin0[0],col++);
//      fprintf(fpp,"1 0 0 %f %d\n",net->xmax0[0],col++);
//      fprintf(fpp,"0 1 0 %f %d\n",net->xmin0[1],col++);
//      fprintf(fpp,"0 1 0 %f %d\n",net->xmax0[1],col++);
//      fprintf(fpp,"0 0 1 %f %d\n",net->ymin0,col++);
//      fprintf(fpp,"0 0 1 %f %d\n",net->ymax0,col++);
//
//	//
//      for(j=0;j<colmax;j++){//initialize
//	sprintf(fname,"tmp/wplane%d.dat",j);
//	fp=fopen(fname,"w");
//	fclose(fp);
//      }
//      FILE *fpw=fopen("tmp/wplanes.dat","w");
//      //      for(i=0;i<colmax;i++){//
//      for(i=0;i<maxnp;i++){//
//	col=sncolnn[i];
//	if(col==0) continue;//no color
//	//	fprintf(stderr,"Plane#%d (%d units;col#=%d) ",i,ncolnn[col],col);
//	//	fprintf(stderr,"Normal Vector (Hoh-Sen vector)=");
//	for(j=0;j<k1;j++){
//	  //	  fprintf(stderr,"%+.5f ",nvmeannn[col*k2+j]);
//	  fprintf(fpp,"%+.5f ",nvmeannn[col*k2+j]);
//	}
//	//	fprintf(stderr,"Distance from the origin (LRF)=%.5f\n",nvmeannn[col*k2+k1]);//??
//	fprintf(fpp,"%+.5f %d\n",nvmeannn[col*k2+k1],i);
//	
//	sprintf(fname,"tmp/wplane%d.dat",i);
//	fp=fopen(fname,"a+");
//	
//	for(u1=0;u1<n_cells;u1++){
//	  if(colornn[u1]==col && netnn->cell[u1].v>0){
//	    FLOAT *w=&(netnn->cell[u1].w[0]);
//	    //	    fprintf(stderr,"u%d",u1);//check
//	    for(j=0;j<k;j++) {
//	      fprintf(fp,"%e ",w[j]);
//	      fprintf(fpw,"%e ",w[j]);
//	    }
//	    
//	    FLOAT y=0;
//	    for(j=0;j<k;j++) y += (-nvmeannn[col*k2+j]/nvmeannn[col*k2+k])*w[j];
//	    y+=(nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);
//	    //	    if(y>1800 || y<400 || col==10){
//	    //	    if(y>1800 || y<400){
//	    //	      fprintf(stderr,"check y=%f w=%f %f u=%d\n",y,w[0],w[1],u1);
//	    //	    }
//	    if(1==1){//for estimated y by the CAN2
//	      FLOAT *x=(FLOAT*)malloc(sizeof(FLOAT)*k1),_yt,_y,_Yt;
//	      for(j=0;j<k;j++) x[j]=w[j];
//	      x[k]=1;
//	      calc_output(&net[nn], x,&_yt,&_y,&_Yt);
//	      y=_Yt;
//	    }
//	    
//	    fprintf(fp,"%e\n",y);
//	    fprintf(fpw,"%e %d\n",y,i);
//	    if(i==0 || col==10){//check
//	      fprintf(stderr,"\nw=");
//	      for(j=0;j<k;j++) fprintf(stderr,"%f ",netnn->cell[u1].w[j]);
//	      fprintf(stderr,"y=%f",y);
//	      
//	      //	      fprintf(stderr,"\nm=");
////	      for(j=0;j<k;j++) fprintf(stderr,"%f ",(-nvmean[col*k2+j]/nvmean[col*k2+k]));
////	      fprintf(stderr,"%f ",(nvmean[col*k2+k1]/nvmean[col*k2+k]));
//	      
//	      
//	      fprintf(stderr,"\nM=");
//	      for(j=0;j<=k;j++) fprintf(stderr,"%f ",netnn->cell[u1].am.M[0][j]);
//	      fprintf(stderr,"v=%f\nn=",netnn->cell[u1].v);
//	      
//	      for(j=0;j<k1;j++) fprintf(stderr,"%f ",nvnn[(u1)*k2+j]);
//	      fprintf(stderr,"d=%f",nvnn[(u1)*k2+k1]);
//	      fprintf(stderr,"\n");
//	    }
//	  }
//	}
//	fclose(fp);
//      }
//      fclose(fpw);
//      fp=fopen("tmp/wplanes.plt","w");
//      //      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
//      //      fprintf(fp,"set view 0,0; set pointsize 1.0;\n");
//      fprintf(fp,"set title \"Do 'cd tmp;gnuplot planes.plt' for review.\"\n");
//      //      fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set zlabel \"z\";set view 0,0; set pointsize 1\n");
//      fprintf(fp,"set xlabel \"x\"; set ylabel \"z\"; set zlabel \"y\";set view 60,350; set pointsize 1\n");
//      //      fprintf(fp,"splot [%f:%f][%f:%f][0:] \"../%s\" using 1:2:3 pointsize 0.1",
//      //      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1 linetype 0",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      fprintf(fp,"splot [%f:%f][][%f:%f] \"../%s\" using 1:3:2 pointsize 0.1 linetype 0",
//	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      for(i=0;i<colmax;i++){
//	col=sncolnn[i];
//	//	if(ncolnn[col]>=2){
//	if(i<maxnp){
//	  //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" w l %d",i,i,ncolnn[col],i+1);
//	  fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"wplane%d(%dunits)\" w l %d",i,i,ncolnn[col],i+1);
//	}
//	else if(ncolnn[col]<=1 || i>maxnp) {
//	  //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"\" pointsize 0.5 linetype 0",i);
//	  fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"\" pointsize 0.5 linetype 0",i);
//	}
//	else {
//	  //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" pointsize 0.5 linetype 0",i,i,ncolnn[col]);
//	  fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"wplane%d(%dunits)\" pointsize 0.5 linetype 0",i,i,ncolnn[col]);
//	}
//      }
//      fprintf(fp,"\npause -1 \"Do 'cd tmp;gnuplot wplanes.plt' for review. Return to continue.\"\n");
//
//      fprintf(fp,"set term tgif\n");
//      fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
//      //
////      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
////      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1",
////	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
////      for(i=0;i<colmax;i++){
////	col=sncolnn[i];
////	fprintf(fp,", \"wplane%d.dat\" using 1:2:3 w p",col);
////	if(ncolnn[col]<=1) fprintf(fp," 0 0 ps 0.1");
////      }
////      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
////      
////      fprintf(fp,"set termi tgif\n");
////      fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
//      fclose(fp);
//      system("cd tmp;xterm -geometry 80x5+0+0 -e gnuplot wplanes.plt&");
//      fprintf(stdout,"### Data of planes are saved in planes.dat\n");
//      fclose(fpp);
//    }
//    {
//      int t;
//      FLOAT _y,_yr,_Y;
//      for(j=0;j<colmax;j++){//initialize
//	sprintf(fname,"tmp/xplane%d.dat",j);
//	fp=fopen(fname,"w");
//	fclose(fp);
//      }
//      FILE *fpx=fopen("tmp/xplanes.dat","w");
//      for (t=givendata->n_train; t<givendata->n_total; t++) {
//	//      for (j=0; j<=k; j++) test->x[t][j] = givendata->x[t][j];
//	//      calc_output(netnn, test->x[t],&test->y[t],&y,&(test->Y[t]));
//	calc_output(netnn, givendata->x[t],&_y,&_yr,&_Y);
//	//	col=color[netnn->c];
//	for(i=0;i<colmax;i++) if(color[netnn->c]==sncolnn[i]) break;
//	sprintf(fname,"tmp/xplane%d.dat",i);
//	fp=fopen(fname,"a+");
//	for(j=0; j<k; j++){
//	  fprintf(fp,"%e ",test->x[t][j]);
//	  fprintf(fpx,"%e ",test->x[t][j]);
//	}
//	fprintf(fp,"%e\n",_Y);
//	fprintf(fpx,"%e %d\n",_Y,i);
//	fclose(fp);
//      }
//      fclose(fpx);
//      fp=fopen("tmp/planes.plt","w");
//      //      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
//      //      fprintf(fp,"set view 0,0; set pointsize 1.0;\n");
//      fprintf(fp,"set title \"Do 'cd tmp;gnuplot planes.plt' for review.\"\n");
//      //      fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set zlabel \"z\";set view 0,0; set pointsize 1\n");
//      fprintf(fp,"set xlabel \"x\"; set ylabel \"z\"; set zlabel \"y\";set view 60,350; set pointsize 1\n");
//      //      fprintf(fp,"splot [%f:%f][%f:%f][0:] \"../%s\" using 1:2:3 pointsize 0.1",
//      //      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1 linetype 0",
//      fprintf(fp,"splot [%f:%f][][%f:%f] \"../%s\" using 1:3:2 pointsize 0.1 linetype 0",
//	    //	    test->xmin0[0],test->xmax0[0],test->xmin0[1],test->xmax0[1],givendata->path);
//	    //	    netnn->xmin0[0],netnn->xmax0[0],netnn->xmin0[1],netnn->xmax0[1],givendata->path);
//	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      for(i=0;i<colmax;i++){
//	col=sncolnn[i];
//	//	if(ncolnn[col]>1) {
//	if(ncolnn[col]>ncolnn[sncolnn[maxnp]]) {
//	  //	  fprintf(fp,", \"xplane%d.dat\" using 1:2:3 t \"xplane%d(%dunits)\" w p %d",i,i,ncolnn[col],i+1);
//	  fprintf(fp,", \"xplane%d.dat\" using 1:3:2 t \"xplane%d(%dunits)\" w p %d",i,i,ncolnn[col],i+1);
//	  //	fprintf(fp,", \"plane%d.dat\" using 1:2:3 w p",col);
//	  //	  if(ncolnn[col]<=1) fprintf(fp," ps 0.1");
//	}
//	else{
//	  //	  fprintf(fp,", \"xplane%d.dat\" using 1:2:3 t \"\" pointsize 0.5 linetype 0",i);
//	  fprintf(fp,", \"xplane%d.dat\" using 1:3:2 t \"\" pointsize 0.5 linetype 0",i);
//	}
//      }
//      fprintf(fp,"\npause -1 \"Do 'cd tmp;gnuplot planes.plt' for review. Return to continue.\"\n");
//      
//      fprintf(fp,"set term tgif\n");
//      fprintf(fp,"set output \"planes.obj\"\nreplot\n");
//      fclose(fp);
//      system("cd tmp;xterm -geometry 80x5-0+0 -e gnuplot planes.plt&");
//    }
//  }
//  //  fclose(fp);
//  free(nv);
//  free(nvmean);
//  free(wmean);
//  //  free(_nvmean);
//  free(_Mmean);
//  free(_wmean);
//  free(sameplanes);
//  free(color);
//  free(ncol);
//  free(sncol);
//}
#define ZERO 1e-20
void search_planes_2d(NET *net, DATA *givendata, DATA *test,FLOAT t1,int maxnp,int maxoptit,int DISP)
{
  int j,nn;
  int u1,u2;
  int k=net->k;//input channels;
  int k1=k+1;
  int k2=k+2;
  int nEns=net[0].nEns;//number of ensemble
  int n_cells=net[0].n_cells;
  FLOAT *nv=(FLOAT *)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//normal (hosen) vectors
  FLOAT *nvmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//normal (hosen) vectors
  FLOAT *_Mmean=(FLOAT*)malloc(sizeof(FLOAT)*k2);//normal (hosen) vectors
  FLOAT *wmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//center of weight vectors
  FLOAT *_wmean=(FLOAT*)malloc(sizeof(FLOAT)*k2);//center of weight vectors
  int *sameplanes=(int*)malloc(sizeof(int)*nEns*n_cells*n_cells);//normal (hosen) vectors
  int *color=(int*)malloc(sizeof(int)*nEns*n_cells);//number of planes
  int *ncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
  int *sncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
  int *nxcol=(int*)malloc(sizeof(int)*nEns*n_cells);//number of data for col
  //  int *nxunt=(int*)calloc(nEns*n_cells,sizeof(int));//number of data for unit
  FLOAT cost1,sint1;
  char fname[20];
  int col=0;
  int it,i;
  FILE *fp;
  if(maxoptit<1) maxoptit=1;
  if(t1>1e-5 && t1<90-1e-5) ;
  else {
    fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90]:");
    scanf1("%lf",&t1);
  }
  if(t1<1)t1=1; else if(t1>89) t1=89;
  cost1=cos(t1*PI/180.);
  sint1=sin(t1*PI/180.);

  //  double *A=(double *)malloc(sizeof(double)*k1*n_cells);//k1=k+1
  //  double *b=(double *)malloc(sizeof(double)*n_cells);
  for(nn=0;nn<nEns;nn++){//calc normal vectors
    FLOAT *nvmeannn=&nvmean[(nn*n_cells)*k2];
    FLOAT *nvnn=&nv[(nn*n_cells)*k2];
    FLOAT *wmeannn=&wmean[(nn*n_cells)*k2];
    int *colornn=&color[nn*n_cells];
    int *ncolnn=&ncol[nn*n_cells];
    int *sncolnn=&sncol[nn*n_cells];
    int *nxcolnn=&nxcol[nn*n_cells];
    //    int *nxuntnn=&nxunt[nn*n_cells];
    NET *netnn=&net[nn];
    for(u1=0;u1<n_cells;u1++){
      FLOAT sumM2=0;for(j=0;j<k;j++) sumM2+=square(netnn->cell[u1].am.M[0][j]);
      //      FLOAT sumMw=0;for(j=0;j<k;j++) sumMw+=(netnn->cell[u1].am.M[0][j])*(netnn->cell[u1].am.w[0][j]);
      FLOAT nvk=nvnn[u1*k2+k]=1./sqrt(1.0+sumM2);
      if(netnn->cell[u1].am.M[0][k]<0) nvk*=-1.;//20111018 modified for slam LRF
      for(j=0;j<k;j++) nvnn[(u1)*k2+j]=-netnn->cell[u1].am.M[0][j]*nvk;
      nvnn[(u1)*k2+k1]=netnn->cell[u1].am.M[0][k]*nvk;//distance from the origin;
      colornn[u1]=0;//no color for each unit
      ncolnn[u1]=0;//number of units for a color
    }

    for(col=1;col<n_cells;col++){
      for(u1=0;u1<n_cells;u1++) if(colornn[u1]==0 && netnn->cell[u1].v>0) break;
      if(u1>=n_cells-1) break;
      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] = nvnn[(u1)*k2+j];
      for(j=0;j<k  ;j++) wmeannn[col*k2+j] = netnn->cell[u1].w[j];

      for(it=0;it<maxoptit;it++){
	AM am; init_AM(&am,k1,1);
	for(u2=0;u2<n_cells;u2++) if(colornn[u2]==col) colornn[u2]=0;
	ncolnn[col]=0;
	for(j=0;j<=k1;j++) _Mmean[j]=_wmean[j] =0;
	FLOAT vsum=0.0;

	FLOAT y1=0, y2=0;
	for(i=0;i<k;i++) y1 += (-nvmeannn[col*k2+i]/nvmeannn[col*k2+k])*wmeannn[col*k2+i];
	y1+=(nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);

	for(u2=0;u2<n_cells;u2++){
	  if(colornn[u2]==0 && netnn->cell[u2].v>0){
	    FLOAT n1n2=0.0;
	    for(j=0;j<=k;j++) n1n2+=(nvmeannn[col*k2+j]*nvnn[(u2)*k2+j]);//
//	    if(n1n2<=-cost1){
//	      fprintf(stderr,"check cost1\n");
//	    }
//	    if(u2==9){
//	      int check;check=1;
//	    }
//	    if(n1n2>=cost1){//condition #1
	    if(n1n2>=cost1 || n1n2<=-cost1){//condition #1
	      FLOAT n1w12=0.0, n2w12=0.0, w12=0.0,w12j;
	      for(j=0;j<k;j++){
		w12j=wmeannn[col*k2+j]-netnn->cell[u2].w[j];
		w12   += (w12j*w12j);
		n1w12 += (nvmeannn[col*k2+j]*w12j);
		n2w12 += (nvnn[(u2)*k2+j]*w12j);
	      }
	      {//for j=k
		y2=0;
		for(i=0;i<k;i++) y2 += (netnn->cell[u2].am.M[0][i]*netnn->cell[u2].w[i]);
		y2 += (netnn->cell[u2].am.M[0][k]);//constant element
		
		w12j=y1-y2;
		w12 += w12j*w12j;
		n1w12 += (nvmeannn[col*k2+k]*w12j);
		n2w12 += (nvnn[(u2)*k2+k]*w12j);
	      }
	      w12=sqrt(w12);
	      //	      if(fabs(n1w12/w12)<=sint1 && fabs(n2w12/w12)<=sint1){//condition #2,#3
	      double resolution=30.;
	      if(fabs(w12)<ZERO||//added 090103
		 fabs(n1w12)<resolution){//condition #2+#3+resolution condition 081230
		 //(fabs(n1w12/w12)<=sint1&&fabs(n2w12/w12)<=sint1&&fabs(n1w12)<resolution)){//condition #2+#3+resolution condition 081230
		 //		 (fabs(n1w12/w12)<=sint1&&fabs(n2w12/w12)<=sint1)){//condition #2+#3+resolution condition 081230
		FLOAT y3=0;
		if(it>1000){//no use??//		if(it>10){//no use??
		  for(i=0;i<k;i++) y3 += (nvmeannn[col*k2+i]/nvmeannn[col*k2+k])*netnn->cell[u2].w[i];
		  y3 += (nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);//constant element
		  if(fabs(y3-y2)>100) continue;
		}
		{
		  am.y[0]=y2; am.x[k]=1.0;
		  for(j=0;j<k;j++) am.x[j]=netnn->cell[u2].w[j];
		  calc_AM(am);
		}

		for(j=0;j<=k;j++) _Mmean[j]+= (netnn->cell[u2].am.M[0][j]*netnn->cell[u2].v);
		for(j=0;j<k  ;j++) _wmean[j]+= (netnn->cell[u2].w[j]*netnn->cell[u2].v);
		vsum+=netnn->cell[u2].v;
		colornn[u2]=col;
		ncolnn[col]++;
	      }
	    }
	  }
	}
	if(ncolnn[col]<=0){
	  ncolnn[col]=1;
	  colornn[u1]=col;
	  break;
	}
	else{
	  //	  if(ncolnn[col]>=k1){//not enough
	  if(ncolnn[col]>k1){//new version
	  //	  if(1==0){//for old version
	    for(j=0;j<=k;j++) _Mmean[j] = am.M[0][j];
	  }
	  else{
	    for(j=0;j<=k;j++) _Mmean[j] /= vsum;
	  }
	  {//??work for k=2 but doesnot for k=1?? why?
	    FLOAT sumM2=0;
	    for(j=0;j<k;j++) sumM2+=square(_Mmean[j]);
	    FLOAT nvk=nvmeannn[col*k2+k]=1./sqrt(1.0+sumM2);
	    if(_Mmean[k]<0) _Mmean[k]*=-1;//20111018 modified for slam LFR
	    for(j=0;j<k;j++) nvmeannn[col*k2+j]=-_Mmean[j]*nvk;
	    nvmeannn[col*k2+k1]=+_Mmean[k]*nvk;//distance from the origin;
	    for(j=0;j<k  ;j++) wmeannn[col*k2+j] =_wmean[j]/vsum; 
	    //	    fprintf(stderr,"it%d col%d wmean=%g\n",it,col,_wmean[0]);//for check
	  }
	}
	free_AM(&am);
      }
    }
    int colmax=col;
    if(maxnp>colmax) maxnp=colmax;


    //sorting below
    for(i=0;i<colmax;i++){
      sncolnn[i]=i;//prepare
      nxcolnn[i]=0;
    }
    for(i=0;i<colmax-1;i++){//sort
      for(j=i+1;j<colmax;j++){
	if(ncolnn[sncolnn[j]]>ncolnn[sncolnn[i]]){
	  int ntmp=sncolnn[j];
	  sncolnn[j]=sncolnn[i];
	  sncolnn[i]=ntmp;
	}
      }
    } //sorting above

    {
      int t;
      FLOAT _y,_yr,_Y;
      for(j=0;j<colmax;j++){//initialize
	sprintf(fname,"tmp/xplane%d.dat",j);
	fp=fopen(fname,"w");
	fclose(fp);
      }
      FILE *fpx=fopen("tmp/xplanes.dat","w");
      int t0=0;
      if(givendata->n_total-givendata->n_train>0) t0=givendata->n_train;
      int t1=givendata->n_total;
      for (t=t0; t<t1; t++) {
	//      for (t=givendata->n_train; t<givendata->n_total; t++) {
	//      for (j=0; j<=k; j++) test->x[t][j] = givendata->x[t][j];
	//      calc_output(netnn, test->x[t],&test->y[t],&y,&(test->Y[t]));
	calc_output(netnn, givendata->x[t],&_y,&_yr,&_Y);
	//	col=color[netnn->c];
	////	nxuntnn[netnn->c]++;
	for(i=0;i<colmax;i++) if(color[netnn->c]==sncolnn[i]) break;
	sprintf(fname,"tmp/xplane%d.dat",i);
	fp=fopen(fname,"a+");
	for(j=0; j<k; j++){
	  fprintf(fp,"%e ",test->x[t][j]);
	  fprintf(fpx,"%e ",test->x[t][j]);
	}
	fprintf(fp,"%e\n",_Y);
	fprintf(fpx,"%e %d %d\n",_Y,i,netnn->c);
	fclose(fp);
      }
      fclose(fpx);
    }
    {    //
      fp=fopen("tmp/wmc.dat","w");
      fprintf(fp,"#w0,w1,M1,M2,M0,PlaneNo.,NumDdatainVoronoi\n");
      for(u1=0;u1<n_cells;u1++){
	//	if(netnn->cell[u1].v<=0) continue;//commentout on 2001.02.19 for can2pfrir.c
	for(j=0;j<k;j++) fprintf(fp,"%+e ",netnn->cell[u1].w[j]);
	for(j=0;j<=k;j++) fprintf(fp,"%+e ",netnn->cell[u1].am.M[0][j]);
	int col=colornn[u1];
	int nx=netnn->cell[u1].v*netnn->vmax;
	if(nx>0) nxcolnn[col]+=nx;
	fprintf(fp,"%d %d\n",colornn[u1],nx);
	//	fprintf(fp,"%d %g\n",colornn[u1],netnn->cell[u1].v*netnn->vmax);
	//fprintf(fp,"%d %d %g\n",colornn[u1],nxuntnn[u1],netnn->cell[u1].v*netnn->vmax);
      }
      fclose(fp);
//      fp=fopen("tmp/wmy.dat","w");
//      for(u1=0;u1<n_cells;u1++){
//	if(netnn->cell[u1].v<=0) continue;
//	for(j=0;j<k;j++) fprintf(fp,"%e ",netnn->cell[u1].w[j]);
//	double y=netnn->cell[u1].am.M[0][k];
//	for(j=0;j<k;j++) y+=netnn->cell[u1].am.M[0][j]*netnn->cell[u1].w[j];
//	fprintf(fp,"%e ",y);
//	fprintf(fp,"%d %d\n",colornn[u1],netnn->cell[u1].v*netnn->vmax);
//      }
//      fclose(fp);
    }
    {
      FILE *fpp=fopen("tmp/planes.dat","w");
      col=colmax;
      for(j=0;j<k;j++){
	for(i=0;i<=k;i++){if(i==j) fprintf(fpp,"1 "); else fprintf(fpp,"0 ");}
	fprintf(fpp,"%f %d 1 0 #nx ny d colnumber numberofw's\n",net->xmin0[j],col++);
	for(i=0;i<=k;i++){if(i==j) fprintf(fpp,"1 "); else fprintf(fpp,"0 ");}
	fprintf(fpp,"%f %d 1 0\n",net->xmax0[j],col++);
      }
      for(i=0;i<k;i++) fprintf(fpp,"0 ");
      fprintf(fpp,"1 %f %d 1 0\n",net->ymin0,col++);
      for(i=0;i<k;i++) fprintf(fpp,"0 ");
      fprintf(fpp,"1 %f %d 1 0\n",net->ymax0,col++);

      for(j=0;j<colmax;j++){//initialize
	sprintf(fname,"tmp/wplane%d.dat",j);
	fp=fopen(fname,"w");
	fclose(fp);
      }
      FILE *fpw=fopen("tmp/wplanes.dat","w");
      //      for(i=0;i<colmax;i++){//
      for(i=0;i<maxnp;i++){//
	col=sncolnn[i];
	if(col==0) continue;//no color
	//	fprintf(stderr,"Plane#%d (%d units;col#=%d) ",i,ncolnn[col],col);
	//	fprintf(stderr,"Normal Vector (Hoh-Sen vector)=");
	for(j=0;j<k1;j++){
	  //	  fprintf(stderr,"%+.5f ",nvmeannn[col*k2+j]);
	  fprintf(fpp,"%+.5f ",nvmeannn[col*k2+j]);
	}
	//	fprintf(stderr,"Distance from the origin (LRF)=%.5f\n",nvmeannn[col*k2+k1]);//??
	fprintf(fpp,"%+.5f %d %d %d %d\n",nvmeannn[col*k2+k1],i,ncolnn[col],nxcolnn[col],col);
	//	fprintf(fpp,"%+.5f %d %d %d\n",nvmeannn[col*k2+k1],i,ncolnn[col],nxcolnn[i]);
	
	sprintf(fname,"tmp/wplane%d.dat",i);//	sprintf(fname,"tmp/wplane%d.dat",i);
	fp=fopen(fname,"a+");
	
	for(u1=0;u1<n_cells;u1++){
	  if(colornn[u1]==col && netnn->cell[u1].v>0){
	    FLOAT *w=&(netnn->cell[u1].w[0]);
	    //	    fprintf(stderr,"u%d",u1);//check
	    for(j=0;j<k;j++) {
	      fprintf(fp,"%e ",w[j]);
	      fprintf(fpw,"%e ",w[j]);
	    }
	    
	    FLOAT y=0;
	    for(j=0;j<k;j++) y += (-nvmeannn[col*k2+j]/nvmeannn[col*k2+k])*w[j];
	    y+=(nvmeannn[col*k2+k1]/nvmeannn[col*k2+k]);
	    if(1==1){//for estimated y by the CAN2
	      FLOAT *x=(FLOAT*)malloc(sizeof(FLOAT)*k1),_yt,_y,_Yt;
	      for(j=0;j<k;j++) x[j]=w[j];
	      x[k]=1;
	      calc_output(&net[nn], x,&_yt,&_y,&_Yt);
	      y=_Yt;
	    }
	    fprintf(fp,"%e\n",y);
	    fprintf(fpw,"%e %d\n",y,i);
//	    if(i==0 || col==10){//check
//	      //	      fprintf(stderr,"\nw=");
//	      //	      for(j=0;j<k;j++) fprintf(stderr,"%f ",netnn->cell[u1].w[j]);
//	      //	      fprintf(stderr,"y=%f",y);
//	      
//	      //	      fprintf(stderr,"\nm=");
////	      for(j=0;j<k;j++) fprintf(stderr,"%f ",(-nvmean[col*k2+j]/nvmean[col*k2+k]));
////	      fprintf(stderr,"%f ",(nvmean[col*k2+k1]/nvmean[col*k2+k]));
//	      
//	      
////	      fprintf(stderr,"\nM=");
////	      for(j=0;j<=k;j++) fprintf(stderr,"%f ",netnn->cell[u1].am.M[0][j]);
////	      fprintf(stderr,"v=%f\nn=",netnn->cell[u1].v);
//	      
////	      for(j=0;j<k1;j++) fprintf(stderr,"%f ",nvnn[(u1)*k2+j]);
//	      //	      fprintf(stderr,"d=%f",nvnn[(u1)*k2+k1]);
//	      //	      fprintf(stderr,"\n");
//	    }

	  }
	}
	fclose(fp);
      }
      fclose(fpw);
      fclose(fpp);
      //fprintf(stderr,"#laptime:search_planes_2d=%5.3f\n",mytimer_lap());//check
      fprintf(stderr,"#time elapsed:search_planes_2d=%5.3f\n",mytimer_total());//check
      if(DISP==1){
	if(k==2){
	  fp=fopen("tmp/wplanes.plt","w");
	  //      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
	  //      fprintf(fp,"set view 0,0; set pointsize 1.0;\n");
	  fprintf(fp,"set title \"Do 'cd tmp;gnuplot planes.plt' for review.\"\n");
	  //      fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set zlabel \"z\";set view 0,0; set pointsize 1\n");
	  fprintf(fp,"set xlabel \"x\"; set ylabel \"z\"; set zlabel \"y\";set view 60,350; set pointsize 1\n");
	  //      fprintf(fp,"splot [%f:%f][%f:%f][0:] \"../%s\" using 1:2:3 pointsize 0.1",
	  //      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1 linetype 0",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  fprintf(fp,"splot [%f:%f][][%f:%f] \"../%s\" using 1:3:2 pointsize 0.1 linetype 0",
		  net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  for(i=0;i<colmax;i++){
	    col=sncolnn[i];
	    //	if(ncolnn[col]>=2){
	    if(i<maxnp){
	      //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" w l %d",i,i,ncolnn[col],i+1);
	      fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"wplane%d(%dunits)\" w l %d",i,i,ncolnn[col],i+1);
	    }
	    else if(ncolnn[col]<=1 || i>maxnp) {
	      //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"\" pointsize 0.5 linetype 0",i);
	      fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"\" pointsize 0.5 linetype 0",i);
	    }
	    else {
	      //	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" pointsize 0.5 linetype 0",i,i,ncolnn[col]);
	      fprintf(fp,", \"wplane%d.dat\" using 1:3:2 t \"wplane%d(%dunits)\" pointsize 0.5 linetype 0",i,i,ncolnn[col]);
	    }
	  }
	  fprintf(fp,"\npause -1 \"Do 'cd tmp;gnuplot wplanes.plt' for review. Return to continue.\"\n");
	  fprintf(fp,"set term tgif\n");
	  fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
	  //
	  //      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
	  //      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1",
	  //	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  //      for(i=0;i<colmax;i++){
	  //	col=sncolnn[i];
	  //	fprintf(fp,", \"wplane%d.dat\" using 1:2:3 w p",col);
	  //	if(ncolnn[col]<=1) fprintf(fp," 0 0 ps 0.1");
	  //      }
	  //      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
	  //      
	  //      fprintf(fp,"set termi tgif\n");
	  //      fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
	  fclose(fp);
	  system("cd tmp;xterm -geometry 80x5+0+0 -e gnuplot wplanes.plt&");
	  fprintf(stdout,"### Data of planes are saved in planes.dat\n");
	}//closing if(k==2)
	else if(k==1){
	  fp=fopen("tmp/wplanes.plt","w");
	  fprintf(fp,"set title \"Do 'gnuplot wplanes.plt' for review.\"\n");
	  fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set pointsize 1\n");
	  //	  fprintf(fp,"plot [%f:%f][%f:%f] \"%s\" using 1:2 pointsize 0.1 linetype 0",
	  //		  net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  fprintf(fp,"plot \"%s\" using 1:2 pointsize 0.1 linetype 0",
		  givendata->path);
	  for(i=0;i<colmax;i++){
	    col=sncolnn[i];
	    if(i<maxnp && ncolnn[col]>0){
	      fprintf(fp,", \"tmp/wplane%d.dat\" using 1:2 t \"wplane%d(%dunits)\" w p pt %d",i,i,ncolnn[col],i+1);
	    }
//	    else if(ncolnn[col]<=1 || i>maxnp) {
//	      fprintf(fp,", \"tmp/wplane%d.dat\" using 1:2 t \"\" pointsize 0.5 linetype 0",i);
//	    }
//	    else {
//	      fprintf(fp,", \"tmp/wplane%d.dat\" using 1:2 t \"wplane%d(%dunits)\" pointsize 0.5 linetype 0",i,i,ncolnn[col]);
//	    }
	  }
	  fprintf(fp,"\npause -1 \"Do 'gnuplot wplanes.plt' for review. Return to continue.\"\n");
	  fprintf(fp,"set term tgif\n");
	  fprintf(fp,"set output \"tmp/wplanes.obj\"\nreplot\n");
	  fclose(fp);
	  system("xterm -geometry 80x5+0+0 -e gnuplot tmp/wplanes.plt&");
	  fprintf(stdout,"### Data of planes are saved in planes.dat\n");
	}//closing if(k==1)
      }
    }
    {
#ifdef ORIGxplane
      int t;
      FLOAT _y,_yr,_Y;
      for(j=0;j<colmax;j++){//initialize
	sprintf(fname,"tmp/xplane%d.dat",j);
	fp=fopen(fname,"w");
	fclose(fp);
      }
      FILE *fpx=fopen("tmp/xplanes.dat","w");
      for (t=givendata->n_train; t<givendata->n_total; t++) {
	//      for (j=0; j<=k; j++) test->x[t][j] = givendata->x[t][j];
	//      calc_output(netnn, test->x[t],&test->y[t],&y,&(test->Y[t]));
	calc_output(netnn, givendata->x[t],&_y,&_yr,&_Y);
	//	col=color[netnn->c];
	for(i=0;i<colmax;i++) if(color[netnn->c]==sncolnn[i]) break;
	sprintf(fname,"tmp/xplane%d.dat",i);
	fp=fopen(fname,"a+");
	for(j=0; j<k; j++){
	  fprintf(fp,"%e ",test->x[t][j]);
	  fprintf(fpx,"%e ",test->x[t][j]);
	}
	fprintf(fp,"%e\n",_Y);
	fprintf(fpx,"%e %d %d\n",_Y,i,netnn->c);
	fclose(fp);
      }
      fclose(fpx);
#endif//#ifdef ORIGxplane

      if(DISP==1){
	if(k==2){
	  fp=fopen("tmp/planes.plt","w");
	  //      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
	  //      fprintf(fp,"set view 0,0; set pointsize 1.0;\n");
	  fprintf(fp,"set title \"Do 'cd tmp;gnuplot planes.plt' for review.\"\n");
	  //      fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set zlabel \"z\";set view 0,0; set pointsize 1\n");
	  fprintf(fp,"set xlabel \"x\"; set ylabel \"z\"; set zlabel \"y\";set view 60,350; set pointsize 1\n");
	  //      fprintf(fp,"splot [%f:%f][%f:%f][0:] \"../%s\" using 1:2:3 pointsize 0.1",
	  //      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1 linetype 0",
	  fprintf(fp,"splot [%f:%f][][%f:%f] \"../%s\" using 1:3:2 pointsize 0.1 linetype 0",
		  //	    test->xmin0[0],test->xmax0[0],test->xmin0[1],test->xmax0[1],givendata->path);
		  //	    netnn->xmin0[0],netnn->xmax0[0],netnn->xmin0[1],netnn->xmax0[1],givendata->path);
		  net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  for(i=0;i<colmax;i++){
	    col=sncolnn[i];
	    //	if(ncolnn[col]>1) {
	    if(ncolnn[col]>ncolnn[sncolnn[maxnp]]) {
	      //	  fprintf(fp,", \"xplane%d.dat\" using 1:2:3 t \"xplane%d(%dunits)\" w p %d",i,i,ncolnn[col],i+1);
	      fprintf(fp,", \"xplane%d.dat\" using 1:3:2 t \"xplane%d(%dunits,%ddata)\" w p pt %d",i,i,ncolnn[col],nxcolnn[i],i+1);
	      //	fprintf(fp,", \"plane%d.dat\" using 1:2:3 w p",col);
	      //	  if(ncolnn[col]<=1) fprintf(fp," ps 0.1");
	    }
	    else{
	      //	  fprintf(fp,", \"xplane%d.dat\" using 1:2:3 t \"\" pointsize 0.5 linetype 0",i);
	      fprintf(fp,", \"xplane%d.dat\" using 1:3:2 t \"\" pointsize 0.5 linetype 0",i);
	    }
	  }
	  fprintf(fp,"\npause -1 \"Do 'cd tmp;gnuplot planes.plt' for review. Return to continue.\"\n");
	  
	  fprintf(fp,"set term tgif\n");
	  fprintf(fp,"set output \"planes.obj\"\nreplot\n");
	  fclose(fp);
	  system("cd tmp;xterm -geometry 80x5-0+0 -e gnuplot planes.plt&");
	}//closing if(k==2)
	else if(k==1){
	  fp=fopen("tmp/planes.plt","w");
	  fprintf(fp,"set title \"Do 'gnuplot planes.plt' for review.\"\n");
	  fprintf(fp,"set xlabel \"x\"; set ylabel \"y\"; set pointsize 1\n");
	  //	  fprintf(fp,"plot [%f:%f][%f:%f] \"%s\" using 1:2 pointsize 0.1 linetype 0",
	  //		  net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
	  fprintf(fp,"plot \"%s\" using 1:2 pointsize 0.1 linetype 0",givendata->path);
	  for(i=0;i<colmax;i++){
	    col=sncolnn[i];
	    //	    if(ncolnn[col]>ncolnn[sncolnn[maxnp]]) {
	    if(ncolnn[col]>0) {
	      fprintf(fp,", \"tmp/xplane%d.dat\" using 1:2 t \"xplane%d(%dunits,%ddata)\" w p pt %d",i,i,ncolnn[col],nxcolnn[i],i+1);
	    }
//	    else{
//	      fprintf(fp,", \"xplane%d.dat\" using 1:2 t \"\" pointsize 0.5 linetype 0",i);
//	    }
	  }
	  fprintf(fp,"\npause -1 \"Do 'gnuplot planes.plt' for review. Return to continue.\"\n");
	  fprintf(fp,"set term tgif\n");
	  fprintf(fp,"set output \"planes.obj\"\nreplot\n");
	  fclose(fp);
	  system("xterm -geometry 80x5-0+0 -e gnuplot tmp/planes.plt&");
	}//closing if(k==1)
      }
    }
  }
  //  fprintf(stderr,"#laptime:search_planes_2dend=%5.3f\n",mytimer_lap());//check
  //  fclose(fp);
  free(nv);
  free(nvmean);
  free(wmean);
  //  free(_nvmean);
  free(_Mmean);
  free(_wmean);
  free(sameplanes);
  free(color);
  free(ncol);
  free(sncol);
}
//void search_planes_2dold(NET *net, DATA *givendata, DATA *test,FLOAT t1,int maxnp)
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
//  FLOAT *wmean=(FLOAT*)malloc(sizeof(FLOAT)*nEns*n_cells*k2);//center of weight vectors
//  int *sameplanes=(int*)malloc(sizeof(int)*nEns*n_cells*n_cells);//normal (hosen) vectors
//  int *color=(int*)malloc(sizeof(int)*nEns*n_cells);//number of planes
//  int *ncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  int *sncol=(int*)malloc(sizeof(int)*nEns*n_cells);//sort of planes
//  FLOAT nvk1;
//  FLOAT cost1,sint1;
//  char fname[20];
//  int col=0;
//  if(t1>1 && t1<90) ;
//  else {
//    //    fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90] and number of planes to extract[<%d]:",n_cells);
//    //    scanf2("%lf%d",&t1,&maxnp);
//    fprintf(stdout,"Allowable angle (degree) for eauivalent plane[0-90]:");
//    scanf1("%lf",&t1);
//  }
//  if(t1<1)t1=1; else if(t1>89) t1=89;
//  cost1=cos(t1*PI/180.);
//  sint1=sin(t1*PI/180.);
//
//  for(nn=0;nn<nEns;nn++){//calc normal vectors
//    FLOAT *nvnn=&nv[(nn*n_cells)*k2];
//    FLOAT *nvmeannn=&nvmean[(nn*n_cells)*k2];
//    FLOAT *wmeannn=&wmean[(nn*n_cells)*k2];
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
//      colornn[u1]=0;//no color for each unit
//      ncolnn[u1]=0;//number of units for a color
//    }
//
//    for(u1=0;u1<n_cells;u1++){
//      if(colornn[u1]!=0) continue;
//      col++;
//      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] = nvnn[(u1)*k2+j];
//      for(j=0;j<k  ;j++) wmeannn[col*k2+j] = netnn->cell[u1].w[j];
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
//	      for(j=0;j<k  ;j++) wmeannn[col*k2+j]+= netnn->cell[u2].w[j];
//	      //	    if(u1<10) fprintf(stderr,"check u1=%d,u2=%d.\n",u1,u2);
//	    }
//	  }
//	}
//      }
//      for(j=0;j<=k1;j++) nvmeannn[col*k2+j] /=(ncolnn[col]);
//      //      fprintf(stderr,"check npnn[%d]=%d\n",u1,npnn[u1]);
//
//
//
//    }
//    int i;
//    int colmax=col;
//    FILE *fp;
//    if(maxnp>colmax) maxnp=colmax;
//    //sorting below
//    for(i=0;i<colmax;i++) sncolnn[i]=i;//prepare
//    for(i=0;i<colmax-1;i++){//sort
//      for(j=i+1;j<colmax;j++){
//	if(ncolnn[sncolnn[j]]>ncolnn[sncolnn[i]]){
//	  int ntmp=sncolnn[j];
//	  sncolnn[j]=sncolnn[i];
//	  sncolnn[i]=ntmp;
//	}
//      }
//    }
//    //sorting above
//    {
//      for(j=0;j<colmax;j++){//initialize
//	sprintf(fname,"tmp/wplane%d.dat",j);
//	fp=fopen(fname,"w");
//	fclose(fp);
//      }
//      for(i=0;i<colmax;i++){//
//	col=sncolnn[i];
//	if(col==0) continue;//no color
//	fprintf(stdout,"Plane#%d: with %d units; col#=%d);",i,ncolnn[col],col);
//	fprintf(stdout,"\nNormal Vector (Hohsen)=");
//	
//	for(j=0;j<k1;j++) fprintf(stdout,"%+.5f ",nvmeannn[col*k2+j]);
//	
//	fprintf(stdout,"Distance from origin=%.5f\n",nvmeannn[col*k2+k1]);
//	
//	sprintf(fname,"tmp/wplane%d.dat",col);
//	fp=fopen(fname,"a+");
//	
//	for(u1=0;u1<n_cells;u1++){
//	  if(colornn[u1]==col){
//	    FLOAT *w=&(netnn->cell[u1].w[0]);
//	    //	    fprintf(stderr,"u%d",u1);//check
//	    for(j=0;j<k;j++){
//	      fprintf(fp,"%e ",w[j]);
//	    }
//	    FLOAT y=0;
//	    for(j=0;j<k;j++) y += (-nvmean[col*k2+j]/nvmean[col*k2+k])*w[j];
//	    y+=(-nvmean[col*k2+k1]/nvmean[col*k2+k]);
//	    //	  {
//	    //	    FLOAT *x=(FLOAT*)malloc(sizeof(FLOAT)*k1),_yt,_y,_Yt;
//	    //	    for(j=0;j<k;j++) x[j]=w[j];
//	    //	    x[k]=1;
//	    //	    calc_output(&net[nn], x,&_yt,&_y,&_Yt);
//	    //	    y=_Yt;
//	    //	  }
//	    
//	    fprintf(fp,"%e\n",y);
//	    if(i==0){//check
//	      fprintf(stderr,"\nw=");
//	      for(j=0;j<k;j++) fprintf(stderr,"%f ",netnn->cell[u1].w[j]);
//	      fprintf(stderr,"y=%f",y);
//	      
//	      fprintf(stderr,"\nm=");
//	      for(j=0;j<k;j++) fprintf(stderr,"%f ",(-nvmean[col*k2+j]/nvmean[col*k2+k]));
//	      fprintf(stderr,"%f ",(-nvmean[col*k2+k1]/nvmean[col*k2+k]));
//	      
//	      
//	      fprintf(stderr,"\nM=");
//	      for(j=0;j<=k;j++) fprintf(stderr,"%f ",netnn->cell[u1].am.M[0][j]);
//	      fprintf(stderr,"v=%f\nn=",netnn->cell[u1].v);
//	      
//	      for(j=0;j<=k1;j++) fprintf(stderr,"%f ",nvnn[(u1)*k2+j]);
//	      fprintf(stderr,"\n");
//	    }
//	  }
//	}
//	fclose(fp);
//      }
//      fp=fopen("tmp/wplanes.plt","w");
//      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
//      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1",
//	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      for(i=0;i<colmax;i++){
//	col=sncolnn[i];
//	//	if(ncolnn[col]>=2){
//	if(i<maxnp){
//	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" w l",col,col,ncolnn[col]);
//	}
//	else {
//	  fprintf(fp,", \"wplane%d.dat\" using 1:2:3 t \"wplane%d(%dunits)\" w p",col,col,ncolnn[col]);
//	}
//      }
//      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//      
//      fprintf(fp,"set termi tgif\n");
//      fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
//      //
//      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1",
//	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      for(i=0;i<colmax;i++){
//	col=sncolnn[i];
//	fprintf(fp,", \"wplane%d.dat\" using 1:2:3 w p",col);
//	if(ncolnn[col]<=1) fprintf(fp," ps 0.1");
//      }
//      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//      
//      fprintf(fp,"set termi tgif\n");
//      fprintf(fp,"set output \"wplanes.obj\"\nreplot\n");
//      fclose(fp);
//      system("cd tmp;xterm -geometry 80x5-0+50 -e gnuplot wplanes.plt&");
//    }
//    {
//      int t;
//      FLOAT _y,_yr,_Y;
//      for(j=0;j<colmax;j++){//initialize
//	sprintf(fname,"tmp/plane%d.dat",j);
//	fp=fopen(fname,"w");
//	fclose(fp);
//      }
//      for (t=givendata->n_train; t<givendata->n_total; t++) {
//	//      for (j=0; j<=k; j++) test->x[t][j] = givendata->x[t][j];
//	//      calc_output(netnn, test->x[t],&test->y[t],&y,&(test->Y[t]));
//	calc_output(netnn, givendata->x[t],&_y,&_yr,&_Y);
//	col=color[netnn->c];
//	sprintf(fname,"tmp/plane%d.dat",col);
//	fp=fopen(fname,"a+");
//	for(j=0; j<k; j++) fprintf(fp,"%e ",test->x[t][j]);
//	fprintf(fp,"%e\n",_Y);
//	fclose(fp);
//      }      
//      fp=fopen("tmp/planes.plt","w");
//      fprintf(fp,"set view 144,181; set pointsize 1.0;\n");
//      fprintf(fp,"splot [%f:%f][%f:%f] \"../%s\" using 1:2:3 pointsize 0.1",
//	    //	    test->xmin0[0],test->xmax0[0],test->xmin0[1],test->xmax0[1],givendata->path);
//	    //	    netnn->xmin0[0],netnn->xmax0[0],netnn->xmin0[1],netnn->xmax0[1],givendata->path);
//	      net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//      for(i=0;i<colmax;i++){
//	col=sncolnn[i];
//	fprintf(fp,", \"plane%d.dat\" using 1:2:3 t \"plane%d(%dunits)\" w p",col,col,ncolnn[col]);
//	//	fprintf(fp,", \"plane%d.dat\" using 1:2:3 w p",col);
//	if(ncolnn[col]<1) fprintf(fp," ps 0.1");
//      }
//      fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//      
//      fprintf(fp,"set termi tgif\n");
//      fprintf(fp,"set output \"planes.obj\"\nreplot\n");
//      fclose(fp);
//      system("cd tmp;xterm -geometry 80x5-0+50 -e gnuplot planes.plt&");
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
void exec_ssp_test_rt(NET *net, DATA *givendata, DATA *test) { //p2p ssprt 
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_test = givendata->n_test;
  FLOAT y;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  FILE *fp;
  int i;
  char line[256];
  int GasosuH,GasosuV;
  FLOAT GakakuH,GakakuV,roll,pitch,yaw;
  FLOAT cH1,cD1,cV1,cH,cD,cV;
  FLOAT *yx,*x1,*x2,dd,dd0;
  FLOAT ddmin,ddminx1,ddminx2,ddminyx,ddmindc;
  int *dc;
  int ViewV=60,ViewH=330,zahyo0,zahyo1;//,zahyohenkan
  FLOAT xxxx[3],yyyy,YYYY;
  FLOAT cr,sr,cp,sp,cy,sy;
  FLOAT JV,JH,JD,XcwH,XcwD,XcwV,t;
  FLOAT *nH,*nD,*nV;
  FLOAT ZjH,ZjD,ZjV;
  FLOAT ZjH0,ZjD0,ZjV0;
  FLOAT XpH,XpD,XpV;
  FLOAT Xpx1,Xpx2,Xpy;
  int jV,jH,ii;
  int atta;
  FLOAT Hmin,Hmax,Vmin,Vmax,dH,dV;
  FLOAT Dmin,Dmax;
  FLOAT *wH,*wD,*wV;
  FLOAT ny,nx1,nx2,*wx1,*wx2,*wy;
  FLOAT cwx1,swx1,scwx1;
  FLOAT cwx2,swx2,scwx2;
  //  FLOAT x1,x2,y;
  printf("cos(phi) for smoothing (2 for y0x), zahyo_orig, zahyo_target(0:Polar,1:Rect):");
  fgets(line,256,stdin);

  sscanf(line,"%lf%d%d",&(net->cost),&zahyo0,&zahyo1);
  //  sscanf(line,"%lf%d",&(net->cost),&zahyohenkan);
  printf("Camera parameters (GasosuH,GasosuV,GakakuH,GakakuV):");
  fgets(line,80,stdin);
  sscanf(line,"%d%d%lf%lf",&GasosuH,&GasosuV,&GakakuH,&GakakuV);

  printf("Camera pose (Horizontal,Depth,Vertical,roll,pitch(tilt),yaw(pan):");
  fgets(line,256,stdin);
  sscanf(line,"%lf%lf%lf%lf%lf%lf",&cH,&cD,&cV,&roll,&pitch,&yaw);

  fp=open_file("rangegot.dat","w");  //write predict.dat
  {//zahyohenkan
//    yx=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x1=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x2=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    dc=(int *)calloc(GasosuH*GasosuV,sizeof(int));

    yx=(FLOAT *)malloc(GasosuH*GasosuV*sizeof(FLOAT));
    x1=(FLOAT *)malloc(GasosuH*GasosuV*sizeof(FLOAT));
    x2=(FLOAT *)malloc(GasosuH*GasosuV*sizeof(FLOAT));
    dc=(int *)malloc(GasosuH*GasosuV*sizeof(int));

    wH=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    wD=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    wV=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    wx1=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    wx2=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    wy =(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    nH =(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    nD =(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    nV =(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
    
    xxxx[2]=1;
    cr=cos(roll/Degree);
    sr=sin(roll/Degree);
    cp=cos(pitch/Degree);
    sp=sin(pitch/Degree);
    cy=cos(yaw/Degree);
    sy=sin(yaw/Degree);

    {
      Hmin=Vmin=1e30;
      Hmax=Vmax=-1e30;
      for(i=0;i<net->n_cells;i++){
	//	if(zahyohenkan==Polar2Polar || zahyohenkan==Polar2Rect){
	if(zahyo0==Polar){
	  wx1[i]=moverange1(net->cell[i].w[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0])/Degree;
	  wx2[i]=moverange1(net->cell[i].w[1], net->xmin[1],net->xmax[1],net->xmin0[1],net->xmax0[1])/Degree;
	  wy[i] =moverange1(net->cell[i].am.M[0][2]
			   +net->cell[i].am.M[0][0]*net->cell[i].w[0]
			   +net->cell[i].am.M[0][1]*net->cell[i].w[1], 
			   net->ymin0,net->ymax,net->ymin0,net->ymax0);
	  wH[i]= wy[i]*cos(wx2[i])*sin(wx1[i]);
	  wD[i]= wy[i]*cos(wx2[i])*cos(wx1[i]);
	  wV[i]= wy[i]*sin(wx2[i]);

	  ny =(net->ymax-net->ymin)/(net->ymax0-net->ymin0);
	  nx1=-((net->xmax[0]-net->xmin[0])/(net->xmax0[0]-net->xmin0[0]))*(net->cell[i].am.M[0][0]);
	  nx2=-((net->xmax[1]-net->xmin[1])/(net->xmax0[1]-net->xmin0[1]))*(net->cell[i].am.M[0][1]);
	  
	  cwx1=cos(wx1[i]);swx1=sin(wx1[i]);scwx1=1./cos(wx1[i]);
	  cwx2=cos(wx2[i]);swx2=sin(wx2[i]);scwx2=1./cos(wx2[i]);
	  nH[i]=(nx1*cwx1*scwx2 + swx1*(ny*y*cwx2 - nx2*swx2))/wy[i];
	  nD[i]=(-(nx1*scwx2*swx1) + cwx1*(ny*wy[i]*cwx2 - nx2*swx2))/wy[i];
	  nV[i]=nx2*cwx2/wy[i]+ny*swx2;
	  //	    nx1 Cos[x1] Sec[x2] + Sin[x1] (ny y Cos[x2] - nx2 Sin[x2])
	  //Out[37]= {----------------------------------------------------------, 
	  //					y
	  // 
	  //     -(nx1 Sec[x2] Sin[x1]) + Cos[x1] (ny y Cos[x2] - nx2 Sin[x2])
	  //>    -------------------------------------------------------------, 
	  //				     y
	  // 
	  //     nx2 Cos[x2]
	  //>    ----------- + ny Sin[xg2]}
	  //	    y
	}
	//	else if(zahyohenkan==Rect2Polar || zahyohenkan==Rect2Rect){
	else if(zahyo0==Rect){
	  wH[i]= moverange1(net->cell[i].w[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]);
	  wV[i]= moverange1(net->cell[i].w[1], net->xmin[1],net->xmax[1],net->xmin0[1],net->xmax0[1]);
	  wD[i]= moverange1(net->cell[i].am.M[0][2]
			   +net->cell[i].am.M[0][0]*net->cell[i].w[0]
			   +net->cell[i].am.M[0][1]*net->cell[i].w[1], 
			   net->ymin0,net->ymax,net->ymin0,net->ymax0);
	  
	  wy[i] =sqrt(wH[i]*wH[i]+wV[i]*wV[i]+wD[i]*wD[i]);
	  wx1[i]=atan2(wH[i],wD[i]);
	  wx2[i]=atan2(wV[i],sqrt(wH[i]*wH[i]+wV[i]*wV[i]));
		       
	  nD[i]=(net->ymax-net->ymin)/(net->ymax0-net->ymin0);
	  nH[i]=-((net->xmax[0]-net->xmin[0])/(net->xmax0[0]-net->xmin0[0]))*(net->cell[i].am.M[0][0]);
	  nV[i]=-((net->xmax[1]-net->xmin[1])/(net->xmax0[1]-net->xmin0[1]))*(net->cell[i].am.M[0][1]);
	}
//	  JH= wH[i]-cH;
//	  JD= wD[i]-cD;
//	  JV= wV[i]-cV;

//	  JH= wH[i]+cH;
//	  JD= wD[i]+cD;
//	  JV= wV[i]+cV;
//
	JH= wH[i];
	JD= wD[i];
	JV= wV[i];
	// YI.PI.RI. = Inverse R.P.Y
//	  XpH=cp*JD*sy + JV*(-(cy*sr) + cr*sp*sy) + JH*(cr*cy + sp*sr*sy);
//	  XpD=cp*cy*JD + JH*(cy*sp*sr - cr*sy) + JV*(cr*cy*sp + sr*sy);
//	  XpV=cp*cr*JV - JD*sp + cp*JH*sr;
	//R.P.Y
	//	    Zjx=cp JV sr + JD (cy sp sr - cr sy) + JH (cr cy + sp sr sy); 
	//	    Zjy=cp cy JD - JV sp + cp JH sy;
	//	    Zjz=cp cr JV + JH (-(cy sr) + cr sp sy) + JD (cr cy sp + sr sy);
	XpH=wH[i];
	XpD=wD[i];
	XpV=wV[i];
	//	if(zahyohenkan==Polar2Polar || zahyohenkan==Rect2Polar){
	if(zahyo1==Polar){
	  Xpy = sqrt(XpH*XpH+XpD*XpD+XpV*XpV);
	  Xpx1=atan2(XpH,XpD)*Degree;
	  Xpx2=atan2(XpV,sqrt(XpH*XpH+XpD*XpD))*Degree;
	}
	//	else if(zahyohenkan==Polar2Rect || zahyohenkan==Rect2Rect){
	else {//if(zahyo1==Rect){
	  Xpy =XpD;
	  Xpx1=XpH;
	  Xpx2=XpV;
	}
	//	Xpx1=moverange1(Xpx1, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
	//	Xpx2=moverange1(Xpx2, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
	if(Xpx1<Hmin) Hmin=Xpx1;
	if(Xpx1>Hmax) Hmax=Xpx1;
	if(Xpx2<Vmin) Vmin=Xpx2;
	if(Xpx2>Vmax) Vmax=Xpx2;
	if(Xpy<Dmin) Dmin=Xpy;
	if(Xpy>Dmax) Dmax=Xpy;
      }
    }
//    Hmax=net->xmax0[0];
//    Hmin=net->xmin0[0];
//    Vmax=net->xmax0[1];
//    Vmin=net->xmin0[1];
    Hmax=Hmax+(Hmax-Hmin)*0.2;
    Hmin=Hmin-(Hmax-Hmin)*0.2;
    Vmax=Vmax+(Vmax-Vmin)*0.2;
    Vmin=Vmin-(Vmax-Hmin)*0.2;
    dH=(Hmax-Hmin)/GasosuH;  dV=(Vmax-Vmin)/GasosuV;
    //    dH=(Hmax-Hmin)/GasosuH;    dV=(Vmax-Vmin)/GasosuV;
    //    dH=tan((GakakuH/2.)/Degree)/GasosuH/2.;  dV=tan((GakakuV/2.)/Degree)/GasosuV/2.;
    
    printf("Hmin,Hmax,Vmin,Vmax=%+6.3f,%+6.3f,%+6.3f,%+6.3f.\n",Hmin,Hmax,Vmin,Vmax);

    //    if(zahyohenkan==Polar2Polar || zahyohenkan==Rect2Polar){
    if(zahyo1==Polar){
      cH1=cH;
      cD1=cD;
      cV1=cV;
    }
    else {//if(zahyo1==Rect){
//	ZjH = cy*sp*sr - cr*sy;
//	ZjD = cp*cy;
//	ZjV = cr*cy*sp + sr*sy;
//
      ZjH = -(cy*sp*sr) - cr*sy; 
      ZjD = cp*cy;
      ZjV = cr*cy*sp - sr*sy;
      //zj={0,1,0}
      //A.zj
      //Out[21]= {-(cy sp sr) - cr sy, cp cy, cr cy sp - sr sy}

    }
    
    for(jV=0;jV<GasosuV;jV++){
      JV=jV*dV+Vmin;   //      ((2.*jV+1.)/GasosuV-1.)*tan((GakakuV/2.)/Degree);
      //Dmax (-(cy sp sr) - cr sy), cp cy Dmax, Dmax (cr cy sp - sr sy)
      //JV=jV*GakakuV/GasosuV;
      for(jH=0;jH<GasosuH;jH++){
	JH=jH*dH+Hmin;   // if(zahyohenkan==Polar2Polar) JH=((2.*jH+1.)/GasosuH-1.)*tan((GakakuH/2.)/Degree);
      //Dmax (-(cy sp sr) - cr sy), cp cy Dmax, Dmax (cr cy sp - sr sy)
	//JH=jH*GakakuH/GasosuH;
	atta=0;
	//	if(zahyohenkan==Polar2Polar || zahyohenkan==Rect2Polar){
	if(zahyo1==Polar){
	  ZjH0=tan(JH/Degree);
	  ZjD0=1.;
	  ZjV0=sqrt(1.+ZjH0*ZjH0)*tan(JV/Degree);
	  
	  ZjH=-cp*ZjV0*sr + ZjD0*(-cy*sp*sr - cr*sy) + ZjH0*(cr*cy - sp*sr*sy); 
	  ZjD=cp*cy*ZjD0 - ZjV0*sp + cp*ZjH0*sy;
	  ZjV=cp*cr*ZjV0 + ZjH0*(cy*sr + cr*sp*sy) + ZjD0*(cr*cy*sp - sr*sy);
	}
	//	else if(zahyohenkan==Polar2Rect || zahyohenkan==Rect2Rect){
	else {//if(zahyo1==Rect){
	  JV-=Dmax*(cr*cy*sp - sr*sy);
	  JH-=Dmax*(-(cy*sp*sr) - cr*sy);
	  cH1=-cp*JV*sr + JH*(cr*cy - sp*sr*sy)    ;
	  cD1=-(JV*sp) + cp*JH*sy                 ;
	  cV1=cp*cr*JV + JH*(cy*sr + cr*sp*sy) ;
	}
	ddmin=1e30;
	dd0=1e30;
	for(i=0;i<net->n_cells;i++){
	  XcwH= cH1-wH[i];
	  XcwD= cD1-wD[i];
	  XcwV= cV1-wV[i];
	  
	  //Solve[{nn.(t Zj + Xcw) == 0}, t] 
	  dd=ZjH*nH[i] + ZjD*nD[i] + ZjV*nV[i]; if(fabs(dd)<1e-20) dd=1e-20;
	  t=(-(nH[i]*XcwH) - nD[i]*XcwD - nV[i]*XcwV)/dd;//Solve[{n.(t Zj + Xcw) == 0}, t] 
	  XpH= t*ZjH + cH1;
	  XpD= t*ZjD + cD1;
	  XpV= t*ZjV + cV1;
	  //	  dd=((d=(Xpx-Xwx))*d)+((d=(Xpy-Xwy))*d)+((d=(Xpy-Xwy))*d);
	  //	  if(zahyohenkan<=1){
	  if(zahyo0==Polar){
	    Xpy = sqrt(XpH*XpH+XpD*XpD+XpV*XpV);
	    Xpx1=atan2(XpH,XpD)*Degree;
	    Xpx2=atan2(XpV,sqrt(XpH*XpH+XpD*XpD))*Degree;
	  }
	  else{
	    Xpy = XpD;//Koten?
	    Xpx1= XpH;//Koten?
	    Xpx2= XpV;//Koten?
	  }
	  if(Xpx1<net->xmin0[0]) Xpx1=net->xmin0[0];
	  if(Xpx1>net->xmax0[0]) Xpx1=net->xmax0[0];
	  if(Xpx2<net->xmin0[1]) Xpx2=net->xmin0[1];
	  if(Xpx2>net->xmax0[1]) Xpx2=net->xmax0[1];
	  //	  if(Xpx1<net->xmin0[0] || Xpx1>net->xmax0[0] || Xpx2<net->xmin0[1] || Xpx2>net->xmax0[1]){
	    //	    fprintf(stdout,"(%f,%f) is out of range\n",Xpx1,Xpx2);
	    //	    continue;//???
	  //	  }
	  xxxx[0]=moverange1(Xpx1, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
	  xxxx[1]=moverange1(Xpx2, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
	  if(fabs(net->cost)>1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
	  else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
//	  if(jV==0 && jH==0){
//	    printf("%+e %+e %+e %+e %+e %+e %+e %+e %+e %f %d %d\n",
//		   Xpx1,Xpx2,Xpy,YYYY,
//		   xxxx[0],xxxx[1],yyyy,
//		   net->cell[net->c].w[0],net->cell[net->c].w[1],
//		   //		   sqrt(((d=(net->cell[i].w[0]-xxxx[0]))*d)+((d=(net->cell[i].w[1]-xxxx[1]))*d)),
//		   sqrt(square(net->cell[i].w[0]-xxxx[0])+square(net->cell[i].w[1]-xxxx[1])),
//		   net->c,i);//check050706
//	  }

	  //	  yx[jV*GasosuH+jH]=YYYY;//for check
	  //	    if(1==0 && jV==0 && jH == 16){
	  //	      //ssprt
	  //	      //-2                   cos(theta) for search neighboring weight vector
	  //	      //50 50 60 60         GasosuH,GasosuV,GakakuH,GakakuV
	  //	      //1 0 0 0 0 0          xc,yc,zc,roll,pitch(tilt),yaw(pan)
	  //	      //0 10 00 10 -30 60 Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin<0 for quit
	  //	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx1,Xpx2,Xpr,YYYY,i,net->c);
	  //	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e i=%d c=%d\n",XpH,XpD,XpV,YYYY,i,net->c);
	  //	      printf("%+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e\n",XpH,XpD,XpV,cH1,cD1,cV1,XcwH,XcwD,XcwV);
	  //	    }
	  //if(jV==1&&jH==0){
	  //printf("???check i,c=%3d,%3d. Y(%+6.0f,%+6.0f)=%+6.0f?%+6.0f,n=%+8.2e,%+8.2e,%+8.2e\n",i,net->c,Xpx1,Xpx2,YYYY,Xpy,nH[i],nD[i],nV[i]);
	  //}
	  if(net->c==i){ 
	    //	      printf("cell %3d w=(%+6.3f,%+6.3f) for (%3d,%3d) xp=(%+6.3f,%+8.3f)<-Xp=(%+6.3f,%+6.3f)\n",
	    //		     i,net->cell[i].w[0],net->cell[i].w[1],
	    //		     jH,jV,xxxx[0],xxxx[1], Xpx1, Xpx2
	    //		     );
	    //	    if(yx[jV*GasosuH+jH]>1){
	    if(atta>0){
	      //	      printf("j[%d,%d]x(%+7.4f,%+7.4f) new r(%d)=%7.2f old r(%d)=%7.2f\n",jV,jH,Xpx1,Xpx2, i,YYYY, dc[jV*GasosuH+jH],yx[jV*GasosuH+jH]);
	      printf("j[%d,%d] newHVD(%g,%g,%g) oldHVD=(%g,%g,%g)\n",jV,jH,Xpx1,Xpx2, YYYY, 
		     x1[jV*GasosuH+jH],
		     x2[jV*GasosuH+jH],
		     yx[jV*GasosuH+jH]);
	    }
	    //	    dd=sqrt((Xpx1-cH1)*(Xpx1-cH1)+(Xpx2-cV1)*(Xpx2-cV1)+(YYYY-cH1)*(YYYY-cH1));
	    //	    if(dd<dd0){
	    if(atta==0 || YYYY<yx[jV*GasosuH+jH]){
	      //	      dd0=dd;
	      //x1[jV*GasosuH+jH]=Xpx1; x2[jV*GasosuH+jH]=Xpx2;
	      //	      x1[jV*GasosuH+jH]=JH;  x2[jV*GasosuH+jH]=JV;
	      yx[jV*GasosuH+jH]=YYYY;
	      x1[jV*GasosuH+jH]=Xpx1;
	      x2[jV*GasosuH+jH]=Xpx2;
	      //	      yx[jV*GasosuH+jH]=cD1; x1[jV*GasosuH+jH]=cH1; x2[jV*GasosuH+jH]=cV1;//for check
	      dc[jV*GasosuH+jH]=net->dc;
	      //	      dd0=dd;
	    }
	    atta++;
	  }
	  else if(atta==0){
	    //	    dd=sqrt(((d=(net->cell[i].w[0]-xxxx[0]))*d)+((d=(net->cell[i].w[1]-xxxx[1]))*d));
	    dd=sqrt(square(net->cell[i].w[0]-xxxx[0])+square(net->cell[i].w[1]-xxxx[1]));
	    //if(dd>0.2*net->xwidth && atta==1) continue;
	    //	    if(dd>0.1*net->xwidth && atta==1) continue;
	    if(dd<ddmin){
	      ddmin=dd;
	  //	  if(Xpx1<net->xmin0[0] || Xpx1>net->xmax0[0] || Xpx2<net->xmin0[1] || Xpx2>net->xmax0[1]){
	      ddminx1=Xpx1;
	      ddminx2=Xpx2;
	      //	      ddminyx=Xpy;//
	      ddminyx=YYYY;

	      ddmindc=net->dc;
	    }
	  }
	}
	if(atta==0){
	  printf("Hen? kana!??? NaiNai j(%d,%d),HVD=%g %g %g>>\n",jV,jH,ddminx1,ddminx2,ddminyx);
	  yx[jV*GasosuH+jH]=ddminyx;
	  if(ddminx1<net->xmin0[0]) ddminx1=net->xmin0[0];
	  if(ddminx1>net->xmax0[0]) ddminx1=net->xmax0[0];
	  if(ddminx2<net->xmin0[1]) ddminx2=net->xmin0[1];
	  if(ddminx2>net->xmax0[1]) ddminx2=net->xmax0[1];
	  //	  printf("HVD=%g %g %g\n",ddminx1,ddminx2,ddminyx);
	  x1[jV*GasosuH+jH]=ddminx1;
	  x2[jV*GasosuH+jH]=ddminx2;
	  dc[jV*GasosuH+jH]=ddmindc;
	  yx[jV*GasosuH+jH]=-1e-30;
	  yx[jV*GasosuH+jH]=ddminyx;
	  //	      yx[jV*GasosuH+jH]=cD1; x1[jV*GasosuH+jH]=cH1; x2[jV*GasosuH+jH]=cV1;//for check
	}
      }
      fprintf(stdout,"."); fflush(stdout);
    }
    for(jV=0;jV<GasosuV;jV++){
      JV=jV*dV+Vmin; 
      for(jH=0;jH<GasosuH;jH++){
	JH=jH*dH+Hmin;
	if(yx[jV*GasosuH+jH]>-1e-30+1){
	  fprintf(fp,"%13.7e %13.7e %13.7e", 
		  //JH,JV,
		  x1[jV*GasosuH+jH],
		  x2[jV*GasosuH+jH],
		  yx[jV*GasosuH+jH]);
	}
	else{
	  fprintf(fp,"%13.7e %13.7e -", 
		  //JH,JV,
		  x1[jV*GasosuH+jH],
		  x2[jV*GasosuH+jH]);
		  //yx[jV*GasosuH+jH]);
	}
	for(ii=-2;ii<=8;ii++){
	  //	  if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
	  if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %13.7e",x1[jV*GasosuH+jH]);
	  else fprintf(fp," -");
	}
	//	fprintf(fp," %13.7e %13.7e", x1[jV*GasosuH+jH],x2[jV*GasosuH+jH]);
	fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
    }
  }
  //  for(ii=0;jH<=14;ii) fprintf(fp,"- ");  fprintf(fp,"%e %e %e\n",cH,cD,cV);
  fclose(fp);
  
  fp=open_file("cameracenter.dat","w");
  fprintf(fp,"%e %e %e\n",cH,cD,cV);
  fclose(fp);
  fp=open_file("weights.dat","w");
  //  if(zahyohenkan<=1){
  if(zahyo0==Polar){
    for(i=0;i<net->n_cells;i++)  fprintf(fp,"%e %e\n",wx1[i]*Degree,wx2[i]*Degree);
  }
  else{
    for(i=0;i<net->n_cells;i++)  fprintf(fp,"%e %e\n",wH[i],wV[i]);
  }
  //  if(zahyohenkan==Polar2Polar)  for(i=0;i<net->n_cells;i++)  fprintf(fp,"%e %e\n",wx1[i],wx2[i]);
  //  else if(zahyohenkan==Polar2Rect)  for(i=0;i<net->n_cells;i++) fprintf(fp,"%e %e\n",wH[i],wV[i]);
  fclose(fp);

  fp=open_file("rangegot.gpl","w");
  
  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 1:2:3\n");
  fprintf(fp,"set term x11\nset view 60,330\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
  fprintf(fp,"set nocontour\nset style data lines\n");
  fprintf(fp,"set label \"X Camera\"  at first %e, first  %e, first  %e\n",cH,cD,cV);
  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p pt 0, \"rangegot.dat\" using 1:3:2 w p pt 1\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"splot \"rangegot.dat\" using 1:3:2 w l 1 0, \"cameracenter.dat\" using 1:3:2 w p pt 2\n");
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"set contour\nset view 60,350\n");
  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3, \"cameracenter.dat\" using 1:2:3 w p pt 2\n\n");
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
  //  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3 w l\npause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
  fprintf(fp,"plot \"rangegot.dat\" u 5:2 t \"no opposite weight (open boundary)\" w p pt 2");
  fprintf(fp,", \"rangegot.dat\" u  6:2 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p pt 2");
  fprintf(fp,", \"rangegot.dat\" u  7:2 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p pt 3");
  fprintf(fp,", \"rangegot.dat\" u  8:2 t \"xi(fu renzoku?) with y0x for t<0\"          w p pt 4");
  fprintf(fp,", \"rangegot.dat\" u  9:2 t \"xi(   renzoku?) with y0x for t>1\"          w p pt 6");
  fprintf(fp,", \"rangegot.dat\" u 10:2 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p pt 0");
  fprintf(fp,", \"rangegot.dat\" u 11:2 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p pt 7");
  fprintf(fp,", \"rangegot.dat\" u 12:2 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p pt 8");
  fprintf(fp,", \"weights.dat\" u 1:2 t \"weight vecotr\"        w p pt 9");
  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l ls 1");
  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
  fclose(fp);
  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");

  //////////////////////////////////////
  {
    //    int Hmin,Hmax,Vmin,Vmax;
    char cmd='r';
    FILE *fp1;
    printf("u,d,r,l,q for up,down,right,left,quit:");
    //    fp = open_pipe(GNUPLOT, "w");
    fp = open_pipe("gnuplot -geometry 510x360-0+0", "w");
    fp1 = open_pipe("gnuplot -geometry 640x420-0+460", "w");
    fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
    fprintf(fp,"set nocontour\nset style data lines\n");
    fprintf(fp,"set label \". Camera\"  at first %e, first  %e, first  %e\n",cH,cD,cV);

    fprintf(fp1,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
    fprintf(fp1,"set nocontour\nset style data lines\n");
    fprintf(fp1,"set label \". Camera\"  at first %e, first  %e, first  %e\n",cH,cD,cV);
    for(i=0;;){
      if(i>0) fgets(line,256,stdin);
      i=1;
      if(line[0]=='q') break;
      if(line[0]=='u' || line[0]=='d' || line[0]=='r' || line[0]=='l') cmd=line[0];
      if(cmd=='u') ViewV+=10;
      else if(cmd=='d') ViewV-=10;
      else if(cmd=='r') ViewH+=10;
      else if(cmd=='l') ViewH-=10;

      ViewH = (ViewH+360)%360;
      ViewV =(ViewV+180)%180;
      fprintf(fp,"set view %d,%d\nsplot \"rangegot.dat\" using 1:3:2 w l ls 1, \"cameracenter.dat\" using 1:3:2 w p pt 2\n",ViewV,ViewH);
      fprintf(fp1,"set view %d,%d\nsplot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p pt 0, \"cameracenter.dat\" using 1:3:2 w p pt 2\n",ViewV,ViewH,net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
      fflush(fp);
      fflush(fp1);
    }
    close_pipe(fp);
    close_pipe(fp1);
  }
  //////////////////////////////////////
  
  MSE2=MSE1=MSE /= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  
  //  printf("=>MSEtest=%f\n",MSE);
  
  test->MSE2 =test->MSE1 =test->MSE  = MSE;
  test->NMSE = NMSE;
  free(yx);
  free(x1);
  free(x2);
  free(dc);
  free(wH);
  free(wD);
  free(wV);
  free(wy);
  free(wx1);
  free(wx2);
  free(nH);
  free(nD);
  free(nV);
  return;
}
//void exec_ssp_test_rt1(NET *net, DATA *givendata, DATA *test) { //p2p
//  //char name[32] = "exec_ssp_test";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_test = givendata->n_test;
//  FLOAT y;
//  FLOAT MSE1=0;
//  FLOAT MSE2=0;
//  FILE *fp;
//  int i;
//  char line[256];
//  int GasosuH,GasosuV;
//  FLOAT GakakuH,GakakuV,roll,pitch,yaw;
//  FLOAT cH1,cD1,cV1,cH,cD,cV;
//  FLOAT *yx,*x1,*x2,dd;
//  int *dc;
//  int ViewV=60,ViewH=330,zahyohenkan;
//  FLOAT xxxx[3],yyyy,YYYY;
//  FLOAT cr,sr,cp,sp,cy,sy;
//  FLOAT JV,JH,JD,yw,Xwx,Xwy,Xwz,XcwH,XcwD,XcwV,t;
//  FLOAT nH,nD,nV;
//  FLOAT ZjH,ZjD,ZjV;
//  FLOAT XpH,XpD,XpV;
//  FLOAT Xpx1,Xpx2,Xpy;
//  int jV,jH,ii;
//  int atta;
//  FLOAT Hmin,Hmax,Vmin,Vmax,dH,dV;
//  FLOAT *wH,*wD,*wV;
//  FLOAT ny,nx1,nx2,wx1,wx2,wy;
//
//  printf("cos(theta)for search neighboring weight vector(2 for y0x), zahyohenkan(0:p2p,1:p2r) :");
//  fgets(line,256,stdin);
//
//  sscanf(line,"%lf%d",&(net->cost),&zahyohenkan);
//  printf("Camera parameters (GasosuH,GasosuV,GakakuH,GakakuV):");
//  fgets(line,80,stdin);
//  sscanf(line,"%d%d%lf%lf",&GasosuH,&GasosuV,&GakakuH,&GakakuV);
//
//  printf("Camera pose (cH,cD,cV,roll,pitch(tilt),yaw(pan):");
//  fgets(line,256,stdin);
//  sscanf(line,"%lf%lf%lf%lf%lf%lf",&cH,&cD,&cV,&roll,&pitch,&yaw);
//
//  fp=open_file("rangegot.dat","w");  //write predict.dat
//  {//zahyohenkan
//    //    r=(FLOAT *)malloc(sizeof(FLOAT)*GasosuH*GasosuV);
//    yx=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x1=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x2=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    dc=(int *)calloc(GasosuH*GasosuV,sizeof(int));
//
//    wH=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
//    wD=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
//    wV=(FLOAT *)malloc(net->n_cells*sizeof(FLOAT));
//    
//    xxxx[2]=1;
//    cr=cos(roll/Degree);
//    sr=sin(roll/Degree);
//    cp=cos(pitch/Degree);
//    sp=sin(pitch/Degree);
//    cy=cos(yaw/Degree);
//    sy=sin(yaw/Degree);
//
//    {
//	Hmin=Vmin=1e30;
//	Hmax=Vmax=-1e30;
//	for(i=0;i<net->n_cells;i++){
//	  yw =net->cell[i].am.M[0][2]
//	    +net->cell[i].am.M[0][0]*net->cell[i].w[0]
//	    +net->cell[i].am.M[0][1]*net->cell[i].w[1];
//	  wx1=moverange1(net->cell[i].w[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]);
//	  wx2=moverange1(net->cell[i].w[1], net->xmin[1],net->xmax[1],net->xmin0[1],net->xmax0[1]);
//	  wy =moverange1(yw, net->ymin0,net->ymax,net->ymin0,net->ymax0);
//	  wH[i]= wy*cos(wx2/Degree)*sin(wx1/Degree);
//	  wD[i]= wy*cos(wx2/Degree)*cos(wx1/Degree);
//	  wV[i]= wy*sin(wx2/Degree);
////
////	  wH[i]=(yw*cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree));
////	  wD[i]=(yw*cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree));
////	  wV[i]=(yw*sin(net->cell[i].am.M[0][1]/Degree))                                     ;
//
//	  JH= wH[i]-cH;
//	  JD= wD[i]-cD;
//	  JV= wV[i]-cV;
//	  // YI.PI.RI. = Inverse R.P.Y
//	  XpH=cp*JD*sy + JV*(-(cy*sr) + cr*sp*sy) + JH*(cr*cy + sp*sr*sy);
//	  XpD=cp*cy*JD + JH*(cy*sp*sr - cr*sy) + JV*(cr*cy*sp + sr*sy);
//	  XpV=cp*cr*JV - JD*sp + cp*JH*sr;
//	  //R.P.Y
//	  //	    Zjx=cp JV sr + JD (cy sp sr - cr sy) + JH (cr cy + sp sr sy); 
//	  //	    Zjy=cp cy JD - JV sp + cp JH sy;
//	  //	    Zjz=cp cr JV + JH (-(cy sr) + cr sp sy) + JD (cr cy sp + sr sy);
//	  if(zahyohenkan==Polar2Polar){
//	    Xpy = sqrt(XpH*XpH+XpD*XpD+XpV*XpV);
//	    Xpx1=atan2(XpH,XpD)*Degree;
//	    Xpx2=atan2(XpV,sqrt(XpH*XpH+XpD*XpD))*Degree;
//	  }
//	  else if(zahyohenkan==Polar2Rect){
//	    Xpy =XpD;
//	    Xpx1=XpH;
//	    Xpx2=XpV;
//	  }
//	  if(Xpx1<Hmin) Hmin=Xpx1;
//	  if(Xpx1>Hmax) Hmax=Xpx1;
//	  if(Xpx2<Vmin) Vmin=Xpx2;
//	  if(Xpx2>Vmax) Vmax=Xpx2;
//	}
//    }
//    dH=(Hmax-Hmin)/GasosuH;  dV=(Vmax-Vmin)/GasosuV;
//    //    dH=(Hmax-Hmin)/GasosuH;    dV=(Vmax-Vmin)/GasosuV;
//    //    dH=tan((GakakuH/2.)/Degree)/GasosuH/2.;  dV=tan((GakakuV/2.)/Degree)/GasosuV/2.;
//
//    printf("Hmin,Hmax,Vmin,Vmax=%+6.3f,%+6.3f,%+6.3f,%+6.3f.\n",Hmin,Hmax,Vmin,Vmax);
//    if(zahyohenkan==Polar2Polar){
//	cH1=cH;
//	cD1=cD;
//	cV1=cV;
//    }
//    else if(zahyohenkan==Polar2Rect){
//	ZjH = cy*sp*sr - cr*sy;
//	ZjD = cp*cy;
//	ZjV = cr*cy*sp + sr*sy;
//    }
//
//    for(jV=0;jV<GasosuV;jV++){
//	JV=jV*dV+Vmin;   //      ((2.*jV+1.)/GasosuV-1.)*tan((GakakuV/2.)/Degree);
//	//JV=jV*GakakuV/GasosuV;
//	for(jH=0;jH<GasosuH;jH++){
//	  JH=jH*dH+Hmin;   // if(zahyohenkan==Polar2Polar) JH=((2.*jH+1.)/GasosuH-1.)*tan((GakakuH/2.)/Degree);
//	  //JH=jH*GakakuH/GasosuH;
//	  atta=0;
//	  if(zahyohenkan==Polar2Polar){
//	    ZjH = cp*JV*sr + cy*sp*sr - cr*sy + JH*(cr*cy + sp*sr*sy);
//	    ZjD = cp*cy - JV*sp + cp*JH*sy; 
//	    ZjV = cp*cr*JV + cr*cy*sp + sr*sy + JH*(-(cy*sr) + cr*sp*sy);
//	  }
//	  //	  Zjx=cp JV sr + JD (cy sp sr - cr sy) + JH (cr cy + sp sr sy); 
//	  //	  Zjy=cp cy JD - JV sp + cp JH sy;
//	  //	  Zjz=cp cr JV + JH (-(cy sr) + cr sp sy) + JD (cr cy sp + sr sy);
//	  
//	  for(i=0;i<net->n_cells;i++){
////	    nH=-cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree) ;//Normal( Hosen) Vector
////	    nD= cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree) ;
////	    nV=-sin(net->cell[i].am.M[0][1]/Degree);
//	    ny =(net->ymax-net->ymin)/(net->ymax0-net->ymin0);
//	    nx1=-((net->xmax[0]-net->xmin[0])/(net->xmax0[0]-net->xmin0[0]))*(net->cell[i].am.M[0][0]/Degree);
//	    nx2=-((net->xmax[1]-net->xmin[1])/(net->xmax0[1]-net->xmin0[1]))*(net->cell[i].am.M[0][1]/Degree);
//	    nH = ny*cos(nx2/Degree)*cos(nx1/Degree);
//	    nD = ny*cos(nx2/Degree)*sin(nx1/Degree);
//	    nV = ny*sin(nx2/Degree);
//
//	    Xwx= wH[i];
//	    Xwy= wD[i];
//	    Xwz= wV[i];
//	    
//	    if(zahyohenkan==Polar2Rect){
//	      cH1=cp*JV*sr + JH*(cr*cy + sp*sr*sy);
//	      cD1=-(JV*sp) + cp*JH*sy;
//	      cV1=cp*cr*JV + JH*(-(cy*sr) + cr*sp*sy);
//	    }
//	    XcwH= cH1-wH[i];
//	    XcwD= cD1-wD[i];
//	    XcwV= cV1-wV[i];
//	    
//	    //Solve[{nn.(t Zj + Xcw) == 0}, t] 
//	    {
//	      dd=ZjH*nH + ZjD*nD + ZjV*nV; if(fabs(dd)<1e-20) dd=1e-20;
//	      t=(-(nH*XcwH) - nD*XcwD - nV*XcwV)/dd;//Solve[{nn.(t Zj + Xcw) == 0}, t] 
//	      XpH= t*ZjH + cH1;
//	      XpD= t*ZjD + cD1;
//	      XpV= t*ZjV + cV1;
//	    }
//	    //	  dd=((d=(Xpx-Xwx))*d)+((d=(Xpy-Xwy))*d)+((d=(Xpy-Xwy))*d);
//	    Xpy = sqrt(XpH*XpH+XpD*XpD+XpV*XpV);
//	    Xpx1=atan2(XpH,XpD)*Degree;
//	    Xpx2=atan2(XpV,sqrt(XpH*XpH+XpD*XpD))*Degree;
//
//	    xxxx[0]=moverange1(Xpx1, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
//	    xxxx[1]=moverange1(Xpx2, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
//	    //{double d;
//	    //	  dd=sqrt(((d=(net->cell[i].am.M[0][0]-xxxx[0]))*d)+((d=(net->cell[i].am.M[0][1]-xxxx[1]))*d));
//	    //	  if(dd>0.2*net->xwidth) continue;
//	    //}
//	    if(fabs(net->cost)>1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
//	    else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
//	    if(jV==0 && jH == 16){
//	      //ssprt
//	      //-2                   cos(theta) for search neighboring weight vector
//	      //50 50 60 60         GasosuH,GasosuV,GakakuH,GakakuV
//	      //1 0 0 0 0 0          xc,yc,zc,roll,pitch(tilt),yaw(pan)
//	      //0 10 00 10 -30 60 Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin<0 for quit
//	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx1,Xpx2,Xpr,YYYY,i,net->c);
//	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e i=%d c=%d\n",XpH,XpD,XpV,YYYY,i,net->c);
//	      printf("%+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e %+13.7e\n",XpH,XpD,XpV,cH1,cD1,cV1,XcwH,XcwD,XcwV);
//	    }
//	    if(net->c==i){ 
//	      //	      printf("cell %3d w=(%+6.3f,%+6.3f) for (%3d,%3d) xp=(%+6.3f,%+8.3f)<-Xp=(%+6.3f,%+6.3f)\n",
//	      //		     i,net->cell[i].w[0],net->cell[i].w[1],
//	      //		     jH,jV,xxxx[0],xxxx[1], Xpx1, Xpx2
//	      //		     );
//	      atta=1;
//	      if(yx[jV*GasosuH+jH]>1){
//		//	      printf("(%d,%d) new r(%d)=%7.2f old r(%d)=%7.2f\n",jV,jH, i,YYYY, dc[jV*GasosuH+jH],yx[jV*GasosuH+jH]);
//	      }
//	      if(yx[jV*GasosuH+jH]<=1 || YYYY<yx[jV*GasosuH+jH]){
//		yx[jV*GasosuH+jH]=YYYY;
//		x1[jV*GasosuH+jH]=xxxx[0];
//		x2[jV*GasosuH+jH]=xxxx[1];
//		dc[jV*GasosuH+jH]=net->dc;
//	      }
//	      //	    break; 
//	    }
//	  }
//	  if(atta==0) {
//	    printf("Hendayo!??? nakatta (%d,%d)\n",jV,jH);
//	  }
//	}
//	fprintf(stdout,"."); fflush(stdout);
//    }
//    for(jV=0;jV<GasosuV;jV++){JV=jV*dV+Vmin; 
//	for(jH=0;jH<GasosuH;jH++){JH=jH*dH+Hmin;
//	  fprintf(fp,"%13.7e %13.7e %13.7e", 
//		  JH,
//		  JV,
//		  yx[jV*GasosuH+jH]);
//	  for(ii=-2;ii<=8;ii++){
//	    //	  if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
//	    if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %d",jH);
//	    else fprintf(fp," -");
//	  }
//	  fprintf(fp," %13.7e %13.7e", x1[jV*GasosuH+jH],x2[jV*GasosuH+jH]);
//	  fprintf(fp,"\n");
//	}
//	fprintf(fp,"\n");
//    }
//  }
//  fclose(fp);
//  
//  fp=open_file("rangegot.gpl","w");
//  
//  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
//  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"set term x11\nset view 60,330\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"set nocontour\nset data style lines\n");
//  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p 0 0, \"rangegot.dat\" using 1:3:2 w p 1 0\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"splot \"rangegot.dat\" using 1:3:2 w l 1 0\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"set contour\nset view 60,350\n");
//  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//  //  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3 w l\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:2 t \"no opposite weight (open boundary)\" w p 2 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:2 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:2 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:2 t \"xi(fu renzoku?) with y0x for t<0\"          w p 1 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:2 t \"xi(   renzoku?) with y0x for t>1\"          w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 10:2 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 11:2 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p 3 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:2 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p 3 8");
//  fprintf(fp,", \"weights.dat\" u 1:2 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l 3 1");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//
//  //////////////////////////////////////
//  {
//    int Hmin,Hmax,Vmin,Vmax;
//    for(;;){
//	printf("Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin=-1 for quit:");
//	fgets(line,256,stdin);
//	sscanf(line,"%d%d%d%d%d%d",&Hmin,&Hmax,&Vmin,&Vmax,&ViewH,&ViewV);
//	if(Hmin==Hmax) break;
//	if(Hmin<0) break;
//	if(Vmin<0) Vmin=0;
//	if(Hmax>GasosuH) Hmax=GasosuH;
//	if(Vmax>GasosuV) Vmax=GasosuV;
//	ViewH = (ViewH+360)%360;
//	ViewV =(ViewV+180)%180;
//
//	fp = open_pipe(GNUPLOT, "w");
//	//      fpd  = open_file("rangedisp.dat", "w");
//	
//	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
//	fprintf(fp,"set nocontour\nset data style lines\n");
//	fprintf(fp,"splot [%d:%d][][%d:%d] \"rangegot.dat\" using 1:3:2 w l 1 0\n",Hmin,Hmax,Vmin,Vmax);
//	//      fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//	//	for(jV=0;jV<GasosuV;jV++){
//	//	  for(jH=0;jH<GasosuH;jH++){
//	//	    fprintf(fp,"%d %d %13.7e", 
//	//		    jH,
//	//		    jV,
//	//		    r[jV*GasosuH+jH]);
//	//	    for(ii=-2;ii<=8;ii++){
//	//	      if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %e",jH);
//	//	      else fprintf(fp," -");
//	//	    }
//	//	    fprintf(fp,"\n");
//	//	  }
//	//	  fprintf(fp,"\n");
//	//	}
//	//	fprintf(fp, "e\n");
//	fflush(fp);
//	//
//	//	fp=open_file("rangedisp.gpl","w");
//	//	
//	//	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//	//	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
//	//	fprintf(fp,"set nocontour\nset data style lines\n");
//	//	fprintf(fp,"splot [%e:%e][][%e:%e] \"rangegot.dat\" using 1:3:2 w p 1 0\n",Hmin,Hmax,Vmin,Vmax);
//	//	fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//	//	fclose(fp);
//	//
//	//	system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//    }
//    fclose(fp);
//  }
//  //////////////////////////////////////
//  
//  MSE2=MSE1=MSE /= n_test;
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//  
//  printf("=>MSEtest=%f\n",MSE);
//  
//  test->MSE2 =test->MSE1 =test->MSE  = MSE;
//  test->NMSE = NMSE;
//  free(yx);
//  free(x1);
//  free(x2);
//  free(dc);
//  free(wH);
//  free(wD);
//  free(wV);
//  return;
//}
//void exec_ssp_test_rt0(NET *net, DATA *givendata, DATA *test) { //p2p
//  //char name[32] = "exec_ssp_test";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_test = givendata->n_test;
//  FLOAT y;
//  FLOAT MSE1=0;
//  FLOAT MSE2=0;
//  FILE *fp;
//  int i;
//  char line[256];
//  int GasosuH,GasosuV;
//  FLOAT GakakuH,GakakuV,roll,pitch,yaw;
//  FLOAT xc1,yc1,zc1,xc,yc,zc;
//  FLOAT *r,*x1,*x2,dd;
//  int *dc;
//  int ViewV=60,ViewH=330,zahyohenkan;
//
//  printf("cos(theta)for search neighboring weight vector(2 for y0x), zahyohenkan(0:p2p,1:p2r) :");
//  fgets(line,256,stdin);
//  sscanf(line,"%lf%d",&(net->cost),&zahyohenkan);
//  printf("Camera parameters (GasosuH,GasosuV,GakakuH,GakakuV):");
//  fgets(line,80,stdin);
//  sscanf(line,"%d%d%lf%lf",&GasosuH,&GasosuV,&GakakuH,&GakakuV);
//
//  printf("Camera pose (xc,yc,zc,roll,pitch(tilt),yaw(pan):");
//  fgets(line,256,stdin);
//  sscanf(line,"%lf%lf%lf%lf%lf%lf",&xc,&yc,&zc,&roll,&pitch,&yaw);
//
//  fp=open_file("rangegot.dat","w");  //write predict.dat
//  {//zahyohenkan
//    FLOAT xxxx[3],yyyy,YYYY;
//    FLOAT cr,sr,cp,sp,cy,sy;
//    FLOAT JV,JH,yw,Nx,Ny,Nz,Xwx,Xwy,Xwz,Xcwx,Xcwy,Xcwz,t;
//    FLOAT Zjx,Zjy,Zjz;
//    FLOAT Xpx,Xpy,Xpz,Xpr;
//    FLOAT Xpx1,Xpx2;
//    int jV,jH,ii;
//    int atta;
//
//    //    r=(FLOAT *)malloc(sizeof(FLOAT)*GasosuH*GasosuV);
//    r=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x1=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    x2=(FLOAT *)calloc(GasosuH*GasosuV,sizeof(FLOAT));
//    dc=(int *)calloc(GasosuH*GasosuV,sizeof(int));
//    
//    xxxx[2]=1;
//    cr=cos(roll/Degree);
//    sr=sin(roll/Degree);
//    cp=cos(pitch/Degree);
//    sp=sin(pitch/Degree);
//    cy=cos(yaw/Degree);
//    sy=sin(yaw/Degree);
//    if(zahyohenkan==Polar2Polar){
//	xc1=xc;
//	yc1=yc;
//	zc1=zc;
//    }
//    else if(zahyohenkan==Polar2Rect){
//	Zjx = cy*sp*sr - cr*sy;
//	Zjy = cp*cy;
//	Zjz = cr*cy*sp + sr*sy;
//    }
//    for(jV=0;jV<GasosuV;jV++){
//	if(zahyohenkan==Polar2Polar) JV=((2.*jV+1.)/GasosuV-1.)*tan((GakakuV/2.)/Degree);
//	//JV=jV*GakakuV/GasosuV;
//	for(jH=0;jH<GasosuH;jH++){
//	  if(zahyohenkan==Polar2Polar) JH=((2.*jH+1.)/GasosuH-1.)*tan((GakakuH/2.)/Degree);
//	  //JH=jH*GakakuH/GasosuH;
//	  atta=0;
//	  if(zahyohenkan==Polar2Polar){
//	    Zjx = cp*JV*sr + cy*sp*sr - cr*sy + JH*(cr*cy + sp*sr*sy);
//	    Zjy = cp*cy - JV*sp + cp*JH*sy; 
//	    Zjz = cp*cr*JV + cr*cy*sp + sr*sy + JH*(-(cy*sr) + cr*sp*sy);
//	  }
////	  Zjx=cp JV sr + JD (cy sp sr - cr sy) + JH (cr cy + sp sr sy); 
////	  Zjy=cp cy JD - JV sp + cp JH sy;
////	  Zjz=cp cr JV + JH (-(cy sr) + cr sp sy) + JD (cr cy sp + sr sy);
//	  
//	  for(i=0;i<net->n_cells;i++){
//	    Nx=-cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree) ;//Normal( Hosen) Vector
//	    Ny= cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree) ;
//	    Nz=-sin(net->cell[i].am.M[0][1]/Degree);
//	    
//	    yw =net->cell[i].am.M[0][2]
//	      +net->cell[i].am.M[0][0]*net->cell[i].w[0]
//	      +net->cell[i].am.M[0][1]*net->cell[i].w[1];
//	    Xwx= (yw*cos(net->cell[i].am.M[0][1]/Degree)*sin(net->cell[i].am.M[0][0]/Degree)) ;
//	    Xwy= (yw*cos(net->cell[i].am.M[0][1]/Degree)*cos(net->cell[i].am.M[0][0]/Degree)) ;
//	    Xwz= (yw*sin(net->cell[i].am.M[0][1]/Degree));
//	    if(zahyohenkan==Polar2Rect){
//	      JV = jV*GakakuV/GasosuV;
//	      JH = jH*GakakuH/GasosuH;
//	      xc1=cp*JV*sr + JH*(cr*cy + sp*sr*sy);
//	      yc1=-(JV*sp) + cp*JH*sy;
//	      zc1=Xcwz=cp*cr*JV + JH*(-(cy*sr) + cr*sp*sy);
//	    }
//	    Xcwx= xc1-Xwx;
//	    Xcwy= yc1-Xwy;
//	    Xcwz= zc1-Xwz;
//
//	    //Solve[{nn.(t Zj + Xcw) == 0}, t] 
//	    {
//	      dd=Zjx*Nx + Zjy*Ny + Zjz*Nz; if(fabs(dd)<1e-20) dd=1e-20;
//	      t=(-(Nx*Xcwx) - Ny*Xcwy - Nz*Xcwz)/dd;//Solve[{nn.(t Zj + Xcw) == 0}, t] 
//	      Xpx= t*Zjx + xc1;
//	      Xpy= t*Zjy + yc1;
//	      Xpz= t*Zjz + zc1;
//	    }
//	    //	  dd=((d=(Xpx-Xwx))*d)+((d=(Xpy-Xwy))*d)+((d=(Xpy-Xwy))*d);
//	    Xpr= sqrt(Xpx*Xpx+Xpy*Xpy+Xpz*Xpz);
//	    Xpx1=atan2(Xpx,Xpy)*Degree;
//	    Xpx2=atan2(Xpz,sqrt(Xpx*Xpx+Xpy*Xpy))*Degree;
//	    xxxx[0]=moverange1(Xpx1, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
//	    xxxx[1]=moverange1(Xpx2, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
//	    //{double d;
//	    //	  dd=sqrt(((d=(net->cell[i].am.M[0][0]-xxxx[0]))*d)+((d=(net->cell[i].am.M[0][1]-xxxx[1]))*d));
//	    //	  if(dd>0.2*net->xwidth) continue;
//	    //}
//	    if(fabs(net->cost)>1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
//	    else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
//	    if(jV==0 && jH == 16){
//	    //ssprt
//	    //-2                   cos(theta) for search neighboring weight vector
//	    //50 50 60 60         GasosuH,GasosuV,GakakuH,GakakuV
//	    //1 0 0 0 0 0          xc,yc,zc,roll,pitch(tilt),yaw(pan)
//	    //0 10 00 10 -30 60 Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin<0 for quit
//	      //	    printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx1,Xpx2,Xpr,YYYY,i,net->c);
//	      printf("%+13.7e %+13.7e %+13.7e %+13.7e hendayo i=%d c=%d\n",Xpx,Xpy,Xpr,YYYY,i,net->c);
//	    }
//	    if(net->c==i){ 
////	      printf("cell %3d w=(%+6.3f,%+6.3f) for (%3d,%3d) xp=(%+6.3f,%+8.3f)<-Xp=(%+6.3f,%+6.3f)\n",
////		     i,net->cell[i].w[0],net->cell[i].w[1],
////		     jH,jV,xxxx[0],xxxx[1], Xpx1, Xpx2
////		     );
//	      atta=1;
//	      if(r[jV*GasosuH+jH]>1){
//		printf("(%d,%d) new r(%d)=%7.2f old r(%d)=%7.2f\n",jV,jH, i,YYYY, dc[jV*GasosuH+jH],r[jV*GasosuH+jH]);
//	      }
//	      if(r[jV*GasosuH+jH]<=1 || YYYY<r[jV*GasosuH+jH]){
//		r[jV*GasosuH+jH]=YYYY;
//		x1[jV*GasosuH+jH]=xxxx[0];
//		x2[jV*GasosuH+jH]=xxxx[1];
//		dc[jV*GasosuH+jH]=net->dc;
//	      }
//	      //	    break; 
//	    }
//	  }
//	  if(atta==0) {
//	    printf("Hendayo!??? nakatta (%d,%d)\n",jV,jH);
//	  }
//	}
//	fprintf(stdout,"."); fflush(stdout);
//    }
//    for(jV=0;jV<GasosuV;jV++){
//	for(jH=0;jH<GasosuH;jH++){
//	  fprintf(fp,"%d %d %13.7e", 
//		  jH,
//		  jV,
//		  r[jV*GasosuH+jH]);
//	  for(ii=-2;ii<=8;ii++){
//	    //	  if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
//	    if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %d",jH);
//	    else fprintf(fp," -");
//	  }
//	  fprintf(fp," %13.7e %13.7e", x1[jV*GasosuH+jH],x2[jV*GasosuH+jH]);
//	  fprintf(fp,"\n");
//	}
//	fprintf(fp,"\n");
//    }
//  }
//  fclose(fp);
//  
//  fp=open_file("rangegot.gpl","w");
//  
//  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
//  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"set term x11\nset view 60,330\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"set nocontour\nset data style lines\n");
//  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p 0 0, \"rangegot.dat\" using 1:3:2 w p 1 0\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"splot \"rangegot.dat\" using 1:3:2 w l 1 0\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//
//  fprintf(fp,"set contour\nset view 60,350\n");
//  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3\n");
//  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//  //  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3 w l\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:2 t \"no opposite weight (open boundary)\" w p 2 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:2 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:2 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:2 t \"xi(fu renzoku?) with y0x for t<0\"          w p 1 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:2 t \"xi(   renzoku?) with y0x for t>1\"          w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 10:2 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 11:2 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p 3 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:2 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p 3 8");
//  fprintf(fp,", \"weights.dat\" u 1:2 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l 3 1");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//
//  //////////////////////////////////////
//  {
//    int Hmin,Hmax,Vmin,Vmax;
//    for(;;){
//	printf("Range of Display (Hmin,Hmax,Vmin,Vmax) and View (Yoko in [0:360], Tate in [0:180]), Hmin=-1 for quit:");
//	fgets(line,256,stdin);
//	sscanf(line,"%d%d%d%d%d%d",&Hmin,&Hmax,&Vmin,&Vmax,&ViewH,&ViewV);
//	if(Hmin==Hmax) break;
//	if(Hmin<0) break;
//	if(Vmin<0) Vmin=0;
//	if(Hmax>GasosuH) Hmax=GasosuH;
//	if(Vmax>GasosuV) Vmax=GasosuV;
//	ViewH = (ViewH+360)%360;
//	ViewV =(ViewV+180)%180;
//
//	fp = open_pipe(GNUPLOT, "w");
//	//      fpd  = open_file("rangedisp.dat", "w");
//	
//	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
//	fprintf(fp,"set nocontour\nset data style lines\n");
//	fprintf(fp,"splot [%d:%d][][%d:%d] \"rangegot.dat\" using 1:3:2 w l 1 0\n",Hmin,Hmax,Vmin,Vmax);
//	//      fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
//	//	for(jV=0;jV<GasosuV;jV++){
//	//	  for(jH=0;jH<GasosuH;jH++){
//	//	    fprintf(fp,"%d %d %13.7e", 
//	//		    jH,
//	//		    jV,
//	//		    r[jV*GasosuH+jH]);
//	//	    for(ii=-2;ii<=8;ii++){
//	//	      if(dc[jV*GasosuH+jH]==ii) fprintf(fp," %e",jH);
//	//	      else fprintf(fp," -");
//	//	    }
//	//	    fprintf(fp,"\n");
//	//	  }
//	//	  fprintf(fp,"\n");
//	//	}
//	//	fprintf(fp, "e\n");
//	fflush(fp);
////
////	fp=open_file("rangedisp.gpl","w");
////	
////	fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
////	fprintf(fp,"set view %d,%d\n",ViewV,ViewH);
////	fprintf(fp,"set nocontour\nset data style lines\n");
////	fprintf(fp,"splot [%e:%e][][%e:%e] \"rangegot.dat\" using 1:3:2 w p 1 0\n",Hmin,Hmax,Vmin,Vmax);
////	fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
////	fclose(fp);
////
////	system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//    }
//    fclose(fp);
//  }
//  //////////////////////////////////////
//
//  MSE2=MSE1=MSE /= n_test;
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//
//  printf("=>MSEtest=%f\n",MSE);
//
//  test->MSE2 =test->MSE1 =test->MSE  = MSE;
//  test->NMSE = NMSE;
//  free(r);
//  free(x1);
//  free(x2);
//  free(dc);
//  return;
//}

void exec_ssp_test_r(NET *net, DATA *givendata, DATA *test) {//sspr
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  //  FLOAT p1,p2,p1_sum=0,p2_sum=0;
  FLOAT e2;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  //  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
  FILE *fp;
  int i;
  char line[256];
  int p2r;//transform the polar coordinate system to rectangle one (1), inverse(-1), non(0)
  FLOAT x0;

  printf("cos(theta) for search neighboring weight vector(-2 for y0x), p2r: ");
  fgets(line,256,stdin);
  sscanf(line,"%lf%d",&(net->cost),&p2r);

  if((fp=open_file("rangegot.dat","w"))==NULL) return;  //write predict.dat
  if(p2r==0){
    // 入力をコピー
    for (t=0; t<n_total; t++) {
      for (k=0; k<=(n_channels); k++) {
	test->x[t][k] = givendata->x[t][k];
	//      test->X[t][k] = givendata->X[t][k];
      }
    }
    test->ymax = givendata->ymax;
    test->ymin = givendata->ymin;
    // 出力および誤差を計算
    x0=test->x[n_train][0];
    for (t=n_train; t<n_total; t++){
      if(fabs(net->cost)>1) calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
      else calc_output_c(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
      test->e[t] = (test->Y[t]-givendata->Y[t]);
      if(test->x[t][0]>x0) fprintf(fp,"\n"); x0=test->x[t][0];
      //    if(t>n_train && t< n_total && test->x[t][0]>test->x[t-1][0]) fprintf(fp,"\n");
      //    if(t>n_train && t< n_total && test->x[t][0]<test->x[t-1][0]) break;
      //    fprintf(fp,"%15.7e %e %e",test->Y[t],test->x[t][0],test->x[t][1]);
      //      fprintf(fp,"%+13.7e %+13.7e %+13.7e",test->Y[t],givendata->x[t][0],givendata->x[t][1]);
      fprintf(fp,"%+13.7e %+13.7e %+13.7e",givendata->x[t][0],givendata->x[t][1],test->Y[t]);
      for(i=-2;i<=8;i++){
	if(net->dc==i) fprintf(fp," %e",givendata->x[t][0]);
	else fprintf(fp," -");
      }
      //    fprintf(fp," %d",net->dc);
      //    if(net->dc==1 || net->dc==5) fprintf(fp," 0");
      //    else if(net->dc==4) fprintf(fp," 1");
      //    else if(net->dc==2 || net->dc==3) fprintf(fp," 2");
      //    else if(net->dc==0) fprintf(fp," 3");
      //    fprintf(fp," %d %d\n",net->dc,t);
      fprintf(fp,"\n");
      MSE += (e2=square(test->e[t]));
      /**/
    }
  }
  else if(p2r>0){//higher resolution
    FLOAT xxxx[3],yyyy,YYYY,dxxxx[2];
    int j,ii;
    xxxx[2]=1;
    dxxxx[0]=(net->xmax[0]-net->xmin[0])/p2r;
    dxxxx[1]=(net->xmax[1]-net->xmin[1])/p2r;
    for (j=0;j<p2r;j++){
      xxxx[1]=j*dxxxx[1]+net->xmin[1];
      for (i=0;i<p2r;i++){
	xxxx[0]=i*dxxxx[0]+net->xmin[0];
	if(fabs(net->cost)>1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
	else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
	//moverange1(xxxx[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]);
	//	fprintf(fp,"%13.7e %13.7e %13.7e",xxxx[0],xxxx[1],YYYY);
	fprintf(fp,"%13.7e %13.7e %13.7e",
		moverange1(xxxx[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]),
		moverange1(xxxx[1], net->xmin[1],net->xmax[1],net->xmin0[1],net->xmax0[1]),
		YYYY);

	for(ii=-2;ii<=8;ii++){
	  //	  if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
	  if(net->dc==ii) fprintf(fp," %e",moverange1(xxxx[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]));
	  else fprintf(fp," -");
	}
	fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
      fprintf(stdout,"."); fflush(stdout);
    }
  }
  fclose(fp);

  fp=open_file("weights.dat","w");
  for(i=0;i<net->n_cells;i++){
    fprintf(fp," %e %e\n",
	    moverange1(net->cell[i].w[0], net->xmin[0],net->xmax[0],net->xmin0[0],net->xmax0[0]),
	    moverange1(net->cell[i].w[1], net->xmin[1],net->xmax[1],net->xmin0[1],net->xmax0[1]));
  }
  fclose(fp);

  fp=open_file("rangegot.gpl","w");

  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 1:2:3\n");

  fprintf(fp,"set term x11\nset view 60,330\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
  fprintf(fp,"set nocontour\nset style data lines\n");
  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p pt 0, \"rangegot.dat\" using 1:3:2 w p pt 1\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");

  fprintf(fp,"splot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p pt 0, \"rangegot.dat\" using 1:3:2 w lp ls 1 pt 0\n",net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");

  fprintf(fp,"set contour\nset view 60,350\n");
  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3\n");
  fprintf(fp,"pause -1 \"Hit Return to continue.\"\n");
  //  fprintf(fp,"splot \"rangegot.dat\" using 1:2:3 w l\npause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
  fprintf(fp,"plot \"rangegot.dat\" u 5:2 t \"no opposite weight (open boundary)\" w p ps 1");
  fprintf(fp,", \"rangegot.dat\" u  6:2 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p ps 2");
  fprintf(fp,", \"rangegot.dat\" u  7:2 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p ps 3");
  fprintf(fp,", \"rangegot.dat\" u  8:2 t \"xi(fu renzoku?) with y0x for t<0\"          w p ps 4");
  fprintf(fp,", \"rangegot.dat\" u  9:2 t \"xi(   renzoku?) with y0x for t>1\"          w p ps 6");
  fprintf(fp,", \"rangegot.dat\" u 10:2 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p ps 0");
  fprintf(fp,", \"rangegot.dat\" u 11:2 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p ps 7");
  fprintf(fp,", \"rangegot.dat\" u 12:2 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p ps 8");
  fprintf(fp,", \"weights.dat\" u 1:2 t \"weight vecotr\"        w p ps 9");
  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l lt 1");
  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
  fclose(fp);
  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");


//  fp=open_file("rangegot1.gpl","w");
//  fprintf(fp,"set data style lines\nset key below\nset pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:3 t \"no opposite weight (open boundary)\" w p 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:3 t \"continuous    y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:3 t \"discontinuous y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:3 t \"discontinuous y0x for t<0\"          w p 2 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:3 t \"continuous    y0x for t>1\"          w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 10:3 t \"continuous    y1x for 0<t<1\"        w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 11:3 t \"discontinuous y0x for t< 0 dyx>2\"   w p 2 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:3 t \"discontinuous y0x for t>10 dyx>2\"   w p 2 8");
//  fprintf(fp,", \"rangegot.dat\" u 13:14 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+200 -e gnuplot rangegot1.gpl&");


  {
    //    int Hmin,Hmax,Vmin,Vmax;
    char cmd='r';
    FILE *fp1;
    FLOAT cH=0,cD=0,cV=0;
   int ViewV=60,ViewH=330;//,zahyohenkan
   
    printf("u,d,r,l,q for up,down,right,left,quit:");
    //    fp = open_pipe(GNUPLOT, "w");
    fp = open_pipe("gnuplot -geometry 510x360-0+0", "w");
    fp1 = open_pipe("gnuplot -geometry 640x420-0+460", "w");
    fprintf(fp,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
    fprintf(fp,"set nocontour\nset style data lines\n");
    fprintf(fp,"set label \". Camera\"  at first %e, first  %e, first  %e\n",cH,cD,cV);

    fprintf(fp1,"set term x11\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
    fprintf(fp1,"set nocontour\nset style data lines\n");
    fprintf(fp1,"set label \". Camera\"  at first %e, first  %e, first  %e\n",cH,cD,cV);
    for(i=0;;){
      if(i>0) fgets(line,256,stdin);
      i=1;
      if(line[0]=='q') break;
      if(line[0]=='u' || line[0]=='d' || line[0]=='r' || line[0]=='l') cmd=line[0];
      if(cmd=='u') ViewV+=10;
      else if(cmd=='d') ViewV-=10;
      else if(cmd=='r') ViewH+=10;
      else if(cmd=='l') ViewH-=10;

      ViewH = (ViewH+360)%360;
      ViewV =(ViewV+180)%180;
      fprintf(fp,"set view %d,%d\nsplot \"rangegot.dat\" using 1:3:2 w l 1 0\n",ViewV,ViewH);
      fprintf(fp1,"set view %d,%d\nsplot [%e:%e][][%e:%e] \"%s\" using 1:3:2 w p pt 0\n",ViewV,ViewH,net->xmin0[0],net->xmax0[0],net->xmin0[1],net->xmax0[1],givendata->path);
      fflush(fp);
      fflush(fp1);
    }
    close_pipe(fp);
    close_pipe(fp1);
  }

  MSE2=MSE1=MSE /= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;

  //  printf("=>MSEtest=%f\n",MSE);

  test->MSE2 =test->MSE1 =test->MSE  = MSE;
  test->NMSE = NMSE;
  return;
}

//void exec_ssp_test_r0(NET *net, DATA *givendata, DATA *test) {
//  //char name[32] = "exec_ssp_test";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_total = givendata->n_total;
//  int n_channels = givendata->k;
//  int n_train = givendata->n_train;
//  int n_test = givendata->n_test;
//  int t=0;
//  int k=0;
//  FLOAT y;
//  //  FLOAT p1,p2,p1_sum=0,p2_sum=0;
//  FLOAT e2;
//  FLOAT MSE1=0;
//  FLOAT MSE2=0;
//  //  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
//  FILE *fp;
//  int i;
//  char line[80];
//  int p2r;//transform the polar coordinate system to rectangle one (1), inverse(-1), non(0)
//  FLOAT x0;
//  printf("cos(theta) for search neighboring weight vector(-2 for y0x), p2r(1,0,-1): ");
//  fgets(line,80,stdin);
//  sscanf(line,"%lf%d",&(net->cost),&p2r);
//
//  if((fp=open_file("rangegot.dat","w"))==NULL) return;  //write predict.dat
//  if(p2r==0){
//    // 入力をコピー
//    for (t=0; t<n_total; t++) {
//	for (k=0; k<=(n_channels); k++) {
//	  test->x[t][k] = givendata->x[t][k];
//	  //      test->X[t][k] = givendata->X[t][k];
//	}
//    }
//    test->ymax = givendata->ymax;
//    test->ymin = givendata->ymin;
//    // 出力および誤差を計算
//    x0=test->x[n_train][0];
//    for (t=n_train; t<n_total; t++){
//	if(net->cost<-1) calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
//	else calc_output_c(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
//	test->e[t] = (test->Y[t]-givendata->Y[t]);
//	if(test->x[t][0]>x0) fprintf(fp,"\n"); x0=test->x[t][0];
//	//    if(t>n_train && t< n_total && test->x[t][0]>test->x[t-1][0]) fprintf(fp,"\n");
//	//    if(t>n_train && t< n_total && test->x[t][0]<test->x[t-1][0]) break;
//	//    fprintf(fp,"%15.7e %e %e",test->Y[t],test->x[t][0],test->x[t][1]);
//	fprintf(fp,"%13.7e %e %e",test->Y[t],givendata->x[t][0],givendata->x[t][1]);
//	for(i=-2;i<=8;i++){
//	  if(net->dc==i) fprintf(fp," %e",givendata->x[t][0]);
//	  else fprintf(fp," -");
//	}
//	//    fprintf(fp," %d",net->dc);
//	//    if(net->dc==1 || net->dc==5) fprintf(fp," 0");
//	//    else if(net->dc==4) fprintf(fp," 1");
//	//    else if(net->dc==2 || net->dc==3) fprintf(fp," 2");
//	//    else if(net->dc==0) fprintf(fp," 3");
//	//    fprintf(fp," %d %d\n",net->dc,t);
//	fprintf(fp,"\n");
//	MSE += (e2=square(test->e[t]));
//	/**/
//    }
//  }
//  else if(p2r==1){
//    FLOAT theta,zzz,lll,xxx,yyy,xxxx[4],yyyy,YYYY;
//    FLOAT xxxmin=1e30,xxxmax=-1e30,dxxx;
//    FLOAT zzzmin=1e30,zzzmax=-1e30,dzzz;
//    int ii,j;
//    for(i=0;i<4;i++){
//	if(i==0){
//	  theta=(90.0-net->xmin0[0])/Degree;
//	  zzz=net->ymax0*cos((90.0-net->xmin0[1])/Degree);
//	}
//	else if(i==1){
//	  theta=(90.0-net->xmin0[0])/Degree;
//	  zzz=net->ymax0*cos((90.0-net->xmax0[1])/Degree);
//	}
//	else if(i==2){
//	  theta=(90.0-net->xmax0[0])/Degree;
//	  zzz=net->ymax0*cos((90.0-net->xmin0[1])/Degree);
//	}
//	else if(i==3){
//	  theta=(90.0-net->xmax0[0])/Degree;
//	  zzz=net->ymax0*cos((90.0-net->xmax0[1])/Degree);
//	}
//	lll=sqrt(net->ymax0*net->ymax0-zzz*zzz);
//	xxx=lll*cos(theta);
//	yyy=lll*sin(theta);
//	if(xxx<xxxmin) xxxmin=xxx;
//	if(xxx>xxxmax) xxxmax=xxx;
//	if(zzz<zzzmin) zzzmin=zzz;
//	if(zzz>zzzmax) zzzmax=zzz;
//    }
//    dxxx=xxxmax-xxxmin;
//    dzzz=zzzmax-zzzmin;
//    xxxx[2]=1;
//    for (j=0;j<100;j++){
//	for (i=0;i<100;i++){
//	  xxx=i*dxxx+xxxmin;
//	  zzz=j*dzzz+zzzmin;
//	  yyy=
//	  xxxx[0]=moverange1(i*dxxx+xxxmin, net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]);
//	  xxxx[1]=moverange1(j*dzzz+zzzmin, net->xmin0[1],net->xmax0[1],net->xmin[1],net->xmax[1]);
//	  if(net->cost<-1) calc_output(net, xxxx,&yyyy,&y,&YYYY);
//	  else calc_output_c(net, xxxx,&yyyy,&y,&YYYY);
//	  fprintf(fp,"%13.7e %13.7e %13.7e",xxxx[0],xxxx[1],YYYY);
//	  for(ii=-2;ii<=8;ii++){
//	    if(net->dc==ii) fprintf(fp," %e",xxxx[0]);
//	    else fprintf(fp," -");
//	  }
//	  fprintf(fp,"\n");
//	}
//	fprintf(fp,"\n");
//    }
//  }
//
//  for(i=0;i<net->n_cells;i++){
//    for(t=0;t<3;t++) fprintf(fp," -");
//    for(t=-2;t<=6;t++) fprintf(fp," -");
//    fprintf(fp," %e %e\n",net->cell[i].w[0], net->cell[i].w[1]);
//  }
//  fclose(fp);
//  fp=open_file("rangegot.gpl","w");
//
//  fprintf(fp,"set term table\nset contour\nset view 0,0\nset nosurface\nset cntrparam levels auto 10\n");
//  fprintf(fp,"set output \"rangegotcnt.tbl\"\nsplot \"rangegot.dat\" using 2:3:1\n");
//  fprintf(fp,"set term x11\nset view 60,350\nset surface\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"splot \"rangegot.dat\" using 2:3:1 w l\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set view 0,0\nset nosurface\nreplot\npause -1 \"Hit Return to continue.\"\n");
//  fprintf(fp,"set pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:3 t \"no opposite weight (open boundary)\" w p 2 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:3 t \"xi(   renzoku?) with y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:3 t \"xi(fu renzoku?) with y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:3 t \"xi(fu renzoku?) with y0x for t<0\"          w p 1 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:3 t \"xi(   renzoku?) with y0x for t>1\"          w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 10:3 t \"xi(   renzoku?) with y1x for 0<t<1\"        w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 11:3 t \"xi(fu renzoku?) with y0x for t< 0 dyx>1\"   w p 3 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:3 t \"xi(fu renzoku?) with y0x for t>10 dyx>1\"   w p 3 8");
//  fprintf(fp,", \"rangegot.dat\" u 13:14 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,", \"rangegotcnt.tbl\" u 1:2 t \"contour\"        w l 3 1");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+50 -e gnuplot rangegot.gpl&");
//
//
////  fp=open_file("rangegot1.gpl","w");
////  fprintf(fp,"set data style lines\nset key below\nset pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
////  fprintf(fp,"plot \"rangegot.dat\" u 5:3 t \"no opposite weight (open boundary)\" w p 1");
////  fprintf(fp,", \"rangegot.dat\" u  6:3 t \"continuous    y0x for |y0x-y1x|<e2\" w p 1 2");
////  fprintf(fp,", \"rangegot.dat\" u  7:3 t \"discontinuous y0x for |y0x-y1x|>e2\" w p 1 3");
////  fprintf(fp,", \"rangegot.dat\" u  8:3 t \"discontinuous y0x for t<0\"          w p 2 4");
////  fprintf(fp,", \"rangegot.dat\" u  9:3 t \"continuous    y0x for t>1\"          w p 1 0");
////  fprintf(fp,", \"rangegot.dat\" u 10:3 t \"continuous    y1x for 0<t<1\"        w p 1 6");
////  fprintf(fp,", \"rangegot.dat\" u 11:3 t \"discontinuous y0x for t< 0 dyx>2\"   w p 2 7");
////  fprintf(fp,", \"rangegot.dat\" u 12:3 t \"discontinuous y0x for t>10 dyx>2\"   w p 2 8");
////  fprintf(fp,", \"rangegot.dat\" u 13:14 t \"weight vecotr\"        w p 3 9");
////  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
////  fclose(fp);
////  system("xterm -geometry 80x5-0+200 -e gnuplot rangegot1.gpl&");
//
//  MSE2=MSE1=MSE /= n_test;
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//
//  printf("=>MSEtest=%f\n",MSE);
//
//  test->MSE2 =test->MSE1 =test->MSE  = MSE;
//  test->NMSE = NMSE;
//  return;
//}

//void exec_ssp_test_r0(NET *net, DATA *givendata, DATA *test) {
//  //char name[32] = "exec_ssp_test";
//  FLOAT MSE = 0.0, NMSE = 0.0;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_total = givendata->n_total;
//  int n_channels = givendata->k;
//  int n_train = givendata->n_train;
//  int n_test = givendata->n_test;
//  int t=0;
//  int k=0;
//  FLOAT y;
//  //  FLOAT p1,p2,p1_sum=0,p2_sum=0;
//  FLOAT e2;
//  FLOAT MSE1=0;
//  FLOAT MSE2=0;
//  //  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
//  FILE *fp;
//  int i;
//
//  // 入力をコピー
//  for (t=0; t<n_total; t++) {
//    for (k=0; k<=(n_channels); k++) {
//	test->x[t][k] = givendata->x[t][k];
//	//      test->X[t][k] = givendata->X[t][k];
//    }
//  }
//  test->ymax = givendata->ymax;
//  test->ymin = givendata->ymin;
//  // 出力および誤差を計算
//  if((fp=open_file("rangegot.dat","w"))==NULL) return;  //write predict.dat
//  for (t=n_train; t<n_total; t++) {
//
//    calc_output_c(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
//    test->e[t] = (test->Y[t]-givendata->Y[t]);
//
//    if(t>n_train && t< n_total && test->x[t][0]>test->x[t-1][0]) fprintf(fp,"\n");
//    //    if(t>n_train && t< n_total && test->x[t][0]<test->x[t-1][0]) break;
//    //    fprintf(fp,"%15.7e %e %e",test->Y[t],test->x[t][0],test->x[t][1]);
//    fprintf(fp,"%15.7e %e %e",test->Y[t],givendata->X[t][0],givendata->X[t][1]);
//    for(i=-2;i<=8;i++){
//	if(net->dc==i) fprintf(fp," %e",givendata->X[t][0]);
//	else fprintf(fp," -");
//    }
//	  //    fprintf(fp," %d",net->dc);
//    //    if(net->dc==1 || net->dc==5) fprintf(fp," 0");
//    //    else if(net->dc==4) fprintf(fp," 1");
//    //    else if(net->dc==2 || net->dc==3) fprintf(fp," 2");
//    //    else if(net->dc==0) fprintf(fp," 3");
//    fprintf(fp," %d %d\n",net->dc,t);
//    MSE += (e2=square(test->e[t]));
//    /**/
//  }
//  for(i=0;i<net->n_cells;i++){
//    for(t=0;t<3;t++) fprintf(fp," - -");
//    for(t=-2;t<=8;t++) fprintf(fp," - -");
//    fprintf(fp," %e %e\n",net->cell[i].w[0], net->cell[i].w[1]);
//  }
//  fclose(fp);
//  fp=open_file("rangegot.gpl","w");
//  fprintf(fp,"set data style lines\nset key below\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"splot \"rangegot.dat\" using 2:3:1\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+0 -e gnuplot rangegot.gpl&");
//
//
//  fp=open_file("rangegot1.gpl","w");
//  fprintf(fp,"set data style lines\nset key below\nset pointsize 0.8\n");// fprintf(fp,"set data style lines\nset key outside\n");
//  fprintf(fp,"plot \"rangegot.dat\" u 5:3 t \"no opposite weight (open boundary)\" w p 1");
//  fprintf(fp,", \"rangegot.dat\" u  6:3 t \"continuous    y0x for |y0x-y1x|<e2\" w p 1 2");
//  fprintf(fp,", \"rangegot.dat\" u  7:3 t \"discontinuous y0x for |y0x-y1x|>e2\" w p 1 3");
//  fprintf(fp,", \"rangegot.dat\" u  8:3 t \"discontinuous y0x for t<0\"          w p 2 4");
//  fprintf(fp,", \"rangegot.dat\" u  9:3 t \"continuous    y0x for t>1\"          w p 1 0");
//  fprintf(fp,", \"rangegot.dat\" u 10:3 t \"continuous    y1x for 0<t<1\"        w p 1 6");
//  fprintf(fp,", \"rangegot.dat\" u 11:3 t \"discontinuous y0x for t< 0 dyx>2\"   w p 2 7");
//  fprintf(fp,", \"rangegot.dat\" u 12:3 t \"discontinuous y0x for t>10 dyx>2\"   w p 2 8");
//  fprintf(fp,", \"rangegot.dat\" u 13:14 t \"weight vecotr\"        w p 3 9");
//  fprintf(fp,"\npause -1 \"Hit Return to continue.\"\n");
//  fclose(fp);
//  system("xterm -geometry 80x5-0+200 -e gnuplot rangegot1.gpl&");
//
//  MSE2=MSE1=MSE /= n_test;
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//
//  printf("=>MSEtest=%f\n",MSE);
//
//  test->MSE2 =test->MSE1 =test->MSE  = MSE;
//  test->NMSE = NMSE;
//  return;
//}

void exec_ssp_test(NET *net, DATA *givendata, DATA *test) {
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  FLOAT p1,p2,p1_sum=0,p2_sum=0,e2;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
  FILE *fp;

  // 入力をコピー
  for (t=0; t<n_total; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  // 出力および訓練誤差を計算20191012
  //  if((fp=open_file("predict-tr.dat","w"))==NULL) return;  
  givendata->MSEtr=0.0;
  for (t=0;t<n_train; t++) {
    calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
    FLOAT e= (test->Y[t]-givendata->Y[t]);
    givendata->MSEtr += e*e;
//    fprintf(fp,"%.7e %d %.7e %.7e %.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var=0,c,v,e2--#1\n",
//	    test->Y[t],t+1,givendata->Y[t],test->y[t],
//	    0.,//net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
//	    net->c,//c
//	    (int)net->cell[net->c].vv,//number of data in V_c
//	    e*e,0.);
  }
  givendata->MSEtr/=n_train;
  //  fclose(fp);
  // 出力および誤差を計算
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  for (t=n_train; t<n_total; t++) {
    calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    //    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t-n_train+1,givendata->Y[t],test->y[t]);

    fprintf(fp,"%.7e %d %.7e %.7e %.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var=0,c,v,e2--#1\n",
	    test->Y[t],t-n_train+1,givendata->Y[t],test->y[t],
	    0.,//net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),0.);

    //    printf("%d)%f %f %f\n",t,test->y[t],y,test->Y[t]);
    MSE += (e2=square(test->e[t]));
    //    printf("%d)MSE%e,e%e tY%e gY%e\n",t,MSE,test->e[t],test->Y[t],givendata->Y[t]);//check060517
    /*bunpu=1*/
    p1=exp(-(square(test->x[t][0]-mu1)+square(test->x[t][1]-0.5))/2./sigma/sigma);
    p1_sum += p1;
    MSE1      += e2*p1; 
    /*bunpu=2*/
    p2=exp(-(square(test->x[t][0]-mu2)+square(test->x[t][1]-0.5))/2./sigma/sigma);
    p2_sum += p2;
    MSE2      += e2*p2; 
    /**/
  }
  fclose(fp);
  MSE /= n_test;
  MSE1/=p1_sum;
  MSE2/=p2_sum;
 
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  //  printf("=>MSEtest=%f\n",MSE);
  test->MSE  = MSE;
  test->MSE1 = MSE1;
  test->MSE2 = MSE2;
  test->NMSE = NMSE;
  {
    FILE *fp;
    fp=fopen("mse.dat","w");
    fprintf(fp,"%e %e %e #MSEtr, MSE, MSE",givendata->MSEtr,MSE,NMSE);
    //    fprintf(fp,"%e %e %e #MSE, NMSE, MSEtr",MSE,NMSE,givendata->MSEtr);
    fclose(fp);
  }
  return;
} //void exec_ssp_test(NET *net, DATA *givendata, DATA *test) {

void exec_ssp_test_ensemble(NET *net, DATA *givendata, DATA *test) {
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  //  FLOAT y;
  FLOAT p1,p2,p1_sum=0,p2_sum=0,e2;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
  FILE *fp;
  FLOAT _yt,_Yt,_ytsum,_Ytsum,_y;
  int j;
  double n_test1=n_test;
  // 入力をコピー
  for (t=0; t<n_total; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  // 出力および誤差を計算
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  for (t=n_train; t<n_total; t++) {
    _ytsum=_Ytsum=0;
    for(j=0;j<net[0].nEns;j++){
      calc_output(&net[j], test->x[t],&_yt,&_y,&_Yt);
      if(net[j].cell[net[j].c].v>1) {//v>1 for excluding calculation of MSE see my_plinn.c NG
	_Yt=givendata->Y[t];
	_yt=givendata->y[t];
	n_test1 -= (1./net[0].nEns);
      }
      _ytsum+=_yt;
      _Ytsum+=_Yt;
    }
    test->y[t]=_ytsum/net[0].nEns;
    test->Y[t]=_Ytsum/net[0].nEns;
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += (e2=square(test->e[t]));
    //    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t-n_train+1,givendata->Y[t],test->y[t]);
    fprintf(fp,"%15.7e %d %15.7e %15.7e %15.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var,c,v,e2#2\n",
	    test->Y[t],t-n_train+1,givendata->Y[t],test->y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),
	    0.);
    //    printf("%d)%f %f %f\n",t,test->y[t],y,test->Y[t]);
    /*bunpu=1*/
    p1=exp(-(square(test->x[t][0]-mu1)+square(test->x[t][1]-0.5))/2./sigma/sigma);
    p1_sum += p1;
    MSE1      += e2*p1; 
    /*bunpu=2*/
    p2=exp(-(square(test->x[t][0]-mu2)+square(test->x[t][1]-0.5))/2./sigma/sigma);
    p2_sum += p2;
    MSE2   += e2*p2; 
    /**/
  }
  fclose(fp);
  //  if(n_test1<=1) fprintf(stdout,"check n_test1=%e, MSE=%e\n",n_test1,MSE);
  MSE /= n_test1;
  MSE1/=p1_sum;
  MSE2/=p2_sum;

  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  //  printf("=>MSEtest=%f\n",MSE);
  test->MSE  = MSE;
  test->MSE1 = MSE1;
  test->MSE2 = MSE2;
  test->NMSE = NMSE;
  {
    FILE *fp;
    fp=fopen("mse.dat","w");
    fprintf(fp,"%e %e %e #MSEtr, MSE, MSE",givendata->MSEtr,MSE,NMSE);
    //    fprintf(fp,"%e %e %e #MSE, NMSE, MSEtr",MSE,NMSE,givendata->MSEtr);
    //    fprintf(fp,"%e %e #MSE, NMSE",MSE,NMSE);
    fclose(fp);
  }
  return;
} //void exec_ssp_test_ensemble(NET *net, DATA *givendata, DATA *test) {

void exec_ssp_test_NIPS04(NET *net, DATA *givendata, DATA *test) {
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  //  FLOAT p1,p2,p1_sum=0,p2_sum=0,e2;
  FLOAT e2;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  //  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
  FILE *fp;
  int i;
  // 入力をコピー
  for (t=0; t<n_total; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  ////////////////////////////予測の分散を求める
  for (i=0; i<net->n_cells;i++){
    net->cell[i].SS=0;
    net->cell[i].vv=0;
  }
  for (t=0; t<n_train; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    net->cell[net->c].SS += square(test->Y[t]-givendata->Y[t]);//square error of training
    net->cell[net->c].vv +=1;
    if (t == n_total-1) break;
  }
  ////////////////////////////
  // 出力および誤差を計算
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat for NIPS04 and IJCNN06
  for (t=n_train; t<n_total; t++){
    calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
#ifdef NEEDKYORI
    FLOAT kyori2testtrain;
    FLOAT kyori2testtraint;
    int tt;
    kyori2testtrain=1e30;
    for(tt=0;tt<n_train;tt++){
      kyori2testtraint=distance2(givendata->X[t],givendata->X[tt],n_channels);
      if(kyori2testtrain>kyori2testtraint) kyori2testtrain=kyori2testtraint;
    }
    fprintf(fp,"%15.7e %d %15.7e %15.7e %15.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var,c,v,e2#3\n",
	    test->Y[t],
	    t-n_train+1,
	    givendata->Y[t],
	    test->y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),
	    kyori2testtrain);
#else
    fprintf(fp,"%15.7e %d %15.7e %15.7e %15.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var,c,v,e2#4\n",
	    test->Y[t],
	    t-n_train+1,
	    givendata->Y[t],
	    test->y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),
	    0.);
#endif
    //    fprintf(fp,"%15.7e %d %15.7e\n",test->Y[t],t-n_train+1,givendata->Y[t]);
    //    printf("%d)%f %f %f\n",t,test->y[t],y,test->Y[t]);
    MSE += (e2=square(test->e[t]));
    /*bunpu=1*/
    MSE1      += net->cell[net->c].SS/net->cell[net->c].vv; //variance
    /**/
  }
  fclose(fp);
  MSE /= n_test;
  MSE1/= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;

  printf("==>MSEtest=%f, predictingMSEtest =%f\n",MSE,MSE1);
  test->MSE  = MSE;
  test->MSE1 = MSE1;
  test->MSE2 = MSE2;
  test->NMSE = NMSE;
  return;
}

void exec_ssp_test_NIPS04_OLD(NET *net, DATA *givendata, DATA *test) {
  //char name[32] = "exec_ssp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  //  FLOAT p1,p2,p1_sum=0,p2_sum=0,e2;
  FLOAT e2;
  FLOAT MSE1=0;
  FLOAT MSE2=0;
  //  FLOAT mu1=0.666666666667,mu2=0.333333333333,sigma=0.1;
  FILE *fp;
  int i;
  // 入力をコピー
  for (t=0; t<n_total; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  ////////////////////////////予測の分散を求める
  for (i=0; i<net->n_cells;i++){
    net->cell[i].SS=0;
    net->cell[i].vv=0;
  }
  for (t=0; t<n_train; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    net->cell[net->c].SS += square(test->Y[t]-givendata->Y[t]);
    net->cell[net->c].vv +=1;
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
  }
  ////////////////////////////
  // 出力および誤差を計算
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  for (t=n_train; t<n_total; t++) {
    calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    fprintf(fp,"%15.7e %15.7e %2d %2d %15.7e %f # #5\n",
	    test->Y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance
	    //	    net->cell[net->c].SS/net->cell[net->c].vv, //variance
	    net->c,//firing cell
	    (int)net->cell[net->c].vv,//firing cell
	    givendata->Y[t],
	    (test->e[t])*(test->e[t]));
    //    fprintf(fp,"%15.7e %d %15.7e\n",test->Y[t],t-n_train+1,givendata->Y[t]);
    //    printf("%d)%f %f %f\n",t,test->y[t],y,test->Y[t]);
    MSE += (e2=square(test->e[t]));
    /*bunpu=1*/
    MSE1      += net->cell[net->c].SS/net->cell[net->c].vv; //variance
    /**/
  }
  fclose(fp);
  MSE /= n_test;
  MSE1/= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;

  printf("==>MSEtest=%f, predictingMSEtest =%f\n",MSE,MSE1);
  test->MSE  = MSE;
  test->MSE1 = MSE1;
  test->MSE2 = MSE2;
  test->NMSE = NMSE;
  return;
}

/*====================================================================*
 * exec_msp_test
 * 多段予測を行う(時系列)
 *====================================================================*/
void exec_msp_test(NET *net, DATA *givendata, DATA *test,double err4terminate,double err4propagate){
  //テストデータの多段予測
  //char name[32] = "exec_msp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int _n_test=0;
  int t=0;
  int k=0;
  FLOAT y;
  FILE *fp;
  int tpD=givendata->tpD;//20131101
  int tr1=givendata->tr1;//20131101
  int tp0=givendata->tp0;//20131101
  int tpG=givendata->tpG;//20131105
// 入力をコピー
//  for (t=0; t<=n_train; t++) {
//    for (k=0; k<=(n_channels); k++) {
//	test->x[t][k] = givendata->x[t][k];
//    }
//  }
  t=n_train;//initial time for prediction
  if(t<n_total) for(k=0;k<=(n_channels);k++) test->x[t][k] = givendata->x[t][k];
  test->x[t][n_channels-1]+=err4propagate;//20150218
  //  test->x[t][n_channels]+=err4propagate;//20150218
//?  if(1==0 && tpG==0){//??
//?    for(k=1;k<tpD;k++){
//?      int tk;
//?      if((tk=t+k)>=n_total) break;
//?      test->x[tk][0]=givendata->x[tk][0];//20131101??<=tpD
//?    }
//?  }
  //  for(k=1;k<=tpD;k++)          test->x[t+k][0]=givendata->x[t+k][0];//20131101?
  //  for (k=0; k<k2;k++) test->x[t][k] = 0;//ARMA
  //    for (k=0; k<k2;k++) test->x[t][k] = givendata->y[t+k2-k];//inchiki
  //  for (k=0; k<k2;k++) test->x[t][k] = t/5000.;//Jikantest
  //  printf("k2=%d\n",k2);
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  if((fp=open_file("msp.dat","w"))==NULL) return;  //write predict.dat
  //  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  //  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  // 出力および誤差を計算
  for (t=n_train; t<n_total; t++) {
//    {//for check
//      int i;
//      //      fprintf(stdout,"check tx[%d]",t);for(i=0;i<n_channels;i++) fprintf(stdout,"%e ",test->x[t][i]);fprintf(stdout,"\n");
//      //      fprintf(stdout,"check gx[%d]",t);for(i=0;i<n_channels;i++) fprintf(stdout,"%e ",givendata->x[t][i]);fprintf(stdout,"\n");
//    }
    calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += square(test->e[t]);
    _n_test++;
    fprintf(fp,"%.7e %d %.7e %.7e %.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var=0,c,v,e2--#1\n",
	    test->Y[t],t-n_train+1,givendata->Y[t],test->y[t],
	    0.,//net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),0.);

#ifdef DEBUG
    // Debug
    printf("t %d | ", t);
    for (k=0; k<=(n_channels); k++) printf("x[%d]=%3.2e, ", k, test->x[t][k]);
    printf("y=%3.2e\n", test->y[t]);
#endif
    //    int tp=t+tpD;
    if(err4terminate>0 && fabs(test->e[t])>err4terminate){
      fprintf(stderr,"#t=%d err|%g|>err4terminate%g\n",t,test->e[t],err4terminate);
      break;
    }
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//
    //    test->x[t+1][0] = test->y[t];
    //    for (k=1; k<n_channels; k++){
    //	test->x[t+1][k] = test->x[t][k-1];
    //    }
    //    test->x[t+tpD+1][0]=test->y[t];
    if(tpG!=0){//use givendata 
      test->x[t+1][0] =givendata->x[t+1][0];
    }
    else{//use prediction for unknown x
      //  t \in [n_train,n_total)
      //givendata->y[t]    =y00[t-n_train +tp0]; for t \in [0,n_total-n_train)
      //givendata->x[t][0] =y00[t-n_train +tp0 -tpD -1];
      //      if(t-n_train>=tr1-tp0 && (t+tpD+1<n_total)) test->x[t+tpD+1][0]=test->y[t];
      if(t-n_train<tr1-tp0) test->x[t+tpD+1][0]=givendata->x[t+tpD+1][0];//use given data if t in [tr0,tr1)
      else if(t+tpD+1<n_total) test->x[t+tpD+1][0]=test->y[t];
    }
    
    //    if((t+tpD+1>=n_train) && (t+tpD+1<n_total)) test->x[t+tpD+1][0]=test->y[t];//>=???20131102
    //    if((t>=n_train+tr1-tp0) && (t+tpD+1<n_total)) test->x[t+tpD+1][0]=test->y[t];//>=???20131102
    //    if(t>n_train+tr1-tp0) test->x[t+tpD+1][0]=test->y[t];//>=???20131102
    //    if(t>=n_train+tr1-tp0) test->x[t+tpD+1][0]=test->y[t];
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
    }
//    for (k=0; k<n_channels; k++){
//      if(k==k2) test->x[t+1][k] = test->y[t];//k2=0
//      else if(k<k2) test->x[t+1][k] = 0; //ARMA
//      else test->x[t+1][k] = test->x[t][k-1];
//      //else if(k<k2) test->x[t+1][k] = (t+1.)/5000.;//Jikantest
//      //else if(k<k2) test->x[t+1][k] = givendata->y[t+1+k2-k];//inchiki
//    }
  }//for (t=n_train; t<n_total; t++) {
  fclose(fp);
  MSE /= _n_test;
  NMSE = MSE/givendata->VARtest;
  //  NMSE = MSE/givendata->VARtrain;//060111
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  {
    FILE *fp;
    fp=fopen("mse.dat","w");
    //    fprintf(fp,"%d %e %e #MSE,NMSE",GlobalTime,MSE,NMSE);
    fprintf(fp,"%e %e %e #MSEtr, MSE, MSE",givendata->MSEtr,MSE,NMSE);
    //    fprintf(fp,"%e %e %e #MSE, NMSE, MSEtr",MSE,NMSE,givendata->MSEtr);
    fclose(fp);
  }
  return;
}//void exec_msp_test(NET *net, DATA *givendata, DATA *test,double err4terminate,double err4propagate){

void exec_msp_test_ensemble(NET *net, DATA *givendata, DATA *test){
  //テストデータのアンサンブル多段予測
  //char name[32] = "exec_msp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  int j;
  FLOAT _yt,_Yt,_ytsum,_Ytsum,_y;
  FILE *fp;
  // 入力をコピー
  t=n_train;
  for (k=0; k<=(n_channels); k++) {
    test->x[t][k] = givendata->x[t][k];
  }
  for (k=0; k<k2;k++) test->x[t][k] = 0;//ARMA
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  // 出力および誤差を計算
  for (t=n_train; t<n_total; t++) {
    _ytsum=_Ytsum=0;
    for(j=0;j<net[0].nEns;j++){
      calc_output(&net[j], test->x[t],&_yt,&_y,&_Yt);
      _ytsum+=_yt;
      _Ytsum+=_Yt;
    }
    test->y[t]=_ytsum/net[0].nEns;
    test->Y[t]=_Ytsum/net[0].nEns;
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += square(test->e[t]);
    //    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t-n_train+1,givendata->Y[t],test->y[t]);
    fprintf(fp,"%15.7e %d %15.7e %15.7e %15.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var,c,v,e2#6\n",
	    test->Y[t],t-n_train+1,givendata->Y[t],test->y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),
	    0.);

#ifdef DEBUG
    // Debug
    printf("t %d | ", t);
    for (k=0; k<=(n_channels); k++) printf("x[%d]=%3.2e, ", k, test->x[t][k]);
    printf("y=%3.2e\n", test->y[t]);
#endif
    if (t == n_total-1) break;

    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    for (k=0; k<n_channels; k++){
      if(k==k2) test->x[t+1][k] = test->y[t];
      else if(k<k2) test->x[t+1][k] = 0; //ARMA
      else test->x[t+1][k] = test->x[t][k-1];
    }
  }
  fclose(fp);
  MSE /= n_test;
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  {
    FILE *fp;
    fp=fopen("mse.dat","w");
    fprintf(fp,"%e %e %e #MSEtr, MSE, MSE",givendata->MSEtr,MSE,NMSE);
    //    fprintf(fp,"%e %e #MSE, NMSE",MSE,NMSE);
    fclose(fp);
  }
  return;
}
void exec_msp_test_Ensemble(NET *net, DATA *givendata, DATA *test){
  //テストデータの多段予測アンサンブル
  //char name[32] = "exec_msp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  int j;
  FLOAT *ytsum,*Ytsum,_y;
  FILE *fp;
  int nEns=net[0].nEns;
  // 入力をコピー
  ytsum=(FLOAT *)malloc(sizeof(FLOAT)*(n_total));
  Ytsum=(FLOAT *)malloc(sizeof(FLOAT)*(n_total));
  
  for (t=n_train; t<n_total; t++) ytsum[t]=Ytsum[t]=0;
  
  t=n_train; for (k=0; k<=(n_channels); k++)  test->x[t][k] = givendata->x[t][k];
  for (k=0; k<k2;k++) test->x[t][k] = 0;//ARMA
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  // 出力および誤差を計算
  
  for(j=0;j<nEns;j++){
    for (t=n_train; t<n_total; t++) {
      calc_output(&net[j], test->x[t],&test->y[t],&_y,&test->Y[t]);
      ytsum[t]+=test->y[t];
      Ytsum[t]+=test->Y[t];
      
#ifdef DEBUG
      // Debug
      printf("t %d | ", t);
      for (k=0; k<=(n_channels); k++) printf("x[%d]=%3.2e, ", k, test->x[t][k]);
      printf("y=%3.2e\n", test->y[t]);
#endif
      if (t == n_total-1) break;
      
      // 次の時刻での入力データを準備
      test->x[t+1][n_channels] = BIAS;//      test->x[t+1][n_channels] = 1.0;
      for (k=0; k<n_channels; k++){
	if(k==k2) test->x[t+1][k] = test->y[t];
	else if(k<k2) test->x[t+1][k] = 0; //ARMA
	else test->x[t+1][k] = test->x[t][k-1];
      }
    }
  }
  if((fp=open_file("predict.dat","w"))==NULL) return;  //write predict.dat
  for (t=n_train; t<n_total; t++){
    test->y[t]=ytsum[t]/nEns;
    test->Y[t]=Ytsum[t]/nEns;
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += square(test->e[t]);
    //    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t-n_train+1,givendata->Y[t],test->y[t]);
    fprintf(fp,"%15.7e %d %15.7e %15.7e %15.7e %2d %2d %.7e %.7e #Y^,t,Y,y,var,c,v,e2#7\n",
	    test->Y[t],t-n_train+1,givendata->Y[t],test->y[t],
	    net->cell[net->c].SS/net->cell[net->c].vv, //variance in V_c
	    net->c,//c
	    (int)net->cell[net->c].vv,//number of data in V_c
	    (test->e[t])*(test->e[t]),
	    0.);
  }
  fclose(fp);
  MSE /= n_test;
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  free(Ytsum);
  free(ytsum);
  return;
}

void exec_msp_test_IJCNN04(NET *net, DATA *givendata, DATA *test){
  //テストデータの多段予測
  //char name[32] = "exec_msp_test";
  //  FLOAT MSE = 0.0, NMSE = 0.0;
  FLOAT MSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  //  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
#define NT 20
  FLOAT mse[NT],e2;
  FLOAT msedp[NT],e2_p;
  int j,t0,nj=0;
  FLOAT *ds=givendata->ijcnn04data->ds;//data from smooth.dat
  FLOAT *dt=givendata->ijcnn04data->dt;//data of data.txt
  FLOAT *dp=givendata->ijcnn04data->dp;//data of prediction write out to predict.dat
  FLOAT *dl=givendata->ijcnn04data->dl;//data of linear approximation with data.txt and data.dat
  int t_t1=givendata->ijcnn04data->t_t1;//
  int _t,td=-(n_train-t_t1)+givendata->i_testblock*1000;
  FLOAT MSEdl;
  FLOAT MSEds;
  FLOAT MSEdp;
  FLOAT MSEsum=0;
  int n_MSEsum=0;
  FLOAT NMSEsum=0;
  FLOAT a,b;
  int  _t0;//real start time for prediction in all data
  int  _t1;//real start time for prediction in all data

  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  //  printf("\nNMSE");
  for(j=NT-1;j>=0;j--){
  //for(j=NT-1;j>=1;j--){
  //  for(j=0;j<NT;j++){
    t0=n_train+j;
    if(t_t1+j>981) continue;
    nj++;
    msedp[j]=mse[j]=0;
    //    printf("\n[%d]y=%5.3e=",t0+td,givendata->y[t0]);
    for (k=0; k<=(n_channels); k++) {
      test->x[t0][k] = givendata->x[t0][k];
      //      printf("%5.3e ",test->x[t0][k]);
    }
    for (k=0; k<k2;k++) test->x[t0][k] = 0;//ARMA
    
    // 出力および誤差を計算
    _t0=t0+td;//real start time for prediction in all data
    _t1=_t0+n_test;//real start time for prediction in all data
    a=dt[_t0-1];
    b=(dt[_t1]-dt[_t0-1])/(n_test+1);
    MSEdp=MSEdl=MSEds=0;
    for (t=t0; t<t0+n_test; t++) {
      _t=t+td;
      calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
      test->e[t] = (test->Y[t]-givendata->Y[t]);
      e2=square(test->e[t]);//smooth_.datに対する誤差
      mse[j] += e2;
#define ShortTermPred  //040917
#ifdef ShortTermPred
      if(j<5) dp[_t]=ds[_t]+(test->Y[k]); else dp[_t]=ds[_t];//>>MSE=4.406258e+02 is of real.dat (num=100)
      //      dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*8/(n_test-1.));//040917 MSE=442.628740 
      //dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*3/(n_test-1.));//040917 MSE=4.449830e+02
      //      dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j/(n_test-1.));//040917 MSE=453.371558 
      //dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*4/(n_test-1.));//040917 MSE=4.437768e+02
      //dp[_t]=ds[_t]+(test->Y[k])*(1.-(double)j/(n_test-1.));//040917  MSE=4.510872e+02
      //      dp[_t]=ds[_t];//040917   MSE=4.442810e+02    //      dp[_t]=ds[_t]+test->Y[k];
#else
      dp[_t]=ds[_t]+(test->Y[k]);
#endif
      //      dp[_t]=ds[_t]+test->Y[t];
      //      dp[_t]=ds[_t]+test->Y[t];
      //      printf("dp[%d]=%e,dt=%e\n",_t,dp[_t],dt[_t]);
      dl[_t]=a+(_t-_t0+1)*b; //linear approximation
      MSEdp +=square(dp[_t]-dt[_t]); //
      MSEdl +=square(dl[_t]-dt[_t]); //linear approximation
      MSEds +=square(ds[_t]-dt[_t]); //smooth approximation
      e2_p=square(dp[_t]-dt[_t]);//data.datに対する誤差
      msedp[j] += e2_p;
      //      printf("--->>>%+15.7e %d \n",dp[_t],_t);//check      
      // 次の時刻での入力データを準備
      //      if(t+1>=n_total) break;//???
      //      if(_t>=5001) break;//???
      if(t_t1+j+t-t0>=1000) break;//???
      test->x[t+1][n_channels] = BIAS;//      test->x[t+1][n_channels] = 1.0;
      for (k=0; k<n_channels; k++){
	if(k==k2) test->x[t+1][k] = test->y[t];
	else if(k<k2) test->x[t+1][k] = 0; //ARMA
	else test->x[t+1][k] = test->x[t][k-1];
      }
    }
    mse[j]/=n_test;
    msedp[j]/=n_test;
    //calc NMSEdl
    {
      MSEdl/=n_test;
      MSEds/=n_test;
      MSEdp/=n_test;
      if(j>=0){
	MSEsum += MSEdp;n_MSEsum++;
	NMSEsum += MSEdp/MSEdl;
      }
      if(MSEdl<1e-30) MSEdl=-MSE;
      //      if((t_t1+j==961) || (j==0 && !(961>=t_t1 && 961<=t_t1+NT-1))) {
      if(j==0) {
	test->MSE=MSEdp;
	givendata->VARtest=MSEdl;//for check
      }
      //      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSE=%6.4f",t_t1+n_test+j-1,MSEdp/MSEdl,MSEdp/MSEds,MSEdp);
      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSEdp=%8.3e MSEmean=%8.3e NMSE_Lmean=%8.3e",t_t1+j,MSEdp/MSEdl,MSEdp/MSEds,MSEdp,MSEsum/n_MSEsum,NMSEsum/n_MSEsum);
      //      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSE=%6.4f",_t0+j,MSEdp/MSEdl,MSEdp/MSEds,MSEdp);

      //      msedp[j] /=(MSEdl);
      //	if(_t0+j>981) printf("?");
      if(t_t1+j>981) printf("?");
	else if(t_t1+j==981) printf("*");
	else printf(" ");
	printf("\n");
    }
    //    printf("NMSE%d-%d=%4.2f ",t_t1+j,t_t1+j+n_test-1,mse[j]/n_test/givendata->VARtest);
  }
  printf("\n");
  //  t=givendata->i_testblock*1000+t_t1;

//  printf("(NMSE  %d-%d..%d-%d:",t,t+n_test-1,t+nj-1,t+n_test-1+nj-1);
//  for(j=0;j<nj;j++){
//    printf("%5.3f",mse[j]/givendata->VARtest);
//    if(t_t1+n_test+j-1>=981) printf("?");
//    else if(t_t1+n_test+j-1==980) printf("*");
//    else printf(" ");
//  }
//
//  printf(")\n(NMSE_L%d-%d..%d-%d:",t,t+n_test-1,t+nj-1,t+n_test-1+nj-1);
//  for(j=0;j<nj;j++){
//    printf("%5.3f",msedp[j]);
//    if(t_t1+n_test+j-1>=981) printf("?");
//    else if(t_t1+n_test+j-1==980) printf("*");
//    else printf(" ");
//  }
  //  printf(")\n");
  //  MSE=mse[0];
  //  MSE /= (n_test);
  //  NMSE = MSE/givendata->VARtest;
  //  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  //  test->MSE = MSE;
  //  test->NMSE = NMSE;
  givendata->ijcnn04data->MSEdp = MSEdp;//ueno(Mon Feb 16 23:07:57 2004)
  givendata->ijcnn04data->MSEdl = MSEdl;
  givendata->ijcnn04data->MSEds = MSEds;
  //  test->MSE = NMSEsum/n_MSEsum;//kuro 
  test->MSE =givendata->ijcnn04data->MSEmean = MSEsum/n_MSEsum;//kuro 
  //givendata->ijcnn04data->MSEdp = MSEsum/n_MSEsum;//kuro 
  return;
}

#ifdef OKASHII
void exec_msp_test_IJCNN04old(NET *net, DATA *givendata, DATA *test){
  //テストデータの多段予測
  //char name[32] = "exec_msp_test";
  //  FLOAT MSE = 0.0, NMSE = 0.0;
  FLOAT MSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
#define NT 30
  FLOAT mse[NT],e2;
  FLOAT msedp[NT],e2_p;
  int j,t0,nj=0;
  FLOAT *ds=givendata->ijcnn04data->ds;//data from smooth.dat
  FLOAT *dt=givendata->ijcnn04data->dt;//data of data.txt
  FLOAT *dp=givendata->ijcnn04data->dp;//data of prediction write out to predict.dat
  FLOAT *dl=givendata->ijcnn04data->dl;//data of linear approximation with data.txt and data.dat
  int t_t1=givendata->ijcnn04data->t_t1;//
  int _t,td=-(n_train-t_t1)+givendata->i_testblock*1000;
  FLOAT MSEdl;
  FLOAT MSEds;
  FLOAT MSEdp;
  FLOAT MSEsum=0;
  int n_MSEsum=0;

  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  //  printf("\nNMSE");
  for(j=NT-1;j>=0;j--){
  //for(j=NT-1;j>=1;j--){
  //  for(j=0;j<NT;j++){
    t0=n_train+j;
    if(t_t1+j>981) continue;
    nj++;
    msedp[j]=mse[j]=0;
    //    printf("\n[%d]y=%5.3e=",t0+td,givendata->y[t0]);
    for (k=0; k<=(n_channels); k++) {
      test->x[t0][k] = givendata->x[t0][k];
      //      printf("%5.3e ",test->x[t0][k]);
    }
    for (k=0; k<k2;k++) test->x[t0][k] = 0;//ARMA
    
    // 出力および誤差を計算
    for (t=t0; t<t0+n_test; t++) {
      calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
      test->e[t] = (test->Y[t]-givendata->Y[t]);
      e2=square(test->e[t]);//smooth_.datに対する誤差
      mse[j] += e2;
      _t=t+td;
#define ShortTermPred  //040917
#ifdef ShortTermPred
      if(j<5) dp[_t]=ds[_t]+(test->Y[k]); else dp[_t]=ds[_t];//>>MSE=4.406258e+02 is of real.dat (num=100)
      //      dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*8/(n_test-1.));//040917 MSE=442.628740 
      //dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*3/(n_test-1.));//040917 MSE=4.449830e+02
      //      dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j/(n_test-1.));//040917 MSE=453.371558 
      //dp[_t]=ds[_t]+(test->Y[k])*exp(-(double)j*4/(n_test-1.));//040917 MSE=4.437768e+02
      //dp[_t]=ds[_t]+(test->Y[k])*(1.-(double)j/(n_test-1.));//040917  MSE=4.510872e+02
      //      dp[_t]=ds[_t];//040917   MSE=4.442810e+02    //      dp[_t]=ds[_t]+test->Y[k];
#else
      dp[_t]=ds[_t]+(test->Y[k]);
#endif
      //      dp[_t]=ds[_t]+test->Y[t];
      //      dp[_t]=ds[_t]+test->Y[t];
      e2_p=square(dp[_t]-dt[_t]);//data.datに対する誤差
      msedp[j] += e2_p;
      //      printf("--->>>%+15.7e %d \n",dp[_t],_t);//check      
      // 次の時刻での入力データを準備
      if(t+1>=n_total) break;//???
      test->x[t+1][n_channels] = BIAS;//      test->x[t+1][n_channels] = 1.0;
      for (k=0; k<n_channels; k++){
	if(k==k2) test->x[t+1][k] = test->y[t];
	else if(k<k2) test->x[t+1][k] = 0; //ARMA
	else test->x[t+1][k] = test->x[t][k-1];
      }
    }
    mse[j]/=n_test;
    msedp[j]/=n_test;
    //calc NMSEdl
    {
      FLOAT a,b;
      int  _t0=t0+td;//real start time for prediction in all data
      int  _t1=_t0+n_test;//real start time for prediction in all data
      int  _t;//
      a=dt[_t0-1];
      b=(dt[_t1]-dt[_t0-1])/(n_test+1);
      MSEdp=MSEdl=MSEds=0;
      for(_t=_t0;_t<_t1;_t++){
	dl[_t]=a+(_t-_t0+1)*b; //linear approximation
	MSEdp +=square(dp[_t]-dt[_t]); //
	MSEdl +=square(dl[_t]-dt[_t]); //linear approximation
	MSEds +=square(ds[_t]-dt[_t]); //smooth approximation
	//	printf("%+15.7e %d %15.7e %15.7e %15.7e,a=%e,b=%e\n",dp[_t],_t,ds[_t],dt[_t],dl[_t],a,b);//check
      }
      MSEdl/=n_test;
      MSEds/=n_test;
      MSEdp/=n_test;
      MSEsum += MSEdp;n_MSEsum++;
      if(MSEdl<1e-30) MSEdl=-MSE;
      //      if((t_t1+j==961) || (j==0 && !(961>=t_t1 && 961<=t_t1+NT-1))) {
      if(j==0) {
	test->MSE=MSEdp;
	givendata->VARtest=MSEdl;//for check
      }
      //      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSE=%6.4f",t_t1+n_test+j-1,MSEdp/MSEdl,MSEdp/MSEds,MSEdp);
      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSEdp=%6.4f MSEsum=%12.6e",t_t1+j,MSEdp/MSEdl,MSEdp/MSEds,MSEdp,MSEsum/n_MSEsum);
      //      printf("%d:NMSE_L=%6.4f NMSE_S=%6.4f MSE=%6.4f",_t0+j,MSEdp/MSEdl,MSEdp/MSEds,MSEdp);

      //      msedp[j] /=(MSEdl);
      //	if(_t0+j>981) printf("?");
      if(t_t1+j>981) printf("?");
	else if(t_t1+j==981) printf("*");
	else printf(" ");
	printf("\n");
    }
    //    printf("NMSE%d-%d=%4.2f ",t_t1+j,t_t1+j+n_test-1,mse[j]/n_test/givendata->VARtest);
  }
  printf("\n");
  //  t=givendata->i_testblock*1000+t_t1;

//  printf("(NMSE  %d-%d..%d-%d:",t,t+n_test-1,t+nj-1,t+n_test-1+nj-1);
//  for(j=0;j<nj;j++){
//    printf("%5.3f",mse[j]/givendata->VARtest);
//    if(t_t1+n_test+j-1>=981) printf("?");
//    else if(t_t1+n_test+j-1==980) printf("*");
//    else printf(" ");
//  }
//
//  printf(")\n(NMSE_L%d-%d..%d-%d:",t,t+n_test-1,t+nj-1,t+n_test-1+nj-1);
//  for(j=0;j<nj;j++){
//    printf("%5.3f",msedp[j]);
//    if(t_t1+n_test+j-1>=981) printf("?");
//    else if(t_t1+n_test+j-1==980) printf("*");
//    else printf(" ");
//  }
  //  printf(")\n");
  //  MSE=mse[0];
  //  MSE /= (n_test);
  //  NMSE = MSE/givendata->VARtest;
  //  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  //  test->MSE = MSE;
  //  test->NMSE = NMSE;
  givendata->ijcnn04data->MSEdp = MSEdp;//ueno(Mon Feb 16 23:07:57 2004)
  givendata->ijcnn04data->MSEdl = MSEdl;
  givendata->ijcnn04data->MSEds = MSEds;
  test->MSE = MSEsum/n_MSEsum;//kuro 
  //givendata->ijcnn04data->MSEdp = MSEsum/n_MSEsum;//kuro 
  return;
}
#endif
void exec_msp_test_IJCNN04_out(NET *net, DATA *givendata, DATA *test){
  //テストデータの多段予測
  //char name[32] = "exec_msp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  FLOAT mse[NT],e2;
  int j,t0,nj=0,t000;
  FILE *fp;
  int t_t1=givendata->ijcnn04data->t_t1;//
  FLOAT *ds=givendata->ijcnn04data->ds;//data from smooth.dat
  //  FLOAT *dr=givendata->ijcnn04data->dr;//data from smooth.dat
  FLOAT *dp=givendata->ijcnn04data->dp;//data of prediction write out to predict.dat
  FLOAT *dt=givendata->ijcnn04data->dt;//data of data.txt
  FLOAT *dl=givendata->ijcnn04data->dl;//data of linear approximation with data.txt and data.dat
  FLOAT *dd=givendata->ijcnn04data->dd;//data of data.dat which is from dataconv0 or dataconv0_
  char *fname_p=givendata->ijcnn04data->fname_p;//predicted data
  char *fname_g=givendata->ijcnn04data->fname_g;//gnuplot 
  char *fname_d=givendata->ijcnn04data->fname_d;//gnuplot 
  FLOAT MSEdl;
  FLOAT MSEds;
  FLOAT MSEdp;

  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  t000=givendata->i_testblock*1000+t_t1;//real start time for prediction in all data
  for(j=0;j>=0;j--){
    t0=n_train+j;
    if(t_t1+j>981) continue;
    nj++;
    mse[j]=0;
    for (k=0; k<=(n_channels); k++) {
      test->x[t0][k] = givendata->x[t0][k];
    }
    for (k=0; k<k2;k++) test->x[t0][k] = 0;//ARMA
    
    // 出力および誤差を計算
    for (t=t0; t<t0+n_test; t++) {
      calc_output(net, test->x[t],&test->y[t],&y,&(test->Y[t]));
      //      printf("%e %d %e\n",test->Y[t],t,test->y[t]);
      test->e[t] = (test->Y[t]-givendata->Y[t]);
      e2=square(test->e[t]);
      mse[j] += e2;
      // 次の時刻での入力データを準備
      if(t+1>=n_total) break;//???if(t_t1+j+1>981) break;//???
      test->x[t+1][n_channels] = BIAS;//      test->x[t+1][n_channels] = 1.0;
      for (k=0; k<n_channels; k++){
	if(k==k2) test->x[t+1][k] = test->y[t];
	else if(k<k2) test->x[t+1][k] = 0; //ARMA
	else test->x[t+1][k] = test->x[t][k-1];
      }
    }
    mse[j]/=n_test;
  }
  //write predict.dat
  {
    //    sprintf(fname_p,"predict%s.dat",givendata->i_testblock);
    sprintf(fname_p,"predict.dat");
    if((fp=open_file(fname_p,"w"))==NULL) return;
    fprintf(fp,"#y^ t err dy^\n");
//    t0=n_train;
//    for(j=0;j<n_test;j++){
//	k=n_train+j;
//	t=t000   +j;
//	dp[t]=ds[t]+test->Y[k];
//    }
  }
  {//calc NMSE, MSEdl, MSEds
    FLOAT a,b;
    int  _t0=t000;//real start time for prediction in all data
    int  _t1=_t0+n_test;//real start time for prediction in all data
    int  _t;//
    a=dt[_t0-1];
    b=(dt[_t1]-dt[_t0-1])/(n_test+1);
    MSEdp=MSEdl=MSEds=0;
    //    for(_t=_t0;_t<_t1;_t++){
    for(j=0;j<n_test;j++){
      k=n_train+j;
      _t=_t0+j;
      dp[_t]=ds[_t]+test->Y[k];
      dl[_t]=a+(_t-_t0+1)*b; //linear approximation
      //      dp[_t] = (dp[_t]+dl[_t])/2.;//????????????????????????
      MSEdp +=square(dp[_t]-dt[_t]);
      MSEds +=square(ds[_t]-dt[_t]); //smooth approximation
      MSEdl +=square(dl[_t]-dt[_t]); //linear approximation
      //	printf("MSEdl[%d]=%e\n",_t,MSEdl);
      fprintf(fp,"%+15.7e %d %15.7e %15.7e %15.7e %15.7e\n",dp[_t],_t,ds[_t],dt[_t],dl[_t],dd[t]);
    }
    MSEdp/=n_test;
    MSEdl/=n_test;
    MSEds/=n_test;
    if(MSEdl<1e-20) {printf("MSEdl is ZERO.\n");MSEdl=-MSE;}
    //    printf("NMSE=%e=%e/%e\n",MSE/MSEdl,MSE,MSEdl);
    printf("--->NMSE_L=%8.6f NMSE_S=%8.6f MSEdp=%e\n",MSEdp/MSEdl,MSEdp/MSEds,MSEdp);
    NMSE=MSEdp/MSEdl;
    test->MSE = MSEdp;
    test->NMSE = MSEdp/MSEdl;
  }
  fclose(fp);
  //write predict.gp
  //  sprintf(fname_g,"predict%d.gp",givendata->i_testblock);
  sprintf(fname_g,"predict.gp");
  if((fp=open_file(fname_g,"w"))==NULL) return;
  fprintf(fp, "set style data linespoints\nplot [%d:%d] \"data.txt\" using 1:2, \"smooth.dat\" using 5:1, \"%s\" using 2:1, \"%s\" using 2:5 t \"Linear Prediction\", \"%s\" using 2:1\npause -1 \"Hit Return to continue.\"\n",givendata->i_testblock*1000+t_t1-10, givendata->i_testblock*1000+t_t1+n_test+10,fname_p,fname_p,fname_d);

//  fprintf(fp, "set data style linespoints\n
//plot [%d:%d] \"%s\" using 2:4 t 'data.txt~', \"smooth.dat\" using 5:1, \"%s\" using 2:1, \"%s\" using 2:5 t \"Linear Prediction\"\n
//pause -1 \"Hit Return to continue.\"\n
//",givendata->i_testblock*1000+t_t1-10, givendata->i_testblock*1000+t_t1+n_test+10,fname_p,fname_p,fname_p);
  fclose(fp);
  //exec gnuplot
#if GPLT >= 1
  {
    char line[512];
    sprintf(line,"xterm -geometry 80x5-0+0 -e gnuplot %s&",fname_g);
    system(line);
  }
#endif
  return;
}

void exec_msp_test_can2_2(NET *net, DATA *givendata, DATA *test){
  //CAN2-2 形式で予測を行う。
  //char name[32] = "exec_msp_test";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int k2 = net->k2;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  //  FLOAT y;
  int K=5,i,tt,i_EKmin,j;
  FLOAT yy,EK,EKmin;
  for (t=n_train-1-K-n_channels; t<n_train; t++) test->y[t]= givendata->y[t-k2];
  for (; t<n_total;t++) test->y[t]= 0;

  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  // 出力および誤差を計算
  for (t=n_train; t<n_total; t++) {
    {
      EKmin=1e30;
      i_EKmin=0;
      for(i=0;i<net->n_cells;i++){
	if(net->cell[i].v<1e-10) continue;
	EK=0;
	for(j=-K;j<0;j++){
	  tt=t+j;
	  for (k=0; k<k2;k++) test->x[tt][k] = test->y[tt+k2-k];
	  for (k=k2; k<n_channels; k++) test->x[tt][k] = test->y[tt-1-k];
	  test->x[tt][n_channels]=BIAS;//	  test->x[tt][n_channels]=1.0;
	  yy=0;
	  for(k=0;k<=n_channels;k++) yy += net->cell[i].am.M[0][k]*test->x[tt][k];
	  EK += square(test->y[tt]-yy);
	}
	if(i_EKmin==0 || EKmin>EK){
	  i_EKmin=i;
	  EKmin=EK;
	}
      }
    }
    printf("[%d]i_EKmin=%d,",t,i_EKmin);
    tt=t;
    for (k=0; k<k2;k++) test->x[tt][k] = test->y[tt+k2-k];
    for (k=k2; k<n_channels; k++) test->x[tt][k] = test->y[tt-1-k];
    test->x[tt][n_channels]=BIAS;//    test->x[tt][n_channels]=1.0;
    yy=0;
    for(k=0;k<=n_channels;k++) yy += net->cell[i_EKmin].am.M[0][k]*test->x[tt][k];
    test->y[t]=yy;
    test->Y[t]=moverange1(yy,net->ymin,net->ymax,net->ymin0,net->ymax0);
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += square(test->e[t]);

    MSE += square(test->e[t]);


    if (t == n_total-1) break;

    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = 1.0;

    for (k=0; k<n_channels; k++){
      if(k==k2) test->x[t+1][k] = test->y[t];
      else if(k<k2) test->x[t+1][k] = 0; //ARMA
      //else if(k<k2) test->x[t+1][k] = (t+1.)/5000.;//Jikantest
      //      else if(k<k2) test->x[t+1][k] = givendata->y[t+1+k2-k];//inchiki
      else test->x[t+1][k] = test->x[t][k-1];
    }
  }
  MSE /= n_test;
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  return;
}


void exec_msp_train(NET *net, DATA *givendata, DATA *test) {
    //全データの多段予測(multistep prediction)
  //char name[32] = "exec_msp_train";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  //  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;
  FLOAT e2;

  // 入力をコピー
  for (t=0; t<1; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  // 出力および誤差を計算
  for (t=0; t<n_train; t++) {
    //  for (t=1; t<n_total; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    e2=square(test->e[t]);
    net->cell[net->c].E += test->e[t];
    net->cell[net->c].S += e2;
    //    test->y[t]=yr;//?????
    MSE += e2;
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    test->x[t+1][0] = test->y[t];
    //    test->X[t+1][0] =  moverange1(test->x[t+1][0],net->xmin,net->xmax,net->xmin0,net->xmax0);
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
      //      test->X[t+1][k] = test->X[t][k-1];
    }
  }
  MSE /= n_train;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtrain;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  net->printf(strx2("PREDTRAIN: MSE=%e, NMSE=%e\n",MSE,NMSE));
  FILE *fp=fopen("mse0.dat","w");fprintf(fp,"%g %g",MSE,NMSE);fclose(fp);
  return;
}

void exec_msp_traintest(NET *net, DATA *givendata, DATA *test) {
    //全データの多段予測(multistep prediction)
  //char name[32] = "exec_msp_train";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  //  int n_train = givendata->n_train;
  //  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;

  // 入力をコピー
  for (t=0; t<1; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  // 出力および誤差を計算
  //  for (t=1; t<n_train; t++) {
  for (t=0; t<n_total; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    //    test->y[t]=yr;//?????
    MSE += square(test->e[t]);
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    test->x[t+1][0] = test->y[t];
    //    test->X[t+1][0] =  moverange1(test->x[t+1][0],net->xmin,net->xmax,net->xmin0,net->xmax0);
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
      //      test->X[t+1][k] = test->X[t][k-1];
    }
  }
  MSE /= n_total;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/(givendata->VARtrain+givendata->VARtest);
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  net->printf(strx2("PREDTRAIN: MSE=%e, NMSE=%e\n",MSE,NMSE));
  FILE *fp=fopen("mse0.dat","w");fprintf(fp,"%g %g",MSE,NMSE);fclose(fp);
  return;
}

void exec_msp_test1(NET *net, DATA *givendata, DATA *test) {
    //全データの多段予測(multistep prediction)
  //char name[32] = "exec_msp_train";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FLOAT y;

  // 入力をコピー
  for (t=0; t<1; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  // 出力および誤差を計算
  for (t=0; t<n_train; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    test->x[t+1][0] = test->y[t];
    //    test->X[t+1][0] =  moverange1(test->x[t+1][0],net->xmin,net->xmax,net->xmin0,net->xmax0);
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
      //      test->X[t+1][k] = test->X[t][k-1];
    }
  }
  //
  for (t=n_train; t<n_total; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    MSE += square(test->e[t]);
    if (t == n_total-1) break;

    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    test->x[t+1][0] = test->y[t];
    //    test->X[t+1][0] = moverange1(test->x[t+1][0],net->xmin,net->xmax,net->xmin0,net->xmax0);
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
      //      test->X[t+1][k] = test->X[t][k-1];
    }
  }

  MSE /= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  net->printf(strx2("PREDTRAIN: MSE=%e, NMSE=%e\n",MSE,NMSE));
  FILE *fp=fopen("mse0.dat","w");fprintf(fp,"%g %g",MSE,NMSE);fclose(fp);
  return;
}  //void exec_msp_test1(NET *net, DATA *givendata, DATA *test) {

//void exec_ssp_train_orig(NET *net, DATA *givendata, DATA *test) {
//  //訓練データの一段予測(single-step prediction)＝学習能力
//  //char name[32] = "exec_ssp_train";
//  FLOAT MSE = 0.0, NMSE;
//  //  FLOAT max = givendata->ymax, min = givendata->ymin;
//  int n_total = givendata->n_total;
//  int n_channels = givendata->k;
//  int n_train = givendata->n_train;
//  //  int n_test = givendata->n_test;
//  int t=0;
//  int k=0;
//  FILE *fp;
//  FLOAT y;
//
//  // 入力をコピー
//  for (t=0; t<n_train; t++) {
//    for (k=0; k<=(n_channels); k++) {
//      test->x[t][k] = givendata->x[t][k];
//      //      test->X[t][k] = givendata->X[t][k];
//    }
//  }
//  test->ymax = givendata->ymax;
//  test->ymin = givendata->ymin;
//
//  // 出力および誤差を計算
//  if((fp=open_file("predict0.dat","w"))==NULL) return;  //prediction of training data
//  for (t=0; t<n_train; t++) {
//    //original data
//    //    for(k=0;k<n_channels;k++) fprintf(fp,"%15.7e ",givendata->X[t][k]);
//    //    fprintf(fp,"%15.7e ",givendata->Y[t]);
//    //
//    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
//    test->e[t] = (test->Y[t]-givendata->Y[t]);
//    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t+1,givendata->Y[t],test->y[t]);
//
//    MSE += square(test->e[t]);
//    if (t == n_total-1) break;
//  }
//  fclose(fp);
//  MSE /= (n_train-1);
//  //  NMSE = MSE/(max-min);
//  NMSE = MSE/givendata->VARtrain;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//  test->MSE = MSE;
//  test->NMSE = NMSE;
//  printf("PREDTRAIN: MSE=%e, NMSE=%e\n",MSE,NMSE);
//  return;
//}
void exec_ssp_train(NET *net, DATA *givendata, DATA *test) {
  //訓練データの一段予測(single-step prediction)＝学習能力
  //char name[32] = "exec_ssp_train";
  FLOAT MSE = 0.0, NMSE;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  //  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  FILE *fp;
  FLOAT y;
  FLOAT e2;

  // 入力をコピー
  for (t=0; t<n_train; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;
  for(t=0;t<net->n_cells;t++) net->cell[t].E=net->cell[t].S=0;
  // 出力および誤差を計算
  if((fp=open_file("predict0.dat","w"))==NULL) return;  //prediction of training data
  for (t=0; t<n_train; t++) {
    //original data
    //    for(k=0;k<n_channels;k++) fprintf(fp,"%15.7e ",givendata->X[t][k]);
    //    fprintf(fp,"%15.7e ",givendata->Y[t]);
    //
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    fprintf(fp,"%15.7e %d %15.7e %15.7e #Y^ t Y yr\n",test->Y[t],t+1,givendata->Y[t],test->y[t]);
    net->cell[net->c].E += test->e[t];
    net->cell[net->c].S += e2;
    e2=square(test->e[t]);
    MSE += e2;
    if (t == n_total-1) break;
  }
  fclose(fp);
  MSE /= (n_train-1);
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtrain;
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  test->NMSE = NMSE;
  net->printf(strx2("PREDTRAIN: MSE=%e, NMSE=%e\n",MSE,NMSE));
  fp=fopen("mse0.dat","w");fprintf(fp,"%g %g",MSE,NMSE);fclose(fp);
  net_save(net,"bestssp.net");//for nips04
  return;
} //void exec_ssp_train(NET *net, DATA *givendata, DATA *test) {
/*====================================================================*
 * pred_out
 * 学習結果をチェックするための条件(時系列?)
 *====================================================================*/
void pred_out(NET *net, DATA *givendata, DATA *test) {
  //char name[32] = "pred_out";
  FLOAT MSE = 0.0, NMSE = 0.0;
  //  FLOAT max = givendata->ymax, min = givendata->ymin;
  int n_total = givendata->n_total;
  int n_channels = givendata->k;
  int n_train = givendata->n_train;
  int n_test = givendata->n_test;
  int t=0;
  int k=0;
  int n1;
  FILE *fp;
  char buff[256];
  FLOAT y;

  printf("File name to save:"); buff[0]=0;scanf1("%s",buff);
  if(buff[0]==0) return;
  printf("How many pred data:"); scanf1("%d",&n1);

  fp=open_file(buff, "w");
  // 入力をコピー
  for (t=0; t<n_total; t++) {
    for (k=0; k<=(n_channels); k++) {
      test->x[t][k] = givendata->x[t][k];
      //      test->X[t][k] = givendata->X[t][k];
    }
  }
  test->ymax = givendata->ymax;
  test->ymin = givendata->ymin;

  // 出力および誤差を計算
  for (t=0; t<n_train+givendata->t0;t++) fprintf(fp,"%e  %6d\n",givendata->y0[t],t);
  for (t=n_train; t<n_train+n1; t++) {
    calc_output(net, test->x[t],&(test->y[t]),&y,&(test->Y[t]));
    test->e[t] = (test->Y[t]-givendata->Y[t]);
    fprintf(fp,"%+e %+5d %+e\n",test->y[t],t+givendata->t0,givendata->Y[t]);

    MSE += square(test->e[t]);
    if (t == n_total-1) break;
    // 次の時刻での入力データを準備
    test->x[t+1][n_channels] = BIAS;//    test->x[t+1][n_channels] = 1.0;
    test->x[t+1][0] = test->y[t];
    //test->X[t+1][0]=moverange1(test->x[t+1][0],net->xmin,net->xmax,net->xmin0,net->xmax0);
    for (k=1; k<n_channels; k++){
      test->x[t+1][k] = test->x[t][k-1];
      //      test->X[t+1][k] = test->X[t][k-1];
    }
  }
  for (t=n_train+givendata->t0+n1; t<n_total+givendata->t0;t++) 
    fprintf(fp,"%e  %6d\n",givendata->y0[t],t);
  fclose(fp);
  MSE /= n_test;
  //  NMSE = MSE/(max-min);
  NMSE = MSE/givendata->VARtest;//???
  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
  test->MSE = MSE;
  //  test->NMSE = NMSE;

  return;
}


/*====================================================================*
 * cond_check_learn
 * 学習結果をチェックするための条件
 *====================================================================*/
int cond_check_learn (int i_times, int max_i_times, int d_times, int mode) {
  char name[32] = "cond_check_learn";
  int ret = 1;

  switch (mode) {
  case ONLINE_MODE: // Online Mode
    //    ret = ((i_times == 1 || (i_times < 10000 && i_times%100 == 0))
    //	   || (i_times%(max_i_times/d_times) == 0));
    //ret = (i_times%(max_i_times/d_times) == 0);
    ret = (i_times == 1 || i_times < 100
	   || (i_times%(max_i_times/d_times) == 0));
    break;
    
  case BATCH_MODE: // Batch Mode
    //    ret = (i_times < 100//)
    ret = (i_times == 1 || i_times < 100
	   || (i_times%(max_i_times/d_times) == 0));
    break;

  default:
    printf("error: unknown learning mode \"%d\" (%s)\n", mode, name);
    exit(-1);
  }
  return(ret);
}

/*====================================================================*
 * check_learn
 * 学習結果をチェックする
 *====================================================================*/
///*
////#define N_GP 5
//#define N_GP 6
//void check_learn(NET *net, DATA *givendata, DATA *test,MSEbank msebank, FILE **gp, int mode)
//{
//  //  int i_times=msebank.i_times[0];
//  int c_times=msebank.c_times[0];
//  FLOAT *MSEtrain=msebank.MSEtrain;
//  FLOAT *MSEssp=msebank.MSEssp;
//  FLOAT *MSEmsp=msebank.MSEmsp;
//  FLOAT *MSE1=msebank.MSE1;
//  FLOAT *MSE2=msebank.MSE2;
//  FLOAT *Entro=msebank.Entro;
//
//  int *MSEtime=msebank.MSEtime;
//
//  //void check_learn(NET *net, DATA *givendata, DATA *test,
//  //		 FLOAT *MSEssp, FLOAT *MSEmsp, int *MSEtime,
//  //		 FILE **gp, int i_times, int c_times, int mode,FLOAT *MSEtrain) {
//  //char name[32] = "check_learn";
//  int n_cells = net->n_cells;
//  FLOAT MIN, AVE, MAX;
//  int i;
//  FILE *gpi;
//  //  FLOAT yr;
//
//  // Calculate MSE
//  //printf("6:Time=%4d\n",GlobalTime);
//  //  exec_msp_test(net, givendata, test0);
//  //  MSEmsp[c_times] = test0->MSE;
//  if(givendata->data_class==TIME_SERIES){
//    exec_msp_test(net, givendata, test);
//    MSEmsp[c_times] = test->MSE;
//  }
////exec_ssp_train(net, givendata, test);MSEssp0[c_times] = test->MSE;
//  MSEtrain[c_times] = net->S[GlobalTime]/givendata->n_train;//Sは学習前のもの
//  exec_ssp_test(net, givendata, test);
//  MSEssp[c_times] = test->MSE;
//  MSE1[c_times] = test->MSE1;
//  MSE2[c_times] = test->MSE2;
//  Entro[c_times] =(1.-net->nentropy)*MSEtrain[c_times];//nentropyは学習前のもの
//  //  Entro[c_times-1] =-(1.-net->nentropy)*log(MSEtrain[c_times-1]);//nentropyは学習前のもの
//
//
//  //  Entro[c_times] =net->nentropy;
//  //  Entro[c_times] =(1-net->nentropy)*MSEtrain[c_times];
//  //  Entro[c_times] =-(1.-net->nentropy)*log(MSEtrain[c_times]);
//  //  MSEssp[c_times] = test->NMSE;
//  //  MSEmsp[c_times] = test0->NMSE;
//  MSEtime[c_times] = (int)GlobalTime;//    MSEtime[c_times] = i_times;//
//  //printf("7:Time=%4d\n",GlobalTime);
//
//  // Display
//  if (mode == 0) {
//    if(givendata->data_class==TIME_SERIES){
//      printf("MSE>%3dssp%8.2e(%8.2e)msp%8.2e(%8.2e)MSEtr%8.2eN%dk%dw%2.1fvm%d:%d\n",
//	     (int)GlobalTime, //	   (int)i_times, //
//	     MSEssp[c_times], MSEssp[c_times]/givendata->VARtest, 
//	     MSEmsp[c_times], MSEmsp[c_times]/givendata->VARtest,
//	     MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
//	     net->n_cells,net->k,net->width,net->vmin,net->vmin2);
//    }
//    else if(givendata->data_class==FUNCTION_APPROX){
////	printf("->%3dMSE012 %8.3e %8.2e %8.2eNMSE0%8.2eMSEtr%8.2eN%dk%dw%2.1fvm%d:%d\n",
////	       MSEtime[c_times],//c_times,//(int)GlobalTime, //	   (int)i_times, //
////	       MSEssp[c_times], test->MSE1,test->MSE2,
////	       MSEssp[c_times]/givendata->VARtest, 
////	       MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
////	       net->n_cells,net->k,net->width,net->vmin,net->vmin2);
//      printf("->%3dMSE%8.3eNMSE0%8.2eMSEtr%8.2eN%dk%dw%2.1fvm%d:%d\n",
//	     MSEtime[c_times],//c_times,//(int)GlobalTime, //	   (int)i_times, //
//	     MSEssp[c_times], 
//	     MSEssp[c_times]/givendata->VARtest, 
//	     MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
//	     net->n_cells,net->k,net->width,net->vmin,net->vmin2);
//    }
//    //    printf(">%3.0f MSE%e;%e NMSE%e;%e N%dk%dw%2.1f\n",
//    //	     (FLOAT)i_times, MSEssp[c_times], MSEmsp[c_times], 
//    //	     MSEssp[c_times]/givendata->VARtest, MSEmsp[c_times]/givendata->VARtest,
//    //	     net->n_cells,net->k,net->width);
//    // Mean Square Error
//    gpi=gp[0];// mse.vs.t
//    disp_mse(net,givendata, test, gp[0], "", "-", "", "", "",msebank);
//
//    if(net->l_mode==BATCH_MODE){
//      // Alpha for Every Unit
//      AVE = 0.0; MIN = MAX = net->cell[0].alpha;
//      for (i=0; i<n_cells; i++) {
//	AVE += net->cell[i].alpha;
//	if (net->cell[i].alpha < MIN) MIN = net->cell[i].alpha;
//	if (net->cell[i].alpha > MAX) MAX = net->cell[i].alpha;
//      }
//      AVE /= n_cells;
//      
//      gpi=gp[1];// alpha.vs.i
//      fprintf(gpi, "set title \"Alpha (MAX%3.2e,AVE%3.2e,MIN%3.2e)\"\n",
//	      MAX, AVE, MIN);
//      fprintf(gpi, "set data style boxes\n");
//      if(n_cells<=1){
//	fprintf(gpi, "###set yrange [%e:%e]\n", 0.,net->cell[0].alpha+0.1);
//	fprintf(gpi, "###set xrange [%e:%e]\n", 0.,1.);
//      }      
//      fprintf(gpi, "plot '-' using 1:2 t \"Alpha\"\n");
//      for (i=0; i<n_cells; i++) {
//	fprintf(gpi, "%d %e\n", i, net->cell[i].alpha);
//      }
//      fprintf(gpi, "e\n");
//      fflush(gpi);
//    }
//      //printf("1:Time=%4d,ctime=%d\n",GlobalTime,c_times);
//      //    printf("gp[%d]:Time=%4d\n",wn,GlobalTime);
//      
//    if(net->l_mode==BATCH_MODE){      // Square Error for Every Unit
//      AVE = 0.0; MIN = MAX = net->cell[0].S;
//      for (i=0; i<n_cells; i++) {
//	AVE += net->cell[i].S;
//	if (net->cell[i].S < MIN) MIN = net->cell[i].S;
//	if (net->cell[i].S > MAX) MAX = net->cell[i].S;
//      }
//      AVE /= n_cells;
//      gpi=gp[2];// S.vs.i
//      fprintf(gpi, "set title \"S(MAX%3.2e,AVE%3.2e,MIN%3.2e)\"\n",
//	      MAX, AVE, MIN);
//      fprintf(gpi, "set data style boxes\n");
//      if(n_cells<=1){
//	fprintf(gpi, "###set yrange [%e:%e]\n", 0.,net->cell[0].S);
//	fprintf(gpi, "###set xrange [%e:%e]\n", 0.,1.);
//      }      
//      fprintf(gpi, "plot '-' using 1:2 t \"S\"\n");
//      for (i=0; i<n_cells; i++) {
//	fprintf(gpi, "%d %e\n", i, net->cell[i].S);
//      }
//      fprintf(gpi, "e\n");
//      fflush(gpi);
//    }
//    if(net->l_mode==BATCH_MODE){      // Square Error for Every Unit
//      //    printf("2:Time=%4d,ctime=%d\n",GlobalTime,c_times);
//      // 
//      AVE = 0.0; MIN = MAX = net->cell[0].v;
//      for (i=0; i<n_cells; i++) {
//	AVE += net->cell[i].v;
//	if (net->cell[i].v < MIN) MIN = net->cell[i].v;
//	if (net->cell[i].v > MAX) MAX = net->cell[i].v;
//      }
//      AVE /= n_cells;
//      
//      gpi=gp[3];
//#define Sv
//#undef Sv
//#ifdef Sv
//      fprintf(gpi, "set title \"S vs. v(MAX%3.2e,AVE%3.2e,MIN%3.2e)\"\n",
//	      MAX, AVE, MIN);
//#else
//      fprintf(gpi, "set title \"Alpha vs. v(MAX%3.2e,AVE%3.2e,MIN%3.2e)\"\n",
//	      MAX, AVE, MIN);
//#endif
//      //    fprintf(gpi, "set data style boxes\n");
//      fprintf(gpi, "set data style points\n");
//      if(n_cells<=1){
//	fprintf(gpi, "###set yrange [%e:%e]\n", 0.,net->cell[0].v);
//	fprintf(gpi, "###set xrange [%e:%e]\n", 0.,1.);
//      }
//#ifdef Sv
//      fprintf(gpi, "plot '-' using 1:2 t \"S-v\"\n");
//#else
//      fprintf(gpi, "plot '-' using 1:2 t \"Alpha-v\"\n");
//#endif
//      for (i=0; i<n_cells; i++) {
//	//      fprintf(gpi, "%d %e\n", i, net->cell[i].v);
//#ifdef Sv
//	fprintf(gpi, "%e %e\n", net->cell[i].v, net->cell[i].S);
//#else
//	fprintf(gpi, "%e %e\n", net->cell[i].v, net->cell[i].alpha);
//#endif
//      }
//      fprintf(gpi, "e\n");
//      fflush(gpi);
//      //    printf("3:Time=%4d,ctime=%d\n",GlobalTime,c_times);
//      //    printf("gp[%d]:Time=%4d\n",wn,GlobalTime);
//    }
//    // kuro
//    
//    gpi=gp[4];//w1,w2
//    //    out_weights2d(net,gpi,gpi,"-");
//    //disp_weights(NET *net, FILE *fpgd0, char *gpl, char *dat, char *obj,char *title, char *type)
//    disp_weights(net, gp[4], "", "-", "","weights", "");
//    fflush(gpi);
//    //    printf("4:Time=%4d,ctime=%d\n",GlobalTime,c_times);
//    //kuro
//    gpi=gp[5];//err.vs.x
//    //    out_err2d(net, givendata, test, gpi, gpi,"-");
//    disp_yhat(net,       //NET *net,
//	      givendata, //DATA *givendata, 
//	      test,      //DATA *result, 
//	      gpi,       //FILE *fpgd0, 
//	      "",      //char *gpl,
//	      "-",       //char *dat, 
//	      0,         //char *obj,
//	      "yhat",    //char *title, 
//	      SSP,        //char *type
//	      0         //int disperr
//	      );
//  }
//  // Write Data File
//  else if (mode == 1) {//結果をファイルに出力
//    //    FILE *fp_dat = NULL, *fp_gpl = NULL;
//    char dat[64], gpl[64], obj[64];
//
//    strcpy(dat, "trans_mse.dat");
//    strcpy(gpl, "trans_mse.gpl");
//    strcpy(obj, "trans_mse.obj");
//
//    disp_mse(net,givendata, test, 0 , gpl, dat, obj, "", "",msebank);
//
//  }
//  else
//
//  return;
//}
//*/

char *check_learn(NET *net, DATA *givendata, DATA *testdata, MSEbank msebank)
//void check_learn(SIMDATA *simdata, NET *net, DATA *givendata, DATA *testdata, MSEbank msebank)
{
  //  int i_times=msebank.i_times[0];
  int c_times=msebank.c_times[0];
  FLOAT *MSEtrain=msebank.MSEtrain;
  FLOAT *MSEssp=msebank.MSEssp;
  FLOAT *MSEmsp=msebank.MSEmsp;
  FLOAT *MSE1=msebank.MSE1;
  FLOAT *MSE2=msebank.MSE2;
  FLOAT *MSEdp=msebank.MSEdp;//ueno(Tue Feb 17 00:22:49 2004)
  FLOAT *MSEdl=msebank.MSEdl;
  FLOAT *MSEds=msebank.MSEds;
  FLOAT *MSEmean=msebank.MSEmean;
  FLOAT *Entro=msebank.Entro;
  int *MSEtime=msebank.MSEtime;
  char *mes="???";
  //calcuration
  {
    //calc train
    //#orig
    MSEtrain[c_times] = net->S[GlobalTime]/givendata->n_train;//
    Entro[c_times] =(1.-net->nentropy)*MSEtrain[c_times];//
    //only modify MSEtrain
    //    MSEtrain[c_times]=givendata->MSEtr;//20191012 check whether is S=MSEtr good or not

    //    Entro[c_times] =(1.-net->nentropy)*(MSEtrain[c_times]-net->sigma2_hat);//
    //calc ssp
    {
      exec_ssp_test(net, givendata, testdata);
      MSEssp[c_times] = testdata->MSE;
      MSE1[c_times] = testdata->MSE1;
      MSE2[c_times] = testdata->MSE2;
      MSEtime[c_times] = (int)GlobalTime;
    }
    if(givendata->data_class==TIME_SERIES){//after ssp for display Y
      if(application==IJCNN04_APPLI) {
	exec_msp_test_IJCNN04(net, givendata, testdata);//test
	MSEdp[c_times] = givendata->ijcnn04data->MSEdp;//ueno(Mon Feb 16 23:14:03 2004)
	MSEdl[c_times] = givendata->ijcnn04data->MSEdl;
	MSEds[c_times] = givendata->ijcnn04data->MSEds;
	MSEmean[c_times] = givendata->ijcnn04data->MSEmean;
      }
      else exec_msp_test(net, givendata, testdata,0,0);//test
      //      exec_msp_test_can2_2(net, givendata, testdata);//test
      MSE1[c_times]=MSEmsp[c_times] = testdata->MSE;
      MSEtime[c_times] = (int)GlobalTime;
    }
    if(givendata->data_class==TIME_SERIES){
      if(application==IJCNN04_APPLI) {
	printf("MSE>%3dMSEdl%8.2eMSEds%8.2eMSEdp%8.2eN%dk%d+%dw%2.1fvm%d:%d\n",
	       (int)GlobalTime,
	       MSEdl[c_times],MSEds[c_times],MSEdp[c_times],
	       net->n_cells,net->k1,net->k2,net->width,net->vmin,net->vmin2);
      }
      else{
	printf("MSE>%3dssp%8.2e(%8.2e)msp%8.2e(%8.2e)MSEtr%8.2eN%dk%d+%dw%2.1fvm%d:%d\n",
	       (int)GlobalTime, //	   (int)i_times, //
	       MSEssp[c_times], MSEssp[c_times]/givendata->VARtest, 
	       MSEmsp[c_times], MSEmsp[c_times]/givendata->VARtest,
	       MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
	       net->n_cells,net->k1,net->k2,net->width,net->vmin,net->vmin2);
      }
    }
    else if(givendata->data_class==FUNCTION_APPROX){
      mes=strx10("%d %8.3e %8.3e %8.2e %8.2e MSEtr MSE NMSEtr NMSE N%dk%dw%2.1fvm%d:%d\n",
                 MSEtime[c_times],//c_times,//(int)GlobalTime, //	   (int)i_times, //
                 MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
                 MSEssp[c_times], 
                 MSEtrain[c_times]/givendata->VARtrain, 
                 MSEssp[c_times]/givendata->VARtest, 
                 net->n_cells,net->k,net->width,net->vmin,net->vmin2);

//      printf("->%3dMSE%8.3eNMSE0%8.2eMSEtr%8.2eN%dk%dw%2.1fvm%d:%d\n",
//	     MSEtime[c_times],//c_times,//(int)GlobalTime, //	   (int)i_times, //
//	     MSEssp[c_times], 
//	     MSEssp[c_times]/givendata->VARtest, 
//	     MSEtrain[c_times],  //	     MSEssp0[c_times]/givendata->VARtrain,
//	     net->n_cells,net->k,net->width,net->vmin,net->vmin2);
    }
    else {
    }
  }

  //display
//#if GPLT >= 1
//  {
//    //    simdata_net_display(simdata, net, GlobalTime);
//    //    if (GlobalTime>=1) simdata_ability_display(simdata, &msebank, GlobalTime);
//    if (givendata->data_class == TIME_SERIES) {
//      //      simdata_error_msp_display(simdata, givendata, testdata, GlobalTime);
//      //      simdata_func_msp_display(simdata, givendata, testdata, GlobalTime);
//    }
//
//    else if (givendata->data_class == FUNCTION_APPROX) {
//      //      simdata_error_ssp_display(simdata, givendata, testdata, GlobalTime);
//      //      simdata_func_ssp_display(simdata, givendata, testdata, GlobalTime);
//    }
//    else {
//    }
//  }
//#endif
//
//  //write out
//#if GPLT >= 1
//  {
//    if (GlobalTime == msebank.max_i_times) {
//      //      simdata_net_write(simdata, net, GlobalTime);
//      //      if (GlobalTime>=1) simdata_ability_write(simdata, &msebank, GlobalTime);
//      if (givendata->data_class == TIME_SERIES) {
////	simdata_error_msp_write(simdata, givendata, testdata, GlobalTime);
////	simdata_func_msp_write(simdata, givendata, testdata, GlobalTime);
//      }
//      else if (givendata->data_class == FUNCTION_APPROX) {
//	simdata_error_ssp_write(simdata, givendata, testdata, GlobalTime);
//	simdata_func_ssp_write(simdata, givendata, testdata, GlobalTime);
//      }
//      else {
//      }
//    }
//  }
//#endif
  return(mes);
}

/*====================================================================*
 * learn_net_base
 * ネットに学習を行う
 *====================================================================*/

void learn_net_base(NET *net, DATA *givendata, DATA *test, int mode) {
  char name[32] = "learn_net_base";
  //  BDATA *bdata = NULL;
  //FILE *gp[N_GP];
  //  FLOAT **w_cache = NULL;
  //  FLOAT *MSEgivendata0 = NULL, *MSEssp = NULL, *MSEmsp = NULL;
  //  FLOAT *MSEgivendata0 = NULL, *MSEssp = NULL, *MSEmsp = NULL;
  FLOAT *MSEssp = NULL, *MSEmsp = NULL,*MSEtrain;
  FLOAT *MSE1, *MSE2;
  FLOAT *MSEdp = NULL, *MSEdl = NULL, *MSEds = NULL;//ueno(Mon Feb 16 23:30:23 2004)
  FLOAT *MSEmean = NULL;//ueno(Wed Feb 25 05:08:10 2004)
  int *MSEtime = NULL;
  int max_i_times = net->i_times;
  int i_times = 0;
  int c_times = 0;
  int d_times = net->d_times;
  //  int n_cells = net->n_cells;
  //  int n_channels = net->k;
  int n_train = givendata->n_train;
  int t;//  int i, t;
  int is=0,im=0,is0=0,ine=0;
  //int lbuff=(102+max_i_times/d_times);//100以下はすべて表示
  int lbuff=(int)(max_i_times+10);//ueno(Fri Feb 13 00:56:01 2004)
  FLOAT *Entro;
  MSEbank msebank;
  char *mes="OK?";
  //  SIMDATA *simdata = NULL;//ueno
  //  int lbuff=(1+max_i_times/d_times);
  // 学習回数毎のMSEを記録
  msebank.MSEssp = MSEssp = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEmsp = MSEmsp = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEtrain= MSEtrain= (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEtime = MSEtime = (int *)malloc(sizeof(int)*lbuff);
  msebank.MSE1 = MSE1 = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSE2 = MSE2 = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEdp = MSEdp = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);//ueno(Mon Feb 16 23:29:46 2004)
  msebank.MSEdl = MSEdl = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEds = MSEds = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.MSEmean = MSEmean = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.Entro = Entro = (FLOAT *)malloc(sizeof(FLOAT)*lbuff);
  msebank.VARtrain = &(givendata->VARtrain);//ueno
  msebank.VARssp = &(givendata->VARtest);
  msebank.VARmsp = &(givendata->VARtest);
  msebank.max_i_times = max_i_times;//ueno
  msebank.c_times =&c_times;
  msebank.i_times =&i_times;
//  if (MSEgivendata0 == NULL) {
//    printf("error: memory allocation error #?. (%s)\n", name);
//    exit(-1);
//  }
  if (MSEssp == NULL || MSEmsp == NULL || MSEtime == NULL) {
    printf("error: memory allocation error #?. (%s)\n", name);
    exit(-1);
  }
  //  for (t=0; t<(10000+max_i_times/d_times); t++) MSEgivendata0[t] = 0.0;
  for (t=0; t<lbuff; t++){
    MSEssp[t] = MSEmsp[t] = MSEtime[t] = 0.0;
    MSE1[t] = MSE2[t] = 0.0;
    MSEdp[t] = MSEdl[t] = MSEds[t] = 0.0;//ueno(Mon Feb 16 23:32:45 2004)
    MSEmean[t] = 0.0;//ueno(Wed Feb 25 05:08:30 2004)
  }

  // GNUPLOTへのプロセスを開く
  //for (i=0; i<N_GP; i++) gp[i] = open_pipe(GNUPLOT, "w");
  //  FILE *fpmse=fopen("mse_t.dat","w");
  switch (mode) {
  case ONLINE_MODE: // Online Mode
    net->printf("> Learning Mode :  ONLINE MODE\n\n");
    //    ReinitTime=-100;
    //    GlobalTime=0;
    //    GlobalTimeMax += max_i_times;
    //    if(GlobalTime==0) simdata = simdata_create("online", net, givendata, givendata->data_class);//ueno
//#if GPLT >= 1  //070130
//      simdata = simdata_create("online", net, givendata, givendata->data_class, givendata->application);//ueno(Tue Feb 17 00:20:13 2004)
//#endif
    mytimer_start();
    while (i_times < max_i_times) {
      ++GlobalTime;
      ++i_times;
      for (t=0; t<n_train; t++) {
	//	printf("x[%d]=%+e %+e y=%+e\n",t,givendata->x[t][0],givendata->x[t][1],givendata->y[t]);//check
	store_vector(net, givendata->x[t], givendata->y[t]);
      }
      // 学習状況のチェック
      if (cond_check_learn(i_times, max_i_times, d_times, mode)) {
	FLOAT yr,y,Y,S=0;
	int tt;
	for(tt=0;tt<n_train;tt++){
	  calc_output(net, givendata->x[tt],&yr,&y,&Y);
	  S += square(Y-givendata->Y[tt]);
	}
	net->S[GlobalTime]=S;
	//check_learn(net, givendata, test, msebank,gp, 0);
        //	check_learn(simdata, net, givendata, test, msebank);//ueno
	if(MSEssp[is]>MSEssp[c_times]) is=c_times;
	if(givendata->data_class==TIME_SERIES && MSEmsp[im]>MSEmsp[c_times]) im=c_times;
	if(MSEtrain[is0]>MSEtrain[c_times]) is0=c_times;
	++c_times;
      }
      //	if(c_times>5) break;
      //	if (i_times == max_i_times)
      //	  check_learn(net, givendata, test,
      //		      MSEssp, MSEmsp, MSEtime,
      //		      gp, i_times, c_times, 1);
      //      if(i_times >= max_i_times) break;
      //      }
      //      check_learn(net, givendata, test, msebank,gp, 1);
      //      printf("Min MSEssp[%d]=%e,MSEmsp[%d]=%e.MSEtrain[%d]=%e",
      //   MSEtime[is],MSEssp[is],MSEtime[im],MSEmsp[im],MSEtime[is0],MSEtrain[is0]);
    }
    break;

  case BATCH_MODE: // Batch Mode

    //    mytimer_start();//fprintf(stderr,"#laptime:start=%5.3f\n",mytimer_lap());//for check

    net->printf("> Learning Mode :  BATCH MODE\n\n");
//#if GPLT >= 1
//    //    if(GlobalTime==0) simdata = simdata_create("batch", net, givendata, givendata->data_class);//ueno
//    simdata = simdata_create("batch", net, givendata, givendata->data_class, givendata->application);//ueno(Tue Feb 17 00:20:27 2004)
//#endif
    net->S[0] = 1e+12;
    //    w_cache = init_batch_wvector(givendata->x, givendata->y, n_cells, n_channels, n_train);
    //    set_w_batch(net, w_cache);
    init_net_batch(net, givendata->x, n_train);
    GlobalTime=0;//kuro050206
    //    ReinitTime=-100;
    //    ReinitTime=0;
    //    GlobalTimeMax += max_i_times;
    mytimer_start();
    for(;;) {    //    while (i_times < max_i_times) {
      ++i_times;
      ++GlobalTime;
      net->GlobalTime=GlobalTime;
      store_vector_batch(net, givendata->x, givendata->y, n_train,0);
      //      if(i_times>=max_i_times) break;//??
      // 学習状況のチェック
      if (cond_check_learn(i_times, max_i_times, d_times, mode)) {
	//check_learn(net, givendata, test, msebank,gp, 0);
	mes=check_learn(net, givendata, test, msebank);//ueno
        //	check_learn(simdata, net, givendata, test, msebank);//ueno
	if(c_times==0 || MSEssp[is]>MSEssp[c_times]){
	  is=c_times;
	  net_save(net,"bestssp.net");//for nips04
	}
	if(givendata->data_class==TIME_SERIES && (c_times==0 || MSEmsp[im]>MSEmsp[c_times])){
	  im=c_times;
	  net_save(net,"bestmsp.net");
	}
	if(c_times==0 || MSEtrain[is0]>MSEtrain[c_times]){
	  is0=c_times;
	  //	  net_save(net,"besttrain.net");
	}
	if(c_times==0 || Entro[ine]>Entro[c_times]){//
	  ine=c_times;
	}
	++c_times;
      }
      if(i_times>=max_i_times) break;
      store_vector_batch(net, givendata->x, givendata->y, n_train,1);
    }
    break;
    
  default: // Oops ,,,
    printf("Error: Undefined Learning Mode \"%d\". (%s)\n", mode, name);
    exit(-1);
  }
  double elapsedtime=mytimer_total();
  {//      if (i_times == max_i_times){
    //check_learn(net, givendata, test, msebank,gp, 1);//??lasttime two times
    char *mes=check_learn(net, givendata, test, msebank);//ueno
//    printf("1\n");
    net->printf(strx1("%s",mes));
    //    net->printf("%s %p",mes,&mes);
//    printf("2\n");
    //    check_learn(simdata, net, givendata, test, msebank);//ueno
    net_save(net,"last.net");
    FILE *fp=fopen("tmp/mse_t.csv","w");
    int i;
    for(i=0;i<max_i_times;i++){
      fprintf(fp,"%.7e %.7e %.7e %.7e\n",MSEtrain[i],MSEssp[i],MSEtrain[i]/givendata->VARtrain,MSEssp[i]/givendata->VARtrain);
      //fprintf(stdout,"%.7e %.7e %.7e %.7e\n",MSEssp[i],MSEtrain[i],MSEssp[i]/givendata->n_train,MSEtrain[i]/givendata->n_train);
    }
    fclose(fp);

    //
    net->printf("last!!!!!!!!!!!\n");
  }
  //    remove_batch_wvector(w_cache, n_cells);
  net->printf(strx1("total time elapsed=%5.3f\n",mytimer_total()));
  //    check_learn(net, givendata, test,
  //		  MSEssp, MSEmsp, MSEtime,
  //		  gp, i_times, c_times-1, 1);
  net->printf("*********\n");
  if(givendata->data_class==TIME_SERIES){
    net->printf(">Minimum MSE(NMSE) of ssp and msp:\n");
    net->printf(strx3(">>%dssp%10.4e(%8.6f)",
	   MSEtime[is],
	   MSEssp[is],
	   MSEssp[is]/givendata->VARtest
	   // Entro[is],
	   //MSEtime[ine],
	   //MSEssp[ine],Entro[ine]
                      ));
    net->printf(strx3("%dmsp%10.4e(%8.6f)",
                      MSEtime[im],
                      MSEmsp[im],
                      MSEmsp[im]/givendata->VARtest));
    net->printf(strx3("%dtrain%10.4e(%8.6f)",
                      MSEtime[is0],
                      MSEtrain[is0],
                      MSEtrain[is0]));
    net->printf(strx5("N%dk%dw%2.1fvm%d:%d\n",net->n_cells,net->k,net->width,net->vmin,net->vmin2));
    int t;
    net->tpH=givendata->n_test;
    for(t=givendata->n_train;t<givendata->n_total;t++){
      if(fabs(givendata->Y[t]-test->Y[t])>net->tpEy){
	net->tpH=t-givendata->n_train;
	break;
      }
    }
    int i=max_i_times;
    sprintf(&net->msg[strlen(net->msg)],"%d(%.3fs) %.3e %.3e #ep(time),MSEtr,MSE k%d N%d T%d,%d seed%lu nop%d H%.0f(Ey%g) ",i,elapsedtime,MSEtrain[i],MSEssp[i],net->k,net->n_cells,net->i_times,net->Tpinv,net->seed,net->nop,net->tpH,net->tpEy);
    //    sprintf(&net->msg[strlen(net->msg)],"%d(%.3fs) %.3e %.3e %.3e #ep(time),MSEtr,MSEtst,MSEmsp seed%d nop%d H%.0f(Ey%g) ",i,elapsedtime,MSEtrain[i],MSEssp[i],MSEmsp[i],net->seed,net->nop,net->tpH,net->tpEy);
  }
  else if(givendata->data_class==FUNCTION_APPROX){
//    printf("\n>Minimum MSE(NMSE) of ssp:\n");
//    printf(">>%dssp%10.4eH%e>>%dssp%10.4eH%e",
//	     MSEtime[is],
//	     MSEssp[is],//	   MSEssp[is]/givendata->VARtest,
//	     Entro[is],
//	     MSEtime[ine],MSEssp[ine],Entro[ine]);
//    printf("N%dk%dw%2.1fvm%d:%d\n",net->n_cells,net->k,net->width,net->vmin,net->vmin2);

    {
      FILE *fp;
      fp=fopen("is.dat","w+");
      fprintf(fp,"%d %15.7e %d %15.7e %d %d %d %15.7e %15.7e\n",
	      //MSEtime[is],MSEssp[is],givendata->n_train,MSEtrain[c_times],net->n_cells2,net->n_cells,givendata->n_test);
	      MSEtime[is],MSEssp[is],givendata->n_train,MSEtrain[c_times],net->n_cells2,net->n_cells,givendata->n_test,MSEssp[c_times],MSEtrain[is]);
      //      printf("?????????????????????\n");
      fclose(fp);
    }
    net->printf("\n>Minimum MSE(NMSE) of ssp:\n");
    net->printf(strx4(">>%dtest%10.4e,%dtrain%10.4e",
                      MSEtime[is],
                      MSEssp[is],//	   MSEssp[is]/givendata->VARtest,
                      MSEtime[is0],
                      MSEtrain[is0]));
    net->printf(strx5("N%dk%dw%2.1fvm%d:%d\n",net->n_cells,net->k,net->width,net->vmin,net->vmin2));
    int i=max_i_times;
    sprintf(&net->msg[strlen(net->msg)],"%d(%.3fs) %.3e %.3e #ep(time),MSEtr,MSEtst ",i,elapsedtime,MSEtrain[i],MSEssp[i]);
  }
  //
  net->printf("To see the best msp or ssp result, do as follows.\n");
  net->printf("nl                          \n");
  net->printf("bestmsp.net       (or bestssp.net)\n");
  net->printf("msp               (or ssp) \n");
  net->printf("*********\n");
  
  free(MSEssp);
  free(MSEtime);
  free(MSEmsp);
  free(MSEtrain);
  free(MSE1);//ueno(Fri Feb 13 00:57:24 2004)
  free(MSE2);
  free(MSEdp);//ueno(Mon Feb 16 23:33:01 2004)
  free(MSEdl);
  free(MSEds);
  free(MSEmean);
  free(Entro);
//#if GPLT >= 1
//  simdata_remove(simdata);
//#endif
  return;
}

/*====================================================================*
 * learn_net
 * ネットの学習を行う
 *====================================================================*/
void learn_net(NET *net, DATA *givendata, DATA *test) {
  char name[32] = "learn_net";
  int mode = net->l_mode;

  net->printf("================================================================\n");
  net->printf("	Learning\n");
  net->printf("================================================================\n");

  switch (mode) {
  case ONLINE_MODE: // オンライン学習
    learn_net_base(net, givendata, test, mode);
    break;

  case BATCH_MODE: // バッチ学習
    learn_net_base(net, givendata, test, mode);
    break;

  default: // エラー処理
    printf("Error: Undefine Learning Mode \"%d\". (%s)\n", mode, name);
    exit(-1);
  }
  return;
}
//void svr_train1(NET *net, DATA *givendata, DATA *test)
//{
//  FILE *fp2;
//  FLOAT r,C,epsilon;
//  int dim;
//  char buff[256];
//  FLOAT max_iterations;
//  int t,k,i;
//
//  printf("* r, epsilon, max_iterations : "); scanf3("%lf%lf%lf", &r,&epsilon,&max_iterations);//
//  C=1000;
//  //  dim=givendata->k-1;
//  sprintf(buff,"%s","svrparam.dat");
//  if((fp2=fopen(buff,"w+"))==NULL){
//    printf("File(%s) Open Eoor\n",buff);
//    fclose(fp2);
//    return;
//  }
//  fprintf(fp2,"@kernel\n");
//  fprintf(fp2,"type radial gamma %f\n",1./r/r);
//  fprintf(fp2,"@parameters\n");
//  fprintf(fp2,"C %f\n",C);
//  fprintf(fp2,"epsilon %f\n",epsilon);
//  fprintf(fp2,"max_iterations %10.0f\n",max_iterations);
//  fprintf(fp2,"regression\n");
//  fclose(fp2);
//  
//  sprintf(buff,"%s","svrtrain.dat");
//  if((fp2=fopen(buff,"w+"))==NULL){
//    printf("File(%s) Open Eoor\n",buff);
//    fclose(fp2);
//    return;
//  }
//  fprintf(fp2,"@examples\n");
//  fprintf(fp2,"format xy\n");
//  fprintf(fp2,"dim %d\n",givendata->k);
//
//  for (t=0; t<givendata->n_train; t++) {
//    for(k=0;k<givendata->k;k++) fprintf(fp2,"%e ", givendata->x[t][k]);
//    fprintf(fp2,"%e\n", givendata->y[t]);
//  }
//  fclose(fp2);
//  //
//  sprintf(buff,"%s","svrtest.dat");
//  if((fp2=fopen(buff,"w+"))==NULL){
//    printf("File(%s) Open Eoor\n",buff);
//    fclose(fp2);
//    return;
//  }
//  //
//  fprintf(fp2,"@examples\n");
//  fprintf(fp2,"format xy\n");
//  fprintf(fp2,"dim %d\n",dim);
//  //  learningsteps=n_train;
//  //learningsteps=2200;
//  for(t=givendata->n_train+1;t<givendata->n_total;t++){
//    for(k=0;k<givendata->k;k++) fprintf(fp2,"%e ", givendata->x[t][k]);
//    fprintf(fp2,"%e\n", givendata->y[t]);
//  }
//  fclose(fp2);
//  system("time ~/bin/mysvm svrparam.dat svrtrain.dat");
//  system("~/bin/predict svrparam.dat svrtrain.dat.svm svrtest.dat");
//  system("./pred2y svrtest.dat.pred svrtest.dat>y.dat");
//#if GPLT >= 1
//  system("xterm -geometry 80x5-0+0 -e gnuplot svry.gp&");
//#endif
//  printf("%5.0f) dim%d C%3.0f r%5.3f epsilon%6.4f###\n",max_iterations,dim,C,r,epsilon);
//}
//void svr_msp_test(NET *net, DATA *givendata, DATA *test)
//{//テストデータの多段予測
//  FILE *fp2,*fp;
//  char buff[256];
//  int t,k;
//  FLOAT MSE=0,NMSE;
//
//  sprintf(buff,"%s","y1.dat");
//  if((fp=fopen(buff,"w+"))==NULL){
//    printf("File Open Eoor\n");
//    fclose(fp);
//    return;
//  }
//
//  // 入力をコピー
//  for (t=givendata->n_train; t<=givendata->n_train; t++) {
//    for (k=0; k<=(givendata->k); k++) {
//      test->X[t][k] = givendata->X[t][k];
//      //      test->X[t][k] = givendata->X[t][k];
//    }
//  }
//  for(t=givendata->n_train+1;t<givendata->n_total;t++){
//    sprintf(buff,"%s","svrpred.dat");
//    if((fp2=fopen(buff,"w+"))==NULL){
//      printf("File(%s) Open Eoor\n",buff);
//      fclose(fp2);
//      return;
//    }
//    fprintf(fp2,"@examples\n");
//    fprintf(fp2,"format xy\n");
//    fprintf(fp2,"dim %d\n",givendata->k);
//    for(k=0;k<givendata->k;k++){
//      fprintf(fp2,"%e ", test->X[t][k]);
//    }
//    fprintf(fp2,"%e\n", givendata->Y[t]);
//    fclose(fp2);
//    system("predict svrparam.dat svrtrain.dat.svm svrpred.dat>/dev/null");
//    sprintf(buff,"%s","svrpred.dat.pred");
//    if((fp2=fopen(buff,"r+"))==NULL){
//      printf("File(%s) Open Eoor\n",buff);
//      fclose(fp2);
//      return;
//    }
//    fgets(buff,256,fp2);//read @examples
//    fgets(buff,256,fp2);//read the data
//    sscanf(buff,"%lf",&(test->Y[t]));
//    fclose(fp2);
//
//    test->e[t]= test->Y[t]-givendata->Y[t];
//    MSE += square(test->e[t]);
//    if (t == givendata->n_total-1) break;
//
//    test->X[t+1][0]=test->Y[t];
//    for (k=1; k<givendata->k; k++){
//      test->X[t+1][k] = test->X[t][k-1];
//    }
//    //    fprintf(fp,"%d %e %e %e\n",i,test->Y[t],test->y[t],test->e[t]);
//  }
//  fclose(fp);
//  MSE/=givendata->n_test;
//  NMSE = MSE/givendata->VARtest;
//  if(MSE>Infty) MSE=Infty; if(NMSE>Infty) NMSE=Infty;
//  test->MSE = MSE;
//  test->NMSE = NMSE;
//
//  printf("\nSVR PREDICTION: MSE=%e NMSE=%e\n", MSE,NMSE);
//}
/*====================================================================*
 * exec_sim
 * シミュレーションを実行する
 *====================================================================*/
void exec_sim(NET *net, DATA *givendata, DATA *test) {
  //  FLOAT i_times,d_times;
  int t;
  //  if(net->i_times>0) SuccessiveLearn=1;
  //  else SuccessiveLearn=0;

  //  printf("* Learning Mode?             : "); scanf1("%d", &(net->l_mode));//
  //  printf("* Learning Mode, Output Precision(1=int)?:"); scanf3("%d%d%d", &(net->l_mode),&net->r1,&net->r2);//
  //  printf("resolution=%f=%d/%d\n",(FLOAT)net->r1/net->r2,net->r1,net->r2);
  //  printf("* Learning Mode (0:online,1:batch)?:"); scanf1("%d", &(net->l_mode));//
  net->printf("* Learning Mode (0:online,1:batch), gamma0, entropy_thresh?:"); 
  scanf3("%d%lf%lf",&(net->l_mode),&(net->gamma0),&(net->nentropy_thresh));//
  net->printf("* Number of Total Iterations?           : "); scanf1("%d", &(net->i_times));//
  //  printf("* Number of Displaying Interval Iterations? : "); scanf1("%d", &(net->d_times));//
  //  printf("* Number of Displaying Interval Iterations? : "); scanf1("%d", &(net->d_times));//
  net->printf("* Displaying Interval, rot_x, rot_z ? : "); scanf3("%d%d%d", &(net->d_times),&(givendata->rot_x),&(givendata->rot_z));//
  //  printf("* gamma0, nentropy_thresh? : "); scanf2("%lf%lf", &(net->gamma0),&(net->nentropy_thresh);//

  //  GlobalTimeMax = net->i_times;  //  GlobalTimeMax += net->i_times;
  GlobalTimeMax += net->i_times; 
//  net->i_times= net->i_times - GlobalTime; 
//  if(net->d_times<=1) net->d_times=1;
//  if(net->d_times>=net->i_times) net->d_times=net->i_times;
  net->printf(strx2("GlobalTime=%d & GlobalTimeMax=%d:",(int)GlobalTime,(int)GlobalTimeMax));

  //  net->S = (FLOAT *)malloc(sizeof(FLOAT)*((int)(net->i_times+1)));
  //  for (t=0; t<(int)(net->i_times+1); t++) net->S[t] = 0.0;
  if(net->S!=NULL) free(net->S);
  net->S = (FLOAT *)malloc(sizeof(FLOAT)*((int)(GlobalTimeMax+2)));
  for (t=0; t<=(int)(GlobalTimeMax); t++) net->S[t] = 0.0;

//  printf("RAND=%d seed=%lu\n",RAND,net->seed);
//#if RAND == RANDOM
//  srandom((unsigned int)net->seed);
//#elif RAND == DRAND48
//  seed48((unsigned short*)&net->seed);
//#elif RAND == MYRAND
//  _randn=(unsigned long)net->seed;
//#elif RAND == ZMTRAND
//  InitMt((unsigned long)net->seed);
//#endif

//
//#ifdef MYRANDOM
//  printf("MYRANDOM seed=%lu\n",net->seed);
//  _randn=net->seed;
//#else
//  printf("noMYRANDOM seed=%lu\n",net->seed);
//  srandom(net->seed);
//#endif
  // ネットの学習を行う
  learn_net(net, givendata, test);//ここでtrans_mse.gplをセーブ

  // 一段予測を行う
  if(1==0){
    exec_ssp_test(net, givendata, test);//check_learnでもやっている?
    exec_plot(givendata, test, "Sinle-Step Prediction",SSP,net);//ここでy_ssp.gplなどをセーブ
    
    if(givendata->data_class==TIME_SERIES){
      // 多段予測を行う
      exec_msp_test(net, givendata, test,0,0);
      exec_plot(givendata, test, "Multi-Step Prediction",MSP,net);
    }
  }
  free(net->S);
  SuccessiveLearn=1;
  return;
}


/* this file ends here,,, */
