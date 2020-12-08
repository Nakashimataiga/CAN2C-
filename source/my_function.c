/*
 *  $id: time-stamp: "2004/02/17 00:48:12" / takamasa ueno exp.$040213kuro
 *
 * my_function.c
 *
 */
#include "my_function.h"
#define GNUPLOT "gnuplot -geometry 400x300"
//#define DEBUG

///*====================================================================*
// * get_function_id
// * 予測対象となる関数を選ぶ
// *
// * data_class   : function id
// * fname : file name （時系列データがファイルにある場合）
// *====================================================================*/
//void get_function_id(int *data_class, char *fname) {
//  printf("================================================================\n");
//  printf("	Select The Objective Function\n");
//  printf("================================================================\n");
//  printf("	<funcition id list>\n");
//  printf("	他の関数はまだ移植してません (^^;\n");
//  printf("	%d : mackey-glass chaos time series\n", MACKEY_GLASS);
//  printf("----------------------------------------------------------------\n");
//  printf("* function id? : "); scanf1("%d", data_class);
//
//  switch (*data_class) {
//  case MACKEY_GLASS: // Makcey-Glass 時系列
//    printf("* data file?   : "); fname[0]=0;scanf1("%s", fname);
//    break;
//
//  default:
//  }
//
//  return;
//}
//
///*====================================================================
// *
// * TDATA（時系列データ構造体）の処理
// *
// *====================================================================*/
//
///*====================================================================*
// * create_time_data
// * 時系列データを初期化
// *
// * t   : total steps (time channels)
// * data_class : function id
// *====================================================================*/
//TDATA* create_time_data(int t, int data_class) {
//  char name[32] = "create_time_data";
//  TDATA *data = NULL;
//
//  data = (TDATA *)malloc(sizeof(TDATA));
//  if (data == NULL) {
//    printf("error: memory allocation error #1. (%s)\n", name);
//    exit(-1);
//  }
//  data->t = t;
//  data->data_class = data_class;
//  data->x = (FLOAT *)malloc(sizeof(FLOAT)*t);
//  if (data->x == NULL) {
//    printf("error: memory allocation error #2. (%s)\n", name);
//    exit(-1);
//  }
//  return(data);
//}
//
///*====================================================================*
// * remove_time_data
// * 時系列データを破棄
// *====================================================================*/
//void remove_time_data(TDATA *tdata) {
//  free(tdata->x);
//  free(tdata);
//  return;
//}
//
///*====================================================================*
// * init_time_data
// * 時系列データを読み込む
// *
// * data_class   : function id
// * fname : file name
// *====================================================================*/
//TDATA* init_time_data(int data_class, char *fname) {
//  char name[32] = "init_time_data";
//  TDATA *data;
//  FILE *fp;
//  char line[512];
//  int size;
//  int t, i;
//
//  // ファイルを開き行数をカウント
//  fp = open_file(fname, "r");
//  size = count_file(fp);
//  close_file(fp);
//
//  data = create_time_data(size, data_class);
//  fp = open_file(fname, "r");
//
//  t = 0;
//
//  switch (data_class) {
//  case MACKEY_GLASS: // Mackey-Glass (ID: 1000)
//    while (!feof(fp)) {
//	fgets(line, 512, fp);
//	if (line[0] == '#') continue;
//	if (sscanf(line, "%d%lf", &i, &(data->x[t])) == 2) ++t;
//    }
//    data->t = t-1;
//    break;
//
//  default: // Oops ,,,
//    printf("error: undefined function \"%d\". (%s)\n", data_class, name);
//    exit(-1);
//  }
//
//  close_file(fp);
//  return(data);
//}


/*====================================================================*
 *
 * DATA（NETで扱うデータ構造体）の処理
 *
 *====================================================================*/

/*====================================================================*
 * show_data_parms
 * データを表示する
 *====================================================================*/
void show_data_parms(DATA *data) {
  char str[32];
  int n_channels = data->k;
  int n_total = data->n_total;
  int t, k;

  printf("================================================================\n");
  printf("	Show The Data-Set Parameters\n");
  printf("================================================================\n");
  printf("> data type           : %s\n", data->type);
  printf("> total data          : %d\n", data->n_total);
  printf("> # input channels    : %d\n", data->k);
  printf("> training data       : %d\n", data->n_train);
  printf("> prediction data     : %d\n", data->n_test);
  printf("> max (output) value  : %e\n", data->ymax);
  printf("> min (output) value  : %e\n", data->ymin);
  printf("> function id.        : %d\n", data->data_class);
  printf("----------------------------------------------------------------\n");
  printf("* show each data set? [y/n] : "); str[0]=0;scanf1("%s", str);

  if (str[0] == 'y') {
    for (t=0; t<n_total; t++) {
      printf("%4d : ",t);
      for (k=0; k<=n_channels; k++) printf("x[%d] = %e, ", k, data->x[t][k]);
      //      for (k=0; k<=n_channels; k++) printf("x[%d] = %e, ", k, data->X[t][k]);
      printf("y = %e", data->y[t]);
      printf("\n");
    }
  }
  return;
}

/*====================================================================*
 * get_data_parms
 * 時系列データのパラメータを得る
 *
 * n_channels : input channels
 * n_train    : training steps
 * n_pred     : predict steps
 * n_total    : total steps
 *====================================================================*/
void get_data_parms(int *n_channels, int *n_train, int *n_test, int n_total) {
  char name[32] = "get_net_parms";

  printf("================================================================\n");
  printf("	Input The Data Set Prameters\n");
  printf("================================================================\n");
  printf("> total steps     : %d\n", n_total);
  printf("* training steps? : "); scanf1("%d", n_train);
  printf("* input channels? : "); scanf1("%d", n_channels);

  // エラー処理
  if (*n_train < 1 || n_total <= *n_train) {
    printf("error: training steps is invalied value. (%s)\n", name);
    exit(-1);
  }
  if (*n_train < *n_channels) {
    printf("error: input channels is invalied value. (%s)\n", name);
    exit(-1);
  }

  *n_test = n_total-(*n_train);
  printf("> training steps : %d\n", *n_train);
  printf("> predict steps  : %d\n", *n_test);
  printf("> input channels : %d\n", *n_channels);
  return;
}

//void create_data1(char *type, int data_class,
//		  int n_channels, int n_train, int n_total,DATA *data,int t0) {
void create_data1(int data_class,
		  int n_channels, int n_train, int n_total,DATA *data,int t0) {
  char name[32] = "create_data";
  int t,n_buff;

  // データセット
  //  data = (DATA *)malloc(sizeof(DATA));
  if (data == NULL) {
    printf("error: memory allocation error #1. (%s)\n", name);
    exit(-1);
  }
  //  strcpy(data->type, type);
  data->MSE  = 0.;
  data->NMSE = 0.;
  data->max = 0.;
  data->min = 0.;
  data->mean = 0.;
  data->n_total = n_total;
  data->n_train = n_train;
  data->block = data->block_begin = data->block_end = 0;//ueno(Tue Feb 10 19:36:46 2004)
  data->t0=t0;
  data->k = n_channels;
  //  data->k1 = n_channels+1;
  data->n_test = (n_total-n_train);
  if(data->n_test<=0) data->n_test=1;
  data->data_class = data_class;
  n_buff=n_total+t0;
  data->xmin = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
  data->xmax = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
  data->xmin0 = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
  data->xmax0 = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
  // 入力データ
  //  data->x = (FLOAT **)malloc(sizeof(FLOAT *)*n_total);
//  if(data->x!=NULL){
//    for (t=0; t<n_total; t++) if(data->x[t]!=NULL) free(data->x[t]);
//    free(data->x);
//  }
//  if(data->X!=NULL){
//    for (t=0; t<n_total; t++) if(data->X[t]!=NULL) free(data->X[t]);
//    free(data->X);
//  }
  data->x = (FLOAT **)malloc(sizeof(FLOAT *)*(n_total+1));
  data->X = (FLOAT **)malloc(sizeof(FLOAT *)*(n_total+1));
  if (data->x == NULL) {
    printf("error: memory allocation error #?. (%s)\n", name);
    exit(-1);
  }
  for (t=0; t<n_total; t++) {
    data->x[t] = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    data->X[t] = (FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    if (data->x[t] == NULL) {
      printf("error: memory allocation error #?-%d. (%s)\n", t, name);
      exit(-1);
    }
  }
  //  for (t=0; t<n_total; t++) data->x[t]=data->x0[t+t0];
  // 出力データ
  //  if(data->y0!=NULL){free(data->y0);}
  //  if(data->Y0!=NULL){free(data->Y0); }
  data->y0 = (FLOAT *)malloc(sizeof(FLOAT)*n_buff);
  data->Y0 = (FLOAT *)malloc(sizeof(FLOAT)*n_buff);
  data->y = &(data->y0[t0]);
  data->Y = &(data->Y0[t0]);
  if (data->y0 == NULL) {
    printf("error: memory allocation error #?. (%s)\n", name);
    exit(-1);
  }
  // 予測誤差データ
  data->e = (FLOAT *)malloc(sizeof(FLOAT)*n_total);
  data->rot_x=60;data->rot_z=30;
  if (data->e == NULL) {
    printf("error: memory allocation error #?. (%s)\n", name);
    exit(-1);
  }
}

/*====================================================================*
 * remove_data
 * データセットを破棄する
 *====================================================================*/
void remove_data(DATA *data) {
  int t;

  for (t=0; t<data->n_total; t++) free(data->x[t]);
  free(data->e);
  free(data->y);
  free(data->x);
  free(data);
  return;
}



/*
 * Input Vector belong Voronoi Area for Batch Learning
 */
//typedef struct {
//  int **id;	// id[i][j] : BDATA.x[i][j][k]に代入されるDATA.x[t][k]のインデックスt
//  int *m;	// m[i]     : w[i]のVoronoi領域に属するx[t]の個数
//  int *m0;	// m[i]     : w[i]のVoronoi領域に属するx[t]の個数
//  int *s;	// s[t]     : x[t]の属するVoronoi領域のインデックスi
//  int n;	// number of cells
//  int k;	// input channels
//  int t;	// training steps
//} BIDATA;

/*====================================================================*
 * create_batch_ivector
 * Voronoi領域に属する入力ベクトルのデータを作成
 *====================================================================*/
//BIDATA* create_batch_ivector(int n_cells, int n_channels, int n_train) {
//  char name[32] = "create_batch_ivector";
//  BIDATA *bidata = NULL;
//  int i;
//
//  bidata = (BIDATA *)malloc(sizeof(BIDATA));
//  if (bidata == NULL) {
//    printf("error: memory allocation error #1. (%s)\n", name);
//    exit(-1);
//  }
//
//  bidata->id = (int **)malloc(sizeof(int *)*n_cells);
//  if (bidata->id == NULL) {
//    printf("error: memory allocation error #1. (%s)\n", name);
//    exit(-1);
//  }
//  for (i=0; i<n_cells; i++) {
//    bidata->id[i] = (int *)malloc(sizeof(int)*(n_train+1));
//    if (bidata->id[i] == NULL) {
//	printf("error: memory allocation error #2-%d. (%s)\n", i, name);
//	exit(-1);
//    }
//  }
//
//  bidata->m0 = (int *)malloc(sizeof(int)*n_cells);
//  bidata->m = (int *)malloc(sizeof(int)*n_cells);
//  if (bidata->m == NULL) {
//    printf("error: memory allocation error #2. (%s)\n", name);
//    exit(-1);
//  }
//  bidata->s = (int *)malloc(sizeof(int)*n_train);
//  if (bidata->s == NULL) {
//    printf("error: memory allocation error #3. (%s)\n", name);
//    exit(-1);
//  }
//  return(bidata);
//}

/*====================================================================*
 * remove_batch_ivector
 * Voronoi領域に属する入力ベクトルのデータを破棄
 *====================================================================*/
//void remove_batch_ivector(BIDATA *bidata) {
//  int n_cells = bidata->n_test;
//  int i;
//  for (i=0; i<n_cells; i++) free(bidata->id[i]);
//  free(bidata->s);
//  free(bidata->m);
//  free(bidata->m0);
//  free(bidata->id);
//  free(bidata);
//  return;
//}

/*====================================================================*
 * init_batch_ivector
 * Voronoi領域に属する入力ベクトルのデータを初期化
 *====================================================================*/
//BIDATA *init_batch_ivector(FLOAT **x, FLOAT *y, FLOAT **w,
//			     int n_cells, int n_channels, 
//			     int n_train, NET *net) {
//  char name[32] = "init_batch_ivector";
//  BIDATA *bidata = NULL;
//  FLOAT length, min_length;
//  int i_min_length;
//  int i, j, t;
//
//  bidata = create_batch_ivector(n_cells, n_channels, n_train);
//  bidata->n = n_cells;
//  bidata->k = n_channels;
//  bidata->n_total = n_train;
//
//  // bidata->m，bidata->sの初期化
//  for (i=0; i<n_cells; i++) bidata->m[i] = 0;
//
//  for (t=0; t<n_train; t++) {
//    min_length = calc_length(x[t], w[0], n_channels);
//    i_min_length = 0;
//
//    for (i=0; i<n_cells; i++) {
//	length = calc_length(x[t], w[i], n_channels);
//	if (length < min_length) {
//	  min_length = length;
//	  i_min_length = i;
//	}
//    }
//    bidata->s[t] = i_min_length;
//    bidata->id[i_min_length][bidata->m[i_min_length]]=t;
//    ++(bidata->m[i_min_length]);
//  }
//
////  // bidata->idの初期化
////  for (i=0; i<n_cells; i++) {
////    j = 0;
////    for (t=0; t<n_train; t++) {
////	if (i == bidata->s[t]) {
////	  bidata->id[i][j] = t;
////	  ++j;
////	}
////    }
////  }
//  // Complete the Input Vectors
//  // 荷重ベクトルが参照する入力ベクトルの個数は，
//  // 最低でも入力次数以上ないと，連想行列Mの値が保証されない
//  {
//    FLOAT *l = NULL,ltemp;
//    int *s = NULL, n1 = (n_channels+1), temp;
//    int n2,k,tt,n0,n10;
//    l = (FLOAT *)malloc(sizeof(FLOAT)*n_train);
//    if (l == NULL) {
//	printf("error: memory allocation error #1. (%s)", name);
//	exit(-1);
//    }
//    s = (int *)malloc(sizeof(int)*n_train);
//    if (s == NULL) {
//	printf("error: memory allocation error #2. (%s)", name);
//	exit(-1);
//    }
//
//    //    n1=(n_channels+1)+net->m/net->n_cells;//kuro: larger for remove noise
//    //    n1=(n_channels+1)+50./(GlobalTime+1);//kuro: larger for remove noise
//    //    n1=(n_channels+1)+20./(GlobalTime+1);//kuro: larger for remove noise
//    //    n1=n_channels+1+20./GlobalTime;//kuro: larger for remove noise
//    //    n1=n_channels;  //MSE9.262654e-09;8.725766e-07
//    //    n1=n_channels+1;//MSE9.006930e-09;1.706565e-06
//    for (i=0; i<n_cells; i++) {
//	//      n1=bidata->m[i]+(n_channels+1.)*(1./GlobalTime-1./GlobalTimeMax);//MSE8.274131e-09***;1.448717e-06
//	n1=bidata->m[i]+(n_channels+1.-bidata->m[i])*(1./GlobalTime-1./GlobalTimeMax);//MSE1.009213e-08;5.891452e-07***
//	//      n1=n_channels+1;
//	bidata->m0[i]=bidata->m[i];
//	if (bidata->m[i] < n1) {//入力次数以上ないとき
//	  n0=bidata->m[i];
//	  n2=n1-n0;
//	  //n2=n1*3-n0;
//	  //	n2=n1*1.5-n0;
//	  n10=n1-n0;
//	  tt=0;
//	  //入力ベクトルと荷重ベクトルの距離による入力ベクトルのソート
//	  for (t=0; t<n_train; t++){
//	    if(bidata->s[t]==i) continue;
//	    l[t] = calc_length(x[t], w[i], n_channels);//距離
//	    s[tt++] = t;//ソート準備
//	  }
//	  // ソート
//	  for (t=0; t<n2; t++) {
//	    for (j=t+1; j<tt; j++) {
//	      if (l[(int)s[t]] > l[(int)s[j]]) {
//		temp = s[j];
//		s[j] = s[t];
//		s[t] = temp;
//	      }
//	    }
//	  }
//	  //kuro
//	  if(GlobalTime>1){//kuro
//	    //近似誤差による入力ベクトルのソート
//	    for(j=0;j<n2;j++){
//	      t=s[j];
//	      if(bidata->s[t]==i) l[t]=0;
//	      else{
//		ltemp =y[t];
//		for(k=0; k<=(n_channels); k++) 
//		  ltemp -= net->cell[i].am.M[0][k]*x[t][k];
//		l[t]=fabs(ltemp); //誤差
//	      }
//	    }
//	    for(j=0;j<n10;j++){//Sort
//	      for(k=j+1;k<n2;k++){
//		if(l[s[j]]>l[s[k]]){
//		  temp=s[j];
//		  s[j]=s[k];
//		  s[k]=temp;
//		}
//	      }
//	    }
//	  }
//	  for (j=0; j<n10; j++) bidata->id[i][j+n0] = s[j];
//	  bidata->m[i] = n1;
//	}
//    }
//    free(s);
//    free(l);
//  }
//
//  return(bidata);
//}

/*====================================================================*
 * create_batch_data
 * データセットを作成
 *====================================================================*/
//BDATA* create_batch_data(int n_cells, int *n_ivectors, int n_channels) {
//  char name[32] = "create_batch_data";
//  BDATA *bdata = NULL;
//  int i, j;
//
//  // メモリ確保
//  bdata = (BDATA *)malloc(sizeof(BDATA));
//  if (bdata == NULL) {
//    printf("error: memory allocation error #1. (%s)\n", name);
//    exit(-1);
//  }
//  bdata->m0 = (int *)malloc(sizeof(int)*n_cells);
//  bdata->m = (int *)malloc(sizeof(int)*n_cells);
//  if (bdata->m == NULL) {
//    printf("error: memory allocation error #2. (%s)\n", name);
//    exit(-1);
//  }
//
//  // x[i][j][k]
//  bdata->x = (FLOAT ***)malloc(sizeof(FLOAT **)*n_cells);
//  if (bdata->x == NULL) {
//    printf("error: memory allocation error #3. (%s)\n", name);
//    exit(-1);
//  }
//  for (i=0; i<n_cells; i++) {
//    bdata->x[i] = (FLOAT **)malloc(sizeof(FLOAT *)*n_ivectors[i]);
//    if (bdata->x[i] == NULL) {
//	printf("error: memory allocation error #4-%d. (%s)\n", i, name);
//	exit(-1);
//    }
//  }
//  for (i=0; i<n_cells; i++) {
//    for (j=0; j<n_ivectors[i]; j++) {
//	bdata->x[i][j] = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels+1));
//	if (bdata->x[i][j] == NULL) {
//	  printf("error: memory allocation error #5-%d-%d. (%s)\n", i, j, name);
//	  exit(-1);
//	}
//    }
//  }
//
//  // y[i][j]
//  bdata->y = (FLOAT **)malloc(sizeof(FLOAT *)*n_cells);
//  if (bdata->y == NULL) {
//    printf("error: memory allocation error #6. (%s)\n", name);
//    exit(-1);
//  }
//  for (i=0; i<n_cells; i++) {
//    bdata->y[i] = (FLOAT *)malloc(sizeof(FLOAT)*n_ivectors[i]);
//    if (bdata->y[i] == NULL) {
//	printf("error: memory allocation error #7-%d. (%s)\n", i, name);
//	exit(-1);
//    }
//  }
//
//  // w[i][k]
//  bdata->w = (FLOAT **)malloc(sizeof(FLOAT *)*n_cells);
//  if (bdata->w == NULL) {
//    printf("error: memory allocation error #8. (%s)\n", name);
//    exit(-1);
//  }
//  for (i=0; i<n_cells; i++) {
//    bdata->w[i] = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels));
//    if (bdata->w[i] == NULL) {
//	printf("error: memory allocation error #9. (%s)\n", name);
//	exit(-1);
//    }
//  }
//
//  /*
//  // x_init[i][k]
//  bdata->x_init = (FLOAT **)malloc(sizeof(FLOAT *)*n_cells);
//  if (bdata->x_init == NULL) {
//    printf("error: memory allocation error #8. (%s)\n", name);
//    exit(-1);
//  }
//  for (i=0; i<n_cells; i++) {
//    bdata->x_init[i] = (FLOAT *)malloc(sizeof(FLOAT)*(n_channels+1));
//    if (bdata->x_init[i] == NULL) {
//	printf("error: memory allocation error #9. (%s)\n", name);
//	exit(-1);
//    }
//  }
//
//  // y_init[i]
//  bdata->y_init = (FLOAT *)malloc(sizeof(FLOAT)*n_cells);
//  if (bdata->y_init == NULL) {
//    printf("error: memory allocation error #10. (%s)\n", name);
//    exit(-1);
//  }
//
//  // t_init[i]
//  bdata->t_init = (int *)malloc(sizeof(int)*n_cells);
//  if (bdata->t_init == NULL) {
//    printf("error: memory allocation error #11. (%s)\n", name);
//    exit(-1);
//  }
//  */
//
//  return(bdata);
//}

/*====================================================================*
 * remove_batch_data
 * バッチ型学習用のデータセットを破棄
 *====================================================================*/
//void remove_batch_data(BDATA *bdata) {
//  int i, j;
//
//  for (i=0; i<bdata->n; i++)
//    for (j=0; j<bdata->m[i]; j++)
//	free(bdata->x[i][j]);
//  for (i=0; i<bdata->n; i++) {
//    //free(bdata->x_init[i]);
//    free(bdata->w[i]);
//    free(bdata->y[i]);
//    free(bdata->x[i]);
//  }
//  //free(bdata->t_init);
//  //free(bdata->y_init);
//  //free(bdata->x_init);
//  free(bdata->w);
//  free(bdata->y);
//  free(bdata->x);
//  free(bdata);
//  return;
//}

/*====================================================================*
 * show_batch_data
 * バッチ型学習用のデータセットを表示
 *====================================================================*/
//void show_batch_data_parms(BDATA *bdata) {
//  //char name[32] = "show_batch_data";
//  FILE *gp, *fp_dat, *fp_gpl;
//  char dat[64], gpl[64], obj[64];
//  int n_cells = bdata->n;
//  int i;
//
//  // Display
//  gp = open_pipe(GNUPLOT, "w");
//  fprintf(gp, "set title \"number of input vectors in each voronoi area\"\n");
//  fprintf(gp, "set data style boxes\n");
//  fprintf(gp, "plot '-' using 1:2\n");
//  for (i=0; i<n_cells; i++) fprintf(gp, "%d %d\n", i, bdata->m[i]);
//  fprintf(gp, "e\n");
//  fflush(gp);
//  usleep(2000000);
//  close_pipe(gp);
//
//  // Write Out
//  strcpy(dat, "dist_ivectors.dat");
//  strcpy(gpl, "dist_ivectors.gpl");
//  strcpy(obj, "dist_ivectors.obj");
//  fp_dat = open_file(dat, "w");
//  fp_gpl = open_file(gpl, "w");
//
//  for (i=0; i<n_cells; i++) fprintf(fp_dat, "%d %d\n", i, bdata->m[i]);
//  fprintf(fp_gpl, "set data style boxes\n");
//  fprintf(fp_gpl, "set title \"number of input vectors in each voronoi area\"\n");
//  // for Output X11
//  fprintf(fp_gpl, "set terminal x11\n");
//  fprintf(fp_gpl, "plot \"%s\" using 1:2\n", dat);
//  fprintf(fp_gpl, "pause -1 \"Hit Return to Save %s\"\n", obj);
//  // for Output Tgif
//  fprintf(fp_gpl, "set terminal tgif\n");
//  fprintf(fp_gpl, "set output \"%s\"\n", obj);
//  fprintf(fp_gpl, "plot \"%s\" using 1:2\n", dat);
//
//  close_file(fp_gpl);
//  close_file(fp_dat);
//  return;
//}

/*====================================================================*
 * init_batch_data
 * バッチ型学習用のデータセットを初期化
 *====================================================================*/
//BDATA *init_batch_data(FLOAT **x, FLOAT *y, FLOAT **w,
//			   int n_cells, int n_channels, 
//			 int n_train,
//			 NET *net) {
//  //char name[32] = "init_batch_data";
//  BDATA *bdata = NULL;
//  BIDATA *bidata = NULL;
//  int i, j, k, t;
//
//  // Set the input data ,,,
//
//  // Voronoi領域に属する入力ベクトルのデータ
//  bidata = init_batch_ivector(x, y, w, n_cells, n_channels, n_train,net);
//
//  // データを作成
//  bdata = create_batch_data(n_cells, bidata->m, n_channels);
//
//  // bdata->n，bdata->m，bdata->kの初期化
//  bdata->n = n_cells;
//  for (i=0; i<n_cells; i++){
//    bdata->m[i] = bidata->m[i];
//    bdata->m0[i] = bidata->m0[i];
//  }
//  bdata->k = n_channels;
//  bdata->t = n_train;
//
//  // bdata->wの初期化
//  for (i=0; i<n_cells; i++) {
//    for (k=0; k<(n_channels); k++) {
//	bdata->w[i][k] = w[i][k];
//    }
//  }
//
//  // bdata->x，bdata->yの初期化
//  for (i=0; i<n_cells; i++) {
//    for (j=0; j<bdata->m[i]; j++) {
//	t = bidata->id[i][j];
//	for (k=0; k<=(n_channels); k++) bdata->x[i][j][k] = x[t][k];
//	bdata->y[i][j] = y[t];
//    }
//  }
//
//
//  remove_batch_ivector(bidata);
//  return(bdata);
//}
//#define isspace(c) (c==' ')
//#define isnump(c) ((c >= '0' && c <= '9')||c=='.'||c=='-'||c=='+'||c==' ')
//void prepare_data(int *n_channels, 
//		    int *train_steps, 
//		    int *pred_steps, 
//		    int *total_steps,
//		    int *data_class,
//		    DATA *train, 
//		    DATA *test, 
//		    DATA *pred,
//		    NET *net)
void calc_statistics(DATA *data,NET *net)
{
  int t,k;
  int n_total=data->n_total;
  int n_train=data->n_train;
  int n_channels=data->k;
  int nn;

  data->MEANtrain=nn=0;
  data->ymin=data->ymax=data->y[0];
  for(k=0;k<n_channels;k++){data->xmin[k]=data->xmax[k]=data->x[0][k];}
  for(t=0;t<n_train;t++){//  for(t=1;t<n_train;t++){
    data->MEANtrain+=data->y[t];
    nn++;
    if(data->y[t]<data->ymin) data->ymin=data->y[t];
    if(data->y[t]>data->ymax) data->ymax=data->y[t];
    for(k=0;k<n_channels;k++){
      if(data->x[t][k]<data->xmin[k]) data->xmin[k]=data->x[t][k];
      if(data->x[t][k]>data->xmax[k]) data->xmax[k]=data->x[t][k];
    }
  }
  nn=n_train;//  nn=n_train-1;
  data->MEANtrain/=nn;
  
  data->VARtrain=0;
  for(t=0;t<n_train;t++) //  for(t=1;t<n_train;t++) 
    data->VARtrain+=square(data->y[t]-data->MEANtrain);
  data->VARtrain/=nn;

  if(n_total>n_train){
    data->MEANtest=nn=0;
    for (t=n_train; t<n_total; t++){
      data->MEANtest+=data->y[t];
      nn++;
    }
    nn=n_total-n_train;
    data->MEANtest/=nn;
    data->VARtest=0;
    for (t=n_train; t<n_total; t++) 
      data->VARtest +=square(data->y[t]-data->MEANtest);
    data->VARtest/=nn;
  }
  else data->VARtest=data->VARtrain;//??//??060111
  if(data->VARtest<=0) data->VARtest=data->VARtrain;//??060111
  net->printf(strx2("ymin=%13.7e;ymax=%13.7e;\n", data->ymin,data->ymax));
  for(k=0;k<n_channels;k++){
    //    printf("x[%d] is in [%f,%f].\n",k,data->xmin[k],data->xmax[k]);
    net->printf(strx4("xmin[%d]=%13.7e;xmax[%d]=%13.7e\n",k,data->xmin[k],k,data->xmax[k]));
  }
  net->printf(strx2("VARtrain=%e. VARtest=%e.\n",data->VARtrain,data->VARtest));
}


void normalize_data(DATA *givendata, NET *net)
{
  int t,k,i;
  int n_total=givendata->n_total;
  //  int n_train=givendata->n_train;
  int n_channels=givendata->k;
  //  int t0=givendata->t0;
  int data_class=givendata->data_class;
  //  int nn;
  
  net->printf("--- tips ---\n");
  //  printf("Original xmin,xmax,ymin,ymax=%f,%f,%f,%f.\n",
  //	 givendata->xmin,givendata->xmax,givendata->ymin,givendata->ymax);
  net->printf("Normalize from [min0,max0] to [min,max].\n");
  net->printf("No normalization if ymin==0 && ymax==0.\n");
  net->printf("Enter ymin0 ymax0 ymin ymax:");
  double _ym0=0,_yM0=0,_ym1=0,_yM1=0,_y0m=0,_y1m=0,_y0M=0,_y1M=0,_yhh=1.;
  double _xm0=0,_xM0=0,_xm1=0,_xM1=0,_x0m=0,_x1m=0,_x0M=0,_x1M=0,_xhh=1.;
  {
    char buff[256],*p=buff;
    fgets(buff,256,stdin);
    for(;;p++){if(*p==0x00 || *p==':') break;}
    if(*p!=':'){
      //      sscanf(buff,"%lf%lf%lf%lf",&net->ymin0,&net->ymax0,&net->ymin,&net->ymax);//scanf4("%lf%lf%lf%lf",&net->ymin0,&net->ymax0,&net->ymin,&net->ymax);
      sscanf(buff,"%lf%lf%lf%lf",&_ym0,&_yM0,&_ym1,&_yM1);
    }
    else{
      //      sscanf(buff,"%lf%lf%lf%lf %lf:%lf:%lf:%lf",&_ym0,&_yM0,&net->ymin,&net->ymax,&_y0m,&_y1m,&_y0M,&_y1M);
      sscanf(buff,"%lf%lf%lf%lf %lf:%lf:%lf:%lf",&_ym0,&_yM0,&_ym1,&_yM1,   &_y0m,&_y1m,&_y0M,&_y1M);
    }
  }
  if(net->init==0){
    net->xmin0=(FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    net->xmax0=(FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    net->xmin=(FLOAT *)malloc(sizeof(FLOAT)*n_channels);
    net->xmax=(FLOAT *)malloc(sizeof(FLOAT)*n_channels);
  }
  if(_yM0-_ym0<1e-20 && _yM1-_ym1<1e-20){
    _ym0=givendata->ymin;_yM0=givendata->ymax;_yhh=_yM0-_ym0;
    net->ymin0=net->ymin=_ym0-_yhh*_y1m-_y0m;//slightly smaller ?
    net->ymax0=net->ymax=_yM0+_yhh*_y1M+_y0M;//slightly bigger ?
    _x0m=_y0m;_x1m=_y1m;_x0M=_y0M;_x1M=_y1M;
    for(k=0;k<n_channels;k++){
      _xm0=givendata->xmin[k];
      _xM0=givendata->xmax[k];
      _yhh=_yM0-_ym0;
      net->xmin0[k]=net->xmin[k]=_xm0-_xhh*_x1m-_x0m;
      net->xmax0[k]=net->xmax[k]=_xM0+_xhh*_x1M+_x0M;
    }
    net->printf(strx4("#Modified-y:%f:%f:%f:%f,",net->ymin0,net->ymax0,net->ymin,net->ymax));
    net->printf(strx4("#Modified-x0:%f:%f:%f:%f\n",net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]));
  }
  else{//
    if(_yM0-_ym0<1e-20){_ym0=givendata->ymin;_yM0=givendata->ymax;}
    _yhh=_yM0-_ym0;
    net->ymin0=_ym0-_yhh*_y1m-_y0m;
    net->ymax0=_yM0+_yhh*_y1M+_y0M;

    if(_yM1-_ym1<1e-20){_ym1=givendata->ymin;_yM1=givendata->ymax;}
    _yhh=_yM1-_ym1;
    net->ymin=_ym1-_yhh*_y1m-_y0m;
    net->ymax=_yM1+_yhh*_y1M+_y0M;
    //    fprintf(stderr,"y:%f:%f:%f:%f ",net->ymin0,net->ymax0,net->ymin,net->ymax);
    if(data_class == TIME_SERIES){
      for(k=0;k<n_channels;k++){
	net->xmin0[k]=net->ymin0;
	net->xmax0[k]=net->ymax0;
	net->xmin[k]=net->ymin;
	net->xmax[k]=net->ymax;
      }
    }
    else{//FUNCTION_APPROX
      for(k=0;k<n_channels;k++){
	net->printf(strx1("Enter x[%d]min0 max0 min max:",k));
	{
	  char buff[256],*p=buff;
	  fgets(buff,256,stdin);
	  for(;;p++){if(*p==0x00 || *p==':') break;}
	  if(*p!=':'){
	    //	    sscanf(buff,"%lf%lf%lf%lf",&(net->xmin0[k]),&(net->xmax0[k]),&(net->xmin[k]),&(net->xmax[k]));}
	    sscanf(buff,"%lf%lf%lf%lf",&_xm0,&_xM0,&_xm1,&_xM1);
	  }
	  else{
	    //	    sscanf(buff,"%lf%lf%lf%lf %lf:%lf:%lf:%lf",&_xm0,&_xM0,&(net->xmin[k]),&(net->xmax[k]),&_x0m,&_x1m,&_x0M,&_x1M);
	    sscanf(buff,"%lf%lf%lf%lf %lf:%lf:%lf:%lf",&_xm0,&_xM0,&_xm1,&_xM1,&_x0m,&_x1m,&_x0M,&_x1M);
	  }
	  //	  _xhh=(_xM0-_xm0);
	  //	  net->xmin0[k]=_xm0-_xhh*_x1m-_x0m;
	  //	  net->xmax0[k]=_xM0+_xhh*_x1M+_x0M;
	}
	if(_xM0-_xm0<1e-20){//net->xmax0[k]-net->xmin0[k]<1e-20){
	  _xm0=givendata->xmin[k];
	  _xM0=givendata->xmax[k];
	}
	_xhh=_xM0-_xm0;
	net->xmin0[k]=_xm0-_xhh*_x1m-_x0m;
	net->xmax0[k]=_xM0+_xhh*_x1M+_x0M;

	if(_xM1-_xm1<1e-20){//net->xmax0[k]-net->xmin0[k]<1e-20){
	  _xm1=givendata->xmin[k];
	  _xM1=givendata->xmax[k];
	}
	_xhh=_xM1-_xm1;
	net->xmin[k]=_xm1-_xhh*_x1m-_x1m;
	net->xmax[k]=_xM1+_xhh*_x1M+_x1M;
      }
    }
    net->printf(strx4("#Modified-y:%f:%f:%f:%f,",net->ymin0,net->ymax0,net->ymin,net->ymax));
    net->printf(strx4("#Modified-x0:%f:%f:%f:%f\n",net->xmin0[0],net->xmax0[0],net->xmin[0],net->xmax[0]));
    for(t=0;t<n_total;t++){
      givendata->y[t]=moverange(givendata->y[t],
				net->ymin0,net->ymax0,
				net->ymin,net->ymax);
      for(k=0;k<n_channels;k++){
	givendata->x[t][k]=moverange(givendata->x[t][k],
				     net->xmin0[k],net->xmax0[k],
				     net->xmin[k],net->xmax[k]);
      }
    }
    net->printf("--- tips for resolution ---\n");
    net->printf("r1 r2 (r1/r2 for integers r1 and r2 is the resolution).");
    net->printf("no digitization if r1=0\n");
    net->printf("r1 r2 r3:");
    //    scanf2("%d%d",&net->r1,&net->r2);
    scanf3("%d%d%lf",&net->r1,&net->r2,&net->r3);

    net->ywidth=net->ymax-net->ymin;
    if(net->r1>0){// resolution
      char buff[126];
      net->r12=((double)net->r1/net->r2)*(net->ywidth)/(net->ymax0-net->ymin0);
      net->R12=((double)net->r1/net->r2);

      double R3=net->r3, r3=R3*(net->ywidth)/(net->ymax0-net->ymin0);

      net->nr=(int)((net->ywidth)/net->r12+1);//少し大きめ
      net->r=(FLOAT *)malloc(sizeof(FLOAT)*(net->nr+1));
      net->R=(FLOAT *)malloc(sizeof(FLOAT)*(net->nr+1));
      if(fabs(net->r3)<1e-20){//no round
	for(i=0;i<=net->nr;i++){
	  sprintf(buff,"%.10e",(double)i*net->r12+net->ymin);
	  sscanf(buff,"%lf",&net->r[i]);
	  sprintf(buff,"%.10e",(double)i*net->R12+net->ymin0);
	  sscanf(buff,"%lf",&net->R[i]);
	  //	printf("???R[%d]=%f\n",i,net->R[i]);
	}
      }
      else{//if(fabs(res->r3)>ZERO) for round at r3
	for(i=0;i<=net->nr;i++){
	  sprintf(buff,"%.10e",(double)((int)((i*net->r12+net->ymin)/r3+0.5))*r3);
	  sscanf(buff,"%lf",&net->r[i]);
	  sprintf(buff,"%.10e",(double)((int)((i*net->R12+net->ymin0)/R3+0.5))*R3);
	  sscanf(buff,"%lf",&net->R[i]);
	  //	printf("???R[%d]=%f\n",i,net->R[i]);
	}
      }
      for(t=0;t<n_total;t++){
	givendata->y[t]=net->r[(int)((givendata->y[t]-net->ymin)/net->r12+0.5)];
	//	printf("???y[%d]=%f\n",t,givendata->y[t]);
      }
      if(data_class==TIME_SERIES){
	for(t=0;t<n_total;t++){
	  for(k=0;k<n_channels;k++){
	    givendata->x[t][k]=net->r[(int)((givendata->x[t][k]-net->ymin)/net->r12+0.5)];
	  }
	}
      }
      else{
	printf("not digitize x for function approximation\n");
      }
    }
//    if(net->r3>0){//net->r3>0 for nonlinear parameter;old version modified at 060211
//      for(t=0;t<n_total;t++){
//	givendata->y[t]=net->ywidth*pow((givendata->y[t]-net->ymin)/net->ywidth,net->r3);
//      }
//    }
  }
  //  if(1==1){//for datacheck nips04 
  {
    FILE *fp;
    if(opendir("./tmp/")==NULL) system("mkdir ./tmp");
    //    fp=fopen("./tmp/nips04.dat","w");
    fp=fopen("./tmp/train+test.dat","w");
    for(t=0;t<n_total;t++){
      for(k=0;k<n_channels;k++){
	fprintf(fp,"%.7e ",givendata->x[t][k]);
      }
      if(t<givendata->n_train) fprintf(fp,"%.7e\n",givendata->y[t]);
      else fprintf(fp,"%.7e\n",givendata->y[t]);
      //     else fprintf(fp,"%.7e #\n",givendata->y[t]);
    }
    fclose(fp);
  }
  //以上 030718m 正規化の試行
}
//#define LINESIZE 5120
//#define LINESIZE 10240
void load_data_IJCNN04(DATA *givendata, DATA *test, NET *net,int *n_train,int *n_total, int n_channels)
{
  int i,j,k;
  char *line=(char *)malloc(sizeof(char)*LINESIZE);//  char line[LINESIZE];
  int i_testblock;//
  FILE *fp;
  int i1,i2,ii,t,tt;
  //  FLOAT dd[5000];
  int t2,t_t1,t_t2,t_d=1000;
  int k2=net->k2;//MAモデル入力次数 not implemented yet
  FLOAT *ds;//data from smooth.dat
  FLOAT *dr;///data from smooth_.dat, or the rest of smooth.dat
  FLOAT *dp;//data of prediction write out to predict.dat
  FLOAT *dt;//data of data.txt and linear approximation is done
  FLOAT *dd;//data of linear approximation data.dat which is from dataconv0 or dataconv0_ with dlast
  FLOAT *dl;//data of linear approximation with data of n_train-1 and n_total+1
  char *fname_s;
  char *fname_r;
  char *fname_p;
  char *fname_g;
  char *fname_t;
  char *fname_d;
  char *data_path;

  printf("Test data are from <i_testblock*1000+t1> to <i_testblock*1000+t2>. "); 
  printf(" For search parameters, set t1=960 and t2=980 (or t2=1000). "); 
  printf(" For final prediction, set t1=980 t2=1000. "); 
  //  printf(" i_testblock(0-4), t1, t2 : "); //
  //  scanf3("%d%d%d", &i_testblock,&t_t1,&t_t2);//
  printf(" i_testblock(0-4), t1 : "); //
  scanf2("%d%d", &i_testblock,&t_t1);//
  givendata->i_testblock=i_testblock;
  t_t2=t_t1+19;
  //  k2=0; //k2 is for ARMA model but has not implemented yet
  *n_train=(980-n_channels)*4+(t_t1-1-n_channels);
  *n_total=*n_train+(t_t2-t_t1+1);//
  t_d=1000-t_t2;

  create_data1(data_class,n_channels,*n_train,(*n_total+t_d),givendata,n_channels);//20131102
  create_data1(data_class,n_channels,*n_train,(*n_total+t_d),test,n_channels);//20131102
  //  create_data1(data_class,n_channels,*n_train,(*n_total+t_d),givendata,n_channels);
  //  create_data1(data_class,n_channels,*n_train,(*n_total+t_d),test,n_channels);

  //givendata->n_total=test->n_total=*n_total;
  //givendata->n_test =test->n_test =*n_total-*n_train;  

  givendata->block = i_testblock;//ueno(Tue Feb 10 19:36:56 2004)
  givendata->block_begin = t_t1+i_testblock*1000;
  givendata->block_end   = t_t2+i_testblock*1000;
  givendata->application = application;//ueno(Mon Feb 16 23:27:06 2004)
  givendata->ijcnn04data=(IJCNN04DATA *)malloc((sizeof(IJCNN04DATA)));
  givendata->ijcnn04data->i_testblock=i_testblock;
  givendata->ijcnn04data->t_t1=t_t1;
  ds=givendata->ijcnn04data->ds;
  dr=givendata->ijcnn04data->dr;
  dp=givendata->ijcnn04data->dp;
  dt=givendata->ijcnn04data->dt;
  dl=givendata->ijcnn04data->dl;
  dd=givendata->ijcnn04data->dd;
  fname_s=givendata->ijcnn04data->fname_s;
  fname_r=givendata->ijcnn04data->fname_r;
  fname_p=givendata->ijcnn04data->fname_p;
  fname_g=givendata->ijcnn04data->fname_g;
  fname_t=givendata->ijcnn04data->fname_t;
  fname_d=givendata->ijcnn04data->fname_d;
  data_path=givendata->ijcnn04data->data_path;

  sprintf(data_path,"%s",givendata->path);
  //  fp = open_file(givendata->path, "r");
  sprintf(fname_p,"%s","predict.dat");
  sprintf(fname_g,"%s","predict.gp");

  //read smooth_.dat
  sprintf(fname_r,"%s/%s",givendata->path,"smooth_.dat");
  fp = open_file(fname_r, "r");
  for(i=1;i<=5000;i++){
    fgets(line,LINESIZE,fp);
    sscanf(line,"%lf",&dr[i]);
    if(feof(fp)) break;
  }
  fclose(fp);
  //read smooth.dat
  {
    sprintf(fname_s,"%s/%s",givendata->path,"smooth.dat");
    if((fp=open_file(fname_s,"r"))==NULL) return;
    for(t=1;t<=5000;t++){
      fgets(line, 512, fp);
      sscanf(line,"%lf",&ds[t]);
    }
    fclose(fp);
  }
  //read data.dat
  {
    sprintf(fname_d,"%s/%s",givendata->path,"data.dat");
    if((fp=open_file(fname_d,"r"))==NULL) return;
    for(t=1;t<=5000;t++){
      fgets(line, 512, fp);
      sscanf(line,"%lf",&dd[t]);
    }
    fclose(fp);
  }
  //read data.txt
  {
    FLOAT a,b;
    sprintf(fname_t,"%s/%s",givendata->path,"data.txt");
    if((fp=open_file(fname_t,"r"))==NULL) return;
    for(t=1;t<=5000;t++) dt[t]=1e30;
    dt[5000]=dd[5000];
    for(t=1;t<=5000;t++){
      line[0]=0;
      fgets(line, 512, fp);
      if(feof(fp)) break;
      if(line[0]==0) break;
      int cr=sscanf(line,"%lf%lf",&a,&b);
      if(cr<2) continue;
      dt[(int)a]=b;
    }
    fclose(fp);
    for(t=1;t<=5000;t++){
      if(dt[t]<1e29) a=dt[t];
      else if(dt[t]>1e29){
	for(j=1;;j++) if(dt[t+j]<1e29) break;
	b=(dt[t+j]-a)/(j+1);
	for(k=1;k<=j;k++){
	  dt[t+k-1]=a+k*b;
	}
      }
      //      printf("--->dt[%d]=%e\n",t,dt[t]);
    }
  }
  //calc AVEdt dl[t]
  {
    FLOAT a,b,MSEdl;
    int  t0=givendata->i_testblock*1000+t_t1;//real start time for prediction in all data
    int  t1=t0+givendata->n_test;//real start time for prediction in all data
    int n_test=givendata->n_test;
    a=dt[t0-1];
    b=(dt[t1]-dt[t0-1])/(n_test+1);
    MSEdl=0;
    for(t=t0;t<t1;t++){
      dl[t]=a+(t-t0+1)*b; //linear approximation
      MSEdl +=square(dl[t]-dt[t]); //linear approximation
    }
    MSEdl/=n_test;
    givendata->ijcnn04data->MSEdl=MSEdl;
  }

  t=0; i=0;
  /////////////////////////////
  for(ii=0,i=0;i<=4;i++){
    i1=i*1000;//original block
    if(i==i_testblock) {// test block 
      t2=t_t1-1;    //last # of the training data in the test block 
      i2=(980-n_channels)*4;//start # of test data in the memory
    }
    else{//not test_block
      t2=980;    //last # of the data in the training block
      i2=ii*(980-n_channels);//start # of training data in the memory
      ii++;
    }

    for(tt=n_channels+1;tt<=t2;tt++){//training data
      t=tt+i1;
      k=tt+i2-(n_channels+1);
      givendata->x[k][n_channels]=BIAS;
      givendata->y[k]=dr[t-k2];
      for (j=0; j<k2; j++) givendata->x[k][j] = dr[t-j];//ARMA
      //for (j=0; j<k2; j++) givendata->x[k][j] = t/5000.;//Jikantest
      for (j=k2; j<n_channels; j++) givendata->x[k][j] = dr[t-1-j];
      //	printf("\nd[%d]=y[%d]=%+e for train",t-k2,k,givendata->y[k]);
    }
    if(i==i_testblock){//test data
      //      for(tt=t_t1;tt<=t_t2;tt++){
      for(tt=t_t1;tt<=1000;tt++){
	t=tt+i1;
	k=tt+i2-(n_channels+1);
	//k=tt+i2-(n_channels);//??
	givendata->x[k][n_channels]=BIAS;
	givendata->y[k]=dr[t];
	for (j=0; j<k2; j++) givendata->x[k][j] = 0;//ARMA
	for (j=k2; j<n_channels; j++) givendata->x[k][j] = dr[t-1-j];//will not be used
	printf("\nd[%d]=y[%d]=%+e for test",t,k,givendata->y[k]);
      }
    }
  }
  
//  for(i=0;i<*n_total;i++){
//    printf("\n[%d]:y=%+11.6fx=",i,givendata->y[i]);
//    //    for(j=0;j<=n_channels;j++) printf("%+11.6f ",givendata->x[i][j]);
//    for(j=0;j<=5;j++) printf("%+11.6f ",givendata->x[i][j]);
//  }
  printf("n_train=%d,n_total=%d\n",*n_train,*n_total);
}

void load_data_RANGEDATA(DATA *givendata, DATA *test, NET *net,int *n_train,int *n_total, int n_channels)
{
  char *line=(char *)malloc(sizeof(char)*LINESIZE); //char line[LINESIZE];
  FILE *fp;
  int t,tt;
  FLOAT xmin,xmax,zmin,zmax,xx,zz,rr;
  int t0;
  FLOAT xmin0,xmax0,zmin0,zmax0,dx,nx,dz,nz;
  FLOAT d,dmin;
  int tmin;
  //  int resolution=50;//51;
  //    int resolution=100;//51;
  //    int resolution=99;//51;
    int resolution=33;//51;
  //int resolution=2;//50;//51;
  //  int SGN=1;  //  int SGN=-1;

  {
    FLOAT d1,d2,d3,d4;
    printf("x1min,x1max,x2min,x2max (range for test):");
    //scanf4("%lf%lf%lf%lf",&xmin,&xmax,&zmin,&zmax);//test
    scanf4("%lf%lf%lf%lf",&d1,&d2,&d3,&d4);//test
    if(fabs(d1-d2)>1e-10){
      xmin=d1;xmax=d2;zmin=d3;zmax=d4;
    }
    else {
      xmin=zmin=-1e-20;
      xmax=zmax=1e20;
    }
  }
  t=tt=0; 
  fp = open_file(givendata->path, "r");
  while(!feof(fp)){
    line[0]=0;
    fgets(line, LINESIZE, fp);
    if(line[0] == 0) break;
    if(feof(fp)) break;
    if(line[0] == '#') continue;
    //    if(!isnump(line[0])) continue;
    int nr=sscanf(line,"%lf%lf%lf",&xx,&zz,&rr);
    if(nr<3) continue;
    if(t==0 || xx<xmin0) xmin0=xx;
    if(t==0 || xx>xmax0) xmax0=xx;
    if(t==0 || zz<zmin0) zmin0=zz;
    if(t==0 || zz>zmax0) zmax0=zz;
    if(xx>=xmin && xx<=xmax && zz>=zmin && zz<=zmax){
      //      printf("t=%d xyz=%e %e %e, xxminmax=%e %e %e %e\n",t,xx,zz,rr,xmin,xmax,zmin,zmax);//debug
      t++;
    }
  }
  printf("\nOriginall x1min,x1max,x2zmin,x2max=%+6.1f %+6.1f %+6.1f %+6.1f\n",xmin0,xmax0,zmin0,zmax0);
  nx=nz=resolution;
  dx=(xmax-xmin)/nx;
  dz=(zmax-zmin)/nz;
  tt=0;
  for(xx=xmin;xx<xmax+dx;xx+=dx){
    for(zz=zmin;zz<zmax+dz;zz+=dz){
      tt++;
    }
  }
  //  tt=1;//??070712
  fclose(fp);
  t0=0;//??
  *n_total=t+tt-t0;
  *n_train=t-t0;
  //  printf("\nOriginall x1min,x1max,x2zmin,x2max=%+6.1f %+6.1f %+6.1f %+6.1f\n",xmin0,xmax0,zmin0,zmax0);
  net->printf(strx2("Number of training and test data=%d and %d.\n",*n_train,*n_total-*n_train));
      
  create_data1(data_class,n_channels,*n_train,*n_total,givendata,t0);
  create_data1(data_class,n_channels,*n_train,*n_total,test,t0);

  fp = open_file(givendata->path, "r");
  t=tt=0; 
  while(!feof(fp)){
    line[0]=0;
    fgets(line, LINESIZE, fp);
    if(line[0] == 0) break;
    if(feof(fp)) break;
    if(!isnump(line[0])) continue;
    int cr=sscanf(line,"%lf%lf%lf",&xx,&zz,&rr);
    if(cr<3) continue;
    if(xx>=xmin && xx<=xmax && zz>=zmin && zz<=zmax){
      givendata->x[t][0]=xx;
      givendata->x[t][1]=zz;
      givendata->y[t]=rr;
      givendata->x[t][n_channels]=BIAS;
//	if(p2r){//polar -> rectangle
//	  FLOAT theta,lll,xxx,yyy,zzz;
//	  theta=(90.0-xx)*PI/180.;
//	  zzz=rr*cos((90.0-zz)*PI/180.);
//	  lll=sqrt(rr*rr-zzz*zzz);
//	  xxx=lll*cos(theta);
//	  yyy=lll*sin(theta);
//	  if(xxx<xxxmin) xxxmin=xxx;
//	  if(xxx>xxxmax) xxxmax=xxx;
//	  if(yyy<yyymin) yyymin=yyy;
//	  if(yyy>yyymax) yyymax=yyy;
//	}
//	else {//rectangle -> polar
//	  FLOAT theta1,theta2,rrr;
//	  rrr=sqrt(xx*xx+zz*zz+rr*rr);
//	  theta1=90.-atan2(xx,rr)*180./PI;
//	  theta2=90.-acos(zz/rrr)*180./PI;
//	  if(theta1<theta1min) theta1min=theta1;
//	  if(theta1>thata1max) theta1max=theta1;
//	  if(theta2<theta2min) theta2min=theta2;
//	  if(theta2>thata2max) theta2max=theta2;
//	}
      //      givendata->y[t]=-rr;
      //      printf("t=%d xyz=%e %e %e, xxminmax=%e %e %e %e\n",t,xx,zz,rr,xmin,xmax,zmin,zmax);//debug
      t++;
    }
  }
  tt=t;
  for(xx=xmin;xx<xmax+dx;xx+=dx){//test data
    //printf("xx,zz=%f,%f\n",xx,zz);
    //          printf("xi=%f\n",xx);
    //    putchar('.');
    for(zz=zmin;zz<zmax+dz;zz+=dz){
      givendata->x[tt][0]=xx;
      givendata->x[tt][1]=zz;
      givendata->x[tt][n_channels]=BIAS;
      //      if(tt-*n_train>10000) printf("%d)xx,zz=%f,%f\n",tt,xx,zz);
      tmin=0;dmin=-1;
      for(t=0;t<*n_train;t++){
	d=square(xx-givendata->x[t][0])+square(zz-givendata->x[t][1]);
	if(dmin<0 || d<dmin){
	  dmin=d;
	  tmin=t;
	}
      }
      givendata->y[tt]=givendata->y[tmin];
      tt++;
    }
  }
  fclose(fp);
  if(1==1){//for checkdata
    FILE *fp;
    int k;
    fp=fopen("rangedataact.dat","w");
    //    for(t=0;t<*n_total;t++){
    for(t=0;t<*n_train;t++){
      for(k=0;k<n_channels;k++){
	fprintf(fp,"%13.7e ",givendata->x[t][k]);
      }
      fprintf(fp,"%13.7e\n",givendata->y[t]);
    }
    fclose(fp);
  }
}

//void load_data_RANGEDATA1(DATA *givendata, DATA *test, NET *net,int *n_train,int *n_total, int n_channels)
//{
//  char line[LINESIZE];
//  FILE *fp;
//  int t,tt;
//  FLOAT xmin,xmax,zmin,zmax,xx,zz,rr;
//  int t0;
//  FLOAT xmin0,xmax0,zmin0,zmax0;
//
//  printf("x1min,x1max,x2min,x2max (range for test):");scanf4("%lf%lf%lf%lf",&xmin,&xmax,&zmin,&zmax);//test
//  t=tt=0; 
//  fp = open_file(givendata->path, "r");
//  while(!feof(fp)){
//    line[0]=0;
//    fgets(line, LINESIZE, fp);
//    if(line[0] == 0) break;
//    if(feof(fp)) break;
//    sscanf(line,"%lf%lf%lf",&xx,&zz,&rr);
//    if(t==0 || xx<xmin0) xmin0=xx;
//    if(t==0 || xx>xmax0) xmax0=xx;
//    if(t==0 || zz<zmin0) zmin0=zz;
//    if(t==0 || zz>zmax0) zmax0=zz;
//    if(xx>=xmin && xx<=xmax && zz>=zmin && zz<=zmax){
//	t++;
//    }
//  }
//
//  t0=0;//??
//  *n_total=t+t-t0;
//  *n_train=t-t0;
//  printf("x1min,x1max,x2zmin,x2max=%+6.1f %+6.1f %+6.1f %+6.1f\n",xmin0,xmax0,zmin0,zmax);
//  printf("Number of training and test data=%d and %d.\n",*n_train,*n_total-*n_train);
//	
//  create_data1(data_class,n_channels,*n_train,*n_total,givendata,t0);
//  create_data1(data_class,n_channels,*n_train,*n_total,test,t0);
//
//  fp = open_file(givendata->path, "r");
//  t=0; tt=*n_train;
//  while(!feof(fp)){
//    line[0]=0;
//    fgets(line, LINESIZE, fp);
//    if(line[0] == 0) break;
//    if(feof(fp)) break;
//    sscanf(line,"%lf%lf%lf",&xx,&zz,&rr);
//    if(xx>=xmin && xx<=xmax && zz>=zmin && zz<=zmax){
//	givendata->x[tt][0]=givendata->x[t][0]=xx;
//	givendata->x[tt][1]=givendata->x[t][1]=zz;
//	givendata->x[tt][2]=givendata->x[t][2]=1.0;
//	givendata->y[tt]=givendata->y[t]=rr;
//	//      givendata->y[t]=-rr;
//	t++;tt++;
//    }
//  }
//  fclose(fp);
//}
void load_data(DATA *givendata, DATA *test, NET *net)
{
  int n_channels;
  int total_steps = 0, train_steps = 0;//, pred_steps = 0;
  char fname[2][256]={0};
  FILE *fp;
  char *line=(char *)malloc(sizeof(char)*LINESIZE);//  char line[LINESIZE];
  int t, i,t1,k;
  int t0,n_total,n_train;
  char *p;
  givendata->tp0=0;
  application=DIRECT_APPLI;
  net->printf("Method(0:time series, 1:funciton approximation,3:ijcnn04,4:RangeData)? "); scanf1("%d", &data_class);
  if(data_class==-1){
    return;
  }
  if(data_class==0){
    net->printf("Data is for Time Series\n");
  }
  else if(data_class==1){
    net->printf("Data is for Function Approximation\n");
  }
  else if(data_class==3){ 
    data_class=TIME_SERIES;
    application=IJCNN04_APPLI;
    net->printf("JobClass is IJCNN04. Data is for Time Series\n");
  }
  else if(data_class==4){ 
    data_class=FUNCTION_APPROX;
    application=RANGE_APPLI;
    net->printf("JobClass is Range Data. Data is for Function Approximation\n");
  }
  //  printf("TimeSeries(0) or funciton approximation(1), and precision (1:int)? "); scanf2("%d%d", &data_class,&net->prec);
  //  net->r1=0;
  //  printf("Number of input channels? : "); scanf3("%d%d%d", &n_channels,&net->r1,&net->r2);
  net->printf("Number of input channels k1 and k2 ? : ");scanf2("%d%d", &net->k1,&net->k2);
  net->k=n_channels=net->k1+net->k2;
  switch (data_class) {
    int ntotal;
  case TIME_SERIES: // 時系列
    net->printf("Data file name ? : "); scanf1("%s", fname[0]);//train
    strcpy(givendata->path, fname[0]);
    //    printf("#datafilename:%s=%s\n",givendata->path,fname[0]);//check 20191114
    t0=n_channels;//??
    if(application==IJCNN04_APPLI) load_data_IJCNN04(givendata,test,net,&n_train,&n_total,n_channels);
    else{//DIRECT_APPLI
      // 以下、ファイルを開き行数をカウント
      t=0; i=0;
      fp = open_file(givendata->path, "r");
      //      while(!feof(fp)){
      while(1){
	fgets(line, LINESIZE, fp);
	if(feof(fp)) break;
	//	if (!isnump(line[0])) continue;//
	if(line[0] == '#') continue;
	++t;
      }
      fclose(fp);
      total_steps = t;
      // 以上、ファイルを開き行数をカウント
      printf("Points in time of training, total and offset, data(Ttrain Ttotal<%d)? : ",total_steps); 
      //      scanf2("%d%d", &train_steps,&ntotal);//pred
      {
	int noffset=0; 
	char buff[256],*p=buff;
	fgets(buff,256,stdin);
	// pattern1 t1 t2 toffset
	// pattern2 tr0-tr1:tp0-tp1
	for(;;p++){if(*p==0x00 || *p==':' || *p=='#') break;}
	if(*p=='#') *p=0;
	if(*p!=':'){
	  sscanf(buff, "%d%d%d", &train_steps,&ntotal,&noffset);//modified 060110
	  //	  scanf3("%d%d%d", &train_steps,&ntotal,&noffset);//modified 060110
	  if(noffset<0){ noffset*=(-1);train_steps-=noffset;ntotal-=noffset;}//modified 061114
	  if(noffset<0 || noffset>total_steps) noffset=0;
	  if(train_steps>total_steps) train_steps=total_steps;
	  if(ntotal>=total_steps) total_steps=ntotal;//061109
	  if(ntotal>=train_steps && ntotal<=total_steps) total_steps=ntotal;
	  n_total=total_steps-t0;
	  n_train=train_steps-t0;
	  givendata->tp0=n_train+n_channels;//??
	  //    pred_steps=n_total-n_train;
	  printf("Training and Test Steps=%d and %d.\n",n_train, n_total-n_train);
	  create_data1(data_class,n_channels,n_train,n_total,givendata,t0);
	  create_data1(data_class,n_channels,n_train,n_total,test,t0);
	  givendata->block_begin = train_steps;//ueno(Tue Feb 10 19:36:56 2004)
	  givendata->block_end   = ntotal;
	  givendata->application = application;//ueno(Mon Feb 16 23:27:06 2004)
	  t=-noffset;//      t=0; 
	  i=0; fp = open_file(fname[i], "r");
	  while(1){
	    fgets(line, LINESIZE, fp);
	    if(feof(fp)) break;
	    if (!isnump(line[0])) continue;//	if(line[0] == '#') continue;
	    if(t<0) {t++;continue;}
	    int nr=sscanf(line, "%lf", &(givendata->y0[t]));
	    if(nr<1) continue;
	    if(++t>=total_steps) break;
	  }
	  fclose(fp);
	  for(t=n_channels;t<total_steps;t++){
	    t1=t-n_channels;
	    givendata->x[t1][n_channels]=BIAS;
	    for (k=0; k<n_channels; k++) {
	      givendata->x[t1][k] = givendata->y0[t-k-1];//bad??
	    }
	  }
	}
	else{ //if(*p!=':'){ //first char non-digit is not ':'
	  int _tr0,_tr1,_tp0,_tp1,_tpD=0,_tpG=0,_tpEy=15,col;
	  col=sscanf(buff,"%d-%d:%d-%d:%d:%d:%d",&_tr0,&_tr1,&_tp0,&_tp1,&_tpD,&_tpG,&_tpEy);//sscanf(buff,"%d-%d:%d-%d",&_tr0,&_tr1,&_tp0,&_tp1);
	  //	  if(col<=4) {_tpD=_tpG=0;_tpEy=15;}
	  //	  if(col<=5) {_tpG=0;} 
	  if(_tp0<n_channels) _tp0=n_channels;
	  if(_tr1>total_steps) _tr1=total_steps;
	  if(_tr0<0) _tr0=0;
	  n_train=_tr1-_tr0-n_channels-_tpD;//n_train=_tr1-_tr0-n_channels;//20131101
	  n_total=n_train+_tp1-_tp0;
	  train_steps=n_train+_tr0;//train_steps=tr1-n_channels;
	  if(_tp1>_tr1) ntotal=_tp1; else ntotal=_tr1;//??
	  if(ntotal>=total_steps) total_steps=ntotal;//061109
	  givendata->tr0=_tr0;
	  givendata->tr1=_tr1;
	  givendata->tp0=_tp0;
	  givendata->tp1=_tp1;
	  givendata->tpD=_tpD;//delay for y
	  givendata->tpG=_tpG;//use x with given data for prediction (for ability check) or not (for prediction)
	  net->tpEy=_tpEy;//delay for y

	  create_data1(data_class,n_channels,n_train,n_total,givendata,t0);
	  create_data1(data_class,n_channels,n_train,n_total,test,t0);

	  givendata->block_begin = train_steps;//ueno(Tue Feb 10 19:36:56 2004)
	  givendata->block_end   = ntotal;
	  givendata->application = application;//ueno(Mon Feb 16 23:27:06 2004)

	  i=0;fp=open_file(fname[i], "r");
	  double *y00=(double *) malloc(sizeof(double)*total_steps);
	  for(t=0;t<total_steps;){
	    fgets(line, LINESIZE, fp);
	    if(feof(fp)) break;
	    if(line[0] == '#') continue;//20140303
	    int nr=sscanf(line, "%lf", &y00[t]);
	    if(nr<1) continue;
	    t++;
	  }
	  fclose(fp);
	  //	  for(t=n_channels;t<n_train;t++){
	  for(t1=0;t1<n_train;t1++){//training
	    t=t1+n_channels;
	    givendata->y[t1]=y00[t+_tr0+_tpD];//??20131031
	    givendata->x[t1][n_channels]=BIAS;
	    for (k=0; k<n_channels; k++) {
	      givendata->x[t1][k] = y00[t-k-1+_tr0];//
	    }
	    //	    if(t1<10) fprintf(stderr,"hen [%d]x%f, y%f\n",t1,givendata->x[t1][0],givendata->y0[t1]);
	  }
	  //	  for(;t<n_total;t++){
	  for(;t1<n_total;t1++){//prediction
	    t=t1-n_train;//t=t1+n_channels;
	    givendata->y[t1]=y00[t+_tp0];
	    givendata->x[t1][n_channels]=BIAS;
	    for (k=0; k<n_channels; k++) {
	      givendata->x[t1][k] = y00[t-k-1+_tp0-_tpD];//?20131031
	    }
	  }
	  free(y00);
	}
      }
    }

    //    total_steps-=n_channels;
    //    givendata_steps-=n_channels;
    //    free(tdata);
    break;
  case FUNCTION_APPROX: //関数近似
    net->printf("Training   Data file name:");fname[0][0]=0;scanf1("%s",fname[0]);//train
    strcpy(givendata->path, fname[0]);
    if(application==RANGE_APPLI) load_data_RANGEDATA(givendata,test,net,&n_train,&n_total,n_channels);
    else {
      net->printf("Prediction Data file name:");fname[1][0]=0;scanf1("%s",fname[1]);//test
      
      t=0; 
      for(i=0;i<2;i++){
	fp = open_file(fname[i], "r");
	while(!feof(fp)){
	  line[0]=0;
	  fgets(line, LINESIZE, fp);
	  if(feof(fp)) break;
	  if(line[0] == 0) break;
	  if (!isnump(line[0])){
	    //	  printf("line=%s.",line);
	    continue;	//if (line[0] == '#') continue;
	  }
	  else{
	    int numexist=1;
	    char *p=&line[0];
	    for(;;){
	      if(*p=='.'||*p=='-'||*p=='+'||*p==' ') {p++;continue;}
	      if((*p >= '0' && *p <= '9')) {break;}
	      else {numexist=0;break;}
	    }
	    if(numexist==0) continue;
	  }
	  ++t;
	  //	  printf("t=%d %s",t,line);
	}
	if(i==0) train_steps = t;
	fclose(fp);
      }
      total_steps = t;
      t0=0;//??
      n_total=total_steps-t0;
      n_train=train_steps-t0;
      net->printf(strx2("Number of training and test data=%d and %d.\n",n_train, n_total-n_train));
      
      create_data1(data_class,n_channels,n_train,n_total,givendata,t0);
      create_data1(data_class,n_channels,n_train,n_total,test,t0);
      t=0;
      for(i=0;i<2;i++){
	fp = open_file(fname[i], "r");
	for(;;){
	  line[0]=0;
	  fgets(line, LINESIZE, fp);
	  if(line[0] == 0) break;
	  if(feof(fp)) break;
	  if (!isnump(line[0])) continue;//	if(line[0] == '#') continue;
	  else{
	    int numexist=1;
	    char *p=&line[0];
	    for(;;){
	      if(*p=='.'||*p=='-'||*p=='+'||*p==' ') {p++;continue;}
	      if((*p >= '0' && *p <= '9')) {break;}
	      else {numexist=0;break;}
	    }
	    if(numexist==0) continue;
	  }

	  givendata->x[t][n_channels]=BIAS;
	  p=line;
	  while(isspace(*p)) p++;
	  for (k=0; k<n_channels; k++){
	    sscanf(p,"%lf",&(givendata->x[t][k]));
	    while(!isspace(*p)) p++;
	    while( isspace(*p)) p++;
	  }
	  sscanf(p,"%lf",&(givendata->y[t]));
	  //if(givendata->y[t] < givendata->min) givendata->min = givendata->y[t];
	  //if(givendata->y[t] > givendata->max) givendata->max = givendata->y[t];
	  if(++t>=total_steps) break;
	}
	fclose(fp);
      }
    }
  }
  //copy x and y to X and Y for Original Scale 
  for(t=0;t<n_total+t0;t++) givendata->Y0[t]=givendata->y0[t];
  for(t=0;t<n_total;t++)for(k=0;k<=n_channels;k++)
    givendata->X[t][k]=givendata->x[t][k];
  //
  calc_statistics(givendata,net);
  normalize_data(givendata,net);
//#define NoAffine
//#ifdef NoAffine
//  for(t=0;t<n_total;t++) givendata->x[t][n_channels]=0.0;
//#endif
  if(application==IJCNN04_APPLI){
    int t_t1,t_t2;
    t_t1=givendata->ijcnn04data->t_t1;
    t_t2=t_t1+19;
    //  k2=0; //k2 is for ARMA model but has not implemented yet
    n_train=(980-n_channels)*4+(t_t1-1-n_channels);
    n_total=n_train+(t_t2-t_t1+1);//
    givendata->n_total=test->n_total=n_total;
    givendata->n_test =test->n_test =n_total-n_train;  
  }

//  net->ymax=givendata->ymax;
//  net->ymin=givendata->ymin;
//  net->xmin=givendata->xmin;
//  net->xmax=givendata->xmax;
  net->ywidth=net->ymax-net->ymin;
  //  net->xwidth=net->xmax-net->xmin;
////  net->xwidth=1.;
////  for(k=0;k<n_channels;k++){
////	net->xwidth *= (net->xmax[k]-net->xmin[k]);
////  }
////  net->xwidth = pow(net->xwidth,1./k);
  net->xwidth=(net->xmax[0]-net->xmin[0]);
  for(k=1;k<n_channels;k++){
    if((net->xmax[k]-net->xmin[k])>net->xwidth) net->xwidth=(net->xmax[k]-net->xmin[k]);
  }

  net->printf(strx4("***ymin,ymax,ymin0,ymax0=%f,%f,%f,%f\n",net->ymin,net->ymax,net->ymin0,net->ymax0));
}
//#define SenkeiHokan(x,x0,y0,x1,y1) ((FLOAT)((x)-(x0))*((y1)-(y0))/((x1)-(x0))+(y0))

/* this file ends here,,, */
