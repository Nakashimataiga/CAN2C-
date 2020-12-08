/*
 * $id: time-stamp: "6月  4, 2003 23:52:36" / takamasa ueno exp.$
 *
 * main.c
 *
 * 2003.05.27
 * ユニットの初期化の方法を変えることで，
 * （store_vector() の if (net->cell[s0].v < v_thresh) と
 *   modify_x() の if (net->cell[m].v >= v_thresh) を思い切ってはずす！）
 * α_i（ユニットの価値）の値が負数になることを一応回避．
 * v_ratio = 1.1 あたりにすると，前よりもイイ学習結果が得られた．
 *
 * sim.c で exec_pred() の 入力データを書き換えるところをバグフィックス
 * （pred->x[t+1][0] = 1.0 に修正）
 *
 * 2003.05.25
 * 自分で書いた my_plinn.c がどーしてもうまくいかないので，
 * 先生の plinn2.c を変数の部分だけ変えてそのまま流用することに．
 *
 * 2003.05.15
 * α_i（ユニットの価値）の計算が明らかにおかしい
 * 書き直す必要がある？
 *
 * 2003.05.07
 * 時間遅れが発生していたがデータセットの中身を変更することで回避．
 *
 * 2003.05.01
 * とりあえずプログラムが完成するが，かなりバギー（汗）
 * オンライン学習の予測結果に時間遅れが生じる．
 *
 * 2003.04.06
 * 競合連想ネット Competitive Associative Networks (CAN2)
 * 先生の plinn2.c をバッチ学習処理にも対応できるよう
 * コーディングし直す予定だが先は長いぞ．
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "my_function.h"
#include "my_plinn.h"
#include "sim.h"
#include "my_misc.h"
#include "am.h"
#include "randoms.h"

/**********************************************************************
 * show_all_pointer
 * データ（構造体）のポインタをチェックする
 * 初期化されていないならNULLを表示する
 **********************************************************************/
//void show_all_pointer(NET *net, TDATA *tdata, DATA *train,
//			DATA *aprx, DATA *test, DATA *pred) {
//  printf("================================================================\n");
//  printf("	Show The Pointer Of Structure Data\n");
//  printf("================================================================\n");
//  printf("> net         : %p\n", net);
//    printf("> tdata       : %p\n", tdata);
//  printf("> train       : %p\n", train);
//  printf("> aprx        : %p\n", aprx);
//  printf("> test        : %p\n", test);
//  printf("> pred        : %p\n", pred);
//  return;
//}

/**********************************************************************
 * remove_all_data
 * 全てのデータ（構造体）を破棄する
 **********************************************************************/
//void remove_all_data(NET *net, TDATA *tdata,
//		       DATA *train, DATA *aprx, DATA *test, DATA *pred) {
//  //remove_net(net);
//  //  remove_time_data(tdata);
//  remove_data(train);
//  remove_data(aprx);
//  remove_data(test);
//  remove_data(pred);
//  return;
//}

/**********************************************************************
 * quit_program
 * プログラムを終了する
 **********************************************************************/
//void quit_program(NET *net, TDATA *tdata,
//		  DATA *train, DATA *aprx, DATA *test, DATA *pred) {
int quit_program(NET *net)
{
  //remove_all_data(net, tdata, train, aprx, test, pred);
  //  printf("\n");
  net->printf("Good Bye ... (\"-\")zzz\n");
  exit(0);
  return(0);
}

/**********************************************************************
 * get_command
 * 実行するコマンドを得る
 **********************************************************************/
void help() {
  printf("================================================================\n");
  printf("	Command\n");
  printf("================================================================\n");
  printf("	ex: execute the simulation\n");
  printf("	ri: reinitialize each data\n");
  printf("	sn: show the network parameters\n");
  printf("	sr: show the traininig data sets\n");
  printf("	sa: show the approximated data sets\n");
  printf("	st: show the test data sets\n");
  printf("	sp: show the predict data sets\n");
  printf("	sb: show the data for batch learning\n");
  printf("	pt: show the pointer of structure data\n");
  printf("     msp: multi-step prediction of test data\n");
  printf("    msp0: multi-step prediction of training data\n");
  printf("    ssp0: single-step prediction of training data\n");
  printf("	he: see this list\n");
  printf("	qu: quit the program\n");
  printf("----------------------------------------------------------------\n");
  return;
}

void get_command(char *buff,NET *net) {
  //  if(net->nop!=1) net->printf("* Command? : "); 
  buff[0]=0;
  //  scanf1("%s",buff);
  fgets(buff,256,stdin);
  return;
}


// 全てのデータが初期化されているか？
//static int b_init_all_data = 0;

/**********************************************************************
 * main
 * メインルーチン
 **********************************************************************/
void normalize_data(DATA *givendata,NET *net);
void calc_statistics(DATA *givendata, NET *net);

NET *net;

int main(int argc,char *argv[]) {
  //////////////////////////////////////////
  // Initialize memories and parameters
  //////////////////////////////////////////
  DATA *givendata= NULL, *test = NULL;//, *test0 = NULL;
  char buff[256];
  net=(NET *)malloc(sizeof(NET)*1);
  givendata=(DATA *)malloc(sizeof(DATA)*1);//
  test =(DATA *)malloc(sizeof(DATA)*1);
  PREDTRAIN=0;
  BIAS=1;
  LINESIZE=10240;
  net->r1=0; net->r3=0;//
  net->init=net->winit=net->Vinit=0;
  net->seed=1;
  net->Tpinv=999;
  char netmsg[999]={0,};
  net->msg=(char *)&netmsg;
  net->DISP=1;
  net->nop=0;
  //  int nop=0;
  int i;
  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"BIAS:",5)==0) sscanf(&argv[i][5],"%lf",&BIAS);
    else if(strncmp(argv[i],"LINESIZE:",9)==0) sscanf(&argv[i][9],"%d",&LINESIZE);
    else if(strncmp(argv[i],"seed:",5)==0) sscanf(&argv[i][5],"%lu",&net->seed);
    else if(strncmp(argv[i],"Tpinv:",6)==0) sscanf(&argv[i][6],"%d",&net->Tpinv);
    else if(strncmp(argv[i],"DISP:",5)==0) sscanf(&argv[i][5],"%d",&net->DISP);
    else if(strncmp(argv[i],"nop:",4)==0) sscanf(&argv[i][4],"%d",&net->nop);
  }
  if(net->nop==0) net->printf=printf1;
  else net->printf=noprintf;
  //////////////////////////////////////////
  if(net->nop==0){
    printf("Called with RAND=%d:",RAND);
    for(i=0;i<argc;i++){
      printf("%s ",argv[i]);
    }
    printf("\n");
  }
  //////////////////////////////////////////
  // set seed of random function
  //////////////////////////////////////////
#if RAND == RANDOM
  srandom((unsigned int)net->seed);
#elif RAND == DRAND48
  seed48((unsigned short*)&net->seed);
#elif RAND == MYRAND
  _randn=(unsigned long)net->seed;
#elif RAND == ZMTRAND
  InitMt((unsigned long)net->seed);
#endif

  //////////////////////////////////////////
  // Load data (training and test)
  //////////////////////////////////////////
  load_data(givendata,test,net);//in my_function.c

  //  if(nop==0) test->printf=givendata->printf=net->printf=printf1;
  //else test->printf=givendata->printf=net->printf=noprintf;

  //////////////////////////////////////////
  // Initialize time
  //////////////////////////////////////////
  ReinitTime=GlobalTime=GlobalTimeMax=0;
  SuccessiveLearn=0;
  //  help();
  /* example of param.dat
1          #0:time-series,1:function,3:ijcnn04,4:range-data,No.1
1 0       #dimensionality k1 k2					
/home/kuro/sotu/2019/chainer/can2py/tmp/train.csv						       	        
/home/kuro/sotu/2019/chainer/can2py/tmp/test.csv				        
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00:0.000000e+00:0.000000e+00:0.000000e+00   #y-normalization[ymin0:ymax0]->[ymin1:ymax1]
in          init-net			                
90          number-of-cells				          	
6           n_compare					
2.000000e-01 3 0     v_thresh vmin vmin2	         		
5.000000e+00           v_ratio				                
2.000000e-01 5.000000e-02 7.000000e-01   width,gamma  #window width                         
ex          execution						
1 5.000000e-02 7.000000e-01   #<0:online,1:batch>, <gamma0>, <entropy_thresh>     
100        i_times number of learning iterations	        
100 50 350    number-of-display, rot_x, rot_z			       
ssp0_      #prediction of training data -> predict0.dat    
qu 
   */
  while (1) {
    get_command(buff,net); 
    //    net->printf("\n");
    if (strncmp(buff, "qu", 2) == 0) {
      return(0);
      //      quit_program(net);
    }
    //////////////////////////////////////////
    // Initialize net
    //////////////////////////////////////////
    else if (strncmp(buff, "in", 2) == 0) {
      net=init_net(net);//should be called after load_data n_channels
    }
    //////////////////////////////////////////
    // Execute Simulation
    //////////////////////////////////////////
    else if (strncmp(buff, "ex", 2) == 0) {
      //      net->printf=printf1;
      exec_sim(net, givendata, test);//strange message20191101
    }
    else if (strncmp(buff, "ssp0", 4) == 0) {
    //訓練データの一段予測(single-step prediction)＝学習能力
      //single step prediction from beggining
      PREDTRAIN=1;
      exec_ssp_train(net, givendata, test);//for ssp0_, ssp00 <- same??
      if(!(buff[4]=='_' || buff[4]=='0')) exec_plot(givendata, test, "SSPofTrainData",SSP,net);
      PREDTRAIN=0;
      printf("%s\n",net->msg);
    }
//    //
//    else if (strncmp(buff, "he", 2) == 0) help();
    else if (buff[0]=='!') system(&buff[1]);
//    else if (strncmp(buff, "in", 2) == 0) {
//      net=init_net(net);//should be called after load_data n_channels
//    }
//    // ネットの保存
//    else if (strncmp(buff, "ns", 2) == 0) {
//      net_save(net,NULL);
//    }
//    // ネット読込
    else if (strncmp(buff, "sM", 2) == 0) {
      //      printf("\n*** Saving nets with N%d b%d***\n",net->n_cells,net->nEns);//check20191114
      save_M(net);
    }
    else if (strncmp(buff, "nls", 3) == 0) {//for onsei
      net=net_loads(net);
      if(net->nop==0) net->printf=printf1;
      else net->printf=noprintf;
      printf("\n***Loaded nets with N%d b%d***\n",net->n_cells,net->nEns);//check20191114
      //      net->printf("##b=%d\n",net->nEns);
    }
    else if (strncmp(buff, "nl", 2) == 0) {
      net_load(net,NULL);
    }
//    else if (strncmp(buff, "swm", 3) == 0) {
//      save_wm(net);
//    }
//    else if (strncmp(buff, "swM", 3) == 0) {
//      save_wM(net);
//    }
//    else if (strncmp(buff, "pMs", 3) == 0) {//for onsei
//      poles_of_M_shrink(net);
//    }
//    else if (strncmp(buff, "ld", 2) == 0) {
//      load_data(givendata,test,net);
//      //      load_data(givendata,test,test0,net);
//    }
//    else if (strncmp(buff, "no", 2) == 0) {
//      normalize_data(givendata,net);
//      //      test0->VAR=givendata->VARtest;
//      //      test->VAR=givendata->VARtest;
//    }
////    else if (strncmp(buff, "ssp00", 5) == 0) {
////    //訓練データの一段予測(single-step prediction)＝学習能力
////      //single step prediction from beggining
////      PREDTRAIN=1;
////      exec_ssp_train(net, givendata, test);
////      PREDTRAIN=0;
////    }
//    else if (strncmp(buff, "ssprt", 5) == 0) {//for rangedata
//      //テストデータの一段予測(single-step prediction)＝学習能力
//      PREDTRAIN=0;
//      exec_ssp_test_rt(net, givendata, test);
//      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//    }
////    else if (strncmp(buff, "sspE", 4) == 0) {
////    //訓練データの一段予測(single-step prediction)＝学習能力
////      //single step prediction from beggining
////      PREDTRAIN=1;
////      exec_ssp_train_E(net, givendata, test);
////      //      exec_plot(givendata, test, "SSPofTrainData",SSP,net);
////      PREDTRAIN=0;
////    }
//    else if (strncmp(buff, "sspn", 4) == 0) {
//    //テストデータの予測 
//      exec_ssp_test_NIPS04(net, givendata, test);
//    }
//    else if (strncmp(buff, "sspr", 4) == 0) {//for range data
//      //テストデータの一段予測(single-step prediction)＝学習能力
//      PREDTRAIN=0;
//      exec_ssp_test_r(net, givendata, test);
//      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "sspe", 4) == 0) {
//      //テストデータの一段予測(single-step prediction)＝学習能力
//      //single step prediction from beggining
//      PREDTRAIN=0;
//      exec_ssp_test_ensemble(net, givendata, test);
//#if GPLT >= 1
//      exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//#endif
//    }
////    else if (strncmp(buff, "ssp_", 4) == 0) {//without display
////      //テストデータの一段予測(single-step prediction)＝学習能力
////      //single step prediction from beggining
////      PREDTRAIN=0;
////      exec_ssp_test(net, givendata, test);
////      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
////      PREDTRAIN=0;
////      fprintf(stderr,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
////    }
//    else if (strncmp(buff, "ssp", 3) == 0) {
//      //テストデータの一段予測(single-step prediction)＝学習能力
//      //single step prediction from beggining
//      PREDTRAIN=0;
//      exec_ssp_test(net, givendata, test);
//      if(buff[3]!='_') exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
//    else if (strncmp(buff, "msp0", 4) == 0) {
//    //訓練データの多段予測(multistep prediction)＝学習予測能力
//      PREDTRAIN=2;
//      exec_msp_train(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTrainData",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "msp2", 4) == 0) {
//    //訓練データの多段予測(multistep prediction)＝学習予測能力
//      PREDTRAIN=3;
//      exec_msp_traintest(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTrainTestData",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "msp1", 4) == 0) {
//    //訓練データの多段予測結果をテストデータの予測に利用
//      PREDTRAIN=0;
//      exec_msp_test1(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData1",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "mspj", 4) == 0) {
//    //テストデータの多段予測(multistep prediction)
//      exec_msp_test_IJCNN04_out(net, givendata, test);
//    }
////    else if (strncmp(buff, "mspr", 4) == 0) {
////	//テストデータの多段予測(multistep prediction)
////	exec_msp_test_real(net, givendata, test);
////	exec_plot(givendata, test, "MSPofTestData",MSP,net);
////	printf(">>MSE---;%e NMSE---;%e\n",test0->MSE,test0->NMSE);
////    }
//    else if (strncmp(buff, "mspe", 4) == 0) {
//    //テストデータの多段予測(multistep prediction)
//      exec_msp_test_ensemble(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData",MSP,net);
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
//    else if (strncmp(buff, "mspE", 4) == 0) {
//    //テストデータの多段予測(multistep prediction)
//      exec_msp_test_Ensemble(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData",MSP,net);
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
    else if (strncmp(buff, "msp", 3) == 0) {
    //テストデータの多段予測(multistep prediction)
      //      int _n_total,n_total;
      //      sscanf(&buff[4],"%d",&_n_total);
      //      n_total=givendata->n_total;
      //      givendata->n_total=_n_total;
      double _et,_ep;
      double err4terminate=0;
      double err4propagate=0;//20150218 
      int nc=sscanf(&buff[4],"%lf%lf",&_et,&_ep);
      if(nc>=2) err4propagate=_ep; 
      if(nc>=1) err4terminate=_et; 
      exec_msp_test(net, givendata, test,err4terminate,err4propagate);
      if(buff[3]!='_') exec_plot(givendata, test, "MSPofTestData",MSP,net);
      
      //      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
      printf("%s\n",net->msg);
      //      givendata->n_total=n_total;
    }
//    // 構造体データのポインタをチェック
//    else if (strncmp(buff, "pt", 2) == 0) {
//      //      show_all_pointer(net, tdata, givendata, aprx, test, test0);
//    }
//    // 真値のデータセットの確認
//    else if (strncmp(buff, "r0", 2) == 0){
//      net->r1=net->r2=0;
//      net->r3=0;
//    }
//    else if (strncmp(buff, "sr", 2) == 0) show_data_parms(givendata);
//    // 学習データセットの確認
//    //    else if (strncmp(buff, "sa", 2) == 0) show_data_parms(givendata0);
//    // 一段予測データセットの確認
//    else if (strncmp(buff, "st", 2) == 0) show_data_parms(test);
//    // 多段予測データセットの確認
//    //    else if (strncmp(buff, "sp", 2) == 0) show_data_parms(test0);
//    // バッチ型学習用のデータセットの確認
//    //else if (strncmp(buff, "sb", 2) == 0) show_batch_data_parms(bdata);
//    // ネットのパラメータを確認
//    else if (strncmp(buff, "sn", 2) == 0) show_net_parms(net);
//    // 重みの表示
//    else if (strncmp(buff, "spr:", 4) == 0){//距離データ中の平面の探索
//      double theta1;
//      int maxnp,maxoptit,DISP=1;
//      int ncol=0;
//      char *p=&buff[4];
//      for(;;p++){
//	if(*p==0) break;
//	if(*p==':') ncol++;
//      }
//      if(ncol>=3) sscanf(&buff[4],"%lf:%d:%d:%d",&theta1,&maxnp,&maxoptit,&DISP);
//      else sscanf(&buff[4],"%lf:%d:%d",&theta1,&maxnp,&maxoptit);
//      search_planes_2d(net,givendata,test,theta1,maxnp,maxoptit,DISP);    
//    }
//    //    else if (strncmp(buff, "sp", 2) == 0) search_planes(net);    // 平面の探索
//    else if (strncmp(buff, "sw", 2) == 0){
//      int nn=10;
//      if(strlen(buff)>3) sscanf(&buff[3],"%d",&nn);
//      show_weights(net,nn);
//    }
//    else if (strncmp(buff, "po", 2) == 0) {
//      pred_out(net, givendata, test);
//    }
//    else if (strncmp(buff, "pMm", 3) == 0) {
//      calc_poles_of_Mmean(net);
//    }
//    else if (strncmp(buff, "pM", 2) == 0) {
//      calc_poles_of_M(net);
//    }
//    else if (strncmp(buff, "wait", 4) == 0) {
//      system("xterm -T \"Waiting CAN2. Hit a Key to Continue.\" -geometry 60x5-0+0 -e /bin/sh -c \"read\";wait");
//      //      system("xterm -T \"Waiting CAN2. Hit a Key to Continue.\" -geometry 40x5-0+0 -e /bin/sh -c \"echo \\\"Hit A Key to continue. \\\"; read\"");
//    }
    // エラー処理
    else continue;

    //    net->printf("\n\n");
    //    printf("\n");
  }
  return 0;
}


/* this file ends here,,, */
