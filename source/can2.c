/*
 * $id: time-stamp: "6��  4, 2003 23:52:36" / takamasa ueno exp.$
 *
 * main.c
 *
 * 2003.05.27
 * ��˥åȤν��������ˡ���Ѥ��뤳�Ȥǡ�
 * ��store_vector() �� if (net->cell[s0].v < v_thresh) ��
 *   modify_x() �� if (net->cell[m].v >= v_thresh) ��פ��ڤäƤϤ�������
 * ��_i�ʥ�˥åȤβ��͡ˤ��ͤ�����ˤʤ뤳�Ȥ�������
 * v_ratio = 1.1 ������ˤ���ȡ������⥤���ؽ���̤�����줿��
 *
 * sim.c �� exec_pred() �� ���ϥǡ�����񤭴�����Ȥ����Х��ե��å���
 * ��pred->x[t+1][0] = 1.0 �˽�����
 *
 * 2003.05.25
 * ��ʬ�ǽ񤤤� my_plinn.c ���ɡ����Ƥ⤦�ޤ������ʤ��Τǡ�
 * ������ plinn2.c ���ѿ�����ʬ�����Ѥ��Ƥ��Τޤ�ή�Ѥ��뤳�Ȥˡ�
 *
 * 2003.05.15
 * ��_i�ʥ�˥åȤβ��͡ˤη׻������餫�ˤ�������
 * ��ľ��ɬ�פ����롩
 *
 * 2003.05.07
 * �����٤줬ȯ�����Ƥ������ǡ������åȤ���Ȥ��ѹ����뤳�Ȥǲ���
 *
 * 2003.05.01
 * �Ȥꤢ�����ץ���ब�������뤬�����ʤ�Х����ʴ���
 * ����饤��ؽ���ͽ¬��̤˻����٤줬�����롥
 *
 * 2003.04.06
 * ����Ϣ�ۥͥå� Competitive Associative Networks (CAN2)
 * ������ plinn2.c ��Хå��ؽ������ˤ��б��Ǥ���褦
 * �����ǥ��󥰤�ľ��ͽ��������Ĺ������
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
 * �ǡ����ʹ�¤�ΡˤΥݥ��󥿤�����å�����
 * ���������Ƥ��ʤ��ʤ�NULL��ɽ������
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
 * ���ƤΥǡ����ʹ�¤�Ρˤ��˴�����
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
 * �ץ�����λ����
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
 * �¹Ԥ��륳�ޥ�ɤ�����
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


// ���ƤΥǡ��������������Ƥ��뤫��
//static int b_init_all_data = 0;

/**********************************************************************
 * main
 * �ᥤ��롼����
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
    //�����ǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
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
//    // �ͥåȤ���¸
//    else if (strncmp(buff, "ns", 2) == 0) {
//      net_save(net,NULL);
//    }
//    // �ͥå��ɹ�
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
////    //�����ǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
////      //single step prediction from beggining
////      PREDTRAIN=1;
////      exec_ssp_train(net, givendata, test);
////      PREDTRAIN=0;
////    }
//    else if (strncmp(buff, "ssprt", 5) == 0) {//for rangedata
//      //�ƥ��ȥǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
//      PREDTRAIN=0;
//      exec_ssp_test_rt(net, givendata, test);
//      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//    }
////    else if (strncmp(buff, "sspE", 4) == 0) {
////    //�����ǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
////      //single step prediction from beggining
////      PREDTRAIN=1;
////      exec_ssp_train_E(net, givendata, test);
////      //      exec_plot(givendata, test, "SSPofTrainData",SSP,net);
////      PREDTRAIN=0;
////    }
//    else if (strncmp(buff, "sspn", 4) == 0) {
//    //�ƥ��ȥǡ�����ͽ¬ 
//      exec_ssp_test_NIPS04(net, givendata, test);
//    }
//    else if (strncmp(buff, "sspr", 4) == 0) {//for range data
//      //�ƥ��ȥǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
//      PREDTRAIN=0;
//      exec_ssp_test_r(net, givendata, test);
//      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "sspe", 4) == 0) {
//      //�ƥ��ȥǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
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
////      //�ƥ��ȥǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
////      //single step prediction from beggining
////      PREDTRAIN=0;
////      exec_ssp_test(net, givendata, test);
////      //      exec_plot(givendata, test, "SSPofTestData",SSP,net);
////      PREDTRAIN=0;
////      fprintf(stderr,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
////    }
//    else if (strncmp(buff, "ssp", 3) == 0) {
//      //�ƥ��ȥǡ����ΰ���ͽ¬(single-step prediction)��ؽ�ǽ��
//      //single step prediction from beggining
//      PREDTRAIN=0;
//      exec_ssp_test(net, givendata, test);
//      if(buff[3]!='_') exec_plot(givendata, test, "SSPofTestData",SSP,net);
//      PREDTRAIN=0;
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
//    else if (strncmp(buff, "msp0", 4) == 0) {
//    //�����ǡ�����¿��ͽ¬(multistep prediction)��ؽ�ͽ¬ǽ��
//      PREDTRAIN=2;
//      exec_msp_train(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTrainData",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "msp2", 4) == 0) {
//    //�����ǡ�����¿��ͽ¬(multistep prediction)��ؽ�ͽ¬ǽ��
//      PREDTRAIN=3;
//      exec_msp_traintest(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTrainTestData",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "msp1", 4) == 0) {
//    //�����ǡ�����¿��ͽ¬��̤�ƥ��ȥǡ�����ͽ¬������
//      PREDTRAIN=0;
//      exec_msp_test1(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData1",MSP,net);
//      PREDTRAIN=0;
//    }
//    else if (strncmp(buff, "mspj", 4) == 0) {
//    //�ƥ��ȥǡ�����¿��ͽ¬(multistep prediction)
//      exec_msp_test_IJCNN04_out(net, givendata, test);
//    }
////    else if (strncmp(buff, "mspr", 4) == 0) {
////	//�ƥ��ȥǡ�����¿��ͽ¬(multistep prediction)
////	exec_msp_test_real(net, givendata, test);
////	exec_plot(givendata, test, "MSPofTestData",MSP,net);
////	printf(">>MSE---;%e NMSE---;%e\n",test0->MSE,test0->NMSE);
////    }
//    else if (strncmp(buff, "mspe", 4) == 0) {
//    //�ƥ��ȥǡ�����¿��ͽ¬(multistep prediction)
//      exec_msp_test_ensemble(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData",MSP,net);
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
//    else if (strncmp(buff, "mspE", 4) == 0) {
//    //�ƥ��ȥǡ�����¿��ͽ¬(multistep prediction)
//      exec_msp_test_Ensemble(net, givendata, test);
//      exec_plot(givendata, test, "MSPofTestData",MSP,net);
//      fprintf(stdout,">>MSE---;%e NMSE---;%e\n",test->MSE,test->NMSE);
//    }
    else if (strncmp(buff, "msp", 3) == 0) {
    //�ƥ��ȥǡ�����¿��ͽ¬(multistep prediction)
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
//    // ��¤�Υǡ����Υݥ��󥿤�����å�
//    else if (strncmp(buff, "pt", 2) == 0) {
//      //      show_all_pointer(net, tdata, givendata, aprx, test, test0);
//    }
//    // ���ͤΥǡ������åȤγ�ǧ
//    else if (strncmp(buff, "r0", 2) == 0){
//      net->r1=net->r2=0;
//      net->r3=0;
//    }
//    else if (strncmp(buff, "sr", 2) == 0) show_data_parms(givendata);
//    // �ؽ��ǡ������åȤγ�ǧ
//    //    else if (strncmp(buff, "sa", 2) == 0) show_data_parms(givendata0);
//    // ����ͽ¬�ǡ������åȤγ�ǧ
//    else if (strncmp(buff, "st", 2) == 0) show_data_parms(test);
//    // ¿��ͽ¬�ǡ������åȤγ�ǧ
//    //    else if (strncmp(buff, "sp", 2) == 0) show_data_parms(test0);
//    // �Хå����ؽ��ѤΥǡ������åȤγ�ǧ
//    //else if (strncmp(buff, "sb", 2) == 0) show_batch_data_parms(bdata);
//    // �ͥåȤΥѥ�᡼�����ǧ
//    else if (strncmp(buff, "sn", 2) == 0) show_net_parms(net);
//    // �Ťߤ�ɽ��
//    else if (strncmp(buff, "spr:", 4) == 0){//��Υ�ǡ������ʿ�̤�õ��
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
//    //    else if (strncmp(buff, "sp", 2) == 0) search_planes(net);    // ʿ�̤�õ��
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
    // ���顼����
    else continue;

    //    net->printf("\n\n");
    //    printf("\n");
  }
  return 0;
}


/* this file ends here,,, */
