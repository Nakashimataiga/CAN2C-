//gcc -g pred2y.c -o pred2y
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define isnum(c) (((c>='0') && (c<='9'))||(c=='-')||(c=='+')||(c==' '))
int main(int argc, char **argv)
{
  //  int i,j;
  FILE *fp1,*fp2,*fp3;
  char buff1[1024],buff2[1024];
  double x1=0,y,yp,x1b,mse=0,num=0,err;
  double psum=0,p;

  //  printf("\n# %0x,%0x,%0x\n",'0','1','9');
  if(argc<3){
    fprintf(stderr,"Usage: %s <predfile> <truefile> \n",argv[0]);
    exit(1);
  }

  if((fp1=fopen(argv[1],"rt"))==NULL){
    fprintf(stderr,"File(%s) Open Eoor\n",argv[1]);
    fclose(fp1);
    return(-1);
  }
  if((fp2=fopen(argv[2],"rt"))==NULL){
    fprintf(stderr,"File(%s) Open Eoor\n",argv[2]);
    fclose(fp1);
    return(-1);
  }
  if(argc>=4){
    if((fp3=fopen(argv[3],"rt"))==NULL){
      fprintf(stderr,"File(%s) Open Eoor\n",argv[2]);
      fclose(fp1);
      return(-1);
    }
  }
  psum=0;
  //  printf("# yp,y,err,mse\n");
  for(num=0;;){
    for(;;){
      fgets(buff1,1024,fp1);
      //      printf("buff1:%s\n",buff1);
      if(feof(fp1)) break;
      if(isnum(buff1[0])) break;
    }
    //    printf("fp1=%0x,fp2=%0x\n",fp1,fp2);
    for(;;){
      fgets(buff2,1024,fp2);
      //      printf("buff2:%s\n",buff2);
      if(feof(fp2)) break;
      if(isnum(buff2[0])) break;
    }
    if(feof(fp1)) break;
    if(feof(fp2)) break;
    sscanf(buff1,"%lf",&yp);
    sscanf(buff2,"%lf",&y);
    p=1;
    psum+=p;
    num++;
    err=yp-y;
    mse +=(err*err)*p;
    if(fabs(x1b-x1)>1e-10) printf("\n");
    x1b=x1;
    fprintf(stdout,"%15.7e %15.7e %+15.7e %+15.7e #yp,y,err,mse\n",yp,y,err,mse/num);
  }
  fclose(fp1);
  fclose(fp2);
  fprintf(stdout,"#MSE=%e for %s (num=%d)\n",mse/psum,argv[2],(int)num);
  fprintf(stdout,"#MSE=%f (num=%d) between %s and %s.\n",mse/psum,(int)num,argv[1],argv[2]);
  //  fprintf(stderr,"# from %s and %s \n",argv[1],argv[2]);
  //  fprintf(stderr,">>MSE=%e is of %s (num=%d)\n",mse/psum,argv[2],(int)num);
  //  printf("#MSE=%e=%f,num=%10.0f ",mse/psum,mse/psum,num);
  //  fprintf(stderr,">>MSE=%e=%f,num=%10.0f ",mse/psum,mse/psum,num);
  fp1=fopen("mse.dat","w");
  fprintf(fp1,"%13.7e %d\n",mse/num,(int)num);
  fclose(fp1);
  return(1);
}
