/*
 * $id: time-stamp: "5月 30, 2003 18:42:35" / takamasa ueno exp.$
 *
 * my_misc.c
 *
 */
#include "my_misc.h"
//double square(double x){return(x*x);}
/*====================================================================*
 * open_file
 * ファイルを開きファイルポインタを返す
 *====================================================================*/
FILE* open_file(char *fname, char *mode) {
  char name[32] = "open_file";
  FILE *fp = NULL;

  if ((fp=fopen(fname, mode)) == NULL) {
    printf("error: cannot open file \"%s\". (%s, mode:%s)\n", fname, name, mode);
    //    exit(-1);
    return(NULL);
  }
  return(fp);
}

/*====================================================================*
 * close_file
 * ファイルを閉じファイルポインタを空にする
 *====================================================================*/
void close_file(FILE *fp) {
  fclose(fp);
  return;
}

/*====================================================================*
 * open_pipe
 * パイプを開きファイルポインタを返す
 *====================================================================*/
FILE *open_pipe(char *cmd, char *mode) {
  char name[32] = "open_pipe";
  FILE *fp = NULL;

  if ((fp=popen(cmd, mode)) == NULL) {
    printf("error: oops, cannot find \"%s\". (%s, mode:%s)\n", cmd, name, mode);
    exit(-1);
  }
  return(fp);
}

/*====================================================================*
 * close_pipe
 * パイプを閉じファイルポインタを空にする
 *====================================================================*/
void close_pipe(FILE *fp) {
  pclose(fp);
  return;
}

/*====================================================================*
 * count_file
 * ファイルの行数を数える
 *====================================================================*/
int count_file(FILE *fp) {
  int count = 0;
  char line[512];

  while (!feof(fp)) {
    fgets(line, 512, fp);
    ++count;
    //printf(line); // for Debug
  }
  return(count);
}

/*====================================================================*
 * calc_length
 * x[dim]とy[dim]との距離を計算する
 *====================================================================*/
FLOAT calc_length(FLOAT *x, FLOAT *y, int dim) {
  FLOAT length = 0.0;
  int i;

  for (i=0; i<dim; i++) length += (FLOAT)square(x[i]-y[i]);
  length = (FLOAT)sqrt(length);
  return(length);
}

FLOAT distance2(FLOAT *x, FLOAT *y, int dim) {
  FLOAT d2 = 0.0;
  int i;
  for (i=0; i<dim; i++) d2 += (x[i]-y[i])*(x[i]-y[i]);

  //  for (i=0; i<dim; i++) d2 += (dd=(x[i]-y[i]))*dd;
  return(d2);
}

//static int my_malloc_total=0;
void *my_malloc(int size,char *mes,int mes_onoff)
{
  char *ptr;
  if((ptr=(char *)malloc(size)) == 0){
    fprintf(stderr,"Could not allocate '%s' (%d bytes)\n",mes,size);
    fprintf(stderr,"malloc_total=%d at %s[%d]\n",my_malloc_total,mes,size);
  }
  else{
    my_malloc_total+=size;
    if(mes_onoff==1){
      printf("malloc_total=%d at %s[%d]\n",my_malloc_total,mes,size);
    }
  }
  return(ptr);
}
void my_free(void *ptr)
{
  my_malloc_total-=sizeof(ptr);
  free(ptr);
}
//struct timeval mytv[3],*mytv0,*mytv1,*mytv2,*mytv3;
double timeval_subtract (struct timeval *x, struct timeval *y) /* from gnu c-lib info.texi */
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  return((x->tv_sec - y->tv_sec)+(x->tv_usec - y->tv_usec)/1000000.);
  
  /* Compute the time remaining to wait.
     `tv_usec' is certainly positive. */
  //  result->tv_sec = x->tv_sec - y->tv_sec;
  //  result->tv_usec = x->tv_usec - y->tv_usec;
  
  /* Return 1 if result is negative. */
  //  return x->tv_sec < y->tv_sec;
}
void mytimer_start()
{
  //  mytv=(timeval *)malloc(sizeof(timeval)*3);
  mytv0=&mytv[0];
  mytv1=&mytv[1];
  mytv2=&mytv[2];
  gettimeofday(mytv0,NULL);
  gettimeofday(mytv2,NULL);
}
double mytimer_lap()
{
  mytv3=mytv1;
  mytv1=mytv2;
  mytv2=mytv3;
  gettimeofday(mytv2,NULL);
  return(timeval_subtract(mytv2,mytv1));
}
double mytimer_total()
{
  gettimeofday(mytv2,NULL);
  return(timeval_subtract(mytv2,mytv0));
}
char *fnbody(char *p,char *pr)
{
  int i,i0,i1,len=strlen(p);
  i=len-1;
  i1=len;
  i0=0;
  for(i=strlen(p);i>=0;i--){
    if(i1==len && p[i]=='.') i1=i;
    if(p[i]=='/') {
      i0=i+1;
      break;
    }
  }
  for(i=0;i<i1-i0;i++) pr[i]=p[i+i0];
  pr[i1-i0]=0;
  return(pr);
}

int printf1(char *buf) {printf("%s",buf);return 0;}
int noprintf(char *buf) {return 0;}
//gettimeofday(&tv1, NULL);
//timeval_subtract(&tv10, &tv1, &tv0);

/* this file ends here,,, */
