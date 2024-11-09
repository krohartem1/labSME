//diaproc.cpp

#include <stdio.h>
#include <math.h>
#include <stdlib.h>    // atof
#include <conio.h>     // getch
#include <string.h>    // strcpy
#include "defglob.h"
#include "diaproc.h"
#define MALBUF 16
#define MAXBUF 80

// ------------- new --------------------
//IAPROC.OBJ : error LNK2001: unresolved external symbol "void __cdecl clear_buffer(int,char *)" (?clear_buffer@@YAXHPAD@Z)
//Debug/CYLINDER.exe : fatal error LNK1120: 1 unresolved externals
//Error executing link.exe.
//Creating browse info file...

DATA_SECTION::DATA_SECTION(void)
{
  clear_buffer(MAX_SECTION_NAME_LEN,name);
}
Real my_square(Real x){return(x*x);}
Real XHeavisideX(Real x){ return( x>=0.0 ? x : 0.0 ); }
Real Sign(Real x){if(x>0.0) return 1.0; else if(x<0.0) return -1.0; else return 0.0;}
Real Min(Real a, Real b){return(a<b ? a : b); }
Real Min(Real a, Real b, Real c)
{
 Real m12 = a<b ? a : b;
 return (m12<c ? m12 : c);
}

Real Min(Real a, Real b, Real c, Real d)
{
 Real m12 = a<b ? a : b;
 Real m34 = c<d ? c : d;
 return (m12<m34 ? m12 : m34);
}

Real Max(Real a, Real b){return(a>b ? a : b); }
Real Max(Real a, Real b, Real c)
{
 Real m12 = a>b ? a : b;
 return (m12>c ? m12 : c);
}

Real Max(Real a, Real b, Real c, Real d)
{
 Real m12 = a>b ? a : b;
 Real m34 = c>d ? c : d;
 return (m12>m34 ? m12 : m34);
}

int is_numeral(char c)
{
  if(  c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||c=='5'||c=='6'||c=='7'||c=='8'
     ||c=='9'||c=='+'||c=='-'||c=='.'||c=='e'||c=='E') return 1;
  else return 0;
}

void get_number(char* s)
{//procedure truncates string s leaving char-s that suit for a number
  char *p = s;
  while(is_numeral(*p)) p++;
  *p='\0';
  return;
}

void remove_comment(char* s)
{//procedure places '\0' just before symbols // are encountered
  char*p = s;
  char c=' ';
  while( !(*p=='\0' || *p=='/'&& c=='/') ) c = *p++;
  if(*p!='\0') *(p-1)='\0';
  return;
}


//Function Hev = Heaviside(x)
//returns 1.0 if x>0 , 0.0 if x<0, and Hev_at_0 if x==0 (by default 0.0)
double Hev(double x, double Hev_at_0){
  if(x>0.0) return 1.0;
  if(x<0.0) return 0.0;
  return Hev_at_0;
}
//Function Hev2 = Heaviside(x*y)
//returns 1.0 if x and y have the same sign and 0.0 otherwise
double Hev2(double x, double y){
  if(x>0.0 && y>0.0 || x<0.0 && y<0.0) return 1.0;
  return 0.0;
}

//Function My_pow(x,y)
//returns x**y when x>0 and returns 0 when x==0 && y>=0 
double My_pow(double x, double y){
  if(x<0){fprintf(stderr,"\nMy_power: x<0 in x**y for double y"); exit(2);}
  if(x==0){
    if(y<0){fprintf(stderr,"\nMy_power: x==0 && y<0 in x**y"); exit(2);}
    return 0.0;
  }
  double exp_pow = y*log(x);
  if(exp_pow>EXP_POW_MAX){ // define EXP_POW_MAX  709.7
    fprintf(stderr,"\nMy_power: overflow in x**y:");
    fprintf(stderr,"\n x=%le, y=%le", x,y);
    exit(2);
  }
  if(exp_pow<EXP_POW_MIN){ // define EXP_POW_MIN  -708.6
    fprintf(stderr,"\nMy_power: underflow in x**y:");
    fprintf(stderr,"\n x=%le, y=%le", x,y);
    return 0.0;
  }
  return exp(exp_pow);
}

//Function My_pow(x,n)
//returns x**n for integer n and double x
double My_pow(double x, int n){
  if(x==0){
    if(n<0){fprintf(stderr,"\nMy_power: x==0 && n<0 in x**n"); exit(2);}
    else return 0.0;
  }
  double res=1.0;
  if(n>0) while(n-- > 0) res *= x;
  if(n<0) while(n++ < 0) res /= x;
  return res;
}

//Function My_pow(k,n)
//returns k**n for integer k and n>=0
int My_pow(int k, int n){
  if(n<0){fprintf(stderr,"\nInteger version of My_power does not allow n<0 in k**n"); exit(2);}
  if(k==0){
    if(n==0){fprintf(stderr,"\n0**0 occured in function My_power"); exit(2);}
    else return 0;
  }
  int res=1; // ??? почему 1,0-double
  if(n>0) while(n-- > 0) res *= k;
  if(n<0) while(n++ < 0) res /= k;
  return res;
}

double sqr(double x){return x*x;}

double odd_sqrt(double x){if(x>=0) return sqrt(x); else return -sqrt(-x); }

//Function aliquot(x,y,eps) for positive x,y,eps
//returns 1 if exists N : fabs(x - N*y) < eps * y
//oterwise returns 0
int aliquot(double x, double y, double eps){
	if(x<0.0 || y<=0.0 || eps<=0.0 ){
	  fprintf(stderr,"\n\nArguments of function aliquot must be positive.");
	  fprintf(stderr,"\n\n Program will be terminated.");
	  pause();
	  exit(2);
	}
  double r = fmod(x,y), eps_y=eps*y;
  if(r<eps_y || y - r <eps_y) return 1;
  else return 0;
}

int freadnextkey(FILE *fc, char* found_key)
{//This utility finds and reads the next key including '=' and then stops.
 //Spaces and tabs inside a key are ignored.
 //found_key must contain '\0' only at its end
  char *p;
  char c=' '; //c - char from the stream.
  //int key_ready = 0;
  p=found_key;
  New_char: 
    if(*p=='\0'){
      fprintf(stderr,"ERROR: too long word %s... in a data file",found_key);
      exit(1);
    }
    c=fgetc(fc);
    if(c==EOF) return EOF;
    if(c=='\n'){p=found_key; goto New_char;}
    if(c==' '|| c=='\t') goto New_char;
    if(c=='#'){ while(c!='\n')c=fgetc(fc); p=found_key; goto New_char; } //skip comment

    if(c=='/'){ *p++=c; c=fgetc(fc);
      if(c=='/'){while(c!='\n'&&c!=EOF)c=fgetc(fc); p=found_key; goto New_char;} //skip comment
    }
    if(c != '='){*p++ = c; goto New_char;}
    *p='\0';
    //fprintf(stdout,"\nFound key %s",found_key);
    //pause();
    return 0;
}//freadnextkey

void ffindkey(FILE *fc, const char* sample_key)
{
  int key_is_found, key_is_wrong;
  char current_key[MAX_KEY_LEN];

  rewind(fc); 

  
  key_is_found=0;
  while(clear_buffer(MAX_KEY_LEN,current_key),//buffer must be cleared before each iteration
        key_is_found==0 && freadnextkey(fc,current_key)!=EOF){
    //freadnextkey(fc,current_key);
    char *p_current_key;
    const char *p_sample_key;
    p_current_key=current_key;
    p_sample_key=sample_key;

    key_is_wrong=0;
    while(*p_sample_key !='\0' && key_is_wrong==0){
      while(*p_sample_key == ' ' || *p_sample_key == '\t') p_sample_key++;
      if(*p_current_key++ != *p_sample_key++) key_is_wrong=1;
    }
    //Now if key_is_wrong==0 then *p_sample_key must be =='\0'.
    //So for matching of the two keys also *p_current_key must be =='\0'.
    //Thus:
    if(*p_current_key !='\0') key_is_wrong=1;
    if(key_is_wrong==0) key_is_found=1;

    if(key_is_wrong==1){//skip a word belonging to this key
       char c=' ';
       while(c==' ' || c=='\t') c=fgetc(fc);//skip blank char-s
       while( !(c==' ' || c=='\t' || c=='\n' ) ) c=fgetc(fc);//skip not blank char-s
    }
  }//while

  if(key_is_found==0){
    fprintf(stderr,"ERROR in a data file: key %s was not found",sample_key);
	pause();
    exit(1);
  }
  return;
}//ffindkey

 void fkeyread_string(FILE *fc, const char* sample_key, char* res, int res_length)
{//in the result string comment is removed
 //,i.e. slash '/' if followed by another '/' is replaced by '\0'
 char c;
 char *p_res;
 clear_buffer(res_length,res);
 //find the key, which is needed:
 ffindkey(fc, sample_key);//now the last character read is '='
 //skip possible blank characters after '=' :
 c=' ';
 while(c==' ' || c=='\t' || c=='\n') c=fgetc(fc);
 //now c == the first character of the result string

 p_res = res;
 while(c!=' ' && c!='\t' && c!='\n'){
   if(*p_res=='\0'){
     fprintf(stderr,"\nWARNING: string %s... in a data file is too long",res);
     break;
   }
   *p_res++ = c;
   c=fgetc(fc);
 }
 *p_res='\0';
 remove_comment(res);
 return;
}// fkeyread_string

 double fkeyread_double(FILE *fc, const char* sample_key, double *d)
{
  char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = strtod(buf, &endptr);
}//fkeyread_double

 int fkeyread_int(FILE *fc, const char* sample_key,int *d)
{
  //char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = atoi(buf);
}//fkeyread_int

 long fkeyread_long(FILE *fc, const char* sample_key, long *d)
{
  //char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = atol(buf);
}//fkeyread_int
//---------------------------------------

 void pause(void)
{
// char a='a';
// printf("\nPAUSE. Continue?(Y/N)\n");
// while( a!='y'&&a!='Y'&&a!='n'&&a!='N' ) a=getch();
// if(a=='n'||a=='N'||a=='Т'||a=='в') exit(1);
 char a='2';
 printf("\nPAUSE. Continue?(0/1)\n");
 while( a!='0'&& a!='1' ) a=getch();
 if(a=='0') exit(1);
 return;
}

 void pause(char* info)
{
// char a='a';
// printf("\nPAUSE. Continue?(Y/N)\n");
// while( a!='y'&&a!='Y'&&a!='n'&&a!='N' ) a=getch();
// if(a=='n'||a=='N'||a=='Т'||a=='в') exit(1);
 char a='2';
 printf("\n%s\nPAUSE. CONTINUE?(0/1)\n",info);
 while( a!='0'&& a!='1' ) a=getch();
 if(a=='0') exit(1);
 return;
}

 void PressAnyKey(void)
{
 printf("\nPress any key to continue\n");
 char a=getch();
 return;
}

 void skiplines(FILE *fc, int k_line)
{ char buf[MAXBUF];
  while(k_line-- >0) fgets(buf,MAXBUF,fc);
}
 double frd1d(FILE *fc, Real *d)
{ char buf[MAXBUF],buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,MAXBUF,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t')) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  return *d = atof(buf1);
}
 
 int frd1i(FILE *fc, int *k)
{
  char buf[MAXBUF], buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,MAXBUF,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  return *k = atoi(buf1);
}
 char* frd1s(FILE *fc, char *s)
{ char buf[MAXBUF];
  char *pb=buf, *pb1=s;
  int i=0;
  fgets(buf,MAXBUF,fc);
//  printf("buf=%s",buf);
  while(i<MALBUF && (*pb==' ' || *pb=='\t')) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
//  printf("2:s=%s",s);
  return s;
}
 
 void enter1Real(char *s, Real *d)
{
  printf("\n"); printf("%s",s); scanf("%lg",d);
}
 
 void enter1Index(char *s, int *i)
{
  printf("\n"); printf("%s ",s); scanf("%i",i);
}

  char* substr(char* s, int n1, int n2)
{
  if(n1>n2){
    strcpy(s,"");
    return s;
  }
  int i;
//  if(n1<0)n1=0;
  if(n2>=strlen(s)) n2=strlen(s)-1;
  char* res=s;
  char* x = s + n1;
  for(i=0;i<n2-n1+1;i++) res[i] = *x++;
  res[n2-n1+1]='\0';
  return strcpy(s,res);
}


  void delblanks(char* s)
{
  int i, j, lens=strlen(s);
  for(i=0;s[i]!='\0';i++){
   if(s[i]==' '){
    for(j=i;j<lens;j++) s[j]=s[j+1];
   }
  }
}
  void blanks2zeros(char* s)
{
  int i, lens=strlen(s);
  for(i=0;i<lens;i++) if(s[i]==' ') s[i]='0';
}

 void clear_buffer(int buffer_length, char * buffer)
{
  int i;
  for(i=0;i<buffer_length-1;i++) buffer[i]=' '; buffer[i]='\0';
  return;
}

//solution of combined equations by Kramer's method for e1*X+g1*Y=q1,
//														e2*X+g2*Y=q2,
void Kramer (double e1, double e2, double g1, double g2, double q1, double q2, double& X, double& Y)
{
	 double Det = e1*g2 - e2*g1;
	 double Det1 = q1*g2 - q2*g1;
	 double Det2 = e1*q2 - e2*q1;
	//	 if(fabs(Det) <= 1e-17){
	//		 printf("\nDeterminant of combined equations is %lg near 0", Det);
	//		 printf("\nDeterminant Det1 = %lg , Det2 = %lg", Det1, Det2);
	//		 pause();
	//	 }	
					
//				pause();
	 X = Det1/Det;
	 Y = Det2/Det;
}