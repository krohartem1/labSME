//diaproc.h
// must be included after  stdio.h
#ifndef _DIAPROC

#define _DIAPROC


#include <stdio.h>
#include "defglob.h"
#define MALBUF 16
#define MAX_KEY_LEN 80
#define MAX_SECTION_NAME_LEN 32
#define LINE_LEN 80

// -------------------new---------------------
struct DATA_SECTION
{
  char name[MAX_SECTION_NAME_LEN];
  long pos0;
  long pos1;
  //Constructor
  DATA_SECTION(void);
};

enum TASK {load, strain33, strain32, zhestkost33, zhestkost32 ,vibr32 ,vibr32_quake};
enum ETAP_TYPE {etap,cycle,endcycle};
struct ETAP
{
  //char name[MAX_SECTION_NAME_LEN];
  ETAP_TYPE type;
  int num;
};

struct CYCLE
{
  int activated;
  long k; //number of iterations
  long i;  //current iteration
  long pos0;//starting position in wzd-file, i.e. position after
};

//--------------------------------------------

 void pause(void);
 void pause(char* info);
 void PressAnyKey(void);
 void skiplines(FILE*
		,int //k_line
	       );
 double frd1d(FILE*, Real*);
 int frd1i(FILE*, int*);
 char* frd1s(FILE*, char*);
 void enter1Real(char*, Real*);
 void enter1Index(char*, int*);
 char* substr(char*, int, int );
 void delblanks(char* s);
 void blanks2zeros(char* s);
  
// ------------------new----------------------
 void Exit_error(const char* msg);
 void Exit_error(const char* msg, int n);
 void Exit_error(const char* msg, double x);
 void Exit_error(const char* msg, char* s);
 void Exit_error(const char* msg_with_int, const char* msg, int n);

 void clear_buffer(int buffer_length, char * buffer);
 int freadnextkey(FILE *fc, char* found_key);
 void ffindkey(FILE *fc, const char* sample_key);
    //This utility finds and reads the next key including '=' and then stops.
    //Spaces and tabs inside a key are ignored.
    //all characters from // until the end of the line are skipped as comments.

 void fkeyread_string(FILE *fc, const char* sample_key, char* res, int res_length);
 double fkeyread_double(FILE *fc, const char* sample_key, double *d);
 int fkeyread_int(FILE *fc, const char* sample_key, int *d);
 long fkeyread_long(FILE *fc, const char* sample_key, long *d);

 int freadnextsection(FILE *fc, DATA_SECTION& sect); //auxiliary utility
 int ffindsection(FILE *fc, const char* sample_name, DATA_SECTION& current_section);
 //int ffindsection(FILE *fc, const char* sample_name);//looks for section with name "sample_name" (0 - success)
 int ffindkey_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key);
 int fkeyread_string_in_section(FILE *fc, DATA_SECTION& sect ,const char* sample_key,char* res, int res_length);
 //three functions reading upon key "sample_key" from DATA_SECTION "sect" (return 0 on success)
 int fkeyread_double_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, double *d);
 int fkeyread_int_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, int *d);
 int fkeyread_long_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, long *d);

 int recognize_etap(ETAP& et, DATA_SECTION& sect);

Real my_square(Real x);
Real XHeavisideX(Real x);
Real Sign(Real x);

//These functions may not be made inline: they are not class-members

Real Min(Real a, Real b);
Real Min(Real a, Real b, Real c);
Real Min(Real a, Real b, Real c, Real d);
Real Max(Real a, Real b);
Real Max(Real a, Real b, Real c);
Real Max(Real a, Real b, Real c, Real d);
int is_numeral(char c);
void get_number(char* s);//procedure truncates string s leaving char-s that suit for a number
void remove_comment(char* s); //in the result string comment is removed
                              //,i.e. slash '/' if followed by another '/' is replaced by '\0'

//Function aliquot(x,y,eps) for positive x,y,eps
//returns 1 if exists N : fabs(x - N*y) < eps*y
//oterwise returns 0
int aliquot(double x, double y, double eps);

//function ceil_2 returns minimum n=2**k, such that n>=x
int ceil_2(double x);

//Function Hev = Heaviside(x)
//returns 1.0 if x>0 , 0.0 if x<0, and second argument if x==0
//by default Heaviside(0) is 0.0
double Hev(double x, double Hev_at_0=0.0);

//Function Hev2 = Heaviside(x*y) 
//returns 1.0 if x and y have the same sign and 0.0 otherwise
double Hev2(double x, double y);

//Function My_power(x,y)
//returns x**y when x>0 and returns 0 when x==0 && y>=0 
double My_pow(double x, double y);

//Function My_power(x,n)
//returns x**n for integer n and double x
double My_pow(double x, int n);

int My_pow(int k, int n);

double sqr(double x);
double odd_sqrt(double x);

//solution of combined equations by Kramer's method
void Kramer (double a1, double a2, double b1, double b2, double h1, 
			 double h2, double& X, double& Y);

#endif
