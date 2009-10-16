#include <stream.h>
#include <math.h>
#include <LEDA/real.h>
#include <LEDA/integer.h>
#include "ABDPY.h"

extern "C" int getpid();
extern "C" long clock();

void exact_out(ostream& s, double a)
{
 s.precision(20);
 s<<a;
 s.precision(6);
}


// random value generation
// return in a double a signed integer with i bits
double rnd(int i)
{ if (i < 25 )
	return  (((int) random()%2)?1:-1) *
	((double) ((int) random()% ((int)exp2(i)) ));
  else
	return  (((int) random()%2)?1:-1) *
	(((double) ((int) random()% ((int)exp2(i-25)) )) * exp2(25) +
	((double) ((int) random()% ((int)exp2(25)) )));
}


// Few 2x2 determinant sign computation
int quadruple_det2x2(double a, double b, double c, double d)
{
long double DET=(long double)a*(long double)d-(long double)b*(long double)c;
if (DET>0.0) return 1;
if (DET<0.0) return -1;
return 0;
}

int double_det2x2(double a, double b, double c, double d)
{ return ( ((a*d)<(b*c)) ? -1 : 1) ; }

int leda_integer_det2x2(double a, double b, double c, double d)
{ return sign( integer(a)*integer(d)-integer(b)*integer(c) ) ; }

int leda_real_det2x2(double a, double b, double c, double d)
{ real D=real(a)*real(d)-real(b)*real(c); return D.sign(1); }


int double_det3x3(double a, double b, double c,
		   double d, double e, double f, double g, double h, double i)
{ if ( a*(e*i-f*h) + b*(f*g-d*i) + c*(d*h-e*g) < 0 ) return -1 ;
  else return 1; }

// Few 3x3 determinant sign computation
int leda_integer_det3x3(double a, double b, double c, 
                         	  double d, double e, double f,
                         	  double g, double h, double i)
{ return sign( integer(a)* (integer(e)* integer(i) -integer(f)* integer(h) )
		    +  integer(b)* (integer(f)* integer(g) -integer(d)* integer(i) )
            +  integer(c)* (integer(d)* integer(h) -integer(e)* integer(g) ) );}

int leda_real_det3x3(double a, double b, double c, 
                         	  double d, double e, double f,
                         	  double g, double h, double i)
{ real D =     real(a)* (real(e)* real(i) - real(f)* real(h) )
		    +  real(b)* (real(f)* real(g) - real(d)* real(i) )
            +  real(c)* (real(d)* real(h) - real(e)* real(g) );
return D.sign(1); }






void generate2(double& a,double& b,double& c,double& d, char type)
{
  double cxx,cyy;
  switch (type) {
  case 'r' :               // random
    a = rnd(53);
    b = rnd(53);
	c = rnd(53);
    d = rnd(53);
	break;

  case 'p' :              // x+y=0
    a = rnd(53);
    b = -a;
	c = rnd(53);
    d = -c;
    break;

  case 'q' :              // aU bV
    cxx = rnd(26);
    cyy = rnd(26);
    a =   rnd(27);
    c =   rnd(27);
    b = cxx * a;
    d = cxx * c;
	a *= cyy;
	c *= cyy;
    break;

  case 'f' :              // y = floor(u*x)      u in [0,1]
    cxx = rnd(52);
    a = rnd(53);
    b = floor(a*cxx/exp2(52));
	c = rnd(53);
    d = floor(c*cxx/exp2(52));
	break;

  case 'c' :              // all entries equal
    a = b = c = d = rnd(53);
	break;
}}

void generate3(double& a,double& b,double& c,double& d,double& e,
			   double& f,double& g,double& h,double& i, char type)
{
  double cxx,cyy;
   switch (type) {
  case 'r' :              // random
    a = rnd(51);
    b = rnd(51);
	c = rnd(51);
    d = rnd(51);
    e = rnd(51);
    f = rnd(51);
    g = rnd(51);
	h = rnd(51);
    i = rnd(51);
	break;

  case 'q' :            // aU bV cU+dV
    cxx = rnd(25);
    cyy = rnd(25);
    a =   rnd(25);
    b =   rnd(25);
    c =   cxx * a + cyy * b;
    d =   rnd(26);
    e =   rnd(26);
    f =   cxx * d + cyy * e;
    g =   rnd(25);
    h =   rnd(25);
    i =   cxx * g + cyy * h;
    cxx = rnd(25);
    cyy = rnd(25);
	a *= cxx;
	d *= cxx;
	g *= cxx;
	b *= cyy;
	e *= cyy;
	h *= cyy;
    break;

  case 'p' :              // x+y+z=0
    a = rnd(50);
    b = rnd(50);
	c = -a -b;
    d = rnd(50);
    e = rnd(50);
    f = -d -e;
    g = rnd(50);
	h = rnd(50);
    i = -g -h;
    break;

  case 'f' :              // z = floor( t*x + u*y )     t,u in [0,1]
    cxx = rnd(51);
	cyy = rnd(51);
    a = rnd(51);
    b = rnd(51);
	c = floor((-a*cxx-b*cyy)/exp2(52));
    d = rnd(51);
    e = rnd(51);
    f = floor((-d*cxx-e*cyy)/exp2(52));
	g = rnd(51);
	h = rnd(51);
	i = floor((-g*cxx-h*cyy)/exp2(52));
	break;

  case 'c' :             // all entries equal
    a = b = c = d = e = f = g = h = i = rnd(51);
	break;
}}

void perturbation2(double& a,double& b,double& c,double& d)
{ a += rnd(2); b += rnd(2); c += rnd(2); d += rnd(2);}

void perturbation3(double& a,double& b,double& c,double& d,double& e,
			   double& f,double& g,double& h,double& i)
{ a += rnd(2); b += rnd(2); c += rnd(2); d += rnd(2); e += rnd(2); f += rnd(2);
	g += rnd(2); h += rnd(2); i += rnd(2); }

void transpose2(double& a,double& b,double& c,double& d)
{ double swap; (void) a; (void) d; swap = b; b = c; c = swap; }

void transpose3(double& a,double& b,double& c,double& d,double& e,
			   double& f,double& g,double& h,double& i)
{ double swap; swap = b; b = d; d = swap; swap = c; c = g; g = swap;
	swap = f; f = h; h = swap; (void) a; (void) e; (void) i; }





void interactive()
{
	int dim;
	double a,b,c,d,e,f,g,h,i;

	cout << " dim ? ";
	cin >> dim;
	while (1) if (dim == 2) {
	   cout << "a = "; cin >> a;
	   cout << "b = "; cin >> b;
	   cout << "c = "; cin >> c;
	   cout << "d = "; cin >> d;
	   cout << endl;
	   cout << " | " ; exact_out(cout,a);
       cout << " " ; exact_out(cout,c);
       cout << " |" << endl;
	   cout << " | " ; exact_out(cout,b);
       cout << " " ; exact_out(cout,d);
       cout << " |" << endl;
	 cout << "double :" << double_det2x2(a,b,c,d) << endl;
	 cout << "quadruple :" << quadruple_det2x2(a,b,c,d) << endl;
	 cout << "integer :" << leda_integer_det2x2(a,b,c,d) << endl;
	 cout << "real:" << leda_real_det2x2(a,b,c,d) << endl;
	 cout << "ABDPY not lazy:" << ABDPY__not_lazy_det2x2(a,b,c,d) << endl;
	 cout << "ABDPY lazy:" << ABDPY_det2x2(a,b,c,d) << endl;
	 cout << "ABDPY secure:" << ABDPY_secure_det2x2(a,b,c,d) <<  endl;
	} else if (dim == 3) {
	   cout << "a = "; cin >> a;
	   cout << "b = "; cin >> b;
	   cout << "c = "; cin >> c;
	   cout << "d = "; cin >> d;
	   cout << "e = "; cin >> e;
	   cout << "f = "; cin >> f;
	   cout << "g = "; cin >> g;
	   cout << "h = "; cin >> h;
	   cout << "i = "; cin >> i;
	   cout << " | " ; exact_out(cout,a);
       cout << " " ; exact_out(cout,d);
       cout << " " ; exact_out(cout,g);
       cout << " |" << endl;
	   cout << " | " ; exact_out(cout,b);
       cout << " " ; exact_out(cout,e);
       cout << " " ; exact_out(cout,h);
       cout << " |" << endl;
	   cout << " | " ; exact_out(cout,c);
       cout << " " ; exact_out(cout,f);
       cout << " " ; exact_out(cout,i);
       cout << " |" << endl;
	 cout << "double :" << double_det3x3(a,b,c,d,e,f,g,h,i) << endl;
	 cout << "integer :" << leda_integer_det3x3(a,b,c,d,e,f,g,h,i) << endl;
	 cout << "real :" << leda_real_det3x3(a,b,c,d,e,f,g,h,i) << endl;
	 cout << "ABDPY not lazy:"<<ABDPY__not_lazy_det3x3(a,b,c,d,e,f,g,h,i)<<endl;
	 cout << "ABDPY lazy:" << ABDPY_det3x3(a,b,c,d,e,f,g,h,i) << endl;
	 cout << "ABDPY secure:" << ABDPY_secure_det3x3(a,b,c,d,e,f,g,h,i) << endl;
	}
}




main(int argc, char **argv)
{ double a,b,c,d,e,f,g,h,i;
  int sign;
  int n,N=80000;
  int input, method, trial, j;
  double time[10];
  long t1,t2;
  int (*FF2)(double, double, double, double);
  int (*FF3)(double, double, double, double, double, double, double, double, double);

 srandom(n=getpid());
 cerr <<"init random with pid ="<<n<<endl;

 if (--argc)
 {
	++argv;
	if (    (**argv>='1') && (**argv<='9') ) n  = atoi(*argv);
	else interactive();
 }
 else n = 10;

  if (!argc){
  cout <<"\\documentstyle{article}" << endl;
  cout <<"\\voffset=-2.5cm" << endl;
  cout <<"\\hoffset=-2.5cm" << endl;
  cout <<"\\textwidth              17.5cm" << endl;
  cout <<"\\textheight             24.2cm" << endl;
  cout <<"\\begin{document}" << endl;
  cout <<"\\title{Test of determinant sign algorithms}" << endl;
  cout <<"\\author{O. Devillers}" << endl;
  cout <<"\\maketitle" << endl;
  cout << endl;
  cout <<"Few determinant sign algorithms" << endl;
  cout <<"are tested with various kinds of input." << endl;
  cout <<"For each kind of input, "<<n<<" differents trials are done." << endl;
  cout <<"For each, the determinant is computed "<<N<<" times," << endl;
  cout <<"and the time are obtained using the {\\tt clock} function."  << endl;
  cout << endl;
  cout <<"Thus the precision on the given running time is about " << endl;
  cout <<16000/(double)N << "$\\mu$s." << endl;
  cout << endl;
  cout <<"The program verify that all methods" << endl;
  cout <<"produces the same result, except" << endl;
  cout <<"the direct floating point calculation" << endl;
  cout <<"which can be wrong due to" << endl;
  cout <<"rounding error." << endl;
  cout << endl;
  cout <<"\\section*{Platform}" << endl;
  cout <<"% ADD HERE YOUR PLATFORM" << endl;
  cout <<"Computer :\\\\" << endl;
  cout <<"Compiler :\\\\" << endl;
  cout <<"Compiler options : " << endl;
  cout << endl;
  cout <<"\\section*{Algorithms}" << endl;
  cout <<"Here are the different algorithms used," << endl;
  cout <<"the original data are of {\\tt double} type." << endl;
  cout << endl;
  cout <<"\\subsection*{Direct computation}" << endl;
  cout <<"Direct computation of the determinant" << endl;
  cout <<"using $ad-bc$ formula in 2D and" << endl;
  cout <<"$a(ei-fh) + b(fg-di) + c(dh-eg)$ in 3D." << endl;
  cout <<"Plus an evaluation of the sign." << endl;
  cout <<"\\subsubsection*{double}" << endl;
  cout <<"If the computation is done using {\\tt double} arithmetic," << endl;
  cout <<"it is fast, but there is rounding errors." << endl;
  cout <<"\\subsubsection*{quadruple}" << endl;
  cout <<"If the computation is done using quadruple precision" << endl;
  cout <<"({\\tt long double} type) the result is exact in 2D." << endl;
  cout <<"\\subsubsection*{leda-integer}" << endl;
  cout <<"LEDA provides exact computation on integer of arbitrary" << endl;
  cout <<"length. The result is exact \\cite{leda}." << endl;
  cout <<"\\subsubsection*{leda-real}" << endl;
  cout <<"LEDA provides a floating point filter which invokes" << endl;
  cout <<"an exact computation only when the rounded computation" << endl;
  cout <<"does not allow a conclusion \\cite{leda}." << endl;
  cout <<"The result is exact." << endl;
  cout <<"\\subsection*{Iterative method}" << endl;
  cout <<"Avnaim, Boissonnat, Devillers, Preparata and Yvinec" << endl;
  cout <<"provides an iterative method which transform the " << endl;
  cout <<"determinant in another one with smaller entries" << endl;
  cout <<"up the entries have sign allowing a conclusion" << endl;
  cout <<"\\cite{iter1,iter2}." << endl;
  cout <<"The algorithm uses double arithmetic to perform" << endl;
  cout <<"fast operations on 53 bits integers, thus there is precondition" << endl;
  cout <<"that the entries must represent an integer and their" << endl;
  cout <<"absolute value must be smaller than $2^{53}$ in 2D" << endl;
  cout <<"and $2^{51}$ in 3D." << endl;
  cout <<"\\subsubsection*{not-lazy}" << endl;
  cout <<"The algorithm is iterated, up to the comparison" << endl;
  cout <<"of the entries and computation of sign of minors" << endl;
  cout <<"in 3D allows a conclusion." << endl;
  cout <<"\\subsubsection*{lazy}" << endl;
  cout <<"The algorithm is combinated with a floating point filter." << endl;
  cout <<"\\subsubsection*{secure}" << endl;
  cout <<"The algorithm is combinated with a floating point filter," << endl;
  cout <<"and the entries are verified to match the preconditions." << endl;
  cout << endl;
  cout <<"\\small" << endl;
  cout <<"\\section*{Input 2D}" << endl;
  cout <<"$a$ is random on $b$ bits means a is evenly distributed" << endl;
  cout <<"among the integers between $-2^b+1$ and $2^b-1$." << endl;
  cout <<"\\subsection*{random}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a& c\\\\ b & d \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d \\mbox{ random on 53 bits} $" << endl;
  cout <<"\\subsection*{$x=-y$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a& c\\\\ -a & -c \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,c \\mbox{ random on 53 bits} $" << endl;
  cout <<"\\subsection*{$x=-y+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a+e_a& c+e_c\\\\ -a+e_b & -c+e_d \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,c \\mbox{ random on 53 bits, } e_a,e_b,e_c,e_d \\mbox{ on 2 bits}$" << endl;
  cout <<"\\subsection*{$x=-y^t$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a& -a\\\\ b & -b \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,c \\mbox{ random on 53 bits} $" << endl;
  cout <<"\\subsection*{$x=-y+\\varepsilon^t$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a+e_a& -a+e_c\\\\ b+e_b & -b+e_d \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b \\mbox{ random on 53 bits, } e_a,e_b,e_c,e_d \\mbox{ on 2 bits}$" << endl;
  cout <<"\\subsection*{$kU,lU$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} ka& la\\\\ kb & lb \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b \\mbox{ random on 26 bits, } k,l \\mbox{ on 27 bits}$" << endl;
  cout <<"\\subsection*{$kU,lU+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} ka+e_a& la+e_c\\\\ kb+e_b & lb+e_d \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b \\mbox{ random on 53 bits, } k,l \\mbox{ on 27 bits, } e_a,e_b,e_c,e_d \\mbox{ on 2 bits}$" << endl;
  cout <<"\\subsection*{$U,\\lfloor\\alpha U\\rfloor$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a & c\\\\ \\lfloor \\alpha a\\rfloor & \\lfloor \\alpha c\\rfloor \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,c \\mbox{ random on 53 bits, } \\alpha \\mbox{ random in [-1,1]}$" << endl;
  cout <<"\\subsection*{$U,\\lfloor\\alpha U\\rfloor^t$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a& \\lfloor \\alpha a\\rfloor\\\\ b & \\lfloor \\alpha b\\rfloor \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b \\mbox{ random on 53 bits, } \\alpha \\mbox{ random in [-1,1]}$" << endl;
  cout <<"\\subsection*{$=$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a& a\\\\ a & a \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a \\mbox{ random on 53 bits} $" << endl;
  cout <<"\\subsection*{$=+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{cc} a+e_a& a+e_c\\\\ a+e_b & a+e_d \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a \\mbox{ random on 53 bits, } e_a,e_b,e_c,e_d \\mbox{ on 2 bits}$" << endl;
  cout << endl;
  cout <<"\\subsection*{Running time in $\\mu$s.}" << endl;
  cout <<"\\begin{tabular}{|l||r|r|r|r|r|r|r|}\\hline" << endl;
  cout <<"input&double&quadruple&leda-integer&leda-real&not-lazy&lazy&secure\\\\\\hline" << endl;
  }

 for ( input = 0 ; input < 11 ; input++ )
 {
	switch (input) {
	  case  0: cout<< "random"; break;
	  case  1: cout<< "$x=-y$"; break;
	  case  2: cout<< "$x=-y+\\varepsilon$"; break;
	  case  3: cout<< "$x=-y^t$"; break;
	  case  4: cout<< "$x=-y+\\varepsilon^t$"; break;
	  case  5: cout<< "$kU,lU$"; break;
	  case  6: cout<< "$kU,lU+\\varepsilon$"; break;
	  case  7: cout<< "$U,\\lfloor\\alpha U\\rfloor$"; break;
	  case  8: cout<< "$U,\\lfloor\\alpha U\\rfloor^t$"; break;
	  case  9: cout<< "$=$"; break;
	  case 10: cout<< "$=+\\varepsilon$"; break;
	}

	for ( method = 0 ; method < 7 ; method++)  time[method] = 0.0;
	for (trial=0; trial < n ; trial++){
	  switch (input) {
	  	case  0: generate2(a,b,c,d,'r'); break;
	  	case  1: generate2(a,b,c,d,'p'); break;
	  	case  2: generate2(a,b,c,d,'p'); perturbation2(a,b,c,d); break;
	  	case  3: generate2(a,b,c,d,'p'); transpose2(a,b,c,d); break;
	  	case  4: generate2(a,b,c,d,'p'); perturbation2(a,b,c,d);transpose2(a,b,c,d); break;
	  	case  5: generate2(a,b,c,d,'q'); break;
	  	case  6: generate2(a,b,c,d,'q'); perturbation2(a,b,c,d); break;
	  	case  7: generate2(a,b,c,d,'f'); break;
	  	case  8: generate2(a,b,c,d,'f'); transpose2(a,b,c,d); break;
	  	case  9: generate2(a,b,c,d,'c'); break;
	  	case 10: generate2(a,b,c,d,'c'); perturbation2(a,b,c,d); break;
	  }
	  sign =  leda_integer_det2x2(a,b,c,d);
	  for ( method = 0 ; method < 7 ; method++) {
		switch (method) {
          case 0 : FF2 = double_det2x2; break;
          case 1 : FF2 = quadruple_det2x2; break;
		  case 2 : FF2 = leda_integer_det2x2; break;
          case 3 : FF2 = leda_real_det2x2;                   N/=100; break;
          case 4 : FF2 = ABDPY__not_lazy_det2x2;             N*=100; break;
          case 5 : FF2 = ABDPY_det2x2; break;
          case 6 : FF2 = ABDPY_secure_det2x2; break;
		}
		if ( method) if (sign != FF2(a,b,c,d) )
		{ cerr << "ERROR det2x2- send the following lines to "
                         << "olivier.devillers@sophia.inria.fr"
                         <<input<<" "<<trial<< " double "
                         << double_det2x2(a,b,c,d) << ", quadruple "
                         << quadruple_det2x2(a,b,c,d) << ", integer "
		                 << leda_integer_det2x2(a,b,c,d) << ", real"
                         << leda_real_det2x2(a,b,c,d) << ", ABDPY"
                         << ABDPY__not_lazy_det2x2(a,b,c,d) << ", ABDPY"
                         << ABDPY_det2x2(a,b,c,d) << ", ABDPY"
                         << ABDPY_secure_det2x2(a,b,c,d) <<  endl;
          cerr<<"a= ";exact_out(cerr,a); cerr<<endl;
          cerr<<"b= ";exact_out(cerr,b); cerr<<endl;
          cerr<<"c= ";exact_out(cerr,c); cerr<<endl;
          cerr<<"d= ";exact_out(cerr,d); cerr<<endl;
		  cerr<<endl<<"send the above lines to olivier.devillers@sophia.inria.fr"<<endl;
		}
		if (!argc)
		{	
			t1 = clock();
			for (j=0; j<N; j++)
							(void) FF2(a,b,c,d);
			t2 = clock()-t1;
		}
		time[method] += t2 / (double) N;
	  } // end for method
	} //end for trial
	if(!argc)for ( method = 0 ; method < 7 ; method++)cout<<"&"<<time[method]/n;
	if(!argc)cout << "\\\\\\hline";
    cout << endl;
  }  // end for input


  if (!argc){
  cout <<"\\end{tabular}" << endl;
  cout << endl;
  cout <<"\\section*{Input 3D}" << endl;
  cout <<"\\subsection*{random}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&d&g\\\\b&e&h\\\\c&f&i \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d,e,f,g,h,i \\mbox{ random on 51 bits} $" << endl;
  cout <<"\\subsection*{$x+y+z=0$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&d&g\\\\b&e&h\\\\-a-b&-d-e&-g-h \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,d,e,g,h \\mbox{ random on 51 bits} $" << endl;
  cout <<"\\subsection*{$x+y+z=0+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a+e_a&d+e_d&g+e_g\\\\b+e_b&e+e_e&h+e_h\\\\-a-b+e_c&-d-e+e_f&-g-h+e_i \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,d,e,g,h \\mbox{ random on 51 bits, }" << endl;
  cout <<"e_a,e_b,e_c,e_d,e_e,e_f,e_g,e_h,e_i \\mbox{ on 2 bits}$" << endl;
  cout <<"\\subsection*{$x+y+z=0^t$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&d&-a-d\\\\b&e&-b-e\\\\c&f&-c-f \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d,e,f \\mbox{ random on 51 bits} $" << endl;
  cout <<"\\subsection*{$x+y+z=0+\\varepsilon^t$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a+e_a&d+e_d&-a-d+e_g\\\\b+e_b&e+e_e&-b-e+e_h\\\\c+e_c&f+e_f&-c-f+e_i \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d,e,f \\mbox{ random on 51 bits, }" << endl;
  cout <<"e_a,e_b,e_c,e_d,e_e,e_f,e_g,e_h,e_i \\mbox{ on 2 bits}$" << endl;
  cout <<"\\subsection*{$kU,lV,mU+nV$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} ka&ld&ma+nd\\\\kb&le&mb+ne\\\\kc&lf&mc+lf \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d,e,f,m,n \\mbox{ random on 25 bits, } k,l \\mbox{ on 26 bits}$" << endl;
  cout <<"\\subsection*{$kU,lV,mU+nV+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} ka+e_a&ld+e_d&ma+nd+e_g\\\\kb+e_b&le+e_e&mb+ne+e_h\\\\kc+e_c&lf+e_f&mc+lf+e_i \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} \\begin{array}{l}a,b,c,d,e,f,m,n \\mbox{ random on 25 bits, } k,l \\mbox{ on 26 bits,}\\\\" << endl;
  cout <<"e_a,e_b,e_c,e_d,e_e,e_f,e_g,e_h,e_i \\mbox{ on 2 bits} \\end{array}$" << endl;
  cout <<"\\subsection*{$U,V,\\lfloor\\alpha U+\\beta V\\rfloor$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&d&g\\\\b&e&h\\\\" << endl;
  cout <<"\\lfloor \\alpha a+\\beta b\\rfloor&\\lfloor \\alpha d+\\beta e\\rfloor&" << endl;
  cout <<"\\lfloor \\alpha g+\\beta h\\rfloor \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,d,e,g,h \\mbox{ random on 51 bits, } \\alpha,\\beta \\mbox{ random in [-1/2,1/2]}$" << endl;
  cout <<"\\subsection*{$U,V,\\lfloor\\alpha U+\\beta V\\rfloor^t$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&d&\\lfloor \\alpha a+\\beta d\\rfloor\\\\" << endl;
  cout <<"b&e&\\lfloor \\alpha b+\\beta e\\rfloor\\\\" << endl;
  cout <<"c&f&\\lfloor \\alpha c+\\beta f\\rfloor \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a,b,c,d,e,f \\mbox{ random on 51 bits, } \\alpha,\\beta \\mbox{ random in [-1/2,1/2]}$" << endl;
  cout <<"\\subsection*{$=$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a&a&a\\\\a&a&a\\\\a&a&a \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a \\mbox{ random on 53 bits} $" << endl;
  cout <<"\\subsection*{$=+\\varepsilon$}" << endl;
  cout <<"$\\left| \\begin{array}{ccc} a+e_a&d+e_d&g+e_g\\\\b+e_b&e+e_e&h+e_h\\\\-a-b+e_c&-d-e+" << endl;
  cout <<"e_f&-g-h+e_i \\end{array} \\right|" << endl;
  cout <<"\\hspace{1cm} a \\mbox{ random on 53 bits, }" << endl;
  cout <<"e_a,e_b,e_c,e_d,e_e,e_f,e_g,e_h,e_i \\mbox{ on 2 bits}$" << endl;
  cout << endl;
  cout <<"\\subsection*{Running time in $\\mu$s.}" << endl;
  cout <<"\\begin{tabular}{|l||r|r|r|r|r|r|}\\hline" << endl;
  cout <<"input&double&leda-integer&leda-real&not-lazy&lazy&secure\\\\\\hline" << endl;
  }

 for ( input = 0 ; input < 11 ; input++ )
 {
	switch (input) {
	  case  0: cout<< "random"; break;
	  case  1: cout<< "$x+y+z=0$"; break;
	  case  2: cout<< "$x+y+z=0+\\varepsilon$"; break;
	  case  3: cout<< "$x+y+z=0^t$"; break;
	  case  4: cout<< "$x+y+z=0+\\varepsilon^t$"; break;
	  case  5: cout<< "$kU,lV,mU+nV$"; break;
	  case  6: cout<< "$kU,lV,mU+nV+\\varepsilon$"; break;
	  case  7: cout<< "$U,V,\\lfloor\\alpha U+\\beta V\\rfloor$"; break;
	  case  8: cout<< "$U,V,\\lfloor\\alpha U+\\beta V\\rfloor^t$"; break;
	  case  9: cout<< "$=$"; break;
	  case 10: cout<< "$=+\\varepsilon$"; break;
	}
	for ( method = 0 ; method < 6 ; method++)  time[method] = 0.0;
	for (trial=0; trial < n ; trial++){
	  	switch (input) {
	  	case  0: generate3(a,b,c,d,e,f,g,h,i,'r');  break;
	  	case  1: generate3(a,b,c,d,e,f,g,h,i,'p');  break;
	  	case  2: generate3(a,b,c,d,e,f,g,h,i,'p'); perturbation3(a,b,c,d,e,f,g,h,i); break;
	  	case  3: generate3(a,b,c,d,e,f,g,h,i,'p'); transpose3(a,b,c,d,e,f,g,h,i); break;
	  	case  4: generate3(a,b,c,d,e,f,g,h,i,'p'); perturbation3(a,b,c,d,e,f,g,h,i);transpose3(a,b,c,d,e,f,g,h,i);  break;
	  	case  5: generate3(a,b,c,d,e,f,g,h,i,'q');  break;
	  	case  6: generate3(a,b,c,d,e,f,g,h,i,'q'); perturbation3(a,b,c,d,e,f,g,h,i); break;
	  	case  7: generate3(a,b,c,d,e,f,g,h,i,'f'); break;
	  	case  8: generate3(a,b,c,d,e,f,g,h,i,'f'); transpose3(a,b,c,d,e,f,g,h,i); break;
	  	case  9: generate3(a,b,c,d,e,f,g,h,i,'c'); break;
	  	case 10: generate3(a,b,c,d,e,f,g,h,i,'c'); perturbation3(a,b,c,d,e,f,g,h,i); break;
	  	}
		sign =  leda_integer_det3x3(a,b,c,d,e,f,g,h,i);
	  for ( method = 0 ; method < 6 ; method++) {
		switch (method) {
          case 0 : FF3 = double_det3x3; break;
		  case 1 : FF3 = leda_integer_det3x3; break;
          case 2 : FF3 = leda_real_det3x3;                 N/=100; break;
          case 3 : FF3 = ABDPY__not_lazy_det3x3;           N*=100; break;
          case 4 : FF3 = ABDPY_det3x3; break;
          case 5 : FF3 = ABDPY_secure_det3x3; break;
		}
		if ( method) if (sign != FF3(a,b,c,d,e,f,g,h,i) )
		{ cerr << "ERROR det3x3- send the following lines to "
                         << "olivier.devillers@sophia.inria.fr"
                         <<input<<" "<<trial<< " double "
                         << double_det3x3(a,b,c,d,e,f,g,h,i) << ", integer "
		                 << leda_integer_det3x3(a,b,c,d,e,f,g,h,i) << ", real "
                         << leda_real_det3x3(a,b,c,d,e,f,g,h,i) << ", ABDPY "
                         <<ABDPY__not_lazy_det3x3(a,b,c,d,e,f,g,h,i)<<", ABDPY "
                         << ABDPY_det3x3(a,b,c,d,e,f,g,h,i) << ", ABDPY "
                         << ABDPY_secure_det3x3(a,b,c,d,e,f,g,h,i) << endl;
          cerr<<"a= ";exact_out(cerr,a); cerr<<endl;
          cerr<<"b= ";exact_out(cerr,b); cerr<<endl;
          cerr<<"c= ";exact_out(cerr,c); cerr<<endl;
          cerr<<"d= ";exact_out(cerr,d); cerr<<endl;
          cerr<<"e= ";exact_out(cerr,e); cerr<<endl;
          cerr<<"f= ";exact_out(cerr,f); cerr<<endl;
          cerr<<"g= ";exact_out(cerr,g); cerr<<endl;
          cerr<<"h= ";exact_out(cerr,h); cerr<<endl;
          cerr<<"i= ";exact_out(cerr,i); cerr<<endl;
		  cerr<<endl<<"send the above lines to olivier.devillers@sophia.inria.fr"<<endl;
        }
		if (!argc)
		{
			t1 = clock();
			for (j=0; j<N; j++)
				(void) FF3(a,b,c,d,e,f,g,h,i);
			t2 = clock()-t1;
		}
		time[method] += t2 / (double) N;
	  } // end for method
	} //end for trial
	if(!argc)for ( method = 0 ; method < 6 ; method++)  cout <<"&" <<  time[method]/n;
	if(!argc)cout << "\\\\\\hline";
    cout << endl;
  }  // end for input


  if (!argc){
  cout <<"\\end{tabular}" << endl;
  cout << endl;
  cout << endl;
  cout <<"\\newcommand{\\etalchar}[1]{$^{#1}$}" << endl;
  cout <<"\\begin{thebibliography}{ABD{\\etalchar{+}}94}" << endl;
  cout <<"\\bibitem[ABD{\\etalchar{+}}94]{iter1}" << endl;
  cout <<"F.~Avnaim, J.-D. Boissonnat, O.~Devillers, F.~Preparata, and M.~Yvinec." << endl;
  cout <<"\\newblock Evaluating signs of determinants using single-precision arithmetic." << endl;
  cout <<"\\newblock Research {Report} 2306, INRIA, BP93, 06902 Sophia-Antipolis, France, 1994." << endl;
  cout <<"\\bibitem[ABD{\\etalchar{+}}95]{iter2}" << endl;
  cout <<"F.~Avnaim, J.-D. Boissonnat, O.~Devillers, F.~Preparata, and M.~Yvinec." << endl;
  cout <<"\\newblock Evaluating signs of determinants using single-precision arithmetic." << endl;
  cout <<"\\newblock In {\\em Proc. 11th Annu. ACM Sympos. Comput. Geom.},1995." << endl;
  cout <<"{\\tt http://www.inria.fr:/prisme/personnel/devillers/anglais/determinant.html}" << endl;
  cout <<"\\bibitem[MBK{\\etalchar{+}}95]{leda}" << endl;
  cout <<"K.~Mehlhorn, C.~Burnikel, J.~K\\\"onnemann, S.~N\\\"aher, S.~Schirra, and C.~Uhrig." << endl;
  cout <<"\\newblock Exact Geometric Computation in LEDA" << endl;
  cout <<"\\newblock In {\\em Proc. 11th Annu. ACM Sympos. Comput. Geom.},1995." << endl;
  cout <<"{\\tt http://www.mpi-sb.mpg.de/guide/staff/uhrig/leda.html}" << endl;
  cout <<"\\end{thebibliography}" << endl;
  cout <<"\\end{document}" << endl;
  }

}
