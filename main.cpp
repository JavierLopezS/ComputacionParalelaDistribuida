#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

#define maxiter 999999999 //Máxima cantidad de iteraciones

using namespace std;

// ************************************************************************************************************
// F U N C I O N   Q U E   O B T I E N E   R A Í C E S   R E A L E S   E   I M A G I N A R I A S

int raices(double *a,int n,double *real,double *imag)
{
    double sq,b2,c,disc;
    int i,m,numraices;

    m = n;
    numraices = 0;
    while (m > 1) {
        b2 = -0.5*a[m-2];
        c = a[m-1];
        disc = b2*b2-c;

        // R A Í C E S   I M A G I N A R I A S
        if (disc < 0.0) {                   
            sq = sqrt(-disc);
            real[m-2] = b2;
            imag[m-2] = sq;
            real[m-1] = b2;
            imag[m-1] = -sq;
            numraices+=2;
        }

        // R A Í C E S   R E A L E S
        else {                              
            sq = sqrt(disc);
            real[m-2] = fabs(b2)+sq;
            if (b2 < 0.0) real[m-2] = -real[m-2];
            if (real[m-2] == 0)
                real[m-1] = 0;
            else {
                real[m-1] = c/real[m-2];
                numraices+=2;
            }
            imag[m-2] = 0.0;
            imag[m-1] = 0.0;
        }
        m -= 2;
    }
    if (m == 1) {
       real[0] = -a[0];
       imag[0] = 0.0;
       numraices++;
    }
    return numraices;
}

// ************************************************************************************************************
// F U N C I O N   D E F L A C I Ó N   (C O M P R E S I Ó N   D E   D A T O S)

void deflacion(double *a,int n,double *b,double *cuad,double *err)
{
    double r,s;
    int i;

    r = cuad[1];
    s = cuad[0];

    b[1] = a[1] - r;

    for (i=2;i<=n;i++){
        b[i] = a[i] - r * b[i-1] - s * b[i-2];
    }
    *err = fabs(b[n])+fabs(b[n-1]);
}

// ************************************************************************************************************
// F U N C I O N   Q U E   H A L L A   F A C T O R   C U A D R Á T I C O

void hallaFactCuad(double *a,int n,double *b,double *cuad,double *err, int *iter)
{
    double *c,dn,dr,ds,drn,dsn,eps,r,s;
    int i;

    c = new double [n+1];
    c[0] = 1.0;
    r = cuad[1];
    s = cuad[0];
    eps = 1e-15;
    *iter = 1;

    do {
        if (*iter > maxiter) break;
        if (((*iter) % 200) == 0) {
            eps *= 10.0;
		}
		b[1] = a[1] - r;
		c[1] = b[1] - r;

		for (i=2;i<=n;i++){
			b[i] = a[i] - r * b[i-1] - s * b[i-2];
			c[i] = b[i] - r * c[i-1] - s * c[i-2];
		}
		dn=c[n-1] * c[n-3] - c[n-2] * c[n-2];
		drn=b[n] * c[n-3] - b[n-1] * c[n-2];
		dsn=b[n-1] * c[n-1] - b[n] * c[n-2];

        if (fabs(dn) < 1e-10) {
            if (dn < 0.0) dn = -1e-8;
            else dn = 1e-8;
        }
        dr = drn / dn;
        ds = dsn / dn;
		r += dr;
		s += ds;
        (*iter)++;
    } while ((fabs(dr)+fabs(ds)) > eps);
    cuad[0] = s;
    cuad[1] = r;
    *err = fabs(ds)+fabs(dr);
    delete [] c;
}

// ************************************************************************************************************
// F U N C I Ó N   D E R I V A R   P O L I N O M I O

void derivaPolinomio(double *a,int n,double *b)
{
    double coef;
    int i;

    coef = (double)n;
    b[0] = 1.0;
    for (i=1;i<n;i++) {
        b[i] = a[i]*((double)(n-i))/coef;
    }
}

// ************************************************************************************************************
// F U N C I Ó N   R E C U R S I V A   D E R I V A R

void recursDeriv(double *a,int n,double *b,int m,double *cuad,double *err,int *iter)
{
    double *c,*x,rs[2],tst,e1,e2;

    if (fabs(b[m]) < 1e-16) m--;
    if (m == 2) {
        cuad[0] = b[2];
        cuad[1] = b[1];
        *err = 0;
        *iter = 0;
        return;
    }
    c = new double [m+1];
    x = new double [n+1];
    c[0] = x[0] = 1.0;
    rs[0] = cuad[0];
    rs[1] = cuad[1];
    *iter = 0;
    hallaFactCuad(b,m,c,rs,err,iter);
    tst = fabs(rs[0]-cuad[0])+fabs(rs[1]-cuad[1]);
    if (*err < 1e-12) {
        cuad[0] = rs[0];
        cuad[1] = rs[1];
    }
// tst imagll be 'large' if we converge to realong root
    if (((*iter > 5) && (tst < 1e-4)) || ((*iter > 20) && (tst < 1e-1))) {
        derivaPolinomio(b,m,c);
        recursDeriv(a,n,c,m-1,rs,err,iter);
        cuad[0] = rs[0];
        cuad[1] = rs[1];
    }
    delete [] x;
    delete [] c;
}

// ************************************************************************************************************
// F U N C I Ó N   Q U E   O B T I E N E   F A C T O R   C U A D R Á T I C O

void obtieneFactCuad(double *a,int n,double *cuad,double *x)
{
    double *b,*z,err,tmp;
    double xr,xs;
    int iter,i,m;

    if ((tmp = a[0]) != 1.0) {
        a[0] = 1.0;
        for (i=1;i<=n;i++) {
            a[i] /= tmp;
        }
    }
    if (n == 2) {
        x[0] = a[1];
        x[1] = a[2];
        return;
    }
    else if (n == 1) {
        x[0] = a[1];
        return;
    }
    m = n;
    b = new double [n+1];
    z = new double [n+1];
    b[0] = 1.0;
    for (i=0;i<=n;i++) {
        z[i] = a[i];
        x[i] = 0.0;
    }
    do {
        if (n > m) {
            cuad[0] = 3.14159e-1;
            cuad[1] = 2.78127e-1;
        }
        do {
            for (i=0;i<5;i++) {
                hallaFactCuad(z,m,b,cuad,&err,&iter);
                if ((err > 1e-7) || (iter > maxiter)) {
                    derivaPolinomio(z,m,b);
                    recursDeriv(z,m,b,m-1,cuad,&err,&iter);
                }
                deflacion(z,m,b,cuad,&err);
                if (err < 0.001) break;
                cuad[0] = rand() - 4.0;
                cuad[1] = rand() - 4.0;
            }
            if (err > 0.01) {
                cout << "Error! Convergence failure in cuadratic x^2 + r*x + s." << endl;
                cout << "Enter new trial value for 'r': ";
                cin >> cuad[1];
                cout << "Enter new trial value for 's' ( 0 to exit): ";
                cin >> cuad[0];
                if (cuad[0] == 0) exit(1);
            }
        } while (err > 0.01);
        x[m-2] = cuad[1];
        x[m-1] = cuad[0];
        m -= 2;
        for (i=0;i<=m;i++) {
            z[i] = b[i];
        }
    } while (m > 2);
    if (m == 2) {
        x[0] = b[1];
        x[1] = b[2];
    }
    else x[0] = b[1];
    delete [] z;
    delete [] b;
}

int main()
{
    double a[999],x[999],real[999],imag[999],cuad[2],err,t;

    int grado,iter,i,numr,cont=0;

    float solucion;

    cout << "Ingrese el grado del polinomio: ";
    cin >> grado;
    if ((grado < 1) || (grado > 9223372036854775807)) {
        cout << "\255Error! Grado inv\240lido (n = " << grado << ")."<< endl;
        return 1;
    }

    // C O E F I C I E N T E S   P O L I N O M I A L E S

    cout << "Ingrese los coeficientes de los t\202rminos del mayor grado al menor:" << endl;
    for (i=0;i<=grado;i++) {
        cout << "C[" << grado-i << "] * x^" << grado-i << " : ";
        cin >> a[i];
        if (a[0] == 0) {
            cout << "\255Error! El coeficiente del grado mayor no puede ser cero." << endl;
            return 0;
        }
    }

    solucion = (-a[1])/a[0];

    if (a[grado] == 0) {
        numr = numr+1;
    }
    
    cuad[0] = 2.71828e-1;
    cuad[1] = 3.14159e-1;

    obtieneFactCuad(a,grado,cuad,x);
    numr = raices(x,grado,real,imag);

    if (grado == 2)
    {
        float discriminante = a[1]*a[1] - 4*a[0]*a[2];

        if (discriminante == 0)
        {
            numr = 1;
        }
    }
    
    for (i=1;i<=grado;i++) {
        
        if (a[i]==0)
            cont=cont+1;
    }
    

    if (grado == cont)
        numr = 1;


    cout << endl << "=== Resultado ====" << endl;


    if(numr == 1)
    {
        
        if(grado==1)
        {
            cout << endl << "Ra\241ces (" << numr << " encontrada):" << endl << endl;
            cout << solucion << endl;
        }
        if(grado==2)
        {
            cout << endl << "Ra\241ces (" << numr << " encontrada):" << endl << endl;
            cout.setf(ios::showpoint|ios::showpos|ios::left|ios::scientific);
            cout.precision(5);
            cout << real[0] << endl;
        }
        else
        {
            cout << endl << "Ra\241ces (" << numr << " encontrada):" << endl << endl;
            cout << 0 << endl;
        }
    }
    else
    {
        cout << endl << "Ra\241ces (" << numr << " encontradas):" << endl << endl;
        cout << "*************************************************" << endl;
        cout << endl<< "R E A L E S:" << endl << endl;
    
        cout.setf(ios::showpoint|ios::showpos|ios::left|ios::scientific);
        cout.precision(5);

        for (i=0;i<grado;i++) {
            if (imag[i] == 0.0)
                cout << real[i] << endl;
        }

        cout << endl << "*************************************************" << endl;
        cout << endl << "I M A G I N A R I A S:" << endl << endl;

        for (i=0;i<grado;i++) {
            if (imag[i] != 0.0)
                cout << real[i] << " (" << imag[i] << ")i" << endl;
        }
    }

    cout << endl << "=== Integrantes ====" << endl << endl;
    cout << "Alex Bidart" << endl;
    cout << "Shu-yi Wong" << endl;
    cout << "Javier L\242pez" << endl;
    
    return 0;
}
