using namespace std;
typedef double Real;
Real norm_2(Real v[],int n);
void newton_sist(void(*effe)(Real*,Real*),void(*Jeffe)(Real*,Real*), int n, Real x[], int* nit, Real toll, int nitmax);
void newtonVERO_sist(void(*effe)(Real*,Real*),void(*Jeffe)(Real*,Real*), int n, Real x[], int *nit, Real toll, int nitmax);
