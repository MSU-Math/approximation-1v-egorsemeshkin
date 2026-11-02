#include "approximation.h"

using namespace std;

static int find_interval(double x_eval, int n, double *x) {
    if (x_eval <= x[0]) return 0;
    if (x_eval >= x[n - 1]) return n - 2;
    int l = 0;
    int r = n - 1;
    while (r - l > 1) {
        int mid = l + (r - l) / 2;
        if (x_eval >= x[mid]) l = mid;
        else r = mid;
    }
    return l;
}

int spline_init(int n, double *x, double *f, double *a, double *extra, double d2f_left, double d2f_right) {
    if (n < 3) return -1;
    double *alpha = extra;
    double *M = extra + n;
    
    M[0] = d2f_left;
    M[n-1] = d2f_right;
    
    double h0 = x[1] - x[0];
    double h1 = x[2] - x[1];
    
    double A1 = h0 / 6.0;
    double B1 = (h0 + h1) / 3.0;
    double C1 = h1 / 6.0;
    double F1 = ((f[2] - f[1]) / h1) - ((f[1] - f[0]) / h0) - A1 * M[0];
    
    alpha[1] = -C1 / B1;
    M[1] = F1 / B1;
    
    for (int i = 2; i < n - 1; ++i) {
        double hi_prev = x[i] - x[i-1];
        double hi = x[i+1] - x[i];
        
        double Ai = hi_prev / 6.0;
        double Bi = (hi_prev + hi) / 3.0;
        double Ci = hi / 6.0;
        double Fi = ((f[i+1] - f[i]) / hi) - ((f[i] - f[i-1]) / hi_prev);
        
        if (i == n - 2) {
            Fi -= Ci * M[n-1];
            Ci = 0.0;
        }
        
        double denom = Bi + Ai * alpha[i-1];
        alpha[i] = -Ci / denom;
        M[i] = (Fi - Ai * M[i-1]) / denom;
    }
    
    for (int i = n - 3; i >= 1; --i) {
        M[i] = alpha[i] * M[i+1] + M[i];
    }
    
    for (int i = 0; i < n - 1; ++i) {
        double h = x[i+1] - x[i];
        a[4*i + 0] = f[i];
        a[4*i + 1] = (f[i+1] - f[i]) / h - h * (2.0 * M[i] + M[i+1]) / 6.0;
        a[4*i + 2] = M[i] / 2.0;
        a[4*i + 3] = (M[i+1] - M[i]) / (6.0 * h);
    }
    
    return 0;
}

double spline_compute(double x_eval, double a_bnd, double b_bnd, int n, double *x, double *a) {
    (void)a_bnd;
    (void)b_bnd;
    int i = find_interval(x_eval, n, x);
    double dx = x_eval - x[i];
    return a[4*i + 0] + dx * (a[4*i + 1] + dx * (a[4*i + 2] + dx * a[4*i + 3]));
}
