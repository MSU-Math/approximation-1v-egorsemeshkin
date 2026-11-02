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

int bessel_init(int n, double *x, double *f, double *a, double *extra, double d2f_left, double d2f_right) {
    if (n < 3) return -1;
    double *df = extra;
    
    for (int i = 1; i < n - 1; ++i) {
        double h_prev = x[i] - x[i-1];
        double h_next = x[i+1] - x[i];
        df[i] = (h_prev / (h_prev + h_next)) * ((f[i+1] - f[i]) / h_next) +
                (h_next / (h_prev + h_next)) * ((f[i] - f[i-1]) / h_prev);
    }
    
    double h0 = x[1] - x[0];
    df[0] = -(h0 / 4.0) * d2f_left + (3.0 / (2.0 * h0)) * (f[1] - f[0]) - 0.5 * df[1];
    
    double hn = x[n-1] - x[n-2];
    df[n-1] = (hn / 4.0) * d2f_right + (3.0 / (2.0 * hn)) * (f[n-1] - f[n-2]) - 0.5 * df[n-2];
    
    for (int i = 0; i < n - 1; ++i) {
        double h = x[i+1] - x[i];
        a[4*i + 0] = f[i];
        a[4*i + 1] = df[i];
        a[4*i + 2] = (3.0 * (f[i+1] - f[i]) / (h * h)) - (2.0 * df[i] + df[i+1]) / h;
        a[4*i + 3] = (2.0 * (f[i] - f[i+1]) / (h * h * h)) + (df[i] + df[i+1]) / (h * h);
    }
    
    return 0;
}

double bessel_compute(double x_eval, double a_bnd, double b_bnd, int n, double *x, double *a) {
    (void)a_bnd;
    (void)b_bnd;
    int i = find_interval(x_eval, n, x);
    double dx = x_eval - x[i];
    return a[4*i + 0] + dx * (a[4*i + 1] + dx * (a[4*i + 2] + dx * a[4*i + 3]));
}
