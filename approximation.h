#ifndef APPROXIMATION_H
#define APPROXIMATION_H

int bessel_init(int n, double *x, double *f, double *a, double *extra, double d2f_left, double d2f_right);
double bessel_compute(double x_eval, double a_bnd, double b_bnd, int n, double *x, double *a);

int spline_init(int n, double *x, double *f, double *a, double *extra, double d2f_left, double d2f_right);
double spline_compute(double x_eval, double a_bnd, double b_bnd, int n, double *x, double *a);

class IFunction {
public:
    virtual double compute(double x) const = 0;
    virtual ~IFunction() = default;
};

class ApproxFunction : public IFunction {
private:
    double a_bnd;
    double b_bnd;
    int n;
    double *x_array;
    double *a_array;
    typedef double (*ComputeMethod)(double, double, double, int, double*, double*);
    ComputeMethod method;

public:
    ApproxFunction(double a, double b, int n_points, double *x, double *coeff, ComputeMethod comp_method)
        : a_bnd(a), b_bnd(b), n(n_points), x_array(x), a_array(coeff), method(comp_method) {}

    double compute(double x) const override {
        return method(x, a_bnd, b_bnd, n, x_array, a_array);
    }
};

#endif
