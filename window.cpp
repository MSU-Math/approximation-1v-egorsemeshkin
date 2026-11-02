#include "window.h"
#include "approximation.h"
#include <QPainter>
#include <QKeyEvent>
#include <cmath>
#include <algorithm>
#include <cstdio>

using namespace std;

Window::Window(QWidget *parent) : QWidget(parent) {
    a = -10.0;
    b = 10.0;
    n = 10;
    func_id = 0;
    display_mode = 0;
    zoom_s = 0;
    p_perturb = 0;
    setFocusPolicy(Qt::StrongFocus);
}

int Window::parse_command_line(int argc, char *argv[]) {
    if (argc != 5) return -1;
    if (sscanf(argv[1], "%lf", &a) != 1 ||
        sscanf(argv[2], "%lf", &b) != 1 ||
        sscanf(argv[3], "%d", &n) != 1 ||
        sscanf(argv[4], "%d", &func_id) != 1) {
        return -2;
    }
    if (b - a < 1.e-6 || n < 3 || func_id < 0 || func_id > 6) return -3;
    update_data();
    return 0;
}

double Window::f_exact(double x) const {
    switch (func_id) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return x * x;
        case 3: return x * x * x;
        case 4: return x * x * x * x;
        case 5: return exp(x);
        case 6: return 1.0 / (25.0 * x * x + 1.0);
        default: return 0.0;
    }
}

double Window::d2f_exact(double x) const {
    switch (func_id) {
        case 0: return 0.0;
        case 1: return 0.0;
        case 2: return 2.0;
        case 3: return 6.0 * x;
        case 4: return 12.0 * x * x;
        case 5: return exp(x);
        case 6: {
            double den = 25.0 * x * x + 1.0;
            return 50.0 * (75.0 * x * x - 1.0) / (den * den * den);
        }
        default: return 0.0;
    }
}

void Window::update_data() {
    switch (func_id) {
        case 0: f_name = "f(x) = 1"; break;
        case 1: f_name = "f(x) = x"; break;
        case 2: f_name = "f(x) = x^2"; break;
        case 3: f_name = "f(x) = x^3"; break;
        case 4: f_name = "f(x) = x^4"; break;
        case 5: f_name = "f(x) = e^x"; break;
        case 6: f_name = "f(x) = 1/(25x^2 + 1)"; break;
    }

    x_arr.assign(n, 0.0);
    f_arr.assign(n, 0.0);
    
    double max_f = 0.0;
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i) {
        x_arr[i] = a + i * step;
        f_arr[i] = f_exact(x_arr[i]);
        if (fabs(f_arr[i]) > max_f) max_f = fabs(f_arr[i]);
    }

    if (p_perturb != 0) {
        f_arr[n / 2] += p_perturb * 0.1 * max_f;
    }

    a_bessel.assign(4 * (n - 1), 0.0);
    extra_bessel.assign(n, 0.0);
    bessel_init(n, &x_arr[0], &f_arr[0], &a_bessel[0], &extra_bessel[0], d2f_exact(a), d2f_exact(b));

    a_spline.assign(4 * (n - 1), 0.0);
    extra_spline.assign(2 * n, 0.0);
    spline_init(n, &x_arr[0], &f_arr[0], &a_spline[0], &extra_spline.data()[0], d2f_exact(a), d2f_exact(b));
}

void Window::keyPressEvent(QKeyEvent *event) {
    switch (event->key()) {
        case Qt::Key_0:
            func_id = (func_id + 1) % 7;
            update_data();
            break;
        case Qt::Key_1:
            display_mode = (display_mode + 1) % 4;
            break;
        case Qt::Key_2:
            zoom_s++;
            break;
        case Qt::Key_3:
            zoom_s--;
            break;
        case Qt::Key_4:
            n *= 2;
            update_data();
            break;
        case Qt::Key_5:
            if (n / 2 >= 3) {
                n /= 2;
                update_data();
            }
            break;
        case Qt::Key_6:
            p_perturb++;
            update_data();
            break;
        case Qt::Key_7:
            p_perturb--;
            update_data();
            break;
        default:
            QWidget::keyPressEvent(event);
            return;
    }
    update();
}

void Window::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    
    double scale_pow = pow(2.0, zoom_s);
    double mid = (a + b) / 2.0;
    double half_width = (b - a) / (2.0 * scale_pow);
    double draw_a = mid - half_width;
    double draw_b = mid + half_width;

    int pixels = width();
    
    ApproxFunction f_bessel(a, b, n, &x_arr[0], &a_bessel[0], bessel_compute);
    ApproxFunction f_spline(a, b, n, &x_arr[0], &a_spline[0], spline_compute);

    vector<double> val_exact(pixels);
    vector<double> val_bessel(pixels);
    vector<double> val_spline(pixels);

    double min_y = 1e300, max_y = -1e300;
    double max_err_bessel = 0.0;
    double max_err_spline = 0.0;

    for (int px = 0; px < pixels; ++px) {
        double math_x = draw_a + px * (draw_b - draw_a) / (pixels - 1);
        val_exact[px] = f_exact(math_x);
        val_spline[px] = f_spline.compute(math_x);

        double err_s = fabs(val_exact[px] - val_spline[px]);
        if (err_s > max_err_spline) max_err_spline = err_s;

        double err_b = 0.0;
        if (n <= 50) {
            val_bessel[px] = f_bessel.compute(math_x);
            err_b = fabs(val_exact[px] - val_bessel[px]);
            if (err_b > max_err_bessel) max_err_bessel = err_b;
        }

        if (display_mode == 0) {
            min_y = min(min_y, val_exact[px]);
            max_y = max(max_y, val_exact[px]);
            if (n <= 50) {
                min_y = min(min_y, val_bessel[px]);
                max_y = max(max_y, val_bessel[px]);
            }
        } else if (display_mode == 1) {
            min_y = min(min_y, min(val_exact[px], val_spline[px]));
            max_y = max(max_y, max(val_exact[px], val_spline[px]));
        } else if (display_mode == 2) {
            min_y = min(min_y, min(val_exact[px], val_spline[px]));
            max_y = max(max_y, max(val_exact[px], val_spline[px]));
            if (n <= 50) {
                min_y = min(min_y, val_bessel[px]);
                max_y = max(max_y, val_bessel[px]);
            }
        } else if (display_mode == 3) {
            min_y = min(min_y, err_s);
            max_y = max(max_y, err_s);
            if (n <= 50) {
                min_y = min(min_y, err_b);
                max_y = max(max_y, err_b);
            }
        }
    }

    double true_max_abs = max(fabs(min_y), fabs(max_y));

    if (max_y - min_y < 1e-12) {
        double pad = (true_max_abs > 1e-12) ? 0.1 * true_max_abs : 1e-6;
        max_y += pad;
        min_y -= pad;
    }

    double delta_y = 0.05 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    painter.save();
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (draw_b - draw_a), -height() / (max_y - min_y));
    painter.translate(-0.5 * (draw_a + draw_b), -0.5 * (min_y + max_y));

    QPen pen_axis(Qt::red);
    pen_axis.setWidth(0);
    painter.setPen(pen_axis);
    painter.drawLine(QPointF(draw_a, 0), QPointF(draw_b, 0));
    painter.drawLine(QPointF(0, max_y), QPointF(0, min_y));

    QPen pen_exact(Qt::black);
    pen_exact.setWidth(0);
    QPen pen_bessel(Qt::blue);
    pen_bessel.setWidth(0);
    QPen pen_spline(Qt::green);
    pen_spline.setWidth(0);

    for (int px = 0; px < pixels - 1; ++px) {
        double x1 = draw_a + px * (draw_b - draw_a) / (pixels - 1);
        double x2 = draw_a + (px + 1) * (draw_b - draw_a) / (pixels - 1);

        if (display_mode != 3) {
            painter.setPen(pen_exact);
            painter.drawLine(QPointF(x1, val_exact[px]), QPointF(x2, val_exact[px+1]));
        }

        if (display_mode == 0 || display_mode == 2) {
            if (n <= 50) {
                painter.setPen(pen_bessel);
                painter.drawLine(QPointF(x1, val_bessel[px]), QPointF(x2, val_bessel[px+1]));
            }
        }

        if (display_mode == 1 || display_mode == 2) {
            painter.setPen(pen_spline);
            painter.drawLine(QPointF(x1, val_spline[px]), QPointF(x2, val_spline[px+1]));
        }

        if (display_mode == 3) {
            if (n <= 50) {
                painter.setPen(pen_bessel);
                painter.drawLine(QPointF(x1, fabs(val_exact[px] - val_bessel[px])), 
                                 QPointF(x2, fabs(val_exact[px+1] - val_bessel[px+1])));
            }
            painter.setPen(pen_spline);
            painter.drawLine(QPointF(x1, fabs(val_exact[px] - val_spline[px])), 
                             QPointF(x2, fabs(val_exact[px+1] - val_spline[px+1])));
        }
    }

    painter.restore();

    painter.setPen(Qt::black);
    char info[256];
    if (n <= 50) {
        snprintf(info, sizeof(info), "k=%d %s | n=%d | zoom=%d | p=%d | max|F|=%.4e | Err_B=%.4e | Err_S=%.4e",
                 func_id, f_name.c_str(), n, zoom_s, p_perturb, true_max_abs, max_err_bessel, max_err_spline);
    } else {
        snprintf(info, sizeof(info), "k=%d %s | n=%d | zoom=%d | p=%d | max|F|=%.4e | Err_S=%.4e",
                 func_id, f_name.c_str(), n, zoom_s, p_perturb, true_max_abs, max_err_spline);
    }
    painter.drawText(10, 20, info);
    
    const char *mode_str = "";
    if (display_mode == 0) mode_str = "Mode: Exact + Bessel";
    else if (display_mode == 1) mode_str = "Mode: Exact + Spline";
    else if (display_mode == 2) mode_str = "Mode: Exact + Bessel + Spline";
    else if (display_mode == 3) mode_str = "Mode: Errors";
    painter.drawText(10, 40, mode_str);

    static int last_func = -1, last_mode = -1, last_n = -1, last_zoom = -999, last_p = -999;
    if (func_id != last_func || display_mode != last_mode || n != last_n || zoom_s != last_zoom || p_perturb != last_p) {
        printf("%s | %s\n", info, mode_str);
        fflush(stdout);
        last_func = func_id; last_mode = display_mode; last_n = n; last_zoom = zoom_s; last_p = p_perturb;
    }
}
