#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <vector>
#include <string>

using namespace std;

class Window : public QWidget {
    Q_OBJECT

private:
    int func_id;
    string f_name;
    double a;
    double b;
    int n;
    int display_mode;
    int zoom_s;
    int p_perturb;

    vector<double> x_arr;
    vector<double> f_arr;
    vector<double> a_bessel;
    vector<double> a_spline;
    vector<double> extra_bessel;
    vector<double> extra_spline;

    double f_exact(double x) const;
    double d2f_exact(double x) const;
    void update_data();

public:
    Window(QWidget *parent = nullptr);
    int parse_command_line(int argc, char *argv[]);

protected:
    void paintEvent(QPaintEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
};

#endif
