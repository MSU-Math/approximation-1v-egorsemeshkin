#include <QApplication>
#include <QMainWindow>
#include "window.h"

using namespace std;

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QMainWindow *window = new QMainWindow;
    Window *graph_area = new Window(window);

    if (graph_area->parse_command_line(argc, argv) != 0) {
        return -1;
    }

    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");
    window->resize(1000, 800);
    window->show();

    return app.exec();
}
