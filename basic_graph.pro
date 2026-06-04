QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h approximation.h
SOURCES       = main.cpp \
                window.cpp \
                bessel.cpp \
                spline.cpp
QT += widgets core gui
