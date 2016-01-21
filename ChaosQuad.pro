#-------------------------------------------------
#
# Project created by QtCreator 2014-06-21T11:18:37
#
#-------------------------------------------------
QMAKE_CXXFLAGS += -std=gnu++11 #\
                  #-O3

LIBS += -lquadmath

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = chaosDetector2
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp \
    odesolver.cpp \
    childsolver.cpp \
    qcustomplot.cpp \
    multithreadsolver.cpp \
    task.cpp \
    declarations.cpp


HEADERS  += mainwindow.h \
    odesolver.h \
    childsolver.h \
    qcustomplot.h \
    multithreadsolver.h \
    task.h \
    declarations.h


FORMS    += mainwindow.ui \
    plotwindow.ui
