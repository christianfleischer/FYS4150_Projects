TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    tridiag.cpp \
    ran2.cpp \
    gaussian_deviate.cpp

LIBS += -llapack -lblas -larmadillo


