TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    gauleg.cpp \
    gauss_laguerre.cpp \
    ran0.cpp

LIBS += -llapack -lblas -larmadillo -fopenmp



