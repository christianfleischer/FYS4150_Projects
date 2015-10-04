TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    initialize_electrons.cpp \
    jacobi_rotations.cpp

LIBS += -llapack -lblas -larmadillo

HEADERS += \
    catch.h

