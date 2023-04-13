QT -= core
QT -= gui

QMAKE_CXXFLAGS -= -std=gnu++98
QMAKE_CXXFLAGS -= -std=gnu++0x
QMAKE_CXXFLAGS -= --std=c++11
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++98
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++0x
QMAKE_CXXFLAGS_RELEASE -= --std=c++11

QMAKE_CXXFLAGS += -std=c++14

TARGET = aacp_template
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

#LIBS += -L/usr/lib64 -lgslcblas
LIBS += -lblas -llapack

SOURCES += main.cpp

HEADERS += \
    la_wrapper.h \
    la_base_obj.h \
    la_operations.h

