QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2

INCLUDEPATH += "/usr/local/include/eigen3/"

HEADERS += \
    filereader.h \
    parameters.h \
    constants.h \
    matrixnumerov.h \
    generalizedmatrixnumerov.h \
    matrixreader.h

SOURCES += \
    filereader.cpp \
    parameters.cpp \
    main.cpp \
    matrixnumerov.cpp \
    generalizedmatrixnumerov.cpp \
    matrixreader.cpp

DISTFILES += \
    parameters.inp
