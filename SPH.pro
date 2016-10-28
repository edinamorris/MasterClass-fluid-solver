TARGET=FluidSim
CONFIG+=c++11
CONFIG-= x86_64 app_bundle

OBJECTS_DIR=obj
MOC_DIR=moc
QT+= opengl core gui

INCLUDEPATH+=./include \
            /usr/local/include

SOURCES+=$$PWD/src/sph_main.cpp \
         $$PWD/src/sph_system.cpp

HEADERS+=$$PWD/include/sph_data.h \
         $$PWD/include/sph_header.h \
         $$PWD/include/sph_system.h \
         $$PWD/include/sph_type.h

OTHER_FILES+=$$PWD/Shader/shader.fs \
             $$PWD/Shader/shader.vs

DESTDIR=./

macx:LIBS+= -framework OpenGL -framework GLUT -L/usr/local/lib -lglew

linux: LIBS+= -lGLU -lGLEW
