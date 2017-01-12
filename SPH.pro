TARGET=FluidSim
CONFIG+=c++11
CONFIG-= x86_64 app_bundle

OBJECTS_DIR=obj
MOC_DIR=moc
QT+= opengl core gui

INCLUDEPATH+=./include \
            /usr/local/include

SOURCES+=$$PWD/src/SphMain.cpp \
         $$PWD/src/SphSystem.cpp \
         $$PWD/src/SphPhaseData.cpp \
         $$PWD/src/Vec3.cpp

HEADERS+=$$PWD/include/SphData.h \
         $$PWD/include/SphHeader.h \
         $$PWD/include/SphSystem.h \
         $$PWD/include/SphType.h \
         $$PWD/include/SphParticle.h \
         $$PWD/include/Vec3.h \
         $$PWD/include/SphPhaseData.h

OTHER_FILES+=$$PWD/Shader/shader.fs \
             $$PWD/Shader/shader.vs

DESTDIR=./

macx:LIBS+= -framework OpenGL -framework GLUT -L/usr/local/lib -lglew

linux: LIBS+= -lGLU -lGLEW -lglut
