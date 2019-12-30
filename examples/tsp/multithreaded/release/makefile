CC=g++
ARC=-m64 -mfpmath=sse -msse2
CFLAGS=-c -std=c++17 -fno-rtti -O3 -I$(HRD_MIPCL) $(ARC) -minline-all-stringops
LDFLAGS=-shared-libgcc
INS=install
MIP_DIR=$(MIPCLDIR)
LDFLAGS+= -L$(MIP_DIR)/lib -lmipcl -lpthread
CFLAGS+=-DMIP_API=
TARGET=tsp
#
PRJ_DIR=$(MIP_DIR)/examples/mipcl/tsp
#
HRD_MIPCL=$(MIP_DIR)/mipcl/headers
#
INSDIR=$(PRJ_DIR)/bin
SRC_PATH=$(PRJ_DIR)/sources
HRD_PATH=$(PRJ_DIR)/sources
#
vpath %.h $(HRD_PATH)
vpath %.cpp $(SRC_PATH)
SRC=tsp.cpp tspPool.cpp flowNet.cpp main.cpp
OBJS=$(SRC:.cpp=.o)
.cpp.o:
	$(CC) $(CFLAGS)  $<
build: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(TARGET)
#dependencies
flowNet.o: flowNet.cpp flowNet.h compare.h
tspPool.o: tspPool.cpp tspPool.h
tsp.o: tsp.cpp tsp.h flowNet.h tspPool.h
main.o: main.cpp tsp.h
#
.PHONY: clean
clean:
	rm -f $(OBJS)
#
install: $(TARGET)
	$(INS) $(TARGET) $(INSDIR)
	rm $(TARGET)
all: build install
