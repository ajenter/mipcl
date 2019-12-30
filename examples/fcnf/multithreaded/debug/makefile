CC=g++
ARC=-m64 -mfpmath=sse -msse2
CFLAGS=-c -std=c++17 -fno-rtti -g -I$(HRD_MIPCL) $(ARC)
INS=install
#
MIP_DIR=$(MIPCLDIR)
LDFLAGS= -L$(MIP_DIR)/lib -lmipcl_dbg
CFLAGS+=-DMIP_API=
TARGET=fcnf_dbg
#
PRJ_DIR=$(MIP_DIR)/examples/mipcl/fcnf
#
HRD_MIPCL=$(MIP_DIR)/mipcl/headers
#
INSDIR=$(PRJ_DIR)/bin
SRC_PATH=$(PRJ_DIR)/sources
HRD_PATH=$(PRJ_DIR)/sources
#
vpath %.h $(HRD_PATH)
vpath %.cpp $(SRC_PATH)
SRC=fcnf.cpp main.cpp
OBJS=$(SRC:.cpp=.o)
.cpp.o:
	$(CC) $(CFLAGS)  $<
build: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)
#dependencies
fcnf.o: fcnf.cpp fcnf.h
main.o: main.cpp fcnf.h
#
.PHONY: clean
clean:
	rm -f $(OBJS)
#
install: $(TARGET)
	$(INS) $(TARGET) $(INSDIR)
	rm $(TARGET)
all: build install
