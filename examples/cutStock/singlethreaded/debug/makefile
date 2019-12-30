CC=g++
ARC=-m64 -mfpmath=sse -msse2
CFLAGS=-c -std=c++17 -fno-rtti -g -I$(HRD_MIPCL) $(ARC) -D__ONE_THREAD_
LDFLAGS=-shared-libgcc
INS=install
#
MIP_DIR=$(MIPCLDIR)
LDFLAGS+= -L$(MIP_DIR)/lib -lmipcl-1_dbg
CFLAGS+=-DMIP_API=
TARGET=cutStock-1_dbg
#
PRJ_DIR=$(MIP_DIR)/examples/mipcl/cutStock
#
HRD_MIPCL=$(MIP_DIR)/mipcl/headers
#
INSDIR=$(PRJ_DIR)/bin
SRC_PATH=$(PRJ_DIR)/sources
HRD_PATH=$(PRJ_DIR)/sources
#
vpath %.h $(HRD_PATH)
vpath %.cpp $(SRC_PATH)
SRC=cutstock.cpp main.cpp
OBJS=$(SRC:.cpp=.o)
.cpp.o:
	$(CC) $(CFLAGS)  $<	
build: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(TARGET)
#dependencies
cutstock.o: cutstock.cpp cutstock.h
main.o: main.cpp cutstock.h
#
.PHONY: clean
clean:
	rm -f $(OBJS)
#
install: $(TARGET)
	$(INS) $(TARGET) $(INSDIR)
	rm $(TARGET)
all: build install
test:
	echo $(MIP_DIR)
