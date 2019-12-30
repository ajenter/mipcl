CC=g++
ARC=-m64 -mfpmath=sse -msse2
CFLAGS=-c -std=c++17 -fno-rtti -O3 -I$(HRD_MIPCL) $(ARC) -minline-all-stringops -D__ONE_THREAD_
INS=install
MIP_DIR=$(MIPCLDIR)
LDFLAGS+= -L$(MIP_DIR)/lib -lmipcl-1
CFLAGS+= -DMIP_API=
TARGET=genassign-1
#
PRJ_DIR=$(MIP_DIR)/examples/mipcl/genAssign
#
HRD_MIPCL=$(MIP_DIR)/mipcl/headers
#
INSDIR=$(PRJ_DIR)/bin
SRC_PATH=$(PRJ_DIR)/sources
HRD_PATH=$(PRJ_DIR)/sources
#
vpath %.h $(HRD_PATH)
vpath %.cpp $(SRC_PATH)
SRC=genAssign.cpp knapsack.cpp main.cpp
OBJS=$(SRC:.cpp=.o)
.cpp.o:
	$(CC) $(CFLAGS)  $<	
build: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(TARGET)
#dependencies
genAssign.o: genAssign.cpp genAssign.h knapsack.h
knapsack.o: knapsack.cpp knapsack.h
main.o: main.cpp genAssign.h
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
