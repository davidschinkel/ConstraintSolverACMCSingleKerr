EXECUTABLES := ACMC_Data
LIBS:= mpfr quadmath
#PREFIXES:=$(HOME)/local
BUILD_DIR=build/
#SRC_DIR=src/
DEP_DIR=dep/
COMMON_OBJS:= allocate_and_free_derivs BondiMass FieldAndBoundEqns FuncAndJacobian
COMMON_OBJS:= $(COMMON_OBJS) Get_Coefficients GetNewResolution IndexRoutines JacobianFD
COMMON_OBJS:= $(COMMON_OBJS) newton ScanAndPrint utilities cast

CXX:=g++
INC:=$(patsubst %,-I%/include,$(PREFIXES))
LIB:=$(patsubst %,-I%/lib,$(PREFIXES))

CXXFLAGS:=-g -MMD $(INC) -O3 -Wextra
CPPFLAGS:=-std=c++0x
CFLAGS:=
LDFLAGS:=-g -fopenmp $(LIB)

LDLIBS:=$(patsubst %,-l%,$(LIBS))

SRCS:= $(wildcard *.c)

DEPS := $(patsubst %.c,$(DEP_DIR)%.d,$(SRCS))

.SUFFIXES:

all: $(patsubst %,%,$(EXECUTABLES))

ACMC_Data: $(BUILD_DIR)main.o $(patsubst %,$(BUILD_DIR)%.o,$(COMMON_OBJS))
	$(CXX) $(LDFLAGS) $< $(patsubst %,$(BUILD_DIR)%.o,$(COMMON_OBJS)) $(LDLIBS) -o $@
		notify-send -t 10000 "Information" "Compiler fertig."

$(BUILD_DIR)%.o: %.c
	mkdir -p dep build ../run/SequenceElements ../run/BondiMass
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) $< -MF $(patsubst %.c,$(DEP_DIR)%.d,$<) -o $@

-include $(DEPS)

clean:
	rm -f -r $(BUILD_DIR) $(patsubst %,%,$(EXECUTABLES))

veryclean: clean
	rm -f $(DEPS)

	