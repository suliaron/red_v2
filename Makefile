# http://scottmcpeak.com/autodepend/autodepend.html
# https://www.gnu.org/software/make/manual/make.html#Rules
# http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#axzz3ODhDRl4y
# http://owen.sj.ca.us/~rk/howto/slides/make/slides/makerecurs.html

CUDA_PATH := /usr/local/cuda
NVCC_PATH := /usr/local/cuda/bin
# NVCC Compiler
NVCC := $(NVCC_PATH)/nvcc
# Linker
LINK := $(NVCC_PATH)/nvcc
#
RM = rm -f

# Options for the nvcc compiler
# NVCC_FLAGS := -Xcompiler -Wall -G -g -O0 -gencode arch=compute_20,code=sm_20 -fmad=false
NVCC_FLAGS := -Xcompiler -Wall -O2 -gencode arch=compute_30,code=sm_30 -fmad=false

# Paths
SRC             := src
RED             := $(SRC)/red_v2
RED_INITIAL     := $(SRC)/red_v2_initial
RED_UTILITY     := $(SRC)/red_v2_utility
RED_BENCHMARK   := $(SRC)/benchmark

INCLUDES        := -Ibase_type -I$(RED_UTILITY) -I$(RED_INITIAL) -I$(CUDA_PATH)/include

RED_LIBRARY     := red.library
BIN             := bin

# List the objects for the executables and library
RED_OBJS := \
integrator.o \
int_euler.o \
int_hermite4.o \
int_rungekutta2.o \
int_rungekutta4.o \
int_rungekutta5.o \
int_rungekutta7.o \
main.o \
nbody.o \
ode.o \
options.o \
parameter.o \
rtbp1D.o \
rtbp2D.o \
rtbp3D.o \
tbp1D.o \
tbp2D.o \
tbp3D.o \
threebody.o

RED_INITIAL_OBJS := \
distribution.o \
main.o \
util_init.o

RED_UTILITY_OBJS := \
file_util.o \
tokenizer.o \
tools.o \
util.o

RED_BENCHMARK_OBJS := \
cpu_base.o \
main.o \
util.o

RED_DEPS := $(RED_OBJS:.o=.d)
RED_INITIAL_DEPS := $(RED_INITIAL_OBJS:.o=.d)
RED_UTILITY_DEPS := $(RED_UTILITY_OBJS:.o=.d)
RED_BENCHMARK_DEPS := $(RED_BENCHMARK_OBJS:.o=.d)

# Targets
all : redutil2 red_v2 red_v2_init red_v2_benchmark

red_v2 : redutil2 $(RED)/red_v2 | $(BIN)

red_v2_init : redutil2 $(RED_INITIAL)/red_v2_init | $(BIN)

red_v2_benchmark : redutil2 $(RED_BENCHMARK)/red_v2_benchmark | $(BIN)

redutil2 : $(RED_UTILITY)/redutil2.a | $(RED_LIBRARY)

-include $(addprefix $(RED)/, $(RED_DEPS))
-include $(addprefix $(RED_INITIAL)/, $(RED_INITIAL_DEPS))
-include $(addprefix $(RED_UTILITY)/, $(RED_UTILITY_DEPS))
-include $(addprefix $(RED_BENCHMARK)/, $(RED_BENCHMARK_DEPS))

# Build rules
$(RED)/red_v2 : $(addprefix $(RED)/, $(RED_OBJS)) | $(BIN)
	$(LINK) $(RED_LIBRARY)/redutil2.a -o $@ $?
	cp $@ $(BIN)/
	
$(RED_INITIAL)/red_v2_init : $(addprefix $(RED_INITIAL)/, $(RED_INITIAL_OBJS)) | $(BIN)
	$(LINK) $(RED_LIBRARY)/redutil.a -o $@ $?
	cp $@ $(BIN)/

$(RED_BENCHMARK)/red_v2_benchmark : $(addprefix $(RED_BENCHMARK)/, $(RED_BENCHMARK_OBJS)) | $(BIN)
	$(LINK) $(RED_LIBRARY)/redutil.a -o $@ $?
	cp $@ $(BIN)/
	
$(RED_UTILITY)/redutil2.a : $(addprefix $(RED_UTILITY)/, $(RED_UTILITY_OBJS)) | $(RED_LIBRARY)
	ar cr $@ $?
	cp $@ $(RED_LIBRARY)/

# compile and generate dependency info
$(RED)/%.o : $(RED)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED)/%.o : $(RED)/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_INITIAL)/%.o : $(RED_INITIAL)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_INITIAL)/%.o : $(RED_INITIAL)/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_BENCHMARK)/%.o : $(RED_BENCHMARK)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_BENCHMARK)/%.o : $(RED_BENCHMARK)/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_UTILITY)/%.o : $(RED_UTILITY)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

# compile and generate dependency info
$(RED_UTILITY)/%.o : $(RED_UTILITY)/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking $(NVCC) Compiler'
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -M -odir "" $(NVCC_FLAGS) $(INCLUDES) -o "$(@:%.o=%.d)" $<
	@echo 'Finished building: $<'
	@echo ''

$(RED_LIBRARY) :
	mkdir $(RED_LIBRARY)

$(BIN) : 
	mkdir $(BIN)

clean:
	-$(RM) $(RED_LIBRARY)/*.a $(RED_UTILITY)/*.o $(RED_UTILITY)/*.d $(RED_INITIAL)/*.o $(RED_INITIAL)/*.d  $(RED_BENCHMARK)/*.o $(RED_BENCHMARK)/*.d $(RED)/*.o $(RED)/*.d $(BIN)/*
