################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: sharpSAT

#g++ -L/opt/apps/gmp/5.0.4/lib/ -L/usr/lib/ -I/opt/apps/gmp/5.0.4/include/ -g -pg -o "sharpSAT" $(OBJS) $(USER_OBJS) $(LIBS)
# Tool invocations
sharpSAT: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L/home/shubham/Downloads/gmp-6.1.2/ -I/home/shubham/Downloads/gmp-6.1.2/ -L/usr/lib/ -L/usr/local/lib -g -pg -o "sharpSAT" $(OBJS) $(USER_OBJS) $(LIBS)
	# g++ -L/home/kgm2/Experimentation/Scripts/INSTALL/gmp/lib/ -I/home/kgm2/Experimentation/Scripts/INSTALL/gmp/include/ -L/usr/lib/ -g -pg -o "sharpSAT" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) sharpSAT
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
