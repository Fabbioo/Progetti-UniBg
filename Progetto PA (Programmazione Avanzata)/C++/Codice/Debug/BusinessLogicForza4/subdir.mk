################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BusinessLogicForza4/Forza4.cpp \
../BusinessLogicForza4/Generatore.cpp \
../BusinessLogicForza4/Gettone.cpp \
../BusinessLogicForza4/GettoneRosso.cpp \
../BusinessLogicForza4/GettoneVerde.cpp \
../BusinessLogicForza4/Griglia.cpp 

CPP_DEPS += \
./BusinessLogicForza4/Forza4.d \
./BusinessLogicForza4/Generatore.d \
./BusinessLogicForza4/Gettone.d \
./BusinessLogicForza4/GettoneRosso.d \
./BusinessLogicForza4/GettoneVerde.d \
./BusinessLogicForza4/Griglia.d 

OBJS += \
./BusinessLogicForza4/Forza4.o \
./BusinessLogicForza4/Generatore.o \
./BusinessLogicForza4/Gettone.o \
./BusinessLogicForza4/GettoneRosso.o \
./BusinessLogicForza4/GettoneVerde.o \
./BusinessLogicForza4/Griglia.o 


# Each subdirectory must supply rules for building sources it contributes
BusinessLogicForza4/%.o: ../BusinessLogicForza4/%.cpp BusinessLogicForza4/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-BusinessLogicForza4

clean-BusinessLogicForza4:
	-$(RM) ./BusinessLogicForza4/Forza4.d ./BusinessLogicForza4/Forza4.o ./BusinessLogicForza4/Generatore.d ./BusinessLogicForza4/Generatore.o ./BusinessLogicForza4/Gettone.d ./BusinessLogicForza4/Gettone.o ./BusinessLogicForza4/GettoneRosso.d ./BusinessLogicForza4/GettoneRosso.o ./BusinessLogicForza4/GettoneVerde.d ./BusinessLogicForza4/GettoneVerde.o ./BusinessLogicForza4/Griglia.d ./BusinessLogicForza4/Griglia.o

.PHONY: clean-BusinessLogicForza4

