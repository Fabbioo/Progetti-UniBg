################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BusinessLogic/Animale.cpp \
../BusinessLogic/Cane.cpp \
../BusinessLogic/EssereVivente.cpp \
../BusinessLogic/Gatto.cpp \
../BusinessLogic/OperazioniVettore.cpp \
../BusinessLogic/Umano.cpp 

CPP_DEPS += \
./BusinessLogic/Animale.d \
./BusinessLogic/Cane.d \
./BusinessLogic/EssereVivente.d \
./BusinessLogic/Gatto.d \
./BusinessLogic/OperazioniVettore.d \
./BusinessLogic/Umano.d 

OBJS += \
./BusinessLogic/Animale.o \
./BusinessLogic/Cane.o \
./BusinessLogic/EssereVivente.o \
./BusinessLogic/Gatto.o \
./BusinessLogic/OperazioniVettore.o \
./BusinessLogic/Umano.o 


# Each subdirectory must supply rules for building sources it contributes
BusinessLogic/%.o: ../BusinessLogic/%.cpp BusinessLogic/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-BusinessLogic

clean-BusinessLogic:
	-$(RM) ./BusinessLogic/Animale.d ./BusinessLogic/Animale.o ./BusinessLogic/Cane.d ./BusinessLogic/Cane.o ./BusinessLogic/EssereVivente.d ./BusinessLogic/EssereVivente.o ./BusinessLogic/Gatto.d ./BusinessLogic/Gatto.o ./BusinessLogic/OperazioniVettore.d ./BusinessLogic/OperazioniVettore.o ./BusinessLogic/Umano.d ./BusinessLogic/Umano.o

.PHONY: clean-BusinessLogic

