################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Implementazioni/Animale.cpp \
../Implementazioni/Cane.cpp \
../Implementazioni/EssereVivente.cpp \
../Implementazioni/Gatto.cpp \
../Implementazioni/OperazioniVettore.cpp \
../Implementazioni/Umano.cpp 

CPP_DEPS += \
./Implementazioni/Animale.d \
./Implementazioni/Cane.d \
./Implementazioni/EssereVivente.d \
./Implementazioni/Gatto.d \
./Implementazioni/OperazioniVettore.d \
./Implementazioni/Umano.d 

OBJS += \
./Implementazioni/Animale.o \
./Implementazioni/Cane.o \
./Implementazioni/EssereVivente.o \
./Implementazioni/Gatto.o \
./Implementazioni/OperazioniVettore.o \
./Implementazioni/Umano.o 


# Each subdirectory must supply rules for building sources it contributes
Implementazioni/%.o: ../Implementazioni/%.cpp Implementazioni/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Implementazioni

clean-Implementazioni:
	-$(RM) ./Implementazioni/Animale.d ./Implementazioni/Animale.o ./Implementazioni/Cane.d ./Implementazioni/Cane.o ./Implementazioni/EssereVivente.d ./Implementazioni/EssereVivente.o ./Implementazioni/Gatto.d ./Implementazioni/Gatto.o ./Implementazioni/OperazioniVettore.d ./Implementazioni/OperazioniVettore.o ./Implementazioni/Umano.d ./Implementazioni/Umano.o

.PHONY: clean-Implementazioni

