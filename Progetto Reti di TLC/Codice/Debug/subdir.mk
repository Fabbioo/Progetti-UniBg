################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../buffer.cpp \
../calendar.cpp \
../easyio.cpp \
../event.cpp \
../main.cpp \
../packet.cpp \
../queue.cpp \
../rand.cpp \
../simulator.cpp \
../stat.cpp 

CPP_DEPS += \
./buffer.d \
./calendar.d \
./easyio.d \
./event.d \
./main.d \
./packet.d \
./queue.d \
./rand.d \
./simulator.d \
./stat.d 

OBJS += \
./buffer.o \
./calendar.o \
./easyio.o \
./event.o \
./main.o \
./packet.o \
./queue.o \
./rand.o \
./simulator.o \
./stat.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean--2e-

clean--2e-:
	-$(RM) ./buffer.d ./buffer.o ./calendar.d ./calendar.o ./easyio.d ./easyio.o ./event.d ./event.o ./main.d ./main.o ./packet.d ./packet.o ./queue.d ./queue.o ./rand.d ./rand.o ./simulator.d ./simulator.o ./stat.d ./stat.o

.PHONY: clean--2e-

