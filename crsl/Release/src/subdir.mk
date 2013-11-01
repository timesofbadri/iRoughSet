################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CRS2.cpp \
../src/IncCRS2.cpp \
../src/crsl.cpp \
../src/parse_args.cpp 

OBJS += \
./src/CRS2.o \
./src/IncCRS2.o \
./src/crsl.o \
./src/parse_args.o 

CPP_DEPS += \
./src/CRS2.d \
./src/IncCRS2.d \
./src/crsl.d \
./src/parse_args.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


