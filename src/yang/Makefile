# Build makefile for rubens (Original makefile = JOE makefile by Vincent Terrapon)
# Author: Ik Jang 
# Version: 1.0
# Date: 09/2015

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep
include $(MUM_HOME)/src/ike/Makefile.dep

#
# 
#
.PHONY: default help clean

default: rubens

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'joe    :  compiles joe (RANS solver)'
	@echo 'rubens :  compiles rubens (QoI-recalculation solver)'
	@echo 'clean  :  clean rubens (included common files!)'
	@echo ' '

##########
# JOE
##########

COMMON = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv2.o Logging.o
JOE = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o
JOE_OBJS:  $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE)) 

JOE_LIBS = $(CLIBS) -lm

# COMMON_OBJS: For common, you can simply make it
COMMON_OBJS:  
	$(MAKE) $(COMMON) common_comb -C $(MUM_HOME)/src/common 

$(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o: $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp -o $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o

$(MUM_HOME)/src/ike/JOE/Scalars.o: $(MUM_HOME)/src/ike/JOE/Scalars.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/Scalars.cpp -o $(MUM_HOME)/src/ike/JOE/Scalars.o

$(MUM_HOME)/src/ike/JOE/JoeWithModels.o: $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp $(JoeWithModels_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp -o $(MUM_HOME)/src/ike/JOE/JoeWithModels.o


##########
# YANG (YANG requires JOE)
##########
JOE_FOR_YANG = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o 

YANG_OBJS = yangWithModels.o

yang: COMMON_OBJS $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_FOR_YANG)) $(addprefix $(MUM_HOME)/src/yang/, $(YANG_OBJS)) yang.o 
	$(CC) -o $@ $(CLD) $(addprefix $(MUM_HOME)/src/common/, $(COMMON)) $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_FOR_YANG)) $(addprefix $(MUM_HOME)/src/yang/, $(YANG_OBJS)) yang.o $(JOE_LIBS)  

$(MUM_HOME)/src/yang/yangWithModels.o: $(MUM_HOME)/src/yang/yangWithModels.cpp $(MUM_HOME)/src/yang/yangWithModels.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/yang/yangWithModels.cpp -o $(MUM_HOME)/src/yang/yangWithModels.o

yang.o: yang.cpp $(MUM_HOME)/src/yang/yangWithModels.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -I$(MUM_HOME)/src/ike/JOE/turbModels -I$(MUM_HOME)/src/ike/JOE/combModels -I$(MUM_HOME)/src/yang -c $<


##########
# ALL
##########

clean:
	rm $(addprefix $(MUM_HOME)/src/common/, $(COMMON)) JOE_OBJS $(addprefix $(MUM_HOME)/src/yang/, $(YANG_OBJS)) yang.o joe yang

