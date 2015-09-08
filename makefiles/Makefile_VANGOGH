# Build makefile for vanGogh (Original makefile = JOE makefile by Vincent Terrapon)
# Author: Ik Jang 
# Version: 1.0
# Date: 08/2015

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

default: vanGogh

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'joe :      compiles joe (RANS solver)'
	@echo 'vanGogh :  compiles vanGogh (Perturbations solver)'
	@echo 'clean :    clean vanGogh (included common files!)'
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

joe: $(addprefix $(MUM_HOME)/src/common/, $(COMMON)) JOE_OBJS joe.o
	$(CC) -o $@ $(CLD) $(OBJS) $(JOE_LIBS) #; mv $@ $(MUM_HOME)/bin

joe.o: joe.cpp $(joe_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<

$(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o: $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp -o $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o

$(MUM_HOME)/src/ike/JOE/Scalars.o: $(MUM_HOME)/src/ike/JOE/Scalars.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/Scalars.cpp -o $(MUM_HOME)/src/ike/JOE/Scalars.o

$(MUM_HOME)/src/ike/JOE/JoeWithModels.o: $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp $(JoeWithModels_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp -o $(MUM_HOME)/src/ike/JOE/JoeWithModels.o


##########
# VANGOGH (VANGOGH requires JOE)
##########
JOE_FOR_VANGOGH = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o 

VANGOGH_OBJS = DiffFilter.o vanGoghWithModels.o

vanGogh: COMMON_OBJS $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_FOR_VANGOGH)) $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) vanGogh.o 
	$(CC) -o $@ $(CLD) $(addprefix $(MUM_HOME)/src/common/, $(COMMON)) $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_FOR_VANGOGH)) $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) vanGogh.o $(JOE_LIBS)  

$(MUM_HOME)/src/vangogh/DiffFilter.o: $(MUM_HOME)/src/vangogh/DiffFilter.cpp $(MUM_HOME)/src/vangogh/DiffFilter.h $(UgpWithCV2_H) $(Param_H) $(myMem_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/vangogh/DiffFilter.cpp -o $(MUM_HOME)/src/vangogh/DiffFilter.o

$(MUM_HOME)/src/vangogh/vanGoghWithModels.o: $(MUM_HOME)/src/vangogh/vanGoghWithModels.cpp $(MUM_HOME)/src/vangogh/vanGoghWithModels.h $(MUM_HOME)/src/vangogh/vanGoghUtils.h $(MUM_HOME)/src/vangogh/DiffFilter.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/vangogh/vanGoghWithModels.cpp -o $(MUM_HOME)/src/vangogh/vanGoghWithModels.o

vanGogh.o: vanGogh.cpp $(MUM_HOME)/src/vangogh/vanGoghWithModels.h $(MUM_HOME)/src/vangogh/vanGoghUtils.h $(MUM_HOME)/src/vangogh/DiffFilter.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -I$(MUM_HOME)/src/ike/JOE/turbModels -I$(MUM_HOME)/src/ike/JOE/combModels -I$(MUM_HOME)/src/vangogh -c $<


##########
# ALL
##########

clean:
	rm $(addprefix $(MUM_HOME)/src/common/, $(COMMON)) JOE_OBJS $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) vanGogh.o joe vanGogh 

