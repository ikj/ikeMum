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

COMMON_JOE = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv2.o Logging.o
JOE_OBJS = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o
OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_OBJS)) joe.o

JOE_LIBS = $(CLIBS) -lm


joe: obj_joe joe.o 
	$(CC) -o $@ $(CLD) $(OBJS) $(JOE_LIBS) #; mv $@ $(MUM_HOME)/bin

obj_joe: 
	$(MAKE) obj -C $(MUM_HOME)/src/joe

joe.o: joe.cpp $(joe_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -c $<


##########
# VANGOGH (VANGOGH requires JOE)
##########
JOE_OBJS_FOR_VANGOGH = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o 

VANGOGH_OBJS = DiffFilter.o vanGoghWithModels.o

vanGogh: COMMON_OBJS_FOR_VANGOGH $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_OBJS_FOR_VANGOGH)) $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) vanGogh.o 
	$(CC) -o $@ $(CLD) $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE_OBJS_FOR_VANGOGH)) $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) vanGogh.o $(JOE_LIBS)  

$(MUM_HOME)/src/vangogh/DiffFilter.o: $(MUM_HOME)/src/vangogh/DiffFilter.cpp $(MUM_HOME)/src/vangogh/DiffFilter.h $(UgpWithCV2_H) $(Param_H) $(myMem_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/vangogh/DiffFilter.cpp -o $(MUM_HOME)/src/vangogh/DiffFilter.o

$(MUM_HOME)/src/vangogh/vanGoghWithModels.o: $(MUM_HOME)/src/vangogh/vanGoghWithModels.cpp $(MUM_HOME)/src/vangogh/vanGoghWithModels.h $(MUM_HOME)/src/vangogh/vanGoghUtils.h $(MUM_HOME)/src/vangogh/DiffFilter.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/vangogh/vanGoghWithModels.cpp -o $(MUM_HOME)/src/vangogh/vanGoghWithModels.o

vanGogh.o: vanGogh.cpp $(MUM_HOME)/src/vangogh/vanGoghWithModels.h $(MUM_HOME)/src/vangogh/vanGoghUtils.h $(MUM_HOME)/src/vangogh/DiffFilter.h $(joe_H) 
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -I$(MUM_HOME)/src/ike/JOE/turbModels -I$(MUM_HOME)/src/ike/JOE/combModels -I$(MUM_HOME)/src/vangogh -c $<

COMMON_OBJS_FOR_VANGOGH: $(addprefix $(MUM_HOME)/src/common/, $(COMMON_JOE)) 
	$(MAKE) $(COMMON) -C $(MUM_HOME)/src/common

## Note: For joe, each object file is compiled seperately since you cannot simply make joe because of ike
#$(MUM_HOME)/src/joe/UgpWithCvCompFlow.o: $(MUM_HOME)/src/joe/UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
#	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp -o $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o
#
#$(MUM_HOME)/src/joe/Scalars.o: $(MUM_HOME)/src/joe/Scalars.cpp $(UgpWithCvCompFlow_H)
#	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/Scalars.cpp -o $(MUM_HOME)/src/ike/JOE/Scalars.o
#
#$(MUM_HOME)/src/joe/JoeWithModels.o: $(MUM_HOME)/src/joe/JoeWithModels.cpp $(JoeWithModels_H)
#	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp -o $(MUM_HOME)/src/ike/JOE/JoeWithModels.o

##########
# ALL
##########
clean:
	rm $(OBJS) $(addprefix $(MUM_HOME)/src/vangogh/, $(VANGOGH_OBJS)) *.o joe vanGogh 

