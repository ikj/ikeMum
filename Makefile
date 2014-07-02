# Build makefile for joe
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(MUM_HOME)/Makefile.in
include $(MUM_HOME)/src/common/Makefile.dep
include $(MUM_HOME)/src/ike/Makefile.dep

ADOLC_INC_DIR = $(ADOLC_HOME)/include/adolc
#ADOLC_LIB = -L$(ADOLC_HOME)/lib64 -ladolc -L$(COLPACK_HOME)/lib/ -lColPack -Wl,--rpath -Wl,$(ADOLC_HOME)/lib64
ADOLC_LIB = -L$(ADOLC_HOME)/lib64 -ladolc  $(COLPACK_HOME)/lib/libColPack.a -lm

COMMON = mpi_tommie.o MpiStuff.o Gp.o MiscUtils.o Ugp.o Param.o UgpWithTools.o CdpFilter.o MshFilter.o UgpWithCv2.o Logging.o
JOE = UgpWithCvCompFlow.o Scalars.o JoeWithModels.o 
ADJOINT = UgpWithCvCompFlowAD.o ScalarsAD.o JoeWithModelsAD.o  
IKE = IkeWithPsALC.o IkeWithModels.o IkeUtilsAD.o # ike.o is not included on purpose since user can have own joe.cpp
adjoint_H = UgpWithCvCompFlowAD.h JoeWithModelsAD.h 
ike_H = IkeUgpWithCvCompFlow.h IkeWithPsALC.h IkeWithModels.h IkeUtilsAD.h MatComprsed.h PetscSolver2.h ADvar.h IkeTurbModel_KOM.h 
OBJS = $(addprefix $(MUM_HOME)/src/common/, $(COMMON))\
       $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE))\
			 $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(ADJOINT))\
			 $(addprefix $(MUM_HOME)/src/ike/, $(IKE))\
       $(MUM_HOME)/src/common/mixing_combustion/CombustionGeneralDefinitions.o\
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/ludcmp.o\
       $(MUM_HOME)/src/common/mixing_combustion/Numerical_Recipes/qrdcmp.o\
       $(MUM_HOME)/src/common/mixing_combustion/ChemtableCartesianLinear_single.o\
       $(MUM_HOME)/src/common/mixing_combustion/Table4D.o\
       $(MUM_HOME)/src/common/mixing_combustion/PressureTable.o\
			 ike.o

.PHONY: default help clean

default: ike

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'ike :      compiles ike (Pseudo-arclength solver)'
	@echo 'clean :    clean joe (included common files!)'
	@echo ' '

JOE_LIBS = $(CLIBS) -lm
#IKE_LIBS = -I${ADOLC_INC} ${ADOLC_LIBOPTS} 

ike: COMMON_OBJS JOE_OBJS ADJOINT_OBJS IKE_OBJS ike.o 
	$(CC) -o $@ $(CLD) $(OBJS) $(JOE_LIBS) $(ADOLC_LIB)

# ike.o: ike.o treated separately because can be changed by user in case directory
ike.o: ike.cpp ikeHyShot.h $(joe_H) $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(adjoint_H)) $(addprefix $(MUM_HOME)/src/ike/, $(ike_H)) $(MUM_HOME)/src/ike/IkeTurbModel_KOM.h $(addprefix $(MUM_HOME)/src/ike/, JoeDebug.h) $(addprefix $(MUM_HOME)/src/ike/JOE/combModels/, CombModel_TripletMixing_ik.h) $(addprefix $(MUM_HOME)/src/ike/JOE/combModels/, CombModel_QuadrupletMixing_ik.h)
	$(CC) $(CFLAGS) -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/common/mixing_combustion -I$(MUM_HOME)/src/ike/JOE -I$(MUM_HOME)/src/ike -c ike.cpp -o ike.o

# COMMON_OBJS: For common, you can simply make it
COMMON_OBJS:  
	$(MAKE) $(COMMON) common_comb -C $(MUM_HOME)/src/common 

# JOE_OBJS: For joe, each object file is compiled seperately since you cannot simply make joe because of ike
JOE_OBJS:  $(addprefix $(MUM_HOME)/src/ike/JOE/, $(JOE)) 

$(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o: $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.cpp -o $(MUM_HOME)/src/ike/JOE/UgpWithCvCompFlow.o

$(MUM_HOME)/src/ike/JOE/Scalars.o: $(MUM_HOME)/src/ike/JOE/Scalars.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/Scalars.cpp -o $(MUM_HOME)/src/ike/JOE/Scalars.o

$(MUM_HOME)/src/ike/JOE/JoeWithModels.o: $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp $(JoeWithModels_H)
	$(CC) $(CFLAGS) -I$(MUM_HOME)/src/common -c $(MUM_HOME)/src/ike/JOE/JoeWithModels.cpp -o $(MUM_HOME)/src/ike/JOE/JoeWithModels.o

# ADJOINT_OBJS: For adjoint, each object file is compiled seperately
ADJOINT_OBJS: $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(ADJOINT)) 

$(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/JoeWithModelsAD.o: $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/JoeWithModelsAD.cpp $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/JoeWithModelsAD.h $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h $(MUM_HOME)/src/ike/IkeUgpWithCvCompFlow.h
	$(CC) $(CFLAGS) -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/JoeWithModelsAD.cpp -o $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/JoeWithModelsAD.o 

$(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.o: $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.cpp $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h $(MUM_HOME)/src/ike/IkeUgpWithCvCompFlow.h $(MUM_HOME)/src/ike/ADvar.h 
	$(CC) $(CFLAGS) -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.cpp -o $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.o 

$(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/ScalarsAD.o: $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/ScalarsAD.cpp $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h $(MUM_HOME)/src/ike/IkeUgpWithCvCompFlow.h  $(MUM_HOME)/src/ike/ADvar.h 
	$(CC) $(CFLAGS) -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/ScalarsAD.cpp -o $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/ScalarsAD.o

# IKE_OBJS:
IKE_OBJS: $(addprefix $(MUM_HOME)/src/ike/, $(IKE)) 

$(MUM_HOME)/src/ike/IkeUtilsAD.o: $(MUM_HOME)/src/ike/IkeUtilsAD.cpp $(MUM_HOME)/src/ike/IkeUtilsAD.h $(MUM_HOME)/src/ike/MatComprsed.h 
	$(CC) $(CFLAGS)  -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/IkeUtilsAD.cpp -o $(MUM_HOME)/src/ike/IkeUtilsAD.o 

$(MUM_HOME)/src/ike/IkeWithModels.o: $(MUM_HOME)/src/ike/IkeWithModels.cpp $(MUM_HOME)/src/ike/IkeWithModels.h $(MUM_HOME)/src/ike/IkeUtilsAD.h $(MUM_HOME)/src/ike/MatComprsed.h $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(adjoint_H)) $(MUM_HOME)/src/ike/IkeUgpWithCvCompFlow.h $(MUM_HOME)/src/ike/IkeTurbModel_KOM.h $(MUM_HOME)/src/ike/ADvar.h  
	$(CC) $(CFLAGS)  -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/IkeWithModels.cpp -o $(MUM_HOME)/src/ike/IkeWithModels.o 

$(MUM_HOME)/src/ike/IkeWithPsALC.o: $(MUM_HOME)/src/ike/IkeWithPsALC.cpp $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(adjoint_H)) $(addprefix $(MUM_HOME)/src/ike/, $(ike_H))  $(MUM_HOME)/src/ike/ADvar.h
	$(CC) $(CFLAGS)  -I$(ADOLC_INC_DIR) -I$(MUM_HOME)/src/common -I$(MUM_HOME)/src/ike/JOE -c $(MUM_HOME)/src/ike/IkeWithPsALC.cpp -o $(MUM_HOME)/src/ike/IkeWithPsALC.o 

clean:
	rm $(OBJS) ike  

cleanike:
	rm $(addprefix $(MUM_HOME)/src/ike/JOE/ADJOINT_FILES/, $(ADJOINT)) $(addprefix $(MUM_HOME)/src/ike/, $(IKE)) ike.o ike


