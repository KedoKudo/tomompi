	
SUPPORT_SRC = filelistclass \
			centeringclass \
			errorlogclass

SUPPORT_OBJS = ./support/filelistclass.o \
			./support/centeringclass.o \
			./support/errorlogclass.o
		
filelistclass: $(filelistclass)
	$(CCC) $(NAPI_INC) $(HDF_INC) $(NEXUSLIB_INC) \
	$(MPI_INC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./support/filelistclass.cpp \
	-o ./support/filelistclass.o
	@echo ' '
	
centeringclass : $(centeringclass) 
	$(CCC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./support/centeringclass.cpp \
	-o ./support/centeringclass.o
	@echo ' '

errorlogclass : $(errorlogclass) 
	$(CCC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./support/errorlogclass.cpp \
	-o ./support/errorlogclass.o
	@echo ' '


support_clean : $(support_clean)
	/bin/rm -f ./support/*.o ./support/*~ 
	


