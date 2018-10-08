	
RECON_SRC = recon_algorithm \
			filteredbackprojection \
			gridrec

RECON_OBJS = ./recon/recon_algorithm.o \
			./recon/filteredbackprojection.o \
			./recon/gridrec.o
			
recon_algorithm : $(recon_algorithm)
	$(CCC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./recon/recon_algorithm.cpp \
	-o ./recon/recon_algorithm.o 
	@echo ' '

filteredbackprojection : $(filteredbackprojection)
	$(CCC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./recon/filteredbackprojection.cpp \
	-o ./recon/filteredbackprojection.o 
	@echo ' '

gridrec : $(gridrec)
	$(CCC) $(UTILITY_INC) $(TOMOMPI_INC) $(OPTFLAGS) \
	-c ./recon/gridrec.cpp \
	-o ./recon/gridrec.o 
	@echo ' '

	
recon_clean : $(recon_clean)
	/bin/rm -f ./recon/*.o ./recon/*~ 
	


