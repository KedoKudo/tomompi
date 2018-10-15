	
MAIN_SRC = tomompi_sinoserver \
			tomompi_client \
			tomompi \
			tomosolo
	
MAIN_PROGS = link_tomompi link_tomosolo

INSTALL_PROGS = install

install: $(install)
	cp ./main/tomompi ./bin/tomompi
	cp ./main/tomosolo ./bin/tomosolo

link_tomompi: $(link_tomompi)
	@echo 'Linking tomompi...'
	$(CCLINKER) -static \
	$(OTHER_LIBS) \
	$(NEXUS_LIBS) $(FFTW_LIBS) \
	$(UTILITY_OBJS) $(NAPI_OBJS) \
	$(FFT_OBJS) $(SUPPORT_OBJS) $(RECON_OBJS) \
	$(HDF_LIBS) $(HDF5_LIBS) $(XML_LIBS) \
	./main/tomompi_sinogramserver.o \
	./main/tomompi_client.o \
	./main/tomompi.o \
	-o ./main/tomompi \
	$(OPTFLAGS)
	@echo '...link completed'
	@echo ' '
#	./compress.o \  not needed on tomo but not on blacklab?

link_tomosolo: $(link_tomosolo)
	@echo 'Linking tomosolo...'
	$(CCLINKER) $(HDF_LIBS) $(HDF5_LIBS) $(XML_LIBS) \
	$(OTHER_LIBS) \
	$(NEXUS_LIBS) $(FFTW_LIBS) \
	$(UTILITY_OBJS) $(NAPI_OBJS) \
	$(FFT_OBJS) $(SUPPORT_OBJS) $(RECON_OBJS) \
	./main/tomosolo.o \
	-o ./main/tomosolo \
	$(OPTFLAGS)
	@echo '...link completed'
	@echo ' '
		
tomompi: $(tomompi)
	$(CCC) $(TOMOMPI_INC) $(NAPI_INC) $(HDF_INC) $(NEXUSLIB_INC) \
	$(MPICH_INC) $(UTILITY_INC) $(OPTFLAGS) \
	-c ./main/tomompi.cpp \
	-o ./main/tomompi.o
	@echo ' '

tomompi_sinoserver: $(tomompi_sinoserver)
	$(CCC) $(TOMOMPI_INC) $(NAPI_INC) $(HDF_INC) $(NEXUSLIB_INC) \
	$(MPI_INC) $(UTILITY_INC) $(OPTFLAGS) \
	-c ./main/tomompi_sinogramserver.cpp \
	-o ./main/tomompi_sinogramserver.o
	@echo ' '

tomompi_client: $(tomompi_client)
	$(CCC) $(TOMOMPI_INC) $(NAPI_INC) $(HDF_INC) $(NEXUSLIB_INC) \
	$(MPI_INC) $(UTILITY_INC) $(OPTFLAGS) \
	-c ./main/tomompi_client.cpp \
	-o ./main/tomompi_client.o
	@echo ' '

tomosolo: $(tomosolo)
	$(CCC) $(TOMOMPI_INC) $(NAPI_INC) $(HDF_INC) $(NEXUSLIB_INC) \
	$(MPICH_INC) $(UTILITY_INC) $(OPTFLAGS) \
	-c ./main/tomosolo.cpp \
	-o ./main/tomosolo.o
	@echo ' '

	
main_clean : $(main_clean)
	/bin/rm -f ./main/*.o ./main/*~ 
	/bin/rm -f ./main/tomompi_sinogramserver
	/bin/rm -f ./main/tomompi_client
	/bin/rm -f ./main/tomompi
	/bin/rm -f ./bin/tomompi
	


