
FFT_CODE = fft_fftw
#		fft_fftw 
#		fft_nr 
#		fft_acml 
#		fft_mkl 

FFT_SRC = $(FFT_CODE)

FFT_OBJS = ./fft/$(FFT_CODE).o

fft_fftw : $(fft_fftw)
	$(CC) $(OPTFLAGS) $(TOMOMPI_INC) -I$(USR_LOCAL_SHARED)/include \
	-c ./fft/fft_fftw.c \
	-o ./fft/fft_fftw.o
	@echo ' '

fft_nr : $(fft_nr)
	$(CC) $(OPTFLAGS) $(TOMOMPI_INC) \
	-c ./fft/fft_nr.c \
	-o ./fft/fft_nr.o
	@echo ' '

fft_acml : $(fft_acml)
	$(CC) $(OPTFLAGS) $(TOMOMPI_INC) \
	-c ./fft/fft_acml.c \
	-o ./fft/fft_acml.o
	@echo ' '

fft_mkl : $(fft_mkl)
	$(CC) $(OPTFLAGS) $(TOMOMPI_INC) \
	-c ./fft/fft_mkl.c \
	-o ./fft/fft_mkl.o
	@echo ' '

fft_clean : $(fft_clean)
	/bin/rm -f ./fft/*.o ./fft/*~
	