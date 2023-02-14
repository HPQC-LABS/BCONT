#  make FC=ifort
ifndef $FC
   FC = nagfor
   FC = f90             
   FC = gfortran      
   FC = ifort
endif
 
#option to choose level of optimization:  make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
#   FFLAGS =  -O2             # debug not specified, so optimize, nagfor: perform all checks except for -C=undefined
#   FFLAGS =  -O2             # debug not specified, so optimize, nagfor: check for undefined variables
  else                           
    FFLAGS =  -q  -C          # debug=1  (basic gtortran debug option)
#   FFLAGS =  -g     -u       # debug=1  (basic debug option)
    ifeq ($(debug),2)
        FFLAGS =  -C -g       # debug=2  (higher-level debug option)
    endif
#                 
endif
 
OBJECTS = bcont.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o bcont.x

# To run the code, execute:
# ./level.x < input.5 > fort.6

