
.SUFFIXES:.C
.C.o:
		$(CC) -c -o $*.o $<

# Your LEDA 3.1.1 directory with the following corrections
# src/numbers/_real.c
#   line 609 (or so) immediately before "adjust_dbl()" insert:
#    "error_bf() = q;"
# 
#   line 728 add before the while loop in int real_rep::sign(const bigfloat&)
#          improve(current_precision);
#   line 745 replace "get_significant()" by "get_exponent()" 
LEDA = YOUR LEDA DIRECTORY

CC = CC -ptr/tmp -O -I$(LEDA)/incl

determinant: main.o ABDPY.o
	$(CC) -o determinant main.o ABDPY.o $(LEDA)/libL.a -lm
