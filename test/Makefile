FC=mpiifort
MPAL_INCLUDE=../build/
MPAL_LIB=../build/

LIB_DIR=../lib/

# test1: test1.f90
# 	$(FC) -I$(MPAL_INCLUDE) -L$(MPAL_INCLUDE) -L$(LIB_DIR) $^ -o $@ -lmpal -lcadnaF -lrpe

test3: test3.f90
	$(FC) -I$(MPAL_INCLUDE) -L$(MPAL_INCLUDE) -L$(LIB_DIR) $^ -o $@ -lmpal -lcadnaF -lrpe -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

.PHONY: clean

clean:
	rm -f test1 test2
