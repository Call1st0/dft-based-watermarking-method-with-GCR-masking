UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
libGcr.so: Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o
	$(CXX) -shared $^ -o $@ -Wl,--whole-archive -llcms2 -Wl,--no-whole-archive -Wl,--whole-archive -llapack -Wl,--no-whole-archive -Wl,--whole-archive -lblas -Wl,--no-whole-archive -Wl,--whole-archive -larmadillo -Wl,--no-whole-archive

Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o : CXXFLAGS+=-fPIC
endif
ifeq ($(UNAME), Darwin)
libGcr.dylib: Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o
	$(CXX) -shared $^ -o $@ -Wl,-force_load /usr/local/lib/liblcms2.a  -Wl,-noall_load -Wl,-force_load /usr/local/opt/lapack/lib/libblas.dylib -Wl,-noall_load -Wl,-force_load /usr/local/opt/lapack/lib/liblapack.dylib -Wl,-noall_load -Wl,-force_load /usr/local/lib/libarmadillo.dylib -Wl,-noall_load

Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o : CXXFLAGS+=-fPIC
endif
