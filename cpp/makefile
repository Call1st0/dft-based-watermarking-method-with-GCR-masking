libGcr.so: Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o
	$(CXX) -shared $^ -o $@ -Wl,--whole-archive -llcms2 -Wl,--no-whole-archive -Wl,--whole-archive -llapack -Wl,--no-whole-archive -Wl,--whole-archive -lblas -Wl,--no-whole-archive -Wl,--whole-archive -larmadillo -Wl,--no-whole-archive

Gamut.o GCR_py_wrap.o lcmswrap.o RevInt.o : CXXFLAGS+=-fPIC
