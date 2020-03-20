CXXFLAGS := -fPIC -O3 -g
LDFLAGS = -shared -llcms2 -larmadillo
CXX := g++
NAME := wmgcr
PREFIX = /usr

lib = lib$(NAME).so
src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)

$(lib): $(obj)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(lib) 

.PHONY: install
install: $(solib)
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	cp $(lib) $(DESTDIR)$(PREFIX)/lib/$(lib)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/$(lib)

