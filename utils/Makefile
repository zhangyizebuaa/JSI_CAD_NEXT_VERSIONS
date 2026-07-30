#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = .
endif

$(TDIR)/main: $(ODIR)/main.o  $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/main $(ODIR)/main.o -L$(LDIR) $(RPATH) -legads -lm

$(ODIR)/main.o:	main.cpp $(IDIR)/egads.h $(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CXX) -c $(COPTS) $(DEFINE) -I$(IDIR) main.cpp -I/home/sx/OpenCASCADE-7.6.0/include/opencascade -I/home/sx/EngSketchPad/src/EGADS/include -I/home/sx/EngSketchPad/src/EGADS/util -I/home/sx/EngSketchPad/src/EGADS/src -o $(ODIR)/main.o

clean:
	-rm $(ODIR)/main.o

cleanall:	clean
	-rm $(TDIR)/main
