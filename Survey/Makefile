LDFLAGS += -Wl,--no-as-needed $(shell root-config --libs) -lMinuit -lGeom -ggdb3
CXXFLAGS += $(shell root-config --cflags) -ggdb3 -I.
CXXFLAGS += -Wextra

DEP_OBJ=Castor.h Point.h
EXE=main

all: $(EXE)

.PHONY: all clean

$(EXE): $(DEP_OBJ)

clean:
	rm -f *~ *.o *.so $(EXE) *.d
