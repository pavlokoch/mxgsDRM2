name := mxgsDRM
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
CPPFLAGS += -O2 -lgsl -lgslcblas -lsqlite3
