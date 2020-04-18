# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch
-include ../MyConfig.mk

#------------------------------------------------------------------------------


ANAHV5O      = normalizeVar.o
ANAHV5S      = normalizeVar.C
ANAHV5       = normalizeVar.exe



OBJS          = $(ANAHV5O) 

PROGRAMS      = $(ANAHV5) 


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

		@echo "$@ done"


$(ANAHV5):      $(ANAHV5O)
		$(LD) $(LDFLAGS) $^ $(LIBS) -lMLP $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"




clean:
		@rm -f $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations \
		   files/*
		@rm -rf cxx_repository
		-@cd RootShower && $(MAKE) distclean
		-@cd rhtml && $(MAKE) distclean
		-@cd RootIDE && $(MAKE) distclean
		-@cd periodic && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

