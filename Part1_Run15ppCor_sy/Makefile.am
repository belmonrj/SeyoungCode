AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libFvtxHigh.la 

noinst_HEADERS = \
  pp_Cor.h \
  pp_filter.h \
  pp_filter_ep.h \
  FvtxHigh_LinkDef.h

libFvtxHigh_la_SOURCES = \
  pp_Cor.C \
  pp_filter.C \
  pp_filter_ep.C \
  FvtxHigh_Dict.C

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib -lTOAD 

libFvtxHigh_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lFROG \
  -lfun4allfuncs \
  -lreactionplane \
  -lDoubleInteractionUtil \
  -lfvtxgeom

testexternals_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libFvtxHigh.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
FvtxHigh_Dict.C: \
  pp_Cor.h \
  pp_filter.h \
  pp_filter_ep.h \
  FvtxHigh_LinkDef.h 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
