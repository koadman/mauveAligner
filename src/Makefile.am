
if DEBUG
D_CXXFLAGS = -Wall -g -DCOMMAND_LINE -D__GNDEBUG__ -O2
LIBGENOME = Genomed
LIBMEMS = Memsd
else
P_CXXFLAGS = -O2 -DCOMMAND_LINE
LIBGENOME = Genome
LIBMEMS = Mems
endif
CXXFLAGS = $(D_CXXFLAGS) $(P_CXXFLAGS) `wx-config --cxxflags` 

bin_PROGRAMS = gril 
gril_SOURCES = \
memsApp.cpp

LIBRARY_CL = -l$(LIBGENOME)  -l$(LIBMEMS) -L@CONFIG_LIBGENOMEDIR@/lib  -L@CONFIG_LIBMEMSDIR@/lib 

gril_LDADD = $(LIBRARY_CL)
gril_LDFLAGS = `wx-config --libs`

INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_WXWINDOWSDIR@/include/ -I$(top_srcdir)/include/gn -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn -I@CONFIG_LIBMEMSDIR@/include

EXTRA_PROGRAMS = grilStatic

grilStatic_SOURCES = memsApp.cpp
grilStatic_LDADD = @CONFIG_LIBMEMSDIR@/lib/libMems.a @CONFIG_LIBGENOMEDIR@/lib/libGenome.a 
grilStatic_LDFLAGS = -static `wx-config --static --libs` 
