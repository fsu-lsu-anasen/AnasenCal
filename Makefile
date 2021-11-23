CC=g++
ROOTCFLAGS=`root-config --cflags`
ROOTGLIBS=`root-config --glibs`
CFLAGS=-std=c++11 -g -Wall $(ROOTCFLAGS)

INCLDIR=./include
SRCDIR=./src
OBJDIR=./objs

CPPFLAGS=-I $(INCLDIR)
LDFLAGS=$(ROOTGLIBS) -lSpectrum

SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

DICT_PAGES=$(INCLDIR)/DataStructs.h $(INCLDIR)/LinkDef_AnasenEvent.h
DICT=$(SRCDIR)/AnasenEvent_dict.cxx
LIB=$(OBJDIR)/AnasenEvent_dict.o
DICTSO=$(OBJDIR)/libAnasenEvent_dict.so

EXE=bin/anasencal

.PHONY: all clean

all: $(EXE)

$(EXE): $(LIB) $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(LIB): $(DICT)
	$(CC) $(CFLAGS) $(CPPFLAGS) -fPIC -I ./ -o $@ -c $^
	$(CC) $(CFLAGS) $(CPPFLAGS) -I ./ -fPIC -shared -o $(DICTSO) $@
	mv $(SRCDIR)/*.pcm ./bin/
	cp ./bin/*.pcm  $(OBJDIR)/

$(DICT): $(DICT_PAGES)
	rootcling -f $@ $^

clean:
	$(RM) $(OBJS) $(EXE) $(LIB) $(DICT) ./bin/*.pcm $(OBJDIR)/*.pcm $(DICTSO)

VPATH=$(SRCDIR)
$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $^

