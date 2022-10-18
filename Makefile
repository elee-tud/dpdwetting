
CC		= mpicxx
INC		= 
LIB		= 
CXXFLAGS	=  -g -std=c++0x -O2 

SRCDIR1 = src
SRCDIR2 = src/extensions src/externalforce src/integrators src/nonbonded src/bonded src/tools
SRCDIR3 = src/analysis src/analysis/droplet src/analysis/capillarybridge src/analysis/polymers src/analysis/general

OBJDIR = build
BINDIR = bin
LIBDIR = lib

VPATH= $(SRCDIR1) $(SRCDIR2) $(SRCDIR3)

SRCS= $(foreach dir, .. $(VPATH), $(wildcard $(dir)/*.cpp))
HEADS:= $(SRCS:.cpp=.hpp) ../src/types.hpp

SRCS := $(notdir $(SRCS))

OBJS := $(SRCS:%.cpp=$(OBJDIR)/%.o)




$(OBJDIR)/%.o: %.cpp %.hpp Makefile
	@echo "COMPILING SOURSE $< INTO OBJECT $@"
	$(CC) $(CXXFLAGS) -c $< -o $@


TARGET = $(BINDIR)/dpdwetting
TARGETLIB = $(LIBDIR)/libdpd.a


all		: $(TARGET) $(TARGETLIB)
	if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi

main	: $(TARGET)
lib		: $(TARGETLIB)


$(TARGET)	: $(OBJS) $(HEADS)
	if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CC) -o $(TARGET) $(OBJS) -I$(INC) -L$(LIB) $(CFLAGS) -lm 

$(TARGETLIB)	: $(OBJS) $(HEADS)
	if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	ar rcu $@ $+
	ranlib $@

clean	: 
	rm -rf $(OBJS) $(TARGET) $(TARGETLIB)
	@rmdir $(OBJDIR) $(BINDIR) $(LIBDIR)

depend	: $(SRCS)
	$(CC) -M $(CFLAGS) $^ > $@


