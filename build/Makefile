
CC		= mpicxx
INC		= 
LIB		= 
CXXFLAGS	=  -g -std=c++0x -O2 

SRCDIR1 = ../src
SRCDIR2 = ../src/extensions ../src/externalforce ../src/integrators ../src/nonbonded ../src/bonded ../src/tools
SRCDIR3 = ../src/analysis ../src/analysis/droplet ../src/analysis/capillarybridge ../src/analysis/polymers ../src/analysis/general

VPATH= $(SRCDIR1) $(SRCDIR2) $(SRCDIR3)

SRCS= $(foreach dir, .. $(VPATH), $(wildcard $(dir)/*.cpp))
HEADS:= $(SRCS:.cpp=.hpp) ../src/types.hpp

SRCS:= $(notdir $(SRCS))

OBJS= $(SRCS:.cpp=.o)



%.o: %.cpp %.hpp Makefile
	$(CC) $(CXXFLAGS) -c $< -o $@


TARGET = ../bin/dpdwetting
TARGETLIB = ../lib/libdpd.a


all		: $(TARGET) $(TARGETLIB)

main	: $(TARGET)
lib		: $(TARGETLIB)


$(TARGET)	: $(OBJS) $(HEADS)
	$(CC) -o $(TARGET) $(OBJS) -I$(INC) -L$(LIB) $(CFLAGS) -lm 
$(TARGETLIB)	: $(OBJS) $(HEADS)
	ar rcu $@ $+
	ranlib $@
	cp $(HEADS) ../include/

clean	: 
	rm -rf $(OBJS) $(TARGET) $(TARGETLIB)

depend	: $(SRCS)
	$(CC) -M $(CFLAGS) $^ > $@


