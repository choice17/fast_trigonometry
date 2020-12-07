CC=gcc
STRIP=strip

Q=@
CFLAGS=-O3 -g -fPIC
TARGET=fast_tri.so

SRCS=fast_tri.c
OBJS=fast_tri.o
INCS_DIR=.
INCS=fast_tri.h

PREFIX=install
DEST=${PREFIX}

CFLAGS+= -I$(INCS_DIR)

all: $(TARGET) install  

strip: $(TARGET)
	$(Q)strip $<

install: $(TARGET)
	install -d ${PREFIX}
	install -d ${PREFIX}/lib
	install -d ${PREFIX}/inc
	install -m 755 $(INCS) ${DEST}/inc
	install -m 755 $(TARGET) ${DEST}/lib

$(TARGET): $(OBJS)
	@printf "    CC $@\n"
	$(Q)$(CC) $(CFLAGS) $^ -o $@

%.o : %.c
	@printf "    CC $@\n"
	$(Q)$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(Q)rm -rf $(OBJS) $(DEST) $(TARGET)

