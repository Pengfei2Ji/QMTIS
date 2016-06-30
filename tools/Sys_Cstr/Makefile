#####################################################
# QMTIS QMTIMS QMTIS QMTIS QMTIS QMTIMS QMTIS QMTIS #
#####################################################
# Version 1.2
# Datae_ 2015-04-08
#####################################################
#
# Tool: Sys_Cstr (System Construction) ##############
#
# Developed by: Pengfei Ji 
#####################################################

# Compilers
CFLAGS = ${FLAGS}
CC     = gcc
#CC     = icc
FLAGS  = -O

ifneq (,$(strip $(findstring alpha,${HOSTTYPE})))
CC      = cc
FLAGS  += -DALPHA -O3 -float -fp_reorder -arch ev56 -tune host -w
endif

ifneq (,$(strip $(findstring sparc,${HOSTTYPE})))
CC      = cc
FLAGS  += -fast -xtarget=ultra3 -xarch=v9b
endif

# where to put executable
#BINDIR = ${HOME}/bin/${HOSTTYPE}
BINDIR = ../../bin/${HOSTTYPE}

############################################################################

# sample constructor
Sys_Cstr: Sys_Cstr.c
#	${CC} ${CFLAGS} -o ${BINDIR}/$@ Sys_Cstr.c -lm
	${CC} ${CFLAGS} -o ${BINDIR}/Sys_Cstr.X Sys_Cstr.c -lm


###########################################################################
#
#	 Rules
#
###########################################################################

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

# all objects depend on headers
${OBJECTS}: ${HEADERS}

.SUFFIXES:

# generic compilation rule
%.o : %.c
	${CC} ${CFLAGS} -c $<

# How to link
${MAKETARGET}: ${OBJECTS} 
	${CC} ${FLAGS} -o $@ ${OBJECTS} ${LIBS}
	mv $@ ${BINDIR}; rm -f $@

# First recursion only set the MAKETARGET Variable
.DEFAULT:
	${MAKE} MAKETARGET='$@' STAGE2

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${MAKETARGET})
	${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

# clean
clean:
	rm -f *.o *~ 

count:
	wc -l Makefile *.[ch] *.csh *.sh *_uv *.perl rav
