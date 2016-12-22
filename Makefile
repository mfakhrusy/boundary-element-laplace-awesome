#   
#	File name		: Makefile
#	Date			: 
#	Version			: 
#	Author			: 
#

DEST			= .

HDRS			= 	global.hpp				\
				math_libs/math_function.hpp		\
				pre_process/initialization.hpp		\
				pre_process/airfoil.hpp			\
				solver/matrix_init.hpp			\
				solver/matrix_solver.hpp		\
				misc/miscellaneous.hpp			\
				misc/post_computation.hpp		\
	
LIBS			=	

INPS			=	

COMPILER		= g++ 

OPTFLAG			= -std=c++11 -O2

MAKEFILE		= Makefile


PROGRAM			= BEM_laplace

SRCS			= main.cpp					\
			  math_libs/math_function.cpp			\
			  pre_process/initialization.cpp		\
			  pre_process/airfoil.cpp			\
			  solver/matrix_init.cpp			\
			  solver/matrix_solver.cpp			\
			  misc/miscellaneous.cpp			\
			  misc/post_computation.cpp			\

OBJS			= $(SRCS:.cpp=.o) 	

.cpp.o:
			$(COMPILER) $(OPTFLAG) -c $*.cpp -o $*.o 

all:			$(PROGRAM)

$(PROGRAM):		$(OBJS) $(LIBS)
			@echo -n "Loading Program $(PROGRAM) ... "
			@$(COMPILER) $(OPTFLAG) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
			@echo "done"

clean:;			@rm -f $(SRCS:.cpp=.o) $(SRCS:.cpp=.il) $(PROGRAM)


