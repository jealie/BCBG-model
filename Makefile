CC      = g++
CFLAGS  = -O3
#CFLAGS  = -g
LDFLAGS = 

all: single_run clean

single_run: single_run.o bcbg2.o run_sim.o cells.o helper_fct.o multichannelsnucleus.o singlechannelnucleus.o
	@echo " [BUILD] $@"
	$(CC) -o $@ $^ $(LDFLAGS)

single_run.o: single_run.cpp helper_fct.hpp constants.hpp bcbg2.hpp run_sim.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

bcbg2.o: bcbg2.cpp constants.hpp bcbg2.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

run_sim.o: run_sim.cpp constants.hpp bcbg2.hpp run_sim.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

cells.o: cells.cpp constants.hpp bcbg2.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

helper_ftc.o: run_sim.o helper_ftc.cpp helper_fct.hpp constants.hpp bcbg2.hpp run_sim.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

singlechannelnucleus.o: singlechannelnucleus.cpp bcbg2.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

multichannelsnucleus.o: multichannelsnucleus.cpp bcbg2.hpp
	@echo " [BUILD] $@"
	$(CC) -c $(CFLAGS) $<

.PHONY: clean cleanest

clean:
		rm *.o

cleanest: clean
		rm bcbg2

