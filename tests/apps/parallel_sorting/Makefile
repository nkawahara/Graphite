TARGET = p_sort
SOURCES = p_sort.cc

# Network Node ful Cores size ( is not worker Core 
CORES ?= 14
#CORES ?= 258
#CORES ?= 130

APP_FLAGS ?= $(CORES)
ENABLE_SM ?= true
CLOCK_SKEW_MANAGEMENT_SCHEME ?= lax_p2p

include ../../Makefile.tests
