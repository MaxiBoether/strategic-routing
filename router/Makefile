ifeq ($(TYPE),RELEASE)
CPPFLAGS+=-Ofast
else
CPPFLAGS+=-g -fstack-protector-strong -fstack-clash-protection -fcf-protection -Wall -Wpedantic -Wextra -std=c++2a -D_GLIBCXX_ASSERTIONS 
ifeq ($(SANITIZE),THREAD)
CPPFLAGS+=-fsanitize=thread
else
CPPFLAGS+=-fsanitize=address
endif

ifeq ($(TYPE),DEBUG)
CPPFLAGS+=-O0 -Wshadow -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wrestrict -Wnull-dereference -Wno-missing-braces -Wno-missing-field-initializers -Winline -Wundef -Wcast-qual -Wredundant-decls -Wunreachable-code -Wstrict-overflow=5
else
CPPFLAGS+=-O2 -D_FORTIFY_SOURCE=2
endif

endif

CPPFLAGS+=-Wall -Wpedantic -Wextra -std=c++17 -fopenmp
SRCDIR=src
BUILDDIR=build
INCDIRS=include lib

LDLIBS+=-lm -ltinyxml2 -lgsl -lgslcblas
CORE_OBJS=$(BUILDDIR)/io.o $(BUILDDIR)/data.o $(BUILDDIR)/psychmod.o
INC=$(addprefix -I ,$(INCDIRS))

ifndef PSYCHMOD
PSYCHMOD=user_equilibrium_2r
endif
CPPFLAGS+=-DPSYCH_MODEL_CLASS=$(PSYCHMOD)

ifndef STRATEGY
STRATEGY=ssotd
endif

ROUTING_STRATEGY=$(BUILDDIR)/$(STRATEGY).o

ifeq ($(STRATEGY),ssotd)
ifndef SSOTD_VARIANT
SSOTD_VARIANT=fulldisjoint
endif
SSOTD_VARIANT_FILE=ssotd_$(SSOTD_VARIANT)
ROUTING_STRATEGY=$(BUILDDIR)/$(SSOTD_VARIANT_FILE).o
endif

ifeq ($(STRATEGY),ssotd)
	ADDITIONALS=$(addprefix $(BUILDDIR)/,ssotd_core.o)
else ifeq ($(STRATEGY),ea)
	ADDITIONALS=$(addprefix $(BUILDDIR)/,ea_io.o ea_islands.o ea_logging.o ea_mutations.o ea_util.o ea_scoring.o ea_crossover.o)
endif

define cc-command
$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
endef


all: $(BUILDDIR)/router
poa: $(BUILDDIR)/poa

$(BUILDDIR)/router: $(BUILDDIR)/router.o $(CORE_OBJS) $(ADDITIONALS) $(ROUTING_STRATEGY)
	$(CXX) $(CPPFLAGS) $(INC) $^ -o $@ $(LDLIBS)

$(BUILDDIR)/poa: $(BUILDDIR)/poa.o $(CORE_OBJS)
	$(CXX) $(CPPFLAGS) $(INC) $^ -o $@ $(LDLIBS)

$(BUILDDIR)/%.o: src/*/%.cpp
	$(cc-command)

clean:
	$(RM) -f $(BUILDDIR)/*.o $(BUILDDIR)/router $(BUILDDIR)/poa

.PHONY: all clean
