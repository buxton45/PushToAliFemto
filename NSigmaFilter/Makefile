LIBS   = $(shell root-config --libs)
CFLAGS = $(shell root-config --cflags) -Wall -I${KLAM_ANALYZE_DIR} -g -O0

MY_LIBS = $(addprefix build/, AliFemtoNSigmaFilter.o)

all: runTest

runTest: ${MY_LIBS} test.C
	g++ test.C ${MY_LIBS} -o runTest $(LIBS) $(CFLAGS)

build/%.o: %.cxx %.h
	g++ -c $< $(CFLAGS) -o $@

clean:
	rm -f runTest ${MY_LIBS}
