CXXFLAGS += -std=c++11

all: target_fuzzer

clean:
	rm -fv *.a *.o *unittest *_fuzzer *_seed_corpus.zip crash-* *.zip

check: all
	./target_fuzzer test_data/*

# Fuzz targets.
target_fuzzer: target_fuzzer.cpp main.a
	${CXX} ${CXXFLAGS} $< main.a ${LIB_FUZZING_ENGINE} -o $@
	zip -q -r target_fuzzer_seed_corpus.zip test_data

# The library itself.
main.a: main.cpp main.h
	${CXX} ${CXXFLAGS} $< -c
	ar ruv main.a main.o