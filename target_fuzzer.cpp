#include <string>
#include <sstream>

#include "main.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    wrong_vector_access(data, size);

    return 0;
}