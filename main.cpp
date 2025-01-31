#include "main.h"

int wrong_vector_access(const uint8_t *vector, size_t size)
{
    int sum = 0;
    for (size_t i = 0; i < size + 1; i++)
    {
        sum += vector[i];
    }

    return sum;
}