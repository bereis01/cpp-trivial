#include "main.h"

int wrong_vector_access(const uint8_t *vector, size_t size)
{
    int vector_size = vector[0];
    int sum = size;
    for (size_t i = 0; i < vector_size; i++)
    {
        sum += vector[i];
    }

    return sum;
}