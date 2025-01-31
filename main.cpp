#include "main.h"

void wrong_vector_access(const uint8_t *vector, size_t size)
{
    int i = 0;
    while (true)
    {
        if (vector[i] == size || vector[i] == 0)
        {
            break;
        }

        i++;
    }
}