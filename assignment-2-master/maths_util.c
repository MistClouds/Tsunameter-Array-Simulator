#include "maths_util.h"

// ---------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------------------------------

/*
    void seed_with_offset(int offset)

    Seeds the C pseudo-random number generator, using an offset 
    to assist in the seeding process. The offset is useful for 
    seeding different processes, where seeding based on time 
    alone may still result in the same seeds for each process;
    hence the same PRNG values being generated. 

    Inspired by: https://stackoverflow.com/questions/23742550/mpi-rand-gives-the-same-constant-numbers-across-all-processes-in-every-run

    Args:
        offset: the integer offset applied to the seed 

*/
void seed_with_offset(int offset)
{
    srand(offset + time(NULL)); 
}

// ---------------------------------------------------------------------------------------------

/*
    float random_float(float min, float max)

    Generates a psuedo-random float between the provided
    minimum and maximum float values. 

    Inspired By (provided in assignment spec): https://stackoverflow.com/questions/13408990/how-to-generate-random-float-number-in-c

    Args:
        min: the minimum generatable random number
        max: the maximum generatable random number

    Returns:
        The pseudo-randomly generated float
*/
float random_float(float min, float max)
{
    float range = max - min;
    return (float) (rand() / (float) (RAND_MAX/range)) + min;
}

// ---------------------------------------------------------------------------------------------

/*
    float average(float* data, const int count)

    Calculates the average of the given array

    Args:
        data: a pointer to a float array whose average we want
        count: the amount of data to include in the average (or length of array)

    Returns:
        The average
*/
float average(float* data, const int count)
{
    float sum = 0.0;
    for(int i = 0; i < count; i++)
    {
        sum += data[i];
    }
    return sum / count;
}

// ---------------------------------------------------------------------------------------------

/*
    bool between(int value, int min, int max)

    Returns true if the given value is between the
    given minimum and maximum values (inclusive).

    Args:
        value: the value to check the 'inbetweeness' of
        min: the minimum value
        max: the maximum value

    Returns:
        True if value is between min and max (inclusive), otherwise false. 
*/
bool between(int value, int min, int max)
{
    return value >= min && value <= max;
}

// ---------------------------------------------------------------------------------------------

/*
    int max(int v1, int v2)

    Returns the higher of the provided integers

    Args:
        v1: the first value
        v2: the second value

    Returns:
        The higher of v1 and v2
*/
int max(int v1, int v2)
{
    return v1 > v2 ? v1 : v2;
}

// ---------------------------------------------------------------------------------------------

/*
    int min(int v1, int v2)

    Returns the lower of the provided integers

    Args:
        v1: the first value
        v2: the second value

    Returns:
        The lower of v1 and v2
*/
int min(int v1, int v2)
{
    return v1 < v2 ? v1 : v2;
}

// ---------------------------------------------------------------------------------------------

/*
    double elapsed_s(struct timespec start, struct timespec end)

    Calculates the elapsed time between the provided start and end timespecs
    in seconds, but with nanosecond level precision. 

    Args:
        start: the start time of the 'elapsed' period we want to calculate
        end: the end time of the 'elapsed' period we want to calculate

    Returns:
        The elapsed time in seconds
*/
double elapsed_s(struct timespec start, struct timespec end)
{
    return ((end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec)) * 1e-9;
}

