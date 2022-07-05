#include <time.h>
#include <pthread.h>

/*
    This file provides all public functions and data required
    to interface with the Satellite Altimeter. See the .c
    file for all private functions and data. 
*/

// guard
#ifndef SATELLITE_INCLUDE
#define SATELLITE_INCLUDE

#include <stdbool.h>

// ---------------------------------------------------------------------------------------------
//    PUBLIC CONSTANTS
// ---------------------------------------------------------------------------------------------

#define SAT_ARRAY_LENGTH 500     // length of the data array in satellite

// ---------------------------------------------------------------------------------------------
//    PUBLIC STRUCTS
// ---------------------------------------------------------------------------------------------

// The SatelliteData struct contains the data about each coord generated and stored by the satellite
typedef struct 
{
    float height;           // height of the water
    time_t timestamp;       // time the data was obtained
    int coordinates[2];     // coordinates of the data obtained

} SatelliteData;

// The Satellite struct holds all data for a Satellite 
typedef struct 
{
    // array which holds the Satellite's data
    SatelliteData data_array[SAT_ARRAY_LENGTH];

    // A flag which if set to 0 allows the simulate_satellite function  
    // to gather data iteratively otherwise causes satellite stops gathering data
    int terminate_satellite;

    // mutex used to protect the satellite data_array, write_index and data_count
    pthread_mutex_t satellite_mutex;

    // the index to write to data_array 
    int write_index;  

    // the amount of data stored in data_array
    int data_count;

    // rows, m, and columns, n, in tsunamter cartesian topology
    int m, n;
} Satellite;

// global satellite
// this variable is not valid until satellite_simulate is called in thread creation  
Satellite g_satellite;

// ---------------------------------------------------------------------------------------------
// PUBLIC FUNCTIONS
// ---------------------------------------------------------------------------------------------

/*
    void terminate_satellite() 

    Terminates the global satellite and stops any instances of 'satellite_simulate' from running. 

*/
void terminate_satellite();

/*
    void* satellite_simulate(void* grid_size_ptr)

    Initialises the Satellite struct and begins simulating it; producing data
    and storing it as required until the terminate_satellite function is called. 

    Args:
        grid_size_ptr: an integer array with two values, the first representing rows in the tsunameter
            cartesian topology and the second columns in the topology 

*/
void* satellite_simulate(void* grid_size_ptr);

// end guard
#endif