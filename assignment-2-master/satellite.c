#include "satellite.h"

#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>

#include "time.h"
#include "maths_util.h"
#include "base_station.h"
#include "simulation.h"

// ---------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------------------------------

/*
    void initialize_satellite(int m, int n)

    Initialises the global satellite struct. The struct should not be used until this
    has been called, and this should not be called again until terminate is called first. 

    Args:
        m: represents the number of rows of tsunameters
        n: represents the number of columns of tsunameters
*/
void initialize_satellite(int m, int n) {
    // initialising variables
    g_satellite.data_count = 0;
    g_satellite.write_index = 0;
    g_satellite.terminate_satellite = 0;
    g_satellite.m = m;
    g_satellite.n = n;

    // seed the satellite PRNG
    srand(time(NULL));

    // create the pthread mutex
    pthread_mutex_init(&g_satellite.satellite_mutex, NULL);
}

// ---------------------------------------------------------------------------------------------

/*
    void terminate_satellite() 

    Terminates the global satellite and stops any instances of 'satellite_simulate' from running. 
*/
void terminate_satellite() { 
    // set the termination signal
    g_satellite.terminate_satellite = 1;

    // destroying satellite mutex
    pthread_mutex_destroy(&g_satellite.satellite_mutex);
}

// ---------------------------------------------------------------------------------------------

/*
    void* satellite_simulate(void* grid_size_ptr)

    Initialises the Satellite struct and begins simulating it; producing data
    and storing it as required until the terminate_satellite function is called. 

    Args:
        grid_size_ptr: an integer array with two values, the first representing rows in the tsunameter
            cartesian topology and the second columns in the topology 

*/
void* satellite_simulate(void* grid_size_ptr) {
    
    int* grid_size = (int*) grid_size_ptr;
    int m = grid_size[0]; // number of rows in tsunameter cartesian topology
    int n = grid_size[1]; // number of columns in tsunameter cartesian topology
    initialize_satellite(m, n); // initiate satellite variables

    // initialise two variables to be used for time
    struct timespec start_time, end_time;

    // beginning of data production, occurs while there is no termination message
    while (g_satellite.terminate_satellite == false) {
        clock_gettime(CLOCK_MONOTONIC, &start_time);  
        
        // lock the mutex  
        pthread_mutex_lock(&g_satellite.satellite_mutex);

        // produce data for each m and n coordinate
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
    	        SatelliteData satellite_data; // initialise variable which holds satellite data

                // save the coordinates in satellite_data
                satellite_data.coordinates[0] = i;
                satellite_data.coordinates[1] = j;

    	        // generate the height for coordinate and add to satellite_data
   	            satellite_data.height = random_float(SAT_MIN_WATER, SAT_MAX_WATER);
    
                // obtain current time and save in satelite_data
                satellite_data.timestamp = time(NULL);

                // add satellite_data to the satellite's data array write_index is used as
                // it holds the position in the array that is to be written to next
                g_satellite.data_array[g_satellite.write_index] = satellite_data;
                
                g_satellite.write_index = (g_satellite.write_index + 1) % SAT_ARRAY_LENGTH; // update the write index to the next element in the array
                
                // Update the amount of data stored in the data array by adding 1, but then capping
                // it to the size of the array. 
                g_satellite.data_count = min(g_satellite.data_count + 1, SAT_ARRAY_LENGTH);

            }
        }

        // unlock the mutex as finished fully updating satellite data
        pthread_mutex_unlock(&g_satellite.satellite_mutex);
        
        // Calculate amount of time we need to sleep until the next cycle
        clock_gettime(CLOCK_MONOTONIC, &end_time);  
        double sleep_time_s = (double)SAT_CYCLE_TIME_S - elapsed_s(start_time, end_time);

        // Only sleep if the sleep time is positive (we aren't late) and we aren't supposed to terminate
        if(sleep_time_s > 0 && !g_satellite.terminate_satellite)
            usleep(sleep_time_s * 1e6); 
    }

    return NULL;
}