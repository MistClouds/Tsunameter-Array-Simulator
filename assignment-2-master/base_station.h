

/*
    This file provides all public funcitons and data required to 
    simulate the base station. See the .c file for all private
    functions and data.  
*/

// guard
#ifndef BASE_STATION_INCLUDE
#define BASE_STATION_INCLUDE

// ---------------------------------------------------------------------------------------------
// PUBLIC FUNCTIONS
// ---------------------------------------------------------------------------------------------

/*
    void station_simulate(int min_tsu_world_rank, int max_tsu_world_rank, int m, int n, int req_iterations)

    Starts simulating the base station. If a >0 number of iterations is provided, then the
    station will perform the given number of iterations before shutting down. Otherwise it will wait
    for a enter to be pressed in the console before terminating. 

    Args:
        m: the number of rows in the cartesian communicator
        n: the number of columns in the cartesian communicator
        min_tsu_world_rank: The starting world rank of the continuous chunk of ranks used for the cartesian communicator of the tsunameter nodes.    
        max_tsu_world_rank: The (inclusive) ending world rank of the continuous chunk of ranks used for the cartesian communicator of the tsunameter nodes.    
        req_iterations: If >0, this represents the number of iterations to be performed by the Base Station before terminating

*/
void station_simulate(int m, int n, int min_tsu_world_rank, int max_tsu_world_rank, int req_iterations);

#endif