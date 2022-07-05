#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h> 
#include <time.h> 
#include <unistd.h> 
#include <getopt.h>

#include "tsunameter.h"
#include "maths_util.h"
#include "base_station.h"
#include "simulation.h"

/*
    This file contains the launching point of the simulation.

    See the documentation above the main function to see how to see what 
    arguments you must/can supply when running the compiled program

    See simulation.h to see all the simulation settings
    (e.g. max and min random water heights)
*/

// See definition for documentation
bool parse_args(int argc, char* argv[], int* m, int* n, float* threshold, int* fault_node, int* iterations);

// See definition for documentation
void print_arg_error_info();


/*
    The main function.

    See the following for information on what arguments to supply:

    Please use format: ./ass2.o <m> <n> <threshold>

    Where:  
            m = integer representing number of sensor network rows.
            n = integer representing number of sensor network columns.
            threshold = float respresenting the water column height threshold for the nodes. Heights are randomly generated between 6000 and 7500, so a good value is 6800.


    Optional Args:
            -i <iterations> = forces simulation run the given integer number of iterations (>0), instead of waiting on a sentinel values (REQUIRED FOR MONARCH). Each iteration is atleast 4 seconds long.
            -f <node> = introduces a fault into the given integer node (from 0 to m*n-1) randomly, 3-6 seconds after starting.


    NOTE:
            You need to supply m*n+1 processes to the program when running it

*/
int main(int argc, char* argv[])
{
    // Initialise MPI with full thread support which is required as the simulation uses many threads
    int given_thread_support;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &given_thread_support); 
    if(given_thread_support != MPI_THREAD_MULTIPLE)
    {
        printf("ERROR: Required level of threading is not supported!! \n");
        MPI_Finalize();
        return -1;
    }

    // Grab current process world rank and world communicator size
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); 

    // Parse input arguments to control the simulation
    float threshold = 0;
    int m, n, fault_node = -1, iterations = -1;
    bool success = parse_args(argc, argv, &m, &n, &threshold, &fault_node, &iterations);
    
    // If parsing did not succeed, then print the argument error message and exit
    if(!success && world_rank == STA_WORLD_RANK)
    {
        print_arg_error_info();
        exit(-1);
    }

    // Wait on a barrier, so that no thread proceeds if we have an error in the args
    MPI_Barrier(MPI_COMM_WORLD);

    // Excluding the base station process (world rank 0) and any extra unnecessarily supplied
    // processes; forms a 2D cartesian communicator of size m*n which represents the WSN.
    // For a non-tsunameter process, this returns a communicator consisting of all
    // non-tsunameter processes. For a tsunameter node process, this returns the cartesian
    // communicator that mimics the m x n organisation of tsunameter nodes.
    MPI_Comm split_comm;
    int min_tsu_world_rank, max_tsu_world_rank; // Store the resulting minimum and max world ranks of the tsunameter sensor nodes
    if(tsu_distribute_sensors(m, n, &split_comm, &min_tsu_world_rank, &max_tsu_world_rank) == -1)
    {
        printf("ERROR: Not enough processes to simulate all nodes, need %3.d got %3.d!\n", m * n + 1, world_size);
        MPI_Finalize();
        return -2;
    }

    // Split the execution path to either simulate a tsunameter or the base station
    // Processes whose world rank isn't the base station one (rank 0) or one of the 
    // tsunameter world ranks (between the min_tsu_world_rank and max_tsu_world_rank)
    // will simply pass this if statement and wait for termination of all other processes
    // which are actually being used. 
    if(between(world_rank, min_tsu_world_rank, max_tsu_world_rank))
    {
        tsu_simulate(split_comm, m, n, threshold, fault_node);
    } 
    else if(world_rank == STA_WORLD_RANK)
    {
        station_simulate(m, n, min_tsu_world_rank, max_tsu_world_rank, iterations);
    }

    // Flush the output and wait on a barrier to try make sure that this printf is always the last thing printed
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank == STA_WORLD_RANK)
        printf("Check the generated report.log file for more info.\n");

    // Ensure we free the created communicator, 
    // regardless of if its the cartesian one or non-cartesian one
    MPI_Comm_free(&split_comm);
   
    MPI_Finalize(); 
    return 0;
}

// -----------------------------------------------------------------------------

/*
    bool parse_args(int argc, char* argv[], int* m, int* n, float* threshold, int* fault_node, int* iterations);

    parses the main function arguments to try and determine output arguments m, n and the threshold; with 
    optional arguments fault_node and iterations.  

    Args:
        argc: the count of arguments. This should be passed directly from the main function argument
        argv: the arguments, this should be passed directly form the main function argument
        m: the parsed output number of rows to use in sensor node grid
        n: the parsed output number of columns to usein the sensor node grid
        threshold: the parsed output water column height threshold used by the 
            tsunameter nodes to determine when to begin a tsunami alert 
        fault_node: the optional parsed output representing the tsunameter sensor node to introduce
            a fault in randomly after 3-6 seconds.
        iterations: the optional parsed output representing the number of iterations 
            that the base station should run for

    Returns:
        True if the supplied arguments were valid, otherwise false if there was an error or invalid format
*/
bool parse_args(int argc, char* argv[], int* m, int* n, float* threshold, int* fault_node, int* iterations)
{
    bool error = false;

    // If we don't have at least 4 arguments, then the args are invalid so 
    // we skip trying to read them and say that there's an error. 
    if(argc >= 4)
    {
        // We must always be supplied m, n and the threshold, so read them in directly
        *m = atoi(argv[1]);
        *n = atoi(argv[2]);
        *threshold = atof(argv[3]);

        // Parse optional arguments such as the fault node and the iteration
        int option = 0;
        while((option = getopt(argc, argv, "f:i:")) != -1)
        {
            switch(option)
            {
                // The i option is for iterations, so parse that
                case 'i':
                {
                    *iterations = atoi(optarg);
                    if(*iterations <= 0)
                        error = true;
                    break;
                }
                // The f option is for the fault node, so parse that
                case 'f':
                {
                    *fault_node = atoi(optarg);
                    if(*fault_node < 0 || *fault_node >= (*m) * (*n))
                        error = true;
                    break;
                }
                // If we have an unknown option, then we have an error
                default:
                {
                    error = true;
                    break;
                }
            }
        }
    } else error = true;
        
    return !error;
}

// -----------------------------------------------------------------------------

/*
    void print_arg_error_info()

    Prints an informational message depicting the argument format 
    that needs to be followed when passed in when starting the program
   
*/
void print_arg_error_info()
{
    printf("INCORRECT ARGUMENT VALUE OR FORMAT!\n\n");
    printf("Please use format: ./ass2.o <m> <n> <threshold>\n\n");
    printf("Where:\n");
    printf("\t m = integer representing number of sensor network rows.\n");
    printf("\t n = integer representing number of sensor network columns.\n");
    printf("\t threshold = float respresenting the water column height threshold for the nodes. Heights are randomly generated between 6000 and 7500, so a good value is 6800.\n");
    printf("\n\n");
    printf("Optional Args:\n");
    printf("\t-i <iterations> = forces simulation run the given integer number of iterations (>0), instead of waiting on a sentinel values (REQUIRED FOR MONARCH). Each iteration is atleast %d seconds long.\n", STA_ITER_TIME_S);
    printf("\t-f <node> = introduces a fault into the given integer node (from 0 to m*n-1) randomly, 3-6 seconds after starting.\n");
    printf("\n\n");
    printf("NOTE:\n");
    printf("\tYou need to supply m*n+1 processes to the program when running it\n");

}