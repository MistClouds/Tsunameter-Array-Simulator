#include <time.h> 
#include <mpi.h> 

/*
    This file provides all public funcitons and data required to 
    simulate and interface with a tsunameter sensor. See the .c
    file for all private functions and data.  
*/

// guard
#ifndef TSUNAMETER_INCLUDE
#define TSUNAMETER_INCLUDE

// ---------------------------------------------------------------------------------------------
//  STRUCTS
// ---------------------------------------------------------------------------------------------

// The TsuReport struct represents the report message sent from the Tsunameter node to 
// the base station whenever enough neighbour nodes agree with a tsunami alert. 
typedef struct {  
    
    // The node number which sent the report (AKA rank within cartesian communicator)
    int reporting_node;  

    // The sma height which triggered the report in the reporting tsunameter
    float trigger_height; 
    
    // The time that the triggering height was sensed
    time_t trigger_timestamp; 
    
    // The sma heights of the tsunameter's neighbour nodes. 
    // This is only filled up to neighbour_count in the TsuInfo 
    // struct for the tsunameter as declared below. 
    float neighbour_heights[4]; 

    // The number of neighbours which had similar heights, hence agreed with the alert
    int agreeing_neighbour_count; 

    // The time that the alert started being processing
    struct timespec alert_timestamp; 

    // The timestamp for the neighbour query being done. 
    // That is, this time is measured right after the node receives
    // responses about all neighbours heights when handling a tsunami alert. 
    struct timespec query_timestamp; 
} TsuReport;

// The TsuIno struct contains node related info about a tsunameter
typedef struct {
    // The node number (AKA rank within cartesian communicator)
    int node;

    // The number of neighbours the tsunameter has.
    // All the provided neighbour arrays are of size 4, but are 
    // only filled from index 0 to neighbour_count. 
    int neighbour_count;

    // The tsunameter's neighbour nodes (AKA their ranks within the cartesian communicator)
    int neighbour_nodes[4];

    // The coordinates of the node within the cartesian communicator
    int coord[2];

    // The ipv4 address of the host device excuting the tsunameter code
    char ipv4[16];

    // The name of the processor of the host device executing the tsunameter code
    char processor[MPI_MAX_PROCESSOR_NAME];
} TsuInfo;

// ---------------------------------------------------------------------------------------------
// PUBLIC FUNCTIONS
// ---------------------------------------------------------------------------------------------

/*
    MPI_Datatype tsu_report_datatype()

    Creates the MPI_Datatype required to send and receive the TsuReport struct
    directly in a message. The user is responsible for freeing this datatype
    after creation.

    Returns:
        The MPI_Datatype which represents the TsuReport struct. 

*/
MPI_Datatype tsu_report_datatype();

/*
    MPI_Datatype tsu_info_datatype()

    Creates the MPI_Datatype required to send and receive the TsuInfo struct
    directly in a message. The user is responsible for freeing this datatype
    after creation.

    Returns:
        The MPI_Datatype which represents the TsuInfo struct. 

*/
MPI_Datatype tsu_info_datatype();

/*
    int tsu_distribute_sensors(const int m, const int n, MPI_Comm* tsu_comm, int* min_tsu_world_rank, int* max_tsu_world_rank)

    Creates a cartesian communicator of size m * n. The communicator skips the process of world rank 0, hence starts
    using consecutive processes of rank 1 and higher up until it has all required m*n processes. The output argument 
    tsu_comm will contain the cartesian communicator iff the calling process ends up being assigned to the communicator. 
    Otherwise, it simply contains a communicator which has all excluded processes. It is up to the user to free the 
    communicator once they are done with it. 

    Args:
        m: the number of rows in the cartesian communicator
        n: the number of columns in the cartesian communicator
        tsu_comm: the output arg containing the cartesian communicator. 
            This is only a valid output if the calling process ends up being used in the communicator. 
        min_tsu_world_rank: the output arg containing the minimum tsunameter world rank. 
                That is, the starting world rank of the continuous chunk of ranks used for the cartesian communicator.    
        max_tsu_world_rank: the output arg containing the maximum tsunameter world rank. 
                That is, the ending (inclusive) world rank of the continuous chunk of ranks used for the cartesian communicator. 

    Returns:
        -1 on error, otherwise 0. It should only error if not enough processes are supplied to create the communicator. 
*/
int tsu_distribute_sensors(const int m, const int n, MPI_Comm* tsu_comm, int* min_tsu_world_rank, int* max_tsu_world_rank);

/*
    void tsu_simulate(const MPI_Comm tsu_cart_comm, const int m, const int n, const float threshold, const int fault_node)

    Starts simulating a tsunameter node, which will operate until a termination signal is sent to it from the base station.  
    The calling process must be in the supplied cartesian communicator, resulting from a call to tsu_distribute_sensors.

    Args:
        tsu_cart_comm: the cartesian communicator. 
            Note that the calling process must be in this communicator.  
        m: the number of rows in the cartesian communicator
        n: the number of columns in the cartesian communicator
        threshold: the water column sma height threshold that will have the tsunameter trigger start a tsunami alert 
        fault_node: the rank of the node within the cartesian communicator which will experience a simulated fault
            within a random time of 3-6 seconds after calling tsu_simulate. Use -1 to have no faulty nodes. 
*/
void tsu_simulate(const MPI_Comm tsu_cart_comm, const int m, const int n, const float threshold, const int fault_node);

// end guard
#endif