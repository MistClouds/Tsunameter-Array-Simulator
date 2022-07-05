#include "tsunameter.h"

#include <pthread.h>
#include <stdbool.h>
#include <stddef.h>
#include <unistd.h> 
#include <string.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <time.h>

#include "maths_util.h"
#include "message_type.h"
#include "simulation.h"

// ---------------------------------------------------------------------------------------------
//    PRIVATE CONSTANTS
// ---------------------------------------------------------------------------------------------

#define HIST_POINTS 5  // Amount of height history to store for the simple moving average 
#define QUERY_RESPONSE_POLL_TIME_MS 5 // The amount of time to sleep between polls when waiting for query responses from neighbours
#define REQ_AGREEING_NEIGHBOURS 2 // Amount of adjacent sensors which need to agree to confirm a tsunami alert and turn it into a report

// ---------------------------------------------------------------------------------------------
//    PRIVATE STRUCTS
// ---------------------------------------------------------------------------------------------

// The Tsunameter struct represents an instance of a tsunameter, hence
// holds all data and information pertaining to that tsunameter.  
typedef struct {
    // The cartesian communicator which the tsunameter is found in
    MPI_Comm cart_comm;
    
    // The rank of the tsunameter within its cartesian communicator
    int rank;

    // A flag representing whether the tsunameter should keep running or not.
    // If false, the tsunameter should shut down. 
    bool runflag;

    // The number of neighbours the tsunameter node has
    int neighbour_count;

    // The ranks the neighbours of the tsunameter, with respect to
    // the cartesian communicator. This array is only filled from
    // index 0 up to neighbour_count. 
    int neighbour_ranks[4];

    // The MPI datatype used to transmit the TsuReport struct through messages
    MPI_Datatype report_type;

    // A mutex used to protect all sensing data such as the last_sense_time, 
    // history, history_index, and height_sma from race conditions. 
    pthread_mutex_t data_mutex;
    
    // The timestamp of the last sensing operation done by the tsunameter
    time_t last_sense_time;

    // The rolling history of sensed water heights
    float history[HIST_POINTS];

    // The history index to write the next sensed 
    // height to in the history array
    int history_index;

    // The calculated simple moving average of the history array. 
    float height_sma;
} Tsunameter; 

// ---------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------------------------------

/*
    MPI_Datatype tsu_report_datatype()

    Creates the MPI_Datatype required to send and receive the TsuReport struct
    directly in a message. The user is responsible for freeing this datatype
    after creation.

    Returns:
        The MPI_Datatype which represents the TsuReport struct. 

*/
MPI_Datatype tsu_report_datatype() 
{
    // The struct we want to create a datatype for...
    /*
        typedef struct {  
            int reporting_node;  
            float trigger_height; 
            time_t trigger_timestamp; 
            float neighbour_heights[4]; 
            int agreeing_neighbour_count; 
            struct timespec alert_timestamp; 
            struct timespec query_timestamp; 
        } TsuReport;
    */

    // Number of members within the struct
    const int blocks = 7;

    // The length of each member in terms of its datatype. This is 1
    // unless the member is an array. For more complicated types such 
    // as time_t and timespec, we treat them as arrays of bytes.   
    const int block_lengths[7] = {
        1, 1, sizeof(time_t), 4, 1, sizeof(struct timespec), sizeof(struct timespec)  
    };

    // Memory offsets of all members within the struct
    const MPI_Aint block_displacements[7] = {
        offsetof(TsuReport, reporting_node),
        offsetof(TsuReport, trigger_height),
        offsetof(TsuReport, trigger_timestamp),
        offsetof(TsuReport, neighbour_heights),
        offsetof(TsuReport, agreeing_neighbour_count),
        offsetof(TsuReport, alert_timestamp),
        offsetof(TsuReport, query_timestamp),
    };

    // Types of each block/member of the struct.  
    // Again note that some members are sent as byte arrays instead. 
    const MPI_Datatype block_types[7] = {
        MPI_INT,
        MPI_FLOAT,
        MPI_BYTE,
        MPI_FLOAT,
        MPI_INT,
        MPI_BYTE,
        MPI_BYTE
    };

    MPI_Datatype report_datatype;
    MPI_Type_create_struct(blocks, block_lengths, block_displacements, block_types, &report_datatype);
    MPI_Type_commit(&report_datatype);
    return report_datatype;
}

// ---------------------------------------------------------------------------------------------

/*
    MPI_Datatype tsu_info_datatype()

    Creates the MPI_Datatype required to send and receive the TsuInfo struct
    directly in a message. The user is responsible for freeing this datatype
    after creation.

    Returns:
        The MPI_Datatype which represents the TsuInfo struct. 

*/
MPI_Datatype tsu_info_datatype()
{
    // The struct we want to create a datatype for...
    /*
        typedef struct {
            int node;
            int neighbour_count;
            int neighbour_nodes[4];
            int coord[2];
            char ipv4[16];
            char processor[MPI_MAX_PROCESSOR_NAME];
        } TsuInfo;
    */
   
    // Number of members within the struct
    const int blocks = 6;

    // The length of each member in terms of its datatype. 
    const int block_lengths[6] = {
        1, 1, 4, 2, 16, MPI_MAX_PROCESSOR_NAME   
    };

    // Memory offsets of all blocks/members within the struct
    const MPI_Aint block_displacements[6] = {
        offsetof(TsuInfo, node),
        offsetof(TsuInfo, neighbour_count),
        offsetof(TsuInfo, neighbour_nodes),
        offsetof(TsuInfo, coord),
        offsetof(TsuInfo, ipv4),
        offsetof(TsuInfo, processor)
    };

    // Types of each block/member of the struct.  
    const MPI_Datatype block_types[6] = {
        MPI_INT,
        MPI_INT,
        MPI_INT,
        MPI_INT,
        MPI_CHAR,
        MPI_CHAR,
    };

    MPI_Datatype info_datatype;
    MPI_Type_create_struct(blocks, block_lengths, block_displacements, block_types, &info_datatype);
    MPI_Type_commit(&info_datatype);
    return info_datatype;
}

// -----------------------------------------------------------------------------

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
int tsu_distribute_sensors(const int m, const int n, MPI_Comm* tsu_comm, int* min_tsu_world_rank, int* max_tsu_world_rank) 
{
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // If the world size is less than the number of ranks we need, return with error
    if(world_size < m * n + 1)
        return -1;

    // Rank 0 is always the base station, so that means that the min sensor world rank is
    // always 1. We need m * n sensor ranks, which means that the max rank is always m * n (inclusive*). 
    *min_tsu_world_rank = 1;
    *max_tsu_world_rank = m * n; 

    // Flag which represents whether the current process's rank is one which should be
    // included in the cartesian communicator. This avoids including the basestation
    // or extra unrequired processes. 
    bool is_tsu_rank = between(world_rank, *min_tsu_world_rank, *max_tsu_world_rank);

    // We need to split the communicator such that we end up with a communicator 
    // consisting of only sensor nodes, and one with the rest of the ranks. 
    MPI_Comm_split(MPI_COMM_WORLD,is_tsu_rank, world_rank, tsu_comm);

    // If we are not a sensor rank, then we exit now and return no error, with tsu_comm as the communicator that doesn't have sensor nodes
    if(!is_tsu_rank)
        return 0;

    // If we are here, we have a communicator of tsunameters only, however we want to re-organise them 
    // into a communicator which has Cartesian topology, as this mimics their layout.
    MPI_Comm tsu_cart_comm;
    const int dims[2] = {n, m};
    const int wrap[2] = {false ,false};  
    MPI_Cart_create(*tsu_comm, 2, dims, wrap, false, &tsu_cart_comm);

    // Update our output communicator to the Cartesian one
    *tsu_comm = tsu_cart_comm;

    return 0;
}

// -----------------------------------------------------------------------------

/*
    void find_neighbour_ranks(const MPI_Comm tsu_cart_comm, int neighbour_ranks[4], int* neighbour_count)

    Finds the ranks of all neighbour nodes in the given 2D cartesian communicator for the calling process
    which is in said communicator.

    Args:
        tsu_cart_comm: the 2D cartesian communicator
        neighbour_ranks: the output array which will be filled with the ranks of the neighbours
        neighbour_count: the  output number of neighbours the node has; also the number of 
            ranks which have been written into the neighbour_ranks array starting
            at index 0.  
*/
void find_neighbour_ranks(const MPI_Comm tsu_cart_comm, int neighbour_ranks[4], int* neighbour_count)
{
    // Fill neighbour test array with the ranks of all neighbour
    int neighbour_test[4] = {-1};
    MPI_Cart_shift(tsu_cart_comm, 0, 1, &neighbour_test[0], &neighbour_test[1]);
    MPI_Cart_shift(tsu_cart_comm, 1, 1, &neighbour_test[2], &neighbour_test[3]);

    // If a neighbour didn't exist, then its value within the neighbour test array
    // will be MPI_PROC_NULL. So we want to go through the test array and add only the real 
    // neighbours to the output neighbour ranks array. 
    *neighbour_count = 0;
    for(int i = 0; i < 4; i++)
        if(neighbour_test[i] != MPI_PROC_NULL)
            neighbour_ranks[(*neighbour_count)++] = neighbour_test[i]; 
}

// -----------------------------------------------------------------------------

/*
    Tsunameter create_tsunameter_context(MPI_Comm tsu_cart_comm)

    Creates and initialises a tsunameter struct, thus creating a new 'context' for a tsunameter.
    The tsunameter created corresponds to the calling process, which is a node in the given
    2D cartesian communicator.  

    Args:
        tsu_cart_comm: the 2D cartesian communicator

    Returns:
        The new tsunameter context. 
*/
Tsunameter create_tsunameter_context(MPI_Comm tsu_cart_comm)
{
    Tsunameter tsunameter = {0};
    tsunameter.runflag = true;
    tsunameter.height_sma = 0.0;
    tsunameter.history_index = 0;
    tsunameter.cart_comm = tsu_cart_comm;
    tsunameter.report_type = tsu_report_datatype();
    MPI_Comm_rank(tsu_cart_comm, &tsunameter.rank);

    // Get the ranks of all the adjacent tsunameters, the 'neighbours'
    find_neighbour_ranks(tsu_cart_comm, tsunameter.neighbour_ranks, &tsunameter.neighbour_count);

    // Initialize the mutex used to protect tsunameter sense data
    pthread_mutex_init(&tsunameter.data_mutex, NULL);

    return tsunameter;
}

// -----------------------------------------------------------------------------

/*
    void destroy_tsunameter_context(Tsunameter* tsunameter)

    Properly terminates a tsunameter context, freeing all initialised resources.

    Args:
        tsunameter: a pointer to an initialised tsunameter context to destroy
*/
void destroy_tsunameter_context(Tsunameter* tsunameter)
{
    pthread_mutex_destroy(&tsunameter->data_mutex);
    MPI_Type_free(&tsunameter->report_type);
}

// -----------------------------------------------------------------------------

/*
    TsuInfo gather_info(Tsunameter* tsunameter)

    Gathers info about a tsunameter node

    Args:
        tsunameter: a pointer to an initialised tsunameter whose info we want gatherd

    Returns:
        A TsuInfo containing info about the given tsunameter
*/
TsuInfo gather_info(Tsunameter* tsunameter)
{
    TsuInfo info = {0};
    info.node = tsunameter->rank; // the node number is simply the tsunameter rank
    info.neighbour_count = tsunameter->neighbour_count;
    MPI_Cart_coords(tsunameter->cart_comm, tsunameter->rank, 2,  info.coord);

    // The neighbour ranks array will need to be copied in
    memcpy(info.neighbour_nodes, tsunameter->neighbour_ranks, sizeof(tsunameter->neighbour_ranks));
    
    // Grab processor name of node by using MPI
    int len = 0;
	MPI_Get_processor_name(info.processor, &len);

    // Grab the ip address of the node by first grabbing the host
    // then finding the ipv4 of the hosts first address. 
    char hostname[512];
    gethostname(hostname, 512);
    struct hostent* host = gethostbyname(hostname);
  
    struct in_addr* first_addr = (struct in_addr*) host->h_addr_list[0];
    char* ipv4 = inet_ntoa(*first_addr);

    // Copy the ipv4 string into the tsu info
    strcpy(info.ipv4, ipv4);

    return info;
}

// -----------------------------------------------------------------------------

/*
    void* alert_procedure(void* tsunameter_ptr)

    Procedure used to handle a tsunami alert from start to finish. That is, this 
    handles querying all neighbours for their heights, comparing the heights, and
    then sending a report to the base station if necessary. This function is designed
    to be ran from another thread. 

    Args:
        tsunameter: a pointer to the initialised tsunameter context
*/
void* alert_procedure(void* tsunameter_ptr)
{
    // Grab the alert time, which is the rough time that the alert procedure was begun.
    // We use a realtime clock so that the timing is consistent between different processes.  
    struct timespec alert_time; 
    clock_gettime(CLOCK_REALTIME, &alert_time); 

    // Grab the tsunameter context passed in to the thread.
    // This is the live context which is being updated, which means
    // that we need to make a copy of the triggering height water 
    // height, so that it doesn't get overwritten by a new one
    // while we are handling the alert procedure. 
    Tsunameter* tsunameter = (Tsunameter*) tsunameter_ptr;

    // To keep things safe, acquire the data mutex when grabbing a copy of the trigger height
    pthread_mutex_lock(&tsunameter->data_mutex);
    float trigger_height = tsunameter->height_sma;
    time_t trigger_time = tsunameter->last_sense_time;
    pthread_mutex_unlock(&tsunameter->data_mutex);

    MPI_Status status[4];
    MPI_Request requests[4];
    float neighbour_heights[4] = {-1};
    
    // Send query messages to all neighbours so that they respond with their sensed heights
    for(int n = 0; n < tsunameter->neighbour_count; n++)
        MPI_Isend(NULL, 0, MPI_BYTE, tsunameter->neighbour_ranks[n], QUERY, tsunameter->cart_comm, &requests[n]);

    // Wait for all sends to be completed
    MPI_Waitall(tsunameter->neighbour_count, requests, status);

    // Wait for response from neighbours
    for(int n = 0; n < tsunameter->neighbour_count; n++)
        MPI_Irecv(&neighbour_heights[n], 1, MPI_FLOAT, tsunameter->neighbour_ranks[n], QUERY_RESPONSE, tsunameter->cart_comm, &requests[n]);
    
    // We need to wait until we have recevied all response from our neighbours. 
    // However, it is possible that a neighbour has experienced a fault or has terminated
    // before they can respond to us, thus leading to getting stuck on the wait. If our neighbour
    // is terminated, or has experienced a fault, then we can assume that we are also terminated
    // or will be terminated soon after the fault detection takes place. So instead of doing a blocking
    // wait for all responses, we periodically test for responses and also for our runflag. If our runflag
    // goes false (we are terminated), we continue with whatever responses we have. This is because we 
    // could be expecting 4 response, but receive 2, which happen to agree. It is important that 
    // a tsunami alert is still sent in this case, even if we are terminating. 
    int response_flag = 0;
    while(!response_flag && tsunameter->runflag)
    {
        // Test for responses, if all are in then the response flag becomes true and we exit the loop. 
        MPI_Testall(tsunameter->neighbour_count, requests, &response_flag, status);

        // Sleep for a few ms so we aren't wasting too much CPU power polling for responses,
        // but also aren't adding a lot of latency to the alert procedure.  
        usleep(QUERY_RESPONSE_POLL_TIME_MS * 1000);
    }

    // Count how many neighbours have agreeing heights with us, 
    // and decide whether to send a report to the base station. 
    int agrees = 0;
    for(int n = 0; n < tsunameter->neighbour_count; n++)
        if(fabs(neighbour_heights[n] - trigger_height) < TSU_HEIGHT_AGREE_TOLERANCE_M) 
            agrees++;

    // Grab the query time, which is the rough time that the neighbour querying was completed.
    // We use a realtime clock so that the timing is consistent between different processes.  
    struct timespec query_time; 
    clock_gettime(CLOCK_REALTIME, &query_time); 

    // If we have have enough agreeing neighbours, then report it to the base station
    if(agrees >= REQ_AGREEING_NEIGHBOURS)
    {
        // Start by actually creating the report to send
        TsuReport message;
        message.reporting_node = tsunameter->rank;
        message.trigger_height = trigger_height;
        message.trigger_timestamp = trigger_time;
        message.alert_timestamp = alert_time;
        message.query_timestamp = query_time;
        message.agreeing_neighbour_count = agrees;
            
        // Copy the heights array over to the report's array
        memcpy(message.neighbour_heights, neighbour_heights, sizeof(message.neighbour_heights));

        // Send the report to the base station
        MPI_Send(&message, 1, tsunameter->report_type, STA_WORLD_RANK, REPORT, MPI_COMM_WORLD);
    }

    return NULL;
}

// -----------------------------------------------------------------------------

/*
    void* station_ping_procedure(void* tsunameter_ptr)

    Procedure used to provide passive fault detection of the tsunameter
    by periodically sending a pinging the base station with the current
    time of the tsunameter. This is designed to run in another thread. 
    It runs until the tsunameter runflag becomes false. 

    Args:
        tsunameter: a pointer to the initialised tsunameter context
*/
void* station_ping_procedure(void* tsunameter_ptr)
{
    Tsunameter* tsunameter = (Tsunameter*) tsunameter_ptr;

    while(tsunameter->runflag)
    {
        // Send ping to base station to let it know we are still alive. 
        // Send the local time of the tsunameter as time_t. 
        // We will send it as a byte array that the other side reconstructs.  
        time_t curr_time = time(NULL);
        MPI_Send(&curr_time, sizeof(curr_time), MPI_BYTE, STA_WORLD_RANK, PING, MPI_COMM_WORLD);

        // Sleep until time for next ping to be sent.
        sleep(TSU_PING_CYCLE_TIME_S);
    }

    return NULL;
}

// -----------------------------------------------------------------------------

/*
    void* listen_procedure(void* tsunameter_ptr)

    Procedure used by the tsunameter to listen to and respond to messages
    from the base station and other tsunameters. This is designed to 
    run in another thread. It runs until the tsunameter runflag becomes false. 

    Args:
        tsunameter: a pointer to the initialised tsunameter context
*/
void* listen_procedure(void* tsunameter_ptr)
{
    // Grab the tsunameter context passed in to the thread.
    // This is the live context which is being updated, so we need to be careful.
    Tsunameter* tsunameter = (Tsunameter*) tsunameter_ptr;

    MPI_Status status[2];
    MPI_Request requests[2];

    // Both these sources we listen to have different communicators, so we need to wait for their messages, 
    // separately. This is done by doing a non-blocking recv on both, and waiting for one to catch a message. 
    const int tsu_index = 0, station_index = 1;
    MPI_Irecv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, QUERY, tsunameter->cart_comm, &requests[tsu_index]);
    MPI_Irecv(NULL, 0, MPI_BYTE, STA_WORLD_RANK, TERMINATE, MPI_COMM_WORLD, &requests[station_index]);

    while(tsunameter->runflag)
    {
        // Wait for any of the two receives to grab something
        int request_index = -1;
        MPI_Waitany(2, requests, &request_index, status);

        // Check if the message is from a tsunameter. 
        // Note that by doing this first we are also proritising responding over terminating, which is a good thing
        if(request_index == tsu_index)
        {
            // This is extra code added here to make it easier to terminate the tsunameter when 
            // its not supposed to and introduce a fault. Think of it as an intentional weakness
            // in the system so that we can break it on command for testing. If we receive a query
            // message from ourselves, then just exit the loop so we terminate the thread. 
            if(status[request_index].MPI_SOURCE == tsunameter->rank)
            {
                // Cancel the base station listen receive, since we are terminating. 
                MPI_Cancel(&requests[station_index]);
                return NULL;
            }

            // If the request index is for the tsunameter request, then we have received 
            // a query message from another tsunameter; so we need to lock the tsunameter
            // data mutex and respond with our height.  
            
            pthread_mutex_lock(&tsunameter->data_mutex);

            MPI_Send(&tsunameter->height_sma, 1, MPI_FLOAT, status[request_index].MPI_SOURCE, QUERY_RESPONSE, tsunameter->cart_comm); 

            pthread_mutex_unlock(&tsunameter->data_mutex);

            // We have received something, so we need to set up a new non-blocking receive for the tsunameter again. 
            MPI_Irecv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, QUERY, tsunameter->cart_comm, &requests[tsu_index]);
        }
        else if(request_index == station_index)
        {
            // If the request index is for the base station request, then we have received 
            // a termination message from the base station. In response, we need to update 
            // the runflag and shut everything down. 

            // Update runflag so that everything terminates
            // Note that don't use a mutex because this is the only thread that writes to it
            // and race conditions are not critical, or have no effect on the operation of
            // the tsunameter. That is, a race condition in the condition check the while loops
            // is indistinguisable from the runflag just being updated at a later time.  
            tsunameter->runflag = false;
        }
    }

    // If we are here then we are terminating.
    // Since we are the ones that update the termination flag, we can guarantee
    // that there is no non-blocking recv waiting for the base station any more
    // However, there is probably one waiting on tsunameters; so we need to 
    // cancel the non-blocking receive. 
    MPI_Cancel(&requests[tsu_index]);

    return NULL;
}

// -----------------------------------------------------------------------------

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
void tsu_simulate(const MPI_Comm tsu_cart_comm, const int m, const int n, const float threshold, const int fault_node) 
{   
    // Create the context data for the tsunameter, which will store all of its data for us
    Tsunameter tsunameter = create_tsunameter_context(tsu_cart_comm);

    // Find all necessary info about the tsunameter, and send it to the base station
    TsuInfo info = gather_info(&tsunameter);
    MPI_Datatype info_datatype = tsu_info_datatype();
    MPI_Send(&info, 1, info_datatype, STA_WORLD_RANK, INFO, MPI_COMM_WORLD);
    MPI_Type_free(&info_datatype);

    int sense_count = 0;
    bool alert_started = false;
    struct timespec cycle_start, cycle_end;

    // Seed PRNG using the rank as an offset so all the tsunameters give different values
    seed_with_offset(tsunameter.rank);
    
    // Start the auxiliary threads that handle tasks such listening procedures and the fault detection ping
    pthread_t alert_thread, listen_thread, station_ping_thread;
    pthread_create(&station_ping_thread, NULL, &station_ping_procedure, (void*)(&tsunameter));
    pthread_create(&listen_thread, NULL, &listen_procedure, (void*)(&tsunameter));

    // This code is not really part of the tsunameter, its just here to allow us to introduce a 
    // fault into a tsunameter for testing. If we are the fault node, then generate a random
    // number between 3 and 6, which will be how long we will wait before terminating ourselves. 
    bool faulty = false;
    float fault_time_s = -1;
    if(fault_node == tsunameter.rank)
    {
        fault_time_s = random_float(3, 6);
        printf("A fault will be introduced at tsunameter node %d in %.3f seconds\n", tsunameter.rank, fault_time_s);
    }

    // Grab start time for the fault introduction    
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time); 

    while(tsunameter.runflag)
    {
        // Get cycle start time for sleep timing later on
        clock_gettime(CLOCK_MONOTONIC, &cycle_start); 

        // Run the tsunameter sensor to measure the water height and log it
        float sensed_height = random_float(TSU_MAX_WATER, TSU_MIN_WATER);
        sense_count++;

        // We need to save the new height to the history and update the simple moving average
        // This data is accessed by multiple threads, hence is protected by mutex
        pthread_mutex_lock(&tsunameter.data_mutex);

        // Add the sensed height to the history in a FIFO/circular manner
        tsunameter.last_sense_time = time(NULL);
        tsunameter.history[tsunameter.history_index] = sensed_height;
        tsunameter.history_index = (tsunameter.history_index + 1) % HIST_POINTS;

        // Calculate simple moving average.
        // If the sense count is less than the amount of history we are storing
        // then only ensure only the part of the history array that is actually
        // filled is used for the calculations.
        const int history_count = (sense_count >= HIST_POINTS) ? HIST_POINTS : sense_count;
        tsunameter.height_sma = average(tsunameter.history, history_count); 
    
        pthread_mutex_unlock(&tsunameter.data_mutex);

        // If the height simple moving average is above the tsunami threshold, 
        // then we want to start the tsunami alert procedure in another thread.
        if(tsunameter.height_sma >= threshold)
        {
            // Wait for any previous threads that are already alerting to terminate before
            // we start another one to handle another alert procedure. 
            if(alert_started)
                pthread_join(alert_thread, NULL);

            // Create a new thread to perform the full alert procedure in the background
            // Here we are passing in a pointer to the tsunameter data context. 
            // WARNING: ensure we do not call this while under protection of the data mutex or
            // we will run into a deadlock !!
            pthread_create(&alert_thread, NULL, &alert_procedure, (void*) (&tsunameter));
            alert_started = true;
        }

        // If we are the faulty node, and we have passed the required amount of time, shut down with no warning
        if(fault_node == tsunameter.rank && elapsed_s(start_time, cycle_end) >= fault_time_s && tsunameter.runflag)
        {
            // We become faulty by updating our runflag to false so that the tsunameter terminates next loop iteration.
            tsunameter.runflag = false;
            faulty = true;

            // All threads should naturally terminate as well, however the listen thread will be
            // blocked waiting for input, so we send ourselves a quick message to unblock it so
            // that it terminates. This is not great to do, but is required to simulate a fault. 
            // Note: this doesn't deadlock because the receive end is non-blocking. 
            MPI_Send(NULL, 0, MPI_BYTE, tsunameter.rank, QUERY, tsunameter.cart_comm);
            printf("Introducing fault...\n");
            break;
        }

        // Get time at the end of that cycle. 
        clock_gettime(CLOCK_MONOTONIC, &cycle_end);  

        // Sleep for whatver time is left to get the next cycle to be TSU_CYCLE_TIME_S seconds from the start of the last. 
        // Only sleep if the sleep time is positive (we aren't late) and we aren't supposed to terminate
        double sleep_time_s = (double)TSU_CYCLE_TIME_S - elapsed_s(cycle_start, cycle_end);
        if(sleep_time_s > 0 && tsunameter.runflag)
            usleep(sleep_time_s * 1e6); 
            
    }

    // Ensure we join any joinable threads we have created
    if(alert_started)
        pthread_join(alert_thread, NULL);

    pthread_join(station_ping_thread, NULL);
    pthread_join(listen_thread, NULL);

    destroy_tsunameter_context(&tsunameter);

    // Send termination confirmation back to base station (if not faulty)
    if(!faulty)
        MPI_Send(NULL, 0, MPI_BYTE, STA_WORLD_RANK, TERMINATE_COMPLETE, MPI_COMM_WORLD);
    
    printf("Tsunameter node %d has terminated\n", tsunameter.rank);
}

// -----------------------------------------------------------------------------

