#include "base_station.h"

#include <math.h>
#include <mpi.h> 
#include <stdbool.h>
#include <stdlib.h> 
#include <time.h>
#include <unistd.h> 
#include <stdio.h> 

#include "message_type.h"
#include "tsunameter.h"
#include "satellite.h"
#include "maths_util.h"
#include "simulation.h"

// ---------------------------------------------------------------------------------------------
//    PRIVATE CONSTANTS
// ---------------------------------------------------------------------------------------------

#define REPORT_ARRAY_SIZE 100 // The size of the pending reports queue in the BaseStation struct below

// ---------------------------------------------------------------------------------------------
//    PRIVATE STRUCTS
// ---------------------------------------------------------------------------------------------

// The PendingReport struct represents a report which has been received and is pending to be processed
typedef struct {
    // The report received from a Tsunameter
    TsuReport report;

    // The time that the report was received
    struct timespec timestamp;
} PendingReport;

// The TsuMetrics struct is used to record specific performance metrics about
// the entire tsunami detection system on a per tsunameter node basis.  
typedef struct {
    // The number of matching reports sent by this node
    int matching_reports;

    // The number of mismatched reports sent by this node
    int mismatched_reports;

    // The total number of messages sent/received by this 
    // node in respect to the generation of its reports. 
    int tsunami_report_messages;

    // The total number of messages sent by this node with
    // the aim of fault detection only. 
    int fault_detection_messages;
} TsuMetrics;

// The BaseStation struct represents an instance of a base station, hence 
// holds all data and information pertaining to the base station. 
typedef struct {
    // A flag representing whether the base station should keep running or not
    // If false, the base station will shut down. 
    bool runflag;

    // A queue of received reports which need to be processed by the base station 
    PendingReport pending_reports[REPORT_ARRAY_SIZE];

    // The next index to write to in the pending reports queue. 
    // Also equivalent to the count of reports in the queue. 
    int write_index;
    
    // A flag which asserts whether the queue has pending
    // reports which need to be processed or not. 
    bool report_available;

    // The mutex which protects all report processing variables such as
    // the pending_reports queue, the report_available flag and write_index
    pthread_mutex_t report_data_mutex;

    // The starting world rank of the continuous chunk of ranks
    // used for the cartesian communicator of the tsunameter nodes.    
    int min_tsu_world_rank;

    // The (inclusive) ending world rank of the continuous chunk of ranks
    // used for the cartesian communicator of the tsunameter nodes.    
    int max_tsu_world_rank;
    
    // The number of rows in the cartesian communicator
    int m;

    // The number of columns in the cartesian communicator
    int n;

    // The number of tsunameter nodes contained the cartesian communicator
    // This is simply a quick reference to m*n
    int node_count;

    // A pointer to a dynamically allocated array of size m*n, 
    // which holds all the last ping times of all tsunameter nodes
    // for use in the passive fault detection system. 
    // Array indices correspond to the tsunameter node numbers/ranks 
    // within their respective cartesian communicator
    time_t* ping_times; 

    // A pointer to a dynamically allocated array of size m*n, 
    // which holds a flag to assert which nodes are faulty and which aren't. 
    // Array indices correspond to the tsunameter node numbers/ranks 
    // within their respective cartesian communicator
    bool* fault_nodes; 
    
    // The mutex which protects all fault detection & ping variables
    // such as the ping_times and fault_nodes arrays. 
    pthread_mutex_t ping_data_mutex;

    // The 'lifetime' count of matching reports which have been processed
    int matching_reports;

    // The 'lifetime' count of mismatched reports which have been processed
    int mismatched_reports;

    // The total number of messages sent in the whole system by both the 
    // base station and tsunameters in regards to sent tsunami reports only. 
    int total_tsunami_report_messages; 

    // The total number of messages sent in the whole system by both the 
    // base station and tsunameters in the aim of fault detection only. 
    int total_fault_detection_messages;
    
    // The total amount of seconds spent between every tsunami alert triggering
    // in the tsunameter and the point that the report is received in the base station. 
    double total_comm_time_s;
    
    // The total amount of seconds spent between every tsunami alert triggering
    // in the tsunameter and the tsunameter getting all responses from its neighbours
    // about their recorded heights. 
    double total_query_comm_time_s;

    // The total amount of seconds spent between every tsunami alert triggering
    // in the tsunameter and the base station actually iterating and processing
    // the received report. This is a metric on processing latency, it doesn't
    // really say anything about network performance. 
    double total_process_time_s; 

    // A pointer to a dynamically allocated array of size m*n, 
    // which holds TsuMetrics structs containing performance metrics
    // for each of the tsunameters.
    // Array indices correspond to the tsunameter node numbers/ranks 
    // within their respective cartesian communicator
    TsuMetrics* node_metrics;

    // A pointer to a dynamically allocated array of size m*n, 
    // which holds TsuInfo structs containing useful information
    // about each of the tsunameter nodes. 
    // Array indices correspond to the tsunameter node numbers/ranks 
    // within their respective cartesian communicator
    TsuInfo* node_info;

} BaseStation;

// This enum simply stores different possible reasons
// that prompt system termination for the base station. 
typedef enum {
    ITERATIONS_REACHED,
    SENTINEL_ASSERTED,
    FAULT_DETECTED
} TerminationReason;

// ---------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------------------------------

/*
    BaseStation create_base_station_context(int m, int n, int min_tsu_world_rank, int max_tsu_world_rank)

    Creates and initialises a BaseStation struct, thus creating a new 'context' for a base station. 
    Only one instance of the base station should ever exist so this should only be called once. 
    The function only returns after communicating with all newly created tsunameters, so calling
    it another time would cause this function to never return. 

    Args:
        m: the number of rows in the cartesian communicator
        n: the number of columns in the cartesian communicator
        min_tsu_world_rank: The starting world rank of the continuous chunk of ranks used for the cartesian communicator of the tsunameter nodes.    
        max_tsu_world_rank: The (inclusive) ending world rank of the continuous chunk of ranks used for the cartesian communicator of the tsunameter nodes.    

    Returns:
        The new base station context. 
*/
BaseStation create_base_station_context(int m, int n, int min_tsu_world_rank, int max_tsu_world_rank)
{
    // Create the station struct, initialising it to 0 initially.
    // There are many members in the struct which should be set to 0
    // we avoid having to do that by manually initialising the struct to 0. 
    BaseStation station = {0};
    station.runflag = true;
    station.m = m;
    station.n = n;
    station.node_count = m * n; 
    station.min_tsu_world_rank = min_tsu_world_rank;
    station.max_tsu_world_rank = max_tsu_world_rank;
    
    // Initialise all the dynamic arrays which store node data
    station.fault_nodes = malloc(sizeof(bool) * station.node_count);
    station.ping_times = malloc(sizeof(time_t) * station.node_count);
    station.node_info = malloc(sizeof(TsuInfo) * station.node_count);
    station.node_metrics = malloc(sizeof(TsuMetrics) * station.node_count);
    
    // Gather all tsunameter node info. The tsunameters send their
    // info to the base station as soon as they start up, so we 
    // should get a message from them now. 
    for(int r = min_tsu_world_rank; r <= max_tsu_world_rank; r++)
    {
        // Change world rank to actual node numbers (AKA cartesian comm ranks)
        // We use the node number as the index in the array to reference the node. 
        int node = r - min_tsu_world_rank;
        
        MPI_Status status;
        MPI_Datatype info_type = tsu_info_datatype();
        MPI_Recv(&station.node_info[node], 1, info_type, r, INFO, MPI_COMM_WORLD,  &status);
        MPI_Type_free(&info_type);
    }

    // Set all default values in the dynamically created arrays
    time_t curr_time = time(NULL);
    TsuMetrics default_metrics = {0}; 
    for(int i = 0; i < station.node_count; i++)
    {
        // Set default metrics, which are all initialised to 0
        station.node_metrics[i] = default_metrics;

        // Set last ping to current time otherwise
        // all nodes will immediately fault as the last
        // ping time would refer back to 1970.
        station.ping_times[i] = curr_time;
        
        // Set not faulty
        station.fault_nodes[i] = false;
    }

    // Init both the required mutexes
    pthread_mutex_init(&station.ping_data_mutex, NULL);
    pthread_mutex_init(&station.report_data_mutex, NULL);

    return station;
}

// ---------------------------------------------------------------------------------------------

/*
    void destroy_base_station_context(BaseStation* base_station)

    Properly terminates a base station context, correctly freeing all
    dynamic memory and other initialised resources. 

    Args:
        station: a pointer to an initialised base station context to destroy
*/
void destroy_base_station_context(BaseStation* station)
{
    pthread_mutex_destroy(&station->ping_data_mutex);
    pthread_mutex_destroy(&station->report_data_mutex);
    free(station->node_metrics);
    free(station->fault_nodes);
    free(station->ping_times);
    free(station->node_info);
}

// ---------------------------------------------------------------------------------------------

/*
    SatelliteData find_corresponding_satellite_reading(int required_coord[2], time_t required_time, bool* data_found) 

    Searches through the satellite data array to find data corresponding to the given coord and required time. 

    Args:
        required_coord: the coordinates that the data must be for
        required_time: the time that the data must be for, within some allowed tolerance
        data_found: an output flag asserting whether data was successfully found for the 
            given requirements. It is possible that the required time is too old for 
            data to exist anymore, hence nothing is found. 

    Returns:
        The corresponding satellite data. If data_found is false, then the return value is undefined. 
*/
SatelliteData find_corresponding_satellite_reading(int required_coord[2], time_t required_time, bool* data_found) { 
    SatelliteData data = {0};
    bool found = false;

    // lock the satellite mutex so that nothing is 
    // written to the data while we are trying to read from it
    pthread_mutex_lock(&g_satellite.satellite_mutex);

    // This is the offset applied to the satellite data array's write index which 
    // determines the index of the data we will check in each loop iteration    
    int read_offset = 0;

    // We will loop through the satellite data array backwards from its writing index, until
    // we either find the corresponding data, look through the whole array, or the data 
    // we are searching is too old to serve useful so we exit early.  
    while(!found)
    {
        // If the read_offset is greater than the amount of data stored in the array
        // then we have covered all the data there is, meaning we have not found
        // satellite data that is suitably within the required time frame.  
        if(read_offset >= g_satellite.data_count)
        {
            // In this case return whatever satellite data we have, but
            // set the data_found flag as false to indicate that we didn't find anything
            *data_found = false;

            // End the loop by pretending we found something
            found =true;
        }

        // Calculate index for next value to check.
        // We iterate backwards from the write index so that we always get the latest data
        // Use a simple circular buffer index calculation to calculate the index, however
        // this time we go backwards, meaning we can end up with negative numbers before
        // the % operator. This make us end up with a negative index, so we add the 
        // array length to the value to always keep it positive. The array length 
        // value gets modulo'd away so it doesn't affect that calculation. 
        // For example: -5 % 10 = -5, but (-5 + 10) % 10 = 5
        int index = (g_satellite.write_index - read_offset - 1 + SAT_ARRAY_LENGTH) % SAT_ARRAY_LENGTH;

        // Grab satellite data at that index
        data = g_satellite.data_array[index];

        // If data is not for our current coord, then continue and check the next index
        if(data.coordinates[0] != required_coord[0] || data.coordinates[1] != required_coord[1])
        {
            read_offset++;
            continue;
        }

        // If we are here, then that data is for the required coord. 

        // Get time difference between the data's time and the required time
        double time_diff_s = difftime(required_time, data.timestamp); // time1 - time2

        if(abs(time_diff_s) <= STA_SAT_TIME_TOLERANCE_S)
        {
            // If the time difference is within the time tolerance of the data, then we accept it.
            // Note that this assumes that the data we are comparing to isn't in the future by more than the time tolerance
            // which should not happen for reasonable values of the STA_SAT_TIME_TOLERANCE_S (at least SAT_CYCLE_TIME_S)
            *data_found = true;
            found = true;
        } 
        else if(time_diff_s >= 2 * STA_SAT_TIME_TOLERANCE_S)
        {
            // If the time difference is positive, meaning that the required time is 
            // in the future compared to the data's timestamp, and its in the future
            // by more than two times the time tolerance, then we break out of the search
            // and assume we did not find any data. We can do this because we iterate 
            // through the array backwards, so the latest data is iterated through 
            // first. If we run into data which is way too old, then the rest of the array
            // will only contain even older data which is not useful. 
            
            // In this case return whatever satellite data we have, but
            // set the data_found flag as false to indicate that we didn't find anything
            *data_found = false;

            // End the loop by pretending we found something
            found =true;
        }
        else 
        {
            // If the data was no good, then we need to find an older one that still works 
            // Since the index is currently at a coordinate data point, and there are
            // m * n coordinate data points per time. We can simply advance by 
            // m * n to jump straight to the next data point for this coord.
            int coord_offset = g_satellite.m * g_satellite.n;
            read_offset += coord_offset;
        } 
    }

    // unlock the satellite mutex
    pthread_mutex_unlock(&g_satellite.satellite_mutex);

    return data;
}   

// ---------------------------------------------------------------------------------------------

/*
    void process_report(BaseStation* station, FILE* log_file, PendingReport pending_report, int iteration)

    Processes a pending report by determining whether its a match or mismatch. Additionaly it logs
    extensive information about the report to the log file, and updates performance metric trackers
    within the station based on the report. 
    Args:
        station: a pointer to the active base station context
        log_file: a pointer to an open file stream to which the report data will be logged to
        pending_report: the pending report being processed
        iteration: the current iteration number of the base station
*/
void process_report(BaseStation* station, FILE* log_file, PendingReport pending_report, int iteration)
{
    time_t curr_time = time(NULL);

    TsuReport report = pending_report.report;

    // Grab node specific data
    TsuInfo node_info = station->node_info[report.reporting_node];
    TsuMetrics* node_metrics = &station->node_metrics[report.reporting_node];

    // Grab the satellite data which best matches the report
    bool data_found, matching = false;
    SatelliteData corresponding_data = find_corresponding_satellite_reading(node_info.coord, report.trigger_timestamp, &data_found);
    
    // Check if the satellite data matches, hence determining if have a true alert or false alert
    matching = data_found && fabs(report.trigger_height - corresponding_data.height) < STA_SAT_HEIGHT_TOLERANCE_M;

    // Update report metrics
    if(matching)
    {
        station->matching_reports++;
        node_metrics->matching_reports++;
    }
    else
    {
        station->mismatched_reports++;
        node_metrics->mismatched_reports++;
    }

    // Update message metrics

    // Each tsunameter sends query messages to its neighbours and gets
    // responses, then sends a single report to the base station. Hence its
    // total number of messages is 2 * neighbours + 1. 
    int total_messages = 2 * node_info.neighbour_count + 1;
    station->total_tsunami_report_messages += total_messages;
    node_metrics->tsunami_report_messages += total_messages;

    // Update timing metrics
    double query_comm_time_s = elapsed_s(report.alert_timestamp, report.query_timestamp);
    station->total_query_comm_time_s += query_comm_time_s;
    
    double report_comm_time_s = elapsed_s(report.alert_timestamp, pending_report.timestamp);
    station->total_comm_time_s += report_comm_time_s;

    // Consider the report as being officially processed now, 
    // so take the timestamp for the processing now. 
    struct timespec report_processing_timestamp;
    clock_gettime(CLOCK_REALTIME, &report_processing_timestamp); 

    double report_processing_time_s = elapsed_s(report.alert_timestamp, report_processing_timestamp);
    station->total_process_time_s += report_processing_time_s;

    // With the report fully processed, log everything to file
    fprintf(log_file, "REPORT\n");
    fprintf(log_file, "****************************************************************\n\n");
    
    // Initial report info
    fprintf(log_file, "Iteration: %d\n", iteration);
    fprintf(log_file, "Logged time: %s\n", ctime(&curr_time));
    fprintf(log_file, "Alert trigger time: %s", ctime(&report.trigger_timestamp));

    // Whether the report was a true or false alert
    if(matching)
        fprintf(log_file, "Alert type: %s\n", "MATCH\n");
    else
        fprintf(log_file, "Alert type: %s\n", "MISMATCH\n"); 

    // Reporting node info
    fprintf(log_file, "%-20s %-10s %-14s %s\n", "Reporting Node", "Coord", "Height (m)", "IPv4");
    fprintf(log_file, "%-20d (%d,%d)   %14f    %-16s (%s)\n\n", node_info.node, node_info.coord[0], node_info.coord[1], report.trigger_height, node_info.ipv4, node_info.processor);
    
    // Neighbour node info
    fprintf(log_file, "%-20s %-10s %-14s %s\n", "Neighbour Node", "Coord", "Height (m)", "IPv4");
    for(int j = 0; j < node_info.neighbour_count; j++){
        TsuInfo neighbour_info = station->node_info[node_info.neighbour_nodes[j]];
        fprintf(log_file, "%-20d (%d,%d)   %14f    %-16s (%s)\n",neighbour_info.node , neighbour_info.coord[0], neighbour_info.coord[1], report.neighbour_heights[j], neighbour_info.ipv4, neighbour_info.processor);
    }

    // Agreeing neighbour node info
    fprintf(log_file,"\nNumber of agreeing neighbour nodes: %d\n", report.agreeing_neighbour_count);
    fprintf(log_file,"Max. tolerance range between node readings (m): %d\n\n", TSU_HEIGHT_AGREE_TOLERANCE_M);

    // Satellite info, only print if we actually got corresponding data from the satellite
    if(data_found)
    {
        fprintf(log_file,"Satellite altimeter reporting coord: (%d,%d)\n", corresponding_data.coordinates[0], corresponding_data.coordinates[1]);
        fprintf(log_file,"Satellite altimeter reporting time: %s", ctime(&corresponding_data.timestamp));
        fprintf(log_file,"Satellite altimeter reporting height: %f\n", corresponding_data.height);
        fprintf(log_file,"Max. tolerance range between altimeter and node readings (m): %d\n\n", STA_SAT_HEIGHT_TOLERANCE_M);
    }
    else
    {
        // Handle case where we didnt get corresponding data. 
        // This really should never happen with correct operation of the whole system. 
        fprintf(log_file,"No satellite altimeter data was found for this coord within the time tolerance of %d seconds\n\n", STA_SAT_TIME_TOLERANCE_S);
    }

    // Log the timing metrics
    fprintf(log_file,"Node alert triggered to neighbour query complete time (seconds): %f\n", query_comm_time_s);
    fprintf(log_file,"Node alert triggered to base station report received time (seconds): %f\n", report_comm_time_s);
    fprintf(log_file,"Node alert triggered to base station report fully processed time (seconds): %f\n", report_processing_time_s);

    // Log message metrics
    fprintf(log_file,"Total messages sent by tsunameter nodes in order to achieve the report: %d\n", total_messages);
    fprintf(log_file, "\n----------------------------------------------------------------\n");
}

// ---------------------------------------------------------------------------------------------

/*
    void log_fault(BaseStation* station, FILE* log_file, time_t ping_time, time_t detection_time, int faulty_node, int iteration)

    Logs information about a detected fault to the given file stream

    Args:
        station: a pointer to the active base station context
        log_file: a pointer to an open file stream to which the statistics data will be logged to
        ping_time: the ping time which lead to raising the fault
        detection_time: the time that the fault was actually detected by the base station
        faulty_node: the node which was detected to be faulty, this is its rank with respect to the 
            cartesian communicator. It is not its world rank. 
        iteration: the current iteration number of the base station
*/
void log_fault(BaseStation* station, FILE* log_file, time_t ping_time, time_t detection_time, int faulty_node, int iteration)
{
    time_t curr_time = time(NULL);

    // Grab the information we know about the faulty node
    TsuInfo node_info = station->node_info[faulty_node];

    // Initial logging info
    fprintf(log_file, "FAULT DETECTED\n");
    fprintf(log_file, "****************************************************************\n\n");
    fprintf(log_file, "Iteration %d\n", iteration);
    fprintf(log_file, "Logged time: %s\n", ctime(&curr_time));

    // Log faulty node specific info
    fprintf(log_file, "Faulty node: %d\n", faulty_node);
    fprintf(log_file, "Node coord: (%d,%d) \n", node_info.coord[0], node_info.coord[1]);
    fprintf(log_file, "Node IPv4: %s (%s)\n\n", node_info.ipv4, node_info.processor);

    // Add ping cycle time to the last ping time to get the time 
    // that we expected the next ping from the faulty node
    struct tm expected_ping_time = *localtime(&ping_time);
    expected_ping_time.tm_sec += TSU_PING_CYCLE_TIME_S; 
    mktime(&expected_ping_time);

    // Log faulty specific timing data
    fprintf(log_file, "Last ping occurred at: %s", ctime(&ping_time));
    fprintf(log_file, "Expected ping at: %s", asctime(&expected_ping_time));
    fprintf(log_file, "Detected fault at: %s", ctime(&detection_time));

    // Log the full ping table, hence ping times for all nodes
    fprintf(log_file, "\nFull Ping Status\n");
    fprintf(log_file, "%-8s %-10s %-32s %s\n", "Node", "Coord", "IPv4/Processor" , "Last Ping");
    for(int i = 0; i < station->node_count; i++)
    {
        TsuInfo info = station->node_info[i];
        fprintf(log_file, "%-8d (%d,%d) %16s/%-20s %s",i , info.coord[0], info.coord[1], info.ipv4, info.processor,  ctime(&station->ping_times[i]));
    }

    fprintf(log_file, "\nPrompting termination of system...\n");
    fprintf(log_file, "\n----------------------------------------------------------------\n");
}

// ---------------------------------------------------------------------------------------------

/*
    void log_final_stats(BaseStation* station, FILE* log_file, int iteration)

    Grabs all the compiled metrics found in the given station and logs a summary
    to the provided file stream
    

    Args:
        station: a pointer to the active base station context
        log_file: a pointer to an open file stream to which the statistics data will be logged to
        iteration: the current iteration number of the base station
*/
void log_final_stats(BaseStation* station, FILE* log_file, int iteration)
{
    fprintf(log_file, "SESSION STATISTICS\n");
    fprintf(log_file, "****************************************************************\n\n");

    // Log 'global' metrics
    const int total_reports = station->mismatched_reports + station->matching_reports;
    fprintf(log_file,"Total Iterations: %d\n", iteration);
    fprintf(log_file,"Total Number of Reports: %d\n", total_reports);
    fprintf(log_file,"Total Number of Matches: %d\n", station->matching_reports);
    fprintf(log_file,"Total Number of Mismatches: %d\n\n", station->mismatched_reports);
    fprintf(log_file,"Avg. Neighbour Query Communication Time (s): %f\n", station->total_query_comm_time_s / total_reports);
    fprintf(log_file,"Avg. Report Communication Time (s): %f\n", station->total_comm_time_s / total_reports);
    fprintf(log_file,"Avg. Report Processing Time (s): %f\n\n", station->total_process_time_s / total_reports);
    fprintf(log_file,"Total Tsunami Detection Messages: %d\n", station->total_tsunami_report_messages);
    fprintf(log_file,"Total Fault Detection Messages: %d\n", station->total_fault_detection_messages);

    // Log the 'per node' metrics
    fprintf(log_file, "%-8s %-10s %-10s %-20s %-35s %-35s\n", "Node", "Coord", "Reports", "Matching Reports" , "Total Tsunami Report Messages", "Total Fault Detection Messages");
    for(int i = 0; i < station->node_count; i++)
    {
        TsuInfo node_info = station->node_info[i];
        TsuMetrics node_metrics = station->node_metrics[i];
        int total_node_reports = node_metrics.matching_reports + node_metrics.mismatched_reports;

        fprintf(log_file, "%-8d (%d,%d)      %-10d %-20d %-35d %-35d\n",i , node_info.coord[0], node_info.coord[1], total_node_reports, node_metrics.matching_reports,  node_metrics.tsunami_report_messages, node_metrics.fault_detection_messages);
    }
    fprintf(log_file, "\n----------------------------------------------------------------\n");
}

// ---------------------------------------------------------------------------------------------

/*
    void terminate_system(BaseStation* station, FILE* log_file, TerminationReason reason, pthread_t satellite_thread, pthread_t listen_thread, pthread_t sentinel_thread, time_t station_start_time)

    Properly and gracefully terminates the full system, also logging the termination to the provided file stream

    Args:
        station: a pointer to the active base station context
        log_file: a pointer to an open file stream to which the statistics data will be logged to
        reason: the reason that we are performing a termination, this will change how we handle the termination
        satellite_thread: the active satellite thread which is simulating the satellite altimeter
        listen_thread: the active listener thread which is running the node_listening_procedure
        sentintel_thrread: the listener thread which is (or isn't) running the sentinel procedure
        station_start_time: a timestamp to the time that the station started being simulated
*/
void terminate_system(BaseStation* station, FILE* log_file, TerminationReason reason, pthread_t satellite_thread, pthread_t listen_thread, pthread_t sentinel_thread, time_t station_start_time )
{
    // Ensure that the runflag is actually off so that all the threads can properly terminate
    station->runflag = false;

    fprintf(log_file, "SYSTEM TERMINATION\n");
    fprintf(log_file, "****************************************************************\n\n");
    fprintf(log_file, "Starting System Termination...\n\n");

    // Terminate the satellite, this is done by first
    // calling terminate_satellite, to signal its simulation to end
    fprintf(log_file, "Sending termination signal to satellite...\n");
    terminate_satellite();

    // We then wait for the satellite thread to join to confirm termination
    pthread_join(satellite_thread, NULL);
    fprintf(log_file, "Satellite thread joined, termination successful\n\n");

    // Next terminate the tsunameters by sending messages to all nodes which aren't known to be faulty.
    fprintf(log_file, "Sending termination signals to nodes...\n");
    for(int r = station->min_tsu_world_rank; r <= station->max_tsu_world_rank; r++)
    {
        // Convert the world ranks of the nodes to their ranks respective to the
        // cartesian communicator, which is how the faulty nodes array is indexed. 
        int node = r - station->min_tsu_world_rank;
        if(!station->fault_nodes[node])
        {
            // If not faulty, send the termination message
            MPI_Send(NULL, 0, MPI_BYTE, r, TERMINATE, MPI_COMM_WORLD);
            fprintf(log_file, "Termination sent to node %d\n", node);
        }
        else
        {
            // If faulty, simply log that we are skipping sending the termination signal
            fprintf(log_file, "Termination skipped for faulty node %d\n", node);
        }
    }

    // Wait for termination confirmations from non-faulty tsunamter nodes, and log the confirmations
    MPI_Status status;
    bool any_confirmation_received = false;
    fprintf(log_file, "\nWaiting for confirmations...\n");
    for(int r = station->min_tsu_world_rank; r <= station->max_tsu_world_rank; r++)
    {
        // Convert the world ranks of the nodes to their ranks respective to the
        // cartesian communicator, which is how the faulty nodes array is indexed. 
        int node = r - station->min_tsu_world_rank;
        if(!station->fault_nodes[node])
        {
            // If not faulty, wait for the termination confirmation message
            any_confirmation_received = true;
            MPI_Recv(NULL, 0, MPI_BYTE, r, TERMINATE_COMPLETE, MPI_COMM_WORLD, &status);
            fprintf(log_file, "Node %d has confirmed termination\n", r - station->min_tsu_world_rank);
        }
    }

    // Next terminate the listener thread. The listener thread performs a non-blocking
    // wait for a message to be sent. If the thread has been blocked ever since the runflag
    // was updated to false, then it will still be blocked and the thread can't be joined. 
    // However, the thread should have been unblocked after detecting the termination complete
    // messages we got sent above. In the incredibly rare case that all nodes were faulty
    // hence no messages we received above, we just cancel the thread instead of joining it.
    fprintf(log_file, "\nTerminating base station listener thread...\n");
    if(any_confirmation_received)
        pthread_join(listen_thread, NULL);
    else
        pthread_cancel(listen_thread);
    fprintf(log_file, "Thread joined, base station listener terminated...\n");

    // Next we terminate the sentinel thread. How this is done will depend on the
    // termination reason. If we terminated due to a sentinel value, then we need
    // to join the thread. If we terminated due to reaching the specified iterations
    // then we don't need to terminate it because it was never running. If we terminated
    // due to a faulty node, then we need to cancel the thread to unblock and kill it. 
    if(reason != ITERATIONS_REACHED)
    {
        fprintf(log_file, "\nTerminating base station sentinel thread...\n");
        if(reason == SENTINEL_ASSERTED)
            pthread_join(sentinel_thread, NULL); 
        else if(reason == FAULT_DETECTED)
        // Note that its possible to terminate due to a fault while also using iterations
        // so the sentinel thread was never created. In that case, we just eat the error
        // the the pthread_cancel will produce since it won't affect operation at all. 
            pthread_cancel(sentinel_thread); 
        fprintf(log_file, "Thread joined, base station sentinel listener terminated...\n");
    }

    // Log the total system run time in seconds
    time_t curr_time = time(NULL);
    fprintf(log_file, "\nAll systems properly terminated after being alive for (s): %f\n",difftime(curr_time, station_start_time));
    fprintf(log_file, "\n----------------------------------------------------------------\n");
}

// ---------------------------------------------------------------------------------------------

/*
    void* node_listen_procedure(void* base_station)

    Procedure used by the base station to listen to messages being sent
    to it from tsunameter nodes. It runs until the station runflag becomes false. 

    Args:
        station_ptr: a pointer to the initialised station context
*/
void* node_listen_procedure(void* station_ptr)
{
    //Note that this is a pointer to the live BaseStation context, 
    // so we need to be careful about what we modify
    BaseStation* station = (BaseStation*) station_ptr;
    
    // Create the report datatype we can use to receive
    // the full TsuReport struct directly in a message
    MPI_Datatype report_type = tsu_report_datatype();

    MPI_Status status;
    while(station->runflag)
    {
        // Listen for any message sent to the base station, these should
        // all be from tsunameters because nothing else sends messages.
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // We differentiate messages by their tag
        if(status.MPI_TAG == REPORT)
        {
            // If the message was a report, then we want to add the report 
            // to the pending reports queue for the main thread to process

            // We need to modify the report data so lock its mutex
            pthread_mutex_lock(&station->report_data_mutex);

            // Before we receive anything, double check that the report array
            // isn't full. If it is, then skip grabbing the report for now
            // until the main thread has processed the reports.
            if(station->write_index == REPORT_ARRAY_SIZE)
                continue;

            // Receive the report
            TsuReport report;
            MPI_Recv(&report, 1, report_type, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
            
            // Create a pending report package which also has the timestamp of the time we received the report. 
            PendingReport pending_report;
            pending_report.report = report;
            clock_gettime(CLOCK_REALTIME, &pending_report.timestamp);

            // Push the pending report to the report queue
            station->pending_reports[station->write_index++] = pending_report;
            station->report_available = true;

            pthread_mutex_unlock(&station->report_data_mutex);
        }
        else if(status.MPI_TAG == PING)
        {
            // If the message was a ping message, then we simply 
            // want to save the ping time that was sent to us.

            // Receive the ping message, instead of defining an
            // MPI_datatype for the time_t, we just receive it 
            // directly as a byte array. This should work fine
            // running on modern machines with the same endianness. 
            time_t ping_time;
            MPI_Recv(&ping_time, sizeof(time_t), MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

            // Lock the ping mutex since we are going to modify the ping times array
            pthread_mutex_lock(&station->ping_data_mutex);

            // The array stores the ping times indexed by the actual
            // sensor ranks, but we only have access to the world ranks
            // So use the min sensor rank to convert it to sensor rank indexed by 0. 
            int node = status.MPI_SOURCE - station->min_tsu_world_rank;

            // Update fault detection metrics
            station->total_fault_detection_messages++;
            station->node_metrics[node].fault_detection_messages++;

            // Update the last ping time for the node
            station->ping_times[node] = ping_time;

            pthread_mutex_unlock(&station->ping_data_mutex);
        }
    }

    // Remember to free the report datatype that 
    // we created at the start of the function
    MPI_Type_free(&report_type);

    return NULL;
}

// ---------------------------------------------------------------------------------------------

/*
    void* sentinel_procedure(void* station_ptr)

    Procedure used by the base station to listen for a sentinel value
    to be entered to the terminal in order to start termination. The
    sentinel value itself is just the enter key being pressed in the 
    terminal. This runs until enter is pressed, so should be cancelled
    if used in a thread which is to be stopped.  

    Args:
        station_ptr: a pointer to the initialised station context
*/
void* sentinel_procedure(void* station_ptr)
{
    BaseStation* station = (BaseStation*) station_ptr;

    // Simply wait for sentinel value, and terminate when we get it
    printf("Press enter at any time to terminate system\n\n");
    getchar();
    printf("System termination started... this may take some time\n");

    // If we are here, then we are terminating the system so set the runflag to false
    station->runflag = false;

    return NULL;
}

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
void station_simulate(int m, int n, int min_tsu_world_rank, int max_tsu_world_rank, int req_iterations)
{
    // Create and initialise the context data for the base station
    BaseStation station = create_base_station_context(m, n, min_tsu_world_rank, max_tsu_world_rank);
    
    // Find the start time of the system to get total run time at the end
    time_t start_time = time(NULL);

    // Start all auxiliary threads, including the satellite thread
    int grid_size[2] = {m, n};
    pthread_t satellite_thread, sentinel_thread, listen_thread;
    pthread_create(&satellite_thread, NULL, &satellite_simulate, &grid_size);
    pthread_create(&listen_thread, NULL, &node_listen_procedure, &station);

    // Assert the use_sentinel flag if we will terminate the system based
    // on a sentinel value instead of a number of iterations. We know that
    // this is the expected behaviour when the given required iterations is negative
    bool use_sentinel = req_iterations <= 0;

    // If we are using a sentinel value to terminate, then create the thread that waits for it
    if(use_sentinel)
        pthread_create(&sentinel_thread, NULL, &sentinel_procedure, &station);
    else
        printf("System will run for %d iterations before terminating\n\n", req_iterations);

    // Open the file which will store all logged data
    FILE* log_file = fopen("report.log", "w");

    // This variable will hold the reason we are terminating, whenever we terminate
    // We initialise it to be either the iterations reached or sentinel value 
    // reason which are the default expected reasons. We update it to the fault
    // detected reason later if necessary.  
    TerminationReason termination_reason = use_sentinel ? SENTINEL_ASSERTED : ITERATIONS_REACHED;

    int iteration = 0;
    struct timespec iteration_start_time, iteration_end_time;
    while(station.runflag)
    {
        // Grab iteration start time for sleep timing later on
        clock_gettime(CLOCK_MONOTONIC, &iteration_start_time); 
        
        // First thing we do in an iteration is process all available reports

        // Acquire mutex for the report data as we want to process pending reports
        // so don't want reports to be added into the queue at the same time. 
        pthread_mutex_lock(&station.report_data_mutex);

        // Only look through the pending report array if the listener thread flagged
        // that we have reports available to us (it actually pushed pending reports)
        if(station.report_available)
        {
            // 'pop' reports from the queue and process them
            for(int i = 0; i < station.write_index; i++)
            {
                PendingReport report = station.pending_reports[i];
                process_report(&station, log_file, report, iteration); 
            }

            // Reset the write index back to 0 to 'clear' the queue 
            station.write_index = 0;

            // Also reset the flag so we don't try to process reports 
            // until the listener thread adds pending reports again.
            station.report_available = false;
        }

        pthread_mutex_unlock(&station.report_data_mutex);

        // Next we want to go through all the ping times from 
        // the tsunameter nodes to check if any are too old.
        // Nodes with substantially old ping times will be 
        // conisdered faulty.  

        // Acquire mutex for the ping data as we want to go through 
        // the ping times and we want to avoid race conditions. 
        pthread_mutex_lock(&station.ping_data_mutex);

        // We compare the ping times to the stations current time, so grab the current time
        time_t curr_time = time(NULL);

        // Go through each ping time
        for(int i = 0; i < station.node_count; i++)
        {
            // If any of the ping times are longer than 1.5 pings worth of time, then we consider it at fault
            // We don't want to check for exactly longer than one pings worth of time because we could just 
            // be checking the ping at the same time that the ping is being sent. And two pings worth of 
            // time adds too much latency for fault detection, so we check for longer than 1.5 times as 
            // a compromise. The allowed ping should be increased if the system happened to have a lot
            // of latency when sending messages. 
            time_t ping_time = station.ping_times[i];
            if( difftime(curr_time, ping_time) > 1.5 * TSU_PING_CYCLE_TIME_S)
            {
                // If we are here, then we identified a fault, so log the fault
                log_fault(&station, log_file, ping_time, curr_time, i, iteration);
                
                printf("\nFault detected, terminating... this make take some time\n");
                
                // Mark the node as faulty
                station.fault_nodes[i] = true; 

                // Update the runflag so that we terminate
                // Also update the termination reason to be a fault detection
                termination_reason = FAULT_DETECTED;
                station.runflag = false; 
            }
            else station.fault_nodes[i] = false; // mark node as not faulty 
        }
        pthread_mutex_unlock(&station.ping_data_mutex);

        // Grab iteration end time for sleep timing
        clock_gettime(CLOCK_MONOTONIC, &iteration_end_time);

        // Sleep for whatver time is left to get the next cycle to be STA_ITER_TIME_S seconds from the start of the last. 
        // Only sleep if the sleep time is positive (we aren't late) and we aren't supposed to terminate
        double sleep_time_s = (double)STA_ITER_TIME_S - elapsed_s(iteration_start_time, iteration_end_time);
        if(sleep_time_s > 0 && station.runflag)
            usleep(sleep_time_s * 1e6); 

        // If this next iteration is the one in which we were asked to terminate, then update the runflag to terminate
        iteration++;
        if(iteration == req_iterations)
        {
            printf("\nIterations reached, terminating... this make take some time\n");
            station.runflag = false;
        }
    }

    // Run the full system termination procedure to gracefully terminate it and log it at the same time
    terminate_system(&station, log_file, termination_reason, satellite_thread, listen_thread, sentinel_thread, start_time);

    // Log overall performance statistics
    log_final_stats(&station, log_file, iteration);
    
    // Close the log file stream
    fclose(log_file);

    // If we are here, then we have terminated so properly destroy the base station context
    destroy_base_station_context(&station);
}

