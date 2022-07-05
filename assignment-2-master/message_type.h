
#ifndef MESSAGE_TYPE_INCLUDE
#define MESSAGE_TYPE_INCLUDE

/*
    This file simply holds the message tag enum which contains
    all the standardised tags which are used in the simulation
    for specific purposes.  
*/

typedef enum
{
    QUERY,                  // Used for prompting a query which requires a response of type QUERY_RESPONSE 
    QUERY_RESPONSE,         // Used for replying to a QUERY message
    INFO,                   // Used for sending information about a node/process
    REPORT,                 // Used for sending a data report
    PING,                   // Used for pinging another process
    TERMINATE,              // Used for signalling termination of another process
    TERMINATE_COMPLETE,     // Used for confirming the termination of a process
} MessageTags;

#endif