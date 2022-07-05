
#ifndef SIMULATION_INCLUDE
#define SIMULATION_INCLUDE

/*
    This file provides all simulation settings used through out all files.
    Hence, the simulation can be quickly modified from this one location. 
*/


#define STA_WORLD_RANK 0 // The base station world rank

#define SAT_CYCLE_TIME_S 2 // The amount of seconds between satellite cycles which generate data
#define TSU_CYCLE_TIME_S 2 // The amount of seconds between tsunameter cycles which generate data
#define STA_ITER_TIME_S 4 // The amount of seconds between base station iterations
#define TSU_PING_CYCLE_TIME_S 2 // The amount of time between tsunameter fault detection pings to the base station

#define TSU_MIN_WATER 6500 // The minimum generatable water height for the tsunameter
#define TSU_MAX_WATER 7500 // The maximum generatable water height for the tsunameter

#define SAT_MIN_WATER 6800 // The minimum generatable water height for the satellite
#define SAT_MAX_WATER 7500 // The maximum generatable water height for the satellite

#define TSU_HEIGHT_AGREE_TOLERANCE_M 100 // The allowed tolerance in height when comparing a tsunameter's heigh with its neighbour nodes
#define STA_SAT_HEIGHT_TOLERANCE_M 200 // The allowed tolerance in height between a satellite reading and tsunameter report for it to be accepted
#define STA_SAT_TIME_TOLERANCE_S 3 // The allowed tolerance in seconds between a satellite reading and tsunameter report for it to be accepted

#endif