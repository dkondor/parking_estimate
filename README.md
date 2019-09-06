# parking_estimate
Simulation code used for paper:  
Kondor D, Zhang H, Tachet R, Santi P, Ratti C (2018).  
Estimating savings in parking demand using shared vehicles for home-work commuting.  
*IEEE Transactions on Intelligent Transportation Systems.*  
https://doi.org/10.1109/TITS.2018.2869085  
https://arxiv.org/abs/1710.04983

## Overview

This repository contains code to simulate demand for parking and cars during hypothetical commuting trips made at randomized times of day. The methodology is described in more detail in our paper. Each cpp file is a separate program that is suitable for a slightly different variation of this estimation. Thus, each file can be compiled and run separately. These programs require a C++11 compliant compiler and the Boost geometry library (https://www.boost.org/doc/libs/1_71_0/libs/geometry/doc/html/index.html). The code was tested and ran on Linux (Ubuntu 16.04) with GCC and Boost version 1.58.

## Inputs

All programs take their main input as a list of home and work location pairs and corresponding travel times. Agents are assumed to have a daily commute between these locations that takes the travel time given. Commute timings are chosen randomly: either as uniformly random in a given interval or according to an empirical distribution of activity.

Common command line parameters are the following:
* -i <file> -- main input file name; format should be home longitude, home latitude, work longitude, work latitude, home->work travel time, work->home travel time; fields separated by blanks, travel time should be in seconds (read from standard input if not given)
* -s <seed> -- seed used for random number generation (current time is used by default)
* -d <days> -- run the simulation for this many days
* -D <filename> <morning_start> <morning_end> <evening_start> <evening_end> -- if given, read empirical distribution of commute activity from the file specified; the format should be timestamp, frequency (separated by blank, timestamp is in seconds); the next four parameters (if given) specify the start and end of the morning and afternoon peak periods, i.e. commutes are limited to these periods (default: morning between 6am and 9am; afternoon between 5pm and 8pm)
* -w <morning_length> <afternoon_length> -- length of morning and afternoon commute periods; only used with uniform commute timings, i.e. if the previous parameter is not present (default is two hours for both)
* -c <lon> <lat> -- coordinates of the center of the area of the simulation; this is important since this is used to convert coordinates to a local Euclidean projection for spatial indexing
* -C -- if present, the previous is not used, and coordinates are expected in meters from a reference point
* -r <distance> -- maximum distance of parking from destination (main simulation parameter, in meters; default is 500m)
* -o <file> <day> -- if present, a detailed record of the given day of the simulation is written in the file specified here
* -0 -- if given, all travel times are treated as zero, the result is a limit on efficiency (due to imbalances in flows)

## Output

Main output is:  
day, n_cars, n_parking, distance  
i.e. the number of days elapsed, the number of cars and parking spaces used respectively and the extra distance traveled between users' trip start and end location and parking. Distance is calculated as Euclidean distance, thus is a lower bound on real distances.

## Program usage and differences

### parkestimate_boost.cpp
Basic implementation, it is assumed that people are driving themselves (as opposed to using self-driving cars) and are either using shared or private cars. All parking spaces are shared. Extra options:
* -p -- if given, everyone uses their own private car; otherwise, everyone uses shared cars
* -O <file> -- instead of running the simulation, just generate trip start times and output the resulting trips to this file
* -I <file> -- instead of generating trip start times, read the previously generated trips from the given file

### parkestimate_boost_sd.cpp
Implementation with self-driving cars in mind, it is expected that maximum distances from parking can be higher, thus the time it takes to find these is explicitely taken into account. Extra option:
* -S <speed> -- (average) speed that self-driving cars are able to travel when going to park or coming from parking to meet a passenger

### parkestimate_boost_sd_pmax.cpp
Similarly, assumes self-driving cars, and also uses a strict limit for the maximum number of parking spaces and cars. If the limit is is reached, cars will travel further then r_max to find parking. Extra options:
* -P <np> -- maximum allowed number of parking spaces
* -C <nc> -- maximum allowed number of cars

### parkestimate_boost_sd_pmax2.cpp
Similarly, assumes self-driving cars, but distributes a fixed number of cars and parking at the beginning of the day and it is not allowed to add more. Exits with an error if a trip cannot be served. Extra options:
* -P <np> -- fixed number of parking spaces
* -C <nc> -- fixed number of cars





