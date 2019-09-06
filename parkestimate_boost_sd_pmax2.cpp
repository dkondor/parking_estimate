/*
 * parkestimate_boost_sd_pmax2.cpp -- do parking space estimate calculations using Boost's spatial index
 * 	addition: include extra time used for reaching the parking space from the final destination
 *  and the start from a parking space
 * 	use heap structures to keep track of these
 * 
 *  use a constant number of cars and parking which are generated at a given ratio of home and work locations
 * 	if these are not sufficient, exit with an error
 * 	otherwise, report the number of extra travel distance due to looking for parking
 * 
 * Copyright 2017 Kondor Dániel <dkondor@mit.edu>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include <vector>
#include <random>
#include <queue>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef std::pair<point, unsigned int> value;


class coordconverter {
	protected:
		double clon;
		double clat;
		double factor;
		coordconverter() { }
	public:
		coordconverter(double clon_, double clat_) {
			clon = clon_;
			clat = clat_;
			factor = cos(M_PI*clat/180.0);
		}
		point operator () (double lon,double lat) {
			int x = (int)round(factor*(lon-clon)*20000000.0/180.0);
			int y = (int)round((lat-clat)*20000000.0/180.0);
			return point(x,y);
		}
};

class tdist_empirical { // class for storing and estimating temporal distribution
	protected:	
		std::vector<unsigned int> times;
		std::vector<unsigned int> cdf;
		std::uniform_int_distribution<unsigned int> dist;
	
	public:
		tdist_empirical() { }

		// read the distribution of frequencies from a given file
		// format is timestamp,frequency
		// cumulative distribution is created on the fly
		// if the minimum and maximum parameters are given, limit the data read in between [tmin,tmax)
		int ReadFreqs(char* fn, unsigned int tmin = 0, unsigned int tmax = 0) {
			unsigned int l = 0;
			times.clear();
			cdf.clear();
			
			FILE* f = fopen(fn,"r");
			if(f == 0) { fprintf(stderr,"tdist::ReadFreqs(): Error opening file %s!\n",fn); return 1; }
			
			while(1) {
				int a;
				unsigned int ts,freq;
				do a = fgetc(f); while(a == ' ' || a == '\t');
				if(a == EOF) break;
				l++;
				if(a == '\n') continue;
				
				ungetc(a,f);
				a = fscanf(f,"%u",&ts);
				if(a != 1) { fprintf(stderr,"tdist::ReadFreqs(): invalid data on input line %u!\n",l); return 2; }
				do a = fgetc(f); while(a == ' ' || a == '\t');
				if(a == '\n' || a == EOF) { fprintf(stderr,"tdist.ReadFreqs(): invalid data on input line %u!\n",l); return 2; }
				ungetc(a,f);
				a = fscanf(f,"%u",&freq);
				if(a != 1) { fprintf(stderr,"tdist::ReadFreqs(): invalid data on input line %u!\n",l); return 2; }
				
				if(tmin > 0 && tmax > 0) if(ts < tmin || ts >= tmax) goto readfreq_endl;
				if(times.size() > 0) if(ts <= times.back()) { fprintf(stderr,"tdist::ReadFreqs(): input not sorted on line %u!\n",l); return 3; }
				times.push_back(ts);
				if(cdf.size() > 0) cdf.push_back(cdf.back() + freq);
				else cdf.push_back(freq);
				
readfreq_endl:
				do a = fgetc(f); while( ! (a == '\n' || a == EOF) );
				if(a == EOF) break;
			}
			
			fclose(f);
			if(times.size() == 0) { fprintf(stderr,"tdist::ReadFreqs(): no data read from the input file %s!\n",fn); return 4; }
			dist = std::uniform_int_distribution<unsigned int>(0,cdf.back()-1);
			
			return 0;
		}

		size_t NRecords() { return times.size(); }
		unsigned int Total() { if(cdf.size() > 0) return cdf.back(); else return 0; }

		unsigned int operator () (std::mt19937& r) {
			if(cdf.size() == 0) throw new std::runtime_error("tdist.GetRandomTS(): no distribution loaded!\n");
			uint x = dist(r);
			unsigned int i = std::upper_bound(cdf.begin(),cdf.end(),x) - cdf.begin();
			return times[i];
		}
};



struct start_event {
	int x;
	int y;
	unsigned int ts;
	int dest_x;
	int dest_y;
	unsigned int ttime;
	unsigned int user_id;
	bool operator < (const start_event& x) const { return ts < x.ts; }
};

struct end_event {
	int x;
	int y;
	unsigned int ts;
	unsigned int user_id;
};

struct end_comparer { // comparer for priority_queue 
	bool operator() (end_event& x, end_event& y) const { return x.ts > y.ts; }
};

// struct to store main results
struct res_struct {
	unsigned int ncars; // number of cars 
	unsigned int nparkspaces; // number of parking spots
	double dist_tot; // total 'extra' distance traveled (i.e. between the start / destination and parking
	double dist_thres; // distance traveled above the set maximum distance threshold
	unsigned int trips_thres; // trip start or end events or where the car has to travel more than the maximum distance threshold
};

/*
 * process a set of commute trips, given in the [seq,end) sequence
 * parameters:
 * 		seq -- (in) iterator to trips to process (forward iterator, dereferencable as start_event
 * 		end -- (in) iterator to the end of sequence or sentinel class (events are processed until seq != end)
 * 		parkspaces_empty -- (in/out) spatial index for empty parking spaces (could be already populated with available parking)
 * 		parkspaces_occupied -- (in/out) spatial index for occupied parking spaces, i.e. parked cars (could be already populated with available cars)
 * 		end_events -- (temp) queue class to be used as work space for storing events (should be empty originally, emptied before returning)
 * 		occupy_events -- (temp) queue class to be used as work space for storing events (should be empty originally, emptied before returning)
 * 		dmax -- (in) maximum distance for comfortable walking -- trips where the parking is further from the destination than this are counted separately
 * 		speed1 -- (in) driving speed (using as the crow flies distance) to calculate the time needed to cover the distance between the parking space
 * 			and trip start / end
 *		out -- (out) output detailed info on all processed events
 * 		res -- (out) store the main result here (only add to the extra distances traveled)
 * return value: 0 if OK, 1 if ran out of cars or parking spaces
 */
template <class it, class se, class index1>
int process_events(it seq, se end, index1& parkspaces_empty, index1& parkspaces_occupied,
		std::priority_queue<end_event, std::vector<end_event>, end_comparer>& end_events,
		std::priority_queue<end_event, std::vector<end_event>, end_comparer>& occupy_events,
		double dmax, double speed1, FILE* out, bool shared, std::vector<point>& end_coords, res_struct& res) {
	
	if( ! (end_events.empty() && occupy_events.empty()) ) throw new std::runtime_error("process_events(): invalid (non-empty) queue objects provided!\n");
	
	std::vector<value> results;
	
	while(1) {
		bool snext = (seq != end);
		bool enext = (end_events.size() > 0);
		bool onext = (occupy_events.size() > 0);
		
		if(snext && enext) if( end_events.top().ts < (*seq).ts ) snext = false;
		if(snext && onext) if( occupy_events.top().ts < (*seq).ts) snext = false;
		
		if(snext) {
			// process start event
			start_event s = *seq;
			// search for an available vehicle
			point p(s.x,s.y);
			unsigned int end_ts = s.ts + s.ttime;
			double dist = 0.0;
			// search for "free" cars around the users's location
			if(shared) {
				bool found = false;
				parkspaces_occupied.query(bgi::nearest(p,1),std::back_inserter(results));
				if(results.size() > 0) found = true; // note: no need to check the distance here, just use the closest one all the time
				if(!found) {
					fprintf(stderr,"process_events(): ran out of cars!\n");
					return 1;
				}
				
				dist = bg::distance(results[0].first,p);
				// remove and add to as empty parking space
				unsigned int tsadd = (unsigned int)round(dist / speed1);
				
				parkspaces_empty.insert(results[0]); // otherwise add it to the available empty parking to be reused later
				if(parkspaces_occupied.remove(results[0]) != 1) {
					throw new std::runtime_error("process_events(): error with remove!\n");
				}
				// add the extra travel time needed to move from the parking space to the user's location
				end_ts += tsadd;
				res.dist_tot += dist;
				if(dist >= dmax) { res.dist_thres += dist; res.trips_thres++; }
			}
			else { // shared == false
				// just add an empty parking space
				results.push_back(std::make_pair(p,0));
				parkspaces_empty.insert(results[0]);
			}
			
			if(out) fprintf(out,"%u\t%d\t%d\tTrue\t%d\t%d\t%f\n",s.ts,s.x,s.y,results[0].first.get<0>(),results[0].first.get<1>(),dist);
			
			// add the end of this trip to the end events to be processed
			end_event e;
			e.x = s.dest_x;
			e.y = s.dest_y;
			e.ts = end_ts;
			if(!shared) e.user_id = s.user_id;
			end_events.push(e);
			
			seq++;
			results.clear();
			continue;
		}
		
		if(enext && onext) if( occupy_events.top().ts < end_events.top().ts ) enext = false;
		
		if(enext) {
			// process trip end event
			const end_event& e = end_events.top();
			point p(e.x,e.y);
			
			// search for free parking spaces around the user's location
			parkspaces_empty.query(bgi::nearest(p,1),std::back_inserter(results));
			bool found = false;
			if(results.size() > 0) found = true;
			if(!found) {
				fprintf(stderr,"process_events(): ran out of parking spots!\n");
				return 1;
			}
			double dist = bg::distance(results[0].first,p);
			
			if(shared) {
				// add to the queue of reserved parking spaces, to be made available later
				end_event o;
				o.x = results[0].first.get<0>();
				o.y = results[0].first.get<1>();
				o.ts = e.ts + (unsigned int)round(dist / speed1);
				occupy_events.push(o);
				if(dist >= dmax) { res.dist_thres += dist; res.trips_thres++; }
			}
			else end_coords[e.user_id] = results[0].first; // shared == false -> just store the end coordinates
				
			
			// remove the empty parking spot
			if(parkspaces_empty.remove(results[0]) != 1) {
				throw new std::runtime_error("process_events(): error with remove!\n");
			}
			res.dist_tot += dist;
			if(out) fprintf(out,"%u\t%d\t%d\tFalse\t%d\t%d\t%f\n",e.ts,e.x,e.y,results[0].first.get<0>(),results[0].first.get<1>(),dist);
			
			end_events.pop();
			results.clear();
			continue;
		}
		
		if(onext) {
			// process event for occupying a reserved parking space (add to the list of available cars
			const end_event& o = occupy_events.top();
			parkspaces_occupied.insert(std::make_pair(point(o.x,o.y),0)); 
			occupy_events.pop();
			continue;
		}
		
		break; // no events left to process
	}
	
	if(seq != end || end_events.size() > 0 || occupy_events.size() > 0) throw new std::runtime_error("process_events(): not all events processed!\n");
	return 0;
}


template <class tdist>
unsigned int do_estimate(std::vector<point>& home_loc, std::vector<point>& work_loc,
			std::vector<std::pair<unsigned int, unsigned int> >& travel_times, std::mt19937& rg, tdist&& morning_dist,
			tdist&& evening_dist, double dmax, bool shared, double rcars, double rpark, double speed1,
			bool travel0, std::vector<res_struct>& res, FILE* out, unsigned int day_out) {
	if(home_loc.size() != work_loc.size() || work_loc.size() != travel_times.size() || home_loc.size() == 0)
		throw new std::runtime_error("do_estimate(): invalid input!\n");
	
	unsigned int nparkspaces = 0;
	unsigned int ncars = 0;
	unsigned int nusers = home_loc.size();
	
	bgi::rtree< value, bgi::rstar<16> > parkspaces_empty;
	bgi::rtree< value, bgi::rstar<16> > parkspaces_occupied;
	std::priority_queue<end_event, std::vector<end_event>, end_comparer> end_events;
	std::priority_queue<end_event, std::vector<end_event>, end_comparer> occupy_events;
	
	// generate parking spaces
	nparkspaces = (unsigned int)ceil(rpark*2*nusers);
	if(shared) ncars = (unsigned int)ceil(rcars*nusers);
	else ncars = nusers;
	if(ncars > nparkspaces) throw std::runtime_error("do_estimate(): ncars > nparkspaces!\n");
	// generate exactly ncars occupied parking and (nparkspaces - ncars) empty parking
	// use random arrays for this
	{
		std::vector<std::pair<unsigned int, unsigned int> > tmp(nusers);
		if(shared) {
			for(unsigned int i=0;i<nusers;i++) tmp[i] = std::make_pair(rg(),i);
			std::sort(tmp.begin(),tmp.end(),[](auto a, auto b) { return a.first < b.first; });
			// first ncars places are occupied with cars
			unsigned int j=0;
			for(;j<ncars;j++) parkspaces_occupied.insert(std::make_pair(home_loc[tmp[j].second],j));
			// next nparkspaces/2 - ncars are empty parking spaces
			for(;j<nparkspaces/2;j++) parkspaces_empty.insert(std::make_pair(home_loc[tmp[j].second],j));
			// generate the remaining empty parking spaces from work locations
			for(unsigned int i=0;i<nusers;i++) tmp[i] = std::make_pair(rg(),i);
			std::sort(tmp.begin(),tmp.end(),[](auto a, auto b) { return a.first < b.first; });
			for(unsigned int i=0;i+j<nparkspaces;i++) parkspaces_empty.insert(std::make_pair(work_loc[tmp[i].second],j+i));
		}
		else {
			// all home locations are assumed to be occupied parking
			unsigned int j=nusers;
			for(unsigned int i=0;i<nusers;i++) tmp[i] = std::make_pair(rg(),i);
			std::sort(tmp.begin(),tmp.end(),[](auto a, auto b) { return a.first < b.first; });
			for(unsigned int i=0;i+j<nparkspaces;i++) parkspaces_empty.insert(std::make_pair(work_loc[tmp[i].second],j+i));
		}
	}
	
	std::vector<point> end_coords;
	if(!shared) {
		end_coords.resize(nusers);
		for(unsigned int i=0;i<nusers;i++) end_coords[i] = home_loc[i];
	}
	
	std::vector<start_event> events(nusers);
	int r = 0;
	
	for(size_t j=0;j<res.size();j++) {
		res[j].nparkspaces = nparkspaces; // these are always the same in this case
		res[j].ncars = ncars;
		res[j].dist_thres = 0.0;
		res[j].dist_tot = 0.0;
		res[j].trips_thres = 0;
		
		// 1. home to work trips
		for(unsigned int i=0;i<nusers;i++) {
			if(shared) {
				events[i].x = home_loc[i].get<0>();
				events[i].y = home_loc[i].get<1>();
			}
			else {
				events[i].x = end_coords[i].get<0>();
				events[i].y = end_coords[i].get<1>();
				events[i].user_id = i;
				// calculate the distance from the users' home location, add to the extra distances
				double dist = bg::distance(end_coords[i],home_loc[i]);
				res[j].dist_tot += dist;
				if(dist >= dmax) { res[j].dist_thres += dist; res[j].trips_thres++; }
			}
			
			events[i].dest_x = work_loc[i].get<0>();
			events[i].dest_y = work_loc[i].get<1>();
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[i].ts = seq;
				events[i].ttime = 0;
			}
			else {
				unsigned int ts1 = morning_dist(rg);
				events[i].ts = ts1;
				events[i].ttime = travel_times[i].first;
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end());
		if(out && day_out == j) r = process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied,
			end_events, occupy_events, dmax, speed1, out, shared, end_coords, res[j]);
		else r = process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied,
			end_events, occupy_events, dmax, speed1, 0, shared, end_coords, res[j]);
		if(r) return j;
		
		// 2. work to home trips
		for(unsigned int i=0;i<nusers;i++) {
			if(shared) {
				events[i].x = work_loc[i].get<0>();
				events[i].y = work_loc[i].get<1>();
			}
			else {
				events[i].x = end_coords[i].get<0>();
				events[i].y = end_coords[i].get<1>();
				events[i].user_id = i;
				// calculate the distance from the users' work location, add to the extra distances
				double dist = bg::distance(end_coords[i],work_loc[i]);
				res[j].dist_tot += dist;
				if(dist >= dmax) { res[j].dist_thres += dist; res[j].trips_thres++; }
			}

			events[i].dest_x = home_loc[i].get<0>();
			events[i].dest_y = home_loc[i].get<1>();
			
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[i].ts = seq;
				events[i].ttime = 0;
			}
			else {
				unsigned int ts1 = evening_dist(rg);
				events[i].ts = ts1;
				events[i].ttime = travel_times[i].second;
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end()); //, [](const user_event a, const user_event b) { return a.ts < b.ts; });
		if(out && day_out ==j) r = process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied,
			end_events, occupy_events, dmax, speed1, out, shared, end_coords, res[j]);
		else r = process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied,
			end_events, occupy_events, dmax, speed1, 0, shared, end_coords, res[j]);
		if(r) return j;
	}
	return res.size();
}



int main(int argc, char** args)
{
	// input: file with home -- work coordinates, number of times to run, number of threads, repetitions for each day
	char* usersfile = 0;
	char* tdistfile = 0;
	unsigned int morning_length = 2*3600;
	unsigned int evening_length = 2*3600;
	unsigned int morning_start = 6*3600;
	unsigned int morning_end = 9*3600;
	unsigned int afternoon_start = 17*3600;
	unsigned int afternoon_end = 20*3600;
	unsigned int seed = time(0);
	unsigned int days = 30; // number of days to run the simulation to see if it "converges"
	double clon = -71.0584775; // Boston city hall coordinates
	double clat = 42.3605468;
	double rcars = 0.7; // ratio of cars and parking to the number of people to include
	double rpark = 0.7;
	char* outfile = 0;
	int day_out = -1;
	bool shared = true; // use shared cars (if false, use private cars)

	double speed1 = 5.5555555555555; // average travel speed locally: 20 km/h -> 5.5555m m/s
	double dmax = 500; // radius to use when looking for nearby cars / parking spots (in meters)
	bool travel0 = false; // travel times are taken as zero, i.e. estimate a limit on efficiency
	
	for(int i=1;i<argc;i++) if(args[i][0] == '-') switch(args[i][1]) {
		case 'i':
			usersfile = args[i + 1];
			i++;
			break;
		case 's':
			seed = atoi(args[i + 1]);
			i++;
			break;
		case 'd':
			days = atoi(args[i+1]);
			break;
		case 'D':
			tdistfile = args[i+1];
			if(argc > i+2 && args[i+2][0] != '-') {
				morning_start = atoi(args[i+2]);
				morning_end = atoi(args[i+3]);
				afternoon_start = atoi(args[i+4]);
				afternoon_end = atoi(args[i+5]);
			}
			break;
		case 'w':
			morning_length = atoi(args[i+1]);
			if(argc > i+2 && args[i+2][0] != '-') evening_length = atoi(args[i+2]);
			else evening_length = morning_length;
			break;
		case 'c':
			clon = atof(args[i+1]);
			clat = atof(args[i+2]);
			break;
		case 'r':
			dmax = atoi(args[i+1]);
			break;
		case 'S':
			speed1 = atof(args[i+1])*1000.0/3600.0;
			break;
		case '0':
			travel0 = true;
			break;
		case 'o':
			outfile = args[i+1];
			if(i+2 < argc && args[i+1][0] != '-' && args[i+2][0] != '-') day_out = atoi(args[i+2]);
			else day_out = -1;
			break;
		case 'C':
			rcars = atof(args[i+1]);
			break;
		case 'P':
			rpark = atof(args[i+1]);
			break;
		case 'p':
			shared = false;
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s !", args[i]);
			break;
	}
	
	if(outfile && day_out == -1) day_out = days-1;
	
	std::mt19937 rg;
	rg.seed(seed);

	tdist_empirical dist_morning;
	tdist_empirical dist_evening;

	if(tdistfile) {
		if(dist_morning.ReadFreqs(tdistfile,morning_start,morning_end)) {
			fprintf(stderr,"Error reading morning commute time distribution from file %s!\n",tdistfile);
			return 1;
		}
		if(dist_evening.ReadFreqs(tdistfile,afternoon_start,afternoon_end)) {
			fprintf(stderr,"Error reading afternoon commute time distribution from file %s!\n",tdistfile);
			return 1;
		}
	}
	
	std::vector<point> home_loc;
	std::vector<point> work_loc;
	std::vector<std::pair<unsigned int, unsigned int> > travel_times;
	
	FILE* inf = stdin;
	if(usersfile) {
		inf = fopen(usersfile,"r");
		if(inf == 0) {
			fprintf(stderr,"Error opening input file %s!\n",usersfile);
			return 1;
		}
	}
	
	unsigned int line = 0;
	coordconverter conv(clon,clat);
	while(1) {
		int a;
		double hlon,hlat,wlon,wlat;
		unsigned int hwtime,whtime;
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == EOF) break;
		line++;
		if(a == '\n') continue;
		
		ungetc(a,inf); a = fscanf(inf,"%lf",&hlon);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == '\n' || a == EOF) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		ungetc(a,inf); a = fscanf(inf,"%lf",&hlat);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == '\n' || a == EOF) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		ungetc(a,inf); a = fscanf(inf,"%lf",&wlon);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == '\n' || a == EOF) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		ungetc(a,inf); a = fscanf(inf,"%lf",&wlat);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == '\n' || a == EOF) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		ungetc(a,inf); a = fscanf(inf,"%u",&hwtime);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		do a = getc(inf); while(a == ' ' || a == '\t');
		if(a == '\n' || a == EOF) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		ungetc(a,inf); a = fscanf(inf,"%u",&whtime);
		if(a != 1) { fprintf(stderr,"Invalid data on input line %u!\n",line); return 1; }
		
		home_loc.push_back(conv(hlon,hlat));
		work_loc.push_back(conv(wlon,wlat));
		travel_times.push_back(std::make_pair(hwtime,whtime));
		
		do a = getc(inf); while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	if(inf != stdin) fclose(inf);
	if(home_loc.size() == 0) {
		fprintf(stderr,"Error: no data read from input!\n");
		return 1;
	}
	
	FILE* out = 0;
	if(outfile) {
		out = fopen(outfile,"w");
		if(out == 0) fprintf(stderr,"Error opening output file %s!\n",outfile);
	}
	std::vector<res_struct> res(days);
	if(tdistfile) do_estimate(home_loc, work_loc, travel_times, rg, dist_morning, dist_evening, dmax, shared, rcars, rpark, speed1, travel0,
		res, out, day_out);
	else do_estimate(home_loc, work_loc, travel_times, rg, std::uniform_int_distribution<unsigned int>(0,morning_length),
		std::uniform_int_distribution<unsigned int>(0,evening_length), dmax, shared, rcars, rpark, speed1, travel0, res, out, day_out);
	if(out) fclose(out);
	
	for(unsigned int i=0;i<days;i++) fprintf(stdout,"%u\t%u\t%u\t%g\t%g\t%u\n",i,res[i].ncars,res[i].nparkspaces,res[i].dist_tot,res[i].dist_thres,res[i].trips_thres);
	
	return 0;
}

