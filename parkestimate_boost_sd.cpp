/*
 * parkestimate_boost_sd.cpp -- do parking space estimate calculations using Boost's spatial index
 * 	addition: include extra time used for reaching the parking space from the final destination
 *  and the start from a parking space
 * 	use heap structures to keep track of these
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


struct start_event {
	int x;
	int y;
	unsigned int ts;
	int dest_x;
	int dest_y;
	unsigned int ttime;
	bool operator < (const start_event& x) const { return ts < x.ts; }
};

struct end_event {
	int x;
	int y;
	unsigned int ts;
};

struct end_comparer { // comparer for priority_queue 
	bool operator() (end_event& x, end_event& y) const { return x.ts > y.ts; }
};

/*
struct user_event {
	int x;
	int y;
	unsigned int ts;
	bool leave;
	int user_id;
	bool operator < (const user_event& x) const {
		if(ts < x.ts) return true;
		if(ts > x.ts) return false;
		if(user_id < x.user_id) return true;
		if(user_id > x.user_id) return false;
		if(leave == true && x.leave == false) return true;
		return false;
	}
};*/

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


// struct to store main results
struct res_struct {
	unsigned int ncars; // number of cars 
	unsigned int nparkspaces; // number of parking spots
	double dist_tot; // total 'extra' distance traveled (i.e. between the start / destination and parking
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
 * 		res -- (in/out) return the result here
 * 		dmax -- (in) maximum distance that parking should be from start / destination of journeys
 * 		speed1 -- (in) driving speed (using as the crow flies distance) to calculate the time needed to cover the distance between the parking space
 * 			and trip start / end
 * 		grace_period -- (in) if a newly 'created' parking space is occupied less than this time (in seconds) it is discarded, i.e. assumed that
 * 			it was unneeded
 */
template <class it, class se, class index1>
void process_events(it seq, se end, index1& parkspaces_empty, index1& parkspaces_occupied,
		std::priority_queue<end_event, std::vector<end_event>, end_comparer>& end_events,
		std::priority_queue<end_event, std::vector<end_event>, end_comparer>& occupy_events,
		res_struct& res, double dmax, double speed1, unsigned int grace_period) {
	
	if( ! (end_events.empty() && occupy_events.empty()) ) throw new std::runtime_error("process_events(): invalid (non-empty) queue objects provided!\n");
	
	std::vector<value> results;
	double dist_tot = 0.0;
	unsigned int ncars = res.ncars;
	unsigned int nparkspaces = res.nparkspaces;
	
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
			// search for "free" cars around the users's location
			parkspaces_occupied.query(bgi::nearest(p,1),std::back_inserter(results));
			bool found = false;
			double dist = 0.0;
			if(results.size() > 0) {
				dist = bg::distance(results[0].first,p);
				if(dist < dmax) found = true;
				else results.clear();
			}
			
			if(found) {
				// remove and add to as empty parking space
				unsigned int tsadd = (unsigned int)round(dist / speed1);
				unsigned int ts_park = results[0].second; // time this parking space was created
				if(grace_period > 0 && ts_park + grace_period > s.ts) {
					nparkspaces--; // if it was only required very recently, just discard it
								// (assuming the car could just spend that time on the road)
					unsigned int tsdiff = s.ts - ts_park;
					if(tsdiff >= tsadd) tsadd = 0;
					else tsadd -= tsdiff; // decrease the travel time for the time, as it could start earlier already
				}
				else parkspaces_empty.insert(results[0]); // otherwise add it to the available empty parking to be reused later
				if(parkspaces_occupied.remove(results[0]) != 1) {
					throw new std::runtime_error("process_events(): error with remove!\n");
				}
				// add the extra travel time needed to move from the parking space to the user's location
				end_ts += tsadd;
				dist_tot += dist; // keep track of extra distance traveled
			}
			else {
				// add a "new" empty parking space
				results.push_back(std::make_pair(p,nparkspaces));
				nparkspaces++;
				ncars++;
				parkspaces_empty.insert(results[0]);
			}
			
			// add the end of this trip to the end events to be processed
			end_event e;
			e.x = s.dest_x;
			e.y = s.dest_y;
			e.ts = end_ts;
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
			double dist = 0.0;
			if(results.size() > 0) {
				dist = bg::distance(results[0].first,p);
				if(dist < dmax) found = true;
				else results.clear();
			}
			
			if(found) {
				// add to the queue of reserved parking spaces, to be made available later
				end_event o;
				o.x = results[0].first.get<0>();
				o.y = results[0].first.get<1>();
				o.ts = e.ts + (unsigned int)round(dist / speed1);
				occupy_events.push(o);
				// remove the empty parking spot and add as occupied
				if(parkspaces_empty.remove(results[0]) != 1) {
					throw new std::runtime_error("process_events(): error with remove!\n");
				}
				dist_tot += dist; // add to the total extra distance traveled
			}
			else {
				// add a new occupied parking space
				results.push_back(std::make_pair(p,e.ts)); // add with the current timestamp, so if the car is immediately needed,
						// the parking space can be discarded
				nparkspaces++;
				parkspaces_occupied.insert(results[0]);
			}
			
			end_events.pop();
			results.clear();
			continue;
		}
		
		if(onext) {
			// process event for occupying a reserved parking space (add to the list of available cars
			const end_event& o = occupy_events.top();
			parkspaces_occupied.insert(std::make_pair(point(o.x,o.y),0)); // add with zero timestamp so that this parking space cannot be discarded
				// it was there already
			occupy_events.pop();
			continue;
		}
		
		break; // no events left to process
	}
	
	if(seq != end || end_events.size() > 0 || occupy_events.size() > 0) throw new std::runtime_error("process_events(): not all events processed!\n");
	res.dist_tot += dist_tot;
	res.nparkspaces = nparkspaces;
	res.ncars = ncars;
}


template <class tdist>
void do_estimate(std::vector<point>& home_loc, std::vector<point>& work_loc,
			std::vector<std::pair<unsigned int, unsigned int> >& travel_times,
			std::mt19937& rg, tdist&& morning_dist, tdist&& evening_dist, double dmax, double speed1, bool travel0, unsigned int grace_period,
			std::vector<res_struct>& res) {
	if(home_loc.size() != work_loc.size() || work_loc.size() != travel_times.size() || home_loc.size() == 0)
		throw new std::runtime_error("do_estimate(): invalid input!\n");
	
	unsigned int nparkspaces = 0;
	unsigned int ncars = 0;
	unsigned int nusers = home_loc.size();
	
	bgi::rtree< value, bgi::rstar<16> > parkspaces_empty;
	bgi::rtree< value, bgi::rstar<16> > parkspaces_occupied;
	std::priority_queue<end_event, std::vector<end_event>, end_comparer> end_events;
	std::priority_queue<end_event, std::vector<end_event>, end_comparer> occupy_events;
	
	std::vector<start_event> events(nusers);
	
	
	for(size_t j=0;j<res.size();j++) {
		res[j].nparkspaces = nparkspaces;
		res[j].ncars = ncars;
		res[j].dist_tot = 0.0;
		
		// 1. home to work trips
		unsigned int tsmax = 0; // minimum timestamp to put in events -- main purpose: it has to be greater than the maximum ts for creating any new
		// parking spot + grace_period in the last run, so parking spaces created in the last run will not get discarded in the current
		for(auto it = parkspaces_occupied.qbegin(bgi::satisfies([](value const&){ return true; })); it != parkspaces_occupied.qend(); ++it) {
			unsigned int ts0 = it->second;
			if(ts0 > tsmax) tsmax = ts0;
		}
		
		for(unsigned int i=0;i<nusers;i++) {
			events[i].x = home_loc[i].get<0>();
			events[i].y = home_loc[i].get<1>();
			
			events[i].dest_x = work_loc[i].get<0>();
			events[i].dest_y = work_loc[i].get<1>();
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[i].ts = seq;
				events[i].ttime = 0;
			}
			else {
				unsigned int ts1 = tsmax + grace_period + morning_dist(rg);
				events[i].ts = ts1;
				events[i].ttime = travel_times[i].first;
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end()); //, [](const user_event a, const user_event b) { return a.ts < b.ts; });
		process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, end_events, occupy_events,
			res[j], dmax, speed1, grace_period);
		
		// 2. work to home trips
		tsmax = 0;
		for(auto it = parkspaces_occupied.qbegin(bgi::satisfies([](value const&){ return true; })); it != parkspaces_occupied.qend(); ++it) {
			unsigned int ts0 = it->second;
			if(ts0 > tsmax) tsmax = ts0;
		}
		
		for(unsigned int i=0;i<nusers;i++) {
			events[i].x = work_loc[i].get<0>();
			events[i].y = work_loc[i].get<1>();

			events[i].dest_x = home_loc[i].get<0>();
			events[i].dest_y = home_loc[i].get<1>();
			
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[i].ts = seq;
				events[i].ttime = 0;
			}
			else {
				unsigned int ts1 = evening_dist(rg) + grace_period + tsmax;
				events[i].ts = ts1;
				events[i].ttime = travel_times[i].second;
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end()); //, [](const user_event a, const user_event b) { return a.ts < b.ts; });
		process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, end_events, occupy_events,
			res[j], dmax, speed1, grace_period);
		
		ncars = res[j].ncars;
		nparkspaces = res[j].nparkspaces;
	}
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
	unsigned int days = 30;
	double clon = -71.0584775; // Boston city hall coordinates
	double clat = 42.3605468;

	double speed1 = 5.5555555555555; // average travel speed locally: 20 km/h -> 5.5555m m/s
	double dmax = 500; // radius to use when looking for nearby cars / parking spots (in meters)
	bool travel0 = false; // travel times are taken as zero, i.e. estimate a limit on efficiency
	unsigned int grace_period = 0; // if a newly 'created' parking space is occupied in less than this time (in seconds)
			// it is discarded, i.e. assumed that it was unneeded

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
		case 'g':
			grace_period = atoi(args[i+1]);
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s !", args[i]);
			break;
	}

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
	
	std::vector<res_struct> res(days);
	if(tdistfile) do_estimate(home_loc, work_loc, travel_times, rg, dist_morning, dist_evening, dmax, speed1, travel0, grace_period, res);
	else do_estimate(home_loc, work_loc, travel_times, rg, std::uniform_int_distribution<unsigned int>(0,morning_length),
		std::uniform_int_distribution<unsigned int>(0,evening_length), dmax, speed1, travel0, grace_period, res);
	
	for(unsigned int i=0;i<days;i++) fprintf(stdout,"%u\t%u\t%u\t%g\n",i,res[i].ncars,res[i].nparkspaces,res[i].dist_tot);
	
	return 0;
}

