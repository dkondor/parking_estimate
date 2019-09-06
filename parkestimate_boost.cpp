/*
 * parkestimate_boost.cpp -- do parking space estimate calculations using Boost's spatial index
 * 
 * Copyright 2017,2018 Kondor Dániel <dkondor@mit.edu>
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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <random>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef std::pair<point, unsigned int> value;


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
};

class coordconverter {
	protected:
		double clon;
		double clat;
		double factor;
		coordconverter():do_nothing(false) { }
	public:
		bool do_nothing;
		coordconverter(double clon_, double clat_) {
			clon = clon_;
			clat = clat_;
			factor = cos(M_PI*clat/180.0);
			do_nothing = false;
		}
		point operator () (double lon,double lat) {
			if(do_nothing) return point((int)lon,(int)lat);
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

template <class it, class se, class index1>
void process_events(it seq, se end, index1& parkspaces_empty, index1& parkspaces_occupied, res_struct& res,
		FILE* out, FILE* dists_out, bool shared, std::vector<point>& end_coords, double dmax, bool use_end = false) {
	
	std::vector<value> results;
	double dist_tot = 0.0;
	unsigned int ncars = res.ncars;
	unsigned int nparkspaces = res.nparkspaces;
	
	for(;seq != end;seq++) {
		user_event e = *seq;
		point p(e.x,e.y);
		if(e.leave) {
			double dist = 0.0;
			if(shared) {
				// search for "free" cars around the users's location
				parkspaces_occupied.query(bgi::nearest(p,1),std::back_inserter(results));
				bool found = false;
				if(results.size() > 0) {
					dist = bg::distance(results[0].first,p);
					if(dist < dmax) { found = true; dist_tot += dist; }
					else results.clear();
				}
				
				if(found) {
					// remove and add to as empty parking space
					parkspaces_empty.insert(results[0]);
					if(parkspaces_occupied.remove(results[0]) != 1) {
						throw new std::runtime_error("process_events(): error with remove!\n");
					}
				}
				else {
					// add a "new" empty parking space
					dist = 0.0;
					results.push_back(std::make_pair(p,nparkspaces));
					nparkspaces++;
					ncars++;
					parkspaces_empty.insert(results[0]);
				}
			}
			else {
				if(use_end) { 
					results.push_back(std::make_pair(end_coords[e.user_id],e.user_id));
					dist = bg::distance(end_coords[e.user_id],p);
				}
				else results.push_back(std::make_pair(p,e.user_id));
				parkspaces_empty.insert(results[0]);
			}
			if(out) {
				fprintf(out,"%u\t%d\t%d\t",e.ts,e.x,e.y);
				if(e.leave) fprintf(out,"True");
				else fprintf(out,"False");
				fprintf(out,"\t%d\t%d\t%d\t%u\n",e.user_id,results[0].first.get<0>(),results[0].first.get<1>(),results[0].second);
			}
			if(dists_out) fprintf(dists_out,"%f\n",dist);
		}
		else {
			// search for free parking spaces around the user's location
			parkspaces_empty.query(bgi::nearest(p,1),std::back_inserter(results));
			bool found = false;
			double dist = 0.0;
			if(results.size() > 0) {
				dist = bg::distance(results[0].first,p);
				if(dist < dmax) { found = true; dist_tot += dist; }
				else results.clear();
			}
			
			if(found) {
				// remove the empty parking spot and add as occupied
				if(shared) parkspaces_occupied.insert(results[0]);
				else end_coords[e.user_id] = results[0].first;
				if(parkspaces_empty.remove(results[0]) != 1) {
					throw new std::runtime_error("process_events(): error with remove!\n");
				}
			}
			else {
				// add a new occupied parking space
				dist = 0.0;
				results.push_back(std::make_pair(p,nparkspaces));
				nparkspaces++;
				if(shared) parkspaces_occupied.insert(results[0]);
				else end_coords[e.user_id] = results[0].first;
			}
			if(out) {
				fprintf(out,"%u\t%d\t%d\t",e.ts,e.x,e.y);
				if(e.leave) fprintf(out,"True");
				else fprintf(out,"False");
				fprintf(out,"\t%d\t%d\t%d\t%u\n",e.user_id,results[0].first.get<0>(),results[0].first.get<1>(),results[0].second);
			}
			if(dists_out) fprintf(dists_out,"%f\n",dist);
		}
		
		results.clear();
	}
	
	res.dist_tot += dist_tot;
	res.nparkspaces = nparkspaces;
	res.ncars = ncars;
}


template <class tdist>
void do_estimate(std::vector<point>& home_loc, std::vector<point>& work_loc,
			std::vector<std::pair<unsigned int, unsigned int> >& travel_times,
			std::mt19937& rg, tdist&& morning_dist, tdist&& evening_dist, double dmax, bool shared, bool travel0,
			std::vector<res_struct>& res, FILE* out, unsigned int day_out, FILE* trips_out, char* dists_out_base) {
	if(home_loc.size() != work_loc.size() || work_loc.size() != travel_times.size() || home_loc.size() == 0)
		throw new std::runtime_error("do_estimate(): invalid input!\n");
	
	char* fntmp = 0;
	if(dists_out_base) {
		fntmp = (char*)malloc(sizeof(char)*(strlen(dists_out_base) + 20));
		if(!fntmp) throw std::runtime_error("do_estimate(): error allocating memory!\n");
	}
	
	unsigned int nparkspaces = 0;
	unsigned int ncars = 0;
	unsigned int nusers = home_loc.size();
	
	std::vector<point> end_coords;
	bgi::rtree< value, bgi::rstar<16> > parkspaces_empty;
	bgi::rtree< value, bgi::rstar<16> > parkspaces_occupied;
	
	if(!shared) {
		end_coords.resize(nusers);
		for(unsigned int i=0;i<nusers;i++) end_coords[i] = home_loc[i];
		ncars = nusers;
		nparkspaces = nusers;
	}
	
	std::vector<user_event> events(2*nusers);
	
	for(size_t j=0;j<res.size();j++) {
		res[j].nparkspaces = nparkspaces;
		res[j].ncars = ncars;
		res[j].dist_tot = 0.0;
		
		FILE* dists_out = 0;
		if(dists_out_base) {
			sprintf(fntmp,"%s_d%u.dat",dists_out_base,(unsigned int)j);
			dists_out = fopen(fntmp,"w");
			if(!dists_out) {
				fprintf(stderr,"do_estimate(): error opening file %s!\n",fntmp);
				throw std::runtime_error("do_estimate(): error opening output file!\n");
			}
		}
		
		// 1. home to work trips
		unsigned int tsmax = 0;
		for(unsigned int i=0;i<nusers;i++) {
			if(shared) {
				events[2*i].x = home_loc[i].get<0>();
				events[2*i].y = home_loc[i].get<1>();
			}
			else {
				events[2*i].x = end_coords[i].get<0>();
				events[2*i].y = end_coords[i].get<1>();
			}
			events[2*i].leave = true;
			events[2*i].user_id = i;
			events[2*i+1].x = work_loc[i].get<0>();
			events[2*i+1].y = work_loc[i].get<1>();
			events[2*i+1].leave = false;
			events[2*i+1].user_id = i;
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[2*i].ts = seq;
				events[2*i+1].ts = seq;
			}
			else {
				unsigned int ts1 = morning_dist(rg);
				events[2*i].ts = ts1;
				events[2*i+1].ts = ts1 + travel_times[i].first;
				if(events[2*i+1].ts > tsmax) tsmax = events[2*i+1].ts;
			}
			if(trips_out) {
				// output whole trip
				fprintf(trips_out,"%d\t%u\t%d\t%d\t%u\t%d\t%d\n",i,events[2*i].ts,home_loc[i].get<0>(),home_loc[i].get<1>(),
					events[2*i+1].ts,work_loc[i].get<0>(),work_loc[i].get<1>());
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end()); //, [](const user_event a, const user_event b) { return a.ts < b.ts; });
		if(out && j == day_out)
			process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, res[j], out, dists_out, shared, end_coords, dmax);
		else process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, res[j], 0, dists_out, shared, end_coords, dmax);
		
		// 2. work to home trips
		for(unsigned int i=0;i<nusers;i++) {
			if(shared) {
				events[2*i].x = work_loc[i].get<0>();
				events[2*i].y = work_loc[i].get<1>();
			}
			else {
				events[2*i].x = end_coords[i].get<0>();
				events[2*i].y = end_coords[i].get<1>();
			}
			events[2*i].leave = true;
			events[2*i].user_id = i;
			events[2*i+1].x = home_loc[i].get<0>();
			events[2*i+1].y = home_loc[i].get<1>();
			events[2*i+1].leave = false;
			events[2*i+1].user_id = i;
			if(travel0) {
				unsigned int seq = rg(); // random "sequence number" for this person
				events[2*i].ts = seq;
				events[2*i+1].ts = seq;
			}
			else {
				unsigned int ts1 = tsmax + evening_dist(rg);
				events[2*i].ts = ts1;
				events[2*i+1].ts = ts1 + travel_times[i].second;
			}
			if(trips_out) {
				// output whole trip
				fprintf(trips_out,"%d\t%u\t%d\t%d\t%u\t%d\t%d\n",i,events[2*i].ts,work_loc[i].get<0>(),work_loc[i].get<1>(),
					events[2*i+1].ts,home_loc[i].get<0>(),home_loc[i].get<1>());
			}
		}
		// sort by timestamp, do the processing
		std::sort(events.begin(), events.end()); //, [](const user_event a, const user_event b) { return a.ts < b.ts; });
		if(out && j == day_out)
			process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, res[j], out, dists_out, shared, end_coords, dmax);
		else process_events(events.begin(), events.end(), parkspaces_empty, parkspaces_occupied, res[j], 0, dists_out, shared, end_coords, dmax);
		
		ncars = res[j].ncars;
		nparkspaces = res[j].nparkspaces;
		if(dists_out) fclose(dists_out);
	}
}


// read in trips from the given input
// format should be: user_id, start_ts, start_lon, start_lat, end_ts, end_lon, end_lat
int read_user_events(std::vector<user_event>& events, FILE* f, coordconverter& c) {
	unsigned int line = 0;
	while(1) {
		int a;
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == EOF) break;
		line++;
		
		user_event start,end;
		double lon,lat;
		point p;
		
		ungetc(a,f);
		a = fscanf(f,"%d",&(start.user_id));
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%u",&(start.ts));
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%lg",&lon);
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%lg",&lat);
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		p = c(lon,lat);
		start.x = p.get<0>();
		start.y = p.get<1>();
		
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%u",&(end.ts));
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%lg",&lon);
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		do a = getc(f); while(a == '\t' || a == ' ');
		if(a == '\n' || a == '\r' || a == EOF) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		ungetc(a,f);
		a = fscanf(f,"%lg",&lat);
		if(a != 1) { fprintf(stderr,"read_user_events(): invalid data on input line %u!\n",line); return 1; }
		p = c(lon,lat);
		end.x = p.get<0>();
		end.y = p.get<1>();
		
		end.user_id = start.user_id;
		start.leave = true;
		end.leave = false;
		
		events.push_back(start);
		events.push_back(end);
		
		do a = getc(f); while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	return 0;
}


int main(int argc, char** args)
{
	// input: file with home -- work coordinates, number of times to run, number of threads, repetitions for each day
	char* usersfile = 0;
	char* tdistfile = 0;
	char* outfile = 0;
	char* dists_out_base = 0; /* base filename to output distance distribution (from trip start / end to the vehicles) */
	int day_out = -1;
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
	bool coord_no_convert = false;

	double dmax = 500; // radius to use when looking for nearby cars / parking spots (in meters)
	bool shared = true; // use car sharing (i.e. everyone takes any available car) vs. use private cars and share only parking spots
	bool travel0 = false; // travel times are taken as zero, i.e. estimate a limit on efficiency

	/*
	 * special usage: just consider one day, with the possibility to write out the sequence of trips to a file or read them from a file
	 */
	bool oneday = false;
	char* trips_outfile = 0;
	char* trips_infile = 0;

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
			i+=2;
			break;
		case 'C':
			coord_no_convert = true;
			break;
		case 'r':
			dmax = atoi(args[i+1]);
			break;
		case 'p':
			shared = false;
			break;
		case '0':
			travel0 = true;
			break;
		case 'o':
			outfile = args[i+1];
			if(i+2 < argc && args[i+1][0] != '-' && args[i+2][0] != '-') day_out = atoi(args[i+2]);
			else day_out = -1;
			break;
		case 'I':
			trips_infile = args[i+1];
			i++;
			oneday = true;
			break;
		case 'O':
			trips_outfile = args[i+1];
			i++;
			oneday = true;
			break;
		case 'b':
			dists_out_base = args[i+1];
			i++;
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
	
	coordconverter conv(clon,clat);
	conv.do_nothing = coord_no_convert;
	
	if(oneday && trips_infile) {
		// just read trips from the input, run one iteration of the simulation
		std::vector<user_event> events;
		FILE* tinf = fopen(trips_infile,"r");
		if(!tinf) {
			fprintf(stderr,"Error opening input file %s!\n",trips_infile);
			return 1;
		}
		int a = read_user_events(events,tinf,conv);
		fclose(tinf);
		if(a != 0) return 1;
		
		std::sort(events.begin(),events.end());
		
		std::vector<point> end_coords;
		bgi::rtree< value, bgi::rstar<16> > parkspaces_empty;
		bgi::rtree< value, bgi::rstar<16> > parkspaces_occupied;
		res_struct res1 = {0,0,0.0};
		
		
		if(!shared) {
			// for each user, fill out their "last known" location as the start of their first trip
			//		(i.e. assume they have their car parked there)
			// note: this part is not optimal in terms of runtime
			int max_user_id = 0;
			for(auto it = events.begin(); it != events.end(); ++it) if(it->user_id > max_user_id) max_user_id = it->user_id;
			end_coords.resize(max_user_id+1);
			for(auto it = events.rbegin(); it != events.rend(); ++it) if(it->leave == true)
				end_coords[it->user_id] = point(it->x,it->y);
		}
		
		process_events(events.begin(),events.end(),parkspaces_empty,parkspaces_occupied,res1,0,0,shared,end_coords,dmax,true);
		fprintf(stdout,"%u\t%u\t%u\t%g\n",0,res1.ncars,res1.nparkspaces,res1.dist_tot);
		return 0;
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
	
	FILE* trips_out = 0;
	if(oneday && trips_outfile) {
		days = 1;
		trips_out = fopen(trips_outfile,"w");
		if(!trips_out) {
			fprintf(stderr,"Error opening output file %s!\n",trips_outfile);
			return 1;
		}
	}
	FILE* out = 0;
	if(outfile) {
		out = fopen(outfile,"w");
		if(out == 0) fprintf(stderr,"Error opening output file %s!\n",outfile);
	}
	std::vector<res_struct> res(days);
	if(tdistfile) do_estimate(home_loc, work_loc, travel_times, rg, dist_morning, dist_evening, dmax, shared, false, res, out, day_out, trips_out, dists_out_base);
	else do_estimate(home_loc, work_loc, travel_times, rg, std::uniform_int_distribution<unsigned int>(0,morning_length),
		std::uniform_int_distribution<unsigned int>(0,evening_length), dmax, shared, travel0, res, out, day_out, trips_out, dists_out_base);
	if(out) fclose(out);
	if(trips_out) fclose(trips_out);
	
	for(unsigned int i=0;i<days;i++) fprintf(stdout,"%u\t%u\t%u\t%g\n",i,res[i].ncars,res[i].nparkspaces,res[i].dist_tot);
	
	return 0;
}

