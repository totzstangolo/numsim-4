/*
Copyright (C) 2012   Malte Brunn

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <sys/time.h>
#include <vector>
//------------------------------------------------------------------------------
#ifndef _ZEITGEIST_HPP
#define _ZEITGEIST_HPP
//------------------------------------------------------------------------------
class ZeitGeist {
public:
	ZeitGeist () {
		ticks.clear();
	}
	void Start () {
		ticks.clear();
		timeval tv;
		gettimeofday(&tv, NULL);
		ticks.push_back(tv);
	}
	unsigned long Stop () {
		timeval tv;
		gettimeofday(&tv, NULL);
		ticks.push_back(tv);
		return Total ();
	}
	unsigned long Step () {
		timeval tv;
		gettimeofday(&tv, NULL);
		ticks.push_back(tv);
		return Interval(ticks.size()-2,1);
	}
	void Reset () {
		ticks.clear();
	}
	unsigned int Size () { return ticks.size(); }
	unsigned long Interval (unsigned int start, unsigned int steps) {
		if (start >= ticks.size()) return 0;
		if (start + steps >= ticks.size())
			return ticks.back().tv_sec * 1000000 + ticks.back().tv_usec
				- ticks[start].tv_sec * 1000000 - ticks[start].tv_usec;
		return ticks[start+steps].tv_sec * 1000000 + ticks[start+steps].tv_usec
			- ticks[start].tv_sec * 1000000 - ticks[start].tv_usec;
	}
	unsigned long Total () {
		return ticks.back().tv_sec * 1000000 + ticks.back().tv_usec
			- ticks.front().tv_sec * 1000000 - ticks.front().tv_usec;
	}
	unsigned long Last () {
		if (ticks.size() <= 1) return 0;
		return ticks.back().tv_sec * 1000000 + ticks.back().tv_usec
			- ticks[ticks.size()-2].tv_sec * 1000000 - ticks[ticks.size()-2].tv_usec;
	}
private:
	std::vector<timeval> ticks;
};
//------------------------------------------------------------------------------
#endif
