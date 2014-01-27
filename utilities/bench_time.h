#ifndef BENCH_TIME_H__
#define BENCH_TIME_H__

#include <sys/times.h>
#include <unistd.h>

inline void bench_time(double &utime, double &stime)
{
	tms t;
	times(&t);
	double ticks = sysconf(_SC_CLK_TCK);
	utime = t.tms_utime / ticks;
	stime = t.tms_stime / ticks;
}

#endif