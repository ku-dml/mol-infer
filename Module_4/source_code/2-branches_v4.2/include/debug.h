// This file gives some convenient ways of debugging.

#ifndef _DEBUG_H
#define _DEBUG_H

#ifndef DEBUG
#define DEBUG false
#endif

#define DBG1(x) if(DEBUG) std::cerr << #x << " = " << x << std::endl<< std::flush
#define DBG1d(x, y) if(DEBUG) std::cerr << #x << "[" << y << "] = " << x[y] << std::endl<< std::flush
#define DBG2d(x, y, z) if(DEBUG) std::cerr << #x << "[" << y << "]" << "[" << z << "] = " << x[y][z] << std::endl<< std::flush
#define DBG1D(x)  if (DEBUG){ cerr << #x << endl<< std::flush; for (const auto it: x) cerr << it <<", "; cerr << std::endl<< std::flush; }
#define DBG2D(x) if (DEBUG) {cerr << #x << endl<< std::flush;  for (const auto it: x){ for (const auto  itt: it) cerr << itt << ", "; cerr << std::endl<< std::flush;} cerr << std::endl << std::flush;}
#define MSG(msg) if(DEBUG) std::cerr << msg << std::endl << std::flush

#define ECHO(x) std::cout << #x << " = " << x << std::endl
#define ECHO1d(x, y) std::cout << #x << "[" << y << "] = " << x[y] << std::endl
#define ECHO2d(x, y, z) if(DEBUG) std::cout << #x << "[" << y << "]" << "[" << z << "] = " << x[y][z] << std::endl
#define ECHO1D(x) cout << #x << endl; for (const auto it: x) cout << it <<", "; cout << std::endl

#define DBGCon(x) if (DEBUG){ std::cerr<< #x << std::endl; for (auto &y:x) cerr << y << std::endl; }

namespace std{
	template <class T1, class T2>
	ostream& operator << (ostream& os, const pair<T1, T2> &p){
		os << "(" << p.first << "," << p.second << ")";
		return os;
	}
}
#endif
