// C++11
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <queue>
#include <map>
#include <set>
#include <random>
using namespace std;
#define rep(i,n) for (int (i)=(0);(i)<(int)(n);++(i))

#define LOCAL

// clock
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

#define LIMIT_SEC 9.8

#define LOCAL_SCALE 2494000000ULL
#define TOPCODER_SCALE 2800000000ULL

class Timer {
private:
    double start, end;
    double limit;

    double get_cycle() {
#ifdef LOCAL
        unsigned long long tm;
        tm = __rdtsc();
        return (double)tm / LOCAL_SCALE;
#else
        unsigned long long low, high;
        __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
        return (double)((low) | (high << 32)) / TOPCODER_SCALE;
#endif
    }

public:

    Timer() {
        start = get_cycle();
        limit = LIMIT_SEC;
    }
    Timer(double limit) :limit(limit) {
        start = get_cycle();
    }
    double get_time() {
        end = get_cycle();
        return end - start;
    }
    bool time_over() {
        if (get_time() > limit) {
            return true;
        }
        return false;
    }
    void set_limit(double l) {
        limit = l;
    }
    void set_start() {
        start = get_cycle();
    }
    double get_limit() const {
        return limit;
    }
};

class XorShift {
public:
    unsigned long x, y, z, w;
    XorShift() {
        x = 123456789; y = 362436069; z = 521288629; w = 88675123;
    }
    XorShift(unsigned long seed) {
        XorShift();
        w = seed;
        for (int i = 0; i < 100; ++i) (*this)();
    }
    unsigned long operator()() {
        unsigned long t = x ^ (x << 11);
        x = y; y = z; z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
};

XorShift rnd;

random_device rd;
mt19937 rand_gen(rd());

struct UnionFind {
    static const int MAX = 1001;
    int par[MAX];
    void init(int N) {
        for (int i=0; i<N; ++i) {
            par[i] = -1;
        }
    }
    void init() {
        for (int i=0; i<MAX; ++i) {
            par[i] = -1;
        }
    }
    int find(int x) {
        if (par[x] < 0) return x;
        else return par[x] = find(par[x]);
    }
    bool unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) return false;
        if (par[x] > par[y]) {
            int tmp = x;
            x = y;
            y = tmp;
        }
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    bool same(int x, int y) {
        return find(x) == find(y);
    }
    // そのグループの木の大きさを返す
    int size(int x) {
        return -par[find(x)];
    }
};


class RoadsAndJunctions {
public:

    int NC;
    int S;
    vector<int> c;
    UnionFind uf;

    vector<int> buildJunctions(int S, vector<int> cities, double junctionCost, double failureProbability) {
        // store number of cities for building the roads
        NC = cities.size() / 2;
        S = this->S;
        c = cities;
        uf.init(NC);
        return {};
    }
    vector<int> buildRoads(vector<int> junctionStatus) {
        // build a road from the single junction to each city
        // (we assume that it was built, but don't check it)
        vector<int> ret;

        while (uf.size(0) < NC) {
            for (int u=0; u<NC; ++u) {
                double d = 1e8;
                int idx = -1;
                for (int v=0; v<NC; ++v) {
                    if (uf.same(u, v)) continue;
                    double tmp = dist(c[u*2], c[u*2+1], c[v*2], c[v*2+1]);
                    if (d > tmp) {
                        d = tmp;
                        idx = v;
                    }
                }
                if (idx != -1) {
                    uf.unite(u, idx);
                    ret.push_back(u);
                    ret.push_back(idx);
                }
            }
        }

        return ret;
    }

    double dist(double ax, double ay, double bx, double by) {
        return sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by));
    }
};
// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

#ifdef LOCAL

int main() {
    RoadsAndJunctions rj;
    int S, C;
    cin >> S >> C;
    vector<int> cities(C);
    getVector(cities);
    double junctionCost, failureProbability;
    cin >> junctionCost >> failureProbability;

    vector<int> ret = rj.buildJunctions(S, cities, junctionCost, failureProbability);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();

    int J;
    cin >> J;
    vector<int> junctionStatus(J);
    getVector(junctionStatus);

    ret = rj.buildRoads(junctionStatus);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}
#endif
