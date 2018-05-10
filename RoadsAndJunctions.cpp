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
#include <fstream>
using namespace std;
#define rep(i,n) for (int (i)=(0);(i)<(int)(n);++(i))

#define LOCAL

// clock
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

#define LIMIT_SEC 9.5

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
    static const int MAX = 1000001;
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

const int MAXE = 100010 * 2;
struct edge { int u, v; double cost; edge(int u, int v, double cost) :u(u), v(v), cost(cost) {} edge() {} };

bool comp(const edge& e1, const edge& e2) {
    return e1.cost < e2.cost;
}

ofstream fout("text.txt");

class RoadsAndJunctions {
public:

    int NC;
    int S;
    vector<int> c;
    UnionFind uf;

    //int E, V;
    edge es[MAXE];
    vector<int> junctions;

    Timer timer;

    vector<int> js;

    vector<int> buildJunctions(int s, vector<int> cities, double junctionCost, double failureProbability) {
        // store number of cities for building the roads

        NC = cities.size() / 2;
        //V = NC;
        S = s;
        c = cities;

        // for (int i=0; i<NC; ++i) {
        //     for (int j=i+1; j<NC; ++j) {
        //         double d = dist(c[i*2], c[i*2+1], c[j*2], c[j*2+1]);
        //         es[E++] = edge(i, j, d);
        //         es[E++] = edge(j, i, d);
        //     }
        // }
        //
        vector<int> ccc = c;
        double score = make_greedy_score(ccc, ccc.size()/2);

        set<pair<int, int>> st;
        for (int i=0; i<c.size(); ++i) {
            st.insert(make_pair(c[i*2], c[i*2+1]));
        }

        int idx = 0;
        while (!timer.time_over()) {

            int x, y;
            if (timer.get_time() < timer.get_limit()/4) {
                 x = rnd()%(S/2), y = rnd()%(S/2);
            }
            else if (timer.get_time() < timer.get_limit()/4*2) {
                x = rnd()%S; y = rnd()%(S/2);
            }
            else if (timer.get_time() < timer.get_limit()/4*3) {
                x = rnd()%(S/2), y = rnd()%S;
            }
            else {
                x = rnd()%S, y = rnd()%S;
            }

            // int x = cities[idx*2], y = cities[idx*2+1];
            // bool tmp = rnd()%2;
            // x += (rnd()%100 + 10) * (tmp ? 1 : -1);
            // if (x < 0) x = 0; else if (x >= S) x = S-1;
            // tmp = rnd()%2;
            // y += (rnd()%100 + 10) * (tmp ? 1 : -1);
            // if (y < 0) y = 0; else if (y >= S) y = S-1;
            // idx = (idx + 1) % NC;

            if (st.count(make_pair(x, y))) continue;
            ccc.push_back(x);
            ccc.push_back(y);
            double tmp_score = make_greedy_score(ccc, ccc.size()/2);
            if (tmp_score+(ccc.size()/2)*junctionCost < score) {
                score = tmp_score+(ccc.size()/2)*junctionCost;
                junctions.push_back(x);
                junctions.push_back(y);
                st.insert(make_pair(x, y));
            }
            else {
                ccc.pop_back();
                ccc.pop_back();
            }
        }

        return junctions;
    }

    int dame = 0;

    vector<int> buildRoads(vector<int> junctionStatus) {
        cerr << junctionStatus.size() << endl;
        for (int i = 0; i < junctionStatus.size(); ++i) {
            c.push_back(junctions[i*2]);
            c.push_back(junctions[i*2+1]);
            if (junctionStatus[i] == 0) dame++;
        }
        js = junctionStatus;
        return make_greedy(c, c.size()/2).first;
    }

    double make_greedy_score(vector<int>& c, int V) {
        double score = 0;
        uf.init(V);
        while (uf.size(0) < V) {
            double d = 1e8;
            int idx1 = -1, idx2 = -1;
            for (int u=0; u<V; ++u) {
                for (int v=0; v<V; ++v) {
                    if (uf.same(u, v)) continue;
                    double tmp = dist(c[u*2], c[u*2+1], c[v*2], c[v*2+1]);
                    if (d > tmp) {
                        d = tmp;
                        idx1 = u;
                        idx2 = v;
                    }
                }
            }
            if (d != 1e8) {
                uf.unite(idx1, idx2);
                score += dist(c[idx1*2], c[idx1*2+1], c[idx2*2], c[idx2*2+1]);
            }
        }
        return score;
    }

    pair<vector<int>, double> make_greedy(vector<int>& c, int V) {
        double score = 0;
        uf.init(V);
        vector<int> ret;
        while (uf.size(0) < V-dame) {
            double d = 1e8;
            int idx1 = -1, idx2 = -1;
            for (int u=0; u<V; ++u) {
                if (u >= NC and js[u-NC] == 0) continue;
                for (int v=0; v<V; ++v) {
                    if (v >= NC and js[v-NC] == 0) continue;
                    if (uf.same(u, v)) continue;
                    double tmp = dist(c[u*2], c[u*2+1], c[v*2], c[v*2+1]);
                    if (d > tmp) {
                        d = tmp;
                        idx1 = u;
                        idx2 = v;
                    }
                }
            }
            if (d != 1e8) {
                uf.unite(idx1, idx2);
                score += dist(c[idx1*2], c[idx1*2+1], c[idx2*2], c[idx2*2+1]);
                ret.push_back(idx1);
                ret.push_back(idx2);
            }
        }
        return {ret, score};
    }

    // pair<vector<int>, double> make_kruskal(int V) {
    //     sort(es, es+V, comp);
    //     double score = 0;
    //     vector<int> ret;
    //     uf.init(V);
    //     for (int i=0; i<E; ++i) {
    //         edge e = es[i];
    //         if (!uf.same(e.u, e.v)) {
    //             uf.unite(e.u, e.v);
    //             ret.push_back(e.u);
    //             ret.push_back(e.v);
    //             score += e.cost;
    //         }
    //     }
    //     return {ret, score};
    // }

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
