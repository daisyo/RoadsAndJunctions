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
#include <bitset>
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
    static const int MAX = 500;
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

// const int MAXE = 100010 * 2;
// struct edge { int u, v; double cost; edge(int u, int v, double cost) :u(u), v(v), cost(cost) {} edge() {} };
//
// bool comp(const edge& e1, const edge& e2) {
//     return e1.cost < e2.cost;
// }
//
// ofstream fout("text.txt");

bool used[1000*1000+2];

int best_state[4][500];
int best_sz[4] = {0};
pair<int, int> pos[4][500];
int pos_sz[4] = {0};

int ccc[500];
int ccc_sz = 0;

int junctions_state[4][500];
int junctions_state_sz[4] = {0};

class RoadsAndJunctions {
public:

    int NC;
    int S;
    vector<int> c;
    UnionFind uf;

    //int E, V;
    //edge es[MAXE];
    vector<int> junctions;

    Timer timer;

    vector<int> js;

    vector<int> buildJunctions(int s, vector<int> cities, double junctionCost, double failureProbability) {
        // store number of cities for building the roads

        //V = NC;
        NC = cities.size() / 2;
        S = s;
        c = cities;
        for (int i=0; i<NC; ++i) {
            ccc[i*2] = c[i*2];
            ccc[i*2+1] = c[i*2+1];
            ccc_sz += 2;
        }

        double score = make_greedy_score2();
        //vector<int> best;

        for (int i=0; i<NC; ++i) {
            int x = cities[i*2], y = cities[i*2+1], idx = 0;
            used[cities[i*2+1]*S+cities[i*2]] = 1;
            if (x < S/2 and y < S/2)       idx = 0;
            else if (x >= S/2 and y < S/2) idx = 1;
            else if (x < S/2 and y >= S/2) idx = 2;
            else                           idx = 3;
            pos[idx][pos_sz[idx]++] = make_pair(x, y);
        }

        double limit = LIMIT_SEC - 3, lim = limit/4.0;;
        long long cnt=0;
        timer.set_limit(limit);

        int idx = 0;
        double current_time, start_time = timer.get_time();
        double start_temp = 10, end_temp = 0.5;
        double diff_temp = (end_temp - start_temp) / lim;

        while (idx < 4 and pos_sz[idx] == 0) idx++;

        while ((current_time = timer.get_time()) < LIMIT_SEC and idx < 4) {
            ++cnt;
            if (current_time >= (idx+1)*lim) {
                start_time = current_time;
                idx++;
                if (pos_sz[idx] == 0) idx++;
                if (idx >= 4) break;
            }

            bool f = (rnd()%10>=4);
            if (junctions_state_sz[idx]/2 >= (2*NC)/4) f = false;

            if (f) {
                int rid = rnd()%pos_sz[idx];
                int gx = pos[idx][rid].first, gy = pos[idx][rid].second;
                bool tmp1 = rnd()%2, tmp2 = rnd()%2;
                gx += (rnd()%30+5)*(tmp1 ? 1 : -1);
                gy += (rnd()%30+5)*(tmp2 ? 1 : -1);
                if (gx < 0) gx = 0; else if (gx >= S) gx = S-1;
                if (gy < 0) gy = 0; else if (gy >= S) gy = S-1;
                if (used[gy*S+gx]) continue;
                ccc[ccc_sz++] = gx;
                ccc[ccc_sz++] = gy;
                double tmp_score = make_greedy_score2();
                bool force = true;
                double change = score - (tmp_score + (junctions_state_sz[idx]/2+1)*junctionCost);
                if (change < 0) force = exp(change / (start_temp + diff_temp * (current_time - start_time))) * UINT32_MAX > rnd();
                else force = true;

                if (force) {
                    score = tmp_score + (junctions_state_sz[idx]/2+1)*junctionCost;
                    junctions_state[idx][junctions_state_sz[idx]++] = gx;
                    junctions_state[idx][junctions_state_sz[idx]++] = gy;
                    used[gy*S+gx] = 1;
                    if (change > 0) {
                        memcpy(best_state[idx], junctions_state[idx], sizeof(junctions_state[idx]));
                        best_sz[idx] = junctions_state_sz[idx];
                    }
                }
                else {
                    ccc_sz -= 2;
                }
            }
            else {
                int JC = ccc_sz/2 - NC;
                if (JC == 0) continue;

                int rid = rnd()%JC;
                int x = ccc[NC*2+rid*2], y = ccc[NC*2+rid*2+1];
                swap(ccc[ccc_sz-2], ccc[NC*2+rid*2]);
                swap(ccc[ccc_sz-1], ccc[NC*2+rid*2+1]);
                ccc_sz -= 2;
                double tmp_score = make_greedy_score2();
                bool force = false;
                double change = score - (tmp_score + (junctions_state_sz[idx]/2+1)*junctionCost);
                if (change < 0) force = exp(change / (start_temp + diff_temp * (current_time - start_time))) * UINT32_MAX > rnd();
                else force = true;

                if (force) {
                    score = tmp_score + (junctions_state_sz[idx]/2+1)*junctionCost;
                    swap(junctions_state[idx][junctions_state_sz[idx]-2], junctions_state[idx][rid*2]);
                    swap(junctions_state[idx][junctions_state_sz[idx]-1], junctions_state[idx][rid*2+1]);
                    junctions_state_sz[idx] -= 2;
                    if (change > 0) {
                        memcpy(best_state[idx], junctions_state[idx], sizeof(junctions_state[idx]));
                        best_sz[idx] = junctions_state_sz[idx];
                    }
                }
                else {
                    ccc[ccc_sz++] = x;
                    ccc[ccc_sz++] = y;
                }
            }
        }

        for (int idx=0; idx<4; ++idx)
        for (int i=0; i<best_sz[idx]; ++i) {
            junctions.emplace_back(best_state[idx][i]);
        }

        timer.set_limit(LIMIT_SEC);
        int jcsz = junctions.size()/2;
        if (jcsz == 0) return junctions;

        long long cnt2 = 0;
        while (!timer.time_over()) {
            ++cnt2;
            int idx = rnd()%jcsz;
            int prev_x = ccc[2*NC+idx*2], prev_y = ccc[2*NC+idx*2+1];
            int& x = ccc[2*NC+idx*2];
            int& y = ccc[2*NC+idx*2+1];
            bool tmp1 = rnd()%2, tmp2 = rnd()%2;
            if (rnd()%2) {
                x += (rnd()%5)*(tmp1 ? 1 : -1);
            }
            else {
                y += (rnd()%5)*(tmp2 ? 1 : -1);
            }
            if (x < 0) x = 0; else if (x >= S) x = S-1;
            if (y < 0) y = 0; else if (y >= S) y = S-1;
            if (used[y*S+x]) {
                x = prev_x;
                y = prev_y;
            }
            double tmp_score = make_greedy_score2();
            if (tmp_score+(junctions.size()/2)*junctionCost < score) {
                score = tmp_score+(junctions.size()/2)*junctionCost;
                junctions[idx*2] = x;
                junctions[idx*2+1] = y;
                used[prev_y*S+prev_x] = 0;
                used[y*S+x] = 1;
            }
            else {
                x = prev_x;
                y = prev_y;
            }
        }

        return junctions;
    }

    int dame = 0;

    vector<int> buildRoads(vector<int> junctionStatus) {
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

    double make_greedy_score2() {
        double score = 0;
        int V = ccc_sz / 2;
        uf.init(V);
        while (uf.size(0) < V) {
            double d = 1e8;
            int idx1 = -1, idx2 = -1;
            for (int u=0; u<V; ++u) {
                for (int v=0; v<V; ++v) {
                    if (uf.same(u, v)) continue;
                    double tmp = dist(ccc[u*2], ccc[u*2+1], ccc[v*2], ccc[v*2+1]);
                    if (d > tmp) {
                        d = tmp;
                        idx1 = u;
                        idx2 = v;
                    }
                }
            }
            if (d != 1e8) {
                uf.unite(idx1, idx2);
                score += dist(ccc[idx1*2], ccc[idx1*2+1], ccc[idx2*2], ccc[idx2*2+1]);
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
