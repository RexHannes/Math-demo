#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

static int cpp_mod_int(const cpp_int &x, int m){
    cpp_int r = x % m;
    int v = (int)r.convert_to<long long>();
    if(v<0) v += m;
    return v;
}

static cpp_int lcm_range_int(int m, int M){
    cpp_int L = 1;
    for(int n=m;n<=M;n++){
        int rem = cpp_mod_int(L, n);
        int g = std::gcd(rem, n);
        L = (L / g) * n;
    }
    return L;
}

static int m_star_ld(int M){
    long double s=0;
    for(int t=M;t>=2;t--){
        s += 1.0L/t;
        if(s >= 1.0L - 1e-18L) return t;
    }
    return 2;
}

struct ModProbeResult{
    bool reachable;
    int target;
    int reach_count;
    vector<pair<int,int>> nonzero;
};

static ModProbeResult probe_modulus(int m, int M, int R, bool keep_nonzero){
    cpp_int L = lcm_range_int(m,M);
    int len = M-m+1;
    vector<int> weights(len);
    int target = cpp_mod_int(-L, R);
    vector<pair<int,int>> nz;
    for(int i=0;i<len;i++){
        int t=m+i;
        cpp_int q = L / t;
        int w = cpp_mod_int(q, R);
        weights[i]=w;
        target += w;
        target %= R;
        if(keep_nonzero && w!=0) nz.push_back({t,w});
    }
    vector<unsigned char> dp0(R,0), dp1(R,0), ndp0(R,0), ndp1(R,0);
    dp0[0]=1;
    for(int i=0;i<len;i++){
        int t=m+i;
        fill(ndp0.begin(), ndp0.end(), 0);
        fill(ndp1.begin(), ndp1.end(), 0);
        int w=weights[i];
        for(int res=0;res<R;res++){
            if(dp0[res]||dp1[res]) ndp0[res]=1;
            if(t!=m && t!=M && dp0[res]) ndp1[(res+w)%R]=1;
        }
        dp0.swap(ndp0); dp1.swap(ndp1);
    }
    int reach_count=0;
    for(int res=0;res<R;res++) if(dp0[res]||dp1[res]) reach_count++;
    bool reachable = dp0[target] || dp1[target];
    return {reachable,target,reach_count,nz};
}

static string nz_to_string(const vector<pair<int,int>> &nz){
    string s;
    for(size_t i=0;i<nz.size();i++){
        if(i) s += ";";
        s += to_string(nz[i].first) + ":" + to_string(nz[i].second);
    }
    return s;
}

int main(int argc, char** argv){
    int Mmin=200, Mmax=450, R=33;
    string prefix="results/modulus_probe";
    bool keep_examples=true;
    for(int i=1;i<argc;i++){
        string a=argv[i];
        auto next=[&](){ return string(argv[++i]); };
        if(a=="--M-min") Mmin=stoi(next());
        else if(a=="--M-max") Mmax=stoi(next());
        else if(a=="--modulus") R=stoi(next());
        else if(a=="--prefix") prefix=next();
        else if(a=="--no-examples") keep_examples=false;
    }
    ofstream summ(prefix + "_summary_by_M.csv");
    summ << "M,mstar,total,blocked,reachable,first_reachable_m,last_reachable_m,example_blocked_nz,example_reachable_nz\n";
    ofstream rec(prefix + "_reachable_records.csv");
    rec << "m,M,target,reach_count,nonzero_weights\n";
    ofstream all(prefix + "_all_records.csv");
    all << "m,M,blocked,target,reach_count,nonzero_weights\n";

    for(int M=Mmin; M<=Mmax; M++){
        int ms=m_star_ld(M);
        int total=0, blocked=0, reachable=0;
        int firstReach=-1,lastReach=-1;
        string exBlocked="", exReach="";
        for(int m=2;m<=ms;m++){
            total++;
            auto pr = probe_modulus(m,M,R,keep_examples);
            bool isBlocked = !pr.reachable;
            if(isBlocked){
                blocked++;
                if(exBlocked.empty()) exBlocked=nz_to_string(pr.nonzero);
            }else{
                reachable++;
                if(firstReach==-1) firstReach=m;
                lastReach=m;
                if(exReach.empty()) exReach=nz_to_string(pr.nonzero);
                rec << m << ',' << M << ',' << pr.target << ',' << pr.reach_count << ',' << '"' << nz_to_string(pr.nonzero) << '"' << "\n";
            }
            all << m << ',' << M << ',' << (isBlocked?1:0) << ',' << pr.target << ',' << pr.reach_count << ',' << '"' << nz_to_string(pr.nonzero) << '"' << "\n";
        }
        summ << M << ',' << ms << ',' << total << ',' << blocked << ',' << reachable << ',' << firstReach << ',' << lastReach << ',' << '"' << exBlocked << '"' << ',' << '"' << exReach << '"' << "\n";
        if(M%10==0 || reachable>0){
            cerr << "M="<<M<<" mstar="<<ms<<" total="<<total<<" blocked="<<blocked<<" reachable="<<reachable<<"\n";
        }
    }
    cerr << "Done. Wrote " << prefix << "_*\n";
    return 0;
}
