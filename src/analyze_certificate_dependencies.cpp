#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;

static vector<int> primes_upto(int n){
    vector<int> primes; vector<bool> is(n+1,true); if(n>=0) is[0]=false; if(n>=1) is[1]=false;
    for(int i=2;i<=n;i++) if(is[i]){ primes.push_back(i); if((long long)i*i<=n) for(long long j=1LL*i*i;j<=n;j+=i) is[(int)j]=false; }
    return primes;
}

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

struct PairData{
    int m, M;
    cpp_int L;
    vector<cpp_int> term; // L/t for t=m..M
    cpp_int target; // sum L/t - L
};

static PairData make_pair_data(int m, int M){
    PairData d; d.m=m; d.M=M; d.L=lcm_range_int(m,M);
    d.target = -d.L;
    for(int t=m;t<=M;t++){
        cpp_int q = d.L / t;
        d.term.push_back(q);
        d.target += q;
    }
    return d;
}

static bool target_reachable_mod(const PairData &d, int r, const vector<char> *forcedPtr=nullptr, string *fail_reason=nullptr, int *at_t=nullptr){
    int len = d.M-d.m+1;
    int target = cpp_mod_int(d.target, r);
    vector<int> weights(len);
    for(int i=0;i<len;i++) weights[i] = cpp_mod_int(d.term[i], r);
    vector<unsigned char> dp0(r,0), dp1(r,0), ndp0(r,0), ndp1(r,0);
    dp0[0]=1;
    for(int i=0;i<len;i++){
        int t = d.m+i;
        fill(ndp0.begin(), ndp0.end(), 0);
        fill(ndp1.begin(), ndp1.end(), 0);
        bool forced = forcedPtr && (*forcedPtr)[i];
        int w = weights[i];
        if(forced){
            // must skip, only if previous not skipped
            for(int res=0; res<r; ++res) if(dp0[res]) ndp1[(res+w)%r]=1;
        }else{
            for(int res=0; res<r; ++res){
                if(dp0[res]||dp1[res]) ndp0[res]=1;
                if(t!=d.m && t!=d.M && dp0[res]) ndp1[(res+w)%r]=1;
            }
        }
        dp0.swap(ndp0); dp1.swap(ndp1);
        bool any=false; for(int res=0;res<r;res++) if(dp0[res]||dp1[res]){any=true;break;}
        if(!any){ if(fail_reason) *fail_reason="no_state_survives"; if(at_t) *at_t=t; return false; }
    }
    bool reach = dp0[target] || dp1[target];
    if(!reach && fail_reason) *fail_reason="target_unreachable";
    return reach;
}

static vector<int> candidate_moduli(int limit, bool rich_squarefree){
    set<int> mods;
    auto primes = primes_upto(limit);
    for(int p: primes) mods.insert(p);
    auto small = primes_upto(min(limit, 97));
    for(int p: small){
        long long x=1LL*p*p;
        while(x<=limit){ mods.insert((int)x); x*=p; }
    }
    vector<int> sf_primes;
    if(rich_squarefree) sf_primes = {2,3,5,7,11,13,17,19,23,29};
    else sf_primes = {2,3,5,7,11,13};
    int n=sf_primes.size();
    for(int mask=1; mask<(1<<n); ++mask){
        long long prod=1;
        bool ok=true;
        for(int i=0;i<n;i++) if(mask&(1<<i)){
            prod *= sf_primes[i];
            if(prod>limit){ok=false; break;}
        }
        if(ok) mods.insert((int)prod);
    }
    vector<int> v(mods.begin(), mods.end());
    sort(v.begin(), v.end());
    return v;
}

static int find_mod_block(const PairData &d, const vector<int> &mods, int &checked){
    checked=0;
    for(int r: mods){
        checked++;
        if(!target_reachable_mod(d,r,nullptr,nullptr,nullptr)) return r;
    }
    return -1;
}

static long long B_j(int j){
    long long L=1;
    auto gcdll = [](long long a,long long b){ while(b){long long t=a%b;a=b;b=t;} return a; };
    for(int a=1;a<=j;a++) L = L / gcdll(L,a) * a;
    long long s=0; for(int a=1;a<=j;a++) s += L/a;
    return s;
}

struct ForcedInfo{
    vector<int> forced_positions;
    map<int, vector<array<int,5>>> reasons; // x -> {j,p,a,x,Bj}
    vector<char> forced_mask;
};

static ForcedInfo forced_positions(int m, int M, int J, const vector<int> &prime_list){
    set<int> forced;
    map<int, vector<array<int,5>>> reasons;
    for(int j=1;j<=J;j++){
        long long Bj = B_j(j);
        int lo = M/(j+1)+1;
        int hi = M/j;
        for(int p: prime_list){
            if(p<max(2,lo)) continue;
            if(p>hi) break;
            if(p<=Bj) continue;
            for(int a=1;a<=j;a++){
                int x=a*p;
                if(m<=x && x<=M){
                    forced.insert(x);
                    reasons[x].push_back({j,p,a,x,(int)Bj});
                }
            }
        }
    }
    ForcedInfo info;
    info.forced_positions.assign(forced.begin(), forced.end());
    info.reasons = std::move(reasons);
    info.forced_mask.assign(M-m+1, 0);
    for(int x: info.forced_positions) info.forced_mask[x-m]=1;
    return info;
}

static pair<int,int> adjacent_forced_pair(const ForcedInfo &fi){
    for(size_t i=1;i<fi.forced_positions.size();i++){
        if(fi.forced_positions[i] == fi.forced_positions[i-1]+1) return {fi.forced_positions[i-1], fi.forced_positions[i]};
    }
    return {-1,-1};
}

static tuple<bool,int,int,int,int> dense_interval(const ForcedInfo &fi, int m, int M, int min_len){
    int n=M-m+1;
    vector<int> pref(n+1,0);
    for(int i=0;i<n;i++) pref[i+1]=pref[i]+(fi.forced_mask[i]?1:0);
    bool found=false; int bestA=-1,bestB=-1,bestCnt=0,bestLen=0;
    for(int a=0;a<n;a++){
        for(int b=a+min_len-1;b<n;b++){
            int len=b-a+1;
            int cnt=pref[b+1]-pref[a];
            int max_noadj=(len+1)/2; // ceil(len/2)
            if(cnt>max_noadj){
                if(!found || len>bestLen || (len==bestLen && cnt>bestCnt)){
                    found=true; bestA=m+a; bestB=m+b; bestCnt=cnt; bestLen=len;
                }
            }
        }
    }
    return {found,bestA,bestB,bestCnt,bestLen};
}

static int m_star_ld(int M){
    long double s=0;
    for(int t=M;t>=2;t--){
        s += 1.0L/t;
        if(s >= 1.0L - 1e-18L) return t;
    }
    return 2;
}

static string json_escape(const string &s){ string o; for(char c:s){ if(c=='"') o+="\\\""; else if(c=='\\') o+="\\\\"; else o+=c;} return o; }

int main(int argc, char** argv){
    int Mmin=10, Mmax=200, primeLimit=1000, mixedLimit=10000, J=12, denseMinLen=3;
    bool rich=true;
    string prefix="/mnt/data/cert_diag";
    for(int i=1;i<argc;i++){
        string a=argv[i]; auto next=[&](){ return string(argv[++i]); };
        if(a=="--M-min") Mmin=stoi(next());
        else if(a=="--M-max") Mmax=stoi(next());
        else if(a=="--prime-limit") primeLimit=stoi(next());
        else if(a=="--mixed-limit") mixedLimit=stoi(next());
        else if(a=="--J") J=stoi(next());
        else if(a=="--dense-min-len") denseMinLen=stoi(next());
        else if(a=="--prefix") prefix=next();
        else if(a=="--not-rich") rich=false;
    }
    vector<int> primesMixed = primes_upto(max(mixedLimit, Mmax+10));
    vector<int> directPrimes; for(int p: primesMixed) if(p<=primeLimit) directPrimes.push_back(p);
    vector<int> allMods = candidate_moduli(mixedLimit, rich);
    vector<int> afterDirectMods; for(int r: allMods) if(!(r<=primeLimit && binary_search(directPrimes.begin(), directPrimes.end(), r))) afterDirectMods.push_back(r);

    ofstream rec(prefix + "_records.csv");
    rec << "m,M,category,modulus,forced_count,adjacent_pair,dense_interval,dense_count,dense_len,forced_modulus\n";
    ofstream summ(prefix + "_summary_by_M.csv");
    summ << "M,total,direct_prime,mixed_or_large_modulus,dense_forced_holes,linked_prime_pair,forced_modular,endpoint_forced,no_state_forced,survivor\n";

    map<string,long long> overall;
    map<string,long long> pairCounter, signatureCounter;
    vector<tuple<int,int,string,string>> examples;
    for(int M=Mmin; M<=Mmax; ++M){
        int ms=m_star_ld(M);
        map<string,int> counts;
        for(int m=2; m<=ms; ++m){
            PairData d = make_pair_data(m,M);
            int checked=0;
            int mod = find_mod_block(d, directPrimes, checked);
            string cat; int forced_count=0; string adj=""; string denseStr=""; int denseCnt=0,denseLen=0; int forcedMod=-1;
            if(mod!=-1){ cat="direct_prime"; }
            else{
                mod = find_mod_block(d, afterDirectMods, checked);
                if(mod!=-1){ cat="mixed_or_large_modulus"; }
                else{
                    ForcedInfo fi = forced_positions(m,M,J,primesMixed);
                    forced_count=fi.forced_positions.size();
                    // endpoint forced
                    if(fi.forced_mask[0] || fi.forced_mask[M-m]){
                        cat="endpoint_forced";
                        adj = fi.forced_mask[0] ? to_string(m) : to_string(M);
                    } else {
                        auto [dfound,A,B,cnt,len] = dense_interval(fi,m,M,denseMinLen);
                        auto ap = adjacent_forced_pair(fi);
                        if(dfound){
                            cat="dense_forced_holes"; denseStr=to_string(A)+"-"+to_string(B); denseCnt=cnt; denseLen=len;
                        } else if(ap.first!=-1){
                            cat="linked_prime_pair"; adj=to_string(ap.first)+"-"+to_string(ap.second);
                            // signatures
                            for(auto &L: fi.reasons[ap.first]) for(auto &R: fi.reasons[ap.second]){
                                string sig=to_string(L[2])+"*"+to_string(L[1])+"+1="+to_string(R[2])+"*"+to_string(R[1])+" | j "+to_string(L[0])+"->"+to_string(R[0]);
                                signatureCounter[sig]++;
                            }
                            pairCounter[adj]++;
                        } else {
                            // Try forced-modular automaton with candidate moduli, which uses forced positions but not immediate adjacent/local contradiction.
                            bool blocked=false; string reason; int at_t=-1;
                            for(int r: allMods){
                                reason.clear(); at_t=-1;
                                if(!target_reachable_mod(d,r,&fi.forced_mask,&reason,&at_t)){
                                    forcedMod=r;
                                    if(reason=="no_state_survives") cat="no_state_forced"; else cat="forced_modular";
                                    blocked=true; break;
                                }
                            }
                            if(!blocked) cat="survivor";
                        }
                    }
                }
            }
            counts[cat]++; overall[cat]++;
            rec << m << ',' << M << ',' << cat << ',' << mod << ',' << forced_count << ',' << adj << ',' << denseStr << ',' << denseCnt << ',' << denseLen << ',' << forcedMod << "\n";
            if((cat=="linked_prime_pair" || cat=="survivor" || cat=="forced_modular") && examples.size()<50){
                examples.push_back({m,M,cat,adj});
            }
        }
        int total=0; for(auto &kv: counts) total+=kv.second;
        summ << M << ',' << total << ',' << counts["direct_prime"] << ',' << counts["mixed_or_large_modulus"] << ',' << counts["dense_forced_holes"] << ',' << counts["linked_prime_pair"] << ',' << counts["forced_modular"] << ',' << counts["endpoint_forced"] << ',' << counts["no_state_forced"] << ',' << counts["survivor"] << "\n";
        if(M%10==0 || counts["survivor"]>0){
            cerr << "M="<<M<<" total="<<total<<" prime="<<counts["direct_prime"]<<" mixed="<<counts["mixed_or_large_modulus"]<<" dense="<<counts["dense_forced_holes"]<<" linked="<<counts["linked_prime_pair"]<<" fmod="<<counts["forced_modular"]<<" surv="<<counts["survivor"]<<"\n";
        }
    }
    ofstream js(prefix + "_report.json");
    js << "{\n";
    js << "  \"M_min\": "<<Mmin<<", \"M_max\": "<<Mmax<<", \"prime_limit\": "<<primeLimit<<", \"mixed_limit\": "<<mixedLimit<<", \"J\": "<<J<<", \"dense_min_len\": "<<denseMinLen<<",\n";
    js << "  \"overall_counts\": {"; bool first=true; long long totalAll=0; for(auto &kv: overall) totalAll+=kv.second; overall["total"]=totalAll; for(auto &kv: overall){ if(!first) js<<", "; first=false; js << "\""<<kv.first<<"\": "<<kv.second; } js << "},\n";
    js << "  \"top_linked_pairs\": {"; first=true; vector<pair<string,long long>> pc(pairCounter.begin(),pairCounter.end()); sort(pc.begin(),pc.end(),[](auto&a,auto&b){return a.second>b.second;}); int c=0; for(auto &kv: pc){ if(c++>=20) break; if(!first) js<<", "; first=false; js << "\""<<kv.first<<"\": "<<kv.second;} js << "},\n";
    js << "  \"top_signatures\": {"; first=true; vector<pair<string,long long>> sc(signatureCounter.begin(),signatureCounter.end()); sort(sc.begin(),sc.end(),[](auto&a,auto&b){return a.second>b.second;}); c=0; for(auto &kv: sc){ if(c++>=20) break; if(!first) js<<", "; first=false; js << "\""<<json_escape(kv.first)<<"\": "<<kv.second;} js << "}\n";
    js << "}\n";
    cerr << "Done. Wrote " << prefix << "_*.\n";
    return 0;
}
