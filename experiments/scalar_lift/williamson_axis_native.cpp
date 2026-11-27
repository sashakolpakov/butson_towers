// Exact axis search for dense symmetric Williamson multiplier charts.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <sys/file.h>
#include <unistd.h>
#include <vector>

using Clock=std::chrono::steady_clock;
using Axis=std::vector<int>;
static std::atomic<bool> stopping(false);
static void interrupt(int){stopping=true;}
static void require(bool value,const std::string& why){if(!value)throw std::runtime_error(why);}
static int rotate(int a,int phase){return a==0?0:1+(a-1+phase)%3;}
static std::string utc(){auto now=std::time(nullptr);std::tm t{};gmtime_r(&now,&t);char b[32];std::strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&t);return b;}
static std::string quote(const std::string& s){std::string v="\"";for(char c:s){if(c=='\\'||c=='\"')v+='\\';if(c=='\n')v+="\\n";else v+=c;}return v+'\"';}
template<class C>static std::string array(const C& c){std::ostringstream out;out<<'[';bool first=true;for(auto v:c){if(!first)out<<',';out<<v;first=false;}return out.str()+']';}
static void atomic_write(const std::filesystem::path& path,const std::string& text){
 auto temporary=path;temporary+=".tmp";
 FILE* f=std::fopen(temporary.c_str(),"wb");require(f!=nullptr,"cannot open checkpoint temporary");
 bool good=std::fwrite(text.data(),1,text.size(),f)==text.size()&&std::fflush(f)==0&&::fsync(::fileno(f))==0;
 int closed=std::fclose(f);require(good&&closed==0,"cannot persist checkpoint temporary");
 require(std::rename(temporary.c_str(),path.c_str())==0,"cannot replace checkpoint");
}

struct Term{int i,j,delta,count;};
struct Problem{
 int p,m,mu,n,action;std::vector<std::vector<Term>> terms;
 std::vector<int> index,phase,representatives;std::vector<std::array<int,4>> targets;
 explicit Problem(const std::string& file){
  std::ifstream in(file);std::string magic;in>>magic;require(magic=="DENSE_WILLIAMSON_AXES_V1","bad input format");
  in>>p>>m>>mu>>n>>action;require(p>2&&m>0&&m%2==0&&(p-1)%m==0&&n==(p-1)/m,"bad chart dimensions");
  require(action==0||(action==1&&(m/2)%3==0),"bad coordinate action");
  int count;in>>count;require(count>0,"no exact row targets");targets.resize(count);for(auto& r:targets)for(int& x:r)in>>x;
  terms.resize(n);for(auto& row:terms){int size;in>>size;require(size>0&&size<=p,"bad coefficient count");row.resize(size);for(auto& t:row){in>>t.i>>t.j>>t.delta>>t.count;require(t.i>=0&&t.i<=t.j&&t.j<n&&t.delta>=0&&t.delta<3&&t.count>0,"bad coefficient");}}
  require(bool(in),"truncated coefficient input");std::string extra;require(!(in>>extra),"extra coefficient input");
  // Rebuild actual residue orbits; no serialized index array is trusted.
  index.assign(p,-1);phase.assign(p,0);std::vector<bool> seen(p,false);seen[0]=true;int next=0;
  for(int r=1;r<p;r++)if(!seen[r]){representatives.push_back(r);int x=r;for(int k=0;k<m;k++){require(!seen[x],"overlapping residue orbit");seen[x]=true;index[x]=next;phase[x]=action?k%3:0;x=int(int64_t(x)*mu%p);}require(x==r,"wrong multiplier order");next++;}
  require(next==n&&std::all_of(seen.begin(),seen.end(),[](bool b){return b;}),"incomplete residue orbit cover");
  for(int x=1;x<p;x++)require(index[x]==index[p-x]&&phase[x]==phase[p-x],"inversion does not fix source states");
  std::vector<std::map<std::array<int,3>,int>> expected(n);
  for(int s=1;s<p;s++){
   std::map<std::array<int,3>,int> actual;
   for(int x=1;x<p;x++){int y=(x+s)%p;if(!y)continue;int i=index[x],j=index[y],delta=(phase[y]-phase[x]+3)%3;if(i>j){std::swap(i,j);delta=(3-delta)%3;}actual[{i,j,delta}]++;}
   int k=index[s];if(s==representatives[k])expected[k]=actual;
   require(actual==expected[k],"expanded shift coefficient mismatch");
  }
  for(int k=0;k<n;k++){
   std::map<std::array<int,3>,int> supplied;int total=0,diag0=0;
   for(auto t:terms[k]){require(!supplied.count({t.i,t.j,t.delta}),"duplicate input coefficient");supplied[{t.i,t.j,t.delta}]=t.count;total+=t.count;if(t.i!=t.j)require(t.count%2==0,"odd offdiagonal coefficient");if(t.i==t.j&&t.delta==0)diag0+=t.count;}
   require(supplied==expected[k]&&total==p-2&&diag0%2==1,"literal coefficient audit failed");
   for(int i=0;i<n;i++)require((supplied[{i,i,1}]+supplied[{i,i,2}])%2==0,"rotated diagonal parity failed");
  }
  for(auto target:targets){int64_t norm=0;for(int v:target)norm+=int64_t(v)*v;require(norm==4*p,"row target norm failed");}
 }
};

struct Evaluation{
 std::vector<int> signs,correlations;std::array<int,4> sums{};
 int64_t correlation_energy=0,row_energy=0,energy=0;bool exact_counts=false,zero=false;
};
static Evaluation evaluate(const Axis& axes,const Problem& q,bool full=true){
 require(int(axes.size())==q.n,"wrong axis vector size");for(int a:axes)require(a>=0&&a<4,"bad axis alphabet");Evaluation e;e.signs.resize(q.n);e.correlations.resize(q.n);
 std::vector<int> diagonal(q.n,0);std::vector<std::vector<Term>> active(q.n);
 for(int k=0;k<q.n;k++){
  int A=0;active[k].reserve(q.terms[k].size()/3);
  for(auto t:q.terms[k])if(axes[t.i]==rotate(axes[t.j],t.delta)){A+=t.count;if(t.i==t.j)diagonal[k]+=t.count;else active[k].push_back(t);}
  require(A%2==1,"unsigned diagonal parity failed");e.signs[k]=1-2*(((A+1)/2)%2);
 }
 e.sums.fill(1);
 for(int i=0;i<q.n;i++){
  int s=e.signs[i];if(q.action==0){for(int j=0;j<4;j++)e.sums[j]+=q.m*s*(1-2*int(j==axes[i]));}
  else if(axes[i]==0){e.sums[0]-=q.m*s;for(int j=1;j<4;j++)e.sums[j]+=q.m*s;}
  else{e.sums[0]+=q.m*s;for(int j=1;j<4;j++)e.sums[j]+=(q.m/3)*s;}
 }
 e.row_energy=INT64_MAX;
 for(auto target:q.targets){int64_t error=0;for(int j=0;j<4;j++){int64_t d=e.sums[j]-target[j];error+=d*d;}e.row_energy=std::min(e.row_energy,error);}
 e.exact_counts=e.row_energy==0;
 if(!full&&!e.exact_counts){e.energy=e.row_energy;return e;}
 for(int k=0;k<q.n;k++){
  int value=diagonal[k]+e.signs[k];for(auto t:active[k])value+=t.count*e.signs[t.i]*e.signs[t.j];
  require(value%4==0,"derived unique signs do not satisfy mod16");e.correlations[k]=value;e.correlation_energy+=int64_t(value)*value;
 }
 e.energy=e.correlation_energy+e.row_energy;e.zero=e.exact_counts&&e.correlation_energy==0;return e;
}

static std::string verify(const Axis& axes,const Evaluation& e,const Problem& q,const std::string& fingerprint){
 require(e.zero,"positive guard called on a nonzero residual");std::vector<std::vector<int>> rows(4,std::vector<int>(q.p,1));
 for(int x=1;x<q.p;x++){int i=q.index[x],a=rotate(axes[i],q.phase[x]);for(int j=0;j<4;j++)rows[j][x]=e.signs[i]*(1-2*int(j==a));}
 std::array<int,4> sums{};for(int j=0;j<4;j++)for(int x=0;x<q.p;x++){require(rows[j][x]==1||rows[j][x]==-1,"non-sign source entry");require(rows[j][x]==rows[j][(q.p-x)%q.p],"source symmetry failed");sums[j]+=rows[j][x];}
 require(sums==e.sums,"independent row sums failed");
 for(int x=1;x<q.p;x++){int prod=1;for(int j=0;j<4;j++)prod*=rows[j][x];require(prod==-1,"product condition failed");}
 std::vector<int> correlations(q.p,0);
 for(int shift=0;shift<q.p;shift++){int value=0;for(int j=0;j<4;j++)for(int x=0;x<q.p;x++)value+=rows[j][x]*rows[j][(x+shift)%q.p];require(value==(shift==0?4*q.p:0),"full integer source correlation failed");correlations[shift]=value;}
 std::ostringstream out;out<<"{\"schema\":\"dense-williamson-axis-witness-v1\",\"p\":"<<q.p<<",\"m\":"<<q.m<<",\"mu\":"<<q.mu<<",\"action\":"<<quote(q.action?"C3":"identity")<<",\"axes\":"<<array(axes)<<",\"signs\":"<<array(e.signs)<<",\"row_sums\":"<<array(sums)<<",\"rows\":[";
 for(int j=0;j<4;j++){if(j)out<<',';out<<array(rows[j]);}out<<"],\"all_integer_correlations\":"<<array(correlations)<<",\"all_expanded_shifts_verified\":true,\"sign_entries\":true,\"symmetric_circulant_rows\":true,\"product_condition\":true,\"fingerprint\":"<<quote(fingerprint)<<",\"new_order_claimed\":false}";return out.str();
}

struct State{
 uint64_t patterns=0,row_pass=0,zeros=0,witnesses=0,accepted=0;double elapsed=0;
 int64_t best=INT64_MAX;Axis current,best_axes;std::mt19937_64 rng;
};
static void save_state(const std::filesystem::path& output,const State& s,const std::string& fingerprint){
 std::ostringstream out;out<<"DENSE_WILLIAMSON_STATE_V1\n"<<fingerprint<<'\n'<<s.patterns<<' '<<s.row_pass<<' '<<s.zeros<<' '<<s.witnesses<<' '<<s.accepted<<' '<<std::setprecision(17)<<s.elapsed<<' '<<s.best<<'\n';
 out<<s.current.size();for(int v:s.current)out<<' '<<v;out<<'\n'<<s.best_axes.size();for(int v:s.best_axes)out<<' '<<v;out<<'\n'<<s.rng<<'\n';atomic_write(output/"checkpoint.txt",out.str());
}
static State load_state(const std::filesystem::path& output,uint64_t seed,const std::string& fingerprint){
 State s;s.rng.seed(seed);std::ifstream in(output/"checkpoint.txt");if(!in)return s;
 std::string magic,old;in>>magic>>old;require(magic=="DENSE_WILLIAMSON_STATE_V1"&&old==fingerprint,"checkpoint input/source fingerprint mismatch");
 in>>s.patterns>>s.row_pass>>s.zeros>>s.witnesses>>s.accepted>>s.elapsed>>s.best;size_t n;in>>n;s.current.resize(n);for(int&v:s.current)in>>v;in>>n;s.best_axes.resize(n);for(int&v:s.best_axes)in>>v;in>>s.rng;require(bool(in),"truncated native state");return s;
}

int main(int argc,char**argv){
 std::filesystem::path output=argc>2?argv[2]:".";
 int reported_p=0,reported_m=0;
 bool owns_output=false;
 try{
  require(argc==14,"usage: engine input output seconds seed mode fingerprint max_patterns max_witnesses T0 Tmin epoch nice");
  std::filesystem::create_directories(output);
  int lock=::open((output/"worker.lock").c_str(),O_CREAT|O_RDWR,0600);
  require(lock>=0&&::flock(lock,LOCK_EX|LOCK_NB)==0,"another native worker owns this output");
  owns_output=true;
  Problem q(argv[1]);reported_p=q.p;reported_m=q.m;double budget=std::stod(argv[3]);uint64_t seed=std::stoull(argv[4]);std::string mode=argv[5],fingerprint=argv[6];uint64_t limit=std::stoull(argv[7]),max_witness=std::stoull(argv[8]);
  double T0=std::stod(argv[9]),Tmin=std::stod(argv[10]);uint64_t epoch=std::stoull(argv[11]);int priority=std::stoi(argv[12]);
  // argv13 is a reserved format version, making positional invocations explicit.
  require(std::string(argv[13])=="1","bad execution format version");require(budget>0&&epoch>0&&T0>0&&Tmin>0,"invalid execution bounds");
  require(mode=="random"||mode=="anneal"||mode=="enumerate","invalid native method");
  if(priority)require(::setpriority(PRIO_PROCESS,0,priority)==0,"cannot set native scheduling priority");
  auto state=load_state(output,seed,fingerprint);const double before=state.elapsed;auto started=Clock::now(),last=started;
  std::signal(SIGTERM,interrupt);std::signal(SIGINT,interrupt);
  uint64_t exhaustive=0;if(mode=="enumerate"){require(q.n<=31,"complete enumeration supports at most31 axes");exhaustive=uint64_t(1)<<(2*q.n);if(!limit||limit>exhaustive)limit=exhaustive;}
  auto elapsed=[&](){return before+std::chrono::duration<double>(Clock::now()-started).count();};
  auto checkpoint=[&](const std::string& phase){state.elapsed=elapsed();save_state(output,state,fingerprint);std::ostringstream out;out<<"{\"schema\":\"dense-williamson-axis-search-v1\",\"p\":"<<q.p<<",\"m\":"<<q.m<<",\"mu\":"<<q.mu<<",\"action\":"<<quote(q.action?"C3":"identity")<<",\"nonzero_orbits\":"<<q.n<<",\"pid\":"<<getpid()<<",\"phase\":"<<quote(phase)<<",\"method\":"<<quote(mode)<<",\"updated_utc\":"<<quote(utc())<<",\"elapsed_seconds\":"<<std::setprecision(17)<<state.elapsed<<",\"budget_seconds\":"<<budget<<",\"fingerprint\":"<<quote(fingerprint)<<",\"totals\":{\"axis_patterns_evaluated\":"<<state.patterns<<",\"exact_row_count_survivors\":"<<state.row_pass<<",\"zero_correlation_candidates\":"<<state.zeros<<",\"verified_witnesses\":"<<state.witnesses<<",\"accepted_moves\":"<<state.accepted<<"},\"best_objective\":"<<state.best<<",\"best_axes\":"<<array(state.best_axes)<<"}\n";atomic_write(output/"status.json",out.str());};
  checkpoint("running");std::string final="budget_exhausted";Evaluation current;if(!state.current.empty())current=evaluate(state.current,q,true);
  while(elapsed()<budget){
   if(stopping||std::filesystem::exists(output/"STOP")){final="stopped";break;}
   if(limit&&state.patterns>=limit){final=(mode=="enumerate"&&limit==exhaustive)?"exhausted":"pattern_limit";break;}
   if(max_witness&&state.witnesses>=max_witness){final="verified_witness";break;}
   Axis proposal(q.n);bool restart=mode!="anneal"||state.current.empty()||state.patterns%epoch==0;
   if(mode=="enumerate"){uint64_t bits=state.patterns;for(int&i:proposal){i=int(bits&3);bits>>=2;}}
   else if(restart){for(int&i:proposal)i=int(state.rng()&3);}
   else{proposal=state.current;int j=int(state.rng()%q.n);proposal[j]=(proposal[j]+1+int(state.rng()%3))%4;}
   auto value=evaluate(proposal,q,mode=="anneal");state.patterns++;state.row_pass+=value.exact_counts;
   if(value.exact_counts&&value.energy<state.best){state.best=value.energy;state.best_axes=proposal;}
   bool accept=restart||mode!="anneal"||value.energy<=current.energy;
   if(!accept){double fraction=double((state.patterns-1)%epoch)/epoch;double T=T0*std::pow(Tmin/T0,fraction);double u=double(state.rng()>>11)*(1.0/9007199254740992.0);accept=u<std::exp(double(current.energy-value.energy)/T);}
   if(accept){state.current=proposal;current=value;state.accepted++;}
   if(value.zero){state.zeros++;std::string witness=verify(proposal,value,q,fingerprint);atomic_write(output/"verified-witness.json",witness+'\n');FILE*f=std::fopen((output/"witnesses.jsonl").c_str(),"ab");require(f!=nullptr,"cannot append witness");require(std::fwrite(witness.data(),1,witness.size(),f)==witness.size(),"cannot write witness");std::fputc('\n',f);std::fflush(f);::fsync(::fileno(f));std::fclose(f);state.witnesses++;checkpoint("running");}
   if(Clock::now()-last>std::chrono::seconds(5)){checkpoint("running");last=Clock::now();}
  }
  if(max_witness&&state.witnesses>=max_witness)final="verified_witness";
  checkpoint(final);std::cout<<"{\"phase\":"<<quote(final)<<",\"patterns\":"<<state.patterns<<",\"witnesses\":"<<state.witnesses<<"}\n";return 0;
 }catch(const std::exception&e){if(owns_output)try{std::string failure="{\"schema\":\"dense-williamson-axis-search-v1\",\"phase\":\"failed\",\"pid\":"+std::to_string(getpid())+",\"p\":"+std::to_string(reported_p)+",\"m\":"+std::to_string(reported_m)+",\"updated_utc\":"+quote(utc())+",\"error\":"+quote(e.what())+"}\n";atomic_write(output/"failure.json",failure);atomic_write(output/"status.json",failure);}catch(...){}std::cerr<<e.what()<<'\n';return 1;}
}
