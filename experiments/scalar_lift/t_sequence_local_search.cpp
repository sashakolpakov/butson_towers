#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Exact incremental periodic-correlation search. Every proposed move preserves
// the signed sums of head and tail orbits. The optional free-head moves
// exchange a balanced head pair with a balanced tail pair. Only integral correlation
// values enter the objective; annealing probabilities guide proposals only.
struct Engine {
    int p,m,mu,n,L,headSum,tailSum,absoluteWeight,parityWeight;
    bool freeHeads=false;
    int minHeads=std::numeric_limits<int>::max(),maxHeads=-1;
    uint64_t conversionProposals=0,conversionAccepted=0;
    std::vector<int16_t> H,K;
    std::vector<int> s,c,F,best,initial;
    std::vector<std::array<int,8>> tabu;
    std::mt19937_64 rng;
    long long energy=0,bestEnergy=std::numeric_limits<long long>::max(),bestEvenEnergy=std::numeric_limits<long long>::max();
    uint64_t proposals=0,accepted=0,restarts=0,iterations=0;
    std::chrono::steady_clock::time_point started;
    std::string output;

    size_t hi(int i,int a,int l=0) const {return (size_t(i)*8+a)*L+l;}
    size_t ki(int i,int j,int a,int b,int l=0) const {return (((size_t(i)*n+j)*8+a)*8+b)*L+l;}
    int pick(int lim){return int(rng()%lim);}
    double elapsed()const{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
    Engine(const std::string& table,uint64_t seed,const std::string& target,int weight=0,int parity=0,bool free=false):absoluteWeight(weight),parityWeight(parity),freeHeads(free),rng(seed),output(target){
        if(weight<0||parity<0)throw std::runtime_error("objective weights must be nonnegative");
        std::ifstream f(table,std::ios::binary); int32_t a[7]; f.read((char*)a,sizeof a);
        if(!f)throw std::runtime_error("table header missing");
        p=a[0];m=a[1];mu=a[2];n=a[3];L=a[4];headSum=a[5];tailSum=a[6];
        if(n<2||L<1||n>2000||L>2000)throw std::runtime_error("invalid table dimensions");
        H.resize(size_t(n)*8*L);K.resize(size_t(n)*n*64*L);
        f.read((char*)H.data(),H.size()*sizeof(int16_t));f.read((char*)K.data(),K.size()*sizeof(int16_t));
        if(!f||f.peek()!=EOF)throw std::runtime_error("table size differs from dimensions");
        s.resize(n);c.resize(L);F.resize(size_t(n)*8*L);tabu.resize(n);
        started=std::chrono::steady_clock::now();
    }
    int headCount()const{return std::count_if(s.begin()+1,s.end(),[](int a){return a/2==0;});}
    void observeHeads(){int h=headCount();minHeads=std::min(minHeads,h);maxHeads=std::max(maxHeads,h);}
    bool conversion(int i,int j,int& a,int& b){
        if(i<1||j<1||i==j||!((s[i]^s[j])&1)||bool(s[i]/2)!=bool(s[j]/2))return false;
        if(s[i]/2){a=s[i]&1;b=s[j]&1;}
        else {a=2*(1+pick(3))+(s[i]&1);b=2*(1+pick(3))+(s[j]&1);}
        return true;
    }
    void rebuild(){
        std::fill(c.begin(),c.end(),0);std::fill(F.begin(),F.end(),0);
        for(int i=0;i<n;++i)for(int a=0;a<8;++a)for(int l=0;l<L;++l)F[hi(i,a,l)]=H[hi(i,a,l)];
        for(int i=0;i<n;++i){
            for(int l=0;l<L;++l)c[l]+=H[hi(i,s[i],l)];
            for(int j=i+1;j<n;++j)for(int l=0;l<L;++l)c[l]+=K[ki(i,j,s[i],s[j],l)];
            for(int j=0;j<n;++j)if(i!=j)for(int a=0;a<8;++a){
                const auto at=ki(i,j,a,s[j]);const auto dest=hi(i,a);
                for(int l=0;l<L;++l)F[dest+l]+=K[at+l];
            }
        }
        energy=0;for(int z:c)energy+=1LL*z*z+1LL*absoluteWeight*std::abs(z)+1LL*parityWeight*(std::abs(z)%2);
    }
    void initialize(int heads){
        int hn=(headSum-1)/m,tn=tailSum/(m/3),tails=n-1-heads;
        if(heads<0||tails<0||(heads+hn)%2||(tails+tn)%2||heads<std::abs(hn)||tails<std::abs(tn))
            throw std::runtime_error("head count incompatible with exact row sums");
        s[0]=0;int k=1;
        for(int i=0;i<(heads+hn)/2;++i)s[k++]=0;
        for(int i=0;i<(heads-hn)/2;++i)s[k++]=1;
        for(int i=0;i<(tails+tn)/2;++i)s[k++]=2*(1+pick(3));
        for(int i=0;i<(tails-tn)/2;++i)s[k++]=2*(1+pick(3))+1;
        std::shuffle(s.begin()+1,s.end(),rng);
        for(auto& t:tabu)t.fill(0);
        rebuild();observeHeads();++restarts;
    }
    long long delta(int i,int a,int j=-1,int b=-1,std::vector<int>* values=nullptr)const{
        const int oldi=s[i],oldj=j>=0?s[j]:0;
        long long change=0;
        if(values)values->resize(L);
        for(int l=0;l<L;++l){
            int d=F[hi(i,a,l)]-F[hi(i,oldi,l)];
            if(j>=0)d+=F[hi(j,b,l)]-F[hi(j,oldj,l)]
                +K[ki(i,j,a,b,l)]-K[ki(i,j,a,oldj,l)]-K[ki(i,j,oldi,b,l)]+K[ki(i,j,oldi,oldj,l)];
            change+=1LL*d*(2*c[l]+d)+1LL*absoluteWeight*(std::abs(c[l]+d)-std::abs(c[l]))
                +1LL*parityWeight*((std::abs(c[l]+d)%2)-(std::abs(c[l])%2));
            if(values)(*values)[l]=d;
        }
        return change;
    }
    void mutate(int i,int a){
        int old=s[i];if(old==a)return;
        for(int k=0;k<n;++k)if(k!=i)for(int b=0;b<8;++b){
            auto dest=hi(k,b),add=ki(k,i,b,a),sub=ki(k,i,b,old);
            for(int l=0;l<L;++l)F[dest+l]+=K[add+l]-K[sub+l];
        }
        s[i]=a;
    }
    void apply(int i,int a,int j=-1,int b=-1){
        int headChange=int(a/2==0)-int(s[i]/2==0);
        if(j>=0)headChange+=int(b/2==0)-int(s[j]/2==0);
        std::vector<int> d;long long change=delta(i,a,j,b,&d);
        mutate(i,a);if(j>=0)mutate(j,b);
        for(int l=0;l<L;++l)c[l]+=d[l];
        energy+=change;++accepted;
        if(freeHeads&&headChange){++conversionAccepted;observeHeads();}
    }
    void assertCounts()const{
        int h=1,t=0;
        for(int i=1;i<n;++i){int sign=1-2*(s[i]&1);if(s[i]/2==0)h+=m*sign;else t+=(m/3)*sign;}
        if(h!=headSum||t!=tailSum||s[0]!=0)throw std::runtime_error("exact row sums changed");
    }
    void checkpoint(bool final=false){
        if(energy<bestEnergy){bestEnergy=energy;best=s;assertCounts();}
        else if(!final)return;
        const auto save=s;const auto oldCorr=c;const auto oldEnergy=energy;
        s=best;rebuild();assertCounts();
        if(energy!=bestEnergy)throw std::runtime_error("checkpoint energy does not replay");
        long long square=0;for(int z:c)square+=1LL*z*z;
        const std::string temp=output+".tmp";std::ofstream f(temp);
        f<<"{\"p\":"<<p<<",\"m\":"<<m<<",\"mu\":"<<mu<<",\"best_energy\":"<<square<<",\"search_cost\":"<<bestEnergy<<",\"absolute_weight\":"<<absoluteWeight<<",\"parity_weight\":"<<parityWeight
          <<",\"free_heads\":"<<(freeHeads?"true":"false")<<",\"head_orbits\":"<<headCount()
          <<",\"minimum_heads_visited\":"<<minHeads<<",\"maximum_heads_visited\":"<<maxHeads
          <<",\"conversion_proposals\":"<<conversionProposals<<",\"conversion_moves_accepted\":"<<conversionAccepted
          <<",\"exact_correlation_zero\":"<<(bestEnergy==0?"true":"false")<<",\"elapsed_seconds\":"<<elapsed()
          <<",\"proposals\":"<<proposals<<",\"accepted\":"<<accepted<<",\"restarts\":"<<restarts<<",\"iterations\":"<<iterations<<",\"states\":[";
        for(int i=0;i<n;++i)f<<(i?",":"")<<best[i];f<<"],\"correlations\":[";
        for(int l=0;l<L;++l)f<<(l?",":"")<<c[l];f<<"],\"row_sums\":["<<headSum<<","<<tailSum<<","<<tailSum<<","<<tailSum<<"]}\n";f.close();
        if(std::rename(temp.c_str(),output.c_str()))throw std::runtime_error("checkpoint rename failed");
        std::cout<<"{\"best_energy\":"<<square<<",\"search_cost\":"<<bestEnergy<<",\"seconds\":"<<elapsed()<<",\"proposals\":"<<proposals<<",\"restarts\":"<<restarts<<"}"<<std::endl;
        s=save;rebuild();
        if(c!=oldCorr||energy!=oldEnergy)throw std::runtime_error("checkpoint restoration changed correlations");
    }
    void checkpointEven(){
        long long square=0;
        for(int z:c){if(std::abs(z)%2)return;square+=1LL*z*z;}
        if(square>=bestEvenEnergy)return;
        bestEvenEnergy=square;assertCounts();
        const std::string name=output+".even.json",temp=name+".tmp";std::ofstream f(temp);
        f<<"{\"p\":"<<p<<",\"m\":"<<m<<",\"mu\":"<<mu<<",\"best_energy\":"<<square
         <<",\"free_heads\":"<<(freeHeads?"true":"false")<<",\"head_orbits\":"<<headCount()
         <<",\"exact_correlation_zero\":"<<(square==0?"true":"false")
         <<",\"all_class_correlations_even\":true,\"elapsed_seconds\":"<<elapsed()<<",\"states\":[";
        for(int i=0;i<n;++i)f<<(i?",":"")<<s[i];f<<"],\"correlations\":[";
        for(int l=0;l<L;++l)f<<(l?",":"")<<c[l];
        f<<"],\"row_sums\":["<<headSum<<","<<tailSum<<","<<tailSum<<","<<tailSum<<"]}\n";f.close();
        if(std::rename(temp.c_str(),name.c_str()))throw std::runtime_error("even checkpoint rename failed");
        std::cout<<"{\"best_all_even_energy\":"<<square<<",\"seconds\":"<<elapsed()<<"}"<<std::endl;
    }
    void deltaChecks(){
        for(int i=0;i<n;++i)s[i]=pick(8);s[0]=0;rebuild();
        for(int trial=0;trial<500;++trial){
            int i=1+pick(n-1),j=1+pick(n-1);if(j==i)j=-1;
            int a=pick(8),b=pick(8);auto old=c;long long e=energy;
            long long change=delta(i,a,j,b);apply(i,a,j,b);auto updated=c;auto updatedF=F;
            rebuild();
            if(c!=updated||F!=updatedF||energy!=e+change)throw std::runtime_error("incremental delta check failed");
        }
        std::cout<<"{\"incremental_random_move_checks\":500,\"passed\":true}"<<std::endl;
    }
    void conversionChecks(int heads){
        initialize(heads);
        int increases=0,decreases=0;
        for(int trial=0;trial<300;trial++){
            int i=-1,j=-1,a=0,b=0;
            for(int pass=0;pass<2&&i<0;pass++){
                bool wantTail=((trial%2)==0)^bool(pass);
                for(int u=1;u<n&&i<0;u++)for(int v=u+1;v<n;v++){
                    if(bool(s[u]/2)!=wantTail)continue;
                    if(conversion(u,v,a,b)){i=u;j=v;break;}
                }
            }
            if(i<0)throw std::runtime_error("no balanced conversion available for the control");
            int before=headCount();long long oldEnergy=energy,change=delta(i,a,j,b);
            apply(i,a,j,b);assertCounts();
            int after=headCount();
            if(after==before+2)increases++;else if(after==before-2)decreases++;
            else throw std::runtime_error("conversion did not change the head count by two");
            auto updated=c;auto updatedF=F;rebuild();
            if(c!=updated||F!=updatedF||energy!=oldEnergy+change)
                throw std::runtime_error("balanced conversion incremental check failed");
        }
        if(!increases||!decreases)throw std::runtime_error("both conversion directions were not checked");
        std::cout<<"{\"balanced_conversion_checks\":300,\"increases\":"<<increases
                 <<",\"decreases\":"<<decreases<<",\"exact_DC_preserved\":true}"<<std::endl;
    }
    void anneal(double seconds,int heads){
        initialize(heads);if(!initial.empty()){s=initial;rebuild();assertCounts();}checkpoint();checkpointEven();double cycleStart=elapsed();uint64_t cycle=0;
        while(elapsed()<seconds&&bestEnergy){
            if(cycle>=100000){
                if(pick(4)==0)initialize(heads);
                else {s=best;rebuild();for(int k=0;k<3+pick(15);++k){int i=1+pick(n-1),j=1+pick(n-1);if(i!=j)apply(i,s[j],j,s[i]);}}
                cycle=0;cycleStart=elapsed();
            }
            ++cycle;++proposals;++iterations;
            int i=1+pick(n-1),j=-1,a=s[i],b=-1;
            if(freeHeads&&pick(5)==0){
                bool found=false;
                for(int attempt=0;attempt<40&&!found;attempt++){
                    i=1+pick(n-1);j=1+pick(n-1);found=conversion(i,j,a,b);
                }
                if(!found)continue;
                ++conversionProposals;
            }
            else if(s[i]/2&&pick(5)<2)a=2*(1+((s[i]/2-1+1+pick(2))%3))+(s[i]&1);
            else {
                j=1+pick(n-1);if(i==j||s[i]==s[j])continue;
                if(s[i]/2&&s[j]/2&&((s[i]^s[j])&1)&&pick(2)){a=s[i]^1;b=s[j]^1;}
                else {a=s[j];b=s[i];}
            }
            long long d=delta(i,a,j,b);
            double phase=double(cycle)/100000.0;
            double top=(bestEnergy<100?24.0:80.0)*(1.0+absoluteWeight/4.0+parityWeight/4.0);
            double T=0.25+top*std::pow(1-phase,3);
            if(d<=0||double(rng()>>11)*(1.0/9007199254740992.0)<std::exp(-double(d)/T)){
                apply(i,a,j,b);checkpoint();checkpointEven();
            }
        }
        checkpoint(true);
    }
    void tabuSearch(double seconds,int heads){
        initialize(heads);if(!initial.empty()){s=initial;rebuild();assertCounts();}checkpoint();checkpointEven();uint64_t stale=0;
        while(elapsed()<seconds&&bestEnergy){
            ++iterations;++stale;
            long long bestDelta=std::numeric_limits<long long>::max();int bi=-1,bj=-1,ba=0,bb=0,ties=0;
            auto consider=[&](int i,int a,int j,int b){
                ++proposals;
                int hc=int(a/2==0)-int(s[i]/2==0);
                if(j>=0)hc+=int(b/2==0)-int(s[j]/2==0);
                if(hc)++conversionProposals;
                long long d=delta(i,a,j,b);
                bool forbidden=tabu[i][a]>int(iterations)||(j>=0&&tabu[j][b]>int(iterations));
                if(forbidden&&energy+d>=bestEnergy)return;
                if(d<bestDelta){bestDelta=d;bi=i;ba=a;bj=j;bb=b;ties=1;}
                else if(d==bestDelta&&pick(++ties)==0){bi=i;ba=a;bj=j;bb=b;}
            };
            for(int i=1;i<n;++i){
                if(s[i]/2)for(int a=2;a<8;++a)if((a&1)==(s[i]&1)&&a!=s[i])consider(i,a,-1,-1);
                for(int j=i+1;j<n;++j){
                    if(s[i]!=s[j])consider(i,s[j],j,s[i]);
                    if(s[i]/2&&s[j]/2&&s[i]/2!=s[j]/2&&((s[i]^s[j])&1))
                        consider(i,s[i]^1,j,s[j]^1);
                    if(freeHeads&&((s[i]^s[j])&1)&&bool(s[i]/2)==bool(s[j]/2)){
                        if(s[i]/2)consider(i,s[i]&1,j,s[j]&1);
                        else for(int a=2;a<8;a++)if((a&1)==(s[i]&1))
                            for(int b=2;b<8;b++)if((b&1)==(s[j]&1))consider(i,a,j,b);
                    }
                }
            }
            if(bi<0||stale>5000){initialize(heads);stale=0;continue;}
            int oi=s[bi],oj=bj>=0?s[bj]:0;
            apply(bi,ba,bj,bb);checkpointEven();int tenure=5+pick(8);
            tabu[bi][oi]=int(iterations)+tenure;if(bj>=0)tabu[bj][oj]=int(iterations)+tenure;
            if(energy<bestEnergy){checkpoint();stale=0;}
        }
        checkpoint(true);
    }
};
int main(int argc,char**argv){try{
    if(argc<6)throw std::runtime_error("usage: binary table output seconds seed heads [anneal|tabu|check]");
    Engine e(argv[1],std::stoull(argv[4]),argv[2],argc>7?std::stoi(argv[7]):0,argc>8?std::stoi(argv[8]):0,argc>10&&std::stoi(argv[10])!=0);std::string mode=argc>6?argv[6]:"tabu";
    if(argc>9&&std::string(argv[9])!="-"){std::ifstream f(argv[9]);e.initial.resize(e.n);for(int& z:e.initial)if(!(f>>z)||z<0||z>7)throw std::runtime_error("invalid initial states");}
    if(mode=="check"){e.deltaChecks();if(e.freeHeads)e.conversionChecks(std::stoi(argv[5]));}
    else if(mode=="anneal")e.anneal(std::stod(argv[3]),std::stoi(argv[5]));
    else if(mode=="tabu")e.tabuSearch(std::stod(argv[3]),std::stoi(argv[5]));
    else throw std::runtime_error("unknown search mode");
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
