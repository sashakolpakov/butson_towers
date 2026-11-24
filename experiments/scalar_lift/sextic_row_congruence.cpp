// Exact balanced-triple moves into, then inside, the full mod-two syndrome.
#include <algorithm>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

struct Entry { uint32_t code; int32_t next; };
struct Table {
    std::vector<uint64_t> keys; std::vector<int32_t> heads; std::vector<Entry> entries; size_t mask;
    explicit Table(size_t count){size_t n=1;while(n<2*count+1)n<<=1;keys.resize(n);heads.assign(n,-1);entries.reserve(count);mask=n-1;}
    size_t bucket(uint64_t key)const{uint64_t z=key;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;size_t i=z&mask;while(heads[i]!=-1&&keys[i]!=key)i=(i+1)&mask;return i;}
    void insert(uint64_t key,uint32_t code){size_t i=bucket(key);keys[i]=key;entries.push_back({code,heads[i]});heads[i]=int32_t(entries.size()-1);}
};

int main(int argc,char**argv){
    if(argc!=4){std::cerr<<"CASE CPU_SECONDS SEED\n";return 2;}
    try{
        std::ifstream in(argv[1]);int n,m,quota;in>>n>>m>>quota;
        if(!in||n>1024||m>63||m<1||quota<3||n-quota<3)throw std::runtime_error("invalid dimensions");
        std::vector<int> target(m),bits(n*m),selected(quota);for(int&x:target)in>>x;
        for(int&x:bits){in>>x;if(x<0||x>1)throw std::runtime_error("nonbinary input");}
        for(int&x:selected)in>>x;if(!in)throw std::runtime_error("truncated case");
        double seconds=std::stod(argv[2]);uint64_t seed=std::stoull(argv[3]);if(seconds<=0||seconds>20)throw std::runtime_error("invalid budget");
        std::vector<bool> membership(n,false);for(int x:selected){if(x<0||x>=n||membership[x])throw std::runtime_error("invalid quota subset");membership[x]=true;}
        std::mt19937_64 rng(seed);std::clock_t start=std::clock();auto elapsed=[&](){return double(std::clock()-start)/CLOCKS_PER_SEC;};
        std::vector<uint64_t> signature(n);for(int i=0;i<n;++i)for(int q=0;q<m;++q)signature[i]|=uint64_t(bits[i*m+q])<<q;
        auto residual_for=[&](const std::vector<int>&s){std::vector<int>r(m);for(int q=0;q<m;++q){r[q]=-target[q];for(int x:s)r[q]+=bits[x*m+q];}return r;};
        auto energy_for=[&](const std::vector<int>&r){int e=0;for(int x:r)e+=x*x;return e;};
        auto syndrome_for=[&](const std::vector<int>&r){uint64_t s=0;for(int q=0;q<m;++q)s|=uint64_t(r[q]&1)<<q;return s;};
        auto residual=residual_for(selected);int initial_energy=energy_for(residual);
        bool initial_sector=syndrome_for(residual)==0;
        auto best_selected=selected,best_residual=residual;int best_energy=initial_sector?initial_energy:std::numeric_limits<int>::max();
        std::vector<int> phase_energy;std::vector<uint64_t> phase_matches;std::vector<bool> phase_complete;
        uint64_t scans=0,total_matches=0;int phase=0;
        while(elapsed()<seconds && best_energy!=0){
            residual=residual_for(selected);uint64_t defect=syndrome_for(residual);
            std::vector<int> outside;std::fill(membership.begin(),membership.end(),false);for(int x:selected)membership[x]=true;for(int i=0;i<n;++i)if(!membership[i])outside.push_back(i);
            size_t tuples=size_t(quota)*(quota-1)*(quota-2)/6;Table table(tuples);bool timed=false;
            for(int a=0;a<quota&&!timed;++a)for(int b=a+1;b<quota&&!timed;++b)for(int c=b+1;c<quota;++c){
                table.insert(signature[selected[a]]^signature[selected[b]]^signature[selected[c]],selected[a]|(selected[b]<<10)|(selected[c]<<20));
                if((table.entries.size()&16383)==0&&elapsed()>=seconds){timed=true;break;}
            }
            int candidate_energy=std::numeric_limits<int>::max();uint32_t chosen_out=0;int chosen_a=0,chosen_b=0,chosen_c=0;uint64_t matches=0,ties=0;
            for(size_t a=0;a<outside.size()&&!timed&&candidate_energy;++a)for(size_t b=a+1;b<outside.size()&&!timed&&candidate_energy;++b)for(size_t c=b+1;c<outside.size()&&candidate_energy;++c){
                uint64_t key=signature[outside[a]]^signature[outside[b]]^signature[outside[c]]^defect;size_t bucket=table.bucket(key);
                for(int e=table.heads[bucket];e!=-1;e=table.entries[e].next){
                    ++matches;uint32_t code=table.entries[e].code;int energy=0;
                    for(int q=0;q<m;++q){int z=residual[q]+bits[outside[a]*m+q]+bits[outside[b]*m+q]+bits[outside[c]*m+q]-bits[(code&1023)*m+q]-bits[((code>>10)&1023)*m+q]-bits[((code>>20)&1023)*m+q];if(z&1)throw std::runtime_error("syndrome hash mismatch");energy+=z*z;}
                    if(energy<candidate_energy || (energy==candidate_energy && rng()%++ties==0)){if(energy<candidate_energy)ties=1;candidate_energy=energy;chosen_out=code;chosen_a=outside[a];chosen_b=outside[b];chosen_c=outside[c];}
                }
                if((++scans&16383)==0&&elapsed()>=seconds){timed=true;break;}
            }
            total_matches+=matches;
            if(candidate_energy==std::numeric_limits<int>::max())break;
            for(int k=0;k<3;++k)membership[(chosen_out>>(10*k))&1023]=false;membership[chosen_a]=membership[chosen_b]=membership[chosen_c]=true;
            selected.clear();for(int i=0;i<n;++i)if(membership[i])selected.push_back(i);
            if(int(selected.size())!=quota)throw std::runtime_error("quota lost");
            residual=residual_for(selected);if(syndrome_for(residual)||energy_for(residual)!=candidate_energy)throw std::runtime_error("candidate replay failed");
            phase_energy.push_back(candidate_energy);phase_matches.push_back(matches);phase_complete.push_back(!timed && candidate_energy!=0);++phase;
            if(candidate_energy<=best_energy){best_energy=candidate_energy;best_selected=selected;best_residual=residual;}
        }
        bool sector=best_energy!=std::numeric_limits<int>::max();if(!sector){best_energy=initial_energy;best_selected=selected;best_residual=residual_for(selected);}
        auto checked=residual_for(best_selected);if(checked!=best_residual)throw std::runtime_error("final residual mismatch");
        int modulus=1;if(sector){modulus=2;while(modulus<1024&&std::all_of(checked.begin(),checked.end(),[&](int x){return x%(2*modulus)==0;}))modulus*=2;}
        std::cout<<"{\"schema\":\"sextic-row-congruence-v1\",\"seed\":"<<seed<<",\"cpu_seconds\":"<<elapsed()<<",\"initial_energy\":"<<initial_energy<<",\"initial_mod2_sector\":"<<(initial_sector?"true":"false")<<",\"best_energy\":"<<best_energy<<",\"all_mod2_passed\":"<<(sector?"true":"false")<<",\"passed_modulus\":"<<modulus<<",\"scanned_triples\":"<<scans<<",\"syndrome_matches\":"<<total_matches<<",\"phases\":[";
        for(size_t j=0;j<phase_energy.size();++j){if(j)std::cout<<',';std::cout<<"{\"energy\":"<<phase_energy[j]<<",\"matches\":"<<phase_matches[j]<<",\"complete\":"<<(phase_complete[j]?"true":"false")<<'}';}
        std::cout<<"],\"selected\":[";for(size_t j=0;j<best_selected.size();++j){if(j)std::cout<<',';std::cout<<best_selected[j];}std::cout<<"],\"residual\":[";for(int q=0;q<m;++q){if(q)std::cout<<',';std::cout<<checked[q];}std::cout<<"]}\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
