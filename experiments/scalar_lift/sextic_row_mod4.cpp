// Compose disjoint parity-kernel triple moves to lift the full syndrome to mod4.
#include <algorithm>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

struct Entry{uint32_t code;int32_t next;};
struct Table{
    std::vector<uint64_t>keys;std::vector<int32_t>heads;std::vector<Entry>entries;size_t mask;
    explicit Table(size_t count){size_t n=1;while(n<2*count+1)n<<=1;keys.resize(n);heads.assign(n,-1);entries.reserve(count);mask=n-1;}
    size_t bucket(uint64_t k)const{uint64_t z=k;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;size_t i=z&mask;while(heads[i]!=-1&&keys[i]!=k)i=(i+1)&mask;return i;}
    void insert(uint64_t k,uint32_t c){size_t i=bucket(k);keys[i]=k;entries.push_back({c,heads[i]});heads[i]=int32_t(entries.size()-1);}
};
struct Move{uint64_t positive=0,negative=0;uint16_t out[3],in[3];uint64_t syndrome()const{return positive|negative;}};
bool disjoint(const Move&a,const Move&b){for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(a.out[i]==b.out[j]||a.in[i]==b.in[j])return false;return true;}

int main(int argc,char**argv){
    if(argc!=4){std::cerr<<"CASE CPU_SECONDS SEED\n";return 2;}
    try{
        std::ifstream input(argv[1]);int n,m,quota;input>>n>>m>>quota;
        if(!input||n>1024||m>63||m<1||quota<3||n-quota<3)throw std::runtime_error("invalid dimensions");
        std::vector<int>target(m),bits(n*m),selected(quota);for(int&x:target)input>>x;for(int&x:bits){input>>x;if(x<0||x>1)throw std::runtime_error("nonbinary");}for(int&x:selected)input>>x;if(!input)throw std::runtime_error("truncated case");
        double seconds=std::stod(argv[2]);uint64_t seed=std::stoull(argv[3]);if(seconds<=0||seconds>10)throw std::runtime_error("invalid CPU limit");
        std::clock_t start=std::clock();auto elapsed=[&](){return double(std::clock()-start)/CLOCKS_PER_SEC;};std::mt19937_64 rng(seed);
        std::vector<bool>membership(n,false);for(int x:selected){if(x<0||x>=n||membership[x])throw std::runtime_error("invalid subset");membership[x]=true;}
        std::vector<int>outside,residual(m);for(int i=0;i<n;++i)if(!membership[i])outside.push_back(i);
        std::vector<uint64_t>signature(n);for(int i=0;i<n;++i)for(int q=0;q<m;++q)signature[i]|=uint64_t(bits[i*m+q])<<q;
        uint64_t wanted=0;int initial_energy=0;for(int q=0;q<m;++q){residual[q]=-target[q];for(int x:selected)residual[q]+=bits[x*m+q];if(residual[q]&1)throw std::runtime_error("input has not passed mod2");initial_energy+=residual[q]*residual[q];wanted|=uint64_t((residual[q]/2)&1)<<q;}
        Table triples(size_t(quota)*(quota-1)*(quota-2)/6);bool complete=true;uint64_t scans=0;
        for(int a=0;a<quota&&complete;++a)for(int b=a+1;b<quota&&complete;++b)for(int c=b+1;c<quota;++c){triples.insert(signature[selected[a]]^signature[selected[b]]^signature[selected[c]],selected[a]|(selected[b]<<10)|(selected[c]<<20));if((triples.entries.size()&16383)==0&&elapsed()>.72*seconds){complete=false;break;}}
        std::vector<Move>moves;moves.reserve(50000);
        for(size_t a=0;a<outside.size()&&complete;++a)for(size_t b=a+1;b<outside.size()&&complete;++b)for(size_t c=b+1;c<outside.size();++c){
            uint64_t key=signature[outside[a]]^signature[outside[b]]^signature[outside[c]];size_t bucket=triples.bucket(key);
            for(int e=triples.heads[bucket];e!=-1;e=triples.entries[e].next){Move move;uint32_t code=triples.entries[e].code;for(int k=0;k<3;++k)move.out[k]=(code>>(10*k))&1023;move.in[0]=outside[a];move.in[1]=outside[b];move.in[2]=outside[c];
                for(int q=0;q<m;++q){int d=0;for(int k=0;k<3;++k)d+=bits[move.in[k]*m+q]-bits[move.out[k]*m+q];if(d&1)throw std::runtime_error("kernel syndrome mismatch");if(d>0)move.positive|=uint64_t(1)<<q;if(d<0)move.negative|=uint64_t(1)<<q;}
                moves.push_back(move);
            }
            if((++scans&16383)==0&&elapsed()>.72*seconds){complete=false;break;}
        }
        std::unordered_multimap<uint64_t,uint32_t>dictionary;dictionary.reserve(2*moves.size()+1);for(uint32_t i=0;i<moves.size();++i)dictionary.emplace(moves[i].syndrome(),i);
        std::ofstream catalogue(std::string(argv[1])+".kernel.tsv");
        if(!catalogue)throw std::runtime_error("cannot save kernel catalogue");
        catalogue<<"MOD4CAT1 "<<n<<' '<<m<<' '<<quota<<' '<<moves.size()<<' '<<complete<<'\n';
        for(const Move&move:moves){for(int k=0;k<3;++k)catalogue<<move.out[k]<<' ';for(int k=0;k<3;++k)catalogue<<move.in[k]<<' ';catalogue<<move.positive<<' '<<move.negative<<'\n';}catalogue.close();if(!catalogue)throw std::runtime_error("kernel catalogue write failed");
        uint64_t basis[64]={0};int rank=0;for(const Move&move:moves){uint64_t v=move.syndrome();for(int b=62;b>=0;--b)if(v>>b&1){if(basis[b])v^=basis[b];else{basis[b]=v;++rank;break;}}}
        uint64_t remainder=wanted;for(int b=62;b>=0;--b)if(remainder>>b&1)remainder^=basis[b];
        int best_energy=wanted==0?initial_energy:std::numeric_limits<int>::max();std::vector<int>best_selected=selected,best_residual=residual;uint64_t samples=0,disjoint_pairs=0,macros=0;int best_count=0;
        std::vector<Move>best_macro;
        auto inspect=[&](const std::vector<uint32_t>&indices){
            for(size_t a=0;a<indices.size();++a)for(size_t b=0;b<a;++b)if(!disjoint(moves[indices[a]],moves[indices[b]]))return;
            std::vector<int>r=residual;for(int q=0;q<m;++q)for(uint32_t i:indices)r[q]+=2*int(moves[i].positive>>q&1)-2*int(moves[i].negative>>q&1);
            int energy=0;for(int x:r){if(x%4)throw std::runtime_error("mod4 macro mismatch");energy+=x*x;}++macros;
            if(energy>=best_energy)return;std::vector<bool>yes=membership;for(uint32_t i:indices)for(int k=0;k<3;++k){yes[moves[i].out[k]]=false;yes[moves[i].in[k]]=true;}
            best_selected.clear();for(int i=0;i<n;++i)if(yes[i])best_selected.push_back(i);if(int(best_selected.size())!=quota)throw std::runtime_error("macro quota lost");best_residual=r;best_energy=energy;best_count=indices.size();best_macro.clear();for(uint32_t i:indices)best_macro.push_back(moves[i]);
        };
        auto one=dictionary.equal_range(wanted);for(auto it=one.first;it!=one.second;++it)inspect({it->second});
        for(uint32_t i=0;i<moves.size()&&best_energy;++i){auto range=dictionary.equal_range(wanted^moves[i].syndrome());for(auto it=range.first;it!=range.second;++it)if(it->second!=i)inspect({i,it->second});if((i&1023)==0&&elapsed()>=seconds)break;}
        while(!moves.empty()&&samples<1000000&&best_energy&&elapsed()<seconds){
            ++samples;uint32_t a=rng()%moves.size(),b=rng()%moves.size();if(a==b||!disjoint(moves[a],moves[b]))continue;++disjoint_pairs;
            auto range=dictionary.equal_range(wanted^moves[a].syndrome()^moves[b].syndrome());for(auto it=range.first;it!=range.second;++it)if(it->second!=a&&it->second!=b)inspect({a,b,it->second});
        }
        bool lifted=best_energy!=std::numeric_limits<int>::max();if(!lifted){best_energy=initial_energy;best_selected=selected;best_residual=residual;}
        std::vector<int>checked(m);for(int q=0;q<m;++q){checked[q]=-target[q];for(int x:best_selected)checked[q]+=bits[x*m+q];}if(checked!=best_residual)throw std::runtime_error("independent macro replay failed");
        int modulus=2;while(modulus<1024&&std::all_of(checked.begin(),checked.end(),[&](int x){return x%(2*modulus)==0;}))modulus*=2;
        std::cout<<"{\"schema\":\"sextic-row-mod4-v1\",\"seed\":"<<seed<<",\"cpu_seconds\":"<<elapsed()<<",\"initial_energy\":"<<initial_energy<<",\"initial_mod2_sector\":true,\"best_energy\":"<<best_energy<<",\"all_mod2_passed\":true,\"passed_modulus\":"<<modulus<<",\"mod4_lift_found\":"<<(lifted?"true":"false")<<",\"catalogue_complete\":"<<(complete?"true":"false")<<",\"kernel_moves\":"<<moves.size()<<",\"half_syndrome_rank\":"<<rank<<",\"target_in_catalogue_syndrome_span\":"<<(remainder==0?"true":"false")<<",\"sampled_pairs\":"<<samples<<",\"disjoint_pairs\":"<<disjoint_pairs<<",\"verified_macro_matches\":"<<macros<<",\"best_macro_moves\":"<<best_count<<",\"selected\":[";
        for(size_t i=0;i<best_selected.size();++i){if(i)std::cout<<',';std::cout<<best_selected[i];}std::cout<<"],\"residual\":[";for(int q=0;q<m;++q){if(q)std::cout<<',';std::cout<<checked[q];}
        std::cout<<"],\"chosen_macro\":[";for(size_t j=0;j<best_macro.size();++j){if(j)std::cout<<',';const Move&move=best_macro[j];std::cout<<"{\"out\":[";for(int k=0;k<3;++k){if(k)std::cout<<',';std::cout<<move.out[k];}std::cout<<"],\"in\":[";for(int k=0;k<3;++k){if(k)std::cout<<',';std::cout<<move.in[k];}std::cout<<"],\"positive_half_mask\":"<<move.positive<<",\"negative_half_mask\":"<<move.negative<<'}';}std::cout<<"]}\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
