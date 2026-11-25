// Stream all solutions of a Boolean quadratic system in 64-assignment words.
// Only a short prefix of equation words is updated at every Gray-code step;
// the remaining words are evaluated on demand. Memory is polynomial in d,m.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct Equation {
    uint64_t base=0, linearHigh=0;
    std::vector<uint64_t> highRows;
    std::vector<uint8_t> crossLow;
};

static inline unsigned parity(uint64_t x) { return __builtin_parityll(x); }

int main(int argc,char**argv) {
    if(argc<4) return 2;
    try {
        std::ifstream in(argv[1]);
        unsigned d=0,m=0;
        in>>d>>m;
        if(!in||d>60||m==0) throw std::runtime_error("invalid dimensions");
        const unsigned low=std::min(6u,d),high=d-low;
        const uint64_t valid=low==6?~0ULL:((1ULL<<(1u<<low))-1);
        const uint64_t lowMask=(1ULL<<low)-1;
        std::vector<uint64_t> bits(low),linearTruth(1u<<low);
        for(unsigned u=0;u<low;u++)
            for(unsigned z=0;z<(1u<<low);z++) if((z>>u)&1) bits[u]|=1ULL<<z;
        for(unsigned mask=1;mask<(1u<<low);mask++) {
            unsigned u=__builtin_ctz(mask);
            linearTruth[mask]=linearTruth[mask^(1u<<u)]^bits[u];
        }
        std::vector<Equation> equations(m);
        for(unsigned ell=0;ell<m;ell++) {
            unsigned constant;uint64_t linear;
            in>>constant>>linear;
            if(!in||constant>1||(d<64&&(linear>>d))) throw std::runtime_error("invalid polynomial");
            std::vector<uint64_t> rows(d);
            for(unsigned u=0;u<d;u++) {
                in>>rows[u];
                if(!in||(rows[u]>>u)) throw std::runtime_error("quadratic rows must be strictly lower triangular");
            }
            auto &e=equations[ell];
            e.base=(constant?valid:0)^linearTruth[linear&lowMask];
            e.linearHigh=linear>>low;
            e.highRows.assign(high,0);e.crossLow.assign(high,0);
            for(unsigned u=0;u<low;u++)
                for(unsigned v=0;v<u;v++) if((rows[u]>>v)&1) e.base^=bits[u]&bits[v];
            for(unsigned u=0;u<high;u++) {
                e.crossLow[u]=rows[low+u]&lowMask;
                for(unsigned v=0;v<u;v++) if((rows[low+u]>>(low+v))&1) {
                    e.highRows[u]|=1ULL<<v;
                    e.highRows[v]|=1ULL<<u;
                }
            }
        }
        const unsigned tracked=std::min(8u,m);
        std::vector<uint64_t> words(tracked),derivative(tracked);
        std::vector<std::vector<uint64_t>> changeWord(high,std::vector<uint64_t>(tracked)),
            changeDerivative(high,std::vector<uint64_t>(tracked));
        for(unsigned ell=0;ell<tracked;ell++) {
            words[ell]=equations[ell].base;
            derivative[ell]=equations[ell].linearHigh;
            for(unsigned u=0;u<high;u++) {
                changeWord[u][ell]=linearTruth[equations[ell].crossLow[u]];
                changeDerivative[u][ell]=equations[ell].highRows[u];
            }
        }
        auto evaluate=[&](unsigned ell,uint64_t h) {
            auto &e=equations[ell];
            uint8_t cross=0;
            unsigned value=parity(e.linearHigh&h);
            for(uint64_t q=h;q;q&=q-1) {
                unsigned u=__builtin_ctzll(q);
                cross^=e.crossLow[u];
                value^=parity(e.highRows[u]&h&((1ULL<<u)-1));
            }
            return e.base^linearTruth[cross]^(value?valid:0);
        };
        const double seconds=std::stod(argv[3]);
        const uint64_t rootLimit=argc>4?std::stoull(argv[4]):100000;
        if(seconds<=0) throw std::runtime_error("positive time limit required");
        auto started=std::chrono::steady_clock::now();
        const uint64_t blocks=1ULL<<high;
        std::vector<uint64_t> roots;
        uint64_t completed=0;
        bool timeLimit=false,countLimit=false;
        for(uint64_t step=0;step<blocks;step++) {
            if((step&16383)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>=seconds) {
                timeLimit=true;break;
            }
            uint64_t h=step^(step>>1);
            if(step) {
                unsigned u=__builtin_ctzll(step);
                const auto &cw=changeWord[u],&cd=changeDerivative[u];
                for(unsigned ell=0;ell<tracked;ell++) {
                    words[ell]^=cw[ell]^(-uint64_t((derivative[ell]>>u)&1));
                    derivative[ell]^=cd[ell];
                }
            }
            uint64_t good=valid;
            for(unsigned ell=0;ell<tracked&&good;ell++) good&=~words[ell];
            for(unsigned ell=tracked;ell<m&&good;ell++) good&=~evaluate(ell,h);
            completed++;
            while(good) {
                unsigned z=__builtin_ctzll(good);good&=good-1;
                roots.push_back(uint64_t(z)|(h<<low));
                if(rootLimit&&roots.size()>=rootLimit) {countLimit=true;break;}
            }
            if(countLimit) break;
        }
        double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        bool complete=!timeLimit&&!countLimit&&completed==blocks;
        std::ofstream out(argv[2]);
        out<<"{\"variables\":"<<d<<",\"equations\":"<<m
           <<",\"complete\":"<<(complete?"true":"false")
           <<",\"time_limit_reached\":"<<(timeLimit?"true":"false")
           <<",\"root_limit_reached\":"<<(countLimit?"true":"false")
           <<",\"assignments_covered\":"<<(completed<<low)
           <<",\"elapsed_seconds\":"<<elapsed<<",\"roots\":[";
        for(size_t i=0;i<roots.size();i++) {if(i)out<<",";out<<roots[i];}
        out<<"]}\n";
        std::cout<<"{\"complete\":"<<(complete?"true":"false")
                 <<",\"assignments_covered\":"<<(completed<<low)
                 <<",\"roots\":"<<roots.size()<<",\"seconds\":"<<elapsed<<"}\n";
        return 0;
    } catch(const std::exception&e) {std::cerr<<e.what()<<"\n";return 1;}
}
