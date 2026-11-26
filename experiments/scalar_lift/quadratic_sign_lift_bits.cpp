// Exact exhaustive Boolean-MQ filter. Polynomial values on small blocks are
// bit sliced. Across equations, high-variable values and low-variable linear
// coefficients are packed into words. Reflected Gray order changes a variable's
// derivative by one Hessian coefficient between consecutive uses.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

template<unsigned Low> int scan(std::istream& in,unsigned d,unsigned m,char**argv,int argc) {
    constexpr unsigned Lanes=Low>6?1u<<(Low-6):1;
    constexpr unsigned PerWord=Low?64/Low:64;
    constexpr uint64_t LowMask=(1ULL<<Low)-1;
    constexpr uint64_t Valid=Low>=6?~0ULL:((1ULL<<(1u<<Low))-1);
    using Word=std::array<uint64_t,Lanes>;
    const unsigned high=d-Low,crossWords=(m+PerWord-1)/PerWord,constantWords=(m+63)/64;
    std::array<Word,Low> bits{};
    std::vector<Word> linearTruth(1u<<Low),base(m);
    for(unsigned u=0;u<Low;u++)
        for(unsigned z=0;z<(1u<<Low);z++) if((z>>u)&1) bits[u][z/64]|=1ULL<<(z%64);
    for(unsigned mask=1;mask<(1u<<Low);mask++) {
        unsigned u=__builtin_ctz(mask);
        for(unsigned w=0;w<Lanes;w++) linearTruth[mask][w]=linearTruth[mask^(1u<<u)][w]^bits[u][w];
    }
    std::vector<uint64_t> cross(high*crossWords),derivative(high*constantWords),
        hessian(high*high*constantWords);
    for(unsigned ell=0;ell<m;ell++) {
        unsigned constant;uint64_t linear;
        in>>constant>>linear;
        if(!in||constant>1||(linear>>d)) throw std::runtime_error("invalid polynomial");
        std::vector<uint64_t> rows(d);
        for(unsigned u=0;u<d;u++) {
            in>>rows[u];
            if(!in||(rows[u]>>u)) throw std::runtime_error("quadratic rows must be strictly lower triangular");
        }
        for(unsigned w=0;w<Lanes;w++) base[ell][w]=(constant?Valid:0)^linearTruth[linear&LowMask][w];
        for(unsigned u=0;u<Low;u++)
            for(unsigned v=0;v<u;v++) if((rows[u]>>v)&1)
                for(unsigned w=0;w<Lanes;w++) base[ell][w]^=bits[u][w]&bits[v][w];
        const unsigned cw=ell/64;
        const uint64_t bit=1ULL<<(ell%64);
        for(unsigned u=0;u<high;u++) {
            cross[u*crossWords+ell/PerWord]|=(rows[Low+u]&LowMask)<<((ell%PerWord)*Low);
            // At u's first flip the sole set lower Gray bit is u-1.
            if(((linear>>(Low+u))&1) ^ (u?((rows[Low+u]>>(Low+u-1))&1):0))
                derivative[u*constantWords+cw]|=bit;
            for(unsigned v=0;v<u;v++) if((rows[Low+u]>>(Low+v))&1) {
                hessian[(u*high+v)*constantWords+cw]|=bit;
                hessian[(v*high+u)*constantWords+cw]|=bit;
            }
        }
    }
    in>>std::ws;
    if(!in.eof()) throw std::runtime_error("unexpected trailing polynomial data");
    const double seconds=std::stod(argv[3]);
    const uint64_t rootLimit=argc>4?std::stoull(argv[4]):100000;
    if(seconds<=0) throw std::runtime_error("positive time limit required");
    auto started=std::chrono::steady_clock::now();
    const uint64_t blocks=1ULL<<high;
    std::vector<uint64_t> roots,coefficients(crossWords),constants(constantWords);
    uint64_t completed=0;
    bool timeLimit=false,countLimit=false;
    for(uint64_t step=0;step<blocks;step++) {
        if((step&16383)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>=seconds) {
            timeLimit=true;break;
        }
        uint64_t h=step^(step>>1);
        if(step) {
            unsigned u=__builtin_ctzll(step);
            for(unsigned w=0;w<crossWords;w++) coefficients[w]^=cross[u*crossWords+w];
            for(unsigned w=0;w<constantWords;w++) constants[w]^=derivative[u*constantWords+w];
            // Between consecutive flips of u, exactly one higher bit v
            // changes, at step+2^u. All lower Gray bits repeat.
            unsigned v=__builtin_ctzll(step+(1ULL<<u));
            if(v<high) for(unsigned w=0;w<constantWords;w++)
                derivative[u*constantWords+w]^=hessian[(u*high+v)*constantWords+w];
        }
        Word good;good.fill(Valid);
        uint64_t any=Valid;
        unsigned ell=0;
        for(unsigned w=0;w<crossWords&&any;w++) {
            uint64_t packed=coefficients[w];
            unsigned stop=std::min(m,ell+PerWord);
            for(;ell<stop&&any;ell++,packed>>=Low) {
                const Word &truth=linearTruth[packed&LowMask],&initial=base[ell];
                uint64_t sign=-((constants[ell/64]>>(ell%64))&1ULL);
                any=0;
                for(unsigned lane=0;lane<Lanes;lane++) {
                    good[lane]&=~(initial[lane]^truth[lane]^sign);
                    any|=good[lane];
                }
            }
        }
        completed++;
        for(unsigned lane=0;lane<Lanes&&!countLimit;lane++) while(good[lane]) {
            unsigned z=__builtin_ctzll(good[lane]);good[lane]&=good[lane]-1;
            roots.push_back(uint64_t(64*lane+z)|(h<<Low));
            if(rootLimit&&roots.size()>=rootLimit) {countLimit=true;break;}
        }
        if(countLimit) break;
    }
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    bool complete=!timeLimit&&!countLimit&&completed==blocks;
    std::ofstream out(argv[2]);
    out<<"{\"variables\":"<<d<<",\"equations\":"<<m<<",\"block_variables\":"<<Low
       <<",\"algorithm\":\"packed-gray-bit-sliced-v2\""
       <<",\"complete\":"<<(complete?"true":"false")
       <<",\"time_limit_reached\":"<<(timeLimit?"true":"false")
       <<",\"root_limit_reached\":"<<(countLimit?"true":"false")
       <<",\"assignments_covered\":"<<(completed<<Low)
       <<",\"elapsed_seconds\":"<<elapsed<<",\"roots\":[";
    for(size_t i=0;i<roots.size();i++) {if(i)out<<",";out<<roots[i];}
    out<<"]}\n";
    std::cout<<"{\"complete\":"<<(complete?"true":"false")
             <<",\"assignments_covered\":"<<(completed<<Low)
             <<",\"roots\":"<<roots.size()<<",\"seconds\":"<<elapsed<<"}\n";
    return 0;
}

int main(int argc,char**argv) {
    if(argc<4) return 2;
    try {
        std::ifstream in(argv[1]);
        unsigned d=0,m=0;
        in>>d>>m;
        if(!in||d>60||m==0) throw std::runtime_error("invalid dimensions");
        unsigned low=argc>5?std::stoul(argv[5]):11;
        if(low<6||low>12) throw std::runtime_error("block variables must be between 6 and 12");
        low=std::min(low,d);
        switch(low) {
            case 0:return scan<0>(in,d,m,argv,argc);
            case 1:return scan<1>(in,d,m,argv,argc);
            case 2:return scan<2>(in,d,m,argv,argc);
            case 3:return scan<3>(in,d,m,argv,argc);
            case 4:return scan<4>(in,d,m,argv,argc);
            case 5:return scan<5>(in,d,m,argv,argc);
            case 6:return scan<6>(in,d,m,argv,argc);
            case 7:return scan<7>(in,d,m,argv,argc);
            case 8:return scan<8>(in,d,m,argv,argc);
            case 9:return scan<9>(in,d,m,argv,argc);
            case 10:return scan<10>(in,d,m,argv,argc);
            case 11:return scan<11>(in,d,m,argv,argc);
            case 12:return scan<12>(in,d,m,argv,argc);
        }
        return 2;
    } catch(const std::exception&e) {std::cerr<<e.what()<<"\n";return 1;}
}
