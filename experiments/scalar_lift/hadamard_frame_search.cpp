// Active exact twisted-T frame search; see docs/HADAMARD_FRAME_SEARCH.md.
// Extracted from the successful H7572 kernel, without research-checkpoint imports.
// Modes: search (SL2 parity axes for p631) and lift (generic finite axis pools).
// Original optimized sign kernel SHA256: 28d3ce519e1ac9066b070b835f6a82badbd1d300e298c15f4e1e9d9f736b1d87
// Original optimized sampler SHA256: d5b1b14202e52e39fddb5cf81d91cdd3a271d4bdf4f3230beed62afb93de09aa
// Exact mod4/mod8 sign lifting for a fixed twisted one-hot axis word.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace twisted_sign_lift {
    inline int parity(uint64_t x){
        return __builtin_popcountll(x)&1;
    } struct Config{
        int headNet=1,tailNet=1,firstNegative=0,firstAxis=0,maxNullity=24;
        bool mod8=true,stopFirst=true,validateMod8=false,collectMod8=false;
    };
    struct Result{
        std::string status;
        int rank=0,nullity=0;
        uint64_t affineAssignments=0,mod8Passed=0,DCpassed=0;
        std::vector<uint64_t>solutions,mod8Masks;
    };
    class Engine{
        public: int p,m,mu,n,L;
        std::vector<int16_t>H,K;
        std::vector<uint8_t>axisWeights;
        mutable std::array<std::vector<std::array<uint64_t,2>>,17>outerProducts;
        std::vector<int>index,phase;
        explicit Engine(const std::string&path){
            std::ifstream in(path,std::ios::binary);
            int32_t h[7];
            in.read((char*)h,sizeof(h));
            p=h[0];
            m=h[1];
            mu=h[2];
            n=h[3];
            L=h[4];
            if(!in||n<2||n>63||L<=0)throw std::runtime_error("bad T coefficient table");
            H.resize(n*8*L);
            K.resize((size_t)n*n*64*L);
            in.read((char*)H.data(),H.size()*2);
            in.read((char*)K.data(),K.size()*2);
            if(!in)throw std::runtime_error("truncated T table");
            index.resize(p);
            phase.resize(p);
            std::vector<int>seen(p);
            seen[0]=1;
            int orbit=0;
            for(int x=1;x<p;x++)if(!seen[x]){
                int y=x;
                orbit++;
                for(int k=0;k<m;k++){
                    if(seen[y])throw std::runtime_error("bad orbit");
                    seen[y]=1;
                    index[y]=orbit;
                    phase[y]=k%3;
                    y=y*mu%p;
                }if(y!=x)throw std::runtime_error("bad multiplier");
            }if(orbit+1!=n)throw std::runtime_error("bad orbit count");
            axisWeights.resize((size_t)n*(n-1)/2*16*L);
            size_t e=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,e++)for(int a=0;a<4;a++)for(int b=0;b<4;b++)for(int l=0;l<L;l++){
                int w=K[((((size_t)i*n+j)*8+2*a)*8+2*b)*L+l];
                if(w<0||w>255)throw std::runtime_error("unsupported unsigned coefficient");
                axisWeights[(e*16+4*a+b)*L+l]=w;
            }
        } bool verifyFull(const std::vector<int>&a,uint64_t x)const{
            std::vector<int>b(4*p);
            for(int r=0;r<p;r++){
                int i=index[r],v=a[i],at=v?1+(v-1+phase[r])%3:0;
                b[4*r+at]=(x>>i&1)?-1:1;
            }for(int s=0;s<p;s++){
                int value=0;
                for(int r=0;r<p;r++)for(int c=0;c<4;c++)value+=b[4*r+c]*b[4*((r+s)%p)+c];
                if(value!=(s==0?p:0))return false;
            }return true;
        } Result lift(const std::vector<int>&axes,const Config&cfg)const{
            Result out;
            if(int(axes.size())!=n||axes[0]||axes[1]!=cfg.firstAxis)throw std::runtime_error("bad axis gauge");
            for(int a:axes)if(a<0||a>3)throw std::runtime_error("bad axis");
            uint64_t head=0,tail=0;
            for(int i=1;i<n;i++)(axes[i]?tail:head)|=1ULL<<i;
            int kh=__builtin_popcountll(head),kt=n-1-kh;
            if((kh-cfg.headNet)%2||(kt-cfg.tailNet)%2){
                out.status="DC parity";
                return out;
            }int hn=(kh-cfg.headNet)/2,tn=(kt-cfg.tailNet)/2;
            if(hn<0||hn>kh||tn<0||tn>kt){
                out.status="DC range";
                return out;
            } std::vector<int>cabs(L);
            std::vector<uint64_t>degree(L);
            std::vector<uint16_t>weights((size_t)n*(n-1)/2*L);
            for(int i=0;i<n;i++)for(int l=0;l<L;l++)cabs[l]+=H[(i*8+2*axes[i])*L+l];
            size_t edge=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,edge+=L)for(int l=0;l<L;l++){
                unsigned w=axisWeights[edge*16+(4*axes[i]+axes[j])*L+l];
                weights[edge+l]=w;
                cabs[l]+=w;
                degree[l]^=((1ULL<<i)|(1ULL<<j))&(-uint64_t(w&1));
            } std::vector<uint64_t>rows;
            for(int l=0;l<L;l++){
                if(cabs[l]%2){
                    out.status="odd unsigned correlation";
                    return out;
                }rows.push_back(degree[l]|(uint64_t((cabs[l]/2)&1)<<n));
            }rows.push_back(1);
            rows.push_back(2|(uint64_t(cfg.firstNegative)<<n));
            rows.push_back(head|(uint64_t(hn&1)<<n));
            rows.push_back(tail|(uint64_t(tn&1)<<n));
            std::vector<int>pivot;
            int at=0;
            for(int col=0;col<n;col++){
                int found=-1;
                for(int i=at;i<int(rows.size());i++)if(rows[i]>>col&1){
                    found=i;
                    break;
                }if(found<0)continue;
                std::swap(rows[at],rows[found]);
                for(int i=0;i<int(rows.size());i++)if(i!=at&&(rows[i]>>col&1))rows[i]^=rows[at];
                pivot.push_back(col);
                at++;
            } uint64_t all=(1ULL<<n)-1;
            for(uint64_t r:rows)if(!(r&all)&&(r>>n&1)){
                out.status="inconsistent mod4";
                return out;
            }out.rank=at;
            uint64_t particular=0;
            for(int i=0;i<at;i++)if(rows[i]>>n&1)particular|=1ULL<<pivot[i];
            std::vector<uint64_t>basis;
            for(int j=0;j<n;j++)if(std::find(pivot.begin(),pivot.end(),j)==pivot.end()){
                uint64_t v=1ULL<<j;
                for(int i=0;i<at;i++)if(rows[i]>>j&1)v|=1ULL<<pivot[i];
                basis.push_back(v);
            }int d=basis.size();
            out.nullity=d;
            if(d>cfg.maxNullity){
                out.status="nullity limit";
                return out;
            }out.affineAssignments=1ULL<<d;
            auto correlation=[&](uint64_t x,int l){
                int cut=0;
                size_t edge=0;
                for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,edge+=L)if(((x>>i)^(x>>j))&1)cut+=weights[edge+l];
                return cabs[l]-2*cut;
            };
            auto accept=[&](uint64_t x){
                if(__builtin_popcountll(x&head)!=hn||__builtin_popcountll(x&tail)!=tn)return false;
                out.DCpassed++;
                for(int l=0;l<L;l++)if(correlation(x,l))return false;
                if(!verifyFull(axes,x))throw std::runtime_error("compressed zero failed expanded verification");
                out.solutions.push_back(x);
                return cfg.stopFirst;
            };
            if(!cfg.mod8||d>24){
                uint64_t x=particular;
                for(uint64_t k=0;k<(1ULL<<d);k++){
                    if(k)x^=basis[__builtin_ctzll(k)];
                    if(accept(x))break;
                }out.status="enumerated mod4";
                return out;
            }  // Compute all d+1 singleton affine values together in bitplanes of cut mod4.
            // Coordinate d represents particular; coordinate u represents particular xor basis[u].
            uint64_t affineBits=(1ULL<<(d+1))-1,kernelBits=(1ULL<<d)-1;
            std::vector<uint64_t>vertexBasis(n),cutLow(L),cutHigh(L),packedQ0(L),packedQ1(L),laplacianBasis(d>16?L*n:0);
            if(d<=16&&outerProducts[d].empty()){
                auto&table=outerProducts[d];
                table.resize(1ULL<<d);
                for(uint64_t mask=0;mask<(1ULL<<d);mask++){
                    int bit=0;
                    for(int u=0;u<d;u++)for(int v=0;v<u;v++,bit++)if((mask>>u&1)&&(mask>>v&1))table[mask][bit/64]|=1ULL<<(bit%64);
                }
            } for(int u=0;u<d;u++)for(uint64_t q=basis[u];q;q&=q-1)vertexBasis[__builtin_ctzll(q)]|=1ULL<<u;
            edge=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,edge+=L){
                uint64_t difference=vertexBasis[i]^vertexBasis[j];
                uint64_t flip=difference^(((((particular>>i)^(particular>>j))&1)!=0)?affineBits:0);
                if(d<=16){
                    auto outer=outerProducts[d][difference];
                    for(int l=0;l<L;l++){
                        unsigned w=weights[edge+l];
                        uint64_t odd=-uint64_t(w&1),twice=-uint64_t((w>>1)&1);
                        cutHigh[l]^=(cutLow[l]&flip&odd)^(flip&twice);
                        cutLow[l]^=flip&odd;
                        packedQ0[l]^=outer[0]&odd;
                        packedQ1[l]^=outer[1]&odd;
                    }
                } else for(int l=0;l<L;l++){
                    unsigned w=weights[edge+l];
                    if(w&1){
                        cutHigh[l]^=cutLow[l]&flip;
                        cutLow[l]^=flip;
                        laplacianBasis[l*n+i]^=difference;
                        laplacianBasis[l*n+j]^=difference;
                    }if(w&2)cutHigh[l]^=flip;
                }
            } std::vector<int>constant(L);
            std::vector<uint64_t>linear(L);
            std::vector<std::vector<uint64_t>>quadratic(L,std::vector<uint64_t>(d));
            for(int l=0;l<L;l++){
                if(cutLow[l]!=((cabs[l]/2&1)?affineBits:0))throw std::runtime_error("bad bitplane mod4 affine values");
                uint64_t singleton=cutHigh[l]^((cabs[l]/4&1)?affineBits:0);
                constant[l]=(singleton>>d)&1;
                linear[l]=(singleton^(constant[l]?affineBits:0))&kernelBits;
                if(d<=16){
                    int bit=0;
                    for(int u=0;u<d;u++)for(int v=0;v<u;v++,bit++)if(((bit<64?packedQ0[l]:packedQ1[l])>>(bit%64))&1){
                        quadratic[l][u]|=1ULL<<v;
                        quadratic[l][v]|=1ULL<<u;
                    }
                } else for(int u=0;u<d;u++){
                    uint64_t row=0;
                    for(uint64_t q=basis[u];q;q&=q-1)row^=laplacianBasis[l*n+__builtin_ctzll(q)];
                    if(row>>u&1)throw std::runtime_error("nonalternating bitplane kernel Laplacian");
                    quadratic[l][u]=row;
                }
            } if(cfg.validateMod8){
                uint64_t cases=std::min<uint64_t>(1ULL<<d,128);
                for(uint64_t k=0;k<cases;k++){
                    uint64_t z=(k*0x9e3779b97f4a7c15ULL)&((1ULL<<d)-1),x=particular;
                    for(uint64_t q=z;q;q&=q-1)x^=basis[__builtin_ctzll(q)];
                    for(int l=0;l<L;l++){
                        int predicted=constant[l]^parity(linear[l]&z);
                        for(int u=0;u<d;u++)if(z>>u&1)predicted^=parity(quadratic[l][u]&z&((1ULL<<u)-1));
                        int value=correlation(x,l);
                        if(value%4||predicted!=((value/4)&1))throw std::runtime_error("quadratic mod8 disagrees with exact correlation");
                    }
                }
            } int low=std::min(6,d),high=d-low,blocks=1<<high;
            uint64_t lowValid=low==6?~0ULL:((1ULL<<(1<<low))-1);
            std::vector<uint64_t>Z(low),linearTruth(1<<low);
            for(int u=0;u<low;u++)for(int z=0;z<(1<<low);z++)if(z>>u&1)Z[u]|=1ULL<<z;
            for(int z=1;z<(1<<low);z++){
                int u=__builtin_ctz(z);
                linearTruth[z]=linearTruth[z^(1<<u)]^Z[u];
            } if(high<=18){
                int truthWords=(blocks+63)/64;
                int split=high/2,as=1<<split,bs=1<<(high-split);
                std::vector<uint64_t>lowBase(L);
                std::vector<uint64_t>highTruth(L*truthWords);
                std::vector<uint8_t>crossA(L*as),crossB(L*bs);
                for(int l=0;l<L;l++){
                    uint64_t word=(constant[l]?lowValid:0)^linearTruth[linear[l]&((1ULL<<low)-1)];
                    for(int u=0;u<low;u++)for(int v=0;v<u;v++)if(quadratic[l][u]>>v&1)word^=Z[u]&Z[v];
                    lowBase[l]=word;
                    // A high-variable truth word has six local variables;
                    // upper variables alter only its linear mask and constant.
                    int smallHigh=std::min(6,high);
                    uint64_t baseHigh=linearTruth[(linear[l]>>low)&((1ULL<<smallHigh)-1)];
                    for(int u=0;u<smallHigh;u++)for(int t=0;t<u;t++)
                        if(quadratic[l][low+u]>>(low+t)&1)baseHigh^=Z[u]&Z[t];
                    std::vector<uint8_t>upperCross(truthWords),upperValue(truthWords);
                    highTruth[l*truthWords]=baseHigh;
                    for(int word=1;word<truthWords;word++){
                        int u=__builtin_ctz(word),prev=word^(1<<u),v=low+6+u;
                        upperCross[word]=upperCross[prev]^uint8_t((quadratic[l][v]>>low)&63);
                        upperValue[word]=upperValue[prev]^int(linear[l]>>v&1)
                            ^parity((quadratic[l][v]>>(low+6))&uint64_t(prev));
                        highTruth[l*truthWords+word]=baseHigh^linearTruth[upperCross[word]]
                            ^(upperValue[word]?~0ULL:0);
                    }
                    for(int h=1;h<as;h++){
                        int u=__builtin_ctz(h),prev=h^(1<<u);
                        crossA[l*as+h]=crossA[l*as+prev]^uint8_t(quadratic[l][low+u]&((1ULL<<low)-1));
                    } for(int h=1;h<bs;h++){
                        int u=__builtin_ctz(h),prev=h^(1<<u);
                        crossB[l*bs+h]=crossB[l*bs+prev]^uint8_t(quadratic[l][low+split+u]&((1ULL<<low)-1));
                    }
                } bool done=false;
                for(int h=0;h<blocks&&!done;h++){
                    uint64_t good=lowValid;
                    for(int l=0;l<L&&good;l++){
                        int cross=crossA[l*as+(h&(as-1))]^crossB[l*bs+(h>>split)];
                        good&=~(lowBase[l]^linearTruth[cross]^((highTruth[l*truthWords+h/64]>>(h%64)&1)?lowValid:0));
                    }while(good){
                        int lowWord=__builtin_ctzll(good);
                        good&=good-1;
                        uint64_t z=uint64_t(lowWord)|(uint64_t(h)<<low),x=particular;
                        for(uint64_t q=z;q;q&=q-1)x^=basis[__builtin_ctzll(q)];
                        out.mod8Passed++;
                        if(cfg.collectMod8)out.mod8Masks.push_back(x);
                        if(accept(x)){
                            done=true;
                            break;
                        }
                    }
                } out.status="enumerated mod8";
                return out;
            } std::vector<uint64_t>lowBase(L);
            std::vector<std::vector<uint8_t>>cross(L,std::vector<uint8_t>(blocks)),highValue(L,std::vector<uint8_t>(blocks));
            for(int l=0;l<L;l++){
                uint64_t word=constant[l]?lowValid:0;
                for(int u=0;u<low;u++){
                    if(linear[l]>>u&1)word^=Z[u];
                    for(int v=0;v<u;v++)if(quadratic[l][u]>>v&1)word^=Z[u]&Z[v];
                }lowBase[l]=word;
                for(int h=1;h<blocks;h++){
                    int u=__builtin_ctz(h),prev=h^(1<<u),v=low+u;
                    cross[l][h]=cross[l][prev]^uint8_t(quadratic[l][v]&((1ULL<<low)-1));
                    highValue[l][h]=highValue[l][prev]^int(linear[l]>>v&1)^parity((quadratic[l][v]>>low)&uint64_t(prev));
                }
            } bool done=false;
            for(int h=0;h<blocks&&!done;h++){
                uint64_t good=lowValid;
                for(int l=0;l<L&&good;l++)good&=~(lowBase[l]^linearTruth[cross[l][h]]^(highValue[l][h]?lowValid:0));
                while(good){
                    int lowWord=__builtin_ctzll(good);
                    good&=good-1;
                    uint64_t z=uint64_t(lowWord)|(uint64_t(h)<<low),x=particular;
                    for(uint64_t q=z;q;q&=q-1)x^=basis[__builtin_ctzll(q)];
                    out.mod8Passed++;
                    if(cfg.collectMod8)out.mod8Masks.push_back(x);
                    if(accept(x)){
                        done=true;
                        break;
                    }
                }
            } out.status="enumerated mod8";
            return out;
        }
    };
}  // namespace twisted_sign_lift
// Uniform SL(2,2^15) parity axes -> native GF2/mod8 sign lifts -> exact T631.
// field_sign_search T_TABLE FIELD_COLUMNS OUTPUT_JSON SECONDS SEED [SAMPLE_ONLY_COUNT=0] [HEAD_COUNTS_CSV=all]
#include <chrono>
#include <cstdio>
#include <iostream>
#include <random>
#include <sstream>
#include <unistd.h>
using ColorWord = unsigned __int128;
std::string color_decimal(ColorWord value){
    if(!value)return "0";
    std::string digits;
    while(value){digits.push_back(char('0'+value%10));value/=10;}
    std::reverse(digits.begin(),digits.end());
    return digits;
}
int search_main(int argc,char**argv){
    if(argc<6)return 2;
    twisted_sign_lift::Engine engine(argv[1]);
    std::ifstream in(argv[2],std::ios::binary);
    int32_t h[7];
    in.read((char*)h,sizeof(h));
    if(!in||h[0]!=0x46503135||(h[1]!=2&&h[1]!=3)||h[3]<1||h[3]>25||engine.n!=2*h[3]+1)return 4;
    if(h[1]==2&&h[3]>16)return 4;
    int degree=h[3];
    uint32_t size=1u<<degree,mask=size-1,polynomial=h[2];
    std::vector<ColorWord>column(4*degree);
    for(auto&value:column){
        uint64_t words[2]={0,0};
        in.read((char*)words,(h[1]==3?2:1)*sizeof(uint64_t));
        value=ColorWord(words[0])|(ColorWord(words[1])<<64);
    }
    if(!in)return 4;
    auto mul=[&](uint32_t a,uint32_t b){
        uint32_t c=0;
        while(b){
            if(b&1)c^=a;
            b>>=1;
            a<<=1;
            if(a&size)a^=polynomial;
        }return c;
    };
    auto power=[&](uint32_t a,uint32_t k){
        uint32_t c=1;
        while(k){
            if(k&1)c=mul(c,a);
            a=mul(a,a);
            k>>=1;
        }return c;
    };
    bool preloadInverse=degree<=17;
    std::vector<uint32_t>inverse(preloadInverse?size:0);
    for(uint32_t a=1;preloadInverse&&a<size;a++){
        inverse[a]=power(a,size-2);
        if(mul(a,inverse[a])!=1)return 5;
    }
    auto invert=[&](uint32_t value){
        if(!value)throw std::runtime_error("zero field inverse");
        if(preloadInverse)return inverse[value];
        uint32_t a=value,b=polynomial;
        uint64_t u=1,v=0;
        while(a!=1){
            if(!a)throw std::runtime_error("field modulus is reducible");
            int shift=(31-__builtin_clz(a))-(31-__builtin_clz(b));
            if(shift<0){std::swap(a,b);std::swap(u,v);shift=-shift;}
            a^=b<<shift;u^=v<<shift;
        }
        while(u>=size){int shift=(63-__builtin_clzll(u))-degree;u^=uint64_t(polynomial)<<shift;}
        uint32_t result=uint32_t(u);
        if(mul(value,result)!=1)throw std::runtime_error("field inverse recurrence failed");
        return result;
    };
    int colorBytes=(4*degree+7)/8;
    std::vector<std::array<ColorWord,256>>lookup(colorBytes);
    for(int byte=0;byte<colorBytes;byte++)for(int v=1;v<256;v++){
        int bit=__builtin_ctz(v);
        lookup[byte][v]=lookup[byte][v^(1<<bit)]^((8*byte+bit<4*degree)?column[8*byte+bit]:0);
    } double seconds=atof(argv[4]);
    uint64_t seed=strtoull(argv[5],nullptr,10),sampleOnly=argc>6?strtoull(argv[6],nullptr,10):0;
    std::mt19937_64 rng(seed);
    std::string output=argv[3];
    twisted_sign_lift::Config cfg;
    cfg.mod8=true;
    cfg.stopFirst=true;
    cfg.firstNegative=h[4];
    cfg.headNet=h[5];
    cfg.tailNet=h[6];
    std::vector<bool>allowHead(engine.n,true);
    std::vector<int>headFilter;
    if(argc>7){
        std::fill(allowHead.begin(),allowHead.end(),false);
        std::stringstream csv(argv[7]);
        std::string token;
        while(std::getline(csv,token,',')){
            size_t used=0;
            int value=std::stoi(token,&used);
            if(used!=token.size()||value<0||value>=engine.n)throw std::runtime_error("invalid head filter");
            allowHead[value]=true;
            headFilter.push_back(value);
        }if(headFilter.empty())throw std::runtime_error("empty head filter");
    } uint64_t headRejected=0;
    uint64_t sampled=0,gauged=0,tested=0,affine=0,mod8=0,dc=0,inconsistent=0,rankLimit=0;
    std::vector<uint64_t>heads(engine.n);
    std::vector<std::string>examples;
    std::vector<int>witness;
    uint64_t witnessNegative=0;
    auto start=std::chrono::steady_clock::now();
    double lastout=0;
    auto report=[&](){
        std::string temporary=output+".tmp";
        std::ofstream out(temporary);
        out<<"{\"schema\":\"field-parity-mod8-sign-search-v1\",\"p\":"<<engine.p<<",\"m\":"<<engine.m<<",\"mu\":"<<engine.mu<<",\"field_degree\":"<<degree<<",\"pid\":"<<getpid()<<",\"seed\":"<<seed<<",\"field_parameters_sampled\":"<<sampled<<",\"head_gauge_passed\":"<<gauged<<",\"axis_patterns_sign_lifted\":"<<tested<<",\"affine_sign_assignments_covered\":"<<affine<<",\"mod8_passed\":"<<mod8<<",\"DC_passed_after_mod8\":"<<dc<<",\"mod4_inconsistent\":"<<inconsistent<<",\"rank_limit_skipped\":"<<rankLimit<<",\"head_filter_rejected\":"<<headRejected<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"head_count_histogram\":[";
        for(int i=0;i<engine.n;i++){
            if(i)out<<",";
            out<<heads[i];
        }out<<"],\"head_filter\":[";
        for(size_t i=0;i<headFilter.size();i++){
            if(i)out<<",";
            out<<headFilter[i];
        }out<<"],\"parameter_examples\":[";
        for(size_t i=0;i<examples.size();i++){
            if(i)out<<",";
            out<<examples[i];
        }out<<"]";
        if(!witness.empty()){
            out<<",\"states\":[";
            for(int i=0;i<engine.n;i++){
                if(i)out<<",";
                out<<2*witness[i]+int(witnessNegative>>i&1);
            }out<<"],\"full_expanded_integer_correlations_verified\":true";
            std::ofstream txt(output+".txt");
            txt<<engine.p<<" "<<engine.m<<" "<<engine.mu<<" "<<1+engine.m*cfg.headNet<<" "<<(engine.m/3)*cfg.tailNet<<" 0\n";
            for(int i=0;i<engine.n;i++)txt<<2*witness[i]+int(witnessNegative>>i&1)<<" ";
            txt<<"\n";
        } out<<",\"T_sequence_verified\":"<<(!witness.empty()?"true":"false")<<",\"new_order_claimed\":false}\n";
        out.close();
        if(std::rename(temporary.c_str(),output.c_str()))throw std::runtime_error("atomic report failed");
    };
    std::cout<<"{\"stage\":\"field_sign_search_started\",\"pid\":"<<getpid()<<",\"seed\":"<<seed<<"}"<<std::endl;
    while(true){
        if(sampled%256==0){
            double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if((sampleOnly&&sampled>=sampleOnly)||(!sampleOnly&&elapsed>=seconds))break;
            if(elapsed-lastout>5){
                report();
                std::cout<<"{\"seconds\":"<<elapsed<<",\"field_parameters\":"<<sampled<<",\"lifts\":"<<tested<<",\"mod8_passed\":"<<mod8<<",\"DC_passed\":"<<dc<<"}"<<std::endl;
                lastout=elapsed;
            }
        } uint32_t p=rng()&mask,q=rng()&mask,w=rng()&mask;
        if(!(p|q))continue;
        uint32_t b0,b1;
        if(p){
            b0=w;
            b1=mul(1^mul(q,w),invert(p));
        }else{
            b0=invert(q);
            b1=w;
        }if((mul(p,b1)^mul(q,b0))!=1)throw std::runtime_error("SL2 determinant failed");
        ColorWord spectral=ColorWord(p)|(ColorWord(q)<<degree)|(ColorWord(b0)<<(2*degree))|(ColorWord(b1)<<(3*degree)),colors=0;
        for(int i=0;i<colorBytes;i++)colors^=lookup[i][(spectral>>(8*i))&255];
        std::vector<int>axes(engine.n);
        int k=0;
        for(int i=1;i<engine.n;i++){
            axes[i]=(colors>>(2*(i-1)))&3;
            k+=axes[i]==0;
        }sampled++;
        heads[k]++;
        if(examples.size()<64){
            examples.push_back("{\"p_field\":"+std::to_string(p)+",\"q_field\":"+std::to_string(q)+",\"w_field\":"+std::to_string(w)+",\"b0\":"+std::to_string(b0)+",\"b1\":"+std::to_string(b1)+",\"spectral\":"+color_decimal(spectral)+",\"color_bits\":"+color_decimal(colors)+"}");
        } if(sampleOnly)continue;
        if(axes[1])continue;
        gauged++;
        if(!allowHead[k]){
            headRejected++;
            continue;
        }cfg.validateMod8=tested<32;
        auto result=engine.lift(axes,cfg);
        tested++;
        affine+=result.affineAssignments;
        mod8+=result.mod8Passed;
        dc+=result.DCpassed;
        if(result.status=="odd unsigned correlation")throw std::runtime_error("field axis failed correlation parity");
        if(result.status=="inconsistent mod4")inconsistent++;
        if(result.status=="nullity limit")rankLimit++;
        if(!result.solutions.empty()){
            witness=axes;
            witnessNegative=result.solutions[0];
            report();
            std::cout<<"{\"candidate\":true,\"full_expanded_integer_correlations_verified\":true,\"path\":\""<<output<<"\"}"<<std::endl;
            return 0;
        }
    } report();
    std::cout<<"{\"candidate\":false,\"field_parameters\":"<<sampled<<",\"axis_patterns_sign_lifted\":"<<tested<<",\"mod8_passed\":"<<mod8<<",\"DC_passed_after_mod8\":"<<dc<<",\"output\":\""<<output<<"\"}"<<std::endl;
    return 0;
}  // Reusable native preparation plus mod8-filtered finite sign enumeration.
// native_sign_lift T_TABLE_BIN AXIS_POOL_BIN OUTPUT_JSON [STOP_FIRST=1] [VALIDATE_MOD8=0] [USE_MOD8=1]
#include <chrono>
#include <iostream>
#include <unistd.h>
int pool_main(int argc,char**argv){
    if(argc<4)return 2;
    twisted_sign_lift::Engine engine(argv[1]);
    std::ifstream pool(argv[2],std::ios::binary);
    int32_t h[9];
    pool.read((char*)h,sizeof(h));
    if(!pool||h[0]!=0x41584c34||h[1]!=engine.p||h[2]!=engine.m||h[3]!=engine.n)return 3;
    int count=h[4];
    twisted_sign_lift::Config cfg;
    cfg.headNet=h[5];
    cfg.tailNet=h[6];
    cfg.firstNegative=h[7];
    cfg.maxNullity=h[8];
    cfg.stopFirst=argc<5||atoi(argv[4])!=0;
    cfg.validateMod8=argc>5&&atoi(argv[5]);
    cfg.mod8=argc<7||atoi(argv[6])!=0;
    cfg.collectMod8=argc>7&&atoi(argv[7]);
    auto started=std::chrono::steady_clock::now();
    uint64_t affine=0,mod8=0,dc=0,solutions=0;
    std::vector<std::string>records;
    std::vector<int>firstStates;
    int completed=0;
    for(int id=0;id<count;id++){
        std::vector<uint8_t>bytes(engine.n);
        pool.read((char*)bytes.data(),bytes.size());
        if(!pool)return 4;
        std::vector<int>axes(bytes.begin(),bytes.end());
        auto r=engine.lift(axes,cfg);
        completed++;
        affine+=r.affineAssignments;
        mod8+=r.mod8Passed;
        dc+=r.DCpassed;
        solutions+=r.solutions.size();
        std::string record="{\"pattern\":"+std::to_string(id)+",\"status\":\""+r.status+"\",\"rank\":"+std::to_string(r.rank)+",\"nullity\":"+std::to_string(r.nullity)+",\"negative_masks\":[";
        for(size_t j=0;j<r.solutions.size();j++){
            if(j)record+=",";
            record+=std::to_string(r.solutions[j]);
        }record+="],\"mod8_negative_masks\":[";
        for(size_t j=0;j<r.mod8Masks.size();j++){
            if(j)record+=",";
            record+=std::to_string(r.mod8Masks[j]);
        }record+="]}";
        records.push_back(record);
        if(!r.solutions.empty()&&firstStates.empty()){
            for(int i=0;i<engine.n;i++)firstStates.push_back(2*axes[i]+int(r.solutions[0]>>i&1));
        } if(cfg.stopFirst&&!r.solutions.empty())break;
    } std::ofstream out(argv[3]);
    out<<"{\"p\":"<<engine.p<<",\"m\":"<<engine.m<<",\"mu\":"<<engine.mu<<",\"pid\":"<<getpid()<<",\"patterns_completed\":"<<completed<<",\"patterns_available\":"<<count<<",\"affine_sign_assignments_covered\":"<<affine<<",\"mod8_passed\":"<<mod8<<",\"DC_passed_after_filter\":"<<dc<<",\"solutions\":"<<solutions<<",\"mod8_enabled\":"<<(cfg.mod8?"true":"false")<<",\"mod8_random_validation_enabled\":"<<(cfg.validateMod8?"true":"false")<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\"patterns\":[";
    for(size_t i=0;i<records.size();i++){
        if(i)out<<",";
        out<<records[i];
    }out<<"]";
    if(!firstStates.empty()){
        out<<",\"states\":[";
        for(size_t i=0;i<firstStates.size();i++){
            if(i)out<<",";
            out<<firstStates[i];
        }out<<"],\"full_expanded_integer_correlations_verified\":true";
        std::ofstream text(std::string(argv[3])+".txt");
        text<<engine.p<<" "<<engine.m<<" "<<engine.mu<<" "<<1+engine.m*cfg.headNet<<" "<<(engine.m/3)*cfg.tailNet<<" 0\n";
        for(int s:firstStates)text<<s<<" ";
        text<<"\n";
    } out<<",\"new_order_claimed\":false}\n";
    out.close();
    std::cout<<"{\"patterns\":"<<completed<<",\"affine_covered\":"<<affine<<",\"mod8_passed\":"<<mod8<<",\"DC_passed_after_filter\":"<<dc<<",\"solutions\":"<<solutions<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\"output\":\""<<argv[3]<<"\"}"<<std::endl;
    return 0;
}
#include <functional>
#include <set>
struct ConditionalSampler {
    int degree,n; uint32_t polynomial,size,mask;
    std::vector<ColorWord>columns;
    std::vector<std::array<ColorWord,256>>lookup;
    explicit ConditionalSampler(const std::string&path){
        std::ifstream input(path,std::ios::binary);int32_t h[7];input.read((char*)h,sizeof(h));
        if(!input||h[0]!=0x46503135||h[1]!=3||h[3]<1||h[3]>25||h[5]!=0||h[6]<=0||h[4])
            throw std::runtime_error("headless sampler requires positive-tail DC and first sign positive");
        degree=h[3];polynomial=h[2];size=1u<<degree;mask=size-1;n=2*degree+1;
        columns.resize(4*degree);
        for(auto&value:columns){uint64_t words[2];input.read((char*)words,16);value=ColorWord(words[0])|(ColorWord(words[1])<<64);}
        if(!input)throw std::runtime_error("truncated headless map");
        lookup.resize((4*degree+7)/8);
        for(int i=0;i<int(lookup.size());i++)for(int v=1;v<256;v++){
            int b=__builtin_ctz(v);lookup[i][v]=lookup[i][v^(1<<b)]^((8*i+b<4*degree)?columns[8*i+b]:0);
        }
    }
    uint32_t mul(uint32_t a,uint32_t b)const{
        uint32_t value=0;while(b){if(b&1)value^=a;b>>=1;a<<=1;if(a&size)a^=polynomial;}return value;
    }
    uint32_t inverse(uint32_t x)const{
        if(!x)throw std::runtime_error("zero inverse");uint32_t a=x,b=polynomial;uint64_t u=1,v=0;
        while(a!=1){if(!a)throw std::runtime_error("reducible field polynomial");int j=(31-__builtin_clz(a))-(31-__builtin_clz(b));if(j<0){std::swap(a,b);std::swap(u,v);j=-j;}a^=b<<j;u^=v<<j;}
        while(u>=size)u^=uint64_t(polynomial)<<((63-__builtin_clzll(u))-degree);
        if(mul(x,uint32_t(u))!=1)throw std::runtime_error("bad field inverse");return uint32_t(u);
    }
    ColorWord map(ColorWord value)const{
        ColorWord out=0;for(int i=0;i<int(lookup.size());i++)out^=lookup[i][(value>>(8*i))&255];return out;
    }
    // Calls visit for EVERY allowed w. False from visit stops at an exact zero.
    bool each(uint32_t p,uint32_t q,const std::function<bool(uint32_t,ColorWord)>&visit,uint64_t&space)const{
        if(!(p|q))throw std::runtime_error("zero SL2 row");
        std::vector<ColorWord>delta(degree);ColorWord constant;
        if(p){
            uint32_t ip=inverse(p),slope=mul(q,ip);
            constant=map(ColorWord(p)|(ColorWord(q)<<degree)|(ColorWord(ip)<<(3*degree)));
            for(int j=0;j<degree;j++)delta[j]=map((ColorWord(1)<<(2*degree+j))|(ColorWord(mul(slope,1u<<j))<<(3*degree)));
        }else{
            constant=map((ColorWord(q)<<degree)|(ColorWord(inverse(q))<<(2*degree)));
            for(int j=0;j<degree;j++)delta[j]=map(ColorWord(1)<<(3*degree+j));
        }
        // First orbit is positive axis1: solve its two affine color equations.
        std::vector<uint64_t>rows(2);std::vector<int>pivot;
        for(int bit=0;bit<2;bit++){
            for(int j=0;j<degree;j++)rows[bit]|=uint64_t((delta[j]>>bit)&1)<<j;
            rows[bit]|=uint64_t(int((constant>>bit)&1)^(bit==0))<<degree;
        }
        for(int col=0;col<degree;col++){
            int at=pivot.size(),found=-1;for(int i=at;i<2;i++)if(rows[i]>>col&1){found=i;break;}
            if(found<0)continue;std::swap(rows[at],rows[found]);for(int i=0;i<2;i++)if(i!=at&&(rows[i]>>col&1))rows[i]^=rows[at];pivot.push_back(col);
        }
        for(uint64_t row:rows)if(!(row&mask)&&(row>>degree&1))return true;
        uint32_t particular=0;for(int i=0;i<int(pivot.size());i++)if(rows[i]>>degree&1)particular|=1u<<pivot[i];
        std::vector<uint32_t>basis;
        for(int col=0;col<degree;col++)if(std::find(pivot.begin(),pivot.end(),col)==pivot.end()){
            uint32_t v=1u<<col;for(int i=0;i<int(pivot.size());i++)if(rows[i]>>col&1)v|=1u<<pivot[i];basis.push_back(v);
        }
        ColorWord base=constant;for(uint32_t z=particular;z;z&=z-1)base^=delta[__builtin_ctz(z)];
        std::vector<ColorWord>changes(basis.size());for(int i=0;i<int(basis.size());i++)for(uint32_t z=basis[i];z;z&=z-1)changes[i]^=delta[__builtin_ctz(z)];
        int d=basis.size(),low=std::min(6,d),high=d-low;space+=1ULL<<d;
        uint64_t valid=low==6?~0ULL:((1ULL<<(1<<low))-1);
        std::vector<uint64_t>Z(low),truth(1<<low),left(n-1),right(n-1);
        for(int i=0;i<low;i++)for(int z=0;z<(1<<low);z++)if(z>>i&1)Z[i]|=1ULL<<z;
        for(int z=1;z<(1<<low);z++){int j=__builtin_ctz(z);truth[z]=truth[z^(1<<j)]^Z[j];}
        for(int i=0;i<n-1;i++){
            int a=0,b=0;for(int j=0;j<low;j++){a|=int((changes[j]>>(2*i))&1)<<j;b|=int((changes[j]>>(2*i+1))&1)<<j;}left[i]=truth[a];right[i]=truth[b];
        }
        ColorWord highColors=base;
        for(uint64_t counter=0;counter<(1ULL<<high);counter++){
            if(counter)highColors^=changes[low+__builtin_ctzll(counter)];
            uint64_t good=valid;ColorWord colors=highColors;
            for(int i=0;i<n-1&&good;i++,colors>>=2){unsigned bits=colors&3;good&=(left[i]^(-uint64_t(bits&1)))|(right[i]^(-uint64_t((bits>>1)&1)));}
            while(good){int a=__builtin_ctzll(good);good&=good-1;uint64_t z=uint64_t(a)|((counter^(counter>>1))<<low);uint32_t w=particular;ColorWord result=highColors;
                for(uint64_t b=z;b;b&=b-1)w^=basis[__builtin_ctzll(b)];for(unsigned b=a;b;b&=b-1)result^=changes[__builtin_ctz(b)];
                if((result&3)!=1)throw std::runtime_error("first-tail gauge failed");
                for(int i=0;i<n-1;i++)if(((result>>(2*i))&3)==0)throw std::runtime_error("zero color escaped filter");
                if(!visit(w,result))return false;
            }
        }
        return true;
    }
};
int headless_main(int argc,char**argv){try{
    if(argc<7)throw std::runtime_error("headless TABLE COLUMNS OUTPUT SECONDS SEED MODE(search|exhaustive)");
    twisted_sign_lift::Engine engine(argv[1]);ConditionalSampler sampler(argv[2]);if(engine.n!=sampler.n)throw std::runtime_error("chart size mismatch");
    std::ifstream input(argv[2],std::ios::binary);int32_t cfgData[7];input.read((char*)cfgData,sizeof(cfgData));twisted_sign_lift::Config cfg;cfg.headNet=0;cfg.firstAxis=1;cfg.tailNet=cfgData[6];cfg.firstNegative=0;cfg.collectMod8=false;cfg.stopFirst=std::string(argv[6])!="exhaustive";
    bool exhaustive=!cfg.stopFirst;if(exhaustive&&sampler.degree>5)throw std::runtime_error("exhaustive control too large");
    std::mt19937_64 rng(strtoull(argv[5],nullptr,10));double limit=atof(argv[4]);auto start=std::chrono::steady_clock::now();uint64_t pairs=0,space=0,axis=0,affine=0,mod8=0,dc=0;std::vector<int>found;std::vector<std::string>records;
    if(std::string(argv[6])=="checkpair"){
        uint32_t p=argc>7?uint32_t(strtoul(argv[7],nullptr,10)):uint32_t(rng()&sampler.mask);
        uint32_t q=argc>8?uint32_t(strtoul(argv[8],nullptr,10)):uint32_t(rng()&sampler.mask);
        if(!(p|q))q=1;
        uint32_t base=p?sampler.inverse(p):sampler.inverse(q),slope=p?sampler.mul(q,base):0;
        auto at=[&](uint32_t w){uint32_t u=p?w:base,v=p?(base^sampler.mul(slope,w)):w;
            return sampler.map(ColorWord(p)|(ColorWord(q)<<sampler.degree)|(ColorWord(u)<<(2*sampler.degree))|(ColorWord(v)<<(3*sampler.degree)));};
        std::set<uint32_t>expected,actual;
        auto timing=std::chrono::steady_clock::now();
        for(uint32_t w=0;w<sampler.size;w++){
            ColorWord colors=at(w);if((colors&3)!=1)continue;bool good=true;
            for(int i=0;i<sampler.n-1;i++)if(((colors>>(2*i))&3)==0){good=false;break;}
            if(good)expected.insert(w);
        }
        double slow=std::chrono::duration<double>(std::chrono::steady_clock::now()-timing).count();timing=std::chrono::steady_clock::now();
        sampler.each(p,q,[&](uint32_t w,ColorWord colors){if(colors!=at(w))throw std::runtime_error("conditional color mismatch");actual.insert(w);return true;},space);
        double fast=std::chrono::duration<double>(std::chrono::steady_clock::now()-timing).count();
        if(expected!=actual)throw std::runtime_error("conditional w set differs from exhaustive enumeration");
        std::ofstream out(argv[3]);out<<"{\"p_field\":"<<p<<",\"q_field\":"<<q<<",\"all_w_checked\":"<<sampler.size<<",\"first_axis_affine_space_size\":"<<space<<",\"headless_w_count\":"<<actual.size()<<",\"complete_w_sets_equal\":true,\"naive_seconds\":"<<slow<<",\"conditional_seconds\":"<<fast<<"}";
        std::cout<<"{\"headless_complete_pair_check\":true,\"points\":"<<sampler.size<<",\"survivors\":"<<actual.size()<<"}"<<std::endl;return 0;
    }
    auto visit=[&](uint32_t w,ColorWord colors){
        std::vector<int>a(engine.n);for(int i=1;i<engine.n;i++)a[i]=(colors>>(2*(i-1)))&3;axis++;cfg.validateMod8=axis<32;auto lift=engine.lift(a,cfg);affine+=lift.affineAssignments;mod8+=lift.mod8Passed;dc+=lift.DCpassed;
        if(exhaustive){std::string s="{\"colors\":"+color_decimal(colors)+",\"w\":"+std::to_string(w)+",\"negative_masks\":[";for(size_t i=0;i<lift.solutions.size();i++){if(i)s+=",";s+=std::to_string(lift.solutions[i]);}records.push_back(s+"]}");}
        if(!lift.solutions.empty()&&found.empty()){for(int i=0;i<engine.n;i++)found.push_back(2*a[i]+int(lift.solutions[0]>>i&1));}
        return found.empty()||exhaustive;
    };
    std::cout<<"{\"stage\":\"headless_started\",\"pid\":"<<getpid()<<"}"<<std::endl;
    if(exhaustive){for(uint32_t p=0;p<sampler.size;p++)for(uint32_t q=0;q<sampler.size;q++)if(p|q){pairs++;sampler.each(p,q,visit,space);}}
    else while(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<limit&&found.empty()){
        uint32_t p=rng()&sampler.mask,q=rng()&sampler.mask;if(!(p|q))continue;pairs++;sampler.each(p,q,visit,space);
    }
    std::ofstream out(argv[3]);out<<"{\"p\":"<<engine.p<<",\"m\":"<<engine.m<<",\"mu\":"<<engine.mu<<",\"pairs\":"<<pairs<<",\"conditional_w_cardinality_sum\":"<<space<<",\"headless_axes\":"<<axis<<",\"sign_kernel_cardinality_sum\":"<<affine<<",\"mod8\":"<<mod8<<",\"DC\":"<<dc<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"headless_first_axis\":1,\"first_sign_negative\":false,\"records\":[";for(size_t i=0;i<records.size();i++){if(i)out<<",";out<<records[i];}out<<"]";
    if(!found.empty()){out<<",\"states\":[";for(size_t i=0;i<found.size();i++){if(i)out<<",";out<<found[i];}out<<"],\"T_sequence_verified\":true";}
    out<<",\"new_order_claimed\":false}";out.close();std::cout<<"{\"pairs\":"<<pairs<<",\"axes\":"<<axis<<",\"output\":\""<<argv[3]<<"\"}"<<std::endl;
    return 0;
}catch(const std::exception&error){std::cerr<<error.what()<<std::endl;return 1;}}

int main(int argc,char**argv){
    try{
        if(argc<2)throw std::runtime_error("usage: hadamard_frame_search search|lift ...");
        std::string mode=argv[1];
        if(mode=="search")return search_main(argc-1,argv+1);
        if(mode=="lift")return pool_main(argc-1,argv+1);
        if(mode=="headless")return headless_main(argc-1,argv+1);
        throw std::runtime_error("unknown mode: "+mode);
    }catch(const std::exception&error){
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
}
