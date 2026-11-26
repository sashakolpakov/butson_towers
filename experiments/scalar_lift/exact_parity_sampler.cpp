// Exact determinant-one sampling over one binary field, with packed colors.
#include <algorithm>
#include <cstdint>

using wide = unsigned __int128;

static unsigned degree64(uint64_t x) { return 63u-__builtin_clzll(x); }
static unsigned degree128(wide x) {
    uint64_t upper=uint64_t(x>>64);
    return upper?64u+degree64(upper):degree64(uint64_t(x));
}
static uint64_t multiply(uint64_t a,uint64_t b,unsigned d,uint64_t modulus) {
    uint64_t result=0;
    while(b) {
        if(b&1)result^=a;
        b>>=1;a<<=1;
        if(a&(uint64_t(1)<<d))a^=modulus;
    }
    return result;
}
static uint64_t inverse(uint64_t a,unsigned d,uint64_t modulus) {
    uint64_t u=a,v=modulus;
    wide g1=1,g2=0;
    while(u!=1) {
        if(!u)return 0;
        int j=int(degree64(u))-int(degree64(v));
        if(j<0){std::swap(u,v);std::swap(g1,g2);j=-j;}
        u^=v<<j;g1^=g2<<j;
    }
    while(g1&&(degree128(g1)>=d))g1^=wide(modulus)<<(degree128(g1)-d);
    return uint64_t(g1);
}
static bool parameters(unsigned d,uint64_t modulus) {
    return d>=1&&d<=63&&(modulus&1)&&(modulus>>d)==1;
}
extern "C" int exact_parity_field(unsigned d,uint64_t modulus,uint64_t a,uint64_t b,
                                   uint64_t* product,uint64_t* reciprocal) {
    if(!parameters(d,modulus)||!product||!reciprocal||(a>>d)||(b>>d))return 1;
    *product=multiply(a,b,d,modulus);
    *reciprocal=a?inverse(a,d,modulus):0;
    return a&&!*reciprocal?2:0;
}
static int colors(unsigned d,uint64_t modulus,const uint64_t* columns,
                  uint64_t p,uint64_t q,uint64_t w,uint64_t* output) {
    if(!(p|q)||(p>>d)||(q>>d)||(w>>d))return 1;
    uint64_t b0=w,b1=0;
    if(p) {
        uint64_t inv=inverse(p,d,modulus);
        if(!inv||multiply(p,inv,d,modulus)!=1)return 2;
        b1=multiply(1^multiply(q,w,d,modulus),inv,d,modulus);
    } else {
        b0=inverse(q,d,modulus);b1=w;
        if(!b0||multiply(q,b0,d,modulus)!=1)return 2;
    }
    const uint64_t values[4]={p,q,b0,b1};
    std::fill(output,output+4,0);
    for(unsigned block=0;block<4;++block)for(uint64_t bits=values[block];bits;bits&=bits-1) {
        unsigned bit=__builtin_ctzll(bits),offset=4*(block*d+bit);
        for(unsigned word=0;word<4;++word)output[word]^=columns[offset+word];
    }
    return 0;
}
extern "C" int exact_parity_colors(unsigned d,uint64_t modulus,const uint64_t* columns,
                                    uint64_t p,uint64_t q,uint64_t w,uint64_t* output) {
    if(!parameters(d,modulus)||!columns||!output)return 1;
    return colors(d,modulus,columns,p,q,w,output);
}
// xoshiro256** 1.0, adapted from the public-domain reference by
// David Blackman and Sebastiano Vigna: https://prng.di.unimi.it/xoshiro256starstar.c
static uint64_t rotate(uint64_t x,unsigned k){return (x<<k)|(x>>(64-k));}
static uint64_t random64(uint64_t* state) {
    uint64_t result=rotate(state[1]*5,7)*9,t=state[1]<<17;
    state[2]^=state[0];state[3]^=state[1];state[1]^=state[2];state[0]^=state[3];
    state[2]^=t;state[3]=rotate(state[3],45);
    return result;
}
extern "C" int exact_parity_sample(unsigned d,uint64_t modulus,const uint64_t* columns,
                                    const uint8_t* allowed_heads,int first_head,
                                    uint64_t max_draws,uint64_t* state,uint64_t* draws,
                                    uint64_t* output) {
    if(!parameters(d,modulus)||!columns||!allowed_heads||!state||!draws||!output)return 2;
    if(!(state[0]|state[1]|state[2]|state[3]))return 7;
    const uint64_t mask=(uint64_t(1)<<d)-1;
    for(uint64_t trial=0;trial<max_draws;++trial) {
        uint64_t p,q;
        do{p=random64(state)&mask;q=random64(state)&mask;}while(!(p|q));
        uint64_t w=random64(state)&mask;
        int error=colors(d,modulus,columns,p,q,w,output);
        if(error)return 2+error;
        ++*draws;
        if(first_head&&(output[0]&3))continue;
        unsigned occupied=0;
        for(unsigned word=0;word<4;++word)
            occupied+=__builtin_popcountll((output[word]|(output[word]>>1))&UINT64_C(0x5555555555555555));
        if(occupied>2*d)return 6;
        if(allowed_heads[2*d-occupied])return 0;
    }
    return 1;
}
