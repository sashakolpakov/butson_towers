#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Project numeric orbit scores onto the four exact sign-count categories.
// scores[i,axis] is the full orbit inner product with +e_axis.
// Categories are +head,-head,+tail,-tail. Each tail category first chooses
// its best coordinate; a square assignment then enforces all four counts.
extern "C" int t_sequence_assign(int n,const double* scores,const int* counts,int* states){
    if(n<1||n>2000||!scores||!counts||!states)return 1;
    int total=0;for(int c=0;c<4;++c){if(counts[c]<0)return 2;total+=counts[c];}
    if(total!=n)return 3;
    std::vector<int> category;category.reserve(n);
    for(int c=0;c<4;++c)for(int j=0;j<counts[c];++j)category.push_back(c);
    std::vector<double> values(size_t(n)*4);
    std::vector<int> choice(size_t(n)*4);
    for(int i=0;i<n;++i){
        for(int axis=0;axis<4;++axis)if(!std::isfinite(scores[4*i+axis]))return 4;
        values[4*i]=scores[4*i];choice[4*i]=0;
        values[4*i+1]=-scores[4*i];choice[4*i+1]=1;
        int plus=1,minus=1;
        for(int axis=2;axis<4;++axis){
            if(scores[4*i+axis]>scores[4*i+plus])plus=axis;
            if(scores[4*i+axis]<scores[4*i+minus])minus=axis;
        }
        values[4*i+2]=scores[4*i+plus];choice[4*i+2]=2*plus;
        values[4*i+3]=-scores[4*i+minus];choice[4*i+3]=2*minus+1;
    }
    std::vector<double> u(n+1),v(n+1),minimum(n+1);
    std::vector<int> matching(n+1),way(n+1);
    std::vector<unsigned char> used(n+1);
    for(int row=1;row<=n;++row){
        matching[0]=row;int column=0;
        std::fill(minimum.begin(),minimum.end(),std::numeric_limits<double>::infinity());
        std::fill(used.begin(),used.end(),0);
        do{
            used[column]=1;int i=matching[column],next=0;
            double delta=std::numeric_limits<double>::infinity();
            for(int j=1;j<=n;++j)if(!used[j]){
                double cost=-values[4*(i-1)+category[j-1]]-u[i]-v[j];
                if(cost<minimum[j]){minimum[j]=cost;way[j]=column;}
                if(minimum[j]<delta){delta=minimum[j];next=j;}
            }
            if(!std::isfinite(delta))return 5;
            for(int j=0;j<=n;++j)if(used[j]){u[matching[j]]+=delta;v[j]-=delta;}
                else minimum[j]-=delta;
            column=next;
        }while(matching[column]);
        do{int previous=way[column];matching[column]=matching[previous];column=previous;}while(column);
    }
    for(int j=1;j<=n;++j){int i=matching[j]-1;states[i]=choice[4*i+category[j-1]];}
    return 0;
}

// Exact numeric projection with no prescribed head-support size. The two
// integer resources are total positive signs and the signed head count.
// Total positive minus negative signs is head_net+tail_net, so fixing these
// two resources enforces both DC sums while allowing every feasible head size.
extern "C" int t_sequence_assign_free_heads(int n,const double* scores,int head_net,int tail_net,int* states){
    if(n<1||n>2000||!scores||!states)return 1;
    const int positive_total=n+head_net+tail_net;
    if(positive_total%2)return 2;
    const int positive=positive_total/2,negative=n-positive;
    if(positive<0||negative<0||std::abs(head_net)>n)return 3;
    const int width=2*n+1,cells=(positive+1)*width;
    const double absent=-std::numeric_limits<double>::infinity();
    std::vector<double> values(size_t(n)*4),previous(cells,absent),current(cells,absent);
    std::vector<int> choice(size_t(n)*4);
    std::vector<unsigned char> decision(size_t(n+1)*cells,255);
    for(int i=0;i<n;++i){
        for(int axis=0;axis<4;++axis)if(!std::isfinite(scores[4*i+axis]))return 4;
        values[4*i]=scores[4*i];choice[4*i]=0;
        values[4*i+1]=-scores[4*i];choice[4*i+1]=1;
        int plus=1,minus=1;
        for(int axis=2;axis<4;++axis){
            if(scores[4*i+axis]>scores[4*i+plus])plus=axis;
            if(scores[4*i+axis]<scores[4*i+minus])minus=axis;
        }
        values[4*i+2]=scores[4*i+plus];choice[4*i+2]=2*plus;
        values[4*i+3]=-scores[4*i+minus];choice[4*i+3]=2*minus+1;
    }
    auto index=[&](int p,int h){return p*width+h+n;};
    previous[index(0,0)]=0;
    const int dp[4]={1,0,1,0},dh[4]={1,-1,0,0};
    for(int i=0;i<n;++i){
        std::fill(current.begin(),current.end(),absent);
        for(int p=std::max(0,i-negative);p<=std::min(i,positive);++p){
            for(int h=-(i-p);h<=p;++h){
                double value=previous[index(p,h)];if(!std::isfinite(value))continue;
                for(int category=0;category<4;++category){
                    const int np=p+dp[category],nh=h+dh[category],remaining=n-i-1;
                    if(np>positive||i+1-np>negative||np+remaining<positive||std::abs(head_net-nh)>remaining)continue;
                    const int next=index(np,nh);double candidate=value+values[4*i+category];
                    if(candidate>current[next]){current[next]=candidate;decision[size_t(i+1)*cells+next]=category;}
                }
            }
        }
        previous.swap(current);
    }
    if(!std::isfinite(previous[index(positive,head_net)]))return 5;
    int p=positive,h=head_net;
    for(int i=n-1;i>=0;--i){
        int category=decision[size_t(i+1)*cells+index(p,h)];if(category<0||category>3)return 6;
        states[i]=choice[4*i+category];p-=dp[category];h-=dh[category];
    }
    return p==0&&h==0?0:7;
}
