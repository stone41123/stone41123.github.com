#include<cstdio>
#include<cstring>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#define LL long long
using namespace std;
inline int read(){
	int x=0,f=1;char ch=' ';
	while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while(ch>='0'&&ch<='9')x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
	return x*f;
}
const int N=2e4+5,mod=1004535809,g=3,gi=334845270;
int n,m,x,S,G,vis[N],H[N],R[N],bin[1<<19];
LL C[N],ans[N];
inline LL ksm(LL a,LL n,int p=mod){
    LL ans=1LL;
    while(n){
        if(n&1)ans=ans*a%p;
        a=a*a%p;
        n>>=1;
    }
    return ans;
}
inline void init(){
    for(int i=1;i<m;++i){
        memset(vis,0,sizeof(vis));
        int flag=1;
        for(int j=1;j<m;++j){
            int tmp=ksm(i,j,m);
            if(vis[tmp]){flag=0;break;}
            vis[tmp]=1;
        }
        if(flag){G=i;break;}
    }
    for(int i=0;i<m;++i)H[ksm(G,i,m)]=i;
    for(int i=0;i<18;++i)bin[1<<i]=i;
}
inline void NTT(LL *a,int n,int f){
    int L=bin[n];
    for(int i=0;i<n;++i)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
    for(int i=0;i<n;++i)if(i<R[i])swap(a[i],a[R[i]]);
    for(int i=1;i<n;i<<=1){
        LL wn=ksm(f==1?g:gi,(mod-1)/(i<<1));
        for(int j=0;j<n;j+=(i<<1)){
            LL w=1LL;
            for(int k=0;k<i;++k,w=w*wn%mod){
                LL x=a[j+k],y=w*a[j+k+i]%mod;
                a[j+k]=(x+y)%mod;a[j+k+i]=(x-y+mod)%mod;
            }
        }
    }
    if(f==-1)for(int i=0;i<n;++i)a[i]=a[i]*ksm(n,mod-2)%mod;
}
inline void poly_ksm(LL *a,LL *ans,int n,int b){
    for(int i=0;i<n;++i)ans[i]=a[i];b--;
    while(b){
        if(b&1){
            NTT(ans,n,1);NTT(a,n,1);
            for(int i=0;i<n;++i)ans[i]=ans[i]*a[i]%mod;
            NTT(ans,n,-1);NTT(a,n,-1);
            for(int i=m-1;i<n;++i)ans[i%(m-1)]=(ans[i%(m-1)]+ans[i])%mod;
            for(int i=m-1;i<n;++i)ans[i]=0;
        }
        NTT(a,n,1);
        for(int i=0;i<n;++i)a[i]=a[i]*a[i]%mod;
        NTT(a,n,-1);
        for(int i=m-1;i<n;++i)a[i%(m-1)]=(a[i%(m-1)]+a[i])%mod;
        for(int i=m-1;i<n;++i)a[i]=0;
        b>>=1;
    }
}
int main(){
    n=read();m=read();x=read();S=read();
    init();
    for(int i=1;i<=S;++i){
        int tmp=read();
        if(tmp)C[H[tmp]%(m-1)]=1;
    }
    int Len;
    for(Len=1;Len<(m<<1);Len<<=1);
    poly_ksm(C,ans,Len,n);
    printf("%lld",ans[H[x]%(m-1)]);
    return 0;
}
