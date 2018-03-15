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
const int N=3e5+5,mod=998244353,g=3,gi=332748118;
int n,m,L,R[N],bin[1<<21];
LL a[N],b[N],tmp[N];
inline LL ksm(LL a,LL n){
    LL ans=1LL;
    while(n){
        if(n&1)ans=ans*a%mod;
        a=a*a%mod;
        n>>=1;
    }
    return ans;
}
inline void NTT(LL *a,int n,int f){
    int L=bin[n];
    for(int i=0;i<n;++i)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
    for(int i=0;i<n;++i)if(i<R[i])swap(a[i],a[R[i]]);
    for(int i=1;i<n;i<<=1){
        LL wn=ksm(f==1?g:gi,(mod-1)/(i<<1));
        for(int j=0;j<n;j+=(i<<1)){
            LL w=1;
            for(int k=0;k<i;++k,w=w*wn%mod){
                LL x=a[j+k],y=w*a[j+k+i]%mod;
                a[j+k]=(x+y)%mod;a[j+k+i]=(x-y+mod)%mod;
            }
        }
    }
    if(f==-1)for(int i=0;i<n;++i)a[i]=a[i]*ksm(n,mod-2)%mod;
}
inline void inv(LL *a,LL *b,int n){
    if(n==1){
        b[0]=ksm(a[0],mod-2);
        return;
    }
    inv(a,b,(n+1)>>1);
    int m=1;
    for(;m<(n<<1);m<<=1);
    for(int i=0;i<n;++i)tmp[i]=a[i];
    for(int i=n;i<m;++i)tmp[i]=0;
    NTT(tmp,m,1);NTT(b,m,1);
    for(int i=0;i<m;++i)
        b[i]=(2*b[i]%mod-tmp[i]*b[i]%mod*b[i]%mod+mod)%mod;
    NTT(b,m,-1);
    for(int i=n;i<m;++i)b[i]=0;
    return;
}
int main(){
    for(int i=0;i<=19;++i)bin[1<<i]=i;
    n=read();
    for(m=1;m<=n;m<<=1)L++;
    for(int i=0;i<n;++i)a[i]=read();
    inv(a,b,m);
    for(int i=0;i<n;++i)printf("%d ",b[i]);
    return 0;
}
