#include<cstdio>
#include<cstring>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#define LL long long
#define max(a,b) a>b?a:b
#define min(a,b) a<b?a:b
using namespace std;
inline int read(){
	int x=0,f=1;char ch=' ';
	while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while(ch>='0'&&ch<='9')x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
	return x*f;
}
const int N=3e6,mod=998244353,g=3,gi=332748118;
int n,m,L,R[N];
LL inv,A[N],B[N];
inline LL ksm(LL a,LL n){
    LL ans=1;
    while(n){
        if(n&1)ans=ans*a%mod;
        a=a*a%mod;
        n>>=1;
    }
    return ans;
}
inline void NTT(LL *a,int f){
    for(int i=0;i<n;++i)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
    for(int i=0;i<n;++i)if(i<R[i])swap(a[i],a[R[i]]);
    for(int i=1;i<n;i<<=1){
        LL wn=ksm((f==1)?g:gi,(mod-1)/(i<<1));
        for(int j=0;j<n;j+=(i<<1)){
            LL w=1;
            for(int k=0;k<i;++k,w=w*wn%mod){
                LL x=a[j+k],y=w*a[j+k+i]%mod;
                a[j+k]=(x+y)%mod;a[j+k+i]=(x-y+mod)%mod;
            }
        }
    }
    if(f==-1)for(int i=0;i<n;++i)a[i]=a[i]*inv%mod;
}
int main(){
    scanf("%d %d",&n,&m);
    for(int i=0;i<=n;++i)A[i]=read(),A[i]=(A[i]+mod)%mod;
    for(int i=0;i<=m;++i)B[i]=read(),B[i]=(B[i]+mod)%mod;
    m+=n;
    for(n=1;n<=m;n<<=1)++L;
    inv=ksm(n,mod-2);
    NTT(A,1);NTT(B,1);
    for(int i=0;i<n;++i)A[i]=A[i]*B[i]%mod;
    NTT(A,-1);
    for(int i=0;i<=m;++i)printf("%d ",A[i]);
    return 0;
}
