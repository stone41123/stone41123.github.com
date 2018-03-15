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
inline void poly_inv(LL *a,LL *b,int n){
    if(n==1){b[0]=ksm(a[0],mod-2);return;}
    poly_inv(a,b,(n+1)>>1);
    int m=1;
    for(;m<(n<<1);m<<=1);
    for(int i=0;i<n;++i)tmp[i]=a[i];
    for(int i=n;i<m;++i)tmp[i]=0;
    NTT(tmp,m,1);NTT(b,m,1);
    for(int i=0;i<m;++i)
        b[i]=(2LL*b[i]%mod-tmp[i]*b[i]%mod*b[i]%mod+mod)%mod;
    NTT(b,m,-1);
    for(int i=n;i<m;++i)b[i]=0;
    return;
}
LL dao_f[N],inv_f[N],ln_f[N],f0[N],inv[N];
inline void poly_ln(LL *f,LL *g,int n){
	int m=1;
	for(;m<(n<<1);m<<=1);
	for(int i=0;i<n;++i)dao_f[i]=f[i+1]*(i+1)%mod;
	for(int i=n;i<m;++i)dao_f[i]=0;
	for(int i=0;i<m;++i)inv_f[i]=0;
	poly_inv(f,inv_f,n);
	NTT(dao_f,m,1);NTT(inv_f,m,1);
	for(int i=0;i<m;++i)g[i]=dao_f[i]*inv_f[i]%mod;
	NTT(g,m,-1);
	for(int i=n-1;i>0;--i)g[i]=g[i-1]*inv[i]%mod;g[0]=0;
	for(int i=n;i<m;++i)g[i]=0;
}
inline void poly_exp(LL *f,LL *g,int n){
	if(n==1){g[0]=1;return;}
	poly_exp(f,g,(n+1)>>1);
	int m=1;
	for(;m<(n<<1);m<<=1);
	poly_ln(g,ln_f,n);
	for(int i=0;i<n;++i)tmp[i]=((i==0)-ln_f[i]+f[i]+mod)%mod;
	for(int i=n;i<m;++i)tmp[i]=0;
	NTT(tmp,m,1);NTT(g,m,1);
	for(int i=0;i<m;++i)
		g[i]=g[i]*tmp[i]%mod;
	NTT(g,m,-1);
	for(int i=n;i<m;++i)g[i]=0;
	return;
}
LL num[N],f[N],ans[N];
int main(){
	// freopen("./data/data10.in","r",stdin);
	// freopen("./data/data10.out","w",stdout);
    for(int i=0;i<=19;++i)bin[1<<i]=i;
    n=read();int mx=read()+1;inv[0]=inv[1]=1;
	for(int i=2;i<mx;++i)inv[i]=(mod-mod/i)*inv[mod%i]%mod;
    for(int i=1;i<=n;++i)num[read()]++;
	for(m=1;m<mx;m<<=1)L++;
	for(int i=1;i<mx;++i)if(num[i])
		for(int j=1;i*j<mx;++j)
			f[i*j]=(f[i*j]+num[i]*inv[j]%mod)%mod;
	poly_exp(f,ans,mx);
	for(int i=1;i<mx;++i)printf("%lld\n",(ans[i]+mod)%mod);
    return 0;
}
