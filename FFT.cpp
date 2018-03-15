#include<cstdio>
#include<cstring>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#define ll long long
#define max(a,b) a>b?a:b
#define min(a,b) a<b?a:b
using namespace std;
inline int read(){
	int x=0,f=1;char ch=' ';
	while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while(ch>='0'&&ch<='9')x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
	return x*f;
}
const int N=3e6;
const double pi=acos(-1);
struct cp{
    double r,i;
    cp(){}
    cp(double _r,double _i):r(_r),i(_i){}
    inline cp operator + (const cp& b) const {return cp(r+b.r,i+b.i);}
    inline cp operator - (const cp& b) const {return cp(r-b.r,i-b.i);}
    inline cp operator * (const cp& b) const {return cp(r*b.r-i*b.i,r*b.i+b.r*i);}
}a[N],b[N];
int n,m,L,R[N],ans[N],Len;
inline void FFT(cp* a,int n,int f){
    for(int i=0;i<n;++i)R[i]=(R[i>>1]>>1)|((i&1)<<(Len-1));
    for(int i=0;i<n;++i)if(i<R[i])swap(a[i],a[R[i]]);
    for(int i=1;i<n;i<<=1){
        cp wn(cos(pi/i),f*sin(pi/i));
        for(int j=0;j<n;j+=(i<<1)){
            cp w(1,0);
            for(int k=0;k<i;++k,w=w*wn){
                cp x=a[j+k],y=w*a[j+k+i];
                a[j+k]=x+y;a[j+k+i]=x-y;
            }
        }
    }
    if(f==-1)for(int i=0;i<n;++i)a[i].r/=n;
}
int main(){
    n=read();m=read();
    for(int i=0;i<=n;++i)a[i].r=read();
    for(int i=0;i<=m;++i)b[i].r=read();
    for(L=1;L<n+m+1;L<<=1,Len++);
    FFT(a,L,1);
    FFT(b,L,1);
    for(int i=0;i<=L;++i)a[i]=a[i]*b[i];
    FFT(a,L,-1);
    for(int i=0;i<n+m+1;++i)ans[i]=(int)(a[i].r+0.3);
    for(int i=0;i<n+m+1;++i)printf("%d ",ans[i]);
    return 0;
}
