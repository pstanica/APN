// All quadratic switching neighbours F + v*Q (n<=12): for each v, nullspace of linear system in the n(n-1)/2
// coefficients of B_Q; reports nullity and writes one representative per coset of (component span) for each v,
// enumerating up to 2^min(excess,6) cosets.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <cstring>
typedef unsigned __int128 u128;
int n,N,K; std::vector<uint32_t> G,pi_; std::vector<uint32_t> sol; int pidx[16][16];
static int rk(std::vector<uint32_t> v){uint32_t b[32]={0};int r=0;for(uint32_t x:v){while(x){int h=31-__builtin_clz(x);if(!b[h]){b[h]=x;r++;break;}x^=b[h];}}return r;}
int main(int argc,char**argv){ n=atoi(argv[1]); N=1<<n; std::string pre=argv[3]; int vstride=argc>4?atoi(argv[4]):1;
 FILE*f=fopen(argv[2],"rb"); G.resize(N); if(fread(G.data(),4,N,f)!=(size_t)N)return 1; fclose(f); for(int i=N-1;i>=0;i--) G[i]^=G[0];
 K=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){pidx[i][j]=pidx[j][i]=K++;}
 pi_.assign(N,0); sol.assign((size_t)N*N,0);
 for(int a=1;a<N;a++){ std::vector<uint32_t> cols(n); for(int i=0;i<n;i++) cols[i]=G[a^(1<<i)]^G[a]^G[1<<i];
   if(rk(cols)!=n-1){fprintf(stderr,"not APN\n");return 1;}
   for(uint32_t b=1;b<(uint32_t)N;b++){bool ok=true;for(auto c:cols) if(__builtin_parity(c&b)){ok=false;break;} if(ok){pi_[a]=b;break;}}
   for(int x=0;x<N;x++) sol[(size_t)a*N+(G[a^x]^G[a]^G[x])]=x; }
 auto formOfBeta=[&](uint32_t beta){u128 r=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){uint32_t p=G[(1<<i)^(1<<j)]^G[1<<i]^G[1<<j]; if(__builtin_parity(p&beta)) r|=(u128)1<<pidx[i][j];} return r;};
 long saved=0; std::vector<long> hist(K+1,0);
 for(int v=1;v<N;v+=vstride){
   u128 basis[128]; memset(basis,0,sizeof basis); int rank=0;
   for(int a=1;a<N;a++){ if(__builtin_parity(pi_[a]&v)) continue; int x=sol[(size_t)a*N+v]; u128 row=0;
     for(int i=0;i<n;i++)for(int j=i+1;j<n;j++) if(((a>>i&1)&(x>>j&1))^((a>>j&1)&(x>>i&1))) row|=(u128)1<<pidx[i][j];
     for(int h=K-1;h>=0&&row;h--) if(row>>h&1){ if(basis[h]) row^=basis[h]; else {basis[h]=row;rank++;break;} } }
   int nul=K-rank; hist[nul]++; int excess=nul-(n-1); if(excess<=0) continue;
   for(int h=0;h<K;h++) if(basis[h]) for(int g=h+1;g<K;g++) if(basis[g]>>h&1) basis[g]^=basis[h];
   std::vector<u128> null; for(int fv=0;fv<K;fv++){ if(basis[fv]) continue; u128 vec=(u128)1<<fv; for(int h=0;h<K;h++) if(basis[h]&&(basis[h]>>fv&1)) vec|=(u128)1<<h; null.push_back(vec);}
   // component span
   u128 cb[128]; memset(cb,0,sizeof cb); for(uint32_t beta=1;beta<(uint32_t)N;beta++) if(!__builtin_parity(beta&v)){u128 r=formOfBeta(beta); for(int h=K-1;h>=0&&r;h--) if(r>>h&1){ if(cb[h]) r^=cb[h]; else {cb[h]=r;break;} } }
   // quotient basis
   std::vector<u128> qb; for(u128 q:null){ u128 r=q; for(int h=K-1;h>=0&&r;h--) if(r>>h&1){ if(cb[h]) r^=cb[h]; else {cb[h]=r; qb.push_back(q); break;} } }
   int e=qb.size(); int lim=e<6?e:6;
   for(int s=1;s<(1<<lim);s++){ u128 q=0; for(int t=0;t<lim;t++) if(s>>t&1) q^=qb[t];
     std::vector<uint32_t> out(N); for(int x=0;x<N;x++){int qv=0; for(int i=0;i<n;i++) if(x>>i&1) for(int j=i+1;j<n;j++) if(x>>j&1) qv^=(int)(q>>pidx[i][j]&1); out[x]=G[x]^(qv?v:0);}
     char fn[256]; sprintf(fn,"%s_v%d_s%d.raw",pre.c_str(),v,s); FILE*o=fopen(fn,"wb"); fwrite(out.data(),4,N,o); fclose(o); saved++; }
 }
 fprintf(stderr,"nullity histogram:"); for(int k=0;k<=K;k++) if(hist[k]) fprintf(stderr," %d:%ld",k,hist[k]); fprintf(stderr,"  saved %ld\n",saved);
}
