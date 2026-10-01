// generic-n invariants for raw tables (uint16/uint32 per point, 2^n entries): APN check, ODS, radical dist, maxsub (target n/2)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <cstring>
#include <algorithm>
int n,N; std::vector<uint32_t> G; std::vector<uint32_t> pi_;
static int rk(std::vector<uint32_t> v){uint32_t b[32]={0};int r=0;for(uint32_t x:v){while(x){int h=31-__builtin_clz(x);if(!b[h]){b[h]=x;r++;break;}x^=b[h];}}return r;}
uint32_t nullvec(std::vector<uint32_t> m){int r=0;int pc[32];for(int c=0;c<n&&r<n;c++){int p=-1;for(int i=r;i<n;i++)if(m[i]>>c&1){p=i;break;}if(p<0)continue;std::swap(m[p],m[r]);for(int i=0;i<n;i++)if(i!=r&&(m[i]>>c&1))m[i]^=m[r];pc[r++]=c;}
 bool isp[32]={0};for(int i=0;i<r;i++)isp[pc[i]]=1;int fc=-1;for(int c=0;c<n;c++)if(!isp[c]){fc=c;break;} if(fc<0)return 0;uint32_t b=1u<<fc;for(int i=0;i<r;i++)if(m[i]>>fc&1)b|=1u<<pc[i];return b;}
std::vector<char> Ein; int best;
void rec(std::vector<uint32_t>&span,const std::vector<uint32_t>&c,int dim,int target){ if(dim>best)best=dim; if(best>=target)return; if(dim+(int)c.size()<=best)return;
 std::vector<char> in(N,0);
 for(size_t i=0;i<c.size();i++){ uint32_t x=c[i]; std::vector<uint32_t> ns=span; for(uint32_t s:span) ns.push_back(s^x); for(uint32_t s:ns) in[s]=1;
  std::vector<uint32_t> rest; for(size_t j=i+1;j<c.size();j++){uint32_t y=c[j]; if(in[y])continue; bool canonical=true; for(uint32_t s:ns) if((y^s)<y){canonical=false;break;} if(!canonical)continue; bool ok=true; for(uint32_t s:ns) if(!Ein[s^y]){ok=false;break;} if(ok)rest.push_back(y);}
  for(uint32_t s:ns) in[s]=0; rec(ns,rest,dim+1,target); if(best>=target)return; } }
int main(int argc,char**argv){ n=atoi(argv[1]); N=1<<n; int dosub=argc>3?atoi(argv[3]):1;
 FILE*f=fopen(argv[2],"rb"); G.resize(N); if(fread(G.data(),4,N,f)!=(size_t)N) return 1; fclose(f);
 pi_.assign(N,0); bool apn=true;
 for(int a=1;a<N;a++){ std::vector<uint32_t> cols(n); for(int i=0;i<n;i++) cols[i]=G[a^(1<<i)]^G[a]^G[1<<i]^G[0]; if(rk(cols)!=n-1){apn=false;break;} pi_[a]=nullvec(cols);}
 printf("APN=%d",apn); if(!apn){printf("\n");return 0;}
 std::map<int,long> ods; std::vector<int> cnt(N);
 for(int a=1;a<N;a++){ std::fill(cnt.begin(),cnt.end(),0); for(int x=0;x<N;x++) cnt[pi_[x]^pi_[x^a]]++; for(int b=0;b<N;b++) ods[cnt[b]]++; }
 printf(" ods{"); for(auto&p:ods) printf("%d:%ld ",p.first,p.second); printf("}");
 std::vector<std::vector<uint32_t>> pol(n,std::vector<uint32_t>(n)); for(int i=0;i<n;i++)for(int j=0;j<n;j++) pol[i][j]=G[(1<<i)^(1<<j)]^G[1<<i]^G[1<<j]^G[0];
 std::map<int,int> rad; Ein.assign(N,0); Ein[0]=1; std::vector<uint32_t> E;
 for(int b=1;b<N;b++){ std::vector<uint32_t> rows(n,0); for(int i=0;i<n;i++)for(int j=0;j<n;j++) if(__builtin_parity(pol[i][j]&b)) rows[i]|=1u<<j; int r=n-rk(rows); rad[r]++; if(r>0){Ein[b]=1;E.push_back(b);} }
 printf(" rad{"); for(auto&p:rad) printf("%d:%d ",p.first,p.second); printf("}");
 if(dosub){ best=0; std::vector<uint32_t> sp{0}; rec(sp,E,0,n/2); printf(" maxsub%s%d (need %d)", best>=n/2?">=":"=",best,n/2);} printf("\n"); }
