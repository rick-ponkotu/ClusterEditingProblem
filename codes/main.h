#define FOR(i,s,t) for(int i=s; i < t; i++)
#define REP(i,n) FOR(i,0,n)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SWAP(type,a,b) do{type tmp=a;a=b;b=tmp;} while(0)

#define dump(x) cerr<<#x<<"="<<x<<", "
#define dumpl(x) cerr<<#x<<"="<<x<<endl
#define dump_vec(a)  for(const auto& e : a) {cerr << e << std::endl; }

typedef std::pair<int, int> edge;

#define DEBUG 1  // 0 -> no output, 1 -> debug
#define RANDOM_SEED 1243
#define NUM_PIVOT 5

void show_sol(const std::vector <edge>& sol);
void write_sol(const std::vector <edge>& sol, const char *fo);
