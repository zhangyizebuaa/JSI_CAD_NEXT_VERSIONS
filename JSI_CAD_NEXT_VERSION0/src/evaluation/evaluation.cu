#include <cstdio>
#include <cstdlib>

#include "evaluation.cuh"
#include "jsi_cad.hpp"
#include "common/common.cuh"
#include "common/profiler.hpp"

#define MAXDEG 25
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

static float BinomialCoefficient(int i, int j) {
#define MAX_HALF_N 26
  static const float bc[((MAX_HALF_N - 2) * (MAX_HALF_N - 1)) / 2 + MAX_HALF_N - 2] = {
      15.0,      20.0,       28.0,       56.0,       70.0,       45.0,       120.0,      210.0,
      252.0,     66.0,       220.0,      495.0,      792.0,      924.0,      91.0,       364.0,
      1001.0,    2002.0,     3003.0,     3432.0,     120.0,      560.0,      1820.0,     4368.0,
      8008.0,    11440.0,    12870.0,    153.0,      816.0,      3060.0,     8568.0,     18564.0,
      31824.0,   43758.0,    48620.0,    190.0,      1140.0,     4845.0,     15504.0,    38760.0,
      77520.0,   125970.0,   167960.0,   184756.0,   231.0,      1540.0,     7315.0,     26334.0,
      74613.0,   170544.0,   319770.0,   497420.0,   646646.0,   705432.0,   276.0,      2024.0,
      10626.0,   42504.0,    134596.0,   346104.0,   735471.0,   1307504.0,  1961256.0,  2496144.0,
      2704156.0, 325.0,      2600.0,     14950.0,    65780.0,    230230.0,   657800.0,   1562275.0,
      3124550.0, 5311735.0,  7726160.0,  9657700.0,  10400600.0, 378.0,      3276.0,     20475.0,
      98280.0,   376740.0,   1184040.0,  3108105.0,  6906900.0,  13123110.0, 21474180.0, 30421755.0,
      37442160.0, 40116600.0, 435.0,      4060.0,     27405.0,    142506.0,   593775.0,   2035800.0,
      5852925.0, 14307150.0, 30045015.0, 54627300.0, 86493225.0, 119759850.0,145422675.0,155117520.0,
      496.0,     4960.0,     35960.0,    201376.0,   906192.0,   3365856.0,  10518300.0, 28048800.0,
      64512240.0,129024480.0,225792840.0,347373600.0,471435600.0,565722720.0,601080390.0,561.0,
      5984.0,    46376.0,    278256.0,   1344904.0,  5379616.0,  18156204.0, 52451256.0, 131128140.0,
      286097760.0,548354040.0,927983760.0,1391975640.0,1855967520.0,2203961430.0,2333606220.0,630.0,
      7140.0,    58905.0,    376992.0,   1947792.0,  8347680.0,  30260340.0, 94143280.0, 254186856.0,
      600805296.0,1251677700.0,2310789600.0,3796297200.0,5567902560.0,7307872110.0,8597496600.0,9075135300.0,
      703.0,     8436.0,     73815.0,    501942.0,   2760681.0,  12620256.0, 48903492.0, 163011640.0,
      472733756.0,1203322288.0,2707475148.0,5414950296.0,9669554100.0,15471286560.0,22239974430.0,28781143380.0,
      33578000610.0,35345263800.0,780.0, 9880.0, 91390.0, 658008.0, 3838380.0, 18643560.0, 76904685.0, 273438880.0,
      847660528.0,2311801440.0,5586853480.0,12033222880.0,23206929840.0,40225345056.0,62852101650.0,88732378800.0,
      113380261800.0,131282408400.0,137846528820.0,861.0,11480.0,111930.0,850668.0,5245786.0,26978328.0,118030185.0,
      445891810.0,1471442973.0,4280561376.0,11058116888.0,25518731280.0,52860229080.0,98672427616.0,166509721602.0,
      254661927156.0,353697121050.0,446775310800.0,513791607420.0,538257874440.0,946.0,13244.0,135751.0,1086008.0,
      7059052.0,38320568.0,177232627.0,708930508.0,2481256778.0,7669339132.0,21090682613.0,51915526432.0,114955808528.0,
      229911617056.0,416714805914.0,686353797976.0,1029530696964.0,1408831480056.0,1761039350070.0,2012616400080.0,
      2104098963720.0,1035.0,15180.0,163185.0,1370754.0,9366819.0,53524680.0,260932815.0,1101716330.0,4076350421.0,
      13340783196.0,38910617655.0,101766230790.0,239877544005.0,511738760544.0,991493848554.0,1749695026860.0,
      2818953098830.0,4154246671960.0,5608233007146.0,6943526580276.0,7890371113950.0,8233430727600.0,1128.0,17296.0,
      194580.0,1712304.0,12271512.0,73629072.0,377348994.0,1677106640.0,6540715896.0,22595200368.0,69668534468.0,
      192928249296.0,482320623240.0,1093260079344.0,2254848913647.0,4244421484512.0,7309837001104.0,11541847896480.0,
      16735679449896.0,22314239266528.0,27385657281648.0,30957699535776.0,32247603683100.0,1225.0,19600.0,230300.0,
      2118760.0,15890700.0,99884400.0,536878650.0,2505433700.0,10272278170.0,37353738800.0,121399651100.0,354860518600.0,
      937845656300.0,2250829575120.0,4923689695575.0,9847379391150.0,18053528883775.0,30405943383200.0,47129212243960.0,
      67327446062800.0,88749815264600.0,108043253365600.0,121548660036300.0,126410606437752.0,1326.0,22100.0,270725.0,
      2598960.0,20358520.0,133784560.0,752538150.0,3679075400.0,15820024220.0,60403728840.0,206379406870.0,635013559600.0,
      1768966344600.0,4481381406320.0,10363194502115.0,21945588357420.0,42671977361650.0,76360380541900.0,125994627894135.0,
      191991813933920.0,270533919634160.0,352870329957600.0,426384982032100.0,477551179875952.0,495918532948104.0};

  int n, half_n, bc_i;
  if (i < 0 || j < 0) return 0.0f;
  if (0 == i || 0 == j) return 1.0f;
  n = i + j;
  if (1 == i || 1 == j) return static_cast<float>(n);
  if (4 == n) return 6.0f;
  if (5 == n) return 10.0f;
  if (n % 2) return BinomialCoefficient(i - 1, j) + BinomialCoefficient(i, j - 1);
  half_n = n >> 1;
  if (half_n > MAX_HALF_N) return BinomialCoefficient(i - 1, j) + BinomialCoefficient(i, j - 1);
  if (i > half_n) i = n - i;
  half_n -= 2;
  bc_i = ((half_n * (half_n + 1)) >> 1) + i - 3;
  return bc[bc_i];
#undef MAX_HALF_N
}

static void EG_EvaluateQuotientRule2(int dim, int der_count, int v_stride, float *v) {
  float F, Fs, Ft, ws, wt, wss, wtt, wst, *f, *x;
  int i, j, n, q, ii, jj, Fn;
  // Same two-step reciprocal as FORMER `EG_EvaluateQuotientRule2` __device__ version.
  F = v[dim];
  F = 1.0f / F;
  if (v_stride > dim + 1) {
    i = ((der_count + 1) * (der_count + 2) >> 1);
    x = v;
    j = dim + 1;
    q = v_stride - j;
    while (i--) {
      jj = j;
      while (jj--) *x++ *= F;
      x += q;
    }
  } else {
    i = (((der_count + 1) * (der_count + 2)) >> 1) * v_stride;
    x = v;
    while (i--) *x++ *= F;
  }

  if (!der_count) return;
  f = v;
  x = v + v_stride;
  ws = -x[dim];
  wt = -x[dim + v_stride];
  j = dim;
  while (j--) {
    F = *f++;
    *x += ws * F;
    x[v_stride] += wt * F;
    x++;
  }
  if (der_count <= 1) return;

  f += (v_stride - dim);
  x = v + 3 * v_stride;
  wss = -x[dim];
  wst = -x[v_stride + dim];
  n = 2 * v_stride;
  wtt = -x[n + dim];
  j = dim;
  while (j--) {
    F = *v++;
    Ft = f[v_stride];
    Fs = *f++;
    *x += wss * F + 2.0f * ws * Fs;
    x[v_stride] += wst * F + wt * Fs + ws * Ft;
    x[n] += wtt * F + 2.0f * wt * Ft;
    x++;
  }

  if (der_count <= 2) return;
  v -= dim;
  f = v + 6 * v_stride;
  for (n = 3; n <= der_count; n++) {
    for (j = 0; j <= n; j++) {
      i = n - j;
      for (ii = 0; ii <= i; ii++) {
        ws = BinomialCoefficient(ii, i - ii);
        for (jj = ii ? 0 : 1; jj <= j; jj++) {
          q = ii + jj;
          Fn = ((q * (q + 1)) / 2 + jj) * v_stride + dim;
          wt = -ws * BinomialCoefficient(jj, j - jj) * v[Fn];
          q = n - q;
          Fn = ((q * (q + 1)) / 2 + j - jj) * v_stride;
          for (q = 0; q < dim; q++) f[q] += wt * v[Fn + q];
        }
      }
      f += v_stride;
    }
  }
}

static int FindSpan(int nKnots, int degree, float u, float *U) {
  int n, low, mid, high;
  if (u <= U[degree]) return degree;
  n = nKnots - degree - 1;
  if (u >= U[n]) return n - 1;
  low = degree;
  high = n;
  mid = (low + high) / 2;
  while ((u < U[mid]) || (u >= U[mid + 1])) {
    if (u < U[mid]) high = mid;
    else low = mid;
    mid = (low + high) / 2;
  }
  return mid;
}

static void DersBasisFuns(int i, int p, float u, float *knot, int der, float **ders) {
  int j, k, j1, j2, r, s1, s2, rk, pk;
  float d, saved, temp, ndu[MAXDEG + 1][MAXDEG + 1];
  float a[2][MAXDEG + 1], left[MAXDEG + 1], right[MAXDEG + 1];
  ndu[0][0] = 1.0f;
  for (j = 1; j <= p; j++) {
    left[j] = u - knot[i + 1 - j];
    right[j] = knot[i + j] - u;
    saved = 0.0f;
    for (r = 0; r < j; r++) {
      ndu[j][r] = right[r + 1] + left[j - r];
      temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }
  for (j = 0; j <= p; j++) ders[0][j] = ndu[j][p];
  for (r = 0; r <= p; r++) {
    s1 = 0;
    s2 = 1;
    a[0][0] = 1.0f;
    for (k = 1; k <= der; k++) {
      d = 0.0f;
      rk = r - k;
      pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
        d = a[s2][0] * ndu[rk][pk];
      }
      j1 = rk >= -1 ? 1 : -rk;
      j2 = (r - 1 <= pk) ? k - 1 : p - r;
      for (j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
        d += a[s2][j] * ndu[rk + j][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
        d += a[s2][k] * ndu[r][pk];
      }
      ders[k][r] = d;
      j = s1;
      s1 = s2;
      s2 = j;
    }
  }
  r = p;
  for (k = 1; k <= der; k++) {
    for (j = 0; j <= p; j++) ders[k][j] *= r;
    r *= p - k;
  }
}

struct EvalContext {
  int flags;
  int flat;
  int rat;
  int degu;
  int degv;
  int nCPu;
  int nKu;
  int nKv;
  int du;
  int dv;
  int der;
  float *U;
  float *Kv;
  float *CP;
  float *w;
};

static EvalContext BuildEvalContext(int *ivec, float *data) {
  EvalContext ctx;
  ctx.flags = ivec[0];
  ctx.flat = ivec[0] & 1;
  ctx.rat = ivec[0] & 2;
  ctx.degu = ivec[1];
  ctx.nCPu = ivec[2];
  ctx.nKu = ivec[3];
  ctx.degv = ivec[4];
  ctx.nKv = ivec[6];
  ctx.der = 2;
  ctx.du = MIN(ctx.der, ctx.degu);
  ctx.dv = MIN(ctx.der, ctx.degv);
  ctx.U = data;
  ctx.Kv = data + ivec[3];
  ctx.CP = data + ivec[3] + ivec[6];
  ctx.w = ctx.CP + 3 * ivec[2] * ivec[5];
  return ctx;
}

static void EvaluateOneUV(const EvalContext &ctx, float U, float V, float *deriv) {
  int i, j, k, l, m, s, spanu, spanv;
  float *NderU[MAXDEG + 1], *NderV[MAXDEG + 1];
  float Nu[MAXDEG + 1][MAXDEG + 1], Nv[MAXDEG + 1][MAXDEG + 1], temp[4 * MAXDEG];
  float v[24];

  for (m = l = 0; l <= ctx.der; l++) {
    for (k = 0; k <= ctx.der - l; k++, m++) {
      deriv[3 * m] = deriv[3 * m + 1] = deriv[3 * m + 2] = 0.0f;
      v[4 * m] = v[4 * m + 1] = v[4 * m + 2] = v[4 * m + 3] = 0.0f;
    }
  }
  if (ctx.flat == 0 && (ctx.flags & 12) != 0) return;
  if (ctx.degu >= MAXDEG || ctx.degv >= MAXDEG) return;
  for (i = 0; i <= ctx.degu; i++) NderU[i] = &Nu[i][0];
  for (i = 0; i <= ctx.degv; i++) NderV[i] = &Nv[i][0];
  spanu = FindSpan(ctx.nKu, ctx.degu, U, ctx.U);
  DersBasisFuns(spanu, ctx.degu, U, ctx.U, ctx.du, NderU);
  spanv = FindSpan(ctx.nKv, ctx.degv, V, ctx.Kv);
  DersBasisFuns(spanv, ctx.degv, V, ctx.Kv, ctx.dv, NderV);

  if (ctx.rat == 0) {
    for (m = l = 0; l <= ctx.dv; l++) {
      for (k = 0; k <= ctx.der - l; k++, m++) {
        if (k > ctx.du) continue;
        for (s = 0; s <= ctx.degv; s++) {
          temp[3 * s] = temp[3 * s + 1] = temp[3 * s + 2] = 0.0f;
          int cp_row = ctx.nCPu * (spanv - ctx.degv + s);
          for (j = 0; j <= ctx.degu; j++) {
            i = spanu - ctx.degu + j + cp_row;
            temp[3 * s] += Nu[k][j] * ctx.CP[3 * i];
            temp[3 * s + 1] += Nu[k][j] * ctx.CP[3 * i + 1];
            temp[3 * s + 2] += Nu[k][j] * ctx.CP[3 * i + 2];
          }
        }
        for (s = 0; s <= ctx.degv; s++) {
          deriv[3 * m] += Nv[l][s] * temp[3 * s];
          deriv[3 * m + 1] += Nv[l][s] * temp[3 * s + 1];
          deriv[3 * m + 2] += Nv[l][s] * temp[3 * s + 2];
        }
      }
    }
    temp[0] = deriv[6];
    temp[1] = deriv[7];
    temp[2] = deriv[8];
    deriv[6] = deriv[9];
    deriv[7] = deriv[10];
    deriv[8] = deriv[11];
    deriv[9] = temp[0];
    deriv[10] = temp[1];
    deriv[11] = temp[2];
    return;
  }

  for (m = l = 0; l <= ctx.dv; l++) {
    for (k = 0; k <= ctx.der - l; k++, m++) {
      if (k > ctx.du) continue;
      for (s = 0; s <= ctx.degv; s++) {
        temp[4 * s] = temp[4 * s + 1] = temp[4 * s + 2] = temp[4 * s + 3] = 0.0f;
        int cp_row = ctx.nCPu * (spanv - ctx.degv + s);
        for (j = 0; j <= ctx.degu; j++) {
          i = spanu - ctx.degu + j + cp_row;
          temp[4 * s] += Nu[k][j] * ctx.w[i] * ctx.CP[3 * i];
          temp[4 * s + 1] += Nu[k][j] * ctx.w[i] * ctx.CP[3 * i + 1];
          temp[4 * s + 2] += Nu[k][j] * ctx.w[i] * ctx.CP[3 * i + 2];
          temp[4 * s + 3] += Nu[k][j] * ctx.w[i];
        }
      }
      for (s = 0; s <= ctx.degv; s++) {
        v[4 * m] += Nv[l][s] * temp[4 * s];
        v[4 * m + 1] += Nv[l][s] * temp[4 * s + 1];
        v[4 * m + 2] += Nv[l][s] * temp[4 * s + 2];
        v[4 * m + 3] += Nv[l][s] * temp[4 * s + 3];
      }
    }
  }

  temp[0] = v[8];
  temp[1] = v[9];
  temp[2] = v[10];
  temp[3] = v[11];
  v[8] = v[12];
  v[9] = v[13];
  v[10] = v[14];
  v[11] = v[15];
  v[12] = temp[0];
  v[13] = temp[1];
  v[14] = temp[2];
  v[15] = temp[3];
  EG_EvaluateQuotientRule2(3, ctx.der, 4, v);
  for (m = l = 0; l <= ctx.der; l++) {
    for (k = 0; k <= ctx.der - l; k++, m++) {
      deriv[3 * m] = v[4 * m];
      deriv[3 * m + 1] = v[4 * m + 1];
      deriv[3 * m + 2] = v[4 * m + 2];
    }
  }
}

void EG_spline2dDeriv_irrational(int *ivec, float *data, float ustep, float vstep, float ustart, float vstart,
                                 int nuv, float *results) {
  CAD_PROFILE_SCOPE("EG_spline2dDeriv_irrational");
  EvalContext ctx = BuildEvalContext(ivec, data);
  for (int idy = 0; idy < nuv; ++idy) {
    for (int idx = 0; idx < nuv; ++idx) {
      int offset = idx + idy * nuv;
      EvaluateOneUV(ctx, ustep * idx + ustart, vstep * idy + vstart, results + offset * 18);
    }
  }
}

void EG_spline2dDeriv_irrational_opt(int *ivec, float *data, float ustep, float vstep, float ustart, float vstart,
                                     int nuv, int ystart, int yend, float *results) {
  CAD_PROFILE_SCOPE("EG_spline2dDeriv_irrational_opt");
  EvalContext ctx = BuildEvalContext(ivec, data);
  int y0 = ystart < 0 ? 0 : ystart;
  int y1 = yend > nuv ? nuv : yend;
  for (int idy = y0; idy < y1; ++idy) {
    for (int idx = 0; idx < nuv; ++idx) {
      int offset = idx + idy * nuv;
      EvaluateOneUV(ctx, ustep * idx + ustart, vstep * idy + vstart, results + offset * 18);
    }
  }
}

void evalWithCuda(int *ivec, float *data, int ndata, float ustep, float vstep, float ustart, float vstart, int nuv,
                  float *results) {
  CAD_PROFILE_SCOPE("evalWithCuda");
  (void)ndata;
  // FORMER `evalWithCuda`: PIPELINE_SIZE strips along v, each strip calls `EG_spline2dDeriv_irrational_opt` on device.
#define PIPELINE_SIZE 4
  int pipeline_split_size = (nuv + PIPELINE_SIZE - 1) / PIPELINE_SIZE;
  for (int i = 0; i < PIPELINE_SIZE; i++) {
    int ystart = pipeline_split_size * i;
    int yend = MIN(nuv, pipeline_split_size * (i + 1));
    if (ystart >= yend) break;
    EG_spline2dDeriv_irrational_opt(ivec, data, ustep, vstep, ustart, vstart, nuv, ystart, yend, results);
  }
#undef PIPELINE_SIZE
}
