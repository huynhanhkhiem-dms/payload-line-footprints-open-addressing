/* hashcost.c -- probe-model vs memory-hierarchy cost of open addressing.
 *
 * Schemes:
 *   up  : independent random-slot probing (with replacement)
 *   dh  : double hashing (full-cycle permutation for power-of-two n)
 *   lp  : linear probing
 *   wtl : walking-the-line baseline (Laarman-style line-local scan, then rehash)
 *   bp  : blocked/bucketized probing (group scan, linear over groups)
 *   eh  : elastic hashing   (Farach-Colton, Krapivin, Kuszmaul, FOCS'24)
 *   fh  : funnel hashing    (idem)
 *
 * Metrics per operation: probes (slot reads), distinct 64B cache lines,
 * distinct 4KiB pages.  Slots are 8 bytes (uint64 key), EMPTY == 0.
 *
 * Build with -DINSTRUMENT for counting; without it for timing.
 * The `up` baseline is deliberately not called classical uniform hashing: FKK
 * define uniform probing as a random permutation, whereas `up` samples slots
 * independently with replacement.  `dh` supplies a no-repeat scattered baseline.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
typedef uint32_t u32;

#define EMPTY 0ULL
#define SLOTB 8            /* bytes per slot */
#define LINEB 64
#define PAGEB 4096
#define SLOTS_PER_LINE (LINEB / SLOTB)   /* 8 */
#define SLOTS_PER_PAGE (PAGEB / SLOTB)   /* 512 */

/* ---------------- hashing ---------------- */
static inline u64 splitmix64(u64 x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
static u64 RUN_SEED = 0;
static inline u64 hash2(u64 x, u64 s) {
    return splitmix64((x ^ RUN_SEED) ^ splitmix64(s + RUN_SEED * 0xD1B54A32D192ED03ULL));
}

/* Multiply-high range reduction: value in [0,m). For a uniform 64-bit input,
 * bucket preimage counts differ by at most one. At m <= 2^24 the resulting
 * per-bucket probability imbalance is below 2^-64 in absolute probability and
 * is negligible relative to sampling variability in this experiment. */
static inline u64 red(u64 h, u64 m) { return (u64)(((__uint128_t)h * (__uint128_t)m) >> 64); }


/* Exact ideal-model expectation for distinct b-slot payload lines touched by
 * an unsuccessful independent random-slot query.  The e empty positions are
 * averaged over a uniformly random e-subset of the n slots. */
static double up_expected_miss_lines(u64 n, u64 e, u64 b) {
    if (!e || !b || n % b) return NAN;
    long double p = 1.0L; /* P(R=0), R empties in one b-slot line */
    for (u64 t = 0; t < e; t++) p *= (long double)(n - b - t) / (long double)(n - t);
    long double ans = 0.0L;
    for (u64 r = 0; r <= b; r++) {
        if (r <= e && e - r <= n - b)
            ans += p * (long double)b / (long double)(b + e - r);
        if (r == b) break;
        long double num1 = (long double)(b - r), den1 = (long double)(r + 1);
        long double num2 = (long double)(e - r);
        long double den2 = (long double)(n - b - e + r + 1);
        if (den2 <= 0 || num2 <= 0) p = 0.0L;
        else p *= (num1 / den1) * (num2 / den2);
    }
    return (double)((long double)(n / b) * ans);
}

/* ---------------- per-operation instrumentation ---------------- */
typedef struct {
    u64 probes;
    u64 lines;
    u64 pages;
} OpCost;

#ifdef INSTRUMENT
/* Per-operation dedup uses generation-stamped open-addressed sets. This keeps
 * counting exact without the quadratic linear scan that would otherwise
 * dominate long unsuccessful uniform-probing queries. */
#define HSET (1u << 20)
static u64 lkey[HSET], pkey[HSET];
static u32 lstamp[HSET], pstamp[HSET], hgen = 1;
static u64 cur_probes, cur_lines, cur_pages;
static int INSTRUMENT_ENABLED = 1;

static inline u64 hset_hash(u64 x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}
static void op_begin(void) {
    cur_probes = cur_lines = cur_pages = 0;
    if (++hgen == 0) { memset(lstamp,0,sizeof lstamp); memset(pstamp,0,sizeof pstamp); hgen = 1; }
}
static inline void add_line(u64 L) {
    u64 i = hset_hash(L) & (HSET - 1), start = i;
    for (;;) {
        if (lstamp[i] != hgen) { lstamp[i] = hgen; lkey[i] = L; cur_lines++; return; }
        if (lkey[i] == L) return;
        i = (i + 1) & (HSET - 1);
        if (i == start) { fprintf(stderr, "line-dedup set exhausted; increase HSET\n"); exit(4); }
    }
}
static inline void add_page(u64 P) {
    u64 i = hset_hash(P) & (HSET - 1), start = i;
    for (;;) {
        if (pstamp[i] != hgen) { pstamp[i] = hgen; pkey[i] = P; cur_pages++; return; }
        if (pkey[i] == P) return;
        i = (i + 1) & (HSET - 1);
        if (i == start) { fprintf(stderr, "page-dedup set exhausted; increase HSET\n"); exit(4); }
    }
}
static inline void touch(u64 s) {
    if (!INSTRUMENT_ENABLED) return;
    cur_probes++;
    add_line(s / SLOTS_PER_LINE);
    add_page(s / SLOTS_PER_PAGE);
}
static inline void touch_run(u64 s, u64 len) {
    if (!INSTRUMENT_ENABLED || !len) return;
    cur_probes += len;
    u64 l0 = s / SLOTS_PER_LINE, l1 = (s + len - 1) / SLOTS_PER_LINE;
    for (u64 L = l0; L <= l1; L++) add_line(L);
    u64 p0 = s / SLOTS_PER_PAGE, p1 = (s + len - 1) / SLOTS_PER_PAGE;
    for (u64 P = p0; P <= p1; P++) add_page(P);
}
static void op_end(OpCost *c) { c->probes = cur_probes; c->lines = cur_lines; c->pages = cur_pages; }
#else
static inline void op_begin(void) {}
static inline void touch(u64 s) { (void)s; }
static inline void touch_run(u64 s, u64 len) { (void)s; (void)len; }
static inline void op_end(OpCost *c) { c->probes = c->lines = c->pages = 0; }
#endif

/* ---------------- statistics accumulator ---------------- */
typedef struct {
    u64 n;
    double sp, sl, spg;      /* sums */
    u64 maxp, maxl;
    u64 *hp;                 /* probe histogram for quantiles */
    u64 hcap;
} Acc;
static void acc_init(Acc *a, u64 cap) {
    memset(a, 0, sizeof *a); a->hcap = cap;
    a->hp = calloc(cap + 1, sizeof(u64));
    if (!a->hp) { fprintf(stderr, "out of memory allocating accumulator\n"); exit(3); }
}
static void acc_free(Acc *a) {
    free(a->hp);
    a->hp = NULL;
}
static void acc_add(Acc *a, OpCost c) {
    a->n++; a->sp += c.probes; a->sl += c.lines; a->spg += c.pages;
    if (c.probes > a->maxp) a->maxp = c.probes;
    if (c.lines  > a->maxl) a->maxl = c.lines;
    u64 b = c.probes < a->hcap ? c.probes : a->hcap;
    a->hp[b]++;
}
static u64 acc_q(Acc *a, double q) {
    u64 target = (u64)(q * (double)a->n), run = 0;
    for (u64 i = 0; i <= a->hcap; i++) { run += a->hp[i]; if (run >= target) return i; }
    return a->hcap;
}

/* ---------------- table ---------------- */
typedef struct {
    u64 *t;      /* slots */
    u64  n;      /* number of slots */
    double delta;
    /* elastic */
    int  k;            /* number of subarrays */
    u64  off[64], sz[64], occ[64];
    u64  cap_arr[64];  /* deepest successful local probe/bucket index (1-based) */
    int  batch;
    /* funnel */
    int  alpha, beta;
    u64  loff[64], lsz[64], lbuckets[64];
    u64  spec_off, spec_sz, Boff, Bsz, Coff, Csz;
    int  llg;          /* ceil(log2 log2 n) style give-up bound */
    int  cbucket;      /* 2 log log n */
    /* blocked */
    int  bsize;
} Table;

static void tbl_alloc(Table *T, u64 n) {
    T->t = aligned_alloc(4096, ((n * 8 + 4095) / 4096) * 4096);
    if (!T->t) { fprintf(stderr, "out of memory allocating table\n"); exit(3); }
    memset(T->t, 0, n * 8);
    T->n = n;
}

/* ============ uniform probing (i.i.d. random slots) ============ */
static int up_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    for (u64 j = 0;; j++) {
        u64 s = red(hash2(key, j * 0x1000193ULL + 7), T->n);
        touch(s);
        if (T->t[s] == EMPTY) { T->t[s] = key; op_end(c); return 1; }
        if (T->t[s] == key)   { op_end(c); return 1; }
        if (j > 200ULL * T->n) { op_end(c); return 0; }
    }
}
static int up_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    for (u64 j = 0;; j++) {
        u64 s = red(hash2(key, j * 0x1000193ULL + 7), T->n);
        touch(s);
        if (T->t[s] == key)   { op_end(c); return 1; }
        if (T->t[s] == EMPTY) { op_end(c); return 0; }
        if (j > 200ULL * T->n) { op_end(c); return 0; }
    }
}

/* ============ double hashing ============ */
static int dh_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    u64 h = hash2(key, 11), s = red(h, T->n);
    u64 step = (hash2(key, 12) | 1ULL) % T->n; if (!step) step = 1;
    for (u64 j = 0; j < T->n; j++) {
        touch(s);
        if (T->t[s] == EMPTY) { T->t[s] = key; op_end(c); return 1; }
        if (T->t[s] == key)   { op_end(c); return 1; }
        s += step; if (s >= T->n) s -= T->n;
    }
    op_end(c); return 0;
}
static int dh_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    u64 h = hash2(key, 11), s = red(h, T->n);
    u64 step = (hash2(key, 12) | 1ULL) % T->n; if (!step) step = 1;
    for (u64 j = 0; j < T->n; j++) {
        touch(s);
        if (T->t[s] == key)   { op_end(c); return 1; }
        if (T->t[s] == EMPTY) { op_end(c); return 0; }
        s += step; if (s >= T->n) s -= T->n;
    }
    op_end(c); return 0;
}

/* ============ linear probing ============ */
/* A linear-probing operation touches one circularly contiguous slot interval.
 * For this scheme only, compute the instrumentation counters analytically
 * rather than inserting every touched line/page into the generic dedup sets.
 * This is exact because n is a power of two and therefore line/page aligned. */
static inline void lp_cost(OpCost *c, u64 start, u64 len, u64 n) {
#ifdef INSTRUMENT
    c->probes = len;
    if (!len) { c->lines = c->pages = 0; return; }
    if (len >= n) {
        c->lines = n / SLOTS_PER_LINE;
        c->pages = n / SLOTS_PER_PAGE;
        return;
    }
    u64 last = start + len - 1;
    if (last < n) {
        c->lines = last / SLOTS_PER_LINE - start / SLOTS_PER_LINE + 1;
        c->pages = last / SLOTS_PER_PAGE - start / SLOTS_PER_PAGE + 1;
    } else {
        u64 end = last - n;
        c->lines = n / SLOTS_PER_LINE - start / SLOTS_PER_LINE + end / SLOTS_PER_LINE + 1;
        c->pages = n / SLOTS_PER_PAGE - start / SLOTS_PER_PAGE + end / SLOTS_PER_PAGE + 1;
    }
#else
    (void)start; (void)len; (void)n;
    c->probes = c->lines = c->pages = 0;
#endif
}
static int lp_insert(Table *T, u64 key, OpCost *c) {
    u64 start = red(hash2(key, 21), T->n), s = start;
    for (u64 j = 0; j < T->n; j++) {
        u64 len = j + 1;
        if (T->t[s] == EMPTY) { T->t[s] = key; lp_cost(c, start, len, T->n); return 1; }
        if (T->t[s] == key)   { lp_cost(c, start, len, T->n); return 1; }
        if (++s == T->n) s = 0;
    }
    lp_cost(c, start, T->n, T->n); return 0;
}
static int lp_find(Table *T, u64 key, OpCost *c) {
    u64 start = red(hash2(key, 21), T->n), s = start;
    for (u64 j = 0; j < T->n; j++) {
        u64 len = j + 1;
        if (T->t[s] == key)   { lp_cost(c, start, len, T->n); return 1; }
        if (T->t[s] == EMPTY) { lp_cost(c, start, len, T->n); return 0; }
        if (++s == T->n) s = 0;
    }
    lp_cost(c, start, T->n, T->n); return 0;
}

/* ============ walking-the-line (cache-line scan + double hashing) ============ */
/* A locality-aware baseline inspired by walking-the-line: scan every slot in
 * one payload line before advancing to another line.  Payload lines themselves
 * follow a full-cycle double-hashing permutation.  With n and the line count
 * powers of two, an odd line step visits every payload line exactly once.
 * Insertion and lookup use the same order, so the first EMPTY slot is a valid
 * absence certificate and the resulting empty-line layout is endogenous. */
static int wtl_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    const u64 B = SLOTS_PER_LINE, nl = T->n / B;
    u64 h1 = hash2(key, 41), h2 = hash2(key, 42);
    u64 line = red(h1, nl);
    u64 step = (h2 | 1ULL) & (nl - 1); if (!step) step = 1;
    u64 off = hash2(key, 43) & (B - 1);
    for (u64 j = 0; j < nl; j++) {
        u64 base = line * B;
        for (u64 r = 0; r < B; r++) {
            u64 pos = base + ((off + r) & (B - 1));
            touch(pos);
            if (T->t[pos] == EMPTY) { T->t[pos] = key; op_end(c); return 1; }
            if (T->t[pos] == key)   { op_end(c); return 1; }
        }
        line = (line + step) & (nl - 1);
    }
    op_end(c); return 0;
}
static int wtl_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    const u64 B = SLOTS_PER_LINE, nl = T->n / B;
    u64 h1 = hash2(key, 41), h2 = hash2(key, 42);
    u64 line = red(h1, nl);
    u64 step = (h2 | 1ULL) & (nl - 1); if (!step) step = 1;
    u64 off = hash2(key, 43) & (B - 1);
    for (u64 j = 0; j < nl; j++) {
        u64 base = line * B;
        for (u64 r = 0; r < B; r++) {
            u64 pos = base + ((off + r) & (B - 1));
            touch(pos);
            if (T->t[pos] == key)   { op_end(c); return 1; }
            if (T->t[pos] == EMPTY) { op_end(c); return 0; }
        }
        line = (line + step) & (nl - 1);
    }
    op_end(c); return 0;
}

/* ============ blocked probing (bucket scan, linear over buckets) ============ */
static int bp_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    u64 B = T->bsize, nb = T->n / B;
    u64 b = red(hash2(key, 31), nb);
    for (u64 j = 0; j < nb; j++) {
        u64 base = b * B;
        for (u64 i = 0; i < B; i++) {
            if (T->t[base + i] == EMPTY) { touch_run(base, i + 1); T->t[base + i] = key; op_end(c); return 1; }
            if (T->t[base + i] == key)   { touch_run(base, i + 1); op_end(c); return 1; }
        }
        touch_run(base, B);
        b++; if (b >= nb) b = 0;
    }
    op_end(c); return 0;
}
static int bp_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    u64 B = T->bsize, nb = T->n / B;
    u64 b = red(hash2(key, 31), nb);
    for (u64 j = 0; j < nb; j++) {
        u64 base = b * B;
        for (u64 i = 0; i < B; i++) {
            if (T->t[base + i] == key)  { touch_run(base, i + 1); op_end(c); return 1; }
            if (T->t[base + i] == EMPTY){ touch_run(base, i + 1); op_end(c); return 0; }
        }
        touch_run(base, B);
        b++; if (b >= nb) b = 0;
    }
    op_end(c); return 0;
}

/* ============ elastic hashing (FKK) ============ */
/* Finite realization of FKK Section II.  For the power-of-two table sizes used
 * here, k=ceil(log2 n) and A_1,...,A_k satisfy |A_{i+1}|=|A_i|/2 +/- 1
 * and sum exactly to n.  Batch B_0 fills A_1 to 0.75|A_1|; batch B_i uses
 * A_i,A_{i+1} with f(eps)=C*min(log2^2(1/eps),log2(1/delta)). */
static double EH_C = 2.0;
static u64 eh_spills = 0;
static u64 eh_bleft = 0;

static void eh_setup(Table *T, double delta) {
    u64 n = T->n, tot = 0;
    int k = (int)ceil(log2((double)n));
    if (k < 2) k = 2;
    if (k > 63) k = 63;
    T->k = k;
    for (int i = 0; i < k; i++) {
        u64 z;
        if (i == k - 1) z = n - tot;
        else {
            int sh = i + 1;
            z = (sh < 64) ? (n >> sh) : 0;
            if (!z) z = 1;
        }
        T->sz[i] = z;
        T->off[i] = tot;
        T->occ[i] = 0;
        T->cap_arr[i] = 0;
        tot += z;
    }
    /* The experiments use n=2^lg; this assertion also checks the finite split. */
    if (tot != n) { fprintf(stderr, "elastic geometry does not sum to n\n"); exit(3); }
    for (int i = 0; i + 1 < k; i++) {
        double d = fabs((double)T->sz[i+1] - 0.5 * (double)T->sz[i]);
        if (d > 1.0000001) { fprintf(stderr, "elastic geometry violates half +/- 1\n"); exit(3); }
    }
    T->batch = 0;
    T->delta = delta;
    eh_spills = 0;
    eh_bleft = 0;
}
static inline u64 eh_probe(Table *T, u64 key, int i, u64 j) {
    return T->off[i] + red(hash2(key, ((u64)(i + 1) << 32) ^ (j + 1)), T->sz[i]);
}
/* f(eps) = c * min(log2^2 eps^-1, log2 delta^-1).  The asymptotic proof
 * leaves the sufficiently-large constant and finite rounding open; this
 * artifact fixes base-2 logs and ceil rounding and reports C in every row. */
static u64 eh_f(Table *T, double eps) {
    double l = log2(1.0 / eps);
    double a = l * l;
    double b = log2(1.0 / T->delta);
    double v = EH_C * (a < b ? a : b);
    if (v < 1) v = 1;
    return (u64)ceil(v);
}
/* Place in the first free local probe of array a0.  A spill is a defensive
 * correctness fallback and is reported; paper-faithful runs are required to
 * finish with zero spills. */
static inline int eh_place_uncapped(Table *T, u64 key, int a0, OpCost *c) {
    for (int t = 0; t < T->k; t++) {
        int a = (a0 + t) % T->k;
        if (T->occ[a] >= T->sz[a]) continue;
        if (t) eh_spills++;
        for (u64 j = 0;; j++) {
            u64 s = eh_probe(T, key, a, j);
            touch(s);
            if (T->t[s] == EMPTY) {
                T->t[s] = key; T->occ[a]++;
                if (j + 1 > T->cap_arr[a]) T->cap_arr[a] = j + 1;
                op_end(c); return 1;
            }
        }
    }
    op_end(c); return 0;
}
/* Equation (1) in FKK Section II, with i 1-based. */
static u64 eh_batch_len(Table *T, int i) {
    u64 a = T->sz[i-1], b = T->sz[i];
    u64 fin = a - (u64)floor(T->delta * (double)a / 2.0);
    u64 start = (u64)ceil(0.75 * (double)a);
    u64 add = fin > start ? fin - start : 0;
    return add + (u64)ceil(0.75 * (double)b);
}

static int eh_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    int bi = T->batch;
    if (bi == 0) {
        if (T->occ[0] < (u64)ceil(0.75 * (double)T->sz[0]))
            return eh_place_uncapped(T, key, 0, c);
        T->batch = 1; bi = 1; eh_bleft = eh_batch_len(T, 1);
    }
    while (bi >= 1 && eh_bleft == 0) {
        T->batch = ++bi;
        if (bi > T->k - 1) break;
        eh_bleft = eh_batch_len(T, bi);
    }
    /* FKK's insertion sequence should end before this branch for the tested
     * parameter range. Keep a counted fallback to prevent undefined behavior. */
    if (bi > T->k - 1) return eh_place_uncapped(T, key, T->k - 1, c);

    int A = bi - 1, B = bi;
    if (eh_bleft) eh_bleft--;
    double e1 = 1.0 - (double)T->occ[A] / (double)T->sz[A];
    double e2 = 1.0 - (double)T->occ[B] / (double)T->sz[B];

    if (e1 > T->delta / 2.0 && e2 > 0.25) {
        u64 f = eh_f(T, e1 > 0 ? e1 : 1e-15);
        for (u64 j = 0; j < f; j++) {
            u64 s = eh_probe(T, key, A, j);
            touch(s);
            if (T->t[s] == EMPTY) {
                T->t[s] = key; T->occ[A]++;
                if (j + 1 > T->cap_arr[A]) T->cap_arr[A] = j + 1;
                op_end(c); return 1;
            }
        }
        /* Do NOT increase the query certificate after an exhausted prefix:
         * cap_arr records resident depth, not insertion work. */
        return eh_place_uncapped(T, key, B, c);
    } else if (e1 <= T->delta / 2.0) {
        return eh_place_uncapped(T, key, B, c);
    } else {
        return eh_place_uncapped(T, key, A, c);
    }
}

/* Exact FKK injection from Lemma 1. i1,j1 are positive (1-based).  Its binary
 * representation is 1 b1 1 b2 ... 1 bq 0 a1...ap, where b encodes j and a
 * encodes i. */
static u64 eh_phi(u64 i1, u64 j1) {
    int p = 64 - __builtin_clzll(i1);
    int q = 64 - __builtin_clzll(j1);
    u64 v = 0;
    for (int bit = q - 1; bit >= 0; bit--) {
        v = (v << 1) | 1ULL;
        v = (v << 1) | ((j1 >> bit) & 1ULL);
    }
    v <<= 1; /* separator 0 */
    for (int bit = p - 1; bit >= 0; bit--) v = (v << 1) | ((i1 >> bit) & 1ULL);
    return v;
}

typedef struct { int i; u64 j; u64 phi; } EHOrd;
static EHOrd *ord = NULL;
static size_t ord_n = 0;
static int ord_cmp(const void *aa, const void *bb) {
    const EHOrd *a = (const EHOrd *)aa, *b = (const EHOrd *)bb;
    if (a->phi < b->phi) return -1;
    if (a->phi > b->phi) return 1;
    if (a->i != b->i) return a->i - b->i;
    return (a->j > b->j) - (a->j < b->j);
}
static void eh_build_order(Table *T) {
    free(ord); ord = NULL; ord_n = 0;
    size_t need = 0;
    for (int i = 0; i < T->k; i++) need += (size_t)T->cap_arr[i];
    ord = need ? malloc(need * sizeof(*ord)) : NULL;
    if (need && !ord) { fprintf(stderr, "out of memory building elastic query order\n"); exit(3); }
    for (int i = 0; i < T->k; i++) {
        for (u64 j = 1; j <= T->cap_arr[i]; j++) {
            ord[ord_n].i = i; ord[ord_n].j = j - 1;
            ord[ord_n].phi = eh_phi((u64)i + 1, j);
            ord_n++;
        }
    }
    qsort(ord, ord_n, sizeof(*ord), ord_cmp);
}
static int eh_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    unsigned char done[64]; memset(done, 0, sizeof done);
    int left = 0;
    for (int i = 0; i < T->k; i++) if (T->cap_arr[i]) left++;
    for (size_t t = 0; t < ord_n && left > 0; t++) {
        int i = ord[t].i; u64 j = ord[t].j;
        if (done[i]) continue;
        u64 s = eh_probe(T, key, i, j);
        touch(s);
        if (T->t[s] == key) { op_end(c); return 1; }
        /* An empty earlier in this key's local probe sequence is a valid
         * insertion-only certificate; otherwise D_i=cap_arr[i] is sufficient. */
        if (T->t[s] == EMPTY || j + 1 == T->cap_arr[i]) { done[i] = 1; left--; }
    }
    op_end(c); return 0;
}

/* ===== blocked elastic variant: one logical candidate is a contiguous bucket
 * of at most b slots.  The last bucket of a level may be shorter, so no access
 * crosses a level boundary. */
static inline void ehb_bucket(Table *T, u64 key, int i, u64 j, u64 b,
                              u64 *base, u64 *len) {
    u64 nb = (T->sz[i] + b - 1) / b;
    u64 bi = red(hash2(key, ((u64)(i + 1) << 32) ^ (j + 1)), nb);
    u64 local = bi * b;
    *base = T->off[i] + local;
    *len = T->sz[i] - local;
    if (*len > b) *len = b;
}
static inline int ehb_place_uncapped(Table *T, u64 key, int a0, u64 b, OpCost *c) {
    for (int t = 0; t < T->k; t++) {
        int a = (a0 + t) % T->k;
        if (T->occ[a] >= T->sz[a]) continue;
        if (t) eh_spills++;
        for (u64 j = 0;; j++) {
            u64 base, len; ehb_bucket(T, key, a, j, b, &base, &len);
            for (u64 q = 0; q < len; q++) {
                if (T->t[base + q] == EMPTY) {
                    touch_run(base, q + 1); T->t[base + q] = key; T->occ[a]++;
                    if (j + 1 > T->cap_arr[a]) T->cap_arr[a] = j + 1;
                    op_end(c); return 1;
                }
            }
            touch_run(base, len);
            if (j > 200ULL * T->sz[a] / (T->sz[a] - T->occ[a] + 1) + 1000) break;
        }
    }
    op_end(c); return 0;
}
static int ehb_insert(Table *T, u64 key, OpCost *c) {
    u64 b = (u64)T->bsize;
    op_begin();
    int bi = T->batch;
    if (bi == 0) {
        if (T->occ[0] < (u64)ceil(0.75 * (double)T->sz[0]))
            return ehb_place_uncapped(T, key, 0, b, c);
        T->batch = 1; bi = 1; eh_bleft = eh_batch_len(T, 1);
    }
    while (bi >= 1 && eh_bleft == 0) {
        T->batch = ++bi;
        if (bi > T->k - 1) break;
        eh_bleft = eh_batch_len(T, bi);
    }
    if (bi > T->k - 1) return ehb_place_uncapped(T, key, T->k - 1, b, c);
    int A = bi - 1, B = bi;
    if (eh_bleft) eh_bleft--;
    double e1 = 1.0 - (double)T->occ[A] / (double)T->sz[A];
    double e2 = 1.0 - (double)T->occ[B] / (double)T->sz[B];
    if (e1 > T->delta / 2.0 && e2 > 0.25) {
        u64 f = eh_f(T, e1 > 0 ? e1 : 1e-15);
        for (u64 j = 0; j < f; j++) {
            u64 base, len; ehb_bucket(T, key, A, j, b, &base, &len);
            for (u64 q = 0; q < len; q++) {
                if (T->t[base + q] == EMPTY) {
                    touch_run(base, q + 1); T->t[base + q] = key; T->occ[A]++;
                    if (j + 1 > T->cap_arr[A]) T->cap_arr[A] = j + 1;
                    op_end(c); return 1;
                }
            }
            touch_run(base, len);
        }
        return ehb_place_uncapped(T, key, B, b, c);
    } else if (e1 <= T->delta / 2.0) return ehb_place_uncapped(T, key, B, b, c);
    else return ehb_place_uncapped(T, key, A, b, c);
}
static int ehb_find(Table *T, u64 key, OpCost *c) {
    u64 b = (u64)T->bsize;
    op_begin();
    unsigned char done[64]; memset(done, 0, sizeof done);
    int left = 0;
    for (int i = 0; i < T->k; i++) if (T->cap_arr[i]) left++;
    for (size_t t = 0; t < ord_n && left > 0; t++) {
        int i = ord[t].i; u64 j = ord[t].j;
        if (done[i]) continue;
        u64 base, len; ehb_bucket(T, key, i, j, b, &base, &len);
        int sawempty = 0;
        for (u64 q = 0; q < len; q++) {
            if (T->t[base + q] == key) { touch_run(base, q + 1); op_end(c); return 1; }
            if (T->t[base + q] == EMPTY) { sawempty = 1; touch_run(base, q + 1); break; }
        }
        if (!sawempty) touch_run(base, len);
        if (sawempty || j + 1 == T->cap_arr[i]) { done[i] = 1; left--; }
    }
    op_end(c); return 0;
}

/* ============ funnel hashing (FKK) ============ */
static int FH_ALIGN = 0;   /* round beta up to a whole 64B line and align levels */

static void fh_setup(Table *T, double delta) {
    u64 n = T->n;
    int alpha = (int)ceil(4.0 * log2(1.0 / delta) + 10.0);
    int beta  = (int)ceil(2.0 * log2(1.0 / delta));
    if (beta < 2) beta = 2;
    if (FH_ALIGN) beta = ((beta + SLOTS_PER_LINE - 1) / SLOTS_PER_LINE) * SLOTS_PER_LINE;
    if (alpha > 60) alpha = 60;
    T->alpha = alpha; T->beta = beta; T->delta = delta;

    u64 spec = (u64)(0.6 * delta * (double)n);     /* in [delta n/2, 3 delta n/4] */
    if (spec < 16) spec = 16;
    u64 main_sz = n - spec;

    /* |A_i| = beta * a_i with a_{i+1} = 3 a_i / 4 ; sum beta*a_i = main_sz */
    double gs = 0, r = 1.0;
    for (int i = 0; i < alpha; i++) { gs += r; r *= 0.75; }
    double a1 = (double)main_sz / (beta * gs);
    u64 used = 0; r = 1.0;
    for (int i = 0; i < alpha; i++) {
        u64 ai = (u64)floor(a1 * r); if (ai < 1) ai = 1;
        u64 z = (u64)ai * (u64)beta;
        if (used + z > main_sz) z = (main_sz > used) ? (main_sz - used) : 0;
        z = (z / beta) * beta;
        T->lsz[i] = z; T->loff[i] = used; T->lbuckets[i] = z / beta;
        used += z; r *= 0.75;
    }
    T->spec_off = used; T->spec_sz = n - used;
    double lg = log2((double)n); if (lg < 2) lg = 2;
    T->llg = (int)ceil(log2(lg));
    T->cbucket = (int)ceil(2.0 * log2(lg));
    if (T->cbucket < 2) T->cbucket = 2;

    /* Finite divisibility convention for A_{alpha+1}.  FKK split it into
     * equal (+/-1) B and C halves and make C out of 2 log log n buckets.
     * Those requirements need not be simultaneously integral.  We round C
     * down to a whole number of buckets and assign the O(log log n) remainder
     * to B, so every table slot remains usable. */
    u64 half = T->spec_sz / 2;
    T->Csz = (half / (u64)T->cbucket) * (u64)T->cbucket;
    if (!T->Csz) T->Csz = (u64)T->cbucket;
    if (T->Csz >= T->spec_sz) T->Csz = (T->spec_sz / (u64)T->cbucket) * (u64)T->cbucket;
    T->Bsz = T->spec_sz - T->Csz;
    T->Boff = T->spec_off;
    T->Coff = T->spec_off + T->Bsz;

    if (delta <= 0.125) {
        u64 lo = (u64)ceil(delta * (double)n / 2.0);
        u64 hi = (u64)floor(3.0 * delta * (double)n / 4.0);
        if (T->spec_sz < lo || T->spec_sz > hi || T->Csz % (u64)T->cbucket) {
            fprintf(stderr, "funnel finite geometry outside documented constraints\n");
            exit(3);
        }
        for (int i = 0; i + 1 < alpha; i++) {
            if (!T->lsz[i] || !T->lsz[i+1]) continue;
            double ai = (double)(T->lsz[i] / (u64)beta);
            double aj = (double)(T->lsz[i+1] / (u64)beta);
            if (fabs(aj - 0.75 * ai) > 1.0000001) {
                fprintf(stderr, "funnel main-array geometry violates 3/4 +/- 1\n");
                exit(3);
            }
        }
    }
}
static int fh_insert(Table *T, u64 key, OpCost *c) {
    op_begin();
    for (int i = 0; i < T->alpha; i++) {
        if (!T->lbuckets[i]) continue;
        u64 b = red(hash2(key, ((u64)(i + 1) << 40) ^ 0xABCD), T->lbuckets[i]);
        u64 base = T->loff[i] + b * (u64)T->beta;
        for (int q = 0; q < T->beta; q++) {
            if (T->t[base + q] == EMPTY) { touch_run(base, q + 1); T->t[base + q] = key; op_end(c); return 1; }
            if (T->t[base + q] == key)   { touch_run(base, q + 1); op_end(c); return 1; }
        }
        touch_run(base, T->beta);
    }
    /* special array: B = uniform probing, give up after log log n attempts */
    for (int j = 0; j < T->llg; j++) {
        u64 s = T->Boff + red(hash2(key, 0xB00000ULL + j), T->Bsz);
        touch(s);
        if (T->t[s] == EMPTY) { T->t[s] = key; op_end(c); return 1; }
        if (T->t[s] == key)   { op_end(c); return 1; }
    }
    /* C: exact FKK interleaving a_1,b_1,a_2,b_2,... .  With prefix-filled
     * buckets, the first empty lies in the emptier bucket (ties favor a). */
    {
        u64 nb = T->Csz / (u64)T->cbucket; if (!nb) nb = 1;
        u64 b1 = red(hash2(key, 0xC1), nb), b2 = red(hash2(key, 0xC2), nb);
        u64 base1 = T->Coff + b1 * T->cbucket, base2 = T->Coff + b2 * T->cbucket;
        for (int q = 0; q < T->cbucket; q++) {
            u64 s1 = base1 + (u64)q; touch(s1);
            if (T->t[s1] == EMPTY) { T->t[s1] = key; op_end(c); return 1; }
            if (T->t[s1] == key)   { op_end(c); return 1; }
            u64 s2 = base2 + (u64)q; touch(s2);
            if (T->t[s2] == EMPTY) { T->t[s2] = key; op_end(c); return 1; }
            if (T->t[s2] == key)   { op_end(c); return 1; }
        }
    }
    op_end(c); return 0;
}
static int fh_find(Table *T, u64 key, OpCost *c) {
    op_begin();
    for (int i = 0; i < T->alpha; i++) {
        if (!T->lbuckets[i]) continue;
        u64 b = red(hash2(key, ((u64)(i + 1) << 40) ^ 0xABCD), T->lbuckets[i]);
        u64 base = T->loff[i] + b * (u64)T->beta;
        int sawempty = 0;
        for (int q = 0; q < T->beta; q++) {
            if (T->t[base + q] == key)   { touch_run(base, q + 1); op_end(c); return 1; }
            if (T->t[base + q] == EMPTY) { sawempty = 1; touch_run(base, q + 1); break; }
        }
        if (sawempty) { op_end(c); return 0; }   /* bucket not full => key absent */
        touch_run(base, T->beta);
    }
    for (int j = 0; j < T->llg; j++) {
        u64 s = T->Boff + red(hash2(key, 0xB00000ULL + j), T->Bsz);
        touch(s);
        if (T->t[s] == key)   { op_end(c); return 1; }
        if (T->t[s] == EMPTY) { op_end(c); return 0; }
    }
    {
        u64 nb = T->Csz / (u64)T->cbucket; if (!nb) nb = 1;
        u64 b1 = red(hash2(key, 0xC1), nb), b2 = red(hash2(key, 0xC2), nb);
        u64 base1 = T->Coff + b1 * T->cbucket, base2 = T->Coff + b2 * T->cbucket;
        for (int q = 0; q < T->cbucket; q++) {
            u64 s1 = base1 + (u64)q; touch(s1);
            if (T->t[s1] == key)   { op_end(c); return 1; }
            if (T->t[s1] == EMPTY) { op_end(c); return 0; }
            u64 s2 = base2 + (u64)q; touch(s2);
            if (T->t[s2] == key)   { op_end(c); return 1; }
            if (T->t[s2] == EMPTY) { op_end(c); return 0; }
        }
    }
    op_end(c); return 0;
}


/* Deterministic payload-line upper bound for the finite Funnel probe schedule.
 * Any contiguous run of r slots can intersect at most ceil((r+b-1)/b)
 * b-slot payload lines, irrespective of alignment.  Funnel probes at most one
 * beta-slot bucket in each main array, llg singleton positions in B, and two
 * cbucket-slot buckets in C. */
static u64 interval_line_bound(u64 r) {
    if (!r) return 0;
    return (r + 2 * SLOTS_PER_LINE - 2) / SLOTS_PER_LINE;
}
static u64 fh_payload_line_bound(const Table *T) {
    u64 main_levels = 0;
    for (int i = 0; i < T->alpha; i++) if (T->lbuckets[i]) main_levels++;
    return main_levels * interval_line_bound((u64)T->beta)
         + (u64)T->llg
         + 2ULL * interval_line_bound((u64)T->cbucket);
}

/* State-aware unsuccessful-query bound. If a main level has no full bucket,
 * every miss terminates in that level.  Sum the worst line span of one bucket
 * in each reachable main level and stop at the first level with zero full
 * buckets. If all main levels contain a full bucket, use the general bound. */
static u64 fh_realized_miss_line_bound(const Table *T) {
    u64 total = 0;
    for (int i = 0; i < T->alpha; i++) {
        if (!T->lbuckets[i]) continue;
        u64 maxspan = 0, fullb = 0;
        for (u64 bb = 0; bb < T->lbuckets[i]; bb++) {
            u64 base = T->loff[i] + bb * (u64)T->beta;
            u64 start_line = base / SLOTS_PER_LINE;
            u64 end_line = (base + (u64)T->beta - 1) / SLOTS_PER_LINE;
            u64 span = end_line - start_line + 1;
            if (span > maxspan) maxspan = span;
            int full = 1;
            for (int q = 0; q < T->beta; q++) {
                if (T->t[base + (u64)q] == EMPTY) { full = 0; break; }
            }
            if (full) fullb++;
        }
        total += maxspan;
        if (fullb == 0) return total;
    }
    return fh_payload_line_bound(T);
}

/* ---------------- driver ---------------- */
typedef int (*InsFn)(Table *, u64, OpCost *);
typedef int (*FndFn)(Table *, u64, OpCost *);

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s scheme log2n delta [bsize] [ehC]\n", argv[0]); return 1; }
    const char *scheme = argv[1];
    int lg = atoi(argv[2]);
    double delta = atof(argv[3]);
    int bsize = (argc > 4) ? atoi(argv[4]) : 8;
    if (argc > 5) EH_C = atof(argv[5]);
    { const char *e = getenv("SEED"); if (e && *e) RUN_SEED = strtoull(e, 0, 10); }

    u64 n = 1ULL << lg;
    Table T; memset(&T, 0, sizeof T);
    tbl_alloc(&T, n);
    T.delta = delta; T.bsize = bsize;

    InsFn ins; FndFn fnd;
    if      (!strcmp(scheme, "up")) { ins = up_insert; fnd = up_find; }
    else if (!strcmp(scheme, "dh")) { ins = dh_insert; fnd = dh_find; }
    else if (!strcmp(scheme, "lp")) { ins = lp_insert; fnd = lp_find; }
    else if (!strcmp(scheme, "wtl")){ ins = wtl_insert; fnd = wtl_find; }
    else if (!strcmp(scheme, "bp")) { ins = bp_insert; fnd = bp_find; }
    else if (!strcmp(scheme, "eh")) { eh_setup(&T, delta); ins = eh_insert; fnd = eh_find; }
    else if (!strcmp(scheme, "ehb")){ eh_setup(&T, delta); ins = ehb_insert; fnd = ehb_find; }
    else if (!strcmp(scheme, "fh")) { fh_setup(&T, delta); ins = fh_insert; fnd = fh_find; }
    else if (!strcmp(scheme, "fha")){ FH_ALIGN = 1; fh_setup(&T, delta); ins = fh_insert; fnd = fh_find; }
    else { fprintf(stderr, "unknown scheme\n"); free(T.t); return 1; }

    u64 m = n - (u64)floor(delta * (double)n);     /* items to insert */

    u64 *keys = malloc(m * sizeof(u64));
    if (!keys) { fprintf(stderr, "out of memory allocating key array\n"); free(T.t); return 3; }
    for (u64 i = 0; i < m; i++) { u64 k = splitmix64(i * 0x5DEECE66DULL + 0x9E3779B9ULL + RUN_SEED * 0x94D049BB133111EBULL); keys[i] = k ? k : 1; }

    Acc ai, apos, aneg; acc_init(&ai, 1 << 16); acc_init(&apos, 1 << 16); acc_init(&aneg, 1 << 16);

#ifdef INSTRUMENT
    int query_only_instrument = 0;
    { const char *e = getenv("QUERY_ONLY_INSTRUMENT"); if (e && *e && strcmp(e, "0")) query_only_instrument = 1; }
    if (query_only_instrument) INSTRUMENT_ENABLED = 0;
#endif

    double t0 = now();
    u64 fails = 0;
    for (u64 i = 0; i < m; i++) {
        OpCost c; if (!ins(&T, keys[i], &c)) fails++;
        acc_add(&ai, c);
    }
    double t_build = now() - t0;

    /* Post-build empty-slot layout diagnostics.  These make concentration
     * visible for line-grouped/self-consistent schemes rather than assuming
     * the uniform-empty-set model used by Theorem 4. */
    u64 layout_lines = (T.n + SLOTS_PER_LINE - 1) / SLOTS_PER_LINE;
    u64 nonfull_lines = 0, empty_lines = 0, empty_slots = 0, max_empty_per_line = 0;
    long double empty_sq = 0.0L;
    for (u64 L = 0; L < layout_lines; L++) {
        u64 base = L * SLOTS_PER_LINE;
        u64 lim = base + SLOTS_PER_LINE; if (lim > T.n) lim = T.n;
        u64 r = 0;
        for (u64 s = base; s < lim; s++) if (T.t[s] == EMPTY) r++;
        empty_slots += r;
        if (r) nonfull_lines++;
        if (r == lim - base) empty_lines++;
        if (r > max_empty_per_line) max_empty_per_line = r;
        empty_sq += (long double)r * (long double)r;
    }
    printf("LAYOUTCSV,%s,%d,%.6f,%llu,%llu,%llu,%llu,%.9Lf\n",
           scheme, lg, delta,
           (unsigned long long)empty_slots, (unsigned long long)nonfull_lines,
           (unsigned long long)empty_lines, (unsigned long long)max_empty_per_line, empty_sq);


    if (!strcmp(scheme, "fh")) {
        const char *fp = getenv("FH_PROFILE");
        if (fp && *fp && strcmp(fp, "0")) {
            for (int i = 0; i < T.alpha; i++) {
                if (!T.lbuckets[i]) continue;
                u64 fullb = 0, nonemptyb = 0, occslots = 0;
                for (u64 bb = 0; bb < T.lbuckets[i]; bb++) {
                    u64 base = T.loff[i] + bb * (u64)T.beta, occb = 0;
                    for (int q = 0; q < T.beta; q++) {
                        if (T.t[base + (u64)q] != EMPTY) occb++;
                    }
                    occslots += occb;
                    if (occb) nonemptyb++;
                    if (occb == (u64)T.beta) fullb++;
                }
                printf("FHPROFILE,%d,%llu,%llu,%llu,%llu\n", i,
                       (unsigned long long)T.lsz[i],
                       (unsigned long long)occslots,
                       (unsigned long long)nonemptyb,
                       (unsigned long long)fullb);
            }
        }
    }

    if (!strcmp(scheme, "eh") || !strcmp(scheme, "ehb")) eh_build_order(&T);
#ifdef INSTRUMENT
    if (query_only_instrument) INSTRUMENT_ENABLED = 1;
#endif

    /* correctness: all keys must be found */
    u64 nq = m < 200000 ? m : 200000;
    { const char *e = getenv("NQ"); if (e) { u64 v = strtoull(e, 0, 10); if (v && v < nq) nq = v; } }
    u64 missing = 0;
    t0 = now();
    for (u64 q = 0; q < nq; q++) {
        u64 idx = splitmix64(q * 0x2545F4914F6CDD1DULL + RUN_SEED) % m;
        OpCost c; if (!fnd(&T, keys[idx], &c)) missing++;
        acc_add(&apos, c);
    }
    double t_pos = now() - t0;

    t0 = now();
    u64 falsepos = 0;
    for (u64 q = 0; q < nq; q++) {
        u64 k = splitmix64(0xDEADBEEF00000000ULL + q + RUN_SEED * 0xBF58476D1CE4E5B9ULL); if (!k) k = 3;
        OpCost c; if (fnd(&T, k, &c)) falsepos++;
        acc_add(&aneg, c);
    }
    double t_neg = now() - t0;

    printf("scheme=%s n=%llu delta=%.5f m=%llu seed=%llu fails=%llu missing=%llu falsepos=%llu spills=%llu\n",
           scheme, (unsigned long long)n, delta, (unsigned long long)m, (unsigned long long)RUN_SEED,
           (unsigned long long)fails, (unsigned long long)missing, (unsigned long long)falsepos, (unsigned long long)eh_spills);
    if (!strcmp(scheme, "up")) {
        u64 e = n - m;
        printf("UPTHEORY,n=%llu,empty=%llu,b=%d,probe_expect=%.9f,line_expect=%.9f\n",
               (unsigned long long)n, (unsigned long long)e, SLOTS_PER_LINE,
               (double)n / (double)e, up_expected_miss_lines(n, e, SLOTS_PER_LINE));
    }
    u64 fh_bound = !strcmp(scheme, "fh") ? fh_payload_line_bound(&T) : 0;
    u64 fh_state_bound = !strcmp(scheme, "fh") ? fh_realized_miss_line_bound(&T) : 0;
    if (!strcmp(scheme, "fh"))
        printf("FHBOUND,general=%llu,state_miss=%llu\n",
               (unsigned long long)fh_bound, (unsigned long long)fh_state_bound);
    printf("CSV,%s,%d,%.6f,%d,%.3f,"
           "%.4f,%.4f,%.4f,%llu,%llu,%llu,"
           "%.4f,%.4f,%llu,%llu,"
           "%.4f,%.4f,%llu,%llu,"
           "%.6f,%.6f,%.6f,%llu,%llu\n",
           scheme, lg, delta, bsize, EH_C,
           ai.sp / ai.n, ai.sl / ai.n, ai.spg / ai.n, (unsigned long long)acc_q(&ai, 0.99), (unsigned long long)ai.maxp, (unsigned long long)ai.maxl,
           apos.sp / apos.n, apos.sl / apos.n, (unsigned long long)acc_q(&apos, 0.99), (unsigned long long)apos.maxp,
           aneg.sp / aneg.n, aneg.sl / aneg.n, (unsigned long long)acc_q(&aneg, 0.99), (unsigned long long)aneg.maxp,
           t_build, t_pos, t_neg, (unsigned long long)m, (unsigned long long)nq);
    printf("CHECKCSV,%s,%d,%.6f,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
           scheme, lg, delta,
           (unsigned long long)fails, (unsigned long long)missing, (unsigned long long)falsepos,
           (unsigned long long)eh_spills, (unsigned long long)ai.maxl,
           (unsigned long long)apos.maxl, (unsigned long long)aneg.maxl,
           (unsigned long long)fh_bound, (unsigned long long)fh_state_bound);
    int status = (fails || missing || falsepos) ? 2 : 0;
    if (!strcmp(scheme, "fh") && (ai.maxl > fh_bound || apos.maxl > fh_bound || aneg.maxl > fh_bound)) {
        fprintf(stderr, "funnel payload-line structural bound violated\n");
        status = 2;
    }
    if (!strcmp(scheme, "fh") && aneg.maxl > fh_state_bound) {
        fprintf(stderr, "funnel state-aware miss bound violated\n");
        status = 2;
    }
    if ((!strcmp(scheme, "eh") || !strcmp(scheme, "ehb")) && eh_spills) {
        fprintf(stderr, "elastic spill occurred; finite realization left the intended FKK batch path\n");
        status = 2;
    }

    acc_free(&ai);
    acc_free(&apos);
    acc_free(&aneg);
    free(keys);
    free(ord); ord = NULL; ord_n = 0;
    free(T.t); T.t = NULL;
    return status;
}
