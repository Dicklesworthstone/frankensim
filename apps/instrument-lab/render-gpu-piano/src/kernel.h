// The contact-free grand_piano step (plan.rs) as ONE persistent kernel, shared by Metal and CUDA.
//
// A prelude (metal.rs / cuda.rs) defines `real` (double-float pairs on hardware without fast FP64, or
// double), R(x), the address-space qualifiers PC (parameters) CG CI (read-only data) SG (state, coherent
// across blocks) OG (outputs) TGM (threadgroup memory) TH (thread), shfl, warp_sum, block_sum,
// BLOCK_SYNC (threadgroup and device memory), BLOCK_SYNC_TG (threadgroup memory only), GRID_SYNC,
// UNROLL / UNROLL_ALL, NW (warps per block), the P_* parameter slots, and KERNEL_SIGNATURE /
// KERNEL_PROLOGUE (tid, blk and the dynamic threadgroup memory `tmem`).
//
// Every block must be resident at once (the host keeps the count within the GPU's budget). Group blocks
// own one or two adjacent bridge groups (GBLK): their string velocities live in threadgroup memory for the
// whole launch and all 32*NW threads share the per-mode work. Warp w owns the block's damper pad w (at
// most NW) with one station per lane: a damper call runs in station space (emulate.rs), the lane holding
// its station's speed u = s.v and accumulated change Delta in registers, and the string velocities take
// S Delta once per step. Board blocks give each warp RPW (1 or 2) board rows. Phases are separated by
// grid barriers, 7 per step with dampers engaged and 3 without (plan.rs):
//   D_b   pads: bridge speed of direction 0, station changes, chain of direction 1;  rows: board after dir 0
//   D_c   pads: bridge speed of direction 1, station changes (then S Delta, or the next chain);  rows: dir 1
//   Y1    groups: string loss, free motion, reaction;  rows: damped board velocity, free board motion
//   P2    groups: end state, finish ledger, contact positions;  rows: board end position
//   P3    rows: board velocity, flow, string loss, ledger;  groups: string loss (commit), station speeds, chain
// Each thread accumulates the modal and damper losses it owns; the block sums are written once.

#define NT (32 * NW)

KERNEL_SIGNATURE {
  KERNEL_PROLOGUE
  const int n = P[P_N], r = P[P_R], ng = P[P_NG], np = P[P_NP], nc = P[P_NC];
  const int steps = P[P_STEPS], damp = P[P_DAMP];
  const int lane = tid & 31, warp = tid >> 5;
  const int GB = P[P_GROUP_BLOCKS], BB = P[P_BOARD_BLOCKS], RPW = P[P_RPW];
  // P_SKIP (a timing diagnostic; results are wrong): bit 0 runs only the barriers; bits 1-4 skip the
  // superpositions, chains, board damper rows, station speeds and S Delta
  const int skip = P[P_SKIP];
  const int is_group = blk < GB && (skip & 1) == 0;
  const int is_board = blk >= GB && (skip & 1) == 0;

  CG real* QQ = C + P[P_QQ]; CG real* QV = C + P[P_QV]; CG real* VQ = C + P[P_VQ]; CG real* VV = C + P[P_VV];
  CG real* BQ = C + P[P_BQ]; CG real* BV = C + P[P_BV]; CG real* OM2 = C + P[P_OM2];
  CG real* AK = C + P[P_A]; CG real* BETA = C + P[P_BETA]; CG real* HAM = C + P[P_HAM];
  CG real* INV = C + P[P_INV]; CG real* CROSS = C + P[P_CROSS]; CG real* MASS = C + P[P_MASS]; CG real* STRAIN = C + P[P_STRAIN];
  const real dur = C[P[P_DUR]];
  CG real* GM = C + P[P_GM]; CG real* R1 = C + P[P_R1]; CG real* R2 = C + P[P_R2]; CG real* R3 = C + P[P_R3];
  CG real* R4 = C + P[P_R4]; CG real* R5 = C + P[P_R5]; CG real* R6 = C + P[P_R6];
  CG real* K1 = C + P[P_K1]; CG real* K2 = C + P[P_K2]; CG real* K3 = C + P[P_K3]; CG real* K4 = C + P[P_K4];
  CG real* K5 = C + P[P_K5]; CG real* K6 = C + P[P_K6]; CG real* K7 = C + P[P_K7];
  CG real* E2 = C + P[P_E2]; CG real* E3 = C + P[P_E3]; CG real* GT = C + P[P_GT]; CG real* CLIFT = C + P[P_CLIFT];
  CG real* ULOSS0 = D + P[P_ULOSS0]; CG real* ULOSS1 = D + P[P_ULOSS1]; CG real* SHAPE = D + P[P_SHAPE];
  CG real* UC0 = D + P[P_UC0]; CG real* UC1 = D + P[P_UC1]; CG real* GUC0 = D + P[P_GUC0]; CG real* GUC1 = D + P[P_GUC1];
  CG real* TM0 = D + P[P_TM0]; CG real* TM1 = D + P[P_TM1]; CG real* GTM0 = D + P[P_GTM0]; CG real* GTM1 = D + P[P_GTM1];
  CG real* LAM0 = D + P[P_LAM0]; CG real* LAM1 = D + P[P_LAM1];
  CG real* RHO0 = D + P[P_RHO0]; CG real* RHO1 = D + P[P_RHO1]; CG real* MU0 = D + P[P_MU0]; CG real* MU1 = D + P[P_MU1];
  CG real* PM0 = D + P[P_P0]; CG real* QM0 = D + P[P_Q0]; CG real* VM0 = D + P[P_V0]; CG real* AM0 = D + P[P_AM0];
  CG real* PM1 = D + P[P_P1]; CG real* QM1 = D + P[P_Q1]; CG real* VM1 = D + P[P_V1]; CG real* AM1 = D + P[P_AM1];
  CI int* GSTART = I + P[P_GSTART]; CI int* GPAD = I + P[P_GPAD]; CI int* GCON = I + P[P_GCON];
  CI int* CONLIST = I + P[P_CONLIST]; CI int* CSTART = I + P[P_CSTART]; CI int* CEND = I + P[P_CEND];
  CI int* PSTART = I + P[P_PSTART]; CI int* PEND = I + P[P_PEND]; CI int* PST = I + P[P_PST]; CI int* PSH = I + P[P_PSH];
  CI int* TOFF = I + P[P_TOFF]; CI int* STPAD = I + P[P_STPAD]; CI int* MPAD = I + P[P_MPAD]; CI int* GBLK = I + P[P_GBLK];
  SG real* Q = S + P[P_SQ]; SG real* V = S + P[P_SV]; SG real* FQ = S + P[P_SFQ]; SG real* FV = S + P[P_SFV];
  SG real* VBT = S + P[P_SVBT]; SG real* SRC = S + P[P_SSRC]; SG real* DD = S + P[P_SDD];
  SG real* REA = S + P[P_SREA]; SG real* REA2 = S + P[P_SREA2]; SG real* SRC2 = S + P[P_SSRC2]; SG real* DD2 = S + P[P_SDD2];
  SG real* A0 = S + P[P_SA0]; SG real* A1 = S + P[P_SA1];
  // threadgroup memory: the block's string velocities, its stations' changes, per-thread scratch, sums
  TGM real* tv = tmem;
  TGM real* tc = tv + P[P_MAXG];
  TGM real* scratch = tc + P[P_MAXS];
  TGM real* sh = scratch + NT;

  real modal = R(0.0), dloss = R(0.0);
  real tot[8];
  // group blocks: group ga (and gb = ga + 1 when `two`), modes [bs, be), the second group's from `split`
  int ga = 0, two = 0, bs = 0, be = 0, split = 0, pb = 0, pe = 0, cb = 0, ce = 0;
  if (is_group) {
    ga = GBLK[blk]; two = GBLK[blk + 1] - ga == 2;
    bs = GSTART[ga]; split = GSTART[ga + 1]; be = GSTART[ga + 1 + two];
    pb = GPAD[ga]; pe = GPAD[ga + 1 + two]; cb = GCON[ga]; ce = GCON[ga + 1 + two];
  }
  const int gb = ga + 1;
  real b0a = R(0.0), b0b = R(0.0);
  // board blocks: this warp's rows j0 (and j1 when RPW = 2), with registers kept between phases
  const int j0 = blk >= GB ? (blk - GB) * NW + warp : r;
  const int j1 = (blk >= GB && RPW == 2) ? j0 + BB * NW : r;
  real dvb0 = R(0.0), qbo0 = R(0.0), qn0 = R(0.0), dvb1 = R(0.0), qbo1 = R(0.0), qn1 = R(0.0);
  // this warp's pad (-1: none, or not engaged) and its stations; lane i < pns owns station i
  const int mp = (is_group && pb + warp < pe && ACT[pb + warp]) ? pb + warp : -1;
  const int pst = mp >= 0 ? PST[mp] : 0, pns = mp >= 0 ? PST[mp + 1] - pst : 0;
  const int mine = lane < pns;
  real su = R(0.0), sdel = R(0.0), sc = R(0.0), sgc = R(0.0), lreg = R(0.0);

// This warp's pad: the chain of direction DIR from its station speeds (plan::ChainOp): station changes
// c = T u and their effect G T u on the station speeds, the board amount moved, the loss cross term
// (kept for the superposition) and the loss u' Lambda u. One shuffle loop serves the three products.
#define CHAIN(DIR) if ((skip & 4) == 0) do {                                                          \
    CG real* t_ = ((DIR) == 0 ? TM0 : TM1) + TOFF[mp] + lane * pns;                                    \
    CG real* g_ = ((DIR) == 0 ? GTM0 : GTM1) + TOFF[mp] + lane * pns;                                  \
    CG real* l_ = ((DIR) == 0 ? LAM0 : LAM1) + TOFF[mp] + lane * pns;                                  \
    real c_ = R(0.0), gc_ = R(0.0), lu_ = R(0.0);                                                      \
    /* the pad's stations, in fours (pns <= 32, so every shuffled lane exists) */                     \
    for (int x0_ = 0; x0_ < pns; x0_ += 4) {                                                           \
      UNROLL_ALL                                                                                       \
      for (int x_ = x0_; x_ < x0_ + 4; x_++) {                                                         \
        const real ux_ = shfl(su, x_);                                                                 \
        if (mine && x_ < pns) { c_ = c_ + t_[x_] * ux_; gc_ = gc_ + g_[x_] * ux_; lu_ = lu_ + l_[x_] * ux_; } \
      }                                                                                                \
    }                                                                                                  \
    sc = c_; sgc = gc_;                                                                                \
    const real mv_ = warp_sum(mine ? ((DIR) == 0 ? MU0 : MU1)[pst + lane] * su : R(0.0));              \
    lreg = warp_sum(mine ? ((DIR) == 0 ? RHO0 : RHO1)[pst + lane] * su : R(0.0));                      \
    const real ls_ = warp_sum(mine ? su * lu_ : R(0.0));                                               \
    if (lane == 0) { if ((DIR) == 0) A0[mp] = mv_; else A1[mp] = mv_; dloss = dloss + ls_; }           \
  } while (0)

// This warp's pad: direction DIR's bridge speed from the board velocity XB and every pad's moved amount;
// its station changes c + beta uc (accumulated) move the station speeds by G T u + beta G uc; its loss.
#define SUPERPOSE(DIR, XB) if ((skip & 2) == 0) do {                                                  \
    CG real* PM_ = (DIR) == 0 ? PM0 : PM1; CG real* QM_ = (DIR) == 0 ? QM0 : QM1;                      \
    SG real* AA_ = (DIR) == 0 ? A0 : A1;                                                               \
    real part_ = R(0.0), padp_ = R(0.0);                                                               \
    UNROLL                                                                                             \
    for (int i_ = lane; i_ < r; i_ += 32) part_ = part_ + PM_[mp * r + i_] * (XB)[i_];                 \
    UNROLL                                                                                             \
    for (int i_ = lane; i_ < np; i_ += 32) padp_ = padp_ + QM_[mp * np + i_] * AA_[i_];                \
    const real beta_ = warp_sum(part_ + padp_);                                                        \
    if (mine) {                                                                                        \
      sdel = sdel + sc + beta_ * ((DIR) == 0 ? UC0 : UC1)[pst + lane];                                 \
      su = su + sgc + beta_ * ((DIR) == 0 ? GUC0 : GUC1)[pst + lane];                                  \
    }                                                                                                  \
    if (lane == 0) dloss = dloss + beta_ * lreg + beta_ * beta_ * ((DIR) == 0 ? ULOSS0[mp] : ULOSS1[mp]); \
  } while (0)

// The station speeds u = s.v of every engaged pad of the block, from the string velocities (each
// station's dot product split over NT / stations threads); the accumulated changes restart at zero.
#define STATIONS() if ((skip & 16) == 0) do {                                                         \
    const int st0_ = PST[pb], tot_ = PST[pe] - st0_;                                                   \
    const int kk_ = tot_ > 0 ? NT / tot_ : 1;                                                          \
    real part_ = R(0.0);                                                                               \
    if (tid < tot_ * kk_) {                                                                            \
      const int s_ = st0_ + tid / kk_, j_ = tid % kk_, p_ = STPAD[s_];                                 \
      if (ACT[p_]) {                                                                                   \
        const int ms_ = PSTART[p_], m_ = PEND[p_] - ms_;                                               \
        CG real* sv_ = SHAPE + PSH[p_] + (s_ - PST[p_]) * m_;                                          \
        real odd_ = R(0.0);                                                                            \
        int k_ = j_;                                                                                   \
        UNROLL                                                                                         \
        for (; k_ + kk_ < m_; k_ += 2 * kk_) {                                                         \
          part_ = part_ + sv_[k_] * tv[ms_ - bs + k_];                                                 \
          odd_ = odd_ + sv_[k_ + kk_] * tv[ms_ - bs + k_ + kk_];                                       \
        }                                                                                              \
        if (k_ < m_) part_ = part_ + sv_[k_] * tv[ms_ - bs + k_];                                      \
        part_ = part_ + odd_;                                                                          \
      }                                                                                                \
    }                                                                                                  \
    scratch[tid] = part_;                                                                              \
    BLOCK_SYNC_TG();                                                                                   \
    if (mine) {                                                                                        \
      real u_ = R(0.0);                                                                                \
      for (int j_ = 0; j_ < kk_; j_++) u_ = u_ + scratch[(pst - st0_ + lane) * kk_ + j_];              \
      su = u_;                                                                                         \
    }                                                                                                  \
    sdel = R(0.0);                                                                                     \
    BLOCK_SYNC_TG();                                                                                   \
  } while (0)

// The accumulated station changes into the block's string velocities: v += S Delta.
#define APPLY() if ((skip & 16) == 0) do {                                                            \
    const int st0_ = PST[pb];                                                                          \
    if (mine) tc[pst - st0_ + lane] = sdel;                                                            \
    BLOCK_SYNC_TG();                                                                                   \
    UNROLL                                                                                             \
    for (int k_ = bs + tid; k_ < be; k_ += NT) {                                                       \
      const int p_ = MPAD[k_];                                                                         \
      if (p_ >= 0 && ACT[p_]) {                                                                        \
        const int ms_ = PSTART[p_], m_ = PEND[p_] - ms_, base_ = PST[p_], ns_ = PST[p_ + 1] - base_;   \
        CG real* col_ = SHAPE + PSH[p_] + (k_ - ms_);                                                  \
        real acc_ = R(0.0), odd_ = R(0.0);                                                             \
        UNROLL                                                                                         \
        for (int i_ = 0; i_ < ns_; i_ += 2) {                                                          \
          acc_ = acc_ + col_[i_ * m_] * tc[base_ - st0_ + i_];                                         \
          if (i_ + 1 < ns_) odd_ = odd_ + col_[(i_ + 1) * m_] * tc[base_ - st0_ + i_ + 1];             \
        }                                                                                              \
        tv[k_ - bs] = tv[k_ - bs] + (acc_ + odd_);                                                     \
      }                                                                                                \
    }                                                                                                  \
    BLOCK_SYNC_TG();                                                                                   \
  } while (0)

// Board row J after direction DIR: V_DIR XB + A_DIR a_DIR, into DEST[J].
#define BOARD_DAMP(DIR, XB, DEST, J) if ((skip & 8) == 0) do {                                        \
    CG real* VM_ = (DIR) == 0 ? VM0 : VM1; CG real* AM_ = (DIR) == 0 ? AM0 : AM1;                      \
    SG real* AA_ = (DIR) == 0 ? A0 : A1;                                                               \
    real part_ = R(0.0), padp_ = R(0.0);                                                               \
    UNROLL                                                                                             \
    for (int i_ = lane; i_ < r; i_ += 32) part_ = part_ + VM_[(J) * r + i_] * (XB)[i_];                \
    UNROLL                                                                                             \
    for (int i_ = lane; i_ < np; i_ += 32) padp_ = padp_ + AM_[(J) * np + i_] * AA_[i_];               \
    const real out_ = warp_sum(part_ + padp_);                                                         \
    if (lane == 0) (DEST)[J] = out_;                                                                   \
  } while (0)

// The string-loss sources of the block's groups' present velocities.
#define SOURCES() do {                                                                                 \
    real sa_ = R(0.0), da_ = R(0.0), sb_ = R(0.0), db_ = R(0.0);                                       \
    UNROLL                                                                                             \
    for (int k_ = bs + tid; k_ < be; k_ += NT) {                                                       \
      const real x_ = tv[k_ - bs];                                                                     \
      const real s_ = CROSS[k_] * x_, d_ = STRAIN[k_] * BETA[k_] * INV[k_] * x_;                       \
      if (k_ < split) { sa_ = sa_ + s_; da_ = da_ + d_; } else { sb_ = sb_ + s_; db_ = db_ + d_; }      \
    }                                                                                                  \
    tot[0] = sa_; tot[1] = da_; tot[2] = sb_; tot[3] = db_;                                            \
    block_sum(tot, 2 + 2 * two, sh, tid);                                                              \
    if (tid == 0) { SRC[ga] = tot[0]; DD[ga] = tot[1]; if (two) { SRC[gb] = tot[2]; DD[gb] = tot[3]; } } \
  } while (0)

// The block's per-mode string-loss pass at bridge speeds BA (first group) and BB_ (second).
#define STRING_LOSS(BA, BB_) do {                                                                      \
    UNROLL                                                                                             \
    for (int k_ = bs + tid; k_ < be; k_ += NT) {                                                       \
      const real b_ = k_ < split ? (BA) : (BB_);                                                       \
      const real old_ = tv[k_ - bs];                                                                   \
      const real mid_ = INV[k_] * old_ + CROSS[k_] * b_;                                               \
      const real rel_ = mid_ - BETA[k_] * b_;                                                          \
      const real nv_ = old_ - dur * (MASS[k_] * mid_ + STRAIN[k_] * rel_);                             \
      modal = modal + R(0.5) * (old_ - nv_) * (old_ + nv_);                                            \
      tv[k_ - bs] = nv_;                                                                               \
    }                                                                                                  \
  } while (0)

// Add this thread's part of the block's groups' row products M[g] . X (row length LEN) to PA and PB.
#define GROUP_DOT(M, X, LEN, PA, PB) do {                                                              \
    UNROLL                                                                                             \
    for (int i_ = tid; i_ < (LEN); i_ += NT) {                                                         \
      const real x_ = (X)[i_];                                                                         \
      PA = PA + (M)[ga * (LEN) + i_] * x_;                                                             \
      if (two) PB = PB + (M)[gb * (LEN) + i_] * x_;                                                    \
    }                                                                                                  \
  } while (0)

  // prologue: the block's velocities into threadgroup memory; the first chains or string-loss sources
  if (is_group) {
    UNROLL
    for (int k = bs + tid; k < be; k += NT) tv[k - bs] = V[k];
    BLOCK_SYNC_TG();
    if (damp) { STATIONS(); if (mp >= 0) CHAIN(0); } else SOURCES();
  }
  GRID_SYNC();

  for (int step = 0; step < steps; step++) {
    if (damp) {
      // D_b (first damper call)
      if (is_group) {
        if (mp >= 0) { SUPERPOSE(0, V + n); CHAIN(1); }
      } else if (is_board) {
        if (j0 < r) BOARD_DAMP(0, V + n, VBT, j0);
        if (j1 < r) BOARD_DAMP(0, V + n, VBT, j1);
      }
      GRID_SYNC();
      // D_c, the step's station changes into the strings, the predicted velocity's string-loss sources
      if (is_group) {
        if (mp >= 0) SUPERPOSE(1, VBT);
        APPLY();
        SOURCES();
      } else if (is_board) {
        if (j0 < r) BOARD_DAMP(1, VBT, V + n, j0);
        if (j1 < r) BOARD_DAMP(1, VBT, V + n, j1);
      }
      GRID_SYNC();
    }
    // Y1: string loss (b = R1 v_b + R2 src), free motion, reactions (b0 = GM q_b)
    if (is_group) {
      real pa = R(0.0), pbb = R(0.0), qa = R(0.0), qb = R(0.0);
      GROUP_DOT(R1, V + n, r, pa, pbb);
      GROUP_DOT(R2, SRC, ng, pa, pbb);
      GROUP_DOT(GM, Q + n, r, qa, qb);
      tot[0] = pa; tot[1] = qa; tot[2] = pbb; tot[3] = qb;
      block_sum(tot, 2 + 2 * two, sh, tid);
      const real ba = tot[0], bbv = tot[2];
      b0a = tot[1]; b0b = tot[3];
      real ra = R(0.0), rb = R(0.0);
      UNROLL
      for (int k = bs + tid; k < be; k += NT) {
        const int second = k >= split;
        const real b = second ? bbv : ba;
        const real b0 = second ? b0b : b0a;
        const real old = tv[k - bs];
        const real mid = INV[k] * old + CROSS[k] * b;
        const real rel = mid - BETA[k] * b;
        const real nv = old - dur * (MASS[k] * mid + STRAIN[k] * rel);
        modal = modal + R(0.5) * (old - nv) * (old + nv);
        const real q = Q[k];
        const real fq = QQ[k] * q + QV[k] * nv;
        tv[k - bs] = nv;
        FQ[k] = fq;
        FV[k] = VQ[k] * q + VV[k] * nv;
        const real x = R(0.5) * AK[k] * (q + fq + R(0.5) * BQ[k] * AK[k] * b0);
        if (second) rb = rb + x; else ra = ra + x;
      }
      tot[0] = ra; tot[1] = rb;
      block_sum(tot, 1 + two, sh, tid);
      if (tid == 0) { REA[ga] = tot[0]; if (two) REA[gb] = tot[1]; }
    } else if (is_board) {
      UNROLL_ALL
      for (int t = 0; t < 2; t++) {
        const int j = t == 0 ? j0 : j1;
        if (j < r) {
          real part = R(0.0);
          UNROLL
          for (int i = lane; i < r; i += 32) part = part + K1[j * r + i] * V[n + i];
          UNROLL
          for (int i = lane; i < ng; i += 32) part = part + K2[j * ng + i] * SRC[i] + K3[j * ng + i] * DD[i];
          const real dvb = warp_sum(part);
          if (t == 0) dvb0 = dvb; else dvb1 = dvb;
          if (lane == 0) {
            const int k = n + j;
            const real vb = V[k];
            modal = modal + R(0.5) * (vb * vb - dvb * dvb);
            const real q = Q[k];
            FQ[k] = QQ[k] * q + QV[k] * dvb;
            FV[k] = VQ[k] * q + VV[k] * dvb;
          }
        }
      }
    }
    GRID_SYNC();
    // P2: end state with the finish ledger (b1 = R3 fq_b + R4 rea), contact positions
    if (is_group) {
      real pa = R(0.0), pbb = R(0.0);
      GROUP_DOT(R3, FQ + n, r, pa, pbb);
      GROUP_DOT(R4, REA, ng, pa, pbb);
      tot[0] = pa; tot[1] = pbb;
      block_sum(tot, 1 + two, sh, tid);
      const real b1a = tot[0], b1b = tot[1];
      const real bbara = R(0.5) * (b0a + b1a), bbarb = R(0.5) * (b0b + b1b);
      real ra = R(0.0), sa = R(0.0), da = R(0.0), rb = R(0.0), sb = R(0.0), db = R(0.0);
      UNROLL
      for (int k = bs + tid; k < be; k += NT) {
        const int second = k >= split;
        const real f = AK[k] * (second ? bbarb : bbara);
        const real nq = FQ[k] + BQ[k] * f;
        const real nv = FV[k] + BV[k] * f;
        const real q = Q[k];
        const real v = tv[k - bs];
        const real rx = AK[k] * R(0.5) * (q + nq);
        modal = modal + f * (nq - q) - R(0.5) * (nv * nv - v * v + OM2[k] * (nq * nq - q * q));
        Q[k] = nq;
        tv[k - bs] = nv;
        const real sx = CROSS[k] * nv, dx = STRAIN[k] * BETA[k] * INV[k] * nv;
        if (second) { rb = rb + rx; sb = sb + sx; db = db + dx; } else { ra = ra + rx; sa = sa + sx; da = da + dx; }
      }
      BLOCK_SYNC();
      tot[0] = ra; tot[1] = sa; tot[2] = da; tot[3] = rb; tot[4] = sb; tot[5] = db;
      block_sum(tot, 3 + 3 * two, sh, tid);
      if (tid == 0) {
        REA2[ga] = tot[0]; SRC2[ga] = tot[1]; DD2[ga] = tot[2];
        if (two) { REA2[gb] = tot[3]; SRC2[gb] = tot[4]; DD2[gb] = tot[5]; }
      }
      for (int ci = cb + warp; ci < ce; ci += NW) {
        const int c = CONLIST[ci];
        real part = R(0.0);
        UNROLL
        for (int k = CSTART[c] + lane; k < CEND[c]; k += 32) part = part + HAM[k] * Q[k];
        const real pos = warp_sum(part) + CLIFT[c] * (CSTART[c] >= split ? b1b : b1a);
        if (lane == 0) store_out(POS, step * nc + c, pos);
      }
    } else if (is_board) {
      UNROLL_ALL
      for (int t = 0; t < 2; t++) {
        const int j = t == 0 ? j0 : j1;
        if (j < r) {
          real part = R(0.0);
          UNROLL
          for (int i = lane; i < r; i += 32) part = part + K4[j * r + i] * FQ[n + i];
          UNROLL
          for (int i = lane; i < ng; i += 32) part = part + K5[j * ng + i] * REA[i];
          const real qn = warp_sum(part);
          const real qbo = Q[n + j];
          if (t == 0) { qn0 = qn; qbo0 = qbo; } else { qn1 = qn; qbo1 = qbo; }
          if (lane == 0) Q[n + j] = qn;
        }
      }
    }
    GRID_SYNC();
    // P3: board velocity, flow, string loss and ledger; the groups' string loss (b2 = R5 fv_b + R6 rea2
    // + R2 src2) and commit; then the station speeds and first chain of the second damper call
    if (is_group) {
      real pa = R(0.0), pbb = R(0.0);
      GROUP_DOT(R5, FV + n, r, pa, pbb);
      GROUP_DOT(R6, REA2, ng, pa, pbb);
      GROUP_DOT(R2, SRC2, ng, pa, pbb);
      tot[0] = pa; tot[1] = pbb;
      block_sum(tot, 1 + two, sh, tid);
      const real b2a = tot[0], b2b = tot[1];
      STRING_LOSS(b2a, b2b);
      BLOCK_SYNC_TG();
      if (damp) { STATIONS(); if (mp >= 0) CHAIN(0); } else SOURCES();
    } else if (is_board) {
      UNROLL_ALL
      for (int t = 0; t < 2; t++) {
        const int j = t == 0 ? j0 : j1;
        if (j < r) {
          real fp = R(0.0), part = R(0.0);
          UNROLL
          for (int i = lane; i < ng; i += 32) {
            const real ra = REA2[i];
            fp = fp + GT[j * ng + i] * ra;
            part = part + K7[j * ng + i] * ra + E2[j * ng + i] * SRC2[i] + E3[j * ng + i] * DD2[i];
          }
          UNROLL
          for (int i = lane; i < r; i += 32) part = part + K6[j * r + i] * FV[n + i];
          const real f = warp_sum(fp);
          const real fin = warp_sum(part);
          if (lane == 0) {
            const int k = n + j;
            const real qn = t == 0 ? qn0 : qn1, qbo = t == 0 ? qbo0 : qbo1, dvb = t == 0 ? dvb0 : dvb1;
            const real nv = FV[k] + BV[k] * f;
            modal = modal + f * (qn - qbo) - R(0.5) * (nv * nv - dvb * dvb + OM2[k] * (qn * qn - qbo * qbo))
              + R(0.5) * (nv * nv - fin * fin);
            V[k] = fin;
            if (!damp) store_out(TRACE, step * r + j, fin);
          }
        }
      }
    }
    GRID_SYNC();
    if (damp) {
      // D_b (second damper call)
      if (is_group) {
        if (mp >= 0) { SUPERPOSE(0, V + n); CHAIN(1); }
      } else if (is_board) {
        if (j0 < r) BOARD_DAMP(0, V + n, VBT, j0);
        if (j1 < r) BOARD_DAMP(0, V + n, VBT, j1);
      }
      GRID_SYNC();
      // D_c; the next step's first damper call continues from these station speeds
      if (is_group) {
        if (mp >= 0) { SUPERPOSE(1, VBT); if (step + 1 < steps) CHAIN(0); }
        if (step + 1 == steps) APPLY();
      } else if (is_board) {
        if (j0 < r) { BOARD_DAMP(1, VBT, V + n, j0); if (lane == 0) store_out(TRACE, step * r + j0, V[n + j0]); }
        if (j1 < r) { BOARD_DAMP(1, VBT, V + n, j1); if (lane == 0) store_out(TRACE, step * r + j1, V[n + j1]); }
      }
      GRID_SYNC();
    }
  }
  if (is_group) {
    BLOCK_SYNC_TG();
    UNROLL
    for (int k = bs + tid; k < be; k += NT) V[k] = tv[k - bs];
  }
  tot[0] = modal; tot[1] = dloss;
  block_sum(tot, 2, sh, tid);
  if (tid == 0) { store_out(LOSS, 2 * blk, tot[0]); store_out(LOSS, 2 * blk + 1, tot[1]); }
}
