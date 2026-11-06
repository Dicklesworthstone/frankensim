//! Sparse Krylov exponential propagator for DG ducts above the dense cap.
//!
//! The dense propagator ([`super::DgDuct::propagator`]) forms `phi0 = e^{dt L}` and
//! `dt phi1` as full matrices, so its build is O(n^3) and every step O(n^2); that is why it
//! stops at [`super::MAX_DG_PROPAGATOR_DOFS`]. The generator `L` of a DG bore with wall memory
//! is very sparse (a banded acoustic block plus one arrow per wall channel), so this kernel
//! never forms an exponential matrix. For the current state `x` it builds an Arnoldi basis
//! `V_m` of `K_m(tau L, x)` and reads
//!
//! ```text
//! e^{tau L} x = beta V_m e^{H_m} e_1,   phi_1(tau L) x = beta V_m phi_1(H_m) e_1,
//! phi_2(tau L) x = beta V_m phi_2(H_m) e_1,
//! ```
//!
//! with all three `phi_j(H_m) e_1` from ONE small exponential of the augmented matrix
//! `[[H, e1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]]` (Sidje 1998). The dense kernel's
//! per-step products are exactly these: `next = e^{dt L} x + sum u c1`, `int = dt phi1 x + sum u c2`,
//! with `c1 = dt phi1(dt L) b`, `c2 = dt^2 phi2(dt L) b` per port. The port columns do not depend on
//! the state, so they are formed once at build.
//!
//! Certificate. The stopping rule is the a-posteriori estimate `beta h_{m+1,m} |e_m^T phi_{j+1}(H_m) e_1|`
//! for `phi_j` (Saad 1992; Sidje 1998), times a safety factor, held under a fraction of the
//! card's tolerance at every step. At build the same estimate is bounded on the port vectors and on
//! fixed probe states, together with the change between `s` and `2s` substeps of the same step. That is an
//! estimate and not an enclosure, as is the dense kernel's halving check. The semigroup is
//! contractive in the energy norm, so per-step errors add at most linearly.
//!
//! Full state and the energy ledger are kept; nothing is reduced.

use super::*;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

/// Each step's estimated relative error is held under this fraction of the certificate ceiling,
/// after the safety factor, so that the independent halving check at build lands under the ceiling.
const STEP_TOL_FRACTION: f64 = 0.25;
/// The a-posteriori estimate is the leading term of the error series; the next term is not bounded.
const ESTIMATE_SAFETY: f64 = 2.0;
/// Arnoldi vectors a state step may use per substep before the build doubles the substeps.
const STATE_MAX_DIM: usize = 80;
/// Arnoldi vectors the one-time port columns may use.
const PORT_MAX_DIM: usize = 400;
/// Largest substep count the build will try (powers of two).
const MAX_SUBSTEPS: u32 = 64;
/// Fixed pseudo-random probe states beside the port vectors in the build certificate.
const PROBE_COUNT: usize = 3;

/// Measured use of the Krylov kernel since it was built.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgKrylovReport {
    /// Total DOFs of the generator.
    pub dofs: usize,
    /// Nonzero entries of the sparse generator.
    pub nonzeros: usize,
    /// Substeps per time step chosen at build (a power of two).
    pub substeps: u32,
    /// State solves done (each is `substeps` Arnoldi runs).
    pub solves: u64,
    /// Mean Arnoldi dimension per run.
    pub mean_dimension: f64,
    /// Largest Arnoldi dimension of any run.
    pub max_dimension: usize,
    /// Largest relative error estimate, after the safety factor, at the build certificate.
    pub build_estimate: f64,
}

#[derive(Debug, Default)]
struct KrylovStats {
    solves: AtomicU64,
    runs: AtomicU64,
    dim_sum: AtomicU64,
    dim_max: AtomicUsize,
    last_dim: AtomicUsize,
}

/// Row-compressed copy of the dense generator with its exact zeros dropped. No value changes.
#[derive(Debug, Clone)]
pub(super) struct SparseGenerator {
    n: usize,
    row_start: Vec<usize>,
    col: Vec<u32>,
    val: Vec<f64>,
}

impl SparseGenerator {
    fn from_dense(l: &[f64], n: usize) -> Result<Self, DuctError> {
        if l.len() != n * n || u32::try_from(n).is_err() {
            return Err(DuctError::BadParameter {
                what: "DG sparse generator needs a square operator with a representable size",
            });
        }
        let mut row_start = Vec::with_capacity(n + 1);
        row_start.push(0);
        let mut col = Vec::new();
        let mut val = Vec::new();
        for i in 0..n {
            for (j, &v) in l[i * n..(i + 1) * n].iter().enumerate() {
                if v != 0.0 {
                    if !v.is_finite() {
                        return Err(DuctError::BadParameter {
                            what: "DG sparse generator entry is not finite",
                        });
                    }
                    col.push(j as u32);
                    val.push(v);
                }
            }
            row_start.push(col.len());
        }
        Ok(Self { n, row_start, col, val })
    }

    /// `y = scale * (L x)`.
    fn apply(&self, x: &[f64], y: &mut [f64], scale: f64) {
        for i in 0..self.n {
            let mut acc = 0.0;
            for k in self.row_start[i]..self.row_start[i + 1] {
                acc += self.val[k] * x[self.col[k] as usize];
            }
            y[i] = scale * acc;
        }
    }
}

/// Arnoldi basis storage, reused between steps (a step is sequential; clones start empty).
#[derive(Debug, Default)]
struct Scratch(Mutex<Vec<f64>>);

impl Clone for Scratch {
    fn clone(&self) -> Self {
        Self::default()
    }
}

/// One accepted solve for one exact state: `next = e^{dt L} x`, `int = dt phi1(dt L) x`.
#[derive(Debug)]
pub(super) struct Solved {
    dt_bits: u64,
    tilde: Vec<f64>,
    pub(super) next: Vec<f64>,
    pub(super) int: Vec<f64>,
}

/// A one-entry cache of the last solve, keyed on the dt and the bits of the state. It is shared
/// by clones of a state (so a trial step on a clone reuses the response solve) and is never part
/// of equality.
#[derive(Debug, Default)]
pub(super) struct SolveCache(Mutex<Option<Arc<Solved>>>);

impl Clone for SolveCache {
    fn clone(&self) -> Self {
        let held = self.0.lock().unwrap_or_else(PoisonError::into_inner).clone();
        Self(Mutex::new(held))
    }
}

impl PartialEq for SolveCache {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

/// Outcome of one Arnoldi run.
struct PhiAction {
    /// `e^{tau L} v`.
    y0: Vec<f64>,
    /// `phi_1(tau L) v`.
    y1: Vec<f64>,
    /// `phi_2(tau L) v`, when asked for.
    y2: Vec<f64>,
    converged: bool,
    /// Relative error estimates (to `beta`) of `y0`, `y1`, `y2` at exit.
    estimate: f64,
}

#[derive(Debug, Clone)]
pub(super) struct KrylovKernel {
    generator: SparseGenerator,
    dt: f64,
    /// The certificate ceiling the propagator was built under.
    tol: f64,
    substeps: u32,
    pub(super) ports: Vec<(PortId, Vec<f64>)>,
    pub(super) port_cols_phi1: Vec<Vec<f64>>,
    pub(super) port_cols_phi2: Vec<Vec<f64>>,
    /// The state-independent mean-pressure impedance of every port pair, row-major,
    /// `sum_i b_row[i] port_cols_phi2[col][i] / dt`, formed once at build.
    pub(super) impedance: Vec<f64>,
    build_estimate: f64,
    stats: Arc<KrylovStats>,
    scratch: Scratch,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = [0.0_f64; 4];
    let chunks = a.len() / 4;
    for c in 0..chunks {
        for k in 0..4 {
            acc[k] += a[4 * c + k] * b[4 * c + k];
        }
    }
    let mut tail = 0.0;
    for i in 4 * chunks..a.len() {
        tail += a[i] * b[i];
    }
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + tail
}

fn norm2(a: &[f64]) -> f64 {
    det::sqrt(dot(a, a))
}

/// `phi_1(H) e_1`, `phi_2(H) e_1`, `phi_3(H) e_1` and `e^H e_1` for the leading `m x m` block of the
/// `(max_dim + 1) x max_dim` Hessenberg `h`, from one exponential of the augmented matrix.
fn small_phi(h: &[f64], ld: usize, m: usize) -> Result<[Vec<f64>; 4], DuctError> {
    let d = m + 3;
    let mut a = vec![0.0; d * d];
    for i in 0..m {
        for j in 0..m {
            a[i * d + j] = h[i * ld + j];
        }
    }
    a[m] = 1.0; // e_1 in column m
    a[m * d + m + 1] = 1.0;
    a[(m + 1) * d + m + 2] = 1.0;
    let shift = series_shift(&a, d, 1.0);
    let (e, _, _, _) = series_phis(&a, d, 1.0, shift)?;
    let col = |c: usize| -> Vec<f64> { (0..m).map(|i| e[i * d + c]).collect() };
    Ok([col(0), col(m), col(m + 1), col(m + 2)])
}

impl KrylovKernel {
    /// Arnoldi run on `v` for the step length `tau`: the three phi-actions and their estimates.
    ///
    /// `tol_rel` is the relative error (to `|v|`) every estimate must stay under after the safety factor;
    /// `want_y2` asks for `phi_2(tau L) v` (the port columns).
    fn exp_action(
        &self,
        v: &[f64],
        tau: f64,
        max_dim: usize,
        tol_rel: f64,
        want_y2: bool,
    ) -> Result<PhiAction, DuctError> {
        let n = self.generator.n;
        let beta = norm2(v);
        if !beta.is_finite() {
            return Err(DuctError::BadParameter {
                what: "DG Krylov input is not finite",
            });
        }
        if beta == 0.0 {
            return Ok(PhiAction {
                y0: vec![0.0; n],
                y1: vec![0.0; n],
                y2: vec![0.0; n],
                converged: true,
                estimate: 0.0,
            });
        }
        let ld = max_dim;
        let mut hess = vec![0.0; (max_dim + 1) * ld];
        let mut basis_guard = self.scratch.0.lock().unwrap_or_else(PoisonError::into_inner);
        let need = (max_dim + 2) * n;
        if basis_guard.len() < need {
            basis_guard.resize(need, 0.0);
        }
        let basis = &mut basis_guard[..need];
        for i in 0..n {
            basis[i] = v[i] / beta;
        }
        let hint = self.stats.last_dim.load(Ordering::Relaxed);
        let mut next_check = hint.saturating_sub(3).max(4);
        let mut m = 0;
        let mut coeffs: Option<[Vec<f64>; 4]> = None;
        let mut estimate = f64::INFINITY;
        let mut converged = false;
        for j in 0..max_dim {
            let (head, tail) = basis.split_at_mut((j + 1) * n);
            let vj = &head[j * n..(j + 1) * n];
            let w = &mut tail[..n]; // the slot of V[j + 1], used as the work vector
            self.generator.apply(vj, w, tau);
            let before = norm2(w);
            for i in 0..=j {
                let vi = &head[i * n..(i + 1) * n];
                let c = dot(vi, w);
                hess[i * ld + j] = c;
                for k in 0..n {
                    w[k] -= c * vi[k];
                }
            }
            let mut after = norm2(w);
            if after < 0.7 * before {
                for i in 0..=j {
                    let vi = &head[i * n..(i + 1) * n];
                    let c = dot(vi, w);
                    hess[i * ld + j] += c;
                    for k in 0..n {
                        w[k] -= c * vi[k];
                    }
                }
                after = norm2(w);
            }
            hess[(j + 1) * ld + j] = after;
            m = j + 1;
            let breakdown = !(after > 1e-14 * before.max(1e-300));
            if m >= next_check || breakdown || m == max_dim {
                let h_next = if breakdown { 0.0 } else { after };
                let c = small_phi(&hess, ld, m)?;
                // leading error terms: beta h |e_m^T phi_{j+1}(H) e_1| (relative: beta factored out)
                let e0 = h_next * c[1][m - 1].abs();
                let e1 = h_next * c[2][m - 1].abs();
                let e2 = h_next * c[3][m - 1].abs();
                estimate = ESTIMATE_SAFETY * if want_y2 { e0.max(e1).max(e2) } else { e0.max(e1) };
                coeffs = Some(c);
                if estimate <= tol_rel || breakdown {
                    converged = estimate <= tol_rel || breakdown;
                    break;
                }
                next_check = m + 2;
            }
            if m == max_dim {
                break;
            }
            for k in 0..n {
                w[k] /= after;
            }
        }
        let Some(c) = coeffs else {
            return Err(DuctError::BadParameter {
                what: "DG Krylov run produced no coefficients",
            });
        };
        let combine = |col: &[f64]| -> Vec<f64> {
            let mut y = vec![0.0; n];
            for (j, &cj) in col.iter().enumerate() {
                let vj = &basis[j * n..(j + 1) * n];
                let s = beta * cj;
                for k in 0..n {
                    y[k] += s * vj[k];
                }
            }
            y
        };
        let y0 = combine(&c[0]);
        let y1 = combine(&c[1]);
        let y2 = if want_y2 { combine(&c[2]) } else { Vec::new() };
        self.stats.runs.fetch_add(1, Ordering::Relaxed);
        self.stats.dim_sum.fetch_add(m as u64, Ordering::Relaxed);
        self.stats.dim_max.fetch_max(m, Ordering::Relaxed);
        self.stats.last_dim.store(m, Ordering::Relaxed);
        Ok(PhiAction { y0, y1, y2, converged, estimate })
    }

    /// `(e^{dt L} x, dt phi_1(dt L) x)` over `substeps` equal substeps, and the largest estimate used.
    fn solve_with(&self, x: &[f64], substeps: u32) -> Result<(Vec<f64>, Vec<f64>, f64), DuctError> {
        let n = self.generator.n;
        let tau = self.dt / f64::from(substeps);
        let tol_rel = STEP_TOL_FRACTION * self.tol;
        let mut cur = x.to_vec();
        let mut int = vec![0.0; n];
        let mut worst = 0.0_f64;
        for _ in 0..substeps {
            let act = self.exp_action(&cur, tau, STATE_MAX_DIM, tol_rel, false)?;
            worst = worst.max(act.estimate);
            if !act.converged {
                return Err(DuctError::PropagatorCertificate {
                    residual: act.estimate,
                    tol: self.tol,
                });
            }
            for i in 0..n {
                int[i] += tau * act.y1[i];
            }
            cur = act.y0;
        }
        self.stats.solves.fetch_add(1, Ordering::Relaxed);
        Ok((cur, int, worst))
    }

    /// The solve for `state`, reusing the cached one when the dt and every bit of the state match.
    pub(super) fn solve_cached(&self, dt: f64, state: &DgState) -> Result<Arc<Solved>, DuctError> {
        {
            let held = state.solve_cache.0.lock().unwrap_or_else(PoisonError::into_inner);
            if let Some(hit) = held.as_ref() {
                if hit.dt_bits == dt.to_bits()
                    && hit.tilde.len() == state.tilde.len()
                    && hit.tilde.iter().zip(&state.tilde).all(|(a, b)| a.to_bits() == b.to_bits())
                {
                    return Ok(Arc::clone(hit));
                }
            }
        }
        let (next, int, _) = self.solve_with(&state.tilde, self.substeps)?;
        let solved = Arc::new(Solved {
            dt_bits: dt.to_bits(),
            tilde: state.tilde.clone(),
            next,
            int,
        });
        *state.solve_cache.0.lock().unwrap_or_else(PoisonError::into_inner) = Some(Arc::clone(&solved));
        Ok(solved)
    }

    pub(super) fn report(&self) -> DgKrylovReport {
        let runs = self.stats.runs.load(Ordering::Relaxed);
        DgKrylovReport {
            dofs: self.generator.n,
            nonzeros: self.generator.val.len(),
            substeps: self.substeps,
            solves: self.stats.solves.load(Ordering::Relaxed),
            mean_dimension: if runs == 0 {
                0.0
            } else {
                self.stats.dim_sum.load(Ordering::Relaxed) as f64 / runs as f64
            },
            max_dimension: self.stats.dim_max.load(Ordering::Relaxed),
            build_estimate: self.build_estimate,
        }
    }
}

/// Deterministic probe states: unit vectors from a fixed splitmix64 stream.
fn probe_states(n: usize) -> Vec<Vec<f64>> {
    let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        seed = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = seed;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    };
    (0..PROBE_COUNT)
        .map(|_| {
            let mut v: Vec<f64> = (0..n).map(|_| (next() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0).collect();
            let norm = norm2(&v);
            for x in &mut v {
                *x /= norm;
            }
            v
        })
        .collect()
}

impl<M: DgMedium> DgDuct<M> {
    /// The sparse Krylov propagator ([`Self::propagator`] above the dense cap).
    pub(super) fn build_krylov_propagator(&self, dt: f64, tol: f64) -> Result<DgPropagator, DuctError> {
        let n = self.total_dofs();
        let generator = SparseGenerator::from_dense(&self.l_tilde, n)?;
        let mut kernel = KrylovKernel {
            generator,
            dt,
            tol,
            substeps: 1,
            ports: self.ports.clone(),
            port_cols_phi1: Vec::new(),
            port_cols_phi2: Vec::new(),
            impedance: Vec::new(),
            build_estimate: 0.0,
            stats: Arc::new(KrylovStats::default()),
            scratch: Scratch::default(),
        };
        let mut probes = probe_states(n);
        for (_, b) in &self.ports {
            let norm = norm2(b);
            if norm > 0.0 {
                probes.push(b.iter().map(|v| v / norm).collect());
            }
        }
        // The smallest power-of-two substep count at which every probe converges within STATE_MAX_DIM
        // vectors, with the independent check: the change between s and 2s substeps stays under tol.
        let mut certificate;
        let mut substeps = 1u32;
        loop {
            certificate = 0.0_f64;
            let mut failed = false;
            for x in &probes {
                match (kernel.solve_with(x, substeps), kernel.solve_with(x, 2 * substeps)) {
                    (Ok((coarse, _, e1)), Ok((fine, _, e2))) => {
                        let diff: f64 = coarse.iter().zip(&fine).map(|(a, b)| (a - b) * (a - b)).sum::<f64>();
                        certificate = certificate.max(e1).max(e2).max(det::sqrt(diff) / norm2(x));
                    }
                    (Err(DuctError::PropagatorCertificate { .. }), _) | (_, Err(DuctError::PropagatorCertificate { .. })) => {
                        failed = true;
                        break;
                    }
                    (Err(e), _) | (_, Err(e)) => return Err(e),
                }
            }
            if !failed {
                break;
            }
            if substeps >= MAX_SUBSTEPS {
                return Err(DuctError::PropagatorCertificate {
                    residual: f64::INFINITY,
                    tol,
                });
            }
            substeps *= 2;
        }
        kernel.substeps = substeps;
        if !(certificate <= tol) {
            return Err(DuctError::PropagatorCertificate {
                residual: certificate,
                tol,
            });
        }
        kernel.build_estimate = certificate;
        // The held-input columns, once: c1 = dt phi_1(dt L) b and c2 = dt^2 phi_2(dt L) b per port.
        let mut cols1 = Vec::with_capacity(self.ports.len());
        let mut cols2 = Vec::with_capacity(self.ports.len());
        for (_, b) in &self.ports {
            let act = kernel.exp_action(b, dt, PORT_MAX_DIM, STEP_TOL_FRACTION * tol, true)?;
            if !act.converged {
                return Err(DuctError::PropagatorCertificate {
                    residual: act.estimate,
                    tol,
                });
            }
            cols1.push(act.y1.iter().map(|v| dt * v).collect::<Vec<f64>>());
            cols2.push(act.y2.iter().map(|v| dt * dt * v).collect::<Vec<f64>>());
            certificate = certificate.max(act.estimate);
        }
        kernel.port_cols_phi1 = cols1;
        kernel.port_cols_phi2 = cols2;
        // The port-pair impedance does not depend on the state: the responses used to re-form it on
        // every call (`count^2 n` multiply-adds), with these same loops in this same order.
        let count = self.ports.len();
        let mut impedance = vec![0.0; count * count];
        for (row, (_, b)) in self.ports.iter().enumerate() {
            for i in 0..n {
                for col in 0..count {
                    impedance[row * count + col] += b[i] * kernel.port_cols_phi2[col][i];
                }
            }
            for col in 0..count {
                impedance[row * count + col] /= dt;
            }
        }
        kernel.impedance = impedance;
        kernel.build_estimate = certificate;
        if !(certificate <= tol) {
            return Err(DuctError::PropagatorCertificate {
                residual: certificate,
                tol,
            });
        }
        // The build probes and the port columns are not use of the kernel: keep only the dimension hint.
        let hint = kernel.stats.last_dim.load(Ordering::Relaxed);
        kernel.stats = Arc::new(KrylovStats::default());
        kernel.stats.last_dim.store(hint, Ordering::Relaxed);
        Ok(DgPropagator {
            operator_identity: self.operator_identity,
            dt,
            tol,
            certificate_residual: certificate,
            field_dofs: self.field_dofs,
            total_dofs: n,
            elements: self.elements.clone(),
            order: self.mesh.order,
            kernel: PropagatorKernel::Krylov(kernel),
        })
    }
}

impl DgPropagator {
    /// One step of the Krylov kernel: the same products as the dense series step, from one
    /// Arnoldi solve of the current state (shared with [`Self::multiport_response`] through the
    /// state's cache) plus the precomputed port columns.
    pub(super) fn step_krylov(
        &self,
        state: &mut DgState,
        port_inputs: &[PortInput],
        kernel: &KrylovKernel,
    ) -> Result<DgStepOutput, DuctError> {
        let ids: Vec<PortId> = kernel.ports.iter().map(|p| p.0).collect();
        let efforts = resolve_efforts(port_inputs, &ids)?;
        let e_prev = state.stored_energy_j();
        let solved = kernel.solve_cached(self.dt, state)?;
        let mut next_tilde = solved.next.clone();
        let mut int_tilde = solved.int.clone();
        for (slot, &u) in efforts.iter().enumerate() {
            if u != 0.0 {
                for i in 0..self.total_dofs {
                    next_tilde[i] += kernel.port_cols_phi1[slot][i] * u;
                    int_tilde[i] += kernel.port_cols_phi2[slot][i] * u;
                }
            }
        }
        self.finish_series_step(state, &kernel.ports, &efforts, e_prev, next_tilde, &int_tilde)
    }

    /// Measured use of the Krylov kernel (dimension, substeps, nonzeros), or `None` on the dense and
    /// lossless kernels.
    #[must_use]
    pub fn krylov_report(&self) -> Option<DgKrylovReport> {
        match &self.kernel {
            PropagatorKernel::Krylov(kernel) => Some(kernel.report()),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests;
