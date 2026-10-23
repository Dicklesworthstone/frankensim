//! Hysteretic (constant loss factor) material loss of an assembled plate,
//! projected ONCE on the complete certified source modes.
//!
//! Wood's loss factor is roughly independent of frequency in the measured
//! audio bands, so a material share c with loss factor `eta_c` adds the loss
//! stiffness `eta_c K_c` to the elastic stiffness `K = sum_c K_c`. This module
//! splits the plate's elastic stiffness into named shares (orthotropic panel
//! bending components, and stiffener groups as differences of two assemblies),
//! forms each share's modal strain energy `E_c = Phi^T K_c Phi` on the same
//! mass-normalized source modes, and combines them with per-mode loss factors
//! into `H = sum_c sqrt(eta_c,a eta_c,b) E_c` (diagonal `sum_c eta_c,a E_c,aa`).
//!
//! The reduction consumes the diagonal: mode a's loss factor
//! `eta_a = H_aa / omega_a^2` and source damping `c_a = H_aa / omega_a` [1/s]
//! at its own frequency (the modal strain-energy approximation). The neglected
//! off-diagonal terms are measured, never dropped silently: per mode the
//! largest normalized coupling `|H_ab| / sqrt(H_aa H_bb)` and, for its nearest
//! lower and upper neighbours, `|H_ab| / |omega_a^2 - omega_b^2|`, the ratio that
//! decides whether the diagonal holds. Where it is not small,
//! [`ModalEnergyShares::cluster_exact`] solves each coupled cluster's complex
//! stiffness `diag(omega^2) + iH` densely, so the diagonal's error is measured
//! per mode rather than bounded by first-order theory.
//!
//! This is an offline projection over the source slice (at most
//! [`MAX_SHARE_MODES`] modes). It creates no runtime dashpot channels; the
//! Kelvin-Voigt tensor image and its 64-connection limit are unchanged.
use super::{AcousticRealizeError, invalid, owner, poll};
use fs_exec::CancelGate;
use fs_la::eigen_complex;
use fs_math::c64::C64;
use fs_plate::{PlateMesh, PlateModel};

/// The complete certified source slice a board reduction admits; not a
/// runtime connection count.
pub const MAX_SHARE_MODES: usize = 512;
/// Upper bound on the named stiffness shares of one projection.
pub const MAX_SHARES: usize = 64;
/// Upper bound on one coupled cluster's dense complex eigenproblem in
/// [`ModalEnergyShares::cluster_exact`] (one LU per eigenvector: cost m^4).
pub const MAX_CLUSTER_MODES: usize = 128;

/// One named share of an assembled plate stiffness.
#[derive(Clone, Copy, Debug)]
pub enum StiffnessShare<'a> {
    /// Panel bending tensors [N m], one row-major 3x3 tensor per mesh triangle
    /// in mesh coordinates (DKT curvature order `xx, yy, 2xy`). Zero tensors
    /// contribute nothing; a tensor need not be definite.
    Panel(&'a [[f64; 9]]),
    /// The stiffness of `plus` minus that of `minus`: two assemblies of the
    /// same mesh, supports and DOF map that differ only in this share (for
    /// example a stiffener group present and absent).
    Difference { plus: &'a PlateModel, minus: &'a PlateModel },
}

/// One orthotropic panel section split into its four material-axis bending
/// components, each rotated to mesh coordinates exactly as fs-plate rotates the
/// whole section (`PlateSection::orthotropic_plane_stress_at_angle`): the L
/// (`D11`), R (`D22`), in-plane shear LR (`D66`) and Poisson (`D12`) parts.
/// Their sum is that section's bending tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicBendingParts {
    pub longitudinal: [f64; 9],
    pub radial: [f64; 9],
    pub shear: [f64; 9],
    pub poisson: [f64; 9],
}

impl OrthotropicBendingParts {
    /// The same constants fs-plate takes: `E_L, E_R, nu_LR, G_LR` [Pa],
    /// thickness [m] and material axis 1 at `angle_rad` from mesh x toward y.
    ///
    /// # Errors
    /// Nonpositive or nonfinite moduli or thickness, a nonfinite angle or a
    /// nonpositive `1 - nu_LR nu_RL` refuse.
    pub fn new(e_l: f64, e_r: f64, nu_lr: f64, g_lr: f64, thickness: f64, angle_rad: f64)
        -> Result<Self, AcousticRealizeError> {
        if [e_l, e_r, g_lr, thickness].iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !nu_lr.is_finite() || !angle_rad.is_finite() {
            return Err(invalid("orthotropic bending parts need positive finite moduli and thickness"));
        }
        let denom = 1.0 - nu_lr * (nu_lr * e_r / e_l);
        if !denom.is_finite() || denom <= 0.0 {
            return Err(invalid("orthotropic bending parts need a positive 1 - nu_LR nu_RL"));
        }
        let h3 = thickness * thickness * thickness / 12.0;
        let (d11, d22, d12, d33) = (e_l * h3 / denom, e_r * h3 / denom, nu_lr * e_r * h3 / denom, g_lr * h3);
        // fs-plate keeps the aligned and the isotropic paths unrotated.
        let aligned = angle_rad == 0.0 || (e_l == e_r && g_lr == e_l / (2.0 * (1.0 + nu_lr)));
        let rotate = |d: [f64; 9]| -> [f64; 9] {
            if aligned { return d; }
            let (s, c) = angle_rad.sin_cos();
            let (cc, ss, cs) = (c * c, s * s, c * s);
            let t = [[cc, ss, cs], [ss, cc, -cs], [-2.0 * cs, 2.0 * cs, cc - ss]];
            let mut out = [0.0; 9];
            for i in 0..3 {
                for j in i..3 {
                    let mut value = 0.0;
                    for a in 0..3 { for b in 0..3 { value += t[a][i] * d[3 * a + b] * t[b][j]; } }
                    out[3 * i + j] = value;
                    out[3 * j + i] = value;
                }
            }
            out
        };
        Ok(Self {
            longitudinal: rotate([d11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            radial: rotate([0.0, 0.0, 0.0, 0.0, d22, 0.0, 0.0, 0.0, 0.0]),
            shear: rotate([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d33]),
            poisson: rotate([0.0, d12, 0.0, d12, 0.0, 0.0, 0.0, 0.0, 0.0]),
        })
    }
}

/// Modal strain energy of every named share on the same source modes:
/// `E_c[a][b] = phi_a^T K_c phi_b` [1/s^2] for mass-normalized `phi`.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalEnergyShares {
    names: Vec<String>,
    omega: Vec<f64>,
    energy: Vec<Vec<f64>>,
}

/// One neighbour's coupling `|H_ab| / |omega_a^2 - omega_b^2|`; infinite for a
/// repeated frequency with nonzero coupling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeighbourCoupling {
    pub mode: usize,
    pub ratio: f64,
}

/// The diagonal hysteretic image consumed by the reduction, with the measured
/// size of the off-diagonal terms it neglects.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticModalLoss {
    /// `eta_a = H_aa / omega_a^2`.
    pub loss_factor: Vec<f64>,
    /// `c_a = H_aa / omega_a` [1/s]: the diagonal viscous term at the mode's
    /// own frequency.
    pub source_damping: Vec<f64>,
    /// `max_{b != a} |H_ab| / sqrt(H_aa H_bb)` (zero for a single mode).
    pub max_coupling: Vec<f64>,
    /// The nearest lower and upper modes in frequency, when they exist.
    pub neighbours: Vec<[Option<NeighbourCoupling>; 2]>,
}

/// The coupling the diagonal image neglects, resolved exactly where it is not
/// small. Modes joined by any pair whose `|H_ab| / |omega_a^2 - omega_b^2|`
/// reaches `pair_ratio` form clusters: connected components over ALL pairs, so
/// every pair between two clusters lies below it. Each cluster's complex
/// stiffness `diag(omega^2) + iH` (mass-normalized modes) is solved densely.
#[derive(Clone, Debug, PartialEq)]
pub struct ClusterExactLoss {
    /// The pair ratio that joins two modes.
    pub pair_ratio: f64,
    /// Every cluster of two or more modes, ascending within and by first mode.
    pub clusters: Vec<Vec<usize>>,
    /// Per mode: the index of its cluster in `clusters`, if any.
    pub cluster_of: Vec<Option<usize>>,
    /// Per mode: `[Re, Im]` [1/s^2] of the cluster eigenvalue matched to it (by
    /// the largest eigenvector overlap); the diagonal entry `[omega_a^2, H_aa]`
    /// for a mode in no cluster.
    pub eigenvalue: Vec<[f64; 2]>,
    /// Per mode: `Im / Re` of that eigenvalue, the exact loss factor.
    pub loss_factor: Vec<f64>,
    /// Per mode: `|x_a|^2 / |x|^2` of the matched eigenvector (one outside
    /// clusters); near one half the cluster's modes are fully mixed.
    pub overlap: Vec<f64>,
    /// Per mode: the sum, over every pair left between clusters, of the
    /// relative change that pair alone makes to the mode's loss factor at second
    /// order (`lambda_a = d_a - H_ab^2 / (d_a - d_b)`, `d = omega^2 + iH_aa`): an
    /// estimate of what the clustering leaves out.
    pub neglected: Vec<f64>,
    /// The largest pair ratio between two modes in different clusters (zero
    /// when there is no such pair).
    pub max_outside_ratio: f64,
}

/// Project every share's modal strain energy once on the supplied source modes.
///
/// # Errors
/// An empty or oversized basis, misaligned or nonfinite modes or frequencies,
/// unnamed or repeated shares, a panel share with the wrong triangle count or
/// nonfinite entries, a difference whose assemblies do not share the model's
/// DOF map, a degenerate element and cancellation refuse.
pub fn project_energy_shares(mesh: &PlateMesh, model: &PlateModel, shares: &[(&str, StiffnessShare<'_>)],
    phi: &[Vec<f64>], omega: &[f64], gate: &CancelGate) -> Result<ModalEnergyShares, AcousticRealizeError> {
    let n = phi.len();
    let free = model.free;
    if n == 0 || n > MAX_SHARE_MODES || omega.len() != n {
        return Err(invalid("hysteretic projection needs 1..=512 modes with one frequency each"));
    }
    if phi.iter().any(|v| v.len() != free || v.iter().any(|x| !x.is_finite()))
        || omega.iter().any(|w| !w.is_finite() || *w <= 0.0) {
        return Err(invalid("hysteretic projection needs finite reduced-DOF modes and positive frequencies"));
    }
    if shares.is_empty() || shares.len() > MAX_SHARES {
        return Err(invalid("hysteretic projection needs 1..=64 named stiffness shares"));
    }
    let mut names: Vec<String> = Vec::with_capacity(shares.len());
    for (name, _) in shares {
        if name.trim().is_empty() || names.iter().any(|seen| seen.as_str() == *name) {
            return Err(invalid("every stiffness share needs a unique nonempty name"));
        }
        names.push((*name).to_string());
    }
    let mut energy = Vec::with_capacity(shares.len());
    for (_, share) in shares {
        poll(gate)?;
        // y_b = K_c phi_b on the reduced DOFs.
        let mut y = vec![vec![0.0; free]; n];
        match share {
            StiffnessShare::Panel(tensors) => {
                if tensors.len() != mesh.tris.len() || tensors.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(invalid("a panel share needs one finite tensor per mesh triangle"));
                }
                for (e, tri) in mesh.tris.iter().enumerate() {
                    poll(gate)?;
                    let d = &tensors[e];
                    if d.iter().all(|v| *v == 0.0) { continue; }
                    let (c, _) = fs_plate::dkt_stiffness(&tri.map(|v| mesh.nodes[v].0),
                        &tri.map(|v| mesh.nodes[v].1), d, e).map_err(owner)?;
                    let dofs: [Option<usize>; 9] = std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]);
                    for (b, out) in y.iter_mut().enumerate() {
                        let local = dofs.map(|v| v.map_or(0.0, |i| phi[b][i]));
                        for (a, row) in dofs.iter().enumerate() {
                            let Some(row) = *row else { continue };
                            out[row] += (0..9).map(|k| c[9 * a + k] * local[k]).sum::<f64>();
                        }
                    }
                }
            }
            StiffnessShare::Difference { plus, minus } => {
                if plus.free != free || minus.free != free || plus.dof_map != model.dof_map
                    || minus.dof_map != model.dof_map {
                    return Err(invalid("a stiffness difference needs two assemblies on the model's own DOF map"));
                }
                let mut less = vec![0.0; free];
                for (b, out) in y.iter_mut().enumerate() {
                    poll(gate)?;
                    plus.k.spmv(&phi[b], out);
                    minus.k.spmv(&phi[b], &mut less);
                    for (o, l) in out.iter_mut().zip(&less) { *o -= l; }
                }
            }
        }
        let mut matrix = vec![0.0; n * n];
        for a in 0..n {
            poll(gate)?;
            for b in a..n {
                let value: f64 = phi[a].iter().zip(&y[b]).map(|(p, q)| p * q).sum();
                if !value.is_finite() { return Err(invalid("hysteretic share energy overflows")); }
                matrix[n * a + b] = value;
                matrix[n * b + a] = value;
            }
        }
        energy.push(matrix);
    }
    Ok(ModalEnergyShares { names, omega: omega.to_vec(), energy })
}

/// [`project_energy_shares`] for an offline caller with no cancellation source
/// (board generation): the projection is never interrupted.
///
/// # Errors
/// Exactly those of [`project_energy_shares`] other than cancellation.
pub fn project_energy_shares_offline(mesh: &PlateMesh, model: &PlateModel, shares: &[(&str, StiffnessShare<'_>)],
    phi: &[Vec<f64>], omega: &[f64]) -> Result<ModalEnergyShares, AcousticRealizeError> {
    project_energy_shares(mesh, model, shares, phi, omega, &CancelGate::new_clock_free())
}

impl ModalEnergyShares {
    /// Share names in projection order.
    #[must_use] pub fn names(&self) -> &[String] { &self.names }
    /// Source angular frequencies [rad/s].
    #[must_use] pub fn omega(&self) -> &[f64] { &self.omega }
    /// Number of source modes.
    #[must_use] pub fn modes(&self) -> usize { self.omega.len() }
    /// `phi_a^T K_c phi_b` [1/s^2] of share `share`.
    #[must_use] pub fn energy(&self, share: usize, a: usize, b: usize) -> f64 {
        self.energy[share][self.omega.len() * a + b]
    }
    /// Share `share`'s fraction of mode a's strain energy, `E_c,aa / omega_a^2`.
    #[must_use] pub fn fraction(&self, share: usize, a: usize) -> f64 {
        self.energy(share, a, a) / (self.omega[a] * self.omega[a])
    }
    /// The sum of every share's fraction of mode a: one when the shares
    /// partition the elastic stiffness that produced the modes.
    #[must_use] pub fn fraction_sum(&self, a: usize) -> f64 {
        (0..self.names.len()).map(|c| self.fraction(c, a)).sum()
    }

    /// Combine the shares with per-mode loss factors: `eta[c][a]` is share c's
    /// loss factor at mode a (zero for a share whose loss is omitted).
    ///
    /// # Errors
    /// A loss table that does not cover every share and mode, negative or
    /// nonfinite loss factors, and a non-dissipative or nonfinite diagonal refuse.
    pub fn hysteretic(&self, eta: &[Vec<f64>]) -> Result<HystereticModalLoss, AcousticRealizeError> {
        let n = self.omega.len();
        if eta.len() != self.names.len() || eta.iter().any(|row| row.len() != n
            || row.iter().any(|v| !v.is_finite() || *v < 0.0)) {
            return Err(invalid("hysteretic loss needs one nonnegative finite loss factor per share and mode"));
        }
        let h = |a: usize, b: usize| self.loss_stiffness(eta, a, b);
        let diagonal: Vec<f64> = (0..n).map(|a| h(a, a)).collect();
        if diagonal.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid("hysteretic material loss must be finite and dissipative on every mode"));
        }
        let mut max_coupling = vec![0.0_f64; n];
        for a in 0..n {
            for b in 0..n {
                if a == b { continue; }
                let off = h(a, b).abs();
                let scale = (diagonal[a] * diagonal[b]).sqrt();
                let ratio = if scale > 0.0 { off / scale } else if off > 0.0 { f64::INFINITY } else { 0.0 };
                max_coupling[a] = max_coupling[a].max(ratio);
            }
        }
        let neighbours = (0..n).map(|a| {
            let coupling = |b: usize| {
                let gap = (self.omega[a] * self.omega[a] - self.omega[b] * self.omega[b]).abs();
                let off = h(a, b).abs();
                NeighbourCoupling { mode: b,
                    ratio: if gap > 0.0 { off / gap } else if off > 0.0 { f64::INFINITY } else { 0.0 } }
            };
            let lower = (0..n).filter(|&b| b != a && self.omega[b] <= self.omega[a])
                .max_by(|&x, &y| self.omega[x].total_cmp(&self.omega[y]).then(x.cmp(&y)));
            let upper = (0..n).filter(|&b| b != a && self.omega[b] > self.omega[a])
                .min_by(|&x, &y| self.omega[x].total_cmp(&self.omega[y]).then(x.cmp(&y)));
            [lower.map(coupling), upper.map(coupling)]
        }).collect();
        Ok(HystereticModalLoss {
            loss_factor: (0..n).map(|a| diagonal[a] / (self.omega[a] * self.omega[a])).collect(),
            source_damping: (0..n).map(|a| diagonal[a] / self.omega[a]).collect(),
            max_coupling, neighbours,
        })
    }

    /// `H_ab = sum_c sqrt(eta_c,a eta_c,b) E_c,ab` [1/s^2]: exact when each
    /// share's loss factor is the same at both modes, a stated choice otherwise.
    fn loss_stiffness(&self, eta: &[Vec<f64>], a: usize, b: usize) -> f64 {
        (0..self.names.len()).map(|c| (eta[c][a] * eta[c][b]).sqrt() * self.energy(c, a, b)).sum()
    }

    /// Measure the diagonal image's error where its neglected coupling is not
    /// small: cluster the modes by `pair_ratio` over all pairs and solve each
    /// cluster's `diag(omega^2) + iH` densely. The eigenvalues come from fs-la's
    /// complex QR on the cluster matrix shifted by its mean `omega^2` and scaled
    /// by its largest remaining entry, so the solver's relative accuracy applies
    /// at the coupling's own scale; eigenvectors come from inverse iteration and
    /// match each mode to the eigenvalue whose eigenvector it dominates.
    ///
    /// # Errors
    /// Those of [`Self::hysteretic`], a nonpositive or nonfinite `pair_ratio`, a
    /// cluster above [`MAX_CLUSTER_MODES`] and a cluster eigenproblem that does
    /// not converge or is singular refuse.
    pub fn cluster_exact(&self, eta: &[Vec<f64>], pair_ratio: f64) -> Result<ClusterExactLoss, AcousticRealizeError> {
        let diagonal = self.hysteretic(eta)?;
        if !pair_ratio.is_finite() || pair_ratio <= 0.0 {
            return Err(invalid("cluster_exact needs a positive finite pair ratio"));
        }
        let n = self.omega.len();
        let w2: Vec<f64> = self.omega.iter().map(|w| w * w).collect();
        let mut h = vec![0.0; n * n];
        for a in 0..n {
            for b in a..n {
                let value = self.loss_stiffness(eta, a, b);
                h[n * a + b] = value;
                h[n * b + a] = value;
            }
        }
        let ratio = |a: usize, b: usize| {
            let (gap, off) = ((w2[a] - w2[b]).abs(), h[n * a + b].abs());
            if gap > 0.0 { off / gap } else if off > 0.0 { f64::INFINITY } else { 0.0 }
        };
        // Connected components over every pair; each root is its component's lowest mode.
        let mut root: Vec<usize> = (0..n).collect();
        let find = |root: &mut Vec<usize>, mut x: usize| {
            while root[x] != x { root[x] = root[root[x]]; x = root[x]; }
            x
        };
        for a in 0..n {
            for b in a + 1..n {
                if ratio(a, b) >= pair_ratio {
                    let (ra, rb) = (find(&mut root, a), find(&mut root, b));
                    if ra != rb { root[ra.max(rb)] = ra.min(rb); }
                }
            }
        }
        let label: Vec<usize> = (0..n).map(|a| find(&mut root, a)).collect();
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for a in 0..n {
            if label[a] == a { clusters.push(vec![a]); } else if let Some(c) = clusters.iter_mut().find(|c| c[0] == label[a]) { c.push(a); }
        }
        clusters.retain(|c| c.len() > 1);
        let mut cluster_of = vec![None; n];
        for (k, cluster) in clusters.iter().enumerate() { for &a in cluster { cluster_of[a] = Some(k); } }
        let mut max_outside_ratio = 0.0_f64;
        let mut neglected = vec![0.0; n];
        for a in 0..n {
            let (re, im) = (w2[a], h[n * a + a]);
            for b in 0..n {
                if b == a || (cluster_of[a].is_some() && cluster_of[a] == cluster_of[b]) { continue; }
                max_outside_ratio = max_outside_ratio.max(ratio(a, b));
                let off = h[n * a + b];
                if off == 0.0 { continue; }
                // lambda_a moves by delta; (im + delta.im)/(re + delta.re) - im/re without cancellation.
                let delta = -(C64::from_re(off * off) / C64::new(re - w2[b], im - h[n * b + b]));
                let change = (delta.im * re - im * delta.re).abs() / (re + delta.re);
                neglected[a] += if im > 0.0 { change / im } else if change > 0.0 { f64::INFINITY } else { 0.0 };
            }
        }
        let mut eigenvalue: Vec<[f64; 2]> = (0..n).map(|a| [w2[a], h[n * a + a]]).collect();
        let mut overlap = vec![1.0; n];
        for cluster in &clusters {
            let m = cluster.len();
            if m > MAX_CLUSTER_MODES {
                return Err(invalid("a coupled cluster exceeds MAX_CLUSTER_MODES; the diagonal check cannot resolve it densely"));
            }
            let mean = cluster.iter().map(|&a| w2[a]).sum::<f64>() / m as f64;
            let scale = cluster.iter().flat_map(|&a| cluster.iter().map(move |&b| (a, b)))
                .map(|(a, b)| if a == b { (w2[a] - mean).abs().max(h[n * a + a].abs()) } else { h[n * a + b].abs() })
                .fold(0.0_f64, f64::max);
            let matrix: Vec<C64> = (0..m * m).map(|ij| {
                let (a, b) = (cluster[ij / m], cluster[ij % m]);
                C64::new(if a == b { (w2[a] - mean) / scale } else { 0.0 }, h[n * a + b] / scale)
            }).collect();
            let values = eigen_complex::eig(&matrix, m).map_err(owner)?;
            let weights = values.iter().map(|&beta| {
                let x = inverse_iteration(&matrix, m, beta)?;
                let norm: f64 = x.iter().map(|v| v.norm_sq()).sum();
                Ok(x.iter().map(|v| v.norm_sq() / norm).collect::<Vec<f64>>())
            }).collect::<Result<Vec<_>, AcousticRealizeError>>()?;
            // Greedy assignment by largest overlap; ties to the lowest indices.
            let (mut open_value, mut open_mode) = (vec![true; m], vec![true; m]);
            for _ in 0..m {
                let mut best = (f64::NEG_INFINITY, 0, 0);
                for k in (0..m).filter(|&k| open_value[k]) {
                    for i in (0..m).filter(|&i| open_mode[i]) {
                        if weights[k][i] > best.0 { best = (weights[k][i], k, i); }
                    }
                }
                let (w, k, i) = best;
                open_value[k] = false;
                open_mode[i] = false;
                eigenvalue[cluster[i]] = [mean + scale * values[k].re, scale * values[k].im];
                overlap[cluster[i]] = w;
            }
        }
        if eigenvalue.iter().any(|e| !e.iter().all(|v| v.is_finite()) || e[0] <= 0.0) {
            return Err(invalid("a cluster eigenvalue is nonfinite or has a nonpositive stiffness"));
        }
        Ok(ClusterExactLoss { pair_ratio, clusters, cluster_of,
            loss_factor: eigenvalue.iter().map(|e| e[1] / e[0]).collect(), eigenvalue, overlap, neglected,
            max_outside_ratio })
    }
}

/// The eigenvector of the O(1) cluster matrix `b` at its eigenvalue `beta`:
/// three steps of inverse iteration from a fixed generic start, shifted 1e-10
/// off `beta` so the factorization exists.
fn inverse_iteration(b: &[C64], m: usize, beta: C64) -> Result<Vec<C64>, AcousticRealizeError> {
    let shift = beta + C64::new(1e-10, 1e-10);
    let shifted: Vec<C64> = (0..m * m).map(|ij| if ij / m == ij % m { b[ij] - shift } else { b[ij] }).collect();
    let lu = eigen_complex::lu_complex(&shifted, m).map_err(owner)?;
    let mut x: Vec<C64> = (0..m).map(|j| C64::new(1.0, 0.618_033_988_749_895 * (j + 1) as f64)).collect();
    for _ in 0..3 {
        lu.solve(&mut x);
        let norm = x.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt();
        if !norm.is_finite() || norm == 0.0 { return Err(invalid("cluster inverse iteration lost its vector")); }
        for v in &mut x { *v = v.scale(1.0 / norm); }
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_plate::{AssemblyOptions, EdgeSupport, PlateSection, SliceOptions, Stiffener};

    const E_L: f64 = 11.9e9;
    const E_R: f64 = 0.93e9;
    const NU: f64 = 0.37;
    const G_LR: f64 = 0.76e9;
    const H: f64 = 0.008;
    const RHO: f64 = 450.0;
    const ANGLE: f64 = 0.61;

    /// A clamped orthotropic panel with one eccentric stiffener, its source
    /// modes, and the same assembly without the stiffener and with its E doubled.
    struct Board { mesh: PlateMesh, model: PlateModel, bare: PlateModel, stiff: PlateModel,
        phi: Vec<Vec<f64>>, omega: Vec<f64> }
    fn board() -> Board {
        let (nx, ny) = (8, 6);
        let mesh = PlateMesh::rectangle(0.6, 0.45, nx, ny);
        let boundary = PlateMesh::rectangle_boundary(nx, ny);
        let section = PlateSection::orthotropic_plane_stress_at_angle(E_L, E_R, NU, G_LR, H, RHO, ANGLE).unwrap();
        let rib = |e: f64| Stiffener { nodes: (0..=nx).map(|i| 3 * (nx + 1) + i).collect(), e, g: 0.6e9,
            area: 0.025 * 0.02, inertia: 0.025 * 0.02_f64.powi(3) / 12.0, torsion: 3.0e-8,
            eccentricity: -0.014, density: 400.0 };
        let options = AssemblyOptions { support: EdgeSupport::Clamped, pretension: 0.0 };
        let assemble = |ribs: &[Stiffener]| fs_plate::assemble(&mesh, &section, &boundary, ribs, &options).unwrap();
        let model = assemble(&[rib(10e9)]);
        let report = fs_plate::modes(&model, (0.0, (std::f64::consts::TAU * 2500.0).powi(2)),
            &SliceOptions::default()).unwrap();
        assert!(report.modes.len() >= 6, "{} modes", report.modes.len());
        let phi = report.modes.iter().map(|pair| pair.phi.clone()).collect();
        let omega = report.modes.iter().map(|pair| pair.lambda.sqrt()).collect();
        Board { bare: assemble(&[]), stiff: assemble(&[rib(20e9)]), mesh, model, phi, omega }
    }
    fn parts(mesh: &PlateMesh) -> Vec<OrthotropicBendingParts> {
        vec![OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, ANGLE).unwrap(); mesh.tris.len()]
    }
    /// Panel L, R, LR, Poisson, rib bending and rib torsion: a partition of K.
    fn shares(b: &Board) -> ModalEnergyShares {
        let p = parts(&b.mesh);
        let pick = |f: fn(&OrthotropicBendingParts) -> [f64; 9]| p.iter().map(f).collect::<Vec<_>>();
        let (l, r, lr, nu) = (pick(|x| x.longitudinal), pick(|x| x.radial), pick(|x| x.shear), pick(|x| x.poisson));
        let gate = CancelGate::new_clock_free();
        // Bending is linear in E: (E doubled) - (E) is the rib's bending share;
        // the rest of (with) - (without) is its torsion.
        project_energy_shares(&b.mesh, &b.model, &[("panel-l", StiffnessShare::Panel(&l)),
            ("panel-r", StiffnessShare::Panel(&r)), ("panel-lr", StiffnessShare::Panel(&lr)),
            ("panel-12", StiffnessShare::Panel(&nu)),
            ("rib", StiffnessShare::Difference { plus: &b.model, minus: &b.bare }),
            ("rib-bending", StiffnessShare::Difference { plus: &b.stiff, minus: &b.model })],
            &b.phi, &b.omega, &gate).unwrap()
    }

    // G0: the four material-axis parts sum to fs-plate's own rotated section.
    #[test]
    fn g0_bending_parts_sum_to_the_fs_plate_section() {
        for angle in [0.0, ANGLE, -1.2] {
            let p = OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, angle).unwrap();
            let d = PlateSection::orthotropic_plane_stress_at_angle(E_L, E_R, NU, G_LR, H, RHO, angle).unwrap().d;
            let scale = d.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            for i in 0..9 {
                let sum = p.longitudinal[i] + p.radial[i] + p.shear[i] + p.poisson[i];
                assert!((sum - d[i]).abs() <= 1e-13 * scale, "angle {angle} entry {i}: {sum:e} vs {:e}", d[i]);
            }
        }
        assert!(OrthotropicBendingParts::new(E_L, E_R, 4.0, G_LR, H, 0.0).is_err());
        assert!(OrthotropicBendingParts::new(E_L, -1.0, NU, G_LR, H, 0.0).is_err());
        assert!(OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, f64::NAN).is_err());
    }

    // G0: the shares partition the stiffness that produced the modes, so every
    // mode's fractions sum to one; a uniform loss factor gives every mode that
    // loss factor; and a source damping is its loss factor times omega.
    #[test]
    fn g0_shares_partition_each_mode_and_a_uniform_loss_factor_is_returned() {
        let b = board();
        let e = shares(&b);
        let n = e.modes();
        // rib = rib-bending + rib-torsion; the partition uses rib-torsion = rib - rib-bending.
        for a in 0..n {
            let sum = e.fraction(0, a) + e.fraction(1, a) + e.fraction(2, a) + e.fraction(3, a) + e.fraction(4, a);
            assert!((sum - 1.0).abs() < 1e-9, "mode {a}: {sum}");
            assert!(e.fraction(5, a) >= 0.0 && e.fraction(5, a) <= e.fraction(4, a) + 1e-12, "mode {a}");
        }
        let eta = 0.013;
        let uniform: Vec<Vec<f64>> = (0..6).map(|c| vec![if c == 5 { 0.0 } else { eta }; n]).collect();
        let loss = e.hysteretic(&uniform).unwrap();
        for a in 0..n {
            assert!((loss.loss_factor[a] / eta - 1.0).abs() < 1e-9, "mode {a}: {}", loss.loss_factor[a]);
            assert!((loss.source_damping[a] / (loss.loss_factor[a] * b.omega[a]) - 1.0).abs() < 1e-12);
            // Uniform loss is proportional to K, which the modes diagonalize.
            assert!(loss.max_coupling[a] < 1e-6, "mode {a}: {}", loss.max_coupling[a]);
        }
    }

    // G1: changing one share's loss factor moves each mode's loss factor by
    // exactly that change times the mode's energy fraction in the share.
    #[test]
    fn g1_one_share_change_moves_each_mode_by_its_energy_fraction() {
        let b = board();
        let e = shares(&b);
        let n = e.modes();
        let base: Vec<Vec<f64>> = vec![vec![0.0124; n], vec![0.0188; n], vec![0.0155; n], vec![0.0; n],
            vec![0.0; n], vec![0.0; n]];
        let mut swapped = base.clone();
        swapped[0] = vec![0.0070; n];
        let (before, after) = (e.hysteretic(&base).unwrap(), e.hysteretic(&swapped).unwrap());
        for a in 0..n {
            let predicted = (0.0070 - 0.0124) * e.fraction(0, a);
            assert!(((after.loss_factor[a] - before.loss_factor[a]) - predicted).abs() < 1e-15 + 1e-12 * predicted.abs(),
                "mode {a}");
            assert!(after.loss_factor[a] < before.loss_factor[a], "mode {a}");
        }
        // Nonuniform loss is not proportional to K: coupling appears and is measured.
        assert!(before.max_coupling.iter().any(|c| *c > 1e-6));
        for (a, pair) in before.neighbours.iter().enumerate() {
            assert!(pair.iter().flatten().all(|c| c.mode != a && c.ratio.is_finite() && c.ratio >= 0.0));
        }
        assert!(before.neighbours[0][0].is_none() && before.neighbours[n - 1][1].is_none());
    }

    /// Two modes 3 Hz apart near 1 kHz split into two shares with different
    /// loss factors, so the loss stiffness couples them by `(0.02 - 0.005) q`.
    fn pair(q: f64) -> ModalEnergyShares {
        let omega = vec![std::f64::consts::TAU * 1000.0, std::f64::consts::TAU * 1003.0];
        let w2 = [omega[0] * omega[0], omega[1] * omega[1]];
        let one = vec![0.6 * w2[0], q, q, 0.3 * w2[1]];
        let two = vec![w2[0] - one[0], -q, -q, w2[1] - one[3]];
        ModalEnergyShares { names: vec!["one".into(), "two".into()], omega, energy: vec![one, two] }
    }
    fn pair_eta() -> Vec<Vec<f64>> { vec![vec![0.02; 2], vec![0.005; 2]] }
    /// Closed-form loss factors of the pair's complex stiffness
    /// `[[a, c], [c, d]]`, `lambda = (a + d)/2 +- sqrt(((a - d)/2)^2 + c^2)`,
    /// ordered so mode 0 takes the eigenvalue whose eigenvector `(c, lambda - a)`
    /// it dominates.
    fn closed_form(e: &ModalEnergyShares) -> [f64; 2] {
        let eta = pair_eta();
        let (h, w2) = (|a, b| e.loss_stiffness(&eta, a, b), |a: usize| e.omega[a] * e.omega[a]);
        let (a, d, c) = (C64::new(w2(0), h(0, 0)), C64::new(w2(1), h(1, 1)), C64::new(0.0, h(0, 1)));
        let half = (a - d).scale(0.5);
        let disc = (half * half + c * c).sqrt();
        let (plus, minus) = ((a + d).scale(0.5) + disc, (a + d).scale(0.5) - disc);
        let on_zero = |l: C64| c.norm_sq() / (c.norm_sq() + (l - a).norm_sq());
        let (mine, other) = if on_zero(plus) >= on_zero(minus) { (plus, minus) } else { (minus, plus) };
        [mine.im / mine.re, other.im / other.re]
    }

    // G0: a joined pair's exact loss factors are the closed-form eigenvalues of
    // its 2x2 complex stiffness, far from the diagonal when the pair ratio is
    // near one; a weak pair below the ratio stays diagonal and the second-order
    // estimate predicts what that leaves out; a smaller ratio joins it exactly.
    #[test]
    fn g0_cluster_exact_matches_the_closed_form_pair_and_estimates_what_it_leaves_out() {
        let eta = pair_eta();
        let strong = pair(1.58e7);
        let (diagonal, exact) = (strong.hysteretic(&eta).unwrap(), strong.cluster_exact(&eta, 0.03).unwrap());
        let r = diagonal.neighbours[0][1].unwrap().ratio;
        assert!(r > 0.5, "{r}");
        assert_eq!(exact.clusters, vec![vec![0, 1]]);
        assert_eq!(exact.cluster_of, vec![Some(0), Some(0)]);
        let closed = closed_form(&strong);
        for a in 0..2 {
            // fs-la's QR is good to about 1e-12 at the scaled cluster's unit scale.
            assert!((exact.loss_factor[a] / closed[a] - 1.0).abs() < 1e-9, "mode {a}: {} vs {}",
                exact.loss_factor[a], closed[a]);
            assert!(exact.overlap[a] > 0.5 && exact.overlap[a] <= 1.0, "mode {a}: {}", exact.overlap[a]);
        }
        assert!((exact.loss_factor[0] / diagonal.loss_factor[0] - 1.0).abs() > 0.05, "{exact:?}");
        assert_eq!((exact.neglected.clone(), exact.max_outside_ratio), (vec![0.0, 0.0], 0.0));

        let weak = pair(1.58e5);
        let (diagonal, exact) = (weak.hysteretic(&eta).unwrap(), weak.cluster_exact(&eta, 0.03).unwrap());
        let r = diagonal.neighbours[0][1].unwrap().ratio;
        assert!(r > 0.003 && r < 0.03, "{r}");
        assert!(exact.clusters.is_empty() && exact.cluster_of == vec![None, None]);
        assert_eq!((exact.loss_factor.clone(), exact.overlap.clone()), (diagonal.loss_factor.clone(), vec![1.0, 1.0]));
        assert_eq!(exact.max_outside_ratio, r);
        let closed = closed_form(&weak);
        for a in 0..2 {
            let left_out = (closed[a] / diagonal.loss_factor[a] - 1.0).abs();
            // A second-order estimate in the pair ratio: its own error is O(r) relative.
            assert!((exact.neglected[a] - left_out).abs() <= 3.0 * r * left_out, "mode {a}: {} vs {left_out}",
                exact.neglected[a]);
        }
        let joined = weak.cluster_exact(&eta, 0.001).unwrap();
        assert_eq!(joined.clusters, vec![vec![0, 1]]);
        for a in 0..2 { assert!((joined.loss_factor[a] / closed[a] - 1.0).abs() < 1e-9, "mode {a}"); }
        for bad in [0.0, -0.03, f64::NAN, f64::INFINITY] { assert!(weak.cluster_exact(&eta, bad).is_err(), "{bad}"); }
    }

    // G0: on the panel with nonuniform loss, the clusters partition the joined
    // modes, every pair between clusters lies below the ratio, each cluster's
    // eigenvalues keep its trace (sum omega^2 and sum H_aa), and a mode in no
    // cluster keeps its diagonal loss factor exactly.
    #[test]
    fn g0_cluster_eigenvalues_keep_each_clusters_trace_and_unjoined_modes_stay_diagonal() {
        let b = board();
        let e = shares(&b);
        let n = e.modes();
        // Panel L, R, LR, no Poisson loss, the whole rib at 0.03; rib-bending overlaps the rib.
        let eta = vec![vec![0.0124; n], vec![0.0188; n], vec![0.0155; n], vec![0.0; n], vec![0.03; n], vec![0.0; n]];
        let diagonal = e.hysteretic(&eta).unwrap();
        let top = diagonal.neighbours.iter().flatten().flatten().map(|c| c.ratio).fold(0.0_f64, f64::max);
        let exact = e.cluster_exact(&eta, top / 2.0).unwrap();
        assert!(!exact.clusters.is_empty(), "the strongest pair joins at half its ratio");
        assert!(exact.max_outside_ratio < exact.pair_ratio);
        let w2 = |a: usize| b.omega[a] * b.omega[a];
        let h = |a: usize| e.loss_stiffness(&eta, a, a);
        for (k, cluster) in exact.clusters.iter().enumerate() {
            assert!(cluster.len() >= 2 && cluster.windows(2).all(|w| w[0] < w[1]));
            assert!(cluster.iter().all(|&a| exact.cluster_of[a] == Some(k)));
            let (re, im): (f64, f64) = (cluster.iter().map(|&a| exact.eigenvalue[a][0]).sum(),
                cluster.iter().map(|&a| exact.eigenvalue[a][1]).sum());
            let (w2_sum, h_sum): (f64, f64) = (cluster.iter().map(|&a| w2(a)).sum(), cluster.iter().map(|&a| h(a)).sum());
            // The QR works at the coupling's scale (at most sum H_aa); 1e-9 of it is three
            // orders above its observed 1e-12 relative accuracy.
            assert!((re - w2_sum).abs() <= 1e-9 * h_sum && (im - h_sum).abs() <= 1e-9 * h_sum,
                "cluster {cluster:?}: {re} vs {w2_sum}, {im} vs {h_sum}");
            for &a in cluster {
                println!("cluster {k} mode {a} {:.2} Hz: diagonal eta {:.5}, exact {:.5}, overlap {:.3}",
                    b.omega[a] / std::f64::consts::TAU, diagonal.loss_factor[a], exact.loss_factor[a], exact.overlap[a]);
            }
        }
        for a in 0..n {
            assert!(exact.overlap[a] > 0.0 && exact.overlap[a] <= 1.0 && exact.neglected[a].is_finite()
                && exact.neglected[a] >= 0.0, "mode {a}");
            if exact.cluster_of[a].is_none() {
                assert_eq!(exact.loss_factor[a], diagonal.loss_factor[a], "mode {a}");
                assert_eq!(exact.eigenvalue[a], [w2(a), h(a)], "mode {a}");
            }
        }
    }

    // G0: malformed inputs refuse before any projection.
    #[test]
    fn g0_malformed_projection_inputs_refuse() {
        let b = board();
        let gate = CancelGate::new_clock_free();
        let p = parts(&b.mesh);
        let l: Vec<[f64; 9]> = p.iter().map(|x| x.longitudinal).collect();
        let panel = [("panel-l", StiffnessShare::Panel(&l))];
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &[], &[], &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega[1..], &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[], &b.phi, &b.omega, &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[("a", StiffnessShare::Panel(&l)),
            ("a", StiffnessShare::Panel(&l))], &b.phi, &b.omega, &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[("short", StiffnessShare::Panel(&l[1..]))],
            &b.phi, &b.omega, &gate).is_err());
        let mut bad = b.omega.clone();
        bad[0] = -1.0;
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &bad, &gate).is_err());
        let e = project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega, &gate).unwrap();
        let n = e.modes();
        assert!(e.hysteretic(&[vec![0.01; n - 1]]).is_err());
        assert!(e.hysteretic(&[vec![-0.01; n]]).is_err());
        assert!(e.hysteretic(&[vec![f64::NAN; n]]).is_err());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega, &cancelled).is_err());
    }
}
