//! The radiation loss of a free shell's own modes (`Body::attach_shell_radiation`, MB 14574 (d), a
//! correction): each rest mode of the shell up to a declared frequency takes the baffled Rayleigh
//! resistance of its own shape on the shell's projected plane, from both faces, times the microphone's
//! own first-order dipole factor (`Mic::unbaffled`), as the explicit low-rank loss the heads' measured
//! shape damping uses (`ShapeDamping`): C = sum_m 2 sigma_m (M phi_m)(M phi_m)'.
//!
//! Modes. K is the linear stiffness about the rest shape, `Body::rest_stiffness_apply` (bending,
//! supports and the von Karman curvature Hessian), the operator the clamp-share analysis (.372)
//! assembles. A body whose thickness, rest height and supports depend on r alone, on the centred
//! square grid, commutes with the reflections x -> -x and y -> -y, so M^-1/2 K M^-1/2 splits into four
//! parity sectors, each solved densely (fs-la `householder_eigh`). The sectors are assembled from the
//! operator applied to the symmetry-adapted basis, so K's asymmetry and the coupling between sectors
//! are measured on every entry, and both are gated.
//!
//! Resistance. With phi mass-normalised (phi' M phi = 1) and the surface velocity phi q', one face's
//! baffled Rayleigh integral radiates the mean power R |q'|^2 / 2, with
//!     R = (omega rho0 / 2 pi) sum_i sum_j phi_i phi_j sin(k r_ij) / r_ij h^4   (the i = j term k h^4),
//! the h^2 node weights the microphone's Rayleigh taps use. The double sum is a convolution on the
//! grid, done by FFT. Both faces radiate; an unbaffled sheet's low-ka cancellation takes the
//! microphone's own factor |H|^2 = (omega tau)^2 / (1 + (omega tau)^2), tau = 4a / (3 pi c). The mode's
//! loss is 2 sigma = 2 |H|^2 R, eta = 2 sigma / omega; the receipt prints, beside it, the value without
//! the dipole factor.
//!
//! Cross terms (MB 14957, a correction). Two loaded modes of one sector also couple through their
//! cross resistance R_mn, and the coupling moves their decay where it is not small against their
//! frequency gap. Their overlap is |2 sigma_mn| / |omega_m - omega_n|, 2 sigma_mn = 2 |H|^2 R_mn taken
//! at either mode's frequency (the larger). The pairs whose overlap exceeds the declared
//! `CLUSTER_OVERLAP` join into clusters, and each cluster takes its full block 2 |H|^2 R at its mean
//! frequency. sin(kr)/r is a positive definite kernel, so every sampled block is positive semidefinite;
//! the block must be so, and symmetric, to `BLOCK_GATE` or the attach refuses. It enters the low-rank
//! form through its eigenvectors, M-orthonormal combinations of its modes. Every other mode keeps its
//! own resistance only, and the receipt prints, beside the overlaps, the share sum_n |R_mn| / R_mm of
//! the cross resistance it drops, over the modes of its sector outside its cluster, at its own
//! frequency. Modes in different sectors do not couple through the kernel.
use super::*;
use fs_la::eigen::householder_eigh;
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Eigenvalues within this fraction of the largest are rigid motions (an unsupported plate's
/// translation and two tilts); they take no load.
const RIGID_FRACTION: f64 = 1e-10;

/// The gate on K's asymmetry (relative Frobenius) and on the largest coupling between parity
/// sectors (relative to K's largest entry), both on M^-1/2 K M^-1/2 in the symmetry-adapted basis.
/// The curvature Hessian's central difference leaves about 2e-15 (.372).
const SYMMETRY_GATE: f64 = 1e-12;

/// The overlap above which two loaded modes of one parity sector share their resistance block,
/// |2 sigma_mn| / |omega_m - omega_n| (MB 14957; declared in the plan line before any run).
pub(crate) const CLUSTER_OVERLAP: f64 = 0.05;

/// A cluster's resistance block must be symmetric to this fraction of its largest entry and its
/// eigenvalues at least minus this fraction of its largest; a negative eigenvalue within it is
/// rounding and is set to zero.
const BLOCK_GATE: f64 = 1e-9;

/// A body's rest modes up to a frequency.
pub(crate) struct RestModes {
    /// (omega rad/s, phi with phi' M phi = 1, parity sector 0..4), ascending in omega.
    pub(crate) modes: Vec<(f64, Vec<f64>, usize)>,
    pub(crate) rigid: usize,
    /// The grid's highest angular frequency [rad/s].
    pub(crate) omega_max: f64,
    pub(crate) asymmetry: f64,
    pub(crate) coupling: f64,
    pub(crate) sectors: [usize; 4],
}

/// The modes of M^-1/2 K M^-1/2, K = `rest_stiffness_apply`, with omega <= `omega_limit`, through
/// the four parity sectors of the centred grid (sector s has parities (+x+y, -x+y, +x-y, -x-y)[s]).
/// Refuses a node off the centred grid, a grid that is not reflection-symmetric, an asymmetric K, a
/// coupling between sectors, a sector too large for the dense solve and a negative eigenvalue.
pub(crate) fn rest_modes(body: &mut Body, omega_limit: f64) -> Result<RestModes, String> {
    let (n, h, name) = (body.n, body.h, body.name.clone());
    let cell = |v: f64| (v / h).round() as i64;
    let mut at = BTreeMap::new();
    for (i, p) in body.xy.iter().enumerate() {
        let (ci, cj) = (cell(p[0]), cell(p[1]));
        if (p[0] - ci as f64 * h).abs() > 1e-9 * h || (p[1] - cj as f64 * h).abs() > 1e-9 * h {
            return Err(format!("{name}: node {i} is off the centred grid"));
        }
        at.insert((ci, cj), i);
    }
    let mut basis: [Vec<Vec<(usize, f64)>>; 4] = Default::default();
    for &(i, j) in at.keys().filter(|(i, j)| *i >= 0 && *j >= 0) {
        for (s, (sx, sy)) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)].into_iter().enumerate() {
            let mut coef = BTreeMap::<usize, f64>::new();
            for (fx, fy) in [(1i64, 1i64), (-1, 1), (1, -1), (-1, -1)] {
                let node = *at.get(&(fx * i, fy * j))
                    .ok_or_else(|| format!("{name}: the grid is not reflection-symmetric at ({i}, {j})"))?;
                let c = (if fx < 0 { sx } else { 1.0 }) * (if fy < 0 { sy } else { 1.0 });
                *coef.entry(node).or_default() += c;
            }
            let norm = coef.values().map(|c| c * c).sum::<f64>().sqrt();
            if norm > 0.0 {
                basis[s].push(coef.into_iter().filter(|(_, c)| *c != 0.0).map(|(k, c)| (k, c / norm)).collect());
            }
        }
    }
    let sizes = [basis[0].len(), basis[1].len(), basis[2].len(), basis[3].len()];
    if sizes.iter().sum::<usize>() != n {
        return Err(format!("{name}: the symmetry-adapted basis has {} vectors for {n} nodes", sizes.iter().sum::<usize>()));
    }
    for (s, &m) in sizes.iter().enumerate() {
        let _ = fs_la::eigen::admit_jacobi_eigh(m)
            .map_err(|e| format!("{name}: parity sector {s} ({m} vectors) is too large for the dense eigensolver: {e}"))?;
    }
    let root: Vec<f64> = body.mass.iter().map(|m| m.sqrt()).collect();
    if root.iter().any(|r| !(*r > 0.0 && r.is_finite())) {
        return Err(format!("{name}: rest modes need positive finite nodal masses"));
    }
    let mut blocks: [Vec<f64>; 4] = Default::default();
    let (mut largest, mut coupling) = (0.0_f64, 0.0_f64);
    let mut x = vec![0.0; n];
    for s in 0..4 {
        let m = sizes[s];
        let mut block = vec![0.0; m * m];
        for q in 0..m {
            for &(i, c) in &basis[s][q] {
                x[i] = c / root[i];
            }
            let mut y = body.rest_stiffness_apply(&x)?;
            for &(i, _) in &basis[s][q] {
                x[i] = 0.0;
            }
            for (v, r) in y.iter_mut().zip(&root) {
                *v /= r;
            }
            for (t, vectors) in basis.iter().enumerate() {
                for (p, u) in vectors.iter().enumerate() {
                    let v: f64 = u.iter().map(|&(i, c)| c * y[i]).sum();
                    largest = largest.max(v.abs());
                    if t == s {
                        block[p * m + q] = v;
                    } else {
                        coupling = coupling.max(v.abs());
                    }
                }
            }
        }
        blocks[s] = block;
    }
    if !(largest > 0.0 && largest.is_finite()) {
        return Err(format!("{name}: the rest stiffness is zero or not finite"));
    }
    let (mut d2, mut n2) = (0.0_f64, 0.0_f64);
    for (s, block) in blocks.iter_mut().enumerate() {
        let m = sizes[s];
        for p in 0..m {
            n2 += block[p * m + p].powi(2);
            for q in (p + 1)..m {
                let (a, b) = (block[p * m + q], block[q * m + p]);
                d2 += 2.0 * (a - b).powi(2);
                n2 += a * a + b * b;
                block[p * m + q] = 0.5 * (a + b);
                block[q * m + p] = 0.5 * (a + b);
            }
        }
    }
    let asymmetry = (d2 / n2).sqrt();
    let coupling = coupling / largest;
    if asymmetry > SYMMETRY_GATE || coupling > SYMMETRY_GATE {
        return Err(format!(
            "{name}: the rest stiffness is not symmetric and reflection-invariant to rounding (asymmetry {asymmetry:.2e}, cross-sector coupling {coupling:.2e}, gate {SYMMETRY_GATE:e})"
        ));
    }
    let solved: Vec<(Vec<f64>, Vec<f64>)> = blocks.iter().zip(sizes).map(|(b, m)| householder_eigh(b, m)).collect();
    let lambda_max = solved.iter().flat_map(|(l, _)| l.iter()).fold(0.0_f64, |a, l| a.max(l.abs()));
    let tol = RIGID_FRACTION * lambda_max;
    let (mut modes, mut rigid) = (Vec::new(), 0);
    for (s, (lam, vectors)) in solved.iter().enumerate() {
        let m = sizes[s];
        for (c, &l) in lam.iter().enumerate() {
            if l < -tol {
                return Err(format!("{name}: the rest state is unstable: eigenvalue {l:e} in parity sector {s}"));
            }
            if l <= tol {
                rigid += 1;
                continue;
            }
            let omega = l.sqrt();
            if omega > omega_limit {
                continue;
            }
            let mut phi = vec![0.0; n];
            for (p, u) in basis[s].iter().enumerate() {
                let yp = vectors[p * m + c];
                for &(i, coef) in u {
                    phi[i] += coef * yp;
                }
            }
            for (v, r) in phi.iter_mut().zip(&root) {
                *v /= r;
            }
            modes.push((omega, phi, s));
        }
    }
    modes.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(RestModes { modes, rigid, omega_max: lambda_max.sqrt(), asymmetry, coupling, sectors: sizes })
}

/// The grid's nodes on a zero-padded square array for the Rayleigh convolution.
pub(crate) struct RayleighGrid {
    size: usize,
    at: Vec<usize>,
    distance: Vec<f64>,
    twiddle: Vec<(f64, f64)>,
}

impl RayleighGrid {
    /// Nodes at integer multiples of `h`; the array is a power of two at least 2 span - 1 wide, so
    /// the circular convolution is the linear one at every node.
    pub(crate) fn new(xy: &[[f64; 2]], h: f64) -> Result<RayleighGrid, String> {
        if xy.is_empty() || !(h > 0.0 && h.is_finite()) {
            return Err("Rayleigh grid: no nodes or no positive finite spacing".into());
        }
        let cells: Vec<(i64, i64)> = xy.iter().map(|p| ((p[0] / h).round() as i64, (p[1] / h).round() as i64)).collect();
        if xy.iter().zip(&cells).any(|(p, c)| (p[0] - c.0 as f64 * h).abs() > 1e-9 * h || (p[1] - c.1 as f64 * h).abs() > 1e-9 * h) {
            return Err("Rayleigh grid: a node is off the square grid".into());
        }
        let (i0, i1) = cells.iter().fold((i64::MAX, i64::MIN), |(a, b), c| (a.min(c.0), b.max(c.0)));
        let (j0, j1) = cells.iter().fold((i64::MAX, i64::MIN), |(a, b), c| (a.min(c.1), b.max(c.1)));
        let span = ((i1 - i0).max(j1 - j0) + 1) as usize;
        let size = (2 * span - 1).next_power_of_two().max(2);
        let at = cells.iter().map(|&(i, j)| (i - i0) as usize + size * (j - j0) as usize).collect();
        let offset = |a: usize| if a <= size / 2 { a as f64 } else { a as f64 - size as f64 };
        let distance = (0..size * size).map(|p| h * (offset(p % size).powi(2) + offset(p / size).powi(2)).sqrt()).collect();
        let twiddle = (0..size / 2).map(|k| {
            let t = -2.0 * PI * k as f64 / size as f64;
            (t.cos(), t.sin())
        }).collect();
        Ok(RayleighGrid { size, at, distance, twiddle })
    }

    /// (G * phi)(x) = sum_y sin(k |x - y|) / |x - y| phi(y) at every node x (the y = x term k phi(x)).
    pub(crate) fn convolve(&self, phi: &[f64], k: f64) -> Vec<f64> {
        assert_eq!(phi.len(), self.at.len(), "one value per node");
        let len = self.size * self.size;
        let mut gr: Vec<f64> = self.distance.iter().map(|&d| if d > 0.0 { (k * d).sin() / d } else { k }).collect();
        let mut gi = vec![0.0; len];
        let (mut pr, mut pi) = (vec![0.0; len], vec![0.0; len]);
        for (&a, &v) in self.at.iter().zip(phi) {
            pr[a] = v;
        }
        self.fft2(&mut gr, &mut gi, false);
        self.fft2(&mut pr, &mut pi, false);
        for p in 0..len {
            let (a, b) = (gr[p] * pr[p] - gi[p] * pi[p], gr[p] * pi[p] + gi[p] * pr[p]);
            pr[p] = a;
            pi[p] = b;
        }
        self.fft2(&mut pr, &mut pi, true);
        let scale = 1.0 / len as f64;
        self.at.iter().map(|&a| pr[a] * scale).collect()
    }

    /// In-place radix-2 transform of one row or column (unnormalised; `inverse` conjugates).
    fn fft(&self, re: &mut [f64], im: &mut [f64], inverse: bool) {
        let n = re.len();
        let mut j = 0usize;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j |= bit;
            if i < j {
                re.swap(i, j);
                im.swap(i, j);
            }
        }
        let mut len = 2;
        while len <= n {
            let (half, step) = (len / 2, n / len);
            for start in (0..n).step_by(len) {
                for k in 0..half {
                    let (c, s) = self.twiddle[k * step];
                    let s = if inverse { -s } else { s };
                    let (a, b) = (start + k, start + k + half);
                    let (xr, xi) = (re[b] * c - im[b] * s, re[b] * s + im[b] * c);
                    re[b] = re[a] - xr;
                    im[b] = im[a] - xi;
                    re[a] += xr;
                    im[a] += xi;
                }
            }
            len <<= 1;
        }
    }

    fn fft2(&self, re: &mut [f64], im: &mut [f64], inverse: bool) {
        let n = self.size;
        for row in 0..n {
            let span = row * n..(row + 1) * n;
            self.fft(&mut re[span.clone()], &mut im[span], inverse);
        }
        let (mut cr, mut ci) = (vec![0.0; n], vec![0.0; n]);
        for col in 0..n {
            for r in 0..n {
                cr[r] = re[r * n + col];
                ci[r] = im[r * n + col];
            }
            self.fft(&mut cr, &mut ci, inverse);
            for r in 0..n {
                re[r * n + col] = cr[r];
                im[r * n + col] = ci[r];
            }
        }
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n == 0 { f64::NAN } else if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

fn t60_text(t: f64) -> String {
    if t.is_finite() { format!("{t:.3} s") } else { "none".into() }
}

/// The eigen-decomposition of a cluster's resistance block (row-major `s` x `s`, rates 1/s), after
/// symmetrising it. Refused unless the block was symmetric to `BLOCK_GATE` of its largest entry, its
/// largest eigenvalue is positive and finite and its smallest is at least -`BLOCK_GATE` of the largest.
/// Returns the eigenvalues ascending, a rounding negative set to zero; the eigenvectors as columns
/// (row-major); the asymmetry; the smallest eigenvalue over the largest, before zeroing; and how many
/// eigenvalues were set to zero.
pub(crate) fn gated_block(mut block: Vec<f64>, s: usize) -> Result<(Vec<f64>, Vec<f64>, f64, f64, usize), String> {
    if s == 0 || block.len() != s * s {
        return Err(format!("a resistance block of {} entries is not {s} x {s}", block.len()));
    }
    let big = block.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if !(big > 0.0 && big.is_finite()) {
        return Err(format!("a resistance block's largest entry is {big:e}"));
    }
    let mut skew = 0.0_f64;
    for p in 0..s {
        for q in (p + 1)..s {
            let (x, y) = (block[p * s + q], block[q * s + p]);
            skew = skew.max((x - y).abs());
            block[p * s + q] = 0.5 * (x + y);
            block[q * s + p] = 0.5 * (x + y);
        }
    }
    let asymmetry = skew / big;
    let (lam, vectors) = householder_eigh(&block, s);
    let top = lam[s - 1];
    let smallest = lam[0] / top;
    if !(asymmetry <= BLOCK_GATE && top > 0.0 && top.is_finite() && smallest >= -BLOCK_GATE) {
        return Err(format!(
            "a resistance block is not symmetric and positive semidefinite to {BLOCK_GATE:e} (asymmetry {asymmetry:.2e}, largest eigenvalue {top:e}, smallest over largest {smallest:.2e})"
        ));
    }
    let zeroed = lam.iter().filter(|l| !(**l > 0.0)).count();
    let lam = lam.into_iter().map(|l| if l > 0.0 { l } else { 0.0 }).collect();
    Ok((lam, vectors, asymmetry, smallest, zeroed))
}

/// Loaded modes of one parity sector that share their resistance block (`radiation_load`).
pub(crate) struct Cluster {
    /// Indices into the loaded modes, ascending.
    pub(crate) members: Vec<usize>,
    /// The mean angular frequency the block is taken at [rad/s].
    pub(crate) omega: f64,
    /// The block's eigenvalues 2 sigma [1/s], ascending, a rounding negative set to zero.
    pub(crate) two_sigma: Vec<f64>,
    /// The block's eigenvectors, M-orthonormal combinations of its modes, in the order of `two_sigma`.
    pub(crate) shapes: Vec<Vec<f64>>,
    /// The block's asymmetry before symmetrising, and its smallest eigenvalue over its largest.
    pub(crate) asymmetry: f64,
    pub(crate) smallest: f64,
    /// The largest overlap between two of its modes.
    pub(crate) overlap: f64,
}

/// The radiation load on a shell's loaded modes (`radiation_load`).
pub(crate) struct RadiationLoad<'a> {
    /// Per loaded mode: its one-face resistance R_mm at its own frequency.
    pub(crate) own: Vec<f64>,
    /// Per loaded mode: its own loss alone, 2 |H|^2 R_mm (R_mm rounded up to zero if below).
    pub(crate) two_sigma: Vec<f64>,
    /// Per loaded mode: its largest overlap with another loaded mode of its sector.
    pub(crate) overlap: Vec<f64>,
    /// Per loaded mode: sum_n |R_mn| / R_mm at its own frequency over the modes of its sector outside
    /// its cluster (the cross resistance the load drops), infinite where R_mm is not positive.
    pub(crate) share: Vec<f64>,
    /// Per loaded mode: its cluster's index in `clusters`, if it is in one.
    pub(crate) cluster: Vec<Option<usize>>,
    pub(crate) clusters: Vec<Cluster>,
    /// The largest overlap of any pair of one sector, and of any pair in no common cluster:
    /// (overlap, mode, mode).
    pub(crate) largest: Option<(f64, usize, usize)>,
    pub(crate) largest_apart: Option<(f64, usize, usize)>,
    /// The load as kept, in mode order: each mode alone with its own loss, each cluster (at its first
    /// mode) with its eigenvectors; only shapes with a positive rate. M-orthonormal.
    pub(crate) shapes: Vec<Cow<'a, [f64]>>,
    pub(crate) rates: Vec<f64>,
    /// Modes whose own resistance rounded to zero or below; cluster eigenvalues set to zero.
    pub(crate) clipped: usize,
    pub(crate) zeroed: usize,
}

/// The radiation load on `modes` (omega rad/s, phi with phi' M phi = 1, parity sector), on `grid` at
/// spacing `h`, in a gas of density `rho` and sound speed `c`, with the dipole factor's time `tau` (module
/// docs): every mode's resistance at its own frequency to each mode of its sector, the pairs' overlaps,
/// the clusters that the overlaps above `threshold` join (+inf leaves every mode alone), and each
/// cluster's block at its mean frequency (`gated_block`). Refuses a threshold that is not positive, a
/// resistance that is not finite and every refusal of `gated_block`.
pub(crate) fn radiation_load<'a>(modes: &[&'a (f64, Vec<f64>, usize)], grid: &RayleighGrid, h: f64, rho: f64, c: f64, tau: f64, threshold: f64) -> Result<RadiationLoad<'a>, String> {
    if !(threshold > 0.0) {
        return Err(format!("the cluster overlap must be positive, +inf for none (got {threshold})"));
    }
    let n = modes.len();
    let h4 = h.powi(4);
    let hz = |omega: f64| omega / (2.0 * PI);
    let dipole = |omega: f64| {
        let wt = omega * tau;
        wt * wt / (1.0 + wt * wt)
    };
    // at[a][b]: R between modes a and b at mode a's frequency, (omega_a rho0 / 2 pi) h^4 phi_b' G phi_a
    let mut at = vec![vec![0.0; n]; n];
    for (a, mode) in modes.iter().enumerate() {
        let (omega, phi, sector) = (mode.0, &mode.1, mode.2);
        let u = grid.convolve(phi, omega / c);
        let factor = omega * rho / (2.0 * PI) * h4;
        for (b, other) in modes.iter().enumerate() {
            if other.2 == sector {
                at[a][b] = factor * dot(&other.1, &u);
            }
        }
        if at[a].iter().any(|r| !r.is_finite()) {
            return Err(format!("a resistance at {:.1} Hz is not finite", hz(omega)));
        }
    }
    let own: Vec<f64> = (0..n).map(|a| at[a][a]).collect();
    let clipped = own.iter().filter(|r| !(**r > 0.0)).count();
    let two_sigma: Vec<f64> = (0..n).map(|a| {
        let r1 = if own[a] > 0.0 { own[a] } else { 0.0 };
        2.0 * dipole(modes[a].0) * r1
    }).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    // union by the smaller index, so each cluster's root is its first mode
    let mut parent: Vec<usize> = (0..n).collect();
    let mut overlap = vec![0.0_f64; n];
    let mut pairs = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            if modes[a].2 != modes[b].2 {
                continue;
            }
            let coupling = (2.0 * dipole(modes[a].0) * at[a][b].abs()).max(2.0 * dipole(modes[b].0) * at[b][a].abs());
            let gap = (modes[a].0 - modes[b].0).abs();
            let o = if coupling == 0.0 { 0.0 } else if gap > 0.0 { coupling / gap } else { f64::INFINITY };
            overlap[a] = overlap[a].max(o);
            overlap[b] = overlap[b].max(o);
            pairs.push((o, a, b));
            if o > threshold {
                let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
    }
    let roots: Vec<usize> = (0..n).map(|i| root(&mut parent, i)).collect();
    let mut size = vec![0usize; n];
    for &r in &roots {
        size[r] += 1;
    }
    let mut cluster = vec![None; n];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut slot = BTreeMap::new();
    for a in 0..n {
        if size[roots[a]] > 1 {
            let next = groups.len();
            let g = *slot.entry(roots[a]).or_insert(next);
            if g == next {
                groups.push(Vec::new());
            }
            groups[g].push(a);
            cluster[a] = Some(g);
        }
    }
    let first = |keep: &dyn Fn(&(f64, usize, usize)) -> bool| {
        pairs.iter().filter(|p| keep(*p)).fold(None, |best: Option<(f64, usize, usize)>, p| match best {
            Some(q) if q.0 >= p.0 => Some(q),
            _ => Some(*p),
        })
    };
    let largest = first(&|_: &(f64, usize, usize)| true);
    let largest_apart = first(&|p: &(f64, usize, usize)| roots[p.1] != roots[p.2]);
    let share: Vec<f64> = (0..n).map(|a| {
        let cross: f64 = (0..n)
            .filter(|&b| b != a && modes[b].2 == modes[a].2 && (cluster[a].is_none() || cluster[b] != cluster[a]))
            .map(|b| at[a][b].abs())
            .sum();
        if own[a] > 0.0 { cross / own[a] } else { f64::INFINITY }
    }).collect();
    let mut clusters = Vec::with_capacity(groups.len());
    for (g, members) in groups.into_iter().enumerate() {
        let s = members.len();
        let omega = members.iter().map(|&i| modes[i].0).sum::<f64>() / s as f64;
        let factor = omega * rho / (2.0 * PI) * h4;
        let scale = 2.0 * dipole(omega);
        let u: Vec<Vec<f64>> = members.iter().map(|&j| grid.convolve(&modes[j].1, omega / c)).collect();
        let mut block = vec![0.0; s * s];
        for (p, &i) in members.iter().enumerate() {
            for (q, uq) in u.iter().enumerate() {
                block[p * s + q] = scale * (factor * dot(&modes[i].1, uq));
            }
        }
        let (lam, vectors, asymmetry, smallest, zeroed) = gated_block(block, s)
            .map_err(|e| format!("the cluster of {s} modes at {:.1} Hz (parity sector {}): {e}", hz(omega), modes[members[0]].2))?;
        let nodes = modes[members[0]].1.len();
        let shapes: Vec<Vec<f64>> = (0..s).map(|col| {
            let mut phi = vec![0.0; nodes];
            for (p, &i) in members.iter().enumerate() {
                let w = vectors[p * s + col];
                for (v, x) in phi.iter_mut().zip(&modes[i].1) {
                    *v += w * x;
                }
            }
            phi
        }).collect();
        let inner = pairs.iter().filter(|p| cluster[p.1] == Some(g) && cluster[p.2] == Some(g)).fold(0.0_f64, |m, p| m.max(p.0));
        clusters.push((Cluster { members, omega, two_sigma: lam, shapes, asymmetry, smallest, overlap: inner }, zeroed));
    }
    let zeroed: usize = clusters.iter().map(|c| c.1).sum();
    let clusters: Vec<Cluster> = clusters.into_iter().map(|c| c.0).collect();
    let (mut shapes, mut rates) = (Vec::new(), Vec::new());
    for a in 0..n {
        match cluster[a] {
            None => {
                if two_sigma[a] > 0.0 {
                    shapes.push(Cow::Borrowed(&modes[a].1[..]));
                    rates.push(two_sigma[a]);
                }
            }
            Some(g) if clusters[g].members[0] == a => {
                for (phi, &r) in clusters[g].shapes.iter().zip(&clusters[g].two_sigma) {
                    if r > 0.0 {
                        shapes.push(Cow::Owned(phi.clone()));
                        rates.push(r);
                    }
                }
            }
            Some(_) => {}
        }
    }
    Ok(RadiationLoad { own, two_sigma, overlap, share, cluster, clusters, largest, largest_apart, shapes, rates, clipped, zeroed })
}

impl Body {
    /// Attach the radiation loss of this free shell's own rest modes (module docs): every mode up to
    /// `loaded_hz` takes 2 sigma = 2 |H|^2 R, R its one-face baffled Rayleigh resistance in `gas`,
    /// |H|^2 the microphone's first-order dipole factor for a sheet of radius `dipole_radius`, and the
    /// modes of one sector whose overlap exceeds `cluster_overlap` share their full resistance block
    /// (`kit` passes `CLUSTER_OVERLAP`; +inf leaves every mode alone); the modes from there to
    /// `listed_hz` are counted and take none. The loss adds to every loss the shell already carries.
    /// Refused on a tensioned head (heads take `attach_head_radiation`), a body that does not radiate,
    /// a second attach, a body with shape damping or a head radiation load, invalid arguments, no
    /// loaded mode, and every refusal of the rest-mode solve and of `radiation_load`. Returns the
    /// receipt lines: a summary, the operator checks and overlaps, octave medians, one line per
    /// cluster, then one line per loaded mode.
    pub fn attach_shell_radiation(&mut self, gas: &GasState, dipole_radius: f64, loaded_hz: f64, listed_hz: f64, cluster_overlap: f64) -> Result<Vec<String>, String> {
        let name = self.name.clone();
        if self.tension > 0.0 || self.n < 9 {
            return Err(format!("{name}: shell radiation applies to free untensioned shells; a head takes its radiation table"));
        }
        if !self.radiates {
            return Err(format!("{name}: shell radiation needs a body the microphones hear"));
        }
        if self.shape_damping.is_some() || self.head_radiation.is_some() {
            return Err(format!("{name}: shell radiation attaches once, and not beside shape damping or a head radiation load"));
        }
        let (rho, c) = (gas.density, gas.sound_speed);
        if !(rho > 0.0 && rho.is_finite() && c > 0.0 && c.is_finite() && dipole_radius > 0.0 && dipole_radius.is_finite()) {
            return Err(format!("{name}: shell radiation needs a positive finite gas density, sound speed and dipole radius"));
        }
        if !(loaded_hz > 0.0 && loaded_hz.is_finite() && listed_hz >= loaded_hz && listed_hz.is_finite()) {
            return Err(format!("{name}: shell radiation needs 0 < loaded_hz <= listed_hz, finite (got {loaded_hz} and {listed_hz})"));
        }
        if !(cluster_overlap > 0.0) {
            return Err(format!("{name}: shell radiation needs a positive cluster overlap, +inf for none (got {cluster_overlap})"));
        }
        let rest = rest_modes(self, 2.0 * PI * listed_hz)?;
        let grid = RayleighGrid::new(&self.xy, self.h)?;
        let tau = 4.0 * dipole_radius / (3.0 * PI * c);
        let ih2 = 1.0 / (self.h * self.h);
        let loaded: Vec<&(f64, Vec<f64>, usize)> = rest.modes.iter().filter(|m| m.0 <= 2.0 * PI * loaded_hz).collect();
        let unloaded = rest.modes.len() - loaded.len();
        if loaded.is_empty() {
            return Err(format!("{name}: no rest mode lies at or below {loaded_hz:.1} Hz"));
        }
        let load = radiation_load(&loaded, &grid, self.h, rho, c, tau, cluster_overlap)
            .map_err(|e| format!("{name}: shell radiation: {e}"))?;
        // sigma of this body's other losses on a shape with phi' M phi = 1
        let other_of = |phi: &[f64]| -> f64 {
            let mut other: f64 = phi.iter().zip(&self.mass).zip(&self.s0).map(|((p, m), s)| m * s * p * p).sum();
            if self.s1 > 0.0 {
                other += self.s1 * ih2 * self.lap.mass_flux_quadratic(phi, &self.mass);
            }
            if let Some(v) = &self.viscosity {
                other += 0.5 * v.c.quad(phi);
            }
            if let Some(m) = &self.modal {
                other += 0.5 * m.c.quad(phi);
            }
            other
        };
        // (omega, sector, 2 sigma alone, eta without the dipole factor, sigma of the other losses, dropped share)
        let rows: Vec<(f64, usize, f64, f64, f64, f64)> = loaded.iter().enumerate()
            .map(|(i, m)| (m.0, m.2, load.two_sigma[i], 2.0 * load.own[i].max(0.0) / m.0, other_of(&m.1[..]), load.share[i]))
            .collect();
        let others: Vec<Vec<f64>> = load.clusters.iter().map(|k| k.shapes.iter().map(|p| other_of(&p[..])).collect()).collect();
        let mphi: Vec<Vec<f64>> = load.shapes.iter().map(|p| p.iter().zip(&self.mass).map(|(p, m)| p * m).collect()).collect();
        let two_sigma = load.rates.clone();
        // Gershgorin over the kept shapes' Gram matrix scaled by sqrt(2 sigma_i 2 sigma_j), as for heads
        let mut bound = 0.0_f64;
        for i in 0..mphi.len() {
            let row: f64 = (0..mphi.len()).map(|j| (two_sigma[i] * two_sigma[j]).sqrt() * dot(&mphi[i], &load.shapes[j]).abs()).sum();
            bound = bound.max(row);
        }
        let largest = two_sigma.iter().cloned().fold(0.0_f64, f64::max);
        let hz = |omega: f64| omega / (2.0 * PI);
        let t60 = |sigma: f64| if sigma > 0.0 { 6.907755 / sigma } else { f64::INFINITY };
        let (worst, worst_hz) = rows.iter().fold((0.0_f64, 0.0), |acc, r| if r.5 > acc.0 { (r.5, hz(r.0)) } else { acc });
        let overlap_text = |o: f64| if o.is_finite() { format!("{o:.4}") } else { "unbounded (equal frequencies)".into() };
        let pair_text = |p: Option<(f64, usize, usize)>| match p {
            Some((o, a, b)) => format!("{} ({:.1} and {:.1} Hz, parity sector {})", overlap_text(o), hz(loaded[a].0), hz(loaded[b].0), loaded[a].2),
            None => "none (no two loaded modes share a sector)".into(),
        };
        let clusters_text = if load.clusters.is_empty() {
            format!("no pair overlaps above {cluster_overlap}, so every mode keeps only its own resistance")
        } else {
            format!(
                "the {} modes joined by overlaps above {cluster_overlap} form {} clusters, each of which takes its full resistance block 2 |H|^2 R at its mean frequency (MB 14957, a correction), symmetric and positive semidefinite to {BLOCK_GATE:e} (largest asymmetry {:.2e}, smallest eigenvalue over largest {:.2e}, {} eigenvalues rounded to zero); the largest overlap of two modes in no common cluster is {}",
                load.clusters.iter().map(|k| k.members.len()).sum::<usize>(), load.clusters.len(),
                load.clusters.iter().map(|k| k.asymmetry).fold(0.0_f64, f64::max),
                load.clusters.iter().map(|k| k.smallest).fold(f64::INFINITY, f64::min),
                load.zeroed, pair_text(load.largest_apart)
            )
        };
        let mut lines = vec![
            format!(
                "shell radiation (MB 14574 (d), a correction): {} rest modes of this shell up to {loaded_hz:.0} Hz each take the baffled Rayleigh resistance of its own shape on the shell's projected plane, from both faces, times the microphone's own first-order dipole factor (omega tau)^2/(1 + (omega tau)^2), tau = 4a/(3 pi c) = {:.4} ms (a {dipole_radius} m, c {c:.2} m/s, rho0 {rho:.4} kg/m^3); an approximation: neither edge diffraction nor the curvature's effect on radiation is resolved; {unloaded} modes from {loaded_hz:.0} Hz to {listed_hz:.0} Hz take none, so their decay times are upper bounds; {} rigid motions take none",
                rows.len(), tau * 1e3, rest.rigid
            ),
            format!(
                "shell radiation modes: K about rest from this body's own operator (bending, supports, curvature Hessian), asymmetry {:.2e} and cross-sector coupling {:.2e} (gate {SYMMETRY_GATE:e}), parity sectors {:?}, highest grid frequency {:.0} Hz; the overlap of two loaded modes of one sector, |2 sigma_mn| / |omega_m - omega_n| with 2 sigma_mn = 2 |H|^2 R_mn at either mode's frequency (the larger), is at most {}; {clusters_text}; outside its cluster, a mode's dropped cross resistance to the other loaded modes of its sector, at its own frequency, is at most {worst:.3} of its own ({worst_hz:.1} Hz), median {:.3}; {} modes rounded to no resistance; rate bound {bound:.3} 1/s against 2 max sigma {largest:.3} 1/s",
                rest.asymmetry, rest.coupling, rest.sectors, hz(rest.omega_max), pair_text(load.largest), median(rows.iter().map(|r| r.5).collect()), load.clipped
            ),
        ];
        let mut bands: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (i, r) in rows.iter().enumerate() {
            bands.entry((hz(r.0) / 1000.0).log2().floor() as i32).or_default().push(i);
        }
        let octaves: Vec<String> = bands.iter().map(|(&b, v)| {
            let lo = 1000.0 * 2f64.powi(b);
            let of = |pick: &dyn Fn(usize) -> f64| t60_text(median(v.iter().map(|&i| pick(i)).collect()));
            format!(
                "{:.0}-{:.0} Hz ({} modes): {} / {} / {}",
                lo, 2.0 * lo, v.len(),
                of(&|i| t60(rows[i].2 / 2.0)),
                of(&|i| t60(rows[i].4)),
                of(&|i| t60(rows[i].4 + rows[i].2 / 2.0)),
            )
        }).collect();
        lines.push(format!("shell radiation by octave, median T60 of each mode's own resistance (clusters below): radiation / this body's other losses / both: {}", octaves.join("; ")));
        let list = |v: Vec<String>| v.join(", ");
        for (g, k) in load.clusters.iter().enumerate() {
            lines.push(format!(
                "shell radiation cluster {g}: modes {} at {} Hz, parity sector {}, overlap up to {}; the block at their mean {:.1} Hz: T60 radiation {}, both {} (each mode alone: radiation {}, both {})",
                list(k.members.iter().map(|i| i.to_string()).collect()),
                list(k.members.iter().map(|&i| format!("{:.1}", hz(rows[i].0))).collect()),
                rows[k.members[0]].1, overlap_text(k.overlap), hz(k.omega),
                list(k.two_sigma.iter().map(|r| t60_text(t60(r / 2.0))).collect()),
                list(k.two_sigma.iter().zip(&others[g]).map(|(r, o)| t60_text(t60(o + r / 2.0))).collect()),
                list(k.members.iter().map(|&i| t60_text(t60(rows[i].2 / 2.0))).collect()),
                list(k.members.iter().map(|&i| t60_text(t60(rows[i].4 + rows[i].2 / 2.0))).collect()),
            ));
        }
        for (i, r) in rows.iter().enumerate() {
            lines.push(format!(
                "shell radiation mode {i}: {:.1} Hz, parity sector {}: eta {:.3e} (both faces baffled, without the dipole factor, {:.3e}); T60 radiation {}, other losses {}, both {}; dropped cross share {:.3}; overlap {}, cluster {}",
                hz(r.0), r.1, r.2 / r.0, r.3, t60_text(t60(r.2 / 2.0)), t60_text(t60(r.4)), t60_text(t60(r.4 + r.2 / 2.0)), r.5,
                overlap_text(load.overlap[i]), load.cluster[i].map_or("none".into(), |g| g.to_string())
            ));
        }
        self.shape_damping = Some(ShapeDamping { mphi, two_sigma, rate_bound: bound, dissipated_j: 0.0 });
        Ok(lines)
    }
}
