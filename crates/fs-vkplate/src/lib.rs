//! # fs-vkplate: struck thin plates and membranes in the time domain
//!
//! One explicit finite-difference solver for every vibrating surface of a drum kit.
//! A body is a thin plate on a square grid (spacing `h`) with
//!
//! `rho t w_tt = T Lap(w) - div div (D curv) + vK(w, F) - 2 rho t (s0 w_t - s1 Lap(w_t)) + p + f_contact`
//!
//! - **membranes** (drum heads): clamped rim, in-plane pretension `T`, optional bending,
//!   Berger tension modulation (the mean stretch raises the tension: pitch glide);
//! - **plates and shallow shells** (cymbals): free edges, bending with variable thickness,
//!   full von Karman coupling through the Airy stress function `F`, including a rest shape
//!   `z0(r)` (shallow shell): `Lap^2 F = -E t [l(w + z0) - l(z0)]`, `l(u) = u_xx u_yy - u_xy^2`.
//!
//! The discrete operators come from discrete ENERGIES (bending energy summed over 3x3
//! curvature stencils, membrane energy over grid edges, Airy energy over the grid), so the
//! linear part is symmetric and the nonlinear force is the exact gradient of the discrete
//! von Karman potential. Free edges are the natural boundary condition of that energy.
//!
//! Time stepping is leapfrog. Structural loss is either the legacy `s0 + s1 beta^2`
//! prescription from two (frequency, T60) points or an explicitly selected, source-bound
//! plane-stress Kelvin--Voigt bending tensor. Its backward-velocity force requires the
//! damping CFL and its physical loss uses centered velocity; see `Body::discrete_energy`.
//! Contacts (stick, felt beater, snare wires, hi-hat plates) are
//! solved IMPLICITLY per contact as monotone scalar equations, so stiff Hertz contacts are
//! stable at audio rates. The enclosed air of a drum couples its two heads (lumped adiabatic
//! cavity); each head's volume velocity sees a small-ka monopole radiation resistance. Sound
//! is the baffled Rayleigh integral of surface acceleration at each microphone.
//!
//! Honest scope: square-grid (staircase) boundaries; moderate rotations (von Karman); preset
//! loss parameters are declared estimates, while sourced bending loss retains its admitted
//! fixed-state/band limits; lumped cavity (no cavity modes); baffled
//! radiation (no diffraction around a cymbal's edge).

use std::collections::BTreeMap;
use std::f64::consts::PI;

pub use fs_material::gas::GasState;
pub use fs_material::plane_stress_viscosity::PlaneStressViscosity;

/// Immutable source-bound elastic data and physical striker geometry for kit builders.
pub mod materials;

/// Baked exterior modal loading of drum heads' non-volume motion.
pub mod radiation;

/// Default environmental datum; all coefficients come from the shared gas law.
/// Custom atmospheres enter through `Instrument::with_gas`.
pub fn standard_gas() -> GasState {
    GasState::try_new_moist_air(293.15, 101_325.0, 0.0)
        .expect("the standard dry-air state is in the material model's domain")
}

#[cfg(test)]
mod material_gas_tests;

#[cfg(test)]
mod material_tests;

#[cfg(test)]
mod radiation_tests;

// ----------------------------------------------------------------------------- sparse
#[derive(Clone, Debug, Default)]
pub struct Csr {
    pub n: usize,
    pub ptr: Vec<usize>,
    pub col: Vec<usize>,
    pub val: Vec<f64>,
    /// The grid interior as stencils: rows that `mul_add` computes without the column indices.
    /// Built from `ptr`/`col`/`val` once, in `from_rows` (nothing edits a matrix afterwards).
    runs: Vec<Run>,
}

/// Consecutive rows `r0..r0 + len` whose entries sit at the same offsets from the diagonal (one
/// grid row of a stencil's interior). Offsets are ascending, so a row sums in its CSR order.
#[derive(Clone, Debug)]
struct Run {
    r0: usize,
    len: usize,
    /// Column minus row, per entry.
    off: Vec<isize>,
    /// One coefficient per entry when every row has the same ones; otherwise every row's, entry-major
    /// (`val[k * len + t]` belongs to row `r0 + t`).
    uniform: bool,
    val: Vec<f64>,
}

/// Shortest run worth a stencil loop.
const MIN_RUN: usize = 4;

impl Csr {
    fn from_rows(rows: &[BTreeMap<usize, f64>]) -> Csr {
        let mut c = Csr { n: rows.len(), ptr: vec![0], col: vec![], val: vec![], runs: vec![] };
        for r in rows {
            for (&j, &v) in r {
                if v != 0.0 {
                    c.col.push(j);
                    c.val.push(v);
                }
            }
            c.ptr.push(c.col.len());
        }
        c.runs = c.find_runs();
        c
    }

    /// Maximal runs of rows sharing one stencil shape with 3, 5 or 17 entries (a string, the
    /// Laplacian, the bending operator).
    fn find_runs(&self) -> Vec<Run> {
        let (ptr, col) = (&self.ptr, &self.col);
        let off = |i: usize| col[ptr[i]..ptr[i + 1]].iter().map(move |&j| j as isize - i as isize);
        let mut runs = vec![];
        let mut i = 0;
        while i < self.n {
            let k = ptr[i + 1] - ptr[i];
            let mut e = i + 1;
            if matches!(k, 3 | 5 | 17) {
                while e < self.n && ptr[e + 1] - ptr[e] == k && off(e).eq(off(i)) {
                    e += 1;
                }
            }
            let len = e - i;
            if len >= MIN_RUN && matches!(k, 3 | 5 | 17) {
                let row = |t: usize| &self.val[ptr[i + t]..ptr[i + t + 1]];
                let uniform = (1..len).all(|t| row(t).iter().zip(row(0)).all(|(a, b)| a.to_bits() == b.to_bits()));
                let val = if uniform { row(0).to_vec() } else { (0..k).flat_map(|q| (0..len).map(move |t| row(t)[q])).collect() };
                runs.push(Run { r0: i, len, off: off(i).collect(), uniform, val });
            }
            i = e;
        }
        runs
    }

    /// y += s * A x
    #[inline]
    pub fn mul_add(&self, x: &[f64], s: f64, y: &mut [f64]) {
        if self.n == 0 {
            return; // no rows (and a `Csr::default()` has no `ptr[0]` for `rows_add` to slice)
        }
        let y = &mut y[..self.n];
        let mut i = 0;
        for r in &self.runs {
            self.rows_add(i, r.r0, x, s, y);
            match r.off.len() {
                3 => r.mul_add::<3>(x, s, y),
                5 => r.mul_add::<5>(x, s, y),
                _ => r.mul_add::<17>(x, s, y),
            }
            i = r.r0 + r.len;
        }
        self.rows_add(i, self.n, x, s, y);
    }

    /// y += s * A x on rows `a..b`, from the CSR arrays.
    #[inline]
    fn rows_add(&self, a: usize, b: usize, x: &[f64], s: f64, y: &mut [f64]) {
        // Interleave independent rows; each accumulator still follows its own ascending columns.
        let mut row = a;
        let mut pairs = y[a..b].chunks_exact_mut(2);
        for pair in pairs.by_ref() {
            let (start, middle, end) = (self.ptr[row], self.ptr[row + 1], self.ptr[row + 2]);
            let (cols0, vals0) = (&self.col[start..middle], &self.val[start..middle]);
            let (cols1, vals1) = (&self.col[middle..end], &self.val[middle..end]);
            let common = cols0.len().min(cols1.len());
            let (mut acc0, mut acc1) = (0.0, 0.0);
            for ((&j0, &v0), (&j1, &v1)) in cols0[..common]
                .iter()
                .zip(&vals0[..common])
                .zip(cols1[..common].iter().zip(&vals1[..common]))
            {
                acc0 += v0 * x[j0];
                acc1 += v1 * x[j1];
            }
            for (&j, &v) in cols0[common..].iter().zip(&vals0[common..]) {
                acc0 += v * x[j];
            }
            for (&j, &v) in cols1[common..].iter().zip(&vals1[common..]) {
                acc1 += v * x[j];
            }
            pair[0] += s * acc0;
            pair[1] += s * acc1;
            row += 2;
        }
        for yi in pairs.into_remainder() {
            let (start, end) = (self.ptr[row], self.ptr[row + 1]);
            let (cols, vals) = (&self.col[start..end], &self.val[start..end]);
            let mut acc = 0.0;
            for (&j, &v) in cols.iter().zip(vals) {
                acc += v * x[j];
            }
            *yi += s * acc;
        }
    }
    pub fn quad(&self, x: &[f64]) -> f64 {
        let mut q = 0.0;
        for i in 0..self.n {
            let mut acc = 0.0;
            for p in self.ptr[i]..self.ptr[i + 1] {
                acc += self.val[p] * x[self.col[p]];
            }
            q += x[i] * acc;
        }
        q
    }
}

impl Run {
    /// y += s * A x on this run's rows: each row's products are summed in the CSR row's order
    /// (ascending column), from 0, then scaled and added, exactly as `Csr::rows_add` does. Rows are
    /// independent, so the compiler may compute several at once (SIMD across rows) without
    /// changing any row's sum.
    #[inline(always)]
    fn mul_add<const K: usize>(&self, x: &[f64], s: f64, y: &mut [f64]) {
        let len = self.len;
        let y = &mut y[self.r0..self.r0 + len];
        let xs: [&[f64]; K] = std::array::from_fn(|k| {
            let a = (self.r0 as isize + self.off[k]) as usize;
            &x[a..a + len]
        });
        if self.uniform {
            let v: [f64; K] = std::array::from_fn(|k| self.val[k]);
            for t in 0..len {
                let mut acc = 0.0;
                for k in 0..K {
                    acc += v[k] * xs[k][t];
                }
                y[t] += s * acc;
            }
        } else {
            let v: [&[f64]; K] = std::array::from_fn(|k| &self.val[k * len..k * len + len]);
            for t in 0..len {
                let mut acc = 0.0;
                for k in 0..K {
                    acc += v[k][t] * xs[k][t];
                }
                y[t] += s * acc;
            }
        }
    }
}

fn add(rows: &mut [BTreeMap<usize, f64>], i: usize, j: usize, v: f64) {
    *rows[i].entry(j).or_insert(0.0) += v;
}

// ------------------------------------------------------------------ banded Cholesky
/// Symmetric positive definite band matrix, factored once (the Airy operator is constant).
#[derive(Clone, Debug)]
pub struct Band {
    n: usize,
    b: usize,
    l: Vec<f64>, // L(i, j) at l[i*(b+1) + (j + b - i)], j in [i-b, i]
}

/// Dot product with eight independent accumulators (a fixed, deterministic summation order that
/// keeps several SIMD lanes and FMA pipes busy).
#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let mut s = [0.0f64; 8];
    let (ca, cb) = (a.chunks_exact(8), b.chunks_exact(8));
    let (ra, rb) = (ca.remainder(), cb.remainder());
    for (x, y) in ca.zip(cb) {
        for k in 0..8 {
            s[k] += x[k] * y[k];
        }
    }
    let mut t = ((s[0] + s[4]) + (s[1] + s[5])) + ((s[2] + s[6]) + (s[3] + s[7]));
    for (x, y) in ra.iter().zip(rb) {
        t += x * y;
    }
    t
}

impl Band {
    pub fn factor(a: &[BTreeMap<usize, f64>]) -> Result<Band, String> {
        let n = a.len();
        let mut b = 0;
        for (i, r) in a.iter().enumerate() {
            for &j in r.keys() {
                b = b.max(i.abs_diff(j));
            }
        }
        let w = b + 1;
        let mut l = vec![0.0; n * w];
        for (i, r) in a.iter().enumerate() {
            for (&j, &v) in r.range(i.saturating_sub(b)..=i) {
                l[i * w + (j + b - i)] = v;
            }
        }
        for i in 0..n {
            for j in i.saturating_sub(b)..=i {
                let mut s = l[i * w + (j + b - i)];
                let k0 = i.saturating_sub(b).max(j.saturating_sub(b));
                for k in k0..j {
                    s -= l[i * w + (k + b - i)] * l[j * w + (k + b - j)];
                }
                if i == j {
                    if s <= 0.0 {
                        return Err(format!("Airy operator not positive definite at row {i}"));
                    }
                    l[i * w + b] = s.sqrt();
                } else {
                    l[i * w + (j + b - i)] = s / l[j * w + b];
                }
            }
        }
        Ok(Band { n, b, l })
    }
    pub fn solve(&self, x: &mut [f64]) {
        let (n, b, w) = (self.n, self.b, self.b + 1);
        for i in 0..n {
            let lo = i.saturating_sub(b);
            let row = &self.l[i * w + (lo + b - i)..i * w + b + 1];
            let s = x[i] - dot(&row[..i - lo], &x[lo..i]);
            x[i] = s / row[i - lo];
        }
        self.backward(x);
    }
    /// `solve` for the right-hand side `rhs(0), rhs(1), ...`, each entry built inside the forward
    /// sweep just before its row is eliminated (`x` is only written). Every row of that sweep waits
    /// on the previous row's division, so the right-hand side's independent loads and arithmetic
    /// fill that wait instead of costing a pass of their own. Same operations, same order, same
    /// bits as filling `x` first and calling `solve`.
    #[inline]
    pub fn solve_with(&self, x: &mut [f64], mut rhs: impl FnMut(usize) -> f64) {
        let (n, b, w) = (self.n, self.b, self.b + 1);
        for i in 0..n {
            let lo = i.saturating_sub(b);
            let row = &self.l[i * w + (lo + b - i)..i * w + b + 1];
            let s = rhs(i) - dot(&row[..i - lo], &x[lo..i]);
            x[i] = s / row[i - lo];
        }
        self.backward(x);
    }
    #[inline]
    fn backward(&self, x: &mut [f64]) {
        let (n, b, w) = (self.n, self.b, self.b + 1);
        // L' x = y by columns of L' (rows of L): each solved x[i] is pushed into the rows above,
        // so both sweeps stream the same contiguous factor once
        for i in (0..n).rev() {
            let lo = i.saturating_sub(b);
            let row = &self.l[i * w + (lo + b - i)..i * w + b + 1];
            let xi = x[i] / row[i - lo];
            x[i] = xi;
            for (xj, lij) in x[lo..i].iter_mut().zip(&row[..i - lo]) {
                *xj -= lij * xi;
            }
        }
    }
}

// ------------------------------------------------------------------------ bodies
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Edge {
    /// Displacement zero outside the domain (a drum head on its bearing edge).
    Clamped,
    /// Natural boundary of the discrete energy (a cymbal's rim and centre hole).
    Free,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Nonlinear {
    None,
    /// Uniform tension increase from the mean stretch (clamped membranes).
    Berger,
    /// Full von Karman / shallow-shell coupling through the Airy stress function (free plates).
    VonKarman,
}

/// A circular plate or membrane. Thickness and rest height are functions of radius.
pub struct PlateSpec<'a> {
    pub radius: f64,
    pub hole: f64,
    pub grid: f64,
    pub edge: Edge,
    pub thickness: &'a dyn Fn(f64) -> f64,
    pub rest_height: &'a dyn Fn(f64) -> f64,
    pub young: f64,
    pub poisson: f64,
    pub density: f64,
    /// In-plane pretension per unit length [N/m].
    pub tension: f64,
    pub bending: bool,
    pub nonlinear: Nonlinear,
    /// Two (frequency Hz, T60 s) points fixing s0 + s1 beta^2.
    pub t60: [(f64, f64); 2],
    /// Felt support to ground over an annulus: (r0, r1, stiffness N/m per m^2, damping N s/m per m^2).
    pub support: Option<(f64, f64, f64, f64)>,
    /// Extra viscous damping to ground over a patch (centre x, y, radius, N s/m per m^2): a pillow or gel.
    pub muffle: Option<(f64, f64, f64, f64)>,
}

#[derive(Clone)]
struct VonKarman {
    eh: f64,
    nb: Vec<[Option<usize>; 9]>, // C, E, W, N, S, NE, NW, SE, SW
    l0: Vec<f64>,              // l(z0) at core
    z0: Vec<f64>,              // rest height per node
    airy: Band,
    f: Vec<f64>,
    q: Vec<f64>,
    curv: Vec<(f64, f64, f64)>,
}

#[derive(Clone)]
pub struct Body {
    pub name: String,
    pub h: f64,
    pub n: usize,
    pub xy: Vec<[f64; 2]>,
    /// Supplied mid-surface height in the body's outward displacement coordinates [m].
    pub rest_z: Vec<f64>,
    grid_index: BTreeMap<(i64, i64), usize>,
    pub mass: Vec<f64>,
    s0: Vec<f64>,
    s1: f64,
    lap: Csr,
    stiff: Csr,
    viscosity: Option<ViscousBending>,
    tension: f64,
    berger: f64,
    vk: Option<VonKarman>,
    pub w: Vec<f64>,
    pub wp: Vec<f64>,
    pub wn: Vec<f64>,
    pub f: Vec<f64>,
    /// Constant external force per node (snare-wire preload).
    pub bias: Vec<f64>,
    tmp: Vec<f64>,
    /// L w at the current step, and L w at the previous one (reused: one Laplacian product per step).
    lw: Vec<f64>,
    lwp: Vec<f64>,
    g: Vec<f64>,
    psi: f64,
    psi_next: f64,
    sav: bool,
    /// sqrt(2 V) at the previous step, and the energy the loss terms removed in the last step
    /// (the most the auxiliary variable may be raised by when it is pulled back to sqrt(2 V)).
    r_prev: f64,
    budget: f64,
    /// Step size the cached update coefficients `cf` belong to.
    kc: f64,
    cf: Vec<[f64; 4]>,
    /// 3-D placement for microphones: origin and outward normal (+1 up, -1 down).
    pub origin: [f64; 3],
    pub facing: f64,
    pub radiates: bool,
    pub area: f64,
    /// Air mass loading the volume motion (one side: 8 rho a^3 / 3 for a baffled piston) [kg].
    pub air_mass: f64,
    /// Exterior radiation resistance and air mass of non-volume head shapes.
    head_radiation: Option<Box<radiation::ModalAirLoad>>,
    radiated_j: f64,
}

/// C = G'G is the assembled curvature-rate damping matrix. Retaining G lets
/// physical dissipation be accumulated as squares, without clipping negative
/// quadratic roundoff. Three G rows belong to each plate node (empty outside
/// the curvature integration points), which external steppers also retain.
#[derive(Clone)]
struct ViscousBending {
    source: PlaneStressViscosity,
    angle: f64,
    c: Csr,
    root: Csr,
    scratch: Vec<f64>,
    rate_bound: f64,
    dissipated_j: f64,
    last_j: f64,
}

fn sigmas(t60: [(f64, f64); 2], beta2: &dyn Fn(f64) -> f64) -> (f64, f64) {
    let (fa, ta) = t60[0];
    let (fb, tb) = t60[1];
    let (sa, sb) = (6.907755 / ta, 6.907755 / tb);
    let (ba, bb) = (beta2(fa), beta2(fb));
    let s1 = if (bb - ba).abs() > 0.0 { ((sb - sa) / (bb - ba)).max(0.0) } else { 0.0 };
    ((sa - s1 * ba).max(0.0), s1)
}

impl Body {
    pub fn plate(name: &str, s: &PlateSpec) -> Result<Body, String> {
        Self::plate_impl(name, s, None, true)
    }

    /// Same operators with no surface loss at all: `s.t60` is ignored, and
    /// support and muffle dashpots remain. For a surface whose loss the caller
    /// declares a stated unknown entered as zero.
    pub fn plate_without_surface_loss(name: &str, s: &PlateSpec) -> Result<Body, String> {
        Self::plate_impl(name, s, None, false)
    }

    /// Same elastic/geometry operators with a sourced reduced plane-stress
    /// Kelvin--Voigt bending law. Replaces the T60 prescription completely;
    /// support and muffle dashpots remain independent. The receiving instrument
    /// must admit the source band and damping CFL at its actual internal clock.
    pub fn plate_with_viscosity(name: &str, s: &PlateSpec,
        source: &PlaneStressViscosity, angle_rad: f64) -> Result<Body, String> {
        if !s.bending { return Err("bending viscosity requires the plate bending operator".into()); }
        Self::plate_impl(name, s, Some((source, angle_rad)), false)
    }

    fn plate_impl(name: &str, s: &PlateSpec,
        viscosity: Option<(&PlaneStressViscosity, f64)>, t60_loss: bool) -> Result<Body, String> {
        let h = s.grid;
        let m = (s.radius / h).ceil() as i64 + 3;
        let inside = |i: i64, j: i64| {
            let r = ((i * i + j * j) as f64).sqrt() * h;
            match s.edge {
                Edge::Clamped => r < s.radius - 1e-9 && r >= s.hole,
                Edge::Free => r <= s.radius + 1e-9 && r >= s.hole - 1e-9,
            }
        };
        let mut core = vec![];
        let mut dom = std::collections::BTreeSet::new();
        for j in -m..=m {
            for i in -m..=m {
                match s.edge {
                    Edge::Clamped => {
                        if inside(i, j) {
                            dom.insert((j, i));
                        }
                    }
                    Edge::Free => {
                        let all = (-1..=1).all(|dj| (-1..=1).all(|di| inside(i + di, j + dj)));
                        if all {
                            core.push((i, j));
                            for dj in -1..=1 {
                                for di in -1..=1 {
                                    dom.insert((j + dj, i + di));
                                }
                            }
                        }
                    }
                }
            }
        }
        // row-major numbering (j outer, i inner) keeps the Airy band narrow
        let mut grid_index = BTreeMap::new();
        let mut xy = vec![];
        for &(j, i) in &dom {
            grid_index.insert((i, j), xy.len());
            xy.push([i as f64 * h, j as f64 * h]);
        }
        let n = xy.len();
        if n < 9 {
            return Err(format!("{name}: grid {h} m leaves only {n} nodes"));
        }
        let at = |i: i64, j: i64| grid_index.get(&(i, j)).copied();
        let r_of = |p: [f64; 2]| (p[0] * p[0] + p[1] * p[1]).sqrt();
        let thick: Vec<f64> = xy.iter().map(|&p| (s.thickness)(r_of(p))).collect();
        let mut mass: Vec<f64> = thick.iter().map(|t| s.density * t * h * h).collect();
        // graph Laplacian (membrane energy T/2 sum (w_i - w_j)^2), Dirichlet outside when clamped
        let mut lrows = vec![BTreeMap::new(); n];
        for (&(i, j), &a) in &grid_index {
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                match at(i + di, j + dj) {
                    Some(b) => {
                        if a < b {
                            add(&mut lrows, a, a, 1.0);
                            add(&mut lrows, b, b, 1.0);
                            add(&mut lrows, a, b, -1.0);
                            add(&mut lrows, b, a, -1.0);
                        }
                    }
                    None => {
                        if s.edge == Edge::Clamped {
                            add(&mut lrows, a, a, 1.0);
                        }
                    }
                }
            }
        }
        let lap = Csr::from_rows(&lrows);
        // bending energy D/2 h^2 [wxx^2 + wyy^2 + 2 nu wxx wyy + 2 (1 - nu) wxy^2] per curvature point
        let mut krows = vec![BTreeMap::new(); n];
        let nu = s.poisson;
        let curv_pts: Vec<(i64, i64)> = match s.edge {
            Edge::Clamped => grid_index.keys().copied().collect(),
            Edge::Free => core.clone(),
        };
        let stencil = |i: i64, j: i64| -> [Vec<(Option<usize>, f64)>; 3] {
            let h2 = h * h;
            [
                vec![(at(i - 1, j), 1.0 / h2), (at(i, j), -2.0 / h2), (at(i + 1, j), 1.0 / h2)],
                vec![(at(i, j - 1), 1.0 / h2), (at(i, j), -2.0 / h2), (at(i, j + 1), 1.0 / h2)],
                vec![
                    (at(i + 1, j + 1), 0.25 / h2),
                    (at(i - 1, j + 1), -0.25 / h2),
                    (at(i + 1, j - 1), -0.25 / h2),
                    (at(i - 1, j - 1), 0.25 / h2),
                ],
            ]
        };
        if s.bending {
            let q = [[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 2.0 * (1.0 - nu)]];
            for &(i, j) in &curv_pts {
                let r = ((i * i + j * j) as f64).sqrt() * h;
                let t = (s.thickness)(r);
                let d = s.young * t * t * t / (12.0 * (1.0 - nu * nu));
                let b = stencil(i, j);
                for a in 0..3 {
                    for c in 0..3 {
                        if q[a][c] == 0.0 {
                            continue;
                        }
                        for &(pa, va) in &b[a] {
                            for &(pc, vc) in &b[c] {
                                if let (Some(pa), Some(pc)) = (pa, pc) {
                                    add(&mut krows, pa, pc, d * h * h * q[a][c] * va * vc);
                                }
                            }
                        }
                    }
                }
            }
        }
        let viscosity = if let Some((source, angle)) = viscosity {
            let mut roots = vec![BTreeMap::new(); 3 * n];
            for &(i, j) in &curv_pts {
                let r = ((i * i + j * j) as f64).sqrt() * h;
                let d = source.bending_matrix_n_m_s((s.thickness)(r), angle)
                    .map_err(|e| e.to_string())?;
                if d.iter().all(|v| *v == 0.0) { continue; }
                let rows: Vec<BTreeMap<usize, f64>> = (0..3).map(|a|
                    (0..3).map(|b| (b, d[3*a+b])).collect()).collect();
                let factor = Band::factor(&rows)?;
                let b = stencil(i, j);
                let node = at(i, j).expect("curvature point belongs to the plate");
                // Engineering curvature is [w_xx,w_yy,2w_xy]. The square
                // root includes sqrt(area)=h, so G'G has units N s/m.
                for column in 0..3 {
                    for a in column..3 {
                        if a - column > factor.b { continue; }
                        let l = factor.l[a * (factor.b + 1) + column + factor.b - a];
                        let scale = h * l * if a == 2 { 2.0 } else { 1.0 };
                        for &(p, value) in &b[a] {
                            if let Some(p) = p { add(&mut roots, 3 * node + column, p, scale * value); }
                        }
                    }
                }
            }
            let mut rows = vec![BTreeMap::new(); n];
            for root in &roots {
                for (&i, &a) in root {
                    for (&j, &b) in root { add(&mut rows, i, j, a * b); }
                }
            }
            if rows.iter().chain(&roots).flat_map(|r| r.values()).any(|v| !v.is_finite()) {
                return Err("assembled bending viscosity is not representable".into());
            }
            let rate_bound = rows.iter().enumerate().map(|(i, row)| row.iter()
                .map(|(&j, value)| value.abs() / mass[i].sqrt() / mass[j].sqrt()).sum::<f64>())
                .fold(0.0_f64, f64::max);
            if !rate_bound.is_finite() { return Err("bending-viscosity CFL bound is not representable".into()); }
            Some(ViscousBending { source: source.clone(), angle, c: Csr::from_rows(&rows),
                root: Csr::from_rows(&roots), scratch: vec![0.0; 3*n], rate_bound,
                dissipated_j: 0.0, last_j: 0.0 })
        } else { None };
        // supports and muffling: springs/dashpots to ground per node
        let mut s0_extra = vec![0.0; n];
        for (k, &p) in xy.iter().enumerate() {
            let r = r_of(p);
            if let Some((r0, r1, kk, cc)) = s.support {
                if r >= r0 && r <= r1 {
                    add(&mut krows, k, k, kk * h * h);
                    s0_extra[k] += cc * h * h / (2.0 * mass[k]);
                }
            }
            if let Some((mx, my, mr, cc)) = s.muffle {
                if ((p[0] - mx).powi(2) + (p[1] - my).powi(2)).sqrt() <= mr {
                    s0_extra[k] += cc * h * h / (2.0 * mass[k]);
                }
            }
        }
        let stiff = Csr::from_rows(&krows);
        // loss from the dispersion of the thinnest region (the bow of a cymbal, the head itself)
        let t_min = thick.iter().cloned().fold(f64::INFINITY, f64::min);
        let rho_t = s.density * t_min;
        let dd = if s.bending { s.young * t_min.powi(3) / (12.0 * (1.0 - nu * nu)) } else { 0.0 };
        let (tt, ddc) = (s.tension / rho_t, dd / rho_t);
        let beta2 = move |f: f64| {
            let om2 = (2.0 * PI * f).powi(2);
            if ddc > 0.0 { (-tt + (tt * tt + 4.0 * ddc * om2).sqrt()) / (2.0 * ddc) } else { om2 / tt }
        };
        let (s0, s1) = if t60_loss { sigmas(s.t60, &beta2) } else { (0.0, 0.0) };
        let s0v: Vec<f64> = (0..n).map(|k| s0 + s0_extra[k]).collect();
        let area = PI * (s.radius * s.radius - s.hole * s.hole);
        let t_mean = thick.iter().sum::<f64>() / n as f64;
        let berger = if s.nonlinear == Nonlinear::Berger { s.young * t_mean / (2.0 * area * (1.0 - nu)) } else { 0.0 };
        let rest_z: Vec<f64> = xy.iter().map(|&p| (s.rest_height)(r_of(p))).collect();
        if rest_z.iter().any(|z| !z.is_finite()) {
            return Err(format!("{name}: nonfinite rest-surface height"));
        }
        let vk = if s.nonlinear == Nonlinear::VonKarman {
            if s.edge != Edge::Free {
                return Err(format!("{name}: von Karman coupling is implemented for free edges (stress-free rim)"));
            }
            let z0 = rest_z.clone();
            let mut cidx = BTreeMap::new();
            let mut sorted = core.clone();
            sorted.sort_by_key(|&(i, j)| (j, i));
            for (idx, &(i, j)) in sorted.iter().enumerate() {
                cidx.insert((i, j), idx);
            }
            let nb: Vec<[Option<usize>; 9]> = sorted
                .iter()
                .map(|&(i, j)| {
                    [at(i, j), at(i + 1, j), at(i - 1, j), at(i, j + 1), at(i, j - 1), at(i + 1, j + 1), at(i - 1, j + 1), at(i + 1, j - 1), at(i - 1, j - 1)]
                })
                .collect();
            let curv = |v: &[f64], s9: &[Option<usize>; 9]| {
                let g = |o: Option<usize>| o.map_or(0.0, |k| v[k]);
                let h2 = h * h;
                let xx = (g(s9[1]) - 2.0 * g(s9[0]) + g(s9[2])) / h2;
                let yy = (g(s9[3]) - 2.0 * g(s9[0]) + g(s9[4])) / h2;
                let xy_ = (g(s9[5]) - g(s9[6]) - g(s9[7]) + g(s9[8])) / (4.0 * h2);
                (xx, yy, xy_)
            };
            let l0: Vec<f64> = nb.iter().map(|s9| {
                let (a, b, c) = curv(&z0, s9);
                a * b - c * c
            }).collect();
            // Airy energy 1/2 sum (Lap_h F)^2 over every node where it is nonzero; F = 0 off the core
            let nc = sorted.len();
            let mut arows = vec![BTreeMap::new(); nc];
            let h2 = h * h;
            for (&(i, j), _) in &grid_index {
                let mut row: Vec<(usize, f64)> = vec![];
                if let Some(&c) = cidx.get(&(i, j)) {
                    row.push((c, -4.0 / h2));
                }
                for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                    if let Some(&c) = cidx.get(&(i + di, j + dj)) {
                        row.push((c, 1.0 / h2));
                    }
                }
                for &(a, va) in &row {
                    for &(b, vb) in &row {
                        add(&mut arows, a, b, va * vb);
                    }
                }
            }
            let airy = Band::factor(&arows)?;
            Some(VonKarman { eh: s.young * t_mean, nb, l0, z0, airy, f: vec![0.0; nc], q: vec![0.0; nc], curv: vec![(0.0, 0.0, 0.0); nc] })
        } else {
            None
        };
        // edge nodes of a free plate outside the true disc carry real mass but a staircase rim
        for m_ in mass.iter_mut() {
            *m_ = m_.max(1e-12);
        }
        Ok(Body {
            name: name.into(), h, n, xy, rest_z, grid_index, mass, s0: s0v, s1, lap, stiff, viscosity, tension: s.tension, berger, vk,
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: true, area, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
        })
    }

    /// A 1-D string (a snare wire) of `n` nodes spaced `h`, both ends fixed beyond the last node.
    pub fn string(name: &str, n: usize, h: f64, mass_per_m: f64, tension: f64, t60: f64, x0: f64, y: f64) -> Body {
        let mut lrows = vec![BTreeMap::new(); n];
        for a in 0..n {
            add(&mut lrows, a, a, 2.0);
            if a + 1 < n {
                add(&mut lrows, a, a + 1, -1.0);
                add(&mut lrows, a + 1, a, -1.0);
            }
        }
        let xy = (0..n).map(|k| [x0 + k as f64 * h, y]).collect();
        Body {
            name: name.into(), h, n, xy, rest_z: vec![0.0; n], grid_index: BTreeMap::new(), mass: vec![mass_per_m * h; n],
            s0: vec![6.907755 / t60; n], s1: 0.0, lap: Csr::from_rows(&lrows), stiff: Csr::from_rows(&vec![BTreeMap::new(); n]), viscosity: None,
            tension: tension / h, berger: 0.0, vk: None, // tension/h: the chain energy T/2 sum (dw)^2/h
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
        }
    }

    /// A lumped rigid mass (stick, beater, hi-hat rod).
    pub fn point(name: &str, mass: f64) -> Body {
        Body {
            name: name.into(), h: 1.0, n: 1, xy: vec![[0.0, 0.0]], rest_z: vec![0.0], grid_index: BTreeMap::new(), mass: vec![mass],
            s0: vec![0.0], s1: 0.0, lap: Csr::from_rows(&[BTreeMap::new()]), stiff: Csr::from_rows(&[BTreeMap::new()]), viscosity: None,
            tension: 0.0, berger: 0.0, vk: None,
            w: vec![0.0], wp: vec![0.0], wn: vec![0.0], f: vec![0.0], bias: vec![0.0], tmp: vec![0.0], lw: vec![0.0], lwp: vec![0.0], g: vec![0.0], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
        }
    }

    /// Attach the exterior load of a clamped full circular head of `radius`
    /// (its nodes centred on the rim's centre) from `table` in `gas`. Mass
    /// and stiffness are this body's own; returns the receipt lines. The
    /// optional compact net-volume resistance [Pa s/m^3] is only compared
    /// with the baked volume shapes in the receipt.
    pub fn attach_head_radiation(&mut self, table: &radiation::RadiationTable, radius: f64, gas: &GasState,
        compact_volume_resistance: Option<f64>) -> Result<Vec<String>, String> {
        // Air loading is the same physics whatever the film-loss source: a
        // sourced bending viscosity composes with it (its force is explicit,
        // the air mass implicit). Shallow shells are not clamped films.
        if self.vk.is_some() || self.n < 2 {
            return Err(format!("{}: head radiation applies to clamped tensioned heads, not shells", self.name));
        }
        let (stiff, lap, tension) = (&self.stiff, &self.lap, self.tension);
        let stiffness = |v: &[f64]| {
            let mut out = vec![0.0; v.len()];
            stiff.mul_add(v, 1.0, &mut out);
            lap.mul_add(v, tension, &mut out);
            out
        };
        let load = radiation::ModalAirLoad::new(table, radius, self.h, &self.xy, &self.mass, &stiffness, gas, compact_volume_resistance)?;
        let receipt = load.receipt.clone();
        self.head_radiation = Some(Box::new(load));
        self.kc = 0.0; // the implicit air mass is prepared with the update coefficients
        Ok(receipt)
    }

    pub fn head_radiation(&self) -> Option<&radiation::ModalAirLoad> {
        self.head_radiation.as_deref()
    }

    /// Energy the head-radiation resistance removed, at centred velocity [J].
    pub fn head_radiated_j(&self) -> f64 {
        self.radiated_j
    }

    pub fn node_at(&self, x: f64, y: f64) -> Option<usize> {
        self.grid_index.get(&((x / self.h).round() as i64, (y / self.h).round() as i64)).copied()
    }

    /// Finite-difference weights of dw/dx and dw/dy at node `i`: central between the two
    /// neighbours, one-sided where the plate ends on one side (empty where it ends on both).
    pub fn slope(&self, i: usize) -> [Vec<(usize, f64)>; 2] {
        let (gi, gj) = ((self.xy[i][0] / self.h).round() as i64, (self.xy[i][1] / self.h).round() as i64);
        let at = |di: i64, dj: i64| self.grid_index.get(&(gi + di, gj + dj)).copied();
        let ih = 1.0 / self.h;
        let d = |p: Option<usize>, m: Option<usize>| match (p, m) {
            (Some(p), Some(m)) => vec![(p, 0.5 * ih), (m, -0.5 * ih)],
            (Some(p), None) => vec![(p, ih), (i, -ih)],
            (None, Some(m)) => vec![(i, ih), (m, -ih)],
            (None, None) => vec![],
        };
        [d(at(1, 0), at(-1, 0)), d(at(0, 1), at(0, -1))]
    }

    /// A point contact uses bilinear grid weights, preserving its position and total force.
    /// A finite felt face uses equally weighted nodes within `radius`. At an incomplete
    /// boundary cell, retain the nearest-node fallback rather than extrapolating off the body.
    pub fn footprint(&self, x: f64, y: f64, radius: f64) -> Vec<(usize, f64)> {
        if radius == 0.0 && (x / self.h).abs() < i64::MAX as f64 && (y / self.h).abs() < i64::MAX as f64 {
            let (gx, gy) = (x / self.h, y / self.h);
            let (i, j) = (gx.floor() as i64, gy.floor() as i64);
            let (tx, ty) = (gx - i as f64, gy - j as f64);
            let corners = [
                ((i, j), (1.0 - tx) * (1.0 - ty)),
                ((i + 1, j), tx * (1.0 - ty)),
                ((i, j + 1), (1.0 - tx) * ty),
                ((i + 1, j + 1), tx * ty),
            ];
            let mut point = Vec::with_capacity(4);
            let mut complete = true;
            for (cell, weight) in corners {
                if weight == 0.0 {
                    continue;
                }
                if let Some(&node) = self.grid_index.get(&cell) {
                    point.push((node, weight));
                } else {
                    complete = false;
                    break;
                }
            }
            if complete && !point.is_empty() {
                // The same weights read the contact displacement and distribute its
                // force, so interpolation and loading remain reciprocal in work.
                return point;
            }
        }
        let mut v: Vec<usize> = (0..self.n)
            .filter(|&k| ((self.xy[k][0] - x).powi(2) + (self.xy[k][1] - y).powi(2)).sqrt() <= radius)
            .collect();
        if v.is_empty() {
            let k = (0..self.n)
                .min_by(|&a, &b| {
                    let da = (self.xy[a][0] - x).powi(2) + (self.xy[a][1] - y).powi(2);
                    let db = (self.xy[b][0] - x).powi(2) + (self.xy[b][1] - y).powi(2);
                    da.partial_cmp(&db).unwrap()
                })
                .unwrap();
            v.push(k);
        }
        let w = 1.0 / v.len() as f64;
        v.into_iter().map(|k| (k, w)).collect()
    }

    pub fn nodes(&self) -> usize {
        self.n
    }

    /// Linear elastic forces and the constant bias into `f` (overwrites).
    fn linear_forces(&mut self) {
        // the one Laplacian product of the step: tension, Berger stretch and loss all reuse it
        for x in self.lw.iter_mut() {
            *x = 0.0;
        }
        self.lap.mul_add(&self.w, 1.0, &mut self.lw);
        let t = self.tension;
        for ((f, b), l) in self.f.iter_mut().zip(&self.bias).zip(&self.lw) {
            *f = b - t * l;
        }
        if !self.stiff.col.is_empty() {
            self.stiff.mul_add(&self.w, -1.0, &mut self.f);
        }
    }

    /// Nonlinear potential V(w) >= 0 and its gradient into `g`; returns V.
    fn nonlinear_gradient(&mut self) -> f64 {
        for x in self.g.iter_mut() {
            *x = 0.0;
        }
        let mut v = 0.0;
        if self.berger > 0.0 {
            // V = berger/4 (w' L w)^2, grad = berger (w' L w) L w
            let sq: f64 = self.w.iter().zip(&self.lw).map(|(a, b)| a * b).sum();
            v += 0.25 * self.berger * sq * sq;
            let c = self.berger * sq;
            for (g, l) in self.g.iter_mut().zip(&self.lw) {
                *g += c * l;
            }
        }
        if let Some(vk) = self.vk.as_mut() {
            let h2 = self.h * self.h;
            let w = &self.w;
            let VonKarman { eh, nb, l0, z0, airy, f, q, curv } = vk;
            let eh = *eh;
            let u = |o: Option<usize>| o.map_or(0.0, |k| w[k] + z0[k]);
            // q = -E t (l(w + z0) - l(z0)) per core node, built row by row inside the Airy solve's
            // forward sweep (it fills that sweep's latency chain); A F = q
            airy.solve_with(f, |c| {
                let s9 = &nb[c];
                let xx = (u(s9[1]) - 2.0 * u(s9[0]) + u(s9[2])) / h2;
                let yy = (u(s9[3]) - 2.0 * u(s9[0]) + u(s9[4])) / h2;
                let xy = (u(s9[5]) - u(s9[6]) - u(s9[7]) + u(s9[8])) / (4.0 * h2);
                let qc = -eh * (xx * yy - xy * xy - l0[c]);
                q[c] = qc;
                curv[c] = (xx, yy, xy);
                qc
            });
            // V = h^2/(2 E t) F' A F with A F = q; the force h^2 (dq/dw)' F is -grad V
            v += h2 / (2.0 * eh) * f.iter().zip(q.iter()).map(|(a, b)| a * b).sum::<f64>();
            for (c, s9) in nb.iter().enumerate() {
                let (xx, yy, xy) = curv[c];
                let fc = f[c]; // h^2 (area) times the 1/h^2 of each stencil weight
                let mut put = |o: Option<usize>, val: f64| {
                    if let Some(k) = o {
                        self.g[k] -= fc * val;
                    }
                };
                put(s9[1], yy);
                put(s9[2], yy);
                put(s9[0], -2.0 * yy - 2.0 * xx);
                put(s9[3], xx);
                put(s9[4], xx);
                let gx = -0.5 * xy;
                put(s9[5], gx);
                put(s9[6], -gx);
                put(s9[7], -gx);
                put(s9[8], gx);
            }
        }
        v.max(0.0)
    }

    /// Leapfrog step from the forces in `f`. The nonlinear potential enters through a scalar
    /// auxiliary variable psi = sqrt(2 V) (Bilbao, Ducceschi and Zama): the update stays explicit
    /// (one Sherman-Morrison rank-one correction) and conserves a discrete energy
    /// 1/2 |dw/dt|^2_M + 1/2 w'Kw + 1/2 psi^2, so the nonlinearity cannot destabilise it.
    /// Per-node update coefficients for step `k`: [1/a, 2m/k^2, (1 - s0 k) m/k^2, 2 s1 m/k] with
    /// a = m (1 + s0 k)/k^2. Mass and loss are fixed once stepping starts.
    fn coefs(&mut self, k: f64) {
        if self.kc == k && self.cf.len() == self.n {
            return;
        }
        let k2 = k * k;
        self.cf = (0..self.n)
            .map(|i| {
                let (m, s0k) = (self.mass[i], self.s0[i] * k);
                [k2 / (m * (1.0 + s0k)), 2.0 * m / k2, (1.0 - s0k) * m / k2, 2.0 * self.s1 * m / k]
            })
            .collect();
        self.kc = k;
        if let Err(e) = self.prepare_head_radiation(k) {
            // Instrument::validate_head_radiation refuses this at build time
            // for the instrument's own clock; another clock cannot reach here.
            panic!("{e}");
        }
    }

    /// Factor the implicit air-mass system at step `k` (with the volume
    /// column when the piston air mass is present).
    fn prepare_head_radiation(&mut self, k: f64) -> Result<(), String> {
        let Some(load) = self.head_radiation.as_deref_mut() else { return Ok(()) };
        let mu = if self.air_mass > 0.0 { self.air_mass / (self.area * self.area * k * k) } else { 0.0 };
        let ia: Vec<f64> = (0..self.n).map(|i| k * k / (self.mass[i] * (1.0 + self.s0[i] * k))).collect();
        let volume = (mu > 0.0).then_some((self.h * self.h, mu));
        load.prepare(k, &ia, volume).map_err(|e| format!("{}: {e}", self.name))
    }

    fn predict(&mut self, k: f64) {
        self.coefs(k);
        let n = self.n;
        if let Some(viscous) = &self.viscosity {
            for i in 0..n { self.tmp[i] = self.w[i] - self.wp[i]; }
            // Backward velocity preserves the existing diagonal contact and
            // SAV compliance. Its additional explicit stability bound is
            // dt² lambda_K + 2 dt lambda_C < 4; see discrete_energy below.
            viscous.c.mul_add(&self.tmp, -1.0 / k, &mut self.f);
        }
        if self.s1 > 0.0 {
            // Lap_h (w - wp) = -(L w - L wp) / h^2 from the two stored products
            let ih2 = 1.0 / (self.h * self.h);
            for ((t, a), b) in self.tmp.iter_mut().zip(&self.lw).zip(&self.lwp) {
                *t = -(a - b) * ih2;
            }
        }
        let nl = self.berger > 0.0 || self.vk.is_some();
        let mut scale = 0.0;
        if nl {
            let v = self.nonlinear_gradient();
            let root = (2.0 * v).sqrt();
            if self.vk.is_some() {
                // psi is advanced along g only, so it drifts from sqrt(2V): on a live grid it wandered
                // between 0.5 and 1.5 times sqrt(2V) within 2 s, scaling the shell's curvature
                // stiffness (which a shallow shell carries entirely in V) by that factor. Pull it back
                // to sqrt(2V) at the half step without adding energy the losses did not remove:
                // lowering psi only removes energy; raising it spends at most the last step's loss.
                let target = if self.r_prev > 0.0 { 0.5 * (root + self.r_prev) } else { root };
                if self.psi > target {
                    self.psi = target;
                } else {
                    self.psi = (self.psi * self.psi + 2.0 * self.budget).sqrt().min(target);
                }
                self.r_prev = root;
            }
            if root > 0.0 {
                let ir = 1.0 / root;
                for x in self.g.iter_mut() {
                    *x *= ir;
                }
                scale = 1.0;
            } else {
                for x in self.g.iter_mut() {
                    *x = 0.0;
                }
            }
        }
        let k2 = k * k;
        let gwp: f64 = if scale > 0.0 { dot(&self.g, &self.wp) } else { 0.0 };
        let sc = self.psi - 0.25 * gwp;
        // air mass on the volume motion: (M_r / A^2) u u' with u = h^2 per node, implicit
        let h2 = self.h * self.h;
        let mu = if self.air_mass > 0.0 { self.air_mass / (self.area * self.area * k2) } else { 0.0 };
        let vol_prev = if mu > 0.0 { self.w.iter().zip(&self.wp).map(|(a, b)| 2.0 * a - b).sum::<f64>() * h2 } else { 0.0 };
        let extra = if mu > 0.0 { mu * h2 * vol_prev } else { 0.0 };
        if let Some(load) = self.head_radiation.as_deref_mut() {
            load.load_rhs(&self.w[..n], &self.wp[..n], k);
        }
        let (mut gy, mut gzg, mut uy, mut uzu, mut gzu) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (cf, w, wp, f, g) = (&self.cf[..n], &self.w[..n], &self.wp[..n], &self.f[..n], &self.g[..n]);
        let (wn, tmp) = (&mut self.wn[..n], &mut self.tmp[..n]);
        let s1 = self.s1 > 0.0;
        let modal: &[f64] = self.head_radiation.as_deref().map_or(&[], |load| &load.rhs);
        for i in 0..n {
            let [ia, cw, cwp, cs1] = cf[i];
            let mut rhs = cw * w[i] - cwp * wp[i] + f[i] + extra;
            if !modal.is_empty() {
                rhs += modal[i];
            }
            if s1 {
                rhs += cs1 * tmp[i];
            }
            if scale > 0.0 {
                rhs -= g[i] * sc;
            }
            let x = rhs * ia;
            wn[i] = x;
            if scale > 0.0 {
                let z = g[i] * ia;
                tmp[i] = z;
                gy += g[i] * x;
                gzg += g[i] * z;
                gzu += h2 * z;
            }
            if mu > 0.0 {
                uy += h2 * x;
                uzu += h2 * h2 * ia;
            }
        }
        if let Some(load) = self.head_radiation.as_deref_mut() {
            // Woodbury with U = [g, V], V = [u?, modal shapes], C = diag(1/4, C_V):
            // with Q = C_V^-1 + V'D^-1V factored, b = V'D^-1 g and y = V'D^-1 r,
            // c_g = (g'D^-1 r - b'Q^-1 y) / (4 + g'D^-1 g - b'Q^-1 b), c_V = Q^-1 (y - b c_g).
            let prepared = load.prepared.as_ref().filter(|p| p.k == k)
                .expect("prepared for this step with the update coefficients");
            let off = usize::from(prepared.with_volume);
            let size = off + load.cols.len();
            if off == 1 {
                (load.y[0], load.b[0]) = (uy, gzu);
            }
            for (r, col) in load.cols.iter().enumerate() {
                load.y[off + r] = dot(col, wn);
                load.b[off + r] = if scale > 0.0 { dot(col, tmp) } else { 0.0 };
            }
            radiation::cholesky_solve(&prepared.chol, size, &mut load.y[..size]);
            let cg = if scale > 0.0 {
                load.w2[..size].copy_from_slice(&load.b[..size]);
                radiation::cholesky_solve(&prepared.chol, size, &mut load.w2[..size]);
                let cg = (gy - dot(&load.b[..size], &load.y[..size])) / (4.0 + gzg - dot(&load.b[..size], &load.w2[..size]));
                for (y, w2) in load.y[..size].iter_mut().zip(&load.w2[..size]) {
                    *y -= w2 * cg;
                }
                cg
            } else {
                0.0
            };
            load.corr.fill(if off == 1 { h2 * load.y[0] } else { 0.0 });
            for (r, col) in load.cols.iter().enumerate() {
                let c = load.y[off + r];
                for (out, u) in load.corr.iter_mut().zip(col) {
                    *out += c * u;
                }
            }
            for i in 0..n {
                let zg = if scale > 0.0 { tmp[i] } else { 0.0 };
                wn[i] -= zg * cg + load.corr[i] * cf[i][0];
            }
            self.sav = scale > 0.0;
            return;
        }
        // Woodbury with U = [g, u], C = diag(1/4, mu): solve (C^-1 + U' D^-1 U) c = U' D^-1 r
        let (cg, cu) = match (scale > 0.0, mu > 0.0) {
            (false, false) => (0.0, 0.0),
            (true, false) => (gy / (4.0 + gzg), 0.0),
            (false, true) => (0.0, uy / (1.0 / mu + uzu)),
            (true, true) => {
                let (a11, a12, a22) = (4.0 + gzg, gzu, 1.0 / mu + uzu);
                let det = a11 * a22 - a12 * a12;
                ((a22 * gy - a12 * uy) / det, (a11 * uy - a12 * gy) / det)
            }
        };
        if cg != 0.0 || cu != 0.0 {
            let hc = h2 * cu;
            for i in 0..n {
                let zg = if scale > 0.0 { tmp[i] } else { 0.0 };
                wn[i] -= zg * cg + hc * cf[i][0];
            }
        }
        self.sav = scale > 0.0;
    }

    /// After contacts have adjusted `wn`: advance the auxiliary variable consistently.
    fn finish(&mut self, k: f64) {
        let mut bending_loss = 0.0;
        if let Some(viscous) = &mut self.viscosity {
            for i in 0..self.n { self.tmp[i] = self.wn[i] - self.wp[i]; }
            viscous.scratch.fill(0.0);
            viscous.root.mul_add(&self.tmp, 1.0, &mut viscous.scratch);
            // dt * ||G v_center||², with v_center=(wn-wp)/(2dt).
            // The squares are nonnegative without clipping or an empirical
            // decay correction. This is not an additional thermal update.
            bending_loss = dot(&viscous.scratch, &viscous.scratch) / (4.0 * k);
            viscous.last_j = bending_loss;
            viscous.dissipated_j += bending_loss;
        }
        if let Some(load) = self.head_radiation.as_deref_mut() {
            let radiated = load.radiated(&self.wn, k);
            self.radiated_j += radiated;
            bending_loss += radiated;
        }
        if self.sav && self.vk.is_some() {
            // energy the loss terms removed over this step (the budget for raising psi next step):
            // s0: sum m s0 (wn - wp)^2 / (2k); s1: s1/(k h^2) sum m (L(w - wp)) . (wn - wp), if positive
            let ih2 = 1.0 / (self.h * self.h);
            let mut d0 = 0.0;
            let mut d1 = 0.0;
            for i in 0..self.n {
                let dv = self.wn[i] - self.wp[i];
                d0 += self.mass[i] * self.s0[i] * dv * dv;
                d1 += self.mass[i] * (self.lw[i] - self.lwp[i]) * dv;
            }
            self.budget = d0 / (2.0 * k) + (self.s1 * ih2 / k * d1).max(0.0) + bending_loss;
        }
        if self.sav {
            let d: f64 = self.g.iter().zip(&self.wn).zip(&self.wp).map(|((g, a), b)| g * (a - b)).sum();
            self.psi_next = self.psi + 0.5 * d;
        } else {
            self.psi_next = self.psi;
        }
    }

    /// Displacement per unit force on node `i` within one step (k^2 / a).
    #[inline]
    fn compliance(&self, i: usize, k: f64) -> f64 {
        if self.kc == k && i < self.cf.len() { self.cf[i][0] } else { k * k / (self.mass[i] * (1.0 + self.s0[i] * k)) }
    }

    fn rotate(&mut self) {
        std::mem::swap(&mut self.wp, &mut self.w);
        std::mem::swap(&mut self.w, &mut self.wn);
        std::mem::swap(&mut self.lwp, &mut self.lw); // L wp of the next step is this step's L w
        self.psi = self.psi_next;
    }

    /// After setting w or wp directly (initial conditions), refresh the stored L wp.
    pub fn sync(&mut self) {
        for x in self.lwp.iter_mut() {
            *x = 0.0;
        }
        self.lap.mul_add(&self.wp, 1.0, &mut self.lwp);
    }

    /// Largest eigenvalue of M^-1 K for the LINEAR stiffness (power iteration on the symmetric
    /// form). The nonlinear part is handled by the auxiliary variable and needs no CFL headroom.
    pub fn max_rate2(&mut self) -> f64 {
        let n = self.n;
        let mut v: Vec<f64> = (0..n).map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5).collect();
        let mut lam = 0.0;
        let saved = self.w.clone();
        for _ in 0..300 {
            for i in 0..n {
                self.w[i] = v[i] / self.mass[i].sqrt();
            }
            let bias = std::mem::replace(&mut self.bias, vec![0.0; n]);
            self.linear_forces();
            self.bias = bias;
            let y: Vec<f64> = (0..n).map(|i| -self.f[i] / self.mass[i].sqrt()).collect();
            let norm = y.iter().map(|a| a * a).sum::<f64>().sqrt();
            let vn = v.iter().map(|a| a * a).sum::<f64>().sqrt();
            lam = norm / vn.max(1e-300);
            if norm == 0.0 {
                break;
            }
            v = y.iter().map(|a| a / norm).collect();
        }
        self.w = saved;
        lam
    }

    /// Largest eigenvalue of M^-1 (K + B), B the stiffness of the nonlinear potential linearised
    /// about the rest state (a shallow shell's curvature stiffness), by power iteration.
    pub fn max_rate2_with_shell(&mut self) -> f64 {
        let n = self.n;
        let saved = (self.w.clone(), self.g.clone());
        let mut v: Vec<f64> = (0..n).map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5).collect();
        let mut lam = 0.0;
        let eps = 1e-9;
        for _ in 0..400 {
            let vn = v.iter().map(|a| a * a).sum::<f64>().sqrt();
            for i in 0..n {
                self.w[i] = v[i] / vn / self.mass[i].sqrt();
            }
            let bias = std::mem::replace(&mut self.bias, vec![0.0; n]);
            self.linear_forces();
            self.bias = bias;
            // B x ~ grad V(eps x) / eps
            let xs = self.w.clone();
            for i in 0..n {
                self.w[i] = eps * xs[i];
            }
            self.nonlinear_gradient();
            let y: Vec<f64> = (0..n).map(|i| (-self.f[i] + self.g[i] / eps) / self.mass[i].sqrt()).collect();
            let norm = y.iter().map(|a| a * a).sum::<f64>().sqrt();
            lam = norm;
            if norm == 0.0 {
                break;
            }
            v = y;
        }
        self.w = saved.0;
        self.g = saved.1;
        lam
    }

    /// Diagnostics between steps: (psi, sqrt(2 V(w))) of the auxiliary variable.
    pub fn sav_state(&mut self) -> (f64, f64) {
        let v = self.nonlinear_gradient();
        (self.psi, (2.0 * v).sqrt())
    }

    /// Legacy s0/s1 prescription (including the minimum support s0), not the
    /// selected tensor law. Use `bending_viscosity` for sourced bending loss.
    pub fn loss(&self) -> (f64, f64) {
        (self.s0.iter().cloned().fold(f64::INFINITY, f64::min), self.s1)
    }

    /// Irreversible curvature-rate loss at centered velocity, in joules.
    /// Supports, contact, air and numerical dissipation are separate accounts.
    pub fn bending_loss_dissipated_j(&self) -> f64 {
        self.viscosity.as_ref().map_or(0.0, |v| v.dissipated_j)
    }
    pub fn bending_viscosity(&self) -> Option<(&PlaneStressViscosity, f64)> {
        self.viscosity.as_ref().map(|v| (&v.source, v.angle))
    }

    /// Staggered structural energy for the explicit curvature-viscosity step.
    /// For a linear isolated body it satisfies delta E = -dt ||G v_center||².
    /// The negative dt/4 v'C v correction is numerical stored energy, never
    /// reported as material heat. Positivity needs the admitted damping CFL.
    /// Coupling/contact/source work belongs to the complete instrument balance.
    pub fn discrete_energy(&self, k: f64) -> f64 {
        let velocity: Vec<f64> = self.w.iter().zip(&self.wp).map(|(w, wp)| (w - wp) / k).collect();
        let mut kw = vec![0.0; self.n];
        self.stiff.mul_add(&self.w, 1.0, &mut kw);
        self.lap.mul_add(&self.w, self.tension, &mut kw);
        let kinetic = self.mass.iter().zip(&velocity).map(|(m, v)| 0.5*m*v*v).sum::<f64>();
        let correction = self.viscosity.as_ref().map_or(0.0, |v| k * v.c.quad(&velocity) / 4.0);
        // Head radiation: the shapes' air kinetic energy, less the same
        // backward-velocity damping correction as the viscous law.
        let air = self.head_radiation.as_deref().map_or(0.0, |load| {
            let p: Vec<f64> = load.cols.iter().map(|c| dot(c, &velocity)).collect();
            let (mut a, mut b) = (vec![0.0; p.len()], vec![0.0; p.len()]);
            radiation::ModalAirLoad::apply(&load.mass, &p, &mut a);
            radiation::ModalAirLoad::apply(&load.damping, &p, &mut b);
            0.5 * dot(&p, &a) - k * dot(&p, &b) / 4.0
        });
        kinetic + air + 0.5 * dot(&self.wp, &kw) + 0.5 * self.psi * self.psi - correction
    }

    /// Kinetic + linear potential + auxiliary (nonlinear) energy at the half step (diagnostics).
    pub fn energy(&self, k: f64) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n {
            let v = (self.w[i] - self.wp[i]) / k;
            e += 0.5 * self.mass[i] * v * v;
        }
        let mid: Vec<f64> = (0..self.n).map(|i| 0.5 * (self.w[i] + self.wp[i])).collect();
        e + 0.5 * self.tension * self.lap.quad(&mid) + 0.5 * self.stiff.quad(&mid) + 0.5 * self.psi * self.psi
    }
}

// --------------------------------------------------------------------- couplings
/// Lumped adiabatic air between heads: p = -(rho c^2 / V) dV, dV = sum of outward volume.
#[derive(Clone)]
pub struct Cavity {
    pub heads: Vec<usize>,
    pub volume: f64,
    pub port: Option<Port>,
}

/// A hole in a head (a kick drum's front port): the air plug in it is an inertance
/// rho L_eff / A driven by the cavity pressure; it vents the cavity and radiates as a monopole.
#[derive(Clone)]
pub struct Port {
    pub area: f64,
    pub length: f64,
    /// Radiation (and flow) resistance on the volume velocity [Pa s / m^3].
    pub resistance: f64,
    /// Opening position and outward normal sign along z (for microphones).
    pub pos: [f64; 3],
    pub facing: f64,
    pub x: f64,
    pub u: f64,
    pub du: f64,
}

impl Port {
    pub fn new(radius: f64, head_thickness: f64, pos: [f64; 3], facing: f64) -> Port {
        let area = PI * radius * radius;
        // two flanged end corrections of 0.85 r each
        let length = head_thickness + 1.7 * radius;
        Port { area, length, resistance: 0.0, pos, facing, x: 0.0, u: 0.0, du: 0.0 }
    }
}

/// Small-ka monopole radiation resistance on the NET outward volume velocity of a group of
/// surfaces (both heads of one drum, which sit far closer together than a wavelength), at
/// frequency `f_ref`. Two heads moving the same way in space cancel and are barely damped.
#[derive(Clone)]
pub struct RadDamp {
    pub bodies: Vec<usize>,
    pub resistance: f64,
}

impl RadDamp {
    pub fn monopole(bodies: Vec<usize>, f_ref: f64, gas: &GasState) -> RadDamp {
        RadDamp { bodies, resistance: gas.density * (2.0 * PI * f_ref).powi(2) / (4.0 * PI * gas.sound_speed) }
    }
}

/// Linear spring/dashpot between two nodes (or one node and ground).
#[derive(Clone)]
pub struct Link {
    pub a: (usize, usize),
    pub b: Option<(usize, usize)>,
    pub k: f64,
    pub c: f64,
}

impl Link {
    /// Assemble a distributed grip with independent axial and angular impedance.
    ///
    /// `nodes` supplies positive quadrature weights (normalized here). Angular
    /// matrices are [xx, xy, yy] for slopes [dw/dx, dw/dy], in N m/rad and
    /// N m s/rad; axial coefficients are [N/m, N s/m]. The weighted affine fit
    /// separates translation, tilt and non-affine deformation. Its energy is
    /// K/2 * (mean(w)-anchor)^2 + K/2 * sum(weight*residual^2)
    /// + theta' * angular_stiffness * theta/2; damping has the same form.
    ///
    /// The returned links are ONE assembled PSD form. Individual pair links can
    /// have negative coefficients and must not be removed or used separately.
    /// This cold assembly uses the ordinary Link force/snapshot path unchanged.
    pub fn affine_mount(
        body_index: usize, body: &Body, nodes: &[(usize, f64)],
        anchor: Option<(usize, usize)>, axial: [f64; 2],
        angular_stiffness: [f64; 3], angular_damping: [f64; 3],
    ) -> Result<Vec<Link>, String> {
        let psd = |q: [f64; 3]| {
            let scale = q.iter().fold(0.0_f64, |s,x| s.max(x.abs()));
            q.iter().all(|x| x.is_finite()) && q[0] >= 0.0 && q[2] >= 0.0
                && (scale == 0.0 || (q[0]/scale)*(q[2]/scale) >= (q[1]/scale).powi(2))
        };
        if nodes.len() < 3 || !axial.iter().all(|x| x.is_finite() && *x >= 0.0)
            || !psd(angular_stiffness) || !psd(angular_damping) {
            return Err("affine mount: need >=3 nodes and finite nonnegative impedances (PSD angular matrices)".into());
        }
        let mut sum = 0.0;
        let mut seen = BTreeMap::new();
        for &(i,w) in nodes {
            if i >= body.n || !w.is_finite() || w <= 0.0
                || !body.xy[i].iter().all(|x| x.is_finite())
                || anchor == Some((body_index,i)) || seen.insert(i,()).is_some() {
                return Err("affine mount: invalid/duplicate attachment node, weight or self anchor".into());
            }
            sum += w;
        }
        if !sum.is_finite() { return Err("affine mount: weight sum overflow".into()); }
        let center: [f64; 2] = std::array::from_fn(|d| nodes.iter().map(|&(i,w)| (w/sum)*body.xy[i][d]).sum());
        let points: Vec<_> = nodes.iter().map(|&(i,w)| {
            (i,w/sum,[body.xy[i][0]-center[0],body.xy[i][1]-center[1]])
        }).collect();
        let mut cov = [0.0;3];
        for &(_,w,[x,y]) in &points {
            cov[0] += w*x*x; cov[1] += w*x*y; cov[2] += w*y*y;
        }
        let det = cov[0]*cov[2]-cov[1]*cov[1];
        if !det.is_finite() || det <= 1e-12*(cov[0]+cov[2]).powi(2) {
            return Err("affine mount: attachment points must span a well-conditioned plane".into());
        }
        let [a,b,d] = [cov[2]/det,-cov[1]/det,cov[0]/det];
        // R = covariance^-1. The pair correction is W X (R A R - K R) X' W.
        let correction = |k: f64, [u,v,w]: [f64;3]| [
            a*a*u+2.0*a*b*v+b*b*w-k*a,
            a*b*u+(a*d+b*b)*v+b*d*w-k*b,
            b*b*u+2.0*b*d*v+d*d*w-k*d,
        ];
        let (hk,hc) = (correction(axial[0],angular_stiffness),correction(axial[1],angular_damping));
        let mut links: Vec<_> = points.iter().map(|&(i,w,_)| Link {
            a:(body_index,i), b:anchor, k:axial[0]*w, c:axial[1]*w,
        }).collect();
        for (j,&(nj,wj,[xj,yj])) in points.iter().enumerate() {
            for &(ni,wi,[xi,yi]) in &points[..j] {
                let coefficient = |h: [f64;3]| -wi*wj*(xi*(h[0]*xj+h[1]*yj)+yi*(h[1]*xj+h[2]*yj));
                links.push(Link { a:(body_index,ni), b:Some((body_index,nj)), k:coefficient(hk), c:coefficient(hc) });
            }
        }
        if links.iter().any(|l| !l.k.is_finite() || !l.c.is_finite()) {
            return Err("affine mount: assembled coefficient overflow".into());
        }
        Ok(links)
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Law {
    /// f = K d^alpha (1 + mu d_dot): Hertz with Hunt-Crossley loss.
    Power { k: f64, alpha: f64, mu: f64 },
}

/// One-sided contact: d = sum_side sign * weight * w - gap; d > 0 is interpenetration.
#[derive(Clone)]
pub struct Contact {
    pub a: Vec<(usize, usize, f64)>, // (body, node, weight) with sign +1
    pub b: Vec<(usize, usize, f64)>, // sign -1
    pub gap: f64,
    pub law: Law,
    pub d_now: f64,
    pub active: bool,
    pub peak_force: f64,
    /// Force of the last step [N] (diagnostics).
    pub force: f64,
    /// The force solve stops at an exact Newton fixed point (see `contact_force`). On for the snare
    /// wires: held under load and warm-started, about one solve in five lands on its own root. Off
    /// elsewhere: it is as exact there, but it would change every other piece's output (by rounding
    /// for most; the open hi-hat's rim chatter decorrelates), and that has not been graded.
    pub fixed_point_exit: bool,
}

impl Contact {
    pub fn new(a: Vec<(usize, usize, f64)>, b: Vec<(usize, usize, f64)>, gap: f64, law: Law) -> Contact {
        Contact { a, b, gap, law, d_now: 0.0, active: true, peak_force: 0.0, force: 0.0, fixed_point_exit: false }
    }
}

/// A scheduled strike: at `time` the striker body (a point mass) is placed `clear` metres from
/// the struck surface and launched toward it at `speed`.
#[derive(Clone)]
pub struct Strike {
    pub time: f64,
    pub striker: usize,
    pub contact: usize,
    pub speed: f64,
    pub clear: f64,
    pub done: bool,
}

/// Piecewise-linear external force program on one node.
#[derive(Clone)]
pub struct Drive {
    pub body: usize,
    pub node: usize,
    pub knots: Vec<(f64, f64)>,
}

impl Drive {
    fn at(&self, t: f64) -> f64 {
        let k = &self.knots;
        if k.is_empty() || t <= k[0].0 {
            return k.first().map_or(0.0, |p| p.1);
        }
        for w in k.windows(2) {
            if t <= w[1].0 {
                let u = (t - w[0].0) / (w[1].0 - w[0].0);
                return w[0].1 + u * (w[1].1 - w[0].1);
            }
        }
        k.last().unwrap().1
    }
}

/// Rayleigh (baffled) pressure microphone.
#[derive(Clone)]
pub struct Mic {
    pub pos: [f64; 3],
    taps: Vec<(usize, usize, usize, f64, f64)>, // body, node, delay, weight lo, weight hi
    ptaps: Vec<(usize, usize, f64, f64)>,      // cavity port, delay, weight lo, weight hi
    ring: Vec<f64>,
    /// Unbaffled thin source of radius a: the far field tends to the baffled (Rayleigh) one for
    /// ka >> 1 and to a dipole (4/3pi) ka times it for ka << 1 (rigid disc). A first-order
    /// high-pass with tau = 4a/(3 pi c) joins the two limits. None = baffled.
    pub unbaffled: Option<f64>,
    /// A rigid floor at this height: surfaces facing away from the microphone still reach it by
    /// the floor reflection (image source; the snare's underside heard by an overhead).
    pub floor: Option<f64>,
    hp: (f64, f64),
}

impl Mic {
    pub fn new(pos: [f64; 3]) -> Mic {
        Mic { pos, taps: vec![], ptaps: vec![], ring: vec![], unbaffled: None, floor: None, hp: (0.0, 0.0) }
    }
}

// --------------------------------------------------------------------- contact solves
/// Secant of U(d)=d_+^(alpha+1)/(alpha+1), and its derivative in the first argument.
/// The Hertz factorization avoids subtracting nearly equal fractional powers.
#[inline(always)]
fn power_secant(d: f64, previous: f64, alpha: f64) -> (f64, f64) {
    if d <= 0.0 && previous <= 0.0 { return (0.0, 0.0); }
    if previous <= 0.0 {
        let den = d - previous;
        if alpha == 1.5 {
            let (root, ratio) = (d.sqrt(), d/den);
            return (0.4*d*root*ratio, 0.4*root*ratio*(2.5-ratio));
        }
        let p = d.powf(alpha);
        return (d*p/((alpha+1.0)*den), p*(alpha*d-(alpha+1.0)*previous)/((alpha+1.0)*den*den));
    }
    if d <= 0.0 {
        let den = previous - d;
        if alpha == 1.5 {
            let (root, ratio) = (previous.sqrt(), previous/den);
            return (0.4*previous*root*ratio, 0.4*root*ratio*ratio);
        }
        let potential = previous.powf(alpha+1.0)/(alpha+1.0);
        return (potential/den, potential/(den*den));
    }
    if alpha == 1.5 {
        let (a, b) = (d.sqrt(), previous.sqrt());
        let m = a.max(b);
        let (a, b) = (a/m, b/m);
        let (a2, b2, sum) = (a*a, b*b, a+b);
        let value = 0.4*m*m*m*(a2*a2+a2*a*b+a2*b2+a*b2*b+b2*b2)/sum;
        let derivative = 0.2*m*(3.0*a2*a+6.0*a2*b+4.0*a*b2+2.0*b2*b)/(sum*sum);
        return (value, derivative);
    }
    let z = (d-previous)/previous;
    let base = previous.powf(alpha);
    if z.abs()*alpha.abs().max(1.0) < 1e-3 {
        let a2 = alpha*(alpha-1.0);
        let a3 = a2*(alpha-2.0);
        let a4 = a3*(alpha-3.0);
        let a5 = a4*(alpha-4.0);
        let value = base*(1.0+z*(0.5*alpha+z*(a2/6.0+z*(a3/24.0+z*(a4/120.0+z*a5/720.0)))));
        let derivative = base/previous*(0.5*alpha+z*(a2/3.0+z*(a3/8.0+z*(a4/30.0+z*a5/144.0))));
        (value, derivative)
    } else {
        let delta = ((alpha+1.0)*z.ln_1p()).exp_m1();
        (base*delta/((alpha+1.0)*z),
         base/previous*((alpha+1.0)*z*(1.0+z).powf(alpha)-delta)/((alpha+1.0)*z*z))
    }
}

#[inline(always)]
fn contact_phi(d: f64, previous: f64, kk: f64, alpha: f64, mu: f64, k: f64) -> (f64, f64) {
    let (elastic, derivative) = power_secant(d, previous, alpha);
    let damping_rate = mu/(2.0*k);
    let damp = 1.0+damping_rate*(d-previous);
    let value = kk*elastic*damp;
    if value <= 0.0 { (0.0, 0.0) }
    else { (value, kk*(derivative*damp+elastic*damping_rate)) }
}

/// Centered discrete-gradient contact, phi = delta U / delta d times the
/// Hunt--Crossley factor at (d_next-d_previous)/(2k), clipped at zero force.
/// Thus zero-loss free masses conserve staggered kinetic plus contact energy.
/// The discrete-potential construction follows Bilbao et al.,
/// https://arxiv.org/abs/1405.2589, section 2.3, equation (11).
/// phi is increasing, so f=phi(d_pred-comp*f) lies in [0,phi(d_pred)]; unloading
/// can still exert force with negative predicted or final penetration.
/// Safeguarded Newton (a
/// bisection step whenever Newton leaves the bracket), from the last force when it lies in the bracket.
/// With `fixed_exit`, a Newton step that lands exactly on the iterate stops the solve there: the
/// correction is below half an ulp, so the iterate is the root to rounding. Without it the bracket
/// update has just put that iterate on the bracket's edge, so the safeguard bisects away from the
/// root and halves back down to it (30-47 iterations instead of about 3).
fn contact_force(d_pred: f64, comp: f64, force: f64, previous: f64, law: Law, fixed_exit: bool, k: f64) -> f64 {
    let Law::Power { k: kk, alpha, mu } = law;
    let phi = |d| contact_phi(d, previous, kk, alpha, mu, k);
    // no penetration, or a NaN one (a state already blown up): no force, and no 60 Newton steps on NaN
    if d_pred.is_nan() || previous.is_nan() || (d_pred <= 0.0 && previous <= 0.0) {
        return 0.0;
    }
    let (mut lo, mut hi) = (0.0, phi(d_pred).0);
    if hi <= 0.0 { return 0.0; }
    // a contact held under load barely changes between steps: start from the last force
    let mut x = if force > 0.0 && force < hi { force } else { 0.5 * hi };
    for _ in 0..60 {
        let (v, dv) = phi(d_pred - comp * x);
        let g = x - v;
        if g > 0.0 { hi = x; } else { lo = x; }
        let mut next = x - g / (1.0 + comp * dv);
        if !((next > lo && next < hi) || (fixed_exit && next == x)) {
            next = 0.5 * (lo + hi);
        }
        let done = (next - x).abs() <= 1e-13 * next.abs().max(1e-30) || hi - lo <= 1e-13 * hi.max(1e-30);
        x = next;
        if done {
            break;
        }
    }
    x
}

/// One iteration of `contact_force`'s loop for a = 1.5, without branches: the same operations in the
/// same order, so the same iterate, bracket and stop. Returns (next iterate, stop).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn newton_step_sqrt(x: f64, lo: &mut f64, hi: &mut f64, dp: f64, cm: f64, dn: f64, kk: f64, mu: f64, fixed_exit: bool, k: f64) -> (f64, bool) {
    let d = dp - cm * x;
    let (v, dv) = contact_phi(d, dn, kk, 1.5, mu, k);
    let g = x - v;
    let up = g > 0.0;
    *hi = if up { x } else { *hi };
    *lo = if up { *lo } else { x };
    let newton = x - g / (1.0 + cm * dv);
    let next = if (newton > *lo && newton < *hi) | (fixed_exit & (newton == x)) { newton } else { 0.5 * (*lo + *hi) };
    let done = ((next - x).abs() <= 1e-13 * next.abs().max(1e-30)) | (*hi - *lo <= 1e-13 * hi.max(1e-30));
    (next, done)
}

/// Contacts whose solves run interleaved. No two of a batch share a node, so each one's predicted
/// penetration, force and update are what solving them one at a time in list order gives.
#[derive(Clone, Default)]
struct ContactQueue {
    /// Per queued contact: its index, d_pred, comp, first iterate, upper bracket, d_previous, K, mu, fixed-point
    /// exit; the force.
    j: Vec<usize>,
    dp: Vec<f64>,
    cm: Vec<f64>,
    x0: Vec<f64>,
    hi: Vec<f64>,
    dn: Vec<f64>,
    kk: Vec<f64>,
    mu: Vec<f64>,
    fx: Vec<bool>,
    f: Vec<f64>,
    /// Per body and node, the batch that last claimed the node.
    mark: Vec<Vec<u32>>,
    batch: u32,
}

/// Solves in flight at once: independent Newton iterations overlap instead of waiting on each
/// other's square roots and divisions.
const LANES: usize = 4;

impl ContactQueue {
    fn new_batch(&mut self) {
        self.batch = self.batch.wrapping_add(1);
        if self.batch == 0 {
            for m in self.mark.iter_mut() {
                m.fill(0);
            }
            self.batch = 1;
        }
    }

    /// `contact_force` for every queued contact (all a = 1.5), LANES at a time; a lane that stops
    /// takes the next contact.
    fn solve(&mut self, k: f64) {
        // a batch with no fixed-point exit (every hi-hat batch) runs the loop without the per-lane flag
        if self.fx.contains(&true) { self.solve_lanes::<true>(k) } else { self.solve_lanes::<false>(k) }
    }

    fn solve_lanes<const FX: bool>(&mut self, k: f64) {
        let n = self.j.len();
        self.f.clear();
        self.f.resize(n, 0.0);
        let mut at = [n; LANES];
        let (mut x, mut lo, mut hi) = ([0.5; LANES], [0.0; LANES], [1.0; LANES]);
        let (mut dp, mut cm, mut dn, mut kk, mut mu) = ([1.0; LANES], [1.0; LANES], [0.0; LANES], [0.0; LANES], [0.0; LANES]);
        let mut fx = [false; LANES];
        let mut it = [0u32; LANES];
        let mut next = 0;
        let mut busy = 0;
        macro_rules! load {
            ($l:expr) => {{
                let (l, t) = ($l, next);
                at[l] = t;
                x[l] = self.x0[t];
                lo[l] = 0.0;
                hi[l] = self.hi[t];
                dp[l] = self.dp[t];
                cm[l] = self.cm[t];
                dn[l] = self.dn[t];
                kk[l] = self.kk[t];
                mu[l] = self.mu[t];
                if FX {
                    fx[l] = self.fx[t];
                }
                it[l] = 0;
                next += 1;
            }};
        }
        for l in 0..LANES {
            if next < n {
                load!(l);
                busy += 1;
            }
        }
        while busy > 0 {
            // an idle lane iterates on harmless values (a halving from 0.5) and is ignored
            let mut done = [false; LANES];
            for l in 0..LANES {
                let (nx, stop) = newton_step_sqrt(x[l], &mut lo[l], &mut hi[l], dp[l], cm[l], dn[l], kk[l], mu[l], FX & fx[l], k);
                x[l] = nx;
                done[l] = stop;
                it[l] += 1;
            }
            for l in 0..LANES {
                if at[l] < n && (done[l] || it[l] == 60) {
                    self.f[at[l]] = x[l];
                    if next < n {
                        load!(l);
                    } else {
                        at[l] = n;
                        (x[l], lo[l], hi[l]) = (0.5, 0.0, 1.0);
                        busy -= 1;
                    }
                }
            }
        }
    }

    /// Solve the batch, apply its forces in list order and start a new batch.
    fn flush(&mut self, bodies: &mut [Body], contacts: &mut [Contact], k: f64) {
        if !self.j.is_empty() {
            self.solve(k);
            for t in 0..self.j.len() {
                let c = &mut contacts[self.j[t]];
                let f = self.f[t];
                apply_force(bodies, c, f, self.dp[t], self.cm[t], k);
            }
            for v in [&mut self.dp, &mut self.cm, &mut self.x0, &mut self.hi, &mut self.dn, &mut self.kk, &mut self.mu, &mut self.f] {
                v.clear();
            }
            self.fx.clear();
            self.j.clear();
        }
        self.new_batch();
    }
}

/// Push the contact's nodes apart by force `f` and record it.
#[inline]
fn apply_force(bodies: &mut [Body], c: &mut Contact, f: f64, d_pred: f64, comp: f64, k: f64) {
    if f > 0.0 {
        for &(b, i, w) in &c.a {
            let cm = bodies[b].compliance(i, k);
            bodies[b].wn[i] -= w * cm * f;
        }
        for &(b, i, w) in &c.b {
            let cm = bodies[b].compliance(i, k);
            bodies[b].wn[i] += w * cm * f;
        }
    }
    c.peak_force = c.peak_force.max(f);
    c.force = f;
    c.d_now = d_pred - comp * f;
}

/// Every active contact in list order, each solved implicitly with the others' updates applied.
/// Runs of contacts on disjoint nodes with the sqrt law (snare wires, hi-hat rims, sticks) are
/// queued and solved interleaved; any other contact, or one that shares a node with the queue,
/// first flushes the queue. The result is bit-identical to solving one contact at a time.
fn solve_contacts(bodies: &mut [Body], contacts: &mut [Contact], q: &mut ContactQueue, k: f64) {
    if q.mark.len() != bodies.len() {
        q.mark.resize_with(bodies.len(), Vec::new);
    }
    for (m, b) in q.mark.iter_mut().zip(bodies.iter()) {
        if m.len() != b.n {
            m.clear();
            m.resize(b.n, 0);
        }
    }
    q.new_batch();
    for jc in 0..contacts.len() {
        if !contacts[jc].active {
            continue;
        }
        let Law::Power { k: kk, alpha, mu } = contacts[jc].law;
        let queued = alpha == 1.5;
        let c = &contacts[jc];
        if !queued || c.a.iter().chain(&c.b).any(|&(b, i, _)| q.mark[b][i] == q.batch) {
            q.flush(bodies, contacts, k);
        }
        let c = &mut contacts[jc];
        let mut d_pred = -c.gap;
        let mut previous = -c.gap;
        let mut comp = 0.0;
        for &(b, i, w) in &c.a {
            d_pred += w * bodies[b].wn[i];
            previous += w * bodies[b].wp[i];
            comp += w * w * bodies[b].compliance(i, k);
        }
        for &(b, i, w) in &c.b {
            d_pred -= w * bodies[b].wn[i];
            previous -= w * bodies[b].wp[i];
            comp += w * w * bodies[b].compliance(i, k);
        }
        // Read the actual previous pose, not a cache another contact or
        // friction solve may have invalidated. The secant spans two steps.
        if !queued || d_pred.is_nan() || previous.is_nan() || (d_pred <= 0.0 && previous <= 0.0) {
            // solved now: the queue is empty, or this contact moves nothing
            let f = contact_force(d_pred, comp, c.force, previous, c.law, c.fixed_point_exit, k);
            apply_force(bodies, c, f, d_pred, comp, k);
            continue;
        }
        for &(b, i, _) in c.a.iter().chain(&c.b) {
            q.mark[b][i] = q.batch;
        }
        let hi = contact_phi(d_pred, previous, kk, alpha, mu, k).0;
        q.j.push(jc);
        q.dp.push(d_pred);
        q.cm.push(comp);
        q.x0.push(if c.force > 0.0 && c.force < hi { c.force } else { 0.5 * hi });
        q.hi.push(hi);
        q.dn.push(previous);
        q.kk.push(kk);
        q.mu.push(mu);
        q.fx.push(c.fixed_point_exit);
    }
    q.flush(bodies, contacts, k);
}

/// Dry (Coulomb) friction between two plate surfaces that a contact presses together (a closed
/// hi-hat's rims). Bending moves a plate's faces sideways by half its thickness times its slope
/// (Kirchhoff: u = -z grad w), so the two faces of a contact slip against each other by
/// s = (t_a / 2) grad w_a + (t_b / 2) grad w_b, the face below a plate moving opposite to the face
/// above it. `sx` and `sy` hold the two slip components as weighted sums of node displacements. The
/// friction force acts against the final slip with at most `mu` times the normal force. It is
/// solved after normal contact using backward slip. This reduces free-mass kinetic energy;
/// it does not certify passivity of the complete coupled oscillator discretization.
#[derive(Clone)]
pub struct Friction {
    pub contact: usize,
    pub sx: Vec<(usize, usize, f64)>,
    pub sy: Vec<(usize, usize, f64)>,
    pub mu: f64,
    /// Friction force of the last step [N] (diagnostics), and whether the faces stuck.
    pub force: f64,
    pub stuck: bool,
}

impl Friction {
    /// The slip of the faces between node `ia` of body `a` (its upper face) and node `ib` of body `b`
    /// (its lower face), with plate thicknesses `ta` and `tb` there.
    pub fn between(bodies: &[Body], contact: usize, a: usize, ia: usize, ta: f64, b: usize, ib: usize, tb: f64, mu: f64) -> Friction {
        let [ax, ay] = bodies[a].slope(ia);
        let [bx, by] = bodies[b].slope(ib);
        // a face at height z above its plate's midplane moves by -z grad w: the upper face of the lower
        // plate a (z = ta/2) by -(ta/2) grad w_a, the lower face of the upper plate b (z = -tb/2) by
        // +(tb/2) grad w_b, so the slip of b's face over a's is (tb/2) grad w_b + (ta/2) grad w_a
        let side = |terms: Vec<(usize, f64)>, body: usize, t: f64| terms.into_iter().map(move |(i, w)| (body, i, 0.5 * t * w));
        let sx = side(ax, a, ta).chain(side(bx, b, tb)).collect();
        let sy = side(ay, a, ta).chain(side(by, b, tb)).collect();
        Friction { contact, sx, sy, mu, force: 0.0, stuck: false }
    }
}

/// Minimize d'f + f'Cf/2 over the Coulomb disk |f| <= cap, for SPD C.
/// Sliding satisfies (C + lambda I) f = -d, so f opposes the final slip
/// d + C f = -lambda f. Radially clipping C^-1 d fails for anisotropic C.
fn friction_force(dx: f64, dy: f64, cxx: f64, cyy: f64, cxy: f64, cap: f64) -> (f64, f64, bool) {
    if cap == 0.0 { return (0.0, 0.0, false); }
    let det = cxx * cyy - cxy * cxy;
    let (mut fx, mut fy) = (-(cyy * dx - cxy * dy) / det, -(cxx * dy - cxy * dx) / det);
    let mag = (fx * fx + fy * fy).sqrt();
    let stuck = mag <= cap;
    if !stuck {
        let (mut lo, mut hi) = (0.0, (dx * dx + dy * dy).sqrt() / cap);
        // Exact for isotropic C, and near zero at the stick/slide boundary.
        // An overflowed unconstrained norm must not seed lambda with NaN.
        let mut lambda = if mag.is_finite() { hi * ((mag - cap) / mag) } else { 0.5 * hi };
        for _ in 0..64 {
            let (a, b) = (cxx + lambda, cyy + lambda);
            let det = a * b - cxy * cxy;
            fx = -(b * dx - cxy * dy) / det;
            fy = -(a * dy - cxy * dx) / det;
            let mag = (fx * fx + fy * fy).sqrt();
            let residual = mag - cap;
            if residual.abs() <= 1e-12 * cap { break; }
            if residual > 0.0 { lo = lambda; } else { hi = lambda; }
            let derivative = -(b * fx * fx - 2.0 * cxy * fx * fy + a * fy * fy) / (det * mag);
            let next = lambda - residual / derivative;
            lambda = if next > lo && next < hi { next } else { 0.5 * (lo + hi) };
        }
    }
    (fx, fy, stuck)
}

fn solve_frictions(bodies: &mut [Body], contacts: &[Contact], frictions: &mut [Friction], k: f64) {
    for fr in frictions.iter_mut() {
        let c = &contacts[fr.contact];
        let n = if c.active { c.force } else { 0.0 };
        if !(n > 0.0) || fr.mu <= 0.0 {
            fr.force = 0.0;
            fr.stuck = false;
            continue;
        }
        // slip over this step as predicted (contacts included), and the slip compliance matrix
        let inc = |t: &[(usize, usize, f64)]| t.iter().map(|&(b, i, w)| w * (bodies[b].wn[i] - bodies[b].w[i])).sum::<f64>();
        let (dx, dy) = (inc(&fr.sx), inc(&fr.sy));
        let (mut cxx, mut cyy, mut cxy) = (0.0, 0.0, 0.0);
        for &(b, i, w) in &fr.sx {
            let cm = bodies[b].compliance(i, k);
            cxx += w * w * cm;
            for &(b2, i2, w2) in &fr.sy {
                if b2 == b && i2 == i {
                    cxy += w * w2 * cm;
                }
            }
        }
        for &(b, i, w) in &fr.sy {
            cyy += w * w * bodies[b].compliance(i, k);
        }
        let det = cxx * cyy - cxy * cxy;
        if !(det > 0.0) {
            fr.force = 0.0;
            fr.stuck = false;
            continue;
        }
        let (fx, fy, stuck) = friction_force(dx, dy, cxx, cyy, cxy, fr.mu * n);
        fr.stuck = stuck;
        for &(b, i, w) in &fr.sx {
            let cm = bodies[b].compliance(i, k);
            bodies[b].wn[i] += cm * w * fx;
        }
        for &(b, i, w) in &fr.sy {
            let cm = bodies[b].compliance(i, k);
            bodies[b].wn[i] += cm * w * fy;
        }
        fr.force = (fx * fx + fy * fy).sqrt();
    }
}

// --------------------------------------------------------------------- instrument
#[derive(Clone)]
pub struct Instrument {
    pub bodies: Vec<Body>,
    pub cavities: Vec<Cavity>,
    pub rads: Vec<RadDamp>,
    pub links: Vec<Link>,
    pub contacts: Vec<Contact>,
    /// Dry friction on the faces some contacts press together, solved after the contacts.
    pub frictions: Vec<Friction>,
    pub strikes: Vec<Strike>,
    pub drives: Vec<Drive>,
    pub mics: Vec<Mic>,
    pub rate: f64,
    pub step_index: u64,
    // One admitted stationary medium for cavities, ports, exterior loading and
    // microphones. Changing it during a trajectory requires physical state and
    // parameter-work transfer, so callers build a new instrument instead.
    gas: GasState,
    materials: Option<std::sync::Arc<materials::KitMaterials>>,
    prepared: bool,
    cq: ContactQueue,
}

impl Instrument {
    pub fn new(rate: f64) -> Instrument {
        Self::with_gas(rate, standard_gas()).expect("instrument sample rate must be positive and finite")
    }
    /// Assemble in a coherent material-resolved, spatially uniform atmosphere.
    /// Temperature, pressure and composition are evaluated by fs-material;
    /// independent density/sound-speed edits are refused before assembly.
    pub fn with_gas(rate: f64, gas: GasState) -> Result<Instrument, String> {
        gas.validate_coherence().map_err(|e| e.to_string())?;
        if !rate.is_finite() || rate <= 0.0 {
            return Err("instrument sample rate must be positive and finite".into());
        }
        Ok(Instrument { bodies: vec![], cavities: vec![], rads: vec![], links: vec![], contacts: vec![], frictions: vec![], strikes: vec![], drives: vec![], mics: vec![], rate, step_index: 0, gas, materials: None, prepared: false, cq: ContactQueue::default() })
    }
    /// A new trajectory retaining the same elastic sources its kit builder uses.
    pub(crate) fn with_materials(rate: f64, gas: GasState, materials: materials::KitMaterials) -> Result<Instrument, String> {
        let mut instrument = Self::with_gas(rate, gas)?;
        instrument.materials = Some(std::sync::Arc::new(materials));
        Ok(instrument)
    }
    pub fn materials(&self) -> Option<&materials::KitMaterials> { self.materials.as_deref() }
    /// The same material snapshot used by every gas-dependent operator.
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    pub fn add(&mut self, b: Body) -> usize {
        self.bodies.push(b);
        self.bodies.len() - 1
    }
    pub fn k(&self) -> f64 {
        1.0 / self.rate
    }
    pub fn time(&self) -> f64 {
        self.step_index as f64 / self.rate
    }

    /// Largest lambda_max(M^-1 K) over the bodies, and which body: leapfrog needs k^2 lambda < 4.
    pub fn max_rate2(&mut self) -> (f64, String) {
        let mut worst = (0.0, String::new());
        for b in self.bodies.iter_mut() {
            if b.n < 2 {
                continue;
            }
            let lam = b.max_rate2();
            if lam > worst.0 {
                worst = (lam, b.name.clone());
            }
        }
        worst
    }
    /// Rest-state timestep estimate, including shell curvature, without altering live state.
    /// Power iteration is an estimate; large-deflection/contact convergence is separate.
    pub fn admission_rate2(&self) -> (f64, String) {
        let mut worst = (0.0, String::new());
        for b in &self.bodies {
            if b.n < 2 { continue; }
            // The existing diagnostic mutates force caches; use a cold image.
            let mut image = b.clone();
            let linear = image.max_rate2();
            let lam = if image.vk.is_some() { linear.max(image.max_rate2_with_shell()) } else { linear };
            if lam > worst.0 { worst = (lam, b.name.clone()); }
        }
        let gas_rate2 = self.gas_admission_rate2(self.rate);
        if gas_rate2 > 0.0 {
            worst.0 += gas_rate2;
            worst.1.push_str(" + gas coupling");
        }
        let viscous_rate2 = self.viscosity_admission_rate2(self.rate);
        if viscous_rate2 > 0.0 {
            worst.0 += viscous_rate2;
            worst.1.push_str(" + material bending viscosity");
        }
        let radiation_rate2 = 2.0 * self.rate * self.bodies.iter().filter_map(|b| b.head_radiation.as_deref())
            .map(|load| load.rate_bound).fold(0.0_f64, f64::max);
        if radiation_rate2 > 0.0 {
            worst.0 += radiation_rate2;
            worst.1.push_str(" + head radiation");
        }
        worst
    }

    /// Factor every head-radiation air-mass system at this instrument's own
    /// step, refusing at build time what stepping could not use.
    pub fn validate_head_radiation(&mut self) -> Result<(), String> {
        let k = self.k();
        for body in &mut self.bodies {
            body.prepare_head_radiation(k)?;
        }
        Ok(())
    }

    fn viscosity_admission_rate2(&self, rate: f64) -> f64 {
        2.0 * rate * self.bodies.iter().filter_map(|b| b.viscosity.as_ref())
            .map(|v| v.rate_bound).fold(0.0_f64, f64::max)
    }

    /// Source-domain admission at the actual internal integration clock.
    /// A prescribed broadband time-domain law must cover DC through Nyquist;
    /// a source query inside a narrow band is not sufficient.
    pub fn validate_bending_loss_clock(&self) -> Result<(), String> {
        let omega = PI * self.rate;
        for body in &self.bodies {
            if let Some(viscous) = &body.viscosity {
                let (lo, hi) = viscous.source.omega_band_rad_s();
                if !omega.is_finite() || self.rate <= 0.0 || lo > 0.0 || hi < omega {
                    return Err(format!("{}: bending-viscosity source must cover 0 through internal Nyquist {} Hz; source band [{lo}, {hi}] rad/s", body.name, self.rate / 2.0));
                }
            }
        }
        if let Some(materials) = self.materials() { materials.validate_bending_loss_clock(self.rate)?; }
        Ok(())
    }

    // Separate the rate-dependent gas bound from the expensive structural
    // eigensolve so clock selection can evaluate every candidate correctly.
    fn gas_admission_rate2(&self, rate: f64) -> f64 {
        // Each cavity is a positive rank-one stiffness in volume displacement.
        // Its sole mass-normalized eigenvalue is K_V * sum(u_i^2 / m_i),
        // including the port coordinate's inverse inertance. Adding these
        // bounds to the structural estimate prevents a light/hot or compressed
        // gas from silently bypassing the timestep admission. Ignoring added
        // gas mass is conservative. Explicit radiation adds 2 D / dt to the
        // leapfrog stability inequality dt^2 K + 2 dt D < 4 M.
        let inverse_volume_mass = |heads: &[usize]| -> f64 {
            heads.iter().map(|&index| {
                let body = &self.bodies[index];
                body.h.powi(4) * body.mass.iter().map(|m| 1.0 / m).sum::<f64>()
            }).sum()
        };
        let mut gas_rate2 = 0.0;
        for cavity in &self.cavities {
            let inverse_port_mass = cavity.port.as_ref()
                .map_or(0.0, |port| port.area / (self.gas.density * port.length));
            gas_rate2 += self.gas.density * self.gas.sound_speed.powi(2) / cavity.volume
                * (inverse_volume_mass(&cavity.heads) + inverse_port_mass);
            if let Some(port) = &cavity.port {
                gas_rate2 += 2.0 * rate * port.resistance * inverse_port_mass;
            }
        }
        for radiation in &self.rads {
            gas_rate2 += 2.0 * rate * radiation.resistance
                * inverse_volume_mass(&radiation.bodies);
        }
        gas_rate2
    }

    /// (k^2 lambda_K + 2 k lambda_C) / 4 including rest curvature and explicit
    /// gas/material damping; excludes contact and large-deflection changes.
    pub fn stability(&mut self) -> (f64, String) {
        let k = self.k();
        let (lam, who) = self.admission_rate2();
        (k * k * lam / 4.0, who)
    }

    fn prepare(&mut self) {
        let k = self.k();
        for m in self.mics.iter_mut() {
            m.taps.clear();
            let mut dmax = 0;
            for (bi, b) in self.bodies.iter().enumerate() {
                if !b.radiates {
                    continue;
                }
                let h2 = b.h * b.h;
                // the direct path, and the floor's image of the microphone when there is a floor
                let mut targets = vec![m.pos];
                if let Some(z) = m.floor {
                    targets.push([m.pos[0], m.pos[1], 2.0 * z - m.pos[2]]);
                }
                for tp in targets {
                    for i in 0..b.n {
                        let p = [b.origin[0] + b.xy[i][0], b.origin[1] + b.xy[i][1], b.origin[2] + b.facing * b.rest_z[i]];
                        let d = [tp[0] - p[0], tp[1] - p[1], tp[2] - p[2]];
                        if d[2] * b.facing <= 0.0 {
                            continue; // the other side of this surface faces away from this path
                        }
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        let delay = r / self.gas.sound_speed * self.rate;
                        let di = delay.floor() as usize;
                        let fr = delay - di as f64;
                        // w is the OUTWARD displacement of the facing side, so no sign flip
                        let wgt = self.gas.density * h2 / (2.0 * PI * r) / (k * k);
                        m.taps.push((bi, i, di, wgt * (1.0 - fr), wgt * fr));
                        dmax = dmax.max(di + 2);
                    }
                }
            }
            m.ptaps.clear();
            for (ci, c) in self.cavities.iter().enumerate() {
                if let Some(q) = &c.port {
                    let d = [m.pos[0] - q.pos[0], m.pos[1] - q.pos[1], m.pos[2] - q.pos[2]];
                    if d[2] * q.facing <= 0.0 {
                        continue;
                    }
                    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    let delay = r / self.gas.sound_speed * self.rate;
                    let di = delay.floor() as usize;
                    let fr = delay - di as f64;
                    let wgt = self.gas.density / (2.0 * PI * r);
                    m.ptaps.push((ci, di, wgt * (1.0 - fr), wgt * fr));
                    dmax = dmax.max(di + 2);
                }
            }
            m.ring = vec![0.0; (dmax + 2).next_power_of_two()];
        }
        self.prepared = true;
    }

    /// One time step; returns each microphone's pressure sample [Pa].
    pub fn step(&mut self, out: &mut [f64]) {
        if !self.prepared {
            self.prepare();
        }
        let k = self.k();
        let t = self.time();
        for b in self.bodies.iter_mut() {
            b.linear_forces();
        }
        for l in &self.links {
            let (ba, na) = l.a;
            let (wa, va) = (self.bodies[ba].w[na], (self.bodies[ba].w[na] - self.bodies[ba].wp[na]) / k);
            let (wb, vb) = match l.b {
                Some((bb, nb)) => (self.bodies[bb].w[nb], (self.bodies[bb].w[nb] - self.bodies[bb].wp[nb]) / k),
                None => (0.0, 0.0),
            };
            let f = -l.k * (wa - wb) - l.c * (va - vb);
            self.bodies[ba].f[na] += f;
            if let Some((bb, nb)) = l.b {
                self.bodies[bb].f[nb] -= f;
            }
        }
        for c in self.cavities.iter_mut() {
            let mut dv = 0.0;
            for &hb in &c.heads {
                let b = &self.bodies[hb];
                dv += b.w.iter().sum::<f64>() * b.h * b.h;
            }
            let p = -self.gas.density * self.gas.sound_speed.powi(2) / c.volume * (dv + c.port.as_ref().map_or(0.0, |q| q.x));
            for &hb in &c.heads {
                let b = &mut self.bodies[hb];
                let fa = p * b.h * b.h;
                for f in b.f.iter_mut() {
                    *f += fa;
                }
            }
            if let Some(q) = c.port.as_mut() {
                // the air plug: rho L / A du/dt = p - R u (symplectic Euler, far below its limit)
                let du = k * (p - q.resistance * q.u) * q.area / (self.gas.density * q.length);
                q.u += du;
                q.du = du / k;
                q.x += k * q.u;
            }
        }
        for r in &self.rads {
            let mut q = 0.0;
            for &bi in &r.bodies {
                let b = &self.bodies[bi];
                q += b.w.iter().zip(&b.wp).map(|(a, c)| a - c).sum::<f64>() / k * b.h * b.h;
            }
            for &bi in &r.bodies {
                let b = &mut self.bodies[bi];
                let fa = -r.resistance * q * b.h * b.h;
                for f in b.f.iter_mut() {
                    *f += fa;
                }
            }
        }
        for d in &self.drives {
            let v = d.at(t);
            self.bodies[d.body].f[d.node] += v;
        }
        for b in self.bodies.iter_mut() {
            b.predict(k);
        }
        // launches
        for s in self.strikes.iter_mut() {
            if s.done || t < s.time {
                continue;
            }
            s.done = true;
            // one stick, one contact: a new launch retires the striker's earlier contacts (a hi-hat or
            // cymbal adds one contact per hit, and a stale one would push on the shell again)
            for (j, c) in self.contacts.iter_mut().enumerate() {
                if j != s.contact && c.a.iter().chain(&c.b).any(|&(b, _, _)| b == s.striker) {
                    c.active = false;
                }
            }
            let c = &mut self.contacts[s.contact];
            let surf: f64 = c.a.iter().map(|&(b, i, w)| w * self.bodies[b].w[i]).sum();
            // striker on side b: d = surf - z - gap, so z starts at surf - gap - clear and moves "down" toward -d
            let z0 = surf - c.gap + s.clear;
            let st = &mut self.bodies[s.striker];
            st.w[0] = z0;
            st.wp[0] = z0 + s.speed * k;
            st.wn[0] = z0 - s.speed * k;
            c.active = true;
            c.d_now = surf - z0 - c.gap;
        }
        // contacts, solved implicitly one at a time (see `solve_contacts`)
        solve_contacts(&mut self.bodies, &mut self.contacts, &mut self.cq, k);
        if !self.frictions.is_empty() {
            solve_frictions(&mut self.bodies, &self.contacts, &mut self.frictions, k);
        }
        for b in self.bodies.iter_mut() {
            b.finish(k);
        }
        // microphones: centred acceleration at this step
        for (mi, m) in self.mics.iter_mut().enumerate() {
            let mask = m.ring.len() - 1;
            let base = self.step_index as usize;
            for &(bi, i, di, lo, hi) in &m.taps {
                let b = &self.bodies[bi];
                let a = b.wn[i] - 2.0 * b.w[i] + b.wp[i];
                m.ring[(base + di) & mask] += lo * a;
                m.ring[(base + di + 1) & mask] += hi * a;
            }
            for &(ci, di, lo, hi) in &m.ptaps {
                let a = self.cavities[ci].port.as_ref().map_or(0.0, |q| q.du);
                m.ring[(base + di) & mask] += lo * a;
                m.ring[(base + di + 1) & mask] += hi * a;
            }
            let slot = base & mask;
            let mut y = m.ring[slot];
            if let Some(a) = m.unbaffled {
                // bilinear first-order high-pass, H(s) = s tau / (1 + s tau)
                let tau = 4.0 * a / (3.0 * PI * self.gas.sound_speed);
                let g = 2.0 * tau * self.rate;
                let (x1, y1) = m.hp;
                let yn = (g * (y - x1) - (1.0 - g) * y1) / (1.0 + g);
                m.hp = (y, yn);
                y = yn;
            }
            out[mi] = y;
            m.ring[slot] = 0.0;
        }
        for b in self.bodies.iter_mut() {
            b.rotate();
        }
        self.step_index += 1;
    }

    /// Everything one time step reads, as plain data: for a stepper outside this crate, the CUDA renderer
    /// in apps/instrument-lab/render-gpu. Exporting changes nothing the steps compute, and `step` stays the
    /// reference that stepper is checked against. It prepares the microphones' taps as the first step would,
    /// so place microphones and bodies before exporting (a later move would no longer re-prepare them).
    pub fn image(&mut self) -> InstrumentImage {
        if !self.prepared {
            self.prepare();
        }
        let k = self.k();
        let bodies = self.bodies.iter_mut().map(|b| b.image(k)).collect();
        let mics = self
            .mics
            .iter()
            .map(|m| MicImage { taps: m.taps.clone(), ptaps: m.ptaps.clone(), ring_len: m.ring.len(), unbaffled: m.unbaffled, hp: m.hp })
            .collect();
        InstrumentImage {
            rate: self.rate,
            gas: self.gas,
            materials: self.materials.clone(),
            step_index: self.step_index,
            bodies,
            cavities: self.cavities.clone(),
            rads: self.rads.clone(),
            links: self.links.clone(),
            contacts: self.contacts.clone(),
            frictions: self.frictions.clone(),
            strikes: self.strikes.clone(),
            drives: self.drives.clone(),
            mics,
        }
    }
}

// --------------------------------------------------------------------- export
/// A sparse matrix as compressed rows: row `i` holds `val[ptr[i]..ptr[i + 1]]` at columns
/// `col[ptr[i]..ptr[i + 1]]`, summed in that order.
#[derive(Clone, Debug)]
pub struct CsrImage {
    pub ptr: Vec<usize>,
    pub col: Vec<usize>,
    pub val: Vec<f64>,
}

/// The von Karman data of a plate: per core node its 9-point stencil in the body's node numbering
/// (centre, E, W, N, S, NE, NW, SE, SW; None off the core), l(z0) there, the rest height per node,
/// E t, and the stress function of the last step (the core numbering of `Body::airy_solve`).
#[derive(Clone, Debug)]
pub struct VonKarmanImage {
    pub eh: f64,
    pub nb: Vec<[Option<usize>; 9]>,
    pub l0: Vec<f64>,
    pub z0: Vec<f64>,
    pub f: Vec<f64>,
}

/// One body's operators, per-node update coefficients and state. `cf[i]` is
/// [k^2 / a, 2 m / k^2, (1 - s0 k) m / k^2, 2 s1 m / k] with a = m (1 + s0 k) / k^2, at the
/// instrument's step k; the first is also node i's compliance (displacement per unit force in a step).
#[derive(Clone, Debug)]
pub struct BodyImage {
    pub name: String,
    pub n: usize,
    pub h: f64,
    pub mass: Vec<f64>,
    pub s0: Vec<f64>,
    pub s1: f64,
    pub cf: Vec<[f64; 4]>,
    pub lap: CsrImage,
    pub stiff: CsrImage,
    /// Source-retaining curvature-rate law, its C=G'G operator and three G rows
    /// per node. The separate joule counter is physical centered-velocity loss.
    pub viscosity: Option<PlaneStressViscosity>,
    pub viscosity_angle_rad: f64,
    pub viscous: CsrImage,
    pub viscous_root: CsrImage,
    pub viscous_dissipated_j: f64,
    pub tension: f64,
    pub berger: f64,
    pub bias: Vec<f64>,
    pub vk: Option<VonKarmanImage>,
    pub area: f64,
    pub air_mass: f64,
    pub w: Vec<f64>,
    pub wp: Vec<f64>,
    /// L w of the previous step (`lap` times `wp`).
    pub lwp: Vec<f64>,
    pub psi: f64,
    pub r_prev: f64,
    pub budget: f64,
}

/// A microphone's prepared taps: (body, node, whole-sample delay, weight on that sample, weight on
/// the next) on surface acceleration, and (cavity, delay, weights) on port flow acceleration.
#[derive(Clone, Debug)]
pub struct MicImage {
    pub taps: Vec<(usize, usize, usize, f64, f64)>,
    pub ptaps: Vec<(usize, usize, f64, f64)>,
    pub ring_len: usize,
    pub unbaffled: Option<f64>,
    pub hp: (f64, f64),
}

#[derive(Clone)]
pub struct InstrumentImage {
    pub rate: f64,
    /// Retained coherent medium for external steppers and observation filters.
    pub gas: GasState,
    /// Exact elastic sources used to build the shared operators, when supplied.
    pub materials: Option<std::sync::Arc<materials::KitMaterials>>,
    pub step_index: u64,
    pub bodies: Vec<BodyImage>,
    pub cavities: Vec<Cavity>,
    pub rads: Vec<RadDamp>,
    pub links: Vec<Link>,
    pub contacts: Vec<Contact>,
    pub frictions: Vec<Friction>,
    pub strikes: Vec<Strike>,
    pub drives: Vec<Drive>,
    pub mics: Vec<MicImage>,
}

impl Csr {
    fn image(&self) -> CsrImage {
        // a matrix with no rows has no ptr[0]
        let ptr = if self.ptr.is_empty() { vec![0; self.n + 1] } else { self.ptr.clone() };
        CsrImage { ptr, col: self.col.clone(), val: self.val.clone() }
    }
}

impl Body {
    fn image(&mut self, k: f64) -> BodyImage {
        self.coefs(k);
        BodyImage {
            name: self.name.clone(),
            n: self.n,
            h: self.h,
            mass: self.mass.clone(),
            s0: self.s0.clone(),
            s1: self.s1,
            cf: self.cf.clone(),
            lap: self.lap.image(),
            stiff: self.stiff.image(),
            viscosity: self.viscosity.as_ref().map(|v| v.source.clone()),
            viscosity_angle_rad: self.viscosity.as_ref().map_or(0.0, |v| v.angle),
            viscous: self.viscosity.as_ref().map_or_else(|| Csr::default().image(), |v| v.c.image()),
            viscous_root: self.viscosity.as_ref().map_or_else(|| Csr::default().image(), |v| v.root.image()),
            viscous_dissipated_j: self.bending_loss_dissipated_j(),
            tension: self.tension,
            berger: self.berger,
            bias: self.bias.clone(),
            vk: self.vk.as_ref().map(|v| VonKarmanImage { eh: v.eh, nb: v.nb.clone(), l0: v.l0.clone(), z0: v.z0.clone(), f: v.f.clone() }),
            area: self.area,
            air_mass: self.air_mass,
            w: self.w.clone(),
            wp: self.wp.clone(),
            lwp: self.lwp.clone(),
            psi: self.psi,
            r_prev: self.r_prev,
            budget: self.budget,
        }
    }

    /// Solve the plate's Airy system A F = x in place (core numbering), with the factor every step
    /// uses. False (and `x` untouched) for a body without von Karman coupling. Column j of A^-1 is
    /// the solve of the j-th unit vector.
    pub fn airy_solve(&self, x: &mut [f64]) -> bool {
        match &self.vk {
            Some(v) if x.len() == v.nb.len() => {
                v.airy.solve(x);
                true
            }
            _ => false,
        }
    }
}

// --------------------------------------------------------------------- decimation
fn decimation_kernel(factor: usize) -> Vec<f64> {
    let taps = 32 * factor + 1;
    let fc = 0.46 / factor as f64; // cycles per input sample: 0.92 x output Nyquist
    let mid = (taps / 2) as f64;
    (0..taps)
        .map(|i| {
            let t = i as f64 - mid;
            let s = if t == 0.0 { 2.0 * fc } else { (2.0 * PI * fc * t).sin() / (PI * t) };
            let wdw = 0.42 - 0.5 * (2.0 * PI * i as f64 / (taps - 1) as f64).cos() + 0.08 * (4.0 * PI * i as f64 / (taps - 1) as f64).cos();
            s * wdw
        })
        .collect()
}

/// The offline filter in causal streaming form: 16 output frames of delay.
/// Phase and history survive arbitrary render blocks and cloned instruments.
#[derive(Clone)]
pub(crate) struct StreamingDecimator {
    kernel: Vec<f64>,
    history: Vec<f64>,
    next: usize,
}
impl StreamingDecimator {
    fn new(factor: usize) -> Self {
        let kernel = if factor == 1 { vec![1.0] } else { decimation_kernel(factor) };
        Self { history: vec![0.0; kernel.len()], kernel, next: 0 }
    }
    fn push(&mut self, value: f64) {
        self.history[self.next] = value;
        self.next = (self.next + 1) % self.history.len();
    }
    fn value(&self) -> f64 {
        let (a,b) = self.history.split_at(self.next);
        let mut value = 0.0;
        for (h,x) in self.kernel.iter().zip(b.iter().chain(a)) { value += h*x; }
        value
    }
    fn clear(&mut self) { self.history.fill(0.0); self.next = 0; }
}

/// Windowed-sinc low-pass then keep every `factor`-th sample.
pub fn decimate(x: &[f64], factor: usize) -> Vec<f64> {
    if factor == 1 { return x.to_vec(); }
    let h = decimation_kernel(factor);
    let taps = h.len();
    let half = taps / 2;
    (0..x.len() / factor)
        .map(|o| {
            let c = o * factor;
            let mut acc = 0.0;
            for (j, hj) in h.iter().enumerate() {
                let idx = c as i64 + j as i64 - half as i64;
                if idx >= 0 && (idx as usize) < x.len() {
                    acc += hj * x[idx as usize];
                }
            }
            acc
        })
        .collect()
}

/// Float32 WAV bytes, interleaved channels.
pub fn wav_f32(channels: &[Vec<f64>], rate: u32) -> Vec<u8> {
    let nch = channels.len();
    let frames = channels.iter().map(|c| c.len()).min().unwrap_or(0);
    let data_len = (frames * nch * 4) as u32;
    let mut v = Vec::with_capacity(44 + data_len as usize);
    v.extend_from_slice(b"RIFF");
    v.extend_from_slice(&(36 + data_len).to_le_bytes());
    v.extend_from_slice(b"WAVEfmt ");
    v.extend_from_slice(&16u32.to_le_bytes());
    v.extend_from_slice(&3u16.to_le_bytes());
    v.extend_from_slice(&(nch as u16).to_le_bytes());
    v.extend_from_slice(&rate.to_le_bytes());
    v.extend_from_slice(&(rate * 4 * nch as u32).to_le_bytes());
    v.extend_from_slice(&((4 * nch) as u16).to_le_bytes());
    v.extend_from_slice(&32u16.to_le_bytes());
    v.extend_from_slice(b"data");
    v.extend_from_slice(&data_len.to_le_bytes());
    for f in 0..frames {
        for c in channels {
            v.extend_from_slice(&(c[f] as f32).to_le_bytes());
        }
    }
    v
}

pub mod kit;
pub mod live;

#[cfg(test)]
mod tests;
