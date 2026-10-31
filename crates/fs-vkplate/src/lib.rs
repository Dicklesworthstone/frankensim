//! # fs-vkplate: struck thin plates and membranes in the time domain
//!
//! One explicit finite-difference solver for every vibrating surface of a drum kit.
//! A body is a thin plate on a square grid (spacing `h`) with
//!
//! `rho t w_tt = T Lap(w) - div div (D curv) + vK(w, F) - 2 rho t s0 w_t + 2 s1 div(rho t grad(w_t)) + p + f_contact`
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

mod annular_membrane;
pub use annular_membrane::MembraneTerm;

pub use fs_material::gas::GasState;
pub use fs_material::plane_stress_viscosity::PlaneStressViscosity;

/// Immutable source-bound elastic data and physical striker geometry for kit builders.
pub mod materials;

/// Baked exterior modal loading of drum heads' motion beyond their compact volume model.
pub mod radiation;

/// Air loss laws read from the gas record: oscillating-cylinder drag and the
/// viscous layer on a plate's faces.
pub mod air_loss;

/// Coiled steel snare wires from a declared geometry, two material cards and the gas.
pub mod coil;

/// Elliptical Hertz contact of two elastic bodies.
pub mod hertz;

/// A head film's generalized-Maxwell law on a plate's bending operator, fitted to the film card's
/// loss-tangent curve (every grid mode carries the film loss).
pub mod film;

/// A cymbal struck with the stick's shaft across its edge: the shaft's effective mass and its
/// elastic contact across the grain from the stick card.
pub mod shaft;

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

#[cfg(test)]
mod air_loss_tests;

#[cfg(test)]
mod bending_quadrature_tests;

#[cfg(test)]
mod film_tests;

/// Four quarter-cell quadrature points, each with three tensor-root rows.
/// Shared with consumers of `BodyImage::viscous_root`.
pub const BENDING_ROOT_ROWS_PER_NODE: usize = 12;

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

    /// Maximal runs of rows sharing one stencil shape with 3, 5, 13 or 17 entries (a string, the
    /// Laplacian, the bending operator).
    fn find_runs(&self) -> Vec<Run> {
        let (ptr, col) = (&self.ptr, &self.col);
        let off = |i: usize| col[ptr[i]..ptr[i + 1]].iter().map(move |&j| j as isize - i as isize);
        let mut runs = vec![];
        let mut i = 0;
        while i < self.n {
            let k = ptr[i + 1] - ptr[i];
            let mut e = i + 1;
            if matches!(k, 3 | 5 | 13 | 17) {
                while e < self.n && ptr[e + 1] - ptr[e] == k && off(e).eq(off(i)) {
                    e += 1;
                }
            }
            let len = e - i;
            if len >= MIN_RUN && matches!(k, 3 | 5 | 13 | 17) {
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
                13 => r.mul_add::<13>(x, s, y),
                17 => r.mul_add::<17>(x, s, y),
                _ => self.rows_add(r.r0, r.r0 + r.len, x, s, y),
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

    /// Correction from M L to the symmetric edge-flux Laplacian L_M,
    /// evaluated on w-wp. Each edge carries the arithmetic-mean nodal mass:
    /// (L_M x)_i = m_i (L x)_i + sum_j L_ij (m_j-m_i)(x_j-x_i)/2.
    /// Thus -L_M is an internal dissipative force even on a tapered plate.
    /// Dirichlet boundary edges keep their existing m_i weight. Constant mass
    /// gives exactly zero correction, preserving the uniform-body operator.
    fn mass_flux_correction(&self, i: usize, w: &[f64], wp: &[f64], mass: &[f64]) -> f64 {
        let xi = w[i] - wp[i];
        let mut correction = 0.0;
        for p in self.ptr[i]..self.ptr[i + 1] {
            let j = self.col[p];
            correction += 0.5 * self.val[p] * (mass[j] - mass[i]) * ((w[j] - wp[j]) - xi);
        }
        correction
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
        // Eliminate two adjacent columns together. Each destination still sees
        // exactly the descending-row subtraction order; only its load/store is
        // shared. Solve the second pivot after the first column updates it.
        let mut end = n;
        while end >= 2 && b > 0 {
            let i = end - 1;
            let j = i - 1;
            let lo = i.saturating_sub(b);
            let lower = j.saturating_sub(b);
            let row = &self.l[i * w + (lo + b - i)..i * w + b + 1];
            let previous = &self.l[j * w + (lower + b - j)..j * w + b + 1];
            let xi = x[i] / row[i - lo];
            x[i] = xi;
            x[j] -= row[j - lo] * xi;
            let xj = x[j] / previous[j - lower];
            x[j] = xj;
            for ((value, li), lj) in x[lo..j].iter_mut()
                .zip(&row[..j - lo])
                .zip(&previous[lo - lower..j - lower]) {
                *value -= li * xi;
                *value -= lj * xj;
            }
            // The earlier column may reach one extra node at the lower edge.
            for (value, lij) in x[lower..lo].iter_mut().zip(&previous[..lo - lower]) {
                *value -= lij * xj;
            }
            end -= 2;
        }
        for i in (0..end).rev() {
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

#[cfg(test)]
mod band_pair_tests {
    use super::*;

    #[test]
    fn adjacent_columns_preserve_order_and_linear_residual() {
        for n in [0_usize, 1, 2, 3, 5, 8, 17, 32] {
            for width in [0_usize, 1, 2, 4, 9, 31] {
                let width = width.min(n.saturating_sub(1));
                let mut a = vec![BTreeMap::new(); n];
                for i in 0..n {
                    for j in i.saturating_sub(width)..i {
                        let value = ((i * 17 + j * 13) % 11) as f64 / 32.0 - 0.125;
                        a[i].insert(j, value);
                        a[j].insert(i, value);
                    }
                }
                for (i, row) in a.iter_mut().enumerate() {
                    row.insert(i, 1.0 + row.values().map(|x| x.abs()).sum::<f64>());
                }
                let factor = Band::factor(&a).unwrap();
                for scale in [1e-9, 1.0, 1e9] {
                    let rhs: Vec<_> = (0..n).map(|i| scale * ((i * 19 % 23) as f64 - 11.0)).collect();
                    let mut expected = rhs.clone();
                    // Retained incumbent order, including both triangular
                    // sweeps, independently checks every output bit.
                    let (b, stride) = (factor.b, factor.b + 1);
                    for i in 0..n {
                        let lo = i.saturating_sub(b);
                        let row = &factor.l[i * stride + (lo + b - i)..i * stride + b + 1];
                        expected[i] = (expected[i] - dot(&row[..i - lo], &expected[lo..i])) / row[i - lo];
                    }
                    for i in (0..n).rev() {
                        let lo = i.saturating_sub(b);
                        let row = &factor.l[i * stride + (lo + b - i)..i * stride + b + 1];
                        let xi = expected[i] / row[i - lo];
                        expected[i] = xi;
                        for (value, lij) in expected[lo..i].iter_mut().zip(&row[..i - lo]) {
                            *value -= lij * xi;
                        }
                    }
                    let mut actual = rhs.clone();
                    factor.solve(&mut actual);
                    assert_eq!(actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                               expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
                    for (i, row) in a.iter().enumerate() {
                        let product: f64 = row.iter().map(|(&j, &v)| v * actual[j]).sum();
                        assert!((product - rhs[i]).abs() <= scale * 1e-12,
                                "n={n} width={width} row={i}: residual {}", product - rhs[i]);
                    }
                }
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
    rest_curv: Vec<(f64, f64, f64)>, // xx, yy, xy of z0 at core
    z0: Vec<f64>,              // rest height per node
    airy: Option<Band>,
    matrix: Csr,
    annulus: Option<annular_membrane::AnnularMembrane>,
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
    modal: Option<ModalLoss>,
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
    /// Exterior air mass, radiation resistance and viscous layer of the head shapes
    /// it retains: the table's m >= 1 set, or every resolved shape once the head's
    /// own loss is modelled (the volume shapes then less their compact share).
    head_radiation: Option<Box<radiation::ModalAirLoad>>,
    radiated_j: f64,
    /// The film's generalized-Maxwell law on the bending operator (`film`), if attached.
    film: Option<Box<film::FilmField>>,
    /// Whether `stiff` also holds support springs (it is then not the bending operator alone).
    support_springs: bool,
    /// Whether `s0`/`s1` hold the legacy fitted T60 surface loss, which lumps every loss path.
    surface_t60: bool,
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

/// A per-mode loss on a body whose modes are known in closed form (a coil's sine modes):
/// C = sum_j 2 sigma_j M phi_j phi_j' M, applied explicitly to the backward velocity like
/// `ViscousBending`, so its admission is the same dt^2 lambda_K + 2 dt lambda_C < 4 with
/// lambda_C = 2 max sigma_j exactly.
#[derive(Clone)]
struct ModalLoss {
    c: Csr,
    rate_bound: f64,
    dissipated_j: f64,
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
        // Average four positive tensor energies, not the twisting curvature
        // before squaring it. Centered Dxy loses alternating cell twists:
        // its squared Fourier symbol has extra cos²(kx h/2) cos²(ky h/2).
        // Cell D+x D+y instead satisfies the discrete integration-by-parts
        // identity; constant-D interior bending is exactly D*Laplacian_h².
        // Symmetric quarter-cell weights retain second-order consistency and
        // the energy's natural free-edge conditions without boundary clamps.
        let stencil = |i: i64, j: i64, di: i64, dj: i64| -> [Vec<(Option<usize>, f64)>; 3] {
            let h2 = h * h;
            let (ci, cj) = (i + di, j + dj);
            [
                vec![(at(i - 1, j), 1.0 / h2), (at(i, j), -2.0 / h2), (at(i + 1, j), 1.0 / h2)],
                vec![(at(i, j - 1), 1.0 / h2), (at(i, j), -2.0 / h2), (at(i, j + 1), 1.0 / h2)],
                vec![
                    (at(ci + 1, cj + 1), 1.0 / h2),
                    (at(ci, cj + 1), -1.0 / h2),
                    (at(ci + 1, cj), -1.0 / h2),
                    (at(ci, cj), 1.0 / h2),
                ],
            ]
        };
        if s.bending {
            let q = [[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 2.0 * (1.0 - nu)]];
            for &(i, j) in &curv_pts {
                let r = ((i * i + j * j) as f64).sqrt() * h;
                let t = (s.thickness)(r);
                let d = s.young * t * t * t / (12.0 * (1.0 - nu * nu));
                for (di, dj) in [(-1, -1), (0, -1), (-1, 0), (0, 0)] {
                    let b = stencil(i, j, di, dj);
                    for a in 0..3 {
                        for c in 0..3 {
                            if q[a][c] == 0.0 {
                                continue;
                            }
                            for &(pa, va) in &b[a] {
                                for &(pc, vc) in &b[c] {
                                    if let (Some(pa), Some(pc)) = (pa, pc) {
                                        add(&mut krows, pa, pc, 0.25 * d * h * h * q[a][c] * va * vc);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let viscosity = if let Some((source, angle)) = viscosity {
            let mut roots = vec![BTreeMap::new(); BENDING_ROOT_ROWS_PER_NODE * n];
            for &(i, j) in &curv_pts {
                let r = ((i * i + j * j) as f64).sqrt() * h;
                let d = source.bending_matrix_n_m_s((s.thickness)(r), angle)
                    .map_err(|e| e.to_string())?;
                if d.iter().all(|v| *v == 0.0) { continue; }
                let rows: Vec<BTreeMap<usize, f64>> = (0..3).map(|a|
                    (0..3).map(|b| (b, d[3*a+b])).collect()).collect();
                let factor = Band::factor(&rows)?;
                let node = at(i, j).expect("curvature point belongs to the plate");
                // Engineering curvature is [w_xx,w_yy,2w_xy]. Each root
                // includes sqrt(quarter area)=h/2, so G'G has units N s/m.
                for (quadrant, (di, dj)) in [(-1, -1), (0, -1), (-1, 0), (0, 0)].into_iter().enumerate() {
                    let b = stencil(i, j, di, dj);
                    for column in 0..3 {
                        for a in column..3 {
                            if a - column > factor.b { continue; }
                            let l = factor.l[a * (factor.b + 1) + column + factor.b - a];
                            let scale = 0.5 * h * l * if a == 2 { 2.0 } else { 1.0 };
                            for &(p, value) in &b[a] {
                                if let Some(p) = p {
                                    add(&mut roots, BENDING_ROOT_ROWS_PER_NODE * node + 3 * quadrant + column, p, scale * value);
                                }
                            }
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
                root: Csr::from_rows(&roots), scratch: vec![0.0; BENDING_ROOT_ROWS_PER_NODE*n], rate_bound,
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
            let rest_curv: Vec<_> = nb.iter().map(|s9| curv(&z0, s9)).collect();
            let l0: Vec<f64> = rest_curv.iter().map(|&(a, b, c)| a * b - c * c).collect();
            // Complementary membrane energy, with F=0 outside the core:
            // 1/(2 E t) [Fxx^2 + Fyy^2 - 2 nu Fxx Fyy + 2(1+nu) Fxy^2].
            // Keep the existing mean E*t normalization of q and V, but put
            // the LOCAL reciprocal thickness inside the stress operator.
            // Mixed derivatives live at cell centres. For constant thickness,
            // sum Fxx Fyy = sum Fxy^2 (zero extension), so this reduces exactly
            // to the former Laplacian-squared matrix, including its boundary.
            let nc = sorted.len();
            let mut arows = vec![BTreeMap::new(); nc];
            let h2 = h * h;
            let compliance = |x: f64, y: f64| -> Result<f64, String> {
                let t = (s.thickness)(x.hypot(y));
                if !t.is_finite() || t <= 0.0 {
                    return Err(format!("{name}: membrane thickness must be finite and positive"));
                }
                Ok(t_mean / t)
            };
            for &(i, j) in grid_index.keys() {
                let mut rows = [Vec::new(), Vec::new()];
                for (axis, offsets) in [[(-1, 0), (0, 0), (1, 0)], [(0, -1), (0, 0), (0, 1)]]
                    .iter()
                    .enumerate()
                {
                    for (&(di, dj), weight) in offsets.iter().zip([1.0, -2.0, 1.0]) {
                        if let Some(&c) = cidx.get(&(i + di, j + dj)) {
                            rows[axis].push((c, weight / h2));
                        }
                    }
                }
                let weight = compliance(i as f64 * h, j as f64 * h)?;
                for (a, ra) in rows.iter().enumerate() {
                    for (b, rb) in rows.iter().enumerate() {
                        let q = if a == b { 1.0 } else { -nu };
                        for &(ia, va) in ra {
                            for &(ib, vb) in rb {
                                add(&mut arows, ia, ib, weight * q * va * vb);
                            }
                        }
                    }
                }
            }
            let cells: std::collections::BTreeSet<_> = sorted
                .iter()
                .flat_map(|&(i, j)| [(i - 1, j - 1), (i - 1, j), (i, j - 1), (i, j)])
                .collect();
            for (i, j) in cells {
                let row: Vec<_> = [(0, 0, 1.0), (1, 0, -1.0), (0, 1, -1.0), (1, 1, 1.0)]
                    .iter()
                    .filter_map(|&(di, dj, v)| cidx.get(&(i + di, j + dj)).map(|&c| (c, v / h2)))
                    .collect();
                let weight = 2.0 * (1.0 + nu) * compliance((i as f64 + 0.5) * h, (j as f64 + 0.5) * h)?;
                for &(a, va) in &row {
                    for &(b, vb) in &row {
                        add(&mut arows, a, b, weight * va * vb);
                    }
                }
            }
            let airy = if s.hole > 0.0 { None } else { Some(Band::factor(&arows)?) };
            let (annulus, matrix) = if s.hole > 0.0 {
                let (a, matrix) = annular_membrane::AnnularMembrane::new(&grid_index, &sorted, s, compliance)?;
                (Some(a), matrix)
            } else {
                (None, Csr::from_rows(&arows))
            };
            let na = matrix.n;
            Some(VonKarman { eh: s.young * t_mean, nb, l0, rest_curv, z0, airy, matrix, annulus, f: vec![0.0; na], q: vec![0.0; na], curv: vec![(0.0, 0.0, 0.0); nc] })
        } else {
            None
        };
        // edge nodes of a free plate outside the true disc carry real mass but a staircase rim
        for m_ in mass.iter_mut() {
            *m_ = m_.max(1e-12);
        }
        Ok(Body {
            name: name.into(), h, n, xy, rest_z, grid_index, mass, s0: s0v, s1, lap, stiff, viscosity, modal: None, tension: s.tension, berger, vk,
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: true, area, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
            film: None, support_springs: s.support.is_some(), surface_t60: s0 > 0.0 || s1 > 0.0,
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
            s0: vec![6.907755 / t60; n], s1: 0.0, lap: Csr::from_rows(&lrows), stiff: Csr::from_rows(&vec![BTreeMap::new(); n]), viscosity: None, modal: None,
            tension: tension / h, berger: 0.0, vk: None, // tension/h: the chain energy T/2 sum (dw)^2/h
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
            film: None, support_springs: false, surface_t60: 6.907755 / t60 > 0.0,
        }
    }

    /// A coiled snare wire: the `string` chain plus the coil's own bending stiffness `bending`
    /// [N m^2], with both ends pinned beyond the last node (w = 0 there, curvature free). The
    /// bending operator B/h^3 lap^2 then shares the chain's sine modes phi_j, with
    /// beta_j^2 = 4/h^2 sin^2(pi j / (2 (n + 1))) and omega_j^2 = (T beta_j^2 + B beta_j^4)/mu,
    /// so the loss is exact per mode: `sigma(omega_j, share_j)` [1/s] at each mode's own discrete
    /// frequency, where share_j = B beta_j^2/(T + B beta_j^2) is the bending part of its potential
    /// energy. No s0 or s1 is used.
    #[allow(clippy::too_many_arguments)]
    pub fn coil(name: &str, n: usize, h: f64, mass_per_m: f64, tension: f64, bending: f64,
        sigma: &dyn Fn(f64, f64) -> Result<f64, String>, x0: f64, y: f64) -> Result<Body, String> {
        if n < 2 || !(h > 0.0 && mass_per_m > 0.0 && tension >= 0.0 && bending > 0.0) {
            return Err(format!("{name}: a coil needs n >= 2 and positive spacing, mass and bending (tension >= 0)"));
        }
        let mut b = Body::string(name, n, h, mass_per_m, tension, f64::INFINITY, x0, y);
        b.s0 = vec![0.0; n];
        b.surface_t60 = false;
        // stiff = (B/h^3) lap lap, lap the chain's (2, -1) stencil with zero ghosts
        let lrows: Vec<BTreeMap<usize, f64>> = (0..n).map(|a| {
            let mut r = BTreeMap::new();
            r.insert(a, 2.0);
            if a > 0 { r.insert(a - 1, -1.0); }
            if a + 1 < n { r.insert(a + 1, -1.0); }
            r
        }).collect();
        let mut krows = vec![BTreeMap::new(); n];
        for (a, row) in lrows.iter().enumerate() {
            for (&mid, &x) in row {
                for (&c, &z) in &lrows[mid] {
                    add(&mut krows, a, c, bending / h.powi(3) * x * z);
                }
            }
        }
        b.stiff = Csr::from_rows(&krows);
        let m = mass_per_m * h;
        let np1 = (n + 1) as f64;
        let mut two_sigma = Vec::with_capacity(n);
        for j in 1..=n {
            let beta2 = 4.0 / (h * h) * (PI * j as f64 / (2.0 * np1)).sin().powi(2);
            let omega = ((tension * beta2 + bending * beta2 * beta2) / mass_per_m).sqrt();
            let share = bending * beta2 / (tension + bending * beta2);
            let s = sigma(omega, share)?;
            if !(s >= 0.0 && s.is_finite()) {
                return Err(format!("{name}: mode {j} at {:.1} Hz has loss {s} 1/s", omega / (2.0 * PI)));
            }
            two_sigma.push(2.0 * s);
        }
        let phi = |j: usize, i: usize| (2.0 / np1).sqrt() * (PI * (j * (i + 1)) as f64 / np1).sin();
        let mut crows = vec![BTreeMap::new(); n];
        for (r, row) in crows.iter_mut().enumerate() {
            for c in 0..n {
                let v: f64 = (1..=n).map(|j| two_sigma[j - 1] * phi(j, r) * phi(j, c)).sum();
                row.insert(c, m * v);
            }
        }
        let rate_bound = two_sigma.iter().cloned().fold(0.0_f64, f64::max);
        b.modal = Some(ModalLoss { c: Csr::from_rows(&crows), rate_bound, dissipated_j: 0.0 });
        Ok(b)
    }

    /// Energy the coil's per-mode loss removed at centred velocity [J].
    pub fn modal_loss_dissipated_j(&self) -> f64 {
        self.modal.as_ref().map_or(0.0, |m| m.dissipated_j)
    }

    /// A lumped rigid mass (stick, beater, hi-hat rod).
    pub fn point(name: &str, mass: f64) -> Body {
        Body {
            name: name.into(), h: 1.0, n: 1, xy: vec![[0.0, 0.0]], rest_z: vec![0.0], grid_index: BTreeMap::new(), mass: vec![mass],
            s0: vec![0.0], s1: 0.0, lap: Csr::from_rows(&[BTreeMap::new()]), stiff: Csr::from_rows(&[BTreeMap::new()]), viscosity: None, modal: None,
            tension: 0.0, berger: 0.0, vk: None,
            w: vec![0.0], wp: vec![0.0], wn: vec![0.0], f: vec![0.0], bias: vec![0.0], tmp: vec![0.0], lw: vec![0.0], lwp: vec![0.0], g: vec![0.0], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0, head_radiation: None, radiated_j: 0.0,
            film: None, support_springs: false, surface_t60: false,
        }
    }

    /// Attach the exterior load of a clamped full circular head of `radius`
    /// (its nodes centred on the rim's centre) from `table` in `gas`. Mass
    /// and stiffness are this body's own; returns the receipt lines. The
    /// optional compact net-volume resistance [Pa s/m^3] is compared with the
    /// baked volume shapes in the receipt; on a head carrying a film law (the
    /// load then retains every resolved shape) the volume shapes (0,n>=2) also
    /// take the baffled load less its share.
    pub fn attach_head_radiation(&mut self, table: &radiation::RadiationTable, radius: f64, gas: &GasState,
        compact_volume_resistance: Option<f64>) -> Result<Vec<String>, String> {
        self.attach_head_radiation_impl(table, radius, gas, compact_volume_resistance, None)
    }

    pub fn attach_head_radiation_with_viscous(
        &mut self,
        table: &radiation::RadiationTable,
        radius: f64,
        gas: &GasState,
        compact_volume_resistance: Option<f64>,
        coincidence_hz: f64,
    ) -> Result<Vec<String>, String> {
        self.attach_head_radiation_impl(
            table,
            radius,
            gas,
            compact_volume_resistance,
            Some(coincidence_hz),
        )
    }

    fn attach_head_radiation_impl(
        &mut self,
        table: &radiation::RadiationTable,
        radius: f64,
        gas: &GasState,
        compact_volume_resistance: Option<f64>,
        coincidence_hz: Option<f64>,
    ) -> Result<Vec<String>, String> {
        // Air loading is the same physics whatever the head's own loss: a
        // sourced bending viscosity or a film law composes with it (its force
        // is explicit, the air mass implicit). Shallow shells are not clamped films.
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
        let viscous = coincidence_hz.map(|fc| (self.mass.iter().sum::<f64>() / self.area, fc));
        // A head whose own loss is modelled (a film law on its bending operator,
        // the viscous layer) has its air mass and viscous layer on every shape
        // the grid resolves, not only the table's: the load then also needs the
        // piston's share of each volume shape.
        let extend = (self.film.is_some() || viscous.is_some()).then(|| {
            if self.air_mass > 0.0 {
                self.air_mass / (self.area * self.area)
            } else {
                0.0
            }
        });
        let mut load = radiation::ModalAirLoad::new_with_losses(
            table,
            radius,
            self.h,
            &self.xy,
            &self.mass,
            &stiffness,
            gas,
            compact_volume_resistance,
            viscous,
            extend,
        )?;
        if let Some(film) = self.film.as_deref() {
            // The load's loaded frequencies use the card's E; the law's E' differs (not iterated).
            // E' rises with frequency, so the retained band's ends bound the difference.
            let (low, high) = load.loaded_hz;
            let (at_low, at_high) = (film.law.storage(low), film.law.storage(high));
            let shift = (at_low.sqrt() - 1.0)
                .abs()
                .max((at_high.sqrt() - 1.0).abs());
            let note = format!(
                "; the film law is this head's own, on its bending operator, not in this load; each shape's load is read at its loaded frequency on the card's E, not the law's E', which is {:+.1} % against it at the lowest retained shape ({low:.0} Hz) and {:+.1} % at the highest ({high:.0} Hz); a shape's frequency differs from the one its load is read at by about half that times its bending share, at most {:.1} % (pure bending)",
                100.0 * (at_low - 1.0),
                100.0 * (at_high - 1.0),
                100.0 * shift
            );
            if let Some(scope) = load
                .receipt
                .iter_mut()
                .find(|line| line.starts_with("head radiation scope:"))
            {
                scope.push_str(&note);
            }
        }
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

    /// Make this plate's bending operator the film's generalized-Maxwell law
    /// `law` (moduli relative to the plate's own E, the card's), on every grid
    /// mode. The tension (a pre-stress) and the Berger in-plane stretching
    /// stay elastic. Returns the receipt lines. A shallow shell, a plate
    /// without bending, one whose stiffness operator also holds support
    /// springs, one with another bending loss or the fitted T60 surface loss
    /// (which already lumps the film's), a second law, or a law attached after
    /// the head radiation refuses.
    pub fn attach_film_maxwell(&mut self, law: film::FilmLaw) -> Result<Vec<String>, String> {
        if self.vk.is_some() || self.n < 2 || self.stiff.col.is_empty() {
            return Err(format!(
                "{}: a film law needs a flat plate's bending operator",
                self.name
            ));
        }
        if self.support_springs {
            return Err(format!(
                "{}: this stiffness operator also holds support springs, not the bending operator alone",
                self.name
            ));
        }
        if self.film.is_some() {
            return Err(format!("{}: a film law attaches once", self.name));
        }
        if self.viscosity.is_some() || self.modal.is_some() || self.surface_t60 {
            return Err(format!(
                "{}: a film law and another bending loss, or the fitted T60 surface loss that lumps it, are exclusive",
                self.name
            ));
        }
        if self.head_radiation.is_some() {
            return Err(format!(
                "{}: attach the film law before the head radiation, so that the load retains every shape the grid resolves",
                self.name
            ));
        }
        // The highest grid mode at the card's E (tension and bending), before the law is attached.
        let top_hz = self.max_rate2().sqrt() / (2.0 * PI);
        let (lo, hi) = law.band_hz;
        let place = if top_hz < lo {
            format!(
                "below the curve's {lo} Hz, where the law continues it, law tan delta {:.4} there",
                law.tan_delta(top_hz)
            )
        } else if top_hz <= hi {
            format!(
                "inside the curve's band, law tan delta {:.4} there",
                law.tan_delta(top_hz)
            )
        } else {
            format!(
                "ABOVE the curve's {hi} Hz: the modes above it take the law's own continuation, not the curve"
            )
        };
        let lines = vec![
            law.receipt(),
            format!(
                "film law on the bending operator: every grid mode carries it, the highest at {top_hz:.0} Hz ({place}); below {lo} Hz the law continues the curve; tension (a pre-stress) and the Berger in-plane stretching stay elastic"
            ),
            "film law time stepping: fs-material's exact linear-ramp branch recursion (GeneralizedMaxwell::step) on the field K q; its heat (dt/tau)|q_mid|^2_K per branch closes the discrete energy, and admission uses the glassy modulus".to_string(),
        ];
        self.film = Some(Box::new(film::FilmField::new(law, self.n)));
        Ok(lines)
    }

    /// The attached film law, if any.
    pub fn film_law(&self) -> Option<&film::FilmLaw> {
        self.film.as_deref().map(|f| &f.law)
    }

    /// Heat the film law's branches have dissipated [J] (the midpoint heat that
    /// closes `discrete_energy`; see `film`).
    pub fn film_heat_j(&self) -> f64 {
        self.film.as_deref().map_or(0.0, |f| f.heat_j)
    }

    /// The film law's heat in the last step [J].
    pub fn film_step_j(&self) -> f64 {
        self.film.as_deref().map_or(0.0, |f| f.last_j)
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
            if let Some(film) = self.film.as_deref_mut() {
                // K w is kept for the film law's branches (`predict`); f gets the
                // same -K w as below, and the law's difference from it in `predict`.
                film.sw.fill(0.0);
                self.stiff.mul_add(&self.w, 1.0, &mut film.sw);
                for (f, s) in self.f.iter_mut().zip(&film.sw) {
                    *f -= s;
                }
            } else {
                self.stiff.mul_add(&self.w, -1.0, &mut self.f);
            }
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
            let VonKarman { eh, nb, rest_curv, airy, f, q, curv, annulus, z0, .. } = vk;
            let eh = *eh;
            if let Some(annulus) = annulus {
                annulus.prepare_load(w, z0);
                annulus.solve_with(f, |c| {
                    let qc = annulus.prepared_load_row(c, eh);
                    q[c] = qc;
                    qc
                });
                v += h2 / (2.0 * eh) * f.iter().zip(q.iter()).map(|(a, b)| a * b).sum::<f64>();
                annulus.gradient(w, z0, f, h2, &mut self.g);
                return v.max(0.0);
            }
            let u = |o: Option<usize>| o.map_or(0.0, |k| w[k]);
            let airy = airy.as_ref().expect("simply-connected plate retains its Airy factor");
            // q = -E t (l(w + z0) - l(z0)) per core node, built row by row inside the Airy solve's
            // forward sweep (it fills that sweep's latency chain); A F = q
            airy.solve_with(f, |c| {
                let s9 = &nb[c];
                let xx = (u(s9[1]) - 2.0 * u(s9[0]) + u(s9[2])) / h2;
                let yy = (u(s9[3]) - 2.0 * u(s9[0]) + u(s9[4])) / h2;
                let xy = (u(s9[5]) - u(s9[6]) - u(s9[7]) + u(s9[8])) / (4.0 * h2);
                // Expand l(w + z0) - l(z0) before evaluating: adding a tiny
                // displacement to z0 can erase it, as can subtracting the
                // two nearly equal curvature products. Keep the rest tangent.
                let (zx, zy, zxy) = rest_curv[c];
                let qc = -eh * (xx * yy - xy * xy + xx * zy + yy * zx - 2.0 * xy * zxy);
                q[c] = qc;
                curv[c] = (xx + zx, yy + zy, xy + zxy);
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
        if let Some(film) = self.film.as_deref_mut() {
            // The branches advance once per step here (not in `linear_forces`,
            // which the admission's power iteration also calls).
            film.advance(&self.w, &self.wp, &mut self.f, k, &self.stiff);
        }
        if let Some(viscous) = &self.viscosity {
            for i in 0..n { self.tmp[i] = self.w[i] - self.wp[i]; }
            // Backward velocity preserves the existing diagonal contact and
            // SAV compliance. Its additional explicit stability bound is
            // dt² lambda_K + 2 dt lambda_C < 4; see discrete_energy below.
            viscous.c.mul_add(&self.tmp, -1.0 / k, &mut self.f);
        }
        if let Some(modal) = &self.modal {
            // the same backward velocity and admission as the viscous law
            for i in 0..n { self.tmp[i] = self.w[i] - self.wp[i]; }
            modal.c.mul_add(&self.tmp, -1.0 / k, &mut self.f);
        }
        if self.s1 > 0.0 {
            // div(m grad velocity), not m Lap(velocity): the latter is not
            // symmetric when thickness varies and can create mechanical energy.
            let ih2 = 1.0 / (self.h * self.h);
            for i in 0..self.n {
                let correction = self.lap.mass_flux_correction(i, &self.w, &self.wp, &self.mass);
                self.tmp[i] = -((self.lw[i] - self.lwp[i]) + correction / self.mass[i]) * ih2;
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
                // Zero is a physical previous potential, not an uninitialized
                // sentinel: keep the same half-step average through zero.
                let target = 0.5 * (root + self.r_prev);
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
        if let Some(modal) = &mut self.modal {
            for i in 0..self.n { self.tmp[i] = self.wn[i] - self.wp[i]; }
            // dt v_center' C v_center, v_center = (wn - wp)/(2 dt); C is positive semidefinite
            let loss = modal.c.quad(&self.tmp) / (4.0 * k);
            modal.dissipated_j += loss;
            bending_loss += loss;
        }
        if let Some(load) = self.head_radiation.as_deref_mut() {
            let radiated = load.radiated(&self.wn, k);
            self.radiated_j += radiated;
            bending_loss += radiated;
        }
        if self.vk.is_some() {
            // energy the loss terms removed over this step (the budget for raising psi next step):
            // Refresh at zero nonlinear gradient too: otherwise a zero-work
            // step retains credit already allocated by the previous step.
            // s0: sum m s0 (wn - wp)^2 / (2k); s1 uses the same symmetric
            // mass-weighted edge flux as predict, paired with centered velocity.
            let ih2 = 1.0 / (self.h * self.h);
            let mut d0 = 0.0;
            let mut d1 = 0.0;
            for i in 0..self.n {
                let dv = self.wn[i] - self.wp[i];
                d0 += self.mass[i] * self.s0[i] * dv * dv;
                let correction = self.lap.mass_flux_correction(i, &self.w, &self.wp, &self.mass);
                d1 += (self.mass[i] * (self.lw[i] - self.lwp[i]) + correction) * dv;
            }
            // Backward s1 velocity can do negative work against centered
            // motion. Clip the signed total, not its separate contributions.
            let loss = d0 / (2.0 * k) + self.s1 * ih2 / k * d1 + bending_loss;
            self.budget = if loss < 0.0 { 0.0 } else { loss }; // preserve NaN
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
    /// A film law's branches relax: the set state has no film history.
    pub fn sync(&mut self) {
        for x in self.lwp.iter_mut() {
            *x = 0.0;
        }
        self.lap.mul_add(&self.wp, 1.0, &mut self.lwp);
        if let Some(film) = self.film.as_deref_mut() {
            film.reset();
        }
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
            if let Some(film) = self.film.as_deref() {
                // A film law's bending stiffness reaches its glassy modulus within a step.
                let excess = film.law.glassy() - 1.0;
                for (f, s) in self.f.iter_mut().zip(&film.sw) {
                    *f -= excess * s;
                }
            }
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
        // Annular strain is quadratic in displacement. Its rest derivative
        // gives the exact tangent with one constant sparse Jacobian, rather
        // than rebuilding the full nonlinear strain at every power iterate.
        // This changes setup only; rendering still evaluates the full energy.
        let mut tangent = self.vk.as_ref().and_then(|vk| vk.annulus.as_ref()
            .map(|a| (a.rest_jacobian(&vk.z0), vec![0.0; vk.f.len()])));
        for _ in 0..400 {
            let vn = v.iter().map(|a| a * a).sum::<f64>().sqrt();
            for i in 0..n {
                self.w[i] = v[i] / vn / self.mass[i].sqrt();
            }
            let bias = std::mem::replace(&mut self.bias, vec![0.0; n]);
            self.linear_forces();
            self.bias = bias;
            let gradient_divisor = if let Some((j, airy)) = tangent.as_mut() {
                let vk = self.vk.as_ref().unwrap();
                vk.annulus.as_ref().unwrap().rest_hessian_apply(j, &self.w,
                    vk.eh, self.h * self.h, airy, &mut self.g);
                1.0
            } else {
                // Simply-connected/other bodies retain their existing path.
                let xs = self.w.clone();
                for i in 0..n {
                    self.w[i] = eps * xs[i];
                }
                self.nonlinear_gradient();
                eps
            };
            let y: Vec<f64> = (0..n).map(|i| (-self.f[i] + self.g[i] / gradient_divisor) / self.mass[i].sqrt()).collect();
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

    /// K x alone (bending, tension, supports; the snare-wire bias excluded), without the
    /// curvature Hessian: the linear part of [`Body::rest_stiffness_apply`]. State is restored.
    pub fn linear_stiffness_apply(&mut self, x: &[f64]) -> Result<Vec<f64>, String> {
        let n = self.n;
        if x.len() != n || x.iter().any(|v| !v.is_finite()) {
            return Err(format!("{}: linear stiffness needs {n} finite nodal values", self.name));
        }
        let saved = (self.w.clone(), self.f.clone(), self.lw.clone());
        let bias = std::mem::replace(&mut self.bias, vec![0.0; n]);
        self.w.copy_from_slice(x);
        self.linear_forces();
        let out: Vec<f64> = self.f.iter().map(|f| -f).collect();
        self.bias = bias;
        (self.w, self.f, self.lw) = saved;
        Ok(out)
    }

    /// (K + B) x about the rest state, for offline modal analysis: the linear forces (bending,
    /// tension, supports; the snare-wire bias excluded) plus the Hessian B of the nonlinear
    /// potential at w = 0, by a central difference of its gradient,
    /// B x = (grad V(eps x) - grad V(-eps x)) / (2 eps), eps = 1e-12 m. The difference cancels the
    /// quadratic term, leaving an O(eps^2) truncation; the Airy load expands l(w + z0) - l(z0)
    /// analytically, so a small step costs no cancellation. Measured on the kit cymbals (.372):
    /// K's asymmetry falls as eps^2 down to its roundoff floor near 2e-15 at 1e-12 m. The same two
    /// owners (linear_forces, nonlinear_gradient); the body's state is restored.
    pub fn rest_stiffness_apply(&mut self, x: &[f64]) -> Result<Vec<f64>, String> {
        self.rest_stiffness_apply_step(x, 1e-12)
    }

    /// [`Body::rest_stiffness_apply`] with an explicit central-difference step `eps` [m], for
    /// attributing its error: truncation scales as eps^2, the Airy solve's roundoff does not.
    pub fn rest_stiffness_apply_step(&mut self, x: &[f64], eps: f64) -> Result<Vec<f64>, String> {
        let n = self.n;
        if !(eps.is_finite() && eps > 0.0) {
            return Err(format!("{}: central-difference step must be positive and finite", self.name));
        }
        if x.len() != n || x.iter().any(|v| !v.is_finite()) {
            return Err(format!("{}: rest stiffness needs {n} finite nodal values", self.name));
        }
        let saved = (self.w.clone(), self.g.clone(), self.f.clone(), self.lw.clone());
        let bias = std::mem::replace(&mut self.bias, vec![0.0; n]);
        self.w.copy_from_slice(x);
        self.linear_forces();
        let mut out: Vec<f64> = self.f.iter().map(|f| -f).collect();
        for sign in [1.0, -1.0] {
            for i in 0..n {
                self.w[i] = sign * eps * x[i];
            }
            self.linear_forces(); // the Berger term reads this w's Laplacian
            self.nonlinear_gradient();
            for i in 0..n {
                out[i] += sign * self.g[i] / (2.0 * eps);
            }
        }
        self.bias = bias;
        (self.w, self.g, self.f, self.lw) = saved;
        Ok(out)
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
    /// With a film law it includes the branches' storage, and the law's midpoint heat is one more loss.
    pub fn discrete_energy(&self, k: f64) -> f64 {
        let velocity: Vec<f64> = self.w.iter().zip(&self.wp).map(|(w, wp)| (w - wp) / k).collect();
        let mut kw = vec![0.0; self.n];
        // A film law keeps e_inf of the card's bending stiffness; its branches hold the rest.
        let bending = self.film.as_deref().map_or(1.0, |f| f.e_inf());
        self.stiff.mul_add(&self.w, bending, &mut kw);
        self.lap.mul_add(&self.w, self.tension, &mut kw);
        let film = self
            .film
            .as_deref()
            .map_or(0.0, |f| f.stored_j(&self.w, &self.wp));
        let kinetic = self.mass.iter().zip(&velocity).map(|(m, v)| 0.5*m*v*v).sum::<f64>();
        let correction = self.viscosity.as_ref().map_or(0.0, |v| k * v.c.quad(&velocity) / 4.0)
            + self.modal.as_ref().map_or(0.0, |m| k * m.c.quad(&velocity) / 4.0);
        // Head radiation: the shapes' air kinetic energy, less the same
        // backward-velocity damping correction as the viscous law.
        let air = self.head_radiation.as_deref().map_or(0.0, |load| {
            let p: Vec<f64> = load.cols.iter().map(|c| dot(c, &velocity)).collect();
            let (mut a, mut b) = (vec![0.0; p.len()], vec![0.0; p.len()]);
            radiation::ModalAirLoad::apply(&load.mass, &p, &mut a);
            radiation::ModalAirLoad::apply(&load.damping, &p, &mut b);
            0.5 * dot(&p, &a) - k * dot(&p, &b) / 4.0
        });
        kinetic + air + film + 0.5 * dot(&self.wp, &kw) + 0.5 * self.psi * self.psi - correction
    }

    /// Kinetic + linear potential + auxiliary (nonlinear) energy at the half step (diagnostics).
    pub fn energy(&self, k: f64) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n {
            let v = (self.w[i] - self.wp[i]) / k;
            e += 0.5 * self.mass[i] * v * v;
        }
        let mid: Vec<f64> = (0..self.n).map(|i| 0.5 * (self.w[i] + self.wp[i])).collect();
        let (bending, film) = self
            .film
            .as_deref()
            .map_or((1.0, 0.0), |f| (f.e_inf(), f.elastic_j()));
        e + 0.5 * self.tension * self.lap.quad(&mid)
            + 0.5 * bending * self.stiff.quad(&mid)
            + film
            + 0.5 * self.psi * self.psi
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

#[cfg(test)]
mod link_clock_tests {
    use super::*;

    #[test]
    fn assembled_signed_links_retain_the_analytical_relative_mode() {
        let mut ins = Instrument::new(48_000.0);
        ins.add(Body::point("a", 1.0));
        ins.add(Body::point("b", 1.0));
        // Together these are one positive spring/dashpot, like signed terms
        // in an affine mount. Separately taking magnitudes gives a false bound.
        ins.links.push(Link { a: (0, 0), b: Some((1, 0)), k: 10.0, c: 3.0 });
        ins.links.push(Link { a: (0, 0), b: Some((1, 0)), k: -7.0, c: -2.0 });
        for rate in [48_000.0, 96_000.0] {
            assert_eq!(ins.link_admission_rate2(rate), 6.0 + 4.0 * rate);
        }
        assert_eq!(ins.admission_rate2().0, 192_006.0);
    }

    #[test]
    fn link_clock_matches_real_explicit_dashpot_instability() {
        let run = |rate| {
            let mut ins = Instrument::new(rate);
            ins.add(Body::point("a", 1.0));
            ins.add(Body::point("b", 1.0));
            ins.links.push(Link { a: (0, 0), b: Some((1, 0)), k: 0.0, c: 96_000.0 });
            ins.bodies[0].w[0] = 1e-6;
            ins.bodies[1].w[0] = -1e-6;
            let ratio = ins.stability().0;
            let mut peak = 0.0_f64;
            for _ in 0..32 {
                ins.step(&mut []);
                peak = peak.max(ins.bodies[0].w[0].abs());
            }
            (ratio, peak)
        };
        let (bad_ratio, bad_peak) = run(48_000.0);
        let (good_ratio, good_peak) = run(192_000.0);
        assert!((bad_ratio - 2.0).abs() < 1e-14);
        assert!((good_ratio - 0.5).abs() < 1e-14);
        assert!(bad_peak > 1.0, "the omitted dashpot really destabilizes the step");
        assert!(good_peak <= 1e-6, "the admitted clock must keep the actual trajectory bounded");
    }

    #[test]
    fn invalid_link_coefficients_cannot_disappear_from_the_bound() {
        let mut ins = Instrument::new(48_000.0);
        ins.add(Body::point("a", 1.0));
        ins.links.push(Link { a: (0, 0), b: None, k: f64::NAN, c: 0.0 });
        assert!(ins.admission_rate2().0.is_infinite());
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Law {
    /// f = K d^alpha (1 + mu d_dot): Hertz with Hunt-Crossley loss.
    Power { k: f64, alpha: f64, mu: f64 },
}

/// A contact whose softer body crushes: the stick's wood across its grain under a cymbal's edge
/// (`shaft`). The contact's `Law::Power` (alpha 1.5) is its elastic Hertz law F = K d^(3/2), whose
/// contact area is `area` d. Loading is Hertz until the mean pressure F/A reaches `pressure`, at the
/// onset d_y = (pressure area / K)^2. Deeper, the mean pressure stays there and the area keeps
/// Hertz's proportion to the approach, F = pressure area d (Thornton 1997, J. Appl. Mech. 64, 383,
/// whose limit is the peak pressure; here it is the mean pressure, which a hardness test reports).
/// Below the deepest approach d_m the contact unloads and reloads along Hertz from a permanent set,
/// with the curvature that keeps the loaded area at the peak force (Thornton's construction):
/// F = K (d_m/d_y)^(1/4) (d - s)^(3/2) with s = d_m - sqrt(d_m d_y). At every approach a deeper
/// memory gives no more force, so the contact's potential only falls as the wood crushes and the
/// discrete-gradient solve stays passive. Only the Hunt-Crossley loss of the law is added to it.
/// The solve's secant spans two steps, d^(n-1) to d^(n+1), so it reads the memory through d^(n-1):
/// a loading step then lies wholly on the crushing line, and the crush is exact to second order
/// (a memory through d^n would put the step's lower half on an unloading curve, a first-order
/// error). Each launch crushes fresh wood.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Crush {
    /// Contact area per unit approach [m].
    pub area: f64,
    /// Mean contact pressure at which the softer body crushes [Pa].
    pub pressure: f64,
    /// The deepest approach so far [m], through the previous step: the contact's memory.
    pub deepest: f64,
}

impl Crush {
    pub fn new(area: f64, pressure: f64) -> Result<Crush, String> {
        if !(area.is_finite() && area > 0.0 && pressure.is_finite() && pressure > 0.0) {
            return Err(format!("crush needs a positive area per approach and crush pressure, got {area} m and {pressure} Pa"));
        }
        Ok(Crush { area, pressure, deepest: 0.0 })
    }

    /// The approach at which a Hertz contact of stiffness `k` [N/m^1.5] begins to crush [m].
    #[must_use]
    pub fn onset(&self, k: f64) -> f64 {
        let root = self.area * self.pressure / k;
        root * root
    }

    /// At the current memory: the unloading curve's stiffness [N/m^1.5], its permanent set [m] and
    /// the approach where it meets the crushing line [m].
    fn branch(&self, k: f64) -> (f64, f64, f64) {
        let onset = self.onset(k);
        let reach = self.deepest.max(onset);
        (k * (reach / onset).sqrt().sqrt(), reach - (reach * onset).sqrt(), reach)
    }

    /// The permanent set the contact keeps once unloaded [m].
    #[must_use]
    pub fn set(&self, k: f64) -> f64 {
        self.branch(k).1
    }

    /// The force at approach `d` at the current memory, without Hunt-Crossley loss [N].
    #[must_use]
    pub fn force(&self, k: f64, d: f64) -> f64 {
        let (stiff, set, reach) = self.branch(k);
        if d > reach {
            self.pressure * self.area * d
        } else {
            let x = (d - set).max(0.0);
            stiff * x * x.sqrt()
        }
    }
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
    /// The force solve may stop when a nonzero-residual Newton correction rounds
    /// to the current iterate (see `contact_force`). Exact zero residuals always
    /// stop, independently of this option. On for the snare
    /// wires: held under load and warm-started, about one solve in five lands on its own root. Off
    /// elsewhere: it is as exact there, but it would change every other piece's output (by rounding
    /// for most; the open hi-hat's rim chatter decorrelates), and that has not been graded.
    pub fixed_point_exit: bool,
    /// The crushing of the softer body, with its memory; None for an elastic contact.
    pub crush: Option<Crush>,
}

impl Contact {
    pub fn new(a: Vec<(usize, usize, f64)>, b: Vec<(usize, usize, f64)>, gap: f64, law: Law) -> Contact {
        Contact { a, b, gap, law, d_now: 0.0, active: true, peak_force: 0.0, force: 0.0, fixed_point_exit: false, crush: None }
    }

    /// The same contact crushing as `crush` describes; its law must be Hertz (alpha 1.5).
    pub fn crushing(mut self, crush: Crush) -> Result<Contact, String> {
        let Law::Power { k, alpha, .. } = self.law;
        if alpha != 1.5 || !(k.is_finite() && k > 0.0) {
            return Err(format!("a crushing contact needs a Hertz law (alpha 1.5, K > 0), got alpha {alpha}, K {k}"));
        }
        self.crush = Some(crush);
        Ok(self)
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
    /// A rigid floor: head and cavity-port radiation both include the delayed,
    /// positive-pressure image path. Sources facing away from the microphone
    /// can still reach it by reflection (the snare's underside heard overhead).
    pub floor: Option<Floor>,
    hp: (f64, f64),
}

/// A rigid plane: the point `point` on it and its unit `normal`, pointing into the air.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Floor {
    point: [f64; 3],
    normal: [f64; 3],
}

/// The plane z = `height`, normal +z: the floor of a drum standing on its axis (z up).
pub fn floor_z(height: f64) -> Floor {
    Floor { point: [0.0, 0.0, height], normal: [0.0, 0.0, 1.0] }
}

/// The plane y = `height`, normal +y: the floor of a drum lying on its side (the kick, y up).
pub fn floor_y(height: f64) -> Floor {
    Floor { point: [0.0, height, 0.0], normal: [0.0, 1.0, 0.0] }
}

impl Floor {
    /// The mirror image of `p`, p - 2 (p.n) n + 2 (q.n) n, evaluated in that order so that an
    /// axis normal reflects exactly as the floor's former height did: 2 h - p along the axis, p
    /// across it.
    fn image(&self, p: [f64; 3]) -> [f64; 3] {
        let n = self.normal;
        let pn = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
        let qn = self.point[0] * n[0] + self.point[1] * n[1] + self.point[2] * n[2];
        std::array::from_fn(|i| p[i] - 2.0 * pn * n[i] + 2.0 * qn * n[i])
    }
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

/// The secant of a crushing contact's potential between `previous` and `d` at its current memory,
/// and its derivative in `d` (`power_secant` for `Crush`). The potential is convex, so the secant
/// increases with `d`.
fn crush_secant(d: f64, previous: f64, kk: f64, crush: &Crush) -> (f64, f64) {
    let (stiff, set, reach) = crush.branch(kk);
    let slope = crush.pressure * crush.area;
    if d <= reach && previous <= reach {
        let (value, derivative) = power_secant(d - set, previous - set, 1.5);
        return (stiff * value, stiff * derivative);
    }
    if d >= reach && previous >= reach {
        return (0.5 * slope * (d + previous), 0.5 * slope);
    }
    // one end on each branch: the unloading curve up to `reach`, the crushing line beyond it
    let (lo, hi) = if d < previous { (d, previous) } else { (previous, d) };
    let below = stiff * power_secant(reach - set, lo - set, 1.5).0 * (reach - lo);
    let above = 0.5 * slope * (hi + reach) * (hi - reach);
    let value = (below + above) / (hi - lo);
    (value, (crush.force(kk, d) - value) / (d - previous))
}

/// `contact_phi` for a crushing contact.
#[inline(always)]
fn crush_phi(d: f64, previous: f64, kk: f64, crush: &Crush, mu: f64, k: f64) -> (f64, f64) {
    let (elastic, derivative) = crush_secant(d, previous, kk, crush);
    let damping_rate = mu/(2.0*k);
    let damp = 1.0+damping_rate*(d-previous);
    let value = elastic*damp;
    if value <= 0.0 { (0.0, 0.0) }
    else { (value, derivative*damp+elastic*damping_rate) }
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
/// An exactly zero residual always stops the solve. With `fixed_exit`, a Newton
/// step that rounds to the iterate also stops the solve there: the
/// correction is below half an ulp, so the iterate is the root to rounding. Without it the bracket
/// update has just put that iterate on the bracket's edge, so the safeguard bisects away from the
/// root and halves back down to it (30-47 iterations instead of about 3).
fn contact_force(d_pred: f64, comp: f64, force: f64, previous: f64, law: Law, fixed_exit: bool, k: f64) -> f64 {
    let Law::Power { k: kk, alpha, mu } = law;
    solve_force(d_pred, comp, force, previous, fixed_exit, |d| contact_phi(d, previous, kk, alpha, mu, k))
}

/// `contact_force` for a crushing contact, at its current memory.
#[allow(clippy::too_many_arguments)]
fn crush_force(d_pred: f64, comp: f64, force: f64, previous: f64, law: Law, crush: &Crush, fixed_exit: bool, k: f64) -> f64 {
    let Law::Power { k: kk, mu, .. } = law;
    solve_force(d_pred, comp, force, previous, fixed_exit, |d| crush_phi(d, previous, kk, crush, mu, k))
}

/// The root of f = phi(d_pred - comp f) for `contact_force` and `crush_force`.
#[inline(always)]
fn solve_force(d_pred: f64, comp: f64, force: f64, previous: f64, fixed_exit: bool, phi: impl Fn(f64) -> (f64, f64)) -> f64 {
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
        // An exact root must be accepted before it becomes a bracket endpoint.
        // Otherwise the strict interior safeguard bisects away from the solution.
        if g == 0.0 { return x; }
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

/// One iteration of `contact_force`'s loop for a = 1.5: the same operations in the
/// same order, so the same iterate, bracket and stop. Returns (next iterate, stop).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn newton_step_sqrt(x: f64, lo: &mut f64, hi: &mut f64, dp: f64, cm: f64, dn: f64, kk: f64, mu: f64, fixed_exit: bool, k: f64) -> (f64, bool) {
    let d = dp - cm * x;
    let (v, dv) = contact_phi(d, dn, kk, 1.5, mu, k);
    let g = x - v;
    if g == 0.0 { return (x, true); }
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
/// Runs of contacts on disjoint nodes with the sqrt law and no crush (snare wires, hi-hat rims,
/// sticks) are queued and solved interleaved; any other contact, or one that shares a node with the queue,
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
        let queued = alpha == 1.5 && contacts[jc].crush.is_none();
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
            let f = match &mut c.crush {
                Some(crush) => {
                    // the memory through the secant's lower end, d^(n-1) (see `Crush`)
                    crush.deepest = crush.deepest.max(previous);
                    crush_force(d_pred, comp, c.force, previous, c.law, crush, c.fixed_point_exit, k)
                }
                None => contact_force(d_pred, comp, c.force, previous, c.law, c.fixed_point_exit, k),
            };
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
    /// Set by `load_from_rest`. It lives with the state, so a clone finishes the ramp.
    load_on: Option<LoadOn>,
}

/// A static load coming on from zero: the steps it takes, and each loaded body's full bias.
#[derive(Clone)]
struct LoadOn {
    steps: u64,
    full: Vec<(usize, Vec<f64>)>,
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
        Ok(Instrument { bodies: vec![], cavities: vec![], rads: vec![], links: vec![], contacts: vec![], frictions: vec![], strikes: vec![], drives: vec![], mics: vec![], rate, step_index: 0, gas, materials: None, prepared: false, cq: ContactQueue::default(), load_on: None })
    }
    /// Bring every body's bias on from zero as a raised cosine over `seconds`, not all at once, so
    /// the piece ends that time at rest under its load. A step rings every mode it reaches. This
    /// leaves a mode of frequency f at most 1 / ((2 f T)^2 - 1) of that, and the load rate has no jump.
    pub(crate) fn load_from_rest(&mut self, seconds: f64) {
        let full: Vec<(usize, Vec<f64>)> = self.bodies.iter().enumerate()
            .filter(|(_, b)| b.bias.iter().any(|&x| x != 0.0))
            .map(|(i, b)| (i, b.bias.clone()))
            .collect();
        self.load_on = Some(LoadOn { steps: (seconds * self.rate).round() as u64, full });
    }
    /// Whether the static load is still coming on. `InstrumentImage` carries the bias as a constant,
    /// so a stepper outside this crate must refuse a piece that says so.
    pub fn load_is_coming_on(&self) -> bool {
        self.load_on.as_ref().is_some_and(|load| self.step_index < load.steps)
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
            let material_loss = self.bodies.iter()
                .any(|b| b.viscosity.as_ref().is_some_and(|v| v.rate_bound > 0.0));
            worst.1.push_str(if material_loss {
                " + material bending viscosity"
            } else if self.bodies.iter().any(|b| b.modal.is_some()) {
                " + coil modal loss"
            } else {
                " + gradient loss"
            });
        }
        let radiation_rate2 = 2.0 * self.rate * self.bodies.iter().filter_map(|b| b.head_radiation.as_deref())
            .map(|load| load.rate_bound).fold(0.0_f64, f64::max);
        if radiation_rate2 > 0.0 {
            worst.0 += radiation_rate2;
            worst.1.push_str(" + head radiation");
        }
        let link_rate2 = self.link_admission_rate2(self.rate);
        if link_rate2 > 0.0 {
            worst.0 += link_rate2;
            worst.1.push_str(" + spring/dashpot links");
        }
        worst
    }

    /// Bound the explicit links' K + 2 rate C in mass-symmetric coordinates.
    /// The backward-velocity dashpots in `step` enter the leapfrog stability
    /// condition at the candidate clock, not the clock at construction.
    fn link_admission_rate2(&self, rate: f64) -> f64 {
        type Node = (usize, usize);
        let mut rows = BTreeMap::<Node, BTreeMap<Node, f64>>::new();
        for link in &self.links {
            let q = link.k + 2.0 * rate * link.c;
            if !q.is_finite() { return f64::INFINITY; }
            *rows.entry(link.a).or_default().entry(link.a).or_default() += q;
            if let Some(b) = link.b {
                *rows.entry(b).or_default().entry(b).or_default() += q;
                *rows.entry(link.a).or_default().entry(b).or_default() -= q;
                *rows.entry(b).or_default().entry(link.a).or_default() -= q;
            }
        }
        // Affine mounts are PSD only after assembling their signed pair
        // terms. Taking absolute values before assembly destroys cancellation
        // between those terms and invents extra stiffness/damping.
        let mut bound = 0.0_f64;
        for (&(a, i), row) in &rows {
            let ma = self.bodies[a].mass[i];
            if !ma.is_finite() || ma <= 0.0 { return f64::INFINITY; }
            let mut sum = 0.0;
            for (&(b, j), &q) in row {
                let mb = self.bodies[b].mass[j];
                if !q.is_finite() || !mb.is_finite() || mb <= 0.0 { return f64::INFINITY; }
                sum += if (a, i) == (b, j) { q.abs() / ma }
                    else { q.abs() / ma.sqrt() / mb.sqrt() };
            }
            bound = bound.max(sum);
        }
        bound
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
        // Both explicit viscous forces need dt^2 lambda_K + 2 dt lambda_C < 4.
        // For the variable-mass gradient loss, bound M^-1 L_M by its absolute
        // row sums. The constant-mass shortcut 8/h^2 is not valid on a taper.
        let bound = self.bodies.iter().map(|b| {
            let mut gradient = 0.0_f64;
            if b.s1 > 0.0 {
                for i in 0..b.n {
                    let (mut diagonal, mut off) = (0.0, 0.0);
                    for p in b.lap.ptr[i]..b.lap.ptr[i + 1] {
                        let j = b.lap.col[p];
                        let a = b.lap.val[p];
                        if i == j { diagonal += b.mass[i] * a; } else {
                            diagonal += 0.5 * a * (b.mass[i] - b.mass[j]);
                            off += (0.5 * a * (b.mass[i] + b.mass[j])).abs();
                        }
                    }
                    gradient = gradient.max((diagonal.abs() + off) / b.mass[i]);
                }
                gradient *= 2.0 * b.s1 / (b.h * b.h);
            }
            gradient + b.viscosity.as_ref().map_or(0.0, |v| v.rate_bound)
                + b.modal.as_ref().map_or(0.0, |m| m.rate_bound)
        }).fold(0.0_f64, f64::max);
        2.0 * rate * bound
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
    /// gas/material damping and spring/dashpot links; excludes contact and
    /// large-deflection changes.
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
            // One geometric image for every radiator, including cavity ports.
            // Rigid pressure reflection has gain +1; distance and delay still
            // belong to each source. Preparing these paths allocates no ray list.
            let targets = [Some(m.pos), m.floor.map(|floor| floor.image(m.pos))];
            for (bi, b) in self.bodies.iter().enumerate() {
                if !b.radiates {
                    continue;
                }
                let h2 = b.h * b.h;
                for tp in targets.into_iter().flatten() {
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
                    for tp in targets.into_iter().flatten() {
                        let d = [tp[0] - q.pos[0], tp[1] - q.pos[1], tp[2] - q.pos[2]];
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
        if let Some(load) = self.load_on.as_ref().filter(|load| self.step_index <= load.steps) {
            let g = 0.5 * (1.0 - (PI * self.step_index as f64 / load.steps as f64).cos());
            for (body, full) in &load.full {
                for (x, f) in self.bodies[*body].bias.iter_mut().zip(full) {
                    *x = f * g;
                }
            }
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
            if let Some(crush) = &mut c.crush {
                crush.deepest = 0.0;
            }
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
    /// in apps/instrument-lab/render-gpu. A head's film law, its radiation load and a load still coming
    /// on are not in it, so such a stepper must refuse them (render-gpu's plan does). Exporting changes nothing the steps compute, and
    /// `step` stays the reference that stepper is checked against. It prepares the microphones' taps as the
    /// first step would, so place microphones and bodies before exporting (a later move would no longer
    /// re-prepare them).
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
    /// Cached rest curvatures (xx, yy, xy), one triple per core stencil.
    pub rest_curv: Vec<(f64, f64, f64)>,
    /// Weak annular geometric-strain rows; includes the hole's three affine modes.
    pub membrane_terms: Option<Vec<Vec<MembraneTerm>>>,
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
    /// Source-retaining curvature-rate law and its C=G'G operator. Each node
    /// has `BENDING_ROOT_ROWS_PER_NODE` tensor-root rows. The separate joule
    /// counter is physical centered-velocity loss.
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
            vk: self.vk.as_ref().map(|v| VonKarmanImage { eh: v.eh, nb: v.nb.clone(), l0: v.l0.clone(), rest_curv: v.rest_curv.clone(), membrane_terms: v.annulus.as_ref().map(|a| a.terms.clone()), z0: v.z0.clone(), f: v.f.clone() }),
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

    /// Solve A F = x with the factor every step uses. Rows are core nodes followed
    /// by three affine hole-trace coordinates for an annulus. This is the same order
    /// as the exported stress and weak-load rows. False for an incompatible length.
    /// False (and `x` untouched) for a body without von Karman coupling. Column j of A^-1 is
    /// the solve of the j-th unit vector.
    pub fn airy_solve(&self, x: &mut [f64]) -> bool {
        match &self.vk {
            Some(v) if x.len() == v.f.len() => {
                if let Some(a) = &v.annulus { a.solve(x); } else if let Some(a) = &v.airy { a.solve(x); }
                true
            }
            _ => false,
        }
    }

    /// Original Airy operator, last solved right-hand side and stress function
    /// in the same core-node order as `airy_solve`. The last evaluation may
    /// precede the current nodal state after an instrument step rotates it.
    /// This borrows the actual source system; it does not reconstruct A from
    /// an inverse or rerun/alter the physical trajectory.
    pub fn airy_system(&self) -> Option<(&Csr, &[f64], &[f64])> {
        self.vk.as_ref().map(|v| (&v.matrix, v.q.as_slice(), v.f.as_slice()))
    }
}

// --------------------------------------------------------------------- decimation
const DECIMATION_RADIUS_FRAMES: usize = 16;

fn decimation_kernel(factor: usize) -> Vec<f64> {
    let taps = 2 * DECIMATION_RADIUS_FRAMES * factor + 1;
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
    decimate_from(x, factor, 0)
}

/// Filter the complete trajectory, emitting samples starting at `start`.
/// Earlier samples remain available to the FIR instead of being replaced by
/// silence when removing settling time. Output phase is relative to `start`;
/// the frame count is `(x.len() - start) / factor`, as for decimating the crop.
/// A start at or beyond the end returns no frames. `factor` must be positive.
pub fn decimate_from(x: &[f64], factor: usize, start: usize) -> Vec<f64> {
    assert!(factor > 0, "decimation factor must be positive");
    if start >= x.len() { return Vec::new(); }
    if factor == 1 { return x[start..].to_vec(); }
    let h = decimation_kernel(factor);
    let taps = h.len();
    let half = taps / 2;
    (0..(x.len() - start) / factor)
        .map(|o| {
            let c = start + o * factor;
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

#[cfg(test)]
mod settling_decimation_tests {
    use super::*;

    #[test]
    fn settled_constant_does_not_acquire_a_startup_transient() {
        for factor in [2, 3, 4, 6, 8] {
            let input = vec![1.0; 100 * factor];
            let output = decimate_from(&input, factor, 32 * factor);
            let dc = decimation_kernel(factor).iter().sum::<f64>();
            assert_eq!(output[0].to_bits(), dc.to_bits());
            assert!((decimate(&input[32 * factor..], factor)[0] - dc).abs() > 0.2);
        }
    }

    #[test]
    fn crop_preserves_continuous_filter_history_phase_and_frame_count() {
        for factor in [2, 3, 4, 6, 8] {
            let input: Vec<_> = (0..200 * factor + 1)
                .map(|i| (0.17 * i as f64).sin() + 0.3 * (0.79 * i as f64).cos())
                .collect();
            for start in [0, 1, 5, 32 * factor, 32 * factor + 1] {
                let output = decimate_from(&input, factor, start);
                assert_eq!(output.len(), (input.len() - start) / factor);
                let mut filter = StreamingDecimator::new(factor);
                let delay = DECIMATION_RADIUS_FRAMES * factor;
                let mut expected = Vec::new();
                for i in 0..input.len() + delay {
                    filter.push(input.get(i).copied().unwrap_or(0.0));
                    if i >= start + delay && (i - start - delay) % factor == 0
                        && expected.len() < output.len() {
                        expected.push(filter.value());
                    }
                }
                assert_eq!(output, expected, "factor {factor}, start {start}");
            }
        }
    }

    #[test]
    fn identity_and_empty_crops_keep_their_boundaries() {
        let input = [-0.0, 0.25, -0.5];
        assert_eq!(decimate_from(&input, 1, 1), input[1..]);
        assert_eq!(decimate(&input, 1)[0].to_bits(), (-0.0_f64).to_bits());
        for factor in [1, 2, 6] {
            assert!(decimate_from(&input, factor, input.len()).is_empty());
            assert!(decimate_from(&input, factor, usize::MAX).is_empty());
            assert!(decimate_from(&[], factor, 0).is_empty());
        }
    }
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
mod membrane_compliance_tests;

#[cfg(test)]
#[test]
fn exact_contact_root_is_not_bisected_away() {
    // Hertz contact at d = 1/4: 8 * d^(3/2) = 1 exactly. The predictor
    // d_pred = 5/4 and unit compliance leave that penetration under force 1.
    // This covers both scalar and interleaved solvers without a tolerance
    // relaxation or reliance on the optional rounded-fixed-point policy.
    for fixed_exit in [false, true] {
        let (mut lo, mut hi) = (0.0, 16.0);
        let (force, done) = newton_step_sqrt(
            1.0, &mut lo, &mut hi, 1.25, 1.0, 0.25, 8.0, 0.0,
            fixed_exit, 0.1,
        );
        assert_eq!((force, done), (1.0, true));
        assert_eq!((lo, hi), (0.0, 16.0));
        assert_eq!(contact_force(
            1.25, 1.0, 1.0, 0.25,
            Law::Power { k: 8.0, alpha: 1.5, mu: 0.0 }, fixed_exit, 0.1,
        ), 1.0);
    }
}

#[cfg(test)]
#[test]
fn g0_crush_is_hertz_to_its_pressure_then_holds_it_and_unloads_from_a_set() {
    // The scale of the shaft on a crash edge: K, the Hertz area per approach, hickory's 84 MPa.
    let (kk, area, pressure) = (9.0e7, 4.4e-3, 84.0e6);
    let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-12 * b.abs(), "{a} vs {b}");
    let mut c = Crush::new(area, pressure).unwrap();
    let onset = c.onset(kk);
    // Hertz until its mean pressure is the crush pressure, then that mean pressure on the Hertz area
    close(kk * onset.powf(1.5) / (area * onset), pressure);
    for x in [0.25, 0.5, 1.0] {
        close(c.force(kk, x * onset), kk * (x * onset).powf(1.5));
    }
    for x in [1.5, 4.0, 9.0] {
        close(c.force(kk, x * onset) / (area * x * onset), pressure);
    }
    assert_eq!(c.set(kk), 0.0);
    // crushed to 4 d_y: set 4 - sqrt(4) = 2 d_y, unloading stiffness K 4^(1/4), continuous at the turn
    c.deepest = 4.0 * onset;
    close(c.set(kk), 2.0 * onset);
    close(c.force(kk, 4.0 * onset), pressure * area * 4.0 * onset);
    close(c.force(kk, 3.0 * onset), kk * 2f64.sqrt() * onset.powf(1.5));
    assert_eq!(c.force(kk, 2.0 * onset), 0.0);
    assert_eq!(c.force(kk, -onset), 0.0);
    // at every approach a deeper memory gives no more force: the potential only falls as it crushes
    let curve = |m: f64| -> Vec<f64> {
        let at = Crush { deepest: m * onset, ..c };
        (1..=48).map(|i| at.force(kk, 0.25 * i as f64 * onset)).collect()
    };
    let mut last = curve(0.0);
    for m in [1.0, 1.05, 1.5, 2.0, 4.0, 11.0] {
        let now = curve(m);
        for (i, (a, b)) in now.iter().zip(&last).enumerate() {
            assert!(*a <= *b * (1.0 + 1e-12), "memory {m} d_y gives more force at {} d_y", 0.25 * (i + 1) as f64);
        }
        last = now;
    }
    assert!(Crush::new(0.0, pressure).is_err() && Crush::new(area, f64::NAN).is_err());
    let felt = Contact::new(vec![], vec![], 0.0, Law::Power { k: kk, alpha: 2.5, mu: 0.0 });
    assert!(felt.crushing(c).is_err());
}

#[cfg(test)]
#[test]
fn g1_a_mass_dropped_on_a_crushing_contact_rebounds_as_thornton_and_never_gains_energy() {
    // A mass on a rigid crushing contact, stepped as solve_contacts steps it: the memory through
    // the secant's lower end, then the discrete-gradient solve. Without Hunt-Crossley loss the
    // staggered energy (at the memory of the step) never rises; an elastic drop returns its speed,
    // and a crushing one crushes and rebounds as the continuous law does: loading takes
    // 2/5 K d_y^(5/2) + p A' (d_m^2 - d_y^2) / 2, unloading returns 2/5 p A' d_m sqrt(d_m d_y).
    // What is left is the turn, between samples: about a k^2 / d_m, 3e-4 at the crash's own 288 kHz.
    // A memory through d^n instead was 1.1e-4 too deep at 1 m/s and 20 MHz (scaled by the step,
    // about 7e-3 at 288 kHz).
    let (kk, area, pressure, mass) = (9.0e7, 4.4e-3, 84.0e6, 0.013);
    let law = Law::Power { k: kk, alpha: 1.5, mu: 0.0 };
    let slope = area * pressure;
    let potential = |c: &Crush, d: f64| {
        let (stiff, set, reach) = c.branch(kk);
        let x = (d.min(reach) - set).max(0.0);
        let line = if d > reach { 0.5 * slope * (d * d - reach * reach) } else { 0.0 };
        0.4 * stiff * x * x * x.sqrt() + line
    };
    let cases = [(2.0e7, 1e-6), (2.88e5, 2e-3)].into_iter().flat_map(|rt| [0.05, 1.0, 4.0].map(move |v| (rt, v)));
    for ((rate, tol), speed) in cases {
        let k = 1.0 / rate;
        let mut c = Crush::new(area, pressure).unwrap();
        let onset = c.onset(kk);
        let energy = |c: &Crush, a: f64, b: f64| 0.5 * mass * ((b - a) * rate).powi(2) + 0.5 * (potential(c, a) + potential(c, b));
        let (mut d0, mut d1, mut f) = (-1e-6 - speed * k, -1e-6, 0.0);
        let start = energy(&c, d0, d1);
        let (mut h, mut rise) = (start, 0.0_f64);
        while !(d1 < -1e-6 && d1 < d0) {
            let d_pred = 2.0 * d1 - d0;
            c.deepest = c.deepest.max(d0);
            let before = energy(&c, d0, d1);
            f = crush_force(d_pred, k * k / mass, f, d0, law, &c, false, k);
            let d2 = d_pred - k * k / mass * f;
            let next = energy(&c, d1, d2);
            rise = rise.max(before - h).max(next - before);
            (h, d0, d1) = (next, d1, d2);
        }
        assert!(rise <= 1e-12 * start, "{speed} m/s: energy rose by {rise:e} J");
        let restitution = (d0 - d1) * rate / speed;
        if speed < 0.08 {
            assert!(c.deepest < onset && c.set(kk) == 0.0);
            assert!((restitution - 1.0).abs() < 1e-9, "{speed} m/s: elastic restitution {restitution}");
            continue;
        }
        let deepest = (onset * onset + (mass * speed * speed - 0.8 * kk * onset.powf(2.5)) / slope).sqrt();
        let expected = (0.4 * slope * deepest * (deepest * onset).sqrt() / (0.5 * mass * speed * speed)).sqrt();
        assert!(c.deepest > 10.0 * onset);
        assert!((c.deepest / deepest - 1.0).abs() < tol, "{rate} Hz, {speed} m/s: crushed to {} vs {deepest}", c.deepest);
        assert!((c.set(kk) / (deepest - (deepest * onset).sqrt()) - 1.0).abs() < tol);
        assert!((restitution / expected - 1.0).abs() < tol, "{rate} Hz, {speed} m/s: restitution {restitution} vs {expected}");
    }
}

#[cfg(test)]
#[test]
fn internal_damping_preserves_momentum_and_dissipates_on_unequal_masses() {
    // Two free nodes joined by one viscous edge. The continuum flux gives
    // equal/opposite forces -c(v0-v1), +c(v0-v1), c=2*s1*mean(m)/h^2.
    // The old M L ordering instead gives [-1,4] here: +2 W and +3 N.
    for masses in [[1.0, 4.0], [2.0, 2.0]] {
        for velocity in [[2.0, 1.0], [1.0, 2.0], [3.0, 3.0]] {
            let mut b = Body::string("viscous edge", 2, 1.0, 1.0, 0.0, f64::INFINITY, 0.0, 0.0);
            b.mass.copy_from_slice(&masses);
            b.s0.fill(0.0);
            b.s1 = 0.5;
            b.lap = Csr::from_rows(&[
                BTreeMap::from([(0, 1.0), (1, -1.0)]),
                BTreeMap::from([(0, -1.0), (1, 1.0)]),
            ]);
            let k = 1.0 / 128.0;
            for i in 0..2 { b.wp[i] = -k * velocity[i]; }
            b.sync();
            b.linear_forces();
            b.predict(k);
            let force: Vec<_> = (0..2).map(|i| masses[i] * (b.wn[i] + b.wp[i]) / (k*k)).collect();
            let c = 0.5 * (masses[0] + masses[1]);
            let expected = -c * (velocity[0] - velocity[1]);
            assert!((force[0] - expected).abs() < 1e-12);
            assert!((force[0] + force[1]).abs() < 1e-12, "internal force changed momentum");
            let power = force[0] * velocity[0] + force[1] * velocity[1];
            assert!((power + c * (velocity[0] - velocity[1]).powi(2)).abs() < 1e-12);
            let mut ins = Instrument::new(1.0 / k);
            ins.add(b);
            // The nonzero eigenvalue of this two-node viscous operator is
            // c*(1/m0+1/m1). Admission must enclose it at every trial clock.
            let exact_rate = c * (1.0 / masses[0] + 1.0 / masses[1]);
            let momentum = masses[0] * velocity[0] + masses[1] * velocity[1];
            let exact_relative = (velocity[0] - velocity[1]) * (-exact_rate * 0.25).exp();
            let mut previous_error = f64::INFINITY;
            for rate in [128.0, 256.0, 512.0] {
                assert!(ins.viscosity_admission_rate2(rate) >= 2.0 * rate * exact_rate);
                // Check the actual time evolution against the analytic viscous
                // dimer, not another implementation of the sparse stencil.
                let mut trial = ins.clone();
                trial.rate = rate;
                for i in 0..2 { trial.bodies[0].wp[i] = -velocity[i] / rate; }
                trial.bodies[0].sync();
                let mut last_energy = 0.5 * (masses[0] * velocity[0].powi(2) + masses[1] * velocity[1].powi(2));
                let mut relative = 0.0;
                for _ in 0..(rate / 4.0) as usize {
                    trial.step(&mut []);
                    let b = &trial.bodies[0];
                    let v = [(b.w[0] - b.wp[0]) * rate, (b.w[1] - b.wp[1]) * rate];
                    assert!((masses[0] * v[0] + masses[1] * v[1] - momentum).abs() < 1e-9);
                    let energy = 0.5 * (masses[0] * v[0].powi(2) + masses[1] * v[1].powi(2));
                    assert!(energy <= last_energy + 1e-10, "viscous trajectory gained kinetic energy");
                    last_energy = energy;
                    relative = v[0] - v[1];
                }
                if velocity[0] != velocity[1] {
                    let error = (relative - exact_relative).abs();
                    assert!(error < 0.55 * previous_error, "backward-velocity loss must converge at first order");
                    previous_error = error;
                } else {
                    assert!(relative.abs() < 1e-10, "rigid translation must remain undamped");
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[test]
fn coil_mode_rings_at_its_discrete_frequency_and_decays_at_its_own_loss() {
    // G1 manufactured: a free pinned coil started in one sine mode stays in it, rings at the
    // leapfrog image of omega_j^2 = (T beta^2 + B beta^4)/mu, decays at the sigma it was given
    // for that mode, and its staggered energy falls by exactly the counted modal loss.
    let (n, h, mu, t, bend) = (31usize, 3.5e-3, 2.0e-3, 0.8, 4.8e-4);
    let law = |omega: f64, share: f64| Ok::<f64, String>(2.0 + 3.0e-4 * omega * share);
    let rate = 400_000.0;
    let k = 1.0 / rate;
    for j in [1usize, 6, 23] {
        let mut b = Body::coil("coil", n, h, mu, t, bend, &law, 0.0, 0.0).unwrap();
        let np1 = (n + 1) as f64;
        let phi: Vec<f64> = (0..n).map(|i| (2.0 / np1).sqrt() * (PI * (j * (i + 1)) as f64 / np1).sin()).collect();
        let beta2 = 4.0 / (h * h) * (PI * j as f64 / (2.0 * np1)).sin().powi(2);
        let omega = ((t * beta2 + bend * beta2 * beta2) / mu).sqrt();
        let sigma = law(omega, bend * beta2 / (t + bend * beta2)).unwrap();
        let amp = 1e-4;
        for i in 0..n { b.w[i] = amp * phi[i]; b.wp[i] = amp * phi[i]; }
        b.sync();
        let mut ins = Instrument::new(rate);
        ins.add(b);
        assert!(ins.stability().0 < 0.5);
        let e0 = ins.bodies[0].discrete_energy(k);
        let steps = (0.05 * rate) as usize;
        let (mut q, mut leak) = (Vec::with_capacity(steps), 0.0_f64);
        for _ in 0..steps {
            ins.step(&mut []);
            let w = &ins.bodies[0].w;
            let qj: f64 = w.iter().zip(&phi).map(|(a, p)| a * p).sum();
            let rest: f64 = w.iter().zip(&phi).map(|(a, p)| (a - qj * p).powi(2)).sum::<f64>().sqrt();
            leak = leak.max(rest);
            q.push(qj);
        }
        assert!(leak < 1e-9 * amp, "mode {j} leaked {leak:e} into other modes");
        // frequency from upward zero crossings, against the leapfrog image of omega
        let ups: Vec<f64> = (1..steps).filter(|&s| q[s - 1] < 0.0 && q[s] >= 0.0)
            .map(|s| (s - 1) as f64 + q[s - 1] / (q[s - 1] - q[s])).collect();
        let f = (ups.len() - 1) as f64 / ((ups[ups.len() - 1] - ups[0]) * k);
        let leapfrog = 2.0 / k * (omega * k / 2.0).asin() / (2.0 * PI);
        assert!((f / leapfrog - 1.0).abs() < 2e-3, "mode {j}: {f} Hz vs {leapfrog} Hz");
        // decay from the per-cycle peaks of |q|
        let mut peaks = vec![];
        for s in 1..steps - 1 {
            if q[s].abs() >= q[s - 1].abs() && q[s].abs() > q[s + 1].abs() && q[s].abs() > 0.0 {
                peaks.push((s as f64 * k, q[s].abs().ln()));
            }
        }
        let m = peaks.len() as f64;
        let (st, sy) = peaks.iter().fold((0.0, 0.0), |(a, b), p| (a + p.0, b + p.1));
        let (mt, my) = (st / m, sy / m);
        let slope = peaks.iter().map(|p| (p.0 - mt) * (p.1 - my)).sum::<f64>()
            / peaks.iter().map(|p| (p.0 - mt).powi(2)).sum::<f64>();
        assert!((-slope / sigma - 1.0).abs() < 0.02, "mode {j}: decay {} vs sigma {sigma}", -slope);
        let e1 = ins.bodies[0].discrete_energy(k);
        let lost = ins.bodies[0].modal_loss_dissipated_j();
        assert!(((e0 - e1) - lost).abs() < 1e-9 * e0, "mode {j}: energy fell {} J, counted {lost} J", e0 - e1);
    }
}
