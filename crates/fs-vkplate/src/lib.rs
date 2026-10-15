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
//! Time stepping is leapfrog with a frequency-dependent loss `s0 + s1 beta^2`, set from two
//! (frequency, T60) points. Contacts (stick, felt beater, snare wires, hi-hat plates) are
//! solved IMPLICITLY per contact as monotone scalar equations, so stiff Hertz contacts are
//! stable at audio rates. The enclosed air of a drum couples its two heads (lumped adiabatic
//! cavity); each head's volume velocity sees a small-ka monopole radiation resistance. Sound
//! is the baffled Rayleigh integral of surface acceleration at each microphone.
//!
//! Honest scope: square-grid (staircase) boundaries; moderate rotations (von Karman); loss
//! parameters are declared estimates, not measured; lumped cavity (no cavity modes); baffled
//! radiation (no diffraction around a cymbal's edge).

use std::collections::BTreeMap;
use std::f64::consts::PI;

pub const RHO_AIR: f64 = 1.204;
pub const C_AIR: f64 = 343.0;

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
        for (yi, p) in y[a..b].iter_mut().zip(self.ptr[a..=b].windows(2)) {
            let (cols, vals) = (&self.col[p[0]..p[1]], &self.val[p[0]..p[1]]);
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

pub struct Body {
    pub name: String,
    pub h: f64,
    pub n: usize,
    pub xy: Vec<[f64; 2]>,
    grid_index: BTreeMap<(i64, i64), usize>,
    pub mass: Vec<f64>,
    s0: Vec<f64>,
    s1: f64,
    lap: Csr,
    stiff: Csr,
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
        let (s0, s1) = sigmas(s.t60, &beta2);
        let s0v: Vec<f64> = (0..n).map(|k| s0 + s0_extra[k]).collect();
        let area = PI * (s.radius * s.radius - s.hole * s.hole);
        let t_mean = thick.iter().sum::<f64>() / n as f64;
        let berger = if s.nonlinear == Nonlinear::Berger { s.young * t_mean / (2.0 * area * (1.0 - nu)) } else { 0.0 };
        let vk = if s.nonlinear == Nonlinear::VonKarman {
            if s.edge != Edge::Free {
                return Err(format!("{name}: von Karman coupling is implemented for free edges (stress-free rim)"));
            }
            let z0: Vec<f64> = xy.iter().map(|&p| (s.rest_height)(r_of(p))).collect();
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
            name: name.into(), h, n, xy, grid_index, mass, s0: s0v, s1, lap, stiff, tension: s.tension, berger, vk,
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: true, area, air_mass: 0.0,
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
            name: name.into(), h, n, xy, grid_index: BTreeMap::new(), mass: vec![mass_per_m * h; n],
            s0: vec![6.907755 / t60; n], s1: 0.0, lap: Csr::from_rows(&lrows), stiff: Csr::from_rows(&vec![BTreeMap::new(); n]),
            tension: tension / h, berger: 0.0, vk: None, // tension/h: the chain energy T/2 sum (dw)^2/h
            w: vec![0.0; n], wp: vec![0.0; n], wn: vec![0.0; n], f: vec![0.0; n], bias: vec![0.0; n], tmp: vec![0.0; n], lw: vec![0.0; n], lwp: vec![0.0; n], g: vec![0.0; n], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0,
        }
    }

    /// A lumped rigid mass (stick, beater, hi-hat rod).
    pub fn point(name: &str, mass: f64) -> Body {
        Body {
            name: name.into(), h: 1.0, n: 1, xy: vec![[0.0, 0.0]], grid_index: BTreeMap::new(), mass: vec![mass],
            s0: vec![0.0], s1: 0.0, lap: Csr::from_rows(&[BTreeMap::new()]), stiff: Csr::from_rows(&[BTreeMap::new()]),
            tension: 0.0, berger: 0.0, vk: None,
            w: vec![0.0], wp: vec![0.0], wn: vec![0.0], f: vec![0.0], bias: vec![0.0], tmp: vec![0.0], lw: vec![0.0], lwp: vec![0.0], g: vec![0.0], psi: 0.0, psi_next: 0.0, sav: false, r_prev: 0.0, budget: 0.0, kc: 0.0, cf: vec![],
            origin: [0.0; 3], facing: 1.0, radiates: false, area: 0.0, air_mass: 0.0,
        }
    }

    pub fn node_at(&self, x: f64, y: f64) -> Option<usize> {
        self.grid_index.get(&((x / self.h).round() as i64, (y / self.h).round() as i64)).copied()
    }

    /// Nodes within `radius` of (x, y), equally weighted (a felt face); at least the nearest node.
    pub fn footprint(&self, x: f64, y: f64, radius: f64) -> Vec<(usize, f64)> {
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
    }

    fn predict(&mut self, k: f64) {
        self.coefs(k);
        let n = self.n;
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
        let (mut gy, mut gzg, mut uy, mut uzu, mut gzu) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (cf, w, wp, f, g) = (&self.cf[..n], &self.w[..n], &self.wp[..n], &self.f[..n], &self.g[..n]);
        let (wn, tmp) = (&mut self.wn[..n], &mut self.tmp[..n]);
        let s1 = self.s1 > 0.0;
        for i in 0..n {
            let [ia, cw, cwp, cs1] = cf[i];
            let mut rhs = cw * w[i] - cwp * wp[i] + f[i] + extra;
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
            self.budget = d0 / (2.0 * k) + (self.s1 * ih2 / k * d1).max(0.0);
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

    pub fn loss(&self) -> (f64, f64) {
        (self.s0.iter().cloned().fold(f64::INFINITY, f64::min), self.s1)
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
pub struct Cavity {
    pub heads: Vec<usize>,
    pub volume: f64,
    pub port: Option<Port>,
}

/// A hole in a head (a kick drum's front port): the air plug in it is an inertance
/// rho L_eff / A driven by the cavity pressure; it vents the cavity and radiates as a monopole.
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
pub struct RadDamp {
    pub bodies: Vec<usize>,
    pub resistance: f64,
}

impl RadDamp {
    pub fn monopole(bodies: Vec<usize>, f_ref: f64) -> RadDamp {
        RadDamp { bodies, resistance: RHO_AIR * (2.0 * PI * f_ref).powi(2) / (4.0 * PI * C_AIR) }
    }
}

/// Linear spring/dashpot between two nodes (or one node and ground).
pub struct Link {
    pub a: (usize, usize),
    pub b: Option<(usize, usize)>,
    pub k: f64,
    pub c: f64,
}

#[derive(Clone, Copy, Debug)]
pub enum Law {
    /// f = K d^alpha (1 + mu d_dot): Hertz with Hunt-Crossley loss.
    Power { k: f64, alpha: f64, mu: f64 },
}

/// One-sided contact: d = sum_side sign * weight * w - gap; d > 0 is interpenetration.
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
pub struct Strike {
    pub time: f64,
    pub striker: usize,
    pub contact: usize,
    pub speed: f64,
    pub clear: f64,
    pub done: bool,
}

/// Piecewise-linear external force program on one node.
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
/// One contact's force: the root of g(f) = f - phi(d_pred - comp f), phi(d) = K d^a (1 + mu (d - d_now)/k)
/// clipped at zero force. g is increasing with a root in [0, d_pred / comp]: safeguarded Newton (a
/// bisection step whenever Newton leaves the bracket), from the last force when it lies in the bracket.
/// With `fixed_exit`, a Newton step that lands exactly on the iterate stops the solve there: the
/// correction is below half an ulp, so the iterate is the root to rounding. Without it the bracket
/// update has just put that iterate on the bracket's edge, so the safeguard bisects away from the
/// root and halves back down to it (30-47 iterations instead of about 3).
fn contact_force(d_pred: f64, comp: f64, force: f64, d_now: f64, law: Law, fixed_exit: bool, k: f64) -> f64 {
    let Law::Power { k: kk, alpha, mu } = law;
    let phi = |d: f64| -> (f64, f64) {
        if d <= 0.0 {
            return (0.0, 0.0);
        }
        let pa1 = if alpha == 1.5 { d.sqrt() } else { d.powf(alpha - 1.0) };
        let pa = pa1 * d;
        let damp = 1.0 + mu * (d - d_now) / k;
        let v = kk * pa * damp;
        if v <= 0.0 { (0.0, 0.0) } else { (v, kk * (alpha * pa1 * damp + pa * mu / k)) }
    };
    if d_pred <= 0.0 {
        return 0.0;
    }
    let (mut lo, mut hi) = (0.0, d_pred / comp);
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
    let pa1 = d.sqrt();
    let pa = pa1 * d;
    let damp = 1.0 + mu * (d - dn) / k;
    let v = kk * pa * damp;
    let dv = kk * (1.5 * pa1 * damp + pa * mu / k);
    // phi is (0, 0) at d <= 0 (where v and dv above are not used) and wherever v <= 0
    let zero = (d <= 0.0) | (v <= 0.0);
    let (v, dv) = if zero { (0.0, 0.0) } else { (v, dv) };
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
#[derive(Default)]
struct ContactQueue {
    /// Per queued contact: its index, d_pred, comp, first iterate, upper bracket, d_now, K, mu, fixed-point
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
        let mut comp = 0.0;
        for &(b, i, w) in &c.a {
            d_pred += w * bodies[b].wn[i];
            comp += w * w * bodies[b].compliance(i, k);
        }
        for &(b, i, w) in &c.b {
            d_pred -= w * bodies[b].wn[i];
            comp += w * w * bodies[b].compliance(i, k);
        }
        if !queued || d_pred <= 0.0 {
            // solved now: the queue is empty, or this contact moves nothing
            let f = contact_force(d_pred, comp, c.force, c.d_now, c.law, c.fixed_point_exit, k);
            apply_force(bodies, c, f, d_pred, comp, k);
            continue;
        }
        for &(b, i, _) in c.a.iter().chain(&c.b) {
            q.mark[b][i] = q.batch;
        }
        let hi = d_pred / comp;
        q.j.push(jc);
        q.dp.push(d_pred);
        q.cm.push(comp);
        q.x0.push(if c.force > 0.0 && c.force < hi { c.force } else { 0.5 * hi });
        q.hi.push(hi);
        q.dn.push(c.d_now);
        q.kk.push(kk);
        q.mu.push(mu);
        q.fx.push(c.fixed_point_exit);
    }
    q.flush(bodies, contacts, k);
}

// --------------------------------------------------------------------- instrument
pub struct Instrument {
    pub bodies: Vec<Body>,
    pub cavities: Vec<Cavity>,
    pub rads: Vec<RadDamp>,
    pub links: Vec<Link>,
    pub contacts: Vec<Contact>,
    pub strikes: Vec<Strike>,
    pub drives: Vec<Drive>,
    pub mics: Vec<Mic>,
    pub rate: f64,
    pub step_index: u64,
    prepared: bool,
    cq: ContactQueue,
}

impl Instrument {
    pub fn new(rate: f64) -> Instrument {
        Instrument { bodies: vec![], cavities: vec![], rads: vec![], links: vec![], contacts: vec![], strikes: vec![], drives: vec![], mics: vec![], rate, step_index: 0, prepared: false, cq: ContactQueue::default() }
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
    /// k^2 lambda_max / 4 at the current rate (below 1 is stable without loss or contact).
    pub fn stability(&mut self) -> (f64, String) {
        let k = self.k();
        let (lam, who) = self.max_rate2();
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
                        let p = [b.origin[0] + b.xy[i][0], b.origin[1] + b.xy[i][1], b.origin[2]];
                        let d = [tp[0] - p[0], tp[1] - p[1], tp[2] - p[2]];
                        if d[2] * b.facing <= 0.0 {
                            continue; // the other side of this surface faces away from this path
                        }
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        let delay = r / C_AIR * self.rate;
                        let di = delay.floor() as usize;
                        let fr = delay - di as f64;
                        // w is the OUTWARD displacement of the facing side, so no sign flip
                        let wgt = RHO_AIR * h2 / (2.0 * PI * r) / (k * k);
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
                    let delay = r / C_AIR * self.rate;
                    let di = delay.floor() as usize;
                    let fr = delay - di as f64;
                    let wgt = RHO_AIR / (2.0 * PI * r);
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
            let p = -RHO_AIR * C_AIR * C_AIR / c.volume * (dv + c.port.as_ref().map_or(0.0, |q| q.x));
            for &hb in &c.heads {
                let b = &mut self.bodies[hb];
                let fa = p * b.h * b.h;
                for f in b.f.iter_mut() {
                    *f += fa;
                }
            }
            if let Some(q) = c.port.as_mut() {
                // the air plug: rho L / A du/dt = p - R u (symplectic Euler, far below its limit)
                let p = -RHO_AIR * C_AIR * C_AIR / c.volume * (dv + q.x);
                let du = k * (p - q.resistance * q.u) * q.area / (RHO_AIR * q.length);
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
                let tau = 4.0 * a / (3.0 * PI * C_AIR);
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
}

// --------------------------------------------------------------------- decimation
/// Windowed-sinc low-pass then keep every `factor`-th sample.
pub fn decimate(x: &[f64], factor: usize) -> Vec<f64> {
    if factor == 1 {
        return x.to_vec();
    }
    let taps = 32 * factor + 1;
    let fc = 0.46 / factor as f64; // cycles per input sample: 0.92 x output Nyquist
    let mid = (taps / 2) as f64;
    let h: Vec<f64> = (0..taps)
        .map(|i| {
            let t = i as f64 - mid;
            let s = if t == 0.0 { 2.0 * fc } else { (2.0 * PI * fc * t).sin() / (PI * t) };
            let wdw = 0.42 - 0.5 * (2.0 * PI * i as f64 / (taps - 1) as f64).cos() + 0.08 * (4.0 * PI * i as f64 / (taps - 1) as f64).cos();
            s * wdw
        })
        .collect();
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
