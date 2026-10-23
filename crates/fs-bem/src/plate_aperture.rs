//! Potential-flow (Laplace) conductance of apertures cut through a rigid plate.
//!
//! The plate fills `0 <= z <= t` and is an infinite rigid baffle. Its
//! apertures are a raster of open square cells of side `cell`, cut straight
//! through it. The air is incompressible and irrotational, with potential 1
//! far below the plate and 0 far above. The conductance `K = Q / delta-phi`
//! [m] is the volume flux per unit potential difference: a thin circular hole
//! of radius `a` has Rayleigh's `K = 2a`, and the opening's acoustic mass is
//! `rho / K`.
//!
//! Each raster cell is split `m x m` into cells of side `h = cell / m`.
//! - Half-spaces: each is exact through its Neumann Green's function
//!   `1/(2 pi r)`. On a face, `phi = phi_far -/+ H w`, where `w` is the upward
//!   velocity and `H` is P0 collocation of `1/(2 pi r)` over the face cells.
//!   The cell integrals are exact: closed form up to eight cells away, and
//!   8x8 Gauss-Legendre beyond. `H` is applied as a convolution on a padded
//!   FFT grid, so it equals the dense collocation matrix to rounding.
//! - Duct (`t > 0`): cell-centred finite volumes on the voxel prism, in
//!   `L = max(1, round(t / h))` layers. Walls are Neumann, and each face
//!   potential sits half a layer beyond its end voxels. The layer operator
//!   is diagonal in the DST-II basis, so a duct solve is `L` shifted
//!   cross-section solves, each by Jacobi-preconditioned CG.
//! - Coupling: with `u = [w_b; -w_t]` and the duct's Dirichlet-to-Neumann
//!   map `Lambda`, `(I + Lambda H2) u = Lambda e`, where `e = [1; 0]`. The
//!   operator is self-adjoint and positive definite in the `H2` inner
//!   product, and CG runs in that inner product.
//! - Thin plate (`t = 0`): `2 H w = 1`, by CG.
//!
//! `K` converges at about first order in `h`, because the velocity is
//! singular at the rim (`d^-1/2` in a thin plate). Every result reports its
//! discretization, so a caller ladders `m` and extrapolates. Nothing is
//! fitted. Not represented: plate curvature, a finite baffle, a closed cavity
//! on either side, viscosity and compressibility.

use std::fmt;

use fs_exec::{CancelGate, TilePool};
use fs_fft::{C64, FftNd};
use fs_math::det;

/// An aperture raster: open square cells of side `cell` [m]. Cell `(i, j)`
/// covers `[i, i + 1) x [j, j + 1)` in units of `cell`.
#[derive(Debug, Clone, PartialEq)]
pub struct ApertureRaster {
    /// Cell side [m].
    pub cell: f64,
    /// Open cells, each listed once, in any order.
    pub open: Vec<(i64, i64)>,
}

/// Declared discretization, solver tolerance and budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateApertureOptions {
    /// Plate thickness [m]. Zero is a thin plate.
    pub thickness: f64,
    /// Each raster cell is split `subdivision x subdivision`. The duct takes
    /// `max(1, round(thickness / h))` layers at the resulting side `h`.
    pub subdivision: usize,
    /// Relative Euclidean residual of the outer solve. Each cross-section
    /// solve of the duct runs to a hundredth of it.
    pub tolerance: f64,
    /// Outer CG iteration limit.
    pub max_iterations: usize,
    /// Iteration limit of each cross-section CG solve.
    pub max_inner_iterations: usize,
    /// More refined face cells than this are refused.
    pub max_cells: usize,
    /// More padded FFT grid points than this are refused.
    pub max_fft_points: usize,
    /// More duct voxels (cells times layers) than this are refused.
    pub max_voxels: usize,
    /// Workers for the FFTs and the duct's layer modes. Results are bitwise
    /// identical at every count.
    pub workers: usize,
}

impl Default for PlateApertureOptions {
    fn default() -> Self {
        Self {
            thickness: 0.0,
            subdivision: 1,
            tolerance: 1e-10,
            max_iterations: 1000,
            max_inner_iterations: 20_000,
            max_cells: 1 << 20,
            max_fft_points: 1 << 24,
            max_voxels: 1 << 24,
            workers: 1,
        }
    }
}

/// A solved conductance and its discretization.
#[derive(Debug, Clone, PartialEq)]
pub struct PlateApertureConductance {
    /// `K = Q / delta-phi` [m]: the flux into the bottom face.
    pub conductance: f64,
    /// The flux out of the top face [m]. It equals `conductance` to the
    /// duct solves' tolerance (exactly, for a thin plate).
    pub top_flux: f64,
    /// Refined cell side `h` [m].
    pub cell_side: f64,
    /// Refined face cells.
    pub cells: usize,
    /// Duct layers (zero for a thin plate).
    pub layers: usize,
    /// Padded FFT grid points.
    pub fft_points: usize,
    /// Outer CG iterations.
    pub iterations: usize,
    /// True relative residual of the outer system, recomputed after the solve.
    pub residual: f64,
    /// Cross-section CG iterations, summed over every duct solve.
    pub inner_iterations: usize,
}

/// Why a conductance was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum PlateApertureError {
    /// A declared input is outside its domain.
    Invalid(&'static str),
    /// A declared budget is smaller than the problem needs.
    Budget {
        /// The budgeted quantity.
        what: &'static str,
        /// What the problem needs.
        required: usize,
        /// What was declared.
        maximum: usize,
    },
    /// A CG solve did not reach its tolerance.
    NotConverged {
        /// Which solve.
        stage: &'static str,
        /// Iterations spent.
        iterations: usize,
        /// Relative residual reached.
        residual: f64,
    },
    /// The FFT executor failed.
    Executor(String),
}

impl fmt::Display for PlateApertureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "invalid plate aperture input: {what}"),
            Self::Budget {
                what,
                required,
                maximum,
            } => write!(f, "{what}: {required} exceeds the declared {maximum}"),
            Self::NotConverged {
                stage,
                iterations,
                residual,
            } => write!(
                f,
                "{stage} CG stopped at relative residual {residual:e} after {iterations} iterations"
            ),
            Self::Executor(message) => write!(f, "FFT executor: {message}"),
        }
    }
}

impl std::error::Error for PlateApertureError {}

type Result<T> = std::result::Result<T, PlateApertureError>;

/// Solve the conductance of `raster` through a plate under `options`.
///
/// # Errors
/// [`PlateApertureError::Invalid`] for an empty raster, a repeated cell or a
/// non-physical option; [`PlateApertureError::Budget`] when a declared
/// budget is exceeded; [`PlateApertureError::NotConverged`] when a CG solve
/// (or the recomputed outer residual, at ten times the tolerance) misses.
pub fn plate_aperture_conductance(
    raster: &ApertureRaster,
    options: &PlateApertureOptions,
) -> Result<PlateApertureConductance> {
    let o = options;
    if !(raster.cell.is_finite() && raster.cell > 0.0) {
        return Err(PlateApertureError::Invalid(
            "cell must be finite and positive",
        ));
    }
    if !(o.thickness.is_finite() && o.thickness >= 0.0) {
        return Err(PlateApertureError::Invalid(
            "thickness must be finite and non-negative",
        ));
    }
    if o.subdivision == 0 || o.workers == 0 || o.max_iterations == 0 {
        return Err(PlateApertureError::Invalid(
            "subdivision, workers and max_iterations must be positive",
        ));
    }
    if !(o.tolerance.is_finite() && o.tolerance > 0.0 && o.tolerance < 1.0) {
        return Err(PlateApertureError::Invalid("tolerance must lie in (0, 1)"));
    }
    if raster.open.is_empty() {
        return Err(PlateApertureError::Invalid("the raster has no open cell"));
    }
    // Indices up to 2^40 refined by up to 2^20 stay inside i64.
    if o.subdivision > 1 << 20 {
        return Err(PlateApertureError::Invalid("subdivision above 2^20"));
    }
    let m = o.subdivision;
    let cells_needed = raster.open.len().saturating_mul(m.saturating_mul(m));
    if cells_needed > o.max_cells {
        return Err(PlateApertureError::Budget {
            what: "refined face cells",
            required: cells_needed,
            maximum: o.max_cells,
        });
    }
    let mut open = raster.open.clone();
    open.sort_unstable();
    if open.windows(2).any(|w| w[0] == w[1]) {
        return Err(PlateApertureError::Invalid("an open cell is listed twice"));
    }
    if open
        .iter()
        .any(|&(i, j)| i.unsigned_abs().max(j.unsigned_abs()) > 1 << 40)
    {
        return Err(PlateApertureError::Invalid(
            "cell index magnitude above 2^40",
        ));
    }
    let mi = m as i64;
    let mut cells = Vec::with_capacity(cells_needed);
    for &(i, j) in &open {
        for a in 0..mi {
            for b in 0..mi {
                cells.push((i * mi + a, j * mi + b));
            }
        }
    }
    let (mut i0, mut i1, mut j0, mut j1) = (i64::MAX, i64::MIN, i64::MAX, i64::MIN);
    for &(i, j) in &cells {
        i0 = i0.min(i);
        i1 = i1.max(i);
        j0 = j0.min(j);
        j1 = j1.max(j);
    }
    let nx = (i1 - i0 + 1) as u64;
    let ny = (j1 - j0 + 1) as u64;
    // A linear convolution of n cells needs a pad of at least 2n - 1.
    let pad = |n: u64| {
        n.checked_mul(2)
            .and_then(|v| (v - 1).checked_next_power_of_two())
    };
    let grid = match (pad(nx), pad(ny)) {
        (Some(px), Some(py)) => px.checked_mul(py).map(|points| (px, py, points)),
        _ => None,
    };
    let (px, py, fft_points) = match grid {
        Some((px, py, points)) if points <= o.max_fft_points as u64 => {
            (px as usize, py as usize, points as usize)
        }
        _ => {
            return Err(PlateApertureError::Budget {
                what: "padded FFT grid points",
                required: grid.map_or(usize::MAX, |g| usize::try_from(g.2).unwrap_or(usize::MAX)),
                maximum: o.max_fft_points,
            });
        }
    };
    let (nx, ny) = (nx as usize, ny as usize);
    let h = raster.cell / m as f64;
    let n = cells.len();
    let layers = if o.thickness == 0.0 {
        0
    } else {
        let layers = (o.thickness / h).round().max(1.0);
        let voxels = layers * n as f64;
        if voxels > o.max_voxels as f64 {
            return Err(PlateApertureError::Budget {
                what: "duct voxels",
                required: voxels.min(usize::MAX as f64) as usize,
                maximum: o.max_voxels,
            });
        }
        layers as usize
    };
    let local: Vec<(usize, usize)> = cells
        .iter()
        .map(|&(i, j)| ((i - i0) as usize, (j - j0) as usize))
        .collect();
    let face = FaceOperator::new(&local, nx, ny, px, py, h, o.workers)?;
    if layers == 0 {
        return thin(&face, n, h, fft_points, o);
    }
    let duct = Duct::new(&local, nx, ny, h, o.thickness, layers, o);
    thick(&face, &duct, n, h, fft_points, o)
}

fn thin(
    face: &FaceOperator,
    n: usize,
    h: f64,
    fft_points: usize,
    o: &PlateApertureOptions,
) -> Result<PlateApertureConductance> {
    let mut grid = vec![C64::default(); face.fft.total()];
    let zeros = vec![0.0; n];
    let mut sink = vec![0.0; n];
    let mut apply = |v: &[f64], out: &mut [f64]| -> Result<()> {
        face.apply(v, &zeros, out, &mut sink, &mut grid)?;
        out.iter_mut().for_each(|x| *x *= 2.0);
        Ok(())
    };
    let b = vec![1.0; n];
    let b_norm = norm(&b);
    let mut x = vec![0.0; n];
    let mut r = b.clone();
    let mut p = r.clone();
    let mut ap = vec![0.0; n];
    let mut rr = dot(&r, &r);
    let mut iterations = 0;
    let mut converged = false;
    for it in 1..=o.max_iterations {
        iterations = it;
        apply(&p, &mut ap)?;
        let alpha = rr / dot(&p, &ap);
        axpy(alpha, &p, &mut x);
        axpy(-alpha, &ap, &mut r);
        let rr_new = dot(&r, &r);
        if det::sqrt(rr_new) <= o.tolerance * b_norm {
            converged = true;
            break;
        }
        let beta = rr_new / rr;
        rr = rr_new;
        for (pi, ri) in p.iter_mut().zip(&r) {
            *pi = ri + beta * *pi;
        }
    }
    apply(&x, &mut ap)?;
    let residual = norm(&b.iter().zip(&ap).map(|(b, a)| b - a).collect::<Vec<_>>()) / b_norm;
    if !converged || residual > 10.0 * o.tolerance {
        return Err(PlateApertureError::NotConverged {
            stage: "thin plate",
            iterations,
            residual,
        });
    }
    let conductance = h * h * x.iter().sum::<f64>();
    Ok(PlateApertureConductance {
        conductance,
        top_flux: conductance,
        cell_side: h,
        cells: n,
        layers: 0,
        fft_points,
        iterations,
        residual,
        inner_iterations: 0,
    })
}

fn thick(
    face: &FaceOperator,
    duct: &Duct,
    n: usize,
    h: f64,
    fft_points: usize,
    o: &PlateApertureOptions,
) -> Result<PlateApertureConductance> {
    let mut grid = vec![C64::default(); face.fft.total()];
    let mut inner = 0;
    // H2 applied to [a; b] (both faces in one complex convolution).
    let h2 = |v: &[f64], out: &mut [f64], grid: &mut [C64]| -> Result<()> {
        let (a, b) = v.split_at(n);
        let (oa, ob) = out.split_at_mut(n);
        face.apply(a, b, oa, ob, grid)
    };
    let lambda = |v: &[f64], out: &mut [f64], inner: &mut usize| -> Result<()> {
        *inner += duct.dirichlet_to_neumann(v, out)?;
        Ok(())
    };
    let mut e = vec![0.0; 2 * n];
    e[..n].fill(1.0);
    let mut b = vec![0.0; 2 * n];
    lambda(&e, &mut b, &mut inner)?;
    let b_norm = norm(&b);
    let mut x = vec![0.0; 2 * n];
    let mut r = b.clone();
    let mut hr = vec![0.0; 2 * n];
    h2(&r, &mut hr, &mut grid)?;
    let mut rho = dot(&r, &hr);
    let mut p = r.clone();
    let mut hp = hr.clone();
    let mut tp = vec![0.0; 2 * n];
    let mut iterations = 0;
    let mut converged = false;
    for it in 1..=o.max_iterations {
        iterations = it;
        lambda(&hp, &mut tp, &mut inner)?;
        axpy(1.0, &p, &mut tp);
        let alpha = rho / dot(&hp, &tp);
        axpy(alpha, &p, &mut x);
        axpy(-alpha, &tp, &mut r);
        if norm(&r) <= o.tolerance * b_norm {
            converged = true;
            break;
        }
        h2(&r, &mut hr, &mut grid)?;
        let rho_new = dot(&r, &hr);
        let beta = rho_new / rho;
        rho = rho_new;
        for i in 0..2 * n {
            p[i] = r[i] + beta * p[i];
            hp[i] = hr[i] + beta * hp[i];
        }
    }
    h2(&x, &mut hp, &mut grid)?;
    lambda(&hp, &mut tp, &mut inner)?;
    let residual = norm(&(0..2 * n).map(|i| b[i] - x[i] - tp[i]).collect::<Vec<_>>()) / b_norm;
    if !converged || residual > 10.0 * o.tolerance {
        return Err(PlateApertureError::NotConverged {
            stage: "coupled duct",
            iterations,
            residual,
        });
    }
    let area = h * h;
    Ok(PlateApertureConductance {
        conductance: area * x[..n].iter().sum::<f64>(),
        top_flux: -area * x[n..].iter().sum::<f64>(),
        cell_side: h,
        cells: n,
        layers: duct.layers,
        fft_points,
        iterations,
        residual,
        inner_iterations: inner,
    })
}

// Eight-point Gauss-Legendre rule mapped to [0, 1] (the table of
// helmholtz.rs, pinned there).
const GL8_X: [f64; 8] = [
    0.019_855_071_751_231_884,
    0.101_666_761_293_186_64,
    0.237_233_795_041_835_5,
    0.408_282_678_752_175_1,
    0.591_717_321_247_824_9,
    0.762_766_204_958_164_5,
    0.898_333_238_706_813_4,
    0.980_144_928_248_768_1,
];
const GL8_W: [f64; 8] = [
    0.050_614_268_145_188_13,
    0.111_190_517_226_687_24,
    0.156_853_322_938_943_66,
    0.181_341_891_689_181,
    0.181_341_891_689_181,
    0.156_853_322_938_943_66,
    0.111_190_517_226_687_24,
    0.050_614_268_145_188_13,
];

/// Offsets up to this many cells use the closed form; farther ones the 8x8
/// rule, whose error there is far below rounding while the closed form's
/// corner cancellation grows with the distance.
const NEAR_CELLS: u64 = 8;

/// `integral over the unit square centred at (dx, dy) of 1/r`, seen from the
/// origin.
fn unit_cell_potential(dx: u64, dy: u64) -> f64 {
    if dx.max(dy) <= NEAR_CELLS {
        unit_cell_closed_form(dx as f64, dy as f64)
    } else {
        unit_cell_quadrature(dx as f64, dy as f64)
    }
}

/// `F(x2, y2) - F(x1, y2) - F(x2, y1) + F(x1, y1)` with the antiderivative
/// `F(x, y) = x asinh(y/|x|) + y asinh(x/|y|)` of `1/r`.
fn unit_cell_closed_form(dx: f64, dy: f64) -> f64 {
    let f = |x: f64, y: f64| {
        let mut out = 0.0;
        if x != 0.0 {
            out += x * asinh(y / x.abs());
        }
        if y != 0.0 {
            out += y * asinh(x / y.abs());
        }
        out
    };
    f(dx + 0.5, dy + 0.5) - f(dx - 0.5, dy + 0.5) - f(dx + 0.5, dy - 0.5) + f(dx - 0.5, dy - 0.5)
}

fn unit_cell_quadrature(dx: f64, dy: f64) -> f64 {
    let mut sum = 0.0;
    for (u, wu) in GL8_X.iter().zip(GL8_W) {
        let x = dx + u - 0.5;
        let mut row = 0.0;
        for (v, wv) in GL8_X.iter().zip(GL8_W) {
            let y = dy + v - 0.5;
            row += wv / det::sqrt(x * x + y * y);
        }
        sum += wu * row;
    }
    sum
}

fn asinh(s: f64) -> f64 {
    let a = s.abs();
    let v = det::ln(a + det::sqrt(a * a + 1.0));
    if s < 0.0 { -v } else { v }
}

/// `H`: P0 collocation of `1/(2 pi r)` over the face cells, as a linear
/// convolution on a zero-padded power-of-two grid.
struct FaceOperator {
    fft: FftNd,
    kernel: Vec<C64>,
    position: Vec<usize>,
    pool: Option<TilePool>,
    gate: CancelGate,
}

impl FaceOperator {
    fn new(
        cells: &[(usize, usize)],
        nx: usize,
        ny: usize,
        px: usize,
        py: usize,
        h: f64,
        workers: usize,
    ) -> Result<Self> {
        let scale = h / (2.0 * core::f64::consts::PI);
        let table: Vec<f64> = (0..nx * ny)
            .map(|k| scale * unit_cell_potential((k / ny) as u64, (k % ny) as u64))
            .collect();
        let fft = FftNd::new(&[px, py]);
        let mut kernel = vec![C64::default(); px * py];
        // A grid residue is an offset of either sign; |offset| < n is all a
        // linear convolution of n cells reaches when the pad is >= 2n - 1.
        for ix in 0..px {
            let ox = if 2 * ix < px { ix } else { px - ix };
            if ox >= nx {
                continue;
            }
            for iy in 0..py {
                let oy = if 2 * iy < py { iy } else { py - iy };
                if oy < ny {
                    kernel[ix * py + iy] = C64::new(table[ox * ny + oy], 0.0);
                }
            }
        }
        let pool = (workers > 1).then(|| TilePool::for_host(workers, 0x504C_4154_4150_4552));
        let gate = CancelGate::new_clock_free();
        let mut face = Self {
            fft,
            kernel: Vec::new(),
            position: cells.iter().map(|&(i, j)| i * py + j).collect(),
            pool,
            gate,
        };
        face.forward(&mut kernel)?;
        face.kernel = kernel;
        Ok(face)
    }

    fn forward(&self, data: &mut [C64]) -> Result<()> {
        match &self.pool {
            None => {
                self.fft.forward(data);
                Ok(())
            }
            Some(pool) => self
                .fft
                .forward_pooled(data, pool, &self.gate)
                .map_err(|e| PlateApertureError::Executor(format!("{e:?}"))),
        }
    }

    fn inverse(&self, data: &mut [C64]) -> Result<()> {
        match &self.pool {
            None => {
                self.fft.inverse(data);
                Ok(())
            }
            Some(pool) => self
                .fft
                .inverse_pooled(data, pool, &self.gate)
                .map_err(|e| PlateApertureError::Executor(format!("{e:?}"))),
        }
    }

    /// `out_a = H a` and `out_b = H b`: the kernel is real and even, so the
    /// real and imaginary parts convolve independently.
    fn apply(
        &self,
        a: &[f64],
        b: &[f64],
        out_a: &mut [f64],
        out_b: &mut [f64],
        grid: &mut [C64],
    ) -> Result<()> {
        grid.fill(C64::default());
        for ((&p, &va), &vb) in self.position.iter().zip(a).zip(b) {
            grid[p] = C64::new(va, vb);
        }
        self.forward(grid)?;
        for (g, k) in grid.iter_mut().zip(&self.kernel) {
            *g = C64::new(g.re * k.re - g.im * k.im, g.re * k.im + g.im * k.re);
        }
        self.inverse(grid)?;
        for ((&p, oa), ob) in self.position.iter().zip(out_a).zip(out_b) {
            *oa = grid[p].re;
            *ob = grid[p].im;
        }
        Ok(())
    }
}

/// One DST-II layer mode `v_k[l] = sin(k pi (l + 1/2) / L)` of the duct.
struct LayerMode {
    /// `(h/dz)^2 * 4 sin^2(k pi / 2L)`: the mode's layer eigenvalue, over dz.
    shift: f64,
    /// `v_k[0]`.
    bottom: f64,
    /// `v_k[L - 1]`.
    top: f64,
    /// `1 / |v_k|^2`.
    inverse_norm: f64,
}

/// Cell-centred finite volumes on the voxel prism, solved mode by mode.
struct Duct {
    neighbours: Vec<[u32; 4]>,
    degree: Vec<f64>,
    modes: Vec<LayerMode>,
    layers: usize,
    beta: f64,
    dz: f64,
    tolerance: f64,
    max_iterations: usize,
    workers: usize,
}

impl Duct {
    fn new(
        cells: &[(usize, usize)],
        nx: usize,
        ny: usize,
        h: f64,
        thickness: f64,
        layers: usize,
        o: &PlateApertureOptions,
    ) -> Self {
        let mut index = vec![u32::MAX; nx * ny];
        for (c, &(i, j)) in cells.iter().enumerate() {
            index[i * ny + j] = c as u32;
        }
        let at = |i: usize, j: usize| index[i * ny + j];
        let neighbours: Vec<[u32; 4]> = cells
            .iter()
            .map(|&(i, j)| {
                [
                    if i + 1 < nx { at(i + 1, j) } else { u32::MAX },
                    if i > 0 { at(i - 1, j) } else { u32::MAX },
                    if j + 1 < ny { at(i, j + 1) } else { u32::MAX },
                    if j > 0 { at(i, j - 1) } else { u32::MAX },
                ]
            })
            .collect();
        let degree = neighbours
            .iter()
            .map(|nb| nb.iter().filter(|&&q| q != u32::MAX).count() as f64)
            .collect();
        let dz = thickness / layers as f64;
        let beta = (h / dz) * (h / dz);
        let modes = (1..=layers)
            .map(|k| {
                let s = det::sin(k as f64 * core::f64::consts::PI / (2 * layers) as f64);
                LayerMode {
                    shift: 4.0 * beta * s * s,
                    bottom: s,
                    top: if k % 2 == 1 { s } else { -s },
                    inverse_norm: if k == layers {
                        1.0 / layers as f64
                    } else {
                        2.0 / layers as f64
                    },
                }
            })
            .collect();
        Self {
            neighbours,
            degree,
            modes,
            layers,
            beta,
            dz,
            tolerance: 0.01 * o.tolerance,
            max_iterations: o.max_inner_iterations,
            workers: o.workers.min(layers),
        }
    }

    /// `[phi_b; phi_t]` to `[w_b; -w_t]`, the outward normal derivative on
    /// both faces. Returns the cross-section CG iterations spent.
    fn dirichlet_to_neumann(&self, phi: &[f64], out: &mut [f64]) -> Result<usize> {
        let n = self.degree.len();
        let (pb, pt) = phi.split_at(n);
        let mut solutions = vec![0.0; n * self.layers];
        let per = self.layers.div_ceil(self.workers);
        let outcomes: Vec<Result<usize>> = std::thread::scope(|scope| {
            let handles: Vec<_> = solutions
                .chunks_mut(n * per)
                .zip(self.modes.chunks(per))
                .map(|(block, modes)| {
                    scope.spawn(move || -> Result<usize> {
                        let mut spent = 0;
                        let mut rhs = vec![0.0; n];
                        for (x, mode) in block.chunks_mut(n).zip(modes) {
                            let weight = 2.0 * self.beta * mode.inverse_norm;
                            for c in 0..n {
                                rhs[c] = weight * (mode.bottom * pb[c] + mode.top * pt[c]);
                            }
                            spent += self.shifted_solve(mode.shift, &rhs, x)?;
                        }
                        Ok(spent)
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("duct mode worker panicked"))
                .collect()
        });
        let mut spent = 0;
        for outcome in outcomes {
            spent += outcome?;
        }
        let (wb, wt) = out.split_at_mut(n);
        let gain = 2.0 / self.dz;
        for c in 0..n {
            let (mut bottom, mut top) = (0.0, 0.0);
            for (k, mode) in self.modes.iter().enumerate() {
                let x = solutions[k * n + c];
                bottom += mode.bottom * x;
                top += mode.top * x;
            }
            wb[c] = gain * (pb[c] - bottom);
            wt[c] = gain * (pt[c] - top);
        }
        Ok(spent)
    }

    /// `(L_perp + shift) x = rhs` by Jacobi-preconditioned CG, from zero.
    fn shifted_solve(&self, shift: f64, rhs: &[f64], x: &mut [f64]) -> Result<usize> {
        let n = rhs.len();
        x.fill(0.0);
        let b_norm = norm(rhs);
        if b_norm == 0.0 {
            return Ok(0);
        }
        let diagonal: Vec<f64> = self.degree.iter().map(|d| d + shift).collect();
        let mut r = rhs.to_vec();
        let mut z: Vec<f64> = r.iter().zip(&diagonal).map(|(r, d)| r / d).collect();
        let mut p = z.clone();
        let mut q = vec![0.0; n];
        let mut rz = dot(&r, &z);
        for it in 1..=self.max_iterations {
            for c in 0..n {
                let mut s = diagonal[c] * p[c];
                for &nb in &self.neighbours[c] {
                    if nb != u32::MAX {
                        s -= p[nb as usize];
                    }
                }
                q[c] = s;
            }
            let alpha = rz / dot(&p, &q);
            axpy(alpha, &p, x);
            axpy(-alpha, &q, &mut r);
            if norm(&r) <= self.tolerance * b_norm {
                return Ok(it);
            }
            for c in 0..n {
                z[c] = r[c] / diagonal[c];
            }
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for c in 0..n {
                p[c] = z[c] + beta * p[c];
            }
        }
        Err(PlateApertureError::NotConverged {
            stage: "duct cross-section",
            iterations: self.max_iterations,
            residual: norm(&r) / b_norm,
        })
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(a: &[f64]) -> f64 {
    det::sqrt(dot(a, a))
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (y, x) in y.iter_mut().zip(x) {
        *y += alpha * x;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options(thickness: f64, subdivision: usize) -> PlateApertureOptions {
        PlateApertureOptions {
            thickness,
            subdivision,
            tolerance: 1e-12,
            ..PlateApertureOptions::default()
        }
    }

    fn l_mask() -> ApertureRaster {
        ApertureRaster {
            cell: 1e-3,
            open: vec![(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (3, 2), (5, 3)],
        }
    }

    #[test]
    fn unit_cell_potential_closed_form_meets_quadrature() {
        // The self integral of 1/r over the unit square from its centre.
        let own = unit_cell_closed_form(0.0, 0.0);
        let exact = 4.0 * (1.0 + 2.0_f64.sqrt()).ln();
        assert!(
            (own - exact).abs() <= 4.0 * f64::EPSILON * exact,
            "{own} {exact}"
        );
        // Both forms agree where the switch is made and beyond it.
        for &(dx, dy) in &[(8.0, 0.0), (8.0, 5.0), (8.0, 8.0), (9.0, 3.0), (12.0, 12.0)] {
            let closed = unit_cell_closed_form(dx, dy);
            let quad = unit_cell_quadrature(dx, dy);
            assert!(
                ((closed - quad) / quad).abs() < 1e-13,
                "{dx} {dy}: {closed} {quad}"
            );
        }
        // Far away the square is a point charge plus its isotropic quadrupole
        // term 1/(24 R^3); the next term is O(R^-5).
        let r: f64 = 40.0;
        let far = unit_cell_quadrature(r, 0.0);
        let asymptotic = 1.0 / r + 1.0 / (24.0 * r * r * r);
        assert!(
            ((far - asymptotic) / far).abs() < 1e-8,
            "{far} {asymptotic}"
        );
    }

    #[test]
    fn face_convolution_equals_dense_collocation() {
        let raster = l_mask();
        let cells: Vec<(usize, usize)> = raster
            .open
            .iter()
            .map(|&(i, j)| (i as usize, j as usize))
            .collect();
        let (nx, ny) = (6, 4);
        let h = 0.5e-3;
        let a: Vec<f64> = (0..cells.len())
            .map(|k| 1.0 + 0.37 * k as f64 - 0.05 * (k * k) as f64)
            .collect();
        let b: Vec<f64> = (0..cells.len()).map(|k| (k as f64 - 2.5).abs()).collect();
        let dense = |v: &[f64], c: usize| -> f64 {
            let (i, j) = cells[c];
            let mut s = 0.0;
            for (q, &(qi, qj)) in cells.iter().enumerate() {
                s += unit_cell_potential(i.abs_diff(qi) as u64, j.abs_diff(qj) as u64) * v[q];
            }
            s * h / (2.0 * core::f64::consts::PI)
        };
        let da: Vec<f64> = (0..cells.len()).map(|c| dense(&a, c)).collect();
        let db: Vec<f64> = (0..cells.len()).map(|c| dense(&b, c)).collect();
        let scale = da.iter().chain(&db).fold(0.0_f64, |m, v| m.max(v.abs()));
        let mut results = Vec::new();
        for workers in [1, 3] {
            let face = FaceOperator::new(&cells, nx, ny, 16, 8, h, workers).unwrap();
            let mut grid = vec![C64::default(); 128];
            let (mut oa, mut ob) = (vec![0.0; cells.len()], vec![0.0; cells.len()]);
            face.apply(&a, &b, &mut oa, &mut ob, &mut grid).unwrap();
            for c in 0..cells.len() {
                assert!(
                    (oa[c] - da[c]).abs() <= 1e-13 * scale,
                    "{c}: {} {}",
                    oa[c],
                    da[c]
                );
                assert!(
                    (ob[c] - db[c]).abs() <= 1e-13 * scale,
                    "{c}: {} {}",
                    ob[c],
                    db[c]
                );
            }
            results.push((oa, ob));
        }
        assert_eq!(results[0], results[1], "pooled FFT changed the bits");
    }

    /// Dense Gaussian elimination with partial pivoting (test oracle only).
    fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
        let n = b.len();
        for col in 0..n {
            let pivot = (col..n)
                .max_by(|&p, &q| a[p][col].abs().total_cmp(&a[q][col].abs()))
                .unwrap();
            a.swap(col, pivot);
            b.swap(col, pivot);
            for row in col + 1..n {
                let f = a[row][col] / a[col][col];
                for k in col..n {
                    a[row][k] -= f * a[col][k];
                }
                b[row] -= f * b[col];
            }
        }
        let mut x = vec![0.0; n];
        for row in (0..n).rev() {
            let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
            x[row] = (b[row] - s) / a[row][row];
        }
        x
    }

    #[test]
    fn duct_map_is_the_finite_volume_schur_complement() {
        let raster = l_mask();
        let cells: Vec<(usize, usize)> = raster
            .open
            .iter()
            .map(|&(i, j)| (i as usize, j as usize))
            .collect();
        let n = cells.len();
        let (h, thickness, layers) = (0.5e-3, 1.4e-3, 3);
        let o = options(thickness, 1);
        let duct = Duct::new(&cells, 6, 4, h, thickness, layers, &o);
        let dz = thickness / layers as f64;
        // The 3-D finite-volume system, voxel (c, l) at c * L + l.
        let size = n * layers;
        let mut a = vec![vec![0.0; size]; size];
        let (c_in, c_z, c_face) = (dz, h * h / dz, 2.0 * h * h / dz);
        for c in 0..n {
            for &nb in &duct.neighbours[c] {
                if nb == u32::MAX {
                    continue;
                }
                for l in 0..layers {
                    a[c * layers + l][c * layers + l] += c_in;
                    a[c * layers + l][nb as usize * layers + l] -= c_in;
                }
            }
            for l in 0..layers - 1 {
                let (p, q) = (c * layers + l, c * layers + l + 1);
                a[p][p] += c_z;
                a[q][q] += c_z;
                a[p][q] -= c_z;
                a[q][p] -= c_z;
            }
            a[c * layers][c * layers] += c_face;
            a[c * layers + layers - 1][c * layers + layers - 1] += c_face;
        }
        let phi: Vec<f64> = (0..2 * n)
            .map(|k| 0.3 + 0.11 * k as f64 - 0.02 * (k * k % 7) as f64)
            .collect();
        let mut rhs = vec![0.0; size];
        for c in 0..n {
            rhs[c * layers] += c_face * phi[c];
            rhs[c * layers + layers - 1] += c_face * phi[n + c];
        }
        let v = dense_solve(a, rhs);
        let mut out = vec![0.0; 2 * n];
        duct.dirichlet_to_neumann(&phi, &mut out).unwrap();
        for c in 0..n {
            let wb = c_face * (phi[c] - v[c * layers]) / (h * h);
            let wt_neg = -c_face * (v[c * layers + layers - 1] - phi[n + c]) / (h * h);
            assert!(
                (out[c] - wb).abs() <= 1e-10 * wb.abs().max(1.0),
                "{c}: {} {wb}",
                out[c]
            );
            assert!(
                (out[n + c] - wt_neg).abs() <= 1e-10 * wt_neg.abs().max(1.0),
                "{c}: {} {wt_neg}",
                out[n + c]
            );
        }
        // A uniform drop is the linear profile: w = 1/t on every cell, on both faces.
        let mut drop = vec![0.0; 2 * n];
        drop[..n].fill(1.0);
        duct.dirichlet_to_neumann(&drop, &mut out).unwrap();
        for c in 0..n {
            assert!(
                (out[c] * thickness - 1.0).abs() < 1e-11,
                "{}",
                out[c] * thickness
            );
            assert!(
                (out[n + c] * thickness + 1.0).abs() < 1e-11,
                "{}",
                out[n + c] * thickness
            );
        }
        // Lambda is symmetric, and a uniform potential carries no flux.
        let other: Vec<f64> = (0..2 * n).map(|k| (k as f64 * 0.7).sin()).collect();
        let (mut lp, mut lo) = (vec![0.0; 2 * n], vec![0.0; 2 * n]);
        duct.dirichlet_to_neumann(&phi, &mut lp).unwrap();
        duct.dirichlet_to_neumann(&other, &mut lo).unwrap();
        let (x, y) = (dot(&other, &lp), dot(&phi, &lo));
        assert!((x - y).abs() <= 1e-10 * x.abs().max(y.abs()), "{x} {y}");
        duct.dirichlet_to_neumann(&vec![1.0; 2 * n], &mut out)
            .unwrap();
        assert!(out.iter().all(|w| w.abs() < 1e-9), "{out:?}");
    }

    #[test]
    fn workers_and_cell_order_leave_the_bits() {
        let raster = l_mask();
        let mut reversed = raster.clone();
        reversed.open.reverse();
        for thickness in [0.0, 1.2e-3] {
            let one = plate_aperture_conductance(&raster, &options(thickness, 2)).unwrap();
            let mut pooled = options(thickness, 2);
            pooled.workers = 3;
            let three = plate_aperture_conductance(&raster, &pooled).unwrap();
            let shuffled = plate_aperture_conductance(&reversed, &options(thickness, 2)).unwrap();
            assert_eq!(one, three, "workers changed the result");
            assert_eq!(one, shuffled, "the listed cell order changed the result");
            assert!(((one.top_flux - one.conductance) / one.conductance).abs() < 1e-10);
        }
    }

    #[test]
    fn refusals() {
        let raster = l_mask();
        let ok = options(1e-3, 1);
        let err = |r: &ApertureRaster, o: &PlateApertureOptions| {
            plate_aperture_conductance(r, o).unwrap_err()
        };
        let mut bad = raster.clone();
        bad.cell = 0.0;
        assert!(matches!(err(&bad, &ok), PlateApertureError::Invalid(_)));
        bad.cell = f64::NAN;
        assert!(matches!(err(&bad, &ok), PlateApertureError::Invalid(_)));
        let mut twice = raster.clone();
        twice.open.push((3, 2));
        assert!(matches!(err(&twice, &ok), PlateApertureError::Invalid(_)));
        let empty = ApertureRaster {
            cell: 1e-3,
            open: Vec::new(),
        };
        assert!(matches!(err(&empty, &ok), PlateApertureError::Invalid(_)));
        for o in [
            PlateApertureOptions {
                thickness: -1e-3,
                ..ok
            },
            PlateApertureOptions {
                thickness: f64::INFINITY,
                ..ok
            },
            PlateApertureOptions {
                subdivision: 0,
                ..ok
            },
            PlateApertureOptions { workers: 0, ..ok },
            PlateApertureOptions {
                tolerance: 0.0,
                ..ok
            },
            PlateApertureOptions {
                max_iterations: 0,
                ..ok
            },
        ] {
            assert!(
                matches!(err(&raster, &o), PlateApertureError::Invalid(_)),
                "{o:?}"
            );
        }
        let cells = PlateApertureOptions {
            subdivision: 2,
            max_cells: 27,
            ..ok
        };
        assert_eq!(
            err(&raster, &cells),
            PlateApertureError::Budget {
                what: "refined face cells",
                required: 28,
                maximum: 27
            }
        );
        let fft = PlateApertureOptions {
            max_fft_points: 127,
            ..ok
        };
        assert_eq!(
            err(&raster, &fft),
            PlateApertureError::Budget {
                what: "padded FFT grid points",
                required: 128,
                maximum: 127
            }
        );
        // 1 mm through 1 mm cells: one layer.
        let voxels = PlateApertureOptions {
            max_voxels: 6,
            ..ok
        };
        assert_eq!(
            err(&raster, &voxels),
            PlateApertureError::Budget {
                what: "duct voxels",
                required: 7,
                maximum: 6
            }
        );
        let short = PlateApertureOptions {
            max_iterations: 1,
            ..ok
        };
        assert!(matches!(
            err(&raster, &short),
            PlateApertureError::NotConverged {
                stage: "coupled duct",
                iterations: 1,
                ..
            }
        ));
        let thin_short = PlateApertureOptions {
            thickness: 0.0,
            max_iterations: 1,
            ..ok
        };
        assert!(matches!(
            err(&raster, &thin_short),
            PlateApertureError::NotConverged {
                stage: "thin plate",
                iterations: 1,
                ..
            }
        ));
        let inner_short = PlateApertureOptions {
            thickness: 4e-3,
            max_inner_iterations: 1,
            ..ok
        };
        assert!(matches!(
            err(&raster, &inner_short),
            PlateApertureError::NotConverged {
                stage: "duct cross-section",
                ..
            }
        ));
    }
}
