//! Causal forced response of a planar jet to the acoustic transverse
//! velocity at a flue exit (beads pexfp.135/.138/.321, ruling 7376): ONE
//! impulse response from the exit's transverse velocity to the jet
//! centreline displacement at the labium, integrated in TIME from the
//! linearised vorticity dynamics of the jet. It is causal by construction
//! and carries the exit receptivity, every spatial mode, the continuous
//! spectrum and the finite-thickness fall-off together; no single mode is
//! transformed. The response is held after its window, which cuts a slowly
//! decaying tail (2.5-3.9 % at w 1-2 for window 90; CONTRACT invariant 14).
//!
//! Units: lengths in the Bickley half-width `b`, speeds in the centre speed
//! `U0`, time in `b / U0`, Reynolds number `Re = U0 b / nu`.
//!
//! Model.
//! - Base flow `U = sech^2 y`, parallel (quasi-parallel: the base flow's own
//!   viscous spreading over the window is neglected; the disturbance is
//!   viscous).
//! - A uniform transverse acoustic velocity `v_ac(t)` with displacement
//!   `xi = int v_ac dt`. The disturbance is the jet carried with the fluid
//!   (`zeta_p = U'' xi`, which induces no flow) plus a deviation `zeta`, `psi`:
//!   `d_t zeta + U d_x zeta = U'' v + lap(zeta) / Re`, `v = -d_x psi`,
//!   `lap psi = -zeta`.
//! - The flue exit `x = 0` is a wall from which the jet issues undisplaced
//!   and with zero transverse velocity (Fletcher & Thwaites' zero
//!   displacement and zero slope, the Kutta condition): `psi(0, y) = xi U`
//!   and Thom's wall vorticity
//!   `zeta(0, y) = -2 (psi_1 - psi_0 - dx v_ac) / dx^2 - xi U''`.
//! - Laterally unbounded: outside `|y| > L_y` the field is harmonic, and
//!   each streamwise sine mode continues with the exact discrete exterior
//!   decay `r^|j|`, `r + 1/r = 2 - lambda_k dx^2`. Downstream, a wall at
//!   `x = L_x` stays beyond the disturbance front (speed at most `U0`) for
//!   the whole window.
//! - Output: the centreline displacement relative to the nozzle at station
//!   `X`, `eta = xi + eta_dev`, with `d_t eta_dev + d_x eta_dev = v(x, 0)`
//!   and `eta_dev(0) = -xi`.
//! - Opt-in, by name ([`ExitReceptivity::Measured`]): the measured exit
//!   receptivity in place of the Kutta seed. The jet issues displaced by
//!   `eta(0, t) = R v_ac(t)`, in phase, with `R = (h / U_j)` in `b / U0`
//!   units (Terrien, Vergez, de la Cuadra and Fabre, arXiv 1402.3143, Eq. 1,
//!   "according to the empirical model proposed by de la Cuadra"; de la
//!   Cuadra, thesis, Stanford 2005, sec. 2.3.5 and Fig. 2.28: V_ac and eta_0
//!   in phase). Fletcher's `1/omega` and `pi` phase are "not observed" there
//!   (sec. 4.3.3), so no carried term enters: `psi(0, y) = -R v_ac U`, the
//!   exit's own transverse velocity stays zero (the zero-slope half kept),
//!   and the output is `eta = eta_dev` with `eta_dev(0) = R v_ac`. Growth and
//!   convection remain this solve's. The Kutta law is the default and its
//!   arithmetic is unchanged.
//! - Opt-in, by name ([`BaseSpreading`]): the base jet's own streamwise development, quasi-parallel.
//!   At each `x` the base is `U = a(x) sech^2(y / beta(x))`, with `beta = b(x) / b0` and the centre
//!   speed `a = beta^(-1/2)` from the momentum flux `a^2 beta` conserved. The deviation keeps the
//!   equations above with that local `U` and `U''`; the base's own transverse velocity and every
//!   streamwise derivative of the base are dropped, the quasi-parallel truncation. The
//!   jet carried with the fluid is a rigid translation of the base, an exact solution of the full
//!   equations, so the deviation carries no term in `xi`. The centreline deviation moves at the
//!   local centre speed: `d_t eta_dev + a d_x eta_dev = v(x, 0)`. The exit is the
//!   `beta = 1` profile, so its arithmetic is the parallel jet's. Past the station the base is held
//!   at the station's: the jet beyond the labium is outside this response, and holding it keeps the
//!   band finite. [`BaseSpreading::Parallel`] is the default and its arithmetic is unchanged.
//!
//! Numerics: third-order upwind-biased streamwise transport (second-order
//! closures at both ends), the discrete second difference of `U` for `U''`
//! (the carried jet is then an exact discrete solution, so no error grows
//! with `xi`), five-point viscous Laplacian, SSP-RK3, and a
//! Poisson solve by a sine transform in `x` (two rows per complex FFT)
//! with one tridiagonal solve per mode in `y`. Rows are independent in
//! both transforms, so the result is bit-identical for every worker count.
//!
//! The step response `s(t)` (unit `v_ac` from `t = 0`) is recorded over
//! the window and held at its last value afterwards. At audio rate, a mouth
//! velocity that is linear between samples gives the displacement
//! `eta_n = sum_k g_k v_(n-k)`, with hat weights
//! `g_k = (1/D) [int_(kD)^((k+1)D) s - int_((k-1)D)^(kD) s]` in seconds;
//! `g_0` multiplies the current, still unknown, sample.

use crate::AeroacError;
use fs_fft::{C64 as Bin, Fft};
use fs_math::c64::C64;
use fs_math::det;

/// Smallest clearance [b] between the disturbance front at the end of the
/// window and the downstream wall.
pub const FRONT_MARGIN: f64 = 8.0;

/// How the acoustic transverse velocity enters the jet at the flue exit.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ExitReceptivity {
    /// Fletcher & Thwaites' zero displacement and zero slope: the jet issues
    /// undisplaced while the fluid is carried by `xi = int v_ac dt`.
    #[default]
    Kutta,
    /// The measured receptivity (module documentation): the jet issues
    /// displaced by `R v_ac`, in phase, with no carried term; `R` is the exit
    /// displacement per exit transverse velocity [b per v/U0], `h / U_j` in
    /// jet units.
    Measured {
        /// `R` [b/U0].
        displacement_per_velocity: f64,
    },
}

impl ExitReceptivity {
    /// The exit shift `s` (`psi(0, y) = s U`, `eta_dev(0) = -s`), the exit's
    /// imposed transverse velocity and the carried displacement added to the
    /// output, at transverse velocity `vac` and its integral `xi`.
    fn exit_values(self, vac: f64, xi: f64) -> (f64, f64, f64) {
        match self {
            Self::Kutta => (xi, vac, xi),
            Self::Measured { displacement_per_velocity } => {
                (-(displacement_per_velocity * vac), 0.0, 0.0)
            }
        }
    }

    fn validate(self) -> Result<(), AeroacError> {
        match self {
            Self::Kutta => Ok(()),
            Self::Measured { displacement_per_velocity: r } => {
                if r > 0.0 && r.is_finite() {
                    Ok(())
                } else {
                    Err(AeroacError::InvalidParameter {
                        what: "forced jet: measured exit receptivity must be positive and finite",
                    })
                }
            }
        }
    }
}

/// The base jet's streamwise development (module documentation). Widths are in the exit
/// half-width `b0` and `x` in `b0`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BaseSpreading {
    /// `U = sech^2 y` at every `x`.
    #[default]
    Parallel,
    /// The laminar Bickley jet at the solve's own Reynolds number, matched to the exit's width and
    /// centre speed: the centreline momentum balance `u_c du_c/dx = -2 nu u_c / b^2` with
    /// `u_c^2 b` conserved gives `beta = (1 + 6 x / Re)^(2/3)`, a virtual origin `Re / 6`
    /// half-widths upstream of the exit.
    LaminarBickley,
    /// A width growing linearly, `beta = 1 + width_slope x`, the centre speed from the momentum
    /// flux conserved.
    Linear {
        /// `db / dx`, positive.
        width_slope: f64,
    },
}

impl BaseSpreading {
    /// `beta = b(x) / b0` at `x` [b0] for a solve at `reynolds` (`U0 b0 / nu`).
    #[must_use]
    pub fn width_ratio(self, x: f64, reynolds: f64) -> f64 {
        match self {
            Self::Parallel => 1.0,
            Self::LaminarBickley => det::pow(1.0 + 6.0 * x / reynolds, 2.0 / 3.0),
            Self::Linear { width_slope } => 1.0 + width_slope * x,
        }
    }

    /// Virtual origin [b0] upstream of the exit, where the law's width extrapolates to zero.
    #[must_use]
    pub fn virtual_origin(self, reynolds: f64) -> Option<f64> {
        match self {
            Self::Parallel => None,
            Self::LaminarBickley => Some(reynolds / 6.0),
            Self::Linear { width_slope } => Some(1.0 / width_slope),
        }
    }

    fn validate(self, reynolds: f64) -> Result<(), AeroacError> {
        let bad = |what| Err(AeroacError::InvalidParameter { what });
        match self {
            Self::Parallel => Ok(()),
            Self::LaminarBickley if !reynolds.is_finite() => {
                bad("forced jet: laminar Bickley spreading needs a finite Reynolds number")
            }
            Self::LaminarBickley => Ok(()),
            Self::Linear { width_slope } if !(width_slope > 0.0 && width_slope.is_finite()) => {
                bad("forced jet: linear spreading needs a positive finite width slope")
            }
            Self::Linear { .. } => Ok(()),
        }
    }
}

/// The spread base on the grid: `U` and the discrete `U''` per row and column, and the centre
/// speed per column.
struct SpreadBase {
    u: Vec<f64>,
    upp: Vec<f64>,
    centre_speed: Vec<f64>,
}

/// Grid, window and viscosity of one forced-jet computation, in jet units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForcedJetResolution {
    /// Grid spacing in both directions [b].
    pub spacing: f64,
    /// Streamwise cells (a power of two); the box is `cells * spacing` long.
    pub streamwise_cells: usize,
    /// Half-height of the vortical band [b]; the exact exterior closes it.
    pub half_height: f64,
    /// Recorded window [b/U0]; the response is held at its last value after it.
    pub window: f64,
    /// Time step over spacing (the centreline Courant number).
    pub courant: f64,
    /// `U0 b / nu`; `f64::INFINITY` selects the inviscid jet.
    pub reynolds: f64,
    /// Worker threads; results do not depend on the count.
    pub workers: usize,
}

impl ForcedJetResolution {
    /// Box length [b].
    #[must_use]
    pub fn length(&self) -> f64 {
        self.streamwise_cells as f64 * self.spacing
    }

    /// Time step [b/U0].
    #[must_use]
    pub fn time_step(&self) -> f64 {
        self.courant * self.spacing
    }

    fn validate(&self) -> Result<(usize, usize), AeroacError> {
        let bad = |what| Err(AeroacError::InvalidParameter { what });
        if !(self.spacing > 0.0 && self.spacing <= 0.5) {
            return bad("forced jet: spacing must lie in (0, 0.5] half-widths");
        }
        let nx = self.streamwise_cells;
        if !(nx >= 32 && nx.is_power_of_two()) {
            return bad("forced jet: streamwise cells must be a power of two, at least 32");
        }
        let rows = 2.0 * self.half_height / self.spacing;
        let ny = rows.round() as usize;
        if !(self.half_height >= 4.0)
            || (rows - ny as f64).abs() > 1e-9 * rows
            || !ny.is_multiple_of(2)
        {
            return bad(
                "forced jet: half-height must be at least 4 b and an even number of spacings across",
            );
        }
        if !(self.window > 0.0) || !self.window.is_finite() {
            return bad("forced jet: window must be positive and finite");
        }
        if !(self.length() >= self.window + FRONT_MARGIN) {
            return bad("forced jet: box must exceed the window by the front margin");
        }
        if !(self.courant > 0.0 && self.courant <= 1.5) {
            return bad("forced jet: Courant number must lie in (0, 1.5]");
        }
        if !(self.reynolds > 0.0) {
            return bad("forced jet: Reynolds number must be positive (infinity is inviscid)");
        }
        let diffusion = 4.0 * self.time_step() / (self.reynolds * self.spacing * self.spacing);
        if !(diffusion <= 2.0) {
            return bad("forced jet: viscous step exceeds the explicit stability limit");
        }
        if self.workers == 0 {
            return bad("forced jet: at least one worker");
        }
        Ok((nx, ny))
    }
}

/// Bickley profile `U = sech^2 y`.
fn profile(y: f64) -> f64 {
    let t = det::tanh(y);
    1.0 - t * t
}

struct Scratch {
    buf: Vec<Bin>,
    tmp: Vec<Bin>,
    deriv: Vec<f64>,
}

impl Scratch {
    fn new(cols: usize) -> Self {
        let n = 2 * (cols + 1);
        Self {
            buf: vec![Bin::default(); n],
            tmp: vec![Bin::default(); n],
            deriv: vec![0.0; cols],
        }
    }
}

/// What one row's right-hand side needs besides its own data.
struct RowContext<'a> {
    z: &'a [f64],
    u: &'a [f64],
    upp: &'a [f64],
    cols: usize,
    rows: usize,
    h: f64,
    /// The exit's imposed transverse velocity (`v_ac` under Kutta).
    wall_v: f64,
    /// The exit shift `s`, `psi(0, y) = s U` (`xi` under Kutta).
    shift: f64,
    /// `1 / (Re h^2)`, absent for the inviscid jet.
    visc: Option<f64>,
    /// The spread base's `U` and `U''` per row and column; absent for the parallel jet, whose
    /// exit profiles `u`, `upp` hold at every column.
    spread: Option<(&'a [f64], &'a [f64])>,
}

impl RowContext<'_> {
    /// Transport, production and diffusion of row `r` from its streamfunction `prow`.
    fn row(&self, r: usize, prow: &[f64], deriv: &mut [f64], drow: &mut [f64]) {
        let (cols, h) = (self.cols, self.h);
        let zr = &self.z[r * cols..(r + 1) * cols];
        // The exit is the unspread profile under every base.
        let psi0 = self.shift * self.u[r];
        let z0 = -2.0 * (prow[0] - psi0 - h * self.wall_v) / (h * h) - self.shift * self.upp[r];
        upwind(zr, z0, h, deriv);
        let half = 0.5 / h;
        match self.spread {
            None => {
                for c in 0..cols {
                    let left = if c == 0 { psi0 } else { prow[c - 1] };
                    let right = if c + 1 == cols { 0.0 } else { prow[c + 1] };
                    drow[c] = -self.u[r] * deriv[c] - self.upp[r] * (right - left) * half;
                }
            }
            Some((su, supp)) => {
                let (ur, ppr) = (
                    &su[r * cols..(r + 1) * cols],
                    &supp[r * cols..(r + 1) * cols],
                );
                for c in 0..cols {
                    let left = if c == 0 { psi0 } else { prow[c - 1] };
                    let right = if c + 1 == cols { 0.0 } else { prow[c + 1] };
                    drow[c] = -ur[c] * deriv[c] - ppr[c] * (right - left) * half;
                }
            }
        }
        if let Some(visc) = self.visc {
            let below = (r > 0).then(|| &self.z[(r - 1) * cols..r * cols]);
            let above = (r + 1 < self.rows).then(|| &self.z[(r + 1) * cols..(r + 2) * cols]);
            for c in 0..cols {
                let zl = if c == 0 { z0 } else { zr[c - 1] };
                let zrt = if c + 1 == cols { 0.0 } else { zr[c + 1] };
                let zd = below.map_or(0.0, |row| row[c]);
                let zu = above.map_or(0.0, |row| row[c]);
                drow[c] += visc * (zl + zrt + zd + zu - 4.0 * zr[c]);
            }
        }
    }
}

/// Sine transform `F_k = scale sum_(i=1)^(n-1) f_i sin(pi k i / n)` of one
/// or two rows in place, through one complex FFT of the odd extension.
fn sine_rows(
    fft: &Fft,
    n: usize,
    a: &mut [f64],
    b: Option<&mut [f64]>,
    sc: &mut Scratch,
    scale: f64,
) {
    let buf = &mut sc.buf;
    buf[0] = Bin::default();
    buf[n] = Bin::default();
    match &b {
        Some(b) => {
            for i in 1..n {
                buf[i] = Bin::new(a[i - 1], b[i - 1]);
                buf[2 * n - i] = Bin::new(-a[i - 1], -b[i - 1]);
            }
        }
        None => {
            for i in 1..n {
                buf[i] = Bin::new(a[i - 1], 0.0);
                buf[2 * n - i] = Bin::new(-a[i - 1], 0.0);
            }
        }
    }
    fft.forward(buf, &mut sc.tmp);
    // Z_k = -2i F^a_k + 2 F^b_k for the packed rows a + i b.
    for k in 1..n {
        a[k - 1] = -0.5 * scale * buf[k].im;
    }
    if let Some(b) = b {
        for k in 1..n {
            b[k - 1] = 0.5 * scale * buf[k].re;
        }
    }
}

/// Run `tasks` over `workers` scoped threads in contiguous groups, each with its own scratch.
fn run_tasks<T: Send>(
    workers: usize,
    cols: usize,
    tasks: Vec<T>,
    f: &(impl Fn(&mut Scratch, T) + Sync),
) {
    let n = tasks.len();
    let w = workers.min(n).max(1);
    if w == 1 {
        let mut sc = Scratch::new(cols);
        for task in tasks {
            f(&mut sc, task);
        }
        return;
    }
    let mut groups: Vec<Vec<T>> = (0..w).map(|_| Vec::new()).collect();
    for (i, task) in tasks.into_iter().enumerate() {
        groups[i * w / n].push(task);
    }
    std::thread::scope(|s| {
        for group in groups {
            s.spawn(move || {
                let mut sc = Scratch::new(cols);
                for task in group {
                    f(&mut sc, task);
                }
            });
        }
    });
}

/// Upwind-biased `d/dx` at interior nodes `i = 1..n-1` of the row `f`
/// (node `i` stored at `f[i-1]`), with the inflow value `f0` at `i = 0`.
fn upwind(f: &[f64], f0: f64, h: f64, out: &mut [f64]) {
    let m = f.len();
    let (s2, s6) = (0.5 / h, 1.0 / (6.0 * h));
    out[0] = (f[1] - f0) * s2;
    out[1] = (2.0 * f[2] + 3.0 * f[1] - 6.0 * f[0] + f0) * s6;
    for (o, w) in out[2..m - 1].iter_mut().zip(f.windows(4)) {
        *o = (2.0 * w[3] + 3.0 * w[2] - 6.0 * w[1] + w[0]) * s6;
    }
    out[m - 1] = (3.0 * f[m - 1] - 4.0 * f[m - 2] + f[m - 3]) * s2;
}

/// Time-domain solver of the forced linear jet (the module documentation
/// states the model). The forcing closure returns `(v_ac, xi)` at a time.
pub struct ForcedJet {
    res: ForcedJetResolution,
    nx: usize,
    cols: usize,
    rows: usize,
    centre: usize,
    u: Vec<f64>,
    upp: Vec<f64>,
    inv_den: Vec<f64>,
    fft: Fft,
    zeta: Vec<f64>,
    eta: Vec<f64>,
    time: f64,
    xi: f64,
    vac: f64,
    exit: ExitReceptivity,
    spreading: BaseSpreading,
    spread: Option<SpreadBase>,
    stage_z: Vec<f64>,
    stage_e: Vec<f64>,
    dz: Vec<f64>,
    de: Vec<f64>,
    psi: Vec<f64>,
}

impl ForcedJet {
    /// Quiescent jet at `t = 0` under the Kutta exit.
    pub fn new(resolution: ForcedJetResolution) -> Result<Self, AeroacError> {
        Self::with_exit(resolution, ExitReceptivity::Kutta)
    }

    /// Quiescent jet at `t = 0` under the named exit receptivity.
    pub fn with_exit(
        resolution: ForcedJetResolution,
        exit: ExitReceptivity,
    ) -> Result<Self, AeroacError> {
        Self::with_base(resolution, exit, BaseSpreading::Parallel, f64::INFINITY)
    }

    /// Quiescent jet at `t = 0` under the named exit receptivity and base development, the base
    /// held past `held_beyond` [b0] (the station). A spread base needs a band of at least 4 of its
    /// held half-widths and a positive finite hold.
    pub fn with_base(
        resolution: ForcedJetResolution,
        exit: ExitReceptivity,
        spreading: BaseSpreading,
        held_beyond: f64,
    ) -> Result<Self, AeroacError> {
        exit.validate()?;
        spreading.validate(resolution.reynolds)?;
        let (nx, ny) = resolution.validate()?;
        if spreading != BaseSpreading::Parallel {
            if !(held_beyond > 0.0 && held_beyond.is_finite()) {
                return Err(AeroacError::InvalidParameter {
                    what: "forced jet: a spread base needs a positive finite hold station",
                });
            }
            let widest = spreading.width_ratio(held_beyond, resolution.reynolds);
            if !(widest.is_finite() && resolution.half_height >= 4.0 * widest) {
                return Err(AeroacError::InvalidParameter {
                    what: "forced jet: half-height must be at least 4 of the spread base's held half-widths",
                });
            }
        }
        let h = resolution.spacing;
        let cols = nx - 1;
        let rows = ny - 1;
        // Row r sits at y = (r + 1 - ny/2) h: exactly symmetric, the centre on
        // y = 0. U'' is the DISCRETE second difference of U, so the jet carried
        // with the fluid (psi = -xi U, zeta = U'' xi) solves the discrete
        // Poisson problem exactly; the analytic U'' would leave a residual
        // xi h^2 U''''/12 that grows with the acoustic displacement.
        let at = |j: usize| profile((j as f64 - (ny / 2) as f64) * h);
        let u: Vec<f64> = (1..=rows).map(at).collect();
        let upp: Vec<f64> = (1..=rows)
            .map(|j| (at(j + 1) - 2.0 * at(j) + at(j - 1)) / (h * h))
            .collect();
        // Thomas factors per sine mode: unit off-diagonals, diagonal
        // -2 + lambda_k, plus the exterior decay r_k at both edges.
        let mut inv_den = vec![0.0; rows * cols];
        for k in 0..cols {
            let lam = 2.0 * det::cos(core::f64::consts::PI * (k + 1) as f64 / nx as f64) - 2.0;
            let mu = -lam;
            let r_k = 0.5 * (2.0 + mu - det::sqrt(mu * (4.0 + mu)));
            let mut prev = 0.0;
            for r in 0..rows {
                let mut diag = -2.0 + lam;
                if r == 0 {
                    diag += r_k;
                }
                if r == rows - 1 {
                    diag += r_k;
                }
                let den = diag - prev;
                prev = 1.0 / den;
                inv_den[r * cols + k] = prev;
            }
        }
        // The spread base at node c + 1 (x = (c + 1) h), held past the station: U and its
        // discrete second difference across rows, as for the exit profile.
        let spread = (spreading != BaseSpreading::Parallel).then(|| {
            let mut base = SpreadBase {
                u: vec![0.0; rows * cols],
                upp: vec![0.0; rows * cols],
                centre_speed: vec![0.0; cols],
            };
            for c in 0..cols {
                let x = ((c + 1) as f64 * h).min(held_beyond);
                let beta = spreading.width_ratio(x, resolution.reynolds);
                let a = 1.0 / det::sqrt(beta);
                let at = |j: usize| a * profile((j as f64 - (ny / 2) as f64) * h / beta);
                for r in 0..rows {
                    let j = r + 1;
                    base.u[r * cols + c] = at(j);
                    base.upp[r * cols + c] = (at(j + 1) - 2.0 * at(j) + at(j - 1)) / (h * h);
                }
                base.centre_speed[c] = a;
            }
            base
        });
        Ok(Self {
            res: resolution,
            nx,
            cols,
            rows,
            centre: ny / 2 - 1,
            u,
            upp,
            inv_den,
            fft: Fft::new(2 * nx),
            zeta: vec![0.0; rows * cols],
            eta: vec![0.0; cols],
            time: 0.0,
            xi: 0.0,
            vac: 0.0,
            exit,
            spreading,
            spread,
            stage_z: vec![0.0; rows * cols],
            stage_e: vec![0.0; cols],
            dz: vec![0.0; rows * cols],
            de: vec![0.0; cols],
            psi: vec![0.0; rows * cols],
        })
    }

    /// The resolution this solver runs.
    #[must_use]
    pub fn resolution(&self) -> &ForcedJetResolution {
        &self.res
    }

    /// Current time [b/U0].
    #[must_use]
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Centreline deviation `eta_dev` at the interior nodes `x = i dx`, `i = 1..`.
    #[must_use]
    pub fn centreline_deviation(&self) -> &[f64] {
        &self.eta
    }

    /// The exit receptivity this solver runs.
    #[must_use]
    pub fn exit(&self) -> ExitReceptivity {
        self.exit
    }

    /// The base development this solver runs.
    #[must_use]
    pub fn spreading(&self) -> BaseSpreading {
        self.spreading
    }

    /// Centreline displacement relative to the nozzle, `xi + eta_dev(x)` [b]
    /// (`eta_dev(x)` under the measured exit, which carries nothing), by cubic
    /// interpolation between nodes.
    pub fn displacement_at(&self, station: f64) -> Result<f64, AeroacError> {
        let h = self.res.spacing;
        let s = station / h;
        if !(s >= 1.0 && s + 2.0 <= self.cols as f64) {
            return Err(AeroacError::InvalidParameter {
                what: "forced jet: station outside the interpolation range",
            });
        }
        let (shift, _, carried) = self.exit.exit_values(self.vac, self.xi);
        let i = s.floor() as usize;
        let f = s - i as f64;
        let val = |n: usize| if n == 0 { -shift } else { self.eta[n - 1] };
        let w = [
            -f * (f - 1.0) * (f - 2.0) / 6.0,
            (f + 1.0) * (f - 1.0) * (f - 2.0) / 2.0,
            -(f + 1.0) * f * (f - 2.0) / 2.0,
            (f + 1.0) * f * (f - 1.0) / 6.0,
        ];
        let dev = w[0] * val(i - 1) + w[1] * val(i) + w[2] * val(i + 1) + w[3] * val(i + 2);
        Ok(carried + dev)
    }

    /// Advance one SSP-RK3 step of `courant * spacing`.
    pub fn step(&mut self, forcing: &impl Fn(f64) -> (f64, f64)) -> Result<(), AeroacError> {
        let dt = self.res.time_step();
        let t = self.time;
        let mut z = core::mem::take(&mut self.zeta);
        let mut e = core::mem::take(&mut self.eta);
        let mut sz = core::mem::take(&mut self.stage_z);
        let mut se = core::mem::take(&mut self.stage_e);
        self.rhs(&z, &e, forcing(t));
        for (s, (a, d)) in sz.iter_mut().zip(z.iter().zip(&self.dz)) {
            *s = a + dt * d;
        }
        for (s, (a, d)) in se.iter_mut().zip(e.iter().zip(&self.de)) {
            *s = a + dt * d;
        }
        self.rhs(&sz, &se, forcing(t + dt));
        for (s, (a, d)) in sz.iter_mut().zip(z.iter().zip(&self.dz)) {
            *s = 0.75 * a + 0.25 * (*s + dt * d);
        }
        for (s, (a, d)) in se.iter_mut().zip(e.iter().zip(&self.de)) {
            *s = 0.75 * a + 0.25 * (*s + dt * d);
        }
        self.rhs(&sz, &se, forcing(t + 0.5 * dt));
        for (a, (s, d)) in z.iter_mut().zip(sz.iter().zip(&self.dz)) {
            *a = *a / 3.0 + 2.0 / 3.0 * (s + dt * d);
        }
        for (a, (s, d)) in e.iter_mut().zip(se.iter().zip(&self.de)) {
            *a = *a / 3.0 + 2.0 / 3.0 * (s + dt * d);
        }
        self.zeta = z;
        self.eta = e;
        self.stage_z = sz;
        self.stage_e = se;
        self.time = t + dt;
        let (vac, xi) = forcing(self.time);
        self.xi = xi;
        self.vac = vac;
        if !(vac.is_finite() && xi.is_finite() && self.eta.iter().all(|v| v.is_finite())) {
            return Err(AeroacError::NonFinite {
                what: "forced jet: state left the finite domain",
            });
        }
        Ok(())
    }

    /// Right-hand side into `self.dz`, `self.de` (`self.psi` is the work field).
    fn rhs(&mut self, z: &[f64], e: &[f64], (vac, xi): (f64, f64)) {
        // Under Kutta these are (xi, v_ac) themselves, so its arithmetic is unchanged.
        let (shift, wall_v, _) = self.exit.exit_values(vac, xi);
        let (nx, cols, rows) = (self.nx, self.cols, self.rows);
        let h = self.res.spacing;
        let (u, upp, fft) = (&self.u, &self.upp, &self.fft);
        let pair = 2 * cols;
        // Assemble h^2 times the Poisson source (the exit value psi_0 = s U
        // moves to the right side) and transform each row pair.
        {
            let tasks: Vec<(usize, &mut [f64])> = self
                .psi
                .chunks_mut(pair)
                .enumerate()
                .map(|(p, c)| (2 * p, c))
                .collect();
            run_tasks(self.res.workers, cols, tasks, &|sc: &mut Scratch,
                                                       (r0, chunk): (
                usize,
                &mut [f64],
            )| {
                for (dr, row) in chunk.chunks_mut(cols).enumerate() {
                    let r = r0 + dr;
                    for (o, zz) in row.iter_mut().zip(&z[r * cols..(r + 1) * cols]) {
                        *o = -h * h * zz;
                    }
                    row[0] -= shift * u[r];
                }
                let (a, b) = chunk.split_at_mut(cols.min(chunk.len()));
                let b = if b.is_empty() { None } else { Some(b) };
                sine_rows(fft, nx, a, b, sc, 1.0);
            });
        }
        // Tridiagonal solve across rows, one sine mode per column.
        let inv = &self.inv_den;
        let psi = &mut self.psi;
        for (p, f) in psi[..cols].iter_mut().zip(&inv[..cols]) {
            *p *= f;
        }
        for r in 1..rows {
            let (head, tail) = psi.split_at_mut(r * cols);
            let prev = &head[(r - 1) * cols..];
            let fac = &inv[r * cols..(r + 1) * cols];
            for ((p, q), f) in tail[..cols].iter_mut().zip(prev).zip(fac) {
                *p = (*p - q) * f;
            }
        }
        for r in (0..rows - 1).rev() {
            let (head, tail) = psi.split_at_mut((r + 1) * cols);
            let fac = &inv[r * cols..(r + 1) * cols];
            for ((p, q), f) in head[r * cols..].iter_mut().zip(&tail[..cols]).zip(fac) {
                *p -= f * q;
            }
        }
        // Back to rows, then transport, production and diffusion per row.
        {
            let tasks: Vec<(usize, &mut [f64], &mut [f64])> = self
                .psi
                .chunks_mut(pair)
                .zip(self.dz.chunks_mut(pair))
                .enumerate()
                .map(|(p, (c, d))| (2 * p, c, d))
                .collect();
            let inverse = 2.0 / nx as f64;
            let ctx = RowContext {
                z,
                u,
                upp,
                cols,
                rows,
                h,
                wall_v,
                shift,
                visc: self
                    .res
                    .reynolds
                    .is_finite()
                    .then(|| 1.0 / (self.res.reynolds * h * h)),
                spread: self
                    .spread
                    .as_ref()
                    .map(|base| (base.u.as_slice(), base.upp.as_slice())),
            };
            run_tasks(self.res.workers, cols, tasks, &|sc: &mut Scratch,
                                                       (r0, chunk, dchunk): (
                usize,
                &mut [f64],
                &mut [f64],
            )| {
                {
                    let (a, b) = chunk.split_at_mut(cols.min(chunk.len()));
                    let b = if b.is_empty() { None } else { Some(b) };
                    sine_rows(fft, nx, a, b, sc, inverse);
                }
                for (dr, (prow, drow)) in
                    chunk.chunks(cols).zip(dchunk.chunks_mut(cols)).enumerate()
                {
                    ctx.row(r0 + dr, prow, &mut sc.deriv, drow);
                }
            });
        }
        // Centreline displacement deviation: unit speed (the local centre speed of a spread
        // base), inflow -s; the exit streamfunction there is s U(0) = s.
        let prow = &self.psi[self.centre * cols..(self.centre + 1) * cols];
        upwind(e, -shift, h, &mut self.de);
        let half = 0.5 / h;
        for c in 0..cols {
            let left = if c == 0 {
                shift * self.u[self.centre]
            } else {
                prow[c - 1]
            };
            let right = if c + 1 == cols { 0.0 } else { prow[c + 1] };
            self.de[c] = match &self.spread {
                None => -self.de[c] - (right - left) * half,
                Some(base) => -base.centre_speed[c] * self.de[c] - (right - left) * half,
            };
        }
    }
}

/// Step response `s(t)` of the centreline displacement at one station to a
/// unit transverse velocity switched on at `t = 0` [b per U0, at times in b/U0].
#[derive(Debug, Clone)]
pub struct JetStepResponse {
    station: f64,
    time_step: f64,
    values: Vec<f64>,
    integrals: Vec<f64>,
    resolution: ForcedJetResolution,
    exit: ExitReceptivity,
    spreading: BaseSpreading,
}

/// Integrate the forced jet's unit-step response at `station` [b] over the window (Kutta exit).
pub fn bickley_step_response(
    station: f64,
    resolution: &ForcedJetResolution,
) -> Result<JetStepResponse, AeroacError> {
    bickley_step_response_with_exit(station, resolution, ExitReceptivity::Kutta)
}

/// Integrate the forced jet's unit-step response at `station` [b] over the window under the
/// named exit receptivity.
pub fn bickley_step_response_with_exit(
    station: f64,
    resolution: &ForcedJetResolution,
    exit: ExitReceptivity,
) -> Result<JetStepResponse, AeroacError> {
    bickley_step_response_with_base(station, resolution, exit, BaseSpreading::Parallel)
}

/// Integrate the forced jet's unit-step response at `station` [b0] over the window under the
/// named exit receptivity and base development, the base held past the station.
pub fn bickley_step_response_with_base(
    station: f64,
    resolution: &ForcedJetResolution,
    exit: ExitReceptivity,
    spreading: BaseSpreading,
) -> Result<JetStepResponse, AeroacError> {
    let mut jet = ForcedJet::with_base(*resolution, exit, spreading, station)?;
    jet.displacement_at(station)?;
    let dt = resolution.time_step();
    let steps = (resolution.window / dt).ceil() as usize;
    let mut values = Vec::with_capacity(steps + 1);
    values.push(0.0);
    let forcing = |t: f64| (1.0, t);
    for _ in 0..steps {
        jet.step(&forcing)?;
        values.push(jet.displacement_at(station)?);
    }
    let mut integrals = Vec::with_capacity(values.len());
    let mut acc = 0.0;
    integrals.push(0.0);
    for w in values.windows(2) {
        acc += 0.5 * dt * (w[0] + w[1]);
        integrals.push(acc);
    }
    Ok(JetStepResponse {
        station,
        time_step: dt,
        values,
        integrals,
        resolution: *resolution,
        exit,
        spreading,
    })
}

/// `(1 - e^-x) / x` and `(1 - (1 + x) e^-x) / x^2` for imaginary `x = i w h`.
fn segment_factors(x: f64) -> (C64, C64) {
    if x.abs() < 0.5 {
        // Taylor sums: sum (-x)^n / (n+1)! and sum (-x)^n (n+1) / (n+2)!.
        let (mut e1, mut e2) = (C64::new(0.0, 0.0), C64::new(0.0, 0.0));
        let mut power = C64::new(1.0, 0.0);
        let step = C64::new(0.0, -x);
        let mut fact = 1.0;
        for n in 0..20 {
            let f1 = fact * (n + 1) as f64;
            let f2 = f1 * (n + 2) as f64;
            e1 = e1 + power.scale(1.0 / f1);
            e2 = e2 + power.scale((n + 1) as f64 / f2);
            fact = f1;
            power = power * step;
        }
        (e1, e2)
    } else {
        let (s, c) = det::sin_cos(x);
        let ex = C64::new(c, -s);
        let ix = C64::new(0.0, x);
        let one = C64::new(1.0, 0.0);
        let e1 = (one - ex) / ix;
        let e2 = (one - (one + ix) * ex) / (ix * ix);
        (e1, e2)
    }
}

impl JetStepResponse {
    /// Station [b].
    #[must_use]
    pub fn station(&self) -> f64 {
        self.station
    }

    /// Recording interval [b/U0].
    #[must_use]
    pub fn time_step(&self) -> f64 {
        self.time_step
    }

    /// Recorded window [b/U0].
    #[must_use]
    pub fn window(&self) -> f64 {
        (self.values.len() - 1) as f64 * self.time_step
    }

    /// `s(k dt)`, `k = 0..`; `s(0) = 0`.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The resolution the response was integrated at.
    #[must_use]
    pub fn resolution(&self) -> &ForcedJetResolution {
        &self.resolution
    }

    /// The exit receptivity the response was integrated under.
    #[must_use]
    pub fn exit(&self) -> ExitReceptivity {
        self.exit
    }

    /// The base development the response was integrated under.
    #[must_use]
    pub fn spreading(&self) -> BaseSpreading {
        self.spreading
    }

    /// `b(X) / b0` of that base at the station.
    #[must_use]
    pub fn station_width_ratio(&self) -> f64 {
        self.spreading
            .width_ratio(self.station, self.resolution.reynolds)
    }

    /// Zero-frequency gain of the held response, `s(window)`.
    #[must_use]
    pub fn held_gain(&self) -> f64 {
        *self.values.last().expect("at least one value")
    }

    /// `s(tau)`, linear between records, zero before `t = 0`, held after the window.
    #[must_use]
    pub fn value_at(&self, tau: f64) -> f64 {
        if !(tau > 0.0) {
            return 0.0;
        }
        let x = tau / self.time_step;
        let k = x.floor() as usize;
        if k + 1 >= self.values.len() {
            return self.held_gain();
        }
        let f = x - k as f64;
        self.values[k] + f * (self.values[k + 1] - self.values[k])
    }

    /// `int_0^tau s`, exact for the linear interpolant and the held tail.
    #[must_use]
    pub fn integral_at(&self, tau: f64) -> f64 {
        if !(tau > 0.0) {
            return 0.0;
        }
        let dt = self.time_step;
        let x = tau / dt;
        let k = x.floor() as usize;
        let last = self.values.len() - 1;
        if k >= last {
            return self.integrals[last] + self.held_gain() * (tau - last as f64 * dt);
        }
        let d = tau - k as f64 * dt;
        let slope = (self.values[k + 1] - self.values[k]) / dt;
        self.integrals[k] + self.values[k] * d + 0.5 * slope * d * d
    }

    /// Transfer `G(w) = int_0^inf h(t) e^(-i w t) dt` of the held response at
    /// Strouhal `w = omega b / U0`, in b/U0 per unit `v/U0` (a delay is a phase lag).
    pub fn transfer(&self, omega: f64) -> Result<C64, AeroacError> {
        if !(omega > 0.0) || !omega.is_finite() {
            return Err(AeroacError::InvalidParameter {
                what: "forced jet transfer: omega must be positive",
            });
        }
        let dt = self.time_step;
        let (e1, e2) = segment_factors(omega * dt);
        let mut sum = C64::new(0.0, 0.0);
        for k in 0..self.values.len() - 1 {
            let (s, c) = det::sin_cos(omega * k as f64 * dt);
            let rot = C64::new(c, -s);
            let ya = self.values[k];
            let slope = (self.values[k + 1] - ya) / dt;
            sum = sum + rot * (e1.scale(ya * dt) + e2.scale(slope * dt * dt));
        }
        let t = self.window();
        let (s, c) = det::sin_cos(omega * t);
        Ok(C64::new(0.0, omega) * sum + C64::new(c, -s).scale(self.held_gain()))
    }

    /// Physical-time endpoint weights [s] of one accepted convective interval.
    ///
    /// `younger_age` and `advance` are in the response clock [b/U0]; the interval
    /// runs from age younger_age+advance to younger_age over physical_duration_s.
    /// For endpoint velocities (old, new) in m/s, its contribution is
    /// `older_weight * old + newer_weight * new` metres. This returns
    /// (newer_weight, older_weight), including the current unknown endpoint.
    ///
    /// This is exact for this owner's piecewise-linear step response and a
    /// physical velocity linear over the interval. It integrates the derivative
    /// h against both linear hats locally; no large accumulated integrals are
    /// subtracted. The held response tail has zero derivative and contributes
    /// no additional impulse. Each interval retains its OWN physical clock scale.
    ///
    /// # Errors
    /// Nonfinite/negative age, nonpositive advance/duration, or nonfinite weights.
    pub fn interval_hat_weights(
        &self,
        younger_age: f64,
        advance: f64,
        physical_duration_s: f64,
    ) -> Result<(f64, f64), AeroacError> {
        if !younger_age.is_finite() || younger_age < 0.0
            || !advance.is_finite() || advance <= 0.0
            || !physical_duration_s.is_finite() || physical_duration_s <= 0.0
        {
            return Err(AeroacError::InvalidParameter {
                what: "forced jet interval needs a finite nonnegative age and positive advance/duration",
            });
        }
        if younger_age >= self.window() { return Ok((0.0, 0.0)); }
        let spacing = self.time_step;
        let end = (younger_age + advance).min(self.window());
        let mut age = younger_age;
        let mut knot = (age / spacing).floor() as usize;
        let (mut newer, mut older) = (0.0, 0.0);
        while age < end && knot + 1 < self.values.len() {
            let next = (((knot + 1) as f64) * spacing).min(end);
            if next > age {
                let impulse = (self.values[knot + 1] - self.values[knot])
                    * ((next - age) / spacing);
                let fraction = 0.5
                    * ((age - younger_age) / advance + (next - younger_age) / advance);
                newer += impulse * (1.0 - fraction);
                older += impulse * fraction;
                age = next;
            }
            knot += 1;
        }
        let seconds_per_clock = physical_duration_s / advance;
        let weights = (seconds_per_clock * newer, seconds_per_clock * older);
        if !weights.0.is_finite() || !weights.1.is_finite() {
            return Err(AeroacError::NonFinite { what: "forced jet interval weights" });
        }
        Ok(weights)
    }

    /// Audio-rate hat weights `g_k` [s], `k = 0..`, for a jet of centre speed
    /// `centre_speed_m_s` and Bickley half-width `half_width_m` sampled every
    /// `sample_period_s`: `eta_n = sum_k g_k v_(n-k)` in metres for mouth
    /// velocities in m/s, with `g_0` the current-sample weight. The weights
    /// end one sample after the window; their sum is the held gain in seconds.
    pub fn hat_weights(
        &self,
        centre_speed_m_s: f64,
        half_width_m: f64,
        sample_period_s: f64,
    ) -> Result<Vec<f64>, AeroacError> {
        self.leading_hat_weights(centre_speed_m_s, half_width_m, sample_period_s, usize::MAX)
            .map(|(_, weights)| weights)
    }

    /// The number of [`Self::hat_weights`] and the first `count` of them (all when there are
    /// fewer), bit for bit. A caller that reads only the leading lags of a slow jet's long
    /// response need not compute the rest.
    pub fn leading_hat_weights(
        &self,
        centre_speed_m_s: f64,
        half_width_m: f64,
        sample_period_s: f64,
        count: usize,
    ) -> Result<(usize, Vec<f64>), AeroacError> {
        if !(centre_speed_m_s > 0.0 && half_width_m > 0.0 && sample_period_s > 0.0)
            || !(centre_speed_m_s * half_width_m * sample_period_s).is_finite()
        {
            return Err(AeroacError::InvalidParameter {
                what: "forced jet weights: speed, width and period must be positive",
            });
        }
        let c = half_width_m / centre_speed_m_s;
        let d = sample_period_s / c;
        let all = ((self.window() / d).ceil() as usize).saturating_add(2);
        let taps = all.min(count);
        let mut out = Vec::with_capacity(taps);
        for k in 0..taps {
            let kf = k as f64;
            let second = self.integral_at((kf + 1.0) * d) - 2.0 * self.integral_at(kf * d)
                + self.integral_at((kf - 1.0) * d);
            out.push(c / d * second);
        }
        if !out.iter().all(|g| g.is_finite()) {
            return Err(AeroacError::NonFinite {
                what: "forced jet weights",
            });
        }
        Ok((all, out))
    }
}
