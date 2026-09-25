//! Quasistatic, incompressible Reynolds-film pressure on a finite-volume graph.
//!
//! This is an orthogonal image to the compressible 1-D gas-film time solver.
//! It has no clock or pressure history. At every caller-supplied configuration,
//! `L(h) p = A B v`, and the positive resisting force is `B^T A p`.
//! Thus mechanical power equals `sum_edges G (p_i-p_j)^2 >= 0`.
//! A caller must evaluate this port *inside* its implicit mechanical equation.
//!
//! Poiseuille conductance is `width * h_+^3 / (12 * viscosity * length)`.
//! A nonpositive CHANNEL aperture closes that passage exactly; a collapsed
//! VOLUME cell, disconnected trapped pocket, or exceeded validity limit refuses.
//! There is no gap floor, pressure clipping, diagonal jitter, or hidden leakage.
//! Reservoir pressures are zero gauge. Suction is retained within the supplied
//! pressure limit. Compressibility, gas inertia, thermal effects, slip, moving
//! topology and fluid acoustic radiation are not part of this image.
/// Compressible mass storage and transport on the same geometric pressure graph.
pub mod isothermal;

use core::fmt;

/// Explicit dense pressure-island ceiling, not a real-time performance claim.
pub const MAX_CELLS: usize = 64;
/// Maximum number of generalized mechanical ports.
pub const MAX_PORTS: usize = 256;
const MAX_CHANNELS: usize = 256;

/// Linearized physical separation `reference_m - closure dot q`.
#[derive(Debug, Clone)]
pub struct GapPort {
    /// Reference separation [m].
    pub reference_m: f64,
    /// Work-conjugate closure row in the caller's coordinate basis.
    pub closure: Vec<f64>,
}

/// One retained fluid volume, with pressure work conjugate to swept volume.
#[derive(Debug, Clone)]
pub struct FilmCell {
    /// Fixed projected cell area [m^2].
    pub area_m2: f64,
    /// Separation and normal motion at the cell's quadrature station.
    pub gap: GapPort,
}

/// One Poiseuille passage between two cells or to zero-gauge ambient pressure.
#[derive(Debug, Clone)]
pub struct FilmChannel {
    /// Upstream cell (the law itself is reciprocal).
    pub from: usize,
    /// Adjacent cell, or `None` for an ambient boundary.
    pub to: Option<usize>,
    /// Face width [m].
    pub width_m: f64,
    /// Distance between pressure samples [m].
    pub length_m: f64,
    /// Explicit hydraulic aperture at the connecting face.
    pub gap: GapPort,
}

/// Caller-declared domain of this low-pressure, thin-gap viscous image.
#[derive(Debug, Clone, Copy)]
pub struct FilmLimits {
    /// Positive minimum retained VOLUME-cell gap [m]. Never a replacement gap.
    pub minimum_cell_gap_m: f64,
    /// Maximum admitted cell/channel gap [m].
    pub maximum_gap_m: f64,
    /// Maximum absolute GAUGE pressure [Pa]; choose small relative to ambient.
    pub maximum_pressure_pa: f64,
}

/// A refusal leaves all caller-owned output slices untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilmError(pub &'static str);
impl fmt::Display for FilmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for FilmError {}

/// Endpoint port diagnostics, not a time-integrated mechanical energy receipt.
#[derive(Debug, Clone, Copy)]
pub struct FilmReport {
    /// Smallest retained cell separation [m].
    pub minimum_gap_m: f64,
    /// Largest absolute cell gauge pressure [Pa].
    pub peak_pressure_pa: f64,
    /// Nonnegative viscous power at the supplied work-conjugate rate [W].
    pub dissipated_power_w: f64,
}

/// Immutable film geometry. Evaluation and analytic tangent use bounded stack
/// scratch, never mutate pressure history, and allocate only at construction.
#[derive(Debug, Clone)]
pub struct ResistiveFilm {
    cells: Vec<FilmCell>,
    channels: Vec<FilmChannel>,
    factors: Vec<f64>,
    ports: usize,
    limits: FilmLimits,
}

fn fail<T>(s: &'static str) -> Result<T, FilmError> {
    Err(FilmError(s))
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

impl ResistiveFilm {
    /// Admit supplied SI geometry and viscosity. Cell area and the transpose
    /// force projection are never normalized independently.
    pub fn new(
        cells: Vec<FilmCell>,
        channels: Vec<FilmChannel>,
        ports: usize,
        viscosity_pa_s: f64,
        limits: FilmLimits,
    ) -> Result<Self, FilmError> {
        if ports == 0
            || ports > MAX_PORTS
            || cells.is_empty()
            || cells.len() > MAX_CELLS
            || channels.is_empty()
            || channels.len() > MAX_CHANNELS
            || !viscosity_pa_s.is_finite()
            || viscosity_pa_s <= 0.0
            || !limits.minimum_cell_gap_m.is_finite()
            || limits.minimum_cell_gap_m <= 0.0
            || !limits.maximum_gap_m.is_finite()
            || limits.maximum_gap_m <= limits.minimum_cell_gap_m
            || !limits.maximum_pressure_pa.is_finite()
            || limits.maximum_pressure_pa <= 0.0
        {
            return fail("invalid film dimensions, viscosity, or validity limits");
        }
        let valid_gap = |g: &GapPort| {
            g.reference_m.is_finite()
                && g.closure.len() == ports
                && g.closure.iter().all(|x| x.is_finite())
        };
        for cell in &cells {
            if !cell.area_m2.is_finite() || cell.area_m2 <= 0.0 || !valid_gap(&cell.gap) {
                return fail("invalid film cell area or closure port");
            }
        }
        let mut factors = Vec::with_capacity(channels.len());
        for edge in &channels {
            if edge.from >= cells.len()
                || edge.to.is_some_and(|j| j >= cells.len() || j == edge.from)
                || !valid_gap(&edge.gap)
                || !edge.width_m.is_finite()
                || edge.width_m <= 0.0
                || !edge.length_m.is_finite()
                || edge.length_m <= 0.0
            {
                return fail("invalid film channel geometry or connectivity");
            }
            let factor = (edge.width_m / edge.length_m) / (12.0 * viscosity_pa_s);
            if !factor.is_finite() || factor <= 0.0 {
                return fail("film conductance scale overflow");
            }
            factors.push(factor);
        }
        // Reference gaps need not be current gaps: attachments may be admitted
        // after motion. The mechanical owner checks its ACTUAL state on attach.
        Ok(Self {
            cells,
            channels,
            factors,
            ports,
            limits,
        })
    }

    /// Number of work-conjugate mechanical coordinates.
    pub fn port_count(&self) -> usize {
        self.ports
    }
    /// Number of retained pressure cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Check cell and channel gap validity without advancing or solving pressure.
    pub fn validate_configuration(&self, q: &[f64]) -> Result<f64, FilmError> {
        if q.len() != self.ports || q.iter().any(|x| !x.is_finite()) {
            return fail("invalid film configuration");
        }
        let mut minimum = f64::INFINITY;
        for c in &self.cells {
            let h = c.gap.reference_m - dot(&c.gap.closure, q);
            if !h.is_finite() || h < self.limits.minimum_cell_gap_m || h > self.limits.maximum_gap_m
            {
                return fail("film volume cell outside declared gap domain");
            }
            minimum = minimum.min(h);
        }
        for edge in &self.channels {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            if !h.is_finite() || h > self.limits.maximum_gap_m {
                return fail("film channel outside declared gap domain");
            }
        }
        Ok(minimum)
    }

    fn solve(&self, q: &[f64], v: &[f64]) -> Result<Solution, FilmError> {
        let minimum = self.validate_configuration(q)?;
        if v.len() != self.ports || v.iter().any(|x| !x.is_finite()) {
            return fail("invalid film velocity/effort");
        }
        let n = self.cells.len();
        let mut s = Solution::new(n, minimum);
        for (k, (edge, factor)) in self.channels.iter().zip(&self.factors).enumerate() {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            let g = if h > 0.0 { factor * h * h * h } else { 0.0 };
            if !g.is_finite() {
                return fail("nonfinite film conductance");
            }
            s.conductance[k] = g;
            let i = edge.from;
            s.lower[i * n + i] += g;
            if let Some(j) = edge.to {
                s.lower[j * n + j] += g;
                s.lower[i * n + j] -= g;
                s.lower[j * n + i] -= g;
            }
        }
        for (i, cell) in self.cells.iter().enumerate() {
            s.pressure[i] = cell.area_m2 * dot(&cell.gap.closure, v);
        }
        s.factor()?;
        solve_factored(&s.lower, &s.scale, n, &mut s.pressure)?;
        if s.pressure[..n]
            .iter()
            .any(|p| p.abs() > self.limits.maximum_pressure_pa)
        {
            return fail("film pressure exceeds declared incompressible domain");
        }
        Ok(s)
    }

    /// Fill the positive resisting generalized force and gauge pressures.
    /// Extra fluid storage is NOT inferred. Return values are work-conjugate to
    /// `velocity`, which may be a discrete-gradient effort rather than endpoint p.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<FilmReport, FilmError> {
        if force.len() != self.ports || pressure.len() != self.cells.len() {
            return fail("invalid film output shape");
        }
        let s = self.solve(q, velocity)?;
        let mut candidate = [0.0; MAX_PORTS];
        self.project(&s.pressure, &mut candidate);
        let mut power = 0.0;
        for (k, edge) in self.channels.iter().enumerate() {
            let dp = s.pressure[edge.from] - edge.to.map_or(0.0, |j| s.pressure[j]);
            power += s.conductance[k] * dp * dp;
        }
        let work = dot(velocity, &candidate[..self.ports]);
        if candidate[..self.ports].iter().any(|x| !x.is_finite())
            || !power.is_finite()
            || !work.is_finite()
            || (work - power).abs() > 1e-9 * (work.abs() + power.abs()).max(1e-30)
        {
            return fail("film pressure solve failed reciprocal work balance");
        }
        force.copy_from_slice(&candidate[..self.ports]);
        pressure.copy_from_slice(&s.pressure[..self.cells.len()]);
        Ok(FilmReport {
            minimum_gap_m: s.minimum_gap,
            peak_pressure_pa: pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs())),
            dissipated_power_w: power,
        })
    }

    /// Exact directional action of the resisting port, with independently
    /// supplied configuration and effort directions. Differentiates hydraulic
    /// apertures too: `L dp = A B dv - (dL) p`. At a shut channel, h_+^3 has
    /// derivative zero. The same pressure factorization is reused for dp.
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        if dq.len() != self.ports
            || dv.len() != self.ports
            || out.len() != self.ports
            || dq.iter().chain(dv).any(|x| !x.is_finite())
        {
            return fail("invalid film tangent shape or direction");
        }
        let s = self.solve(q, velocity)?;
        let mut rhs = [0.0; MAX_CELLS];
        for (i, c) in self.cells.iter().enumerate() {
            rhs[i] = c.area_m2 * dot(&c.gap.closure, dv);
        }
        for (edge, factor) in self.channels.iter().zip(&self.factors) {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            let dh = -dot(&edge.gap.closure, dq);
            let dg = if h > 0.0 {
                3.0 * factor * h * h * dh
            } else {
                0.0
            };
            let change = dg * (s.pressure[edge.from] - edge.to.map_or(0.0, |j| s.pressure[j]));
            rhs[edge.from] -= change;
            if let Some(j) = edge.to {
                rhs[j] += change;
            }
        }
        solve_factored(&s.lower, &s.scale, self.cells.len(), &mut rhs)?;
        let mut candidate = [0.0; MAX_PORTS];
        self.project(&rhs, &mut candidate);
        if candidate[..self.ports].iter().any(|x| !x.is_finite()) {
            return fail("nonfinite film force tangent");
        }
        out.copy_from_slice(&candidate[..self.ports]);
        Ok(())
    }

    fn project(&self, p: &[f64; MAX_CELLS], f: &mut [f64; MAX_PORTS]) {
        for (i, cell) in self.cells.iter().enumerate() {
            let force = cell.area_m2 * p[i];
            for (j, weight) in cell.gap.closure.iter().enumerate() {
                f[j] += weight * force;
            }
        }
    }
}

// Local pressure-island elimination only, not a replacement mechanical solver.
// Diagonal equilibration preserves the operator. No pivot repair or gap floor.
struct Solution {
    lower: [f64; MAX_CELLS * MAX_CELLS],
    scale: [f64; MAX_CELLS],
    pressure: [f64; MAX_CELLS],
    conductance: [f64; MAX_CHANNELS],
    n: usize,
    minimum_gap: f64,
}
impl Solution {
    fn new(n: usize, minimum_gap: f64) -> Self {
        Self {
            lower: [0.0; MAX_CELLS * MAX_CELLS],
            scale: [0.0; MAX_CELLS],
            pressure: [0.0; MAX_CELLS],
            conductance: [0.0; MAX_CHANNELS],
            n,
            minimum_gap,
        }
    }
    fn factor(&mut self) -> Result<(), FilmError> {
        let n = self.n;
        for i in 0..n {
            let d = self.lower[i * n + i];
            if !d.is_finite() || d <= 0.0 {
                return fail("undrained or nonfinite film pressure cell");
            }
            self.scale[i] = d.sqrt();
        }
        for i in 0..n {
            for j in 0..n {
                self.lower[i * n + j] = (self.lower[i * n + j] / self.scale[i]) / self.scale[j];
            }
        }
        for i in 0..n {
            for j in 0..=i {
                let mut sum = self.lower[i * n + j];
                for k in 0..j {
                    sum -= self.lower[i * n + k] * self.lower[j * n + k];
                }
                if i == j {
                    if !sum.is_finite() || sum <= 64.0 * f64::EPSILON {
                        return fail("singular or unresolved trapped film pressure island");
                    }
                    self.lower[i * n + i] = sum.sqrt();
                } else {
                    self.lower[i * n + j] = sum / self.lower[j * n + j];
                }
            }
        }
        Ok(())
    }
}
fn solve_factored(
    l: &[f64; MAX_CELLS * MAX_CELLS],
    scale: &[f64; MAX_CELLS],
    n: usize,
    rhs: &mut [f64; MAX_CELLS],
) -> Result<(), FilmError> {
    for i in 0..n {
        let mut x = rhs[i] / scale[i];
        for j in 0..i {
            x -= l[i * n + j] * rhs[j];
        }
        rhs[i] = x / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut x = rhs[i];
        for j in i + 1..n {
            x -= l[j * n + i] * rhs[j];
        }
        rhs[i] = x / l[i * n + i];
    }
    for i in 0..n {
        rhs[i] /= scale[i];
    }
    if rhs[..n].iter().any(|x| !x.is_finite()) {
        return fail("nonfinite film pressure solution");
    }
    Ok(())
}
