//! Quasistatic, incompressible Reynolds-film pressure on a finite-volume graph.
//!
//! This is an orthogonal image to the compressible 1-D gas-film time solver.
//! It has no clock or pressure history. At every caller-supplied configuration,
//! `L(h) psi = A B v`, and the positive resisting force is `B^T A p`.
//! With constant viscosity, `psi=p`. The optional Barus pressure-viscosity law
//! uses `psi=(1-exp(-alpha*p))/alpha`, without another pressure iteration.
//! Mechanical power is `sum_edges G (psi_i-psi_j) (p_i-p_j) >= 0`.
//! A caller must evaluate this port *inside* its implicit mechanical equation.
//!
//! Poiseuille conductance is `width * h_+^3 / (12 * viscosity * length)`.
//! Each channel may have its own positive viscosity, held fixed during a solve.
//! A channel coefficient is uniform across its gap; for serial segments with
//! the same aperture it is the LENGTH-weighted arithmetic viscosity, since
//! viscous resistance, not fluidity, adds along the passage.
//! A nonpositive CHANNEL aperture closes that passage exactly; a collapsed
//! VOLUME cell, disconnected trapped pocket, or exceeded validity limit refuses.
//! There is no gap floor, pressure clipping, diagonal jitter, or hidden leakage.
//! Reservoir pressures are zero gauge. Suction is retained within the supplied
//! pressure limit. Compressibility, gas inertia, thermal effects, slip, moving
//! topology and fluid acoustic radiation are not part of this image.
/// Mass-conserving liquid-film rupture, starvation and reformation on this graph.
pub mod cavitation;
/// Compressible mass storage and transport on the same geometric pressure graph.
pub mod isothermal;
/// Exact pressure-flow reduction for prescribed through-gap viscosity layers.
pub mod stratified;
/// Work-conjugate sliding walls, Couette entrainment and shear on this same graph.
pub mod shear;

/// Regularized power-law pressure flow on the original Reynolds graph.
pub mod power_law;

use core::fmt;
use std::sync::Arc;

/// Explicit dense pressure-island ceiling, not a real-time performance claim.
pub const MAX_CELLS: usize = 64;
/// Maximum number of generalized mechanical ports.
pub const MAX_PORTS: usize = 256;
/// Maximum number of hydraulic passages in one pressure island.
pub const MAX_CHANNELS: usize = 256;

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
    /// Orientation of positive volume flow out of this cell.
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

/// Caller-declared domain of this incompressible, thin-gap viscous image.
#[derive(Debug, Clone, Copy)]
pub struct FilmLimits {
    /// Positive minimum retained VOLUME-cell gap [m]. Never a replacement gap.
    pub minimum_cell_gap_m: f64,
    /// Maximum admitted cell/channel gap [m].
    pub maximum_gap_m: f64,
    /// Maximum absolute GAUGE pressure [Pa]. The caller independently bounds
    /// liquid incompressibility and any pressure-viscosity fit over this interval.
    /// A gas in this pressure-eliminated image additionally needs excursions
    /// small relative to its absolute reservoir pressure.
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

/// One passage's flow and irreversible work from the same pressure solve.
/// This has no independent time, mass inventory or heat-partition assumption.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FilmChannelReport {
    /// Signed volume flow from `from` to `to` (or ambient) [m^3/s].
    pub volume_flow_m3_s: f64,
    /// Nonnegative `Q * (p_from - p_to)` [W], including the Barus law.
    /// A closed passage contributes exactly zero flow and power.
    pub dissipated_power_w: f64,
}

/// Immutable film geometry and transport law. Evaluation and tangent use bounded stack
/// scratch, never mutate pressure history, and allocate only at construction.
/// Clones and viscosity recompilations share the admitted geometry; only transport
/// coefficients are copied or rebuilt, independent of closure-port storage size.
#[derive(Debug, Clone)]
pub struct ResistiveFilm {
    cells: Arc<[FilmCell]>,
    channels: Arc<[FilmChannel]>,
    factors: Vec<f64>,
    channel_viscosities_pa_s: Vec<f64>,
    pressure_viscosity_coefficient_per_pa: f64,
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
            channel_viscosities_pa_s: vec![viscosity_pa_s; channels.len()],
            cells: cells.into(),
            channels: channels.into(),
            factors,
            pressure_viscosity_coefficient_per_pa: 0.0,
            ports,
            limits,
        })
    }

    /// Common dynamic viscosity at zero gauge pressure [Pa s], only when every
    /// channel has exactly the same value. A heterogeneous film has no single
    /// viscosity: this returns `None` rather than inventing an effective mean.
    #[must_use]
    pub fn viscosity_pa_s(&self) -> Option<f64> {
        let first = self.channel_viscosities_pa_s[0];
        self.channel_viscosities_pa_s
            .iter()
            .all(|&mu| mu == first)
            .then_some(first)
    }

    /// Exact zero-gauge dynamic viscosities [Pa s] in original channel order.
    /// Geometry, gap changes and pressure work use these same coefficients.
    #[must_use]
    pub fn channel_viscosities_pa_s(&self) -> &[f64] {
        &self.channel_viscosities_pa_s
    }

    /// Exponential pressure-viscosity coefficient [1/Pa], zero for a constant law.
    #[must_use]
    pub const fn pressure_viscosity_coefficient_per_pa(&self) -> f64 {
        self.pressure_viscosity_coefficient_per_pa
    }

    /// Retain this graph with the Barus law `mu(p)=mu_0 exp(alpha*p)`.
    ///
    /// `p` is GAUGE pressure relative to every zero-gauge reservoir; the existing
    /// channel viscosity is `mu_0` there. The caller supplies a nonnegative coefficient
    /// and a pressure range supporting this constitutive approximation. This
    /// does not infer a material source, elastic deformation, liquid density,
    /// cavitation or thermal evolution. It remains the classical isothermal
    /// incompressible Reynolds reduction, not a general high-pressure fluid law.
    ///
    /// Integrating fluidity exactly along each passage gives
    /// `Q=G_0 (psi_i-psi_j)`, where `psi=(1-exp(-alpha*p))/alpha`.
    /// The original Cholesky solve therefore remains valid. If `alpha*psi>=1`,
    /// the requested flow has no finite pressure in this model and refuses.
    /// Zero alpha restores the original constant-viscosity arithmetic exactly.
    /// Gas-storage and cavitating owners require their own derivation and refuse
    /// this pressure-dependent image.
    pub fn with_barus_viscosity(
        &self,
        pressure_viscosity_coefficient_per_pa: f64,
    ) -> Result<Self, FilmError> {
        if !pressure_viscosity_coefficient_per_pa.is_finite()
            || pressure_viscosity_coefficient_per_pa < 0.0
            || !(pressure_viscosity_coefficient_per_pa * self.limits.maximum_pressure_pa)
                .is_finite()
        {
            return fail("invalid Barus pressure-viscosity coefficient or pressure interval");
        }
        let mut result = self.clone();
        result.pressure_viscosity_coefficient_per_pa = pressure_viscosity_coefficient_per_pa;
        Ok(result)
    }

    /// Recompile a UNIFORM reference viscosity on this exact graph and pressure law.
    ///
    /// Cells, apertures, force ports, pressure-viscosity coefficient and limits
    /// remain unchanged. This stateless resistance has no stored energy to remap.
    /// The caller must hold
    /// the returned coefficient fixed during its mechanical solve and account
    /// for the resulting viscous work with that same operator. This does not
    /// update a compressible film's gas storage or infer thermal expansion.
    pub fn with_viscosity(&self, viscosity_pa_s: f64) -> Result<Self, FilmError> {
        let viscosities = [viscosity_pa_s; MAX_CHANNELS];
        self.with_channel_viscosities(&viscosities[..self.channels.len()])
    }

    /// Compile one positive dynamic viscosity [Pa s] per ORIGINAL channel.
    ///
    /// The coefficients are constant across each channel's gap and length at
    /// zero gauge pressure; an existing common Barus coefficient is retained.
    /// For a fixed aperture with viscosity varying along the flow direction,
    /// supply the length-weighted arithmetic mean, not a harmonic mean. Gapwise
    /// thermal profiles require their own momentum integration and are not
    /// represented by substituting a cell or wall viscosity.
    /// [`stratified::StratifiedFilm`] supplies that integration for prescribed
    /// piecewise-constant layers while retaining their original viscosities.
    ///
    /// Original geometry, limits and force ports remain exact and share their
    /// admitted storage without another topology/closure validation or copy.
    /// All conductance factors are compiled once here; pressure, tangent and work use the same
    /// existing solve with no added evaluation allocation or iteration. The
    /// caller owns source/temperature admission and holds these values fixed
    /// inside each mechanical solve. This supplies no heat or mass evolution.
    /// Refuses a wrong count, nonpositive/nonfinite viscosity or unrepresentable
    /// conductance without changing the original operator.
    pub fn with_channel_viscosities(&self, viscosities_pa_s: &[f64]) -> Result<Self, FilmError> {
        if viscosities_pa_s.len() != self.channels.len() {
            return fail("film viscosity count must match the original channels");
        }
        let mut factors = Vec::with_capacity(self.channels.len());
        for (edge, &viscosity) in self.channels.iter().zip(viscosities_pa_s) {
            if !viscosity.is_finite() || viscosity <= 0.0 {
                return fail("film channel viscosity must be finite and positive");
            }
            let factor = (edge.width_m / edge.length_m) / (12.0 * viscosity);
            if !factor.is_finite() || factor <= 0.0 {
                return fail("film conductance scale overflow");
            }
            factors.push(factor);
        }
        Ok(Self {
            cells: Arc::clone(&self.cells),
            channels: Arc::clone(&self.channels),
            factors,
            channel_viscosities_pa_s: viscosities_pa_s.to_vec(),
            pressure_viscosity_coefficient_per_pa: self.pressure_viscosity_coefficient_per_pa,
            ports: self.ports,
            limits: self.limits,
        })
    }

    /// Number of work-conjugate mechanical coordinates.
    #[must_use]
    pub fn port_count(&self) -> usize {
        self.ports
    }
    /// Number of retained pressure cells.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    /// Number of retained hydraulic passages, including currently shut channels.
    #[must_use]
    pub fn channel_count(&self) -> usize {
        self.channels.len()
    }

    /// Check cell and channel gap validity without advancing or solving pressure.
    pub fn validate_configuration(&self, q: &[f64]) -> Result<f64, FilmError> {
        if q.len() != self.ports || q.iter().any(|x| !x.is_finite()) {
            return fail("invalid film configuration");
        }
        let mut minimum = f64::INFINITY;
        for c in self.cells.iter() {
            let h = c.gap.reference_m - dot(&c.gap.closure, q);
            if !h.is_finite() || h < self.limits.minimum_cell_gap_m || h > self.limits.maximum_gap_m
            {
                return fail("film volume cell outside declared gap domain");
            }
            minimum = minimum.min(h);
        }
        for edge in self.channels.iter() {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            if !h.is_finite() || h > self.limits.maximum_gap_m {
                return fail("film channel outside declared gap domain");
            }
        }
        Ok(minimum)
    }

    fn solve(&self, q: &[f64], v: &[f64]) -> Result<Solution, FilmError> {
        self.solve_with_entrainment::<false>(q, v, &[])
    }

    // The sliding-wall owner supplies oriented volume entrainment and its
    // conjugate wall tractions. Matrix assembly and elimination remain here.
    // The ordinary pressure image specializes away the additional RHS path.
    fn solve_with_entrainment<const ENTRAINMENT: bool>(
        &self,
        q: &[f64],
        v: &[f64],
        entrainment: &[f64],
    ) -> Result<Solution, FilmError> {
        if ENTRAINMENT
            && (entrainment.len() != self.channels.len()
                || entrainment.iter().any(|flow| !flow.is_finite()))
        {
            return fail("invalid oriented film entrainment");
        }
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
        if ENTRAINMENT {
            for (edge, &flow) in self.channels.iter().zip(entrainment) {
                if flow != 0.0 {
                    s.pressure[edge.from] -= flow;
                    if let Some(j) = edge.to {
                        s.pressure[j] += flow;
                    }
                }
            }
        }
        s.factor()?;
        solve_factored(&s.lower, &s.scale, n, &mut s.pressure)?;
        s.potential = s.pressure;
        let alpha = self.pressure_viscosity_coefficient_per_pa;
        if alpha != 0.0 {
            for (p, psi) in s.pressure[..n].iter_mut().zip(&s.potential[..n]) {
                let scaled = alpha * psi;
                if !scaled.is_finite() || scaled >= 1.0 {
                    return fail("Barus film flow exceeds its finite-pressure limit");
                }
                // ln_1p retains the constant-viscosity limit. An underflowed
                // alpha*psi has no representable pressure-viscosity correction.
                *p = if scaled == 0.0 {
                    *psi
                } else {
                    -(-scaled).ln_1p() / alpha
                };
                if !p.is_finite() {
                    return fail("Barus film pressure is unrepresentable");
                }
            }
            // Each passage retains its own reference viscosity. Pressure is
            // monotone along a passage, so checking both ends bounds its full
            // Barus viscosity; do not test an unrelated cell's coefficient.
            for (edge, &mu) in self.channels.iter().zip(&self.channel_viscosities_pa_s) {
                for psi in [
                    s.potential[edge.from],
                    edge.to.map_or(0.0, |j| s.potential[j]),
                ] {
                    let viscosity = mu / (1.0 - alpha * psi);
                    if !viscosity.is_finite() || viscosity <= 0.0 {
                        return fail("Barus film viscosity is unrepresentable");
                    }
                }
            }
        }
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
        self.evaluate::<false>(q, velocity, force, pressure, &mut [])
    }

    /// Evaluate force, pressure and every passage's flow/work in one solve.
    ///
    /// Reports follow the original channel order. Their powers sum to the
    /// same resisting-port power; a thermal owner can route selected passages
    /// once to explicit finite stores. For pressure-dependent viscosity the
    /// power uses physical pressure drop, not the squared transformed potential.
    /// No wall/fluid partition, fluid storage, or advection model is inferred.
    /// All output slices remain untouched if any admission or work check fails.
    pub fn evaluate_with_channels_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
    ) -> Result<FilmReport, FilmError> {
        self.evaluate::<true>(q, velocity, force, pressure, channels)
    }

    // Specialization removes channel-report scratch/writes from ordinary force
    // and Newton evaluations. Both images use precisely the same arithmetic.
    fn evaluate<const CHANNEL_REPORTS: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
    ) -> Result<FilmReport, FilmError> {
        if force.len() != self.ports || pressure.len() != self.cells.len()
            || (CHANNEL_REPORTS && channels.len() != self.channels.len()) {
            return fail("invalid film output shape");
        }
        let s = self.solve(q, velocity)?;
        let mut candidate = [0.0; MAX_PORTS];
        let mut channel_candidate = [FilmChannelReport::default(); MAX_CHANNELS];
        self.project(&s.pressure, &mut candidate);
        let mut power = 0.0;
        for (k, edge) in self.channels.iter().enumerate() {
            let dp = s.pressure[edge.from] - edge.to.map_or(0.0, |j| s.pressure[j]);
            let dpsi = s.potential[edge.from] - edge.to.map_or(0.0, |j| s.potential[j]);
            let flow = s.conductance[k] * dpsi;
            let loss = flow * dp;
            if CHANNEL_REPORTS {
                if !flow.is_finite() || !loss.is_finite() || loss < 0.0 {
                    return fail("nonfinite or negative film channel work");
                }
                channel_candidate[k] = FilmChannelReport {
                    volume_flow_m3_s: flow,
                    dissipated_power_w: loss,
                };
            }
            power += loss;
        }
        let work = dot(velocity, &candidate[..self.ports]);
        if candidate[..self.ports].iter().any(|x| !x.is_finite())
            || !power.is_finite()
            || power < 0.0
            || !work.is_finite()
            || (work - power).abs() > 1e-9 * (work.abs() + power.abs()).max(1e-30)
        {
            return fail("film pressure solve failed reciprocal work balance");
        }
        force.copy_from_slice(&candidate[..self.ports]);
        pressure.copy_from_slice(&s.pressure[..self.cells.len()]);
        if CHANNEL_REPORTS {
            channels.copy_from_slice(&channel_candidate[..self.channels.len()]);
        }
        Ok(FilmReport {
            minimum_gap_m: s.minimum_gap,
            peak_pressure_pa: pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs())),
            dissipated_power_w: power,
        })
    }

    /// Exact directional action of the resisting port, with independently
    /// supplied configuration and effort directions. Differentiates hydraulic
    /// apertures too: `L dpsi = A B dv - (dL) psi`, followed by
    /// `dp_i=dpsi_i/(1-alpha*psi_i)`. At a shut channel, h_+^3 has derivative
    /// zero. The same pressure factorization is reused. For pressure-dependent
    /// viscosity the complete mechanical tangent is generally NONSYMMETRIC;
    /// it must not be replaced by a symmetric constant damping matrix.
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
        self.tangent_from_solution(q, &s, dq, dv, out)
    }

    // Reuse a solution admitted by this exact operator. The caller has
    // already checked direction/output shapes. Constitutive continuations
    // retain the original tangent arithmetic without a second pressure solve.
    fn tangent_from_solution(
        &self,
        q: &[f64],
        s: &Solution,
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
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
            let change = dg * (s.potential[edge.from] - edge.to.map_or(0.0, |j| s.potential[j]));
            rhs[edge.from] -= change;
            if let Some(j) = edge.to {
                rhs[j] += change;
            }
        }
        solve_factored(&s.lower, &s.scale, self.cells.len(), &mut rhs)?;
        let alpha = self.pressure_viscosity_coefficient_per_pa;
        if alpha != 0.0 {
            for (dp, psi) in rhs[..self.cells.len()].iter_mut().zip(&s.potential) {
                *dp /= 1.0 - alpha * psi;
            }
        }
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
    potential: [f64; MAX_CELLS],
    conductance: [f64; MAX_CHANNELS],
    n: usize,
    minimum_gap: f64,
}
impl Solution {
    #[allow(clippy::large_stack_arrays)] // Fixed-size stack scratch for bounded MAX_CELLS elimination.
    fn new(n: usize, minimum_gap: f64) -> Self {
        Self {
            lower: [0.0; MAX_CELLS * MAX_CELLS],
            scale: [0.0; MAX_CELLS],
            pressure: [0.0; MAX_CELLS],
            potential: [0.0; MAX_CELLS],
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

#[cfg(test)]
mod viscosity_tests;
