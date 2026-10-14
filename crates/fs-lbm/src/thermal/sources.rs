//! Prescribed volumetric heating, external forces and independently measured
//! wall exchange on the existing thermal slabs. No second transport solver.
//!
//! A source is a signed temperature increment over ONE lattice step. Two
//! isotropic half-kicks surround the existing flow/heat step; the first kick
//! therefore participates in buoyancy immediately. Sources are frozen over the
//! step, not an implicit reaction law. Symmetric source splitting alone does
//! not establish second-order accuracy of the partitioned flow/heat solver.

pub mod transport;

use super::{ThermalForcingError, ThermalLbm};
use crate::core2::{Cell, CollisionModel2};
use crate::d3q19::ThermalLbm3;

/// Explicit limits for one checked, transactional source-coupled step.
#[derive(Debug, Clone, Copy)]
pub struct ThermalSourceLimits {
    /// Maximum total cells, including the two wall planes. The checked path
    /// temporarily holds a second complete solver state, plus transport scratch.
    pub max_cells: usize,
    /// Admitted temperature interval in the solver's temperature coordinates.
    /// Checked on the input, both source kicks and the transported candidate.
    pub temperature_range: [f64; 2],
    /// Absolute tolerance on the SUMMED temperature balance, in lattice units.
    pub balance_tolerance: f64,
}

/// Signed, summed thermal balance for one accepted step (fluid cells only).
/// Multiply each value by the same physical cell heat-capacity/temperature
/// scale to obtain energy. Positive wall exchange means heat INTO the fluid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalSourceBalance {
    /// Prescribed source increment, summed over fluid cells.
    pub source: f64,
    /// Change of stored temperature, summed cell-by-cell.
    pub stored_change: f64,
    /// Incoming-minus-outgoing populations on the bottom wall links.
    pub bottom_wall: f64,
    /// Incoming-minus-outgoing populations on the top wall links.
    pub top_wall: f64,
    /// `stored_change - source - bottom_wall - top_wall`.
    /// Independently measured wall links, not a residual relabelled as heat.
    pub defect: f64,
}

/// Refusal leaves the original physical state unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalSourceError {
    /// Invalid limits, layout, source or numerical state.
    Invalid(&'static str),
    /// Physical material/source conversion was refused.
    Material(crate::material::MaterialLbmError),
    /// Existing external-force admission rejected the candidate.
    Force(ThermalForcingError),
    /// A candidate temperature is outside the caller's interval.
    Temperature {
        /// Linear cell index.
        cell: usize,
        /// Candidate temperature in lattice coordinates.
        value: f64,
    },
    /// The independent thermal balance exceeds the requested tolerance.
    Balance {
        /// Observed balance residual.
        defect: f64,
    },
}
impl core::fmt::Display for ThermalSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "thermal source: {what}"),
            Self::Material(error) => write!(f, "thermal source: {error}"),
            Self::Force(error) => write!(f, "thermal source: {error}"),
            Self::Temperature { cell, value } => write!(f, "thermal source temperature {value} at cell {cell} is outside its admitted range"),
            Self::Balance { defect } => write!(f, "thermal source balance defect {defect} exceeds tolerance"),
        }
    }
}
impl std::error::Error for ThermalSourceError {}
impl From<ThermalForcingError> for ThermalSourceError {
    fn from(error: ThermalForcingError) -> Self { Self::Force(error) }
}

/// Physical limits for source-coupled material slabs.
#[derive(Debug, Clone, Copy)]
pub struct HeatSourceLimits {
    /// Maximum cells in either grid, including its two wall planes.
    /// This counts grid cells, not bytes; see [`ThermalSourceLimits`].
    pub max_cells: usize,
    /// Positive absolute-temperature interval [K] admitted by the caller's
    /// frozen-property material model and Boussinesq approximation.
    pub temperature_range_k: [f64; 2],
    /// Absolute energy-balance tolerance for one step [J].
    pub balance_tolerance_j: f64,
}

/// Physical heat balance for the represented 3-D volume, or for the explicit
/// extrusion depth of a 2-D slab. Positive wall energy enters the fluid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatBalance {
    /// Authored volumetric source integrated over cells and time [J].
    pub source_j: f64,
    /// Change in stored sensible thermal energy [J].
    pub stored_change_j: f64,
    /// Bottom-wall heat input, measured from population exchange [J].
    pub bottom_wall_j: f64,
    /// Top-wall heat input, measured from population exchange [J].
    pub top_wall_j: f64,
    /// Independently recomputed stored-minus-source-minus-wall defect [J].
    pub defect_j: f64,
}

impl From<crate::material::MaterialLbmError> for ThermalSourceError {
    fn from(error: crate::material::MaterialLbmError) -> Self { Self::Material(error) }
}

// Private adapters retain each solver's complete physical state. The source
// weights come from the existing D2Q5 constant or the D3Q7 constructor
// at unit wall temperature, rather than a duplicate quadrature table.
pub(crate) trait Slab<const D: usize, const Q: usize>: Sized {
    fn layout(&self) -> Result<(usize, usize), ThermalSourceError>;
    fn populations(&self) -> &[[f64; Q]];
    fn populations_mut(&mut self) -> &mut [[f64; Q]];
    fn snapshot(&self) -> (Self, [f64; Q]);
    fn walls(&self) -> (f64, f64);
    fn incoming() -> (usize, usize);
    fn flow_valid(&self) -> bool;
    fn advance(&mut self, forces: &[[f64; D]]) -> Result<(), ThermalForcingError>;
    fn momentum_relaxation(&self) -> &[f64];
    fn momentum_relaxation_mut(&mut self) -> &mut [f64];
    fn advance_local(&mut self, forces: &[[f64; D]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError>;
}

fn invalid(what: &'static str) -> ThermalSourceError { ThermalSourceError::Invalid(what) }

fn check_temperatures<const Q: usize>(pops: &[[f64; Q]], range: [f64; 2]) -> Result<(), ThermalSourceError> {
    for (cell, pop) in pops.iter().enumerate() {
        let value: f64 = pop.iter().sum();
        if !pop.iter().all(|p| p.is_finite()) || !value.is_finite() || value < range[0] || value > range[1] {
            return Err(ThermalSourceError::Temperature { cell, value });
        }
    }
    Ok(())
}

fn kick<const Q: usize>(pops: &mut [[f64; Q]], heat: &[f64], weights: &[f64; Q], plane: usize) {
    let end = pops.len() - plane;
    for (pop, &source) in pops[plane..end].iter_mut().zip(&heat[plane..end]) {
        // Preserve the exact old path for zero sources, including signed zero.
        if source != 0.0 {
            for (value, &weight) in pop.iter_mut().zip(weights) {
                *value += weight * (0.5 * source);
            }
        }
    }
}

pub(crate) fn source_step_mapped<S: Slab<D, Q>, const D: usize, const Q: usize, R>(
    state: &mut S, heat: &[f64], forces: &[[f64; D]], limits: ThermalSourceLimits,
    finish: impl FnOnce(ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<R, ThermalSourceError> {
    source_step_transport_mapped::<S, D, Q, R>(
        state, heat, forces, limits, None, |_, balance| finish(balance),
    )
}

pub(crate) fn source_step_transport_mapped<S: Slab<D, Q>, const D: usize, const Q: usize, R>(
    state: &mut S, heat: &[f64], forces: &[[f64; D]], limits: ThermalSourceLimits,
    transport: Option<transport::RelaxationField<'_>>,
    finish: impl FnOnce(&S, ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<R, ThermalSourceError> {
    let (cells, plane) = state.layout()?;
    if cells > limits.max_cells || !limits.balance_tolerance.is_finite() || limits.balance_tolerance < 0.0
        || !limits.temperature_range.iter().all(|x| x.is_finite())
        || limits.temperature_range[0] > limits.temperature_range[1]
    {
        return Err(invalid("invalid temperature/balance limits or cell budget exceeded"));
    }
    if heat.len() != cells || forces.len() != cells { return Err(invalid("one source and force vector per cell, including walls, is required")); }
    if !heat.iter().all(|x| x.is_finite()) || !forces.iter().flatten().all(|x| x.is_finite()) {
        return Err(invalid("source and force fields must be finite"));
    }
    if heat[..plane].iter().chain(&heat[cells-plane..]).any(|x| *x != 0.0) {
        return Err(invalid("wall source entries must be zero; wall temperatures are prescribed"));
    }
    if !state.flow_valid() { return Err(invalid("invalid thermal/flow state or unsupported slab topology")); }
    check_temperatures(state.populations(), limits.temperature_range)?;
    let (bottom, top) = state.walls();
    if ![bottom, top].iter().all(|t| t.is_finite() && *t >= limits.temperature_range[0] && *t <= limits.temperature_range[1]) {
        return Err(invalid("wall temperatures are outside the admitted interval"));
    }
    let source: f64 = heat[plane..cells-plane].iter().sum();
    if !source.is_finite() { return Err(invalid("summed source is not representable")); }
    if let Some(field) = transport { field.admit(cells)?; }
    let (mut candidate, weights) = state.snapshot();
    if let Some(field) = transport {
        candidate.momentum_relaxation_mut().copy_from_slice(field.momentum);
    }
    kick(candidate.populations_mut(), heat, &weights, plane);
    check_temperatures(candidate.populations(), limits.temperature_range)?;
    match transport {
        Some(field) => candidate.advance_local(forces, field.thermal)?,
        None => candidate.advance(forces)?,
    }
    if !candidate.flow_valid() { return Err(invalid("transport produced an inadmissible flow state")); }
    check_temperatures(candidate.populations(), limits.temperature_range)?;
    // Anti-bounce-back gives incoming = 2*w*Twall - outgoing. Therefore the
    // independently measured link exchange is 2*(incoming - w*Twall).
    // Measure BEFORE the second source kick, which is not wall transport.
    let (qb, qt) = S::incoming();
    let bottom_wall: f64 = candidate.populations()[plane..2*plane].iter()
        .map(|p| 2.0 * (p[qb] - weights[qb] * bottom)).sum();
    let top_wall: f64 = candidate.populations()[cells-2*plane..cells-plane].iter()
        .map(|p| 2.0 * (p[qt] - weights[qt] * top)).sum();
    kick(candidate.populations_mut(), heat, &weights, plane);
    check_temperatures(candidate.populations(), limits.temperature_range)?;
    let stored_change: f64 = candidate.populations()[plane..cells-plane].iter()
        .zip(&state.populations()[plane..cells-plane])
        .map(|(after, before)| after.iter().zip(before).map(|(a,b)| a-b).sum::<f64>()).sum();
    let defect = stored_change - source - bottom_wall - top_wall;
    if !defect.is_finite() || defect.abs() > limits.balance_tolerance {
        return Err(ThermalSourceError::Balance { defect });
    }
    // Coefficients are borrowed for ONE step, like sources. Preserve background
    // rates so a later legacy step cannot combine stale local viscosity with
    // default thermal diffusivity. Populations and actual forces are retained.
    if transport.is_some() {
        candidate.momentum_relaxation_mut().copy_from_slice(state.momentum_relaxation());
    }
    let result = finish(&candidate, ThermalSourceBalance { source, stored_change, bottom_wall, top_wall, defect })?;
    *state = candidate;
    Ok(result)
}

fn flow_valid<const Q: usize, const D: usize>(f: &[[f64; Q]], tau: &[f64], fext: &[[f64; D]], tau_g: f64) -> bool {
    tau_g.is_finite() && tau_g > 0.5
        && tau.iter().all(|t| t.is_finite() && *t > 0.5)
        && fext.iter().flatten().all(|x| x.is_finite())
        && f.iter().all(|p| {
            let rho: f64 = p.iter().sum();
            rho.is_finite() && rho > 0.0 && p.iter().all(|x| x.is_finite())
        })
}

impl Slab<2, 5> for ThermalLbm {
    fn layout(&self) -> Result<(usize, usize), ThermalSourceError> {
        let g = &self.grid;
        let n = g.nx.checked_mul(g.ny).ok_or(invalid("cell count overflows"))?;
        if g.nx == 0 || g.ny < 4 || [g.f.len(), g.flags.len(), g.tau.len(), g.fext.len(), self.gpop.len()].iter().any(|l| *l != n) {
            return Err(invalid("invalid 2D slab storage"));
        }
        Ok((n, g.nx))
    }
    fn populations(&self) -> &[[f64; 5]] { &self.gpop }
    fn populations_mut(&mut self) -> &mut [[f64; 5]] { &mut self.gpop }
    fn snapshot(&self) -> (Self, [f64; 5]) {
        (Self {
            grid: self.grid.clone(), gpop: self.gpop.clone(), tau_g: self.tau_g,
            t_bottom: self.t_bottom, t_top: self.t_top, gbeta: self.gbeta,
            t_ref: self.t_ref, scratch: Vec::new(),
        }, super::W5)
    }
    fn walls(&self) -> (f64, f64) { (self.t_bottom, self.t_top) }
    fn incoming() -> (usize, usize) { (2, 4) }
    fn flow_valid(&self) -> bool {
        let g = &self.grid;
        flow_valid(&g.f, &g.tau, &g.fext, self.tau_g) && self.gbeta.is_finite() && self.t_ref.is_finite()
            && g.periodic_x && !g.periodic_y && matches!(g.collision, CollisionModel2::Bgk)
            && g.g.iter().all(|x| x.is_finite())
            && g.flags.iter().enumerate().all(|(i, flag)| *flag == if i < g.nx || i >= g.f.len()-g.nx { Cell::Wall } else { Cell::Fluid })
    }
    fn advance(&mut self, forces: &[[f64; 2]]) -> Result<(), ThermalForcingError> { self.step_with_external_forces(forces) }
    fn momentum_relaxation(&self) -> &[f64] { &self.grid.tau }
    fn momentum_relaxation_mut(&mut self) -> &mut [f64] { &mut self.grid.tau }
    fn advance_local(&mut self, forces: &[[f64; 2]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError> {
        self.step_with_local_thermal_relaxation(forces, thermal_tau)
    }
}
impl Slab<3, 7> for ThermalLbm3 {
    fn layout(&self) -> Result<(usize, usize), ThermalSourceError> {
        let g = &self.grid;
        let plane = g.nx.checked_mul(g.ny).ok_or(invalid("wall plane overflows"))?;
        let n = plane.checked_mul(g.nz).ok_or(invalid("cell count overflows"))?;
        if plane == 0 || g.nz < 4 || [g.f.len(), g.tau.len(), g.fext.len(), self.gpop.len()].iter().any(|l| *l != n) {
            return Err(invalid("invalid 3D slab storage"));
        }
        Ok((n, plane))
    }
    fn populations(&self) -> &[[f64; 7]] { &self.gpop }
    fn populations_mut(&mut self) -> &mut [[f64; 7]] { &mut self.gpop }
    fn snapshot(&self) -> (Self, [f64; 7]) {
        let g = &self.grid;
        let mut s = Self::slab(g.nx, g.ny, g.nz-2, 1.0, self.tau_g, self.gbeta);
        let weights = s.gpop[0];
        s.grid.f.clone_from(&g.f); s.grid.tau.clone_from(&g.tau); s.grid.fext.clone_from(&g.fext);
        s.gpop.clone_from(&self.gpop);
        s.t_bottom = self.t_bottom; s.t_top = self.t_top; s.t_ref = self.t_ref;
        (s, weights)
    }
    fn walls(&self) -> (f64, f64) { (self.t_bottom, self.t_top) }
    fn incoming() -> (usize, usize) { (5, 6) }
    fn flow_valid(&self) -> bool {
        flow_valid(&self.grid.f, &self.grid.tau, &self.grid.fext, self.tau_g) && self.gbeta.is_finite() && self.t_ref.is_finite()
    }
    fn advance(&mut self, forces: &[[f64; 3]]) -> Result<(), ThermalForcingError> { self.step_with_external_forces(forces) }
    fn momentum_relaxation(&self) -> &[f64] { &self.grid.tau }
    fn momentum_relaxation_mut(&mut self) -> &mut [f64] { &mut self.grid.tau }
    fn advance_local(&mut self, forces: &[[f64; 3]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError> {
        self.step_with_local_thermal_relaxation(forces, thermal_tau)
    }
}

impl ThermalLbm {
    /// Apply prescribed signed temperature increments and external forces for
    /// one step. Includes walls in each input array; wall heating must be zero.
    /// Both temperature half-kicks, flow, heat transport and the measured wall
    /// balance must pass before the candidate replaces the original state.
    /// This checked path supports the constructor's periodic BGK slab, not
    /// edited obstacle masks or alternative collision operators.
    pub fn step_with_sources(&mut self, heat_per_step: &[f64], forces: &[[f64; 2]], limits: ThermalSourceLimits) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_mapped::<Self, 2, 5, _>(self, heat_per_step, forces, limits, Ok)
    }
}
impl ThermalLbm3 {
    /// Three-dimensional counterpart of [`ThermalLbm::step_with_sources`].
    /// Uses the original D3Q7 equilibrium weights and the same source and
    /// balance implementation. Source rates are prescribed, not reaction laws.
    pub fn step_with_sources(&mut self, heat_per_step: &[f64], forces: &[[f64; 3]], limits: ThermalSourceLimits) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_mapped::<Self, 3, 7, _>(self, heat_per_step, forces, limits, Ok)
    }
}
