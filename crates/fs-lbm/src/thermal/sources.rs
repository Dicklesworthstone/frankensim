//! Prescribed volumetric heating, external forces and independently measured
//! wall exchange on the existing thermal slabs. No second transport solver.
//!
//! A source is a signed temperature increment over ONE lattice step. Two
//! equilibrium half-kicks surround the existing flow/heat step; the first
//! kick therefore participates in buoyancy immediately. Each kick preserves
//! thermal nonequilibrium and rebinds the advective-flux history at the current
//! velocity. Sources are frozen over the step, not an implicit reaction law.
//! Symmetric source splitting alone does not establish second-order accuracy
//! of the partitioned flow/heat solver.

pub mod sgs3;
pub mod transport;

use super::{ThermalForcingError, ThermalLbm};
use super::boundary::{ThermalWall, ThermalWallError, ThermalWallExchange, ThermalWallField};
use crate::core2::{Cell, CollisionModel2};
use crate::d3q19::{CollisionModel3, ThermalLbm3};

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
    /// The selected thermal boundary declaration is invalid.
    Wall(ThermalWallError),
    /// A candidate temperature is outside the caller's interval.
    Temperature {
        /// Linear cell index.
        cell: usize,
        /// Candidate temperature in lattice coordinates.
        value: f64,
    },
    /// A resolved halfway wall trace leaves the admitted scalar interval.
    WallTemperature {
        /// Bottom = 0, top = 1.
        wall: usize,
        /// Smallest individual trace on this wall.
        minimum: f64,
        /// Largest individual trace on this wall.
        maximum: f64,
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
            Self::Wall(error) => write!(f, "thermal source: {error}"),
            Self::Temperature { cell, value } => write!(f, "thermal source temperature {value} at cell {cell} is outside its admitted range"),
            Self::WallTemperature { wall, minimum, maximum } => write!(f,
                "thermal source wall {wall} trace [{minimum}, {maximum}] is outside its admitted range"),
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
    fn source_half_kick(&mut self, heat: &[f64], plane: usize);
    fn snapshot(&self) -> (Self, [f64; Q]);
    // Original coordinate anchors, retained for physical-unit consistency checks.
    fn walls(&self) -> (f64, f64);
    fn thermal_wall_fields(&self) -> [ThermalWallField<'_>; 2];
    fn wall_exchange(&self) -> Option<[ThermalWallExchange; 2]>;
    fn flow_valid(&self) -> bool;
    fn advance(&mut self, forces: &[[f64; D]]) -> Result<(), ThermalForcingError>;
    fn momentum_relaxation(&self) -> &[f64];
    fn momentum_relaxation_mut(&mut self) -> &mut [f64];
    fn momentum_relaxation_valid(&self, tau: &[f64]) -> bool;
    fn advance_local(&mut self, forces: &[[f64; D]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError>;
    fn advance_local_resolved<E: From<ThermalForcingError>>(
        &mut self, forces: &[[f64; D]], thermal_tau: &[f64],
        resolve: impl FnOnce(&Self) -> Result<Vec<[f64; D]>, E>,
    ) -> Result<(), E>;
    fn advance_tensor(&mut self, forces: &[[f64; D]], diffusion: super::tensor::DiffusivityField<'_, D>) -> Result<(), ThermalForcingError>;
}

fn invalid(what: &'static str) -> ThermalSourceError { ThermalSourceError::Invalid(what) }

fn check_temperatures<const Q: usize>(
    pops: &[[f64; Q]], plane: usize, range: [f64; 2],
) -> Result<(), ThermalSourceError> {
    // Ghost wall planes carry no thermal storage. Flux/film surfaces have
    // resolved traces, not the constructor's unused ghost temperatures.
    for (offset, pop) in pops[plane..pops.len()-plane].iter().enumerate() {
        let value: f64 = pop.iter().sum();
        if !pop.iter().all(|p| p.is_finite()) || !value.is_finite() || value < range[0] || value > range[1] {
            return Err(ThermalSourceError::Temperature { cell: offset + plane, value });
        }
    }
    Ok(())
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
    source_step_impl::<S, D, Q, R>(state, heat, forces, limits, transport, None, finish)
}

pub(crate) fn source_step_tensor_mapped<S: Slab<D, Q>, const D: usize, const Q: usize, R>(
    state: &mut S, heat: &[f64], forces: &[[f64; D]], limits: ThermalSourceLimits,
    diffusion: super::tensor::DiffusivityField<'_, D>,
    finish: impl FnOnce(&S, ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<R, ThermalSourceError> {
    source_step_impl::<S, D, Q, R>(state, heat, forces, limits, None, Some(diffusion), finish)
}

fn source_step_impl<S: Slab<D, Q>, const D: usize, const Q: usize, R>(
    state: &mut S, heat: &[f64], forces: &[[f64; D]], limits: ThermalSourceLimits,
    transport: Option<transport::RelaxationField<'_>>,
    diffusion: Option<super::tensor::DiffusivityField<'_, D>>,
    finish: impl FnOnce(&S, ThermalSourceBalance) -> Result<R, ThermalSourceError>,
) -> Result<R, ThermalSourceError> {
    source_step_with_advance::<S, D, Q, R, ThermalSourceError>(
        state, heat, forces, limits, transport, diffusion,
        |candidate| {
            if let Some(field) = diffusion {
                candidate.advance_tensor(forces, field)?;
            } else {
                match transport {
                    Some(field) => candidate.advance_local(forces, field.thermal)?,
                    None => candidate.advance(forces)?,
                }
            }
            Ok(())
        },
        finish,
    )
}

// The numerical owner may resolve a force after its one thermal transport.
// Keep all source, wall, storage, endpoint and publication gates here, and
// propagate constitutive errors without erasing their original source type.
#[allow(clippy::too_many_arguments)]
fn source_step_with_advance<
    S: Slab<D, Q>, const D: usize, const Q: usize, R, E: From<ThermalSourceError>,
>(
    state: &mut S, heat: &[f64], forces: &[[f64; D]], limits: ThermalSourceLimits,
    transport: Option<transport::RelaxationField<'_>>,
    diffusion: Option<super::tensor::DiffusivityField<'_, D>>,
    advance: impl FnOnce(&mut S) -> Result<(), E>,
    finish: impl FnOnce(&S, ThermalSourceBalance) -> Result<R, E>,
) -> Result<R, E> {
    let (cells, plane) = state.layout()?;
    if cells > limits.max_cells || !limits.balance_tolerance.is_finite() || limits.balance_tolerance < 0.0
        || !limits.temperature_range.iter().all(|x| x.is_finite())
        || limits.temperature_range[0] > limits.temperature_range[1]
    {
        return Err(invalid("invalid temperature/balance limits or cell budget exceeded").into());
    }
    if heat.len() != cells || forces.len() != cells { return Err(invalid("one source and force vector per cell, including walls, is required").into()); }
    if !heat.iter().all(|x| x.is_finite()) || !forces.iter().flatten().all(|x| x.is_finite()) {
        return Err(invalid("source and force fields must be finite").into());
    }
    if heat[..plane].iter().chain(&heat[cells-plane..]).any(|x| *x != 0.0) {
        return Err(invalid("ghost wall source entries must be zero; use the thermal wall law for surface heat").into());
    }
    if !state.flow_valid() { return Err(invalid("invalid thermal/flow state or unsupported slab topology").into()); }
    check_temperatures(state.populations(), plane, limits.temperature_range)?;
    for field in state.thermal_wall_fields() {
        field.validate(plane).map_err(ThermalSourceError::Wall)?;
        for wall in field.laws() {
            if let ThermalWall::Temperature(temperature) = wall
                && (temperature < limits.temperature_range[0] || temperature > limits.temperature_range[1])
            {
                return Err(invalid("prescribed wall temperature is outside the admitted interval").into());
            }
        }
    }
    let source: f64 = heat[plane..cells-plane].iter().sum();
    if !source.is_finite() { return Err(invalid("summed source is not representable").into()); }
    if let Some(field) = transport {
        field.admit(cells)?;
        if !state.momentum_relaxation_valid(field.momentum) {
            return Err(invalid("local transport is incompatible with the selected collision rates").into());
        }
    }
    if let Some(field) = diffusion { field.admit(cells)?; }
    let (mut candidate, _) = state.snapshot();
    if let Some(field) = transport {
        candidate.momentum_relaxation_mut().copy_from_slice(field.momentum);
    }
    candidate.source_half_kick(heat, plane);
    check_temperatures(candidate.populations(), plane, limits.temperature_range)?;
    advance(&mut candidate)?;
    if !candidate.flow_valid() { return Err(invalid("transport produced an inadmissible flow state").into()); }
    check_temperatures(candidate.populations(), plane, limits.temperature_range)?;
    // The streaming owner measures incoming minus outgoing on each actual
    // reflected pair. This remains valid for flux/film/insulated walls and
    // does not reconstruct exterior heat from the volume balance residual.
    // Its trace is sampled before the second source kick, which is not wall heat.
    let wall_exchange = candidate.wall_exchange().ok_or(invalid("missing thermal wall exchange"))?;
    for (wall, exchange) in wall_exchange.iter().enumerate() {
        if !exchange.mean_temperature.is_finite()
            || !exchange.min_temperature.is_finite() || !exchange.max_temperature.is_finite()
            || exchange.min_temperature < limits.temperature_range[0]
            || exchange.max_temperature > limits.temperature_range[1]
        {
            return Err(ThermalSourceError::WallTemperature {
                wall, minimum: exchange.min_temperature, maximum: exchange.max_temperature,
            }.into());
        }
    }
    let [bottom_wall, top_wall] = wall_exchange.map(|exchange| exchange.heat_in);
    candidate.source_half_kick(heat, plane);
    check_temperatures(candidate.populations(), plane, limits.temperature_range)?;
    let stored_change: f64 = candidate.populations()[plane..cells-plane].iter()
        .zip(&state.populations()[plane..cells-plane])
        .map(|(after, before)| after.iter().zip(before).map(|(a,b)| a-b).sum::<f64>()).sum();
    let defect = stored_change - source - bottom_wall - top_wall;
    if !defect.is_finite() || defect.abs() > limits.balance_tolerance {
        return Err(ThermalSourceError::Balance { defect }.into());
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

// Both lattices use the same even/odd TRT relaxation law. Reuse the native
// checked model's admission rather than a second formula for the odd rate.
fn trt_rates_valid(tau: &[f64], magic: f64) -> bool {
    tau.iter().all(|&tau| CollisionModel3::Trt { tau, magic }.validate().is_ok())
}

impl Slab<2, 5> for ThermalLbm {
    fn layout(&self) -> Result<(usize, usize), ThermalSourceError> {
        let g = &self.grid;
        let n = g.nx.checked_mul(g.ny).ok_or(invalid("cell count overflows"))?;
        if g.nx == 0 || g.ny < 4 || [g.f.len(), g.flags.len(), g.tau.len(), g.fext.len(), self.gpop.len()].iter().any(|l| *l != n)
            || self.spatial_wall_boundaries.as_ref().is_some_and(|walls| walls.iter().any(|wall| wall.len() != g.nx))
        {
            return Err(invalid("invalid 2D slab storage"));
        }
        Ok((n, g.nx))
    }
    fn populations(&self) -> &[[f64; 5]] { &self.gpop }
    fn source_half_kick(&mut self, heat: &[f64], plane: usize) {
        let end = self.gpop.len() - plane;
        let grid = &self.grid;
        super::advection::D2Q5.source_half_kick(
            &mut self.gpop, heat, plane..end, &mut self.advection_history,
            |cell| grid.moments(cell).u,
        );
    }
    fn snapshot(&self) -> (Self, [f64; 5]) {
        (Self {
            grid: self.grid.clone(), gpop: self.gpop.clone(), tau_g: self.tau_g,
            t_bottom: self.t_bottom, t_top: self.t_top, gbeta: self.gbeta,
            t_ref: self.t_ref, scratch: Vec::new(),
            buoyancy_direction: self.buoyancy_direction,
            buoyancy_quadrature: self.buoyancy_quadrature,
            coriolis: self.coriolis,
            wall_boundaries: self.wall_boundaries, wall_exchange: self.wall_exchange,
            spatial_wall_boundaries: self.spatial_wall_boundaries.clone(),
            advection_history: self.advection_history.clone(),
        }, super::W5)
    }
    fn walls(&self) -> (f64, f64) { (self.t_bottom, self.t_top) }
    fn thermal_wall_fields(&self) -> [ThermalWallField<'_>; 2] { self.thermal_wall_fields() }
    fn wall_exchange(&self) -> Option<[ThermalWallExchange; 2]> { self.wall_exchange }
    fn flow_valid(&self) -> bool {
        let g = &self.grid;
        flow_valid(&g.f, &g.tau, &g.fext, self.tau_g) && self.gbeta.is_finite() && self.t_ref.is_finite()
            && g.periodic_x && !g.periodic_y && self.momentum_relaxation_valid(&g.tau)
            && g.g.iter().all(|x| x.is_finite())
            && g.flags.iter().enumerate().all(|(i, flag)| *flag == if i < g.nx || i >= g.f.len()-g.nx { Cell::Wall } else { Cell::Fluid })
    }
    fn advance(&mut self, forces: &[[f64; 2]]) -> Result<(), ThermalForcingError> { self.step_with_external_forces(forces) }
    fn momentum_relaxation(&self) -> &[f64] { &self.grid.tau }
    fn momentum_relaxation_mut(&mut self) -> &mut [f64] { &mut self.grid.tau }
    fn momentum_relaxation_valid(&self, tau: &[f64]) -> bool {
        match self.grid.collision {
            CollisionModel2::Bgk => true,
            CollisionModel2::Trt { magic } => trt_rates_valid(tau, magic),
            _ => false,
        }
    }
    fn advance_local(&mut self, forces: &[[f64; 2]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError> {
        self.step_with_local_thermal_relaxation(forces, thermal_tau)
    }
    fn advance_local_resolved<E: From<ThermalForcingError>>(
        &mut self, forces: &[[f64; 2]], thermal_tau: &[f64],
        resolve: impl FnOnce(&Self) -> Result<Vec<[f64; 2]>, E>,
    ) -> Result<(), E> {
        self.step_with_resolved_external_forces(forces, thermal_tau, resolve)
    }
    fn advance_tensor(&mut self, forces: &[[f64; 2]], diffusion: super::tensor::DiffusivityField<'_, 2>) -> Result<(), ThermalForcingError> {
        self.step_with_tensor_diffusion(forces, diffusion)
    }
}
impl Slab<3, 7> for ThermalLbm3 {
    fn layout(&self) -> Result<(usize, usize), ThermalSourceError> {
        let g = &self.grid;
        let plane = g.nx.checked_mul(g.ny).ok_or(invalid("wall plane overflows"))?;
        let n = plane.checked_mul(g.nz).ok_or(invalid("cell count overflows"))?;
        if plane == 0 || g.nz < 4 || [g.f.len(), g.tau.len(), g.fext.len(), self.gpop.len()].iter().any(|l| *l != n)
            || self.spatial_wall_boundaries.as_ref().is_some_and(|walls| walls.iter().any(|wall| wall.len() != plane))
        {
            return Err(invalid("invalid 3D slab storage"));
        }
        Ok((n, plane))
    }
    fn populations(&self) -> &[[f64; 7]] { &self.gpop }
    fn source_half_kick(&mut self, heat: &[f64], plane: usize) {
        let end = self.gpop.len() - plane;
        let grid = &self.grid;
        super::advection::D3Q7.source_half_kick(
            &mut self.gpop, heat, plane..end, &mut self.advection_history,
            |cell| grid.moments3(cell).1,
        );
    }
    fn snapshot(&self) -> (Self, [f64; 7]) {
        let g = &self.grid;
        let mut s = Self::slab(g.nx, g.ny, g.nz-2, 1.0, self.tau_g, self.gbeta);
        let weights = s.gpop[0];
        s.grid.f.clone_from(&g.f); s.grid.tau.clone_from(&g.tau); s.grid.fext.clone_from(&g.fext);
        s.grid.trt_magic = g.trt_magic();
        s.grid.central_moment_bulk.clone_from(&g.central_moment_bulk);
        s.grid.copy_plate_state_from(g);
        s.gpop.clone_from(&self.gpop);
        s.advection_history.clone_from(&self.advection_history);
        s.t_bottom = self.t_bottom; s.t_top = self.t_top; s.t_ref = self.t_ref;
        s.buoyancy_direction = self.buoyancy_direction;
        s.buoyancy_quadrature = self.buoyancy_quadrature;
        s.coriolis = self.coriolis;
        s.wall_boundaries = self.wall_boundaries; s.wall_exchange = self.wall_exchange;
        s.spatial_wall_boundaries.clone_from(&self.spatial_wall_boundaries);
        (s, weights)
    }
    fn walls(&self) -> (f64, f64) { (self.t_bottom, self.t_top) }
    fn thermal_wall_fields(&self) -> [ThermalWallField<'_>; 2] { self.thermal_wall_fields() }
    fn wall_exchange(&self) -> Option<[ThermalWallExchange; 2]> { self.wall_exchange }
    fn flow_valid(&self) -> bool {
        flow_valid(&self.grid.f, &self.grid.tau, &self.grid.fext, self.tau_g) && self.gbeta.is_finite() && self.t_ref.is_finite()
            && self.momentum_relaxation_valid(&self.grid.tau)
    }
    fn advance(&mut self, forces: &[[f64; 3]]) -> Result<(), ThermalForcingError> { self.step_with_external_forces(forces) }
    fn momentum_relaxation(&self) -> &[f64] { &self.grid.tau }
    fn momentum_relaxation_mut(&mut self) -> &mut [f64] { &mut self.grid.tau }
    fn momentum_relaxation_valid(&self, tau: &[f64]) -> bool {
        if self.grid.has_independent_bulk_collision() {
            return self.grid.bulk_field_length_valid(tau.len())
                && tau.iter().enumerate().all(|(cell, &value)| self.grid.collision_with_tau_at(cell, value).validate().is_ok());
        }
        self.grid.trt_magic().is_none_or(|magic| trt_rates_valid(tau, magic))
    }
    fn advance_local(&mut self, forces: &[[f64; 3]], thermal_tau: &[f64]) -> Result<(), ThermalForcingError> {
        self.step_with_local_thermal_relaxation(forces, thermal_tau)
    }
    fn advance_local_resolved<E: From<ThermalForcingError>>(
        &mut self, forces: &[[f64; 3]], thermal_tau: &[f64],
        resolve: impl FnOnce(&Self) -> Result<Vec<[f64; 3]>, E>,
    ) -> Result<(), E> {
        self.step_with_resolved_external_forces(forces, thermal_tau, resolve)
    }
    fn advance_tensor(&mut self, forces: &[[f64; 3]], diffusion: super::tensor::DiffusivityField<'_, 3>) -> Result<(), ThermalForcingError> {
        self.step_with_tensor_diffusion(forces, diffusion)
    }
}

impl ThermalLbm {
    /// Apply prescribed signed temperature increments and external forces for
    /// one step. Includes walls in each input array; wall heating must be zero.
    /// Both temperature half-kicks, flow, heat transport and the measured wall
    /// balance must pass before the candidate replaces the original state.
    /// This checked path supports the periodic slab with BGK or TRT collision.
    /// Edited obstacle masks and other collision operators remain unsupported.
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

#[cfg(test)]
mod history_tests {
    use super::*;

    fn reject_then_advance<S: Slab<D, Q>, const D: usize, const Q: usize>(
        mut actual: S, mut reference: S, force: [f64; D],
    ) {
        let (cells, plane) = actual.layout().unwrap();
        let forces = vec![force; cells];
        let zeros = vec![0.0; cells];
        let limits = ThermalSourceLimits {
            max_cells: cells,
            temperature_range: [-2.0, 2.0],
            balance_tolerance: 1e-10,
        };
        // Warm two independently constructed states through the ordinary
        // transport path. Using snapshot() to construct the oracle would
        // hide a missing history copy in the very operation being tested.
        for _ in 0..5 {
            actual.advance(&forces).unwrap();
            reference.advance(&forces).unwrap();
        }
        let before = actual.populations().to_vec();
        let before_walls = actual.thermal_walls();
        let before_exchange = actual.wall_exchange();
        let mut heat = vec![0.0; cells];
        for cell in plane..cells - plane {
            heat[cell] = 0.0002 * ((cell % 7) as f64 - 3.0);
        }
        // Reject only after BOTH source kicks, thermal/flow advancement and
        // the independent balance have succeeded. The candidate has already
        // changed its history as well as its visible populations.
        let rejected = source_step_transport_mapped::<S, D, Q, ()>(
            &mut actual, &heat, &forces, limits, None,
            |_, _| Err(invalid("test endpoint rejection")),
        );
        assert_eq!(rejected, Err(invalid("test endpoint rejection")));
        assert_eq!(actual.populations(), before.as_slice());
        assert_eq!(actual.thermal_walls(), before_walls);
        assert_eq!(actual.wall_exchange(), before_exchange);
        // The next accepted source steps must retain the same warmed
        // derivative as the ordinary solver, not restart it in each clone.
        // This catches both rollback leaks and omission from snapshot().
        for _ in 0..4 {
            reference.advance(&forces).unwrap();
            source_step_mapped::<S, D, Q, ()>(
                &mut actual, &zeros, &forces, limits, |_| Ok(()),
            ).unwrap();
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(actual.wall_exchange(), reference.wall_exchange());
        }
    }

    #[test]
    fn g0_rejected_slab_candidate_preserves_history_for_next_accepted_step() {
        let mut two = ThermalLbm::slab(7, 10, 0.9, 0.8, 1e-5);
        let mut two_reference = ThermalLbm::slab(7, 10, 0.9, 0.8, 1e-5);
        two.perturb(0.012);
        two_reference.perturb(0.012);
        reject_then_advance::<_, 2, 5>(two, two_reference, [1e-5, -2e-6]);

        let mut three = ThermalLbm3::slab(5, 3, 10, 0.9, 0.8, 1e-5);
        let mut three_reference = ThermalLbm3::slab(5, 3, 10, 0.9, 0.8, 1e-5);
        three.perturb(0.012);
        three_reference.perturb(0.012);
        reject_then_advance::<_, 3, 7>(three, three_reference, [1e-5, -2e-6, 3e-6]);
    }

    #[test]
    fn g0_flux_and_film_laws_survive_source_snapshots_and_endpoint_rejection() {
        let walls = [
            ThermalWall::HeatFlux(0.0002),
            ThermalWall::Convective { transfer: 0.04, ambient: 0.0 },
        ];
        let mut two = ThermalLbm::slab(7, 10, 0.9, 0.8, 1e-5);
        let mut two_reference = ThermalLbm::slab(7, 10, 0.9, 0.8, 1e-5);
        two.set_thermal_walls(walls).unwrap();
        two_reference.set_thermal_walls(walls).unwrap();
        reject_then_advance::<_, 2, 5>(two, two_reference, [1e-5, -2e-6]);

        let mut three = ThermalLbm3::slab(5, 3, 10, 0.9, 0.8, 1e-5);
        let mut three_reference = ThermalLbm3::slab(5, 3, 10, 0.9, 0.8, 1e-5);
        three.set_thermal_walls(walls).unwrap();
        three_reference.set_thermal_walls(walls).unwrap();
        reject_then_advance::<_, 3, 7>(three, three_reference, [1e-5, -2e-6, 3e-6]);
    }

    #[test]
    fn g4_local_wall_trace_refuses_even_when_its_mean_and_fluid_are_in_range() {
        let limits = ThermalSourceLimits {
            max_cells: 8, temperature_range: [0.0, 1.0], balance_tolerance: 1e-12,
        };
        let walls = [ThermalWall::Adiabatic, ThermalWall::Temperature(0.0)];
        let mut two = ThermalLbm::slab(2, 2, 0.9, 0.8, 0.0);
        two.set_thermal_walls(walls).unwrap();
        let cell = two.grid.idx(0, 1);
        // Retain cell temperature while changing its nonequilibrium normal
        // flux. After collision one insulated wall trace is negative, although
        // both fluid temperatures and the mean lower-wall trace remain in [0,1].
        two.gpop[cell][4] += 0.8;
        two.gpop[cell][0] -= 0.8;
        let before = two.gpop.clone();
        let refusal = two.step_with_sources(&[0.0; 8], &[[0.0; 2]; 8], limits);
        assert!(matches!(refusal,
            Err(ThermalSourceError::WallTemperature { wall: 0, minimum, maximum })
                if minimum < 0.0 && maximum <= 1.0));
        assert_eq!(two.gpop, before);
        assert_eq!(two.last_wall_exchange(), None);

        let mut three = ThermalLbm3::slab(2, 1, 2, 0.9, 0.8, 0.0);
        three.set_thermal_walls(walls).unwrap();
        let cell = three.grid.idx(0, 0, 1);
        three.gpop[cell][6] += 0.6;
        three.gpop[cell][0] -= 0.6;
        let before = three.gpop.clone();
        let refusal = three.step_with_sources(&[0.0; 8], &[[0.0; 3]; 8], limits);
        assert!(matches!(refusal,
            Err(ThermalSourceError::WallTemperature { wall: 0, minimum, maximum })
                if minimum < 0.0 && maximum <= 1.0));
        assert_eq!(three.gpop, before);
        assert_eq!(three.last_wall_exchange(), None);
    }
}
