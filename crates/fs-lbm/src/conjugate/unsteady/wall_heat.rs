//! Equilibrium wall transport refreshed by the original coupled time marcher.
//! No momentum solver, thermal integrator, accepted state or clock is duplicated.

use fs_exec::CancelGate;
use super::{Boussinesq, ConjugateMarch, UnsteadyConfig, march_conjugate_inner};
use crate::conjugate::{ChtError, EnergyConfig, FlowField, FluidProperties,
    FvBoundary, RadiationConfig, SimpleConfig, SolidMaterial, ThermalFace,
    ThermalSetup, Turbulence, VoxelDomain, poll};
use crate::conjugate::simple::Solver;
use crate::conjugate::turbulence::thermal_wall::{BoundaryWallHeat,
    BoundaryWallHeatExchange, PreparedWallHeatBoundaries, WallHeatConfig,
    WallHeatFace, WallHeatFaceReport};
use crate::conjugate::turbulence::thermal_wall::boundary::prepare_wall_heat_boundary_fields;
use crate::conjugate::turbulence::thermal_wall::coupling::WallHeatFlow;

/// Explicit equilibrium wall transport on a native LVEL flow timeline.
#[derive(Clone, Debug, PartialEq)]
pub struct WallHeatMarchConfig {
    /// Stationary internal solid/fluid walls, retaining original contacts.
    pub internal_walls: Vec<WallHeatFace>,
    /// Exterior walls; original thermal faces must delegate as Adiabatic.
    pub boundary_walls: Vec<BoundaryWallHeat>,
    /// Physical turbulent Prandtl number and bounded quadrature controls.
    pub transport: WallHeatConfig,
}

/// Bounded per-step summary; full per-face coefficients are retained only for
/// the final accepted endpoint, avoiding O(steps * wall_faces) history storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatStepRecord {
    /// Native endpoint time, s.
    pub time_s: f64,
    /// Signed net exterior wall heat leaving the domain, W. Already counted
    /// once in the corresponding native energy record's source/outflow terms.
    pub boundary_outflow_w: f64,
    /// Imposed inward wall-flux watts included in the native source term.
    /// Unlike authored volume sources, these are not scaled by power_schedule.
    pub prescribed_flux_input_w: f64,
    /// Minimum/maximum exterior fluid-side surface temperature, K; None when
    /// the selection consists solely of internal conjugate walls.
    pub surface_temperature_range_k: Option<(f64, f64)>,
    /// Largest y+ among the internal and exterior wall samples at this step.
    pub maximum_y_plus: f64,
    /// Largest fluid-side wall heat coefficient at this step, W/(m^2 K).
    pub maximum_heat_transfer_w_m2_k: f64,
    /// Actual integrand evaluations used after exact exterior-query caching.
    pub quadrature_evaluations: usize,
}

/// Native coupled trajectory and wall observations from those SAME steps.
#[derive(Clone, Debug, PartialEq)]
pub struct WallHeatMarch {
    /// Original momentum/temperature solution, averages and energy closures.
    pub march: ConjugateMarch,
    /// One small wall-transport summary per accepted native step.
    pub wall_records: Vec<WallHeatStepRecord>,
    /// Signed dt-weighted exterior heat, using the native endpoint quadrature.
    /// This is a disaggregation, not an additional energy loss.
    pub boundary_outflow_j: f64,
    /// Internal wall coefficients used by the final energy step.
    pub final_internal_faces: Vec<WallHeatFaceReport>,
    /// Exterior coefficients, temperatures and watts from the final step.
    pub final_boundary_faces: Vec<BoundaryWallHeatExchange>,
    /// Actual tangent-based turbulent conductivity consumed by the final
    /// energy step, W/(m K). Use this for thermal output, not the unchanged
    /// native flow report's legacy secant-based `eddy_conductivity` helper.
    pub final_eddy_conductivity_w_m_k: Vec<f64>,
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput { field: "unsteady.wall_heat", reason: reason.into() }
}

pub(super) struct Coupling<'a> {
    domain: &'a VoxelDomain,
    fluid: &'a FluidProperties,
    faces: &'a [FvBoundary; 6],
    config: &'a WallHeatMarchConfig,
    velocity: Vec<[f64; 3]>,
    records: Vec<WallHeatStepRecord>,
    boundary_j: f64,
    compensation_j: f64,
    final_internal: Vec<WallHeatFaceReport>,
    final_boundary: Vec<BoundaryWallHeatExchange>,
    final_eddy: Vec<f64>,
}

impl<'a> Coupling<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(domain: &'a VoxelDomain, fluid: &'a FluidProperties, flow: &'a SimpleConfig,
        setup: &ThermalSetup, time: &UnsteadyConfig, config: &'a WallHeatMarchConfig,
        gate: &CancelGate) -> Result<Self, ChtError>
    {
        poll(gate)?;
        // Reject invalid time controls before sizing any step history. The
        // original march still owns full flow and thermal admission.
        super::admit_unsteady(domain, flow, time)?;
        config.transport.validate()?;
        if flow.turbulence != Turbulence::Lvel {
            return Err(invalid("the equilibrium wall-heat march requires explicit LVEL momentum, not SGS or an implicit turbulence switch"));
        }
        if config.internal_walls.is_empty() && config.boundary_walls.is_empty() {
            return Err(invalid("select at least one physical wall"));
        }
        if config.internal_walls.len() > domain.cell_count().saturating_mul(3)
            || config.boundary_walls.len() > domain.cell_count().saturating_mul(6)
            || !setup.eddy_conductivity_w_m_k.is_empty()
        {
            return Err(invalid("wall selection exceeds the face budget or thermal transport already has an owner"));
        }
        Ok(Self { domain, fluid, faces: &flow.faces, config,
            velocity: vec![[0.0; 3]; domain.cell_count()],
            records: Vec::with_capacity(time.steps), boundary_j: 0.0,
            compensation_j: 0.0, final_internal: Vec::new(), final_boundary: Vec::new(),
            final_eddy: Vec::new() })
    }

    /// Borrow the accepted momentum viscosity and the very same projected face
    /// flux passed to energy. Reuse only the cell-velocity scratch, never a
    /// temperature, old wall correction, or a duplicated momentum report.
    pub(super) fn prepare(&mut self, solver: &Solver<'_>, flux: &FlowField,
        setup: &ThermalSetup, gate: &CancelGate) -> Result<PreparedWallHeatBoundaries, ChtError>
    {
        for cell in 0..self.domain.cell_count() {
            if cell % 64 == 0 { poll(gate)?; }
            self.velocity[cell] = if self.domain.is_fluid(cell) {
                solver.cell_velocity(cell)
            } else { [0.0; 3] };
        }
        prepare_wall_heat_boundary_fields(self.domain, self.fluid,
            WallHeatFlow { flux, velocity: &self.velocity, eddy: solver.tangent_eddy_viscosity() },
            setup, self.faces, &self.config.internal_walls, &self.config.boundary_walls,
            self.config.transport, gate)
    }

    pub(super) fn record(&mut self, prepared: PreparedWallHeatBoundaries,
        temperature: &[f64], time_s: f64, dt: f64, gate: &CancelGate) -> Result<(), ChtError>
    {
        let exchanges = prepared.boundary_exchanges(temperature, gate)?;
        let mut outward = 0.0;
        let mut imposed = 0.0;
        let mut range: Option<(f64, f64)> = None;
        let mut max_y: f64 = 0.0;
        let mut max_h: f64 = 0.0;
        for (index, face) in prepared.internal_faces.iter().enumerate() {
            if index % 64 == 0 { poll(gate)?; }
            max_y = max_y.max(face.transfer.y_plus);
            max_h = max_h.max(face.transfer.heat_transfer_w_m2_k);
        }
        for (index, face) in exchanges.iter().enumerate() {
            if index % 64 == 0 { poll(gate)?; }
            outward += face.outward_w;
            if matches!(face.coefficients.boundary.condition, ThermalFace::HeatFlux(_)) {
                imposed -= face.outward_w;
            }
            let t = face.surface_temperature_k;
            range = Some(range.map_or((t, t), |(lo, hi)| (lo.min(t), hi.max(t))));
            max_y = max_y.max(face.coefficients.transfer.y_plus);
            max_h = max_h.max(face.coefficients.transfer.heat_transfer_w_m2_k);
        }
        let increment = dt * outward - self.compensation_j;
        let total = self.boundary_j + increment;
        let compensation = (total - self.boundary_j) - increment;
        if ![outward, imposed, increment, total, compensation].into_iter().all(f64::is_finite) {
            return Err(invalid("integrated boundary heat is not representable"));
        }
        poll(gate)?;
        self.records.push(WallHeatStepRecord { time_s, boundary_outflow_w: outward,
            prescribed_flux_input_w: imposed, surface_temperature_range_k: range,
            maximum_y_plus: max_y, maximum_heat_transfer_w_m2_k: max_h,
            quadrature_evaluations: prepared.quadrature_evaluations });
        self.boundary_j = total;
        self.compensation_j = compensation;
        self.final_internal = prepared.internal_faces;
        self.final_boundary = exchanges;
        self.final_eddy = prepared.thermal.eddy_conductivity_w_m_k;
        Ok(())
    }

    fn finish(self, march: ConjugateMarch) -> WallHeatMarch {
        WallHeatMarch { march, wall_records: self.records, boundary_outflow_j: self.boundary_j,
            final_internal_faces: self.final_internal, final_boundary_faces: self.final_boundary,
            final_eddy_conductivity_w_m_k: self.final_eddy }
    }
}

/// Advance native LVEL flow and temperature with fresh wall heat transport.
/// After each converged and projected momentum step, prepare from the original
/// thermal setup with that step's authored source/radiation contributions, then
/// run the original backward-Euler energy step. Internal contacts, exterior
/// laws and tangent bulk conductivity are each applied once. Sources follow
/// power_schedule; independently prescribed exterior heat flux does not.
///
/// Momentum may use BDF2; the flow/thermal coupling remains first order, with
/// the original lagged Boussinesq force and radiation treatment. This is an
/// evolving EQUILIBRIUM wall model, not a model of separation, roughness or
/// turbulent thermal memory. No existing march opts into it automatically.
///
/// # Errors
/// Original flow/energy admission and convergence errors, incompatible wall
/// selection/transport, nonphysical wall temperatures, or cancellation. The
/// complete result is returned only on success; caller-owned inputs are unchanged.
#[allow(clippy::too_many_arguments)]
pub fn march_conjugate_with_wall_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    walls: &WallHeatMarchConfig,
    buoyancy: Option<&Boussinesq>,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<WallHeatMarch, ChtError> {
    let mut coupling = Coupling::new(domain, fluid, config, setup, unsteady, walls, gate)?;
    let march = march_conjugate_inner(domain, fluid, solids, setup, config, unsteady,
        None, Some(&mut coupling), buoyancy, initial_temperature, inlet_schedule,
        power_schedule, energy, radiation, gate)?;
    poll(gate)?;
    Ok(coupling.finish(march))
}

#[cfg(test)]
mod tests;
