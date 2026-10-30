//! Prepared source-viscosity feedback without rebuilding inertia or erasing fluid history.
//!
//! At fixed density, radius and dimensionless quadrature, cylinder spectral
//! masses do not depend on viscosity: c_j/lambda_j is constant. The ORIGINAL
//! stored-energy matrix therefore survives a temperature-dependent viscosity.
//! Only positive physical resistance factors change. The native prepared
//! Gonzalez owner advances all coordinates together, using its exact analytic
//! tangent and optional independent-memory condensation. This is not a new
//! integrator, mode retuning, diagonal damping approximation or history reset.

use super::*;
use fs_matdb::EvaluationDecision;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::liquid::{DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_PROPERTY};
use fs_material::state_point::{DENSITY_PROPERTY, MaterialStatePointError,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY};
use fs_phs::{PreparedStepError, PreparedStepRecord, StepWorkspace};
use fs_qty::Temperature;

/// Complete native candidate and physical losses, not an accepted clock update.
/// The original model still owns every state coordinate and source section.
pub struct PreparedFrameFluidStep {
    pub state: Vec<f64>,
    pub port_outputs: Vec<f64>,
    pub native: PreparedStepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub material_port_dissipation_j: Vec<f64>,
    pub fluid_member_dissipation_j: Vec<f64>,
    pub nodal_work_j: f64,
    pub energy_residual_j: f64,
    /// Original native (condensed, full dense) solve counts, not a timing claim.
    pub linear_solve_counts: (usize, usize),
}

/// Reusable native scratch and the SAME original frame/fluid Hamiltonian.
/// Only this scratch mutates during a proposal; no physical state is retained
/// here. Coefficient/source queries and returned records may allocate. Native
/// stepping and resistance/Hessian callbacks allocate nothing after preparation.
pub struct PreparedFluidFrame<'f, 's, F> {
    model: &'f MaterialFluidFrame<'s, F>,
    lossless: PortHamiltonian,
    workspace: StepWorkspace,
    input: Vec<f64>,
    next: Vec<f64>,
    output: Vec<f64>,
    unit_ratios: Vec<f64>,
    material_resistances: Vec<f64>,
}

impl<'s, F> MaterialFluidFrame<'s, F> {
    /// Cache the existing implicit solver and exact fluid-memory elimination.
    /// Every physical state remains present. The native owner falls back to
    /// full dense LU if its exact elimination checks refuse the current layout.
    /// Original state/dense budgets bound each matrix, not total workspace bytes.
    pub fn prepare(&self, gate: &CancelGate) -> Result<PreparedFluidFrame<'_, 's, F>> {
        poll(gate)?;
        let n = self.system.state_dim();
        let m = self.basis.modes().len();
        let q = self.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("original frame/fluid quadratic storage declaration is absent"))?;
        let (j, _, g) = self.system.structure();
        // All original losses are supplied below as their physical PSD factors,
        // so the frozen dense R must NOT also be added to the state equation.
        let storage = QuadraticStorage::new(q.to_vec(), n)?;
        let lossless = PortHamiltonian::new(n, m, j.to_vec(), vec![0.0; n*n],
            g.to_vec(), Box::new(storage))?;
        poll(gate)?;
        let mut workspace = StepWorkspace::new(&lossless)?;
        let pairs: Vec<_> = (2*m..n.saturating_sub(1)).step_by(2)
            .map(|i| [i, i+1]).collect();
        workspace.set_condensed_pairs(&pairs)?;
        poll(gate)?;
        Ok(PreparedFluidFrame {
            model: self, lossless, workspace, input: vec![0.0; m],
            next: vec![0.0; n], output: vec![0.0; m],
            unit_ratios: vec![1.0; self.members.len()],
            material_resistances: self.material_ports.iter().map(|p| p.resistance).collect(),
        })
    }
}

impl<'f, 's, F> PreparedFluidFrame<'f, 's, F> {
    pub const fn model(&self) -> &'f MaterialFluidFrame<'s, F> { self.model }
    /// Actual native dense border after exact independent-memory elimination.
    pub fn condensed_dimension(&self) -> Option<usize> { self.workspace.condensed_dimension() }

    /// Repeated frozen-source stepping through the original physical model.
    pub fn step_candidate(&mut self, state: &[f64], nodal_force: &[f64], dt: f64,
        energy_tolerance_j: f64, gate: &CancelGate) -> Result<PreparedFrameFluidStep> {
        // Borrow fields separately: no per-step clone of the baseline factors.
        solve(self.model, &self.lossless, &mut self.workspace, &mut self.input,
            &mut self.next, &mut self.output, &self.unit_ratios, &self.material_resistances, state, nodal_force,
            dt, energy_tolerance_j, gate)
    }

    /// Crate-private seam for the original layer/caloric owner. These are
    /// source-derived port resistances, never user-authored material overrides.
    pub(crate) fn step_with_material_resistances(&mut self, material_resistances: &[f64],
        state: &[f64], nodal_force: &[f64], dt: f64, energy_tolerance_j: f64,
        gate: &CancelGate) -> Result<PreparedFrameFluidStep> {
        solve(self.model, &self.lossless, &mut self.workspace, &mut self.input,
            &mut self.next, &mut self.output, &self.unit_ratios, material_resistances,
            state, nodal_force, dt, energy_tolerance_j, gate)
    }

    // Only the pinned liquid source owner below supplies non-unit factors.
    // In particular, no public raw damping override acquires source authority.
    fn step_with_ratios(&mut self, ratios: &[f64], state: &[f64], nodal_force: &[f64],
        dt: f64, energy_tolerance_j: f64, gate: &CancelGate) -> Result<PreparedFrameFluidStep> {
        solve(self.model, &self.lossless, &mut self.workspace, &mut self.input,
            &mut self.next, &mut self.output, ratios, &self.material_resistances, state, nodal_force,
            dt, energy_tolerance_j, gate)
    }
}

#[allow(clippy::too_many_arguments)]
fn solve<F>(model: &MaterialFluidFrame<'_, F>, system: &PortHamiltonian,
    work: &mut StepWorkspace, input: &mut [f64], next: &mut [f64], output: &mut [f64],
    ratios: &[f64], material_resistances: &[f64], state: &[f64], nodal_force: &[f64], dt: f64,
    energy_tolerance_j: f64, gate: &CancelGate) -> Result<PreparedFrameFluidStep> {
    poll(gate)?;
    positive(dt)?;
    positive(energy_tolerance_j)?;
    model.admit_state(state)?;
    if ratios.len() != model.members.len() || ratios.iter().any(|r| !r.is_finite() || *r <= 0.0) {
        return Err(invalid("one positive finite source viscosity ratio per original member required"));
    }
    if material_resistances.len() != model.material_ports.len()
        || material_resistances.iter().any(|c| !c.is_finite() || *c < 0.0) {
        return Err(invalid("complete finite nonnegative original material-port resistances required"));
    }
    for station in &model.stations {
        poll(gate)?;
        let weight = positive(station.weight_m*ratios[station.member])?;
        let flow = &model.members[station.member].flow;
        positive(weight*flow.low_tail_resistance_per_length_n_s_m2())?;
        for branch in flow.branches() { positive(weight*branch.resistance_per_length_n_s_m2)?; }
    }
    if let Some(damping) = model.basis.material_damping() { damping.admit_timestep(dt)?; }
    model.basis.project_force_into(nodal_force, [0.0; 2], input, || gate.is_requested())?;
    let m = model.basis.modes().len();
    let n = system.state_dim();
    let q = system.storage().exact_quadratic_matrix()
        .ok_or_else(|| invalid("prepared frame/fluid storage declaration is absent"))?;
    let force = |_: &[f64], effort: &[f64], out: &mut [f64]| {
        resistance_action(model, ratios, material_resistances, effort, out)
    };
    let tangent = |_: &[f64], _: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        resistance_action(model, ratios, material_resistances, de, out)
    };
    let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
        // The declared original Q is diagonal except for the loaded modal
        // momentum block. Preserve that FULL block without n*n zero products.
        for i in 0..m {
            out[i] = q[i*n+i]*direction[i];
            out[m+i] = (0..m).map(|j| q[(m+i)*n+m+j]*direction[m+j]).sum();
        }
        for i in 2*m..n { out[i] = q[i*n+i]*direction[i]; }
        out.iter().all(|v| v.is_finite())
    };
    let attempt = work.step_into_dissipative_controlled(system, state, input, dt, next, output,
        &force, Some(&hessian), Some(&tangent), || gate.is_requested());
    poll(gate)?;
    let native = attempt.map_err(|e| match e {
        PreparedStepError::Cancelled => FrameFluidError::Cancelled,
        PreparedStepError::Solver(e) => FrameFluidError::Native(e),
    })?;
    model.admit_state(next)?;
    let effort = work.last_step_discrete_effort()
        .ok_or_else(|| invalid("native prepared frame effort is absent"))?;
    let mut material_port_dissipation_j = Vec::with_capacity(model.material_ports.len());
    for (port, &resistance) in model.material_ports.iter().zip(material_resistances) {
        poll(gate)?;
        let v = dot(&port.projection, &effort[m..2*m])?;
        material_port_dissipation_j.push(finite(dt*resistance*v*v)?);
    }
    let mut fluid_member_dissipation_j = vec![0.0; model.members.len()];
    for station in &model.stations {
        poll(gate)?;
        let v = dot(&station.projection, &effort[m..2*m])?;
        let flow = &model.members[station.member].flow;
        let weight = finite(station.weight_m*ratios[station.member])?;
        let mut power = weight*flow.low_tail_resistance_per_length_n_s_m2()*v*v;
        for (k, branch) in flow.branches().iter().enumerate() {
            let slip = v-effort[2*m+station.first_memory+k];
            power += weight*branch.resistance_per_length_n_s_m2*slip*slip;
        }
        fluid_member_dissipation_j[station.member] = finite(
            fluid_member_dissipation_j[station.member]+dt*power)?;
    }
    let mut nodal_work_j = 0.0;
    for (dof, f) in nodal_force.iter().enumerate() {
        if dof%128 == 0 { poll(gate)?; }
        let du: f64 = model.basis.modes().iter().enumerate()
            .map(|(i,mode)| mode.displacement[dof]*(next[i]-state[i])).sum();
        nodal_work_j = finite(nodal_work_j+f*du)?;
    }
    let stored_before_j = model.system.hamiltonian(state);
    let stored_after_j = model.system.hamiltonian(next);
    let loss = finite(material_port_dissipation_j.iter().chain(&fluid_member_dissipation_j).sum())?;
    let energy_residual_j = finite(stored_after_j-stored_before_j+loss-nodal_work_j)?;
    for residual_j in [energy_residual_j, finite(native.supplied-nodal_work_j)?,
        finite(native.dissipated-loss)?, finite(native.delta_h-(stored_after_j-stored_before_j))?] {
        if residual_j.abs() > energy_tolerance_j {
            return Err(FrameFluidError::Balance {residual_j, tolerance_j: energy_tolerance_j});
        }
    }
    poll(gate)?;
    Ok(PreparedFrameFluidStep { state: next.to_vec(), port_outputs: output.to_vec(), native,
        stored_before_j, stored_after_j, material_port_dissipation_j,
        fluid_member_dissipation_j, nodal_work_j, energy_residual_j,
        linear_solve_counts: work.linear_solve_counts() })
}

// One physical factor action serves both nonlinear residual and exact analytic
// tangent. Each station evaluates its body velocity once, accumulates all slip
// forces, then projects their SUM once: O(stations*(m+branches)), not O(n*n).
fn resistance_action<F>(model: &MaterialFluidFrame<'_, F>, ratios: &[f64],
    material_resistances: &[f64], effort: &[f64], out: &mut [f64]) -> bool {
    let m = model.basis.modes().len();
    out.fill(0.0);
    for (port, &resistance) in model.material_ports.iter().zip(material_resistances) {
        let v: f64 = port.projection.iter().zip(&effort[m..2*m]).map(|(b,e)| b*e).sum();
        let f = resistance*v;
        for i in 0..m { out[m+i] += port.projection[i]*f; }
    }
    for station in &model.stations {
        let flow = &model.members[station.member].flow;
        let v: f64 = station.projection.iter().zip(&effort[m..2*m]).map(|(b,e)| b*e).sum();
        let weight = station.weight_m*ratios[station.member];
        let mut reaction = weight*flow.low_tail_resistance_per_length_n_s_m2()*v;
        for (k, branch) in flow.branches().iter().enumerate() {
            let index = 2*m+station.first_memory+k;
            let f = weight*branch.resistance_per_length_n_s_m2*(v-effort[index]);
            out[index] -= f;
            reaction += f;
        }
        for i in 0..m { out[m+i] += station.projection[i]*reaction; }
    }
    out.iter().all(|v| v.is_finite())
}

#[derive(Debug)]
pub enum LiquidFrameError {
    Mechanical(FrameFluidError),
    Source(MaterialStatePointError),
    Invalid(&'static str),
}
impl core::fmt::Display for LiquidFrameError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "liquid frame: {self:?}") }
}
impl std::error::Error for LiquidFrameError {}
impl From<FrameFluidError> for LiquidFrameError {
    fn from(e: FrameFluidError) -> Self { Self::Mechanical(e) }
}
impl From<MaterialStatePointError> for LiquidFrameError {
    fn from(e: MaterialStatePointError) -> Self { Self::Source(e) }
}

/// Existing fixed-density liquid laws on the original fluid-memory coordinates.
/// Temperatures are prescribed here. The finite thermal owner, not this class,
/// must account for heat and determine temperature in a closed simulation.
pub struct LiquidViscosityFrame<'f, 's, 'a> {
    prepared: PreparedFluidFrame<'f, 's, ResolvedLiquidState>,
    laws: Vec<LiquidTemperatureLaw<'a>>,
}

pub struct LiquidFrameStep {
    pub mechanical: PreparedFrameFluidStep,
    /// Exact original source claims evaluated at this proposal's temperatures.
    pub liquids: Vec<ResolvedLiquidState>,
}

impl<'f, 's, 'a> LiquidViscosityFrame<'f, 's, 'a> {
    /// Original source and order must match EVERY immersed member. Scalar,
    /// validity-wide constant density is required; equality of two samples does
    /// not establish constant inertia. No density tolerance or mass remapping.
    pub fn new(model: &'f MaterialFluidFrame<'s, ResolvedLiquidState>,
        laws: Vec<LiquidTemperatureLaw<'a>>, gate: &CancelGate)
        -> core::result::Result<Self, LiquidFrameError> {
        poll(gate)?;
        if laws.len() != model.members.len() {
            return Err(LiquidFrameError::Invalid("one original liquid law per loaded member required"));
        }
        for (member, law) in model.members.iter().zip(&laws) {
            poll(gate)?;
            if member.flow.fluid() != law.reference_state() {
                return Err(LiquidFrameError::Invalid("liquid law must retain the exact original prepared source bundle"));
            }
            let density = law.reference_state().material().property(DENSITY_PROPERTY)
                .ok_or(LiquidFrameError::Invalid("original liquid density source is absent"))?;
            if density.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
                return Err(LiquidFrameError::Invalid("viscosity feedback requires an explicitly constant-density source"));
            }
            let viscosity = match law.viscosity_basis() {
                LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_PROPERTY,
                LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_PROPERTY,
            };
            for name in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY,
                THERMAL_CONDUCTIVITY_PROPERTY, viscosity] {
                let property = law.reference_state().material().property(name)
                    .ok_or(LiquidFrameError::Invalid("original required transport source is absent"))?;
                if !property.is_static_at_fixed_external_conditions() {
                    return Err(LiquidFrameError::Invalid("fixed-context Newtonian transport cannot freeze unowned rate or frequency dependence"));
                }
            }
        }
        let prepared = model.prepare(gate)?;
        Ok(Self {prepared, laws})
    }
    pub fn model(&self) -> &'f MaterialFluidFrame<'s, ResolvedLiquidState> { self.prepared.model() }
    pub fn laws(&self) -> &[LiquidTemperatureLaw<'a>] { &self.laws }
    pub fn condensed_dimension(&self) -> Option<usize> { self.prepared.condensed_dimension() }

    /// Re-evaluate the SAME claims, retaining density, momentum and all memory.
    /// Missing source support refuses without selecting a replacement claim.
    pub fn liquids_at_temperatures(&self, temperatures_k: &[f64], gate: &CancelGate)
        -> core::result::Result<Vec<ResolvedLiquidState>, LiquidFrameError> {
        poll(gate)?;
        if temperatures_k.len() != self.laws.len() {
            return Err(LiquidFrameError::Invalid("complete liquid temperature vector required"));
        }
        let mut states = Vec::with_capacity(self.laws.len());
        for (law, &t) in self.laws.iter().zip(temperatures_k) {
            poll(gate)?;
            let state = law.at_temperature(Temperature::new(t))?;
            if state.density_kg_m3().to_bits() != law.reference_state().density_kg_m3().to_bits() {
                return Err(LiquidFrameError::Invalid("liquid temperature update changed original fluid inertia"));
            }
            states.push(state);
        }
        poll(gate)?;
        Ok(states)
    }

    /// Hold the proposed source viscosities over ONE original implicit step.
    /// No eigenvectors, inverse mass, Bessel quadrature or fluid coordinates are
    /// rebuilt. Warm native algebra retains its original iteration/budget gates.
    pub fn step_at_temperatures(&mut self, state: &[f64], nodal_force: &[f64],
        temperatures_k: &[f64], dt: f64, energy_tolerance_j: f64, gate: &CancelGate)
        -> core::result::Result<LiquidFrameStep, LiquidFrameError> {
        let liquids = self.liquids_at_temperatures(temperatures_k, gate)?;
        let ratios: Vec<_> = liquids.iter().zip(self.model().members())
            .map(|(next, original)| positive(next.dynamic_viscosity_pa_s()
                / original.flow().transport().dynamic_viscosity_pa_s()))
            .collect::<Result<_>>()?;
        let mechanical = self.prepared.step_with_ratios(&ratios, state, nodal_force,
            dt, energy_tolerance_j, gate)?;
        poll(gate)?;
        Ok(LiquidFrameStep {mechanical, liquids})
    }
}

#[cfg(test)]
pub(crate) mod tests;
