//! Pressure-constrained homogeneous chemistry on the original extent solver.
//!
//! At every trial, conserve H1=H0+Q, retain p1=p0 and derive V1 from the same
//! composition-dependent EOS. Solve xi_r=dt*V1*rate_r(T1,p0,n1), including that
//! evolving volume in both residuals and Jacobian probes. No outer pressure
//! iteration, second kinetics integrator or extra chemical heat source exists.
//! Reference: https://www.cantera.org/3.2/reference/reactors/constant-pressure-reactor.html

/// Same-endpoint thermal contacts with bounded total chemical work.
pub mod thermal;

use super::*;
use boundary::{BoundaryGasReactionStep, GasReactionBoundary};
use fs_thermal::enthalpy::advection::species::close_generated_simplex;

/// Internal storage constraint shared by the original coupled mechanism.
#[derive(Clone, Copy)]
pub(super) enum Constraint {
    Boundary(GasReactionBoundary),
    Pressure { heat_j: f64 },
}
impl Constraint {
    pub(super) fn validate(self) -> Result<(), GasReactionError> {
        match self {
            Self::Boundary(boundary) => boundary.validate(),
            Self::Pressure { heat_j } if heat_j.is_finite() => Ok(()),
            _ => Err(invalid("finite constant-pressure reaction heat required")),
        }
    }

    pub(super) fn reconstruct(
        self, law: &GasReaction, gas: &ClosedGasParcel,
        changes: &[f64], amounts: &[f64], controls: GasReactionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ClosedGasParcel, usize), GasReactionError> {
        let heat_j = match self {
            Self::Boundary(boundary) => return law.reconstruct_inventory_at_boundary(
                gas, changes, amounts, boundary, controls, cancelled,
            ),
            Self::Pressure { heat_j } => heat_j,
        };
        poll(cancelled)?;
        let mut next = gas.clone();
        if changes.iter().all(|&change| change == 0.0) {
            let iterations = if heat_j == 0.0 { 0 }
                else { next.add_heat_at_constant_pressure(heat_j, controls.gas, cancelled)?.iterations };
            return Ok((next, iterations));
        }
        // The existing extent/inventory owner has already formed n0 + S xi.
        // Close only generated simplex roundoff through the original thermal
        // composition helper. Never add a concentration floor or clip amounts.
        let total: f64 = amounts.iter().sum();
        if !positive(total) {
            return Err(invalid("isobaric reaction produced no represented total amount"));
        }
        let mut fractions: Vec<_> = amounts.iter().map(|n| n / total).collect();
        close_generated_simplex(&mut fractions)
            .map_err(|_| invalid("isobaric reaction composition exceeds normalization roundoff"))?;
        if fractions.iter().zip(amounts).any(|(&x, &n)|
            (n > 0.0 && !positive(x)) || (n == 0.0 && x != 0.0))
        {
            return Err(invalid("isobaric composition lost or invented a species"));
        }
        let composition = Composition::new(CompositionBasis::MoleFraction, fractions)
            .map_err(|_| invalid("isobaric composition refused by original quantity owner"))?;
        let model = law.catalog.with_composition(composition)?;
        let report = next.replace_model_with_heat_at_constant_pressure(
            model, heat_j, controls.gas, cancelled,
        )?;
        Ok((next, report.energy.iterations))
    }

    pub(super) fn boundary(
        self, initial: &ClosedGasParcel, endpoint: &ClosedGasParcel,
    ) -> Result<GasReactionBoundary, GasReactionError> {
        match self {
            Self::Boundary(boundary) => Ok(boundary),
            Self::Pressure { heat_j } => {
                let pressure = initial.state().state().pressure;
                if endpoint.state().state().pressure.to_bits() != pressure.to_bits() {
                    return Err(invalid("isobaric chemistry did not retain the original pressure"));
                }
                let boundary = GasReactionBoundary {
                    final_volume_m3: endpoint.volume_m3(), heat_j,
                    work_by_gas_j: pressure * (endpoint.volume_m3() - initial.volume_m3()),
                };
                boundary.validate()?;
                Ok(boundary)
            }
        }
    }
}

/// Original chemical/heat/work equations plus the independent enthalpy balance.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricGasReactionStep {
    pub pressure_pa: f64,
    /// Actual kinetic equations, species inventory, Delta U, p Delta V and entropy.
    pub boundary: BoundaryGasReactionStep,
    /// Total formation-inclusive enthalpy change [J].
    pub enthalpy_change_j: f64,
    /// Delta H minus supplied external heat [J].
    pub enthalpy_residual_j: f64,
}

/// Immutable input, complete candidate output; publication stays with the caller.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricGasReactionCandidate {
    gas: ClosedGasParcel,
    step: IsobaricGasReactionStep,
}
impl IsobaricGasReactionCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub const fn step(&self) -> &IsobaricGasReactionStep { &self.step }
    pub fn into_parts(self) -> (ClosedGasParcel, IsobaricGasReactionStep) { (self.gas, self.step) }
}

impl GasReaction {
    /// Solve all original mechanism channels at ONE constant-pressure endpoint.
    ///
    /// Requires `with_implicit_mechanism`. The initial parcel pressure is the
    /// quasistatic load; heat is signed whole-step external energy, not reaction
    /// enthalpy. Initial mass/species are shared, and the endpoint volume enters
    /// every backward-Euler rate equation. Existing source scales, reverse-rate
    /// thermodynamics, reachable absent species, bounded Newton work, inventory
    /// and heat-entropy admission remain authoritative. Complete gas/kinetic
    /// reports are retained; residual/Jacobian probes omit kinetic receipts.
    ///
    /// Fixed-volume subcycling/scalar policies refuse rather than silently
    /// evolving a different physical constraint. This is first-order closed,
    /// homogeneous ideal-gas evolution, not a flame, shock, constant-pressure
    /// flow reactor, piston trajectory or unconditional stiff-stability claim.
    /// Enclosing thermal reservoirs must account for their own entropy transfer.
    pub fn isobaric_implicit_candidate(
        &self, gas: &ClosedGasParcel, duration_s: f64, heat_j: f64,
        controls: GasReactionSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricGasReactionCandidate, GasReactionError> {
        let policy = self.isobaric_policy(duration_s, heat_j, controls, &mut cancelled)?;
        self.isobaric_candidate_with_policy(gas, duration_s, heat_j, controls, policy, &mut cancelled)
    }

    // Shared admission and numerical owner. A thermal caller may only tighten
    // the borrowed policy's work limit, never change the physical constraint.
    fn isobaric_policy(
        &self, duration_s: f64, heat_j: f64, controls: GasReactionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<&ImplicitMechanismSolve, GasReactionError> {
        poll(cancelled)?;
        controls.validate(duration_s)?;
        Constraint::Pressure { heat_j }.validate()?;
        if self.subcycling.is_some() || self.implicit.is_some() {
            return Err(invalid("constant-pressure chemistry requires one coupled implicit interval"));
        }
        self.implicit_mechanism.as_ref()
            .ok_or_else(|| invalid("constant-pressure chemistry requires the original implicit mechanism"))
    }

    fn isobaric_candidate_with_policy(
        &self, gas: &ClosedGasParcel, duration_s: f64, heat_j: f64,
        controls: GasReactionSolve, policy: &ImplicitMechanismSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IsobaricGasReactionCandidate, GasReactionError> {
        let constraint = Constraint::Pressure { heat_j };
        let candidate = candidate_with_constraint(self, gas, duration_s, constraint, true,
            controls, policy, cancelled)?;
        let (next, reaction) = candidate.into_parts();
        let boundary = constraint.boundary(gas, &next)?;
        let enthalpy = gas.mass_kg() * (next.state().specific_enthalpy_j_kg()
            - gas.state().specific_enthalpy_j_kg());
        let residual = enthalpy - heat_j;
        if !enthalpy.is_finite() || !residual.is_finite()
            || residual.abs() > gas.mass_kg() * controls.gas.energy_tolerance_j_kg
        {
            return Err(invalid("constant-pressure chemistry did not close total enthalpy"));
        }
        let heat_entropy = heat_j / next.state().state().temperature;
        let energy_residual = (reaction.stored_energy_change_j + boundary.work_by_gas_j) - heat_j;
        let entropy_residual = reaction.entropy_change_j_k - heat_entropy;
        poll(cancelled)?;
        Ok(IsobaricGasReactionCandidate {
            gas: next,
            step: IsobaricGasReactionStep {
                pressure_pa: gas.state().state().pressure,
                boundary: BoundaryGasReactionStep {
                    reaction, initial_volume_m3: gas.volume_m3(), boundary,
                    energy_residual_j: energy_residual, endpoint_heat_entropy_j_k: heat_entropy,
                    entropy_residual_j_k: entropy_residual,
                },
                enthalpy_change_j: enthalpy, enthalpy_residual_j: residual,
            },
        })
    }
}
