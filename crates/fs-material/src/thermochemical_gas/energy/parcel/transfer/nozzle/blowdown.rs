//! Adiabatic depletion of a finite, well-mixed, fixed-volume rocket reservoir.
//!
//! With frozen composition, no heat or shaft work and negligible tank kinetic
//! energy, d(m u) = h dm implies m du = R T dm and ds = 0 for the remaining gas.
//! This remains true for temperature-dependent heat capacities. Reconstruct the
//! remaining inventory on that SAME material isentrope and export its lost
//! internal energy once as exhaust stagnation enthalpy. Jet kinetic energy and
//! pressure work are already included; neither is another tank-energy sink.
//!
//! The physical-time update uses implicit midpoint quadrature of nozzle mass
//! rate along this exact thermodynamic path. It is second order for smooth
//! admitted states, not a resolved tank/jet solver or a combustion model.

use super::super::super::super::entropy::{GasEntropyError, GasEntropySolve};
use super::super::{GasSpeciesTransfer, inventory};
use super::rocket::{RocketNozzle, RocketNozzleFlow, RocketNozzleSolve};
use super::{ClosedGasParcel, GasEnergyError, poll, positive, updated_point};
use fs_qty::Temperature;

/// Source-domain, entropy-inversion and extensive-inventory tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasBlowdownSolve {
    /// The same material's reversible density/entropy inversion.
    pub entropy: GasEntropySolve,
    /// Absolute mass discrepancy allowed in total and individual species [kg].
    pub mass_tolerance_kg: f64,
}
impl GasBlowdownSolve {
    fn validate(self) -> Result<(), GasBlowdownError> {
        self.entropy.validate()?;
        if !positive(self.mass_tolerance_kg) {
            return Err(GasBlowdownError::Invalid(
                "positive finite blowdown mass tolerance required",
            ));
        }
        Ok(())
    }
}

/// All failures leave the original reservoir untouched.
#[derive(Clone, Debug, PartialEq)]
pub enum GasBlowdownError {
    /// Invalid input, unrepresentable inventory or a physical step too large.
    Invalid(&'static str),
    /// Original thermodynamic, transport, nozzle or cancellation refusal.
    Energy(GasEnergyError),
    /// Source domain, entropy-inversion or cancellation refusal.
    Entropy(GasEntropyError),
    /// The implicit mass quadrature did not close within the work budget.
    NotConverged {
        iterations: usize,
        mass_residual_kg: f64,
    },
}
impl core::fmt::Display for GasBlowdownError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite gas reservoir blowdown: {self:?}")
    }
}
impl std::error::Error for GasBlowdownError {}
impl From<GasEnergyError> for GasBlowdownError {
    fn from(error: GasEnergyError) -> Self {
        Self::Energy(error)
    }
}
impl From<GasEntropyError> for GasBlowdownError {
    fn from(error: GasEntropyError) -> Self {
        Self::Entropy(error)
    }
}

/// Actual reservoir inventory and exported stagnation-energy accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct GasBlowdownStep {
    /// Mass leaving the finite reservoir [kg].
    pub expelled_mass_kg: f64,
    /// Integral of stagnation enthalpy carried by that exhaust [J]. Signed
    /// formation/reference energy is retained, so this need not be positive.
    /// This is minus the reservoir internal-energy change, not an independent
    /// calorimetric measurement or an extra nozzle kinetic-energy contribution.
    pub exhausted_enthalpy_j: f64,
    /// Mass-average exhaust stagnation enthalpy [J/kg]; None for zero withdrawal.
    pub mean_stagnation_enthalpy_j_kg: Option<f64>,
    /// Remaining-mass specific entropy change, bounded by the solver [J/(kg K)].
    pub specific_entropy_residual_j_kg_k: f64,
    /// Actual remaining plus expelled mass minus initial mass [kg].
    pub mass_residual_kg: f64,
    /// Largest absolute original-species inventory discrepancy [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Exhausted masses in original canonical species order.
    pub species: Vec<GasSpeciesTransfer>,
    /// Interior entropy-inversion iterations.
    pub iterations: usize,
}

/// Tentative endpoint and exhaust flux; publish only after coupled acceptance.
#[derive(Clone, Debug, PartialEq)]
pub struct GasBlowdownCandidate {
    endpoint: ClosedGasParcel,
    step: GasBlowdownStep,
}
impl GasBlowdownCandidate {
    /// Fully resolved remaining gas with original composition and source laws.
    #[must_use]
    pub const fn endpoint(&self) -> &ClosedGasParcel {
        &self.endpoint
    }
    /// Exported species and energy from this same withdrawal.
    #[must_use]
    pub const fn step(&self) -> &GasBlowdownStep {
        &self.step
    }
    /// Move the accepted reservoir and accounting without resampling materials.
    #[must_use]
    pub fn into_parts(self) -> (ClosedGasParcel, GasBlowdownStep) {
        (self.endpoint, self.step)
    }
}

fn remaining(
    old: &ClosedGasParcel,
    expelled: f64,
    controls: GasEntropySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, f64, usize), GasBlowdownError> {
    let mass = old.mass_kg() - expelled;
    if !positive(mass) || mass == old.mass_kg() {
        return Err(GasBlowdownError::Invalid(
            "blowdown must retain a representable positive reservoir mass",
        ));
    }
    let solved = old.model().reconstruct_density_entropy_change(
        old.point(),
        mass / old.volume_m3(),
        0.0,
        controls,
        &mut *cancelled,
    )?;
    let residual = solved.residual_j_kg_k();
    let iterations = solved.iterations();
    let mut parcel = old.parcel.clone();
    parcel.point = updated_point(
        old.point(),
        solved.state().state().temperature,
        solved.state().state().pressure,
    )?;
    parcel.state = solved.into_state();
    parcel.mass_kg = mass;
    Ok((
        ClosedGasParcel::from_isochoric(parcel),
        residual,
        iterations,
    ))
}

impl ClosedGasParcel {
    /// Withdraw a prescribed finite mass adiabatically to an external exhaust.
    /// This thermodynamic operation supplies no rate or downstream-pressure law;
    /// `rocket_blowdown` supplies both. All original species fractions remain
    /// unchanged in a well-mixed reservoir. Source limits, vacuum endpoints,
    /// lost trace-species inventory and cancellation refuse atomically.
    ///
    /// Remaining gas follows ds=0; it does NOT remain isothermal, and its energy
    /// loss is not dm times an unchanged initial enthalpy for a finite withdrawal.
    pub fn exhaust_mass_adiabatic(
        &self,
        mass_kg: f64,
        controls: GasBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasBlowdownCandidate, GasBlowdownError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !mass_kg.is_finite() || mass_kg < 0.0 || mass_kg >= self.mass_kg() {
            return Err(GasBlowdownError::Invalid(
                "finite nonnegative exhaust mass smaller than reservoir inventory required",
            ));
        }
        if mass_kg == 0.0 {
            return Ok(GasBlowdownCandidate {
                endpoint: self.clone(),
                step: GasBlowdownStep {
                    expelled_mass_kg: 0.0,
                    exhausted_enthalpy_j: 0.0,
                    mean_stagnation_enthalpy_j_kg: None,
                    specific_entropy_residual_j_kg_k: 0.0,
                    mass_residual_kg: 0.0,
                    maximum_species_mass_residual_kg: 0.0,
                    species: vec![],
                    iterations: 0,
                },
            });
        }
        let (endpoint, entropy_residual, iterations) =
            remaining(self, mass_kg, controls.entropy, &mut cancelled)?;
        let (t0, t1) = (
            self.state().state().temperature,
            endpoint.state().state().temperature,
        );
        let r = self.state().state().specific_gas_constant;
        let dh = self
            .model()
            .thermodynamics()
            .specific_enthalpy_change(Temperature::new(t0), Temperature::new(t1))
            .map_err(|error| {
                GasEnergyError::from(crate::thermochemical_gas::failure(
                    "blowdown enthalpy",
                    error,
                ))
            })?
            .value();
        let du = dh - r * (t1 - t0);
        let exhausted =
            mass_kg * self.state().specific_internal_energy_j_kg() - endpoint.mass_kg() * du;
        let mean_h =
            self.state().specific_internal_energy_j_kg() - (endpoint.mass_kg() / mass_kg) * du;
        // dU/dm=h. The finite withdrawal's mean h must lie between the endpoint
        // values. Compare relative to h0 so formation energy cannot mask a defect.
        let mean_minus_h0 = -(endpoint.mass_kg() / mass_kg) * du - r * t0;
        let energy_slack =
            4.0 * t0 * controls.entropy.entropy_tolerance_j_kg_k * (self.mass_kg() / mass_kg)
                + 128.0 * f64::EPSILON * r * t0;
        if !du.is_finite()
            || du >= 0.0
            || t1 >= t0
            || endpoint.state().state().pressure >= self.state().state().pressure
            || !exhausted.is_finite()
            || !mean_h.is_finite()
            || !mean_minus_h0.is_finite()
            || !energy_slack.is_finite()
            || mean_minus_h0 > energy_slack
            || mean_minus_h0 < dh - energy_slack
        {
            return Err(GasBlowdownError::Invalid(
                "remaining gas or integrated exhaust enthalpy is inconsistent with adiabatic depletion",
            ));
        }
        let incoming = inventory(self.model(), mass_kg)?;
        let before = inventory(self.model(), self.mass_kg())?;
        let after = inventory(endpoint.model(), endpoint.mass_kg())?;
        let mass_residual = (endpoint.mass_kg() - self.mass_kg()) + mass_kg;
        let mut maximum: f64 = 0.0;
        for (species, amount) in &incoming {
            let residual = (after[species] - before[species]) + amount;
            if !residual.is_finite() {
                return Err(GasBlowdownError::Invalid(
                    "exhaust species accounting overflow",
                ));
            }
            maximum = maximum.max(residual.abs());
        }
        if mass_residual.abs() > controls.mass_tolerance_kg || maximum > controls.mass_tolerance_kg
        {
            return Err(GasBlowdownError::Invalid(
                "exhaust did not close reservoir mass and species inventories",
            ));
        }
        poll(&mut cancelled)?;
        Ok(GasBlowdownCandidate {
            endpoint,
            step: GasBlowdownStep {
                expelled_mass_kg: mass_kg,
                exhausted_enthalpy_j: exhausted,
                mean_stagnation_enthalpy_j_kg: Some(mean_h),
                specific_entropy_residual_j_kg_k: entropy_residual,
                mass_residual_kg: mass_residual,
                maximum_species_mass_residual_kg: maximum,
                species: incoming
                    .into_iter()
                    .map(|(species, mass_kg)| GasSpeciesTransfer { species, mass_kg })
                    .collect(),
                iterations,
            },
        })
    }
}

/// Physical-step bounds for midpoint nozzle-rate integration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketBlowdownSolve {
    /// Converging-diverging nozzle roots at the same original gas composition.
    pub nozzle: RocketNozzleSolve,
    /// Reservoir entropy inversion and absolute mass quadrature tolerance.
    pub reservoir: GasBlowdownSolve,
    /// Largest mass fraction expelled in one physical step, strictly in (0,1).
    /// This is a step bound, not a truncation-error estimate.
    pub maximum_donor_fraction: f64,
    /// At most 128 midpoint mass-rate iterations, each bounded by the source solvers.
    pub maximum_iterations: usize,
}

/// A finite-reservoir rocket step with mass, enthalpy and thrust from one path.
#[derive(Clone, Debug, PartialEq)]
pub struct RocketBlowdownCandidate {
    reservoir: GasBlowdownCandidate,
    midpoint_flow: RocketNozzleFlow,
    /// Requested physical duration [s].
    pub duration_s: f64,
    /// Second-order midpoint quadrature of axial thrust [N s].
    /// The composing momentum solver must consume this impulse only once.
    pub axial_impulse_n_s: f64,
    /// Expelled mass minus duration times the midpoint nozzle rate [kg].
    pub mass_quadrature_residual_kg: f64,
    /// Number of midpoint rate iterations; source roots have separate budgets.
    pub iterations: usize,
}
impl RocketBlowdownCandidate {
    /// Exact thermodynamic depletion and its exported extensive quantities.
    #[must_use]
    pub const fn reservoir(&self) -> &GasBlowdownCandidate {
        &self.reservoir
    }
    /// Nozzle state at the midpoint of the remaining-mass interval.
    #[must_use]
    pub const fn midpoint_flow(&self) -> &RocketNozzleFlow {
        &self.midpoint_flow
    }
    /// Move the reservoir, its accounting and the midpoint nozzle state.
    #[must_use]
    pub fn into_parts(self) -> (ClosedGasParcel, GasBlowdownStep, RocketNozzleFlow) {
        let (endpoint, step) = self.reservoir.into_parts();
        (endpoint, step, self.midpoint_flow)
    }
}

impl ClosedGasParcel {
    /// Advance a finite adiabatic reservoir through an ideal rocket nozzle.
    /// Vacuum ambient is permitted. Initial, midpoint AND endpoint nozzle
    /// states must remain inside the nozzle's admitted expansion regime.
    ///
    /// Solve dm = dt*mdot(m0-dm/2) by bounded fixed-point iteration, reconstructing
    /// every trial on the original isentrope. This is second order in time for
    /// smooth admitted flows; the fixed mass-fraction cap does not certify a
    /// desired time error. Refine time steps for transient pressure and impulse.
    /// No mass clipping, source extrapolation or partial update is performed.
    pub fn rocket_blowdown(
        &self,
        ambient_pressure_pa: f64,
        nozzle: RocketNozzle,
        duration_s: f64,
        controls: RocketBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RocketBlowdownCandidate, GasBlowdownError> {
        controls.reservoir.validate()?;
        if !positive(duration_s)
            || !positive(controls.maximum_donor_fraction)
            || controls.maximum_donor_fraction >= 1.0
            || !(1..=128).contains(&controls.maximum_iterations)
        {
            return Err(GasBlowdownError::Invalid(
                "positive finite duration, donor fraction in (0,1) and 1..=128 midpoint iterations required",
            ));
        }
        let initial =
            self.rocket_nozzle_flow(ambient_pressure_pa, nozzle, controls.nozzle, &mut cancelled)?;
        let mut mass = duration_s * initial.mass_flow_kg_s;
        let mut last_residual = f64::INFINITY;
        for iterations in 1..=controls.maximum_iterations {
            if !positive(mass)
                || mass / self.mass_kg() > controls.maximum_donor_fraction
                || mass <= controls.reservoir.mass_tolerance_kg
            {
                return Err(GasBlowdownError::Invalid(
                    "rocket blowdown timestep exceeds the donor bound or mass resolution",
                ));
            }
            poll(&mut cancelled)?;
            let (midpoint, _, _) =
                remaining(self, 0.5 * mass, controls.reservoir.entropy, &mut cancelled)?;
            let flow = midpoint.rocket_nozzle_flow(
                ambient_pressure_pa,
                nozzle,
                controls.nozzle,
                &mut cancelled,
            )?;
            let target = duration_s * flow.mass_flow_kg_s;
            last_residual = mass - target;
            if !last_residual.is_finite() {
                return Err(GasBlowdownError::Invalid(
                    "rocket midpoint mass rate is unrepresentable",
                ));
            }
            if last_residual.abs() <= controls.reservoir.mass_tolerance_kg {
                let reservoir =
                    self.exhaust_mass_adiabatic(mass, controls.reservoir, &mut cancelled)?;
                // A valid midpoint does not permit crossing a source or nozzle
                // applicability boundary at the accepted end of the interval.
                reservoir.endpoint().rocket_nozzle_flow(
                    ambient_pressure_pa,
                    nozzle,
                    controls.nozzle,
                    &mut cancelled,
                )?;
                let impulse = duration_s * flow.thrust_n;
                if !positive(impulse) {
                    return Err(GasBlowdownError::Invalid(
                        "rocket thrust impulse is unrepresentable",
                    ));
                }
                poll(&mut cancelled)?;
                return Ok(RocketBlowdownCandidate {
                    reservoir,
                    midpoint_flow: flow,
                    duration_s,
                    axial_impulse_n_s: impulse,
                    mass_quadrature_residual_kg: last_residual,
                    iterations,
                });
            }
            mass = target;
        }
        Err(GasBlowdownError::NotConverged {
            iterations: controls.maximum_iterations,
            mass_residual_kg: last_residual,
        })
    }
}

#[cfg(test)]
mod tests;
