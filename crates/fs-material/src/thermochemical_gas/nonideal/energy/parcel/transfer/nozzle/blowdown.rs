//! Adiabatic depletion of one finite, well-mixed Peng--Robinson reservoir.
//!
//! At fixed tank volume, d(m u) = h dm and the source identity
//! du = T ds + p d(rho)/rho^2 imply ds = 0 for the remaining gas.
//! Reconstruct that original NASA/PR isentrope instead of holding pressure,
//! temperature or upstream enthalpy constant during a finite withdrawal.
//! Export U_initial - U_final once as exhaust stagnation enthalpy; pressure
//! work and jet kinetic energy are already included in this stream account.
//!
//! A bounded implicit-midpoint mass-rate solve supplies physical-time nozzle
//! blowdown. It is second order for smooth admitted flow, not a certified time
//! error bound. The original homogeneous stable branch must remain applicable;
//! there is no phase flash, wall heat, reaction, tank momentum or finite receiver.

use super::{GasNozzle, IsochoricRealGasParcel, RealGasNozzleFlow, RealGasNozzleStepSolve};
use crate::thermochemical_gas::energy::{
    entropy::GasEntropySolve, parcel::transfer::{GasSpeciesTransfer, inventory},
    poll, updated_point,
};
use crate::thermochemical_gas::positive;

/// The same source-entropy and extensive-inventory controls and reports as the
/// ideal-reference parcel. Their use here includes the original PR departures.
pub use crate::thermochemical_gas::energy::parcel::transfer::nozzle::blowdown::{
    GasBlowdownError as RealGasBlowdownError,
    GasBlowdownSolve as RealGasBlowdownSolve,
    GasBlowdownStep as RealGasBlowdownStep,
};

/// Fully reconstructed reservoir and exported species/enthalpy from one path.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasBlowdownCandidate {
    endpoint: IsochoricRealGasParcel,
    step: RealGasBlowdownStep,
}
impl RealGasBlowdownCandidate {
    #[must_use]
    pub const fn endpoint(&self) -> &IsochoricRealGasParcel { &self.endpoint }
    #[must_use]
    pub const fn step(&self) -> &RealGasBlowdownStep { &self.step }
    #[must_use]
    pub fn into_parts(self) -> (IsochoricRealGasParcel, RealGasBlowdownStep) {
        (self.endpoint, self.step)
    }
}

fn admit(
    parcel: &IsochoricRealGasParcel,
    controls: RealGasBlowdownSolve,
) -> Result<(), RealGasBlowdownError> {
    controls.entropy.validate()?;
    let temperature = parcel.state().temperature_k();
    if !positive(controls.mass_tolerance_kg)
        || !positive(parcel.mass_kg() * controls.entropy.entropy_tolerance_j_kg_k)
        || temperature < controls.entropy.minimum_temperature_k
        || temperature > controls.entropy.maximum_temperature_k
    {
        return Err(RealGasBlowdownError::Invalid(
            "positive real-gas inventory tolerances and an entropy bracket containing the reservoir required",
        ));
    }
    Ok(())
}

// Reuse one scratch parcel in the time solve. No source/composition clone or
// species-account allocation is needed for each midpoint rate iteration.
fn remaining(
    original: &IsochoricRealGasParcel,
    expelled: f64,
    controls: GasEntropySolve,
    target: &mut IsochoricRealGasParcel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, usize), RealGasBlowdownError> {
    let mass = original.mass_kg() - expelled;
    let density = mass / original.volume_m3();
    if !positive(mass) || !positive(density)
        || mass >= original.mass_kg() || density >= original.density_kg_m3()
    {
        return Err(RealGasBlowdownError::Invalid(
            "real-gas blowdown must retain a represented positive lower mass and density",
        ));
    }
    let solved = original.model().reconstruct_density_entropy_change(
        original.point(), original.density_kg_m3(), density, 0.0,
        controls, &mut *cancelled,
    )?;
    let residual = solved.residual_j_kg_k();
    let iterations = solved.iterations();
    let point = updated_point(
        original.point(), solved.state().temperature_k(), solved.state().pressure_pa(),
    )?;
    target.point = point;
    target.state = solved.into_state();
    target.mass_kg = mass;
    target.density_kg_m3 = density;
    Ok((residual, iterations))
}

impl IsochoricRealGasParcel {
    /// Withdraw finite mass adiabatically at the current fixed tank volume.
    ///
    /// Species fractions and all original source bindings remain unchanged.
    /// The complete remaining material follows ds=0 at its actual new density;
    /// finite exhaust carries the integral of changing stagnation enthalpy,
    /// not expelled mass times the initial enthalpy. The source's reference
    /// energy is retained in reports but cancels before the physical admission.
    ///
    /// This operation supplies no flow rate or back-pressure rule. Use
    /// `nozzle_blowdown` for a pressure-driven physical interval. Positive
    /// withdrawals must exceed the absolute mass tolerance. No vacuum endpoint,
    /// source extrapolation or partial mutation of the original is allowed.
    pub fn exhaust_mass_adiabatic(
        &self,
        mass_kg: f64,
        controls: RealGasBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasBlowdownCandidate, RealGasBlowdownError> {
        admit(self, controls)?;
        poll(&mut cancelled)?;
        if !mass_kg.is_finite() || mass_kg < 0.0 || mass_kg >= self.mass_kg() {
            return Err(RealGasBlowdownError::Invalid(
                "finite nonnegative real-gas exhaust smaller than reservoir inventory required",
            ));
        }
        if mass_kg == 0.0 {
            return Ok(RealGasBlowdownCandidate {
                endpoint: self.clone(),
                step: RealGasBlowdownStep {
                    expelled_mass_kg: 0.0,
                    exhausted_enthalpy_j: 0.0,
                    mean_stagnation_enthalpy_j_kg: None,
                    specific_entropy_residual_j_kg_k: 0.0,
                    mass_residual_kg: 0.0,
                    maximum_species_mass_residual_kg: 0.0,
                    species: Vec::new(),
                    iterations: 0,
                },
            });
        }
        if mass_kg <= controls.mass_tolerance_kg {
            return Err(RealGasBlowdownError::Invalid(
                "positive real-gas exhaust must exceed its absolute mass tolerance",
            ));
        }
        let mut endpoint = self.clone();
        let (entropy_residual, iterations) =
            remaining(self, mass_kg, controls.entropy, &mut endpoint, &mut cancelled)?;
        let t0 = self.state().temperature_k();
        let t1 = endpoint.state().temperature_k();
        let du = self.model().specific_internal_energy_change_between_densities(
            self.density_kg_m3(), t0, endpoint.density_kg_m3(), t1,
        )?;
        let dh = self.model().specific_enthalpy_change_between_densities(
            self.density_kg_m3(), t0, endpoint.density_kg_m3(), t1,
        )?;
        let u0 = self.state().specific_internal_energy_j_kg();
        let flow_work = self.state().pressure_pa() / self.density_kg_m3();
        let exhausted = mass_kg * u0 - endpoint.mass_kg() * du;
        let mean_increment = -(endpoint.mass_kg() / mass_kg) * du;
        let mean_h = u0 + mean_increment;
        // Along the admitted stable isentrope, dh = dp/rho and du = p drho/rho^2.
        // Thus h decreases during depletion and its mass-average lies between
        // the endpoint values. Compare relative increments, never large datums.
        let mean_minus_h0 = mean_increment - flow_work;
        let energy_slack = 4.0 * t0.max(t1)
            * controls.entropy.entropy_tolerance_j_kg_k * (self.mass_kg() / mass_kg)
            + 128.0 * f64::EPSILON * (flow_work + du.abs() + dh.abs());
        if !du.is_finite() || du >= 0.0 || !dh.is_finite() || dh >= 0.0
            || endpoint.state().pressure_pa() >= self.state().pressure_pa()
            || !exhausted.is_finite() || !mean_h.is_finite()
            || !mean_minus_h0.is_finite() || !energy_slack.is_finite()
            || mean_minus_h0 > energy_slack || mean_minus_h0 < dh - energy_slack
        {
            return Err(RealGasBlowdownError::Invalid(
                "real-gas depletion or integrated exhaust enthalpy violates the admitted isentrope",
            ));
        }
        let expelled = inventory(self.model().source(), mass_kg)?;
        let before = inventory(self.model().source(), self.mass_kg())?;
        let after = inventory(endpoint.model().source(), endpoint.mass_kg())?;
        let mass_residual = (endpoint.mass_kg() - self.mass_kg()) + mass_kg;
        let mut maximum_species: f64 = 0.0;
        for (species, amount) in &expelled {
            let residual = (after[species] - before[species]) + amount;
            if !residual.is_finite() {
                return Err(RealGasBlowdownError::Invalid(
                    "real-gas exhaust species accounting is unrepresentable",
                ));
            }
            maximum_species = maximum_species.max(residual.abs());
        }
        if mass_residual.abs() > controls.mass_tolerance_kg
            || maximum_species > controls.mass_tolerance_kg
        {
            return Err(RealGasBlowdownError::Invalid(
                "real-gas exhaust did not close total and species inventories",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RealGasBlowdownCandidate {
            endpoint,
            step: RealGasBlowdownStep {
                expelled_mass_kg: mass_kg,
                exhausted_enthalpy_j: exhausted,
                mean_stagnation_enthalpy_j_kg: Some(mean_h),
                specific_entropy_residual_j_kg_k: entropy_residual,
                mass_residual_kg: mass_residual,
                maximum_species_mass_residual_kg: maximum_species,
                species: expelled.into_iter()
                    .map(|(species, mass_kg)| GasSpeciesTransfer { species, mass_kg }).collect(),
                iterations,
            },
        })
    }
}

/// Source and physical-work bounds for a finite nozzle blowdown interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasNozzleBlowdownSolve {
    /// Original PR nozzle closure. Its maximum donor fraction bounds the total
    /// withdrawal; the paired-transfer energy controls are validated but unused.
    pub nozzle: RealGasNozzleStepSolve,
    /// Reservoir entropy reconstruction and absolute mass quadrature tolerance.
    pub reservoir: RealGasBlowdownSolve,
    /// At most 128 implicit-midpoint mass-rate iterations.
    pub maximum_iterations: usize,
}

/// One complete pressure-driven physical interval, retaining the original source.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleBlowdownCandidate {
    reservoir: RealGasBlowdownCandidate,
    midpoint_flow: RealGasNozzleFlow,
    /// Requested physical duration [s].
    pub duration_s: f64,
    /// Expelled mass minus duration times the midpoint nozzle rate [kg].
    pub mass_quadrature_residual_kg: f64,
    /// Midpoint rate iterations; inner source roots have separate budgets.
    pub iterations: usize,
}
impl RealGasNozzleBlowdownCandidate {
    #[must_use]
    pub const fn reservoir(&self) -> &RealGasBlowdownCandidate { &self.reservoir }
    #[must_use]
    pub const fn midpoint_flow(&self) -> &RealGasNozzleFlow { &self.midpoint_flow }
    #[must_use]
    pub fn into_parts(self) -> (IsochoricRealGasParcel, RealGasBlowdownStep, RealGasNozzleFlow) {
        let (endpoint, step) = self.reservoir.into_parts();
        (endpoint, step, self.midpoint_flow)
    }
}

impl IsochoricRealGasParcel {
    /// Advance a finite adiabatic tank through its converging real-gas nozzle.
    ///
    /// Solve dm = dt*mdot(m0-dm/2) on the original reservoir isentrope. Pressure,
    /// temperature, choking and mass rate therefore evolve together. Initial,
    /// midpoint and final nozzle states must all remain source-admitted. Vacuum
    /// is a valid external back pressure, never a zero-pressure EOS query.
    ///
    /// The fixed-point work and donor-fraction bounds are explicit. They do not
    /// certify a time error: refine physical intervals for transient pressure.
    /// An interval crossing back-pressure equilibrium or a source/phase domain
    /// refuses and must be shortened. It is never clipped or silently stopped
    /// early. Equal pressures or a closed nozzle return exact identity.
    pub fn nozzle_blowdown(
        &self,
        back_pressure_pa: f64,
        nozzle: GasNozzle,
        duration_s: f64,
        controls: RealGasNozzleBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasNozzleBlowdownCandidate, RealGasBlowdownError> {
        admit(self, controls.reservoir)?;
        controls.nozzle.validate()?;
        if !positive(duration_s) || !(1..=128).contains(&controls.maximum_iterations) {
            return Err(RealGasBlowdownError::Invalid(
                "positive real-gas blowdown duration and 1..=128 midpoint iterations required",
            ));
        }
        let initial = self.nozzle_flow_to(
            back_pressure_pa, nozzle, controls.nozzle, &mut cancelled,
        )?;
        if initial.mass_flow_kg_s == 0.0 {
            let reservoir = self.exhaust_mass_adiabatic(
                0.0, controls.reservoir, &mut cancelled,
            )?;
            poll(&mut cancelled)?;
            return Ok(RealGasNozzleBlowdownCandidate {
                reservoir, midpoint_flow: initial, duration_s,
                mass_quadrature_residual_kg: 0.0, iterations: 0,
            });
        }
        let mut mass = duration_s * initial.mass_flow_kg_s;
        let mut midpoint = self.clone();
        let mut last_residual = f64::INFINITY;
        for iterations in 1..=controls.maximum_iterations {
            if !positive(mass) || mass <= controls.reservoir.mass_tolerance_kg
                || mass / self.mass_kg() > controls.nozzle.maximum_donor_fraction
            {
                return Err(RealGasBlowdownError::Invalid(
                    "real-gas blowdown interval exceeds donor bounds or mass resolution",
                ));
            }
            poll(&mut cancelled)?;
            remaining(self, 0.5 * mass, controls.reservoir.entropy,
                &mut midpoint, &mut cancelled)?;
            let flow = midpoint.nozzle_flow_to(
                back_pressure_pa, nozzle, controls.nozzle, &mut cancelled,
            )?;
            let target = duration_s * flow.mass_flow_kg_s;
            last_residual = mass - target;
            if !last_residual.is_finite() {
                return Err(RealGasBlowdownError::Invalid(
                    "real-gas midpoint mass quadrature is unrepresentable",
                ));
            }
            if last_residual.abs() <= controls.reservoir.mass_tolerance_kg {
                let reservoir = self.exhaust_mass_adiabatic(
                    mass, controls.reservoir, &mut cancelled,
                )?;
                if reservoir.endpoint().state().pressure_pa() < back_pressure_pa {
                    return Err(RealGasBlowdownError::Invalid(
                        "real-gas blowdown crosses back-pressure equilibrium; reduce the interval",
                    ));
                }
                reservoir.endpoint().nozzle_flow_to(
                    back_pressure_pa, nozzle, controls.nozzle, &mut cancelled,
                )?;
                poll(&mut cancelled)?;
                return Ok(RealGasNozzleBlowdownCandidate {
                    reservoir, midpoint_flow: flow, duration_s,
                    mass_quadrature_residual_kg: last_residual, iterations,
                });
            }
            mass = target;
        }
        Err(RealGasBlowdownError::NotConverged {
            iterations: controls.maximum_iterations,
            mass_residual_kg: last_residual,
        })
    }
}

#[cfg(test)]
mod tests;
