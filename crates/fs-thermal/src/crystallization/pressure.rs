//! Fixed-mass crystallization calorimetry at a declared constant absolute pressure.
//!
//! Original phase densities give `v=(1-X)/rho_parent+X/rho_crystal` at the actual
//! kinetic MASS fraction. Heat supplies the existing enthalpy kernel unchanged.
//! This module accounts for `delta_U=delta_H-p*delta_V`; it never adds pressure
//! work or latent heat as another enthalpy source. Here `u=h-pv` DEFINES the
//! internal-energy convention, rather than independently verifying an EOS.
//!
//! Pressure must be the fixed condition at which both original calorics and
//! densities apply. A caller-supplied number does not change their source
//! pressure, predict a pressure-dependent melting point, or validate source
//! applicability. Only quasistatic pressure-volume work is included; kinetic,
//! deviatoric/elastic, surface, chemical and other work need their own owners.
//! No moving mesh, fluid continuity solve or constrained-volume response is added.

/// Continuing accepted phase inventory and actual volume under this pressure load.
#[path = "pressure/parcel.rs"]
pub mod parcel;

/// Implicit conservative heat exchange between two original phase inventories.
#[path = "pressure/exchange.rs"]
pub mod exchange;

use super::{
    AdaptiveCrystallizationStep, CrystallizationAdaptivity, CrystallizationError,
    CrystallizationLaw, CrystallizationSolve, CrystallizationState, EquivalentAgeRemelting,
    check_cancelled, nonnegative, positive,
};

/// Absolute arithmetic/admission budgets, separate from material and time error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricCrystallizationControls {
    pub caloric: CrystallizationSolve,
    /// Noncancelling heat, pressure-work, split-arithmetic and caloric defects, J.
    pub energy_tolerance_j: f64,
    /// Difference between published volume change and stable source change, m3.
    /// Enforced even at zero pressure; a zero work term cannot hide lost volume.
    pub volume_tolerance_m3: f64,
}

impl Default for IsobaricCrystallizationControls {
    fn default() -> Self {
        Self {
            caloric: CrystallizationSolve::default(),
            energy_tolerance_j: 1e-7,
            volume_tolerance_m3: 1e-15,
        }
    }
}

impl IsobaricCrystallizationControls {
    pub fn validate(self) -> Result<(), CrystallizationError> {
        self.caloric.validate()?;
        if !nonnegative(self.energy_tolerance_j) || !nonnegative(self.volume_tolerance_m3) {
            return Err(invalid(
                "isobaric tolerances must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

/// One complete original-state proposal and its fixed-mass pressure-work account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricCrystallizationStep {
    pub initial: CrystallizationState,
    pub state: CrystallizationState,
    /// Zero only for an instantaneous heat transaction.
    pub duration_s: f64,
    pub mass_kg: f64,
    pub pressure_pa: f64,
    pub supplied_heat_j: f64,
    pub initial_volume_m3: f64,
    pub final_volume_m3: f64,
    /// ACTUAL subtraction of the two published volumes, used for boundary work.
    pub volume_change_m3: f64,
    /// Stable original-density difference, including thermal and phase changes.
    pub source_volume_change_m3: f64,
    /// Represented `m*(h_final-h_initial)`, with no absolute-datum subtraction.
    pub enthalpy_change_j: f64,
    /// `delta_H - p*source_delta_V`; no independent internal-energy EOS is inferred.
    pub internal_energy_change_j: f64,
    /// Work OUT of the parcel, `p*(V_final-V_initial)`, using published volumes.
    pub boundary_work_j: f64,
    /// Signed `delta_U + boundary_work - supplied_heat`, summed with compensation.
    pub energy_residual_j: f64,
    /// Sum of absolute heat/volume/work-split and endpoint caloric defects.
    pub absolute_energy_residual_j: f64,
    pub caloric_residual_j_kg: f64,
    /// The original adaptive result, unchanged. Absent for instantaneous heat.
    pub adaptive: Option<AdaptiveCrystallizationStep>,
}

/// Prepared access to BOTH identified original caloric/density branches.
#[derive(Clone, Copy)]
pub struct IsobaricCrystallization<'law, 'source> {
    law: &'law CrystallizationLaw<'source>,
    pressure_pa: f64,
}

impl<'law, 'source> IsobaricCrystallization<'law, 'source> {
    /// The caller admits the common fixed source pressure. Endpoint checks
    /// reject missing/invalid density data; every actual state is queried too.
    pub fn new(
        law: &'law CrystallizationLaw<'source>,
        pressure_pa: f64,
    ) -> Result<Self, CrystallizationError> {
        if !nonnegative(pressure_pa) {
            return Err(invalid(
                "isobaric absolute pressure must be finite and nonnegative",
            ));
        }
        if law.parent_caloric_identity().is_none() || law.crystalline_caloric_identity().is_none() {
            return Err(invalid(
                "isobaric crystallization requires both original source identities",
            ));
        }
        let prepared = Self { law, pressure_pa };
        for temperature in law.config().temperature_range_k {
            prepared.densities(temperature, &mut || false)?;
        }
        Ok(prepared)
    }

    pub const fn law(&self) -> &'law CrystallizationLaw<'source> {
        self.law
    }

    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }

    /// Free mixture volume from original densities at the actual accepted X.
    pub fn volume_m3(
        &self,
        current: CrystallizationState,
        mass_kg: f64,
        solve: CrystallizationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, CrystallizationError> {
        check_cancelled(&mut cancelled)?;
        check_mass(mass_kg)?;
        self.law.validate_state(current, solve)?;
        let density = self.densities(current.temperature_k(), &mut cancelled)?;
        let volume = density.volume(current.crystalline_mass_fraction(), mass_kg)?;
        check_cancelled(&mut cancelled)?;
        Ok(volume)
    }

    /// Zero-time heat transfer on the original caloric/remelting owner.
    pub fn apply_heat(
        &self,
        current: CrystallizationState,
        mass_kg: f64,
        heat_j: f64,
        controls: IsobaricCrystallizationControls,
        remelting: Option<EquivalentAgeRemelting>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricCrystallizationStep, CrystallizationError> {
        let input = self.heat_input(current, mass_kg, heat_j, controls, &mut cancelled)?;
        let step = self.law.apply_specific_enthalpy(
            current,
            input.target_h,
            controls.caloric,
            remelting,
            &mut cancelled,
        )?;
        self.account(
            input,
            step.state,
            0.0,
            step.caloric_residual_j_kg,
            None,
            controls,
            &mut cancelled,
        )
    }

    /// Finite-time heat transfer delegates to the SAME adaptive phase kernel.
    /// A zero supplied heat still permits adiabatic crystallization/contraction.
    /// Pressure work is accounted afterwards, never fed back as another heater.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_heat(
        &self,
        current: CrystallizationState,
        mass_kg: f64,
        heat_j: f64,
        duration_s: f64,
        controls: IsobaricCrystallizationControls,
        adaptivity: CrystallizationAdaptivity,
        remelting: Option<EquivalentAgeRemelting>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricCrystallizationStep, CrystallizationError> {
        let input = self.heat_input(current, mass_kg, heat_j, controls, &mut cancelled)?;
        let step = self.law.advance_specific_enthalpy_adaptive(
            current,
            input.target_h,
            duration_s,
            controls.caloric,
            adaptivity,
            remelting,
            &mut cancelled,
        )?;
        self.account(
            input,
            step.state,
            duration_s,
            step.caloric_residual_j_kg,
            Some(step),
            controls,
            &mut cancelled,
        )
    }

    fn densities(
        &self,
        temperature: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Densities, CrystallizationError> {
        check_cancelled(cancelled)?;
        self.law.check_temperature(temperature)?;
        let parent = self
            .law
            .amorphous
            .density_kg_m3(temperature)?
            .ok_or(invalid("original parent-phase density unavailable"))?;
        check_cancelled(cancelled)?;
        let crystal = self
            .law
            .crystalline
            .density_kg_m3(temperature)?
            .ok_or(invalid("original crystalline-phase density unavailable"))?;
        if [parent, crystal]
            .into_iter()
            .any(|rho| !positive(rho) || !positive(rho.recip()))
        {
            return Err(invalid(
                "original phase density and reciprocal volume must be finite and positive",
            ));
        }
        check_cancelled(cancelled)?;
        Ok(Densities { parent, crystal })
    }

    fn heat_input(
        &self,
        current: CrystallizationState,
        mass_kg: f64,
        heat_j: f64,
        controls: IsobaricCrystallizationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeatInput, CrystallizationError> {
        check_cancelled(cancelled)?;
        controls.validate()?;
        check_mass(mass_kg)?;
        self.law.validate_state(current, controls.caloric)?;
        let dh = finite(heat_j / mass_kg)?;
        if heat_j != 0.0 && dh == 0.0 {
            return Err(invalid("isobaric specific heat increment underflows"));
        }
        let target_h = if heat_j == 0.0 {
            current.specific_enthalpy_j_kg()
        } else {
            finite(current.specific_enthalpy_j_kg() + dh)?
        };
        let density = self.densities(current.temperature_k(), cancelled)?;
        let volume = density.volume(current.crystalline_mass_fraction(), mass_kg)?;
        let caloric_residual = finite(
            self.law.specific_enthalpy_j_kg(
                current.temperature_k(),
                current.crystalline_mass_fraction(),
            )? - current.specific_enthalpy_j_kg(),
        )?;
        Ok(HeatInput {
            current,
            mass_kg,
            heat_j,
            target_h,
            density,
            volume,
            caloric_residual,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn account(
        &self,
        input: HeatInput,
        state: CrystallizationState,
        duration_s: f64,
        caloric_residual_j_kg: f64,
        adaptive: Option<AdaptiveCrystallizationStep>,
        controls: IsobaricCrystallizationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IsobaricCrystallizationStep, CrystallizationError> {
        check_cancelled(cancelled)?;
        let final_density = if state.temperature_k() == input.current.temperature_k() {
            input.density
        } else {
            self.densities(state.temperature_k(), cancelled)?
        };
        let final_volume_m3 =
            final_density.volume(state.crystalline_mass_fraction(), input.mass_kg)?;
        let volume_change_m3 = finite(final_volume_m3 - input.volume)?;
        let source_volume_change_m3 = product(
            input.mass_kg,
            density_volume_change(
                input.density,
                final_density,
                input.current.crystalline_mass_fraction(),
                state.crystalline_mass_fraction(),
            )?,
        )?;
        let volume_defect = finite(volume_change_m3 - source_volume_change_m3)?;
        if volume_defect.abs() > controls.volume_tolerance_m3 {
            return Err(invalid(
                "isobaric volume representation exceeds its tolerance",
            ));
        }
        let enthalpy_change_j = product(
            input.mass_kg,
            finite(state.specific_enthalpy_j_kg() - input.current.specific_enthalpy_j_kg())?,
        )?;
        let source_work = product(self.pressure_pa, source_volume_change_m3)?;
        let boundary_work_j = product(self.pressure_pa, volume_change_m3)?;
        let pressure_defect = product(self.pressure_pa, volume_defect)?;
        let internal_energy_change_j = finite(enthalpy_change_j - source_work)?;
        let heat_defect = finite(enthalpy_change_j - input.heat_j)?;
        let split_defect = sum(&[internal_energy_change_j, source_work, -enthalpy_change_j])?;
        let work_defect = sum(&[boundary_work_j, -source_work, -pressure_defect])?;
        let energy_residual_j = sum(&[internal_energy_change_j, boundary_work_j, -input.heat_j])?;
        let absolute_energy_residual_j = sum(&[
            heat_defect.abs(),
            pressure_defect.abs(),
            split_defect.abs(),
            work_defect.abs(),
            product(input.mass_kg, input.caloric_residual.abs())?,
            product(input.mass_kg, caloric_residual_j_kg.abs())?,
        ])?
        .max(energy_residual_j.abs());
        if absolute_energy_residual_j > controls.energy_tolerance_j {
            return Err(invalid("isobaric energy residual exceeds its tolerance"));
        }
        check_cancelled(cancelled)?;
        Ok(IsobaricCrystallizationStep {
            initial: input.current,
            state,
            duration_s,
            mass_kg: input.mass_kg,
            pressure_pa: self.pressure_pa,
            supplied_heat_j: input.heat_j,
            initial_volume_m3: input.volume,
            final_volume_m3,
            volume_change_m3,
            source_volume_change_m3,
            enthalpy_change_j,
            internal_energy_change_j,
            boundary_work_j,
            energy_residual_j,
            absolute_energy_residual_j,
            caloric_residual_j_kg,
            adaptive,
        })
    }
}

#[derive(Clone, Copy)]
struct Densities {
    parent: f64,
    crystal: f64,
}

impl Densities {
    fn volume(self, fraction: f64, mass: f64) -> Result<f64, CrystallizationError> {
        let vp = product(1.0 - fraction, self.parent.recip())?;
        let vc = product(fraction, self.crystal.recip())?;
        let volume = product(mass, sum(&[vp, vc])?)?;
        if !positive(volume) {
            return Err(invalid(
                "isobaric mixture volume must be finite and positive",
            ));
        }
        Ok(volume)
    }
}

struct HeatInput {
    current: CrystallizationState,
    mass_kg: f64,
    heat_j: f64,
    target_h: f64,
    density: Densities,
    volume: f64,
    caloric_residual: f64,
}

fn density_volume_change(
    initial: Densities,
    final_density: Densities,
    initial_x: f64,
    final_x: f64,
) -> Result<f64, CrystallizationError> {
    // Exact symmetric product difference for v=vp+X*(vc-vp). This is an
    // endpoint decomposition, not separate physical heat or work mechanisms.
    let mean_x = if initial_x == final_x {
        initial_x
    } else {
        sum(&[product(0.5, initial_x)?, product(0.5, final_x)?])?
    };
    // Averaging the complementary fractions separately retains a very small
    // parent fraction even when the rounded mean crystal fraction equals one.
    let mean_parent = if initial_x == final_x {
        1.0 - initial_x
    } else {
        sum(&[product(0.5, 1.0 - initial_x)?, product(0.5, 1.0 - final_x)?])?
    };
    let dx = final_x - initial_x;
    let parent_change = if mean_parent == 0.0 {
        0.0
    } else {
        reciprocal_change(initial.parent, final_density.parent)?
    };
    let crystal_change = if mean_x == 0.0 {
        0.0
    } else {
        reciprocal_change(initial.crystal, final_density.crystal)?
    };
    let phase_change = if dx == 0.0 {
        0.0
    } else {
        let old_gap = reciprocal_change(initial.parent, initial.crystal)?;
        let new_gap = reciprocal_change(final_density.parent, final_density.crystal)?;
        let mean_gap = if old_gap == new_gap {
            old_gap
        } else {
            sum(&[product(0.5, old_gap)?, product(0.5, new_gap)?])?
        };
        product(dx, mean_gap)?
    };
    sum(&[
        product(mean_parent, parent_change)?,
        product(mean_x, crystal_change)?,
        phase_change,
    ])
}

/// Difference of reciprocal densities, without first rounding either volume.
fn reciprocal_change(initial: f64, final_density: f64) -> Result<f64, CrystallizationError> {
    let difference = initial - final_density;
    if difference == 0.0 {
        return Ok(0.0);
    }
    let change = finite((difference / initial.max(final_density)) / initial.min(final_density))?;
    if change == 0.0 {
        return Err(invalid("original density volume difference underflows"));
    }
    Ok(change)
}

fn check_mass(mass: f64) -> Result<(), CrystallizationError> {
    if positive(mass) {
        Ok(())
    } else {
        Err(invalid("isobaric mass must be finite and positive"))
    }
}

fn product(a: f64, b: f64) -> Result<f64, CrystallizationError> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(invalid("required isobaric product underflows"));
    }
    Ok(value)
}

/// Fixed-order Neumaier sum retains a small residual after two large accounts cancel.
fn sum(values: &[f64]) -> Result<f64, CrystallizationError> {
    let mut total: f64 = 0.0;
    let mut correction = 0.0;
    for &value in values {
        let next = finite(total + value)?;
        let remainder = if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        correction = finite(correction + remainder)?;
        total = next;
    }
    finite(total + correction)
}

fn finite(value: f64) -> Result<f64, CrystallizationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite isobaric crystallization arithmetic"))
    }
}

fn invalid(what: &'static str) -> CrystallizationError {
    CrystallizationError::Invalid(what)
}

#[cfg(test)]
#[path = "pressure/tests.rs"]
mod tests;
