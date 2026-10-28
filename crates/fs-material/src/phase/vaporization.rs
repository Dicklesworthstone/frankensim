//! Saturated liquid/vapor water at one constant absolute pressure.
//!
//! The original IAPWS saturation source owns temperature, both phase densities,
//! and caloric endpoints. Total specific enthalpy is the independent coordinate;
//! vapor quality is a mass fraction, never a solid/liquid phase label. Specific
//! volumes add. This is an equilibrium, spatially uniform constitutive law, not
//! a nucleation, evaporation-rate, pressure-control or moving-interface solver.
//!
//! Only the closed coexistence interval is admitted: neither subcooled liquid
//! nor superheated vapor is extrapolated. At the critical point the two phases
//! coincide and quality is undefined. A thermal consumer using this law models
//! a quasistatic external pressure load; it must not interpret its enthalpy as
//! the internal energy of a rigid, sealed volume. Heat-exchange geometry remains
//! the transport owner's explicit reference geometry.

/// Steady adiabatic pressure-drop flashing on the same original water source.
pub mod flash;

use fs_qty::Pressure;
use fs_thermochem::water::{SaturatedWaterState, WaterSaturationError, WaterSaturationSolve};

/// Immutable source-backed coexistence chart at one pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricWaterVaporization {
    saturation: SaturatedWaterState,
    enthalpy_interval_j_kg: f64,
}

/// A validated liquid/vapor state on one immutable saturation chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterVaporizationState {
    law: IsobaricWaterVaporization,
    specific_enthalpy_j_kg: f64,
    vapor_mass_fraction: f64,
}

impl WaterVaporizationState {
    /// Retained original pressure chart for subsequent heating or condensation.
    /// Reusing it avoids another source evaluation or pressure inversion.
    pub const fn law(self) -> IsobaricWaterVaporization {
        self.law
    }

    /// Absolute equilibrium temperature [K].
    pub fn temperature_k(self) -> f64 {
        self.law.saturation.temperature().value()
    }
    /// Constant absolute pressure [Pa].
    pub fn pressure_pa(self) -> f64 {
        self.law.saturation.pressure().value()
    }
    /// Total specific enthalpy, including vaporization [J/kg].
    pub const fn specific_enthalpy_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    /// Actual vapor mass fraction, including exactly pure endpoints.
    pub const fn vapor_mass_fraction(self) -> f64 {
        self.vapor_mass_fraction
    }
    /// Actual liquid mass fraction; no solid inventory is implied.
    pub fn liquid_mass_fraction(self) -> f64 {
        1.0 - self.vapor_mass_fraction
    }
    /// Additive phase specific volume [m3/kg], not an averaged density.
    pub fn specific_volume_m3_kg(self) -> f64 {
        let source = self.law.saturation;
        if self.vapor_mass_fraction == 1.0 {
            1.0 / source.vapor_density().value()
        } else {
            1.0 / source.liquid_density().value()
                + self.vapor_mass_fraction * source.specific_volume_gap().value()
        }
    }
    /// Bulk density from the actual liquid/vapor inventory [kg/m3].
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    /// Original-source mixture internal energy [J/kg]. The source datum is retained.
    pub fn specific_internal_energy_j_kg(self) -> f64 {
        let source = self.law.saturation;
        if self.vapor_mass_fraction == 1.0 {
            source.vapor_specific_internal_energy().value()
        } else {
            source.liquid_specific_internal_energy().value()
                + self.vapor_mass_fraction * source.latent_internal_energy().value()
        }
    }
    /// Original-source mixture entropy [J/(kg K)], including latent entropy.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        let source = self.law.saturation;
        if self.vapor_mass_fraction == 1.0 {
            source.vapor_specific_entropy().value()
        } else {
            source.liquid_specific_entropy().value()
                + self.vapor_mass_fraction * source.latent_entropy().value()
        }
    }
}

/// Complete immutable heat proposal; a failure changes neither input state nor law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterVaporizationHeatStep {
    /// Original accepted state.
    pub initial: WaterVaporizationState,
    /// Proposed state at the same pressure and fixed total mass.
    pub state: WaterVaporizationState,
    /// Invariant total mass [kg].
    pub mass_kg: f64,
    /// Requested heat into the specimen [J]. Cooling is negative.
    pub supplied_heat_j: f64,
    /// Original [liquid, vapor] phase inventories [kg], separately rounded.
    pub initial_phase_masses_kg: [f64; 2],
    /// Proposed [liquid, vapor] phase inventories [kg], separately rounded.
    /// Their displayed sum can differ from the invariant total mass by roundoff.
    pub phase_masses_kg: [f64; 2],
    /// Original and proposed total volumes [m3], separately rounded.
    pub volumes_m3: [f64; 2],
    /// Volume change from the actual represented vapor-mass change and source
    /// specific-volume gap [m3]. Subtracting displayed volumes can lose this change.
    pub volume_change_m3: f64,
    /// Actual m*(h_final-h_initial) [J].
    pub enthalpy_change_j: f64,
    /// Work delivered to the pressure load, p*delta_V [J].
    pub boundary_work_j: f64,
    /// Source latent internal energy times actual vapor-mass change [J].
    /// Evaluated independently of supplied heat and the enthalpy increment.
    pub internal_energy_change_j: f64,
    /// Source latent entropy times actual vapor-mass change [J/K].
    pub entropy_change_j_k: f64,
    /// Requested heat minus represented enthalpy change [J].
    pub heat_residual_j: f64,
    /// Requested heat minus source internal-energy change and outward work [J].
    /// Includes source endpoint/gap rounding and represented phase-mass changes.
    pub first_law_residual_j: f64,
}

/// Explicit refusal of source, domain, representation or caller accuracy limits.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterVaporizationError {
    /// Original saturation source refused the pressure query.
    Source(WaterSaturationError),
    /// Coexistence has no positive representable latent/volume interval.
    CollapsedCoexistence,
    /// No subcooling or superheating branch is provided.
    OutsideCoexistence,
    /// The proposed state belongs to a different source pressure/chart.
    StateMismatch,
    /// Invalid finite mass, heat or absolute accuracy budget.
    InvalidInput(&'static str),
    /// A required physical increment or endpoint is not representable.
    Unrepresentable(&'static str),
    /// At least one actual physical residual exceeded the caller's J budget.
    EnergyBalance {
        heat_residual_j: f64,
        first_law_residual_j: f64,
    },
    /// Caller cancellation; no proposal is published.
    Cancelled,
}

impl std::fmt::Display for WaterVaporizationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "isobaric water vaporization refusal: {self:?}")
    }
}
impl std::error::Error for WaterVaporizationError {}

impl IsobaricWaterVaporization {
    /// Bind one validated original saturation state, excluding critical collapse.
    /// The represented source h endpoints define the exact admitted interval.
    /// Their difference can differ by roundoff from the source's directly
    /// evaluated latent heat; the heat-step first-law gate retains that difference.
    pub fn from_saturation(
        saturation: SaturatedWaterState,
    ) -> Result<Self, WaterVaporizationError> {
        let interval = saturation.vapor_specific_enthalpy().value()
            - saturation.liquid_specific_enthalpy().value();
        if !interval.is_finite()
            || interval <= 0.0
            || saturation.latent_enthalpy().value() <= 0.0
            || saturation.specific_volume_gap().value() <= 0.0
        {
            return Err(WaterVaporizationError::CollapsedCoexistence);
        }
        Ok(Self {
            saturation,
            enthalpy_interval_j_kg: interval,
        })
    }

    /// Use the original bounded saturation-pressure inversion; no second root solver.
    /// The chart retains the evaluated source pressure within the caller's Pa
    /// residual budget. It does not replace that pressure with the requested value.
    pub fn at_pressure(
        pressure: Pressure,
        controls: WaterSaturationSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WaterVaporizationError> {
        let source = SaturatedWaterState::at_pressure(pressure, controls, cancelled)
            .map_err(WaterVaporizationError::Source)?;
        Self::from_saturation(source.into_state())
    }

    /// Original source, including directly evaluated latent gaps.
    pub const fn saturation(self) -> SaturatedWaterState {
        self.saturation
    }

    /// Resolve total specific enthalpy to real liquid/vapor inventory without allocation.
    pub fn state_at_specific_enthalpy(
        self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<WaterVaporizationState, WaterVaporizationError> {
        let liquid = self.saturation.liquid_specific_enthalpy().value();
        let vapor = self.saturation.vapor_specific_enthalpy().value();
        if !specific_enthalpy_j_kg.is_finite()
            || specific_enthalpy_j_kg < liquid
            || specific_enthalpy_j_kg > vapor
        {
            return Err(WaterVaporizationError::OutsideCoexistence);
        }
        let fraction = if specific_enthalpy_j_kg == vapor {
            1.0
        } else {
            (specific_enthalpy_j_kg - liquid) / self.enthalpy_interval_j_kg
        };
        if specific_enthalpy_j_kg > liquid
            && specific_enthalpy_j_kg < vapor
            && (fraction == 0.0 || fraction == 1.0)
        {
            return Err(WaterVaporizationError::Unrepresentable("vapor quality"));
        }
        let state = WaterVaporizationState {
            law: self,
            specific_enthalpy_j_kg,
            vapor_mass_fraction: fraction,
        };
        if !state.bulk_density_kg_m3().is_finite()
            || state.bulk_density_kg_m3() <= 0.0
            || !state.specific_internal_energy_j_kg().is_finite()
            || !state.specific_entropy_j_kg_k().is_finite()
        {
            return Err(WaterVaporizationError::Unrepresentable("mixture state"));
        }
        Ok(state)
    }

    /// Propose reversible equilibrium boiling/condensation of a fixed mass.
    ///
    /// Heat changes H once. The actual represented vapor inventory supplies
    /// independent delta_U and p*delta_V from the original source gaps. Both
    /// Q-delta_H and Q-delta_U-p*delta_V must fit the explicit absolute J budget;
    /// zero requests strict representational closure. No source uncertainty or
    /// time-integration error is claimed. A completely lost nonzero heat/phase
    /// increment always refuses. No accepted state is mutated, even on cancellation.
    pub fn advance_heat(
        self,
        current: WaterVaporizationState,
        mass_kg: f64,
        supplied_heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterVaporizationHeatStep, WaterVaporizationError> {
        use WaterVaporizationError as Error;
        if cancelled() {
            return Err(Error::Cancelled);
        }
        if current.law != self {
            return Err(Error::StateMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !supplied_heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(Error::InvalidInput(
                "positive mass, finite heat and nonnegative J budget required",
            ));
        }
        let dh_requested = supplied_heat_j / mass_kg;
        if !dh_requested.is_finite() || (supplied_heat_j != 0.0 && dh_requested == 0.0) {
            return Err(Error::Unrepresentable("specific heat increment"));
        }
        let liquid_h = self.saturation.liquid_specific_enthalpy().value();
        let vapor_h = self.saturation.vapor_specific_enthalpy().value();
        let liquid_heat = mass_kg * (liquid_h - current.specific_enthalpy_j_kg);
        let vapor_heat = mass_kg * (vapor_h - current.specific_enthalpy_j_kg);
        if supplied_heat_j < liquid_heat || supplied_heat_j > vapor_heat {
            return Err(Error::OutsideCoexistence);
        }
        // Recognize an exact request for a source endpoint before Q/m and h+dh
        // round it outside coexistence or leave a fictitious residual phase.
        // This does not clip excess heat; both actual caloric gates still apply.
        let target_h = if supplied_heat_j != 0.0 && supplied_heat_j == liquid_heat {
            liquid_h
        } else if supplied_heat_j != 0.0 && supplied_heat_j == vapor_heat {
            vapor_h
        } else {
            current.specific_enthalpy_j_kg + dh_requested
        };
        let state = self.state_at_specific_enthalpy(target_h)?;
        let enthalpy_change_j =
            mass_kg * (state.specific_enthalpy_j_kg - current.specific_enthalpy_j_kg);
        let masses = |state: WaterVaporizationState| {
            [
                mass_kg * state.liquid_mass_fraction(),
                mass_kg * state.vapor_mass_fraction(),
            ]
        };
        let initial_phase_masses_kg = masses(current);
        let phase_masses_kg = masses(state);
        let vapor_mass_change = phase_masses_kg[1] - initial_phase_masses_kg[1];
        let volumes_m3 = [
            mass_kg * current.specific_volume_m3_kg(),
            mass_kg * state.specific_volume_m3_kg(),
        ];
        let volume_change_m3 = vapor_mass_change * self.saturation.specific_volume_gap().value();
        let boundary_work_j = self.saturation.pressure().value() * volume_change_m3;
        let internal_energy_change_j =
            vapor_mass_change * self.saturation.latent_internal_energy().value();
        let entropy_change_j_k = vapor_mass_change * self.saturation.latent_entropy().value();
        let heat_residual_j = supplied_heat_j - enthalpy_change_j;
        let first_law_residual_j = (supplied_heat_j - internal_energy_change_j) - boundary_work_j;
        if ![
            enthalpy_change_j,
            volume_change_m3,
            boundary_work_j,
            internal_energy_change_j,
            entropy_change_j_k,
            heat_residual_j,
            first_law_residual_j,
        ]
        .into_iter()
        .all(f64::is_finite)
            || !volumes_m3.into_iter().all(|v| v.is_finite() && v > 0.0)
            || (supplied_heat_j != 0.0
                && (enthalpy_change_j == 0.0
                    || vapor_mass_change == 0.0
                    || volume_change_m3 == 0.0
                    || boundary_work_j == 0.0
                    || internal_energy_change_j == 0.0))
        {
            return Err(Error::Unrepresentable(
                "heat, phase inventory or pressure work",
            ));
        }
        for (state, masses) in [(current, initial_phase_masses_kg), (state, phase_masses_kg)] {
            if (state.liquid_mass_fraction() > 0.0 && masses[0] == 0.0)
                || (state.vapor_mass_fraction() > 0.0 && masses[1] == 0.0)
            {
                return Err(Error::Unrepresentable("phase mass"));
            }
        }
        if heat_residual_j.abs() > energy_tolerance_j
            || first_law_residual_j.abs() > energy_tolerance_j
        {
            return Err(Error::EnergyBalance {
                heat_residual_j,
                first_law_residual_j,
            });
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        Ok(WaterVaporizationHeatStep {
            initial: current,
            state,
            mass_kg,
            supplied_heat_j,
            initial_phase_masses_kg,
            phase_masses_kg,
            volumes_m3,
            volume_change_m3,
            enthalpy_change_j,
            boundary_work_j,
            internal_energy_change_j,
            entropy_change_j_k,
            heat_residual_j,
            first_law_residual_j,
        })
    }
}

#[cfg(test)]
mod tests;
