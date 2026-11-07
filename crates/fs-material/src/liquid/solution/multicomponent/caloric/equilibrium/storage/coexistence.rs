//! Two original equilibrium branches joined by verified liquid coexistence.
//!
//! The coexistence owner admits the actual common temperature, every chemical
//! potential, conserved feed and original latent H/S. The plateau is an exact
//! isothermal enthalpy interval, with dT/dh=0, not an apparent-Cp smoothing.
//! The lower-enthalpy coexistence endpoint is h=0: the cold sensible branch
//! consequently has negative relative enthalpy. This keeps a small latent heat
//! out of an unrelated large cold-reference offset.
//!
//! ThermalLaw returns its defining caloric integral dh/T. The separate source
//! entropy query retains the original phase entropy, including the measured
//! coexistence closure residual; no source value is fitted to a heat budget.
//! The native entropy test does not turn a numerically admitted common tangent
//! into a global Gibbs-minimum, kinetic or phase morphology certificate.

/// Declared bulk conduction on the original sensible/latent storage path.
pub mod conduction;

use super::{
    IsobaricHeatCapacityLaw, MulticomponentEquilibriumCaloricStateV1,
    MulticomponentLiquidCaloricError, MulticomponentLiquidEquilibriumThermalStorage, Result,
    Temperature, poll, positive,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use fs_thermochem::solution::multicomponent::caloric::equilibrium::coexistence::{
    MulticomponentCaloricCoexistenceControlV1, MulticomponentCaloricCoexistenceV1,
    MulticomponentCoexistenceStateV1,
};
use std::sync::Arc;

/// One continuous same-feed thermal path with a genuine common-tangent plateau.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidCoexistenceThermalStorage {
    branches: [MulticomponentLiquidEquilibriumThermalStorage; 2],
    coexistence: MulticomponentCaloricCoexistenceV1,
    enthalpy_bounds: [f64; 2],
    latent_h_j_kg: f64,
    maximum_slope: f64,
}

impl MulticomponentLiquidEquilibriumThermalStorage {
    /// Join two adjacent original phase branches using their actual endpoint
    /// chemical potentials and caloric storage. The upper temperature of self
    /// must equal the lower temperature of upper. Both laws must retain the
    /// SAME material, pure Cp, component masses, pressure and conserved feed.
    /// General phase discovery is not used to replace an incompatible seam.
    /// The joined path has its own coexistence-relative enthalpy coordinate;
    /// preparing either branch for streams does not propagate its flow binding.
    pub fn join_coexistence(
        &self,
        upper: &Self,
        controls: MulticomponentCaloricCoexistenceControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidCoexistenceThermalStorage> {
        poll(&mut cancelled)?;
        if (!Arc::ptr_eq(&self.owner, &upper.owner) && self.owner != upper.owner)
            || self.component_moles() != upper.component_moles()
            || self.temperature_domain()[1] != upper.temperature_domain()[0]
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "coexistence storage needs adjacent branches of the same original source and feed",
            ));
        }
        let coexistence = self.owner.model.coexistence_plateau(
            &self.endpoints[1],
            &upper.endpoints[0],
            controls,
            &mut cancelled,
        )?;
        let mass = self.mass_kg;
        let latent_h_j_kg = coexistence.latent_enthalpy_j() / mass;
        let lower_span = self
            .owner
            .model
            .equilibrium_enthalpy_change_j(&self.endpoints[0], &self.endpoints[1])?
            / mass;
        let upper_span = self
            .owner
            .model
            .equilibrium_enthalpy_change_j(&upper.endpoints[0], &upper.endpoints[1])?
            / mass;
        // Each branch may carry a different stream or closed-storage datum.
        // Only original endpoint differences define this joined heat coordinate.
        let enthalpy_bounds = [-lower_span, latent_h_j_kg + upper_span];
        // Both storages normalize the same conserved feed. Use this path's
        // original lower mass even when their phase summations differ by ulps.
        let upper_slope = (upper.maximum_slope * (mass / upper.mass_kg).next_up()).next_up();
        let maximum_slope = self.maximum_slope.max(upper_slope);
        if !positive(latent_h_j_kg)
            || !positive(lower_span)
            || !positive(upper_span)
            || !positive(maximum_slope)
            || !enthalpy_bounds[1].is_finite()
            || enthalpy_bounds[1] <= latent_h_j_kg
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable coexistence storage interval or temperature slope",
            ));
        }
        let coordinate_error = ((enthalpy_bounds[1] - latent_h_j_kg) - upper_span).abs() * mass;
        let coordinate_allowance = self
            .controls
            .thermal
            .energy_tolerance_j
            .min(upper.controls.thermal.energy_tolerance_j)
            + 64.0 * f64::EPSILON * upper_span * mass;
        if !coordinate_error.is_finite() || coordinate_error > coordinate_allowance {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "upper sensible span loses its original heat accuracy beside latent storage",
            ));
        }
        let storage = MulticomponentLiquidCoexistenceThermalStorage {
            branches: [self.clone(), upper.clone()],
            coexistence,
            enthalpy_bounds,
            latent_h_j_kg,
            maximum_slope,
        };
        poll(&mut cancelled)?;
        Ok(storage)
    }
}

impl MulticomponentLiquidCoexistenceThermalStorage {
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        [
            self.branches[0].temperature_domain()[0],
            self.branches[1].temperature_domain()[1],
        ]
    }
    /// Inclusive relative h domain [J/kg]; the lower coexistence endpoint is zero.
    pub const fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        self.enthalpy_bounds
    }
    pub const fn coexistence_specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [0.0, self.latent_h_j_kg]
    }
    pub const fn coexistence(&self) -> &MulticomponentCaloricCoexistenceV1 {
        &self.coexistence
    }
    pub fn mass_kg(&self) -> f64 {
        self.branches[0].mass_kg()
    }
    pub fn component_moles(&self) -> &[f64] {
        self.coexistence.component_moles()
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        self.branches[0].sources()
    }

    /// Original branch Cp away from coexistence. The common temperature has
    /// a finite enthalpy interval, not a finite single-valued heat capacity.
    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        let state = self.branch_at_temperature(temperature, &mut cancelled)?;
        let capacity = state.isobaric_heat_capacity_j_per_k() / self.mass_kg();
        if !positive(capacity) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable coexistence-branch specific heat",
            ));
        }
        Ok(capacity)
    }

    /// Original enthalpy away from the transition. At its common temperature,
    /// callers must choose within coexistence_specific_enthalpy_domain_j_kg;
    /// this query refuses instead of inventing a phase fraction.
    pub fn specific_enthalpy_j_kg(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        let state = self.branch_at_temperature(temperature, &mut cancelled)?;
        if temperature == self.temperature_domain()[0] {
            return Ok(self.enthalpy_bounds[0]);
        }
        if temperature == self.temperature_domain()[1] {
            return Ok(self.enthalpy_bounds[1]);
        }
        let upper = temperature.value() > self.coexistence.temperature().value();
        let anchor = if upper {
            self.coexistence.upper()
        } else {
            self.coexistence.lower()
        };
        let increment = self.branches[0]
            .owner
            .model
            .equilibrium_enthalpy_change_j(anchor, &state)?
            / self.mass_kg();
        let h = if upper {
            self.latent_h_j_kg + increment
        } else {
            increment
        };
        self.admit_enthalpy(h)?;
        if (upper && h <= self.latent_h_j_kg) || (!upper && h >= 0.0) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "nonzero sensible heat is unrepresentable beside coexistence",
            ));
        }
        poll(&mut cancelled)?;
        Ok(h)
    }

    /// Original sensible inverse or the exact latent lever rule. At either
    /// plateau corner the selected one-sided Newton tangent is zero; there is
    /// no unique two-sided derivative at that corner.
    pub fn caloric_state_at_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidCoexistenceThermalState<'_>> {
        poll(&mut cancelled)?;
        self.admit_enthalpy(h)?;
        let point = if h < 0.0 {
            CaloricPoint::Branch {
                state: self.branch_at_enthalpy(h, false, &mut cancelled)?,
                upper: false,
            }
        } else if h > self.latent_h_j_kg {
            CaloricPoint::Branch {
                state: self.branch_at_enthalpy(h, true, &mut cancelled)?,
                upper: true,
            }
        } else {
            let heat = if h == self.latent_h_j_kg {
                self.coexistence.latent_enthalpy_j()
            } else {
                h * self.mass_kg()
            };
            if h > 0.0 && !positive(heat) {
                return Err(MulticomponentLiquidCaloricError::Invalid(
                    "positive original latent heat is unrepresentable after mass normalization",
                ));
            }
            CaloricPoint::Coexisting(
                self.coexistence
                    .state_at_enthalpy_increment(heat, &mut cancelled)?,
            )
        };
        let response = MulticomponentLiquidCoexistenceThermalState {
            point,
            mass_kg: self.mass_kg(),
        };
        let slope = response.temperature_slope_kg_k_j();
        if !slope.is_finite() || (response.branch_state().is_some() && slope <= 0.0) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable original coexistence-branch temperature slope",
            ));
        }
        poll(&mut cancelled)?;
        Ok(response)
    }

    /// Actual source phase entropy, separately from ThermalLaw's heat integral.
    /// Their difference is exactly the admitted portion of the original latent
    /// Clausius residual, up to floating-point normalization and summation.
    pub fn specific_source_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        self.entropy_change(initial_h, final_h, true, &mut cancelled)
    }

    fn branch_at_temperature(
        &self,
        temperature: Temperature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(cancelled)?;
        let [lo, hi] = self.temperature_domain();
        if !temperature.value().is_finite()
            || temperature.value() < lo.value()
            || temperature.value() > hi.value()
            || temperature == self.coexistence.temperature()
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "temperature is outside the joined domain or has set-valued coexistence enthalpy",
            ));
        }
        if temperature == lo {
            return Ok(self.branches[0].endpoints[0].clone());
        }
        if temperature == hi {
            return Ok(self.branches[1].endpoints[1].clone());
        }
        let upper = temperature.value() > self.coexistence.temperature().value();
        let branch = &self.branches[usize::from(upper)];
        let anchor = if upper {
            self.coexistence.upper()
        } else {
            self.coexistence.lower()
        };
        Ok(branch.owner.model.equilibrium_at_temperature(
            anchor,
            temperature,
            branch.controls,
            cancelled,
        )?)
    }

    fn branch_at_enthalpy(
        &self,
        h: f64,
        upper: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(cancelled)?;
        let branch = &self.branches[usize::from(upper)];
        if h == self.enthalpy_bounds[usize::from(upper)] {
            return Ok(branch.endpoints[usize::from(upper)].clone());
        }
        let anchor = if upper {
            self.coexistence.upper()
        } else {
            self.coexistence.lower()
        };
        let offset = if upper { self.latent_h_j_kg } else { 0.0 };
        if h == offset {
            return Ok(anchor.clone());
        }
        let heat = (h - offset) * self.mass_kg();
        if !heat.is_finite() || heat == 0.0 {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "nonzero original sensible heat beside coexistence is unrepresentable",
            ));
        }
        let balance =
            branch
                .owner
                .model
                .add_heat_equilibrium(anchor, heat, branch.controls, cancelled)?;
        Ok(balance.final_state().clone())
    }

    fn entropy_change(
        &self,
        initial_h: f64,
        final_h: f64,
        source: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(cancelled)?;
        self.admit_enthalpy(initial_h)?;
        self.admit_enthalpy(final_h)?;
        if initial_h == final_h {
            return Ok(0.0);
        }
        let (lo, hi, sign) = if initial_h < final_h {
            (initial_h, final_h, 1.0)
        } else {
            (final_h, initial_h, -1.0)
        };
        let mut entropy = 0.0;
        for (upper, start, stop) in [
            (false, lo.min(0.0), hi.min(0.0)),
            (true, lo.max(self.latent_h_j_kg), hi.max(self.latent_h_j_kg)),
        ] {
            if start == stop {
                continue;
            }
            let before = self.branch_at_enthalpy(start, upper, cancelled)?;
            let after = self.branch_at_enthalpy(stop, upper, cancelled)?;
            let increment = self.branches[0]
                .owner
                .model
                .equilibrium_entropy_change_j_per_k(&before, &after)?
                / self.mass_kg();
            if !positive(increment) {
                return Err(MulticomponentLiquidCaloricError::Invalid(
                    "positive original sensible entropy beside coexistence is unresolved",
                ));
            }
            entropy += increment;
        }
        let latent_start = lo.max(0.0);
        let latent_stop = hi.min(self.latent_h_j_kg);
        if latent_stop > latent_start {
            let width = latent_stop - latent_start;
            entropy += if source {
                (width / self.latent_h_j_kg)
                    * (self.coexistence.latent_entropy_j_per_k() / self.mass_kg())
            } else {
                width / self.coexistence.temperature().value()
            };
        }
        if !positive(entropy) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "positive coexistence path entropy is unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(sign * entropy)
    }

    fn admit_enthalpy(&self, h: f64) -> Result<()> {
        if !h.is_finite() || h < self.enthalpy_bounds[0] || h > self.enthalpy_bounds[1] {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "enthalpy outside the original coexistence storage domain",
            ));
        }
        Ok(())
    }
}

impl ThermalLaw for MulticomponentLiquidCoexistenceThermalStorage {
    fn temperature_k(&self, h: f64) -> core::result::Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| state.temperature().value())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> core::result::Result<Option<f64>, DepositionError> {
        self.entropy_change(initial_h, final_h, false, &mut || false)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

#[derive(Clone, Debug, PartialEq)]
enum CaloricPoint<'a> {
    Branch {
        state: MulticomponentEquilibriumCaloricStateV1,
        upper: bool,
    },
    Coexisting(MulticomponentCoexistenceStateV1<'a>),
}

/// Temperature, tangent and actual phase inventory from the same heat coordinate.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidCoexistenceThermalState<'a> {
    point: CaloricPoint<'a>,
    mass_kg: f64,
}
impl<'a> MulticomponentLiquidCoexistenceThermalState<'a> {
    pub fn temperature(&self) -> Temperature {
        match &self.point {
            CaloricPoint::Branch { state, .. } => state.temperature(),
            CaloricPoint::Coexisting(state) => state.temperature(),
        }
    }
    pub fn temperature_slope_kg_k_j(&self) -> f64 {
        match &self.point {
            CaloricPoint::Branch { state, .. } => {
                self.mass_kg / state.isobaric_heat_capacity_j_per_k()
            }
            CaloricPoint::Coexisting(_) => 0.0,
        }
    }
    /// Fraction of the higher-enthalpy endpoint ensemble. Both ensembles are
    /// liquids; this is not a solid/liquid or vapor quality interpretation.
    pub fn upper_branch_mass_fraction(&self) -> f64 {
        match &self.point {
            CaloricPoint::Branch { upper, .. } => {
                if *upper {
                    1.0
                } else {
                    0.0
                }
            }
            CaloricPoint::Coexisting(state) => state.upper_ensemble_fraction(),
        }
    }
    pub fn branch_state(&self) -> Option<&MulticomponentEquilibriumCaloricStateV1> {
        match &self.point {
            CaloricPoint::Branch { state, .. } => Some(state),
            CaloricPoint::Coexisting(_) => None,
        }
    }
    pub fn coexistence_state(&self) -> Option<&MulticomponentCoexistenceStateV1<'a>> {
        match &self.point {
            CaloricPoint::Coexisting(state) => Some(state),
            CaloricPoint::Branch { .. } => None,
        }
    }
}
