//! Source-integrated isotropic thermoelasticity around the original J2 owner.
//!
//! Let f(T) = integral(alpha(theta) dtheta), a(T) = alpha(T) I and
//! e = eps - eps_p - f(T) I. With fixed C and H, the Helmholtz contribution
//! is psi = e:C:e/2 + H p^2/2, its entropy is s = a(T):sigma, and its
//! internal energy is psi + T s. The native backward-Euler work is
//! sigma_1:delta(eps-f I), with its original physical and numerical losses.
//! Hence reversible caloric transfer is sigma_1:I delta(f) - delta(T s).
//! Constant alpha recovers -T_0 alpha delta(tr(sigma)). This endpoint rule
//! closes the discrete first law; time refinement is still required.
//! Numerical loss remains separate and is never deposited as physical heat.

use super::*;
use crate::tensor::{contract, trace};
use crate::{MaterialAdmissibility, Tangent6};
use fs_evidence::ModelCard;

/// Native J2 law evaluated at the actual temperature, with its sourced
/// stress-free strain offset. It uses the original return map and tangent.
#[derive(Clone, Debug)]
pub struct ThermalJ2 {
    source: ResolvedJ2Plasticity,
    expansion_strain: f64,
}
impl ThermalJ2 {
    /// Complete original property bundle at this law's temperature.
    pub const fn resolved(&self) -> &ResolvedJ2Plasticity {
        &self.source
    }

    fn mechanical_strain(&self, strain: &Voigt) -> Voigt {
        core::array::from_fn(|i| strain[i] - if i < 3 { self.expansion_strain } else { 0.0 })
    }
}
impl SmallStrainLaw for ThermalJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.source.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.source.stress(&self.mechanical_strain(strain), state)
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        self.source.tangent(&self.mechanical_strain(strain), state)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        self.source
            .update_state(&self.mechanical_strain(strain), state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.source.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.source.card()
    }
}

/// One native J2 increment and the additional reversible thermal accounts.
#[derive(Clone, Debug)]
pub struct ThermoplasticIncrement {
    /// Native mechanical strain eps-integral(alpha dT) and J2 work/storage.
    pub mechanical: EnergyIncrement,
    /// Endpoint stress contracted with the TOTAL strain increment [J/m3].
    pub mechanical_work_j_m3: f64,
    /// Signed reversible exchange into the stress-free caloric chart [J/m3].
    pub reversible_heat_j_m3: f64,
    /// Change in T alpha(T) tr(sigma), separate from elastic free energy [J/m3].
    pub thermoelastic_storage_change_j_m3: f64,
}
impl ThermoplasticIncrement {
    /// Plastic heat plus reversible exchange, supplied once to heat transport.
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.mechanical.energy.physical_dissipation_j_m3 + self.reversible_heat_j_m3
    }
}

impl ThermoplasticSource<'_> {
    /// Sourced expansion coefficient at the reference temperature [1/K].
    /// It is zero when the caller explicitly selects the zero-expansion reduction.
    /// Evolving values remain available in each resolved material snapshot.
    pub const fn expansion_per_k(&self) -> f64 {
        self.expansion_per_k
    }

    /// Native elastic/plastic strain argument for a TOTAL strain and actual T.
    /// Both total and mechanical strains must remain inside the declared bound.
    pub fn mechanical_strain(
        &self,
        temperature_k: f64,
        strain: &Voigt,
    ) -> Result<Voigt, FeedbackError> {
        self.mechanical_strain_at(temperature_k, strain)
    }

    pub(super) fn mechanical_strain_at<E>(
        &self,
        temperature_k: f64,
        strain: &Voigt,
    ) -> Result<Voigt, FeedbackError<E>> {
        let expansion =
            self.expansion_strain_between(self.reference_temperature_k, temperature_k)?;
        let result = core::array::from_fn(|i| strain[i] - if i < 3 { expansion } else { 0.0 });
        let limit = self.source.law().elastic.strain_limit;
        for value in [strain, &result] {
            let norm_squared = contract(value, value);
            if value.iter().any(|v| !v.is_finite())
                || !norm_squared.is_finite()
                || norm_squared > limit * limit
            {
                return Err(FeedbackError::Invalid(
                    "total or thermally offset strain exceeds the small-strain bound",
                ));
            }
        }
        Ok(result)
    }

    /// Resolve the actual caloric state into the same native spatial law,
    /// including its material-derived stress-free strain.
    pub fn resolve_expanded_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalJ2), FeedbackError> {
        let (state, source) = self.resolve_enthalpy(specific_enthalpy_j_kg)?;
        let expansion_strain =
            self.expansion_strain_between(self.reference_temperature_k, state.temperature_k())?;
        Ok((
            state,
            ThermalJ2 {
                source,
                expansion_strain,
            },
        ))
    }

    /// Audit a total-strain increment against the original source and caloric
    /// states. Accepted history remains unchanged. No thermal solver is added.
    pub fn advance_with_energy(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ThermoplasticIncrement, FeedbackError> {
        let old = self
            .caloric
            .curve()
            .state_at_specific_enthalpy(old_h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let new = self
            .caloric
            .curve()
            .state_at_specific_enthalpy(new_h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        Ok(self
            .increment_at_temperature(
                old.temperature_k(),
                new.temperature_k(),
                old_strain,
                old_state,
                strain,
                tolerances,
            )?
            .1)
    }

    pub(super) fn increment_at_temperature<E>(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(ResolvedJ2Plasticity, ThermoplasticIncrement), FeedbackError<E>> {
        let source = self.resolve_temperature(temperature_k)?;
        let old_mechanical = self.mechanical_strain_at(old_temperature_k, old_strain)?;
        let mechanical_strain = self.mechanical_strain_at(temperature_k, strain)?;
        let mechanical = source
            .law()
            .advance_with_energy(&old_mechanical, old_state, &mechanical_strain, tolerances)
            .map_err(FeedbackError::Mechanical)?;
        let old_elastic = core::array::from_fn(|i| old_mechanical[i] - old_state.plastic_strain[i]);
        // Recover the OLD elastic stress directly. The endpoint yield offset
        // must not return-map the already committed state a second time.
        let old_stress = source.law().elastic.stress(&old_elastic, &());
        let delta_strain = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let mechanical_work_j_m3 = finite(contract(&mechanical.stress_pa, &delta_strain))?;
        let (old_alpha, new_alpha) = if self.expansion_knots.is_empty() {
            (self.expansion_per_k, self.expansion_per_k)
        } else {
            // Read both endpoint coefficients from the original pinned claims.
            // This refreshes source data, never the old committed plastic history.
            let old_source = self.resolve_temperature(old_temperature_k)?;
            let coefficient = |law: &ResolvedJ2Plasticity| {
                law.material()
                    .property(crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
                    .expect("the admitted variable expansion source is present")
                    .value_si()
            };
            (coefficient(&old_source), coefficient(&source))
        };
        let old_entropy = product(old_alpha, trace(&old_stress))?;
        let new_entropy = product(new_alpha, trace(&mechanical.stress_pa))?;
        let thermoelastic_storage_change_j_m3 = finite(
            product(temperature_k, new_entropy)? - product(old_temperature_k, old_entropy)?,
        )?;
        let reversible_heat_j_m3 = if self.expansion_knots.is_empty() {
            // Keep the original constant-coefficient reduction and arithmetic.
            product(-old_temperature_k, finite(new_entropy - old_entropy)?)?
        } else {
            // Integrate THIS interval, avoiding cancellation of two large
            // primitives relative to the original stress-free temperature.
            let delta_expansion =
                self.expansion_strain_between(old_temperature_k, temperature_k)?;
            finite(
                product(trace(&mechanical.stress_pa), delta_expansion)?
                    - thermoelastic_storage_change_j_m3,
            )?
        };
        let residual = finite(
            mechanical_work_j_m3
                - mechanical.energy.stored_change_j_m3()
                - thermoelastic_storage_change_j_m3
                - mechanical.energy.physical_dissipation_j_m3
                - reversible_heat_j_m3
                - mechanical.energy.numerical_dissipation_j_m3,
        )?;
        if residual.abs() > tolerances.energy_j_m3 {
            return Err(FeedbackError::Mechanical(EnergyError::EnergyResidual {
                residual_j_m3: residual,
                tolerance_j_m3: tolerances.energy_j_m3,
            }));
        }
        Ok((
            source,
            ThermoplasticIncrement {
                mechanical,
                mechanical_work_j_m3,
                reversible_heat_j_m3,
                thermoelastic_storage_change_j_m3,
            },
        ))
    }

    fn expansion_strain_between<E>(
        &self,
        from_k: f64,
        to_k: f64,
    ) -> Result<f64, FeedbackError<E>> {
        let knots = self.caloric.curve().knots();
        let lower = knots[0].temperature_k;
        let upper = knots[knots.len() - 1].temperature_k;
        if [from_k, to_k]
            .iter()
            .any(|&t| !t.is_finite() || t < lower || t > upper)
        {
            return Err(FeedbackError::Invalid(
                "thermal strain path is outside the original caloric support",
            ));
        }
        if self.expansion_knots.is_empty() {
            product(self.expansion_per_k, to_k - from_k)
        } else {
            finite(
                crate::state_point::integrate_piecewise_linear(
                    &self.expansion_knots,
                    from_k,
                    to_k,
                )
                .map_err(|e| FeedbackError::Source(e.into()))?,
            )
        }
    }

    pub(super) fn admit_thermoelastic_capacity(&self) -> Result<(), FeedbackError> {
        if !self.expansion_knots.is_empty() {
            return self.admit_variable_expansion_capacity();
        }
        if self.expansion_per_k == 0.0 {
            return Ok(());
        }
        let law = self.source.law();
        let bulk = law.elastic.youngs / (3.0 * (1.0 - 2.0 * law.elastic.poisson));
        let correction_per_k =
            9.0 * bulk * self.expansion_per_k * self.expansion_per_k / self.density_kg_m3;
        if !positive(correction_per_k) {
            return Err(FeedbackError::Invalid(
                "thermoelastic heat capacity correction is not representable",
            ));
        }
        let cp_name = &self.caloric.receipts()[0].property;
        let knots = self.caloric.curve().knots();
        for knot in knots {
            let source = self.resolve_temperature(knot.temperature_k)?;
            let cp = source
                .material()
                .property(cp_name)
                .ok_or(FeedbackError::Invalid(
                    "missing original heat capacity source",
                ))?
                .value_si();
            if !positive(cp - correction_per_k * knot.temperature_k) {
                return Err(FeedbackError::Invalid(
                    "sourced isochoric heat capacity must be positive over the complete thermal span",
                ));
            }
        }
        // The admitted thermal owner uses piecewise-linear h(T), with its own
        // declared interpolation error. Its actual segment slopes must also
        // stay positive after reversible thermoelastic coupling.
        for pair in knots.windows(2) {
            let capacity = (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                / (pair[1].temperature_k - pair[0].temperature_k);
            if !positive(capacity - correction_per_k * pair[1].temperature_k) {
                return Err(FeedbackError::Invalid(
                    "caloric chart isochoric slope must be positive on every segment",
                ));
            }
        }
        Ok(())
    }

    fn admit_variable_expansion_capacity(&self) -> Result<(), FeedbackError> {
        let law = self.source.law();
        let bulk = law.elastic.youngs / (3.0 * (1.0 - 2.0 * law.elastic.poisson));
        let maximum_alpha = self
            .expansion_knots
            .iter()
            .map(|(_, alpha)| alpha.abs())
            .fold(0.0_f64, f64::max);
        let mut maximum_slope = 0.0_f64;
        for pair in self.expansion_knots.windows(2) {
            let width = pair[1].0 - pair[0].0;
            if !positive(width) {
                return Err(FeedbackError::Invalid(
                    "expansion source temperatures must be strictly increasing",
                ));
            }
            let slope = finite((pair[1].1 - pair[0].1) / width)?;
            maximum_slope = maximum_slope.max(slope.abs());
        }

        // At fixed total strain and plastic history:
        // c_eff = rho c_caloric + T alpha'(T) tr(sigma) - 9 K T alpha(T)^2.
        // Native energy::admit_energy_history bounds the ACTUAL elastic tensor
        // e by ||e|| <= strain_limit under tensor::contract: diagonal squares
        // plus twice the shear squares. Thus |tr(e)| <= sqrt(3) ||e|| and
        // |tr(sigma)| = 3 K |tr(e)| <= 3 sqrt(3) K strain_limit. No assumption
        // about an exactly zero floating-point plastic trace is required.
        let maximum_trace_stress = product(
            product(3.0 * fs_math::det::sqrt(3.0), bulk)?,
            law.elastic.strain_limit,
        )?;
        let correction_per_k = finite(
            product(9.0 * bulk, product(maximum_alpha, maximum_alpha)?)?
                + product(maximum_slope, maximum_trace_stress)?,
        )? / self.density_kg_m3;
        let knots = self.caloric.curve().knots();
        let correction = product(knots[knots.len() - 1].temperature_k, correction_per_k)?;
        // Piecewise-linear alpha has its largest |alpha| at a knot; each
        // segment's derivative is constant. The full-span bound is sufficient,
        // deliberately conservative, and includes both one-sided knot limits.
        let mut minimum_capacity = self.minimum_source_heat_capacity_j_kg_k;
        for pair in knots.windows(2) {
            let capacity = (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                / (pair[1].temperature_k - pair[0].temperature_k);
            minimum_capacity = minimum_capacity.min(capacity);
        }
        if !positive(minimum_capacity) || !positive(minimum_capacity - correction) {
            return Err(FeedbackError::Invalid(
                "variable expansion cannot ensure positive effective heat capacity over the admitted strain and temperature span",
            ));
        }
        Ok(())
    }
}
