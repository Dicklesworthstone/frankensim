//! Source-integrated isotropic thermoelasticity around the original J2 owner.
//!
//! Let f(T) = integral(alpha(theta) dtheta), a(T) = alpha(T) I and
//! e = eps - eps_p - f(T) I. The Helmholtz contribution is
//! psi = e:C(T):e/2 + H(T) p^2/2 + H_k dev(eps_p):dev(eps_p)/3.
//! The optional Prager modulus H_k is fixed over the admitted thermal chart
//! and contributes storage but no temperature derivative. With fixed nu and affine E/H its entropy is
//! s = a(T):sigma - (E'/E) psi_elastic - H' p^2/2; internal energy is psi + T s.
//! The native backward-Euler return map freezes coefficients at T_1. Retain its
//! original physical/numerical losses and add the frozen-history coefficient
//! change Dpsi_c = psi(T_1,e_0,p_0)-psi(T_0,e_0,p_0) to its stored-energy change.
//! Reversible caloric transfer is sigma_1:I delta(f) - Dpsi_c - delta(T s).
//! Fixed moduli/constant alpha retain -T_0 alpha delta(tr(sigma)) exactly.
//! This endpoint rule closes the discrete first law; time refinement is required.
//! Numerical loss remains separate and is never deposited as physical heat.
//! Selected linear overstress uses eta(T_1)/dt in the original return map and
//! physical dissipation. Eta stores no free energy and adds no entropy term.

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
    rate: Option<super::super::rate::LinearOverstressJ2>,
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
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.stress(&strain, state),
            None => self.source.stress(&strain, state),
        }
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.tangent(&strain, state),
            None => self.source.tangent(&strain, state),
        }
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.update_state(&strain, state),
            None => self.source.update_state(&strain, state),
        }
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.rate.as_ref().map_or_else(|| self.source.admissibility(), SmallStrainLaw::admissibility)
    }
    fn card(&self) -> ModelCard {
        self.rate.as_ref().map_or_else(|| self.source.card(), SmallStrainLaw::card)
    }
}

/// One native J2 increment and the additional reversible thermal accounts.
#[derive(Clone, Debug)]
pub struct ThermoplasticIncrement {
    /// Native mechanical strain eps-integral(alpha dT) and endpoint-modulus
    /// J2 work/storage. Both native energy endpoints use the NEW moduli.
    pub mechanical: EnergyIncrement,
    /// Frozen old strain/history free-energy change from sourced E(T)/H(T).
    /// Adding this to the native stored change gives the actual two-temperature
    /// free-energy change [J/m3]. This is neither plastic nor numerical loss.
    pub coefficient_storage_change_j_m3: f64,
    /// Endpoint stress contracted with the TOTAL strain increment [J/m3].
    pub mechanical_work_j_m3: f64,
    /// Signed reversible exchange into the stress-free caloric chart [J/m3].
    pub reversible_heat_j_m3: f64,
    /// Change in T times mechanical entropy, separate from free energy [J/m3].
    pub thermoelastic_storage_change_j_m3: f64,
}
impl ThermoplasticIncrement {
    /// Actual elastic/hardening free-energy change between the two temperatures.
    /// The native endpoint-modulus account remains independently inspectable.
    pub fn stored_change_j_m3(&self) -> f64 {
        self.mechanical.energy.stored_change_j_m3() + self.coefficient_storage_change_j_m3
    }

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
    /// Rate-dependent sources require `resolve_expanded_enthalpy_for_step`.
    pub fn resolve_expanded_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalJ2), FeedbackError> {
        self.resolve_expanded_for_interval(specific_enthalpy_j_kg, None)
    }

    /// Resolve source coefficients and the return map for this physical interval.
    /// Thermal substeps and nonlinear iterations do not create extra plastic time.
    pub fn resolve_expanded_enthalpy_for_step(&self, specific_enthalpy_j_kg: f64, duration_s: f64)
        -> Result<(EquilibriumPhaseState, ThermalJ2), FeedbackError> {
        if !positive(duration_s) {
            return Err(FeedbackError::Invalid("positive physical constitutive interval required"));
        }
        self.resolve_expanded_for_interval(specific_enthalpy_j_kg, Some(duration_s))
    }

    fn resolve_expanded_for_interval(&self, specific_enthalpy_j_kg: f64, duration_s: Option<f64>)
        -> Result<(EquilibriumPhaseState, ThermalJ2), FeedbackError> {
        let (state, source) = self.resolve_enthalpy(specific_enthalpy_j_kg)?;
        let rate = self.overstress_law(&source, duration_s)?;
        let expansion_strain =
            self.expansion_strain_between(self.reference_temperature_k, state.temperature_k())?;
        Ok((
            state,
            ThermalJ2 {
                source,
                expansion_strain,
                rate,
            },
        ))
    }

    /// Audit a total-strain increment against the original source and caloric
    /// states. Accepted history remains unchanged. No thermal solver is added.
    /// Rate-dependent sources require `advance_with_energy_for_step`.
    pub fn advance_with_energy(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ThermoplasticIncrement, FeedbackError> {
        self.advance_for_interval(old_h_j_kg, new_h_j_kg, old_strain, old_state,
            strain, tolerances, None)
    }

    /// Audit yield/viscous dissipation and reciprocal heat at the same physical dt
    /// used by the native spatial constitutive response.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_energy_for_step(&self, old_h_j_kg: f64, new_h_j_kg: f64,
        old_strain: &Voigt, old_state: &J2State, strain: &Voigt,
        tolerances: EnergyTolerances, duration_s: f64)
        -> Result<ThermoplasticIncrement, FeedbackError> {
        if !positive(duration_s) {
            return Err(FeedbackError::Invalid("positive physical constitutive interval required"));
        }
        self.advance_for_interval(old_h_j_kg, new_h_j_kg, old_strain, old_state,
            strain, tolerances, Some(duration_s))
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_for_interval(&self, old_h_j_kg: f64, new_h_j_kg: f64,
        old_strain: &Voigt, old_state: &J2State, strain: &Voigt,
        tolerances: EnergyTolerances, duration_s: Option<f64>)
        -> Result<ThermoplasticIncrement, FeedbackError> {
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
                duration_s,
            )?
            .1)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn increment_at_temperature<E>(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: Option<f64>,
    ) -> Result<(ResolvedJ2Plasticity, ThermoplasticIncrement), FeedbackError<E>> {
        let source = self.resolve_temperature(temperature_k)?;
        let old_mechanical = self.mechanical_strain_at(old_temperature_k, old_strain)?;
        let mechanical_strain = self.mechanical_strain_at(temperature_k, strain)?;
        let mechanical = match self.overstress_law(&source, duration_s)? {
            Some(rate) => rate.advance_with_energy(&old_mechanical, old_state, &mechanical_strain, tolerances),
            None => source.law().advance_with_energy(&old_mechanical, old_state, &mechanical_strain, tolerances),
        }.map_err(FeedbackError::Mechanical)?;
        let old_elastic = core::array::from_fn(|i| old_mechanical[i] - old_state.plastic_strain[i]);
        let variable_moduli = self.youngs_temperature_slope_pa_k != 0.0
            || self.hardening_temperature_slope_pa_k != 0.0;
        let old_source = if variable_moduli || !self.expansion_knots.is_empty() {
            Some(self.resolve_temperature(old_temperature_k)?)
        } else {
            None
        };
        // Recover OLD elastic stress at its actual temperature. Refreshing the
        // source must not return-map the committed plastic history a second time.
        let old_law = old_source.as_ref().unwrap_or(&source);
        let old_stress = old_law.law().elastic.stress(&old_elastic, &());
        let delta_strain = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let mechanical_work_j_m3 = finite(contract(&mechanical.stress_pa, &delta_strain))?;
        let (old_alpha, new_alpha) = if self.expansion_knots.is_empty() {
            (self.expansion_per_k, self.expansion_per_k)
        } else {
            let coefficient = |law: &ResolvedJ2Plasticity| {
                law.material()
                    .property(crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
                    .expect("the admitted variable expansion source is present")
                    .value_si()
            };
            (coefficient(old_law), coefficient(&source))
        };
        let mut old_entropy = product(old_alpha, trace(&old_stress))?;
        let mut new_entropy = product(new_alpha, trace(&mechanical.stress_pa))?;
        let coefficient_storage_change_j_m3 = if variable_moduli {
            let youngs = source.law().elastic.youngs;
            // Fixed nu makes C(T)=E(T)*C_unit. Native initial/final elastic
            // energies therefore give the exact tensor contraction with C'.
            let elastic_derivative = self.youngs_temperature_slope_pa_k / youngs;
            let old_p_squared = product(old_state.alpha, old_state.alpha)?;
            let new_p_squared = product(mechanical.state.alpha, mechanical.state.alpha)?;
            old_entropy = finite(
                old_entropy
                    - product(elastic_derivative, mechanical.energy.initial.elastic_j_m3)?
                    - product(0.5 * self.hardening_temperature_slope_pa_k, old_p_squared)?,
            )?;
            new_entropy = finite(
                new_entropy
                    - product(
                        elastic_derivative,
                        mechanical.energy.final_state.elastic_j_m3,
                    )?
                    - product(0.5 * self.hardening_temperature_slope_pa_k, new_p_squared)?,
            )?;
            // Use actual pinned endpoint coefficients. The correction changes
            // native frozen-coefficient delta(psi) into true endpoint delta(psi).
            finite(
                product(
                    (youngs - old_law.law().elastic.youngs) / youngs,
                    mechanical.energy.initial.elastic_j_m3,
                )? + product(
                    0.5 * (source.law().hardening - old_law.law().hardening),
                    old_p_squared,
                )?,
            )?
        } else {
            0.0
        };
        let thermoelastic_storage_change_j_m3 = finite(
            product(temperature_k, new_entropy)? - product(old_temperature_k, old_entropy)?,
        )?;
        let reversible_heat_j_m3 = if self.expansion_knots.is_empty() && !variable_moduli {
            // Keep the original constant-coefficient reduction and arithmetic.
            product(-old_temperature_k, finite(new_entropy - old_entropy)?)?
        } else {
            // Integrate THIS interval, avoiding cancellation of two large
            // primitives relative to the original stress-free temperature.
            let delta_expansion =
                self.expansion_strain_between(old_temperature_k, temperature_k)?;
            finite(
                product(trace(&mechanical.stress_pa), delta_expansion)?
                    - coefficient_storage_change_j_m3
                    - thermoelastic_storage_change_j_m3,
            )?
        };
        let stored_change_j_m3 =
            finite(mechanical.energy.stored_change_j_m3() + coefficient_storage_change_j_m3)?;
        let residual = finite(
            mechanical_work_j_m3
                - stored_change_j_m3
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
                coefficient_storage_change_j_m3,
                mechanical_work_j_m3,
                reversible_heat_j_m3,
                thermoelastic_storage_change_j_m3,
            },
        ))
    }

    fn overstress_law<E>(&self, source: &ResolvedJ2Plasticity, duration_s: Option<f64>)
        -> Result<Option<super::super::rate::LinearOverstressJ2>, FeedbackError<E>> {
        use super::super::rate::{LinearOverstressJ2, OVERSTRESS_VISCOSITY_PROPERTY};
        if !self.linear_overstress { return Ok(None); }
        let duration_s = duration_s.ok_or(FeedbackError::Invalid(
            "source-defined overstress requires the physical step duration"))?;
        let viscosity = source.material().property(OVERSTRESS_VISCOSITY_PROPERTY)
            .ok_or(FeedbackError::Invalid("resolved source omitted overstress viscosity"))?.value_si();
        let law = LinearOverstressJ2::new(source.law().clone(), viscosity, duration_s)
            .map_err(|e| FeedbackError::Source(e.into()))?;
        Ok(Some(law))
    }

    fn expansion_strain_between<E>(&self, from_k: f64, to_k: f64) -> Result<f64, FeedbackError<E>> {
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
                crate::state_point::integrate_piecewise_linear(&self.expansion_knots, from_k, to_k)
                    .map_err(|e| FeedbackError::Source(e.into()))?,
            )
        }
    }

    pub(super) fn admit_thermoelastic_capacity(&self) -> Result<(), FeedbackError> {
        if !self.expansion_knots.is_empty() || self.youngs_temperature_slope_pa_k != 0.0 {
            return self.admit_variable_thermoelastic_capacity();
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

    fn admit_variable_thermoelastic_capacity(&self) -> Result<(), FeedbackError> {
        let law = self.source.law();
        let bulk_denominator = 3.0 * (1.0 - 2.0 * law.elastic.poisson);
        let bulk = self.maximum_youngs_modulus_pa / bulk_denominator;
        let bulk_slope = self.youngs_temperature_slope_pa_k / bulk_denominator;
        let maximum_alpha = self
            .expansion_knots
            .iter()
            .map(|(_, alpha)| alpha.abs())
            .fold(self.expansion_per_k.abs(), f64::max);
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
        // c_eff = rho c_caloric + T [alpha' tr(sigma)
        //          + 2 alpha tr(C' e) - 9 K alpha^2].
        // E/H are globally affine with fixed nu, so C''=0 and H''=0;
        // no unbounded hardening-history term is hidden in this capacity gate.
        // Native energy::admit_energy_history bounds the ACTUAL elastic tensor
        // e by ||e|| <= strain_limit under tensor::contract: diagonal squares
        // plus twice the shear squares. Thus |tr(e)| <= sqrt(3) ||e|| and
        // |tr(sigma)| = 3 K |tr(e)| <= 3 sqrt(3) K strain_limit. No assumption
        // about an exactly zero floating-point plastic trace is required.
        let maximum_trace_stress = product(
            product(3.0 * fs_math::det::sqrt(3.0), bulk)?,
            law.elastic.strain_limit,
        )?;
        let maximum_trace_stress_derivative = product(
            product(3.0 * fs_math::det::sqrt(3.0), bulk_slope.abs())?,
            law.elastic.strain_limit,
        )?;
        let correction_per_k = finite(
            product(9.0 * bulk, product(maximum_alpha, maximum_alpha)?)?
                + product(maximum_slope, maximum_trace_stress)?
                + product(2.0 * maximum_alpha, maximum_trace_stress_derivative)?,
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
                "sourced moduli and expansion cannot ensure positive effective heat capacity over the admitted strain and temperature span",
            ));
        }
        Ok(())
    }
}
