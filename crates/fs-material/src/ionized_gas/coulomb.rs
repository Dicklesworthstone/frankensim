//! Correlated charged-carrier velocities with original background drag anchors.

use super::{
    ELEMENTARY_CHARGE_C, IonizedGasConductivity, IonizedGasError, IonizedGasJouleHeating, poll,
    positive,
};
use fs_math::det::sqrt;
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombModelV1, CoulombPairV1};

const MAXIMUM_CARRIERS: usize = 128;
const BALANCE_RELATIVE_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct CoulombDrag {
    model: CoulombModelV1,
    pairs: Vec<CoulombPairV1>,
    // Heavy/electron friction and signed correction to electrons, per |E|².
    thermal_coefficients: [f64; 3],
    force_balance_relative_residual: f64,
    power_balance_relative_residual: f64,
}

impl IonizedGasConductivity {
    /// Add mutual Coulomb drag to this snapshot's original positive background
    /// anchors. The source bank, populations and volume must match exactly;
    /// all pair moments are evaluated at this snapshot's actual Th and Te.
    ///
    /// For every carrier, `a_s*u_s + sum_t K_st*(u_s-u_t) = n_s*q_s*E`,
    /// where `a_s=n_s*m_s*nu_background_s`. Every unordered charged pair is
    /// required, including ion-ion pairs. The symmetric system is factored once
    /// at construction, using the existing Cholesky owner and diagonal scaling.
    /// At most 128 carriers are admitted; original force and electrical-power
    /// relative residuals must be at most 1e-10. These are numerical admission
    /// controls, not empirical uncertainty. No neutral anchor is replaced or
    /// discarded, and a pure-Coulomb floating common velocity is not inferred.
    pub fn with_coulomb_drag(
        mut self,
        model: CoulombModelV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasError> {
        poll(&mut cancelled)?;
        let invalid = IonizedGasError::Invalid;
        let n = self.carriers.len();
        if self.coulomb.is_some()
            || !(2..=MAXIMUM_CARRIERS).contains(&n)
            || !model.is_bound_to_populations(&self.state)
        {
            return Err(invalid(
                "Coulomb drag requires an original unmodified anchored snapshot and matching populations",
            ));
        }
        let mut pairs = Vec::with_capacity(n * (n - 1) / 2);
        model
            .evaluate_with_pairs(
                self.state.heavy_temperature(),
                self.state.electron_temperature(),
                |pair| pairs.push(pair),
                &mut cancelled,
            )
            .map_err(source_error)?;
        if pairs.len() != n * (n - 1) / 2 {
            return Err(invalid("complete unordered charged-pair coverage required"));
        }
        pairs.sort_by_key(|p| (p.first_index(), p.second_index()));
        let mut original_indices = vec![usize::MAX; self.state.amounts().len()];
        for (i, source) in self.state.thermodynamics().models().iter().enumerate() {
            if let Ok(carrier) = self
                .carriers
                .binary_search_by(|carrier| carrier.species().cmp(source.species()))
            {
                original_indices[i] = carrier;
            }
        }
        let mut background_drag = Vec::with_capacity(n);
        let mut forcing = Vec::with_capacity(n);
        let mut matrix = vec![0.0; n * n];
        for (i, carrier) in self.carriers.iter().enumerate() {
            poll(&mut cancelled)?;
            let drag =
                carrier.number_density_per_m3 * carrier.mass_kg * carrier.relaxation.rate_per_s();
            let force = carrier.number_density_per_m3
                * ((carrier.charge_number as f64) * ELEMENTARY_CHARGE_C);
            if !positive(drag) || !force.is_finite() || force == 0.0 {
                return Err(invalid(
                    "background anchor or electric force is unrepresentable",
                ));
            }
            background_drag.push(drag);
            forcing.push(force);
            matrix[i * n + i] = drag;
        }
        let mut pair_carriers = Vec::with_capacity(pairs.len());
        for pair in &pairs {
            poll(&mut cancelled)?;
            let i = *original_indices
                .get(pair.first_index())
                .ok_or(invalid("Coulomb pair source index is absent"))?;
            let j = *original_indices
                .get(pair.second_index())
                .ok_or(invalid("Coulomb pair source index is absent"))?;
            let k = pair.friction_coefficient_kg_m3_s();
            if i >= n || j >= n || i == j || !positive(k) || matrix[i * n + j] != 0.0 {
                return Err(invalid(
                    "Coulomb pair is duplicate, absent or unrepresentable",
                ));
            }
            matrix[i * n + i] += k;
            matrix[j * n + j] += k;
            matrix[i * n + j] = -k;
            matrix[j * n + i] = -k;
            pair_carriers.push([i, j]);
        }
        let mut scales = Vec::with_capacity(n);
        for i in 0..n {
            let diagonal = matrix[i * n + i];
            if !positive(diagonal) {
                return Err(invalid("Coulomb friction diagonal is unrepresentable"));
            }
            scales.push(sqrt(diagonal));
        }
        let mut velocity = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let rhs = forcing[i] / scales[i];
            if !rhs.is_finite() || rhs == 0.0 {
                return Err(invalid("scaled electric force is unrepresentable"));
            }
            velocity.push(rhs);
            for j in 0..i {
                let entry = matrix[i * n + j] / scales[i] / scales[j];
                if !entry.is_finite() || entry == 0.0 {
                    return Err(invalid("scaled Coulomb pair friction is unrepresentable"));
                }
                matrix[i * n + j] = entry;
                matrix[j * n + i] = entry;
            }
            matrix[i * n + i] = 1.0;
        }
        poll(&mut cancelled)?;
        let factor = fs_la::factor::cholesky(&matrix, n)
            .map_err(|_| invalid("anchored Coulomb friction lost positive definiteness"))?;
        factor.solve(&mut velocity);
        poll(&mut cancelled)?;
        for (v, scale) in velocity.iter_mut().zip(&scales) {
            *v /= scale;
            if !v.is_finite() {
                return Err(invalid("Coulomb-coupled velocity is unrepresentable"));
            }
        }

        let mut force = vec![0.0; n];
        let mut gross_force = vec![0.0; n];
        let mut dissipation = Sum::default();
        let mut conductivity = Sum::default();
        let mut thermal_coefficients = [Sum::default(); 3];
        let mut maximum_field = f64::INFINITY;
        for i in 0..n {
            poll(&mut cancelled)?;
            force[i] = background_drag[i] * velocity[i];
            gross_force[i] = forcing[i].abs() + force[i].abs();
            dissipation.add(force[i] * velocity[i]);
            thermal_coefficients[usize::from(self.carriers[i].is_electron)]
                .add(force[i] * velocity[i]);
            conductivity.add(forcing[i] * velocity[i]);
            let carrier = &mut self.carriers[i];
            carrier.drift_mobility_m2_v_s = velocity[i];
            carrier.conductivity_s_m = forcing[i] * velocity[i];
            carrier.maximum_electric_field_v_m = if velocity[i] == 0.0 {
                f64::INFINITY
            } else {
                sqrt(super::BOLTZMANN_J_K * carrier.temperature_k / carrier.mass_kg)
                    * self.limits.maximum_drift_thermal_ratio
                    / velocity[i].abs()
            };
            maximum_field = maximum_field.min(carrier.maximum_electric_field_v_m);
        }
        for (pair, &[i, j]) in pairs.iter().zip(&pair_carriers) {
            poll(&mut cancelled)?;
            let relative = velocity[i] - velocity[j];
            let drag = pair.friction_coefficient_kg_m3_s() * relative;
            force[i] += drag;
            force[j] -= drag;
            gross_force[i] += drag.abs();
            gross_force[j] += drag.abs();
            dissipation.add(drag * relative);
            if relative != 0.0 {
                // The response is linear at this frozen state, but a 1 V/m
                // computational basis need not lie inside its physical domain.
                // Resolve the source's thermal fractions at an admitted probe
                // speed, then cache their coefficients for the actual field.
                let probe = pair.maximum_relative_speed_squared_m2_s2().min(1.0);
                let heating = pair.drift_heating(probe).map_err(source_error)?;
                let total = heating.total_power_density_w_m3();
                if !positive(total) {
                    return Err(invalid("Coulomb thermal fraction is unrepresentable"));
                }
                let dissipated = drag * relative;
                thermal_coefficients[usize::from(self.carriers[i].is_electron)]
                    .add(dissipated * (heating.first_friction_power_density_w_m3() / total));
                thermal_coefficients[usize::from(self.carriers[j].is_electron)]
                    .add(dissipated * (heating.second_friction_power_density_w_m3() / total));
                let correction =
                    dissipated * (heating.thermal_exchange_correction_to_first_w_m3() / total);
                if self.carriers[i].is_electron {
                    thermal_coefficients[2].add(correction);
                } else if self.carriers[j].is_electron {
                    thermal_coefficients[2].add(-correction);
                }
                let bound = sqrt(pair.maximum_relative_speed_squared_m2_s2()) / relative.abs();
                // Leave room for the actual vector-field arithmetic to reach
                // the source's closed relative-drift boundary without overshoot.
                maximum_field = maximum_field.min(bound * (1.0 - 16.0 * f64::EPSILON));
            }
        }
        let mut force_residual = 0.0_f64;
        for i in 0..n {
            if !force[i].is_finite() || !positive(gross_force[i]) {
                return Err(invalid("original Coulomb force scale is unrepresentable"));
            }
            let residual = (force[i] - forcing[i]).abs() / gross_force[i];
            if !residual.is_finite() || residual > BALANCE_RELATIVE_TOLERANCE {
                return Err(invalid("original Coulomb momentum balance failed"));
            }
            force_residual = force_residual.max(residual);
        }
        let power_residual = relative_power_error(conductivity.value, dissipation.value)?;
        let thermal_coefficients = thermal_coefficients.map(|sum| sum.value);
        if thermal_coefficients.iter().any(|value| !value.is_finite())
            || thermal_coefficients[0] < 0.0
            || thermal_coefficients[1] < 0.0
        {
            return Err(invalid("Coulomb thermal response is unrepresentable"));
        }
        relative_power_error(
            conductivity.value,
            thermal_coefficients[0] + thermal_coefficients[1],
        )?;
        if !positive(conductivity.value) || !positive(maximum_field) {
            return Err(invalid(
                "Coulomb conductivity or weak-drift field is unrepresentable",
            ));
        }
        self.conductivity_s_m = conductivity.value;
        self.maximum_electric_field_v_m = maximum_field;
        self.coulomb = Some(CoulombDrag {
            model,
            pairs,
            thermal_coefficients,
            force_balance_relative_residual: force_residual,
            power_balance_relative_residual: power_residual,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original source model; its frozen populations also bind thermal evolution.
    pub fn coulomb_model(&self) -> Option<&CoulombModelV1> {
        self.coulomb.as_ref().map(|drag| &drag.model)
    }
    /// Complete pair moments at this snapshot's actual two temperatures.
    pub fn coulomb_pairs(&self) -> &[CoulombPairV1] {
        self.coulomb.as_ref().map_or(&[], |drag| &drag.pairs)
    }
    /// Largest original per-carrier force residual, normalized by gross force.
    pub fn coulomb_force_balance_relative_residual(&self) -> Option<f64> {
        self.coulomb
            .as_ref()
            .map(|drag| drag.force_balance_relative_residual)
    }
    /// Relative difference between summed electrical work and drag dissipation.
    pub fn coulomb_power_balance_relative_residual(&self) -> Option<f64> {
        self.coulomb
            .as_ref()
            .map(|drag| drag.power_balance_relative_residual)
    }

    pub(super) fn coulomb_joule_heating_w_m3(
        &self,
        field: [f64; 2],
    ) -> Result<IonizedGasJouleHeating, IonizedGasError> {
        let drag = self.coulomb.as_ref().ok_or(IonizedGasError::Invalid(
            "Coulomb heating requires its original coupled response",
        ))?;
        let [x, y] = field;
        let [heavy_friction, electron_friction, correction] = drag
            .thermal_coefficients
            .map(|coefficient| x * (coefficient * x) + y * (coefficient * y));
        let heating = IonizedGasJouleHeating {
            electron_power_density_w_m3: electron_friction + correction,
            heavy_power_density_w_m3: heavy_friction - correction,
            electron_friction_power_density_w_m3: electron_friction,
            heavy_friction_power_density_w_m3: heavy_friction,
            thermal_exchange_correction_to_electrons_w_m3: correction,
        };
        let electrical = x * (self.conductivity_s_m * x) + y * (self.conductivity_s_m * y);
        if !heating.electron_power_density_w_m3.is_finite()
            || !heating.heavy_power_density_w_m3.is_finite()
            || (field != [0.0; 2] && electrical == 0.0)
        {
            return Err(IonizedGasError::Invalid(
                "Coulomb thermal powers are unrepresentable",
            ));
        }
        relative_power_error(electrical, heating.total_power_density_w_m3())?;
        Ok(heating)
    }
}

fn source_error(error: CoulombErrorV1) -> IonizedGasError {
    match error {
        CoulombErrorV1::Cancelled => IonizedGasError::Cancelled,
        error => IonizedGasError::Coulomb(error),
    }
}

fn relative_power_error(work: f64, heat: f64) -> Result<f64, IonizedGasError> {
    let residual = if work == 0.0 && heat == 0.0 {
        0.0
    } else {
        (work - heat).abs() / work.abs().max(heat.abs())
    };
    if !work.is_finite()
        || !heat.is_finite()
        || work < 0.0
        || heat < 0.0
        || !residual.is_finite()
        || residual > BALANCE_RELATIVE_TOLERANCE
    {
        return Err(IonizedGasError::Invalid(
            "electrical work and Coulomb friction heat disagree",
        ));
    }
    Ok(residual)
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.correction;
        let next = self.value + increment;
        self.correction = (next - self.value) - increment;
        self.value = next;
    }
}
