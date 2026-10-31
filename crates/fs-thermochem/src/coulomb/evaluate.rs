use super::*;

impl CoulombModelV1 {
    pub(super) fn evaluate_inner(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        mut visitor: impl FnMut(CoulombPairV1),
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoulombSummaryV1, CoulombErrorV1> {
        poll(&mut cancelled)?;
        let (th, te) = (heavy_temperature.value(), electron_temperature.value());
        if !positive(th) || !positive(te) {
            return Err(CoulombErrorV1::Invalid(
                "positive finite actual heavy/electron temperatures",
            ));
        }
        // Preserve original source support, including the neutral populations
        // whose shared heavy temperature is part of this exact state binding.
        for (index, (model, amount)) in self
            .initial_state
            .thermodynamics()
            .models()
            .iter()
            .zip(self.initial_state.amounts())
            .enumerate()
        {
            poll(&mut cancelled)?;
            if amount.value() == 0.0 {
                continue;
            }
            let temperature = if self.initial_state.electron_index() == Some(index) {
                electron_temperature
            } else {
                heavy_temperature
            };
            model
                .evaluate_properties(temperature)
                .map_err(|source| CoulombErrorV1::Species {
                    species: model.species().clone(),
                    source,
                })?;
        }
        let mut summary = CoulombSummaryV1 {
            pair_count: 0,
            electron_heavy_thermal: 0.0,
            debye_length: None,
            minimum_logarithm: None,
            debye_occupation: None,
            maximum_coupling: 0.0,
            maximum_phase_space: 0.0,
            maximum_relativistic: 0.0,
            maximum_quantum: 0.0,
        };
        if self.carriers.is_empty() {
            poll(&mut cancelled)?;
            return Ok(summary);
        }
        let log_kb = ln(BOLTZMANN_J_PER_K);
        let log_eps = ln(VACUUM_PERMITTIVITY_F_M);
        let log_h = ln(PLANCK_J_S);
        let log_2pi = ln(2.0 * core::f64::consts::PI);
        let log_4pi = ln(4.0 * core::f64::consts::PI);
        let log_th = ln(th);
        let log_te = ln(te);
        let temperature_log = |carrier: &Carrier| if carrier.electron { log_te } else { log_th };
        let mut log_screening = f64::NEG_INFINITY;
        let mut log_density = f64::NEG_INFINITY;
        let mut largest_log_charge = f64::NEG_INFINITY;
        let mut coldest_log_temperature = f64::INFINITY;
        for carrier in self.carriers.iter() {
            poll(&mut cancelled)?;
            let log_t = temperature_log(carrier);
            log_density = log_add(log_density, carrier.log_density);
            largest_log_charge = largest_log_charge.max(carrier.log_charge);
            coldest_log_temperature = coldest_log_temperature.min(log_t);
            if self.screening == CoulombScreeningV1::StaticAllChargedDebye || carrier.electron {
                log_screening = log_add(
                    log_screening,
                    carrier.log_density + 2.0 * carrier.log_charge - log_eps - log_kb - log_t,
                );
            }
            let log_phase = carrier.log_density + 3.0 * log_h
                - 1.5 * (log_2pi + carrier.log_mass + log_kb + log_t);
            let log_relativistic = log_kb + log_t - carrier.log_mass - 2.0 * ln(SPEED_OF_LIGHT_M_S);
            if log_phase > ln(self.limits.maximum_thermal_phase_space_density) {
                return Err(CoulombErrorV1::Domain(
                    "charged species exceeds classical phase-space limit",
                ));
            }
            if log_relativistic > ln(self.limits.maximum_thermal_energy_rest_energy_ratio) {
                return Err(CoulombErrorV1::Domain(
                    "charged species exceeds nonrelativistic thermal-energy limit",
                ));
            }
            summary.maximum_phase_space = summary.maximum_phase_space.max(exp(log_phase));
            summary.maximum_relativistic = summary.maximum_relativistic.max(exp(log_relativistic));
        }
        if !log_screening.is_finite() {
            return Err(CoulombErrorV1::Arithmetic("Debye screening population"));
        }
        let log_debye = -0.5 * log_screening;
        let log_occupation = log_4pi - ln(3.0) + log_density + 3.0 * log_debye;
        if log_occupation < ln(self.limits.minimum_particles_per_debye_sphere) {
            return Err(CoulombErrorV1::Domain(
                "too few charged particles per Debye sphere",
            ));
        }
        let log_coupling = 2.0 * largest_log_charge - log_4pi - log_eps
            + (log_4pi - ln(3.0) + log_density) / 3.0
            - log_kb
            - coldest_log_temperature;
        if log_coupling > ln(self.limits.maximum_coulomb_coupling) {
            return Err(CoulombErrorV1::Domain(
                "charged state exceeds weak Coulomb-coupling limit",
            ));
        }
        summary.debye_length = Some(from_log(log_debye, "Debye length")?);
        summary.debye_occupation = Some(from_log(log_occupation, "Debye occupation")?);
        summary.maximum_coupling = exp(log_coupling);
        let mut thermal_sum = 0.0_f64;
        let mut thermal_correction = 0.0_f64;
        let mut minimum_logarithm = f64::INFINITY;
        for (first_slot, first) in self.carriers.iter().enumerate() {
            for second in &self.carriers[first_slot + 1..] {
                poll(&mut cancelled)?;
                let log_first_variance = log_kb + temperature_log(first) - first.log_mass;
                let log_second_variance = log_kb + temperature_log(second) - second.log_mass;
                let log_variance = log_add(log_first_variance, log_second_variance);
                let variance = from_log(log_variance, "relative thermal variance")?;
                let first_variance = from_log(log_first_variance, "first thermal variance")?;
                let second_variance = from_log(log_second_variance, "second thermal variance")?;
                let mass_ratio = first.mass.min(second.mass) / first.mass.max(second.mass);
                let log_mass_sum = first.log_mass.max(second.log_mass) + ln(1.0 + mass_ratio);
                let log_reduced_mass = first.log_mass.min(second.log_mass) - ln(1.0 + mass_ratio);
                let first_mass_fraction = if first.mass <= second.mass {
                    mass_ratio / (1.0 + mass_ratio)
                } else {
                    1.0 / (1.0 + mass_ratio)
                };
                let second_mass_fraction = if first.mass <= second.mass {
                    1.0 / (1.0 + mass_ratio)
                } else {
                    mass_ratio / (1.0 + mass_ratio)
                };
                let log_c = first.log_charge + second.log_charge - log_4pi - log_eps;
                let log_b90 = log_c - log_reduced_mass - log_variance;
                let logarithm = log_debye - log_b90;
                if !logarithm.is_finite() || logarithm < self.limits.minimum_coulomb_logarithm {
                    return Err(CoulombErrorV1::Domain(
                        "pair Coulomb logarithm below declared leading-log limit",
                    ));
                }
                // bq/b90 = h*sqrt(S)/(4*pi*C); reduced mass cancels.
                let log_quantum = log_h - log_4pi + 0.5 * log_variance - log_c;
                if log_quantum > ln(self.limits.maximum_quantum_to_classical_ratio) {
                    return Err(CoulombErrorV1::Domain(
                        "quantum diffraction exceeds classical pair-scattering limit",
                    ));
                }
                summary.maximum_quantum = summary.maximum_quantum.max(exp(log_quantum));
                let log_k = ln(4.0 * sqrt(2.0 * core::f64::consts::PI) / 3.0)
                    + first.log_density
                    + second.log_density
                    + 2.0 * log_c
                    + ln(logarithm)
                    - log_reduced_mass
                    - 1.5 * log_variance;
                let friction = from_log(log_k, "pair Coulomb momentum coefficient")?;
                let thermal = from_log(
                    ln(3.0) + log_kb + log_k - log_mass_sum,
                    "pair Coulomb thermal conductance",
                )?;
                let first_t = if first.electron { te } else { th };
                let second_t = if second.electron { te } else { th };
                let beta = if first_t >= second_t {
                    first_mass_fraction
                        * (first_variance / variance)
                        * ((first_t - second_t) / first_t)
                } else {
                    -second_mass_fraction
                        * (second_variance / variance)
                        * ((second_t - first_t) / second_t)
                };
                let maximum_speed = self.limits.maximum_relative_drift_thermal_ratio
                    * sqrt(first_variance.min(second_variance));
                let maximum_speed_squared = maximum_speed * maximum_speed;
                if !positive(maximum_speed_squared) || !beta.is_finite() {
                    return Err(CoulombErrorV1::Arithmetic(
                        "pair drift bound or thermal correction",
                    ));
                }
                let pair = CoulombPairV1 {
                    first_index: first.index,
                    second_index: second.index,
                    friction,
                    thermal,
                    coulomb_logarithm: logarithm,
                    relative_variance: variance,
                    maximum_speed_squared,
                    first_mass_fraction,
                    second_mass_fraction,
                    thermal_drift_coefficient: 1.5 * beta,
                };
                if first.electron || second.electron {
                    let sum = thermal_sum + thermal;
                    thermal_correction += if thermal_sum.abs() >= thermal.abs() {
                        (thermal_sum - sum) + thermal
                    } else {
                        (thermal - sum) + thermal_sum
                    };
                    thermal_sum = sum;
                }
                minimum_logarithm = minimum_logarithm.min(logarithm);
                summary.pair_count += 1;
                visitor(pair);
            }
        }
        summary.electron_heavy_thermal = thermal_sum + thermal_correction;
        if !summary.electron_heavy_thermal.is_finite() {
            return Err(CoulombErrorV1::Arithmetic(
                "total electron-heavy thermal conductance",
            ));
        }
        summary.minimum_logarithm = (summary.pair_count > 0).then_some(minimum_logarithm);
        poll(&mut cancelled)?;
        Ok(summary)
    }
}
