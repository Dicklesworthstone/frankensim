use super::*;
mod collisions;
use collisions::collision_error;

const CV: f64 = 1.5 * R;

#[derive(Clone, Copy)]
pub(super) struct Point {
    pub time: f64,
    pub extent: f64,
    pub heavy_temperature: f64,
    pub electron_temperature: f64,
    pub exchange: f64,
    pub coulomb_exchange: f64,
    pub electron_temperature_extent_derivative: f64,
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    channel: &'a AtomicIonizationV1,
    powers: [f64; 2],
    controls: AtomicIonizationSolve,
    amounts: [f64; 3],
    indices: [usize; 3],
    electron_sensible: f64,
    heavy_limits: [f64; 2],
    exchange: Option<PreparedExchange<'a>>,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        channel: &'a AtomicIonizationV1,
        collisions: Option<&'a ElectronHeavyCollisionSources>,
        powers: [f64; 2],
        controls: AtomicIonizationSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AtomicIonizationError> {
        if initial.thermodynamics() != channel.thermodynamics()
            || initial.electron_index() != Some(channel.electron_index())
            || initial.relative_charge_residual() > controls.relative_charge_tolerance
        {
            return Err(AtomicIonizationError::Invalid(
                "state and original channel bank disagree",
            ));
        }
        let indices = [
            channel.neutral_index(),
            channel.ion_index(),
            channel.electron_index(),
        ];
        let amounts = indices.map(|i| initial.amounts()[i].value());
        if amounts[2] == 0.0 && powers[1] != 0.0 {
            return Err(AtomicIonizationError::Invalid(
                "electron power without electrons",
            ));
        }
        for (i, model) in initial.thermodynamics().models().iter().enumerate() {
            if initial.amounts()[i].value() == 0.0 && !indices.contains(&i) {
                continue;
            }
            let electron = i == indices[2];
            let (lo, hi, t) = if electron {
                (
                    controls.minimum_electron_temperature_k,
                    controls.maximum_electron_temperature_k,
                    initial.electron_temperature().value(),
                )
            } else {
                (
                    controls.minimum_heavy_temperature_k,
                    controls.maximum_heavy_temperature_k,
                    initial.heavy_temperature().value(),
                )
            };
            if !(lo..=hi).contains(&t) {
                return Err(AtomicIonizationError::OutsideDomain);
            }
            admit_caloric_interval(
                model,
                Temperature::new(lo),
                Temperature::new(hi),
                &mut *cancelled,
            )
            .map_err(kinetic_error)?;
        }
        let mut result = Self {
            initial,
            channel,
            powers,
            controls,
            amounts,
            indices,
            electron_sensible: CV * amounts[2] * initial.electron_temperature().value(),
            heavy_limits: [0.0; 2],
            exchange: collisions
                .map(|sources| PreparedExchange::new(initial, sources, indices[0], cancelled))
                .transpose()
                .map_err(collision_error)?,
        };
        result.heavy_limits = [
            result
                .heavy_calorics(controls.minimum_heavy_temperature_k, cancelled)?
                .0,
            result
                .heavy_calorics(controls.maximum_heavy_temperature_k, cancelled)?
                .0,
        ];
        if !result.electron_sensible.is_finite() {
            return Err(AtomicIonizationError::Arithmetic(
                "initial electron sensible energy",
            ));
        }
        // Admit the initial collision/temperature regime even for a zero-duration call.
        channel
            .evaluate(
                initial.heavy_temperature(),
                initial.electron_temperature(),
                amounts.map(|n| n / initial.volume().value()),
                &mut *cancelled,
            )
            .map_err(channel_error)?;
        Ok(result)
    }

    fn heavy_calorics(
        &self,
        t: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), AtomicIonizationError> {
        super::super::energy::source_pool_calorics(self.initial, false, t, cancelled).map_err(
            |error| match error {
                super::super::energy::ElectronHeavyEnergyError::Cancelled => {
                    AtomicIonizationError::Cancelled
                }
                error => AtomicIonizationError::Caloric(error),
            },
        )
    }

    fn heavy_temperature(
        &self,
        time: f64,
        guess: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AtomicIonizationError> {
        let target = self.powers[0] * time;
        if !target.is_finite() {
            return Err(AtomicIonizationError::Arithmetic("heavy power integral"));
        }
        if target == 0.0 {
            return Ok(self.initial.heavy_temperature().value());
        }
        if !(self.heavy_limits[0]..=self.heavy_limits[1]).contains(&target) {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let mut lo = self.controls.minimum_heavy_temperature_k;
        let mut hi = self.controls.maximum_heavy_temperature_k;
        let mut t = guess;
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            let (energy, cv) = self.heavy_calorics(t, cancelled)?;
            let residual = energy - target;
            let tolerance = (self.controls.energy_tolerance_j / 32.0).min(
                cv * (self.controls.temperature_tolerance_k + self.controls.relative_tolerance * t)
                    / 128.0,
            );
            if residual.abs() <= tolerance {
                return Ok(t);
            }
            if residual > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
            let next = t - residual / cv;
            t = if next > lo && next < hi {
                next
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(AtomicIonizationError::Nonconvergence)
    }

    fn amounts_at(&self, extent: f64) -> [f64; 3] {
        [
            self.amounts[0] - extent,
            self.amounts[1] + extent,
            self.amounts[2] + extent,
        ]
    }

    fn point(&self, time: f64, extent: f64, heavy: f64) -> Result<Point, AtomicIonizationError> {
        let amounts = self.amounts_at(extent);
        if amounts.iter().any(|n| !n.is_finite() || *n < 0.0) {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let te = self.electron_temperature(time, extent);
        if !te.is_finite()
            || !(self.controls.minimum_electron_temperature_k
                ..=self.controls.maximum_electron_temperature_k)
                .contains(&te)
        {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        Ok(Point {
            time,
            extent,
            heavy_temperature: heavy,
            electron_temperature: te,
            exchange: 0.0,
            coulomb_exchange: 0.0,
            electron_temperature_extent_derivative: if amounts[2] == 0.0 {
                0.0
            } else {
                -(self.channel.ionization_energy_j_mol() + CV * te) / (CV * amounts[2])
            },
        })
    }

    fn electron_temperature(&self, time: f64, extent: f64) -> f64 {
        let initial = self.initial.electron_temperature().value();
        let amount = self.amounts[2] + extent;
        if amount == 0.0 && self.amounts[2] == 0.0 && extent == 0.0 {
            return initial;
        }
        let work = self.powers[1] * time;
        let chi = self.channel.ionization_energy_j_mol();
        // Preserve a small physical increment, including the exact unchanged
        // initial temperature at a closed source endpoint. For large cooling,
        // the direct energy quotient avoids cancellation against initial T.
        let increment = (work - (chi + CV * initial) * extent) / (CV * amount);
        if increment.is_finite() && increment.abs() <= 0.5 * initial {
            initial + increment
        } else {
            (self.electron_sensible + work - chi * extent) / (CV * amount)
        }
    }

    fn bracket_point(
        &self,
        time: f64,
        mut extent: f64,
        other: f64,
        heavy: f64,
        temperatures: [f64; 2],
    ) -> Result<Point, AtomicIonizationError> {
        let lower = extent < other;
        for _ in 0..16 {
            let te = self.electron_temperature(time, extent);
            if (temperatures[0]..=temperatures[1]).contains(&te) {
                return self.point(time, extent, heavy);
            }
            if !te.is_finite()
                || extent == other
                || (lower && te < temperatures[0])
                || (!lower && te > temperatures[1])
            {
                return Err(AtomicIonizationError::OutsideDomain);
            }
            // Only move an endpoint if reconstructing its actual energy state
            // falls outside the admitted temperature interval. An extent ULP
            // alone is insufficient near zero: scale by the electron inventory
            // and the inverse energy derivative, without clipping temperature.
            let boundary = temperatures[usize::from(lower)];
            let correction = (CV * (self.amounts[2] + extent)
                / (self.channel.ionization_energy_j_mol() + CV * boundary))
                * (te - boundary).abs();
            let roundoff = 8.0 * f64::EPSILON * self.amounts[2].max(extent.abs());
            let shift = correction.max(roundoff);
            let candidate = if lower {
                (extent + shift).max(extent.next_up())
            } else {
                (extent - shift).min(extent.next_down())
            };
            let next = if candidate > extent.min(other) && candidate < extent.max(other) {
                candidate
            } else {
                extent + 0.5 * (other - extent)
            };
            if next == extent || !next.is_finite() {
                return Err(AtomicIonizationError::OutsideDomain);
            }
            extent = next;
        }
        Err(AtomicIonizationError::OutsideDomain)
    }

    fn cheap_bracket_sign(
        &self,
        point: &Point,
        start: &Point,
        lower: bool,
    ) -> Result<bool, AtomicIonizationError> {
        if (lower && point.extent > start.extent) || (!lower && point.extent < start.extent) {
            return Ok(false);
        }
        let concentrations = self
            .amounts_at(point.extent)
            .map(|n| n / self.initial.volume().value());
        let [neutral, ion, electron] = concentrations;
        if electron == 0.0 || (neutral == 0.0 && ion == 0.0) {
            return Ok(true);
        }
        let affinity = if neutral == 0.0 {
            f64::NEG_INFINITY
        } else if ion == 0.0 {
            f64::INFINITY
        } else {
            let log_kc = self
                .channel
                .log_equilibrium_constant(Temperature::new(point.electron_temperature))
                .map_err(channel_error)?;
            // Use the same concentration normalization and operation order as
            // the source channel's signed affinity, including near balance.
            log_kc + fs_math::det::ln(neutral) - fs_math::det::ln(ion) - fs_math::det::ln(electron)
        };
        // Detailed balance certifies q's sign without Maxwellian quadrature.
        // If both terms of xi-xi_start-dt*V*q have the needed sign, its
        // magnitude cannot change the bracket admission.
        Ok(if lower {
            affinity >= 0.0
        } else {
            affinity <= 0.0
        })
    }

    fn residual(
        &self,
        point: &Point,
        start: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), AtomicIonizationError> {
        let volume = self.initial.volume().value();
        let amounts = self.amounts_at(point.extent);
        let evaluated = self
            .channel
            .evaluate(
                Temperature::new(point.heavy_temperature),
                Temperature::new(point.electron_temperature),
                amounts.map(|n| n / volume),
                &mut *cancelled,
            )
            .map_err(channel_error)?;
        let dt = point.time - start.time;
        let residual =
            point.extent - start.extent - dt * (volume * evaluated.progress_rate_mol_m3_s);
        let dc = evaluated.concentration_derivatives_per_s;
        let temperature_derivative = point.electron_temperature_extent_derivative;
        let derivative = 1.0
            - dt * (-dc[0]
                + dc[1]
                + dc[2]
                + volume
                    * evaluated.electron_temperature_derivative_mol_m3_s_k
                    * temperature_derivative);
        if !residual.is_finite() || !derivative.is_finite() {
            return Err(AtomicIonizationError::Arithmetic(
                "constrained extent residual or derivative",
            ));
        }
        Ok((residual, derivative))
    }

    pub fn step(
        &self,
        start: &Point,
        end: f64,
        coupling: ExchangeConductance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, AtomicIonizationError> {
        if self.exchange.is_some() && self.amounts[2] > 0.0 {
            return self.step_with_exchange(start, end, coupling, cancelled);
        }
        let th = self.heavy_temperature(end, start.heavy_temperature, cancelled)?;
        if self.amounts[2] == 0.0 {
            // Electron impact and its three-body reverse both require electrons.
            return self.point(end, 0.0, th);
        }
        let chi = self.channel.ionization_energy_j_mol();
        let energy = self.electron_sensible + self.powers[1] * end;
        if !energy.is_finite() {
            return Err(AtomicIonizationError::Arithmetic("electron power integral"));
        }
        let initial_te = self.initial.electron_temperature().value();
        let lower_energy =
            (CV * self.amounts[2] * (initial_te - self.controls.maximum_electron_temperature_k)
                + self.powers[1] * end)
                / (chi + CV * self.controls.maximum_electron_temperature_k);
        // Intersect with the actual stationary-target domain at this heavy T.
        // An inadmissibly cold bracket endpoint must not reject a supported
        // interior solution.
        let mut minimum_te = self.controls.minimum_electron_temperature_k;
        for &i in &self.indices[..2] {
            let mass = self.initial.thermodynamics().models()[i]
                .molar_mass()
                .value()
                / fs_thermochem::electron_collision::AVOGADRO_PER_MOL;
            minimum_te = minimum_te.max(
                (th * (fs_thermochem::electron_collision::ELECTRON_MASS_KG / mass)
                    / self
                        .channel
                        .forward_rate()
                        .maximum_target_thermal_speed_ratio()
                        .powi(2))
                .next_up(),
            );
        }
        if minimum_te > self.controls.maximum_electron_temperature_k {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let upper_energy = (CV * self.amounts[2] * (initial_te - minimum_te)
            + self.powers[1] * end)
            / (chi + CV * minimum_te);
        // Keep closed physical endpoints, including an unchanged state exactly
        // at its source limit. Repair only actual reconstruction roundoff below.
        let mut lo = (-self.amounts[1])
            .max((-self.amounts[2]).next_up())
            .max(lower_energy);
        let mut hi = self.amounts[0].min(upper_energy);
        if lo > hi || !lo.is_finite() || !hi.is_finite() {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let temperatures = [minimum_te, self.controls.maximum_electron_temperature_k];
        let lower = self.bracket_point(end, lo, hi, th, temperatures)?;
        lo = lower.extent;
        let upper = self.bracket_point(end, hi, lo, th, temperatures)?;
        hi = upper.extent;
        if (!self.cheap_bracket_sign(&lower, start, true)?
            && self.residual(&lower, start, cancelled)?.0 > 0.0)
            || (!self.cheap_bracket_sign(&upper, start, false)?
                && self.residual(&upper, start, cancelled)?.0 < 0.0)
        {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let mut extent = if start.extent >= lo && start.extent <= hi {
            start.extent
        } else {
            lo + 0.5 * (hi - lo)
        };
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            let point = self.point(end, extent, th)?;
            let (residual, derivative) = self.residual(&point, start, cancelled)?;
            let n = self.amounts_at(extent);
            let dte = (chi + CV * point.electron_temperature) / (CV * n[2]);
            let amount_scale = self.controls.absolute_amount_tolerance_mol
                + self.controls.relative_tolerance
                    * n.iter().copied().fold(f64::INFINITY, f64::min);
            let temperature_scale = self.controls.temperature_tolerance_k
                + self.controls.relative_tolerance * point.electron_temperature;
            let tolerance = (amount_scale / 128.0).min(temperature_scale / (128.0 * dte));
            // For stiff channels, residuals naturally amplify tiny extent errors.
            // Require a small Newton correction as well as a bounded bracket;
            // the derivative magnitude converts the residual back to mol.
            if residual.abs() <= tolerance * derivative.abs().max(1.0) {
                return Ok(point);
            }
            if residual > 0.0 {
                hi = extent;
            } else {
                lo = extent;
            }
            let next = extent - residual / derivative;
            let candidate = if derivative > 0.0 && next > lo && next < hi {
                next
            } else {
                lo + 0.5 * (hi - lo)
            };
            if candidate == extent {
                return Err(AtomicIonizationError::Nonconvergence);
            }
            extent = candidate;
        }
        Err(AtomicIonizationError::Nonconvergence)
    }

    pub fn error(&self, full: &Point, half: &Point) -> f64 {
        let a = self.amounts_at(full.extent);
        let b = self.amounts_at(half.extent);
        let mut error = 0.0_f64;
        for (a, b) in a.into_iter().zip(b) {
            error = error.max(
                (a - b).abs()
                    / (self.controls.absolute_amount_tolerance_mol
                        + self.controls.relative_tolerance * a.abs().max(b.abs())),
            );
        }
        error = error.max(
            (full.electron_temperature - half.electron_temperature).abs()
                / (self.controls.temperature_tolerance_k
                    + self.controls.relative_tolerance
                        * full
                            .electron_temperature
                            .abs()
                            .max(half.electron_temperature.abs())),
        );
        error = error.max(
            (full.heavy_temperature - half.heavy_temperature).abs()
                / (self.controls.temperature_tolerance_k
                    + self.controls.relative_tolerance
                        * full
                            .heavy_temperature
                            .abs()
                            .max(half.heavy_temperature.abs())),
        );
        for (a, b) in [
            (full.exchange, half.exchange),
            (full.coulomb_exchange, half.coulomb_exchange),
            (
                full.exchange - full.coulomb_exchange,
                half.exchange - half.coulomb_exchange,
            ),
        ] {
            error = error.max(
                (a - b).abs()
                    / (self.controls.energy_tolerance_j
                        + self.controls.relative_tolerance * a.abs().max(b.abs())),
            );
        }
        error
    }

    pub fn entropy_supply(&self, dt: f64, point: &Point) -> f64 {
        dt * (self.powers[0] / point.heavy_temperature
            + self.powers[1] / point.electron_temperature)
    }

    fn amount(&self, index: usize, extent: f64) -> f64 {
        let sign = if index == self.indices[0] {
            -1.0
        } else if index == self.indices[1] || index == self.indices[2] {
            1.0
        } else {
            0.0
        };
        self.initial.amounts()[index].value() + sign * extent
    }

    pub fn entropy_change(
        &self,
        from: &Point,
        to: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AtomicIonizationError> {
        let mut sum = 0.0;
        let mut correction = 0.0;
        let log_volume = fs_math::det::ln(self.initial.volume().value());
        for (i, model) in self.initial.thermodynamics().models().iter().enumerate() {
            let n0 = self.amount(i, from.extent);
            let n1 = self.amount(i, to.extent);
            if n0 == 0.0 && n1 == 0.0 {
                continue;
            }
            poll(cancelled)?;
            let (t0, t1) = if i == self.indices[2] {
                (from.electron_temperature, to.electron_temperature)
            } else {
                (from.heavy_temperature, to.heavy_temperature)
            };
            let source_error = |source| AtomicIonizationError::Species {
                species: model.species().clone(),
                source,
            };
            let thermal = model
                .molar_entropy_change(Temperature::new(t0), Temperature::new(t1))
                .map_err(source_error)?
                .value()
                - R * log_ratio(t1, t0);
            let dn = n1 - n0;
            let composition = if dn == 0.0 {
                0.0
            } else {
                let reference = model
                    .evaluate_properties(Temperature::new(t0))
                    .map_err(source_error)?
                    .s()
                    .value()
                    - R * (fs_math::det::ln(R) + fs_math::det::ln(t0)
                        - log_volume
                        - fs_math::det::ln(model.convention().reference_pressure().value()));
                let mixing = if n0 > 0.0 && n1 > 0.0 {
                    dn * fs_math::det::ln(n0) + n1 * log_ratio(n1, n0)
                } else if n1 > 0.0 {
                    n1 * fs_math::det::ln(n1)
                } else {
                    -n0 * fs_math::det::ln(n0)
                };
                dn * reference - R * mixing
            };
            let increment = n1 * thermal + composition - correction;
            let next = sum + increment;
            correction = (next - sum) - increment;
            sum = next;
        }
        if !sum.is_finite() {
            return Err(AtomicIonizationError::Arithmetic(
                "two-pool mixture entropy",
            ));
        }
        Ok(sum)
    }

    pub fn publish(
        &self,
        point: &Point,
        balance: &mut AtomicIonizationBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, AtomicIonizationError> {
        let state = self.state_at(point, cancelled)?;
        let heavy_change = self.heavy_calorics(point.heavy_temperature, cancelled)?.0;
        let electron = &self.initial.thermodynamics().models()[self.indices[2]];
        let delta_u = electron
            .molar_internal_energy_change(
                self.initial.electron_temperature(),
                state.electron_temperature(),
            )
            .map_err(|source| AtomicIonizationError::Species {
                species: electron.species().clone(),
                source,
            })?
            .value();
        let electron_change =
            self.amounts[2] * delta_u + point.extent * CV * point.electron_temperature;
        balance.electron_deposited_energy_j = self.powers[1] * point.time;
        balance.heavy_deposited_energy_j = self.powers[0] * point.time;
        balance.chemical_storage_change_j = self.channel.ionization_energy_j_mol() * point.extent;
        balance.electron_to_heavy_energy_j = point.exchange;
        balance.electron_to_heavy_coulomb_energy_j = point.coulomb_exchange;
        balance.electron_to_heavy_neutral_energy_j = point.exchange - point.coulomb_exchange;
        balance.electron_sensible_energy_residual_j =
            electron_change + balance.chemical_storage_change_j + point.exchange
                - balance.electron_deposited_energy_j;
        balance.heavy_sensible_energy_residual_j =
            heavy_change - point.exchange - balance.heavy_deposited_energy_j;
        // Reference-invariant original-source increment, not subtraction of
        // two large absolute parcel energies containing arbitrary formation zeros.
        balance.energy_residual_j =
            balance.electron_sensible_energy_residual_j + balance.heavy_sensible_energy_residual_j;
        if [
            balance.electron_sensible_energy_residual_j,
            balance.heavy_sensible_energy_residual_j,
            balance.energy_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > self.controls.energy_tolerance_j)
        {
            return Err(AtomicIonizationError::Balance("source energy budget"));
        }
        Ok(state)
    }
}

// Stable ln(a/b) using the same deterministic Kahan log1p identity as source
// entropy differences. No ratio overflow or loss of a tiny temperature change.
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    sign * if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        fs_math::det::ln(hi) - fs_math::det::ln(lo)
    }
}
