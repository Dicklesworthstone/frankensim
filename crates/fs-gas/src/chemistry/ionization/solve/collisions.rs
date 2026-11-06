//! Simultaneous reaction and thermal exchange at actual endpoint populations.
//!
//! A bounded scalar conductance iteration surrounds the original constrained
//! reaction solve. For each positive trial G, Ce=3R/2*(ne0+xi), alpha=h*G, and
//! A=Ue0+Pe*t-chi*xi-Qstart give Te=(A+alpha*Th)/(Ce+alpha).
//! The inner thermal residual has derivative Ch+alpha*Ce/(Ce+alpha)>0.
//! These fixed-G derivatives are used only inside that constrained subsolve;
//! no monotonicity of the population-dependent source conductance is assumed.
//! Publication requires Q-Qstart=h*G(endpoint)*(Te-Th) within the declared
//! energy/temperature controls. No nonlinear trial allocates a gas state.

use super::*;
use fs_thermochem::equilibrium::composition::MAX_EQUILIBRIUM_SPECIES_V1;

enum ThermalTrial {
    Admitted(Point),
    Hot,
    Cold,
}

impl Problem<'_> {
    pub(super) fn state_at(
        &self,
        point: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, AtomicIonizationError> {
        let len = self.initial.amounts().len();
        let mut amounts = [Amount::new(0.0); MAX_EQUILIBRIUM_SPECIES_V1];
        amounts[..len].copy_from_slice(self.initial.amounts());
        for (i, n) in self.indices.into_iter().zip(self.amounts_at(point.extent)) {
            amounts[i] = Amount::new(n);
        }
        FrozenIonizedGasStateV1::new(
            self.initial.thermodynamics(),
            &amounts[..len],
            self.initial.volume(),
            Temperature::new(point.heavy_temperature),
            Temperature::new(point.electron_temperature),
            self.controls.relative_charge_tolerance,
            &mut *cancelled,
        )
        .map_err(|error| match error {
            FrozenIonizedGasStateErrorV1::Cancelled => AtomicIonizationError::Cancelled,
            error => AtomicIonizationError::State(error),
        })
    }

    pub(in crate::chemistry::ionization) fn coupling(
        &self,
        point: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExchangeConductance, AtomicIonizationError> {
        let Some(exchange) = &self.exchange else {
            return Ok(ExchangeConductance::default());
        };
        let len = self.initial.amounts().len();
        let mut amounts = [Amount::new(0.0); MAX_EQUILIBRIUM_SPECIES_V1];
        amounts[..len].copy_from_slice(self.initial.amounts());
        for (index, amount) in self.indices.into_iter().zip(self.amounts_at(point.extent)) {
            amounts[index] = Amount::new(amount);
        }
        // Retain the former frozen-state carrier-scaled charge admission,
        // using the same normalization and compensated summation. Neutral
        // background cannot dilute an imbalance as recombination removes
        // carriers. Only the final accepted state needs allocated ownership.
        let charges = self.initial.thermodynamics().charge_vector();
        let mut scale = 0.0_f64;
        for (index, amount) in amounts[..len].iter().enumerate() {
            if charges.get(index) != Some(0) {
                scale = scale.max(amount.value());
            }
        }
        if scale > 0.0 {
            let mut sums = [0.0_f64; 2];
            let mut corrections = [0.0_f64; 2];
            for (index, amount) in amounts[..len].iter().enumerate() {
                poll(cancelled)?;
                let charge = charges.get(index).expect("canonical charge axis") as f64;
                if charge == 0.0 {
                    continue;
                }
                let q = (amount.value() / scale) * charge;
                for (slot, value) in [q, q.abs()].into_iter().enumerate() {
                    let sum = sums[slot] + value;
                    corrections[slot] += if sums[slot].abs() >= value.abs() {
                        (sums[slot] - sum) + value
                    } else {
                        (value - sum) + sums[slot]
                    };
                    sums[slot] = sum;
                }
            }
            let relative = (sums[0] + corrections[0]).abs() / (sums[1] + corrections[1]);
            if !relative.is_finite() || relative > self.controls.relative_charge_tolerance {
                return Err(AtomicIonizationError::State(
                    FrozenIonizedGasStateErrorV1::ChargeImbalance {
                        relative_residual: relative,
                    },
                ));
            }
        }
        exchange
            .coupling_at_point(
                &amounts[..len],
                self.initial.volume(),
                point.heavy_temperature,
                point.electron_temperature,
                cancelled,
            )
            .map_err(collision_error)
    }

    fn thermal_trial(
        &self,
        start: &Point,
        end: f64,
        extent: f64,
        conductance: f64,
        minimum_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermalTrial, AtomicIonizationError> {
        let amounts = self.amounts_at(extent);
        if amounts.iter().any(|n| !n.is_finite() || *n < 0.0) || amounts[2] == 0.0 {
            return Err(AtomicIonizationError::OutsideDomain);
        }
        let ce = CV * amounts[2];
        let alpha = (end - start.time) * conductance;
        let denominator = ce + alpha;
        if !positive(ce) || !alpha.is_finite() || alpha < 0.0 || !positive(denominator) {
            return Err(AtomicIonizationError::Arithmetic(
                "coupled thermal coefficient",
            ));
        }
        let fraction = alpha / denominator;
        let electron_work = self.powers[1] * end;
        let heavy_target = self.powers[0] * end + start.exchange;
        let chi = self.channel.ionization_energy_j_mol();
        let initial_te = self.initial.electron_temperature().value();
        // An energy difference avoids subtracting absolute source formation
        // energies and preserves an unchanged source endpoint exactly.
        let exchange_at = |th: f64| {
            fraction
                * (CV * self.amounts[2] * (initial_te - th) + electron_work
                    - (chi + CV * th) * extent
                    - start.exchange)
        };
        if !electron_work.is_finite() || !heavy_target.is_finite() {
            return Err(AtomicIonizationError::Arithmetic("held coupled pool work"));
        }
        let mut lo = self.controls.minimum_heavy_temperature_k;
        let mut hi = self.controls.maximum_heavy_temperature_k;
        let rlo = self.heavy_limits[0] - heavy_target - exchange_at(lo);
        let rhi = self.heavy_limits[1] - heavy_target - exchange_at(hi);
        if !rlo.is_finite() || !rhi.is_finite() {
            return Err(AtomicIonizationError::Arithmetic("coupled caloric bracket"));
        }
        if rlo > 0.0 {
            return Ok(ThermalTrial::Cold);
        }
        if rhi < 0.0 {
            return Ok(ThermalTrial::Hot);
        }
        let mut th = if (lo..=hi).contains(&start.heavy_temperature) {
            start.heavy_temperature
        } else {
            lo + 0.5 * (hi - lo)
        };
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            let (energy, ch) = self.heavy_calorics(th, cancelled)?;
            let transfer = exchange_at(th);
            let residual = energy - heavy_target - transfer;
            let derivative = ch + fraction * ce;
            let temperature_scale =
                self.controls.temperature_tolerance_k + self.controls.relative_tolerance * th.abs();
            let tolerance = (self.controls.energy_tolerance_j / 128.0)
                .min(derivative * temperature_scale / 256.0);
            if !residual.is_finite() || !positive(derivative) {
                return Err(AtomicIonizationError::Arithmetic(
                    "coupled caloric residual",
                ));
            }
            if residual.abs() <= tolerance {
                let change = (electron_work - start.exchange - (chi + CV * initial_te) * extent
                    + alpha * (th - initial_te))
                    / denominator;
                let te = if change.is_finite() && change.abs() <= 0.5 * initial_te {
                    initial_te + change
                } else {
                    (self.electron_sensible + electron_work - chi * extent - start.exchange
                        + alpha * th)
                        / denominator
                };
                if te > self.controls.maximum_electron_temperature_k {
                    return Ok(ThermalTrial::Hot);
                }
                if !te.is_finite()
                    || te < self.controls.minimum_electron_temperature_k
                    || te < minimum_ratio * th
                {
                    return Ok(ThermalTrial::Cold);
                }
                // Implicit differentiation includes both changing electron
                // heat capacity and the thermal response of the heavy pool.
                let effective_capacity = ce + (alpha / (ch + alpha)) * ch;
                let dte = -(chi + CV * te) / effective_capacity;
                let point = Point {
                    time: end,
                    extent,
                    heavy_temperature: th,
                    electron_temperature: te,
                    exchange: start.exchange + transfer,
                    // The outer source-consistency solve assigns the actual
                    // neutral/Coulomb partition only after it converges.
                    coulomb_exchange: start.coulomb_exchange,
                    electron_temperature_extent_derivative: dte,
                };
                if !dte.is_finite()
                    || !point.exchange.is_finite()
                    || !point.coulomb_exchange.is_finite()
                {
                    return Err(AtomicIonizationError::Arithmetic(
                        "coupled thermal response",
                    ));
                }
                return Ok(ThermalTrial::Admitted(point));
            }
            if residual > 0.0 {
                hi = th;
            } else {
                lo = th;
            }
            let next = th - residual / derivative;
            let next = if next > lo && next < hi {
                next
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == th {
                return Err(AtomicIonizationError::Nonconvergence);
            }
            th = next;
        }
        Err(AtomicIonizationError::Nonconvergence)
    }

    fn reaction_converged(&self, point: &Point, residual: f64, derivative: f64) -> bool {
        let amounts = self.amounts_at(point.extent);
        let amount_scale = self.controls.absolute_amount_tolerance_mol
            + self.controls.relative_tolerance * amounts.into_iter().fold(f64::INFINITY, f64::min);
        let temperature_scale = self.controls.temperature_tolerance_k
            + self.controls.relative_tolerance * point.electron_temperature;
        let tolerance = (amount_scale / 128.0)
            .min(temperature_scale / (128.0 * point.electron_temperature_extent_derivative.abs()));
        residual.abs() <= tolerance * derivative.abs().max(1.0)
    }

    /// Close the original endpoint collision law around the constrained
    /// reaction/energy subsolve. The bracket, when available, consists only of
    /// observed opposite source residuals; its orientation is not assumed.
    /// Before finding one, a bounded positive secant or source update supplies
    /// the next conductance. Failure to converge requests adaptive refinement.
    pub(super) fn step_with_exchange(
        &self,
        start: &Point,
        end: f64,
        coupling: ExchangeConductance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Point, ExchangeConductance), AtomicIonizationError> {
        let dt = end - start.time;
        if !positive(dt) {
            return Err(AtomicIonizationError::Arithmetic(
                "positive finite coupled time step",
            ));
        }
        let mut conductance = coupling.total();
        let mut extent_guess = start.extent;
        let mut previous: Option<(f64, f64)> = None;
        let mut bracket: Option<[(f64, f64); 2]> = None;
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            let mut point =
                self.step_at_conductance(start, end, conductance, extent_guess, cancelled)?;
            let actual = self.coupling(&point, cancelled)?;
            let increment = point.exchange - start.exchange;
            let difference = point.electron_temperature - point.heavy_temperature;
            let expected = (dt * actual.total()) * difference;
            let residual = increment - expected;
            let ce = CV * self.amounts_at(point.extent)[2];
            let ch = self.heavy_calorics(point.heavy_temperature, cancelled)?.1;
            let temperature_scale = self.controls.temperature_tolerance_k
                + self.controls.relative_tolerance
                    * point.heavy_temperature.min(point.electron_temperature);
            let tolerance = (self.controls.energy_tolerance_j / 16.0)
                .min(ch.min(ce) * temperature_scale / 64.0);
            if !expected.is_finite() || !residual.is_finite() || !positive(tolerance) {
                return Err(AtomicIonizationError::Arithmetic(
                    "endpoint collision residual or tolerance",
                ));
            }
            if residual.abs() <= tolerance {
                if (increment > tolerance && difference < 0.0)
                    || (increment < -tolerance && difference > 0.0)
                {
                    return Err(AtomicIonizationError::Balance(
                        "endpoint heat transfer against temperature difference",
                    ));
                }
                point.coulomb_exchange = start.coulomb_exchange
                    + if actual.total() == 0.0 {
                        0.0
                    } else {
                        increment * (actual.coulomb / actual.total())
                    };
                if !point.coulomb_exchange.is_finite() {
                    return Err(AtomicIonizationError::Arithmetic(
                        "endpoint Coulomb transfer",
                    ));
                }
                return Ok((point, actual));
            }
            let source_residual = conductance - actual.total();
            if source_residual == 0.0 {
                return Err(AtomicIonizationError::Nonconvergence);
            }
            if let Some([lower, upper]) = bracket {
                bracket = if source_residual.is_sign_negative() == lower.1.is_sign_negative() {
                    Some([(conductance, source_residual), upper])
                } else {
                    Some([lower, (conductance, source_residual)])
                };
            } else if let Some(old) = previous
                && old.1.is_sign_negative() != source_residual.is_sign_negative()
            {
                let current = (conductance, source_residual);
                bracket = Some(if old.0 < conductance {
                    [old, current]
                } else {
                    [current, old]
                });
            }
            let secant = previous.map(|(old_g, old_residual)| {
                conductance
                    - source_residual * (conductance - old_g) / (source_residual - old_residual)
            });
            let reducing = previous
                .is_none_or(|(_, old_residual)| source_residual.abs() < 0.75 * old_residual.abs());
            let next = if let Some([lower, upper]) = bracket {
                match secant {
                    Some(next) if reducing && next > lower.0 && next < upper.0 => next,
                    _ => 0.5 * lower.0 + 0.5 * upper.0,
                }
            } else {
                match secant {
                    Some(next)
                        if next.is_finite()
                            && next >= 0.0
                            && next <= 4.0 * conductance.max(actual.total()) =>
                    {
                        next
                    }
                    _ => actual.total(),
                }
            };
            if !next.is_finite() || next < 0.0 || next == conductance {
                return Err(AtomicIonizationError::Nonconvergence);
            }
            previous = Some((conductance, source_residual));
            extent_guess = point.extent;
            conductance = next;
        }
        Err(AtomicIonizationError::Nonconvergence)
    }

    fn step_at_conductance(
        &self,
        start: &Point,
        end: f64,
        conductance: f64,
        initial_extent: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, AtomicIonizationError> {
        let mut minimum_ratio = self
            .exchange
            .as_ref()
            .expect("collision branch")
            .minimum_electron_heavy_ratio();
        for &i in &self.indices[..2] {
            let mass = self.initial.thermodynamics().models()[i]
                .molar_mass()
                .value()
                / fs_thermochem::electron_collision::AVOGADRO_PER_MOL;
            minimum_ratio = minimum_ratio.max(
                (fs_thermochem::electron_collision::ELECTRON_MASS_KG / mass)
                    / self
                        .channel
                        .forward_rate()
                        .maximum_target_thermal_speed_ratio()
                        .powi(2),
            );
        }
        // Both constrained temperatures decrease with extent. Find a supported
        // interior point without evaluating the impact source at unphysical
        // simplex endpoints. Source boundaries remain closed, without clipping.
        let physical_lo = (-self.amounts[1]).max((-self.amounts[2]).next_up());
        let physical_hi = self.amounts[0];
        let mut lo = physical_lo;
        let mut hi = physical_hi;
        let mut extent = initial_extent.max(lo).min(hi);
        let mut seed = None;
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            match self.thermal_trial(start, end, extent, conductance, minimum_ratio, cancelled)? {
                ThermalTrial::Admitted(point) => {
                    seed = Some(point);
                    break;
                }
                ThermalTrial::Hot => lo = extent,
                ThermalTrial::Cold => hi = extent,
            }
            let next = lo + 0.5 * (hi - lo);
            if next == extent || lo >= hi {
                break;
            }
            extent = next;
        }
        let mut point = seed.ok_or(AtomicIonizationError::OutsideDomain)?;
        let (mut residual, mut derivative) = self.residual(&point, start, cancelled)?;
        if self.reaction_converged(&point, residual, derivative) {
            return Ok(point);
        }

        // Search only the needed side until a sign-changing reaction bracket
        // is found. Most steps need one endpoint, avoiding a separate search
        // for both precise temperature-domain boundaries on every iteration.
        let upper = residual < 0.0;
        let mut admitted = point.extent;
        let mut boundary = if upper { physical_hi } else { physical_lo };
        let mut trial = boundary;
        let mut opposite = None;
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            match self.thermal_trial(start, end, trial, conductance, minimum_ratio, cancelled)? {
                ThermalTrial::Admitted(candidate) => {
                    // Detailed balance can certify the endpoint sign without
                    // an unnecessary Maxwellian impact quadrature.
                    if self.cheap_bracket_sign(&candidate, start, !upper)? {
                        opposite = Some(candidate.extent);
                        break;
                    }
                    let (value, slope) = self.residual(&candidate, start, cancelled)?;
                    if self.reaction_converged(&candidate, value, slope) {
                        return Ok(candidate);
                    }
                    if (upper && value >= 0.0) || (!upper && value <= 0.0) {
                        opposite = Some(candidate.extent);
                        break;
                    }
                    admitted = trial;
                }
                ThermalTrial::Hot if !upper => boundary = trial,
                ThermalTrial::Cold if upper => boundary = trial,
                _ => return Err(AtomicIonizationError::OutsideDomain),
            }
            let next = admitted + 0.5 * (boundary - admitted);
            if next == admitted || next == boundary {
                break;
            }
            trial = next;
        }
        let opposite = opposite.ok_or(AtomicIonizationError::OutsideDomain)?;
        lo = point.extent.min(opposite);
        hi = point.extent.max(opposite);
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            if self.reaction_converged(&point, residual, derivative) {
                return Ok(point);
            }
            if residual > 0.0 {
                hi = point.extent;
            } else {
                lo = point.extent;
            }
            let newton = point.extent - residual / derivative;
            let candidate = if derivative > 0.0 && newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if candidate == point.extent {
                return Err(AtomicIonizationError::Nonconvergence);
            }
            point = match self.thermal_trial(
                start,
                end,
                candidate,
                conductance,
                minimum_ratio,
                cancelled,
            )? {
                ThermalTrial::Admitted(point) => point,
                _ => return Err(AtomicIonizationError::OutsideDomain),
            };
            (residual, derivative) = self.residual(&point, start, cancelled)?;
        }
        Err(AtomicIonizationError::Nonconvergence)
    }
}

pub(super) fn collision_error(
    error: super::super::super::energy::ElectronHeavyEnergyError,
) -> AtomicIonizationError {
    use super::super::super::energy::ElectronHeavyEnergyError;
    match error {
        ElectronHeavyEnergyError::Cancelled => AtomicIonizationError::Cancelled,
        ElectronHeavyEnergyError::OutsideDomain => AtomicIonizationError::OutsideDomain,
        error => AtomicIonizationError::Caloric(error),
    }
}
