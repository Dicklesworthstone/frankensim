//! Caloric-only stations and bounded back-pressure regime selection.
use super::*;

#[derive(Clone, Copy, Debug)]
pub(super) struct Station {
    pub temperature: f64,
    pub pressure: f64,
    pub velocity: f64,
    pub sound_speed: f64,
}
impl Station {
    pub fn flux(self, r: f64) -> f64 {
        (self.pressure / (r * self.temperature)) * self.velocity
    }
    pub fn mach(self) -> f64 {
        self.velocity / self.sound_speed
    }
}

/// A single isentrope, with a separately retained total-enthalpy temperature.
struct Isentrope<'a> {
    parcel: &'a ClosedGasParcel,
    stagnation_temperature: f64,
    anchor_temperature: f64,
    anchor_pressure: f64,
}
impl Isentrope<'_> {
    fn sample(&self, t: f64) -> Result<Station, GasEnergyError> {
        let thermo = self.parcel.model().thermodynamics();
        let r = self.parcel.state().state().specific_gas_constant;
        let properties = thermo
            .evaluate_properties(Temperature::new(t), Pressure::new(self.anchor_pressure))
            .map_err(source_error)?;
        let (cp, cv) = (
            properties.mass_specific_cp().value(),
            properties.mass_specific_cv().value(),
        );
        let gamma = cp / cv;
        if !positive(cp)
            || !positive(cv)
            || !gamma.is_finite()
            || gamma <= 1.0
            || gamma > 5.0 / 3.0 + 1e-9
        {
            return Err(invalid(
                "rocket caloric sample has inadmissible heat capacities",
            ));
        }
        let dh = thermo
            .specific_enthalpy_change(
                Temperature::new(self.stagnation_temperature),
                Temperature::new(t),
            )
            .map_err(source_error)?
            .value();
        let ds = thermo
            .specific_entropy_change(
                Temperature::new(self.anchor_temperature),
                Pressure::new(self.anchor_pressure),
                Temperature::new(t),
                Pressure::new(self.anchor_pressure),
            )
            .map_err(source_error)?
            .value();
        let pressure = self.anchor_pressure * det::exp(ds / r);
        let speed_squared = -2.0 * dh;
        if !positive(pressure) || !speed_squared.is_finite() || speed_squared < 0.0 {
            return Err(invalid(
                "rocket expansion pressure or kinetic energy is unrepresentable",
            ));
        }
        Ok(Station {
            temperature: t,
            pressure,
            velocity: det::sqrt(speed_squared),
            sound_speed: det::sqrt(gamma * r * t),
        })
    }

    fn mass_flux(
        &self,
        target: f64,
        lo: f64,
        hi: f64,
        controls: RocketNozzleSolve,
        iterations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Station, GasEnergyError> {
        if !positive(target) {
            return Err(invalid("rocket mass flux is unrepresentable"));
        }
        let r = self.parcel.state().state().specific_gas_constant;
        let t = root(
            lo,
            hi,
            0.25 * controls.relative_mass_flow_tolerance,
            controls.nozzle.maximum_iterations,
            iterations,
            cancelled,
            |t| Ok(self.sample(t)?.flux(r) / target - 1.0),
        )?;
        self.sample(t)
    }
}

/// Safeguarded monotone bracket with a bounded, polled iteration budget.
/// Either sign orientation is allowed. The caller supplies a physical residual.
pub(super) fn root(
    mut lo: f64,
    mut hi: f64,
    tolerance: f64,
    maximum_iterations: usize,
    total_iterations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
    mut residual: impl FnMut(f64) -> Result<f64, GasEnergyError>,
) -> Result<f64, GasEnergyError> {
    poll(cancelled)?;
    let mut lower = residual(lo)?;
    poll(cancelled)?;
    let upper = residual(hi)?;
    if !lower.is_finite() || !upper.is_finite() {
        return Err(invalid("nozzle root has a nonfinite residual"));
    }
    if lower.abs() <= tolerance {
        return Ok(lo);
    }
    if upper.abs() <= tolerance {
        return Ok(hi);
    }
    if lower.is_sign_negative() == upper.is_sign_negative() {
        return Err(invalid(
            "nozzle root is outside its admitted physical bracket",
        ));
    }
    for _ in 0..maximum_iterations {
        let mid = lo + 0.5 * (hi - lo);
        if mid <= lo || mid >= hi {
            return Err(invalid("nozzle root exhausted floating-point progress"));
        }
        poll(cancelled)?;
        let value = residual(mid)?;
        *total_iterations += 1;
        if !value.is_finite() {
            return Err(invalid("nozzle root produced a nonfinite residual"));
        }
        if value.abs() <= tolerance {
            return Ok(mid);
        }
        if value.is_sign_negative() == lower.is_sign_negative() {
            lo = mid;
            lower = value;
        } else {
            hi = mid;
        }
    }
    Err(invalid("nozzle root exhausted its iteration budget"))
}

fn unchoked(
    initial: &Isentrope<'_>,
    minimum_temperatures: [f64; 2],
    ambient: f64,
    nozzle: RocketNozzle,
    controls: RocketNozzleSolve,
    mut iterations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RocketNozzleFlow, GasEnergyError> {
    if ambient <= 0.0 {
        return Err(invalid(
            "admitted temperature range does not reach a vacuum nozzle's sonic state",
        ));
    }
    let t0 = initial.stagnation_temperature;
    let r = initial.parcel.state().state().specific_gas_constant;
    let temperature = root(
        minimum_temperatures[1],
        t0,
        0.25 * controls.relative_mass_flow_tolerance,
        controls.nozzle.maximum_iterations,
        &mut iterations,
        cancelled,
        |t| Ok(initial.sample(t)?.pressure / ambient - 1.0),
    )?;
    let exit = initial.sample(temperature)?;
    let rate = exit.flux(r) * nozzle.exit_area_m2;
    let throat = initial.mass_flux(
        rate / nozzle.throat_area_m2,
        minimum_temperatures[0],
        t0,
        controls,
        &mut iterations,
        cancelled,
    )?;
    assemble(
        initial.parcel,
        ambient,
        nozzle,
        controls,
        throat,
        exit,
        rate,
        RocketNozzleRegime::UnchokedSubsonic,
        None,
        iterations,
        cancelled,
    )
}

pub(super) fn flow(
    parcel: &ClosedGasParcel,
    ambient: f64,
    nozzle: RocketNozzle,
    controls: RocketNozzleSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RocketNozzleFlow, GasEnergyError> {
    let state = parcel.state().state();
    let (t0, p0, r) = (
        state.temperature,
        state.pressure,
        state.specific_gas_constant,
    );
    if t0 < controls.nozzle.minimum_temperature_k || t0 > controls.nozzle.maximum_temperature_k {
        return Err(invalid(
            "rocket reservoir temperature is outside the admitted source bracket",
        ));
    }
    let initial = Isentrope {
        parcel,
        stagnation_temperature: t0,
        anchor_temperature: t0,
        anchor_pressure: p0,
    };
    let mut iterations = 0;
    let minimum_temperature = controls.nozzle.minimum_temperature_k;
    let lower = initial.sample(minimum_temperature)?;
    if 0.5 * (lower.sound_speed * lower.sound_speed - lower.velocity * lower.velocity)
        > controls.nozzle.energy_tolerance_j_kg
    {
        // The entire admitted branch is subsonic. A small pressure drop can
        // have source-valid throat/exit states even when a hypothetical sonic
        // temperature is below the source's gas-phase validity boundary.
        return unchoked(
            &initial,
            [minimum_temperature; 2],
            ambient,
            nozzle,
            controls,
            iterations,
            cancelled,
        );
    }
    let sonic_temperature = root(
        controls.nozzle.minimum_temperature_k,
        t0,
        controls.nozzle.energy_tolerance_j_kg,
        controls.nozzle.maximum_iterations,
        &mut iterations,
        cancelled,
        |t| {
            let at = initial.sample(t)?;
            Ok(0.5 * (at.sound_speed * at.sound_speed - at.velocity * at.velocity))
        },
    )?;
    let sonic = initial.sample(sonic_temperature)?;
    let choked_rate = sonic.flux(r) * nozzle.throat_area_m2;
    let target_exit_flux = choked_rate / nozzle.exit_area_m2;
    // The limiting downstream subsonic branch is the zero-strength shock at
    // the throat. Its pressure cannot be lower than the sonic pressure.
    // Vacuum/low-back-pressure calls therefore need no extra subsonic root.
    // No hypothetical transport state is resolved here.
    if ambient >= sonic.pressure {
        let critical_exit = initial.mass_flux(
            target_exit_flux,
            sonic_temperature,
            t0,
            controls,
            &mut iterations,
            cancelled,
        )?;
        if ambient >= critical_exit.pressure {
            return unchoked(
                &initial,
                [sonic_temperature, critical_exit.temperature],
                ambient,
                nozzle,
                controls,
                iterations,
                cancelled,
            );
        }
    }

    let supersonic_exit = if nozzle.exit_area_m2 == nozzle.throat_area_m2 {
        sonic
    } else {
        initial.mass_flux(
            target_exit_flux,
            controls.nozzle.minimum_temperature_k,
            sonic_temperature,
            controls,
            &mut iterations,
            cancelled,
        )?
    };
    if ambient <= supersonic_exit.pressure {
        return assemble(
            parcel,
            ambient,
            nozzle,
            controls,
            sonic,
            supersonic_exit,
            choked_rate,
            RocketNozzleRegime::ChokedSupersonic,
            None,
            iterations,
            cancelled,
        );
    }

    let exit_shock = shock::calculate(
        parcel,
        supersonic_exit,
        t0,
        controls.nozzle,
        &mut iterations,
        cancelled,
    )?;
    if ambient < exit_shock.downstream.pressure {
        return assemble(
            parcel,
            ambient,
            nozzle,
            controls,
            sonic,
            supersonic_exit,
            choked_rate,
            RocketNozzleRegime::ChokedExternalCompression,
            None,
            iterations,
            cancelled,
        );
    }

    // Parametrize shock position by its upstream static temperature. Temperature
    // increases monotonically from the supersonic exit back to the sonic throat.
    // The actual area follows from the same conserved choked mass flow.
    let mut lo = supersonic_exit.temperature;
    let mut hi = sonic.temperature;
    let pressure_tolerance = controls.relative_mass_flow_tolerance * p0;
    let mut shock = exit_shock;
    let mut exit = shock.downstream;
    let mut accepted = (exit.pressure - ambient).abs() <= 0.25 * pressure_tolerance;
    for _ in 0..controls.nozzle.maximum_iterations {
        if accepted {
            break;
        }
        let temperature = lo + 0.5 * (hi - lo);
        if temperature <= lo || temperature >= hi {
            return Err(invalid(
                "internal nozzle shock exhausted floating-point progress",
            ));
        }
        poll(cancelled)?;
        let upstream = initial.sample(temperature)?;
        shock = shock::calculate(
            parcel,
            upstream,
            t0,
            controls.nozzle,
            &mut iterations,
            cancelled,
        )?;
        let downstream = Isentrope {
            parcel,
            stagnation_temperature: t0,
            anchor_temperature: shock.downstream.temperature,
            anchor_pressure: shock.downstream.pressure,
        };
        exit = downstream.mass_flux(
            target_exit_flux,
            shock.downstream.temperature,
            t0,
            controls,
            &mut iterations,
            cancelled,
        )?;
        iterations += 1;
        let residual = exit.pressure - ambient;
        accepted = residual.abs() <= 0.25 * pressure_tolerance;
        if residual < 0.0 {
            lo = temperature;
        } else {
            hi = temperature;
        }
    }
    if !accepted {
        return Err(invalid(
            "internal nozzle shock exhausted its back-pressure iteration budget",
        ));
    }
    let area = choked_rate / shock.upstream.flux(r);
    if !area.is_finite()
        || area < nozzle.throat_area_m2 * (1.0 - controls.relative_mass_flow_tolerance)
        || area > nozzle.exit_area_m2 * (1.0 + controls.relative_mass_flow_tolerance)
    {
        return Err(invalid(
            "normal shock area lies outside the physical diverging nozzle",
        ));
    }
    assemble(
        parcel,
        ambient,
        nozzle,
        controls,
        sonic,
        exit,
        choked_rate,
        RocketNozzleRegime::ChokedNormalShock,
        Some((area, shock)),
        iterations,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble(
    parcel: &ClosedGasParcel,
    ambient: f64,
    nozzle: RocketNozzle,
    controls: RocketNozzleSolve,
    throat_at: Station,
    exit_at: Station,
    rate: f64,
    regime: RocketNozzleRegime,
    shock: Option<(f64, shock::RawShock)>,
    iterations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RocketNozzleFlow, GasEnergyError> {
    poll(cancelled)?;
    let throat = parcel.model().evaluate(&updated_point(
        parcel.point(),
        throat_at.temperature,
        throat_at.pressure,
    )?)?;
    poll(cancelled)?;
    let exit = parcel.model().evaluate(&updated_point(
        parcel.point(),
        exit_at.temperature,
        exit_at.pressure,
    )?)?;
    let normal_shock = if let Some((area, jump)) = shock {
        Some(RocketNormalShock {
            area_m2: area,
            flow: shock::resolve(parcel, jump, controls.nozzle, cancelled)?,
        })
    } else {
        None
    };
    let entropy_increase = normal_shock
        .as_ref()
        .map_or(0.0, |s| s.flow.entropy_increase_j_kg_k);
    let state = parcel.state().state();
    let thermo = parcel.model().thermodynamics();
    let energy_residual = thermo
        .specific_enthalpy_change(
            Temperature::new(state.temperature),
            Temperature::new(exit_at.temperature),
        )
        .map_err(source_error)?
        .value()
        + 0.5 * exit_at.velocity * exit_at.velocity;
    let entropy_residual = state.temperature
        * (thermo
            .specific_entropy_change(
                Temperature::new(state.temperature),
                Pressure::new(state.pressure),
                Temperature::new(exit_at.temperature),
                Pressure::new(exit_at.pressure),
            )
            .map_err(source_error)?
            .value()
            - entropy_increase);
    let throat_rate = throat.state().density * throat_at.velocity * nozzle.throat_area_m2;
    let exit_rate = exit.state().density * exit_at.velocity * nozzle.exit_area_m2;
    let mass_residual = exit_rate / throat_rate - 1.0;
    let (throat_mach, exit_mach) = (
        throat_at.velocity / throat.state().sound_speed,
        exit_at.velocity / exit.state().sound_speed,
    );
    let matches_pressure = matches!(
        regime,
        RocketNozzleRegime::UnchokedSubsonic | RocketNozzleRegime::ChokedNormalShock
    );
    let pressure_residual = matches_pressure.then_some(exit_at.pressure - ambient);
    let mach_tolerance = 4.0 * controls.nozzle.energy_tolerance_j_kg
        / (throat.state().sound_speed * throat.state().sound_speed);
    if !positive(rate)
        || !positive(throat_rate)
        || !positive(exit_rate)
        || !mass_residual.is_finite()
        || mass_residual.abs() > controls.relative_mass_flow_tolerance
        || (throat_rate / rate - 1.0).abs() > controls.relative_mass_flow_tolerance
        || (exit_rate / rate - 1.0).abs() > controls.relative_mass_flow_tolerance
        || !positive(throat_mach)
        || !positive(exit_mach)
        || (matches_pressure && exit_mach > 1.0 + mach_tolerance)
        || (!matches_pressure && exit_mach < 1.0 - mach_tolerance)
        || (regime == RocketNozzleRegime::UnchokedSubsonic && throat_mach > 1.0 + mach_tolerance)
        || (regime != RocketNozzleRegime::UnchokedSubsonic
            && (throat_mach - 1.0).abs() > mach_tolerance)
        || !energy_residual.is_finite()
        || energy_residual.abs() > 2.0 * controls.nozzle.energy_tolerance_j_kg
        || !entropy_residual.is_finite()
        || entropy_residual.abs() > 4.0 * controls.nozzle.energy_tolerance_j_kg
        || pressure_residual
            .is_some_and(|p| p.abs() > controls.relative_mass_flow_tolerance * state.pressure)
    {
        return Err(invalid(
            "rocket exit failed mass, energy, entropy, pressure or Mach admission",
        ));
    }
    let momentum = rate * exit_at.velocity;
    let pressure_thrust = (exit_at.pressure - ambient) * nozzle.exit_area_m2;
    let thrust = momentum + pressure_thrust;
    let effective = thrust / rate;
    let characteristic = state.pressure * nozzle.throat_area_m2 / rate;
    let specific_impulse = effective / STANDARD_GRAVITY_M_S2;
    if !positive(momentum)
        || !pressure_thrust.is_finite()
        || !positive(thrust)
        || !positive(effective)
        || !positive(characteristic)
        || !positive(specific_impulse)
    {
        return Err(invalid("rocket thrust or performance is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(RocketNozzleFlow {
        throat,
        exit,
        regime,
        throat_mach_number: throat_mach,
        normal_shock,
        mass_flow_kg_s: rate,
        exit_velocity_m_s: exit_at.velocity,
        exit_mach_number: exit_mach,
        momentum_thrust_n: momentum,
        pressure_thrust_n: pressure_thrust,
        thrust_n: thrust,
        effective_exhaust_velocity_m_s: effective,
        characteristic_velocity_m_s: characteristic,
        specific_impulse_s: specific_impulse,
        relative_mass_flow_residual: mass_residual,
        energy_residual_j_kg: energy_residual,
        entropy_energy_residual_j_kg: entropy_residual,
        entropy_increase_j_kg_k: entropy_increase,
        back_pressure_residual_pa: pressure_residual,
        iterations,
    })
}
