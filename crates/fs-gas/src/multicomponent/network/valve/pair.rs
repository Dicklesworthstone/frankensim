//! Inventory-conserving finite transfer with an implicit endpoint mass rate.
use super::*;
use crate::dynamics::wall::shared::CaloricLaw;
use crate::multicomponent::wall::LocalCaloric;

pub(super) struct Update {
    pub(super) states: [[f64; MAX_COMPONENTS]; 2],
    pub(super) exchange: ValveExchange,
    primitives: [Primitive; 2],
}

pub(super) fn solve(
    gas: &MixtureGas,
    plenums: [&Plenum; 2],
    orifice: Orifice,
    dt: f64,
    controls: ValveControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Update, ValveError> {
    poll(cancelled)?;
    let n = gas.species_count();
    let initial = [plenums[0].state(), plenums[1].state()];
    let volumes = [plenums[0].volume_m3(), plenums[1].volume_m3()];
    let mut unchanged = [[0.0; MAX_COMPONENTS]; 2];
    for side in 0..2 {
        unchanged[side][..n].copy_from_slice(initial[side].partial_densities);
        unchanged[side][n + 1] = initial[side].energy;
    }
    let primitives = [gas.primitive(initial[0], &mut *cancelled)?, gas.primitive(initial[1], &mut *cancelled)?];
    let rate = flow::admitted(gas, initial, primitives, orifice, cancelled)?;
    let Some(donor) = rate.donor else {
        return Ok(Update { states: unchanged, exchange: ValveExchange::zero(n), primitives });
    };
    let sign = if donor == 0 { 1.0 } else { -1.0 };
    let p = primitives[donor];
    let mass = p.density * volumes[donor];
    let requested = dt * rate.mass_kg_s.abs();
    if !positive(dt) || !positive(mass) || !positive(requested) {
        return Err(ValveError::Invalid("unrepresentable valve interval or requested transfer"));
    }
    let (mut lo, mut hi) = (0.0, mass);
    let mut x = requested.min(0.5 * mass);
    let mut hit_domain = false;
    for iteration in 1..=controls.maximum_iterations {
        poll(cancelled)?;
        if !positive(x) || x <= lo || x >= hi {
            return Err(if hit_domain { ValveError::CaloricDomain } else { ValveError::IterationLimit });
        }
        let candidate = candidate(gas, initial, volumes, donor, p, x, cancelled);
        let mut update = match candidate {
            Ok(value) => value,
            Err(ValveError::CaloricDomain) | Err(ValveError::Flow(FlowError::EnergyDomain)) => {
                hit_domain = true;
                hi = x;
                x = lo.midpoint(hi);
                continue;
            }
            Err(error) => return Err(error),
        };
        // An over-equilibrated trial is an upper mass bound, not a request to
        // solve an unused reverse nozzle (which may leave the other source path).
        if update.primitives[donor].pressure <= update.primitives[1 - donor].pressure {
            hi = x;
            x = lo.midpoint(hi);
            continue;
        }
        let endpoint = flow::admitted(gas,
            [cell(&update.states[0][..n + 2], n), cell(&update.states[1][..n + 2], n)],
            update.primitives, orifice, cancelled)?;
        // A pressure-actuated aperture can close before the tanks equilibrate.
        // Such a trial is an upper transfer bound, never an accepted tiny leak
        // justified only by the absolute mass tolerance.
        if endpoint.donor.is_none() {
            hi = x;
            x = lo.midpoint(hi);
            continue;
        }
        let actual_mass = sign * update.exchange.mass_to_second_kg();
        let expected_mass = dt * sign * endpoint.mass_kg_s;
        let residual = actual_mass - expected_mass;
        if !residual.is_finite() { return Err(ValveError::Invalid("implicit valve mass equation overflow")); }
        let tolerance = controls.absolute_mass_tolerance_kg
            + controls.relative_mass_tolerance * actual_mass.abs().max(expected_mass.abs());
        if expected_mass >= 0.0 && residual.abs() <= tolerance {
            update.exchange.iterations = iteration;
            update.exchange.choked_endpoints = usize::from(endpoint.choked);
            update.exchange.maximum_mass_residual_kg = residual.abs();
            poll(cancelled)?;
            return Ok(update);
        }
        if residual < 0.0 { lo = x; } else { hi = x; }
        x = lo.midpoint(hi);
    }
    Err(if hit_domain { ValveError::CaloricDomain } else { ValveError::IterationLimit })
}

#[allow(clippy::too_many_arguments)]
fn candidate(
    gas: &MixtureGas,
    initial: [Cell<'_>; 2],
    volumes: [f64; 2],
    donor: usize,
    p: Primitive,
    withdrawn: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Update, ValveError> {
    let n = gas.species_count();
    let receiver = 1 - donor;
    let volume = volumes[donor];
    let mass = p.density * volume;
    let fraction = withdrawn / mass;
    let remaining = 1.0 - fraction;
    if !positive(fraction) || !positive(remaining) { return Err(ValveError::CaloricDomain); }
    let law = LocalCaloric::new(gas, initial[donor].partial_densities)?;
    let temperature = blowdown_temperature(&law, p, remaining, cancelled)?;
    // Integrate d(m*u)=h*dm along the original donor isentrope. This is the
    // actual finite enthalpy transfer, not withdrawn mass times a stale h.
    // Evaluate differences before formation-energy rounding can obscure them.
    let cooling = law.internal_energy_change(temperature, p.temperature, cancelled)?;
    let energy = withdrawn * (initial[donor].energy / p.density)
        + (mass - withdrawn) * cooling;
    if !energy.is_finite() || energy < 0.0 { return Err(ValveError::Invalid("valve blowdown energy")); }
    let sign = if donor == 0 { 1.0 } else { -1.0 };
    let mut result = Update { states: [[0.0; MAX_COMPONENTS]; 2], exchange: ValveExchange::zero(n), primitives: [Primitive::default(); 2] };
    for k in 0..n {
        let transfer = (initial[donor].partial_densities[k] * fraction) * volume;
        if initial[donor].partial_densities[k] > 0.0 && transfer == 0.0 {
            return Err(ValveError::Invalid("nonzero valve species transfer underflowed"));
        }
        result.exchange.species_to_second_kg[k] = sign * transfer;
        result.states[donor][k] = initial[donor].partial_densities[k] - transfer / volume;
        result.states[receiver][k] = initial[receiver].partial_densities[k] + transfer / volumes[receiver];
    }
    result.exchange.energy_to_second_j = sign * energy;
    result.states[donor][n + 1] = initial[donor].energy - energy / volume;
    result.states[receiver][n + 1] = initial[receiver].energy + energy / volumes[receiver];
    for (side, state) in result.states.iter().enumerate() {
        if !gas.admissible(&state[..n + 2]) { return Err(ValveError::CaloricDomain); }
        result.primitives[side] = gas.decode(&state[..n + 2], None, cancelled)?;
    }
    Ok(result)
}

fn blowdown_temperature(
    law: &LocalCaloric<'_>,
    initial: Primitive,
    remaining: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ValveError> {
    let domain = law.domain();
    let (mut lo, mut hi) = (domain.minimum_temperature, initial.temperature);
    let entropy = |t: f64, cancelled: &mut _| -> Result<f64, ValveError> {
        let pressure = (initial.pressure * remaining) * (t / initial.temperature);
        if !positive(pressure) { return Err(ValveError::CaloricDomain); }
        Ok(law.entropy_change(initial.temperature, initial.pressure, t, pressure, cancelled)?)
    };
    // Along a fixed-volume blowdown p=rho*R*T and dS/dT=Cv/T, not Cp/T.
    let cv0 = law.heat_capacities(initial.temperature, initial.pressure, cancelled)?.cv;
    let tolerance = domain.absolute_energy_tolerance / initial.temperature
        + domain.relative_energy_tolerance * cv0;
    if entropy(lo, cancelled)? > tolerance { return Err(ValveError::CaloricDomain); }
    let mut t = hi;
    for _ in 0..domain.maximum_iterations {
        poll(cancelled)?;
        let residual = entropy(t, cancelled)?;
        if residual.abs() <= tolerance { return Ok(t); }
        if residual < 0.0 { lo = t; } else { hi = t; }
        let cv = law.heat_capacities(t, initial.pressure, cancelled)?.cv;
        if !positive(cv) { return Err(ValveError::Invalid("valve blowdown heat capacity")); }
        let trial = t - residual * (t / cv);
        t = if trial > lo && trial < hi { trial } else { lo.midpoint(hi) };
    }
    Err(ValveError::Flow(FlowError::TemperatureIterations))
}
