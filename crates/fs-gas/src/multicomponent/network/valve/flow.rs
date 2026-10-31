//! One shared source-owned nozzle solve and donor composition per valve.
use super::*;
use crate::dynamics::characteristic::shared::nozzle;
use crate::dynamics::wall::shared::WallError;
use crate::multicomponent::wall::LocalCaloric;

pub(super) fn evaluate(
    gas: &MixtureGas,
    states: [Cell<'_>; 2],
    orifice: Orifice,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ValveFlow, ValveError> {
    poll(&mut cancelled)?;
    if states.iter().any(|q| q.momentum != 0.0) {
        return Err(ValveError::Invalid("valve endpoints must be stationary compact volumes"));
    }
    let p = [gas.primitive(states[0], &mut cancelled)?, gas.primitive(states[1], &mut cancelled)?];
    admitted(gas, states, p, orifice, &mut cancelled)
}

pub(super) fn admitted(
    gas: &MixtureGas,
    states: [Cell<'_>; 2],
    p: [Primitive; 2],
    orifice: Orifice,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ValveFlow, ValveError> {
    let n = gas.species_count();
    let mut area = orifice.effective_area().value();
    if area == 0.0 || p[0].pressure == p[1].pressure {
        poll(cancelled)?;
        return Ok(ValveFlow::zero(n));
    }
    let donor = usize::from(p[1].pressure > p[0].pressure);
    if (donor == 1 && orifice.direction == ValveDirection::FirstToSecond)
        || (donor == 0 && orifice.direction == ValveDirection::SecondToFirst)
    {
        poll(cancelled)?;
        return Ok(ValveFlow::zero(n));
    }
    if let Some(opening) = orifice.pressure_opening {
        let fraction = opening.fraction(Pressure::new(p[donor].pressure - p[1 - donor].pressure))?;
        if fraction == 0.0 {
            poll(cancelled)?;
            return Ok(ValveFlow::zero(n));
        }
        area *= fraction;
        if area == 0.0 {
            return Err(ValveError::Invalid("nonzero pressure-actuated orifice area underflowed"));
        }
    }
    let law = LocalCaloric::new(gas, states[donor].partial_densities)?;
    let (exit, choked) = nozzle::state(&law, p[donor], p[1 - donor].pressure, cancelled)
        .map_err(|error| match error {
            WallError::Source(e) => ValveError::from(e),
            WallError::Cancelled => ValveError::Cancelled,
            WallError::Invalid(reason) => ValveError::Invalid(reason),
            WallError::CaloricPath(reason) => ValveError::Flow(FlowError::CharacteristicPath(reason)),
            WallError::Iterations => ValveError::Flow(FlowError::TemperatureIterations),
        })?;
    let magnitude = (area * exit.density) * exit.velocity;
    if !positive(magnitude) { return Err(ValveError::Invalid("unrepresentable nonzero valve flow")); }
    let sign = if donor == 0 { 1.0 } else { -1.0 };
    let mut out = ValveFlow::zero(n);
    for k in 0..n {
        out.species_kg_s[k] = sign * magnitude * (states[donor].partial_densities[k] / p[donor].density);
        if states[donor].partial_densities[k] > 0.0 && out.species_kg_s[k] == 0.0 {
            return Err(ValveError::Invalid("nonzero valve species rate underflowed"));
        }
    }
    out.mass_kg_s = sum(out.species_kg_s[..n].iter().copied());
    // Stagnation enthalpy already includes the eventual dissipated jet energy.
    // Adding a second kinetic-energy/heat-release source would double count it.
    out.energy_w = out.mass_kg_s
        * (states[donor].energy / p[donor].density + p[donor].pressure / p[donor].density);
    out.donor = Some(donor);
    out.exit = Some(exit);
    out.choked = choked;
    if !out.energy_w.is_finite() || out.mass_kg_s == 0.0 || !out.mass_kg_s.is_finite()
        || !out.species_kg_s[..n].iter().all(|x| x.is_finite())
    { return Err(ValveError::Invalid("valve species or enthalpy rate overflow")); }
    poll(cancelled)?;
    Ok(out)
}
