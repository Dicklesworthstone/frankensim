//! Coupled impulse and work, using the existing exact slip-decay map once.

use super::{RecoilLaw, SlidingSupport, WallExchangeError, nonnegative, positive};

pub(super) struct Motion {
    pub momentum: f64,
    pub support_velocity_m_s: f64,
    pub gas_impulse_ns: f64,
    pub work_to_gas_j: f64,
    pub dissipation_j: f64,
    pub support_kinetic_energy_change_j: f64,
    pub momentum_residual_ns: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    momentum: f64,
    density: f64,
    velocity: f64,
    volume: f64,
    dt: f64,
    support: SlidingSupport,
    law: RecoilLaw,
) -> Result<Motion, WallExchangeError> {
    let mass = density * volume;
    let slip = velocity - support.velocity_m_s;
    if !positive(density)
        || !positive(volume)
        || !positive(mass)
        || !positive(support.mass_kg)
        || !support.velocity_m_s.is_finite()
        || !slip.is_finite()
        || !nonnegative(dt)
        || !nonnegative(law.linear_rate_per_s)
        || !nonnegative(law.quadratic_per_m)
    {
        return Err(WallExchangeError::Invalid(
            "finite-support drag mass, rate or velocity",
        ));
    }
    let reduced_mass = if mass <= support.mass_kg {
        mass / (1.0 + mass / support.mass_kg)
    } else {
        support.mass_kg / (1.0 + support.mass_kg / mass)
    };
    let relative_duration = dt * (mass / reduced_mass);
    if !positive(reduced_mass) || !nonnegative(relative_duration) {
        return Err(WallExchangeError::Invalid(
            "finite-support relative relaxation interval",
        ));
    }
    let (factor, removed) = super::super::drag::combined_decay(
        law.linear_rate_per_s,
        law.quadratic_per_m,
        slip.abs(),
        relative_duration,
    );
    let gas_weight = reduced_mass / mass;
    let support_weight = reduced_mass / support.mass_kg;
    let next_momentum = if slip == 0.0 || removed == 0.0 {
        momentum
    } else if removed <= 0.5 {
        momentum - (density * gas_weight) * slip * removed
    } else {
        let center = if gas_weight <= 0.5 {
            velocity - gas_weight * slip
        } else {
            support.velocity_m_s + support_weight * slip
        };
        density * (center + gas_weight * slip * factor)
    };
    // Retain actual accepted gas impulse, then give the support its opposite.
    // Endpoint rounding is audited rather than hidden in an energy correction.
    let impulse = (next_momentum - momentum) * volume;
    let support_delta = -impulse / support.mass_kg;
    let next_support = support.velocity_m_s + support_delta;
    let actual_support_impulse = support.mass_kg * (next_support - support.velocity_m_s);
    let out = Motion {
        momentum: next_momentum,
        support_velocity_m_s: next_support,
        gas_impulse_ns: impulse,
        work_to_gas_j: impulse * (support.velocity_m_s + 0.5 * support_delta),
        dissipation_j: -impulse * (slip + 0.5 * (impulse / reduced_mass)),
        support_kinetic_energy_change_j: actual_support_impulse
            * (0.5 * next_support + 0.5 * support.velocity_m_s),
        momentum_residual_ns: impulse + actual_support_impulse,
    };
    if ![
        out.momentum,
        out.support_velocity_m_s,
        out.gas_impulse_ns,
        out.work_to_gas_j,
        out.dissipation_j,
        out.support_kinetic_energy_change_j,
        out.momentum_residual_ns,
    ]
    .into_iter()
    .all(f64::is_finite)
        || out.dissipation_j < 0.0
    {
        return Err(WallExchangeError::Invalid(
            "finite-support mechanical exchange",
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
