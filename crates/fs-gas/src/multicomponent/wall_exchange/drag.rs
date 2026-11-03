//! Exact frozen relative drag, with the same stationary arithmetic as before.

use super::{WallDrag, WallExchangeError, det, nonnegative};

pub(super) struct Update {
    pub momentum: f64,
    pub energy: f64,
    pub work_j: f64,
    pub heating_j: f64,
    pub energy_residual_j: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    law: WallDrag,
    initial_momentum: f64,
    initial_energy: f64,
    density: f64,
    velocity: f64,
    volume: f64,
    dt: f64,
) -> Result<Update, WallExchangeError> {
    let wall = law.moving_velocity().unwrap_or(0.0);
    let slip = velocity - wall;
    if !slip.is_finite() {
        return Err(WallExchangeError::Invalid("wall relative velocity"));
    }
    // Keep both the surviving and removed fractions: expm1 and d/(1+d)
    // retain weak resolved impulses that subtracting a near-unity factor loses.
    let (factor, removed) = match law {
        WallDrag::None => (1.0, 0.0),
        WallDrag::Linear { rate_per_s } | WallDrag::MovingLinear { rate_per_s, .. } => {
            let exponent = rate_per_s * dt;
            if !nonnegative(exponent) {
                return Err(WallExchangeError::Invalid("wall drag exponent"));
            }
            (
                det::exp(-exponent),
                if wall == 0.0 {
                    0.0
                } else {
                    -det::expm1(-exponent)
                },
            )
        }
        WallDrag::Darcy {
            friction_factor,
            hydraulic_diameter_m,
        }
        | WallDrag::MovingDarcy {
            friction_factor,
            hydraulic_diameter_m,
            ..
        } => {
            let damping = 0.5 * friction_factor * dt * (slip.abs() / hydraulic_diameter_m);
            if !nonnegative(damping) {
                return Err(WallExchangeError::Invalid("wall Darcy increment"));
            }
            (1.0 / (1.0 + damping), damping / (1.0 + damping))
        }
    };
    if wall == 0.0 {
        let momentum = initial_momentum * factor;
        let old_kinetic = 0.5 * initial_momentum * (initial_momentum / density);
        let kinetic = 0.5 * momentum * (momentum / density);
        return Ok(Update {
            momentum,
            energy: initial_energy,
            work_j: 0.0,
            heating_j: (old_kinetic - kinetic) * volume,
            energy_residual_j: 0.0,
        });
    }
    // At no slip, zero duration or zero coupling, retain the original state
    // exactly, including its caloric datum and represented momentum.
    let momentum = if slip == 0.0 || removed == 0.0 {
        initial_momentum
    } else if removed <= 0.5 {
        initial_momentum - (density * slip) * removed
    } else {
        density * (wall + slip * factor)
    };
    let change = momentum - initial_momentum;
    let work_density = wall * change;
    let energy = initial_energy + work_density;
    // Actual accepted impulse determines the slip loss. This avoids a
    // subtraction of boosted bulk kinetic energies and never deposits a second
    // heat source. dK + Q_slip = U_wall * dp for fixed gas amount.
    let heating_j = -change * (slip + 0.5 * (change / density)) * volume;
    let update = Update {
        momentum,
        energy,
        work_j: work_density * volume,
        heating_j,
        energy_residual_j: ((energy - initial_energy) - work_density) * volume,
    };
    if ![
        update.momentum,
        update.energy,
        update.work_j,
        update.heating_j,
        update.energy_residual_j,
    ]
    .into_iter()
    .all(f64::is_finite)
        || update.heating_j < 0.0
    {
        return Err(WallExchangeError::Invalid(
            "moving wall mechanical exchange",
        ));
    }
    Ok(update)
}

#[cfg(test)]
mod tests;
