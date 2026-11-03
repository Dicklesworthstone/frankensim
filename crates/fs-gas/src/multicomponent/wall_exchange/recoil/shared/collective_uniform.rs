//! Exact nonlinear recoil for an invariant common-velocity gas group.
//!
//! Equal frozen acceleration laws and equal entering velocities keep every
//! active gas velocity equal. Their summed mass therefore follows the original
//! exact two-body map; each cell receives its own mass share of work and heat.

use super::{
    RecoilLaw, SharedRecoilCell, SlidingSupport, Stage, Sum, WallExchangeError, mechanics, poll,
    positive,
};

pub(super) fn advance(
    cells: &[SharedRecoilCell],
    states: &mut [Stage],
    dt: f64,
    support: SlidingSupport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<SlidingSupport>, WallExchangeError> {
    let mut group = None;
    let mut total_mass = Sum::default();
    for cell in cells {
        poll(cancelled)?;
        if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
            continue;
        }
        // Linear groups already have the more general exact collective mode.
        // No tolerance merges different physical states or resistance laws.
        if cell.quadratic_per_m == 0.0 {
            return Ok(None);
        }
        let descriptor = (
            cell.momentum_kg_m2_s / cell.density_kg_m3,
            cell.linear_rate_per_s,
            cell.quadratic_per_m,
        );
        if group.is_some_and(|previous| previous != descriptor) {
            return Ok(None);
        }
        group = Some(descriptor);
        total_mass.add(cell.density_kg_m3 * cell.volume_m3);
    }
    let Some((velocity, linear, quadratic)) = group else {
        return Ok(None);
    };
    let mass = total_mass.value();
    if !positive(mass) {
        // Retain the general owner when one aggregate mass is unrepresentable.
        return Ok(None);
    }
    // This unit-density mechanical coordinate makes the map's momentum equal
    // to velocity and reproduces the actual summed inertial mass through its
    // volume parameter. No synthetic volume or gas state is ever published.
    let motion = mechanics::advance(
        velocity,
        1.0,
        velocity,
        mass,
        dt,
        support,
        RecoilLaw {
            linear_rate_per_s: linear,
            quadratic_per_m: quadratic,
            gas_heat_fraction: 1.0,
        },
    )?;
    let delta = motion.momentum - velocity;
    let support_delta = -motion.gas_impulse_ns / support.mass_kg;
    let raw_work = motion.gas_impulse_ns * support_delta;
    let frame_work = if raw_work.is_finite() {
        0.5 * raw_work
    } else {
        (0.5 * motion.gas_impulse_ns) * support_delta
    };
    let slip = velocity - support.velocity_m_s;
    for (cell, state) in cells.iter().zip(states) {
        poll(cancelled)?;
        if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
            continue;
        }
        let momentum = if delta == 0.0 {
            cell.momentum_kg_m2_s
        } else if delta.abs() <= 0.5 * slip.abs() {
            cell.momentum_kg_m2_s + cell.density_kg_m3 * delta
        } else {
            cell.density_kg_m3 * motion.momentum
        };
        let cell_mass = cell.density_kg_m3 * cell.volume_m3;
        let heat = mass_share(motion.dissipation_j, cell_mass, mass);
        let work = mass_share(frame_work, cell_mass, mass);
        if ![momentum, heat, work].into_iter().all(f64::is_finite) || heat < 0.0 {
            return Err(WallExchangeError::Invalid(
                "collective uniform recoil endpoint",
            ));
        }
        state.momentum = momentum;
        state.heat.add(heat);
        state.work.add(work);
    }
    // The outer owner restores U0 times each ACTUAL cell impulse and audits
    // actual gas/support endpoints. No heat or energy residual is repaired.
    poll(cancelled)?;
    Ok(Some(SlidingSupport {
        velocity_m_s: motion.support_velocity_m_s,
        ..support
    }))
}

fn mass_share(total: f64, mass: f64, combined_mass: f64) -> f64 {
    let fraction = mass / combined_mass;
    let specific = total / combined_mass;
    if fraction < f64::MIN_POSITIVE && specific.is_finite() {
        // A very small mass share can still carry a resolved part of a large
        // extensive exchange; divide that exchange before multiplying mass.
        // If that quotient overflows, the nonzero mass fraction retains the
        // finite extensive share without materializing a specific energy.
        mass * specific
    } else {
        fraction * total
    }
}
