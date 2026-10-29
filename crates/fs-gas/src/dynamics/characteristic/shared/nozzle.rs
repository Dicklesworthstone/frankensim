//! A compact stagnation reservoir discharging through a converging throat.
//! The existing source-owned acceleration and sonic limit remain authoritative.

use super::{enthalpy_difference, inlet_limit, poll};
use crate::dynamics::wall::shared::{CaloricLaw, WallError};
use crate::dynamics::{Primitive, positive};
use fs_math::det;

/// Positive exit speed is directed from the donor toward the receiver. This
/// closure is for two compact volumes, not a duct characteristic boundary.
pub(crate) fn state<L: CaloricLaw>(
    gas: &L,
    donor: Primitive,
    backpressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Primitive, bool), WallError<L::Error>> {
    poll(cancelled)?;
    if !positive(backpressure) || backpressure > donor.pressure || donor.velocity != 0.0 {
        return Err(WallError::Invalid("stationary nozzle donor or backpressure"));
    }
    if backpressure == donor.pressure {
        return Ok((donor, false));
    }
    let (speed, limit, sonic) = inlet_limit(gas, donor, cancelled)?;
    if backpressure <= limit.pressure {
        if !sonic && backpressure < limit.pressure {
            return Err(WallError::CaloricPath("nozzle exceeds cold source path"));
        }
        return Ok((Primitive { velocity: speed, ..limit }, sonic));
    }
    // Solve the original isentrope at the actual backpressure. The derivative
    // is Cp(T)/T, so no finite-difference Jacobian or frozen gamma is needed.
    let mut lo = limit.temperature;
    let mut hi = donor.temperature;
    let mut t = lo.midpoint(hi);
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let entropy = gas.entropy_change(
            donor.temperature, donor.pressure, t, backpressure, cancelled,
        ).map_err(WallError::Source)?;
        let cp = gas.heat_capacities(t, backpressure, cancelled)
            .map_err(WallError::Source)?.cp;
        if !positive(cp) || !entropy.is_finite() {
            return Err(WallError::CaloricPath("nozzle isentropic source"));
        }
        let tolerance = gas.domain().absolute_energy_tolerance / donor.temperature
            + gas.domain().relative_energy_tolerance * cp;
        if entropy.abs() <= tolerance {
            let kinetic = 2.0 * enthalpy_difference(gas, t, donor.temperature, cancelled)?;
            if !positive(kinetic) {
                return Err(WallError::CaloricPath("nozzle pressure drop is unresolved"));
            }
            let velocity = det::sqrt(kinetic);
            let density = (backpressure / gas.gas_constant()) / t;
            let exit = gas.primitive(density, velocity, backpressure, cancelled)
                .map_err(WallError::Source)?;
            if velocity > exit.sound_speed * (1.0 + 128.0 * f64::EPSILON) {
                return Err(WallError::CaloricPath("nonclassical subsonic nozzle path"));
            }
            return Ok((exit, false));
        }
        if entropy < 0.0 { lo = t; } else { hi = t; }
        let next = t - entropy * (t / cp);
        t = if next > lo && next < hi { next } else { lo.midpoint(hi) };
    }
    Err(WallError::Iterations)
}
