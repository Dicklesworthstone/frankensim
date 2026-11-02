//! A compact stagnation reservoir discharging through a converging throat.
//! The existing source-owned acceleration and sonic limit remain authoritative.

use super::{inlet_limit, poll};
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
        return Err(WallError::Invalid(
            "stationary nozzle donor or backpressure",
        ));
    }
    if backpressure == donor.pressure {
        return Ok((donor, false));
    }
    // At fixed temperature the original source supplies R*ln(p0/pb), with
    // its cancellation-safe pressure difference. Resolve the pressure DROP,
    // rather than accepting an absolute entropy residual of order Cp, which
    // can swamp the entire driving potential near pressure equilibrium.
    let entropy_drop = gas
        .entropy_change(
            donor.temperature,
            donor.pressure,
            donor.temperature,
            backpressure,
            cancelled,
        )
        .map_err(WallError::Source)?;
    let r = gas.gas_constant();
    let cp0 = gas
        .heat_capacities(donor.temperature, donor.pressure, cancelled)
        .map_err(WallError::Source)?
        .cp;
    if !positive(entropy_drop) || !positive(r) || !positive(cp0) {
        return Err(WallError::CaloricPath(
            "unrepresentable nozzle pressure drop",
        ));
    }
    let logarithmic_drop = entropy_drop / r;
    let pressure_ratio = backpressure / donor.pressure;
    let weak_cooling_bound = donor.temperature * ((donor.pressure - backpressure) / donor.pressure);
    // Positive Cv gives Cp>R, hence T_exit>T0*(pb/p0) and
    // u^2=2*integral(Cp dT) <= 2*R*T0*ln(p0/pb). If the latter is below
    // R*T0*(pb/p0)<c_exit^2, the requested state is certainly subsonic.
    // Common near-equilibrium valves need no preliminary sonic-root solve.
    let maximum_cooling = if 2.0 * logarithmic_drop < pressure_ratio
        && donor.temperature - weak_cooling_bound >= gas.domain().minimum_temperature
    {
        weak_cooling_bound
    } else {
        let (speed, limit, sonic) = inlet_limit(gas, donor, cancelled)?;
        if backpressure <= limit.pressure {
            if !sonic && backpressure < limit.pressure {
                return Err(WallError::CaloricPath("nozzle exceeds cold source path"));
            }
            return Ok((
                Primitive {
                    velocity: speed,
                    ..limit
                },
                sonic,
            ));
        }
        donor.temperature - limit.temperature
    };
    // z = (T0-T)/Delta s_p factors out the weak driving pressure. Averaging
    // the source's exact Cp/T increment over its represented temperature span
    // retains the derivative limit even if T0-Delta T rounds back to T0.
    let mut lo = 0.0;
    let mut hi = maximum_cooling / entropy_drop;
    if !positive(hi) {
        return Err(WallError::CaloricPath(
            "nozzle has no supported cooling interval",
        ));
    }
    let initial = donor.temperature / cp0;
    let mut z = if initial > lo && initial < hi {
        initial
    } else {
        lo.midpoint(hi)
    };
    let tolerance = gas.domain().relative_energy_tolerance + 128.0 * f64::EPSILON;
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let cooling = z * entropy_drop;
        let t = donor.temperature - cooling;
        if !positive(cooling) || t < gas.domain().minimum_temperature {
            return Err(WallError::CaloricPath("nozzle exceeds cold source path"));
        }
        let represented_cooling = donor.temperature - t;
        let entropy_secant = if represented_cooling == 0.0 {
            cp0 / donor.temperature
        } else {
            gas.entropy_change(
                t,
                donor.pressure,
                donor.temperature,
                donor.pressure,
                cancelled,
            )
            .map_err(WallError::Source)?
                / represented_cooling
        };
        let cp = gas
            .heat_capacities(t, backpressure, cancelled)
            .map_err(WallError::Source)?
            .cp;
        let residual = z * entropy_secant - 1.0;
        if !positive(cp) || !positive(entropy_secant) || !residual.is_finite() {
            return Err(WallError::CaloricPath("nozzle isentropic source"));
        }
        if residual.abs() <= tolerance {
            let enthalpy_secant = if represented_cooling == 0.0 {
                cp0
            } else {
                gas.internal_energy_change(t, donor.temperature, cancelled)
                    .map_err(WallError::Source)?
                    / represented_cooling
                    + r
            };
            let kinetic = 2.0 * cooling * enthalpy_secant;
            if !positive(kinetic) {
                return Err(WallError::CaloricPath("nozzle pressure drop is unresolved"));
            }
            let velocity = det::sqrt(kinetic);
            let density = (backpressure / gas.gas_constant()) / t;
            let exit = gas
                .primitive(density, velocity, backpressure, cancelled)
                .map_err(WallError::Source)?;
            if velocity > exit.sound_speed * (1.0 + 128.0 * f64::EPSILON) {
                return Err(WallError::CaloricPath("nonclassical subsonic nozzle path"));
            }
            return Ok((exit, false));
        }
        if residual < 0.0 {
            lo = z;
        } else {
            hi = z;
        }
        let next = z - residual * (t / cp);
        z = if next > lo && next < hi {
            next
        } else {
            lo.midpoint(hi)
        };
    }
    Err(WallError::Iterations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::wall::shared::HeatCapacities;
    use crate::dynamics::{CaloricDomain, Error, FrozenGas};
    use std::cell::Cell;

    // Count calls into the actual source-owning gas. No substitute caloric
    // values or numerical oracle enter the production closure in this probe.
    struct Counted<'a> {
        gas: &'a FrozenGas,
        queries: Cell<usize>,
    }
    impl Counted<'_> {
        fn query(&self) {
            self.queries.set(self.queries.get() + 1);
        }
    }
    impl CaloricLaw for Counted<'_> {
        type Error = Error;
        fn domain(&self) -> CaloricDomain {
            self.gas.domain()
        }
        fn gas_constant(&self) -> f64 {
            self.gas.gas_constant()
        }
        fn heat_capacities(
            &self,
            t: f64,
            p: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<HeatCapacities, Error> {
            self.query();
            self.gas.heat_capacities(t, p, cancelled)
        }
        fn internal_energy_change(
            &self,
            initial: f64,
            final_t: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<f64, Error> {
            self.query();
            self.gas.internal_energy_change(initial, final_t, cancelled)
        }
        fn entropy_change(
            &self,
            initial_t: f64,
            initial_p: f64,
            final_t: f64,
            final_p: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<f64, Error> {
            self.query();
            self.gas
                .entropy_change(initial_t, initial_p, final_t, final_p, cancelled)
        }
        fn sound_speed(
            &self,
            density: f64,
            pressure: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<f64, Error> {
            self.query();
            self.gas.sound_speed(density, pressure, cancelled)
        }
        fn primitive(
            &self,
            density: f64,
            velocity: f64,
            pressure: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Primitive, Error> {
            self.query();
            CaloricLaw::primitive(self.gas, density, velocity, pressure, cancelled)
        }
        fn next_seam(
            &self,
            start: f64,
            end: f64,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<f64, Error> {
            self.query();
            self.gas.next_seam(start, end, cancelled)
        }
    }

    #[test]
    fn weak_pressure_solve_avoids_the_preliminary_sonic_source_queries() {
        let source = crate::dynamics::tests::gas(0.4, 0.0);
        let donor = source.from_rho_u_p(1.0, 0.0, 10.0).unwrap();
        let counted = Counted {
            gas: &source,
            queries: Cell::new(0),
        };
        let (speed, limit, sonic) = inlet_limit(&counted, donor, &mut || false).unwrap();
        assert!(sonic);
        let preliminary_queries = counted.queries.replace(0);
        let (exit, choked) = state(&counted, donor, 9.9999, &mut || false).unwrap();
        let complete_queries = counted.queries.get();
        assert!(!choked && exit.velocity > 0.0);
        assert!(
            complete_queries < preliminary_queries,
            "full weak solve={complete_queries}, former mandatory preliminary solve={preliminary_queries}"
        );
        eprintln!(
            "source queries: full weak solve={complete_queries}, preliminary sonic solve={preliminary_queries}"
        );
        // Strong-drop choking still returns the original sonic owner's state.
        assert_eq!(
            state(&counted, donor, 0.5 * limit.pressure, &mut || false).unwrap(),
            (
                Primitive {
                    velocity: speed,
                    ..limit
                },
                true
            )
        );
    }
}
