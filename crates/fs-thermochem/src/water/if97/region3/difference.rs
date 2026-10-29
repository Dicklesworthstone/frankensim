//! Cancellation-resistant changes between original Region 3 source states.
//!
//! The same Helmholtz coefficients as [`super::If97DenseState`] are differenced
//! analytically. Powers retain the input density/temperature increments before
//! reduced coordinates round; the logarithm uses the stable entropy log-ratio
//! primitive. Thus small heat or
//! pressure changes do not subtract two megajoule-scale reference energies.
//! The endpoints have already passed the source's domain and stability checks.
//! This operation asserts no admissible path between them: in particular a
//! liquid/vapor endpoint pair is not a homogeneous intermediate-density state.

use super::{CRITICAL_DENSITY_KG_M3, If97DenseState, LOG_COEFFICIENT, TERMS};
use crate::water::MassSpecificVolume;
use crate::water::if97::{CRITICAL_TEMPERATURE_K, R_SPECIFIC_J_KG_K};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use fs_qty::Pressure;

/// Original-source endpoint changes, in the direction `to - from`.
///
/// Stable energy differences are not replacements for the retained absolute
/// source values. A consumer that stores absolute energy must also audit their
/// representability, its pressure closure, and its physical energy balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseStateDifference {
    internal_energy_j_kg: f64,
    enthalpy_j_kg: f64,
    entropy_j_kg_k: f64,
    volume_m3_kg: f64,
    pressure_pa: f64,
    flow_work_j_kg: f64,
}

impl DenseStateDifference {
    pub const fn specific_internal_energy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    pub const fn specific_enthalpy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    pub const fn specific_entropy_change(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    pub const fn specific_volume_change(self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.volume_m3_kg)
    }
    pub const fn pressure_change(self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
    /// Change in the original `p/rho`, already included in the enthalpy change.
    pub const fn flow_work_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.flow_work_j_kg)
    }
}

/// Evaluate six endpoint changes with fixed work, no allocation and no solve.
///
/// Exact zero and reversal follow from a canonical endpoint ordering. The
/// arithmetic approximates the original binary64 coefficient model, not a
/// separate fitted equation or a quantified physical uncertainty interval.
pub fn difference(from: If97DenseState, to: If97DenseState) -> DenseStateDifference {
    let (ra, ta) = (from.density().value(), from.temperature().value());
    let (rb, tb) = (to.density().value(), to.temperature().value());
    let (a, b, sign) = if ta < tb || (ta == tb && ra <= rb) {
        (from, to, 1.0)
    } else {
        (to, from, -1.0)
    };
    let (ra, ta) = (a.density().value(), a.temperature().value());
    let (rb, tb) = (b.density().value(), b.temperature().value());
    let drho = rb - ra;
    let dt = tb - ta;
    let delta = Powers::<13>::new(
        ra / CRITICAL_DENSITY_KG_M3,
        rb / CRITICAL_DENSITY_KG_M3,
        drho / CRITICAL_DENSITY_KG_M3,
    );
    let tau_a = CRITICAL_TEMPERATURE_K / ta;
    let tau = Powers::<27>::new(tau_a, CRITICAL_TEMPERATURE_K / tb, -tau_a * (dt / tb));
    let mut u = Sum::default();
    let mut pv = Sum::default();
    let mut p = Sum::default();
    let mut s = Sum::default();
    pv.add(LOG_COEFFICIENT * dt);
    p.add(LOG_COEFFICIENT * (tb * delta.change[1] + delta.a[1] * dt));
    s.add(-LOG_COEFFICIENT * crate::entropy::log_ratio(rb, ra));
    for &(i, j, n) in &TERMS {
        let caloric = if j == 0 {
            tb * delta.change[i] + delta.a[i] * dt
        } else {
            CRITICAL_TEMPERATURE_K * product_change(&delta, i, &tau, j - 1)
        };
        if j != 0 {
            u.add(n * j as f64 * caloric);
        }
        if i != 0 {
            pv.add(n * i as f64 * caloric);
            let pressure = if j == 0 {
                tb * delta.change[i + 1] + delta.a[i + 1] * dt
            } else {
                CRITICAL_TEMPERATURE_K * product_change(&delta, i + 1, &tau, j - 1)
            };
            p.add(n * i as f64 * pressure);
        }
        if j != 1 {
            s.add(n * (j as f64 - 1.0) * product_change(&delta, i, &tau, j));
        }
    }
    let scale = sign * R_SPECIFIC_J_KG_K;
    let du = scale * u.value();
    let dpv = scale * pv.value();
    DenseStateDifference {
        internal_energy_j_kg: du,
        enthalpy_j_kg: du + dpv,
        entropy_j_kg_k: scale * s.value(),
        volume_m3_kg: -sign * (drho / ra) / rb,
        pressure_pa: scale * CRITICAL_DENSITY_KG_M3 * p.value(),
        flow_work_j_kg: dpv,
    }
}

struct Powers<const N: usize> {
    a: [f64; N],
    b: [f64; N],
    change: [f64; N],
}
impl<const N: usize> Powers<N> {
    fn new(a: f64, b: f64, change: f64) -> Self {
        let mut powers = Self {
            a: [1.0; N],
            b: [1.0; N],
            change: [0.0; N],
        };
        for i in 1..N {
            powers.a[i] = powers.a[i - 1] * a;
            powers.b[i] = powers.b[i - 1] * b;
            powers.change[i] = b * powers.change[i - 1] + powers.a[i - 1] * change;
        }
        powers
    }
}

fn product_change<const A: usize, const B: usize>(
    x: &Powers<A>,
    i: usize,
    y: &Powers<B>,
    j: usize,
) -> f64 {
    x.b[i] * y.change[j] + y.a[j] * x.change[i]
}

#[derive(Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }
    fn value(self) -> f64 {
        self.total + self.correction
    }
}

#[cfg(test)]
mod tests;
