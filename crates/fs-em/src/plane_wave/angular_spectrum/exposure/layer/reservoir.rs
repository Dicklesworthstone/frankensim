//! Original instantaneous and Debye reservoirs for settled harmonic fields.
//!
//! tau P_dot + P = Delta_epsilon E implies
//! E dot P_dot = (tau/Delta_epsilon)|P_dot|²
//!               + d(|P|²/(2 Delta_epsilon))/dt.
//! Heat is formed from each reservoir's actual phased current, not from the
//! magnitudes of per-frequency complex permittivity or boundary flux loss.

use super::*;
use crate::plane_wave::{DielectricResponse, ohmic_conductivity, shared_ohmic_conductivity};
use fs_material::field_storage::{
    dielectric::DebyePole, linear_electric_permittivity, linear_magnetic_permeability,
};

pub(super) struct LossLaw<'m> {
    pub epsilon_infinity: f64,
    pub permeability: f64,
    pub conductivity: f64,
    pub poles: &'m [DebyePole],
}

impl<'m> LossLaw<'m> {
    pub fn new(response: DielectricResponse<'m>) -> Result<Self, AngularSpectrumError> {
        let (epsilon_infinity, conductivity, poles): (f64, f64, &'m [DebyePole]) = match response {
            DielectricResponse::Lossless(material) => (
                linear_electric_permittivity(material).map_err(ReflectionError::from)?,
                0.0,
                &[],
            ),
            DielectricResponse::Ohmic(conductor) => (
                linear_electric_permittivity(conductor.material()).map_err(ReflectionError::from)?,
                ohmic_conductivity(conductor)?,
                &[],
            ),
            DielectricResponse::Debye(dielectric) => (
                dielectric.high_frequency_permittivity_f_m(),
                0.0,
                dielectric.poles(),
            ),
            DielectricResponse::DebyeWithOhmic { dielectric, conductor } => (
                dielectric.high_frequency_permittivity_f_m(),
                shared_ohmic_conductivity(dielectric, conductor)?,
                dielectric.poles(),
            ),
            DielectricResponse::OpticalConductor(_) | DielectricResponse::OpticalMaterial(_) => {
                return Err(bad(
                    "finite-window film heat needs an original causal material reservoir",
                ));
            }
        };
        let permeability = linear_magnetic_permeability(response.material())
            .map_err(ReflectionError::from)?;
        Ok(Self { epsilon_infinity, permeability, conductivity, poles })
    }
}

#[derive(Clone, Copy)]
pub(super) struct PoleTransfer {
    // P/sqrt(Delta_epsilon), so its squared real field stores twice the energy.
    pub storage: C64,
    // sqrt(tau/Delta_epsilon) P_dot, whose squared real field is heat.
    pub heat: C64,
}

impl PoleTransfer {
    pub fn new(pole: DebyePole, omega: f64) -> Result<Self, AngularSpectrumError> {
        let delta = pole.permittivity_increment_f_m;
        let tau = pole.relaxation_time_s;
        if [delta, tau, omega].iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(bad("settled Debye response needs positive source poles and frequency"));
        }
        let root_delta = det::sqrt(delta);
        let root_tau = det::sqrt(tau);
        let x = omega * tau;
        let (storage, heat) = if x <= 1.0 {
            let denominator = 1.0 + x * x;
            let storage = C64::new(
                root_delta / denominator,
                -multiply_three(root_delta, omega, tau) / denominator,
            );
            let heat_im = multiply_three(root_delta, root_tau, omega) / denominator;
            (storage, C64::new(heat_im * x, heat_im))
        } else {
            // This remains finite when omega*tau overflows. Avoid forming a
            // tiny reciprocal before its multiplication by a large source
            // amplitude, which can erase a representable quadrature component.
            let inverse = (1.0 / omega) / tau;
            let denominator = 1.0 + inverse * inverse;
            let storage_im = -divide_product(root_delta, omega, tau) / denominator;
            let heat_re = (root_delta / root_tau) / denominator;
            (
                C64::new(-storage_im * inverse, storage_im),
                C64::new(heat_re, divide_product(heat_re, omega, tau)),
            )
        };
        if !finite(storage) || !finite(heat) {
            return Err(bad("settled Debye transfer is unrepresentable"));
        }
        Ok(Self { storage, heat })
    }
}

// All operands are positive source coefficients. Combining the largest and
// smallest first avoids spurious overflow/underflow in a representable product.
fn multiply_three(a: f64, b: f64, c: f64) -> f64 {
    let mut values = [a, b, c];
    values.sort_unstable_by(f64::total_cmp);
    (values[0] * values[2]) * values[1]
}

fn divide_product(value: f64, a: f64, b: f64) -> f64 {
    let (small, large) = (a.min(b), a.max(b));
    let first = value / small;
    if first.is_finite() { first / large } else { (value / large) / small }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pole(delta: f64, tau: f64) -> DebyePole {
        DebyePole {
            permittivity_increment_f_m: delta,
            relaxation_time_s: tau,
            capacitance_f: 1.0,
            resistance_ohm: tau,
        }
    }

    #[test]
    fn original_pole_transfer_retains_ode_phase_and_positive_heat() {
        let pole = pole(0.7, 0.4);
        for omega in [1e-8, 2.0, 1e8] {
            let response = PoleTransfer::new(pole, omega).unwrap();
            let p = response.storage.scale(det::sqrt(0.7));
            let current = C64::new(-p.im, p.re).scale(omega);
            let ode = p + current.scale(0.4);
            assert!((ode.re / 0.7 - 1.0).abs() < 5e-15);
            assert!(ode.im.abs() < 2e-16);
            let heat = response.heat.norm_sq();
            let expected = (current.norm_sq() / 0.7) * 0.4;
            assert!((heat / expected - 1.0).abs() < 5e-15);
        }
    }

    #[test]
    fn extreme_frequency_products_do_not_erase_representable_reservoir_components() {
        let low = PoleTransfer::new(pole(1e300, 1e-200), 1e-200).unwrap();
        assert!((low.storage.im / -1e-250 - 1.0).abs() < 2e-14);
        assert!((low.heat.im / 1e-150 - 1.0).abs() < 2e-14);
        let high = PoleTransfer::new(pole(1e300, 1e200), 1e200).unwrap();
        assert!((high.storage.im / -1e-250 - 1.0).abs() < 2e-14);
        assert!((high.heat.re / 1e50 - 1.0).abs() < 2e-14);
    }
}
