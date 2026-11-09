//! Passive scalar termination of the original pressure/velocity two-port.
//!
//! Use z*(1-S22)+r*(1+S22), with scaled backing impedance z and reference
//! impedance r. Forming a rounded reflection first can lose the small wall
//! motion of a stiff finite support. Rear pressure and velocity are independent
//! retained coefficients; neither is recovered by subtracting opposite waves.

use super::super::{Result, divided, finite, from_log, product_ratio};
use super::{C64, Dd, linear};
use fs_math::det;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Termination {
    pub(crate) incoming_right: C64,
    pub(crate) reflected: C64,
    pub(crate) pressure: C64,
    pub(crate) velocity: C64,
    resistance: f64,
}

impl Termination {
    /// `None` is an exact clamp. `Some(0)` is the exact traction-free limit.
    /// Finite impedances use p=Z*v in the original exp(+i omega t) convention.
    pub(crate) fn new(scattering: [C64; 4], reference: f64, backing: Option<C64>) -> Result<Self> {
        if !reference.is_finite()
            || reference <= 0.0
            || scattering
                .into_iter()
                .any(|v| !finite(v) || !v.abs().is_finite())
            || backing.is_some_and(|z| !finite(z) || z.re < 0.0)
        {
            return Err("supported layer needs a finite passive backing and positive reference");
        }
        let (z, r, scale, resistance) = if let Some(backing) = backing {
            let scale = reference.max(backing.re).max(backing.im.abs());
            let r = reference / scale;
            if r == 0.0 {
                return Err("backing/reference impedance contrast is unrepresentable");
            }
            (divided(backing, scale), r, scale, backing.re)
        } else {
            (C64::ONE, 0.0, reference, 0.0)
        };
        let [sr, si, zr, zi, r] =
            [scattering[3].re, scattering[3].im, z.re, z.im, r].map(Dd::from_f64);
        let minus = Dd::ONE - sr;
        let plus = Dd::ONE + sr;
        let denominator = C64::new(
            (zr * minus + zi * si + r * plus).to_f64(),
            (zi * minus - zr * si + r * si).to_f64(),
        );
        if !finite(denominator) || denominator == C64::ZERO {
            return Err("supported viscoelastic round trip is singular or unrepresentable");
        }
        let amplitude = scattering[2] / denominator;
        let incoming_right = linear([z, C64::from_re(-r.to_f64())], [amplitude; 2]);
        let reflected = linear([scattering[0], scattering[1]], [C64::ONE, incoming_right]);
        let pressure = (z * amplitude).scale(2.0);
        let velocity = if backing.is_none() {
            C64::ZERO
        } else {
            let component = |value: f64| -> Result<f64> {
                Ok(product_ratio(2.0, 1.0, value.abs(), scale)?.copysign(value))
            };
            C64::new(component(amplitude.re)?, component(amplitude.im)?)
        };
        if [incoming_right, reflected, pressure, velocity]
            .into_iter()
            .any(|v| !finite(v) || !v.abs().is_finite())
        {
            return Err("supported viscoelastic boundary coefficients are unrepresentable");
        }
        Ok(Self {
            incoming_right,
            reflected,
            pressure,
            velocity,
            resistance,
        })
    }

    /// Positive support work, independent of the reflected-power deficit.
    pub(crate) fn absorbed_power(self, velocity: C64, area: f64) -> Result<f64> {
        if self.resistance == 0.0 {
            return Ok(0.0);
        }
        from_log(
            det::ln(self.resistance) + 2.0 * det::ln(velocity.abs()) + det::ln(area) - det::ln(2.0),
        )
    }
}

pub(crate) fn traveling_power(pressure: C64, impedance: f64, area: f64) -> Result<f64> {
    from_log(2.0 * det::ln(pressure.abs()) + det::ln(area) - det::ln(2.0) - det::ln(impedance))
}
