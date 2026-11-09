//! Coherent counter-propagating acoustic modes in a finite homogeneous segment.
//! Endpoint amplitudes use only decaying propagation, including in opaque cells.

use super::{BulkWaveAbsorption, BulkWaveField, admit, finite, invalid, propagate};
use crate::MaterialError;
use fs_math::{c64::C64, dd::Dd, det};

/// Two coherent native acoustic modes over x in [0,length]. Both supplied
/// velocities are signed along +x: u=uf+ub, while p, traction and T' are
/// proportional to uf-ub. Their phases matter; summing separate traveling-wave
/// losses would discard standing-wave interference. No reflection coefficient,
/// boundary condition or independent thermal mode is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkCoherentWaveSegment {
    absorption: BulkWaveAbsorption,
    length_m: f64,
    forward_at_front: C64,
    backward_at_back: C64,
}

/// Coherent cell absorption with signed +x acoustic powers at its two faces.
/// Left minus right is the power lost inside the cell, including energy entering
/// from the right. Only viscous_heat_w is direct caloric heating; conductive
/// entropy and its T0-weighted availability loss retain their distinct meaning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkCoherentCellAbsorption {
    pub left_acoustic_power_w: f64,
    pub right_acoustic_power_w: f64,
    pub viscous_heat_w: f64,
    pub conductive_entropy_w_k: f64,
    pub thermal_availability_loss_w: f64,
    pub acoustic_power_loss_w: f64,
    /// Left minus right minus independently integrated constitutive losses.
    pub balance_residual_w: f64,
}

impl BulkWaveAbsorption {
    /// Bind coherent peak velocity phasors at their respective incident faces.
    /// The forward wave is specified at x=0 and the backward wave at x=length;
    /// both velocities use the same +x sign convention. Length is finite and
    /// nonnegative. Queries allocate nothing and never reevaluate the material.
    pub fn coherent_segment(
        &self,
        length_m: f64,
        forward_at_front_m_s: C64,
        backward_at_back_m_s: C64,
    ) -> Result<BulkCoherentWaveSegment, MaterialError> {
        admit(forward_at_front_m_s, length_m)?;
        admit(backward_at_back_m_s, length_m)?;
        Ok(BulkCoherentWaveSegment {
            absorption: *self,
            length_m,
            forward_at_front: forward_at_front_m_s,
            backward_at_back: backward_at_back_m_s,
        })
    }
}

impl BulkCoherentWaveSegment {
    /// Original frequency, wave and source-coherent equilibrium.
    pub const fn absorption(&self) -> &BulkWaveAbsorption {
        &self.absorption
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }

    /// Observe coherent fields, signed power flux and nonnegative local losses.
    /// Pressure/temperature nodes and velocity nodes dissipate differently:
    /// grad(u)=ik(uf-ub), while grad(T')=ik(T'/uf)(uf+ub).
    pub fn at_distance(&self, distance_m: f64) -> Result<BulkWaveField, MaterialError> {
        self.admit_distance(distance_m)?;
        let k = self.absorption.wave.wavenumber_rad_m();
        self.field(
            propagate(self.forward_at_front, k, distance_m)?,
            propagate(self.backward_at_back, k, self.length_m - distance_m)?,
        )
    }

    fn field(&self, forward: C64, backward: C64) -> Result<BulkWaveField, MaterialError> {
        let prepared = &self.absorption;
        let k = prepared.wave.wavenumber_rad_m();
        let ik = C64::new(-k.im, k.re);
        let t0 = prepared.equilibrium.temperature_k();
        let kappa = prepared
            .equilibrium
            .absorption()
            .state()
            .transport()
            .thermal_conductivity_w_m_k();
        let velocity = forward + backward;
        let difference = forward - backward;
        let pressure = prepared.wave.pressure_impedance_pa_s_m() * difference;
        let traction = prepared.wave.traction_impedance_pa_s_m() * difference;
        let temperature = prepared.temperature_over_velocity * difference;
        let gradient_t = ik * (prepared.temperature_over_velocity * velocity);
        let heat_flux = -gradient_t.scale(kappa);
        let viscous = (ik * difference)
            .scale(det::sqrt(0.5 * prepared.longitudinal_viscosity))
            .norm_sq();
        let conductive = gradient_t.scale(det::sqrt(0.5 * kappa / t0)).norm_sq();
        let entropy = conductive / t0;
        let mechanical = 0.5 * (traction * velocity.conj()).re;
        let thermal = (0.5 / t0) * (temperature * heat_flux.conj()).re;
        if ![velocity, pressure, traction, temperature, heat_flux]
            .into_iter()
            .all(finite)
            || ![
                viscous,
                conductive,
                entropy,
                mechanical,
                thermal,
                mechanical + thermal,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "coherent bulk fields or physical work are unrepresentable",
            ));
        }
        Ok(BulkWaveField {
            velocity_m_s: velocity,
            pressure_pa: pressure,
            normal_traction_pa: traction,
            temperature_perturbation_k: temperature,
            conductive_heat_flux_w_m2: heat_flux,
            mechanical_power_flux_w_m2: mechanical,
            thermal_availability_flux_w_m2: thermal,
            viscous_heating_w_m3: viscous,
            conductive_entropy_w_m3_k: entropy,
            thermal_availability_loss_w_m3: conductive,
        })
    }

    /// Integrate coherent constitutive losses exactly within the admitted
    /// constant-medium model. A positive Gram form retains interference without
    /// subtracting nearly equal powers or growing exponentials. Stable small-cell
    /// coefficients resolve loss near standing-wave nodes below ordinary d-c
    /// roundoff. The explicit physical area is applied once. No quadrature runs.
    pub fn absorption_between(
        &self,
        start_m: f64,
        end_m: f64,
        area_m2: f64,
    ) -> Result<BulkCoherentCellAbsorption, MaterialError> {
        self.admit_distance(start_m)?;
        self.admit_distance(end_m)?;
        if end_m < start_m || !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "coherent absorption needs ordered depths and positive finite area",
            ));
        }
        let prepared = &self.absorption;
        let k = prepared.wave.wavenumber_rad_m();
        let width = end_m - start_m;
        let front = propagate(self.forward_at_front, k, start_m)?;
        let back = propagate(self.backward_at_back, k, self.length_m - end_m)?;
        let left = self.field(front, propagate(back, k, width)?)?;
        let right = self.field(propagate(front, k, width)?, back)?;
        let (plus, minus) = gram_weights(k, width)?;
        let ik = C64::new(-k.im, k.re);
        let sum_gradient = ik * (front + back);
        let difference_gradient = ik * (front - back);
        let t0 = prepared.equilibrium.temperature_k();
        let kappa = prepared
            .equilibrium
            .absorption()
            .state()
            .transport()
            .thermal_conductivity_w_m_k();
        let viscous_scale = det::sqrt(0.5 * prepared.longitudinal_viscosity);
        let thermal_scale = det::sqrt(0.5 * kappa / t0);
        let viscous = positive_form(
            difference_gradient.scale(viscous_scale),
            sum_gradient.scale(viscous_scale),
            plus,
            minus,
        ) * area_m2;
        let conductive = positive_form(
            (prepared.temperature_over_velocity * sum_gradient).scale(thermal_scale),
            (prepared.temperature_over_velocity * difference_gradient).scale(thermal_scale),
            plus,
            minus,
        ) * area_m2;
        let entropy = conductive / t0;
        let absorbed = viscous + conductive;
        let left_power = left.acoustic_power_flux_w_m2() * area_m2;
        let right_power = right.acoustic_power_flux_w_m2() * area_m2;
        let residual =
            (Dd::from_f64(left_power) - Dd::from_f64(right_power) - Dd::from_f64(absorbed))
                .to_f64();
        if ![
            viscous,
            conductive,
            entropy,
            absorbed,
            left_power,
            right_power,
            residual,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid("coherent bulk cell absorption is unrepresentable"));
        }
        Ok(BulkCoherentCellAbsorption {
            left_acoustic_power_w: left_power,
            right_acoustic_power_w: right_power,
            viscous_heat_w: viscous,
            conductive_entropy_w_k: entropy,
            thermal_availability_loss_w: conductive,
            acoustic_power_loss_w: absorbed,
            balance_residual_w: residual,
        })
    }

    fn admit_distance(&self, distance: f64) -> Result<(), MaterialError> {
        if distance.is_finite() && distance >= 0.0 && distance <= self.length_m {
            Ok(())
        } else {
            Err(invalid(
                "coherent bulk observation must lie inside its finite segment",
            ))
        }
    }
}

fn positive_form(a: C64, b: C64, plus: f64, minus: f64) -> f64 {
    a.scale(det::sqrt(plus)).norm_sq() + b.scale(det::sqrt(minus)).norm_sq()
}

/// With Fa at the left and Bb at the right, the integral of |uf +/- ub|²
/// is (d+c)|Fa +/- Bb|²/2 + (d-c)|Fa -/+ Bb|²/2. Here
/// d=(1-exp(-2 alpha w))/(2 alpha), c=w exp(-alpha w) sinc(Re(k) w).
fn gram_weights(k: C64, width: f64) -> Result<(f64, f64), MaterialError> {
    if width == 0.0 {
        return Ok((0.0, 0.0));
    }
    let a = k.im * width;
    let decay = det::exp(-a);
    let d = if a == 0.0 {
        width
    } else if a < 0.5 {
        width * (-det::expm1(-2.0 * a) / (2.0 * a))
    } else {
        (-0.5 * det::expm1(-2.0 * a)) / k.im
    };
    let b = k.re * width;
    let (plus, minus) = if a < 0.5 && b < 0.5 {
        // d-c = w exp(-a) [(sinhc(a)-1) + (1-sinc(b))]. Each series
        // starts at order x², avoiding catastrophic cancellation at a node.
        // Terms through x^14 leave <2e-19 absolute remainder for |x|<0.5.
        let a2 = a * a;
        let b2 = b * b;
        let sinhc_minus_one = a2
            * (1.0 / 6.0
                + a2 * (1.0 / 120.0
                    + a2 * (1.0 / 5040.0
                        + a2 * (1.0 / 362880.0
                            + a2 * (1.0 / 39916800.0
                                + a2 * (1.0 / 6227020800.0 + a2 / 1307674368000.0))))));
        let one_minus_sinc = b2
            * (1.0 / 6.0
                + b2 * (-1.0 / 120.0
                    + b2 * (1.0 / 5040.0
                        + b2 * (-1.0 / 362880.0
                            + b2 * (1.0 / 39916800.0
                                + b2 * (-1.0 / 6227020800.0 + b2 / 1307674368000.0))))));
        let minus = ((0.5 * width) * decay) * (sinhc_minus_one + one_minus_sinc);
        (d - minus, minus)
    } else {
        let cross = if decay == 0.0 {
            0.0
        } else if !b.is_finite() {
            return Err(invalid("surviving coherent wave phase is unrepresentable"));
        } else {
            (width * decay) * if b == 0.0 { 1.0 } else { det::sin(b) / b }
        };
        (0.5 * d + 0.5 * cross, 0.5 * d - 0.5 * cross)
    };
    if [plus, minus].into_iter().all(|v| v.is_finite() && v >= 0.0) {
        Ok((plus, minus))
    } else {
        Err(invalid(
            "coherent absorption Gram coefficients are unrepresentable",
        ))
    }
}

#[cfg(test)]
mod tests;
