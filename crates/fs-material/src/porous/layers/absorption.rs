//! Physical absorption in each original equivalent-fluid layer.
//!
//! With exp(+i omega t), the existing complex density rho and bulk modulus K
//! imply the local mean loss
//! `q = omega/2 * (-Im(rho) |v|² - Im(1/K) |p|²)` [W/m³].
//! This follows by differentiating `I_n = Re(p conj(v_n))/2` in the original
//! momentum/continuity equations. At oblique incidence the tangential velocity
//! `v_t = k_t p/(omega rho)` MUST participate in the viscous loss. Equivalent
//! fluid/source conventions: COMSOL, About the Poroacoustics Models,
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.182.html>.
//!
//! The SAME bounded R/T propagation supplies both faces of every layer. Two
//! inward-decaying waves give analytic integrals of |p|² and |v_n|² without
//! growing exponentials. A positive Gram decomposition retains thin-layer
//! standing-wave losses, and logarithmic amplitude scaling retains weak heat
//! below the resolution of `1-R-T`. No power residual is redistributed as heat.
//! Rigid frame, infinite planar layers and the original source support remain;
//! the viscous/thermal mechanisms are not a fluid/solid temperature partition.

use super::transmission::{LogAmplitude, log_scaled};
use super::{
    LayerTransfer, MAX_LAYERS, RigidBackedLayers, TransmittingLayers, divided, finite,
    magnitude_scale,
};
use crate::MaterialError;
use fs_math::{c64::C64, dd::Dd, det};

/// Caller-owned allowance for layer propagation and its analytic loss integral.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LayerAbsorptionBudget {
    /// Declared layers, including zero-thickness identities; at most MAX_LAYERS.
    /// Zero is valid for an empty stack.
    pub max_layers: usize,
}

/// Work in one declared layer, in the original exposed-to-rear order.
/// The fractions use incident normal power through the physical face. The
/// viscous term includes BOTH normal and tangential particle velocity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LayerAbsorbedPower {
    pub viscous_absorbed_fraction: f64,
    pub thermal_absorbed_fraction: f64,
    pub absorbed_power_w: f64,
}
impl LayerAbsorbedPower {
    pub fn absorbed_fraction(&self) -> f64 {
        self.viscous_absorbed_fraction + self.thermal_absorbed_fraction
    }
}

/// Independent R/T and positive local loss accounts at one incidence angle.
/// No coefficient is clipped or renormalized. Zero incident watts retain the
/// coefficients while all powers vanish. Preparation and this query are not
/// a transient thermal feedback solve or an empirical heating calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerAbsorption {
    pub incident_power_w: f64,
    pub reflected_power_w: f64,
    pub transmitted_power_w: f64,
    pub absorbed_power_w: f64,
    pub reflected_fraction: f64,
    pub transmitted_fraction: f64,
    pub absorbed_fraction: f64,
    pub balance_residual_w: f64,
    layers: [LayerAbsorbedPower; MAX_LAYERS],
    layer_count: usize,
}
impl LayerAbsorption {
    /// One row per original layer, including exact-zero lossless/identity rows.
    pub fn layers(&self) -> &[LayerAbsorbedPower] {
        &self.layers[..self.layer_count]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum LayerAbsorptionError {
    Material(MaterialError),
    Invalid(&'static str),
    Cancelled,
}
impl From<MaterialError> for LayerAbsorptionError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for LayerAbsorptionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "layer absorption refused: {self:?}")
    }
}
impl std::error::Error for LayerAbsorptionError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LayerAbsorptionError> {
    if cancelled() {
        Err(LayerAbsorptionError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit(
    layers: usize,
    cosine: f64,
    incident_power_w: f64,
    budget: LayerAbsorptionBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LayerAbsorptionError> {
    poll(cancelled)?;
    if budget.max_layers > MAX_LAYERS || layers > budget.max_layers {
        return Err(LayerAbsorptionError::Invalid(
            "layer absorption work budget exceeded",
        ));
    }
    if !cosine.is_finite()
        || cosine <= 0.0
        || cosine > 1.0
        || !incident_power_w.is_finite()
        || incident_power_w < 0.0
    {
        return Err(LayerAbsorptionError::Invalid(
            "layer absorption needs cosine in (0,1] and finite nonnegative incident watts",
        ));
    }
    Ok(())
}

impl RigidBackedLayers<'_> {
    /// Integrate each layer's original dynamic-density and bulk-modulus loss.
    /// Uses one original reflection pass plus a linear analytic integration
    /// pass, fixed stack storage and no additional material/source queries.
    /// Cancellation is polled per layer in BOTH passes and before publication.
    pub fn layer_absorption(
        &self,
        cosine: f64,
        incident_power_w: f64,
        budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerAbsorption, LayerAbsorptionError> {
        admit(
            self.layer_count(),
            cosine,
            incident_power_w,
            budget,
            &mut cancelled,
        )?;
        let mut transfers = [LayerTransfer::default(); MAX_LAYERS];
        let response = self.reflection_with_layers::<LayerAbsorptionError>(cosine, |i, step| {
            poll(&mut cancelled)?;
            transfers[i] = step;
            Ok(())
        })?;
        partition(
            self,
            cosine,
            incident_power_w,
            [response.pressure_reflection.norm_sq(), 0.0],
            &transfers,
            &mut cancelled,
        )
    }
}

impl TransmittingLayers<'_> {
    /// The same local positive loss integrals with the actual outgoing-fluid
    /// termination. The original R/T amplitudes and receiving power weights
    /// are unchanged; this report resolves loss independently of their residual.
    pub fn layer_absorption(
        &self,
        cosine: f64,
        incident_power_w: f64,
        budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerAbsorption, LayerAbsorptionError> {
        let prepared = self.prepared_layers();
        admit(
            prepared.layer_count(),
            cosine,
            incident_power_w,
            budget,
            &mut cancelled,
        )?;
        let mut transfers = [LayerTransfer::default(); MAX_LAYERS];
        let response = self.scattering_with_layers::<LayerAbsorptionError>(cosine, |i, step| {
            poll(&mut cancelled)?;
            transfers[i] = step;
            Ok(())
        })?;
        partition(
            prepared,
            cosine,
            incident_power_w,
            [response.reflected_fraction, response.transmitted_fraction],
            &transfers,
            &mut cancelled,
        )
    }
}

fn partition(
    stack: &RigidBackedLayers<'_>,
    cosine: f64,
    incident_power_w: f64,
    outgoing: [f64; 2],
    transfers: &[LayerTransfer; MAX_LAYERS],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<LayerAbsorption, LayerAbsorptionError> {
    let mut layers = [LayerAbsorbedPower::default(); MAX_LAYERS];
    let mut total = Dd::ZERO;
    // The incident pressure is one: p_front + Z_ref*v_front = 2. Only
    // magnitude is needed; a common phase cancels in every work integral.
    let front = transfers[0];
    let mut amplitude = LogAmplitude::new(if stack.layer_count() == 0 {
        0.0
    } else {
        core::f64::consts::LN_2 - det::ln((front.front_pressure + front.front_velocity).abs())
    });
    for (i, layer) in stack.layers.iter().enumerate() {
        poll(cancelled)?;
        let step = transfers[i];
        let d = layer.source.thickness_m();
        if d == 0.0 {
            continue;
        }
        let front_log = amplitude.log;
        let kn = step.normal_wavenumber;
        amplitude.add(if kn.im == 0.0 {
            -log_scaled(16.0, step.normalization)
        } else {
            kn.im * d - log_scaled(8.0, step.normalization)
        });
        if let Some(response) = layer.porous {
            let [viscous, thermal] = layer_loss(
                stack,
                cosine,
                d,
                &response,
                step,
                [front_log, amplitude.log],
            )?;
            let absorbed = viscous + thermal;
            let power = incident_power_w * absorbed;
            if !power.is_finite() || !absorbed.is_finite() {
                return Err(LayerAbsorptionError::Invalid(
                    "layer absorbed power is unrepresentable",
                ));
            }
            layers[i] = LayerAbsorbedPower {
                viscous_absorbed_fraction: viscous,
                thermal_absorbed_fraction: thermal,
                absorbed_power_w: power,
            };
            total = total + Dd::from_f64(viscous) + Dd::from_f64(thermal);
        }
    }
    let absorbed = total.to_f64();
    let [reflected, transmitted] = outgoing;
    let fraction_residual =
        (Dd::from_f64(1.0) - Dd::from_f64(reflected) - Dd::from_f64(transmitted) - total).to_f64();
    if !absorbed.is_finite()
        || !fraction_residual.is_finite()
        || fraction_residual.abs() > 1024.0 * f64::EPSILON * (stack.layer_count() as f64 + 1.0)
    {
        return Err(LayerAbsorptionError::Invalid(
            "local absorption disagrees with boundary power",
        ));
    }
    let powers = [reflected, transmitted, absorbed].map(|fraction| incident_power_w * fraction);
    let residual = (Dd::from_f64(incident_power_w)
        - Dd::from_f64(powers[0])
        - Dd::from_f64(powers[1])
        - Dd::from_f64(powers[2]))
    .to_f64();
    if !powers.into_iter().all(f64::is_finite) || !residual.is_finite() {
        return Err(LayerAbsorptionError::Invalid(
            "layer power account is unrepresentable",
        ));
    }
    poll(cancelled)?;
    Ok(LayerAbsorption {
        incident_power_w,
        reflected_power_w: powers[0],
        transmitted_power_w: powers[1],
        absorbed_power_w: powers[2],
        reflected_fraction: reflected,
        transmitted_fraction: transmitted,
        absorbed_fraction: absorbed,
        balance_residual_w: residual,
        layers,
        layer_count: stack.layer_count(),
    })
}

fn layer_loss(
    stack: &RigidBackedLayers<'_>,
    cosine: f64,
    thickness: f64,
    response: &super::PorousFluidResponse,
    step: LayerTransfer,
    amplitude_logs: [f64; 2],
) -> Result<[f64; 2], LayerAbsorptionError> {
    let kn = step.normal_wavenumber;
    let (zc, k) = response.characteristic();
    let ratio = divided(zc * (k / kn), stack.z0).scale(cosine);
    let [log_p, log_q] = field_integral_logs(step, thickness, ratio, amplitude_logs)?;
    let density = response.dynamic_density_kg_m3();
    let inverse_density = density.recip();
    let compliance = response.dynamic_bulk_modulus_pa().recip();
    if !finite(inverse_density)
        || !finite(compliance)
        || inverse_density.im < 0.0
        || compliance.im > 0.0
    {
        return Err(LayerAbsorptionError::Invalid(
            "local passive constitutive loss is unrepresentable",
        ));
    }
    let log_omega = det::ln(stack.omega);
    let log_reference_impedance = det::ln(stack.z0) - det::ln(cosine);
    let kt = stack.k0 * det::sqrt((1.0 - cosine) * (1.0 + cosine));
    let normal = log_omega + det::ln(-density.im) - log_reference_impedance + log_q;
    let tangential = 2.0 * det::ln(kt) - log_omega
        + det::ln(inverse_density.im)
        + log_reference_impedance
        + log_p;
    let thermal = log_omega + det::ln(-compliance.im) + log_reference_impedance + log_p;
    let result = [det::exp(log_add(normal, tangential)), det::exp(thermal)];
    if result.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(LayerAbsorptionError::Invalid(
            "integrated constitutive loss is unrepresentable",
        ));
    }
    Ok(result)
}

// Integrals of |p|² and |Z_ref*v_n|² for unit incident pressure, retained
// logarithmically until the actual constitutive loss supplies its coefficient.
fn field_integral_logs(
    step: LayerTransfer,
    thickness: f64,
    ratio: C64,
    amplitude_logs: [f64; 2],
) -> Result<[f64; 2], LayerAbsorptionError> {
    let common_log = amplitude_logs[0].max(amplitude_logs[1]);
    if common_log == f64::NEG_INFINITY {
        return Ok([f64::NEG_INFINITY; 2]);
    }
    let inverse = ratio.recip();
    if !finite(ratio) || !finite(inverse) || !common_log.is_finite() {
        return Err(LayerAbsorptionError::Invalid(
            "local wave amplitudes are unrepresentable",
        ));
    }
    // Scale pressure and velocity separately: using the inverse-impedance
    // scale for pressure would erase a representable O(cos(theta)) grazing
    // pressure before its O(1/cos(theta)) power coefficient can participate.
    let pressure_reference = magnitude_scale(ratio).max(1.0);
    let velocity_reference = magnitude_scale(inverse).max(1.0);
    let r = divided(ratio, pressure_reference);
    let y = divided(inverse, velocity_reference);
    let p_back = divided(step.back_pressure, pressure_reference);
    let q_back = divided(step.back_velocity, velocity_reference);
    let kn = step.normal_wavenumber;
    let x = -kn.im * thickness;
    let phase = kn.re * thickness;
    if x < 0.1 && phase.abs() < 0.1 {
        // At a thin rigid-backed layer, v_n is O(d), while the two traveling
        // wave velocities are O(1). Work from their difference would erase the
        // O(d³) viscous integral. Center fields and parity remove that subtraction:
        // integral |p|² = I_cos |p_mid|² + I_sin |r*q_mid|², and similarly for q.
        let half = kn.scale(0.5 * thickness);
        let (sin_re, cos_re) = det::sin_cos(half.re);
        let sinh_im = 0.5 * (det::expm1(half.im) - det::expm1(-half.im));
        let cosh_im = 0.5 * (det::exp(half.im) + det::exp(-half.im));
        let cos = C64::new(cos_re * cosh_im, -sin_re * sinh_im);
        let i_sin = C64::new(-cos_re * sinh_im, sin_re * cosh_im);
        let r_q_back = r * step.back_velocity;
        let y_p_back = y * step.back_pressure;
        let p = cos * p_back + i_sin * r_q_back;
        let q = cos * q_back + i_sin * y_p_back;
        let r_q = cos * r_q_back + i_sin * p_back;
        let y_p = cos * y_p_back + i_sin * q_back;
        let x2 = x * x;
        let p2 = phase * phase;
        let sinhc_minus_one = x2
            * (1.0 / 6.0
                + x2 * (1.0 / 120.0
                    + x2 * (1.0 / 5040.0 + x2 * (1.0 / 362880.0 + x2 / 39916800.0))));
        let one_minus_sinc = p2
            * (1.0 / 6.0
                - p2 * (1.0 / 120.0
                    - p2 * (1.0 / 5040.0 - p2 * (1.0 / 362880.0 - p2 / 39916800.0))));
        let i_cos = thickness * (1.0 + 0.5 * sinhc_minus_one - 0.5 * one_minus_sinc);
        let i_sin = (0.5 * thickness) * (sinhc_minus_one + one_minus_sinc);
        return Ok([
            2.0 * (amplitude_logs[1] + det::ln(pressure_reference))
                + log_add(
                    det::ln(i_cos) + 2.0 * det::ln(p.abs()),
                    det::ln(i_sin) + 2.0 * det::ln(r_q.abs()),
                ),
            2.0 * (amplitude_logs[1] + det::ln(velocity_reference))
                + log_add(
                    det::ln(i_cos) + 2.0 * det::ln(q.abs()),
                    det::ln(i_sin) + 2.0 * det::ln(y_p.abs()),
                ),
        ]);
    }
    let p_front = divided(step.front_pressure, pressure_reference);
    let q_front = divided(step.front_velocity, velocity_reference);
    let front_scale = det::exp(amplitude_logs[0] - common_log) * 0.5;
    let back_scale = det::exp(amplitude_logs[1] - common_log) * 0.5;
    let a = (p_front + r * step.front_velocity).scale(front_scale);
    let b = (p_back - r * step.back_velocity).scale(back_scale);
    let c = (y * step.front_pressure + q_front).scale(front_scale);
    let e = (y * step.back_pressure - q_back).scale(back_scale);
    let gram = wave_gram(kn, thickness)?;
    // The backward propagation's bounded factor is exp(-i*k_n*d)/8.
    // Align the right-face modes with the left-face phase only when their
    // interference survives; an opaque layer needs no unrepresentable phase.
    let rotation = if gram[0] == gram[1] || back_scale == 0.0 {
        C64::ONE
    } else {
        let (sine, cosine) = det::sin_cos(-phase);
        C64::new(cosine, sine)
    };
    Ok([
        2.0 * (common_log + det::ln(pressure_reference)) + mode_integral_log(a, b * rotation, gram),
        2.0 * (common_log + det::ln(velocity_reference))
            + mode_integral_log(c, -(e * rotation), gram),
    ])
}

// For p=A exp(-ikz)+B exp(-ik(d-z)), integral |p|² is
// g_plus |A+B|² + g_minus |A-B|², g_±=(J±L)/2,
// J=(1-exp(-2ad))/(2a), L=exp(-ad) sin(bd)/b, k=b-ia.
// Both Gram eigenvalues are nonnegative. Optically thin layers instead use
// the midpoint identity above, so near-equal J and L are never subtracted.
fn wave_gram(kn: C64, d: f64) -> Result<[f64; 2], LayerAbsorptionError> {
    let a = -kn.im;
    let x = a * d;
    let phase = kn.re * d;
    let attenuation = det::exp(-x);
    let j = if x < 1.0e-4 {
        d * if x == 0.0 {
            1.0
        } else {
            -det::expm1(-2.0 * x) / (2.0 * x)
        }
    } else {
        (0.5 * -det::expm1(-2.0 * x)) / a
    };
    let l = if attenuation == 0.0 {
        0.0
    } else if phase.abs() < 1.0e-4 {
        (d * attenuation)
            * if phase == 0.0 {
                1.0
            } else {
                det::sin_cos(phase).0 / phase
            }
    } else {
        attenuation * (det::sin_cos(phase).0 / kn.re)
    };
    let result = [0.5 * j + 0.5 * l, 0.5 * j - 0.5 * l];
    if result.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(LayerAbsorptionError::Invalid(
            "local wave integral is unrepresentable",
        ));
    }
    Ok(result)
}

fn log_add(a: f64, b: f64) -> f64 {
    let larger = a.max(b);
    if larger == f64::NEG_INFINITY {
        larger
    } else {
        larger + det::ln(1.0 + det::exp(a.min(b) - larger))
    }
}

fn mode_integral_log(a: C64, b: C64, gram: [f64; 2]) -> f64 {
    log_add(
        det::ln(gram[0]) + 2.0 * det::ln((a + b).abs()),
        det::ln(gram[1]) + 2.0 * det::ln((a - b).abs()),
    )
}
