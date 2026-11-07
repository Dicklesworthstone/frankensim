//! Original angle-dependent magnetic storage and its symmetric component work.
use super::saturation::MagneticLaw;
use super::{SynchronousMotorParameters, SynchronousOperatingPoint, invalid};
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::{dd::Dd, det};

#[derive(Clone, Copy)]
pub(super) struct SynchronousStorage {
    parameters: SynchronousMotorParameters,
    inverse: [f64; 4], // Ld, Lq, Jr, Js
    pole_pairs: f64,
    pm_reluctance_flux: f64,
}

fn complete(values: &[f64]) -> bool {
    values.len() == 5 && values.iter().all(|v| v.is_finite())
}

pub(super) fn to_dq(sine: f64, cosine: f64, value: [f64; 2]) -> [f64; 2] {
    [
        cosine.mul_add(value[0], sine * value[1]),
        (-sine).mul_add(value[0], cosine * value[1]),
    ]
}

pub(super) fn to_ab(sine: f64, cosine: f64, value: [f64; 2]) -> [f64; 2] {
    [
        cosine.mul_add(value[0], -sine * value[1]),
        sine.mul_add(value[0], cosine * value[1]),
    ]
}

pub(super) fn operating_point(
    p: SynchronousMotorParameters,
    law: &MagneticLaw,
    state: &[f64; 5],
) -> Result<SynchronousOperatingPoint, PhsError> {
    let (s, c) = det::sin_cos(f64::from(p.pole_pairs) * state[2]);
    let flux_dq_wb = to_dq(s, c, [state[0], state[1]]);
    let current_dq_a = law.current_dq(p, flux_dq_wb)?;
    let differential_inductance_dq_h = law.differential_inductance_dq(p, flux_dq_wb)?;
    let current_alpha_beta_a = to_ab(s, c, current_dq_a);
    let phase_a = det::sqrt(2.0 / 3.0) * current_alpha_beta_a[0];
    let quadrature = det::sqrt(0.5) * current_alpha_beta_a[1];
    let phase_current_a = [
        phase_a,
        -0.5 * phase_a + quadrature,
        -0.5 * phase_a - quadrature,
    ];
    let shaft_speed_rad_s = [
        state[3] / p.rotor_inertia_kg_m2,
        state[4] / p.stator_inertia_kg_m2,
    ];
    let electrical_speed = f64::from(p.pole_pairs) * (shaft_speed_rad_s[0] - shaft_speed_rad_s[1]);
    let electromagnetic_torque_nm = match law {
        MagneticLaw::Linear => {
            (f64::from(p.pole_pairs) * current_dq_a[1])
                * (p.d_inductance_h - p.q_inductance_h)
                    .mul_add(current_dq_a[0], p.pm_flux_linkage_wb)
        }
        MagneticLaw::Curves { .. } => {
            f64::from(p.pole_pairs)
                * flux_dq_wb[0].mul_add(current_dq_a[1], -flux_dq_wb[1] * current_dq_a[0])
        }
    };
    let speed_voltage_dq_v = [
        -electrical_speed * flux_dq_wb[1],
        electrical_speed * flux_dq_wb[0],
    ];
    if !complete(state)
        || current_alpha_beta_a
            .iter()
            .chain(&current_dq_a)
            .chain(&phase_current_a)
            .chain(&flux_dq_wb)
            .chain(&differential_inductance_dq_h)
            .chain(&shaft_speed_rad_s)
            .chain(&speed_voltage_dq_v)
            .chain(&[electromagnetic_torque_nm])
            .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "synchronous operating current, torque or speed voltage is unrepresentable",
        ));
    }
    Ok(SynchronousOperatingPoint {
        current_alpha_beta_a,
        current_dq_a,
        phase_current_a,
        flux_dq_wb,
        differential_inductance_dq_h,
        shaft_speed_rad_s,
        electromagnetic_torque_nm,
        speed_voltage_dq_v,
    })
}

struct Pair {
    sine: f64,
    cosine: f64,
    half_angle: f64,
    sine_half: f64,
    cosine_half: f64,
    sinc_half: f64,
    mean: [f64; 2],
    delta: [f64; 2],
    mean_reluctivity: [f64; 2],
}

impl SynchronousStorage {
    pub(super) fn new(parameters: SynchronousMotorParameters) -> Result<Self, PhsError> {
        let values = [
            parameters.d_inductance_h,
            parameters.q_inductance_h,
            parameters.rotor_inertia_kg_m2,
            parameters.stator_inertia_kg_m2,
        ];
        let inverse = values.map(|value| 1.0 / value);
        if values
            .iter()
            .chain(&inverse)
            .any(|v| !v.is_finite() || *v <= 0.0)
            || parameters.pole_pairs == 0
            || !parameters.pm_flux_linkage_wb.is_finite()
        {
            return Err(invalid(
                "synchronous machine needs finite positive representable Ld/Lq/inertias, pole pairs and finite PM flux",
            ));
        }
        let pm_reluctance_flux = parameters.pm_flux_linkage_wb * inverse[0];
        if !pm_reluctance_flux.is_finite() {
            return Err(invalid(
                "synchronous PM constitutive coefficient is unrepresentable",
            ));
        }
        Ok(Self {
            parameters,
            inverse,
            pole_pairs: f64::from(parameters.pole_pairs),
            pm_reluctance_flux,
        })
    }

    fn pair(&self, a: &[f64], b: &[f64]) -> Option<Pair> {
        if !complete(a) || !complete(b) {
            return None;
        }
        let middle_angle = self.pole_pairs * f64::midpoint(a[2], b[2]);
        let half_angle = (0.5 * self.pole_pairs) * (b[2] - a[2]);
        let (sine, cosine) = det::sin_cos(middle_angle);
        let (sine_half, cosine_half) = det::sin_cos(half_angle);
        let mean = to_dq(
            sine,
            cosine,
            [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])],
        );
        let delta = to_dq(sine, cosine, [b[0] - a[0], b[1] - a[1]]);
        let sinc_half = if half_angle == 0.0 {
            1.0
        } else {
            sine_half / half_angle
        };
        let cc = cosine_half * cosine_half;
        let ss = sine_half * sine_half;
        // Half the two endpoint reluctivity tensors, in their common middle
        // frame. Positive weighted sums retain either small principal value.
        let mean_reluctivity = [
            self.inverse[0].mul_add(cc, self.inverse[1] * ss),
            self.inverse[0].mul_add(ss, self.inverse[1] * cc),
        ];
        if [
            middle_angle,
            half_angle,
            sine,
            cosine,
            sine_half,
            cosine_half,
            sinc_half,
        ]
        .iter()
        .chain(&mean)
        .chain(&delta)
        .chain(&mean_reluctivity)
        .any(|v| !v.is_finite())
        {
            return None;
        }
        Some(Pair {
            sine,
            cosine,
            half_angle,
            sine_half,
            cosine_half,
            sinc_half,
            mean,
            delta,
            mean_reluctivity,
        })
    }

    fn magnetic_effort(&self, pair: &Pair) -> ([f64; 2], f64) {
        let dq = [
            pair.mean_reluctivity[0]
                .mul_add(pair.mean[0], -self.pm_reluctance_flux * pair.cosine_half),
            pair.mean_reluctivity[1] * pair.mean[1],
        ];
        let product_mean = pair.mean[0].mul_add(pair.mean[1], 0.25 * pair.delta[0] * pair.delta[1]);
        let angle = self.pole_pairs
            * pair.sinc_half
            * ((self.inverse[0] - self.inverse[1]) * pair.cosine_half * product_mean
                - self.pm_reluctance_flux * pair.mean[1]);
        (dq, angle)
    }

    // Independent constitutive energy difference, using differences of squares
    // of the ORIGINAL d/q flux deviations and momenta. No accepted torque,
    // port work or dissipation is used to manufacture the energy ledger.
    fn energy_change(&self, a: &[f64], b: &[f64], pair: &Pair) -> Dd {
        let mean_d = pair
            .cosine_half
            .mul_add(pair.mean[0], 0.5 * pair.sine_half * pair.delta[1])
            - self.parameters.pm_flux_linkage_wb;
        let mean_q = pair
            .cosine_half
            .mul_add(pair.mean[1], -0.5 * pair.sine_half * pair.delta[0]);
        let delta_d = pair
            .cosine_half
            .mul_add(pair.delta[0], 2.0 * pair.sine_half * pair.mean[1]);
        let delta_q = pair
            .cosine_half
            .mul_add(pair.delta[1], -2.0 * pair.sine_half * pair.mean[0]);
        let mut change =
            Dd::from_f64(self.inverse[0]) * Dd::from_f64(mean_d) * Dd::from_f64(delta_d)
                + Dd::from_f64(self.inverse[1]) * Dd::from_f64(mean_q) * Dd::from_f64(delta_q);
        for i in 3..5 {
            let mean = Dd::from_f64(0.5 * a[i]) + Dd::from_f64(0.5 * b[i]);
            change = change + Dd::from_f64(self.inverse[i - 1]) * mean * Dd::from_f64(b[i] - a[i]);
        }
        change
    }
}

impl Storage for SynchronousStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(*self))
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }

    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !complete(x) {
            return f64::NAN;
        }
        let (s, c) = det::sin_cos(self.pole_pairs * x[2]);
        let flux = to_dq(s, c, [x[0], x[1]]);
        let d = flux[0] - self.parameters.pm_flux_linkage_wb;
        0.5 * (d * (d * self.inverse[0])
            + flux[1] * (flux[1] * self.inverse[1])
            + x[3] * (x[3] * self.inverse[2])
            + x[4] * (x[4] * self.inverse[3]))
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        // Coincidence is the analytic gradient of the same original storage.
        self.evaluate(x, x, x, out);
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let Some(pair) = self.pair(a, b) else {
            return Some(f64::NAN);
        };
        let mut gradient = [0.0; 5];
        self.gradient(mid, &mut gradient);
        let mut remainder = self.energy_change(a, b, &pair);
        for i in 0..5 {
            remainder = remainder - Dd::from_f64(gradient[i]) * Dd::from_f64(b[i] - a[i]);
        }
        Some(remainder.to_f64())
    }
}

impl ConstitutiveDiscreteGradient for SynchronousStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], _mid: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != 5 {
            return;
        }
        let Some(pair) = self.pair(a, b) else {
            return;
        };
        // Symmetric coordinate split: average the exact flux gradient at both
        // angles, and the exact angle divided difference at both flux endpoints.
        // sin(h)/h supplies its continuous zero-angle limit without subtracting
        // nearly equal magnetic energies. Kinetic work remains independent.
        let (dq, angle) = self.magnetic_effort(&pair);
        out[..2].copy_from_slice(&to_ab(pair.sine, pair.cosine, dq));
        out[2] = angle;
        out[3] = f64::midpoint(a[3], b[3]) * self.inverse[2];
        out[4] = f64::midpoint(a[4], b[4]) * self.inverse[3];
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !complete(direction) || out.len() != 5 {
            return false;
        }
        let Some(pair) = self.pair(a, b) else {
            return false;
        };
        let dangle = (0.5 * self.pole_pairs) * direction[2];
        // Both the middle angle and half endpoint-angle difference have this
        // SAME end-state derivative. No midpoint-Hessian substitution occurs.
        let flux_direction = to_dq(pair.sine, pair.cosine, [direction[0], direction[1]]);
        let dmean = [
            0.5 * flux_direction[0] + dangle * pair.mean[1],
            0.5 * flux_direction[1] - dangle * pair.mean[0],
        ];
        let ddelta = [
            flux_direction[0] + dangle * pair.delta[1],
            flux_direction[1] - dangle * pair.delta[0],
        ];
        let dreluctivity =
            2.0 * (self.inverse[1] - self.inverse[0]) * pair.cosine_half * pair.sine_half * dangle;
        let (dq, _) = self.magnetic_effort(&pair);
        let ddq = [
            pair.mean_reluctivity[0] * dmean[0]
                + dreluctivity * pair.mean[0]
                + self.pm_reluctance_flux * pair.sine_half * dangle
                - dq[1] * dangle,
            pair.mean_reluctivity[1] * dmean[1] - dreluctivity * pair.mean[1] + dq[0] * dangle,
        ];
        out[..2].copy_from_slice(&to_ab(pair.sine, pair.cosine, ddq));
        let product = pair.mean[0] * pair.mean[1] + 0.25 * pair.delta[0] * pair.delta[1];
        let dproduct = dmean[0] * pair.mean[1]
            + pair.mean[0] * dmean[1]
            + 0.25 * (ddelta[0] * pair.delta[1] + pair.delta[0] * ddelta[1]);
        let inverse_difference = self.inverse[0] - self.inverse[1];
        let angular = inverse_difference * pair.cosine_half * product
            - self.pm_reluctance_flux * pair.mean[1];
        let dangular = inverse_difference
            * (pair.cosine_half * dproduct - pair.sine_half * dangle * product)
            - self.pm_reluctance_flux * dmean[1];
        let sinc_derivative = if pair.half_angle.abs() < 0.01 {
            // Taylor evaluation of the SAME sinc derivative near coincidence;
            // the omitted term here is below floating-point resolution.
            let square = pair.half_angle * pair.half_angle;
            pair.half_angle
                * (-1.0 / 3.0 + square * (1.0 / 30.0 + square * (-1.0 / 840.0 + square / 45360.0)))
        } else {
            (pair.cosine_half - pair.sinc_half) / pair.half_angle
        };
        out[2] = self.pole_pairs * (sinc_derivative * dangle * angular + pair.sinc_half * dangular);
        out[3] = 0.5 * self.inverse[2] * direction[3];
        out[4] = 0.5 * self.inverse[3] * direction[4];
        out.iter().all(|v| v.is_finite())
    }
}

impl ContinuousStorageTangent for SynchronousStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        // Symmetry and consistency give D_2 DG(x,x) = H''(x)/2. Reusing the
        // complete analytic action keeps the continuous and discrete laws tied.
        if !self.end_jacobian_action(state, state, direction, out) {
            return false;
        }
        for value in out.iter_mut() {
            *value *= 2.0;
        }
        out.iter().all(|v| v.is_finite())
    }
}
