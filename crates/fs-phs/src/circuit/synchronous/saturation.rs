//! Independent original d/q flux laws, exact component work and its tangent.
//!
//! The d law acts on excess flux lambda_d - psi_pm; the q law acts on lambda_q.
//! These are caller-authored reversible laws, with no inferred cross-saturation.
use super::{
    SynchronousMotorParameters, invalid,
    storage::{to_ab, to_dq},
};
use crate::circuit::inductance::InductorFluxCurve;
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::{dd::Dd, det};

#[derive(Clone)]
pub(super) enum MagneticLaw {
    Linear,
    Curves {
        d: InductorFluxCurve,
        q: InductorFluxCurve,
    },
}

impl MagneticLaw {
    pub(super) fn flux_dq(
        &self,
        p: SynchronousMotorParameters,
        current: [f64; 2],
    ) -> Result<[f64; 2], PhsError> {
        let flux = match self {
            Self::Linear => [
                p.d_inductance_h.mul_add(current[0], p.pm_flux_linkage_wb),
                p.q_inductance_h * current[1],
            ],
            Self::Curves { d, q } => [
                d.flux_wb(current[0]).map_err(|_| {
                    invalid("synchronous d current is outside its original flux law")
                })? + p.pm_flux_linkage_wb,
                q.flux_wb(current[1]).map_err(|_| {
                    invalid("synchronous q current is outside its original flux law")
                })?,
            ],
        };
        if current.iter().chain(&flux).any(|v| !v.is_finite()) {
            return Err(invalid(
                "synchronous current-to-flux relation is unrepresentable",
            ));
        }
        Ok(flux)
    }

    pub(super) fn current_dq(
        &self,
        p: SynchronousMotorParameters,
        flux: [f64; 2],
    ) -> Result<[f64; 2], PhsError> {
        match self {
            Self::Linear => Ok([
                (flux[0] - p.pm_flux_linkage_wb) / p.d_inductance_h,
                flux[1] / p.q_inductance_h,
            ]),
            Self::Curves { d, q } => Ok([
                d.current_a(flux[0] - p.pm_flux_linkage_wb)
                    .map_err(|_| invalid("synchronous d flux is outside its original law"))?,
                q.current_a(flux[1])
                    .map_err(|_| invalid("synchronous q flux is outside its original law"))?,
            ]),
        }
    }

    pub(super) fn differential_inductance_dq(
        &self,
        p: SynchronousMotorParameters,
        flux: [f64; 2],
    ) -> Result<[f64; 2], PhsError> {
        match self {
            Self::Linear => Ok([p.d_inductance_h, p.q_inductance_h]),
            Self::Curves { d, q } => Ok([
                d.differential_inductance_h(flux[0] - p.pm_flux_linkage_wb)
                    .map_err(|_| invalid("synchronous d flux is outside its original law"))?,
                q.differential_inductance_h(flux[1])
                    .map_err(|_| invalid("synchronous q flux is outside its original law"))?,
            ]),
        }
    }
}

#[derive(Clone)]
pub(super) struct SaturatedStorage {
    parameters: SynchronousMotorParameters,
    curves: [InductorFluxCurve; 2],
    inverse_inertia: [f64; 2],
    pole_pairs: f64,
}

fn complete(values: &[f64]) -> bool {
    values.len() == 5 && values.iter().all(|v| v.is_finite())
}

// Each corner is original physical d/q flux, indexed [fixed flux][fixed angle].
// A shared middle frame keeps the small angular increment explicit.
struct Pair {
    sine: f64,
    cosine: f64,
    half_angle: f64,
    sine_half: f64,
    cosine_half: f64,
    sinc_half: f64,
    middle: [[f64; 2]; 2],
    corner: [[[f64; 2]; 2]; 2],
}

impl Pair {
    fn to_angle(&self, angle: usize, stationary: [f64; 2]) -> [f64; 2] {
        let middle = to_dq(self.sine, self.cosine, stationary);
        if angle == 0 {
            to_ab(self.sine_half, self.cosine_half, middle)
        } else {
            to_dq(self.sine_half, self.cosine_half, middle)
        }
    }

    fn from_angle(&self, angle: usize, dq: [f64; 2]) -> [f64; 2] {
        let middle = if angle == 0 {
            to_dq(self.sine_half, self.cosine_half, dq)
        } else {
            to_ab(self.sine_half, self.cosine_half, dq)
        };
        to_ab(self.sine, self.cosine, middle)
    }

    fn sinc_derivative(&self) -> f64 {
        if self.half_angle.abs() < 0.01 {
            let square = self.half_angle * self.half_angle;
            self.half_angle
                * (-1.0 / 3.0 + square * (1.0 / 30.0 + square * (-1.0 / 840.0 + square / 45360.0)))
        } else {
            (self.cosine_half - self.sinc_half) / self.half_angle
        }
    }
}

impl SaturatedStorage {
    pub(super) fn new(
        parameters: SynchronousMotorParameters,
        d: InductorFluxCurve,
        q: InductorFluxCurve,
    ) -> Result<Self, PhsError> {
        let inertia = [
            parameters.rotor_inertia_kg_m2,
            parameters.stator_inertia_kg_m2,
        ];
        let inverse_inertia = inertia.map(|v| 1.0 / v);
        if inertia
            .iter()
            .chain(&inverse_inertia)
            .any(|v| !v.is_finite() || *v <= 0.0)
            || parameters.pole_pairs == 0
            || !parameters.pm_flux_linkage_wb.is_finite()
        {
            return Err(invalid(
                "saturated synchronous machine needs positive representable inertias, pole pairs and finite PM flux",
            ));
        }
        Ok(Self {
            parameters,
            curves: [d, q],
            inverse_inertia,
            pole_pairs: f64::from(parameters.pole_pairs),
        })
    }

    fn physical_flux(&self, state: &[f64]) -> Option<[f64; 2]> {
        if !complete(state) {
            return None;
        }
        let (s, c) = det::sin_cos(self.pole_pairs * state[2]);
        let flux = to_dq(s, c, [state[0], state[1]]);
        flux.iter().all(|v| v.is_finite()).then_some(flux)
    }

    fn coordinate(&self, flux: [f64; 2]) -> [f64; 2] {
        [flux[0] - self.parameters.pm_flux_linkage_wb, flux[1]]
    }

    fn pair(&self, a: &[f64], b: &[f64]) -> Option<Pair> {
        if !complete(a) || !complete(b) {
            return None;
        }
        let middle_angle = self.pole_pairs * f64::midpoint(a[2], b[2]);
        let half_angle = (0.5 * self.pole_pairs) * (b[2] - a[2]);
        let (sine, cosine) = det::sin_cos(middle_angle);
        let (sine_half, cosine_half) = det::sin_cos(half_angle);
        let sinc_half = if half_angle == 0.0 {
            1.0
        } else {
            sine_half / half_angle
        };
        let middle = [
            to_dq(sine, cosine, [a[0], a[1]]),
            to_dq(sine, cosine, [b[0], b[1]]),
        ];
        let corner = middle.map(|flux| {
            [
                to_ab(sine_half, cosine_half, flux),
                to_dq(sine_half, cosine_half, flux),
            ]
        });
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
        .chain(middle.iter().flatten())
        .chain(corner.iter().flatten().flatten())
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
            middle,
            corner,
        })
    }

    fn chord(&self, a: [f64; 2], b: [f64; 2]) -> Option<[f64; 2]> {
        let (a, b) = (self.coordinate(a), self.coordinate(b));
        Some([
            self.curves[0].discrete_current_a(a[0], b[0]).ok()?,
            self.curves[1].discrete_current_a(a[1], b[1]).ok()?,
        ])
    }

    fn chord_action(
        &self,
        a: [f64; 2],
        b: [f64; 2],
        da: [f64; 2],
        db: [f64; 2],
    ) -> Option<[f64; 2]> {
        let (a, b) = (self.coordinate(a), self.coordinate(b));
        let mut action = [0.0; 2];
        for k in 0..2 {
            let start = self.curves[k]
                .discrete_current_end_derivative(b[k], a[k])
                .ok()?;
            let end = self.curves[k]
                .discrete_current_end_derivative(a[k], b[k])
                .ok()?;
            action[k] = start.mul_add(da[k], end * db[k]);
        }
        Some(action)
    }

    fn effort(&self, a: &[f64], b: &[f64]) -> Option<[f64; 5]> {
        let pair = self.pair(a, b)?;
        let mut electrical = [[0.0; 2]; 2];
        let mut angular = [0.0; 2];
        for k in 0..2 {
            electrical[k] = pair.from_angle(k, self.chord(pair.corner[0][k], pair.corner[1][k])?);
            let chord = self.chord(pair.corner[k][0], pair.corner[k][1])?;
            angular[k] = pair.middle[k][1].mul_add(chord[0], -pair.middle[k][0] * chord[1]);
        }
        let result = [
            f64::midpoint(electrical[0][0], electrical[1][0]),
            f64::midpoint(electrical[0][1], electrical[1][1]),
            self.pole_pairs * pair.sinc_half * f64::midpoint(angular[0], angular[1]),
            f64::midpoint(a[3], b[3]) * self.inverse_inertia[0],
            f64::midpoint(a[4], b[4]) * self.inverse_inertia[1],
        ];
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    fn action(&self, a: &[f64], b: &[f64], direction: &[f64]) -> Option<[f64; 5]> {
        if !complete(direction) {
            return None;
        }
        let pair = self.pair(a, b)?;
        let eta = self.pole_pairs * direction[2];
        let gamma = 0.5 * eta;
        let flux_direction = [direction[0], direction[1]];
        let mut corner_action = [[[0.0; 2]; 2]; 2];
        corner_action[1][0] = pair.to_angle(0, flux_direction);
        let at_end = pair.to_angle(1, flux_direction);
        for k in 0..2 {
            corner_action[k][1] = [
                eta.mul_add(pair.corner[k][1][1], if k == 1 { at_end[0] } else { 0.0 }),
                (-eta).mul_add(pair.corner[k][1][0], if k == 1 { at_end[1] } else { 0.0 }),
            ];
        }
        let middle_direction = to_dq(pair.sine, pair.cosine, flux_direction);
        let mut electrical_action = [[0.0; 2]; 2];
        let mut angular = [0.0; 2];
        let mut angular_action = [0.0; 2];
        for k in 0..2 {
            let chord = self.chord(pair.corner[0][k], pair.corner[1][k])?;
            let mut dchord = self.chord_action(
                pair.corner[0][k],
                pair.corner[1][k],
                corner_action[0][k],
                corner_action[1][k],
            )?;
            if k == 1 {
                dchord[0] -= eta * chord[1];
                dchord[1] += eta * chord[0];
            }
            electrical_action[k] = pair.from_angle(k, dchord);
            let chord = self.chord(pair.corner[k][0], pair.corner[k][1])?;
            let dchord = self.chord_action(
                pair.corner[k][0],
                pair.corner[k][1],
                corner_action[k][0],
                corner_action[k][1],
            )?;
            let middle = pair.middle[k];
            let dmiddle = [
                gamma.mul_add(middle[1], if k == 1 { middle_direction[0] } else { 0.0 }),
                (-gamma).mul_add(middle[0], if k == 1 { middle_direction[1] } else { 0.0 }),
            ];
            angular[k] = middle[1].mul_add(chord[0], -middle[0] * chord[1]);
            angular_action[k] = dmiddle[1] * chord[0] + middle[1] * dchord[0]
                - dmiddle[0] * chord[1]
                - middle[0] * dchord[1];
        }
        let result = [
            f64::midpoint(electrical_action[0][0], electrical_action[1][0]),
            f64::midpoint(electrical_action[0][1], electrical_action[1][1]),
            self.pole_pairs
                * (pair.sinc_derivative() * gamma * f64::midpoint(angular[0], angular[1])
                    + pair.sinc_half * f64::midpoint(angular_action[0], angular_action[1])),
            0.5 * self.inverse_inertia[0] * direction[3],
            0.5 * self.inverse_inertia[1] * direction[4],
        ];
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    fn energy_change(&self, a: &[f64], b: &[f64]) -> Option<Dd> {
        let a_flux = self.coordinate(self.physical_flux(a)?);
        let b_flux = self.coordinate(self.physical_flux(b)?);
        let mut change = Dd::from_f64(self.curves[0].energy_change_j(a_flux[0], b_flux[0]).ok()?)
            + Dd::from_f64(self.curves[1].energy_change_j(a_flux[1], b_flux[1]).ok()?);
        for k in 0..2 {
            let i = k + 3;
            let mean = Dd::from_f64(0.5 * a[i]) + Dd::from_f64(0.5 * b[i]);
            change =
                change + Dd::from_f64(self.inverse_inertia[k]) * mean * Dd::from_f64(b[i] - a[i]);
        }
        Some(change)
    }
}

impl Storage for SaturatedStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(self.clone()))
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }

    fn hamiltonian(&self, state: &[f64]) -> f64 {
        let Some(flux) = self.physical_flux(state).map(|v| self.coordinate(v)) else {
            return f64::NAN;
        };
        let Ok(d) = self.curves[0].stored_energy_j(flux[0]) else {
            return f64::NAN;
        };
        let Ok(q) = self.curves[1].stored_energy_j(flux[1]) else {
            return f64::NAN;
        };
        d + q
            + 0.5
                * (state[3] * (state[3] * self.inverse_inertia[0])
                    + state[4] * (state[4] * self.inverse_inertia[1]))
    }

    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != 5 {
            return;
        }
        let Some(flux) = self.physical_flux(state) else {
            return;
        };
        let coordinate = self.coordinate(flux);
        let (Ok(d), Ok(q)) = (
            self.curves[0].current_a(coordinate[0]),
            self.curves[1].current_a(coordinate[1]),
        ) else {
            return;
        };
        let (s, c) = det::sin_cos(self.pole_pairs * state[2]);
        out[..2].copy_from_slice(&to_ab(s, c, [d, q]));
        out[2] = self.pole_pairs * flux[1].mul_add(d, -flux[0] * q);
        out[3] = state[3] * self.inverse_inertia[0];
        out[4] = state[4] * self.inverse_inertia[1];
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], middle: &[f64]) -> Option<f64> {
        let Some(mut remainder) = self.energy_change(a, b) else {
            return Some(f64::NAN);
        };
        let mut gradient = [0.0; 5];
        self.gradient(middle, &mut gradient);
        for k in 0..5 {
            remainder = remainder - Dd::from_f64(gradient[k]) * Dd::from_f64(b[k] - a[k]);
        }
        Some(remainder.to_f64())
    }
}

impl ConstitutiveDiscreteGradient for SaturatedStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], _middle: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() == 5 {
            if let Some(effort) = self.effort(a, b) {
                out.copy_from_slice(&effort);
            }
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if out.len() != 5 {
            return false;
        }
        let Some(action) = self.action(a, b, direction) else {
            return false;
        };
        out.copy_from_slice(&action);
        true
    }
}

impl ContinuousStorageTangent for SaturatedStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        let Some(flux) = self.physical_flux(state).map(|v| self.coordinate(v)) else {
            return false;
        };
        if self
            .curves
            .iter()
            .zip(flux)
            .any(|(curve, value)| curve.scalar.continuous_slope(value).is_none())
        {
            return false;
        }
        // At a differentiable source point, symmetry gives H'' = 2 D_2 DG.
        // Newton chords use the scalar source's one-sided endpoint convention;
        // a true continuous tangent must instead refuse a kink or support edge.
        if !self.end_jacobian_action(state, state, direction, out) {
            return false;
        }
        for value in out.iter_mut() {
            *value *= 2.0;
        }
        out.iter().all(|v| v.is_finite())
    }
}
