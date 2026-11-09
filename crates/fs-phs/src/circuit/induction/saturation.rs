//! Isotropic saturation on the original coupled winding fluxes.
//!
//! For positive leakages ls,lr, let a=lr/(ls+lr), b=ls/(ls+lr),
//! ell=ls*lr/(ls+lr), d=lambda_s-R*lambda_r and y=a*lambda_s+b*R*lambda_r.
//! Eliminating the common magnetizing flux gives y=psi_m+ell*i_m. Transforming
//! each original piecewise-affine knot therefore gives the exact I_m(|y|) law,
//! without an inner nonlinear solve. The magnetic energy is
//! |d|^2/(2*(ls+lr)) + integral_0^|y| I_m(u) du. This is exactly the original
//! leakage energy plus integral I_m dpsi_m, not a secant-inductance energy.
use super::{
    InductionMotorParameters, invalid,
    storage::{InductionStorage, rotate},
};
use crate::circuit::inductance::InductorFluxCurve;
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::{dd::Dd, det};

#[derive(Clone)]
pub(super) struct SaturatedStorage {
    leakage: InductionStorage,
    effective: InductorFluxCurve,
    weights: [f64; 2],
    pole_pairs: f64,
}

fn complete(values: &[f64]) -> bool {
    values.len() == 7 && values.iter().all(|v| v.is_finite())
}
fn norm(x: [f64; 2]) -> f64 {
    det::hypot(x[0], x[1])
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}

struct Pair {
    sine: f64,
    cosine: f64,
    half_angle: f64,
    sine_half: f64,
    cosine_half: f64,
    sinc_half: f64,
    rotor_middle: [[f64; 2]; 2],
    rotor_corner: [[[f64; 2]; 2]; 2],
    // [fixed flux endpoint][fixed angle endpoint]
    corner: [[[f64; 2]; 2]; 2],
}

impl Pair {
    fn rotate_at(&self, angle: usize, x: [f64; 2]) -> [f64; 2] {
        rotate(
            if angle == 0 {
                -self.sine_half
            } else {
                self.sine_half
            },
            self.cosine_half,
            rotate(self.sine, self.cosine, x),
        )
    }
    fn inverse_at(&self, angle: usize, x: [f64; 2]) -> [f64; 2] {
        rotate(
            if angle == 0 {
                self.sine_half
            } else {
                -self.sine_half
            },
            self.cosine_half,
            rotate(-self.sine, self.cosine, x),
        )
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
        p: InductionMotorParameters,
        curve: &InductorFluxCurve,
    ) -> Result<Self, PhsError> {
        if curve
            .knots()
            .iter()
            .zip(curve.knots().iter().rev())
            .any(|(a, b)| a.0 != -b.0 || a.1 != -b.1)
        {
            return Err(invalid(
                "isotropic induction magnetizing knots must be odd-symmetric",
            ));
        }
        let ls = p.stator_leakage_inductance_h;
        let lr = p.rotor_leakage_inductance_h;
        let scale = ls.max(lr);
        let denominator = ls / scale + lr / scale;
        let weights = [(lr / scale) / denominator, (ls / scale) / denominator];
        let leakage = InductionStorage::reduced_leakage(p, (1.0 / scale) / denominator)?;
        let ell = ls * weights[0];
        if !ell.is_finite() || ell <= 0.0 || weights.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid(
                "induction magnetizing leakage reduction is unrepresentable",
            ));
        }
        let effective = InductorFluxCurve::try_new(
            curve
                .knots()
                .iter()
                .map(|&(flux, current)| (ell.mul_add(current, flux), current))
                .collect(),
        )
        .map_err(|_| invalid("induction transformed magnetizing law is unrepresentable"))?;
        Ok(Self {
            leakage,
            effective,
            weights,
            pole_pairs: f64::from(p.pole_pairs),
        })
    }

    fn weighted(&self, s: [f64; 2], r: [f64; 2]) -> [f64; 2] {
        std::array::from_fn(|i| self.weights[0].mul_add(s[i], self.weights[1] * r[i]))
    }

    fn coordinate(&self, state: &[f64]) -> Option<([f64; 2], [f64; 2])> {
        if !complete(state) {
            return None;
        }
        let (s, c) = det::sin_cos(self.pole_pairs * state[4]);
        let rotor = rotate(s, c, [state[2], state[3]]);
        let y = self.weighted([state[0], state[1]], rotor);
        y.iter()
            .chain(&rotor)
            .all(|v| v.is_finite())
            .then_some((y, rotor))
    }

    fn pair(&self, a: &[f64], b: &[f64]) -> Option<Pair> {
        if !complete(a) || !complete(b) {
            return None;
        }
        let angle = self.pole_pairs * f64::midpoint(a[4], b[4]);
        let half_angle = (0.5 * self.pole_pairs) * (b[4] - a[4]);
        let (sine, cosine) = det::sin_cos(angle);
        let (sine_half, cosine_half) = det::sin_cos(half_angle);
        let sinc_half = if half_angle == 0.0 {
            1.0
        } else {
            sine_half / half_angle
        };
        let rotor_middle = [
            rotate(sine, cosine, [a[2], a[3]]),
            rotate(sine, cosine, [b[2], b[3]]),
        ];
        let rotor_corner = rotor_middle.map(|r| {
            [
                rotate(-sine_half, cosine_half, r),
                rotate(sine_half, cosine_half, r),
            ]
        });
        let stator = [[a[0], a[1]], [b[0], b[1]]];
        let corner = std::array::from_fn(|i| rotor_corner[i].map(|r| self.weighted(stator[i], r)));
        if [
            angle,
            half_angle,
            sine,
            cosine,
            sine_half,
            cosine_half,
            sinc_half,
        ]
        .iter()
        .chain(rotor_middle.iter().flatten())
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
            rotor_middle,
            rotor_corner,
            corner,
        })
    }

    // Radial discrete gradient: I_bar(ra,rb)*(a+b)/(ra+rb). The scalar kernel
    // integrates every original segment; equal radii and the origin have their
    // continuous limits. No division by rb-ra or endpoint-energy subtraction.
    fn chord(&self, a: [f64; 2], b: [f64; 2]) -> Option<[f64; 2]> {
        let (ra, rb) = (norm(a), norm(b));
        let mean_radius = f64::midpoint(ra, rb);
        let current = self.effective.discrete_current_a(ra, rb).ok()?;
        let result = if mean_radius == 0.0 {
            [0.0; 2]
        } else {
            std::array::from_fn(|i| current * (f64::midpoint(a[i], b[i]) / mean_radius))
        };
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    fn chord_action(
        &self,
        a: [f64; 2],
        b: [f64; 2],
        da: [f64; 2],
        db: [f64; 2],
    ) -> Option<[f64; 2]> {
        let (ra, rb) = (norm(a), norm(b));
        let mean_radius = f64::midpoint(ra, rb);
        let result = if mean_radius == 0.0 {
            let half_slope = 0.5 / self.effective.differential_inductance_h(0.0).ok()?;
            std::array::from_fn(|i| half_slope.mul_add(da[i], half_slope * db[i]))
        } else {
            let current = self.effective.discrete_current_a(ra, rb).ok()?;
            let c = (0.5 * current) / mean_radius;
            // At a zero endpoint its radial coefficient vanishes identically;
            // differentiating an arbitrary direction of the zero norm is wrong.
            let dra = if ra == 0.0 {
                0.0
            } else {
                (a[0] / ra).mul_add(da[0], (a[1] / ra) * da[1])
            };
            let drb = if rb == 0.0 {
                0.0
            } else {
                (b[0] / rb).mul_add(db[0], (b[1] / rb) * db[1])
            };
            let start = self
                .effective
                .discrete_current_end_derivative(rb, ra)
                .ok()?;
            let end = self
                .effective
                .discrete_current_end_derivative(ra, rb)
                .ok()?;
            let correction = (start - c).mul_add(dra, (end - c) * drb);
            std::array::from_fn(|i| {
                c.mul_add(da[i], c * db[i]) + (f64::midpoint(a[i], b[i]) / mean_radius) * correction
            })
        };
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    fn magnetic_gradient(&self, state: &[f64]) -> Option<[f64; 7]> {
        let (y, rotor) = self.coordinate(state)?;
        let current = self.chord(y, y)?;
        let (s, c) = det::sin_cos(self.pole_pairs * state[4]);
        let back = rotate(-s, c, current);
        Some([
            self.weights[0] * current[0],
            self.weights[0] * current[1],
            self.weights[1] * back[0],
            self.weights[1] * back[1],
            self.weights[1] * self.pole_pairs * cross(rotor, current),
            0.0,
            0.0,
        ])
    }

    fn effort(&self, a: &[f64], b: &[f64]) -> Option<[f64; 7]> {
        let pair = self.pair(a, b)?;
        let mut stator = [[0.0; 2]; 2];
        let mut rotor = [[0.0; 2]; 2];
        let mut angular = [0.0; 2];
        for k in 0..2 {
            stator[k] = self.chord(pair.corner[0][k], pair.corner[1][k])?;
            rotor[k] = pair.inverse_at(k, stator[k]);
            let angular_chord = self.chord(pair.corner[k][0], pair.corner[k][1])?;
            angular[k] = cross(pair.rotor_middle[k], angular_chord);
        }
        let result = [
            self.weights[0] * f64::midpoint(stator[0][0], stator[1][0]),
            self.weights[0] * f64::midpoint(stator[0][1], stator[1][1]),
            self.weights[1] * f64::midpoint(rotor[0][0], rotor[1][0]),
            self.weights[1] * f64::midpoint(rotor[0][1], rotor[1][1]),
            self.weights[1]
                * self.pole_pairs
                * pair.sinc_half
                * f64::midpoint(angular[0], angular[1]),
            0.0,
            0.0,
        ];
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    fn action(&self, a: &[f64], b: &[f64], direction: &[f64]) -> Option<[f64; 7]> {
        if !complete(direction) {
            return None;
        }
        let pair = self.pair(a, b)?;
        let eta = self.pole_pairs * direction[4];
        let gamma = 0.5 * eta;
        let ds = [direction[0], direction[1]];
        let dr = [direction[2], direction[3]];
        let mut corner_action = [[[0.0; 2]; 2]; 2];
        for j in 0..2 {
            corner_action[1][j] = self.weighted(ds, pair.rotate_at(j, dr));
        }
        for k in 0..2 {
            corner_action[k][1][0] -= self.weights[1] * eta * pair.rotor_corner[k][1][1];
            corner_action[k][1][1] += self.weights[1] * eta * pair.rotor_corner[k][1][0];
        }
        let middle_dr = rotate(pair.sine, pair.cosine, dr);
        let mut stator = [[0.0; 2]; 2];
        let mut rotor = [[0.0; 2]; 2];
        let mut angular = [0.0; 2];
        let mut angular_action = [0.0; 2];
        for k in 0..2 {
            let chord = self.chord(pair.corner[0][k], pair.corner[1][k])?;
            stator[k] = self.chord_action(
                pair.corner[0][k],
                pair.corner[1][k],
                corner_action[0][k],
                corner_action[1][k],
            )?;
            let mut dchord = stator[k];
            if k == 1 {
                dchord[0] += eta * chord[1];
                dchord[1] -= eta * chord[0];
            }
            rotor[k] = pair.inverse_at(k, dchord);
            let chord = self.chord(pair.corner[k][0], pair.corner[k][1])?;
            let dchord = self.chord_action(
                pair.corner[k][0],
                pair.corner[k][1],
                corner_action[k][0],
                corner_action[k][1],
            )?;
            let middle = pair.rotor_middle[k];
            let dmiddle = [
                -gamma * middle[1] + if k == 1 { middle_dr[0] } else { 0.0 },
                gamma * middle[0] + if k == 1 { middle_dr[1] } else { 0.0 },
            ];
            angular[k] = cross(middle, chord);
            angular_action[k] = cross(dmiddle, chord) + cross(middle, dchord);
        }
        let result = [
            self.weights[0] * f64::midpoint(stator[0][0], stator[1][0]),
            self.weights[0] * f64::midpoint(stator[0][1], stator[1][1]),
            self.weights[1] * f64::midpoint(rotor[0][0], rotor[1][0]),
            self.weights[1] * f64::midpoint(rotor[0][1], rotor[1][1]),
            self.weights[1]
                * self.pole_pairs
                * (pair.sinc_derivative() * gamma * f64::midpoint(angular[0], angular[1])
                    + pair.sinc_half * f64::midpoint(angular_action[0], angular_action[1])),
            0.0,
            0.0,
        ];
        result.iter().all(|v| v.is_finite()).then_some(result)
    }

    // Independent original endpoint work; the accepted split efforts and the
    // electrical/shaft work and heat reports are not inputs. Trig polarization
    // retains a tiny rotation beside a large flux without subtracting energies.
    fn magnetic_energy_change(&self, a: &[f64], b: &[f64]) -> Option<Dd> {
        let pair = self.pair(a, b)?;
        let chord = self.chord(pair.corner[0][0], pair.corner[1][1])?;
        let dr = rotate(pair.sine, pair.cosine, [b[2] - a[2], b[3] - a[3]]);
        let mr = [
            f64::midpoint(pair.rotor_middle[0][0], pair.rotor_middle[1][0]),
            f64::midpoint(pair.rotor_middle[0][1], pair.rotor_middle[1][1]),
        ];
        let delta_rotor = [
            pair.cosine_half * dr[0] - 2.0 * pair.sine_half * mr[1],
            pair.cosine_half * dr[1] + 2.0 * pair.sine_half * mr[0],
        ];
        let mut change = Dd::from_f64(0.0);
        for k in 0..2 {
            let delta = Dd::from_f64(self.weights[0]) * Dd::from_f64(b[k] - a[k])
                + Dd::from_f64(self.weights[1]) * Dd::from_f64(delta_rotor[k]);
            change = change + Dd::from_f64(chord[k]) * delta;
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
        let Some((y, _)) = self.coordinate(state) else {
            return f64::NAN;
        };
        let Ok(magnetic) = self.effective.stored_energy_j(norm(y)) else {
            return f64::NAN;
        };
        self.leakage.hamiltonian(state) + magnetic
    }
    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != 7 {
            return;
        }
        let Some(magnetic) = self.magnetic_gradient(state) else {
            return;
        };
        self.leakage.gradient(state, out);
        for k in 0..7 {
            out[k] += magnetic[k];
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let (Some(mut change), Some(gradient)) = (
            self.magnetic_energy_change(a, b),
            self.magnetic_gradient(mid),
        ) else {
            return Some(f64::NAN);
        };
        change = change + Dd::from_f64(self.leakage.midpoint_energy_remainder(a, b, mid)?);
        for k in 0..7 {
            change = change - Dd::from_f64(gradient[k]) * Dd::from_f64(b[k] - a[k]);
        }
        Some(change.to_f64())
    }
}

impl ConstitutiveDiscreteGradient for SaturatedStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], mid: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != 7 {
            return;
        }
        let Some(magnetic) = self.effort(a, b) else {
            return;
        };
        self.leakage.evaluate(a, b, mid, out);
        for k in 0..7 {
            out[k] += magnetic[k];
        }
    }
    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if out.len() != 7 {
            return false;
        }
        let Some(magnetic) = self.action(a, b, direction) else {
            return false;
        };
        if !self.leakage.end_jacobian_action(a, b, direction, out) {
            return false;
        }
        for k in 0..7 {
            out[k] += magnetic[k];
        }
        out.iter().all(|v| v.is_finite())
    }
}

impl ContinuousStorageTangent for SaturatedStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        let Some((y, _)) = self.coordinate(state) else {
            return false;
        };
        if self.effective.scalar.continuous_slope(norm(y)).is_none() {
            return false;
        }
        // A symmetric consistent chord has H''=2 D_2 DG at a differentiable
        // point. A source kink/support boundary is not a continuous tangent.
        if !self.end_jacobian_action(state, state, direction, out) {
            return false;
        }
        for v in out.iter_mut() {
            *v *= 2.0;
        }
        out.iter().all(|v| v.is_finite())
    }
}
