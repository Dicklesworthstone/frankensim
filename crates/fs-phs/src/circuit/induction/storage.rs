//! One original reciprocal winding energy; no separately fitted torque law.
use super::{InductionMotorParameters, invalid};
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::{dd::Dd, det};

#[derive(Clone, Copy)]
pub(super) struct InductionStorage {
    // Hmag = (kr |s|² + ks |r|² + km |s-Rr|²)/2. This positive
    // decomposition avoids subtracting two large nearly equal self/mutual terms.
    kr: f64,
    ks: f64,
    km: f64,
    pole_pairs: f64,
    inverse_inertia: [f64; 2],
}

pub(super) fn rotate(s: f64, c: f64, x: [f64; 2]) -> [f64; 2] {
    [c.mul_add(x[0], -s * x[1]), s.mul_add(x[0], c * x[1])]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}
fn complete(x: &[f64]) -> bool {
    x.len() == 7 && x.iter().all(|v| v.is_finite())
}

struct Pair {
    sine: f64,
    cosine: f64,
    half_angle: f64,
    sine_half: f64,
    cosine_half: f64,
    sinc_half: f64,
    ms: [f64; 2],
    mr: [f64; 2],
    ds: [f64; 2],
    dr: [f64; 2],
    rotated_mr: [f64; 2],
    rotated_dr: [f64; 2],
    inverse_ms: [f64; 2],
}

impl InductionStorage {
    // The exact positive leakage part of the saturated winding elimination.
    // This reuses its original flux/angle work, kinetic energy and tangent.
    pub(super) fn reduced_leakage(
        p: InductionMotorParameters,
        inverse_sum: f64,
    ) -> Result<Self, PhsError> {
        let mut storage = Self::new(p)?;
        if !inverse_sum.is_finite() || inverse_sum <= 0.0 {
            return Err(invalid("induction reduced leakage is unrepresentable"));
        }
        storage.kr = 0.0;
        storage.ks = 0.0;
        storage.km = inverse_sum;
        Ok(storage)
    }

    pub(super) fn new(p: InductionMotorParameters) -> Result<Self, PhsError> {
        let ls = p.stator_leakage_inductance_h;
        let lr = p.rotor_leakage_inductance_h;
        let lm = p.magnetizing_inductance_h;
        let inverse_inertia = [1.0 / p.rotor_inertia_kg_m2, 1.0 / p.stator_inertia_kg_m2];
        if [ls, lr, p.rotor_inertia_kg_m2, p.stator_inertia_kg_m2]
            .iter()
            .chain(&inverse_inertia)
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !lm.is_finite()
            || lm < 0.0
            || p.pole_pairs == 0
        {
            return Err(invalid(
                "induction needs positive finite leakage/inertias and pole pairs, and nonnegative magnetizing inductance",
            ));
        }
        // Scale before forming Lls Llr + Lm(Lls+Llr), avoiding a determinant
        // cancellation and avoidable overflow of the unscaled products.
        let scale = ls.max(lr).max(lm);
        let (s, r, m) = (ls / scale, lr / scale, lm / scale);
        let determinant = s.mul_add(r, m * (s + r));
        let coefficients = [r, s, m].map(|v| (v / determinant) / scale);
        if !determinant.is_finite()
            || determinant <= 0.0
            || coefficients.iter().any(|v| !v.is_finite())
            || coefficients[0] <= 0.0
            || coefficients[1] <= 0.0
            || (lm > 0.0 && coefficients[2] <= 0.0)
            || !(ls + lm).is_finite()
            || !(lr + lm).is_finite()
        {
            return Err(invalid(
                "induction original winding constitutive coefficients are unrepresentable",
            ));
        }
        Ok(Self {
            kr: coefficients[0],
            ks: coefficients[1],
            km: coefficients[2],
            pole_pairs: f64::from(p.pole_pairs),
            inverse_inertia,
        })
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
        let ms = [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])];
        let mr = [f64::midpoint(a[2], b[2]), f64::midpoint(a[3], b[3])];
        let ds = [b[0] - a[0], b[1] - a[1]];
        let dr = [b[2] - a[2], b[3] - a[3]];
        let rotated_mr = rotate(sine, cosine, mr);
        let rotated_dr = rotate(sine, cosine, dr);
        let inverse_ms = rotate(-sine, cosine, ms);
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
        .chain(&ms)
        .chain(&mr)
        .chain(&ds)
        .chain(&dr)
        .chain(&rotated_mr)
        .chain(&rotated_dr)
        .chain(&inverse_ms)
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
            ms,
            mr,
            ds,
            dr,
            rotated_mr,
            rotated_dr,
            inverse_ms,
        })
    }

    fn angular_product(pair: &Pair) -> f64 {
        cross(pair.ms, pair.rotated_mr) + 0.25 * cross(pair.ds, pair.rotated_dr)
    }

    // Polarize the ORIGINAL positive energy decomposition independently of
    // accepted effort, work and heat. Rotated endpoint means/differences use
    // trig identities instead of subtracting nearly equal stored energies.
    fn energy_change(&self, a: &[f64], b: &[f64], pair: &Pair) -> Dd {
        let mut energy = Dd::from_f64(0.0);
        for i in 0..2 {
            let j = 1 - i;
            let sign = if i == 0 { -1.0 } else { 1.0 };
            let mean_difference = pair.ms[i]
                - pair.cosine_half * pair.rotated_mr[i]
                - 0.5 * pair.sine_half * sign * pair.rotated_dr[j];
            let delta_difference = pair.ds[i]
                - pair.cosine_half * pair.rotated_dr[i]
                - 2.0 * pair.sine_half * sign * pair.rotated_mr[j];
            energy = energy
                + Dd::from_f64(self.kr) * Dd::from_f64(pair.ms[i]) * Dd::from_f64(pair.ds[i])
                + Dd::from_f64(self.ks) * Dd::from_f64(pair.mr[i]) * Dd::from_f64(pair.dr[i])
                + Dd::from_f64(self.km)
                    * Dd::from_f64(mean_difference)
                    * Dd::from_f64(delta_difference)
                + Dd::from_f64(self.inverse_inertia[i])
                    * (Dd::from_f64(0.5 * a[5 + i]) + Dd::from_f64(0.5 * b[5 + i]))
                    * Dd::from_f64(b[5 + i] - a[5 + i]);
        }
        energy
    }
}

impl Storage for InductionStorage {
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
        let (s, c) = det::sin_cos(self.pole_pairs * x[4]);
        let r = rotate(s, c, [x[2], x[3]]);
        let mut energy = 0.0;
        for i in 0..2 {
            let difference = x[i] - r[i];
            energy += x[i] * (self.kr * x[i])
                + x[i + 2] * (self.ks * x[i + 2])
                + difference * (self.km * difference)
                + x[i + 5] * (self.inverse_inertia[i] * x[i + 5]);
        }
        0.5 * energy
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.evaluate(x, x, x, out);
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let Some(pair) = self.pair(a, b) else {
            return Some(f64::NAN);
        };
        let mut gradient = [0.0; 7];
        self.gradient(mid, &mut gradient);
        let mut remainder = self.energy_change(a, b, &pair);
        for i in 0..7 {
            remainder = remainder - Dd::from_f64(gradient[i]) * Dd::from_f64(b[i] - a[i]);
        }
        Some(remainder.to_f64())
    }
}

impl ConstitutiveDiscreteGradient for InductionStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], _mid: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if out.len() != 7 {
            return;
        }
        let Some(pair) = self.pair(a, b) else {
            return;
        };
        // Symmetric split between all four fluxes and the mechanical angle:
        // average quadratic flux chords at both angles and exact angular work
        // at both flux endpoints. sinc gives the continuous zero-angle limit.
        for i in 0..2 {
            out[i] = self.kr * pair.ms[i]
                + self.km * (-pair.cosine_half).mul_add(pair.rotated_mr[i], pair.ms[i]);
            out[i + 2] = self.ks * pair.mr[i]
                + self.km * (-pair.cosine_half).mul_add(pair.inverse_ms[i], pair.mr[i]);
            out[i + 5] = f64::midpoint(a[i + 5], b[i + 5]) * self.inverse_inertia[i];
        }
        out[4] = self.pole_pairs * self.km * pair.sinc_half * Self::angular_product(&pair);
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !complete(direction) || out.len() != 7 {
            return false;
        }
        let Some(pair) = self.pair(a, b) else {
            return false;
        };
        let dangle = (0.5 * self.pole_pairs) * direction[4];
        let dms = [0.5 * direction[0], 0.5 * direction[1]];
        let dmr = [0.5 * direction[2], 0.5 * direction[3]];
        let dds = [direction[0], direction[1]];
        let ddr = [direction[2], direction[3]];
        let mut dv = rotate(pair.sine, pair.cosine, dmr);
        let mut dw = rotate(pair.sine, pair.cosine, ddr);
        let mut du = rotate(-pair.sine, pair.cosine, dms);
        for i in 0..2 {
            let j = 1 - i;
            let sign = if i == 0 { -1.0 } else { 1.0 };
            dv[i] += sign * dangle * pair.rotated_mr[j];
            dw[i] += sign * dangle * pair.rotated_dr[j];
            du[i] -= sign * dangle * pair.inverse_ms[j];
            out[i] = self.kr * dms[i]
                + self.km
                    * (dms[i] - pair.cosine_half * dv[i]
                        + pair.sine_half * dangle * pair.rotated_mr[i]);
            out[i + 2] = self.ks * dmr[i]
                + self.km
                    * (dmr[i] - pair.cosine_half * du[i]
                        + pair.sine_half * dangle * pair.inverse_ms[i]);
            out[i + 5] = 0.5 * self.inverse_inertia[i] * direction[i + 5];
        }
        let product = Self::angular_product(&pair);
        let dproduct = cross(dms, pair.rotated_mr)
            + cross(pair.ms, dv)
            + 0.25 * (cross(dds, pair.rotated_dr) + cross(pair.ds, dw));
        let sinc_derivative = if pair.half_angle.abs() < 0.01 {
            let square = pair.half_angle * pair.half_angle;
            pair.half_angle
                * (-1.0 / 3.0 + square * (1.0 / 30.0 + square * (-1.0 / 840.0 + square / 45360.0)))
        } else {
            (pair.cosine_half - pair.sinc_half) / pair.half_angle
        };
        out[4] = self.pole_pairs
            * self.km
            * (sinc_derivative * dangle * product + pair.sinc_half * dproduct);
        out.iter().all(|v| v.is_finite())
    }
}

impl ContinuousStorageTangent for InductionStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        if !self.end_jacobian_action(state, state, direction, out) {
            return false;
        }
        for v in out.iter_mut() {
            *v *= 2.0;
        }
        out.iter().all(|v| v.is_finite())
    }
}
