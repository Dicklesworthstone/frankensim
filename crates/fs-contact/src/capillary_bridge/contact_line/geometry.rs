//! Fixed-volume circular geometry and analytic energy derivatives.
use super::{ContactLineError, ContactLineProfile, SourceContactLineBridge, poll};
use crate::capillary_bridge::BridgeProfile;
use fs_exec::CancelGate;
use fs_math::det;

pub(super) struct Evaluation {
    pub profile: ContactLineProfile,
    hh: f64,
    hb: f64,
    bb: f64,
    normal_geometry_error_n: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct Gradient {
    pub value: f64,
    pub dh: f64,
    pub db: f64,
    pub error: f64,
}

impl SourceContactLineBridge {
    pub(super) fn span_at_angle(&self, h: f64, theta: f64) -> Result<f64, ContactLineError> {
        self.original.check_gap(h)?;
        if !theta.is_finite()
            || theta < self.config.minimum_angle_radians
            || theta > self.config.maximum_angle_radians
        {
            return Err(ContactLineError::Domain("admitted acute contact angle"));
        }
        let (k, _, _, _, _) = shape(theta);
        let span = self.original.geometry.area / h + h * k;
        if !span.is_finite() || span <= 0.0 {
            return Err(ContactLineError::Domain("representable wet footprint"));
        }
        Ok(span)
    }

    pub(super) fn at_angle(
        &self,
        h: f64,
        b: f64,
        theta: f64,
    ) -> Result<Evaluation, ContactLineError> {
        let implied_span = self.span_at_angle(h, theta)?;
        if (theta - self.original.geometry.theta).abs()
            > self.config.maximum_angle_departure_radians
        {
            return Err(ContactLineError::Domain(
                "linear contact-line law angle departure",
            ));
        }
        let span_residual = implied_span - b;
        let tolerance = self.config.geometry_absolute_tolerance_m
            + self.config.geometry_relative_tolerance * b.abs();
        if !b.is_finite() || b <= 0.0 || span_residual.abs() > tolerance {
            return Err(ContactLineError::NotConverged);
        }
        let (k, dk, sine, cosine, alpha) = shape(theta);
        let neck = b - h * cosine / (1.0 + sine);
        if !neck.is_finite() || neck < self.config.minimum_neck_width_m {
            return Err(ContactLineError::Domain(
                "connected resolved neck; no pinch-off model",
            ));
        }
        let g = self.original.geometry;
        let scale = 2.0 * g.sigma_width;
        let theta_b = 1.0 / (h * dk);
        let theta_h = ((g.area / h) / h - k) / (h * dk);
        let energy = scale * (h * (alpha / cosine) - g.cosine * b);
        let uh = scale * (cosine * b / h + sine);
        let ub = scale * (cosine - g.cosine);
        let hh = scale * ((cosine - sine * b / h) * theta_h - (cosine * b / h) / h);
        let hb = -scale * sine * theta_h;
        let bb = -scale * sine * theta_b;
        let profile = ContactLineProfile {
            gap_m: h,
            wetted_span_m: b,
            contact_angle_radians: theta,
            neck_width_m: neck,
            liquid_volume_m3: g.volume,
            liquid_gas_area_m2: 2.0 * h * (alpha / cosine) * g.width,
            wetted_area_m2: 2.0 * b * g.width,
            liquid_minus_gas_pressure_pa: -2.0 * (g.sigma_width / g.width) * cosine / h,
            surface_free_energy_j: energy,
            normal_force_n: -uh,
            footprint_energy_gradient_n: ub,
            geometric_span_residual_m: span_residual,
            estimated_volume_error_m3: g.width * h * span_residual.abs(),
            estimated_surface_energy_error_j: scale * cosine * span_residual.abs(),
        };
        if [
            energy,
            uh,
            ub,
            hh,
            hb,
            bb,
            profile.liquid_gas_area_m2,
            profile.wetted_area_m2,
            profile.liquid_minus_gas_pressure_pa,
            profile.estimated_volume_error_m3,
            profile.estimated_surface_energy_error_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || bb <= 0.0
        {
            return Err(ContactLineError::Domain(
                "finite convex footprint energy geometry",
            ));
        }
        self.original.check_profile(&BridgeProfile {
            gap_m: h,
            attached: true,
            liquid_volume_m3: g.volume,
            detached_cap_volumes_m3: None,
            wetted_span_m: b,
            neck_width_m: Some(neck),
            liquid_gas_area_m2: profile.liquid_gas_area_m2,
            wetted_area_m2: profile.wetted_area_m2,
            liquid_minus_gas_pressure_pa: profile.liquid_minus_gas_pressure_pa,
            surface_free_energy_j: energy,
            normal_force_n: -uh,
            force_gap_derivative_n_per_m: -hh,
        })?;
        // Inverse geometry perturbs theta while the canonical footprint is
        // held fixed. Include the explicit b contribution when bounding the
        // angle-only force error; the total mixed Hessian alone can cancel it.
        let normal_geometry_error_n = (hb.abs() + scale * cosine / h) * span_residual.abs();
        Ok(Evaluation {
            profile,
            hh,
            hb,
            bb,
            normal_geometry_error_n,
        })
    }

    pub(super) fn at(
        &self,
        h: f64,
        b: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        poll(gate)?;
        if *evaluations >= self.config.maximum_evaluations {
            return Err(ContactLineError::Budget);
        }
        *evaluations += 1;
        let mut lower = self.config.minimum_angle_radians;
        let mut upper = self.config.maximum_angle_radians;
        let low_span = self.span_at_angle(h, upper)?;
        let high_span = self.span_at_angle(h, lower)?;
        let tolerance = self.config.geometry_absolute_tolerance_m
            + self.config.geometry_relative_tolerance * b.abs();
        if !b.is_finite() || b < low_span - tolerance || b > high_span + tolerance {
            return Err(ContactLineError::Domain(
                "footprint outside admitted circular branch",
            ));
        }
        let mut theta = 0.5 * (lower + upper);
        for _ in 0..self.config.maximum_geometry_iterations {
            poll(gate)?;
            let (k, dk, _, _, _) = shape(theta);
            let residual = self.original.geometry.area / h + h * k - b;
            if residual.abs() <= tolerance {
                return self.at_angle(h, b, theta);
            }
            if residual > 0.0 {
                lower = theta;
            } else {
                upper = theta;
            }
            let candidate = theta - residual / (h * dk);
            let next = if candidate > lower && candidate < upper {
                candidate
            } else {
                lower + 0.5 * (upper - lower)
            };
            if next == theta {
                break;
            }
            theta = next;
        }
        Err(ContactLineError::NotConverged)
    }

    pub(super) fn footprint_gradient(
        &self,
        h0: f64,
        h1: f64,
        b0: f64,
        b1: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let left = self.integral(false, h0, b0, b1, evaluations, gate)?;
        let right = if h0 == h1 {
            left
        } else {
            self.integral(false, h1, b0, b1, evaluations, gate)?
        };
        Ok(Gradient {
            value: 0.5 * (left.value + right.value),
            db: 0.5 * (left.db + right.db),
            dh: 0.5 * right.dh,
            error: 0.5 * (left.error + right.error),
        })
    }

    pub(super) fn gap_gradient(
        &self,
        h0: f64,
        h1: f64,
        b0: f64,
        b1: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let left = self.integral(true, b0, h0, h1, evaluations, gate)?;
        let right = if b0 == b1 {
            left
        } else {
            self.integral(true, b1, h0, h1, evaluations, gate)?
        };
        Ok(Gradient {
            value: 0.5 * (left.value + right.value),
            dh: 0.5 * (left.dh + right.dh),
            db: 0.5 * right.db,
            error: 0.5 * (left.error + right.error),
        })
    }

    fn integral(
        &self,
        along_gap: bool,
        fixed: f64,
        from: f64,
        to: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        if from == to {
            let p = if along_gap {
                self.at(from, fixed, evaluations, gate)?
            } else {
                self.at(fixed, from, evaluations, gate)?
            };
            return Ok(sample(&p, along_gap, 0.5));
        }
        let coarse = self.quadrature(along_gap, fixed, from, to, &GAUSS8, evaluations, gate)?;
        let mut fine = self.quadrature(along_gap, fixed, from, to, &GAUSS16, evaluations, gate)?;
        fine.error += (fine.value - coarse.value).abs() + 64.0 * f64::EPSILON * fine.value.abs();
        Ok(fine)
    }

    fn quadrature(
        &self,
        along_gap: bool,
        fixed: f64,
        from: f64,
        to: f64,
        nodes: &[(f64, f64)],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Gradient, ContactLineError> {
        let mut sum = Gradient::default();
        for &(x, w) in nodes {
            for sign in [-1.0, 1.0] {
                let t = 0.5 * (1.0 + sign * x);
                let value = from + t * (to - from);
                let p = if along_gap {
                    self.at(value, fixed, evaluations, gate)?
                } else {
                    self.at(fixed, value, evaluations, gate)?
                };
                let q = sample(&p, along_gap, t);
                sum.value += 0.5 * w * q.value;
                sum.dh += 0.5 * w * q.dh;
                sum.db += 0.5 * w * q.db;
                sum.error += 0.5 * w * q.error;
            }
        }
        Ok(sum)
    }
}

fn sample(p: &Evaluation, along_gap: bool, t: f64) -> Gradient {
    let residual = p.profile.geometric_span_residual_m.abs();
    if along_gap {
        Gradient {
            value: -p.profile.normal_force_n,
            dh: t * p.hh,
            db: p.hb,
            error: p.normal_geometry_error_n + 32.0 * f64::EPSILON * p.profile.normal_force_n.abs(),
        }
    } else {
        Gradient {
            value: p.profile.footprint_energy_gradient_n,
            dh: p.hb,
            db: t * p.bb,
            error: p.bb * residual
                + 32.0 * f64::EPSILON * p.profile.footprint_energy_gradient_n.abs(),
        }
    }
}

fn shape(theta: f64) -> (f64, f64, f64, f64, f64) {
    let sine = det::sin(theta);
    let cosine = det::cos(theta);
    let alpha = det::atan2(cosine, sine);
    let numerator = if alpha < 0.01 {
        let a2 = alpha * alpha;
        alpha
            * a2
            * (2.0 / 3.0
                + a2 * (-2.0 / 15.0
                    + a2 * (4.0 / 315.0 + a2 * (-2.0 / 2835.0 + a2 * 4.0 / 155925.0))))
    } else {
        alpha - sine * cosine
    };
    let k = numerator / (2.0 * cosine * cosine);
    (k, -1.0 + 2.0 * k * sine / cosine, sine, cosine, alpha)
}

const GAUSS8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_26),
];
const GAUSS16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_54),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_73),
    (0.755_404_408_355_003, 0.124_628_971_255_533_87),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_095),
];
