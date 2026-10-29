//! Original fixed-volume concave circular menisci with independent plate spans.
use super::{
    AsymmetricContactLineProfile, ContactLineError, SourceAsymmetricContactLineBridge, poll,
};
use fs_exec::CancelGate;
use fs_math::det;

pub(super) struct Evaluation {
    pub profile: AsymmetricContactLineProfile,
    /// Derivatives in [gap, lower full span, upper full span].
    pub gradient: [f64; 3],
    pub hessian: [[f64; 3]; 3],
    pub gradient_error: [f64; 3],
}

impl SourceAsymmetricContactLineBridge {
    pub(super) fn from_angles(
        &self,
        h: f64,
        theta: [f64; 2],
    ) -> Result<Evaluation, ContactLineError> {
        self.check_angles(theta)?;
        self.check_gap(h)?;
        let phi = core::f64::consts::FRAC_PI_2 - 0.5 * (theta[0] + theta[1]);
        let alpha = 0.5 * (theta[0] - theta[1]);
        let d = h * det::sin(alpha) / det::cos(alpha);
        let p = h * h + d * d;
        let s = (self.area + p * segment(phi).0) / h;
        self.at_phi(h, [s - d, s + d], phi)
    }

    pub(super) fn at(
        &self,
        h: f64,
        spans: [f64; 2],
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        poll(gate)?;
        let c = self.config.contact_line;
        if *evaluations >= c.maximum_evaluations {
            return Err(ContactLineError::Budget);
        }
        *evaluations += 1;
        self.check_gap(h)?;
        if spans
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x > self.limits.maximum_wetted_span_m)
        {
            return Err(ContactLineError::Domain(
                "positive admitted per-surface wet footprint",
            ));
        }
        let d = 0.5 * (spans[1] - spans[0]);
        let s = 0.5 * (spans[0] + spans[1]);
        let p = h * h + d * d;
        let alpha = det::atan2(d, h);
        let mut low = 0.0_f64;
        let mut high = core::f64::consts::FRAC_PI_2;
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let young = self.wetting[i].contact_angle_radians();
            let min = c
                .minimum_angle_radians
                .max(young - c.maximum_angle_departure_radians);
            let max = c
                .maximum_angle_radians
                .min(young + c.maximum_angle_departure_radians);
            low = low.max(core::f64::consts::FRAC_PI_2 + sign * alpha - max);
            high = high.min(core::f64::consts::FRAC_PI_2 + sign * alpha - min);
        }
        if !(0.0 < low && low < high && high < core::f64::consts::FRAC_PI_2) {
            return Err(ContactLineError::Domain(
                "common concave acute-angle circular branch",
            ));
        }
        let area_excess = h * s - self.area;
        let tolerance = self.area_tolerance(h, spans);
        if !p.is_finite()
            || !area_excess.is_finite()
            || area_excess < p * segment(low).0 - tolerance
            || area_excess > p * segment(high).0 + tolerance
        {
            return Err(ContactLineError::Domain(
                "fixed-volume footprint outside circular branch",
            ));
        }
        let mut phi = 0.5 * (low + high);
        for _ in 0..c.maximum_geometry_iterations {
            poll(gate)?;
            let (q, dq) = segment(phi);
            let residual = p * q - area_excess;
            if residual.abs() <= tolerance {
                return self.at_phi(h, spans, phi);
            }
            if residual < 0.0 {
                low = phi;
            } else {
                high = phi;
            }
            let candidate = phi - residual / (p * dq);
            let next = if candidate > low && candidate < high {
                candidate
            } else {
                low + 0.5 * (high - low)
            };
            if next == phi {
                break;
            }
            phi = next;
        }
        Err(ContactLineError::NotConverged)
    }

    fn at_phi(&self, h: f64, spans: [f64; 2], phi: f64) -> Result<Evaluation, ContactLineError> {
        self.check_gap(h)?;
        if spans
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x > self.limits.maximum_wetted_span_m)
        {
            return Err(ContactLineError::Domain(
                "positive admitted per-surface wet footprint",
            ));
        }
        let d = 0.5 * (spans[1] - spans[0]);
        let s = 0.5 * (spans[0] + spans[1]);
        let length = h.hypot(d);
        let p = length * length;
        let alpha = det::atan2(d, h);
        let theta = [
            core::f64::consts::FRAC_PI_2 - phi + alpha,
            core::f64::consts::FRAC_PI_2 - phi - alpha,
        ];
        self.check_angles(theta)?;
        let (q, dq) = segment(phi);
        let residual = p * q - (h * s - self.area);
        if residual.abs() > self.area_tolerance(h, spans) {
            return Err(ContactLineError::NotConverged);
        }
        let sin = det::sin(phi);
        let cos = det::cos(phi);
        let radius = length / (2.0 * sin);
        let half_sine = det::sin(0.5 * phi);
        let neck = s - (2.0 * h * half_sine * half_sine + d * d / (length + h)) / sin;
        let neck_height = 0.5 * (h - d * cos / sin);
        if !neck.is_finite()
            || neck < self.config.contact_line.minimum_neck_width_m
            || !(neck_height > 0.0 && neck_height < h)
        {
            return Err(ContactLineError::Domain(
                "connected interior neck; no pinch-off model",
            ));
        }
        let dp = [2.0 * h, -d, d];
        let ds = [s, 0.5 * h, 0.5 * h];
        let dphi: [f64; 3] = core::array::from_fn(|j| (ds[j] - q * dp[j]) / (p * dq));
        let dalpha = [-d / p, -0.5 * h / p, 0.5 * h / p];
        let theta_j: [[_; 3]; 2] = core::array::from_fn(|i| {
            core::array::from_fn(|j| -dphi[j] + if i == 0 { dalpha[j] } else { -dalpha[j] })
        });
        let gamma = self.gamma_width;
        let c = sin * s + cos * h;
        let uh = 2.0 * gamma * c / length;
        let ub = core::array::from_fn(|i| gamma * (det::cos(theta[i]) - self.young_cosines[i]));
        let mut hessian = [[0.0; 3]; 3];
        for j in 0..3 {
            let dc = (cos * s - sin * h) * dphi[j] + if j == 0 { cos } else { 0.5 * sin };
            hessian[0][j] = 2.0 * gamma / length * (dc - c * dp[j] / (2.0 * p));
            for i in 0..2 {
                hessian[i + 1][j] = -gamma * det::sin(theta[i]) * theta_j[i][j];
            }
        }
        let determinant = hessian[1][1] * hessian[2][2] - hessian[1][2] * hessian[2][1];
        if hessian.iter().flatten().any(|v| !v.is_finite())
            || hessian[1][1] <= 0.0
            || hessian[2][2] <= 0.0
            || !determinant.is_finite()
            || determinant <= 0.0
        {
            return Err(ContactLineError::Domain(
                "positive reciprocal footprint energy Hessian",
            ));
        }
        let free_length = 2.0 * length * phi / sin;
        let energy = gamma
            * (free_length - self.young_cosines[0] * spans[0] - self.young_cosines[1] * spans[1]);
        let pressure = -self.gamma / radius;
        let forces = [
            -(pressure * spans[0] * self.width - 2.0 * gamma * det::sin(theta[0])),
            pressure * spans[1] * self.width - 2.0 * gamma * det::sin(theta[1]),
        ];
        let phi_error = residual.abs() / (p * dq);
        // The signed surface energy subtracts two independently sourced
        // wetting terms from the free surface. Near cancellation its residual
        // alone understates arithmetic uncertainty, even with an exact root.
        let energy_scale = gamma
            * (free_length.abs()
                + (self.young_cosines[0] * spans[0]).abs()
                + (self.young_cosines[1] * spans[1]).abs());
        let energy_error = (2.0 * gamma * length * (sin - phi * cos) / (sin * sin)).abs()
            * phi_error
            + 64.0 * f64::EPSILON * energy_scale;
        let gradient = [uh, ub[0], ub[1]];
        let gradient_error = [
            (2.0 * gamma / length * (cos * s - sin * h)).abs() * phi_error,
            gamma * det::sin(theta[0]) * phi_error,
            gamma * det::sin(theta[1]) * phi_error,
        ];
        let gradient_error =
            core::array::from_fn(|i| gradient_error[i] + 64.0 * f64::EPSILON * gradient[i].abs());
        let profile = AsymmetricContactLineProfile {
            gap_m: h,
            wetted_spans_m: spans,
            contact_angles_radians: theta,
            meniscus_radius_m: radius,
            neck_width_m: neck,
            neck_height_m: neck_height,
            liquid_volume_m3: self.volume,
            liquid_gas_area_m2: free_length * self.width,
            wetted_areas_m2: spans.map(|b| b * self.width),
            liquid_minus_gas_pressure_pa: pressure,
            surface_free_energy_j: energy,
            normal_force_n: -uh,
            plate_normal_forces_n: forces,
            footprint_energy_gradients_n: ub,
            geometric_area_residual_m2: residual,
            estimated_volume_error_m3: self.width * residual.abs(),
            estimated_surface_energy_error_j: energy_error,
        };
        if gradient
            .iter()
            .chain(gradient_error.iter())
            .any(|v| !v.is_finite())
            || !energy.is_finite()
            || !energy_error.is_finite()
            || !radius.is_finite()
            || !profile.liquid_gas_area_m2.is_finite()
            || !pressure.is_finite()
        {
            return Err(ContactLineError::Domain(
                "finite circular geometry and original surface work",
            ));
        }
        super::super::limit(
            uh.abs(),
            self.limits.maximum_abs_force_n,
            "asymmetric bridge normal force",
        )?;
        super::super::limit(
            energy.abs(),
            self.limits.maximum_abs_surface_energy_j,
            "asymmetric surface energy",
        )?;
        let scale = h.max(spans[0]).max(spans[1]);
        let bond =
            self.liquid.density_kg_m3() * self.limits.gravity_m_s2 * scale * scale / self.gamma;
        super::super::limit(
            bond,
            self.limits.maximum_bond_number,
            "asymmetric bridge Bond number",
        )?;
        Ok(Evaluation {
            profile,
            gradient,
            hessian,
            gradient_error,
        })
    }

    fn check_angles(&self, angles: [f64; 2]) -> Result<(), ContactLineError> {
        let c = self.config.contact_line;
        for i in 0..2 {
            if !angles[i].is_finite()
                || angles[i] < c.minimum_angle_radians
                || angles[i] > c.maximum_angle_radians
                || (angles[i] - self.wetting[i].contact_angle_radians()).abs()
                    > c.maximum_angle_departure_radians
            {
                return Err(ContactLineError::Domain(
                    "each original surface's admitted acute dynamic angle",
                ));
            }
        }
        Ok(())
    }
    fn area_tolerance(&self, h: f64, spans: [f64; 2]) -> f64 {
        let c = self.config.contact_line;
        let scale = h.max(spans[0]).max(spans[1]);
        c.geometry_absolute_tolerance_m * scale
            + c.geometry_relative_tolerance * self.area.max(h * 0.5 * (spans[0] + spans[1]))
    }
}

fn segment(phi: f64) -> (f64, f64) {
    let sin = det::sin(phi);
    let cos = det::cos(phi);
    let numerator = if phi < 0.01 {
        let x = phi * phi;
        phi * x
            * (2.0 / 3.0
                + x * (-2.0 / 15.0 + x * (4.0 / 315.0 + x * (-2.0 / 2835.0 + x * 4.0 / 155925.0))))
    } else {
        phi - sin * cos
    };
    let q = numerator / (2.0 * sin * sin);
    (q, 1.0 - 2.0 * q * cos / sin)
}
