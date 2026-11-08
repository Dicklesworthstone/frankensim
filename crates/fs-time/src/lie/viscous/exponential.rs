//! Prepared exact frozen-attitude linear drag; the outer free rotation keeps
//! its own numerical defect. No physical heat is inferred from an energy gap.

use super::{
    InertiaTensor3, LieStepError, Mat3, Stage, Vector, add, apply, compensated, components, dd,
    dot, energy_change, invalid, relative, scalar, subtract,
};
use fs_la::eigen::jacobi_eigh;
use fs_math::{dd::Dd, det};

#[derive(Debug, Clone)]
pub(super) struct ExponentialRelaxation {
    velocity_modes: Mat3,
    momentum_modes: Mat3,
    response: [Response; 3],
    half_s: f64,
    rate_bound: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub(in crate::lie) struct Response {
    pub(in crate::lie) decay_change: f64,
    pub(in crate::lie) response: f64,
    pub(in crate::lie) integrated_response: f64,
    // Two squares integrate each mode's heat, avoiding cancellation between
    // initial/forcing terms, including exact equilibrium and weak relaxation.
    pub(in crate::lie) heat_mean: [f64; 2],
    pub(in crate::lie) heat_deviation: [f64; 2],
}

impl ExponentialRelaxation {
    pub(super) fn prepare(
        inertia: InertiaTensor3,
        damping: Mat3,
        half_s: f64,
    ) -> Result<Self, LieStepError> {
        let axes = inertia.principal_to_body().matrix();
        let roots = components(inertia.principal_moments()).map(f64::sqrt);
        let mut velocity_principal = axes;
        let mut momentum_principal = axes;
        for i in 0..3 {
            for (j, root) in roots.iter().enumerate() {
                velocity_principal.m[3 * i + j] /= root;
                momentum_principal.m[3 * i + j] *= root;
            }
        }
        let mut normalized = [0.0; 9];
        for i in 0..3 {
            let column_i = column(velocity_principal, i);
            for j in i..3 {
                let entry = scalar(
                    dot(column_i, apply(damping, column(velocity_principal, j))?).to_f64(),
                    "exponential viscous normalized tensor",
                )?;
                normalized[3 * i + j] = entry;
                normalized[3 * j + i] = entry;
            }
        }
        let scale = normalized.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(invalid("resolved exponential viscous mobility"));
        }
        let scaled = normalized.map(|x| x / scale);
        let (_, vectors) = jacobi_eigh(&scaled, 3);
        let basis = Mat3 {
            m: core::array::from_fn(|i| vectors[i]),
        };
        let mut rates = [0.0; 3];
        for (mode, rate) in rates.iter_mut().enumerate() {
            let v = column(basis, mode);
            let av = apply(Mat3 { m: scaled }, v)?;
            let eigenvalue = scalar(dot(v, av).to_f64(), "exponential viscous eigenvalue")?;
            if eigenvalue < 0.0 {
                return Err(invalid("nonnegative exponential viscous spectrum"));
            }
            for i in 0..3 {
                if (av[i] - dd(eigenvalue) * v[i]).to_f64().abs() > 4096.0 * f64::EPSILON {
                    return Err(invalid("resolved exponential viscous eigenpair"));
                }
            }
            for j in 0..3 {
                let expected = if mode == j { 1.0 } else { 0.0 };
                if (dot(v, column(basis, j)).to_f64() - expected).abs() > 4096.0 * f64::EPSILON {
                    return Err(invalid("orthogonal exponential viscous modes"));
                }
            }
            *rate = scalar(scale * eigenvalue, "exponential viscous relaxation rate")?;
        }
        let velocity_modes = product(velocity_principal, basis)?;
        let momentum_modes = product(momentum_principal, basis)?;
        let mut response = [Response::default(); 3];
        for i in 0..3 {
            response[i] = Response::prepare(rates[i], half_s)?;
        }
        let rate_bound = scalar(
            rates.iter().sum::<f64>() * (1.0 + 256.0 * f64::EPSILON),
            "exponential viscous relaxation bound",
        )?;
        Ok(Self {
            velocity_modes,
            momentum_modes,
            response,
            half_s,
            rate_bound,
        })
    }

    pub(super) const fn rate_bound(&self) -> f64 {
        self.rate_bound
    }

    pub(super) fn stage(
        &self,
        inertia: InertiaTensor3,
        damping: Mat3,
        momentum: Vector,
        drive: Vector,
        bath: Option<Vector>,
    ) -> Result<Stage, LieStepError> {
        // y = Q^T I^-1/2 L and beta = Q^T I^1/2 Omega; relative
        // modes satisfy x' = g - lambda x during this fixed-attitude stage.
        let before = apply(self.velocity_modes.transpose(), momentum)?;
        let bath_modes = bath
            .map(|spin| apply(self.momentum_modes.transpose(), spin))
            .transpose()?;
        let slip = relative(before, bath_modes);
        let source = apply(self.velocity_modes.transpose(), drive)?;
        let mut increment = [Dd::ZERO; 3];
        let mut integrated_slip = [Dd::ZERO; 3];
        let mut heat = Dd::ZERO;
        for i in 0..3 {
            let r = self.response[i];
            increment[i] = dd(r.decay_change) * slip[i] + dd(r.response) * source[i];
            integrated_slip[i] = dd(r.response) * slip[i] + dd(r.integrated_response) * source[i];
            let mean = dd(r.heat_mean[0]) * slip[i] + dd(r.heat_mean[1]) * source[i];
            let deviation = dd(r.heat_deviation[0]) * slip[i] + dd(r.heat_deviation[1]) * source[i];
            heat = heat + mean * mean + deviation * deviation;
        }
        // Apply only the increment: reconstructing the whole momentum would
        // destroy resolved small losses beside a much larger undamped spin.
        let actual = apply(self.momentum_modes, increment)?;
        let after = add(momentum, actual);
        for component in after {
            compensated::checked(component, "exponential viscous momentum")?;
        }
        let integrated_slip = apply(self.velocity_modes, integrated_slip)?;
        let integrated_spin = bath.map_or(integrated_slip, |spin| {
            add(integrated_slip, spin.map(|x| dd(self.half_s) * x))
        });
        let drive_impulse = drive.map(|x| dd(self.half_s) * x);
        let drag_impulse = apply(damping, integrated_slip)?.map(|x| dd(-1.0) * x);
        let drive_work = dot(drive, integrated_spin);
        let bath_work = bath.map_or(Dd::ZERO, |spin| dot(drag_impulse, spin));
        scalar(heat.to_f64(), "exponential viscous heat")?;
        let change = energy_change(inertia, momentum, after)?;
        let work_residual = scalar(
            (change - drive_work - bath_work + heat).to_f64(),
            "exponential viscous work residual",
        )?;
        let mut work_scale = change
            .to_f64()
            .abs()
            .max(drive_work.to_f64().abs())
            .max(bath_work.to_f64().abs())
            .max(heat.to_f64().abs());
        let actual = subtract(after, momentum);
        for i in 0..3 {
            let residual = scalar(
                (actual[i] - drive_impulse[i] - drag_impulse[i]).to_f64(),
                "exponential viscous impulse residual",
            )?;
            let impulse_scale = actual[i]
                .to_f64()
                .abs()
                .max(drive_impulse[i].to_f64().abs())
                .max(drag_impulse[i].to_f64().abs());
            if residual.abs() > 16384.0 * f64::EPSILON * impulse_scale {
                return Err(invalid("resolved exponential viscous impulse"));
            }
            work_scale = work_scale.max((drive[i] * integrated_spin[i]).to_f64().abs());
            if let Some(spin) = bath {
                work_scale = work_scale.max((drag_impulse[i] * spin[i]).to_f64().abs());
            }
        }
        if work_residual.abs() > 32768.0 * f64::EPSILON * work_scale {
            return Err(invalid("resolved exponential viscous work"));
        }
        Ok(Stage {
            momentum: after,
            drive_impulse,
            drag_impulse,
            drive_work,
            bath_work,
            heat,
        })
    }
}

impl Response {
    pub(in crate::lie) fn prepare(rate: f64, h: f64) -> Result<Self, LieStepError> {
        if rate == 0.0 {
            return Ok(Self {
                response: h,
                integrated_response: 0.5 * h * h,
                ..Self::default()
            });
        }
        let x = scalar(h * rate, "exponential viscous relaxation number")?;
        let decay_change = det::expm1(-x);
        let f = if x == 0.0 { 1.0 } else { -decay_change / x };
        // v = phi_2(-x). V = variance(exp(-x*t), t in [0,1])/x^2.
        // Their series avoid subtracting nearly equal numbers in weak drag.
        let (v, root_variance, small) = if x <= 0.5 {
            let mut v = 0.5;
            let mut term_v = v;
            let mut bracket = 1.0 / 12.0;
            let mut term_b = bracket;
            for k in 0..24 {
                term_v *= -x / (k + 3) as f64;
                term_b *= -x * (k + 2) as f64 / ((k + 1) * (k + 4)) as f64;
                v += term_v;
                bracket += term_b;
            }
            (v, (f * bracket).sqrt(), true)
        } else {
            let v = (1.0 - f) / x;
            let twice = if x > 350.0 {
                0.5 / x
            } else {
                -det::expm1(-2.0 * x) / (2.0 * x)
            };
            (v, (twice - f * f).sqrt(), false)
        };
        let root_x = x.sqrt();
        let heat_deviation = if small {
            [-x * root_x * root_variance, h * (root_x * root_variance)]
        } else {
            [-root_x * root_variance, h * (root_variance / root_x)]
        };
        let response = Self {
            decay_change,
            response: h * f,
            integrated_response: h * (h * v),
            heat_mean: [root_x * f, root_x * (h * v)],
            heat_deviation,
        };
        for coefficient in [
            response.decay_change,
            response.response,
            response.integrated_response,
            response.heat_mean[0],
            response.heat_mean[1],
            response.heat_deviation[0],
            response.heat_deviation[1],
        ] {
            scalar(coefficient, "resolved exponential viscous coefficient")?;
        }
        Ok(response)
    }
}

fn column(matrix: Mat3, index: usize) -> Vector {
    [
        dd(matrix.m[index]),
        dd(matrix.m[3 + index]),
        dd(matrix.m[6 + index]),
    ]
}

fn product(left: Mat3, right: Mat3) -> Result<Mat3, LieStepError> {
    let mut result = Mat3::zero();
    for j in 0..3 {
        let value = apply(left, column(right, j))?;
        for (i, entry) in value.iter().enumerate() {
            result.m[3 * i + j] = scalar(entry.to_f64(), "exponential viscous basis")?;
        }
    }
    Ok(result)
}
