//! Exact linear exchange while both attitudes are frozen. The nonlinear free
//! rotations remain in the existing pair splitting, with their own defects.

use super::super::exponential::Response;
use super::{
    Dd, HeldTorque, LieStepError, Mat3, PairStage, So3, Vector, ViscousPairStepper, apply,
    cholesky, compensated, dd, dot, drive_at, energy_change, invalid, scalar, slip, subtract,
    to_space, velocity,
};
use fs_la::eigen::jacobi_eigh;

#[derive(Debug, Clone)]
pub(super) struct PreparedExchange {
    relative_modes: Mat3,
    dual_modes: Mat3,
    response: [Response; 3],
    phi3: [f64; 3],
}

impl PreparedExchange {
    pub(super) fn prepare(
        owner: &ViscousPairStepper,
        rotations: [So3; 2],
    ) -> Result<Self, LieStepError> {
        // I_ref K = L L^T, y = Q^T L^-1 r. Over normalized stage time,
        // y' = g - Lambda y with Lambda = h/I_ref Q^T L^T D L Q.
        // This works with semidefinite D and never inverts the damping tensor.
        let mobility = owner.relative_mobility(rotations)?;
        let factor =
            cholesky(&mobility, 3).map_err(|_| invalid("resolved exponential pair mobility"))?;
        let lower = Mat3 {
            m: core::array::from_fn(|i| factor.l(i / 3, i % 3)),
        };
        let mut matrix = Mat3::zero();
        for i in 0..3 {
            for j in i..3 {
                let physical = dot(
                    column(lower, i),
                    apply(owner.damping_first_body, column(lower, j))?,
                );
                let scaled =
                    compensated::divide(dd(owner.half_s) * physical, dd(owner.reference_moment))?;
                let entry = scalar(scaled.to_f64(), "exponential pair relaxation tensor")?;
                if physical.to_f64() != 0.0 && entry == 0.0 {
                    return Err(invalid("resolved exponential pair relaxation tensor"));
                }
                matrix.m[3 * i + j] = entry;
                matrix.m[3 * j + i] = entry;
            }
        }
        let scale = matrix.max_abs();
        if scale == 0.0 {
            return Err(invalid("resolved exponential pair damping"));
        }
        let scaled = Mat3 {
            m: matrix.m.map(|v| v / scale),
        };
        let (_, vectors) = jacobi_eigh(&scaled.m, 3);
        let modes = Mat3 {
            m: core::array::from_fn(|i| vectors[i]),
        };
        let mut response = [Response::default(); 3];
        let mut phi3 = [0.0; 3];
        let mut relative_modes = Mat3::zero();
        let mut dual_modes = Mat3::zero();
        for i in 0..3 {
            let mode = column(modes, i);
            let applied = apply(scaled, mode)?;
            let eigenvalue = scalar(dot(mode, applied).to_f64(), "exponential pair eigenvalue")?;
            if eigenvalue < 0.0 {
                return Err(invalid("nonnegative exponential pair spectrum"));
            }
            for j in 0..3 {
                if (applied[j] - dd(eigenvalue) * mode[j]).to_f64().abs() > 4096.0 * f64::EPSILON
                    || (dot(mode, column(modes, j)).to_f64() - if i == j { 1.0 } else { 0.0 }).abs()
                        > 4096.0 * f64::EPSILON
                {
                    return Err(invalid("resolved orthogonal exponential pair spectrum"));
                }
            }
            let number = scalar(scale * eigenvalue, "exponential pair relaxation number")?;
            if eigenvalue != 0.0 && number == 0.0 {
                return Err(invalid("resolved exponential pair relaxation number"));
            }
            response[i] = Response::prepare(number, 1.0)?;
            phi3[i] = third_response(number, response[i].integrated_response);
            let relative = apply(lower, mode)?;
            let mut dual = relative.map(Dd::to_f64);
            factor.solve(&mut dual);
            for j in 0..3 {
                relative_modes.m[3 * j + i] =
                    scalar(relative[j].to_f64(), "exponential pair mode")?;
                dual_modes.m[3 * j + i] = scalar(dual[j], "exponential pair dual mode")?;
            }
        }
        // These are reciprocal coordinate transforms, not independent modal
        // approximations of the two original bodies. Refuse unresolved mobility.
        for i in 0..3 {
            for j in 0..3 {
                let product = dot(column(dual_modes, i), column(relative_modes, j)).to_f64();
                if (product - if i == j { 1.0 } else { 0.0 }).abs() > 8192.0 * f64::EPSILON {
                    return Err(invalid("resolved exponential pair duality"));
                }
            }
        }
        Ok(Self {
            relative_modes,
            dual_modes,
            response,
            phi3,
        })
    }
}

impl ViscousPairStepper {
    pub(super) fn exponential_stage(
        &self,
        rotations: [So3; 2],
        momentum: [Vector; 2],
        drives: [HeldTorque; 2],
    ) -> Result<PairStage, LieStepError> {
        let fresh;
        let prepared = if let Some(cached) = &self.cached_exponential {
            cached
        } else {
            fresh = PreparedExchange::prepare(self, rotations)?;
            &fresh
        };
        let before = [
            velocity(self.inertia[0], momentum[0])?,
            velocity(self.inertia[1], momentum[1])?,
        ];
        let mut drive_impulse = [[Dd::ZERO; 3]; 2];
        let mut drive_change = [[Dd::ZERO; 3]; 2];
        for body in 0..2 {
            let drive = drive_at(rotations[body], drives[body])?;
            for i in 0..3 {
                drive_impulse[body][i] = compensated::checked(
                    dd(self.half_s) * drive[i],
                    "exponential pair drive impulse",
                )?;
            }
            drive_change[body] = velocity(self.inertia[body], drive_impulse[body])?;
        }
        let relative = slip(rotations, before)?;
        let relative_drive = slip(rotations, drive_change)?;
        let initial = apply(prepared.dual_modes.transpose(), relative)?;
        let source = apply(prepared.dual_modes.transpose(), relative_drive)?;
        let mut change = [Dd::ZERO; 3];
        let mut integral = [Dd::ZERO; 3];
        let mut weighted = [Dd::ZERO; 3];
        let mut heat = Dd::ZERO;
        let root_reference = dd(self.reference_moment.sqrt());
        for i in 0..3 {
            let r = prepared.response[i];
            change[i] = dd(r.decay_change) * initial[i] + dd(r.response) * source[i];
            integral[i] = dd(r.response) * initial[i] + dd(r.integrated_response) * source[i];
            weighted[i] = dd(r.integrated_response) * initial[i] + dd(prepared.phi3[i]) * source[i];
            let mean = dd(r.heat_mean[0]) * initial[i] + dd(r.heat_mean[1]) * source[i];
            let deviation =
                dd(r.heat_deviation[0]) * initial[i] + dd(r.heat_deviation[1]) * source[i];
            // Weight before squaring: raw rate squares can overflow or become
            // subnormal even while the physical inertia-weighted heat is finite.
            let weighted_mean = root_reference * mean;
            let weighted_deviation = root_reference * deviation;
            heat = heat + weighted_mean * weighted_mean + weighted_deviation * weighted_deviation;
        }
        if scalar(heat.to_f64(), "exponential pair heat")? < 0.0 {
            return Err(invalid("nonnegative exponential pair heat"));
        }
        // j is the coupling impulse on body 0. weighted_j = integral(j(t) dt)/h
        // supplies the exact drive/transfer work, including drives on both bodies.
        let impulse = apply(
            self.damping_first_body,
            apply(prepared.relative_modes, integral)?,
        )?
        .map(|v| dd(-self.half_s) * v);
        let weighted_impulse = apply(
            self.damping_first_body,
            apply(prepared.relative_modes, weighted)?,
        )?
        .map(|v| dd(-self.half_s) * v);
        let coupling_impulse_space = to_space(rotations[0], impulse)?;
        let weighted_space = to_space(rotations[0], weighted_impulse)?;
        let coupling = [
            impulse,
            compensated::spatial_vector_to_body(rotations[1], coupling_impulse_space)?.map(|v| -v),
        ];
        let weighted = [
            weighted_impulse,
            compensated::spatial_vector_to_body(rotations[1], weighted_space)?.map(|v| -v),
        ];
        let mut after = momentum;
        let mut actual_change = [[Dd::ZERO; 3]; 2];
        let mut coupling_change = [[Dd::ZERO; 3]; 2];
        let mut drive_work = [Dd::ZERO; 2];
        let mut coupling_work = [Dd::ZERO; 2];
        let mut work_scale = heat.to_f64().abs();
        for body in 0..2 {
            for i in 0..3 {
                after[body][i] = compensated::checked(
                    momentum[body][i] + drive_impulse[body][i] + coupling[body][i],
                    "exponential pair kicked momentum",
                )?;
            }
            let actual = subtract(after[body], momentum[body]);
            actual_change[body] = velocity(self.inertia[body], actual)?;
            coupling_change[body] = velocity(self.inertia[body], coupling[body])?;
            let weighted_change = velocity(self.inertia[body], weighted[body])?;
            let average = core::array::from_fn(|i| {
                before[body][i] + dd(0.5) * drive_change[body][i] + weighted_change[i]
            });
            let transfer_rate =
                core::array::from_fn(|i| before[body][i] + dd(0.5) * coupling_change[body][i]);
            let time_moment = subtract(coupling[body], weighted[body]);
            drive_work[body] = dot(drive_impulse[body], average);
            // Integrate omega(t).j'(t) using p(t)=p0+t*drive+j(t).
            // The quadratic impulse primitive is exact for constant full inertia.
            coupling_work[body] =
                dot(coupling[body], transfer_rate) + dot(drive_change[body], time_moment);
            let kinetic = energy_change(self.inertia[body], momentum[body], after[body])?;
            let mut body_scale = kinetic
                .to_f64()
                .abs()
                .max(drive_work[body].to_f64().abs())
                .max(coupling_work[body].to_f64().abs());
            for i in 0..3 {
                for work in [
                    drive_impulse[body][i] * average[i],
                    coupling[body][i] * transfer_rate[i],
                    drive_change[body][i] * time_moment[i],
                ] {
                    body_scale = body_scale
                        .max(scalar(work.to_f64(), "exponential pair component work")?.abs());
                }
            }
            let residual = scalar(
                (kinetic - drive_work[body] - coupling_work[body]).to_f64(),
                "exponential pair body work",
            )?;
            if residual.abs() > 32768.0 * f64::EPSILON * body_scale {
                return Err(invalid("resolved exponential pair body work"));
            }
            work_scale = work_scale.max(body_scale);
        }
        let actual_relative = slip(rotations, actual_change)?;
        let coupling_relative = slip(rotations, coupling_change)?;
        let expected_relative = apply(prepared.relative_modes, change)?;
        let rate_scale = actual_relative
            .iter()
            .chain(relative_drive.iter())
            .chain(coupling_relative.iter())
            .map(|v| v.to_f64().abs())
            .fold(0.0_f64, f64::max);
        for i in 0..3 {
            let residual = scalar(
                (actual_relative[i] - expected_relative[i]).to_f64(),
                "exponential pair slip response",
            )?;
            if residual.abs() > 32768.0 * f64::EPSILON * rate_scale {
                return Err(invalid("resolved exponential pair slip response"));
            }
        }
        let residual = scalar(
            (coupling_work[0] + coupling_work[1] + heat).to_f64(),
            "exponential pair reciprocal work",
        )?;
        if residual.abs() > 32768.0 * f64::EPSILON * work_scale {
            return Err(invalid("resolved exponential pair reciprocal work"));
        }
        Ok(PairStage {
            momentum: after,
            drive_impulse,
            coupling_impulse_space,
            drive_work,
            coupling_work,
            heat,
        })
    }
}

fn column(matrix: Mat3, index: usize) -> Vector {
    [
        dd(matrix.m[index]),
        dd(matrix.m[3 + index]),
        dd(matrix.m[6 + index]),
    ]
}

// phi_3(-x) integrates the already-integrated slip; the series retains the
// 1/6 limit without cancellation when physical damping is weak.
fn third_response(x: f64, phi2: f64) -> f64 {
    if x <= 0.5 {
        let mut term = 1.0 / 6.0;
        let mut value = term;
        for k in 0..24 {
            term *= -x / (k + 4) as f64;
            value += term;
        }
        value
    } else {
        (0.5 - phi2) / x
    }
}
