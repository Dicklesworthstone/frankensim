//! Same-time distributed contact/flow solve on all retained material modes.
//! The existing exact or midpoint structural response remains the only linear
//! owner. Trials change neither modal state nor hereditary material history.
use super::closure::{ModalPlateClosure, displacements};
use super::{
    AcousticRealizeError, ApertureDrive, MAX_MODES, ModalAperture, invalid, nonlinear, uniform,
};
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;

#[derive(Clone, Copy)]
struct SpatialTrial {
    velocity: [f64; MAX_MODES],
    dp: f64,
    residual: [f64; MAX_MODES],
    scale: [f64; MAX_MODES],
    flow_residual: f64,
    flow_scale: f64,
}
impl SpatialTrial {
    fn norm(&self, n: usize) -> f64 {
        self.residual[..n].iter().fold(0.0, |sum, x| sum.hypot(*x))
    }
}
impl ModalAperture {
    /// Replace the effective uniform slit/contact with actual nodal geometry.
    /// Uses every original signed DKT mode, one triangle contact quadrature and
    /// the first blade's original physical slit edges. Independently meshed
    /// blades must cover the same simple planar footprint and slit in the same
    /// frame. Original P1 traces and their transpose preserve contact work;
    /// source and interface quadrature resolutions remain independent controls.
    /// Initial motion, full viscosity and exact material memory are preserved.
    /// Pressure and swept flow retain their original conjugate face projections.
    ///
    /// # Errors
    /// Missing/noncovering or ambiguous sources, exhausted interpolation budget,
    /// repeated/late attachment, invalid clearances or law, excessive local
    /// penetration, slope or complete stored energy.
    pub fn with_plate_closure(
        mut self,
        spec: PlateClosureSpec,
    ) -> Result<Self, AcousticRealizeError> {
        if self.samples != 0 || self.closure.is_some() {
            return Err(invalid(
                "modal spatial closure requires one attachment before stepping",
            ));
        }
        let closure = ModalPlateClosure::from_plates(&self.material_plates, spec)?;
        closure.validate_state(&self.state)?;
        if let Some(law) = self.viscous_slit {
            let q = displacements(&self.state);
            closure.viscous_flow(law, &q[..self.state.len()], self.spec.density_kg_m3, 0.0)?;
        }
        self.closure = Some(closure.into());
        self.validate_state(&self.state)?;
        let energy = self.stored_energy_j();
        if !energy.is_finite() || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid(
                "modal spatial contact exceeds the complete energy allowance",
            ));
        }
        Ok(self)
    }
    /// Actual source-bound geometry; absent means the original effective slit.
    pub fn plate_closure(&self) -> Option<&ModalPlateClosure> {
        self.closure.as_deref()
    }

    fn force_velocity(
        &self,
        force: &[f64],
        mobility: &[f64],
        out: &mut [f64],
    ) -> Result<(), AcousticRealizeError> {
        if let Some(exact) = &self.exact {
            exact.force_velocity(force, out);
        } else if let Some(c) = self.viscosity.as_ref().filter(|c| c.is_active()) {
            for (v, f) in out.iter_mut().zip(force) {
                *v = 0.5 * self.spec.time_step_s * f;
            }
            c.solve_free(out)?;
        } else {
            for ((v, f), b) in out.iter_mut().zip(force).zip(mobility) {
                *v = b * f;
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("modal spatial force response overflowed"));
        }
        Ok(())
    }

    pub(super) fn solve_profiled(
        &self,
        drive: ApertureDrive,
        free: &[f64],
        mobility: &[f64],
        areas: &[f64],
        closure: &ModalPlateClosure,
    ) -> Result<uniform::UniformTrial, AcousticRealizeError> {
        let n = free.len();
        let before = displacements(&self.state);
        let contact_step = closure.contact_step(&before[..n])?;
        let dt = self.spec.time_step_s;
        let z = self.spec.impedance_pa_s_m3;
        let mut pressure = [0.0; MAX_MODES];
        self.force_velocity(areas, mobility, &mut pressure[..n])?;
        let pressure_mobility: f64 = areas.iter().zip(&pressure).map(|(a, p)| a * p).sum();
        let denominator = 1.0 + z * pressure_mobility;
        if !pressure_mobility.is_finite()
            || pressure_mobility < 0.0
            || !denominator.is_finite()
            || denominator <= 0.0
        {
            return Err(invalid("modal spatial swept response is unresolved"));
        }
        let evaluations = std::cell::Cell::new(0);
        let evaluate = |velocity: [f64; MAX_MODES]| -> Result<SpatialTrial, AcousticRealizeError> {
            evaluations.set(evaluations.get() + 1);
            if evaluations.get() > self.spec.limits.max_evaluations {
                return Err(invalid(
                    "modal spatial nonlinear evaluation budget exhausted",
                ));
            }
            let mut after = before;
            let mut middle = before;
            for i in 0..n {
                after[i] += dt * velocity[i];
                middle[i] = f64::midpoint(before[i], after[i]);
            }
            let contact = contact_step.response(&after[..n], &velocity[..n])?;
            let mut base = [0.0; MAX_MODES];
            self.force_velocity(&contact.forces_n[..n], mobility, &mut base[..n])?;
            let mut moving = drive.body_flow_m3_s;
            for i in 0..n {
                base[i] += free[i];
                moving -= areas[i] * base[i];
            }
            let source =
                (drive.upstream_pressure_pa - 2.0 * drive.incoming_pressure_pa - z * moving)
                    / denominator;
            let dp = match self.viscous_slit {
                Some(law) => crate::bernoulli_aperture::viscous::profiled_pressure_drop(
                    z / denominator,
                    source,
                    |dp| closure.viscous_flow(law, &middle[..n], self.spec.density_kg_m3, dp),
                )?,
                None => crate::bernoulli_aperture::moving::fixed_opening_pressure_drop(
                    self.spec.width_m,
                    closure.open_area_m2(&middle[..n])? / self.spec.width_m,
                    self.spec.density_kg_m3,
                    z / denominator,
                    source,
                )?,
            };
            let mut residual = [0.0; MAX_MODES];
            let mut scale = [0.0; MAX_MODES];
            for i in 0..n {
                let predicted = base[i] - pressure[i] * dp;
                residual[i] = velocity[i] - predicted;
                scale[i] = velocity[i].abs() + predicted.abs() + free[i].abs();
            }
            // dp balanced the PREDICTED swept flow; the frame reports this
            // iterate, so accepting it leaves -sum a_i r_i in the flow balance.
            // ModalFlow certifies that against relative_tolerance times the
            // |jet|+|swept|+|body|+|wave| scale with no absolute term, which the
            // per-mode test (absolute velocity_tolerance_m_s) does not imply.
            // Hold half of that budget here, on the smaller scale without |wave|.
            let jet = match self.viscous_slit {
                Some(law) => closure
                    .viscous_flow(law, &middle[..n], self.spec.density_kg_m3, dp)?
                    .flow_m3_s,
                None => crate::bernoulli_aperture::moving::volume_flow(
                    self.spec.width_m,
                    closure.open_area_m2(&middle[..n])? / self.spec.width_m,
                    dp,
                    self.spec.density_kg_m3,
                ),
            };
            let mut flow_residual = 0.0;
            let mut swept = 0.0;
            for i in 0..n {
                flow_residual += areas[i] * residual[i];
                swept -= areas[i] * velocity[i];
            }
            let flow_scale = jet.abs() + swept.abs() + drive.body_flow_m3_s.abs();
            if !dp.is_finite()
                || !flow_residual.is_finite()
                || !flow_scale.is_finite()
                || residual[..n]
                    .iter()
                    .chain(&scale[..n])
                    .any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "modal spatial trial left its finite physical domain",
                ));
            }
            Ok(SpatialTrial {
                velocity,
                dp,
                residual,
                scale,
                flow_residual,
                flow_scale,
            })
        };
        let resolved = |t: &SpatialTrial| {
            t.flow_residual.abs() <= 0.5 * self.spec.limits.relative_tolerance * t.flow_scale
                && (0..n).all(|i| {
                    t.residual[i].abs()
                        <= self.spec.limits.velocity_tolerance_m_s
                            + self.spec.limits.relative_tolerance * t.scale[i]
                })
        };
        let mut velocity = [0.0; MAX_MODES];
        velocity[..n].copy_from_slice(free);
        let mut trial = evaluate(velocity)?;
        loop {
            if resolved(&trial) {
                return Ok(uniform::UniformTrial {
                    velocity: trial.velocity,
                    dp: trial.dp,
                    residual: trial.norm(n),
                    scale: trial.scale[..n].iter().sum(),
                    evaluations: evaluations.get(),
                });
            }
            // Every Jacobian probe and backtrack consumes the same explicit
            // nonlinear budget. All coordinates are solved together; there is
            // no lagged contact force or fallback to a mean-gap model.
            let mut jacobian = [0.0; MAX_MODES * MAX_MODES];
            for column in 0..n {
                let delta = 1e-6 * (trial.velocity[column].abs() + 0.01);
                let mut left = trial.velocity;
                let mut right = trial.velocity;
                left[column] -= delta;
                right[column] += delta;
                let width = right[column] - left[column];
                if !width.is_finite() || width <= 0.0 {
                    return Err(invalid("modal spatial Jacobian increment is unresolved"));
                }
                let a = evaluate(left)?;
                let b = evaluate(right)?;
                for row in 0..n {
                    jacobian[n * row + column] = (b.residual[row] - a.residual[row]) / width;
                }
            }
            let mut direction = [0.0; MAX_MODES];
            for row in 0..n {
                let scale = jacobian[n * row..n * (row + 1)]
                    .iter()
                    .fold(0.0_f64, |s, x| s.max(x.abs()));
                if !scale.is_finite() || scale == 0.0 {
                    return Err(invalid("modal spatial Jacobian has an unresolved row"));
                }
                for value in &mut jacobian[n * row..n * (row + 1)] {
                    *value /= scale;
                }
                direction[row] = -trial.residual[row] / scale;
            }
            let rhs = direction;
            let factor = fs_la::factor::lu(&jacobian[..n * n], n).map_err(nonlinear)?;
            factor.solve(&mut direction[..n]);
            for row in 0..n {
                let terms = jacobian[n * row..n * (row + 1)].iter().zip(&direction[..n]);
                let product: f64 = terms.clone().map(|(a, b)| a * b).sum();
                let scale = rhs[row].abs() + terms.map(|(a, b)| (a * b).abs()).sum::<f64>();
                if !product.is_finite()
                    || !scale.is_finite()
                    || (product - rhs[row]).abs()
                        > 512.0 * f64::EPSILON * n as f64 * scale.max(f64::MIN_POSITIVE)
                {
                    return Err(invalid(
                        "modal spatial Newton direction has an unresolved linear residual",
                    ));
                }
            }
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..16 {
                let mut velocity = trial.velocity;
                for i in 0..n {
                    velocity[i] += fraction * direction[i];
                }
                let next = evaluate(velocity)?;
                if resolved(&next) || next.norm(n) < trial.norm(n) {
                    accepted = Some(next);
                    break;
                }
                fraction *= 0.5;
            }
            trial = accepted.ok_or_else(|| {
                invalid("modal spatial solve failed to reduce all momentum residuals")
            })?;
        }
    }
}
