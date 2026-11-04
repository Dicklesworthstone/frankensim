//! Exponentially fitted normal-wall drag on the original staggered wave owner.
//!
//! In energy coordinates, x' = -C div(y), y' = -C grad(x) - R y + f.
//! R damps only the normal component. Freeze the gradient at midpoint x,
//! integrate each face exactly, and use its interval mean in the cell equation.
//! Eliminating faces leaves two scalar tridiagonals after a 2x2 rotation.
//! This retains the clamped-normal fourth-sound limit without a drag time-step
//! restriction. The authored Darcy closure still assumes a quasisteady
//! transverse velocity profile; this is not a resolved viscous boundary layer.

use super::*;
use fs_math::det;

#[derive(Clone, Copy, Debug)]
pub(super) struct ViscousWaveChange {
    pub body_work_j: f64,
    /// Signed wall impulse on the fluid and sum of absolute local impulses.
    pub wall_impulse_n_s: [f64; 2],
    pub relative_residual: f64,
}

#[derive(Clone, Copy, Debug)]
struct ViscousOperator {
    normal: [f64; 2],
    superfluid: [f64; 2],
    // Lambda = C*dt/dx; these columns map component velocities to cell fluxes.
    normal_flux: [f64; 2],
    superfluid_flux: [f64; 2],
    lambda: [f64; 2],
    // Orthonormal rows diagonalizing Lambda*phi2(-dt*R)*Lambda.
    cell_rows: [[f64; 2]; 2],
    normal_rotated_flux: [f64; 2],
    superfluid_rotated_flux: [f64; 2],
    r: [f64; 2],
    density: [f64; 2],
    density_root: [f64; 2],
    decay: f64,
    loss: f64,
    phi1: f64,
    phi2: f64,
    drag_forced: f64,
    heat_initial_root: f64,
    heat_shift: f64,
    heat_forced_root: f64,
}

impl ViscousOperator {
    fn new(channel: &TwoFluidChannel, dt: f64) -> Result<Self, ChannelError> {
        let density = [
            channel.material.normal_density_kg_m3(),
            channel.material.superfluid_density_kg_m3(),
        ];
        let density_root = density.map(f64::sqrt);
        let normal = channel.modes.face_to_modes(
            channel.material.normal_fraction() / density_root[0],
            1.0 / density_root[0],
        );
        let norm = normal[0].hypot(normal[1]);
        let normal = normal.map(|value| value / norm);
        // The modal kinetic rotation has positive determinant. This second
        // column therefore points along positive superfluid velocity.
        let superfluid = [normal[1], -normal[0]];
        let lambda = channel.modes.speed.map(|c| c * dt / channel.dx);
        let normal_flux = [lambda[0] * normal[0], lambda[1] * normal[1]];
        let superfluid_flux = [lambda[0] * superfluid[0], lambda[1] * superfluid[1]];
        let s = channel.normal_wall_rate_s_inv * dt;
        let decay = det::exp(-s);
        let loss = -det::expm1(-s);
        let (phi1, phi2, drag_forced, heat_shift, heat_forced) = if s < 0.125 {
            // phi2(-s) = 1/2 - s/6 + s^2/24 - ...
            // Evaluate it before phi1 so the small drag impulse never comes
            // from subtracting two nearly equal unit coefficients.
            let mut term = 0.5;
            let mut phi2 = term;
            for denominator in 3..=13 {
                term *= -s / f64::from(denominator);
                phi2 += term;
            }
            let drag_forced = s * phi2;
            let square = s * s;
            let shift = 0.5
                + square
                    * (-1.0 / 24.0
                        + square
                            * (1.0 / 240.0
                                + square
                                    * (-17.0 / 40320.0
                                        + square
                                            * (31.0 / 725760.0 - square * 691.0 / 159667200.0))));
            let heat_forced = s
                * (1.0 / 12.0
                    + square
                        * (-1.0 / 120.0
                            + square
                                * (17.0 / 20160.0
                                    + square
                                        * (-31.0 / 362880.0
                                            + square
                                                * (691.0 / 79833600.0
                                                    - square * 5461.0 / 6227020800.0)))));
            (1.0 - drag_forced, phi2, drag_forced, shift, heat_forced)
        } else {
            // The reciprocal ordering also handles an overflowed rate*dt:
            // its exponentially small transient is already unrepresentable.
            let inverse = if s.is_finite() {
                1.0 / s
            } else {
                (1.0 / channel.normal_wall_rate_s_inv) / dt
            };
            let phi1 = loss * inverse;
            let drag_forced = 1.0 - phi1;
            let shift = det::tanh(0.5 * s) * inverse;
            (
                phi1,
                drag_forced * inverse,
                drag_forced,
                shift,
                (1.0 - 2.0 * shift) * inverse,
            )
        };
        let heat_initial = -0.5 * det::expm1(-2.0 * s);

        // This dimensionless matrix is bounded by the admitted Courant step;
        // physical sound speeds need not be squared to find its rotation.
        let diagonal = [
            phi2 * normal_flux[0] * normal_flux[0] + 0.5 * superfluid_flux[0] * superfluid_flux[0],
            phi2 * normal_flux[1] * normal_flux[1] + 0.5 * superfluid_flux[1] * superfluid_flux[1],
        ];
        let off =
            phi2 * normal_flux[0] * normal_flux[1] + 0.5 * superfluid_flux[0] * superfluid_flux[1];
        let difference = diagonal[0] - diagonal[1];
        let (cos, sin) = if off == 0.0 {
            if difference >= 0.0 {
                (1.0, 0.0)
            } else {
                (0.0, 1.0)
            }
        } else {
            let gap = difference.hypot(2.0 * off);
            if difference >= 0.0 {
                let cos = (0.5 * (1.0 + difference / gap)).sqrt();
                (cos, off / (gap * cos))
            } else {
                let sin = (0.5 * (1.0 - difference / gap)).sqrt().copysign(off);
                (off / (gap * sin), sin)
            }
        };
        let rotation_norm = cos.hypot(sin);
        let row = [cos / rotation_norm, sin / rotation_norm];
        let cell_rows = [row, [-row[1], row[0]]];
        let normal_rotated_flux = cell_rows.map(|q| dot(q, normal_flux));
        let superfluid_rotated_flux = cell_rows.map(|q| dot(q, superfluid_flux));
        let r = [0, 1].map(|i| {
            0.5 * (phi2 * normal_rotated_flux[i] * normal_rotated_flux[i]
                + 0.5 * superfluid_rotated_flux[i] * superfluid_rotated_flux[i])
        });
        let cross = phi2 * normal_rotated_flux[0] * normal_rotated_flux[1]
            + 0.5 * superfluid_rotated_flux[0] * superfluid_rotated_flux[1];
        if !positive(norm)
            || !positive(rotation_norm)
            || ![
                decay,
                loss,
                phi1,
                phi2,
                drag_forced,
                heat_initial,
                heat_shift,
                heat_forced,
                r[0],
                r[1],
                cross,
            ]
            .into_iter()
            .all(f64::is_finite)
            || heat_initial < 0.0
            || heat_forced < 0.0
            || cross.abs() > 1024.0 * f64::EPSILON * (diagonal[0] + diagonal[1])
        {
            return Err(ChannelError::Invalid(
                "normal-wall exponential wave operator is not representable",
            ));
        }
        Ok(Self {
            normal,
            superfluid,
            normal_flux,
            superfluid_flux,
            lambda,
            cell_rows,
            normal_rotated_flux,
            superfluid_rotated_flux,
            r,
            density,
            density_root,
            decay,
            loss,
            phi1,
            phi2,
            drag_forced,
            heat_initial_root: heat_initial.sqrt(),
            heat_shift,
            heat_forced_root: heat_forced.sqrt(),
        })
    }

    fn components(&self, y: [f64; 2]) -> [f64; 2] {
        [dot(self.normal, y), dot(self.superfluid, y)]
    }

    // Direct component forces avoid manufacturing an undamped superfluid force
    // by subtracting large barycentric/relative normal-force contributions.
    fn force_increment(&self, force: BodyForce<'_>, face: usize, dt: f64) -> [f64; 2] {
        let value = match force {
            BodyForce::None => return [0.0; 2],
            BodyForce::Acceleration(a) => {
                return self.density_root.map(|root| root * (a * dt));
            }
            BodyForce::Uniform(value) => value,
            BodyForce::Faces(values) => values[face],
        };
        [
            self.density_root[0] * (value.normal_n_m3 / self.density[0] * dt),
            self.density_root[1] * (value.superfluid_n_m3 / self.density[1] * dt),
        ]
    }

    fn base_flux(
        &self,
        pending: &[Field; 2],
        mode: usize,
        face: usize,
        dt: f64,
        force: BodyForce<'_>,
    ) -> f64 {
        let y = [pending[0].faces[face], pending[1].faces[face]];
        if face == 0 || face + 1 == pending[0].faces.len() {
            return dot(
                self.cell_rows[mode],
                [self.lambda[0] * y[0], self.lambda[1] * y[1]],
            );
        }
        let old = self.components(y);
        let body = self.force_increment(force, face, dt);
        self.normal_rotated_flux[mode] * (self.phi1 * old[0] + self.phi2 * body[0])
            + self.superfluid_rotated_flux[mode] * (old[1] + 0.5 * body[1])
    }
}

impl TwoFluidChannel {
    pub(super) fn wave_cell_rotation(&self, dt: f64) -> Result<[[f64; 2]; 2], ChannelError> {
        if self.normal_wall_rate_s_inv == 0.0 {
            Ok([[1.0, 0.0], [0.0, 1.0]])
        } else {
            Ok(ViscousOperator::new(self, dt)?.cell_rows)
        }
    }

    pub(super) fn viscous_wave_step(
        &mut self,
        dt: f64,
        force: BodyForce<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscousWaveChange, ChannelError> {
        let operator = ViscousOperator::new(self, dt)?;
        let n = self.cells();
        let mut residual: f64 = 0.0;
        for mode in 0..2 {
            let q = operator.cell_rows[mode];
            let old =
                |i: usize| q[0] * self.accepted[0].cells[i] + q[1] * self.accepted[1].cells[i];
            let r = operator.r[mode];
            let work = &mut self.workspace[mode];
            super::step::factor_midpoint(work, r, dt, cancelled)?;
            let mut left_flux = operator.base_flux(&self.pending, mode, 0, dt, force);
            for i in 0..n {
                if i % 128 == 0 {
                    poll(cancelled)?;
                }
                let center = old(i);
                let mut laplacian = 0.0;
                if i > 0 {
                    laplacian += old(i - 1) - center;
                }
                if i + 1 < n {
                    laplacian += old(i + 1) - center;
                }
                let right_flux = operator.base_flux(&self.pending, mode, i + 1, dt, force);
                work.rhs[i] = center + r * laplacian - (right_flux - left_flux);
                left_flux = right_flux;
            }
            work.factor.solve_factored_into(
                &work.rhs,
                &mut self.pending[mode].cells,
                &mut *cancelled,
            )?;
            let mut error: f64 = 0.0;
            let mut scale: f64 = 0.0;
            for i in 0..n {
                if i % 128 == 0 {
                    poll(cancelled)?;
                }
                let mut applied = work.diagonal[i] * self.pending[mode].cells[i];
                let mut row_scale = applied.abs() + work.rhs[i].abs();
                if i > 0 {
                    applied -= r * self.pending[mode].cells[i - 1];
                    row_scale += (r * self.pending[mode].cells[i - 1]).abs();
                }
                if i + 1 < n {
                    applied -= r * self.pending[mode].cells[i + 1];
                    row_scale += (r * self.pending[mode].cells[i + 1]).abs();
                }
                error = error.max((applied - work.rhs[i]).abs());
                scale = scale.max(row_scale);
            }
            let relative = if scale == 0.0 { error } else { error / scale };
            if !relative.is_finite() || relative > 4096.0 * f64::EPSILON {
                return Err(ChannelError::Invalid(
                    "normal-wall wave solve has an unresolved residual",
                ));
            }
            residual = residual.max(relative);
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let rotated = [self.pending[0].cells[i], self.pending[1].cells[i]];
            self.pending[0].cells[i] =
                operator.cell_rows[0][0] * rotated[0] + operator.cell_rows[1][0] * rotated[1];
            self.pending[1].cells[i] =
                operator.cell_rows[0][1] * rotated[0] + operator.cell_rows[1][1] * rotated[1];
        }

        let mut body_work = 0.0;
        let mut work_correction = 0.0;
        let mut wall_impulse = 0.0;
        let mut impulse_correction = 0.0;
        let mut impulse_scale = 0.0;
        for i in 1..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let old = operator.components([self.pending[0].faces[i], self.pending[1].faces[i]]);
            let gradient = [0, 1].map(|m| {
                0.5 * ((self.accepted[m].cells[i] - self.accepted[m].cells[i - 1])
                    + (self.pending[m].cells[i] - self.pending[m].cells[i - 1]))
            });
            let body = operator.force_increment(force, i, dt);
            let increment = [
                body[0] - dot(operator.normal_flux, gradient),
                body[1] - dot(operator.superfluid_flux, gradient),
            ];
            let next = [
                operator.decay * old[0] + operator.phi1 * increment[0],
                old[1] + increment[1],
            ];
            let mean = [
                operator.phi1 * old[0] + operator.phi2 * increment[0],
                old[1] + 0.5 * increment[1],
            ];
            let y = [0, 1].map(|m| operator.normal[m] * next[0] + operator.superfluid[m] * next[1]);
            // For z'=-rate*z+g, s=rate*dt, w=dt*g, t=tanh(s/2):
            // integral(rate*z^2 dt) =
            // (1-exp(-2s))/2 * (z0+t*w/s)^2 + (s-2t)*w^2/s^2.
            // This positive two-square form avoids subtracting kinetic stores
            // and remains accurate both for tiny drag and forced stiff flow.
            let first = operator.heat_initial_root * (old[0] + operator.heat_shift * increment[0]);
            let second = operator.heat_forced_root * increment[0];
            let heat = self.cell_volume * (first * first + second * second);
            let work = self.cell_volume * dot(body, mean);
            let impulse = -self.cell_volume
                * operator.density_root[0]
                * (operator.loss * old[0] + operator.drag_forced * increment[0]);
            if ![y[0], y[1], heat, work, impulse]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(ChannelError::Domain {
                    field: "normal-wall wave/source update",
                    index: i,
                });
            }
            self.pending[0].faces[i] = y[0];
            self.pending[1].faces[i] = y[1];
            self.pending_viscous_heat_j[i - 1] += 0.5 * heat;
            self.pending_viscous_heat_j[i] += 0.5 * heat;
            accumulate(&mut body_work, &mut work_correction, work);
            accumulate(&mut wall_impulse, &mut impulse_correction, impulse);
            impulse_scale += impulse.abs();
        }
        poll(cancelled)?;
        if ![body_work, wall_impulse, impulse_scale]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(ChannelError::Invalid(
                "normal-wall work or impulse is not representable",
            ));
        }
        Ok(ViscousWaveChange {
            body_work_j: body_work,
            wall_impulse_n_s: [wall_impulse, impulse_scale],
            relative_residual: residual,
        })
    }
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
