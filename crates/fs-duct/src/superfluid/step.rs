use super::*;
use fs_material::superfluid::{TwoFluidVelocity, VortexMutualFriction};

impl TwoFluidChannel {
    pub(super) fn stage_step(
        &mut self,
        dt: f64,
        boundary: BoundaryVelocities,
        friction: MutualFriction,
        force: BodyForce<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelStep, ChannelError> {
        self.stage_step_with_limits(dt, boundary, friction, force, true, cancelled)
    }

    // A nonpublishing affine pressure query may use an artificial reference
    // velocity. Only the final physical proposal is subject to amplitude limits;
    // both paths retain the same original operator, finite checks and balances.
    pub(super) fn stage_step_with_limits(
        &mut self,
        dt: f64,
        boundary: BoundaryVelocities,
        friction: MutualFriction,
        force: BodyForce<'_>,
        physical_limits: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelStep, ChannelError> {
        poll(cancelled)?;
        let next_time = self.time_s + dt;
        if !positive(dt)
            || !next_time.is_finite()
            || next_time <= self.time_s
            || dt > self.maximum_step_s() * (1.0 + 8.0 * f64::EPSILON)
        {
            return Err(ChannelError::Invalid(
                "step must advance a finite clock and satisfy the declared wave Courant <= 1 policy",
            ));
        }
        let friction = match friction {
            MutualFriction::None => friction,
            MutualFriction::Linear { gamma_kg_m3_s: k }
            | MutualFriction::GorterMellink { a_m_s_kg: k } => {
                if !k.is_finite() || k < 0.0 {
                    return Err(ChannelError::Invalid(
                        "mutual-friction coefficient must be nonnegative and finite",
                    ));
                }
                friction
            }
            MutualFriction::IsotropicVortexTangle {
                circulation_quantum_m2_s,
                line_density_per_m2,
                dissipative_b,
            } => {
                let vortices = VortexMutualFriction::try_new(
                    circulation_quantum_m2_s,
                    line_density_per_m2,
                    dissipative_b,
                    0.0,
                )?;
                MutualFriction::Linear {
                    gamma_kg_m3_s: self.material.isotropic_vortex_drag_kg_m3_s(vortices)?,
                }
            }
        };
        let endpoints = [boundary.left, boundary.right];
        let mut port_modes = [[0.0; 2]; 2];
        for (i, end) in endpoints.into_iter().enumerate() {
            if [end.normal_m_s, end.superfluid_m_s].into_iter().any(|v| {
                !v.is_finite()
                    || (physical_limits
                        && (v.abs() > self.limits.max_component_speed_m_s
                            || v.abs() > 0.1 * self.modes.speed[1]))
            }) {
                return Err(ChannelError::Domain {
                    field: "boundary velocity",
                    index: i,
                });
            }
            let v = self.material.normal_fraction() * end.normal_m_s
                + self.material.superfluid_fraction() * end.superfluid_m_s;
            port_modes[i] = self
                .modes
                .face_to_modes(v, end.normal_m_s - end.superfluid_m_s);
        }
        let n = self.cells();
        force.validate(n + 1, &self.material, cancelled)?;
        let before = energy(&self.accepted, self.cell_volume, cancelled)?;
        for i in 0..=n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            for mode in 0..2 {
                if i < n {
                    self.pending[mode].cells[i] = self.accepted[mode].cells[i];
                }
                self.pending[mode].faces[i] = self.accepted[mode].faces[i];
            }
            if i < n {
                self.pending_heat_j[i] = 0.0;
            }
        }
        for mode in 0..2 {
            self.pending[mode].faces[0] = port_modes[0][mode];
            self.pending[mode].faces[n] = port_modes[1][mode];
        }
        self.friction_half_step(friction, dt / 2.0, cancelled)?;
        // Keep the exact force kicks adjacent to the wave update. Friction
        // outside this affine midpoint block preserves static force balance.
        let mut source_work = self.force_half_step(force, dt / 2.0, cancelled)?;
        let mut work = 0.0;
        let mut residual: f64 = 0.0;
        for mode in 0..2 {
            residual = residual.max(midpoint(
                &self.accepted[mode].cells,
                &mut self.pending[mode],
                &mut self.workspace[mode],
                self.modes.speed[mode] * dt / self.dx,
                dt,
                cancelled,
            )?);
            let old = &self.accepted[mode].cells;
            let new = &self.pending[mode].cells;
            work += dt
                * self.area_m2
                * self.modes.speed[mode]
                * (0.5 * (old[0] + new[0]) * port_modes[0][mode]
                    - 0.5 * (old[n - 1] + new[n - 1]) * port_modes[1][mode]);
        }
        source_work += self.force_half_step(force, dt / 2.0, cancelled)?;
        if force != BodyForce::None {
            self.validate_pending_with_limits(physical_limits, cancelled)?;
        }
        self.friction_half_step(friction, dt / 2.0, cancelled)?;
        self.validate_pending_with_limits(physical_limits, cancelled)?;
        let after = energy(&self.pending, self.cell_volume, cancelled)?;
        let mut heat = 0.0;
        let mut correction = 0.0;
        for (i, &value) in self.pending_heat_j.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            accumulate(&mut heat, &mut correction, value);
        }
        let total_heat = self.total_heat_j + heat;
        let defect = (after - before) - work - source_work + heat;
        let scale = before + after + work.abs() + source_work.abs() + heat;
        let tolerance = 32768.0 * f64::EPSILON * scale;
        if ![heat, total_heat, work, source_work, defect, tolerance]
            .into_iter()
            .all(f64::is_finite)
            || heat < 0.0
        {
            return Err(ChannelError::Invalid(
                "wave work or friction heat is not representable",
            ));
        }
        if defect.abs() > tolerance {
            return Err(ChannelError::Balance {
                defect_j: defect,
                tolerance_j: tolerance,
            });
        }
        let left = self.modes.face(&self.material, port_modes[0]);
        let right = self.modes.face(&self.material, port_modes[1]);
        let boundary_mass = dt * self.area_m2 * (left.mass_flux_kg_m2_s - right.mass_flux_kg_m2_s);
        let boundary_entropy =
            dt * self.area_m2 * (left.entropy_flux_w_m2_k - right.entropy_flux_w_m2_k);
        let boundary_heat = dt
            * self.area_m2
            * (left.counterflow_heat_flux_w_m2 - right.counterflow_heat_flux_w_m2);
        let mut pressure_impulse = [0.0; 2];
        for (end, cell) in [0, n - 1].into_iter().enumerate() {
            let old = self.modes.cell(
                &self.material,
                [self.accepted[0].cells[cell], self.accepted[1].cells[cell]],
            );
            let new = self.modes.cell(
                &self.material,
                [self.pending[0].cells[cell], self.pending[1].cells[cell]],
            );
            pressure_impulse[end] = dt
                * self.area_m2
                * (0.5 * old.pressure_perturbation_pa + 0.5 * new.pressure_perturbation_pa)
                * if end == 0 { 1.0 } else { -1.0 };
        }
        if ![
            boundary_mass,
            boundary_entropy,
            boundary_heat,
            pressure_impulse[0],
            pressure_impulse[1],
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(ChannelError::Invalid(
                "boundary transfer is not representable",
            ));
        }
        poll(cancelled)?;
        let momentum = self.audit_resolved_momentum(force, dt, pressure_impulse, cancelled)?;
        Ok(ChannelStep {
            time_s: next_time,
            wave_energy_before_j: before,
            wave_energy_after_j: after,
            boundary_wave_work_j: work,
            body_force_work_j: source_work,
            boundary_pressure_impulse_n_s: pressure_impulse,
            resolved_momentum_before_kg_m_s: momentum[0],
            resolved_momentum_after_kg_m_s: momentum[1],
            interior_body_impulse_n_s: momentum[2],
            resolved_momentum_defect_n_s: momentum[3],
            mutual_friction_heat_j: heat,
            wave_balance_defect_j: defect,
            boundary_mass_kg: boundary_mass,
            boundary_entropy_j_k: boundary_entropy,
            boundary_counterflow_heat_j: boundary_heat,
            maximum_relative_linear_residual: residual,
        })
    }

    fn friction_half_step(
        &mut self,
        friction: MutualFriction,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        if friction == MutualFriction::None {
            return Ok(());
        }
        for i in 1..self.cells() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let sample = self.modes.face(
                &self.material,
                [self.pending[0].faces[i], self.pending[1].faces[i]],
            );
            let velocities = TwoFluidVelocity {
                normal_m_s: [sample.normal_velocity_m_s, 0.0, 0.0],
                superfluid_m_s: [sample.superfluid_velocity_m_s, 0.0, 0.0],
            };
            let change = match friction {
                MutualFriction::None => unreachable!(),
                MutualFriction::IsotropicVortexTangle { .. } => {
                    unreachable!("vortex tangle drag is resolved once before staging")
                }
                MutualFriction::Linear { gamma_kg_m3_s } => {
                    self.material
                        .relax_linear_friction(velocities, gamma_kg_m3_s, dt)?
                }
                MutualFriction::GorterMellink { a_m_s_kg } => self
                    .material
                    .relax_gorter_mellink(velocities, a_m_s_kg, dt)?,
            };
            let v = self.material.barycentric_velocity_m_s(change.velocity)?[0];
            let w = self.material.counterflow_velocity_m_s(change.velocity)?[0];
            let next = self.modes.face_to_modes(v, w);
            self.pending[0].faces[i] = next[0];
            self.pending[1].faces[i] = next[1];
            let heat = change.dissipated_energy_j_m3 * self.cell_volume * 0.5;
            self.pending_heat_j[i - 1] += heat;
            self.pending_heat_j[i] += heat;
        }
        poll(cancelled)
    }
}

/// Original staggered equations, after eliminating the interior face velocities:
/// (I-rL) P_new = (I+rL) P_old - lambda D V_old, r=lambda²/4.
/// Boundary velocities are prescribed and have no phantom boundary inertia.
fn midpoint(
    old: &[f64],
    proposed: &mut Field,
    work: &mut SolveWorkspace,
    lambda: f64,
    dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ChannelError> {
    let n = old.len();
    let r = 0.25 * lambda * lambda;
    if work.factor_step != Some(dt) {
        work.factor_step = None;
        work.boundary_inverse = [None; 2];
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            work.diagonal[i] = 1.0 + if i == 0 || i + 1 == n { r } else { 2.0 * r };
            if i + 1 < n {
                work.off_diagonal[i] = -r;
            }
        }
        work.factor.factor(
            &work.off_diagonal,
            &work.diagonal,
            &work.off_diagonal,
            &mut *cancelled,
        )?;
        work.factor_step = Some(dt);
    }
    for i in 0..n {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let mut laplacian = 0.0;
        if i > 0 {
            laplacian += old[i - 1] - old[i];
        }
        if i + 1 < n {
            laplacian += old[i + 1] - old[i];
        }
        work.rhs[i] = old[i] + r * laplacian - lambda * (proposed.faces[i + 1] - proposed.faces[i]);
    }
    work.factor
        .solve_factored_into(&work.rhs, &mut proposed.cells, &mut *cancelled)?;
    let mut error: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for i in 0..n {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let mut applied = work.diagonal[i] * proposed.cells[i];
        let mut row_scale = applied.abs() + work.rhs[i].abs();
        if i > 0 {
            applied -= r * proposed.cells[i - 1];
            row_scale += (r * proposed.cells[i - 1]).abs();
        }
        if i + 1 < n {
            applied -= r * proposed.cells[i + 1];
            row_scale += (r * proposed.cells[i + 1]).abs();
        }
        error = error.max((applied - work.rhs[i]).abs());
        scale = scale.max(row_scale);
        if i > 0 {
            proposed.faces[i] -= 0.5
                * lambda
                * ((old[i] - old[i - 1]) + (proposed.cells[i] - proposed.cells[i - 1]));
        }
    }
    let relative = if scale == 0.0 { error } else { error / scale };
    if !relative.is_finite() || relative > 4096.0 * f64::EPSILON {
        return Err(ChannelError::Invalid(
            "two-fluid midpoint solve has an unresolved residual",
        ));
    }
    Ok(relative)
}
