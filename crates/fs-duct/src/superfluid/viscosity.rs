//! Exponentially fitted normal-wall dynamics on the original staggered wave owner.
//!
//! In energy coordinates, x' = -C div(y), y' = -C grad(x) - R y + f.
//! R damps only the normal component. Freeze the gradient at midpoint x,
//! integrate each face exactly, and use its interval mean in the cell equation.
//! Eliminating faces leaves two scalar tridiagonals after a 2x2 rotation.
//! This retains the clamped-normal fourth-sound limit without a drag time-step
//! restriction. Scalar Darcy drag assumes a quasisteady transverse profile.
//! The optional Stokes spectrum instead stores positive-weight normal modes;
//! their weighted phi2 changes only the same 2x2 cell operator, while their
//! history and extra kinetic energy remain in the original channel owner.

use super::*;
use fs_material::superfluid::{
    NormalViscousFactors, NormalViscousSpectrum, PreparedNormalViscousSpectrum,
};

/// Interior-face transverse velocities, scaled by sqrt(rho_n). The original
/// face field stores their weighted section mean. Keeping independent accepted
/// and pending memory makes response probes and rejected heat callbacks atomic.
#[derive(Debug)]
pub(super) struct NormalProfileMemory {
    pub spectrum: NormalViscousSpectrum,
    accepted: Vec<f64>,
    pending: Vec<f64>,
    mean_history: Vec<f64>,
    density_root: f64,
    prepared: Option<(f64, PreparedNormalViscousSpectrum)>,
}

impl NormalProfileMemory {
    pub fn new(spectrum: NormalViscousSpectrum, cells: usize) -> Result<Self, ChannelError> {
        let count = spectrum.modes().len();
        let entries = (cells - 1)
            .checked_mul(count)
            .filter(|n| *n <= isize::MAX as usize / core::mem::size_of::<f64>())
            .ok_or(ChannelError::Capacity)?;
        Ok(Self {
            spectrum,
            accepted: filled(entries)?,
            pending: filled(entries)?,
            mean_history: filled(cells - 1)?,
            density_root: spectrum.material().normal_density_kg_m3().sqrt(),
            prepared: None,
        })
    }

    pub fn initialize_pending(
        &mut self,
        normal: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        let count = self.spectrum.modes().len();
        for (face, row) in self.pending.chunks_exact_mut(count).enumerate() {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            row.fill(self.density_root * normal[face + 1]);
        }
        poll(cancelled)
    }

    pub fn initialize_steady_pending(
        &mut self,
        fields: &[Field; 2],
        modes: Modes,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        let coefficients = self.spectrum.modes();
        if coefficients
            .iter()
            .any(|mode| !positive(mode.decay_rate_s_inv))
        {
            return Err(ChannelError::Invalid(
                "steady normal profile requires positive viscous relaxation rates",
            ));
        }
        let slowest = coefficients
            .iter()
            .map(|mode| mode.decay_rate_s_inv)
            .fold(f64::INFINITY, f64::min);
        let mut normalization = 0.0;
        let mut correction = 0.0;
        for mode in coefficients {
            accumulate(
                &mut normalization,
                &mut correction,
                mode.weight * (slowest / mode.decay_rate_s_inv),
            );
        }
        if !positive(normalization) {
            return Err(ChannelError::Invalid(
                "steady normal profile normalization is unresolved",
            ));
        }
        for (face, row) in self
            .pending
            .chunks_exact_mut(coefficients.len())
            .enumerate()
        {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            let velocity = modes.face(
                self.spectrum.material(),
                [fields[0].faces[face + 1], fields[1].faces[face + 1]],
            );
            let mean = self.density_root * velocity.normal_velocity_m_s;
            for (value, mode) in row.iter_mut().zip(coefficients) {
                *value = mean * ((slowest / mode.decay_rate_s_inv) / normalization);
            }
        }
        poll(cancelled)
    }

    pub fn copy_pending(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        for (accepted, pending) in self.accepted.chunks(128).zip(self.pending.chunks_mut(128)) {
            poll(cancelled)?;
            pending.copy_from_slice(accepted);
        }
        Ok(())
    }

    pub fn publish(&mut self) {
        core::mem::swap(&mut self.accepted, &mut self.pending);
    }

    pub fn lift_mean(&mut self, face: usize, delta_normal_m_s: f64) {
        let count = self.spectrum.modes().len();
        let offset = (face - 1) * count;
        let increment = self.density_root * delta_normal_m_s;
        for value in &mut self.pending[offset..offset + count] {
            *value += increment;
        }
    }

    pub fn mode_velocity(&self, face: usize, mode: usize) -> Option<f64> {
        let count = self.spectrum.modes().len();
        if face == 0 || face > self.mean_history.len() || mode >= count {
            return None;
        }
        Some(self.accepted[(face - 1) * count + mode] / self.density_root)
    }

    pub fn excess_energy(
        &self,
        pending: bool,
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ChannelError> {
        let data = if pending {
            &self.pending
        } else {
            &self.accepted
        };
        let coefficients = self.spectrum.modes();
        let mut total = 0.0;
        let mut correction = 0.0;
        for (face, row) in data.chunks_exact(coefficients.len()).enumerate() {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            // Center about one mode before taking the weighted mean. Uniform
            // states have exactly zero variance, even for a large bulk speed.
            let anchor = row[0];
            let mut mean_offset = 0.0;
            let mut mean_correction = 0.0;
            for (&value, mode) in row.iter().zip(coefficients) {
                accumulate(
                    &mut mean_offset,
                    &mut mean_correction,
                    mode.weight * (value - anchor),
                );
            }
            for (&value, mode) in row.iter().zip(coefficients) {
                let deviation = (value - anchor) - mean_offset;
                accumulate(
                    &mut total,
                    &mut correction,
                    (0.5 * volume * mode.weight * deviation) * deviation,
                );
            }
        }
        if total.is_finite() && total >= 0.0 {
            Ok(total)
        } else {
            Err(ChannelError::Invalid(
                "normal profile kinetic energy is not representable",
            ))
        }
    }

    pub fn validate_pending(
        &self,
        fields: &[Field; 2],
        modes: Modes,
        speed_limit: f64,
        physical_limits: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        let coefficients = self.spectrum.modes();
        for (face, row) in self.pending.chunks_exact(coefficients.len()).enumerate() {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            let mut mean = 0.0;
            let mut mean_correction = 0.0;
            let mut scale = 0.0;
            let mut speed_squared = 0.0;
            for (&value, mode) in row.iter().zip(coefficients) {
                let velocity = value / self.density_root;
                accumulate(&mut mean, &mut mean_correction, mode.weight * value);
                scale += mode.weight * value.abs();
                speed_squared += (mode.weight * velocity) * velocity;
            }
            let physical = modes.face(
                self.spectrum.material(),
                [fields[0].faces[face + 1], fields[1].faces[face + 1]],
            );
            let target = self.density_root * physical.normal_velocity_m_s;
            // A clamped normal mean may vanish by cancellation while the
            // superfluid remains mobile. Its reconstruction error is set by
            // the uncancelled kinetic coordinates, not by tiny v_n alone.
            scale +=
                target.abs() + fields[0].faces[face + 1].abs() + fields[1].faces[face + 1].abs();
            if ![mean, scale, speed_squared, target]
                .into_iter()
                .all(f64::is_finite)
                || speed_squared < 0.0
                || (physical_limits && speed_squared.sqrt() > speed_limit)
            {
                return Err(ChannelError::Domain {
                    field: "normal profile RMS speed",
                    index: face + 1,
                });
            }
            if (mean - target).abs() > 16384.0 * f64::EPSILON * scale {
                return Err(ChannelError::Domain {
                    field: "normal profile mean",
                    index: face + 1,
                });
            }
        }
        poll(cancelled)
    }

    fn prepare(&mut self, dt: f64) -> Result<(), ChannelError> {
        if self.prepared.as_ref().is_none_or(|(step, _)| *step != dt) {
            self.prepared = Some((dt, self.spectrum.prepare(dt)?));
        }
        Ok(())
    }

    fn prepared(&self, dt: f64) -> Result<&PreparedNormalViscousSpectrum, ChannelError> {
        self.prepared
            .as_ref()
            .filter(|(step, _)| *step == dt)
            .map(|(_, factors)| factors)
            .ok_or(ChannelError::Invalid(
                "normal profile factors are not prepared for this interval",
            ))
    }

    fn phi2(&self, dt: f64) -> Result<f64, ChannelError> {
        let mut mean = 0.0;
        let mut correction = 0.0;
        for mode in self.prepared(dt)?.modes() {
            accumulate(&mut mean, &mut correction, mode.weight * mode.factors.phi2);
        }
        Ok(mean)
    }

    fn prepare_history(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        let prepared = &self
            .prepared
            .as_ref()
            .ok_or(ChannelError::Invalid("normal profile factors are missing"))?
            .1;
        let coefficients = prepared.modes();
        for (face, (row, history)) in self
            .pending
            .chunks_exact(coefficients.len())
            .zip(&mut self.mean_history)
            .enumerate()
        {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            *history = 0.0;
            let mut correction = 0.0;
            for (&value, mode) in row.iter().zip(coefficients) {
                accumulate(
                    history,
                    &mut correction,
                    mode.weight * mode.factors.phi1 * value,
                );
            }
        }
        poll(cancelled)
    }

    // Next mean, interval mean, heat density, signed impulse loss in normalized
    // coordinates, and its uncancelled scale. All per-mode factors are prepared
    // before entering the spatial loops.
    fn advance_face(&mut self, face: usize, increment: f64) -> Result<[f64; 5], ChannelError> {
        let prepared = &self
            .prepared
            .as_ref()
            .ok_or(ChannelError::Invalid("normal profile factors are missing"))?
            .1;
        let coefficients = prepared.modes();
        let offset = (face - 1) * coefficients.len();
        let row = &mut self.pending[offset..offset + coefficients.len()];
        let mut totals = [0.0; 5];
        let mut corrections = [0.0; 5];
        for (value, mode) in row.iter_mut().zip(coefficients) {
            let f = mode.factors;
            let next = f.decay * *value + f.phi1 * increment;
            let mean = f.phi1 * *value + f.phi2 * increment;
            let first = f.heat_initial_root * (*value + f.heat_shift * increment);
            let second = f.heat_forced_root * increment;
            let heat = (mode.weight * first) * first + (mode.weight * second) * second;
            let loss = mode.weight * (f.loss * *value + f.drag_forced * increment);
            let values = [
                mode.weight * next,
                mode.weight * mean,
                heat,
                loss,
                loss.abs(),
            ];
            if !values.into_iter().all(f64::is_finite) || !next.is_finite() {
                return Err(ChannelError::Domain {
                    field: "normal profile update",
                    index: face,
                });
            }
            for i in 0..5 {
                accumulate(&mut totals[i], &mut corrections[i], values[i]);
            }
            *value = next;
        }
        Ok(totals)
    }
}

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
        let factors = NormalViscousFactors::try_new(channel.normal_wall_rate_s_inv, dt)?;
        let decay = factors.decay;
        let loss = factors.loss;
        let phi1 = factors.phi1;
        let phi2 = match &channel.normal_profile {
            Some(profile) => profile.phi2(dt)?,
            None => factors.phi2,
        };
        let drag_forced = factors.drag_forced;
        let heat_initial_root = factors.heat_initial_root;
        let heat_shift = factors.heat_shift;
        let heat_forced_root = factors.heat_forced_root;

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
                heat_initial_root,
                heat_shift,
                heat_forced_root,
                r[0],
                r[1],
                cross,
            ]
            .into_iter()
            .all(f64::is_finite)
            || heat_initial_root < 0.0
            || heat_forced_root < 0.0
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
            heat_initial_root,
            heat_shift,
            heat_forced_root,
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
        profile: Option<&NormalProfileMemory>,
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
        let normal_history = match profile {
            Some(profile) => profile.mean_history[face - 1],
            None => self.phi1 * old[0],
        };
        self.normal_rotated_flux[mode] * (normal_history + self.phi2 * body[0])
            + self.superfluid_rotated_flux[mode] * (old[1] + 0.5 * body[1])
    }
}

impl TwoFluidChannel {
    pub(super) fn wave_cell_rotation(&self, dt: f64) -> Result<[[f64; 2]; 2], ChannelError> {
        if !self.uses_viscous_wave_operator() {
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
        if let Some(profile) = &mut self.normal_profile {
            profile.prepare(dt)?;
            profile.prepare_history(cancelled)?;
        }
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
            let mut left_flux = operator.base_flux(
                &self.pending,
                self.normal_profile.as_ref(),
                mode,
                0,
                dt,
                force,
            );
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
                let right_flux = operator.base_flux(
                    &self.pending,
                    self.normal_profile.as_ref(),
                    mode,
                    i + 1,
                    dt,
                    force,
                );
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
            let normal = if let Some(profile) = &mut self.normal_profile {
                profile.advance_face(i, increment[0])?
            } else {
                // The material owner supplies the same positive two-square
                // heat integral for the scalar Darcy closure and every
                // transverse Stokes mode. No kinetic-store subtraction.
                let first =
                    operator.heat_initial_root * (old[0] + operator.heat_shift * increment[0]);
                let second = operator.heat_forced_root * increment[0];
                let loss = operator.loss * old[0] + operator.drag_forced * increment[0];
                [
                    operator.decay * old[0] + operator.phi1 * increment[0],
                    operator.phi1 * old[0] + operator.phi2 * increment[0],
                    first * first + second * second,
                    loss,
                    loss.abs(),
                ]
            };
            let next = [normal[0], old[1] + increment[1]];
            let mean = [normal[1], old[1] + 0.5 * increment[1]];
            let y = [0, 1].map(|m| operator.normal[m] * next[0] + operator.superfluid[m] * next[1]);
            let heat = self.cell_volume * normal[2];
            let work = self.cell_volume * dot(body, mean);
            let impulse = -self.cell_volume * operator.density_root[0] * normal[3];
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
            impulse_scale += self.cell_volume * operator.density_root[0] * normal[4];
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
