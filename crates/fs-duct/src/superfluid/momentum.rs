//! Physical body impulses and closed-end half-cell reactions on the same wave owner.

use super::*;

/// Force on each component per TOTAL mixture volume [N/m³].
/// These are explicit external sources, not inferred mutual friction or viscosity.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct ComponentForceDensity {
    pub normal_n_m3: f64,
    pub superfluid_n_m3: f64,
}

/// Constant-in-time forces during one interval. Spatial fields contain N+1 faces.
/// `Acceleration(a)` applies the SAME acceleration to both components, using
/// their original reference densities; it adds no buoyancy or stratification.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub enum BodyForce<'a> {
    #[default]
    None,
    Acceleration(f64),
    Uniform(ComponentForceDensity),
    Faces(&'a [ComponentForceDensity]),
}

impl BodyForce<'_> {
    pub(super) fn validate(
        self,
        faces: usize,
        material: &TwoFluidState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        if let Self::Faces(values) = self {
            if values.len() != faces {
                return Err(ChannelError::Invalid(
                    "body-force field needs N+1 face values",
                ));
            }
        }
        let count = if matches!(self, Self::Faces(_)) {
            faces
        } else {
            1
        };
        for i in 0..count {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let (bulk, relative, total) = self.components(i, material);
            if ![bulk, relative, total].into_iter().all(f64::is_finite) {
                return Err(ChannelError::Domain {
                    field: "body force",
                    index: i,
                });
            }
        }
        Ok(())
    }

    // Barycentric acceleration, relative acceleration, and total force density.
    fn components(self, i: usize, material: &TwoFluidState) -> (f64, f64, f64) {
        let value = match self {
            Self::None => return (0.0, 0.0, 0.0),
            Self::Acceleration(a) => return (a, 0.0, material.density_kg_m3() * a),
            Self::Uniform(value) => value,
            Self::Faces(values) => values[i],
        };
        let total = value.normal_n_m3 + value.superfluid_n_m3;
        (
            total / material.density_kg_m3(),
            value.normal_n_m3 / material.normal_density_kg_m3()
                - value.superfluid_n_m3 / material.superfluid_density_kg_m3(),
            total,
        )
    }

    fn impulse_scale(self, i: usize, material: &TwoFluidState, dt: f64, volume: f64) -> f64 {
        let value = match self {
            Self::None => return 0.0,
            Self::Acceleration(a) => return (material.density_kg_m3() * a * dt * volume).abs(),
            Self::Uniform(value) => value,
            Self::Faces(values) => values[i],
        };
        (value.normal_n_m3 * dt * volume).abs() + (value.superfluid_n_m3 * dt * volume).abs()
    }
}

/// First-order longitudinal momentum on a closed, stationary reference domain.
/// Pressure loads are perturbations about an explicitly balanced reference
/// pressure; nonlinear advective flux and gravitational stratification are absent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedChannelStep {
    pub channel: ChannelStep,
    pub momentum_before_kg_m_s: f64,
    pub momentum_after_kg_m_s: f64,
    /// Full-domain force integral, including both endpoint half dual volumes.
    pub body_impulse_n_s: f64,
    /// Signed perturbation wall impulses ON THE FLUID. Wall reactions are opposite.
    pub wall_impulse_n_s: [f64; 2],
    pub momentum_defect_n_s: f64,
}

impl TwoFluidChannel {
    /// Momentum of the dynamically resolved interior faces. With stationary
    /// closed ends this equals the full-domain trapezoidal momentum integral.
    /// A moving mechanical termination must add its endpoint half-cell momentum.
    pub fn resolved_momentum_kg_m_s(&self) -> Result<f64, ChannelError> {
        Ok(self.momentum_of(&self.accepted, &mut || false)?.0)
    }

    /// Drive both components with authored physical forces, keeping both ends
    /// stationary and impermeable. An already moving end refuses: closing it
    /// needs a separately modeled impulse, not discarded endpoint momentum.
    pub fn step_forced_closed(
        &mut self,
        dt_s: f64,
        force: BodyForce<'_>,
        friction: MutualFriction,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedChannelStep, ChannelError> {
        let report = self.stage_forced_closed(dt_s, force, friction, &mut cancelled)?;
        self.publish_step(report.channel);
        Ok(report)
    }

    /// Prepare the original heat owner's private candidate before publishing
    /// force, wave, momentum and heat together. The callback must not mutate an
    /// external accepted owner; only install its returned candidate on success.
    pub fn step_forced_closed_with_heat<T, E>(
        &mut self,
        dt_s: f64,
        force: BodyForce<'_>,
        friction: MutualFriction,
        mut cancelled: impl FnMut() -> bool,
        prepare_heat: impl FnOnce(&ForcedChannelStep, &[f64]) -> Result<T, E>,
    ) -> Result<(ForcedChannelStep, T), HeatStepError<E>> {
        let report = self
            .stage_forced_closed(dt_s, force, friction, &mut cancelled)
            .map_err(HeatStepError::Channel)?;
        let candidate = prepare_heat(&report, &self.pending_heat_j).map_err(HeatStepError::Heat)?;
        poll(&mut cancelled).map_err(HeatStepError::Channel)?;
        self.publish_step(report.channel);
        Ok((report, candidate))
    }

    fn stage_forced_closed(
        &mut self,
        dt: f64,
        force: BodyForce<'_>,
        friction: MutualFriction,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForcedChannelStep, ChannelError> {
        poll(cancelled)?;
        let n = self.cells();
        if self
            .accepted
            .iter()
            .any(|mode| mode.faces[0] != 0.0 || mode.faces[n] != 0.0)
        {
            return Err(ChannelError::Invalid(
                "closed forcing requires already stationary end velocities",
            ));
        }
        let (before, before_scale) = self.momentum_of(&self.accepted, cancelled)?;
        let channel =
            self.stage_step(dt, BoundaryVelocities::CLOSED, friction, force, cancelled)?;
        let (after, after_scale) = self.momentum_of(&self.pending, cancelled)?;
        let mut body = 0.0;
        let mut correction = 0.0;
        let mut scale = before_scale + after_scale;
        let mut wall = channel.boundary_pressure_impulse_n_s;
        for i in 0..=n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let (_, _, density) = force.components(i, &self.material);
            let volume = self.cell_volume * if i == 0 || i == n { 0.5 } else { 1.0 };
            let impulse = density * dt * volume;
            accumulate(&mut body, &mut correction, impulse);
            // Equal-and-opposite component impulses cancel physically, but
            // their uncancelled scale still bounds floating-point summation.
            scale += force.impulse_scale(i, &self.material, dt, volume);
            if i == 0 {
                wall[0] -= impulse;
            }
            if i == n {
                wall[1] -= impulse;
            }
        }
        scale += wall[0].abs() + wall[1].abs();
        let defect = (after - before) - body - (wall[0] + wall[1]);
        let tolerance = 65536.0 * f64::EPSILON * scale;
        if ![body, wall[0], wall[1], defect, tolerance]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(ChannelError::Invalid(
                "closed momentum transfer is not representable",
            ));
        }
        if defect.abs() > tolerance {
            return Err(ChannelError::MomentumBalance {
                defect_n_s: defect,
                tolerance_n_s: tolerance,
            });
        }
        poll(cancelled)?;
        Ok(ForcedChannelStep {
            channel,
            momentum_before_kg_m_s: before,
            momentum_after_kg_m_s: after,
            body_impulse_n_s: body,
            wall_impulse_n_s: wall,
            momentum_defect_n_s: defect,
        })
    }

    pub(super) fn force_half_step(
        &mut self,
        force: BodyForce<'_>,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ChannelError> {
        if force == BodyForce::None {
            return Ok(0.0);
        }
        let mut work = 0.0;
        let mut correction = 0.0;
        for i in 1..self.cells() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let (a, relative_a, _) = force.components(i, &self.material);
            let dv = a * dt;
            let dw = relative_a * dt;
            if dv == 0.0 && dw == 0.0 {
                continue;
            }
            let old = self.modes.face(
                &self.material,
                [self.pending[0].faces[i], self.pending[1].faces[i]],
            );
            let v = old.barycentric_velocity_m_s + dv;
            let w = old.counterflow_velocity_m_s + dw;
            let speeds = [
                v + self.material.superfluid_fraction() * w,
                v - self.material.normal_fraction() * w,
            ];
            // These half-kick velocities are algebraic auxiliaries of the
            // affine midpoint solve. A pressure-balanced state may have zero
            // physical velocity even when either auxiliary is large. Apply
            // physical limits after the complete forced-wave block instead.
            if ![dv, dw, v, w].into_iter().all(f64::is_finite)
                || speeds.into_iter().any(|u| !u.is_finite())
            {
                return Err(ChannelError::Domain {
                    field: "forced component velocity",
                    index: i,
                });
            }
            // Exact finite impulse work, including its quadratic term. Do not
            // subtract two large kinetic energies or count this work as heat.
            let energy = self.cell_volume
                * (self.material.density_kg_m3()
                    * dv
                    * (0.5 * old.barycentric_velocity_m_s + 0.5 * v)
                    + self.material.reduced_density_kg_m3()
                        * dw
                        * (0.5 * old.counterflow_velocity_m_s + 0.5 * w));
            if !energy.is_finite() {
                return Err(ChannelError::Invalid(
                    "body-force work is not representable",
                ));
            }
            accumulate(&mut work, &mut correction, energy);
            let next = self.modes.face_to_modes(v, w);
            self.pending[0].faces[i] = next[0];
            self.pending[1].faces[i] = next[1];
        }
        poll(cancelled)?;
        if !work.is_finite() {
            return Err(ChannelError::Invalid(
                "body-force work is not representable",
            ));
        }
        Ok(work)
    }

    fn momentum_of(
        &self,
        fields: &[Field; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ChannelError> {
        let mut total = 0.0;
        let mut correction = 0.0;
        let mut scale = 0.0;
        for i in 1..self.cells() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let face = self
                .modes
                .face(&self.material, [fields[0].faces[i], fields[1].faces[i]]);
            let value = self.cell_volume * face.mass_flux_kg_m2_s;
            accumulate(&mut total, &mut correction, value);
            scale += self.cell_volume
                * (self.material.normal_density_kg_m3() * face.normal_velocity_m_s).abs()
                + self.cell_volume
                    * (self.material.superfluid_density_kg_m3() * face.superfluid_velocity_m_s)
                        .abs();
        }
        if !total.is_finite() || !scale.is_finite() {
            return Err(ChannelError::Invalid(
                "channel momentum is not representable",
            ));
        }
        Ok((total, scale))
    }
}
