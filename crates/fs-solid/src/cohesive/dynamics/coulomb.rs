//! Rigid static/kinetic friction on the original damaged interface area.
//!
//! Two scaled tangent-traction unknowns per selected patch join the ORIGINAL
//! midpoint Newton/FGMRES solve. The shared tangent-disk map enforces Coulomb
//! sticking/sliding without a friction spring or a velocity regularization.
//! Capacity uses the original analytic integral of damage times compression
//! over the complete cohesive path. Static capture which slides is re-solved
//! with the sourced kinetic coefficient; every solve shares one work budget.
//! This is a discrete straight-path model. Resolved direction reversals refuse;
//! event times, redundant reaction uniqueness and trajectory accuracy are not
//! inferred from nonlinear convergence.

use super::{
    CohesiveMotion, CohesiveSystem, DynamicsConfig, DynamicsError, Error, finite, positive,
};
use crate::cohesive::CohesivePatch;
use fs_exec::Cx;
use fs_material::cohesive::friction::project_tangent_disk;
use fs_material::cohesive::source::{
    RESIDUAL_KINETIC_FRICTION_PROPERTY, RESIDUAL_STATIC_FRICTION_PROPERTY,
};
use fs_qty::{Dims, QuantitySpec};
use std::sync::Arc;

/// Explicit numerical conditioning and physical-motion resolution gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombConfig {
    /// Converts the two additional metre-valued coordinates to traction [Pa].
    /// Also scales the natural map; it adds no physical spring or stored energy.
    pub traction_scale_pa_per_m: f64,
    /// Absolute accepted net-slip resolution for static capture [m].
    pub slip_tolerance_m: f64,
    /// Resolved endpoint velocity threshold for branch occupancy [m/s].
    pub velocity_tolerance_m_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoulombActivation {
    AtRest,
    Sliding,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoulombRegime {
    Inactive,
    Sticking,
    Sliding,
}

/// Accepted next-step occupancy and previous mean reaction (a solve predictor,
/// not elastic storage or an extra slip/damage history).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombState {
    pub occupancy: CoulombRegime,
    pub last_traction_pa: [f64; 2],
}

#[derive(Clone, Copy)]
struct Model {
    static_mu: f64,
    kinetic_mu: f64,
    coordinate: Option<usize>,
}

#[derive(Clone)]
pub(super) struct CoulombBinding {
    models: Arc<[Option<Model>]>,
    pub(super) states: Arc<[CoulombState]>,
    // An instantaneous change invalidates inherited sliding occupancy only on
    // affected tangential traces. Reclassify using the next interval's actual
    // incoming velocity, including prescribed motion; keep prior reactions.
    impulse_recheck: Option<Arc<[bool]>>,
    config: CoulombConfig,
    extra: usize,
}

impl Model {
    fn from_patch(patch: &CohesivePatch<'_>) -> Result<Self, Error> {
        let source = patch.source.ok_or(Error::Invalid(
            "Coulomb friction needs the original sourced bond",
        ))?;
        let coefficient = |name| {
            let property = source.resolved().property(name).ok_or(Error::Invalid(
                "original bond bundle lacks static/kinetic residual friction",
            ))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !property.is_static_at_fixed_external_conditions()
            {
                return Err(Error::Invalid(
                    "Coulomb coefficients must be dimensionless and fixed at the original query",
                ));
            }
            finite(property.value_si())
        };
        let static_mu = coefficient(RESIDUAL_STATIC_FRICTION_PROPERTY)?;
        let kinetic_mu = coefficient(RESIDUAL_KINETIC_FRICTION_PROPERTY)?;
        if kinetic_mu < 0.0 || static_mu < kinetic_mu {
            return Err(Error::Invalid("require 0 <= kinetic mu <= static mu"));
        }
        Ok(Self {
            static_mu,
            kinetic_mu,
            coordinate: None,
        })
    }
}

impl<'m, 'a> CohesiveMotion<'m, 'a> {
    /// Select original sourced static/kinetic residual friction once. `None`
    /// leaves a patch frictionless. This is mutually exclusive with the optional
    /// velocity-regularized owner. Initial closed-contact occupancy is explicit.
    pub fn with_damaged_area_coulomb(
        self,
        activation: &[Option<CoulombActivation>],
        config: CoulombConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, DynamicsError> {
        self.with_damaged_area_coulomb_controlled(activation, config, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate counterpart; cancellation returns no partially adopted motion.
    pub fn with_damaged_area_coulomb_controlled(
        mut self,
        activation: &[Option<CoulombActivation>],
        config: CoulombConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DynamicsError> {
        poll(cancelled)?;
        if self.friction.is_some()
            || self.coulomb.is_some()
            || activation.len() != self.state.points.len()
        {
            return Err(Error::Invalid(
                "select one friction owner with one choice per original patch",
            )
            .into());
        }
        positive(config.traction_scale_pa_per_m)?;
        positive(config.slip_tolerance_m)?;
        positive(config.velocity_tolerance_m_s)?;
        let mut models = Vec::with_capacity(activation.len());
        let mut states = Vec::with_capacity(activation.len());
        for (p, choice) in activation.iter().enumerate() {
            poll(cancelled)?;
            models.push(
                choice
                    .map(|_| Model::from_patch(&self.state.system.patches[p]))
                    .transpose()?,
            );
            let point = self.state.points[p];
            let occupancy = if choice.is_none() || point.damage() == 0.0 || point.jump_m()[0] >= 0.0
            {
                CoulombRegime::Inactive
            } else if *choice == Some(CoulombActivation::Sliding) {
                CoulombRegime::Sliding
            } else {
                CoulombRegime::Sticking
            };
            states.push(CoulombState {
                occupancy,
                last_traction_pa: [0.0; 2],
            });
        }
        let binding = CoulombBinding::new(models, states.into(), config)?;
        binding.admit_basis(self.state.system, self.config)?;
        poll(cancelled)?;
        self.coulomb = Some(binding);
        Ok(self)
    }

    /// Same original patch order, including inactive/unselected entries.
    #[must_use]
    pub fn damaged_area_coulomb_states(&self) -> Option<&[CoulombState]> {
        self.coulomb.as_ref().map(|binding| binding.states.as_ref())
    }
}

impl CoulombBinding {
    fn new(
        mut models: Vec<Option<Model>>,
        states: Arc<[CoulombState]>,
        config: CoulombConfig,
    ) -> Result<Self, Error> {
        let mut extra = 0usize;
        for model in models.iter_mut().flatten() {
            if model.static_mu != 0.0 {
                model.coordinate = Some(extra);
                extra = extra
                    .checked_add(2)
                    .ok_or(Error::Budget("Coulomb coordinates"))?;
            }
        }
        Ok(Self {
            models: models.into(),
            states,
            impulse_recheck: None,
            config,
            extra,
        })
    }
    pub(super) fn dimension(&self, n: usize) -> Result<usize, Error> {
        n.checked_add(self.extra)
            .ok_or(Error::Budget("Coulomb dimension"))
    }
    pub(super) fn admit_basis(
        &self,
        system: &CohesiveSystem<'_>,
        config: DynamicsConfig,
    ) -> Result<(), Error> {
        let n = self.dimension(system.n)?;
        let r = config.solve.newton.linear_restart;
        let entries = r
            .checked_add(1)
            .and_then(|a| {
                n.checked_mul(2)
                    .and_then(|b| b.checked_add(r))
                    .and_then(|b| a.checked_mul(b))
            })
            .ok_or(Error::Budget("augmented Coulomb Krylov basis"))?;
        if entries > system.budget.max_krylov_basis_entries {
            return Err(Error::Budget("augmented Coulomb Krylov basis"));
        }
        Ok(())
    }
    pub(super) fn append_predictor(
        &self,
        motion: &CohesiveMotion<'_, '_>,
        x: &mut Vec<f64>,
        offsets: &[[f64; 3]],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let n = x.len();
        x.resize(self.dimension(n)?, 0.0);
        for (p, model) in self.models.iter().enumerate() {
            poll(cancelled)?;
            if let Some(model) = model
                && let Some(c) = model.coordinate
            {
                let mut traction = self.states[p].last_traction_pa;
                let incoming = local_velocity(motion, p, offsets[p], dt, &motion.velocity_m_s)?;
                if incoming[0].hypot(incoming[1]) > self.config.velocity_tolerance_m_s {
                    let point = motion.state.points[p];
                    let jump = motion.state.system.patches[p].jump(&x[..n], offsets[p])?;
                    let old = point.jump_m();
                    let slip = [finite(jump[1] - old[1])?, finite(jump[2] - old[2])?];
                    let distance = finite(slip[0].hypot(slip[1]))?;
                    if distance > 0.0 {
                        // Seed resolved motion on the sourced kinetic disk at
                        // the incoming damaged compression. Starting at zero
                        // reaction can trap a tiny-slip natural map's line
                        // search on the disk interior. This is only a solve
                        // guess: the original path pressure, equations and
                        // physical replay still determine every accepted force.
                        let pressure = finite(point.damage() * (-point.traction_pa()[0]).max(0.0))?;
                        let capacity = finite(model.kinetic_mu * pressure)?;
                        for d in 0..2 {
                            traction[d] = finite(capacity * (slip[d] / distance))?;
                        }
                    }
                }
                for d in 0..2 {
                    x[n + c + d] = finite(traction[d] / self.config.traction_scale_pa_per_m)?;
                }
            }
        }
        Ok(())
    }
    pub(super) fn needs_pressure(&self, p: usize) -> bool {
        self.models[p].is_some_and(|m| m.coordinate.is_some())
    }
    pub(super) fn coordinate(&self, p: usize) -> Option<usize> {
        self.models[p].and_then(|m| m.coordinate)
    }
    pub(super) fn traction(&self, p: usize, coordinates: &[f64]) -> Result<[f64; 2], Error> {
        let c = self
            .coordinate(p)
            .ok_or(Error::Invalid("missing selected Coulomb coordinates"))?;
        Ok([
            finite(coordinates[c] * self.config.traction_scale_pa_per_m)?,
            finite(coordinates[c + 1] * self.config.traction_scale_pa_per_m)?,
        ])
    }
    pub(super) fn coefficients(
        &self,
        motion: &CohesiveMotion<'_, '_>,
        offsets: &[[f64; 3]],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, Error> {
        let mut result = Vec::with_capacity(self.models.len());
        for (p, model) in self.models.iter().enumerate() {
            poll(cancelled)?;
            let mu = if let Some(model) = model {
                let v = local_velocity(motion, p, offsets[p], dt, &motion.velocity_m_s)?;
                let inherited_sliding = self.states[p].occupancy == CoulombRegime::Sliding
                    && !self.impulse_recheck.as_ref().is_some_and(|flags| flags[p]);
                if inherited_sliding || v[0].hypot(v[1]) > self.config.velocity_tolerance_m_s {
                    model.kinetic_mu
                } else {
                    model.static_mu
                }
            } else {
                0.0
            };
            result.push(mu);
        }
        Ok(result)
    }
    pub(super) fn switch_sliding(
        &self,
        motion: &CohesiveMotion<'_, '_>,
        x: &[f64],
        offsets: &[[f64; 3]],
        dt: f64,
        endpoint_velocity: &[f64],
        mu: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Error> {
        let mut changed = false;
        for (p, model) in self.models.iter().enumerate() {
            poll(cancelled)?;
            if let Some(model) = model
                && mu[p] != model.kinetic_mu
            {
                let jump = motion.state.system.patches[p].jump(x, offsets[p])?;
                let old = motion.state.points[p].jump_m();
                let v1 = local_velocity(motion, p, offsets[p], dt, endpoint_velocity)?;
                if (jump[1] - old[1]).hypot(jump[2] - old[2]) > self.config.slip_tolerance_m
                    || v1[0].hypot(v1[1]) > self.config.velocity_tolerance_m_s
                {
                    mu[p] = model.kinetic_mu;
                    changed = true;
                }
            }
        }
        Ok(changed)
    }
    pub(super) fn residual(
        &self,
        patch: u64,
        force: [f64; 2],
        slip: [f64; 2],
        pressure: f64,
        mu: f64,
    ) -> Result<[f64; 2], Error> {
        let cap = finite(mu * pressure)?;
        let k = self.config.traction_scale_pa_per_m;
        let q = [
            finite(force[0] + k * slip[0])?,
            finite(force[1] + k * slip[1])?,
        ];
        let projected =
            project_tangent_disk(q, cap).map_err(|error| Error::Material { patch, error })?;
        // Keep the natural map continuous through the disk boundary: scaling
        // only its exterior creates a merit barrier to static breakaway.
        // The accepted replay separately checks the actual kinetic force and
        // work, including directional errors near zero slip.
        Ok([
            finite(force[0] - projected[0])?,
            finite(force[1] - projected[1])?,
        ])
    }

    /// Exact local derivative of the original natural map, before patch area.
    /// The first block differentiates the full endpoint jump; the second uses
    /// the additional metre-valued coordinates z, with traction = k z. Pressure
    /// and its gradient belong to the same original damaged-compression path.
    /// Neither the conditioning scale k nor this derivative adds a spring.
    ///
    /// The disk interior has zero traction-coordinate derivative and -k times
    /// the tangential jump derivative. Outside, differentiate both the disk
    /// direction and its pressure-dependent radius. The actual disk boundary
    /// has no classical derivative and requests the existing numerical fallback.
    pub(super) fn natural_map_tangent(
        &self,
        patch: u64,
        force: [f64; 2],
        slip: [f64; 2],
        pressure: f64,
        pressure_gradient: [f64; 3],
        mu: f64,
    ) -> Result<Option<([[f64; 3]; 2], [[f64; 2]; 2])>, Error> {
        let capacity = finite(mu * pressure)?;
        let k = self.config.traction_scale_pa_per_m;
        let q = [
            finite(force[0] + k * slip[0])?,
            finite(force[1] + k * slip[1])?,
        ];
        let length = q[0].hypot(q[1]);
        if !length.is_finite() {
            return Ok(None);
        }
        // Retain the scalar owner's domain admission. A derivative refusal
        // must not introduce a different physical projection or repair inputs.
        project_tangent_disk(q, capacity).map_err(|error| Error::Material { patch, error })?;
        if pressure_gradient.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        if capacity == 0.0 {
            // For smooth nonnegative P at P=0, dP is zero and the projected
            // force is o(|dg|), even at q=0. A nonzero dP there is a one-sided
            // contact derivative; an underflowed positive radius is unresolved.
            if mu == 0.0 || (pressure == 0.0 && pressure_gradient.iter().all(|&v| v == 0.0)) {
                return Ok(Some(([[0.0; 3]; 2], [[k, 0.0], [0.0, k]])));
            }
            return Ok(None);
        }
        if length == capacity {
            return Ok(None);
        }
        if length < capacity {
            return Ok(Some(([[0.0, -k, 0.0], [0.0, 0.0, -k]], [[0.0; 2]; 2])));
        }
        let normal = q.map(|v| v / length);
        let transverse = [-normal[1], normal[0]];
        let ratio = capacity / length;
        let mut jump = [[0.0; 3]; 2];
        let mut coordinates = [[0.0; 2]; 2];
        for i in 0..2 {
            for (j, &gradient) in pressure_gradient.iter().enumerate() {
                jump[i][j] = -normal[i] * (mu * gradient);
            }
            for j in 0..2 {
                // In two dimensions, t t^T = I - n n^T. This form preserves
                // a small transverse derivative near an axis without 1-n_i².
                let projection = ratio * transverse[i] * transverse[j];
                jump[i][j + 1] -= k * projection;
                coordinates[i][j] = k * (if i == j { 1.0 } else { 0.0 } - projection);
            }
        }
        Ok(jump
            .iter()
            .flatten()
            .chain(coordinates.iter().flatten())
            .all(|v| v.is_finite())
            .then_some((jump, coordinates)))
    }

    pub(super) fn audit(
        &self,
        motion: &CohesiveMotion<'_, '_>,
        p: usize,
        offset: [f64; 3],
        dt: f64,
        endpoint_velocity: &[f64],
        point: fs_material::cohesive::CohesivePoint<'_>,
        force: [f64; 2],
        slip: [f64; 2],
        pressure: f64,
        capture_mu: f64,
    ) -> Result<(f64, CoulombState), Error> {
        let Some(model) = self.models[p] else {
            return Ok((
                0.0,
                CoulombState {
                    occupancy: CoulombRegime::Inactive,
                    last_traction_pa: [0.0; 2],
                },
            ));
        };
        let patch = &motion.state.system.patches[p];
        let tolerance = motion.config.solve.force_tolerance_n / patch.area_m2;
        let distance = slip[0].hypot(slip[1]);
        let v0 = local_velocity(motion, p, offset, dt, &motion.velocity_m_s)?;
        let v1 = local_velocity(motion, p, offset, dt, endpoint_velocity)?;
        let tol = self.config.velocity_tolerance_m_s;
        if pressure > 0.0
            && v0[0].hypot(v0[1]) > tol
            && v1[0].hypot(v1[1]) > tol
            && v0[0] * v1[0] + v0[1] * v1[1] < -tol * tol
        {
            return Err(Error::UnresolvedFrictionReversal { patch: patch.id });
        }
        // A resolved moving path still dissipates when dt makes its increment
        // smaller than the displacement resolution. Refinement cannot turn
        // physical kinetic work into zero-heat static capture.
        let static_step = distance <= self.config.slip_tolerance_m
            && v0[0].hypot(v0[1]) <= tol
            && v1[0].hypot(v1[1]) <= tol;
        let loss = if pressure == 0.0 {
            if force[0].hypot(force[1]) > tolerance {
                return Err(Error::Invalid("open Coulomb path carries friction"));
            }
            0.0
        } else if static_step {
            if force[0].hypot(force[1]) > capture_mu * pressure + tolerance {
                return Err(Error::Invalid(
                    "static Coulomb reaction exceeds capture cone",
                ));
            }
            0.0
        } else {
            positive(distance)?;
            let cap = finite(model.kinetic_mu * pressure)?;
            for d in 0..2 {
                if (force[d] - cap * (slip[d] / distance)).abs() > tolerance {
                    return Err(Error::Invalid(
                        "Coulomb sliding effort fails original kinetic law",
                    ));
                }
            }
            finite(cap * distance)?
        };
        let work = finite(force[0] * slip[0] + force[1] * slip[1])?;
        if (work - loss).abs() > motion.config.solve.local_energy_tolerance_j_m2 {
            return Err(Error::CohesiveEnergyDefect(work - loss));
        }
        let occupancy = if point.jump_m()[0] >= 0.0 || point.damage() == 0.0 {
            CoulombRegime::Inactive
        } else if v1[0].hypot(v1[1]) > tol {
            CoulombRegime::Sliding
        } else {
            CoulombRegime::Sticking
        };
        Ok((
            loss,
            CoulombState {
                occupancy,
                last_traction_pa: force,
            },
        ))
    }
    pub(super) fn accepted(&self, states: Vec<CoulombState>) -> Self {
        Self {
            states: states.into(),
            impulse_recheck: None,
            ..self.clone()
        }
    }

    /// Retain history and sources, but reconsider affected occupancy after an
    /// impulse. The prescribed rate is only known at the next finite interval,
    /// so classifying from the free velocity alone would be physically wrong.
    pub(super) fn after_impulse(
        &self,
        before: &CohesiveMotion<'_, '_>,
        velocity: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, super::impulse::ImpulseError> {
        use super::impulse::ImpulseError;
        let mut successor = self.clone();
        for (p, model) in self.models.iter().enumerate() {
            if cancelled() {
                return Err(ImpulseError::Cancelled);
            }
            if model.is_none() {
                continue;
            }
            let patch = &before.state.system.patches[p];
            let trace = |v: &[f64]| {
                patch.jump(v, [0.0; 3]).map_err(|_| {
                    ImpulseError::Invalid("non-finite retained contact velocity trace")
                })
            };
            let old = trace(&before.velocity_m_s)?;
            let new = trace(velocity)?;
            if old[1] != new[1] || old[2] != new[2] {
                let flags = successor
                    .impulse_recheck
                    .get_or_insert_with(|| vec![false; self.models.len()].into());
                Arc::make_mut(flags)[p] = true;
            }
        }
        Ok(successor)
    }
    /// Refinement may finish well inside the ORIGINAL velocity-resolution band.
    /// This tests an already audited native candidate; it changes no velocity,
    /// friction coefficient, history, reaction or physical acceptance tolerance.
    pub(super) fn resolved_stop(
        &self,
        before: &CohesiveMotion<'_, '_>,
        after: &CohesiveMotion<'_, '_>,
        dt: f64,
        patch_id: u64,
    ) -> Result<bool, Error> {
        let p = before
            .state
            .system
            .patches
            .iter()
            .position(|p| p.id == patch_id)
            .ok_or(Error::Invalid("reversal patch outside original system"))?;
        let free = before.state.system.patches[p].jump(&after.state.displacements, [0.0; 3])?;
        let jump = after.state.points[p].jump_m();
        let offset = [jump[0] - free[0], jump[1] - free[1], jump[2] - free[2]];
        let v = local_velocity(before, p, offset, dt, &after.velocity_m_s)?;
        // Leave a numerical margin from the native capture threshold so an
        // unresolved residual velocity does not flip branch from roundoff on
        // the next implicit step. Its actual kinetic energy remains accounted.
        Ok(v[0].hypot(v[1]) <= 0.25 * self.config.velocity_tolerance_m_s)
    }
    pub(super) fn temperature_successor(
        &self,
        target: &CohesiveSystem<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut models = Vec::with_capacity(self.models.len());
        for (old, patch) in self.models.iter().zip(&target.patches) {
            poll(cancelled)?;
            models.push(old.map(|_| Model::from_patch(patch)).transpose()?);
        }
        let mut successor = Self::new(models, self.states.clone(), self.config)?;
        successor.impulse_recheck = self.impulse_recheck.clone();
        Ok(successor)
    }
}

fn local_velocity(
    motion: &CohesiveMotion<'_, '_>,
    p: usize,
    offset: [f64; 3],
    dt: f64,
    velocity: &[f64],
) -> Result<[f64; 2], Error> {
    let patch = &motion.state.system.patches[p];
    let free0 = patch.jump(&motion.state.displacements, [0.0; 3])?;
    let old = motion.state.points[p].jump_m();
    let mut rate = [0.0; 3];
    for d in 0..3 {
        rate[d] = finite((offset[d] - (old[d] - free0[d])) / dt)?;
    }
    let value = patch.jump(velocity, rate)?;
    Ok([value[1], value[2]])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tangent_tests {
    use super::*;

    fn binding(k: f64) -> CoulombBinding {
        CoulombBinding {
            models: Vec::new().into(),
            states: Vec::new().into(),
            impulse_recheck: None,
            config: CoulombConfig {
                traction_scale_pa_per_m: k,
                slip_tolerance_m: 1e-10,
                velocity_tolerance_m_s: 1e-8,
            },
            extra: 0,
        }
    }
    fn close(a: f64, b: f64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
            "{a:e} != {b:e}"
        );
    }
    fn action(
        blocks: &([[f64; 3]; 2], [[f64; 2]; 2]),
        jump: [f64; 3],
        coordinates: [f64; 2],
    ) -> [f64; 2] {
        core::array::from_fn(|i| {
            (0..3).map(|j| blocks.0[i][j] * jump[j]).sum::<f64>()
                + (0..2).map(|j| blocks.1[i][j] * coordinates[j]).sum::<f64>()
        })
    }

    #[test]
    fn natural_map_blocks_match_independent_residual_differences() {
        for (k, force, slip, pressure, gradient, mu) in [
            (7.0, [0.1, -0.2], [0.01, 0.02], 3.0, [-2.5, 0.9, -0.4], 0.7),
            (7.0, [4.3, -1.7], [0.03, -0.07], 3.0, [-2.5, 0.9, -0.4], 0.4),
            (37.0, [-4.0, 2.0], [-0.2, 0.05], 2.0, [1.7, -0.5, 0.8], 0.5),
            (0.2, [1.001, 0.0], [0.0; 2], 2.0, [0.3, 0.2, -0.1], 0.5),
        ] {
            let model = binding(k);
            let original = model.residual(9, force, slip, pressure, mu).unwrap();
            let blocks = model
                .natural_map_tangent(9, force, slip, pressure, gradient, mu)
                .unwrap()
                .unwrap();
            let residual = |x: [f64; 5]| {
                let p =
                    pressure + (0..3).map(|j| gradient[j] * x[j]).sum::<f64>() + 0.3 * x[0] * x[1];
                model
                    .residual(
                        9,
                        [force[0] + k * x[3], force[1] + k * x[4]],
                        [slip[0] + x[1], slip[1] + x[2]],
                        p,
                        mu,
                    )
                    .unwrap()
            };
            for j in 0..5 {
                let h = 1e-6;
                let mut plus = [0.0; 5];
                let mut minus = [0.0; 5];
                plus[j] = h;
                minus[j] = -h;
                let a = residual(plus);
                let b = residual(minus);
                for i in 0..2 {
                    let expected = if j < 3 {
                        blocks.0[i][j]
                    } else {
                        blocks.1[i][j - 3]
                    };
                    close((a[i] - b[i]) / (2.0 * h), expected, 2e-8);
                }
            }
            let direction = [-0.6, 0.2, -0.4, 0.1, 0.3];
            let h = 1e-6;
            let a = residual(direction.map(|v| h * v));
            let b = residual(direction.map(|v| -h * v));
            let expected = action(&blocks, [-0.6, 0.2, -0.4], [0.1, 0.3]);
            for i in 0..2 {
                close((a[i] - b[i]) / (2.0 * h), expected[i], 2e-8);
            }
            assert_eq!(
                model.residual(9, force, slip, pressure, mu).unwrap(),
                original
            );
        }
    }

    #[test]
    fn disk_interior_is_exact_constraint_without_a_friction_spring() {
        let model = binding(13.0);
        let tangent = model
            .natural_map_tangent(1, [0.2, -0.1], [0.01, -0.02], 4.0, [3.0, -2.0, 1.0], 0.7)
            .unwrap()
            .unwrap();
        assert_eq!(tangent.0, [[0.0, -13.0, 0.0], [0.0, 0.0, -13.0]]);
        assert_eq!(tangent.1, [[0.0; 2]; 2]);
    }

    #[test]
    fn pressure_and_traction_blocks_rotate_with_the_same_tangent_frame() {
        let model = binding(7.0);
        let angle: f64 = 0.37;
        let (s, c) = angle.sin_cos();
        let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
        let rotate_jump = |v: [f64; 3]| {
            let t = rotate([v[1], v[2]]);
            [v[0], t[0], t[1]]
        };
        let force = [3.0, -4.0];
        let slip = [0.4, 0.2];
        let gradient = [-2.5, 0.9, -0.4];
        let old = model
            .natural_map_tangent(1, force, slip, 2.0, gradient, 0.4)
            .unwrap()
            .unwrap();
        let new = model
            .natural_map_tangent(
                1,
                rotate(force),
                rotate(slip),
                2.0,
                rotate_jump(gradient),
                0.4,
            )
            .unwrap()
            .unwrap();
        for j in 0..5 {
            let mut jump = [0.0; 3];
            let mut coordinates = [0.0; 2];
            if j < 3 {
                jump[j] = 1.0;
            } else {
                coordinates[j - 3] = 1.0;
            }
            let expected = rotate(action(&old, jump, coordinates));
            let actual = action(&new, rotate_jump(jump), rotate(coordinates));
            for i in 0..2 {
                close(actual[i], expected[i], 2e-13);
            }
        }
    }

    #[test]
    fn smooth_zero_radius_includes_zero_force_and_zero_slip() {
        let model = binding(7.0);
        for force in [[0.0; 2], [2.0, -3.0]] {
            for (pressure, gradient, mu) in [(3.0, [-2.0, 0.5, 0.3], 0.0), (0.0, [0.0; 3], 0.4)] {
                let blocks = model
                    .natural_map_tangent(1, force, [0.0; 2], pressure, gradient, mu)
                    .unwrap()
                    .unwrap();
                assert_eq!(blocks, ([[0.0; 3]; 2], [[7.0, 0.0], [0.0, 7.0]]));
            }
            // A smooth nonnegative pressure can leave zero quadratically. Its
            // projected effort is second order even when its direction flips.
            for j in 0..5 {
                let h = 1e-7;
                let evaluate = |sign: f64| {
                    let mut x = [0.0; 5];
                    x[j] = sign * h;
                    let p = x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2];
                    model
                        .residual(
                            1,
                            [force[0] + 7.0 * x[3], force[1] + 7.0 * x[4]],
                            [x[1], x[2]],
                            p,
                            0.4,
                        )
                        .unwrap()
                };
                let a = evaluate(1.0);
                let b = evaluate(-1.0);
                for i in 0..2 {
                    let expected = if j == i + 3 { 7.0 } else { 0.0 };
                    close((a[i] - b[i]) / (2.0 * h), expected, 2e-7);
                }
            }
        }
    }

    #[test]
    fn disk_corners_and_unrepresentable_derivatives_leave_original_residual_unchanged() {
        let model = binding(7.0);
        let huge = [0.9 * f64::MAX; 2];
        assert!(huge.iter().all(|v| v.is_finite()));
        assert!(!huge[0].hypot(huge[1]).is_finite());
        assert!(
            model
                .natural_map_tangent(1, huge, [0.0; 2], 1.0, [0.0; 3], 0.5)
                .unwrap()
                .is_none()
        );
        let before = model.residual(1, [3.0, 4.0], [0.0; 2], 10.0, 0.5).unwrap();
        assert!(
            model
                .natural_map_tangent(1, [3.0, 4.0], [0.0; 2], 10.0, [1.0; 3], 0.5)
                .unwrap()
                .is_none()
        );
        for force in [[0.0; 2], [2.0, -3.0]] {
            assert!(
                model
                    .natural_map_tangent(1, force, [0.0; 2], 0.0, [1.0, 0.0, 0.0], 0.4)
                    .unwrap()
                    .is_none()
            );
        }
        for gradient in [[f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]] {
            assert!(
                model
                    .natural_map_tangent(1, [2.0, 0.0], [0.0; 2], 1.0, gradient, 0.5)
                    .unwrap()
                    .is_none()
            );
        }
        assert!(
            model
                .natural_map_tangent(1, [2.0, 0.0], [0.0; 2], 1e-200, [1e200, 0.0, 0.0], 1e200)
                .unwrap()
                .is_none()
        );
        assert!(
            model
                .natural_map_tangent(1, [2.0, 0.0], [0.0; 2], 1e-300, [0.0; 3], 1e-300)
                .unwrap()
                .is_none()
        );
        assert_eq!(
            model.residual(1, [3.0, 4.0], [0.0; 2], 10.0, 0.5).unwrap(),
            before
        );
        assert!(
            model
                .residual(1, [f64::NAN, 0.0], [0.0; 2], 1.0, 0.5)
                .is_err()
        );
        assert!(
            model
                .natural_map_tangent(1, [f64::NAN, 0.0], [0.0; 2], 1.0, [0.0; 3], 0.5)
                .is_err()
        );
    }
}
