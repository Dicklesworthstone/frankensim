//! Exact finite-sphere interactions with both native sources free to recoil.
//!
//! The spherical mass distributions may overlap or share a center; their
//! Newtonian interaction remains finite. This gravity model supplies neither
//! contact impulses nor deformation and does not authorize material overlap
//! for a mechanical contact owner. All native masses and inertias are retained.

mod pair;

use super::super::adaptive::{NativeStep, advance_native};
use super::super::{GravityAdvanceError, GravityAdvanceOptions, GravityAdvanceReceipt};
use super::{
    GravityError, MutualGravity, MutualGravityLoad, MutualGravityStepReceipt, Sum, VectorSum,
    binary_parts, poll, same_frame, scaled_positive, step_central, storage,
};
use crate::{MassProperties, RigidBodyState};

/// Uniform-density spheres in the original supplied body order.
///
/// Each declared positive radius must agree with all three native principal
/// moments `I = 2 m R^2 / 5` to a 128-epsilon relative roundoff allowance.
/// The source mass distribution is therefore consistent with its existing rigid
/// inertia, without introducing a second mass, pose or momentum owner.
///
/// The original frame, G and pair budget remain binding. For this finite model,
/// a zero separation cutoff allows coincident centers; a positive cutoff keeps
/// its strict exclusion meaning. Overlap is a gravity regime, not contact or
/// material interpenetration resolution. Self-energy within each sphere is
/// excluded from the pair interaction diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MutualSphericalGravity<'a> {
    monopole: MutualGravity,
    body_radii_m: &'a [f64],
}

impl MutualGravity {
    /// Select exact gravity between finite uniform-density spherical sources.
    ///
    /// The borrowed radii and original pair budget are admitted before any
    /// allocation. Native body counts, masses and uniform-sphere inertias are
    /// checked by each load/step under bounded per-body cancellation polling.
    /// A zero cutoff admits the finite coincident-center field in this lane.
    pub fn uniform_spheres<'a>(
        self,
        body_radii_m: &'a [f64],
    ) -> Result<MutualSphericalGravity<'a>, GravityError> {
        self.pair_count(body_radii_m.len())?;
        if body_radii_m
            .iter()
            .any(|&radius| !radius.is_finite() || radius <= 0.0)
        {
            return Err(GravityError::Invalid(
                "mutual spheres need positive finite source radii",
            ));
        }
        Ok(MutualSphericalGravity {
            monopole: self,
            body_radii_m,
        })
    }
}

impl MutualSphericalGravity<'_> {
    /// Evaluate the exact pair potential and equal/opposite COM forces.
    ///
    /// Exterior, partial-overlap and full-containment regimes share one
    /// conservative pair law, continuous at the regime boundaries. Spherical
    /// symmetry supplies no gravity torque about either COM. Both original
    /// source masses contribute and each unordered pair is counted once.
    /// Work is O(N^2), scratch O(N), with fixed body/pair reduction order and
    /// cancellation polled before allocation, per body/pair and at publication.
    pub fn load(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MutualGravityLoad, GravityError> {
        poll(&mut cancelled)?;
        same_frame(frame, self.monopole.frame)?;
        if states.len() != properties.len() || states.len() != self.body_radii_m.len() {
            return Err(GravityError::Invalid(
                "mutual spheres need matching body, property and radius counts",
            ));
        }
        let pair_count = self.monopole.pair_count(states.len())?;
        let mut forces = storage(states.len())?;
        for ((&state, &mass), &radius) in states.iter().zip(properties).zip(self.body_radii_m) {
            poll(&mut cancelled)?;
            state.validate()?;
            admit_uniform_inertia(mass, radius)?;
            forces.push(VectorSum::default());
        }
        let mut potential = Sum::default();
        for i in 0..states.len() {
            for j in i + 1..states.len() {
                poll(&mut cancelled)?;
                let displacement = states[j]
                    .pose()
                    .position_world()
                    .sub(states[i].pose().position_world());
                if self.monopole.minimum_separation_m > 0.0 {
                    self.admit_radius(displacement.stable_norm("mutual sphere separation")?)?;
                }
                let (force, energy) = pair::load(
                    self.monopole.gravitational_constant,
                    [properties[i].mass(), properties[j].mass()],
                    [self.body_radii_m[i], self.body_radii_m[j]],
                    displacement,
                )?;
                forces[i].add(force)?;
                forces[j].add(force.scale(-1.0))?;
                potential.add(energy)?;
            }
        }
        let mut forces_world = storage(states.len())?;
        for sum in forces {
            poll(&mut cancelled)?;
            forces_world.push(sum.value()?);
        }
        let potential_energy_j = potential.value()?;
        poll(&mut cancelled)?;
        Ok(MutualGravityLoad {
            forces_world,
            potential_energy_j,
            pair_count,
        })
    }

    /// Advance all finite sources with the original native central-gravity map.
    ///
    /// Simultaneous exact pair forces drive common half-kicks around the native
    /// free translation/rotation drift. Both sources recoil; their original
    /// isotropic spin has no gravity torque. No energy or momentum is projected.
    /// The fixed resolved map is second order and symmetric; its finite-step
    /// energy defect is not zero. Linear and orbital-plus-spin angular momenta
    /// are conserved to roundoff.
    ///
    /// Zero cutoff admits overlap, containment and center crossing without a
    /// point singularity. A positive cutoff still excludes both endpoints and
    /// every straight relative drift segment; it does not certify the true
    /// curved trajectory. Direct pair work retains both native displacements
    /// and avoids subtracting the absolute binding-energy diagnostics. Inputs
    /// are immutable and any refusal or cancellation returns no partial state.
    pub fn step(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        duration_seconds: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MutualGravityStepReceipt, GravityError> {
        let admission = (self.monopole.minimum_separation_m > 0.0)
            .then_some(|_, _, radius| self.admit_radius(radius));
        step_central(
            states,
            properties,
            duration_seconds,
            |states, cancelled| self.load(frame, states, properties, cancelled),
            admission,
            |i, j, before, after| {
                if self.monopole.minimum_separation_m > 0.0 {
                    // Preserve the typed, refinable domain refusal before any
                    // work arithmetic at an excluded endpoint. Initial loads
                    // already admitted the original configuration.
                    self.admit_radius(
                        after[1]
                            .sub(after[0])
                            .stable_norm("mutual sphere separation")?,
                    )?;
                }
                pair::potential_change(
                    self.monopole.gravitational_constant,
                    [properties[i].mass(), properties[j].mass()],
                    [self.body_radii_m[i], self.body_radii_m[j]],
                    before,
                    after,
                )
            },
            cancelled,
        )
    }

    /// Advance a bounded interval using the existing native adaptive driver.
    ///
    /// Accepted states use two native half steps, with finite-sphere forces and
    /// direct work in every trial. The original pair budget bounds each pass;
    /// `maximum_trials` bounds timestep refinement. Invalid inputs, arithmetic
    /// refusal, budget exhaustion and cancellation publish no partial interval.
    /// A later typed positive-cutoff refusal may trigger bounded refinement.
    /// Local-error estimates do not certify the true path or its global error,
    /// and the variable-step sequence is not claimed symplectic or reversible.
    pub fn advance(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        duration_seconds: f64,
        options: GravityAdvanceOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GravityAdvanceReceipt, GravityAdvanceError> {
        options.validate(duration_seconds)?;
        self.load(frame, states, properties, &mut cancelled)?;
        let mut initial = storage(states.len())?;
        for &state in states {
            poll(&mut cancelled)?;
            initial.push(state);
        }
        advance_native(
            initial,
            properties,
            duration_seconds,
            options,
            |states, duration, cancelled| {
                self.step(frame, states, properties, duration, cancelled)
                    .map(NativeStep::from)
            },
            cancelled,
        )
    }

    fn admit_radius(self, radius: f64) -> Result<(), GravityError> {
        if self.monopole.minimum_separation_m > 0.0 {
            self.monopole.admit_radius(radius)?;
        }
        Ok(())
    }
}

fn admit_uniform_inertia(properties: MassProperties, radius: f64) -> Result<(), GravityError> {
    properties.validate()?;
    // Avoid forming m*R or R^2: native moments can be representable even if
    // either intermediate overflows or underflows at the declared SI scales.
    let (mass, em) = binary_parts(properties.mass());
    let (radius, er) = binary_parts(radius);
    let expected = scaled_positive((0.4 * mass * radius) * radius, em + 2 * er)?;
    let moments = properties.principal_inertia_body();
    if [moments.x, moments.y, moments.z]
        .into_iter()
        .any(|moment| (moment / expected - 1.0).abs() > 128.0 * f64::EPSILON)
    {
        return Err(GravityError::Invalid(
            "native inertia does not match the declared uniform-density sphere",
        ));
    }
    Ok(())
}
