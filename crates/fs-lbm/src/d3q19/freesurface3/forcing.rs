//! Held local forces in the original Guo population frame.

use super::*;
use crate::d3q19::{RigidFrame3, W3, guo_projection3};

impl RigidFrame3 {
    /// Prepare the original sparse VOF owner for one prescribed rigid-frame step.
    ///
    /// All quantities use lattice units (`dx = dt = 1`) and moving-frame axes.
    /// Native cell centres are `[x+0.5, y+0.5, z+0.5]`. Fluid and wall velocities
    /// are relative to this frame. `acceleration` is physical body acceleration
    /// (including gravity), already resolved in the frame. Supply midpoint
    /// motion coefficients and refresh before EVERY subsequent native step.
    ///
    /// This replaces the COMPLETE held force field, including initial gravity.
    /// The shared dense-grid midpoint solve supplies translation, centrifugal,
    /// Euler and Coriolis terms; each wet cell's current density converts the
    /// resulting acceleration to force density. Pressure, capillarity and wall
    /// traction remain in the original free-surface operator.
    ///
    /// The existing Guo source-frame change preserves canonical density exactly,
    /// VOF fill/mass and physical velocity/stress to population roundoff. It
    /// applies no impulse and advances no time. New wet cells inherit the native
    /// neighbor-interpolated force frame during conversion; the next preparation
    /// evaluates their load at their actual position and accepted velocity.
    ///
    /// Only active wet cells require storage. Cancellation is polled per tile
    /// and before the atomic force/population publication. Geometry and mass
    /// routes do not change. This is a prescribed closed vessel with an imposed
    /// gas pressure, not finite-inertia vessel/gas dynamics or an energy-exact
    /// full-flow scheme. The original low-Mach and resolution limits still apply.
    /// Historical `LegacyGlobal` replay is refused because it does not use the
    /// current physical Guo frame for every gas-link reconstruction.
    /// A separately selected `FreeSurface3::rebind_body_acceleration` owner is
    /// also refused, including an owned zero acceleration. That mode owns its
    /// density-dependent endpoint source and cannot become a held force field
    /// implicitly. This entry never combines independent acceleration owners.
    pub fn prepare_free_surface_step(
        self,
        surface: &mut FreeSurface3,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.prepare_free_surface_step_with_acceleration(
            surface,
            acceleration,
            |_, _, _| Ok(acceleration),
            gate,
        )
    }

    // Constitutive owners can supply already-resolved local body acceleration
    // without duplicating tile traversal, frame inertia, or the Guo transform.
    // The nominal value is admitted even for an empty surface. Callback refusal
    // publishes no part of the new force field or its population representation.
    pub(crate) fn prepare_free_surface_step_with_acceleration<E>(
        self,
        surface: &mut FreeSurface3,
        nominal_acceleration: [f64; 3],
        mut acceleration_at: impl FnMut(u64, usize, f64) -> Result<[f64; 3], E>,
        gate: &CancelGate,
    ) -> Result<(), E>
    where
        E: From<FreeSurfaceError3>,
    {
        check_cancel(gate)?;
        let denominator = self
            .midpoint_denominator(nominal_acceleration)
            .map_err(FreeSurfaceError3::Transport)?;
        if surface.redistribution == ConversionRedistribution3::LegacyGlobal {
            return Err(FreeSurfaceError3::Transport(
                "rigid-frame flow requires current free-surface force semantics",
            )
            .into());
        }
        if surface.body_acceleration().is_some() {
            return Err(FreeSurfaceError3::Transport(
                "rigid-frame forcing cannot replace an independent body-acceleration owner",
            )
            .into());
        }
        if !surface.collision_model().supports_body_force() {
            return Err(FreeSurfaceError3::Transport(
                "rigid-frame free surface requires a forcing-compatible collision model",
            )
            .into());
        }
        let mut forces = BTreeMap::new();
        for (&key, tile) in &surface.tiles {
            check_cancel(gate)?;
            let slot = surface.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let (rho, raw) = surface.grid.cell_macros(slot, lane);
                let previous = surface.force_at_lane(key, lane);
                let velocity = core::array::from_fn(|axis| raw[axis] + 0.5 * previous[axis] / rho);
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let position = [x, y, z].map(|coordinate| coordinate as f64 + 0.5);
                let acceleration = acceleration_at(key, lane, rho)?;
                let effective =
                    self.effective_acceleration(position, velocity, acceleration, denominator);
                let next = effective.map(|value| rho * value);
                if !rho.is_finite()
                    || rho <= 0.0
                    || !previous
                        .into_iter()
                        .chain(velocity)
                        .chain(next)
                        .all(f64::is_finite)
                {
                    return Err(FreeSurfaceError3::Transport(
                        "unrepresentable rigid-frame free-surface cell",
                    )
                    .into());
                }
                forces.insert((key, lane), next);
            }
        }
        surface.rebind_force_field(forces, gate)?;
        Ok(())
    }

    /// Enclose source preparation and the native timestep in one transaction.
    /// The SI owner delegates here so momentum-preserving conversion uses this
    /// candidate directly rather than cloning an already-private sparse grid.
    pub(crate) fn step_free_surface(
        self,
        surface: &mut FreeSurface3,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let mut candidate = surface.clone();
        self.prepare_free_surface_step(&mut candidate, acceleration, gate)?;
        candidate.step_inner(gate, None, None, None)?;
        check_cancel(gate)?;
        *surface = candidate;
        Ok(())
    }
}

#[cfg(test)]
mod acceleration_tests;

pub(super) fn body_force(rho: f64, acceleration: [f64; 3]) -> Result<[f64; 3], FreeSurfaceError3> {
    let force = acceleration.map(|value| rho * value);
    if !rho.is_finite()
        || rho <= 0.0
        || force
            .iter()
            .zip(acceleration)
            .any(|(&value, a)| !value.is_finite() || a != 0.0 && value == 0.0)
    {
        return Err(FreeSurfaceError3::Transport(
            "body acceleration requires positive density and representable rho*a",
        ));
    }
    Ok(force)
}

/// Preserve the existing explicit physical boundary velocity while evaluating
/// its even Guo pressure source in the endpoint density-dependent force frame.
/// Missing gas populations obey f_q=base_q-rho*A_q; all other A_q are zero.
/// Thus rho=sum(base)/(1+sum(A)), without an iteration or a force reframe.
pub(super) fn close_body_pressure(
    populations: &mut [f64; Q3],
    coefficients: [f64; Q3],
) -> Result<(), FreeSurfaceError3> {
    let denominator = 1.0 + coefficients.iter().sum::<f64>();
    let rho = populations.iter().sum::<f64>() / denominator;
    if coefficients.iter().any(|value| !value.is_finite())
        || !denominator.is_finite()
        || denominator <= 0.0
        || !rho.is_finite()
        || rho <= 0.0
    {
        return Err(FreeSurfaceError3::Transport(
            "density-dependent gas pressure has no positive representable closure",
        ));
    }
    for (population, coefficient) in populations.iter_mut().zip(coefficients) {
        *population -= rho * coefficient;
    }
    let actual: f64 = populations.iter().sum();
    let allowance = 64.0 * f64::EPSILON * populations.iter().map(|value| value.abs()).sum::<f64>();
    if !allowance.is_finite() || (actual - rho).abs() > allowance {
        return Err(FreeSurfaceError3::Transport(
            "density-dependent pressure closure exceeds population roundoff",
        ));
    }
    preserve_density(populations, rho)?;
    if populations
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(FreeSurfaceError3::Transport(
            "body acceleration needs positive finite endpoint populations",
        ));
    }
    Ok(())
}

pub(super) fn source(q: usize, velocity: [f64; 3], force: [f64; 3]) -> f64 {
    let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
    let eu = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
    W3[q] * guo_projection3(e, velocity, eu, force)
}

// Preserve the actual canonical population reduction, not only its analytic
// zeroth moment. Even a one-ulp density change can turn a full interface into
// a fictitious partial cell. The rest correction changes no momentum/stress;
// try its adjacent representations before closing the final ordered add.
pub(super) fn preserve_density(
    populations: &mut [f64; Q3],
    rho: f64,
) -> Result<(), FreeSurfaceError3> {
    populations[0] += rho - populations.iter().sum::<f64>();
    let corrected_rest = populations[0];
    let actual = populations.iter().sum::<f64>();
    if actual == rho && rho.is_finite() {
        return Ok(());
    }
    let below = actual < rho;
    for _ in 0..8 {
        populations[0] = if below {
            populations[0].next_up()
        } else {
            populations[0].next_down()
        };
        let actual = populations.iter().sum::<f64>();
        if actual == rho && rho.is_finite() {
            return Ok(());
        }
        if !actual.is_finite() || (actual < rho) != below {
            break;
        }
    }
    // Rounded intermediate additions can skip the target for every nearby
    // rest value. Restore that correction and close the LAST add instead.
    // For physical positive populations prefix is near rho, so Sterbenz's
    // subtraction is exact; verify the resulting reduction even outside that
    // case. This moving-population correction affects momentum/stress
    // only at the explicitly bounded population-reduction roundoff scale.
    populations[0] = corrected_rest;
    let prefix = populations[..Q3 - 1].iter().sum::<f64>();
    let last = rho - prefix;
    let magnitude = populations.iter().map(|value| value.abs()).sum::<f64>();
    if !last.is_finite()
        || !magnitude.is_finite()
        || (last - populations[Q3 - 1]).abs() > 64.0 * f64::EPSILON * magnitude
    {
        return Err(FreeSurfaceError3::Transport(
            "force-frame density closure exceeds population roundoff",
        ));
    }
    populations[Q3 - 1] = last;
    if populations.iter().sum::<f64>() != rho {
        return Err(FreeSurfaceError3::Transport(
            "force-frame density is not exactly representable",
        ));
    }
    Ok(())
}

pub(super) fn reframe(
    mut populations: [f64; Q3],
    physical_velocity: [f64; 3],
    previous: [f64; 3],
    next: [f64; 3],
) -> Result<[f64; Q3], FreeSurfaceError3> {
    if previous == next {
        return Ok(populations);
    }
    let rho = populations.iter().sum::<f64>();
    let delta = core::array::from_fn(|d| previous[d] - next[d]);
    for (q, value) in populations.iter_mut().enumerate() {
        *value += 0.5 * source(q, physical_velocity, delta);
    }
    // The analytic source has zero mass. Close the same ordered reduction
    // used by native density and geometry, never changing VOF mass or heat.
    preserve_density(&mut populations, rho)?;
    if populations.iter().any(|v| !v.is_finite()) {
        return Err(FreeSurfaceError3::Transport(
            "unrepresentable force-frame change",
        ));
    }
    Ok(populations)
}

/// Change physical velocity while retaining density, pressure and physical
/// nonequilibrium stress in the same Guo frame. The external mass/momentum
/// owner supplies the impulse and accounts for its kinetic work.
pub(super) fn retarget_velocity(
    mut populations: [f64; Q3],
    previous: [f64; 3],
    next: [f64; 3],
    force: [f64; 3],
) -> Result<[f64; Q3], FreeSurfaceError3> {
    if previous == next {
        return Ok(populations);
    }
    let rho = populations.iter().sum::<f64>();
    let old_eq = equilibrium3(rho, previous);
    let new_eq = equilibrium3(rho, next);
    for (q, value) in populations.iter_mut().enumerate() {
        *value +=
            new_eq[q] - old_eq[q] + 0.5 * (source(q, previous, force) - source(q, next, force));
    }
    preserve_density(&mut populations, rho)?;
    if populations
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(FreeSurfaceError3::Transport(
            "inflow velocity needs positive finite populations",
        ));
    }
    Ok(populations)
}

impl FreeSurface3 {
    /// Explicit lattice body acceleration [dx/dt²], when selected. None is the
    /// historical uniform/reference-density or independently held local force.
    /// Some([0; 3]) remains an acceleration-owned zero load, not unknown forcing.
    #[must_use]
    pub const fn body_acceleration(&self) -> Option<[f64; 3]> {
        self.body_acceleration
    }

    /// Select or change a piecewise-held body acceleration with F=rho*a.
    ///
    /// This is a change of the accepted Guo representation, with no elapsed
    /// time or physical impulse. Each cell retains its canonical density,
    /// tracked mass, physical velocity and force-corrected stress. The next
    /// collision supplies F_n, and streaming publishes F_next=rho_next*a
    /// WITHOUT reframing populations: endpoint physical momentum includes the
    /// new half source. Gas pressure and interface birth use that same law.
    ///
    /// An existing independent local-force field (including thermal buoyancy)
    /// cannot be replaced implicitly. LegacyGlobal replay also refuses; its
    /// historical raw-velocity pressure/conversion convention is unchanged.
    /// Invalid loads, inadmissible populations or cancellation leave all
    /// accepted populations, force history and acceleration unchanged.
    pub fn rebind_body_acceleration(
        &mut self,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.rebind_body_acceleration_checked(acceleration, || check_cancel(gate))
    }

    /// Same candidate admission under an existing Cx/owner cancellation scope.
    pub(crate) fn rebind_body_acceleration_checked(
        &mut self,
        acceleration: [f64; 3],
        mut checkpoint: impl FnMut() -> Result<(), FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        checkpoint()?;
        if acceleration.iter().any(|value| !value.is_finite()) {
            return Err(FreeSurfaceError3::Transport(
                "body acceleration must be finite",
            ));
        }
        if self.redistribution == ConversionRedistribution3::LegacyGlobal {
            return Err(FreeSurfaceError3::Transport(
                "body acceleration cannot use historical LegacyGlobal force semantics",
            ));
        }
        if self.body_acceleration.is_none() && self.local_force.is_some() {
            return Err(FreeSurfaceError3::Transport(
                "body acceleration cannot replace an independent local-force owner",
            ));
        }
        let mut forces = BTreeMap::new();
        for (&key, tile) in &self.tiles {
            checkpoint()?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let populations = self.grid.populations(slot, lane);
                if populations
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.0)
                {
                    return Err(FreeSurfaceError3::Transport(
                        "body acceleration requires positive finite accepted populations",
                    ));
                }
                forces.insert(
                    (key, lane),
                    body_force(populations.iter().sum(), acceleration)?,
                );
            }
        }
        self.rebind_force_field_inner(forces, true, checkpoint)?;
        self.body_acceleration = Some(acceleration);
        Ok(())
    }

    pub(super) fn admit_body_acceleration_mode(&self) -> Result<(), FreeSurfaceError3> {
        if self.body_acceleration.is_some()
            && self.redistribution == ConversionRedistribution3::LegacyGlobal
        {
            return Err(FreeSurfaceError3::Transport(
                "body acceleration cannot use historical LegacyGlobal force semantics",
            ));
        }
        Ok(())
    }

    pub(crate) fn reference_force(&self) -> [f64; 3] {
        self.grid.collision().1
    }

    pub(super) fn force_at_lane(&self, key: u64, lane: usize) -> [f64; 3] {
        self.local_force
            .as_ref()
            .map_or_else(|| self.grid.collision().1, |forces| forces[&(key, lane)])
    }

    /// Actual force density in lattice units at an accepted wet cell. Its
    /// half-source correction is already included by `sample_at`; adding it
    /// again would count the same force twice. Gas/wall/outside returns None.
    #[must_use]
    pub fn force_density_at(&self, x: usize, y: usize, z: usize) -> Option<[f64; 3]> {
        let (slot, lane) = self.grid.resolve_source(
            i64::try_from(x).ok()?,
            i64::try_from(y).ok()?,
            i64::try_from(z).ok()?,
        )?;
        let key = self.grid.active_keys()[slot];
        (self.tiles.get(&key)?.cells[lane] != Cell3::Gas).then(|| self.force_at_lane(key, lane))
    }

    /// Change only the source frame, preserving canonical density exactly and
    /// physical velocity/force-corrected second moment to population roundoff.
    /// The caller owns the enclosing candidate.
    /// No elapsed time, impulse, heat, mass route or mechanical work is added.
    pub(crate) fn rebind_force_field(
        &mut self,
        forces: BTreeMap<(u64, usize), [f64; 3]>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        if self.body_acceleration.is_some() {
            return Err(FreeSurfaceError3::Transport(
                "independent local forces cannot replace the body-acceleration owner",
            ));
        }
        self.rebind_force_field_inner(forces, false, || check_cancel(gate))
    }

    fn rebind_force_field_inner(
        &mut self,
        forces: BTreeMap<(u64, usize), [f64; 3]>,
        retain_uniform_history: bool,
        mut checkpoint: impl FnMut() -> Result<(), FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        checkpoint()?;
        let wet = self
            .tiles
            .values()
            .map(|tile| {
                tile.cells
                    .iter()
                    .filter(|&&cell| cell != Cell3::Gas)
                    .count()
            })
            .sum::<usize>();
        if forces.len() != wet {
            return Err(FreeSurfaceError3::Transport(
                "local force must cover each wet cell once",
            ));
        }
        let uniform = self.grid.collision().1;
        let mut replacements = Vec::new();
        for (&key, tile) in &self.tiles {
            checkpoint()?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let next = *forces
                    .get(&(key, lane))
                    .ok_or(FreeSurfaceError3::Transport(
                        "missing local force at a wet cell",
                    ))?;
                if next.iter().any(|f| !f.is_finite()) {
                    return Err(FreeSurfaceError3::Transport("nonfinite local force"));
                }
                let previous = self.force_at_lane(key, lane);
                // Preserve the original population bits when the force does
                // not change, including the constant-density gravity limit.
                if previous == next {
                    continue;
                }
                let populations = self.grid.populations(slot, lane);
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                if !rho.is_finite() || rho <= 0.0 {
                    return Err(FreeSurfaceError3::Transport(
                        "invalid density in force-frame change",
                    ));
                }
                let velocity = core::array::from_fn(|d| raw[d] + 0.5 * previous[d] / rho);
                let populations = reframe(populations, velocity, previous, next)?;
                if retain_uniform_history && populations.iter().any(|value| *value <= 0.0) {
                    return Err(FreeSurfaceError3::Transport(
                        "body acceleration needs positive reframed populations",
                    ));
                }
                replacements.push((slot, lane, populations));
            }
        }
        checkpoint()?;
        for (slot, lane, populations) in replacements {
            self.grid.set_populations(slot, lane, populations);
        }
        self.local_force =
            if !retain_uniform_history && forces.values().all(|&force| force == uniform) {
                None
            } else {
                Some(forces)
            };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Guo's independently known first and second moments: F and uF + Fu.
    // The physical moments therefore add half that source to stored f.
    fn physical_moments(
        populations: &[f64; Q3],
        velocity: [f64; 3],
        force: [f64; 3],
    ) -> ([f64; 3], [[f64; 3]; 3]) {
        let mut momentum = force.map(|v| 0.5 * v);
        let mut stress = core::array::from_fn(|i| {
            core::array::from_fn(|j| 0.5 * (velocity[i] * force[j] + velocity[j] * force[i]))
        });
        for (q, &population) in populations.iter().enumerate() {
            let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
            for i in 0..3 {
                momentum[i] += population * e[i];
                for j in 0..3 {
                    stress[i][j] += population * e[i] * e[j];
                }
            }
        }
        (momentum, stress)
    }

    #[test]
    fn g1_inflow_velocity_retarget_keeps_density_and_physical_nonequilibrium() {
        for rho in [0.4, 0.8, 1.3, 2.5] {
            let previous = [0.02, -0.01, 0.03];
            let next = [-0.01, 0.025, 0.015];
            let force = [0.0003, -0.0002, 0.0001];
            let mut physical = equilibrium3(rho, previous);
            for (q, value) in physical.iter_mut().enumerate() {
                *value += 0.001 * rho * W3[q] * f64::from(E3[q].0 * E3[q].1);
            }
            let initial = reframe(physical, previous, [0.0; 3], force).unwrap();
            let actual_rho = initial.iter().sum::<f64>();
            let after = retarget_velocity(initial, previous, next, force).unwrap();
            assert_eq!(actual_rho.to_bits(), after.iter().sum::<f64>().to_bits());
            let (momentum, stress) = physical_moments(&after, next, force);
            let (_, old_stress) = physical_moments(&initial, previous, force);
            for a in 0..3 {
                assert!((momentum[a] - actual_rho * next[a]).abs() < 2e-14 * rho);
                for b in 0..3 {
                    let expected_change =
                        actual_rho * (next[a] * next[b] - previous[a] * previous[b]);
                    assert!(
                        (stress[a][b] - old_stress[a][b] - expected_change).abs() < 2e-14 * rho
                    );
                }
            }
            assert_eq!(
                retarget_velocity(initial, previous, previous, force).unwrap(),
                initial
            );
        }
    }

    #[test]
    fn g1_force_reframe_preserves_exact_density_and_physical_moments() {
        let mut seed = 0x6a09_e667_f3bc_c909_u64;
        let mut uniform = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut tail_closures = 0;
        for _ in 0..512 {
            let nominal_density = 0.9 + 0.2 * uniform();
            let velocity = core::array::from_fn(|_| 0.1 * (uniform() - 0.5));
            let previous = core::array::from_fn(|_| 0.002 * (uniform() - 0.5));
            let next = core::array::from_fn(|_| 0.002 * (uniform() - 0.5));
            let initial = equilibrium3(nominal_density, velocity);
            let rho = initial.iter().sum::<f64>();
            let first = reframe(initial, velocity, [0.0; 3], previous).unwrap();
            let after = reframe(first, velocity, previous, next).unwrap();
            assert!(first.iter().chain(&after).all(|&value| value > 0.0));
            assert_eq!(first.iter().sum::<f64>().to_bits(), rho.to_bits());
            assert_eq!(after.iter().sum::<f64>().to_bits(), rho.to_bits());
            let raw_last = first[Q3 - 1]
                + 0.5
                    * source(
                        Q3 - 1,
                        velocity,
                        core::array::from_fn(|d| previous[d] - next[d]),
                    );
            tail_closures += usize::from(after[Q3 - 1].to_bits() != raw_last.to_bits());
            let (initial_momentum, initial_stress) = physical_moments(&initial, velocity, [0.0; 3]);
            let (momentum, stress) = physical_moments(&after, velocity, next);
            for i in 0..3 {
                assert!((momentum[i] - initial_momentum[i]).abs() <= 8.0 * f64::EPSILON * rho);
                for j in 0..3 {
                    assert!(
                        (stress[i][j] - initial_stress[i][j]).abs() <= 8.0 * f64::EPSILON * rho
                    );
                }
            }
            // A true full interface remains full, and every actual partial
            // volume keeps the same bit pattern when its stored mass is held.
            for fraction in [0.0, 1e-8, 0.2, 0.999999, 1.0] {
                let mass = fraction * rho;
                assert_eq!(
                    (mass / rho).to_bits(),
                    (mass / after.iter().sum::<f64>()).to_bits()
                );
            }
            let unchanged = reframe(after, velocity, next, next).unwrap();
            assert_eq!(unchanged.map(f64::to_bits), after.map(f64::to_bits));
        }
        assert!(
            tail_closures > 0,
            "exercise unattainable rest-only reductions"
        );
        eprintln!(
            "force-frame density: 512 positive low-Mach states, {tail_closures} final-add closures"
        );
    }
}
