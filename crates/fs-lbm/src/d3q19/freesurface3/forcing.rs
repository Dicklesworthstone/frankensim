//! Held local forces in the original Guo population frame.

use super::*;
use crate::d3q19::{Coriolis3, RigidFrame3, W3, guo_projection3};

/// Prescribed acceleration in the fixed Cartesian lattice coordinates.
///
/// Every wet cell uses its own density in `F = rho * a(x)`. The descriptor
/// covers the whole domain, including cells activated by advancing liquid,
/// without a dense field or a wet-only map. It is independent of velocity:
/// Velocity-dependent rotation uses `Coriolis3` or the complete `RigidFrame3`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BodyAccelerationField3 {
    /// Spatially uniform acceleration [dx/dt²]; preserves the scalar arithmetic.
    Uniform([f64; 3]),
    /// `acceleration + gradient * (x - origin)` at cell centre `x = i + 1/2`.
    /// The origin stays fixed when a piston physically relocates liquid.
    Affine {
        /// Reference position in lattice coordinates [dx].
        origin: [f64; 3],
        /// Acceleration at the reference position [dx/dt²].
        acceleration: [f64; 3],
        /// Row-major `gradient[component][coordinate]` [1/dt²].
        gradient: [[f64; 3]; 3],
    },
}

impl BodyAccelerationField3 {
    fn admit(&self) -> Result<(), FreeSurfaceError3> {
        let finite = match self {
            Self::Uniform(acceleration) => acceleration.iter().all(|v| v.is_finite()),
            Self::Affine {
                origin,
                acceleration,
                gradient,
            } => origin
                .iter()
                .chain(acceleration)
                .chain(gradient.iter().flatten())
                .all(|v| v.is_finite()),
        };
        if !finite {
            return Err(FreeSurfaceError3::Transport(
                "body-acceleration field coefficients must be finite",
            ));
        }
        Ok(())
    }

    // Parameters are admitted once at binding. Sampling adds no source query,
    // allocation or temporal averaging; the original cell order is retained.
    pub(super) fn at_cell(&self, key: u64, lane: usize) -> Result<[f64; 3], FreeSurfaceError3> {
        let (origin, acceleration, gradient) = match self {
            Self::Uniform(acceleration) => return Ok(*acceleration),
            Self::Affine {
                origin,
                acceleration,
                gradient,
            } => (origin, acceleration, gradient),
        };
        let (x, y, z) = FreeSurface3::global_of(key, lane);
        let displacement: [f64; 3] =
            core::array::from_fn(|axis| [x, y, z][axis] as f64 + 0.5 - origin[axis]);
        let mut value = *acceleration;
        for axis in 0..3 {
            for coordinate in 0..3 {
                // A zero gradient retains the uniform path's bits, including
                // signed zeros, without multiplying unused displacements.
                if gradient[axis][coordinate] != 0.0 {
                    value[axis] += gradient[axis][coordinate] * displacement[coordinate];
                }
            }
        }
        if value.iter().any(|v| !v.is_finite()) {
            return Err(FreeSurfaceError3::Transport(
                "body acceleration is unrepresentable at a wet cell position",
            ));
        }
        Ok(value)
    }
}

impl Coriolis3 {
    /// Prepare local rotating-flow forcing on the original sparse wet cells.
    ///
    /// The source is `acceleration - 2*omega cross u`, with centrifugal
    /// potential already absorbed into effective acceleration. All quantities
    /// use lattice units and the existing relative velocity axes. Refresh
    /// before each native step, including after wet-cell conversion.
    ///
    /// Replaces the complete held force, without an instantaneous impulse.
    /// Canonical density, fill, tracked mass, physical velocity and corrected
    /// stress retain the same guarantees as [`RigidFrame3::prepare_free_surface_step`].
    /// Historical replay and independent body-acceleration ownership refuse.
    /// Material-owned gravity is explicitly replaced by this complete source.
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

    pub(crate) fn prepare_free_surface_step_with_acceleration<E>(
        self,
        surface: &mut FreeSurface3,
        nominal_acceleration: [f64; 3],
        acceleration_at: impl FnMut(u64, usize, f64) -> Result<[f64; 3], E>,
        gate: &CancelGate,
    ) -> Result<(), E>
    where
        E: From<FreeSurfaceError3>,
    {
        self.frame().prepare_free_surface_source::<true, E>(
            surface,
            nominal_acceleration,
            acceleration_at,
            gate,
        )
    }

    pub(crate) fn step_free_surface(
        self,
        surface: &mut FreeSurface3,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.frame()
            .step_free_surface_source::<true>(surface, acceleration, gate)
    }
}

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
    /// The material adapter's uniform gravity transfers to this complete source
    /// without a second gravity load. A separately selected native acceleration
    /// still refuses, including an owned zero or affine acceleration. This
    /// entry never combines independent acceleration owners.
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
        acceleration_at: impl FnMut(u64, usize, f64) -> Result<[f64; 3], E>,
        gate: &CancelGate,
    ) -> Result<(), E>
    where
        E: From<FreeSurfaceError3>,
    {
        self.prepare_free_surface_source::<false, E>(
            surface,
            nominal_acceleration,
            acceleration_at,
            gate,
        )
    }

    fn prepare_free_surface_source<const CORIOLIS_ONLY: bool, E>(
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
        if surface.body_acceleration_field().is_some() && !surface.owns_material_gravity() {
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
                let effective = self.source_acceleration::<CORIOLIS_ONLY>(
                    position,
                    velocity,
                    acceleration,
                    denominator,
                );
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
        // The complete frame source already includes gravity and uses actual
        // cell density. Publish its ownership only after the Guo reframe passes.
        surface.rebind_force_field_inner(forces, false, || check_cancel(gate))?;
        surface.body_acceleration = None;
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
        self.step_free_surface_source::<false>(surface, acceleration, gate)
    }

    fn step_free_surface_source<const CORIOLIS_ONLY: bool>(
        self,
        surface: &mut FreeSurface3,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let mut candidate = surface.clone();
        self.prepare_free_surface_source::<CORIOLIS_ONLY, FreeSurfaceError3>(
            &mut candidate,
            acceleration,
            |_, _, _| Ok(acceleration),
            gate,
        )?;
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
) -> Result<[f64; Q3], FreeSurfaceError3>{
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
    /// Uniform lattice body acceleration [dx/dt²], when selected. An affine
    /// owner returns None; use [`Self::body_acceleration_field`] to distinguish
    /// that owner from historical/reference-density or independently held force.
    /// Some([0; 3]) remains an acceleration-owned zero load, not unknown forcing.
    #[must_use]
    pub const fn body_acceleration(&self) -> Option<[f64; 3]> {
        match self.body_acceleration {
            Some(BodyAccelerationField3::Uniform(acceleration)) => Some(acceleration),
            _ => None,
        }
    }

    /// Complete accepted density-dependent acceleration owner, if selected.
    /// Ordinary timesteps hold this field; endpoint advances publish its next
    /// value only together with accepted populations, mass and force history.
    #[must_use]
    pub const fn body_acceleration_field(&self) -> Option<BodyAccelerationField3> {
        self.body_acceleration
    }

    pub(super) fn acceleration_at_lane(
        &self,
        key: u64,
        lane: usize,
    ) -> Result<Option<[f64; 3]>, FreeSurfaceError3> {
        self.body_acceleration
            .as_ref()
            .map(|field| field.at_cell(key, lane))
            .transpose()
    }

    /// Advance one lattice timestep to a prescribed body-acceleration endpoint.
    ///
    /// The accepted force history supplies `F_n = rho_n * a_n` to collision.
    /// Streaming then retains `F_next = rho_next * next_acceleration` in the
    /// physical velocity and gas-pressure source. These two endpoints give
    /// trapezoidal body impulse, including when density changes during the
    /// interval. There is no zero-time force reframe before or after the step:
    /// that would discard the physical half impulse of the changed load.
    /// Equal acceleration endpoints reproduce [`Self::step_with_gate`] exactly.
    ///
    /// First select the initial acceleration with
    /// [`Self::rebind_body_acceleration`]. This entry cannot infer an acceleration
    /// from an independently held force density. The new acceleration becomes
    /// the retained value only after the complete native step, conversion and
    /// cancellation checks succeed. Ordinary steps thereafter hold that value.
    /// Native low-Mach, capillary and boundary-resolution limits still apply;
    /// this does not evolve a vessel or add a mechanical-energy/heat account.
    pub fn advance_body_acceleration(
        &mut self,
        next_acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.advance_body_acceleration_field(
            BodyAccelerationField3::Uniform(next_acceleration),
            gate,
        )
    }

    /// Advance to a uniform or affine acceleration endpoint on the same native
    /// stream. Collision uses the previous local `rho*a(x)`; gas reconstruction,
    /// accepted velocity and newly wet cells use the next field at their actual
    /// cell centres and their own next densities. The descriptor remains fixed
    /// in world/grid coordinates, not attached to individual parcels.
    ///
    /// This shares the scalar entry's ownership, atomicity and numerical scope.
    /// It adds no interpolation, dense storage, per-cell material callback or
    /// extra population buffer. An affine law is not a velocity-dependent
    /// inertial force or a thermodynamic density/thermal-expansion model.
    pub fn advance_body_acceleration_field(
        &mut self,
        next: BodyAccelerationField3,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.advance_body_acceleration_with_pressure(next, gate, None)
    }

    /// The same endpoint-acceleration step with prescribed gas-link pressure.
    ///
    /// Pressure and its density-excursion bound have the same units, sampling
    /// and physical scope as [`Self::step_with_gas_link_pressure`]. Gas links
    /// close their even Guo moment using the actual next density and the NEXT
    /// acceleration, while collision still uses the retained OLD force. The
    /// pressure field describes this one frozen boundary step. Any pressure,
    /// flow, conversion or cancellation refusal preserves both force epochs
    /// and every accepted physical state.
    pub fn advance_body_acceleration_with_gas_link_pressure(
        &mut self,
        next_acceleration: [f64; 3],
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink3) -> Result<f64, FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        self.advance_body_acceleration_field_with_gas_link_pressure(
            BodyAccelerationField3::Uniform(next_acceleration),
            gate,
            maximum_density_excursion,
            pressure,
        )
    }

    /// Compose a spatial acceleration endpoint with the original gas-link load.
    /// Pressure assembly samples each interface's next acceleration once, then
    /// closes its gas links in that local actual-density force frame. A pressure or
    /// spatial-field refusal leaves the complete accepted trajectory intact.
    pub fn advance_body_acceleration_field_with_gas_link_pressure(
        &mut self,
        next: BodyAccelerationField3,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink3) -> Result<f64, FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let boundary = pressure::GasPressure3::new(
            maximum_density_excursion,
            pressure::GasPressureField3::Link(&pressure),
        )?;
        self.advance_body_acceleration_with_pressure(next, gate, Some(&boundary))
    }

    fn advance_body_acceleration_with_pressure(
        &mut self,
        next: BodyAccelerationField3,
        gate: &CancelGate,
        pressure: Option<&pressure::GasPressure3<'_>>,
    ) -> Result<(), FreeSurfaceError3> {
        let mut candidate = self.body_acceleration_endpoint_candidate(next, gate)?;
        candidate.step_inner(gate, None, None, pressure)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(())
    }

    // A native boundary owner can stage this SAME source epoch on its one
    // private candidate. The returned populations still carry F_n; callers
    // must run the native stream/conversion before publishing the candidate.
    pub(super) fn body_acceleration_endpoint_candidate(
        &self,
        next: BodyAccelerationField3,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceError3> {
        check_cancel(gate)?;
        self.admit_body_acceleration_mode()?;
        if self.body_acceleration.is_none() {
            return Err(FreeSurfaceError3::Transport(
                "acceleration evolution requires an explicit initial body-acceleration owner",
            ));
        }
        next.admit()?;
        let mut candidate = self.clone();
        // Do not rebind local_force: it is the accepted collision source F_n.
        // The existing stream already uses this descriptor only to close the
        // endpoint pressure, refresh F_next, and initialize newly wet cells.
        candidate.body_acceleration = Some(next);
        if candidate.material_gravity.is_some() {
            candidate.material_gravity = match next {
                BodyAccelerationField3::Uniform(acceleration) => Some(acceleration),
                BodyAccelerationField3::Affine { .. } => None,
            };
        }
        Ok(candidate)
    }

    // Collision admits the OLD source. The endpoint half source can make the
    // actual next physical velocity inadmissible even with positive, low-speed
    // raw populations. Run after conversion so new and retargeted wet cells
    // receive the same final admission, without altering their momentum.
    pub(super) fn admit_body_acceleration_state(
        &self,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        if self.body_acceleration.is_none() {
            return Ok(());
        }
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active acceleration tile");
            for (lane, cell) in tile.cells.iter().enumerate() {
                if *cell == Cell3::Gas {
                    continue;
                }
                if self
                    .grid
                    .populations(slot, lane)
                    .iter()
                    .any(|v| !v.is_finite() || *v <= 0.0)
                {
                    return Err(FreeSurfaceError3::Transport(
                        "body acceleration needs positive finite final wet populations",
                    ));
                }
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                let force = self.force_at_lane(key, lane);
                let velocity: [f64; 3] =
                    core::array::from_fn(|d| raw[d] + 0.5 * force[d] / rho);
                let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
                if !rho.is_finite()
                    || rho <= 0.0
                    || !speed.is_finite()
                    || speed / CS2_3.sqrt() >= crate::MACH_LIMIT
                {
                    return Err(FreeSurfaceError3::Transport(
                        "body-acceleration endpoint exceeds the physical lattice Mach domain",
                    ));
                }
            }
        }
        Ok(())
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

    /// Select or replace a uniform/affine field without advancing physical time.
    /// Reuses the original Guo reframe: density and mass remain exact, while
    /// physical velocity and stress are retained to population roundoff.
    /// Subsequent steps evaluate `rho*a(x)` at each actual cell position,
    /// including newly wet cells and cells relocated by a piston remap.
    /// Independent local-force and historical replay owners remain refused.
    pub fn rebind_body_acceleration_field(
        &mut self,
        field: BodyAccelerationField3,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.rebind_body_acceleration_field_checked(field, || check_cancel(gate))
    }

    /// Same candidate admission under an existing Cx/owner cancellation scope.
    pub(crate) fn rebind_body_acceleration_checked(
        &mut self,
        acceleration: [f64; 3],
        checkpoint: impl FnMut() -> Result<(), FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        self.rebind_body_acceleration_field_checked(
            BodyAccelerationField3::Uniform(acceleration),
            checkpoint,
        )
    }

    fn rebind_body_acceleration_field_checked(
        &mut self,
        field: BodyAccelerationField3,
        mut checkpoint: impl FnMut() -> Result<(), FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        checkpoint()?;
        field.admit()?;
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
                    body_force(populations.iter().sum(), field.at_cell(key, lane)?)?,
                );
            }
        }
        self.rebind_force_field_inner(forces, true, checkpoint)?;
        self.body_acceleration = Some(field);
        self.material_gravity = None;
        Ok(())
    }

    /// Select the material adapter's uniform gravity. Unlike an independently
    /// authored acceleration, this owner can hand off to a complete thermal or
    /// frame force law that already includes that same physical gravity.
    pub(crate) fn rebind_material_gravity(
        &mut self,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.rebind_material_gravity_checked(acceleration, || check_cancel(gate))
    }

    pub(crate) fn rebind_material_gravity_checked(
        &mut self,
        acceleration: [f64; 3],
        checkpoint: impl FnMut() -> Result<(), FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        self.rebind_body_acceleration_checked(acceleration, checkpoint)?;
        self.material_gravity = Some(acceleration);
        Ok(())
    }

    pub(crate) const fn material_gravity(&self) -> Option<[f64; 3]> {
        self.material_gravity
    }

    fn owns_material_gravity(&self) -> bool {
        self.material_gravity.is_some_and(|acceleration| {
            self.body_acceleration == Some(BodyAccelerationField3::Uniform(acceleration))
        })
    }

    /// Actual carrier-density gravity at one wet cell, without velocity/stress
    /// reconstruction. Coupled force maps use this same acceleration convention.
    pub(crate) fn material_gravity_force_at(
        &self,
        key: u64,
        lane: usize,
    ) -> Result<[f64; 3], FreeSurfaceError3> {
        let acceleration = self.material_gravity.ok_or(FreeSurfaceError3::Transport(
            "coupled material forcing requires its original uniform gravity",
        ))?;
        let tile = self.tiles.get(&key).ok_or(FreeSurfaceError3::Transport(
            "material gravity needs an active wet cell",
        ))?;
        if lane >= TILE_CELLS || tile.cells[lane] == Cell3::Gas {
            return Err(FreeSurfaceError3::Transport(
                "material gravity needs a wet cell",
            ));
        }
        let slot = self.grid.slot_of(key).expect("active material cell");
        body_force(self.grid.populations(slot, lane).iter().sum(), acceleration)
    }

    pub(crate) fn material_gravity_forces(
        &self,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), [f64; 3]>, FreeSurfaceError3> {
        check_cancel(gate)?;
        let mut forces = BTreeMap::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Gas {
                    forces.insert((key, lane), self.material_gravity_force_at(key, lane)?);
                }
            }
        }
        Ok(forces)
    }

    /// Bind a COMPLETE material force field on an enclosing private candidate.
    /// Its gravity component uses current rho, and additional physical forces
    /// stay force densities. A nonzero extra force selects the existing held
    /// force split. When extras vanish, restore native density-aware gravity so
    /// zero buoyancy/Marangoni/Lorentz force retains the original trajectory.
    /// No impulse, mass, heat or extra sparse force map is introduced.
    pub(crate) fn rebind_material_force_field(
        &mut self,
        forces: BTreeMap<(u64, usize), [f64; 3]>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let gravity = self.material_gravity.ok_or(FreeSurfaceError3::Transport(
            "complete material forces cannot replace independent acceleration",
        ))?;
        if self.body_acceleration.is_some() && !self.owns_material_gravity() {
            return Err(FreeSurfaceError3::Transport(
                "complete material forces cannot replace independent acceleration",
            ));
        }
        let mut gravity_only = true;
        for (index, (&(key, lane), &force)) in forces.iter().enumerate() {
            if index % TILE_CELLS == 0 {
                check_cancel(gate)?;
            }
            gravity_only &= force == self.material_gravity_force_at(key, lane)?;
        }
        // Original admission checks complete wet-cell coverage and reframes
        // atomically. Pure gravity keeps the required full endpoint history.
        self.rebind_force_field_inner(forces, gravity_only, || check_cancel(gate))?;
        self.body_acceleration = gravity_only.then_some(BodyAccelerationField3::Uniform(gravity));
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
