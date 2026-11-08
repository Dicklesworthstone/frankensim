//! Uniform physical body acceleration through the retained full reference mass.
//!
//! A constant translation field has zero director acceleration. Applying the
//! original augmented mass to that field includes graded shell weight and every
//! payload's COM force and offset couple. This is prescribed loading in the fixed
//! Cartesian frame, not self-gravity, an accelerating frame or follower gravity.

use super::{
    MovingLoadError, MovingShellModel, MovingShellMotion, MovingShellStep, ShellAssembly,
    ShellLoadProjection, ShellPointWrench, ShellSupportMotion, ShellTimeConfig, finite, poll,
};

impl MovingShellModel {
    // Three rigid-translation columns of the SAME full assembled inertia.
    // Preparation is bounded by the existing dense shell matrix; endpoints only
    // visit these three coefficients per row. Payload attachment refreshes this
    // once after all payloads have augmented the full mass.
    pub(super) fn prepare_body_acceleration_map(
        assembly: &ShellAssembly,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MovingLoadError> {
        let n = assembly.full_mass.dimension();
        let mut map = Vec::with_capacity(n);
        for row in assembly.full_mass.values().chunks_exact(n) {
            poll(cancelled)?;
            let mut value = [0.0; 3];
            for node in row.chunks_exact(6) {
                for axis in 0..3 {
                    value[axis] += node[axis];
                }
            }
            finite(&value)?;
            map.push(value);
        }
        Ok(map)
    }

    /// Superpose point wrenches with uniform physical body acceleration [m/s²].
    ///
    /// The acceleration acts on the entire retained reference inventory,
    /// including attached payloads, in the model's fixed Cartesian frame. The
    /// full mass maps the repeated translational acceleration and zero director
    /// acceleration to nodal forces and couples before support restriction.
    /// Offset payload moments and restrained-node weight therefore survive.
    ///
    /// The three mass columns are prepared once; endpoint work is linear in the
    /// number of DOFs, with no density resampling or extra point-load entries.
    /// Do not also supply the same shell or payload weight through point ports.
    /// Zero acceleration retains the original point projection exactly.
    pub fn project_with_body_acceleration(
        &self,
        loads: &[ShellPointWrench],
        acceleration_m_s2: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellLoadProjection, MovingLoadError> {
        poll(&mut cancelled)?;
        finite(&acceleration_m_s2)?;
        let mut load = self.project(loads, &mut cancelled)?;
        if acceleration_m_s2 == [0.0; 3] {
            return Ok(load);
        }
        for (value, row) in load.full.iter_mut().zip(&self.body_acceleration_map) {
            poll(&mut cancelled)?;
            *value += row
                .iter()
                .zip(acceleration_m_s2)
                .map(|(m, a)| m * a)
                .sum::<f64>();
        }
        finite(&load.full)?;
        for (value, &row) in load.reduced.iter_mut().zip(&self.assembly.free_dofs) {
            *value = load.full[row];
        }
        load.resultant = self.resultant(&load.full)?;
        poll(&mut cancelled)?;
        Ok(load)
    }

    /// Start native motion with explicit initial body acceleration and wrenches.
    ///
    /// Initial equilibrium, reactions and all later impulse events retain the
    /// combined full load. Supply the next endpoint acceleration explicitly to
    /// [`MovingShellMotion::advance_with_body_acceleration`]; no gravity is
    /// inferred or automatically carried to a newly supplied endpoint.
    pub fn start_with_body_acceleration<'a>(
        &'a self,
        displacement: &[f64],
        velocity: &[f64],
        initial_loads: &[ShellPointWrench],
        acceleration_m_s2: [f64; 3],
        config: ShellTimeConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellMotion<'a>, MovingLoadError> {
        self.start_with_supports(
            displacement,
            velocity,
            initial_loads,
            config,
            None,
            Some(acceleration_m_s2),
            cancelled,
        )
    }

    /// Start prescribed-pin motion with explicit body acceleration and wrenches.
    ///
    /// `supports` retains the original support-node order. Body loading and
    /// prescribed support acceleration are independent physical inputs; the
    /// original full-mass condensation includes both before solving free rows.
    /// Use [`MovingShellMotion::advance_driven_with_body_acceleration`] for the
    /// following endpoints. No support actuator or gravity is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn start_driven_with_body_acceleration<'a>(
        &'a self,
        displacement: &[f64],
        velocity: &[f64],
        supports: ShellSupportMotion,
        initial_loads: &[ShellPointWrench],
        acceleration_m_s2: [f64; 3],
        config: ShellTimeConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellMotion<'a>, MovingLoadError> {
        self.start_with_supports(
            displacement,
            velocity,
            initial_loads,
            config,
            Some(supports),
            Some(acceleration_m_s2),
            cancelled,
        )
    }
}

impl MovingShellMotion<'_> {
    /// Advance with explicit endpoint body acceleration and point wrenches.
    ///
    /// The old combined endpoint force remains held through impulse events and
    /// enters the same trapezoidal impulse/work account as the new force. The
    /// original Newmark and mass factors are reused. Resolve time-varying fields
    /// through timestep refinement; this is uniform prescribed acceleration,
    /// without a spatially varying field or configuration-dependent gravity.
    ///
    /// A failed projection, solve or cancellation leaves the accepted state,
    /// held force and reactions unchanged. Calling the original `advance`
    /// instead supplies zero body acceleration at that new endpoint.
    pub fn advance_with_body_acceleration(
        &mut self,
        time_s: f64,
        loads: &[ShellPointWrench],
        acceleration_m_s2: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellStep, MovingLoadError> {
        self.advance_with_supports(time_s, loads, None, Some(acceleration_m_s2), cancelled)
    }

    /// Advance driven pins with explicit body and support endpoint acceleration.
    ///
    /// The two accelerations have their own physical meanings. Support work and
    /// body-load work retain the original separate accounting; neither is
    /// counted twice. All timestamp, horizon and atomic-publication rules of
    /// [`Self::advance_driven`] apply to the combined load.
    pub fn advance_driven_with_body_acceleration(
        &mut self,
        time_s: f64,
        loads: &[ShellPointWrench],
        body_acceleration_m_s2: [f64; 3],
        support_acceleration_m_s2: [[f64; 3]; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellStep, MovingLoadError> {
        self.advance_with_supports(
            time_s,
            loads,
            Some(support_acceleration_m_s2),
            Some(body_acceleration_m_s2),
            cancelled,
        )
    }
}
