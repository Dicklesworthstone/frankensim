//! Axial laminar flow between two fixed concentric cylindrical walls.
//!
//! Annular cell areas carry pressure, inertia and through-flow. Each shear
//! conductance is the integral cylindrical resistance between neighbouring
//! velocity centres (or a centre and the actual wall):
//! g = 2*pi*mu*L / ln(r_right/r_left). Thus steady zero-pressure axial Couette
//! flow retains its logarithmic profile, even near a small inner cylinder.
//! The original shared tridiagonal momentum owner supplies time advancement,
//! signed pressure/wall work, reaction impulse and nonnegative viscous heat.
//!
//! Geometry is evaluated with the radial gap, avoiding subtraction of squared
//! radii. A compensated logarithm retains narrow-gap shear when adding the
//! relative gap to one would round away. Both are cold assembly operations;
//! the existing factorization and trial field are reused during every step.
//!
//! The model is incompressible, fully developed and Newtonian. Both wall speeds
//! are AXIAL; cylinder rotation, normal motion, eccentricity, entrance effects,
//! thermal evolution and transition are not supplied. A zero inner radius is
//! a CircularChannel, whose axis is symmetry rather than a material wall.
//!
//! Continuum references: axial annular Couette and Poiseuille solutions of
//! rho*u_t = dp/L + mu/r * (r*u_r)_r with prescribed axial wall velocities.
//! <https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2020/08/Viscous-Laminar-Flows-Lesson-2-Internal-Flows-Handout-1.pdf>

use fs_math::det;

use super::{
    AcousticRealizeError, ChannelGasUpdate, GasState, MAX_TRANSVERSE_CELLS, MomentumProfile,
    SlitChannelLoad, SlitChannelStep, SlitExit, invalid,
};

/// Concentric cylindrical walls, their axial length, and independent resolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnnularChannelSpec {
    /// Distance between supplied pressure planes [m].
    pub length_m: f64,
    /// Positive radius of the inner material wall [m].
    pub inner_radius_m: f64,
    /// Outer wall radius, strictly greater than the inner radius [m].
    pub outer_radius_m: f64,
    /// Uniform cells across the radial gap, in 2..=4096.
    pub radial_cells: usize,
    /// Backward-Euler period [s], independently refinable.
    pub sample_period_s: f64,
    /// Original static-pressure or explicit uniform kinetic-head outlet.
    pub exit: SlitExit,
}

/// Original material gas and complete annular momentum field.
#[derive(Debug, Clone)]
pub struct AnnularChannel {
    gas: GasState,
    spec: AnnularChannelSpec,
    profile: MomentumProfile,
}

impl AnnularChannel {
    /// Build a quiescent annulus from one coherent gas; no transport overrides.
    pub fn from_gas(
        gas: &GasState,
        spec: AnnularChannelSpec,
    ) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence()
            .map_err(|_| invalid("annular channel needs one coherent material gas"))?;
        let profile = Self::profile_from_transport(gas.density, gas.dynamic_viscosity, spec)?;
        Ok(Self {
            gas: *gas,
            spec,
            profile,
        })
    }

    /// One geometry assembly for gas and generic material-transport owners.
    pub(super) fn profile_from_transport(
        density: f64,
        viscosity: f64,
        spec: AnnularChannelSpec,
    ) -> Result<MomentumProfile, AcousticRealizeError> {
        if ![
            density,
            viscosity,
            spec.length_m,
            spec.inner_radius_m,
            spec.outer_radius_m,
            spec.sample_period_s,
        ]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
            || spec.outer_radius_m <= spec.inner_radius_m
            || !(2..=MAX_TRANSVERSE_CELLS).contains(&spec.radial_cells)
        {
            return Err(invalid(
                "annular channel needs positive transport, length, time, ordered radii and 2..=4096 cells",
            ));
        }
        let n = spec.radial_cells;
        let gap = spec.outer_radius_m - spec.inner_radius_m;
        let dr = gap / n as f64;
        let shear = core::f64::consts::TAU * viscosity * spec.length_m;
        if ![dr, shear].iter().all(|x| x.is_finite() && *x > 0.0) {
            return Err(invalid(
                "annular channel coefficients cannot be represented",
            ));
        }
        // pi*(r_right^2-r_left^2), evaluated without cancellation in a thin gap.
        let ring_width = core::f64::consts::TAU * dr;
        let areas = (0..n)
            .map(|i| ring_width * (spec.inner_radius_m + (i as f64 + 0.5) * dr))
            .collect();
        let mut conductances = Vec::with_capacity(n + 1);
        for face in 0..=n {
            let (left, distance) = if face == 0 {
                (spec.inner_radius_m, 0.5 * dr)
            } else if face == n {
                (spec.inner_radius_m + (n as f64 - 0.5) * dr, 0.5 * dr)
            } else {
                (spec.inner_radius_m + (face as f64 - 0.5) * dr, dr)
            };
            conductances.push(shear / log_ratio_from_gap(left, distance));
        }
        // Both boundaries are material walls; a rounded-away inner conductance
        // must not silently become the circular owner's symmetry condition.
        if conductances.iter().any(|g| !g.is_finite() || *g <= 0.0) {
            return Err(invalid("annular wall shear cannot be represented"));
        }
        MomentumProfile::new(
            density,
            spec.length_m,
            spec.sample_period_s,
            spec.exit,
            areas,
            conductances,
        )
    }

    /// Original continuous pressure/flow field, with every annular cell retained.
    /// The existing nonlinear outlet refuses this linear harmonic view.
    pub fn continuous_operator(
        &self,
    ) -> Result<super::harmonic::ChannelContinuousOperator<'_>, AcousticRealizeError> {
        self.profile.continuous_operator()
    }

    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    #[must_use]
    pub const fn spec(&self) -> &AnnularChannelSpec {
        &self.spec
    }
    /// Accepted velocities at radial centres, inner to outer [m/s].
    #[must_use]
    pub fn velocities_m_s(&self) -> &[f64] {
        self.profile.velocities()
    }
    /// Historical exchange; a zero-duration material change does not rewrite it.
    #[must_use]
    pub const fn last_step(&self) -> SlitChannelStep {
        self.profile.last_step()
    }
    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.profile.stored_energy_j()
    }

    /// Replace the complete gas at fixed geometry while preserving each cell's
    /// axial momentum. The returned signed kinetic parameter work belongs to the
    /// prescribing material owner; it is neither wall heat nor pressure work.
    /// Refusal or cancellation leaves the original source and field unchanged.
    pub fn rebind_gas_controlled(
        &mut self,
        gas: &GasState,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelGasUpdate, AcousticRealizeError> {
        super::gas::admit(gas, maximum_absolute_parameter_work_j, cancelled)?;
        if gas == &self.gas {
            super::gas::poll(cancelled)?;
            return Ok(ChannelGasUpdate::unchanged(self.stored_energy_j()));
        }
        let mut candidate = Self::from_gas(gas, self.spec)?;
        let update = candidate.profile.transfer_material_history(
            &self.profile,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        super::gas::poll(cancelled)?;
        *self = candidate;
        Ok(update)
    }

    /// Next-step affine pressure port with stationary walls.
    pub fn pressure_load(&self) -> Result<SlitChannelLoad, AcousticRealizeError> {
        self.profile.pressure_load()
    }
    /// Same pressure port with held axial wall speeds in [inner, outer] order.
    /// Wall forcing is in the offset; accepted wall work remains separate.
    pub fn pressure_load_with_wall_velocities(
        &self,
        walls_m_s: [f64; 2],
    ) -> Result<SlitChannelLoad, AcousticRealizeError> {
        self.profile.pressure_load_with_walls(walls_m_s)
    }
    /// Added-mass part rho*L/(area*dt) of the next-step pressure impedance.
    #[must_use]
    pub fn inertial_impedance_pa_s_m3(&self) -> f64 {
        self.profile.inertial_impedance()
    }
    /// Unpublished candidate at the requested signed endpoint through-flow.
    pub fn preview_at_flow(&self, flow: f64) -> Result<Self, AcousticRealizeError> {
        self.preview_at_flow_with_wall_velocities(flow, [0.0; 2])
    }
    /// Unpublished candidate with the same wall-driven pressure condensation.
    /// A coupled owner must accept both pressure work and wall work.
    pub fn preview_at_flow_with_wall_velocities(
        &self,
        flow: f64,
        walls_m_s: [f64; 2],
    ) -> Result<Self, AcousticRealizeError> {
        let pressure = self
            .pressure_load_with_wall_velocities(walls_m_s)?
            .pressure_drop_pa(flow)?;
        let mut trial = self.clone();
        trial.step_with_wall_velocities(pressure, walls_m_s)?;
        Ok(trial)
    }
    /// Stationary-wall steady response of this annular discretization.
    pub fn steady_bulk_speed_m_s(&self, pressure: f64) -> Result<f64, AcousticRealizeError> {
        self.profile.steady_bulk_speed(pressure)
    }
    /// Advance pressure-driven motion with stationary walls.
    pub fn step(&mut self, pressure: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        self.profile.step(pressure)
    }
    /// Advance with independent axial [inner, outer] wall speeds.
    /// Both walls have physical traction. The shared step retains signed work,
    /// impulses on the fluid, relative-shear heat and kinetic/numerical accounts.
    /// There is no allocation or factorization in this operation.
    pub fn step_with_wall_velocities(
        &mut self,
        pressure: f64,
        walls_m_s: [f64; 2],
    ) -> Result<SlitChannelStep, AcousticRealizeError> {
        self.profile.step_with_walls(pressure, walls_m_s)
    }
}

// ln((radius+distance)/radius) without losing the small geometric increment.
// The compensation corrects the rounded addition before calling the existing
// deterministic logarithm. A very wide ratio uses a difference of logarithms
// only when the ratio itself overflows; the result remains bounded and positive.
fn log_ratio_from_gap(radius: f64, distance: f64) -> f64 {
    let relative = distance / radius;
    if relative.is_finite() {
        let sum = 1.0 + relative;
        if sum == 1.0 {
            relative
        } else if relative < 0.5 {
            det::ln(sum) * (relative / (sum - 1.0))
        } else {
            det::ln(sum)
        }
    } else {
        det::ln(radius + distance) - det::ln(radius)
    }
}
