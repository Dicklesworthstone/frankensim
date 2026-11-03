//! SI initial flow and held fluid loads on the original material-scaled owner.

use super::*;

impl<'a, F: FluidTransport + ?Sized, const N: usize, S> RetainedFluidBodies<'a, F, N, S> {
    /// Prepare unit-reference-density fluid with a uniform physical velocity [m/s].
    ///
    /// Available only before the first accepted fluid/body interval. This is an
    /// initial condition, not an impulse applied to an evolved fluid. Body pose,
    /// velocity and momentum remain their authored initial values. A mismatch
    /// with body motion or stationary walls creates a physical startup transient;
    /// it is not a solved steady flow or a no-slip interior prescription.
    ///
    /// The source-owned scaling admits the speed. The native initializer uses
    /// the current Guo force frame, so a previously selected fluid force cannot
    /// add an unwanted half-force velocity. Cancellation changes no accepted state.
    pub fn initialize_uniform_flow_si(
        &mut self,
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        poll(gate)?;
        if self.accepted_steps != 0 {
            return Err(Error::Invalid(
                "fluid initialization requires the initial accepted time",
            ));
        }
        let velocity_m_s = convert(velocity_m_s, 1.0)?;
        let speed = velocity_m_s[0]
            .hypot(velocity_m_s[1])
            .hypot(velocity_m_s[2]);
        if speed > self.maximum_speed_m_s {
            return Err(Error::Invalid(
                "initial flow exceeds the source-scaled scenario speed bound",
            ));
        }
        let velocity = convert(
            velocity_m_s,
            self.scaling.timestep_s() / self.scaling.cell_size_m(),
        )?;
        self.grid
            .initialize_uniform_equilibrium(1.0, velocity, gate)
            .map_err(Error::Flow)
    }

    /// Hold new velocity [m/s] and pressure-perturbation [Pa] reservoir targets.
    ///
    /// Faces follow `Face3::ALL`. Their kinds and all wall/periodic data stay
    /// fixed. Conversion reuses the original source density, dx/dt and scenario
    /// speed bound, including after an admitted viscosity-only source rebind.
    /// All six targets are admitted before the native owner publishes any.
    ///
    /// This advances no time and changes no fluid population, force frame,
    /// body state, geometry or accepted-step count. The next accepted native
    /// interval applies the held targets; its shared `open` receipt owns actual
    /// reservoir mass and linear/angular momentum exchange exactly once.
    /// The caller samples time-varying controls at native tick boundaries.
    /// No waveform interpolation, finite-reservoir or pressure-work model is
    /// implied. Invalid targets or cancellation retain the previous controls.
    pub fn rebind_open_boundaries_si(
        &mut self,
        faces: [FaceBoundarySi; 6],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        poll(gate)?;
        let boundaries =
            boundary::lattice_boundaries(faces, &self.scaling, self.maximum_speed_m_s)?;
        self.grid
            .rebind_open_boundaries(boundaries, gate)
            .map_err(Error::Flow)
    }

    /// Hold a spatial inlet velocity profile in the fixed world/grid frame.
    ///
    /// The callback receives physical cell centers [m] in deterministic z/y/x
    /// order and returns velocity [m/s] or the original profile source error.
    /// It runs once per active fluid face node during this call. The caller owns
    /// the field's time and samples new traces at native tick boundaries.
    ///
    /// Native mean-normal control and reconstruction are reused. The FINAL
    /// corrected trace must satisfy the source-scaled scenario speed bound;
    /// samples and corrected values also retain the original low-Mach bound.
    /// Two plane buffers are reused, without a volume-sized converted field.
    ///
    /// This advances no populations, body state, force frame or accepted time.
    /// The following interval's shared `open` receipt owns actual reservoir
    /// mass and linear/angular momentum. Refusal or cancellation retains the
    /// entire previous trace. The returned diagnostics are IN LATTICE UNITS,
    /// with explicit `_lu` field names; no fallible conversion follows publication.
    pub fn set_velocity_profile_si(
        &mut self,
        face: Face3,
        flux: InletFluxSi,
        gate: &CancelGate,
        sample: impl FnMut([f64; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, Error> {
        boundary::set_velocity_profile_si(
            &mut self.grid,
            &self.scaling,
            self.maximum_speed_m_s,
            face,
            flux,
            gate,
            sample,
        )
    }

    /// Clear a held spatial trace and restore its current uniform face target.
    /// Cancellation preserves the profile. Populations and body history change
    /// only when a later native interval applies that target.
    pub fn clear_velocity_profile(&mut self, face: Face3, gate: &CancelGate) -> Result<(), Error> {
        if gate.is_requested() {
            return Err(Error::Profile(BoundaryProfileError3::Cancelled));
        }
        self.grid.clear_velocity_profile(face);
        Ok(())
    }

    /// Replace an externally authored uniform fluid volume force [N/m^3].
    ///
    /// This is a force on represented fluid only. The load is held for
    /// following native intervals; a rebind itself advances no time and applies
    /// no impulse or work. Original Guo rebinding preserves actual fluid mass,
    /// momentum and force-corrected stress. Body states, material units, viscosity
    /// and accepted clocks stay intact. An unchanged force preserves populations.
    ///
    /// Conversion is F_lattice = f_SI dt^2/(rho_reference dx), using the same
    /// source density as all fluid/body momentum conversions. Native dry-body
    /// correction already removes the force on fictitious interior fluid. Do not
    /// apply that subtraction or the fluid impulse again as an authored body load.
    /// The reference hydrostatic pressure is already removed by dry-body reduced
    /// gravity; do not add that gravity a second time through this external port.
    ///
    /// A reduced affine-pressure scenario needs a distinct solid traction too:
    /// for `p = p_periodic - f.x`, author `V_original * f` as the dry body's
    /// center force as well as supplying `f` here. Periodic pressure perturbations
    /// do not recover that affine traction. This physical solid load is separate
    /// from the fictitious-fluid subtraction and is not a duplicate fluid deposit.
    ///
    /// This supplies a piecewise-held forcing history. It does not implement an
    /// inlet/outlet pressure boundary, variable-density flow, or a closed finite
    /// reservoir. The native step still enforces its actual Mach and balance gates.
    pub fn rebind_uniform_fluid_force_si(
        &mut self,
        force_density_n_m3: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        poll(gate)?;
        let force = convert(force_density_n_m3, self.fluid_force_scale()?)?;
        self.grid
            .rebind_uniform_force(force, gate)
            .map_err(Error::Flow)
    }

    /// Replace a spatial fluid volume-force field [N/m^3] on the existing owner.
    ///
    /// The dense field covers the complete grid in x-fastest order:
    /// `(z * ny + y) * nx + x`. Every value must convert to a finite lattice
    /// force. Static-solid slots are ignored by the native sampler; partially
    /// occupied body cells remain fluid cells. The field is ADDED to the held
    /// uniform force; a zero field switches off only the spatial contribution.
    ///
    /// This prepares a piecewise-held load for following simultaneous fluid/body
    /// steps. The Guo transform preserves actual fluid density, momentum and full
    /// force-corrected stress without advancing time or changing body states.
    /// Dry-body owners retain their existing midpoint interior-force and torque
    /// subtraction. Authored solid loads are separate; this field does not infer
    /// a solid pressure traction or a hydrostatic buoyancy correction.
    ///
    /// Cancellation, a malformed field, conversion failure or native refusal
    /// leaves accepted populations, forces, body history and clocks intact.
    /// Sampling allocates only the native force proposal and reuses population
    /// scratch; no second lattice-sized converted field is constructed.
    pub fn rebind_fluid_force_field_si(
        &mut self,
        force_density_n_m3: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        poll(gate)?;
        let [nx, ny, nz] = self.grid.dimensions();
        if force_density_n_m3.len() != nx * ny * nz {
            return Err(Error::Invalid(
                "fluid force field must cover the whole grid",
            ));
        }
        let scale = self.fluid_force_scale()?;
        for chunk in force_density_n_m3.chunks(64) {
            poll(gate)?;
            for &force in chunk {
                convert(force, scale)?;
            }
        }
        self.grid
            .rebind_force_field_with_gate(
                |[x, y, z]| force_density_n_m3[(z * ny + y) * nx + x].map(|f| f * scale),
                gate,
            )
            .map_err(Error::Flow)
    }

    /// Observe the held external fluid volume force [N/m^3].
    pub fn uniform_fluid_force_density_n_m3(&self) -> Result<[f64; 3], Error> {
        convert(
            self.grid.uniform_force_density(),
            1.0 / self.fluid_force_scale()?,
        )
    }

    fn fluid_force_scale(&self) -> Result<f64, Error> {
        let dt = self.scaling.timestep_s();
        positive(
            (dt / self.scaling.cell_size_m()) * (dt / self.scaling.transport().density_kg_m3()),
        )
    }
}
