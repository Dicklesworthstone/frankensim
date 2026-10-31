//! Source-backed ordinary solids with implicit interior-fluid momentum correction.
//!
//! This one-body lane retains the ORIGINAL dry mass and complete inertia at
//! every density ratio. The native resolved flow corrects its reaction using
//! actual streamed fluid momentum inside the sampled moving geometry. No
//! displaced tensor is subtracted, no epsilon mass is introduced, and neutral
//! buoyancy does not suppress the body's original momentum or fluid feedback.
//!
//! A fixed reference fluid uses reduced pressure: (M_s-rho_f V)g is applied
//! once. The original material volume supplies V, independently of quadrature.
//! The fluid is initially quiescent, single phase and incompressible in the
//! native low-Mach sense. No free surface, thermal expansion, near-contact film
//! or exact-rigid-interior guarantee is implied. Inspect native sampled-volume,
//! actual interior-mass/slip and kinetic diagnostics and refine dx/dt/quadrature.

use fs_lbm::d3q19::{DryBodyCouplingConfig3, DryBodyExchange3};

use super::*;

/// ONE ordinary sphere/ellipsoid on the existing source-scaled D3Q19 owner.
/// Its private native container is reused for units/source identity only; this
/// lane always dispatches the explicit dry-body implicit correction.
pub struct ResolvedImmersedBody<'a, F: FluidTransport + ?Sized> {
    native: RetainedFluidBodies<'a, F, 1, &'a ResolvedBody>,
    volume_m3: f64,
    displaced_mass_kg: f64,
    reduced_gravity_force_n: [f64; 3],
    gravity_m_s2: [f64; 3],
}

impl<'a, F: FluidTransport + ?Sized> ResolvedImmersedBody<'a, F> {
    /// Bind one original uniform spherical solid, of ANY positive density, to
    /// its matching native surface. Initial body motion is retained exactly;
    /// the surrounding/retained fluid is initially quiescent.
    ///
    /// The existing geometry owner checks the original region and full inertia.
    /// Native admission/iteration occurs on each step; construction does not
    /// promise convergence at an arbitrary density ratio or resolution.
    #[allow(clippy::too_many_arguments)]
    pub fn sphere(
        fluid: &'a F,
        domain: ImmersedDomain,
        material: &'a ResolvedBody,
        radius_m: f64,
        initial: BodyStateSi,
        gravity_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let source = AdditionalMaterialSphere::new(material, radius_m)?;
        Self::with_shape(
            fluid,
            domain,
            material,
            initial,
            gravity_m_s2,
            VolumeMoments::sphere(source.radius_m())?.volume_m3(),
            |scaling| NativeEllipsoid3::sphere(source.radius_m() / scaling.cell_size_m()),
            gate,
        )
    }

    /// Bind one original uniform ellipsoid, including its complete rotated
    /// body-frame tensor. Shape and material axes use the same original rotation.
    #[allow(clippy::too_many_arguments)]
    pub fn ellipsoid(
        fluid: &'a F,
        domain: ImmersedDomain,
        material: &'a ResolvedBody,
        semiaxes_m: [f64; 3],
        body_from_shape: So3,
        initial: BodyStateSi,
        gravity_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let source = AdditionalMaterialEllipsoid::new(material, semiaxes_m, body_from_shape)?;
        Self::with_shape(
            fluid,
            domain,
            material,
            initial,
            gravity_m_s2,
            VolumeMoments::ellipsoid(vector(semiaxes_m))?.volume_m3(),
            |scaling| source.lattice_shape(scaling.cell_size_m()),
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn with_shape(
        fluid: &'a F,
        domain: ImmersedDomain,
        material: &'a ResolvedBody,
        initial: BodyStateSi,
        gravity_m_s2: [f64; 3],
        volume_m3: f64,
        shape: impl FnOnce(&MaterialScaling<'a, F>) -> Result<NativeEllipsoid3, Error>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let volume_m3 = positive(volume_m3)?;
        let gravity_m_s2 = convert(gravity_m_s2, 1.0)?;
        // The source-derived dry spatial inertia passes unchanged through the
        // original material scaling/construction path. No density matching or
        // excess-mass branch exists here.
        let native = RetainedFluidBodies::new_with_surfaces(
            fluid,
            domain,
            [material],
            [material],
            [initial],
            |scaling| Ok(CarrierSurfaces::Ellipsoids([shape(scaling)?])),
            gate,
        )?;
        let displaced_mass_kg = positive(native.scaling.transport().density_kg_m3() * volume_m3)?;
        let reduced_mass = finite(material.inertia().mass() - displaced_mass_kg)?;
        let reduced_gravity_force_n = convert(gravity_m_s2, reduced_mass)?;
        poll(gate)?;
        Ok(Self {
            native,
            volume_m3,
            displaced_mass_kg,
            reduced_gravity_force_n,
            gravity_m_s2,
        })
    }

    /// Original body, region bindings and density/query provenance.
    pub const fn material(&self) -> &'a ResolvedBody {
        self.native.sources[0]
    }

    /// Original fluid snapshot and source-owned lattice scaling.
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        self.native.scaling()
    }

    /// Analytic ORIGINAL material volume; never replaced by sampled occupancy.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// Reference fluid density times original material volume.
    pub const fn displaced_mass_kg(&self) -> f64 {
        self.displaced_mass_kg
    }

    /// Held physical gravity. No gravity source is also installed in the grid.
    pub const fn gravity_m_s2(&self) -> [f64; 3] {
        self.gravity_m_s2
    }

    /// (M_dry-rho_f V)g [N], including exact neutral cancellation when the
    /// original source values agree. This changes LOAD only, never inertia.
    pub const fn reduced_gravity_force_n(&self) -> [f64; 3] {
        self.reduced_gravity_force_n
    }

    /// Accepted authoritative fluid; independent state/force writes are excluded.
    pub const fn grid(&self) -> &BoundaryGrid3 {
        self.native.grid()
    }

    /// Original DRY inertia and accepted kinematics, in lattice units.
    pub const fn body_lattice(&self) -> &MovingPartialSaturationBody3 {
        &self.native.bodies[0]
    }

    /// Actual accepted pose and original DRY spin momentum in SI units.
    /// A conversion refusal is observational and never rolls back a prior step.
    pub fn state_si(&self) -> Result<BodyStateSi, Error> {
        self.native.state_si(0)
    }

    /// Number of complete accepted physical intervals.
    pub const fn accepted_steps(&self) -> u64 {
        self.native.accepted_steps()
    }

    /// Reuse the original viscosity update and source identity checks. Changing
    /// reference density is refused: it would alter actual fluid mass, units
    /// and reduced-pressure gravity without conservative transfer.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.native.rebind_fluid(fluid, gate)
    }

    /// Advance authored WORLD forces [N] and center couples [N m], with reduced
    /// gravity added once. Do not include weight or buoyancy again in authored.
    ///
    /// The native solver preserves positive dry inertia, solves actual streamed
    /// interior momentum feedback, bounds sampled-volume changes, and admits
    /// independent linear/angular balances before publishing. Its original
    /// complete receipt is returned IN LATTICE UNITS. Use the explicit unit
    /// converters below; these do not relabel correction work as dissipation.
    ///
    /// Refusal preserves fluid, dry body and step count. No fallible conversion
    /// or diagnostic follows the native atomic publication.
    pub fn step_forces_si(
        &mut self,
        authored: [f64; 6],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodyExchange3, Error> {
        poll(gate)?;
        let next = self
            .native
            .accepted_steps
            .checked_add(1)
            .ok_or(Error::Invalid("step count overflowed"))?;
        let mut loads = authored;
        for (axis, load) in loads.iter_mut().enumerate().take(3) {
            *load = finite(*load + self.reduced_gravity_force_n[axis])?;
        }
        let dt = self.native.scaling.timestep_s();
        let mut impulse = [0.0; 6];
        for axis in 0..6 {
            let unit = if axis < 3 {
                self.native.impulse_unit
            } else {
                self.native.angular_impulse_unit
            };
            impulse[axis] = finite(loads[axis] * (dt / unit))?;
        }
        let origin = convert(moment_origin_m, 1.0 / self.native.scaling.cell_size_m())?;
        let CarrierSurfaces::Ellipsoids(shapes) = &self.native.surfaces else {
            return Err(Error::Invalid(
                "dry-body source requires its original resolved shape",
            ));
        };
        let receipt = self.native.grid.step_with_dry_partial_saturation_body(
            &mut self.native.samplers[0],
            &shapes[0],
            &mut self.native.bodies[0],
            impulse,
            origin,
            controls,
            gate,
        )?;
        self.native.accepted_steps = next;
        Ok(receipt)
    }

    /// Convert native population mass, including the actual interior inventory.
    pub fn mass_kg(&self, lattice: f64) -> Result<f64, Error> {
        self.native.mass_kg(lattice)
    }

    /// Convert signed native linear impulse/momentum to N s.
    pub fn impulse_ns(&self, lattice: f64) -> Result<f64, Error> {
        self.native.impulse_ns(lattice)
    }

    /// Convert signed native angular impulse/momentum to N m s.
    pub fn angular_impulse_nms(&self, lattice: f64) -> Result<f64, Error> {
        self.native.angular_impulse_nms(lattice)
    }

    /// Convert a signed native work/kinetic account to J without relabeling it.
    pub fn energy_j(&self, lattice: f64) -> Result<f64, Error> {
        self.native.energy_j(lattice)
    }

    /// Reduced gravitational potential change along this accepted midpoint
    /// translation. This observational conversion can fail without mutating any
    /// state. The source/body volume is fixed; no free-surface potential is used.
    pub fn gravity_potential_change_j(&self, receipt: &DryBodyExchange3) -> Result<f64, Error> {
        let travel = convert(
            receipt.native.coupling.velocity_midpoint,
            self.native.scaling.cell_size_m(),
        )?;
        finite(
            -(0..3)
                .map(|i| self.reduced_gravity_force_n[i] * travel[i])
                .sum::<f64>(),
        )
    }
}
