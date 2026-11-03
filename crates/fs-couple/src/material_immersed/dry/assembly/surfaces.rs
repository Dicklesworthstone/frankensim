//! Material inertia on the shared moving sphere/surface film solve.

use super::*;
use fs_lbm::d3q19::{
    DryLubricatedSurfacesExchange3, MovingFiniteLubricationPlane3, MovingLubricationSurface3,
};

mod shared;
pub use shared::{MaterialFilmCarrierSi, MaterialFilmPatchSi, ResolvedImmersedSharedSurfaces};

/// A prescribed reservoir or an independent material-backed finite carrier.
///
/// These are local planar missing-film models, not resolved no-slip boundaries.
/// The caller supplies the matching resolved boundary and the applicability of
/// the planar approximation. Source material supplies the finite carrier's full
/// inertia; an infinite plane does not infer its mass or validate a finite face.
/// Equal material references in different slots describe independent carriers.
#[derive(Clone, Copy, Debug)]
pub enum MaterialFilmSurfaceSi<'a> {
    /// Prescribed translation and axial spin, with a separate actuator-work port.
    Prescribed {
        /// Plane geometry, translation and explicit film scales.
        plane: PlaneFilmSi,
        /// Signed spin about the plane normal through its translating point [rad/s].
        axial_spin_rad_s: f64,
    },
    /// A freely translating/rotating carrier receiving the actual film wrench.
    Finite {
        /// Original positive mass and complete inertia about its center of mass.
        material: &'a ResolvedBody,
        /// Carrier state at the assembly's CURRENT accepted physical time.
        state: BodyStateSi,
        /// Attached plane point relative to the center of mass, in body axes [m].
        point_body_m: [f64; 3],
        /// Body-fixed normal pointing into the spheres' permitted half-space.
        normal_body: [f64; 3],
        /// Wall-normal missing-resistance cutoff [m].
        cutoff_gap_m: f64,
        /// Positive wall-normal gap regularization, no larger than cutoff [m].
        roughness_gap_m: f64,
        /// Optional tangent/parallel-spin [cutoff, roughness] gaps [m].
        tangential_gaps_m: Option<[f64; 2]>,
    },
}

/// Original material spheres, resolved fluid and mixed moving film surfaces.
///
/// All constituents share one accepted clock and the original simultaneous
/// resistance solve. Finite carriers recoil and their attached planes follow
/// their native rigid motion. No surface is independently mutable. Viscosity
/// always comes from the assembly's current fluid snapshot.
pub struct ResolvedImmersedSurfaces<'a, F: FluidTransport + ?Sized, const N: usize> {
    assembly: ResolvedImmersedBodies<'a, F, N>,
    sources: Vec<MaterialFilmSurfaceSi<'a>>,
    surfaces: Vec<MovingLubricationSurface3>,
    surface_impulses: Vec<[f64; 6]>,
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedBodies<'a, F, N> {
    /// Attach up to 64 prescribed or finite material surfaces at the current time.
    ///
    /// Original fluid, sphere poses/momenta, reduced gravity and accepted history
    /// are retained. Finite carriers use the same SI/lattice mass, full inertia
    /// and momentum conversion as the immersed solids. They have no sampled
    /// fluid interior, so no buoyancy, displaced mass, or gravity is inferred;
    /// supply their complete external force/couple explicitly on each step.
    ///
    /// Native geometry, clearance, speed, film/mode and conservative sweep
    /// admission remains authoritative. This supplies neither contact nor a
    /// resolved moving-wall boundary or thermal destination for film loss.
    pub fn with_moving_surfaces(
        self,
        sources: &[MaterialFilmSurfaceSi<'a>],
        gate: &CancelGate,
    ) -> Result<ResolvedImmersedSurfaces<'a, F, N>, Error> {
        poll(gate)?;
        if !matches!(self.native.surfaces, CarrierSurfaces::Spheres { .. }) {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        }
        if sources.len() > 64 {
            return Err(Error::Invalid("at most 64 moving film surfaces"));
        }
        let native = &self.native;
        let dx = native.scaling.cell_size_m();
        let dt = native.scaling.timestep_s();
        let mut surfaces = Vec::with_capacity(sources.len());
        for &source in sources {
            poll(gate)?;
            surfaces.push(match source {
                MaterialFilmSurfaceSi::Prescribed {
                    plane,
                    axial_spin_rad_s,
                } => MovingLubricationSurface3::Prescribed(
                    plane
                        .lattice(&native.scaling)?
                        .with_axial_spin(finite(axial_spin_rad_s * dt)?)?,
                ),
                MaterialFilmSurfaceSi::Finite {
                    material,
                    state,
                    point_body_m,
                    normal_body,
                    cutoff_gap_m,
                    roughness_gap_m,
                    tangential_gaps_m,
                } => {
                    let body = body_to_lattice(
                        &material.inertia(),
                        state,
                        dx,
                        dt,
                        native.mass_unit,
                        native.inertia_unit,
                        native.angular_impulse_unit,
                    )?;
                    let mut surface = MovingFiniteLubricationPlane3::new(
                        body,
                        convert(point_body_m, 1.0 / dx)?,
                        normal_body,
                        finite(cutoff_gap_m / dx)?,
                        finite(roughness_gap_m / dx)?,
                    )?;
                    if let Some([cutoff, roughness]) = tangential_gaps_m {
                        surface = surface
                            .with_tangential(finite(cutoff / dx)?, finite(roughness / dx)?)?;
                    }
                    MovingLubricationSurface3::Finite(surface)
                }
            });
        }
        poll(gate)?;
        Ok(ResolvedImmersedSurfaces {
            assembly: self,
            sources: sources.to_vec(),
            surfaces,
            surface_impulses: vec![[0.0; 6]; sources.len()],
        })
    }
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedSurfaces<'a, F, N> {
    /// Original source/scaling and accepted dry-sphere/fluid history, read-only.
    pub const fn assembly(&self) -> &ResolvedImmersedBodies<'a, F, N> {
        &self.assembly
    }

    /// Immutable inputs, including attachment-time states, in stable surface order.
    pub fn source_surfaces_si(&self) -> &[MaterialFilmSurfaceSi<'a>] {
        &self.sources
    }

    /// Authoritative accepted native surface geometry/momenta in lattice units.
    pub fn surfaces_lattice(&self) -> &[MovingLubricationSurface3] {
        &self.surfaces
    }

    /// Accepted finite-carrier state in SI units; prescribed reservoirs return None.
    pub fn surface_body_state_si(&self, index: usize) -> Result<Option<BodyStateSi>, Error> {
        match self.surface(index)? {
            MovingLubricationSurface3::Prescribed(_) => Ok(None),
            MovingLubricationSurface3::Finite(surface) => body_to_si(
                &surface.body(),
                self.assembly.scaling().cell_size_m(),
                self.assembly.scaling().timestep_s(),
                self.assembly.native.angular_impulse_unit,
            )
            .map(Some),
        }
    }

    /// Accepted world point [m] and unit normal for the actual attached surface.
    pub fn surface_geometry_si(&self, index: usize) -> Result<([f64; 3], [f64; 3]), Error> {
        let (point, normal) = match self.surface(index)? {
            MovingLubricationSurface3::Prescribed(surface) => (surface.point(), surface.normal()),
            MovingLubricationSurface3::Finite(surface) => (surface.point(), surface.normal()),
        };
        Ok((
            convert(point, self.assembly.scaling().cell_size_m())?,
            normal,
        ))
    }

    /// Actual rigid surface velocity [m/s] at a physical world point [m].
    /// The original local low-Mach check applies; an infinite rotating plane
    /// has no global speed bound. This observation does not advance state.
    pub fn surface_velocity_at_si(
        &self,
        index: usize,
        point_m: [f64; 3],
    ) -> Result<[f64; 3], Error> {
        let scaling = self.assembly.scaling();
        let point = convert(point_m, 1.0 / scaling.cell_size_m())?;
        let velocity = match self.surface(index)? {
            MovingLubricationSurface3::Prescribed(surface) => surface.velocity_at(point)?,
            MovingLubricationSurface3::Finite(surface) => surface.velocity_at(point)?,
        };
        convert(velocity, scaling.cell_size_m() / scaling.timestep_s())
    }

    fn surface(&self, index: usize) -> Result<&MovingLubricationSurface3, Error> {
        self.surfaces
            .get(index)
            .ok_or(Error::Invalid("moving film surface index out of range"))
    }

    /// Refresh bulk and every pair/surface film viscosity from one fluid source.
    /// Density changes retain the original conservative-transfer refusal.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.assembly.rebind_fluid(fluid, gate)
    }

    /// Prepare the original initial flow without changing surface or body states.
    pub fn initialize_uniform_flow_si(
        &mut self,
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly.initialize_uniform_flow_si(velocity_m_s, gate)
    }

    /// Hold new SI open-reservoir targets without changing accepted fluid,
    /// sphere or moving-surface state. The original assembly owns admission
    /// and the following common step's reservoir momentum account.
    pub fn rebind_open_boundaries_si(
        &mut self,
        faces: [FaceBoundarySi; 6],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly.rebind_open_boundaries_si(faces, gate)
    }

    /// Hold physical cell-center inlet velocities [m/s] with the original
    /// source-scaled speed admission and native reservoir momentum accounting.
    /// Returned profile diagnostics retain their explicit lattice units.
    pub fn set_velocity_profile_si(
        &mut self,
        face: Face3,
        flux: InletFluxSi,
        gate: &CancelGate,
        sample: impl FnMut([f64; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, Error> {
        self.assembly.set_velocity_profile_si(face, flux, gate, sample)
    }

    /// Restore the held uniform face target without advancing accepted state.
    pub fn clear_velocity_profile(&mut self, face: Face3, gate: &CancelGate) -> Result<(), Error> {
        self.assembly.clear_velocity_profile(face, gate)
    }

    /// Hold a uniform external fluid volume force [N/m^3], with no immediate impulse.
    pub fn rebind_uniform_fluid_force_si(
        &mut self,
        force_density_n_m3: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_uniform_fluid_force_si(force_density_n_m3, gate)
    }

    /// Hold the original x-fastest spatial fluid force field [N/m^3].
    pub fn rebind_fluid_force_field_si(
        &mut self,
        force_density_n_m3: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_fluid_force_field_si(force_density_n_m3, gate)
    }

    /// Advance actual fluid, dry spheres and mixed film surfaces together.
    ///
    /// Both load arrays contain WORLD center forces [N] followed by center
    /// couples [N m]. Sphere loads exclude their already-owned reduced weight.
    /// Surface loads are complete external loads, in stable descriptor order;
    /// prescribed surfaces require exact zeros. Finite carriers receive their
    /// film wrench once INSIDE the common solve, not through these authored loads.
    ///
    /// The unchanged native report is IN LATTICE UNITS. Its prescribed film
    /// work belongs to actuator ports; finite-carrier film work changes kinetic
    /// storage. The aggregate surface work is therefore not all reservoir work.
    /// Use `assembly()` conversions for momentum/work observations.
    ///
    /// Every conversion and the next clock value is checked before the original
    /// atomic publication. Refusal leaves all populations, poses, momenta and
    /// accepted time intact; only reusable scratch may change.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forces_si(
        &mut self,
        sphere_loads: [[f64; 6]; N],
        surface_loads: &[[f64; 6]],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedSurfacesExchange3<N>, Error> {
        poll(gate)?;
        if surface_loads.len() != self.surfaces.len() {
            return Err(Error::Invalid("surface/load count mismatch"));
        }
        let loads = self
            .assembly
            .prepare_step(sphere_loads, moment_origin_m, gate)?;
        let native = &mut self.assembly.native;
        for (index, authored) in surface_loads.iter().enumerate() {
            poll(gate)?;
            if matches!(
                self.surfaces[index],
                MovingLubricationSurface3::Prescribed(_)
            ) && authored.iter().any(|&x| x != 0.0)
            {
                return Err(Error::Invalid(
                    "prescribed surface cannot receive an authored load",
                ));
            }
            for (axis, &load) in authored.iter().enumerate() {
                let unit = if axis < 3 {
                    native.impulse_unit
                } else {
                    native.angular_impulse_unit
                };
                self.surface_impulses[index][axis] =
                    finite(load * (native.scaling.timestep_s() / unit))?;
            }
        }
        let CarrierSurfaces::Spheres { radii, film } = native.surfaces else {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        };
        let report = native.grid.step_with_dry_lubricated_spheres_and_surfaces(
            &mut native.samplers,
            radii,
            &mut native.bodies,
            loads.impulses,
            film,
            &mut self.surfaces,
            &self.surface_impulses,
            loads.origin,
            controls,
            gate,
        )?;
        native.accepted_steps = loads.next;
        Ok(report)
    }
}
