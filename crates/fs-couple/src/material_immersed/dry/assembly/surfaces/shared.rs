//! Several body-fixed film patches on each original material carrier.

use super::*;
use fs_lbm::d3q19::{
    DryLubricatedCarriersExchange3, LubricationPlanePatch3, MovingCarrierSurface3,
};

/// One independent physical carrier, with its original complete dry inertia.
///
/// All patches naming its slot share this state, one external wrench and one
/// kinetic-energy account. Equal material references in DIFFERENT slots still
/// describe independent physical bodies; source identity is not runtime identity.
#[derive(Clone, Copy, Debug)]
pub struct MaterialFilmCarrierSi<'a> {
    /// Original positive mass and full center-of-mass inertia.
    pub material: &'a ResolvedBody,
    /// Carrier state at the assembly's current accepted physical time.
    pub state: BodyStateSi,
}

/// A prescribed film reservoir or a patch on an explicitly identified carrier.
///
/// Geometry is the existing local planar missing-resistance model, not a finite
/// face mesh or resolved no-slip wall. Film scales remain explicit; every patch
/// gets viscosity from the assembly's single current fluid snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialFilmPatchSi {
    /// Prescribed translation and axial spin with its own actuator-work port.
    Prescribed {
        /// World geometry, translation and explicit film scales in SI units.
        plane: PlaneFilmSi,
        /// Signed spin about the plane normal through its translating point [rad/s].
        axial_spin_rad_s: f64,
    },
    /// Body-fixed geometry referring to one `MaterialFilmCarrierSi` slot.
    Attached {
        /// Index into the separately declared, stable carrier list.
        carrier: usize,
        /// Offset from the carrier's center of mass, in body axes [m].
        point_body_m: [f64; 3],
        /// Body-fixed normal into the spheres' permitted half-space.
        normal_body: [f64; 3],
        /// Wall-normal missing-resistance cutoff [m].
        cutoff_gap_m: f64,
        /// Positive gap regularization, no larger than the cutoff [m].
        roughness_gap_m: f64,
        /// Optional tangent/parallel-spin [cutoff, roughness] gaps [m].
        tangential_gaps_m: Option<[f64; 2]>,
    },
}

/// Dry spheres, actual fluid and shared finite carriers on one accepted clock.
///
/// Multiple attached patches participate in one simultaneous native resistance
/// solve. Each carrier receives their summed wrench and advances only once.
/// Carrier momentum/kinetic storage and per-patch film exchange stay distinct.
pub struct ResolvedImmersedSharedSurfaces<'a, F: FluidTransport + ?Sized, const N: usize> {
    assembly: ResolvedImmersedBodies<'a, F, N>,
    carrier_sources: Vec<MaterialFilmCarrierSi<'a>>,
    patch_sources: Vec<MaterialFilmPatchSi>,
    carriers: Vec<MovingPartialSaturationBody3>,
    surfaces: Vec<MovingCarrierSurface3>,
    carrier_impulses: Vec<[f64; 6]>,
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedBodies<'a, F, N> {
    /// Attach up to 64 material carriers and 64 prescribed/body-fixed patches.
    ///
    /// Every carrier must be referenced by at least one attached patch. Its
    /// material mass, complete inertia and SI momentum are converted once, and
    /// all accepted sphere/fluid history survives attachment. Carrier states
    /// belong to the current accepted time, even when attachment follows steps.
    ///
    /// Carriers have no sampled fluid interior or inferred displaced volume;
    /// supply their complete external loads, including any gravity, explicitly.
    /// Native geometry, film, speed, conservative sweep and coupling admission
    /// remain authoritative. This adds neither resolved moving walls nor contact,
    /// a thermal receiver, or arbitrary-face lubrication.
    pub fn with_shared_surfaces(
        self,
        carrier_sources: &[MaterialFilmCarrierSi<'a>],
        patch_sources: &[MaterialFilmPatchSi],
        gate: &CancelGate,
    ) -> Result<ResolvedImmersedSharedSurfaces<'a, F, N>, Error> {
        poll(gate)?;
        if !matches!(self.native.surfaces, CarrierSurfaces::Spheres { .. }) {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        }
        if carrier_sources.len() > 64 || patch_sources.len() > 64 {
            return Err(Error::Invalid("at most 64 film carriers and 64 patches"));
        }
        let mut referenced = [false; 64];
        for source in patch_sources {
            poll(gate)?;
            if let MaterialFilmPatchSi::Attached { carrier, .. } = *source {
                if carrier >= carrier_sources.len() {
                    return Err(Error::Invalid("film patch carrier index out of range"));
                }
                referenced[carrier] = true;
            }
        }
        if referenced[..carrier_sources.len()].contains(&false) {
            return Err(Error::Invalid("every film carrier needs an attached patch"));
        }
        let native = &self.native;
        let dx = native.scaling.cell_size_m();
        let dt = native.scaling.timestep_s();
        let mut carriers = Vec::with_capacity(carrier_sources.len());
        for source in carrier_sources {
            poll(gate)?;
            carriers.push(body_to_lattice(
                &source.material.inertia(),
                source.state,
                dx,
                dt,
                native.mass_unit,
                native.inertia_unit,
                native.angular_impulse_unit,
            )?);
        }
        let mut surfaces = Vec::with_capacity(patch_sources.len());
        for &source in patch_sources {
            poll(gate)?;
            let surface = match source {
                MaterialFilmPatchSi::Prescribed {
                    plane,
                    axial_spin_rad_s,
                } => MovingCarrierSurface3::Prescribed(
                    plane
                        .lattice(&native.scaling)?
                        .with_axial_spin(finite(axial_spin_rad_s * dt)?)?,
                ),
                MaterialFilmPatchSi::Attached {
                    carrier,
                    point_body_m,
                    normal_body,
                    cutoff_gap_m,
                    roughness_gap_m,
                    tangential_gaps_m,
                } => {
                    let mut patch_body = LubricationPlanePatch3::new(
                        convert(point_body_m, 1.0 / dx)?,
                        normal_body,
                        finite(cutoff_gap_m / dx)?,
                        finite(roughness_gap_m / dx)?,
                    )?;
                    if let Some([cutoff, roughness]) = tangential_gaps_m {
                        patch_body = patch_body
                            .with_tangential(finite(cutoff / dx)?, finite(roughness / dx)?)?;
                    }
                    MovingCarrierSurface3::Attached {
                        carrier,
                        patch_body,
                    }
                }
            };
            let geometry = surface.geometry(&carriers)?;
            surface.velocity_at(&carriers, geometry.point())?;
            surfaces.push(surface);
        }
        poll(gate)?;
        Ok(ResolvedImmersedSharedSurfaces {
            assembly: self,
            carrier_sources: carrier_sources.to_vec(),
            patch_sources: patch_sources.to_vec(),
            carriers,
            surfaces,
            carrier_impulses: vec![[0.0; 6]; carrier_sources.len()],
        })
    }
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedSharedSurfaces<'a, F, N> {
    /// Original sources/scaling and accepted dry-sphere/fluid state, read-only.
    pub const fn assembly(&self) -> &ResolvedImmersedBodies<'a, F, N> {
        &self.assembly
    }
    /// Original independent carriers and their attachment-time SI states.
    pub fn source_carriers_si(&self) -> &[MaterialFilmCarrierSi<'a>] {
        &self.carrier_sources
    }
    /// Original patch geometry and explicit carrier indices, in stable order.
    pub fn source_patches_si(&self) -> &[MaterialFilmPatchSi] {
        &self.patch_sources
    }
    /// Each actual accepted carrier state appears exactly once, in lattice units.
    pub fn carriers_lattice(&self) -> &[MovingPartialSaturationBody3] {
        &self.carriers
    }
    /// Prescribed planes and body-fixed patch bindings, without mutable aliases.
    pub fn surfaces_lattice(&self) -> &[MovingCarrierSurface3] {
        &self.surfaces
    }
    /// Accepted state of one unique finite carrier in physical units.
    pub fn carrier_state_si(&self, index: usize) -> Result<BodyStateSi, Error> {
        let body = self
            .carriers
            .get(index)
            .ok_or(Error::Invalid("film carrier index out of range"))?;
        body_to_si(
            body,
            self.assembly.scaling().cell_size_m(),
            self.assembly.scaling().timestep_s(),
            self.assembly.native.angular_impulse_unit,
        )
    }
    /// Actual world point [m] and unit normal, using the shared accepted pose.
    pub fn surface_geometry_si(&self, index: usize) -> Result<([f64; 3], [f64; 3]), Error> {
        let geometry = self.surface(index)?.geometry(&self.carriers)?;
        Ok((
            convert(geometry.point(), self.assembly.scaling().cell_size_m())?,
            geometry.normal(),
        ))
    }
    /// Actual rigid velocity [m/s] at a world point [m]; local native speed limits apply.
    pub fn surface_velocity_at_si(
        &self,
        index: usize,
        point_m: [f64; 3],
    ) -> Result<[f64; 3], Error> {
        let scaling = self.assembly.scaling();
        let velocity = self.surface(index)?.velocity_at(
            &self.carriers,
            convert(point_m, 1.0 / scaling.cell_size_m())?,
        )?;
        convert(velocity, scaling.cell_size_m() / scaling.timestep_s())
    }
    fn surface(&self, index: usize) -> Result<&MovingCarrierSurface3, Error> {
        self.surfaces
            .get(index)
            .ok_or(Error::Invalid("film patch index out of range"))
    }
    /// Refresh bulk and all film viscosity together, preserving accepted history.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.assembly.rebind_fluid(fluid, gate)
    }
    /// Prepare the original initial flow without changing carriers or patches.
    pub fn initialize_uniform_flow_si(
        &mut self,
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly.initialize_uniform_flow_si(velocity_m_s, gate)
    }
    /// Hold new SI open-reservoir targets without changing accepted fluid,
    /// spheres, shared carriers or patches. The original assembly owns admission
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

    /// Hold physical cell-center pressure perturbations [Pa] through the original
    /// material scaling and native reservoir momentum owner. Diagnostics use LU.
    pub fn set_pressure_profile_si(
        &mut self,
        face: Face3,
        gate: &CancelGate,
        sample: impl FnMut([f64; 3]) -> Result<f64, BoundaryProfileError3>,
    ) -> Result<PressureProfileReport3, Error> {
        self.assembly.set_pressure_profile_si(face, gate, sample)
    }

    /// Restore the held uniform pressure target without advancing accepted state.
    pub fn clear_pressure_profile(&mut self, face: Face3, gate: &CancelGate) -> Result<(), Error> {
        self.assembly.clear_pressure_profile(face, gate)
    }
    /// Hold an external uniform fluid volume force [N/m^3], without an immediate impulse.
    pub fn rebind_uniform_fluid_force_si(
        &mut self,
        force_density_n_m3: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_uniform_fluid_force_si(force_density_n_m3, gate)
    }
    /// Hold an external x-fastest spatial fluid force field [N/m^3].
    pub fn rebind_fluid_force_field_si(
        &mut self,
        force_density_n_m3: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_fluid_force_field_si(force_density_n_m3, gate)
    }

    /// Advance all bodies, patches and populations through one original solve.
    ///
    /// Loads are WORLD center forces [N] then center couples [N m], with exactly
    /// one load per unique carrier, irrespective of its patch count. Sphere
    /// loads exclude already-owned reduced gravity; carrier loads are complete.
    /// Film reactions are applied inside the shared solve, never authored again.
    ///
    /// The unchanged report is in LATTICE UNITS. Per-carrier trajectory storage
    /// occurs once; per-patch work is its contribution to that same carrier or
    /// a prescribed actuator port. Use `assembly()` for observation conversions.
    /// Every fallible conversion and the next clock check precedes publication.
    /// Refusal preserves all accepted states; reusable trial scratch may change.
    pub fn step_forces_si(
        &mut self,
        sphere_loads: [[f64; 6]; N],
        carrier_loads: &[[f64; 6]],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedCarriersExchange3<N>, Error> {
        poll(gate)?;
        if carrier_loads.len() != self.carriers.len() {
            return Err(Error::Invalid("carrier/load count mismatch"));
        }
        let loads = self
            .assembly
            .prepare_step(sphere_loads, moment_origin_m, gate)?;
        let native = &mut self.assembly.native;
        for (index, authored) in carrier_loads.iter().enumerate() {
            poll(gate)?;
            for (axis, &load) in authored.iter().enumerate() {
                let unit = if axis < 3 {
                    native.impulse_unit
                } else {
                    native.angular_impulse_unit
                };
                self.carrier_impulses[index][axis] =
                    finite(load * (native.scaling.timestep_s() / unit))?;
            }
        }
        let CarrierSurfaces::Spheres { radii, film } = native.surfaces else {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        };
        let report = native.grid.step_with_dry_lubricated_spheres_and_carriers(
            &mut native.samplers,
            radii,
            &mut native.bodies,
            loads.impulses,
            film,
            &mut self.carriers,
            &self.carrier_impulses,
            &mut self.surfaces,
            loads.origin,
            controls,
            gate,
        )?;
        native.accepted_steps = loads.next;
        Ok(report)
    }
}
