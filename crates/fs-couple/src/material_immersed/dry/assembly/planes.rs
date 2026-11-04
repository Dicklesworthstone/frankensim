//! Prescribed film surfaces sharing the material assembly's accepted clock.

use super::*;
use fs_lbm::d3q19::{DryLubricatedPlanesExchange3, LubricationPlane3};

/// Original geometry and prescribed translation of one missing-film surface.
///
/// The normal points into the spheres' permitted fluid half-space. A plane
/// supplies only the native missing-film law and its reaction/work port. It
/// does not create or translate a resolved no-slip fluid boundary. The scenario
/// must supply the matching resolved boundary separately; a stationary plane
/// can accompany an existing stationary domain wall. Finite-mass plane recoil,
/// tilting rotation, contact and a thermal destination for film loss are not supplied.
/// Axial spin is supplied separately by `with_spinning_prescribed_planes`, keeping
/// the original translating-plane descriptor and callers source-compatible.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneFilmSi {
    /// Point on the plane at attachment time [m].
    pub point_m: [f64; 3],
    /// Nonzero world normal, normalized by the original native owner.
    pub normal: [f64; 3],
    /// Constant prescribed world translation [m/s].
    pub velocity_m_s: [f64; 3],
    /// Wall-normal missing-resistance cutoff [m], independently calibrated.
    pub cutoff_gap_m: f64,
    /// Positive wall-normal gap regularization, no larger than cutoff [m].
    pub roughness_gap_m: f64,
    /// Optional tangent/parallel-spin [cutoff, roughness] gaps [m]. Native
    /// asymptotic admission requires cutoff <= 0.02 of every sphere radius.
    pub tangential_gaps_m: Option<[f64; 2]>,
}

impl PlaneFilmSi {
    pub(super) fn lattice<F: FluidTransport + ?Sized>(
        self,
        scaling: &MaterialScaling<'_, F>,
    ) -> Result<LubricationPlane3, Error> {
        let dx = scaling.cell_size_m();
        let mut plane = LubricationPlane3::new(
            convert(self.point_m, 1.0 / dx)?,
            self.normal,
            convert(self.velocity_m_s, scaling.timestep_s() / dx)?,
            finite(self.cutoff_gap_m / dx)?,
            finite(self.roughness_gap_m / dx)?,
        )?;
        if let Some([cutoff, roughness]) = self.tangential_gaps_m {
            plane = plane.with_tangential(finite(cutoff / dx)?, finite(roughness / dx)?)?;
        }
        Ok(plane)
    }
}

/// Original material dry spheres and prescribed films on one accepted clock.
///
/// Created by [`ResolvedImmersedBodies::with_prescribed_planes`]. Source-owned
/// viscosity drives sphere-pair and plane resistance inside every original
/// implicit dry-body candidate. Plane positions advance exactly once with an
/// accepted fluid/body interval. Rejected nonlinear probes never advance them.
/// Individual plane reactions remain separate from authored body loads and
/// resolved-wall exchange. There is no independent mutable assembly access.
pub struct ResolvedImmersedPlanes<'a, F: FluidTransport + ?Sized, const N: usize> {
    assembly: ResolvedImmersedBodies<'a, F, N>,
    source_planes: Vec<PlaneFilmSi>,
    source_plane_axial_spins_rad_s: Vec<f64>,
    planes: Vec<LubricationPlane3>,
    next_planes: Vec<LubricationPlane3>,
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedBodies<'a, F, N> {
    /// Attach up to 64 prescribed surfaces to this spherical material assembly.
    /// Supplied points describe the CURRENT accepted time, including when this
    /// is called after earlier sphere-only steps. The original dry source,
    /// populations, pose/momentum, gravity and accepted-step count are retained.
    ///
    /// Viscosity comes from the assembly's fluid; no independent scalar is
    /// accepted. Existing total-film/mode budgets include pairs AND planes.
    /// Native duplicate-plane, positive-gap, sweep and tangent-asymptotic
    /// admission remains authoritative on each step. Empty planes preserve the
    /// sphere-only result exactly. Ellipsoidal film attachment refuses.
    pub fn with_prescribed_planes(
        self,
        planes: &[PlaneFilmSi],
        gate: &CancelGate,
    ) -> Result<ResolvedImmersedPlanes<'a, F, N>, Error> {
        self.attach_prescribed_planes(planes, None, gate)
    }

    /// Attach prescribed planar surfaces with signed angular speeds [rad/s]
    /// about their own normals through the supplied points. Each spin axis
    /// translates with its original `PlaneFilmSi::velocity_m_s`.
    ///
    /// Rates correspond one-to-one with `planes` and retain their original SI
    /// values. The native rate is `omega * timestep`, with no length factor.
    /// All normal/tangent resistance and work remain in the original simultaneous
    /// fluid/body solve. The matching resolved wall must carry the same motion.
    /// This adds no finite axial torsional law, tilted wall rotation or recoil.
    /// Zero rates recover the existing prescribed-plane path exactly.
    pub fn with_spinning_prescribed_planes(
        self,
        planes: &[PlaneFilmSi],
        axial_spins_rad_s: &[f64],
        gate: &CancelGate,
    ) -> Result<ResolvedImmersedPlanes<'a, F, N>, Error> {
        self.attach_prescribed_planes(planes, Some(axial_spins_rad_s), gate)
    }

    fn attach_prescribed_planes(
        self,
        planes: &[PlaneFilmSi],
        axial_spins_rad_s: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<ResolvedImmersedPlanes<'a, F, N>, Error> {
        poll(gate)?;
        if axial_spins_rad_s.is_some_and(|rates| rates.len() != planes.len()) {
            return Err(Error::Invalid("plane/spin source count mismatch"));
        }
        if !matches!(self.native.surfaces, CarrierSurfaces::Spheres { .. }) {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        }
        if planes.len() > 64 {
            return Err(Error::Invalid("at most 64 prescribed film planes"));
        }
        let mut native = Vec::with_capacity(planes.len());
        for (index, &plane) in planes.iter().enumerate() {
            poll(gate)?;
            let spin = axial_spins_rad_s.map_or(0.0, |rates| rates[index]);
            let lattice_spin = finite(spin * self.scaling().timestep_s())?;
            native.push(
                plane
                    .lattice(self.scaling())?
                    .with_axial_spin(lattice_spin)?,
            );
        }
        poll(gate)?;
        Ok(ResolvedImmersedPlanes {
            assembly: self,
            source_planes: planes.to_vec(),
            source_plane_axial_spins_rad_s: axial_spins_rad_s
                .map_or_else(|| vec![0.0; planes.len()], <[f64]>::to_vec),
            next_planes: native.clone(),
            planes: native,
        })
    }
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedPlanes<'a, F, N> {
    /// Original source/scaling and accepted physical dry states, read-only.
    pub const fn assembly(&self) -> &ResolvedImmersedBodies<'a, F, N> {
        &self.assembly
    }

    /// Original SI inputs in stable plane order, including their attachment pose.
    pub fn source_planes_si(&self) -> &[PlaneFilmSi] {
        &self.source_planes
    }

    /// Original signed SI angular-speed inputs in the stable source-plane order.
    /// The original translating-plane constructor supplies exact zero rates.
    pub fn source_plane_axial_spins_rad_s(&self) -> &[f64] {
        &self.source_plane_axial_spins_rad_s
    }

    /// Authoritative accepted plane geometry and prescribed rate in lattice units.
    pub fn planes_lattice(&self) -> &[LubricationPlane3] {
        &self.planes
    }

    /// Observe an accepted plane in SI units. Observation conversion is separate
    /// from stepping and cannot turn a successful native publication into an error.
    pub fn plane_si(&self, index: usize) -> Result<PlaneFilmSi, Error> {
        let plane = self
            .planes
            .get(index)
            .ok_or(Error::Invalid("prescribed plane index out of range"))?;
        let scaling = self.assembly.scaling();
        Ok(PlaneFilmSi {
            point_m: convert(plane.point(), scaling.cell_size_m())?,
            normal: plane.normal(),
            velocity_m_s: convert(
                plane.velocity(),
                scaling.cell_size_m() / scaling.timestep_s(),
            )?,
            ..self.source_planes[index]
        })
    }

    /// Signed accepted plane rotation about its own normal [rad/s].
    /// Observation conversion remains separate from native state publication.
    pub fn plane_axial_spin_rad_s(&self, index: usize) -> Result<f64, Error> {
        let plane = self
            .planes
            .get(index)
            .ok_or(Error::Invalid("prescribed plane index out of range"))?;
        finite(plane.axial_spin() / self.assembly.scaling().timestep_s())
    }

    /// Actual local rigid wall velocity [m/s] at a physical world point [m].
    /// The native low-Mach admission applies at this point; an infinite rotating
    /// plane has no global speed bound. This does not construct a resolved wall.
    pub fn plane_velocity_at_si(&self, index: usize, point_m: [f64; 3]) -> Result<[f64; 3], Error> {
        let plane = self
            .planes
            .get(index)
            .ok_or(Error::Invalid("prescribed plane index out of range"))?;
        let scaling = self.assembly.scaling();
        let point = convert(point_m, 1.0 / scaling.cell_size_m())?;
        convert(
            plane.velocity_at(point)?,
            scaling.cell_size_m() / scaling.timestep_s(),
        )
    }

    /// Update bulk and ALL film viscosities together while retaining the accepted
    /// plane/fluid/body history. Original density-changing rebind refusal applies.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.assembly.rebind_fluid(fluid, gate)
    }

    /// Prepare physical fluid flow [m/s] at the initial accepted time, retaining
    /// original dry-body states and prescribed plane geometry/rates.
    pub fn initialize_uniform_flow_si(
        &mut self,
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly.initialize_uniform_flow_si(velocity_m_s, gate)
    }

    /// Hold new SI open-reservoir targets without changing accepted fluid,
    /// dry-body or prescribed-plane state. The original assembly owns admission
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

    /// Hold an external fluid volume force [N/m^3] with no immediate impulse or
    /// work. Original fluid/body/plane history and its common clock are retained.
    pub fn rebind_uniform_fluid_force_si(
        &mut self,
        force_density_n_m3: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_uniform_fluid_force_si(force_density_n_m3, gate)
    }

    /// Hold a dense x-fastest spatial fluid volume force [N/m^3], added to the
    /// uniform force. Preserves fluid/body/plane state and the accepted clock.
    pub fn rebind_fluid_force_field_si(
        &mut self,
        force_density_n_m3: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.assembly
            .rebind_fluid_force_field_si(force_density_n_m3, gate)
    }

    /// Held external fluid volume force [N/m^3].
    pub fn uniform_fluid_force_density_n_m3(&self) -> Result<[f64; 3], Error> {
        self.assembly.uniform_fluid_force_density_n_m3()
    }

    /// Advance physical dry spheres, actual fluid and prescribed planes together.
    /// Authored WORLD center forces/couples use N and N m. The original reduced
    /// weight is added once by the same load owner as the sphere-only entry.
    ///
    /// Returns the unchanged native receipt IN LATTICE UNITS. Its plane impulses,
    /// angular impulses (including reaction couples) and work are RECEIVED by
    /// the planes. Their driving reservoirs supply the opposite to this system.
    /// Use `assembly()` unit conversions for observations. Do not apply these
    /// film reactions to the fluid populations or authored body loads again.
    ///
    /// Both current and next prescribed geometries are prepared before the sole
    /// native fluid/body publication. Any failure preserves all accepted states
    /// and accepted time. Only reusable plane/sampler/solver scratch can change.
    pub fn step_forces_si(
        &mut self,
        authored: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedPlanesExchange3<N>, Error> {
        self.step_forces_si_with_fluid_endpoint(authored, moment_origin_m, controls, gate, None)
    }

    /// Advance the original prescribed-plane transaction with the selected SI
    /// fluid-force field scaled at the next endpoint. The uniform force remains
    /// held. Source density, dx/dt, authored loads and reduced weight retain the
    /// same conversion; the native dry solve subtracts each fictitious-interior
    /// source trapezoid once. Film reaction/work retain their original physical owners.
    ///
    /// The native receipt remains in lattice units. Every conversion and
    /// next prescribed-plane position is admitted before the common publication of
    /// fluid, bodies, surface states, endpoint source and accepted clock.
    /// Refusal preserves all accepted history; only reusable scratch may change.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forces_si_with_endpoint_fluid_force_scale(
        &mut self,
        authored: [[f64; 6]; N],
        next_scale: f64,
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedPlanesExchange3<N>, Error> {
        self.step_forces_si_with_fluid_endpoint(
            authored,
            moment_origin_m,
            controls,
            gate,
            Some(next_scale),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_forces_si_with_fluid_endpoint(
        &mut self,
        authored: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
        endpoint_force_scale: Option<f64>,
    ) -> Result<DryLubricatedPlanesExchange3<N>, Error> {
        let loads = self
            .assembly
            .prepare_step(authored, moment_origin_m, gate)?;
        for (next, &plane) in self.next_planes.iter_mut().zip(&self.planes) {
            poll(gate)?;
            *next = plane.translated(1.0)?;
        }
        let native = &mut self.assembly.native;
        let CarrierSurfaces::Spheres { radii, film } = native.surfaces else {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        };
        let report = match endpoint_force_scale {
            None => native.grid.step_with_dry_lubricated_spheres_and_planes(
                &mut native.samplers,
                radii,
                &mut native.bodies,
                loads.impulses,
                film,
                &self.planes,
                loads.origin,
                controls,
                gate,
            )?,
            Some(scale) => native
                .grid
                .step_with_dry_lubricated_spheres_and_planes_and_endpoint_force_scale(
                    &mut native.samplers,
                    radii,
                    &mut native.bodies,
                    loads.impulses,
                    film,
                    &self.planes,
                    scale,
                    loads.origin,
                    controls,
                    gate,
                )?,
        };
        core::mem::swap(&mut self.planes, &mut self.next_planes);
        native.accepted_steps = loads.next;
        Ok(report)
    }
}
