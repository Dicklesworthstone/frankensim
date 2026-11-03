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
/// rotation, contact and a thermal destination for film loss are not supplied.
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
    fn lattice<F: FluidTransport + ?Sized>(
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
        poll(gate)?;
        if !matches!(self.native.surfaces, CarrierSurfaces::Spheres { .. }) {
            return Err(Error::Invalid("plane films require spherical dry bodies"));
        }
        if planes.len() > 64 {
            return Err(Error::Invalid("at most 64 prescribed film planes"));
        }
        let mut native = Vec::with_capacity(planes.len());
        for &plane in planes {
            poll(gate)?;
            native.push(plane.lattice(self.scaling())?);
        }
        poll(gate)?;
        Ok(ResolvedImmersedPlanes {
            assembly: self,
            source_planes: planes.to_vec(),
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
        let report = native.grid.step_with_dry_lubricated_spheres_and_planes(
            &mut native.samplers,
            radii,
            &mut native.bodies,
            loads.impulses,
            film,
            &self.planes,
            loads.origin,
            controls,
            gate,
        )?;
        core::mem::swap(&mut self.planes, &mut self.next_planes);
        native.accepted_steps = loads.next;
        Ok(report)
    }
}
