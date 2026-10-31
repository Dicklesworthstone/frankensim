//! Ordinary heavy solids on the original resolved-fluid moving-body solver.
//!
//! The native partial-saturation operator retains interior fluid. A uniform
//! solid therefore contributes EXCESS inertia: M_s-rho_f V and I_s-rho_f I_V.
//! Constant reference-fluid gravity is absorbed into reduced pressure; only
//! (M_s-rho_f V)*g acts on the carrier. No second fluid gravity or extra
//! Archimedes force is applied. Density and viscosity share one source snapshot.
//!
//! This is the heavy-solid, fully immersed, initially resting approximation.
//! Positive excess inertia is essential to this native solve. Equal/lower solid
//! densities refuse, without epsilon masses. Interior rigidity is approximated
//! by the existing partial-saturation law and sampled geometry, not imposed as
//! an exact constraint. Interior fluid must approach rigid motion under the
//! caller's spatial/time refinement; this adapter supplies no calibrated drag,
//! arbitrary-density stability, moving free surface or finite thermal storage.
//!
//! The continuum subtraction is the fictitious-fluid formulation in Glowinski,
//! Pan, Hesla & Joseph, Int. J. Multiphase Flow 25 (1999), 755-794:
//! <https://doi.org/10.1016/S0301-9322(98)00048-2>. Their exact rigid-interior
//! constraint is distinct from this native finite-resolution cell law.

use super::*;
use fs_ga::Mat3;
use fs_lbm::d3q19::BodyFrameEllipsoid3;
use fs_mbd::articulated::SpatialInertia;

/// Original ordinary solid and the two explicitly separated inertias [SI].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplacedInertia {
    /// Original material-defined solid, never relabelled as additional mass.
    pub physical: SpatialInertia,
    /// Reference fluid occupying the SAME analytic volume and full tensor.
    pub displaced: SpatialInertia,
    /// Positive difference supplied to the native retained-fluid carrier.
    pub excess: SpatialInertia,
    pub volume_m3: f64,
}

impl DisplacedInertia {
    fn new(material: &ResolvedBody, density: f64) -> Result<Self, Error> {
        let geometry = material.regions()[0].geometry();
        let physical = material.inertia();
        let mut fluid_tensor = geometry.inertia_m5();
        for entry in &mut fluid_tensor.m {
            *entry = finite(*entry * density)?;
        }
        let displaced = SpatialInertia::new(
            positive(density * geometry.volume_m3())?,
            geometry.centroid_m(),
            fluid_tensor,
        )
        .map_err(MaterialInertiaError::from)?;
        let mass = physical.mass() - displaced.mass();
        if !mass.is_finite() || mass <= 0.0 {
            return Err(Error::Invalid(
                "submerged settling requires solid density strictly above reference fluid density",
            ));
        }
        let excess = SpatialInertia::new(
            mass,
            physical.center_of_mass(),
            Mat3 {
                m: core::array::from_fn(|i| {
                    physical.inertia_com().m[i] - displaced.inertia_com().m[i]
                }),
            },
        )
        .map_err(MaterialInertiaError::from)?;
        Ok(Self {
            physical,
            displaced,
            excess,
            volume_m3: geometry.volume_m3(),
        })
    }
}

/// Ordinary uniform solid sphere with its original source and geometry.
#[derive(Clone, Copy, Debug)]
pub struct SubmergedSphere<'a>(AdditionalMaterialSphere<'a>);
impl<'a> SubmergedSphere<'a> {
    pub fn new(material: &'a ResolvedBody, radius_m: f64) -> Result<Self, Error> {
        Ok(Self(AdditionalMaterialSphere::new(material, radius_m)?))
    }
    pub const fn material(self) -> &'a ResolvedBody {
        self.0.material()
    }
    pub const fn radius_m(self) -> f64 {
        self.0.radius_m()
    }
}

/// Ordinary uniform solid ellipsoid, including its full material-frame tensor.
#[derive(Clone, Copy, Debug)]
pub struct SubmergedEllipsoid<'a>(AdditionalMaterialEllipsoid<'a>);
impl<'a> SubmergedEllipsoid<'a> {
    pub fn new(
        material: &'a ResolvedBody,
        semiaxes_m: [f64; 3],
        body_from_shape: So3,
    ) -> Result<Self, Error> {
        Ok(Self(AdditionalMaterialEllipsoid::new(
            material,
            semiaxes_m,
            body_from_shape,
        )?))
    }
    pub const fn material(self) -> &'a ResolvedBody {
        self.0.material()
    }
    pub const fn semiaxes_m(self) -> [f64; 3] {
        self.0.semiaxes_m()
    }
    pub const fn body_from_shape(self) -> So3 {
        self.0.body_from_shape()
    }
}

pub type SubmergedMaterialSpheres<'a, F, const N: usize> =
    SubmergedMaterialBodies<'a, F, N, SubmergedSphere<'a>>;
pub type SubmergedMaterialEllipsoids<'a, F, const N: usize> =
    SubmergedMaterialBodies<'a, F, N, SubmergedEllipsoid<'a>>;

/// One native fluid/moving-body owner; no separate particle trajectory.
///
/// Reduced-pressure gravity is fixed at construction. Native populations keep
/// their normalized reference-fluid mass; excess carrier and original physical
/// inertia remain separately accessible. Only read-only native access is given.
pub struct SubmergedMaterialBodies<'a, F: FluidTransport + ?Sized, const N: usize, S> {
    native: RetainedFluidBodies<'a, F, N, S>,
    inertias: [DisplacedInertia; N],
    gravity_m_s2: [f64; 3],
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> SubmergedMaterialSpheres<'a, F, N> {
    /// Release ordinary homogeneous spheres from rest in quiescent fluid.
    /// Initial velocities AND angular momenta must be exactly zero; center and
    /// orientation are retained. Shapes, films and fluid obey native admission.
    pub fn new(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [SubmergedSphere<'a>; N],
        initial: [BodyStateSi; N],
        gravity_m_s2: [f64; 3],
        film: FilmGeometry,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let (scaling, inertias) = prepare(
            fluid,
            domain,
            sources.map(SubmergedSphere::material),
            initial,
            gravity_m_s2,
            gate,
        )?;
        let native = RetainedFluidBodies::new_with_scaling(
            scaling,
            domain,
            sources,
            inertias.map(|s| s.excess),
            initial,
            |scaling| {
                let dx = scaling.cell_size_m();
                let mut model = SphereLubrication3::new(
                    scaling.flow().viscosity,
                    film.cutoff_gap_m / dx,
                    film.roughness_gap_m / dx,
                    film.maximum_pairs,
                )?
                .with_resistance_mode_limit(film.maximum_resistance_modes)?;
                if let Some([cutoff, roughness]) = film.tangential_gaps_m {
                    model = model.with_tangential_pairs(cutoff / dx, roughness / dx)?;
                }
                Ok(CarrierSurfaces::Spheres {
                    radii: sources.map(|s| s.radius_m() / dx),
                    film: model,
                })
            },
            gate,
        )?;
        Ok(Self {
            native,
            inertias,
            gravity_m_s2,
        })
    }
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> SubmergedMaterialEllipsoids<'a, F, N> {
    /// Release ordinary homogeneous ellipsoids from rest. Principal shape axes
    /// may be rotated in the body frame; inertia subtraction retains all terms.
    /// Existing enclosing-sphere clearance applies; no ellipsoid film is added.
    pub fn new(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [SubmergedEllipsoid<'a>; N],
        initial: [BodyStateSi; N],
        gravity_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let (scaling, inertias) = prepare(
            fluid,
            domain,
            sources.map(SubmergedEllipsoid::material),
            initial,
            gravity_m_s2,
            gate,
        )?;
        let native = RetainedFluidBodies::new_with_scaling(
            scaling,
            domain,
            sources,
            inertias.map(|s| s.excess),
            initial,
            |scaling| {
                let mut shapes = Vec::with_capacity(N);
                for source in sources {
                    poll(gate)?;
                    shapes.push(source.0.lattice_shape(scaling.cell_size_m())?);
                }
                Ok(CarrierSurfaces::Ellipsoids(
                    shapes
                        .try_into()
                        .map_err(|_| Error::Invalid("ellipsoid count"))?,
                ))
            },
            gate,
        )?;
        Ok(Self {
            native,
            inertias,
            gravity_m_s2,
        })
    }
}

fn prepare<'a, F: FluidTransport + ?Sized, const N: usize>(
    fluid: &'a F,
    domain: ImmersedDomain,
    materials: [&ResolvedBody; N],
    initial: [BodyStateSi; N],
    gravity: [f64; 3],
    gate: &CancelGate,
) -> Result<(MaterialScaling<'a, F>, [DisplacedInertia; N]), Error> {
    poll(gate)?;
    if !(1..=64).contains(&N)
        || !gravity.into_iter().all(f64::is_finite)
        || initial
            .iter()
            .any(|s| s.velocity_m_s != [0.0; 3] || s.angular_momentum_body_kg_m2_s != [0.0; 3])
    {
        return Err(Error::Invalid(
            "settling requires 1..=64 initially resting bodies and finite reduced-pressure gravity",
        ));
    }
    let scaling = plan_scaling_from_material(
        fluid,
        domain.cell_size_m,
        domain.timestep_s,
        domain.maximum_speed_m_s,
    )?;
    let mut inertias = Vec::with_capacity(N);
    for material in materials {
        poll(gate)?;
        inertias.push(DisplacedInertia::new(
            material,
            scaling.transport().density_kg_m3(),
        )?);
    }
    Ok((
        scaling,
        inertias
            .try_into()
            .map_err(|_| Error::Invalid("body count"))?,
    ))
}

/// Separately reconstructed ordinary-body momentum and represented work [SI].
///
/// The native impulse is the ACTUAL population reaction. The reconstructed dry
/// hydrodynamic impulse additionally includes displaced-inertia acceleration
/// and buoyancy; it is conditional on rigid fictitious interior fluid, not a
/// second measured population exchange. Films are reported separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubmergedBodyExchange {
    pub physical_linear_momentum_change_ns: [f64; 3],
    pub physical_spin_momentum_change_nms: [f64; 3],
    pub resolved_carrier_impulse_ns: [f64; 3],
    pub film_impulse_ns: [f64; 3],
    pub displaced_inertia_impulse_ns: [f64; 3],
    pub buoyancy_impulse_ns: [f64; 3],
    pub reconstructed_hydrodynamic_impulse_ns: [f64; 3],
    pub reconstructed_hydrodynamic_spin_impulse_nms: [f64; 3],
    /// Dry momentum minus reconstructed hydrodynamics, weight and authored load.
    pub dry_momentum_residual: [f64; 6],
    /// Delta of -M_excess g dot x, about the same reduced-pressure datum.
    pub reduced_gravity_potential_change_j: f64,
    /// Delta of EXCESS carrier kinetic energy plus reduced gravity potential.
    pub represented_body_energy_change_j: f64,
    /// Native resolved-fluid and film work on this EXCESS carrier.
    pub represented_hydrodynamic_work_j: f64,
    /// Work of caller-authored forces/couples, excluding gravity.
    pub authored_work_j: f64,
    /// Original native free-rotation energy error, never converted to heat.
    pub free_rotation_energy_change_j: f64,
    /// Represented energy minus hydro/authored work and free-rotation error.
    pub represented_work_residual_j: f64,
}

/// An accepted native step, with immutable data for optional SI observations.
///
/// Publication is atomic in the original solver. Converting an already accepted
/// receipt to extreme SI values can refuse without changing that accepted state;
/// these observations are deliberately separate from step acceptance.
#[derive(Clone, Debug)]
pub struct SubmergedStep<const N: usize> {
    pub native: MovingPartialSaturationBodiesExchange3<N>,
    before: [MovingPartialSaturationBody3; N],
    after: [MovingPartialSaturationBody3; N],
    inertias: [DisplacedInertia; N],
    authored: [[f64; 6]; N],
    gravity: [f64; 3],
    dx: f64,
    dt: f64,
    impulse_unit: f64,
    angular_impulse_unit: f64,
    energy_unit: f64,
}

impl<const N: usize> SubmergedStep<N> {
    pub fn body_exchange_si(&self, index: usize) -> Result<SubmergedBodyExchange, Error> {
        let split = self
            .inertias
            .get(index)
            .ok_or(Error::Invalid("body index out of range"))?;
        let a = self.before[index];
        let b = self.after[index];
        let dv = convert(
            core::array::from_fn(|i| b.velocity[i] - a.velocity[i]),
            self.dx / self.dt,
        )?;
        let travel = convert(core::array::from_fn(|i| b.center[i] - a.center[i]), self.dx)?;
        let dry_dp = convert(dv, split.physical.mass())?;
        let displaced_dp = convert(dv, split.displaced.mass())?;
        let buoyancy = convert(self.gravity, -split.displaced.mass() * self.dt)?;
        let weight = convert(self.gravity, split.physical.mass() * self.dt)?;
        let dry_dl = angular_change(a, b, split.physical, self.dt)?;
        let displaced_dl = angular_change(a, b, split.displaced, self.dt)?;
        let hydro = self.native.coupling.bodies[index].hydrodynamic;
        let resolved = convert(hydro.body_impulse, self.impulse_unit)?;
        let film = convert(
            self.native.coupling.lubrication.body_impulses[index],
            self.impulse_unit,
        )?;
        let hydro_spin = convert(
            hydro.body_angular_impulse_about_center,
            self.angular_impulse_unit,
        )?;
        let film_spin = convert(
            self.native.coupling.lubrication.body_spin_impulses[index],
            self.angular_impulse_unit,
        )?;
        let reconstructed = checked(core::array::from_fn(|i| {
            resolved[i] + film[i] + displaced_dp[i] + buoyancy[i]
        }))?;
        let reconstructed_spin = checked(core::array::from_fn(|i| {
            hydro_spin[i] + film_spin[i] + displaced_dl[i]
        }))?;
        let residual = checked(core::array::from_fn(|i| {
            if i < 3 {
                dry_dp[i] - reconstructed[i] - weight[i] - self.authored[index][i] * self.dt
            } else {
                dry_dl[i - 3] - reconstructed_spin[i - 3] - self.authored[index][i] * self.dt
            }
        }))?;
        let potential = finite(-split.excess.mass() * dot(self.gravity, travel))?;
        let motion = self.native.bodies[index];
        let body_energy = finite(motion.kinetic_energy_change * self.energy_unit + potential)?;
        let hydro_work = finite(
            (hydro.body_work + self.native.coupling.lubrication.body_work[index])
                * self.energy_unit,
        )?;
        let mid_omega = convert(
            self.native.coupling.bodies[index].angular_velocity_midpoint,
            1.0 / self.dt,
        )?;
        let authored_work = finite(
            dot(
                [
                    self.authored[index][0],
                    self.authored[index][1],
                    self.authored[index][2],
                ],
                travel,
            ) + dot(
                [
                    self.authored[index][3],
                    self.authored[index][4],
                    self.authored[index][5],
                ],
                mid_omega,
            ) * self.dt,
        )?;
        let free_error = finite(motion.free_rotation_energy_change * self.energy_unit)?;
        let work_residual = finite(body_energy - hydro_work - authored_work - free_error)?;
        Ok(SubmergedBodyExchange {
            physical_linear_momentum_change_ns: dry_dp,
            physical_spin_momentum_change_nms: dry_dl,
            resolved_carrier_impulse_ns: resolved,
            film_impulse_ns: film,
            displaced_inertia_impulse_ns: displaced_dp,
            buoyancy_impulse_ns: buoyancy,
            reconstructed_hydrodynamic_impulse_ns: reconstructed,
            reconstructed_hydrodynamic_spin_impulse_nms: reconstructed_spin,
            dry_momentum_residual: residual,
            reduced_gravity_potential_change_j: potential,
            represented_body_energy_change_j: body_energy,
            represented_hydrodynamic_work_j: hydro_work,
            authored_work_j: authored_work,
            free_rotation_energy_change_j: free_error,
            represented_work_residual_j: work_residual,
        })
    }
}

/// Optional current-state diagnostic, not a constraint or post-commit refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InteriorFluidDiagnostic {
    pub sampled_volume_m3: f64,
    pub relative_volume_error: f64,
    pub sampled_fluid_mass_kg: f64,
    /// Mass-weighted RMS of actual native fluid velocity minus rigid point speed.
    pub rms_slip_m_s: f64,
    pub maximum_slip_m_s: f64,
    /// Weighted interior fluid momentum minus analytic rho_f V * body velocity.
    /// Includes compressibility, occupancy quadrature and velocity error.
    pub linear_momentum_mismatch_ns: [f64; 3],
}

impl<'a, F: FluidTransport + ?Sized, const N: usize, S> SubmergedMaterialBodies<'a, F, N, S> {
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        self.native.scaling()
    }
    pub const fn sources(&self) -> &[S; N] {
        self.native.sources()
    }
    pub const fn inertias(&self) -> &[DisplacedInertia; N] {
        &self.inertias
    }
    pub const fn grid(&self) -> &BoundaryGrid3 {
        self.native.grid()
    }
    /// The authoritative EXCESS carriers; not ordinary dry-body spin momenta.
    pub const fn excess_carriers_lattice(&self) -> &[MovingPartialSaturationBody3; N] {
        self.native.bodies_lattice()
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.native.accepted_steps()
    }
    pub const fn gravity_m_s2(&self) -> [f64; 3] {
        self.gravity_m_s2
    }

    /// Actual center/rotation/velocity with ORIGINAL dry-body angular momentum.
    pub fn state_si(&self, index: usize) -> Result<BodyStateSi, Error> {
        let mut state = self.native.state_si(index)?;
        state.angular_momentum_body_kg_m2_s = array(
            self.inertias[index]
                .physical
                .inertia_com()
                .apply(vector(omega_body(
                    self.native.bodies[index],
                    self.scaling().timestep_s(),
                )?)),
        );
        checked(state.angular_momentum_body_kg_m2_s)?;
        Ok(state)
    }

    /// Refresh source viscosity on the same history. The existing owner refuses
    /// any density change, preserving the original displaced inertia and units.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.native.rebind_fluid(fluid, gate)
    }

    /// Apply reduced weight plus authored WORLD center forces/couples once.
    ///
    /// The grid has no gravity force: reference hydrostatic pressure is removed.
    /// Callers must not include weight/buoyancy again in `authored`. The same
    /// native moving-body/film collision owns all reciprocal impulses and time.
    /// Refusal preserves fluid, carrier pose/momentum and accepted step count.
    /// Geometry accuracy remains visible in the original sampled-volume receipt.
    pub fn step_forces_si(
        &mut self,
        authored: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        gate: &CancelGate,
    ) -> Result<SubmergedStep<N>, Error> {
        poll(gate)?;
        let mut loads = authored;
        for (load, split) in loads.iter_mut().zip(&self.inertias) {
            for axis in 0..3 {
                load[axis] = finite(load[axis] + split.excess.mass() * self.gravity_m_s2[axis])?;
            }
        }
        let before = self.native.bodies;
        let native = self.native.step_forces_si(loads, moment_origin_m, gate)?;
        // No fallible work follows the native atomic publication.
        Ok(SubmergedStep {
            native,
            before,
            after: self.native.bodies,
            inertias: self.inertias,
            authored,
            gravity: self.gravity_m_s2,
            dx: self.native.scaling.cell_size_m(),
            dt: self.native.scaling.timestep_s(),
            impulse_unit: self.native.impulse_unit,
            angular_impulse_unit: self.native.angular_impulse_unit,
            energy_unit: self.native.energy_unit,
        })
    }

    /// Resample the ACCEPTED geometry and measure its fictitious interior fluid.
    /// Bounded O(N*cells*samples^3) work, reusing existing scratch; invoke on a
    /// resolution ladder or diagnostic cadence, not necessarily every tick.
    /// Reports actual error; neither clamps velocities nor repairs volume/mass.
    /// Cancellation changes only sampler scratch and no physical history.
    pub fn interior_diagnostics(
        &mut self,
        gate: &CancelGate,
    ) -> Result<[InteriorFluidDiagnostic; N], Error> {
        match self.native.surfaces {
            CarrierSurfaces::Spheres { radii, .. } => {
                let mut shapes = Vec::with_capacity(N);
                for radius in radii {
                    shapes.push(BodyFrameEllipsoid3::sphere(radius)?);
                }
                self.sample_interiors(
                    core::array::from_fn(|i| &shapes[i] as &dyn BodyFrameShape3),
                    gate,
                )
            }
            CarrierSurfaces::Ellipsoids(shapes) => self.sample_interiors(
                core::array::from_fn(|i| &shapes[i] as &dyn BodyFrameShape3),
                gate,
            ),
        }
    }

    fn sample_interiors(
        &mut self,
        shapes: [&dyn BodyFrameShape3; N],
        gate: &CancelGate,
    ) -> Result<[InteriorFluidDiagnostic; N], Error> {
        let [nx, ny, nz] = self.grid().dimensions();
        let dx = self.scaling().cell_size_m();
        let velocity_unit = dx / self.scaling().timestep_s();
        let mut output = Vec::with_capacity(N);
        for (index, shape) in shapes.into_iter().enumerate() {
            poll(gate)?;
            let body = self.native.bodies[index];
            // Public sampler accepts a sized wrapper around the original shape.
            let shape = BorrowedShape(shape);
            let volume =
                self.native.samplers[index].sample(body.center, body.orientation, &shape, gate)?;
            if volume <= 0.0 {
                return Err(Error::Invalid("solid is unresolved by interior quadrature"));
            }
            let omega = body.angular_velocity_world()?;
            let fractions = self.native.samplers[index].fractions();
            let mut mass = 0.0;
            let mut slip2 = 0.0;
            let mut maximum = 0.0_f64;
            let mut momentum = [0.0; 3];
            for z in 0..nz {
                for y in 0..ny {
                    poll(gate)?;
                    for x in 0..nx {
                        let epsilon = fractions[(z * ny + y) * nx + x];
                        if epsilon == 0.0 {
                            continue;
                        }
                        if self.native.grid.is_solid(x, y, z) {
                            return Err(Error::Invalid(
                                "sampled interior overlaps stationary solid",
                            ));
                        }
                        let amount = epsilon * self.native.grid.density(x, y, z);
                        let u = self.native.grid.velocity(x, y, z);
                        let arm = vector([
                            x as f64 + 0.5 - body.center[0],
                            y as f64 + 0.5 - body.center[1],
                            z as f64 + 0.5 - body.center[2],
                        ]);
                        let rigid = vector(body.velocity) + vector(omega).cross(arm);
                        let slip = vector(u) - rigid;
                        let squared = finite(slip.dot(slip))?;
                        mass += amount;
                        slip2 += amount * squared;
                        maximum = maximum.max(squared);
                        for axis in 0..3 {
                            momentum[axis] += amount * u[axis];
                        }
                    }
                }
            }
            positive(mass)?;
            let volume_m3 = positive(volume * dx * dx * dx)?;
            output.push(InteriorFluidDiagnostic {
                sampled_volume_m3: volume_m3,
                relative_volume_error: finite(
                    (volume_m3 / self.inertias[index].volume_m3 - 1.0).abs(),
                )?,
                sampled_fluid_mass_kg: positive(mass * self.native.mass_unit)?,
                rms_slip_m_s: finite((slip2 / mass).sqrt() * velocity_unit)?,
                maximum_slip_m_s: finite(maximum.sqrt() * velocity_unit)?,
                linear_momentum_mismatch_ns: checked(core::array::from_fn(|i| {
                    momentum[i] * self.native.impulse_unit
                        - self.inertias[index].displaced.mass() * body.velocity[i] * velocity_unit
                }))?,
            });
        }
        poll(gate)?;
        output
            .try_into()
            .map_err(|_| Error::Invalid("interior diagnostic count"))
    }
}

struct BorrowedShape<'a>(&'a dyn BodyFrameShape3);
impl BodyFrameShape3 for BorrowedShape<'_> {
    fn bounding_radius(&self) -> f64 {
        self.0.bounding_radius()
    }
    fn contains(&self, point: Vec3) -> bool {
        self.0.contains(point)
    }
}

fn omega_body(body: MovingPartialSaturationBody3, dt: f64) -> Result<[f64; 3], Error> {
    let omega = body
        .inertia
        .velocity(body.angular_momentum_body)
        .map_err(|source| Error::Native(PartialSaturationError3::Motion { source }))?;
    convert(array(omega), 1.0 / dt)
}
fn angular_change(
    a: MovingPartialSaturationBody3,
    b: MovingPartialSaturationBody3,
    inertia: SpatialInertia,
    dt: f64,
) -> Result<[f64; 3], Error> {
    let world = |state: MovingPartialSaturationBody3| -> Result<Vec3, Error> {
        state
            .orientation
            .rotate(inertia.inertia_com().apply(vector(omega_body(state, dt)?)))
            .map_err(|source| {
                Error::Native(PartialSaturationError3::Motion {
                    source: source.into(),
                })
            })
    };
    checked(array(world(b)? - world(a)?))
}
fn checked<const N: usize>(v: [f64; N]) -> Result<[f64; N], Error> {
    convert(v, 1.0)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
