//! Source-owned Stokes flow with physical traction boundaries.

use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologySample, MaterialFluxSystem,
    rheology::RheologyField,
};
use crate::ns::FluxSolution;
use crate::ns::boundary::{
    BoundaryWrench, StokesBoundary, StokesPressureReference, StokesRheologyControls,
    StokesRheologySolution, StokesSolution,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Physical effort per metre of out-of-plane depth, about an explicit origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialBoundaryWrench {
    /// Force per unit depth [N/m].
    pub force_n_per_m: [f64; 2],
    /// Moment per unit depth [(N m)/m = N], positive counterclockwise.
    pub moment_n: f64,
}

impl MaterialBoundaryWrench {
    /// Mechanical power per unit depth [W/m] for an authored rigid virtual
    /// motion: translation [m/s] at the moment origin and angular rate [1/s].
    /// This evaluates conjugate work; it does not advance body/mesh motion or
    /// include transported kinetic energy.
    ///
    /// # Errors
    /// Refuses nonfinite motion or unrepresentable power.
    pub fn rigid_power_w_per_m(
        self,
        translation: [f64; 2],
        angular_rate: f64,
    ) -> Result<f64, MaterialFluxError> {
        let power = self.force_n_per_m[0] * translation[0]
            + self.force_n_per_m[1] * translation[1]
            + self.moment_n * angular_rate;
        if !translation.iter().all(|v| v.is_finite())
            || !angular_rate.is_finite()
            || !power.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "nonfinite physical boundary power",
            ));
        }
        Ok(power)
    }
}

/// Actual discrete pressure/viscous reactions in SI, with momentum transport
/// kept separate. Two-dimensional forces are N/m, moments N, and powers W/m.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBoundaryLoad {
    /// Origin [m] for every moment.
    pub origin_m: [f64; 2],
    /// Total discrete stress effort exerted by fluid on exterior.
    pub stress: MaterialBoundaryWrench,
    /// Supplied open stress, with sign reversed to act on exterior.
    pub imposed: MaterialBoundaryWrench,
    /// Separately declared open-boundary backflow resistance on exterior.
    pub backflow: MaterialBoundaryWrench,
    /// Outward advective momentum and its moment; not a mechanical wall load.
    pub outward_momentum: MaterialBoundaryWrench,
    /// Selected component references: None retains its traction datum;
    /// Some(offset) adds that authored [Pa] offset to its mean-zero pressure.
    pub pressure_offsets_pa: Vec<(usize, Option<f64>)>,
}

pub(super) fn physical_boundary_load(
    flow: &StokesSolution,
    density: f64,
    edges: &[usize],
    origin_m: [f64; 2],
    pressure_offsets_pa: &[Option<f64>],
) -> Result<MaterialBoundaryLoad, MaterialFluxError> {
    let offsets: Vec<_> = pressure_offsets_pa
        .iter()
        .map(|value| value.map(|p| p / density))
        .collect();
    let load = flow
        .reactions
        .resultant(edges, origin_m, &offsets)
        .map_err(MaterialFluxError::Boundary)?;
    let physical = |value: BoundaryWrench| {
        let force_n_per_m = value.force.map(|f| density * f);
        let moment_n = density * value.moment;
        if !force_n_per_m.iter().all(|f| f.is_finite()) || !moment_n.is_finite() {
            Err(MaterialFluxError::Invalid(
                "unrepresentable physical boundary effort",
            ))
        } else {
            Ok(MaterialBoundaryWrench {
                force_n_per_m,
                moment_n,
            })
        }
    };
    Ok(MaterialBoundaryLoad {
        origin_m,
        stress: physical(load.stress)?,
        imposed: physical(load.imposed)?,
        backflow: physical(load.backflow)?,
        outward_momentum: physical(load.outward_momentum)?,
        pressure_offsets_pa: load
            .pressure_offsets
            .iter()
            .map(|(component, _)| (*component, pressure_offsets_pa[*component]))
            .collect(),
    })
}

/// Exactly one physical velocity or total-traction owner per boundary edge.
/// Geometry is in metres. Callbacks must be pure and fixed during the solve.
#[derive(Clone, Copy)]
pub enum MaterialStokesBoundary<'a> {
    /// Prescribed velocity, with the existing strong-normal/SIP-tangent rule.
    Velocity {
        /// Mesh boundary edge index.
        edge: usize,
        /// Physical velocity [m/s].
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
    /// Physical total stress vector (2 mu D(u) - p I)n.
    TotalTraction {
        /// Mesh boundary edge index.
        edge: usize,
        /// Traction [Pa], evaluated at position [m] with outward unit normal.
        value: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    },
}

/// Accepted linear Stokes flow in SI units, with the exact transport snapshot.
pub struct MaterialStokesSolution {
    flow: StokesSolution,
    transport: FluidTransportState,
    pressure_pa: Vec<f64>,
    kinetic_energy_j_per_m: f64,
}

impl MaterialStokesSolution {
    /// Signed outward mass flux [kg/(m s)] and its first spatial moment [kg/s]
    /// about origin_m. Both come from the original retained convective
    /// advector, exactly as the endpoint momentum-transport receipt does.
    /// Stokes-only solutions return zero. No source callback is queried.
    ///
    /// # Errors
    /// Refuses invalid selected edges/origin or physical conversion overflow.
    pub fn outward_mass_flux_moments(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
    ) -> Result<(f64, [f64; 2]), MaterialFluxError> {
        let (flux, moment) = self.flow.reactions
            .volume_flux_moments(edges, origin_m)
            .map_err(MaterialFluxError::Boundary)?;
        let rho = self.transport.density_kg_m3();
        let mass_flux = rho * flux;
        let first_moment = moment.map(|v| rho * v);
        if !mass_flux.is_finite() || !first_moment.iter().all(|v| v.is_finite()) {
            return Err(MaterialFluxError::Invalid("physical mass flux moment overflow"));
        }
        Ok((mass_flux, first_moment))
    }

    /// Evaluate selected boundary reactions using retained endpoint coefficients
    /// and source density. Mean-zero components require explicit offsets [Pa];
    /// traction-datum components require None. Velocity edges are not inferred
    /// to be impermeable walls. No material or boundary callback is queried.
    ///
    /// # Errors
    /// Refuses invalid selections/references or unrepresentable physical loads.
    pub fn boundary_loads(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, MaterialFluxError> {
        physical_boundary_load(
            &self.flow,
            self.transport.density_kg_m3(),
            edges,
            origin_m,
            pressure_offsets_pa,
        )
    }
    /// Original BDM1/P0 field and independently checked equation residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.field
    }
    /// Physical velocity moments [m/s], ready for the next step.
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        &self.flow.field.x[..self.flow.field.n_u]
    }
    /// Pressure [Pa], retaining each natural boundary's declared datum.
    /// Mean-zero components have no known relative pressure offset. This is
    /// mechanical incompressible pressure, not a thermodynamic EOS state.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// Whether pressure uses traction datums, mean-zero gauges or both.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.pressure_reference
    }
    /// Local datum in the original triangle order.
    #[must_use]
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.flow.pressure_reference_per_cell.get(cell).copied()
    }
    /// Immutable source response used for viscosity, inertia and SI conversion.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    /// Physical kinetic energy per metre of out-of-plane depth [J/m].
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }
}

/// Accepted source-backed steady flow and its exact final material responses.
pub struct MaterialStokesRheologySolution<G> {
    flow: StokesRheologySolution,
    pressure_pa: Vec<f64>,
    kinetic_energy_j_per_m: f64,
    samples: Vec<MaterialFluxRheologySample<G>>,
}

impl<G> MaterialStokesRheologySolution<G> {
    /// Physical selected-edge reactions at the accepted constitutive field.
    /// Explicit mean-zero offsets [Pa] and retained traction datums follow
    /// [`MaterialStokesSolution::boundary_loads`]. No source requery occurs.
    ///
    /// # Errors
    /// Refuses invalid selections/references or unrepresentable physical loads.
    pub fn boundary_loads(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, MaterialFluxError> {
        let density = self
            .samples
            .first()
            .ok_or(MaterialFluxError::Invalid("missing endpoint material"))?
            .transport()
            .density_kg_m3();
        physical_boundary_load(
            &self.flow.solution,
            density,
            edges,
            origin_m,
            pressure_offsets_pa,
        )
    }
    /// Primitive velocity/kinematic-pressure coefficients and independently
    /// reassembled original and row-scaled nonlinear endpoint residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.solution.field
    }
    /// Physical velocity coefficients [m/s].
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        let x = &self.flow.solution.field.x;
        &x[..x.len() - self.pressure_pa.len()]
    }
    /// Physical pressure [Pa] in the retained traction datum, or mean-zero
    /// gauge for all-velocity boundaries; not an EOS pressure declaration.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// Declared pressure interpretation, with no hidden shift at open ends.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.solution.pressure_reference
    }
    /// Local datum for a componentwise pressure field, in original cell order.
    #[must_use]
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.flow
            .solution
            .pressure_reference_per_cell
            .get(cell)
            .copied()
    }
    /// Physical kinetic energy per metre of out-of-plane depth.
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }
    /// Exact final source owners, in triangle order, without a reporting requery.
    #[must_use]
    pub fn samples(&self) -> &[MaterialFluxRheologySample<G>] {
        &self.samples
    }
    /// Completed linearized constitutive/flow solves.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.flow.iterations
    }
    /// Final undamped relative velocity correction.
    #[must_use]
    pub const fn relative_update(&self) -> f64 {
        self.flow.relative_update
    }
    /// Larger of the two independently reassembled nonlinear residuals.
    #[must_use]
    pub const fn nonlinear_residual(&self) -> f64 {
        self.flow.nonlinear_residual
    }
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Advance linear incompressible Stokes by one backward-Euler step, with
    /// endpoint force density [N/m³], velocity [m/s] and TOTAL traction [Pa].
    ///
    /// Time is in seconds. All callbacks specify the new time level and can
    /// capture time-varying pressure/stress. Viscosity comes only from the
    /// bound material. Its same frozen density applies to inertia, `mu/rho`,
    /// force, traction, reported pressure and kinetic energy at both time
    /// levels; changing density requires a separate conservative state transfer.
    ///
    /// Finite-step inertia supports entirely traction-driven components without
    /// artificial grounding. Sealed components retain their own zero-mean
    /// pressure and endpoint flux check. This is linear fixed-mesh Stokes;
    /// it adds no convection, backflow law, thermal evolution or EOS pressure.
    ///
    /// # Errors
    /// Preserves checked native admission and convergence failures, refusing
    /// unrepresentable physical outputs without modifying previous_velocity.
    pub fn bdf1_stokes_mixed_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<MaterialStokesSolution, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let scaled_tractions: Vec<_> = boundaries
            .iter()
            .map(|boundary| {
                move |point: [f64; 2], normal: [f64; 2]| match boundary {
                    MaterialStokesBoundary::TotalTraction { value, .. } => {
                        let physical = value(point, normal);
                        [physical[0] / density, physical[1] / density]
                    }
                    MaterialStokesBoundary::Velocity { .. } => {
                        unreachable!("only traction owners request stress samples")
                    }
                }
            })
            .collect();
        let primitive_boundaries: Vec<_> = boundaries
            .iter()
            .zip(&scaled_tractions)
            .map(|(boundary, scaled)| match boundary {
                MaterialStokesBoundary::Velocity { edge, value } => StokesBoundary::Velocity {
                    edge: *edge,
                    value: *value,
                },
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    StokesBoundary::TotalTraction {
                        edge: *edge,
                        value: scaled,
                    }
                }
            })
            .collect();
        let flow = self
            .system
            .bdf1_stokes_mixed_step(
                self.params,
                &|point| {
                    let physical = force_density(point);
                    [physical[0] / density, physical[1] / density]
                },
                &primitive_boundaries,
                previous_velocity,
                timestep_s,
            )
            .map_err(MaterialFluxError::Boundary)?;
        self.material_stokes_solution(flow)
    }

    // One SI recovery for linear and nonlinear mixed-boundary endpoints.
    pub(super) fn material_stokes_solution(
        &self,
        flow: StokesSolution,
    ) -> Result<MaterialStokesSolution, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let pressure_pa: Vec<_> = flow.pressure.iter().map(|p| density * p).collect();
        let kinetic_energy_j_per_m = density * self.system.kinetic_energy(&flow.field.x);
        if pressure_pa.iter().any(|value| !value.is_finite()) || !kinetic_energy_j_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "physical pressure or energy is not representable",
            ));
        }
        Ok(MaterialStokesSolution {
            flow,
            transport: self.transport,
            pressure_pa,
            kinetic_energy_j_per_m,
        })
    }

    /// Solve steady generalized-Newtonian mixed Stokes with physical body
    /// force [N/m³], velocity [m/s] and TOTAL traction [Pa].
    ///
    /// Each constitutive response must retain the reference density. Dynamic
    /// viscosity, force and traction are divided by that SAME density; pressure
    /// and kinetic energy are converted back together. A traction pressure
    /// datum is retained rather than shifted to zero mean.
    ///
    /// resolve receives the actual cell centroid and shear invariant [1/s].
    /// Captured temperature/composition remain fixed during this solve. One
    /// cellwise response supplies volume and both one-sided interface traces.
    /// The shared BDF1 material-query owner caches unchanged cell/rate queries
    /// and retains the exact source owners used by the accepted final residual.
    /// It adds no fit, viscosity floor, thermal evolution or viscoelastic state.
    ///
    /// # Errors
    /// Preserves original material refusals and rejects changed density,
    /// invalid boundaries/conversions, nonfinite states and unconverged solves.
    /// Refusal publishes neither fallback properties nor a partial flow, and
    /// never mutates initial_velocity.
    pub fn solve_stokes_mixed_rheology<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesRheologySolution<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        // Homogeneous closure storage avoids a boxed callback per edge.
        // Admission still happens before any force, traction or source query.
        let scaled_tractions: Vec<_> = boundaries
            .iter()
            .map(|boundary| {
                move |point: [f64; 2], normal: [f64; 2]| match boundary {
                    MaterialStokesBoundary::TotalTraction { value, .. } => {
                        let physical = value(point, normal);
                        [physical[0] / density, physical[1] / density]
                    }
                    MaterialStokesBoundary::Velocity { .. } => {
                        unreachable!("only traction owners request stress samples")
                    }
                }
            })
            .collect();
        let primitive_boundaries: Vec<_> = boundaries
            .iter()
            .zip(&scaled_tractions)
            .map(|(boundary, scaled)| match boundary {
                MaterialStokesBoundary::Velocity { edge, value } => StokesBoundary::Velocity {
                    edge: *edge,
                    value: *value,
                },
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    StokesBoundary::TotalTraction {
                        edge: *edge,
                        value: scaled,
                    }
                }
            })
            .collect();
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = self
            .system
            .solve_stokes_mixed_rheology_checked(
                self.params,
                &|cell, _, rate| field.viscosity(cell, rate),
                &|point| {
                    let physical = force_density(point);
                    [physical[0] / density, physical[1] / density]
                },
                &primitive_boundaries,
                initial_velocity,
                controls,
                &|| field.is_valid(),
            )
            .map_err(MaterialFluxError::Boundary);
        let (flow, samples) = field.finish(result)?;
        let pressure_pa: Vec<_> = flow.solution.pressure.iter().map(|p| density * p).collect();
        let kinetic_energy_j_per_m = density * self.system.kinetic_energy(&flow.solution.field.x);
        if pressure_pa.iter().any(|value| !value.is_finite()) || !kinetic_energy_j_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "physical pressure or energy is not representable",
            ));
        }
        Ok(MaterialStokesRheologySolution {
            flow,
            pressure_pa,
            kinetic_energy_j_per_m,
            samples,
        })
    }
}
