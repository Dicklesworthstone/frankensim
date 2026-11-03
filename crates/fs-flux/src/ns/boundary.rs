//! Stokes flow with individually owned velocity and total-traction edges.
//!
//! The natural datum is t = (2 nu D(u) - p I)n, with kinematic pressure
//! and traction (both divided by the same constant reference density).
//! Leaving a boundary velocity unknown also retains its continuity equation:
//! a traction boundary fixes that connected component's pressure datum, so no
//! pressure pin is added there. Sealed components retain independent gauges.

use super::{FluxParams, FluxSolution, FluxSystem};
use crate::bdm::edge_gauss_pub;
use fs_sparse::Csr;
use std::cell::Cell;

#[path = "boundary/reaction.rs"]
pub(super) mod reaction;
pub use reaction::{
    BoundaryLoad, BoundaryReactionOwner, BoundaryReactions, BoundaryWrench,
    DiscreteBoundaryReaction,
};
#[path = "boundary/transport.rs"]
mod transport;
pub use transport::BoundaryTransport;

#[path = "boundary/rheology.rs"]
mod rheology;
pub use rheology::{StokesRheologyControls, StokesRheologySolution};
#[path = "boundary/navier_stokes.rs"]
mod navier_stokes;
#[path = "boundary/transient.rs"]
mod transient;
pub use navier_stokes::{
    NavierStokesBoundary, NavierStokesControls, NavierStokesStep, OpenBoundaryPolicy,
};

/// Exactly one condition must own each boundary edge of the supplied mesh.
/// Interior edges cannot be assigned. Callbacks must be pure.
#[derive(Clone, Copy)]
pub enum StokesBoundary<'a> {
    /// Prescribed velocity: normal BDM moments are strong, tangential data SIP.
    Velocity {
        /// Index in TriMesh::edges.
        edge: usize,
        /// Velocity [m/s] at the boundary quadrature point.
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
    /// Natural total traction, including pressure and viscous stress.
    TotalTraction {
        /// Index in TriMesh::edges.
        edge: usize,
        /// Kinematic traction [m²/s²], evaluated with the outward unit normal.
        /// A pressure-only load is -pressure * normal; a fully developed
        /// channel end also has tangential viscous traction.
        value: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    },
}

/// Meaning of the reported kinematic pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StokesPressureReference {
    /// All boundary velocities are prescribed; each component has zero area mean.
    MeanZero,
    /// The total-traction data fixes the pressure datum; output is not shifted.
    /// This is not thermodynamic/EOS absolute pressure.
    TractionDatum,
    /// Some components retain traction datums while others have mean-zero
    /// gauges. Inspect each cell's reference before comparing pressures.
    Componentwise,
}

/// Checked mixed-boundary admission, assembly or solve failure.
#[derive(Debug, Clone, PartialEq)]
pub enum StokesBoundaryError {
    /// Invalid viscosity, solver parameters, geometry size or arithmetic.
    InvalidInput(&'static str),
    /// An entry names an edge outside the mesh.
    EdgeOutOfRange { edge: usize },
    /// An entry attempts to impose a condition on an interior edge.
    InteriorEdge { edge: usize },
    /// More than one entry owns the same edge.
    DuplicateEdge { edge: usize },
    /// A boundary edge has no condition.
    MissingEdge { edge: usize },
    /// Boundary data is non-finite at an assembly quadrature point.
    NonFiniteBoundary { edge: usize },
    /// The body-force callback returned a non-finite value.
    NonFiniteForce,
    /// Viscosity was not positive and finite at an actual assembly sample.
    InvalidViscosity {
        /// Cell owning this one-sided constitutive sample.
        cell: usize,
        /// Physical sample position [m].
        point: [f64; 2],
        /// Refused kinematic viscosity [m²/s].
        value: f64,
    },
    /// Prescribed velocities violate the all-Dirichlet net-flux compatibility.
    IncompatibleFlux {
        /// Connected fluid component (ordered by its first triangle).
        component: usize,
        /// Signed outward integral of prescribed velocity.
        net_flux: f64,
        /// Roundoff allowance for the assembled boundary quadrature.
        tolerance: f64,
    },
    /// Pure traction leaves rigid velocity modes; no hidden grounding is added.
    NoVelocityBoundary,
    /// One disconnected steady component has neither a velocity boundary nor
    /// positive porous drag. Another component cannot remove its rigid modes.
    UnanchoredComponent {
        /// Connected fluid component (ordered by its first triangle).
        component: usize,
    },
    /// The constitutive/momentum iteration refused its candidate.
    Rheology(crate::ns::FluxRheologyError),
    /// The assembled linear system could not be factored.
    SingularSystem,
    /// An explicitly recomputed residual did not meet the requested tolerance.
    NotConverged {
        /// Linear iterations performed.
        iterations: usize,
        /// True Euclidean ||b-Ax|| / max(||b||, f64::MIN_POSITIVE).
        relative_residual: f64,
        /// Relative residual after normalizing each equation by its row norm.
        scaled_relative_residual: f64,
    },
}

impl std::fmt::Display for StokesBoundaryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "mixed-boundary Stokes refused: {self:?}")
    }
}
impl std::error::Error for StokesBoundaryError {}

#[derive(Clone, Copy)]
enum Samples {
    Velocity([[f64; 2]; 2]),
    Traction([[f64; 2]; 2]),
}

/// Samples are retained once so the strong moments and SIP use identical data.
pub(super) struct AdmittedBoundary {
    samples: Vec<Option<Samples>>,
    pub(super) has_traction: Vec<bool>,
    open_policies: Option<Vec<Option<OpenBoundaryPolicy>>>,
    reservoir_velocities: Option<Vec<Option<[[f64; 2]; 2]>>>,
}

impl AdmittedBoundary {
    pub(super) fn open_policy(&self, edge: usize) -> Option<OpenBoundaryPolicy> {
        self.open_policies
            .as_ref()
            .and_then(|policies| policies[edge])
    }

    pub(super) fn reservoir_velocity(&self, edge: usize) -> Option<[[f64; 2]; 2]> {
        self.reservoir_velocities
            .as_ref()
            .and_then(|velocities| velocities[edge])
    }

    pub(super) fn is_traction(&self, edge: usize) -> bool {
        matches!(self.samples[edge], Some(Samples::Traction(_)))
    }

    pub(super) fn velocity(&self, edge: usize, q: usize) -> [f64; 2] {
        match self.samples[edge] {
            Some(Samples::Velocity(values)) => values[q],
            _ => unreachable!("only admitted velocity edges request velocity samples"),
        }
    }

    pub(super) fn traction(&self, edge: usize, q: usize) -> [f64; 2] {
        match self.samples[edge] {
            Some(Samples::Traction(values)) => values[q],
            _ => unreachable!("only admitted traction edges request traction samples"),
        }
    }
}

/// A linear mixed-flow operator bound to its mesh and pressure interpretation.
/// BDF1 includes mass/history; the Oseen entry also freezes convection.
/// The read-only matrix/RHS can also feed the existing discrete transpose solve.
pub struct StokesAssembly<'a, 'm> {
    system: &'a FluxSystem<'m>,
    matrix: Csr,
    rhs: Vec<f64>,
    params: FluxParams,
    reference: StokesPressureReference,
    has_traction: Vec<bool>,
    reactions: reaction::ReactionOperator,
}

/// The admitted discrete field together with pressure in its declared datum.
pub struct StokesSolution {
    /// Velocity DOFs followed by raw P0 pressure coefficients.
    pub field: FluxSolution,
    /// Per-cell kinematic pressure [m²/s²], in pressure_reference.
    pub pressure: Vec<f64>,
    /// Whether an output gauge shift was applied.
    pub pressure_reference: StokesPressureReference,
    /// Actual reference for each cell, always MeanZero or TractionDatum.
    /// Independent mean-zero components have no known pressure difference.
    pub pressure_reference_per_cell: Vec<StokesPressureReference>,
    /// Original discrete boundary reaction functionals at this accepted field.
    /// Stress exchange is separate from outward advective momentum; mechanics
    /// evaluation requires explicit offsets for selected mean-zero components.
    pub reactions: BoundaryReactions,
}

impl StokesAssembly<'_, '_> {
    /// Assembled operator, including strong normal rows and any pressure pin.
    #[must_use]
    pub fn matrix(&self) -> &Csr {
        &self.matrix
    }

    /// Body, velocity and natural total-traction load.
    #[must_use]
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }

    /// Pressure datum selected by the admitted boundary partition.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.reference
    }

    /// Local pressure interpretation; `None` when cell is outside this mesh.
    #[must_use]
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.system.pressure_component.get(cell).map(|&component| {
            if self.has_traction[component] {
                StokesPressureReference::TractionDatum
            } else {
                StokesPressureReference::MeanZero
            }
        })
    }

    /// Require both original and row-scaled residuals without shifting
    /// traction pressure or hiding constraint error in large viscous loads.
    ///
    /// # Errors
    /// Refuses singular systems, non-finite results or exhausted convergence.
    pub fn solve(&self) -> Result<StokesSolution, StokesBoundaryError> {
        let mut field = self
            .system
            .solve_linear_result(&self.matrix, &self.rhs, self.params)?;
        let mut residual = vec![0.0; self.rhs.len()];
        let (relative_residual, scaled_relative_residual) =
            super::sparse::residuals(&self.matrix, &self.rhs, &field.x, &mut residual);
        field.rel_residual = relative_residual;
        field.scaled_rel_residual = scaled_relative_residual;
        if !field.x.iter().all(|value| value.is_finite()) || !field.converged(self.params.tol) {
            return Err(StokesBoundaryError::NotConverged {
                iterations: field.iters,
                relative_residual,
                scaled_relative_residual,
            });
        }
        self.report_solution(field)
    }

    // Share the physical pressure interpretation after either checked linear
    // convergence or checked nonlinear endpoint convergence.
    fn report_solution(&self, field: FluxSolution) -> Result<StokesSolution, StokesBoundaryError> {
        let pressure = self.system.pressure_in_gauges(&field.x, &self.has_traction);
        if !pressure.iter().all(|value| value.is_finite()) {
            return Err(StokesBoundaryError::InvalidInput(
                "non-finite reported pressure",
            ));
        }
        let reactions =
            self.reactions
                .evaluate(self.system, &field, &pressure, &self.has_traction)?;
        Ok(StokesSolution {
            field,
            pressure,
            pressure_reference: self.reference,
            pressure_reference_per_cell: (0..self.system.mesh.tris.len())
                .map(|cell| {
                    self.pressure_reference_at(cell)
                        .expect("cell belongs to mesh")
                })
                .collect(),
            reactions,
        })
    }
}

impl<'m> FluxSystem<'m> {
    /// Assemble steady -div(2 nu D(u)) + grad(p) = f, div(u) = 0.
    ///
    /// Every boundary edge must have exactly one velocity or total-traction
    /// owner. At least one velocity edge must anchor each connected domain.
    /// Natural edges retain both normal BDM moments and every continuity row.
    /// If all edges prescribe velocity, compatibility is checked before the
    /// existing pressure pin is used. That case is numerically identical to
    /// assemble_variable_viscosity with constant params.nu and the same data.
    ///
    /// Viscosity is kinematic [m²/s], body force is acceleration [m/s²], and
    /// total traction is stress divided by the constant reference density.
    /// There is deliberately no convection, backflow or transient argument.
    ///
    /// # Errors
    /// Refuses invalid/incomplete edge ownership before invoking callbacks,
    /// unanchored domains, non-finite data, incompatible component-local
    /// all-velocity flux, invalid parameters or unrepresentable assembly.
    pub fn assemble_stokes_mixed<'a>(
        &'a self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Stokes parameters",
            ));
        }
        self.assemble_stokes_mixed_variable_viscosity(params, &|_, _| params.nu, force, boundaries)
    }

    /// Assemble steady variable-viscosity `-div(2 nu D(u)) + grad(p) = f`
    /// with individually owned velocity and total-traction boundaries.
    ///
    /// `viscosity(cell, point)` returns positive finite kinematic viscosity
    /// [m²/s]. It must be pure and use the supplied cell for one-sided values
    /// at material interfaces; `params.nu` is unused. Dynamic viscosity must
    /// be divided by one constant reference density, as must supplied physical
    /// body force and traction. The traction callback supplies the TOTAL
    /// `(2 nu D(u) - p I)n`; it is not multiplied by viscosity again.
    ///
    /// This uses the same symmetric cell stresses, one-sided face traces and
    /// harmonic viscosity/normal-height SIP weighting as
    /// [`Self::assemble_variable_viscosity`]. Natural boundaries preserve
    /// their pressure datum and every continuity equation. There is no
    /// convection, transient or implicit rheological update in this entry.
    ///
    /// # Errors
    /// Refuses the same ownership, flux, geometry and solver failures as
    /// [`Self::assemble_stokes_mixed`]. The first invalid viscosity sample is
    /// retained with its cell and position, further viscosity queries stop,
    /// and no assembled operator or candidate flow is returned.
    pub fn assemble_stokes_mixed_variable_viscosity<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        let admitted = self.admit_stokes_problem(params, boundaries)?;
        self.assemble_stokes_admitted_viscosity(params, viscosity, force, &admitted)
    }

    fn assemble_stokes_admitted_viscosity<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        admitted: &AdmittedBoundary,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.assemble_stokes_admitted_viscosity_and_drag(params, viscosity, force, admitted, None)
    }

    pub(super) fn assemble_stokes_admitted_viscosity_and_drag<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        admitted: &AdmittedBoundary,
        darcy: Option<&[super::porous::DarcyDrag]>,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.assemble_stokes_admitted_operator(
            params, viscosity, force, admitted, darcy, 0.0, None, None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_stokes_admitted_operator<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        admitted: &AdmittedBoundary,
        darcy: Option<&[super::porous::DarcyDrag]>,
        inv_dt: f64,
        previous_velocity: Option<&[f64]>,
        advecting_velocity: Option<&[f64]>,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        let bad_viscosity = Cell::new(None);
        let checked_viscosity = |cell, point| {
            if bad_viscosity.get().is_some() {
                return 1.0;
            }
            let value = viscosity(cell, point);
            if value.is_finite() && value > 0.0 {
                value
            } else {
                bad_viscosity.set(Some((cell, point, value)));
                // Only finish the private assembly. The refused coefficient
                // never becomes a published matrix or a linear solve.
                1.0
            }
        };
        let bad_force = Cell::new(false);
        let checked_force = |point| {
            let value = force(point);
            if value.iter().all(|entry| entry.is_finite()) {
                value
            } else {
                bad_force.set(true);
                [0.0; 2]
            }
        };
        let mut reactions = reaction::ReactionAssembly::new(self, admitted, advecting_velocity);
        let (matrix, rhs) = self.assemble_operator(
            params,
            Some(&checked_viscosity),
            &checked_force,
            &|_| [0.0; 2],
            advecting_velocity,
            inv_dt,
            previous_velocity,
            Some(admitted),
            darcy,
            Some(&mut reactions),
        );
        if let Some((cell, point, value)) = bad_viscosity.get() {
            return Err(StokesBoundaryError::InvalidViscosity { cell, point, value });
        }
        if bad_force.get() {
            return Err(StokesBoundaryError::NonFiniteForce);
        }
        if !rhs.iter().all(|value| value.is_finite())
            || (0..matrix.nrows())
                .any(|row| !matrix.row(row).1.iter().all(|value| value.is_finite()))
        {
            return Err(StokesBoundaryError::InvalidInput(
                "non-finite assembled operator or load",
            ));
        }
        Ok(StokesAssembly {
            system: self,
            matrix,
            rhs,
            params,
            reference: if admitted.has_traction.iter().all(|&value| value) {
                StokesPressureReference::TractionDatum
            } else if admitted.has_traction.iter().all(|&value| !value) {
                StokesPressureReference::MeanZero
            } else {
                StokesPressureReference::Componentwise
            },
            has_traction: admitted.has_traction.clone(),
            reactions: reactions.finish(),
        })
    }

    fn admit_stokes_problem(
        &self,
        params: FluxParams,
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        self.admit_stokes_problem_with_velocity_anchors(params, boundaries, &[])
    }

    pub(super) fn admit_stokes_problem_with_velocity_anchors(
        &self,
        params: FluxParams,
        boundaries: &[StokesBoundary<'_>],
        anchors_velocity: &[bool],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        if !(params.sigma.is_finite()
            && params.sigma > 0.0
            && params.tol.is_finite()
            && params.tol > 0.0)
            || params.max_iters == 0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Stokes parameters",
            ));
        }
        if self.mesh.tris.is_empty() {
            return Err(StokesBoundaryError::InvalidInput("empty fluid mesh"));
        }
        self.admit_stokes_boundaries_with_velocity_anchors(boundaries, anchors_velocity)
    }

    pub(super) fn admit_stokes_boundaries(
        &self,
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        self.admit_stokes_boundaries_with_velocity_anchors(boundaries, &[])
    }

    fn admit_stokes_boundaries_with_velocity_anchors(
        &self,
        boundaries: &[StokesBoundary<'_>],
        anchors_velocity: &[bool],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        let mesh = self.mesh;
        let mut owners = vec![None; mesh.edges.len()];
        let components = self.pressure_roots.len();
        let mut has_velocity = vec![false; components];
        let mut has_traction = vec![false; components];
        let mut boundary_count = vec![0_usize; components];
        // Complete the ownership audit before querying any user callback.
        for boundary in boundaries {
            let edge = match *boundary {
                StokesBoundary::Velocity { edge, .. }
                | StokesBoundary::TotalTraction { edge, .. } => edge,
            };
            if edge >= mesh.edges.len() {
                return Err(StokesBoundaryError::EdgeOutOfRange { edge });
            }
            if !mesh.is_boundary(edge) {
                return Err(StokesBoundaryError::InteriorEdge { edge });
            }
            if owners[edge].replace(*boundary).is_some() {
                return Err(StokesBoundaryError::DuplicateEdge { edge });
            }
            let component = self.pressure_component[mesh.edges[edge].tris.0];
            boundary_count[component] += 1;
            match boundary {
                StokesBoundary::Velocity { .. } => has_velocity[component] = true,
                StokesBoundary::TotalTraction { .. } => has_traction[component] = true,
            }
        }
        for (edge, owner) in owners.iter().enumerate() {
            if mesh.is_boundary(edge) && owner.is_none() {
                return Err(StokesBoundaryError::MissingEdge { edge });
            }
        }
        for (component, &velocity) in has_velocity.iter().enumerate() {
            if !velocity && !anchors_velocity.get(component).copied().unwrap_or(false) {
                return Err(if components == 1 {
                    StokesBoundaryError::NoVelocityBoundary
                } else {
                    StokesBoundaryError::UnanchoredComponent { component }
                });
            }
        }
        let mut samples = vec![None; mesh.edges.len()];
        let mut net_flux = vec![0.0_f64; components];
        let mut absolute_flux = vec![0.0_f64; components];
        for (edge, owner) in owners.into_iter().enumerate() {
            let Some(owner) = owner else { continue };
            let geometry = &mesh.edges[edge];
            let component = self.pressure_component[geometry.tris.0];
            let points = edge_gauss_pub(mesh.verts[geometry.verts.0], mesh.verts[geometry.verts.1]);
            let mut values = [[0.0; 2]; 2];
            for (q, (point, weight)) in points.into_iter().enumerate() {
                values[q] = match owner {
                    StokesBoundary::Velocity { value, .. } => value(point),
                    StokesBoundary::TotalTraction { value, .. } => value(point, geometry.normal),
                };
                if !values[q].iter().all(|value| value.is_finite()) {
                    return Err(StokesBoundaryError::NonFiniteBoundary { edge });
                }
                if matches!(owner, StokesBoundary::Velocity { .. }) {
                    let flux = weight
                        * (values[q][0] * geometry.normal[0] + values[q][1] * geometry.normal[1]);
                    net_flux[component] += flux;
                    absolute_flux[component] += flux.abs();
                }
            }
            samples[edge] = Some(match owner {
                StokesBoundary::Velocity { .. } => Samples::Velocity(values),
                StokesBoundary::TotalTraction { .. } => Samples::Traction(values),
            });
        }
        for component in 0..components {
            let net_flux = net_flux[component];
            let absolute_flux = absolute_flux[component];
            if !net_flux.is_finite() || !absolute_flux.is_finite() {
                return Err(StokesBoundaryError::InvalidInput("boundary flux overflow"));
            }
            if has_traction[component] {
                continue;
            }
            let tolerance = 256.0
                * f64::EPSILON
                * boundary_count[component] as f64
                * absolute_flux.max(f64::MIN_POSITIVE);
            if net_flux.abs() > tolerance {
                return Err(StokesBoundaryError::IncompatibleFlux {
                    component,
                    net_flux,
                    tolerance,
                });
            }
        }
        Ok(AdmittedBoundary {
            samples,
            has_traction,
            open_policies: None,
            reservoir_velocities: None,
        })
    }
}
