//! Steady Stokes with individually owned velocity and total-traction edges.
//!
//! The natural datum is t = (2 nu D(u) - p I)n, with kinematic pressure
//! and traction (both divided by the same constant reference density).
//! Leaving a boundary velocity unknown also retains its continuity equation:
//! a traction boundary fixes the pressure datum, so no pressure pin is added.

use super::{FluxParams, FluxSolution, FluxSystem};
use crate::bdm::edge_gauss_pub;
use fs_sparse::Csr;
use std::cell::Cell;

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
    /// All boundary velocities are prescribed; output has zero area mean.
    MeanZero,
    /// The total-traction data fixes the pressure datum; output is not shifted.
    /// This is not thermodynamic/EOS absolute pressure.
    TractionDatum,
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
    /// Prescribed velocities violate the all-Dirichlet net-flux compatibility.
    IncompatibleFlux {
        /// Signed outward integral of prescribed velocity.
        net_flux: f64,
        /// Roundoff allowance for the assembled boundary quadrature.
        tolerance: f64,
    },
    /// Pure traction leaves rigid velocity modes; no hidden grounding is added.
    NoVelocityBoundary,
    /// This entry point admits one connected fluid domain.
    DisconnectedMesh,
    /// The assembled linear system could not be factored.
    SingularSystem,
    /// The explicitly recomputed residual did not meet the requested tolerance.
    NotConverged {
        /// Linear iterations performed.
        iterations: usize,
        /// True Euclidean ||b-Ax|| / max(||b||, 1e-150).
        relative_residual: f64,
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
    pub(super) has_traction: bool,
}

impl AdmittedBoundary {
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

/// A steady linear operator bound to its mesh and pressure interpretation.
/// The read-only matrix/RHS can also feed the existing discrete transpose solve.
pub struct StokesAssembly<'a, 'm> {
    system: &'a FluxSystem<'m>,
    matrix: Csr,
    rhs: Vec<f64>,
    params: FluxParams,
    reference: StokesPressureReference,
}

/// The admitted discrete field together with pressure in its declared datum.
pub struct StokesSolution {
    /// Velocity DOFs followed by raw P0 pressure coefficients.
    pub field: FluxSolution,
    /// Per-cell kinematic pressure [m²/s²], in pressure_reference.
    pub pressure: Vec<f64>,
    /// Whether an output gauge shift was applied.
    pub pressure_reference: StokesPressureReference,
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

    /// Solve and require a true residual, without changing traction pressure.
    ///
    /// # Errors
    /// Refuses singular systems, non-finite results or exhausted convergence.
    pub fn solve(&self) -> Result<StokesSolution, StokesBoundaryError> {
        let mut field = self
            .system
            .solve_linear_result(&self.matrix, &self.rhs, self.params)?;
        let mut residual = vec![0.0; self.rhs.len()];
        self.matrix.spmv(&field.x, &mut residual);
        let mut error = 0.0_f64;
        let mut load = 0.0_f64;
        for (actual, target) in residual.into_iter().zip(&self.rhs) {
            error = error.hypot(actual - target);
            load = load.hypot(*target);
        }
        let relative_residual = error / load.max(1e-150);
        if !field.x.iter().all(|value| value.is_finite())
            || !relative_residual.is_finite()
            || relative_residual > self.params.tol
        {
            return Err(StokesBoundaryError::NotConverged {
                iterations: field.iters,
                relative_residual,
            });
        }
        field.rel_residual = relative_residual;
        let pressure = match self.reference {
            StokesPressureReference::MeanZero => self.system.pressure(&field.x),
            StokesPressureReference::TractionDatum => field.x[self.system.n_u..].to_vec(),
        };
        if !pressure.iter().all(|value| value.is_finite()) {
            return Err(StokesBoundaryError::InvalidInput("non-finite reported pressure"));
        }
        Ok(StokesSolution {
            field,
            pressure,
            pressure_reference: self.reference,
        })
    }
}

impl<'m> FluxSystem<'m> {
    /// Assemble steady -div(2 nu D(u)) + grad(p) = f, div(u) = 0.
    ///
    /// Every boundary edge must have exactly one velocity or total-traction
    /// owner. At least one velocity edge must anchor this connected domain.
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
    /// disconnected or unanchored domains, non-finite data, incompatible
    /// all-velocity flux, invalid parameters or unrepresentable assembly.
    pub fn assemble_stokes_mixed<'a>(
        &'a self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        if !(params.nu.is_finite() && params.nu > 0.0
            && params.sigma.is_finite() && params.sigma > 0.0
            && params.tol.is_finite() && params.tol > 0.0)
            || params.max_iters == 0
        {
            return Err(StokesBoundaryError::InvalidInput("invalid Stokes parameters"));
        }
        if self.mesh.tris.is_empty() {
            return Err(StokesBoundaryError::InvalidInput("empty fluid mesh"));
        }
        if self.n > 1500 && params.max_iters < 100 {
            return Err(StokesBoundaryError::InvalidInput(
                "iterative Stokes needs a budget of at least one 100-step GMRES cycle",
            ));
        }
        let admitted = self.admit_stokes_boundaries(boundaries)?;
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
        let viscosity = |_: usize, _: [f64; 2]| params.nu;
        let (matrix, rhs) = self.assemble_operator(
            params,
            Some(&viscosity),
            &checked_force,
            &|_| [0.0; 2],
            None,
            0.0,
            None,
            Some(&admitted),
        );
        if bad_force.get() {
            return Err(StokesBoundaryError::NonFiniteForce);
        }
        if !rhs.iter().all(|value| value.is_finite())
            || (0..matrix.nrows()).any(|row| {
                !matrix.row(row).1.iter().all(|value| value.is_finite())
            })
        {
            return Err(StokesBoundaryError::InvalidInput("non-finite assembled operator or load"));
        }
        Ok(StokesAssembly {
            system: self,
            matrix,
            rhs,
            params,
            reference: if admitted.has_traction {
                StokesPressureReference::TractionDatum
            } else {
                StokesPressureReference::MeanZero
            },
        })
    }

    fn admit_stokes_boundaries(
        &self,
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        let mesh = self.mesh;
        let mut owners = vec![None; mesh.edges.len()];
        let mut has_velocity = false;
        let mut has_traction = false;
        // Complete the ownership audit before querying any user callback.
        for boundary in boundaries {
            let edge = match *boundary {
                StokesBoundary::Velocity { edge, .. } => {
                    has_velocity = true;
                    edge
                }
                StokesBoundary::TotalTraction { edge, .. } => {
                    has_traction = true;
                    edge
                }
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
        }
        for (edge, owner) in owners.iter().enumerate() {
            if mesh.is_boundary(edge) && owner.is_none() {
                return Err(StokesBoundaryError::MissingEdge { edge });
            }
        }
        if !has_velocity {
            return Err(StokesBoundaryError::NoVelocityBoundary);
        }
        let mut visited = vec![false; mesh.tris.len()];
        let mut queue = vec![0];
        visited[0] = true;
        while let Some(cell) = queue.pop() {
            for &(edge, _) in &mesh.tri_edges[cell] {
                let adjacent = mesh.edges[edge].tris;
                let other = if adjacent.0 == cell { adjacent.1 } else { adjacent.0 };
                if other != usize::MAX && !visited[other] {
                    visited[other] = true;
                    queue.push(other);
                }
            }
        }
        if visited.iter().any(|seen| !seen) {
            return Err(StokesBoundaryError::DisconnectedMesh);
        }
        let mut samples = vec![None; mesh.edges.len()];
        let mut net_flux = 0.0_f64;
        let mut absolute_flux = 0.0_f64;
        for (edge, owner) in owners.into_iter().enumerate() {
            let Some(owner) = owner else { continue };
            let geometry = &mesh.edges[edge];
            let points = edge_gauss_pub(
                mesh.verts[geometry.verts.0],
                mesh.verts[geometry.verts.1],
            );
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
                    net_flux += flux;
                    absolute_flux += flux.abs();
                }
            }
            samples[edge] = Some(match owner {
                StokesBoundary::Velocity { .. } => Samples::Velocity(values),
                StokesBoundary::TotalTraction { .. } => Samples::Traction(values),
            });
        }
        if !net_flux.is_finite() || !absolute_flux.is_finite() {
            return Err(StokesBoundaryError::InvalidInput("boundary flux overflow"));
        }
        if !has_traction {
            let tolerance = 256.0
                * f64::EPSILON
                * boundaries.len() as f64
                * absolute_flux.max(f64::MIN_POSITIVE);
            if net_flux.abs() > tolerance {
                return Err(StokesBoundaryError::IncompatibleFlux { net_flux, tolerance });
            }
        }
        Ok(AdmittedBoundary { samples, has_traction })
    }
}
