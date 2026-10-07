//! Steady conjugate heat transfer (CHT) on voxel domains: D3Q19 airflow plus
//! one conservative finite-volume energy equation over fluid AND solid cells
//! (bead frankensim-rc-root-q61wp.34, "fluids for the wedge").
//!
//! # What this computes
//!
//! A forced-convection problem with one-way coupling: the incompressible
//! steady flow does not depend on temperature (no buoyancy), so the pipeline
//! is
//!
//! 1. **Flow.** [`lbm_duct_flow`] drives the existing [`crate::BoundaryGrid3`]
//!    (velocity inlet on `x-min`, pressure outlet on `x-max`, halfway
//!    bounce-back walls and voxel obstacles) to a steady state measured by
//!    the relative velocity change between checks. [`simple_flow`] solves
//!    the same steady flow by finite volumes on the staggered faces (SIMPLEC)
//!    with no lattice constraint, returning divergence-free face fluxes
//!    directly; [`fv_natural_convection`] couples it to the energy equation
//!    through the Boussinesq force, including pressure openings. Any other producer of
//!    cell or face velocities (an analytic profile, another solver) enters
//!    through [`FlowField::from_cell_velocities`] or
//!    [`FlowField::from_face_velocity`].
//! 2. **Projection.** Face fluxes interpolated from cell velocities are made
//!    discretely divergence-free by one SPD pressure-correction solve on the
//!    fluid cells ([`FlowField::from_cell_velocities`]); the returned
//!    [`ProjectionReport`] retains the divergence before and after and the
//!    size of the correction. Fluxes across fluid/solid faces are exactly
//!    zero by construction.
//! 3. **Energy.** [`solve_energy`] assembles
//!    `div(rho c_p u T) - div(k grad T) = q'''` on every cell. Diffusion uses
//!    the harmonic-mean face conductivity, which is the exact two-point flux
//!    for piecewise-constant conductivity meeting at a cell face, so the
//!    fluid/solid interface conditions (continuity of temperature AND normal
//!    heat flux) hold by construction even at the ~10^4 aluminium/air
//!    contrast that defeats explicit thermal-LBM conjugate schemes. Convection
//!    uses Patankar's power-law scheme (or first-order upwind) in its
//!    conservative form, so the discrete total energy flux telescopes exactly:
//!    the [`EnergyBalance`] closes to solver tolerance whatever the residual
//!    divergence of the supplied flux field.
//!
//! The linear system is an M-matrix (non-positive off-diagonals, weakly
//! diagonally dominant for a divergence-free flux field), so the discrete
//! solution obeys a maximum principle. It is solved by ILU(0)-preconditioned
//! BiCGStab on the Jacobi-row-scaled system; the reported residual is the
//! recomputed true residual, never the recurrence estimate.
//!
//! # Units
//!
//! SI throughout: lengths in metres, velocities m/s, volumetric fluxes
//! m^3/s, conductivity W/(m K), power W, temperature K. Lattice quantities
//! appear only inside [`lbm_duct_flow`] and its report.
//!
//! # Determinism
//!
//! Single-threaded, fixed traversal order (`x` fastest, then `y`, then `z`;
//! faces in [`Face3::ALL`] order), deterministic COO assembly and sequential
//! reductions. Results are bit-reproducible for identical inputs on the same
//! ISA/toolchain profile. No cross-ISA identity is claimed.
//!
//! # Cancellation
//!
//! Every public driver takes a [`CancelGate`] and polls it between LBM step
//! batches and between Krylov iterations; a tripped gate returns
//! [`ChtError::Cancelled`] and publishes nothing.
//!
//! # No-claim boundaries
//!
//! - Steady, constant-property convection: forced, or Boussinesq natural
//!   and mixed; turbulence only through the algebraic LVEL closure
//!   ([`turbulence`]); radiation only as surface emission to the
//!   surroundings; no temperature-dependent properties. Fluid properties are
//!   frozen at the declared state.
//! - Voxel (staircase) geometry at the declared resolution; curved walls are
//!   represented to `O(dx)`. No mesh-convergence claim is made by a single
//!   run; refinement ladders are the caller's evidence.
//! - The power-law convection scheme is first-order where the cell Péclet
//!   number exceeds about 10; results are Estimated numerical evidence, not
//!   certified enclosures.
//! - The LBM flow is a weakly compressible approximation; its divergence is
//!   removed by the projection, whose correction magnitude is retained.

use fs_exec::CancelGate;

use crate::d3q19::Face3;

mod buoyant;
mod domain;
mod energy;
mod flow;
mod krylov;
mod natural;
mod radiation;
mod simple;
mod transient;
pub mod turbulence;
mod unsteady;

pub use buoyant::{
    FvBuoyancyConfig, FvNaturalConvection, FvNaturalConvectionReport, fv_natural_convection,
};
pub use domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
pub use energy::{
    CellSink, CompactComponent, ContactResistance, ConvectionScheme, EnergyBalance, EnergyConfig,
    EnergyReport, EnergySolution, JunctionSolution, ThermalFace, ThermalSetup, solve_energy,
};
pub use flow::{
    FlowFace, FlowField, LbmCollisionChoice, LbmFlow, LbmFlowConfig, LbmFlowReport,
    ProjectionReport, lbm_duct_flow,
};
pub use natural::{BuoyancyConfig, NaturalConvection, NaturalConvectionReport, natural_convection};
pub use radiation::{
    ExposedFace, RadiationConfig, RadiationReport, STEFAN_BOLTZMANN, SurfaceExchange,
    escape_factors, radiated_power, radiative_sinks, solve_energy_radiating,
};
pub use simple::{
    AMG_REBUILD_SWEEPS, FacePatch, FanCurve, FanInlet, FlowResistance, FvBoundary, FvFlow,
    InternalFan, PressureSolver, SimpleConfig, SimpleReport, TimeScheme, Turbulence, simple_flow,
};
pub use transient::{TransientConfig, TransientRecord, TransientSolution, march_energy};
pub use unsteady::{
    Boussinesq, ConjugateMarch, FlowStepRecord, UnsteadyConfig, UnsteadyFlow, march_conjugate,
    simple_unsteady,
};

/// Semantics version of the conjugate pipeline: covers voxel indexing, face
/// flux layout, projection, energy discretization, solver policy, and the
/// LBM driver's steady criterion. Bump on any change that can move results.
pub const CHT_SEMANTICS_VERSION: u32 = 1;

/// Structured refusal from any stage of the conjugate pipeline.
#[derive(Debug, Clone, PartialEq)]
pub enum ChtError {
    /// A domain dimension is zero, a cell count overflows, or the spacing is
    /// not finite and positive.
    InvalidDomain {
        /// Human-readable reason.
        reason: String,
    },
    /// A material, fluid property, boundary value, or source is not finite or
    /// outside its physical range.
    InvalidInput {
        /// Which input was refused.
        field: &'static str,
        /// Human-readable reason.
        reason: String,
    },
    /// A voxel names a solid material that the material table does not have.
    UnknownMaterial {
        /// Linear cell index.
        cell: usize,
        /// The missing material index.
        material: u16,
    },
    /// The flow field carries volumetric flux through a domain face whose
    /// thermal rule cannot account for advected energy.
    FlowThroughClosedFace {
        /// The offending domain face.
        face: Face3,
        /// Net outward volumetric flux through that face, m^3/s.
        net_flux_m3_s: f64,
    },
    /// The lattice cannot represent the declared flow at this resolution.
    LatticeResolution {
        /// Relaxation time the declared state requires.
        tau: f64,
        /// Cell Reynolds number `U dx / nu`.
        cell_reynolds: f64,
        /// Human-readable remedy.
        remedy: String,
    },
    /// The LBM state became non-finite.
    FlowDiverged {
        /// Step at which the non-finite state was observed.
        step: usize,
    },
    /// The LBM run exhausted its step budget before the steady criterion.
    FlowNotSteady {
        /// Steps taken.
        steps: usize,
        /// Last relative velocity change between checks.
        last_change: f64,
        /// Declared steady tolerance.
        tolerance: f64,
    },
    /// The Krylov solve did not reach the requested true residual.
    SolverNotConverged {
        /// Which system (`"projection"` or `"energy"`).
        system: &'static str,
        /// Iterations performed.
        iterations: usize,
        /// Recomputed relative residual of the row-scaled system.
        relative_residual: f64,
        /// Requested tolerance.
        tolerance: f64,
    },
    /// ILU(0) met a zero pivot (should not occur for the admitted M-matrices).
    PreconditionerBreakdown {
        /// Which system.
        system: &'static str,
        /// Failing row.
        row: usize,
    },
    /// The tile pool refused a pass (worker fault or admission).
    Executor {
        /// Rendered pool diagnostic.
        detail: String,
    },
    /// The cancel gate tripped; nothing was published.
    Cancelled,
}

impl core::fmt::Display for ChtError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDomain { reason } => write!(f, "invalid CHT domain: {reason}"),
            Self::InvalidInput { field, reason } => {
                write!(f, "invalid CHT input `{field}`: {reason}")
            }
            Self::UnknownMaterial { cell, material } => {
                write!(
                    f,
                    "cell {cell} names solid material {material}, which is not declared"
                )
            }
            Self::FlowThroughClosedFace {
                face,
                net_flux_m3_s,
            } => write!(
                f,
                "flow crosses closed thermal face {face:?} (net outward {net_flux_m3_s:e} m^3/s); declare Inflow/Outflow there"
            ),
            Self::LatticeResolution {
                tau,
                cell_reynolds,
                remedy,
            } => write!(
                f,
                "lattice cannot represent the flow: tau = {tau}, cell Reynolds = {cell_reynolds}; {remedy}"
            ),
            Self::FlowDiverged { step } => write!(f, "LBM flow became non-finite at step {step}"),
            Self::FlowNotSteady {
                steps,
                last_change,
                tolerance,
            } => write!(
                f,
                "flow not steady after {steps} steps/iterations (residual or change {last_change:e} > {tolerance:e})"
            ),
            Self::SolverNotConverged {
                system,
                iterations,
                relative_residual,
                tolerance,
            } => write!(
                f,
                "{system} solve stopped after {iterations} iterations at relative residual {relative_residual:e} (tolerance {tolerance:e})"
            ),
            Self::PreconditionerBreakdown { system, row } => {
                write!(f, "{system} ILU(0) zero pivot at row {row}")
            }
            Self::Executor { detail } => write!(f, "tile pool refused an LBM pass: {detail}"),
            Self::Cancelled => write!(f, "conjugate heat-transfer run cancelled"),
        }
    }
}

impl std::error::Error for ChtError {}

fn poll(gate: &CancelGate) -> Result<(), ChtError> {
    if gate.is_requested() {
        Err(ChtError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite_positive(field: &'static str, value: f64) -> Result<(), ChtError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(ChtError::InvalidInput {
            field,
            reason: format!("must be finite and positive, got {value}"),
        })
    }
}

fn finite(field: &'static str, value: f64) -> Result<(), ChtError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ChtError::InvalidInput {
            field,
            reason: format!("must be finite, got {value}"),
        })
    }
}
