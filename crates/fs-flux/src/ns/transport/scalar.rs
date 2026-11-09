//! Passive scalar transport on the original BDM1 velocity and triangle mesh.
//!
//! The equation is dc/dt + div(u c - alpha grad(c)) = source. P1 discontinuous
//! scalars use conservative upwind advection and symmetric interior-penalty
//! diffusion, including both streams of a reversing affine normal trace.
//! Three values per triangle follow that triangle's original vertex order.
//! No centroid-orthogonality assumption or two-point diffusion is introduced.
//!
//! Temperature uses alpha=k/(rho Cp), source=Q/(rho Cp), with explicitly fixed
//! rho and Cp. Multiply area-integrated temperature budgets by rho Cp and
//! physical depth for energy. This passive equation supplies no EOS, variable
//! heat capacity, buoyancy, phase change or automatic viscous/SGS heating.
//! DG/BDF transport is not a positivity or maximum-principle guarantee.

use crate::TriMesh;
use crate::bdm::tri_quad;
use crate::ns::FluxSystem;
use crate::ns::boundary::bdf::StokesBdfStep;
use fs_math::dd::Dd;
use fs_solver::{FgmresState, SolverRunProgress, op::CsrOp};
use fs_sparse::precond::{Ilu0, IluControlledError, Precond, ilu0_controlled};
use fs_sparse::{Coo, Csr};

pub(in crate::ns) mod assembly;
/// Opt-in conservative scalar bounds using original-operator flux correction.
pub mod bounded;

/// Each original boundary edge needs exactly one owner. Callbacks must be pure.
#[derive(Clone, Copy)]
pub enum ScalarBoundary<'a> {
    /// Weak diffusion value and incoming advective donor. Outflow uses the
    /// actual interior scalar. Both use the same two-point affine boundary trace.
    Value {
        edge: usize,
        value: &'a dyn Fn([f64; 2]) -> f64,
    },
    /// Outward diffusive flux -alpha grad(c).n, in scalar*m/s. Zero is insulating.
    /// Every incoming stream needs an explicit donor even at zero NET flow.
    Flux {
        edge: usize,
        outward: &'a dyn Fn([f64; 2]) -> f64,
        inflow: Option<&'a dyn Fn([f64; 2]) -> f64>,
    },
}

/// Geometry/operator admission. The penalty uses harmonic alpha/h with the
/// original triangle normal heights; a value at least 12 is required for P1.
#[derive(Clone, Copy, Debug)]
pub struct ScalarTransportControls {
    pub penalty: f64,
    /// Numerical continuity allowance in [0,1e-12], never a flow projection.
    pub relative_volume_divergence_tolerance: f64,
}
impl Default for ScalarTransportControls {
    fn default() -> Self {
        Self {
            penalty: 12.0,
            relative_volume_divergence_tolerance: 1e-12,
        }
    }
}

/// A bounded original sparse Krylov solve. Both original and row-scaled
/// residuals are independently recomputed before an endpoint is returned.
#[derive(Clone, Copy, Debug)]
pub struct ScalarSolveControls {
    pub relative_tolerance: f64,
    pub maximum_iterations: usize,
}
impl Default for ScalarSolveControls {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            maximum_iterations: 400,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ScalarTransportError {
    Invalid(&'static str),
    Volume(super::FluxTransportError),
    VolumeDivergence { cell: usize, relative_residual: f64 },
    MissingInflow { edge: usize },
    Arithmetic(&'static str),
    LinearSolve { iterations: usize, residual: f64 },
    Cancelled,
}
impl core::fmt::Display for ScalarTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "BDM scalar transport: {self:?}")
    }
}
impl std::error::Error for ScalarTransportError {}
impl From<super::FluxTransportError> for ScalarTransportError {
    fn from(error: super::FluxTransportError) -> Self {
        Self::Volume(error)
    }
}
type Result<T> = core::result::Result<T, ScalarTransportError>;

/// Numerical flux leaving one boundary, in scalar*m²/s per out-of-plane depth.
/// Diffusive value boundaries include the same consistent SIP penalty flux
/// used in the equations; it is not merely a sampled gradient diagnostic.
#[derive(Clone, Copy, Debug)]
pub struct ScalarBoundaryRate {
    pub edge: usize,
    pub outward_advection: f64,
    pub outward_diffusion: f64,
}

/// One accepted endpoint. Histories and the prepared operator remain unchanged.
pub struct ScalarTransportStep {
    pub values: Vec<f64>,
    pub iterations: usize,
    pub relative_residual: f64,
    /// Integral of the endpoint scalar [scalar*m²].
    pub integral: f64,
    /// Integral of the actual BDF derivative [scalar*m²/s]. For BDF2 this
    /// uses BOTH histories and is not (new-old)/dt.
    pub storage_rate: f64,
    pub source_rate: f64,
    pub boundaries: Vec<ScalarBoundaryRate>,
    /// Storage + outward advection + outward diffusion - source. Numerical
    /// roundoff/solve defect is reported, never repaired or deposited as heat.
    pub balance_residual: f64,
}

struct BoundaryFlux {
    edge: usize,
    cell: usize,
    advection: [f64; 3],
    diffusion: [f64; 3],
    advective_load: f64,
    diffusive_load: f64,
    // Complete native weak boundary load, resolved onto the SAME two original
    // Gauss values. Some SIP/trace-interpolation coefficients can be negative.
    reservoir_coefficients: [[f64; 2]; 3],
    reservoir_values: [f64; 2],
    prescribed_flux: bool,
}

/// Frozen endpoint velocity, diffusivity and boundary traces. Reprepare when
/// any changes. Assembly is O(cells+edges); storage and solves remain sparse.
pub struct PreparedScalarTransport<'m> {
    mesh: &'m TriMesh,
    operator: Csr,
    boundary_load: Vec<f64>,
    boundary_fluxes: Vec<BoundaryFlux>,
    diffusivity_m2_s: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(ScalarTransportError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64, name: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ScalarTransportError::Arithmetic(name))
    }
}
fn sum(values: impl Iterator<Item = f64>) -> f64 {
    values
        .fold(Dd::from_f64(0.0), |total, x| total + Dd::from_f64(x))
        .to_f64()
}

struct ScalarPreconditioner(Option<Ilu0>);
impl Precond for ScalarPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        if let Some(factor) = &self.0 {
            factor.apply(r, z);
        } else {
            z.copy_from_slice(r);
        }
    }
}

impl PreparedScalarTransport<'_> {
    #[must_use]
    pub fn diffusivity_m2_s(&self) -> &[f64] {
        &self.diffusivity_m2_s
    }

    /// BDF1 with current endpoint source [scalar/s]. No explicit diffusion CFL
    /// restriction is imposed; temporal accuracy still requires refinement.
    pub fn bdf1_step(
        &self,
        previous: &[f64],
        timestep_s: f64,
        source: &dyn Fn(usize, [f64; 2]) -> f64,
        controls: ScalarSolveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ScalarTransportStep> {
        self.bdf2_variable_step(previous, None, timestep_s, source, controls, cancelled)
    }

    /// Original native unequal-step BDF2 history; no older state is exactly
    /// BDF1 startup. The new/old interval ratio must lie in (0,2]. Velocity,
    /// diffusivity and boundary data are evaluated at the supplied endpoint;
    /// they are not extrapolated here. Sources are sampled on the original
    /// degree-four cell quadrature. Cancellation/refusal returns no partial step.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_variable_step(
        &self,
        previous: &[f64],
        older: Option<(&[f64], f64)>,
        timestep_s: f64,
        source: &dyn Fn(usize, [f64; 2]) -> f64,
        controls: ScalarSolveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScalarTransportStep> {
        poll(&mut cancelled)?;
        let n = self.operator.nrows();
        if previous.len() != n
            || previous.iter().any(|v| !v.is_finite())
            || older.is_some_and(|(v, _)| v.len() != n || v.iter().any(|x| !x.is_finite()))
            || !controls.relative_tolerance.is_finite()
            || 0.1 * controls.relative_tolerance <= 0.0
        {
            return Err(ScalarTransportError::Invalid(
                "finite original P1 histories and positive tolerance required",
            ));
        }
        let time = StokesBdfStep::new(timestep_s, older.map(|(_, dt)| dt))
            .map_err(ScalarTransportError::Invalid)?;
        let inv_dt = time.effective_timestep_s.recip();
        let history: Vec<_> = previous
            .iter()
            .enumerate()
            .map(|(i, &value)| time.history(value, older.map_or(value, |(v, _)| v[i])))
            .collect();
        if history.iter().any(|v| !v.is_finite()) {
            return Err(ScalarTransportError::Arithmetic("BDF scalar history"));
        }
        let mut matrix = Coo::new(n, n);
        for row in 0..n {
            poll(&mut cancelled)?;
            let (columns, values) = self.operator.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                matrix.push(row, column, value);
            }
        }
        let mut rhs = self.boundary_load.clone();
        let mut sources = Vec::with_capacity(self.mesh.tris.len());
        for (cell, tri) in self.mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let area = self.mesh.areas[cell];
            for i in 0..3 {
                for j in 0..3 {
                    let mass = finite(
                        area * inv_dt / 12.0 * if i == j { 2.0 } else { 1.0 },
                        "scalar mass",
                    )?;
                    matrix.push(3 * cell + i, 3 * cell + j, mass);
                    rhs[3 * cell + i] += mass * history[3 * cell + j];
                }
            }
            let points = tri.map(|vertex| self.mesh.verts[vertex]);
            let basis = assembly::ScalarBasis::new(points, area);
            let mut integrated = Dd::from_f64(0.0);
            for (point, weight) in tri_quad(points, area) {
                let value = finite(source(cell, point), "scalar source")?;
                integrated = integrated + Dd::from_f64(weight * value);
                let phi = basis.value(point);
                for i in 0..3 {
                    rhs[3 * cell + i] += weight * value * phi[i];
                }
            }
            sources.push(finite(integrated.to_f64(), "integrated scalar source")?);
        }
        let matrix = matrix.assemble();
        let mut scaled = Coo::new(n, n);
        let mut scaled_rhs = Vec::with_capacity(n);
        for row in 0..n {
            poll(&mut cancelled)?;
            let (columns, values) = matrix.row(row);
            let scale = values.iter().fold(0.0_f64, |s, &v| s.max(v.abs()));
            if !scale.is_finite() || scale <= 0.0 || values.iter().any(|v| !v.is_finite()) {
                return Err(ScalarTransportError::Arithmetic("scalar operator"));
            }
            for (&column, &value) in columns.iter().zip(values) {
                scaled.push(row, column, value / scale);
            }
            scaled_rhs.push(finite(rhs[row] / scale, "scalar load")?);
        }
        let scaled = scaled.assemble();
        let factor = match ilu0_controlled(&scaled, || poll(&mut cancelled)) {
            Ok(factor) => Some(factor),
            Err(IluControlledError::Checkpoint(error)) => return Err(error),
            Err(IluControlledError::Breakdown(_)) => None,
            Err(IluControlledError::Capacity) => {
                return Err(ScalarTransportError::Invalid("scalar factor capacity"));
            }
        };
        let preconditioner = ScalarPreconditioner(factor);
        let operator = CsrOp::general(scaled);
        let mut state = FgmresState::new(&scaled_rhs, 40.min(n).max(1));
        let mut work = vec![0.0; n];
        let residual = |state: &[f64], work: &mut [f64]| {
            let (original, scaled) = crate::ns::sparse::residuals(&matrix, &rhs, state, work);
            if original.is_finite() && scaled.is_finite() {
                original.max(scaled)
            } else {
                f64::INFINITY
            }
        };
        let mut defect = residual(&state.x, &mut work);
        let mut inner_tolerance = 0.1 * controls.relative_tolerance;
        while defect > controls.relative_tolerance && state.iters < controls.maximum_iterations {
            poll(&mut cancelled)?;
            state.restart = 40
                .min(n)
                .min(controls.maximum_iterations - state.iters)
                .max(1);
            let before = state.iters;
            let report = state.run_controlled(
                &operator,
                &preconditioner,
                &scaled_rhs,
                inner_tolerance,
                1,
                &mut cancelled,
            );
            match report.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => return Err(ScalarTransportError::Cancelled),
                _ => return Err(ScalarTransportError::Arithmetic("scalar solver callback")),
            }
            defect = residual(&state.x, &mut work);
            if state.iters == before {
                break;
            }
            if report.report.converged_euclidean() && defect > controls.relative_tolerance {
                inner_tolerance = inner_tolerance
                    .min(0.5 * state.rel_residual() * controls.relative_tolerance / defect);
            }
        }
        if !defect.is_finite()
            || defect > controls.relative_tolerance
            || state.x.iter().any(|x| !x.is_finite())
        {
            return Err(ScalarTransportError::LinearSolve {
                iterations: state.iters,
                residual: defect,
            });
        }
        let mut boundaries = Vec::with_capacity(self.boundary_fluxes.len());
        for face in &self.boundary_fluxes {
            poll(&mut cancelled)?;
            let values = &state.x[3 * face.cell..3 * face.cell + 3];
            boundaries.push(ScalarBoundaryRate {
                edge: face.edge,
                outward_advection: finite(
                    sum(face.advection.iter().zip(values).map(|(a, c)| a * c))
                        + face.advective_load,
                    "scalar advection rate",
                )?,
                outward_diffusion: finite(
                    sum(face.diffusion.iter().zip(values).map(|(a, c)| a * c))
                        + face.diffusive_load,
                    "scalar diffusion rate",
                )?,
            });
        }
        let integral = finite(
            sum(self.mesh.areas.iter().enumerate().map(|(cell, area)| {
                area / 3.0 * sum(state.x[3 * cell..3 * cell + 3].iter().copied())
            })),
            "scalar integral",
        )?;
        let storage_rate = finite(
            sum(self.mesh.areas.iter().enumerate().map(|(cell, area)| {
                area / 3.0
                    * sum((0..3).map(|i| (state.x[3 * cell + i] - history[3 * cell + i]) * inv_dt))
            })),
            "scalar storage rate",
        )?;
        let source_rate = finite(sum(sources.into_iter()), "scalar source rate")?;
        let balance_residual = finite(
            sum(std::iter::once(storage_rate)
                .chain(std::iter::once(-source_rate))
                .chain(
                    boundaries
                        .iter()
                        .flat_map(|face| [face.outward_advection, face.outward_diffusion]),
                )),
            "scalar balance",
        )?;
        poll(&mut cancelled)?;
        Ok(ScalarTransportStep {
            values: state.x,
            iterations: state.iters,
            relative_residual: defect,
            integral,
            storage_rate,
            source_rate,
            boundaries,
            balance_residual,
        })
    }
}
