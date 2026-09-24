//! Real steady-conduction consumer of the tetrahedral continuum verifier.
//!
//! This adapter consumes `fs_conduction::ConductionProblem` without inventing
//! another thermal solver or replacing its material/boundary model. The primal
//! and unit-source homogeneous-boundary dual both use the existing FEM solver;
//! the verifier independently encloses their discretization AND algebraic errors.
//! An existing P1 field can also be bounded without running another primal solve.
//!
//! Admitted class: linear SPD tensor conductivity (possibly assigned per element),
//! element-constant sources, face-constant Neumann flux and Robin coefficient,
//! and affine Dirichlet/reference data. Unsupported inputs refuse before solves;
//! no averaging, tensor isotropization or frozen nonlinear coefficient can mint
//! a bound. Bounded-temperature material tables refuse because this adapter has
//! no maximum-temperature argument establishing their continuum validity.
//!
//! The interval is conditional on the exact declared polyhedral domain and
//! nominal coefficients, not a CAD, model-discrepancy or material-uncertainty
//! certificate. It is a VOLUME MEAN, never a nodal/point maximum. No evidence
//! colour or capability maturity is promoted. As in `fs-conduction::solve`, the
//! boundary partition must have been built for this same mesh.

use std::collections::BTreeSet;

use fs_conduction::{
    ConductionError, ConductionProblem, ConductionSolution, ScalarField, SolveConfig,
    TemperatureSpan, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_exec::Cx;

use crate::tet::{self, BoundaryCondition, BoundaryFace, ConductivityTensor, FluxBudget,
    MeanBound, TetError, TensorTetProblem};

/// Solvers retain their own tolerances; these are not substituted for an error bound.
#[derive(Debug, Clone, Default)]
pub struct MeanSolveConfig {
    pub primal: SolveConfig,
    pub dual: SolveConfig,
    pub flux: FluxBudget,
}

/// An actual FEM temperature field, actual unit-source dual, and continuum mean enclosure.
#[derive(Debug)]
pub struct MeanTemperatureSolution {
    pub primal: ConductionSolution,
    pub dual: ConductionSolution,
    pub bound: MeanBound,
}

/// A bound on an existing field with an actual same-operator unit-source dual.
/// No primal solve is performed or claimed. The supplied field may be an
/// incomplete iterate; its algebraic error stays in the bound and correction.
#[derive(Debug)]
pub struct MeanFieldBound {
    pub dual: ConductionSolution,
    pub bound: MeanBound,
}

/// Original solver refusals are preserved; an unsupported model yields no bound.
#[derive(Debug)]
pub enum ConductionBoundError {
    Conduction(ConductionError),
    Verification(TetError),
}
impl From<ConductionError> for ConductionBoundError {
    fn from(error: ConductionError) -> Self { Self::Conduction(error) }
}
impl From<TetError> for ConductionBoundError {
    fn from(error: TetError) -> Self { Self::Verification(error) }
}
impl std::fmt::Display for ConductionBoundError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Conduction(error) => write!(f, "conduction mean solve: {error}"),
            Self::Verification(error) => write!(f, "conduction mean bound: {error}"),
        }
    }
}
impl std::error::Error for ConductionBoundError {}

type Result<T> = std::result::Result<T, ConductionBoundError>;

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| TetError::Cancelled)?;
    Ok(())
}
fn unsupported(what: &'static str) -> ConductionBoundError {
    TetError::Unsupported(what).into()
}

struct Admitted {
    tets: Vec<[usize; 4]>,
    conductivity: Vec<ConductivityTensor>,
    source: Vec<f64>,
    boundary: Vec<BoundaryFace>,
}
impl Admitted {
    fn problem<'a>(&'a self, vertices: &'a [[f64; 3]]) -> TensorTetProblem<'a> {
        TensorTetProblem { vertices, tets: &self.tets, conductivity: &self.conductivity,
            source: &self.source, boundary: &self.boundary }
    }
}

fn admit(cx: &Cx<'_>, problem: ConductionProblem<'_>, budget: FluxBudget) -> Result<Admitted> {
    poll(cx)?;
    let mesh = problem.mesh;
    let n = mesh.vertex_count();
    let ne = mesh.element_count();
    if ne > budget.max_cells || n > budget.max_cells.saturating_mul(4)
        || budget.max_iterations > 1_000_000 {
        return Err(TetError::Budget.into());
    }
    problem.source.validate("bounded volumetric source", n)?;
    if let Some(materials) = problem.element_materials { materials.validate_for(mesh)?; }
    for condition in problem.boundary.conditions() { poll(cx)?; condition.validate(n)?; }
    let mut tets = Vec::with_capacity(ne);
    let mut conductivity = Vec::with_capacity(ne);
    let mut source = Vec::with_capacity(ne);
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        poll(cx)?;
        let model = match problem.element_materials {
            Some(materials) => materials.model_for(e)?,
            None => problem.material,
        };
        if model.is_temperature_dependent() {
            return Err(unsupported("temperature-dependent conductivity needs a nonlinear residual bound"));
        }
        if !matches!(model.temperature_span(), TemperatureSpan::Unbounded) {
            return Err(unsupported("bounded material validity needs a separate continuum temperature-range proof"));
        }
        // Preserve every coefficient and its mesh-frame orientation. Outward
        // SPD/inverse admission below is stricter than an ordinary Cholesky test.
        let tensor = model.tensor_at(0.0)?;
        let vertices = tet.map(|i| i as usize);
        let values = vertices.map(|i| problem.source.at(i));
        if values.iter().any(|&v| v != values[0]) {
            return Err(unsupported("nonconstant element source requires a data-oscillation or higher-order flux term"));
        }
        tets.push(vertices);
        conductivity.push(tensor);
        source.push(values[0]);
    }
    let mut boundary = Vec::with_capacity(mesh.boundary().len());
    for (slot, face) in mesh.boundary().iter().enumerate() {
        poll(cx)?;
        let vertices = face.vertices.map(|v| v as usize);
        let condition = match problem.boundary.condition_for(slot) {
            Some(ThermalBc::Dirichlet { temperature }) =>
                BoundaryCondition::Dirichlet(vertices.map(|v| temperature.at(v))),
            Some(ThermalBc::Neumann { outward_flux }) => {
                let flux = vertices.map(|v| outward_flux.at(v));
                if flux.iter().any(|&v| v != flux[0]) {
                    return Err(unsupported("nonconstant face Neumann flux requires a higher-order normal trace"));
                }
                BoundaryCondition::Neumann(flux[0])
            }
            Some(ThermalBc::Robin { htc, t_ref }) => {
                let h = vertices.map(|v| htc.at(v));
                if h.iter().any(|&v| v != h[0]) {
                    return Err(unsupported("nonconstant face Robin coefficient cannot be replaced by its mean"));
                }
                BoundaryCondition::Robin { h: h[0], reference: vertices.map(|v| t_ref.at(v)) }
            }
            // This state is produced only by an explicit adiabatic remainder.
            None => BoundaryCondition::Neumann(0.0),
        };
        boundary.push(BoundaryFace { vertices, condition });
    }
    let admitted = Admitted { tets, conductivity, source, boundary };
    admitted.problem(mesh.positions()).validate_conductivity(budget, || cx.checkpoint().is_ok())?;
    Ok(admitted)
}

fn dual_boundary(cx: &Cx<'_>, problem: ConductionProblem<'_>) -> Result<ThermalBoundary> {
    let mut faces = vec![BTreeSet::new(); problem.boundary.conditions().len()];
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        poll(cx)?;
        if let Some(region) = problem.boundary.region_for(slot) { faces[region].insert(face.face); }
    }
    let mut builder = ThermalBoundaryBuilder::new(problem.mesh);
    for (region, condition) in problem.boundary.conditions().iter().enumerate() {
        poll(cx)?;
        let dual = match condition {
            ThermalBc::Dirichlet { .. } => ThermalBc::dirichlet(0.0)?,
            ThermalBc::Neumann { .. } => ThermalBc::adiabatic(),
            ThermalBc::Robin { htc, .. } => ThermalBc::Robin {
                htc: htc.clone(), t_ref: ScalarField::Uniform(0.0),
            },
        };
        builder = builder.region(&problem.boundary.region_names()[region],
            |face| faces[region].contains(&face.face), dual)?;
    }
    if problem.boundary.adiabatic_remainder_faces() != 0 {
        builder = builder.adiabatic_remainder();
    }
    let dual = builder.finish()?;
    if dual.adiabatic_remainder_faces() != problem.boundary.adiabatic_remainder_faces() {
        return Err(TetError::Invalid("dual boundary partition does not reproduce the primal").into());
    }
    Ok(dual)
}

fn bound_field(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, admitted: &Admitted,
    temperature: &[f64], dual_config: SolveConfig, flux: FluxBudget,
) -> Result<MeanFieldBound> {
    poll(cx)?;
    if temperature.len() != problem.mesh.vertex_count() {
        return Err(TetError::Invalid("candidate temperature length").into());
    }
    for chunk in temperature.chunks(256) {
        poll(cx)?;
        if chunk.iter().any(|t| !t.is_finite()) {
            return Err(TetError::Invalid("finite candidate temperature required").into());
        }
    }
    for &(vertex, value) in problem.boundary.dirichlet() {
        poll(cx)?;
        if temperature.get(vertex) != Some(&value) {
            return Err(TetError::Invalid("candidate does not match prescribed Dirichlet trace").into());
        }
    }
    let boundary = dual_boundary(cx, problem)?;
    let unit_source = ScalarField::Uniform(1.0);
    let dual = fs_conduction::solve(cx, ConductionProblem {
        boundary: &boundary, source: &unit_source, ..problem
    }, dual_config)?;
    poll(cx)?;
    let bound = tet::tensor_mean_bound(&admitted.problem(problem.mesh.positions()),
        temperature, &dual.temperature, flux, || cx.checkpoint().is_ok())?;
    Ok(MeanFieldBound { dual, bound })
}

/// Bound a supplied P1 temperature field without repeating its primal solve.
/// The field need not satisfy the discrete equations: the complete residual
/// correction and equilibrated majorant include its remaining algebraic error.
/// Length, finite values and the exact Dirichlet trace are checked before the
/// dual solve. The caller binds the nominal problem, not a trusted solve report.
/// No mesh, contact, nonlinear or uncertainty assumption is relaxed.
pub fn bound_temperature_mean(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
    dual_config: SolveConfig, flux: FluxBudget,
) -> Result<MeanFieldBound> {
    let admitted = admit(cx, problem, flux)?;
    bound_field(cx, problem, &admitted, temperature, dual_config, flux)
}

/// Solve the original thermal problem and its unit-volume-source adjoint, then
/// enclose the exact-domain average temperature. Both solve reports and their
/// material provenance remain available. The bound applies to the nominal
/// admitted LINEAR tensor operator; it cannot certify a coupled nonlinear loop.
///
/// The caller supplies independent primal, dual and flux budgets. Cancellation
/// uses the same Cx for both existing FEM solves and the flux/goal integration.
/// An exhausted solve or verifier returns an error, never a partial certificate.
pub fn solve_with_mean_bound(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, config: MeanSolveConfig,
) -> Result<MeanTemperatureSolution> {
    let admitted = admit(cx, problem, config.flux)?;
    let primal = fs_conduction::solve(cx, problem, config.primal)?;
    let result = bound_field(cx, problem, &admitted, &primal.temperature, config.dual, config.flux)?;
    Ok(MeanTemperatureSolution { primal, dual: result.dual, bound: result.bound })
}

#[cfg(test)]
mod tests;
#[cfg(test)]
#[path = "conduction/tensor_tests.rs"]
mod tensor_tests;
