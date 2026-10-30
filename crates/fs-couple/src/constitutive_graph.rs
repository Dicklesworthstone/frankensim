//! Close material-law feedback outside the local constitutive DAG.
//!
//! Equations are data: dimensionally checked linear balances of native graph
//! inputs and outputs. Constitutive nonlinearity and memory remain in the law
//! nodes. The shared Newton--Krylov owner handles iteration; exact chain-rule
//! actions use the SAME committed history as every trial. Nothing publishes
//! history until the caller accepts the returned candidate. Balance equations
//! do not establish passivity, conservation or calibration of arbitrary laws.

/// Original material equations coupled to finite enthalpy and heat transport.
pub mod enthalpy;

/// The same material equations coupled to native spatial nodal enthalpy.
pub mod spatial;

use fs_exec::Cx;
use fs_material::graph::{ConstitutiveGraph, GraphError, GraphOutput};
use fs_qty::Dims;
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
    SolverRunProgress,
};
use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet};

/// A named input or output, in the graph's original node/port vocabulary.
pub type PortKey = (String, String);

/// A normalized solver coordinate `x` means `input = scale_si * x`.
#[derive(Debug, Clone)]
pub struct MaterialUnknown {
    /// Original unconnected graph input to solve.
    pub input: PortKey,
    /// Required SI base dimensions of that input.
    pub dims: Dims,
    /// Positive physical input value represented by one normalized coordinate.
    pub scale_si: f64,
}

/// An equation term can close an external feedback edge without putting a
/// cycle inside the graph. Coefficients are dimensionless; differing physical
/// dimensions need a real conversion/geometry law node, not a numeric trick.
#[derive(Debug, Clone)]
pub enum MaterialValue {
    /// Original external input, either fixed or solved.
    Input(PortKey),
    /// Native node output evaluated from the same committed history.
    Output(PortKey),
}

/// `sum(coefficient * value) = target_si`, normalized by `scale_si`.
#[derive(Debug, Clone)]
pub struct MaterialEquation {
    /// Dimensionless coefficients multiplying original physical values.
    pub terms: Vec<(f64, MaterialValue)>,
    /// Common SI dimensions of every term and the target.
    pub dims: Dims,
    /// Prescribed physical right-hand side in those SI dimensions.
    pub target_si: f64,
    /// Positive physical residual corresponding to one normalized residual.
    pub scale_si: f64,
}

/// Explicit size, evaluation and nonlinear work limits for one proposal.
#[derive(Debug, Clone, Copy)]
pub struct MaterialGraphSolve {
    /// Maximum number of independent external inputs and balance equations.
    pub maximum_unknowns: usize,
    /// Maximum total terms across all equations.
    pub maximum_terms: usize,
    /// Value and Jacobian calls, including the independent final admission.
    pub maximum_evaluations: usize,
    /// Maximum Newton iterations in the native nonlinear solver.
    pub maximum_iterations: usize,
    /// Every final equation must satisfy this dimensionless absolute tolerance,
    /// even when Newton's relative criterion alone would accept it.
    pub equation_tolerance: f64,
    /// Existing Newton, linear-solve and globalization controls on scaled values.
    pub newton: NewtonKrylovConfig,
}

/// Refusals preserve the original law or solver failure and accepted history.
#[derive(Debug)]
pub enum MaterialGraphError {
    /// Invalid dimensions, bindings, scales, controls or numerical values.
    Invalid(&'static str),
    /// Original graph/source/history refusal, without a replacement law.
    Graph(GraphError),
    /// Native nonlinear solver configuration or arithmetic refusal.
    Newton(NewtonError),
    /// Native solver exhausted its convergence path; no history was accepted.
    NotConverged(NewtonReport),
    /// Fresh physical equations failed the independent normalized endpoint gate.
    EquationResidual {
        /// Largest absolute normalized equation defect.
        maximum: f64,
        /// Caller-supplied absolute normalized acceptance limit.
        tolerance: f64,
    },
    /// Native solver stopped before completing its admitted iteration work.
    Interrupted(SolverRunProgress),
    /// A value/Jacobian/final-admission call exceeded the evaluation allowance.
    EvaluationBudget,
    /// Shared execution context cancelled before a candidate was returned.
    Cancelled,
}
impl core::fmt::Display for MaterialGraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constitutive material solve: {self:?}")
    }
}
impl std::error::Error for MaterialGraphError {}
impl From<GraphError> for MaterialGraphError {
    fn from(error: GraphError) -> Self {
        match error {
            GraphError::Cancelled => Self::Cancelled,
            other => Self::Graph(other),
        }
    }
}

/// A complete accepted numerical proposal, with original law state and outputs.
#[derive(Debug)]
pub struct MaterialGraphCandidate {
    /// Normalized solution, in the original unknown order.
    pub coordinates: Vec<f64>,
    /// Complete physical fixed and solved boundary inputs in SI units.
    pub external: BTreeMap<PortKey, f64>,
    /// Original node outputs and proposed aggregate history, accepted together.
    pub graph: GraphOutput,
    /// Largest fresh normalized equation residual; not a model-error bound.
    pub maximum_scaled_residual: f64,
    /// Report from the original Newton--Krylov owner.
    pub nonlinear: NewtonReport,
    /// Total value/Jacobian calls, including final physical admission.
    pub evaluations: usize,
}

/// Public residual/Jacobian adapter. All callbacks are pure with respect to
/// physical history; a source refusal latches and cannot become convergence.
/// Use `solve` to check the original failure before interpreting a solver
/// report, or check `take_failure` after directly using `NonlinearProblem`.
pub struct MaterialGraphResidual<'a, 'scope> {
    cx: &'a Cx<'scope>,
    graph: &'a ConstitutiveGraph,
    version: u64,
    state: &'a [f64],
    fixed: BTreeMap<PortKey, f64>,
    unknowns: Vec<MaterialUnknown>,
    equations: Vec<MaterialEquation>,
    controls: MaterialGraphSolve,
    evaluations: Cell<usize>,
    failed: Cell<bool>,
    failure: RefCell<Option<MaterialGraphError>>,
}

impl<'a, 'scope> MaterialGraphResidual<'a, 'scope> {
    /// Admit a square, dimensionally consistent feedback block over the original
    /// graph and committed state schema. An empty block admits prescribed-input
    /// material evolution: the native solver still evaluates the complete graph
    /// and the final source/history proposal, with zero nonlinear iterations.
    /// Fixed values and unknowns may name only external inputs; constitutive
    /// history remains borrowed and unchanged.
    pub fn new(
        cx: &'a Cx<'scope>,
        graph: &'a ConstitutiveGraph,
        version: u64,
        state: &'a [f64],
        fixed: BTreeMap<PortKey, f64>,
        unknowns: Vec<MaterialUnknown>,
        equations: Vec<MaterialEquation>,
        controls: MaterialGraphSolve,
    ) -> Result<Self, MaterialGraphError> {
        cx.checkpoint().map_err(|_| MaterialGraphError::Cancelled)?;
        controls
            .newton
            .validate()
            .map_err(MaterialGraphError::Newton)?;
        if unknowns.len() != equations.len()
            || unknowns.len() > controls.maximum_unknowns
            || controls.maximum_evaluations == 0
            || controls.maximum_iterations == 0
            || !positive(controls.equation_tolerance)
        {
            return Err(MaterialGraphError::Invalid(
                "square system and positive finite work/tolerance controls required",
            ));
        }
        graph.state_schema().decode(version, state)?;
        let mut keys = BTreeSet::new();
        for unknown in &unknowns {
            if !positive(unknown.scale_si)
                || !keys.insert(unknown.input.clone())
                || fixed.contains_key(&unknown.input)
                || graph.external_input_dims(&unknown.input.0, &unknown.input.1)
                    != Some(unknown.dims)
            {
                return Err(MaterialGraphError::Invalid(
                    "unknown must be a unique, dimensioned, unfixed external input with a positive scale",
                ));
            }
        }
        for ((node, port), value) in &fixed {
            if !value.is_finite() || graph.external_input_dims(node, port).is_none() {
                return Err(MaterialGraphError::Invalid(
                    "fixed values must name finite external graph inputs",
                ));
            }
        }
        let mut terms = 0usize;
        for equation in &equations {
            terms = terms
                .checked_add(equation.terms.len())
                .ok_or(MaterialGraphError::Invalid("equation term budget overflow"))?;
            if terms > controls.maximum_terms
                || equation.terms.is_empty()
                || !positive(equation.scale_si)
                || !equation.target_si.is_finite()
            {
                return Err(MaterialGraphError::Invalid(
                    "equation terms, scale, target or budget invalid",
                ));
            }
            for (coefficient, value) in &equation.terms {
                let dims = match value {
                    MaterialValue::Input((node, port)) => graph.external_input_dims(node, port),
                    MaterialValue::Output((node, port)) => graph.output_dims(node, port),
                };
                if !coefficient.is_finite() || *coefficient == 0.0 || dims != Some(equation.dims) {
                    return Err(MaterialGraphError::Invalid(
                        "equation term dimensions or coefficient invalid",
                    ));
                }
            }
        }
        Ok(Self {
            cx,
            graph,
            version,
            state,
            fixed,
            unknowns,
            equations,
            controls,
            evaluations: Cell::new(0),
            failed: Cell::new(false),
            failure: RefCell::new(None),
        })
    }

    /// Whether a solver callback has permanently refused this proposal.
    pub fn has_failed(&self) -> bool {
        self.failed.get()
    }
    /// Consumed value/Jacobian evaluations, including unsuccessful proposals.
    /// A containing physical solve can bound all nested work with this count.
    pub fn evaluations(&self) -> usize {
        self.evaluations.get()
    }
    /// Taking the first error does not clear the permanent failure latch.
    pub fn take_failure(&self) -> Option<MaterialGraphError> {
        self.failure.borrow_mut().take()
    }

    fn external(
        &self,
        x: &[f64],
        direction: bool,
    ) -> Result<BTreeMap<PortKey, f64>, MaterialGraphError> {
        if x.len() != self.unknowns.len() {
            return Err(MaterialGraphError::Invalid("coordinate dimension mismatch"));
        }
        let mut external = if direction {
            BTreeMap::new()
        } else {
            self.fixed.clone()
        };
        for (unknown, value) in self.unknowns.iter().zip(x) {
            let value = unknown.scale_si * value;
            if !value.is_finite() {
                return Err(MaterialGraphError::Invalid("nonfinite scaled input"));
            }
            external.insert(unknown.input.clone(), value);
        }
        Ok(external)
    }

    fn equations(
        &self,
        external: &BTreeMap<PortKey, f64>,
        outputs: &BTreeMap<PortKey, f64>,
        direction: bool,
    ) -> Result<Vec<f64>, MaterialGraphError> {
        self.equations
            .iter()
            .map(|equation| {
                let mut sum = if direction { 0.0 } else { -equation.target_si };
                for (coefficient, value) in &equation.terms {
                    let value = match value {
                        MaterialValue::Input(key) if direction => {
                            external.get(key).copied().unwrap_or(0.0)
                        }
                        MaterialValue::Input(key) => *external
                            .get(key)
                            .ok_or(MaterialGraphError::Invalid("unfed equation input"))?,
                        MaterialValue::Output(key) => *outputs
                            .get(key)
                            .ok_or(MaterialGraphError::Invalid("missing equation output"))?,
                    };
                    sum += coefficient * value;
                }
                let value = sum / equation.scale_si;
                if !value.is_finite() {
                    return Err(MaterialGraphError::Invalid(
                        "unrepresentable equation residual",
                    ));
                }
                Ok(value)
            })
            .collect()
    }

    fn begin(&self) -> Result<(), MaterialGraphError> {
        self.cx
            .checkpoint()
            .map_err(|_| MaterialGraphError::Cancelled)?;
        if self.failed.get() {
            return Err(MaterialGraphError::Invalid(
                "material solve previously refused",
            ));
        }
        if self.evaluations.get() >= self.controls.maximum_evaluations {
            return Err(MaterialGraphError::EvaluationBudget);
        }
        self.evaluations.set(self.evaluations.get() + 1);
        Ok(())
    }

    fn values(
        &self,
        x: &[f64],
    ) -> Result<(BTreeMap<PortKey, f64>, GraphOutput, Vec<f64>), MaterialGraphError> {
        self.begin()?;
        let external = self.external(x, false)?;
        let value = self
            .graph
            .execute_cancellable(self.version, self.state, &external, || {
                self.cx.checkpoint().is_err()
            })?;
        let residual = self.equations(&external, &value.outputs, false)?;
        Ok((external, value, residual))
    }

    fn publish(&self, result: Result<Vec<f64>, MaterialGraphError>, output: &mut [f64]) {
        output.fill(f64::NAN);
        let result = if output.len() == self.dimension() {
            result
        } else {
            Err(MaterialGraphError::Invalid(
                "solver output dimension mismatch",
            ))
        };
        match result {
            Ok(value) => output.copy_from_slice(&value),
            Err(error) => {
                if !self.failed.replace(true) {
                    *self.failure.borrow_mut() = Some(error);
                }
            }
        }
    }

    /// Solve and independently admit the final physical residual. Neither
    /// failure nor any number of Newton/Jacobian probes advances law memory.
    pub fn solve(
        &self,
        initial_coordinates: Vec<f64>,
    ) -> Result<MaterialGraphCandidate, MaterialGraphError> {
        let constructed = NewtonKrylovState::new(self, initial_coordinates, self.controls.newton);
        if let Some(error) = self.take_failure() {
            return Err(error);
        }
        let mut solver = constructed.map_err(MaterialGraphError::Newton)?;
        let run = solver.run_controlled(self, self.controls.maximum_iterations, &mut || {
            self.has_failed() || self.cx.checkpoint().is_err()
        });
        if let Some(error) = self.take_failure() {
            return Err(error);
        }
        self.cx
            .checkpoint()
            .map_err(|_| MaterialGraphError::Cancelled)?;
        if run.progress != SolverRunProgress::Complete {
            return Err(MaterialGraphError::Interrupted(run.progress));
        }
        if !run.report.converged {
            return Err(MaterialGraphError::NotConverged(run.report));
        }
        let (external, graph, residual) = self.values(&solver.x)?;
        let maximum = residual.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum > self.controls.equation_tolerance {
            return Err(MaterialGraphError::EquationResidual {
                maximum,
                tolerance: self.controls.equation_tolerance,
            });
        }
        self.cx
            .checkpoint()
            .map_err(|_| MaterialGraphError::Cancelled)?;
        Ok(MaterialGraphCandidate {
            coordinates: solver.x,
            external,
            graph,
            maximum_scaled_residual: maximum,
            nonlinear: run.report,
            evaluations: self.evaluations.get(),
        })
    }
}

impl NonlinearProblem for MaterialGraphResidual<'_, '_> {
    fn dimension(&self) -> usize {
        self.unknowns.len()
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        self.publish(self.values(x).map(|(_, _, residual)| residual), output);
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        let result = (|| {
            self.begin()?;
            let external = self.external(x, false)?;
            let direction = self.external(direction, true)?;
            let value =
                self.graph
                    .linearize(self.version, self.state, &external, &direction, || {
                        self.cx.checkpoint().is_err()
                    })?;
            self.equations(&direction, &value.output_direction, true)
        })();
        self.publish(result, output);
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
