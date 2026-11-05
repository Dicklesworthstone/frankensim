//! Entropy of the original stationary P1 enthalpy step.
//!
//! This is a discrete thermal account, not a continuum entropy enclosure. The
//! original caloric charts supply integral dh/T, including latent plateaus;
//! bonded films supply their original sensible caloric entropy. Bulk, surface and contact transport
//! use the SAME assembled operator as the accepted energy equation.
//!
//! A positive material tensor alone does not make every P1 mesh thermodynamically
//! compatible. Signed discrete production is retained, and the optional step
//! refuses negative production beyond the caller's absolute allowance. It never
//! clips coefficients, changes the mesh or labels numerical defects as heating.
//! See Limache and Idelsohn (2016), doi:10.1186/s40323-016-0066-8.
//!
//! Each active Robin triangle retains its actual positive affine bath field.
//! The original h(x), T(x) and B(x) determine its exact rational bath entropy;
//! no patch mean substitutes. Photon entropy requires its own integration.
//! Prescribed heat flux and volumetric watts remain unclassified inputs.

use std::fmt;

mod robin;

use fs_exec::Cx;

use super::{
    ASSEMBLY_TILE, AssembledSystem, DofMap, EnthalpyBackwardEuler, EnthalpyEndpointEvaluation,
    EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext, finite, poll,
};
use crate::{
    ConductionProblem, ThermalInterfaces,
    assemble::{assemble_source_load, robin_face_matrix},
    bc::ThermalBc,
};

/// Signed integrated entropy contributions in J/K over one physical interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnthalpyEntropyBalance {
    /// Actual original-chart entropy, plus the sensible films' entropy.
    pub thermal_storage_change_j_k: f64,
    /// Prescribed-temperature reactions at their actual time-quadrature temperatures.
    /// Includes fixed-node storage when a boundary temperature changes.
    pub dirichlet_entropy_input_j_k: f64,
    /// Integral of inward Robin heat divided by the actual affine bath B(x).
    pub robin_entropy_input_j_k: f64,
    /// Signed discrete bulk, surface and contact generation on this state.
    /// General P1 meshes can make this negative even with positive conductivity.
    pub conduction_generation_j_k: f64,
    /// Signed generation of the original consistent Robin face operator.
    /// No diagonal boundary approximation is substituted to enforce a sign.
    pub robin_generation_j_k: f64,
    /// Original integrated volumetric watts divided by receiving nodal T.
    /// This does not assign entropy to an unspecified heat/work source.
    pub source_at_node_j_k: f64,
    /// Signed prescribed heat input divided by receiving nodal T.
    /// An external reservoir temperature is not inferred from a heat flux.
    pub neumann_at_node_j_k: f64,
    /// True if any volumetric source or prescribed flux is nonzero.
    pub has_unclassified_inputs: bool,
    /// Exact state entropy minus native time-weighted energy increments
    /// divided by their actual stage temperatures.
    /// This is time discretization, not physical conductive generation.
    pub temporal_defect_j_k: f64,
    /// Time-weighted free-row stage residuals divided by actual stage temperature.
    pub algebraic_residual_j_k: f64,
    /// Remaining operator/reduction roundoff after ALL accounts above.
    /// No physical source or balancing correction is created from this value.
    pub balance_residual_j_k: f64,
}

/// Refusal of physical evaluation or the explicitly requested entropy gate.
#[derive(Clone, Debug, PartialEq)]
pub enum EnthalpyEntropyError {
    /// Original constitutive, source, work, energy or cancellation refusal.
    Physical(EnthalpyError),
    /// A measured discrete entropy account is below the caller's allowance.
    NegativeProduction {
        account: &'static str,
        value_j_k: f64,
        tolerance_j_k: f64,
    },
    /// The independent entropy identity does not close within the allowance.
    Balance {
        residual_j_k: f64,
        tolerance_j_k: f64,
    },
}

impl From<EnthalpyError> for EnthalpyEntropyError {
    fn from(error: EnthalpyError) -> Self {
        Self::Physical(error)
    }
}
impl fmt::Display for EnthalpyEntropyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Physical(error) => error.fmt(f),
            Self::NegativeProduction {
                account,
                value_j_k,
                tolerance_j_k,
            } => write!(
                f,
                "{account} entropy {value_j_k} J/K is below -{tolerance_j_k} J/K"
            ),
            Self::Balance {
                residual_j_k,
                tolerance_j_k,
            } => write!(
                f,
                "entropy residual {residual_j_k} J/K exceeds {tolerance_j_k} J/K"
            ),
        }
    }
}
impl std::error::Error for EnthalpyEntropyError {}

impl EnthalpyEntropyBalance {
    /// Check this state's transport signs, thermal storage-minus-input
    /// inequality and independent closure. This does not prove nodal
    /// monotonicity for other states or classify an external heat/work owner.
    pub fn admit(&self, tolerance_j_k: f64) -> Result<(), EnthalpyEntropyError> {
        admit_tolerance(tolerance_j_k)?;
        for value in [
            self.thermal_storage_change_j_k,
            self.dirichlet_entropy_input_j_k,
            self.robin_entropy_input_j_k,
            self.conduction_generation_j_k,
            self.robin_generation_j_k,
            self.source_at_node_j_k,
            self.neumann_at_node_j_k,
            self.temporal_defect_j_k,
            self.algebraic_residual_j_k,
            self.balance_residual_j_k,
        ] {
            finite(value)?;
        }
        let available = finite(
            self.thermal_storage_change_j_k
                - self.dirichlet_entropy_input_j_k
                - self.robin_entropy_input_j_k
                - self.source_at_node_j_k
                - self.neumann_at_node_j_k,
        )?;
        for (account, value) in [
            ("bulk/surface/contact", self.conduction_generation_j_k),
            ("Robin", self.robin_generation_j_k),
            ("thermal storage minus inputs", available),
        ] {
            finite(value)?;
            if value < -tolerance_j_k {
                return Err(EnthalpyEntropyError::NegativeProduction {
                    account,
                    value_j_k: value,
                    tolerance_j_k,
                });
            }
        }
        if finite(self.balance_residual_j_k)?.abs() > tolerance_j_k {
            return Err(EnthalpyEntropyError::Balance {
                residual_j_k: self.balance_residual_j_k,
                tolerance_j_k,
            });
        }
        Ok(())
    }
}

fn admit_tolerance(tolerance: f64) -> Result<(), EnthalpyError> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(EnthalpyError::InvalidInput(
            "entropy tolerance must be finite and nonnegative",
        ));
    }
    Ok(())
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance the native backward-Euler solve and admit its discrete thermal
    /// entropy before returning either state or account. All existing Newton,
    /// energy, material, boundary and cancellation gates remain in force.
    ///
    /// Entropy is evaluated once on the final assembled operator, not during
    /// Newton iterations. Additional work is O(vertices + nonzeros + faces),
    /// plus original caloric segment queries, with one nodal source buffer.
    /// Unclassified inputs remain explicit; acceptance is not a whole-system
    /// second-law certificate or a spatial/time accuracy guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_entropy(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        entropy_tolerance_j_k: f64,
    ) -> Result<(EnthalpyStepSolution, EnthalpyEntropyBalance), EnthalpyEntropyError> {
        admit_tolerance(entropy_tolerance_j_k)?;
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        admit_boundary(cx, problem)?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation: None,
        };
        let (h, report) = context.solve_newton(config, None)?;
        let (solution, entropy) =
            context.finish_and_entropy(h, report, config.energy_tolerance_j, true)?;
        let entropy = entropy.expect("requested native entropy account");
        entropy.admit(entropy_tolerance_j_k)?;
        poll(cx, solution.temperature.len())?;
        Ok((solution, entropy))
    }

    /// Measure an arbitrary nonradiative endpoint with the original physical
    /// operator. This performs no solve and applies no acceptance gate: signed
    /// production and the actual local residual remain available to a coupled
    /// owner, including a candidate that must subsequently be refused.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_entropy(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<(EnthalpyEndpointEvaluation, EnthalpyEntropyBalance), EnthalpyError> {
        let (endpoint, entropy) = self.evaluate_endpoint_and_entropy(
            cx,
            problem,
            interfaces,
            old_specific_h,
            candidate_specific_h,
            dt_s,
            config,
            None,
            None,
            true,
        )?;
        Ok((
            endpoint,
            entropy.expect("requested endpoint entropy account"),
        ))
    }
}

// Every source field has already passed native shape admission. The exact
// boundary partition must also match before indexing its face assignments.
pub(super) fn admit_boundary(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
) -> Result<(), EnthalpyError> {
    if problem.boundary.face_count() != problem.mesh.boundary().len() {
        return Err(EnthalpyError::InvalidInput(
            "entropy boundary partition must match the original mesh",
        ));
    }
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % ASSEMBLY_TILE == 0 {
            poll(cx, slot)?;
        }
        if let Some(ThermalBc::Robin { htc, t_ref }) = problem.boundary.condition_for(slot) {
            let vertices = face.vertices.map(|v| v as usize);
            if vertices.iter().all(|&v| htc.at(v) == 0.0) {
                continue;
            }
            let temperatures = vertices.map(|v| t_ref.at(v));
            if temperatures.iter().any(|t| !t.is_finite() || *t <= 0.0) {
                return Err(EnthalpyError::InvalidInput(
                    "entropy requires positive actual Robin bath temperatures",
                ));
            }
        }
    }
    Ok(())
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), EnthalpyError> {
        let term = finite(term)?;
        let next = finite(self.value + term)?;
        let lost = if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, EnthalpyError> {
        finite(self.value + self.correction)
    }
}

// -Kij*(Ti-Tj)^2/(Ti*Tj), without squaring a dimensional temperature.
fn pair_generation(entry: f64, ti: f64, tj: f64) -> Result<f64, EnthalpyError> {
    if entry == 0.0 || ti == tj {
        return Ok(0.0);
    }
    let delta = ti - tj;
    finite((-entry * (delta / ti)) * (delta / tj))
}

pub(super) fn evaluate(
    context: &StepContext<'_, '_, '_, '_>,
    h: &[f64],
    system: &AssembledSystem,
    temperature: &[f64],
    dofs: &DofMap,
) -> Result<EnthalpyEntropyBalance, EnthalpyError> {
    if context.pointwise_radiation.is_some() {
        return Err(EnthalpyError::InvalidInput(
            "conductive entropy does not assign photon entropy to radiation",
        ));
    }
    admit_boundary(context.cx, context.problem)?;
    let model = context.storage;
    let source_load = assemble_source_load(context.cx, model.mesh, context.problem.source)?;
    let (mut state, mut temporal, mut residual, mut fixed, mut source) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    let (mut conduction, mut robin, mut bath, mut neumann) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    let mut has_unclassified = false;
    for (vertex, &ti) in temperature.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(context.cx, vertex)?;
        }
        if !ti.is_finite() || ti <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "entropy requires positive absolute temperatures",
            ));
        }
        let curve = model.curve_for_vertex(vertex);
        let mut ds = finite(
            model.masses[vertex]
                * curve
                    .specific_entropy_change_j_kg_k(context.old[vertex], h[vertex])
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?,
        )?;
        ds = finite(ds + model.film_entropy_change_j_k(vertex, context.old[vertex], ti)?)?;
        let de = model.storage_change_j(vertex, h[vertex], context.old[vertex], ti)?;
        state.add(ds)?;
        temporal.add(finite(ds - de / ti)?)?;
        source.add(finite(source_load[vertex] / ti)?)?;
        has_unclassified |= context.problem.source.at(vertex) != 0.0;
        let (columns, entries) = system.operator.row(vertex);
        let mut transport = 0.0;
        for (index, (&column, &entry)) in columns.iter().zip(entries).enumerate() {
            if index != 0 && index % ASSEMBLY_TILE == 0 {
                poll(context.cx, index)?;
            }
            transport = finite(transport + finite(entry * temperature[column])?)?;
            if column != vertex {
                conduction.add(pair_generation(0.5 * entry, ti, temperature[column])?)?;
            }
        }
        let row = finite(de + context.dt * finite(transport - system.load[vertex])?)?;
        if dofs.slot_of(vertex).is_some() {
            residual.add(finite(row / ti)?)?;
        } else {
            fixed.add(finite(row / ti)?)?;
        }
    }

    for (slot, face) in model.mesh.boundary().iter().enumerate() {
        if slot % ASSEMBLY_TILE == 0 {
            poll(context.cx, slot)?;
        }
        let vertices = face.vertices.map(|v| v as usize);
        match context.problem.boundary.condition_for(slot) {
            Some(ThermalBc::Neumann { outward_flux }) => {
                for (a, &va) in vertices.iter().enumerate() {
                    let mut out = 0.0_f64;
                    for (b, &vb) in vertices.iter().enumerate() {
                        let weight = if a == b {
                            face.area / 6.0
                        } else {
                            face.area / 12.0
                        };
                        out = weight.mul_add(outward_flux.at(vb), out);
                        has_unclassified |= outward_flux.at(vb) != 0.0;
                    }
                    neumann.add(finite(-out / temperature[va])?)?;
                }
            }
            Some(ThermalBc::Robin { htc, t_ref }) => {
                if vertices.iter().all(|&v| htc.at(v) == 0.0) {
                    continue;
                }
                let baths = vertices.map(|v| t_ref.at(v));
                let bath_reference = baths.into_iter().fold(0.0_f64, f64::max);
                let transfer = vertices.map(|v| htc.at(v));
                let temperatures = vertices.map(|v| temperature[v]);
                let block = robin_face_matrix(face.area, transfer);
                let delta = std::array::from_fn::<_, 3, _>(|a| temperatures[a] - baths[a]);
                let correction = robin::correction(
                    context.cx, face.area, transfer, temperatures, baths,
                )?;
                for (a, &va) in vertices.iter().enumerate() {
                    let mut out = 0.0_f64;
                    for b in 0..3 {
                        out = block[a][b].mul_add(delta[b], out);
                    }
                    bath.add(finite(-out / bath_reference)?)?;
                    // The reference cancels algebraically with the exact
                    // rational correction below; the original B(x) is retained.
                    robin.add(finite((out / temperature[va])
                        * ((temperature[va] - bath_reference) / bath_reference))?)?;
                    for b in (a + 1)..3 {
                        // Remove the Robin off-diagonal from the full sparse
                        // operator pair sum, leaving original bulk/contact K.
                        conduction.add(-pair_generation(
                            block[a][b],
                            temperature[va],
                            temperature[vertices[b]],
                        )?)?;
                    }
                }
                bath.add(-correction)?;
                robin.add(correction)?;
            }
            Some(ThermalBc::Dirichlet { .. }) | None => {}
        }
    }
    let mut result = EnthalpyEntropyBalance {
        thermal_storage_change_j_k: state.total()?,
        dirichlet_entropy_input_j_k: fixed.total()?,
        robin_entropy_input_j_k: finite(context.dt * bath.total()?)?,
        conduction_generation_j_k: finite(context.dt * conduction.total()?)?,
        robin_generation_j_k: finite(context.dt * robin.total()?)?,
        source_at_node_j_k: finite(context.dt * source.total()?)?,
        neumann_at_node_j_k: finite(context.dt * neumann.total()?)?,
        has_unclassified_inputs: has_unclassified,
        temporal_defect_j_k: temporal.total()?,
        algebraic_residual_j_k: residual.total()?,
        balance_residual_j_k: 0.0,
    };
    let mut closure = Sum::default();
    for term in [
        result.thermal_storage_change_j_k,
        -result.dirichlet_entropy_input_j_k,
        -result.robin_entropy_input_j_k,
        -result.conduction_generation_j_k,
        -result.robin_generation_j_k,
        -result.source_at_node_j_k,
        -result.neumann_at_node_j_k,
        -result.temporal_defect_j_k,
        -result.algebraic_residual_j_k,
    ] {
        closure.add(term)?;
    }
    result.balance_residual_j_k = closure.total()?;
    poll(context.cx, h.len())?;
    Ok(result)
}
