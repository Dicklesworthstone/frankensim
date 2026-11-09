//! Original-law nonlinear conduction on the native full-tensor face topology.
//!
//! Traces remain physical Kelvin temperatures. Cell unknowns are temperature
//! or the original caloric energy coordinate. Conductivity laws retain their
//! exact Kelvin interval integrals or explicitly select cellwise K(e). The
//! original endpoint tensors, dT/de and optional analytic dK/de give the Newton
//! derivative. Cells are statically condensed, leaving only face increments
//! in the native sparse solve.

use fs_exec::CancelGate;
use fs_sparse::Coo;

use super::{
    ChtError, ConductivityTensor, EnergyConfig, Prepared, PseudoStep, SolidMaterial,
    TensorConductionSolution, ThermalFace, ThermalSetup, VoxelDomain, audit_endpoint_with_contacts,
    bicgstab_ilu0, finite, finite_positive, invalid, poll, prepare_geometry, residual_norms,
    scale_rows, validate_sensible_storage,
};

#[cfg(feature = "enthalpy-coupling")]
mod caloric;
mod coordinate;
mod local;
mod material;
mod sdirk;
mod transient;
#[cfg(feature = "enthalpy-coupling")]
pub(in crate::conjugate) use caloric::{
    CaloricAdaptiveForcing, CaloricAdaptiveTrial, advance_caloric_adaptive_with, embedded_ratio,
    advance_radiating_caloric_prepared, advance_radiating_sdirk2_prepared,
};
#[cfg(feature = "enthalpy-coupling")]
pub use caloric::{
    TensorCaloricAdaptiveConfig, TensorCaloricAdaptiveError, TensorCaloricAdaptiveReport,
    TensorCaloricAdaptiveSolution, TensorCaloricForcing, TensorCaloricForcingSample,
    TensorCaloricEndpoint, TensorCaloricLaw, TensorCaloricSdirk2Step, TensorCaloricStep,
    advance_caloric_tensor_conduction_adaptive,
    TensorCaloricStorage, advance_caloric_tensor_conduction,
    advance_caloric_tensor_conduction_sdirk2,
    TensorCaloricRadiationReport, advance_caloric_tensor_conduction_radiating,
    advance_caloric_tensor_conduction_radiating_sdirk2,
    TensorCaloricRadiatingAdaptiveSolution, advance_caloric_tensor_conduction_radiating_adaptive,
};
use coordinate::{CellCoordinate, KelvinCoordinate};
use super::contact::{ContactHeat, PreparedSources, joint_heat};
use local::{KirchhoffCell, KirchhoffRecovery};
pub use material::{RotatedSolidConductivity, RotatedSourceConductivity, TensorConductivityLaw};
#[cfg(feature = "enthalpy-coupling")]
pub use material::{
    IsotropicPhaseConductivity, TensorCaloricConductivity, TensorCaloricConductivityPoint,
};
pub use sdirk::{TensorNonlinearSdirk2Step, advance_nonlinear_tensor_conduction_sdirk2};
pub use transient::{TensorNonlinearStep, advance_nonlinear_tensor_conduction};

/// Bounded nonlinear and algebraic work for one original-law spatial endpoint,
/// or shared over both stages of a physical SDIRK2 interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorNonlinearConfig {
    /// Native linear tolerance and total Krylov budget over ALL Newton solves.
    /// The same tolerance admits every original physical heat balance.
    /// The convection scheme is unused by this conduction-only entry.
    pub linear: EnergyConfig,
    /// Maximum Newton corrections, including zero-iteration linear solves.
    pub max_newton_iterations: usize,
    /// Total rejected line-search candidates over the complete endpoint.
    /// This is not reset for each Newton correction.
    pub max_backtracks: usize,
}

impl Default for TensorNonlinearConfig {
    fn default() -> Self {
        Self {
            linear: EnergyConfig::default(),
            max_newton_iterations: 32,
            max_backtracks: 128,
        }
    }
}

/// One accepted original-law endpoint; no trial material history is published.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorNonlinearSolution {
    /// Actual cell/face temperatures, integrated face heat, native physical
    /// ledger and original local heat acceptance. Its linear residual belongs
    /// to the last Newton correction, not a frozen-temperature surrogate.
    pub solution: TensorConductionSolution,
    /// Number of Newton corrections used.
    pub newton_iterations: usize,
    /// Total rejected original-law candidates.
    pub backtracks: usize,
    /// Largest absolute original free-cell or free-trace equation defect, W,
    /// before the native arithmetic roundoff allowance.
    pub max_equation_defect_w: f64,
}

/// Solve steady full-tensor conduction with one immutable original law per cell.
///
/// The laws field is the explicit conductivity owner. Material indices still
/// define occupancy and native material-pair contacts; bind native temperature
/// tables through RotatedSolidConductivity to retain their original source and
/// declared principal directions. A constant ConductivityTensor also
/// implements the law. No law is inferred from a cell label or an average T.
///
/// Each cell-to-face interval uses the law's exact tensor integral. For
/// K(T)=k(T)K0 this is the scalar Kirchhoff transform of native full-tensor HMM;
/// for grid-diagonal K(T) it is the integrated directional half-cell law.
/// Arbitrary full K(T) uses the documented nonsymmetric columnwise primitive
/// extension. Its smooth local flux consistency and exact derivative do not
/// establish a general discrete entropy/maximum principle or unconditional
/// nonlinear convergence. Unsupported or unconverged endpoints are refused.
///
/// The caller supplies one positive cell-temperature seed. Prescribed cells
/// take their authored temperature; perfect interfaces start at the two-cell
/// mean, fixed-resistance contact traces at their respective cell temperatures,
/// and prescribed boundary traces at the authored Kelvin value. Source contacts
/// start from their original declared surface-temperature query, so their bulk
/// centres need not lie in the interface source box. Every initial cell-to-face
/// interval must still be supported by the original conductivity law.
///
/// Newton corrections use original endpoint tensors, including all mixed
/// conductivity components, without an extra K'(T) term. Local cell increments
/// are eliminated before assembling the face-only system. A bounded line
/// search reevaluates every original interval and requires decrease of the
/// maximum absolute equation defect. Final acceptance independently checks
/// local cell, interface, contact and boundary heat using the shared native
/// physical ledger. Both source support and positive Kelvin values hold at
/// the returned endpoint. All inputs remain unchanged on any outcome.
///
/// SourceFaceContact binds its original positive R(T_surface_A,T_surface_B)
/// directly to the two physical traces. Every trial retains the selected source,
/// fixed external conditions and A/B surface order. Its full analytic resistance
/// chain rule enters the SAME sparse Newton system; no outer frozen-R iteration
/// or extra contact heat is introduced. Curve/source-cell corners use original
/// one-sided primal tangents and confer no unique-adjoint claim. Physical
/// source receipts and surface heat are returned with the solution.
/// TensorNonlinearConfig owns all iteration limits; scalar SourceContactControls
/// do not govern this direct Newton path. Other geometry, fixed/film/flux/sink
/// scope matches frozen tensor conduction. Advection, eddy additions and compact
/// components retain their existing owners.
///
/// # Errors
/// Invalid ownership, unsupported source intervals, nonphysical/unrepresentable
/// fields, an unanchored steady problem, exhausted nonlinear/linear/backtracking
/// budgets, a failed sparse solve, or cancellation returns no partial solution.
#[allow(clippy::too_many_arguments)]
pub fn solve_nonlinear_tensor_conduction(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    initial_temperature_k: &[f64],
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorNonlinearSolution, ChtError> {
    solve_inner(
        domain,
        solids,
        setup,
        laws,
        initial_temperature_k,
        config,
        None,
        None,
        gate,
    )
}

#[derive(Clone, Copy)]
struct CoordinateValue {
    value: f64,
    temperature_slope: f64,
}

struct Endpoint {
    solution: TensorConductionSolution,
    // None for Kelvin, avoiding duplicate cell storage in all existing paths.
    coordinates: Option<Vec<CoordinateValue>>,
    trace: Vec<f64>,
    // Source-order heat and full two-trace derivatives from this SAME trial.
    source_contact_heat: Vec<ContactHeat>,
    cell_residual: Vec<f64>,
    trace_residual: Vec<f64>,
    merit: f64,
}

impl Endpoint {
    fn value(&self, cell: usize) -> f64 {
        self.coordinates
            .as_ref()
            .map_or(self.solution.temperature[cell], |values| values[cell].value)
    }

    fn temperature_slope(&self, cell: usize) -> f64 {
        self.coordinates
            .as_ref()
            .map_or(1.0, |values| values[cell].temperature_slope)
    }
}

struct Direction {
    cell: Vec<f64>,
    trace: Vec<f64>,
    iterations: usize,
    nonzeros: usize,
    relative_residual: f64,
}

/// The temperature primitive remains the default. A caloric constitutive law
/// explicitly selects cellwise K(e) without assigning an enthalpy to a face.
#[derive(Clone, Copy)]
enum CellLaw {
    Kirchhoff,
    #[cfg(feature = "enthalpy-coupling")]
    Caloric(material::TensorCaloricConductivity),
}

fn cell_law<C: CellCoordinate>(
    coordinate: &C,
    law: &dyn TensorConductivityLaw,
    cell: usize,
    value: f64,
    temperature: f64,
) -> Result<CellLaw, ChtError> {
    #[cfg(feature = "enthalpy-coupling")]
    if let Some(point) = coordinate.conductivity_point(cell, value, temperature)?
        && let Some(response) = law.caloric_conductivity(point)?
    {
        return Ok(CellLaw::Caloric(response));
    }
    #[cfg(not(feature = "enthalpy-coupling"))]
    let _ = (coordinate, law, cell, value, temperature);
    Ok(CellLaw::Kirchhoff)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_inner(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    initial_temperature_k: &[f64],
    config: &TensorNonlinearConfig,
    step: Option<&PseudoStep<'_>>,
    carried_power_w: Option<&[f64]>,
    gate: &CancelGate,
) -> Result<TensorNonlinearSolution, ChtError> {
    poll(gate)?;
    validate_sensible_storage(domain, step, gate)?;
    solve_coordinate(
        domain,
        solids,
        setup,
        laws,
        initial_temperature_k,
        config,
        &KelvinCoordinate { step },
        carried_power_w,
        gate,
    )
    .map(|(solution, _)| solution)
}

/// Shared original-law solve: Kelvin and caloric coordinates use exactly the
/// same face topology, heat residuals, sparse solve and work budgets.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_coordinate<C: CellCoordinate>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    initial_values: &[f64],
    config: &TensorNonlinearConfig,
    coordinate: &C,
    carried_power_w: Option<&[f64]>,
    gate: &CancelGate,
) -> Result<(TensorNonlinearSolution, Option<Vec<CoordinateValue>>), ChtError> {
    poll(gate)?;
    let cells = domain.cell_count();
    if laws.len() != cells || initial_values.len() != cells {
        return Err(invalid(
            "nonlinear tensor conduction needs one original law and seed per cell",
        ));
    }
    if let Some(carry) = carried_power_w {
        if coordinate.reference_temperatures().is_none() || carry.len() != cells {
            return Err(invalid(
                "stage heat carry needs original storage and one value per cell",
            ));
        }
        for (cell, &value) in carry.iter().enumerate() {
            if cell.is_multiple_of(4096) {
                poll(gate)?;
            }
            finite("thermal.tensor_conduction.carried_power_w", value)?;
        }
    }
    let mut values = initial_values.to_vec();
    for (cell, value) in values.iter_mut().enumerate() {
        if cell.is_multiple_of(4096) {
            poll(gate)?;
        }
        *value = coordinate.normalize(cell, *value)?;
    }
    // Bounds/duplicates are still admitted by the shared geometry owner.
    // This pre-assignment merely selects the correct original-law query.
    for &(cell, t) in &setup.fixed_temperature {
        if cell >= cells {
            return Err(invalid("fixed-temperature cell is outside the domain"));
        }
        finite_positive("thermal.fixed_temperature", t)?;
        values[cell] = coordinate.prescribed_value(cell, t)?;
    }
    let temperature: std::borrow::Cow<'_, [f64]> = if C::IS_TEMPERATURE {
        std::borrow::Cow::Borrowed(&values)
    } else {
        std::borrow::Cow::Owned(
            values
                .iter()
                .enumerate()
                .map(|(cell, &value)| {
                    if cell.is_multiple_of(256) {
                        poll(gate)?;
                    }
                    coordinate
                        .temperature_and_slope(cell, value)
                        .map(|state| state[0])
                })
                .collect::<Result<Vec<_>, ChtError>>()?,
        )
    };
    let mut point = Vec::<ConductivityTensor>::with_capacity(cells);
    let mut geometry = Vec::with_capacity(cells);
    for cell in 0..cells {
        if cell.is_multiple_of(512) {
            poll(gate)?;
        }
        finite_positive("thermal.tensor_conduction.initial_temperature_k", temperature[cell])?;
        point.push(match cell_law(
            coordinate, laws[cell], cell, values[cell], temperature[cell],
        )? {
            CellLaw::Kirchhoff => laws[cell].conductivity_at_temperature(temperature[cell])?,
            #[cfg(feature = "enthalpy-coupling")]
            CellLaw::Caloric(response) => response.conductivity(),
        });
        geometry.push(
            KirchhoffCell::new(domain.widths(cell))
                .ok_or_else(|| invalid("nonlinear tensor cell geometry is unrepresentable"))?,
        );
    }
    let sources = PreparedSources::new(domain, setup, gate)?;
    let data = prepare_geometry(
        domain,
        solids,
        setup,
        &point,
        &config.linear,
        coordinate.reference_temperatures(),
        Some(&sources),
        gate,
    )?;
    drop(point);
    let mut trace = initial_trace(domain, &data, &temperature, gate)?;
    sources.seed_traces(&data, &mut trace, gate)?;
    drop(temperature);
    let mut endpoint = evaluate(
        domain,
        setup,
        laws,
        &geometry,
        &data,
        &sources,
        coordinate,
        carried_power_w,
        values,
        trace,
        gate,
    )?;
    let (mut corrections, mut backtracks, mut iterations) = (0, 0, 0);
    let (mut nonzeros, mut last_linear_residual) = (0, 0.0);
    loop {
        poll(gate)?;
        if endpoint.solution.report.max_heat_residual <= config.linear.tolerance {
            endpoint.solution.report.unknowns = data.known.len();
            endpoint.solution.report.nonzeros = nonzeros;
            endpoint.solution.report.iterations = iterations;
            endpoint.solution.report.relative_residual = last_linear_residual;
            return Ok((TensorNonlinearSolution {
                solution: endpoint.solution,
                newton_iterations: corrections,
                backtracks,
                max_equation_defect_w: endpoint.merit,
            }, endpoint.coordinates));
        }
        if corrections == config.max_newton_iterations {
            return Err(not_converged(
                "nonlinear tensor endpoint",
                iterations,
                endpoint.solution.report.max_heat_residual,
                config.linear.tolerance,
            ));
        }
        let remaining = config
            .linear
            .max_iterations
            .checked_sub(iterations)
            .ok_or_else(|| invalid("nonlinear tensor solve exceeded its shared Krylov budget"))?;
        let direction = linearize(
            domain,
            setup,
            laws,
            &geometry,
            &data,
            coordinate,
            &endpoint,
            &config.linear,
            remaining,
            gate,
        )?;
        iterations += direction.iterations;
        corrections += 1;
        nonzeros = direction.nonzeros;
        last_linear_residual = direction.relative_residual;
        let mut fraction = 1.0_f64;
        loop {
            poll(gate)?;
            let values = (0..cells)
                .map(|cell| {
                    if data.fixed[cell].is_some() {
                        endpoint.value(cell)
                    } else {
                        fraction.mul_add(direction.cell[cell], endpoint.value(cell))
                    }
                })
                .collect();
            let trace = endpoint
                .trace
                .iter()
                .enumerate()
                .map(|(id, &t)| {
                    data.known[id].unwrap_or_else(|| fraction.mul_add(direction.trace[id], t))
                })
                .collect();
            match evaluate(
                domain,
                setup,
                laws,
                &geometry,
                &data,
                &sources,
                coordinate,
                carried_power_w,
                values,
                trace,
                gate,
            ) {
                Ok(candidate)
                    if candidate.solution.report.max_heat_residual <= config.linear.tolerance
                        || candidate.merit < endpoint.merit * (1.0 - 1.0e-4 * fraction) =>
                {
                    endpoint = candidate;
                    break;
                }
                // Source work and constitutive refusals cannot be repaired by
                // trying another state, and must retain their original cause.
                Err(error @ ChtError::InvalidInput {
                    field: "thermal.tensor_caloric.original_source_refusal"
                        | "gas.conduction.original_source_refusal",
                    ..
                }) => return Err(error),
                // A trial may leave source support or positive Kelvin. Only
                // that candidate is rejected; original history is immutable.
                Ok(_) | Err(ChtError::InvalidInput { .. }) => {}
                Err(error) => return Err(error),
            }
            if backtracks == config.max_backtracks {
                return Err(not_converged(
                    "nonlinear tensor line search",
                    iterations,
                    endpoint.solution.report.max_heat_residual,
                    config.linear.tolerance,
                ));
            }
            backtracks += 1;
            fraction *= 0.5;
            if fraction == 0.0 {
                return Err(not_converged(
                    "nonlinear tensor line search",
                    iterations,
                    endpoint.solution.report.max_heat_residual,
                    config.linear.tolerance,
                ));
            }
        }
    }
}

fn not_converged(
    system: &'static str,
    iterations: usize,
    relative_residual: f64,
    tolerance: f64,
) -> ChtError {
    ChtError::SolverNotConverged {
        system,
        iterations,
        relative_residual,
        tolerance,
    }
}

fn initial_trace(
    domain: &VoxelDomain,
    data: &Prepared,
    temperature: &[f64],
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    let mut trace = vec![0.0; data.known.len()];
    for c in 0..domain.cell_count() {
        if c.is_multiple_of(512) {
            poll(gate)?;
        }
        for f in 0..6 {
            let id = data.face_ids[c][f];
            if trace[id] > 0.0 {
                continue;
            }
            trace[id] = if let Some(t) = data.known[id] {
                t
            } else if let Some(n) = domain.neighbor(c, f)
                && data.face_ids[n][f ^ 1] == id
            {
                0.5 * temperature[c] + 0.5 * temperature[n]
            } else {
                temperature[c]
            };
            finite_positive(
                "thermal.tensor_conduction.initial_face_temperature_k",
                trace[id],
            )?;
        }
    }
    Ok(trace)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn evaluate<C: CellCoordinate>(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    geometry: &[KirchhoffCell],
    data: &Prepared,
    sources: &PreparedSources<'_>,
    coordinate: &C,
    carried_power_w: Option<&[f64]>,
    values: Vec<f64>,
    trace: Vec<f64>,
    gate: &CancelGate,
) -> Result<Endpoint, ChtError> {
    let cells = domain.cell_count();
    let (temperature, coordinates) = if C::IS_TEMPERATURE {
        (values, None)
    } else {
        let mut temperature = Vec::with_capacity(cells);
        let mut coordinates = Vec::with_capacity(cells);
        for (cell, value) in values.into_iter().enumerate() {
            if cell.is_multiple_of(512) {
                poll(gate)?;
            }
            let value = coordinate.normalize(cell, value)?;
            let [t, slope] = coordinate.temperature_and_slope(cell, value)?;
            finite_positive("thermal.tensor_conduction.temperature_k", t)?;
            finite("thermal.tensor_conduction.temperature_slope", slope)?;
            if slope < 0.0 {
                return Err(invalid(
                    "the original caloric temperature slope must be nonnegative",
                ));
            }
            temperature.push(t);
            coordinates.push(CoordinateValue {
                value,
                temperature_slope: slope,
            });
        }
        (temperature, Some(coordinates))
    };
    for (id, &t) in trace.iter().enumerate() {
        if id.is_multiple_of(4096) {
            poll(gate)?;
        }
        finite_positive("thermal.tensor_conduction.face_temperature_k", t)?;
    }
    let mut faces = Vec::with_capacity(cells);
    let mut heat = Vec::with_capacity(cells);
    let mut arithmetic = Vec::with_capacity(cells);
    let mut cell_residual = vec![0.0; cells];
    let mut trace_residual = vec![0.0; trace.len()];
    for c in 0..cells {
        if c.is_multiple_of(256) {
            poll(gate)?;
        }
        let t = temperature[c];
        finite_positive("thermal.tensor_conduction.temperature_k", t)?;
        let face = data.face_ids[c].map(|id| trace[id]);
        let value = coordinates.as_ref().map_or(t, |values| values[c].value);
        let means = match cell_law(coordinate, laws[c], c, value, t)? {
            CellLaw::Kirchhoff => {
                let mut means = [[[0.0; 3]; 3]; 6];
                for f in 0..6 {
                    means[f] = laws[c].mean_conductivity(t, face[f])?.as_array();
                }
                means
            }
            #[cfg(feature = "enthalpy-coupling")]
            CellLaw::Caloric(response) => [response.conductivity().as_array(); 6],
        };
        let actual = geometry[c]
            .heat(t, face, means)
            .ok_or_else(|| invalid("original tensor integral heat is unrepresentable"))?;
        if data.fixed[c].is_none() {
            let sink = data.sinks[c];
            cell_residual[c] = actual.outward_w.iter().sum::<f64>()
                + sink.conductance * (t - data.reference)
                - sink.relative_load
                - setup.power_w.get(c).copied().unwrap_or(0.0)
                - carried_power_w.map_or(0.0, |carry| carry[c]);
            cell_residual[c] += coordinate.storage_rate(c, value)?.power_w;
        }
        for (f, &q) in actual.outward_w.iter().enumerate() {
            let id = data.face_ids[c][f];
            if data.known[id].is_some() {
                continue;
            }
            trace_residual[id] -= q;
            if domain.neighbor(c, f).is_none() {
                let area = domain.face_area(c, f / 2);
                match setup.faces[f] {
                    ThermalFace::HeatFlux(into) => trace_residual[id] -= area * into,
                    ThermalFace::Convective { h, ambient } => {
                        trace_residual[id] += area * h * (face[f] - ambient)
                    }
                    _ => {}
                }
            }
        }
        faces.push(face);
        heat.push(actual.outward_w);
        arithmetic.push(actual.arithmetic_w);
    }
    let (source_contact_heat, source_contact_states) =
        sources.evaluate(data, &trace, &heat, gate)?;
    for (index, joint) in data.joints.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        let a = data.face_ids[joint.cell][joint.face];
        let b = data.face_ids[joint.neighbor][joint.face ^ 1];
        let response = joint_heat(joint, [trace[a], trace[b]], &source_contact_heat)?;
        trace_residual[a] += response.heat_w;
        trace_residual[b] -= response.heat_w;
    }
    let mut merit = 0.0_f64;
    for &value in cell_residual.iter().chain(&trace_residual) {
        if !value.is_finite() {
            return Err(invalid(
                "original tensor equation defect is unrepresentable",
            ));
        }
        merit = merit.max(value.abs());
    }
    let mut solution = audit_endpoint_with_contacts(
        domain,
        setup,
        data,
        |cell, temperature| coordinate.storage_rate(cell,
            coordinates.as_ref().map_or(temperature, |values| values[cell].value)),
        |joint, temperatures| joint_heat(joint, temperatures, &source_contact_heat)
            .map(|response| (response.heat_w, response.arithmetic_w)),
        carried_power_w,
        temperature,
        faces,
        heat,
        arithmetic,
        gate,
    )?;
    solution.source_contact_states = source_contact_states;
    Ok(Endpoint {
        solution,
        coordinates,
        trace,
        source_contact_heat,
        cell_residual,
        trace_residual,
        merit,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn linearize<C: CellCoordinate>(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    geometry: &[KirchhoffCell],
    data: &Prepared,
    coordinate: &C,
    endpoint: &Endpoint,
    config: &EnergyConfig,
    remaining_iterations: usize,
    gate: &CancelGate,
) -> Result<Direction, ChtError> {
    let unknowns = data.known.len();
    let mut coo = Coo::new(unknowns, unknowns);
    let mut rhs: Vec<_> = endpoint.trace_residual.iter().map(|r| -r).collect();
    let mut recovery = Vec::<Option<KirchhoffRecovery>>::with_capacity(domain.cell_count());
    for c in 0..domain.cell_count() {
        if c.is_multiple_of(256) {
            poll(gate)?;
        }
        let (tangent, compose_temperature) = match cell_law(
            coordinate,
            laws[c],
            c,
            endpoint.value(c),
            endpoint.solution.temperature[c],
        )? {
            CellLaw::Kirchhoff => {
                let center = laws[c]
                    .conductivity_at_temperature(endpoint.solution.temperature[c])?
                    .as_array();
                let mut faces = [[[0.0; 3]; 3]; 6];
                for (f, face) in faces.iter_mut().enumerate() {
                    *face = laws[c]
                        .conductivity_at_temperature(endpoint.trace[data.face_ids[c][f]])?
                        .as_array();
                }
                (geometry[c].tangent(center, faces), !C::IS_TEMPERATURE)
            }
            #[cfg(feature = "enthalpy-coupling")]
            CellLaw::Caloric(response) => {
                let tensor = response.conductivity().as_array();
                let tangent = if data.fixed[c].is_some() {
                    geometry[c].tangent(tensor, [tensor; 6])
                } else {
                    geometry[c].caloric_tangent(
                        tensor,
                        response.derivative_kg_w_m_k_j(),
                        endpoint.solution.temperature[c],
                        data.face_ids[c].map(|id| endpoint.trace[id]),
                        endpoint.temperature_slope(c),
                    )
                };
                (tangent, false)
            }
        };
        let tangent = tangent
            .ok_or_else(|| invalid("original tensor endpoint tangent is unrepresentable"))?;
        let (matrix, load) = if data.fixed[c].is_some() {
            recovery.push(None);
            (tangent.face.map(|row| row.map(|value| -value)), [0.0; 6])
        } else {
            let slope = endpoint.temperature_slope(c);
            let added = data.sinks[c].conductance * slope + coordinate.storage_tangent(c);
            let tangent = if compose_temperature {
                tangent.with_cell_temperature_slope(slope)
                    .ok_or_else(|| invalid("original caloric tensor tangent is unrepresentable"))?
            } else {
                tangent
            };
            let condensed = tangent
                .condense(added, endpoint.cell_residual[c])
                .ok_or_else(|| invalid("nonlinear cell condensation is unrepresentable"))?;
            recovery.push(Some(condensed.recovery()));
            (condensed.matrix, condensed.load)
        };
        for f in 0..6 {
            let i = data.face_ids[c][f];
            if data.known[i].is_some() {
                continue;
            }
            rhs[i] += load[f];
            for g in 0..6 {
                let j = data.face_ids[c][g];
                if data.known[j].is_none() && matrix[f][g] != 0.0 {
                    coo.push(i, j, matrix[f][g]);
                }
            }
            if domain.neighbor(c, f).is_none()
                && let ThermalFace::Convective { h, .. } = setup.faces[f]
            {
                let conductance = domain.face_area(c, f / 2) * h;
                finite_positive("thermal.tensor_conduction.film_conductance", conductance)?;
                coo.push(i, i, conductance);
            }
        }
    }
    for (index, joint) in data.joints.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        let a = data.face_ids[joint.cell][joint.face];
        let b = data.face_ids[joint.neighbor][joint.face ^ 1];
        let response = joint_heat(
            joint, [endpoint.trace[a], endpoint.trace[b]], &endpoint.source_contact_heat,
        )?;
        let [da, db] = response.derivative_w_k;
        coo.push(a, a, da);
        coo.push(a, b, db);
        coo.push(b, b, -db);
        coo.push(b, a, -da);
    }
    for (id, known) in data.known.iter().enumerate() {
        if id.is_multiple_of(4096) {
            poll(gate)?;
        }
        if known.is_some() {
            coo.push(id, id, 1.0);
            rhs[id] = 0.0;
        }
    }
    let matrix = scale_rows(&coo, &mut rhs);
    drop(coo);
    for (id, value) in rhs.iter().enumerate() {
        if id.is_multiple_of(4096) {
            poll(gate)?;
        }
        if !value.is_finite()
            || matrix.row(id).1.iter().any(|a| !a.is_finite())
            || matrix.get(id, id) <= 0.0
        {
            return Err(invalid("original tensor Newton system is unrepresentable"));
        }
    }
    let mut trace = vec![0.0; unknowns];
    let outcome = bicgstab_ilu0(
        "nonlinear tensor Newton correction",
        &matrix,
        &rhs,
        &mut trace,
        config.tolerance,
        remaining_iterations,
        gate,
    )?;
    // Constrained trace increments are exactly zero in the physical update.
    for (id, known) in data.known.iter().enumerate() {
        if known.is_some() {
            trace[id] = 0.0;
        }
    }
    let (residual, scale) = residual_norms(&matrix, &rhs, &trace);
    let relative_residual = if scale > 0.0 {
        residual / scale
    } else {
        residual
    };
    if !relative_residual.is_finite() || relative_residual > config.tolerance {
        return Err(not_converged(
            "nonlinear tensor true linear residual",
            outcome.iterations,
            relative_residual,
            config.tolerance,
        ));
    }
    let mut cell = Vec::with_capacity(recovery.len());
    for (c, saved) in recovery.into_iter().enumerate() {
        if c.is_multiple_of(512) {
            poll(gate)?;
        }
        let increment = match saved {
            None => 0.0,
            Some(saved) => saved
                .cell_increment(data.face_ids[c].map(|id| trace[id]))
                .ok_or_else(|| invalid("nonlinear cell correction is unrepresentable"))?,
        };
        cell.push(increment);
    }
    Ok(Direction {
        cell,
        trace,
        iterations: outcome.iterations,
        nonzeros: matrix.nnz(),
        relative_residual,
    })
}
