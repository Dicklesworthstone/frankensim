//! Full-tensor conduction on the native rectilinear voxel grid.
//!
//! The hybrid finite-volume energy retains cell and face temperatures. Cell
//! unknowns are eliminated locally, leaving one shared trace at perfect
//! interfaces or two traces at a finite thermal contact. The SPD tensor is
//! used in full: a tangential temperature gradient may carry normal heat.
//! This is an opt-in conduction solve, with no advective/compact/source-law
//! inference and no discrete maximum-principle claim for rotated anisotropy.
//!
//! See Droniou et al., M3AS 20 (2010), doi:10.1142/S0218202510004222,
//! equations 2.17--2.21; `local` specializes their energy to Cartesian cells.

use std::collections::BTreeMap;
use fs_exec::CancelGate;
use fs_sparse::Coo;
use super::energy::{PseudoStep, SourceContactState, heat_residual_ratio};
use super::krylov::{bicgstab_ilu0, residual_norms};
use super::{
    ChtError, EnergyBalance, EnergyConfig, SolidMaterial, ThermalFace, ThermalSetup,
    Voxel, VoxelDomain, finite, finite_positive, poll, scale_rows,
};

mod contact;
mod local;
mod nonlinear;
mod sdirk;
mod transient;
use local::LocalOperator;
#[cfg(feature = "enthalpy-coupling")]
pub(super) use nonlinear::{
    CaloricAdaptiveForcing, CaloricAdaptiveTrial, advance_caloric_adaptive_with, embedded_ratio,
    advance_radiating_caloric_prepared, advance_radiating_sdirk2_prepared,
};
#[cfg(feature = "enthalpy-coupling")]
pub(super) use sdirk::stage_clock as caloric_stage_clock;
pub use nonlinear::{
    RotatedSolidConductivity, RotatedSourceConductivity, TensorConductivityLaw, TensorNonlinearConfig,
    TensorNonlinearSdirk2Step, TensorNonlinearSolution, TensorNonlinearStep,
    advance_nonlinear_tensor_conduction, advance_nonlinear_tensor_conduction_sdirk2,
    solve_nonlinear_tensor_conduction,
};
#[cfg(feature = "enthalpy-coupling")]
pub use nonlinear::{
    IsotropicPhaseConductivity, TensorCaloricConductivity, TensorCaloricConductivityPoint,
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
pub use sdirk::{
    TENSOR_SDIRK2_GAMMA, TensorConductionSdirk2Step, advance_tensor_conduction_sdirk2,
};
pub use transient::{TensorConductionStep, advance_tensor_conduction};

/// Frozen physical conductivity in grid Cartesian coordinates, W/(m K).
/// All six independent components are retained. The constructor refuses
/// nonsymmetric, non-finite and numerically non-positive-definite matrices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConductivityTensor([[f64; 3]; 3]);

impl ConductivityTensor {
    /// Admit one symmetric positive-definite tensor without projecting it.
    ///
    /// # Errors
    /// A non-finite, nonsymmetric or non-positive LDL pivot is refused.
    pub fn new(value: [[f64; 3]; 3]) -> Result<Self, ChtError> {
        if !local::is_spd(value) {
            return Err(invalid("conductivity needs a finite symmetric positive-definite 3x3 tensor"));
        }
        Ok(Self(value))
    }

    /// Conductivity in the original Cartesian frame, W/(m K).
    #[must_use]
    pub const fn as_array(self) -> [[f64; 3]; 3] { self.0 }
}

/// Algebraic and independently evaluated physical evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorConductionReport {
    /// Face-trace unknowns, including constrained boundary traces.
    pub unknowns: usize,
    /// Stored nonzeros in the condensed face system.
    pub nonzeros: usize,
    /// Krylov iterations, including any same-matrix residual refinement.
    pub iterations: usize,
    /// Recomputed row-scaled linear residual. For frozen conduction this is
    /// evaluated at the returned trace; for a nonlinear solve it belongs to
    /// the last Newton correction. Original heat acceptance is reported below.
    pub relative_residual: f64,
    /// Worst local cell, interface/contact or boundary heat defect, after
    /// the native arithmetic roundoff allowance, relative to local throughput.
    pub max_heat_residual: f64,
    /// Actual heat leaving solid cells toward fluid cells, W.
    /// This observer sums physical face fluxes, never scalar conductances.
    pub solid_to_fluid_heat_w: f64,
    /// Native heat ledger. Storage, if supplied by a physical step owner,
    /// is excluded here and included in the local endpoint acceptance.
    pub balance: EnergyBalance,
}

/// Cell temperatures, both contact traces and actual physical face heat.
/// Arrays use `Face3::ALL` order. Perfect-interface traces are identical;
/// finite contacts retain both temperatures and their heat-carrying jump.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorConductionSolution {
    /// One temperature per voxel, K.
    pub temperature: Vec<f64>,
    /// Face-plane temperature on each side of each cell, K.
    pub face_temperature_k: Vec<[f64; 6]>,
    /// Conductive heat leaving each cell face, W, including thermostat faces.
    pub outward_heat_w: Vec<[f64; 6]>,
    /// Original source-contact receipts and actual surface heat, in authored
    /// A/B order. Empty for fixed-resistance or perfect contacts.
    pub source_contact_states: Vec<SourceContactState>,
    /// Solver, conservation and heat-transfer evidence.
    pub report: TensorConductionReport,
}

impl TensorConductionSolution {
    /// Actual conductive heat leaving solids into fluid voxels, W.
    #[must_use]
    pub fn solid_to_fluid_heat_w(&self) -> f64 { self.report.solid_to_fluid_heat_w }
}

/// Solve steady conduction with one explicitly resolved frozen tensor per cell.
///
/// Uses the native graded geometry, sources, fixed cells, linear cell sinks,
/// temperature/flux/adiabatic/film boundaries, material-pair contacts and
/// explicit face contacts. Material constants are not substituted for the
/// supplied field. Native k(T) tables, source contacts, eddy additions, compact
/// components and open-flow boundaries require their existing owners and are
/// refused. `config.scheme` is unused because this API has no advection.
///
/// The hybrid energy is affine-consistent and positive definite once anchored;
/// it reduces to native half-cell/harmonic resistances for diagonal tensors.
/// General rotated anisotropy is not an M-matrix and has no maximum-principle
/// guarantee. Returned cell and trace temperatures must remain positive Kelvin.
/// Up to eight same-matrix refinements share `config.max_iterations`; both
/// the true algebraic residual and every physical local balance must pass.
///
/// # Errors
/// Refuses incompatible owners, invalid/unrepresentable data, an unanchored
/// steady problem, non-positive returned temperatures, cancellation or an
/// exhausted algebraic/physical residual budget.
pub fn solve_tensor_conduction(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    tensors: &[ConductivityTensor],
    config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    solve_tensor_conduction_inner(domain, solids, setup, tensors, config, None, gate)
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput { field: "thermal.tensor_conduction", reason: reason.into() }
}

#[derive(Clone, Copy, Default)]
struct Sink {
    conductance: f64,
    relative_load: f64,
    absolute_load: f64,
    count: usize,
}

struct Joint {
    cell: usize,
    neighbor: usize,
    face: usize,
    conductance: f64,
    source: Option<contact::SourceJoint>,
}

struct Prepared {
    reference: f64,
    fixed: Vec<Option<f64>>,
    sinks: Vec<Sink>,
    face_ids: Vec<[usize; 6]>,
    known: Vec<Option<f64>>,
    joints: Vec<Joint>,
}

fn operator(
    domain: &VoxelDomain, tensors: &[ConductivityTensor], cell: usize,
) -> Result<LocalOperator, ChtError> {
    LocalOperator::new(tensors[cell].0, domain.widths(cell))
        .ok_or_else(|| invalid("cell geometry or tensor half-cell coefficients are unrepresentable"))
}

fn validate_sensible_storage(
    domain: &VoxelDomain,
    step: Option<&PseudoStep<'_>>,
    gate: &CancelGate,
) -> Result<(), ChtError> {
    if let Some(step) = step {
        let cells = domain.cell_count();
        if step.coefficient_w_k.len() != cells
            || step.previous.len() != cells
            || !step.junction_coefficient_w_k.is_empty()
        {
            return Err(invalid(
                "storage needs one original temperature/coefficient per cell and no junctions",
            ));
        }
        for c in 0..cells {
            if c.is_multiple_of(4096) {
                poll(gate)?;
            }
            finite_positive("thermal.storage.coefficient_w_k", step.coefficient_w_k[c])?;
        }
    }
    Ok(())
}

fn prepare(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    tensors: &[ConductivityTensor],
    config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>,
    gate: &CancelGate,
) -> Result<Prepared, ChtError> {
    validate_sensible_storage(domain, step, gate)?;
    prepare_geometry(
        domain,
        solids,
        setup,
        tensors,
        config,
        step.map(|s| s.previous),
        None,
        gate,
    )
}

/// Geometry admission only needs the original positive temperature reference
/// and knowledge that physical storage anchors every cell. Its energy
/// coordinate and storage tangent stay with the selected transient owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn prepare_geometry(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    tensors: &[ConductivityTensor],
    config: &EnergyConfig,
    storage_reference_k: Option<&[f64]>,
    sources: Option<&contact::PreparedSources<'_>>,
    gate: &CancelGate,
) -> Result<Prepared, ChtError> {
    poll(gate)?;
    let cells = domain.cell_count();
    cells.checked_mul(36).ok_or_else(|| invalid("tensor stencil size overflows"))?;
    domain.check_materials(solids.len())?;
    finite_positive("thermal.tensor_conduction.tolerance", config.tolerance)?;
    if tensors.len() != cells {
        return Err(invalid("needs exactly one prepared conductivity tensor per cell"));
    }
    if (!setup.source_face_contacts.is_empty() && sources.is_none())
        || !setup.eddy_conductivity_w_m_k.is_empty()
        || !setup.compact_components.is_empty()
        || !setup.compact_heat_storage.is_empty()
    {
        return Err(invalid("source contacts require the nonlinear tensor owner; eddy additions and compact components are unsupported"));
    }
    if !setup.power_w.is_empty() && setup.power_w.len() != cells {
        return Err(invalid("power needs zero entries or one entry per cell"));
    }
    for (index, &power) in setup.power_w.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        finite("thermal.power_w", power)?;
    }
    let mut fixed = vec![None; cells];
    for (index, &(cell, temperature)) in setup.fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        finite_positive("thermal.fixed_temperature", temperature)?;
        if cell >= cells { return Err(invalid("fixed-temperature cell is outside the domain")); }
        if fixed[cell].replace(temperature).is_some() {
            return Err(invalid("a cell has more than one fixed-temperature owner"));
        }
    }
    let mut anchored = !setup.fixed_temperature.is_empty();
    for rule in setup.faces {
        match rule {
            ThermalFace::Adiabatic => {}
            ThermalFace::HeatFlux(q) => finite("thermal.face.heat_flux", q)?,
            ThermalFace::Temperature(t) => {
                finite_positive("thermal.face.temperature", t)?;
                anchored = true;
            }
            ThermalFace::Convective { h, ambient } => {
                finite_positive("thermal.face.h", h)?;
                finite_positive("thermal.face.ambient", ambient)?;
                anchored = true;
            }
            ThermalFace::Inflow { .. } | ThermalFace::Outflow { .. } =>
                return Err(invalid("open-flow boundaries require the advective energy owner")),
        }
    }
    if let Some(previous) = storage_reference_k {
        if previous.len() != cells {
            return Err(invalid(
                "storage needs one original temperature reference per cell",
            ));
        }
        for (c, &temperature) in previous.iter().enumerate() {
            if c.is_multiple_of(4096) {
                poll(gate)?;
            }
            finite_positive("thermal.storage.previous_k", temperature)?;
        }
        anchored = true;
    }
    let reference = setup.fixed_temperature.first().map(|p| p.1)
        .or_else(|| setup.faces.iter().find_map(|rule| match *rule {
            ThermalFace::Temperature(t) => Some(t),
            ThermalFace::Convective { ambient, .. } => Some(ambient),
            _ => None,
        }))
        .or_else(|| setup.cell_sinks.first().map(|s| s.temperature_k))
        .or_else(|| storage_reference_k.and_then(|s| s.first().copied()))
        .unwrap_or(0.0);
    let mut sinks = vec![Sink::default(); cells];
    for (index, sink) in setup.cell_sinks.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        if sink.cell >= cells { return Err(invalid("sink cell is outside the domain")); }
        finite("thermal.cell_sinks.conductance_w_k", sink.conductance_w_k)?;
        finite_positive("thermal.cell_sinks.temperature_k", sink.temperature_k)?;
        if sink.conductance_w_k < 0.0 { return Err(invalid("sink conductance must be non-negative")); }
        let entry = &mut sinks[sink.cell];
        entry.conductance += sink.conductance_w_k;
        entry.relative_load += sink.conductance_w_k * (sink.temperature_k - reference);
        entry.absolute_load += sink.conductance_w_k * sink.temperature_k;
        entry.count += 1;
        anchored |= sink.conductance_w_k > 0.0;
        if !entry.conductance.is_finite() || !entry.relative_load.is_finite()
            || !entry.absolute_load.is_finite()
        { return Err(invalid("accumulated sink load is unrepresentable")); }
    }
    if !anchored {
        return Err(invalid("steady conduction needs a temperature, positive film/sink or fixed-cell anchor"));
    }
    let mut pairs = BTreeMap::new();
    for (index, contact) in setup.contacts.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        let (a, b) = contact.materials;
        let r = contact.resistance_m2_k_w;
        if a == b || usize::from(a) >= solids.len() || usize::from(b) >= solids.len()
            || !r.is_finite() || r < 0.0
        { return Err(invalid("material contacts need distinct declared solids and finite non-negative resistance")); }
        if pairs.insert((a.min(b), a.max(b)), r).is_some() {
            return Err(invalid("duplicate material-pair contact owner"));
        }
    }
    let explicit = super::energy::face_contact::prepare_conduction(domain, setup, gate)?;
    if let Some(sources) = sources {
        for (index, (a, b)) in sources.cell_pairs().enumerate() {
            if index.is_multiple_of(256) { poll(gate)?; }
            let pair_owned = match (domain.voxel_at(a), domain.voxel_at(b)) {
                (Voxel::Solid(a), Voxel::Solid(b)) => pairs.contains_key(&(a.min(b), a.max(b))),
                _ => false,
            };
            if explicit.contains_key(&(a, b)) || pair_owned {
                return Err(invalid("a source contact cannot also have a fixed face or material-pair owner"));
            }
        }
    }
    let mut face_ids = vec![[usize::MAX; 6]; cells];
    let mut known = Vec::new();
    let mut joints = Vec::new();
    for c in 0..cells {
        if c.is_multiple_of(512) { poll(gate)?; }
        operator(domain, tensors, c)?;
        for f in 0..6 {
            if face_ids[c][f] != usize::MAX { continue; }
            face_ids[c][f] = known.len();
            let prescribed = if domain.neighbor(c, f).is_none() {
                match setup.faces[f] { ThermalFace::Temperature(t) => Some(t), _ => None }
            } else { None };
            known.push(prescribed);
            if let Some(n) = domain.neighbor(c, f) {
                let source = sources.and_then(|sources| sources.face(c, n));
                let resistance = source.map(|(_, resistance)| resistance).or_else(|| {
                    explicit.get(&(c.min(n), c.max(n))).copied().or_else(|| {
                        match (domain.voxel_at(c), domain.voxel_at(n)) {
                            (Voxel::Solid(a), Voxel::Solid(b)) => pairs.get(&(a.min(b), a.max(b))).copied(),
                            _ => None,
                        }
                    })
                }).unwrap_or(0.0);
                if resistance > 0.0 {
                    let conductance = domain.face_area(c, f / 2) / resistance;
                    finite_positive("thermal.tensor_conduction.contact_conductance", conductance)?;
                    face_ids[n][f ^ 1] = known.len();
                    known.push(None);
                    joints.push(Joint {
                        cell: c, neighbor: n, face: f, conductance,
                        source: source.map(|(owner, _)| owner),
                    });
                } else {
                    face_ids[n][f ^ 1] = face_ids[c][f];
                }
            }
        }
    }
    Ok(Prepared { reference, fixed, sinks, face_ids, known, joints })
}

fn cell_terms(
    op: &LocalOperator, cell: usize, setup: &ThermalSetup, data: &Prepared,
    step: Option<&PseudoStep<'_>>, carried_power_w: Option<&[f64]>,
) -> Result<(f64, f64), ChtError> {
    let mut added = data.sinks[cell].conductance;
    let mut load = data.sinks[cell].relative_load
        + setup.power_w.get(cell).copied().unwrap_or(0.0)
        + carried_power_w.map_or(0.0, |carry| carry[cell]);
    if let Some(step) = step {
        added += step.coefficient_w_k[cell];
        load += step.coefficient_w_k[cell] * (step.previous[cell] - data.reference);
    }
    if !added.is_finite() || !(op.cell_diagonal + added).is_finite() || !load.is_finite() {
        return Err(invalid("cell source, sink or storage load is unrepresentable"));
    }
    Ok((added, load))
}

/// Shared spatial owner: physical time drivers supply only free-cell storage
/// and the original old state. Fixed-cell storage belongs to their ledger.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
pub(crate) fn solve_tensor_conduction_inner(
    domain: &VoxelDomain, solids: &[SolidMaterial], setup: &ThermalSetup,
    tensors: &[ConductivityTensor], config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>, gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    solve_tensor_conduction_stage(domain, solids, setup, tensors, config, step, None, gate)
}

/// A physical stage may retain heat from an earlier stage. This immutable
/// carry changes the storage equation only; it is never generated heat.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn solve_tensor_conduction_stage(
    domain: &VoxelDomain, solids: &[SolidMaterial], setup: &ThermalSetup,
    tensors: &[ConductivityTensor], config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>, carried_power_w: Option<&[f64]>, gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    poll(gate)?;
    if solids.iter().any(SolidMaterial::has_temperature_dependent_conductivity) {
        return Err(invalid("frozen tensors cannot also own native temperature-dependent conductivity"));
    }
    let data = prepare(domain, solids, setup, tensors, config, step, gate)?;
    if let Some(carry) = carried_power_w {
        if step.is_none() || carry.len() != domain.cell_count() {
            return Err(invalid("stage heat carry needs original storage and one value per cell"));
        }
        for (cell, &value) in carry.iter().enumerate() {
            if cell.is_multiple_of(4096) { poll(gate)?; }
            finite("thermal.tensor_conduction.carried_power_w", value)?;
        }
    }
    let unknowns = data.known.len();
    let mut coo = Coo::new(unknowns, unknowns);
    let mut rhs = vec![0.0; unknowns];
    for c in 0..domain.cell_count() {
        if c.is_multiple_of(512) { poll(gate)?; }
        let op = operator(domain, tensors, c)?;
        let (matrix, face_load) = if let Some(t) = data.fixed[c] {
            (op.matrix, op.cell_link.map(|b| b * (t - data.reference)))
        } else {
            let (added, load) = cell_terms(&op, c, setup, &data, step, carried_power_w)?;
            let diagonal = op.cell_diagonal + added;
            (op.schur(added), op.cell_link.map(|b| b / diagonal * load))
        };
        for f in 0..6 {
            let i = data.face_ids[c][f];
            if data.known[i].is_some() { continue; }
            rhs[i] += face_load[f];
            for g in 0..6 {
                let j = data.face_ids[c][g];
                if let Some(t) = data.known[j] {
                    rhs[i] -= matrix[f][g] * (t - data.reference);
                } else if matrix[f][g] != 0.0 {
                    coo.push(i, j, matrix[f][g]);
                }
            }
            if domain.neighbor(c, f).is_none() {
                let area = domain.face_area(c, f / 2);
                match setup.faces[f] {
                    ThermalFace::HeatFlux(q) => rhs[i] += area * q,
                    ThermalFace::Convective { h, ambient } => {
                        let conductance = area * h;
                        finite_positive("thermal.tensor_conduction.film_conductance", conductance)?;
                        coo.push(i, i, conductance);
                        rhs[i] += conductance * (ambient - data.reference);
                    }
                    _ => {}
                }
            }
        }
    }
    for (index, joint) in data.joints.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        let i = data.face_ids[joint.cell][joint.face];
        let j = data.face_ids[joint.neighbor][joint.face ^ 1];
        let g = joint.conductance;
        coo.push(i, i, g); coo.push(i, j, -g);
        coo.push(j, j, g); coo.push(j, i, -g);
    }
    for (i, known) in data.known.iter().enumerate() {
        if i.is_multiple_of(4096) { poll(gate)?; }
        if let Some(t) = known {
            coo.push(i, i, 1.0);
            rhs[i] = t - data.reference;
        }
    }
    let matrix = scale_rows(&coo, &mut rhs);
    drop(coo);
    for (i, value) in rhs.iter().enumerate() {
        if i.is_multiple_of(4096) { poll(gate)?; }
        if !value.is_finite() || matrix.row(i).1.iter().any(|a| !a.is_finite())
            || matrix.get(i, i) <= 0.0
        { return Err(invalid("condensed face system is unrepresentable")); }
    }
    let mut trace = vec![0.0; unknowns];
    let mut iterations = 0;
    let mut target = config.tolerance;
    let mut last_residual = f64::INFINITY;
    for _ in 0..8 {
        poll(gate)?;
        let outcome = bicgstab_ilu0(
            "tensor conduction", &matrix, &rhs, &mut trace, target,
            config.max_iterations - iterations, gate,
        )?;
        iterations += outcome.iterations;
        let mut solution = recover(domain, setup, tensors, &data, step, carried_power_w, &mut trace, gate)?;
        let (residual, scale) = residual_norms(&matrix, &rhs, &trace);
        let relative = if scale > 0.0 { residual / scale } else { residual };
        if !relative.is_finite() {
            return Err(invalid("true tensor residual norm is unrepresentable"));
        }
        last_residual = relative.max(solution.report.max_heat_residual);
        if last_residual.is_finite() && last_residual <= config.tolerance {
            solution.report.unknowns = unknowns;
            solution.report.nonzeros = matrix.nnz();
            solution.report.iterations = iterations;
            solution.report.relative_residual = relative;
            poll(gate)?;
            return Ok(solution);
        }
        if iterations == config.max_iterations { break; }
        target = 0.1 * target.min(outcome.relative_residual);
        if target == 0.0 || !target.is_finite() { break; }
    }
    Err(ChtError::SolverNotConverged {
        system: "tensor conduction physical heat", iterations,
        relative_residual: last_residual, tolerance: config.tolerance,
    })
}

#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn recover(
    domain: &VoxelDomain, setup: &ThermalSetup, tensors: &[ConductivityTensor],
    data: &Prepared, step: Option<&PseudoStep<'_>>, carried_power_w: Option<&[f64]>,
    trace: &mut [f64], gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    let cells = domain.cell_count();
    let mut surface = Vec::with_capacity(trace.len());
    for (i, value) in trace.iter_mut().enumerate() {
        if i.is_multiple_of(4096) { poll(gate)?; }
        let t = data.known[i].unwrap_or(data.reference + *value);
        finite_positive("thermal.tensor_conduction.face_temperature_k", t)?;
        *value = t - data.reference;
        surface.push(t);
    }
    let mut temperature = Vec::with_capacity(cells);
    let mut face_temperature_k = Vec::with_capacity(cells);
    let mut outward_heat_w = Vec::with_capacity(cells);
    let mut arithmetic = Vec::with_capacity(cells);
    for c in 0..cells {
        if c.is_multiple_of(512) { poll(gate)?; }
        let op = operator(domain, tensors, c)?;
        let faces = data.face_ids[c].map(|id| surface[id]);
        let t = if let Some(t) = data.fixed[c] { t } else {
            let (added, load) = cell_terms(&op, c, setup, data, step, carried_power_w)?;
            let weighted = op.cell_link.iter().zip(data.face_ids[c])
                .fold(load, |sum, (b, id)| b.mul_add(trace[id], sum));
            data.reference + weighted / (op.cell_diagonal + added)
        };
        finite_positive("thermal.tensor_conduction.temperature_k", t)?;
        temperature.push(t);
        face_temperature_k.push(faces);
        outward_heat_w.push(op.outward_heat(t, faces));
        arithmetic.push(op.heat_arithmetic(t, faces));
    }
    audit_endpoint(domain, setup, data, step, carried_power_w, temperature,
        face_temperature_k, outward_heat_w, arithmetic, gate)
}

/// One physical heat owner for frozen and original-law tensor endpoints.
/// The nonlinear caller supplies its actual integral heat, without recovering
/// a different temperature from a frozen surrogate.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn audit_endpoint(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    data: &Prepared,
    step: Option<&PseudoStep<'_>>,
    carried_power_w: Option<&[f64]>,
    temperature: Vec<f64>,
    face_temperature_k: Vec<[f64; 6]>,
    outward_heat_w: Vec<[f64; 6]>,
    arithmetic: Vec<[f64; 6]>,
    gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    audit_endpoint_with_storage(
        domain,
        setup,
        data,
        |cell, temperature| {
            Ok(step.map_or(StorageRate::default(), |step| StorageRate {
                power_w: step.coefficient_w_k[cell] * (temperature - step.previous[cell]),
                arithmetic_w: step.coefficient_w_k[cell]
                    * (temperature.abs() + step.previous[cell].abs()),
            }))
        },
        carried_power_w,
        temperature,
        face_temperature_k,
        outward_heat_w,
        arithmetic,
        gate,
    )
}

/// Actual physical storage rate and its arithmetic scale, independent of the
/// coordinate used by the nonlinear solve. Enthalpy supplies centered energy
/// change here, never an apparent heat capacity or a datum-sized allowance.
#[derive(Clone, Copy, Default)]
struct StorageRate {
    power_w: f64,
    arithmetic_w: f64,
}

#[allow(clippy::too_many_arguments)]
fn audit_endpoint_with_storage(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    data: &Prepared,
    storage: impl Fn(usize, f64) -> Result<StorageRate, ChtError>,
    carried_power_w: Option<&[f64]>,
    temperature: Vec<f64>,
    face_temperature_k: Vec<[f64; 6]>,
    outward_heat_w: Vec<[f64; 6]>,
    arithmetic: Vec<[f64; 6]>,
    gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    audit_endpoint_with_contacts(
        domain, setup, data, storage,
        |joint, [a, b]| {
            if joint.source.is_some() {
                return Err(invalid("source contact audit requires its original physical response"));
            }
            Ok((joint.conductance * (a - b), joint.conductance * (a.abs() + b.abs())))
        },
        carried_power_w, temperature, face_temperature_k, outward_heat_w, arithmetic, gate,
    )
}

/// One physical heat audit for fixed contacts and cached ORIGINAL source laws.
/// The contact closure supplies constitutive heat and its arithmetic scale;
/// actual bulk face heat is checked independently against that same jump law.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit_endpoint_with_contacts(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    data: &Prepared,
    storage: impl Fn(usize, f64) -> Result<StorageRate, ChtError>,
    contact: impl Fn(&Joint, [f64; 2]) -> Result<(f64, f64), ChtError>,
    carried_power_w: Option<&[f64]>,
    temperature: Vec<f64>,
    face_temperature_k: Vec<[f64; 6]>,
    outward_heat_w: Vec<[f64; 6]>,
    arithmetic: Vec<[f64; 6]>,
    gate: &CancelGate,
) -> Result<TensorConductionSolution, ChtError> {
    let cells = domain.cell_count();
    let mut maximum = 0.0_f64;
    let (mut source, mut fixed, mut boundary, mut sinks, mut scale) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut solid_fluid = 0.0;
    for c in 0..cells {
        if c.is_multiple_of(512) { poll(gate)?; }
        let t = temperature[c];
        let sink = data.sinks[c];
        let sink_heat = sink.conductance * (t - data.reference) - sink.relative_load;
        let power = setup.power_w.get(c).copied().unwrap_or(0.0);
        let outward: f64 = outward_heat_w[c].iter().sum();
        let mut residual = outward + sink_heat - power;
        let mut throughput = outward_heat_w[c].iter().map(|q| q.abs()).sum::<f64>()
            + sink_heat.abs() + power.abs();
        let mut rounding = arithmetic[c].iter().sum::<f64>()
            + sink.conductance * t.abs() + sink.absolute_load + power.abs();
        sinks += sink_heat;
        scale += sink_heat.abs();
        if data.fixed[c].is_some() {
            fixed += outward + sink_heat;
        } else {
            source += power;
            scale += power.abs();
            let stored = storage(c, t)?;
            finite("thermal.tensor_conduction.storage_power_w", stored.power_w)?;
            finite("thermal.tensor_conduction.storage_arithmetic_w", stored.arithmetic_w)?;
            if stored.arithmetic_w < 0.0 {
                return Err(invalid("physical storage arithmetic scale must be nonnegative"));
            }
            residual += stored.power_w;
            throughput += stored.power_w.abs();
            rounding += stored.arithmetic_w;
            if let Some(carry) = carried_power_w {
                residual -= carry[c];
                throughput += carry[c].abs();
                rounding += carry[c].abs();
            }
            maximum = maximum.max(heat_residual_ratio(
                residual, throughput, rounding, sink.count,
            )?);
        }
        for f in 0..6 {
            let q = outward_heat_w[c][f];
            if let Some(n) = domain.neighbor(c, f) {
                if matches!(domain.voxel_at(c), Voxel::Solid(_)) && domain.is_fluid(n) {
                    solid_fluid += q;
                }
                if c < n {
                    let other = outward_heat_w[n][f ^ 1];
                    maximum = maximum.max(heat_residual_ratio(
                        q + other, q.abs() + other.abs(),
                        arithmetic[c][f] + arithmetic[n][f ^ 1], 0,
                    )?);
                }
            } else {
                // Match the native fixed-cell ledger: its exterior heat is
                // served by its thermostat, outside the free-cell heat balance.
                if data.fixed[c].is_some() { fixed -= q; } else { boundary += q; }
                scale += q.abs();
                let area = domain.face_area(c, f / 2);
                let (expected, bath_arithmetic) = match setup.faces[f] {
                    ThermalFace::Adiabatic => (0.0, 0.0),
                    ThermalFace::HeatFlux(into) => (-area * into, (area * into).abs()),
                    ThermalFace::Convective { h, ambient } => (
                        area * h * (face_temperature_k[c][f] - ambient),
                        area * h * (face_temperature_k[c][f].abs() + ambient.abs()),
                    ),
                    ThermalFace::Temperature(_) => continue,
                    ThermalFace::Inflow { .. } | ThermalFace::Outflow { .. } => unreachable!(),
                };
                maximum = maximum.max(heat_residual_ratio(
                    q - expected, q.abs() + expected.abs(),
                    arithmetic[c][f] + bath_arithmetic, 0,
                )?);
            }
        }
    }
    for (index, joint) in data.joints.iter().enumerate() {
        if index.is_multiple_of(4096) { poll(gate)?; }
        let a = joint.cell; let b = joint.neighbor; let f = joint.face;
        let ta = face_temperature_k[a][f];
        let tb = face_temperature_k[b][f ^ 1];
        let (jump_heat, contact_arithmetic) = contact(joint, [ta, tb])?;
        let q = outward_heat_w[a][f];
        maximum = maximum.max(heat_residual_ratio(
            q - jump_heat, q.abs() + jump_heat.abs(),
            arithmetic[a][f] + contact_arithmetic, 0,
        )?);
    }
    scale += fixed.abs();
    let residual_w = source + fixed - boundary - sinks;
    if [source, fixed, boundary, sinks, scale, residual_w, solid_fluid]
        .iter().any(|v| !v.is_finite())
    { return Err(invalid("physical tensor heat ledger is unrepresentable")); }
    let balance = EnergyBalance {
        source_w: source, fixed_cell_injection_w: fixed,
        boundary_outflow_w: boundary, advective_outflow_w: 0.0,
        sink_outflow_w: sinks, residual_w,
        relative_residual: if scale > 0.0 { residual_w.abs() / scale } else { 0.0 },
    };
    Ok(TensorConductionSolution {
        temperature, face_temperature_k, outward_heat_w,
        source_contact_states: Vec::new(),
        report: TensorConductionReport {
            unknowns: 0, nonzeros: 0, iterations: 0, relative_residual: 0.0,
            max_heat_residual: maximum, solid_to_fluid_heat_w: solid_fluid, balance,
        },
    })
}
