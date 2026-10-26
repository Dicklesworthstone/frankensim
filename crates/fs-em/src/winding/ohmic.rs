//! Material Ohmic loss on the existing field-derived winding/circuit owner.
//!
//! Each represented conductor uses the SAME material card and signed unit-current
//! profile as the magnetic solve. Integrating rho_e * J_unit^2 over its actual
//! extruded mesh volume gives the winding's intrinsic series resistance. It is
//! stamped into the original flux-coordinate dissipation block; the complete
//! mutual-inductance storage, Kirchhoff structure, ports and graph resistors stay
//! with fs-phs. Every heat receipt comes from that one native midpoint step.
//!
//! This is a fixed-geometry, disjoint stranded-conductor reduction. End resistance,
//! skin/proximity/eddy loss, movement and variable magnetic storage are not inferred.
//! Changing the material data before preparation changes both field and resistance;
//! during evolution only temperature may change, with permeability held constant.
//! Explicit material replacement at an accepted time uses `rebind_material_candidate`:
//! conserved flux/charge, refreshed fields/current/resistance and external work.

mod rebind;
pub use rebind::{WindingMaterialRebind, WindingMaterialRebindControls};

use super::{WindingError, circuit::WindingCircuit};
use crate::planar::MagneticError;
use fs_material::{
    conductor::ResolvedConductorState,
    field_storage::{FieldStorageError, admit_fixed_magnetic_refresh},
    state_point::ResolvedScalarProperty,
};
use fs_phs::{
    DescriptorPortHamiltonian, PhsError, QuadraticStorage, StepRecord,
    circuit::{Branch, CircuitDae, CircuitError},
};

/// One original mesh triangle's complete material/electrical source.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingConductor {
    /// Original triangle index, not winding-local order.
    pub triangle: usize,
    /// Complete original source snapshot, including magnetic/thermal properties.
    pub conductor: ResolvedConductorState,
}

/// Absolute physical accounting limits; not a time-discretization certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OhmicWindingControls {
    /// Agreement of independently reconstructed material and circuit power [W].
    pub power_tolerance_w: f64,
    /// Agreement of material heat, native dissipation, storage and port work [J].
    pub energy_tolerance_j: f64,
}

/// Refusals never alter accepted fluxes, charges, sources or time.
#[derive(Debug)]
pub enum OhmicWindingError {
    /// Incomplete correspondence, unsupported source/model or bad arithmetic.
    Invalid(&'static str),
    /// Existing magnetic owner refusal.
    Winding(WindingError),
    /// Existing Kirchhoff/constraint/step owner refusal.
    Circuit(CircuitError),
    /// Existing descriptor structure/storage admission refusal.
    Phs(PhsError),
    /// Original source permeability admission refusal.
    Magnetic(FieldStorageError),
    /// Independently reconstructed power/energy does not close.
    Balance(&'static str),
    /// Cancellation before publication.
    Cancelled,
}
impl From<WindingError> for OhmicWindingError {
    fn from(error: WindingError) -> Self {
        match error {
            WindingError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Winding(other),
        }
    }
}
impl From<CircuitError> for OhmicWindingError {
    fn from(error: CircuitError) -> Self { Self::Circuit(error) }
}
impl From<PhsError> for OhmicWindingError {
    fn from(error: PhsError) -> Self { Self::Phs(error) }
}
impl From<FieldStorageError> for OhmicWindingError {
    fn from(error: FieldStorageError) -> Self { Self::Magnetic(error) }
}
impl core::fmt::Display for OhmicWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material winding: {self:?}")
    }
}
impl std::error::Error for OhmicWindingError {}

/// Cloneable accepted electrical history; immutable magnetic geometry stays shared.
/// Only fluxes/charges are stored energy. Algebraic coordinates in `state` are the
/// last native solution and are reinitialized against each candidate's held ports.
#[derive(Clone)]
pub struct OhmicWindingCircuit<'c, 'm> {
    magnetic: &'c WindingCircuit<'m>,
    conductors: Vec<WindingConductor>,
    conductor_windings: Vec<usize>,
    resistance_ohm: Vec<f64>,
    state: Vec<f64>,
    time_s: f64,
}

/// Spatial electrical/heat quantities at the owning scheme's common current:
/// native midpoint for linear windings, endpoint for nonlinear BE windings.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingCellHeat {
    /// Original complete source, with actual resolved physical query and receipts.
    pub source: WindingConductor,
    /// Owning winding, in original unit-profile order.
    pub winding: usize,
    /// Actual extruded triangle volume [m3].
    pub volume_m3: f64,
    /// Signed out-of-plane current density [A/m2].
    pub current_density_a_m2: f64,
    /// Signed Ohmic electric field rho_e*J [V/m].
    pub electric_field_v_m: f64,
    /// Volume-integrated irreversible power [W].
    pub power_w: f64,
    /// Heat over the physical step [J].
    pub heat_j: f64,
}

/// An explicit graph resistor keeps its authored R and separate heat receiver.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingGraphResistorHeat {
    /// Index in the original full graph, not resistor-only order.
    pub branch: usize,
    /// Signed first-to-second terminal voltage [V].
    pub voltage_v: f64,
    /// Signed first-to-second terminal current [A].
    pub current_a: f64,
    /// Actual graph-resistor power [W].
    pub power_w: f64,
    /// Heat over the physical step [J].
    pub heat_j: f64,
}

/// One original descriptor solve with a complete spatial/graph heat partition.
#[derive(Clone, Debug)]
pub struct OhmicWindingStep {
    pub native: StepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub midpoint_currents_a: Vec<f64>,
    pub cells: Vec<WindingCellHeat>,
    pub graph_resistors: Vec<WindingGraphResistorHeat>,
    /// Sum of intrinsic winding and every explicit graph-resistor heat [J].
    pub heat_j: f64,
    /// Storage change plus total heat minus native source work [J].
    pub energy_residual_j: f64,
}

/// Proposal only: publish with every receiving thermal candidate atomically.
pub struct OhmicWindingCandidate<'c, 'm> {
    circuit: OhmicWindingCircuit<'c, 'm>,
    step: OhmicWindingStep,
}
impl<'c, 'm> OhmicWindingCandidate<'c, 'm> {
    pub fn circuit(&self) -> &OhmicWindingCircuit<'c, 'm> { &self.circuit }
    pub fn step(&self) -> &OhmicWindingStep { &self.step }
    pub fn into_parts(self) -> (OhmicWindingCircuit<'c, 'm>, OhmicWindingStep) {
        (self.circuit, self.step)
    }
}

fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), OhmicWindingError> {
    if cancelled() { Err(OhmicWindingError::Cancelled) } else { Ok(()) }
}
fn admit_electrical_law(property: &ResolvedScalarProperty) -> Result<(), OhmicWindingError> {
    // Use the same canonical source-axis and applicability policy as magnetic
    // storage and thermal material owners. Pressure/composition data remain
    // genuine fixed external conditions; dynamic field/rate laws still refuse.
    if !property.is_static_at_fixed_external_conditions() {
        return Err(OhmicWindingError::Invalid(
            "winding loss requires a static Ohmic scalar with only fixed external or owned temperature support",
        ));
    }
    Ok(())
}

fn bind(
    magnetic: &WindingCircuit<'_>,
    sources: Vec<WindingConductor>,
    previous: Option<&[WindingConductor]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingConductor>, Vec<usize>, Vec<f64>), OhmicWindingError> {
    use OhmicWindingError::Invalid;
    let basis = magnetic.basis();
    let problem = basis.problem();
    if magnetic.graph().branches.iter().any(|&(a, b, _)| a == b) {
        return Err(Invalid("invalid winding conductor capacity or self-loop circuit branch"));
    }
    bind_sources(problem.mesh(), basis.profiles_per_m2(), sources, previous, |source| {
        if problem.recoil_materials()[source.triangle].is_some() {
            if &problem.materials()[source.triangle] != source.conductor.material() {
                return Err(Invalid("conducting recoil material requires its exact original fixed-state source bundle"));
            }
        } else {
            admit_fixed_magnetic_refresh(&problem.materials()[source.triangle], source.conductor.material())
                .map_err(|e| match e {
                    FieldStorageError::Invalid(why) => OhmicWindingError::Invalid(why),
                    other => OhmicWindingError::Magnetic(other),
                })?;
        }
        Ok(())
    }, cancelled)
}

// Shared by constant-inductance and nonlinear-field circuits. Only magnetic
// source admission differs; spatial ownership, source continuity and rho*J^2
// resistance are computed once here.
pub(super) fn bind_sources(
    mesh: &fs_rep_mesh::TriComplex2,
    profiles: &[Vec<f64>],
    mut sources: Vec<WindingConductor>,
    previous: Option<&[WindingConductor]>,
    mut admit_material: impl FnMut(&WindingConductor) -> Result<(), OhmicWindingError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingConductor>, Vec<usize>, Vec<f64>), OhmicWindingError> {
    use OhmicWindingError::Invalid;
    let cells = mesh.faces().len();
    if sources.len() > cells || profiles.iter().any(|p| p.len() != cells) {
        return Err(Invalid("invalid winding conductor or profile capacity"));
    }
    let mut owners = vec![None; cells];
    for (winding, profile) in profiles.iter().enumerate() {
        for (triangle, &current) in profile.iter().enumerate() {
            poll(cancelled)?;
            if current != 0.0 && owners[triangle].replace(winding).is_some() {
                return Err(Invalid("overlapping current-carrying windings require an explicit conductor fill/mixture model"));
            }
        }
    }
    if sources.len() != owners.iter().filter(|o| o.is_some()).count() {
        return Err(Invalid("every current-carrying triangle requires exactly one complete conductor source"));
    }
    sources.sort_by_key(|s| s.triangle);
    if previous.is_some_and(|p| p.len() != sources.len()) {
        return Err(Invalid("winding source correspondence changed"));
    }
    let mut windings = Vec::with_capacity(sources.len());
    let mut resistance = vec![0.0; profiles.len()];
    for (index, source) in sources.iter().enumerate() {
        poll(cancelled)?;
        if source.triangle >= cells || index > 0 && sources[index - 1].triangle == source.triangle {
            return Err(Invalid("duplicate or out-of-range winding conductor triangle"));
        }
        let winding = owners[source.triangle].ok_or(Invalid("conductor binding names an inactive triangle"))?;
        if previous.is_some_and(|p| p[index].triangle != source.triangle
            || p[index].conductor.source_basis() != source.conductor.source_basis()
            || {
                let before = p[index].conductor.material();
                let after = source.conductor.material();
                before.card_identity() != after.card_identity()
                    || before.material() != after.material()
                    || before.axis_quantities() != after.axis_quantities()
                    || before.properties().len() != after.properties().len()
                    || before.properties().iter().zip(after.properties()).any(|(a,b)|
                        a.requirement() != b.requirement() || a.source_axes() != b.source_axes()
                            || a.answer().receipt.selected != b.answer().receipt.selected)
            })
        {
            return Err(Invalid("winding electrical source basis or triangle correspondence changed"));
        }
        admit_material(source)?;
        admit_electrical_law(source.conductor.source_property())?;
        let j = profiles[winding][source.triangle];
        let volume = mesh.face_measure(source.triangle).ok_or(Invalid("winding triangle has no physical mesh volume"))?;
        let contribution = source.conductor.resistivity_ohm_m() * j * j * volume;
        if !positive(volume) || !positive(contribution) {
            return Err(Invalid("spatial winding resistance is not representable and positive"));
        }
        resistance[winding] += contribution;
        if !positive(resistance[winding]) { return Err(Invalid("summed winding resistance is not representable")); }
        windings.push(winding);
    }
    if resistance.iter().any(|&r| !positive(r)) { return Err(Invalid("every winding needs represented finite positive resistance")); }
    Ok((sources, windings, resistance))
}

// Preserve the original full coupled magnetic Q and every J/G/R graph stamp.
// R on flux coordinates acts on their conjugate currents L^-1*lambda, not on
// lambda_i/L_ii. Stamping a shunt conductance would implement different physics.
fn lower(magnetic: &WindingCircuit<'_>, resistance: &[f64]) -> Result<CircuitDae, OhmicWindingError> {
    let original = magnetic.dae();
    let n = original.system.state_dim();
    let nd = original.system.differential_dim();
    let matrix = original.system.storage().exact_quadratic_matrix().ok_or(
        OhmicWindingError::Invalid("native winding storage must declare its complete quadratic form"),
    )?;
    let mut q = vec![0.0; nd * nd];
    for row in 0..nd { q[row * nd..(row + 1) * nd].copy_from_slice(&matrix[row * n..row * n + nd]); }
    let (j, r, g) = original.system.structure();
    let mut r = r.to_vec();
    for (&coordinate, &ohms) in magnetic.flux_coordinates().iter().zip(resistance) {
        r[coordinate * n + coordinate] += ohms;
    }
    let system = DescriptorPortHamiltonian::new(n, nd, original.ports,
        j.to_vec(), r, g.to_vec(), Box::new(QuadraticStorage::new(q, nd)?))?;
    Ok(CircuitDae {
        system,
        flux_index: original.flux_index.clone(),
        charge_index: original.charge_index.clone(),
        node_potential_index: original.node_potential_index.clone(),
        ports: original.ports,
    })
}

impl<'c, 'm> OhmicWindingCircuit<'c, 'm> {
    /// Bind every active triangle to the same prepared magnetic material and
    /// initialize original winding-order currents and capacitor-branch charges.
    /// Intrinsic loss covers only the represented extruded conductor geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn new(magnetic: &'c WindingCircuit<'m>, sources: Vec<WindingConductor>,
        time_s: f64, currents_a: &[f64], charges_c: &[f64], ports: &[f64],
        cancelled: &mut impl FnMut() -> bool) -> Result<Self, OhmicWindingError>
    {
        poll(cancelled)?;
        if !time_s.is_finite() || time_s < 0.0 { return Err(OhmicWindingError::Invalid("finite nonnegative initial time required")); }
        let (conductors, conductor_windings, resistance_ohm) = bind(magnetic, sources, None, cancelled)?;
        let dae = lower(magnetic, &resistance_ohm)?;
        let initial = magnetic.consistent_initial_state(currents_a, charges_c, ports, cancelled)?;
        let state = dae.consistent_initial_state(&initial, ports)?;
        magnetic.admit_state(&state, cancelled)?;
        let h = dae.system.hamiltonian(&state);
        if !h.is_finite() || h < 0.0 || state.iter().any(|v| !v.is_finite()) {
            return Err(OhmicWindingError::Invalid("initial electrical state is not representable"));
        }
        poll(cancelled)?;
        Ok(Self { magnetic, conductors, conductor_windings, resistance_ohm, state, time_s })
    }
    pub fn magnetic(&self) -> &'c WindingCircuit<'m> { self.magnetic }
    /// Sources in stable original triangle order.
    pub fn conductors(&self) -> &[WindingConductor] { &self.conductors }
    /// Owning winding in parallel with `conductors()`.
    pub fn conductor_windings(&self) -> &[usize] { &self.conductor_windings }
    pub fn resistance_ohm(&self) -> &[f64] { &self.resistance_ohm }
    pub fn state(&self) -> &[f64] { &self.state }
    pub fn time_s(&self) -> f64 { self.time_s }
    pub fn currents_a(&self) -> Result<Vec<f64>, OhmicWindingError> { Ok(self.magnetic.currents_a(&self.state)?) }

    /// Evaluate refreshed material resistance on this SAME accepted history.
    /// Electrical midpoint integration is owned by the original descriptor;
    /// resistance is caller-resolved and held for this trial. Coupled thermal
    /// iteration supplies its actual candidate temperature and deposits this
    /// exact heat once. No accepted state is mutated, even on late cancellation.
    pub fn step_to_candidate(&self, target_s: f64, sources: Vec<WindingConductor>,
        ports: &[f64], controls: OhmicWindingControls,
        cancelled: &mut impl FnMut() -> bool) -> Result<OhmicWindingCandidate<'c, 'm>, OhmicWindingError>
    {
        use OhmicWindingError::{Invalid, Balance};
        poll(cancelled)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite() || !positive(dt)
            || !positive(controls.power_tolerance_w) || !positive(controls.energy_tolerance_j)
            || ports.len() != self.magnetic.dae().ports || ports.iter().any(|v| !v.is_finite())
        {
            return Err(Invalid("positive representable interval/tolerances and complete finite ports required"));
        }
        let (conductors, conductor_windings, resistance_ohm) = bind(self.magnetic, sources, Some(&self.conductors), cancelled)?;
        let dae = lower(self.magnetic, &resistance_ohm)?;
        let nd = dae.system.differential_dim();
        let mut seed = vec![0.0; self.state.len()];
        seed[..nd].copy_from_slice(&self.state[..nd]);
        let x0 = dae.consistent_initial_state(&seed, ports)?;
        for i in 0..nd {
            if (x0[i] - self.state[i]).abs() > 64.0 * f64::EPSILON * self.state[i].abs() {
                return Err(Invalid("constraint initialization changed conserved winding/capacitor history"));
            }
        }
        poll(cancelled)?;
        // Native flux coordinates remain delta_lambda=L*i. Permanent linkage
        // is a separate fixed offset; the same whole-field validator used by
        // field sampling checks both actual endpoints before publication.
        self.magnetic.admit_state(&x0, cancelled)?;
        let (native, _) = dae.step_audited(&x0, ports, dt)?;
        self.magnetic.admit_state(&native.x, cancelled)?;
        let stored_before_j = dae.system.hamiltonian(&x0);
        let stored_after_j = dae.system.hamiltonian(&native.x);
        if native.x.iter().chain(&native.y).any(|v| !v.is_finite())
            || [stored_before_j, stored_after_j, native.delta_h, native.dissipated, native.supplied,
                native.solver_residual].iter().any(|v| !v.is_finite())
            || stored_before_j < 0.0 || stored_after_j < 0.0 || native.dissipated < 0.0
        { return Err(Balance("nonfinite or negative native storage/loss")); }
        let midpoint: Vec<_> = x0.iter().zip(&native.x).map(|(&a, &b)| f64::midpoint(a, b)).collect();
        let midpoint_currents_a = self.magnetic.currents_a(&midpoint)?;
        let (cells, winding_heat) = spatial_heat(
            self.magnetic.basis().problem().mesh(), self.magnetic.basis().profiles_per_m2(),
            &conductors, &conductor_windings, &resistance_ohm, &midpoint_currents_a,
            dt, controls.power_tolerance_w, cancelled,
        )?;
        let (graph_resistors, graph_heat) = graph_heat(
            self.magnetic.graph(), &dae, &midpoint, dt, controls.power_tolerance_w, cancelled,
        )?;
        let heat_j = winding_heat + graph_heat;
        let energy_residual_j = stored_after_j - stored_before_j + heat_j - native.supplied;
        if !heat_j.is_finite() || !energy_residual_j.is_finite()
            || (heat_j - native.dissipated).abs() > controls.energy_tolerance_j
            || ((stored_after_j - stored_before_j) - native.delta_h).abs() > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j
        { return Err(Balance("spatial and graph heat, original storage and source work do not close")); }
        let circuit = Self { magnetic: self.magnetic, conductors, conductor_windings, resistance_ohm,
            state: native.x.clone(), time_s: target_s };
        let step = OhmicWindingStep { native, stored_before_j, stored_after_j,
            midpoint_currents_a, cells, graph_resistors, heat_j, energy_residual_j };
        poll(cancelled)?;
        Ok(OhmicWindingCandidate { circuit, step })
    }
}

// Both winding models use the same spatial work/heat owner, with the current
// selected by their declared time scheme (midpoint or backward Euler).
#[allow(clippy::too_many_arguments)]
pub(super) fn spatial_heat(
    mesh: &fs_rep_mesh::TriComplex2,
    profiles: &[Vec<f64>],
    conductors: &[WindingConductor],
    conductor_windings: &[usize],
    resistance_ohm: &[f64],
    currents_a: &[f64],
    dt: f64,
    power_tolerance_w: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingCellHeat>, f64), OhmicWindingError> {
    use OhmicWindingError::{Invalid, Balance};
        let mut cells = Vec::with_capacity(conductors.len());
        let mut winding_power = vec![0.0; resistance_ohm.len()];
        let mut heat_j = 0.0;
        for (source, &winding) in conductors.iter().zip(conductor_windings) {
            poll(cancelled)?;
            let volume_m3 = mesh.face_measure(source.triangle).ok_or(Invalid("winding triangle has no physical mesh volume"))?;
            let current_density_a_m2 = profiles[winding][source.triangle] * currents_a[winding];
            let electric_field_v_m = source.conductor.resistivity_ohm_m() * current_density_a_m2;
            let power_w = electric_field_v_m * current_density_a_m2 * volume_m3;
            let cell_heat = power_w * dt;
            if [current_density_a_m2, electric_field_v_m, power_w, cell_heat].iter().any(|v| !v.is_finite())
                || power_w < 0.0 || current_density_a_m2 != 0.0 && power_w == 0.0
                || power_w != 0.0 && cell_heat == 0.0
            { return Err(Balance("spatial Ohmic power/heat is unrepresentable")); }
            winding_power[winding] += power_w;
            heat_j += cell_heat;
            cells.push(WindingCellHeat { source: source.clone(), winding, volume_m3,
                current_density_a_m2, electric_field_v_m, power_w, heat_j: cell_heat });
        }
        for winding in 0..resistance_ohm.len() {
            let i = currents_a[winding];
            let terminal_power = i * resistance_ohm[winding] * i;
            if !terminal_power.is_finite() || !winding_power[winding].is_finite()
                || (terminal_power - winding_power[winding]).abs() > power_tolerance_w
            { return Err(Balance("integrated material heat disagrees with intrinsic winding series work")); }
        }
    Ok((cells, heat_j))
}

pub(super) fn graph_heat(
    graph: &fs_phs::circuit::CircuitGraph,
    dae: &CircuitDae,
    midpoint: &[f64],
    dt: f64,
    power_tolerance_w: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingGraphResistorHeat>, f64), OhmicWindingError> {
    use OhmicWindingError::Balance;
    let mut heat_j = 0.0;
        let voltage = |node: usize| if node == 0 { 0.0 } else { midpoint[dae.node_potential_index[node - 1]] };
        let mut graph_resistors = Vec::new();
        for (branch, &(a, b, law)) in graph.branches.iter().enumerate() {
            poll(cancelled)?;
            if let Branch::Resistor { ohms } = law {
                let voltage_v = voltage(a) - voltage(b);
                let current_a = voltage_v / ohms;
                let power_w = voltage_v * current_a;
                let q = power_w * dt;
                if [voltage_v, current_a, power_w, q].iter().any(|v| !v.is_finite())
                    || power_w < 0.0 || voltage_v != 0.0 && power_w == 0.0 || power_w != 0.0 && q == 0.0
                    || (current_a * ohms * current_a - power_w).abs() > power_tolerance_w
                { return Err(Balance("original graph resistor work/heat is unrepresentable")); }
                heat_j += q;
                graph_resistors.push(WindingGraphResistorHeat { branch, voltage_v, current_a, power_w, heat_j: q });
            }
        }
    Ok((graph_resistors, heat_j))
}
