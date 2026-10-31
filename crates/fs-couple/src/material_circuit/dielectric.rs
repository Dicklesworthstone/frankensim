//! Source-backed dielectric absorption on the existing Kirchhoff descriptor.
//!
//! One original capacitor supplies instantaneous storage in parallel with
//! source-derived series RC polarization arms. Linear C_infinity is the default;
//! [`nonlinear`] admits an original reversible instantaneous D(E) instead.
//! The descriptor owns all electric history; this module adds no time integrator.
//! Polarization heat comes from the SAME native midpoint efforts as circuit
//! dissipation and is never called bulk Ohmic heat.

use std::sync::Arc;

use fs_material::field_storage::dielectric::{MAX_DIELECTRIC_POLES, ResolvedDebyeCapacitor};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitGraph, assemble_circuit},
};

use super::capacitance::dielectric::MaterialDielectricCapacitor;
use super::{MaterialCircuitError, MaterialCircuitSolve, admit_graph, extended, poll, positive};

/// Source-backed nonlinear instantaneous response with the same Debye memory.
pub mod nonlinear;
/// Conservative material refresh without erasing polarization history.
pub mod rebind;

struct Prepared {
    source: ResolvedDebyeCapacitor,
    graph: CircuitGraph,
    dae: CircuitDae,
    branch: usize,
    instantaneous_charge_index: usize,
    polarization_charge_indices: Vec<usize>,
    loss_branches: Vec<usize>,
    material_resistors: Vec<super::MaterialResistor>,
    reversible_instantaneous: Option<MaterialDielectricCapacitor>,
}

/// Fixed source/geometry and an accepted native circuit state. The compiled
/// descriptor is shared across candidates, so stepping does not reassemble it.
/// Original R/L/source/transformer branches retain their caller-supplied values.
#[derive(Clone)]
pub struct DebyeCircuit {
    prepared: Arc<Prepared>,
    state: Vec<f64>,
}

/// One original polarization process and its physical conversion to heat.
#[derive(Clone, Debug, PartialEq)]
pub struct PolarizationHeat {
    /// Original source process index.
    pub pole: usize,
    /// Actual resistor index in the expanded native graph.
    pub resistor_branch: usize,
    /// Same midpoint resistor effort used by the native dissipation [V].
    pub voltage_drop_v: f64,
    /// Current toward the process's polarization capacitor [A].
    pub current_a: f64,
    /// Irreversible polarization energy released during the interval [J].
    pub heat_j: f64,
}

/// Actual native loss in one original, non-polarization resistor.
#[derive(Clone, Debug, PartialEq)]
pub struct DebyeResistorHeat {
    pub branch: usize,
    pub voltage_drop_v: f64,
    pub current_a: f64,
    pub heat_j: f64,
}

/// Physical accounts from one native step. Polarization heat is already part
/// of `native.dissipated`, not an additional debit to circuit energy.
pub struct DebyeCircuitStep {
    /// Original accepted descriptor step and its complete work/loss account.
    pub native: StepRecord,
    /// The dielectric's instantaneous plus polarization storage endpoints [J].
    pub dielectric_energy_j: [f64; 2],
    /// Storage increment from charge changes, retaining work below the spacing
    /// of the endpoint energies under a large pre-existing bias [J].
    pub dielectric_energy_change_j: f64,
    /// Work supplied to this dielectric at its original terminal pair [J].
    pub dielectric_terminal_work_j: f64,
    /// Source-ordered physical heat contributions.
    pub polarization: Vec<PolarizationHeat>,
    /// Sum of those contributions, already included in native dissipation [J].
    pub polarization_heat_j: f64,
    /// Separate loss in the other original authored circuit resistors [J].
    pub other_resistor_heat_j: f64,
    /// Original branch identities for explicit thermal destinations.
    pub other_resistors: Vec<DebyeResistorHeat>,
    /// Dielectric storage change plus polarization heat minus terminal work [J].
    pub dielectric_energy_residual_j: f64,
    /// Whole-circuit storage change plus all heat minus external supply [J].
    pub circuit_energy_residual_j: f64,
}

/// A failed or cancelled proposal never changes the circuit or its memory.
pub struct DebyeCircuitCandidate {
    circuit: DebyeCircuit,
    step: DebyeCircuitStep,
}
impl DebyeCircuitCandidate {
    /// Proposed native charge/flux history with the same immutable source.
    pub const fn circuit(&self) -> &DebyeCircuit {
        &self.circuit
    }
    /// Original solve and independently recovered physical accounts.
    pub const fn step(&self) -> &DebyeCircuitStep {
        &self.step
    }
    /// Accept together with the physical receivers of the separately named heat.
    pub fn into_parts(self) -> (DebyeCircuit, DebyeCircuitStep) {
        (self.circuit, self.step)
    }
}

fn invalid(message: &'static str) -> MaterialCircuitError {
    MaterialCircuitError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, MaterialCircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("dielectric circuit arithmetic is unrepresentable"))
    }
}

fn accumulate(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), MaterialCircuitError> {
    let term = finite(term)?;
    let next = finite(*sum + term)?;
    *correction = finite(
        *correction
            + if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

/// Charge/flux increments times midpoint voltage/current give the quadratic
/// storage change without subtracting large endpoint energies. The correction
/// also retains the exact midpoint when adjacent represented states straddle it.
fn accumulate_storage_change(
    sum: &mut f64,
    correction: &mut f64,
    before: f64,
    after: f64,
    coefficient: f64,
) -> Result<(), MaterialCircuitError> {
    let change = finite(after - before)?;
    let midpoint = f64::midpoint(before, after);
    accumulate(sum, correction, change * (midpoint / coefficient))?;
    let offset = 0.5 * (before - midpoint) + 0.5 * (after - midpoint);
    accumulate(sum, correction, change * (offset / coefficient))
}

fn dielectric_storage_energy(
    source: &ResolvedDebyeCapacitor,
    instantaneous: Option<&MaterialDielectricCapacitor>,
    instantaneous_charge_c: f64,
    polarization_charges_c: &[f64],
) -> Result<f64, MaterialCircuitError> {
    match instantaneous {
        Some(law) => finite(
            law.stored_energy_j(instantaneous_charge_c)?
                + source.stored_energy_j(0.0, polarization_charges_c)?,
        ),
        None => Ok(source.stored_energy_j(instantaneous_charge_c, polarization_charges_c)?),
    }
}

impl DebyeCircuit {
    /// Replace ONE capacitor placeholder with its original dielectric model and
    /// append its internal nodes/branches. All original graph indices stay fixed.
    /// `charges_c` contains original capacitors in graph order; the selected
    /// capacitor's entry is its INSTANTANEOUS charge, not its total charge.
    /// Supply each polarization charge separately, in the source's pole order.
    ///
    /// At equilibrium under V, initialize q_inf=C_inf V and q_k=C_k V. A
    /// previously shorted capacitor may instead have q_inf=0 with nonzero q_k;
    /// the native circuit then predicts voltage recovery without a fitted gain.
    /// All capacitor charges and inductor flux linkages enter the original
    /// consistent-initialization solve unchanged. Incompatible ideal clamps
    /// refuse instead of creating an unaccounted charging impulse.
    ///
    /// The source and geometry remain frozen. This does not supply leakage,
    /// electrode ESR/ESL, breakdown, thermal feedback or parameter-work updates.
    /// Original other resistors are authored circuit elements; their aggregate
    /// loss stays separate from the source-bound dielectric's physical heat.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: CircuitGraph,
        branch: usize,
        source: ResolvedDebyeCapacitor,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        polarization_charges_c: Vec<f64>,
        initial_inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        Self::assemble_with_instantaneous(
            graph,
            branch,
            source,
            None,
            fluxes_wb,
            charges_c,
            polarization_charges_c,
            initial_inputs,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_with_instantaneous(
        mut graph: CircuitGraph,
        branch: usize,
        source: ResolvedDebyeCapacitor,
        reversible_instantaneous: Option<MaterialDielectricCapacitor>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        polarization_charges_c: Vec<f64>,
        initial_inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(&mut cancelled)?;
        let &(a, b, original) = graph
            .branches
            .get(branch)
            .ok_or_else(|| invalid("dielectric binding must name an original capacitor branch"))?;
        let poles = source.poles().len();
        if !matches!(original, Branch::Capacitor { .. })
            || a == b
            || a >= graph.node_count
            || b >= graph.node_count
            || polarization_charges_c.len() != poles
            || graph
                .node_count
                .checked_add(poles)
                .is_none_or(|n| n > super::MAX_MATERIAL_CIRCUIT_NODES)
            || graph
                .branches
                .len()
                .checked_add(2 * poles)
                .is_none_or(|n| n > super::MAX_MATERIAL_CIRCUIT_BRANCHES)
        {
            return Err(invalid(
                "bounded capacitor graph and complete original polarization history required",
            ));
        }
        let original_capacitors = graph
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Capacitor { .. }))
            .count();
        if charges_c.len() != original_capacitors {
            return Err(invalid(
                "one instantaneous charge per original capacitor required",
            ));
        }
        let capacitor_slot = graph.branches[..branch]
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Capacitor { .. }))
            .count();
        dielectric_storage_energy(
            &source,
            reversible_instantaneous.as_ref(),
            charges_c[capacitor_slot],
            &polarization_charges_c,
        )?;
        graph.branches[branch].2 = Branch::Capacitor {
            farads: source.instantaneous_capacitance_f(),
        };
        let mut loss_branches = Vec::with_capacity(poles);
        for pole in source.poles() {
            poll(&mut cancelled)?;
            let node = graph.node_count;
            graph.node_count += 1;
            loss_branches.push(graph.branches.len());
            graph.branches.push((
                a,
                node,
                Branch::Resistor {
                    ohms: pole.resistance_ohm,
                },
            ));
            graph.branches.push((
                node,
                b,
                Branch::Capacitor {
                    farads: pole.capacitance_f,
                },
            ));
        }
        admit_graph(&graph)?;
        let dae = if let Some(law) = &reversible_instantaneous {
            fs_phs::circuit::capacitance::assemble_circuit_with_capacitor_curves(
                &graph,
                &[fs_phs::circuit::capacitance::CapacitorCurve {
                    branch,
                    curve: law.curve().clone(),
                }],
                &[],
            )?
        } else {
            assemble_circuit(&graph)?
        };
        let mut all_charges = charges_c;
        all_charges.extend(polarization_charges_c);
        let initial = extended(&dae, &fluxes_wb, &all_charges, initial_inputs)?;
        let state = dae.consistent_initial_state(&initial, initial_inputs)?;
        let h = dae.system.hamiltonian(&state);
        if !h.is_finite() || h < 0.0 || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable initial dielectric circuit state"));
        }
        let instantaneous_charge_index = dae.charge_index[capacitor_slot];
        let polarization_charge_indices = dae.charge_index[original_capacitors..].to_vec();
        poll(&mut cancelled)?;
        Ok(Self {
            prepared: Arc::new(Prepared {
                source,
                graph,
                dae,
                branch,
                instantaneous_charge_index,
                polarization_charge_indices,
                loss_branches,
                material_resistors: Vec::new(),
                reversible_instantaneous,
            }),
            state,
        })
    }

    /// Original pole law, zero-field background and fixed dielectric geometry.
    /// With [`Self::reversible_instantaneous`], that source supplies the poles
    /// and matching low-field baseline; the retained nonlinear law supplies
    /// the actual instantaneous voltage/energy. Its finite-field response is
    /// not the linear source's frequency-domain admittance.
    pub fn source(&self) -> &ResolvedDebyeCapacitor {
        &self.prepared.source
    }
    /// Expanded native graph, preserving original branch indices.
    pub fn graph(&self) -> &CircuitGraph {
        &self.prepared.graph
    }
    /// Prepared original Kirchhoff descriptor and its coordinate maps.
    pub fn dae(&self) -> &CircuitDae {
        &self.prepared.dae
    }
    /// Accepted native differential and algebraic coordinates.
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    /// Instantaneous charge q_infinity alone [C]; polarization stays separate.
    pub fn instantaneous_charge_c(&self) -> f64 {
        self.state[self.prepared.instantaneous_charge_index]
    }
    /// Retained polarization charges in original source order [C].
    pub fn polarization_charges_c(&self) -> Vec<f64> {
        self.prepared
            .polarization_charge_indices
            .iter()
            .map(|&i| self.state[i])
            .collect()
    }
    /// Physical original terminal voltage from the instantaneous storage [V].
    pub fn terminal_voltage_v(&self) -> f64 {
        match self.reversible_instantaneous() {
            Some(law) => law
                .voltage_v(self.instantaneous_charge_c())
                .unwrap_or(f64::NAN),
            None => self.instantaneous_charge_c() / self.source().instantaneous_capacitance_f(),
        }
    }
    /// Instantaneous and all polarization storage, counted exactly once [J].
    pub fn dielectric_energy_j(&self) -> Result<f64, MaterialCircuitError> {
        self.energy(&self.state)
    }
    fn energy(&self, state: &[f64]) -> Result<f64, MaterialCircuitError> {
        let mut charges = [0.0; MAX_DIELECTRIC_POLES];
        for (charge, &index) in charges
            .iter_mut()
            .zip(&self.prepared.polarization_charge_indices)
        {
            *charge = state[index];
        }
        dielectric_storage_energy(
            self.source(),
            self.reversible_instantaneous(),
            state[self.prepared.instantaneous_charge_index],
            &charges[..self.source().poles().len()],
        )
    }

    fn dielectric_energy_change(
        &self,
        before: &[f64],
        after: &[f64],
    ) -> Result<f64, MaterialCircuitError> {
        let (mut sum, mut correction) = (0.0, 0.0);
        let slot = self.prepared.instantaneous_charge_index;
        if let Some(law) = self.reversible_instantaneous() {
            accumulate(
                &mut sum,
                &mut correction,
                law.curve().energy_change_j(before[slot], after[slot])?,
            )?;
        } else {
            accumulate_storage_change(
                &mut sum,
                &mut correction,
                before[slot],
                after[slot],
                self.source().instantaneous_capacitance_f(),
            )?;
        }
        for (&slot, pole) in self
            .prepared
            .polarization_charge_indices
            .iter()
            .zip(self.source().poles())
        {
            accumulate_storage_change(
                &mut sum,
                &mut correction,
                before[slot],
                after[slot],
                pole.capacitance_f,
            )?;
        }
        finite(sum + correction)
    }

    fn circuit_energy_change(
        &self,
        before: &[f64],
        after: &[f64],
    ) -> Result<f64, MaterialCircuitError> {
        let (mut sum, mut correction) = (0.0, 0.0);
        let (mut capacitors, mut inductors) = (0, 0);
        for &(_, _, branch) in &self.graph().branches {
            let (slot, coefficient) = match branch {
                Branch::Capacitor { farads } => {
                    let slot = self.dae().charge_index[capacitors];
                    capacitors += 1;
                    (slot, farads)
                }
                Branch::Inductor { henries } => {
                    let slot = self.dae().flux_index[inductors];
                    inductors += 1;
                    (slot, henries)
                }
                _ => continue,
            };
            if let Some(law) = self
                .reversible_instantaneous()
                .filter(|_| slot == self.prepared.instantaneous_charge_index)
            {
                accumulate(
                    &mut sum,
                    &mut correction,
                    law.curve().energy_change_j(before[slot], after[slot])?,
                )?;
            } else {
                accumulate_storage_change(
                    &mut sum,
                    &mut correction,
                    before[slot],
                    after[slot],
                    coefficient,
                )?;
            }
        }
        finite(sum + correction)
    }

    /// Advance once through `CircuitDae::step_audited`. Coefficients stay prepared
    /// and charge history is never reinitialized. Held port changes must be
    /// compatible with the existing charges/fluxes; ideal impulses are refused.
    /// Native midpoint voltage drops recover every polarization resistor's heat
    /// independently of the storage residual. Its sum and dielectric terminal
    /// work must close, as must the complete circuit's source/storage/loss balance.
    /// The supplied tolerances do not certify time accuracy: resolve the shortest
    /// relaxation time and check step refinement for physical transients.
    pub fn step_candidate(
        &self,
        inputs: &[f64],
        duration_s: f64,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeCircuitCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.power_tolerance_w)
        {
            return Err(invalid(
                "positive dielectric interval and physical balance tolerances required",
            ));
        }
        let dae = self.dae();
        let x0 = dae.consistent_initial_state(&self.state, inputs)?;
        poll(&mut cancelled)?;
        let (native, _) = dae.step_audited(&x0, inputs, duration_s)?;
        if native.x.iter().chain(&native.y).any(|v| !v.is_finite())
            || [
                native.delta_h,
                native.dissipated,
                native.supplied,
                native.solver_residual,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || native.dissipated < 0.0
        {
            return Err(invalid("invalid native dielectric circuit step"));
        }
        if self.reversible_instantaneous().is_some() {
            // The same frozen source must admit the physical endpoint without
            // projecting a charge/flux or hiding an ideal-clamp impulse.
            dae.consistent_initial_state(&native.x, inputs)?;
        }
        let voltage = |node| {
            if node == 0 {
                0.0
            } else {
                let slot = dae.node_potential_index[node - 1];
                f64::midpoint(x0[slot], native.x[slot])
            }
        };
        let mut polarization = Vec::with_capacity(self.source().poles().len());
        let mut polarization_heat_j = 0.0;
        let mut other_resistor_heat_j = 0.0;
        let mut other_resistors = Vec::new();
        for (branch, &(a, b, element)) in self.graph().branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let Branch::Resistor { ohms } = element else {
                continue;
            };
            let drop = finite(voltage(a) - voltage(b))?;
            let current = finite(drop / ohms)?;
            let power = finite(drop * current)?;
            let heat = finite(duration_s * power)?;
            if power < 0.0
                || (drop != 0.0 && (current == 0.0 || power == 0.0))
                || (power != 0.0 && heat == 0.0)
            {
                return Err(invalid(
                    "nonzero dielectric or circuit resistor heat disappeared",
                ));
            }
            if let Some(pole) = self
                .prepared
                .loss_branches
                .iter()
                .position(|&index| index == branch)
            {
                let charge = self.prepared.polarization_charge_indices[pole];
                let charge_current = finite((native.x[charge] - x0[charge]) / duration_s)?;
                if finite(drop * charge_current - power)?.abs() > controls.power_tolerance_w {
                    return Err(MaterialCircuitError::Balance(
                        "polarization charge work and resistor heat disagree",
                    ));
                }
                polarization_heat_j = finite(polarization_heat_j + heat)?;
                polarization.push(PolarizationHeat {
                    pole,
                    resistor_branch: branch,
                    voltage_drop_v: drop,
                    current_a: current,
                    heat_j: heat,
                });
            } else {
                other_resistor_heat_j = finite(other_resistor_heat_j + heat)?;
                other_resistors.push(DebyeResistorHeat {
                    branch,
                    voltage_drop_v: drop,
                    current_a: current,
                    heat_j: heat,
                });
            }
        }
        let (a, b, _) = self.graph().branches[self.prepared.branch];
        let (mut charge_change, mut charge_correction) = (0.0, 0.0);
        for slot in std::iter::once(self.prepared.instantaneous_charge_index)
            .chain(self.prepared.polarization_charge_indices.iter().copied())
        {
            accumulate(
                &mut charge_change,
                &mut charge_correction,
                native.x[slot] - x0[slot],
            )?;
        }
        let charge_change = finite(charge_change + charge_correction)?;
        let dielectric_terminal_work_j = finite((voltage(a) - voltage(b)) * charge_change)?;
        let dielectric_energy_j = [self.energy(&x0)?, self.energy(&native.x)?];
        let dielectric_energy_change_j = self.dielectric_energy_change(&x0, &native.x)?;
        let dielectric_energy_residual_j =
            finite(dielectric_energy_change_j + polarization_heat_j - dielectric_terminal_work_j)?;
        let h0 = finite(dae.system.hamiltonian(&x0))?;
        let h1 = finite(dae.system.hamiltonian(&native.x))?;
        let circuit_energy_change_j = self.circuit_energy_change(&x0, &native.x)?;
        let heat = finite(polarization_heat_j + other_resistor_heat_j)?;
        let circuit_energy_residual_j = finite(native.delta_h + heat - native.supplied)?;
        if h0 < 0.0
            || h1 < 0.0
            || finite(heat - native.dissipated)?.abs() > controls.energy_tolerance_j
            || finite(circuit_energy_change_j - native.delta_h)?.abs() > controls.energy_tolerance_j
            || dielectric_energy_residual_j.abs() > controls.energy_tolerance_j
            || circuit_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(MaterialCircuitError::Balance(
                "dielectric, other resistor and native circuit energy do not close",
            ));
        }
        poll(&mut cancelled)?;
        let circuit = Self {
            prepared: self.prepared.clone(),
            state: native.x.clone(),
        };
        Ok(DebyeCircuitCandidate {
            circuit,
            step: DebyeCircuitStep {
                native,
                dielectric_energy_j,
                dielectric_energy_change_j,
                dielectric_terminal_work_j,
                polarization,
                polarization_heat_j,
                other_resistor_heat_j,
                other_resistors,
                dielectric_energy_residual_j,
                circuit_energy_residual_j,
            },
        })
    }
}
