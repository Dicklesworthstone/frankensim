//! Fixed-clock exact flow of the original material-bound Kirchhoff circuit.
//!
//! Preparation eliminates only algebraic multipliers, by one native LU of the
//! original index-one block. Every flux, charge and mechanical coordinate is
//! retained. The native quadratic-flow owner advances those coordinates; native
//! source-square observers independently integrate each original resistor and
//! the signed physical source outputs, including direct input feedthrough.
//! No heat is inferred from an energy deficit and no endpoint is rescaled.
//!
//! Material, topology, resistor state and duration are held after preparation.
//! A singular algebraic block refuses explicitly; the ordinary descriptor API
//! remains available for compatible higher-index circuits. Exact means the
//! linear held-input model up to floating-point error, not an exact continuum
//! rod, an outward-rounded domain certificate or a real-time deadline.

use super::{
    MaterialCircuitError, MaterialCircuitSolve, PiezoelectricCircuit, PiezoelectricCircuitStorage,
    bind, invalid, native, poll, positive,
};
use fs_la::factor::lu;
use fs_material::conductor::OhmicDrive;
use fs_math::{
    det,
    eft::{two_prod, two_sum},
};
use fs_phs::circuit::CircuitDae;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedSourceLoss, split_quadratic,
};
use fs_phs::{PortHamiltonian, PreparedStepError, ResistanceWitness, Storage};

type Result<T> = core::result::Result<T, MaterialCircuitError>;

/// One original resistor's observations over the accepted continuous interval.
/// Mean power includes current variance: it is not mean voltage times mean
/// current. The complete source binding is retained by model().circuit().
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactPiezoelectricResistorHeat {
    /// Index in the original CircuitGraph branches.
    pub branch: usize,
    pub mean_voltage_v: f64,
    pub mean_current_a: f64,
    pub rms_current_a: f64,
    pub mean_power_w: f64,
    pub heat_j: f64,
}

/// Independent physical accounting for one accepted fixed-clock interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactPiezoelectricCircuitStep {
    pub duration_s: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub piezoelectric_energy_j: [f64; 2],
    pub other_electrical_energy_j: [f64; 2],
    pub electrical_supply_j: f64,
    pub mechanical_supply_j: f64,
    pub heat_j: f64,
    pub energy_residual_j: f64,
    /// Conservative computed bound on the complete stored energy during the
    /// interval. It can refuse a physically valid window with unrelated large
    /// electrical storage; it is not a directed-rounding certificate.
    pub interval_energy_bound_j: f64,
}

/// Reusable actual circuit, with immutable material/topology and one clock.
/// Warm steps allocate no owner buffers, rebuild no matrices and perform no
/// factorization, exponential or Newton solve. Original constitutive admission,
/// physical source work and resistor heat remain authoritative. On any refusal
/// or cancellation, the model and all previously published reports are intact.
pub struct PreparedPiezoelectricCircuit {
    model: PiezoelectricCircuit,
    descriptor: CircuitDae,
    storage: PiezoelectricCircuitStorage,
    flow: PreparedQuadraticFlow,
    losses: Vec<PreparedSourceLoss>,
    // Each group observes at most d signed original electrical outputs. Its
    // squared integral is NEVER used as physical dissipation.
    source_outputs: Vec<PreparedSourceLoss>,
    algebraic_effort: Vec<f64>,
    algebraic_input: Vec<f64>,
    root_conductances: Vec<f64>,
    duration_s: f64,
    state: Vec<f64>,
    full_state: Vec<f64>,
    inputs: Vec<f64>,
    endpoint: Vec<f64>,
    full_candidate: Vec<f64>,
    effort: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    observation_scratch: Vec<f64>,
    output_candidate: Vec<f64>,
    mean_outputs: Vec<f64>,
    heat_candidate: Vec<ExactPiezoelectricResistorHeat>,
    resistor_heat: Vec<ExactPiezoelectricResistorHeat>,
    last_step: Option<ExactPiezoelectricCircuitStep>,
}

impl PreparedPiezoelectricCircuit {
    /// Prepare the original descriptor at the specified initial electrical
    /// inputs. Force does not constrain the initial conserved coordinates.
    /// Dense byte/work allowances are divided across the propagation and
    /// observation caches, so each cache cannot consume the entire allowance.
    /// Existing descriptor/state limits also bound the cold Kirchhoff solve.
    pub fn new(
        model: PiezoelectricCircuit,
        initial_electrical_ports: &[f64],
        duration_s: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !positive(duration_s) {
            return Err(invalid(
                "prepared piezoelectric circuit requires a positive finite duration",
            ));
        }
        let (_, _, original) = bind(&model.circuit.graph, model.circuit.resistors.clone(), &[])?;
        let (descriptor, storage) = model.assemble(&original)?;
        let initial = model.extended(&original, &descriptor, initial_electrical_ports)?;
        let (n, d, m) = (
            descriptor.system.state_dim(),
            descriptor.system.differential_dim(),
            descriptor.ports,
        );
        let mut inputs = initial_electrical_ports.to_vec();
        inputs.push(0.0);
        model.admit(&descriptor, &initial, &inputs)?;
        poll(&mut cancelled)?;
        let (algebraic_effort, algebraic_input) = eliminate(&descriptor, &mut cancelled)?;
        let (j, r, input_map, loss_rows, direct_loss, root_conductances) = lower(
            &model,
            &descriptor,
            &algebraic_effort,
            &algebraic_input,
            &mut cancelled,
        )?;
        let witness = ResistanceWitness::try_from_source_rows(d, r.clone(), loss_rows.clone())
            .map_err(native)?;
        let system = PortHamiltonian::new(d, m, j, r, input_map, Box::new(storage.clone()))
            .map_err(native)?
            .with_resistance_witness(witness)
            .map_err(native)?;
        let split = split_quadratic(&system, split_budget).map_err(native)?;
        let groups = (m - 1).div_ceil(d);
        let count = 1 + model.circuit.resistors.len() + groups;
        let cache_budget = DenseExponentialBudget {
            max_dimension: dense_budget.max_dimension,
            max_dense_bytes: dense_budget.max_dense_bytes / count,
            max_multiply_terms: dense_budget.max_multiply_terms / count as u64,
        };
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &[],
            &[],
            cache_budget,
            &mut cancelled,
        )
        .map_err(prepared_error)?;
        // Use the identical operation order as the native flow's generator.
        let (j, r, input_map) = system.structure();
        let q = &storage.quadratic;
        let mut generator = vec![0.0; d * d];
        for i in 0..d {
            poll(&mut cancelled)?;
            for k in 0..d {
                for column in 0..d {
                    generator[i * d + column] += (j[i * d + k] - r[i * d + k]) * q[k * d + column];
                }
            }
        }
        let loss_state = multiply(&loss_rows, q, model.circuit.resistors.len(), d, d);
        let mut losses = Vec::with_capacity(model.circuit.resistors.len());
        for index in 0..model.circuit.resistors.len() {
            losses.push(
                PreparedSourceLoss::new_with_feedthrough_controlled(
                    d,
                    m,
                    &generator,
                    input_map,
                    &loss_state[index * d..(index + 1) * d],
                    &direct_loss[index * m..(index + 1) * m],
                    duration_s,
                    cache_budget,
                    &mut cancelled,
                )
                .map_err(prepared_error)?,
            );
        }
        // y = G_d^T e + G_a^T (L e + N u), in the original port order.
        // Integrating S*x + D*u together retains a small physical source
        // current between large cancelling node/source terms.
        let (_, _, g) = descriptor.system.structure();
        let mut output_effort = vec![0.0; (m - 1) * d];
        let mut output_input = vec![0.0; (m - 1) * m];
        for port in 0..m - 1 {
            poll(&mut cancelled)?;
            for column in 0..d {
                output_effort[port * d + column] = g[column * m + port]
                    + (d..n)
                        .map(|i| g[i * m + port] * algebraic_effort[(i - d) * d + column])
                        .sum::<f64>();
            }
            for source in 0..m {
                output_input[port * m + source] = (d..n)
                    .map(|i| g[i * m + port] * algebraic_input[(i - d) * m + source])
                    .sum();
            }
        }
        let output_state = multiply(&output_effort, q, m - 1, d, d);
        let mut source_outputs = Vec::with_capacity(groups);
        for first in (0..m - 1).step_by(d) {
            let last = (first + d).min(m - 1);
            source_outputs.push(
                PreparedSourceLoss::new_with_feedthrough_controlled(
                    d,
                    m,
                    &generator,
                    input_map,
                    &output_state[first * d..last * d],
                    &output_input[first * m..last * m],
                    duration_s,
                    cache_budget,
                    &mut cancelled,
                )
                .map_err(prepared_error)?,
            );
        }
        let state = initial[..d].to_vec();
        let mut full_state = vec![0.0; n];
        let mut effort = vec![0.0; n];
        reconstruct(
            &descriptor,
            &storage,
            &algebraic_effort,
            &algebraic_input,
            &state,
            &inputs,
            &mut full_state,
            &mut effort,
            &mut cancelled,
        )?;
        let resistor_heat: Vec<_> = model
            .circuit
            .resistors
            .iter()
            .map(|source| ExactPiezoelectricResistorHeat {
                branch: source.branch,
                mean_voltage_v: 0.0,
                mean_current_a: 0.0,
                rms_current_a: 0.0,
                mean_power_w: 0.0,
                heat_j: 0.0,
            })
            .collect();
        poll(&mut cancelled)?;
        Ok(Self {
            model,
            descriptor,
            storage,
            flow,
            losses,
            source_outputs,
            algebraic_effort,
            algebraic_input,
            root_conductances,
            duration_s,
            state,
            full_state,
            inputs,
            endpoint: vec![0.0; d],
            full_candidate: vec![0.0; n],
            effort,
            augmented: vec![0.0; d + m + 1],
            balanced: vec![0.0; d],
            observation_scratch: vec![0.0; 2 * d],
            output_candidate: vec![0.0; m - 1],
            mean_outputs: vec![0.0; m - 1],
            heat_candidate: resistor_heat.clone(),
            resistor_heat,
            last_step: None,
        })
    }

    pub const fn model(&self) -> &PiezoelectricCircuit {
        &self.model
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Return every conserved coordinate and its original source bindings.
    pub fn into_model(self) -> PiezoelectricCircuit {
        self.model
    }
    /// Complete original descriptor layout, including Kirchhoff multipliers.
    /// Before the first step it uses the preparation's initial electrical inputs.
    pub fn descriptor_state(&self) -> &[f64] {
        &self.full_state
    }
    pub const fn last_step(&self) -> Option<&ExactPiezoelectricCircuitStep> {
        self.last_step.as_ref()
    }
    /// Last accepted interval's heat, in the original stable resistor order.
    pub fn resistor_heat(&self) -> Option<&[ExactPiezoelectricResistorHeat]> {
        self.last_step
            .as_ref()
            .map(|_| self.resistor_heat.as_slice())
    }
    /// Actual trajectory-mean electrical port outputs, in original port order.
    pub fn mean_electrical_outputs(&self) -> Option<&[f64]> {
        self.last_step
            .as_ref()
            .map(|_| self.mean_outputs.as_slice())
    }

    /// Advance one held source/force interval and publish only after all checks.
    /// The complete-energy envelope checks strain and field over the window,
    /// including excursions that endpoint-only checks could miss.
    pub fn step(
        &mut self,
        electrical_ports: &[f64],
        force_n: f64,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactPiezoelectricCircuitStep> {
        poll(&mut cancelled)?;
        let (d, m) = (self.state.len(), self.inputs.len());
        if electrical_ports.len() != m - 1
            || electrical_ports.iter().any(|v| !v.is_finite())
            || !force_n.is_finite()
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.power_tolerance_w)
        {
            return Err(invalid(
                "prepared circuit needs complete finite sources and positive balance tolerances",
            ));
        }
        self.inputs[..m - 1].copy_from_slice(electrical_ports);
        self.inputs[m - 1] = force_n;
        let norm = self
            .flow
            .interval_coordinate_norm_bound(&self.state, &self.inputs)
            .map_err(native)?;
        let rod = &self.model.rod;
        let limits = rod.limits();
        let extension_bound = norm / det::sqrt(rod.short_circuit_stiffness_n_m());
        let voltage_bound = norm / det::sqrt(rod.axially_clamped_capacitance_f());
        let interval_energy_bound_j = 0.5 * norm * norm;
        if !extension_bound.is_finite()
            || !voltage_bound.is_finite()
            || !interval_energy_bound_j.is_finite()
            || extension_bound / rod.length_m() > limits.max_abs_strain
            || voltage_bound / rod.length_m() > limits.max_abs_electric_field_v_m
        {
            return Err(invalid(
                "complete prepared circuit window exceeds the piezoelectric small-signal envelope",
            ));
        }
        self.flow
            .evaluate_into(
                &self.state,
                &self.inputs,
                &mut self.augmented,
                &mut self.balanced,
                &mut self.endpoint,
                &mut [],
            )
            .map_err(native)?;
        poll(&mut cancelled)?;
        reconstruct(
            &self.descriptor,
            &self.storage,
            &self.algebraic_effort,
            &self.algebraic_input,
            &self.endpoint,
            &self.inputs,
            &mut self.full_candidate,
            &mut self.effort,
            &mut cancelled,
        )?;
        let piezoelectric_energy_j = [
            self.storage.component_energy(&self.state),
            self.storage.component_energy(&self.endpoint),
        ];
        let other_electrical_energy_j = [
            self.storage.other.hamiltonian(&self.state),
            self.storage.other.hamiltonian(&self.endpoint),
        ];
        let stored_before_j = piezoelectric_energy_j[0] + other_electrical_energy_j[0];
        let stored_after_j = piezoelectric_energy_j[1] + other_electrical_energy_j[1];
        let mechanical_supply_j =
            force_n * (self.endpoint[self.storage.extension] - self.state[self.storage.extension]);
        let (mut electrical_supply_j, mut work_tail) = (0.0, 0.0);
        for port in 0..m - 1 {
            let integral = self.source_outputs[port / d]
                .signed_integral_into_controlled(
                    port % d,
                    &self.state,
                    &self.inputs,
                    &mut self.observation_scratch,
                    &mut cancelled,
                )
                .map_err(prepared_error)?;
            self.output_candidate[port] = integral / self.duration_s;
            add_product(
                &mut electrical_supply_j,
                &mut work_tail,
                electrical_ports[port],
                integral,
            );
        }
        electrical_supply_j += work_tail;
        let mut heat_j = 0.0;
        for (index, loss) in self.losses.iter().enumerate() {
            let source = &self.model.circuit.resistors[index];
            let joules = loss
                .evaluate_into_controlled(
                    &self.state,
                    &self.inputs,
                    &mut self.observation_scratch,
                    &mut cancelled,
                )
                .map_err(prepared_error)?;
            let mean_voltage_v = loss
                .signed_integral_into_controlled(
                    0,
                    &self.state,
                    &self.inputs,
                    &mut self.observation_scratch,
                    &mut cancelled,
                )
                .map_err(prepared_error)?
                / self.root_conductances[index]
                / self.duration_s;
            let mean_current_a = source
                .conductor
                .current_a_for_drive(OhmicDrive::Voltage(mean_voltage_v))?;
            let mean_power_w = joules / self.duration_s;
            let rms_current_a = det::sqrt(mean_power_w / source.conductor.resistance_ohm());
            let material_power = source.conductor.joule_power_w(rms_current_a)?;
            let mean_current_power = source.conductor.joule_power_w(mean_current_a)?;
            if [
                joules,
                mean_voltage_v,
                mean_current_a,
                mean_power_w,
                rms_current_a,
                material_power,
                mean_current_power,
            ]
            .iter()
            .any(|v| !v.is_finite())
                || joules < 0.0
                || joules != 0.0 && mean_power_w == 0.0
                || mean_power_w != 0.0 && rms_current_a == 0.0
                || mean_current_power != 0.0 && joules == 0.0
                || (material_power - mean_power_w).abs() > controls.power_tolerance_w
                || mean_current_power - mean_power_w > controls.power_tolerance_w
            {
                return Err(MaterialCircuitError::Balance(
                    "prepared resistor integral and original material heat disagree",
                ));
            }
            heat_j += joules;
            self.heat_candidate[index] = ExactPiezoelectricResistorHeat {
                branch: source.branch,
                mean_voltage_v,
                mean_current_a,
                rms_current_a,
                mean_power_w,
                heat_j: joules,
            };
        }
        let energy_residual_j =
            stored_after_j - stored_before_j + heat_j - electrical_supply_j - mechanical_supply_j;
        if [
            stored_before_j,
            stored_after_j,
            mechanical_supply_j,
            electrical_supply_j,
            heat_j,
            energy_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || self.output_candidate.iter().any(|v| !v.is_finite())
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(MaterialCircuitError::Balance(
                "prepared piezoelectric circuit work and physical heat do not close",
            ));
        }
        let record = ExactPiezoelectricCircuitStep {
            duration_s: self.duration_s,
            stored_before_j,
            stored_after_j,
            piezoelectric_energy_j,
            other_electrical_energy_j,
            electrical_supply_j,
            mechanical_supply_j,
            heat_j,
            energy_residual_j,
            interval_energy_bound_j,
        };
        poll(&mut cancelled)?;
        // No fallible operation or callback follows publication.
        for (value, &index) in self
            .model
            .circuit
            .fluxes_wb
            .iter_mut()
            .zip(&self.descriptor.flux_index)
        {
            *value = self.endpoint[index];
        }
        for (value, &index) in self
            .model
            .circuit
            .charges_c
            .iter_mut()
            .zip(&self.descriptor.charge_index)
        {
            *value = self.endpoint[index];
        }
        self.model.extension_m = self.endpoint[self.storage.extension];
        self.model.momentum_kg_m_s = self.endpoint[self.storage.momentum];
        self.state.copy_from_slice(&self.endpoint);
        self.full_state.copy_from_slice(&self.full_candidate);
        self.resistor_heat.copy_from_slice(&self.heat_candidate);
        self.mean_outputs.copy_from_slice(&self.output_candidate);
        self.last_step = Some(record);
        Ok(record)
    }
}

fn prepared_error(error: PreparedStepError) -> MaterialCircuitError {
    match error {
        PreparedStepError::Cancelled => MaterialCircuitError::Cancelled,
        PreparedStepError::Solver(error) => native(error),
    }
}

// E xdot = M e + G u, M = J-R, e = [grad H(z), lambda].
// The original algebraic equations give lambda = L grad H + N u.
// LU acts only during preparation, without choosing a nullspace or discarding
// any conserved coordinate. Its residual gate is componentwise and unit-free.
fn eliminate(
    dae: &CircuitDae,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let (n, d, m) = (
        dae.system.state_dim(),
        dae.system.differential_dim(),
        dae.ports,
    );
    let a = n - d;
    let (j, r, g) = dae.system.structure();
    let mut block = vec![0.0; a * a];
    for i in 0..a {
        poll(cancelled)?;
        for k in 0..a {
            block[i * a + k] = j[(d + i) * n + d + k] - r[(d + i) * n + d + k];
        }
    }
    let factor = lu(&block, a).map_err(|_| {
        invalid("prepared circuit requires an invertible index-one algebraic block")
    })?;
    let (mut l, mut input) = (vec![0.0; a * d], vec![0.0; a * m]);
    let mut rhs = vec![0.0; a];
    let mut solution = vec![0.0; a];
    for column in 0..d + m {
        poll(cancelled)?;
        for i in 0..a {
            rhs[i] = if column < d {
                -(j[(d + i) * n + column] - r[(d + i) * n + column])
            } else {
                -g[(d + i) * m + column - d]
            };
        }
        solution.copy_from_slice(&rhs);
        factor.solve(&mut solution);
        for i in 0..a {
            let mut residual = -rhs[i];
            let mut scale = rhs[i].abs();
            for k in 0..a {
                let term = block[i * a + k] * solution[k];
                residual += term;
                scale += term.abs();
            }
            if !roundoff(residual, scale, a + 1) {
                return Err(invalid(
                    "prepared Kirchhoff elimination failed its original algebraic rows",
                ));
            }
            if column < d {
                l[i * d + column] = solution[i];
            } else {
                input[i * m + column - d] = solution[i];
            }
        }
    }
    Ok((l, input))
}

type Lowered = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

// The reduced R retains the original physical resistor squares. The native
// ResistanceWitness and PortHamiltonian admissions check its PSD factor.
// The skew part must reconstruct the unchanged Schur operator within its
// componentwise arithmetic scale; no eigenvalue clipping repairs a bad map.
fn lower(
    model: &PiezoelectricCircuit,
    dae: &CircuitDae,
    l: &[f64],
    input: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Lowered> {
    let (n, d, m) = (
        dae.system.state_dim(),
        dae.system.differential_dim(),
        dae.ports,
    );
    let (j, r, g) = dae.system.structure();
    let mut k = vec![0.0; d * d];
    let mut scales = vec![0.0; d * d];
    let mut b = vec![0.0; d * m];
    for i in 0..d {
        poll(cancelled)?;
        for column in 0..d {
            k[i * d + column] = j[i * n + column] - r[i * n + column];
            scales[i * d + column] = k[i * d + column].abs();
            for a in d..n {
                let term = (j[i * n + a] - r[i * n + a]) * l[(a - d) * d + column];
                k[i * d + column] += term;
                scales[i * d + column] += term.abs();
            }
        }
        for port in 0..m {
            b[i * m + port] = g[i * m + port];
            for a in d..n {
                b[i * m + port] += (j[i * n + a] - r[i * n + a]) * input[(a - d) * m + port];
            }
        }
    }
    let mut rows = vec![0.0; model.circuit.resistors.len() * d];
    let mut direct = vec![0.0; model.circuit.resistors.len() * m];
    let mut roots = Vec::with_capacity(model.circuit.resistors.len());
    for (index, source) in model.circuit.resistors.iter().enumerate() {
        poll(cancelled)?;
        let root = det::sqrt(1.0 / source.conductor.resistance_ohm());
        if !positive(root) {
            return Err(invalid("prepared resistor source scale is unrepresentable"));
        }
        roots.push(root);
        let (a, b, _) = model.circuit.graph.branches[source.branch];
        for (node, sign) in [(a, 1.0), (b, -1.0)] {
            if node == 0 {
                continue;
            }
            let row = dae.node_potential_index[node - 1] - d;
            for column in 0..d {
                rows[index * d + column] += sign * root * l[row * d + column];
            }
            for port in 0..m {
                direct[index * m + port] += sign * root * input[row * m + port];
            }
        }
    }
    let mut reduced_j = vec![0.0; d * d];
    let mut reduced_r = vec![0.0; d * d];
    for i in 0..d {
        poll(cancelled)?;
        for column in i..d {
            let value = (0..roots.len())
                .map(|row| rows[row * d + i] * rows[row * d + column])
                .sum::<f64>();
            reduced_r[i * d + column] = value;
            reduced_r[column * d + i] = value;
            if i != column {
                let skew = 0.5 * k[i * d + column] - 0.5 * k[column * d + i];
                reduced_j[i * d + column] = skew;
                reduced_j[column * d + i] = -skew;
            }
        }
    }
    for i in 0..d {
        for column in 0..d {
            let scale =
                scales[i * d + column] + scales[column * d + i] + reduced_r[i * d + column].abs();
            if !roundoff(
                reduced_j[i * d + column] - reduced_r[i * d + column] - k[i * d + column],
                scale,
                n + 1,
            ) {
                return Err(invalid(
                    "physical resistor squares do not reconstruct the original Kirchhoff flow",
                ));
            }
        }
    }
    if b.iter().chain(&rows).chain(&direct).any(|v| !v.is_finite()) {
        return Err(invalid(
            "prepared circuit input or resistor maps are unrepresentable",
        ));
    }
    Ok((reduced_j, reduced_r, b, rows, direct, roots))
}

fn multiply(a: &[f64], b: &[f64], rows: usize, inner: usize, columns: usize) -> Vec<f64> {
    let mut result = vec![0.0; rows * columns];
    for i in 0..rows {
        for k in 0..inner {
            for j in 0..columns {
                result[i * columns + j] += a[i * inner + k] * b[k * columns + j];
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn reconstruct(
    dae: &CircuitDae,
    storage: &PiezoelectricCircuitStorage,
    l: &[f64],
    input: &[f64],
    state: &[f64],
    sources: &[f64],
    full: &mut [f64],
    effort: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let (n, d, m) = (
        dae.system.state_dim(),
        dae.system.differential_dim(),
        dae.ports,
    );
    full[..d].copy_from_slice(state);
    storage.gradient(state, &mut effort[..d]);
    if !storage.hamiltonian(state).is_finite() || effort[..d].iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "prepared circuit left the original piezoelectric material domain",
        ));
    }
    for a in 0..n - d {
        poll(cancelled)?;
        let (mut value, mut tail) = (0.0, 0.0);
        for i in 0..d {
            add_product(&mut value, &mut tail, l[a * d + i], effort[i]);
        }
        for port in 0..m {
            add_product(&mut value, &mut tail, input[a * m + port], sources[port]);
        }
        full[d + a] = value + tail;
    }
    effort[d..].copy_from_slice(&full[d..]);
    let (j, r, g) = dae.system.structure();
    for row in d..n {
        poll(cancelled)?;
        let (mut residual, mut tail, mut scale) = (0.0, 0.0, 0.0);
        for column in 0..n {
            let coefficient = j[row * n + column] - r[row * n + column];
            add_product(&mut residual, &mut tail, coefficient, effort[column]);
            scale += (coefficient * effort[column]).abs();
        }
        for port in 0..m {
            add_product(&mut residual, &mut tail, g[row * m + port], sources[port]);
            scale += (g[row * m + port] * sources[port]).abs();
        }
        if !roundoff(residual + tail, scale, n + m) {
            return Err(MaterialCircuitError::Balance(
                "prepared endpoint does not satisfy the original Kirchhoff rows",
            ));
        }
    }
    if full.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "prepared circuit algebraic state is unrepresentable",
        ));
    }
    Ok(())
}

fn roundoff(residual: f64, scale: f64, terms: usize) -> bool {
    residual.is_finite()
        && scale.is_finite()
        && residual.abs() <= 256.0 * f64::EPSILON * terms as f64 * scale
}

fn add_product(sum: &mut f64, tail: &mut f64, a: f64, b: f64) {
    let (product, product_error) = two_prod(a, b);
    let (next, error) = two_sum(*sum, product);
    *sum = next;
    *tail += product_error + error;
}
