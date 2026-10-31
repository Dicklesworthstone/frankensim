//! Original circuit constraints closed against the resolved winding impedance.
//!
//! A positive reference resistor and a compensating current port replace each
//! field-owned winding only while the native descriptor is factored. If H maps
//! those injected currents to terminal voltages and G is the reference diagonal
//! conductance, the physical currents solve `[(I-HG) Z + H] i = v_open`.
//! Replaying `j = G Z i - i` cancels the reference termination exactly. Its heat
//! is excluded from every physical result. No circuit incidence, MNA solver,
//! electromagnetic field operator or empirical loss law is reimplemented here.

use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;
use fs_phs::{
    circuit::{
        Branch, CircuitError, CircuitGraph, MutualInductance,
        assemble_circuit_with_mutual_inductance,
    },
    harmonic::{HarmonicBudget, HarmonicError, HarmonicModel, MAX_HARMONIC_DIMENSION},
};

use super::{HarmonicWindingError, HarmonicWindingFrequency, HarmonicWindingResponse};

/// Bound on the caller's graph before any native dense circuit allocation.
pub const MAX_HARMONIC_CIRCUIT_BRANCHES: usize = 512;
/// Original source control indices and physical source branches are both bounded.
pub const MAX_HARMONIC_CIRCUIT_SOURCES: usize = 128;
/// Field terminal closure remains a small dense system.
pub const MAX_HARMONIC_CIRCUIT_WINDINGS: usize = 16;

/// Native descriptor budgets and independently admitted terminal/work defects.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicCircuitControls {
    pub native: HarmonicBudget,
    /// Additional network, interface and result numerical payload [f64 slots].
    /// Existing native and spatial-field budgets apply separately.
    pub max_interface_scalars: usize,
    /// Numerical reference resistance is this positive factor times |Z_ii|.
    /// It cancels from the physical circuit; changing it cannot add damping.
    pub reference_resistance_scale: f64,
    pub relative_terminal_residual_tolerance: f64,
    pub voltage_tolerance_v: f64,
    pub current_tolerance_a: f64,
    /// Absolute real AND reactive circuit/field power allowance [W or var].
    pub power_tolerance_w: f64,
}

#[derive(Debug)]
pub enum HarmonicCircuitError {
    Invalid(&'static str),
    Circuit(CircuitError),
    Native(HarmonicError),
    Field(HarmonicWindingError),
    TerminalSingular { index: usize },
    TerminalResidual { relative: f64, tolerance: f64 },
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for HarmonicCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "harmonic winding circuit: {self:?}")
    }
}
impl std::error::Error for HarmonicCircuitError {}
impl From<CircuitError> for HarmonicCircuitError {
    fn from(error: CircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<HarmonicError> for HarmonicCircuitError {
    fn from(error: HarmonicError) -> Self {
        match error {
            HarmonicError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl From<HarmonicWindingError> for HarmonicCircuitError {
    fn from(error: HarmonicWindingError) -> Self {
        match error {
            HarmonicWindingError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<LuComplexError> for HarmonicCircuitError {
    fn from(error: LuComplexError) -> Self {
        match error {
            LuComplexError::Cancelled => Self::Cancelled,
            LuComplexError::Singular { index } => Self::TerminalSingular { index },
        }
    }
}
type Result<T> = core::result::Result<T, HarmonicCircuitError>;

/// Original input/output units and orientation, retained through port expansion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicSourceKind {
    /// Input is first-minus-second terminal voltage; output is delivered current.
    Voltage,
    /// Input injects current at the first terminal; output is its voltage.
    Current,
}

/// Original graph branch, including field windings and source branches.
#[derive(Clone, Debug)]
pub struct HarmonicCircuitBranch {
    pub branch: usize,
    pub voltage_v: C64,
    /// Passive first-to-second current. Current-source branches therefore use -u.
    pub current_a: C64,
    /// `conj(I) V / 2`; sources may absorb or supply signed power.
    pub absorbed_power_va: C64,
}

/// One original source branch and its actual delivered cycle-average power.
#[derive(Clone, Debug)]
pub struct HarmonicCircuitSource {
    pub branch: usize,
    pub port: usize,
    pub kind: HarmonicSourceKind,
    /// Original peak voltage or injected-current input.
    pub input: C64,
    /// Power-conjugate delivered current or terminal voltage.
    pub output: C64,
    /// Signed supplied complex power; inductive loads have negative imaginary part.
    pub supplied_power_va: C64,
}

/// Original ideal-transformer reactions, without fictitious magnetic storage.
#[derive(Clone, Debug)]
pub struct HarmonicCircuitTransformer {
    pub transformer: usize,
    pub primary_voltage_v: C64,
    pub secondary_voltage_v: C64,
    pub primary_current_a: C64,
    pub secondary_current_a: C64,
    pub power_residual_va: C64,
}

/// Physical network and spatial-field response at the same frequency and phase.
/// Reference resistors and compensating ports do not appear as physical sources.
#[derive(Clone, Debug)]
pub struct HarmonicCircuitResponse {
    pub field: HarmonicWindingResponse,
    /// Includes ground at index zero, then the original graph's node order.
    pub node_voltages_v: Vec<C64>,
    pub branches: Vec<HarmonicCircuitBranch>,
    pub sources: Vec<HarmonicCircuitSource>,
    pub transformers: Vec<HarmonicCircuitTransformer>,
    /// Outputs summed by the original source-port index. Gaps retain zero.
    pub source_outputs: Vec<C64>,
    /// Original non-field capacitive and inductive cycle-average storage.
    pub mean_network_stored_energy_j: f64,
    pub mean_resistor_power_w: f64,
    /// Sum of independently reconstructed spatial winding/conductor heat.
    pub mean_field_power_w: f64,
    pub supplied_power_va: C64,
    pub absorbed_power_va: C64,
    /// Physical absorption minus original external source supply.
    pub power_residual_va: C64,
    pub relative_terminal_residual: f64,
    pub relative_network_residual: f64,
    pub maximum_terminal_voltage_defect_v: f64,
    pub maximum_terminal_current_defect_a: f64,
    pub maximum_kirchhoff_current_defect_a: f64,
}

fn invalid(message: &'static str) -> HarmonicCircuitError {
    HarmonicCircuitError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(HarmonicCircuitError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable harmonic circuit arithmetic"))
    }
}
fn phasor(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn power(voltage: C64, current: C64) -> Result<C64> {
    phasor((current.conj() * voltage).scale(0.5))
}

struct Source {
    branch: usize,
    port: usize,
    kind: HarmonicSourceKind,
}

// Native assembly allocates before HarmonicModel checks its budget. Admit the
// identical native dimension/payload/work expression BEFORE entering assembly.
fn admit_native_budget(n: usize, ports: usize, budget: HarmonicBudget) -> Result<()> {
    let square = n
        .checked_mul(n)
        .ok_or_else(|| invalid("native dimension overflow"))?;
    let bytes = square
        .checked_mul(88)
        .and_then(|v| {
            n.checked_mul(ports)
                .and_then(|p| p.checked_mul(16))
                .and_then(|p| v.checked_add(p))
        })
        .and_then(|v| n.checked_mul(256).and_then(|p| v.checked_add(p)))
        .and_then(|v| ports.checked_mul(64).and_then(|p| v.checked_add(p)))
        .ok_or_else(|| invalid("native payload overflow"))?;
    let work = (square as u64)
        .checked_mul(n as u64 + ports as u64 + 8)
        .ok_or_else(|| invalid("native work overflow"))?;
    if n == 0
        || n > MAX_HARMONIC_DIMENSION
        || n > budget.max_dimension
        || bytes > budget.max_dense_bytes
        || work > budget.max_dense_work
        || !positive(budget.relative_residual_tolerance)
        || budget.relative_residual_tolerance >= 1.0
    {
        return Err(invalid(
            "native circuit dimension, payload, work or accuracy budget refused",
        ));
    }
    Ok(())
}

fn terminal_solve(
    matrix: &[C64],
    rhs: &[C64],
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64)> {
    let n = rhs.len();
    let mut scaled = matrix.to_vec();
    let mut values = rhs.to_vec();
    for (row, input) in scaled.chunks_exact_mut(n).zip(&mut values) {
        poll(cancelled)?;
        let scale = row.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        if !positive(scale) {
            return Err(invalid("terminal row has no representable response"));
        }
        for z in row {
            *z = phasor(C64::new(z.re / scale, z.im / scale))?;
        }
        *input = phasor(C64::new(input.re / scale, input.im / scale))?;
    }
    let mut columns = vec![0.0; n];
    for column in 0..n {
        poll(cancelled)?;
        let scale = (0..n)
            .map(|row| scaled[row * n + column].abs())
            .fold(0.0_f64, f64::max);
        if !positive(scale) {
            return Err(invalid("terminal column has no representable response"));
        }
        columns[column] = scale;
        for row in 0..n {
            let z = &mut scaled[row * n + column];
            *z = phasor(C64::new(z.re / scale, z.im / scale))?;
        }
    }
    let factor = lu_complex_with_cancel(&scaled, n, &mut *cancelled)?;
    factor.solve_with_cancel(&mut values, &mut *cancelled)?;
    for (value, scale) in values.iter_mut().zip(columns) {
        *value = phasor(C64::new(value.re / scale, value.im / scale))?;
    }
    let mut maximum = 0.0_f64;
    for (row, &input) in rhs.iter().enumerate() {
        poll(cancelled)?;
        let mut defect = -input;
        let mut scale = input.abs();
        for (column, &value) in values.iter().enumerate() {
            let coefficient = matrix[row * n + column];
            defect = phasor(defect + coefficient * value)?;
            scale = finite(scale + coefficient.abs() * value.abs())?;
        }
        maximum = maximum.max(if scale == 0.0 {
            0.0
        } else {
            defect.abs() / scale
        });
    }
    if maximum > tolerance {
        return Err(HarmonicCircuitError::TerminalResidual {
            relative: maximum,
            tolerance,
        });
    }
    Ok((values, maximum))
}

impl HarmonicWindingFrequency<'_, '_, '_> {
    /// Drive every field winding from the original RLC, transformer and source
    /// graph. `winding_branches` is in the field model's ORIGINAL winding order;
    /// every entry replaces one distinct inductor branch, whose numeric L is a
    /// placeholder. Its resistance, mutual field and declared end storage come
    /// exclusively from this frequency's source-bound field owner.
    ///
    /// `external_mutual` may couple remaining ordinary inductors only. External
    /// graph resistors remain separate loads; do not repeat the represented
    /// winding copper there. Sources retain original indices and exp(-i omega t)
    /// peak phasors. Shared source controls of the same physical kind are allowed;
    /// one input cannot simultaneously mean volts and amperes.
    ///
    /// One original descriptor factorization serves at most 16 terminal unit
    /// responses and two drives. Only the terminal Schur matrix is additionally
    /// factored, and the stored field unit responses are reused. True singular
    /// networks, exhausted budgets and failed power/constraint checks refuse.
    /// This is a fixed-bias small-signal AC response, not a transient, nonlinear
    /// large-signal saturation cycle, moving geometry or temperature update.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn circuit_response(
        &self,
        graph: &CircuitGraph,
        winding_branches: &[usize],
        external_mutual: &[MutualInductance],
        source_phasors: &[C64],
        controls: HarmonicCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicCircuitResponse> {
        poll(cancelled)?;
        let count = self.unit_responses().len();
        if count == 0
            || count > MAX_HARMONIC_CIRCUIT_WINDINGS
            || winding_branches.len() != count
            || self.impedance_ohm().len() != count * count
            || graph.branches.len() > MAX_HARMONIC_CIRCUIT_BRANCHES
            || graph.node_count < 2
            || graph.node_count > MAX_HARMONIC_DIMENSION
            || graph.transformers.len() > MAX_HARMONIC_CIRCUIT_BRANCHES
            || external_mutual.len()
                > graph
                    .branches
                    .len()
                    .saturating_mul(graph.branches.len().saturating_sub(1))
                    / 2
            || !positive(controls.reference_resistance_scale)
            || !positive(controls.relative_terminal_residual_tolerance)
            || controls.relative_terminal_residual_tolerance >= 1.0
            || !positive(controls.voltage_tolerance_v)
            || !positive(controls.current_tolerance_a)
            || !positive(controls.power_tolerance_w)
            || source_phasors.len() > MAX_HARMONIC_CIRCUIT_SOURCES
            || source_phasors
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid(
                "bounded graph, complete winding terminals and finite circuit controls required",
            ));
        }
        let mut winding_at = vec![None; graph.branches.len()];
        let mut reference = Vec::with_capacity(count);
        for (winding, &branch) in winding_branches.iter().enumerate() {
            poll(cancelled)?;
            let &(a, b, element) = graph
                .branches
                .get(branch)
                .ok_or_else(|| invalid("field winding branch is outside original graph"))?;
            if a == b
                || a >= graph.node_count
                || b >= graph.node_count
                || !matches!(element, Branch::Inductor { .. })
                || winding_at[branch].is_some()
            {
                return Err(invalid(
                    "each field winding must replace a distinct oriented inductor branch",
                ));
            }
            winding_at[branch] = Some(winding);
            let r = self.impedance_ohm()[winding * count + winding].abs()
                * controls.reference_resistance_scale;
            if !positive(r) || !positive(1.0 / r) {
                return Err(invalid(
                    "numerical reference termination is unrepresentable",
                ));
            }
            reference.push(r);
        }
        for link in external_mutual {
            poll(cancelled)?;
            if [link.branches.0, link.branches.1]
                .iter()
                .any(|&branch| winding_at.get(branch).is_none_or(Option::is_some))
            {
                return Err(invalid(
                    "external mutual links cannot duplicate field-owned winding storage",
                ));
            }
        }
        let mut sources = Vec::new();
        let mut source_at = vec![None; graph.branches.len()];
        let mut source_kinds = vec![None; MAX_HARMONIC_CIRCUIT_SOURCES];
        let (mut ports, mut inductors, mut capacitors, mut voltage_sources) = (0, 0, 0, 0);
        for (branch, &(_, _, element)) in graph.branches.iter().enumerate() {
            poll(cancelled)?;
            let (port, kind) = match element {
                Branch::Inductor { .. } => {
                    if winding_at[branch].is_none() {
                        inductors += 1;
                    }
                    continue;
                }
                Branch::Capacitor { .. } => {
                    capacitors += 1;
                    continue;
                }
                Branch::Resistor { .. } => continue,
                Branch::VoltageSource { port } => {
                    voltage_sources += 1;
                    (port, HarmonicSourceKind::Voltage)
                }
                Branch::CurrentSource { port } => (port, HarmonicSourceKind::Current),
            };
            if port >= MAX_HARMONIC_CIRCUIT_SOURCES
                || sources.len() >= MAX_HARMONIC_CIRCUIT_SOURCES
                || source_kinds[port].is_some_and(|previous| previous != kind)
            {
                return Err(invalid(
                    "source capacity exceeded or one control mixes voltage and current units",
                ));
            }
            source_kinds[port] = Some(kind);
            ports = ports.max(port + 1);
            source_at[branch] = Some(sources.len());
            sources.push(Source { branch, port, kind });
        }
        if source_phasors.len() != ports
            || source_phasors
                .iter()
                .enumerate()
                .any(|(port, &value)| source_kinds[port].is_none() && value != C64::ZERO)
        {
            return Err(invalid(
                "finite inputs must retain exactly the original source-port mapping",
            ));
        }
        let n = inductors + 2 * capacitors + graph.node_count - 1
            + voltage_sources
            + graph.transformers.len();
        let native_ports = sources.len() + count;
        admit_native_budget(n, native_ports, controls.native)?;
        let payload = 16 * n
            + 24 * graph.branches.len()
            + 24 * sources.len()
            + 8 * graph.node_count
            + 16 * count * count
            + 32 * count
            + 10 * graph.transformers.len()
            + 8 * external_mutual.len()
            + 16 * native_ports
            + 4 * ports
            + MAX_HARMONIC_CIRCUIT_SOURCES;
        if payload > controls.max_interface_scalars {
            return Err(invalid(
                "harmonic circuit interface payload exceeds caller budget",
            ));
        }
        let mut lowered = graph.clone();
        for source in &sources {
            let port = source_at[source.branch].unwrap();
            lowered.branches[source.branch].2 = match source.kind {
                HarmonicSourceKind::Voltage => Branch::VoltageSource { port },
                HarmonicSourceKind::Current => Branch::CurrentSource { port },
            };
        }
        for (winding, &branch) in winding_branches.iter().enumerate() {
            lowered.branches[branch].2 = Branch::Resistor {
                ohms: reference[winding],
            };
            let (a, b, _) = lowered.branches[branch];
            lowered.branches.push((
                a,
                b,
                Branch::CurrentSource {
                    port: sources.len() + winding,
                },
            ));
        }
        let dae = assemble_circuit_with_mutual_inductance(&lowered, external_mutual)?;
        let native = HarmonicModel::from_descriptor(&dae.system, controls.native, &mut *cancelled)?;
        let factor = native.factor(self.omega_rad_s(), &mut *cancelled)?;
        let mut input = vec![C64::ZERO; native_ports];
        for (index, source) in sources.iter().enumerate() {
            input[index] = source_phasors[source.port];
        }
        let open = factor.response(&input, &mut *cancelled)?;
        let mut terminal = vec![C64::ZERO; count * count];
        let mut unit = vec![C64::ZERO; native_ports];
        for column in 0..count {
            poll(cancelled)?;
            unit[sources.len() + column] = C64::ONE;
            let response = factor.response(&unit, &mut *cancelled)?;
            for row in 0..count {
                terminal[row * count + column] = response.output[sources.len() + row];
            }
            unit[sources.len() + column] = C64::ZERO;
        }
        let impedance = self.impedance_ohm();
        let mut schur = terminal.clone();
        for row in 0..count {
            poll(cancelled)?;
            for inner in 0..count {
                let value = terminal[row * count + inner];
                let multiplier = phasor(
                    C64::from_re(if row == inner { 1.0 } else { 0.0 })
                        - C64::new(value.re / reference[inner], value.im / reference[inner]),
                )?;
                for column in 0..count {
                    schur[row * count + column] = phasor(
                        schur[row * count + column]
                            + multiplier * impedance[inner * count + column],
                    )?;
                }
            }
        }
        let (currents, relative_terminal_residual) = terminal_solve(
            &schur,
            &open.output[sources.len()..],
            controls.relative_terminal_residual_tolerance,
            cancelled,
        )?;
        let mut voltages = vec![C64::ZERO; count];
        for row in 0..count {
            for column in 0..count {
                voltages[row] =
                    phasor(voltages[row] + impedance[row * count + column] * currents[column])?;
            }
            input[sources.len() + row] = phasor(
                C64::new(
                    voltages[row].re / reference[row],
                    voltages[row].im / reference[row],
                ) - currents[row],
            )?;
        }
        // Final native replay audits the full original descriptor equations;
        // the field response independently reconstructs spatial heat and storage.
        let response = factor.response(&input, &mut *cancelled)?;
        let field = self.response(&currents, cancelled)?;
        let mut nodes = vec![C64::ZERO; graph.node_count];
        for (node, &index) in dae.node_potential_index.iter().enumerate() {
            nodes[node + 1] = response.effort[index];
        }
        let voltage = |a: usize, b: usize| nodes[a] - nodes[b];
        let mut maximum_voltage = 0.0_f64;
        let mut maximum_current = 0.0_f64;
        let mut reference_heat = 0.0;
        for (winding, &branch) in winding_branches.iter().enumerate() {
            poll(cancelled)?;
            let (a, b, _) = graph.branches[branch];
            let v = phasor(voltage(a, b))?;
            let reference_current = phasor(C64::new(
                v.re / reference[winding],
                v.im / reference[winding],
            ))?;
            let terminal_current = phasor(reference_current - input[sources.len() + winding])?;
            maximum_voltage = maximum_voltage.max((v - field.terminal_voltages_v[winding]).abs());
            maximum_current =
                maximum_current.max((terminal_current - field.currents_a[winding]).abs());
            reference_heat = finite(reference_heat + power(v, reference_current)?.re)?;
        }
        if maximum_voltage > controls.voltage_tolerance_v
            || maximum_current > controls.current_tolerance_a
        {
            return Err(HarmonicCircuitError::Balance(
                "field terminal relation and native circuit replay disagree",
            ));
        }
        let mut branches = Vec::with_capacity(graph.branches.len());
        let mut kcl = vec![C64::ZERO; graph.node_count];
        let (mut flux_slot, mut charge_slot) = (0, 0);
        let mut resistor_heat = 0.0;
        let mut network_absorption = C64::ZERO;
        for (branch, &(a, b, element)) in graph.branches.iter().enumerate() {
            poll(cancelled)?;
            let v = phasor(voltage(a, b))?;
            let current = if let Some(winding) = winding_at[branch] {
                currents[winding]
            } else {
                match element {
                    Branch::Resistor { ohms } => phasor(C64::new(v.re / ohms, v.im / ohms))?,
                    Branch::Inductor { .. } => {
                        let current = response.effort[dae.flux_index[flux_slot]];
                        flux_slot += 1;
                        current
                    }
                    Branch::Capacitor { .. } => {
                        let charge = response.state[dae.charge_index[charge_slot]];
                        charge_slot += 1;
                        phasor(C64::new(0.0, -self.omega_rad_s()) * charge)?
                    }
                    Branch::VoltageSource { .. } => -response.output[source_at[branch].unwrap()],
                    Branch::CurrentSource { port } => -source_phasors[port],
                }
            };
            let absorbed = power(v, current)?;
            if matches!(element, Branch::Resistor { .. }) {
                resistor_heat = finite(resistor_heat + absorbed.re)?;
            }
            if source_at[branch].is_none() && winding_at[branch].is_none() {
                network_absorption = phasor(network_absorption + absorbed)?;
            }
            kcl[a] = phasor(kcl[a] + current)?;
            kcl[b] = phasor(kcl[b] - current)?;
            branches.push(HarmonicCircuitBranch {
                branch,
                voltage_v: v,
                current_a: current,
                absorbed_power_va: absorbed,
            });
        }
        let mut transformers = Vec::with_capacity(graph.transformers.len());
        // The native CircuitDae's published state convention ends in i_T.
        // Read those original reactions; no transformer or node equations are solved here.
        let offset = dae.system.state_dim() - graph.transformers.len();
        for (index, transformer) in graph.transformers.iter().enumerate() {
            poll(cancelled)?;
            let first = phasor(response.effort[offset + index])?;
            let second = phasor(first.scale(-transformer.ratio))?;
            let vp = phasor(voltage(transformer.primary.0, transformer.primary.1))?;
            let vs = phasor(voltage(transformer.secondary.0, transformer.secondary.1))?;
            let defect = phasor(power(vp, first)? + power(vs, second)?)?;
            if defect.re.abs() > controls.power_tolerance_w
                || defect.im.abs() > controls.power_tolerance_w
                || (vp - vs.scale(transformer.ratio)).abs() > controls.voltage_tolerance_v
            {
                return Err(HarmonicCircuitError::Balance(
                    "original ideal-transformer voltage or power identity failed",
                ));
            }
            for ((a, b), current) in [
                (transformer.primary, first),
                (transformer.secondary, second),
            ] {
                kcl[a] = phasor(kcl[a] + current)?;
                kcl[b] = phasor(kcl[b] - current)?;
            }
            transformers.push(HarmonicCircuitTransformer {
                transformer: index,
                primary_voltage_v: vp,
                secondary_voltage_v: vs,
                primary_current_a: first,
                secondary_current_a: second,
                power_residual_va: defect,
            });
        }
        let maximum_kcl = kcl.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        if maximum_kcl > controls.current_tolerance_a {
            return Err(HarmonicCircuitError::Balance(
                "original branch and transformer currents fail Kirchhoff balance",
            ));
        }
        let mut source_outputs = vec![C64::ZERO; ports];
        let mut source_reports = Vec::with_capacity(sources.len());
        let mut supplied = C64::ZERO;
        for (index, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            let drive = source_phasors[source.port];
            let output = response.output[index];
            let branch = &branches[source.branch];
            let (voltage_defect, current_defect) = match source.kind {
                HarmonicSourceKind::Voltage => (
                    (drive - branch.voltage_v).abs(),
                    (output + branch.current_a).abs(),
                ),
                HarmonicSourceKind::Current => (
                    (output - branch.voltage_v).abs(),
                    (drive + branch.current_a).abs(),
                ),
            };
            if voltage_defect > controls.voltage_tolerance_v
                || current_defect > controls.current_tolerance_a
            {
                return Err(HarmonicCircuitError::Balance(
                    "original source voltage or current constraint failed",
                ));
            }
            // Real port power is symmetric; reactive power must preserve which
            // member is voltage and which is delivered current.
            let supply = match source.kind {
                HarmonicSourceKind::Voltage => power(drive, output)?,
                HarmonicSourceKind::Current => power(output, drive)?,
            };
            let defect = phasor(supply + branches[source.branch].absorbed_power_va)?;
            if defect.re.abs() > controls.power_tolerance_w
                || defect.im.abs() > controls.power_tolerance_w
            {
                return Err(HarmonicCircuitError::Balance(
                    "original source branch and native port power disagree",
                ));
            }
            supplied = phasor(supplied + supply)?;
            source_outputs[source.port] = phasor(source_outputs[source.port] + output)?;
            source_reports.push(HarmonicCircuitSource {
                branch: source.branch,
                port: source.port,
                kind: source.kind,
                input: drive,
                output,
                supplied_power_va: supply,
            });
        }
        let field_heat = finite(field.mean_dissipated_power_w())?;
        let absorbed = phasor(network_absorption + C64::new(field_heat, field.reactive_power_var))?;
        let residual = phasor(absorbed - supplied)?;
        let network_heat_defect =
            finite(response.dissipated_power_w - reference_heat - resistor_heat)?;
        if resistor_heat < 0.0
            || field_heat < 0.0
            || response.mean_stored_energy_j < 0.0
            || residual.re.abs() > controls.power_tolerance_w
            || residual.im.abs() > controls.power_tolerance_w
            || network_heat_defect.abs() > controls.power_tolerance_w
            || response.power_residual_w.abs() > controls.power_tolerance_w
        {
            return Err(HarmonicCircuitError::Balance(
                "original sources, physical network and spatial field power disagree",
            ));
        }
        poll(cancelled)?;
        Ok(HarmonicCircuitResponse {
            field,
            node_voltages_v: nodes,
            branches,
            sources: source_reports,
            transformers,
            source_outputs,
            mean_network_stored_energy_j: response.mean_stored_energy_j,
            mean_resistor_power_w: resistor_heat,
            mean_field_power_w: field_heat,
            supplied_power_va: supplied,
            absorbed_power_va: absorbed,
            power_residual_va: residual,
            relative_terminal_residual,
            relative_network_residual: response.relative_residual,
            maximum_terminal_voltage_defect_v: maximum_voltage,
            maximum_terminal_current_defect_a: maximum_current,
            maximum_kirchhoff_current_defect_a: maximum_kcl,
        })
    }
}
