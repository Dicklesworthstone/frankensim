//! Finite-energy capacitor terminals on the original bidirectional converter.
//!
//! Each bus connects through its own resistor to a retained external capacitor.
//! Its voltage sags or rises with transferred charge; it is never replaced by
//! an ideal voltage command. Two optional held currents charge those external
//! capacitors. Zero currents close the complete five-state electrical system.
//!
//! This is a linear capacitor/supercapacitor equivalent with frozen C and ESR,
//! not battery chemistry, an SOC/OCV curve, leakage or a protection controller.
//! Signed capacitor voltage and regenerative current remain physical coordinates;
//! no clipping or energy repair hides a trajectory outside caller validity.

use super::{
    AveragedBuckBoost, BuckBoostParameters, CONDUCTION_CHANNELS, ConductionHeat, ConverterStep,
    invalid, poll, positive,
};
use crate::linear_exact::prepared::PreparedQuadraticFlow;
use crate::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedSourceLoss, split_quadratic,
};
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage,
    storage_energy_change,
};
use fs_math::{Dd, det};

/// Two finite external capacitors, in input/output bus order. Their charges are
/// additional continuing storage, not replacements for the original bus charge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapacitorTerminals {
    pub capacitance_f: [f64; 2],
    /// Physical positive series resistance between each capacitor and its bus.
    pub resistance_ohm: [f64; 2],
}

/// Independent physical accounts on one complete five-state trajectory.
#[derive(Clone, Debug)]
pub struct CapacitorFedStep {
    /// Original three-state converter and its five conduction heat owners.
    pub converter: ConverterStep,
    /// Complete converter + external capacitor boundary. Dissipation includes
    /// the same five converter heats plus both terminal heats, never extra heat.
    pub total: PreparedStepRecord,
    pub terminal_storage_change_j: [f64; 2],
    /// Signed work leaving each external capacitor through its series resistor.
    /// This is an independent integral of V_cap * I_to_bus, not an energy deficit.
    pub terminal_outflow_work_j: [f64; 2],
    /// Signed work of the held external charging current, I_ext * integral V_cap.
    pub external_work_j: [f64; 2],
    pub terminal_resistor_heat_j: [f64; 2],
}

/// State order: original flux, original input/output bus charges, then input/
/// output terminal-capacitor charges. Both external port currents enter storage.
pub struct CapacitorFedBuckBoost {
    converter: AveragedBuckBoost,
    terminals: CapacitorTerminals,
    system: PortHamiltonian,
}

impl AveragedBuckBoost {
    /// Retain both external capacitor charges and their reciprocal resistive
    /// connections in the same original pHS. With g=1/R each connection stamps
    /// g*[1,-1]*[1,-1]^T, giving positive g*(V_cap-V_bus)^2 dissipation.
    /// No source voltage is frozen or sampled at a step boundary.
    pub fn with_capacitor_terminals(
        self,
        terminals: CapacitorTerminals,
    ) -> Result<CapacitorFedBuckBoost, PhsError> {
        if terminals
            .capacitance_f
            .iter()
            .chain(&terminals.resistance_ohm)
            .any(|&v| !positive(v) || !positive(1.0 / v))
        {
            return Err(invalid(
                "capacitor terminals require positive representable C and R",
            ));
        }
        let (old_j, old_r, _) = self.system.structure();
        let old_q = self
            .system
            .storage()
            .quadratic_part()
            .ok_or_else(|| invalid("converter requires its original quadratic storage"))?;
        let (mut j, mut r, mut q, mut g) =
            (vec![0.0; 25], vec![0.0; 25], vec![0.0; 25], vec![0.0; 10]);
        for i in 0..3 {
            for k in 0..3 {
                j[i * 5 + k] = old_j[i * 3 + k];
                r[i * 5 + k] = old_r[i * 3 + k];
                q[i * 5 + k] = old_q.matrix[i * 3 + k];
            }
        }
        for p in 0..2 {
            let (b, s) = (p + 1, p + 3);
            let conductance = 1.0 / terminals.resistance_ohm[p];
            r[b * 5 + b] += conductance;
            r[s * 5 + s] += conductance;
            r[b * 5 + s] -= conductance;
            r[s * 5 + b] -= conductance;
            q[s * 5 + s] = 1.0 / terminals.capacitance_f[p];
            g[s * 2 + p] = 1.0;
        }
        let system = PortHamiltonian::new(5, 2, j, r, g, Box::new(QuadraticStorage::new(q, 5)?))?;
        Ok(CapacitorFedBuckBoost {
            converter: self,
            terminals,
            system,
        })
    }
}

impl CapacitorFedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn converter_system(&self) -> &PortHamiltonian {
        self.converter.system()
    }
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.converter.parameters()
    }
    pub const fn terminals(&self) -> CapacitorTerminals {
        self.terminals
    }

    /// Original charge/flux coordinates from physical current and four voltages.
    pub fn state(
        &self,
        current_a: f64,
        bus_voltage_v: [f64; 2],
        terminal_voltage_v: [f64; 2],
    ) -> Result<[f64; 5], PhsError> {
        let base = self.converter.state(current_a, bus_voltage_v)?;
        let state = [
            base[0],
            base[1],
            base[2],
            self.terminals.capacitance_f[0] * terminal_voltage_v[0],
            self.terminals.capacitance_f[1] * terminal_voltage_v[1],
        ];
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid(
                "finite capacitor state or stored energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    /// Cache exact evolution, signed physical work and positive heat on the
    /// SAME five-state generator. The dense allowance applies to one 8-coordinate
    /// quadratic-flow cache and each of three 10-coordinate positive observers.
    /// Warm calls allocate nothing and perform no factorization or exponential.
    pub fn prepare_controlled(
        self,
        duration_s: f64,
        budget: LinearSplitBudget,
        observation_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCapacitorFedBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s) {
            return Err(invalid("finite capacitor step requires positive duration").into());
        }
        let split = split_quadratic(&self.system, budget)?;
        poll(&mut cancelled)?;
        let q = split.quadratic_matrix();
        let (j, r, input) = self.system.structure();
        let generator: [f64; 25] = core::array::from_fn(|slot| {
            let (i, k) = (slot / 5, slot % 5);
            (0..5)
                .map(|a| (j[i * 5 + a] - r[i * 5 + a]) * q[a * 5 + k])
                .sum()
        });
        // Independent signed observations: bus input power, capacitor outflow
        // power. Their difference is checked against a separately positive
        // resistor-square observation; no balance remainder supplies any heat.
        let mut forms = vec![vec![0.0; 25]; 4];
        let mut voltages = vec![vec![0.0; 5]; 4];
        for p in 0..2 {
            let (b, s) = (p + 1, p + 3);
            let conductance = 1.0 / self.terminals.resistance_ohm[p];
            let (qb, qs) = (q[b * 5 + b], q[s * 5 + s]);
            let cross = 0.5 * conductance * qb * qs;
            forms[p][b * 5 + b] = -conductance * qb * qb;
            forms[p][b * 5 + s] = cross;
            forms[p][s * 5 + b] = cross;
            forms[p + 2][s * 5 + s] = conductance * qs * qs;
            forms[p + 2][b * 5 + s] = -cross;
            forms[p + 2][s * 5 + b] = -cross;
            voltages[p][b] = qb;
            voltages[p + 2][s] = qs;
        }
        let flow = PreparedQuadraticFlow::new_controlled(
            &self.system,
            &split,
            duration_s,
            &forms,
            &voltages,
            observation_budget,
            &mut cancelled,
        )?;
        let total_r: f64 = self.converter.effective_resistance_ohm.iter().sum();
        let fractions = self
            .converter
            .effective_resistance_ohm
            .map(|v| if total_r > 0.0 { v / total_r } else { 0.0 });
        let prepare_heat = |row: &[f64; 5], cancel: &mut _| {
            PreparedSourceLoss::new_controlled(
                5,
                2,
                &generator,
                input,
                row,
                duration_s,
                observation_budget,
                cancel,
            )
        };
        let mut winding = [0.0; 5];
        winding[0] = det::sqrt(total_r) * q[0];
        let mut terminal = [[0.0; 5]; 2];
        for p in 0..2 {
            let (b, s) = (p + 1, p + 3);
            let root = det::sqrt(1.0 / self.terminals.resistance_ohm[p]);
            terminal[p][b] = -root * q[b * 5 + b];
            terminal[p][s] = root * q[s * 5 + s];
        }
        let heat = [
            prepare_heat(&winding, &mut cancelled)?,
            prepare_heat(&terminal[0], &mut cancelled)?,
            prepare_heat(&terminal[1], &mut cancelled)?,
        ];
        poll(&mut cancelled)?;
        Ok(PreparedCapacitorFedBuckBoost {
            model: self,
            flow,
            heat,
            fractions,
            duration_s,
        })
    }
}

pub struct PreparedCapacitorFedBuckBoost {
    model: CapacitorFedBuckBoost,
    flow: PreparedQuadraticFlow,
    heat: [PreparedSourceLoss; 3],
    fractions: [f64; CONDUCTION_CHANNELS],
    duration_s: f64,
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}

impl PreparedCapacitorFedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn converter_system(&self) -> &PortHamiltonian {
        self.model.converter_system()
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn terminals(&self) -> CapacitorTerminals {
        self.model.terminals
    }

    /// Advance continuing capacitor/inductor storage. Inputs are held currents
    /// INTO the external capacitors. The four means are input/output bus voltage,
    /// then input/output external-capacitor voltage. No endpoint voltage supplies
    /// interval work. Cancellation or any failed physical balance leaves BOTH
    /// caller buffers unchanged; a failed trial is safe to retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 5],
        external_current_a: &[f64; 2],
        next: &mut [f64; 5],
        average_voltage_v: &mut [f64; 4],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitorFedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let stationary = self.heat[0].has_exactly_zero_derivative_controlled(
            initial,
            external_current_a,
            &mut cancelled,
        )?;
        let (mut candidate, mut integrals) = ([0.0; 5], [0.0; 8]);
        let mut heats = [0.0; 3];
        for (p, observer) in self.heat.iter().enumerate() {
            heats[p] = observer.evaluate_into_controlled(
                initial,
                external_current_a,
                &mut [0.0; 10],
                &mut cancelled,
            )?;
        }
        if stationary {
            candidate = *initial;
            let mut effort = [0.0; 5];
            self.system().storage().gradient(initial, &mut effort);
            for p in 0..2 {
                let current =
                    (effort[p + 3] - effort[p + 1]) / self.model.terminals.resistance_ohm[p];
                integrals[p] = self.duration_s * effort[p + 1] * current;
                integrals[p + 2] = self.duration_s * effort[p + 3] * current;
                integrals[p + 4] = self.duration_s * effort[p + 1];
                integrals[p + 6] = self.duration_s * effort[p + 3];
            }
        } else {
            self.flow.evaluate_into(
                initial,
                external_current_a,
                &mut [0.0; 8],
                &mut [0.0; 5],
                &mut candidate,
                &mut integrals,
            )?;
        }
        // An isolated bus/capacitor pair can be at exact equilibrium while the
        // other pair or the inductor evolves. Preserve only a proved zero-flow
        // spectator: its bridge is disconnected, its source is zero, the
        // independent positive heat is zero, and the projected ORIGINAL
        // generator has an exact-zero derivative witness. No tolerance or
        // threshold turns a weak physical transfer into equilibrium.
        let modulation = self.model.converter.modulation();
        let connection = [
            modulation.input_high_fraction,
            modulation.output_high_fraction,
        ];
        for p in 0..2 {
            if connection[p] == 0.0 && external_current_a[p] == 0.0 && heats[p + 1] == 0.0 {
                let mut projected = [0.0; 5];
                projected[p + 1] = initial[p + 1];
                projected[p + 3] = initial[p + 3];
                let mut effort = [0.0; 5];
                self.system().storage().gradient(&projected, &mut effort);
                if effort[p + 1] == effort[p + 3]
                    && self.heat[0].has_exactly_zero_derivative_controlled(
                        &projected,
                        &[0.0; 2],
                        &mut cancelled,
                    )?
                {
                    candidate[p + 1] = initial[p + 1];
                    candidate[p + 3] = initial[p + 3];
                    integrals[p] = 0.0;
                    integrals[p + 2] = 0.0;
                    integrals[p + 4] = self.duration_s * effort[p + 1];
                    integrals[p + 6] = self.duration_s * effort[p + 3];
                }
            }
        }
        poll(&mut cancelled)?;
        let original = [initial[0], initial[1], initial[2]];
        let endpoint = [candidate[0], candidate[1], candidate[2]];
        let converter_change = storage_energy_change(
            self.converter_system().storage(),
            &original,
            &endpoint,
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let total_change = storage_energy_change(
            self.system().storage(),
            initial,
            &candidate,
            &mut [0.0; 5],
            &mut [0.0; 5],
        );
        let terminal_storage_change_j = core::array::from_fn(|p| {
            let before = Dd::from_f64(initial[p + 3]);
            let after = Dd::from_f64(candidate[p + 3]);
            ((after - before) * (after * Dd::from_f64(0.5) + before * Dd::from_f64(0.5))
                / Dd::from_f64(self.model.terminals.capacitance_f[p]))
            .to_f64()
        });
        let terminal_work_j = [integrals[0], integrals[1]];
        let terminal_outflow_work_j = [integrals[2], integrals[3]];
        let external_work_j = core::array::from_fn(|p| external_current_a[p] * integrals[p + 6]);
        let terminal_resistor_heat_j = [heats[1], heats[2]];
        let joules = self.fractions.map(|f| f * heats[0]);
        let voltage: [f64; 4] = core::array::from_fn(|p| integrals[p + 4] / self.duration_s);
        let supplied = sum(external_work_j);
        let dissipated = sum(heats);
        let converter_supplied = sum(terminal_work_j);
        let mut residuals = [0.0; 8];
        residuals[0] = sum([total_change, dissipated, -supplied]);
        residuals[1] = sum([converter_change, heats[0], -converter_supplied]);
        residuals[2] = sum([
            converter_change,
            terminal_storage_change_j[0],
            terminal_storage_change_j[1],
            -total_change,
        ]);
        residuals[3] = sum(joules.into_iter().chain([-heats[0]]));
        for p in 0..2 {
            residuals[4 + p] = sum([
                terminal_storage_change_j[p],
                terminal_outflow_work_j[p],
                -external_work_j[p],
            ]);
            residuals[6 + p] = sum([
                terminal_outflow_work_j[p],
                -terminal_work_j[p],
                -terminal_resistor_heat_j[p],
            ]);
        }
        // Each boundary has its OWN transfer scale. Neither untouched stored
        // energy nor a much stronger unrelated terminal can admit a weak,
        // unresolved transfer. General signed observations may refuse at such
        // ill-conditioned states; no balance remainder replaces their work.
        let scale = sum(external_work_j
            .into_iter()
            .chain(terminal_work_j)
            .chain(terminal_outflow_work_j)
            .map(f64::abs))
        .max(dissipated)
        .max(total_change.abs())
        .max(converter_change.abs())
        .max(sum(terminal_storage_change_j.map(f64::abs)));
        let mut scales = [scale; 8];
        scales[1] = sum(terminal_work_j.map(f64::abs))
            .max(heats[0])
            .max(converter_change.abs());
        scales[2] = sum([
            converter_change.abs(),
            terminal_storage_change_j[0].abs(),
            terminal_storage_change_j[1].abs(),
            total_change.abs(),
        ]);
        scales[3] = heats[0];
        for p in 0..2 {
            scales[4 + p] = sum([
                terminal_storage_change_j[p].abs(),
                terminal_outflow_work_j[p].abs(),
                external_work_j[p].abs(),
            ]);
            scales[6 + p] = sum([
                terminal_outflow_work_j[p].abs(),
                terminal_work_j[p].abs(),
                terminal_resistor_heat_j[p],
            ]);
        }
        if scales.iter().any(|v| !v.is_finite())
            || !self.system().hamiltonian(initial).is_finite()
            || !self.system().hamiltonian(&candidate).is_finite()
            || candidate
                .iter()
                .chain(&integrals)
                .chain(&voltage)
                .chain(&joules)
                .chain(&heats)
                .chain(&external_work_j)
                .chain(&terminal_storage_change_j)
                .any(|v| !v.is_finite())
            || heats.iter().chain(&joules).any(|v| *v < 0.0)
            || residuals
                .iter()
                .zip(scales)
                .any(|(v, scale)| !v.is_finite() || v.abs() > 1e-12 * scale)
        {
            return Err(
                invalid("finite capacitor converter work, storage and heat do not close").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        average_voltage_v.copy_from_slice(&voltage);
        let record = |delta_h, dissipated, supplied| PreparedStepRecord {
            delta_h,
            dissipated,
            supplied,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        };
        Ok(CapacitorFedStep {
            converter: ConverterStep {
                native: record(converter_change, heats[0], converter_supplied),
                terminal_work_j,
                conduction: ConductionHeat { joules },
            },
            total: record(total_change, dissipated, supplied),
            terminal_storage_change_j,
            terminal_outflow_work_j,
            external_work_j,
            terminal_resistor_heat_j,
        })
    }
}
