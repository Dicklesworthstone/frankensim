//! Finite battery storage with a reciprocal current/voltage bus terminal.
//!
//! State is `[remaining chemical charge, polarization charge, bus charge]` in
//! coulombs. The source-supplied OCV/SOC curve owns chemical energy through
//! `H_chem = integral OCV(q) dq`. A series resistance, one polarization RC pair
//! and a retained bus capacitor give voltage sag, relaxation and regeneration.
//! Current entering the bus is positive; a load draws negative port current.
//!
//! The same effort drives chemical depletion, polarization and bus charging:
//! `I_cell = (OCV - V_pol - V_bus)/R_series`. Thus the continuous power identity
//! is `H_dot = V_bus*I_port - R_series*I_cell² - V_pol²/R_pol`. The source
//! squares are retained for native composition. The componentwise discrete
//! gradient and its analytic derivative preserve this identity in the native
//! time-step quadrature without constructing or correcting energy after a step.
//!
//! The OCV curve is a frozen, reversible, caller-owned equivalent circuit law,
//! not a chemistry fit. SOC is bounded: exhaustion/overcharge refuses, never
//! clips charge or invents an ideal supply. Protection events, thermal drift,
//! hysteresis, aging and charge-rate capacity effects require separate laws.
//! Native discrete heat is not an exact continuous RMS-current integral; time
//! refinement remains necessary for resolved transient voltage and losses.

use std::sync::Arc;

use crate::{
    ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, PortHamiltonian,
    PreparedStepError, PreparedStepRecord, ResistanceWitness, SolveBudget, StepWorkspace, Storage,
};
use fs_math::{dd::Dd, det};

fn invalid(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}

#[derive(Clone, Debug)]
struct CurveData {
    capacity_c: f64,
    knots: Vec<(f64, f64)>,
    slopes: Vec<f64>,
    energy: Vec<f64>,
}

/// Bounded piecewise-linear open-circuit voltage versus remaining charge.
/// Flat voltage plateaus are permitted. Clones share immutable source data.
#[derive(Clone, Debug)]
pub struct BatteryOcvCurve(Arc<CurveData>);

impl BatteryOcvCurve {
    /// `capacity_c` is usable capacity in coulombs (3600 C = 1 Ah). Knots are
    /// `(SOC, volts)`, start at SOC 0, end at SOC 1, and strictly increase in SOC.
    /// Voltage must be finite, positive and nondecreasing. Discrete-gradient
    /// Newton tangents use the right segment at knots and the last at full
    /// charge. Continuous tangents require a two-sided derivative: unequal-slope
    /// knots and finite support endpoints refuse; collinear interior knots pass.
    pub fn new(capacity_c: f64, soc_voltage: &[(f64, f64)]) -> Result<Self, PhsError> {
        if !capacity_c.is_finite()
            || capacity_c <= 0.0
            || !(2..=4096).contains(&soc_voltage.len())
            || soc_voltage[0].0 != 0.0
            || soc_voltage[soc_voltage.len() - 1].0 != 1.0
        {
            return Err(invalid(
                "battery requires positive capacity and 2..4096 SOC knots covering [0,1]",
            ));
        }
        let mut knots: Vec<(f64, f64)> = Vec::with_capacity(soc_voltage.len());
        let mut slopes = Vec::with_capacity(soc_voltage.len() - 1);
        let mut energy = Vec::with_capacity(soc_voltage.len());
        let mut accumulated = Dd::ZERO;
        for &(soc, voltage) in soc_voltage {
            let charge = soc * capacity_c;
            if !soc.is_finite()
                || !(0.0..=1.0).contains(&soc)
                || !voltage.is_finite()
                || voltage <= 0.0
                || !charge.is_finite()
            {
                return Err(invalid(
                    "battery OCV knots must have finite SOC and positive voltage",
                ));
            }
            if let Some(&(previous_q, previous_v)) = knots.last() {
                let width = charge - previous_q;
                let slope = (voltage - previous_v) / width;
                if width <= 0.0 || !slope.is_finite() || slope < 0.0 {
                    return Err(invalid(
                        "battery charge knots must increase and OCV must not decrease",
                    ));
                }
                slopes.push(slope);
                accumulated = accumulated
                    + Dd::from_f64(width) * Dd::from_f64(f64::midpoint(previous_v, voltage));
            }
            if !accumulated.to_f64().is_finite() {
                return Err(invalid("battery chemical energy must be representable"));
            }
            energy.push(accumulated.to_f64());
            knots.push((charge, voltage));
        }
        Ok(Self(Arc::new(CurveData {
            capacity_c,
            knots,
            slopes,
            energy,
        })))
    }

    pub fn capacity_c(&self) -> f64 {
        self.0.capacity_c
    }

    fn segment(&self, charge: f64) -> Result<usize, PhsError> {
        if !charge.is_finite() || !(0.0..=self.capacity_c()).contains(&charge) {
            return Err(invalid("battery charge is outside supplied SOC support"));
        }
        Ok(self
            .0
            .knots
            .partition_point(|&(q, _)| q <= charge)
            .saturating_sub(1)
            .min(self.0.slopes.len() - 1))
    }

    fn voltage_in(&self, charge: f64, segment: usize) -> f64 {
        let (q, v) = self.0.knots[segment];
        v + self.0.slopes[segment] * (charge - q)
    }

    fn continuous_slope(&self, charge: f64) -> Option<f64> {
        let i = self.segment(charge).ok()?;
        if charge <= 0.0
            || charge >= self.capacity_c()
            || (i > 0 && charge == self.0.knots[i].0 && self.0.slopes[i - 1] != self.0.slopes[i])
        {
            return None;
        }
        Some(self.0.slopes[i])
    }

    pub fn voltage_v(&self, charge_c: f64) -> Result<f64, PhsError> {
        Ok(self.voltage_in(charge_c, self.segment(charge_c)?))
    }

    /// Remaining reversible chemical energy above the empty-state reference.
    /// In general this is not charge times terminal voltage.
    pub fn energy_j(&self, charge_c: f64) -> Result<f64, PhsError> {
        let i = self.segment(charge_c)?;
        let (q, v) = self.0.knots[i];
        Ok((Dd::from_f64(self.0.energy[i])
            + Dd::from_f64(charge_c - q)
                * Dd::from_f64(f64::midpoint(v, self.voltage_in(charge_c, i))))
        .to_f64())
    }

    fn intervals(
        &self,
        a: f64,
        b: f64,
        mut visit: impl FnMut(f64, f64, usize),
    ) -> Result<(), PhsError> {
        self.segment(a)?;
        self.segment(b)?;
        let (mut left, right) = (a.min(b), a.max(b));
        let mut i = self.segment(left)?;
        while left < right {
            let end = right.min(self.0.knots[i + 1].0);
            visit(left, end, i);
            left = end;
            i += 1;
        }
        Ok(())
    }

    /// Exact integral of the supplied piecewise-linear source between two
    /// charge states. Direct local integration avoids subtracting stored energy.
    pub fn energy_change_j(&self, a: f64, b: f64) -> Result<f64, PhsError> {
        let mut integral = Dd::ZERO;
        self.intervals(a, b, |left, right, i| {
            integral = integral
                + Dd::from_f64(right - left)
                    * Dd::from_f64(f64::midpoint(
                        self.voltage_in(left, i),
                        self.voltage_in(right, i),
                    ));
        })?;
        Ok(if b < a {
            -integral.to_f64()
        } else {
            integral.to_f64()
        })
    }

    fn mean_voltage(&self, a: f64, b: f64) -> Result<f64, PhsError> {
        let i = self.segment(a)?;
        if i == self.segment(b)? {
            return Ok(self.voltage_in(f64::midpoint(a, b), i));
        }
        Ok(self.energy_change_j(a, b)? / (b - a))
    }

    fn mean_tangent(&self, a: f64, b: f64) -> Result<f64, PhsError> {
        let i = self.segment(a)?;
        if a == b {
            return Ok(0.5 * self.0.slopes[i]);
        }
        // Integral t*V'(a+t*(b-a)) dt; avoids subtracting nearly equal OCVs.
        let width = (b - a).abs();
        let mut tangent = Dd::ZERO;
        self.intervals(a, b, |left, right, segment| {
            let t = f64::midpoint((left - a) / (b - a), (right - a) / (b - a));
            tangent = tangent
                + Dd::from_f64(self.0.slopes[segment])
                    * Dd::from_f64(t)
                    * Dd::from_f64((right - left) / width);
        })?;
        Ok(tangent.to_f64())
    }

    fn remainder(&self, a: f64, b: f64, midpoint: f64) -> Result<f64, PhsError> {
        let midpoint_v = self.voltage_v(midpoint)?;
        let mut remainder = Dd::ZERO;
        self.intervals(a, b, |left, right, i| {
            let average = Dd::from_f64(0.5) * Dd::from_f64(self.voltage_in(left, i))
                + Dd::from_f64(0.5) * Dd::from_f64(self.voltage_in(right, i));
            remainder =
                remainder + Dd::from_f64(right - left) * (average - Dd::from_f64(midpoint_v));
        })?;
        Ok(if b < a {
            -remainder.to_f64()
        } else {
            remainder.to_f64()
        })
    }
}

/// Frozen battery equivalent circuit parameters in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryParameters {
    pub series_resistance_ohm: f64,
    pub polarization_resistance_ohm: f64,
    pub polarization_capacitance_f: f64,
    /// Physical bus storage: keeps terminal voltage a retained dynamic effort.
    pub bus_capacitance_f: f64,
}

#[derive(Clone)]
struct BatteryStorage {
    curve: BatteryOcvCurve,
    parameters: BatteryParameters,
}

impl Storage for BatteryStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(self.clone()))
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if x.len() != 3 {
            return f64::NAN;
        }
        let Ok(chemical) = self.curve.energy_j(x[0]) else {
            return f64::NAN;
        };
        (Dd::from_f64(chemical)
            + Dd::from_f64(0.5) * Dd::from_f64(x[1]) * Dd::from_f64(x[1])
                / Dd::from_f64(self.parameters.polarization_capacitance_f)
            + Dd::from_f64(0.5) * Dd::from_f64(x[2]) * Dd::from_f64(x[2])
                / Dd::from_f64(self.parameters.bus_capacitance_f))
        .to_f64()
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if x.len() != 3 || out.len() != 3 {
            out.fill(f64::NAN);
            return;
        }
        out[0] = self.curve.voltage_v(x[0]).unwrap_or(f64::NAN);
        out[1] = x[1] / self.parameters.polarization_capacitance_f;
        out[2] = x[2] / self.parameters.bus_capacitance_f;
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
    fn coordinate_scale(&self, coordinate: usize) -> Option<f64> {
        match coordinate {
            0 => Some(self.curve.capacity_c()),
            1 => Some(self.parameters.polarization_capacitance_f * self.curve.0.knots.last()?.1),
            2 => Some(self.parameters.bus_capacitance_f * self.curve.0.knots.last()?.1),
            _ => None,
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        if a.len() != 3 || b.len() != 3 || midpoint.len() != 3 {
            return Some(f64::NAN);
        }
        Some(
            self.curve
                .remainder(a[0], b[0], midpoint[0])
                .unwrap_or(f64::NAN),
        )
    }
}

impl ContinuousStorageTangent for BatteryStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        if state.len() != 3
            || direction.len() != 3
            || out.len() != 3
            || !state.iter().chain(direction).all(|v| v.is_finite())
        {
            return false;
        }
        let Some(slope) = self.curve.continuous_slope(state[0]) else {
            return false;
        };
        out[0] = slope * direction[0];
        out[1] = direction[1] / self.parameters.polarization_capacitance_f;
        out[2] = direction[2] / self.parameters.bus_capacitance_f;
        out.iter().all(|v| v.is_finite())
    }
}

impl ConstitutiveDiscreteGradient for BatteryStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        if a.len() != 3 || b.len() != 3 || midpoint.len() != 3 || out.len() != 3 {
            out.fill(f64::NAN);
            return;
        }
        out[0] = self.curve.mean_voltage(a[0], b[0]).unwrap_or(f64::NAN);
        out[1] = midpoint[1] / self.parameters.polarization_capacitance_f;
        out[2] = midpoint[2] / self.parameters.bus_capacitance_f;
    }
    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if a.len() != 3
            || b.len() != 3
            || direction.len() != 3
            || out.len() != 3
            || !a.iter().chain(b).chain(direction).all(|v| v.is_finite())
        {
            return false;
        }
        let Ok(tangent) = self.curve.mean_tangent(a[0], b[0]) else {
            return false;
        };
        out[0] = tangent * direction[0];
        out[1] = 0.5 * direction[1] / self.parameters.polarization_capacitance_f;
        out[2] = 0.5 * direction[2] / self.parameters.bus_capacitance_f;
        out.iter().all(|v| v.is_finite())
    }
}

/// Instantaneous original-source quantities; heat rates are nonnegative watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryObservation {
    pub state_of_charge: f64,
    pub open_circuit_voltage_v: f64,
    pub terminal_voltage_v: f64,
    pub polarization_voltage_v: f64,
    /// Positive discharges chemical storage, negative regenerates it.
    pub cell_current_a: f64,
    pub series_heat_w: f64,
    pub polarization_heat_w: f64,
}

/// Original component work over one discrete-gradient interval. These are the
/// native time quadrature's signed work and nonnegative branch heats, not exact
/// continuous RMS integrals. The chemical effort integrates every crossed OCV
/// segment; evaluating instantaneous OCV at the midpoint is not equivalent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryInterval {
    pub chemical_energy_change_j: f64,
    pub polarization_energy_change_j: f64,
    pub bus_energy_change_j: f64,
    /// Positive work enters the battery's physical bus terminal.
    pub terminal_work_j: f64,
    pub series_heat_j: f64,
    pub polarization_heat_j: f64,
    pub discrete_ocv_v: f64,
    pub discrete_polarization_voltage_v: f64,
    pub discrete_terminal_voltage_v: f64,
    /// Positive discharges chemical storage; negative charges it.
    pub discrete_cell_current_a: f64,
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}

impl BatteryInterval {
    pub fn storage_change_j(&self) -> f64 {
        sum([
            self.chemical_energy_change_j,
            self.polarization_energy_change_j,
            self.bus_energy_change_j,
        ])
    }
    pub fn heat_j(&self) -> f64 {
        sum([self.series_heat_j, self.polarization_heat_j])
    }
    /// Independent component balance; a diagnostic of the accepted numerical
    /// solve, not a correction to storage, source work or either heat channel.
    pub fn balance_residual_j(&self) -> f64 {
        sum([
            self.chemical_energy_change_j,
            self.polarization_energy_change_j,
            self.bus_energy_change_j,
            self.series_heat_j,
            self.polarization_heat_j,
            -self.terminal_work_j,
        ])
    }
}

/// Accepted original three-state solve and its independent physical accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryStep {
    pub native: PreparedStepRecord,
    pub battery: BatteryInterval,
}

/// A finite source with one native current-input/voltage-output bus port.
/// Consume with `into_system` for reciprocal connection to a motor or converter,
/// or reuse `StepWorkspace` for prescribed external current histories.
pub struct BatterySource {
    system: PortHamiltonian,
    curve: BatteryOcvCurve,
    parameters: BatteryParameters,
}

impl BatterySource {
    pub fn new(curve: BatteryOcvCurve, parameters: BatteryParameters) -> Result<Self, PhsError> {
        let p = parameters;
        let values = [
            p.series_resistance_ohm,
            p.polarization_resistance_ohm,
            p.polarization_capacitance_f,
            p.bus_capacitance_f,
        ];
        if !values
            .iter()
            .all(|v| v.is_finite() && *v > 0.0 && v.recip().is_finite())
        {
            return Err(invalid(
                "battery resistances and capacitances must be finite and positive",
            ));
        }
        let g = p.series_resistance_ohm.recip();
        let gp = p.polarization_resistance_ohm.recip();
        let r = vec![g, -g, -g, -g, g + gp, g, -g, g, g];
        let root_g = det::sqrt(g);
        let witness = ResistanceWitness::try_from_source_rows(
            3,
            r.clone(),
            vec![root_g, -root_g, -root_g, 0.0, det::sqrt(gp), 0.0],
        )?;
        let storage = BatteryStorage {
            curve: curve.clone(),
            parameters,
        };
        let system = PortHamiltonian::new(
            3,
            1,
            vec![0.0; 9],
            r,
            vec![0.0, 0.0, 1.0],
            Box::new(storage),
        )?
        .with_resistance_witness(witness)?;
        Ok(Self {
            system,
            curve,
            parameters,
        })
    }

    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }
    pub fn curve(&self) -> &BatteryOcvCurve {
        &self.curve
    }
    pub fn parameters(&self) -> BatteryParameters {
        self.parameters
    }

    // Cold composite preparation retains the same immutable constitutive law.
    pub(crate) fn share_system(&mut self) -> PortHamiltonian {
        self.system.share_storage()
    }

    /// Initialize retained charge without assuming a precharged or empty bus.
    pub fn state(&self, soc: f64, polarization_v: f64, bus_v: f64) -> Result<[f64; 3], PhsError> {
        let state = [
            soc * self.curve.capacity_c(),
            polarization_v * self.parameters.polarization_capacitance_f,
            bus_v * self.parameters.bus_capacitance_f,
        ];
        if !soc.is_finite()
            || !(0.0..=1.0).contains(&soc)
            || !state.iter().all(|v| v.is_finite())
            || !self.system.storage().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "battery initial state must be finite and inside SOC support",
            ));
        }
        Ok(state)
    }

    pub fn observe(&self, state: &[f64; 3]) -> Result<BatteryObservation, PhsError> {
        if !state.iter().all(|v| v.is_finite()) {
            return Err(invalid("battery observation requires finite state"));
        }
        let ocv = self.curve.voltage_v(state[0])?;
        let vp = state[1] / self.parameters.polarization_capacitance_f;
        let vb = state[2] / self.parameters.bus_capacitance_f;
        let current = ((Dd::from_f64(ocv) - Dd::from_f64(vp) - Dd::from_f64(vb))
            / Dd::from_f64(self.parameters.series_resistance_ohm))
        .to_f64();
        let series_heat = self.parameters.series_resistance_ohm * current * current;
        let polarization_heat = (vp / self.parameters.polarization_resistance_ohm) * vp;
        if ![vp, vb, current, series_heat, polarization_heat]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(invalid("battery observation is not representable"));
        }
        Ok(BatteryObservation {
            state_of_charge: state[0] / self.curve.capacity_c(),
            open_circuit_voltage_v: ocv,
            terminal_voltage_v: vb,
            polarization_voltage_v: vp,
            cell_current_a: current,
            series_heat_w: series_heat,
            polarization_heat_w: polarization_heat,
        })
    }

    /// Account for this source's part of an accepted native interval. This
    /// also serves monolithic electrical compositions without another solve.
    /// Every effort is the original componentwise discrete gradient, including
    /// crossed OCV knots. The supplied terminal current is its discrete current
    /// entering this battery (held for standalone stepping, solved for a wire).
    /// Arbitrary endpoint pairs are not thereby certified as a trajectory:
    /// `balance_residual_j` exposes their physical equation defect.
    pub fn interval(
        &self,
        initial: &[f64; 3],
        next: &[f64; 3],
        terminal_current_a: f64,
        duration_s: f64,
    ) -> Result<BatteryInterval, PhsError> {
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !terminal_current_a.is_finite()
            || !self.system.hamiltonian(initial).is_finite()
            || !self.system.hamiltonian(next).is_finite()
        {
            return Err(invalid(
                "battery interval requires finite supported states, current and positive duration",
            ));
        }
        let p = self.parameters;
        let ocv = self.curve.mean_voltage(initial[0], next[0])?;
        let vp = f64::midpoint(initial[1], next[1]) / p.polarization_capacitance_f;
        let vb = f64::midpoint(initial[2], next[2]) / p.bus_capacitance_f;
        let current = ((Dd::from_f64(ocv) - Dd::from_f64(vp) - Dd::from_f64(vb))
            / Dd::from_f64(p.series_resistance_ohm))
        .to_f64();
        let capacitor_change = |index: usize, capacitance: f64| {
            let a = Dd::from_f64(initial[index]);
            let b = Dd::from_f64(next[index]);
            ((b - a) * (a * Dd::from_f64(0.5) + b * Dd::from_f64(0.5)) / Dd::from_f64(capacitance))
                .to_f64()
        };
        let dt = Dd::from_f64(duration_s);
        let interval = BatteryInterval {
            chemical_energy_change_j: self.curve.energy_change_j(initial[0], next[0])?,
            polarization_energy_change_j: capacitor_change(1, p.polarization_capacitance_f),
            bus_energy_change_j: capacitor_change(2, p.bus_capacitance_f),
            terminal_work_j: (dt * Dd::from_f64(vb) * Dd::from_f64(terminal_current_a)).to_f64(),
            series_heat_j: (Dd::from_f64(p.series_resistance_ohm)
                * Dd::from_f64(current)
                * Dd::from_f64(current)
                * dt)
                .to_f64(),
            polarization_heat_j: (Dd::from_f64(vp) / Dd::from_f64(p.polarization_resistance_ohm)
                * Dd::from_f64(vp)
                * dt)
                .to_f64(),
            discrete_ocv_v: ocv,
            discrete_polarization_voltage_v: vp,
            discrete_terminal_voltage_v: vb,
            discrete_cell_current_a: current,
        };
        if ![
            interval.chemical_energy_change_j,
            interval.polarization_energy_change_j,
            interval.bus_energy_change_j,
            interval.terminal_work_j,
            interval.series_heat_j,
            interval.polarization_heat_j,
            ocv,
            vp,
            vb,
            current,
            interval.storage_change_j(),
            interval.heat_j(),
            interval.balance_residual_j(),
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(invalid("battery interval work or heat is unrepresentable"));
        }
        Ok(interval)
    }

    /// Prepare the existing three-state solver once. No new integration rule
    /// or allocation is introduced during a successful warm step.
    pub fn prepare_controlled(
        self,
        solve: SolveBudget,
        backtracking_halvings: u8,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedBatterySource, PreparedStepError> {
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(solve)?;
        workspace.set_backtracking_limit(backtracking_halvings)?;
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        Ok(PreparedBatterySource {
            model: self,
            workspace,
        })
    }
}

/// Reusable native solve scratch; all accepted chemical and capacitor charges
/// remain in caller-owned state. A refused/cancelled call publishes no state.
pub struct PreparedBatterySource {
    model: BatterySource,
    workspace: StepWorkspace,
}

impl PreparedBatterySource {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn parameters(&self) -> BatteryParameters {
        self.model.parameters()
    }
    pub fn curve(&self) -> &BatteryOcvCurve {
        self.model.curve()
    }

    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: f64,
        duration_s: f64,
        next: &mut [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryStep, PreparedStepError> {
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        let (mut candidate, mut voltage) = ([0.0; 3], [0.0]);
        let native = self.workspace.step_into_controlled(
            &self.model.system,
            initial,
            &[terminal_current_a],
            duration_s,
            &mut candidate,
            &mut voltage,
            &mut cancelled,
        )?;
        let battery = self
            .model
            .interval(initial, &candidate, terminal_current_a, duration_s)?;
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        *next = candidate;
        Ok(BatteryStep { native, battery })
    }
}
