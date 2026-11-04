//! A non-inverting, bidirectional, averaged four-switch buck/boost power stage.
//!
//! The two synchronous half bridges connect one inductor to two capacitive
//! buses. With high-side connection fractions `a,b`, the original equations are
//! `L di/dt = a*v_in - b*v_out - R*i`, `C_in dv_in/dt = I_in - a*i`,
//! and `C_out dv_out/dt = I_out + b*i`. Both terminal currents enter the stage.
//! The SAME fractions appear in the voltage and current laws, so conversion
//! exchanges power without manufacturing it. Buck operation uses `(a,b)=(d,1)`;
//! boost uses `(1,1-d)`. Neither direction is diode-clamped at zero current.
//!
//! This is the forced synchronous, small-switching-ripple averaged lane, not a
//! transistor or PWM-event solver. Ideal complementary gates have no dead time,
//! shoot-through, gate power, switching loss, reverse recovery or diode mode.
//! Winding and channel resistances are frozen caller-admitted physical values.
//! Their losses describe the squared AVERAGED inductor current; unresolved PWM
//! ripple would add RMS loss. No efficiency multiplier or independent loss map
//! is applied. Temperature/source validity and the averaging regime belong to
//! the caller, not to an inferred semiconductor card.
//!
//! Topology and averaging limits: TI SLVA535B (2018),
//! <https://www.ti.com/lit/an/slva535b/slva535b.pdf>; conduction versus switching
//! losses: Analog Devices AN-140,
//! <https://www.analog.com/en/resources/app-notes/an-140.html>.

use crate::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use crate::{PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage};

/// Five distinct heat owners, in the same order as [`ConductionHeat::joules`].
pub const CONDUCTION_CHANNELS: usize = 5;

/// Frozen linear constitutive data, in SI units, for one physical power stage.
/// Numerical admission is not material provenance or a device safe-area claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BuckBoostParameters {
    pub inductance_h: f64,
    /// Input and output bus capacitance, in that order [F].
    pub bus_capacitance_f: [f64; 2],
    pub winding_resistance_ohm: f64,
    /// Input half-bridge high-side and low-side on-resistances [ohm].
    pub input_on_resistance_ohm: [f64; 2],
    /// Output half-bridge high-side and low-side on-resistances [ohm].
    pub output_on_resistance_ohm: [f64; 2],
}

/// Duty-averaged connection to each positive bus, held over a solve interval.
/// Each low-side fraction is its high-side complement; the two legs can have
/// different fractions. Fractions zero and one are valid static connections.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeModulation {
    pub input_high_fraction: f64,
    pub output_high_fraction: f64,
}

/// Independent physical conduction-loss allocation on one native trajectory.
/// No entry includes the energy delivered to a load or another storage element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionHeat {
    /// Winding, input high, input low, output high, output low [J]. All five
    /// paths carry the same averaged current; their effective resistances
    /// partition the original positive integral of `R*i(t)^2`.
    pub joules: [f64; CONDUCTION_CHANNELS],
}

/// The native record plus terminal work and its sole physical heat partition.
/// `native.dissipated` and `conduction` describe the SAME heat, never two debits.
#[derive(Clone, Debug)]
pub struct ConverterStep {
    pub native: PreparedStepRecord,
    /// Signed electrical work entering the input and output buses [J]. A load
    /// draws negative work; reverse conversion exchanges these signs naturally.
    pub terminal_work_j: [f64; 2],
    pub conduction: ConductionHeat,
}

/// Native three-state model: `[inductor flux, input charge, output charge]`.
/// Port inputs are currents INTO the two buses; outputs are their voltages.
/// The immutable original pHS can be composed with other physical owners.
pub struct AveragedBuckBoost {
    parameters: BuckBoostParameters,
    modulation: BridgeModulation,
    effective_resistance_ohm: [f64; CONDUCTION_CHANNELS],
    system: PortHamiltonian,
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl AveragedBuckBoost {
    /// Assemble the full original flux/charge model. All component energies,
    /// both source ports, and all five conduction paths are retained. No static
    /// conversion-ratio approximation replaces the inductor or either bus.
    ///
    /// A caller changing duty keeps the SAME flux/charges: only the skew
    /// exchange and duty-weighted resistance change. Storage is independent of
    /// modulation, so the ideal gate command does not introduce parameter work.
    /// Gate-drive energy is outside this model, not relabeled as zero physical
    /// semiconductor loss. Changing L or C at runtime needs its own energy law.
    ///
    /// # Errors
    /// Nonpositive/unrepresentable storage, negative/nonfinite resistance,
    /// out-of-range duty, overflow, or original skew/PSD admission refusal.
    pub fn new(
        parameters: BuckBoostParameters,
        modulation: BridgeModulation,
    ) -> Result<Self, PhsError> {
        let storage = [
            parameters.inductance_h,
            parameters.bus_capacitance_f[0],
            parameters.bus_capacitance_f[1],
        ];
        if storage
            .iter()
            .any(|&value| !positive(value) || !positive(1.0 / value))
        {
            return Err(invalid(
                "converter requires finite positive representable L and bus C",
            ));
        }
        let resistance = [
            parameters.winding_resistance_ohm,
            parameters.input_on_resistance_ohm[0],
            parameters.input_on_resistance_ohm[1],
            parameters.output_on_resistance_ohm[0],
            parameters.output_on_resistance_ohm[1],
        ];
        if resistance
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "converter conduction requires nonnegative finite resistances",
            ));
        }
        let a = modulation.input_high_fraction;
        let b = modulation.output_high_fraction;
        if ![a, b]
            .iter()
            .all(|value| value.is_finite() && (0.0..=1.0).contains(value))
        {
            return Err(invalid(
                "converter high-side connection fractions must lie in [0,1]",
            ));
        }
        let effective = [
            resistance[0],
            a * resistance[1],
            (1.0 - a) * resistance[2],
            b * resistance[3],
            (1.0 - b) * resistance[4],
        ];
        let total: f64 = effective.iter().sum();
        if !total.is_finite() {
            return Err(invalid(
                "converter total conduction resistance is unrepresentable",
            ));
        }
        let q = vec![
            1.0 / storage[0],
            0.0,
            0.0,
            0.0,
            1.0 / storage[1],
            0.0,
            0.0,
            0.0,
            1.0 / storage[2],
        ];
        let system = PortHamiltonian::new(
            3,
            2,
            vec![0.0, a, -b, -a, 0.0, 0.0, b, 0.0, 0.0],
            vec![total, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            Box::new(QuadraticStorage::new(q, 3)?),
        )?;
        Ok(Self {
            parameters,
            modulation,
            effective_resistance_ohm: effective,
            system,
        })
    }

    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn modulation(&self) -> BridgeModulation {
        self.modulation
    }
    pub const fn effective_resistance_ohm(&self) -> &[f64; CONDUCTION_CHANNELS] {
        &self.effective_resistance_ohm
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Transfer this SAME owner into native interconnection, conductance loading
    /// or multirate evolution. Further owners' losses must remain separate from
    /// this stage's five channels; a composite total is not converter heat.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Express physical initial current and bus voltages as original flux/charge
    /// coordinates. Signs are retained; no zero-current clipping or projection.
    pub fn state(&self, current_a: f64, bus_voltage_v: [f64; 2]) -> Result<[f64; 3], PhsError> {
        if [current_a, bus_voltage_v[0], bus_voltage_v[1]]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "converter initial current and bus voltages must be finite",
            ));
        }
        let state = [
            self.parameters.inductance_h * current_a,
            self.parameters.bus_capacitance_f[0] * bus_voltage_v[0],
            self.parameters.bus_capacitance_f[1] * bus_voltage_v[1],
        ];
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid(
                "converter initial flux, charge or energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    /// Prepare the original exact quadratic held-duty/held-current flow once.
    /// The fixed three-state preparation uses the existing dimension/work caps.
    /// Repeated successful steps perform no owner allocation or exponential.
    /// They include independently integrated conduction heat and bus voltages.
    /// Cancellation is polled around bounded native preparation and publication.
    pub fn prepare_controlled(
        self,
        duration_s: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAveragedBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        let total: f64 = self.effective_resistance_ohm.iter().sum();
        let fractions = self
            .effective_resistance_ohm
            .map(|value| if total > 0.0 { value / total } else { 0.0 });
        let native = PreparedExactZoh::new_owned(self.system, duration_s, budget)?;
        poll(&mut cancelled)?;
        Ok(PreparedAveragedBuckBoost { native, fractions })
    }
}

/// Fixed-duty native propagator. Caller state is published only after the same
/// original physical-energy gate, complete heat/work evaluation and final poll.
/// Both terminal currents and duty are held over the entire prepared duration.
pub struct PreparedAveragedBuckBoost {
    native: PreparedExactZoh<'static>,
    fractions: [f64; CONDUCTION_CHANNELS],
}

impl PreparedAveragedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        self.native.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }

    /// Advance all flux/charge coordinates on the native cached trajectory.
    /// The native total heat is partitioned using the actual duty-weighted
    /// resistances; `terminal_work_j[k] = I[k] * integral(v[k] dt)`.
    /// In particular, endpoint bus voltage is never used as interval work.
    /// No accepted state/history belongs to this cache. A refused/cancelled
    /// trial leaves both output buffers unchanged and is safe to retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConverterStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut candidate = [0.0; 3];
        let mut voltage = [0.0; 2];
        let native =
            self.native
                .step_into(initial, terminal_current_a, &mut candidate, &mut voltage)?;
        let terminal_work_j =
            core::array::from_fn(|k| (self.native.duration() * voltage[k]) * terminal_current_a[k]);
        let joules = self.fractions.map(|fraction| fraction * native.dissipated);
        if terminal_work_j
            .iter()
            .chain(&joules)
            .any(|value| !value.is_finite())
            || joules.iter().any(|value| *value < 0.0)
        {
            return Err(
                invalid("converter terminal work or conduction heat is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        average_bus_voltage_v.copy_from_slice(&voltage);
        Ok(ConverterStep {
            native,
            terminal_work_j,
            conduction: ConductionHeat { joules },
        })
    }
}
