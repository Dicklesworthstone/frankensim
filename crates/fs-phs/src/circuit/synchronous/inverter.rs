//! Finite DC capacitor and averaged complementary three-phase bridge feeding
//! the SAME native synchronous machine, including modulation changes.
//!
//! A high-side fraction d_j in [0,1] and complementary low-side fraction 1-d_j
//! select each physical leg. In the original power-invariant stator axes,
//! `v_alpha_beta = Clarke(d)*Vdc - R_bridge*i_alpha_beta`, while
//! `Qdc' = I_source - Clarke(d).i_alpha_beta`. The reciprocal current is signed
//! and is never reconstructed by dividing power by a possibly zero bus voltage.
//! Each phase has `r_j=d_j*Rhigh_j+(1-d_j)*Rlow_j`; its six positive channel
//! heats describe that SAME original resistance, not extra energy debits.
//!
//! The capacitor, winding fluxes, mechanical angle and both shaft momenta are
//! solved simultaneously by the existing component discrete-gradient owner.
//! Each held-duty interval refreshes the reciprocal structure and its checked
//! resistance witness in preallocated storage. No flux, charge or angle is
//! projected at a command edge. Duty changes do not change the Hamiltonian.
//!
//! This is the forced-synchronous small-ripple AVERAGED lane: squared averaged
//! phase current omits unresolved PWM ripple heating. Ideal complementary gates
//! supply no dead time, diode/recovery behavior, shoot-through, switching/gate
//! heat or device safe-area model. Resistance is frozen caller-supplied data.
//! Initial and accepted DC charge must be nonnegative; reverse-bus semiconductor
//! operation refuses instead of silently introducing diode physics or clipping.
//! Endpoint admission is not a continuous-time bus-positivity certificate.
//!
//! Physical bridge states and DC current signs: TI SBOA587, table 1,
//! <https://www.ti.com/document-viewer/lit/html/SBOA587/GUID-915C0393-540E-4DBA-B00B-2A6B52D87277>.

use super::{
    SynchronousMotor, SynchronousMotorParameters, SynchronousOperatingPoint, invalid, poll,
};
use crate::{
    InterconnectionBudget, PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord,
    QuadraticStorage, ResistanceWitness, SolveBudget, StepWorkspace, interconnect_mapped,
};
use fs_math::det;

/// Original finite bus and per-phase high/low channel resistances in SI units.
/// Winding resistance and machine data stay in the supplied SynchronousMotor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InverterParameters {
    pub dc_capacitance_f: f64,
    /// a,b,c high-side channel resistances [ohm].
    pub high_on_resistance_ohm: [f64; 3],
    /// a,b,c low-side channel resistances [ohm].
    pub low_on_resistance_ohm: [f64; 3],
}

/// Native discrete work over one held-duty/source interval. All current,
/// voltage, heat and impulse observations use the accepted constitutive effort;
/// they have the existing second-order time-discretization error.
#[derive(Clone, Debug)]
pub struct InverterStep {
    pub native: PreparedStepRecord,
    /// Signed external source-current work entering the finite DC capacitor [J].
    pub source_work_j: f64,
    /// Signed external rotor/stator work using absolute shaft speeds [J].
    pub shaft_work_j: [f64; 2],
    /// Signed work entering the machine winding terminals, AFTER bridge drops.
    pub machine_electrical_work_j: f64,
    /// Signed work leaving the DC capacitor into the bridge [J]; negative is
    /// regeneration into the finite bus. Includes physical channel conduction.
    pub bridge_work_j: f64,
    /// Balanced a,b,c discrete winding currents [A].
    pub phase_current_a: [f64; 3],
    /// Actual phase-to-floating-neutral terminal voltages [V], including each
    /// channel drop AND its neutral shift. Not ideal pole voltages d_j*Vdc.
    pub phase_voltage_v: [f64; 3],
    /// Signed current leaving the DC capacitor into the bridge [A].
    pub dc_bridge_current_a: f64,
    pub copper_heat_j: f64,
    pub bearing_heat_j: f64,
    /// a,b,c phase heat, each `[high-side, low-side]` [J]. These channels,
    /// winding copper and bearing heat partition the SAME native dissipation.
    pub channel_heat_j: [[f64; 2]; 3],
    /// Electromagnetic impulse ON rotor [N m s]; stator gets its negative.
    pub electromagnetic_angular_impulse_nm_s: f64,
    /// Internal bearing impulse ON rotor [N m s]; stator gets its negative.
    pub bearing_angular_impulse_nm_s: f64,
    pub shaft_angle_increment_rad: [f64; 2],
}

#[derive(Clone, Copy)]
struct BridgeStamp {
    high_fraction: [f64; 3],
    modulation: [f64; 2],
    channel_resistance: [[f64; 2]; 3],
    phase_resistance: [f64; 3],
    resistance: [f64; 36],
    root: [f64; 36],
}

impl BridgeStamp {
    fn new(
        p: InverterParameters,
        machine: SynchronousMotorParameters,
        base_resistance: &[f64; 36],
        scales: [f64; 2],
        high_fraction: [f64; 3],
    ) -> Result<Self, PhsError> {
        if high_fraction
            .iter()
            .any(|d| !d.is_finite() || !(0.0..=1.0).contains(d))
        {
            return Err(invalid(
                "inverter requires three finite physical high-side fractions in [0,1]",
            ));
        }
        let [a, b] = scales; // sqrt(2/3), sqrt(1/2), original power-invariant Clarke map
        let modulation = [
            a * (high_fraction[0] - f64::midpoint(high_fraction[1], high_fraction[2])),
            b * (high_fraction[1] - high_fraction[2]),
        ];
        let channel_resistance: [[f64; 2]; 3] = std::array::from_fn(|j| {
            [
                high_fraction[j] * p.high_on_resistance_ohm[j],
                (1.0 - high_fraction[j]) * p.low_on_resistance_ohm[j],
            ]
        });
        let phase_resistance = channel_resistance.map(|r| r[0] + r[1]);
        let [ra, rb, rc] = phase_resistance;
        let mut resistance = *base_resistance;
        // C diag(r_phase) C^T on the two ORIGINAL winding-flux coordinates.
        // Form the small off-diagonal from the source difference directly.
        resistance[7] += a * a * ra + (0.25 * a * a) * rb + (0.25 * a * a) * rc;
        resistance[14] += b * b * rb + b * b * rc;
        resistance[8] = (0.5 * a * b) * (rc - rb);
        resistance[13] = resistance[8];
        let mut root = [0.0; 36];
        let winding = machine.winding_resistance_ohm;
        let scale = phase_resistance.iter().fold(winding, |s, r| s.max(*r));
        if scale > 0.0 {
            // Positive source determinant of the ORIGINAL 2x2 electrical R.
            // This avoids subtracting Raa*Rbb-Rab^2 and preserves exact rank-one
            // loss. Normalization avoids squaring the dimensional resistance.
            let [x, y, z] = phase_resistance.map(|r| r / scale);
            let w = winding / scale;
            let determinant = w * w
                + w * (a * a * x + (0.25 * a * a + b * b) * (y + z))
                + (a * a * b * b) * (x * y + x * z + y * z);
            if resistance[7] <= 0.0 || !determinant.is_finite() {
                return Err(invalid(
                    "inverter electrical resistance factor is unrepresentable",
                ));
            }
            root[1] = det::sqrt(resistance[7]);
            root[2] = resistance[8] / root[1];
            root[8] = det::sqrt(scale) * det::sqrt(determinant / (resistance[7] / scale));
        }
        root[16] = det::sqrt(machine.bearing_drag_nm_s_per_rad);
        root[17] = -root[16];
        if resistance
            .iter()
            .chain(&root)
            .chain(&modulation)
            .chain(&phase_resistance)
            .chain(channel_resistance.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "inverter held structure or physical losses are unrepresentable",
            ));
        }
        crate::require_symmetric(&resistance, 6, "inverter original resistance")?;
        // The SAME native componentwise admission used by every source witness.
        // No replaced matrix, eigenvalue clipping, small-source cutoff or QR.
        crate::linear_exact::positive_heat::check_resistance_root(&resistance, &root, 6)?;
        Ok(Self {
            high_fraction,
            modulation,
            channel_resistance,
            phase_resistance,
            resistance,
            root,
        })
    }
}

/// Native six-state owner `[Qdc, lambda_alpha, lambda_beta, theta, p_r, p_s]`.
/// External inputs are `[I_source, tau_r, tau_s]`, outputs `[Vdc, omega_r, omega_s]`.
/// Construction selects the physical all-low zero vector. Modulation changes
/// are commands supplied to the prepared owner, never new stored coordinates.
pub struct AveragedSynchronousInverter {
    parameters: InverterParameters,
    machine_parameters: SynchronousMotorParameters,
    system: PortHamiltonian,
    base_resistance: [f64; 36],
    scales: [f64; 2],
    held: BridgeStamp,
}

impl AveragedSynchronousInverter {
    pub fn new(
        machine: SynchronousMotor,
        parameters: InverterParameters,
    ) -> Result<Self, PhsError> {
        let c = parameters.dc_capacitance_f;
        if !c.is_finite()
            || c <= 0.0
            || !(1.0 / c).is_finite()
            || 1.0 / c <= 0.0
            || parameters
                .high_on_resistance_ohm
                .iter()
                .chain(&parameters.low_on_resistance_ohm)
                .any(|r| !r.is_finite() || *r < 0.0)
        {
            return Err(invalid(
                "inverter requires positive representable DC capacitance and finite nonnegative channel resistances",
            ));
        }
        let machine_parameters = machine.parameters();
        // Two current ports share ONE capacitor: source current stays external;
        // the second closes reciprocally through the actual inverter vector.
        let capacitor = PortHamiltonian::new(
            1,
            2,
            vec![0.0],
            vec![0.0],
            vec![1.0, 1.0],
            Box::new(QuadraticStorage::new(vec![1.0 / c], 1)?),
        )?;
        let mut system = interconnect_mapped(
            capacitor,
            machine.into_system(),
            &[1],
            &[0, 1],
            &[0.0; 2],
            InterconnectionBudget {
                max_states: 6,
                max_total_ports: 6,
                max_dense_bytes: 1 << 16,
                max_coupling_multiply_terms: 128,
            },
        )?;
        let mut base_resistance = [0.0; 36];
        base_resistance.copy_from_slice(&system.r);
        let scales = [det::sqrt(2.0 / 3.0), det::sqrt(0.5)];
        let held = BridgeStamp::new(
            parameters,
            machine_parameters,
            &base_resistance,
            scales,
            [0.0; 3],
        )?;
        crate::require_psd(&held.resistance, 6, "inverter original resistance")?;
        system.r.copy_from_slice(&held.resistance);
        let witness = ResistanceWitness::try_from_source_rows(
            6,
            held.resistance.to_vec(),
            held.root.to_vec(),
        )?;
        system = system.with_resistance_witness(witness)?;
        Ok(Self {
            parameters,
            machine_parameters,
            system,
            base_resistance,
            scales,
            held,
        })
    }

    pub const fn parameters(&self) -> InverterParameters {
        self.parameters
    }
    pub const fn machine_parameters(&self) -> SynchronousMotorParameters {
        self.machine_parameters
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Attach a nonnegative finite bus voltage [V] to the SAME five machine
    /// coordinates previously initialized by `SynchronousMotor::state`.
    pub fn state(
        &self,
        bus_voltage_v: f64,
        machine_state: &[f64; 5],
    ) -> Result<[f64; 6], PhsError> {
        let mut state = [0.0; 6];
        state[0] = self.parameters.dc_capacitance_f * bus_voltage_v;
        state[1..].copy_from_slice(machine_state);
        if !bus_voltage_v.is_finite()
            || bus_voltage_v < 0.0
            || state.iter().any(|v| !v.is_finite())
            || !self.system.hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "inverter initial nonnegative bus charge or machine energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    fn admit_bus(&self, state: &[f64; 6]) -> Result<(), PhsError> {
        if !state[0].is_finite()
            || state[0] < 0.0
            || !(state[0] / self.parameters.dc_capacitance_f).is_finite()
        {
            Err(invalid(
                "inverter requires a finite nonnegative DC bus charge and voltage",
            ))
        } else {
            Ok(())
        }
    }

    fn apply_stamp(&mut self, stamp: &BridgeStamp) {
        // Private constructor guarantees these fixed dimensions and the witness.
        // Every stamp passed here has already passed the original native gate.
        self.system.j[1] = -stamp.modulation[0];
        self.system.j[2] = -stamp.modulation[1];
        self.system.j[6] = stamp.modulation[0];
        self.system.j[12] = stamp.modulation[1];
        self.system.r.copy_from_slice(&stamp.resistance);
        let witness = self
            .system
            .resistance_witness
            .as_mut()
            .expect("inverter retains its source resistance witness");
        witness.matrix.copy_from_slice(&stamp.resistance);
        witness.root.copy_from_slice(&stamp.root);
    }

    /// Prepare native six-state scratch once. Existing explicit residual/update
    /// and 0..16 backtracking limits apply unchanged at every modulation sample.
    pub fn prepare_controlled(
        self,
        budget: SolveBudget,
        backtracking_halvings: u8,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedSynchronousInverter, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(budget)?;
        workspace.set_backtracking_limit(backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedSynchronousInverter {
            model: self,
            workspace,
        })
    }
}

/// Prepared averaged bridge. Accepted state belongs to caller; the last accepted
/// held modulation and its matching native matrices/witness belong to this cache.
/// No allocation, interconnection rebuild or exponential occurs on warm steps.
pub struct PreparedSynchronousInverter {
    model: AveragedSynchronousInverter,
    workspace: StepWorkspace,
}

impl PreparedSynchronousInverter {
    pub fn system(&self) -> &PortHamiltonian {
        &self.model.system
    }
    pub const fn parameters(&self) -> InverterParameters {
        self.model.parameters
    }
    pub const fn machine_parameters(&self) -> SynchronousMotorParameters {
        self.model.machine_parameters
    }
    pub const fn high_fraction(&self) -> [f64; 3] {
        self.model.held.high_fraction
    }

    pub fn machine_operating_point(
        &self,
        state: &[f64; 6],
    ) -> Result<SynchronousOperatingPoint, PhsError> {
        self.model.admit_bus(state)?;
        super::storage::operating_point(
            self.model.machine_parameters,
            &[state[1], state[2], state[3], state[4], state[5]],
        )
    }

    /// Hold actual high-side fractions and external source current/shaft torques
    /// for this interval. A common duty shift leaves the ideal voltage/current
    /// map unchanged. It can change REAL heat when high/low resistances differ.
    ///
    /// The native workspace rereads J/R, refreshes its complete dependency rows
    /// and refactors the actual analytic Jacobian at every call. The witness is
    /// updated with the SAME new R; no prepared exact-flow cache is reused.
    /// Negative accepted bus charge, solver/report refusal or final cancellation
    /// restores the previous J/R/witness and publishes neither output buffer.
    #[allow(clippy::too_many_arguments)]
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 6],
        high_fraction: [f64; 3],
        inputs: &[f64; 3],
        duration_s: f64,
        next: &mut [f64; 6],
        outputs: &mut [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<InverterStep, PreparedStepError> {
        poll(&mut cancelled)?;
        self.model.admit_bus(initial)?;
        let stamp = BridgeStamp::new(
            self.model.parameters,
            self.model.machine_parameters,
            &self.model.base_resistance,
            self.model.scales,
            high_fraction,
        )?;
        poll(&mut cancelled)?;
        self.model.apply_stamp(&stamp);
        let (mut candidate, mut output) = ([0.0; 6], [0.0; 3]);
        let result: Result<InverterStep, PreparedStepError> = (|| {
            let native = self.workspace.step_into_controlled(
                &self.model.system,
                initial,
                inputs,
                duration_s,
                &mut candidate,
                &mut output,
                &mut cancelled,
            )?;
            self.model.admit_bus(&candidate)?;
            let e = self
                .workspace
                .last_step_discrete_effort()
                .ok_or_else(|| invalid("inverter accepted constitutive effort is unavailable"))?;
            let phase_a = self.model.scales[0] * e[1];
            let quadrature = self.model.scales[1] * e[2];
            let phase_current_a = [
                phase_a,
                -0.5 * phase_a + quadrature,
                -0.5 * phase_a - quadrature,
            ];
            let dc_bridge_current_a = stamp.modulation[0] * e[1] + stamp.modulation[1] * e[2];
            let drop: [f64; 3] =
                std::array::from_fn(|j| stamp.phase_resistance[j] * phase_current_a[j]);
            let mean_drop = drop.iter().map(|v| v / 3.0).sum::<f64>();
            let mean_duty = high_fraction.iter().map(|d| d / 3.0).sum::<f64>();
            let phase_voltage_v: [f64; 3] = std::array::from_fn(|j| {
                (high_fraction[j] - mean_duty) * e[0] - drop[j] + mean_drop
            });
            let machine_electrical_work_j = duration_s
                * phase_voltage_v
                    .iter()
                    .zip(phase_current_a)
                    .map(|(v, i)| v * i)
                    .sum::<f64>();
            let bridge_work_j = duration_s * e[0] * dc_bridge_current_a;
            let channel_heat_j: [[f64; 2]; 3] = std::array::from_fn(|j| {
                std::array::from_fn(|k| {
                    duration_s
                        * phase_current_a[j]
                        * stamp.channel_resistance[j][k]
                        * phase_current_a[j]
                })
            });
            let p = self.model.machine_parameters;
            let relative_speed = e[4] - e[5];
            let copper_heat_j = duration_s
                * (e[1] * p.winding_resistance_ohm * e[1] + e[2] * p.winding_resistance_ohm * e[2]);
            let bearing_heat_j =
                duration_s * relative_speed * p.bearing_drag_nm_s_per_rad * relative_speed;
            let electromagnetic_angular_impulse_nm_s = -duration_s * e[3];
            let bearing_angular_impulse_nm_s =
                -duration_s * p.bearing_drag_nm_s_per_rad * relative_speed;
            let source_work_j = duration_s * inputs[0] * output[0];
            let shaft_angle_increment_rad = [duration_s * output[1], duration_s * output[2]];
            let shaft_work_j = [
                shaft_angle_increment_rad[0] * inputs[1],
                shaft_angle_increment_rad[1] * inputs[2],
            ];
            if [
                source_work_j,
                machine_electrical_work_j,
                bridge_work_j,
                dc_bridge_current_a,
                copper_heat_j,
                bearing_heat_j,
                electromagnetic_angular_impulse_nm_s,
                bearing_angular_impulse_nm_s,
            ]
            .iter()
            .chain(&shaft_work_j)
            .chain(&shaft_angle_increment_rad)
            .chain(&phase_current_a)
            .chain(&phase_voltage_v)
            .chain(channel_heat_j.iter().flatten())
            .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "inverter signed work, current or physical heat is unrepresentable",
                )
                .into());
            }
            poll(&mut cancelled)?;
            Ok(InverterStep {
                native,
                source_work_j,
                shaft_work_j,
                machine_electrical_work_j,
                bridge_work_j,
                phase_current_a,
                phase_voltage_v,
                dc_bridge_current_a,
                copper_heat_j,
                bearing_heat_j,
                channel_heat_j,
                electromagnetic_angular_impulse_nm_s,
                bearing_angular_impulse_nm_s,
                shaft_angle_increment_rad,
            })
        })();
        match result {
            Ok(record) => {
                self.model.held = stamp;
                next.copy_from_slice(&candidate);
                outputs.copy_from_slice(&output);
                Ok(record)
            }
            Err(error) => {
                let previous = self.model.held;
                self.model.apply_stamp(&previous);
                Err(error)
            }
        }
    }
}
