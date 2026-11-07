//! Reciprocal permanent-magnet and synchronous-reluctance motor/generator.
//!
//! The five original states are `[lambda_alpha, lambda_beta, theta, p_r, p_s]`:
//! two stator-frame flux linkages, relative MECHANICAL rotor angle, and absolute
//! rotor/stator angular momenta. Electrical angle is `phi = pole_pairs*theta`.
//! With `lambda_dq = rotation(-phi)*lambda_alpha_beta`, magnetic storage is
//! `(lambda_d-psi_pm)^2/(2 Ld) + lambda_q^2/(2 Lq)`. Its angle derivative gives
//! `tau = p*(psi_pm*i_q + (Ld-Lq)*i_d*i_q)` on the rotor and its negative on the
//! stator. The SAME energy law therefore owns PM torque, saliency torque and
//! motional voltage; none is a separately fitted torque map.
//!
//! All alpha/beta AND d/q variables use POWER-INVARIANT normalization. For
//! positive phase order a,b,c, alpha is sqrt(2/3)*(a-(b+c)/2), beta is
//! (b-c)/sqrt(2), and d is aligned with the rotor PM axis. Positive q is a
//! counterclockwise quarter turn from d. Thus electrical power is v_d*i_d +
//! v_q*i_q, with NO 3/2 multiplier. A conventional phase-peak PM linkage must
//! be multiplied by sqrt(3/2) before supplying `pm_flux_linkage_wb`.
//!
//! The winding axes move with the stator; only relative shaft motion generates
//! speed voltage. Both finite inertias, equal/opposite bearing torque and all
//! signed ports remain in one native nonlinear port-Hamiltonian solve. The
//! source-owned symmetric discrete gradient preserves magnetic and kinetic
//! component work, with an analytic Newton tangent and no new integrator.
//!
//! [`SynchronousMotor::with_flux_curves`] instead admits bounded, reversible
//! independent d/q saturation laws. Their original energy and current own
//! `tau = p*(lambda_d*i_q-lambda_q*i_d)` and the same reciprocal speed voltage.
//!
//! This is a balanced three-wire, sinusoidal, frozen-parameter, coaxial machine.
//! Zero PM linkage admits synchronous reluctance; equal Ld/Lq removes saliency.
//! Cross-saturation, zero-sequence current, cogging, iron losses,
//! demagnetization, thermal evolution and device switching need their own laws.
//! Parameter admission is numerical, not a field or dynamometer validation.
//!
//! Reference: MIT 6.685, Notes 7, section 3.2, equations 1–19 (the reference uses
//! amplitude-invariant 3/2; this owner explicitly uses power invariance):
//! <https://ocw.mit.edu/courses/6-685-electric-machines-fall-2013/c5edfb0f28f7c1b1b4d1588f179074c6_MIT6_685F13_chapter7.pdf>.

use super::inductance::InductorFluxCurve;
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, ResistanceWitness,
    SolveBudget, StepWorkspace, Storage,
};
use fs_math::det;

mod storage;
use storage::SynchronousStorage;
mod saturation;
use saturation::{MagneticLaw, SaturatedStorage};

/// Finite DC bus, physical leg modulation and averaged channel conduction.
pub mod inverter;
/// Original-law fixed-d-current torque allocation under current/source limits.
pub mod reference;

/// Explicit frozen machine data, in SI units, with no inferred source card.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SynchronousMotorParameters {
    /// Linear d inductance, or source-derived zero-excess-flux reference when
    /// returned by a machine constructed with `with_flux_curves`.
    pub d_inductance_h: f64,
    /// Linear q inductance, or source-derived zero-flux reference when returned
    /// by a machine constructed with `with_flux_curves`.
    pub q_inductance_h: f64,
    /// Per-phase resistance, also the power-invariant alpha/beta resistance.
    pub winding_resistance_ohm: f64,
    /// POWER-INVARIANT d-axis PM flux linkage [Wb]. Signed magnet orientation;
    /// zero is a synchronous-reluctance machine. Not a phase-peak linkage.
    pub pm_flux_linkage_wb: f64,
    pub pole_pairs: u32,
    pub rotor_inertia_kg_m2: f64,
    pub stator_inertia_kg_m2: f64,
    /// Internal rotor-to-stator viscous drag [N m s/rad], not drag to ground.
    pub bearing_drag_nm_s_per_rad: f64,
}

/// Instantaneous observables from the SAME nonlinear storage, without advancing
/// any state. Electrical quantities below are power-invariant except abc current.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SynchronousOperatingPoint {
    pub current_alpha_beta_a: [f64; 2],
    pub current_dq_a: [f64; 2],
    /// Balanced a,b,c winding currents [A], with no zero-sequence channel.
    pub phase_current_a: [f64; 3],
    pub flux_dq_wb: [f64; 2],
    /// Local original-law differential d/q inductance [H], not secant L.
    /// At a source knot this follows the curve's one-sided convention: right
    /// segment at an interior knot, left segment at the final support endpoint.
    pub differential_inductance_dq_h: [f64; 2],
    pub shaft_speed_rad_s: [f64; 2],
    /// Electromagnetic torque ON the rotor; the stator receives its negative.
    pub electromagnetic_torque_nm: f64,
    /// Speed-voltage terms `[-omega_e*lambda_q, omega_e*lambda_d]` [V].
    /// Together with `R*i + d(lambda_dq)/dt`, these are the actual d/q voltages;
    /// they include PM back EMF and the reciprocal saliency exchange.
    pub speed_voltage_dq_v: [f64; 2],
}

/// Nonlinear native owner. Inputs are `[v_alpha, v_beta, tau_rotor, tau_stator]`;
/// outputs are `[i_alpha, i_beta, omega_rotor, omega_stator]`.
pub struct SynchronousMotor {
    parameters: SynchronousMotorParameters,
    magnetic_law: MagneticLaw,
    system: PortHamiltonian,
}

/// One accepted native discrete-gradient interval. Work and heat are its
/// second-order, energy-consistent DISCRETE quantities, not exact trajectory
/// integrals. Timestep refinement remains necessary for transient accuracy.
#[derive(Clone, Debug)]
pub struct SynchronousStep {
    pub native: PreparedStepRecord,
    /// Signed work entering the alpha and beta winding ports [J]. Negative
    /// total electrical work is regeneration, without a current sign clamp.
    pub electrical_work_j: [f64; 2],
    /// Signed external rotor/stator work [J], using absolute shaft speeds.
    pub shaft_work_j: [f64; 2],
    /// `dt*R*(i_alpha_bar^2+i_beta_bar^2)` from the accepted constitutive effort.
    pub copper_heat_j: f64,
    /// `dt*b*(omega_rotor_bar-omega_stator_bar)^2` from that SAME effort.
    /// Both heat entries describe native dissipation, never additional debits.
    pub bearing_heat_j: f64,
    /// Discrete electromagnetic impulse ON the rotor [N m s]; stator opposite.
    pub electromagnetic_angular_impulse_nm_s: f64,
    /// Discrete internal bearing impulse ON the rotor [N m s]; stator opposite.
    pub bearing_angular_impulse_nm_s: f64,
    /// Absolute rotor/stator angle increments [rad], consistent with quadratic
    /// kinetic storage. Their difference advances the relative-angle state.
    pub shaft_angle_increment_rad: [f64; 2],
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

/// Balanced three-wire winding voltage from physical phase potentials a,b,c.
/// Any common-mode potential is projected out because this model has no neutral
/// or zero-sequence current. Do not use this to discard an actual neutral path.
/// Phase power equals `v_alpha*i_alpha + v_beta*i_beta` for this owner's returned
/// balanced `phase_current_a`. No amplitude-invariant 3/2 factor is applied.
pub fn phase_voltage_to_alpha_beta(phase_voltage_v: [f64; 3]) -> Result<[f64; 2], PhsError> {
    let result = [
        det::sqrt(2.0 / 3.0)
            * (phase_voltage_v[0] - f64::midpoint(phase_voltage_v[1], phase_voltage_v[2])),
        det::sqrt(0.5) * (phase_voltage_v[1] - phase_voltage_v[2]),
    ];
    if phase_voltage_v
        .iter()
        .chain(&result)
        .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "synchronous phase voltage transform is unrepresentable",
        ));
    }
    Ok(result)
}

impl SynchronousMotor {
    /// Assemble one original five-state storage and constant reciprocal J/R/G.
    /// Angle and both momenta are retained: no imposed-speed or fixed-stator
    /// approximation hides torque reaction, acceleration or generated voltage.
    pub fn new(parameters: SynchronousMotorParameters) -> Result<Self, PhsError> {
        let storage = SynchronousStorage::new(parameters)?;
        Self::assemble(parameters, MagneticLaw::Linear, Box::new(storage))
    }

    /// Admit independent, bounded original current-versus-flux laws. The d
    /// curve's coordinate is `lambda_d - psi_pm`; the q curve's is `lambda_q`.
    /// Their stored energies replace the two quadratic magnetic terms while
    /// preserving the original five states and both finite shaft inertias.
    ///
    /// These curves are authoritative: the supplied `d_inductance_h` and
    /// `q_inductance_h` entries are replaced by their differential inductance
    /// at zero curve flux in returned parameter metadata. Use the operating
    /// point's actual differential inductances for a saturated state.
    ///
    /// Every source query, including Newton probes and the discrete gradient's
    /// four flux/angle corners, must remain within the original closed support.
    /// No extrapolation, fitted torque law, cross-saturation or hysteresis is
    /// inferred. Warm source queries and stepping allocate nothing.
    pub fn with_flux_curves(
        mut parameters: SynchronousMotorParameters,
        d: InductorFluxCurve,
        q: InductorFluxCurve,
    ) -> Result<Self, PhsError> {
        parameters.d_inductance_h = d
            .differential_inductance_h(0.0)
            .map_err(|_| invalid("synchronous d flux law has no zero-flux inductance"))?;
        parameters.q_inductance_h = q
            .differential_inductance_h(0.0)
            .map_err(|_| invalid("synchronous q flux law has no zero-flux inductance"))?;
        let storage = SaturatedStorage::new(parameters, d.clone(), q.clone())?;
        Self::assemble(parameters, MagneticLaw::Curves { d, q }, Box::new(storage))
    }

    fn assemble(
        parameters: SynchronousMotorParameters,
        magnetic_law: MagneticLaw,
        storage: Box<dyn Storage>,
    ) -> Result<Self, PhsError> {
        let r = parameters.winding_resistance_ohm;
        let b = parameters.bearing_drag_nm_s_per_rad;
        if [r, b].iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid(
                "synchronous winding and bearing losses must be finite nonnegative",
            ));
        }
        let mut j = vec![0.0; 25];
        j[13] = 1.0;
        j[14] = -1.0;
        j[17] = -1.0;
        j[22] = 1.0;
        let mut resistance = vec![0.0; 25];
        resistance[0] = r;
        resistance[6] = r;
        resistance[18] = b;
        resistance[19] = -b;
        resistance[23] = -b;
        resistance[24] = b;
        let witness = ResistanceWitness::try_from_source_rows(
            5,
            resistance.clone(),
            vec![
                det::sqrt(r),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                det::sqrt(r),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                det::sqrt(b),
                -det::sqrt(b),
            ],
        )?;
        let mut g = vec![0.0; 20];
        g[0] = 1.0;
        g[5] = 1.0;
        g[14] = 1.0;
        g[19] = 1.0;
        let system = PortHamiltonian::new(5, 4, j, resistance, g, storage)?
            .with_resistance_witness(witness)?;
        Ok(Self {
            parameters,
            magnetic_law,
            system,
        })
    }

    pub const fn parameters(&self) -> SynchronousMotorParameters {
        self.parameters
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Transfer all original states into native circuit/mechanical composition.
    /// For a DC converter bus and held ideal inverter vector m_alpha_beta,
    /// `interconnect_mapped(bus, motor, &[bus_port], &[0,1], &m, budget)` gives
    /// v_alpha_beta=m*Vbus and bus load current=−m.dot(i_alpha_beta), including
    /// regeneration. The composite nonlinear solve retains finite bus charge,
    /// armature flux, angle and both momenta; no lagged torque solve is added.
    /// Inverter gate validity/losses remain the actual converter owner's task.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Initialize power-invariant d/q current [A], relative MECHANICAL angle
    /// [rad], and absolute rotor/stator speeds [rad/s]. PM flux is retained even
    /// at zero current; resetting it to zero would create a different machine.
    pub fn state(
        &self,
        current_dq_a: [f64; 2],
        relative_angle_rad: f64,
        shaft_speed_rad_s: [f64; 2],
    ) -> Result<[f64; 5], PhsError> {
        let p = self.parameters;
        let (s, c) = det::sin_cos(f64::from(p.pole_pairs) * relative_angle_rad);
        let [d, q] = self.magnetic_law.flux_dq(p, current_dq_a)?;
        let state = [
            c.mul_add(d, -s * q),
            s.mul_add(d, c * q),
            relative_angle_rad,
            p.rotor_inertia_kg_m2 * shaft_speed_rad_s[0],
            p.stator_inertia_kg_m2 * shaft_speed_rad_s[1],
        ];
        if current_dq_a
            .iter()
            .chain(&shaft_speed_rad_s)
            .chain(&state)
            .any(|v| !v.is_finite())
            || !self.system.hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "synchronous initial current, angle, momentum or energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    /// Rotate power-invariant commanded d/q voltage into the stator axes using
    /// this machine's pole pairs and relative mechanical angle convention.
    /// Holding this returned alpha/beta voltage over a step does NOT continually
    /// rotate it as the rotor moves; a sampled controller must issue new values.
    pub fn voltage_from_dq(
        &self,
        relative_angle_rad: f64,
        voltage_dq_v: [f64; 2],
    ) -> Result<[f64; 2], PhsError> {
        let (s, c) = det::sin_cos(f64::from(self.parameters.pole_pairs) * relative_angle_rad);
        let voltage = [
            c.mul_add(voltage_dq_v[0], -s * voltage_dq_v[1]),
            s.mul_add(voltage_dq_v[0], c * voltage_dq_v[1]),
        ];
        if !relative_angle_rad.is_finite()
            || voltage_dq_v.iter().chain(&voltage).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "synchronous d/q voltage transform is unrepresentable",
            ));
        }
        Ok(voltage)
    }

    pub fn operating_point(&self, state: &[f64; 5]) -> Result<SynchronousOperatingPoint, PhsError> {
        storage::operating_point(self.parameters, &self.magnetic_law, state)
    }

    /// Allocate native five-state Newton scratch once, with explicit existing
    /// residual/update and 0..16 backtracking limits. The stored law supplies its
    /// own exact discrete-current/torque Jacobian; warm calls allocate nothing.
    /// Duration and held voltage/torque are supplied separately at each step.
    pub fn prepare_controlled(
        self,
        budget: SolveBudget,
        backtracking_halvings: u8,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedSynchronousMotor, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(budget)?;
        workspace.set_backtracking_limit(backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedSynchronousMotor {
            model: self,
            workspace,
        })
    }
}

/// Reusable native implicit solve; all accepted state/history stays with caller.
pub struct PreparedSynchronousMotor {
    model: SynchronousMotor,
    workspace: StepWorkspace,
}

impl PreparedSynchronousMotor {
    pub fn system(&self) -> &PortHamiltonian {
        &self.model.system
    }
    pub const fn parameters(&self) -> SynchronousMotorParameters {
        self.model.parameters
    }
    pub fn operating_point(&self, state: &[f64; 5]) -> Result<SynchronousOperatingPoint, PhsError> {
        self.model.operating_point(state)
    }

    /// Simultaneously advance winding flux, relative angle and both momenta.
    /// `inputs=[v_alpha,v_beta,tau_r,tau_s]` is held over the nonnegative interval.
    /// `outputs` receives the work-conjugate discrete currents and absolute
    /// shaft speeds. These are not endpoint currents or an exact-flow average.
    /// Native residual admission, finite physical reporting and final polling
    /// precede publication. A refused/cancelled trial is safe to retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 5],
        inputs: &[f64; 4],
        duration_s: f64,
        next: &mut [f64; 5],
        outputs: &mut [f64; 4],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SynchronousStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let (mut candidate, mut effort_ports) = ([0.0; 5], [0.0; 4]);
        let native = self.workspace.step_into_controlled(
            &self.model.system,
            initial,
            inputs,
            duration_s,
            &mut candidate,
            &mut effort_ports,
            &mut cancelled,
        )?;
        let effort = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("synchronous accepted constitutive effort is unavailable"))?;
        let p = self.model.parameters;
        let electrical_work_j = [
            duration_s * effort_ports[0] * inputs[0],
            duration_s * effort_ports[1] * inputs[1],
        ];
        let shaft_angle_increment_rad =
            [duration_s * effort_ports[2], duration_s * effort_ports[3]];
        let shaft_work_j = [
            shaft_angle_increment_rad[0] * inputs[2],
            shaft_angle_increment_rad[1] * inputs[3],
        ];
        let relative_speed = effort[3] - effort[4];
        let copper_heat_j = duration_s
            * (effort[0] * p.winding_resistance_ohm * effort[0]
                + effort[1] * p.winding_resistance_ohm * effort[1]);
        let bearing_heat_j =
            duration_s * relative_speed * p.bearing_drag_nm_s_per_rad * relative_speed;
        let electromagnetic_angular_impulse_nm_s = -duration_s * effort[2];
        let bearing_angular_impulse_nm_s =
            -duration_s * p.bearing_drag_nm_s_per_rad * relative_speed;
        if electrical_work_j
            .iter()
            .chain(&shaft_work_j)
            .chain(&shaft_angle_increment_rad)
            .chain(&[
                copper_heat_j,
                bearing_heat_j,
                electromagnetic_angular_impulse_nm_s,
                bearing_angular_impulse_nm_s,
            ])
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("synchronous work, heat or angular impulse is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        outputs.copy_from_slice(&effort_ports);
        Ok(SynchronousStep {
            native,
            electrical_work_j,
            shaft_work_j,
            copper_heat_j,
            bearing_heat_j,
            electromagnetic_angular_impulse_nm_s,
            bearing_angular_impulse_nm_s,
            shaft_angle_increment_rad,
        })
    }
}
