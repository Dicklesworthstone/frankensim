//! Finite chemical energy, DC/DC conversion and a reciprocal voltage motor.
//!
//! The original six electrical states precede the motor's flux and two absolute
//! momenta. A native Dirac connection sends the actual output-bus voltage into
//! the winding and draws the SAME signed armature current from that bus. Battery
//! sag, back EMF, acceleration, auxiliary demand and regeneration therefore
//! participate in one nine-state solve. No requested shaft-power map or held
//! interface voltage substitutes for the electrical/mechanical trajectory.
//!
//! Converter modulation and constitutive parameters remain frozen. Commands
//! are external rotor/stator torques and optional regulated auxiliary demand.
//! Both inertias and the equal/opposite electromagnetic and bearing reaction
//! remain physical. A supported stator requires an authored reaction; it is not
//! silently pinned. This is the existing averaged converter and commutated
//! linear-winding motor, without new switching, magnetic or control models.

use super::{
    BatteryFedBuckBoost, BatteryFedInterval, PowerLoadControls, RegulatedPowerLoad, admit_balances,
    invalid, load_point, poll, positive, product, sum,
};
use crate::circuit::motor::{DcMotor, DcMotorParameters};
use crate::{
    InterconnectionBudget, NonlinearPortPower, PhsError, PortHamiltonian, PreparedStepError,
    PreparedStepRecord, StepWorkspace, interconnect_mapped,
};
use fs_math::dd::Dd;

/// One actual nine-state native solve and its physical constituent accounts.
/// No separate six-state or three-state solver certificate is manufactured.
#[derive(Clone, Debug)]
pub struct BatteryMotorStep {
    pub native: PreparedStepRecord,
    pub electrical: BatteryFedInterval,
    pub magnetic_energy_change_j: f64,
    pub kinetic_energy_change_j: [f64; 2],
    /// Signed work entering the motor winding. Negative is generation; this
    /// cancels `electrical.external_work_j` inside the complete drive boundary.
    pub motor_electrical_work_j: f64,
    /// Signed external work into rotor and stator, using absolute shaft speeds.
    pub shaft_work_j: [f64; 2],
    pub copper_heat_j: f64,
    pub bearing_heat_j: f64,
    pub discrete_motor_current_a: f64,
    pub discrete_shaft_speed_rad_s: [f64; 2],
    /// Torque on the rotor; the stator receives its exact opposite.
    pub discrete_electromagnetic_torque_nm: f64,
    pub discrete_bearing_torque_nm: f64,
    pub discrete_back_emf_v: f64,
    /// Momentum equation diagnostic, in its own units, not an energy correction.
    pub angular_momentum_residual_nms: f64,
}

impl BatteryMotorStep {
    /// Existing eight electrical heat channels, then motor copper and bearing.
    /// Auxiliary-load export and external shaft work are never routed as heat.
    pub fn heat_joules(&self) -> [f64; 10] {
        let e = self.electrical.heat_joules();
        [
            e[0],
            e[1],
            e[2],
            e[3],
            e[4],
            e[5],
            e[6],
            e[7],
            self.copper_heat_j,
            self.bearing_heat_j,
        ]
    }
    pub fn motor_storage_change_j(&self) -> f64 {
        sum([
            self.magnetic_energy_change_j,
            self.kinetic_energy_change_j[0],
            self.kinetic_energy_change_j[1],
        ])
    }
    fn storage_changes(&self) -> [f64; 7] {
        let b = self.electrical.battery;
        [
            b.chemical_energy_change_j,
            b.polarization_energy_change_j,
            b.bus_energy_change_j,
            self.electrical.converter_storage_change_j,
            self.magnetic_energy_change_j,
            self.kinetic_energy_change_j[0],
            self.kinetic_energy_change_j[1],
        ]
    }
    pub fn balance_residual_j(&self) -> f64 {
        sum(self
            .storage_changes()
            .into_iter()
            .chain(self.heat_joules())
            .chain([
                -self.shaft_work_j[0],
                -self.shaft_work_j[1],
                self.electrical.load_work_j,
            ]))
    }
}

/// State: `[battery3, converter3, motor_flux, rotor_momentum, stator_momentum]`.
/// Source/storage and both finite inertias remain the original owners' laws.
pub struct BatteryMotorDrive {
    electrical: BatteryFedBuckBoost,
    motor_parameters: DcMotorParameters,
    system: PortHamiltonian,
}

impl BatteryFedBuckBoost {
    /// Close output port 2 against the existing motor voltage port 0. The
    /// reciprocal map has V_motor=V_output and I_into_output=-I_motor.
    pub fn with_dc_motor(mut self, motor: DcMotor) -> Result<BatteryMotorDrive, PhsError> {
        let motor_parameters = motor.parameters();
        let system = interconnect_mapped(
            self.system.share_storage(),
            motor.into_system(),
            &[2],
            &[0],
            &[1.0],
            InterconnectionBudget {
                max_states: 9,
                max_total_ports: 6,
                max_dense_bytes: 1 << 19,
                max_coupling_multiply_terms: 256,
            },
        )?;
        Ok(BatteryMotorDrive {
            electrical: self,
            motor_parameters,
            system,
        })
    }
}

impl BatteryMotorDrive {
    /// Remaining native ports: battery/source-input currents, rotor torque,
    /// stator torque. The prepared finite-source path holds both currents zero.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }
    pub fn electrical(&self) -> &BatteryFedBuckBoost {
        &self.electrical
    }
    pub fn motor_parameters(&self) -> DcMotorParameters {
        self.motor_parameters
    }

    /// Concatenate coordinates initialized by the two ORIGINAL owners. No
    /// current, angle, momentum, capacitor voltage or stored energy is reset.
    pub fn state(
        &self,
        electrical_state: &[f64; 6],
        motor_state: &[f64; 3],
    ) -> Result<[f64; 9], PhsError> {
        let mut state = [0.0; 9];
        state[..6].copy_from_slice(electrical_state);
        state[6..].copy_from_slice(motor_state);
        if !state.iter().all(|v| v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid(
                "battery-motor drive requires finite supported original states",
            ));
        }
        Ok(state)
    }

    pub fn prepare_controlled(
        self,
        duration_s: f64,
        controls: PowerLoadControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedBatteryMotorDrive, PreparedStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid(
                "battery-motor drive requires positive duration and finite energy allowances",
            )
            .into());
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedBatteryMotorDrive {
            model: self,
            workspace,
            duration_s,
            controls,
        })
    }
}

/// Nine-state fixed scratch. A failed physical account or cancellation leaves
/// the entire caller-owned state untouched; successful warm calls allocate none.
pub struct PreparedBatteryMotorDrive {
    model: BatteryMotorDrive,
    workspace: StepWorkspace,
    duration_s: f64,
    controls: PowerLoadControls,
}

fn quadratic_change(before: f64, after: f64, storage: f64) -> f64 {
    let a = Dd::from_f64(before);
    let b = Dd::from_f64(after);
    ((b - a) * (Dd::from_f64(0.5) * a + Dd::from_f64(0.5) * b) / Dd::from_f64(storage)).to_f64()
}

impl PreparedBatteryMotorDrive {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn electrical(&self) -> &BatteryFedBuckBoost {
        self.model.electrical()
    }
    pub fn motor_parameters(&self) -> DcMotorParameters {
        self.model.motor_parameters()
    }
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Joint physical interval with held external torques and auxiliary load.
    /// Regeneration is signed motor current through the same native connection,
    /// with no efficiency multiplier, power command or current clipping.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 9],
        shaft_torque_nm: &[f64; 2],
        auxiliary_load: Option<RegulatedPowerLoad>,
        next: &mut [f64; 9],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryMotorStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let c_out = self.model.electrical.parameters().bus_capacitance_f[1];
        load_point(auxiliary_load, initial[5] / c_out)?;
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if effort.len() != 9 || out.len() != 9 {
                return None;
            }
            let point = load_point(auxiliary_load, effort[5]).ok()?;
            out.fill(0.0);
            out[5] = point.current_a;
            Some(NonlinearPortPower {
                supplied_power_w: -point.delivered_power_w,
                dissipated_power_w: 0.0,
            })
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if effort.len() != 9 || de.len() != 9 || out.len() != 9 {
                return false;
            }
            let Ok(point) = load_point(auxiliary_load, effort[5]) else {
                return false;
            };
            out.fill(0.0);
            out[5] = point.incremental_conductance_s * de[5];
            out[5].is_finite()
        };
        let dt = self.duration_s;
        let (mut candidate, mut outputs) = ([0.0; 9], [0.0; 4]);
        let native = self.workspace.step_into_supplied_controlled(
            &self.model.system,
            initial,
            &[0.0, 0.0, shaft_torque_nm[0], shaft_torque_nm[1]],
            dt,
            &mut candidate,
            &mut outputs,
            &port,
            None,
            Some(&tangent),
            &mut cancelled,
        )?;
        let e = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("battery-motor drive requires the accepted native effort"))?;
        let electrical = self.model.electrical.account_interval(
            &core::array::from_fn(|i| initial[i]),
            &core::array::from_fn(|i| candidate[i]),
            &core::array::from_fn(|i| e[i]),
            -e[6],
            auxiliary_load,
            dt,
        )?;
        electrical.admit(self.controls)?;
        let p = self.model.motor_parameters;
        let current = e[6];
        let speed = [e[7], e[8]];
        let relative_speed = sum([speed[0], -speed[1]]);
        let electromagnetic_torque = p.coupling_nm_per_a * current;
        let bearing_torque = -p.bearing_drag_nm_s_per_rad * relative_speed;
        let emf = p.coupling_nm_per_a * relative_speed;
        let step = BatteryMotorStep {
            native,
            electrical,
            magnetic_energy_change_j: quadratic_change(initial[6], candidate[6], p.inductance_h),
            kinetic_energy_change_j: [
                quadratic_change(initial[7], candidate[7], p.rotor_inertia_kg_m2),
                quadratic_change(initial[8], candidate[8], p.stator_inertia_kg_m2),
            ],
            motor_electrical_work_j: product(dt, e[5], current),
            shaft_work_j: [
                product(dt, shaft_torque_nm[0], speed[0]),
                product(dt, shaft_torque_nm[1], speed[1]),
            ],
            copper_heat_j: (Dd::from_f64(p.winding_resistance_ohm)
                * Dd::from_f64(current)
                * Dd::from_f64(current)
                * Dd::from_f64(dt))
            .to_f64(),
            bearing_heat_j: (Dd::from_f64(p.bearing_drag_nm_s_per_rad)
                * Dd::from_f64(relative_speed)
                * Dd::from_f64(relative_speed)
                * Dd::from_f64(dt))
            .to_f64(),
            discrete_motor_current_a: current,
            discrete_shaft_speed_rad_s: speed,
            discrete_electromagnetic_torque_nm: electromagnetic_torque,
            discrete_bearing_torque_nm: bearing_torque,
            discrete_back_emf_v: emf,
            angular_momentum_residual_nms: (Dd::from_f64(candidate[7]) - Dd::from_f64(initial[7])
                + Dd::from_f64(candidate[8])
                - Dd::from_f64(initial[8])
                - Dd::from_f64(dt)
                    * (Dd::from_f64(shaft_torque_nm[0]) + Dd::from_f64(shaft_torque_nm[1])))
            .to_f64(),
        };
        let em_work = product(dt, emf, current);
        let magnetic_scale = sum([
            step.magnetic_energy_change_j.abs(),
            step.copper_heat_j,
            em_work.abs(),
            step.motor_electrical_work_j.abs(),
        ]);
        let mechanical_scale = sum(step
            .kinetic_energy_change_j
            .map(f64::abs)
            .into_iter()
            .chain(step.shaft_work_j.map(f64::abs))
            .chain([step.bearing_heat_j, em_work.abs()]));
        let shaft_work = sum(step.shaft_work_j);
        let heat = sum(step.heat_joules());
        let stored = sum(step.storage_changes());
        let [source_scale, converter_scale, _, output_scale] = step.electrical.scales();
        let scale = source_scale
            .max(converter_scale)
            .max(output_scale)
            .max(magnetic_scale)
            .max(mechanical_scale)
            .max(heat);
        let values = [
            step.motor_electrical_work_j,
            step.copper_heat_j,
            step.bearing_heat_j,
            current,
            speed[0],
            speed[1],
            electromagnetic_torque,
            bearing_torque,
            emf,
            step.angular_momentum_residual_nms,
            stored,
            heat,
        ];
        if values
            .iter()
            .chain(&step.storage_changes())
            .chain(&step.shaft_work_j)
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("battery-motor drive work, heat or momentum is unrepresentable").into(),
            );
        }
        admit_balances(
            self.controls,
            &[
                (
                    sum([
                        step.electrical.external_work_j,
                        step.motor_electrical_work_j,
                    ]),
                    sum([
                        step.electrical.external_work_j.abs(),
                        step.motor_electrical_work_j.abs(),
                    ]),
                ),
                (
                    sum([
                        step.magnetic_energy_change_j,
                        step.copper_heat_j,
                        em_work,
                        -step.motor_electrical_work_j,
                    ]),
                    magnetic_scale,
                ),
                (
                    sum([
                        step.kinetic_energy_change_j[0],
                        step.kinetic_energy_change_j[1],
                        step.bearing_heat_j,
                        -em_work,
                        -step.shaft_work_j[0],
                        -step.shaft_work_j[1],
                    ]),
                    mechanical_scale,
                ),
                (sum([stored, -native.delta_h]), scale),
                (sum([heat, -native.dissipated]), heat),
                (
                    sum([native.supplied, -shaft_work, step.electrical.load_work_j]),
                    sum(step
                        .shaft_work_j
                        .map(f64::abs)
                        .into_iter()
                        .chain([step.electrical.load_work_j])),
                ),
                (step.balance_residual_j(), scale),
            ],
            "battery, converter, winding and shafts do not close on the same trajectory",
        )?;
        poll(&mut cancelled)?;
        *next = candidate;
        Ok(step)
    }
}
