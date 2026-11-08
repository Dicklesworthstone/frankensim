//! A finite mechanical load driven through an elastic, damped shaft.
//!
//! The original nine-state battery/inverter/machine and the shaft twist and load
//! momentum evolve in ONE native solve. With twist `z = theta_rotor-theta_load`
//! relative to the stress-free phase, `z' = omega_rotor-omega_load` and the
//! transmitted torque is `T = k*z + c*z'`. The rotor receives `-T`, the load `T`;
//! spring storage is `k*z*z/2` and physical shaft heat is `c*z'*z'`.
//!
//! Load motion therefore changes rotor back EMF, winding current, DC bus sag
//! and battery charge during the same interval. There is no sampled load torque
//! or separately advanced mechanical subsystem. All three external torques use
//! their absolute shaft speeds; stator recoil remains the original machine's.
//!
//! This is a coaxial lumped shaft with fixed coefficients and finite inertia,
//! not a rigid constraint, gearbox, distributed shaft, backlash or fracture law.
//! The original averaged bridge and battery support bounds remain unchanged.

use super::{BatteryInverterControls, BatterySynchronousInverter, DriveStamp, sum};
use crate::circuit::battery::BatteryInterval;
use crate::circuit::synchronous::inverter::{BridgeStamp, InverterInterval};
use crate::circuit::synchronous::{SynchronousOperatingPoint, invalid, poll, reference, storage};
use crate::{
    InterconnectionBudget, PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord,
    QuadraticStorage, ResistanceWitness, StepWorkspace, interconnect_mapped, storage_energy_change,
};
use fs_math::{dd::Dd, det};

/// Original coaxial load inertia and passive shaft coefficients in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticShaftParameters {
    pub load_inertia_kg_m2: f64,
    pub stiffness_nm_per_rad: f64,
    pub damping_nm_s_per_rad: f64,
}

/// One complete eleven-state solve and separate physical component accounts.
#[derive(Clone, Debug)]
pub struct ElasticShaftStep {
    pub native: PreparedStepRecord,
    pub battery: BatteryInterval,
    /// The rotor work here includes the external torque MINUS shaft reaction.
    /// It is the actual inverter constituent boundary, not total external work.
    pub inverter: InverterInterval,
    pub inverter_storage_change_j: f64,
    pub spring_energy_change_j: f64,
    pub load_kinetic_energy_change_j: f64,
    pub lead_current_a: f64,
    pub lead_heat_j: f64,
    /// Positive torque acts on the load and negatively on the motor rotor.
    pub shaft_torque_nm: f64,
    pub shaft_heat_j: f64,
    /// Signed work leaving the rotor and entering the elastic/damped shaft.
    pub rotor_shaft_work_j: f64,
    /// Signed work leaving the shaft and entering the load inertia.
    pub load_shaft_work_j: f64,
    /// External rotor, stator and load torque work at ABSOLUTE shaft speeds.
    pub external_shaft_work_j: [f64; 3],
    /// Original work-conjugate battery and inverter bus efforts.
    pub discrete_bus_voltage_v: [f64; 2],
    /// Work-conjugate rotor, stator and load angular velocities [rad/s].
    pub discrete_shaft_speed_rad_s: [f64; 3],
}

impl ElasticShaftStep {
    /// Battery series/polarization, lead, eight inverter losses, shaft damping.
    pub fn heat_joules(&self) -> [f64; 12] {
        let mut heat = [0.0; 12];
        heat[..3].copy_from_slice(&[
            self.battery.series_heat_j,
            self.battery.polarization_heat_j,
            self.lead_heat_j,
        ]);
        heat[3..11].copy_from_slice(&self.inverter.heat_joules());
        heat[11] = self.shaft_heat_j;
        heat
    }

    pub fn shaft_balance_residual_j(&self) -> f64 {
        sum([
            self.spring_energy_change_j,
            self.shaft_heat_j,
            self.load_shaft_work_j,
            -self.rotor_shaft_work_j,
        ])
    }

    pub fn balance_residual_j(&self) -> f64 {
        sum(self.heat_joules().into_iter().chain([
            self.battery.storage_change_j(),
            self.inverter_storage_change_j,
            self.spring_energy_change_j,
            self.load_kinetic_energy_change_j,
            -self.external_shaft_work_j[0],
            -self.external_shaft_work_j[1],
            -self.external_shaft_work_j[2],
        ]))
    }
}

#[derive(Clone, Copy)]
struct ShaftStamp {
    drive: DriveStamp,
    resistance: [f64; 121],
    root: [f64; 121],
}

impl ShaftStamp {
    fn new(drive: DriveStamp, damping: f64) -> Result<Self, PhsError> {
        let mut resistance = [0.0; 121];
        let mut root = [0.0; 121];
        for i in 0..9 {
            for j in 0..9 {
                resistance[i * 11 + j] = drive.resistance[i * 9 + j];
                root[i * 11 + j] = drive.root[i * 9 + j];
            }
        }
        resistance[7 * 11 + 7] += damping;
        resistance[10 * 11 + 10] = damping;
        resistance[7 * 11 + 10] = -damping;
        resistance[10 * 11 + 7] = -damping;
        // One retained source square c*(omega_rotor-omega_load)^2. Original
        // battery, bridge, copper and rotor/stator bearing rows stay untouched.
        root[9 * 11 + 7] = det::sqrt(damping);
        root[9 * 11 + 10] = -root[9 * 11 + 7];
        if resistance.iter().chain(&root).any(|v| !v.is_finite()) {
            return Err(invalid(
                "elastic drive resistance or loss rows are unrepresentable",
            ));
        }
        crate::linear_exact::positive_heat::check_resistance_root(&resistance, &root, 11)?;
        Ok(Self {
            drive,
            resistance,
            root,
        })
    }
}

/// State order: the original battery drive's nine coordinates, shaft twist
/// relative to its stress-free phase [rad], then load angular momentum [N m s].
pub struct ElasticShaftDrive {
    drive: BatterySynchronousInverter,
    parameters: ElasticShaftParameters,
    system: PortHamiltonian,
    held: ShaftStamp,
}

impl BatterySynchronousInverter {
    /// Retain a finite load inertia and its elastic/viscous reaction in the
    /// original native solve. Zero stiffness/damping are admitted passive limits;
    /// infinite stiffness/inertia do not silently replace a missing constraint.
    pub fn with_elastic_shaft(
        mut self,
        parameters: ElasticShaftParameters,
    ) -> Result<ElasticShaftDrive, PhsError> {
        let p = parameters;
        if !p.load_inertia_kg_m2.is_finite()
            || p.load_inertia_kg_m2 <= 0.0
            || !p.load_inertia_kg_m2.recip().is_finite()
            || p.load_inertia_kg_m2.recip() <= 0.0
            || !p.stiffness_nm_per_rad.is_finite()
            || p.stiffness_nm_per_rad < 0.0
            || !p.damping_nm_s_per_rad.is_finite()
            || p.damping_nm_s_per_rad < 0.0
        {
            return Err(invalid(
                "elastic drive requires positive finite load inertia and nonnegative finite shaft coefficients",
            ));
        }
        let held = ShaftStamp::new(self.held, p.damping_nm_s_per_rad)?;
        let load = PortHamiltonian::new(
            2,
            1,
            vec![0.0, -1.0, 1.0, 0.0],
            vec![0.0; 4],
            vec![0.0, 1.0],
            Box::new(QuadraticStorage::new(
                vec![
                    p.stiffness_nm_per_rad,
                    0.0,
                    0.0,
                    p.load_inertia_kg_m2.recip(),
                ],
                2,
            )?),
        )?;
        let mut system = interconnect_mapped(
            self.system.share_storage(),
            load,
            &[],
            &[],
            &[],
            InterconnectionBudget {
                max_states: 11,
                max_total_ports: 5,
                max_dense_bytes: 1 << 18,
                max_coupling_multiply_terms: 256,
            },
        )?;
        // Native Dirac link: twist' gains omega_rotor and rotor momentum'
        // gains minus the spring effort. The load's J supplies -omega_load/+T.
        system.j[9 * 11 + 7] = 1.0;
        system.j[7 * 11 + 9] = -1.0;
        system.r.copy_from_slice(&held.resistance);
        system = system.with_resistance_witness(ResistanceWitness::try_from_source_rows(
            11,
            held.resistance.to_vec(),
            held.root.to_vec(),
        )?)?;
        Ok(ElasticShaftDrive {
            drive: self,
            parameters,
            system,
            held,
        })
    }
}

impl ElasticShaftDrive {
    /// Additive native ports: battery current, inverter current, then external
    /// rotor/stator/load torques. Prepared stepping holds both currents at zero.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn drive(&self) -> &BatterySynchronousInverter {
        &self.drive
    }
    pub fn parameters(&self) -> ElasticShaftParameters {
        self.parameters
    }

    /// Initialize every original drive coordinate, shaft strain and load speed
    /// independently. Preload and relative motion retain their original energy.
    pub fn state(
        &self,
        drive_state: &[f64; 9],
        twist_rad: f64,
        load_speed_rad_s: f64,
    ) -> Result<[f64; 11], PhsError> {
        let mut state = [0.0; 11];
        state[..9].copy_from_slice(drive_state);
        state[9] = twist_rad;
        state[10] = self.parameters.load_inertia_kg_m2 * load_speed_rad_s;
        if !load_speed_rad_s.is_finite() {
            return Err(invalid("elastic drive requires finite initial load speed"));
        }
        self.admit(&state)?;
        Ok(state)
    }

    fn admit(&self, state: &[f64; 11]) -> Result<(), PhsError> {
        self.drive.admit(&std::array::from_fn(|i| state[i]))?;
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(state).is_finite() {
            return Err(invalid(
                "elastic drive requires finite supported original storage",
            ));
        }
        Ok(())
    }

    fn apply_stamp(&mut self, stamp: &ShaftStamp) {
        self.drive.apply_stamp(&stamp.drive);
        for i in 0..9 {
            for j in 0..9 {
                self.system.j[i * 11 + j] = self.drive.system.j[i * 9 + j];
            }
        }
        self.system.r.copy_from_slice(&stamp.resistance);
        let witness = self
            .system
            .resistance_witness
            .as_mut()
            .expect("elastic drive retains the original source loss witness");
        witness.matrix.copy_from_slice(&stamp.resistance);
        witness.root.copy_from_slice(&stamp.root);
    }

    pub fn prepare_controlled(
        self,
        controls: BatteryInverterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElasticShaftDrive, PreparedStepError> {
        poll(&mut cancelled)?;
        if !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(
                invalid("elastic drive requires finite nonnegative work allowances").into(),
            );
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedElasticShaftDrive {
            model: self,
            workspace,
            controls,
        })
    }
}

/// Fixed eleven-state native scratch; successful warm modulation changes need
/// no allocation, new interconnection, lagged load torque or separate time step.
pub struct PreparedElasticShaftDrive {
    model: ElasticShaftDrive,
    workspace: StepWorkspace,
    controls: BatteryInverterControls,
}

fn quadratic_change(a: f64, b: f64, coefficient: f64) -> f64 {
    let (a, b) = (Dd::from_f64(a), Dd::from_f64(b));
    ((b - a) * (a * Dd::from_f64(0.5) + b * Dd::from_f64(0.5)) * Dd::from_f64(coefficient)).to_f64()
}

impl PreparedElasticShaftDrive {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn drive(&self) -> &BatterySynchronousInverter {
        self.model.drive()
    }
    pub fn parameters(&self) -> ElasticShaftParameters {
        self.model.parameters()
    }
    pub fn high_fraction(&self) -> [f64; 3] {
        self.model.held.drive.bridge.high_fraction
    }

    pub fn machine_operating_point(
        &self,
        state: &[f64; 11],
    ) -> Result<SynchronousOperatingPoint, PhsError> {
        self.model.admit(state)?;
        let inverter = &self.model.drive.inverter;
        storage::operating_point(
            inverter.machine_parameters,
            &inverter.machine_magnetic_law,
            &[state[4], state[5], state[6], state[7], state[8]],
        )
    }

    pub fn current_reference_for_torque_controlled(
        &self,
        requested_torque_nm: f64,
        d_current_a: f64,
        current_limit_a: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<reference::TorqueCurrentReference, PreparedStepError> {
        let inverter = &self.model.drive.inverter;
        reference::allocate(
            inverter.machine_parameters,
            &inverter.machine_magnetic_law,
            requested_torque_nm,
            d_current_a,
            current_limit_a,
            cancelled,
        )
    }

    /// Hold duty and external rotor/stator/load torques through the same native
    /// interval. Internal spring/damper torques use the accepted constitutive
    /// effort. Refusal or cancellation restores the whole held stamp and leaves
    /// caller state untouched; no state is projected at a duty or load change.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 11],
        high_fraction: [f64; 3],
        external_torque_nm: [f64; 3],
        duration_s: f64,
        next: &mut [f64; 11],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElasticShaftStep, PreparedStepError> {
        poll(&mut cancelled)?;
        self.model.admit(initial)?;
        let drive = &self.model.drive;
        let inverter = &drive.inverter;
        let bridge = BridgeStamp::new(
            inverter.parameters,
            inverter.machine_parameters,
            &inverter.base_resistance,
            inverter.scales,
            high_fraction,
        )?;
        let drive_stamp = DriveStamp::new(&drive.battery, bridge, drive.lead_resistance_ohm)?;
        let stamp = ShaftStamp::new(drive_stamp, self.model.parameters.damping_nm_s_per_rad)?;
        poll(&mut cancelled)?;
        self.model.apply_stamp(&stamp);
        let (mut candidate, mut output) = ([0.0; 11], [0.0; 5]);
        let result: Result<ElasticShaftStep, PreparedStepError> = (|| {
            let native = self.workspace.step_into_controlled(
                &self.model.system,
                initial,
                &[
                    0.0,
                    0.0,
                    external_torque_nm[0],
                    external_torque_nm[1],
                    external_torque_nm[2],
                ],
                duration_s,
                &mut candidate,
                &mut output,
                &mut cancelled,
            )?;
            self.model.admit(&candidate)?;
            let e = self
                .workspace
                .last_step_discrete_effort()
                .ok_or_else(|| invalid("elastic drive requires accepted constitutive efforts"))?;
            let drive = &self.model.drive;
            let p = self.model.parameters;
            let relative_speed = e[7] - e[10];
            let shaft_torque_nm = e[9] + p.damping_nm_s_per_rad * relative_speed;
            let lead_current_a = ((Dd::from_f64(e[2]) - Dd::from_f64(e[3]))
                / Dd::from_f64(drive.lead_resistance_ohm))
            .to_f64();
            let battery = drive.battery.interval(
                &[initial[0], initial[1], initial[2]],
                &[candidate[0], candidate[1], candidate[2]],
                -lead_current_a,
                duration_s,
            )?;
            let inverter = drive.inverter.interval(
                &stamp.drive.bridge,
                &[e[3], e[4], e[5], e[6], e[7], e[8]],
                &[
                    lead_current_a,
                    external_torque_nm[0] - shaft_torque_nm,
                    external_torque_nm[1],
                ],
                duration_s,
            )?;
            let inverter_storage_change_j = storage_energy_change(
                drive.inverter.system.storage(),
                &initial[3..9],
                &candidate[3..9],
                &mut [0.0; 6],
                &mut [0.0; 6],
            );
            let spring_energy_change_j =
                quadratic_change(initial[9], candidate[9], p.stiffness_nm_per_rad);
            let load_kinetic_energy_change_j =
                quadratic_change(initial[10], candidate[10], p.load_inertia_kg_m2.recip());
            let lead_heat_j =
                duration_s * lead_current_a * drive.lead_resistance_ohm * lead_current_a;
            let shaft_heat_j =
                duration_s * relative_speed * p.damping_nm_s_per_rad * relative_speed;
            let rotor_shaft_work_j = duration_s * shaft_torque_nm * e[7];
            let load_shaft_work_j = duration_s * shaft_torque_nm * e[10];
            let speeds = [e[7], e[8], e[10]];
            let external_shaft_work_j =
                std::array::from_fn(|i| duration_s * external_torque_nm[i] * speeds[i]);
            let step = ElasticShaftStep {
                native,
                battery,
                inverter,
                inverter_storage_change_j,
                spring_energy_change_j,
                load_kinetic_energy_change_j,
                lead_current_a,
                lead_heat_j,
                shaft_torque_nm,
                shaft_heat_j,
                rotor_shaft_work_j,
                load_shaft_work_j,
                external_shaft_work_j,
                discrete_bus_voltage_v: [e[2], e[3]],
                discrete_shaft_speed_rad_s: speeds,
            };
            let inverter_heat = sum(inverter.heat_joules());
            let heat = sum(step.heat_joules());
            let source_scale = sum([
                battery.chemical_energy_change_j.abs(),
                battery.polarization_energy_change_j.abs(),
                battery.bus_energy_change_j.abs(),
                battery.terminal_work_j.abs(),
                battery.heat_j(),
            ]);
            let inverter_scale = sum([
                inverter_storage_change_j.abs(),
                inverter_heat,
                inverter.source_work_j.abs(),
                inverter.shaft_work_j[0].abs(),
                inverter.shaft_work_j[1].abs(),
            ]);
            let lead_scale = sum([
                battery.terminal_work_j.abs(),
                inverter.source_work_j.abs(),
                lead_heat_j,
            ]);
            let shaft_scale = sum([
                spring_energy_change_j.abs(),
                shaft_heat_j,
                rotor_shaft_work_j.abs(),
                load_shaft_work_j.abs(),
            ]);
            let load_scale = sum([
                load_kinetic_energy_change_j.abs(),
                load_shaft_work_j.abs(),
                external_shaft_work_j[2].abs(),
            ]);
            let external_scale = sum(external_shaft_work_j.map(f64::abs));
            let change = sum([
                battery.storage_change_j(),
                inverter_storage_change_j,
                spring_energy_change_j,
                load_kinetic_energy_change_j,
            ]);
            let scale = source_scale
                .max(inverter_scale)
                .max(shaft_scale)
                .max(load_scale)
                .max(heat)
                .max(external_scale);
            let balances = [
                (battery.balance_residual_j(), source_scale),
                (
                    sum([
                        inverter_storage_change_j,
                        inverter_heat,
                        -inverter.source_work_j,
                        -inverter.shaft_work_j[0],
                        -inverter.shaft_work_j[1],
                    ]),
                    inverter_scale,
                ),
                (
                    sum([battery.terminal_work_j, inverter.source_work_j, lead_heat_j]),
                    lead_scale,
                ),
                (step.shaft_balance_residual_j(), shaft_scale),
                (
                    sum([
                        load_kinetic_energy_change_j,
                        -load_shaft_work_j,
                        -external_shaft_work_j[2],
                    ]),
                    load_scale,
                ),
                (
                    sum([
                        inverter.shaft_work_j[0],
                        rotor_shaft_work_j,
                        -external_shaft_work_j[0],
                    ]),
                    inverter_scale.max(shaft_scale),
                ),
                (sum([change, -native.delta_h]), scale),
                (sum([heat, -native.dissipated]), heat),
                (
                    sum([native.supplied, -sum(external_shaft_work_j)]),
                    external_scale,
                ),
                (step.balance_residual_j(), scale),
            ];
            if [
                inverter_storage_change_j,
                spring_energy_change_j,
                load_kinetic_energy_change_j,
                lead_current_a,
                shaft_torque_nm,
                rotor_shaft_work_j,
                load_shaft_work_j,
            ]
            .iter()
            .chain(&step.heat_joules())
            .chain(&external_shaft_work_j)
            .chain(&step.discrete_bus_voltage_v)
            .chain(&speeds)
            .any(|v| !v.is_finite())
                || balances.iter().any(|&(residual, scale)| {
                    let allowance = self.controls.absolute_energy_tolerance_j
                        + self.controls.relative_energy_tolerance * scale;
                    !scale.is_finite()
                        || !allowance.is_finite()
                        || !residual.is_finite()
                        || residual.abs() > allowance
                })
            {
                return Err(
                    invalid("battery, inverter, elastic shaft and load work do not close").into(),
                );
            }
            poll(&mut cancelled)?;
            Ok(step)
        })();
        match result {
            Ok(step) => {
                self.model.held = stamp;
                self.model.drive.held = stamp.drive;
                self.model.drive.inverter.held = stamp.drive.bridge;
                *next = candidate;
                Ok(step)
            }
            Err(error) => {
                let previous = self.model.held;
                self.model.apply_stamp(&previous);
                Err(error)
            }
        }
    }
}
