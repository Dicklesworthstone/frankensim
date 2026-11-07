//! Reversible nonlinear armature storage on the original reciprocal DC motor.
//!
//! The retained state is `[flux linkage, rotor momentum, stator momentum]`.
//! Current is the derivative of the source-bounded integral `H = integral I dλ`.
//! Saturation therefore changes current, differential inductance and magnetic
//! energy together. Neither `λ I / 2` nor a frozen secant L replaces that law.
//!
//! The existing scalar circuit storage owns exact finite magnetic work and its
//! analytic discrete-current derivative. The original prepared pHS solver uses
//! that work-conjugate current in BOTH winding voltage and reciprocal shaft
//! torque. Independent rotor/stator inertia and equal/opposite bearing reaction
//! remain present. No new nonlinear solver, current clipping or state projection
//! is introduced; native composition retains the same constitutive owner.
//!
//! Reported source work and positive heat use the native DISCRETE quadrature.
//! In particular R*I_bar²*dt is not an exact continuous RMS-current integral.
//! Nonlinear current I_bar is the magnetic-energy secant over flux change, not
//! generally the arithmetic mean of endpoint currents or a time average.
//! Time refinement and source interpolation accuracy remain separate duties.
//! Frozen k and reversible scalar I(λ) add no iron hysteresis, magnetic angle
//! dependence, commutation, brush loss, demagnetization or thermal parameter law.

use super::{DcMotor, DcMotorParameters, invalid, poll};
use crate::circuit::{
    CircuitError, capacitance::storage_with_inductor_curve, inductance::InductorFluxCurve,
};
use crate::{
    PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage, SolveBudget,
    StepWorkspace,
};
use fs_math::dd::Dd;

/// Frozen mechanical, reciprocal-coupling and loss data. The separate original
/// flux/current curve is the complete winding storage; there is no placeholder L.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearDcMotorParameters {
    pub winding_resistance_ohm: f64,
    pub rotor_inertia_kg_m2: f64,
    pub stator_inertia_kg_m2: f64,
    /// Signed torque/current and back-EMF/relative-speed coefficient in SI units.
    pub coupling_nm_per_a: f64,
    /// Internal rotor-to-stator viscous drag, not a torque to ground.
    pub bearing_drag_nm_s_per_rad: f64,
}

/// Explicit native Newton bounds and independent physical-balance allowances.
/// The fixed three-state/three-port workspace and bounded source knot count
/// bound storage and each callback. No tolerance is raised to accept a trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearMotorControls {
    pub solve: SolveBudget,
    /// Existing native globalization, 0..=16. Out-of-source probes still refuse.
    pub backtracking_halvings: u8,
    pub absolute_energy_tolerance_j: f64,
    pub absolute_momentum_tolerance_nms: f64,
    /// Relative to each boundary's actual transfers, excluding unchanged energy.
    pub relative_balance_tolerance: f64,
}
impl Default for NonlinearMotorControls {
    fn default() -> Self {
        Self {
            solve: SolveBudget::default(),
            backtracking_halvings: 0,
            absolute_energy_tolerance_j: 0.0,
            absolute_momentum_tolerance_nms: 0.0,
            relative_balance_tolerance: 1e-12,
        }
    }
}

/// Three original flux/momentum states and one exact bounded winding source.
pub struct NonlinearDcMotor {
    parameters: NonlinearDcMotorParameters,
    curve: InductorFluxCurve,
    system: PortHamiltonian,
}

/// Native discrete work, storage and both positive heat owners. None of these
/// accounts is inferred by subtracting total heat from a supplied-work deficit.
#[derive(Clone, Debug)]
pub struct NonlinearMotorStep {
    pub native: PreparedStepRecord,
    pub magnetic_energy_change_j: f64,
    /// Independent absolute rotor/stator kinetic-energy changes, in that order.
    pub kinetic_energy_change_j: [f64; 2],
    /// Signed electrical work into the armature; negative denotes generation.
    pub electrical_work_j: f64,
    /// Signed work of the externally supplied rotor/stator torques.
    pub shaft_work_j: [f64; 2],
    /// Reciprocal internal electromagnetic work on rotor/stator respectively.
    pub electromagnetic_work_j: [f64; 2],
    /// Internal bearing work on rotor/stator respectively; their sum is -heat.
    pub bearing_work_j: [f64; 2],
    pub copper_heat_j: f64,
    pub bearing_heat_j: f64,
    pub discrete_current_a: f64,
    pub discrete_shaft_speed_rad_s: [f64; 2],
    /// Torque on rotor; the stator receives its exact opposite.
    pub discrete_electromagnetic_torque_nm: f64,
    /// Bearing torque on rotor; the stator receives its exact opposite.
    pub discrete_bearing_torque_nm: f64,
    pub discrete_back_emf_v: f64,
    /// Discrete angular increments paired with the held external torques.
    pub shaft_angle_increment_rad: [f64; 2],
    pub angular_momentum_residual_nms: f64,
}

impl NonlinearDcMotor {
    /// Admit the SAME motor J/R/G and source-square resistance witness as the
    /// linear owner, then retain the complete original nonlinear winding law.
    /// Its zero-flux differential L supplies only that existing constructor's
    /// structural admission; it is never used for nonlinear current or energy.
    pub fn new(
        parameters: NonlinearDcMotorParameters,
        curve: InductorFluxCurve,
    ) -> Result<Self, CircuitError> {
        let p = parameters;
        let mut system = DcMotor::new(DcMotorParameters {
            inductance_h: curve.differential_inductance_h(0.0)?,
            winding_resistance_ohm: p.winding_resistance_ohm,
            rotor_inertia_kg_m2: p.rotor_inertia_kg_m2,
            stator_inertia_kg_m2: p.stator_inertia_kg_m2,
            coupling_nm_per_a: p.coupling_nm_per_a,
            bearing_drag_nm_s_per_rad: p.bearing_drag_nm_s_per_rad,
        })?
        .into_system();
        system.storage = storage_with_inductor_curve(
            QuadraticStorage::new(
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0 / p.rotor_inertia_kg_m2,
                    0.0,
                    0.0,
                    0.0,
                    1.0 / p.stator_inertia_kg_m2,
                ],
                3,
            )?,
            0,
            &curve,
        )?;
        Ok(Self {
            parameters,
            curve,
            system,
        })
    }

    pub const fn parameters(&self) -> NonlinearDcMotorParameters {
        self.parameters
    }
    pub const fn flux_curve(&self) -> &InductorFluxCurve {
        &self.curve
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Move the complete original bounded storage into native interconnection.
    /// Bus-voltage/current coupling, reaction torque and constitutive analytic
    /// actions survive composition; nonlinear sources require a nonlinear step.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Initialize original flux from source-supported current and both absolute
    /// shaft momenta. No constitutive extrapolation or fixed-stator substitution.
    pub fn state(
        &self,
        current_a: f64,
        shaft_speed_rad_s: [f64; 2],
    ) -> Result<[f64; 3], CircuitError> {
        let flux = self.curve.flux_wb(current_a)?;
        let state = [
            flux,
            self.parameters.rotor_inertia_kg_m2 * shaft_speed_rad_s[0],
            self.parameters.stator_inertia_kg_m2 * shaft_speed_rad_s[1],
        ];
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(
                invalid("nonlinear motor state or stored energy is unrepresentable").into(),
            );
        }
        Ok(state)
    }

    /// Reuse the native prepared discrete-gradient solver and its exact source
    /// derivative. Cold preparation allocates fixed three-coordinate scratch;
    /// successful warm calls add no allocation in this owner or its callbacks.
    pub fn prepare_controlled(
        self,
        duration_s: f64,
        controls: NonlinearMotorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedNonlinearDcMotor, PreparedStepError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || [
                controls.absolute_energy_tolerance_j,
                controls.absolute_momentum_tolerance_nms,
                controls.relative_balance_tolerance,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "nonlinear motor needs positive duration and finite nonnegative allowances",
            )
            .into());
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedNonlinearDcMotor {
            motor: self,
            workspace,
            duration_s,
            controls,
        })
    }
}

pub struct PreparedNonlinearDcMotor {
    motor: NonlinearDcMotor,
    workspace: StepWorkspace,
    duration_s: f64,
    controls: NonlinearMotorControls,
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}
fn product(a: f64, b: f64, c: f64) -> f64 {
    (Dd::from_f64(a) * Dd::from_f64(b) * Dd::from_f64(c)).to_f64()
}

impl PreparedNonlinearDcMotor {
    pub fn system(&self) -> &PortHamiltonian {
        self.motor.system()
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn parameters(&self) -> NonlinearDcMotorParameters {
        self.motor.parameters
    }
    pub const fn flux_curve(&self) -> &InductorFluxCurve {
        &self.motor.curve
    }

    /// Hold voltage and external rotor/stator torque over the interval.
    /// `discrete_outputs` is [work-conjugate current, midpoint rotor speed,
    /// midpoint stator speed]. In general it is NOT the time-average output.
    /// Curve support applies to every original native probe. Any source,
    /// iteration, physical balance or cancellation refusal preserves BOTH caller
    /// output buffers and leaves this cache safe for deterministic retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        inputs: &[f64; 3],
        next: &mut [f64; 3],
        discrete_outputs: &mut [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearMotorStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let (mut candidate, mut effort) = ([0.0; 3], [0.0; 3]);
        let native = self.workspace.step_into_controlled(
            &self.motor.system,
            initial,
            inputs,
            self.duration_s,
            &mut candidate,
            &mut effort,
            &mut cancelled,
        )?;
        let p = self.motor.parameters;
        let h = self.duration_s;
        // This is the original source integral, independent of supplied work
        // and the native whole-system stored-energy diagnostic.
        let magnetic_energy_change_j =
            self.motor
                .curve
                .energy_change_j(initial[0], candidate[0])
                .map_err(|_| invalid("nonlinear motor magnetic work is outside source support"))?;
        let kinetic_energy_change_j: [f64; 2] = core::array::from_fn(|i| {
            let before = Dd::from_f64(initial[i + 1]);
            let after = Dd::from_f64(candidate[i + 1]);
            let inertia = [p.rotor_inertia_kg_m2, p.stator_inertia_kg_m2][i];
            ((after - before) * (Dd::from_f64(0.5) * before + Dd::from_f64(0.5) * after)
                / Dd::from_f64(inertia))
            .to_f64()
        });
        let relative_speed = sum([effort[1], -effort[2]]);
        let discrete_electromagnetic_torque_nm = p.coupling_nm_per_a * effort[0];
        let discrete_bearing_torque_nm = -p.bearing_drag_nm_s_per_rad * relative_speed;
        let discrete_back_emf_v = p.coupling_nm_per_a * relative_speed;
        let copper_heat_j = product(h * p.winding_resistance_ohm, effort[0], effort[0]);
        let bearing_heat_j = product(
            h * p.bearing_drag_nm_s_per_rad,
            relative_speed,
            relative_speed,
        );
        let electrical_work_j = product(h, inputs[0], effort[0]);
        let shaft_angle_increment_rad: [f64; 2] = [h * effort[1], h * effort[2]];
        let shaft_work_j: [f64; 2] =
            core::array::from_fn(|i| product(h, inputs[i + 1], effort[i + 1]));
        let electromagnetic_work_j: [f64; 2] = core::array::from_fn(|i| {
            product(
                if i == 0 { h } else { -h },
                discrete_electromagnetic_torque_nm,
                effort[i + 1],
            )
        });
        let bearing_work_j: [f64; 2] = core::array::from_fn(|i| {
            product(
                if i == 0 { h } else { -h },
                discrete_bearing_torque_nm,
                effort[i + 1],
            )
        });
        let electromagnetic_transfer = product(h, discrete_back_emf_v, effort[0]);
        let angular_momentum_residual_nms = (Dd::from_f64(candidate[1]) - Dd::from_f64(initial[1])
            + Dd::from_f64(candidate[2])
            - Dd::from_f64(initial[2])
            - Dd::from_f64(h) * (Dd::from_f64(inputs[1]) + Dd::from_f64(inputs[2])))
        .to_f64();
        let heat = sum([copper_heat_j, bearing_heat_j]);
        let supplied = sum([electrical_work_j, shaft_work_j[0], shaft_work_j[1]]);
        let stored = sum([
            magnetic_energy_change_j,
            kinetic_energy_change_j[0],
            kinetic_energy_change_j[1],
        ]);
        let accounts = [
            (
                sum([
                    magnetic_energy_change_j,
                    copper_heat_j,
                    electromagnetic_transfer,
                    -electrical_work_j,
                ]),
                sum([
                    magnetic_energy_change_j.abs(),
                    copper_heat_j,
                    electromagnetic_transfer.abs(),
                    electrical_work_j.abs(),
                ]),
            ),
            (
                sum([
                    kinetic_energy_change_j[0],
                    -shaft_work_j[0],
                    -electromagnetic_work_j[0],
                    -bearing_work_j[0],
                ]),
                sum([
                    kinetic_energy_change_j[0].abs(),
                    shaft_work_j[0].abs(),
                    electromagnetic_work_j[0].abs(),
                    bearing_work_j[0].abs(),
                ]),
            ),
            (
                sum([
                    kinetic_energy_change_j[1],
                    -shaft_work_j[1],
                    -electromagnetic_work_j[1],
                    -bearing_work_j[1],
                ]),
                sum([
                    kinetic_energy_change_j[1].abs(),
                    shaft_work_j[1].abs(),
                    electromagnetic_work_j[1].abs(),
                    bearing_work_j[1].abs(),
                ]),
            ),
            (
                sum([
                    electromagnetic_work_j[0],
                    electromagnetic_work_j[1],
                    -electromagnetic_transfer,
                ]),
                sum(electromagnetic_work_j.map(f64::abs)).max(electromagnetic_transfer.abs()),
            ),
            (
                sum([bearing_work_j[0], bearing_work_j[1], bearing_heat_j]),
                sum(bearing_work_j.map(f64::abs)).max(bearing_heat_j),
            ),
            (
                sum([stored, -native.delta_h]),
                sum([
                    magnetic_energy_change_j.abs(),
                    kinetic_energy_change_j[0].abs(),
                    kinetic_energy_change_j[1].abs(),
                ]),
            ),
            (
                sum([supplied, -native.supplied]),
                sum([
                    electrical_work_j.abs(),
                    shaft_work_j[0].abs(),
                    shaft_work_j[1].abs(),
                ]),
            ),
            (sum([heat, -native.dissipated]), heat),
            (
                sum([stored, heat, -supplied]),
                sum([
                    magnetic_energy_change_j.abs(),
                    kinetic_energy_change_j[0].abs(),
                    kinetic_energy_change_j[1].abs(),
                    heat,
                    electrical_work_j.abs(),
                    shaft_work_j[0].abs(),
                    shaft_work_j[1].abs(),
                ]),
            ),
        ];
        let momentum_scale = sum([
            (candidate[1] - initial[1]).abs(),
            (candidate[2] - initial[2]).abs(),
            (h * inputs[1]).abs(),
            (h * inputs[2]).abs(),
        ]);
        let momentum_tolerance = self.controls.absolute_momentum_tolerance_nms
            + self.controls.relative_balance_tolerance * momentum_scale;
        if [
            magnetic_energy_change_j,
            relative_speed,
            discrete_electromagnetic_torque_nm,
            discrete_bearing_torque_nm,
            discrete_back_emf_v,
            copper_heat_j,
            bearing_heat_j,
            electrical_work_j,
            angular_momentum_residual_nms,
            momentum_scale,
            momentum_tolerance,
        ]
        .iter()
        .chain(&kinetic_energy_change_j)
        .chain(&shaft_work_j)
        .chain(&electromagnetic_work_j)
        .chain(&bearing_work_j)
        .chain(&shaft_angle_increment_rad)
        .any(|v| !v.is_finite())
            || copper_heat_j < 0.0
            || bearing_heat_j < 0.0
            || angular_momentum_residual_nms.abs() > momentum_tolerance
            || accounts.iter().any(|&(residual, scale)| {
                let tolerance = self.controls.absolute_energy_tolerance_j
                    + self.controls.relative_balance_tolerance * scale;
                !residual.is_finite()
                    || !scale.is_finite()
                    || !tolerance.is_finite()
                    || residual.abs() > tolerance
            })
        {
            return Err(invalid(
                "nonlinear motor magnetic, shaft, heat and momentum accounts do not close",
            )
            .into());
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        discrete_outputs.copy_from_slice(&effort);
        Ok(NonlinearMotorStep {
            native,
            magnetic_energy_change_j,
            kinetic_energy_change_j,
            electrical_work_j,
            shaft_work_j,
            electromagnetic_work_j,
            bearing_work_j,
            copper_heat_j,
            bearing_heat_j,
            discrete_current_a: effort[0],
            discrete_shaft_speed_rad_s: [effort[1], effort[2]],
            discrete_electromagnetic_torque_nm,
            discrete_bearing_torque_nm,
            discrete_back_emf_v,
            shaft_angle_increment_rad,
            angular_momentum_residual_nms,
        })
    }
}
