//! G1: the complete finite-source drive against independently integrated SI ODEs.
//! The reference writes KCL, winding KVL and both shaft momentum equations
//! directly; it uses neither a port matrix nor a production work/heat helper.
use fs_phs::SolveBudget;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::power_load::PowerLoadControls;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};

const CAPACITY: f64 = 80.0;
const LEAD: f64 = 0.08;
const TORQUES: [f64; 2] = [-0.2, 0.03];
const END: f64 = 0.12;
const BATTERY: BatteryParameters = BatteryParameters {
    series_resistance_ohm: 0.1,
    polarization_resistance_ohm: 0.2,
    polarization_capacitance_f: 0.5,
    bus_capacitance_f: 0.03,
};
const CONVERTER: BuckBoostParameters = BuckBoostParameters {
    inductance_h: 0.04,
    bus_capacitance_f: [0.04, 0.05],
    winding_resistance_ohm: 0.07,
    input_on_resistance_ohm: [0.02, 0.03],
    output_on_resistance_ohm: [0.02, 0.03],
};
const MODULATION: BridgeModulation = BridgeModulation {
    input_high_fraction: 0.55,
    output_high_fraction: 0.85,
};
const MOTOR: DcMotorParameters = DcMotorParameters {
    inductance_h: 0.02,
    winding_resistance_ohm: 0.2,
    rotor_inertia_kg_m2: 0.08,
    stator_inertia_kg_m2: 0.4,
    coupling_nm_per_a: 0.4,
    bearing_drag_nm_s_per_rad: 0.01,
};

fn initial() -> [f64; 9] {
    [
        CAPACITY * 0.65,
        0.5 * 0.1,
        0.03 * 12.1,
        0.04 * 0.5,
        0.04 * 11.8,
        0.05 * 6.4,
        0.02 * 0.1,
        0.08 * 4.0,
        0.4 * -0.4,
    ]
}
fn model() -> BatteryMotorDrive {
    let battery = BatterySource::new(
        BatteryOcvCurve::new(CAPACITY, &[(0.0, 10.0), (1.0, 14.0)]).unwrap(),
        BATTERY,
    )
    .unwrap();
    AveragedBuckBoost::new(CONVERTER, MODULATION)
        .unwrap()
        .with_battery_source(battery, LEAD)
        .unwrap()
        .with_dc_motor(DcMotor::new(MOTOR).unwrap())
        .unwrap()
}

fn derivative(x: &[f64; 9]) -> [f64; 9] {
    let ocv = 10.0 + 4.0 * x[0] / CAPACITY;
    let polarization = x[1] / BATTERY.polarization_capacitance_f;
    let battery_v = x[2] / BATTERY.bus_capacitance_f;
    let converter_i = x[3] / CONVERTER.inductance_h;
    let input_v = x[4] / CONVERTER.bus_capacitance_f[0];
    let output_v = x[5] / CONVERTER.bus_capacitance_f[1];
    let motor_i = x[6] / MOTOR.inductance_h;
    let relative_speed = x[7] / MOTOR.rotor_inertia_kg_m2 - x[8] / MOTOR.stator_inertia_kg_m2;
    let cell_i = (ocv - polarization - battery_v) / BATTERY.series_resistance_ohm;
    let lead_i = (battery_v - input_v) / LEAD;
    let a = MODULATION.input_high_fraction;
    let b = MODULATION.output_high_fraction;
    let converter_r = CONVERTER.winding_resistance_ohm
        + a * CONVERTER.input_on_resistance_ohm[0]
        + (1.0 - a) * CONVERTER.input_on_resistance_ohm[1]
        + b * CONVERTER.output_on_resistance_ohm[0]
        + (1.0 - b) * CONVERTER.output_on_resistance_ohm[1];
    let internal_torque =
        MOTOR.coupling_nm_per_a * motor_i - MOTOR.bearing_drag_nm_s_per_rad * relative_speed;
    [
        -cell_i,
        cell_i - polarization / BATTERY.polarization_resistance_ohm,
        cell_i - lead_i,
        a * input_v - b * output_v - converter_r * converter_i,
        lead_i - a * converter_i,
        b * converter_i - motor_i,
        output_v
            - MOTOR.winding_resistance_ohm * motor_i
            - MOTOR.coupling_nm_per_a * relative_speed,
        internal_torque + TORQUES[0],
        -internal_torque + TORQUES[1],
    ]
}

fn reference(steps: usize) -> [f64; 9] {
    let dt = END / steps as f64;
    let mut x = initial();
    for _ in 0..steps {
        let a = derivative(&x);
        let b = derivative(&core::array::from_fn(|i| x[i] + dt * 0.5 * a[i]));
        let c = derivative(&core::array::from_fn(|i| x[i] + dt * 0.5 * b[i]));
        let d = derivative(&core::array::from_fn(|i| x[i] + dt * c[i]));
        x = core::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn g1_reciprocal_finite_drive_matches_independent_dynamics_with_second_order_refinement() {
    let fine_reference = reference(32768);
    let coarser_reference = reference(16384);
    for i in 0..9 {
        assert!((fine_reference[i] - coarser_reference[i]).abs() < 1e-10);
    }
    let controls = PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-11,
    };
    let mut errors = Vec::new();
    for steps in [32, 64, 128] {
        let mut solver = model()
            .prepare_controlled(END / steps as f64, controls, || false)
            .unwrap();
        let mut x = initial();
        for _ in 0..steps {
            let mut next = [0.0; 9];
            solver
                .step_controlled(&x, &TORQUES, None, &mut next, || false)
                .unwrap();
            x = next;
        }
        // Each component uses one fixed SI unit of its coordinate as its error
        // scale (coulombs, webers or angular momentum), shared by all refinements.
        let error = x
            .iter()
            .zip(fine_reference)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        let momentum_change = (x[7] - initial()[7]) + (x[8] - initial()[8]);
        assert!((momentum_change - END * (TORQUES[0] + TORQUES[1])).abs() < 2e-11);
    }
    assert!(errors[2] < 2e-4, "full physical state error: {errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "second-order refinement: {errors:?}"
        );
    }
}
