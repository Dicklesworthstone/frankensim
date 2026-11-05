//! Independent reciprocal rod, circuit and four-cell Kelvin heat equations.
//! Coefficients are authored numerical fixtures, not measured active materials.
use super::*;
use fs_couple::circuit_enthalpy::piezoelectric::{self as heated, PiezoelectricEnthalpyInput};
use fs_couple::circuit_enthalpy::thermoelectric::piezoelectric as joined;
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_matdb::{ClaimSet, MaterialStateId};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::{
    COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, LongitudinalPiezoelectric,
    PERMITTIVITY_33_T_PROPERTY, STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS,
};

fn rod(d: f64) -> PiezoelectricRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, d),
        (PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, 2.0),
        ("density", Density::DIMS, 3.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored reciprocal thermocouple actuator".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic ceramic".into(),
            phase: "solid".into(),
            process: "uniform longitudinal fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = LongitudinalPiezoelectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PiezoelectricRod::new(
        material,
        2.0,
        0.5,
        SmallSignalLimits {
            max_abs_strain: 0.05,
            max_abs_electric_field_v_m: 10.0,
        },
    )
    .unwrap()
}

fn actuator(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    d: f64,
    state: [f64; 3],
    current: f64,
) -> PiezoelectricCircuit {
    PiezoelectricCircuit::new(
        CircuitGraph {
            node_count: 5,
            branches: vec![
                (1, 2, Branch::VoltageSource { port: 0 }),
                (2, 3, Branch::Resistor { ohms: f64::NAN }),
                (3, 4, Branch::VoltageSource { port: 1 }),
                (4, 0, Branch::Resistor { ohms: f64::NAN }),
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::CurrentSource { port: 2 }),
            ],
            transformers: vec![],
        },
        material.resolve(thermal, || false).unwrap().resistors,
        vec![],
        vec![state[2]],
        4,
        rod(d),
        state[0],
        state[1],
        &[0.0, 0.0, current],
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9, "{a:.16e} != {b:.16e}");
}
fn state(model: &PiezoelectricCircuit) -> [f64; 3] {
    [model.extension_m(), model.velocity_m_s(), model.charge_c()]
}
fn energy(state: [f64; 3], d: f64) -> f64 {
    let [x, v, q] = state;
    let theta = 2.5 * d;
    let capacitance = 0.5 - 2.5 * d * d;
    0.5 * v * v + 1.25 * x * x + (q - theta * x).powi(2) / (2.0 * capacitance)
}

// Eliminate the native midpoint equations independently for this fixture:
// xdot=v, vdot=F-2.5*x+theta*V, q=C*V+theta*x,
// qdot=I_external-(V-E_thermocouple)/5. No alternate time integrator is used.
fn midpoint(
    old: [f64; 3],
    d: f64,
    dt: f64,
    voltage: f64,
    current: f64,
    force: f64,
) -> ([f64; 3], f64, f64) {
    let [x, v, q] = old;
    let theta = 2.5 * d;
    let capacitance = 0.5 - 2.5 * d * d;
    let stiffness = 2.5 + 4.0 / (dt * dt);
    let a = (force + 4.0 * x / (dt * dt) + 2.0 * v / dt) / stiffness;
    let b = theta / stiffness;
    let middle_voltage = (q - theta * a + 0.5 * dt * (current + voltage / 5.0))
        / (capacitance + theta * b + dt / 10.0);
    let middle_x = a + b * middle_voltage;
    let leg_current = (middle_voltage - voltage) / 5.0;
    (
        [
            2.0 * middle_x - x,
            4.0 * (middle_x - x) / dt - v,
            q + dt * (current - leg_current),
        ],
        middle_voltage,
        leg_current,
    )
}

#[test]
fn g1_finite_gradient_drives_piezo_motion_with_complete_kelvin_heat_and_polarity() {
    let run = |sign: f64, d: f64| {
        let (s1, s2, slope) = (sign * 0.02, -sign * 0.01, sign * 1e-4);
        let (mut cards, mut curves) = fixture(s1, s2, 3.0, 1.0);
        (cards[2], curves[2]) = source_value(
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, s1 - 150.0 * slope), (2000.0, s1 + 1650.0 * slope)],
                dims: SEEBECK_DIMS,
            },
            InterpolationPolicy::LinearInside,
            2.0,
            1.0,
        );
        let context = point(350.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = resistors();
        let links = links();
        let material = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let mut thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let mut model = actuator(&material, &thermal, d, [0.0; 3], 0.0);
        let initial_energy = thermal.energy_j();
        let mut last_current = 0.0;
        for index in 1..=6 {
            let target = f64::from(index) * 0.1;
            let dt = target - thermal.time_s();
            let before = state(&model);
            let old_t = [0, 1, 2, 3].map(|i| thermal.temperature_k(CellId(i)).unwrap());
            let next = joined::step_candidate(
                &material,
                &thermal,
                target,
                &model,
                &[0.0; 3],
                0.0,
                &legs(),
                &BASES,
                config(),
                &[0.0; 4],
                || false,
            )
            .unwrap();
            assert_eq!(next.thermal().time_s(), target);
            let t = [0, 1, 2, 3].map(|i| next.thermal().temperature_k(CellId(i)).unwrap());
            let [a, b, c, e] = t;
            let left = s1 + slope * (0.5 * (a + c) - 350.0);
            let right = s1 + slope * (0.5 * (c + b) - 350.0);
            let first_voltage = (s1 + slope * (0.5 * (a + b) - 350.0)) * (b - a);
            let second_voltage = s2 * (a - b);
            let (expected, middle_voltage, current) =
                midpoint(before, d, dt, first_voltage + second_voltage, 0.0, 0.0);
            for (actual, expected) in state(next.model()).into_iter().zip(expected) {
                close(actual, expected);
            }
            close(
                0.5 * (model.voltage_v().unwrap() + next.model().voltage_v().unwrap()),
                middle_voltage,
            );
            let report = next.report();
            let samples = &report.coupled.source.legs;
            close(samples[0].transfer.current_a, current);
            close(samples[1].transfer.current_a, current);
            close(samples[0].transfer.voltage_drop_v, first_voltage);
            close(samples[1].transfer.voltage_drop_v, second_voltage);
            close(samples[0].transfer.terminal_heat_w[0], -left * a * current);
            close(samples[0].transfer.terminal_heat_w[1], right * b * current);
            close(
                samples[0].transfer.bulk_heat_w,
                c * (left - right) * current,
            );
            // Four independent finite-cell BE balances, including both legs'
            // junction contributions, original Fourier links and bulk Joule heat.
            let powers = [
                c + e - 2.0 * a + (s2 - left) * a * current,
                c + e - 2.0 * b + (right - s2) * b * current,
                a + b - 2.0 * c + 2.0 * current * current + c * (left - right) * current,
                a + b - 2.0 * e + 3.0 * current * current,
            ];
            for ((next_t, old_t), power) in t.into_iter().zip(old_t).zip(powers) {
                near(10.0 * (next_t - old_t), dt * power);
            }
            close(
                report.coupled.deposited_joule_j,
                5.0 * dt * current * current,
            );
            close(
                report.coupled.source.reversible_electrical_work_j,
                -dt * current * (first_voltage + second_voltage),
            );
            close(
                report.coupled.source.reversible_heat_j,
                -report.coupled.source.reversible_electrical_work_j,
            );
            close(report.external_electrical_work_j, 0.0);
            close(report.external_mechanical_work_j, 0.0);
            near(report.external_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            close(report.coupled.electrical.parameter_work_j, 0.0);
            close(
                report.coupled.electrical.stored_after_j,
                energy(expected, d),
            );
            assert!(
                (next.thermal().energy_j() + energy(expected, d) - initial_energy).abs() < 2e-5
            );
            assert_eq!(
                samples[0].material.material()[0].card_identity(),
                cards[2].content_hash()
            );
            assert_eq!(model.rod(), next.model().rod());
            assert_eq!(state(&model), before);
            assert_eq!(
                [0, 1, 2, 3].map(|i| thermal.temperature_k(CellId(i)).unwrap()),
                old_t
            );
            last_current = current;
            (thermal, model, _) = next.into_parts();
        }
        assert!(sign * model.charge_c() < 0.0 && sign * d * model.extension_m() < 0.0);
        assert!(model.velocity_m_s().abs() > 1e-4 && thermal.energy_j() < initial_energy);
        [
            last_current,
            model.charge_c(),
            model.extension_m(),
            model.velocity_m_s(),
            thermal.temperature_k(CellId(0)).unwrap(),
            thermal.temperature_k(CellId(1)).unwrap(),
            thermal.temperature_k(CellId(2)).unwrap(),
            thermal.temperature_k(CellId(3)).unwrap(),
        ]
    };
    let original = run(1.0, 0.2);
    let reverse_thermal = run(-1.0, 0.2);
    let reverse_poling = run(1.0, -0.2);
    for i in 0..4 {
        close(original[i], -reverse_thermal[i]);
    }
    for i in 4..8 {
        close(original[i], reverse_thermal[i]);
        close(original[i], reverse_poling[i]);
    }
    for i in 0..2 {
        close(original[i], reverse_poling[i]);
    }
    for i in 2..4 {
        close(original[i], -reverse_poling[i]);
    }
}

#[test]
fn g1_mechanical_work_and_authored_current_remain_external_while_kelvin_heat_is_internal() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(350.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistors();
    let links = links();
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [350.0; 4]);
    for drive in [0.0, 0.08] {
        let model = actuator(&material, &thermal, 0.2, [0.0; 3], drive);
        let next = joined::step_candidate(
            &material,
            &thermal,
            0.2,
            &model,
            &[0.0, 0.0, drive],
            0.3,
            &legs(),
            &BASES,
            config(),
            &[0.0; 4],
            || false,
        )
        .unwrap();
        let [a, b] = [0, 1].map(|i| next.thermal().temperature_k(CellId(i)).unwrap());
        let (expected, voltage, current) = midpoint([0.0; 3], 0.2, 0.2, 0.03 * (b - a), drive, 0.3);
        for (actual, expected) in state(next.model()).into_iter().zip(expected) {
            close(actual, expected);
        }
        let report = next.report();
        close(report.external_electrical_work_j, 0.2 * drive * voltage);
        close(report.external_mechanical_work_j, 0.3 * expected[0]);
        near(report.external_heat_j, 0.0);
        near(
            next.thermal().energy_j() - thermal.energy_j() + energy(expected, 0.2),
            report.external_electrical_work_j + report.external_mechanical_work_j,
        );
        assert!(report.external_mechanical_work_j > 0.0 && (a - b) * current < 0.0);
        if drive == 0.0 {
            assert!(current < 0.0 && a > b);
        }
        close(
            report.coupled.source.reversible_heat_j,
            -report.coupled.source.reversible_electrical_work_j,
        );
        close(report.coupled.deposited_joule_j, current * current);
    }
}

#[test]
fn g3_zero_seebeck_matches_existing_piezoelectric_enthalpy_and_late_refusal_is_retryable() {
    let (cards, curves) = fixture(0.0, 0.0, 3.0, 1.0);
    let context = point(350.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistors();
    let links = links();
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let model = actuator(&material, &thermal, 0.2, [0.01, 0.02, 0.05], 0.04);
    let prepared = material
        .prepare_with_bases(&thermal, &BASES, || false)
        .unwrap();
    let destinations = prepared.resistor_destinations();
    let original = heated::step_candidate(
        &thermal,
        0.1,
        PiezoelectricEnthalpyInput {
            model: &model,
            electrical_ports: &[0.0, 0.0, 0.04],
            force_n: 0.3,
            destinations: &destinations,
        },
        config(),
        &[0.0; 4],
        |state| prepared.resolve(state, || false),
        || false,
    )
    .unwrap();
    let mut polls = 0;
    let full = joined::step_candidate(
        &material,
        &thermal,
        0.1,
        &model,
        &[0.0, 0.0, 0.04],
        0.3,
        &legs(),
        &BASES,
        config(),
        &[0.0; 4],
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(full.thermal().cells(), original.thermal().cells());
    assert_eq!(full.model(), original.model());
    close(full.report().coupled.source.reversible_heat_j, 0.0);
    close(
        full.report().coupled.source.reversible_electrical_work_j,
        0.0,
    );
    let before = model.clone();
    let before_cells = thermal.cells().clone();
    let mut cancelled = 0;
    assert!(matches!(
        joined::step_candidate(
            &material,
            &thermal,
            0.1,
            &model,
            &[0.0, 0.0, 0.04],
            0.3,
            &legs(),
            &BASES,
            config(),
            &[0.0; 4],
            || {
                cancelled += 1;
                cancelled == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut bounded = config();
    bounded.maximum_circuit_trials = 1;
    assert!(matches!(
        joined::step_candidate(
            &material,
            &thermal,
            0.1,
            &model,
            &[0.0, 0.0, 0.04],
            0.3,
            &legs(),
            &BASES,
            bounded,
            &[0.0; 4],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(
        joined::step_candidate(
            &material,
            &thermal,
            0.1,
            &model,
            &[0.0, 0.0, 0.04],
            1e6,
            &legs(),
            &BASES,
            config(),
            &[0.0; 4],
            || false
        )
        .is_err()
    );
    assert_eq!(model, before);
    assert_eq!(thermal.cells(), &before_cells);
    assert_eq!(thermal.time_s(), 0.0);
    let retry = joined::step_candidate(
        &material,
        &thermal,
        0.1,
        &model,
        &[0.0, 0.0, 0.04],
        0.3,
        &legs(),
        &BASES,
        config(),
        &[0.0; 4],
        || false,
    )
    .unwrap();
    assert_eq!(retry.model(), full.model());
    assert_eq!(retry.thermal().cells(), full.thermal().cells());
}
