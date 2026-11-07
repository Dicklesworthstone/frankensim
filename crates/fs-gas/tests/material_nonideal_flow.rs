//! Original-source real-gas fields: compression, wave propagation and atomicity.
#![cfg(all(feature = "material", feature = "experimental-euler"))]

// Source data only. The shock and acoustic numbers below come from the
// independent 75-digit original PR oracle used by material_nonideal_shock;
// no production shock solution is used as this field's expected answer.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use fs_gas::dynamics::nonideal::{
    RealGas, RealGasBoundaries, RealGasBoundary, RealGasDomain, RealGasEuler1d, RealGasFlowError,
};
use fs_gas::dynamics::{Conserved, Mesh1d, Reconstruction, StepControl, StepReport};
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Length, QuantitySpec, Temperature, Time, Velocity};

const DENSE_SOUND_SPEED: f64 = 434.547_071_783_336_5;

#[path = "material_nonideal_flow/forcing.rs"]
mod forcing;
#[path = "material_nonideal_flow/nozzle.rs"]
mod nozzle;
#[path = "material_nonideal_flow/reconstruction.rs"]
mod reconstruction;
#[path = "material_nonideal_flow/thermal.rs"]
mod thermal;

fn point(temperature: f64) -> QueryPoint {
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 1e5)
        .unwrap()
        .with("lot", 23.0)
        .unwrap()
}

fn domain() -> RealGasDomain {
    RealGasDomain {
        density_bounds_kg_m3: [0.005, 50.0],
        energy: GasEnergySolve {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 1000.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
    }
}

fn bind(model: ThermochemicalPengRobinsonModel, domain: RealGasDomain) -> RealGas {
    RealGas::new(model, point(400.0), Density::new(10.0), domain, || false).unwrap()
}

fn gas() -> RealGas {
    bind(common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]), domain())
}

fn state(gas: &RealGas, density: f64, velocity: f64, temperature: f64) -> Conserved {
    gas.conserved(
        Density::new(density),
        Velocity::new(velocity),
        Temperature::new(temperature),
        || false,
    )
    .unwrap()
}

fn controls() -> StepControl {
    StepControl {
        cfl: 0.4,
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 12,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

fn check_balance(report: StepReport) {
    assert!(report.maximum_cfl <= controls().cfl * (1.0 + 64.0 * f64::EPSILON));
    for (before, after, left, right, residual) in [
        (
            report.before.mass,
            report.after.mass,
            report.left_boundary_change.mass,
            report.right_boundary_change.mass,
            report.conservation_residual.mass,
        ),
        (
            report.before.momentum,
            report.after.momentum,
            report.left_boundary_change.momentum,
            report.right_boundary_change.momentum,
            report.conservation_residual.momentum,
        ),
        (
            report.before.energy,
            report.after.energy,
            report.left_boundary_change.energy,
            report.right_boundary_change.energy,
            report.conservation_residual.energy,
        ),
    ] {
        let scale = 1.0 + before.abs() + after.abs() + left.abs() + right.abs();
        near(residual, 0.0, 512.0 * f64::EPSILON * scale);
    }
}

fn run_until(flow: &mut RealGasEuler1d, finish: f64) {
    let mut steps = 0;
    while flow.time().value() < finish {
        let report = flow
            .step(Time::new(finish - flow.time().value()), controls(), || {
                false
            })
            .unwrap();
        check_balance(report);
        steps += 1;
        assert!(steps <= 4096, "bounded field integration did not finish");
    }
    assert_eq!(flow.time().value(), finish);
}

fn uniform_mesh_field(gas: RealGas, states: Vec<Conserved>) -> RealGasEuler1d {
    RealGasEuler1d::new(
        gas,
        Length::new(1.0 / states.len() as f64),
        0.07,
        states,
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap()
}

#[test]
fn one_common_original_energy_reference_retains_density_calorics_and_source_axes() {
    let gas = gas();
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    assert!((gas.energy_reference() - shifted.energy_reference()).abs() > 1e20);
    for (rho, velocity, temperature) in [
        // This velocity creates a tiny nonzero kinetic subtraction remainder
        // at the exact reference. The retained conservative value is unchanged.
        (10.0, 11.011, 400.0),
        (5.0, -90.0, 350.0),
        (15.0, 120.0, 450.0),
        (20.0, 0.0, 400.0),
    ] {
        let q = state(&gas, rho, velocity, temperature);
        assert_eq!(q, state(&shifted, rho, velocity, temperature));
        let p = gas.primitive(q, || false).unwrap();
        assert_eq!(p, shifted.primitive(q, || false).unwrap());
        let original = gas
            .model()
            .evaluate_density_temperature(&point(temperature), rho)
            .unwrap();
        near(p.density, rho, 0.0);
        near(p.velocity, velocity, 2e-14 * velocity.abs().max(1.0));
        near(p.temperature, temperature, 2e-9);
        near(p.pressure, original.pressure_pa(), 3e-5);
        near(p.sound_speed, original.sound_speed_m_s(), 3e-9);
        near(
            gas.absolute_energy_density(q).unwrap(),
            rho * (original.specific_internal_energy_j_kg() + 0.5 * velocity * velocity),
            3e-7,
        );
        let resolved = gas.thermodynamics(q, || false).unwrap();
        assert_eq!(resolved.constitutive_law(), gas.model().mechanics());
        for component in resolved.source().components() {
            assert_eq!(component.query().axes().get("lot"), Some(&23.0));
            assert_eq!(component.query().axes()["pressure"], resolved.pressure_pa());
            assert_eq!(
                component.query().axes()["temperature"],
                resolved.temperature_k()
            );
        }
    }
    // A density-only change at fixed T has nonzero original residual energy.
    // A fresh e(T_ref,rho_cell) reference per cell would incorrectly erase it.
    let compressed = state(&gas, 20.0, 0.0, 400.0);
    assert!((compressed.energy / compressed.density).abs() > 100.0);
    let q = state(&gas, 10.0, 0.0, 400.0);
    let p = gas.primitive(q, || false).unwrap();
    near(p.pressure, 1_405_488.503_926_915, 1e-6);
    near(p.sound_speed, DENSE_SOUND_SPEED, 1e-10);
}

#[test]
fn source_datum_shifts_leave_complete_compressive_field_history_unchanged() {
    let plain = gas();
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    let states = (0..12)
        .map(|i| state(&plain, 10.0, if i < 6 { 45.0 } else { -45.0 }, 400.0))
        .collect::<Vec<_>>();
    let mut a = uniform_mesh_field(plain, states.clone());
    let mut b = uniform_mesh_field(shifted, states);
    for _ in 0..5 {
        let ra = a.step(Time::new(2e-6), controls(), || false).unwrap();
        let rb = b.step(Time::new(2e-6), controls(), || false).unwrap();
        check_balance(ra);
        assert_eq!(ra, rb);
        assert_eq!(a.states(), b.states());
        assert_eq!(a.primitives(), b.primitives());
        assert_eq!(a.time(), b.time());
    }
}

#[test]
fn nonuniform_mesh_preserves_free_stream_and_closed_compression_balances() {
    let gas = gas();
    let positions = vec![0.0, 0.08, 0.25, 0.5, 0.75, 0.92, 1.0];
    let volumes = positions.windows(2).map(|x| 0.07 * (x[1] - x[0])).collect();
    let mesh = Mesh1d::new(positions, vec![0.07; 7], volumes).unwrap();
    let q = state(&gas, 10.0, 11.011, 400.0);
    let mut uniform = RealGasEuler1d::new_on_mesh(
        gas.clone(),
        mesh.clone(),
        vec![q; 6],
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap();
    run_until(&mut uniform, 3e-5);
    assert!(uniform.states().iter().all(|&actual| actual == q));
    let initial = (0..6)
        .map(|i| state(&gas, 10.0, if i < 3 { 60.0 } else { -60.0 }, 400.0))
        .collect();
    let mut compression = RealGasEuler1d::new_on_mesh(
        gas,
        mesh,
        initial,
        RealGasBoundaries::Ends {
            left: RealGasBoundary::Reflecting,
            right: RealGasBoundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let before = compression.totals(|| false).unwrap();
    for _ in 0..8 {
        let report = compression
            .step(Time::new(2e-6), controls(), || false)
            .unwrap();
        check_balance(report);
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert_eq!(report.left_boundary_change.mass, 0.0);
        assert_eq!(report.right_boundary_change.energy, 0.0);
    }
    let after = compression.totals(|| false).unwrap();
    near(after.mass, before.mass, 3e-14);
    near(after.energy, before.energy, 3e-10);
    assert!(compression.primitives()[2].temperature > 400.0);
    assert!(compression.primitives()[2].pressure > 1_405_488.503_926_915);
    let internal_increase: f64 = compression
        .states()
        .iter()
        .zip(compression.mesh().volumes())
        .map(|(q, volume)| (q.energy - 0.5 * q.momentum * q.momentum / q.density) * volume)
        .sum();
    assert!(internal_increase > 0.0);
}

#[test]
fn dense_translating_shock_refines_toward_independent_real_gas_jump() {
    let mut errors = Vec::new();
    let finish = 1e-4;
    let speed = 150.0;
    let shock_position = 0.5 + speed * finish;
    for n in [40, 80, 160] {
        let gas = gas();
        let left = state(&gas, 10.0, 869.094_143_566_673 + speed, 400.0);
        let right = state(
            &gas,
            27.210_526_175_865_127,
            319.396_301_986_079_2 + speed,
            632.774_712_247_970_2,
        );
        let dx = 1.0 / n as f64;
        let states = (0..n)
            .map(|i| if i < n / 2 { left } else { right })
            .collect();
        let mut flow = RealGasEuler1d::new(
            gas,
            Length::new(dx),
            0.07,
            states,
            RealGasBoundaries::Ends {
                left: RealGasBoundary::Exterior(left),
                right: RealGasBoundary::Exterior(right),
            },
            || false,
        )
        .unwrap();
        run_until(&mut flow, finish);
        let error: f64 = flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                let left_fraction = ((shock_position - i as f64 * dx) / dx).clamp(0.0, 1.0);
                let exact = left_fraction * left.density + (1.0 - left_fraction) * right.density;
                (q.density - exact).abs() * dx
            })
            .sum();
        assert!(error.is_finite());
        assert!(flow.primitives().iter().all(|p| p.temperature > 300.0));
        errors.push(error);
    }
    assert!(errors[1] < 0.9 * errors[0], "shock L1 errors: {errors:?}");
    assert!(errors[2] < 0.9 * errors[1], "shock L1 errors: {errors:?}");
    assert!(errors[2] < 0.75, "shock L1 errors: {errors:?}");
}

#[test]
fn weak_compression_wave_refines_to_original_adiabatic_sound_speed() {
    let mut speed_errors = Vec::new();
    let finish = 0.1 / DENSE_SOUND_SPEED;
    for n in [24, 48, 96] {
        let gas = gas();
        let reference = gas.reference_state();
        let temperature_density_slope = reference.temperature_k()
            * reference.root().pressure_derivatives().0
            / (100.0 * reference.specific_heat_cv_j_kg_k());
        let k = core::f64::consts::TAU;
        let half_cell_phase = core::f64::consts::PI / n as f64;
        let average_factor = half_cell_phase.sin() / half_cell_phase;
        let states = (0..n)
            .map(|i| {
                let phase = k * (i as f64 + 0.5) / n as f64;
                let drho = 1e-4 * phase.sin() * average_factor;
                state(
                    &gas,
                    10.0 + drho,
                    DENSE_SOUND_SPEED * drho / 10.0,
                    400.0 + temperature_density_slope * drho,
                )
            })
            .collect();
        let mut flow = uniform_mesh_field(gas, states);
        run_until(&mut flow, finish);
        let (sine, cosine) =
            flow.states()
                .iter()
                .enumerate()
                .fold((0.0, 0.0), |(sine, cosine), (i, q)| {
                    let phase = k * (i as f64 + 0.5) / n as f64;
                    (
                        sine + (q.density - 10.0) * phase.sin(),
                        cosine + (q.density - 10.0) * phase.cos(),
                    )
                });
        assert!(sine > 0.0 && cosine < 0.0);
        let measured_speed = (-cosine).atan2(sine) / (k * finish);
        speed_errors.push((measured_speed / DENSE_SOUND_SPEED - 1.0).abs());
    }
    assert!(
        speed_errors[1] < 0.6 * speed_errors[0],
        "phase errors: {speed_errors:?}"
    );
    assert!(
        speed_errors[2] < 0.6 * speed_errors[1],
        "phase errors: {speed_errors:?}"
    );
    assert!(speed_errors[2] < 0.002, "phase errors: {speed_errors:?}");
}

#[test]
fn cancellation_source_refusal_and_boundary_edits_publish_no_partial_field() {
    let gas = gas();
    let states = (0..8)
        .map(|i| state(&gas, 10.0, if i < 4 { 50.0 } else { -50.0 }, 400.0))
        .collect::<Vec<_>>();
    let mut baseline = uniform_mesh_field(gas.clone(), states.clone());
    let mut polls = 0;
    let expected = baseline
        .step(Time::new(1e-5), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 20);
    for stop in [1, polls / 2, polls] {
        let mut flow = uniform_mesh_field(gas.clone(), states.clone());
        let primitive = flow.primitives().to_vec();
        let mut calls = 0;
        let result = flow.step(Time::new(1e-5), controls(), || {
            calls += 1;
            calls == stop
        });
        assert_eq!(result, Err(RealGasFlowError::Cancelled));
        assert_eq!(flow.states(), states);
        assert_eq!(flow.primitives(), primitive);
        assert_eq!(flow.time().value(), 0.0);
        let retry = flow.step(Time::new(1e-5), controls(), || false).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(flow.states(), baseline.states());
        assert_eq!(flow.primitives(), baseline.primitives());
    }
    let mut flow = uniform_mesh_field(gas, states.clone());
    let invalid = Conserved {
        density: 100.0,
        ..states[0]
    };
    assert!(
        flow.set_boundaries(
            RealGasBoundaries::Ends {
                left: RealGasBoundary::Exterior(invalid),
                right: RealGasBoundary::Transmissive,
            },
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.boundaries(), RealGasBoundaries::Periodic);
    assert!(
        flow.step(
            Time::new(1e-5),
            StepControl {
                cfl: 0.5,
                ..controls()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time().value(), 0.0);

    let mut narrow = domain();
    narrow.energy.minimum_temperature_k = 399.999;
    narrow.energy.maximum_temperature_k = 400.001;
    let gas = bind(common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]), narrow);
    let states = vec![
        state(&gas, 10.0, 500.0, 400.0),
        state(&gas, 10.0, -500.0, 400.0),
    ];
    let mut flow = uniform_mesh_field(gas, states.clone());
    let bounded = StepControl {
        maximum_halvings: 0,
        ..controls()
    };
    assert!(matches!(
        flow.step(Time::new(1e-5), bounded, || false),
        Err(RealGasFlowError::StepRetries { attempts: 1, .. })
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time().value(), 0.0);
}
