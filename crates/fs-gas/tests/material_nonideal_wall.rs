//! Physical real-gas wall pressure, independent jump/fan anchors and field use.
#![cfg(all(feature = "material", feature = "experimental-euler"))]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use fs_gas::dynamics::nonideal::{
    RealGas, RealGasBoundaries, RealGasBoundary, RealGasDomain, RealGasEuler1d, RealGasFlowError,
    RealGasWallControl,
};
use fs_gas::dynamics::{Conserved, Reconstruction, StepControl};
use fs_gas::material::nonideal_shock::RealGasShockControl;
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Length, QuantitySpec, Temperature, Time, Velocity};

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

fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 300.0,
        maximum_temperature_k: 1000.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn gas(datums: [f64; 2]) -> RealGas {
    RealGas::new(
        common::model(datums, datums, [200.0, 2000.0]),
        point(800.0),
        Density::new(10.0),
        RealGasDomain {
            density_bounds_kg_m3: [0.005, 50.0],
            energy: energy(),
        },
        || false,
    )
    .unwrap()
}

fn controls() -> RealGasWallControl {
    RealGasWallControl {
        compression: RealGasShockControl {
            energy: energy(),
            compression_ratio_bracket: [1.0 + 1e-8, 3.1],
            mass_flux_tolerance_kg_m2_s: 1e-8,
            momentum_tolerance_pa: 1e-4,
            jump_energy_tolerance_j_kg: 1e-5,
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_density_iterations: 128,
        },
        minimum_rarefaction_density_kg_m3: 2.0,
        velocity_tolerance_m_s: 2e-6,
        integral_tolerance_m_s: 1e-7,
        maximum_quadrature_depth: 10,
    }
}

fn state(gas: &RealGas, velocity: f64, temperature: f64) -> Conserved {
    gas.conserved(
        Density::new(10.0),
        Velocity::new(velocity),
        Temperature::new(temperature),
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

#[test]
fn compression_pressure_matches_independent_dense_pr_hugoniot_at_either_wall() {
    let gas = gas([0.0; 2]);
    // Existing independent 75-digit PR normal-shock anchor: subtract its two
    // shock-frame velocities to put the downstream gas at a stationary wall.
    let speed = 869.094_143_566_673 - 319.396_301_986_079_2;
    let right = gas
        .reflecting_wall_flux(state(&gas, speed, 400.0), 1.0, controls(), || false)
        .unwrap();
    let left = gas
        .reflecting_wall_flux(state(&gas, -speed, 400.0), -1.0, controls(), || false)
        .unwrap();
    near(right.momentum, 6_182_880.252_416_264, 8e-4);
    assert_eq!(left, right);
    assert_eq!(right.mass, 0.0);
    assert_eq!(right.energy, 0.0);
    assert!(!right.used_hll_fallback);
    assert!(right.signal_speed >= speed + 434.547_071_783_336_5);
}

#[test]
fn dense_rarefaction_has_positive_native_wall_pressure_when_mirror_flux_is_tensile() {
    let gas = gas([0.0; 2]);
    let inside = state(&gas, -600.0, 800.0);
    let exterior = state(&gas, 600.0, 800.0);
    let mirror = gas.hllc_flux(inside, exterior, || false).unwrap();
    assert!(mirror.momentum < 0.0);
    let physical = gas
        .reflecting_wall_flux(inside, 1.0, controls(), || false)
        .unwrap();
    // Independent PR calculation: molar s=3R ln(T)+R ln(v-b)
    // + a'(T)/(2sqrt(2)b) ln[(v+(1+sqrt(2))b)/(v+(1-sqrt(2))b)].
    // Solve its original isentrope and integrate a dln(rho) with scipy's
    // independent adaptive Gauss--Kronrod quadrature, then a Brent velocity root.
    // rho_wall=3.4938816520658844, T_wall=560.3754865199095 K.
    near(physical.momentum, 704_086.609_244_884_3, 0.02);
    assert_eq!(physical.mass, 0.0);
    assert_eq!(physical.energy, 0.0);
    assert!(!physical.used_hll_fallback);
    near(physical.signal_speed, 1226.318_447_565_742, 1e-8);
    let reflected = gas
        .reflecting_wall_flux(exterior, -1.0, controls(), || false)
        .unwrap();
    assert_eq!(physical, reflected);
}

#[test]
fn wall_pressure_is_invariant_to_large_original_energy_and_entropy_datums() {
    let original = gas([0.0; 2]);
    let shifted = gas([2e17, -7e16]);
    for (velocity, temperature) in [(-600.0, 800.0), (549.697_841_580_593_8, 400.0)] {
        let expected = original
            .reflecting_wall_flux(
                state(&original, velocity, temperature),
                1.0,
                controls(),
                || false,
            )
            .unwrap();
        let actual = shifted
            .reflecting_wall_flux(
                state(&shifted, velocity, temperature),
                1.0,
                controls(),
                || false,
            )
            .unwrap();
        near(actual.momentum, expected.momentum, 2e-3);
        near(actual.signal_speed, expected.signal_speed, 1e-7);
    }
}

#[test]
fn equilibrium_and_weak_waves_have_the_original_acoustic_impedance() {
    let gas = gas([0.0; 2]);
    let rest = state(&gas, 0.0, 800.0);
    let p = gas.primitive(rest, || false).unwrap();
    let flux = gas
        .reflecting_wall_flux(rest, 1.0, controls(), || false)
        .unwrap();
    assert_eq!(flux.momentum, p.pressure);
    assert_eq!(flux.signal_speed, p.sound_speed);
    for speed in [-1e-6, 1e-6] {
        let flux = gas
            .reflecting_wall_flux(state(&gas, speed, 800.0), 1.0, controls(), || false)
            .unwrap();
        near(
            flux.momentum - p.pressure,
            p.density * p.sound_speed * speed,
            1e-8,
        );
    }
}

#[test]
fn nonzero_acoustic_waves_cannot_leave_the_original_temperature_path() {
    let gas = gas([0.0; 2]);
    for (temperature, outward_speed) in [(1000.0, 1e-6), (300.0, -1e-6)] {
        // A tolerance-sized velocity is still a nonzero thermodynamic wave.
        // Only the direction INTO the declared source interval is supported.
        assert!(
            gas.reflecting_wall_flux(
                state(&gas, outward_speed, temperature),
                1.0,
                controls(),
                || false,
            )
            .is_err()
        );
        gas.reflecting_wall_flux(
            state(&gas, -outward_speed, temperature),
            1.0,
            controls(),
            || false,
        )
        .unwrap();
        let rest = state(&gas, 0.0, temperature);
        let pressure = gas.primitive(rest, || false).unwrap().pressure;
        assert_eq!(
            gas.reflecting_wall_flux(rest, 1.0, controls(), || false)
                .unwrap()
                .momentum,
            pressure,
        );
    }
}

#[test]
fn acoustic_and_resolved_wall_waves_join_without_a_compression_gap() {
    let gas = gas([0.0; 2]);
    let primitive = gas.primitive(state(&gas, 0.0, 800.0), || false).unwrap();
    let mut wall = controls();
    wall.compression.compression_ratio_bracket[1] = 1.2;
    let resolution = 32.0 * f64::EPSILON.sqrt() * primitive.sound_speed;
    let mut speeds = vec![
        1e-6,
        1.99e-6,
        2e-6,
        2.01e-6,
        3e-6,
        6e-6,
        0.9 * resolution,
        0.99 * resolution,
        resolution,
        1.01 * resolution,
        1.1 * resolution,
        2.0 * resolution,
        10.0 * resolution,
    ];
    speeds.sort_by(f64::total_cmp);
    for normal in [-1.0, 1.0] {
        for sign in [-1.0, 1.0] {
            let mut previous = primitive.pressure;
            for &magnitude in &speeds {
                let speed = sign * magnitude;
                let flux = gas
                    .reflecting_wall_flux(state(&gas, normal * speed, 800.0), normal, wall, || {
                        false
                    })
                    .unwrap();
                assert!(sign * (flux.momentum - previous) > 0.0);
                // Independent weak-wave impedance, allowing its quadratic
                // nonlinear pressure term and the authored shock residual.
                near(
                    flux.momentum,
                    primitive.pressure + primitive.density * primitive.sound_speed * speed,
                    5e-4 + 4.0 * primitive.density * speed * speed,
                );
                assert_eq!(flux.mass, 0.0);
                assert_eq!(flux.energy, 0.0);
                previous = flux.momentum;
            }
        }
    }
    // A caller's finite compression bracket must not leave a hole between
    // equilibrium and its lower endpoint when a weak shock is resolved.
    wall.compression.compression_ratio_bracket[0] = 1.05;
    let flux = gas
        .reflecting_wall_flux(state(&gas, 0.002, 800.0), 1.0, wall, || false)
        .unwrap();
    near(
        flux.momentum,
        primitive.pressure + primitive.density * primitive.sound_speed * 0.002,
        1e-3,
    );
}

#[test]
fn physical_wall_field_uses_pressure_traction_and_conserves_closed_mass_and_energy() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, 30.0, 800.0);
    let mut wall = controls();
    wall.compression.compression_ratio_bracket[1] = 1.2;
    let mut field = RealGasEuler1d::new(
        gas,
        Length::new(0.02),
        0.5,
        vec![q; 8],
        RealGasBoundaries::Ends {
            left: RealGasBoundary::PhysicalReflecting(wall),
            right: RealGasBoundary::PhysicalReflecting(wall),
        },
        || false,
    )
    .unwrap();
    let before = field.totals(|| false).unwrap();
    let step = StepControl {
        cfl: 0.4,
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 4,
    };
    for _ in 0..4 {
        let report = field.step(Time::new(1e-7), step, || false).unwrap();
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(report.boundary_change.momentum < 0.0);
        near(report.conservation_residual.momentum, 0.0, 2e-12);
    }
    let after = field.totals(|| false).unwrap();
    near(after.mass, before.mass, 2e-13);
    near(after.energy, before.energy, 1e-9);
    assert!(after.momentum < before.momentum);
    assert!(field.primitives()[0].pressure < field.primitives()[7].pressure);
}

#[test]
fn unsupported_wave_paths_and_mid_solve_cancellation_return_no_wall_flux() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, -600.0, 800.0);
    let mut short = controls();
    short.minimum_rarefaction_density_kg_m3 = 8.0;
    assert!(gas.reflecting_wall_flux(q, 1.0, short, || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        gas.reflecting_wall_flux(q, 1.0, controls(), || {
            calls += 1;
            calls == 20
        }),
        Err(RealGasFlowError::Cancelled)
    ));
    // Source pressure and accepted conservative data have not been modified.
    assert_eq!(gas.primitive(q, || false).unwrap().velocity, -600.0);
}
