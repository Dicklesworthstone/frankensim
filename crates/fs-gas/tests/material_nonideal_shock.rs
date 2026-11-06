//! G0/G1/G3/G4: original real-gas shock compression and conservative admission.
#![cfg(feature = "material")]

// Reuse only declared source data. The dense oracle below independently solves
// both jump equations in 75-digit Decimal arithmetic from PR's p(T,v), residual
// energy and entropy formulas; it does not call the production density solver.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use fs_gas::material::nonideal_shock::{
    RealGasNormalShock, RealGasPistonCompression, RealGasShockControl, RealGasShockError as E,
    normal_shock, piston_compression,
};
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Velocity};

fn original() -> ThermochemicalPengRobinsonModel {
    common::model([0.0; 2], [0.0; 2], [200.0, 2000.0])
}

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

fn controls() -> RealGasShockControl {
    RealGasShockControl {
        energy: GasEnergySolve {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 1000.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        compression_ratio_bracket: [1.05, 3.1],
        mass_flux_tolerance_kg_m2_s: 1e-8,
        momentum_tolerance_pa: 1e-4,
        jump_energy_tolerance_j_kg: 1e-5,
        entropy_tolerance_j_kg_k: 1e-8,
        maximum_density_iterations: 128,
    }
}

fn shock(
    model: &ThermochemicalPengRobinsonModel,
    density: f64,
    control: RealGasShockControl,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasNormalShock, E> {
    let context = point(400.0);
    let before = model
        .evaluate_density_temperature(&context, density)
        .unwrap();
    normal_shock(
        model,
        &context,
        Density::new(density),
        Velocity::new(2.0 * before.sound_speed_m_s()),
        control,
        cancelled,
    )
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

fn piston(
    model: &ThermochemicalPengRobinsonModel,
    density: f64,
    upstream_velocity: f64,
    piston_velocity: f64,
    control: RealGasShockControl,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasPistonCompression, E> {
    piston_compression(
        model,
        &point(400.0),
        Density::new(density),
        Velocity::new(upstream_velocity),
        Velocity::new(piston_velocity),
        control,
        cancelled,
    )
}

#[test]
fn dense_compression_matches_independent_pr_hugoniot_and_all_three_flux_balances() {
    let model = original();
    let before = model.clone();
    let result = shock(&model, 10.0, controls(), || false).unwrap();
    let inlet = result.upstream();
    let outlet = result.downstream();
    near(inlet.pressure_pa(), 1_405_488.503_926_915, 1e-6);
    near(inlet.sound_speed_m_s(), 434.547_071_783_336_5, 1e-10);
    near(outlet.density_kg_m3(), 27.210_526_175_865_127, 3e-9);
    near(outlet.temperature_k(), 632.774_712_247_970_2, 5e-8);
    near(outlet.pressure_pa(), 6_182_880.252_416_264, 5e-4);
    near(result.downstream_speed.value(), 319.396_301_986_079_2, 5e-8);
    near(result.downstream_mach, 0.568_258_525_953_416_6, 1e-10);
    near(
        result.internal_energy_increase_j_kg,
        239_980.326_490_340_5,
        1e-4,
    );
    near(result.entropy_change_j_kg_k, 118.500_402_140_820_74, 2e-7);
    assert_eq!(result.upstream_mach, 2.0);
    assert!(outlet.density_kg_m3() / inlet.density_kg_m3() < 2.75);
    // Calorically perfect gamma=4/3 would give rho2/rho1=2.8. The same EOS
    // supplies pressure, energy and sound speed, so this departure survives.
    assert!((outlet.density_kg_m3() / inlet.density_kg_m3() - 2.8).abs() > 0.07);
    assert_eq!(model, before);
    assert_eq!(outlet.constitutive_law(), model.mechanics());
    for component in outlet.source().components() {
        assert_eq!(component.query().axes().get("lot"), Some(&23.0));
        near(
            component.query().axes()["pressure"],
            outlet.pressure_pa(),
            0.0,
        );
        near(
            component.query().axes()["temperature"],
            outlet.temperature_k(),
            0.0,
        );
    }
    let u1 = result.upstream_speed.value();
    let u2 = result.downstream_speed.value();
    let rho1 = inlet.density_kg_m3();
    let rho2 = outlet.density_kg_m3();
    near(
        rho2 * u2 - rho1 * u1,
        0.0,
        controls().mass_flux_tolerance_kg_m2_s,
    );
    near(
        (outlet.pressure_pa() + rho2 * u2 * u2) - (inlet.pressure_pa() + rho1 * u1 * u1),
        0.0,
        controls().momentum_tolerance_pa,
    );
    // The authored zero datum also permits this independent absolute h check.
    near(
        outlet.specific_enthalpy_j_kg() + 0.5 * u2 * u2
            - inlet.specific_enthalpy_j_kg()
            - 0.5 * u1 * u1,
        0.0,
        controls().jump_energy_tolerance_j_kg,
    );
    assert!(result.hugoniot_residual_j_kg.abs() <= controls().jump_energy_tolerance_j_kg);
    assert!(
        result.stagnation_enthalpy_residual_j_kg.abs() <= controls().jump_energy_tolerance_j_kg
    );
    assert!(result.density_iterations <= controls().maximum_density_iterations);
    assert!(result.caloric_iterations > result.density_iterations);
}

#[test]
fn dilute_original_law_converges_to_the_analytic_normal_shock_limit() {
    let model = original();
    let gamma = 4.0 / 3.0;
    let mach_squared = 4.0;
    let density_ratio = (gamma + 1.0) * mach_squared / ((gamma - 1.0) * mach_squared + 2.0);
    let pressure_ratio = (2.0 * gamma * mach_squared - (gamma - 1.0)) / (gamma + 1.0);
    let temperature = 400.0 * pressure_ratio / density_ratio;
    let mut errors = Vec::new();
    for density in [0.08, 0.04, 0.02] {
        let result = shock(&model, density, controls(), || false).unwrap();
        let r = result.downstream().density_kg_m3() / density;
        let p = result.downstream().pressure_pa() / result.upstream().pressure_pa();
        let error = (r / density_ratio - 1.0).abs()
            + (p / pressure_ratio - 1.0).abs()
            + (result.downstream().temperature_k() / temperature - 1.0).abs();
        errors.push(error);
        assert!(result.entropy_change_j_kg_k > 0.0);
        assert!(result.downstream_mach < 1.0);
        near(
            result.mass_flux_residual_kg_m2_s,
            0.0,
            controls().mass_flux_tolerance_kg_m2_s,
        );
    }
    assert!(errors[1] < 0.6 * errors[0]);
    assert!(errors[2] < 0.6 * errors[1]);
    assert!(errors[2] < 1.5e-4);
}

#[test]
fn independent_species_energy_and_entropy_datums_leave_the_jump_unchanged() {
    let plain = original();
    let shifted = common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]);
    let a = shock(&plain, 10.0, controls(), || false).unwrap();
    let b = shock(&shifted, 10.0, controls(), || false).unwrap();
    assert_eq!(
        a.downstream().temperature_k(),
        b.downstream().temperature_k()
    );
    assert_eq!(
        a.downstream().density_kg_m3(),
        b.downstream().density_kg_m3()
    );
    assert_eq!(a.downstream().pressure_pa(), b.downstream().pressure_pa());
    assert_eq!(a.downstream_speed, b.downstream_speed);
    assert_eq!(a.downstream_mach, b.downstream_mach);
    assert_eq!(
        a.internal_energy_increase_j_kg,
        b.internal_energy_increase_j_kg
    );
    assert_eq!(a.entropy_change_j_kg_k, b.entropy_change_j_kg_k);
    assert_eq!(a.hugoniot_residual_j_kg, b.hugoniot_residual_j_kg);
    assert_eq!(
        a.stagnation_enthalpy_residual_j_kg,
        b.stagnation_enthalpy_residual_j_kg
    );
    assert_eq!(a.density_iterations, b.density_iterations);
    assert_eq!(a.caloric_iterations, b.caloric_iterations);
}

#[test]
fn identity_subsonic_unbracketed_and_unsupported_shocks_refuse() {
    let model = original();
    let context = point(400.0);
    let upstream = model.evaluate_density_temperature(&context, 10.0).unwrap();
    for speed in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        0.9 * upstream.sound_speed_m_s(),
        upstream.sound_speed_m_s(),
    ] {
        assert!(matches!(
            normal_shock(
                &model,
                &context,
                Density::new(10.0),
                Velocity::new(speed),
                controls(),
                || false
            ),
            Err(E::Invalid(_))
        ));
    }
    for bracket in [[1.0, 3.1], [0.9, 3.1], [1.1, 1.1], [1.1, f64::INFINITY]] {
        assert!(matches!(
            shock(
                &model,
                10.0,
                RealGasShockControl {
                    compression_ratio_bracket: bracket,
                    ..controls()
                },
                || false
            ),
            Err(E::Invalid(_))
        ));
    }
    for bracket in [[1.05, 1.1], [3.0, 3.1]] {
        assert!(matches!(
            shock(
                &model,
                10.0,
                RealGasShockControl {
                    compression_ratio_bracket: bracket,
                    ..controls()
                },
                || false
            ),
            Err(E::OutsideBracket { .. })
        ));
    }
    assert!(matches!(
        shock(
            &model,
            10.0,
            RealGasShockControl {
                maximum_density_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(E::NotConverged { .. })
    ));
    let narrow = RealGasShockControl {
        energy: GasEnergySolve {
            maximum_temperature_k: 450.0,
            ..controls().energy
        },
        ..controls()
    };
    assert!(matches!(
        shock(&model, 10.0, narrow, || false),
        Err(E::Energy(GasEnergyError::OutsideBracket { .. }))
    ));
    // The authored transport source ends at 10 MPa. A mechanical/caloric
    // candidate beyond that source support must not reuse upstream transport.
    assert!(matches!(
        shock(&model, 20.0, controls(), || false),
        Err(E::Energy(GasEnergyError::Material(_)))
    ));
}

#[test]
fn cancellation_never_returns_a_partial_jump_and_retry_is_deterministic() {
    let model = original();
    let before = model.clone();
    let mut polls = 0;
    let expected = shock(&model, 10.0, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for cancelled_at in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            shock(&model, 10.0, controls(), || {
                count += 1;
                count == cancelled_at
            }),
            Err(E::Cancelled)
        );
        assert_eq!(model, before);
    }
    assert_eq!(shock(&model, 10.0, controls(), || false).unwrap(), expected);
}

#[test]
fn prescribed_piston_matches_independent_dense_pr_oracle_and_moving_front_fluxes() {
    let model = original();
    let result = piston(&model, 10.0, 0.0, 500.0, controls(), || false).unwrap();
    let jump = &result.shock;
    let inlet = jump.upstream();
    let outlet = jump.downstream();
    // An independent 75-digit Decimal solve used prescribed Delta=500 m/s,
    // rho1=10 kg/m3, T1=400 K, explicit PR pressure/internal-energy/entropy,
    // caloric Newton and 220 density bisections. It did not invert the existing
    // Mach-2 shock oracle or call production code.
    near(outlet.density_kg_m3(), 25.512_202_655_197_167, 3e-9);
    near(outlet.temperature_k(), 604.703_814_348_865_1, 5e-8);
    near(outlet.pressure_pa(), 5_517_122.941_388_578, 5e-4);
    near(
        result.shock_normal_velocity.value(),
        822.326_887_492_332_6,
        5e-8,
    );
    near(jump.downstream_speed.value(), 322.326_887_492_332_6, 5e-8);
    near(jump.upstream_mach, 1.892_377_007_898_333_2, 1e-10);
    near(jump.downstream_mach, 0.588_895_735_103_771_7, 1e-10);
    near(
        jump.internal_energy_increase_j_kg,
        210_458.017_079_614_2,
        1e-4,
    );
    near(jump.entropy_change_j_kg_k, 93.873_213_385_056_56, 2e-7);
    assert_eq!(result.upstream_normal_velocity.value(), 0.0);
    assert_eq!(result.piston_normal_velocity.value(), 500.0);
    let front = result.shock_normal_velocity.value();
    let upstream_relative = result.upstream_normal_velocity.value() - front;
    let downstream_relative = result.piston_normal_velocity.value() - front;
    let rho1 = inlet.density_kg_m3();
    let rho2 = outlet.density_kg_m3();
    // Along the declared laboratory normal, both gas velocities relative to
    // the moving front and the conserved mass flux are negative.
    assert!(upstream_relative < 0.0 && downstream_relative < 0.0);
    near(
        rho1 * upstream_relative,
        -jump.mass_flux_kg_m2_s,
        controls().mass_flux_tolerance_kg_m2_s,
    );
    near(
        rho2 * downstream_relative,
        -jump.mass_flux_kg_m2_s,
        controls().mass_flux_tolerance_kg_m2_s,
    );
    near(
        outlet.pressure_pa() + rho2 * downstream_relative.powi(2)
            - inlet.pressure_pa()
            - rho1 * upstream_relative.powi(2),
        0.0,
        controls().momentum_tolerance_pa,
    );
    near(
        outlet.specific_enthalpy_j_kg() + 0.5 * downstream_relative.powi(2)
            - inlet.specific_enthalpy_j_kg()
            - 0.5 * upstream_relative.powi(2),
        0.0,
        controls().jump_energy_tolerance_j_kg,
    );
    assert!(jump.hugoniot_residual_j_kg.abs() <= controls().jump_energy_tolerance_j_kg);
    assert!(jump.stagnation_enthalpy_residual_j_kg.abs() <= controls().jump_energy_tolerance_j_kg);
    assert_eq!(outlet.constitutive_law(), model.mechanics());
    for component in outlet.source().components() {
        assert_eq!(component.query().axes().get("lot"), Some(&23.0));
        assert_eq!(
            component.query().axes()["temperature"],
            outlet.temperature_k()
        );
        assert_eq!(component.query().axes()["pressure"], outlet.pressure_pa());
    }
}

#[test]
fn dilute_piston_converges_to_the_analytic_prescribed_velocity_relation() {
    let model = original();
    let gamma = 4.0 / 3.0;
    let ideal_sound = (gamma * 8.314_462_618_153_24 * 400.0 / common::MIXTURE_MASS).sqrt();
    // Delta/c1 = 2*(M_s^2-1)/[(gamma+1)*M_s], solved independently for M_s.
    let k = (gamma + 1.0) * 500.0 / (2.0 * ideal_sound);
    let mach = 0.5 * (k + (k * k + 4.0).sqrt());
    let m2 = mach * mach;
    let density_ratio = (gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0);
    let pressure_ratio = (2.0 * gamma * m2 - (gamma - 1.0)) / (gamma + 1.0);
    let temperature = 400.0 * pressure_ratio / density_ratio;
    let mut errors = Vec::new();
    for density in [0.08, 0.04, 0.02] {
        let result = piston(&model, density, 0.0, 500.0, controls(), || false).unwrap();
        let jump = &result.shock;
        let r = jump.downstream().density_kg_m3() / density;
        let p = jump.downstream().pressure_pa() / jump.upstream().pressure_pa();
        let error = (r / density_ratio - 1.0).abs()
            + (p / pressure_ratio - 1.0).abs()
            + (jump.downstream().temperature_k() / temperature - 1.0).abs();
        errors.push(error);
        assert!(jump.upstream_mach > 1.0 && jump.downstream_mach < 1.0);
        assert!(jump.entropy_change_j_kg_k > 0.0);
        assert!(result.shock_normal_velocity.value() > 500.0);
    }
    assert!(errors[1] < 0.6 * errors[0]);
    assert!(errors[2] < 0.6 * errors[1]);
    assert!(errors[2] < 1.5e-4);
}

#[test]
fn piston_galilean_shifts_preserve_the_entire_source_bound_jump() {
    let model = original();
    let baseline = piston(&model, 10.0, 0.0, 500.0, controls(), || false).unwrap();
    for shift in [-2000.0, 1e6] {
        let shifted = piston(&model, 10.0, shift, shift + 500.0, controls(), || false).unwrap();
        assert_eq!(shifted.shock, baseline.shock);
        assert_eq!(shifted.upstream_normal_velocity.value(), shift);
        assert_eq!(shifted.piston_normal_velocity.value(), shift + 500.0);
        near(
            shifted.shock_normal_velocity.value() - shift,
            baseline.shock_normal_velocity.value(),
            1e-8,
        );
    }
}

#[test]
fn piston_invalid_motion_unsupported_branches_and_collapsed_fronts_refuse() {
    let model = original();
    for (upstream, wall) in [
        (0.0, 0.0),
        (500.0, 0.0),
        (f64::NAN, 500.0),
        (0.0, f64::INFINITY),
        (-f64::MAX, f64::MAX),
    ] {
        assert!(matches!(
            piston(&model, 10.0, upstream, wall, controls(), || false),
            Err(E::Invalid(_))
        ));
    }
    for bracket in [[1.05, 1.1], [3.0, 3.1]] {
        assert!(matches!(
            piston(
                &model,
                10.0,
                0.0,
                500.0,
                RealGasShockControl {
                    compression_ratio_bracket: bracket,
                    ..controls()
                },
                || false,
            ),
            Err(E::OutsideBracket { .. })
        ));
    }
    assert!(matches!(
        piston(
            &model,
            10.0,
            0.0,
            500.0,
            RealGasShockControl {
                maximum_density_iterations: 1,
                ..controls()
            },
            || false,
        ),
        Err(E::NotConverged { .. })
    ));
    let narrow = RealGasShockControl {
        energy: GasEnergySolve {
            maximum_temperature_k: 450.0,
            ..controls().energy
        },
        ..controls()
    };
    assert!(matches!(
        piston(&model, 10.0, 0.0, 500.0, narrow, || false),
        Err(E::Energy(GasEnergyError::OutsideBracket { .. }))
    ));
    assert!(matches!(
        piston(&model, 20.0, 0.0, 500.0, controls(), || false),
        Err(E::Energy(GasEnergyError::Material(_)))
    ));
    // Delta=1024 remains exactly representable, and the physical front has
    // u2 about 314.4 m/s and T2 about 986.7 K. A common offset of 2^62 makes
    // adding that front separation round to the piston speed itself.
    let wall = 2.0_f64.powi(62);
    let large_frame = RealGasShockControl {
        compression_ratio_bracket: [1.1, 4.6],
        energy: GasEnergySolve {
            maximum_temperature_k: 1200.0,
            ..controls().energy
        },
        ..controls()
    };
    assert!(matches!(
        piston(&model, 0.04, wall - 1024.0, wall, large_frame, || false),
        Err(E::Invalid(
            "laboratory shock speed must be finite and representably exceed the piston velocity"
        ))
    ));
}

#[test]
fn piston_cancellation_is_transactional_through_final_laboratory_admission() {
    let model = original();
    let before = model.clone();
    let mut polls = 0;
    let expected = piston(&model, 10.0, 0.0, 500.0, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for cancelled_at in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            piston(&model, 10.0, 0.0, 500.0, controls(), || {
                count += 1;
                count == cancelled_at
            }),
            Err(E::Cancelled)
        );
        assert_eq!(model, before);
    }
    assert_eq!(
        piston(&model, 10.0, 0.0, 500.0, controls(), || false).unwrap(),
        expected
    );
}
