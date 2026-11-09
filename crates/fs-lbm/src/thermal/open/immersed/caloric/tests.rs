use super::*;
use crate::thermal::boundary::FiniteThermalWall;
use crate::thermal::open::immersed::ImmersedThermalBody3;
use crate::thermal::open::immersed::tests::{fixture, fluid_heat, near, step};
use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, HeatCapacityKnot, SolidLiquidPhase};

fn knot(h: f64, t: f64, liquid: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: 900.0,
    }
}
fn curve(latent: bool) -> Arc<EquilibriumEnthalpyPhaseCurve> {
    Arc::new(
        if latent {
            EquilibriumEnthalpyPhaseCurve::try_new(
                ContentHash([91; 32]),
                vec![
                    knot(0.0, 300.0, 0.0),
                    knot(40_000.0, 340.0, 0.0),
                    knot(80_000.0, 340.0, 1.0),
                    knot(160_000.0, 420.0, 1.0),
                ],
            )
        } else {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ContentHash([92; 32]),
                SolidLiquidPhase::Solid,
                vec![knot(0.0, 300.0, 0.0), knot(120_000.0, 420.0, 0.0)],
            )
        }
        .unwrap(),
    )
}
fn controls() -> ImmersedCaloricControls3 {
    ImmersedCaloricControls3 {
        maximum_iterations: 64,
        maximum_source_knots: 16,
        absolute_energy_tolerance_j: 1e-10,
    }
}
fn binding(latent: bool) -> ImmersedCaloricBinding3 {
    ImmersedCaloricBinding3 {
        mass_kg: 0.004,
        initial_specific_enthalpy_j_kg: if latent { 60_000.0 } else { 56_000.0 },
        specific_enthalpy_range_j_kg: [0.0, if latent { 160_000.0 } else { 120_000.0 }],
        temperature_origin_k: 300.0,
        temperature_scale_k: 40.0,
        fluid_cell_heat_capacity_j_k: 1.0,
    }
}
fn caloric_fixture(latent: bool, source: f64) -> ThermalOpenField3 {
    let mut field = fixture(&[27], source, false);
    if latent {
        field
            .set_immersed_thermal_bodies(
                &[ImmersedThermalBody3 {
                    solid_cells: &[27],
                    solid: FiniteThermalWall {
                        temperature: 1.0,
                        heat_capacity: 4.0,
                        resistance: 0.7,
                        heat_per_step: source,
                    },
                }],
                64 * 1024,
            )
            .unwrap();
    }
    field
        .set_immersed_body_caloric(0, curve(latent), binding(latent), controls())
        .unwrap();
    field
}

#[test]
fn source_sensible_storage_matches_original_constant_capacity_and_heat_accounts() {
    let mut material = caloric_fixture(false, 0.03);
    let mut original = fixture(&[27], 0.03, false);
    for _ in 0..8 {
        step(&mut material);
        step(&mut original);
        let actual = material.immersed_body(0).unwrap();
        let reference = original.immersed_body(0).unwrap();
        near(actual.solid.temperature, reference.solid.temperature, 1e-12);
        let caloric = actual.caloric.unwrap();
        assert_eq!(caloric.mass_kg, 0.004);
        near(
            (caloric.phase.temperature_k() - 300.0) / 40.0,
            actual.solid.temperature,
            0.0,
        );
        near(actual.last_exchange.unwrap().defect(), 0.0, 3e-12);
    }
    assert_eq!(
        material.immersed_body_caloric_curve(0),
        Some(curve(false).as_ref())
    );

    // The ORIGINAL exact Cp(T) constructor retains Cp=1000+10(T-300).
    // In native theta, h=40000*theta+8000*theta^2. Its full implicit
    // balance is a quadratic, so neither endpoint/frozen Cp can pass this.
    let nonlinear = Arc::new(
        EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
            ContentHash([93; 32]),
            SolidLiquidPhase::Solid,
            0.0,
            &[
                HeatCapacityKnot {
                    temperature_k: 300.0,
                    specific_heat_capacity_j_kg_k: 1000.0,
                    bulk_density_kg_m3: 900.0,
                },
                HeatCapacityKnot {
                    temperature_k: 420.0,
                    specific_heat_capacity_j_kg_k: 2200.0,
                    bulk_density_kg_m3: 900.0,
                },
            ],
        )
        .unwrap(),
    );
    let mut field = fixture(&[27], 0.0, false);
    field
        .set_immersed_thermal_bodies(
            &[ImmersedThermalBody3 {
                solid_cells: &[27],
                solid: FiniteThermalWall {
                    temperature: 1.0,
                    heat_capacity: 4.0,
                    resistance: 0.7,
                    heat_per_step: 0.0,
                },
            }],
            64 * 1024,
        )
        .unwrap();
    field
        .set_immersed_body_caloric(
            0,
            nonlinear,
            ImmersedCaloricBinding3 {
                initial_specific_enthalpy_j_kg: 48_000.0,
                specific_enthalpy_range_j_kg: [0.0, 192_000.0],
                ..binding(false)
            },
            controls(),
        )
        .unwrap();
    step(&mut field);
    let k: f64 = 6.0 / 4.7;
    let b = 4.0 + k;
    let c = 4.8 + 0.2 * k;
    let expected = 2.0 * c / (b + (b * b + 3.2 * c).sqrt());
    let actual = field.immersed_body(0).unwrap();
    near(actual.solid.temperature, expected, 1e-12);
    near(
        actual.caloric.unwrap().phase.specific_enthalpy_j_kg(),
        40_000.0 * expected + 8_000.0 * expected * expected,
        1e-7,
    );
    near(actual.last_exchange.unwrap().defect(), 0.0, 3e-12);
}

#[test]
fn original_reflections_freeze_and_melt_latent_store_without_temperature_drift() {
    for source in [0.0, 2.0] {
        let mut field = caloric_fixture(true, source);
        let before = fluid_heat(&field);
        let balance = step(&mut field);
        let body = field.immersed_body(0).unwrap();
        let material = body.caloric.unwrap();
        let independent_heat = (6.0 / 4.7) * (1.0 - 0.2);
        let expected_h = 60_000.0 + (source - independent_heat) / 0.0001;
        assert_eq!(body.solid.temperature, 1.0);
        assert_eq!(material.phase.temperature_k(), 340.0);
        near(material.phase.specific_enthalpy_j_kg(), expected_h, 2e-8);
        near(
            material.phase.liquid_mass_fraction(),
            (expected_h - 40_000.0) / 40_000.0,
            1e-12,
        );
        let exchange = body.last_exchange.unwrap();
        near(exchange.heat_in, independent_heat, 1e-14);
        near(
            fluid_heat(&field) - before + exchange.stored_change,
            source,
            1e-12,
        );
        near(balance.combined_defect, 0.0, 1e-12);
    }
    let mut field = caloric_fixture(true, 10.0);
    step(&mut field);
    let body = field.immersed_body(0).unwrap();
    assert!(body.caloric.unwrap().phase.specific_enthalpy_j_kg() > 80_000.0);
    assert_eq!(body.caloric.unwrap().phase.liquid_mass_fraction(), 1.0);
    assert!(body.solid.temperature > 1.0);
    near(body.last_exchange.unwrap().defect(), 0.0, 3e-12);
}

#[test]
fn source_domain_and_late_cancellation_preserve_native_and_caloric_state_for_exact_retry() {
    let mut field = caloric_fixture(true, 1000.0);
    let old = field.immersed_body(0);
    let populations = field.populations().to_vec();
    let clock = field.completed_steps();
    assert!(matches!(
        field.begin_step::<ThermalOpenError3>(&[0.0; 60], None, || Ok(())),
        Err(ThermalOpenError3::Caloric {
            source: ImmersedCaloricError3::OutsideEnthalpyRange { .. },
            ..
        })
    ));
    assert_eq!(field.immersed_body(0), old);
    assert_eq!(field.populations(), populations);
    assert_eq!(field.completed_steps(), clock);

    field.set_immersed_body_heat(0, 0.0).unwrap();
    let old = field.immersed_body(0);
    let candidate = field
        .begin_step::<ThermalOpenError3>(&[0.0; 60], None, || Ok(()))
        .unwrap();
    assert_eq!(field.immersed_body(0), old);
    assert!(matches!(
        field.finish_step::<ThermalOpenError3>(&[0.0; 60], candidate, || Err(
            ThermalOpenError3::Cancelled
        )),
        Err(ThermalOpenError3::Cancelled)
    ));
    assert_eq!(field.immersed_body(0), old);
    assert_eq!(field.populations(), populations);
    assert_eq!(field.completed_steps(), clock);
    let mut reference = caloric_fixture(true, 0.0);
    assert_eq!(step(&mut field), step(&mut reference));
    assert_eq!(field.populations(), reference.populations());
    assert_eq!(field.immersed_body(0), reference.immersed_body(0));
}

#[test]
fn binding_refuses_mismatched_storage_unsupported_graph_and_post_step_replacement() {
    let mut field = fixture(&[27], 0.0, false);
    let old = field.immersed_body(0);
    let mut wrong = binding(false);
    wrong.mass_kg = -1.0;
    assert!(
        field
            .set_immersed_body_caloric(0, curve(false), wrong, controls())
            .is_err()
    );
    wrong = binding(false);
    wrong.initial_specific_enthalpy_j_kg = 40_000.0;
    assert!(
        field
            .set_immersed_body_caloric(0, curve(false), wrong, controls())
            .is_err()
    );
    assert_eq!(field.immersed_body(0), old);
    let solid = old.unwrap().solid;
    field
        .set_immersed_thermal_network(
            &[ImmersedThermalBody3 {
                solid_cells: &[27],
                solid,
            }],
            &[],
            64 * 1024,
            crate::thermal::boundary::conduction::WallConductionControls::default(),
        )
        .unwrap();
    let old = field.immersed_body(0);
    assert!(
        field
            .set_immersed_body_caloric(0, curve(false), binding(false), controls())
            .is_err()
    );
    assert_eq!(field.immersed_body(0), old);
    let mut field = caloric_fixture(false, 0.0);
    assert!(
        field
            .set_immersed_body_caloric(0, curve(false), binding(false), controls())
            .is_err()
    );
    let mut field = fixture(&[27], 0.0, false);
    step(&mut field);
    let old = field.immersed_body(0);
    assert!(
        field
            .set_immersed_body_caloric(0, curve(false), binding(false), controls())
            .is_err()
    );
    assert_eq!(field.immersed_body(0), old);

    let mut field = fixture(&[27, 28], 0.0, false);
    let solid = field.immersed_body(0).unwrap().solid;
    field
        .set_immersed_thermal_bodies(
            &[
                ImmersedThermalBody3 {
                    solid_cells: &[27],
                    solid,
                },
                ImmersedThermalBody3 {
                    solid_cells: &[28],
                    solid,
                },
            ],
            64 * 1024,
        )
        .unwrap();
    field
        .set_immersed_body_caloric(0, curve(false), binding(false), controls())
        .unwrap();
    let old = [field.immersed_body(0), field.immersed_body(1)];
    let populations = field.populations().to_vec();
    // Each alternate map could reproduce this body's original native T and
    // fit its source/thermal bounds by itself. Sharing one fluid gives it a
    // different physical meaning from the first body's already bound map.
    for mismatch in [
        ImmersedCaloricBinding3 {
            temperature_origin_k: 301.0,
            initial_specific_enthalpy_j_kg: 57_000.0,
            ..binding(false)
        },
        ImmersedCaloricBinding3 {
            temperature_scale_k: 20.0,
            initial_specific_enthalpy_j_kg: 28_000.0,
            specific_enthalpy_range_j_kg: [0.0, 60_000.0],
            ..binding(false)
        },
        ImmersedCaloricBinding3 {
            fluid_cell_heat_capacity_j_k: 2.0,
            ..binding(false)
        },
    ] {
        assert_eq!(
            field.set_immersed_body_caloric(1, curve(false), mismatch, controls()),
            Err(ThermalOpenError3::Invalid(
                "caloric bodies sharing native fluid must share temperature and energy normalization",
            ))
        );
        assert_eq!([field.immersed_body(0), field.immersed_body(1)], old);
        assert_eq!(field.populations(), populations);
        assert!(field.immersed_body_caloric_curve(1).is_none());
    }
    field
        .set_immersed_body_caloric(1, curve(false), binding(false), controls())
        .unwrap();
    assert!(field.immersed_body(1).unwrap().caloric.is_some());
}
