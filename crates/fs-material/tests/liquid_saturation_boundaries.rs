//! Hydraulic boundary behavior on the native nonlinear liquid column.
//! The constitutive fixture is synthetic; references use independent m=1/2
//! closed forms and exact finite-volume inventories.

use fs_material::porous::wet::liquid::{
    LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidSupport,
    transient::{
        LiquidSaturationBoundary as Boundary, LiquidSaturationBoundaryStepInput,
        LiquidSaturationCell, LiquidSaturationColumn, LiquidSaturationImplicitControls,
        LiquidSaturationImplicitError, LiquidSaturationStepInput,
    },
};

fn cell(permeability: f64) -> LiquidSaturationCell {
    LiquidSaturationCell {
        thickness_m: 0.01,
        porosity: 0.4,
        intrinsic_permeability_m2: permeability,
        law: WetPoreLiquidLaw::try_new(
            0.1,
            0.1,
            0.5,
            0.5,
            10_000.0,
            WetPoreLiquidSupport {
                liquid_saturation: (0.1, 0.9),
                temperature_k: (280.0, 330.0),
                liquid_pressure_pa: (50_000.0, 200_000.0),
                maximum_superficial_speed_m_s: 0.1,
            },
            "synthetic hydraulic-boundary fixture; not measured material data",
        )
        .unwrap(),
    }
}

fn column(saturations: &[f64]) -> LiquidSaturationColumn {
    LiquidSaturationColumn::new(
        saturations.iter().map(|_| cell(2e-10)).collect(),
        saturations.to_vec(),
    )
    .unwrap()
}

fn input(boundaries: [Boundary; 2], duration_s: f64) -> LiquidSaturationBoundaryStepInput {
    LiquidSaturationBoundaryStepInput {
        duration_s,
        temperature_k: 300.0,
        gas_pressure_pa: 150_000.0,
        dynamic_viscosity_pa_s: 1e-3,
        body_force_pressure_gradient_pa_m: 0.0,
        boundaries,
    }
}

fn controls() -> LiquidSaturationImplicitControls {
    LiquidSaturationImplicitControls {
        absolute_saturation_residual_tolerance: 1e-11,
        relative_saturation_residual_tolerance: 0.0,
        volume_balance_tolerance_m: 1e-13,
        maximum_iterations: 128,
        maximum_line_search_trials: 48,
    }
}

fn pressure(s: f64) -> f64 {
    let effective = (s - 0.1) / 0.8;
    150_000.0 - 10_000.0 * (1.0 / (effective * effective) - 1.0).sqrt()
}

fn mobility(s: f64) -> f64 {
    let effective = (s - 0.1) / 0.8;
    let integral = 1.0 - (1.0 - effective * effective).sqrt();
    2e-10 * effective.sqrt() * integral * integral / 1e-3
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

#[test]
fn sealed_nonlinear_column_redistributes_to_inventory_determined_equilibrium() {
    let original = column(&[0.4, 0.6]);
    let saved = original.clone();
    let initial_volume = original.liquid_volume_per_area_m().unwrap();
    let mut state = original.clone();
    let mut previous_spread = 0.2;
    for _ in 0..6 {
        let (candidate, report) = state
            .step_candidate_implicit_with_boundaries(
                input([Boundary::Sealed; 2], 2.0),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(candidate.step().face_flux_m_s[0], 0.0);
        assert_eq!(candidate.step().face_flux_m_s[2], 0.0);
        assert_eq!(candidate.step().boundary_volume_per_area_m, 0.0);
        assert!(candidate.state().saturations()[0] > state.saturations()[0]);
        assert!(candidate.state().saturations()[1] < state.saturations()[1]);
        let spread = candidate.state().saturations()[1] - candidate.state().saturations()[0];
        assert!(spread < previous_spread);
        previous_spread = spread;
        assert!(report.saturation_residual_ratio <= 1.0);
        assert!(report.volume_balance_residual_m.abs() <= controls().volume_balance_tolerance_m);
        close(
            candidate.state().liquid_volume_per_area_m().unwrap(),
            initial_volume,
            controls().volume_balance_tolerance_m,
        );
        state = candidate.into_parts().0;
    }
    // Equal calibrated materials and equal storage volumes have common
    // equilibrium saturation set by the original inventory, not an ambient bath.
    for &s in state.saturations() {
        close(s, 0.5, 2e-7);
    }
    assert_eq!(state.elapsed_s(), 12.0);
    assert_eq!(original, saved);

    let explicit = original
        .step_candidate_with_boundaries(input([Boundary::Sealed; 2], 0.001), 1e-13, || false)
        .unwrap();
    close(
        explicit.state().liquid_volume_per_area_m().unwrap(),
        initial_volume,
        1e-13,
    );
    assert!(explicit.step().face_flux_m_s[1] < 0.0);
}

#[test]
fn signed_prescribed_flux_integrates_uptake_drainage_and_throughflow() {
    let original = column(&[0.5]);
    for direction in [1.0, -1.0] {
        let flux = [direction * 8e-6, -direction * 4e-6];
        let forcing = input(
            flux.map(|superficial_flux_m_s| Boundary::Flux {
                superficial_flux_m_s,
            }),
            10.0,
        );
        let (candidate, report) = original
            .step_candidate_implicit_with_boundaries(forcing, controls(), || false)
            .unwrap();
        let expected_transfer = forcing.duration_s * (flux[0] - flux[1]);
        let expected_saturation = 0.5 + expected_transfer / 0.004;
        close(
            candidate.state().saturations()[0],
            expected_saturation,
            2e-15,
        );
        assert_eq!(candidate.step().face_flux_m_s, flux.to_vec());
        close(
            candidate.step().boundary_volume_per_area_m,
            expected_transfer,
            0.0,
        );
        close(
            0.004 * (candidate.state().saturations()[0] - 0.5),
            expected_transfer,
            1e-18,
        );
        assert_eq!(report.iterations, 1);
        let explicit = original
            .step_candidate_with_boundaries(forcing, 1e-13, || false)
            .unwrap();
        close(
            explicit.state().saturations()[0],
            expected_saturation,
            2e-15,
        );
    }
    let q = Boundary::Flux {
        superficial_flux_m_s: 8e-6,
    };
    let (through, report) = original
        .step_candidate_implicit_with_boundaries(input([q; 2], 1e6), controls(), || false)
        .unwrap();
    assert_eq!(through.state().saturations(), original.saturations());
    assert_eq!(through.step().face_flux_m_s, vec![8e-6; 2]);
    assert_eq!(through.step().boundary_volume_per_area_m, 0.0);
    assert_eq!(report.iterations, 0);
}

#[test]
fn resistive_pressure_boundary_matches_nonlinear_manufactured_endpoint_and_mirror() {
    let original = column(&[0.5]);
    let target = 0.54;
    let dt = 0.2;
    let body_gradient = 150.0;
    let resistance = 3e6;
    let q = 0.004 * (target - 0.5) / dt;
    let half_cell_resistance = 0.005 / mobility(target);
    let material_pressure = pressure(target) + q * half_cell_resistance - body_gradient * 0.005;
    let external_pressure = material_pressure + q * resistance;
    let pressure_boundary = Boundary::Pressure {
        liquid_pressure_pa: external_pressure,
        surface_resistance_pa_s_m: resistance,
    };
    for left in [true, false] {
        let boundaries = if left {
            [pressure_boundary, Boundary::Sealed]
        } else {
            [Boundary::Sealed, pressure_boundary]
        };
        let mut forcing = input(boundaries, dt);
        forcing.body_force_pressure_gradient_pa_m =
            if left { body_gradient } else { -body_gradient };
        let (candidate, report) = original
            .step_candidate_implicit_with_boundaries(forcing, controls(), || false)
            .unwrap();
        close(candidate.state().saturations()[0], target, 2e-10);
        let face = if left { 0 } else { 1 };
        close(
            candidate.step().face_flux_m_s[face],
            if left { q } else { -q },
            2e-12,
        );
        assert_eq!(candidate.step().face_flux_m_s[1 - face], 0.0);
        let s = candidate.state().saturations()[0];
        let independently_evaluated = (external_pressure - pressure(s) + body_gradient * 0.005)
            / (0.005 / mobility(s) + resistance);
        close(
            candidate.step().face_flux_m_s[face].abs(),
            independently_evaluated,
            2e-16,
        );
        close(
            candidate.step().boundary_volume_per_area_m,
            0.004 * (s - 0.5),
            controls().volume_balance_tolerance_m,
        );
        assert!(report.saturation_residual_ratio <= 1.0);
    }
}

#[test]
fn zero_resistance_pressure_specialization_preserves_existing_candidates() {
    let original = column(&[0.5, 0.55]);
    let LiquidCapillaryPressure::Finite(pc) = original.states()[0].capillary_pressure() else {
        panic!("fixture must be mobile");
    };
    let legacy = LiquidSaturationStepInput {
        duration_s: 0.01,
        temperature_k: 300.0,
        gas_pressure_pa: 150_000.0,
        dynamic_viscosity_pa_s: 1e-3,
        body_force_pressure_gradient_pa_m: 150.0,
        boundary_liquid_pressure_pa: [150_000.0 - pc + 500.0, pressure(0.55)],
    };
    assert_eq!(
        original.step_candidate(legacy, || false).unwrap(),
        original
            .step_candidate_with_boundaries(legacy.into(), 1e-13, || false)
            .unwrap(),
    );
    assert_eq!(
        original
            .step_candidate_implicit(legacy, controls(), || false)
            .unwrap(),
        original
            .step_candidate_implicit_with_boundaries(legacy.into(), controls(), || false)
            .unwrap(),
    );
}

#[test]
fn material_side_pressure_support_and_exact_seals_cannot_be_bypassed_by_flux() {
    let original = column(&[0.5]);
    let saved = original.clone();
    let high_flux = input(
        [
            Boundary::Flux {
                superficial_flux_m_s: 0.05,
            },
            Boundary::Sealed,
        ],
        1e-5,
    );
    // The imposed speed is below 0.1 m/s and the cell state is admitted, but
    // its required material-side pressure exceeds the source's 200 kPa limit.
    assert!(matches!(
        original.step_candidate_implicit_with_boundaries(high_flux, controls(), || false),
        Err(LiquidSaturationImplicitError::Material(_))
    ));
    assert!(
        original
            .step_candidate_with_boundaries(high_flux, 1e-13, || false)
            .is_err()
    );

    let unsupported_pressure = Boundary::Pressure {
        liquid_pressure_pa: 300_000.0,
        surface_resistance_pa_s_m: 0.0,
    };
    assert!(
        original
            .step_candidate_implicit_with_boundaries(
                input([unsupported_pressure, Boundary::Sealed], 0.01),
                controls(),
                || false,
            )
            .is_err()
    );
    // The external bath is outside the porous source's pressure range, while
    // the actual material-side pressure behind the declared contact is inside.
    let protected = Boundary::Pressure {
        liquid_pressure_pa: 300_000.0,
        surface_resistance_pa_s_m: 1e9,
    };
    let (candidate, _) = original
        .step_candidate_implicit_with_boundaries(
            input([protected, Boundary::Sealed], 0.01),
            controls(),
            || false,
        )
        .unwrap();
    let material_face = 300_000.0 - candidate.step().face_flux_m_s[0] * 1e9;
    assert!((50_000.0..=200_000.0).contains(&material_face));
    assert!(candidate.state().saturations()[0] > 0.5);

    for (s, k) in [(0.1, 2e-10), (0.5, 0.0)] {
        let sealed = LiquidSaturationColumn::new(vec![cell(k)], vec![s]).unwrap();
        let (stationary, _) = sealed
            .step_candidate_implicit_with_boundaries(
                input([Boundary::Sealed; 2], 1e4),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(stationary.state().saturations(), sealed.saturations());
        assert_eq!(stationary.step().face_flux_m_s, vec![0.0; 2]);
        assert!(
            sealed
                .step_candidate_implicit_with_boundaries(
                    input(
                        [
                            Boundary::Flux {
                                superficial_flux_m_s: 1e-9
                            },
                            Boundary::Sealed
                        ],
                        1.0
                    ),
                    controls(),
                    || false,
                )
                .is_err()
        );
    }
    // Original coordinates are admitted; drainage changes mobility enough that
    // the required endpoint face pressure becomes unsupported. Refusal is late.
    assert!(
        original
            .step_candidate_with_boundaries(
                input(
                    [
                        Boundary::Flux {
                            superficial_flux_m_s: -0.001
                        },
                        Boundary::Sealed
                    ],
                    1.0
                ),
                1e-13,
                || false,
            )
            .is_err()
    );
    assert_eq!(original, saved);
}

#[test]
fn explicit_volume_budget_and_late_implicit_cancellation_are_transactional() {
    let original = column(&[0.5]);
    let saved = original.clone();
    let forcing = input(
        [
            Boundary::Flux {
                superficial_flux_m_s: 8e-6,
            },
            Boundary::Flux {
                superficial_flux_m_s: -4e-6,
            },
        ],
        10.0,
    );
    let admitted = original
        .step_candidate_with_boundaries(forcing, 1e-13, || false)
        .unwrap();
    assert!(admitted.step().balance_residual_m.abs() > 1e-30);
    assert!(
        original
            .step_candidate_with_boundaries(forcing, 1e-30, || false)
            .is_err()
    );
    let mut polls = 0;
    let result = original.step_candidate_implicit_with_boundaries(forcing, controls(), || {
        polls += 1;
        polls > 20
    });
    assert!(matches!(
        result,
        Err(LiquidSaturationImplicitError::Cancelled)
    ));
    assert!(polls > 20);
    assert_eq!(original, saved);
}
