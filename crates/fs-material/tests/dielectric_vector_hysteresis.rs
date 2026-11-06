//! Electric vector history and frozen-history variational quantities.
//! Coefficients below are scaled analytic fixtures, not measured ferroelectrics.
#![cfg(feature = "experimental-dielectric-hysteresis")]

use fs_material::field_storage::dielectric::hysteresis::{
    DielectricHysteresisBranch, ScalarDielectricHysteresis,
    vector::{
        ResolvedVectorDielectricHysteresis, VectorDielectricHysteresis,
        VectorDielectricHysteresisControls, VectorDielectricHysteresisError,
        VectorDielectricHysteresisState, VectorDielectricHysteresisStep,
    },
};

type Vector = [f64; 3];

fn branch(chi: f64, pinning: f64, saturation: f64) -> DielectricHysteresisBranch {
    DielectricHysteresisBranch {
        susceptibility_f_m: chi,
        pinning_field_v_m: pinning,
        saturation_polarization_c_m2: saturation,
    }
}
fn law() -> VectorDielectricHysteresis {
    VectorDielectricHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 1.0)]).unwrap()
}
fn controls() -> VectorDielectricHysteresisControls {
    VectorDielectricHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_v_m: 1e-12,
        absolute_displacement_tolerance_c_m2: 1e-12,
        work_tolerance_j_m3: 2e-11,
        maximum_branch_iterations: 100,
        maximum_inverse_iterations: 60,
        maximum_backtracks: 40,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:e}"
    );
}
fn vector_near(actual: Vector, expected: Vector, tolerance: f64) {
    for i in 0..3 {
        near(actual[i], expected[i], tolerance);
    }
}
fn dot(a: Vector, b: Vector) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: Vector) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn sub(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[i] - b[i])
}
fn step(
    law: &VectorDielectricHysteresis,
    previous: &VectorDielectricHysteresisState,
    field: Vector,
) -> VectorDielectricHysteresisStep {
    law.candidate_field(previous, field, controls(), &mut || false)
        .unwrap()
}

#[test]
fn collinear_electric_history_matches_scalar_memory_and_physical_heat() {
    let scalar = ScalarDielectricHysteresis::new(
        0.01,
        30.0,
        vec![
            branch(0.2, 2.0, 1.0),
            branch(0.07, 4.0, 0.4),
            branch(0.1, 0.0, 0.6),
        ],
    )
    .unwrap();
    let vector = VectorDielectricHysteresis::from_isotropic_scalar(&scalar).unwrap();
    assert_eq!(vector.background_permittivity_f_m(), 0.01);
    assert_eq!(vector.field_limit_v_m(), 30.0);
    let mut a = scalar.initial_state();
    let mut b = vector.initial_state();
    for field in [
        0.0, 1.0, 3.0, 12.0, 2.0, -3.0, -12.0, 4.0, 0.0, 8.0, 1.0, 8.0,
    ] {
        let x = scalar.candidate_field(&a, field).unwrap();
        let y = step(&vector, &b, [field, 0.0, 0.0]);
        vector_near(
            y.state.displacement_c_m2(),
            [x.state.displacement_c_m2(), 0.0, 0.0],
            2e-12,
        );
        near(
            y.state.stored_energy_j_m3(),
            x.state.stored_energy_j_m3(),
            2e-12,
        );
        near(y.irreversible_heat_j_m3, x.irreversible_heat_j_m3, 2e-12);
        for (p, scalar_p) in y
            .state
            .polarizations_c_m2()
            .iter()
            .zip(x.state.polarizations_c_m2())
        {
            vector_near(*p, [*scalar_p, 0.0, 0.0], 2e-12);
        }
        near(
            y.energy_change_j_m3 + y.irreversible_heat_j_m3,
            x.electric_work_j_m3,
            3e-12,
        );
        assert!(y.numerical_dissipation_j_m3 >= 0.0);
        near(
            y.endpoint_work_j_m3,
            x.electric_work_j_m3 + y.numerical_dissipation_j_m3,
            2e-11,
        );
        a = x.state;
        b = y.state;
    }
}

#[test]
fn rotation_changes_vector_coordinates_without_changing_memory_losses() {
    let material = law();
    // A proper orthogonal rotation mixes all axes; componentwise scalar play
    // would pass a mere axis permutation but fails this frame transformation.
    let r = [
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    ];
    let rotate = |v: Vector| r.map(|row| dot(row, v));
    let mut a = material.initial_state();
    let mut b = material.initial_state();
    for field in [
        [3.0, 0.0, 0.0],
        [3.0, 4.0, 1.0],
        [2.0, 9.0, -1.0],
        [-1.0, 4.0, 2.0],
        [-7.0, -5.0, -2.0],
        [0.0; 3],
        [4.0, 6.0, 3.0],
    ] {
        let original = step(&material, &a, field);
        let rotated = step(&material, &b, rotate(field));
        vector_near(
            rotated.state.displacement_c_m2(),
            rotate(original.state.displacement_c_m2()),
            2e-10,
        );
        for (&p, &q) in original
            .state
            .polarizations_c_m2()
            .iter()
            .zip(rotated.state.polarizations_c_m2())
        {
            vector_near(q, rotate(p), 2e-10);
        }
        near(
            rotated.state.stored_energy_j_m3(),
            original.state.stored_energy_j_m3(),
            2e-10,
        );
        near(
            rotated.irreversible_heat_j_m3,
            original.irreversible_heat_j_m3,
            2e-10,
        );
        near(
            rotated.numerical_dissipation_j_m3,
            original.numerical_dissipation_j_m3,
            2e-10,
        );
        for i in 0..3 {
            for j in 0..3 {
                let mut expected = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        expected += r[i][k] * original.incremental_permittivity_f_m[k][l] * r[j][l];
                    }
                }
                near(rotated.incremental_permittivity_f_m[i][j], expected, 2e-10);
            }
        }
        a = original.state;
        b = rotated.state;
    }
}

#[test]
fn saturated_rotating_response_matches_independent_angular_stationarity_and_inverse() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    let snapshot = old.clone();
    let field = [3.0, 8.0, 0.0];
    let response = step(&material, &old, field);
    // Minimize 2.5 - 3 cos(theta) - 8 sin(theta)
    //          + 2 sqrt(1.36 - 1.2 cos(theta)) on the saturation circle.
    // This independent angular derivative uses neither the native shrinkage
    // formula nor its Lagrange-multiplier equation.
    let (mut lo, mut hi) = (0.0, core::f64::consts::FRAC_PI_2);
    for _ in 0..80 {
        let angle = f64::midpoint(lo, hi);
        let derivative = 3.0 * angle.sin() - 8.0 * angle.cos()
            + 1.2 * angle.sin() / (1.36 - 1.2 * angle.cos()).sqrt();
        if derivative < 0.0 {
            lo = angle;
        } else {
            hi = angle;
        }
    }
    let angle = f64::midpoint(lo, hi);
    let expected = [angle.cos(), angle.sin(), 0.0];
    vector_near(response.state.polarizations_c_m2()[0], expected, 2e-11);
    near(
        response.irreversible_heat_j_m3,
        2.0 * norm(sub(expected, [0.6, 0.0, 0.0])),
        4e-11,
    );
    let inverse = material
        .candidate_displacement(
            &old,
            response.state.displacement_c_m2(),
            controls(),
            &mut || false,
        )
        .unwrap();
    vector_near(inverse.state.field_v_m(), field, 3e-9);
    vector_near(
        inverse.state.displacement_c_m2(),
        response.state.displacement_c_m2(),
        3e-11,
    );
    near(
        inverse.irreversible_heat_j_m3,
        response.irreversible_heat_j_m3,
        3e-9,
    );
    assert_eq!(old, snapshot);
}

#[test]
fn frozen_electric_potentials_and_tangents_are_work_conjugate() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    let snapshot = old.clone();
    let initial_coenergy = dot(old.field_v_m(), old.displacement_c_m2()) - old.stored_energy_j_m3();
    // Pinned, switching, and saturated states, each away from a kink.
    for field in [[3.0, 0.0, 0.0], [3.0, 3.0, 0.5], [3.0, 8.0, 0.5]] {
        let response = step(&material, &old, field);
        let p = response.state.polarizations_c_m2()[0];
        let direct_coenergy = 0.005 * dot(field, field) + dot(field, p)
            - 2.5 * dot(p, p)
            - 2.0 * norm(sub(p, old.polarizations_c_m2()[0]));
        near(
            response.incremental_coenergy_j_m3().unwrap(),
            direct_coenergy,
            3e-13,
        );
        near(
            response.incremental_coenergy_change_j_m3() + initial_coenergy,
            direct_coenergy,
            3e-13,
        );
        let inverse = response.incremental_inverse_permittivity_m_f().unwrap();
        for j in 0..3 {
            for i in 0..3 {
                let product: f64 = (0..3)
                    .map(|k| response.incremental_permittivity_f_m[i][k] * inverse[k][j])
                    .sum();
                near(product, if i == j { 1.0 } else { 0.0 }, 3e-13);
                assert_eq!(inverse[i][j], inverse[j][i]);
            }
            let delta = 1e-5;
            let mut plus = field;
            let mut minus = field;
            plus[j] += delta;
            minus[j] -= delta;
            let high = step(&material, &old, plus);
            let low = step(&material, &old, minus);
            let coenergy_change = material
                .incremental_coenergy_difference_j_m3(&old, &low.state, &high.state)
                .unwrap();
            near(
                coenergy_change / (2.0 * delta),
                response.state.displacement_c_m2()[j],
                2e-8,
            );
            for i in 0..3 {
                near(
                    (high.state.displacement_c_m2()[i] - low.state.displacement_c_m2()[i])
                        / (2.0 * delta),
                    response.incremental_permittivity_f_m[i][j],
                    2e-7,
                );
            }
            let mut dp = response.state.displacement_c_m2();
            let mut dm = dp;
            dp[j] += delta;
            dm[j] -= delta;
            let high = material
                .candidate_displacement(&old, dp, controls(), &mut || false)
                .unwrap();
            let low = material
                .candidate_displacement(&old, dm, controls(), &mut || false)
                .unwrap();
            near(
                (high.incremental_potential_j_m3().unwrap()
                    - low.incremental_potential_j_m3().unwrap())
                    / (2.0 * delta),
                field[j],
                3e-7,
            );
        }
    }
    assert_eq!(old, snapshot);
}

#[test]
fn pairwise_coenergy_keeps_small_trial_changes_beneath_unrelated_field_energy() {
    let material = VectorDielectricHysteresis::new(1.0, 2e8, vec![]).unwrap();
    let old = material.initial_state();
    let from = step(&material, &old, [1e8, 0.1, 0.0]);
    let to = step(&material, &old, [1e8, 0.2, 0.0]);
    assert_eq!(
        from.incremental_coenergy_j_m3().unwrap(),
        to.incremental_coenergy_j_m3().unwrap()
    );
    assert_eq!(
        from.incremental_coenergy_change_j_m3(),
        to.incremental_coenergy_change_j_m3()
    );
    near(
        material
            .incremental_coenergy_difference_j_m3(&old, &from.state, &to.state)
            .unwrap(),
        0.015,
        1e-17,
    );
    near(
        material
            .incremental_coenergy_difference_j_m3(&old, &to.state, &from.state)
            .unwrap(),
        -0.015,
        1e-17,
    );
    assert_eq!(from.irreversible_heat_j_m3, 0.0);
    assert_eq!(to.irreversible_heat_j_m3, 0.0);
    assert!(to.numerical_dissipation_j_m3 > 0.0);
}

#[test]
fn rotating_loss_converges_to_an_independent_continuous_loop_area() {
    let material =
        VectorDielectricHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 10.0)]).unwrap();
    // Steady continuous polarization lags E by asin(Ec/Eamp):
    // |P|=chi sqrt(Eamp²-Ec²), Q/cycle=2 pi Ec |P|.
    let exact = 2.0 * core::f64::consts::PI * 2.0 * 0.2 * 21.0_f64.sqrt();
    let mut errors = Vec::new();
    let mut numerical = Vec::new();
    for count in [64, 128, 256] {
        let mut state = material.initial_state();
        let (mut heat, mut remainder, mut work) = (0.0, 0.0, 0.0);
        let mut before = 0.0;
        for i in 0..=8 * count {
            let angle = 2.0 * core::f64::consts::PI * i as f64 / count as f64;
            let trial = step(
                &material,
                &state,
                [5.0 * angle.cos(), 5.0 * angle.sin(), 0.0],
            );
            if i == 7 * count {
                before = trial.state.stored_energy_j_m3();
            }
            if i > 7 * count {
                heat += trial.irreversible_heat_j_m3;
                remainder += trial.numerical_dissipation_j_m3;
                work += trial.endpoint_work_j_m3;
            }
            state = trial.state;
        }
        assert!(state.polarizations_c_m2()[0][1] < -0.1);
        near(
            work,
            heat + remainder + state.stored_energy_j_m3() - before,
            2e-10,
        );
        errors.push((heat - exact).abs());
        numerical.push(remainder);
    }
    for pair in errors.windows(2) {
        assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    for pair in numerical.windows(2) {
        assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{numerical:?}");
    }
}

#[test]
fn vector_domains_foreign_memory_budgets_and_cancellation_are_atomic() {
    let material = law();
    let old = step(&material, &material.initial_state(), [5.0, 0.0, 0.0]).state;
    let snapshot = old.clone();
    for field in [
        [f64::NAN, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [25.0, 25.0, 0.0],
    ] {
        let error = material
            .candidate_field(&old, field, controls(), &mut || false)
            .unwrap_err();
        let text = error.to_string();
        assert!(
            !text.contains("magnetic") && !text.contains("permeability") && !text.contains("flux")
        );
    }
    let other = VectorDielectricHysteresis::new(0.02, 30.0, vec![]).unwrap();
    assert!(
        other
            .candidate_field(&old, [0.0; 3], controls(), &mut || false)
            .is_err()
    );
    assert!(
        material
            .incremental_coenergy_difference_j_m3(&old, &old, &other.initial_state())
            .is_err()
    );
    let mut bounded = controls();
    bounded.maximum_branch_iterations = 0;
    assert!(matches!(
        material.candidate_field(&old, [3.0, 8.0, 0.0], bounded, &mut || false),
        Err(VectorDielectricHysteresisError::BranchNonconvergence { .. })
    ));
    bounded = controls();
    bounded.maximum_inverse_iterations = 0;
    assert!(matches!(
        material.candidate_displacement(&old, [0.0; 3], bounded, &mut || false),
        Err(VectorDielectricHysteresisError::DisplacementNonconvergence { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        material.candidate_field(&old, [3.0, 8.0, 0.0], controls(), &mut || {
            polls += 1;
            polls >= 3
        }),
        Err(VectorDielectricHysteresisError::Cancelled)
    ));
    assert!(polls >= 3);
    assert_eq!(old, snapshot);
    assert_eq!(
        step(&material, &old, [3.0, 8.0, 0.0]),
        step(&material, &old, [3.0, 8.0, 0.0])
    );
}

#[test]
fn vector_source_refresh_preserves_original_fit_query_and_prepared_history() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::{
        field_storage::{
            ELECTRIC_PERMITTIVITY_DIMS,
            dielectric::{
                hysteresis::source::{
                    DielectricHysteresisBranchSources, DielectricHysteresisSources,
                    ResolvedDielectricHysteresis,
                },
                reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
            },
        },
        state_point::MaterialPropertySelection,
    };
    use fs_qty::{QuantitySpec, Temperature};
    use std::sync::Arc;

    let roles = DielectricHysteresisSources {
        background_permittivity: PropertyKey::new("epsilon_background", ELECTRIC_PERMITTIVITY_DIMS),
        field_limit: PropertyKey::new("field_support", ELECTRIC_FIELD_DIMS),
        branches: vec![DielectricHysteresisBranchSources {
            susceptibility: PropertyKey::new("susceptibility", ELECTRIC_PERMITTIVITY_DIMS),
            pinning_field: PropertyKey::new("pinning", ELECTRIC_FIELD_DIMS),
            saturation_polarization: PropertyKey::new("saturation", ELECTRIC_DISPLACEMENT_DIMS),
        }],
    };
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (roles.background_permittivity.clone(), 0.01),
        (roles.field_limit.clone(), 30.0),
        (roles.branches[0].susceptibility.clone(), 0.2),
        (roles.branches[0].pinning_field.clone(), 2.0),
        (roles.branches[0].saturation_polarization.clone(), 1.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    600.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored uniaxial electric fit; isotropy explicitly selected".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "fixed isotropic switching fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let point = |temperature| {
        QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                temperature,
            )
            .unwrap()
    };
    let scalar = ResolvedDielectricHysteresis::resolve(
        card.clone(),
        &point(300.0),
        roles.clone(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = ResolvedVectorDielectricHysteresis::from_isotropic_source(scalar).unwrap();
    assert_eq!(source.material().card_identity(), card.content_hash());
    assert_eq!(source.sources(), &roles);
    let old = step(source.law(), &source.law().initial_state(), [3.0, 8.0, 0.0]).state;
    let warmer = source.at_query_point(&point(400.0)).unwrap();
    assert_eq!(warmer.material().query(), point(400.0));
    assert_eq!(
        step(source.law(), &old, [0.0; 3]),
        step(warmer.law(), &old, [0.0; 3])
    );
    assert!(source.at_query_point(&point(601.0)).is_err());
}
