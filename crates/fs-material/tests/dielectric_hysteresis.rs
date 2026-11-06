#![cfg(feature = "experimental-dielectric-hysteresis")]

use fs_material::field_storage::dielectric::hysteresis::{
    DielectricHysteresisBranch, ScalarDielectricHysteresis, UniformHystereticCapacitor,
};

fn law() -> ScalarDielectricHysteresis {
    // Scaled analytical fixture, not a fitted physical ferroelectric specimen.
    ScalarDielectricHysteresis::new(
        0.01,
        30.0,
        vec![DielectricHysteresisBranch {
            susceptibility_f_m: 0.2,
            pinning_field_v_m: 2.0,
            saturation_polarization_c_m2: 1.0,
        }],
    )
    .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}

#[test]
fn exact_switching_work_crosses_pinning_and_saturation_without_endpoint_trapezoids() {
    let law = law();
    let initial = law.initial_state();
    let pinned = law.candidate_field(&initial, 1.0).unwrap();
    near(pinned.state.displacement_c_m2(), 0.01, 1e-15);
    assert_eq!(pinned.irreversible_heat_j_m3, 0.0);
    let switched = law.candidate_field(&initial, 3.0).unwrap();
    near(switched.state.polarizations_c_m2()[0], 0.2, 1e-15);
    near(switched.state.stored_energy_j_m3(), 0.145, 1e-15);
    near(switched.irreversible_heat_j_m3, 0.4, 1e-15);
    // Independent integration: pinned dD/dE=.01 on [0,2], then .21 on [2,3].
    near(switched.electric_work_j_m3, 0.01 * 2.0 + 0.21 * 2.5, 1e-15);
    let saturated = law.candidate_field(&initial, 10.0).unwrap();
    near(saturated.state.displacement_c_m2(), 1.1, 1e-15);
    near(saturated.electric_work_j_m3, 5.0, 1e-14);
    near(saturated.incremental_permittivity_f_m, 0.01, 0.0);
    assert_eq!(initial, law.initial_state());
}

#[test]
fn major_and_minor_loops_retain_remanence_and_account_for_exact_loop_area() {
    let law = law();
    let top = law
        .candidate_field(&law.initial_state(), 10.0)
        .unwrap()
        .state;
    let zero = law.candidate_field(&top, 0.0).unwrap();
    near(zero.state.polarizations_c_m2()[0], 0.4, 1e-15);
    let mid = law.candidate_field(&zero.state, 5.0).unwrap().state;
    let low = law.candidate_field(&mid, 1.0).unwrap().state;
    let returned = law.candidate_field(&low, 5.0).unwrap().state;
    assert_eq!(returned, mid);
    let down = law.candidate_field(&top, -10.0).unwrap();
    let up = law.candidate_field(&down.state, 10.0).unwrap();
    assert_eq!(up.state, top);
    // A saturated branch switches by 2 Ps in each direction: 4 Ec Ps.
    near(
        down.irreversible_heat_j_m3 + up.irreversible_heat_j_m3,
        8.0,
        1e-14,
    );
    near(down.electric_work_j_m3 + up.electric_work_j_m3, 8.0, 1e-14);
}

#[test]
fn charge_inverse_and_tangent_use_the_original_accepted_polarization() {
    let law = law();
    let mut state = law.initial_state();
    for field in [1.0, 3.0, 10.0, 2.0, -4.0, -10.0, 0.0, 5.0] {
        let forward = law.candidate_field(&state, field).unwrap();
        let inverse = law
            .candidate_displacement(&state, forward.state.displacement_c_m2())
            .unwrap();
        near(inverse.state.field_v_m(), field, 3e-13);
        near(
            inverse.irreversible_heat_j_m3,
            forward.irreversible_heat_j_m3,
            3e-13,
        );
        let h = 1e-6;
        let plus = law.candidate_field(&state, field + h).unwrap();
        let minus = law.candidate_field(&state, field - h).unwrap();
        near(
            (plus.state.displacement_c_m2() - minus.state.displacement_c_m2()) / (2.0 * h),
            forward.incremental_permittivity_f_m,
            2e-9,
        );
        state = forward.state;
    }
}

#[test]
fn monotone_subdivision_preserves_state_work_and_heat() {
    let law = law();
    let start = law
        .candidate_field(&law.initial_state(), 10.0)
        .unwrap()
        .state;
    let whole = law.candidate_field(&start, -10.0).unwrap();
    let mut state = start;
    let (mut work, mut heat) = (0.0, 0.0);
    for k in 1..=160 {
        let step = law
            .candidate_field(&state, 10.0 - f64::from(k) / 8.0)
            .unwrap();
        work += step.electric_work_j_m3;
        heat += step.irreversible_heat_j_m3;
        state = step.state;
    }
    assert_eq!(state, whole.state);
    near(work, whole.electric_work_j_m3, 2e-13);
    near(heat, whole.irreversible_heat_j_m3, 2e-13);
}

#[test]
fn terminal_geometry_preserves_vdq_work_and_linear_capacitance_limit() {
    let capacitor = UniformHystereticCapacitor::new(law(), 0.3, 0.02).unwrap();
    let start = capacitor.initial_state();
    let step = capacitor.candidate_voltage(&start, 0.2).unwrap();
    near(step.charge_c, 0.33, 1e-15);
    near(step.electric_work_j, 0.03, 1e-15);
    near(step.switching_heat_j, 0.012, 1e-15);
    near(
        step.mean_voltage_v * step.charge_c,
        step.electric_work_j,
        1e-15,
    );
    let inverse = capacitor.candidate_charge(&start, step.charge_c).unwrap();
    near(inverse.voltage_v, step.voltage_v, 1e-14);
    let held = capacitor.hold_candidate(&step.material.state).unwrap();
    assert_eq!(held.material.state, step.material.state);
    assert_eq!(held.switching_heat_j, 0.0);
    assert_eq!(held.electric_work_j, 0.0);
    let linear = UniformHystereticCapacitor::new(
        ScalarDielectricHysteresis::new(2e-9, 1e8, vec![]).unwrap(),
        0.02,
        1e-3,
    )
    .unwrap();
    let step = linear
        .candidate_voltage(&linear.initial_state(), 100.0)
        .unwrap();
    near(step.incremental_capacitance_f, 4e-8, 1e-23);
    near(step.stored_energy_j, 0.5 * 4e-8 * 100.0 * 100.0, 1e-18);
    assert_eq!(step.switching_heat_j, 0.0);
}

#[test]
fn invalid_fields_and_foreign_history_refuse_without_mutation() {
    let law = law();
    let state = law
        .candidate_field(&law.initial_state(), 5.0)
        .unwrap()
        .state;
    let snapshot = state.clone();
    for field in [f64::NAN, f64::INFINITY, 30.1] {
        assert!(law.candidate_field(&state, field).is_err());
    }
    assert!(law.candidate_displacement(&state, 100.0).is_err());
    let other = ScalarDielectricHysteresis::new(0.02, 30.0, vec![]).unwrap();
    assert!(other.candidate_field(&state, 1.0).is_err());
    assert_eq!(state, snapshot);
    assert!(UniformHystereticCapacitor::new(law.clone(), 0.0, 1.0).is_err());
    assert!(UniformHystereticCapacitor::new(law, 1.0, f64::INFINITY).is_err());
}

mod sources {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::{
        field_storage::{
            ELECTRIC_PERMITTIVITY_DIMS,
            dielectric::{
                hysteresis::source::*,
                reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
            },
        },
        state_point::MaterialPropertySelection,
    };
    use fs_qty::{Dims, QuantitySpec, Temperature};
    use std::sync::Arc;

    fn roles() -> DielectricHysteresisSources {
        DielectricHysteresisSources {
            background_permittivity: PropertyKey::new("background", ELECTRIC_PERMITTIVITY_DIMS),
            field_limit: PropertyKey::new("field_limit", ELECTRIC_FIELD_DIMS),
            branches: vec![DielectricHysteresisBranchSources {
                susceptibility: PropertyKey::new("slope", ELECTRIC_PERMITTIVITY_DIMS),
                pinning_field: PropertyKey::new("pinning", ELECTRIC_FIELD_DIMS),
                saturation_polarization: PropertyKey::new("saturation", ELECTRIC_DISPLACEMENT_DIMS),
            }],
        }
    }
    fn point(t: f64) -> QueryPoint {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
    }
    fn card(change: impl Fn(&mut PropertyClaim)) -> Arc<MaterialCard> {
        let roles = roles();
        let mut claims = ClaimSet::new();
        for (key, value) in [
            (roles.background_permittivity, 0.01),
            (roles.field_limit, 30.0),
            (roles.branches[0].susceptibility.clone(), 0.2),
            (roles.branches[0].pinning_field.clone(), 2.0),
            (roles.branches[0].saturation_polarization.clone(), 1.0),
        ] {
            let dims = key.dims();
            let mut claim = PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    600.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored dielectric play fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            };
            change(&mut claim);
            claims.insert_claim(claim).unwrap();
        }
        Arc::new(
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic".into(),
                    phase: "solid".into(),
                    process: "explicit switching fit".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
        )
    }
    fn resolve(
        card: Arc<MaterialCard>,
        roles: DielectricHysteresisSources,
    ) -> Result<ResolvedDielectricHysteresis, fs_material::field_storage::FieldStorageError> {
        ResolvedDielectricHysteresis::resolve(
            card,
            &point(300.0),
            roles,
            MaterialPropertySelection::SingleClaimOnly,
        )
    }

    #[test]
    fn original_card_drives_switching_and_refresh_preserves_history() {
        let card = card(|_| {});
        let fit = resolve(card.clone(), roles()).unwrap();
        assert_eq!(fit.material().card_identity(), card.content_hash());
        let state = fit
            .law()
            .candidate_field(&fit.law().initial_state(), 10.0)
            .unwrap()
            .state;
        let warmer = fit.at_query_point(&point(400.0)).unwrap();
        assert_eq!(
            fit.law().candidate_field(&state, 0.0).unwrap(),
            warmer.law().candidate_field(&state, 0.0).unwrap()
        );
        assert!(fit.at_query_point(&point(601.0)).is_err());
        let capacitor = UniformHystereticCapacitor::from_source(fit, 0.3, 0.02).unwrap();
        assert_eq!(
            capacitor.source().unwrap().material().card_identity(),
            card.content_hash()
        );
        near(
            capacitor
                .candidate_voltage(&capacitor.initial_state(), 0.2)
                .unwrap()
                .switching_heat_j,
            0.012,
            1e-15,
        );
        let refreshed = capacitor.at_source_query(&point(400.0)).unwrap();
        assert_eq!(refreshed.area_m2(), capacitor.area_m2());
        assert_eq!(refreshed.gap_m(), capacitor.gap_m());
        assert_eq!(refreshed.source().unwrap().material().query(), point(400.0));
        assert_eq!(
            refreshed.candidate_voltage(&state, 0.0).unwrap(),
            capacitor.candidate_voltage(&state, 0.0).unwrap()
        );
        assert!(capacitor.at_source_query(&point(601.0)).is_err());
    }

    #[test]
    fn thermal_ownership_requires_explicit_coverage_for_every_constant_coefficient() {
        let fit = resolve(card(|_| {}), roles()).unwrap();
        fit.admit_temperature_interval("T", [250.0, 600.0]).unwrap();
        for interval in [[249.0, 600.0], [250.0, 601.0], [300.0, 300.0], [0.0, 600.0]] {
            assert!(fit.admit_temperature_interval("T", interval).is_err());
        }
        assert!(
            fit.admit_temperature_interval("frequency", [250.0, 600.0])
                .is_err()
        );
        let missing = resolve(
            card(|claim| {
                if claim.key.name() == "pinning" {
                    claim.validity = ValidityDomain::unconstrained();
                }
            }),
            roles(),
        )
        .unwrap();
        assert!(
            missing
                .admit_temperature_interval("T", [290.0, 310.0])
                .is_err()
        );
    }

    #[test]
    fn relative_permittivity_dynamic_support_and_changing_fits_are_not_frozen() {
        let mut wrong = roles();
        wrong.background_permittivity = PropertyKey::new("background", Dims::NONE);
        assert!(resolve(card(|_| {}), wrong).is_err());
        assert!(
            resolve(
                card(|claim| {
                    claim.validity = claim.validity.clone().with_quantity(
                        "frequency",
                        QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                        0.0,
                        100.0,
                    );
                }),
                roles()
            )
            .is_err()
        );
        assert!(
            resolve(
                card(|claim| {
                    if claim.key.name() == "pinning" {
                        claim.value = PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: vec![(250.0, 2.0), (600.0, 1.0)],
                            dims: ELECTRIC_FIELD_DIMS,
                        };
                        claim.interpolation = InterpolationPolicy::LinearInside;
                    }
                }),
                roles()
            )
            .is_err()
        );
    }
}
