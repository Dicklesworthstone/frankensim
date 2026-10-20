use fs_cheb::Cheb1;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, ChebProfile, Combination, Compat, Environment, Interp,
    LoadCase, Physics, Scenario, TimeSignal, ValidationBudget, ValidationError, parse_ir, write_ir,
};

const MASS: Dims = Dims([0, 1, -1, 0, 0, 0]);

fn inlet(name: &str, signal: TimeSignal) -> BoundaryCondition {
    BoundaryCondition {
        region: name.into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Signal(signal)),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    }
}
fn constant(name: &str, value: f64) -> BoundaryCondition {
    inlet(name, TimeSignal::Constant(QtyAny::new(value, MASS)))
}
fn ramp(from: f64, to: f64) -> TimeSignal {
    TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 2.0,
        from: QtyAny::new(from, MASS),
        to: QtyAny::new(to, MASS),
    }
}
fn table(times: Vec<f64>, values: Vec<f64>, interp: Interp) -> TimeSignal {
    TimeSignal::Table {
        times,
        values,
        dims: MASS,
        interp,
    }
}
fn scenario(bcs: Vec<BoundaryCondition>) -> Scenario {
    let mut scenario = Scenario::new("all-time-flux", 7, Environment::earth_lab());
    scenario.base_bcs = bcs;
    scenario
}
fn assert_admitted(s: &Scenario) {
    assert!(s.validate().is_empty(), "{:?}", s.validate());
    let wire = write_ir(s);
    let decoded = parse_ir(&wire).unwrap();
    assert!(decoded.scenario().validate().is_empty());
    assert_eq!(write_ir(decoded.scenario()), wire);
}

#[test]
fn opposite_transient_histories_admit_through_zero_crossings_and_clamped_tails() {
    let families = [
        (ramp(-2.0, 3.0), ramp(2.0, -3.0)),
        (
            table(vec![-1.0, 0.0, 2.0], vec![-2.0, 1.0, 3.0], Interp::Linear),
            table(vec![-1.0, 0.0, 2.0], vec![2.0, -1.0, -3.0], Interp::Linear),
        ),
        (
            table(vec![-1.0, 0.0, 2.0], vec![-2.0, 1.0, 3.0], Interp::Hold),
            table(vec![-1.0, 0.0, 2.0], vec![2.0, -1.0, -3.0], Interp::Hold),
        ),
        (
            TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 2.0, vec![0.5, 2.0, -0.75, 0.25]),
                dims: MASS,
            }),
            TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 2.0, vec![-0.5, -2.0, 0.75, -0.25]),
                dims: MASS,
            }),
        ),
    ];
    for (a, b) in families {
        let mut s = scenario(vec![
            inlet("in", a),
            constant("offset-in", 3.0),
            inlet("out", b),
            constant("offset-out", -3.0),
        ]);
        assert_admitted(&s);
        for time in [-f64::MAX, -1.0, -0.0, 0.0, 0.3, 1.0, 2.0, f64::MAX] {
            let fluxes: Vec<_> = s
                .base_bcs
                .iter()
                .map(|bc| bc.mass_flow_at(time).unwrap().unwrap())
                .collect();
            let gross: f64 = fluxes.iter().map(|v| v.abs()).sum();
            assert!(fluxes.iter().sum::<f64>().abs() <= 1e-9 * gross);
        }
        s.base_bcs.reverse();
        assert_admitted(&s);
    }
}

#[test]
fn different_linear_and_hold_grids_are_proved_on_the_union_partition() {
    for scale in [1.0, 1e-200, 1e100] {
        assert_admitted(&scenario(vec![
            inlet("in", ramp(scale, 3.0 * scale)),
            inlet(
                "out",
                table(
                    vec![0.0, 0.5, 1.0, 2.0],
                    vec![-scale, -1.5 * scale, -2.0 * scale, -3.0 * scale],
                    Interp::Linear,
                ),
            ),
        ]));
    }
    assert_admitted(&scenario(vec![
        inlet(
            "in",
            table(vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 3.0], Interp::Hold),
        ),
        inlet(
            "out",
            table(
                vec![0.0, 0.5, 1.0, 1.5, 2.0],
                vec![-1.0, -1.0, -2.0, -2.0, -3.0],
                Interp::Hold,
            ),
        ),
    ]));
}

#[test]
fn held_time_zero_has_one_physical_value_for_both_zero_encodings() {
    let a = table(vec![-1.0, -0.0, 1.0], vec![1.0, 2.0, 3.0], Interp::Hold);
    let b = table(vec![-1.0, 0.0, 1.0], vec![-1.0, -2.0, -3.0], Interp::Hold);
    for time in [-0.0, 0.0] {
        assert_eq!(a.eval(time).unwrap().value, 2.0);
        assert_eq!(b.eval(time).unwrap().value, -2.0);
    }
    assert_admitted(&scenario(vec![inlet("in", a), inlet("out", b)]));
}

#[test]
fn smooth_constant_reduction_is_a_proof_not_a_sampled_guess() {
    assert_admitted(&scenario(vec![
        inlet(
            "smooth",
            TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(1.0, 9.0, vec![2.0]),
                dims: MASS,
            }),
        ),
        constant("out", -1.0),
    ]));
    assert_admitted(&scenario(vec![inlet("zero", ramp(0.0, 0.0))]));
}

#[test]
fn opposite_smooth_histories_keep_extreme_finite_domains() {
    for (start, end) in [
        (0.0, f64::MAX),
        (-f64::MAX, f64::MAX),
        (-f64::from_bits(1), f64::from_bits(1)),
    ] {
        let smooth = |sign| {
            TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(start, end, vec![sign, 2.0 * sign, -sign]),
                dims: MASS,
            })
        };
        assert_admitted(&scenario(vec![
            inlet("in", smooth(1.0)),
            inlet("out", smooth(-1.0)),
        ]));
    }
}

#[test]
fn opposite_spectral_pair_composes_with_a_different_grid_affine_balance() {
    let smooth = |sign| {
        TimeSignal::Chebfun(ChebProfile {
            cheb: Cheb1::from_coeffs(0.0, 2.0, vec![sign, 2.0 * sign, -sign]),
            dims: MASS,
        })
    };
    assert_admitted(&scenario(vec![
        inlet("ramp", ramp(1.0, 3.0)),
        inlet("spectral-in", smooth(1.0)),
        inlet(
            "table",
            table(vec![0.0, 1.0, 2.0], vec![-1.0, -2.0, -3.0], Interp::Linear),
        ),
        inlet("spectral-out", smooth(-1.0)),
    ]));
}

#[test]
fn unmatched_domains_and_interior_relative_imbalance_stay_refused() {
    let a = TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(0.0, 1.0, vec![0.0, 1.0]),
        dims: MASS,
    });
    let b = TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(0.0, 2.0, vec![0.0, -1.0]),
        dims: MASS,
    });
    assert!(
        !scenario(vec![inlet("in", a), inlet("out", b)])
            .validate()
            .is_empty()
    );
    let s = scenario(vec![
        inlet("a", ramp(1.0, -1.0)),
        inlet("b", ramp(-1.0 + 1e-9, 1.0 + 1e-9)),
    ]);
    assert!(
        s.validate()
            .iter()
            .any(|v| v.code == "flux-certification-unavailable")
    );
    let s = scenario(vec![
        inlet("a", ramp(0.0, 1.0)),
        inlet("b", table(vec![0.0, 2.0], vec![0.0, -1.0], Interp::Hold)),
    ]);
    assert!(
        s.validate()
            .iter()
            .any(|v| v.code == "flux-certification-unavailable")
    );
}

#[test]
fn rounded_coefficient_cancellation_cannot_erase_an_interior_flux() {
    let scale = 2.0f64.powi(100);
    let smooth = |sign| {
        TimeSignal::Chebfun(ChebProfile {
            cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![-sign * (5.0 / 256.0) * scale, sign * scale]),
            dims: MASS,
        })
    };
    let s = scenario(vec![
        inlet("large-in", smooth(1.0)),
        constant("small-unbalanced", 1.0),
        inlet("large-out", smooth(-1.0)),
    ]);
    // The large opposite histories hide +1 at every fixed smooth checkpoint.
    // At their shared root, between checkpoints, it is the entire gross flux.
    let flux: f64 = s
        .base_bcs
        .iter()
        .map(|bc| bc.mass_flow_at(5.0 / 512.0).unwrap().unwrap())
        .sum();
    assert_eq!(flux, 1.0);
    let findings = s.validate();
    assert!(
        findings
            .iter()
            .any(|v| v.code == "flux-certification-unavailable"),
        "{findings:?}"
    );
    assert!(findings.iter().all(|v| v.code != "flux-imbalance"));
}

#[test]
fn scaled_combinations_are_checked_and_zero_outlets_do_not_absorb_flux() {
    let mut s = scenario(vec![
        constant("base-in", 1.0),
        constant("base-out", -1.0 + 1e-9),
    ]);
    s.cases.push(LoadCase {
        name: "trim".into(),
        bcs: vec![constant("trim", -1e-9)],
    });
    s.cases.push(LoadCase {
        name: "relief".into(),
        bcs: vec![BoundaryCondition {
            region: "relief".into(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::PressureOutlet,
            value: Some(BcValue::Uniform(QtyAny::new(
                0.0,
                Dims([-1, 1, -2, 0, 0, 0]),
            ))),
            compatibility: None,
            frame: 0,
        }],
    });
    assert_admitted(&s);
    s.combinations.push(Combination {
        name: "amplified".into(),
        terms: vec![("trim".into(), 100.0), ("relief".into(), 0.0)],
    });
    let findings = s.validate();
    assert!(
        findings
            .iter()
            .any(|v| v.code == "flux-imbalance" && v.what.contains("combination row 0")),
        "{findings:?}"
    );
    s.combinations[0].terms[1].1 = 1.0;
    assert_admitted(&s);
}

#[test]
fn combination_scaling_checks_overflow_and_signed_transient_terms() {
    let mut s = scenario(vec![]);
    s.cases.push(LoadCase {
        name: "balanced".into(),
        bcs: vec![inlet("in", ramp(-1.0, 2.0)), inlet("out", ramp(1.0, -2.0))],
    });
    s.combinations.push(Combination {
        name: "reverse".into(),
        terms: vec![("balanced".into(), -3.0)],
    });
    assert_admitted(&s);
    s.cases[0].bcs = vec![
        constant("in", f64::MAX / 4.0),
        constant("out", -f64::MAX / 4.0),
    ];
    s.combinations[0].terms[0].1 = 8.0;
    assert!(
        s.validate()
            .iter()
            .any(|v| v.code == "flux-aggregation-nonfinite" && v.what.contains("combination"))
    );
}

#[test]
fn all_time_proof_and_expanded_combinations_are_charged_before_admission() {
    let mut s = scenario(vec![]);
    s.cases.push(LoadCase {
        name: "balanced".into(),
        bcs: vec![inlet("in", ramp(1.0, 2.0)), inlet("out", ramp(-1.0, -2.0))],
    });
    let base = s.validation_plan(ValidationBudget::default()).unwrap();
    s.combinations.push(Combination {
        name: "twice".into(),
        terms: vec![("balanced".into(), 2.0)],
    });
    let plan = s.validation_plan(ValidationBudget::default()).unwrap();
    assert!(plan.flux_checkpoints > base.flux_checkpoints);
    assert!(plan.planned_work > base.planned_work);
    let mut budget = ValidationBudget::default();
    budget.max_work = plan.planned_work - 1;
    assert!(matches!(
        s.validation_plan(budget),
        Err(ValidationError::WorkExceeded { .. })
    ));
}
