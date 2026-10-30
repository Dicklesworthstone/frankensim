//! G1/G3/G4: actual full-matrix slip closes original speed/path queries.
#[path = "trajectory.rs"]
mod trajectory;

#[path = "driven.rs"]
mod driven;

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::motion::{
    MotionQueryError, MotionQuerySolve,
};
use fs_contact::interface_binding::thermal_friction::{FrictionQueryAxis, FrictionQueryAxisKind};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn thermal_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point() -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        ("T", 300.0),
        ("Ta", 290.0),
        ("Tb", 310.0),
        ("ambient", 295.0),
    ] {
        p = p.with_quantity(name, thermal_quantity(), value).unwrap();
    }
    p.with_quantity("p", pressure_quantity(), 1e5)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.5)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap()
}

fn endpoint() -> QueryPoint {
    // Motion entries are deliberately stale placeholders. The source checkpoint
    // supplies the accepted origin; the closure owns their new physical values.
    point()
        .with("T", 320.0)
        .unwrap()
        .with("Ta", 310.0)
        .unwrap()
        .with("Tb", 330.0)
        .unwrap()
        .with("p", 2e5)
        .unwrap()
        .with("speed", 99.0)
        .unwrap()
        .with("path", 99.0)
        .unwrap()
}

fn motion_curve(name: &str, axis: &str, dims: Dims, values: [f64; 2]) -> PropertyClaim {
    let mut claim = curve(name, values);
    claim.value = PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: dims,
        knots: vec![(0.0, values[0]), (2.0, values[1])],
        dims: Dims::NONE,
    };
    claim.validity = ValidityDomain::unconstrained().with_quantity(
        axis,
        QuantitySpec::dimensional(dims),
        0.0,
        2.0,
    );
    claim
}

fn motion_card(x: [f64; 2], y: [f64; 2]) -> InterfaceSystemCard {
    let constant = |name| {
        let mut claim = curve(name, [1.2; 2]);
        claim.validity =
            ValidityDomain::unconstrained().with_quantity("T", thermal_quantity(), 280.0, 320.0);
        claim
    };
    card(vec![
        constant(STATIC),
        constant(STATIC_Y),
        motion_curve(KINETIC, "speed", SPEED, x),
        motion_curve(KINETIC_Y, "path", LENGTH, y),
    ])
}

fn motion_law(card: &InterfaceSystemCard, axis: [f64; 2]) -> SourceThermalEllipticCoulomb {
    let principal = [[STATIC, KINETIC], [STATIC_Y, KINETIC_Y]].map(|names| {
        let owned = [
            ("Ta", FrictionQueryAxisKind::Temperature),
            ("Tb", FrictionQueryAxisKind::Temperature),
            ("speed", FrictionQueryAxisKind::SlidingSpeed),
            ("path", FrictionQueryAxisKind::AccumulatedDistance),
        ]
        .into_iter()
        .map(|(name, kind)| FrictionQueryAxis {
            name: name.into(),
            kind,
        })
        .collect();
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
            point(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: Some("p".into()),
            },
            owned,
        )
        .unwrap()
    });
    SourceThermalEllipticCoulomb::from_principal_sources(principal, axis).unwrap()
}

fn solve(bracket: [f64; 2]) -> MotionQuerySolve {
    MotionQuerySolve {
        distance_bracket_m: bracket,
        absolute_distance_tolerance_m: 2e-13,
        relative_distance_tolerance: 1e-13,
        maximum_evaluations: 96,
    }
}

fn drive(free: [f64; 2]) -> TangentialStepInput {
    TangentialStepInput {
        dt_s: 0.5,
        nominal_pressure_pa: Some(2e5),
        ..input(320.0, free)
    }
}

fn norm(v: [f64; 2]) -> f64 {
    fs_math::det::hypot(v[0], v[1])
}

#[test]
fn independent_full_matrix_roots_refresh_actual_speed_and_original_path_after_reversal() {
    let card = motion_card([0.3, 0.5], [0.2, 0.4]);
    let axis = [0.6, 0.8];
    let law = motion_law(&card, axis);
    let initial = law.initial_state();
    let mut state = initial.clone();
    let matrix = [[0.04, -0.015], [-0.015, 0.02]];
    let compliance = SymmetricTangentialCompliance::new(matrix).unwrap();
    let world = |p: [f64; 2]| {
        [
            axis[0] * p[0] - axis[1] * p[1],
            axis[1] * p[0] + axis[0] * p[1],
        ]
    };
    for (length, path, sign) in [(0.5, 0.75, 1.0), (0.25, 1.0, -1.0)] {
        let speed = length / 0.5;
        let caps = [10.0 * (0.3 + 0.1 * speed), 10.0 * (0.2 + 0.1 * path)];
        let direction = [0.6 / caps[0], 0.8 / caps[1]];
        let multiplier = length / norm(direction);
        let slip = world(direction.map(|x| sign * multiplier * x));
        let force = world([sign * 0.6 * caps[0], sign * 0.8 * caps[1]]);
        let free = [
            slip[0] + matrix[0][0] * force[0] + matrix[0][1] * force[1],
            slip[1] + matrix[1][0] * force[0] + matrix[1][1] * force[1],
        ];
        let report = law
            .advance_motion_consistent(
                &state,
                &drive(free),
                &endpoint(),
                EllipticFrictionMode::Kinetic,
                compliance,
                solve([0.0, 1.0]),
                || false,
            )
            .unwrap();
        assert!(report.distance_residual_m <= 0.0);
        assert!(report.distance_residual_m.abs() <= 3e-13);
        close(report.actual_distance_m, length);
        close(report.queried_distance_m, length);
        close(report.step.dissipated_j, multiplier);
        assert!(report.evaluations > 2 && report.evaluations <= 96);
        for i in 0..2 {
            close(report.step.traction_n[i], force[i]);
            close(report.step.slip_increment_m[i], slip[i]);
            let source = &report.step.state.friction_sources()[i];
            let query = source.query();
            close(query.axes()["speed"], speed);
            close(query.axes()["path"], path);
            for (name, value) in [
                ("T", 320.0),
                ("Ta", 310.0),
                ("Tb", 330.0),
                ("p", 2e5),
                ("ambient", 295.0),
            ] {
                assert_eq!(query.axes()[name], value);
            }
            for property in source.properties() {
                let original = initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        assert_eq!(
            report.step.state.constitutive_state().principal_axis(),
            axis
        );
        state = report.step.state;
    }
}

#[test]
fn matrix_slip_can_exceed_free_motion_and_prescribed_or_unowned_queries_need_one_evaluation() {
    let law = motion_law(&motion_card([1.0; 2], [0.1; 2]), [1.0, 0.0]);
    let slip = [0.1, 1.0];
    let work = norm([10.0 * slip[0], slip[1]]);
    let force = [100.0 * slip[0] / work, slip[1] / work];
    let c = 0.05;
    let free = [
        slip[0] + c * (force[0] - force[1]),
        slip[1] - c * (force[0] - force[1]),
    ];
    assert!(norm(slip) > norm(free));
    let mut input = drive(free);
    input.dt_s = 1.0;
    let matrix = SymmetricTangentialCompliance::new([[c, -c], [-c, c]]).unwrap();
    let result = law
        .advance_motion_consistent(
            &law.initial_state(),
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            matrix,
            solve([0.0, 1.5]),
            || false,
        )
        .unwrap();
    for i in 0..2 {
        close(result.step.slip_increment_m[i], slip[i]);
        close(result.step.traction_n[i], force[i]);
    }
    close(result.actual_distance_m, norm(slip));
    assert!(matches!(
        law.advance_motion_consistent(
            &law.initial_state(),
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            matrix,
            solve([0.0, norm(free)]),
            || false
        ),
        Err(MotionQueryError::Unbracketed { .. })
    ));

    let prescribed = law
        .advance_motion_consistent(
            &law.initial_state(),
            &drive([0.3, 0.4]),
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            SymmetricTangentialCompliance::new([[0.0; 2]; 2]).unwrap(),
            solve([0.5; 2]),
            || false,
        )
        .unwrap();
    assert_eq!(prescribed.evaluations, 1);
    assert_eq!(prescribed.distance_residual_m, 0.0);
    assert_eq!(
        prescribed.step.state.friction_sources()[0].query().axes()["speed"],
        1.0
    );
    assert_eq!(
        prescribed.step.state.friction_sources()[0].query().axes()["path"],
        0.75
    );

    let no_motion = super::law(&directional_card(), [1.0, 0.0]);
    let old = no_motion.initial_state();
    let x = super::input(300.0, [0.3, 0.4]);
    let matrix = SymmetricTangentialCompliance::new([[0.01, 0.004], [0.004, 0.02]]).unwrap();
    let expected = no_motion
        .advance_at_query_with_compliance(
            &old,
            &x,
            &query(300.0),
            EllipticFrictionMode::Kinetic,
            matrix,
        )
        .unwrap();
    let actual = no_motion
        .advance_motion_consistent(
            &old,
            &x,
            &query(300.0),
            EllipticFrictionMode::Kinetic,
            matrix,
            solve([0.0, 1.0]),
            || false,
        )
        .unwrap();
    assert_eq!(actual.evaluations, 1);
    assert_eq!(actual.step, expected);
    assert_eq!(actual.queried_distance_m, actual.actual_distance_m);
}

#[test]
fn static_and_kinetic_branches_keep_their_actual_motion_and_zero_load_path() {
    let law = motion_law(&motion_card([0.3, 0.5], [0.2, 0.4]), [1.0, 0.0]);
    let state = law.initial_state();
    let compliance = SymmetricTangentialCompliance::new([[0.01, 0.004], [0.004, 0.02]]).unwrap();
    let input = drive([0.05, 0.02]); // C * [5 N, 0].
    let rest = law
        .advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::StaticOrSliding,
            compliance,
            solve([0.0, 0.1]),
            || false,
        )
        .unwrap();
    assert_eq!(rest.evaluations, 1);
    assert_eq!(rest.step.regime, FrictionRegime::Sticking);
    assert_eq!(rest.actual_distance_m, 0.0);
    close(rest.step.traction_n[0], 5.0);
    assert_eq!(
        rest.step.state.friction_sources()[0].query().axes()["speed"],
        0.0
    );
    assert_eq!(
        rest.step.state.friction_sources()[0].query().axes()["path"],
        0.25
    );
    let moving = law
        .advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 0.1]),
            || false,
        )
        .unwrap();
    assert_eq!(moving.step.regime, FrictionRegime::Sliding);
    assert!(moving.actual_distance_m > 0.0);
    close(
        moving.step.state.friction_sources()[0].query().axes()["speed"],
        moving.actual_distance_m / input.dt_s,
    );
    assert_eq!(
        moving.step.state.constitutive_law().static_coefficients(),
        rest.step.state.constitutive_law().static_coefficients()
    );
    let no_load = TangentialStepInput {
        normal_n: 0.0,
        ..drive([0.03, 0.04])
    };
    let free = law
        .advance_motion_consistent(
            &state,
            &no_load,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 0.1]),
            || false,
        )
        .unwrap();
    assert_eq!(free.evaluations, 1);
    assert_eq!(free.step.dissipated_j, 0.0);
    assert_eq!(free.step.slip_increment_m, [0.03, 0.04]);
    close(
        free.step.state.friction_sources()[0].query().axes()["path"],
        0.3,
    );
}

#[test]
fn bracket_source_budget_and_mid_solve_cancellation_refuse_without_partial_history() {
    let law = motion_law(&motion_card([0.3, 0.5], [0.2, 0.4]), [1.0, 0.0]);
    let state = law.initial_state();
    let original = state.clone();
    let compliance = SymmetricTangentialCompliance::new([[0.01, 0.004], [0.004, 0.02]]).unwrap();
    let input = drive([0.3, 0.4]);
    let mut budget = solve([0.0, 1.0]);
    budget.maximum_evaluations = 1;
    assert!(matches!(
        law.advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            budget,
            || false
        ),
        Err(MotionQueryError::EvaluationBudget {
            evaluations: 1,
            maximum: 1
        })
    ));
    assert!(matches!(
        law.advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 0.01]),
            || false
        ),
        Err(MotionQueryError::Unbracketed { .. })
    ));
    assert!(matches!(
        law.advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 2.0]),
            || false
        ),
        Err(MotionQueryError::Source(_))
    ));
    assert!(matches!(
        law.advance_motion_consistent(
            &state,
            &input,
            &endpoint().with("ambient", 296.0).unwrap(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 1.0]),
            || false
        ),
        Err(MotionQueryError::Source(_))
    ));
    for stop_at in [1, 3, 4] {
        let mut polls = 0;
        assert_eq!(
            law.advance_motion_consistent(
                &state,
                &input,
                &endpoint(),
                EllipticFrictionMode::Kinetic,
                compliance,
                solve([0.0, 1.0]),
                || {
                    polls += 1;
                    polls == stop_at
                }
            ),
            Err(MotionQueryError::Cancelled)
        );
    }
    let mut limited = input;
    limited.solve.maximum_force_n = 0.1;
    assert!(matches!(
        law.advance_motion_consistent(
            &state,
            &limited,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 1.0]),
            || false
        ),
        Err(MotionQueryError::Source(_))
    ));
    assert_eq!(state, original);
    let a = law
        .advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 1.0]),
            || false,
        )
        .unwrap();
    let b = law
        .advance_motion_consistent(
            &state,
            &input,
            &endpoint(),
            EllipticFrictionMode::Kinetic,
            compliance,
            solve([0.0, 1.0]),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
}
