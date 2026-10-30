//! G1/G3: original directional source claims, texture and physical work together.
#[path = "elliptic/motion.rs"]
mod motion;

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::EllipticFrictionMode;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::SymmetricTangentialCompliance;
use fs_tribo::tangential::transport::{FramedTangentialState, TangentialTransportConfig};

const STATIC_Y: &str = "second-principal-static";
const KINETIC_Y: &str = "second-principal-kinetic";

#[test]
fn full_compliance_refreshes_original_sources_at_an_independent_elliptic_solution() {
    let card = directional_card();
    let law = law(&card, [0.6, 0.8]);
    let initial = law.initial_state();
    let mut state = initial.clone();
    let matrix = [[0.01, 0.003], [0.003, 0.02]];
    let compliance = SymmetricTangentialCompliance::new(matrix).unwrap();
    let world = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    // Construct a boundary force and its associated slip independently, then
    // add the actual matrix recoil to obtain the released displacement.
    for (temperature, caps) in [(280.0, [6.0, 1.0]), (320.0, [2.0, 5.0])] {
        let force = world([0.6 * caps[0], 0.8 * caps[1]]);
        let slip = world([0.04 * 0.6 / caps[0], 0.04 * 0.8 / caps[1]]);
        let free = [
            slip[0] + matrix[0][0] * force[0] + matrix[0][1] * force[1],
            slip[1] + matrix[1][0] * force[0] + matrix[1][1] * force[1],
        ];
        let step = law
            .advance_at_query_with_compliance(
                &state,
                &input(temperature, free),
                &query(temperature),
                EllipticFrictionMode::Kinetic,
                compliance,
            )
            .unwrap();
        for axis in 0..2 {
            close(step.traction_n[axis], force[axis]);
            close(step.slip_increment_m[axis], slip[axis]);
            close(step.compatibility_residual_m[axis], 0.0);
        }
        close(step.dissipated_j, 0.04);
        for (axis, source) in step.state.friction_sources().iter().enumerate() {
            assert_eq!(source.state_coordinate("T"), Some(temperature));
            for property in source.properties() {
                let original = initial.friction_sources()[axis]
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
        state = step.state;
    }
}

#[test]
fn matrix_source_branch_and_refusals_preserve_the_original_checkpoint() {
    let card = directional_card();
    let law = law(&card, [1.0, 0.0]);
    let state = law.initial_state();
    let original = state.clone();
    let compliance = SymmetricTangentialCompliance::new([[0.01, 0.004], [0.004, 0.02]]).unwrap();
    let mut drive = input(300.0, [0.05, 0.02]);
    let advance = |drive: &TangentialStepInput, point: &QueryPoint, mode| {
        law.advance_at_query_with_compliance(&state, drive, point, mode, compliance)
    };
    let rest = advance(&drive, &query(300.0), EllipticFrictionMode::StaticOrSliding).unwrap();
    assert_eq!(rest.regime, FrictionRegime::Sticking);
    close(rest.traction_n[0], 5.0);
    close(rest.traction_n[1], 0.0);
    let moving = advance(&drive, &query(300.0), EllipticFrictionMode::Kinetic).unwrap();
    assert_eq!(moving.regime, FrictionRegime::Sliding);
    assert!(moving.dissipated_j > 0.0);
    assert_eq!(
        moving.state.constitutive_law().static_coefficients(),
        rest.state.constitutive_law().static_coefficients()
    );
    drive.solve.maximum_force_n = 1.0;
    assert!(advance(&drive, &query(300.0), EllipticFrictionMode::StaticOrSliding).is_err());
    drive = input(321.0, [0.05, 0.02]);
    assert!(advance(&drive, &query(321.0), EllipticFrictionMode::Kinetic).is_err());
    drive = input(300.0, [0.05, 0.02]);
    drive.compliance_m_per_n = 0.01;
    assert!(matches!(
        advance(&drive, &query(300.0), EllipticFrictionMode::Kinetic),
        Err(SourceThermalCoulombError::Invalid(
            "matrix owns the sole directional body compliance"
        ))
    ));
    assert_eq!(state, original);
}

fn directional_card() -> InterfaceSystemCard {
    card(vec![
        curve(STATIC, [0.8, 0.4]),
        curve(KINETIC, [0.6, 0.2]),
        curve(STATIC_Y, [0.5, 0.9]),
        curve(KINETIC_Y, [0.1, 0.5]),
    ])
}

fn principal(
    card: &InterfaceSystemCard,
    point: QueryPoint,
    mapped: ThermalFrictionAxes,
) -> [SourceThermalCoulomb; 2] {
    [[STATIC, KINETIC], [STATIC_Y, KINETIC_Y]].map(|names| {
        SourceThermalCoulomb::resolve(
            bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
            point.clone(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            mapped.clone(),
        )
        .unwrap()
    })
}

fn law(card: &InterfaceSystemCard, direction: [f64; 2]) -> SourceThermalEllipticCoulomb {
    SourceThermalEllipticCoulomb::from_principal_sources(
        principal(card, query(300.0), axes()),
        direction,
    )
    .unwrap()
}

fn expected_force(caps: [f64; 2], slip: [f64; 2]) -> ([f64; 2], f64) {
    let work = (caps[0] * slip[0]).hypot(caps[1] * slip[1]);
    (
        [
            caps[0] * caps[0] * slip[0] / work,
            caps[1] * caps[1] * slip[1] / work,
        ],
        work,
    )
}

#[test]
fn temperature_refresh_changes_both_principal_laws_without_erasing_source_or_slip() {
    let card = directional_card();
    let law = law(&card, [1.0, 0.0]);
    let mut state = law.initial_state();
    for (temperature, slip, caps) in [
        (280.0, [0.003, 0.004], [6.0, 1.0]),
        (320.0, [-0.006, -0.008], [2.0, 5.0]),
    ] {
        let step = law
            .advance_sourced(&state, &input(temperature, slip))
            .unwrap();
        let (force, work) = expected_force(caps, slip);
        for i in 0..2 {
            close(step.traction_n[i], force[i]);
        }
        close(step.dissipated_j, work);
        close(step.closure_residual_j, 0.0);
        assert_eq!(step.stored_change_j, 0.0);
        for (axis, source) in step.state.friction_sources().iter().enumerate() {
            assert_eq!(source.state_coordinate("T"), Some(temperature));
            for property in source.properties() {
                let original = law.principal_sources()[axis]
                    .initial_source()
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
        state = step.state;
    }
    close(state.constitutive_state().slip_m()[0], -0.003);
    close(state.constitutive_state().slip_m()[1], -0.004);
    let rest = law.advance(&state, &input(300.0, [0.0; 2])).unwrap();
    assert_eq!(
        rest.state.constitutive_state().slip_m(),
        state.constitutive_state().slip_m()
    );
    assert_eq!(rest.dissipated_j, 0.0);
    assert_eq!(rest.state.constitutive_state().principal_axis(), [1.0, 0.0]);
}

#[test]
fn original_texture_and_static_ellipse_control_breakaway_and_survive_frame_changes() {
    let card = directional_card();
    let law = law(&card, [1.0, 0.0]);
    let rotated = self::law(&card, [0.0, 1.0]);
    let mut x = input(300.0, [0.065, 0.0]);
    x.compliance_m_per_n = 0.01;
    let sliding = law.advance(&law.initial_state(), &x).unwrap();
    let sticking = rotated.advance(&rotated.initial_state(), &x).unwrap();
    assert_eq!(sliding.regime, FrictionRegime::Sliding);
    close(sliding.traction_n[0], 4.0);
    assert_eq!(sticking.regime, FrictionRegime::Sticking);
    close(sticking.traction_n[0], 6.5);
    assert_eq!(sticking.dissipated_j, 0.0);
    let frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let retained = FramedTangentialState::new(&law, sliding.state.clone(), frame).unwrap();
    let moved = retained
        .transport(
            &law,
            frame.rotated(std::f64::consts::FRAC_PI_2).unwrap(),
            0.0,
            TangentialTransportConfig {
                absolute_energy_tolerance_j: 0.0,
                relative_energy_tolerance: 0.0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        moved.checkpoint.state().friction_sources(),
        sliding.state.friction_sources()
    );
    let warm = law
        .advance(moved.checkpoint.state(), &input(320.0, [0.0, -0.01]))
        .unwrap();
    close(warm.traction_n[0], 0.0);
    close(warm.traction_n[1], -2.0);
    close(warm.dissipated_j, 0.02);
    close(warm.state.constitutive_state().slip_m()[0], 0.0);
    close(warm.state.constitutive_state().slip_m()[1], -0.035);
}

#[test]
fn original_pressure_sources_refresh_both_axes_and_keep_their_order() {
    let card = card(vec![
        pressure_surface(STATIC, 0.0),
        pressure_surface(KINETIC, 0.2),
        pressure_surface(STATIC_Y, 0.1),
        pressure_surface(KINETIC_Y, 0.3),
    ]);
    let point = query(300.0)
        .with_quantity("p", pressure_quantity(), 1e5)
        .unwrap();
    let sources = principal(
        &card,
        point,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: Some("p".into()),
        },
    );
    let law = SourceThermalEllipticCoulomb::from_principal_sources(sources, [1.0, 0.0]).unwrap();
    let mut state = law.initial_state();
    for (p, caps) in [(1e5, [4.0, 3.0]), (3e5, [6.0, 5.0])] {
        let mut x = input(300.0, [0.01, -0.02]);
        x.nominal_pressure_pa = Some(p);
        let step = law.advance(&state, &x).unwrap();
        let (force, work) = expected_force(caps, x.slip_increment_m);
        for i in 0..2 {
            close(step.traction_n[i], force[i]);
        }
        close(step.dissipated_j, work);
        assert!(
            step.state
                .friction_sources()
                .iter()
                .all(|s| s.state_coordinate("p") == Some(p))
        );
        state = step.state;
    }
    assert!(law.advance(&state, &input(300.0, [0.01, 0.0])).is_err());
}

#[test]
fn mismatched_sources_and_single_axis_domain_failures_cannot_publish_partial_history() {
    let card = directional_card();
    let sources = principal(&card, query(300.0), axes());
    let mut different = principal(&card, query(301.0), axes());
    different[0] = sources[0].clone();
    assert!(SourceThermalEllipticCoulomb::from_principal_sources(different, [1.0, 0.0]).is_err());
    let original =
        SourceThermalEllipticCoulomb::from_principal_sources(sources.clone(), [1.0, 0.0]).unwrap();
    let swapped = SourceThermalEllipticCoulomb::from_principal_sources(
        [sources[1].clone(), sources[0].clone()],
        [1.0, 0.0],
    )
    .unwrap();
    let state = original
        .advance(&original.initial_state(), &input(300.0, [0.01, 0.02]))
        .unwrap()
        .state;
    assert!(
        swapped
            .advance_sourced(&state, &input(300.0, [0.01, 0.02]))
            .is_err()
    );
    let mut limited = curve(KINETIC_Y, [0.1, 0.5]);
    limited.validity = ValidityDomain::unconstrained().with("T", 280.0, 310.0);
    let restricted_card = super::card(vec![
        curve(STATIC, [0.8, 0.4]),
        curve(KINETIC, [0.6, 0.2]),
        curve(STATIC_Y, [0.5, 0.9]),
        limited,
    ]);
    let restricted = law(&restricted_card, [1.0, 0.0]);
    let accepted = restricted
        .advance(&restricted.initial_state(), &input(300.0, [0.01, 0.02]))
        .unwrap()
        .state;
    let saved = accepted.clone();
    assert!(matches!(
        restricted.advance_sourced(&accepted, &input(320.0, [0.01, 0.02])),
        Err(SourceThermalCoulombError::Material(_))
    ));
    assert_eq!(accepted, saved);
    let next = input(305.0, [0.02, -0.01]);
    assert_eq!(
        restricted.advance(&accepted, &next).unwrap(),
        restricted.advance(&saved, &next).unwrap()
    );
    let mut budget = next;
    budget.solve.maximum_force_n = 0.1;
    assert!(restricted.advance(&accepted, &budget).is_err());
    assert_eq!(accepted, saved);
}
