//! Native material queries and complete source/mechanical acceptance together.

use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceFrictionConditions, SourceImpulseCoefficientRole,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    EllipticFrictionMode, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solid::dynamics::frictional_impact::elliptic::source::{
    SolidEllipticSourceBinding, SolidEllipticSourceError, SolidImpactGeometry,
    SolidImpactKinematics,
};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn point() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
        )
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.25)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap()
}

fn material(axis: [f64; 2]) -> (InterfaceSystemCard, SourceThermalEllipticCoulomb) {
    let mut claims = ClaimSet::new();
    for (name, coordinate, dims, range, values) in [
        (
            "static-x",
            "T",
            Temperature::DIMS,
            [280.0, 320.0],
            [0.7, 0.5],
        ),
        (
            "static-y",
            "T",
            Temperature::DIMS,
            [280.0, 320.0],
            [0.8, 0.6],
        ),
        ("kinetic-x", "speed", SPEED, [0.0, 2.0], [0.2, 0.4]),
        ("kinetic-y", "path", LENGTH, [0.0, 2.0], [0.1, 0.3]),
    ] {
        let validity = if coordinate == "T" {
            ValidityDomain::unconstrained().with("T", range[0], range[1])
        } else {
            ValidityDomain::unconstrained().with_quantity(
                coordinate,
                QuantitySpec::dimensional(dims),
                range[0],
                range[1],
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: coordinate.into(),
                    abscissa_dims: dims,
                    knots: vec![(range[0], values[0]), (range[1], values[1])],
                    dims: Dims::NONE,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: format!(
                        "synthetic solid source-friction fixture {name}; not measured data"
                    ),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic impact fixture".into(),
            revision: 0,
        },
        texture_frame: format!("texture-{name}"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic fixture".into(),
            history: "accepted sliding specimen".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let principal = [["static-x", "kinetic-x"], ["static-y", "kinetic-y"]].map(|names| {
        let requirement = |name| {
            ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
                .unwrap()
        };
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            point(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: Some("p".into()),
            },
            vec![
                FrictionQueryAxis {
                    name: "speed".into(),
                    kind: FrictionQueryAxisKind::SlidingSpeed,
                },
                FrictionQueryAxis {
                    name: "path".into(),
                    kind: FrictionQueryAxisKind::AccumulatedDistance,
                },
            ],
        )
        .unwrap()
    });
    let law = SourceThermalEllipticCoulomb::from_principal_sources(principal, axis).unwrap();
    (card, law)
}

fn accepted(law: &SourceThermalEllipticCoulomb) -> SourceThermalEllipticState {
    // An actual finite sliding interval precedes the impact. Its independent
    // 0.0025 m travel at 0.25 m/s is consistent with the owned query update.
    law.advance_at_query(
        &law.initial_state(),
        &TangentialStepInput {
            normal_n: 10.0,
            slip_increment_m: [0.0025, 0.0],
            dt_s: 0.01,
            compliance_m_per_n: 0.0,
            contact_temperature_k: Some(300.0),
            nominal_pressure_pa: Some(1e5),
            solve: TangentialSolveConfig {
                max_iterations: 64,
                maximum_force_n: 100.0,
                force_absolute_tolerance_n: 1e-12,
                force_relative_tolerance: 1e-12,
            },
        },
        &point().with("path", 0.2525).unwrap(),
        EllipticFrictionMode::Kinetic,
    )
    .unwrap()
    .state
}

fn same_physical_history(actual: &SourceThermalEllipticState, old: &SourceThermalEllipticState) {
    // A coefficient refresh changes the checkpoint's private law owner, while
    // the actual physical history must stay exactly at the accepted endpoint.
    let (actual, old) = (actual.constitutive_state(), old.constitutive_state());
    assert_eq!(actual.slip_m(), old.slip_m());
    assert_eq!(actual.principal_axis(), old.principal_axis());
    assert_eq!(actual.contact_temperature_k(), old.contact_temperature_k());
}

fn conditions() -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: point()
            .with("T", 320.0)
            .unwrap()
            .with("p", 2e5)
            .unwrap()
            .with("speed", 99.0)
            .unwrap()
            .with("path", 99.0)
            .unwrap(),
        contact_temperature_k: Some(320.0),
        nominal_pressure_pa: Some(2e5),
    }
}

fn geometry(c: SolidFrictionalContact<'_>) -> SolidImpactGeometry<'_> {
    SolidImpactGeometry {
        id: c.id,
        velocity_weights: c.velocity_weights,
        gap_m: c.gap_m,
        activation_tolerance_m: c.activation_tolerance_m,
        prescribed_surface_velocity_mps: c.prescribed_surface_velocity_mps,
        normal_law: c.normal_law,
        friction_coefficient: (),
    }
}

#[test]
fn actual_relative_speed_explicit_roles_and_pinned_material_history_drive_both_solid_owners() {
    with_cx(false, |cx| {
        let (card, law) = material([1.0, 0.0]);
        let history = accepted(&law);
        assert!(history.constitutive_state().slip_m()[0] > 0.0);
        assert_eq!(
            history.constitutive_state().contact_temperature_k(),
            Some(300.0)
        );
        let original = history.clone();
        let previous = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
        let current = TangentFrame::new([0.0, 0.0, 1.0], [0.0, -1.0, 0.0]).unwrap();
        let rotated = law
            .transport_state(
                &history,
                TangentialHistoryRotation::between(previous, current, 0.0).unwrap(),
            )
            .unwrap();
        assert_eq!(rotated.constitutive_state().principal_axis(), [0.0, 1.0]);
        let good = conditions();
        let outside = SourceFrictionConditions {
            query: good.query.clone().with("T", 999.0).unwrap(),
            contact_temperature_k: Some(999.0),
            ..good.clone()
        };
        let identity: [f64; 36] = core::array::from_fn(|i| if i / 6 == i % 6 { 1.0 } else { 0.0 });
        let pencil = assembly(&[0.0; 36], &identity, 6);
        let initial = [-0.8, 1.4, 0.3, -1.0, 40.0, 30.0];
        let contacts = [
            geometry(contact(
                40,
                [&identity[0..6], &identity[6..12], &identity[12..18]],
                [0.2, 0.4, 0.3],
                0.0,
            )),
            SolidImpactGeometry {
                gap_m: 1.0,
                ..geometry(contact(
                    7,
                    [&identity[18..24], &identity[24..30], &identity[30..36]],
                    [0.0; 3],
                    0.0,
                ))
            },
        ];
        // Both admitted histories are already expressed in the contact rows'
        // tangent coordinates. Using the law's initial axis instead of the
        // transported checkpoint would give the wrong rotated friction limit.
        for (history, role, mu) in [
            (&original, SourceImpulseCoefficientRole::Kinetic, 0.3),
            (&original, SourceImpulseCoefficientRole::Static, 0.5),
            (&rotated, SourceImpulseCoefficientRole::Kinetic, 0.12525),
            (&rotated, SourceImpulseCoefficientRole::Static, 0.6),
        ] {
            let sources = [
                SolidEllipticSourceBinding {
                    contact_id: 40,
                    law: &law,
                    conditions: &good,
                    coefficient_role: role,
                },
                SolidEllipticSourceBinding {
                    contact_id: 7,
                    law: &law,
                    conditions: &outside,
                    coefficient_role: role,
                },
            ];
            let mut dense =
                ExactLinearSolid::prepare(&pencil, 0.02, &[0.0; 6], &initial, dense_controls(), cx)
                    .unwrap();
            let mut spectral = CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 6],
                &identity,
                0.02,
                &[0.0; 6],
                &initial,
                spectral_controls(),
                cx,
            )
            .unwrap();
            let (mut hd, mut hs) = (
                [history.clone(), original.clone()],
                [history.clone(), original.clone()],
            );
            let d = dense
                .apply_sourced_elliptic_impacts(&contacts, &mut hd, &sources, impact_controls(), cx)
                .unwrap();
            let s = spectral
                .apply_sourced_elliptic_impacts(&contacts, &mut hs, &sources, impact_controls(), cx)
                .unwrap();
            assert_eq!(hd, hs);
            assert_eq!(hd[1], original); // Inactive source is outside its domain and is never queried.
            for (velocity, report) in [(dense.velocity_mps(), &d), (spectral.velocity_mps(), &s)] {
                near(velocity[0], 0.2, 3e-9);
                near(velocity[1], 1.4 - mu, 3e-9);
                near(velocity[2], 0.3, 3e-9);
                assert_eq!(&velocity[3..], &initial[3..]);
                near(report.contact.law.impulses[0][0], 1.0, 3e-9);
                near(report.contact.law.impulses[0][1], -mu, 3e-9);
                assert_eq!(report.contact.law.impulses[1], [0.0; 3]);
                near(
                    report.contact.relative_work_j,
                    -0.5 - mu + 0.5 * mu * mu,
                    3e-9,
                );
                near(
                    report.contact.prescribed_surface_work_j,
                    0.2 - 0.4 * mu,
                    3e-9,
                );
                assert_eq!(report.impulse.energy.dissipated_j, 0.0);
                assert_eq!(report.impulse.time_s, 0.0);
            }
            let next = &hd[0];
            same_physical_history(next, history);
            for (index, source) in next.friction_sources().iter().enumerate() {
                near(source.state_coordinate("speed").unwrap(), 1.0, 1e-14);
                assert_eq!(source.state_coordinate("T"), Some(320.0));
                assert_eq!(source.state_coordinate("p"), Some(2e5));
                assert_eq!(
                    source.state_coordinate("path"),
                    original.friction_sources()[index].state_coordinate("path")
                );
                for property in source.properties() {
                    let previous = original.friction_sources()[index]
                        .property(property.requirement().name())
                        .unwrap();
                    assert_eq!(
                        property.answer().receipt.selected,
                        previous.answer().receipt.selected
                    );
                    card.claims()
                        .verify_receipt(&property.answer().receipt)
                        .unwrap();
                }
            }
            for (actual, expected) in next
                .constitutive_law()
                .static_coefficients()
                .into_iter()
                .zip([0.5, 0.6])
            {
                near(actual, expected, 1e-14);
            }
            for (actual, expected) in next
                .constitutive_law()
                .kinetic_coefficients()
                .into_iter()
                .zip([0.3, 0.12525])
            {
                near(actual, expected, 1e-14);
            }
        }
    });
}

// No Clone requirement: only returned candidates may change these accepted states.
#[derive(Debug, PartialEq)]
struct History {
    contact_id: u64,
    accepted: u32,
}

fn refresh(old: &History, query: &SolidImpactKinematics) -> (EllipticFriction, History) {
    assert_eq!(old.contact_id, query.contact_id);
    (
        EllipticFriction {
            coefficients_mu: [0.2, 0.4],
            principal_axis: [1.0, 0.0],
        },
        History {
            contact_id: old.contact_id,
            accepted: old.accepted + 1,
        },
    )
}

#[test]
fn callback_order_indices_candidate_failure_and_invalid_activation_are_transactional() {
    with_cx(false, |cx| {
        let identity: [f64; 36] = core::array::from_fn(|i| if i / 6 == i % 6 { 1.0 } else { 0.0 });
        let pencil = assembly(&[0.0; 36], &identity, 6);
        let make = || {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 6],
                &identity,
                0.02,
                &[0.0; 6],
                &[-1.0, 0.8, 0.0, -0.5, 0.0, 0.9],
                spectral_controls(),
                cx,
            )
            .unwrap()
        };
        let (mut solid, mut permuted) = (make(), make());
        let contacts = [
            geometry(contact(
                40,
                [&identity[0..6], &identity[6..12], &identity[12..18]],
                [0.0; 3],
                0.0,
            )),
            geometry(contact(
                7,
                [&identity[18..24], &identity[24..30], &identity[30..36]],
                [0.0; 3],
                0.0,
            )),
        ];
        let mut history = [
            History {
                contact_id: 40,
                accepted: 0,
            },
            History {
                contact_id: 7,
                accepted: 0,
            },
        ];
        let before = bits(&solid);
        let mut queries = vec![];
        let failure = solid
            .apply_elliptic_impacts_with(
                &contacts,
                &mut history,
                impact_controls(),
                |old, query| {
                    queries.push((query.contact_id, query.original_index));
                    if query.contact_id == 40 {
                        Err("unavailable second source")
                    } else {
                        Ok(refresh(old, query))
                    }
                },
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            failure,
            SolidEllipticSourceError::Source {
                contact_id: 40,
                original_index: 0,
                source: "unavailable second source"
            }
        ));
        assert_eq!(queries, [(7, 1), (40, 0)]);
        assert_eq!(
            history.iter().map(|h| h.accepted).collect::<Vec<_>>(),
            [0, 0]
        );
        assert_eq!(bits(&solid), before);
        for (gap, tolerance) in [
            (f64::NAN, 0.0),
            (f64::INFINITY, 0.0),
            (-1.0, 0.0),
            (1.0, f64::NAN),
            (1.0, f64::INFINITY),
            (1.0, -1.0),
        ] {
            let mut invalid = contacts;
            invalid[0].gap_m = gap;
            invalid[0].activation_tolerance_m = tolerance;
            let mut calls = 0;
            assert!(
                solid
                    .apply_elliptic_impacts_with(
                        &invalid,
                        &mut history,
                        impact_controls(),
                        |old, query| {
                            calls += 1;
                            Ok::<_, &'static str>(refresh(old, query))
                        },
                        cx
                    )
                    .is_err()
            );
            assert_eq!(calls, 0);
            assert_eq!(bits(&solid), before);
        }
        let mut other = [
            History {
                contact_id: 7,
                accepted: 0,
            },
            History {
                contact_id: 40,
                accepted: 0,
            },
        ];
        solid
            .apply_elliptic_impacts_with(
                &contacts,
                &mut history,
                impact_controls(),
                |old, query| Ok::<_, &'static str>(refresh(old, query)),
                cx,
            )
            .unwrap();
        permuted
            .apply_elliptic_impacts_with(
                &[contacts[1], contacts[0]],
                &mut other,
                impact_controls(),
                |old, query| Ok::<_, &'static str>(refresh(old, query)),
                cx,
            )
            .unwrap();
        assert_eq!(bits(&solid), bits(&permuted));
        assert_eq!(history[0], other[1]);
        assert_eq!(history[1], other[0]);
        assert!(history.iter().all(|h| h.accepted == 1));
    });
}

#[test]
fn material_refusal_and_late_native_refusal_publish_neither_history_nor_motion() {
    with_cx(false, |cx| {
        let (_, law) = material([1.0, 0.0]);
        let old = accepted(&law);
        let good = conditions();
        let bad = SourceFrictionConditions {
            nominal_pressure_pa: Some(1e5),
            ..good.clone()
        };
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        let make = || {
            let mut controls = spectral_controls();
            controls.maximum_relative_operator_error = 1e-3;
            controls.maximum_relative_state_error = 1e-3;
            controls.maximum_work_error_j = 1e-3;
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 3],
                &[1.0001, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                0.05,
                &[0.0; 3],
                &[-0.2, 0.3, 0.4],
                controls,
                cx,
            )
            .unwrap()
        };
        let (mut solid, mut reference) = (make(), make());
        let contacts = [geometry(contact(
            5,
            [&identity[0..3], &identity[3..6], &identity[6..9]],
            [0.0; 3],
            0.0,
        ))];
        let binding = SolidEllipticSourceBinding {
            contact_id: 5,
            law: &law,
            conditions: &good,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        };
        let mut history = [old.clone()];
        let before = bits(&solid);
        let mut controls = impact_controls();
        controls.maximum_impulse_error_kg_m_s = 1e-5;
        controls.maximum_work_error_j = 1e-3;
        assert!(matches!(
            solid.apply_sourced_elliptic_impacts(
                &contacts,
                &mut history,
                &[SolidEllipticSourceBinding {
                    conditions: &bad,
                    ..binding
                }],
                controls,
                cx
            ),
            Err(SolidEllipticSourceError::Source {
                contact_id: 5,
                original_index: 0,
                ..
            })
        ));
        assert_eq!(history[0], old);
        assert!(matches!(
            solid.apply_sourced_elliptic_impacts(
                &contacts,
                &mut history,
                &[SolidEllipticSourceBinding {
                    contact_id: 9,
                    ..binding
                }],
                controls,
                cx
            ),
            Err(SolidEllipticSourceError::Impact(
                SolidFrictionalImpactError::Invalid(_)
            ))
        ));
        assert!(matches!(
            solid.apply_sourced_elliptic_impacts(&contacts, &mut history, &[binding], controls, cx),
            Err(SolidEllipticSourceError::Impact(
                SolidFrictionalImpactError::Speed { .. }
            ))
        ));
        assert_eq!(bits(&solid), before);
        assert_eq!(history[0], old);
        controls.maximum_contact_speed_error_mps = 1e-3;
        assert!(matches!(
            solid.apply_sourced_elliptic_impacts(&contacts, &mut history, &[binding], controls, cx),
            Err(SolidEllipticSourceError::Impact(
                SolidFrictionalImpactError::EllipticContact(
                    EllipticFrictionalContactError::ResidualNotSatisfied(_)
                )
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.apply_sourced_elliptic_impacts(
                    &contacts,
                    &mut history,
                    &[binding],
                    controls,
                    cancelled
                ),
                Err(SolidEllipticSourceError::Impact(
                    SolidFrictionalImpactError::Cancelled
                ))
            ));
        });
        assert_eq!(bits(&solid), before);
        assert_eq!(history[0], old);
        controls.contact.velocity_tolerance = 1e-3;
        controls.contact.complementarity_tolerance = 1e-4;
        controls.contact.energy_tolerance = 1e-4;
        let mut expected_history = [old.clone()];
        solid
            .apply_sourced_elliptic_impacts(&contacts, &mut history, &[binding], controls, cx)
            .unwrap();
        reference
            .apply_sourced_elliptic_impacts(
                &contacts,
                &mut expected_history,
                &[binding],
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(bits(&solid), bits(&reference));
        assert_eq!(history, expected_history);
        same_physical_history(&history[0], &old);
        for source in history[0].friction_sources() {
            assert_eq!(source.state_coordinate("T"), Some(320.0));
            near(source.state_coordinate("speed").unwrap(), 0.5, 1e-14);
        }
    });
}
