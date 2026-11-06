//! Actual incoming slip, original interface identity and atomic source history.

#[path = "source/sustained.rs"]
mod sustained;
#[path = "source/time.rs"]
mod time;

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
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::state_point::{ScalarAdmissibility, ScalarPropertyRequirement};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::source::{
    BrittleFaceFrictionSource, BrittleFaceSourcedFrictionError, BrittleVertexFrictionSource,
};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn query() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
        )
        .unwrap()
        .with_quantity(
            "speed",
            QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])),
            0.25,
        )
        .unwrap()
        .with_quantity(
            "path",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.25,
        )
        .unwrap()
}

fn law(material: &MaterialStateId) -> SourceThermalEllipticCoulomb {
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
        (
            "kinetic-x",
            "speed",
            Dims([1, 0, -1, 0, 0, 0]),
            [0.0, 2.0],
            [0.2, 0.4],
        ),
        (
            "kinetic-y",
            "path",
            Dims([1, 0, 0, 0, 0, 0]),
            [0.0, 2.0],
            [0.1, 0.3],
        ),
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
                    source: format!("manufactured brittle-face friction {name}; not measured data"),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = || SurfaceSpec {
        material: material.clone(),
        texture_frame: "declared face tangents".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface(),
        surface(),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "manufactured".into(),
            history: "accepted sliding".into(),
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
            query(),
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
    SourceThermalEllipticCoulomb::from_principal_sources(principal, [1.0, 0.0]).unwrap()
}

fn history(law: &SourceThermalEllipticCoulomb) -> SourceThermalEllipticState {
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
        &query().with("path", 0.2525).unwrap(),
        EllipticFrictionMode::Kinetic,
    )
    .unwrap()
    .state
}

fn conditions() -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query()
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

#[test]
fn original_interface_queries_actual_slip_and_publishes_sources_with_the_native_event() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let materials = sources(&geometry);
        let solid = failed(&geometry, &materials, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let mut velocity = vec![[0.0; 3]; target.nodes_m().len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { -1.0 } else { 1.0 };
            for &node in tet {
                velocity[node] = [sign * 0.25, 0.0, -sign * SPEED];
            }
        }
        let motion = BrittleTetMotion::try_from_solid(&split.solid, &velocity, 2.0, 2, cx).unwrap();
        let law = law(materials[0].material().material());
        let history = history(&law);
        let original_history = history.clone();
        let conditions = conditions();
        let vertex = BrittleVertexFrictionSource {
            law: &law,
            state: &history,
            conditions: &conditions,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        };
        let binding = BrittleFaceFrictionSource {
            face_id: 951,
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [vertex; 3],
        };
        let hit = motion
            .impact_released_faces_sourced_frictional(
                &[0.0],
                &[binding],
                &[],
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
        retained(&motion, &hit.impact.motion);
        assert_eq!(hit.samples.len(), 3);
        let pn = SPEED / 8.0;
        for (index, sample) in hit.samples.iter().enumerate() {
            assert_eq!(sample.row, hit.impact.rows[index]);
            close(sample.sample.queried_speed_m_s, 0.5, 1e-14);
            close(sample.sample.coefficients_mu[0], 0.25, 1e-14);
            close(sample.sample.coefficients_mu[1], 0.12525, 1e-14);
            close(hit.impact.contact.impulses[index][0], pn, 1e-10);
            close(hit.impact.contact.impulses[index][1], -0.25 * pn, 1e-10);
            close(
                hit.impact.contact.relative_velocities_after[index][1],
                0.5 - 16.0 * 0.25 * pn,
                1e-9,
            );
            let accepted = sample.sample.state.constitutive_state();
            assert_eq!(accepted.slip_m(), history.constitutive_state().slip_m());
            assert_eq!(
                accepted.principal_axis(),
                history.constitutive_state().principal_axis()
            );
            assert_eq!(
                accepted.contact_temperature_k(),
                history.constitutive_state().contact_temperature_k()
            );
            for source in sample.sample.state.friction_sources() {
                close(source.state_coordinate("speed").unwrap(), 0.5, 1e-14);
                close(source.state_coordinate("path").unwrap(), 0.2525, 1e-14);
            }
        }
        assert_eq!(history, original_history);

        // Distinct physical slips on the same face must refresh each local
        // source separately, without a face-average or caller placeholder.
        let face = split.solid.released_faces()[0];
        let mut varying = velocity.clone();
        for vertex in 0..3 {
            varying[face.paired_nodes[1][vertex]][0] =
                varying[face.paired_nodes[0][vertex]][0] + 0.25 * (vertex + 1) as f64;
        }
        let varying_motion =
            BrittleTetMotion::try_from_solid(&split.solid, &varying, 2.0, 2, cx).unwrap();
        let varying_hit = varying_motion
            .impact_released_faces_sourced_frictional(
                &[0.0],
                &[binding],
                &[],
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
        for sample in &varying_hit.samples {
            let speed = 0.25 * (sample.row.vertex + 1) as f64;
            close(sample.sample.queried_speed_m_s, speed, 1e-14);
            close(sample.sample.coefficients_mu[0], 0.2 + 0.1 * speed, 1e-14);
        }

        // A late source failure follows successful candidates for the first
        // vertices. Neither those histories nor the body may be published.
        let mut bad_conditions = conditions.clone();
        bad_conditions.contact_temperature_k = Some(500.0);
        bad_conditions.query = conditions.query.clone().with("T", 500.0).unwrap();
        let mut bad = binding;
        bad.vertices[2].conditions = &bad_conditions;
        assert!(
            matches!(motion.impact_released_faces_sourced_frictional(&[0.0], &[bad], &[], ORIGIN, controls(), cx), Err(BrittleFaceSourcedFrictionError::Source { row, .. }) if row.vertex == 2)
        );
        let all_held: Vec<_> = (0..velocity.len() * 3).collect();
        assert!(
            motion
                .impact_released_faces_sourced_frictional(
                    &[0.0],
                    &[binding],
                    &all_held,
                    ORIGIN,
                    controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(history, original_history);
        let retry = motion
            .impact_released_faces_sourced_frictional(
                &[0.0],
                &[binding],
                &[],
                ORIGIN,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            retry.impact.motion.velocity_mps(),
            hit.impact.motion.velocity_mps()
        );
        assert_eq!(retry.samples[0].sample, hit.samples[0].sample);

        let mut unrelated = materials[0].material().material().clone();
        unrelated.revision += 1;
        let wrong_law = self::law(&unrelated);
        let wrong_history = wrong_law.initial_state();
        let wrong = BrittleFaceFrictionSource {
            vertices: [BrittleVertexFrictionSource {
                law: &wrong_law,
                state: &wrong_history,
                ..vertex
            }; 3],
            ..binding
        };
        assert!(matches!(
            motion.impact_released_faces_sourced_frictional(
                &[0.0],
                &[wrong],
                &[],
                ORIGIN,
                controls(),
                cx
            ),
            Err(BrittleFaceSourcedFrictionError::Binding { .. })
        ));

        // One open vertex has deliberately invalid external conditions. It is
        // retained geometrically but must not query or refresh its source.
        let face = split.solid.released_faces()[0];
        let mut opened = split.solid.solution().displacement_m.clone();
        opened[face.paired_nodes[1][2]][2] += 1.0 / 1024.0;
        let open = split
            .solid
            .advance(&prescribed(&opened), &[], settings(), work_controls(), cx)
            .unwrap()
            .solid;
        let open_motion = BrittleTetMotion::try_from_solid(&open, &velocity, 2.0, 2, cx).unwrap();
        let active = open_motion
            .impact_released_faces_sourced_frictional(&[0.0], &[bad], &[], ORIGIN, controls(), cx)
            .unwrap();
        assert_eq!(active.samples.len(), 2);
        assert!(active.samples.iter().all(|sample| sample.row.vertex != 2));
        assert_eq!(history, original_history);
    });
}
