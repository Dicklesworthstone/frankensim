//! Local source laws on exact matching traces, sharing the existing contact owners.

#[path = "pointwise_matching_contact/transient.rs"]
mod transient;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::interface::temperature::{
    TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair, ScalarField,
    ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;

fn context<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 37,
                    kernel_id: 822,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
}

fn mesh(bent: bool) -> ConductionMesh {
    let mut positions = vec![
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., -1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ];
    let mut tets = vec![[0, 1, 2, 3], [4, 5, 6, 7]];
    if bent {
        let rotated: Vec<_> = positions.iter().map(|p| [p[2] + 3.0, p[1], p[0]]).collect();
        positions.extend(rotated);
        tets.extend([[8, 9, 10, 11], [12, 13, 14, 15]]);
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}

fn kelvin() -> fs_qty::QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("TA", kelvin(), 320.)
        .unwrap()
        .with_quantity("TB", kelvin(), 320.)
        .unwrap()
        .with("normal_pressure", 1e6)
        .unwrap()
}

fn resistance(a: f64, b: f64) -> f64 {
    0.5 + (a - 320.) / 1024. + 3. * (b - 320.) / 1024.
}

fn source_card(reverse: bool) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                dims,
            ),
            value: PropertyValue::Surface2D {
                axes: ["TA".into(), "TB".into()],
                axis_dims: [fs_qty::Temperature::DIMS; 2],
                samples: [[256., 256.], [384., 256.], [256., 384.], [384., 384.]]
                    .into_iter()
                    .map(|p| (p, resistance(p[0], p[1])))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("TA", kelvin(), 256., 384.)
                .with_quantity("TB", kelvin(), 256., 384.)
                .with("normal_pressure", 1e6, 1e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured pointwise contact law; not measured material data".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "fixed reference normal".into(),
    };
    let names = if reverse { ["B", "A"] } else { ["A", "B"] };
    InterfaceSystemCard::assemble(
        side(names[0]),
        side(names[1]),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "test".into(),
            history: "unaged".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn law<'c>(
    mesh: &ConductionMesh,
    card: &'c InterfaceSystemCard,
    reverse: bool,
    pointwise: bool,
) -> TemperatureContact<'c> {
    let mut pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    if reverse {
        pairs = pairs
            .into_iter()
            .map(|p| InterfaceFacePair {
                side_a: p.side_b,
                side_b: p.side_a,
            })
            .collect();
    }
    let axes = if reverse { ["TB", "TA"] } else { ["TA", "TB"] };
    if pointwise {
        TemperatureContact::new_pointwise(
            "bond",
            pairs,
            card,
            &point(),
            TemperatureContactCoordinates::TwoSided(axes.map(String::from)),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap()
    } else {
        TemperatureContact::new_two_sided(
            "bond",
            pairs,
            card,
            &point(),
            axes,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap()
    }
}

fn contacts<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    boundary: &ThermalBoundary,
    card: &'c InterfaceSystemCard,
    reverse: bool,
) -> TemperatureContacts<'m, 'c> {
    let pairs = mesh.element_count() / 2;
    TemperatureContacts::with_matching_quadrature(
        cx,
        mesh,
        boundary,
        vec![],
        vec![law(mesh, card, reverse, true)],
        pairs,
        3 * pairs,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.),
        "{actual:e} != {expected:e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Independent physical triangle rule in the known original vertex order;
// no production geometry, frozen resistance or query result supplies a value.
fn independent_contact_load(t: &[f64]) -> Vec<f64> {
    let mut result = vec![0.; 8];
    for q in 0..3 {
        let basis: [f64; 3] = std::array::from_fn(|i| if i == q { 2. / 3. } else { 1. / 6. });
        let a = (0..3).map(|i| basis[i] * t[i]).sum::<f64>();
        let b = (0..3).map(|i| basis[i] * t[i + 4]).sum::<f64>();
        let heat = (a - b) / resistance(a, b) / 6.;
        for i in 0..3 {
            result[i] += basis[i] * heat;
            result[i + 4] -= basis[i] * heat;
        }
    }
    result
}

#[test]
fn local_hot_spots_change_heat_even_when_face_mean_jump_is_zero() {
    context(|_, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let t = [380., 290., 290., 320., 320., 320., 320., 320.];
        let linear = set.linearize_temperature(cx, &boundary, &t, 1).unwrap();
        let expected = independent_contact_load(&t);
        for (&a, &b) in linear.nodal_outflow_w().iter().zip(&expected) {
            near(a, b, 2e-12);
        }
        near(linear.nodal_outflow_w().iter().sum(), 0., 2e-12);
        assert!(
            dot(&t, linear.nodal_outflow_w()) > 0.,
            "positive thermal dissipation"
        );
        let physical = linear.physical_contacts();
        assert!(physical.states().is_empty());
        assert_eq!(physical.quadrature_states().len(), 3);
        assert!(!physical.frozen_interfaces().surface_is_nonmatching("bond"));
        let flux = physical.frozen_interfaces().fluxes(&t).unwrap();
        near(
            flux[0].heat_rate_a_to_b_w,
            expected[..4].iter().sum(),
            2e-12,
        );
        assert!(
            flux[0].heat_rate_a_to_b_w < -0.1,
            "a face mean would miss net transfer"
        );
        for state in physical.quadrature_states() {
            near(
                state.resistance.value_m2_k_per_w(),
                resistance(state.side_a_temperature_k, state.side_b_temperature_k),
                2e-14,
            );
            assert_eq!(state.resistance.card_identity(), card.content_hash());
        }
        let old = TemperatureContacts::new(
            cx,
            &mesh,
            &boundary,
            vec![],
            vec![law(&mesh, &card, false, false)],
            1,
        )
        .unwrap();
        let legacy = old.evaluate(cx, &boundary, &t).unwrap();
        assert_eq!(legacy.states().len(), 1);
        assert!(legacy.quadrature_states().is_empty());
        near(
            legacy.frozen_interfaces().fluxes(&t).unwrap()[0].heat_rate_a_to_b_w,
            0.,
            1e-12,
        );
    });
}

#[test]
fn affine_trace_fixture_has_exact_constant_physical_heat() {
    context(|_, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let t: Vec<_> = mesh
            .positions()
            .iter()
            .enumerate()
            .map(|(v, p)| {
                if v < 4 {
                    344. + 19. * p[0]
                } else {
                    312. + 15. * p[0]
                }
            })
            .collect();
        // R=.5+x/16, jump=32+4x, so q=64 W/m2 everywhere.
        let linear = set.linearize_temperature(cx, &boundary, &t, 1).unwrap();
        for (&a, b) in linear.nodal_outflow_w().iter().zip([
            64. / 6.,
            64. / 6.,
            64. / 6.,
            0.,
            -64. / 6.,
            -64. / 6.,
            -64. / 6.,
            0.,
        ]) {
            near(a, b, 2e-12);
        }
        near(
            linear
                .physical_contacts()
                .frozen_interfaces()
                .fluxes(&t)
                .unwrap()[0]
                .heat_rate_a_to_b_w,
            32.,
            2e-12,
        );
    });
}

#[test]
fn complete_tangent_transpose_controls_and_side_reversal_share_the_quadrature() {
    context(|_, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let t = [371., 289., 313., 320., 309., 335., 321., 320.];
        let direction = [0.2, -0.3, 0.7, 0.5, -0.4, 0.8, 0.1, -0.2];
        let dual = [0.8, -0.2, 0.4, 0.1, -0.3, 0.9, 0.7, -0.5];
        let linear = set.linearize_temperature(cx, &boundary, &t, 1).unwrap();
        let tangent = linear.apply(cx, &direction).unwrap();
        let transpose = linear.apply_transpose(cx, &dual).unwrap();
        near(dot(&dual, &tangent), dot(&direction, &transpose), 2e-12);
        let step = 0.001;
        let plus: Vec<_> = t.iter().zip(direction).map(|(t, d)| t + step * d).collect();
        let minus: Vec<_> = t.iter().zip(direction).map(|(t, d)| t - step * d).collect();
        let a = independent_contact_load(&plus);
        let b = independent_contact_load(&minus);
        for i in 0..8 {
            near(tangent[i], (a[i] - b[i]) / (2. * step), 3e-9);
        }
        let frozen = linear.physical_contacts().frozen_interfaces();
        let gradient = frozen
            .matching_resistance_scale_pullback(cx, &t, &dual, 1)
            .unwrap();
        near(
            gradient[0].derivative,
            dot(&dual, linear.nodal_outflow_w()),
            2e-12,
        );
        assert!(gradient[0].mapped);
        assert_eq!(gradient[0].face_pairs, 1);
        assert!(
            frozen
                .nonmatching_log_resistance_pullback(cx, "bond", &t, &dual)
                .unwrap()
                .is_none()
        );
        let reverse_card = source_card(true);
        let reversed = contacts(cx, &mesh, &boundary, &reverse_card, true)
            .linearize_temperature(cx, &boundary, &t, 1)
            .unwrap();
        for (&a, &b) in tangent
            .iter()
            .zip(reversed.apply(cx, &direction).unwrap().iter())
        {
            near(a, b, 2e-12);
        }
        let heat = frozen.fluxes(&t).unwrap()[0].heat_rate_a_to_b_w;
        near(
            reversed
                .physical_contacts()
                .frozen_interfaces()
                .fluxes(&t)
                .unwrap()[0]
                .heat_rate_a_to_b_w,
            -heat,
            2e-12,
        );
    });
}

#[test]
fn paired_vertex_support_and_its_derivative_domain_cannot_hide_behind_samples() {
    context(|gate, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let bad = [400., 280., 280., 320., 320., 320., 320., 320.];
        assert!(
            set.evaluate(cx, &boundary, &bad).is_err(),
            "all interior samples are supported but the hot corner is not"
        );
        let edge = [384., 280., 280., 320., 320., 320., 320., 320.];
        assert!(set.evaluate(cx, &boundary, &edge).is_ok());
        assert!(
            set.linearize_temperature(cx, &boundary, &edge, 1).is_err(),
            "corner lacks two-sided source support"
        );
        gate.request();
        assert!(set.evaluate(cx, &boundary, &[320.; 8]).is_err());
        assert!(
            set.linearize_temperature(cx, &boundary, &[320.; 8], 1)
                .is_err()
        );
    });
}

#[test]
fn bent_matching_surfaces_need_no_planar_clipping_and_count_each_face_once() {
    context(|_, cx| {
        let mesh = mesh(true);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let bind = |faces, points| {
            TemperatureContacts::with_matching_quadrature(
                cx,
                &mesh,
                &boundary,
                vec![],
                vec![law(&mesh, &card, false, true)],
                faces,
                points,
            )
        };
        assert!(bind(1, 6).is_err());
        assert!(bind(2, 5).is_err());
        assert!(
            TemperatureContacts::new(
                cx,
                &mesh,
                &boundary,
                vec![],
                vec![law(&mesh, &card, false, true)],
                2
            )
            .is_err()
        );
        let set = bind(2, 6).unwrap();
        let linear = set
            .linearize_temperature(cx, &boundary, &[330.; 16], 2)
            .unwrap();
        assert_eq!(linear.physical_contacts().quadrature_states().len(), 6);
        assert!(
            !linear
                .physical_contacts()
                .frozen_interfaces()
                .surface_is_nonmatching("bond")
        );
        near(
            linear
                .physical_contacts()
                .quadrature_states()
                .iter()
                .map(|s| s.weight_m2)
                .sum(),
            1.,
            2e-14,
        );
        assert!(
            set.linearize_temperature(cx, &boundary, &[330.; 16], 1)
                .is_err()
        );
    });
}
