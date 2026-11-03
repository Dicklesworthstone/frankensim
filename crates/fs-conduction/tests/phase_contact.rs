//! Original two-material phase contact laws on the native spatial h owner.
//! Manufactured sources and analytical balances, not measured interface data.
#[path = "contact_material_state/common.rs"]
mod common;
#[path = "phase_contact/history.rs"]
mod history;
mod support;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::NonmatchingOptions;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_conduction::interface::temperature::{
    MaterialContactState, TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::*;
use fs_exec::Cx;
use fs_matdb::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::with_cx;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b} (tolerance {tolerance})"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn chart(card: &MaterialCard, melt: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, melt - 100.0, 0.0),
            (0.0, melt, 0.0),
            (100.0, melt, 1.0),
            (300.0, melt + 200.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap()
}
fn point() -> QueryPoint {
    common::point()
        .with("aT", 400.0)
        .unwrap()
        .with("bT", 500.0)
        .unwrap()
        .with("T", 450.0)
        .unwrap()
        .with("fa", 0.25)
        .unwrap()
        .with("fb", 0.6)
        .unwrap()
}
fn source(
    cards: &[MaterialCard; 2],
    scale: f64,
    mixed: bool,
    reverse: bool,
    a_limit: f64,
) -> InterfaceSystemCard {
    let mut claim = common::claim("fa", Dims::NONE, scale);
    let (axes, axis_dims, points) = if mixed {
        (
            ["aT".into(), "fb".into()],
            [Temperature::DIMS, Dims::NONE],
            [[300.0, 0.0], [800.0, 0.0], [300.0, 1.0], [800.0, 1.0]],
        )
    } else {
        (
            ["fa".into(), "fb".into()],
            [Dims::NONE; 2],
            [[0.0, 0.0], [a_limit, 0.0], [0.0, 1.0], [a_limit, 1.0]],
        )
    };
    claim.value = PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples: points
            .into_iter()
            .map(|p| {
                (
                    p,
                    scale
                        * if mixed {
                            1.0 + 0.01 * (p[0] - 350.0) + p[1]
                        } else {
                            1.0 + 2.0 * p[0] + p[1]
                        },
                )
            })
            .collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    };
    claim.validity =
        claim
            .validity
            .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.0, a_limit);
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    let order = if reverse { [1, 0] } else { [0, 1] };
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[order[0]].id().clone(),
            texture_frame: "original A".into(),
        },
        SurfaceSpec {
            material: cards[order[1]].id().clone(),
            texture_frame: "original B".into(),
        },
        SystemContext {
            medium: "declared contact".into(),
            third_body: None,
            environment: "test".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn tetrahedra() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ],
    )
    .unwrap()
}
fn slabs(right: [usize; 3]) -> (ConductionMesh, usize) {
    let (left, mut positions) = box_grid([1, 1, 1], [1.0; 3]);
    let split = positions.len();
    let (right, points) = box_grid(right, [1.0; 3]);
    let mut tets = left.tets;
    tets.extend(right.tets.into_iter().map(|t| t.map(|v| v + split as u32)));
    positions.extend(points.into_iter().map(|p| [p[0] + 1.0, p[1], p[2]]));
    (
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap(),
        split,
    )
}
fn sides(mesh: &ConductionMesh, split: usize, axis: usize, coordinate: f64) -> [Vec<usize>; 2] {
    let mut result = [vec![], vec![]];
    for (slot, face) in mesh.boundary().iter().enumerate() {
        if on_box_face(face.centroid[axis], coordinate) {
            let side = usize::from(face.vertices[0] as usize >= split);
            assert!(
                face.vertices
                    .iter()
                    .all(|&v| usize::from(v as usize >= split) == side)
            );
            result[side].push(slot);
        }
    }
    result
}
fn engine<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    split: usize,
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    let assignment: Vec<_> = mesh
        .complex()
        .tets
        .iter()
        .map(|tet| usize::from(tet[0] as usize >= split))
        .collect();
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: &curves[0],
                reference_density_kg_m3: 1.0,
            },
            ReferenceEnthalpyMaterial {
                curve: &curves[1],
                reference_density_kg_m3: 2.0,
            },
        ],
        &assignment,
        EnthalpyBudget {
            max_vertices: mesh.vertex_count(),
            max_elements: mesh.element_count(),
        },
    )
    .unwrap()
}
#[allow(clippy::too_many_arguments)]
fn contacts<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    boundary: &ThermalBoundary,
    cards: &'c [MaterialCard; 2],
    curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    source: &'c InterfaceSystemCard,
    point: &'c QueryPoint,
    sides: [Vec<usize>; 2],
    mode: usize,
    reverse: bool,
) -> TemperatureContacts<'m, 'c> {
    let [mut a, mut b] = sides;
    if reverse {
        std::mem::swap(&mut a, &mut b);
    }
    let axes = if reverse { ["bT", "aT"] } else { ["aT", "bT"] };
    let coordinates = TemperatureContactCoordinates::TwoSided(axes.map(String::from));
    let law = if mode == 2 {
        TemperatureContact::new_nonmatching(
            "joint",
            a,
            b,
            source,
            point,
            coordinates,
            SelectionPolicy::SingleClaimOnly,
            NonmatchingOptions {
                plane_tolerance_m: 1e-12,
                coverage_relative_tolerance: 1e-10,
                max_pair_tests: 100_000,
                max_overlap_triangles: 10_000,
            },
        )
        .unwrap()
    } else {
        let pairs = ThermalInterfaces::coincident_face_pairs(mesh)
            .unwrap()
            .into_iter()
            .map(|p| {
                if a.contains(&p.side_a) {
                    p
                } else {
                    InterfaceFacePair {
                        side_a: p.side_b,
                        side_b: p.side_a,
                    }
                }
            })
            .collect();
        if mode == 0 {
            TemperatureContact::new_two_sided(
                "joint",
                pairs,
                source,
                point,
                axes,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap()
        } else {
            TemperatureContact::new_pointwise(
                "joint",
                pairs,
                source,
                point,
                coordinates,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap()
        }
    };
    let order = if reverse { [1, 0] } else { [0, 1] };
    let law = law
        .with_material_states(MaterialContactState {
            point: ThermalInterfacePoint {
                interface: source,
                materials: [&cards[order[0]], &cards[order[1]]],
                context: point,
                temperature_axes: axes,
                liquid_fraction_axes: if reverse {
                    [Some("fb"), Some("fa")]
                } else {
                    [Some("fa"), Some("fb")]
                },
                mean_temperature_axis: Some("T"),
            },
            curves: [&curves[order[0]], &curves[order[1]]],
        })
        .unwrap();
    TemperatureContacts::with_nonmatching(
        cx,
        mesh,
        boundary,
        vec![],
        vec![],
        vec![law],
        10_000,
        30_000,
    )
    .unwrap()
}
fn config() -> ContactEnthalpyNewtonConfig {
    ContactEnthalpyNewtonConfig {
        step: EnthalpyStepConfig {
            newton: NewtonKrylovConfig {
                relative_tolerance: 0.0,
                absolute_tolerance: 1e-11,
                linear_restart: 16,
                max_linear_cycles: 32,
                forcing_minimum: 1e-12,
                forcing_maximum: 1e-10,
                ..Default::default()
            },
            max_newton_iterations: 32,
            energy_tolerance_j: 1e-9,
        },
        max_residual_evaluations: 256,
        max_face_pairs: 10_000,
        residual_tolerance_j: 1e-9,
    }
}

#[test]
fn local_latent_fraction_laws_retain_both_materials_and_exact_trace_chain_rule() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let engine = engine(cx, &mesh, &curves, 4);
        let source = source(&cards, 1.0, false, false, 1.0);
        let point = point();
        let h = [20., 40., 60., 30., 30., 50., 70., 60.];
        let direction = [1., 1., 1., 0., 0., 0., 0., 0.];
        let dual = [0.3, -0.4, 0.5, 0.2, -0.6, 0.7, -0.8, 0.1];
        for mode in 0..3 {
            let contacts = contacts(
                cx,
                &mesh,
                &boundary,
                &cards,
                &curves,
                &source,
                &point,
                sides(&mesh, 4, 2, 0.),
                mode,
                false,
            );
            assert!(
                contacts
                    .evaluate(
                        cx,
                        &boundary,
                        &[400., 400., 400., 400., 500., 500., 500., 500.]
                    )
                    .is_err()
            );
            let linear = contacts
                .linearize_enthalpy(cx, &boundary, &engine, &h, 10000)
                .unwrap();
            let samples: Vec<(f64, f64, f64)> = if mode == 0 {
                vec![(1. / 3., 1. / 3., 0.5)]
            } else {
                vec![
                    (1. / 6., 1. / 6., 1. / 6.),
                    (2. / 3., 1. / 6., 1. / 6.),
                    (1. / 6., 2. / 3., 1. / 6.),
                ]
            };
            let expected: f64 = samples
                .iter()
                .map(|&(x, y, w)| -100. * w / (1.7 + 0.6 * x + 1.2 * y))
                .sum();
            let expected_derivative: f64 = samples
                .iter()
                .map(|&(x, y, w)| 100. * w * 0.02 / (1.7 + 0.6 * x + 1.2 * y).powi(2))
                .sum();
            near(linear.nodal_outflow_w()[..4].iter().sum(), expected, 2e-11);
            near(linear.nodal_outflow_w().iter().sum(), 0., 2e-11);
            let jx = linear.apply_enthalpy(cx, &direction).unwrap();
            near(jx[..4].iter().sum(), expected_derivative, 2e-12);
            near(jx.iter().sum(), 0., 2e-12);
            let jtz = linear.apply_enthalpy_transpose(cx, &dual).unwrap();
            near(dot(&dual, &jx), dot(&direction, &jtz), 2e-13);
            assert!(linear.apply(cx, &direction).is_err());
            let states = linear.physical_contacts().material_node_states();
            assert_eq!(states.len(), 6);
            for state in states {
                assert_eq!(
                    state.phase.material_card_identity(),
                    cards[usize::from(state.vertex >= 4)].content_hash()
                );
                assert_eq!(state.phase.specific_enthalpy_j_kg(), h[state.vertex]);
                assert_eq!(
                    state.phase.temperature_k(),
                    if state.vertex < 4 { 400. } else { 500. }
                );
            }
            for epsilon in [1e-3, 5e-4] {
                let shifted = |sign: f64| {
                    let hs: Vec<_> = h
                        .iter()
                        .zip(direction)
                        .map(|(h, d)| h + sign * epsilon * d)
                        .collect();
                    let physical = contacts
                        .evaluate_enthalpy(cx, &boundary, &engine, &hs)
                        .unwrap();
                    physical
                        .frozen_interfaces()
                        .fluxes(linear.temperature())
                        .unwrap()[0]
                        .heat_rate_a_to_b_w
                };
                near(
                    (shifted(1.) - shifted(-1.)) / (2. * epsilon),
                    expected_derivative,
                    2e-8,
                );
            }
        }
    });
}

#[test]
fn native_phase_contact_solve_and_adjoint_match_an_independent_latent_quadratic() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let engine = engine(cx, &mesh, &curves, 4);
        let source_card = source(&cards, 1., false, false, 1.);
        let point = point();
        let k = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &k,
            element_materials: None,
            source: &source,
        };
        let old = [25., 25., 25., 30., 60., 60., 60., 65.];
        let dt = 0.05;
        // Each triangle node has m_A=1/24, m_B=1/12 and one third of Q.
        // x(2.1+0.015x)=400dt, while h_B=60-x/2. Bulk T stays 400/500 K.
        let x = 800. * dt / (2.1 + (2.1_f64.powi(2) + 4. * 0.015 * 400. * dt).sqrt());
        let expected_scale = -400. * dt / (2.1 + 0.03 * x);
        for mode in 0..3 {
            let contacts = contacts(
                cx,
                &mesh,
                &boundary,
                &cards,
                &curves,
                &source_card,
                &point,
                sides(&mesh, 4, 2, 0.),
                mode,
                false,
            );
            let accepted = contacts
                .advance_enthalpy_newton(cx, &engine, problem, &old, dt, config())
                .unwrap();
            for v in 0..3 {
                near(accepted.conduction.specific_enthalpy_j_kg[v], 25. + x, 2e-8);
                near(
                    accepted.conduction.specific_enthalpy_j_kg[v + 4],
                    60. - 0.5 * x,
                    2e-8,
                );
            }
            near(accepted.conduction.specific_enthalpy_j_kg[3], old[3], 2e-10);
            near(accepted.conduction.specific_enthalpy_j_kg[7], old[7], 2e-10);
            assert_eq!(
                accepted.conduction.temperature,
                [400., 400., 400., 400., 500., 500., 500., 500.]
            );
            near(accepted.physical.stored_energy_change_j, 0., 1e-9);
            near(accepted.physical.energy_residual_j, 0., 1e-9);
            assert!(accepted.physical.residual_norm_j <= config().residual_tolerance_j);
            let entropy: f64 = accepted
                .conduction
                .specific_enthalpy_j_kg
                .iter()
                .zip(old)
                .zip(engine.reference_nodal_masses_kg())
                .zip(&accepted.conduction.temperature)
                .map(|(((&new, old), &mass), &temperature)| mass * (new - old) / temperature)
                .sum();
            near(entropy, x / 8. * (1. / 400. - 1. / 500.), 2e-11);
            assert!(entropy > 0.);
            let step = contacts
                .linearize_accepted_enthalpy(
                    cx,
                    &engine,
                    problem,
                    &old,
                    dt,
                    config().step,
                    config().residual_tolerance_j,
                    accepted,
                    10000,
                )
                .unwrap();
            let mut seed = [0.; 8];
            seed[0] = 1.;
            assert_eq!(step.temperature_pullback(cx, &seed).unwrap(), vec![0.; 8]);
            let gradient = step
                .pullback_with_resistance_controls(
                    cx,
                    &seed,
                    LinearConfig {
                        tolerance: 1e-11,
                        max_iterations: 256,
                        restart: 8,
                    },
                    10000,
                )
                .unwrap();
            near(
                gradient.resistance_multipliers[0].derivative,
                expected_scale,
                2e-7,
            );
            let direction = [0.3, -0.4, 0.7, 0.2, -0.6, 0.1, 0.5, -0.2];
            let dual = [-0.2, 0.9, 0.1, -0.5, 0.4, -0.3, 0.8, 0.6];
            let jx = step
                .transport()
                .transport()
                .apply_jacobian(cx, &direction)
                .unwrap();
            let jtz = step
                .transport()
                .transport()
                .apply_jacobian_transpose(cx, &dual)
                .unwrap();
            near(dot(&dual, &jx), dot(&direction, &jtz), 2e-12);
            let epsilon = 1e-3;
            let shifted = |sign: f64| {
                let mut h = old;
                h[0] += sign * epsilon;
                contacts
                    .advance_enthalpy_newton(cx, &engine, problem, &h, dt, config())
                    .unwrap()
                    .conduction
                    .specific_enthalpy_j_kg[0]
            };
            near(
                gradient.storage.transport.previous_specific_enthalpy[0],
                (shifted(1.) - shifted(-1.)) / (2. * epsilon),
                3e-7,
            );
        }
    });
}

#[test]
fn mixed_sensible_and_latent_coordinates_keep_temperature_and_fraction_derivatives() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let engine = engine(cx, &mesh, &curves, 4);
        let source = source(&cards, 1., true, false, 1.);
        let point = point();
        let contacts = contacts(
            cx,
            &mesh,
            &boundary,
            &cards,
            &curves,
            &source,
            &point,
            sides(&mesh, 4, 2, 0.),
            1,
            false,
        );
        let h = [-25., -25., -25., -25., 60., 60., 60., 60.];
        let linear = contacts
            .linearize_enthalpy(cx, &boundary, &engine, &h, 10000)
            .unwrap();
        let r = 1.85_f64;
        let jump = -125.;
        near(
            linear.nodal_outflow_w()[..4].iter().sum(),
            0.5 * jump / r,
            2e-11,
        );
        for (direction, expected) in [
            (
                [1., 1., 1., 1., 0., 0., 0., 0.],
                0.5 * (1. / r - jump * 0.01 / r.powi(2)),
            ),
            (
                [0., 0., 0., 0., 1., 1., 1., 1.],
                -0.5 * jump * 0.01 / r.powi(2),
            ),
        ] {
            let actual = linear.apply_enthalpy(cx, &direction).unwrap();
            near(actual[..4].iter().sum(), expected, 2e-12);
        }
    });
}

#[test]
fn independent_opposite_remeshing_and_order_reversal_preserve_phase_contact_power() {
    with_cx(|cx| {
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let point = point();
        for subdivisions in [[1, 1, 1], [1, 2, 3]] {
            let (mesh, split) = slabs(subdivisions);
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let engine = engine(cx, &mesh, &curves, split);
            let h: Vec<_> = (0..mesh.vertex_count())
                .map(|v| if v < split { 25. } else { 60. })
                .collect();
            let direction: Vec<_> = (0..mesh.vertex_count())
                .map(|v| if v < split { 1. } else { 0. })
                .collect();
            let mut forward: Option<Vec<f64>> = None;
            for reverse in [false, true] {
                let card = source(&cards, 1., false, reverse, 1.);
                let contacts = contacts(
                    cx,
                    &mesh,
                    &boundary,
                    &cards,
                    &curves,
                    &card,
                    &point,
                    sides(&mesh, split, 0, 1.),
                    2,
                    reverse,
                );
                let linear = contacts
                    .linearize_enthalpy(cx, &boundary, &engine, &h, 10000)
                    .unwrap();
                near(
                    linear.nodal_outflow_w()[..split].iter().sum(),
                    -100. / 2.1,
                    3e-10,
                );
                near(linear.nodal_outflow_w().iter().sum(), 0., 3e-10);
                let derivative = linear.apply_enthalpy(cx, &direction).unwrap();
                near(
                    derivative[..split].iter().sum(),
                    2. / 2.1_f64.powi(2),
                    3e-11,
                );
                let flux = linear
                    .physical_contacts()
                    .frozen_interfaces()
                    .fluxes(linear.temperature())
                    .unwrap();
                near(
                    flux[0].heat_rate_a_to_b_w,
                    if reverse { 100. / 2.1 } else { -100. / 2.1 },
                    3e-10,
                );
                if let Some(expected) = &forward {
                    for (&a, &b) in linear.nodal_outflow_w().iter().zip(expected) {
                        near(a, b, 3e-10);
                    }
                } else {
                    forward = Some(linear.nodal_outflow_w().to_vec());
                }
            }
        }
    });
}

#[test]
fn original_chart_identity_source_support_and_cancellation_refuse_without_mutation() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let engine = engine(cx, &mesh, &curves, 4);
        let card = source(&cards, 1., false, false, 0.5);
        let point = point();
        let contacts = contacts(
            cx,
            &mesh,
            &boundary,
            &cards,
            &curves,
            &card,
            &point,
            sides(&mesh, 4, 2, 0.),
            1,
            false,
        );
        let h = [25., 25., 25., 30., 60., 60., 60., 65.];
        let before = h;
        let mut outside = h;
        outside[..3].fill(60.);
        assert!(
            contacts
                .evaluate_enthalpy(cx, &boundary, &engine, &outside)
                .is_err()
        );
        let copies = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let wrong_engine = crate::engine(cx, &mesh, &copies, 4);
        assert!(
            contacts
                .evaluate_enthalpy(cx, &boundary, &wrong_engine, &h)
                .is_err()
        );
        let mut corner = h;
        corner[..3].fill(0.);
        assert!(
            contacts
                .evaluate_enthalpy(cx, &boundary, &engine, &corner)
                .is_ok()
        );
        assert!(
            contacts
                .linearize_enthalpy(cx, &boundary, &engine, &corner, 10000)
                .is_err()
        );
        support::with_cancelled_cx(|cancelled| {
            assert!(
                contacts
                    .evaluate_enthalpy(cancelled, &boundary, &engine, &h)
                    .is_err()
            )
        });
        assert_eq!(h, before);
        let first = contacts
            .evaluate_enthalpy(cx, &boundary, &engine, &h)
            .unwrap();
        let retry = contacts
            .evaluate_enthalpy(cx, &boundary, &engine, &h)
            .unwrap();
        assert_eq!(first.quadrature_states(), retry.quadrature_states());
        assert_eq!(first.material_node_states(), retry.material_node_states());
    });
}

#[test]
fn original_source_cell_tangent_allows_melting_to_begin_at_zero_fraction() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let engine = engine(cx, &mesh, &curves, 4);
        let card = source(&cards, 1., false, false, 1.);
        let point = point();
        let contacts = contacts(
            cx,
            &mesh,
            &boundary,
            &cards,
            &curves,
            &card,
            &point,
            sides(&mesh, 4, 2, 0.),
            1,
            false,
        );
        let k = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &k,
            element_materials: None,
            source: &source,
        };
        let old = [0., 0., 0., 0., 60., 60., 60., 60.];
        let dt = 0.005;
        assert!(
            contacts
                .linearize_enthalpy(cx, &boundary, &engine, &old, 10000)
                .is_err()
        );
        let accepted = contacts
            .advance_enthalpy_newton(cx, &engine, problem, &old, dt, config())
            .unwrap();
        let x = 800. * dt / (1.6 + (1.6_f64.powi(2) + 4. * 0.015 * 400. * dt).sqrt());
        for v in 0..3 {
            near(accepted.conduction.specific_enthalpy_j_kg[v], x, 2e-8);
            near(
                accepted.conduction.specific_enthalpy_j_kg[v + 4],
                60. - 0.5 * x,
                2e-8,
            );
            assert!(accepted.conduction.liquid_mass_fraction[v] > 0.);
        }
        near(accepted.physical.energy_residual_j, 0., 1e-9);
    });
}

#[test]
fn strict_local_h_actions_refuse_caloric_corners_for_temperature_only_and_fixed_contacts() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.), chart(&cards[1], 500.)];
        let engine = engine(cx, &mesh, &curves, 4);
        let card = source(&cards, 1., false, false, 1.);
        let point = point();
        for variable in [false, true] {
            let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
            let (fixed, changing) = if variable {
                (
                    vec![],
                    vec![
                        TemperatureContact::new_two_sided(
                            "joint",
                            pairs,
                            &card,
                            &point,
                            ["aT", "bT"],
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .unwrap(),
                    ],
                )
            } else {
                (
                    vec![
                        InterfaceSurface::new(
                            "joint",
                            pairs,
                            InterfaceResistance::from_card(
                                "joint",
                                &card,
                                &point,
                                SelectionPolicy::SingleClaimOnly,
                            )
                            .unwrap(),
                        )
                        .unwrap(),
                    ],
                    vec![],
                )
            };
            let contacts =
                TemperatureContacts::new(cx, &mesh, &boundary, fixed, changing, 10000).unwrap();
            let corner = [0., 0., 0., 25., 60., 60., 60., 60.];
            assert!(
                contacts
                    .linearize_enthalpy(cx, &boundary, &engine, &corner, 10000)
                    .is_err()
            );
            // Nonparticipating apex corners do not invalidate the local contact
            // response: every contact trace lies strictly inside its plateau.
            let interior = [25., 25., 25., 0., 60., 60., 60., 0.];
            let linear = contacts
                .linearize_enthalpy(cx, &boundary, &engine, &interior, 10000)
                .unwrap();
            assert_eq!(linear.apply_enthalpy(cx, &[1.; 8]).unwrap(), vec![0.; 8]);
        }
    });
}
