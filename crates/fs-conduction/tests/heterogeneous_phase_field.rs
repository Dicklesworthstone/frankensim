//! Original phase sources retain material ownership through spatial transport.
//! All constitutive values are manufactured fixtures, not measured material data.
mod support;

use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError,
};
use fs_conduction::transient::enthalpy::material::phase::MaterialPhaseTransport;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyError, EnthalpyStepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionMesh, ConductivityModel, InterfaceResistance, InterfaceSurface, LinearConfig,
    ScalarField, ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const DT: f64 = 0.01;
const DENSITY: [f64; 2] = [2.0, 6.0];
const OLD: [f64; 8] = [-20.0, 80.0, 20.0, 140.0, 980.0, 1080.0, 1020.0, 1140.0];
const VALUES: [[f64; 3]; 2] = [[2.0, 8.0, 3.0], [4.0, 1.0, 7.0]];

fn state_id(material: usize) -> MaterialStateId {
    MaterialStateId {
        chemistry: format!("manufactured phase material {material}"),
        phase: "declared equilibrium chart".into(),
        process: "numerical fixture".into(),
        revision: 0,
    }
}

fn cards() -> [MaterialCard; 2] {
    std::array::from_fn(|material| {
        let values = VALUES[material];
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("k", KD),
                value: PropertyValue::Curve {
                    abscissa: "f".into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, values[0]), (0.5, values[1]), (1.0, values[2])],
                    dims: KD,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        250.0,
                        1000.0,
                    )
                    .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured heterogeneous phase-source fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        MaterialCard::assemble(state_id(material), claims, vec![]).unwrap()
    })
}

fn charts(cards: &[MaterialCard; 2]) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    std::array::from_fn(|material| {
        EquilibriumEnthalpyPhaseCurve::try_new(
            cards[material].content_hash(),
            [
                (-100.0, 300.0, 0.0),
                (0.0, 400.0, 0.0),
                (100.0, 400.0, 1.0),
                (500.0, 800.0, 1.0),
            ]
            .into_iter()
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h + 1000.0 * material as f64,
                temperature_k: t + 50.0 * material as f64,
                liquid_mass_fraction: f,
                // Current bulk density must not overwrite reference storage.
                bulk_density_kg_m3: 1000.0,
            })
            .collect(),
        )
        .unwrap()
    })
}

fn source<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

fn sources<'a>(
    cards: &[MaterialCard; 2],
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
) -> [PreparedPhaseConductivity<'a>; 2] {
    [
        source(&cards[0], &curves[0], [-100.0, 500.0]),
        source(&cards[1], &curves[1], [900.0, 1500.0]),
    ]
}

fn materials<'a>(
    cx: &fs_exec::Cx<'_>,
    sources: &'a [PreparedPhaseConductivity<'a>; 2],
    density: [f64; 2],
) -> [MaterialPhaseTransport<'a>; 2] {
    std::array::from_fn(|id| {
        MaterialPhaseTransport::try_new_isotropic(cx, &sources[id], density[id]).unwrap()
    })
}

fn mesh(contact: bool) -> ConductionMesh {
    let points = if contact {
        vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ]
    } else {
        vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [2., 0., 0.],
            [3., 0., 0.],
            [2., 1., 0.],
            [2., 0., 1.],
        ]
    };
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        points,
    )
    .unwrap()
}

fn single_mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    )
    .unwrap()
}

fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 8,
        max_elements: 2,
    }
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 16,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 40,
        energy_tolerance_j: 1e-10,
    }
}

fn contact(mesh: &ConductionMesh, boundary: &ThermalBoundary) -> ThermalInterfaces {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: 0.5,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite interface resistance".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |id| SurfaceSpec {
        material: state_id(id),
        texture_frame: "fixed planar reference face".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface(0),
        surface(1),
        SystemContext {
            medium: "declared interface".into(),
            third_body: None,
            environment: "numerical fixture".into(),
            history: "unchanged".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let resistance = InterfaceResistance::from_card(
        "joint",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    let mut pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    assert_eq!(pairs.len(), 1);
    for pair in &mut pairs {
        if mesh.boundary()[pair.side_a].element != 0 {
            std::mem::swap(&mut pair.side_a, &mut pair.side_b);
        }
    }
    ThermalInterfaces::new(
        mesh,
        boundary,
        vec![InterfaceSurface::new("joint", pairs, resistance).unwrap()],
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn temperature(h: f64, material: usize) -> f64 {
    let h = h - 1000.0 * material as f64;
    let t = if h < 0.0 {
        400.0 + h
    } else if h <= 100.0 {
        400.0
    } else {
        300.0 + h
    };
    t + 50.0 * material as f64
}

// Independent truncated-power tetrahedral integral, for distinct nodal h.
// It deliberately does not use the production B-spline/Gauss moment rule.
fn average_k(h: &[f64], material: usize) -> f64 {
    let h: Vec<_> = h.iter().map(|x| x - 1000.0 * material as f64).collect();
    let hinge = |cut: f64| {
        (0..4)
            .map(|i| {
                let denominator = (0..4)
                    .filter(|&j| i != j)
                    .map(|j| h[i] - h[j])
                    .product::<f64>();
                assert_ne!(denominator, 0.0);
                (h[i] - cut).max(0.0).powi(4) / (4.0 * denominator)
            })
            .sum::<f64>()
    };
    let values = VALUES[material];
    let left = (values[1] - values[0]) / 50.0;
    let right = (values[2] - values[1]) / 50.0;
    values[0] + left * hinge(0.0) + (right - left) * hinge(50.0) - right * hinge(100.0)
}

fn contact_residual(h: &[f64], old: &[f64]) -> [f64; 8] {
    let mut result = [0.0; 8];
    let t: Vec<_> = h
        .iter()
        .enumerate()
        .map(|(i, &x)| temperature(x, i / 4))
        .collect();
    for material in 0..2 {
        let offset = 4 * material;
        let k = average_k(&h[offset..offset + 4], material);
        let laplace = [
            3.0 * t[offset] - t[offset + 1] - t[offset + 2] - t[offset + 3],
            t[offset + 1] - t[offset],
            t[offset + 2] - t[offset],
            t[offset + 3] - t[offset],
        ];
        for i in 0..4 {
            result[offset + i] = DENSITY[material] / 24.0 * (h[offset + i] - old[offset + i])
                + DT * k * laplace[i] / 6.0;
        }
    }
    // The matching triangle has area 1/2 and R''=1/2. Its consistent
    // transfer matrix is (1+delta_ab)/12, with opposite signs on the traces.
    for a in 0..3 {
        let outward = (0..3)
            .map(|b| {
                let coefficient = if a == b { 1.0 / 6.0 } else { 1.0 / 12.0 };
                coefficient * (t[b] - t[4 + b])
            })
            .sum::<f64>();
        result[a] += DT * outward;
        result[4 + a] -= DT * outward;
    }
    result
}

#[test]
fn source_assignment_preserves_datums_density_and_independent_uniform_owners() {
    let mesh = mesh(false);
    let single = single_mesh();
    let boundary = boundary(&mesh);
    let single_boundary = self::boundary(&single);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let ids = [1, 0];
    let old = [980., 1080., 1020., 1140., -20., 80., 20., 140.];
    let q = ScalarField::Nodal(vec![6., -2., 9., 4., 8., -3., 12., 5.]);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let storage = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &ids,
            budget(),
        )
        .unwrap();
        assert_eq!(storage.element_material_ids(), ids);
        assert_eq!(storage.vertex_material_ids(), [1, 1, 1, 1, 0, 0, 0, 0]);
        for material in 0..2 {
            assert!(std::ptr::eq(
                storage
                    .phase_transport_material(material)
                    .unwrap()
                    .conductivity_source(),
                &sources[material],
            ));
            assert!(std::ptr::eq(
                storage.reference_material(material).unwrap().curve,
                &curves[material],
            ));
        }
        assert!(storage.phase_transport_material(2).is_none());
        for vertex in 0..8 {
            let id = ids[vertex / 4];
            assert!(std::ptr::eq(
                storage.phase_curve_at(vertex).unwrap(),
                &curves[id]
            ));
            near(
                storage.reference_nodal_masses_kg()[vertex],
                DENSITY[id] / 24.0,
                1e-15,
            );
        }
        let problem = storage.phase_transport_problem(&boundary, &q).unwrap();
        let got = storage
            .advance(cx, problem, None, &old, DT, config())
            .unwrap();
        let flux = storage.phase_transport_heat_flux(cx, &old).unwrap();
        for element in 0..2 {
            let id = ids[element];
            let h = &old[4 * element..4 * element + 4];
            let q_values = if element == 0 {
                vec![6., -2., 9., 4.]
            } else {
                vec![8., -3., 12., 5.]
            };
            let single_source = ScalarField::Nodal(q_values);
            let original = materials[id]
                .prepare_uniform(
                    cx,
                    &single,
                    EnthalpyBudget {
                        max_vertices: 4,
                        max_elements: 1,
                    },
                )
                .unwrap();
            let expected = original
                .advance(
                    cx,
                    materials[id].problem(&single, &single_boundary, &single_source),
                    None,
                    h,
                    DT,
                    config(),
                )
                .unwrap();
            for i in 0..4 {
                near(
                    got.specific_enthalpy_j_kg[4 * element + i],
                    expected.specific_enthalpy_j_kg[i],
                    2e-8,
                );
            }
            for axis in 0..3 {
                let gradient = temperature(h[axis + 1], id) - temperature(h[0], id);
                near(flux[element][axis], -average_k(h, id) * gradient, 2e-11);
            }
        }
        near(got.stored_energy_change_j, DT * 39.0 / 24.0, 1e-10);
        assert!(storage.phase_transport_heat_flux(cx, &old[..7]).is_err());
    });
}

#[test]
fn separate_material_traces_exchange_heat_with_consistent_latent_tangents_and_adjoints() {
    let mesh = mesh(true);
    let boundary = boundary(&mesh);
    let interfaces = contact(&mesh, &boundary);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let zero = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let storage = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            budget(),
        )
        .unwrap();
        let problem = storage.phase_transport_problem(&boundary, &zero).unwrap();
        // Coincident material traces cannot silently become welded or adiabatic.
        assert!(
            storage
                .advance(cx, problem, None, &OLD, DT, config())
                .is_err()
        );
        let response = storage
            .linearize_step(cx, problem, Some(&interfaces), &OLD, DT, config())
            .unwrap();
        let h = &response.primal().specific_enthalpy_j_kg;
        for value in contact_residual(h, &OLD) {
            near(value, 0.0, 2e-10);
        }
        let first_energy = (0..4)
            .map(|i| DENSITY[0] / 24.0 * (h[i] - OLD[i]))
            .sum::<f64>();
        assert!(first_energy > 0.0);
        assert_eq!(response.primal().contact_fluxes.len(), 1);
        let transfer = response.primal().contact_fluxes[0].heat_rate_a_to_b_w;
        assert!(transfer < 0.0);
        near(first_energy, -DT * transfer, 1e-10);
        near(response.primal().stored_energy_change_j, 0.0, 1e-10);
        near(response.primal().energy_residual_j, 0.0, 1e-10);
        for column in 0..8 {
            let mut direction = [0.0; 8];
            direction[column] = 1.0;
            let action = response.transport().apply_jacobian(cx, &direction).unwrap();
            let mut plus = h.clone();
            plus[column] += 0.01;
            let mut minus = h.clone();
            minus[column] -= 0.01;
            let rp = storage
                .evaluate_endpoint(
                    cx,
                    problem,
                    Some(&interfaces),
                    &OLD,
                    &plus,
                    DT,
                    config(),
                    None,
                )
                .unwrap();
            let rm = storage
                .evaluate_endpoint(
                    cx,
                    problem,
                    Some(&interfaces),
                    &OLD,
                    &minus,
                    DT,
                    config(),
                    None,
                )
                .unwrap();
            let ip = contact_residual(&plus, &OLD);
            let im = contact_residual(&minus, &OLD);
            for row in 0..8 {
                near(
                    action[row],
                    (rp.residual_j[row] - rm.residual_j[row]) / 0.02,
                    2e-8,
                );
                near(action[row], (ip[row] - im[row]) / 0.02, 2e-8);
            }
            if [1, 2, 5, 6].contains(&column) {
                assert!(action[4 * (column / 4)].abs() > 1e-7);
            }
        }
        let weights = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, -0.2, 0.4];
        let seed = response
            .transport()
            .temperature_pullback(cx, &weights)
            .unwrap();
        for i in [1, 2, 5, 6] {
            assert_eq!(seed[i], 0.0);
        }
        let gradient = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 128,
                    restart: 8,
                },
            )
            .unwrap();
        let objective = |temperature: &[f64]| {
            temperature
                .iter()
                .zip(weights)
                .map(|(t, w)| (t - 400.0) * w)
                .sum::<f64>()
        };
        for column in [0, 1, 4, 5] {
            let perturbed = |delta| {
                let mut old = OLD;
                old[column] += delta;
                objective(
                    &storage
                        .advance(cx, problem, Some(&interfaces), &old, DT, config())
                        .unwrap()
                        .temperature,
                )
            };
            let e = 0.02;
            let fd = (-perturbed(2. * e) + 8. * perturbed(e) - 8. * perturbed(-e)
                + perturbed(-2. * e))
                / (12. * e);
            near(
                gradient.transport.previous_specific_enthalpy[column],
                fd,
                2e-8,
            );
        }
        for material in 0..2 {
            let perturbed = |delta| {
                let mut density = DENSITY;
                density[material] += delta;
                let changed = self::materials(cx, &sources, density);
                let solver = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                    cx,
                    &mesh,
                    &changed,
                    &[0, 1],
                    budget(),
                )
                .unwrap();
                objective(
                    &solver
                        .advance(
                            cx,
                            solver.phase_transport_problem(&boundary, &zero).unwrap(),
                            Some(&interfaces),
                            &OLD,
                            DT,
                            config(),
                        )
                        .unwrap()
                        .temperature,
                )
            };
            let e = 0.001;
            let fd = (-perturbed(2. * e) + 8. * perturbed(e) - 8. * perturbed(-e)
                + perturbed(-2. * e))
                / (12. * e);
            near(gradient.material_reference_density[material], fd, 2e-7);
        }
    });
}

#[test]
fn native_adaptive_heating_preserves_each_materials_single_latent_account() {
    let mesh = mesh(false);
    let boundary = boundary(&mesh);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let q = ScalarField::Uniform(100.0);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let storage = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            budget(),
        )
        .unwrap();
        let masses = storage.reference_nodal_masses_kg().to_vec();
        let old = [-20., -20., -20., -20., 980., 980., 980., 980.];
        let mut state = storage
            .adaptive_state(
                cx,
                &old,
                0.0,
                0.2,
                AdaptiveEnthalpyConfig {
                    enthalpy_absolute_tolerance_j_kg: 1e-5,
                    temperature_absolute_tolerance_k: 1e-5,
                    liquid_fraction_absolute_tolerance: 1e-7,
                    relative_tolerance: 0.0,
                    minimum_step_s: 1e-6,
                    maximum_step_s: 0.4,
                    endpoint: config(),
                },
            )
            .unwrap();
        let report = state
            .advance_to(
                cx,
                storage.phase_transport_problem(&boundary, &q).unwrap(),
                None,
                3.0,
                100,
            )
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        for i in 0..8 {
            let (h, t, f) = if i < 4 {
                (130.0, 430.0, 1.0)
            } else {
                (1030.0, 450.0, 0.3)
            };
            near(state.specific_enthalpy_j_kg()[i], h, 1e-8);
            near(state.temperature()[i], t, 1e-8);
            near(state.liquid_mass_fraction()[i], f, 1e-10);
        }
        let energy = state
            .specific_enthalpy_j_kg()
            .iter()
            .zip(old)
            .zip(&masses)
            .map(|((&new, old), mass)| mass * (new - old))
            .sum::<f64>();
        near(energy, 100.0, 1e-9);
        assert_eq!(storage.reference_nodal_masses_kg(), masses.as_slice());
    });
}

#[test]
fn material_specific_domain_corner_and_assignment_guards_remain_explicit() {
    let mesh = mesh(false);
    let boundary = boundary(&mesh);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let zero = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        assert!(matches!(
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                cx,
                &mesh,
                &materials,
                &[0, 2],
                budget(),
            ),
            Err(HeterogeneousEnthalpyError::UnknownMaterial {
                element: 1,
                material: 2,
                ..
            })
        ));
        assert!(matches!(
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                cx,
                &mesh,
                &materials,
                &[0, 1],
                EnthalpyBudget {
                    max_vertices: 8,
                    max_elements: 1
                },
            ),
            Err(HeterogeneousEnthalpyError::Enthalpy(
                EnthalpyError::Budget { .. }
            ))
        ));
        let shared = ConductionMesh::new(
            TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]),
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [0., 0., -1.],
            ],
        )
        .unwrap();
        assert!(matches!(
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                cx,
                &shared,
                &materials,
                &[0, 1],
                budget(),
            ),
            Err(HeterogeneousEnthalpyError::SharedVertex { .. })
        ));
        let storage = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &[1, 0],
            budget(),
        )
        .unwrap();
        let problem = storage.phase_transport_problem(&boundary, &zero).unwrap();
        let kink = [1050., 1050., 1050., 1050., 20., 20., 20., 20.];
        assert!(matches!(
            storage.linearize_step(cx, problem, None, &kink, DT, config()),
            Err(EnthalpyAdjointError::MaterialKink)
        ));
        let unrelated = ConductivityModel::isotropic_declared(1.0).unwrap();
        let mut replaced = problem;
        replaced.material = &unrelated;
        assert!(matches!(
            storage.advance(cx, replaced, None, &kink, DT, config()),
            Err(EnthalpyError::InvalidInput(_))
        ));
        with_cancelled_cx(|cancelled| {
            assert!(
                HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                    cancelled,
                    &mesh,
                    &materials,
                    &[0, 1],
                    budget(),
                )
                .is_err()
            );
            assert!(storage.phase_transport_heat_flux(cancelled, &kink).is_err());
        });
        assert!(matches!(
            storage.with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(1.0, 2.0).unwrap(); 2]
            ),
            Err(EnthalpyError::InvalidInput(_))
        ));
        let clipped = [
            source(&cards[0], &curves[0], [-100.0, 500.0]),
            source(&cards[1], &curves[1], [1010.0, 1090.0]),
        ];
        let clipped_materials = self::materials(cx, &clipped, DENSITY);
        let clipped_storage = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &clipped_materials,
            &[0, 1],
            budget(),
        )
        .unwrap();
        let valid = [20., 20., 20., 20., 1050., 1050., 1050., 1050.];
        let mut bad = valid;
        bad[5] = 1005.0;
        let problem = clipped_storage
            .phase_transport_problem(&boundary, &zero)
            .unwrap();
        assert!(matches!(
            clipped_storage.evaluate_endpoint(cx, problem, None, &valid, &bad, DT, config(), None,),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 5, .. })
        ));
        assert!(matches!(
            clipped_storage.advance(cx, problem, None, &bad, DT, config()),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 5, .. })
        ));
        assert_eq!(valid[5], 1050.0);
    });
}
