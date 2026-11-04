//! G1/G3 checks of directional phase sources in the existing coupled owners.
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    InterfaceResistance, InterfaceSurface, ThermalBoundary, ThermalInterfaces,
};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};

const K_A: [[f64; 3]; 3] = [[2.0, 8.0, 3.0], [11.0, 4.0, 9.0], [5.0, 3.0, 12.0]];
const K_B: [[f64; 3]; 3] = [[7.0, 2.0, 10.0], [1.5, 6.0, 4.0], [9.0, 11.0, 3.0]];
// Rows are the principal directions in the mesh frame.
const AXES_A: [[f64; 3]; 3] = [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
const AXES_B: [[f64; 3]; 3] = [[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]];
const KEYS: [&str; 3] = ["k0", "k1", "k2"];

fn contact(mesh: &ConductionMesh, boundary: &ThermalBoundary) -> ThermalInterfaces {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: 2.0,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured directional phase contact".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid/liquid".into(),
            process: "manufactured fixture".into(),
            revision: 0,
        },
        texture_frame: "reference".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("first orthotropic material"),
        surface("second orthotropic material"),
        SystemContext {
            medium: "declared finite interface".into(),
            third_body: None,
            environment: "closed coupon".into(),
            history: "fixed".into(),
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
    if mesh.boundary()[pairs[0].side_a].element != 0 {
        let pair = &mut pairs[0];
        std::mem::swap(&mut pair.side_a, &mut pair.side_b);
    }
    ThermalInterfaces::new(
        mesh,
        boundary,
        vec![InterfaceSurface::new("joint", pairs, resistance).unwrap()],
    )
    .unwrap()
}

// Independent exact simplex hinge integral. Shift away each material's enthalpy
// datum before taking powers; the distinct nodal fields below cross phase knots.
fn mean_conductivity(values: [f64; 3], h: &[f64], plateau: [f64; 2]) -> f64 {
    let nodal: [f64; 4] = std::array::from_fn(|i| h[i] - plateau[0]);
    let hinge = |cut: f64| {
        (0..4)
            .map(|i| {
                let denominator = (0..4)
                    .filter(|&j| j != i)
                    .map(|j| nodal[i] - nodal[j])
                    .product::<f64>();
                assert_ne!(denominator, 0.0);
                (nodal[i] - cut).max(0.0).powi(4) / (4.0 * denominator)
            })
            .sum::<f64>()
    };
    let half_latent = (plateau[1] - plateau[0]) / 2.0;
    let low = (values[1] - values[0]) / half_latent;
    let high = (values[2] - values[1]) / half_latent;
    values[0] + low * hinge(0.0) + (high - low) * hinge(half_latent)
        - high * hinge(2.0 * half_latent)
}

#[test]
fn heterogeneous_directional_phase_bodies_exchange_heat_with_one_closed_energy_account() {
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap();
    let cards = [orthotropic_card(K_A), orthotropic_card(K_B)];
    let curves = [
        caloric(&cards[0]),
        EquilibriumEnthalpyPhaseCurve::try_new(
            cards[1].content_hash(),
            [
                (900.0, 350.0, 0.0),
                (1000.0, 450.0, 0.0),
                (1200.0, 450.0, 1.0),
                (1550.0, 800.0, 1.0),
            ]
            .into_iter()
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 700.0,
            })
            .collect(),
        )
        .unwrap(),
    ];
    let ranges = [
        [[-90.0, 450.0], [-80.0, 420.0], [-70.0, 400.0]],
        [[920.0, 1490.0], [930.0, 1470.0], [940.0, 1450.0]],
    ];
    let sources: [[PreparedPhaseConductivity<'_>; 3]; 2] = std::array::from_fn(|m| {
        std::array::from_fn(|axis| {
            prepare_direction(
                &cards[m],
                &curves[m],
                KEYS[axis],
                ranges[m][axis],
                150_000.0,
            )
        })
    });
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let contact = contact(&mesh, &boundary);
    let zero = ScalarField::Uniform(0.0);
    let initial = [-20.0, 80.0, 20.0, 140.0, 960.0, 1080.0, 1160.0, 1330.0];
    with_cx(|cx| {
        let materials = [
            MaterialPhaseTransport::try_new_orthotropic(cx, sources[0].each_ref(), AXES_A, 2.0)
                .unwrap(),
            MaterialPhaseTransport::try_new_orthotropic(cx, sources[1].each_ref(), AXES_B, 3.0)
                .unwrap(),
        ];
        assert_eq!(materials[0].enthalpy_bounds_j_kg(), [-70.0, 400.0]);
        assert_eq!(materials[1].enthalpy_bounds_j_kg(), [940.0, 1450.0]);
        let engine = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let masses = engine.reference_nodal_masses_kg().to_vec();
        for vertex in 0..8 {
            let m = vertex / 4;
            near(masses[vertex], (2.0 + m as f64) / 24.0, 1e-15);
            assert!(std::ptr::eq(
                engine.phase_curve_at(vertex).unwrap(),
                &curves[m]
            ));
        }
        for m in 0..2 {
            let retained = engine.phase_transport_material(m).unwrap();
            assert!(std::ptr::eq(retained, &materials[m]));
            for axis in 0..3 {
                assert!(std::ptr::eq(
                    retained.conductivity_sources()[axis],
                    &sources[m][axis]
                ));
            }
        }
        let problem = engine.phase_transport_problem(&boundary, &zero).unwrap();
        let mut policy = config();
        policy.newton.linear_restart = 8;
        let mut h = initial.to_vec();
        for _ in 0..4 {
            let next = engine
                .advance(cx, problem, Some(&contact), &h, DT, policy)
                .unwrap();
            let energy: [f64; 2] = std::array::from_fn(|m| {
                (4 * m..4 * m + 4)
                    .map(|i| masses[i] * (next.specific_enthalpy_j_kg[i] - h[i]))
                    .sum()
            });
            // The face is area 1/2, with R''=2; a P1 jump integrates to its
            // three-node mean. Internal bulk conduction cannot alter either sum.
            let contact_power = (0..3)
                .map(|i| (next.temperature[i] - next.temperature[i + 4]) / 12.0)
                .sum::<f64>();
            assert!(energy[0] > 0.0 && energy[1] < 0.0);
            near(energy[0], -DT * contact_power, 2e-10);
            near(energy[1], DT * contact_power, 2e-10);
            near(energy.iter().sum(), 0.0, 2e-10);
            near(next.stored_energy_change_j, 0.0, 2e-10);
            near(next.energy_residual_j, 0.0, 2e-10);
            assert_eq!(
                (
                    next.source_w,
                    next.neumann_out_w,
                    next.robin_out_w,
                    next.dirichlet_in_w
                ),
                (0.0, 0.0, 0.0, 0.0)
            );
            near(
                next.contact_fluxes[0].heat_rate_a_to_b_w,
                contact_power,
                1e-11,
            );

            let flux = engine
                .phase_transport_heat_flux(cx, &next.specific_enthalpy_j_kg)
                .unwrap();
            for m in 0..2 {
                let offset = 4 * m;
                let mut gradient: [f64; 3] = std::array::from_fn(|axis| {
                    next.temperature[offset + axis + 1] - next.temperature[offset]
                });
                if m == 0 {
                    gradient[2] = -gradient[2];
                }
                let axes = [AXES_A, AXES_B][m];
                let conductivity = [K_A, K_B][m].map(|values| {
                    mean_conductivity(
                        values,
                        &next.specific_enthalpy_j_kg[offset..offset + 4],
                        [[0.0, 100.0], [1000.0, 1200.0]][m],
                    )
                });
                for row in 0..3 {
                    let expected: f64 = (0..3)
                        .map(|axis| {
                            let along = (0..3).map(|j| axes[axis][j] * gradient[j]).sum::<f64>();
                            -conductivity[axis] * axes[axis][row] * along
                        })
                        .sum();
                    near(flux[m][row], expected, 2e-9);
                }
            }
            h = next.specific_enthalpy_j_kg;
        }
        assert_eq!(engine.reference_nodal_masses_kg(), masses);
        let mut unsupported = initial;
        unsupported[4] = 939.0; // Only the third directional source excludes this h.
        assert!(matches!(
            engine.advance(cx, problem, Some(&contact), &unsupported, DT, policy),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 4, .. })
        ));
    });
}

#[test]
fn accepted_orthotropic_conductivity_scale_adjoint_matches_resolved_source_scaling() {
    let mesh = mesh();
    let card = orthotropic_card(K_A);
    let curve = caloric(&card);
    let sources = KEYS.map(|key| prepare_direction(&card, &curve, key, [-100.0, 500.0], 150_000.0));
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::nodal("manufactured deposition", 4, Q.to_vec()).unwrap();
    let weights = [0.7, -0.4, 0.2, 1.1];
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_orthotropic(cx, sources.each_ref(), AXES_A, 2.0)
            .unwrap();
        let engine = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &heat);
        let accepted = engine
            .advance(cx, problem, None, &OLD, DT, config())
            .unwrap();
        let response = engine
            .linearize_accepted(cx, problem, None, &OLD, DT, config(), accepted)
            .unwrap();
        let seed = response.temperature_pullback(cx, &weights).unwrap();
        let dual = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 64,
                    restart: 4,
                },
            )
            .unwrap();
        let gradient = response
            .conductivity_scale_pullback(cx, &dual.adjoint, 1)
            .unwrap();
        assert_eq!(gradient.len(), 1);
        assert!(gradient[0].abs() > 1e-5);
        // Independently reconstruct the original material claim in all three
        // directions. Caloric values, reference mass, loads and axes stay fixed.
        let objective = |scale: f64| {
            let scaled_card = orthotropic_card(K_A.map(|row| row.map(|value| value * scale)));
            let scaled_curve = caloric(&scaled_card);
            let scaled_sources = KEYS.map(|key| {
                prepare_direction(&scaled_card, &scaled_curve, key, [-100.0, 500.0], 150_000.0)
            });
            let scaled_law = MaterialPhaseTransport::try_new_orthotropic(
                cx,
                scaled_sources.each_ref(),
                AXES_A,
                2.0,
            )
            .unwrap();
            let scaled_engine = scaled_law.prepare_uniform(cx, &mesh, budget()).unwrap();
            let next = scaled_engine
                .advance(
                    cx,
                    scaled_law.problem(&mesh, &boundary, &heat),
                    None,
                    &OLD,
                    DT,
                    config(),
                )
                .unwrap();
            next.temperature
                .iter()
                .zip(weights)
                .map(|(t, weight)| weight * (t - 400.0))
                .sum::<f64>()
        };
        let epsilon = 0.001;
        let finite_difference = (-objective(1.0 + 2.0 * epsilon) + 8.0 * objective(1.0 + epsilon)
            - 8.0 * objective(1.0 - epsilon)
            + objective(1.0 - 2.0 * epsilon))
            / (12.0 * epsilon);
        near(
            gradient[0],
            finite_difference,
            1e-7 * (1.0 + finite_difference.abs()),
        );
    });
}

#[test]
fn phase_history_admits_all_directional_bounds_with_a_nodal_reference_inventory() {
    use fs_conduction::bc::history::ThermalBoundaryHistory;
    use fs_conduction::bc::history::march::HeatSourceHistory;
    use fs_conduction::bc::history::march::enthalpy::{
        EnthalpyHistoryBudget, EnthalpyHistoryRun, PhaseEnthalpyHistoryProblem,
    };
    use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyError;
    use fs_scenario::scenario::Environment;

    let mesh = mesh();
    let card = orthotropic_card(K_A);
    let curve = caloric(&card);
    let ranges = [[-100.0, 500.0], [10.0, 140.0], [20.0, 130.0]];
    let sources: [PreparedPhaseConductivity<'_>; 3] = std::array::from_fn(|axis| {
        prepare_direction(&card, &curve, KEYS[axis], ranges[axis], 150_000.0)
    });
    with_cx(|cx| {
        let materials =
            [
                MaterialPhaseTransport::try_new_orthotropic(cx, sources.each_ref(), AXES_A, 2.0)
                    .unwrap(),
            ];
        let boundary = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            vec![],
            0.0,
        )
        .unwrap();
        let create = |initial| {
            EnthalpyHistoryRun::from_phase_materials_with_nodal_reference_density(
                cx,
                PhaseEnthalpyHistoryProblem {
                    boundary: &boundary,
                    interfaces: None,
                    source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)),
                },
                &materials,
                &[0],
                &[1.0, 2.0, 3.0, 4.0],
                &[initial; 4],
                EnthalpyHistoryBudget {
                    spatial: budget(),
                    max_steps: 1,
                },
            )
        };
        for initial in [15.0, 135.0] {
            assert!(matches!(
                create(initial),
                Err(HeterogeneousEnthalpyError::Enthalpy(
                    EnthalpyError::SourceEnthalpyOutsideDomain {
                        vertex: 0,
                        lower_j_kg: 20.0,
                        upper_j_kg: 130.0,
                        ..
                    }
                ))
            ));
        }
        let run = create(50.0).unwrap();
        assert_eq!(run.specific_enthalpy_j_kg(), &[50.0; 4]);
        assert_eq!(run.accepted_steps(), 0);
        assert!(run.storage().phase_transport_material(0).is_none());
        assert_eq!(
            run.storage().phase_conductivity_enthalpy_bounds_j_kg(0),
            Some([20.0, 130.0])
        );
        let retained = run.storage().phase_conductivity_sources(0).unwrap();
        for axis in 0..3 {
            assert!(std::ptr::eq(retained[axis], &sources[axis]));
        }
    });
}

#[test]
fn assigned_scalar_and_mixed_dispatch_preserve_endpoints_and_adjoints() {
    #[derive(Clone, Copy)]
    enum Interpretation {
        Isotropic,
        EqualPrincipal,
        Directional,
    }
    use Interpretation::{Directional, EqualPrincipal, Isotropic};

    let positions = [-1.0, 1.0]
        .into_iter()
        .flat_map(|z| {
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, z],
            ]
        })
        .collect();
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        positions,
    )
    .unwrap();
    let cards = [orthotropic_card(K_A), orthotropic_card(K_B)];
    let curves = cards.each_ref().map(caloric);
    let sources: [[PreparedPhaseConductivity<'_>; 3]; 2] =
        std::array::from_fn(|material| {
            KEYS.map(|key| {
                prepare_direction(
                    &cards[material],
                    &curves[material],
                    key,
                    [-100.0, 500.0],
                    150_000.0,
                )
            })
        });
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let contact = contact(&mesh, &boundary);
    let zero = ScalarField::Uniform(0.0);
    let old = [-20.0, 80.0, 20.0, 140.0, -10.0, 30.0, 90.0, 180.0];
    let direction = [0.4, -0.7, 0.2, 0.6, -0.5, 0.3, 0.8, -0.1];
    let weights = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, -0.2, 0.8];

    with_cx(|cx| {
        let evaluate = |interpretations: [Interpretation; 2]| {
            let materials: [MaterialPhaseTransport<'_>; 2] =
                std::array::from_fn(|material| {
                    let density = 2.0 + material as f64;
                    let axes = [AXES_A, AXES_B][material];
                    match interpretations[material] {
                        Isotropic => MaterialPhaseTransport::try_new_isotropic(
                            cx,
                            &sources[material][0],
                            density,
                        ),
                        EqualPrincipal => MaterialPhaseTransport::try_new_orthotropic(
                            cx,
                            [&sources[material][0]; 3],
                            axes,
                            density,
                        ),
                        Directional => MaterialPhaseTransport::try_new_orthotropic(
                            cx,
                            sources[material].each_ref(),
                            axes,
                            density,
                        ),
                    }
                    .unwrap()
                });
            let engine = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                cx,
                &mesh,
                &materials,
                &[0, 1],
                EnthalpyBudget {
                    max_vertices: 8,
                    max_elements: 2,
                },
            )
            .unwrap();
            let problem = engine.phase_transport_problem(&boundary, &zero).unwrap();
            let mut policy = config();
            policy.newton.linear_restart = 8;
            let response = engine
                .linearize_step(cx, problem, Some(&contact), &old, DT, policy)
                .unwrap();
            let thermal = response.transport();
            let action = thermal.apply_jacobian(cx, &direction).unwrap();
            let seed = thermal.temperature_pullback(cx, &weights).unwrap();
            let dual = thermal
                .pullback(
                    cx,
                    &seed,
                    LinearConfig {
                        tolerance: 1e-12,
                        max_iterations: 64,
                        restart: 8,
                    },
                )
                .unwrap();
            let scales = thermal
                .conductivity_scale_pullback(cx, &dual.adjoint, 2)
                .unwrap();
            [
                response.primal().specific_enthalpy_j_kg.clone(),
                action,
                dual.previous_specific_enthalpy,
                scales,
            ]
        };

        // Equal principal sources must recover the declared scalar physics.
        let scalar = evaluate([Isotropic, Isotropic]);
        let equal = evaluate([EqualPrincipal, EqualPrincipal]);
        let mixed = evaluate([Isotropic, Directional]);
        let directional = evaluate([EqualPrincipal, Directional]);
        for (actual, expected) in [(&scalar, &equal), (&mixed, &directional)] {
            for (actual, expected) in actual.iter().zip(expected.iter()) {
                assert_eq!(actual.len(), expected.len());
                for (&actual, &expected) in actual.iter().zip(expected.iter()) {
                    near(actual, expected, 3e-9);
                }
            }
        }
        assert!(
            scalar[0]
                .iter()
                .zip(&mixed[0])
                .any(|(&scalar, &mixed)| (scalar - mixed).abs() > 1e-3),
            "the directional second material must change the physical endpoint"
        );
    });
}
