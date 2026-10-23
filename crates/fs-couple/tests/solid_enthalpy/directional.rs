//! Directional source data reaches the same enthalpy, FEM and prestress owners.
use super::*;
use fs_material::state_point::{
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
};

const IDENTITY: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const QUARTER_TURN: [[f64; 3]; 3] = [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]];

// No scalar E, nu or alpha exists. Density is supported ONLY at the original
// temperature: subsequent mechanics must consume the conserved reference mass.
fn directional_material(
    young: [f64; 3],
    alpha: [f64; 3],
    expansion_upper_k: f64,
) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, first, last, upper| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., first), (500., last)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.,
                        upper,
                    )
                    .with("lot", 7., 7.),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic directional solid".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add("density", Density::DIMS, 6., 2., 300.);
    for (name, value) in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.into_iter().zip(young) {
        add(name, Pressure::DIMS, value, 0.6 * value, 500.);
    }
    for name in ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES {
        add(name, Pressure::DIMS, 500., 300., 500.);
    }
    for name in ORTHOTROPIC_POISSON_RATIO_PROPERTIES {
        add(name, Dims::NONE, 0., 0., 500.);
    }
    for (name, value) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .zip(alpha)
    {
        add(
            name,
            INVERSE_TEMPERATURE_DIMS,
            value,
            value,
            expansion_upper_k,
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored directional specimen".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.,
                temperature_k: 300.,
                liquid_mass_fraction: 0.,
                bulk_density_kg_m3: 6.,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 2000.,
                temperature_k: 500.,
                liquid_mass_fraction: 0.,
                bulk_density_kg_m3: 2.,
            },
        ],
    )
    .unwrap();
    (card, curve)
}

fn oriented<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    rotation: [[f64; 3]; 3],
) -> TetCellMaterial<'a> {
    TetCellMaterial {
        stiffness_model: ThermalStiffnessModel::OrientedOrthotropic {
            principal_to_mesh: rotation,
            orientation_identity: ContentHash([8; 32]),
        },
        expansion_model: ThermalExpansionModel::OrientedOrthotropic {
            principal_to_mesh: rotation,
            principal_frame: ContentHash([9; 32]),
        },
        ..source(30, card, curve, context)
    }
}

#[test]
fn rotated_directional_expansion_drives_actual_motion_with_reference_only_density() {
    with_cx(|cx| {
        let (card, curve) = directional_material([1000., 2000., 3000.], [1e-5, 2e-5, -5e-6], 500.);
        let thermal = initial(&curve, &[(30, 1000.)], 1.);
        let context = point(300.);
        let q = std::f64::consts::FRAC_1_SQRT_2;
        let rotation = [[q, -q, 0.], [q, q, 0.], [0., 0., 1.]];
        let model = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![oriented(&card, &curve, &context, rotation)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let loads = model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let next = model
            .advance_to(&thermal, 1., thermal_config(), &loads, &FREE, config(), cx)
            .unwrap();
        let expected = [
            [0., 0., 0.],
            [0.0015, 0., 0.],
            [-0.001, 0.0015, 0.],
            [0., 0., -0.0005],
        ];
        for (actual, expected) in next
            .equilibrium
            .solution
            .displacement
            .displacement_m()
            .iter()
            .zip(expected)
        {
            for (&a, e) in actual.iter().zip(expected) {
                near(a, e, 1e-11);
            }
        }
        for &stress in next.equilibrium.solution.stress.elements()[0].stress_mandel_pa() {
            near(stress, 0., 1e-8);
        }
        let sample = &next.state.samples()[0];
        assert_eq!(sample.elastic.properties().len(), 9);
        assert!(sample.elastic.property("density").is_none());
        assert_eq!(sample.reference_elastic.properties().len(), 10);
        assert_eq!(sample.expansion.reference().properties().len(), 3);
        assert_eq!(sample.expansion.current().properties().len(), 3);
        assert!(sample.expansion.transform().is_some());
        assert_eq!(
            sample
                .expansion
                .thermal_strain()
                .elastic_material_state_identity(),
            next.state.materials()[0].material_state_identity
        );
        near(
            next.state
                .problem(&FREE)
                .assemble(cx)
                .unwrap()
                .total_mass_kg,
            1.,
            1e-12,
        );
        near(next.thermal_audit.external_heat_j, 1000., 1e-8);
        near(next.thermal_audit.stored_energy_change_j, 1000., 1e-8);
        assert_eq!(next.state.mesh_frame(), model.mesh_frame());
        assert_eq!(thermal.time_s(), 0.);
    });
}

#[test]
fn data_and_orientation_changes_reach_restrained_stress_and_the_same_vibration_tangent() {
    with_cx(|cx| {
        for (young, alpha) in [
            ([1000., 2000., 3000.], [1e-5, 2e-5, -5e-6]),
            ([2000., 1500., 3000.], [2e-5, 1e-5, -5e-6]),
        ] {
            let (card, curve) = directional_material(young, alpha, 500.);
            let thermal = initial(&curve, &[(30, 1000.)], 1.);
            let context = point(300.);
            let model = SolidEnthalpy::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![oriented(&card, &curve, &context, QUARTER_TURN)],
                TetAssemblyBudget::standard(),
                0.1,
                config().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = model
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let fixed: Vec<_> = (0..12).collect();
            let clamped = model
                .advance_to(&thermal, 1., thermal_config(), &loads, &fixed, config(), cx)
                .unwrap();
            let expected = [
                -0.8 * young[1] * alpha[1] * 100.,
                -0.8 * young[0] * alpha[0] * 100.,
                -0.8 * young[2] * alpha[2] * 100.,
                0.,
                0.,
                0.,
            ];
            for (&actual, expected) in clamped.equilibrium.solution.stress.elements()[0]
                .stress_mandel_pa()
                .iter()
                .zip(expected)
            {
                near(actual, expected, 1e-9);
            }
            // Fix x while allowing transverse thermal motion. The same accepted
            // source tensor and eigenstrain then feed the existing prestress pencil.
            let fixed = [0, 1, 2, 3, 6, 8, 9];
            let loaded = model
                .advance_to(&thermal, 1., thermal_config(), &loads, &fixed, config(), cx)
                .unwrap();
            let problem = loaded.state.problem(&fixed);
            let tangent = problem
                .assemble_thermal_tangent_with_nodal_forces(
                    &loaded.equilibrium.solution.displacement,
                    loaded.state.thermal_strains(),
                    loads.nodal_forces_n(),
                    config().mesh_frame,
                    0.1,
                    cx,
                )
                .unwrap();
            let transverse = tangent
                .pencil()
                .free_dofs
                .iter()
                .position(|&d| d == 4)
                .unwrap();
            near(
                tangent.geometric_stiffness().get(transverse, transverse),
                expected[0] / 6.,
                1e-9,
            );
            near(tangent.pencil().total_mass_kg, 1., 1e-12);
            assert_eq!(
                tangent.stress().identity(),
                loaded.equilibrium.solution.stress.identity()
            );
        }
    });
}

#[test]
fn expansion_and_stiffness_axes_are_independent_and_the_mesh_frame_cannot_change() {
    with_cx(|cx| {
        let (card, curve) = directional_material([1000., 2000., 3000.], [1e-5, 2e-5, -5e-6], 500.);
        let thermal = initial(&curve, &[(30, 1000.)], 1.);
        let context = point(300.);
        let mut cell = oriented(&card, &curve, &context, IDENTITY);
        cell.expansion_model = ThermalExpansionModel::OrientedOrthotropic {
            principal_to_mesh: QUARTER_TURN,
            principal_frame: ContentHash([9; 32]),
        };
        let model = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![cell],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let loads = model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let next = model
            .advance_to(&thermal, 1., thermal_config(), &loads, &fixed, config(), cx)
            .unwrap();
        let stress = next.equilibrium.solution.stress.elements()[0].stress_mandel_pa();
        near(stress[0], -800. * 0.002, 1e-9);
        near(stress[1], -1600. * 0.001, 1e-9);
        near(stress[2], -2400. * (-0.0005), 1e-9);
        assert!(matches!(
            next.state.solve(
                &loads,
                &fixed,
                SolidEquilibriumConfig {
                    mesh_frame: ContentHash([2; 32]),
                    ..config()
                },
                cx
            ),
            Err(SolidEnthalpyError::Invalid(_))
        ));
        assert!(
            SolidEnthalpy::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![oriented(&card, &curve, &context, IDENTITY)],
                TetAssemblyBudget::standard(),
                0.1,
                ContentHash([0; 32]),
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn intermediate_directional_source_exit_refuses_even_after_recooling() {
    with_cx(|cx| {
        let (card, curve) = directional_material([1000., 2000., 3000.], [1e-5, 2e-5, -5e-6], 350.);
        let seed = initial(&curve, &[(30, 1000.)], 1.);
        let cell = seed.cells()[&CellId(30)].clone();
        let cooling = Boundary {
            source_w: -1000.,
            ..cell.boundary
        };
        let thermal = DepositionSimulation::new(
            vec![
                TimedEvent {
                    time_s: 0.,
                    event: Event::Birth(cell),
                },
                TimedEvent {
                    time_s: 1.,
                    event: Event::Boundary(CellId(30), cooling),
                },
            ],
            &[(MaterialId(0), &curve as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 2,
                max_substeps: 100,
            },
        )
        .unwrap()
        .advance_to(0., thermal_config(), || Ok(()))
        .unwrap()
        .0;
        let context = point(300.);
        let model = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![oriented(&card, &curve, &context, IDENTITY)],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let loads = model
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let independently_cooled = thermal
            .advance_to(2., thermal_config(), || Ok(()))
            .unwrap()
            .0;
        near(
            independently_cooled.temperature_k(CellId(30)).unwrap(),
            300.,
            1e-10,
        );
        assert!(matches!(
            model.advance_to(&thermal, 2., thermal_config(), &loads, &FREE, config(), cx),
            Err(SolidEnthalpyError::Material(_))
        ));
        let supported = model
            .advance_to(
                &thermal,
                0.25,
                thermal_config(),
                &loads,
                &FREE,
                config(),
                cx,
            )
            .unwrap();
        near(
            supported.thermal.temperature_k(CellId(30)).unwrap(),
            325.,
            1e-10,
        );
        assert_eq!(thermal.time_s(), 0.);
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg, 0.);
    });
}
