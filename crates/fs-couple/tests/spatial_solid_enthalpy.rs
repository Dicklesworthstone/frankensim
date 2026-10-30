//! Independent P1 diffusion and affine thermoelastic oracles on the SAME mesh.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::bc::{ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::material::{CONDUCTIVITY_DIMS, LiquidMassFractionConductivity};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_couple::solid_enthalpy::spatial::{
    SpatialMaterialSource, SpatialSolidEnthalpy, SpatialSolidEnthalpyError,
};
use fs_couple::solid_enthalpy::{
    SolidEquilibriumConfig, ThermalExpansionModel, ThermalStiffnessModel,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_solid::loads3::VolumeLoad3;
use fs_solver::NewtonKrylovConfig;

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const FREE: [usize; 6] = [0, 1, 2, 4, 5, 8];

#[path = "spatial_solid_enthalpy/conductivity.rs"]
mod conductivity;
#[path = "spatial_solid_enthalpy/motion.rs"]
mod motion;
#[path = "spatial_solid_enthalpy/radiation.rs"]
mod radiation;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), NODES.to_vec()).unwrap()
}
fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
}
fn thermal_config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}
fn solid_config() -> SolidEquilibriumConfig {
    SolidEquilibriumConfig {
        static_solve: TetStaticSolveConfig::default(),
        mesh_frame: ContentHash([51; 32]),
        maximum_nodal_displacement_m: 0.1,
    }
}
fn native<'m, 'c>(
    mesh: &'m ConductionMesh,
    curve: &'c EquilibriumEnthalpyPhaseCurve,
    density: f64,
    cx: &Cx<'_>,
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        &[ReferenceEnthalpyMaterial {
            curve,
            reference_density_kg_m3: density,
        }],
        &[0],
        EnthalpyBudget {
            max_vertices: 4,
            max_elements: 1,
        },
    )
    .unwrap()
}
fn source<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
) -> SpatialMaterialSource<'a> {
    SpatialMaterialSource {
        curve,
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        reference_specific_enthalpy_j_kg: 0.,
        temperature_interval_k: [300., 500.],
        conductivity: PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
        conductivity_selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
        stiffness_model: ThermalStiffnessModel::Isotropic,
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
        expansion_model: ThermalExpansionModel::Isotropic,
    }
}
fn bind<'a, 'm: 'a, 'c: 'a>(
    native: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    source: SpatialMaterialSource<'a>,
    cx: &Cx<'_>,
) -> Result<SpatialSolidEnthalpy<'a, 'm, 'c>, SpatialSolidEnthalpyError> {
    SpatialSolidEnthalpy::try_new(
        native,
        vec![source],
        TetAssemblyBudget::standard(),
        solid_config().mesh_frame,
        0.1,
        cx,
    )
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

// rho(T) intentionally falls, while the reference inventory remains 6 kg/m3.
// alpha(T) doubles each 100 K, so free strain requires its integral, not alpha(T)*dT.
fn fixture(
    k: f64,
    young: f64,
    alpha: f64,
    cp: f64,
    mechanical_upper: f64,
    interior_young: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    let mut add = |key: PropertyKey, knots: Vec<(f64, f64)>, high: f64| {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.,
                    high,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic P1 identity fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(heat.clone(), vec![(300., cp), (500., cp)], 500.);
    add(density.clone(), vec![(300., 6.), (500., 2.)], 500.);
    add(
        PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
        vec![(300., k), (500., k)],
        500.,
    );
    add(
        PropertyKey::new("young_modulus", Pressure::DIMS),
        vec![(300., young), (400., interior_young), (500., 0.5 * young)],
        mechanical_upper,
    );
    add(
        PropertyKey::new("poisson_ratio", Dims([0; 6])),
        vec![(300., 0.25), (500., 0.25)],
        500.,
    );
    add(
        PropertyKey::new(
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
        ),
        vec![(300., alpha), (500., 3. * alpha)],
        500.,
    );
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "source substitution fixture".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(300.),
            upper: &point(500.),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}

#[test]
fn source_only_substitution_changes_native_diffusion_affine_motion_and_constrained_stress() {
    with_cx(|cx| {
        for (k, young, alpha, cp) in [(2., 1000., 1e-5, 10.), (8., 2000., 2e-5, 20.)] {
            let (card, caloric) = fixture(k, young, alpha, cp, 500., 0.75 * young);
            let mesh = mesh();
            let native = native(&mesh, caloric.curve(), 6., cx);
            let context = point(300.);
            let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
            let geometry = model.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let old = [100. * cp, 120. * cp, 80. * cp, 100. * cp];
            let original_bits = old.map(f64::to_bits);
            let boundary = boundary(&mesh);
            let input = ScalarField::uniform("heat source", 60.).unwrap();
            let result = model
                .advance(
                    &old,
                    1.,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx,
                )
                .unwrap();
            // Unit-tet P1: K=k/6 G'G, m_i=rho/24. (0,1,-1,0)
            // is an exact diffusion eigenvector; uniform source raises the mean.
            let rise = 100. + 60. / (6. * cp);
            let attenuation = 1. / (1. + 4. * k / (6. * cp));
            for (actual, perturbation) in result.thermal.temperature.iter().zip([0., 20., -20., 0.])
            {
                near(*actual, 300. + rise + perturbation * attenuation, 2e-8);
            }
            near(result.thermal.stored_energy_change_j, 10., 1e-8);
            near(result.thermal.energy_residual_j, 0., 1e-8);
            let free_strain = alpha * (rise + 0.005 * rise * rise);
            for (actual, reference) in result.equilibrium.deformed_mesh.nodes_m().iter().zip(NODES)
            {
                for c in 0..3 {
                    near(actual[c], reference[c] * (1. + free_strain), 2e-9);
                }
            }
            for stress in result.equilibrium.solution.stress.elements()[0].stress_mandel_pa() {
                near(*stress, 0., 1e-7);
            }
            near(result.state.materials()[0].density_kg_m3, 6., 0.);
            assert!(result.state.samples()[0].phase.bulk_density_kg_m3() < 6.);
            assert!(
                model
                    .conductivity()
                    .receipts()
                    .iter()
                    .all(|r| card.claims().claim(r.selected).is_some())
            );
            let all_fixed: Vec<_> = (0..12).collect();
            let clamped = model
                .advance(
                    &old,
                    1.,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &all_fixed,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let current_young = young * (1. - 0.0025 * rise);
            let stress = -current_young * free_strain / (1. - 2. * 0.25);
            for &actual in
                &clamped.equilibrium.solution.stress.elements()[0].stress_mandel_pa()[..3]
            {
                near(actual, stress, 1e-7);
            }
            assert_eq!(clamped.equilibrium.deformed_mesh.nodes_m(), &NODES);
            assert_eq!(old.map(f64::to_bits), original_bits);
            assert_eq!(native.reference_nodal_masses_kg(), &[0.25; 4]);
            assert!(core::ptr::eq(model.thermal().mesh(), &mesh));
        }
    });
}

#[test]
fn enthalpy_centroid_is_explicit_and_is_not_mean_nodal_temperature() {
    with_cx(|cx| {
        let (card, _) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Solid,
            [(0., 300.), (1000., 400.), (3000., 500.)]
                .into_iter()
                .map(|(h, t)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h,
                    temperature_k: t,
                    liquid_mass_fraction: 0.,
                    bulk_density_kg_m3: 6.,
                })
                .collect(),
        )
        .unwrap();
        let mesh = mesh();
        let native = native(&mesh, &curve, 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, &curve, &context), cx).unwrap();
        let state = model.resolve(&[0., 0., 0., 3000.], cx).unwrap();
        near(state.samples()[0].phase.specific_enthalpy_j_kg(), 750., 0.);
        near(state.samples()[0].phase.temperature_k(), 375., 0.);
        assert_ne!(state.samples()[0].phase.temperature_k(), 350.); // mean nodal T
    });
}

#[test]
fn every_corner_refuses_melt_and_source_escape_even_when_centroid_remains_admissible() {
    with_cx(|cx| {
        let (card, _) = fixture(2., 1000., 1e-5, 10., 500., 750.);
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
                    specific_enthalpy_j_kg: 1000.,
                    temperature_k: 400.,
                    liquid_mass_fraction: 0.,
                    bulk_density_kg_m3: 4.,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1200.,
                    temperature_k: 400.,
                    liquid_mass_fraction: 1.,
                    bulk_density_kg_m3: 3.,
                },
            ],
        )
        .unwrap();
        let mesh = mesh();
        let native = native(&mesh, &curve, 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, &curve, &context), cx).unwrap();
        assert!(matches!(
            model.resolve(&[0., 0., 0., 1100.], cx),
            Err(SpatialSolidEnthalpyError::Vertex { vertex: 3, .. })
        ));
        let mut restricted = source(&card, &curve, &context);
        restricted.temperature_interval_k = [300., 350.];
        let restricted = bind(&native, restricted, cx).unwrap();
        assert!(matches!(
            restricted.resolve(&[0., 0., 0., 600.], cx),
            Err(SpatialSolidEnthalpyError::Vertex { vertex: 3, .. })
        ));
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let initial = [950.; 4];
        let bits = initial.map(f64::to_bits);
        // Uniform heat crosses into the native latent plateau. The native heat
        // owner can solve it; this solid candidate must refuse its partial melt.
        let hot = ScalarField::uniform("heat source", 600.).unwrap();
        assert!(matches!(
            model.advance(
                &initial,
                1.,
                &boundary(&mesh),
                None,
                &hot,
                thermal_config(),
                &loads,
                &FREE,
                solid_config(),
                cx
            ),
            Err(SpatialSolidEnthalpyError::Vertex { .. })
        ));
        assert_eq!(initial.map(f64::to_bits), bits);
    });
}

#[test]
fn native_chart_inventory_and_continuous_source_support_are_admission_requirements() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let (other, foreign) = fixture(4., 1000., 1e-5, 10., 500., 750.);
        let mesh = mesh();
        let native_good = native(&mesh, caloric.curve(), 6., cx);
        let native_wrong_mass = native(&mesh, caloric.curve(), 3., cx);
        let context = point(300.);
        assert!(
            bind(
                &native_wrong_mass,
                source(&card, caloric.curve(), &context),
                cx
            )
            .is_err()
        );
        assert!(bind(&native_good, source(&other, foreign.curve(), &context), cx).is_err());
        assert!(bind(&native_good, source(&other, caloric.curve(), &context), cx).is_err());
        let native_scaled = native(&mesh, caloric.curve(), 6., cx)
            .with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(2., 2.).unwrap()],
            )
            .unwrap();
        assert!(bind(&native_scaled, source(&card, caloric.curve(), &context), cx).is_err());
        // Valid reference and two endpoint constants are insufficient when the
        // source ends early, or a source's interior knot has an invalid modulus.
        for (high, interior) in [(350., 750.), (500., -1.)] {
            let (card, caloric) = fixture(2., 1000., 1e-5, 10., high, interior);
            let native = native(&mesh, caloric.curve(), 6., cx);
            assert!(bind(&native, source(&card, caloric.curve(), &context), cx).is_err());
        }
    });
}

#[test]
fn a_refused_mechanical_endpoint_does_not_publish_or_modify_the_heat_history() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let old = [500.; 4];
        let bits = old.map(f64::to_bits);
        let input = ScalarField::uniform("heat source", 60.).unwrap();
        let tiny = SolidEquilibriumConfig {
            maximum_nodal_displacement_m: 1e-9,
            ..solid_config()
        };
        assert!(
            model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    tiny,
                    cx
                )
                .is_err()
        );
        assert_eq!(old.map(f64::to_bits), bits);
        // The same original state remains usable after refusal.
        assert!(
            model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx
                )
                .is_ok()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                model.advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    &cancelled
                ),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
            assert!(matches!(
                bind(
                    &native,
                    source(&card, caloric.curve(), &context),
                    &cancelled
                ),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
        });
        assert_eq!(old.map(f64::to_bits), bits);
    });
}

#[test]
fn heterogeneous_assignment_order_drives_distinct_elements_without_a_second_heat_store() {
    with_cx(|cx| {
        let (first, first_curve) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let (second, second_curve) = fixture(8., 2000., 3e-5, 20., 500., 1500.);
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|p| [p[0] + 3., p[1], p[2]]));
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
            nodes,
        )
        .unwrap();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: first_curve.curve(),
                    reference_density_kg_m3: 6.,
                },
                ReferenceEnthalpyMaterial {
                    curve: second_curve.curve(),
                    reference_density_kg_m3: 6.,
                },
            ],
            &[1, 0],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let context = point(300.);
        let model = SpatialSolidEnthalpy::try_new(
            &native,
            vec![
                source(&first, first_curve.curve(), &context),
                source(&second, second_curve.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            solid_config().mesh_frame,
            0.1,
            cx,
        )
        .unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut fixed = FREE.to_vec();
        fixed.extend(FREE.map(|dof| dof + 12));
        let input =
            ScalarField::nodal("local heating", 8, vec![60., 60., 60., 60., 0., 0., 0., 0.])
                .unwrap();
        let result = model
            .advance(
                &[0.; 8],
                1.,
                &boundary(&mesh),
                None,
                &input,
                thermal_config(),
                &loads,
                &fixed,
                solid_config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            result
                .state
                .samples()
                .iter()
                .map(|s| s.material_index)
                .collect::<Vec<_>>(),
            vec![1, 0]
        );
        assert_eq!(
            result.state.samples()[0].elastic.card_identity(),
            second.content_hash()
        );
        assert_eq!(
            result.state.samples()[1].elastic.card_identity(),
            first.content_hash()
        );
        for &t in &result.thermal.temperature[..4] {
            near(t, 300.5, 1e-8);
        }
        for &t in &result.thermal.temperature[4..] {
            near(t, 300., 1e-8);
        }
        let strain = 3e-5 * (0.5 + 0.005 * 0.5 * 0.5);
        near(
            result.equilibrium.deformed_mesh.nodes_m()[1][0],
            1. + strain,
            1e-9,
        );
        for vertex in 4..8 {
            assert_eq!(
                result.equilibrium.deformed_mesh.nodes_m()[vertex],
                mesh.positions()[vertex]
            );
        }
        near(result.thermal.stored_energy_change_j, 10., 1e-8);
    });
}

#[test]
fn changing_only_conductivity_propagates_through_spatial_heat_into_actual_thermal_stress() {
    with_cx(|cx| {
        // Equal-volume mirrored tetrahedra. The odd apex mode (-1,0,0,0,1)
        // has K/m eigenvalue 12*k/rho; the three shared-face nodes stay at 400 K.
        let mut nodes = NODES.to_vec();
        nodes.push([2. / 3.; 3]);
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(5, vec![[0, 1, 2, 3], [4, 1, 2, 3]]),
            nodes,
        )
        .unwrap();
        let context = point(300.);
        let mut responses = Vec::new();
        for k in [1., 8.] {
            let (card, caloric) = fixture(k, 1000., 1e-5, 10., 500., 750.);
            let native = HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &mesh,
                &[ReferenceEnthalpyMaterial {
                    curve: caloric.curve(),
                    reference_density_kg_m3: 6.,
                }],
                &[0, 0],
                EnthalpyBudget {
                    max_vertices: 5,
                    max_elements: 2,
                },
            )
            .unwrap();
            let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
            let geometry = model.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let old = [800., 1000., 1000., 1000., 1200.];
            let result = model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &ScalarField::Uniform(0.),
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let attenuation = 1. / (1. + 12. * k / 60.);
            near(
                result.thermal.temperature[0],
                400. - 20. * attenuation,
                2e-8,
            );
            near(
                result.thermal.temperature[4],
                400. + 20. * attenuation,
                2e-8,
            );
            for &t in &result.thermal.temperature[1..4] {
                near(t, 400., 2e-8);
            }
            near(
                result.state.samples()[0].phase.temperature_k(),
                400. - 5. * attenuation,
                2e-8,
            );
            near(
                result.state.samples()[1].phase.temperature_k(),
                400. + 5. * attenuation,
                2e-8,
            );
            near(result.thermal.stored_energy_change_j, 0., 1e-8);
            let stress: Vec<_> = result
                .equilibrium
                .solution
                .stress
                .elements()
                .iter()
                .flat_map(|element| *element.stress_mandel_pa())
                .collect();
            assert!(stress.iter().any(|s: &f64| s.abs() > 1e-5));
            responses.push((stress, result.equilibrium.deformed_mesh.nodes_m().to_vec()));
        }
        assert!(
            responses[0]
                .0
                .iter()
                .zip(&responses[1].0)
                .any(|(a, b)| (a - b).abs() > 1e-5)
        );
        assert!(
            responses[0]
                .1
                .iter()
                .flatten()
                .zip(responses[1].1.iter().flatten())
                .any(|(a, b)| (a - b).abs() > 1e-6)
        );
    });
}
