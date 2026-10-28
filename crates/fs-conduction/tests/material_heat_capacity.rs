//! Material substitutions reach actual transient storage and diffusion.
mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{
    HeatSourceHistory, ThermalHistoryCapacity, ThermalHistoryMethod, ThermalHistoryProblem,
    ThermalHistoryRun,
};
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::transient::{MaterialHeatCapacity, MaterialHeatCapacitySource};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ElementMaterials,
    LinearConfig, MaterialId, MaterialTable, ScalarField, ThermalBc, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use support::{with_cancelled_cx, with_cx};

fn card(rho: f64, cp: f64, k: f64, varying: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, rho, Density::DIMS),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            cp,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        ),
        (THERMAL_CONDUCTIVITY_PROPERTY, k, THERMAL_CONDUCTIVITY_DIMS),
    ] {
        let curve = varying && name == SPECIFIC_HEAT_CAPACITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curve {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, value), (320.0, 2.0 * value), (360.0, value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 360.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured thermal material; not measured data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn capacity(card: &MaterialCard, high: f64) -> MaterialHeatCapacity {
    MaterialHeatCapacity::from_card(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        280.0,
        high,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn mesh(two: bool) -> ConductionMesh {
    let mut points = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let mut tets = vec![[0, 1, 2, 3]];
    if two {
        points.extend([[2., 0., 0.], [3., 0., 0.], [2., 1., 0.], [2., 0., 1.]]);
        tets.push([4, 5, 6, 7]);
    }
    ConductionMesh::new(TetComplex::from_tets(points.len(), tets), points).unwrap()
}

fn controls() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1e-11,
            max_iterations: 512,
            restart: 16,
        },
        energy_tolerance_j: 1e-8,
    }
}

const CATALOG_SELECTORS: [&str; 2] = ["engineering-reference-model", "source-temperatures-matched"];

fn catalog_keys() -> [PropertyKey; 3] {
    [
        PropertyKey::new("density", Density::DIMS),
        PropertyKey::new("specific-heat-capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
    ]
}

fn catalog_validity(temperature: Option<[f64; 2]>) -> ValidityDomain {
    let mut domain = ValidityDomain::unconstrained();
    for (axis, value) in CATALOG_SELECTORS.into_iter().zip([1.0, 0.0]) {
        domain = domain.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value, value);
    }
    if let Some([low, high]) = temperature {
        domain = domain.with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            low,
            high,
        );
    }
    domain
}

fn catalog_point() -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            297.65,
        )
        .unwrap();
    for (axis, value) in CATALOG_SELECTORS.into_iter().zip([1.0, 0.0]) {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn catalog_card(values: [f64; 3], validity: ValidityDomain) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, value) in catalog_keys().into_iter().zip(values) {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured coefficients with the catalog's exact source schema"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic catalog-schema thermal material".into(),
            phase: "solid".into(),
            process: "numerical regression; no measured-material claim".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

#[test]
fn g1_catalog_source_keys_and_explicit_selectors_drive_actual_transient_fem() {
    let mesh = mesh(false);
    let keys = catalog_keys();
    let point = catalog_point();
    let source = MaterialHeatCapacitySource {
        density: &keys[0],
        specific_heat: &keys[1],
        temperature_axis: "temperature",
        fixed_source_selectors: &CATALOG_SELECTORS,
    };
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(1.5);
    for [rho, cp, conductivity] in [[3.0, 5.0, 2.0], [5.0, 6.0, 1.5]] {
        let card = catalog_card(
            [rho, cp, conductivity],
            catalog_validity(Some([297.15, 298.15])),
        );
        let storage = MaterialHeatCapacity::from_card_at_query_points(
            &card,
            source,
            &point,
            [297.15, 298.15],
            MaterialPropertySelection::PreferObservationBacked,
        )
        .unwrap();
        let k = storage
            .isotropic_conductivity_with_property(
                &card,
                &keys[2],
                SelectionPolicy::PreferObservationBacked,
            )
            .unwrap();
        assert_eq!(storage.value_j_per_m3_k(), rho * cp);
        assert!(
            storage
                .reference()
                .property("specific_heat_capacity")
                .is_none()
        );
        for state in [
            storage.reference(),
            &storage.endpoints()[0],
            &storage.endpoints()[1],
        ] {
            assert_eq!(state.axis_quantities(), point.axis_quantities());
            let query = state.query();
            assert_eq!(query.axes().get("engineering-reference-model"), Some(&1.0));
            assert_eq!(query.axes().get("source-temperatures-matched"), Some(&0.0));
            assert!(!query.axes().contains_key("T"));
            for property in state.properties() {
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        for receipt in k.receipts() {
            card.claims().verify_receipt(receipt).unwrap();
            assert_eq!(receipt.axis_quantities, *point.axis_quantities());
        }
        with_cx(|cx| {
            let engine = BackwardEuler::uniform_material(cx, &mesh, &storage).unwrap();
            let solved = engine
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &k,
                        source: &heat,
                        element_materials: None,
                    },
                    None,
                    &[297.65, 297.85, 297.45, 297.65],
                    0.5,
                    controls(),
                )
                .unwrap();
            // Independent P1 mode decay plus uniform caloric heating on V=1/6.
            let amplitude = 0.2 / (1.0 + 4.0 * conductivity * 0.5 / (rho * cp));
            let mean = 297.65 + 1.5 * 0.5 / (rho * cp);
            for (&actual, expected) in
                solved
                    .temperature
                    .iter()
                    .zip([mean, mean + amplitude, mean - amplitude, mean])
            {
                assert!((actual - expected).abs() < 1e-9);
            }
            assert!((solved.stored_energy_change_j - 0.125).abs() < 1e-8);
            assert_eq!(
                engine.material_capacities()[0].reference().card_identity(),
                card.content_hash()
            );
        });
        let other = catalog_card(
            [rho * 2.0, cp, conductivity],
            catalog_validity(Some([297.15, 298.15])),
        );
        assert!(
            storage
                .isotropic_conductivity_with_property(
                    &other,
                    &keys[2],
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn g0_catalog_source_admission_preserves_temperature_and_selector_limits() {
    let keys = catalog_keys();
    let point = catalog_point();
    let source = MaterialHeatCapacitySource {
        density: &keys[0],
        specific_heat: &keys[1],
        temperature_axis: "temperature",
        fixed_source_selectors: &CATALOG_SELECTORS,
    };
    let normal = catalog_validity(Some([297.15, 298.15]));
    let material = catalog_card([3.0, 5.0, 2.0], normal.clone());
    let resolve = |card: &MaterialCard, source, point: &QueryPoint, span| {
        MaterialHeatCapacity::from_card_at_query_points(
            card,
            source,
            point,
            span,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    assert!(resolve(&material, source, &point, [297.15, 298.15]).is_ok());
    for selectors in [
        &[][..],
        &CATALOG_SELECTORS[..1],
        &[CATALOG_SELECTORS[0], CATALOG_SELECTORS[0]][..],
    ] {
        assert!(
            resolve(
                &material,
                MaterialHeatCapacitySource {
                    fixed_source_selectors: selectors,
                    ..source
                },
                &point,
                [297.15, 298.15]
            )
            .is_err()
        );
    }
    for invalid_point in [
        point
            .clone()
            .with("engineering-reference-model", 0.0)
            .unwrap(),
        point
            .clone()
            .with("source-temperatures-matched", 1.0)
            .unwrap(),
        point
            .clone()
            .with_quantity("temperature", QuantitySpec::dimensional(Dims::NONE), 297.65)
            .unwrap(),
    ] {
        assert!(resolve(&material, source, &invalid_point, [297.15, 298.15]).is_err());
    }
    for validity in [
        catalog_validity(None),
        catalog_validity(Some([297.65, 297.65])),
        normal.clone().with_quantity(
            CATALOG_SELECTORS[0],
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            1.0,
        ),
    ] {
        let card = catalog_card([3.0, 5.0, 2.0], validity);
        assert!(resolve(&card, source, &point, [297.15, 298.15]).is_err());
    }
    assert!(resolve(&material, source, &point, [297.15, 298.16]).is_err());
    assert!(
        resolve(
            &material,
            MaterialHeatCapacitySource {
                specific_heat: &keys[2],
                ..source
            },
            &point,
            [297.15, 298.15]
        )
        .is_err()
    );
    let stress = QuantitySpec::semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static));
    let stressed = catalog_card(
        [3.0, 5.0, 2.0],
        normal.with_quantity("stress", stress, 1.0, 1.0),
    );
    let point = point.with_quantity("stress", stress, 1.0).unwrap();
    assert!(resolve(&stressed, source, &point, [297.15, 298.15]).is_err());
}

#[test]
fn g1_theta_material_swap_changes_diffusion_heating_and_the_exact_source_adjoint() {
    use fs_conduction::transient::{
        FinalStateFunctional, HeatCapacity, TransientConfig,
        march_with_capacity_and_source_scale_gradient,
    };
    let mesh = mesh(false);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let functional = FinalStateFunctional::probe(4, 1).unwrap();
    let heat = ScalarField::Uniform(30.0);
    for material in [card(3.0, 5.0, 2.0, false), card(5.0, 6.0, 1.5, false)] {
        let storage = capacity(&material, 360.0);
        let k = storage
            .isotropic_conductivity(&material, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        let value = storage.value_j_per_m3_k();
        let conductivity = k.tensor_at(300.0).unwrap()[0][0];
        for theta in [0.5, 1.0] {
            let config = TransientConfig::new(theta, 0.2, controls().linear).unwrap();
            with_cx(|cx| {
                let (solved, gradient) = march_with_capacity_and_source_scale_gradient(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &k,
                        source: &heat,
                        element_materials: None,
                    },
                    HeatCapacity::MaterialUniform(&storage),
                    &config,
                    &[300.0, 310.0, 290.0, 300.0],
                    &functional,
                    3,
                )
                .unwrap();
                // Independent P1 eigenmode plus uniform caloric balance.
                let z = 4.0 * conductivity * 0.2 / value;
                let ratio = (1.0 - (1.0 - theta) * z) / (1.0 + theta * z);
                let amplitude = 10.0 * ratio.powi(3);
                let rise = 30.0 * 0.6 / value;
                assert!((solved.temperature[0] - 300.0 - rise).abs() < 1e-8);
                assert!((solved.temperature[1] - 300.0 - rise - amplitude).abs() < 1e-8);
                assert!((solved.temperature[2] - 300.0 - rise + amplitude).abs() < 1e-8);
                assert!((gradient - rise).abs() < 1e-9);
            });
        }
    }
}

#[test]
fn g1_theta_heterogeneous_capacity_and_adjoint_use_the_same_material_ids() {
    use fs_conduction::transient::{
        FinalStateFunctional, HeatCapacity, TransientConfig,
        march_with_capacity_and_source_scale_gradient,
    };
    let mesh = mesh(true);
    let a = card(2.0, 10.0, 2.0, false);
    let b = card(3.0, 20.0, 4.0, false);
    let ca = capacity(&a, 360.0);
    let cb = capacity(&b, 360.0);
    let ka = ca
        .isotropic_conductivity(&a, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let kb = cb
        .isotropic_conductivity(&b, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let assigned = ElementMaterials::new(
        MaterialTable::new([(MaterialId(7), ka.clone()), (MaterialId(2), kb)]).unwrap(),
        vec![MaterialId(7), MaterialId(2)],
    )
    .unwrap();
    let rows = [(MaterialId(2), cb), (MaterialId(7), ca)];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(30.0);
    let config = TransientConfig::crank_nicolson(0.2, controls().linear).unwrap();
    let functional =
        FinalStateFunctional::new(vec![1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0]).unwrap();
    with_cx(|cx| {
        let (solved, gradient) = march_with_capacity_and_source_scale_gradient(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &ka,
                source: &heat,
                element_materials: Some(&assigned),
            },
            HeatCapacity::MaterialById(&rows),
            &config,
            &[300.0; 8],
            &functional,
            3,
        )
        .unwrap();
        for (temperatures, rise) in solved.temperature.chunks_exact(4).zip([0.9, 0.3]) {
            assert!(
                temperatures
                    .iter()
                    .all(|&t| (t - 300.0 - rise).abs() < 1e-8)
            );
        }
        assert!((gradient - 1.5).abs() < 1e-9);
    });
}

#[test]
fn g1_source_adjoint_does_not_scale_robin_or_neumann_boundary_heat() {
    use fs_conduction::transient::{
        FinalStateFunctional, HeatCapacity, TransientConfig, TransientProblem,
        VolumetricHeatCapacity, march_with_capacity, march_with_capacity_and_source_scale_gradient,
        source_scale_gradient,
    };
    let mesh = mesh(false);
    let material = card(2.0, 10.0, 2.0, false);
    let storage = capacity(&material, 360.0);
    let k = storage
        .isotropic_conductivity(&material, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let declared_k = ConductivityModel::isotropic_declared(2.0).unwrap();
    let functional = FinalStateFunctional::probe(4, 1).unwrap();
    let config = TransientConfig::crank_nicolson(0.1, controls().linear).unwrap();
    for bc in [
        ThermalBc::robin(2.0, 320.0).unwrap(),
        ThermalBc::neumann(0.5).unwrap(),
    ] {
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("environment", |_| true, bc)
            .unwrap()
            .finish()
            .unwrap();
        with_cx(|cx| {
            let zero = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &k,
                source: &zero,
                element_materials: None,
            };
            let (_, zero_gradient) = march_with_capacity_and_source_scale_gradient(
                cx,
                problem,
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &[300.0; 4],
                &functional,
                3,
            )
            .unwrap();
            assert_eq!(zero_gradient, 0.0);
            // The correction also reaches the original declared-capacity API.
            assert_eq!(
                source_scale_gradient(
                    cx,
                    TransientProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &declared_k,
                        source: &zero,
                        capacity: VolumetricHeatCapacity::declared(20.0).unwrap(),
                    },
                    &config,
                    &functional,
                    3,
                )
                .unwrap(),
                0.0
            );
            let heat = ScalarField::Uniform(30.0);
            let (_, gradient) = march_with_capacity_and_source_scale_gradient(
                cx,
                ConductionProblem {
                    source: &heat,
                    ..problem
                },
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &[300.0; 4],
                &functional,
                3,
            )
            .unwrap();
            let objective = |scale| {
                let source = ScalarField::Uniform(30.0 * scale);
                let solved = march_with_capacity(
                    cx,
                    ConductionProblem {
                        source: &source,
                        ..problem
                    },
                    HeatCapacity::MaterialUniform(&storage),
                    &config,
                    &[300.0; 4],
                    3,
                )
                .unwrap();
                functional.evaluate(&solved.temperature)
            };
            let difference = (objective(1.01) - objective(0.99)) / 0.02;
            assert!((gradient - difference).abs() < 1e-7);
            assert!(gradient > 0.0);
        });
    }
}

#[test]
fn g4_theta_checks_every_source_endpoint_before_returning_a_field_or_gradient() {
    use fs_conduction::transient::{
        FinalStateFunctional, HeatCapacity, TransientConfig, march_with_capacity,
        march_with_capacity_and_source_scale_gradient,
    };
    let mesh = mesh(false);
    let material = card(2.0, 10.0, 2.0, false);
    let storage = capacity(&material, 300.5);
    let wide = capacity(&material, 360.0);
    let k = ConductivityModel::isotropic_declared(2.0).unwrap();
    let narrow_k = ConductivityModel::isotropic(
        fs_conduction::ConductivityTable::declared_curve(vec![(280.0, 2.0), (300.2, 2.0)]).unwrap(),
    );
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(30.0);
    let initial = [300.0; 4];
    let functional = FinalStateFunctional::probe(4, 0).unwrap();
    let config = TransientConfig::crank_nicolson(0.25, controls().linear).unwrap();
    with_cx(|cx| {
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &k,
            source: &heat,
            element_materials: None,
        };
        assert!(
            march_with_capacity(
                cx,
                problem,
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &initial,
                1,
            )
            .is_ok()
        );
        assert!(matches!(
            march_with_capacity_and_source_scale_gradient(
                cx,
                problem,
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &initial,
                &functional,
                2,
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
        assert!(matches!(
            march_with_capacity(
                cx,
                ConductionProblem {
                    material: &narrow_k,
                    ..problem
                },
                HeatCapacity::MaterialUniform(&wide),
                &config,
                &initial,
                1,
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
        let prescribed = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "hot",
                |face| face.centroid[0] == 0.0,
                ThermalBc::dirichlet(301.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        assert!(matches!(
            march_with_capacity(
                cx,
                ConductionProblem {
                    boundary: &prescribed,
                    ..problem
                },
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &initial,
                0,
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
        assert_eq!(initial, [300.0; 4]);
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            march_with_capacity(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &k,
                    source: &heat,
                    element_materials: None,
                },
                HeatCapacity::MaterialUniform(&storage),
                &config,
                &initial,
                1,
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g1_one_card_swap_changes_storage_and_actual_diffusion_rate() {
    let mesh = mesh(false);
    let original = card(3.0, 5.0, 2.0, false);
    let swapped = card(5.0, 6.0, 1.5, false);
    with_cx(|cx| {
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        for card in [&original, &swapped] {
            let source = capacity(card, 360.0);
            let k = source
                .isotropic_conductivity(card, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            let engine = BackwardEuler::uniform_material(cx, &mesh, &source).unwrap();
            let zero = ScalarField::Uniform(0.0);
            let solved = engine
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &k,
                        source: &zero,
                        element_materials: None,
                    },
                    None,
                    &[300.0, 310.0, 290.0, 300.0],
                    0.5,
                    controls(),
                )
                .unwrap();
            // On a right unit tet, the [0,1,-1,0] P1 mode has decay rate
            // 4*k/(rho*cp): K*v=k*V*v and lumped C=(rho*cp)*V/4.
            let coefficient = k.tensor_at(300.0).unwrap()[0][0];
            let amplitude = 10.0 / (1.0 + 4.0 * coefficient * 0.5 / source.value_j_per_m3_k());
            assert!((solved.temperature[1] - 300.0 - amplitude).abs() < 1e-9);
            assert!((solved.temperature[2] - 300.0 + amplitude).abs() < 1e-9);
            assert!(solved.stored_energy_change_j.abs() < 1e-9);
            assert_eq!(
                engine.material_capacities()[0].reference().card_identity(),
                card.content_hash()
            );
            for endpoint in source.endpoints() {
                for property in endpoint.properties() {
                    card.claims()
                        .verify_receipt(&property.answer().receipt)
                        .unwrap();
                    assert_eq!(
                        property.answer().receipt.selected,
                        source
                            .reference()
                            .property(property.requirement().name())
                            .unwrap()
                            .answer()
                            .receipt
                            .selected
                    );
                }
            }
        }
    });
    assert!(
        capacity(&original, 360.0)
            .isotropic_conductivity(&swapped, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
}

#[test]
fn g1_heterogeneous_history_uses_each_materials_heat_storage_and_retains_sources() {
    let mesh = mesh(true);
    let a = card(2.0, 10.0, 2.0, false);
    let b = card(3.0, 20.0, 4.0, false);
    let ca = capacity(&a, 360.0);
    let cb = capacity(&b, 360.0);
    let ka = ca
        .isotropic_conductivity(&a, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let kb = cb
        .isotropic_conductivity(&b, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let assigned = ElementMaterials::new(
        MaterialTable::new([(MaterialId(7), ka.clone()), (MaterialId(2), kb)]).unwrap(),
        vec![MaterialId(7), MaterialId(2)],
    )
    .unwrap();
    let rows = [(MaterialId(2), cb), (MaterialId(7), ca)];
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            Vec::new(),
            0.0,
        )
        .unwrap();
        for method in [
            ThermalHistoryMethod::Linear,
            ThermalHistoryMethod::Nonlinear(NonlinearStepConfig::default()),
        ] {
            let mut run = ThermalHistoryRun::new_with_capacity(
                cx,
                ThermalHistoryProblem {
                    boundary: &history,
                    material: &ka,
                    element_materials: Some(&assigned),
                    interfaces: None,
                    source: HeatSourceHistory::Fixed(ScalarField::Uniform(30.0)),
                },
                ThermalHistoryCapacity::MaterialById(&rows),
                &[300.0; 8],
                2,
                method,
            )
            .unwrap();
            assert_eq!(run.material_capacities().len(), 2);
            let solved = run.advance_to(cx, 0.5, controls()).unwrap();
            for (temperatures, expected) in solved.temperature.chunks_exact(4).zip([300.75, 300.25])
            {
                assert!(temperatures.iter().all(|&t| (t - expected).abs() < 1e-9));
            }
            // Two tets each V=1/6, source=30 W/m³ for half a second.
            assert!((solved.stored_energy_change_j - 5.0).abs() < 1e-8);
        }
    });
}

#[test]
fn g4_out_of_source_domain_does_not_publish_a_history_or_bypass_coupled_admission() {
    let mesh = mesh(false);
    let card = card(2.0, 10.0, 2.0, false);
    let source = capacity(&card, 300.5);
    let k = source
        .isotropic_conductivity(&card, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            Vec::new(),
            0.0,
        )
        .unwrap();
        for method in [
            ThermalHistoryMethod::Linear,
            ThermalHistoryMethod::Nonlinear(NonlinearStepConfig::default()),
        ] {
            let problem = || ThermalHistoryProblem {
                boundary: &history,
                material: &k,
                element_materials: None,
                interfaces: None,
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(30.0)),
            };
            assert!(
                ThermalHistoryRun::new_with_capacity(
                    cx,
                    problem(),
                    ThermalHistoryCapacity::MaterialUniform(&source),
                    &[300.6; 4],
                    2,
                    method
                )
                .is_err()
            );
            let mut run = ThermalHistoryRun::new_with_capacity(
                cx,
                problem(),
                ThermalHistoryCapacity::MaterialUniform(&source),
                &[300.0; 4],
                2,
                method,
            )
            .unwrap();
            assert!(run.advance_to(cx, 0.5, controls()).is_err());
            assert_eq!(run.time_s(), 0.0);
            assert_eq!(run.accepted_steps(), 0);
            assert_eq!(run.temperature(), &[300.0; 4]);
            assert!(run.advance_to(cx, 0.1, controls()).is_ok());
        }
        let engine = BackwardEuler::uniform_material(cx, &mesh, &source).unwrap();
        // A different, wider-range conductivity cannot widen caloric support.
        let wide = ConductivityModel::isotropic_declared(2.0).unwrap();
        let zero = ScalarField::Uniform(0.0);
        let boundary = history.at(cx, 0.0).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &wide,
            source: &zero,
            element_materials: None,
        };
        assert!(matches!(
            engine.evaluate_prescribed_endpoint(
                cx,
                problem,
                None,
                &[300.0; 4],
                &[300.0, 300.6, 300.0, 300.0],
                0.1,
                controls(),
                None
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("all", |_| true, ThermalBc::dirichlet(301.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        assert!(
            engine
                .advance_prescribed(
                    cx,
                    ConductionProblem {
                        boundary: &boundary,
                        ..problem
                    },
                    None,
                    &[300.0; 4],
                    0.1,
                    controls()
                )
                .is_err()
        );
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            BackwardEuler::uniform_material(cx, &mesh, &source),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g0_varying_nonpositive_and_out_of_support_sources_are_not_frozen_or_replaced() {
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    for card in [
        card(2.0, 10.0, 2.0, true),
        card(0.0, 10.0, 2.0, false),
        card(2.0, -1.0, 2.0, false),
    ] {
        assert!(
            MaterialHeatCapacity::from_card(
                &card,
                &point,
                280.0,
                360.0,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let card = card(2.0, 10.0, 2.0, false);
    assert!(
        MaterialHeatCapacity::from_card(
            &card,
            &point,
            280.0,
            361.0,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}
