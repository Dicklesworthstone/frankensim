//! Material substitutions reach actual transient storage and diffusion.
mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{
    HeatSourceHistory, ThermalHistoryCapacity, ThermalHistoryMethod, ThermalHistoryProblem,
    ThermalHistoryRun,
};
use fs_conduction::transient::MaterialHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
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
use fs_qty::{Density, Temperature};
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
