//! Manufactured source tensors exercise the real coupled owner; these are not
//! specimen data. Independent shear projection, momentum and caloric references.
use super::*;
use fs_couple::friction_heat::third_body::rigid::tensor::*;
use fs_couple::friction_heat::third_body::rigid::thermal::RigidLayerThermalTransport;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    PLASTIC_SHEAR_VISCOSITY_PROPERTY, SHEAR_LAYER_YIELD_PROPERTY,
};
use fs_material::visco_tensor::{
    TensorProjectionBudget,
    source::{
        ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
    },
};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Temperature, Time};

const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("omega", 1000.0)
        .unwrap()
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    material: ResolvedTensorMaterial,
}
fn fixture(softening: bool) -> Fixture {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims: Dims, value: f64, curve: Option<Vec<(f64, f64)>>| {
        let interpolation = if curve.is_some() {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        };
        let value = curve.map_or(PropertyValue::Scalar { value, dims }, |knots| {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            }
        });
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, 200_000.0),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured tensor layer coupling".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (name, value) in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
        .into_iter()
        .zip([30.0, 24.0, 18.0])
    {
        add(name, Pressure::DIMS, value, None);
    }
    for (name, value) in ORTHOTROPIC_POISSON_RATIO_PROPERTIES
        .into_iter()
        .zip([0.1, 0.2, 0.15])
    {
        add(name, Dims::NONE, value, None);
    }
    for (i, (name, value)) in ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES
        .into_iter()
        .zip([9.0, 7.0, 4.0])
        .enumerate()
    {
        add(
            name,
            Pressure::DIMS,
            value,
            if softening && i == 0 {
                Some(vec![(280.0, 9.2), (300.0, 9.0), (320.0, 8.8)])
            } else {
                None
            },
        );
    }
    add("small_strain_limit", Dims::NONE, 0.1, None);
    add("specific_heat_capacity", CP_DIMS, 1.0, None);
    add("density", Density::DIMS, CAP[2], None);
    add(SHEAR_LAYER_YIELD_PROPERTY, Pressure::DIMS, 1000.0, None);
    add(
        PLASTIC_SHEAR_VISCOSITY_PROPERTY,
        DynViscosity::DIMS,
        0.3,
        None,
    );
    add(
        "branch_tau",
        Time::DIMS,
        0.2,
        Some(vec![(280.0, 0.1), (300.0, 0.2), (320.0, 0.3)]),
    );
    let factor = [0.0, 0.0, 0.0, 2.0, 0.0, 1.0];
    let mut names = Vec::new();
    for i in 0..6 {
        for j in i..6 {
            let name = format!("branch_{i}_{j}");
            add(&name, Pressure::DIMS, factor[i] * factor[j], None);
            names.push(name);
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured tensor".into(),
            phase: "solid".into(),
            process: "native shear test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = ResolvedTensorMaterial::resolve(
        &card,
        &point(300.0),
        TensorPronyProperties {
            elasticity: TensorElasticity::Orthotropic,
            strain_limit: "small_strain_limit".into(),
            branches: vec![TensorBranchProperties {
                stiffness_upper: names.try_into().unwrap(),
                relaxation_time: "branch_tau".into(),
            }],
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0,
        1e-8,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        material,
    }
}
fn bind(f: &Fixture, angle: f64) -> Result<TensorLayerMaterial<'_>, TensorLayerSourceError> {
    let (s, c) = angle.sin_cos();
    TensorLayerMaterial::try_new(
        f.material.clone(),
        &f.card,
        &f.caloric,
        "T",
        [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]],
        geometry(),
        MaterialPropertySelection::SingleClaimOnly,
        1000.0,
        TensorProjectionBudget {
            max_coordinates: 2,
            max_scalar_terms: 10000,
        },
    )
}
fn inventory(f: &Fixture) -> DepositionSimulation<'_> {
    let original = thermal();
    let mut events = Vec::new();
    for cell in original.cells().values() {
        let mut cell = cell.clone();
        if cell.id == LAYER {
            cell.material_id = MaterialId(1);
        }
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        });
    }
    for link in original.links().values() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(*link),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &CALORIC as &dyn ThermalLaw),
            (MaterialId(1), f.caloric.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn rotated(k: [[f64; 2]; 2], angle: f64) -> [[f64; 2]; 2] {
    let (s, c) = angle.sin_cos();
    let r = [[c, -s], [s, c]];
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..2)
                .map(|a| (0..2).map(|b| r[i][a] * k[a][b] * r[j][b]).sum::<f64>())
                .sum()
        })
    })
}
fn product(k: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [
        k[0][0] * x[0] + k[0][1] * x[1],
        k[1][0] * x[0] + k[1][1] * x[1],
    ]
}
fn solve(a: [[f64; 2]; 2], b: [f64; 2]) -> [f64; 2] {
    let d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    [
        (a[1][1] * b[0] - a[0][1] * b[1]) / d,
        (a[0][0] * b[1] - a[1][0] * b[0]) / d,
    ]
}

#[test]
fn g1_original_tensor_source_changes_actual_face_momenta_memory_and_finite_heat() {
    let f = fixture(false);
    let angle = 0.4;
    let material = bind(&f, angle).unwrap();
    let heat = inventory(&f);
    let route = ThirdBodyHeatRoute::from_links(&heat, LAYER, [7, 8], K, [0.0; 2]).unwrap();
    let (_, law) = material.at_temperature(300.0).unwrap();
    let expected_k = rotated([[9.0, 0.0], [0.0, 4.0]], angle);
    let expected_arm = rotated([[2.0, 1.0], [1.0, 0.5]], angle);
    for i in 0..2 {
        for j in 0..2 {
            close(law.stiffness_tensor_n_m()[i][j], expected_k[i][j], 2e-13);
            close(
                law.tensor_relaxation().unwrap()[0].stiffness.matrix_n_m()[i][j],
                expected_arm[i][j],
                2e-13,
            );
        }
    }
    let history = RigidLayerHistory::new(&law, geometry(), 0.0).unwrap();
    let old = bodies(false, Vec3::new(0.0, 0.08, 0.0));
    let dt = 0.05;
    let power = 4.0;
    let candidate = step_tensor_heated_rigid_layer_candidate(
        &old,
        properties(false),
        &history,
        &heat,
        &route,
        &material,
        dt,
        10.0,
        law_solve(),
        heat_controls(),
        [0.0; 2],
        &[0.0; 3],
        |_, source| {
            Ok::<_, &str>(RigidLayerThermalTransport {
                conductivity_w_m_k: K,
                area_resistance_m2_k_w: [0.0; 2],
                transport: StepCoefficients {
                    sources_w: BTreeMap::from([(LAYER, power)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: source.clone(),
                },
            })
        },
        || false,
    )
    .unwrap();
    let report = candidate.report();
    let theta = report.heat.layer_temperature_k - 300.0;
    let tau = 0.2 + 0.005 * theta;
    let w = 1.0 / (1.0 + dt / (2.0 * tau));
    let effective =
        std::array::from_fn(|i| std::array::from_fn(|j| expected_k[i][j] + w * expected_arm[i][j]));
    let a = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                1.0 + dt * dt * (17.0 / 6.0) * effective[i][j] / 4.0
            } else {
                dt * dt * (17.0 / 6.0) * effective[i][j] / 4.0
            }
        })
    });
    let force = solve(a, product(effective, [0.08, 0.0]).map(|v| 0.5 * dt * v));
    for i in 0..2 {
        close(report.mechanical.constitutive.traction_n[i], force[i], 3e-9);
    }
    assert!(
        force[1].abs() > 0.001,
        "source orientation must cause transverse force"
    );
    let slip = [
        dt * 0.08 - dt * dt * (17.0 / 6.0) * force[0] / 2.0,
        -dt * dt * (17.0 / 6.0) * force[1] / 2.0,
    ];
    let mid = slip.map(|v| 0.5 * w * v);
    let memory_force = product(expected_arm, mid);
    let dissipation = dt / tau * (mid[0] * memory_force[0] + mid[1] * memory_force[1]);
    close(report.heat.dissipated_j, dissipation, 3e-10);
    assert!(dissipation > 0.0);
    let g = [0.8, 2.0 / 3.0];
    let capacity = CAP[2]
        + dt * g[0] * CAP[0] / (CAP[0] + dt * g[0])
        + dt * g[1] * CAP[1] / (CAP[1] + dt * g[1]);
    close(capacity * theta, dt * power + dissipation, 4e-10);
    close(report.combined_energy_residual_j, 0.0, 5e-10);
    assert!(report.thermoelastic.is_none());
    close(report.heat.reversible_heat_j, 0.0, 0.0);
    close(
        report.source.material.law().branches()[0].relaxation_time_s(),
        tau,
        1e-14,
    );
    assert_eq!(
        report.source.material.state().card_identity(),
        f.card.content_hash()
    );
    assert_eq!(report.source.plastic.card_identity(), f.card.content_hash());
    for side in 0..2 {
        let sign = if side == 0 { 1.0 } else { -1.0 };
        vector_close(
            candidate.bodies()[side]
                .linear_momentum_world()
                .sub(old[side].linear_momentum_world()),
            Vec3::new(0.0, sign * dt * force[0], sign * dt * force[1]),
            3e-10,
        );
    }
    vector_close(
        report.mechanical.angular_momentum_change_nms,
        Vec3::ZERO,
        3e-10,
    );
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(history.time_s(), 0.0);
}

#[test]
fn g4_tensor_softening_and_foreign_orientation_refuse_before_native_state_changes() {
    let soft = fixture(true);
    assert!(bind(&soft, 0.4).is_err());
    let f = fixture(false);
    let material = bind(&f, 0.4).unwrap();
    assert!(material.at_temperature(321.0).is_err());
    let heat = inventory(&f);
    let route = ThirdBodyHeatRoute::from_links(&heat, LAYER, [7, 8], K, [0.0; 2]).unwrap();
    let foreign = bind(&f, 0.0).unwrap();
    let (_, law) = foreign.at_temperature(300.0).unwrap();
    let history = RigidLayerHistory::new(&law, geometry(), 0.0).unwrap();
    let old = bodies(false, Vec3::new(0.0, 0.08, 0.0));
    let calls = Cell::new(0);
    let result = step_tensor_heated_rigid_layer_candidate(
        &old,
        properties(false),
        &history,
        &heat,
        &route,
        &material,
        0.05,
        10.0,
        law_solve(),
        heat_controls(),
        [0.0; 2],
        &[0.0; 3],
        |_, _| {
            calls.set(calls.get() + 1);
            Err::<RigidLayerThermalTransport<()>, _>("should not query transport")
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(HeatedRigidLayerError::Source(
            TensorLayerSourceError::Invalid(_)
        ))
    ));
    assert_eq!(calls.get(), 0);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(history.time_s(), 0.0);
}
