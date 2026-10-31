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
        ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity,
        TensorKelvinViscosityProperties, TensorPronyProperties,
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
    fixture_with_kelvin(softening, false)
}
fn fixture_with_kelvin(softening: bool, kelvin: bool) -> Fixture {
    fixture_with_thermal(softening, kelvin, false)
}
fn fixture_with_thermal(softening: bool, kelvin: bool, affine: bool) -> Fixture {
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
        add(name, Pressure::DIMS, value, affine.then(||
            vec![(280.0,1.05*value),(300.0,value),(320.0,0.95*value)]));
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
            if affine {
                Some(vec![(280.0,1.05*value),(300.0,value),(320.0,0.95*value)])
            } else if softening && i == 0 {
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
            let value = factor[i]*factor[j] + if affine && i==5 && j==5 {2.0} else {0.0};
            let slope = if affine { match (i,j) {
                (3,3) => -0.01, (3,5) => 0.006, (5,5) => 0.02, _ => 0.0,
            }} else {0.0};
            add(&name, Pressure::DIMS, value, affine.then(||
                vec![(280.0,value-20.0*slope),(300.0,value),(320.0,value+20.0*slope)]));
            names.push(name);
        }
    }
    let viscosity_names: [String; 21] = std::array::from_fn(|i| format!("viscosity_{i:02}"));
    if kelvin {
        let factor = [0.0, 0.0, 0.0, 0.2, 0.0, -0.1];
        let mut entry = 0;
        for i in 0..6 {
            for j in i..6 {
                let value = factor[i] * factor[j];
                add(
                    &viscosity_names[entry],
                    DynViscosity::DIMS,
                    value,
                    Some(vec![
                        (280.0, 0.5 * value),
                        (300.0, value),
                        (320.0, 1.5 * value),
                    ]),
                );
                entry += 1;
            }
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
    let mut material = ResolvedTensorMaterial::resolve(
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
    if kelvin {
        material = material
            .with_kelvin_viscosity(
                &card,
                TensorKelvinViscosityProperties {
                    viscosity_upper: viscosity_names,
                },
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
    }
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

fn young_shear_fixture(specific_heat_capacity: f64) -> Fixture {
    let original = fixture_with_thermal(false, true, true);
    let mut claims = ClaimSet::new();
    for (_, value) in original.card.claims().claims_ordered() {
        let mut claim = value.clone();
        if claim.key.name() == "specific_heat_capacity" {
            claim.value = PropertyValue::Scalar {
                value: specific_heat_capacity,
                dims: CP_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    for (name, value, slope) in [
        ("young_modulus", 12.0_f64, -0.012_f64),
        ("shear_modulus", 4.8, 0.006),
    ] {
        let mut claim = original
            .card
            .claims()
            .claims_ordered()
            .next()
            .unwrap()
            .1
            .clone();
        claim.key = PropertyKey::new(name, Pressure::DIMS);
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [280.0, 300.0, 320.0]
                .into_iter()
                .map(|t| (t, slope.mul_add(t - 300.0, value)))
                .collect(),
            dims: Pressure::DIMS,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured independent Young/shear tensor".into(),
            phase: "solid".into(),
            process: "native shear test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut selectors = original.material.properties().clone();
    selectors.elasticity = TensorElasticity::IsotropicYoungShear;
    let material = ResolvedTensorMaterial::resolve(
        &card,
        &point(300.0),
        selectors,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_kelvin_viscosity(
        &card,
        TensorKelvinViscosityProperties {
            viscosity_upper: std::array::from_fn(|i| format!("viscosity_{i:02}")),
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
        -20.0 * specific_heat_capacity,
        1e-8,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        material,
    }
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
    assert_tensor_source_run(false);
}

#[test]
fn g1_kelvin_source_changes_actual_rigid_momenta_and_heat_at_the_same_trial_temperature() {
    assert_tensor_source_run(true);
}

fn assert_tensor_source_run(kelvin: bool) {
    let f = fixture_with_kelvin(false, kelvin);
    let angle = 0.4;
    let material = bind(&f, angle).unwrap();
    let heat = inventory(&f);
    let route = ThirdBodyHeatRoute::from_links(&heat, LAYER, [7, 8], K, [0.0; 2]).unwrap();
    let (_, law) = material.at_temperature(300.0).unwrap();
    let expected_k = rotated([[9.0, 0.0], [0.0, 4.0]], angle);
    let expected_arm = rotated([[2.0, 1.0], [1.0, 0.5]], angle);
    let expected_c = rotated([[0.02, -0.01], [-0.01, 0.005]], angle);
    for i in 0..2 {
        for j in 0..2 {
            close(law.stiffness_tensor_n_m()[i][j], expected_k[i][j], 2e-13);
            close(
                law.tensor_relaxation().unwrap()[0].stiffness.matrix_n_m()[i][j],
                expected_arm[i][j],
                2e-13,
            );
            if kelvin {
                close(
                    law.tensor_kelvin_viscosity().unwrap().matrix_n_s_m()[i][j],
                    expected_c[i][j],
                    2e-15,
                );
            }
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
    let viscosity_scale = if kelvin { 1.0 + 0.025 * theta } else { 0.0 };
    let effective = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            expected_k[i][j]
                + w * expected_arm[i][j]
                + 2.0 * viscosity_scale * expected_c[i][j] / dt
        })
    });
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
    let viscous_force = product(expected_c, slip).map(|v| viscosity_scale * v / dt);
    let kelvin_heat = slip[0] * viscous_force[0] + slip[1] * viscous_force[1];
    let dissipation =
        dt / tau * (mid[0] * memory_force[0] + mid[1] * memory_force[1]) + kelvin_heat;
    close(report.heat.dissipated_j, dissipation, 3e-10);
    assert!(dissipation > 0.0);
    if kelvin {
        assert!(kelvin_heat > 0.0);
        let source = report.source.material.kelvin_source().unwrap();
        close(
            source.state().state_coordinate("T").unwrap(),
            report.heat.layer_temperature_k,
            1e-14,
        );
        close(
            report.source.material.law().kv_pa_s().unwrap()[3][3],
            0.04 * viscosity_scale,
            1e-15,
        );
        assert_eq!(source.state().card_identity(), f.card.content_hash());
        assert_eq!(source.state().properties().len(), 21);
    }
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
fn g1_affine_source_tensors_drive_original_momenta_entropy_memory_and_finite_heat() {
    assert_thermal_source_run(false, false);
}

#[test]
fn g1_independent_young_shear_reaches_original_rigid_entropy_memory_and_finite_heat() {
    assert_thermal_source_run(true, false);
}

#[test]
fn g1_independent_orthotropic_shear_drives_continuing_rigid_entropy_memory_and_finite_heat() {
    assert_thermal_source_run(false, true);
}

fn assert_thermal_source_run(young_shear: bool, independent_orthotropic: bool) {
    let f = if young_shear {
        young_shear_fixture(1.0)
    } else if independent_orthotropic {
        fixture_with_kelvin(true, true)
    } else {
        fixture_with_thermal(false, true, true)
    };
    for angle in [0.4, -0.55] {
        let material = bind(&f, angle).unwrap();
        let mut heat = inventory(&f);
        let route = ThirdBodyHeatRoute::from_links(&heat, LAYER, [7, 8], K, [0.0; 2]).unwrap();
        let (_, law) = material.at_temperature(300.0).unwrap();
        let mut history = RigidLayerHistory::new(&law, geometry(), 0.0).unwrap();
        let mut moving = bodies(false,Vec3::new(0.0,0.08,0.0));
        let k0 = if young_shear { [[4.8,0.0],[0.0,4.8]] }
            else { rotated([[9.0,0.0],[0.0,4.0]],angle) };
        let b0 = rotated([[2.0,1.0],[1.0,if independent_orthotropic {0.5} else {1.5}]],angle);
        let dk = if young_shear { [[0.006,0.0],[0.0,0.006]] }
            else if independent_orthotropic { rotated([[-0.01,0.0],[0.0,0.0]],angle) }
            else { k0.map(|row| row.map(|v| -0.0025*v)) };
        let db = if independent_orthotropic { [[0.0;2];2] }
            else { rotated([[-0.005,0.003],[0.003,0.01]],angle) };
        let c0 = rotated([[0.02,-0.01],[-0.01,0.005]],angle);
        for i in 0..2 { for j in 0..2 {
            close(material.thermal_slopes().equilibrium_n_m_k[i][j],dk[i][j],2e-15);
            close(material.thermal_slopes().relaxation_n_m_k[0][i][j],db[i][j],2e-15);
        }}
        let at = |base: [[f64;2];2], slope: [[f64;2];2], t: f64| {
            std::array::from_fn(|i| std::array::from_fn(|j| base[i][j]+(t-300.0)*slope[i][j]))
        };
        let quadratic = |k: [[f64;2];2], x: [f64;2]| {
            k[0][0]*x[0]*x[0]+2.0*k[0][1]*x[0]*x[1]+k[1][1]*x[1]*x[1]
        };
        let energy = |state: &fs_tribo::layer_contact::ThirdBodyLayerState,t: f64| {
            0.5*quadratic(at(k0,dk,t),state.elastic_slip_m())
                +0.5*quadratic(at(b0,db,t),state.relaxation_slip_m()[0])
        };
        let entropy = |state: &fs_tribo::layer_contact::ThirdBodyLayerState| {
            -0.5*quadratic(dk,state.elastic_slip_m())
                -0.5*quadratic(db,state.relaxation_slip_m()[0])
        };
        let dt = 0.05;
        let power = 4.0;
        for step in 1..=3 {
            let old_temperature = history.law().temperature_kelvin();
            let old_state = history.deformation().clone();
            let candidate = step_tensor_heated_rigid_layer_candidate(&moving,properties(false),
                &history,&heat,&route,&material,step as f64*dt,10.0,law_solve(),heat_controls(),
                [0.0;2],&[0.0;3],|_,source| Ok::<_,&str>(RigidLayerThermalTransport {
                    conductivity_w_m_k: K, area_resistance_m2_k_w: [0.0;2],
                    transport: StepCoefficients { sources_w: BTreeMap::from([(LAYER,power)]),
                        conductances_w_k: BTreeMap::new(), payload: source.clone() },
                }),||false).unwrap();
            let report = candidate.report();
            let account = report.thermoelastic.unwrap();
            let temperature = report.heat.layer_temperature_k;
            let new_state = candidate.history().deformation();
            let k = at(k0,dk,temperature);
            let b = at(b0,db,temperature);
            let c = c0.map(|row| row.map(|v| v*(1.0+0.025*(temperature-300.0))));
            let tau = 0.2+0.005*(temperature-300.0);
            let endpoint_slopes = material.thermal_slopes_at_temperature(temperature).unwrap();
            for i in 0..2 { for j in 0..2 {
                close(candidate.history().law().stiffness_tensor_n_m()[i][j],k[i][j],3e-13);
                close(candidate.history().law().tensor_relaxation().unwrap()[0]
                    .stiffness.matrix_n_m()[i][j],b[i][j],3e-13);
                close(endpoint_slopes.equilibrium_n_m_k[i][j],dk[i][j],2e-15);
                close(endpoint_slopes.relaxation_n_m_k[0][i][j],db[i][j],2e-15);
            }}
            // The first held-force solution is independently inverted at the
            // ORIGINAL endpoint trial temperature, including full cross terms.
            if step == 1 {
                let w = 1.0/(1.0+dt/(2.0*tau));
                let effective = std::array::from_fn(|i| std::array::from_fn(|j|
                    k[i][j]+w*b[i][j]+2.0*c[i][j]/dt));
                let a = std::array::from_fn(|i| std::array::from_fn(|j|
                    if i==j {1.0+dt*dt*(17.0/6.0)*effective[i][j]/4.0}
                    else {dt*dt*(17.0/6.0)*effective[i][j]/4.0}));
                let force = solve(a,product(effective,[0.08,0.0]).map(|v| 0.5*dt*v));
                for i in 0..2 {close(report.mechanical.constitutive.traction_n[i],force[i],3e-9);}
                assert!(force[1].abs()> if young_shear {1e-5} else {0.001});
            } else {
                assert_ne!(account.coefficient_free_energy_change_j,0.0);
            }
            let old_s = entropy(&old_state);
            let new_s = entropy(new_state);
            let coefficient = energy(&old_state,temperature)-energy(&old_state,old_temperature);
            let delta_ts = temperature*new_s-old_temperature*old_s;
            close(account.coefficient_free_energy_change_j,coefficient,2e-14);
            close(account.mechanical_entropy_change_j_k,new_s-old_s,2e-16);
            close(account.entropy_energy_change_j,delta_ts,2e-14);
            close(account.internal_energy_change_j,
                energy(new_state,temperature)+temperature*new_s
                    -energy(&old_state,old_temperature)-old_temperature*old_s,2e-14);
            close(account.reversible_heat_j,-coefficient-delta_ts,2e-14);
            assert_eq!(new_state.plastic_slip_m(),old_state.plastic_slip_m());
            let mid = std::array::from_fn(|i| 0.5*(new_state.relaxation_slip_m()[0][i]
                +old_state.relaxation_slip_m()[0][i]));
            let delta_e = std::array::from_fn(|i| new_state.elastic_slip_m()[i]-old_state.elastic_slip_m()[i]);
            let dissipation = dt/tau*quadratic(b,mid)+quadratic(c,delta_e)/dt;
            close(report.heat.dissipated_j,dissipation,2e-13);
            close(report.heat.thermal.stored_energy_change_j,
                dt*power+dissipation-coefficient-delta_ts,5e-10);
            close(report.combined_energy_residual_j,0.0,5e-10);
            close(report.source.material.state().state_coordinate("T").unwrap(),temperature,1e-14);
            close(report.source.material.kelvin_source().unwrap().state()
                .state_coordinate("T").unwrap(),temperature,1e-14);
            if young_shear {
                assert!(report.source.material.state().property("poisson_ratio").is_none());
                for (name,value,slope) in [("young_modulus",12.0,-0.012),
                    ("shear_modulus",4.8,0.006)] {
                    let source = report.source.material.state().property(name).unwrap();
                    close(source.value_si(),value+slope*(temperature-300.0),2e-14);
                    assert_eq!(source.answer().receipt.selected,
                        f.material.state().property(name).unwrap().answer().receipt.selected);
                }
            } else if independent_orthotropic {
                let name = ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0];
                let source = report.source.material.state().property(name).unwrap();
                close(source.value_si(),9.0-0.01*(temperature-300.0),2e-14);
                assert_eq!(source.answer().receipt.selected,
                    f.material.state().property(name).unwrap().answer().receipt.selected);
            }
            assert_eq!(history.deformation(),&old_state);
            let (next_bodies,next_history,next_heat,_) = candidate.into_parts();
            moving = next_bodies; history = next_history; heat = next_heat;
        }
    }
}

#[test]
fn g4_nonlinear_layer_requires_positive_capacity_over_the_complete_source_strain_domain() {
    let f = young_shear_fixture(0.001);
    let before = f.material.clone();
    assert!(matches!(
        bind(&f, 0.4),
        Err(TensorLayerSourceError::Invalid(
            "tensor layer cannot establish positive total fixed-strain heat capacity over its original domain"
        ))
    ));
    assert_eq!(f.material, before);
    let sufficient = young_shear_fixture(1.0);
    let material = bind(&sufficient, 0.4).unwrap();
    assert!(material.thermal_slopes_at_temperature(321.0).is_err());
    assert!(material.thermal_slopes_at_temperature(f64::NAN).is_err());
}

#[test]
fn g4_independent_equilibrium_path_keeps_source_domain_and_foreign_orientation_refusals() {
    let soft = fixture(true);
    let supported = bind(&soft, 0.4).unwrap();
    let expected_slope = rotated([[-0.01, 0.0], [0.0, 0.0]], 0.4);
    for temperature in [280.0, 300.0, 320.0] {
        let slopes = supported.thermal_slopes_at_temperature(temperature).unwrap();
        for i in 0..2 { for j in 0..2 {
            close(slopes.equilibrium_n_m_k[i][j], expected_slope[i][j], 2e-15);
            close(slopes.relaxation_n_m_k[0][i][j], 0.0, 0.0);
        }}
    }
    assert!(supported.at_temperature(321.0).is_err());
    assert!(supported.thermal_slopes_at_temperature(321.0).is_err());
    let f = fixture(false);
    let material = bind(&f, 0.4).unwrap();
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
