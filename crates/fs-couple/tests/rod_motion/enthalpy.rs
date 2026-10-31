//! Manufactured, independently differentiated section/finite-cell references.
//! These data are numerical fixtures, not measured material calibrations.
use super::*;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::rod_motion::enthalpy::{
    RodLayerCell, RodThermoelasticSolve, RodThermoelasticSystem,
};
use fs_couple::rod_motion::{PreparedRodLoads, RodPointForce};
use fs_matdb::SelectionPolicy;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{QuantitySpec, Temperature};
use fs_solid::rod::material_field::thermal::source::{
    RodLayerSource, RodMaterialSource, RodSegmentSource,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const LENGTH: f64 = 0.8;
const RADIUS: f64 = 0.05;
const RHO: f64 = 3.0;
const YOUNG: f64 = 1200.0;
const SHEAR: f64 = 480.0;

#[derive(Clone, Copy)]
struct Spec {
    young_scale: f64,
    young_slope: f64,
    shear_scale: f64,
    shear_slope: f64,
    alpha: f64,
    cp: f64,
    conductivity: f64,
    mechanical_upper_k: f64,
    kinked_shear: bool,
}
impl Default for Spec {
    fn default() -> Self {
        Self {
            young_scale: 1.0,
            young_slope: -1.5,
            shear_scale: 1.0,
            shear_slope: -0.4,
            alpha: 3e-4,
            cp: 0.1,
            conductivity: 0.01,
            mechanical_upper_k: 500.0,
            kinked_shear: false,
        }
    }
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn material(index: usize, spec: Spec) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, value, slope, dims) in [
        ("density", RHO, 0.0, Density::DIMS),
        (
            "young_modulus",
            YOUNG * spec.young_scale,
            spec.young_slope,
            Pressure::DIMS,
        ),
        (
            "shear_modulus",
            SHEAR * spec.shear_scale,
            spec.shear_slope,
            Pressure::DIMS,
        ),
        (
            "linear_thermal_expansion_coefficient",
            spec.alpha,
            0.0,
            INVERSE_TEMPERATURE_DIMS,
        ),
        (
            "specific_heat_capacity",
            spec.cp,
            0.0,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        ),
        (
            "thermal_conductivity",
            spec.conductivity,
            0.0,
            THERMAL_CONDUCTIVITY_DIMS,
        ),
    ] {
        let upper = if matches!(
            name,
            "young_modulus" | "shear_modulus" | "linear_thermal_expansion_coefficient"
        ) {
            spec.mechanical_upper_k
        } else {
            500.0
        };
        let curve = slope != 0.0 || (name == "shear_modulus" && spec.kinked_shear);
        let mut knots = vec![
            (200.0, value - 100.0 * slope),
            (upper, value + (upper - 300.0) * slope),
        ];
        if name == "shear_modulus" && spec.kinked_shear {
            knots.insert(1, (300.0, value + 1.0));
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curve {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots,
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    upper,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured reciprocal rod calorimetry".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("manufactured rod layer {index}"),
            phase: "solid".into(),
            process: "finite-cell thermoelastic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-7,
    )
    .unwrap();
    (card, curve)
}
fn geometry(layer: usize, layers: usize) -> (f64, f64) {
    let r0 = RADIUS * layer as f64 / layers as f64;
    let r1 = RADIUS * (layer + 1) as f64 / layers as f64;
    (
        core::f64::consts::PI * (r1 * r1 - r0 * r0),
        core::f64::consts::PI * (r1.powi(4) - r0.powi(4)) / 4.0,
    )
}
fn destinations(layers: usize, segments: usize) -> Vec<RodLayerCell> {
    (0..segments)
        .flat_map(|segment| {
            (0..layers).map(move |layer| RodLayerCell {
                segment,
                layer,
                cell: CellId((segment * layers + layer) as u64),
            })
        })
        .collect()
}
fn source<'a>(cards: &[&'a MaterialCard], segments: usize) -> RodMaterialSource<'a> {
    RodMaterialSource::try_new(
        LENGTH,
        (0..segments)
            .map(|_| RodSegmentSource {
                inner_radius_m: 0.0,
                shear_correction: 0.8,
                layers: cards
                    .iter()
                    .enumerate()
                    .map(|(i, card)| RodLayerSource {
                        outer_radius_m: RADIUS * (i + 1) as f64 / cards.len() as f64,
                        card,
                        reference_point: point(300.0),
                        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
                        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
                    })
                    .collect(),
            })
            .collect(),
        "shear_modulus",
        0.1,
        &CancelGate::new(),
    )
    .unwrap()
}
fn fixture<R>(
    specs: &[Spec],
    segments: usize,
    temperatures: &[f64],
    powers: &[f64],
    conduction: bool,
    run: impl FnOnce(&RodMaterialSource<'_>, DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R,
) -> R {
    let sources: Vec<_> = specs
        .iter()
        .enumerate()
        .map(|(i, s)| material(i, *s))
        .collect();
    let cards: Vec<_> = sources.iter().map(|s| &s.0).collect();
    let source = source(&cards, segments);
    let context = point(300.0);
    let routes = destinations(specs.len(), segments);
    assert_eq!(temperatures.len(), routes.len());
    assert_eq!(powers.len(), routes.len());
    let l = LENGTH / segments as f64;
    let mut events: Vec<_> = routes
        .iter()
        .enumerate()
        .map(|(i, route)| {
            let (a, _) = geometry(route.layer, specs.len());
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: route.cell,
                    segment_id: route.segment as u64,
                    body_id: 0,
                    material_id: MaterialId(route.layer as u64),
                    mass_kg: RHO * a * l,
                    specific_h_j_kg: (temperatures[i] - 200.0) * specs[route.layer].cp,
                    geometry: CellGeometry {
                        start_m: [route.segment as f64 * l, 0.0, 0.0],
                        end_m: [(route.segment + 1) as f64 * l, 0.0, 0.0],
                        volume_m3: a * l,
                        surface_area_m2: 2.0 * a + 2.0 * core::f64::consts::PI * RADIUS * l,
                    },
                    boundary: Boundary {
                        ambient_k: temperatures[i],
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: powers[i],
                    },
                }),
            }
        })
        .collect();
    if conduction {
        assert_eq!(routes.len(), 2);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 17,
                a: routes[0].cell,
                b: routes[1].cell,
                area_m2: 0.01,
                conductance_w_k: 0.0,
            }),
        });
    }
    let laws: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, s)| (MaterialId(i as u64), s.1.curve() as &dyn ThermalLaw))
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: routes.len(),
            max_links: 1,
            max_events: routes.len() + 1,
            max_substeps: 1024,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells: Vec<_> = routes
        .iter()
        .map(|r| MaterialCell {
            cell: r.cell,
            material_id: MaterialId(r.layer as u64),
            curve: sources[r.layer].1.curve(),
            point: ThermalMaterialPoint {
                material: &sources[r.layer].0,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: conduction.then_some((&key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect();
    let links = if conduction {
        vec![ConductionLink {
            link: 17,
            distance_a_m: 0.02,
            distance_b_m: 0.03,
        }]
    } else {
        vec![]
    };
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &links,
    };
    run(&source, thermal, &network)
}
fn solve() -> RodThermoelasticSolve {
    RodThermoelasticSolve {
        maximum_cells: 8,
        maximum_material_trials: 256,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-8,
        cell_energy_tolerance_j: 2e-11,
        thermal_energy_tolerance_j: 2e-11,
        coupling_energy_tolerance_j: 5e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-13,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn axial(
    source: &RodMaterialSource<'_>,
    temperatures: &[f64],
    strain: f64,
    clamped: bool,
) -> RodMotion {
    let gate = CancelGate::new();
    let rod = source
        .resolve_at(temperatures, &gate)
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [
                    clamped.then_some(LENGTH * (1.0 + strain)),
                    Some(0.0),
                    Some(0.0),
                ],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        rod.into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] = LENGTH * (1.0 + strain);
    RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, &gate).unwrap()
}
fn total_heat(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum()
}
fn independent_potential(
    spec: Spec,
    a: f64,
    i: f64,
    t: f64,
    gamma: [f64; 3],
    curvature: [f64; 3],
) -> (f64, f64, f64, f64) {
    let e = gamma[0] - spec.alpha * (t - 300.0);
    let young = YOUNG * spec.young_scale + spec.young_slope * (t - 300.0);
    let shear = SHEAR * spec.shear_scale + spec.shear_slope * (t - 300.0);
    let bend = i * (curvature[1] * curvature[1] + curvature[2] * curvature[2]);
    let transverse = 0.8 * a * (gamma[1] * gamma[1] + gamma[2] * gamma[2])
        + 2.0 * i * curvature[0] * curvature[0];
    let psi = 0.5 * (young * (a * e * e + bend) + shear * transverse);
    let entropy = a * (young * spec.alpha * e - 0.5 * spec.young_slope * e * e)
        - 0.5 * (spec.young_slope * bend + spec.shear_slope * transverse);
    let capacity =
        -t * a * (young * spec.alpha * spec.alpha - 2.0 * spec.young_slope * spec.alpha * e);
    (psi, entropy, psi + t * entropy, capacity)
}

#[test]
fn complete_source_energy_entropy_and_capacity_match_independent_six_strain_derivatives() {
    let specs = [
        Spec::default(),
        Spec {
            young_scale: 1.8,
            shear_scale: 1.5,
            young_slope: 0.6,
            shear_slope: 0.3,
            alpha: -1e-4,
            ..Spec::default()
        },
    ];
    let materials = [material(0, specs[0]), material(1, specs[1])];
    let source = source(&[&materials[0].0, &materials[1].0], 1);
    let gate = CancelGate::new();
    let laws = source.thermoelastic_laws(&gate).unwrap();
    let gamma = [0.012, 0.007, -0.005];
    let k = [0.08, -0.06, 0.04];
    let temperatures = [320.0, 275.0];
    let state = laws[0].at_strain(&temperatures, gamma, k, 0.1).unwrap();
    for layer in 0..2 {
        let (a, i) = geometry(layer, 2);
        let expected = independent_potential(specs[layer], a, i, temperatures[layer], gamma, k);
        let actual = state.layers[layer];
        near(actual.free_energy_j_m, expected.0, 2e-15);
        near(actual.entropy_j_m_k, expected.1, 2e-17);
        near(actual.internal_energy_j_m, expected.2, 2e-15);
        near(actual.heat_capacity_j_m_k, expected.3, 2e-17);
        let h = 1e-3;
        let mut lo = temperatures;
        let mut hi = temperatures;
        lo[layer] -= h;
        hi[layer] += h;
        let a = laws[0].at_strain(&lo, gamma, k, 0.1).unwrap().layers[layer];
        let b = laws[0].at_strain(&hi, gamma, k, 0.1).unwrap().layers[layer];
        near(
            -(b.free_energy_j_m - a.free_energy_j_m) / (2.0 * h),
            actual.entropy_j_m_k,
            2e-13,
        );
        near(
            (b.internal_energy_j_m - a.internal_energy_j_m) / (2.0 * h),
            actual.heat_capacity_j_m_k,
            2e-13,
        );
    }
    // The full source potential equals the original native section energy at
    // actual nonzero axial, transverse, twist and bending strains.
    let model = source
        .resolve_at(&temperatures, &gate)
        .unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate)
        .unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1] = [LENGTH * 1.01, 0.004, -0.003];
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(0.02, -0.015, 0.01))).unwrap();
    let strain = model.segment_strain_states(&pose, &gate).unwrap()[0];
    assert!(
        strain
            .gamma
            .iter()
            .chain(&strain.curvature)
            .all(|x| x.abs() > 1e-5)
    );
    let native = laws[0]
        .at_strain(&temperatures, strain.gamma, strain.curvature, 0.1)
        .unwrap();
    near(
        LENGTH * native.layers.iter().map(|s| s.free_energy_j_m).sum::<f64>(),
        strain.strain_energy_j,
        2e-14,
    );
    assert!(laws[0].at_strain(&[501.0, 300.0], gamma, k, 0.1).is_err());
    assert!(
        laws[0]
            .at_strain(&temperatures, gamma, [3.0, 0.0, 0.0], 0.1)
            .is_err()
    );
    let kink = material(
        2,
        Spec {
            kinked_shear: true,
            ..Spec::default()
        },
    );
    assert!(
        self::source(&[&kink.0], 1)
            .thermoelastic_laws(&gate)
            .is_err()
    );
}

#[test]
fn clamped_finite_heat_solves_internal_energy_and_refreshes_source_stress() {
    let gate = CancelGate::new();
    let dt = 0.01;
    let strain = 0.01;
    let power = 2.0;
    let mut endpoints = Vec::new();
    for spec in [
        Spec::default(),
        Spec {
            young_scale: 2.0,
            alpha: 1e-4,
            ..Spec::default()
        },
    ] {
        fixture(
            &[spec],
            1,
            &[300.0],
            &[power],
            false,
            |source, mut thermal, network| {
                let motion = axial(source, &[300.0], strain, true);
                let original_pose = motion.pose().clone();
                let inertia = motion.model().inertia().to_vec();
                let mut owner = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &thermal,
                    network,
                    &destinations(1, 1),
                    1e-11,
                    &gate,
                )
                .unwrap();
                let (a, i) = geometry(0, 1);
                let c = RHO * a * LENGTH * spec.cp;
                let before = independent_potential(spec, a, i, 300.0, [strain, 0.0, 0.0], [0.0; 3]);
                let residual = |t: f64| {
                    c * (t - 300.0)
                        + LENGTH
                            * (independent_potential(spec, a, i, t, [strain, 0.0, 0.0], [0.0; 3]).2
                                - before.2)
                        - power * dt
                };
                let (mut lo, mut hi) = (300.0, 400.0);
                for _ in 0..70 {
                    let mid = 0.5 * (lo + hi);
                    if residual(mid) > 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                let expected = 0.5 * (lo + hi);
                let report = owner
                    .step(
                        &mut thermal,
                        dt,
                        &[RodNodalWrench::default(); 2],
                        controls(dt),
                        solve(),
                        &[power * dt],
                        &gate,
                    )
                    .unwrap();
                let t = thermal.temperature_k(CellId(0)).unwrap();
                near(t, expected, 2e-8);
                let actual = report.layers[0];
                near(actual.frozen_strain_work_j, 0.0, 1e-14);
                near(
                    actual.reversible_heat_j,
                    c * (t - 300.0) - power * dt,
                    2e-11,
                );
                assert!(
                    (actual.reversible_heat_j + actual.coefficient_free_energy_change_j).abs()
                        > 1e-5
                );
                near(report.total_energy_residual_j, 0.0, 5e-11);
                near(report.entropy.accounting_residual_j_k, 0.0, 2e-13);
                assert!(report.entropy.has_unclassified_external_sources);
                assert_eq!(owner.mechanical().pose(), &original_pose);
                assert_eq!(owner.mechanical().model().inertia(), inertia);
                assert_eq!(owner.mechanical().accepted_steps(), 1);
                assert_eq!(owner.mechanical().time_s(), thermal.time_s());
                let model = owner.mechanical().model();
                let force = model.evaluate(owner.mechanical().pose(), &gate).unwrap();
                let e = YOUNG * spec.young_scale + spec.young_slope * (t - 300.0);
                near(
                    force.force_n[1][0],
                    -e * a * (strain - spec.alpha * (t - 300.0)),
                    1e-11,
                );
                endpoints.push(t);
            },
        );
    }
    assert!((endpoints[0] - endpoints[1]).abs() > 0.2);
}

#[test]
fn native_axial_motion_and_finite_heat_match_independent_split_without_numerical_heating() {
    let spec = Spec {
        alpha: 0.0,
        cp: 0.003,
        ..Spec::default()
    };
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.001;
            let motion = axial(source, &[300.0], 0.01, false);
            let mut owner = RodThermoelasticSystem::try_new(
                motion,
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let (a, _) = geometry(0, 1);
            let mass = 0.5 * RHO * a * LENGTH;
            let capacity = RHO * a * LENGTH * spec.cp;
            let kprime = spec.young_slope * a / LENGTH;
            let mut q = 0.008;
            let mut v = 0.0;
            let mut t = 300.0;
            let total0 = owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal);
            let mut native_defect = 0.0;
            let mut heat_sum = 0.0;
            for _ in 0..40 {
                let k = (YOUNG + spec.young_slope * (t - 300.0)) * a / LENGTH;
                let q1 = q + dt * v - 0.5 * dt * dt * k * q / mass;
                let v1 = v - 0.5 * dt * k * (q + q1) / mass;
                let reversible = 0.5 * t * kprime * (q1 * q1 - q * q);
                let t1 = t + reversible / capacity;
                let target = owner.mechanical().time_s() + dt;
                let report = owner
                    .step(
                        &mut thermal,
                        target,
                        &[RodNodalWrench::default(); 2],
                        controls(dt * 1.001),
                        solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                near(
                    owner.mechanical().pose().positions_m[1][0] - LENGTH,
                    q1,
                    2e-12,
                );
                near(
                    owner.mechanical().velocities_m_s().unwrap()[1][0],
                    v1,
                    2e-11,
                );
                near(thermal.temperature_k(CellId(0)).unwrap(), t1, 2e-7);
                near(report.reversible_heat_j, reversible, 2e-12);
                near(report.external_thermal_heat_j, 0.0, 2e-12);
                near(report.coupling_energy_residual_j, 0.0, 2e-11);
                near(report.entropy.accounting_residual_j_k, 0.0, 2e-13);
                assert!(!report.entropy.has_unclassified_external_sources);
                native_defect += report.mechanical.integration_defect_j;
                heat_sum += report.reversible_heat_j;
                q = q1;
                v = v1;
                t = t1;
            }
            near(
                owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal) - total0,
                native_defect,
                5e-10,
            );
            near(
                capacity * (thermal.temperature_k(CellId(0)).unwrap() - 300.0),
                heat_sum,
                5e-10,
            );
            assert!(native_defect.abs() > 1e-8);
            assert!(heat_sum.abs() > 100.0 * native_defect.abs());
            assert!(t > 301.0);
        },
    );
}

#[test]
fn prepared_world_loads_preserve_native_motion_and_finite_heat() {
    fixture(
        &[Spec::default()],
        1,
        &[300.0],
        &[0.5],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.001;
            let motion = axial(source, &[300.0], 0.01, false);
            let mut other_thermal = thermal.clone();
            let routes = destinations(1, 1);
            let mut world = RodThermoelasticSystem::try_new(
                motion.clone(), source, &thermal, network, &routes, 1e-11, &gate,
            )
            .unwrap();
            let mut spatial = RodThermoelasticSystem::try_new(
                motion, source, &other_thermal, network, &routes, 1e-11, &gate,
            )
            .unwrap();
            let loads = [
                RodNodalWrench::default(),
                RodNodalWrench {
                    force_n: [0.03, 0.0, 0.0],
                    torque_nm: [0.002, 0.0, 0.0],
                },
            ];
            let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
            let a = world
                .step(&mut thermal, dt, &loads, controls(dt), solve(), &[0.0], &gate)
                .unwrap();
            let b = spatial
                .step_spatial(
                    &mut other_thermal, dt, &prepared, controls(dt), solve(), &[0.0], &gate,
                )
                .unwrap();
            assert_eq!(world.mechanical().pose(), spatial.mechanical().pose());
            assert_eq!(
                world.mechanical().linear_momenta_kg_m_s(),
                spatial.mechanical().linear_momenta_kg_m_s()
            );
            assert_eq!(
                world.mechanical().body_spins_kg_m2_s(),
                spatial.mechanical().body_spins_kg_m2_s()
            );
            assert_eq!(world.layers(), spatial.layers());
            assert_eq!(thermal.cells(), other_thermal.cells());
            assert_eq!(
                a.mechanical.source_work_quadrature_j,
                b.mechanical.source_work_quadrature_j
            );
            assert_eq!(a.mechanical.support_impulses, b.mechanical.support_impulses);
            assert_eq!(a.total_energy_residual_j, b.total_energy_residual_j);
        },
    );
}

#[test]
fn rotating_material_loads_and_finite_heat_keep_physical_work_and_atomic_retry() {
    let spec = Spec {
        alpha: 0.0,
        mechanical_upper_k: 300.5,
        ..Spec::default()
    };
    fixture(&[spec], 1, &[300.0], &[0.5], false, |source, mut thermal, network| {
        let gate = CancelGate::new();
        let model = Arc::new(
            source.resolve_at(&[300.0], &gate).unwrap()
                .into_motion_model(RodRootBoundary::Free, &gate).unwrap(),
        );
        let mut pose = model.initial_pose();
        pose.frames
            .fill(So3::exp(So3Tangent::new(Vec3::new(0.4, 0.0, 0.0))).unwrap());
        let original = RodMotion::new(
            model, pose, &[[0.0; 3]; 2], &[[2.0, 0.0, 0.0]; 2], 0.0, &gate,
        )
        .unwrap();
        let original_thermal = thermal.clone();
        let mut owner = RodThermoelasticSystem::try_new(
            original.clone(), source, &thermal, network, &destinations(1, 1), 1e-11, &gate,
        )
        .unwrap();
        let original_layers = owner.layers().to_vec();
        let offset = [0.0, 0.04, 0.0];
        let force = [0.0, 0.0, 0.02];
        let couple = [0.0, 0.0002, 0.0];
        let prepared = PreparedRodLoads::new(
            &[RodNodalWrench::default(); 2],
            &[couple; 2],
            &[
                RodPointForce { node: 0, offset_body_m: offset, force_world_n: force },
                RodPointForce { node: 1, offset_body_m: offset, force_world_n: force },
            ],
            &gate,
        )
        .unwrap();
        // Mechanics admits this trial. The later finite heat leaves the source
        // interval; neither physical history nor its momenta may be published.
        original
            .step_spatial_candidate(0.01, &prepared, controls(0.01), &gate)
            .unwrap();
        assert!(
            owner.step_spatial(
                &mut thermal, 0.01, &prepared, controls(0.01), solve(), &[0.0], &gate,
            ).is_err()
        );
        assert_eq!(thermal.cells(), original_thermal.cells());
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(owner.mechanical().pose(), original.pose());
        assert_eq!(
            owner.mechanical().linear_momenta_kg_m_s(), original.linear_momenta_kg_m_s()
        );
        assert_eq!(
            owner.mechanical().body_spins_kg_m2_s(), original.body_spins_kg_m2_s()
        );
        assert_eq!(owner.mechanical().time_s(), 0.0);
        assert_eq!(owner.mechanical().accepted_steps(), 0);
        assert_eq!(owner.layers(), original_layers);

        let dt = 0.001;
        let (frozen, native) = original
            .step_spatial_candidate(dt, &prepared, controls(dt), &gate).unwrap();
        let report = owner.step_spatial(
            &mut thermal, dt, &prepared, controls(dt), solve(), &[0.0], &gate,
        ).unwrap();
        assert_eq!(owner.mechanical().pose(), frozen.pose());
        assert_eq!(
            owner.mechanical().linear_momenta_kg_m_s(), frozen.linear_momenta_kg_m_s()
        );
        assert_eq!(
            owner.mechanical().body_spins_kg_m2_s(), frozen.body_spins_kg_m2_s()
        );
        assert_eq!(
            report.mechanical.source_work_quadrature_j, native.source_work_quadrature_j
        );
        assert_eq!(
            report.mechanical.integration_defect_j, native.integration_defect_j
        );

        // Reconstruct torque about the fixed world origin from physical point
        // positions and body couples, independently of the prepared reducer.
        let vector = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
        let total_moment = |state: &RodMotion| {
            let mut torque = Vec3::default();
            for node in 0..2 {
                let frame = state.pose().frames[node];
                let point = vector(state.pose().positions_m[node])
                    + frame.rotate(vector(offset)).unwrap();
                torque = torque + point.cross(vector(force))
                    + frame.rotate(vector(couple)).unwrap();
            }
            torque
        };
        let expected = (total_moment(&original) + total_moment(owner.mechanical())).scale(0.5 * dt);
        let stale = total_moment(&original).scale(dt);
        assert!((expected - stale).norm() > 1e-10);
        let expected = [expected.x, expected.y, expected.z];
        for c in 0..3 {
            near(report.mechanical.external_angular_impulse_nm_s[c], expected[c], 2e-12);
            near(
                report.mechanical.after.angular_momentum_kg_m2_s[c]
                    - report.mechanical.before.angular_momentum_kg_m2_s[c],
                expected[c], 2e-12,
            );
            near(
                report.mechanical.after.linear_momentum_kg_m_s[c]
                    - report.mechanical.before.linear_momentum_kg_m_s[c],
                2.0 * force[c] * dt, 2e-12,
            );
        }

        let temperature = thermal.temperature_k(CellId(0)).unwrap();
        assert!((300.0..300.5).contains(&temperature));
        let (area, second_moment) = geometry(0, 1);
        let strain = owner.mechanical().model()
            .segment_strain_states(owner.mechanical().pose(), &gate).unwrap()[0];
        let material_energy = LENGTH * independent_potential(
            spec, area, second_moment, temperature, strain.gamma, strain.curvature,
        ).2;
        let mut kinetic = 0.0;
        for node in 0..2 {
            let inertia = &owner.mechanical().model().inertia()[node];
            for c in 0..3 {
                let linear = owner.mechanical().linear_momenta_kg_m_s()[node][c];
                let spin = owner.mechanical().body_spins_kg_m2_s()[node][c];
                kinetic += 0.5 * linear * linear / inertia.mass_kg()
                    + 0.5 * spin * spin / inertia.principal_kg_m2()[c];
            }
        }
        // Initial strain/entropy vanish. Reconstruct complete actual storage;
        // retain the native defect separately rather than injecting it as heat.
        let change = kinetic + material_energy + total_heat(&thermal)
            - native.before.kinetic_energy_j - total_heat(&original_thermal);
        near(
            change - 0.5 * dt - native.source_work_quadrature_j,
            native.integration_defect_j, 5e-10,
        );
        near(
            RHO * area * LENGTH * spec.cp * (temperature - 300.0) - 0.5 * dt,
            report.reversible_heat_j, 2e-11,
        );
        near(report.coupling_energy_residual_j, 0.0, 5e-10);
    });
}

#[test]
fn nonuniform_layer_and_segment_conduction_changes_native_torsion_and_closes_finite_heat() {
    let gate = CancelGate::new();
    let old = [310.0, 350.0];
    let dt = 0.02;
    let axial_strain = 0.012;
    let twist = 0.05;
    for (layers, segments) in [(2, 1), (1, 2)] {
        let specs: Vec<_> = (0..layers)
            .map(|layer| Spec {
                shear_scale: 1.0 + layer as f64,
                shear_slope: -0.4 - 0.2 * layer as f64,
                cp: 0.1 + 0.05 * layer as f64,
                ..Spec::default()
            })
            .collect();
        fixture(
            &specs,
            segments,
            &old,
            &[0.0; 2],
            true,
            |source, mut thermal, network| {
                let l = LENGTH / segments as f64;
                let supports: Vec<_> = (1..=segments)
                    .map(|node| RodSupport {
                        node,
                        position_m: [
                            Some(node as f64 * l * (1.0 + axial_strain)),
                            Some(0.0),
                            Some(0.0),
                        ],
                        frame: Some(
                            So3::exp(So3Tangent::new(Vec3::new(
                                node as f64 * l * twist,
                                0.0,
                                0.0,
                            )))
                            .unwrap(),
                        ),
                    })
                    .collect();
                let rod = source
                    .resolve_at(&old, &gate)
                    .unwrap()
                    .with_supports(&supports, &gate)
                    .unwrap();
                let model = Arc::new(
                    rod.into_motion_model(RodRootBoundary::Clamped, &gate)
                        .unwrap(),
                );
                let pose = model.initial_pose();
                let zeros = vec![[0.0; 3]; segments + 1];
                let motion =
                    RodMotion::new(model, pose.clone(), &zeros, &zeros, 0.0, &gate).unwrap();
                let mut owner = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &thermal,
                    network,
                    &destinations(layers, segments),
                    1e-11,
                    &gate,
                )
                .unwrap();
                let source_spec = [specs[0], specs[1 % layers]];
                let geometry_at = |index: usize| geometry(index % layers, layers);
                let c: [f64; 2] = core::array::from_fn(|index| {
                    RHO * geometry_at(index).0 * l * source_spec[index].cp
                });
                let potential = |index: usize, t: f64| {
                    let (a, i) = geometry_at(index);
                    independent_potential(
                        source_spec[index],
                        a,
                        i,
                        t,
                        [axial_strain, 0.0, 0.0],
                        [twist, 0.0, 0.0],
                    )
                };
                let u0 = [l * potential(0, old[0]).2, l * potential(1, old[1]).2];
                let conductance = 0.01
                    / (0.02 / source_spec[0].conductivity + 0.03 / source_spec[1].conductivity);
                let mut expected = old;
                // Independent two-cell implicit balance of H + Psi + TS.
                for _ in 0..16 {
                    let p = [potential(0, expected[0]), potential(1, expected[1])];
                    let q = dt * conductance * (expected[1] - expected[0]);
                    let residual = [
                        c[0] * (expected[0] - old[0]) + l * p[0].2 - u0[0] - q,
                        c[1] * (expected[1] - old[1]) + l * p[1].2 - u0[1] + q,
                    ];
                    let diagonal = [
                        c[0] + l * p[0].3 + dt * conductance,
                        c[1] + l * p[1].3 + dt * conductance,
                    ];
                    let off = -dt * conductance;
                    let determinant = diagonal[0] * diagonal[1] - off * off;
                    expected[0] -= (diagonal[1] * residual[0] - off * residual[1]) / determinant;
                    expected[1] -= (diagonal[0] * residual[1] - off * residual[0]) / determinant;
                }
                let total0 = owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal);
                let report = owner
                    .step(
                        &mut thermal,
                        dt,
                        &vec![RodNodalWrench::default(); segments + 1],
                        controls(dt),
                        solve(),
                        &[0.0; 2],
                        &gate,
                    )
                    .unwrap();
                for index in 0..2 {
                    near(
                        thermal.temperature_k(CellId(index as u64)).unwrap(),
                        expected[index],
                        2e-7,
                    );
                    near(
                        report.layers[index].states[1].internal_energy_j,
                        l * potential(index, expected[index]).2,
                        2e-10,
                    );
                }
                assert!(expected[0] > old[0] && expected[1] < old[1]);
                near(
                    owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal),
                    total0,
                    5e-10,
                );
                near(report.external_thermal_heat_j, 0.0, 3e-11);
                near(report.total_energy_residual_j, 0.0, 3e-11);
                assert!(report.entropy.conduction_generation_j_k > 0.0);
                near(report.entropy.accounting_residual_j_k, 0.0, 2e-13);
                assert_eq!(owner.mechanical().pose(), &pose);
                let expected_gj: f64 = (0..layers)
                    .map(|layer| {
                        let index = (segments - 1) * layers + layer;
                        let spec = specs[layer];
                        2.0 * geometry(layer, layers).1
                            * (SHEAR * spec.shear_scale
                                + spec.shear_slope * (expected[index] - 300.0))
                    })
                    .sum();
                let response = owner
                    .mechanical()
                    .model()
                    .evaluate(owner.mechanical().pose(), &gate)
                    .unwrap();
                near(
                    response.torque_body_nm[segments][0],
                    -expected_gj * twist,
                    2e-12,
                );
            },
        );
    }
}

fn continuous_axial_reference(spec: Spec, duration: f64) -> [f64; 3] {
    let a = geometry(0, 1).0;
    let mass = 0.5 * RHO * a * LENGTH;
    let capacity = RHO * a * LENGTH * spec.cp;
    let derivative = |y: [f64; 3]| {
        let [q, v, t] = y;
        [
            v,
            -(YOUNG * spec.young_scale + spec.young_slope * (t - 300.0)) * a / LENGTH * q / mass,
            t * spec.young_slope * a / LENGTH * q * v / capacity,
        ]
    };
    let mut y = [0.008, 0.0, 300.0];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        y = core::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn reciprocal_split_converges_at_first_order_to_continuous_thermodynamics() {
    let spec = Spec {
        alpha: 0.0,
        young_slope: -3.0,
        cp: 0.0003,
        ..Spec::default()
    };
    let reference = continuous_axial_reference(spec, 0.04);
    let mut errors = Vec::new();
    for count in [40, 80, 160] {
        fixture(
            &[spec],
            1,
            &[300.0],
            &[0.0],
            false,
            |source, mut thermal, network| {
                let gate = CancelGate::new();
                let dt = 0.04 / count as f64;
                let motion = axial(source, &[300.0], 0.01, false);
                let mut owner = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &thermal,
                    network,
                    &destinations(1, 1),
                    1e-11,
                    &gate,
                )
                .unwrap();
                for _ in 0..count {
                    let target = owner.mechanical().time_s() + dt;
                    owner
                        .step(
                            &mut thermal,
                            target,
                            &[RodNodalWrench::default(); 2],
                            controls(dt * 1.001),
                            solve(),
                            &[0.0],
                            &gate,
                        )
                        .unwrap();
                }
                let actual = [
                    owner.mechanical().pose().positions_m[1][0] - LENGTH,
                    owner.mechanical().velocities_m_s().unwrap()[1][0],
                    thermal.temperature_k(CellId(0)).unwrap(),
                ];
                errors.push(
                    (actual[0] - reference[0]).abs() / 0.008
                        + (actual[1] - reference[1]).abs() / 0.3
                        + (actual[2] - reference[2]).abs() / 300.0,
                );
            },
        );
    }
    for pair in errors.windows(2) {
        assert!((1.85..2.15).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn source_capacity_and_late_domain_refusals_leave_both_original_histories_unchanged() {
    let gate = CancelGate::new();
    let spec = Spec {
        mechanical_upper_k: 305.0,
        ..Spec::default()
    };
    fixture(
        &[spec],
        1,
        &[300.0],
        &[2.0],
        false,
        |source, mut thermal, network| {
            let motion = axial(source, &[300.0], 0.01, true);
            let original = motion.clone();
            let original_thermal = thermal.clone();
            let mut owner = RodThermoelasticSystem::try_new(
                motion,
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let original_layers = owner.layers().to_vec();
            // The source is valid initially. Only the attempted finite heating
            // endpoint leaves E/G/alpha support; no accepted owner may advance.
            assert!(
                owner
                    .step(
                        &mut thermal,
                        0.01,
                        &[RodNodalWrench::default(); 2],
                        controls(0.01),
                        solve(),
                        &[0.0],
                        &gate
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), original_thermal.cells());
            assert_eq!(thermal.time_s(), original_thermal.time_s());
            assert_eq!(owner.mechanical().pose(), original.pose());
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                original.linear_momenta_kg_m_s()
            );
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                original.body_spins_kg_m2_s()
            );
            assert_eq!(owner.mechanical().accepted_steps(), 0);
            assert_eq!(owner.mechanical().time_s(), 0.0);
            assert_eq!(owner.layers(), original_layers);
            let cancelled = CancelGate::new();
            cancelled.request();
            assert!(
                owner
                    .step(
                        &mut thermal,
                        0.001,
                        &[RodNodalWrench::default(); 2],
                        controls(0.001),
                        solve(),
                        &[0.0],
                        &cancelled
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), original_thermal.cells());
            assert_eq!(owner.layers(), original_layers);
            let mut fresh_thermal = original_thermal;
            let mut fresh = RodThermoelasticSystem::try_new(
                original,
                source,
                &fresh_thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let retry = owner
                .step(
                    &mut thermal,
                    0.001,
                    &[RodNodalWrench::default(); 2],
                    controls(0.001),
                    solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let uninterrupted = fresh
                .step(
                    &mut fresh_thermal,
                    0.001,
                    &[RodNodalWrench::default(); 2],
                    controls(0.001),
                    solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert_eq!(thermal.cells(), fresh_thermal.cells());
            assert_eq!(owner.mechanical().pose(), fresh.mechanical().pose());
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                fresh.mechanical().linear_momenta_kg_m_s()
            );
            assert_eq!(owner.layers(), fresh.layers());
            assert_eq!(
                retry.total_energy_residual_j,
                uninterrupted.total_energy_residual_j
            );
        },
    );
    fixture(
        &[Spec {
            cp: 1e-7,
            ..Spec::default()
        }],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, thermal, network| {
            let motion = axial(source, &[300.0], 0.01, true);
            assert!(
                RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &thermal,
                    network,
                    &destinations(1, 1),
                    1e-11,
                    &gate
                )
                .is_err()
            );
        },
    );
    fixture(
        &[Spec::default(); 2],
        1,
        &[300.0; 2],
        &[0.0; 2],
        false,
        |source, thermal, network| {
            let motion = axial(source, &[300.0; 2], 0.01, true);
            let mut routes = destinations(2, 1);
            routes[1].cell = routes[0].cell;
            assert!(
                RodThermoelasticSystem::try_new(
                    motion, source, &thermal, network, &routes, 1e-11, &gate
                )
                .is_err()
            );
        },
    );
}
