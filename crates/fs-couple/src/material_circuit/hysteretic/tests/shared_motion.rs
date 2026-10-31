//! One sourced core memory drives all winding ports and the original body.
//! Scaled authored materials exercise causality/work, not measured alloy fidelity.
use super::*;
use crate::material_circuit::hysteretic::shared::{
    SharedHystereticCircuit, SharedHystereticCircuitInitial,
    motion::{SharedHystereticGapMotion, SharedHystereticGapMotionSolve},
};
use crate::material_circuit::motion::support::AxialSupport;
use fs_ga::Vec3;
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor,
        hysteresis::shared::{SharedHystereticCore, SharedHystereticWinding},
    },
    section::{
        ConcentricElasticSection, ConcentricMaterialLayer, thermal::ConcentricThermoelasticSection,
    },
    state_point::{
        INVERSE_TEMPERATURE_DIMS, ScalarAdmissibility, ScalarPropertyRequirement,
        integrate_isotropic_thermal_expansion, resolve_material_state_point,
    },
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::{Dims, Pressure};

const GAP: f64 = 0.7;
const MASS: f64 = 2.5;
const RATIOS: [f64; 2] = [2.0, -1.0];
const LEAKAGE: [f64; 2] = [0.3, 0.4];
const RESISTANCE: [f64; 2] = [0.7, 1.1];
const PORTS: [f64; 3] = [1.2, 17.0, -0.9];
const INITIAL: [f64; 2] = [0.9, -0.7];

fn card(properties: &[(&str, Dims, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored shared hysteretic motion fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored common-pole fixture".into(),
            phase: "solid".into(),
            process: "analytic test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn body(mass: f64) -> ResolvedBody {
    let card = card(&[("density", Density::DIMS, mass)]);
    resolve_body(
        &[RegionBinding {
            key: "guided common pole",
            material: &card,
            point: &QueryPoint::new(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        }],
        1,
        || false,
    )
    .unwrap()
}

fn shared(core: UniformHystereticCore) -> SharedHystereticCore {
    SharedHystereticCore::new(
        core,
        RATIOS
            .into_iter()
            .zip(LEAKAGE)
            .map(|(turns_ratio, leakage_h)| SharedHystereticWinding {
                turns_ratio,
                leakage_h,
            })
            .collect(),
    )
    .unwrap()
}

fn with_gap(core: SharedHystereticCore, limits: [f64; 2]) -> SharedHystereticCore {
    let medium = card(&[(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        0.8,
    )]);
    core.with_gap(
        ResolvedInductor::resolve(
            &medium,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.5,
            GAP,
            1,
        )
        .unwrap(),
        limits,
    )
    .unwrap()
}

fn linear_core() -> SharedHystereticCore {
    with_gap(
        shared(
            UniformHystereticCore::new(
                ScalarHysteresis::new(0.1, 30.0, vec![]).unwrap(),
                1.0,
                1.0,
                1.0,
            )
            .unwrap(),
        ),
        [0.2, 2.0],
    )
}

fn sourced(pinning: f64) -> SharedHystereticCore {
    let (card, _) = source::core_material(pinning, 310.0);
    with_gap(shared(source::sourced_core(&card)), [0.2, 2.0])
}

fn graph(resistive: bool) -> (CircuitGraph, Vec<usize>) {
    if resistive {
        (
            CircuitGraph {
                node_count: 5,
                branches: vec![
                    (1, 0, Branch::VoltageSource { port: 0 }),
                    (1, 2, Branch::Resistor { ohms: f64::NAN }),
                    (2, 0, Branch::Inductor { henries: f64::NAN }),
                    (3, 0, Branch::VoltageSource { port: 2 }),
                    (3, 4, Branch::Resistor { ohms: f64::NAN }),
                    (4, 0, Branch::Inductor { henries: f64::NAN }),
                ],
                transformers: vec![],
            },
            vec![2, 5],
        )
    } else {
        (
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (1, 0, Branch::VoltageSource { port: 0 }),
                    (1, 0, Branch::Inductor { henries: f64::NAN }),
                    (2, 0, Branch::VoltageSource { port: 2 }),
                    (2, 0, Branch::Inductor { henries: f64::NAN }),
                ],
                transformers: vec![],
            },
            vec![1, 3],
        )
    }
}

fn resistors() -> Vec<MaterialResistor> {
    RESISTANCE
        .into_iter()
        .zip([1, 4])
        .map(|(value, branch)| resistor(branch, &material(value, 0.0).0))
        .collect()
}

fn circuit(
    core: SharedHystereticCore,
    linkage: [f64; 2],
    inputs: &[f64],
    resistive: bool,
) -> SharedHystereticCircuit {
    let magnetic = core
        .candidate_linkages(&core.initial_state(), &linkage)
        .unwrap()
        .state;
    let (graph, branches) = graph(resistive);
    SharedHystereticCircuit::new(
        graph,
        branches,
        core,
        if resistive { resistors() } else { vec![] },
        SharedHystereticCircuitInitial {
            magnetic,
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: inputs.to_vec(),
            time_s: 0.0,
        },
    )
    .unwrap()
}

fn policy() -> SharedHystereticGapMotionSolve {
    let mut electrical = controls();
    electrical.newton.linear_restart = 4;
    SharedHystereticGapMotionSolve {
        circuit: electrical,
        difference_step_m: 1e-6,
        displacement_tolerance_m: 1e-9,
        energy_tolerance_j: 1e-8,
    }
}

// Independent full inductance-matrix inverse and physical ODE. This does not
// call the shared load-line update, force or path-integration implementation.
fn rate(x: [f64; 4]) -> [f64; 4] {
    let permeance = 1.0 / (10.0 + x[2] / 0.4);
    let a = LEAKAGE[0] + 4.0 * permeance;
    let d = LEAKAGE[1] + permeance;
    let mutual = -2.0 * permeance;
    let determinant = a * d - mutual * mutual;
    let currents = [
        (d * x[0] - mutual * x[1]) / determinant,
        (a * x[1] - mutual * x[0]) / determinant,
    ];
    let flux = permeance * (2.0 * currents[0] - currents[1]);
    [
        PORTS[0] - RESISTANCE[0] * currents[0],
        PORTS[2] - RESISTANCE[1] * currents[1],
        x[3] / MASS,
        0.4 - 0.5 * flux * flux / 0.4,
    ]
}

fn rk4() -> [f64; 4] {
    let mut x = [INITIAL[0], INITIAL[1], GAP, MASS * 0.25];
    let dt = 0.2 / 2000.0;
    for _ in 0..2000 {
        let k1 = rate(x);
        let k2 = rate(std::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]));
        let k3 = rate(std::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]));
        let k4 = rate(std::array::from_fn(|i| x[i] + dt * k3[i]));
        x = std::array::from_fn(|i| x[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    x
}

#[test]
fn linear_moving_limit_refines_to_the_independent_mutual_inductance_ode() {
    let reference = rk4();
    let mut errors = Vec::new();
    for ticks in [5, 10] {
        let mut model = SharedHystereticGapMotion::new(
            circuit(linear_core(), INITIAL, &PORTS, true),
            body(MASS),
            None,
            0.25,
        )
        .unwrap();
        for tick in 1..=ticks {
            let next = model
                .step_candidate(
                    0.2 * f64::from(tick) / f64::from(ticks),
                    resistors(),
                    &PORTS,
                    0.4,
                    policy(),
                    || false,
                )
                .unwrap();
            let step = next.step();
            assert!(step.magnetic.gap.as_ref().unwrap().mean_force_n < 0.0);
            near(step.energy_residual_j, 0.0, 1e-8);
            near(step.magnetic.core.core_heat_j, 0.0, 0.0);
            near(
                step.kinetic_energy_j[1] - step.kinetic_energy_j[0],
                step.mechanical_supply_j + step.field_work_on_body_j,
                1e-8,
            );
            assert_eq!(next.model().circuit().branch_indices(), &[2, 5]);
            assert_eq!(step.circuit_trials, step.mechanical_trials + 3);
            assert_eq!(
                step.circuit
                    .resistors
                    .iter()
                    .map(|r| r.material.branch)
                    .collect::<Vec<_>>(),
                vec![1, 4]
            );
            model = next.into_parts().0;
        }
        let linkages = model.circuit().magnetic_state().flux_linkages_wb();
        let actual = [
            linkages[0],
            linkages[1],
            model.gap_m(),
            model.momentum_kg_m_s(),
        ];
        errors.push(
            actual
                .iter()
                .zip(reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(errors[1] < 0.4 * errors[0], "refinement errors: {errors:?}");
    assert!(errors[1] < 3e-5, "refined ODE error: {}", errors[1]);
}

#[test]
fn zero_voltage_motion_conserves_total_linkage_and_exchanges_actual_field_work() {
    let inputs = [0.0, 17.0, 0.0];
    let model = SharedHystereticGapMotion::new(
        circuit(linear_core(), INITIAL, &inputs, false),
        body(MASS),
        None,
        0.2,
    )
    .unwrap();
    let next = model
        .step_candidate(0.05, vec![], &inputs, 0.0, policy(), || false)
        .unwrap();
    let step = next.step();
    for (&actual, expected) in next
        .model()
        .circuit()
        .magnetic_state()
        .flux_linkages_wb()
        .iter()
        .zip(INITIAL)
    {
        near(actual, expected, 3e-10);
    }
    assert!(next.model().gap_m() > GAP);
    assert!(next.model().velocity_m_s() < model.velocity_m_s());
    assert!(
        (next.model().circuit().magnetic_state().currents_a()[0]
            - model.circuit().magnetic_state().currents_a()[0])
            .abs()
            > 1e-5
    );
    near(step.electrical_supply_j, 0.0, 1e-12);
    near(step.mechanical_supply_j, 0.0, 0.0);
    near(step.heat_j, 0.0, 0.0);
    near(step.stored_change_j, 0.0, 1e-8);
    near(
        step.magnetic.stored_change_j,
        -step.field_work_on_body_j,
        1e-8,
    );
    assert!(step.magnetic.stored_change_j > 0.0);
}

#[test]
fn original_pinning_and_body_density_change_motion_with_one_common_heat() {
    let inputs = [5.0, 17.0, -3.0];
    let mut outcomes = Vec::new();
    for (pinning, mass) in [(1.0, MASS), (2.0, MASS), (1.0, 2.0 * MASS)] {
        let model = SharedHystereticGapMotion::new(
            circuit(sourced(pinning), [0.0; 2], &inputs, true),
            body(mass),
            None,
            0.2,
        )
        .unwrap();
        let before = model.clone();
        let next = model
            .step_candidate(0.1, resistors(), &inputs, 0.0, policy(), || false)
            .unwrap();
        let step = next.step();
        let polarization = step.magnetic.state.core_state().polarizations_t()[0];
        assert!(polarization > 0.0);
        near(
            step.magnetic.core.core_heat_j,
            pinning * polarization,
            2e-10,
        );
        near(
            step.heat_j,
            step.circuit.heat_j + pinning * polarization,
            2e-10,
        );
        near(
            step.stored_change_j + step.heat_j,
            step.electrical_supply_j,
            1e-8,
        );
        near(
            step.kinetic_energy_j[1] - step.kinetic_energy_j[0],
            step.field_work_on_body_j,
            1e-8,
        );
        assert!(step.field_work_on_body_j < 0.0);
        assert!(next.model().velocity_m_s() < 0.2);
        assert_eq!(model, before);
        assert_eq!(next.model().body(), model.body());
        near(next.model().generalized_mass_kg(), mass, 1e-12);
        assert!(next.model().with_core_source_query(&point(320.0)).is_err());
        let refreshed = next.model().with_core_source_query(&point(305.0)).unwrap();
        assert_eq!(
            refreshed.circuit().magnetic_state(),
            next.model().circuit().magnetic_state()
        );
        near(
            refreshed.momentum_kg_m_s(),
            next.model().momentum_kg_m_s(),
            0.0,
        );
        outcomes.push((next.model().velocity_m_s(), step.magnetic.core.core_heat_j));
    }
    assert!((outcomes[0].0 - outcomes[1].0).abs() > 1e-8);
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 1e-6);
    assert!(outcomes[2].0 > outcomes[0].0);
}

fn support() -> AxialSupport {
    let card = card(&[
        ("density", Density::DIMS, 3.0),
        ("young_modulus", Pressure::DIMS, 100.0),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            0.0,
        ),
    ]);
    let point = point(300.0);
    let material = resolve_material_state_point(
        &card,
        &point,
        &[
            ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "young_modulus",
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = integrate_isotropic_thermal_expansion(
        &card,
        &point,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    AxialSupport::try_new(
        ConcentricThermoelasticSection::try_new(
            ConcentricElasticSection::try_new(
                0.0,
                vec![ConcentricMaterialLayer {
                    outer_radius_m: (0.02 / core::f64::consts::PI).sqrt(),
                    material,
                }],
            )
            .unwrap(),
            vec![expansion],
        )
        .unwrap(),
        2.0,
        0.68,
        0.1,
    )
    .unwrap()
}

#[test]
fn original_axial_support_mass_and_elastic_energy_enter_the_same_native_motion() {
    let member = support();
    let inputs = [0.0, 17.0, 0.0];
    let model = SharedHystereticGapMotion::new(
        circuit(linear_core(), [0.0; 2], &inputs, false),
        body(MASS),
        Some(member.clone()),
        0.1,
    )
    .unwrap();
    let mass = MASS + 3.0 * 0.02 * 2.0 / 3.0;
    let stiffness = 100.0 * 0.02 / 2.0;
    let (dt, force, velocity) = (0.05, 0.03, 0.1);
    let ratio = dt * dt * stiffness / (4.0 * mass);
    let v1 =
        ((1.0 - ratio) * velocity + dt * (force - stiffness * (GAP - 0.68)) / mass) / (1.0 + ratio);
    let g1 = GAP + 0.5 * dt * (velocity + v1);
    let next = model
        .step_candidate(dt, vec![], &inputs, force, policy(), || false)
        .unwrap();
    near(model.generalized_mass_kg(), mass, 1e-12);
    near(next.model().gap_m(), g1, 3e-10);
    near(next.model().velocity_m_s(), v1, 3e-10);
    near(
        next.step().elastic_energy_j[1],
        0.5 * stiffness * (g1 - 0.68).powi(2),
        3e-10,
    );
    near(next.step().stored_change_j, force * (g1 - GAP), 1e-8);
    near(next.step().heat_j, 0.0, 0.0);
    assert_eq!(next.model().support(), Some(&member));
}

#[test]
fn missing_gap_limits_trial_budget_and_late_cancel_never_publish_partial_motion() {
    let missing = shared(core(false));
    assert!(
        SharedHystereticGapMotion::new(
            circuit(missing, [0.0; 2], &PORTS, true),
            body(MASS),
            None,
            0.0,
        )
        .is_err()
    );
    let model = SharedHystereticGapMotion::new(
        circuit(linear_core(), INITIAL, &PORTS, true),
        body(MASS),
        None,
        0.25,
    )
    .unwrap();
    let saved = model.clone();
    let mut polls = 0;
    let first = model
        .step_candidate(0.04, resistors(), &PORTS, 0.4, policy(), || {
            polls += 1;
            false
        })
        .unwrap();
    let retry = model
        .step_candidate(0.04, resistors(), &PORTS, 0.4, policy(), || false)
        .unwrap();
    assert_eq!(first.model(), retry.model());
    let mut count = 0;
    assert!(matches!(
        model.step_candidate(0.04, resistors(), &PORTS, 0.4, policy(), || {
            count += 1;
            count == polls
        }),
        Err(HystereticCircuitError::Cancelled)
    ));
    let mut limited = policy();
    limited.circuit.maximum_circuit_trials = 4;
    assert!(matches!(
        model.step_candidate(0.04, resistors(), &PORTS, 0.4, limited, || false),
        Err(HystereticCircuitError::TrialBudget)
    ));
    assert_eq!(model, saved);

    let constrained = with_gap(shared(core(false)), [0.699, 0.701]);
    let bounded = SharedHystereticGapMotion::new(
        circuit(constrained, [0.0; 2], &[0.0; 3], false),
        body(MASS),
        None,
        1.0,
    )
    .unwrap();
    let saved = bounded.clone();
    assert!(
        bounded
            .step_candidate(0.1, vec![], &[0.0; 3], 0.0, policy(), || false)
            .is_err()
    );
    assert_eq!(bounded, saved);
}
