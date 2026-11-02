//! G0/G1 momentum, retained-history relaxation and distinct physical/numerical
//! energy checks for native finite Hex8 inertia. Coefficients are synthetic.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::hex::dynamics::FiniteThermoplasticHexMotion;
use fs_solid::hyper3::plastic::hex::{FiniteHexFormulation, FiniteThermoplasticHexSolid};
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticControls;
use fs_solid::plastic3::IncrementSettings;
use fs_solid::plastic3::hex::HexAssemblyBudget;

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const HEX: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const MU: f64 = 8000.0;
const BULK: f64 = 16000.0;
const YIELD: f64 = 200.0;
const HARDENING: f64 = 4000.0;
const ETA: f64 = 8000.0;
const ACTIVATION: f64 = 3000.0;
const DENSITY: [f64; 1] = [6.0];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 178,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 500.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: BULK,
                yield_stress_pa: YIELD,
                hardening_modulus_pa: HARDENING,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: 0.0,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: 0.0,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 3.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}
fn initial(
    law: AffineThermalFiniteJ2,
    formulation: FiniteHexFormulation,
    cx: &Cx<'_>,
) -> FiniteThermoplasticHexSolid<'static> {
    FiniteThermoplasticHexSolid::try_new(
        &NODES,
        &HEX,
        vec![law],
        &[300.0],
        formulation,
        HexAssemblyBudget {
            maximum_tangent_entries: 768,
            ..HexAssemblyBudget::default()
        },
        cx,
    )
    .unwrap()
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 2e-6,
            maximum_endpoint_work_defect_j: 200.0,
            yield_tolerance_pa: 2e-6,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn norm2(v: [f64; 3]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

#[test]
fn completely_free_hex_obeys_exact_uniform_backward_euler_momentum_with_original_mass() {
    with_cx(|cx| {
        let velocity0 = [0.2, -0.1, 0.05];
        let dt = 0.02;
        let mass = DENSITY[0]; // Reference volume is exactly one cubic metre.
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            // A static-only matrix budget cannot silently admit the added mass
            // blocks. The subsequent free-motion solve uses the exact full bound.
            let static_only = FiniteThermoplasticHexSolid::try_new(
                &NODES,
                &HEX,
                vec![family()],
                &[300.0],
                formulation,
                HexAssemblyBudget {
                    maximum_tangent_entries: 576,
                    ..HexAssemblyBudget::default()
                },
                cx,
            )
            .unwrap();
            let refusal =
                FiniteThermoplasticHexMotion::try_new(static_only, &DENSITY, &[velocity0; 8], cx)
                    .unwrap_err();
            assert!(refusal.to_string().contains("inertial tangent entries"));
            for acceleration in [[0.0; 3], [1.0, -0.5, 0.25]] {
                let law = family();
                let mut motion = FiniteThermoplasticHexMotion::try_new(
                    initial(law.clone(), formulation, cx),
                    &DENSITY,
                    &[velocity0; 8],
                    cx,
                )
                .unwrap();
                // For an affine Hex8, each consistent-mass row sums to rho*V/8.
                // These nodal loads therefore give the same uniform acceleration
                // at every free coordinate, with no elastic restraint or strain.
                let forces: Vec<_> = (0..24)
                    .map(|dof| (dof, mass * acceleration[dof % 3] / 8.0))
                    .collect();
                for interval in 1..=3 {
                    let count = f64::from(interval);
                    let old_v = std::array::from_fn(|i| {
                        velocity0[i] + (count - 1.0) * dt * acceleration[i]
                    });
                    let new_v =
                        std::array::from_fn(|i| velocity0[i] + count * dt * acceleration[i]);
                    let new_u: [f64; 3] = std::array::from_fn(|i| {
                        count * dt * velocity0[i]
                            + 0.5 * count * (count + 1.0) * dt * dt * acceleration[i]
                    });
                    let saved_u = motion.solid().displacement_m().to_vec();
                    let saved_v = motion.velocity_mps().to_vec();
                    let step = motion
                        .advance(&[300.0], &[], &forces, dt, settings(), controls(), cx)
                        .unwrap();
                    assert!(step.equilibrium.reactions_n.is_empty());
                    for (u, v) in step
                        .motion
                        .solid()
                        .displacement_m()
                        .iter()
                        .zip(step.motion.velocity_mps())
                    {
                        for i in 0..3 {
                            near(u[i], new_u[i], 2e-10, 2e-9);
                            near(v[i], new_v[i], 2e-8, 2e-9);
                        }
                    }
                    near(step.motion.reference_mass_kg(), mass, 2e-13, 0.0);
                    near(step.inertia.reference_mass_kg, mass, 2e-13, 0.0);
                    assert_eq!(step.motion.reference_density_kg_m3(), DENSITY);
                    assert_eq!(step.motion.solid().families()[0], law);
                    assert_eq!(step.motion.solid().formulation(), formulation);
                    assert!(
                        step.motion
                            .solid()
                            .quadrature_states()
                            .iter()
                            .all(|state| state == &law.initial_state())
                    );
                    near(step.motion.solid().free_energy_j(), 0.0, 1e-10, 0.0);
                    assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                    assert_eq!(step.material_energy.creep_dissipation_j, 0.0);
                    near(step.material_energy.caloric_heat_j, 0.0, 1e-10, 0.0);
                    let before = 0.5 * mass * norm2(old_v);
                    let after = 0.5 * mass * norm2(new_v);
                    let loss = 0.5 * mass * dt * dt * norm2(acceleration);
                    let work = mass * dt * (0..3).map(|i| acceleration[i] * new_v[i]).sum::<f64>();
                    near(step.inertia.kinetic_before_j, before, 3e-9, 2e-9);
                    near(step.inertia.kinetic_after_j, after, 3e-9, 2e-9);
                    near(step.inertia.kinetic_change_j, after - before, 3e-9, 2e-9);
                    near(step.inertia.numerical_loss_j, loss, 3e-9, 2e-9);
                    near(step.inertia.inertial_work_j, work, 3e-9, 2e-9);
                    near(
                        step.material_energy.dead_load_endpoint_work_j,
                        work,
                        3e-9,
                        2e-9,
                    );
                    near(
                        step.material_energy.support_endpoint_work_j,
                        0.0,
                        1e-12,
                        0.0,
                    );
                    near(
                        step.material_energy.internal_endpoint_work_j,
                        0.0,
                        1e-10,
                        0.0,
                    );
                    near(step.balance_error_j, 0.0, 3e-9, 0.0);
                    assert_eq!(motion.solid().displacement_m(), saved_u);
                    assert_eq!(motion.velocity_mps(), saved_v);
                    motion = step.motion;
                }
            }
        }
    });
}

#[test]
fn held_yielded_hex_relaxes_creep_and_arrests_motion_without_heating_numerical_kinetic_loss() {
    with_cx(|cx| {
        let law = family()
            .with_norton_creep(
                FiniteJ2ThermalOverstress {
                    reference_viscosity_pa_s: ETA,
                    stress_scale_pa: 400.0,
                    exponent: 1.0,
                    activation_temperature_k: ACTIVATION,
                },
                3.0,
            )
            .unwrap();
        let (a, initial_dt) = (0.16_f64, 0.2);
        let stretch = [a.exp(), (-0.5 * a).exp(), (-0.5 * a).exp()];
        let fixed: Vec<_> = NODES
            .iter()
            .flat_map(|x| (0..3).map(move |i| (stretch[i] - 1.0) * x[i]))
            .enumerate()
            .collect();
        // Isochoric coaxial loading: a=dp+dc+q/(3mu),
        // q=Y+H*dp and dc=dt*q/eta. This solves both mechanisms together.
        let first_q =
            (YIELD + HARDENING * a) / (1.0 + HARDENING / (3.0 * MU) + HARDENING * initial_dt / ETA);
        let plastic = (first_q - YIELD) / HARDENING;
        let first_creep = initial_dt * first_q / ETA;
        let initial_velocity = [0.2, -0.1, 0.05];
        let arrest_loss = 0.5 * DENSITY[0] * norm2(initial_velocity);
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let loaded = initial(law.clone(), formulation, cx)
                .advance(
                    &[300.0],
                    &fixed,
                    &[],
                    initial_dt,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
            for state in loaded.quadrature_states() {
                near(state.alpha, plastic, 2e-11, 2e-10);
                near(state.creep_alpha, first_creep, 2e-11, 2e-10);
            }
            let mut motion =
                FiniteThermoplasticHexMotion::try_new(loaded, &DENSITY, &[initial_velocity; 8], cx)
                    .unwrap();
            let (mut old_q, mut old_creep) = (first_q, first_creep);
            for (index, dt) in [0.1, 0.2, 0.3].into_iter().enumerate() {
                let temperature = 350.0;
                let viscosity = ETA * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
                let q = old_q / (1.0 + 3.0 * MU * dt / viscosity);
                let dc = dt * q / viscosity;
                let creep = old_creep + dc;
                let heat = q * dc;
                let free = q * q / (6.0 * MU) + 0.5 * HARDENING * plastic * plastic;
                let old_free = old_q * old_q / (6.0 * MU) + 0.5 * HARDENING * plastic * plastic;
                let saved_state = motion.solid().quadrature_states().to_vec();
                let saved_u = motion.solid().displacement_m().to_vec();
                let saved_v = motion.velocity_mps().to_vec();
                if index == 0 {
                    let mut refused = controls();
                    refused.maximum_numerical_kinetic_loss_j = 0.0;
                    assert!(
                        motion
                            .advance(&[temperature], &fixed, &[], dt, settings(), refused, cx)
                            .is_err()
                    );
                    assert_eq!(motion.solid().quadrature_states(), saved_state);
                    assert_eq!(motion.solid().displacement_m(), saved_u);
                    assert_eq!(motion.velocity_mps(), saved_v);
                    assert_eq!(motion.solid().temperatures_k(), &[300.0]);
                    assert_eq!(motion.solid().families()[0], law);
                }
                let step = motion
                    .advance(&[temperature], &fixed, &[], dt, settings(), controls(), cx)
                    .unwrap();
                for ((state, p), old) in step
                    .motion
                    .solid()
                    .quadrature_states()
                    .iter()
                    .zip(step.motion.solid().quadrature_first_piola_pa())
                    .zip(&saved_state)
                {
                    assert_eq!(state.alpha, old.alpha);
                    near(state.creep_alpha, creep, 2e-11, 2e-10);
                    near(
                        state.inverse_plastic_cauchy_green[0],
                        (-2.0 * (plastic + creep)).exp(),
                        2e-10,
                        2e-9,
                    );
                    near(p[0], (2.0 * q / 3.0) / stretch[0], 2e-6, 2e-9);
                    near(p[4], (-q / 3.0) / stretch[1], 2e-6, 2e-9);
                    near(p[8], (-q / 3.0) / stretch[2], 2e-6, 2e-9);
                }
                assert_eq!(step.motion.solid().displacement_m(), saved_u);
                assert!(
                    step.motion
                        .velocity_mps()
                        .iter()
                        .flatten()
                        .all(|v| *v == 0.0)
                );
                near(step.motion.reference_mass_kg(), DENSITY[0], 2e-13, 0.0);
                near(step.motion.solid().free_energy_j(), free, 2e-6, 2e-9);
                near(step.motion.solid().internal_energy_j(), free, 2e-6, 2e-9);
                assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                near(step.material_energy.creep_dissipation_j, heat, 2e-6, 2e-9);
                near(step.material_energy.caloric_heat_j, heat, 2e-6, 2e-9);
                near(step.material_energy.reversible_heat_j, 0.0, 1e-10, 0.0);
                near(
                    step.material_energy.endpoint_work_defect_j,
                    -(free - old_free) - heat,
                    2e-6,
                    2e-9,
                );
                near(
                    step.material_energy.projection_work_defect_j,
                    0.0,
                    1e-9,
                    0.0,
                );
                near(
                    step.material_energy.internal_endpoint_work_j,
                    0.0,
                    1e-9,
                    0.0,
                );
                near(step.material_energy.support_endpoint_work_j, 0.0, 1e-9, 0.0);
                let expected_arrest = if index == 0 { arrest_loss } else { 0.0 };
                near(step.inertia.numerical_loss_j, expected_arrest, 2e-12, 2e-10);
                near(step.inertia.kinetic_before_j, expected_arrest, 2e-12, 2e-10);
                near(step.inertia.kinetic_after_j, 0.0, 1e-12, 0.0);
                near(
                    step.inertia.kinetic_change_j,
                    -expected_arrest,
                    2e-12,
                    2e-10,
                );
                near(step.inertia.inertial_work_j, 0.0, 1e-12, 0.0);
                near(step.balance_error_j, 0.0, 2e-6, 0.0);
                assert_eq!(motion.solid().quadrature_states(), saved_state);
                assert_eq!(motion.velocity_mps(), saved_v);
                motion = step.motion;
                old_q = q;
                old_creep = creep;
            }
            assert!(old_q < 0.1 * first_q && old_creep > first_creep);
        }
    });
}
