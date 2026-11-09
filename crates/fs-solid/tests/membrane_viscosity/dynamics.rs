//! Independent nonlinear material-loss checks; all coefficients are synthetic.
use super::*;
use fs_solid::membrane::dynamics::{
    ConsistentMassBudget, MembraneMotion,
    viscosity::{DampedMembraneMotion, MembraneDampingBudget, MembraneDampingError},
};

fn budget() -> MembraneDampingBudget {
    MembraneDampingBudget {
        max_source_rows: 6,
        max_mass_solves_per_step: 28,
        max_step_s: 0.05,
        max_energy_defect_j: 1.0,
    }
}

fn initial_motion<'a>(
    patch: &'a MembranePatch,
    displacement: &[[f64; 3]],
    velocity: &[[f64; 3]],
    free: &[usize],
    consistent: bool,
) -> MembraneMotion<'a> {
    let fixed: Vec<_> = (0..3 * displacement.len())
        .filter(|dof| !free.contains(dof))
        .map(|dof| (dof, displacement[dof / 3][dof % 3]))
        .collect();
    if consistent {
        MembraneMotion::new_with_consistent_inertia(
            patch,
            &fixed,
            displacement,
            velocity,
            0.0,
            ConsistentMassBudget {
                max_nodes_per_axis: 3,
            },
            || false,
        )
        .unwrap()
    } else {
        MembraneMotion::new(patch, &fixed, displacement, velocity, 0.0, || false).unwrap()
    }
}

fn scalar_motion(patch: &MembranePatch, z: f64, v: f64, consistent: bool) -> MembraneMotion<'_> {
    initial_motion(
        patch,
        &[[0.0; 3], [0.08, 0.0, z], [0.0; 3]],
        &[[0.0; 3], [0.0, 0.0, v], [0.0; 3]],
        &[5],
        consistent,
    )
}

fn scalar_energy(z: f64) -> f64 {
    300.0 * (0.0832 + 0.5 * z * z).powi(2)
}

struct ScalarStep {
    z: f64,
    v: f64,
    heat: f64,
    work: f64,
    defect: f64,
}

// Direct reduction of the authored sheet: U=300(.0832+z^2/2)^2 and
// C(z)=20*z^2. No membrane force, tangent, factor or time-step API is used.
fn scalar_step(z0: f64, v0: f64, mass: f64, dt: f64, load: f64) -> ScalarStep {
    let force = |z: f64| 49.92 * z + 300.0 * z.powi(3);
    let loss_half = |z: f64, v: f64| {
        let exponent = -20.0 * z * z * dt / mass;
        (
            v * (0.5 * exponent).exp(),
            -0.5 * mass * v * v * exponent.exp_m1(),
        )
    };
    let (va, qa) = loss_half(z0, v0);
    let half = va + 0.5 * dt * (load - force(z0)) / mass;
    let z = z0 + dt * half;
    let vb = half + 0.5 * dt * (load - force(z)) / mass;
    let (v, qb) = loss_half(z, vb);
    let work = load * (z - z0);
    ScalarStep {
        z,
        v,
        heat: qa + qb,
        work,
        defect: scalar_energy(z) - scalar_energy(z0) + 0.5 * mass * (vb - va) * (vb + va) - work,
    }
}

#[test]
fn a_flat_prestretched_sheet_develops_material_loss_only_after_transverse_motion() {
    let (p, laws) = patch([20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0], 0.002, 0.0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0; 3]];
    for consistent in [false, true] {
        let original = scalar_motion(&p, 0.0, 1.0, consistent);
        let before = viscosity
            .evaluate(original.displacement_m(), original.velocity_m_s(), || false)
            .unwrap();
        assert_eq!(before.power_w, 0.0);
        assert!(before.resisting_force_n.iter().flatten().all(|f| *f == 0.0));
        let mass = if consistent { 0.5 } else { 1.0 };
        let expected = scalar_step(0.0, 1.0, mass, 0.04, 0.7);
        let mut motion =
            DampedMembraneMotion::new(&viscosity, original, budget(), || false).unwrap();
        let step = motion.step(0.04, &loads, || false).unwrap();
        close(motion.motion().displacement_m()[1][2], expected.z);
        close(motion.motion().velocity_m_s()[1][2], expected.v);
        close(step.dissipated_j, expected.heat);
        close(step.energy_defect_j, expected.defect);
        assert!(step.dissipated_j > 1e-4);
        close(motion.element_dissipated_j()[0], step.dissipated_j);
        close(motion.cumulative_dissipated_j(), step.dissipated_j);
        assert_eq!(step.mass_solves == 0, !consistent);
        assert!(step.mass_solves <= 28);
    }
}

#[test]
fn finite_amplitude_heat_work_and_support_impulse_match_the_scalar_source_law() {
    let (p, laws) = patch([20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0], 0.002, 0.0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0; 3]];
    for consistent in [false, true] {
        let mass = if consistent { 0.5 } else { 1.0 };
        let mut z = 0.2;
        let mut v = -0.3;
        let mut heat = 0.0;
        let mut motion = DampedMembraneMotion::new(
            &viscosity,
            scalar_motion(&p, z, v, consistent),
            budget(),
            || false,
        )
        .unwrap();
        for _ in 0..4 {
            let expected = scalar_step(z, v, mass, 0.04, 0.7);
            let step = motion.step(0.04, &loads, || false).unwrap();
            close(motion.motion().displacement_m()[1][2], expected.z);
            close(motion.motion().velocity_m_s()[1][2], expected.v);
            close(step.strain_energy_j, scalar_energy(expected.z));
            close(step.kinetic_energy_j, 0.5 * mass * expected.v * expected.v);
            close(step.dead_load_work_j, expected.work);
            close(step.dissipated_j, expected.heat);
            close(step.energy_defect_j, expected.defect);
            let energy_change = scalar_energy(expected.z) - scalar_energy(z)
                + 0.5 * mass * (expected.v - v) * (expected.v + v);
            close(
                energy_change + step.dissipated_j - step.dead_load_work_j,
                step.energy_defect_j,
            );
            // The full reference row mass is 1 even when M_ff=.5. The other
            // .5 must appear through prescribed/free consistent-mass coupling.
            close(step.support_impulse_n_s[2], expected.v - v - 0.04 * 0.7);
            close(step.support_impulse_n_s[0], 0.0);
            close(step.support_impulse_n_s[1], 0.0);
            assert!(step.dissipated_j > 0.0);
            heat += expected.heat;
            close(motion.cumulative_dissipated_j(), heat);
            close(motion.element_dissipated_j()[0], expected.heat);
            assert_eq!(motion.motion().displacement_m()[0], [0.0; 3]);
            assert_eq!(motion.motion().displacement_m()[2], [0.0; 3]);
            assert_eq!(motion.motion().displacement_m()[1][0], 0.08);
            assert_eq!(motion.motion().displacement_m()[1][1], 0.0);
            z = expected.z;
            v = expected.v;
        }
    }
}

#[test]
fn pressure_drives_the_native_damped_membrane_without_becoming_material_heat() {
    let (p, laws) = patch([20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0], 0.002, 0.0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0; 3]];
    let dt = 0.04;
    for consistent in [false, true] {
        let mass = if consistent { 0.5 } else { 1.0 };
        for pressure in [-2.0, 2.0] {
            let mut z = 0.2;
            let mut v = -0.3;
            let mut motion = DampedMembraneMotion::new(
                &viscosity,
                scalar_motion(&p, z, v, consistent),
                budget(),
                || false,
            )
            .unwrap();
            for _ in 0..4 {
                // x_b=(1.08,0,z), x_c=(0,1,0): each triangle vertex
                // receives p*(-z,0,1.08)/6. Only b_z is free.
                let expected = scalar_step(z, v, mass, dt, 0.7 + 0.18 * pressure);
                let step = motion
                    .step_with_pressure(dt, &loads, pressure, || false)
                    .unwrap();
                close(motion.motion().displacement_m()[1][2], expected.z);
                close(motion.motion().velocity_m_s()[1][2], expected.v);
                close(step.pressure_work_j, 0.18 * pressure * (expected.z - z));
                close(step.dead_load_work_j, 0.7 * (expected.z - z));
                close(step.dissipated_j, expected.heat);
                close(step.energy_defect_j, expected.defect);
                close(
                    scalar_energy(expected.z) - scalar_energy(z)
                        + 0.5 * mass * (expected.v - v) * (expected.v + v)
                        + step.dissipated_j
                        - step.dead_load_work_j
                        - step.pressure_work_j,
                    step.energy_defect_j,
                );
                // Includes pressure on the two supported vertices and M_cf.
                close(
                    step.support_impulse_n_s[2],
                    expected.v - v - dt * (0.7 + 0.54 * pressure),
                );
                close(
                    step.support_impulse_n_s[0],
                    0.25 * dt * pressure * (z + expected.z),
                );
                close(step.support_impulse_n_s[1], 0.0);
                z = expected.z;
                v = expected.v;
            }
            let before = motion.motion().displacement_m().to_vec();
            let before_heat = motion.cumulative_dissipated_j();
            assert!(
                motion
                    .step_with_pressure(dt, &loads, f64::NAN, || false)
                    .is_err()
            );
            assert_eq!(motion.motion().displacement_m(), before);
            assert_eq!(motion.cumulative_dissipated_j(), before_heat);
        }
    }
}

// Independent two-coordinate continuous equations from F0=(1.08,.12,z),
// F1=(-.04,.94,w), A=.5, h=.002, and the full authored viscosity tensor.
// The last coordinate integrates physical heat, independently of energy balance.
fn coupled_rhs(x: [f64; 5]) -> [f64; 5] {
    let [z, w, v, t, _] = x;
    let strain = [0.0904 + 0.5 * z * z, -0.0574 + 0.5 * w * w, 0.0696 + z * w];
    let stress = [1200.0 * strain[0], 600.0 * strain[1], 250.0 * strain[2]];
    let rate = [z * v, w * t, w * v + z * t];
    let viscous = [
        40_000.0 * rate[0] + 10_000.0 * rate[1] + 4000.0 * rate[2],
        10_000.0 * rate[0] + 30_000.0 * rate[1] + 2000.0 * rate[2],
        4000.0 * rate[0] + 2000.0 * rate[1] + 20_000.0 * rate[2],
    ];
    [
        v,
        t,
        0.7 - 0.5 * (stress[0] * z + stress[2] * w) - 0.001 * (viscous[0] * z + viscous[2] * w),
        -0.4 - 0.5 * (stress[2] * z + stress[1] * w) - 0.001 * (viscous[2] * z + viscous[1] * w),
        0.001 * (0..3).map(|i| rate[i] * viscous[i]).sum::<f64>(),
    ]
}

fn coupled_reference() -> [f64; 5] {
    let dt = 0.12 / 4096.0;
    let mut x = [0.19, 0.25, 0.4, -0.2, 0.0];
    for _ in 0..4096 {
        let a = coupled_rhs(x);
        let b = coupled_rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = coupled_rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = coupled_rhs(std::array::from_fn(|i| x[i] + dt * c[i]));
        x = std::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn coupled_nonlinear_motion_and_material_heat_converge_at_second_order() {
    let (p, laws) = patch(
        [40_000.0, 30_000.0, 20_000.0, 10_000.0, 4000.0, 2000.0],
        0.002,
        0.0,
    );
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0, 0.0, -0.4]];
    let reference = coupled_reference();
    let mut errors = Vec::new();
    for steps in [12, 24, 48] {
        let original = initial_motion(
            &p,
            &deformed(),
            &[[0.0; 3], [0.0, 0.0, 0.4], [0.0, 0.0, -0.2]],
            &[5, 8],
            false,
        );
        let mut motion =
            DampedMembraneMotion::new(&viscosity, original, budget(), || false).unwrap();
        for _ in 0..steps {
            motion.step(0.12 / steps as f64, &loads, || false).unwrap();
        }
        let u = motion.motion().displacement_m();
        let v = motion.motion().velocity_m_s();
        // A fixed .1-second observation scale combines position and velocity.
        let state_error = ((u[1][2] - reference[0]).powi(2)
            + (u[2][2] - reference[1]).powi(2)
            + 0.01 * (v[1][2] - reference[2]).powi(2)
            + 0.01 * (v[2][2] - reference[3]).powi(2))
        .sqrt();
        errors.push((
            state_error,
            (motion.cumulative_dissipated_j() - reference[4]).abs(),
        ));
    }
    for pair in errors.windows(2) {
        assert!(pair[0].0 > 3.5 * pair[1].0, "state refinement: {errors:?}");
        assert!(pair[0].1 > 3.5 * pair[1].1, "heat refinement: {errors:?}");
    }
    assert!(errors[2].0 < 5e-5 && errors[2].1 < 5e-5, "{errors:?}");
}

fn momenta(motion: &MembraneMotion<'_>) -> ([f64; 3], [f64; 3]) {
    let momentum = motion
        .patch()
        .reference_mass_apply(motion.inertia_image(), motion.velocity_m_s(), || false)
        .unwrap();
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (node, p) in momentum.iter().enumerate() {
        let u = motion.displacement_m()[node];
        let reference = motion.patch().mesh().nodes[node];
        let x = [reference[0] + u[0], reference[1] + u[1], u[2]];
        for axis in 0..3 {
            linear[axis] += p[axis];
        }
        angular[0] += x[1] * p[2] - x[2] * p[1];
        angular[1] += x[2] * p[0] - x[0] * p[2];
        angular[2] += x[0] * p[1] - x[1] * p[0];
    }
    (linear, angular)
}

#[test]
fn internal_material_loss_preserves_unforced_linear_and_angular_momentum() {
    let (p, laws) = patch(
        [40_000.0, 30_000.0, 20_000.0, 10_000.0, 4000.0, 2000.0],
        0.002,
        0.37,
    );
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let velocity = [[0.2, -0.3, 0.1], [-0.4, 0.1, 0.8], [0.3, 0.6, -0.2]];
    for consistent in [false, true] {
        let original = initial_motion(
            &p,
            &deformed(),
            &velocity,
            &(0..9).collect::<Vec<_>>(),
            consistent,
        );
        let (linear, angular) = momenta(&original);
        let initial_energy = original.strain_energy_j() + original.kinetic_energy_j();
        let mut motion =
            DampedMembraneMotion::new(&viscosity, original, budget(), || false).unwrap();
        let mut defects = 0.0;
        for _ in 0..24 {
            let step = motion.step(0.001, &[[0.0; 3]; 3], || false).unwrap();
            defects += step.energy_defect_j;
            let (after_linear, after_angular) = momenta(motion.motion());
            for axis in 0..3 {
                close(after_linear[axis], linear[axis]);
                close(after_angular[axis], angular[axis]);
                close(step.support_impulse_n_s[axis], 0.0);
            }
            assert!(step.dissipated_j >= 0.0);
            close(step.dead_load_work_j, 0.0);
        }
        assert!(motion.cumulative_dissipated_j() > 1e-3);
        let final_energy = motion.motion().strain_energy_j() + motion.motion().kinetic_energy_j();
        assert!(
            (final_energy - initial_energy + motion.cumulative_dissipated_j() - defects).abs()
                < 2e-11
        );
    }
}

#[derive(Debug, PartialEq)]
struct Accepted {
    displacement: Vec<[f64; 3]>,
    velocity: Vec<[f64; 3]>,
    time: f64,
    strain: f64,
    kinetic: f64,
}

fn accepted(motion: &MembraneMotion<'_>) -> Accepted {
    Accepted {
        displacement: motion.displacement_m().to_vec(),
        velocity: motion.velocity_m_s().to_vec(),
        time: motion.time_s(),
        strain: motion.strain_energy_j(),
        kinetic: motion.kinetic_energy_j(),
    }
}

#[test]
fn cancellation_at_each_stage_preserves_the_last_accepted_state_and_heat() {
    let (p, laws) = patch([20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0], 0.002, 0.0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0; 3]];
    let mut motion = DampedMembraneMotion::new(
        &viscosity,
        scalar_motion(&p, 0.2, -0.3, true),
        budget(),
        || false,
    )
    .unwrap();
    motion.step(0.01, &loads, || false).unwrap();
    let state = accepted(motion.motion());
    let element_heat = motion.element_dissipated_j().to_vec();
    let total_heat = motion.cumulative_dissipated_j();
    assert!(total_heat > 0.0);
    // Discover the successful path's polls, then stop at EVERY one. This covers
    // both loss stages, the elastic trial and the final publication boundary.
    let mut probe =
        DampedMembraneMotion::new(&viscosity, motion.motion().clone(), budget(), || false).unwrap();
    let mut count = 0;
    probe
        .step(0.01, &loads, || {
            count += 1;
            false
        })
        .unwrap();
    assert!(count > 4);
    for stop_at in 1..=count {
        let mut polls = 0;
        let result = motion.step(0.01, &loads, || {
            polls += 1;
            polls == stop_at
        });
        assert!(
            result.is_err(),
            "uncaught cancellation poll {stop_at}/{count}"
        );
        assert_eq!(accepted(motion.motion()), state);
        assert_eq!(motion.element_dissipated_j(), element_heat.as_slice());
        assert_eq!(motion.cumulative_dissipated_j(), total_heat);
    }
    motion.step(0.01, &loads, || false).unwrap();
    assert!(motion.cumulative_dissipated_j() > total_heat);
}

#[test]
fn source_work_time_and_energy_refusals_publish_no_partial_motion_or_heat() {
    let values = [20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0];
    let (p, laws) = patch(values, 0.002, 0.0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let (foreign, _) = patch(values, 0.002, 0.0);
    assert!(
        DampedMembraneMotion::new(
            &viscosity,
            scalar_motion(&foreign, 0.2, -0.3, false),
            budget(),
            || false,
        )
        .is_err()
    );
    assert!(
        DampedMembraneMotion::new(
            &viscosity,
            scalar_motion(&p, 0.2, -0.3, false),
            MembraneDampingBudget {
                max_source_rows: 0,
                ..budget()
            },
            || false,
        )
        .is_err()
    );
    assert!(
        DampedMembraneMotion::new(
            &viscosity,
            scalar_motion(&p, 0.2, -0.3, true),
            MembraneDampingBudget {
                max_mass_solves_per_step: 0,
                ..budget()
            },
            || false,
        )
        .is_err()
    );
    let mut motion = DampedMembraneMotion::new(
        &viscosity,
        scalar_motion(&p, 0.2, -0.3, false),
        MembraneDampingBudget {
            max_energy_defect_j: 1e-14,
            ..budget()
        },
        || false,
    )
    .unwrap();
    let state = accepted(motion.motion());
    let loads = [[0.0; 3], [0.0, 0.0, 0.7], [0.0; 3]];
    for dt in [0.0, -0.01, 0.051, f64::NAN] {
        assert!(motion.step(dt, &loads, || false).is_err());
        assert_eq!(accepted(motion.motion()), state);
        assert_eq!(motion.cumulative_dissipated_j(), 0.0);
    }
    match motion.step(0.04, &loads, || false) {
        Err(MembraneDampingError::EnergyBudget {
            defect_j,
            maximum_j,
        }) => {
            assert!(defect_j.abs() > maximum_j);
            assert_eq!(maximum_j, 1e-14);
        }
        result => panic!("expected the independent energy-defect gate, got {result:?}"),
    }
    assert_eq!(accepted(motion.motion()), state);
    assert_eq!(motion.cumulative_dissipated_j(), 0.0);
    assert!(motion.element_dissipated_j().iter().all(|q| *q == 0.0));
}
