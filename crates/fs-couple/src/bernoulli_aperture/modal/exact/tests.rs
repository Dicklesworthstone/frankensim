//! Independent physical-momentum equations and power quadrature for regional memory.
use super::*;
use crate::bernoulli_aperture::modal::ApertureMode;
use crate::bernoulli_aperture::paired::PairedApertureLimits;

fn spec() -> ModalApertureSpec {
    ModalApertureSpec {
        rest_opening_m: 4e-4,
        width_m: 0.012,
        modes: vec![
            ApertureMode {
                mass_kg: 2e-5,
                stiffness_n_m: 600.0,
                damping_n_s_m: 0.02,
                pressure_area_m2: 1.2e-4,
                gap_projection: 1.0,
                force_projection: [1.0, 0.0],
            },
            ApertureMode {
                mass_kg: 1.3e-5,
                stiffness_n_m: 1700.0,
                damping_n_s_m: 0.013,
                pressure_area_m2: -3e-5,
                gap_projection: -0.3,
                force_projection: [0.0, 1.0],
            },
        ],
        density_kg_m3: 1.2,
        impedance_pa_s_m3: 1e6,
        time_step_s: 1.0 / 48_000.0,
        max_steps: 1024,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 1.0,
            max_evaluations: 128,
            velocity_tolerance_m_s: 1e-13,
            relative_tolerance: 1e-11,
        },
    }
}
fn arms() -> Vec<CoupledArm> {
    vec![
        CoupledArm {
            coupling: vec![5.0, -7.0],
            relaxation_time_s: 9e-5,
        },
        CoupledArm {
            coupling: vec![0.8, 6.0],
            relaxation_time_s: 3e-4,
        },
    ]
}
fn initial() -> Vec<BladeState> {
    vec![
        BladeState {
            displacement_m: 2e-5,
            velocity_m_s: 0.13,
        },
        BladeState {
            displacement_m: -7e-6,
            velocity_m_s: -0.08,
        },
    ]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-18),
        "{a} != {b}"
    );
}

// This oracle uses physical positions, momenta and forces. It does not use
// the production energy normalization, increment rows or loss Gramian.
fn trajectory(
    spec: &ModalApertureSpec,
    arms: &[CoupledArm],
    old: &[BladeState],
    memory: &[f64],
    held: &[f64],
    time: f64,
) -> Vec<f64> {
    let n = old.len();
    let mut c = vec![0.0; n*n];
    for i in 0..n { c[i*n+i] = spec.modes[i].damping_n_s_m; }
    trajectory_with_c(spec, &c, arms, old, memory, held, time)
}
fn trajectory_with_c(
    spec: &ModalApertureSpec,
    c: &[f64],
    arms: &[CoupledArm],
    old: &[BladeState],
    memory: &[f64],
    held: &[f64],
    time: f64,
) -> Vec<f64> {
    let n = old.len();
    let size = 3 * n + arms.len();
    let mut a = vec![0.0; size * size];
    let mut x = vec![0.0; size];
    for i in 0..n {
        let mode = spec.modes[i];
        x[i] = old[i].displacement_m;
        x[n + i] = mode.mass_kg * old[i].velocity_m_s;
        x[2 * n + i] = held[i];
        a[i * size + n + i] = 1.0 / mode.mass_kg;
        a[(n + i) * size + i] = -mode.stiffness_n_m;
        for j in 0..n {
            a[(n + i) * size + n + j] = -c[i*n+j] / spec.modes[j].mass_kg;
        }
        a[(n + i) * size + 2 * n + i] = 1.0;
    }
    for (j, arm) in arms.iter().enumerate() {
        let at = 3 * n + j;
        x[at] = memory[j];
        a[at * size + at] = -1.0 / arm.relaxation_time_s;
        for i in 0..n {
            a[(n + i) * size + at] = arm.coupling[i];
            a[at * size + i] = arm.coupling[i] / arm.relaxation_time_s;
            for k in 0..n {
                a[(n + i) * size + k] -= arm.coupling[i] * arm.coupling[k];
            }
        }
    }
    for value in &mut a {
        *value *= time;
    }
    let map = exponential(&a, size).unwrap();
    (0..size)
        .map(|i| (0..size).map(|j| map[i * size + j] * x[j]).sum())
        .collect()
}
fn energy(
    spec: &ModalApertureSpec,
    arms: &[CoupledArm],
    state: &[BladeState],
    memory: &[f64],
) -> f64 {
    let modes: f64 = spec
        .modes
        .iter()
        .zip(state)
        .map(|(m, s)| {
            0.5 * m.mass_kg * s.velocity_m_s * s.velocity_m_s
                + 0.5 * m.stiffness_n_m * s.displacement_m * s.displacement_m
        })
        .sum();
    modes
        + arms
            .iter()
            .zip(memory)
            .map(|(arm, z)| {
                let strain = arm
                    .coupling
                    .iter()
                    .zip(state)
                    .map(|(c, s)| c * s.displacement_m)
                    .sum::<f64>()
                    - z;
                0.5 * strain * strain
            })
            .sum::<f64>()
}

#[test]
fn g1_coupled_memory_matches_physical_momentum_and_independent_power() {
    let spec = spec();
    let arms = arms();
    let old = initial();
    let memory = vec![vec![3e-5, -1e-5]];
    let held = [0.0015, -0.0008];
    let dt = spec.time_step_s;
    let exact = ExactModes::with_coupled_arms(&spec, arms.clone()).unwrap();
    let mut velocity = vec![0.0; 2];
    exact
        .free_velocity(&old, &held, dt, Some(&memory), &mut velocity)
        .unwrap();
    let end = exact.endpoint(&old, &held, &velocity, dt, Some(&memory));
    let physical = trajectory(&spec, &arms, &old, &memory[0], &held, dt);
    for i in 0..2 {
        near(end.state[i].displacement_m, physical[i], 2e-10);
        near(
            end.state[i].velocity_m_s,
            physical[2 + i] / spec.modes[i].mass_kg,
            2e-10,
        );
        near(end.memory[0][i], physical[6 + i], 2e-10);
    }
    let panels = 64;
    let mut viscosity = 0.0;
    let mut relaxation = 0.0;
    let mut mean_q = [0.0; 2];
    let mut mean_force = [0.0; 2];
    for panel in 0..=panels {
        let time = dt * panel as f64 / panels as f64;
        let x = trajectory(&spec, &arms, &old, &memory[0], &held, time);
        let weight = if panel == 0 || panel == panels {
            1.0
        } else if panel % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for i in 0..2 {
            let v = x[2 + i] / spec.modes[i].mass_kg;
            viscosity += weight * spec.modes[i].damping_n_s_m * v * v;
            mean_q[i] += weight * x[i];
        }
        for (j, arm) in arms.iter().enumerate() {
            let strain = arm.coupling.iter().zip(&x).map(|(c, q)| c * q).sum::<f64>() - x[6 + j];
            relaxation += weight * strain * strain / arm.relaxation_time_s;
            for i in 0..2 {
                mean_force[i] += weight * arm.coupling[i] * strain;
            }
        }
    }
    let factor = dt / (3 * panels) as f64;
    near(end.viscous_loss_j, viscosity * factor, 2e-7);
    near(end.relaxation_loss_j, relaxation * factor, 2e-7);
    assert!(end.viscous_loss_j > 0.0 && end.relaxation_loss_j > 0.0);
    for i in 0..2 {
        near(end.mean_displacement[i], mean_q[i] * factor / dt, 2e-7);
        near(end.relaxation_force_n[i], mean_force[i] * factor / dt, 2e-7);
        let mode = spec.modes[i];
        near(
            mode.mass_kg * (end.state[i].velocity_m_s - old[i].velocity_m_s) / dt,
            held[i]
                - mode.stiffness_n_m * end.mean_displacement[i]
                - mode.damping_n_s_m * velocity[i]
                - end.relaxation_force_n[i],
            2e-10,
        );
    }
    let work: f64 = held
        .iter()
        .zip(old.iter().zip(&end.state))
        .map(|(f, (a, b))| f * (b.displacement_m - a.displacement_m))
        .sum();
    near(
        energy(&spec, &arms, &old, &memory[0]) + work,
        energy(&spec, &arms, &end.state, &end.memory[0])
            + end.viscous_loss_j
            + end.relaxation_loss_j,
        2e-10,
    );
    // The condensed pressure/contact maps must include all cross-mode arm forces.
    for (force, response) in [
        (
            spec.modes
                .iter()
                .map(|m| m.pressure_area_m2)
                .collect::<Vec<_>>(),
            exact.pressure_response(),
        ),
        (
            spec.modes
                .iter()
                .map(|m| m.gap_projection)
                .collect::<Vec<_>>(),
            exact.contact_response(),
        ),
    ] {
        let zero = vec![BladeState::default(); 2];
        let x = trajectory(&spec, &arms, &zero, &[0.0; 2], &force, dt);
        for i in 0..2 {
            near(response[i], x[i] / dt, 2e-10);
        }
    }
}

#[test]
fn g3_coupled_memory_respects_signed_physical_basis_congruence() {
    let original = spec();
    let arms = arms();
    let old = initial();
    let mut transformed = original.clone();
    let scales = [-2.0, 0.5];
    let mut transformed_arms = arms.clone();
    let mut transformed_old = old.clone();
    let held = [0.0015, -0.0008];
    let mut transformed_held = held;
    for i in 0..2 {
        let m = &mut transformed.modes[i];
        let scale = scales[i];
        m.mass_kg *= scale * scale;
        m.stiffness_n_m *= scale * scale;
        m.damping_n_s_m *= scale * scale;
        m.pressure_area_m2 *= scale;
        m.gap_projection *= scale;
        for f in &mut m.force_projection {
            *f *= scale;
        }
        transformed_old[i].displacement_m /= scale;
        transformed_old[i].velocity_m_s /= scale;
        transformed_held[i] *= scale;
        for arm in &mut transformed_arms {
            arm.coupling[i] *= scale;
        }
    }
    let a = ExactModes::with_coupled_arms(&original, arms).unwrap();
    let b = ExactModes::with_coupled_arms(&transformed, transformed_arms).unwrap();
    let memory = vec![vec![3e-5, -1e-5]];
    let dt = original.time_step_s;
    let (mut va, mut vb) = (vec![0.0; 2], vec![0.0; 2]);
    a.free_velocity(&old, &held, dt, Some(&memory), &mut va)
        .unwrap();
    b.free_velocity(
        &transformed_old,
        &transformed_held,
        dt,
        Some(&memory),
        &mut vb,
    )
    .unwrap();
    let x = a.endpoint(&old, &held, &va, dt, Some(&memory));
    let y = b.endpoint(&transformed_old, &transformed_held, &vb, dt, Some(&memory));
    for i in 0..2 {
        near(
            x.state[i].displacement_m,
            scales[i] * y.state[i].displacement_m,
            2e-12,
        );
        near(
            x.state[i].velocity_m_s,
            scales[i] * y.state[i].velocity_m_s,
            2e-12,
        );
        near(
            scales[i] * x.relaxation_force_n[i],
            y.relaxation_force_n[i],
            2e-12,
        );
        near(x.memory[0][i], y.memory[0][i], 2e-12);
    }
    near(x.viscous_loss_j, y.viscous_loss_j, 2e-12);
    near(x.relaxation_loss_j, y.relaxation_loss_j, 2e-12);
}

#[test]
fn g0_complete_memory_capacity_shape_and_empty_history_are_explicit() {
    let spec = spec();
    let old = initial();
    let dt = spec.time_step_s;
    assert_eq!(CoupledModes::dimensions(2, 122).unwrap(), (128, 132));
    assert!(CoupledModes::dimensions(2, 123).is_err());
    assert!(CoupledModes::dimensions(usize::MAX, 0).is_err());
    assert!(
        ExactModes::with_coupled_arms(
            &spec,
            vec![CoupledArm {
                coupling: vec![1.0],
                relaxation_time_s: 1e-4
            }]
        )
        .is_err()
    );
    assert!(
        ExactModes::with_coupled_arms(
            &spec,
            vec![CoupledArm {
                coupling: vec![1.0, f64::NAN],
                relaxation_time_s: 1e-4
            }]
        )
        .is_err()
    );
    assert!(
        ExactModes::with_coupled_arms(
            &spec,
            vec![CoupledArm {
                coupling: vec![1.0, 2.0],
                relaxation_time_s: 0.0
            }]
        )
        .is_err()
    );
    let exact = ExactModes::with_coupled_arms(&spec, arms()).unwrap();
    let mut out = vec![0.0; 2];
    assert!(
        exact
            .free_velocity(&old, &[0.0; 2], dt, None, &mut out)
            .is_err()
    );
    for invalid_memory in [
        vec![],
        vec![vec![0.0]],
        vec![vec![0.0, f64::NAN]],
        vec![vec![0.0; 2]; 2],
    ] {
        assert!(
            exact
                .free_velocity(&old, &[0.0; 2], dt, Some(&invalid_memory), &mut out)
                .is_err()
        );
    }
    let empty = ExactModes::with_coupled_arms(&spec, Vec::new()).unwrap();
    let legacy = ExactModes::new(&spec, None).unwrap();
    let empty_memory = vec![Vec::new()];
    let held = [0.0015, -0.0008];
    let mut reference = vec![0.0; 2];
    empty
        .free_velocity(&old, &held, dt, Some(&empty_memory), &mut out)
        .unwrap();
    legacy
        .free_velocity(&old, &held, dt, None, &mut reference)
        .unwrap();
    let a = empty.endpoint(&old, &held, &out, dt, Some(&empty_memory));
    let b = legacy.endpoint(&old, &held, &reference, dt, None);
    assert_eq!(a.memory, empty_memory);
    assert_eq!(a.relaxation_loss_j, 0.0);
    for i in 0..2 {
        near(a.state[i].displacement_m, b.state[i].displacement_m, 2e-10);
        near(a.state[i].velocity_m_s, b.state[i].velocity_m_s, 2e-10);
        near(a.mean_displacement[i], b.mean_displacement[i], 2e-10);
    }
    near(a.viscous_loss_j, b.viscous_loss_j, 2e-10);
}

#[test]
fn g1_full_kelvin_and_tensor_memory_share_momentum_and_actual_heat() {
    let mut spec = spec();
    for mode in &mut spec.modes { mode.damping_n_s_m = 0.0; }
    let c = [0.02, -0.008, -0.008, 0.013];
    let arms = arms();
    let old = initial();
    let memory = vec![vec![3e-5, -1e-5]];
    let held = [0.0015, -0.0008];
    let dt = spec.time_step_s;
    let exact = ExactModes::with_coupled_memory(&spec, &c, arms.clone()).unwrap();
    let mut mean_velocity = vec![0.0; 2];
    exact.free_velocity(&old, &held, dt, Some(&memory), &mut mean_velocity).unwrap();
    let end = exact.endpoint(&old, &held, &mean_velocity, dt, Some(&memory));
    let physical = trajectory_with_c(&spec, &c, &arms, &old, &memory[0], &held, dt);
    for i in 0..2 {
        near(end.state[i].displacement_m, physical[i], 2e-10);
        near(end.state[i].velocity_m_s, physical[2+i]/spec.modes[i].mass_kg, 2e-10);
        near(end.memory[0][i], physical[6+i], 2e-10);
        let reaction = (0..2).map(|j| c[i*2+j]*mean_velocity[j]).sum::<f64>();
        near(spec.modes[i].mass_kg*(end.state[i].velocity_m_s-old[i].velocity_m_s)/dt,
            held[i]-spec.modes[i].stiffness_n_m*end.mean_displacement[i]
                -reaction-end.relaxation_force_n[i], 2e-10);
    }
    // Independent physical trajectory quadrature; neither loss is a residual.
    let panels = 128;
    let mut viscous = 0.0;
    let mut relaxation = 0.0;
    for k in 0..=panels {
        let x = trajectory_with_c(&spec, &c, &arms, &old, &memory[0], &held,
            dt*k as f64/panels as f64);
        let weight = if k == 0 || k == panels { 1.0 } else if k%2 == 0 { 2.0 } else { 4.0 };
        let v = [x[2]/spec.modes[0].mass_kg, x[3]/spec.modes[1].mass_kg];
        viscous += weight*(v[0]*(c[0]*v[0]+c[1]*v[1])+v[1]*(c[2]*v[0]+c[3]*v[1]));
        for (j, arm) in arms.iter().enumerate() {
            let elastic = arm.coupling.iter().zip(&x).map(|(a,q)| a*q).sum::<f64>()-x[6+j];
            relaxation += weight*elastic*elastic/arm.relaxation_time_s;
        }
    }
    let scale = dt/(3*panels) as f64;
    near(end.viscous_loss_j, viscous*scale, 2e-8);
    near(end.relaxation_loss_j, relaxation*scale, 2e-8);
    assert!(end.viscous_loss_j > 0.0 && end.relaxation_loss_j > 0.0);
    let work = held.iter().zip(&end.state).zip(&old)
        .map(|((f,new),old)| f*(new.displacement_m-old.displacement_m)).sum::<f64>();
    near(energy(&spec,&arms,&old,&memory[0])+work,
        energy(&spec,&arms,&end.state,&end.memory[0])+end.viscous_loss_j+end.relaxation_loss_j,
        2e-10);
}
