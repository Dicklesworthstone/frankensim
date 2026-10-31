use super::*;
use crate::modal_acoustic_time::ModalAcousticState;

fn law() -> MeanStrainHistory {
    let wire = WireSpan {
        endpoints_m: [[0.0, 0.0], [0.5, 0.0]],
        linear_density_kg_m: 0.002,
        tension_n: 20.0,
        bending_stiffness_n_m2: 1e-4,
        damping_per_s: vec![0.0; 2],
    };
    let body = wire
        .stretching_body_two_planes(
            vec![ModalAcousticState::default(); 4],
            StringStretching {
                axial_rigidity_n: 20_000.0,
                maximum_slope: 0.2,
            },
        )
        .unwrap();
    let BodyPotential::String(potential) = body.potential else {
        panic!("original string owner")
    };
    MeanStrainHistory {
        potential,
        length_m: 0.5,
        material: StrainRelaxation::new(
            0,
            vec![
                StrainRelaxationArm {
                    stiffness_j: 2500.0,
                    relaxation_time_s: 0.001,
                },
                StrainRelaxationArm {
                    stiffness_j: 1000.0,
                    relaxation_time_s: 0.003,
                },
            ],
        )
        .unwrap(),
    }
}
fn strain(q: &[f64]) -> f64 {
    q.iter()
        .enumerate()
        .map(|(i, q)| {
            let k = ((i % 2) + 1) as f64 * core::f64::consts::PI / 0.5;
            q * q * k * k / (2.0 * 0.002 * 0.5)
        })
        .sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_mean_strain_history_closes_work_and_independent_arm_heat() {
    let law = law();
    let a = [2e-5, -1e-5, 1e-5, 0.5e-5];
    let b = [2.2e-5, -1.1e-5, 0.9e-5, 0.6e-5];
    let ea = strain(&a);
    let eb = strain(&b);
    close(law.strain(&a).unwrap(), ea, 1e-14);
    let memory = [ea - 2e-6, ea + 1e-6];
    let saved = memory;
    let mut next = [0.0; 2];
    let mut force = [0.0; 4];
    let mut heat = [0.0; 2];
    let dt = 0.001;
    let loss = law
        .stage_with_dissipation(&a, &b, &memory, dt, &mut next, &mut force, &mut heat)
        .unwrap();
    let mut mean_stress = 0.0;
    for (index, arm) in law.material.arms().iter().enumerate() {
        let u = dt / arm.relaxation_time_s;
        let e = (-u).exp();
        let ramp = (eb - ea) / u;
        let transient = ea - memory[index] - ramp;
        let expected_heat = arm.stiffness_j
            * (ramp * ramp * u
                + 2.0 * ramp * transient * (1.0 - e)
                + 0.5 * transient * transient * (1.0 - e * e));
        close(heat[index], expected_heat, 2e-11);
        close(next[index], eb - ramp - transient * e, 2e-12);
        mean_stress += arm.stiffness_j * (ramp + transient * (1.0 - e) / u);
    }
    for i in 0..4 {
        let k = ((i % 2) + 1) as f64 * core::f64::consts::PI / 0.5;
        let gradient = k * k / (0.002 * 0.5) * 0.5 * (a[i] + b[i]);
        close(force[i], -mean_stress * gradient, 2e-12);
    }
    close(loss, heat.iter().sum(), 1e-14);
    let work: f64 = force
        .iter()
        .zip(a.iter().zip(b))
        .map(|(f, (a, b))| f * (b - a))
        .sum();
    let change = law.energy(&b, &next).unwrap() - law.energy(&a, &memory).unwrap();
    assert!((change + loss + work).abs() < 1e-12 * (change.abs() + loss + work.abs()));
    assert_eq!(
        memory, saved,
        "a constitutive trial cannot advance accepted memory"
    );
}

#[test]
fn held_strain_relaxes_with_original_times_and_keeps_equilibrium_energy_separate() {
    let law = law();
    let q = [2e-5, 0.0, 1e-5, 0.0];
    let epsilon = strain(&q);
    let memory = [epsilon - 2e-6, epsilon - 1e-6];
    let mut next = [0.0; 2];
    let mut f = [0.0; 4];
    let dt = 0.002;
    let before = law.energy(&q, &memory).unwrap();
    let loss = law.stage(&q, &q, &memory, dt, &mut next, &mut f).unwrap();
    for (i, arm) in law.material.arms().iter().enumerate() {
        close(
            epsilon - next[i],
            (epsilon - memory[i]) * (-dt / arm.relaxation_time_s).exp(),
            1e-12,
        );
    }
    close(law.energy(&q, &next).unwrap() + loss, before, 1e-12);
    let relaxed = law.relaxed_memory(&q).unwrap();
    assert_eq!(law.energy(&q, &relaxed).unwrap(), 0.0);
    assert!(law.potential.observe(&q).unwrap().stretching_energy_j > 0.0);
    assert!(
        law.energy(&q, &[1.0; 2]).is_err(),
        "compressive total tension must refuse"
    );
}

#[test]
fn hereditary_tangent_bound_includes_both_planes_and_initial_stress() {
    let law = law();
    let q = [2e-5, -1e-5, 1e-5, 0.5e-5];
    let memory = law.relaxed_memory(&q).unwrap();
    let bound = law.stiffness_upper_bound_per_s2(&memory).unwrap();
    let direction = [0.2, -0.3, 0.1, 0.4];
    let h = 1e-10;
    let plus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q + h * d).collect();
    let minus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q - h * d).collect();
    let mut fp = [0.0; 4];
    let mut fm = [0.0; 4];
    law.reaction(&plus, &memory, &mut fp).unwrap();
    law.reaction(&minus, &memory, &mut fm).unwrap();
    let curvature: f64 = fp
        .iter()
        .zip(fm)
        .zip(direction)
        .map(|((p, m), d)| -d * (p - m) / (2.0 * h))
        .sum();
    assert!(curvature.abs() <= bound * direction.iter().map(|d| d * d).sum::<f64>());
    assert!(law.stiffness_upper_bound_per_s2(&[-1.0; 2]).unwrap() > bound);
    let rotated = [-q[2], -q[3], q[0], q[1]];
    close(
        law.strain(&rotated).unwrap(),
        law.strain(&q).unwrap(),
        1e-14,
    );
}
