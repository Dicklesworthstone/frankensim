//! Complete structural tones: physical resonance, mutual inertia and actual work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, PreparedHarmonicResponse};
use fs_solid::dynamics::{
    ExactHarmonicSolid, ExactLinearSolid, ExactLinearSolidControls, ExactLinearSolidError,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

fn with_cx(gate: &CancelGate, run: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 63,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    });
}
fn controls() -> ExactLinearSolidControls {
    ExactLinearSolidControls {
        maximum_free_dofs: 32,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_exponential_multiply_terms: 10_000_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 1000,
    }
}
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    let (mut stiffness, mut mass) = (Coo::new(n, n), Coo::new(n, n));
    for i in 0..n {
        for j in 0..n {
            stiffness.push(i, j, k[i * n + j]);
            mass.push(i, j, m[i * n + j]);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn real_at(z: C64, omega: f64, time: f64) -> f64 {
    let (sine, cosine) = (omega * time).sin_cos();
    z.re * cosine + z.im * sine
}
fn rotate(solid: &ExactHarmonicSolid, cosine: &mut [f64], sine: &mut [f64], tones: usize) {
    for i in 0..cosine.len() {
        (cosine[i], sine[i]) = solid
            .advance_force_phase(i % tones, cosine[i], sine[i])
            .unwrap();
    }
}
fn state_bits(solid: &ExactHarmonicSolid) -> Vec<u64> {
    solid
        .displacement_m()
        .iter()
        .chain(solid.velocity_mps())
        .copied()
        .chain([
            solid.time_s(),
            solid.stored_energy_j(),
            solid.kinetic_energy_j(),
            solid.cumulative_work_j(),
            solid.cumulative_dissipated_energy_j(),
        ])
        .map(f64::to_bits)
        .chain([solid.accepted_steps()])
        .collect()
}

#[test]
fn resonant_force_keeps_the_physical_secular_transient() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mass, omega, force) = (2.0, 3.0, 1.6);
        let pencil = assembly(&[mass * omega * omega], &[mass], 1);
        let mut solid =
            ExactHarmonicSolid::prepare(&pencil, &[omega], 0.17, &[0.0], &[0.0], controls(), cx)
                .unwrap();
        assert_eq!(solid.angular_frequencies_rad_s(), &[omega]);
        let (mut cosine, mut sine) = ([force], [0.0]);
        for _ in 0..40 {
            let step = solid.advance(&cosine, &sine, cx).unwrap();
            let time = step.time_s;
            let (sn, cs) = (omega * time).sin_cos();
            let q = force * time * sn / (2.0 * mass * omega);
            let v = force * (sn + omega * time * cs) / (2.0 * mass * omega);
            near(solid.displacement_m()[0], q, 3e-10);
            near(solid.velocity_mps()[0], v, 3e-10);
            near(step.energy.dissipated_j, 0.0, 1e-14);
            near(
                solid.cumulative_work_j(),
                0.5 * mass * (v * v + omega * omega * q * q),
                3e-9,
            );
            rotate(&solid, &mut cosine, &mut sine, 1);
        }
    });
}

#[test]
fn full_consistent_mass_and_nonproportional_damping_match_the_physical_pencil() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mass = [2.0, 0.6, 0.6, 1.4];
        let stiffness = [12.0, -2.0, -2.0, 9.0];
        let damping = [0.8, 0.15, 0.15, 0.5];
        let force = [C64::new(1.1, 0.4), C64::new(-0.2, 0.7)];
        let omega = 2.3;
        let reference = PreparedHarmonicResponse::try_new(
            &mass,
            &damping,
            &stiffness,
            2,
            omega,
            HarmonicResponseControls {
                maximum_coordinates: 2,
                maximum_payload_bytes: 1 << 20,
                maximum_relative_residual: 1e-12,
            },
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        let pencil = assembly(&stiffness, &mass, 2);
        let viscosity = assembly(&[0.0; 4], &damping, 2).mass;
        let period = std::f64::consts::TAU / omega;
        let mut solid = ExactHarmonicSolid::prepare_viscous(
            &pencil,
            &viscosity,
            &[omega],
            period / 32.0,
            &reference
                .displacement
                .iter()
                .map(|z| z.re)
                .collect::<Vec<_>>(),
            &reference.velocity.iter().map(|z| z.re).collect::<Vec<_>>(),
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([force[0].re, force[1].re], [force[0].im, force[1].im]);
        for _ in 0..32 {
            let step = solid.advance(&cosine, &sine, cx).unwrap();
            for i in 0..2 {
                near(
                    solid.displacement_m()[i],
                    real_at(reference.displacement[i], omega, step.time_s),
                    2e-10,
                );
                near(
                    solid.velocity_mps()[i],
                    real_at(reference.velocity[i], omega, step.time_s),
                    2e-10,
                );
            }
            rotate(&solid, &mut cosine, &mut sine, 1);
        }
        near(
            solid.cumulative_work_j(),
            reference.average_input_power_w * period,
            2e-9,
        );
        near(
            solid.cumulative_dissipated_energy_j(),
            reference.average_viscous_power_w * period,
            2e-9,
        );
    });
}

#[test]
fn beating_work_and_heat_keep_cross_terms_and_are_interval_partition_invariant() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let frequencies = [3.7, 4.2];
        let force = [C64::new(1.0, 0.2), C64::new(0.8, -0.4)];
        let stiffness = 25.0;
        let damping = 0.6;
        let duration = 1.3;
        let response: Vec<_> = frequencies
            .iter()
            .zip(force)
            .map(|(&w, f)| f / C64::new(stiffness - w * w, -damping * w))
            .collect();
        let velocity: Vec<_> = frequencies
            .iter()
            .zip(&response)
            .map(|(&w, &q)| C64::new(w * q.im, -w * q.re))
            .collect();
        let decay = damping / 2.0;
        let natural = (stiffness - decay * decay).sqrt();
        let a = -response.iter().map(|q| q.re).sum::<f64>();
        let b = (-velocity.iter().map(|v| v.re).sum::<f64>() + decay * a) / natural;
        let trajectory = |time: f64| {
            let (sn, cs) = (natural * time).sin_cos();
            let envelope = (-decay * time).exp();
            let mut q = envelope * (a * cs + b * sn);
            let mut v =
                envelope * ((natural * b - decay * a) * cs - (natural * a + decay * b) * sn);
            let mut f = 0.0;
            for j in 0..2 {
                q += real_at(response[j], frequencies[j], time);
                v += real_at(velocity[j], frequencies[j], time);
                f += real_at(force[j], frequencies[j], time);
            }
            (q, v, f)
        };
        // Independent analytic trajectory quadrature, not an energy difference
        // or sum of isolated steady-state cycle averages.
        let slices = 4096;
        let mut work = 0.0;
        let mut heat = 0.0;
        for i in 0..=slices {
            let (_, v, f) = trajectory(duration * i as f64 / slices as f64);
            let weight = if i == 0 || i == slices {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            work += weight * f * v;
            heat += weight * damping * v * v;
        }
        work *= duration / (3 * slices) as f64;
        heat *= duration / (3 * slices) as f64;
        let pencil = assembly(&[stiffness], &[1.0], 1);
        let viscosity = assembly(&[0.0], &[damping], 1).mass;
        let prepare = |dt| {
            ExactHarmonicSolid::prepare_viscous(
                &pencil,
                &viscosity,
                &frequencies,
                dt,
                &[0.0],
                &[0.0],
                controls(),
                cx,
            )
            .unwrap()
        };
        let mut whole = prepare(duration);
        let cosine = [force[0].re, force[1].re];
        let sine = [force[0].im, force[1].im];
        whole.advance(&cosine, &sine, cx).unwrap();
        let (q, v, _) = trajectory(duration);
        near(whole.displacement_m()[0], q, 2e-10);
        near(whole.velocity_mps()[0], v, 2e-10);
        near(whole.cumulative_work_j(), work, 2e-10);
        near(whole.cumulative_dissipated_energy_j(), heat, 2e-10);
        let mut divided = prepare(duration / 13.0);
        let (mut cosine, mut sine) = (cosine, sine);
        for _ in 0..13 {
            divided.advance(&cosine, &sine, cx).unwrap();
            rotate(&divided, &mut cosine, &mut sine, 2);
        }
        near(divided.displacement_m()[0], q, 2e-10);
        near(divided.velocity_mps()[0], v, 2e-10);
        near(divided.cumulative_work_j(), work, 2e-10);
        near(divided.cumulative_dissipated_energy_j(), heat, 2e-10);
    });
}

#[test]
fn coherent_duplicate_tones_cancel_and_dc_matches_the_held_owner() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let pencil = assembly(&[12.0, -2.0, -2.0, 9.0], &[2.0, 0.6, 0.6, 1.4], 2);
        let viscosity = assembly(&[0.0; 4], &[0.8, 0.15, 0.15, 0.5], 2).mass;
        let mut held = ExactLinearSolid::prepare_viscous(
            &pencil,
            &viscosity,
            0.07,
            &[0.03, -0.08],
            &[0.2, 0.1],
            controls(),
            cx,
        )
        .unwrap();
        let mut tones = ExactHarmonicSolid::prepare_viscous(
            &pencil,
            &viscosity,
            &[0.0, 3.1, 3.1],
            0.07,
            &[0.03, -0.08],
            &[0.2, 0.1],
            controls(),
            cx,
        )
        .unwrap();
        // A sine coefficient at DC does not exert any force. Equal-frequency
        // opposite pairs cancel coherently before work or heat is measured.
        let cosine = [0.4, 1.0, -1.0, -0.3, 0.8, -0.8];
        let sine = [5.0, 0.2, -0.2, -9.0, 0.1, -0.1];
        for _ in 0..8 {
            held.advance(&[0.4, -0.3], cx).unwrap();
            tones.advance(&cosine, &sine, cx).unwrap();
            for i in 0..2 {
                near(tones.displacement_m()[i], held.displacement_m()[i], 1e-10);
                near(tones.velocity_mps()[i], held.velocity_mps()[i], 1e-10);
            }
            near(tones.cumulative_work_j(), held.cumulative_work_j(), 1e-10);
            near(
                tones.cumulative_dissipated_energy_j(),
                held.cumulative_dissipated_energy_j(),
                1e-10,
            );
        }
        let rigid = assembly(&[0.0], &[2.0], 1);
        let mut free =
            ExactHarmonicSolid::prepare(&rigid, &[0.0], 0.4, &[0.3], &[0.2], controls(), cx)
                .unwrap();
        free.advance(&[1.4], &[0.0], cx).unwrap();
        near(
            free.displacement_m()[0],
            0.3 + 0.2 * 0.4 + 0.5 * 0.7 * 0.4 * 0.4,
            1e-12,
        );
        near(free.velocity_mps()[0], 0.2 + 0.7 * 0.4, 1e-12);
    });
}

#[test]
fn refusals_preserve_the_accepted_state_and_source_phase() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<ExactLinearSolid>();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let pencil = assembly(&[9.0], &[1.0], 1);
        for frequencies in [&[][..], &[-1.0][..], &[f64::NAN][..]] {
            assert!(
                ExactHarmonicSolid::prepare(
                    &pencil,
                    frequencies,
                    0.2,
                    &[0.0],
                    &[0.0],
                    controls(),
                    cx,
                )
                .is_err()
            );
        }
        let mut small = controls();
        small.maximum_preparation_bytes = 0;
        assert!(matches!(
            ExactHarmonicSolid::prepare(&pencil, &[3.0], 0.2, &[0.0], &[0.0], small, cx,),
            Err(ExactLinearSolidError::Budget { .. })
        ));
        assert!(matches!(
            ExactHarmonicSolid::prepare(&pencil, &[3.0; 127], 0.2, &[0.0], &[0.0], controls(), cx,),
            Err(ExactLinearSolidError::Budget { .. })
        ));
        let mut solid =
            ExactHarmonicSolid::prepare(&pencil, &[3.0], 0.2, &[0.0], &[0.0], controls(), cx)
                .unwrap();
        solid.advance(&[1.0], &[0.0], cx).unwrap();
        let before = state_bits(&solid);
        let next_phase = solid.advance_force_phase(0, 1.0, 0.0).unwrap();
        for (cosine, sine) in [(&[][..], &[0.0][..]), (&[f64::NAN][..], &[0.0][..])] {
            assert!(solid.advance(cosine, sine, cx).is_err());
            assert_eq!(state_bits(&solid), before);
        }
        assert!(solid.advance_force_phase(1, 1.0, 0.0).is_err());
        gate.request();
        assert!(matches!(
            solid.advance(&[1.0], &[0.0], cx),
            Err(ExactLinearSolidError::Cancelled)
        ));
        assert_eq!(state_bits(&solid), before);
        assert_eq!(solid.advance_force_phase(0, 1.0, 0.0).unwrap(), next_phase);
    });
}
