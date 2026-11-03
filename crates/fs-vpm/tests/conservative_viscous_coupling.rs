//! Independent coupled-ODE checks of the native conservative viscous owner.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_vpm::conservative::{simulate_viscous_with_cx, simulate_with_cx};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::{VortexParticle, VpmBudget};

fn context<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        ))
    })
}

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn coupled_unequal_area_pair_matches_independent_closed_form_at_second_order() {
    let cloud = [
        VortexParticle::new([0.5, 0.0], 0.9),
        VortexParticle::new([-0.5, 0.0], -0.3),
    ];
    let areas = [0.7, 1.3];
    let (nu, width, core, duration): (f64, f64, f64, f64) = (0.2, 1.1, 0.15, 0.7);
    let total = 0.6;
    let velocity_kernel = 1.0 / (std::f64::consts::TAU * (1.0 + core * core));
    let omega = total * velocity_kernel;
    // Pair separation and total circulation remain constant in the complete
    // coupled ODE. Its circulation difference relaxes exponentially while
    // the separation rotates, giving this independent analytic trajectory.
    let rate = 4.0 * nu / (std::f64::consts::PI * width.powi(4))
        * (-1.0 / (width * width)).exp()
        * (areas[0] + areas[1]);
    let equilibrium = areas[0] * total / (areas[0] + areas[1]);
    let equilibrium_difference = total - 2.0 * equilibrium;
    let transient_difference = cloud[1].circulation - cloud[0].circulation - equilibrium_difference;
    let theta = omega * duration;
    let (sine, cosine) = theta.sin_cos();
    let decay = (-rate * duration).exp();
    let denominator = rate * rate + omega * omega;
    let integral_cosine = (rate * (1.0 - decay * cosine) + omega * decay * sine) / denominator;
    let integral_sine = (omega * (1.0 - decay * cosine) - rate * decay * sine) / denominator;
    let midpoint = [
        -0.5 * velocity_kernel
            * (equilibrium_difference * (1.0 - cosine) / omega
                + transient_difference * integral_sine),
        0.5 * velocity_kernel
            * (equilibrium_difference * sine / omega + transient_difference * integral_cosine),
    ];
    let expected = [
        [midpoint[0] + 0.5 * cosine, midpoint[1] + 0.5 * sine],
        [midpoint[0] - 0.5 * cosine, midpoint[1] - 0.5 * sine],
    ];
    let gamma = equilibrium + (cloud[0].circulation - equilibrium) * decay;
    context(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let errors = [8, 16, 32].map(|steps| {
            let run = simulate_viscous_with_cx(
                cx,
                &cloud,
                &areas,
                duration / (steps as f64),
                steps,
                core,
                PseDiffusion::new(nu, width).unwrap(),
                caps(2, steps),
            )
            .unwrap();
            close(run.particles[0].circulation, gamma, 3e-15);
            close(run.particles[1].circulation, total - gamma, 3e-15);
            assert_eq!(run.pair_evaluations, 6 * steps as u64);
            run.particles
                .iter()
                .zip(expected)
                .flat_map(|(a, b)| [a.pos[0] - b[0], a.pos[1] - b[1]])
                .map(|e| e * e)
                .sum::<f64>()
                .sqrt()
        });
        assert!(errors[2] < 2e-6, "errors={errors:?}");
        for pair in errors.windows(2) {
            assert!(
                (3.9..4.1).contains(&(pair[0] / pair[1])),
                "errors={errors:?}"
            );
        }
    });
}

#[test]
fn uniform_sampled_vorticity_is_an_exact_diffusion_identity() {
    let cloud = [
        VortexParticle::new([-0.4, 0.2], 0.25),
        VortexParticle::new([0.6, -0.1], 0.5),
        VortexParticle::new([0.1, 0.9], 1.0),
    ];
    context(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let inviscid = simulate_with_cx(cx, &cloud, 0.01, 30, 0.2, caps(3, 30)).unwrap();
        let viscous = simulate_viscous_with_cx(
            cx,
            &cloud,
            &[1.0, 2.0, 4.0],
            0.01,
            30,
            0.2,
            PseDiffusion::new(100.0, 0.5).unwrap(),
            caps(3, 30),
        )
        .unwrap();
        assert_eq!(viscous.particles, inviscid.particles);
    });
}
