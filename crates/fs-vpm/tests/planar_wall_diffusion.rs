//! G0/G1/G3/G4: stress-free-wall PSE against independent Dirichlet heat
//! solutions, exact self-image decay and the unequal-area two-site generator.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::{self, PseDiffusion, PseKernel};
use fs_vpm::planar_wall::{self, PlanarWall2};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, VpmRun, total_circulation};
use std::f64::consts::PI;
use std::mem::size_of;

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5354_5245_5353,
                kernel_id: 0x5053_4557,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn budget(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

fn wall() -> PlanarWall2 {
    PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap()
}

fn diffuse(
    p: &[VortexParticle],
    areas: &[f64],
    duration: f64,
    steps: usize,
    law: PseDiffusion,
) -> VpmRun {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        planar_wall::diffuse_with_cx(
            cx,
            p,
            areas,
            duration / steps as f64,
            steps,
            law,
            wall(),
            budget(p.len(), steps),
        )
    })
    .unwrap()
}

fn coefficient(law: PseDiffusion, distance: f64) -> f64 {
    let q = distance / law.kernel_width_m();
    let (normalization, shape) = match law.kernel() {
        PseKernel::Gaussian => (4.0 / PI, (-q * q).exp()),
        PseKernel::CompactWendland => (
            1008.0 / (5.0 * PI),
            if q >= 1.0 {
                0.0
            } else {
                (1.0 - q).powi(4) * (1.0 + 4.0 * q)
            },
        ),
    };
    normalization * law.kinematic_viscosity_m2_s() * shape / law.kernel_width_m().powi(4)
}

fn laws(nu: f64, width: f64) -> [PseDiffusion; 2] {
    [
        PseDiffusion::new(nu, width).unwrap(),
        PseDiffusion::compact_wendland(nu, width).unwrap(),
    ]
}

fn enstrophy(p: &[VortexParticle], areas: &[f64]) -> f64 {
    p.iter()
        .zip(areas)
        .map(|(p, area)| 0.5 * p.circulation * p.circulation / area)
        .sum()
}

fn same_bits(a: &[VortexParticle], b: &[VortexParticle]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

#[test]
fn g0_self_image_diffuses_bulk_circulation_and_counts_its_work() {
    for law in laws(0.12, 1.1) {
        for gamma in [-0.7, 0.7] {
            let p = [VortexParticle::new([0.4, 0.15], gamma)];
            let area = 0.6;
            let duration = 0.35;
            let result = diffuse(&p, &[area], duration, 7, law);
            let rate = 2.0 * area * coefficient(law, 0.3);
            let exact = gamma * (-rate * duration).exp();
            assert!((result.particles[0].circulation - exact).abs() < 3e-14);
            assert_eq!(result.particles[0].pos, p[0].pos);
            assert!(result.particles[0].circulation.abs() < gamma.abs());
            assert_eq!(result.pair_evaluations, 14);
            assert_eq!(result.work_units, 22);
            assert_eq!(result.peak_live_bytes, size_of::<VortexParticle>());
        }
    }
}

fn two_site_exact(p: &[VortexParticle; 2], areas: [f64; 2], law: PseDiffusion, t: f64) -> [f64; 2] {
    let dx = p[0].pos[0] - p[1].pos[0];
    let real = coefficient(law, dx.hypot(p[0].pos[1] - p[1].pos[1]));
    let image = coefficient(law, dx.hypot(p[0].pos[1] + p[1].pos[1]));
    let a = -areas[1] * (real + image) - 2.0 * areas[0] * coefficient(law, 2.0 * p[0].pos[1]);
    let d = -areas[0] * (real + image) - 2.0 * areas[1] * coefficient(law, 2.0 * p[1].pos[1]);
    let b = areas[0] * (real - image);
    let c = areas[1] * (real - image);
    // Closed-form exponential of the independently assembled odd-extension
    // 2x2 generator, not a composition of the implementation's pair maps.
    let center = 0.5 * (a + d);
    let difference = 0.5 * (a - d);
    let q = (difference * difference + b * c).sqrt();
    let scale = (center * t).exp();
    let co = (q * t).cosh();
    let si = (q * t).sinh() / q;
    [
        scale * ((co + si * difference) * p[0].circulation + si * b * p[1].circulation),
        scale * (si * c * p[0].circulation + (co - si * difference) * p[1].circulation),
    ]
}

#[test]
fn g1_unequal_area_operator_has_second_order_time_convergence_and_dissipates() {
    let p = [
        VortexParticle::new([-0.3, 0.12], 0.6),
        VortexParticle::new([0.4, 0.38], 0.1),
    ];
    let areas = [0.2, 0.7];
    for law in laws(0.2, 1.4) {
        let exact = two_site_exact(&p, areas, law, 0.45);
        let mut errors = Vec::new();
        for steps in [1, 2, 4] {
            let run = diffuse(&p, &areas, 0.45, steps, law);
            errors.push(
                run.particles
                    .iter()
                    .zip(exact)
                    .map(|(p, exact)| (p.circulation - exact).abs())
                    .fold(0.0_f64, f64::max),
            );
            assert!(run.particles.iter().all(|p| p.circulation >= 0.0));
            assert!(total_circulation(&run.particles) < total_circulation(&p));
            assert!(enstrophy(&run.particles, &areas) < enstrophy(&p, &areas));
        }
        assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
        let signed = [
            VortexParticle::new(p[0].pos, -0.5),
            VortexParticle::new(p[1].pos, 0.8),
        ];
        let stiff = diffuse(&signed, &areas, 1000.0, 1, law);
        assert!(enstrophy(&stiff.particles, &areas) <= enstrophy(&signed, &areas));
        for (p, area) in stiff.particles.iter().zip(areas) {
            assert!((-2.5..=0.8 / 0.7).contains(&(p.circulation / area)));
        }
    }
}

fn heat_error(n: i32) -> f64 {
    let dx = 2.4 / f64::from(n);
    let nu = 0.13;
    let duration = 0.12;
    let omega = |x: f64, y: f64, t: f64| {
        let width2 = 0.18 + 4.0 * nu * t;
        ((-(x * x + (y - 0.5).powi(2)) / width2).exp()
            - (-(x * x + (y + 0.5).powi(2)) / width2).exp())
            / (PI * width2)
    };
    let mut p = Vec::new();
    for j in 0..n {
        for i in -n..=n {
            let pos = [f64::from(i) * dx, (f64::from(j) + 0.5) * dx];
            p.push(VortexParticle::new(
                pos,
                dx * dx * omega(pos[0], pos[1], 0.0),
            ));
        }
    }
    let areas = vec![dx * dx; p.len()];
    let run = diffuse(
        &p,
        &areas,
        duration,
        8,
        PseDiffusion::new(nu, 1.3 * dx).unwrap(),
    );
    let mut error = 0.0;
    let mut norm = 0.0;
    for particle in &run.particles {
        let exact = dx * dx * omega(particle.pos[0], particle.pos[1], duration);
        error += (particle.circulation - exact).abs();
        norm += exact.abs();
        assert!(particle.circulation >= 0.0);
    }
    assert!(total_circulation(&run.particles) < total_circulation(&p));
    error / norm
}

#[test]
fn g1_odd_gaussian_heat_solution_improves_with_cloud_and_kernel_refinement() {
    let coarse = heat_error(9);
    let fine = heat_error(18);
    assert!(fine < 0.32 * coarse, "coarse={coarse}, fine={fine}");
    assert!(fine < 0.020, "fine={fine}");
}

#[test]
fn g0_near_wall_kernel_difference_retains_weak_exchange() {
    let h = 1e-10;
    let p = [
        VortexParticle::new([0.0, h], 0.0),
        VortexParticle::new([0.5, h], 1.0),
    ];
    let areas = [0.25, 0.75];
    let dt = 1e-6;
    for law in laws(0.1, 1.0) {
        let derivative = match law.kernel() {
            PseKernel::Gaussian => 4.0 / PI * 0.1 * (-0.25_f64).exp(),
            PseKernel::CompactWendland => 1008.0 / (5.0 * PI) * 0.1 * 10.0 * 0.5_f64.powi(3),
        };
        let leading_transfer = dt * areas[0] * derivative * 4.0 * h * h;
        let run = diffuse(&p, &areas, dt, 1, law);
        let actual = run.particles[0].circulation;
        assert!(actual > 0.0);
        assert!((actual / leading_transfer - 1.0).abs() < 1e-3);
    }
}

#[test]
fn g1_viscous_image_vortex_has_analytic_decay_and_second_order_translation() {
    let h = 0.3;
    let core = 0.08;
    let gamma = 0.8;
    let area = 0.8;
    let p = [VortexParticle::new([-0.2, h], gamma)];
    let law = PseDiffusion::new(0.16, 0.9).unwrap();
    let duration = 1.2;
    let rate = 2.0 * area * coefficient(law, 2.0 * h);
    let velocity_per_gamma = h / (PI * (4.0 * h * h + core * core));
    let exact_x = p[0].pos[0] + velocity_per_gamma * gamma * -(-rate * duration).exp_m1() / rate;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let run = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
            planar_wall::simulate_viscous_with_cx(
                cx,
                &p,
                &[area],
                duration / steps as f64,
                steps,
                core,
                law,
                wall(),
                budget(1, steps),
            )
        })
        .unwrap();
        assert!((run.particles[0].pos[1] - h).abs() < 1e-15);
        assert!((run.particles[0].circulation - gamma * (-rate * duration).exp()).abs() < 2e-14);
        assert_eq!(run.pair_evaluations, 12 * steps as u64);
        errors.push((run.particles[0].pos[0] - exact_x).abs());
    }
    assert!(errors[1] < 0.26 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.26 * errors[1], "{errors:?}");
}

#[test]
fn g3_zero_viscosity_and_inactive_compact_images_retain_original_owners() {
    let p = [
        VortexParticle::new([-0.2, 2.0], 0.2),
        VortexParticle::new([0.3, 2.1], -0.1),
    ];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let old =
            planar_wall::simulate_with_cx(cx, &p, 0.01, 3, 0.1, wall(), budget(2, 3)).unwrap();
        let zero = planar_wall::simulate_viscous_with_cx(
            cx,
            &p,
            &[0.3, 0.7],
            0.01,
            3,
            0.1,
            PseDiffusion::new(0.0, 1.0).unwrap(),
            wall(),
            budget(2, 3),
        )
        .unwrap();
        assert_eq!(old, zero);
        let law = PseDiffusion::compact_wendland(0.1, 1.0).unwrap();
        let free =
            diffusion::diffuse_with_cx(cx, &p, &[0.3, 0.7], 0.1, 3, law, budget(2, 3)).unwrap();
        let bounded =
            planar_wall::diffuse_with_cx(cx, &p, &[0.3, 0.7], 0.1, 3, law, wall(), budget(2, 3))
                .unwrap();
        same_bits(&free.particles, &bounded.particles);
    });
}

#[test]
fn g4_wall_domains_work_caps_and_mid_sweep_refusals_publish_no_state() {
    let p = [
        VortexParticle::new([-0.3, 0.2], 0.2),
        VortexParticle::new([0.4, 0.5], -0.1),
    ];
    let saved = p;
    let law = PseDiffusion::new(0.1, 1.0).unwrap();
    let active = CancelGate::new_clock_free();
    with_cx(&active, Budget::INFINITE, |cx| {
        let refused = planar_wall::diffuse_with_cx(
            cx,
            &p,
            &[0.3, 0.7],
            0.1,
            1,
            law,
            wall(),
            VpmBudget::new(2, 1, 7, usize::MAX),
        );
        assert_eq!(
            refused,
            Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: 8,
                maximum: 7,
            })
        );
        for y in [0.0, -0.1] {
            let outside = [VortexParticle::new([0.0, y], 0.2)];
            assert!(
                planar_wall::diffuse_with_cx(
                    cx,
                    &outside,
                    &[1.0],
                    0.0,
                    0,
                    law,
                    wall(),
                    budget(1, 0),
                )
                .is_err()
            );
        }
        assert!(
            planar_wall::simulate_viscous_with_cx(
                cx,
                &p,
                &[0.3, 0.7],
                -0.1,
                1,
                0.1,
                law,
                wall(),
                budget(2, 1),
            )
            .is_err()
        );
    });
    let cloud: Vec<_> = (0..17)
        .map(|i| VortexParticle::new([0.05 * f64::from(i), 0.2], f64::from(i) - 8.0))
        .collect();
    let original = cloud.clone();
    let refused = with_cx(&active, Budget::INFINITE.with_poll_quota(4), |cx| {
        planar_wall::diffuse_with_cx(cx, &cloud, &[0.1; 17], 0.1, 1, law, wall(), budget(17, 1))
    });
    assert!(matches!(
        refused,
        Err(VpmError::ExecutionBudget(
            BudgetRefusal::PollsExhausted { .. }
        ))
    ));
    same_bits(&cloud, &original);
    same_bits(&p, &saved);
}
