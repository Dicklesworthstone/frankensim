use super::*;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_feec::fixtures::single_tet;
use fs_math::c64::C64;
use fs_vfit::mimo::delayed::DelayedMimoLoadLimits;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::mimo::{DelayedMimoModel, MimoModel};
use fs_vfit::model::PoleTerm;
use fs_vfit::passivity::CertificateClass;

mod quadratic;

const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 4.0,
};
const M: f64 = 1.0 / (60.0 * 2.0 * 16.0);
const K: f64 = 0.25;
const C: f64 = 0.02 / 12.0;
const B: f64 = 1.0 / 6.0;
const R: f64 = 6.0;
const D: f64 = 3.0;
const A: f64 = 5.0;

fn constant<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}
fn limits() -> CausalTimeLimits {
    CausalTimeLimits {
        interior: ExactTimeLimits::default(),
        maximum_exterior_ports: 2,
        interface_relative_tolerance: 1e-9,
    }
}
fn load_limits(dt: f64) -> MimoLoadLimits {
    MimoLoadLimits {
        band_rad_s: [0.1, core::f64::consts::PI / dt],
        sampled_passivity_points: 65,
        max_states: 8,
        max_memory_bytes: 1 << 20,
    }
}
fn scalar_load(dt: f64) -> CausalExteriorLoad {
    CausalExteriorLoad::Rational(
        PreparedMimoLoad::new(
            &MimoModel {
                ports: 1,
                poles: vec![PoleTerm::Real {
                    pole: -A,
                    residue: 0.0,
                }],
                residues: vec![vec![C64::from_re(R)]],
                d: vec![D],
                e: vec![0.0],
                reciprocal: true,
            },
            dt,
            load_limits(dt),
        )
        .unwrap(),
    )
}
fn scalar_ports() -> [ExactVelocityPort; 2] {
    [
        ExactVelocityPort {
            profile: vec![([0, 1, 2], [1.0; 3])],
        },
        ExactVelocityPort {
            profile: vec![([0, 1, 3], [1.0; 3])],
        },
    ]
}
fn walls() -> [ResistiveBoundary; 1] {
    [ResistiveBoundary {
        faces: vec![[0, 2, 3]],
        admittance_m_pa_s: 0.02,
    }]
}
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}"
    );
}

// Independent one-free-coordinate tetrahedron. For phi0=1-x-y-z after
// releasing [1,2,3], original continuum moments give M=1/(60 rho c²),
// K=1/(2 rho), C=Y/12, b=1/6, and normalized volume-source load 1/4.
// This oracle uses closed-form damped oscillator evolution, not an FE matrix,
// original propagator, port-response cache or rational runtime.
fn oscillator(q: f64, p: f64, force: f64, dt: f64) -> (f64, f64) {
    let alpha = C / (2.0 * M);
    let omega = (K / M - alpha * alpha).sqrt();
    let decay = (-alpha * dt).exp();
    let (sin, cos) = (omega * dt).sin_cos();
    let x = q - force / K;
    (
        force / K + decay * (x * cos + (p + alpha * x) * sin / omega),
        decay * (p * cos - (alpha * p + K / M * x) * sin / omega),
    )
}
fn scalar_step(
    q: f64,
    p: f64,
    x: f64,
    feed: f64,
    volume: f64,
    dt: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let base_force = B * feed + 0.25 * volume;
    let (base_q, _) = oscillator(q, p, base_force, dt);
    let (unit_q, _) = oscillator(0.0, 0.0, B, dt);
    let phi = -(-A * dt).exp_m1() / (A * dt);
    let exterior_history = R * x * phi;
    let direct = D + R / A * (1.0 - phi);
    let velocity = (B * (base_q - q) / dt - exterior_history) / (B * unit_q / dt + direct);
    let force = base_force - B * velocity;
    let (next_q, next_p) = oscillator(q, p, force, dt);
    let next_x = (-A * dt).exp() * x - (-A * dt).exp_m1() / A * velocity;
    // Independent physical Robin heat, integrated on the same exact analytic
    // trajectory. Simpson's 512 subintervals are well below the tested tolerance.
    let panels = 512;
    let h = dt / f64::from(panels);
    let mut heat = p * p + next_p * next_p;
    for i in 1..panels {
        let (_, pressure) = oscillator(q, p, force, f64::from(i) * h);
        heat += if i % 2 == 0 { 2.0 } else { 4.0 } * pressure * pressure;
    }
    heat *= C * h / 3.0;
    (
        next_q,
        next_p,
        next_x,
        velocity,
        exterior_history + direct * velocity,
        heat,
    )
}

#[test]
fn exact_original_tetrahedron_and_causal_pole_match_analytic_work_and_history() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    let volume = model
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let dt = 0.02;
    let mut coupled = CausalInteriorTime::new(
        &model,
        dt,
        &walls(),
        &scalar_ports(),
        &[1],
        &[&volume],
        &[[1, 2, 3]],
        scalar_load(dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    coupled
        .initialize_field(&[0.003, 0.0, 0.0, 0.0], &[0.02, 0.0, 0.0, 0.0], &mut || {
            false
        })
        .unwrap();
    let (mut q, mut p, mut x) = (0.003, 0.02, 0.0);
    let mut volume_extracted = false;
    let mut returning_work = false;
    for step in 0..80 {
        let feed = if step < 32 { 0.02 } else { 0.0 };
        let flow = if step < 16 { -0.001 } else { 0.0 };
        let (qn, pn, xn, velocity, force, heat) = scalar_step(q, p, x, feed, flow, dt);
        let result = coupled.step(&[feed], &[flow], &mut || false).unwrap();
        close(coupled.state().pressure_potential_pa_s()[0], qn, 3e-10);
        close(coupled.state().pressure_pa()[0], pn, 5e-9);
        for i in 1..4 {
            assert_eq!(coupled.state().pressure_pa()[i], 0.0);
            assert_eq!(coupled.state().pressure_potential_pa_s()[i], 0.0);
        }
        let CausalExteriorLoad::Rational(load) = coupled.exterior() else {
            panic!()
        };
        close(load.pole_states()[0].re, xn, 3e-11);
        assert_eq!(load.pole_states()[0].im, 0.0);
        close(result.outward_velocity_m_s[0], velocity, 3e-10);
        close(result.interior_mean_force_n[0], force, 3e-9);
        close(result.exterior_mean_force_n[0], force, 3e-9);
        let feed_work = feed * B * (qn - q);
        let volume_work = 0.25 * flow * (qn - q);
        close(result.work.boundary_drive_work_j, feed_work, 2e-12);
        close(result.work.volume_source_work_j, volume_work, 2e-12);
        close(result.work.source_work_j, feed_work + volume_work, 2e-12);
        close(result.work.boundary_dissipation_j, heat, 2e-12);
        close(
            result.work.exterior_exchange_work_j,
            dt * velocity * force,
            2e-12,
        );
        close(
            result.work.final_fluid_energy_j,
            0.5 * (K * qn * qn + M * pn * pn),
            2e-12,
        );
        close(result.work.interface_work_defect_j, 0.0, 2e-12);
        close(result.work.balance_residual_j(), 0.0, 2e-12);
        volume_extracted |= result.work.volume_source_work_j < 0.0;
        returning_work |= result.work.exterior_exchange_work_j < 0.0;
        (q, p, x) = (qn, pn, xn);
    }
    assert!(volume_extracted);
    // This passive R/(s+a)+D source may return stored energy transiently;
    // the signed report is deliberately not a nonnegative heat account.
    assert!(returning_work);
    assert_eq!(coupled.accepted_steps(), 80);
    assert_eq!(
        coupled.exterior().reports().1.class,
        CertificateClass::GridOnly
    );
}

#[test]
fn held_flow_exchange_converges_to_the_continuously_coupled_original_system() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    let feed = 0.02;
    let total = 0.4;
    // Independent continuous three-state ODE. D>0 permits algebraic elimination
    // only in this oracle: v=(b p-R x)/D. Production never inverts D and admits
    // singular instantaneous terms through its H+D_mean closure.
    let derivative = |y: [f64; 3]| {
        let velocity = (B * y[1] - R * y[2]) / D;
        [
            y[1],
            (B * feed - B * velocity - K * y[0] - C * y[1]) / M,
            -A * y[2] + velocity,
        ]
    };
    let add = |y: [f64; 3], k: [f64; 3], scale: f64| std::array::from_fn(|i| y[i] + scale * k[i]);
    let mut reference = [0.0; 3];
    let h = total / 16_384.0;
    for _ in 0..16_384 {
        let k1 = derivative(reference);
        let k2 = derivative(add(reference, k1, 0.5 * h));
        let k3 = derivative(add(reference, k2, 0.5 * h));
        let k4 = derivative(add(reference, k3, h));
        reference = std::array::from_fn(|i| {
            reference[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = total / f64::from(count);
        let mut coupled = CausalInteriorTime::new(
            &model,
            dt,
            &walls(),
            &scalar_ports(),
            &[1],
            &[],
            &[[1, 2, 3]],
            scalar_load(dt),
            limits(),
            &mut || false,
        )
        .unwrap();
        for _ in 0..count {
            coupled.step(&[feed], &[], &mut || false).unwrap();
        }
        let CausalExteriorLoad::Rational(load) = coupled.exterior() else {
            panic!()
        };
        let dq = coupled.state().pressure_potential_pa_s()[0] - reference[0];
        let dp = coupled.state().pressure_pa()[0] - reference[1];
        let dx = load.pole_states()[0].re - reference[2];
        errors.push((K * dq * dq + M * dp * dp + R * dx * dx).sqrt());
    }
    assert!(errors[0] > 5e-6, "exchange axis must be live: {errors:?}");
    assert!(
        errors[1] < 0.27 * errors[0] && errors[2] < 0.27 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 4.2e-7, "{errors:?}");
}

fn delayed_source(dt: f64) -> DelayedMimoModel {
    DelayedMimoModel {
        delays: vec![0.0, 1.375 * dt, 1.375 * dt, 0.0],
        model: MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real {
                pole: -3.0,
                residue: 0.0,
            }],
            residues: vec![vec![
                C64::from_re(0.4),
                C64::from_re(0.04),
                C64::from_re(0.04),
                C64::from_re(0.3),
            ]],
            d: vec![2.0, 0.2, 0.2, 3.0],
            e: vec![0.0; 4],
            reciprocal: true,
        },
    }
}
fn delayed_load(source: &DelayedMimoModel, dt: f64) -> CausalExteriorLoad {
    CausalExteriorLoad::Delayed(
        PreparedDelayedMimoLoad::new(
            source,
            dt,
            DelayedMimoLoadLimits {
                source: load_limits(dt),
                max_delayed_states: 32,
                max_history_values: 128,
            },
        )
        .unwrap(),
    )
}
// Integral from 0 to t of each independent continuous source unit-step
// response. Differences at actual fractional arrival times reconstruct the
// interval mean without the runtime's split-arrival recurrence or phi functions.
fn step_primitive(source: &MimoModel, entry: usize, time: f64) -> f64 {
    if time <= 0.0 {
        return 0.0;
    }
    let PoleTerm::Real { pole, .. } = source.poles[0] else {
        panic!()
    };
    source.d[entry] * time
        + source.residues[0][entry].re * ((pole * time).exp_m1() / pole - time) / pole
}

#[test]
fn full_p2_field_retains_signed_modes_mutual_reaction_and_fractional_arrivals() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = constant(&geometry);
    let ports = [
        ExactVelocityPort {
            profile: vec![([0, 1, 2], [1.0, 0.0, -0.4])],
        },
        ExactVelocityPort {
            profile: vec![([0, 1, 3], [1.0; 3])],
        },
        ExactVelocityPort {
            profile: vec![([0, 2, 3], [1.0, -1.0, 0.0])],
        },
    ];
    let dt = 0.01;
    let source = delayed_source(dt);
    let mut coupled = CausalInteriorTime::new(
        &model,
        dt,
        &[],
        &ports,
        &[2, 1],
        &[],
        &[],
        delayed_load(&source, dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(coupled.exterior_ports(), [2, 1]);
    assert_eq!(coupled.driven_ports(), [0]);
    assert_eq!(coupled.state().pressure_pa().len(), 10);
    let mut past: Vec<Vec<f64>> = Vec::new();
    let mut signed_mode_live = false;
    let mut delay_live = false;
    for step in 0..48 {
        let feed = if step < 20 {
            0.02 * (0.37 * f64::from(step)).cos()
        } else {
            0.0
        };
        let result = coupled.step(&[feed], &[], &mut || false).unwrap();
        past.push(result.outward_velocity_m_s.clone());
        signed_mode_live |= result.outward_velocity_m_s[0].abs() > 1e-7;
        for a in 0..2 {
            let mut expected = 0.0;
            let mut undelayed = 0.0;
            for b in 0..2 {
                let entry = 2 * a + b;
                for k in 0..=step as usize {
                    let jump = past[k][b] - if k == 0 { 0.0 } else { past[k - 1][b] };
                    let time = (step as usize - k) as f64 * dt;
                    let start = time - source.delays[entry];
                    expected += jump
                        * (step_primitive(&source.model, entry, start + dt)
                            - step_primitive(&source.model, entry, start))
                        / dt;
                    undelayed += jump
                        * (step_primitive(&source.model, entry, time + dt)
                            - step_primitive(&source.model, entry, time))
                        / dt;
                }
            }
            close(result.exterior_mean_force_n[a], expected, 3e-10);
            close(result.interior_mean_force_n[a], expected, 3e-10);
            delay_live |= (expected - undelayed).abs() > 1e-8;
        }
        let scale = result.work.initial_fluid_energy_j
            + result.work.final_fluid_energy_j
            + result.work.source_work_j.abs()
            + result.work.exterior_exchange_work_j.abs();
        close(
            result.work.balance_residual_j(),
            0.0,
            2e-9 * scale.max(1e-12),
        );
    }
    assert!(signed_mode_live && delay_live);
    let CausalExteriorLoad::Delayed(load) = coupled.exterior() else {
        panic!()
    };
    assert_eq!(load.source_model(), &source);
    assert!(load.pole_states().iter().any(|x| x.re != 0.0));
    assert!(load.input_history().iter().any(|x| *x != 0.0));
    assert_eq!(load.reports().0.class, CertificateClass::GridOnly);
    assert_eq!(load.reports().1.class, CertificateClass::GridOnly);
}

#[test]
fn refused_steps_keep_both_original_histories_and_clock_unchanged() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    let dt = 0.02;
    let mut coupled = CausalInteriorTime::new(
        &model,
        dt,
        &walls(),
        &scalar_ports(),
        &[1],
        &[],
        &[[1, 2, 3]],
        scalar_load(dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    coupled.step(&[0.02], &[], &mut || false).unwrap();
    let field = (
        coupled.state().pressure_potential_pa_s().to_vec(),
        coupled.state().pressure_pa().to_vec(),
    );
    let CausalExteriorLoad::Rational(load) = coupled.exterior() else {
        panic!()
    };
    let poles = load.pole_states().to_vec();
    assert!(coupled.step(&[f64::NAN], &[], &mut || false).is_err());
    assert!(coupled.step(&[], &[], &mut || false).is_err());
    assert!(coupled.step(&[0.02], &[0.0], &mut || false).is_err());
    assert!(coupled.step(&[0.02], &[], &mut || true).is_err());
    let mut checkpoints = 0;
    assert!(
        coupled
            .step(&[0.02], &[], &mut || {
                checkpoints += 1;
                checkpoints > 3
            })
            .is_err()
    );
    assert!(
        coupled
            .initialize_field(&field.0, &field.1, &mut || false)
            .is_err()
    );
    // Corrupt only the private test factor to show that a failed original
    // pressure constraint cannot advance even valid exterior pole history.
    coupled.closure = lu(&[0.1], 1).unwrap();
    assert!(coupled.step(&[0.03], &[], &mut || false).is_err());
    assert_eq!(coupled.state().pressure_potential_pa_s(), field.0);
    assert_eq!(coupled.state().pressure_pa(), field.1);
    assert_eq!(coupled.accepted_steps(), 1);
    let CausalExteriorLoad::Rational(load) = coupled.exterior() else {
        panic!()
    };
    assert_eq!(load.pole_states(), poles);
}

#[test]
fn original_clock_port_support_release_and_resource_admission_are_preserved() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    let dt = 0.02;
    let prepare = |ports: &[ExactVelocityPort],
                   selected: &[usize],
                   release: &[[u32; 3]],
                   load: CausalExteriorLoad,
                   budget: CausalTimeLimits| {
        CausalInteriorTime::new(
            &model,
            dt,
            &[],
            ports,
            selected,
            &[],
            release,
            load,
            budget,
            &mut || false,
        )
    };
    assert!(prepare(&scalar_ports(), &[1], &[], scalar_load(dt * 2.0), limits()).is_err());
    assert!(prepare(&scalar_ports(), &[2], &[], scalar_load(dt), limits()).is_err());
    assert!(prepare(&scalar_ports(), &[1, 1], &[], scalar_load(dt), limits()).is_err());
    assert!(
        prepare(
            &scalar_ports(),
            &[1],
            &[[0, 1, 3]],
            scalar_load(dt),
            limits()
        )
        .is_err()
    );
    let mut overlapping = scalar_ports();
    overlapping[0] = overlapping[1].clone();
    assert!(prepare(&overlapping, &[1], &[], scalar_load(dt), limits()).is_err());
    overlapping[1].profile[0].1 = [0.0; 3];
    assert!(prepare(&overlapping, &[1], &[], scalar_load(dt), limits()).is_err());
    let mut budget = limits();
    budget.maximum_exterior_ports = 0;
    assert!(prepare(&scalar_ports(), &[1], &[], scalar_load(dt), budget).is_err());
    budget = limits();
    budget.interior.max_dofs = 3;
    assert!(prepare(&scalar_ports(), &[1], &[], scalar_load(dt), budget).is_err());
    let foreign = constant(&geometry);
    let source = foreign
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    assert!(
        CausalInteriorTime::new(
            &model,
            dt,
            &[],
            &scalar_ports(),
            &[1],
            &[&source],
            &[],
            scalar_load(dt),
            limits(),
            &mut || false
        )
        .is_err()
    );
}
