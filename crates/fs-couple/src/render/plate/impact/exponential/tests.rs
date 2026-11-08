use super::*;
use crate::modal_acoustic_time::ModalAcousticState;
use super::super::{linear::wire::WireSpan, string::StringStretching};

fn config(dt_s: f64) -> ImpactConfig {
    ImpactConfig { dt_s, max_steps: 10_000, maximum_energy_j: 1e3,
        energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-9,
        maximum_generalized_force: 1e4 }
}
fn linear(omega: f64, q: f64) -> ImpactBody {
    ImpactBody { potential: BodyPotential::Linear(vec![omega]),
        initial: vec![ModalAcousticState {
            displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: 0.0,
        }], damping_per_s: vec![0.0] }
}
fn string(amplitude: f64) -> ImpactBody {
    WireSpan { endpoints_m: [[0.0, 0.0], [1.0, 0.0]], linear_density_kg_m: 0.02,
        tension_n: 30.0, bending_stiffness_n_m2: 0.001, damping_per_s: vec![0.0] }
        .stretching_body(vec![ModalAcousticState {
            displacement_m_sqrt_kg: amplitude, velocity_m_sqrt_kg_per_s: 0.0,
        }], StringStretching { axial_rigidity_n: 2000.0, maximum_slope: 0.2 }).unwrap()
}

#[test]
fn physical_bending_modes_keep_their_continuous_phase_and_held_force_work() {
    let gate = CancelGate::new_clock_free();
    let (count, length, mass, tension, ei): (usize, f64, f64, f64, f64) =
        (24, 0.65, 0.003, 80.0, 0.006);
    let wire = WireSpan { endpoints_m: [[0.0, 0.0], [length, 0.0]],
        linear_density_kg_m: mass, tension_n: tension,
        bending_stiffness_n_m2: ei, damping_per_s: vec![0.0; count] };
    let amplitude = 1e-8;
    let mut initial = vec![ModalAcousticState::default(); count];
    initial[count - 1].displacement_m_sqrt_kg = amplitude;
    let body = wire.stretching_body(initial, StringStretching {
        axial_rigidity_n: 0.0, maximum_slope: 0.2,
    }).unwrap();
    let clock = config(1.0 / 48_000.0);
    assert!(matches!(ImpactSystem::new(vec![body.clone()], vec![], vec![], vec![], clock),
        Err(ImpactError::TimeResolution(_))));
    let mut run = ExponentialImpactSystem::new(vec![body], vec![], vec![], clock, &gate).unwrap();
    // Independently substitute sin(k*x) into mu*w_tt = T*w_xx - EI*w_xxxx.
    let k = count as f64 * core::f64::consts::PI / length;
    let omega = ((tension * k * k + ei * k.powi(4)) / mass).sqrt();
    let force = 0.001;
    let mut input = vec![0.0; count]; input[count - 1] = force;
    let mut supplied = 0.0;
    for tick in 1..=600 {
        let frame = run.step(&input, &gate).unwrap(); supplied += frame.supplied_work_j;
        let t = tick as f64 * clock.dt_s;
        let (sine, cosine) = (omega * t).sin_cos();
        let q = amplitude * cosine + force / (omega * omega) * (1.0 - cosine);
        let p = -amplitude * omega * sine + force / omega * sine;
        let scale = omega * amplitude + force / omega;
        assert!(((run.state()[2 * count - 2] - q) * omega).hypot(
            run.state()[2 * count - 1] - p) < 3e-9 * scale);
        assert!(run.state()[..2 * count - 2].iter().all(|value| *value == 0.0));
    }
    let expected_work = force * (run.state()[2 * count - 2] - amplitude);
    assert!((supplied - expected_work).abs() < 1e-18);
}

#[test]
fn volume_coupling_and_spatial_drag_share_one_exact_damped_spectrum() {
    let gate = CancelGate::new_clock_free();
    let omega = core::f64::consts::TAU * 4200.0;
    let stiffness = 0.35 * omega * omega;
    let drag = 31.0;
    let amplitude = 1e-4;
    let clock = config(1.0 / 48_000.0);
    let mut run = ExponentialImpactSystem::new_with_dampers(
        vec![linear(omega, amplitude), linear(omega, -amplitude)], vec![],
        vec![VolumeSpring { bulk_modulus_pa: stiffness, volume_m3: 1.0, areas: vec![1.0, -1.0] }],
        vec![ViscousDamper { weights: vec![1.0, -1.0], damping_n_s_m: drag }], clock, &gate,
    ).unwrap();
    // Antisymmetric coordinate obeys q''+2*c*q'+(w²+2*k)*q=0.
    let frequency_squared = omega * omega + 2.0 * stiffness;
    let damped = (frequency_squared - drag * drag).sqrt();
    let initial_energy = frequency_squared * amplitude * amplitude;
    let mut heat = 0.0;
    for tick in 1..=257 {
        let frame = run.step(&[0.0; 2], &gate).unwrap();
        heat += frame.dissipated_energy_j;
        let t = tick as f64 * clock.dt_s;
        let (sine, cosine) = (damped * t).sin_cos();
        let q = amplitude * (-drag * t).exp() * (cosine + drag / damped * sine);
        let p = -amplitude * (-drag * t).exp() * frequency_squared / damped * sine;
        assert!(((run.state()[0] - q) * omega).hypot(run.state()[1] - p) < 3e-9 * amplitude * omega);
        assert!((run.state()[0] + run.state()[2]).abs() < 1e-12 * amplitude);
        assert!((run.state()[1] + run.state()[3]).abs() < 1e-12 * amplitude * omega);
        assert!(frame.dissipated_energy_j >= 0.0);
    }
    assert!((run.stored_energy_j() + heat - initial_energy).abs() < 1e-9 * initial_energy);
}

#[test]
fn nonlinear_stretching_keeps_energy_and_converges_to_the_duffing_period() {
    let gate = CancelGate::new_clock_free();
    let amplitude: f64 = 0.003;
    let k = core::f64::consts::PI;
    let omega_squared = (30.0 * k * k + 0.001 * k.powi(4)) / 0.02;
    // The sine-field integral gives beta=EA*k^4/(2*mu^2*L).
    let beta = 2000.0 * k.powi(4) / (2.0 * 0.02_f64.powi(2));
    // Independent energy integral for q''+w²*q+beta*q³=0, q(0)=A.
    let cells = 1024;
    let step = core::f64::consts::FRAC_PI_2 / cells as f64;
    let integrand = |theta: f64|
        (omega_squared + 0.5 * beta * amplitude * amplitude * (1.0 + theta.sin().powi(2))).sqrt().recip();
    let mut integral = integrand(0.0) + integrand(core::f64::consts::FRAC_PI_2);
    for j in 1..cells { integral += (if j % 2 == 0 { 2.0 } else { 4.0 }) * integrand(j as f64 * step); }
    let quarter_period = integral * step / 3.0;
    let mut errors = Vec::new();
    for divisions in [16, 32, 64] {
        let dt = quarter_period / divisions as f64;
        let mut run = ExponentialImpactSystem::new(vec![string(amplitude)], vec![], vec![], config(dt), &gate).unwrap();
        let energy = run.stored_energy_j();
        let mut previous = amplitude;
        let mut crossing = None;
        for tick in 1..=2 * divisions {
            run.step(&[0.0], &gate).unwrap();
            let q = run.state()[0];
            assert!((run.stored_energy_j() - energy).abs() < 1e-9 * energy);
            if q <= 0.0 {
                crossing = Some(((tick - 1) as f64 + previous / (previous - q)) * dt);
                break;
            }
            previous = q;
        }
        errors.push((crossing.unwrap() - quarter_period).abs());
    }
    assert!(errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0, "{errors:?}");
    assert!(errors[2] < 3e-5 * quarter_period, "{errors:?}");
}

#[test]
fn rejected_nonlinear_step_and_cancellation_preserve_accepted_state_and_retry() {
    let gate = CancelGate::new_clock_free();
    let make = || ExponentialImpactSystem::new(vec![string(0.003)], vec![], vec![], config(1.0 / 2000.0), &gate).unwrap();
    let mut run = make();
    let before = run.state().to_vec(); let energy = run.stored_energy_j();
    run.set_iteration_limit(0).unwrap();
    assert!(run.step(&[0.0], &gate).is_err());
    assert_eq!(run.state(), before); assert_eq!(run.stored_energy_j(), energy); assert_eq!(run.samples(), 0);
    let stop = CancelGate::new_clock_free(); stop.request();
    assert!(matches!(run.step(&[0.0], &stop), Err(ImpactError::Cancelled)));
    run.set_iteration_limit(50).unwrap();
    let mut fresh = make();
    let actual = run.step(&[0.0], &gate).unwrap();
    let expected = fresh.step(&[0.0], &gate).unwrap();
    assert_eq!(run.state(), fresh.state());
    assert_eq!(actual.stored_energy_j, expected.stored_energy_j);
    assert_eq!(actual.balance_residual_j, expected.balance_residual_j);
}

#[test]
fn complete_coupled_band_capacity_and_unhandled_loss_refuse() {
    let gate = CancelGate::new_clock_free();
    let clock = config(1.0 / 48_000.0);
    let too_fast = (0.95 * core::f64::consts::PI / clock.dt_s).powi(2);
    assert!(ExponentialImpactSystem::new(vec![linear(1.0, 0.0)], vec![],
        vec![VolumeSpring { bulk_modulus_pa: too_fast, volume_m3: 1.0, areas: vec![1.0] }],
        clock, &gate).is_err());
    assert!(ExponentialImpactSystem::new(vec![linear(1.0, 0.0); 52], vec![], vec![], clock, &gate).is_err());
    let contact = Obstacle::new(vec![1.0], 1, 1, vec![0.0], vec![1.0], 100.0, 1.5,
        "explicit elastic contact with nonlinear drag".into()).unwrap().with_internal_loss(0.1).unwrap();
    assert!(ExponentialImpactSystem::new(vec![linear(1.0, 0.0)], vec![contact], vec![], clock, &gate).is_err());
}


#[test]
fn elastic_contact_remainder_converges_without_losing_contact_energy() {
    let gate = CancelGate::new_clock_free();
    let omega = core::f64::consts::TAU * 2000.0;
    let stiffness = 0.4 * omega * omega;
    let gap = -2e-5;
    let horizon = 0.0013;
    let equilibrium = stiffness * gap / (omega * omega + stiffness);
    let loaded = (omega * omega + stiffness).sqrt();
    let expected_q = equilibrium * (1.0 - (loaded * horizon).cos());
    let expected_p = equilibrium * loaded * (loaded * horizon).sin();
    let mut errors = Vec::new();
    for ticks in [64, 128, 256] {
        let contact = Obstacle::new(vec![1.0], 1, 1, vec![gap], vec![1.0],
            stiffness, 1.0, "declared precompressed elastic contact".into()).unwrap();
        let mut run = ExponentialImpactSystem::new(vec![linear(omega, 0.0)], vec![contact],
            vec![], config(horizon / ticks as f64), &gate).unwrap();
        let initial_energy = 0.5 * stiffness * gap * gap;
        for _ in 0..ticks {
            let frame = run.step(&[0.0], &gate).unwrap();
            assert!(run.state()[0] > gap); // the independent active-branch solution applies
            assert!((frame.stored_energy_j - initial_energy).abs() < 1e-9 * initial_energy);
        }
        errors.push(((run.state()[0] - expected_q) * loaded).hypot(run.state()[1] - expected_p)
            / (equilibrium.abs() * loaded));
    }
    assert!(errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5, "{errors:?}");
    assert!(errors[2] < 0.003, "{errors:?}");
}
