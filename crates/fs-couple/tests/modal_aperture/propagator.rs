//! Free structural poles, full material damping and nonlinear work on the
//! actual modal slit owner. Analytic controls do not call its exponential.
use super::*;
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_exec::CancelGate;

fn exact(s: ModalApertureSpec, state: Vec<BladeState>) -> ModalAperture {
    ModalAperture::new_with_propagator(s, state, contact(), BladePropagator::ExactZoh).unwrap()
}

fn free_spec(angular_step: f64, damping_ratio: f64) -> ModalApertureSpec {
    let mut s = spec();
    s.modes.truncate(1);
    s.time_step_s = 1.0 / 48_000.0;
    s.max_steps = 10_000;
    let m = &mut s.modes[0];
    let omega = angular_step / s.time_step_s;
    m.stiffness_n_m = m.mass_kg * omega * omega;
    m.damping_n_s_m = 2.0 * m.mass_kg * omega * damping_ratio;
    m.pressure_area_m2 = 0.0;
    s
}

#[test]
fn exact_free_modes_retain_phase_at_audio_rate_and_midpoint_refuses_the_warped_band() {
    let displacement = 1e-6;
    for angular_step in [0.1, 0.5, 1.0, 2.0] {
        let s = free_spec(angular_step, 0.0);
        let old = vec![BladeState {
            displacement_m: displacement,
            velocity_m_s: 0.0,
        }];
        assert!(ModalAperture::new(s.clone(), old.clone(), contact()).is_err());
        let omega = angular_step / s.time_step_s;
        let mut a = exact(s, old);
        let initial = a.stored_energy_j();
        for sample in 1..=2_000 {
            let frame = a.step(drive(0.0), [0.0; 2]).unwrap();
            let phase = sample as f64 * angular_step;
            let measured_phase =
                (-frame.modes[0].velocity_m_s / omega).atan2(frame.modes[0].displacement_m);
            let error = (measured_phase - phase).sin().abs();
            assert!(
                error < 1e-9,
                "omega*dt={angular_step}, sample={sample}, phase error={error}"
            );
            assert!((frame.modes[0].displacement_m - displacement * phase.cos()).abs() < 1e-14);
            assert!((frame.aperture.stored_energy_j - initial).abs() < 1e-11 * initial);
            assert!(frame.work_certificate.residual <= frame.work_certificate.effective_budget);
        }
    }
}

#[test]
fn exact_damped_modes_retain_original_decay_and_independent_material_loss() {
    let s = free_spec(1.0, 0.025);
    let (dt, amplitude) = (s.time_step_s, 1e-6);
    let omega = 1.0 / dt;
    let alpha = 0.025 * omega;
    let damped = (omega * omega - alpha * alpha).sqrt();
    let mut a = exact(
        s,
        vec![BladeState {
            displacement_m: amplitude,
            velocity_m_s: 0.0,
        }],
    );
    let initial = a.stored_energy_j();
    let mut loss = 0.0;
    for sample in 1..=128 {
        let f = a.step(drive(0.0), [0.0; 2]).unwrap();
        let t = sample as f64 * dt;
        let decay = (-alpha * t).exp();
        let q = amplitude * decay * ((damped * t).cos() + alpha / damped * (damped * t).sin());
        let v = -amplitude * decay * omega * omega / damped * (damped * t).sin();
        assert!((f.modes[0].displacement_m - q).abs() < 1e-12 * amplitude * decay);
        assert!((f.modes[0].velocity_m_s - v).abs() < 1e-12 * amplitude * omega * decay);
        assert!(f.aperture.dissipated_energy_j >= 0.0);
        loss += f.aperture.dissipated_energy_j;
    }
    assert!((initial - a.stored_energy_j() - loss).abs() < 1e-12 * initial);
}

#[test]
fn full_nonproportional_damping_matches_independent_analytic_normal_coordinates() {
    // Equal undamped frequencies, distinct physical masses and a known dense
    // damping basis provide an independent exact solution with signed C_ij.
    let mut s = spec();
    s.time_step_s = 1.0 / 48_000.0;
    let omega = 24_000.0;
    let u = [
        [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
        [2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
        [-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    ];
    let rates = [0.02 * omega, 0.08 * omega, 0.15 * omega];
    let roots: Vec<_> = s.modes.iter().map(|m| m.mass_kg.sqrt()).collect();
    for m in &mut s.modes {
        m.stiffness_n_m = m.mass_kg * omega * omega;
        m.damping_n_s_m = 0.0;
        m.pressure_area_m2 = 0.0;
    }
    let mut c = vec![0.0; 9];
    for i in 0..3 {
        for j in i..3 {
            let value =
                roots[i] * roots[j] * (0..3).map(|k| u[i][k] * rates[k] * u[j][k]).sum::<f64>();
            c[3 * i + j] = value;
            c[3 * j + i] = value;
        }
    }
    let q0 = [1e-6, -0.3e-6, 0.8e-6];
    let normal: Vec<_> = (0..3)
        .map(|k| (0..3).map(|i| u[i][k] * roots[i] * q0[i]).sum::<f64>())
        .collect();
    let old = q0
        .map(|q| BladeState {
            displacement_m: q,
            velocity_m_s: 0.0,
        })
        .to_vec();
    let dt = s.time_step_s;
    let mut a = exact(s, old)
        .with_viscosity(c.clone(), &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(a.viscosity_matrix_n_s_m(), Some(c.as_slice()));
    assert_eq!(a.blade_propagator(), BladePropagator::ExactZoh);
    let initial = a.stored_energy_j();
    let mut loss = 0.0;
    for sample in 1..=128 {
        let f = a.step(drive(0.0), [0.0; 2]).unwrap();
        let t = sample as f64 * dt;
        for i in 0..3 {
            let mut q = 0.0;
            let mut v = 0.0;
            for k in 0..3 {
                let alpha = 0.5 * rates[k];
                let wd = (omega * omega - alpha * alpha).sqrt();
                let scale = u[i][k] * normal[k] * (-alpha * t).exp() / roots[i];
                q += scale * ((wd * t).cos() + alpha / wd * (wd * t).sin());
                v -= scale * omega * omega / wd * (wd * t).sin();
            }
            assert!((f.modes[i].displacement_m - q).abs() < 2e-18);
            assert!((f.modes[i].velocity_m_s - v).abs() < 5e-14);
            assert!(f.momentum_residual_n[i].abs() < 1e-11);
        }
        loss += f.aperture.dissipated_energy_j;
        assert!(f.aperture.dissipated_energy_j >= 0.0);
    }
    assert!((initial - a.stored_energy_j() - loss).abs() < 2e-12 * initial);
}

#[test]
fn exact_full_c_closure_reopening_and_signed_force_ports_close_the_same_work_account() {
    let mut s = spec();
    s.time_step_s = 1.0 / 48_000.0;
    s.rest_opening_m = 4e-5;
    for m in &mut s.modes {
        m.damping_n_s_m = 0.0;
    }
    let mut a = exact(s, vec![BladeState::default(); 3])
        .with_viscosity(
            vec![0.03, 0.005, -0.002, 0.005, 0.02, 0.003, -0.002, 0.003, 0.01],
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    a.set_viscous_slit(
        fs_couple::bernoulli_aperture::viscous::ViscousSlit::new(1.81e-5, 0.001).unwrap(),
    )
    .unwrap();
    let mut closed = false;
    let mut reopened = false;
    for sample in 0..1_024 {
        let before = a.modal_states().to_vec();
        let forces = [0.001, -0.0003];
        let f = a
            .step(drive(if sample < 512 { 2_000.0 } else { -100.0 }), forces)
            .unwrap();
        closed |= f.aperture.state.opening_m < 0.0;
        reopened |= closed && sample >= 512 && f.aperture.state.opening_m > 0.0;
        let increments: Vec<_> = before
            .iter()
            .zip(&f.modes)
            .map(|(x, y)| y.displacement_m - x.displacement_m)
            .collect();
        let swept = -a
            .spec()
            .modes
            .iter()
            .zip(&increments)
            .map(|(m, q)| m.pressure_area_m2 * q)
            .sum::<f64>();
        assert!((a.spec().time_step_s * f.aperture.swept_flow_m3_s - swept).abs() < 1e-21);
        let work: f64 = a
            .spec()
            .modes
            .iter()
            .zip(&increments)
            .map(|(m, q)| {
                (m.force_projection[0] * forces[0] + m.force_projection[1] * forces[1]) * q
            })
            .sum();
        assert!((f.aperture.mechanical_work_j - work).abs() < 1e-21);
        assert!(f.momentum_residual_n.iter().all(|r| r.abs() < 1e-8));
        assert!(f.work_certificate.residual <= f.work_certificate.effective_budget);
        assert!(f.flow_certificate.residual <= f.flow_certificate.effective_budget);
    }
    assert!(
        closed && reopened,
        "the original contact must engage and release"
    );
    let state = a.modal_states().to_vec();
    let energy = a.stored_energy_j();
    assert!(a.step(drive(f64::NAN), [0.0; 2]).is_err());
    assert_eq!(a.modal_states(), state);
    assert_eq!(a.stored_energy_j(), energy);
}
