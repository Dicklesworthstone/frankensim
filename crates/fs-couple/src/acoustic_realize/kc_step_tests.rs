//! G1/G3 integration of the complete KC image with the original exponential
//! discrete-gradient owner. Manufactured mechanics are numerical references,
//! not measured material or instrument validation.
use super::*;

fn string(frequency_hz: f64, modes: usize) -> PrestressedString {
    let length_m: f64 = 0.5;
    let density: f64 = 0.006;
    PrestressedString {
        length_m,
        tension_n: density * (2.0 * length_m * frequency_hz).powi(2),
        lin_density_kg_m: density,
        axial_stiffness_n: 20_000.0,
        width_m: 0.001,
        n_modes: modes,
        damping_ratio: 0.0,
        // This explicit complete authored loss model makes the analytic
        // fixtures lossless, including the otherwise present cylinder drag.
        rayleigh: Some(RayleighParams { alpha_per_s: 0.0, beta_s: 0.0 }),
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: false,
        tension_loss_factor: None,
        bending_loss_factor: None,
    }
}

fn mass_scale(string: &PrestressedString) -> f64 {
    (string.lin_density_kg_m * string.length_m / 2.0).sqrt()
}

// Independent tension-plus-bending dispersion and physical KC integral.
fn omega(string: &PrestressedString, index: usize) -> f64 {
    let k = index as f64 * core::f64::consts::PI / string.length_m;
    ((string.tension_n * k * k + string.bending_stiffness_n_m2 * k.powi(4))
        / string.lin_density_kg_m).sqrt()
}

fn physical_string_energy(string: &PrestressedString, state: &[f64]) -> f64 {
    let mut energy = 0.0;
    let mut squared_slope = 0.0;
    for mode in 0..string.n_modes {
        let k = (mode + 1) as f64 * core::f64::consts::PI / string.length_m;
        let q = state[2 * mode];
        let p = state[2 * mode + 1];
        energy += 0.5 * (p * p + omega(string, mode + 1).powi(2) * q * q);
        squared_slope += (k * q / mass_scale(string)).powi(2);
    }
    energy + string.axial_stiffness_n * string.length_m * squared_slope.powi(2) / 32.0
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale.max(f64::MIN_POSITIVE),
        "{actual:.17e} != {expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn kc_high_frequency_and_inharmonic_partials_retain_the_linear_phase() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let dt = 1.0 / 48_000.0;
    for (frequency, modes, bending) in [
        (1_000.0, 1, 0.0),
        (5_000.0, 1, 0.0),
        (10_000.0, 1, 0.0),
        (1_000.0, 4, 0.5),
    ] {
        let mut string = string(frequency, modes);
        string.bending_stiffness_n_m2 = bending;
        let mut member = kc_string_member(
            &string, None, None, &[], &gas, 1.0, 1.0, dt,
        ).unwrap();
        assert!(matches!(&member.flow, KcFlow::Exact(_)));
        let amplitudes: Vec<_> = (0..modes)
            .map(|mode| mass_scale(&string) * 1e-9 / (mode + 1) as f64).collect();
        for (mode, &amplitude) in amplitudes.iter().enumerate() {
            member.x[2 * mode] = amplitude;
        }
        if bending > 0.0 {
            assert!(omega(&string, 4) / (4.0 * omega(&string, 1)) - 1.0 > 0.01);
        }
        let mut texture = TextureDrive::try_new(None).unwrap();
        for step in 1..=32 {
            let (pressure, physical) = step_kc_member(
                &mut member, None, &mut texture, &[], &string, &gas, dt,
            ).unwrap();
            assert!(pressure.is_finite());
            for (mode, &amplitude) in amplitudes.iter().enumerate() {
                let w = omega(&string, mode + 1);
                let phase = w * f64::from(step) * dt;
                close(member.x[2 * mode], amplitude * phase.cos(), amplitude);
                close(member.x[2 * mode + 1], -w * amplitude * phase.sin(), w * amplitude);
                close(physical[mode], member.x[2 * mode] / mass_scale(&string), 1e-9);
                if step == 1 {
                    assert!(fs_phs::linear_exact::admit_midpoint_band(w, dt).is_err());
                    let observed_phase = (-member.x[2 * mode + 1] / w)
                        .atan2(member.x[2 * mode]);
                    let cents = 1200.0 * (observed_phase / (w * dt)).log2();
                    assert!(cents.abs() < 0.01, "{frequency} Hz mode {mode}: {cents} cents");
                }
            }
        }
    }
}

fn elliptic_k(parameter: f64) -> f64 {
    let (mut a, mut b) = (1.0, (1.0 - parameter).sqrt());
    for _ in 0..32 {
        let arithmetic = 0.5 * (a + b);
        let geometric = (a * b).sqrt();
        a = arithmetic;
        b = geometric;
        if (a - b).abs() <= 2.0 * f64::EPSILON * a {
            break;
        }
    }
    core::f64::consts::PI / (2.0 * a)
}

fn descending_zero_fraction(q0: f64, p0: f64, q1: f64, p1: f64, dt: f64) -> f64 {
    // Fourth-order Hermite event interpolation makes the measured period
    // sensitive to the second-order stepper, not a linear crossing estimate.
    let (mut low, mut high) = (0.0, 1.0);
    for _ in 0..54 {
        let u = 0.5 * (low + high);
        let value = (2.0 * u * u * u - 3.0 * u * u + 1.0) * q0
            + (u * u * u - 2.0 * u * u + u) * dt * p0
            + (-2.0 * u * u * u + 3.0 * u * u) * q1
            + (u * u * u - u * u) * dt * p1;
        if value > 0.0 { low = u; } else { high = u; }
    }
    0.5 * (low + high)
}

#[test]
fn finite_amplitude_kc_period_converges_to_the_duffing_elliptic_oracle() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut string = string(1_000.0, 1);
    let amplitude: f64 = 0.002;
    let w = omega(&string, 1);
    let gamma = 0.5 * w * w / (amplitude * amplitude);
    // Physical single-mode equation: Q'' + w^2 Q + gamma Q^3 = 0.
    string.axial_stiffness_n = gamma * 4.0 * string.lin_density_kg_m
        * string.length_m.powi(4) / core::f64::consts::PI.powi(4);
    let capital_omega = (w * w + gamma * amplitude * amplitude).sqrt();
    let parameter = gamma * amplitude * amplitude / (2.0 * capital_omega * capital_omega);
    let reference_period = 4.0 * elliptic_k(parameter) / capital_omega;
    let mut errors = Vec::new();
    for rate in [48_000.0, 96_000.0, 192_000.0, 384_000.0] {
        let dt = 1.0 / rate;
        let mut member = kc_string_member(
            &string, None, None, &[], &gas, 1.0, 1.0, dt,
        ).unwrap();
        member.x[0] = amplitude * mass_scale(&string);
        let energy = physical_string_energy(&string, &member.x);
        let mut texture = TextureDrive::try_new(None).unwrap();
        let mut period = None;
        for step in 0..1000 {
            let (q0, p0) = (member.x[0], member.x[1]);
            step_kc_member(&mut member, None, &mut texture, &[], &string, &gas, dt).unwrap();
            close(physical_string_energy(&string, &member.x), energy, energy);
            if member.x[0] <= 0.0 {
                let fraction = descending_zero_fraction(q0, p0, member.x[0], member.x[1], dt);
                period = Some(4.0 * (f64::from(step) + fraction) * dt);
                break;
            }
        }
        let period = period.expect("the finite-amplitude oscillator must cross zero");
        errors.push((period / reference_period - 1.0).abs());
        if rate == 384_000.0 {
            let cents = 1200.0 * (reference_period / period).log2();
            assert!(cents.abs() < 0.1, "finest nonlinear period error {cents} cents");
        }
    }
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.35 * pair[0], "period refinement errors {errors:?}");
    }
}

#[test]
fn kc_exact_contact_preserves_the_original_potential_and_held_force_work() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut string = string(800.0, 1);
    string.rayleigh = Some(RayleighParams { alpha_per_s: 4.0, beta_s: 0.0 });
    let obstacle = UnilateralObstacle {
        stations: vec![0.5],
        gaps_m: vec![5e-6],
        measures_m: vec![0.02],
        stiffness: 5e5,
        alpha: 1.0,
        mu_kinetic: 0.0,
        mu_static: 0.0,
        friction_regularization_m_s: None,
        internal_loss: 0.0,
        provenance: "manufactured unilateral elastic contact".into(),
    };
    let dt = 1.0 / 96_000.0;
    let mut member = kc_string_member(
        &string, None, None, std::slice::from_ref(&obstacle), &gas, 1.0, 1.0, dt,
    ).unwrap();
    assert!(matches!(&member.flow, KcFlow::Exact(_)));
    member.x[0] = mass_scale(&string) * 20e-6;
    let physical_energy = |x: &[f64]| {
        let penetration = (x[0] / mass_scale(&string) - obstacle.gaps_m[0]).max(0.0);
        physical_string_energy(&string, x)
            + 0.5 * obstacle.stiffness * obstacle.measures_m[0] * penetration * penetration
    };
    let initial = physical_energy(&member.x);
    close(member.flow.system().hamiltonian(&member.x), initial, initial);
    assert!(initial > 1.1 * physical_string_energy(&string, &member.x));
    let (mut heat, mut work) = (0.0, 0.0);
    let (mut contacted, mut separated) = (false, false);
    let force = [0.05];
    for _ in 0..160 {
        let record = member.flow.step_into(
            &member.x, &force, dt, &mut member.next, &mut member.average_output,
        ).unwrap();
        close(record.delta_h, physical_energy(&member.next) - physical_energy(&member.x), initial);
        close(record.supplied, force[0] * (member.next[0] - member.x[0]), initial);
        close(record.supplied, dt * force[0] * member.average_output[0], initial);
        assert!(record.dissipated >= 0.0);
        heat += record.dissipated;
        work += record.supplied;
        core::mem::swap(&mut member.x, &mut member.next);
        contacted |= member.x[0] / mass_scale(&string) > obstacle.gaps_m[0];
        separated |= member.x[0] / mass_scale(&string) < obstacle.gaps_m[0];
    }
    assert!(contacted && separated);
    assert!(heat > 0.0);
    close(physical_energy(&member.x) - initial + heat, work, initial);
}

#[test]
fn kc_exact_flow_keeps_every_bending_memory_coordinate_and_its_energy() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut string = string(300.0, 2);
    string.rayleigh = None;
    string.bending_stiffness_n_m2 = 0.7;
    string.relaxation_bending = Some(fs_scenario::PronyBending {
        branches: vec![
            fs_scenario::BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 0.5, relaxation_time_s: 0.0005,
            },
            fs_scenario::BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 0.3, relaxation_time_s: 0.005,
            },
        ],
        omega_band_rad_s: (0.0, 20_000.0),
        material_state_identity: None,
        source_properties: None,
    });
    validate_relaxation_modes(&string, 48_000).unwrap();
    let dt = 1.0 / 48_000.0;
    let mut member = kc_string_member(
        &string, Some(Pluck { station_frac: 0.37, height_m: 1e-4 }),
        None, &[], &gas, 1.0, 1.0, dt,
    ).unwrap();
    assert!(matches!(&member.flow, KcFlow::Exact(_)));
    assert_eq!(member.x.len(), 8);
    let initial_memory = member.x[4..].to_vec();
    assert!(initial_memory.iter().all(|z| *z != 0.0));
    let physical_energy = |x: &[f64]| {
        let mut energy = physical_string_energy(&string, x);
        for (branch, stiffness) in [0.5, 0.3].into_iter().enumerate() {
            for mode in 0..2 {
                let k = (mode + 1) as f64 * core::f64::consts::PI / string.length_m;
                let root = (stiffness * k.powi(4) / string.lin_density_kg_m).sqrt();
                energy += 0.5 * (root * x[2 * mode] - x[(2 + branch) * 2 + mode]).powi(2);
            }
        }
        energy
    };
    let initial = physical_energy(&member.x);
    close(member.flow.system().hamiltonian(&member.x), initial, initial);
    let force = [0.002, -0.001];
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..32 {
        let record = member.flow.step_into(
            &member.x, &force, dt, &mut member.next, &mut member.average_output,
        ).unwrap();
        let physical_work: f64 = force.iter().enumerate()
            .map(|(mode, f)| f * (member.next[2 * mode] - member.x[2 * mode])).sum();
        close(record.supplied, physical_work, initial);
        heat += record.dissipated;
        work += record.supplied;
        core::mem::swap(&mut member.x, &mut member.next);
        close(member.flow.system().hamiltonian(&member.x), physical_energy(&member.x), initial);
    }
    assert!(heat > 0.0);
    assert!(member.x[4..].iter().zip(initial_memory).any(|(a, b)| *a != b));
    close(physical_energy(&member.x) - initial + heat, work, initial);
}

#[test]
fn exact_kc_refuses_retained_modes_outside_the_sampled_observer_band() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let dt = 1.0 / 48_000.0;
    let upper_hz = 24_000.0 * ModalAcousticTimeBudget::audible_reference().nyquist_guard_fraction;
    let above = string(upper_hz * 1.001, 1);
    assert!(matches!(
        kc_string_member(&above, None, None, &[], &gas, 1.0, 1.0, dt),
        Err(AcousticRealizeError::InvalidDescription { .. })
    ));
    let mut bending = string(1000.0, 1);
    let k = core::f64::consts::PI / bending.length_m;
    bending.bending_stiffness_n_m2 = bending.lin_density_kg_m
        * (core::f64::consts::TAU * upper_hz * 1.1).powi(2) / k.powi(4);
    assert!(kc_string_member(&bending, None, None, &[], &gas, 1.0, 1.0, dt).is_err());
}

#[test]
fn rejected_kc_trials_keep_mechanics_observer_and_texture_for_retry() {
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut string = string(1000.0, 1);
    string.axial_stiffness_n = 5e6;
    let dt = 1.0 / 48_000.0;
    let mut member = kc_string_member(
        &string, None, None, &[], &gas, 1.0, 1.0, dt,
    ).unwrap();
    member.x[0] = 0.002 * mass_scale(&string);
    let before = member.x.clone();
    let bow = BowStroke {
        station_frac: 0.37, normal_force_n: 0.5, velocity_m_s: 0.2,
        mu_static: 0.8, mu_dynamic: 0.4, stribeck_m_s: 0.01,
    };
    let spec = ContactTexture {
        rms_height_m: 1e-6, hurst_exponent: 0.7, min_cycles: 1, max_cycles: 4,
        phase_seed: 7, track_length_m: 0.01, tangent_stiffness_n_m: 1000.0,
    };
    let mut texture = TextureDrive::try_new(Some(spec)).unwrap();
    assert!(step_kc_member(
        &mut member, Some(bow), &mut texture, &[], &string, &gas, 2.0 * dt,
    ).is_err());
    assert_eq!(member.x, before);
    assert_eq!(texture.inner.as_ref().unwrap().path_m, 0.0);

    let KcFlow::Exact(flow) = &mut member.flow else { panic!("small KC model must use exact flow"); };
    flow.set_max_iterations(0).unwrap();
    assert!(matches!(
        step_kc_member(&mut member, Some(bow), &mut texture, &[], &string, &gas, dt),
        Err(AcousticRealizeError::Nonlinear(_))
    ));
    assert_eq!(member.x, before);
    assert!(member.observer.previous_acceleration.is_none());
    assert_eq!(texture.inner.as_ref().unwrap().path_m, 0.0);

    let KcFlow::Exact(flow) = &mut member.flow else { unreachable!() };
    flow.set_max_iterations(50).unwrap();
    let weight = member.observer.weights[0];
    member.observer.weights[0].0 = f64::MAX;
    assert!(matches!(
        step_kc_member(&mut member, Some(bow), &mut texture, &[], &string, &gas, dt),
        Err(AcousticRealizeError::InvalidDescription { .. })
    ));
    assert_eq!(member.x, before);
    assert!(member.observer.previous_acceleration.is_none());
    assert_eq!(texture.inner.as_ref().unwrap().path_m, 0.0);
    member.observer.weights[0] = weight;

    let mut fresh = kc_string_member(
        &string, None, None, &[], &gas, 1.0, 1.0, dt,
    ).unwrap();
    fresh.x.copy_from_slice(&before);
    let mut fresh_texture = TextureDrive::try_new(Some(spec)).unwrap();
    let expected = step_kc_member(
        &mut fresh, Some(bow), &mut fresh_texture, &[], &string, &gas, dt,
    ).unwrap();
    let actual = step_kc_member(
        &mut member, Some(bow), &mut texture, &[], &string, &gas, dt,
    ).unwrap();
    assert_eq!(member.x, fresh.x);
    assert_eq!(actual, expected);
    assert_eq!(member.observer.previous_acceleration, fresh.observer.previous_acceleration);
    assert_eq!(texture.inner.as_ref().unwrap().path_m, fresh_texture.inner.as_ref().unwrap().path_m);
    assert!(texture.inner.as_ref().unwrap().path_m > 0.0);
}
