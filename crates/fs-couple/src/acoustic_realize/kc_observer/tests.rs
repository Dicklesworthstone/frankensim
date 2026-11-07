//! Independent physical force derivatives and actual renderer integration.
use super::super::{
    AmbientGas, GasState, PrestressedString, RayleighParams, TextureDrive, UnilateralObstacle,
    gas_state, kc_string_member, step_kc_member,
};
use super::*;

fn string() -> PrestressedString {
    PrestressedString {
        length_m: 0.1,
        tension_n: 100.0,
        lin_density_kg_m: 0.004,
        axial_stiffness_n: 20_000.0,
        width_m: 0.001,
        n_modes: 2,
        damping_ratio: 0.0,
        rayleigh: Some(RayleighParams {
            alpha_per_s: 12.0,
            beta_s: 2e-7,
        }),
        bending_stiffness_n_m2: 0.01,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: false,
        tension_loss_factor: None,
        bending_loss_factor: None,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 8e-12 * scale.max(f64::MIN_POSITIVE),
        "{actual:.17e} != {expected:.17e}, scale={scale:.17e}"
    );
}

fn mass(s: &PrestressedString) -> f64 {
    s.lin_density_kg_m * s.length_m / 2.0
}

fn waves(s: &PrestressedString) -> [f64; 2] {
    let k = core::f64::consts::PI / s.length_m;
    [k, 2.0 * k]
}

fn omega_sq(s: &PrestressedString, k: f64) -> f64 {
    (s.tension_n * k.powi(2) + s.bending_stiffness_n_m2 * k.powi(4)) / s.lin_density_kg_m
}

fn pressure(s: &PrestressedString, gas: &GasState, a: [f64; 2], j: [f64; 2]) -> f64 {
    let pi = core::f64::consts::PI;
    let scale = gas.density * s.width_m / (4.0 * pi * mass(s).sqrt());
    scale
        * (2.0 * s.length_m / pi * a[0] - s.length_m.powi(2) / (2.0 * pi * gas.sound_speed) * j[1])
}

fn obstacle(alpha: f64) -> UnilateralObstacle {
    UnilateralObstacle {
        stations: vec![0.31, 0.73],
        gaps_m: vec![1e-4, 0.003],
        measures_m: vec![0.017, 0.026],
        stiffness: 7e6,
        alpha,
        mu_kinetic: 0.0,
        mu_static: 0.0,
        friction_regularization_m_s: None,
        internal_loss: 0.0,
        provenance: "manufactured distributed elastic contact".into(),
    }
}

#[test]
fn shared_stretching_bending_damping_and_contact_reach_endpoint_jerk() {
    let s = string();
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let wave = waves(&s);
    let q = [mass(&s).sqrt() * 3e-4, mass(&s).sqrt() * 2e-4];
    let v = [mass(&s).sqrt() * -0.05, mass(&s).sqrt() * 0.08];
    let squared_slope: f64 = (0..2)
        .map(|i| wave[i].powi(2) * q[i].powi(2) / mass(&s))
        .sum();
    let slope_rate: f64 = (0..2)
        .map(|i| 2.0 * wave[i].powi(2) * q[i] * v[i] / mass(&s))
        .sum();
    let gamma = s.axial_stiffness_n / (4.0 * s.lin_density_kg_m);
    // Physical T+EI dispersion and d=alpha+beta*omega^2, independent of the
    // pHS matrices and storage's Hessian implementation.
    let damping = wave.map(|k| 12.0 + 2e-7 * omega_sq(&s, k));
    for alpha in [1.0, 1.5, 2.0] {
        let ob = obstacle(alpha);
        let mut member = kc_string_member(
            &s,
            None,
            None,
            core::slice::from_ref(&ob),
            &gas,
            1.0,
            1.0,
            1.0 / 48_000.0,
        )
        .unwrap();
        member.x.copy_from_slice(&[q[0], v[0], q[1], v[1]]);
        // Rebinding a prescribed held load changes its smooth endpoint force,
        // not an artificial derivative spanning the two different inputs.
        for force in [[0.04, -0.03], [-0.09, 0.06]] {
            let mut a = [0.0; 2];
            let mut j = [0.0; 2];
            for i in 0..2 {
                a[i] = force[i]
                    - omega_sq(&s, wave[i]) * q[i]
                    - damping[i] * v[i]
                    - gamma * wave[i].powi(2) * squared_slope * q[i];
                j[i] = -omega_sq(&s, wave[i]) * v[i]
                    - gamma * wave[i].powi(2) * (squared_slope * v[i] + slope_rate * q[i]);
            }
            let mut active_points = 0;
            for point in 0..ob.stations.len() {
                let row = [1.0, 2.0].map(|mode| {
                    (mode * core::f64::consts::PI * ob.stations[point]).sin() / mass(&s).sqrt()
                });
                let penetration = row[0] * q[0] + row[1] * q[1] - ob.gaps_m[point];
                if penetration <= 0.0 {
                    continue;
                }
                active_points += 1;
                let rate = row[0] * v[0] + row[1] * v[1];
                let force = ob.measures_m[point] * ob.stiffness * penetration.powf(alpha);
                let derivative = ob.measures_m[point]
                    * ob.stiffness
                    * alpha
                    * penetration.powf(alpha - 1.0)
                    * rate;
                for i in 0..2 {
                    a[i] -= row[i] * force;
                    j[i] -= row[i] * derivative;
                }
            }
            assert_eq!(active_points, 1, "the second point must stay inactive");
            for i in 0..2 {
                j[i] -= damping[i] * a[i];
            }
            let actual = member
                .observation
                .observe(&member.observer, member.flow.system(), &member.x, &force)
                .unwrap();
            for i in 0..2 {
                close(member.observation.state_rate[2 * i + 1], a[i], a[i].abs());
                close(member.observation.state_second[2 * i + 1], j[i], j[i].abs());
            }
            let expected = pressure(&s, &gas, a, j);
            close(actual, expected, expected.abs());
            assert!(member.observer.previous_acceleration.is_none());
        }
    }
}

#[test]
fn nonlinear_renderer_keeps_upper_even_partial_amplitude_and_phase() {
    let mut s = string();
    s.length_m = 0.002;
    s.lin_density_kg_m = 1e-6;
    s.width_m = 1e-5;
    s.tension_n = s.lin_density_kg_m * (3200.0 * s.length_m).powi(2);
    s.axial_stiffness_n = 0.0;
    s.bending_stiffness_n_m2 = 0.0;
    s.rayleigh = Some(RayleighParams {
        alpha_per_s: 0.0,
        beta_s: 0.0,
    });
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let dt = 1.0 / 8000.0;
    let mut member = kc_string_member(&s, None, None, &[], &gas, 1.0, 1.0, dt).unwrap();
    let amplitude = 1e-10;
    member.x[2] = amplitude;
    let omega = 2.0 * core::f64::consts::PI * 3200.0;
    let peak = pressure(&s, &gas, [0.0; 2], [0.0, amplitude * omega.powi(3)]);
    let mut texture = TextureDrive::try_new(None).unwrap();
    let (mut in_phase, mut quadrature, mut error) = (0.0, 0.0, 0.0_f64);
    // Forty complete cycles, f_2 at 80% of Nyquist, k_air L < 0.12.
    // The old acceleration difference scales the amplitude by sinc(0.4*pi)
    // and lags by 0.4*pi even when the mechanical propagation is exact.
    for sample in 1..=100 {
        let phase = omega * f64::from(sample) * dt;
        let (actual, _) =
            step_kc_member(&mut member, None, &mut texture, &[], &s, &gas, dt).unwrap();
        let normalized = actual / peak;
        in_phase += normalized * phase.sin() / 50.0;
        quadrature += normalized * phase.cos() / 50.0;
        error = error.max((normalized - phase.sin()).abs());
    }
    assert!((in_phase - 1.0).abs() < 2e-10, "gain {in_phase}");
    assert!(quadrature.abs() < 2e-10, "quadrature {quadrature}");
    assert!(error < 2e-10, "endpoint waveform error {error}");
    assert!(member.observer.previous_acceleration.is_none());
}

#[test]
fn bending_memory_radiates_while_modal_velocity_and_acceleration_are_zero() {
    use fs_scenario::acoustic::{BendingRelaxationBranch, PronyBending};
    let mut s = string();
    s.rayleigh = None;
    s.relaxation_bending = Some(PronyBending {
        branches: vec![
            BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 0.006,
                relaxation_time_s: 0.0007,
            },
            BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 0.003,
                relaxation_time_s: 0.004,
            },
        ],
        omega_band_rad_s: (0.0, 30_000.0),
        material_state_identity: None,
        source_properties: None,
    });
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut member = kc_string_member(&s, None, None, &[], &gas, 1.0, 1.0, 1.0 / 48_000.0).unwrap();
    let wave = waves(&s);
    let q = [mass(&s).sqrt() * 1e-4, mass(&s).sqrt() * 2e-4];
    member.x[0] = q[0];
    member.x[2] = q[1];
    let squared_slope: f64 = (0..2)
        .map(|i| wave[i].powi(2) * q[i].powi(2) / mass(&s))
        .sum();
    let gamma = s.axial_stiffness_n / (4.0 * s.lin_density_kg_m);
    let mut force = [0.0; 2];
    let mut jerk = [0.0; 2];
    for i in 0..2 {
        force[i] = omega_sq(&s, wave[i]) * q[i] + gamma * wave[i].powi(2) * squared_slope * q[i];
        for (arm, branch) in s
            .relaxation_bending
            .as_ref()
            .unwrap()
            .branches
            .iter()
            .enumerate()
        {
            let root =
                (branch.relaxing_stiffness_n_m2 * wave[i].powi(4) / s.lin_density_kg_m).sqrt();
            let memory = 4 + 2 * arm + i;
            member.x[memory] = 0.2 * arm as f64 * root * q[i];
            let strain = root * q[i] - member.x[memory];
            let memory_rate = strain / branch.relaxation_time_s;
            force[i] += root * strain;
            jerk[i] += root * memory_rate;
        }
    }
    let actual = member
        .observation
        .observe(&member.observer, member.flow.system(), &member.x, &force)
        .unwrap();
    for i in 0..2 {
        close(
            member.observation.state_rate[2 * i + 1],
            0.0,
            force[i].abs(),
        );
        close(
            member.observation.state_second[2 * i + 1],
            jerk[i],
            jerk[i].abs(),
        );
    }
    // The even mode emits a finite signed first moment solely because the
    // original Maxwell memories are relaxing. -omega^2*v-d*a would miss it.
    let expected = pressure(&s, &gas, [0.0; 2], jerk);
    assert!(expected.abs() > 1e-6);
    close(actual, expected, expected.abs());
}

#[test]
fn a_linear_contact_kink_refuses_jerk_but_quiet_contact_is_observable() {
    let s = string();
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let mut ob = obstacle(1.0);
    ob.gaps_m[0] = 0.0;
    let mut member =
        kc_string_member(&s, None, None, &[ob], &gas, 1.0, 1.0, 1.0 / 48_000.0).unwrap();
    let quiet = member
        .observation
        .observe(&member.observer, member.flow.system(), &member.x, &[0.0; 2])
        .unwrap();
    assert_eq!(quiet, 0.0);
    member.x[1] = 1e-3;
    assert!(matches!(
        member
            .observation
            .observe(&member.observer, member.flow.system(), &member.x, &[0.0; 2],),
        Err(AcousticRealizeError::InvalidDescription {
            what: "KC endpoint radiation requires a finite smooth storage tangent",
        })
    ));
    assert!(member.observer.previous_acceleration.is_none());
}
