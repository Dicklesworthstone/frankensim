use super::*;
use fs_phs::harmonic::resonance::{ResonanceControls, search_resonance};

struct Cavity {
    mass: Vec<f64>,
    damping: Vec<f64>,
    stiffness: Vec<f64>,
    coupling: Vec<f64>,
}
impl Cavity {
    // Independent consistent P1 unit-tetrahedron matrices, rho=kappa=0.2.
    // Affine dilation has B_i=1/8 and physical compression loading 0.3.
    fn new() -> Self {
        let mut mass = vec![0.0; 25];
        let mut damping = vec![0.0; 25];
        let mut stiffness = vec![0.0; 25];
        mass[0] = 1.0;
        stiffness[0] = 0.7;
        let gradients = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        for i in 0..4 {
            damping[i + 1] = -0.125;
            damping[(i + 1) * 5] = 0.125;
            for j in 0..4 {
                mass[(i + 1) * 5 + j + 1] = if i == j { 1.0 / 12.0 } else { 1.0 / 24.0 };
                stiffness[(i + 1) * 5 + j + 1] = (0..3)
                    .map(|k| gradients[i][k] * gradients[j][k])
                    .sum::<f64>()
                    / 1.2;
            }
        }
        Self {
            mass,
            damping,
            stiffness,
            coupling: vec![0.125; 4],
        }
    }
    fn source(&self) -> Source<'_> {
        Source {
            structural: 1,
            coordinates: 5,
            mass: &self.mass,
            damping: &self.damping,
            stiffness: &self.stiffness,
            coupling: &self.coupling,
        }
    }
}
fn branch() -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: vec![1.0, 0.0, 0.0, 0.0, 0.0],
        stiffness: 1.25,
        relaxation_time_s: 1.0,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-10 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn scalar(omega: f64, force: C64, flow: C64) -> (C64, C64, C64) {
    let s = C64::new(0.0, -omega);
    let stiffness = s * s + C64::ONE + s.scale(1.25) / (C64::ONE + s);
    let u = (force + flow.scale(0.6) / s) / stiffness;
    (
        u,
        u.scale(-0.6) + flow.scale(1.2) / s,
        u * s.scale(1.25) / (C64::ONE + s),
    )
}

#[test]
fn full_p1_memory_response_matches_original_cubic_and_physical_accounts() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let force = C64::new(0.7, -0.2);
    let flow = C64::new(-0.02, 0.05);
    for omega in [0.12, 0.5, 1.0, 1.4, 2.0, 8.0] {
        let result = prepared
            .factor(omega, &mut || false)
            .unwrap()
            .solve(&[force], &[flow.scale(0.25); 4], &mut || false)
            .unwrap();
        let (u, p, stress) = scalar(omega, force, flow);
        let s = C64::new(0.0, -omega);
        close(result.structural_displacement[0], u, 3e-10);
        close(result.structural_velocity[0], s * u, 3e-10);
        close(result.fluid_force[0], p.scale(0.5), 3e-10);
        for i in 0..4 {
            close(result.pressure_pa[i], p, 3e-10);
            close(result.pressure_potential_pa_s[i], p / s, 3e-10);
            close(result.displaced_flow_m3_s[i], (s * u).scale(0.125), 3e-10);
        }
        close(result.memory[0].stress, stress, 3e-10);
        let energy = result.spatial_energy;
        near(
            energy.structural_kinetic_j,
            0.25 * omega * omega * u.norm_sq(),
        );
        near(energy.structural_strain_j, 0.175 * u.norm_sq());
        near(energy.fluid_compression_j, 5.0 / 24.0 * p.norm_sq());
        near(energy.fluid_kinetic_j, 0.0);
        near(result.relaxation_energy_j, stress.norm_sq() / 5.0);
        near(
            result.average_relaxation_dissipation_w,
            stress.norm_sq() / 2.5,
        );
        near(result.average_wall_dissipation_w, 0.0);
        near(result.average_structural_dissipation_w, 0.0);
        for defect in [
            result.power_balance_residual_w,
            result.structural_power_balance_residual_w(),
            result.acoustic_power_balance_residual_w(),
            result.relaxation_power_balance_residual_w(),
            result.memory[0].power_balance_residual_w(),
        ] {
            near(defect, 0.0);
        }
        assert!(result.relative_residual < 1e-10);
    }
}

#[test]
fn frequency_tangent_retains_both_source_slopes_and_pressure_kinematics() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let (omega, f, q, df, dq) = (
        1.37,
        C64::new(0.7, -0.3),
        C64::new(0.02, 0.01),
        C64::new(-0.08, 0.06),
        C64::new(-0.01, 0.04),
    );
    let factor = prepared.factor(omega, &mut || false).unwrap();
    let result = factor
        .solve_with_frequency_derivative(
            &[f],
            &[q.scale(0.25); 4],
            Some((&[df], &[dq.scale(0.25); 4])),
            &mut || false,
        )
        .unwrap();
    let s = C64::new(0.0, -omega);
    let ds = C64::new(0.0, -1.0);
    let one_s = C64::ONE + s;
    let t = s * s + C64::ONE + s.scale(1.25) / one_s;
    let dt = ds * (s.scale(2.0) + C64::from_re(1.25) / (one_s * one_s));
    let load = f + q.scale(0.6) / s;
    let dload = df + dq.scale(0.6) / s - q.scale(0.6) * ds / (s * s);
    let du = (dload * t - load * dt) / (t * t);
    let (u, p, _) = scalar(omega, f, q);
    let dp = du.scale(-0.6) + dq.scale(1.2) / s - q.scale(1.2) * ds / (s * s);
    let dstress = (ds * u / (one_s * one_s) + s * du / one_s).scale(1.25);
    close(
        result.structural_displacement_derivative_per_rad_s[0],
        du,
        5e-10,
    );
    close(
        result.structural_velocity_derivative_per_rad_s[0],
        ds * u + s * du,
        5e-10,
    );
    close(
        result.relaxation_stress_derivative_per_rad_s[0],
        dstress,
        5e-10,
    );
    close(
        result.fluid_force_derivative_per_rad_s[0],
        dp.scale(0.5),
        5e-10,
    );
    for i in 0..4 {
        close(result.pressure_derivative_pa_per_rad_s[i], dp, 5e-10);
        close(
            result.pressure_potential_derivative_pa_s_per_rad_s[i],
            dp / s - p * ds / (s * s),
            5e-10,
        );
        close(
            result.displaced_flow_derivative_m3_s_per_rad_s[i],
            (ds * u + s * du).scale(0.125),
            5e-10,
        );
    }
    let h = 1e-5;
    let plus = scalar(omega + h, f + df.scale(h), q + dq.scale(h));
    let minus = scalar(omega - h, f - df.scale(h), q - dq.scale(h));
    close((plus.0 - minus.0).scale(0.5 / h), du, 2e-8);
    close((plus.1 - minus.1).scale(0.5 / h), dp, 2e-8);
    assert!(result.relative_derivative_residual < 1e-10);
}

#[test]
fn potential_memory_is_causal_wall_admittance_with_separate_acoustic_heat() {
    let cavity = Cavity::new();
    let branches = [MaxwellDampingBranch {
        strain: vec![0.0, 0.25, 0.25, 0.25, 0.25],
        stiffness: 0.8,
        relaxation_time_s: 0.4,
    }];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let (f, q) = (C64::new(0.6, 0.1), C64::new(0.07, -0.03));
    for omega in [0.2, 0.9, 1.5, 7.0] {
        let s = C64::new(0.0, -omega);
        let a = s.scale(5.0 / 6.0) + C64::from_re(0.8) / (s + C64::from_re(2.5));
        let u = (f * a + q.scale(0.5)) / ((s * s + C64::from_re(0.7)) * a + s.scale(0.25));
        let p = (q - (s * u).scale(0.5)) / a;
        let stress = p.scale(0.8) / (s + C64::from_re(2.5));
        let result = prepared
            .factor(omega, &mut || false)
            .unwrap()
            .solve(&[f], &[q.scale(0.25); 4], &mut || false)
            .unwrap();
        close(result.structural_displacement[0], u, 5e-10);
        for &value in &result.pressure_pa {
            close(value, p, 5e-10);
        }
        close(result.memory[0].stress, stress, 5e-10);
        near(result.average_structural_to_relaxation_power_w, 0.0);
        near(
            result.average_acoustic_to_relaxation_power_w,
            result.average_relaxation_dissipation_w,
        );
        near(result.average_wall_dissipation_w, 0.0);
        near(result.power_balance_residual_w, 0.0);
        near(result.structural_power_balance_residual_w(), 0.0);
        near(result.acoustic_power_balance_residual_w(), 0.0);
    }
}

#[test]
fn driven_peak_and_actual_bandwidth_use_complete_memory_response() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let result = search_resonance(
        (0.8, 2.4),
        ResonanceControls {
            absolute_frequency_tolerance_rad_s: 1e-10,
            relative_frequency_tolerance: 1e-10,
            maximum_evaluations: 256,
            measure_half_power_bandwidth: true,
        },
        || false,
        |omega, _, cancelled| {
            let response = prepared
                .factor(omega, cancelled)?
                .solve_with_frequency_derivative(&[C64::ONE], &[C64::ZERO; 4], None, cancelled)?;
            Ok::<_, Error>((
                response.response.structural_displacement[0],
                response.structural_displacement_derivative_per_rad_s[0],
                response
                    .relative_derivative_residual
                    .max(response.response.relative_residual),
            ))
        },
    )
    .unwrap();
    // Independent stationary polynomial for |s²+1+1.25*s/(1+s)|², x=omega².
    let polynomial = |x: f64| 2.0 * x.powi(4) + 1.5 * x.powi(3) - 7.5 * x * x - 4.9375 * x + 2.0625;
    let (mut low, mut high) = (1.0, 4.0);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if polynomial(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    near(
        result.peak.angular_frequency_rad_s,
        (0.5 * (low + high)).sqrt(),
    );
    assert!((result.peak.angular_frequency_rad_s - 31.0_f64.sqrt() / 4.0).abs() > 0.01);
    let half = result.half_power.unwrap();
    let amplitude_squared = result.observation.norm_sq();
    for frequency in [
        half.lower.angular_frequency_rad_s,
        half.upper.angular_frequency_rad_s,
    ] {
        near(
            scalar(frequency, C64::ONE, C64::ZERO).0.norm_sq(),
            0.5 * amplitude_squared,
        );
    }
    near(
        half.response_quality_factor,
        result.peak.angular_frequency_rad_s / half.bandwidth_rad_s,
    );
}

#[test]
fn manufactured_nonuniform_fields_retain_mutual_memory_and_gyroscopic_rows() {
    let base = Cavity::new();
    let n = 6;
    let mut mass = vec![0.0; n * n];
    let mut damping = vec![0.0; n * n];
    let mut stiffness = vec![0.0; n * n];
    let mut coupling = vec![0.0; 8];
    for (target, values) in [
        (&mut mass, [1.2, 0.1, 0.1, 0.9]),
        (&mut damping, [0.4, 0.7, -0.5, 0.3]),
        (&mut stiffness, [5.0, 0.4, 0.4, 3.0]),
    ] {
        for i in 0..2 {
            for j in 0..2 {
                target[i * n + j] = values[i * 2 + j];
            }
        }
    }
    for i in 0..4 {
        for j in 0..4 {
            mass[(i + 2) * n + j + 2] = base.mass[(i + 1) * 5 + j + 1];
            stiffness[(i + 2) * n + j + 2] = base.stiffness[(i + 1) * 5 + j + 1];
        }
        damping[(i + 2) * n + i + 2] = 0.03;
        for j in 0..2 {
            let b = if j == 0 {
                0.125
            } else {
                0.03 * (i as f64 - 1.5)
            };
            coupling[i * 2 + j] = b;
            damping[(i + 2) * n + j] = b;
            damping[j * n + i + 2] = -b;
        }
    }
    let branches = [
        MaxwellDampingBranch {
            strain: vec![1.0, 0.4, 0.2, -0.1, 0.3, -0.2],
            stiffness: 0.8,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![-0.3, 0.6, -0.2, 0.4, 0.1, 0.3],
            stiffness: 0.5,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![0.0; n],
            stiffness: 2.0,
            relaxation_time_s: 0.1,
        },
    ];
    let source = Source {
        structural: 2,
        coordinates: n,
        mass: &mass,
        damping: &damping,
        stiffness: &stiffness,
        coupling: &coupling,
    };
    let prepared = Prepared::new(source, &branches, Limits::default(), &mut || false).unwrap();
    let omega = 2.3;
    let s = C64::new(0.0, -omega);
    let expected: Vec<_> = (0..n)
        .map(|i| C64::new(0.2 - 0.05 * i as f64, 0.01 * (i * i) as f64 - 0.08))
        .collect();
    let stress: Vec<_> = branches
        .iter()
        .map(|branch| {
            let strain = branch
                .strain
                .iter()
                .zip(&expected)
                .fold(C64::ZERO, |sum, (&b, &q)| sum + q.scale(b));
            strain * s.scale(branch.stiffness)
                / (s + C64::from_re(branch.relaxation_time_s.recip()))
        })
        .collect();
    let input: Vec<_> = (0..n)
        .map(|i| {
            let spatial = (0..n).fold(C64::ZERO, |sum, j| {
                sum + expected[j]
                    * (C64::from_re(stiffness[i * n + j])
                        + s.scale(damping[i * n + j])
                        + s * s.scale(mass[i * n + j]))
            });
            branches
                .iter()
                .zip(&stress)
                .fold(spatial, |sum, (b, &stress)| sum + stress.scale(b.strain[i]))
        })
        .collect();
    let result = prepared
        .factor(omega, &mut || false)
        .unwrap()
        .solve(&input[..2], &input[2..], &mut || false)
        .unwrap();
    for (actual, &expected) in result
        .structural_displacement
        .iter()
        .chain(&result.pressure_potential_pa_s)
        .zip(&expected)
    {
        close(*actual, expected, 3e-10);
    }
    for (actual, &expected) in result
        .structural_velocity
        .iter()
        .chain(&result.pressure_pa)
        .zip(&expected)
    {
        close(*actual, s * expected, 3e-10);
    }
    for (actual, &expected) in result.memory.iter().zip(&stress) {
        close(actual.stress, expected, 3e-10);
        near(actual.power_balance_residual_w(), 0.0);
    }
    let velocity = [s * expected[0], s * expected[1]];
    let structural_heat = 0.5
        * (0.4 * velocity[0].norm_sq()
            + 0.3 * velocity[1].norm_sq()
            + 0.2 * (velocity[0].conj() * velocity[1]).re);
    near(result.average_structural_dissipation_w, structural_heat);
    for defect in [
        result.power_balance_residual_w,
        result.structural_power_balance_residual_w(),
        result.acoustic_power_balance_residual_w(),
        result.relaxation_power_balance_residual_w(),
    ] {
        near(defect, 0.0);
    }
    assert_eq!(result.memory[2].stress, C64::ZERO);
}

#[test]
fn zero_memory_limit_refusals_cancellation_and_factor_reuse() {
    let cavity = Cavity::new();
    let prepared = Prepared::new(cavity.source(), &[], Limits::default(), &mut || false).unwrap();
    let factor = prepared.factor(0.6, &mut || false).unwrap();
    let response = factor
        .solve(&[C64::ONE], &[C64::ZERO; 4], &mut || false)
        .unwrap();
    close(
        response.structural_displacement[0],
        C64::from_re(1.0 / (1.0 - 0.36)),
        1e-10,
    );
    assert!(response.memory.is_empty());
    near(response.relaxation_energy_j, 0.0);
    assert!(matches!(
        factor.solve(&[C64::ONE], &[C64::ZERO; 4], &mut || true),
        Err(Error::Cancelled)
    ));
    assert!(factor.solve(&[], &[C64::ZERO; 4], &mut || false).is_err());
    let mut calls = 0;
    assert!(
        factor
            .solve(&[C64::ONE], &[C64::ZERO; 4], &mut || {
                calls += 1;
                calls > 20
            })
            .is_err()
    );
    let replay = factor
        .solve(&[C64::ONE], &[C64::ZERO; 4], &mut || false)
        .unwrap();
    assert_eq!(
        response.structural_displacement,
        replay.structural_displacement
    );
    let limits = Limits {
        maximum_branches: 0,
        ..Limits::default()
    };
    assert!(matches!(
        Prepared::new(cavity.source(), &[branch()], limits, &mut || false),
        Err(Error::Limit { .. })
    ));
    let limits = Limits {
        maximum_adapter_bytes: 1,
        ..Limits::default()
    };
    assert!(matches!(
        Prepared::new(cavity.source(), &[], limits, &mut || false),
        Err(Error::Limit { .. })
    ));
    let singular = Source {
        structural: 1,
        coordinates: 2,
        mass: &[1.0, 0.0, 0.0, 1.0],
        damping: &[0.0; 4],
        stiffness: &[1.0, 0.0, 0.0, 0.0],
        coupling: &[0.0],
    };
    let singular = Prepared::new(singular, &[], Limits::default(), &mut || false).unwrap();
    assert!(matches!(
        singular.factor(1.0, &mut || false),
        Err(Error::Harmonic(HarmonicError::Singular { .. }))
    ));
}

#[path = "tests/curvature.rs"]
mod curvature;

#[path = "tests/parameter.rs"]
mod parameter;
#[path = "tests/parameter_coupled.rs"]
mod parameter_coupled;

#[path = "tests/power.rs"]
mod power;

#[path = "tests/transfer.rs"]
mod transfer;
