//! G1/G3 physical parameter identities, original-memory recovery and fixed gauge.
use super::super::super::super::{CoupledPencil, MechanicalHarmonicLimits, Wall};
use super::super::super::MechanicalPoleOptions;
use super::super::MechanicalPoleResponseOptions;
use super::*;
use fs_sparse::{Coo, Csr};

fn scalar(value: f64) -> Csr {
    let mut a = Coo::new(1, 1);
    a.push(0, 0, value);
    a.assemble()
}
fn limits() -> MechanicalPoleLimits {
    MechanicalPoleLimits {
        max_evaluations: 128,
        max_payload_bytes: 2 << 20,
        physical_residual_tolerance: 2e-11,
    }
}
fn owner(reference: f64, k: f64, r: f64, sigma: f64, arms: Vec<(f64, f64)>) -> CoupledPencil {
    let fluid = [scalar(0.0), scalar(1.0), scalar(0.0), scalar(0.0)];
    CoupledPencil::prepare(
        reference,
        fluid.each_ref(),
        &[0],
        vec![Wall {
            dofs: vec![0],
            mass: scalar(1.0),
            surface_density: sigma,
            stiffness: k,
            resistance: r,
            arms,
        }],
        MechanicalHarmonicLimits::default(),
        &mut || false,
    )
    .unwrap()
}
fn pole(owner: &CoupledPencil, s: C64, seed: &[C64]) -> super::super::MechanicalResonancePole {
    owner
        .resonance_pole(
            s,
            seed,
            &MechanicalPoleOptions {
                spectral_scale: 3.0,
                step_tolerance: 2e-13,
                residual_tolerance: 2e-13,
                max_iterations: 50,
                max_dimension: 16,
                max_dense_bytes: 2 << 20,
                max_work: 1_000_000_000,
            },
            limits(),
            &mut || false,
        )
        .unwrap()
}
fn prepare<'a>(owner: &'a CoupledPencil, s: C64, seed: &[C64]) -> PreparedCoupledResonance<'a> {
    owner
        .prepare_resonance_response(
            &pole(owner, s, seed),
            MechanicalPoleResponseOptions {
                spectral_scale: 3.0,
                pole_residual_tolerance: 2e-12,
                pole_step_tolerance: 2e-12,
                linear_residual_tolerance: 2e-12,
                max_dimension: 16,
                max_vectors: 1,
                max_dense_bytes: 2 << 20,
                max_work: 1_000_000_000,
            },
            limits(),
            &mut || false,
        )
        .unwrap()
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 3e-9 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}
fn direction(
    relaxation: &[MechanicalRelaxationParameterDirection],
) -> MechanicalWallParameterDirection<'_> {
    MechanicalWallParameterDirection {
        wall: 0,
        surface_density_kg_m2: 0.0,
        stiffness_pa_m: 0.0,
        resistance_pa_s_m: 0.0,
        relaxation,
    }
}

#[test]
fn original_mass_stiffness_and_resistance_move_the_loaded_pole() {
    let s = C64::new(-0.1, (4.0_f64 - 0.01).sqrt());
    for reference in [0.7, 7.0] {
        let owner = owner(reference, 3.0, 0.2, 1.0, vec![]);
        let prepared = prepare(&owner, s, &[C64::from_re(-1.0), C64::ONE, s]);
        for (sigma, k, r) in [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.3, -0.2, 0.4),
        ] {
            let d = MechanicalWallParameterDirection {
                surface_density_kg_m2: sigma,
                stiffness_pa_m: k,
                resistance_pa_s_m: r,
                ..direction(&[])
            };
            let result = prepared
                .wall_parameter_derivative(&[d], limits(), &mut || false)
                .unwrap();
            let expected = -(s * s.scale(sigma) + s.scale(r) + C64::from_re(k))
                / (s.scale(2.0) + C64::from_re(0.2));
            close(result.rate_per_s_per_parameter, expected);
            close(
                result.coordinate_derivative[result.normalization_pivot],
                C64::ZERO,
            );
            close(
                result.normalization_value,
                prepared.mode().coordinates[result.normalization_pivot],
            );
            let w = prepared.mode().walls[0].displacement_m[0];
            let dw = result.walls[0].displacement_m_per_parameter[0];
            close(
                result.walls[0].normal_velocity_m_s_per_parameter[0],
                s * dw + expected * w,
            );
            close(result.pressure_pa_per_parameter[0], -dw);
            assert!(result.physical_derivative_residual < 2e-11);
        }
    }
}

#[test]
fn original_maxwell_stiffness_time_and_coordinate_chain_match_independent_law() {
    let rate = 3.0;
    let ka: f64 = 2.0;
    for frequency in [2.0, 20.0] {
        let s = C64::new(-0.5, frequency);
        let q = s + C64::from_re(rate);
        let dynamic = s * s + (s / q).scale(ka);
        let r = -dynamic.im / s.im;
        let k = -(dynamic + s.scale(r)).re - 1.0;
        let c = ka.sqrt();
        let e = (s / q).scale(c);
        let source = owner(0.7, k, r, 1.0, vec![(c, rate)]);
        let prepared = prepare(&source, s, &[C64::from_re(-1.0), C64::ONE, s, e]);
        for (dk, dtau) in [(1.0, 0.0), (0.0, 1.0), (-0.4, 0.3)] {
            let arms = [MechanicalRelaxationParameterDirection {
                stiffness_pa_m: dk,
                relaxation_time_s: dtau,
            }];
            let result = prepared
                .wall_parameter_derivative(&[direction(&arms)], limits(), &mut || false)
                .unwrap();
            // Independent scalar elimination is confined to this oracle:
            // L(s)=1+K+R*s+sigma*s^2+ka*s/(s+rate).
            let partial = (s / q).scale(dk) + (s / q / q).scale(ka * rate * rate * dtau);
            let spectral = s.scale(2.0) + C64::from_re(r) + C64::from_re(ka * rate) / q / q;
            let ds = -partial / spectral;
            close(result.rate_per_s_per_parameter, ds);
            let w = prepared.mode().walls[0].displacement_m[0];
            let dw = result.walls[0].displacement_m_per_parameter[0];
            let de = result.walls[0].relaxation_extension_per_parameter[0][0];
            let dz = result.walls[0].relaxation_memory_per_parameter[0][0];
            close(dz, dw.scale(c) + w.scale(0.5 * dk / c) - de);
            let expected_memory = (w / q).scale(c * rate);
            let expected_dmemory = (dw / q).scale(c * rate)
                + (w / q).scale(0.5 * dk * rate / c - c * rate * rate * dtau)
                - (expected_memory / q) * (ds - C64::from_re(rate * rate * dtau));
            close(dz, expected_dmemory);
        }
        // Rebuild the COMPLETE original pencil on both sides, then align only
        // the explicit physical pivot. No numerical source scaling is matched.
        let dk = 0.25;
        let dtau = 0.03;
        let arms = [MechanicalRelaxationParameterDirection {
            stiffness_pa_m: dk,
            relaxation_time_s: dtau,
        }];
        let analytic = prepared
            .wall_parameter_derivative(&[direction(&arms)], limits(), &mut || false)
            .unwrap();
        let h = 1e-5;
        let samples: [_; 2] = core::array::from_fn(|i| {
            let parameter = if i == 0 { -h } else { h };
            let other = owner(
                1.9,
                k,
                r,
                1.0,
                vec![(
                    (ka + parameter * dk).sqrt(),
                    1.0 / (1.0 / rate + parameter * dtau),
                )],
            );
            let mut mode = pole(&other, s, &prepared.mode().coordinates);
            let scale =
                analytic.normalization_value / mode.coordinates[analytic.normalization_pivot];
            for value in &mut mode.coordinates {
                *value = *value * scale;
            }
            mode
        });
        close(
            (samples[1].rate_per_s - samples[0].rate_per_s).scale(0.5 / h),
            analytic.rate_per_s_per_parameter,
        );
        for i in 0..analytic.coordinate_derivative.len() {
            close(
                (samples[1].coordinates[i] - samples[0].coordinates[i]).scale(0.5 / h),
                analytic.coordinate_derivative[i],
            );
        }
    }
}

#[test]
fn dark_equal_time_memory_mode_retains_its_material_derivative() {
    let source = owner(0.7, 3.0, 0.2, 1.0, vec![(2.0, 3.0), (2.0, 3.0)]);
    let s = C64::from_re(-3.0);
    let prepared = prepare(
        &source,
        s,
        &[
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ONE,
            C64::from_re(-1.0),
        ],
    );
    let time = [MechanicalRelaxationParameterDirection {
        stiffness_pa_m: 0.0,
        relaxation_time_s: 1.0,
    }; 2];
    let result = prepared
        .wall_parameter_derivative(&[direction(&time)], limits(), &mut || false)
        .unwrap();
    close(result.rate_per_s_per_parameter, C64::from_re(9.0));
    assert!(result.coordinate_derivative.iter().all(|x| x.abs() < 2e-11));
    let stiffness = [
        MechanicalRelaxationParameterDirection {
            stiffness_pa_m: 1.0,
            relaxation_time_s: 0.0,
        },
        MechanicalRelaxationParameterDirection::default(),
    ];
    let result = prepared
        .wall_parameter_derivative(&[direction(&stiffness)], limits(), &mut || false)
        .unwrap();
    close(result.rate_per_s_per_parameter, C64::ZERO);
    for arm in 0..2 {
        close(
            result.walls[0].relaxation_memory_per_parameter[arm][0],
            -result.walls[0].relaxation_extension_per_parameter[arm][0],
        );
    }
}

#[test]
fn invalid_directions_payload_and_cancellation_publish_no_tangent() {
    // k=2^-1074 and d(k)/dtheta=k must give d(sqrt(k))/dtheta=sqrt(k)/2,
    // even though k/2 itself is below the smallest representable number.
    let tiny = f64::from_bits(1);
    let coupling = tiny.sqrt();
    let (dc, _) = branch_direction(
        coupling,
        1.0,
        MechanicalRelaxationParameterDirection {
            stiffness_pa_m: tiny,
            relaxation_time_s: 0.0,
        },
    )
    .unwrap();
    assert_eq!(dc, coupling * 0.5);
    let s = C64::new(-0.1, (4.0_f64 - 0.01).sqrt());
    let source = owner(0.7, 3.0, 0.2, 1.0, vec![]);
    let prepared = prepare(&source, s, &[C64::from_re(-1.0), C64::ONE, s]);
    let valid = MechanicalWallParameterDirection {
        stiffness_pa_m: 1.0,
        ..direction(&[])
    };
    assert!(
        prepared
            .wall_parameter_derivative(&[valid, valid], limits(), &mut || false)
            .is_err()
    );
    for invalid in [
        MechanicalWallParameterDirection { wall: 1, ..valid },
        MechanicalWallParameterDirection {
            stiffness_pa_m: f64::NAN,
            ..valid
        },
    ] {
        assert!(
            prepared
                .wall_parameter_derivative(&[invalid], limits(), &mut || false)
                .is_err()
        );
    }
    let arm = [MechanicalRelaxationParameterDirection::default()];
    assert!(
        prepared
            .wall_parameter_derivative(&[direction(&arm)], limits(), &mut || false)
            .is_err()
    );
    let mut capped = limits();
    capped.max_payload_bytes = 1;
    assert!(matches!(
        prepared.wall_parameter_derivative(&[valid], capped, &mut || false),
        Err(MechanicalPoleError::Budget { .. })
    ));
    let mut calls = 0;
    assert!(matches!(
        prepared.wall_parameter_derivative(&[valid], limits(), &mut || {
            calls += 1;
            calls == 12
        }),
        Err(MechanicalPoleError::Cancelled)
    ));
    assert!(
        prepared
            .wall_parameter_derivative(&[valid], limits(), &mut || false)
            .is_ok()
    );
    let zero = prepared
        .wall_parameter_derivative(&[], limits(), &mut || false)
        .unwrap();
    assert_eq!(zero.rate_per_s_per_parameter, C64::ZERO);
    assert!(zero.coordinate_derivative.iter().all(|&x| x == C64::ZERO));
}
