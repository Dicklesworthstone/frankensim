//! Independent bulk-group checks on the original vector Maxwell carrier.

use super::*;
use super::super::super::kernel as principal;

const TEST_MU: f64 = 1.256_637_062_12e-6;

fn group_input(axis: [f64; 3], incidence: Incidence) -> Input {
    Input {
        omega: core::f64::consts::TAU * C0 / 550e-9,
        ni: 1.0,
        no: 1.5,
        ne: 1.9,
        mui: TEST_MU,
        mut_: TEST_MU,
        axis,
        incidence,
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
    }
}

fn evaluated(
    data: Input,
    group_indices: [f64; 2],
) -> (OrientedCrystalResponse, [Option<[f64; 3]>; 2]) {
    let carrier = solve(data).unwrap();
    let velocities = mode_velocities(
        &carrier,
        data.omega,
        [data.no, data.ne],
        group_indices,
        unit(data.axis).unwrap(),
    )
    .unwrap();
    (carrier, velocities)
}

fn near_group(actual: [f64; 3], expected: [f64; 3]) {
    for j in 0..3 {
        assert!(
            (actual[j] - expected[j]).abs() <= 5e-12 * C0,
            "group component {j}: {actual:?} != {expected:?}"
        );
    }
}

// Independent dispersive Maxwell energy density of one actual solved mode.
// For positive, lossless fixtures v_group=S/u, with
// u=(E*.[d(omega epsilon)/domega].E + mu |H|^2)/4.
// This oracle uses the dimensional E/H fields and constitutive tensor; it
// does not use the production wavevector metric or spectral contraction.
fn maxwell_energy_velocity(
    data: Input,
    mode: OrientedCrystalMode,
    group_indices: [f64; 2],
) -> [f64; 3] {
    let axis_length = data.axis.iter().map(|a| a * a).sum::<f64>().sqrt();
    let axis = data.axis.map(|a| a / axis_length);
    let e = mode.electric_v_m;
    let h = mode.magnetic_a_m;
    let projected_e = (0..3).fold(C64::ZERO, |sum, j| sum + e[j].scale(axis[j]));
    let electric_norm = e.iter().map(|z| z.norm_sq()).sum::<f64>();
    let magnetic_norm = h.iter().map(|z| z.norm_sq()).sum::<f64>();
    let ordinary = data.no * (2.0 * group_indices[0] - data.no) / (data.mut_ * C0 * C0);
    let extraordinary =
        data.ne * (2.0 * group_indices[1] - data.ne) / (data.mut_ * C0 * C0);
    let energy = 0.25
        * (ordinary * electric_norm
            + (extraordinary - ordinary) * projected_e.norm_sq()
            + data.mut_ * magnetic_norm);
    assert!(energy.is_finite() && energy > 0.0);
    [
        (e[1] * h[2].conj() - e[2] * h[1].conj()).re,
        (e[2] * h[0].conj() - e[0] * h[2].conj()).re,
        (e[0] * h[1].conj() - e[1] * h[0].conj()).re,
    ]
    .map(|component| 0.5 * component / energy)
}

#[test]
fn group_velocity_matches_dispersive_maxwell_energy_and_axis_reversal() {
    let data = group_input([0.3, 0.5, -0.7], Incidence::Cosine(0.6));
    let group_indices = [1.8, 2.3];
    let (carrier, velocities) = evaluated(data, group_indices);
    for m in 0..2 {
        assert_eq!(carrier.transmitted[m].wave_kind, UniaxialWaveKind::Propagating);
        near_group(
            velocities[m].unwrap(),
            maxwell_energy_velocity(data, carrier.transmitted[m], group_indices),
        );
    }
    assert!(velocities[0].unwrap()[1].abs() < 1e-5);
    assert!(velocities[1].unwrap()[1].abs() > 1e5);
    assert_eq!(carrier.transmitted[1].phase_direction.unwrap()[1], 0.0);

    let reversed = Input {
        axis: data.axis.map(|a| -a),
        ..data
    };
    let (_, reversed_velocities) = evaluated(reversed, group_indices);
    for m in 0..2 {
        near_group(reversed_velocities[m].unwrap(), velocities[m].unwrap());
    }
}

#[test]
fn oriented_group_reduces_to_principal_modes_and_the_isotropic_speed() {
    let data = group_input([0.6, 0.0, 0.8], Incidence::Cosine(0.6));
    let group_indices = [1.8, 2.3];
    let (_, velocities) = evaluated(data, group_indices);
    for (m, mode) in [principal::UniaxialMode::Ordinary, principal::UniaxialMode::Extraordinary]
        .into_iter()
        .enumerate()
    {
        let carrier = principal::solve(principal::Input {
            omega_rad_s: data.omega,
            incident_index: data.ni,
            incident_mu_h_m: data.mui,
            ordinary_index: data.no,
            extraordinary_index: data.ne,
            crystal_mu_h_m: data.mut_,
            optic_axis: [data.axis[0], data.axis[2]],
            incidence: data.incidence,
            mode,
        })
        .unwrap();
        let velocity = principal::group_velocity(
            carrier,
            data.omega,
            [data.no, data.ne],
            group_indices,
            [data.axis[0], data.axis[2]],
        )
        .unwrap()
        .unwrap();
        near_group(velocities[m].unwrap(), [velocity[0], 0.0, velocity[1]]);
    }

    let isotropic = Input {
        no: 1.6,
        ne: 1.6,
        axis: [1.0, 2.0, 3.0],
        ..data
    };
    let (_, velocities) = evaluated(isotropic, [1.9, 1.9]);
    let qx = isotropic.ni * (1.0_f64 - 0.6 * 0.6).sqrt();
    let qz = (isotropic.no * isotropic.no - qx * qx).sqrt();
    let expected = [
        C0 / 1.9 * qx / isotropic.no,
        0.0,
        C0 / 1.9 * qz / isotropic.no,
    ];
    for velocity in velocities {
        near_group(velocity.unwrap(), expected);
    }
}

#[test]
fn unexcited_evanescent_and_zero_denominator_modes_keep_distinct_carriers() {
    let base = group_input([0.6, 0.0, 0.8], Incidence::Cosine(0.7));
    for (electric, absent) in [
        ([C64::ZERO, C64::ONE], 1),
        ([C64::ONE, C64::ZERO], 0),
    ] {
        let (carrier, velocities) = evaluated(Input { electric, ..base }, [1.8, 2.3]);
        assert_eq!(carrier.transmitted[absent].wave_kind, UniaxialWaveKind::Propagating);
        assert_eq!(carrier.transmitted[absent].electric_v_m, [C64::ZERO; 3]);
        assert_eq!(carrier.transmitted[absent].energy_direction, None);
        assert_eq!(velocities[absent], None);
        assert!(velocities[1 - absent].is_some());
    }

    for (tangent, absent) in [(1.2, [true, false]), (1.8, [true, true])] {
        let data = Input {
            ni: 2.0,
            no: 1.0,
            ne: 1.7,
            axis: [0.0, 1.0, 0.0],
            incidence: Incidence::TangentialIndex(tangent),
            ..base
        };
        let (carrier, velocities) = evaluated(data, [1.2, 2.0]);
        for m in 0..2 {
            assert_eq!(velocities[m].is_none(), absent[m]);
            if absent[m] {
                assert_eq!(carrier.transmitted[m].wave_kind, UniaxialWaveKind::Evanescent);
                assert!(carrier.transmitted[m].normal_wavenumber_rad_m.im < 0.0);
            }
        }
    }

    let (carrier, zero) = evaluated(base, [0.0, 0.0]);
    assert_eq!(zero, [None, None]);
    assert!(carrier.transmitted.iter().all(|mode| {
        mode.wave_kind == UniaxialWaveKind::Propagating && mode.energy_direction.is_some()
    }));
    let (_, positive) = evaluated(base, [1.8, 2.3]);
    let (_, signed) = evaluated(base, [-1.8, -2.3]);
    for m in 0..2 {
        near_group(signed[m].unwrap(), positive[m].unwrap().map(|v| -v));
    }
}

#[test]
fn exact_critical_modes_keep_tangent_group_motion_even_with_negative_phase_kz() {
    for (axis, ni, no, ne, tangent, mode, group_indices) in [
        ([0.3, 0.8, 0.6], 2.0, 1.0, 1.7, 1.0_f64, 0, [1.2, 2.0]),
        // f=(no^2+ne^2)/2=25, so p=5 is exactly critical and qz=-4.8.
        ([1.0, 0.0, 1.0], 6.0, 1.0, 7.0, 5.0, 1, [1.2, 8.0]),
    ] {
        for sign in [-1.0, 1.0] {
            for (p, expected_kind) in [
                (tangent.next_down(), UniaxialWaveKind::Propagating),
                (tangent, UniaxialWaveKind::Critical),
                (tangent.next_up(), UniaxialWaveKind::Evanescent),
            ] {
                let data = Input {
                    ni,
                    no,
                    ne,
                    ..group_input(axis, Incidence::TangentialIndex(sign * p))
                };
                let (carrier, velocities) = evaluated(data, group_indices);
                let transmitted = carrier.transmitted[mode];
                assert_eq!(transmitted.wave_kind, expected_kind);
                if expected_kind == UniaxialWaveKind::Evanescent {
                    assert_eq!(velocities[mode], None);
                    assert!(transmitted.normal_wavenumber_rad_m.im < 0.0);
                    continue;
                }
                let velocity = velocities[mode].unwrap();
                if expected_kind == UniaxialWaveKind::Critical {
                    assert_eq!(transmitted.normal_power_w_m2, 0.0);
                    assert_eq!(velocity[2], 0.0);
                } else {
                    assert!(transmitted.normal_power_w_m2 > 0.0);
                    assert!(velocity[2] > 0.0);
                }
                assert!(velocity[0] * sign > 0.0);
                near_group(velocity, maxwell_energy_velocity(data, transmitted, group_indices));
                if mode == 1 {
                    assert!(sign * transmitted.normal_wavenumber_rad_m.re < 0.0);
                }
            }
        }
    }
}

#[test]
fn unequal_spectral_slopes_at_phase_crossing_require_the_original_eigenbasis() {
    let data = Input {
        no: 1.5,
        ne: 1.5,
        ..group_input([1.0, 1.0, 0.0], Incidence::Cosine(1.0))
    };
    let carrier = solve(data).unwrap();
    // At equal phase indices the existing steady solver chooses TE/TM.
    // These fields mix the two spectral branches for this out-of-plane axis.
    assert_eq!(carrier.transmitted[0].electric_v_m[0], C64::ZERO);
    assert_eq!(carrier.transmitted[1].electric_v_m[1], C64::ZERO);
    assert!(mode_velocities(
        &carrier, data.omega, [data.no, data.ne], [1.6, 2.1], unit(data.axis).unwrap(),
    )
    .is_err());
    let (_, equal_slopes) = evaluated(data, [1.8, 1.8]);
    for velocity in equal_slopes {
        near_group(velocity.unwrap(), [0.0, 0.0, C0 / 1.8]);
    }

    // In the native principal plane TE/TM are the actual spectral branches,
    // including an equal-phase crossing. Along the optic axis both modes
    // depend on the original ordinary slope; across it they remain distinct.
    for (axis, extraordinary_group) in [
        ([1.0, 0.0, 0.0], 2.1),
        ([0.0, 0.0, 1.0], 1.6),
    ] {
        let (_, velocities) = evaluated(Input { axis, ..data }, [1.6, 2.1]);
        near_group(velocities[0].unwrap(), [0.0, 0.0, C0 / 1.6]);
        near_group(velocities[1].unwrap(), [0.0, 0.0, C0 / extraordinary_group]);
    }

    // At a 45-degree principal-axis crossing the extraordinary group index
    // is the arithmetic mean. Preserve the signed one-ulp residual rather
    // than rounding the two ng/n contributions into a false zero.
    let cancelling = [
        f64::from_bits(1.5_f64.to_bits() + 1),
        -f64::from_bits(1.5_f64.to_bits() + 2),
    ];
    let (_, velocities) = evaluated(Input { axis: [1.0, 0.0, 1.0], ..data }, cancelling);
    let velocity = velocities[1].unwrap();
    let expected = C0 / (0.5 * (cancelling[0] + cancelling[1]));
    assert!(velocity[2].is_finite());
    assert_eq!(velocity[0], 0.0);
    assert_eq!(velocity[1], 0.0);
    assert!(
        (velocity[2] - expected).abs() <= 2e-14 * expected.abs(),
        "signed extraordinary cancellation: {velocity:?}, expected {expected:e}"
    );
}
