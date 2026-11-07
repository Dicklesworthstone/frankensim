use super::*;

const MU: f64 = 1.256_637_062_12e-6;
fn input(axis: [f64; 3], cosine: f64) -> Input {
    Input {
        omega: core::f64::consts::TAU * C0 / 550e-9,
        ni: 1.0,
        no: 1.5,
        ne: 1.9,
        mui: MU,
        mut_: MU,
        axis,
        incidence: Incidence::Cosine(cosine),
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
    }
}
fn near(x: f64, y: f64, t: f64) {
    assert!((x - y).abs() <= t, "{x:.17e} != {y:.17e}");
}
fn complex_near(x: C64, y: C64, t: f64) {
    assert!((x - y).abs() <= t, "{x:?} != {y:?}");
}

// Reconstruct Maxwell curls from the dimensional material tensor and actual
// output fields. This reference does not use the production root, inverse
// tensor, dispersion discriminant, mode basis or boundary-system equations.
fn maxwell(data: Input, r: OrientedCrystalResponse) {
    let a = unit(data.axis).unwrap();
    let epsilon_o = data.no * data.no / (data.mut_ * C0 * C0);
    let epsilon_e = data.ne * data.ne / (data.mut_ * C0 * C0);
    let mut et = [C64::ZERO; 3];
    let mut ht = [C64::ZERO; 3];
    for mode in r.transmitted {
        let k = [
            C64::from_re(r.tangential_wavenumber_rad_m),
            C64::ZERO,
            mode.normal_wavenumber_rad_m,
        ];
        let projection = (0..3).fold(C64::ZERO, |s, j| s + mode.electric_v_m[j].scale(a[j]));
        let displacement: Vector = core::array::from_fn(|j| {
            mode.electric_v_m[j].scale(epsilon_o) + projection.scale((epsilon_e - epsilon_o) * a[j])
        });
        let faraday = cross(k, mode.electric_v_m);
        let ampere = cross(k, mode.magnetic_a_m);
        let electric_scale = mode
            .electric_v_m
            .iter()
            .map(|z| z.abs())
            .fold(1.0_f64, f64::max);
        for j in 0..3 {
            complex_near(
                faraday[j],
                mode.magnetic_a_m[j].scale(data.omega * data.mut_),
                4e-8 * electric_scale,
            );
            complex_near(
                ampere[j],
                displacement[j].scale(-data.omega),
                1e-9 * electric_scale,
            );
            et[j] = et[j] + mode.electric_v_m[j];
            ht[j] = ht[j] + mode.magnetic_a_m[j];
        }
        let flux = power(mode.electric_v_m, mode.magnetic_a_m);
        near(
            flux[2],
            mode.normal_power_w_m2,
            1e-17 * electric_scale * electric_scale,
        );
        if mode.wave_kind == UniaxialWaveKind::Propagating && mode.normal_power_w_m2 > 0.0 {
            let direction = unit(flux).unwrap();
            for j in 0..3 {
                near(mode.energy_direction.unwrap()[j], direction[j], 2e-14);
            }
            assert!(mode.normal_wavenumber_rad_m.im == 0.0);
        } else if mode.wave_kind == UniaxialWaveKind::Evanescent {
            assert!(mode.normal_wavenumber_rad_m.im < 0.0);
            assert_eq!(mode.phase_direction, None);
            assert_eq!(mode.energy_direction, None);
        }
    }
    for j in 0..2 {
        complex_near(
            r.incident_electric_v_m[j] + r.reflected_electric_v_m[j],
            et[j],
            3e-14,
        );
        complex_near(
            r.incident_magnetic_a_m[j] + r.reflected_magnetic_a_m[j],
            ht[j],
            3e-16,
        );
    }
    let incident = power(r.incident_electric_v_m, r.incident_magnetic_a_m)[2];
    let reflected = -power(r.reflected_electric_v_m, r.reflected_magnetic_a_m)[2];
    near(incident, r.incident_normal_power_w_m2, 2e-17);
    near(reflected / incident, r.reflected_power_fraction, 3e-14);
    near(
        power(et, ht)[2] / incident,
        r.transmitted_power_fraction,
        3e-14,
    );
    near(
        r.transmitted
            .iter()
            .map(|m| m.normal_power_w_m2)
            .sum::<f64>()
            / incident,
        r.transmitted_power_fraction,
        3e-14,
    );
    near(
        r.reflected_power_fraction + r.transmitted_power_fraction,
        1.0,
        5e-14,
    );
}

#[test]
fn arbitrary_axis_preserves_vector_maxwell_and_normal_power() {
    for axis in [
        [1.0, 2.0, 3.0],
        [-2.0, 1.0, 0.3],
        [0.0, 1.0, 0.0],
        [0.3, 0.4, -0.5],
    ] {
        for cosine in [0.1, 0.5, 0.99, 1.0] {
            for [no, ne] in [[1.5, 1.9], [1.9, 1.5]] {
                let data = Input {
                    no,
                    ne,
                    ..input(axis, cosine)
                };
                maxwell(data, solve(data).unwrap());
            }
        }
    }
}

#[test]
fn azimuthal_optic_axis_generates_cross_polarization_and_two_energy_rays() {
    let mut data = input([1.0, 1.0, 1.0], 0.8);
    data.electric = [C64::ZERO, C64::ONE];
    let r = solve(data).unwrap();
    assert!(r.reflected_electric_v_m[0].abs() > 0.01);
    assert!(r.transmitted.iter().all(|m| m.normal_power_w_m2 > 0.0));
    let eo = r.transmitted[0].energy_direction.unwrap();
    let ee = r.transmitted[1].energy_direction.unwrap();
    assert!(ee[1].abs() > 0.01);
    near(eo[1], 0.0, 1e-15);
    maxwell(data, r);
}

#[test]
fn principal_section_matches_the_existing_independent_decoupled_solver() {
    use crate::plane_wave::uniaxial::kernel;
    for axis in [[0.0, 0.0, 1.0], [0.7, 0.0, -0.2], [1.0, 0.0, 0.0]] {
        for cosine in [0.3, 0.9, 1.0] {
            let mut data = input(axis, cosine);
            data.mut_ = MU * 1.3;
            let r = solve(data).unwrap();
            for (column, mode) in [
                (1, kernel::UniaxialMode::Ordinary),
                (0, kernel::UniaxialMode::Extraordinary),
            ] {
                let expected = kernel::solve(kernel::Input {
                    omega_rad_s: data.omega,
                    incident_index: data.ni,
                    incident_mu_h_m: data.mui,
                    ordinary_index: data.no,
                    extraordinary_index: data.ne,
                    crystal_mu_h_m: data.mut_,
                    optic_axis: [axis[0], axis[2]],
                    incidence: kernel::Incidence::Cosine(cosine),
                    mode,
                })
                .unwrap();
                complex_near(
                    r.reflected_electric_v_m[column],
                    expected.reflected_tangential_e * data.electric[column],
                    2e-14,
                );
            }
            maxwell(data, r);
        }
    }
}

#[test]
fn isotropic_and_optic_axis_degeneracies_are_finite_and_conserve_energy() {
    for axis in [[0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 0.0]] {
        let mut data = input(axis, 1.0);
        data.ne = data.no;
        let r = solve(data).unwrap();
        let reflection = (data.ni - data.no) / (data.ni + data.no);
        for j in 0..2 {
            complex_near(
                r.reflected_electric_v_m[j],
                data.electric[j].scale(reflection),
                1e-15,
            );
        }
        maxwell(data, r);
    }
    let data = input([0.0, 0.0, 1.0], 1.0);
    maxwell(data, solve(data).unwrap());
}

#[test]
fn one_or_both_eigenmodes_can_be_evanescent_without_losing_interface_power() {
    for cosine in [0.7, 0.5, 0.1] {
        let mut data = input([1.0, 2.0, 0.7], cosine);
        data.ni = 2.0;
        data.no = 1.0;
        data.ne = 1.7;
        let r = solve(data).unwrap();
        assert!(
            r.transmitted
                .iter()
                .any(|m| m.wave_kind == UniaxialWaveKind::Evanescent)
        );
        maxwell(data, r);
    }
}

#[test]
fn reversing_axis_and_coherent_source_scaling_preserve_physics() {
    let data = input([0.3, 0.5, -0.7], 0.6);
    let r = solve(data).unwrap();
    let flipped = solve(Input {
        axis: data.axis.map(|v| -v),
        ..data
    })
    .unwrap();
    for j in 0..3 {
        complex_near(
            r.reflected_electric_v_m[j],
            flipped.reflected_electric_v_m[j],
            0.0,
        );
    }
    let factor = C64::new(-1.3, 0.9);
    let scaled = solve(Input {
        electric: data.electric.map(|z| z * factor),
        ..data
    })
    .unwrap();
    near(
        r.reflected_power_fraction,
        scaled.reflected_power_fraction,
        1e-14,
    );
    for j in 0..3 {
        complex_near(
            scaled.reflected_electric_v_m[j],
            r.reflected_electric_v_m[j] * factor,
            3e-15,
        );
        for m in 0..2 {
            complex_near(
                scaled.transmitted[m].electric_v_m[j],
                r.transmitted[m].electric_v_m[j] * factor,
                3e-15,
            );
        }
    }
}

#[test]
fn invalid_and_zero_incident_fields_refuse() {
    let d = input([1.0, 2.0, 3.0], 0.5);
    for cosine in [0.0, -0.5, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            solve(Input {
                incidence: Incidence::Cosine(cosine),
                ..d
            })
            .is_err()
        );
    }
    for axis in [[0.0; 3], [f64::NAN, 0.0, 1.0], [f64::INFINITY, 0.0, 1.0]] {
        assert!(solve(Input { axis, ..d }).is_err());
    }
    assert!(
        solve(Input {
            electric: [C64::ZERO; 2],
            ..d
        })
        .is_err()
    );
}

#[test]
fn signed_adjacent_critical_indices_keep_both_eigenmodes_finite() {
    for (axis, mode, critical) in [([0.3, 0.8, 0.6], 0, 1.0_f64), ([0.0, 1.0, 0.0], 1, 1.7_f64)] {
        let mut data = input(axis, 0.5);
        data.ni = 2.0;
        data.no = 1.0;
        data.ne = 1.7;
        for (p, expected) in [
            (critical.next_down(), UniaxialWaveKind::Propagating),
            (critical, UniaxialWaveKind::Critical),
            (critical.next_up(), UniaxialWaveKind::Evanescent),
        ] {
            for sign in [-1.0, 1.0] {
                data.incidence = Incidence::TangentialIndex(sign * p);
                let r = solve(data).unwrap();
                assert_eq!(r.transmitted[mode].wave_kind, expected);
                maxwell(data, r);
            }
        }
    }
}
