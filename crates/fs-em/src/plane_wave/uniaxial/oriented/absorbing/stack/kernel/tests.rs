//! G0/G1/G3 controls: independent optics, native-film parity and layer invariance.
use super::*;
use crate::plane_wave::uniaxial::oriented::kernel::unit;
use fs_math::det;

const MU: f64 = 1.256_637_062_12e-6;

fn fixture(ni: f64, nt: f64, rows: &[([C64; 2], f64, [f64; 3])]) -> Frequency {
    let mut entrance = 0.0;
    let layers = rows
        .iter()
        .map(|&(indices, thickness, axis)| {
            let exit = entrance + thickness;
            let layer = ResolvedLayer {
                indices,
                mu: MU,
                axis: unit(axis).unwrap(),
                thickness,
                entrance,
                exit,
                loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (MU * C0 * C0)),
            };
            entrance = exit;
            layer
        })
        .collect();
    Frequency {
        omega: C0,
        ni,
        nt,
        mui: MU,
        mut_: MU,
        layers,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}

fn zclose(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}; tolerance {tolerance:e}"
    );
}

fn answer(data: &Frequency, p: f64, electric: [C64; 2]) -> AbsorbingStackResponse {
    solve(
        data,
        Incidence::TangentialIndex(p),
        electric,
        [1.0, 0.0],
        &mut || false,
    )
    .unwrap()
}

// A short-stack scalar characteristic-matrix oracle. Its growing inverse
// phase is acceptable only in this bounded reference, never in production.
fn characteristic(data: &Frequency, p: f64, tm: bool) -> (C64, C64, f64) {
    let admittance = |n: C64, mu: f64| {
        let mut q = (n * n - C64::from_re(p * p)).sqrt();
        if q.im > 0.0 {
            q = -q;
        }
        if tm {
            n * n / q.scale(mu / MU)
        } else {
            q.scale(MU / mu)
        }
    };
    let mut matrix = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    for layer in &data.layers {
        let n = layer.indices[0];
        let mut q = (n * n - C64::from_re(p * p)).sqrt();
        if q.im > 0.0 {
            q = -q;
        }
        let y = admittance(n, layer.mu);
        let phase = C64::new(
            det::cos(q.re * layer.thickness),
            -det::sin(q.re * layer.thickness),
        )
        .scale(det::exp(q.im * layer.thickness));
        let cosine = (phase + C64::ONE / phase).scale(0.5);
        let i_sine = (C64::ONE / phase - phase).scale(0.5);
        let next = [[cosine, i_sine / y], [i_sine * y, cosine]];
        matrix = core::array::from_fn(|j| {
            core::array::from_fn(|k| matrix[j][0] * next[0][k] + matrix[j][1] * next[1][k])
        });
    }
    let yi = admittance(C64::from_re(data.ni), data.mui);
    let yt = admittance(C64::from_re(data.nt), data.mut_);
    let a = matrix[0][0] + matrix[0][1] * yt;
    let b = matrix[1][0] + matrix[1][1] * yt;
    (
        (a * yi - b) / (a * yi + b),
        yi.scale(2.0) / (a * yi + b),
        yt.re / yi.re,
    )
}

#[test]
fn single_layer_sweep_agrees_with_native_absorbing_film_including_lossless_and_evanescent_exit() {
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    for indices in [
        [C64::new(1.5, -0.07), C64::new(2.1, -0.23)],
        [C64::from_re(1.5), C64::from_re(2.1)],
    ] {
        for axis in [[0.4, 0.7, 0.6], [0.6, 0.8, 0.0]] {
            let data = fixture(1.2, 1.0, &[(indices, 0.9, axis)]);
            for p in [0.0, 0.7, 1.05] {
                let layer = data.layers[0];
                let stack = answer(&data, p, electric);
                let slab = crate::plane_wave::uniaxial::oriented::absorbing::slab::kernel::solve(
                    Input {
                        omega: data.omega,
                        ni: data.ni,
                        indices,
                        mui: data.mui,
                        mut_: layer.mu,
                        axis: layer.axis,
                        incidence: Incidence::TangentialIndex(p),
                        electric,
                        loss_f_m: layer.loss_f_m,
                    },
                    data.nt,
                    data.mut_,
                    layer.thickness,
                )
                .unwrap();
                for j in 0..3 {
                    zclose(
                        stack.reflected_electric_v_m[j],
                        slab.reflected_electric_v_m[j],
                        6e-14,
                    );
                    zclose(
                        stack.transmitted_electric_v_m[j],
                        slab.transmitted_electric_v_m[j],
                        6e-14,
                    );
                }
                close(
                    stack.reflected_power_fraction,
                    slab.reflected_power_fraction,
                    5e-14,
                );
                close(
                    stack.transmitted_power_fraction,
                    slab.transmitted_power_fraction,
                    5e-14,
                );
                close(
                    stack.absorbed_power_fraction,
                    slab.absorbed_power_fraction,
                    5e-14,
                );
                let a = stack
                    .field_at(StackRegion::Layer(0), [0.2, -0.1, 0.37])
                    .unwrap();
                let b = slab
                    .field_at(
                        crate::plane_wave::uniaxial::oriented::slab::SlabRegion::Crystal,
                        [0.2, -0.1, 0.37],
                    )
                    .unwrap();
                for j in 0..3 {
                    zclose(a.field.electric_v_m[j], b.field.electric_v_m[j], 1e-13);
                    zclose(a.field.magnetic_a_m[j], b.field.magnetic_a_m[j], 1e-15);
                }
                close(
                    a.absorbed_power_density_w_m3,
                    b.absorbed_power_density_w_m3,
                    3e-17,
                );
            }
        }
    }
}

#[test]
fn absorbing_isotropic_multilayers_match_independent_complex_characteristic_matrices() {
    let mut data = fixture(
        1.3,
        1.2,
        &[
            ([C64::new(1.6, -0.07); 2], 0.7, [0.2, 0.5, 0.8]),
            ([C64::new(2.1, -0.21); 2], 0.43, [1.0, 0.0, 0.0]),
            ([C64::new(1.4, -0.02); 2], 0.29, [0.0, 1.0, 0.0]),
        ],
    );
    data.layers[1].mu *= 1.4;
    data.layers[1].loss_f_m = data.layers[1].loss_f_m.map(|v| v / 1.4);
    data.mut_ *= 0.9;
    for p in [-0.8, 0.0, 0.65, 1.15] {
        for tm in [false, true] {
            let j = if tm { 0 } else { 1 };
            let mut electric = [C64::ZERO; 2];
            electric[j] = C64::new(0.6, -0.3);
            let response = answer(&data, p, electric);
            let (r, t, ratio) = characteristic(&data, p, tm);
            zclose(response.reflected_electric_v_m[j], electric[j] * r, 4e-14);
            zclose(response.transmitted_electric_v_m[j], electric[j] * t, 4e-14);
            close(response.reflected_power_fraction, r.norm_sq(), 5e-14);
            close(
                response.transmitted_power_fraction,
                ratio * t.norm_sq(),
                5e-14,
            );
            close(
                response.absorbed_power_fraction,
                1.0 - r.norm_sq() - ratio * t.norm_sq(),
                5e-14,
            );
            let mut heat = 0.0;
            for layer in &response.layers {
                let direct = layer.absorbed_power_w(1.0, 0.0, layer.thickness_m).unwrap();
                let entering = layer.face(false).unwrap().poynting_w_m2[2];
                let leaving = layer.face(true).unwrap().poynting_w_m2[2];
                assert!(direct > 0.0);
                close(direct, entering - leaving, 2e-16);
                heat += direct;
            }
            close(
                heat,
                response.incident_normal_power_w_m2 * response.absorbed_power_fraction,
                3e-16,
            );
        }
    }
}

#[test]
fn independent_tilted_axes_preserve_all_faces_heat_subdivision_and_global_azimuth() {
    let a = [C64::new(1.3, -0.08), C64::new(2.1, -0.19)];
    let b = [C64::new(1.7, -0.06), C64::new(2.3, -0.27)];
    let c = [C64::from_re(1.4), C64::from_re(1.8)];
    let axes = [[0.3, 0.7, 0.5], [0.6, -0.2, 0.5], [-0.4, 0.5, 0.6]];
    let rows = [(a, 0.8, axes[0]), (b, 0.6, axes[1]), (c, 0.3, axes[2])];
    let data = fixture(1.8, 1.4, &rows);
    let split = fixture(
        1.8,
        1.4,
        &[
            (a, 0.3, axes[0]),
            (a, 0.5, axes[0]),
            (b, 0.6, axes[1]),
            (c, 0.3, axes[2]),
        ],
    );
    let electric = [C64::new(0.7, 0.3), C64::new(-0.2, 0.6)];
    let response = answer(&data, 0.9, electric);
    let refined = answer(&split, 0.9, electric);
    for j in 0..3 {
        zclose(
            response.reflected_electric_v_m[j],
            refined.reflected_electric_v_m[j],
            5e-14,
        );
        zclose(
            response.transmitted_electric_v_m[j],
            refined.transmitted_electric_v_m[j],
            5e-14,
        );
    }
    assert!(
        (response.layers[0].forward[1].normal_wavenumber_rad_m
            + response.layers[0].backward[1].normal_wavenumber_rad_m)
            .abs()
            > 0.1
    );
    let rotation = |v: [f64; 3]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]];
    let rotate = |v: Vector| {
        [
            v[0].scale(0.6) - v[1].scale(0.8),
            v[0].scale(0.8) + v[1].scale(0.6),
            v[2],
        ]
    };
    let rotated_rows = rows.map(|(indices, depth, axis)| (indices, depth, rotation(axis)));
    let rotated_data = fixture(1.8, 1.4, &rotated_rows);
    let rotated_e = rotate([electric[0], electric[1], C64::ZERO]);
    let turned = solve(
        &rotated_data,
        Incidence::TangentialIndex(0.9),
        [rotated_e[0], rotated_e[1]],
        [0.6, 0.8],
        &mut || false,
    )
    .unwrap();
    for (z, original, subdivided) in [(0.15, 0, 0), (0.6, 0, 1), (1.1, 1, 2), (1.55, 2, 3)] {
        let field = response
            .field_at(StackRegion::Layer(original), [0.2, 0.0, z])
            .unwrap();
        let split_field = refined
            .field_at(StackRegion::Layer(subdivided), [0.2, 0.0, z])
            .unwrap();
        let rotated_field = turned
            .field_at(StackRegion::Layer(original), [0.12, 0.16, z])
            .unwrap();
        for j in 0..3 {
            zclose(
                field.field.electric_v_m[j],
                split_field.field.electric_v_m[j],
                7e-14,
            );
            zclose(
                rotated_field.field.electric_v_m[j],
                rotate(field.field.electric_v_m)[j],
                9e-14,
            );
            zclose(
                rotated_field.field.magnetic_a_m[j],
                rotate(field.field.magnetic_a_m)[j],
                1e-15,
            );
        }
        close(
            field.absorbed_power_density_w_m3,
            split_field.absorbed_power_density_w_m3,
            6e-17,
        );
        close(
            field.absorbed_power_density_w_m3,
            rotated_field.absorbed_power_density_w_m3,
            6e-17,
        );
    }
    for result in [&response, &refined, &turned] {
        let mut left = result
            .plain_field_at(StackRegion::TotalIncidentSide, [0.0; 3])
            .unwrap();
        let mut total_heat = 0.0;
        for layer in &result.layers {
            let front = layer.face(false).unwrap();
            let back = layer.face(true).unwrap();
            for j in 0..2 {
                zclose(left.electric_v_m[j], front.electric_v_m[j], 8e-14);
                zclose(left.magnetic_a_m[j], front.magnetic_a_m[j], 1e-15);
            }
            let heat = layer.absorbed_power_w(0.7, 0.0, layer.thickness_m).unwrap();
            close(
                heat,
                0.7 * (front.poynting_w_m2[2] - back.poynting_w_m2[2]),
                5e-16,
            );
            total_heat += heat;
            left = back;
        }
        let exit = result
            .plain_field_at(StackRegion::Transmitted, [0.0, 0.0, result.thickness_m])
            .unwrap();
        for j in 0..2 {
            zclose(left.electric_v_m[j], exit.electric_v_m[j], 8e-14);
            zclose(left.magnetic_a_m[j], exit.magnetic_a_m[j], 1e-15);
        }
        close(
            total_heat,
            0.7 * result.incident_normal_power_w_m2 * result.absorbed_power_fraction,
            7e-16,
        );
    }
}

#[test]
fn opaque_and_evanescent_stacks_preserve_passivity_grazing_ports_and_bounded_cancellation() {
    let indices = [C64::new(1.5, -0.2), C64::new(2.0, -0.5)];
    let rows: Vec<_> = (0..32)
        .map(|j| {
            (
                indices,
                1e4,
                if j % 2 == 0 {
                    [0.4, 0.7, 0.6]
                } else {
                    [-0.3, 0.1, 0.8]
                },
            )
        })
        .collect();
    let data = fixture(1.2, 1.0, &rows);
    let electric = [C64::new(0.6, 0.3), C64::ONE];
    let response = answer(&data, 0.7, electric);
    let first = data.layers[0];
    let halfspace = crate::plane_wave::uniaxial::oriented::absorbing::kernel::solve(Input {
        omega: data.omega,
        ni: data.ni,
        indices,
        mui: data.mui,
        mut_: first.mu,
        axis: first.axis,
        incidence: Incidence::TangentialIndex(0.7),
        electric,
        loss_f_m: first.loss_f_m,
    })
    .unwrap();
    close(
        response.reflected_power_fraction,
        halfspace.reflected_power_fraction,
        4e-14,
    );
    close(
        response.absorbed_power_fraction,
        halfspace.transmitted_power_fraction,
        4e-14,
    );
    assert_eq!(response.transmitted_power_fraction, 0.0);
    for (j, layer) in response.layers.iter().enumerate() {
        let field = response
            .field_at(
                StackRegion::Layer(j),
                [0.0, 0.0, layer.entrance_z_m + 0.5 * layer.thickness_m],
            )
            .unwrap();
        assert!(field.field.electric_v_m.iter().all(|&v| finite(v)));
        assert!(field.absorbed_power_density_w_m3 >= 0.0);
    }
    for depth in [0.2, 5.0, 40.0] {
        let barrier = fixture(
            2.0,
            2.1,
            &[
                (
                    [C64::from_re(1.0), C64::from_re(1.2)],
                    depth,
                    [0.6, 0.4, 0.3],
                ),
                (
                    [C64::from_re(1.1), C64::from_re(1.25)],
                    depth,
                    [0.2, -0.3, 0.8],
                ),
            ],
        );
        let result = answer(&barrier, 1.5, electric);
        assert!(result.transmitted_power_fraction > 0.0);
        assert_eq!(result.absorbed_power_fraction, 0.0);
        close(
            result.reflected_power_fraction + result.transmitted_power_fraction,
            1.0,
            2e-13,
        );
    }
    let grazing = fixture(
        1.2,
        1.2,
        &[
            ([C64::new(1.7, -0.1); 2], 0.4, [1.0, 2.0, 3.0]),
            ([C64::new(2.1, -0.2); 2], 0.6, [3.0, 2.0, 1.0]),
        ],
    );
    let result = solve(
        &grazing,
        Incidence::Cosine(1e-9),
        [C64::ZERO, C64::ONE],
        [1.0, 0.0],
        &mut || false,
    )
    .unwrap();
    close(result.incident_normal_wavenumber_rad_m, 1.2e-9, 1e-23);
    close(result.exit_normal_wavenumber_rad_m.re, 1.2e-9, 1e-23);
    assert_eq!(result.exit_wave_kind, UniaxialWaveKind::Propagating);
    for stop in [1, 4, 37] {
        let mut calls = 0;
        let result = solve(
            &data,
            Incidence::TangentialIndex(0.7),
            electric,
            [1.0, 0.0],
            &mut || {
                calls += 1;
                calls == stop
            },
        );
        assert!(matches!(
            result,
            Err(crate::plane_wave::angular_spectrum::AngularSpectrumError::Cancelled)
        ));
    }
    let small = workspace_bytes(10, 2).unwrap();
    let large = workspace_bytes(20, 2).unwrap();
    assert!(large > small && large < 2 * small);
    assert!(workspace_bytes(usize::MAX, 2).is_err());
}
