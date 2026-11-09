//! G1/G3/G4: original layer quadrature, coherent beats and relaxed static bias.
use super::*;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::Mandel6;
use fs_plate::shell::viscosity::harmonic::{
    TensorShellMaterialBudget, TensorShellMaterialSpectrumField,
};

const FREQUENCIES: [f64; 2] = [2.4, 3.17];
const BIAS: [f64; 6] = [0.0001, 0.0005, 0.0002, -0.0003, 0.0001, 0.0];
fn motions() -> [[C64; 6]; 2] {
    let q = [
        C64::new(2e-4, 1e-4),
        C64::new(-3e-4, 0.5e-4),
        C64::new(1.3e-4, -0.9e-4),
        C64::new(-2.1e-4, 1.7e-4),
        C64::new(0.8e-4, -0.3e-4),
        C64::ZERO,
    ];
    [q, q.map(|z| z * C64::new(0.35, -0.6))]
}
fn field_budget() -> TensorShellMaterialBudget {
    TensorShellMaterialBudget {
        maximum_points: 24,
        maximum_payload_bytes: 1 << 22,
        maximum_work: 2_000_000_000,
    }
}
fn dissipative_laminate() -> LaminateSection {
    let law = |scale: f64| {
        let base = prony();
        let branch: Mandel6 = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    scale * [1.0, 3.0, 7.0, 2.0, 4.0, 5.0][i]
                } else {
                    0.0
                }
            })
        });
        let kelvin: Mandel6 = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    0.1 + scale * (i + 1) as f64 * 0.04
                } else {
                    0.0
                }
            })
        });
        TensorProny::new_mandel(
            *base.c_eq(),
            vec![
                base.branches()[0].clone(),
                PronyBranch::new(branch, 1.4).unwrap(),
            ],
            0.03,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(kelvin)
        .unwrap()
        .with_band((0.0, 20.0))
        .unwrap()
    };
    LaminateSection::try_new(
        vec![
            layer(LayerLaw::Prony(law(0.4)), 0.08, 4.0, 0.31),
            layer(LayerLaw::Prony(law(0.9)), 0.02, 9.0, -0.47),
        ],
        ReferenceSurface::Offset(0.03),
    )
    .unwrap()
}
fn values(field: &TensorShellMaterialSpectrumField) -> [f64; 8] {
    [
        field.start_stored_energy_j,
        field.end_stored_energy_j,
        field.stored_energy_change_j,
        field.kelvin_heat_j,
        field.maxwell_heat_j,
        field.material_heat_j,
        field.stress_work_j,
        field.energy_balance_residual_j,
    ]
}
fn compare(actual: [f64; 8], expected: [f64; 8], tolerance: f64) {
    let scale = expected.iter().map(|v| v.abs()).fold(1e-30, f64::max);
    for (a, b) in actual.into_iter().zip(expected) {
        assert!(
            (a - b).abs() <= tolerance * scale,
            "{a:e} != {b:e}; scale={scale:e}"
        );
    }
}

// Independent explicit CST field and public DKT curvature, using the physical
// free node's translation/axial-rotation convention. No prepared shell rows or
// material observer feed this reference.
fn section_strain(q: &[C64; 6], bary: [f64; 3]) -> [C64; 6] {
    let (rows, _) = dkt_curvature_quadrature(&[0.0, 0.4, 0.0], &[0.0, 0.0, 0.3], 0).unwrap();
    let midpoint = rows.map(|r| {
        core::array::from_fn::<_, 3, _>(|i| {
            q[2].scale(r[i][6]) + q[3].scale(r[i][8]) - q[4].scale(r[i][7])
        })
    });
    let curvature: [C64; 3] = core::array::from_fn(|i| {
        (0..3).fold(C64::ZERO, |sum, v| {
            sum + (midpoint[v][i] + midpoint[(v + 2) % 3][i] - midpoint[(v + 1) % 3][i])
                .scale(bary[v])
        })
    });
    [
        C64::ZERO,
        q[1].scale(1.0 / 0.3),
        q[0].scale(1.0 / 0.3),
        curvature[0],
        curvature[1],
        curvature[2],
    ]
}
fn in_plane(q: &[C64; 6], bary: [f64; 3], z: f64, angle: f64) -> [C64; 3] {
    let section = section_strain(q, bary);
    let e: [C64; 3] = core::array::from_fn(|i| section[i] + section[i + 3].scale(z));
    let (s, c) = angle.sin_cos();
    [
        e[0].scale(c * c) + e[1].scale(s * s) + e[2].scale(c * s),
        e[0].scale(s * s) + e[1].scale(c * c) - e[2].scale(c * s),
        e[0].scale(-2.0 * c * s) + e[1].scale(2.0 * c * s) + e[2].scale(c * c - s * s),
    ]
}

// The fixture is orthotropic in each layer. Its complete complex normal
// constraint has the independent scalar solution ezz=-(Czx*ex+Czy*ey)/Czz.
// Nonproportional branches/Kelvin terms make that ratio complex and dispersive.
fn full_strain(law: &TensorProny, omega: f64, e: [C64; 3]) -> [C64; 6] {
    let mut normal = core::array::from_fn::<_, 3, _>(|j| {
        C64::new(
            law.c_eq()[2][j],
            -omega * law.kv_pa_s().map_or(0.0, |k| k[2][j]),
        )
    });
    for branch in law.branches() {
        let x = omega * branch.tau_s;
        let weight = C64::new(x * x, -x).scale(1.0 / (1.0 + x * x));
        for (j, c) in normal.iter_mut().enumerate() {
            *c = *c + weight.scale(branch.c_pa[2][j]);
        }
    }
    [
        e[0],
        e[1],
        -(normal[0] * e[0] + normal[1] * e[1]) / normal[2],
        e[2].scale(core::f64::consts::FRAC_1_SQRT_2),
        C64::ZERO,
        C64::ZERO,
    ]
}
fn apply(matrix: &Mandel6, value: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * value[j]).sum())
}
fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn real_fields(law: &TensorProny, strains: &[[C64; 6]; 2], bias: &[f64; 6], t: f64) -> [f64; 4] {
    let mut strain = *bias;
    let mut rate = [0.0; 6];
    for (&w, e) in FREQUENCIES.iter().zip(strains) {
        let (s, c) = (w * t).sin_cos();
        for i in 0..6 {
            strain[i] += e[i].re * c + e[i].im * s;
            rate[i] += w * (e[i].im * c - e[i].re * s);
        }
    }
    let mut stress = apply(law.c_eq(), &strain);
    let mut stored = 0.5 * dot(&strain, &stress);
    let viscous = apply(law.kv_pa_s().unwrap(), &rate);
    let kelvin = dot(&rate, &viscous);
    for i in 0..6 {
        stress[i] += viscous[i];
    }
    let mut maxwell = 0.0;
    for branch in law.branches() {
        let mut spring = [0.0; 6];
        for (&w, e) in FREQUENCIES.iter().zip(strains) {
            let (s, c) = (w * t).sin_cos();
            let x = w * branch.tau_s;
            for i in 0..6 {
                spring[i] += ((x * x * e[i].re + x * e[i].im) * c
                    + (x * x * e[i].im - x * e[i].re) * s)
                    / (1.0 + x * x);
            }
        }
        let branch_stress = apply(&branch.c_pa, &spring);
        let square = dot(&spring, &branch_stress);
        stored += 0.5 * square;
        maxwell += square / branch.tau_s;
        for i in 0..6 {
            stress[i] += branch_stress[i];
        }
    }
    [stored, kelvin, maxwell, dot(&stress, &rate)]
}

// A DIFFERENT positive spatial rule: triangle vertices+centroid (degree two),
// layer faces+center (Simpson), followed by real-time Simpson integration.
fn reference(
    section: &LaminateSection,
    q: &[[C64; 6]; 2],
    bias: &[f64; 6],
    window: HarmonicWindow,
) -> [f64; 8] {
    let triangle = [
        ([1.0, 0.0, 0.0], 1.0 / 12.0),
        ([0.0, 1.0, 0.0], 1.0 / 12.0),
        ([0.0, 0.0, 1.0], 1.0 / 12.0),
        ([1.0 / 3.0; 3], 0.75),
    ];
    let mut total = [0.0; 8];
    for (layer, geometry) in section.layers().iter().zip(section.layer_integrals()) {
        let LayerLaw::Prony(law) = &layer.law else {
            unreachable!()
        };
        for (bary, area_weight) in triangle {
            for (fraction, height_weight) in [(0.0, 1.0 / 6.0), (0.5, 2.0 / 3.0), (1.0, 1.0 / 6.0)]
            {
                let z = geometry.z_bottom_m + fraction * layer.thickness_m;
                let strains = core::array::from_fn(|i| {
                    full_strain(
                        law,
                        FREQUENCIES[i],
                        in_plane(&q[i], bary, z, layer.orientation_rad),
                    )
                });
                let bias = full_strain(
                    law,
                    0.0,
                    in_plane(&bias.map(C64::from_re), bary, z, layer.orientation_rad),
                )
                .map(|v| v.re);
                let begin = real_fields(law, &strains, &bias, window.start_s())[0];
                let end =
                    real_fields(law, &strains, &bias, window.start_s() + window.duration_s())[0];
                let n = 2048;
                let mut integral = [0.0; 3];
                for i in 0..=n {
                    let weight = if i == 0 || i == n {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    let fields = real_fields(
                        law,
                        &strains,
                        &bias,
                        window.start_s() + window.duration_s() * i as f64 / n as f64,
                    );
                    for j in 0..3 {
                        integral[j] += weight * fields[j + 1];
                    }
                }
                let [kelvin, maxwell, work] =
                    integral.map(|v| v * window.duration_s() / (3.0 * n as f64));
                let volume = 0.06 * area_weight * layer.thickness_m * height_weight;
                for (total, value) in total.iter_mut().zip([
                    begin,
                    end,
                    end - begin,
                    kelvin,
                    maxwell,
                    kelvin + maxwell,
                    work,
                    work - (end - begin) - kelvin - maxwell,
                ]) {
                    *total += volume * value;
                }
            }
        }
    }
    total
}

#[test]
fn offset_rotated_layers_and_relaxed_bias_match_independent_real_field_integrals() {
    with_cx(|cx| {
        let section = dissipative_laminate();
        let shell = source(section.clone(), CstDktMass::Consistent, cx).unwrap();
        let observer = shell
            .prepare_material_spectrum_observer(&FREQUENCIES, field_budget(), cx)
            .unwrap();
        assert!(core::ptr::eq(observer.source(), &shell));
        let q = motions();
        let displacements = q.each_ref().map(|v| v.as_slice());
        let window = HarmonicWindow::new(0.37, 1.31).unwrap();
        let mut heat = Vec::new();
        for bias in [[0.0; 6], BIAS] {
            let field = observer
                .observe_window_with_bias(&displacements, &bias, window, cx)
                .unwrap();
            compare(values(&field), reference(&section, &q, &bias, window), 3e-9);
            assert_eq!(field.points.len(), 12);
            let mut volumes = [0.0; 2];
            for point in &field.points {
                volumes[point.layer] += point.volume_m3;
                let layer = &section.layers()[point.layer];
                for (tone, strain) in point.in_plane_strain.iter().enumerate() {
                    let expected = in_plane(
                        &q[tone],
                        point.triangle_barycentric,
                        point.thickness_position_m,
                        layer.orientation_rad,
                    );
                    let scale = expected.iter().map(|z| z.abs()).fold(1e-20, f64::max);
                    for (actual, expected) in strain.iter().zip(expected) {
                        assert!((*actual - expected).abs() < 3e-12 * scale);
                    }
                }
            }
            real_near(volumes[0], 0.06 * 0.08, 3e-13);
            real_near(volumes[1], 0.06 * 0.02, 3e-13);
            heat.push(field.material_heat_j);
        }
        assert_eq!(
            heat[0], heat[1],
            "a relaxed static strain creates no dashpot heat"
        );
        let first = observer.observe_window(&displacements, window, cx).unwrap();
        let split = 0.61;
        let left = observer
            .observe_window(
                &displacements,
                HarmonicWindow::new(window.start_s(), split).unwrap(),
                cx,
            )
            .unwrap();
        let right = observer
            .observe_window(
                &displacements,
                HarmonicWindow::new(window.start_s() + split, window.duration_s() - split).unwrap(),
                cx,
            )
            .unwrap();
        for (a, b) in [
            (
                first.material_heat_j,
                left.material_heat_j + right.material_heat_j,
            ),
            (
                first.stress_work_j,
                left.stress_work_j + right.stress_work_j,
            ),
            (
                first.stored_energy_change_j,
                left.stored_energy_change_j + right.stored_energy_change_j,
            ),
        ] {
            real_near(a, b, 3e-11);
        }
    });
}

#[test]
fn one_tone_period_matches_native_mean_fields_and_empty_spectrum_keeps_static_storage() {
    with_cx(|cx| {
        let shell = source(laminate(), CstDktMass::Consistent, cx).unwrap();
        let w = 2.4;
        let q = motions()[0];
        let prepared = shell.prepare_harmonic_response(w, controls(), cx).unwrap();
        let original = prepared
            .prepare_material_observer(field_budget(), cx)
            .unwrap()
            .observe(&q, cx)
            .unwrap();
        let spectrum = shell
            .prepare_material_spectrum_observer(&[w], field_budget(), cx)
            .unwrap();
        let duration = 2.0 * core::f64::consts::PI / w;
        let observed = spectrum
            .observe_window(&[&q], HarmonicWindow::new(0.13, duration).unwrap(), cx)
            .unwrap();
        real_near(
            observed.material_heat_j,
            duration * original.mean_material_heat_w,
            3e-11,
        );
        real_near(
            observed.kelvin_heat_j,
            duration * original.mean_kelvin_heat_w,
            3e-11,
        );
        real_near(
            observed.maxwell_heat_j,
            duration * original.mean_maxwell_heat_w,
            3e-11,
        );
        assert!(observed.stored_energy_change_j.abs() < 3e-11 * observed.start_stored_energy_j);
        let empty = shell
            .prepare_material_spectrum_observer(&[], field_budget(), cx)
            .unwrap();
        let window = HarmonicWindow::new(0.13, 0.7).unwrap();
        let zero = empty.observe_window(&[], window, cx).unwrap();
        assert_eq!(values(&zero), [0.0; 8]);
        let biased = empty
            .observe_window_with_bias(&[], &BIAS, window, cx)
            .unwrap();
        assert!(biased.start_stored_energy_j > 0.0);
        assert_eq!(biased.start_stored_energy_j, biased.end_stored_energy_j);
        assert_eq!(biased.stored_energy_change_j, 0.0);
        assert_eq!(biased.material_heat_j, 0.0);
        assert_eq!(biased.stress_work_j, 0.0);
    });
}

#[test]
fn coherent_line_cancellation_clock_shift_and_drilling_keep_physical_material_semantics() {
    with_cx(|cx| {
        let shell = source(dissipative_laminate(), CstDktMass::Consistent, cx).unwrap();
        let q = motions();
        let refs = q.each_ref().map(|v| v.as_slice());
        let window = HarmonicWindow::new(0.37, 1.31).unwrap();
        let observer = shell
            .prepare_material_spectrum_observer(&FREQUENCIES, field_budget(), cx)
            .unwrap();
        let original = observer.observe_window(&refs, window, cx).unwrap();
        let shift = 4.3;
        let shifted: [[C64; 6]; 2] = core::array::from_fn(|i| {
            let (s, c) = (FREQUENCIES[i] * shift).sin_cos();
            q[i].map(|z| z * C64::new(c, s))
        });
        let shifted_refs = shifted.each_ref().map(|v| v.as_slice());
        let moved = observer
            .observe_window(
                &shifted_refs,
                HarmonicWindow::new(window.start_s() + shift, window.duration_s()).unwrap(),
                cx,
            )
            .unwrap();
        compare(values(&moved), values(&original), 3e-11);
        let bank = shell
            .prepare_material_spectrum_observer(&[2.4; 3], field_budget(), cx)
            .unwrap();
        let huge = q[0].map(|z| z.scale(1e12));
        let minus = huge.map(|z| -z);
        let cancelled = bank
            .observe_window(&[&huge, &q[0], &minus], window, cx)
            .unwrap();
        let single = shell
            .prepare_material_spectrum_observer(&[2.4], field_budget(), cx)
            .unwrap()
            .observe_window(&[&q[0]], window, cx)
            .unwrap();
        compare(values(&cancelled), values(&single), 3e-10);
        let mut drilling = [C64::ZERO; 6];
        drilling[5] = C64::new(0.1, -0.2);
        let pure = observer
            .observe_window(&[&drilling, &drilling], window, cx)
            .unwrap();
        assert!(shell.model().k.get(5, 5) > 0.0);
        assert_eq!(values(&pure), [0.0; 8]);
    });
}

fn peak(strain: [C64; 6]) -> f64 {
    let a: f64 = strain.iter().map(|z| z.re * z.re).sum();
    let b: f64 = strain.iter().map(|z| z.im * z.im).sum();
    let c: f64 = strain.iter().map(|z| z.re * z.im).sum();
    (0.5 * (a + b + ((a - b) * (a - b) + 4.0 * c * c).sqrt())).sqrt()
}

#[test]
fn full_layer_faces_control_admission_before_interior_quadrature_and_refusals_preserve_reuse() {
    with_cx(|cx| {
        let section = dissipative_laminate();
        let shell = source(section.clone(), CstDktMass::Consistent, cx).unwrap();
        for budget in [
            TensorShellMaterialBudget {
                maximum_points: 11,
                ..field_budget()
            },
            TensorShellMaterialBudget {
                maximum_payload_bytes: 0,
                ..field_budget()
            },
            TensorShellMaterialBudget {
                maximum_work: 0,
                ..field_budget()
            },
        ] {
            assert!(
                shell
                    .prepare_material_spectrum_observer(&[2.4], budget, cx)
                    .is_err()
            );
        }
        let observer = shell
            .prepare_material_spectrum_observer(&[2.4], field_budget(), cx)
            .unwrap();
        let window = HarmonicWindow::new(0.37, 1.31).unwrap();
        let mut q = [C64::ZERO; 6];
        q[3] = C64::ONE;
        let (mut boundary, mut interior) = (0.0_f64, 0.0_f64);
        for (layer, geometry) in section.layers().iter().zip(section.layer_integrals()) {
            let LayerLaw::Prony(law) = &layer.law else {
                unreachable!()
            };
            for bary in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
                for z in [geometry.z_bottom_m, geometry.z_top_m] {
                    boundary = boundary.max(peak(full_strain(
                        law,
                        2.4,
                        in_plane(&q, bary, z, layer.orientation_rad),
                    )));
                }
            }
            let center = 0.5 * (geometry.z_bottom_m + geometry.z_top_m);
            for bary in [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] {
                for sign in [-1.0, 1.0] {
                    let z = center + sign * layer.thickness_m / (2.0 * 3.0_f64.sqrt());
                    interior = interior.max(peak(full_strain(
                        law,
                        2.4,
                        in_plane(&q, bary, z, layer.orientation_rad),
                    )));
                }
            }
        }
        assert!(
            boundary > 1.01 * interior,
            "fixture must hide a face violation from interior samples"
        );
        let scale = 0.03 / (boundary * interior).sqrt();
        assert!(scale * interior < 0.03 && scale * boundary > 0.03);
        assert!(
            observer
                .observe_window(&[&q.map(|z| z.scale(scale))], window, cx)
                .is_err()
        );
        let safe = motions()[0];
        let accepted = observer.observe_window(&[&safe], window, cx).unwrap();
        assert!(observer.observe_window(&[&safe[..5]], window, cx).is_err());
        assert!(
            observer
                .observe_window_with_bias(&[&safe], &[0.0; 5], window, cx)
                .is_err()
        );
        assert!(
            observer
                .observe_window_with_bias(&[&safe], &[0.1; 6], window, cx)
                .is_err()
        );
        let shaped =
            HarmonicWindow::with_profile(0.37, 1.31, HarmonicWindowProfile::Triangular).unwrap();
        assert!(observer.observe_window(&[&safe], shaped, cx).is_err());
        assert_eq!(
            accepted,
            observer.observe_window(&[&safe], window, cx).unwrap()
        );
    });
}

#[test]
fn cancellation_keeps_the_complete_source_available_to_a_fresh_scope() {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let key = fs_exec::StreamKey {
            seed: 59,
            kernel_id: 38,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(
            &gate,
            arena,
            key,
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let shell = source(dissipative_laminate(), CstDktMass::Consistent, &cx).unwrap();
        let observer = shell
            .prepare_material_spectrum_observer(&FREQUENCIES, field_budget(), &cx)
            .unwrap();
        let q = motions();
        let displacements = q.each_ref().map(|v| v.as_slice());
        let window = HarmonicWindow::new(0.37, 1.31).unwrap();
        let accepted = observer
            .observe_window_with_bias(&displacements, &BIAS, window, &cx)
            .unwrap();
        gate.request();
        assert!(matches!(
            observer.observe_window_with_bias(&displacements, &BIAS, window, &cx),
            Err(TensorShellHarmonicError::Plate(
                fs_plate::PlateError::Cancelled
            ))
        ));
        assert!(matches!(
            shell.prepare_material_spectrum_observer(&FREQUENCIES, field_budget(), &cx),
            Err(TensorShellHarmonicError::Plate(
                fs_plate::PlateError::Cancelled
            ))
        ));
        let fresh = CancelGate::new_clock_free();
        let cx = Cx::new(
            &fresh,
            arena,
            key,
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        assert_eq!(
            accepted,
            observer
                .observe_window_with_bias(&displacements, &BIAS, window, &cx)
                .unwrap()
        );
    });
}
