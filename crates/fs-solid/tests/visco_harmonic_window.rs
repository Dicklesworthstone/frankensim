//! Physical interval heat/work on original tensor-Prony solid quadrature.
//! Material coefficients and displacements are synthetic numerical fixtures.

use fs_exec::Cx;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseError;
use fs_solid::linear::Formulation;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

const TET: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const HEX: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];

fn with_cx(run: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 81,
                    kernel_id: 820,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn tensor(diagonal: f64, v: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| v[i] * v[j] + if i == j { diagonal } else { 0.0 })
    })
}

fn law() -> TensorProny {
    TensorProny::new_mandel(
        tensor(80.0, [2.0, -1.0, 0.5, 1.5, -0.4, 0.2]),
        vec![
            PronyBranch::new(tensor(0.0, [1.0, 2.0, -0.5, 0.8, -0.4, 0.3]), 0.08).unwrap(),
            PronyBranch::new(tensor(4.0, [-0.2, 0.3, 1.0, -0.6, 0.7, 0.8]), 0.7).unwrap(),
        ],
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(tensor(0.04, [0.2, -0.1, 0.4, 0.05, -0.3, 0.2]))
    .unwrap()
    .with_band((0.1, 100.0))
    .unwrap()
}

fn solid(mesh: SolidMeshRef<'_>, law: TensorProny, cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.0,
        material_state_identity: [81; 32],
    };
    ViscoElements::new(
        mesh,
        TetViscoField::Uniform(&material),
        &[],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}

fn tet_mesh() -> SolidMeshRef<'static> {
    SolidMeshRef::Tet4 {
        nodes_m: &TET,
        elements: &[[0, 1, 2, 3]],
    }
}

fn hex_mesh() -> SolidMeshRef<'static> {
    SolidMeshRef::Hex8 {
        nodes_m: &HEX,
        elements: &[[0, 1, 2, 3, 4, 5, 6, 7]],
        formulation: Formulation::BBar,
        stiffness_points_per_dir: 2,
    }
}

fn budget() -> ViscoHarmonicMaterialBudget {
    ViscoHarmonicMaterialBudget {
        maximum_points: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_work: 10_000_000,
    }
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-11 * scale.max(1e-20), "{a:e} != {b:e}");
}

fn strain() -> [C64; 6] {
    core::array::from_fn(|i| C64::new((i as f64 - 2.0) * 0.001, (3.0 - i as f64) * 0.0007))
}

fn affine(nodes: &[[f64; 3]], e: [C64; 6]) -> Vec<C64> {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    nodes
        .iter()
        .flat_map(|p| {
            [
                e[0].scale(p[0]) + e[3].scale(s * p[1]) + e[5].scale(s * p[2]),
                e[3].scale(s * p[0]) + e[1].scale(p[1]) + e[4].scale(s * p[2]),
                e[5].scale(s * p[0]) + e[4].scale(s * p[1]) + e[2].scale(p[2]),
            ]
        })
        .collect()
}

fn apply(a: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| a[i][j] * x[j]).sum())
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Independent real-field oracle. The periodic Maxwell spring strain follows
// tau*q_dot + q = tau*epsilon_dot, without any harmonic observation routines.
fn instantaneous(law: &TensorProny, e: &[C64; 6], w: f64, t: f64) -> [f64; 4] {
    let (sin, cos) = (w * t).sin_cos();
    let strain = e.map(|z| z.re * cos + z.im * sin);
    let rate = e.map(|z| w * (-z.re * sin + z.im * cos));
    let mut stress = apply(law.c_eq(), &strain);
    let mut energy = 0.5 * dot(&strain, &stress);
    let viscous = apply(law.kv_pa_s().unwrap(), &rate);
    let kelvin = dot(&rate, &viscous);
    for i in 0..6 {
        stress[i] += viscous[i];
    }
    let mut maxwell = 0.0;
    for branch in law.branches() {
        let x = w * branch.tau_s;
        let (a, b) = (x * x / (1.0 + x * x), x / (1.0 + x * x));
        let q = e.map(|z| (a * z.re + b * z.im) * cos + (a * z.im - b * z.re) * sin);
        let branch_stress = apply(&branch.c_pa, &q);
        let spring = dot(&q, &branch_stress);
        energy += 0.5 * spring;
        maxwell += spring / branch.tau_s;
        for i in 0..6 {
            stress[i] += branch_stress[i];
        }
    }
    [energy, kelvin, maxwell, dot(&stress, &rate)]
}

fn integrate(law: &TensorProny, e: &[C64; 6], w: f64, window: HarmonicWindow) -> [f64; 3] {
    let n = 4096;
    let step = window.duration_s() / n as f64;
    let mut sum = [0.0; 3];
    for k in 0..=n {
        let values = instantaneous(law, e, w, window.start_s() + step * k as f64);
        let weight = if k == 0 || k == n {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for i in 0..3 {
            sum[i] += weight * values[i + 1];
        }
    }
    sum.map(|v| v * step / 3.0)
}

#[test]
fn affine_tetrahedron_recovers_real_interval_heat_work_and_storage() {
    with_cx(|cx, _| {
        let law = law();
        let source = solid(tet_mesh(), law.clone(), cx);
        let w = 7.0;
        let e = strain();
        let displacement = affine(&TET, e);
        let observer = source
            .prepare_harmonic_material_observer(w, budget(), cx)
            .unwrap();
        for window in [
            HarmonicWindow::new(0.137, 0.41).unwrap(),
            HarmonicWindow::new(-0.29, 0.00003).unwrap(),
        ] {
            let field = observer.observe_window(&displacement, window, cx).unwrap();
            let volume = 1.0 / 6.0;
            let start = volume * instantaneous(&law, &e, w, window.start_s())[0];
            let end =
                volume * instantaneous(&law, &e, w, window.start_s() + window.duration_s())[0];
            let [kelvin, maxwell, work] = integrate(&law, &e, w, window).map(|x| x * volume);
            let scale = start.abs().max(end.abs()).max(work.abs());
            assert_eq!(field.window, window);
            assert_eq!(field.points.len(), 1);
            let point = &field.points[0];
            assert_eq!((point.element, point.quadrature_point), (0, 0));
            near(point.volume_m3, volume, volume);
            assert!(point.compatible_strain_mandel.is_none());
            for i in 0..6 {
                near((point.strain_mandel[i] - e[i]).abs(), 0.0, e[i].abs());
            }
            near(field.start_stored_energy_j, start, scale);
            near(field.end_stored_energy_j, end, scale);
            near(field.stored_energy_change_j, end - start, scale);
            near(field.kelvin_heat_j, kelvin, kelvin);
            near(field.maxwell_heat_j, maxwell, maxwell);
            near(field.material_heat_j, kelvin + maxwell, kelvin + maxwell);
            near(field.stress_work_j, work, scale);
            near(field.energy_balance_residual_j, 0.0, scale);
            assert_eq!(
                field,
                observer.observe_window(&displacement, window, cx).unwrap()
            );
        }
        let period = core::f64::consts::TAU / w;
        let whole = observer
            .observe_window(
                &displacement,
                HarmonicWindow::new(0.17, period).unwrap(),
                cx,
            )
            .unwrap();
        let means = observer.observe(&displacement, cx).unwrap();
        near(
            whole.material_heat_j,
            period * means.mean_material_heat_w,
            whole.material_heat_j,
        );
        near(
            whole.stress_work_j,
            period * means.mean_stress_power_w,
            whole.material_heat_j,
        );
        near(
            whole.stored_energy_change_j,
            0.0,
            whole.start_stored_energy_j,
        );
        let left = observer
            .observe_window(&displacement, HarmonicWindow::new(0.17, 0.13).unwrap(), cx)
            .unwrap();
        let right = observer
            .observe_window(
                &displacement,
                HarmonicWindow::new(0.30, period - 0.13).unwrap(),
                cx,
            )
            .unwrap();
        near(
            left.material_heat_j + right.material_heat_j,
            whole.material_heat_j,
            whole.material_heat_j,
        );
        near(
            left.stress_work_j + right.stress_work_j,
            whole.stress_work_j,
            whole.material_heat_j,
        );
        near(
            left.stored_energy_change_j + right.stored_energy_change_j,
            whole.stored_energy_change_j,
            whole.start_stored_energy_j,
        );
    });
}

fn bbar_strain(x: f64, y: f64, amplitude: f64, imaginary: f64) -> ([C64; 6], [C64; 6]) {
    let mut raw = [C64::ZERO; 6];
    raw[0] = C64::new(amplitude * y, imaginary);
    raw[3] = C64::from_re(amplitude * x * core::f64::consts::FRAC_1_SQRT_2);
    let mut projected = raw;
    for z in &mut projected[..3] {
        z.re += amplitude * (0.5 - y) / 3.0;
    }
    (projected, raw)
}

#[test]
fn bbar_keeps_all_original_points_and_integrates_the_constitutive_field() {
    with_cx(|cx, _| {
        let law = law();
        let source = solid(hex_mesh(), law.clone(), cx);
        let w = 4.0;
        let (amplitude, imaginary) = (0.003, 0.001);
        let displacement: Vec<_> = HEX
            .iter()
            .flat_map(|p| {
                [
                    C64::new(amplitude * p[0] * p[1], imaginary * p[0]),
                    C64::ZERO,
                    C64::ZERO,
                ]
            })
            .collect();
        let observer = source
            .prepare_harmonic_material_observer(w, budget(), cx)
            .unwrap();
        let window = HarmonicWindow::new(-0.12, 0.39).unwrap();
        let field = observer.observe_window(&displacement, window, cx).unwrap();
        assert_eq!(field.points.len(), 8);
        let mut start = 0.0;
        let mut end = 0.0;
        let mut integrated = [0.0; 3];
        let g = 1.0 / 3.0_f64.sqrt();
        for x in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
            for y in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
                let (e, _) = bbar_strain(x, y, amplitude, imaginary);
                // Both z points have the same strain, each carrying V/8.
                start += 0.25 * instantaneous(&law, &e, w, window.start_s())[0];
                end += 0.25 * instantaneous(&law, &e, w, window.start_s() + window.duration_s())[0];
                for (sum, value) in integrated.iter_mut().zip(integrate(&law, &e, w, window)) {
                    *sum += 0.25 * value;
                }
            }
        }
        for point in &field.points {
            let raw = point.compatible_strain_mandel.unwrap();
            let e = point.strain_mandel;
            near(point.volume_m3, 0.125, 0.125);
            near(e[0].re + e[1].re + e[2].re, amplitude * 0.5, amplitude);
            near(raw[0].re + 3.0 * e[1].re, amplitude * 0.5, amplitude);
            assert!(point.material.kelvin_heat_j_m3 >= 0.0);
            assert!(point.material.maxwell_heat_j_m3 >= 0.0);
        }
        near(field.start_stored_energy_j, start, start);
        near(field.end_stored_energy_j, end, end);
        near(field.stored_energy_change_j, end - start, start.max(end));
        near(field.kelvin_heat_j, integrated[0], integrated[0]);
        near(field.maxwell_heat_j, integrated[1], integrated[1]);
        near(field.stress_work_j, integrated[2], start.max(end));
        near(field.energy_balance_residual_j, 0.0, start.max(end));
    });
}

#[test]
fn window_admission_preserves_cycle_means_and_rejects_nonphysical_exposure() {
    with_cx(|cx, gate| {
        let source = solid(tet_mesh(), law(), cx);
        let displacement = affine(&TET, strain());
        let window = HarmonicWindow::new(0.1, 0.4).unwrap();
        let small = source
            .prepare_harmonic_material_observer(
                4.0,
                ViscoHarmonicMaterialBudget {
                    maximum_work: 300_000,
                    ..budget()
                },
                cx,
            )
            .unwrap();
        let means = small.observe(&displacement, cx).unwrap();
        assert!(matches!(
            small.observe_window(&displacement, window, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Budget {
                    quantity: "harmonic material window observation work",
                    ..
                }
            ))
        ));
        assert_eq!(means, small.observe(&displacement, cx).unwrap());
        let observer = source
            .prepare_harmonic_material_observer(4.0, budget(), cx)
            .unwrap();
        let good = observer.observe_window(&displacement, window, cx).unwrap();
        let triangular =
            HarmonicWindow::with_profile(0.1, 0.4, HarmonicWindowProfile::Triangular).unwrap();
        assert!(
            observer
                .observe_window(&displacement, triangular, cx)
                .is_err()
        );
        assert!(observer.observe_window(&[], window, cx).is_err());
        assert!(
            observer
                .observe_window(&[C64::from_re(f64::NAN); 12], window, cx)
                .is_err()
        );
        let mut excessive = [C64::ZERO; 6];
        excessive[0] = C64::new(0.019, 0.019);
        assert!(
            observer
                .observe_window(&affine(&TET, excessive), window, cx)
                .is_err()
        );
        assert!(
            source
                .prepare_harmonic_material_observer(101.0, budget(), cx)
                .is_err()
        );
        assert_eq!(
            good,
            observer.observe_window(&displacement, window, cx).unwrap()
        );
        gate.request();
        assert!(observer.observe_window(&displacement, window, cx).is_err());
    });
}

#[test]
fn compatible_strain_limit_cannot_be_hidden_by_bbar_projection() {
    with_cx(|cx, _| {
        let law = law();
        let amplitude = 0.0215;
        let g = 1.0 / 3.0_f64.sqrt();
        let mut raw_exceeds = false;
        for x in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
            for y in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
                let (projected, raw) = bbar_strain(x, y, amplitude, 0.0);
                law.validate_harmonic_mandel_strain(&projected).unwrap();
                raw_exceeds |= law.validate_harmonic_mandel_strain(&raw).is_err();
            }
        }
        assert!(raw_exceeds);
        let source = solid(hex_mesh(), law, cx);
        let observer = source
            .prepare_harmonic_material_observer(4.0, budget(), cx)
            .unwrap();
        let displacement: Vec<_> = HEX
            .iter()
            .flat_map(|p| [C64::from_re(amplitude * p[0] * p[1]), C64::ZERO, C64::ZERO])
            .collect();
        assert!(
            observer
                .observe_window(&displacement, HarmonicWindow::new(0.0, 0.1).unwrap(), cx)
                .is_err()
        );
    });
}
