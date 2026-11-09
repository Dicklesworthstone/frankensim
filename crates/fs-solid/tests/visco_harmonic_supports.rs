//! Original tensor-Prony base excitation with independent mechanical/heat oracles.
//! The declared tensors are manufactured numerical fixtures, not measured data.

use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::visco::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, HarmonicResponseError};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn scope(run: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 84,
                    kernel_id: 96,
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

fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}

fn solid(cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law: TensorProny::new_mandel(
            diagonal([6.0; 6]),
            vec![PronyBranch::new(diagonal([18.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 0.1).unwrap()],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal([1.2, 0.0, 0.0, 0.0, 0.0, 0.0]))
        .unwrap()
        .with_band((0.1, 30.0))
        .unwrap(),
        density_kg_m3: 60.0,
        material_state_identity: [84; 32],
    };
    // Original labels 0 and 3 remain represented. All other coordinates are
    // homogeneous stationary supports. Label 0 will subsequently be driven.
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&material),
        &[1, 2, 4, 5, 6, 7, 8, 9, 10, 11],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 100,
        },
        cx,
    )
    .unwrap()
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 16,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-11,
    }
}

fn observation_budget() -> ViscoHarmonicMaterialBudget {
    ViscoHarmonicMaterialBudget {
        maximum_points: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_work: 10_000_000,
    }
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-10 * scale.max(1e-20), "{a:e} != {b:e}");
}

fn complex_near(a: C64, b: C64) {
    near((a - b).abs(), 0.0, b.abs());
}

// Direct P1 integrals: V=1/6, M=[[1,.5],[.5,1]],
// K_eq=[[2,-1],[-1,1]], K_branch=3*[[1,-1],[-1,1]],
// C_kelvin=.2*[[1,-1],[-1,1]]. No production assembly is used here.
fn pencil(omega: f64) -> [C64; 3] {
    let x = 0.1 * omega;
    let storage = 1.0 + 3.0 * x * x / (1.0 + x * x);
    let loss = 0.2 * omega + 3.0 * x / (1.0 + x * x);
    [
        C64::new(1.0 + storage - omega * omega, -loss),
        C64::new(-storage - 0.5 * omega * omega, loss),
        C64::new(storage - omega * omega, -loss),
    ]
}

#[test]
fn dispersive_base_excitation_retains_cross_mass_reaction_and_original_heat() {
    scope(|cx, _| {
        let source = solid(cx);
        for omega in [0.2, 1.0, 5.0, 17.0] {
            let prepared = source
                .prepare_supported_harmonic_response(&[0], omega, controls(), cx)
                .unwrap();
            assert!(std::ptr::eq(prepared.source(), &source));
            assert_eq!(prepared.represented_dofs(), &[0, 3]);
            assert_eq!(prepared.support_dofs(), &[0]);
            assert_eq!(prepared.free_dofs(), &[3]);
            assert_eq!(prepared.omega_rad_s(), omega);
            let support = C64::new(0.0001, -0.00007);
            let force = [C64::new(0.00003, -0.00002), C64::new(0.0002, 0.00004)];
            let result = prepared.solve(&force, &[support], cx).unwrap();
            let [a00, a01, a11] = pencil(omega);
            let free = (force[1] - a01 * support) / a11;
            let reaction = a00 * support + a01 * free - force[0];
            complex_near(result.displacement_m[0], support);
            complex_near(result.displacement_m[1], free);
            complex_near(result.support_reaction_n[0], reaction);
            let velocity = [support, free].map(|q| C64::new(q.im, -q.re).scale(omega));
            let applied_power =
                0.5 * (force[0].conj() * velocity[0] + force[1].conj() * velocity[1]).re;
            let support_power = 0.5 * (reaction.conj() * velocity[0]).re;
            let loss = -a11.im;
            let heat = 0.5 * omega * loss * (free - support).norm_sq();
            near(result.average_input_power_w, applied_power, heat);
            near(result.average_support_power_w, support_power, heat);
            near(result.average_material_heat_w, heat, heat);
            near(result.power_balance_residual_w, 0.0, heat);
            let observed = source
                .prepare_harmonic_material_observer(omega, observation_budget(), cx)
                .unwrap()
                .observe(&result.displacement_m, cx)
                .unwrap();
            near(observed.mean_material_heat_w, heat, heat);
            let replay = prepared.solve(&force, &[support], cx).unwrap();
            assert_eq!(result.displacement_m, replay.displacement_m);
            assert_eq!(
                result.average_material_heat_w.to_bits(),
                replay.average_material_heat_w.to_bits()
            );
        }
    });
}

#[test]
fn empty_and_all_support_maps_preserve_original_force_and_motion_limits() {
    scope(|cx, _| {
        let source = solid(cx);
        let omega = 2.3;
        let force = [C64::new(0.0003, -0.0001), C64::new(0.0001, 0.0002)];
        let native = source
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap()
            .solve(&force, cx)
            .unwrap();
        let empty = source
            .prepare_supported_harmonic_response(&[], omega, controls(), cx)
            .unwrap();
        let response = empty.solve(&force, &[], cx).unwrap();
        assert!(response.support_reaction_n.is_empty());
        for (&a, &b) in response.displacement_m.iter().zip(&native.displacement_m) {
            complex_near(a, b);
        }
        near(
            response.average_material_heat_w,
            native.average_material_heat_w,
            native.average_material_heat_w,
        );
        let all = source
            .prepare_supported_harmonic_response(&[0, 3], omega, controls(), cx)
            .unwrap();
        assert!(all.free_dofs().is_empty());
        let q = [C64::new(0.0002, -0.0001), C64::new(0.0001, 0.00015)];
        let response = all.solve(&force, &q, cx).unwrap();
        let [a00, a01, a11] = pencil(omega);
        complex_near(
            response.support_reaction_n[0],
            a00 * q[0] + a01 * q[1] - force[0],
        );
        complex_near(
            response.support_reaction_n[1],
            a01 * q[0] + a11 * q[1] - force[1],
        );
        assert_eq!(response.relative_residual, 0.0);
        near(
            response.power_balance_residual_w,
            0.0,
            response.average_material_heat_w,
        );
    });
}

fn real(z: C64, omega: f64, time: f64) -> f64 {
    z.re * (omega * time).cos() + z.im * (omega * time).sin()
}

#[test]
fn close_support_tones_feed_original_coherent_heat_on_the_physical_clock() {
    scope(|cx, _| {
        let source = solid(cx);
        let frequencies = [5.0, 5.03];
        let supports = [C64::new(0.001, 0.0002), C64::new(0.0008, 0.0001)];
        let responses: Vec<_> = frequencies
            .iter()
            .zip(supports)
            .map(|(&omega, support)| {
                source
                    .prepare_supported_harmonic_response(&[0], omega, controls(), cx)
                    .unwrap()
                    .solve(&[C64::ZERO; 2], &[support], cx)
                    .unwrap()
            })
            .collect();
        let displacements: Vec<_> = responses
            .iter()
            .map(|r| r.displacement_m.as_slice())
            .collect();
        let start = 0.137;
        let dt = 0.31;
        let observer = source
            .prepare_harmonic_material_spectrum_observer(&frequencies, observation_budget(), cx)
            .unwrap();
        let observed = observer
            .observe_window(&displacements, HarmonicWindow::new(start, dt).unwrap(), cx)
            .unwrap();
        let fields = |time| {
            let (mut e, mut edot, mut spring, mut u0, mut v0) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for (i, &omega) in frequencies.iter().enumerate() {
                let strain = responses[i].displacement_m[1] - responses[i].displacement_m[0];
                let x = 0.1 * omega;
                let branch = strain * C64::new(x * x / (1.0 + x * x), -x / (1.0 + x * x));
                e += real(strain, omega, time);
                edot += real(C64::new(strain.im, -strain.re).scale(omega), omega, time);
                spring += real(branch, omega, time);
                u0 += real(supports[i], omega, time);
                v0 += real(
                    C64::new(supports[i].im, -supports[i].re).scale(omega),
                    omega,
                    time,
                );
            }
            let stored = (3.0 * (e * e + u0 * u0) + 9.0 * spring * spring) / 6.0;
            let kelvin = 1.2 * edot * edot / 6.0;
            let maxwell = 180.0 * spring * spring / 6.0;
            let work = ((6.0 * e + 18.0 * spring + 1.2 * edot) * edot + 6.0 * u0 * v0) / 6.0;
            [stored, kelvin, maxwell, work]
        };
        let steps = 4096;
        let mut integrals = [0.0; 3];
        for j in 0..=steps {
            let weight = if j == 0 || j == steps {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let values = fields(start + dt * j as f64 / steps as f64);
            for k in 0..3 {
                integrals[k] += weight * values[k + 1];
            }
        }
        for value in &mut integrals {
            *value *= dt / (3.0 * steps as f64);
        }
        let scale = integrals.iter().map(|v| v.abs()).sum::<f64>();
        near(observed.start_stored_energy_j, fields(start)[0], scale);
        near(observed.end_stored_energy_j, fields(start + dt)[0], scale);
        near(observed.kelvin_heat_j, integrals[0], scale);
        near(observed.maxwell_heat_j, integrals[1], scale);
        near(observed.stress_work_j, integrals[2], scale);
        near(observed.energy_balance_residual_j, 0.0, scale);
        let independent_heat: f64 = frequencies
            .iter()
            .enumerate()
            .map(|(i, &omega)| {
                source
                    .prepare_harmonic_material_observer(omega, observation_budget(), cx)
                    .unwrap()
                    .observe_window(
                        &responses[i].displacement_m,
                        HarmonicWindow::new(start, dt).unwrap(),
                        cx,
                    )
                    .unwrap()
                    .material_heat_j
            })
            .sum();
        assert!(observed.material_heat_j > 1.8 * independent_heat);
    });
}

#[test]
fn source_support_budget_and_strain_refusals_leave_the_factor_reusable() {
    scope(|cx, gate| {
        let source = solid(cx);
        for supports in [&[3, 0][..], &[0, 0], &[1]] {
            assert!(
                source
                    .prepare_supported_harmonic_response(supports, 1.0, controls(), cx)
                    .is_err()
            );
        }
        for omega in [0.0, -1.0, f64::NAN, 0.05, 40.0] {
            assert!(
                source
                    .prepare_supported_harmonic_response(&[0], omega, controls(), cx)
                    .is_err()
            );
        }
        let bytes = 344 * 4 + 544 * 2;
        let budget = HarmonicResponseControls {
            maximum_payload_bytes: bytes - 1,
            ..controls()
        };
        assert!(matches!(
            source.prepare_supported_harmonic_response(&[0], 1.0, budget, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        let prepared = source
            .prepare_supported_harmonic_response(
                &[0, 3],
                1.0,
                HarmonicResponseControls {
                    maximum_payload_bytes: bytes,
                    ..controls()
                },
                cx,
            )
            .unwrap();
        let force = [C64::ZERO; 2];
        assert!(
            prepared
                .solve(&force, &[C64::new(0.02, 0.02), C64::ZERO], cx)
                .is_err()
        );
        assert!(
            prepared
                .solve(&force, &[C64::new(f64::NAN, 0.0), C64::ZERO], cx)
                .is_err()
        );
        assert!(prepared.solve(&force, &[], cx).is_err());
        let small = [C64::new(0.0001, 0.0001), C64::ZERO];
        assert!(prepared.solve(&force, &small, cx).is_ok());
        gate.request();
        assert!(matches!(
            prepared.solve(&force, &small, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Cancelled
            ))
        ));
    });
}
