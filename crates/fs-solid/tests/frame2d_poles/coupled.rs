//! G1 original-equation poles and work, G3 port-coordinate covariance,
//! G0 unsupported-law/budget refusal, G4 cancellation. Shares native fixtures.

use super::{BAND, NO_RAYLEIGH, RADIUS, guided, with_cx};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::harmonic::{BeynOptions, PoleOptions, ResonanceRegion};
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, QuadraticStorage, Storage};
use fs_solid::frame2d::{
    FrameMassFormulation,
    coupled::{FrameCoupledControls, FrameCoupledError, FrameDynamicPortLoad, FramePortDynamics},
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

fn controls() -> FrameCoupledControls {
    FrameCoupledControls {
        maximum_coordinates: 12,
        maximum_ports: 4,
        maximum_states: 16,
        maximum_payload_bytes: 2 << 20,
        maximum_work: 1 << 36,
        maximum_relative_residual: 1e-10,
    }
}
fn options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 16.0,
        step_tolerance: 1e-13,
        residual_tolerance: 1e-13,
        max_iterations: 40,
        max_dimension: 16,
        max_dense_bytes: 2 << 20,
        max_work: 100_000_000,
    }
}
fn contour(n: usize) -> BeynOptions {
    BeynOptions {
        nodes: 128,
        probe_columns: n,
        rank_relative_tolerance: 1e-10,
        count_tolerance: 1e-6,
        quadrature_tolerance: 1e-6,
        contour_min_reciprocal_condition: 1e-14,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 120,
        max_dimension: 16,
        max_dense_bytes: 2 << 20,
        max_work: 2_000_000_000,
        refinement: options(),
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "actual {a:?}, expected {b:?}");
}
fn scalar(mass: f64, stiffness: f64) -> Frame2d {
    let mut f = Frame2d::new(Geometry::Linear);
    f.add_node(0.0, 0.0);
    f.add_node(1.0, 0.0);
    f.fix(0, [true; 3]);
    f.fix(1, [false, true, true]);
    f.add_mass(1, mass);
    f.add_element(0, 1, 5, 0.0, &|| SectionModel::Elastic {
        ea: stiffness,
        ei: 1.0,
    })
    .unwrap();
    f
}
const TIP: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn relaxation_and_hidden_state() -> PortHamiltonian {
    // x_dot = v - 1.5*x, opposing force .25*x. The second source coordinate
    // is a real independent internal relaxation, with no external port output.
    PortHamiltonian::new(
        2,
        1,
        vec![0.0; 4],
        vec![6.0, 0.0, 0.0, 5.0],
        vec![1.0, 0.0],
        Box::new(QuadraticStorage::new(vec![0.25, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap()
}

#[test]
fn dynamic_memory_moves_poles_and_preserves_the_hidden_relaxation() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = scalar(1.0, 4.0);
        let accepted = format!("{frame:?}");
        let load = relaxation_and_hidden_state();
        let model = frame
            .prepare_coupled_resonance(
                Rayleigh { a0: 4.5, a1: 0.0 },
                FrameDynamicPortLoad {
                    force_map: &TIP,
                    dynamics: FramePortDynamics::Ode(&load),
                },
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(model.free_dofs(), [3]);
        assert_eq!(model.load_state_range(), 2..4);
        assert_eq!(model.load_algebraic_range(), 4..4);
        // Eliminating the first source memory independently gives
        // (s+1.5)*(s*s+4.5*s+4)+.25*s = (s+1)*(s+2)*(s+3).
        // The isolated source pole -1.5 is NOT the coupled mode. A zero-output
        // source state still contributes its independent physical pole -5.
        let spectrum = model
            .resonance_spectrum(
                ResonanceRegion {
                    center: C64::from_re(-3.0),
                    real_radius: 2.75,
                    imaginary_radius: 1.5,
                },
                &contour(4),
                cx,
            )
            .unwrap();
        assert_eq!(
            (
                spectrum.argument_count,
                spectrum.resolved_rank,
                spectrum.poles.len()
            ),
            (4, 4, 4)
        );
        for wanted in [-1.0, -2.0, -3.0, -5.0] {
            let p = spectrum
                .poles
                .iter()
                .min_by(|a, b| {
                    (a.laplace_pole_per_s - C64::from_re(wanted))
                        .abs()
                        .total_cmp(&(b.laplace_pole_per_s - C64::from_re(wanted)).abs())
                })
                .unwrap();
            close(p.laplace_pole_per_s, C64::from_re(wanted), 2e-9);
            assert!(p.relative_residual <= controls().maximum_relative_residual);
            assert!(p.decay_power_residual_w.abs() < 1e-8 * p.dissipated_power_w);
            if wanted == -5.0 {
                assert!(p.port_output[0].abs() < 1e-10);
                assert!(p.stored_energy_j > 0.0);
            } else {
                close(
                    p.state[2],
                    p.state[0] * p.laplace_pole_per_s / (p.laplace_pole_per_s + C64::from_re(1.5)),
                    1e-10,
                );
            }
        }
        // A numerical dark-mode vector may retain a roundoff-sized frame block.
        // Both the shared pole observer and original q/v check must scale the
        // homogeneous row coefficients by the complete modal amplitude.
        let dark = model
            .resonance_pole(
                C64::from_re(-5.0),
                &[
                    C64::from_re(1e-13),
                    C64::from_re(-1e-13),
                    C64::from_re(1e-13),
                    C64::ONE,
                ],
                &options(),
                cx,
            )
            .unwrap();
        assert!(dark.state[0].abs() > 0.0 && dark.port_output[0].abs() < 1e-12);
        assert!(dark.relative_residual < controls().maximum_relative_residual);
        // Actual AC state, storage and heat follow the SAME complete operator.
        for omega in [0.5, 2.0, 8.0] {
            let s = C64::new(0.0, -omega);
            let force = C64::new(0.7, -0.2);
            let displacement = force
                / (s * s
                    + s.scale(4.5)
                    + C64::from_re(4.0)
                    + s.scale(0.25) / (s + C64::from_re(1.5)));
            let velocity = s * displacement;
            let memory = velocity / (s + C64::from_re(1.5));
            let response = model
                .factor(omega, cx)
                .unwrap()
                .response(&[force], cx)
                .unwrap();
            close(response.state[0], displacement, 1e-12);
            close(response.output[0], velocity, 1e-12);
            close(response.state[2], memory, 1e-12);
            assert_eq!(response.state[3], C64::ZERO);
            let energy =
                (4.0 * displacement.norm_sq() + velocity.norm_sq() + 0.25 * memory.norm_sq()) / 4.0;
            let heat = (4.5 * velocity.norm_sq() + 0.375 * memory.norm_sq()) / 2.0;
            assert!((response.mean_stored_energy_j - energy).abs() < 1e-12);
            assert!((response.dissipated_power_w - heat).abs() < 1e-12);
            assert!(response.power_residual_w.abs() < 1e-12);
        }
        assert_eq!(format!("{frame:?}"), accepted);
    });
}

fn inertance_and_drag(mass: f64, drag: f64) -> DescriptorPortHamiltonian {
    // Source state [z, lambda, mu]. z_dot=lambda, 0=-z/m+v,
    // 0=-mu/drag+v, opposing effort=lambda+mu. Only z stores energy.
    DescriptorPortHamiltonian::new(
        3,
        1,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / drag],
        vec![0.0, 1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![1.0 / mass], 1).unwrap()),
    )
    .unwrap()
}

#[test]
fn descriptor_inertance_and_drag_keep_original_constraints_and_physical_work() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = scalar(2.0, 8.0);
        let load = inertance_and_drag(3.0, 4.0);
        let model = frame
            .prepare_coupled_resonance(
                Rayleigh { a0: 0.5, a1: 0.0 },
                FrameDynamicPortLoad {
                    force_map: &TIP,
                    dynamics: FramePortDynamics::Descriptor(&load),
                },
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(model.dimension(), 5);
        assert_eq!(model.load_algebraic_range(), 3..5);
        let expected = C64::new(-0.5, -(1.35_f64).sqrt());
        let seed = [
            C64::ONE,
            expected.scale(2.0),
            expected.scale(3.0),
            (expected * expected).scale(3.0),
            expected.scale(4.0),
        ];
        let pole = model
            .resonance_pole(expected + C64::new(0.01, 0.02), &seed, &options(), cx)
            .unwrap();
        close(pole.laplace_pole_per_s, expected, 1e-11);
        let v = pole.effort[1];
        close(pole.state[2], v.scale(3.0), 1e-11);
        close(pole.state[3], expected * v.scale(3.0), 1e-11);
        close(pole.state[4], v.scale(4.0), 1e-11);
        let energy = (8.0 * pole.state[0].norm_sq() + 5.0 * v.norm_sq()) / 2.0;
        let loss = 5.0 * v.norm_sq();
        assert!((pole.stored_energy_j - energy).abs() < 1e-12);
        assert!((pole.dissipated_power_w - loss).abs() < 1e-12);
        assert!(pole.decay_power_residual_w.abs() < 1e-10 * loss);
        for omega in [0.5, 2.0, 7.0] {
            let s = C64::new(0.0, -omega);
            let expected_q = C64::ONE / ((s * s).scale(5.0) + s.scale(5.0) + C64::from_re(8.0));
            let response = model
                .factor(omega, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], expected_q, 1e-12);
            close(response.state[2], (s * expected_q).scale(3.0), 1e-12);
            assert!(response.power_residual_w.abs() < 1e-12);
        }
    });
}

#[test]
fn energy_null_rigid_position_is_retained_with_a_descriptor_resistance() {
    with_cx(&CancelGate::new(), |cx| {
        let mut frame = Frame2d::new(Geometry::Linear);
        frame.add_node(0.0, 0.0);
        frame.add_mass(0, 2.0);
        frame.fix(0, [false, true, true]);
        let load = DescriptorPortHamiltonian::new(
            1,
            0,
            1,
            vec![0.0],
            vec![0.5],
            vec![1.0],
            Box::new(QuadraticStorage::new(vec![], 0).unwrap()),
        )
        .unwrap();
        let model = frame
            .prepare_coupled_resonance(
                NO_RAYLEIGH,
                FrameDynamicPortLoad {
                    force_map: &[1.0, 0.0, 0.0],
                    dynamics: FramePortDynamics::Descriptor(&load),
                },
                controls(),
                cx,
            )
            .unwrap();
        let zero = model
            .resonance_pole(C64::ZERO, &[C64::ONE, C64::ZERO, C64::ZERO], &options(), cx)
            .unwrap();
        close(zero.laplace_pole_per_s, C64::ZERO, 1e-14);
        assert_eq!(zero.stored_energy_j, 0.0);
        assert_eq!(zero.dissipated_power_w, 0.0);
        let drag = model
            .resonance_pole(
                C64::from_re(-1.0),
                &[C64::ONE, C64::from_re(-2.0), C64::from_re(-2.0)],
                &options(),
                cx,
            )
            .unwrap();
        close(drag.laplace_pole_per_s, C64::from_re(-1.0), 1e-12);
        assert_eq!(model.load_algebraic_range(), 2..3);
    });
}

fn two_coordinates() -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    for x in [0.0, 1.0, 2.0] {
        frame.add_node(x, 0.0);
    }
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame.fix(2, [false, true, true]);
    for (i, j) in [(0, 1), (1, 2)] {
        frame
            .add_element(i, j, 5, 6.0, &|| SectionModel::Elastic { ea: 8.0, ei: 1.0 })
            .unwrap();
    }
    frame
}

#[test]
fn full_native_mutual_mass_and_multiport_storage_are_covariant() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = two_coordinates();
        let b = [[1.0, -0.5], [0.25, 2.0]];
        let gl = [[1.0, 0.25], [0.5, 2.0]];
        let ql = [[3.0, 1.0], [1.0, 2.0]];
        let m = [[4.0, 1.0], [1.0, 2.0]];
        let mut k = [[16.0, -8.0], [-8.0, 8.0]];
        let mut t = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for p in 0..2 {
                    t[i][j] += gl[i][p] * b[j][p];
                }
            }
        }
        for i in 0..2 {
            for j in 0..2 {
                for a in 0..2 {
                    for z in 0..2 {
                        k[i][j] += t[a][i] * ql[a][z] * t[z][j];
                    }
                }
            }
        }
        let omega = 1.7;
        let s = C64::new(0.0, -omega);
        let mut d = [[C64::ZERO; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                d[i][j] = (s * s + s.scale(0.1)).scale(m[i][j]) + C64::from_re(k[i][j]);
            }
        }
        let force = [C64::new(1.0, 0.3), C64::new(-0.2, 0.4)];
        let det = d[0][0] * d[1][1] - d[0][1] * d[1][0];
        let expected = [
            (force[0] * d[1][1] - force[1] * d[0][1]) / det,
            (d[0][0] * force[1] - d[1][0] * force[0]) / det,
        ];
        for port_scales in [[1.0, 1.0], [1e4, -1e-3]] {
            let mut map = vec![0.0; 18];
            // A supported force row changes reactions only; no motion enters it.
            map[0] = 1000.0;
            for i in 0..2 {
                for p in 0..2 {
                    map[(3 + 3 * i) * 2 + p] = b[i][p] * port_scales[p];
                }
            }
            let g = vec![
                gl[0][0] / port_scales[0],
                gl[0][1] / port_scales[1],
                gl[1][0] / port_scales[0],
                gl[1][1] / port_scales[1],
            ];
            let load = PortHamiltonian::new(
                2,
                2,
                vec![0.0; 4],
                vec![0.0; 4],
                g,
                Box::new(QuadraticStorage::new(vec![3.0, 1.0, 1.0, 2.0], 2).unwrap()),
            )
            .unwrap();
            let model = frame
                .prepare_coupled_resonance(
                    Rayleigh { a0: 0.1, a1: 0.0 },
                    FrameDynamicPortLoad {
                        force_map: &map,
                        dynamics: FramePortDynamics::Ode(&load),
                    },
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(model.free_dofs(), [3, 6]);
            let response = model
                .factor(omega, cx)
                .unwrap()
                .response(&force, cx)
                .unwrap();
            for i in 0..2 {
                close(response.state[i], expected[i], 2e-11);
                close(response.output[i], s * expected[i], 2e-11);
                close(
                    response.state[2 + i],
                    (s * expected[0]).scale(m[i][0]) + (s * expected[1]).scale(m[i][1]),
                    2e-11,
                );
                close(
                    response.state[4 + i],
                    expected[0].scale(t[i][0]) + expected[1].scale(t[i][1]),
                    2e-11,
                );
            }
            assert!(response.power_residual_w.abs() < 1e-11);
            // Independent 2x2 generalized characteristic equation for oscillatory modes.
            let det_m = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            let det_k = k[0][0] * k[1][1] - k[0][1] * k[1][0];
            let cross = k[0][0] * m[1][1] + k[1][1] * m[0][0] - 2.0 * k[0][1] * m[0][1];
            for lambda in [
                (cross - (cross * cross - 4.0 * det_m * det_k).sqrt()) / (2.0 * det_m),
                (cross + (cross * cross - 4.0 * det_m * det_k).sqrt()) / (2.0 * det_m),
            ] {
                let rate = C64::new(-0.05, -(lambda - 0.0025).sqrt());
                let x = [-(k[0][1] - lambda * m[0][1]), k[0][0] - lambda * m[0][0]];
                let seed = [
                    C64::from_re(x[0]),
                    C64::from_re(x[1]),
                    rate.scale(m[0][0] * x[0] + m[0][1] * x[1]),
                    rate.scale(m[1][0] * x[0] + m[1][1] * x[1]),
                    C64::from_re(t[0][0] * x[0] + t[0][1] * x[1]),
                    C64::from_re(t[1][0] * x[0] + t[1][1] * x[1]),
                ];
                let pole = model
                    .resonance_pole(rate + C64::new(0.001, 0.002), &seed, &options(), cx)
                    .unwrap();
                close(pole.laplace_pole_per_s, rate, 2e-9);
            }
        }
    });
}

#[test]
fn source_material_and_authored_loss_keep_one_owner_and_ac_band_gate() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = guided(2.0);
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let load = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![0.0],
            vec![1.0],
            Box::new(QuadraticStorage::new(vec![2.0 * area], 1).unwrap()),
        )
        .unwrap();
        let model = frame
            .prepare_coupled_resonance(
                Rayleigh { a0: 0.5, a1: 0.0 },
                FrameDynamicPortLoad {
                    force_map: &TIP,
                    dynamics: FramePortDynamics::Ode(&load),
                },
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(model.material_omega_band_rad_s(), Some(BAND));
        // M=2A, total K=10A, material C=2A plus authored C=A.
        let rate = C64::new(-0.75, -(5.0_f64 - 0.75 * 0.75).sqrt());
        let pole = model
            .resonance_pole(
                rate + C64::new(0.001, 0.002),
                &[C64::ONE, rate.scale(2.0 * area), C64::ONE],
                &options(),
                cx,
            )
            .unwrap();
        close(pole.laplace_pole_per_s, rate, 2e-10);
        let velocity = pole.effort[1];
        assert!((pole.dissipated_power_w - 3.0 * area * velocity.norm_sq()).abs() < 1e-11);
        assert!(model.factor(1.0, cx).is_ok());
        for omega in [0.25, 4.0] {
            assert!(matches!(
                model.factor(omega, cx),
                Err(FrameCoupledError::Invalid(_))
            ));
        }
    });
}

#[test]
fn resource_bad_map_nonlinear_law_and_cancellation_refuse_without_state_changes() {
    let frame = scalar(1.0, 4.0);
    let accepted = format!("{frame:?}");
    let load = relaxation_and_hidden_state();
    with_cx(&CancelGate::new(), |cx| {
        for cap in [
            FrameCoupledControls {
                maximum_coordinates: 5,
                ..controls()
            },
            FrameCoupledControls {
                maximum_ports: 0,
                ..controls()
            },
            FrameCoupledControls {
                maximum_states: 3,
                ..controls()
            },
            FrameCoupledControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
            FrameCoupledControls {
                maximum_work: 1,
                ..controls()
            },
        ] {
            assert!(matches!(
                frame.prepare_coupled_resonance(
                    NO_RAYLEIGH,
                    FrameDynamicPortLoad {
                        force_map: &TIP,
                        dynamics: FramePortDynamics::Ode(&load)
                    },
                    cap,
                    cx
                ),
                Err(FrameCoupledError::Budget { .. })
            ));
        }
        assert!(matches!(
            frame.prepare_coupled_resonance(
                NO_RAYLEIGH,
                FrameDynamicPortLoad {
                    force_map: &[1.0],
                    dynamics: FramePortDynamics::Ode(&load)
                },
                controls(),
                cx
            ),
            Err(FrameCoupledError::Invalid(_))
        ));
        struct Undeclared;
        impl Storage for Undeclared {
            fn hamiltonian(&self, x: &[f64]) -> f64 {
                x[0].powi(4) / 4.0
            }
            fn gradient(&self, x: &[f64], out: &mut [f64]) {
                out[0] = x[0].powi(3);
            }
        }
        let nonlinear =
            PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(Undeclared))
                .unwrap();
        assert!(matches!(
            frame.prepare_coupled_resonance(
                NO_RAYLEIGH,
                FrameDynamicPortLoad {
                    force_map: &TIP,
                    dynamics: FramePortDynamics::Ode(&nonlinear)
                },
                controls(),
                cx
            ),
            Err(FrameCoupledError::Invalid(_))
        ));
    });
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            frame.prepare_coupled_resonance(
                NO_RAYLEIGH,
                FrameDynamicPortLoad {
                    force_map: &TIP,
                    dynamics: FramePortDynamics::Ode(&load)
                },
                controls(),
                cx
            ),
            Err(FrameCoupledError::Cancelled)
        ));
    });
    assert_eq!(format!("{frame:?}"), accepted);
}
