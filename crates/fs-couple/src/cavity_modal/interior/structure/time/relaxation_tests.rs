//! Independent original-coordinate trajectories and transactional runtime tests.

mod numerical {
    use super::super::*;
    use fs_math::{c64::C64, det};
    use fs_phs::RelaxationBranch;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
            "actual {actual:.16e}, expected {expected:.16e}"
        );
    }
    fn exponential(s: C64, t: f64) -> C64 {
        C64::new(det::cos(s.im * t), det::sin(s.im * t)).scale(det::exp(s.re * t))
    }
    fn roots() -> [C64; 3] {
        [
            C64::from_re(-0.5),
            C64::new(-0.25, det::sqrt(31.0) / 4.0),
            C64::new(-0.25, -det::sqrt(31.0) / 4.0),
        ]
    }
    fn arm() -> RelaxationBranch {
        RelaxationBranch {
            projection: vec![1.0, 0.0, 0.0, 0.0, 0.0],
            stiffness: 1.25,
            relaxation_time_s: 1.0,
        }
    }
    struct Tetrahedron {
        mass: Vec<f64>,
        damping: Vec<f64>,
        stiffness: Vec<f64>,
        coupling: Vec<f64>,
    }
    impl Tetrahedron {
        fn new() -> Self {
            // Exact P1 matrices for unit-coordinate tetrahedron, V=1/6,
            // rho=bulk=.2. Uniform pressure is an invariant of affine dilation.
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
            for a in 0..4 {
                damping[a + 1] = -0.125;
                damping[(a + 1) * 5] = 0.125;
                for b in 0..4 {
                    mass[(a + 1) * 5 + b + 1] = (if a == b { 2.0 } else { 1.0 }) / 24.0;
                    stiffness[(a + 1) * 5 + b + 1] = (0..3)
                        .map(|j| gradients[a][j] * gradients[b][j])
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
        fn source(&self) -> kernel::Source<'_> {
            kernel::Source {
                structural: 1,
                coordinates: 5,
                mass: &self.mass,
                damping: &self.damping,
                stiffness: &self.stiffness,
                coupling: &self.coupling,
            }
        }
    }
    fn squared_envelope_integral(amplitude: C64, s: C64, duration: f64) -> f64 {
        // Integral Re(A exp(s t))^2 from direct trigonometric integration.
        let scalar = (det::exp(2.0 * s.re * duration) - 1.0) / (2.0 * s.re);
        let oscillating = (exponential(s.scale(2.0), duration) - C64::ONE) / s.scale(2.0);
        0.5 * amplitude.norm_sq() * scalar
            + 0.5 * (amplitude.re * amplitude.re - amplitude.im * amplitude.im) * oscillating.re
            - amplitude.re * amplitude.im * oscillating.im
    }

    #[test]
    fn all_three_exact_cubic_modes_evolve_original_fields_and_independent_heat() {
        let source = Tetrahedron::new();
        let duration = 0.27;
        let prepared = kernel::Prepared::new(
            source.source(),
            vec![arm()],
            duration,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        for s in roots() {
            // P(s)=(s+.5)(s^2+.5s+2); no production eigensolver in this oracle.
            let psi = C64::from_re(-0.6) / s;
            let sigma = s.scale(1.25) / (C64::ONE + s);
            let mut q = vec![psi.re; 5];
            q[0] = 1.0;
            let mut v = vec![-0.6; 5];
            v[0] = s.re;
            let mut state = prepared
                .state(0.0, &q, &v, &[sigma.re], &mut || false)
                .unwrap();
            let mut workspace = prepared.workspace();
            for step in 1..=3 {
                let t = step as f64 * duration;
                let frame = prepared
                    .step_into(&mut state, &[0.0], &[0.0; 4], &mut workspace, &mut || false)
                    .unwrap();
                let e = exponential(s, t);
                close(state.canonical[0], e.re, 7e-11);
                close(state.velocity[0], (s * e).re, 7e-11);
                for a in 1..5 {
                    close(state.canonical[a], (psi * e).re, 7e-11);
                    close(state.velocity[a], -0.6 * e.re, 7e-11);
                }
                close(state.stress[0], (sigma * e).re, 7e-11);
                close(
                    state.viscous_strain[0],
                    state.canonical[0] - state.stress[0] / 1.25,
                    2e-12,
                );
                let loss =
                    squared_envelope_integral(sigma * exponential(s, t - duration), s, duration)
                        / 1.25;
                close(frame.memory[0].dissipation_j, loss, 2e-11);
                close(frame.work.relaxation_dissipation_j, loss, 2e-11);
                close(frame.fluid_force[0], -0.3 * e.re, 2e-11);
                close(
                    frame.work.final_energy.spatial.fluid_compression_j,
                    0.15 * e.re * e.re,
                    2e-11,
                );
                close(frame.work.balance_residual_j(), 0.0, 3e-11);
                assert!(frame.work.balance_tolerance_j < 1e-8);
            }
        }
    }

    #[test]
    fn held_force_matches_independent_cubic_partial_fractions_and_work() {
        let source = Tetrahedron::new();
        let dt = 0.63;
        let force = 0.8;
        let prepared = kernel::Prepared::new(
            source.source(),
            vec![arm()],
            dt,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        let mut state = prepared
            .state(0.0, &[0.0; 5], &[0.0; 5], &[0.0], &mut || false)
            .unwrap();
        let mut workspace = prepared.workspace();
        let frame = prepared
            .step_into(&mut state, &[force], &[0.0; 4], &mut workspace, &mut || {
                false
            })
            .unwrap();
        let (mut q, mut v, mut integral, mut stress) =
            (C64::ONE, C64::ZERO, C64::from_re(dt), C64::ZERO);
        let mut stress_amplitudes = [C64::ZERO; 3];
        for (i, s) in roots().into_iter().enumerate() {
            // Original transfer (1+s)/(s*P(s)), P=s^3+s^2+2.25s+1.
            let derivative = s * s.scale(3.0) + s.scale(2.0) + C64::from_re(2.25);
            let residue = (C64::ONE + s) / (s * derivative);
            let e = exponential(s, dt);
            q = q + residue * e;
            v = v + residue * s * e;
            integral = integral + residue * (e - C64::ONE) / s;
            stress_amplitudes[i] = C64::from_re(1.25 * force) / derivative;
            stress = stress + stress_amplitudes[i] * e;
        }
        close(state.canonical[0], force * q.re, 6e-11);
        close(state.velocity[0], force * v.re, 6e-11);
        close(state.stress[0], stress.re, 6e-11);
        for a in 1..5 {
            close(state.velocity[a], -0.6 * force * q.re, 6e-11);
            close(state.canonical[a], -0.6 * force * integral.re, 6e-11);
        }
        close(
            frame.work.structural_input_work_j,
            force * force * q.re,
            3e-11,
        );
        let mut heat = C64::ZERO;
        for (i, s) in roots().into_iter().enumerate() {
            for (j, t) in roots().into_iter().enumerate() {
                heat = heat
                    + stress_amplitudes[i]
                        * stress_amplitudes[j]
                        * (exponential(s + t, dt) - C64::ONE)
                        / (s + t);
            }
        }
        close(frame.memory[0].dissipation_j, heat.re / 1.25, 3e-11);
        assert!(frame.work.fluid_to_structure_work_j < 0.0);
    }

    #[test]
    fn dark_history_relaxes_without_displacement_and_retains_its_heat() {
        let source = Tetrahedron::new();
        let dt = 0.4;
        let branch = RelaxationBranch {
            projection: vec![0.0; 5],
            stiffness: 2.0,
            relaxation_time_s: 0.7,
        };
        let prepared = kernel::Prepared::new(
            source.source(),
            vec![branch],
            dt,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        let mut state = prepared
            .state(1.0, &[0.0; 5], &[0.0; 5], &[0.9], &mut || false)
            .unwrap();
        let mut workspace = prepared.workspace();
        let frame = prepared
            .step_into(&mut state, &[0.0], &[0.0; 4], &mut workspace, &mut || false)
            .unwrap();
        let ratio = det::exp(-dt / 0.7);
        close(state.stress[0], 0.9 * ratio, 2e-12);
        close(
            frame.memory[0].dissipation_j,
            0.9 * 0.9 / 4.0 * (1.0 - ratio * ratio),
            2e-12,
        );
        assert!(state.velocity.iter().all(|v| *v == 0.0));
        assert!(state.canonical[..5].iter().all(|v| *v == 0.0));
        assert_eq!(frame.memory[0].structural_work_j, 0.0);
        assert_eq!(frame.memory[0].acoustic_work_j, 0.0);
    }

    #[test]
    fn mixed_and_wall_histories_match_independent_original_equation_rk4() {
        let mass = [2.0, 0.0, 0.0, 3.0];
        let damping = [0.3, -0.7, 0.7, 0.2];
        let stiffness = [4.0, 0.0, 0.0, 5.0];
        let coupling = [0.7];
        let source = kernel::Source {
            structural: 1,
            coordinates: 2,
            mass: &mass,
            damping: &damping,
            stiffness: &stiffness,
            coupling: &coupling,
        };
        let maps = [[0.0, 1.0], [0.5, -0.4]];
        let ks = [1.1, 0.8];
        let taus = [0.3, 1.2];
        let branches = (0..2)
            .map(|j| RelaxationBranch {
                projection: maps[j].to_vec(),
                stiffness: ks[j],
                relaxation_time_s: taus[j],
            })
            .collect();
        let dt = 0.61;
        let input = [0.35, -0.18];
        let prepared = kernel::Prepared::new(
            source,
            branches,
            dt,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        let mut state = prepared
            .state(0.0, &[0.2, -0.1], &[0.3, -0.4], &[0.15, -0.2], &mut || {
                false
            })
            .unwrap();
        let mut workspace = prepared.workspace();
        let frame = prepared
            .step_into(
                &mut state,
                &input[..1],
                &input[1..],
                &mut workspace,
                &mut || false,
            )
            .unwrap();
        // Independent physical q,v,sigma ODE, with direct loss/work quadratures.
        // No canonical pHS, mass factor, energy transform or exponential here.
        let derivative = |y: &[f64; 12]| {
            let mut d = [0.0; 12];
            d[0] = y[2];
            d[1] = y[3];
            for i in 0..2 {
                d[2 + i] = (input[i]
                    - stiffness[i * 2 + i] * y[i]
                    - damping[i * 2] * y[2]
                    - damping[i * 2 + 1] * y[3]
                    - maps[0][i] * y[4]
                    - maps[1][i] * y[5])
                    / mass[i * 2 + i];
            }
            for j in 0..2 {
                d[4 + j] = ks[j] * (maps[j][0] * y[2] + maps[j][1] * y[3]) - y[4 + j] / taus[j];
                d[6 + 3 * j] = y[4 + j] * y[4 + j] / (ks[j] * taus[j]);
                d[7 + 3 * j] = y[4 + j] * maps[j][0] * y[2];
                d[8 + 3 * j] = y[4 + j] * maps[j][1] * y[3];
            }
            d
        };
        let mut y = [0.0; 12];
        y[..6].copy_from_slice(&[0.2, -0.1, 0.3, -0.4, 0.15, -0.2]);
        let h = dt / 2048.0;
        for _ in 0..2048 {
            let a = derivative(&y);
            let mut z = [0.0; 12];
            for i in 0..12 {
                z[i] = y[i] + 0.5 * h * a[i];
            }
            let b = derivative(&z);
            for i in 0..12 {
                z[i] = y[i] + 0.5 * h * b[i];
            }
            let c = derivative(&z);
            for i in 0..12 {
                z[i] = y[i] + h * c[i];
            }
            let d = derivative(&z);
            for i in 0..12 {
                y[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
            }
        }
        for i in 0..2 {
            close(state.canonical[i], y[i], 2e-11);
            close(state.velocity[i], y[2 + i], 2e-11);
        }
        for j in 0..2 {
            close(state.stress[j], y[4 + j], 2e-11);
            close(frame.memory[j].dissipation_j, y[6 + 3 * j], 2e-11);
            close(frame.memory[j].structural_work_j, y[7 + 3 * j], 2e-11);
            close(frame.memory[j].acoustic_work_j, y[8 + 3 * j], 2e-11);
        }
        assert_eq!(frame.memory[0].structural_work_j, 0.0);
        assert!(frame.memory[0].acoustic_work_j.abs() > 1e-4);
        assert!(frame.memory[1].structural_work_j.abs() > 1e-4);
        assert!(frame.memory[1].acoustic_work_j.abs() > 1e-4);
    }

    #[test]
    fn cancellation_invalid_input_and_clock_refusals_preserve_complete_state() {
        let source = Tetrahedron::new();
        let prepared = kernel::Prepared::new(
            source.source(),
            vec![arm()],
            0.2,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        let initial = prepared
            .state(3.0, &[0.0; 5], &[0.0; 5], &[0.3], &mut || false)
            .unwrap();
        let mut state = initial.clone();
        let mut workspace = prepared.workspace();
        let mut calls = 0;
        prepared
            .step_into(&mut state, &[0.1], &[0.0; 4], &mut workspace, &mut || {
                calls += 1;
                false
            })
            .unwrap();
        let mut state = initial.clone();
        let mut index = 0;
        let result = prepared.step_into(&mut state, &[0.1], &[0.0; 4], &mut workspace, &mut || {
            index += 1;
            index == calls
        });
        assert!(matches!(
            result,
            Err(StructuralInteriorRelaxationTimeError::Cancelled)
        ));
        assert_eq!(state.time_s, initial.time_s);
        assert_eq!(state.canonical, initial.canonical);
        assert_eq!(state.stress, initial.stress);
        assert_eq!(state.velocity, initial.velocity);
        assert!(
            prepared
                .step_into(
                    &mut state,
                    &[f64::NAN],
                    &[0.0; 4],
                    &mut workspace,
                    &mut || false
                )
                .is_err()
        );
        assert_eq!(state.canonical, initial.canonical);
        let mut late = prepared
            .state(1e30, &[0.0; 5], &[0.0; 5], &[0.3], &mut || false)
            .unwrap();
        let before = late.canonical.clone();
        assert!(
            prepared
                .step_into(&mut late, &[0.0], &[0.0; 4], &mut workspace, &mut || false)
                .is_err()
        );
        assert_eq!(late.time_s, 1e30);
        assert_eq!(late.canonical, before);
    }

    #[test]
    fn zero_branch_limit_preserves_original_compression_and_damped_forced_motion() {
        let mut source = Tetrahedron::new();
        source.damping[0] = 0.1;
        let dt = 0.43;
        let force = 0.7;
        let prepared = kernel::Prepared::new(
            source.source(),
            vec![],
            dt,
            StructuralInteriorRelaxationTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
        let mut state = prepared
            .state(0.0, &[0.0; 5], &[0.0; 5], &[], &mut || false)
            .unwrap();
        let mut workspace = prepared.workspace();
        let frame = prepared
            .step_into(&mut state, &[force], &[0.0; 4], &mut workspace, &mut || {
                false
            })
            .unwrap();
        // Uniform compression adds .3 to relaxed structural stiffness .7.
        let decay = 0.05;
        let frequency = det::sqrt(1.0 - decay * decay);
        let envelope = det::exp(-decay * dt);
        let q = force
            * (1.0
                - envelope
                    * (det::cos(frequency * dt) + decay / frequency * det::sin(frequency * dt)));
        let v = force / frequency * envelope * det::sin(frequency * dt);
        let integral = force * dt - v - 0.1 * q;
        close(state.canonical[0], q, 3e-11);
        close(state.velocity[0], v, 3e-11);
        for a in 1..5 {
            close(state.velocity[a], -0.6 * q, 3e-11);
            close(state.canonical[a], -0.6 * integral, 3e-11);
        }
        close(frame.work.structural_input_work_j, force * q, 3e-11);
        close(frame.work.fluid_to_structure_work_j, -0.15 * q * q, 3e-11);
        assert!(frame.memory.is_empty());
        assert!(state.stress.is_empty());
        assert_eq!(frame.work.relaxation_energy_change_j, 0.0);
        assert_eq!(frame.work.relaxation_dissipation_j, 0.0);
    }

    #[test]
    fn admission_bounds_and_active_resistance_remain_native_refusals() {
        let source = Tetrahedron::new();
        let limits = StructuralInteriorRelaxationTimeLimits {
            maximum_adapter_bytes: 0,
            ..Default::default()
        };
        assert!(matches!(
            kernel::Prepared::new(source.source(), vec![arm()], 0.2, limits, &mut || false),
            Err(StructuralInteriorRelaxationTimeError::Limit {
                quantity: "adapter bytes",
                ..
            })
        ));
        let mut limits = StructuralInteriorRelaxationTimeLimits::default();
        limits.maximum_branches = 0;
        assert!(
            kernel::Prepared::new(source.source(), vec![arm()], 0.2, limits, &mut || false)
                .is_err()
        );
        assert!(
            kernel::admit(usize::MAX, usize::MAX, 0.2, Default::default(), &mut || {
                false
            })
            .is_err()
        );
        assert!(matches!(
            kernel::Prepared::new(
                source.source(),
                vec![arm()],
                0.2,
                Default::default(),
                &mut || true
            ),
            Err(StructuralInteriorRelaxationTimeError::Cancelled)
        ));
        let mut active = Tetrahedron::new();
        active.damping[0] = -0.2;
        assert!(matches!(
            kernel::Prepared::new(
                active.source(),
                vec![arm()],
                0.2,
                Default::default(),
                &mut || false
            ),
            Err(StructuralInteriorRelaxationTimeError::Dynamics(_))
        ));
    }
}

mod original_fields {
    use super::super::*;
    use crate::cavity_modal::interior::structure::{NormalDisplacementTrace, StructuralMatrices};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
    use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
    use fs_feec::fixtures::single_tet;
    use fs_modal::harmonic::HarmonicResponseControls;
    use fs_sparse::Coo;

    #[test]
    fn zero_memory_matches_existing_full_field_time_and_source_binding_is_enforced() {
        let (complex, points) = single_tet();
        let medium = Medium {
            density_kg_m3: 0.2,
            sound_speed_m_s: 1.0,
        };
        let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
        let interior = HeterogeneousHelmholtzFe::new(
            &geometry,
            3,
            &mut |_, _| Ok::<_, AcousticFeError>(medium),
            &mut || false,
        )
        .unwrap();
        let traces: Vec<NormalDisplacementTrace> = geometry
            .boundary_faces()
            .into_iter()
            .map(|face| {
                let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
                (
                    face,
                    vertices.map(|p| (0..3).map(|axis| (p[axis] - 0.25) * normal[axis]).sum()),
                )
            })
            .collect();
        let scalar = |v| {
            let mut matrix = Coo::new(1, 1);
            matrix.push(0, 0, v);
            matrix.assemble()
        };
        let mass = scalar(1.0);
        let stiffness = scalar(0.7);
        let damping = scalar(0.1);
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 94,
                    kernel_id: 23,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let controls = HarmonicResponseControls {
                maximum_coordinates: 128,
                maximum_payload_bytes: 16 * 1024 * 1024,
                maximum_relative_residual: 1e-10,
            };
            let model = StructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: Some(&damping),
                },
                &[traces],
                &[],
                controls,
                &cx,
            )
            .unwrap();
            let dt = 0.3;
            let old = model
                .prepare_exact_time(dt, Default::default(), &cx)
                .unwrap();
            let new = model
                .prepare_relaxation_exact_time(&[], dt, Default::default(), &cx)
                .unwrap();
            let old_state = old
                .state(0.0, &[0.2], &[-0.1], &[0.05; 4], &[-0.12; 4], &cx)
                .unwrap();
            let new_state = new
                .relaxed_state(0.0, &[0.2], &[-0.1], &[0.05; 4], &[-0.12; 4], &cx)
                .unwrap();
            let a = old.step(&old_state, &[0.3], &[0.02; 4], &cx).unwrap();
            let b = new.step(&new_state, &[0.3], &[0.02; 4], &cx).unwrap();
            for (a, b) in a
                .state
                .structural_displacement()
                .iter()
                .chain(a.state.structural_velocity())
                .chain(a.state.pressure_potential_pa_s())
                .chain(a.state.pressure_pa())
                .zip(
                    b.state
                        .structural_displacement()
                        .iter()
                        .chain(b.state.structural_velocity())
                        .chain(b.state.pressure_potential_pa_s())
                        .chain(b.state.pressure_pa()),
                )
            {
                assert!((a - b).abs() < 3e-10);
            }
            assert!(
                (a.frame.work.structural_input_work_j - b.frame.work.structural_input_work_j).abs()
                    < 3e-10
            );
            assert!(
                (a.frame.work.fluid_to_structure_work_j - b.frame.work.fluid_to_structure_work_j)
                    .abs()
                    < 3e-10
            );
            assert_eq!(b.frame.work.relaxation_dissipation_j, 0.0);
            let branches = [MaxwellDampingBranch {
                strain: vec![1.0, 0.0, 0.0, 0.0, 0.0],
                stiffness: 1.25,
                relaxation_time_s: 1.0,
            }];
            let other_branches = branches.clone();
            let original = model
                .prepare_relaxation_exact_time(&branches, dt, Default::default(), &cx)
                .unwrap();
            let other = model
                .prepare_relaxation_exact_time(&other_branches, dt, Default::default(), &cx)
                .unwrap();
            let mut state = original
                .state(0.0, &[0.0], &[0.0], &[0.0; 4], &[0.0; 4], &[0.4], &cx)
                .unwrap();
            let mut workspace = other.workspace();
            assert!(
                other
                    .step_into(&mut state, &[0.0], &[0.0; 4], &mut workspace, &cx)
                    .is_err()
            );
            assert_eq!(state.time_s(), 0.0);
            assert_eq!(state.branch_stresses(), &[0.4]);
            let mut workspace = original.workspace();
            original
                .step_into(&mut state, &[0.0], &[0.0; 4], &mut workspace, &cx)
                .unwrap();
            assert!(state.structural_displacement()[0].abs() > 1e-4);
        });
    }
}
