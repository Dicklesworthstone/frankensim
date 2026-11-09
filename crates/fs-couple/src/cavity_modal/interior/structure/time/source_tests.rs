//! Physical startup, original hereditary heat and full P2 sparse source loads.

use super::*;
use crate::cavity_modal::interior::structure::sparse::{
    PreparedSparseStructuralInterior, SparseStructuralTimeLimits,
};
use crate::cavity_modal::interior::structure::time::relaxation::{
    MaxwellDampingBranch, StructuralInteriorRelaxationTimeLimits,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{HeterogeneousError, PressureSourceLoadLimits};
use std::cell::Cell;

#[test]
fn held_volume_sources_reproduce_exact_startup_and_negative_decompression_work() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let calls = Cell::new(0);
    let volume = interior
        .volume_velocity_source(
            4,
            &mut |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, HeterogeneousError>(6.0)
            },
            &mut || false,
        )
        .unwrap();
    let sampled = calls.get();
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[7.0], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: None,
            },
            &[Vec::new()],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let h = 0.2;
        let original = model
            .prepare_exact_time(h, StructuralInteriorTimeLimits::default(), cx)
            .unwrap();
        let driven = original
            .with_sources(&[], &[&volume], PressureSourceLoadLimits::default(), cx)
            .unwrap();
        let initial = original
            .state(0.0, &[0.0], &[0.0], &[0.0; 10], &[0.0; 10], cx)
            .unwrap();
        let mut state = initial.clone();
        let mut workspace = driven.workspace();
        let q = 0.02;
        // beta*V*p_dot=Q, with V=1/6, so p=Q*t/(beta*V).
        let compression = 6.0 * medium.density_kg_m3 * medium.sound_speed_m_s.powi(2);
        let first = driven
            .step_into(&mut state, &[0.0], &[], &[q], &mut workspace, cx)
            .unwrap();
        let first_work = first.work.acoustic_input_work_j;
        close(first_work, 0.5 * compression * (q * h).powi(2), 2e-11);
        close(
            first.work.final_energy.fluid_compression_j,
            first_work,
            2e-11,
        );
        close(first.work.balance_residual_j(), 0.0, 2e-11);
        for (i, (&p, &psi)) in state
            .pressure_pa()
            .iter()
            .zip(state.pressure_potential_pa_s())
            .enumerate()
        {
            close(p, if i < 4 { compression * q * h } else { 0.0 }, 2e-10);
            close(
                psi,
                if i < 4 {
                    compression * q * h * h / 2.0
                } else {
                    0.0
                },
                2e-10,
            );
        }
        let second_q = -0.5 * q;
        let second = driven
            .step_into(&mut state, &[0.0], &[], &[second_q], &mut workspace, cx)
            .unwrap();
        let second_work = second.work.acoustic_input_work_j;
        let expected_work = second_q * compression * h * h * (q + 0.5 * second_q);
        close(second_work, expected_work, 2e-11);
        assert!(
            second_work < 0.0,
            "the source can absorb stored compression energy"
        );
        close(
            second.work.final_energy.fluid_compression_j,
            first_work + second_work,
            2e-11,
        );
        for (i, &p) in state.pressure_pa().iter().enumerate() {
            close(
                p,
                if i < 4 {
                    compression * (q + second_q) * h
                } else {
                    0.0
                },
                2e-10,
            );
        }
        close(state.time_s(), 2.0 * h, 1e-15);
        assert_eq!(calls.get(), sampled);

        // Admission fails before accepted state, history or clock can change.
        let accepted = state.clone();
        assert!(
            driven
                .step_into(&mut state, &[0.0], &[], &[f64::NAN], &mut workspace, cx)
                .is_err()
        );
        assert_eq!(state.time_s(), accepted.time_s());
        assert_eq!(state.pressure_pa(), accepted.pressure_pa());
        assert_eq!(
            state.pressure_potential_pa_s(),
            accepted.pressure_potential_pa_s()
        );
        let mut limits = PressureSourceLoadLimits::default();
        limits.maximum_products_per_projection = 9;
        assert!(matches!(
            original.with_sources(&[], &[&volume], limits, cx),
            Err(HeterogeneousError::Limit { .. })
        ));
        limits = PressureSourceLoadLimits::default();
        limits.maximum_payload_bytes = 159; // 10 cached + 10 scratch doubles need 160.
        assert!(matches!(
            original.with_sources(&[], &[&volume], limits, cx),
            Err(HeterogeneousError::Limit { .. })
        ));
        let other = uniform(&geometry, medium);
        let foreign = other
            .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
                false
            })
            .unwrap();
        assert!(
            original
                .with_sources(&[], &[&foreign], PressureSourceLoadLimits::default(), cx)
                .is_err()
        );
        let candidate = driven.step(&initial, &[0.0], &[], &[q], cx).unwrap();
        close(
            candidate.frame.work.acoustic_input_work_j,
            first_work,
            2e-11,
        );
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn original_injection_drives_structure_without_losing_active_or_dark_maxwell_histories() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let volume = interior
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    // Original affine dilation; the uniform pressure shape remains invariant.
    let traces: Vec<NormalDisplacementTrace> = geometry
        .boundary_faces()
        .into_iter()
        .map(|face| {
            let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
            (
                face,
                vertices.map(|point| (0..3).map(|d| (point[d] - 0.25) * normal[d]).sum()),
            )
        })
        .collect();
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[3.0], 1);
    let damping = dense(&[0.4], 1);
    let branches = [
        MaxwellDampingBranch {
            strain: vec![1.0, 0.0, 0.0, 0.0, 0.0],
            stiffness: 1.25,
            relaxation_time_s: 0.3,
        },
        MaxwellDampingBranch {
            strain: vec![0.0; 5],
            stiffness: 2.0,
            relaxation_time_s: 0.17,
        },
    ];
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[traces],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let h = 0.21;
        let original = model
            .prepare_relaxation_exact_time(
                &branches,
                h,
                StructuralInteriorRelaxationTimeLimits::default(),
                cx,
            )
            .unwrap();
        let driven = original
            .with_sources(&[], &[&volume], PressureSourceLoadLimits::default(), cx)
            .unwrap();
        let sigma0 = 0.4;
        let initial = original
            .state(
                0.0,
                &[0.0],
                &[0.0],
                &[0.0; 4],
                &[0.0; 4],
                &[0.0, sigma0],
                cx,
            )
            .unwrap();
        let mut state = initial.clone();
        let mut workspace = driven.workspace();
        let q = 0.03;
        let frame = driven
            .step_into(&mut state, &[0.0], &[], &[q], &mut workspace, cx)
            .unwrap();
        let source_work = q * frame.integrated_pressure_pa_s.iter().sum::<f64>() / 4.0;
        close(frame.work.acoustic_input_work_j, source_work, 2e-11);
        close(frame.work.structural_input_work_j, 0.0, 1e-15);
        close(frame.work.balance_residual_j(), 0.0, 2e-10);
        assert!(state.structural_displacement()[0] > 1e-6);
        assert!(state.branch_stresses()[0].abs() > 1e-6);
        assert!(frame.memory[0].dissipation_j > 0.0);
        // The spatially dark source history decays even while the acoustic
        // source drives the mechanical branch. Its heat has an independent ODE.
        let decay = det::exp(-h / branches[1].relaxation_time_s);
        close(state.branch_stresses()[1], sigma0 * decay, 2e-11);
        close(
            frame.memory[1].dissipation_j,
            sigma0 * sigma0 / (2.0 * branches[1].stiffness) * (1.0 - decay * decay),
            2e-11,
        );
        close(frame.memory[1].structural_work_j, 0.0, 2e-12);
        close(frame.memory[1].acoustic_work_j, 0.0, 2e-12);

        // Analytic original P1 volume moments: each coordinate receives Q/4.
        let explicit = original.step(&initial, &[0.0], &[q / 4.0; 4], cx).unwrap();
        for (&a, &b) in state
            .pressure_pa()
            .iter()
            .chain(state.structural_velocity())
            .chain(state.branch_stresses())
            .chain(state.branch_viscous_strains())
            .zip(
                explicit
                    .state
                    .pressure_pa()
                    .iter()
                    .chain(explicit.state.structural_velocity())
                    .chain(explicit.state.branch_stresses())
                    .chain(explicit.state.branch_viscous_strains()),
            )
        {
            close(a, b, 2e-10);
        }
        let cloned_branches = branches.clone();
        let foreign = model
            .prepare_relaxation_exact_time(
                &cloned_branches,
                h,
                StructuralInteriorRelaxationTimeLimits::default(),
                cx,
            )
            .unwrap();
        let foreign_driven = foreign
            .with_sources(&[], &[&volume], PressureSourceLoadLimits::default(), cx)
            .unwrap();
        let mut wrong_workspace = foreign_driven.workspace();
        let accepted = state.clone();
        assert!(
            driven
                .step_into(&mut state, &[0.0], &[], &[q], &mut wrong_workspace, cx)
                .is_err()
        );
        assert_eq!(state.time_s(), accepted.time_s());
        assert_eq!(state.branch_stresses(), accepted.branch_stresses());
        assert_eq!(
            state.branch_viscous_strains(),
            accepted.branch_viscous_strains()
        );
        assert_eq!(state.pressure_pa(), accepted.pressure_pa());
    });
}

#[test]
fn zero_net_p2_injection_reaches_sparse_full_fields_and_exact_source_work() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let boundary = interior
        .quadratic_velocity_source(&[([0, 1, 2], [0.0, 0.0, 0.0, 1.0, -1.0, 0.0])], &mut || {
            false
        })
        .unwrap();
    assert!(boundary.integrated_shape_m2().abs() < 2e-14);
    assert!(boundary.has_nonzero_projection());
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[7.0], 1);
    let damping = dense(&[0.2], 1);
    let traces = [vec![([0, 1, 2], [0.7, -0.3, 0.9])]];
    let walls = [ResistiveBoundary {
        faces: vec![[0, 1, 3]],
        admittance_m_pa_s: 0.03,
    }];
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let sparse = PreparedSparseStructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &traces,
            &walls,
            0.015,
            SparseStructuralTimeLimits::default(),
            &cx,
        )
        .unwrap();
        let driven = sparse
            .with_sources(&[&boundary], &[], PressureSourceLoadLimits::default(), &cx)
            .unwrap();
        let initial = sparse.zero_state();
        let amplitude = 0.02;
        let step = driven
            .step(&initial, &[0.0], &[amplitude], &[], &cx)
            .unwrap();
        // Independent P2 moments of 4*l0*l1-4*l1*l2 on area A=1/2.
        let mut load = vec![amplitude / 30.0, 0.0, -amplitude / 30.0, 0.0];
        for &edge in &complex.edges {
            load.push(match edge {
                [0, 1] => amplitude / 90.0,
                [1, 2] => -amplitude / 90.0,
                _ => 0.0,
            });
        }
        let explicit = sparse.step(&initial, &[0.0], &load, &cx).unwrap();
        for (&a, &b) in step
            .state
            .pressure_pa()
            .iter()
            .chain(step.state.structural_velocity())
            .zip(
                explicit
                    .state
                    .pressure_pa()
                    .iter()
                    .chain(explicit.state.structural_velocity()),
            )
        {
            close(a, b, 2e-11);
        }
        assert!(
            step.state
                .pressure_pa()
                .iter()
                .skip(4)
                .any(|p| p.abs() > 1e-7)
        );
        assert!(step.state.structural_velocity()[0].abs() > 1e-8);
        // Initial pressure is zero: p_mid=p_end/2 is the actual time midpoint,
        // not the one-half factor used for harmonic peak amplitudes.
        let work = 0.015
            * load
                .iter()
                .zip(step.state.pressure_pa())
                .map(|(&b, &p)| b * p / 2.0)
                .sum::<f64>();
        close(step.work.acoustic_input_work_j, work, 2e-12);
        assert!(work > 0.0);
        close(step.work.balance_residual_j(), 0.0, 2e-11);
        close(step.work.interface_balance_residual_j(), 0.0, 2e-11);
        assert_eq!(initial.time_s(), 0.0);
        assert!(initial.pressure_pa().iter().all(|&p| p == 0.0));

        let original = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &traces,
            &walls,
            controls(),
            &cx,
        )
        .unwrap();
        let exact = original
            .prepare_exact_time(0.015, StructuralInteriorTimeLimits::default(), &cx)
            .unwrap();
        let exact_driven = exact
            .with_sources(&[&boundary], &[], PressureSourceLoadLimits::default(), &cx)
            .unwrap();
        let mut accepted = exact
            .state(0.0, &[0.0], &[0.0], &[0.0; 10], &[0.0; 10], &cx)
            .unwrap();
        let mut workspace = exact_driven.workspace();
        let frame = exact_driven
            .step_into(
                &mut accepted,
                &[0.0],
                &[amplitude],
                &[],
                &mut workspace,
                &cx,
            )
            .unwrap();
        let actual_work = load
            .iter()
            .zip(&frame.integrated_pressure_pa_s)
            .map(|(&b, &integral)| b * integral)
            .sum::<f64>();
        close(frame.work.acoustic_input_work_j, actual_work, 2e-12);
        let saved = accepted.clone();
        gate.request();
        assert!(
            exact_driven
                .step_into(
                    &mut accepted,
                    &[0.0],
                    &[amplitude],
                    &[],
                    &mut workspace,
                    &cx
                )
                .is_err()
        );
        assert_eq!(accepted.time_s(), saved.time_s());
        assert_eq!(accepted.pressure_pa(), saved.pressure_pa());
        assert_eq!(accepted.structural_velocity(), saved.structural_velocity());
        assert!(
            driven
                .step(&initial, &[0.0], &[amplitude], &[], &cx)
                .is_err()
        );
        assert_eq!(initial.time_s(), 0.0);
    });
}
