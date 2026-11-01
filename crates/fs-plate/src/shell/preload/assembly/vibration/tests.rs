use super::super::tests::{Part, budget as tie_budget, options as preload_options, seam};
use super::super::{ShellPreloadComponent, equilibrate_shell_assembly};
use super::*;
use crate::shell::ShellMesh;
use crate::shell::junction::{GlueLine, TieKind, TieOffset};
use crate::shell::preload::tests::reference_newton;
use fs_exec::CancelGate;
use fs_sparse::{Coo, DirectOrdering, LdltOptions, SymbolicLdlt};

pub(super) fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 65,
                kernel_id: 65,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

pub(super) fn modal_options() -> SliceOptions {
    SliceOptions {
        max_lanczos: 512,
        max_restarts: 16,
        // The inverse-operator Ritz metric is distinct from the original
        // shell pencil's residual; the wide stiffness spectrum needs this
        // tighter iteration target to meet the unchanged 1e-8 physical gate.
        ritz_tol: 1e-14,
        mass_diagonal_equilibration: true,
        ldlt: LdltOptions {
            max_factor_bytes: Some(8 * 1024 * 1024),
            ..LdltOptions::default()
        },
        ..SliceOptions::default()
    }
}

pub(super) fn output_budget() -> AssemblyVibrationBudget {
    AssemblyVibrationBudget {
        maximum_modes: 16,
        maximum_damping_entries: 256,
    }
}

// Independent reference for these small fixtures and their empty window gaps;
// all returned production modes and residuals come from the sparse slicer.
pub(super) fn dense_reference(body: &PreloadedShellAssembly) -> Vec<ModePair> {
    let model = &body.assembly.model;
    fs_modal::eigh_gen_dense(&model.k.to_dense(), &model.m.to_dense(), model.free).unwrap()
}

pub(super) fn plan(body: &PreloadedShellAssembly, count: usize) -> WindowPlan {
    let modes = dense_reference(body);
    let split = count / 2;
    let middle = 0.5 * (modes[split - 1].lambda + modes[split].lambda);
    let high = 0.5 * (modes[count - 1].lambda + modes[count].lambda);
    WindowPlan {
        windows: vec![(0., middle), (middle, high)],
        threads: 2,
        maximum_windows: 2,
        maximum_bytes: 64 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}

pub(super) fn solve(parts: &[ShellPreloadComponent<'_>]) -> PreloadedShellAssembly {
    equilibrate_shell_assembly(
        parts,
        &[],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap()
}

pub(super) fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1. + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

fn check_source_residuals(body: &PreloadedShellAssembly, vibration: &AssemblyVibration<'_>) {
    let model = &body.assembly.model;
    let mass = SymbolicLdlt::analyze(&model.m, DirectOrdering::Amd)
        .unwrap()
        .factor(&model.m, &LdltOptions::default())
        .unwrap();
    let mut kphi = vec![0.; model.free];
    let mut mphi = kphi.clone();
    for mode in vibration.modes() {
        model.k.spmv(&mode.phi, &mut kphi);
        model.m.spmv(&mode.phi, &mut mphi);
        near(
            mode.phi.iter().zip(&mphi).map(|(a, b)| a * b).sum(),
            1.,
            1e-12,
        );
        let residual: Vec<_> = kphi
            .iter()
            .zip(&mphi)
            .map(|(k, m)| k - mode.lambda * m)
            .collect();
        let solved = mass.solve(&residual);
        let observed = residual
            .iter()
            .zip(solved)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .sqrt();
        near(observed, mode.residual, 1e-9);
        assert!(observed / mode.lambda < 1e-8);
        assert_eq!(
            mode.interval,
            (mode.lambda - mode.residual, mode.lambda + mode.residual)
        );
    }
}

#[test]
fn tensile_prestress_and_loaded_shape_survive_the_certified_vibration_consumer() {
    let mut part = Part::rectangle(0, 4, 2, 0., true);
    let unloaded = solve(&[part.input()]);
    // Trapezoid-integrated axial traction of total 100 N on the free edge.
    for (node, force) in [(4, 25.), (9, 50.), (14, 25.)] {
        part.loads[6 * node] = force;
    }
    let loaded = solve(&[part.input()]);
    assert!(loaded.displacements[0].iter().any(|x| x.abs() > 1e-7));
    let windows = plan(&loaded, 4);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let vibration =
            eigen_windows(&loaded, &windows, &modal_options(), output_budget(), cx).unwrap();
        let reference = dense_reference(&unloaded);
        assert_eq!(vibration.mode_count(), 4);
        // Resolve the physical shift even after both numerical intervals are
        // included, using an independent reference rather than comparing the
        // sparse method against itself on the very soft unloaded pencil.
        assert!(vibration.modes().next().unwrap().interval.0 > 1.01 * reference[0].interval.1);
        assert_eq!(
            vibration.displacements(&vec![0.; 4], cx).unwrap(),
            loaded.displacements
        );
        assert!(
            vibration
                .perturbation_displacements(&vec![0.; 4], cx)
                .unwrap()
                .iter()
                .flatten()
                .all(|x| *x == 0.)
        );
        check_source_residuals(&loaded, &vibration);
        // The consumer calls the original owner on the EXACT accepted pencil;
        // neither a new eigenproblem nor a second prestress term is introduced.
        let direct = slice_windows(
            &loaded.assembly.model.k,
            &loaded.assembly.model.m,
            &windows,
            &modal_options(),
            cx,
        )
        .unwrap();
        for (a, b) in vibration
            .modes()
            .zip(direct.reports.iter().flat_map(|r| &r.modes))
        {
            assert_eq!(a.lambda.to_bits(), b.lambda.to_bits());
            assert_eq!(a.phi, b.phi);
            assert_eq!(a.interval, b.interval);
        }
        // Deliberate falsifier: rebuild the deformed mesh as stress free. Its
        // spectrum must differ, even though its displayed shape is identical.
        let reset = Part {
            mesh: ShellMesh::new(
                part.mesh
                    .nodes
                    .iter()
                    .enumerate()
                    .map(|(node, p)| {
                        std::array::from_fn(|c| p[c] + loaded.displacements[0][6 * node + c])
                    })
                    .collect(),
                part.mesh.tris.clone(),
            )
            .unwrap(),
            sections: part.sections.clone(),
            fixed: part.fixed.clone(),
            thickness: part.thickness.clone(),
            loads: vec![0.; part.loads.len()],
        };
        let reset = solve(&[reset.input()]);
        let reset_modes = dense_reference(&reset);
        assert!(vibration.modes().next().unwrap().interval.0 > 1.01 * reset_modes[0].interval.1);
    });
}

#[test]
fn loaded_tied_body_matches_monolithic_modes_and_preserves_window_counts() {
    let mut whole = Part::rectangle(0, 4, 2, 0., true);
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 2, 0., false);
    right.loads[6 * 3 + 2] = -0.25;
    right.loads[6 * 5 + 2] = -0.5;
    whole.loads[6 * 7 + 2] = -0.25;
    whole.loads[6 * 9 + 2] = -0.5;
    let monolithic = solve(&[whole.input()]);
    let tied = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(2, TieOffset::Midsurface)],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let windows = plan(&monolithic, 4);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let a = eigen_windows(&tied, &windows, &modal_options(), output_budget(), cx).unwrap();
        let b =
            eigen_windows(&monolithic, &windows, &modal_options(), output_budget(), cx).unwrap();
        assert_eq!(a.mode_count(), b.mode_count());
        for (a, b) in a.modes().zip(b.modes()) {
            near(a.lambda, b.lambda, 1e-10);
        }
        check_source_residuals(&tied, &a);
        let reports = &a.spectrum().reports;
        assert_eq!(reports[0].below_high, reports[1].below_low);
        assert_eq!(a.mode_count(), reports[1].below_high - reports[0].below_low);
        assert_eq!(a.spectrum().receipt.pencil_symbolic_analyses, 1);
        for r in reports {
            assert_eq!(r.expected, r.modes.len());
        }
        let serial = WindowPlan {
            threads: 1,
            ..windows.clone()
        };
        let serial = eigen_windows(&tied, &serial, &modal_options(), output_budget(), cx).unwrap();
        for (a, b) in a.modes().zip(serial.modes()) {
            assert_eq!(a.lambda.to_bits(), b.lambda.to_bits());
            assert_eq!(a.phi, b.phi);
            assert_eq!(a.interval, b.interval);
        }
    });
}

#[test]
fn nonmatching_slave_loads_and_physical_modal_motion_are_work_conjugate() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 3, 0., false);
    right.loads[6 * 8 + 2] = -0.1;
    let body = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Midsurface)],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let windows = plan(&body, 4);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let mut forces: Vec<Vec<f64>> = [&left, &right]
            .iter()
            .map(|p| {
                (0..p.loads.len())
                    .map(|i| ((i * 7 % 11) as f64 - 5.) * 0.1)
                    .collect()
            })
            .collect();
        forces[1][6 * 3 + 2] += 1.; // Explicit load on an eliminated slave.
        let loads = vibration
            .project_loads(&forces.iter().map(Vec::as_slice).collect::<Vec<_>>(), cx)
            .unwrap();
        let q: Vec<_> = (0..vibration.mode_count())
            .map(|i| (i + 1) as f64 * 1e-5)
            .collect();
        let physical = vibration.perturbation_displacements(&q, cx).unwrap();
        let physical_work: f64 = forces
            .iter()
            .flatten()
            .zip(physical.iter().flatten())
            .map(|(a, b)| a * b)
            .sum();
        let modal_work: f64 = loads.iter().zip(&q).map(|(a, b)| a * b).sum();
        near(physical_work, modal_work, 1e-12);
        let total = vibration.displacements(&q, cx).unwrap();
        for ((total, base), increment) in total
            .iter()
            .flatten()
            .zip(body.displacements.iter().flatten())
            .zip(physical.iter().flatten())
        {
            near(*total - base, *increment, 1e-13);
        }
        let mut first = vec![0.; vibration.mode_count()];
        first[0] = 1.;
        assert_eq!(
            vibration.expand_mode(0, cx).unwrap(),
            vibration.perturbation_displacements(&first, cx).unwrap()
        );
        assert!(vibration.expand_mode(vibration.mode_count(), cx).is_err());
        assert!(vibration.project_loads(&[&forces[0]], cx).is_err());
        assert!(vibration.displacements(&[f64::NAN], cx).is_err());
    });
}

#[test]
fn compliant_glue_keeps_cross_mode_damping_and_exact_physical_dissipation() {
    let mut left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 2, 0., true);
    left.loads[6 * 5 + 2] = -0.01;
    let mut tie = seam(2, TieOffset::Midsurface);
    tie.kind = TieKind::Compliant(GlueLine {
        normal_n_m3: 1e5,
        shear_n_m3: [1e5, 2e5],
        viscosity_n_s_m3: Some([1e3, 2e3, 3e3]),
        width_m: 0.01,
    });
    let body = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[tie],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let windows = plan(&body, 6);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let count = vibration.mode_count();
        let c = vibration.glue_damping();
        assert_eq!(c.len(), count * count);
        assert!((0..count).any(|i| (0..i).any(|j| c[i * count + j].abs() > 1e-8)));
        let direct = body
            .assembly
            .project_glue_damping(&vibration.modes().map(|m| m.phi.clone()).collect::<Vec<_>>())
            .unwrap();
        for (a, b) in c.iter().zip(direct) {
            near(*a, b, 1e-12);
        }
        let velocity: Vec<_> = (0..count).map(|i| (i + 1) as f64 * 0.001).collect();
        let mut retained = vec![0.; body.assembly.model.free];
        for (mode, v) in vibration.modes().zip(&velocity) {
            for (q, phi) in retained.iter_mut().zip(&mode.phi) {
                *q += v * phi;
            }
        }
        let physical = body.assembly.glue_dissipation_rate(&retained).unwrap();
        let mut projected = 0.;
        for i in 0..count {
            for j in 0..count {
                assert_eq!(c[i * count + j], c[j * count + i]);
                projected += velocity[i] * c[i * count + j] * velocity[j];
            }
        }
        assert!(physical > 0.);
        near(physical, projected, 1e-12);
    });
}

#[test]
fn output_caps_cancellation_and_nonpositive_modes_publish_no_partial_report() {
    let part = Part::rectangle(0, 2, 2, 0., true);
    let mut body = solve(&[part.input()]);
    let windows = plan(&body, 4);
    let at_cap = AssemblyVibrationBudget {
        maximum_modes: 4,
        maximum_damping_entries: 16,
    };
    let access_gate = CancelGate::new_clock_free();
    with_cx(&access_gate, |cx| {
        let accepted = eigen_windows(&body, &windows, &modal_options(), at_cap, cx).unwrap();
        assert_eq!(accepted.mode_count(), 4);
        let too_many = AssemblyVibrationBudget {
            maximum_modes: 3,
            ..at_cap
        };
        assert!(
            eigen_windows(&body, &windows, &modal_options(), too_many, cx)
                .unwrap_err()
                .contains("mode cap")
        );
        let too_large = AssemblyVibrationBudget {
            maximum_damping_entries: 15,
            ..at_cap
        };
        assert!(
            eigen_windows(&body, &windows, &modal_options(), too_large, cx)
                .unwrap_err()
                .contains("projection budget")
        );
        access_gate.request();
        let zero = vec![0.; accepted.mode_count()];
        assert!(
            accepted
                .project_loads(&[part.loads.as_slice()], cx)
                .unwrap_err()
                .contains("cancelled")
        );
        assert!(
            accepted
                .perturbation_displacements(&zero, cx)
                .unwrap_err()
                .contains("cancelled")
        );
        assert!(
            accepted
                .displacements(&zero, cx)
                .unwrap_err()
                .contains("cancelled")
        );
        assert!(
            accepted
                .expand_mode(0, cx)
                .unwrap_err()
                .contains("cancelled")
        );
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let error = with_cx(&cancelled, |cx| {
        eigen_windows(&body, &windows, &modal_options(), at_cap, cx)
    })
    .unwrap_err();
    assert!(error.contains("cancelled"));
    // No mode clipping can turn a wrongly altered unstable pencil into audible
    // vibrations. The original equilibrium constructor itself refuses this K.
    let mut negative = Coo::new(body.assembly.model.free, body.assembly.model.free);
    for i in 0..body.assembly.model.free {
        let (columns, values) = body.assembly.model.k.row(i);
        for (&j, &value) in columns.iter().zip(values) {
            negative.push(i, j, -value);
        }
    }
    body.assembly.model.k = negative.assemble();
    let negative_windows = WindowPlan {
        windows: vec![(-windows.windows[1].1, 0.)],
        ..windows
    };
    let error = with_cx(&CancelGate::new_clock_free(), |cx| {
        eigen_windows(&body, &negative_windows, &modal_options(), at_cap, cx)
    })
    .unwrap_err();
    assert!(error.contains("nonpositive"), "{error}");
}
