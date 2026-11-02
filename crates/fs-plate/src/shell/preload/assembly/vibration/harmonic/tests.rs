use super::super::super::equilibrate_shell_assembly;
use super::super::super::tests::{Part, budget as tie_budget, options as preload_options, seam};
use super::super::tests::{dense_reference, modal_options, output_budget, plan, solve, with_cx};
use super::super::{WindowPlan, eigen_windows};
use super::*;
use crate::shell::junction::{GlueLine, TieKind, TieOffset};
use crate::shell::preload::tests::reference_newton;
use fs_exec::CancelGate;

fn budget() -> AssemblyHarmonicBudget {
    AssemblyHarmonicBudget {
        maximum_modes: 16,
        maximum_projection_terms: 100_000_000,
        maximum_projection_bytes: 2 * 1024 * 1024,
    }
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 16,
        maximum_payload_bytes: 2 * 1024 * 1024,
        maximum_relative_residual: 1e-10,
    }
}

fn relative(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale.max(f64::MIN_POSITIVE),
        "{a} != {b}; scale {scale}"
    );
}

#[test]
fn g1_loaded_shell_resonance_matches_analytic_one_mode_and_keeps_its_mean_shape() {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    let unloaded = solve(&[part.input()]);
    for (node, force) in [(2, 25.), (5, 50.), (8, 25.)] {
        part.loads[6 * node] = force;
    }
    part.loads[6 * 8 + 2] = -0.01;
    let loaded = solve(&[part.input()]);
    let reference = dense_reference(&loaded);
    let windows = WindowPlan {
        windows: vec![(0., 0.5 * (reference[0].lambda + reference[1].lambda))],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 64 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&loaded, &windows, &modal_options(), output_budget(), cx).unwrap();
        let model = AssemblyHarmonicModel::new(&vibration, budget(), cx).unwrap();
        let (m, k) = (model.mass_matrix()[0], model.stiffness_matrix()[0]);
        let natural = (k / m).sqrt();
        assert!(natural.powi(2) > 1.01 * dense_reference(&unloaded)[0].lambda);
        let mean = loaded.displacements.clone();
        let shape = vibration.expand_mode(0, cx).unwrap();
        let mut loads = vec![C64::ZERO; part.loads.len()];
        loads[6 * 8 + 2] = C64::new(0.02, -0.003);
        let force = loads[6 * 8 + 2].scale(shape[0][6 * 8 + 2]);
        let mut below = C64::ZERO;
        for ratio in [0., 0.8, 0.999, 1.001, 1.2] {
            let w = natural * ratio;
            let prepared = model.prepare(w, controls(), cx).unwrap();
            let actual = prepared.solve(&[&loads], cx).unwrap();
            let expected = force.scale(1. / (k - w * w * m));
            let velocity = C64::new(0., -w) * expected;
            relative(
                (actual.modal.displacement[0] - expected).abs(),
                0.,
                expected.abs(),
                2e-10,
            );
            for j in 0..loads.len() {
                relative(
                    (actual.perturbation_displacements[0][j] - expected.scale(shape[0][j])).abs(),
                    0.,
                    expected.abs() * shape[0][j].abs(),
                    2e-10,
                );
                relative(
                    (actual.velocities[0][j] - velocity.scale(shape[0][j])).abs(),
                    0.,
                    velocity.abs() * shape[0][j].abs(),
                    2e-10,
                );
            }
            assert_eq!(loaded.displacements, mean);
            assert_eq!(vibration.displacements(&[0.], cx).unwrap(), mean);
            assert_eq!(actual.modal.average_viscous_power_w, 0.);
            if ratio == 0.999 {
                below = actual.modal.displacement[0];
            }
            if ratio == 1.001 {
                assert!((below / actual.modal.displacement[0]).re < -0.99);
                assert!(
                    actual.modal.displacement[0].abs() > 400. * (force / C64::from_re(k)).abs()
                );
            }
        }
        let zero = vec![C64::ZERO; loads.len()];
        let rest = model
            .prepare(0.5 * natural, controls(), cx)
            .unwrap()
            .solve(&[&zero], cx)
            .unwrap();
        assert!(
            rest.perturbation_displacements
                .iter()
                .flatten()
                .all(|z| *z == C64::ZERO)
        );
    });
}

#[test]
fn g1_full_glue_response_obeys_original_projected_equations_and_physical_cycle_power() {
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
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration = eigen_windows(
            &body,
            &plan(&body, 6),
            &modal_options(),
            output_budget(),
            cx,
        )
        .unwrap();
        let model = AssemblyHarmonicModel::new(&vibration, budget(), cx).unwrap();
        let n = vibration.mode_count();
        assert!((0..n).any(|i| (0..i).any(|j| model.damping_matrix()[i * n + j].abs() > 1e-8)));
        let w = 1.002 * vibration.modes().next().unwrap().lambda.sqrt();
        let mut forces = [
            vec![C64::ZERO; left.loads.len()],
            vec![C64::ZERO; right.loads.len()],
        ];
        forces[0][6 * 8 + 2] = C64::new(0.02, -0.003);
        forces[1][6 * 8 + 2] = C64::new(-0.01, 0.004);
        let loads = [&forces[0][..], &forces[1][..]];
        let prepared = model.prepare(w, controls(), cx).unwrap();
        let actual = prepared.solve(&loads, cx).unwrap();
        let again = prepared.solve(&loads, cx).unwrap();
        assert_eq!(actual.modal.displacement, again.modal.displacement);
        let qre: Vec<_> = actual.modal.displacement.iter().map(|v| v.re).collect();
        let qim: Vec<_> = actual.modal.displacement.iter().map(|v| v.im).collect();
        let qre = vibration.retained_perturbation(&qre, cx).unwrap();
        let qim = vibration.retained_perturbation(&qim, cx).unwrap();
        let vre: Vec<_> = qim.iter().map(|q| w * q).collect();
        let vim: Vec<_> = qre.iter().map(|q| -w * q).collect();
        let native = &body.assembly.model;
        let apply = |matrix: &fs_sparse::Csr, vector: &[f64]| {
            let mut out = vec![0.; vector.len()];
            matrix.spmv(vector, &mut out);
            out
        };
        let kre = apply(&native.k, &qre);
        let kim = apply(&native.k, &qim);
        let mre = apply(&native.m, &qre);
        let mim = apply(&native.m, &qim);
        let cre = apply(body.assembly.glue_damping.as_ref().unwrap(), &vre);
        let cim = apply(body.assembly.glue_damping.as_ref().unwrap(), &vim);
        let real: Vec<Vec<_>> = forces
            .iter()
            .map(|f| f.iter().map(|z| z.re).collect())
            .collect();
        let imag: Vec<Vec<_>> = forces
            .iter()
            .map(|f| f.iter().map(|z| z.im).collect())
            .collect();
        let real: Vec<_> = real.iter().map(Vec::as_slice).collect();
        let imag: Vec<_> = imag.iter().map(Vec::as_slice).collect();
        let fre = body.assembly.project_loads(&real).unwrap();
        let fim = body.assembly.project_loads(&imag).unwrap();
        for mode in vibration.modes() {
            for (k, m, c, f) in [(&kre, &mre, &cre, &fre), (&kim, &mim, &cim, &fim)] {
                let mut residual = 0.;
                let mut scale = 0.;
                for i in 0..native.free {
                    residual += mode.phi[i] * (k[i] - w * w * m[i] + c[i] - f[i]);
                    scale += mode.phi[i].abs()
                        * (k[i].abs() + w * w * m[i].abs() + c[i].abs() + f[i].abs());
                }
                relative(residual, 0., scale, 1e-10);
            }
        }
        let heat = 0.5
            * (body.assembly.glue_dissipation_rate(&vre).unwrap()
                + body.assembly.glue_dissipation_rate(&vim).unwrap());
        assert!(heat > 0.);
        relative(actual.modal.average_viscous_power_w, heat, heat, 2e-10);
        relative(actual.average_component_power_w, heat, heat, 2e-10);
        // Falsifier: discarding cross-mode glue damping changes the response.
        let mut diagonal = model.damping_matrix().to_vec();
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    diagonal[i * n + j] = 0.;
                }
            }
        }
        let modal_re = vibration.project_loads(&real, cx).unwrap();
        let modal_im = vibration.project_loads(&imag, cx).unwrap();
        let force: Vec<_> = modal_re
            .into_iter()
            .zip(modal_im)
            .map(|(re, im)| C64::new(re, im))
            .collect();
        let wrong = PreparedHarmonicResponse::try_new(
            model.mass_matrix(),
            &diagonal,
            model.stiffness_matrix(),
            n,
            w,
            controls(),
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        let difference = wrong
            .displacement
            .iter()
            .zip(&actual.modal.displacement)
            .map(|(&a, &b)| (a - b).abs())
            .fold(0., f64::max);
        let magnitude = actual
            .modal
            .displacement
            .iter()
            .map(|z| z.abs())
            .fold(0., f64::max);
        assert!(difference > 1e-5 * magnitude);
    });
}

#[test]
fn g0_caps_and_g4_cancellation_preserve_the_loaded_body_and_prepared_factor() {
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
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let vibration = eigen_windows(
            &body,
            &plan(&body, 4),
            &modal_options(),
            output_budget(),
            cx,
        )
        .unwrap();
        for (modes, terms, bytes) in [
            (3, u64::MAX, usize::MAX),
            (4, 0, usize::MAX),
            (4, u64::MAX, 1),
        ] {
            assert!(
                AssemblyHarmonicModel::new(
                    &vibration,
                    AssemblyHarmonicBudget {
                        maximum_modes: modes,
                        maximum_projection_terms: terms,
                        maximum_projection_bytes: bytes,
                    },
                    cx
                )
                .is_err()
            );
        }
        let model = AssemblyHarmonicModel::new(&vibration, budget(), cx).unwrap();
        let w = 0.8 * vibration.modes().next().unwrap().lambda.sqrt();
        let prepared = model.prepare(w, controls(), cx).unwrap();
        let zero = vec![C64::ZERO; left.loads.len()];
        let mut slave = vec![C64::ZERO; right.loads.len()];
        slave[6 * 3 + 2] = C64::new(0.1, -0.04);
        let before = body.displacements.clone();
        let expected = prepared.solve(&[&zero, &slave], cx).unwrap();
        assert!(
            expected
                .perturbation_displacements
                .iter()
                .flatten()
                .any(|z| z.abs() > 0.)
        );
        assert!(prepared.solve(&[&zero], cx).is_err());
        slave[0] = C64::new(f64::NAN, 0.);
        assert!(prepared.solve(&[&zero, &slave], cx).is_err());
        slave[0] = C64::ZERO;
        assert_eq!(
            prepared
                .solve(&[&zero, &slave], cx)
                .unwrap()
                .modal
                .displacement,
            expected.modal.displacement
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, |other| {
            assert!(prepared.solve(&[&zero, &slave], other).is_err());
            assert!(model.prepare(w, controls(), other).is_err());
        });
        assert_eq!(body.displacements, before);
        assert_eq!(
            prepared
                .solve(&[&zero, &slave], cx)
                .unwrap()
                .modal
                .displacement,
            expected.modal.displacement
        );
    });
}
