//! Relaxed static strain reaches the original quadrature and compatible field.
use super::*;

fn real_strain(tones: &[(f64, [C64; 6])], time: f64) -> [f64; 6] {
    core::array::from_fn(|i| {
        tones
            .iter()
            .map(|(w, e)| {
                let (s, c) = (w * time).sin_cos();
                e[i].re * c + e[i].im * s
            })
            .sum()
    })
}
fn biased_storage(law: &TensorProny, tones: &[(f64, [C64; 6])], bias: &[f64; 6], time: f64) -> f64 {
    let stress = apply(law.c_eq(), bias);
    instantaneous_spectrum(law, tones, time)[0]
        + dot(&stress, &real_strain(tones, time))
        + 0.5 * dot(&stress, bias)
}

#[test]
fn affine_bias_keeps_original_support_order_and_independent_material_work() {
    with_cx(|cx, _| {
        let law = law().with_band((0.0, 100.0)).unwrap();
        let source = solid(tet_mesh(), law.clone(), cx);
        let tones = [
            (7.0, strain()),
            (8.7, strain().map(|z| z * C64::new(0.6, -0.2))),
        ];
        let frequencies = tones.map(|(w, _)| w);
        let displacement = tones.map(|(_, e)| affine(&TET, e));
        let input = displacement.each_ref().map(Vec::as_slice);
        let bias = [0.003, -0.002, 0.001, 0.0004, -0.0002, 0.0003];
        let static_displacement: Vec<_> = affine(&TET, bias.map(C64::from_re))
            .iter()
            .map(|z| z.re)
            .collect();
        let observer = source
            .prepare_harmonic_material_spectrum_observer(&frequencies, spectrum_budget(), cx)
            .unwrap();
        let window = HarmonicWindow::new(0.137, 0.41).unwrap();
        let field = observer
            .observe_window_with_bias(&input, &static_displacement, window, cx)
            .unwrap();
        let ac = observer.observe_window(&input, window, cx).unwrap();
        let volume = 1.0 / 6.0;
        let end = window.start_s() + window.duration_s();
        let start_energy = volume * biased_storage(&law, &tones, &bias, window.start_s());
        let end_energy = volume * biased_storage(&law, &tones, &bias, end);
        let begin_strain = real_strain(&tones, window.start_s());
        let end_strain = real_strain(&tones, end);
        let increment = core::array::from_fn(|i| end_strain[i] - begin_strain[i]);
        let expected_work = volume
            * (integrate_spectrum(&law, &tones, window)[2]
                + dot(&apply(law.c_eq(), &bias), &increment));
        near(field.start_stored_energy_j, start_energy, start_energy);
        near(field.end_stored_energy_j, end_energy, end_energy);
        near(
            field.stored_energy_change_j,
            end_energy - start_energy,
            start_energy.max(end_energy),
        );
        near(
            field.stress_work_j,
            expected_work,
            start_energy.max(end_energy),
        );
        near(
            field.energy_balance_residual_j,
            0.0,
            start_energy.max(end_energy),
        );
        assert_eq!(field.material_heat_j, ac.material_heat_j);
        for (a, b) in field.points[0].bias_mandel.iter().zip(bias) {
            near(*a, b, 0.003);
        }
        assert!(field.points[0].compatible_bias_mandel.is_none());

        // A retained moving coordinate is addressed by represented order, not
        // by indexing the complete mesh vector with a reduced row number.
        let material = TetViscoMaterial {
            law,
            density_kg_m3: 60.0,
            material_state_identity: [81; 32],
        };
        let supported = ViscoElements::new(
            tet_mesh(),
            TetViscoField::Uniform(&material),
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 1000,
            },
            cx,
        )
        .unwrap();
        let constrained = supported
            .prepare_harmonic_material_spectrum_observer(&[7.0], spectrum_budget(), cx)
            .unwrap();
        let q = C64::new(0.001, -0.0007);
        let reduced = constrained
            .observe_window_with_bias(&[&[q]], &[0.003], window, cx)
            .unwrap();
        let mut full_q = [C64::ZERO; 12];
        full_q[11] = q;
        let mut full_bias = [0.0; 12];
        full_bias[11] = 0.003;
        let full = source
            .prepare_harmonic_material_spectrum_observer(&[7.0], spectrum_budget(), cx)
            .unwrap()
            .observe_window_with_bias(&[&full_q], &full_bias, window, cx)
            .unwrap();
        assert_eq!(reduced, full);
    });
}

#[test]
fn bbar_static_bias_retains_compatible_strain_and_cannot_hide_its_domain_violation() {
    with_cx(|cx, _| {
        let law = law().with_band((0.0, 100.0)).unwrap();
        let source = solid(hex_mesh(), law.clone(), cx);
        let frequencies = [4.0, 6.3];
        let amplitudes = [(0.003, 0.001), (0.002, -0.0007)];
        let displacement = amplitudes.map(|(a, b)| {
            HEX.iter()
                .flat_map(|p| [C64::new(a * p[0] * p[1], b * p[0]), C64::ZERO, C64::ZERO])
                .collect::<Vec<_>>()
        });
        let input = displacement.each_ref().map(Vec::as_slice);
        let static_displacement = |a| {
            HEX.iter()
                .flat_map(|p| [a * p[0] * p[1], 0.0, 0.0])
                .collect::<Vec<_>>()
        };
        let bias_amplitude = 0.004;
        let observer = source
            .prepare_harmonic_material_spectrum_observer(&frequencies, spectrum_budget(), cx)
            .unwrap();
        let window = HarmonicWindow::new(-0.12, 0.39).unwrap();
        let field = observer
            .observe_window_with_bias(&input, &static_displacement(bias_amplitude), window, cx)
            .unwrap();
        let ac = observer.observe_window(&input, window, cx).unwrap();
        assert_eq!(field.material_heat_j, ac.material_heat_j);
        let g = 1.0 / 3.0_f64.sqrt();
        let (mut start, mut end) = (0.0, 0.0);
        for x in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
            for y in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
                let tones = core::array::from_fn::<_, 2, _>(|i| {
                    (
                        frequencies[i],
                        bbar_strain(x, y, amplitudes[i].0, amplitudes[i].1).0,
                    )
                });
                let bias = bbar_strain(x, y, bias_amplitude, 0.0).0.map(|z| z.re);
                start += 0.25 * biased_storage(&law, &tones, &bias, window.start_s());
                end += 0.25
                    * biased_storage(&law, &tones, &bias, window.start_s() + window.duration_s());
            }
        }
        near(field.start_stored_energy_j, start, start);
        near(field.end_stored_energy_j, end, end);
        near(field.stored_energy_change_j, end - start, start.max(end));
        near(field.energy_balance_residual_j, 0.0, start.max(end));
        for point in &field.points {
            let raw = point.compatible_bias_mandel.unwrap();
            near(
                point.bias_mandel[..3].iter().sum(),
                0.5 * bias_amplitude,
                bias_amplitude,
            );
            near(
                raw[0] + 3.0 * point.bias_mandel[1],
                0.5 * bias_amplitude,
                bias_amplitude,
            );
        }
        let zero = vec![C64::ZERO; 24];
        assert!(
            observer
                .observe_window_with_bias(&[&zero, &zero], &static_displacement(0.0215), window, cx)
                .is_err()
        );
        assert_eq!(
            field,
            observer
                .observe_window_with_bias(&input, &static_displacement(bias_amplitude), window, cx)
                .unwrap()
        );
    });
}

#[test]
fn zero_strain_bias_preserves_ac_only_band_and_refusals_leave_observer_reusable() {
    with_cx(|cx, gate| {
        let source = solid(tet_mesh(), law(), cx);
        let observer = source
            .prepare_harmonic_material_spectrum_observer(&[7.0], spectrum_budget(), cx)
            .unwrap();
        let q = affine(&TET, strain());
        let window = HarmonicWindow::new(0.13, 0.2).unwrap();
        let ac = observer.observe_window(&[&q], window, cx).unwrap();
        assert_eq!(
            ac,
            observer
                .observe_window_with_bias(&[&q], &[0.0; 12], window, cx)
                .unwrap()
        );
        let rigid: Vec<_> = TET.iter().flat_map(|_| [1e80, -3e80, 4e80]).collect();
        assert_eq!(
            ac,
            observer
                .observe_window_with_bias(&[&q], &rigid, window, cx)
                .unwrap()
        );
        let bias: Vec<_> = affine(&TET, [C64::from_re(0.001); 6])
            .iter()
            .map(|z| z.re)
            .collect();
        assert!(
            observer
                .observe_window_with_bias(&[&q], &bias, window, cx)
                .is_err()
        );
        assert!(
            observer
                .observe_window_with_bias(&[&q], &[0.0; 11], window, cx)
                .is_err()
        );
        assert!(
            observer
                .observe_window_with_bias(&[&q], &[f64::NAN; 12], window, cx)
                .is_err()
        );
        assert_eq!(ac, observer.observe_window(&[&q], window, cx).unwrap());
        gate.request();
        assert!(
            observer
                .observe_window_with_bias(&[&q], &[0.0; 12], window, cx)
                .is_err()
        );
    });
}
