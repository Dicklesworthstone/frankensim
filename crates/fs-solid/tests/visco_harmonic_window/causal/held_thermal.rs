//! Original held thermal samples and their B-bar projection remain one history.

use super::*;
use fs_blake3::ContentHash;
use fs_solid::linear3::{
    QuadraticThermalStrainField, TetThermalStrainField, TetThermalStrainState,
};

fn sample(source: &ViscoElements, free: [f64; 6], ceiling: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        310.0,
        300.0,
        free,
        source
            .reference_material(0)
            .unwrap()
            .material_state_identity,
        ContentHash([83; 32]),
        ceiling,
    )
    .unwrap()
}

#[test]
fn free_expansion_preserves_ac_only_material_support_and_admits_the_entire_waveform() {
    with_cx(|cx, gate| {
        let source = solid(tet_mesh(), law(), cx); // Original band excludes DC.
        let free = [0.0008, 0.0, 0.0, 0.0, 0.0, 0.0];
        let thermal = sample(&source, free, 0.001);
        let bias = real_affine(&TET, free);
        let mut original = source.state(cx).unwrap();
        original
            .step_internal_force_with_thermal_strain(
                &bias,
                0.1,
                TetThermalStrainField::Uniform(&thermal),
                cx,
            )
            .unwrap();
        assert_eq!(original.energy(), ViscoEnergy::default());
        let accepted = original.bind_harmonic_clock(Dd::ZERO).unwrap();
        let cold = source
            .state(cx)
            .unwrap()
            .bind_harmonic_clock(Dd::ZERO)
            .unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[7.0], causal_budget(), cx)
            .unwrap();
        let mut e = [C64::ZERO; 6];
        e[0] = C64::new(0.0, 0.0001);
        let tone = affine(&TET, e);
        let next = observer
            .advance_window(&accepted, &[&tone], &bias, 0.2, cx)
            .unwrap();
        let reference = observer
            .advance_window(&cold, &[&tone], &[0.0; 12], 0.2, cx)
            .unwrap();
        let scale = reference
            .field
            .material_heat_j
            .max(reference.field.end_stored_energy_j);
        near(
            next.field.material_heat_j,
            reference.field.material_heat_j,
            scale,
        );
        near(
            next.field.stress_work_j,
            reference.field.stress_work_j,
            scale,
        );
        near(
            next.field.end_stored_energy_j,
            reference.field.end_stored_energy_j,
            scale,
        );
        assert!(next.field.material_heat_j > 0.0);
        assert_eq!(next.field.points[0].bias_mandel, [0.0; 6]);
        assert_eq!(
            next.history.history().thermal_strains(),
            Some(core::slice::from_ref(&thermal))
        );
        assert_eq!(next.history.history().eigenstrain_work_j(), 0.0);
        assert_eq!(accepted.history().energy(), ViscoEnergy::default());

        // Starting total/mechanical strains fit, but a later phase exceeds
        // the actual thermal source ceiling while the material norm still fits.
        e[0] = C64::new(0.0, 0.0004);
        let excessive = affine(&TET, e);
        assert!(
            observer
                .advance_window(&accepted, &[&excessive], &bias, 0.2, cx)
                .is_err()
        );
        observer
            .advance_window(&cold, &[&excessive], &[0.0; 12], 0.2, cx)
            .unwrap();
        assert_eq!(
            accepted.history().thermal_strains(),
            Some(core::slice::from_ref(&thermal))
        );
        assert_eq!(accepted.history().displacement_m(), bias);
        gate.request();
        assert!(
            observer
                .advance_window(&accepted, &[&tone], &bias, 0.2, cx)
                .is_err()
        );
        assert_eq!(accepted.history().energy(), ViscoEnergy::default());
        assert_eq!(
            accepted.history().thermal_strains(),
            Some(core::slice::from_ref(&thermal))
        );
    });
}

#[test]
fn nonuniform_bbar_expansion_keeps_raw_samples_projected_memories_and_original_work() {
    with_cx(|cx, _| {
        let law = law().with_band((0.0, 100.0)).unwrap();
        let source = solid(hex_mesh(), law.clone(), cx);
        let samples: Vec<_> = (0..8)
            .map(|q| {
                let q = q as f64;
                sample(
                    &source,
                    [
                        0.001 + 0.00008 * q,
                        -0.0003 + 0.00002 * q,
                        0.0002,
                        0.00003 * q,
                        0.0,
                        -0.00001 * q,
                    ],
                    0.005,
                )
            })
            .collect();
        let ramp = 0.1;
        let mut original = source.state(cx).unwrap();
        original
            .step_internal_force_with_thermal_strain(
                &[0.0; 24],
                ramp,
                QuadraticThermalStrainField::AtQuadrature {
                    states: &samples,
                    points_per_element: 8,
                },
                cx,
            )
            .unwrap();
        let old = original.energy();
        let eigenstrain_work = original.eigenstrain_work_j();
        assert!(eigenstrain_work > 0.0);
        let accepted = original.bind_harmonic_clock(Dd::ZERO).unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[], causal_budget(), cx)
            .unwrap();
        let h = 0.2;
        let next = observer
            .advance_window(&accepted, &[], &[0.0; 24], h, cx)
            .unwrap();
        let mean = samples
            .iter()
            .map(|s| {
                let f = s.free_strain_mandel();
                (f[0] + f[1] + f[2]) / 3.0
            })
            .sum::<f64>()
            / 8.0;
        let mut heat = 0.0;
        for (point, sample) in next.field.points.iter().zip(&samples) {
            let raw = sample.free_strain_mandel();
            let local_mean = (raw[0] + raw[1] + raw[2]) / 3.0;
            let projected: [f64; 6] = core::array::from_fn(|i| {
                if i < 3 {
                    raw[i] + mean - local_mean
                } else {
                    raw[i]
                }
            });
            for i in 0..6 {
                near(point.bias_mandel[i], -projected[i], 0.005);
                near(point.compatible_bias_mandel.unwrap()[i], -raw[i], 0.005);
            }
            near(point.volume_m3, 0.125, 1.0);
            // Independent real original spring solution for the accepted ramp,
            // followed by natural decay at fixed total/free displacement.
            let local_heat = law
                .branches()
                .iter()
                .map(|b| {
                    let gain = b.tau_s / ramp * -(-ramp / b.tau_s).exp_m1();
                    let q = projected.map(|v| -v * gain);
                    0.5 * dot(&q, &apply(&b.c_pa, &q)) * -(-2.0 * h / b.tau_s).exp_m1()
                })
                .sum::<f64>();
            near(point.material.maxwell_heat_j_m3, local_heat, local_heat);
            heat += local_heat / 8.0;
        }
        near(next.field.material_heat_j, heat, heat);
        near(next.field.stored_energy_change_j, -heat, heat);
        near(next.field.stress_work_j, 0.0, old.stored_j);
        assert_eq!(next.field.kelvin_heat_j, 0.0);
        near(
            next.history.history().energy().stored_j,
            old.stored_j - heat,
            old.stored_j,
        );
        assert_eq!(
            next.history.history().thermal_strains(),
            Some(samples.as_slice())
        );
        assert_eq!(
            next.history.history().eigenstrain_work_j(),
            eigenstrain_work
        );
        assert_eq!(accepted.history().energy(), old);
        assert_eq!(
            accepted.history().thermal_strains(),
            Some(samples.as_slice())
        );
        let (mut resumed, _) = next.history.into_parts();
        let held = resumed.step_internal_force(&[0.0; 24], 0.07, cx).unwrap();
        assert_eq!(held.step_eigenstrain_work_j, 0.0);
        assert_eq!(resumed.eigenstrain_work_j(), eigenstrain_work);
        assert_eq!(resumed.thermal_strains(), Some(samples.as_slice()));
        assert!(held.step_dissipated_j > 0.0);
    });
}
