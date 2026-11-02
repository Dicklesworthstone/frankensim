//! Actual continuous mechanics/filter response to physical attack and release.
use super::*;

fn gain(sample: f64) -> f64 {
    if sample < 20.0 {
        sample * 5.0
    } else if sample < 40.0 {
        (40.0 - sample) * 5.0
    } else {
        0.0
    }
}

#[test]
fn affine_attack_release_follows_continuous_filter_and_preserves_mechanical_bits() {
    let physical = mode(1900.0, 0.03);
    let mut observed = observer(
        vec![physical],
        ModalAcousticTimeBudget::audible_reference(),
        2,
    );
    let plan = observed
        .prepare_affine_tones(&[700.0], ModalToneBudget::realtime_reference())
        .unwrap();
    assert!(plan.is_affine());
    let mut mechanical = observed.model.clone();
    let mut workspace = ModalAcousticWorkspace::new(&mechanical);
    let mut reference = [0.0; 3];
    let mut pressures = Vec::new();
    let dt = 1.0 / f64::from(RATE);
    let mut carrier = (1.0, 0.0);
    for sample in 0..80 {
        let g0 = gain(sample as f64);
        let g1 = gain((sample + 1) as f64);
        let (c0, s0, c1, s1) = (
            [g0 * carrier.0],
            [g0 * carrier.1],
            [g1 * carrier.0],
            [g1 * carrier.1],
        );
        let actual = observed
            .step_affine_tones(&plan, &c0, &s0, &c1, &s1)
            .unwrap();
        let frame = mechanical
            .step_tones_affine_into(&plan.mechanics, &c0, &s0, &c1, &s1, &mut workspace)
            .unwrap();
        assert_eq!(observed.model.states(), mechanical.states());
        assert_eq!(observed.last_modal_frame(), frame);
        // Independent original ODE, with 128 RK4 substeps in each sample.
        // q'=v, v'=F(t)-2*zeta*w*v-w^2*q, z'=pole*z+v.
        let derivative = |fraction: f64, x: [f64; 3]| {
            let t = (sample as f64 + fraction) * dt;
            let force = (g0 + (g1 - g0) * fraction) * (700.0 * t).cos();
            [
                x[1],
                force
                    - 2.0 * physical.damping_ratio * physical.angular_frequency_rad_s * x[1]
                    - physical.angular_frequency_rad_s.powi(2) * x[0],
                FILTER_POLE * x[2] + x[1],
            ]
        };
        let h = dt / 128.0;
        for sub in 0..128 {
            let t = sub as f64 / 128.0;
            let a = derivative(t, reference);
            let b = derivative(
                t + 0.5 / 128.0,
                core::array::from_fn(|i| reference[i] + 0.5 * h * a[i]),
            );
            let c = derivative(
                t + 0.5 / 128.0,
                core::array::from_fn(|i| reference[i] + 0.5 * h * b[i]),
            );
            let d = derivative(
                t + 1.0 / 128.0,
                core::array::from_fn(|i| reference[i] + h * c[i]),
            );
            for i in 0..3 {
                reference[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
        pressures.push(DIRECT * reference[1] + FILTER_RESIDUE * reference[2]);
        let expected = if sample < 2 {
            0.0
        } else {
            pressures[sample - 2]
        };
        close(actual, expected);
        close(
            observed.model.states()[0].displacement_m_sqrt_kg,
            reference[0],
        );
        close(
            observed.model.states()[0].velocity_m_sqrt_kg_per_s,
            reference[1],
        );
        close(observed.filters[0].state[0], reference[2]);
        if sample == 0 {
            let endpoint_hold = reference[1] * (det::exp(FILTER_POLE * dt) - 1.0) / FILTER_POLE;
            assert!((endpoint_hold - reference[2]).abs() > 1e-9);
        }
        carrier = plan.advance_force_phase(0, carrier.0, carrier.1).unwrap();
    }
    // Turning off the load leaves the original resonant/filter state alive.
    assert!(observed.model.states()[0].velocity_m_sqrt_kg_per_s.abs() > 1e-5);
}

#[test]
fn affine_refusals_preserve_warm_mechanics_filter_delay_and_frame() {
    for mechanical_refusal in [false, true] {
        let budget = ModalAcousticTimeBudget {
            maximum_abs_pressure_pa: if mechanical_refusal { 1e6 } else { 1e-4 },
            maximum_abs_velocity_m_sqrt_kg_per_s: if mechanical_refusal { 0.001 } else { 1e6 },
            ..ModalAcousticTimeBudget::audible_reference()
        };
        let modes = vec![mode(1200.0, 0.02), mode(1700.0, 0.04)];
        let mut observed = observer(
            modes.clone(),
            budget,
            if mechanical_refusal { 3 } else { 0 },
        );
        let mut reference = observer(modes, budget, if mechanical_refusal { 3 } else { 0 });
        let plan = observed
            .prepare_affine_tones(&[740.3], ModalToneBudget::realtime_reference())
            .unwrap();
        let other = reference
            .prepare_affine_tones(&[740.3], ModalToneBudget::realtime_reference())
            .unwrap();
        for _ in 0..8 {
            observed
                .step_affine_tones(&plan, &[0.001; 2], &[0.0; 2], &[0.002; 2], &[0.0; 2])
                .unwrap();
            reference
                .step_affine_tones(&other, &[0.001; 2], &[0.0; 2], &[0.002; 2], &[0.0; 2])
                .unwrap();
        }
        let before = (
            observed.model.states().to_vec(),
            observed.last_modal_frame().clone(),
            observed
                .filters
                .iter()
                .map(|f| f.state.clone())
                .collect::<Vec<_>>(),
            observed.delay.history.clone(),
            observed.delay.cursor,
        );
        assert!(
            observed
                .step_affine_tones(&plan, &[0.0; 2], &[0.0; 2], &[0.01, 1e6], &[0.0; 2])
                .is_err()
        );
        assert!(
            observed
                .step_affine_tones(&other, &[0.0; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2])
                .is_err()
        );
        assert!(observed.step_tones(&plan, &[0.0; 2], &[0.0; 2]).is_err());
        assert!(
            observed
                .step_affine_tones(&plan, &[0.0; 2], &[0.0; 2], &[f64::NAN; 2], &[0.0; 2])
                .is_err()
        );
        let after = (
            observed.model.states().to_vec(),
            observed.last_modal_frame().clone(),
            observed
                .filters
                .iter()
                .map(|f| f.state.clone())
                .collect::<Vec<_>>(),
            observed.delay.history.clone(),
            observed.delay.cursor,
        );
        assert_eq!(before, after);
        for _ in 0..8 {
            assert_eq!(
                observed
                    .step_affine_tones(&plan, &[0.001; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2])
                    .unwrap(),
                reference
                    .step_affine_tones(&other, &[0.001; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2])
                    .unwrap()
            );
        }
    }
}

#[test]
fn affine_dc_matches_existing_continuous_ramp_and_real_feedthrough() {
    let mut affine = observer(
        vec![mode(1200.0, 0.02)],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    let mut ramp = observer(
        vec![mode(1200.0, 0.02)],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    let plan = affine
        .prepare_affine_tones(&[0.0], ModalToneBudget::realtime_reference())
        .unwrap();
    ramp.prepare_linear_forcing().unwrap();
    for sample in 0..32 {
        let f0 = sample as f64 / 32.0;
        let f1 = (sample + 1) as f64 / 32.0;
        close(
            affine
                .step_affine_tones(&plan, &[f0], &[9.0], &[f1], &[-7.0])
                .unwrap(),
            ramp.step_linear(&[f0], &[f1]).unwrap(),
        );
    }
    let mut direct = observer(
        vec![mode(0.0, 0.0)],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    direct.constant = Some(vec![C64::from_re(-2.1)]);
    direct.filters.clear();
    let plan = direct
        .prepare_affine_tones(&[600.0], ModalToneBudget::realtime_reference())
        .unwrap();
    let pressure = direct
        .step_affine_tones(&plan, &[0.0], &[0.0], &[2.0], &[0.7])
        .unwrap();
    assert_eq!(
        pressure,
        -2.1 * direct.model.states()[0].velocity_m_sqrt_kg_per_s
    );
    direct.constant = Some(vec![C64::new(-2.1, 0.1)]);
    assert!(
        direct
            .prepare_affine_tones(&[600.0], ModalToneBudget::realtime_reference())
            .is_err()
    );
}
