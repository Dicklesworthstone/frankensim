//! Retained original solid memory against independent real Maxwell histories.

mod held_thermal;

use super::*;
use fs_math::dd::Dd;

fn causal_budget() -> ViscoHarmonicMaterialBudget {
    ViscoHarmonicMaterialBudget {
        maximum_points: 8,
        maximum_payload_bytes: 64 << 20,
        maximum_work: 10_000_000_000,
    }
}
fn real_affine(nodes: &[[f64; 3]], strain: [f64; 6]) -> Vec<f64> {
    affine(nodes, strain.map(C64::from_re))
        .iter()
        .map(|z| z.re)
        .collect()
}

// Original Maxwell ODE solution: the forced particular strain plus the
// exponentially decaying difference from its ACTUAL initial spring strain.
// These real functions do not call the harmonic material or window owners.
fn spring(tones: &[(f64, [C64; 6])], tau: f64, time: f64) -> [f64; 6] {
    let mut q = [0.0; 6];
    for &(omega, e) in tones {
        let x = omega * tau;
        let (a, b) = (x * x / (1.0 + x * x), x / (1.0 + x * x));
        let (sin, cos) = (omega * time).sin_cos();
        for i in 0..6 {
            q[i] += (a * e[i].re + b * e[i].im) * cos + (a * e[i].im - b * e[i].re) * sin;
        }
    }
    q
}
fn causal_instantaneous(
    law: &TensorProny,
    tones: &[(f64, [C64; 6])],
    bias: &[f64; 6],
    initial_q: &[[f64; 6]],
    time: f64,
) -> ([f64; 4], Vec<[f64; 6]>, [f64; 6]) {
    let mut e = *bias;
    let mut rate = [0.0; 6];
    for &(omega, strain) in tones {
        let (sin, cos) = (omega * time).sin_cos();
        for i in 0..6 {
            e[i] += strain[i].re * cos + strain[i].im * sin;
            rate[i] += omega * (-strain[i].re * sin + strain[i].im * cos);
        }
    }
    let mut stress = apply(law.c_eq(), &e);
    let mut stored = 0.5 * dot(&e, &stress);
    let viscous = apply(law.kv_pa_s().unwrap(), &rate);
    let kelvin = dot(&rate, &viscous);
    for i in 0..6 {
        stress[i] += viscous[i];
    }
    let mut maxwell = 0.0;
    let mut memories = Vec::new();
    for (branch, q0) in law.branches().iter().zip(initial_q) {
        let particular = spring(tones, branch.tau_s, time);
        let start = spring(tones, branch.tau_s, 0.0);
        let q = core::array::from_fn(|i| {
            particular[i] + (q0[i] - start[i]) * (-time / branch.tau_s).exp()
        });
        let s = apply(&branch.c_pa, &q);
        let quadratic = dot(&q, &s);
        stored += 0.5 * quadratic;
        maxwell += quadratic / branch.tau_s;
        for i in 0..6 {
            stress[i] += s[i];
        }
        memories.push(q);
    }
    ([stored, kelvin, maxwell, dot(&stress, &rate)], memories, e)
}
fn causal_integral(
    law: &TensorProny,
    tones: &[(f64, [C64; 6])],
    bias: &[f64; 6],
    initial_q: &[[f64; 6]],
    duration: f64,
) -> [f64; 3] {
    let n = 4096;
    let step = duration / n as f64;
    let mut totals = [0.0; 3];
    for j in 0..=n {
        let values = causal_instantaneous(law, tones, bias, initial_q, j as f64 * step).0;
        let weight = if j == 0 || j == n {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for i in 0..3 {
            totals[i] += weight * values[i + 1];
        }
    }
    totals.map(|v| v * step / 3.0)
}

#[test]
fn original_tet_memory_and_prescribed_supports_match_real_ode_work_heat_and_continuation() {
    with_cx(|cx, _| {
        let law = law().with_band((0.0, 100.0)).unwrap();
        let material = TetViscoMaterial {
            law: law.clone(),
            density_kg_m3: 60.0,
            material_state_identity: [81; 32],
        };
        let source = ViscoElements::new(
            tet_mesh(),
            TetViscoField::Uniform(&material),
            &[0, 1, 2, 4, 5, 8],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 1000,
            },
            cx,
        )
        .unwrap();
        assert!(source.equilibrium().free_dofs.len() < 12);
        let tones = [
            (7.0, strain()),
            (7.18, strain().map(|z| C64::new(-0.3 * z.re, 0.4 * z.im))),
        ];
        let bias = core::array::from_fn(|i| 0.0002 * (i as f64 - 2.0));
        let e0 = core::array::from_fn(|i| bias[i] + tones[0].1[i].re + tones[1].1[i].re);
        let u0 = real_affine(&TET, e0);
        let mut state = source.state(cx).unwrap();
        let ramp_time = 0.2;
        state.step_internal_force(&u0, ramp_time, cx).unwrap();
        let original_energy = state.energy();
        let initial_q: Vec<_> = law
            .branches()
            .iter()
            .map(|b| e0.map(|e| e * (b.tau_s / ramp_time) * -(-ramp_time / b.tau_s).exp_m1()))
            .collect();
        let accepted = state.bind_harmonic_clock(Dd::ZERO).unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[7.0, 7.18], causal_budget(), cx)
            .unwrap();
        let banks: Vec<_> = tones.iter().map(|(_, e)| affine(&TET, *e)).collect();
        let input: Vec<_> = banks.iter().map(Vec::as_slice).collect();
        let u_bias = real_affine(&TET, bias);
        let whole = observer
            .advance_window(&accepted, &input, &u_bias, 0.375, cx)
            .unwrap();
        let field = &whole.field;
        let start = causal_instantaneous(&law, &tones, &bias, &initial_q, 0.0).0[0] / 6.0;
        let (last, q_end, e_end) = causal_instantaneous(&law, &tones, &bias, &initial_q, 0.375);
        let end = last[0] / 6.0;
        let [kelvin, maxwell, work] =
            causal_integral(&law, &tones, &bias, &initial_q, 0.375).map(|v| v / 6.0);
        let scale = start.max(end).max(work.abs());
        near(field.start_stored_energy_j, start, scale);
        near(field.end_stored_energy_j, end, scale);
        near(field.stored_energy_change_j, end - start, scale);
        near(field.kelvin_heat_j, kelvin, kelvin);
        near(field.maxwell_heat_j, maxwell, maxwell);
        near(field.stress_work_j, work, scale);
        near(field.energy_balance_residual_j, 0.0, scale);
        near(
            whole.history.history().energy().work_j,
            original_energy.work_j + work,
            scale,
        );
        near(
            whole.history.history().energy().dissipated_j,
            original_energy.dissipated_j + kelvin + maxwell,
            scale,
        );
        assert_eq!(accepted.history().energy(), original_energy);
        assert_eq!(accepted.history().displacement_m(), u0);

        let first = observer
            .advance_window(&accepted, &input, &u_bias, 0.125, cx)
            .unwrap();
        let second = observer
            .advance_window(&first.history, &input, &u_bias, 0.25, cx)
            .unwrap();
        assert_eq!(second.history.time_s_dd(), whole.history.time_s_dd());
        near(
            first.field.material_heat_j + second.field.material_heat_j,
            field.material_heat_j,
            scale,
        );
        near(
            first.field.stress_work_j + second.field.stress_work_j,
            field.stress_work_j,
            scale,
        );
        near(second.field.end_stored_energy_j, end, scale);

        // The returned object resumes the ORIGINAL affine history owner with
        // its accepted branch strains, rather than recreating a settled state.
        let (mut continued, time) = whole.history.into_parts();
        assert_eq!(time, Dd::from_f64(0.375));
        let endpoint = continued.displacement_m().to_vec();
        let held = continued.step_internal_force(&endpoint, 0.17, cx).unwrap();
        let elastic = 0.5 * dot(&e_end, &apply(law.c_eq(), &e_end));
        let held_storage = law
            .branches()
            .iter()
            .zip(q_end)
            .fold(elastic, |sum, (b, q)| {
                let q = q.map(|v| v * (-0.17 / b.tau_s).exp());
                sum + 0.5 * dot(&q, &apply(&b.c_pa, &q))
            })
            / 6.0;
        near(held.energy.stored_j, held_storage, scale);
        near(held.step_work_j, 0.0, scale);
        near(held.step_dissipated_j, end - held_storage, scale);
    });
}

#[test]
fn zero_ac_relaxation_preserves_memory_and_small_steps_on_large_dd_clocks() {
    with_cx(|cx, _| {
        let law = law(); // Zero strain bias requires no fictitious DC band query.
        let source = solid(tet_mesh(), law.clone(), cx);
        let e = strain().map(|z| z.re);
        let mut state = source.state(cx).unwrap();
        state
            .step_internal_force(&real_affine(&TET, e), 0.1, cx)
            .unwrap();
        state.step_internal_force(&[0.0; 12], 0.1, cx).unwrap();
        let start = Dd::from_f64(2.0_f64.powi(80));
        let old_energy = state.energy();
        let accepted = state.bind_harmonic_clock(start).unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[], causal_budget(), cx)
            .unwrap();
        let dt = 2.0_f64.powi(-10);
        let next = observer
            .advance_window(&accepted, &[], &[0.0; 12], dt, cx)
            .unwrap();
        assert_eq!(next.history.time_s_dd().hi, start.hi);
        assert_eq!(next.history.time_s_dd().lo, dt);
        assert_eq!(next.field.start_time_s_dd, start);
        assert_eq!(next.field.end_time_s_dd, start + Dd::from_f64(dt));
        assert!(next.field.points[0].strain_mandel.is_empty());
        assert_eq!(next.field.kelvin_heat_j, 0.0);
        assert_eq!(next.field.stress_work_j, 0.0);
        let heat = law
            .branches()
            .iter()
            .map(|b| {
                let a = (-0.1 / b.tau_s).exp();
                let gain = (b.tau_s / 0.1) * (1.0 - a);
                let q = e.map(|v| v * gain * (a - 1.0));
                0.5 * dot(&q, &apply(&b.c_pa, &q)) * -(-2.0 * dt / b.tau_s).exp_m1()
            })
            .sum::<f64>()
            / 6.0;
        near(next.field.maxwell_heat_j, heat, heat);
        near(next.field.stored_energy_change_j, -heat, heat);
        near(
            next.history.history().energy().stored_j,
            old_energy.stored_j - heat,
            old_energy.stored_j,
        );
        let second = observer
            .advance_window(&next.history, &[], &[0.0; 12], 0.1, cx)
            .unwrap();
        assert_eq!(
            second.history.time_s_dd(),
            start + Dd::from_f64(dt) + Dd::from_f64(0.1)
        );
        assert_eq!(accepted.history().energy(), old_energy);
    });
}

#[test]
fn original_bbar_points_carry_distinct_causal_memories_and_raw_waveform_limits() {
    with_cx(|cx, _| {
        let law = law();
        let source = solid(hex_mesh(), law.clone(), cx);
        let accepted = source
            .state(cx)
            .unwrap()
            .bind_harmonic_clock(Dd::ZERO)
            .unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[7.0], causal_budget(), cx)
            .unwrap();
        let amplitude = 0.004;
        let u: Vec<_> = HEX
            .iter()
            .flat_map(|p| [C64::new(0.0, amplitude * p[0] * p[1]), C64::ZERO, C64::ZERO])
            .collect();
        let next = observer
            .advance_window(&accepted, &[&u], &[0.0; 24], 0.2, cx)
            .unwrap();
        assert_eq!(next.field.points.len(), 8);
        let g = 1.0 / 3.0_f64.sqrt();
        let zero_q = vec![[0.0; 6]; law.branches().len()];
        let mut expected = [0.0; 3];
        for x in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
            for y in [(1.0 - g) * 0.5, (1.0 + g) * 0.5] {
                let e = bbar_strain(x, y, amplitude, 0.0)
                    .0
                    .map(|z| C64::new(0.0, z.re));
                let local = causal_integral(&law, &[(7.0, e)], &[0.0; 6], &zero_q, 0.2);
                for i in 0..3 {
                    expected[i] += 0.25 * local[i];
                }
            }
        }
        near(next.field.kelvin_heat_j, expected[0], expected[0]);
        near(next.field.maxwell_heat_j, expected[1], expected[1]);
        near(next.field.stress_work_j, expected[2], expected[2]);
        assert!(
            next.field
                .points
                .iter()
                .all(|p| p.compatible_strain_mandel.is_some())
        );
        let excessive: Vec<_> = HEX
            .iter()
            .flat_map(|p| [C64::new(0.0, 0.0215 * p[0] * p[1]), C64::ZERO, C64::ZERO])
            .collect();
        assert!(
            observer
                .advance_window(&accepted, &[&excessive], &[0.0; 24], 0.2, cx)
                .is_err()
        );
        assert_eq!(accepted.history().energy(), ViscoEnergy::default());
    });
}

#[test]
fn identity_continuity_budgets_and_cancellation_refuse_without_changing_accepted_history() {
    with_cx(|cx, gate| {
        let source = solid(tet_mesh(), law(), cx);
        let other = solid(tet_mesh(), law(), cx);
        let accepted = source
            .state(cx)
            .unwrap()
            .bind_harmonic_clock(Dd::ZERO)
            .unwrap();
        let foreign = other
            .state(cx)
            .unwrap()
            .bind_harmonic_clock(Dd::ZERO)
            .unwrap();
        let observer = source
            .prepare_harmonic_material_causal_observer(&[7.0], causal_budget(), cx)
            .unwrap();
        let u = affine(&TET, strain().map(|z| C64::new(0.0, z.im)));
        let good = observer
            .advance_window(&accepted, &[&u], &[0.0; 12], 0.1, cx)
            .unwrap();
        assert!(
            observer
                .advance_window(&foreign, &[&u], &[0.0; 12], 0.1, cx)
                .is_err()
        );
        assert!(
            observer
                .advance_window(&accepted, &[&u[..11]], &[0.0; 12], 0.1, cx)
                .is_err()
        );
        assert!(
            observer
                .advance_window(&accepted, &[&u], &[0.0; 11], 0.1, cx)
                .is_err()
        );
        for duration in [0.0, -0.1, f64::INFINITY, f64::NAN] {
            assert!(
                observer
                    .advance_window(&accepted, &[&u], &[0.0; 12], duration, cx)
                    .is_err()
            );
        }
        let discontinuous = affine(&TET, strain());
        assert!(
            observer
                .advance_window(&accepted, &[&discontinuous], &[0.0; 12], 0.1, cx)
                .is_err()
        );
        for denied in [
            ViscoHarmonicMaterialBudget {
                maximum_points: 0,
                ..causal_budget()
            },
            ViscoHarmonicMaterialBudget {
                maximum_payload_bytes: 1,
                ..causal_budget()
            },
            ViscoHarmonicMaterialBudget {
                maximum_work: 1,
                ..causal_budget()
            },
        ] {
            assert!(
                source
                    .prepare_harmonic_material_causal_observer(&[7.0], denied, cx)
                    .is_err()
            );
        }
        gate.request();
        assert!(matches!(
            observer.advance_window(&accepted, &[&u], &[0.0; 12], 0.1, cx),
            Err(ViscoHarmonicError::Element(ViscoElementError::Elastic(
                fs_solid::linear3::TetElasticError::Cancelled
            )))
        ));
        assert_eq!(accepted.history().energy(), ViscoEnergy::default());
        assert_eq!(accepted.history().displacement_m(), [0.0; 12]);
        with_cx(|fresh, _| {
            assert_eq!(
                good.field,
                observer
                    .advance_window(&accepted, &[&u], &[0.0; 12], 0.1, fresh)
                    .unwrap()
                    .field
            );
        });
    });
}
