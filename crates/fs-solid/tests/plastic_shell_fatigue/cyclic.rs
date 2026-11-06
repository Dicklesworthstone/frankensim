//! Native nonlinear cyclic hardening survives shell reduction and fatigue.

use super::*;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, ResolvedChabocheBranches, resolve_chaboche_branches,
};
use fs_material::plastic::reduction::UniaxialVoce;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY, resolve_voce_j2,
};

fn cyclic_card<const N: usize>(branch_scale: f64) -> (MaterialCard, [BackstressProperties; N]) {
    let properties = core::array::from_fn(|i| {
        BackstressProperties::new(
            format!("backstress_{i}_modulus"),
            format!("backstress_{i}_recovery"),
        )
    });
    let mut values = vec![
        (
            VOCE_SATURATION_INCREMENT_PROPERTY.to_string(),
            Pressure::DIMS,
            0.8,
        ),
        (VOCE_SATURATION_RATE_PROPERTY.to_string(), Dims::NONE, 120.0),
    ];
    for (i, names) in properties.iter().enumerate() {
        values.push((
            names.modulus.clone(),
            Pressure::DIMS,
            branch_scale * 80.0 / (i + 1) as f64,
        ));
        values.push((names.recovery.clone(), Dims::NONE, 25.0 * (i + 1) as f64));
    }
    let extra: Vec<_> = values
        .iter()
        .map(|(name, dims, value)| (name.as_str(), *dims, *value))
        .collect();
    (card_with_properties(0, temperature(), &extra), properties)
}

fn cyclic_law<const N: usize>(
    card: &MaterialCard,
    properties: [BackstressProperties; N],
) -> ResolvedChabocheBranches<N> {
    resolve_chaboche_branches(
        card,
        &point(temperature()),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
        properties,
    )
    .unwrap()
}

#[test]
fn voce_shell_morrow_matches_exact_free_sided_cyclic_history() {
    with_cx(&CancelGate::new(), |cx| {
        let (card, _) = cyclic_card::<0>(1.0);
        let law = resolve_voce_j2(
            &card,
            &point(temperature()),
            MaterialPropertySelection::SingleClaimOnly,
            HardeningPropertyBasis::MixedPlasticStrainModuli,
            0.05,
        )
        .unwrap();
        // Borrow the original source owner; each point still retains full state.
        let section = section(&law);
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let factors = [0.006, -0.006, 0.006, -0.006, 0.006];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut g = gauge(&law, 1, StrainLifeMeanStress::Morrow);
        g.curve = sourced_curve(&card);
        let mut direct =
            PairedStrainLifeMiner::new(g.curve.clone(), g.mean_stress, limits().rainflow).unwrap();
        direct
            .push(StrainStressSample {
                mechanical_strain: 0.0,
                normal_stress_pa: 0.0,
            })
            .unwrap();
        let mut run = ShellFatiguePath::new(&p, vec![g], config(), cx).unwrap();
        let axial = UniaxialVoce::new(law.law().clone()).unwrap();
        let mut history = axial.initial_state();
        for strain in factors {
            let response = axial.response(strain, &history);
            run.advance_next(cx).unwrap();
            let sample = run.samples()[0];
            for component in 0..6 {
                close(
                    sample.total_mechanical_strain[component],
                    response.strain[component],
                    2e-10,
                );
                close(
                    sample.stress_pa[component],
                    if component == 0 {
                        response.stress_pa
                    } else {
                        0.0
                    },
                    2e-8,
                );
            }
            for state in &run.solution().section_states {
                for point in &state.points {
                    let expected = response.state.tensor_history();
                    for component in 0..6 {
                        close(
                            point.history.plastic_strain[component],
                            expected.plastic_strain[component],
                            2e-10,
                        );
                    }
                    close(point.history.alpha, expected.alpha, 2e-10);
                }
            }
            direct
                .push(StrainStressSample {
                    mechanical_strain: strain,
                    normal_stress_pa: response.stress_pa,
                })
                .unwrap();
            history = response.state;
        }
        assert!(history.equivalent_plastic_strain() > 0.015);
        assert_eq!(run.material_source(0), Some(law.material()));
        let native = p.solve(cx).unwrap();
        assert_eq!(run.solution().section_states, native.section_states);
        close(
            run.finish(cx).unwrap()[0].summary.miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-10,
        );
    });
}

fn check_chaboche_cycle<const N: usize>(cx: &Cx<'_>) {
    let (card, properties) = cyclic_card::<N>(1.0);
    let law = cyclic_law(&card, properties);
    let section = section(law.clone());
    let mesh = mesh();
    let fixed = membrane_prescribed();
    let factors = [0.006, -0.006, 0.006, -0.006, 0.006];
    let p = problem(&mesh, &section, &fixed, &[], &factors);
    let mut g = gauge(&law, 2, StrainLifeMeanStress::SmithWatsonTopper);
    g.curve = sourced_curve(&card);
    let mut direct =
        PairedStrainLifeMiner::new(g.curve.clone(), g.mean_stress, limits().rainflow).unwrap();
    direct
        .push(StrainStressSample {
            mechanical_strain: 0.0,
            normal_stress_pa: 0.0,
        })
        .unwrap();
    let mut run = ShellFatiguePath::new(&p, vec![g.clone()], config(), cx).unwrap();
    let axial = UniaxialChaboche::new(law.law().clone()).unwrap();
    let mut history = axial.initial_state();
    for strain in factors {
        let response = axial.response(strain, &history).unwrap();
        run.advance_next(cx).unwrap();
        let sample = run.samples()[0];
        for component in 0..6 {
            close(
                sample.total_mechanical_strain[component],
                response.strain[component],
                2e-10,
            );
            close(
                sample.stress_pa[component],
                if component == 0 {
                    response.stress_pa
                } else {
                    0.0
                },
                2e-8,
            );
        }
        for state in &run.solution().section_states {
            for point in &state.points {
                let expected = response.state.tensor_history();
                close(point.history.plastic.alpha, expected.plastic.alpha, 2e-10);
                for component in 0..6 {
                    close(
                        point.history.plastic.plastic_strain[component],
                        expected.plastic.plastic_strain[component],
                        2e-10,
                    );
                    for branch in 0..N {
                        close(
                            point.history.backstress_pa[branch][component],
                            expected.backstress_pa[branch][component],
                            2e-8,
                        );
                    }
                }
            }
        }
        let translation: f64 = response.state.backstress_pa().iter().sum();
        assert!(translation.abs() > 0.05);
        // Fatigue sees physical stress, never stress minus yield translation.
        assert!((sample.stress_pa[0] - (response.stress_pa - translation)).abs() > 0.05);
        direct
            .push(StrainStressSample {
                mechanical_strain: strain,
                normal_stress_pa: response.stress_pa,
            })
            .unwrap();
        history = response.state;
    }
    assert_eq!(run.material_source(0), Some(law.material()));
    let native = p.solve(cx).unwrap();
    assert_eq!(run.solution().section_states, native.section_states);
    close(
        run.finish(cx).unwrap()[0].summary.miner_damage,
        direct.finish().unwrap().miner_damage,
        1e-10,
    );
    // Matching labels, modulus and temperature do not erase branch provenance.
    let (changed_card, _) = cyclic_card::<N>(1.1);
    g.curve = sourced_curve(&changed_card);
    assert!(matches!(
        ShellFatiguePath::new(&p, vec![g], config(), cx),
        Err(ShellFatigueError::SourceMismatch { .. })
    ));
}

#[test]
fn one_two_and_eight_backstress_shells_preserve_native_tensor_histories_and_swt() {
    with_cx(&CancelGate::new(), |cx| {
        check_chaboche_cycle::<1>(cx);
        check_chaboche_cycle::<2>(cx);
        check_chaboche_cycle::<8>(cx);
    });
}

#[test]
fn chaboche_dynamic_continuation_and_late_gauge_refusal_preserve_every_branch() {
    with_cx(&CancelGate::new(), |cx| {
        let (card, properties) = cyclic_card::<2>(1.0);
        let law = cyclic_law(&card, properties);
        let section = section(law.clone());
        let mesh = mesh();
        let fixed = elastic_field(&mesh, 0.006, 0.0);
        let factors = [1.0, -1.0, 1.0, -1.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1.0), cx)
            .unwrap();
        let wide = gauge(&law, 0, StrainLifeMeanStress::Uncorrected);
        let mut narrow = gauge(&law, 2, StrainLifeMeanStress::Uncorrected);
        narrow.curve = curve(1000.0);
        let mut refused =
            ShellFatiguePath::new(&p, vec![wide.clone(), narrow], config(), cx).unwrap();
        let mut run = ShellFatiguePath::new(&p, vec![wide], config(), cx).unwrap();
        refused.advance_next(cx).unwrap();
        run.advance_next(cx).unwrap();
        let saved = refused.clone();
        let dt = 0.1;
        assert!(matches!(
            refused.advance_next_backward_euler(&reference, dt, cx),
            Err(ShellFatigueError::Gauge { gauge: 1, .. })
        ));
        assert_eq!(refused.solution(), saved.solution());
        assert_eq!(refused.samples(), saved.samples());
        assert_eq!(refused.velocity(), saved.velocity());
        assert_eq!(refused.elapsed_time_s(), saved.elapsed_time_s());
        assert_eq!(refused.next_step(), saved.next_step());
        for gauge in 0..2 {
            assert_eq!(refused.plane_summaries(gauge), saved.plane_summaries(gauge));
        }
        let mut native = problem(&mesh, &section, &fixed, &[], &[1.0])
            .solve(cx)
            .unwrap();
        let mut velocity = vec![[0.0; 6]; 3];
        assert_eq!(run.solution(), &native);
        assert!(
            native.section_states[0].points[0]
                .history
                .backstress_pa
                .iter()
                .all(|branch| branch[0].abs() > 0.01)
        );
        let mut resumed = None;
        for (index, factor) in factors.iter().enumerate().skip(1) {
            let prescribed: Vec<_> = fixed
                .iter()
                .map(|&(dof, value)| (dof, factor * value))
                .collect();
            let endpoint = problem(&mesh, &section, &prescribed, &[], &[1.0]);
            let mut step = endpoint
                .advance_backward_euler(Some(&native), &velocity, &reference, dt, cx)
                .unwrap();
            let energy = run
                .advance_next_backward_euler(&reference, dt, cx)
                .unwrap()
                .unwrap();
            assert_eq!(energy, step.inertia);
            step.solution.increments[0].load_factor = *factor;
            assert_eq!(run.solution(), &step.solution);
            assert_eq!(run.velocity(), step.velocity);
            native = step.solution;
            velocity = step.velocity;
            if index == 1 {
                let before = run.clone();
                run.preview_finish(cx).unwrap();
                assert_eq!(run.solution(), before.solution());
                resumed = Some(before);
            } else if let Some(restarted) = &mut resumed {
                assert_eq!(
                    restarted
                        .advance_next_backward_euler(&reference, dt, cx)
                        .unwrap(),
                    Some(energy)
                );
                assert_eq!(restarted.solution(), run.solution());
                assert_eq!(restarted.samples(), run.samples());
            }
        }
        close(run.elapsed_time_s(), 3.0 * dt, 1e-15);
        assert_eq!(
            run.finish(cx).unwrap(),
            resumed.unwrap().finish(cx).unwrap()
        );
    });
}
