//! Independent coupled equilibrium, kinetics and spatial displacement limits.
use super::super::trapping::{TrappingMembraneDiffusionState, TrappingMembraneSolve};
use super::*;

const AREA: f64 = 0.01;
const L: f64 = 2e-4;
const D: f64 = 2e-10;
const S: f64 = 1e-5;
const P: f64 = 80000.0;
const BYTES: usize = 32 * 1024 * 1024;
fn mobile(n: usize, q: f64, sealed: bool, solubility: f64) -> MembraneDiffusionState {
    let layer = MembraneDiffusionLayer {
        thickness_m: L,
        cells: n,
        diffusivity: MembraneDiffusivity::try_new(D).unwrap(),
        solubility: MolarHenrySolubility::try_new(solubility).unwrap(),
    };
    let seal = sealed.then(|| MolarPermeance::try_new(0.0).unwrap());
    MembraneDiffusion::new_with_surfaces(AREA, &[layer], &[], [seal; 2], || false)
        .unwrap()
        .initial_state(&vec![q; n], || false)
        .unwrap()
}
fn channel(
    id: &str,
    n: usize,
    q: f64,
    sealed: bool,
    capacity: f64,
    affinity: f64,
    rate: f64,
    loading: f64,
) -> CompetitiveMembraneChannel {
    CompetitiveMembraneChannel {
        species: SpeciesId::new(id).unwrap(),
        mobile: mobile(n, q, sealed, S),
        traps: vec![Some(
            LangmuirMembraneTrap::try_new(capacity, affinity, rate).unwrap(),
        )],
        initial_trapped_moles: vec![loading; n],
    }
}
fn controls() -> CompetitiveMembraneSolve {
    CompetitiveMembraneSolve {
        balance: MembraneDiffusionSolve {
            absolute_amount_tolerance_mol: 1e-19,
            relative_amount_tolerance: 1e-11,
        },
        ..Default::default()
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}
fn sealed() -> CompetitiveMembraneDiffusionState {
    CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", 1, P, true, 2.0, 2e-5, 0.5, 0.2e-6),
            channel("B", 1, 0.5 * P, true, 2.0, 5e-5, 0.2, 0.1e-6),
        ],
        || false,
    )
    .unwrap()
}

#[test]
fn finite_bath_competition_matches_an_independent_common_vacancy_root() {
    let state = sealed();
    let saved = state.clone();
    let dt = 10.0;
    let mut work = state.workspace(BYTES).unwrap();
    let candidate = state
        .step_candidate(&[[0.0; 2]; 2], dt, controls(), &mut work, || false)
        .unwrap();
    let total = [
        state.species_amount_mol(0).unwrap(),
        state.species_amount_mol(1).unwrap(),
    ];
    let old = [0.2e-6, 0.1e-6];
    let alpha = [0.5 * 2e-5, 0.2 * 5e-5];
    let kd = [0.5, 0.2];
    let cm = AREA * L * S;
    let capacity = AREA * L * 2.0;
    let n_at = |v: f64, i: usize| {
        (old[i] + (dt * alpha[i] * v / cm) * total[i]) / (1.0 + dt * kd[i] + dt * alpha[i] * v / cm)
    };
    let (mut lo, mut hi) = (0.0, capacity);
    for _ in 0..160 {
        let v = 0.5 * (lo + hi);
        if v + n_at(v, 0) + n_at(v, 1) > capacity {
            hi = v;
        } else {
            lo = v;
        }
    }
    let vacancy = 0.5 * (lo + hi);
    for i in 0..2 {
        let expected_n = n_at(vacancy, i);
        let expected_q = (total[i] - expected_n) / cm;
        near(
            candidate.state().trapped_moles(i).unwrap()[0],
            expected_n,
            1e-17,
        );
        near(
            candidate.state().mobile(i).unwrap().partial_pressures_pa()[0],
            expected_q,
            1e-6,
        );
        near(
            candidate.state().species_amount_mol(i).unwrap(),
            total[i],
            1e-18,
        );
        assert_eq!(
            candidate.step().species[i].1.boundary_transfers_mol,
            [0.0; 2]
        );
    }
    near(candidate.step().vacancy_after_mol[0], vacancy, 1e-17);
    assert!(candidate.step().linear_sweeps > candidate.step().nonlinear_iterations);
    assert_eq!(state, saved);
}

#[test]
fn scalar_and_nonadsorbing_limits_reuse_original_native_transport() {
    let input = channel("A", 8, 0.25 * P, false, 4.0, 2e-5, 0.2, 1e-8);
    let original = TrappingMembraneDiffusionState::new(
        input.mobile.clone(),
        &input.traps,
        &input.initial_trapped_moles,
        || false,
    )
    .unwrap();
    let state = CompetitiveMembraneDiffusionState::new(vec![input], || false).unwrap();
    let mut work = state.workspace(BYTES).unwrap();
    let mut original_work = original.workspace().unwrap();
    let got = state
        .step_candidate(&[[P, 0.0]], 2.0, controls(), &mut work, || false)
        .unwrap();
    let expected = original
        .step_candidate(
            [P, 0.0],
            2.0,
            TrappingMembraneSolve {
                balance: controls().balance,
                ..Default::default()
            },
            &mut original_work,
            || false,
        )
        .unwrap();
    for (&a, &b) in got
        .state()
        .mobile(0)
        .unwrap()
        .partial_pressures_pa()
        .iter()
        .zip(expected.state().mobile().partial_pressures_pa())
    {
        near(a, b, 1e-6);
    }
    for (&a, &b) in got
        .state()
        .trapped_moles(0)
        .unwrap()
        .iter()
        .zip(expected.state().trapped_moles())
    {
        near(a, b, 1e-17);
    }
    let initial = [
        mobile(8, 0.25 * P, false, S),
        mobile(8, 0.5 * P, false, 3.0 * S),
    ];
    let inputs = initial
        .iter()
        .enumerate()
        .map(|(i, m)| CompetitiveMembraneChannel {
            species: SpeciesId::new(if i == 0 { "A" } else { "B" }).unwrap(),
            mobile: m.clone(),
            traps: vec![None],
            initial_trapped_moles: vec![0.0; 8],
        })
        .collect();
    let state = CompetitiveMembraneDiffusionState::new(inputs, || false).unwrap();
    let mut work = state.workspace(BYTES).unwrap();
    let boundaries = [[P, 0.0], [0.0, P]];
    let got = state
        .step_candidate(&boundaries, 2.0, controls(), &mut work, || false)
        .unwrap();
    for i in 0..2 {
        let mut w = MembraneDiffusionWorkspace::new(8).unwrap();
        let expected = initial[i]
            .step_candidate(boundaries[i], 2.0, controls().balance, &mut w, || false)
            .unwrap();
        for (&a, &b) in got
            .state()
            .mobile(i)
            .unwrap()
            .partial_pressures_pa()
            .iter()
            .zip(expected.state().partial_pressures_pa())
        {
            near(a, b, 1e-7);
        }
        assert_eq!(got.state().trapped_moles(i).unwrap(), &[0.0; 8]);
    }
}

#[test]
fn shared_equilibrium_and_species_permutation_have_no_false_diffusion() {
    let n = 7;
    let cap = 2.0 * AREA * L / n as f64;
    let inputs = vec![
        channel("B", n, P, false, 2.0, 2.0 / P, 0.5, 0.5 * cap),
        channel("A", n, P, false, 2.0, 1.0 / P, 0.5, 0.25 * cap),
    ];
    let state = CompetitiveMembraneDiffusionState::new(inputs.clone(), || false).unwrap();
    let mut reverse = inputs;
    reverse.reverse();
    assert_eq!(
        state,
        CompetitiveMembraneDiffusionState::new(reverse, || false).unwrap()
    );
    let mut work = state.workspace(BYTES).unwrap();
    let result = state
        .step_candidate(&[[P; 2]; 2], 0.5, controls(), &mut work, || false)
        .unwrap();
    for s in 0..2 {
        assert_eq!(
            result.state().mobile(s).unwrap().partial_pressures_pa(),
            &vec![P; n]
        );
        assert_eq!(
            result.step().species[s]
                .1
                .transport
                .face_molar_fluxes_mol_m2_s,
            vec![0.0; n + 1]
        );
        assert!(
            result.step().species[s]
                .1
                .gross_adsorption_mol
                .iter()
                .all(|x| *x > 0.0)
        );
    }
    assert_eq!(result.step().nonlinear_iterations, 0);
}

#[test]
fn invading_strong_species_displaces_retained_weak_species_through_original_faces() {
    let n = 16;
    let cap = 4.0 * AREA * L / n as f64;
    let weak_p = P * 0.1;
    let mut state = CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", n, weak_p, false, 4.0, 4.0 / weak_p, 0.2, 0.8 * cap),
            channel("B", n, 0.0, false, 4.0, 20.0 / P, 0.2, 0.0),
        ],
        || false,
    )
    .unwrap();
    let initial_weak = state.species_amount_mol(0).unwrap();
    let mut weak_out = 0.0;
    let mut strong_in = 0.0;
    let mut work = state.workspace(BYTES).unwrap();
    for _ in 0..80 {
        let result = state
            .step_candidate(&[[weak_p; 2], [P, 0.0]], 0.5, controls(), &mut work, || {
                false
            })
            .unwrap();
        let a = &result.step().species[0].1;
        let b = &result.step().species[1].1;
        weak_out += a.boundary_transfers_mol[1] - a.boundary_transfers_mol[0];
        strong_in += b.boundary_transfers_mol[0] - b.boundary_transfers_mol[1];
        for i in 0..n {
            assert!(
                result.state().trapped_moles(0).unwrap()[i]
                    + result.state().trapped_moles(1).unwrap()[i]
                    <= cap * (1.0 + 1e-14)
            );
        }
        state = result.into_parts().0;
    }
    assert!(weak_out > 1e-8);
    assert!(strong_in > 1e-7);
    assert!(state.trapped_moles(0).unwrap()[0] < 0.5 * cap);
    assert!(state.trapped_moles(1).unwrap()[0] > 0.3 * cap);
    near(
        state.species_amount_mol(0).unwrap() + weak_out,
        initial_weak,
        1e-16,
    );
    near(state.species_amount_mol(1).unwrap(), strong_in, 1e-16);
}

#[test]
fn finite_reservoirs_relax_to_shared_competitive_equilibrium_conserving_each_species() {
    let membrane = CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", 6, 0.0, false, 2.0, 2e-5, 0.5, 0.0),
            channel("B", 6, 0.0, false, 2.0, 5e-5, 0.2, 0.0),
        ],
        || false,
    )
    .unwrap();
    let cm = AREA * L * S;
    let bath = MolarPressureCapacity::try_new(cm).unwrap();
    let mut state = FiniteCompetitiveMembraneDiffusionState::new(
        membrane,
        [bath; 2],
        vec![[P, 0.0], [0.0, 2.0 * P]],
    )
    .unwrap();
    let totals = [
        state.species_amount_mol(0).unwrap(),
        state.species_amount_mol(1).unwrap(),
    ];
    let affinity = [2e-5, 5e-5];
    let ct = 2.0 * AREA * L;
    let q = |v: f64, i: usize| totals[i] / (3.0 * cm + ct * affinity[i] * v);
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..160 {
        let v = 0.5 * (lo + hi);
        if v * (1.0 + affinity[0] * q(v, 0) + affinity[1] * q(v, 1)) > 1.0 {
            hi = v;
        } else {
            lo = v;
        }
    }
    let vacancy_fraction = 0.5 * (lo + hi);
    let mut work = state.workspace(BYTES).unwrap();
    for _ in 0..160 {
        state = state
            .step_candidate(20.0, controls(), &mut work, || false)
            .unwrap()
            .into_parts()
            .0;
    }
    for s in 0..2 {
        near(state.species_amount_mol(s).unwrap(), totals[s], 2e-16);
        for pressure in state.reservoir_partial_pressures_pa()[s].into_iter().chain(
            state
                .membrane()
                .mobile(s)
                .unwrap()
                .partial_pressures_pa()
                .iter()
                .copied(),
        ) {
            near(pressure, q(vacancy_fraction, s), 0.02);
        }
    }
}

#[test]
fn competitive_time_refinement_matches_independent_riccati_total_loading() {
    let cm = AREA * L * S;
    let ct = 2.0 * AREA * L;
    let initial_total = cm * P * 1.5;
    let alpha = 0.5 * 2e-5;
    let a = alpha / cm;
    let b = -(a * (initial_total + ct) + 0.5);
    let c = a * initial_total * ct;
    let root = (b * b - 4.0 * a * c).sqrt();
    let low = 2.0 * c / (-b + root);
    let high = (-b + root) / (2.0 * a);
    let r = (low / high) * (a * (low - high) * 2.0).exp();
    let exact = (low - r * high) / (1.0 - r);
    let run = |steps: usize| {
        let mut state = CompetitiveMembraneDiffusionState::new(
            vec![
                channel("A", 1, P, true, 2.0, 2e-5, 0.5, 0.0),
                channel("B", 1, 0.5 * P, true, 2.0, 2e-5, 0.5, 0.0),
            ],
            || false,
        )
        .unwrap();
        let mut work = state.workspace(BYTES).unwrap();
        for _ in 0..steps {
            state = state
                .step_candidate(
                    &[[0.0; 2]; 2],
                    2.0 / steps as f64,
                    controls(),
                    &mut work,
                    || false,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        (state.trapped_moles(0).unwrap()[0] + state.trapped_moles(1).unwrap()[0] - exact).abs()
    };
    let ratio = run(64) / run(128);
    assert!(ratio > 1.85 && ratio < 2.1, "{ratio}");
}

#[test]
fn capacity_work_and_late_cancellation_refuse_without_publishing_any_species() {
    let state = sealed();
    let saved = state.clone();
    assert!(state.workspace(16).is_err());
    let mut work = state.workspace(BYTES).unwrap();
    let limited = CompetitiveMembraneSolve {
        maximum_linear_sweeps: 1,
        ..controls()
    };
    assert!(matches!(
        state.step_candidate(&[[0.0; 2]; 2], 10.0, limited, &mut work, || false),
        Err(CompetitiveMembraneError::LinearWork { .. })
    ));
    let mut count = 0;
    let reference = state
        .step_candidate(&[[0.0; 2]; 2], 10.0, controls(), &mut work, || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        assert!(matches!(
            state.step_candidate(&[[0.0; 2]; 2], 10.0, controls(), &mut work, || {
                calls += 1;
                calls == stop
            }),
            Err(CompetitiveMembraneError::Cancelled)
        ));
        assert_eq!(state, saved);
    }
    assert_eq!(
        reference,
        state
            .step_candidate(&[[0.0; 2]; 2], 10.0, controls(), &mut work, || false)
            .unwrap()
    );
    let overfilled = vec![
        channel("A", 1, P, true, 2.0, 2e-5, 0.5, 3e-6),
        channel("B", 1, P, true, 2.0, 2e-5, 0.5, 3e-6),
    ];
    assert!(CompetitiveMembraneDiffusionState::new(overfilled, || false).is_err());
    let unequal = vec![
        channel("A", 1, P, true, 2.0, 2e-5, 0.5, 0.0),
        channel("B", 1, P, true, 3.0, 2e-5, 0.5, 0.0),
    ];
    assert!(CompetitiveMembraneDiffusionState::new(unequal, || false).is_err());
}

#[test]
fn sub_ulp_competitive_uptake_cannot_hide_behind_stored_mass_conservation() {
    let mut inert = channel("B", 1, 0.0, true, 1e18, 1.0, 0.0, 0.0);
    inert.traps[0] = None;
    let state = CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", 1, 500.0, true, 1e18, 0.004, 1e-18, 1e12),
            inert,
        ],
        || false,
    )
    .unwrap();
    let mut work = state.workspace(BYTES).unwrap();
    assert!(
        matches!(state.step_candidate(&[[0.0;2];2],1.0,controls(),&mut work,||false),
        Err(CompetitiveMembraneError::Diffusion(MembraneDiffusionError::AmountBalance {residual_mol,tolerance_mol}))
            if residual_mol.abs()>1e-7&&tolerance_mol<1e-18)
    );
}
