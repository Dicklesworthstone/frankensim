//! Continuous-time references test the estimator, not another implicit solve.
use super::*;

const AREA: f64 = 0.01;
const L: f64 = 2e-4;
const D: f64 = 2e-10;
const S: f64 = 1e-5;
const P: f64 = 8e4;
const BYTES: usize = 16 * 1024 * 1024;

fn channel(
    id: &str,
    cells: usize,
    pressure: f64,
    sealed: bool,
    trap: bool,
    loading: f64,
) -> CompetitiveMembraneChannel {
    let layer = MembraneDiffusionLayer {
        thickness_m: L,
        cells,
        diffusivity: MembraneDiffusivity::try_new(D).unwrap(),
        solubility: MolarHenrySolubility::try_new(S).unwrap(),
    };
    let seal = sealed.then(|| MolarPermeance::try_new(0.0).unwrap());
    CompetitiveMembraneChannel {
        species: SpeciesId::new(id).unwrap(),
        mobile: MembraneDiffusion::new_with_surfaces(AREA, &[layer], &[], [seal; 2], || false)
            .unwrap()
            .initial_state(&vec![pressure; cells], || false)
            .unwrap(),
        traps: vec![trap.then(|| LangmuirMembraneTrap::try_new(2.0, 2e-5, 0.5).unwrap())],
        initial_trapped_moles: vec![loading; cells],
    }
}
fn controls(step: f64, pressure: f64, amount: f64) -> CompetitiveMembraneAdvanceControls {
    CompetitiveMembraneAdvanceControls {
        temporal: MembraneDiffusionAdvanceControls {
            step: MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: 1e-19,
                relative_amount_tolerance: 1e-11,
            },
            initial_step_s: step,
            minimum_step_s: 1e-7,
            maximum_step_s: step,
            maximum_evaluations: 20000,
            absolute_pressure_tolerance_pa: pressure,
            absolute_amount_tolerance_mol: amount,
            relative_tolerance: 0.0,
            amount_balance_tolerance_mol: 1e-14,
        },
        maximum_iterations: 64,
        maximum_backtracks: 32,
        maximum_linear_sweeps: 256,
        linear_relative_tolerance: 0.05,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, allowance {tolerance}"
    );
}

#[test]
fn g1_competitive_loading_accuracy_refines_toward_exact_riccati_dynamics() {
    let state = CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", 1, P, true, true, 0.0),
            channel("B", 1, 0.5 * P, true, true, 0.0),
        ],
        || false,
    )
    .unwrap();
    let cm = AREA * L * S;
    let ct = AREA * L * 2.0;
    let total = cm * 1.5 * P;
    let a = 0.5 * 2e-5 / cm;
    let b = -(a * (total + ct) + 0.5);
    let c = a * total * ct;
    let root = (b * b - 4.0 * a * c).sqrt();
    let low = 2.0 * c / (-b + root);
    let high = (-b + root) / (2.0 * a);
    let r = (low / high) * (a * (low - high) * 2.0).exp();
    let exact = (low - r * high) / (1.0 - r);
    let run = |amount| {
        let mut work = state.adaptive_workspace(BYTES).unwrap();
        state
            .advance_candidate(
                &[[0.0; 2]; 2],
                2.0,
                controls(2.0, 1e6, amount),
                &mut work,
                || false,
            )
            .unwrap()
    };
    let coarse = run(2e-8);
    let fine = run(2e-10);
    let error = |c: &CompetitiveMembraneAdvanceCandidate| {
        (c.state().trapped_moles(0).unwrap()[0] + c.state().trapped_moles(1).unwrap()[0] - exact)
            .abs()
    };
    assert!(error(&fine) < 0.2 * error(&coarse));
    assert!(
        fine.report().species[0].1.transport.evaluations
            > coarse.report().species[0].1.transport.evaluations
    );
    assert!(fine.report().species[0].1.transport.rejected_intervals > 0);
    for (s, (_, row)) in fine.report().species.iter().enumerate() {
        near(
            fine.state().species_amount_mol(s).unwrap(),
            state.species_amount_mol(s).unwrap(),
            1e-17,
        );
        assert_eq!(row.transport.face_amounts_mol, vec![0.0; 2]);
        assert_eq!(row.transport.interval_s, [0.0, 2.0]);
        assert!(row.transport.amount_residual_mol.abs() < 1e-17);
    }
    assert_eq!(state.elapsed_s(), 0.0);
}

#[test]
fn g1_nonfirst_species_controls_finite_bath_accuracy_and_matches_exact_eigenmodes() {
    let state = CompetitiveMembraneDiffusionState::new(
        vec![
            channel("A", 1, 1.0, false, false, 0.0),
            channel("B", 1, 0.0, false, false, 0.0),
        ],
        || false,
    )
    .unwrap();
    let cm = AREA * L * S;
    let state = FiniteCompetitiveMembraneDiffusionState::new(
        state,
        [MolarPressureCapacity::try_new(cm).unwrap(); 2],
        vec![[1.0; 2], [P, 0.0]],
    )
    .unwrap();
    let duration = 200.0;
    let k = 2.0 * D / (L * L);
    let a = (-k * duration).exp();
    let b = (-3.0 * k * duration).exp();
    let exact = [
        P / 3.0 + P * a / 2.0 + P * b / 6.0,
        P * (1.0 - b) / 3.0,
        P / 3.0 - P * a / 2.0 + P * b / 6.0,
    ];
    let run = |pressure, amount| {
        state
            .advance_candidate(
                duration,
                controls(100.0, pressure, amount),
                &mut state.adaptive_workspace(BYTES).unwrap(),
                || false,
            )
            .unwrap()
    };
    let coarse = run(500.0, 1e-8);
    let fine = run(5.0, 1e-10);
    let error =
        |c: &CompetitiveMembraneAdvanceCandidate<FiniteCompetitiveMembraneDiffusionState>| {
            let q = c.state().reservoir_partial_pressures_pa()[1];
            let got = [
                q[0],
                c.state()
                    .membrane()
                    .mobile(1)
                    .unwrap()
                    .partial_pressures_pa()[0],
                q[1],
            ];
            got.into_iter()
                .zip(exact)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0_f64, f64::max)
        };
    assert!(
        error(&fine) < 0.2 * error(&coarse),
        "{} -> {}",
        error(&coarse),
        error(&fine)
    );
    assert!(fine.report().species[1].1.transport.rejected_intervals > 0);
    for s in 0..2 {
        let report = &fine.report().species[s].1.transport;
        near(
            fine.state().species_amount_mol(s).unwrap(),
            state.species_amount_mol(s).unwrap(),
            1e-17,
        );
        assert!(
            report
                .finite_reservoirs
                .as_ref()
                .unwrap()
                .closed_amount_residual_mol
                .abs()
                < 1e-17
        );
        assert_eq!(
            report.evaluations,
            fine.report().species[0].1.transport.evaluations
        );
    }
    assert_eq!(fine.state().reservoir_partial_pressures_pa()[0], [1.0; 2]);
    assert_eq!(
        fine.state()
            .membrane()
            .mobile(0)
            .unwrap()
            .partial_pressures_pa(),
        &[1.0]
    );
}

#[test]
fn g3_shared_site_displacement_keeps_species_face_balances_and_late_refusal_is_atomic() {
    let n = 4;
    let capacity = 2.0 * AREA * L / n as f64;
    let weak_pressure = 0.1 * P;
    let mut weak = channel("A", n, weak_pressure, false, true, 0.8 * capacity);
    weak.traps[0] = Some(LangmuirMembraneTrap::try_new(2.0, 4.0 / weak_pressure, 0.2).unwrap());
    let mut strong = channel("B", n, 0.0, false, true, 0.0);
    strong.traps[0] = Some(LangmuirMembraneTrap::try_new(2.0, 20.0 / P, 0.2).unwrap());
    let state = CompetitiveMembraneDiffusionState::new(vec![weak, strong], || false).unwrap();
    let saved = state.clone();
    let boundaries = [[weak_pressure; 2], [P, 0.0]];
    let c = controls(5.0, 100.0, 1e-9);
    let mut work = state.adaptive_workspace(BYTES).unwrap();
    let mut polls = 0;
    let got = state
        .advance_candidate(&boundaries, 40.0, c, &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(got.state().trapped_moles(0).unwrap()[0] < 0.5 * capacity);
    assert!(got.state().trapped_moles(1).unwrap()[0] > 0.3 * capacity);
    let a = &got.report().species[0].1.transport;
    let b = &got.report().species[1].1.transport;
    assert!(a.face_amounts_mol[0] - a.face_amounts_mol[n] < 0.0);
    assert!(b.face_amounts_mol[0] - b.face_amounts_mol[n] > 0.0);
    for (s, (_, r)) in got.report().species.iter().enumerate() {
        near(
            got.state().species_amount_mol(s).unwrap() - state.species_amount_mol(s).unwrap(),
            r.transport.face_amounts_mol[0] - r.transport.face_amounts_mol[n],
            1e-16,
        );
        assert!(
            r.transport
                .gross_face_amounts_mol
                .iter()
                .zip(&r.transport.face_amounts_mol)
                .all(|(gross, net)| *gross >= net.abs())
        );
    }
    let mut limited = controls(0.001, 1e6, 1.0);
    limited.temporal.maximum_evaluations = 3;
    assert!(
        matches!(state.advance_candidate(&boundaries,40.0,limited,&mut work,||false),
        Err(Error::WorkLimit {time_s,evaluations:3}) if time_s > 0.0)
    );
    let mut seen = 0;
    assert!(matches!(
        state.advance_candidate(&boundaries, 40.0, c, &mut work, || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(state, saved);
    assert_eq!(
        got,
        state
            .advance_candidate(&boundaries, 40.0, c, &mut work, || false)
            .unwrap()
    );
    assert!(state.adaptive_workspace(1).is_err());
}
