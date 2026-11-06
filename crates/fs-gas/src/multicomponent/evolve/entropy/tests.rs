//! Native conservative motion and original-source entropy, without molecular fluxes.

use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn controls() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    }
}

fn field(
    gas: MixtureGas,
    count: usize,
    boundaries: Boundaries,
    initial: impl Fn(usize) -> (Vec<f64>, f64, f64),
) -> MixtureEuler1d {
    let states = (0..count).map(|i| {
        let (rho, velocity, pressure) = initial(i);
        gas.conserved(&rho, velocity, pressure).unwrap()
    }).collect();
    MixtureEuler1d::new(
        gas, Length::new(1.0 / count as f64), 0.7, states, boundaries, || false,
    ).unwrap()
}

fn temperature(flow: &MixtureEuler1d, i: usize) -> f64 {
    flow.gas().primitive(flow.cell(i).unwrap(), || false).unwrap().temperature
}

// Rows are original [Cp/R constant, Cp/R linear coefficient, M/R].
// Integrate the source calorics and ideal partial pressures independently.
fn source_entropy(flow: &MixtureEuler1d, rows: &[[f64; 3]]) -> f64 {
    (0..flow.cell_count()).map(|i| {
        let t = temperature(flow, i);
        let state = flow.cell(i).unwrap();
        flow.mesh().volumes()[i] * rows.iter().enumerate().filter_map(|(k, row)| {
            let rho = state.partial_densities[k];
            (rho > 0.0).then(|| {
                let [a, b, mass] = *row;
                (rho / mass) * ((a - 1.0) * t.ln() + b * t - (rho / mass).ln())
            })
        }).sum::<f64>()
    }).sum()
}

fn two_cell() -> MixtureEuler1d {
    field(
        gas(vec![source("A", R, 3.5, 0.0, 0.0)]),
        2, Boundaries::Periodic,
        |i| if i == 0 { (vec![1.0], 0.0, 1.0) } else { (vec![0.125], 0.0, 0.1) },
    )
}

#[test]
fn g1_two_cell_pressure_density_mode_matches_independent_conservative_rk2_and_entropy() {
    let mut flow = two_cell();
    let before = source_entropy(&flow, &[[3.5, 0.0, 1.0]]);
    let (report, entropy) = flow.step_entropy_consistent(
        Time::new(0.02), controls(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    assert_eq!(report.elapsed.value(), 0.02);
    assert_eq!(flow.time().value(), 0.02);
    assert_eq!(report.halvings, 0);
    assert_eq!(report.hll_fallback_faces, 0);
    assert_eq!(report.limited_cells, 0);

    // On a two-cell periodic grid at u=0 the symmetric EC momentum fluxes
    // cancel. Species and energy obey dqL/dt=-a*(qL-qR)/dx, dqR=-dqL,
    // a=max(sqrt(1.4*T)). Applying the independent SSP-RK2 Butcher formula gives:
    let expected = [
        [0.9605590146647074, 2.3961789258868715, 0.9988268442874513],
        [0.16444098533529256, 0.35100857411312836, 0.8548225999982363],
    ];
    for (i, [rho, energy, t]) in expected.into_iter().enumerate() {
        let cell = flow.cell(i).unwrap();
        assert!((cell.partial_densities[0] - rho).abs() < 3e-12);
        assert!((cell.energy - energy).abs() < 3e-12);
        assert!(cell.momentum.abs() < 3e-13);
        assert!((temperature(&flow, i) - t).abs() < 3e-12);
    }
    let change = source_entropy(&flow, &[[3.5, 0.0, 1.0]]) - before;
    assert!((change - 0.027299800511537853).abs() < 3e-12);
    assert!((entropy.source_entropy_change_j_k - change).abs() < 3e-12);
    assert!((entropy.entropy_change_j_k - change).abs() < 3e-12);
    assert!((entropy.minimum_step_entropy_change_j_k - change).abs() < 3e-12);
    assert!(entropy.maximum_relative_species_residual < 1e-12);
    assert!(report.conservation_residual.energy.abs() < 3e-13);
    assert_eq!(report.boundary_change.energy, 0.0);
    assert_eq!(report.boundary_change.mass(), 0.0);
}

fn shock(formation: f64) -> MixtureEuler1d {
    field(
        gas(vec![
            source("A", R, 3.5, 0.2, formation),
            source("B", 2.0 * R, 4.0, 0.4, -2.0 * formation),
            source("C", R, 4.0, 0.0, 3.0 * formation),
        ]),
        32,
        Boundaries::Ends { left: Boundary::Reflecting, right: Boundary::Reflecting },
        |i| if i < 16 {
            (vec![0.8, 0.2, 0.0], 0.0, 2.0)
        } else {
            (vec![0.1, 1.1, 0.0], 0.0, 0.2)
        },
    )
}

#[test]
fn g1_variable_caloric_mixture_pressure_waves_keep_closed_entropy_energy_and_source_gauges() {
    let rows = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0], [4.0, 0.0, 1.0]];
    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut flow = shock(formation);
        let initial = flow.states.clone();
        let before = source_entropy(&flow, &rows);
        let original = flow.totals();
        let mut accounted = 0.0;
        for _ in 0..6 {
            let (step, entropy) = flow.step_entropy_consistent(
                Time::new(0.002), controls(), MixtureEntropyControl::default(), || false,
            ).unwrap();
            assert!(entropy.entropy_change_j_k >= -1e-10);
            assert_eq!(step.boundary_change.mass(), 0.0);
            assert_eq!(step.boundary_change.energy, 0.0);
            assert!(step.conservation_residual.energy.abs() < 2e-11);
            assert!(step.conservation_residual.momentum.abs() < 2e-11);
            accounted += entropy.source_entropy_change_j_k;
        }
        assert_ne!(flow.states, initial);
        assert!(flow.states.chunks_exact(5).any(|q| q[3].abs() > 1e-6));
        let after = flow.totals();
        assert!((after.energy - original.energy).abs() < 2e-11);
        for k in 0..3 {
            assert!((after.species_masses[k] - original.species_masses[k]).abs() < 2e-12);
        }
        assert_eq!(after.species_masses[2], 0.0);
        let measured = source_entropy(&flow, &rows) - before;
        assert!(measured > 0.0);
        assert!((accounted - measured).abs() < 3e-10);
        for i in 0..flow.cell_count() {
            assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0);
            assert!(temperature(&flow, i) >= 0.001);
        }
        if let Some((states, time, delta)) = &reference {
            assert_eq!(&flow.states, states);
            assert_eq!(&flow.time(), time);
            assert!((accounted - delta).abs() < 3e-12);
        } else {
            reference = Some((flow.states.clone(), flow.time(), accounted));
        }
    }
}

#[test]
fn g0_uniform_motion_preserves_the_source_state_and_common_absence() {
    let mut flow = field(
        gas(vec![
            source("A", R, 3.5, 0.2, 0.0),
            source("B", 2.0 * R, 4.0, 0.4, 0.0),
            source("C", R, 4.0, 0.0, 0.0),
        ]),
        8, Boundaries::Periodic,
        |_| (vec![0.8, 0.2, 0.0], 0.7, 2.0),
    );
    let original = flow.states.clone();
    let (_, entropy) = flow.step_entropy_consistent(
        Time::new(0.01), controls(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    assert_eq!(flow.states, original);
    assert_eq!(entropy.entropy_change_j_k, 0.0);
    assert_eq!(entropy.source_entropy_change_j_k, 0.0);
}

#[test]
fn g4_native_entropy_step_cancellation_and_unsupported_open_boundaries_are_atomic() {
    let mut reference = shock(0.0);
    let mut polls = 0;
    let (_, expected) = reference.step_entropy_consistent(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        || { polls += 1; false },
    ).unwrap();
    let mut flow = shock(0.0);
    let original = flow.states.clone();
    let mut visits = 0;
    let error = flow.step_entropy_consistent(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        || { visits += 1; visits >= polls },
    ).unwrap_err();
    assert_eq!(error, Error::Cancelled);
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);
    let (_, actual) = flow.step_entropy_consistent(
        Time::new(0.002), controls(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    assert!(flow.states.iter().zip(&reference.states).all(|(a,b)| (a-b).abs() < 3e-11));
    assert!((actual.entropy_change_j_k-expected.entropy_change_j_k).abs() < 3e-11);
    let accepted = flow.states.clone();
    let accepted_time = flow.time();
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    assert!(matches!(
        flow.step_entropy_consistent(
            Time::new(0.002), controls(), MixtureEntropyControl::default(), || false,
        ),
        Err(Error::EntropyDomain(_)),
    ));
    assert_eq!(flow.states, accepted);
    assert_eq!(flow.time(), accepted_time);
}
