//! Evaluated MF8/MT458 data through the original material-temperature history.
//! Source curves are manufactured; the independent ODE tests the admitted
//! equations and first-order coupling, not an evaluated reactor or isotope.

use super::*;
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupLimits, BurnupModel};
use crate::endf::fission_energy::temperature::tests::assembled_fixture;

const EV_J: f64 = 1.602_176_634e-19;
// A small manufactured thermal mass makes feedback visible with O(100)
// expected nuclei. The original HOT caloric law still owns temperature.
const CAPACITY_J_K: f64 = 2.5 * EV_J;
const MASS_KG: f64 = CAPACITY_J_K / 4.0;

fn source_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-27,
        },
        max_composition_feedback_terms: 0,
    }
}

fn step_controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-27,
    }
}

fn system(shutdown: Option<f64>, maximum_k: f64) -> (DepositionSimulation<'static>, BurnupHeating) {
    let mut data = assembled_fixture();
    data.temperature_max_k = maximum_k;
    let groups = Arc::clone(&data.groups);
    let inventory = BurnupInventory::new(
        Arc::new(BurnupModel::new(data, source_budget().source.limits).unwrap()),
        vec![100.0, 0.0, 0.0, 0.0],
    )
    .unwrap();
    let changes = shutdown
        .into_iter()
        .map(|time_s| NeutronFluxChange {
            time_s,
            flux: GroupFlux::new(Arc::clone(&groups), vec![0.0]).unwrap(),
        })
        .collect();
    let source = BurnupCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory,
        fluxes: NeutronFluxSchedule::new(GroupFlux::new(groups, vec![2.0]).unwrap(), changes)
            .unwrap(),
        composition_feedback: vec![],
    };
    let initial = thermal(
        vec![cell(1, 1, MASS_KG, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = BurnupHeating::new(&initial, vec![source], source_budget()).unwrap();
    (initial, heating)
}

// Independent continuous equations, without calling the table, collapse or
// nuclear propagator. Coordinates: parent, radioactive product, stable product,
// stable daughter, temperature, deposited/radiation/escaped eV, fissions.
fn derivative(y: &[f64; 9]) -> [f64; 9] {
    let x = (y[4] - 300.0) / 600.0;
    let fissions = 2.0 * (0.1 + 0.2 * x) * y[0];
    let production = 2.0 * (0.1 + 0.25 * x) * y[0];
    let stable_production = 2.0 * (0.1 + 0.15 * x) * y[0];
    let decay = 0.3 * y[1];
    let deposited = 2.0 * (0.6 + 1.4 * x) * y[0] + 4.0 * decay;
    [
        -fissions,
        production - decay,
        stable_production,
        decay,
        deposited / 2.5,
        deposited,
        // The fixture routes both 2-eV prompt neutrons and 1-eV prompt
        // photons to radiation; only the daughter has an escaped channel.
        3.0 * fissions + 3.0 * decay,
        6.0 * decay,
        fissions,
    ]
}

fn ode_reference(duration: f64) -> [f64; 9] {
    let mut y = [100.0, 0.0, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 0.0];
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let a = derivative(&y);
        let b = derivative(&std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(&std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(&std::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..9 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn g1_evaluated_fission_heats_the_original_material_and_refines_temperature_feedback() {
    let (initial, heating) = system(None, 900.0);
    let reference = ode_reference(2.0);
    let mut errors = Vec::new();
    // The explicit material/source split is first order. Resolve its physical
    // feedback error before comparing to the continuous ODE's energy target.
    for dt in [0.02, 0.01, 0.005] {
        let (heated, evolved, audit) = heating
            .advance_to(&initial, 2.0, step_controls(dt), || Ok(()))
            .unwrap();
        let populations = evolved.source(CellId(1)).unwrap().inventory.populations();
        let temperature = heated.temperature_k(CellId(1)).unwrap();
        let released = evolved.audit().released;
        assert!(temperature > 350.0 && temperature < 900.0);
        close(evolved.audit().fissions, 100.0 - populations[0], 2e-10);
        close(
            populations.iter().sum(),
            100.0 + evolved.audit().fissions,
            3e-10,
        );
        close(
            temperature,
            300.0 + released.deposited_j / CAPACITY_J_K,
            2e-9,
        );
        close(
            heated.energy_j() - initial.energy_j(),
            released.deposited_j,
            1e-27,
        );
        close(audit.external_heat_j, released.deposited_j, 1e-27);
        assert_eq!(heated.cells()[&CellId(1)].mass_kg, MASS_KG);
        assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 0.0);
        let error = populations
            .iter()
            .zip(&reference[..4])
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            + (temperature - reference[4]).abs() / 100.0;
        errors.push(error);
        // The independently integrated outward accounts remain distinct from
        // material heating; their convergence follows the same physical step.
        if dt == 0.005 {
            close(released.deposited_j / EV_J, reference[5], 0.25);
            close(released.radiation_j / EV_J, reference[6], 0.15);
            close(released.escaped_j / EV_J, reference[7], 0.15);
        }
    }
    assert!(
        errors
            .windows(2)
            .all(|p| (1.7..2.3).contains(&(p[0] / p[1]))),
        "first-order feedback errors {errors:?}"
    );
    assert!(errors[2] < 0.08, "feedback errors {errors:?}");
}

#[test]
fn g1_shutdown_and_restart_preserve_temperature_series_daughter_heat_once() {
    let (initial, heating) = system(Some(0.5), 900.0);
    let (at_shutdown, stopped, _) = heating
        .advance_to(&initial, 0.5, step_controls(0.05), || Ok(()))
        .unwrap();
    let before = stopped
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .to_vec();
    let lifetime_before = stopped
        .source(CellId(1))
        .unwrap()
        .inventory
        .cumulative_released();
    let (cooled, final_state, _) = stopped
        .advance_to(&at_shutdown, 2.0, step_controls(0.05), || Ok(()))
        .unwrap();
    let after = final_state
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations();
    let remaining = before[1] * (-0.3_f64 * 1.5).exp();
    let decays = before[1] - remaining;
    close(after[0], before[0], 1e-12);
    close(after[1], remaining, 2e-10);
    close(after[2], before[2], 1e-12);
    close(after[3] - before[3], decays, 2e-10);
    // NuclearHeating::audit describes the last advance, whereas the original
    // inventory keeps lifetime counters. Shutdown must have no new fissions.
    close(final_state.audit().fissions, 0.0, 1e-12);
    close(
        final_state.audit().released.deposited_j,
        4.0 * EV_J * decays,
        1e-27,
    );
    close(
        final_state.audit().released.radiation_j,
        3.0 * EV_J * decays,
        1e-27,
    );
    close(
        final_state.audit().released.escaped_j,
        6.0 * EV_J * decays,
        1e-27,
    );
    close(
        cooled.temperature_k(CellId(1)).unwrap() - at_shutdown.temperature_k(CellId(1)).unwrap(),
        4.0 * EV_J * decays / CAPACITY_J_K,
        2e-9,
    );
    let lifetime_after = final_state
        .source(CellId(1))
        .unwrap()
        .inventory
        .cumulative_released();
    close(
        lifetime_after.deposited_j - lifetime_before.deposited_j,
        4.0 * EV_J * decays,
        1e-27,
    );
    let (direct_heat, direct, _) = heating
        .advance_to(&initial, 2.0, step_controls(0.05), || Ok(()))
        .unwrap();
    for (&actual, &expected) in after
        .iter()
        .zip(direct.source(CellId(1)).unwrap().inventory.populations())
    {
        close(actual, expected, 2e-10);
    }
    close(cooled.energy_j(), direct_heat.energy_j(), 1e-27);
}

#[test]
fn g4_temperature_domain_and_cancellation_discard_the_whole_coupled_trial() {
    let (initial, limited) = system(None, 301.0);
    let energy = initial.energy_j();
    assert!(
        limited
            .advance_to(&initial, 0.1, step_controls(0.1), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.energy_j(), energy);
    assert_eq!(
        limited.source(CellId(1)).unwrap().inventory.populations(),
        &[100.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(limited.audit().fissions, 0.0);

    let (initial, heating) = system(None, 900.0);
    let mut total = 0;
    let expected = heating
        .advance_to(&initial, 0.2, step_controls(0.05), || {
            total += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, total / 2, total] {
        let mut calls = 0;
        assert!(matches!(
            heating.advance_to(&initial, 0.2, step_controls(0.05), || {
                calls += 1;
                if calls == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(FissionThermalError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(heating.audit().fissions, 0.0);
        assert_eq!(
            heating.source(CellId(1)).unwrap().inventory.populations(),
            &[100.0, 0.0, 0.0, 0.0]
        );
    }
    let retried = heating
        .advance_to(&initial, 0.2, step_controls(0.05), || Ok(()))
        .unwrap();
    assert_eq!(
        retried.0.energy_j().to_bits(),
        expected.0.energy_j().to_bits()
    );
    assert_eq!(
        retried.1.source(CellId(1)).unwrap().inventory.populations(),
        expected
            .1
            .source(CellId(1))
            .unwrap()
            .inventory
            .populations()
    );
}
