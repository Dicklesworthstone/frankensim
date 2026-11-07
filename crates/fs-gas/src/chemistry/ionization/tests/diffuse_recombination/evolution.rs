//! Independent evolving escape/thermal references and transactional boundaries.
use super::super::{CHI, EV, NA, R, id};
use super::{SIGMA, fixture, parcel};
use crate::chemistry::collisional_radiative::{
    CollisionalRadiativeSolve, CollisionalRadiativeThermal,
};
use crate::chemistry::diffuse_recombination::DiffuseRecombinationError;
use crate::chemistry::diffuse_recombination::evolution::*;
use crate::chemistry::energy::ElectronHeavyCollisionSources;
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Power, Temperature, Time, Volume};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
    MaxwellianRateControlsV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;

const CV: f64 = 1.5 * R;

fn controls(relative: f64, duration: f64) -> DiffuseRecombinationEvolutionSolve {
    let mut local =
        CollisionalRadiativeSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    local.relative_tolerance = relative;
    local.initial_step_s = duration.max(1e-8) / 8.0;
    local.maximum_step_s = duration.max(1e-8);
    local.energy_tolerance_j = 1e-6;
    let mut solve = DiffuseRecombinationEvolutionSolve::new(local);
    solve.radiation.maximum_cells = 4;
    solve
}
fn area() -> Area {
    Area::new(SIGMA * NA / 2.0)
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:.16e} != {expected:.16e}; difference={:.16e}",
        actual - expected
    );
}
fn audit(result: &DiffuseRecombinationEvolutionResult) {
    let mut gas = 0.0;
    let mut external = 0.0;
    let mut net_rad = 0.0;
    let mut net_extent = 0.0;
    for cell in result.cells() {
        let start = cell.initial_state();
        let end = cell.state();
        let before = start.amounts()[1].value();
        let after = end.amounts()[1].value();
        let sensible = CV
            * (before
                * (end.electron_temperature().value() - start.electron_temperature().value())
                + (after - before) * end.electron_temperature().value());
        let chemical = CHI * (after - before);
        let heavy = CV * (end.heavy_temperature().value() - start.heavy_temperature().value());
        gas += sensible + chemical + heavy;
        let b = cell.balance();
        external += b.electron_deposited_energy_j + b.heavy_deposited_energy_j;
        net_rad += cell.radiation().net_gas_energy_j;
        net_extent += cell.radiation().net_ionization_extent_mol;
        assert!(b.energy_residual_j.abs() <= result.controls().chemistry.energy_tolerance_j);
        assert!(
            b.electron_sensible_energy_residual_j.abs()
                <= result.controls().chemistry.energy_tolerance_j
        );
        assert!(
            b.heavy_sensible_energy_residual_j.abs()
                <= result.controls().chemistry.energy_tolerance_j
        );
        assert_eq!(end.amounts()[1], end.amounts()[2]);
        close(end.amounts()[0].value() + after, 1.0, 0.0, 2e-15);
        assert_eq!(end.thermodynamics(), start.thermodynamics());
        assert_eq!(end.volume(), start.volume());
        close(end.mass().value(), start.mass().value(), 0.0, 2e-13);
    }
    let b = result.balance();
    close(gas + b.escaped_energy_j - external, 0.0, 4e-6, 0.0);
    close(net_rad + b.escaped_energy_j, 0.0, 1e-6, 0.0);
    close(net_extent * NA + b.escaped_photons, 0.0, NA * 2e-12, 0.0);
    assert!(b.energy_residual_j.abs() <= result.controls().chemistry.energy_tolerance_j);
    assert!(b.maximum_normalized_local_error <= 1.0);
    assert!(b.maximum_normalized_spectral_estimate <= 0.25);
    assert!(b.maximum_normalized_angular_estimate <= 0.25);
    assert!(b.maximum_normalized_tail_bound <= 0.25);
    assert!(b.radiation_work <= result.controls().maximum_radiation_work);
    assert_eq!(result.interfaces()[0].right_photons, 0.0);
    assert_eq!(result.interfaces().last().unwrap().left_photons, 0.0);
}

#[test]
fn g1_live_escape_changes_shielding_and_temperature_against_independent_ode() {
    let (state, source) = parcel(0.2, 10000.0, 1.0, 0.2);
    let duration = 1e-4;
    let result = advance_diffuse_recombination(
        &[DiffuseRecombinationEvolutionCell::new(&state, &source)],
        area(),
        Time::new(duration),
        controls(1e-8, duration),
        || false,
    )
    .unwrap();
    // Independent SciPy 1.17 DOP853 (rtol2.5e-13) and Radau (rtol1e-11)
    // agree. Integrate x,Ue,gross captures,self-absorptions,captured energy,
    // escaped energy with analytic complete constant-sigma Milne moments and
    // beta(tau)=(1/2-E3(tau))/tau, tau=2*(1-x). Thus current opacity and
    // spectrum change continuously. No production rate/transfer/time code was
    // used. Finite source support is >1000 thermal energies here.
    let cell = &result.cells()[0];
    close(
        cell.state().amounts()[1].value(),
        0.135_570_823_375_338_8,
        2e-9,
        2e-6,
    );
    close(
        cell.state().electron_temperature().value(),
        11000.934_483_724_479,
        2e-4,
        3e-7,
    );
    close(
        cell.balance().radiative_recombination_extent_mol,
        0.236_198_677_370_041_83,
        1e-9,
        3e-6,
    );
    close(
        cell.radiation().self_absorbed_photons / NA,
        0.171_769_500_745_380_58,
        1e-9,
        3e-6,
    );
    close(
        result.balance().escaped_energy_j,
        90877.299_474_883_13,
        0.005,
        3e-6,
    );
    // Preferential capture warms the remaining electrons while their sensible
    // energy and the full gas energy decrease. Recycling reduces net capture.
    assert!(cell.state().electron_temperature() > state.electron_temperature());
    assert!(cell.balance().electron_sensible_energy_change_j < 0.0);
    assert!(cell.radiation().net_gas_energy_j < 0.0);
    assert!(
        cell.balance().radiative_recombination_extent_mol > 3.0 * -cell.reaction_extent().value()
    );
    assert!(result.balance().rejected_steps > 0);
    audit(&result);
}

#[test]
fn g1_nonlocal_radiation_creates_first_electrons_with_physical_birth_temperature() {
    let (hot, source) = parcel(0.2, 30000.0, 1.0, 0.2);
    let duration = 2e-5;
    let mut endpoints = Vec::new();
    for stored_te in [1000.0, 45000.0] {
        let (neutral, neutral_source) = parcel(0.0, stored_te, 1.0, 0.2);
        let result = advance_diffuse_recombination(
            &[
                DiffuseRecombinationEvolutionCell::new(&hot, &source),
                DiffuseRecombinationEvolutionCell::new(&neutral, &neutral_source),
            ],
            area(),
            Time::new(duration),
            controls(1e-8, duration),
            || false,
        )
        .unwrap();
        // Independent analytic-Milne DOP853 ODE (rtol2.5e-13). For constant
        // sigma each i->j absorption fraction is
        // [1/2-E3(tau_i)-E3(tau_j)+E3(tau_i+tau_j)]/(2*tau_i).
        // Net x and Ue use this cross-slab transfer and each beta(tau_i).
        // Zero initial electrons imply zero emission, not a seeded Te-dependent
        // emissivity. Source excess energy fixes the receiver's birth Te.
        for (cell, x, te, captures) in [
            (
                &result.cells()[0],
                0.186_292_325_547_690_15,
                30196.026_251_933_17,
                0.049_058_711_613_199_94,
            ),
            (
                &result.cells()[1],
                0.006_454_860_865_140_347,
                27337.630_048_585_754,
                1.931_047_864_523_505_3e-5,
            ),
        ] {
            close(cell.state().amounts()[1].value(), x, 1e-9, 3e-6);
            close(cell.state().electron_temperature().value(), te, 2e-4, 3e-7);
            close(
                cell.balance().radiative_recombination_extent_mol,
                captures,
                1e-11,
                5e-6,
            );
        }
        close(
            result.balance().escaped_energy_j,
            11988.581_677_307_639,
            0.005,
            3e-6,
        );
        assert!(result.cells()[1].radiation().nonlocal_absorbed_photons > 0.0);
        assert!(result.cells()[1].radiation().net_electron_energy_j > 0.0);
        assert!(result.cells()[1].reaction_extent().value() > 0.0);
        endpoints.push(result.cells()[1].state().clone());
        audit(&result);
    }
    assert_eq!(endpoints[0], endpoints[1]);
}

fn thermal() -> CollisionalRadiativeThermal {
    let momentum = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "analytic constant momentum-transfer cross section, not measured H data".into(),
            license: "analytic fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(1e-20),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(500.0 * EV),
                cross_section: Area::new(1e-20),
            },
        ],
        Area::new(1e-20),
    )
    .unwrap();
    let mut inputs =
        CollisionalRadiativeThermal::new(Temperature::new(1000.0), Temperature::new(50000.0));
    inputs.collisions = Some(ElectronHeavyCollisionSources {
        momentum_sources: vec![momentum],
        coulomb: None,
        maximum_charged_mole_fraction: 0.5,
        maximum_target_thermal_speed_ratio: 0.1,
        rate_controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    });
    inputs.electron_power = Power::new(3e6);
    inputs.heavy_power = Power::new(2e7);
    inputs
}

#[test]
fn g1_g3_diffuse_collisions_both_thermal_pools_and_work_match_ode_and_second_order() {
    let (original, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        original.amounts(),
        Volume::new(1000.0),
        original.heavy_temperature(),
        original.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    let thermal = thermal();
    let cells = [DiffuseRecombinationEvolutionCell::new(&state, &source).with_thermal(&thermal)];
    let duration = 1e-4;
    let result = advance_diffuse_recombination(
        &cells,
        area(),
        Time::new(duration),
        controls(1e-8, duration),
        || false,
    )
    .unwrap();
    // Independent DOP853/Radau agree with simultaneous analytic ramp-impact,
    // statistical three-body balance, complete Milne, E3 escape and neutral
    // exchange: km=(4/3)*sigma_m*sqrt(8*kT/pi/me),
    // G=3R*x*(1-x)*NA/V*me/(mp+me)*km. Evolve x,Ue,Uh and all integrals.
    let cell = &result.cells()[0];
    let b = cell.balance();
    close(
        cell.state().amounts()[1].value(),
        0.099_997_752_371_083_49,
        1e-10,
        1e-7,
    );
    close(
        cell.state().electron_temperature().value(),
        8751.421_682_702_401,
        0.005,
        1e-7,
    );
    close(
        cell.state().heavy_temperature().value(),
        6309.004_878_804_287,
        0.001,
        1e-7,
    );
    close(
        b.collisional_extent_mol,
        2.449_497_195_566_299_2e-6,
        1e-11,
        1e-5,
    );
    close(
        b.radiative_recombination_extent_mol,
        1.832_889_629_027_612_7e-5,
        1e-11,
        1e-5,
    );
    close(
        cell.radiation().self_absorbed_photons / NA,
        1.363_177_017_820_194e-5,
        1e-11,
        1e-5,
    );
    close(
        result.balance().escaped_energy_j,
        6.568_709_207_988_448,
        1e-5,
        1e-5,
    );
    close(
        b.electron_to_heavy_neutral_energy_j,
        1853.814_270_467_824_1,
        0.005,
        1e-6,
    );
    close(b.electron_deposited_energy_j, 300.0, 1e-10, 0.0);
    close(b.heavy_deposited_energy_j, 2000.0, 1e-10, 0.0);
    assert_eq!(b.electron_to_heavy_coulomb_energy_j, 0.0);
    audit(&result);
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut c = controls(0.1, duration);
        c.chemistry.initial_step_s = duration / steps as f64;
        c.chemistry.maximum_step_s = c.chemistry.initial_step_s;
        c.radiation.relative_tolerance = 1e-8;
        let refined =
            advance_diffuse_recombination(&cells, area(), Time::new(duration), c, || false)
                .unwrap();
        assert_eq!(refined.balance().rejected_steps, 0);
        errors.push(
            (refined.cells()[0].state().electron_temperature().value() - 8751.421_682_702_401)
                .abs(),
        );
        audit(&refined);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "refinement errors={errors:?}, ratio={ratio}"
        );
    }
}

#[test]
fn g1_opaque_self_recycling_keeps_finite_live_cooling_and_pair_loss() {
    let (state, source) = parcel(0.2, 10000.0, 1.0, 0.2);
    let duration = 1e14;
    let result = advance_diffuse_recombination(
        &[DiffuseRecombinationEvolutionCell::new(&state, &source)],
        Area::new(SIGMA * NA * 0.8 / 1e20),
        Time::new(duration),
        controls(1e-8, duration),
        || false,
    )
    .unwrap();
    // Independent analytic-Milne ODE with beta=1/(2*tau), exact to far beyond
    // binary64 at tau>=1e20. The long physical exposure resolves net losses
    // while >10^17 mol of repeated recycling cancels before state mutation.
    let cell = &result.cells()[0];
    close(
        cell.state().amounts()[1].value(),
        0.198_235_575_076_013_03,
        1e-9,
        2e-7,
    );
    close(
        cell.state().electron_temperature().value(),
        10022.119_689_317_416,
        0.0001,
        1e-7,
    );
    close(
        result.balance().escaped_energy_j,
        2480.379_230_742_697,
        0.005,
        3e-6,
    );
    close(
        cell.balance().radiative_recombination_extent_mol,
        3.532_741_342_114_400_6e17,
        0.0,
        3e-6,
    );
    assert!(cell.balance().emitted_radiation_energy_j > 1e23);
    assert!(cell.radiation().net_electron_energy_j < -100.0);
    assert!(cell.reaction_extent().value() < -0.001);
    audit(&result);
}

#[test]
fn g0_g4_live_dark_empty_zero_time_and_transactional_resource_refusals() {
    let (neutral, source) = parcel(0.0, 10000.0, 1.0, 0.2);
    let c = controls(1e-7, 1e-5);
    let dark = advance_diffuse_recombination(
        &[DiffuseRecombinationEvolutionCell::new(&neutral, &source)],
        area(),
        Time::new(1e10),
        c,
        || false,
    )
    .unwrap();
    assert_eq!(dark.cells()[0].state(), &neutral);
    assert_eq!(dark.balance().accepted_steps, 0);
    assert_eq!(dark.balance().radiation_evaluations, 0);
    let empty = advance_diffuse_recombination(&[], area(), Time::new(1e10), c, || false).unwrap();
    assert!(empty.cells().is_empty());
    assert_eq!(empty.interfaces().len(), 1);
    assert_eq!(empty.balance().escaped_photons, 0.0);
    let (state, source) = parcel(0.2, 10000.0, 1.0, 0.2);
    let original = state.clone();
    let cells = [DiffuseRecombinationEvolutionCell::new(&state, &source)];
    let zero = advance_diffuse_recombination(&cells, area(), Time::new(0.0), c, || false).unwrap();
    assert_eq!(zero.cells()[0].state(), &state);
    assert_eq!(zero.balance().escaped_energy_j, 0.0);
    let mut budget = c;
    budget.radiation.maximum_cells = usize::MAX;
    assert_eq!(
        advance_diffuse_recombination(&cells, area(), Time::new(1e-5), budget, || false),
        Err(DiffuseRecombinationEvolutionError::Invalid(
            "radiation cell limit"
        ))
    );
    budget = c;
    budget.maximum_radiation_evaluations = 1;
    assert_eq!(
        advance_diffuse_recombination(&cells, area(), Time::new(1e-5), budget, || false),
        Err(DiffuseRecombinationEvolutionError::RadiationBudget)
    );
    budget = c;
    budget.maximum_radiation_work = 1;
    assert_eq!(
        advance_diffuse_recombination(&cells, area(), Time::new(1e-5), budget, || false),
        Err(DiffuseRecombinationEvolutionError::Radiation(
            DiffuseRecombinationError::WorkBudget
        ))
    );
    budget = c;
    budget.chemistry.maximum_steps = 1;
    assert_eq!(
        advance_diffuse_recombination(&cells, area(), Time::new(1e-5), budget, || false),
        Err(DiffuseRecombinationEvolutionError::StepBudget)
    );
    let mut polls = 0;
    assert_eq!(
        advance_diffuse_recombination(&cells, area(), Time::new(1e-5), c, || {
            polls += 1;
            polls == 500
        }),
        Err(DiffuseRecombinationEvolutionError::Cancelled)
    );
    let (_, foreign) = fixture(0.2, 10000.0, 0.0, 1e6, 0.0);
    assert!(matches!(
        advance_diffuse_recombination(
            &[DiffuseRecombinationEvolutionCell::new(&state, &foreign)],
            area(),
            Time::new(1e-5),
            c,
            || false
        ),
        Err(DiffuseRecombinationEvolutionError::Cell { .. })
    ));
    assert_eq!(state, original);
}
