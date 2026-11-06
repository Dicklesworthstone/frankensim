//! Coupled moving-volume alpha distributions and independent physical accounts.

use super::*;
use crate::fusion::volume::{FusionVolumeHistory, FusionVolumeKnot};
use fs_qty::Volume;

fn history(knots: &[(f64, f64)]) -> FusionVolumeHistory {
    FusionVolumeHistory::new(
        knots
            .iter()
            .map(|&(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )
    .unwrap()
}

fn complete_history() -> FusionVolumeHistory {
    history(&[(0.0, 1.0), (0.04, 0.6), (0.06, 0.6), (0.1, 1.0)])
}

fn moving(
    initial: &KineticFusionState,
    alpha: &KineticAlphaTransport,
    thermal: &FusionNetworkTransport,
    geometry: &FusionVolumeHistory,
    c: FusionNetworkSolve,
) -> Result<KineticFusionResult, FusionNetworkError> {
    advance_kinetic_fusion_with_volume(initial, DT, heating(), thermal, alpha, geometry, c, || {
        false
    })
}

fn check_values(values: &[(&str, f64, f64)], relative: f64, absolute: f64) {
    for &(name, actual, expected) in values {
        assert!(
            (actual - expected).abs() <= absolute + relative * expected.abs(),
            "{name}: actual={actual:.17e}, expected={expected:.17e}"
        );
    }
}

fn check_profile(state: &KineticFusionState, expected: &[(usize, f64)]) {
    for &(index, amount) in expected {
        check_values(
            &[("alpha bin", state.alpha_populations_mol()[index], amount)],
            8e-4,
            2e-17,
        );
    }
}

// Independent same-grid reference: SciPy Radau at rtol=2e-12, with direct
// number/thermal-energy/birth/heat/work ODEs and exact breaks at volume knots.
// Coulomb Psi is evaluated by an independent incomplete-gamma implementation.
// Physical centers and widths are s(t) times their original values, and each
// work derivative is -(2/3)*(d ln V/dt)*U_s. There are no production finite-
// fuel maps, geometric splitting steps, or uniformization polynomials.
// Tightening reference rtol by 100 changed any energy account by <3.2e-9 J.
// Agreement tests time/coupling behavior on 128 cells, not continuum grid error.
#[test]
fn g1_compress_hold_expand_alphas_match_independent_three_store_reference() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    let geometry = complete_history();
    let result = moving(
        &initial,
        &alpha,
        &FusionNetworkTransport::default(),
        &geometry,
        controls(),
    )
    .unwrap();
    let state = result.state();
    let b = result.balance();
    let a = result.alpha_balance();
    check_values(
        &[
            (
                "ion kT",
                state.background_heavy_temperature().value() / KELVIN_PER_KEV,
                10.010009761894604,
            ),
            (
                "electron kT",
                state.electron_temperature().value() / KELVIN_PER_KEV,
                5.565815689594611,
            ),
            ("pressure", state.pressure().value(), 259914.1100612465),
        ],
        3e-6,
        1e-8,
    );
    close(
        b.reaction_extents_mol[0],
        1.1068520170675633e-7,
        4e-6,
        1e-20,
    );
    check_values(
        &[
            (
                "ion storage",
                state.thermal_heavy_energy().value(),
                240246.35054814638,
            ),
            (
                "electron storage",
                state.thermal_electron_energy().value(),
                133761.29770528624,
            ),
            (
                "alpha storage",
                state.alpha_energy().value(),
                15863.516838437132,
            ),
            ("alpha birth", a.born_alpha_energy_j, 38080.865178298816),
            ("neutrons", b.neutron_energy_escaped_j, 150460.8709928943),
            ("nuclear Q", b.nuclear_energy_released_j, 187844.34551800543),
            (
                "reactant kinetic debit",
                b.reacting_ion_energy_removed_j,
                697.3906531878583,
            ),
            (
                "alpha ion heat",
                a.deposited_heavy_energy_j,
                568.7986689481922,
            ),
            (
                "alpha electron heat",
                a.deposited_electron_energy_j,
                17027.1331198708,
            ),
            ("ion work", b.boundary_work_on_heavy_j, 48.44743238832956),
            (
                "electron work",
                b.boundary_work_on_electrons_j,
                -3429.082964585027,
            ),
            (
                "alpha work",
                a.boundary_work_on_alphas_j,
                -4621.416551042802,
            ),
        ],
        4e-5,
        3e-3,
    );
    close(distribution_moments(state).2, 410539.6162087517, 2e-4, 1e-3);
    check_profile(
        state,
        &[
            (24, 4.711388585044238e-15),
            (40, 6.511429760776018e-11),
            (48, 1.0268249796203787e-9),
            (56, 3.67563537721166e-9),
            (64, 3.6952773340455906e-9),
            (72, 2.0981219138169387e-9),
            (80, 1.2739419588472776e-9),
            (88, 9.029794965462814e-10),
            (96, 6.51947164957507e-10),
            (104, 4.78766478551778e-10),
            (107, 3.9472206792033526e-10),
            (108, 6.800475422070372e-12),
        ],
    );
    assert_eq!(b.electron_to_heavy_coulomb_energy_j, 0.0);
    assert_eq!(b.photon_energy_escaped_j, 0.0);
    assert_eq!(b.charged_product_heavy_energy_j, 0.0);
    assert_eq!(b.charged_product_electron_energy_j, 0.0);
    // A closed volume cycle need not have zero net work when births and heat
    // exchanges alter the three stores along the path.
    assert!(a.boundary_work_on_alphas_j < 0.0);
    close(
        result.total_boundary_work().value(),
        b.boundary_work_on_heavy_j + b.boundary_work_on_electrons_j + a.boundary_work_on_alphas_j,
        3e-14,
        1e-9,
    );
    assert_eq!(result.volume_history(), Some(&geometry));
    assert_eq!(state.volume(), initial.volume());
    assert_eq!(state.energy_grid(), initial.energy_grid());
    assert_eq!(state.elapsed_time(), Time::new(0.1));
    check_accounts(&result);
}

#[test]
fn g1_compressed_alpha_burn_with_exchange_and_radiation_matches_independent_ode() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    let thermal = thermal_transport(&alpha, 0.1);
    let geometry = complete_history();
    let result = moving(&initial, &alpha, &thermal, &geometry, controls()).unwrap();
    let state = result.state();
    let b = result.balance();
    let a = result.alpha_balance();
    check_values(
        &[
            (
                "ion kT",
                state.background_heavy_temperature().value() / KELVIN_PER_KEV,
                8.481601070204968,
            ),
            (
                "electron kT",
                state.electron_temperature().value() / KELVIN_PER_KEV,
                6.845280108623761,
            ),
            ("pressure", state.pressure().value(), 254821.50889110912),
        ],
        3e-6,
        1e-8,
    );
    close(b.reaction_extents_mol[0], 8.849906019784525e-8, 4e-6, 1e-20);
    check_values(
        &[
            (
                "ion storage",
                state.thermal_heavy_energy().value(),
                203618.07616059814,
            ),
            (
                "electron storage",
                state.thermal_electron_energy().value(),
                164510.21764832956,
            ),
            (
                "alpha storage",
                state.alpha_energy().value(),
                14103.969527735971,
            ),
            ("alpha birth", a.born_alpha_energy_j, 30430.93222499716),
            ("neutrons", b.neutron_energy_escaped_j, 120282.18173402172),
            ("nuclear Q", b.nuclear_energy_released_j, 150192.1466057028),
            (
                "reactant kinetic debit",
                b.reacting_ion_energy_removed_j,
                520.9673533166704,
            ),
            (
                "alpha ion heat",
                a.deposited_heavy_energy_j,
                453.7727555394035,
            ),
            (
                "alpha electron heat",
                a.deposited_electron_energy_j,
                11845.379898255715,
            ),
            (
                "thermal exchange",
                b.electron_to_heavy_coulomb_energy_j,
                -44747.71992916535,
            ),
            ("photons", b.photon_energy_escaped_j, 2214.606181540517),
            ("ion work", b.boundary_work_on_heavy_j, 8106.495587540961),
            (
                "electron work",
                b.boundary_work_on_electrons_j,
                -10031.5235475506,
            ),
            (
                "alpha work",
                a.boundary_work_on_alphas_j,
                -4027.8100434662206,
            ),
        ],
        4e-5,
        3e-3,
    );
    close(
        distribution_moments(state).2,
        319112.87983129313,
        2e-4,
        1e-3,
    );
    check_profile(
        state,
        &[
            (24, 4.145807150301631e-18),
            (40, 8.201921606356777e-13),
            (48, 6.40194312275464e-11),
            (56, 1.096911270418937e-9),
            (64, 3.8280749905827765e-9),
            (72, 3.01504016535425e-9),
            (80, 1.3335518208807623e-9),
            (88, 8.319174110121226e-10),
            (96, 5.533456756771944e-10),
            (104, 3.7365204795752494e-10),
            (107, 2.731418192502248e-10),
            (108, 2.012275516301807e-12),
        ],
    );
    assert!(b.electron_to_heavy_coulomb_energy_j < 0.0);
    assert!(b.photon_energy_escaped_j > 0.0);
    assert_eq!(result.thermal_transport(), &thermal);
    assert_eq!(result.alpha_transport(), &alpha);
    check_accounts(&result);
}

#[test]
fn g3_constant_volume_history_is_identical_to_the_rigid_kinetic_owner() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let thermal = thermal_transport(&alpha, 0.1);
    let geometry = history(&[(0.0, 1.0), (0.004, 1.0)]);
    let held = moving(&initial, &alpha, &thermal, &geometry, controls()).unwrap();
    let rigid = advance_kinetic_fusion(
        &initial,
        DT,
        heating(),
        &thermal,
        &alpha,
        geometry.duration(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(held.state(), rigid.state());
    assert_eq!(held.balance(), rigid.balance());
    assert_eq!(held.alpha_balance(), rigid.alpha_balance());
    assert_eq!(held.total_boundary_work().value(), 0.0);
    assert_eq!(held.alpha_balance().boundary_work_on_alphas_j, 0.0);
    assert_eq!(held.volume_history(), Some(&geometry));
    assert!(rigid.volume_history().is_none());
    check_accounts(&held);
}

#[test]
fn g3_restart_at_compressed_volume_retains_physical_grid_and_whole_interval_accounts() {
    let (initial, alpha) = scenario(0.0, &grid(5000.0));
    let thermal = FusionNetworkTransport::default();
    let first_geometry = history(&[(0.0, 1.0), (0.04, 0.6)]);
    let first = moving(&initial, &alpha, &thermal, &first_geometry, controls()).unwrap();
    assert_eq!(first.state().volume(), Volume::new(0.6));
    assert!(first.state().alpha_amount().value() > 0.0);
    let scale = (1.0_f64 / 0.6).powf(2.0 / 3.0);
    for (physical, original) in first
        .state()
        .energy_grid()
        .energies()
        .iter()
        .zip(initial.energy_grid().energies())
    {
        close(physical.value(), original.value() * scale, 3e-14, 0.0);
    }
    let second_geometry = history(&[(0.0, 0.6), (0.02, 0.6), (0.06, 1.0)]);
    let second = moving(
        first.state(),
        first.alpha_transport(),
        &thermal,
        &second_geometry,
        controls(),
    )
    .unwrap();
    let state = second.state();
    let b1 = first.balance();
    let b2 = second.balance();
    let a1 = first.alpha_balance();
    let a2 = second.alpha_balance();
    close(state.elapsed_time().value(), 0.1, 0.0, 2e-16);
    assert_eq!(second.elapsed_time(), Time::new(0.06));
    assert_eq!(second.initial_state(), first.state());
    assert_eq!(second.alpha_transport(), &alpha);
    for (physical, original) in state
        .energy_grid()
        .energies()
        .iter()
        .zip(initial.energy_grid().energies())
    {
        close(physical.value(), original.value(), 4e-14, 0.0);
    }
    check_values(
        &[
            (
                "restart ion kT",
                state.background_heavy_temperature().value() / KELVIN_PER_KEV,
                10.010009761894604,
            ),
            (
                "restart electron kT",
                state.electron_temperature().value() / KELVIN_PER_KEV,
                5.565815689594611,
            ),
        ],
        3e-6,
        1e-8,
    );
    close(
        b1.reaction_extents_mol[0] + b2.reaction_extents_mol[0],
        1.1068520170675633e-7,
        4e-6,
        1e-20,
    );
    check_values(
        &[
            (
                "restart alpha storage",
                state.alpha_energy().value(),
                15863.516838437132,
            ),
            (
                "restart alpha birth",
                a1.born_alpha_energy_j + a2.born_alpha_energy_j,
                38080.865178298816,
            ),
            (
                "restart alpha ion heat",
                a1.deposited_heavy_energy_j + a2.deposited_heavy_energy_j,
                568.7986689481922,
            ),
            (
                "restart alpha electron heat",
                a1.deposited_electron_energy_j + a2.deposited_electron_energy_j,
                17027.1331198708,
            ),
            (
                "restart ion work",
                b1.boundary_work_on_heavy_j + b2.boundary_work_on_heavy_j,
                48.44743238832956,
            ),
            (
                "restart electron work",
                b1.boundary_work_on_electrons_j + b2.boundary_work_on_electrons_j,
                -3429.082964585027,
            ),
            (
                "restart alpha work",
                a1.boundary_work_on_alphas_j + a2.boundary_work_on_alphas_j,
                -4621.416551042802,
            ),
        ],
        4e-5,
        3e-3,
    );
    close(distribution_moments(state).2, 410539.6162087517, 2e-4, 1e-3);
    assert_eq!(
        state.thermal_amounts()[state.plasma().indices()[4]].value(),
        0.0
    );
    check_accounts(&first);
    check_accounts(&second);
}

#[test]
fn g4_intermediate_grid_support_volume_mismatch_and_late_cancel_refuse_atomically() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let unchanged = initial.clone();
    let thermal = FusionNetworkTransport::default();
    let mismatch = history(&[(0.0, 1.1), (0.001, 1.1)]);
    assert!(matches!(
        moving(&initial, &alpha, &thermal, &mismatch, controls()),
        Err(FusionNetworkError::Invalid(_))
    ));

    let mut c = controls();
    c.minimum_step_s = 1e-7;
    c.maximum_steps = 2000;
    let capped = alpha_transport(
        TestParticleCoulombV1::new(
            alpha.collisions().background().clone(),
            initial.plasma().indices()[4],
            0.002,
        )
        .unwrap(),
    );
    let excessive_compression = history(&[(0.0, 1.0), (0.0002, 0.125), (0.0004, 1.0)]);
    // Both endpoints have the original 5 MeV grid. The interior exceeds the
    // declared nonrelativistic 7.45 MeV ceiling and cannot be skipped.
    assert_eq!(
        advance_kinetic_fusion_with_volume(
            &initial,
            FusionNetworkChannels::NONE,
            heating(),
            &thermal,
            &capped,
            &excessive_compression,
            c,
            || false,
        )
        .unwrap_err(),
        FusionNetworkError::OutsideDomain
    );

    let (empty, birth_source) = scenario(0.0, &grid(5000.0));
    let birth_source =
        KineticAlphaTransport::new(birth_source.collisions().clone(), 0.05, 0.25, 0.99, 100_000)
            .unwrap();
    let excessive_expansion = history(&[(0.0, 1.0), (0.0002, 2.0), (0.0004, 1.0)]);
    // The grid ceiling falls to 3.15 MeV inside this closed geometry cycle,
    // below a DT alpha birth. No missing birth is clipped into an edge cell.
    // Deliberately loose time tolerances permit one trial per segment: its
    // source midpoints still cover the birth, while the exact turnaround knot
    // does not. Endpoint source admission must catch this even without a
    // time-error rejection or a later unsupported midpoint.
    let mut knot_controls = c;
    knot_controls.relative_tolerance = 0.1;
    knot_controls.absolute_amount_tolerance_mol = 1e-6;
    knot_controls.absolute_temperature_tolerance_k = 1e8;
    knot_controls.energy_tolerance_j = 1e6;
    knot_controls.initial_step_s = 0.0002;
    knot_controls.minimum_step_s = 0.0002;
    knot_controls.maximum_step_s = 0.0002;
    assert_eq!(
        moving(
            &empty,
            &birth_source,
            &thermal,
            &excessive_expansion,
            knot_controls
        )
        .unwrap_err(),
        FusionNetworkError::OutsideDomain
    );
    assert!(empty.alpha_populations_mol().iter().all(|n| *n == 0.0));

    let short = history(&[(0.0, 1.0), (0.0001, 0.8), (0.0002, 1.0)]);
    let mut calls = 0;
    advance_kinetic_fusion_with_volume(
        &initial,
        FusionNetworkChannels::NONE,
        heating(),
        &thermal,
        &alpha,
        &short,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_calls = 0;
    assert_eq!(
        advance_kinetic_fusion_with_volume(
            &initial,
            FusionNetworkChannels::NONE,
            heating(),
            &thermal,
            &alpha,
            &short,
            controls(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            },
        ),
        Err(FusionNetworkError::Cancelled)
    );
    assert_eq!(initial, unchanged);
}
