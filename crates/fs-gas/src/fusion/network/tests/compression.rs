//! Independent thermodynamic and coupled-burn checks for prescribed geometry.

use super::*;
use crate::fusion::volume::FusionVolumeKnot;
use fs_qty::Volume;
use fs_thermochem::fusion::plasma::DtFusionPlasmaErrorV1;

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

fn pool_energies(state: &FrozenIonizedGasStateV1) -> [f64; 2] {
    [
        state.heavy_internal_energy().value(),
        state.electron_internal_energy().value(),
    ]
}

#[test]
fn g1_adiabatic_volume_motion_preserves_sources_populations_and_reverses() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let original = state.clone();
    let transport = FusionNetworkTransport::default();
    let no_time = history(&[(0.0, 0.025)]);
    let unchanged = advance_fusion_network_with_volume(
        &state,
        &plasma,
        FusionNetworkChannels::NONE,
        heating(),
        &transport,
        &no_time,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(unchanged.state(), &state);
    assert_eq!(unchanged.elapsed_time(), Time::new(0.0));
    assert_eq!(unchanged.volume_history(), Some(&no_time));
    assert_eq!(unchanged.balance().boundary_work_on_heavy_j, 0.0);
    assert_eq!(unchanged.balance().boundary_work_on_electrons_j, 0.0);

    // An eightfold volume compression gives T*=4 and p*=32 for either
    // monatomic pool, including unequal initial electron/ion temperatures.
    let compress = history(&[(0.0, 0.025), (0.1, 0.003125)]);
    let compressed = advance_fusion_network_with_volume(
        &state,
        &plasma,
        FusionNetworkChannels::NONE,
        heating(),
        &transport,
        &compress,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(compressed.state().amounts(), state.amounts());
    assert_eq!(compressed.state().volume(), Volume::new(0.003125));
    assert_eq!(compressed.plasma(), &plasma);
    assert_eq!(compressed.transport(), &transport);
    assert_eq!(compressed.volume_history(), Some(&compress));
    close(
        compressed.state().heavy_temperature().value(),
        4.0 * state.heavy_temperature().value(),
        2e-12,
        1e-4,
    );
    close(
        compressed.state().electron_temperature().value(),
        4.0 * state.electron_temperature().value(),
        2e-12,
        1e-4,
    );
    close(
        compressed.state().pressure().value(),
        32.0 * state.pressure().value(),
        3e-12,
        1e-6,
    );
    let before = pool_energies(&state);
    let after = pool_energies(compressed.state());
    let compression_work = [
        compressed.balance().boundary_work_on_heavy_j,
        compressed.balance().boundary_work_on_electrons_j,
    ];
    for i in 0..2 {
        close(after[i], 4.0 * before[i], 3e-12, 1e-8);
        close(compression_work[i], 3.0 * before[i], 3e-12, 1e-8);
    }
    check_accounts(&compressed);

    let expand = history(&[(0.0, 0.003125), (0.1, 0.025)]);
    let returned = advance_fusion_network_with_volume(
        compressed.state(),
        &plasma,
        FusionNetworkChannels::NONE,
        heating(),
        &transport,
        &expand,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(returned.state().amounts(), state.amounts());
    assert_eq!(returned.state().volume(), state.volume());
    let expansion_work = [
        returned.balance().boundary_work_on_heavy_j,
        returned.balance().boundary_work_on_electrons_j,
    ];
    for (i, energy) in pool_energies(returned.state()).into_iter().enumerate() {
        assert!(expansion_work[i] < 0.0);
        close(energy, before[i], 5e-12, 1e-8);
        close(expansion_work[i], -compression_work[i], 5e-12, 1e-8);
    }
    check_accounts(&returned);
    assert_eq!(state, original);

    // A millionfold expansion must retain its positive thermal survivor,
    // even though negative work nearly cancels the original thermal energy.
    // Disabled channels permit the lower original caloric support here.
    let (hot_plasma, hot) = fixture::plasma(MIXED_DENSITY, 0.025, [80.0, 20.0]);
    let mut broad = controls();
    broad.minimum_heavy_temperature_k = 1e4;
    let deep = history(&[(0.0, 0.025), (0.1, 25000.0)]);
    let diluted = advance_fusion_network_with_volume(
        &hot,
        &hot_plasma,
        FusionNetworkChannels::NONE,
        heating(),
        &transport,
        &deep,
        broad,
        || false,
    )
    .unwrap();
    assert_eq!(diluted.state().amounts(), hot.amounts());
    for (actual, initial) in pool_energies(diluted.state())
        .into_iter()
        .zip(pool_energies(&hot))
    {
        assert!(actual > 0.0);
        close(actual, 1e-4 * initial, 2e-10, 1e-8);
    }
    close(
        diluted.state().heavy_temperature().value(),
        1e-4 * hot.heavy_temperature().value(),
        2e-10,
        1e-5,
    );
    close(
        diluted.state().electron_temperature().value(),
        1e-4 * hot.electron_temperature().value(),
        2e-10,
        1e-5,
    );
    check_accounts(&diluted);
}

#[test]
fn g1_compression_with_external_powers_matches_exact_energy_and_refines_quadratically() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let path = history(&[(0.0, 0.025), (0.1, 0.0125)]);
    let held = FusionNetworkHeating {
        heavy_power: Power::new(12500.0),
        electron_power: Power::new(-8000.0),
        ..heating()
    };
    // Direct solution of dU/dt=a*U+P, a=-(2/3) dln(V)/dt:
    // U(t)=U0 exp(a*t)+P expm1(a*t)/a; W_on=U(t)-U0-P*t.
    // These independently evaluated constants do not use the source splitting.
    let expected_energy = [64535.54873398812, 17676.30332411597];
    let expected_temperature = [283290349.52529293, 69675615.67844918];
    let expected_work = [23631.677042488125, 6700.305064215969];
    let advance = |c| {
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::NONE,
            held,
            &FusionNetworkTransport::default(),
            &path,
            c,
            || false,
        )
        .unwrap()
    };
    let resolved = advance(controls());
    for (actual, expected) in pool_energies(resolved.state())
        .into_iter()
        .zip(expected_energy)
    {
        close(actual, expected, 3e-8, 1e-5);
    }
    for (actual, expected) in [
        resolved.state().heavy_temperature().value(),
        resolved.state().electron_temperature().value(),
    ]
    .into_iter()
    .zip(expected_temperature)
    {
        close(actual, expected, 3e-8, 1e-3);
    }
    let b = resolved.balance();
    close(b.boundary_work_on_heavy_j, expected_work[0], 1e-7, 1e-5);
    close(b.boundary_work_on_electrons_j, expected_work[1], 1e-7, 1e-5);
    assert_eq!(b.external_heavy_energy_j, 1250.0);
    assert_eq!(b.external_electron_energy_j, -800.0);
    assert_eq!(resolved.state().amounts(), state.amounts());
    check_accounts(&resolved);

    let mut errors = Vec::new();
    for intervals in [4, 8, 16] {
        let mut c = controls();
        c.relative_tolerance = 0.1;
        c.absolute_amount_tolerance_mol = 1e-3;
        c.absolute_temperature_tolerance_k = 1e6;
        c.energy_tolerance_j = 1e6;
        c.initial_step_s = 0.1 / f64::from(intervals);
        c.minimum_step_s = c.initial_step_s;
        c.maximum_step_s = c.initial_step_s;
        let result = advance(c);
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            pool_energies(result.state())
                .into_iter()
                .zip(expected_energy)
                .map(|(actual, expected)| (actual - expected).abs())
                .sum::<f64>(),
        );
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.8..4.2).contains(&ratio),
            "errors={errors:?}, ratio={ratio}"
        );
    }
}

// Independent direct-SI RK4 integrates extensive populations, both thermal
// energies, each product account and -(2/3) U_s dlnV/dt for each pool's work.
// It stops at each knot and evaluates densities as N_s/V(t), without production
// positive reaction maps, pressure-work maps or thermal/source routines.
// 512/1024 steps per segment agree within 5.5e-12 in temperatures and 4.7e-11
// in extents; the net electron work differs by 5.1e-7 J after cancellation.
// Sources: Bosch-Hale (1992), Eq.14/Table VII, doi:10.1088/0029-5515/32/4/I07;
// separate ion/electron PdV: SAND2004-5996, section 10, printed page 99.
const MOVING: Reference = Reference {
    amounts: [
        8.307908979405144e-7,
        1.2413743871962945e-5,
        2.8679112640692567e-6,
        1.6608905348214246e-6,
        4.5368823463202775e-7,
    ],
    extents: [
        3.8454489486396026e-8,
        4.5044599974802036e-10,
        4.223860014210533e-10,
        9.897835217005494e-11,
    ],
    temperature: [142073989.57519463, 143368404.42695433],
    products: [
        [13243.864951049769, 52288.54500417072],
        [36.95000035420831, 108.46622636616011],
        [42.393143479903706, 125.05415150522774],
        [35.58444092623409, 140.70543613199538],
    ],
    energy: [
        65742.9717721965,
        278.5915817877223,
        52397.01123053688,
        10901.938678030188,
        2722.6134454171493,
        -15872.10142259311,
        2442.9071444453816,
    ],
};

#[test]
fn g1_compress_hold_expand_matches_independent_coupled_nuclear_thermal_and_work_history() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let transport = transport(&state);
    let path = history(&[(0.0, 0.025), (0.04, 0.0125), (0.06, 0.0125), (0.1, 0.02)]);
    let mut c = controls();
    // The first proposed interval spans both interior knots; the owner must
    // respect the prescribed hold and turnaround before error refinement.
    c.initial_step_s = 0.1;
    let result = advance_fusion_network_with_volume(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        &path,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.state().volume(), Volume::new(0.02));
    assert_eq!(result.elapsed_time(), Time::new(0.1));
    assert_eq!(result.volume_history(), Some(&path));
    assert_eq!(result.plasma(), &plasma);
    assert_eq!(result.transport(), &transport);
    check_reference(&result, &MOVING);
    check_accounts(&result);
    for (actual, expected) in pool_energies(result.state())
        .into_iter()
        .zip([32296.52563536544, 36371.74037242438])
    {
        close(actual, expected, 3e-8, 1e-5);
    }
    let b = result.balance();
    close(b.boundary_work_on_heavy_j, 6070.733502829125, 3e-7, 1e-5);
    close(
        b.boundary_work_on_electrons_j,
        264.60915634645653,
        3e-7,
        1e-5,
    );
    close(b.neutron_amount_mol, 3.890493548614405e-8, 3e-7, 1e-20);
    assert_eq!(b.external_heavy_energy_j, 0.0);
    assert_eq!(b.external_electron_energy_j, 0.0);
    // Check the physical equation directly, independently of report residuals.
    let before: f64 = pool_energies(&state).into_iter().sum();
    let after: f64 = pool_energies(result.state()).into_iter().sum();
    close(
        after - before + b.neutron_energy_escaped_j + b.photon_energy_escaped_j,
        b.nuclear_energy_released_j + b.boundary_work_on_heavy_j + b.boundary_work_on_electrons_j,
        2e-12,
        1e-6,
    );
}

#[test]
fn g3_constant_volume_history_agrees_with_rigid_network_and_performs_no_work() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let path = history(&[(0.0, 0.025), (0.2, 0.025)]);
    let transport = transport(&state);
    let rigid = mixed(0.025, 0.2, controls());
    let held = advance_fusion_network_with_volume(
        &state,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        &path,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(rigid.volume_history(), None);
    assert_eq!(held.volume_history(), Some(&path));
    assert_eq!(held.balance().boundary_work_on_heavy_j, 0.0);
    assert_eq!(held.balance().boundary_work_on_electrons_j, 0.0);
    assert_eq!(held.state().volume(), rigid.state().volume());
    for (actual, expected) in held.state().amounts().iter().zip(rigid.state().amounts()) {
        close(actual.value(), expected.value(), 3e-9, 1e-20);
    }
    for (actual, expected) in pool_energies(held.state())
        .into_iter()
        .zip(pool_energies(rigid.state()))
    {
        close(actual, expected, 3e-9, 1e-7);
    }
    for (actual, expected) in held
        .balance()
        .reaction_extents_mol
        .into_iter()
        .zip(rigid.balance().reaction_extents_mol)
    {
        close(actual, expected, 3e-8, 1e-20);
    }
    check_reference(&held, &B);
    check_accounts(&held);
}

#[test]
fn g4_volume_mismatch_intermediate_physical_domain_and_cancellation_refuse_transactionally() {
    let (plasma, state) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let original = state.clone();
    let original_plasma = plasma.clone();
    let thermal = transport(&state);
    let original_thermal = thermal.clone();
    // A comparison midpoint must be representable. Otherwise step doubling
    // could compare a zero-length half plus the full step and report zero error.
    let unresolved = history(&[(0.0, 0.025), (f64::from_bits(1), 0.0125)]);
    assert_eq!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::NONE,
            heating(),
            &FusionNetworkTransport::default(),
            &unresolved,
            controls(),
            || false,
        )
        .unwrap_err(),
        FusionNetworkError::StepBudget
    );
    let mismatch = history(&[(0.0, 0.03), (0.1, 0.015)]);
    assert!(matches!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            mixed_heating(),
            &thermal,
            &mismatch,
            controls(),
            || false,
        ),
        Err(FusionNetworkError::Invalid(_))
    ));

    // Both authored endpoint states would be admitted, but compression first
    // raises kBT_e/(m_e*c^2) above the declared .01 nonrelativistic screen.
    // Returning to the original volume cannot erase that intermediate refusal.
    let excursion = history(&[(0.0, 0.025), (0.05, 0.0125), (0.1, 0.025)]);
    let mut restricted = controls();
    restricted
        .plasma_limits
        .maximum_thermal_energy_rest_energy_ratio = 0.01;
    restricted.minimum_step_s = 1e-8;
    assert!(matches!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::NONE,
            heating(),
            &FusionNetworkTransport::default(),
            &excursion,
            restricted,
            || false,
        ),
        Err(FusionNetworkError::OutsideDomain)
            | Err(FusionNetworkError::Plasma(DtFusionPlasmaErrorV1::Domain(_)))
    ));

    let path = history(&[(0.0, 0.025), (0.1, 0.0125)]);
    let original_path = path.clone();
    // A coarse cooling trial may overdraw a thermal pool. That is a retryable
    // physical-domain crossing, not a corrupted arithmetic/account error.
    let mut cooled_controls = controls();
    cooled_controls.minimum_step_s = 1e-8;
    assert_eq!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::NONE,
            FusionNetworkHeating {
                heavy_power: Power::new(-1e10),
                ..heating()
            },
            &FusionNetworkTransport::default(),
            &path,
            cooled_controls,
            || false,
        )
        .unwrap_err(),
        FusionNetworkError::OutsideDomain
    );
    assert_eq!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            mixed_heating(),
            &thermal,
            &path,
            controls(),
            || true,
        )
        .unwrap_err(),
        FusionNetworkError::Cancelled
    );
    let mut polls = 0;
    assert_eq!(
        advance_fusion_network_with_volume(
            &state,
            &plasma,
            FusionNetworkChannels::ALL,
            mixed_heating(),
            &thermal,
            &path,
            controls(),
            || {
                polls += 1;
                polls >= 500
            },
        )
        .unwrap_err(),
        FusionNetworkError::Cancelled
    );
    assert!(polls >= 500);
    assert_eq!(state, original);
    assert_eq!(plasma, original_plasma);
    assert_eq!(thermal, original_thermal);
    assert_eq!(path, original_path);
}
