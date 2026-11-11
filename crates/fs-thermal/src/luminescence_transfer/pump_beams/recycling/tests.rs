use super::*;

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

fn model(energy: f64, radiative: f64) -> SingletTripletDecay {
    SingletTripletDecay::try_new(
        [energy, 0.75 * energy],
        [radiative, 0.0],
        [0.0; 2],
        [0.0; 2],
    )
    .unwrap()
}

fn path(cell: usize, depth: f64, photon: f64) -> SaturablePumpPath<1> {
    SaturablePumpPath {
        photon_energies_j: [photon],
        segments: vec![SaturablePumpSegment {
            cell_index: cell,
            optical_depth_per_ground: [[depth, 0.0]],
            passive_optical_depth: [0.0],
        }],
    }
}

fn emission(source: usize, path: SaturablePumpPath<1>) -> SaturableEmissionPath<1> {
    SaturableEmissionPath {
        source_cell: source,
        source_state: 0,
        photon_fractions: [1.0],
        path,
        terminal: EmissionTerminal::Escape,
    }
}

fn limits() -> SaturableRecyclingLimits {
    let mut limits = SaturableRecyclingLimits::default();
    limits.pump.column.relative_tolerance = 2e-4;
    limits.pump.column.maximum_substeps = 2048;
    limits.pump.column.maximum_cell_evaluations = 500_000;
    limits.annihilation.relative_tolerance = 2e-4;
    limits.annihilation.maximum_substeps = 2048;
    limits.maximum_feedback_iterations = 128;
    limits
}

fn prepare(
    models: &[SingletTripletDecay],
    external: &[SaturablePumpPath<1>],
    emitted: &[SaturableEmissionPath<1>],
) -> SaturablePhotonRecycling<1> {
    SaturablePhotonRecycling::try_new(
        models,
        external,
        emitted,
        64.0 * f64::EPSILON,
        8,
        8,
        32,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_self_return_matches_independent_exponential_integral_decay_and_photon_counts() {
    let count = 1e12;
    let energy = 4e-19;
    let network = prepare(
        &[model(energy, 1000.0)],
        &[],
        &[emission(0, path(0, 10.0 / count, energy))],
    );
    let initial = [[0.8 * count, 0.2 * count, 0.0]];
    let step = network
        .advance(&initial, &[], 0.5, &[None], limits(), || false)
        .unwrap();

    // For a resonant self-return ray, G=N-S and S'=-k*S*exp(-a*G).
    // Therefore E1(a*S(t))-E1(a*S0)=k*t*exp(-a*N), while the
    // integrated radiative count is exp(a*N)/a*(exp(-a*S)-exp(-a*S0)).
    // These values were evaluated independently at 85 decimal digits, with
    // E1 inversion by bisection; no time-stepper or production routine enters.
    let singlets = 0.172_908_123_031_005_94 * count;
    let radiative = 92.758_040_254_932_08 * count;
    let escaped = 0.027_091_876_968_994_058 * count;
    close(step.cells[0].populations[0], count - singlets, count, 5e-4);
    close(step.cells[0].populations[1], singlets, singlets, 5e-4);
    close(
        step.cells[0].radiative_events[0],
        radiative,
        radiative,
        5e-4,
    );
    close(
        step.cells[0].pump_events[0],
        radiative - escaped,
        radiative,
        5e-4,
    );
    close(
        step.emitted_energy_j[0][0],
        radiative * energy,
        radiative * energy,
        5e-4,
    );
    close(
        step.escaped_energy_j[0][0],
        escaped * energy,
        escaped * energy,
        5e-4,
    );
    close(
        step.emission_paths[0].incident_photons[0],
        radiative,
        radiative,
        5e-4,
    );
    close(
        step.emission_paths[0].absorbed_photons[0],
        radiative - escaped,
        radiative,
        5e-4,
    );
    close(
        step.emission_paths[0].transmitted_photons[0],
        escaped,
        escaped,
        5e-4,
    );
    close(step.cells[0].populations.iter().sum(), count, count, 2e-10);
    assert_eq!(step.cells[0].populations[2], 0.0);
    assert_eq!(step.source_relaxation_heat_j, [[0.0; 2]]);
    assert_eq!(step.unresolved_energy_j, [[0.0; 2]]);
    assert!(step.residual_j.abs() < 2e-10 * count * energy);
    eprintln!(
        "finite self-return: substeps={}, evaluations={}, S/N={:.10e}, emissions/N={:.10e}, residual_j={:.10e}",
        step.substeps,
        step.cell_evaluations,
        step.cells[0].populations[1] / count,
        step.cells[0].radiative_events[0] / count,
        step.residual_j,
    );
}

#[test]
fn g1_emitted_fluence_saturates_a_finite_receiver_with_separate_stokes_heat() {
    let source_count = 1e12;
    let receiver_count = 0.5e12;
    let initial_source = 0.8 * source_count;
    let source_energy = 5e-19;
    let photon_energy = 4e-19;
    let receiver_energy = 3e-19;
    let rate: f64 = 2.0;
    let duration: f64 = 0.9;
    let depth = 4.0 / receiver_count;
    let models = [model(source_energy, rate), model(receiver_energy, 0.0)];
    let network = prepare(&models, &[], &[emission(0, path(1, depth, photon_energy))]);
    let initial = [
        [source_count - initial_source, initial_source, 0.0],
        [receiver_count, 0.0, 0.0],
    ];
    let step = network
        .advance(&initial, &[], duration, &[None; 2], limits(), || false)
        .unwrap();
    let remaining_source = initial_source * (-rate * duration).exp();
    let fluence = initial_source * -(-rate * duration).exp_m1();
    // dG/dF=-(1-exp(-a*G)) gives exp(a*G)-1=(exp(a*G0)-1)*exp(-a*F).
    let remaining_ground =
        ((depth * receiver_count).exp_m1() * (-depth * fluence).exp()).ln_1p() / depth;
    let absorbed = receiver_count - remaining_ground;
    let escaped = fluence - absorbed;
    close(
        step.cells[0].populations[1],
        remaining_source,
        source_count,
        5e-4,
    );
    close(
        step.cells[1].populations[0],
        remaining_ground,
        receiver_count,
        5e-4,
    );
    close(step.cells[1].populations[1], absorbed, receiver_count, 5e-4);
    assert!(step.cells[1].populations[1] < receiver_count);
    assert!(step.cells[1].populations[1] > 0.9 * receiver_count);
    close(
        step.emission_paths[0].incident_photons[0],
        fluence,
        fluence,
        5e-4,
    );
    close(
        step.emission_paths[0].absorbed_photons[0],
        absorbed,
        absorbed,
        5e-4,
    );
    close(
        step.emission_paths[0].transmitted_photons[0],
        escaped,
        escaped,
        5e-4,
    );
    close(
        step.source_relaxation_heat_j[0][0],
        fluence * (source_energy - photon_energy),
        fluence * (source_energy - photon_energy),
        5e-4,
    );
    close(
        step.cells[1].pump_relaxation_heat_j[0],
        absorbed * (photon_energy - receiver_energy),
        absorbed * (photon_energy - receiver_energy),
        5e-4,
    );
    close(
        step.escaped_energy_j[0][0],
        escaped * photon_energy,
        escaped * photon_energy,
        5e-4,
    );
    for (cell, count) in step.cells.iter().zip([source_count, receiver_count]) {
        close(cell.populations.iter().sum(), count, count, 2e-10);
    }
    assert!(step.residual_j.abs() < 2e-10 * initial_source * source_energy);
}

#[test]
fn g1_nearly_closed_return_keeps_independent_pumped_equilibrium_and_net_energy() {
    let count = 1e12;
    let energy = 4e-19;
    let pump_energy = 4.5e-19;
    let rate = 1e6;
    let duration = 1.0;
    let initial = [[0.8 * count, 0.2 * count, 0.0]];
    let ground: f64 = initial[0][0];
    let singlets = initial[0][1];
    let external_depth = 1.0 / count;
    let return_depth = 20.0 / count;
    let external_fraction = -(-external_depth * ground).exp_m1();
    let escape_fraction = (-return_depth * ground).exp();
    let external_rate = rate * singlets * escape_fraction / external_fraction;
    let network = prepare(
        &[model(energy, rate)],
        &[path(0, external_depth, pump_energy)],
        &[emission(0, path(0, return_depth, energy))],
    );
    let step = network
        .advance(
            &initial,
            &[[external_rate]],
            duration,
            &[None],
            limits(),
            || false,
        )
        .unwrap();
    let emitted = rate * singlets * duration;
    let escaped = emitted * escape_fraction;
    for state in 0..3 {
        close(
            step.cells[0].populations[state],
            initial[0][state],
            count,
            2e-8,
        );
    }
    close(step.cells[0].radiative_events[0], emitted, emitted, 2e-8);
    close(
        step.emission_paths[0].transmitted_photons[0],
        escaped,
        escaped,
        2e-8,
    );
    close(
        step.external_beams[0].absorbed_photons[0],
        escaped,
        escaped,
        2e-8,
    );
    close(
        step.cells[0].pump_relaxation_heat_j[0],
        escaped * (pump_energy - energy),
        escaped * (pump_energy - energy),
        2e-8,
    );
    // The residual scale is original storage plus external energy, despite
    // 200,000 radiative/re-excitation events per original emitter.
    let external_scale = singlets * energy + external_rate * duration * pump_energy;
    assert!(step.residual_j.abs() < 2e-8 * external_scale);
    eprintln!(
        "finite near-closed equilibrium: substeps={}, evaluations={}, radiative_events_per_emitter={:.10e}, residual_j={:.10e}",
        step.substeps,
        step.cell_evaluations,
        step.cells[0].radiative_events[0] / count,
        step.residual_j,
    );
}

#[test]
fn g4_no_active_return_retains_native_history_and_failed_trials_leave_retry_identical() {
    let energy = 4e-19;
    let original = [model(energy, 2.0)];
    let external = [path(0, 3e-12, energy)];
    let emitted = [emission(
        0,
        SaturablePumpPath {
            photon_energies_j: [energy],
            segments: vec![],
        },
    )];
    let network = prepare(&original, &external, &emitted);
    let initial = [[7e11, 3e11, 0.0]];
    let rates = [[1e12]];
    let duration = 0.3;
    let controls = limits();
    let old = SaturablePumpBeams::try_new(&original, &external, 8, 8, 32, || false)
        .unwrap()
        .advance(&initial, &rates, duration, controls.pump, || false)
        .unwrap();
    let step = network
        .advance(&initial, &rates, duration, &[None], controls, || false)
        .unwrap();
    for state in 0..3 {
        close(
            step.cells[0].populations[state],
            old.cells[0].populations[state],
            1e12,
            5e-10,
        );
    }
    close(
        step.cells[0].radiative_events[0],
        old.cells[0].radiative_events[0],
        1e12,
        5e-10,
    );
    close(
        step.external_beams[0].absorbed_photons[0],
        old.beams[0].absorbed_photons[0],
        1e12,
        5e-10,
    );
    assert!(
        network
            .advance(&initial, &rates, duration, &[None], controls, || true)
            .is_err()
    );
    let mut bounded = controls;
    bounded.pump.column.maximum_cell_evaluations = 1;
    assert!(
        network
            .advance(&initial, &rates, duration, &[None], bounded, || false)
            .is_err()
    );
    let retry = network
        .advance(&initial, &rates, duration, &[None], controls, || false)
        .unwrap();
    assert_eq!(retry.cells, step.cells);
    assert_eq!(retry.external_beams, step.external_beams);
    assert_eq!(retry.emission_paths, step.emission_paths);
    assert_eq!(retry.residual_j.to_bits(), step.residual_j.to_bits());
    assert_eq!(initial, [[7e11, 3e11, 0.0]]);
}

#[test]
fn g1_recycling_and_triplet_pairs_share_an_independent_stationary_heat_balance() {
    let count = 1e12;
    let singlet_energy = 4e-19;
    let triplet_energy = 3e-19;
    let pump_energy = 3.5e-19;
    let radiative_rate = 50.0;
    let yield_per_pair = 0.8;
    let duration = 0.2;
    let initial = [[0.5 * count, 0.2 * count, 0.3 * count]];
    let ground: f64 = initial[0][0];
    let singlets = initial[0][1];
    let triplets = initial[0][2];
    let return_depth = 2.0 / count;
    let external_depth = 1.0 / count;
    let escape_fraction = (-return_depth * ground).exp();
    // Stationary singlets require phi*p=kS*S*escape. Stationary triplets
    // require external absorption=2*p. These independent flow equalities
    // also give G'=0, including the pair's 2-phi returned ground emitters.
    let pair_rate = radiative_rate * singlets * escape_fraction / yield_per_pair;
    let law = TripletAnnihilation::try_new(pair_rate / (triplets * triplets), 1.0, yield_per_pair)
        .unwrap();
    let external_rate = 2.0 * pair_rate / -(-external_depth * ground).exp_m1();
    let mut triplet_path = path(0, 0.0, pump_energy);
    triplet_path.segments[0].optical_depth_per_ground = [[0.0, external_depth]];
    let original = SingletTripletDecay::try_new(
        [singlet_energy, triplet_energy],
        [radiative_rate, 0.0],
        [0.0; 2],
        [0.0; 2],
    )
    .unwrap();
    let network = prepare(
        &[original],
        &[triplet_path],
        &[emission(0, path(0, return_depth, singlet_energy))],
    );
    let step = network
        .advance(
            &initial,
            &[[external_rate]],
            duration,
            &[Some(law)],
            limits(),
            || false,
        )
        .unwrap();
    let pairs = pair_rate * duration;
    let radiative = radiative_rate * singlets * duration;
    let escaped = radiative * escape_fraction;
    let recycled = radiative * (1.0 - escape_fraction);
    let cell = &step.cells[0];
    for state in 0..3 {
        close(cell.populations[state], initial[0][state], count, 2e-8);
    }
    close(cell.annihilation_events, pairs, pairs, 2e-8);
    close(
        cell.annihilation_singlets,
        yield_per_pair * pairs,
        pairs,
        2e-8,
    );
    close(cell.radiative_events[0], radiative, radiative, 2e-8);
    close(cell.pump_events[0], recycled, recycled, 2e-8);
    close(cell.pump_events[1], 2.0 * pairs, 2.0 * pairs, 2e-8);
    let pair_heat = (2.0 * triplet_energy - yield_per_pair * singlet_energy) * pairs;
    let pump_heat = (pump_energy - triplet_energy) * 2.0 * pairs;
    close(cell.annihilation_heat_j, pair_heat, pair_heat, 2e-8);
    close(cell.pump_relaxation_heat_j[1], pump_heat, pump_heat, 2e-8);
    close(
        step.external_beams[0].absorbed_photons[0],
        2.0 * pairs,
        2.0 * pairs,
        2e-8,
    );
    close(
        step.escaped_energy_j[0][0],
        escaped * singlet_energy,
        escaped * singlet_energy,
        2e-8,
    );
    close(
        step.external_absorbed_pump_energy_j[0][1],
        2.0 * pairs * pump_energy,
        2.0 * pairs * pump_energy,
        2e-8,
    );
    close(
        step.recycled_absorbed_pump_energy_j[0][0],
        recycled * singlet_energy,
        recycled * singlet_energy,
        2e-8,
    );
    assert_eq!(cell.nonradiative_events, [0.0; 2]);
    assert_eq!(cell.crossing_events, [0.0; 2]);
    assert_eq!(cell.reverse_crossing_heat_absorbed_j, 0.0);
    let scale = singlets * singlet_energy
        + triplets * triplet_energy
        + external_rate * duration * pump_energy;
    assert!(step.residual_j.abs() < 2e-10 * scale);
    eprintln!(
        "finite pair/recycling equilibrium: substeps={}, evaluations={}, pair_events={:.10e}, residual_j={:.10e}",
        step.substeps, step.cell_evaluations, cell.annihilation_events, step.residual_j,
    );
}

#[test]
fn g1_sixty_five_radiating_cells_recycle_a_nonstationary_ring_on_original_inventories() {
    let cells = 65usize;
    let count = 1e12;
    let energy = 4e-19;
    let models = vec![model(energy, 2.0); cells];
    let emitted: Vec<_> = (0..cells)
        .map(|source| emission(source, path((source + 1) % cells, 2.0 / count, energy)))
        .collect();
    let network = SaturablePhotonRecycling::try_new(
        &models,
        &[],
        &emitted,
        64.0 * f64::EPSILON,
        cells,
        cells,
        cells,
        || false,
    )
    .unwrap();
    let initial = vec![[0.8 * count, 0.2 * count, 0.0]; cells];
    let step = network
        .advance(&initial, &[], 0.2, &vec![None; cells], limits(), || false)
        .unwrap();

    // Every receiver has the same history by ring symmetry, so the independent
    // self-return E1 solution applies to EACH original finite inventory.
    // Here a*N=2, S0/N=.2 and k*t=.4. E1 differences were evaluated with
    // 110-digit fixed-point integer series and scalar bisection:
    // E1(x)-E1(x0)=ln(x0/x)+sum_k[((-x0)^k-(-x)^k)/(k*k!)].
    // Euler's constant cancels. No native kinetic or temporal solver enters.
    let singlets = 0.184_712_975_307_005_3 * count;
    let radiative = 0.076_886_505_698_846_88 * count;
    let escaped = 0.015_287_024_692_994_692 * count;
    assert_eq!(step.cells.len(), cells);
    assert_eq!(step.emission_paths.len(), cells);
    // A real transient, with appreciable loss, cannot use a frozen equilibrium
    // to satisfy this test merely because all source cells happen to be equal.
    assert!(step.cells[0].populations[1] < 0.19 * count);
    for index in 0..cells {
        let cell = &step.cells[index];
        close(cell.populations[0], count - singlets, count, 5e-4);
        close(cell.populations[1], singlets, singlets, 5e-4);
        close(cell.populations.iter().sum(), count, count, 2e-10);
        assert_eq!(cell.populations[2], 0.0);
        close(cell.radiative_events[0], radiative, radiative, 5e-4);
        close(cell.pump_events[0], radiative - escaped, radiative, 5e-4);
        close(
            step.emission_paths[index].incident_photons[0],
            radiative,
            radiative,
            5e-4,
        );
        close(
            step.emission_paths[index].transmitted_photons[0],
            escaped,
            escaped,
            5e-4,
        );
        close(
            step.escaped_energy_j[index][0],
            escaped * energy,
            escaped * energy,
            5e-4,
        );
    }
    let total_escaped: f64 = step
        .escaped_energy_j
        .iter()
        .map(|cell| cell[0] + cell[1])
        .sum();
    let expected_escaped = cells as f64 * escaped * energy;
    close(total_escaped, expected_escaped, expected_escaped, 5e-4);
    let original_energy = cells as f64 * 0.2 * count * energy;
    assert!(step.residual_j.abs() < 2e-10 * original_energy);
    eprintln!(
        "finite 65-source transient ring: substeps={}, evaluations={}, S/N={:.10e}, escaped_j={:.10e}, residual_j={:.10e}",
        step.substeps,
        step.cell_evaluations,
        step.cells[0].populations[1] / count,
        total_escaped,
        step.residual_j,
    );
}

/// One singlet cell of excitation energy `energy` that radiates on one
/// escaping emission path of the same photon energy (no Stokes shift).
fn planted_network(energy: f64) -> SaturablePhotonRecycling<1> {
    prepare(
        &[model(energy, 1e3)],
        &[],
        &[emission(
            0,
            SaturablePumpPath {
                photon_energies_j: [energy],
                segments: vec![],
            },
        )],
    )
}

/// A history of planted_network() built directly, not computed: from a full
/// singlet, `events` radiative events, `incident` photons carried on the
/// emission path and escaped whole, and the carried bound `bound`.
fn planted_history(energy: f64, events: f64, incident: f64, bound: f64) -> History<1> {
    let cell = SaturableSingletTripletStep {
        populations: [events, 1.0 - events, 0.0],
        pump_events: [0.0; 2],
        absorbed_pump_energy_j: [0.0; 2],
        pump_relaxation_heat_j: [0.0; 2],
        radiative_events: [events, 0.0],
        nonradiative_events: [0.0; 2],
        crossing_events: [0.0; 2],
        stored_energy_j: [(1.0 - events) * energy, 0.0],
        radiative_decay_energy_j: [events * energy, 0.0],
        nonradiative_heat_j: [0.0; 2],
        crossing_heat_j: 0.0,
        reverse_crossing_heat_absorbed_j: 0.0,
        annihilation_events: 0.0,
        annihilation_singlets: 0.0,
        annihilation_heat_j: 0.0,
        population_residual: 0.0,
        residual_j: 0.0,
    };
    let beam = PumpBeamStep {
        incident_photons: [incident],
        transmitted_photons: [incident],
        absorbed_photons: [0.0],
        passive_absorbed_photons: [0.0],
        incident_energy_j: [incident * energy],
        transmitted_energy_j: [incident * energy],
        absorbed_energy_j: [0.0],
        passive_absorbed_energy_j: [0.0],
        photon_residual: [0.0],
    };
    History {
        beams: PumpBeamsStep {
            cells: vec![cell],
            passive_heat_j: vec![0.0],
            beams: vec![beam],
            local_pumps: vec![],
            residual_j: 0.0,
            substeps: 1,
            cell_evaluations: 0,
            coupling_iterations: 0,
            max_energy_discrepancy_j: 0.0,
        },
        external_absorbed: vec![[0.0; 2]],
        recycled_absorbed: vec![[0.0; 2]],
        external_events: vec![[0.0; 2]],
        recycled_events: vec![[0.0; 2]],
        external_relaxation: vec![[0.0; 2]],
        recycled_relaxation: vec![[0.0; 2]],
        external_passive: vec![0.0],
        recycled_passive: vec![0.0],
        source_allowance_j: bound,
    }
}

/// Planted negative for the composed source closure. The History is built
/// directly, not computed: one singlet cell whose directly outgoing emission
/// path carried more photons than its radiative events emitted. A composed
/// mismatch 10 % over the bound the history carries is refused with exactly
/// that residual and that bound; the same mismatch under the bound passes.
#[test]
fn g4_composed_source_closure_refuses_ten_percent_over_the_carried_bound() {
    let energy = 2e-19;
    let network = planted_network(energy);
    let events = 0.5;
    let incident = 0.5 * (1.0 + 1e-9);
    let history = |bound: f64| planted_history(energy, events, incident, bound);
    let initial = [[0.0, 1.0, 0.0]];
    let tolerance = 1e-12;
    let mismatch = (events - incident).abs() * energy;
    let over = mismatch / 1.1;
    // The single-step value (tolerance x the initial singlet storage) stays
    // below the carried bound, so the carried bound is the one tested.
    assert!(tolerance * energy < over);
    match network.admit_composed_source(&initial, &history(over), tolerance, &mut || false) {
        Err(SingletTripletError::EnergyBalance {
            residual_j,
            allowance_j,
        }) => {
            assert_eq!(residual_j, mismatch);
            assert_eq!(allowance_j, over);
        }
        other => panic!("expected the composed source refusal, got {other:?}"),
    }
    network
        .admit_composed_source(&initial, &history(mismatch / 0.9), tolerance, &mut || false)
        .unwrap();
}

/// Planted negative for the total energy balance. The History is built
/// directly, not computed: the planted cell with no source mismatch, whose
/// escaped energy exceeds what it emitted by `excess`. Every term is a power
/// of two times the photon energy 2^-62 J, so the balance and its scale are
/// exact. A balance 10 % over tolerance x scale + the carried bound is refused
/// with exactly that residual and that allowance. The same balance under the
/// sum is admitted, although the single-step value alone would refuse it.
#[test]
fn g4_total_balance_refuses_ten_percent_over_the_single_step_value_plus_the_carried_bound() {
    let energy = 0.5_f64.powi(62);
    let network = planted_network(energy);
    let events = 0.5;
    let excess = energy * 0.5_f64.powi(36);
    let initial = [[0.0, 1.0, 0.0]];
    let tolerance = 1e-12;
    let source_relaxation = [[0.0; 2]];
    let escaped = [[events * energy + excess, 0.0]];
    let unresolved = [[0.0; 2]];
    // The scale's nonzero terms: the initial singlet, the stored remainder and
    // the escaped energy.
    let single = tolerance * (2.0 * energy + excess);
    let over = excess / 1.1 - single;
    assert!(0.0 < over && single < excess);
    let balance = |bound: f64| {
        network.admit_total_balance(
            &initial,
            &planted_history(energy, events, events, bound),
            &source_relaxation,
            &escaped,
            &unresolved,
            tolerance,
            &mut || false,
        )
    };
    match balance(over) {
        Err(SingletTripletError::EnergyBalance {
            residual_j,
            allowance_j,
        }) => {
            assert_eq!(residual_j, -excess);
            assert_eq!(allowance_j, single + over);
        }
        other => panic!("expected the total balance refusal, got {other:?}"),
    }
    assert_eq!(balance(excess / 0.9 - single).unwrap(), -excess);
}
