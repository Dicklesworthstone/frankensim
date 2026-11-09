use super::*;

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

fn controls() -> SaturableRecyclingLimits {
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
    model: SingletTripletDecay,
    external: &[SaturablePumpPath<1>],
    emitted: SaturableEmissionPath<1>,
) -> SaturablePhotonRecycling<1> {
    SaturablePhotonRecycling::try_new(
        &[model],
        external,
        &[emitted],
        64.0 * f64::EPSILON,
        1,
        2,
        2,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_local_pump_matches_exact_finite_two_state_transient_and_integrated_events() {
    let count = 1e12;
    let energy = 4e-19;
    let photon_energy = 5e-19;
    let pump_rate: f64 = 5.0;
    let radiative_rate: f64 = 2.0;
    let nonradiative_rate: f64 = 1.0;
    let duration: f64 = 0.4;
    let initial = [[0.85 * count, 0.15 * count, 0.0]];
    let model = SingletTripletDecay::try_new(
        [energy, 0.75 * energy],
        [radiative_rate, 0.0],
        [nonradiative_rate, 0.0],
        [0.0; 2],
    )
    .unwrap();
    let local = LocalPhotonPump {
        rates_s: [pump_rate, 0.0],
        photon_energies_j: [photon_energy, 0.0],
    };
    let beams = SaturablePumpBeams::<1>::try_new(&[model], &[], 1, 1, 1, || false).unwrap();
    let recycling = prepare(
        model,
        &[],
        SaturableEmissionPath {
            source_cell: 0,
            source_state: 0,
            photon_fractions: [1.0],
            path: SaturablePumpPath {
                photon_energies_j: [energy],
                segments: vec![],
            },
            terminal: EmissionTerminal::Escape,
        },
    );
    let direct = beams
        .advance_with_local_pumps(&initial, &[], &[local], duration, controls().pump, || false)
        .unwrap();
    let step = recycling
        .advance_with_local_pumps(
            &initial,
            &[],
            &[local],
            duration,
            &[None],
            controls(),
            || false,
        )
        .unwrap();

    // S'=W*(N-S)-(kr+knr)*S is a scalar linear equation. Its exact
    // population-time integrals count repeated pumping after either decay;
    // no native kinetic, optical, or temporal solver supplies this reference.
    let relaxation = pump_rate + radiative_rate + nonradiative_rate;
    let equilibrium = pump_rate * count / relaxation;
    let transient_integral = -(-relaxation * duration).exp_m1() / relaxation;
    let singlets = equilibrium + (initial[0][1] - equilibrium) * (-relaxation * duration).exp();
    let singlet_residence =
        equilibrium * duration + (initial[0][1] - equilibrium) * transient_integral;
    let ground_residence = count * duration - singlet_residence;
    let pumped = pump_rate * ground_residence;
    let radiative = radiative_rate * singlet_residence;
    let nonradiative = nonradiative_rate * singlet_residence;
    assert!(pumped > count);
    for (cell, local_step) in [
        (&direct.cells[0], &direct.local_pumps[0]),
        (&step.cells[0], &step.local_pumps[0]),
    ] {
        close(cell.populations[0], count - singlets, count, 3e-10);
        close(cell.populations[1], singlets, count, 3e-10);
        assert_eq!(cell.populations[2], 0.0);
        assert!(cell.populations.iter().all(|&value| value >= 0.0));
        close(cell.populations.iter().sum(), count, count, 3e-10);
        close(cell.pump_events[0], pumped, pumped, 3e-10);
        close(cell.radiative_events[0], radiative, radiative, 3e-10);
        close(
            cell.nonradiative_events[0],
            nonradiative,
            nonradiative,
            3e-10,
        );
        close(local_step.events[0], pumped, pumped, 3e-10);
        close(
            local_step.absorbed_energy_j[0],
            pumped * photon_energy,
            pumped * photon_energy,
            3e-10,
        );
        close(
            local_step.relaxation_heat_j[0],
            pumped * (photon_energy - energy),
            pumped * photon_energy,
            3e-10,
        );
        close(
            cell.absorbed_pump_energy_j[0],
            local_step.absorbed_energy_j[0],
            pumped * photon_energy,
            3e-10,
        );
        close(
            cell.pump_relaxation_heat_j[0],
            local_step.relaxation_heat_j[0],
            pumped * photon_energy,
            3e-10,
        );
        close(
            cell.nonradiative_heat_j[0],
            nonradiative * energy,
            nonradiative * energy,
            3e-10,
        );
        assert_eq!(local_step.events[1], 0.0);
        assert_eq!(local_step.absorbed_energy_j[1], 0.0);
        assert_eq!(local_step.relaxation_heat_j[1], 0.0);
    }
    assert_eq!(step.external_pump_events, [[0.0; 2]]);
    assert_eq!(step.recycled_pump_events, [[0.0; 2]]);
    assert_eq!(step.external_absorbed_pump_energy_j, [[0.0; 2]]);
    assert_eq!(step.recycled_absorbed_pump_energy_j, [[0.0; 2]]);
    close(
        step.escaped_energy_j[0][0],
        radiative * energy,
        radiative * energy,
        3e-10,
    );
    let energy_scale = initial[0][1] * energy + pumped * photon_energy;
    assert!(direct.residual_j.abs() < 2e-10 * energy_scale);
    assert!(step.residual_j.abs() < 2e-10 * energy_scale);

    // Inactive authored local fields preserve the existing no-local route.
    let old = beams
        .advance(&initial, &[], duration, controls().pump, || false)
        .unwrap();
    let zero = beams
        .advance_with_local_pumps(
            &initial,
            &[],
            &[LocalPhotonPump::default()],
            duration,
            controls().pump,
            || false,
        )
        .unwrap();
    assert!(old.local_pumps.is_empty());
    assert_eq!(old.cells, zero.cells);
    assert_eq!(old.beams, zero.beams);
    assert_eq!(old.residual_j.to_bits(), zero.residual_j.to_bits());
    assert!(zero.local_pumps.is_empty());
    assert!(
        beams
            .advance_with_local_pumps(
                &initial,
                &[],
                &[local; 2],
                duration,
                controls().pump,
                || false,
            )
            .is_err()
    );
    let underenergy = LocalPhotonPump {
        photon_energies_j: [0.5 * energy, 0.0],
        ..local
    };
    assert!(
        recycling
            .advance_with_local_pumps(
                &initial,
                &[],
                &[underenergy],
                duration,
                &[None],
                controls(),
                || false,
            )
            .is_err()
    );
    assert_eq!(initial, [[0.85 * count, 0.15 * count, 0.0]]);
    eprintln!(
        "local finite transient: events/N={:.10e}, S/N={:.10e}, residual_j={:.10e}",
        step.local_pumps[0].events[0] / count,
        step.cells[0].populations[1] / count,
        step.residual_j,
    );
}

#[test]
fn g1_local_beam_recycling_and_pairs_share_one_independent_stationary_inventory() {
    let count = 1e12;
    let singlet_energy = 4e-19;
    let triplet_energy = 3e-19;
    let beam_energy = 5e-19;
    let duration = 0.2;
    let radiative_rate = 50.0;
    let yield_per_pair = 0.8;
    let initial = [[0.5 * count, 0.2 * count, 0.3 * count]];
    let ground: f64 = initial[0][0];
    let singlets = initial[0][1];
    let triplets = initial[0][2];
    let local = LocalPhotonPump {
        rates_s: [1.0, 2.0],
        photon_energies_j: [4.5e-19, 3.25e-19],
    };
    let local_rates = local.rates_s.map(|rate| rate * ground);
    let destination_fractions = [0.2, 0.8];
    let return_depth = 2.0 / count;
    let beam_depth = 1.0 / count;
    let escape_fraction = (-return_depth * ground).exp();
    let radiative = radiative_rate * singlets;
    // S'=Ls+fs*A+phi*p-R*escape=0 and T'=Lt+ft*A-2*p=0.
    // Solving these two independent flow equations fixes a positive external
    // absorption A and pair rate p before either transport owner is called.
    let beam_absorption =
        (radiative * escape_fraction - local_rates[0] - 0.5 * yield_per_pair * local_rates[1])
            / (destination_fractions[0] + 0.5 * yield_per_pair * destination_fractions[1]);
    let pair_rate = 0.5 * (local_rates[1] + destination_fractions[1] * beam_absorption);
    assert!(beam_absorption > 0.0 && pair_rate > 0.0);
    let incident = beam_absorption / -(-beam_depth * ground).exp_m1();
    let law = TripletAnnihilation::try_new(pair_rate / (triplets * triplets), 1.0, yield_per_pair)
        .unwrap();
    let model = SingletTripletDecay::try_new(
        [singlet_energy, triplet_energy],
        [radiative_rate, 0.0],
        [0.0; 2],
        [0.0; 2],
    )
    .unwrap();
    let external = SaturablePumpPath {
        photon_energies_j: [beam_energy],
        segments: vec![SaturablePumpSegment {
            cell_index: 0,
            optical_depth_per_ground: [destination_fractions.map(|fraction| fraction * beam_depth)],
            passive_optical_depth: [0.0],
        }],
    };
    let emitted = SaturableEmissionPath {
        source_cell: 0,
        source_state: 0,
        photon_fractions: [1.0],
        path: SaturablePumpPath {
            photon_energies_j: [singlet_energy],
            segments: vec![SaturablePumpSegment {
                cell_index: 0,
                optical_depth_per_ground: [[return_depth, 0.0]],
                passive_optical_depth: [0.0],
            }],
        },
        terminal: EmissionTerminal::Escape,
    };
    let network = prepare(model, &[external], emitted);
    let step = network
        .advance_with_local_pumps(
            &initial,
            &[[incident]],
            &[local],
            duration,
            &[Some(law)],
            controls(),
            || false,
        )
        .unwrap();
    let cell = &step.cells[0];
    let local_events = local_rates.map(|rate| rate * duration);
    let beam_events = destination_fractions.map(|fraction| fraction * beam_absorption * duration);
    let recycled_events = [radiative * (1.0 - escape_fraction) * duration, 0.0];
    let pair_events = pair_rate * duration;
    let escaped_energy = radiative * escape_fraction * duration * singlet_energy;
    let state_energies = [singlet_energy, triplet_energy];
    for state in 0..3 {
        close(cell.populations[state], initial[0][state], count, 2e-8);
        assert!(cell.populations[state] >= 0.0);
    }
    close(cell.populations.iter().sum(), count, count, 2e-10);
    for target in 0..2 {
        let local_energy = local_events[target] * local.photon_energies_j[target];
        let external_energy = beam_events[target] * beam_energy;
        let recycled_energy = recycled_events[target] * singlet_energy;
        let local_heat =
            local_events[target] * (local.photon_energies_j[target] - state_energies[target]);
        let beam_heat = beam_events[target] * (beam_energy - state_energies[target]);
        close(
            step.local_pumps[0].events[target],
            local_events[target],
            count,
            2e-8,
        );
        close(
            step.local_pumps[0].absorbed_energy_j[target],
            local_energy,
            count * beam_energy,
            2e-8,
        );
        close(
            step.local_pumps[0].relaxation_heat_j[target],
            local_heat,
            count * beam_energy,
            2e-8,
        );
        close(
            step.external_pump_events[0][target],
            beam_events[target],
            count,
            2e-8,
        );
        close(
            step.recycled_pump_events[0][target],
            recycled_events[target],
            count,
            2e-8,
        );
        close(
            step.external_absorbed_pump_energy_j[0][target],
            external_energy,
            count * beam_energy,
            2e-8,
        );
        close(
            step.recycled_absorbed_pump_energy_j[0][target],
            recycled_energy,
            count * beam_energy,
            2e-8,
        );
        close(
            step.external_pump_relaxation_heat_j[0][target],
            beam_heat,
            count * beam_energy,
            2e-8,
        );
        close(
            cell.pump_events[target],
            local_events[target] + beam_events[target] + recycled_events[target],
            count,
            2e-8,
        );
        close(
            cell.absorbed_pump_energy_j[target],
            local_energy + external_energy + recycled_energy,
            count * beam_energy,
            2e-8,
        );
        close(
            cell.pump_relaxation_heat_j[target],
            local_heat + beam_heat,
            count * beam_energy,
            2e-8,
        );
    }
    close(cell.annihilation_events, pair_events, pair_events, 2e-8);
    close(
        cell.annihilation_singlets,
        yield_per_pair * pair_events,
        pair_events,
        2e-8,
    );
    let pair_heat = (2.0 * triplet_energy - yield_per_pair * singlet_energy) * pair_events;
    close(cell.annihilation_heat_j, pair_heat, pair_heat, 2e-8);
    close(
        cell.radiative_events[0],
        radiative * duration,
        radiative * duration,
        2e-8,
    );
    close(
        step.external_beams[0].incident_photons[0],
        incident * duration,
        incident * duration,
        2e-8,
    );
    close(
        step.external_beams[0].absorbed_photons[0],
        beam_absorption * duration,
        incident * duration,
        2e-8,
    );
    close(
        step.external_beams[0].transmitted_photons[0],
        (incident - beam_absorption) * duration,
        incident * duration,
        2e-8,
    );
    close(
        step.emission_paths[0].incident_photons[0],
        radiative * duration,
        radiative * duration,
        2e-8,
    );
    close(
        step.emission_paths[0].absorbed_photons[0],
        recycled_events[0],
        radiative * duration,
        2e-8,
    );
    close(
        step.escaped_energy_j[0][0],
        escaped_energy,
        escaped_energy,
        2e-8,
    );
    assert_eq!(step.recycled_pump_relaxation_heat_j, [[0.0; 2]]);
    assert_eq!(step.source_relaxation_heat_j, [[0.0; 2]]);
    assert_eq!(step.passive_heat_j, [0.0]);
    assert_eq!(cell.nonradiative_events, [0.0; 2]);
    assert_eq!(cell.crossing_events, [0.0; 2]);
    let local_input: f64 = (0..2)
        .map(|target| local_events[target] * local.photon_energies_j[target])
        .sum();
    let heat: f64 = cell.pump_relaxation_heat_j.iter().sum::<f64>() + cell.annihilation_heat_j;
    let net_input = beam_absorption * duration * beam_energy + local_input;
    close(
        heat + step.escaped_energy_j[0][0],
        net_input,
        net_input,
        2e-8,
    );
    let scale = singlets * singlet_energy
        + triplets * triplet_energy
        + incident * duration * beam_energy
        + local_input;
    assert!(step.residual_j.abs() < 2e-10 * scale);
    eprintln!(
        "local/beam/recycling/pair equilibrium: substeps={}, evaluations={}, local_events={:?}, pair_events={:.10e}, residual_j={:.10e}",
        step.substeps,
        step.cell_evaluations,
        step.local_pumps[0].events,
        cell.annihilation_events,
        step.residual_j,
    );
}

#[test]
fn g1_local_pair_equilibrium_uses_the_direct_nonlinear_beam_owner() {
    let count = 1e12;
    let energies = [4e-19, 3e-19];
    let initial = [[0.5 * count, 0.2 * count, 0.3 * count]];
    let duration = 0.25;
    let radiative_rate = 2.0;
    let yield_per_pair = 0.5;
    let pair_rate = 0.4 * count;
    // At this authored inventory, the exact stationary balances are
    // Lt=2*p and Ls+phi*p=kr*S; the ground equation is their negative sum.
    let local = LocalPhotonPump {
        rates_s: [
            (radiative_rate * initial[0][1] - yield_per_pair * pair_rate) / initial[0][0],
            2.0 * pair_rate / initial[0][0],
        ],
        photon_energies_j: [4.5e-19, 3.5e-19],
    };
    let law = TripletAnnihilation::try_new(
        pair_rate / (initial[0][2] * initial[0][2]),
        1.0,
        yield_per_pair,
    )
    .unwrap();
    let model =
        SingletTripletDecay::try_new(energies, [radiative_rate, 0.0], [0.0; 2], [0.0; 2]).unwrap();
    let network = SaturablePumpBeams::<1>::try_new(&[model], &[], 1, 1, 1, || false).unwrap();
    let limits = controls();
    let step = network
        .advance_with_local_pumps_and_annihilation(
            &initial,
            &[],
            &[local],
            duration,
            &[Some(law)],
            limits.pump,
            limits.annihilation,
            || false,
        )
        .unwrap();
    let cell = &step.cells[0];
    let pairs = pair_rate * duration;
    let radiative = radiative_rate * initial[0][1] * duration;
    let local_events = [radiative - yield_per_pair * pairs, 2.0 * pairs];
    for state in 0..3 {
        close(cell.populations[state], initial[0][state], count, 2e-8);
        assert!(cell.populations[state] >= 0.0);
    }
    close(cell.populations.iter().sum(), count, count, 2e-10);
    for target in 0..2 {
        let absorbed = local_events[target] * local.photon_energies_j[target];
        let relaxation =
            local_events[target] * (local.photon_energies_j[target] - energies[target]);
        close(cell.pump_events[target], local_events[target], count, 2e-8);
        close(
            step.local_pumps[0].events[target],
            local_events[target],
            count,
            2e-8,
        );
        close(
            step.local_pumps[0].absorbed_energy_j[target],
            absorbed,
            count * energies[0],
            2e-8,
        );
        close(
            step.local_pumps[0].relaxation_heat_j[target],
            relaxation,
            count * energies[0],
            2e-8,
        );
        close(
            cell.absorbed_pump_energy_j[target],
            absorbed,
            count * energies[0],
            2e-8,
        );
        close(
            cell.pump_relaxation_heat_j[target],
            relaxation,
            count * energies[0],
            2e-8,
        );
    }
    close(cell.annihilation_events, pairs, pairs, 2e-8);
    close(
        cell.annihilation_singlets,
        yield_per_pair * pairs,
        pairs,
        2e-8,
    );
    close(
        cell.annihilation_heat_j,
        (2.0 * energies[1] - yield_per_pair * energies[0]) * pairs,
        count * energies[0],
        2e-8,
    );
    close(cell.radiative_events[0], radiative, radiative, 2e-8);
    assert!(step.beams.is_empty());
    assert_eq!(step.passive_heat_j, [0.0]);
    assert_eq!(cell.nonradiative_events, [0.0; 2]);
    let local_input: f64 = (0..2)
        .map(|target| local_events[target] * local.photon_energies_j[target])
        .sum();
    let physical_output = cell.radiative_decay_energy_j.iter().sum::<f64>()
        + cell.pump_relaxation_heat_j.iter().sum::<f64>()
        + cell.annihilation_heat_j;
    close(physical_output, local_input, local_input, 2e-8);
    let scale = initial[0][1] * energies[0] + initial[0][2] * energies[1] + local_input;
    assert!(step.residual_j.abs() < 2e-10 * scale);
    eprintln!(
        "local/direct-pair equilibrium: substeps={}, evaluations={}, pair_events={:.10e}, residual_j={:.10e}",
        step.substeps, step.cell_evaluations, cell.annihilation_events, step.residual_j,
    );
}
