//! Independent analytic checks for shared counterpropagating and returning pumps.

use super::{
    LocalPhotonPump, PumpBeamsLimits, SaturablePumpBeams, SaturablePumpPath,
    SaturablePumpSegment,
};
use super::recycling::{
    EmissionTerminal, SaturableEmissionPath, SaturablePhotonRecycling,
    SaturableRecyclingLimits,
};
use crate::luminescence::transfer::{
    SingletTripletDecay, pump_column::PumpColumnLimits,
};

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn opposed_paths(
    coefficient: f64,
    split_targets: bool,
) -> [SaturablePumpPath<1>; 2] {
    std::array::from_fn(|beam| {
        let mut depth = [0.0; 2];
        depth[if split_targets { beam } else { 0 }] = coefficient;
        let indices = if beam == 0 { [0, 1] } else { [1, 0] };
        SaturablePumpPath {
            photon_energies_j: [if split_targets && beam == 1 { 4e-19 } else { 5e-19 }],
            segments: indices
                .map(|cell_index| SaturablePumpSegment {
                    cell_index,
                    optical_depth_per_ground: [depth],
                    passive_optical_depth: [0.0],
                })
                .into_iter()
                .collect(),
        }
    })
}

#[test]
fn g3_opposed_one_shot_pumps_share_one_inventory_and_match_exact_fluence() {
    let count = 1e12;
    let coefficient = 1.0 / count;
    let model =
        SingletTripletDecay::try_new([4e-19, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2])
            .unwrap();
    let paths = opposed_paths(coefficient, true);
    let plan = SaturablePumpBeams::try_new(&[model; 2], &paths, 2, 2, 4, || false)
        .unwrap();
    let step = plan
        .advance(
            &[[count, 0.0, 0.0]; 2],
            &[[2.5 * count]; 2],
            1.0,
            PumpBeamsLimits {
                column: PumpColumnLimits {
                    relative_tolerance: 2e-5,
                    maximum_cell_evaluations: 500_000,
                    ..PumpColumnLimits::default()
                },
                maximum_coupling_iterations: 128,
            },
            || false,
        )
        .unwrap();

    // Each beam traverses the same total ground opacity, regardless of its
    // direction. Thus dGtotal/dt = -Rtotal*(1-exp(-a*Gtotal)). With no return
    // and initial optical depth 2, total incident normalized fluence 5 gives
    // uout = ln(1+exp(-2)*(exp(5)-1)). The constants below were independently
    // evaluated at 70 decimal digits, without the transport/kinetic solver.
    let outgoing_each = 1.521_074_135_570_247_3 / coefficient;
    let remaining_fraction = 0.021_074_135_570_247_29;
    let events_each = count * (1.0 - remaining_fraction);
    assert_eq!(step.cells.len(), 2);
    assert_eq!(step.beams.len(), 2);
    for beam in &step.beams {
        close(beam.transmitted_photons[0], outgoing_each, 2e-5);
        close(beam.absorbed_photons[0], events_each, 2e-5);
        assert_eq!(beam.passive_absorbed_photons[0], 0.0);
        assert!(
            beam.photon_residual[0].abs()
                <= 2e-12 * beam.incident_photons[0]
        );
    }
    for cell in &step.cells {
        close(cell.populations[0] / count, remaining_fraction, 2e-5);
        close(cell.populations.iter().sum(), count, 2e-12);
        assert!(cell.populations.into_iter().all(|x| x >= 0.0 && x <= count));
    }
    // Equal photon rates and equal total cross sections preserve mirrored
    // ground populations even though the beams populate different manifolds.
    close(step.cells[0].populations[0], step.cells[1].populations[0], 2e-12);
    close(step.cells[0].populations[1], step.cells[1].populations[2], 2e-12);
    close(step.cells[0].populations[2], step.cells[1].populations[1], 2e-12);
    for target in 0..2 {
        let events = step.cells.iter().map(|cell| cell.pump_events[target]).sum();
        let population = step.cells.iter().map(|cell| cell.populations[target + 1]).sum();
        close(events, events_each, 2e-5);
        close(population, events_each, 2e-5);
        close(events, step.beams[target].absorbed_photons[0], 2e-12);
        let deposited = step.cells.iter().map(|cell| cell.absorbed_pump_energy_j[target]).sum();
        close(
            deposited,
            step.beams[target].absorbed_energy_j[0],
            2e-12,
        );
    }
    let heat: f64 = step.cells.iter().map(|cell| {
        cell.pump_relaxation_heat_j.iter().sum::<f64>()
    }).sum();
    close(heat, 2.0 * events_each * 1e-19, 2e-5);
    assert!(step.passive_heat_j.iter().all(|x| *x == 0.0));
    assert!(step.residual_j.abs() < 2e-17);
    assert!(step.substeps > 1);
    eprintln!(
        "opposed one-shot: substeps={}, cell_evaluations={}, coupling_iterations={}, residual_j={:.9e}",
        step.substeps, step.cell_evaluations, step.coupling_iterations, step.residual_j
    );
}

#[test]
fn g1_opposed_recovering_pumps_keep_exact_steady_state_across_a_million_lifetimes() {
    let count = 1e12;
    let coefficient = 1.0 / count;
    let return_rate = 1e9;
    let duration = 1e-3;
    let photon_rate = 1.5 * return_rate / coefficient;
    let model =
        SingletTripletDecay::try_new(
            [4e-19, 3e-19], [return_rate, 0.0], [0.0; 2], [0.0; 2],
        ).unwrap();
    let paths = opposed_paths(coefficient, false);
    let plan = SaturablePumpBeams::try_new(&[model; 2], &paths, 2, 2, 4, || false)
        .unwrap();

    // For normalized intensities s=a*R/k, opposite rays obey
    // ds+/dτ=-s+/(1+s++s-) and ds-/dτ=+s-/(1+s++s-).
    // Their product is constant, so both fractional transmissions equal T
    // with -ln(T)+(sL+sR)*(1-T)=τ0. Here τ0=2 and sL=sR=1.5;
    // T=W(3*exp(1))/3 and mean G/N=-ln(T)/2. Both mirrored
    // cells have this exact mean steady population, independently of any
    // cell iteration or temporal integrator.
    let transmission = 0.539_214_155_592_024_8;
    let ground_fraction = 0.308_821_233_388_037_15;
    let initial = [count * ground_fraction, count * (1.0 - ground_fraction), 0.0];
    let events = return_rate * initial[1] * duration;
    let step = plan
        .advance(
            &[initial; 2],
            &[[photon_rate]; 2],
            duration,
            PumpBeamsLimits {
                column: PumpColumnLimits {
                    maximum_substeps: 4,
                    maximum_cell_evaluations: 512,
                    relative_tolerance: 1e-10,
                    ..PumpColumnLimits::default()
                },
                maximum_coupling_iterations: 32,
            },
            || false,
        )
        .unwrap();
    for cell in &step.cells {
        close(cell.populations[0], initial[0], 5e-12);
        close(cell.populations[1], initial[1], 5e-12);
        assert_eq!(cell.populations[2], 0.0);
        close(cell.populations.iter().sum(), count, 2e-12);
        close(cell.pump_events[0], events, 5e-12);
        close(cell.radiative_events[0], events, 5e-12);
        close(cell.radiative_decay_energy_j[0], events * 4e-19, 5e-12);
        close(cell.pump_relaxation_heat_j[0], events * 1e-19, 5e-12);
        assert_eq!(cell.pump_events[1], 0.0);
        assert!(cell.pump_events[0] > 1e5 * count);
    }
    for beam in &step.beams {
        close(
            beam.transmitted_photons[0] / beam.incident_photons[0],
            transmission,
            5e-12,
        );
        close(beam.absorbed_photons[0], events, 5e-12);
        assert_eq!(beam.passive_absorbed_photons[0], 0.0);
        assert!(
            beam.photon_residual[0].abs()
                <= 2e-12 * beam.incident_photons[0]
        );
    }
    assert!(step.residual_j.abs() < 1e-12);
    assert!(step.substeps <= 4);
    assert!(step.cell_evaluations <= 512);
    eprintln!(
        "opposed steady recovery: events_per_emitter={:.9e}, substeps={}, cell_evaluations={}, coupling_iterations={}, residual_j={:.9e}",
        events / count, step.substeps, step.cell_evaluations, step.coupling_iterations, step.residual_j
    );
}

#[test]
fn g1_returning_beam_matches_independent_fluence_and_ordered_passive_absorption() {
    let count: f64 = 1e12;
    let coefficient = 1.0 / count;
    let photon_energy = 5e-19;
    let state_energy = 4e-19;
    let passive_depth: f64 = 0.7;
    let photon_rate = 2.0 * count;
    let duration = 0.75;
    let model = SingletTripletDecay::try_new(
        [state_energy, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2],
    ).unwrap();
    let path = SaturablePumpPath {
        photon_energies_j: [photon_energy],
        segments: [
            (0, coefficient, 0.0),
            (1, 0.0, passive_depth),
            (0, coefficient, 0.0),
        ].map(|(cell_index, active, passive)| SaturablePumpSegment {
            cell_index,
            optical_depth_per_ground: [[active, 0.0]],
            passive_optical_depth: [passive],
        }).into_iter().collect(),
    };
    let network = SaturablePumpBeams::try_new(
        &[model; 2], &[path], 2, 1, 3, || false,
    ).unwrap();
    assert!(network.topological.is_none());
    let initial = [[count, 0.0, 0.0], [0.0; 3]];
    let step = network.advance(
        &initial,
        &[[photon_rate]],
        duration,
        PumpBeamsLimits {
            column: PumpColumnLimits {
                relative_tolerance: 2e-5,
                maximum_cell_evaluations: 500_000,
                ..PumpColumnLimits::default()
            },
            maximum_coupling_iterations: 128,
        },
        || false,
    ).unwrap();

    // A -> passive B -> A sees the SAME changing ground population twice.
    // With x=exp(-a*G), r=exp(-p), and u=a*incident_rate*t,
    // dx/du=x*(1-x)*(1+r*x). Integrating its partial fractions gives
    // F(x)=ln(x)-[ln(1-x)+r*ln(1+r*x)]/(1+r), F(x)-F(x0)=u.
    // This scalar inverse contains no native optical/kinetic stage.
    let survival = (-passive_depth).exp();
    let x0 = (-coefficient * count).exp();
    let primitive = |x: f64| {
        x.ln() - ((1.0 - x).ln() + survival * (1.0 + survival * x).ln())
            / (1.0 + survival)
    };
    let target = primitive(x0) + coefficient * photon_rate * duration;
    let mut low = x0;
    let mut high = 1.0;
    for _ in 0..100 {
        let middle = low + 0.5 * (high - low);
        if primitive(middle) < target {
            low = middle;
        } else {
            high = middle;
        }
    }
    let x = low + 0.5 * (high - low);
    let ground = -x.ln() / coefficient;
    let absorbed = count - ground;
    // Passive photons are I*(1-r)*integral(x dt), not attenuation after
    // moving both A segments in front of B. Its independent integral is exact.
    let passive_primitive = |value: f64| {
        ((1.0 + survival * value) / (1.0 - value)).ln()
    };
    let passive = (1.0 - survival) / (coefficient * (1.0 + survival))
        * (passive_primitive(x) - passive_primitive(x0));
    let incident = photon_rate * duration;
    let transmitted = incident - absorbed - passive;
    assert!(ground > 0.0 && ground < 0.5 * count);
    assert!(passive > 0.0 && transmitted > 0.0);
    close(step.cells[0].populations[0], ground, 3e-5);
    close(step.cells[0].populations[1], absorbed, 3e-5);
    assert_eq!(step.cells[0].populations[2], 0.0);
    close(step.cells[0].populations.iter().sum(), count, 2e-12);
    assert!(step.cells[0].populations.iter().all(|&value| value >= 0.0));
    assert_eq!(step.cells[1].populations, [0.0; 3]);
    close(step.cells[0].pump_events[0], absorbed, 3e-5);
    close(step.beams[0].absorbed_photons[0], absorbed, 3e-5);
    close(step.beams[0].passive_absorbed_photons[0], passive, 3e-5);
    close(step.beams[0].transmitted_photons[0], transmitted, 3e-5);
    close(step.beams[0].incident_photons[0], incident, 2e-12);
    close(step.cells[0].absorbed_pump_energy_j[0], absorbed * photon_energy, 3e-5);
    close(
        step.cells[0].pump_relaxation_heat_j[0],
        absorbed * (photon_energy - state_energy),
        3e-5,
    );
    assert_eq!(step.passive_heat_j[0], 0.0);
    close(step.passive_heat_j[1], passive * photon_energy, 3e-5);
    assert_eq!(step.cells[0].radiative_events, [0.0; 2]);
    assert!(step.beams[0].photon_residual[0].abs() < 2e-12 * incident);
    assert!(step.residual_j.abs() < 2e-12 * incident * photon_energy);
    assert!(step.substeps > 1);
    assert!(step.coupling_iterations > 0);
    assert_eq!(initial, [[count, 0.0, 0.0], [0.0; 3]]);
}

#[test]
fn g1_emitted_return_through_passive_cell_shares_one_stationary_source_inventory() {
    let count: f64 = 1e12;
    let state_energy = 4e-19;
    let local_photon_energy = 5e-19;
    let radiative_rate = 2.0;
    let duration = 2.0;
    let initial = [[0.4 * count, 0.6 * count, 0.0], [0.0; 3]];
    let ground = initial[0][0];
    let singlets = initial[0][1];
    let depth_first = 0.7 / count;
    let depth_return = 1.1 / count;
    let passive_depth: f64 = 0.6;
    let first_survival = (-depth_first * ground).exp();
    let passive_survival = (-passive_depth).exp();
    let return_survival = (-depth_return * ground).exp();
    // Ordered Beer factors fix the three dispositions before solving any
    // kinetics. B is between both source-cell visits and receives only light
    // that survived the FIRST source traversal.
    let returned_fraction = (1.0 - first_survival)
        + first_survival * passive_survival * (1.0 - return_survival);
    let passive_fraction = first_survival * (1.0 - passive_survival);
    let escape_fraction = first_survival * passive_survival * return_survival;
    close(returned_fraction + passive_fraction + escape_fraction, 1.0, 2e-15);
    let emission_rate = radiative_rate * singlets;
    let local_rate = emission_rate * (passive_fraction + escape_fraction) / ground;
    let model = SingletTripletDecay::try_new(
        [state_energy, 3e-19], [radiative_rate, 0.0], [0.0; 2], [0.0; 2],
    ).unwrap();
    let passive_model = SingletTripletDecay::try_new(
        [state_energy, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2],
    ).unwrap();
    let emitted = SaturableEmissionPath {
        source_cell: 0,
        source_state: 0,
        photon_fractions: [1.0],
        path: SaturablePumpPath {
            photon_energies_j: [state_energy],
            segments: [
                (0, depth_first, 0.0),
                (1, 0.0, passive_depth),
                (0, depth_return, 0.0),
            ].map(|(cell_index, active, passive)| SaturablePumpSegment {
                cell_index,
                optical_depth_per_ground: [[active, 0.0]],
                passive_optical_depth: [passive],
            }).into_iter().collect(),
        },
        terminal: EmissionTerminal::Escape,
    };
    let network = SaturablePhotonRecycling::try_new(
        &[model, passive_model], &[], &[emitted], 64.0 * f64::EPSILON,
        2, 1, 3, || false,
    ).unwrap();
    let step = network.advance_with_local_pumps(
        &initial,
        &[],
        &[
            LocalPhotonPump {
                rates_s: [local_rate, 0.0],
                photon_energies_j: [local_photon_energy, 0.0],
            },
            LocalPhotonPump::default(),
        ],
        duration,
        &[None; 2],
        SaturableRecyclingLimits::default(),
        || false,
    ).unwrap();
    let emitted = emission_rate * duration;
    let recycled = emitted * returned_fraction;
    let local_events = emitted * (passive_fraction + escape_fraction);
    let passive_heat = emitted * passive_fraction * state_energy;
    let escaped_energy = emitted * escape_fraction * state_energy;
    assert!(emitted > 2.0 * count);
    assert_eq!(step.cells.len(), 2);
    for state in 0..2 {
        close(step.cells[0].populations[state], initial[0][state], 3e-12);
    }
    assert_eq!(step.cells[0].populations[2], 0.0);
    assert_eq!(step.cells[1].populations, [0.0; 3]);
    close(step.cells[0].pump_events[0], emitted, 3e-12);
    close(step.cells[0].radiative_events[0], emitted, 3e-12);
    close(step.local_pumps[0].events[0], local_events, 3e-12);
    close(step.recycled_pump_events[0][0], recycled, 3e-12);
    close(
        step.local_pumps[0].absorbed_energy_j[0],
        local_events * local_photon_energy,
        3e-12,
    );
    close(
        step.cells[0].absorbed_pump_energy_j[0],
        local_events * local_photon_energy + recycled * state_energy,
        3e-12,
    );
    let local_heat = local_events * (local_photon_energy - state_energy);
    close(step.cells[0].pump_relaxation_heat_j[0], local_heat, 3e-12);
    assert_eq!(step.recycled_pump_relaxation_heat_j, vec![[0.0; 2]; 2]);
    assert_eq!(step.source_relaxation_heat_j, vec![[0.0; 2]; 2]);
    assert_eq!(step.external_absorbed_pump_energy_j, vec![[0.0; 2]; 2]);
    assert!(step.external_beams.is_empty());
    close(step.emission_paths[0].incident_photons[0], emitted, 3e-12);
    close(step.emission_paths[0].absorbed_photons[0], recycled, 3e-12);
    close(
        step.emission_paths[0].passive_absorbed_photons[0],
        emitted * passive_fraction,
        3e-12,
    );
    close(step.passive_heat_j[1], passive_heat, 3e-12);
    assert_eq!(step.passive_heat_j[0], 0.0);
    close(step.recycled_passive_heat_j[1], passive_heat, 3e-12);
    close(step.escaped_energy_j[0][0], escaped_energy, 3e-12);
    assert_eq!(step.unresolved_energy_j, vec![[0.0; 2]; 2]);
    close(
        local_heat + passive_heat + escaped_energy,
        local_events * local_photon_energy,
        3e-12,
    );
    assert!(step.residual_j.abs() < 2e-12 * local_events * local_photon_energy);
    assert_eq!(step.substeps, 1);
    assert_eq!(step.cell_evaluations, 2);
}
