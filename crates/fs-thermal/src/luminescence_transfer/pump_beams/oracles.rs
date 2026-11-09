//! Independent analytic checks for shared counterpropagating pumps.

use super::{
    PumpBeamsLimits, SaturablePumpBeams, SaturablePumpPath, SaturablePumpSegment,
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
