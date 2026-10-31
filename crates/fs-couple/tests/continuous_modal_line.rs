//! Exact coupled reflection, original work, transactionality and an independent
//! three-rung clock axis. Spatial/radial authority lives in fs-duct tests.
use core::f64::consts::TAU;
use fs_couple::mm_line::continuous::ContinuousModalLine;
use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::ModalEnd;
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_phs::linear_exact::LinearSplitBudget;

fn line(duration: f64) -> ContinuousModalLine {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let profile = RadiusProfile::exponential(0.025, 10.0, 0.12).unwrap();
    ContinuousModalLine::new(
        &profile,
        &gas,
        ModeCount(2),
        3,
        ModalEnd::Rigid,
        duration,
        LinearSplitBudget {
            max_modes: 124,
            max_omega_dt: 1.0,
        },
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + floor,
        "actual={actual:e}, expected={expected:e}"
    );
}
fn snapshot(line: &ContinuousModalLine) -> Vec<f64> {
    let mut values = line.state().to_vec();
    values.extend_from_slice(line.incoming());
    values.extend_from_slice(line.mean_pressures());
    for port in line.port_accounts() {
        values.extend([
            port.boundary_work,
            port.incident_energy,
            port.outgoing_energy,
            port.outgoing_wave,
        ]);
    }
    values
}
fn state_distance(line: &ContinuousModalLine, a: &[f64], b: &[f64]) -> f64 {
    let source = line.source();
    let n = source.coordinates();
    let delta: Vec<_> = a.iter().zip(b).map(|(x, y)| x - y).collect();
    let mut squared = 0.0;
    for i in 0..n {
        for j in 0..n {
            squared += delta[i] * source.stiffness()[i * n + j] * delta[j]
                + delta[n + i] * source.inverse_mass()[i * n + j] * delta[n + j];
        }
    }
    assert!(squared >= 0.0);
    det::sqrt(squared)
}

#[test]
fn coherent_all_mode_pulse_closes_original_work_and_preserves_within_step_variation() {
    let mut region = line(1.0 / 96_000.0);
    assert_eq!(region.state().len(), 16);
    let (mut incident_total, mut outgoing_total) = (0.0, 0.0);
    let mut free = [0.0; 2];
    let mut direct = [0.0; 4];
    let mut variation_live = false;
    let mut cross_mode_live = false;
    for step in 0..256 {
        let incident = if step < 72 {
            [
                80.0 * det::sin(0.11 * (step as f64 + 0.5)),
                60.0 * det::cos(0.07 * (step as f64 + 0.5)),
            ]
        } else {
            [0.0, 0.0]
        };
        region
            .reflection_relation_into(&mut free, &mut direct)
            .unwrap();
        cross_mode_live |= direct[1].abs() > 1.0e-7;
        let predicted = [
            free[0] + direct[0] * incident[0] + direct[1] * incident[1],
            free[1] + direct[2] * incident[0] + direct[3] * incident[1],
        ];
        let before = region.stored_energy_j();
        let record = region.push(&incident).unwrap();
        let after = region.stored_energy_j();
        assert_eq!(record.physical_heat, 0.0);
        assert_eq!(record.external_work, 0.0);
        close(record.delta_h, after - before, 2.0e-8, 2.0e-16);
        close(record.delta_h, record.boundary_work, 2.0e-8, 2.0e-16);
        close(
            record.incident_energy - record.outgoing_energy,
            record.boundary_work,
            2.0e-8,
            2.0e-16,
        );
        for mode in 0..2 {
            close(region.incoming()[mode], predicted[mode], 2.0e-8, 2.0e-7);
            let port = region.port_accounts()[mode];
            let z0 = region.source().inlet_reference_impedance();
            let expected_incident = region.duration() * incident[mode] * incident[mode] / z0;
            let squared_mean = region.duration() * region.incoming()[mode].powi(2) / z0;
            close(port.incident_energy, expected_incident, 1.0e-12, 1.0e-20);
            assert!(port.outgoing_energy >= squared_mean * (1.0 - 2.0e-8) - 1.0e-18);
            variation_live |= port.outgoing_energy - squared_mean > 1.0e-12;
        }
        close(
            region
                .port_accounts()
                .iter()
                .map(|p| p.outgoing_energy)
                .sum(),
            record.outgoing_energy,
            2.0e-12,
            1.0e-20,
        );
        incident_total += record.incident_energy;
        outgoing_total += record.outgoing_energy;
    }
    assert!(
        cross_mode_live,
        "geometry must produce current-window cross-mode response"
    );
    assert!(
        variation_live,
        "integrated reflected energy must not collapse to mean-wave squared"
    );
    close(
        region.stored_energy_j() + outgoing_total,
        incident_total,
        3.0e-8,
        1.0e-14,
    );
}

#[test]
fn exact_held_wave_is_invariant_under_three_clock_subdivisions() {
    let duration = 1.0 / 24_000.0;
    let incident = [70.0, -55.0];
    let mut results = Vec::new();
    for subdivisions in [1, 2, 4] {
        let mut region = line(duration / subdivisions as f64);
        let mut outgoing = [0.0; 2];
        let mut means = [0.0; 2];
        for _ in 0..subdivisions {
            region.push(&incident).unwrap();
            for port in 0..2 {
                outgoing[port] += region.port_accounts()[port].outgoing_energy;
                means[port] += region.incoming()[port] / subdivisions as f64;
            }
        }
        results.push((region.state().to_vec(), outgoing, means));
    }
    let metric = line(duration);
    for other in &results[1..] {
        assert!(state_distance(&metric, &results[0].0, &other.0) < 1.0e-9);
        for port in 0..2 {
            close(results[0].1[port], other.1[port], 3.0e-8, 1.0e-15);
            close(results[0].2[port], other.2[port], 3.0e-8, 1.0e-6);
        }
    }
}

#[test]
fn held_waveform_clock_ladder_resolves_nonconstant_input_without_pole_retuning() {
    // Midpoint sampling only defines the piecewise-held incident waveform.
    // Every interval's coupled system evolution remains exact ZOH.
    let duration = 2.0e-4;
    let mut states = Vec::new();
    for intervals in [32, 64, 128] {
        let dt = duration / intervals as f64;
        let mut region = line(dt);
        for i in 0..intervals {
            let t = (i as f64 + 0.5) * dt;
            region
                .push(&[
                    80.0 * det::cos(TAU * 1000.0 * t),
                    60.0 * det::sin(TAU * 6000.0 * t),
                ])
                .unwrap();
        }
        states.push(region.state().to_vec());
    }
    let metric = line(duration / 128.0);
    let coarse = state_distance(&metric, &states[0], &states[1]);
    let fine = state_distance(&metric, &states[1], &states[2]);
    assert!(coarse > 1.0e-6, "clock changes must be live: {coarse:e}");
    assert!(
        fine < 0.35 * coarse && fine > 0.15 * coarse,
        "coarse={coarse:e}, fine={fine:e}"
    );
}

#[test]
fn malformed_and_overflowed_wave_candidates_preserve_all_accepted_history() {
    let mut region = line(1.0 / 96_000.0);
    region.push(&[65.0, -40.0]).unwrap();
    let before = snapshot(&region);
    for input in [&[1.0][..], &[f64::NAN, 0.0][..], &[f64::MAX, 0.0][..]] {
        assert!(region.push(input).is_err());
        assert_eq!(snapshot(&region), before);
    }
    let mut free = [17.0; 2];
    let mut malformed = [23.0; 3];
    assert!(
        region
            .reflection_relation_into(&mut free, &mut malformed)
            .is_err()
    );
    assert_eq!(free, [17.0; 2]);
    assert_eq!(malformed, [23.0; 3]);
    assert_eq!(snapshot(&region), before);
    // Refusals did not poison numerical scratch or the next valid transaction.
    region.push(&[-15.0, 20.0]).unwrap();
}
