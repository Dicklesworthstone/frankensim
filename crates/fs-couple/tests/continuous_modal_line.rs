//! Exact coupled reflection, original work, transactionality and an independent
//! three-rung clock axis. Spatial/radial authority lives in fs-duct tests.
use core::f64::consts::TAU;
use fs_couple::mm_line::continuous::{ContinuousModalLine, ModalAdmittanceLoad};
use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::ModalEnd;
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PortHamiltonian, QuadraticStorage};

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

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 124,
        max_omega_dt: 1.0,
    }
}

fn harmonic(system: &PortHamiltonian, omega: f64) -> Vec<C64> {
    let model = HarmonicModel::from_ode(
        system,
        HarmonicBudget {
            max_dimension: 124,
            max_dense_bytes: 16 * 1024 * 1024,
            max_dense_work: 4_000_000,
            relative_residual_tolerance: 2.0e-10,
        },
        || false,
    )
    .unwrap();
    model
        .factor(omega, || false)
        .unwrap()
        .response(&[C64::ONE], || false)
        .unwrap()
        .output
}

fn coupled_load(radius: f64, gas: GasState, coupling: &[f64]) -> ModalAdmittanceLoad {
    // Authored complete physical acoustic oscillator, not a fit or material
    // claim: volume displacement, acoustic momentum, and original port maps.
    let modes = coupling.len();
    let mut g = vec![0.0; 2 * modes];
    g[modes..].copy_from_slice(coupling);
    let system = PortHamiltonian::new(
        2,
        modes,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 12_000.0],
        g,
        Box::new(QuadraticStorage::new(vec![1.8e8, 0.0, 0.0, 1.0 / 37.0], 2).unwrap()),
    )
    .unwrap();
    ModalAdmittanceLoad::new(system, radius, ModeCount(modes), gas).unwrap()
}

#[test]
fn bidirectional_modal_pulses_keep_each_face_reference_and_all_wave_work() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let profile = RadiusProfile::exponential(0.025, 8.0, 0.08).unwrap();
    let mut region = ContinuousModalLine::new(
        &profile,
        &gas,
        ModeCount(2),
        4,
        ModalEnd::Ports,
        1.0 / 96_000.0,
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(region.port_count(), 4);
    assert_eq!(region.state().len(), 20);
    assert!(region.load_state().is_empty());
    let zi = region.source().inlet_reference_impedance();
    let zm = region.source().mouth_reference_impedance();
    assert!((zi / zm - 1.0).abs() > 0.5);
    let (mut incoming, mut outgoing) = (0.0, 0.0);
    for step in 0..120 {
        let a = if step < 45 {
            [50.0, -30.0, 35.0, 45.0]
        } else {
            [0.0; 4]
        };
        let mut free = [0.0; 4];
        let mut direct = [0.0; 16];
        region
            .reflection_relation_into(&mut free, &mut direct)
            .unwrap();
        let record = region.push(&a).unwrap();
        for port in 0..4 {
            let predicted = free[port] + (0..4).map(|j| direct[port * 4 + j] * a[j]).sum::<f64>();
            close(region.incoming()[port], predicted, 3.0e-8, 1.0e-6);
            let z = if port < 2 { zi } else { zm };
            close(
                region.port_accounts()[port].incident_energy,
                region.duration() * a[port] * a[port] / z,
                1.0e-12,
                1.0e-20,
            );
        }
        assert_eq!(record.physical_heat, 0.0);
        assert_eq!(record.external_work, 0.0);
        incoming += record.incident_energy;
        outgoing += record.outgoing_energy;
    }
    close(
        region.stored_energy_j() + outgoing,
        incoming,
        3.0e-8,
        1.0e-14,
    );
}

#[test]
fn original_dynamic_load_impedance_converges_under_three_spatial_rungs() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let (radius, length) = (0.025, 0.12);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    for hz in [200.0, 600.0] {
        let omega = TAU * hz;
        let phase = omega * length / gas.sound_speed;
        // Independent acoustic network law and the analytic cylinder's TMM.
        let zl = C64::new(12_000.0, 1.8e8 / omega - 37.0 * omega);
        let c = C64::from_re(det::cos(phase));
        let s = C64::new(0.0, -det::sin(phase));
        let reference = C64::from_re(z0);
        let expected = (zl * c + s * reference) / (c + s * zl / reference);
        let mut errors = Vec::new();
        for elements in [12, 24, 48] {
            let load = coupled_load(radius, gas, &[1.0]);
            let region = ContinuousModalLine::new_loaded(
                &profile,
                &gas,
                ModeCount(1),
                elements,
                load,
                1.0 / 96_000.0,
                budget(),
                || false,
            )
            .unwrap();
            assert_eq!(region.state().len(), 2 * (elements + 1) + 2);
            assert_eq!(region.load_state().len(), 2);
            assert_eq!(region.system().port_dim(), 1);
            errors.push((harmonic(region.system(), omega)[0] - expected).abs() / z0);
        }
        assert!(errors[0] > 1.0e-7, "live spatial defect: {errors:?}");
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.35 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.004, "hz={hz}, errors={errors:?}");
    }
}

#[test]
fn full_coupled_load_memory_heat_and_reflection_survive_exact_clock_subdivision() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let profile = RadiusProfile::exponential(0.025, 0.0, 0.04).unwrap();
    // Longer than the round trip, so the reflected higher mode actually
    // interrogates the coupled termination instead of only the inlet.
    let duration = 1.0 / 2400.0;
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.025_f64.powi(2));
    let mut results = Vec::new();
    for subdivision in [1, 2, 4] {
        let load = coupled_load(0.025, gas, &[1.0, 0.6])
            .with_direct_conductance(vec![0.15 / z0, 0.08 / z0, 0.08 / z0, 0.2 / z0])
            .unwrap();
        let mut region = ContinuousModalLine::new_loaded(
            &profile,
            &gas,
            ModeCount(2),
            3,
            load,
            duration / subdivision as f64,
            budget(),
            || false,
        )
        .unwrap();
        let (mut outgoing, mut heat, mut incoming) = (0.0, 0.0, 0.0);
        for _ in 0..subdivision {
            let record = region.push(&[75.0, 0.0]).unwrap();
            outgoing += record.outgoing_energy;
            heat += record.physical_heat;
            incoming += record.incident_energy;
        }
        assert!(heat > 1.0e-13, "actual load resistance must dissipate");
        assert!(region.load_state().iter().all(|x| x.abs() > 1.0e-14));
        assert!(
            region.port_accounts()[1].outgoing_energy > 1.0e-14,
            "uniform duct's cross-mode reflection must come from original mutual load"
        );
        close(
            region.stored_energy_j() + heat + outgoing,
            incoming,
            3.0e-8,
            1.0e-14,
        );
        let before = snapshot(&region);
        assert!(region.push(&[f64::MAX, 0.0]).is_err());
        assert_eq!(snapshot(&region), before);
        results.push((region.state().to_vec(), outgoing, heat));
    }
    for other in &results[1..] {
        for (&a, &b) in results[0].0.iter().zip(&other.0) {
            close(a, b, 2.0e-7, 1.0e-12);
        }
        close(results[0].1, other.1, 3.0e-8, 1.0e-15);
        close(results[0].2, other.2, 3.0e-8, 1.0e-15);
    }
}

#[test]
fn original_direct_load_has_no_fake_state_and_preserves_matching_and_refusals() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let radius = 0.025;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.04).unwrap();
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let direct =
        ModalAdmittanceLoad::from_conductance(vec![1.0 / z0], radius, ModeCount(1), gas).unwrap();
    let mut region = ContinuousModalLine::new_loaded(
        &profile,
        &gas,
        ModeCount(1),
        8,
        direct,
        1.0 / 96_000.0,
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(region.state().len(), 18);
    assert!(region.load_state().is_empty());
    assert_eq!(region.port_count(), 1);
    // The original two physical source slots remain; only inlet is incident.
    assert_eq!(region.system().port_dim(), 2);
    let model = HarmonicModel::from_ode(
        region.system(),
        HarmonicBudget {
            max_dimension: 124,
            max_dense_bytes: 16 * 1024 * 1024,
            max_dense_work: 4_000_000,
            relative_residual_tolerance: 2.0e-10,
        },
        || false,
    )
    .unwrap();
    let output = model
        .factor(TAU * 600.0, || false)
        .unwrap()
        .response(&[C64::ONE, C64::ZERO], || false)
        .unwrap()
        .output;
    assert!((output[0] - C64::from_re(z0)).abs() / z0 < 0.001);
    let (mut heat, mut incident, mut reflected) = (0.0, 0.0, 0.0);
    for _ in 0..64 {
        let frame = region.push(&[40.0]).unwrap();
        heat += frame.physical_heat;
        incident += frame.incident_energy;
        reflected += frame.outgoing_energy;
        assert_eq!(frame.external_work, 0.0);
    }
    assert!(heat > 1.0e-7);
    close(
        region.stored_energy_j() + heat + reflected,
        incident,
        3.0e-8,
        1.0e-14,
    );
    let wrong =
        ModalAdmittanceLoad::from_conductance(vec![1.0 / z0], 2.0 * radius, ModeCount(1), gas)
            .unwrap();
    assert!(
        ContinuousModalLine::new_loaded(
            &profile,
            &gas,
            ModeCount(1),
            8,
            wrong,
            1.0 / 96_000.0,
            budget(),
            || false
        )
        .is_err()
    );
    assert!(
        ModalAdmittanceLoad::from_conductance(vec![-1.0 / z0], radius, ModeCount(1), gas).is_err()
    );
    assert!(
        ContinuousModalLine::new_loaded(
            &profile,
            &gas,
            ModeCount(2),
            3,
            coupled_load(radius, gas, &[1.0, 0.6]),
            1.0 / 96_000.0,
            LinearSplitBudget {
                max_modes: 16,
                max_omega_dt: 1.0
            },
            || false
        )
        .is_err()
    );
}
