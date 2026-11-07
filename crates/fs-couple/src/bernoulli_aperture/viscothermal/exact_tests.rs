//! Actual viscous/thermal state owners and the graph using their mean ports.
use super::*;
use fs_material::gas::GasSpec;

const DT: f64 = 1.0 / 96_000.0;
const BAND: [f64; 2] = [50.0, 1500.0];

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn selection() -> ViscothermalSectionSpec {
    ViscothermalSectionSpec {
        section: 0,
        cells: 16,
        minimum_frequency_hz: BAND[0],
        maximum_frequency_hz: BAND[1],
    }
}
fn graph(gas: GasState, fractional: bool) -> TubeNetworkSpec {
    TubeNetworkSpec {
        nodes: vec![
            NetworkNode::Inlet,
            NetworkNode::Termination { reflection: 0.0 },
        ],
        sections: vec![TubeSection {
            nodes: [0, 1],
            length_m: 0.25,
            radius_m: 0.007,
            max_length_error_m: if fractional { 1e-12 } else { 0.002 },
        }],
        sound_speed_m_s: gas.sound_speed,
        max_wave_memory_bytes: 8 << 20,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn sampled_fit_report_is_the_actual_selected_viscous_and_thermal_transfer() {
    let gas = gas();
    let gate = CancelGate::new_clock_free();
    let midpoint = WideTubeLoss::new(0.007, &gas, DT, BAND, &gate).unwrap();
    let exact = WideTubeLoss::new_with_sampling(
        0.007,
        &gas,
        DT,
        BAND,
        RelaxationSampling::ExactHeldFlow,
        &gate,
    )
    .unwrap();
    assert_eq!(midpoint.terms, exact.terms); // Same physical Foster identification.
    // The shared analogue fitting error can dominate either aggregate report.
    // Compare the selected sampled response itself, then reconstruct both
    // complete report maxima below.
    let midpoint_z = midpoint.excess_at(BAND[1], true).unwrap().0;
    let exact_z = exact.excess_at(BAND[1], true).unwrap().0;
    assert!((midpoint_z - exact_z).abs() > 1e-6 * midpoint_z.abs());
    let length = 0.02;
    for loss in [&midpoint, &exact] {
        let series = loss.series(length).unwrap();
        let shunt = loss.shunt(length).unwrap();
        assert_eq!(series.sampling(), loss.sampling());
        assert_eq!(shunt.sampling(), loss.sampling());
        assert_eq!(series.base().inertance_pa_s2_m3, 0.0);
        assert_eq!(shunt.compliance_m3_pa(), 0.0);
        let (mut complex_error, mut real_error) = (0.0_f64, 0.0_f64);
        for i in 0..=64 {
            let f = BAND[0] * det::pow(BAND[1] / BAND[0], f64::from(i) / 64.0);
            let target = C64::new(
                det::sqrt(core::f64::consts::TAU * f),
                -det::sqrt(core::f64::consts::TAU * f),
            );
            let (z, y) = loss.excess_at(f, true).unwrap();
            close(
                series.sampled_impedance(f, DT).unwrap().conj(),
                z.scale(length),
                2e-13,
            );
            close(
                shunt.sampled_admittance(f, DT).unwrap().conj(),
                y.scale(length),
                2e-13,
            );
            for discrete in [false, true] {
                let value = loss
                    .excess_at(f, discrete)
                    .unwrap()
                    .0
                    .scale(1.0 / loss.series_gain);
                complex_error = complex_error.max((value - target).abs() / target.abs());
                real_error = real_error.max((value.re - target.re).abs() / target.re);
            }
        }
        assert!((loss.max_complex_error() - complex_error).abs() < 2e-14);
        assert!((loss.max_real_loss_error() - real_error).abs() < 2e-14);
        assert!(complex_error <= 0.05 && real_error <= 0.05);
    }
}

#[test]
fn exact_distributed_loss_matches_actual_scattering_and_closes_network_work() {
    let gas = gas();
    let gate = CancelGate::new_clock_free();
    let frequency = 500.0;
    for fractional in [false, true] {
        let policy = fractional.then_some(FractionalTubePolicy {
            maximum_frequency_hz: BAND[1],
            maximum_summed_phase_error_rad: 0.05,
        });
        let options = ViscothermalTimeOptions {
            sampling: RelaxationSampling::ExactHeldFlow,
            fractional: policy,
        };
        let (lowered, gases, reports) = with_regional_viscothermal_sections_options(
            graph(gas, fractional),
            vec![gas],
            DT,
            &[selection()],
            options,
            &gate,
        )
        .unwrap();
        assert_eq!(
            reports[0].loss.sampling(),
            RelaxationSampling::ExactHeldFlow
        );
        assert!(reports[0].max_scattering_error <= 0.03);
        let (mut network, represented, details) = if let Some(policy) = policy {
            let (network, represented, details) = lowered
                .realize_fractional_waveguide(&gases, DT, policy, &gate)
                .unwrap();
            (network, represented, Some(details))
        } else {
            let (network, represented) = lowered.realize_waveguide(&gases, DT, &gate).unwrap();
            (network, represented, None)
        };
        // Compose the declarations actually consumed by the graph, including
        // each real delay/allpass and each selected mean pressure/flow load.
        // Use exp(+i omega t), independently conjugate to the cold duct report.
        let zc = represented[0].impedance_pa_s_m3;
        let mut transfer = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
        for (i, section) in lowered.sections.iter().enumerate() {
            let omega = core::f64::consts::TAU * frequency * DT;
            let lag = omega * represented[i].one_way_samples as f64
                + details
                    .as_ref()
                    .map_or(0.0, |d| d[i].transit.phase_lag_rad(omega).unwrap());
            let (sin, cos) = lag.sin_cos();
            let c = C64::from_re(cos);
            let s = C64::new(0.0, sin);
            transfer = mul(transfer, [c, s, s, c]);
            match lowered.nodes[section.nodes[1]] {
                NetworkNode::Series { load } => {
                    assert_eq!(load.sampling(), RelaxationSampling::ExactHeldFlow);
                    transfer = mul(
                        transfer,
                        [
                            C64::ONE,
                            load.sampled_impedance(frequency, DT)
                                .unwrap()
                                .scale(1.0 / zc),
                            C64::ZERO,
                            C64::ONE,
                        ],
                    );
                }
                NetworkNode::ShuntAdmittance { load } => {
                    assert_eq!(load.sampling(), RelaxationSampling::ExactHeldFlow);
                    transfer = mul(
                        transfer,
                        [
                            C64::ONE,
                            C64::ZERO,
                            load.sampled_admittance(frequency, DT).unwrap().scale(zc),
                            C64::ONE,
                        ],
                    );
                }
                NetworkNode::Junction | NetworkNode::Termination { .. } => {}
                _ => panic!("unexpected lowered physical node"),
            }
        }
        let predicted = scattering(transfer)[2];
        let (warm, measured) = (49_152, 192 * 64);
        let (mut output, mut input) = (C64::ZERO, C64::ZERO);
        let (mut work, mut heat) = (0.0, 0.0);
        for sample in 0..warm + measured {
            let (sin, cos) = (core::f64::consts::TAU * f64::from(sample) / 192.0).sin_cos();
            let before = network.stored_energy_j();
            let frame = network.step(cos).unwrap();
            let scale = before
                + frame.stored_energy_j
                + frame.inlet_work_j.abs()
                + frame.terminal_loss_j
                + frame.interior_loss_j;
            assert!(frame.balance_residual_j().abs() <= 3e-11 * scale.max(1e-30));
            assert!(frame.interior_loss_j >= 0.0);
            work += frame.inlet_work_j;
            heat += frame.terminal_loss_j + frame.interior_loss_j;
            if sample >= warm {
                let demodulate = C64::new(cos, -sin);
                output = output + demodulate.scale(network.node_frame(1).unwrap().pressure_pa);
                input = input + demodulate.scale(cos);
            }
        }
        // Remaining initial low-pole transients are much smaller than the
        // predeclared scattering error gate; this is a finite-duration tone.
        close(output / input, predicted, 2e-6);
        assert!(heat > 0.0);
        assert!((network.stored_energy_j() + heat - work).abs() < 3e-10 * work.abs());
    }
}

#[test]
fn selected_time_options_keep_defaults_geometry_and_existing_refusals() {
    let gas = gas();
    let gate = CancelGate::new_clock_free();
    let original = graph(gas, false);
    let (legacy, legacy_gases, legacy_reports) =
        with_regional_viscothermal_sections(original.clone(), vec![gas], DT, &[selection()], &gate)
            .unwrap();
    let (selected, selected_gases, selected_reports) = with_regional_viscothermal_sections_options(
        original.clone(),
        vec![gas],
        DT,
        &[selection()],
        ViscothermalTimeOptions::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(legacy, selected);
    assert_eq!(legacy_gases, selected_gases);
    assert_eq!(
        legacy_reports[0].max_scattering_error,
        selected_reports[0].max_scattering_error
    );
    let options = ViscothermalTimeOptions {
        sampling: RelaxationSampling::ExactHeldFlow,
        fractional: None,
    };
    let (empty, gases, reports) = with_regional_viscothermal_sections_options(
        original.clone(),
        vec![gas],
        DT,
        &[],
        options,
        &gate,
    )
    .unwrap();
    assert_eq!(empty, original);
    assert_eq!(gases, [gas]);
    assert!(reports.is_empty());
    for (radius, dt) in [(0.0001, DT), (0.007, 1.0 / 20_000.0)] {
        assert!(
            WideTubeLoss::new_with_sampling(
                radius,
                &gas,
                dt,
                BAND,
                RelaxationSampling::ExactHeldFlow,
                &gate
            )
            .is_err()
        );
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        with_regional_viscothermal_sections_options(
            original,
            vec![gas],
            DT,
            &[selection()],
            options,
            &cancelled
        )
        .is_err()
    );
}
