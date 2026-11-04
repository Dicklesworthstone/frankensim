//! G1/G3: physical transient Stokes flow, native ports and source admission.

use super::*;
use crate::graph::tests::with_cx;
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};
use fs_qty::MassFlowRate;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap()
}

fn spec() -> LaminarPipeSpec {
    LaminarPipeSpec {
        length: Length::new(0.4),
        radius: Length::new(0.001),
        maximum_reynolds: 1000.0,
        band_rad_s: [0.01, 900.0],
        maximum_relative_wall_error: 1e-4,
    }
}

fn make(cx: &Cx<'_>, mu: f64) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(cx, fluid(mu), spec(), &ZkDiffusiveBasis::new(16).unwrap()).unwrap()
}

fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}, scale {scale:e}"
    );
}

// Independent start-up Stokes solution, using the first sixteen J0 zeros.
// At t/tau >= 0.01 the omitted flow tail is below 1e-16. The integrated
// expression uses the exact Rayleigh sum sum(j_0n^-6)=1/192, so its early
// history does not inherit a truncated steady-flow or acceleration moment.
fn startup(t_over_tau: f64) -> (f64, f64, f64) {
    const ZEROS: [f64; 16] = [
        2.404_825_557_695_772_4,
        5.520_078_110_286_311,
        8.653_727_912_911_013,
        11.791_534_439_014_281,
        14.930_917_708_487_787,
        18.071_063_967_910_924,
        21.211_636_629_879_258,
        24.352_471_530_749_302,
        27.493_479_132_040_253,
        30.634_606_468_431_976,
        33.775_820_213_573_57,
        36.917_098_353_664_045,
        40.058_425_764_628_24,
        43.199_791_713_176_73,
        46.341_188_371_661_815,
        49.482_609_897_397_815,
    ];
    let decay: f64 = ZEROS
        .iter()
        .map(|z| (-z * z * t_over_tau).exp() / z.powi(4))
        .sum();
    let integrated_decay: f64 = ZEROS
        .iter()
        .map(|z| (-z * z * t_over_tau).exp() / z.powi(6))
        .sum();
    // Orthogonality of the radial J0 modes gives the actual integral of u².
    let energy_decay: f64 = ZEROS
        .iter()
        .map(|z| {
            let e = (-z * z * t_over_tau).exp();
            (2.0 * e - e * e) / z.powi(6)
        })
        .sum();
    (
        1.0 - 32.0 * decay,
        1.0 - 1.0 / (6.0 * t_over_tau) + 32.0 * integrated_decay / t_over_tau,
        4.0 / 3.0 - 256.0 * energy_decay,
    )
}

#[test]
fn laminar_pipe_retains_source_and_parabolic_profile_storage_at_steady_flow() {
    struct Counted {
        calls: std::cell::Cell<usize>,
        transport: FluidTransportState,
    }
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.transport)
        }
    }
    with_cx(|cx| {
        let source = Counted {
            calls: std::cell::Cell::new(0),
            transport: fluid(0.001),
        };
        let pipe =
            LaminarPipe::try_new(cx, source, spec(), &ZkDiffusiveBasis::new(16).unwrap()).unwrap();
        let q = 1e-7;
        let mut state = pipe.steady_state(VolumetricFlowRate::new(q), 3.0).unwrap();
        let poiseuille = pipe
            .steady_law()
            .evaluate(&pipe.transport(), MassFlowRate::new(1000.0 * q))
            .unwrap();
        let expected_energy =
            (2.0 / 3.0) * 1000.0 * 0.4 * q * q / (core::f64::consts::PI * 0.001_f64.powi(2));
        near(
            state.kinetic_energy_j(),
            expected_energy,
            expected_energy,
            3e-12,
        );
        near(
            pipe.steady_resistance_pa_s_m3() * q,
            poiseuille.pressure_drop_pa,
            poiseuille.pressure_drop_pa,
            3e-15,
        );
        let report = pipe
            .prepare(cx, 0.02)
            .unwrap()
            .advance(
                cx,
                &mut state,
                [
                    Pressure::new(poiseuille.pressure_drop_pa),
                    Pressure::new(0.0),
                ],
            )
            .unwrap();
        near(report.final_volume_flow_m3_s, q, q, 2e-11);
        near(report.average_volume_flow_m3_s, q, q, 2e-11);
        near(
            report.native.dissipated,
            0.02 * poiseuille.dissipated_power_w,
            0.02 * poiseuille.dissipated_power_w,
            2e-11,
        );
        near(
            report.native.balance_residual(),
            0.0,
            expected_energy,
            2e-12,
        );
        assert_eq!(pipe.source().calls.get(), 1);
        assert_eq!(state.accepted_time_s(), 3.02);
    });
}

#[test]
fn laminar_pipe_startup_and_reversal_match_independent_radial_stokes_series() {
    with_cx(|cx| {
        let pipe = make(cx, 0.001);
        let q_infinity = 100.0 * core::f64::consts::PI * 0.001_f64.powi(4) / (8.0 * 0.001 * 0.4);
        for duration in [0.01, 0.05, 0.2, 1.0] {
            let (endpoint, average, profile_energy) = startup(duration); // tau=rho*r²/mu=1s.
            let mut prepared = pipe.prepare(cx, duration).unwrap();
            let mut rest = pipe
                .steady_state(VolumetricFlowRate::new(0.0), 0.0)
                .unwrap();
            let (history, mobility) = prepared.average_flow_relation(&rest).unwrap();
            assert_eq!(history, 0.0);
            let report = prepared
                .advance(cx, &mut rest, [Pressure::new(100.0), Pressure::new(0.0)])
                .unwrap();
            near(
                report.final_volume_flow_m3_s,
                q_infinity * endpoint,
                q_infinity,
                2e-9,
            );
            near(
                report.average_volume_flow_m3_s,
                q_infinity * average,
                q_infinity,
                2e-9,
            );
            let plug_energy = 0.5 * pipe.mean_inertance_pa_s2_m3() * q_infinity * q_infinity;
            near(
                rest.kinetic_energy_j(),
                plug_energy * profile_energy,
                plug_energy,
                2e-9,
            );
            near(
                report.average_volume_flow_m3_s,
                100.0 * mobility,
                q_infinity,
                3e-13,
            );
            near(
                report.displaced_volume_m3,
                duration * q_infinity * average,
                duration * q_infinity,
                2e-9,
            );
            assert!(report.native.dissipated > 0.0 && report.native.delta_h > 0.0);
            near(
                report.native.supplied,
                100.0 * report.displaced_volume_m3,
                report.native.supplied,
                3e-13,
            );
            assert!(report.native.dissipated < report.native.supplied);

            let mut reverse = pipe
                .steady_state(VolumetricFlowRate::new(q_infinity), 0.0)
                .unwrap();
            let report = prepared
                .advance(
                    cx,
                    &mut reverse,
                    [Pressure::new(-100.0), Pressure::new(0.0)],
                )
                .unwrap();
            near(
                report.final_volume_flow_m3_s,
                q_infinity * (1.0 - 2.0 * endpoint),
                q_infinity,
                4e-9,
            );
            near(
                report.average_volume_flow_m3_s,
                q_infinity * (1.0 - 2.0 * average),
                q_infinity,
                4e-9,
            );
            assert!(report.native.dissipated > 0.0);
        }
        // A quasi-steady resistance plus plug inertance has the wrong decay.
        assert!((1.0 - (-8.0_f64 * 0.2).exp() - startup(0.2).0).abs() > 0.09);
    });
}

#[test]
fn laminar_pipe_native_harmonic_ports_resolve_womersley_phase_and_dissipation() {
    with_cx(|cx| {
        for mu in [0.001, 0.004] {
            let pipe = make(cx, mu);
            let model = HarmonicModel::from_ode(
                pipe.system(),
                HarmonicBudget {
                    max_dimension: 32,
                    max_dense_bytes: 1_000_000,
                    max_dense_work: 1_000_000,
                    relative_residual_tolerance: 1e-11,
                },
                || false,
            )
            .unwrap();
            for omega in [0.01, 0.1, 1.0, 10.0, 100.0, 900.0] {
                let actual = model
                    .factor(omega, || false)
                    .unwrap()
                    .response(&[C64::from_re(100.0), C64::ZERO], || false)
                    .unwrap();
                let rv = (omega * 1000.0 * 0.001_f64.powi(2) / mu).sqrt();
                // Native harmonic uses exp(-i wt); wall oracle uses exp(+i wt).
                let exact_z = (C64::new(0.0, omega) / fs_phs::zk_one_minus_f(rv).unwrap())
                    .scale(pipe.mean_inertance_pa_s2_m3())
                    .conj();
                let expected_flow = C64::from_re(100.0) / exact_z;
                assert!((actual.output[0] - expected_flow).abs() / expected_flow.abs() < 1e-4);
                let pole_flow = C64::from_re(100.0) / pipe.impedance_at(omega).unwrap().conj();
                assert!((actual.output[0] - pole_flow).abs() / pole_flow.abs() < 2e-10);
                assert_eq!(actual.output[1], -actual.output[0]);
                near(
                    actual.dissipated_power_w,
                    50.0 * actual.output[0].re,
                    actual.dissipated_power_w,
                    2e-10,
                );
                assert!(actual.dissipated_power_w > 0.0);
            }
        }
    });
}

#[test]
fn laminar_pipe_composes_with_finite_pressure_storage_on_original_phs() {
    with_cx(|cx| {
        let pipe = make(cx, 0.001);
        let mut initial = pipe
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap()
            .coordinates()
            .to_vec();
        let capacity: f64 = 1e-8; // dV/dp [m³/Pa], authored fixed compliance.
        let pressure = 100.0;
        initial.push(capacity.sqrt() * pressure);
        let reservoir = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![0.0],
            vec![-1.0 / capacity.sqrt()],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap();
        let (_original_source, native_pipe) = pipe.into_native();
        let joint = fs_phs::interconnect(native_pipe, reservoir, &[(0, 0)]).unwrap();
        let prepared = PreparedExactZoh::new_owned(
            joint,
            0.05,
            LinearSplitBudget {
                max_modes: 32,
                max_omega_dt: 1.0,
            },
        )
        .unwrap();
        let report = prepared.step(&initial, &[0.0]).unwrap(); // outlet held at zero gauge.
        let volume_change = capacity.sqrt() * (report.x.last().unwrap() - initial.last().unwrap());
        near(
            volume_change,
            0.05 * report.y[0],
            capacity * pressure,
            2e-11,
        );
        assert!(volume_change < 0.0 && report.dissipated > 0.0);
        assert!(report.x[0] > 0.0);
        assert!(report.x[1..17].iter().all(|value| *value > 0.0));
        assert_eq!(report.supplied, 0.0);
        near(
            report.delta_h,
            -report.dissipated,
            0.5 * capacity * pressure * pressure,
            2e-12,
        );
    });
}

#[test]
fn laminar_pipe_history_regime_and_pressure_gauge_are_admitted_atomically() {
    with_cx(|cx| {
        let pipe = make(cx, 0.001);
        let other = make(cx, 0.001);
        let mut prepared = pipe.prepare(cx, 0.05).unwrap();
        let mut foreign = other
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        assert!(
            prepared
                .advance(cx, &mut foreign, [Pressure::new(4.0), Pressure::new(0.0)])
                .is_err()
        );
        let mut state = pipe
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let original = state.clone();
        assert!(
            prepared
                .advance(cx, &mut state, [Pressure::new(1e8), Pressure::new(0.0)])
                .is_err()
        );
        assert_eq!(state.coordinates(), original.coordinates());
        assert_eq!(state.accepted_time_s(), original.accepted_time_s());
        assert!(
            prepared
                .advance(
                    cx,
                    &mut state,
                    [Pressure::new(f64::NAN), Pressure::new(0.0)]
                )
                .is_err()
        );
        assert_eq!(state.coordinates(), original.coordinates());
        let mut shifted = original;
        let first = prepared
            .advance(cx, &mut state, [Pressure::new(7.0), Pressure::new(3.0)])
            .unwrap();
        let second = prepared
            .advance(
                cx,
                &mut shifted,
                [Pressure::new(100_007.0), Pressure::new(100_003.0)],
            )
            .unwrap();
        assert_eq!(state.coordinates(), shifted.coordinates());
        assert_eq!(first.native, second.native);
        assert!(pipe.impedance_at(901.0).is_err());
        assert!(
            LaminarPipe::try_new(cx, fluid(0.001), spec(), &ZkDiffusiveBasis::new(2).unwrap())
                .is_err()
        );
    });
}
