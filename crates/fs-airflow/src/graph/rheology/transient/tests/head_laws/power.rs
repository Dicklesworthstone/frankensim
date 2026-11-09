//! Accepted original power values share the native head query and transaction.

use super::super::super::head::power::{FlowHeadPowerReceipt, PumpLossRouting, PumpPowerState};
use super::*;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug)]
struct Powered {
    original: Arc<Cubic>,
    queries: AtomicU64,
    corrupt_partition: bool,
}
impl FlowHeadLaw for Powered {
    fn flow_interval_m3_s(&self) -> [f64; 2] {
        self.original.support
    }
    fn evaluate(
        &self,
        flow: VolumetricFlowRate,
    ) -> Result<FlowHeadResponse, Box<dyn std::error::Error + Send + Sync>> {
        self.queries.fetch_add(1, Ordering::Relaxed);
        let mut response = self.original.evaluate(flow)?;
        let shaft = 5.0 + 2000.0 * flow.value();
        let pump_loss = shaft - response.pressure.value() * flow.value();
        response.power = Some(FlowHeadPowerReceipt {
            state: PumpPowerState {
                shaft_input_w: shaft,
                electrical_input_w: shaft + 2.0,
                pump_loss_w: pump_loss,
                drive_loss_w: 2.0,
                reservoir_heat_w: 0.25 * pump_loss
                    + 1.0
                    + if self.corrupt_partition { 1.0 } else { 0.0 },
                surroundings_heat_w: 0.75 * pump_loss + 1.0,
            },
            routing: PumpLossRouting {
                pump_to_reservoir: 0.25,
                drive_to_reservoir: 0.5,
            },
        });
        Ok(response)
    }
}
fn powered(corrupt_partition: bool) -> Arc<Powered> {
    Arc::new(Powered {
        original: Arc::new(Cubic {
            shutoff: 800.0,
            linear: 2e5,
            cubic: 4e10,
            support: [0.0, 0.001],
        }),
        queries: AtomicU64::new(0),
        corrupt_partition,
    })
}

#[test]
fn g1_shutoff_and_flow_retain_actual_input_losses_without_an_unmetered_second_query() {
    with_cx(|cx| {
        let model = single();
        for base in [-800.0, 0.0] {
            let original = powered(false);
            let source = FlowHeadSource::new(original.clone());
            let mut state = model.rest_state(&[None; 2], cx).unwrap();
            let duration = 0.02;
            let report = model
                .step_with_head_laws(
                    &mut state,
                    &mut model.workspace().unwrap(),
                    &[boundary(0, base), boundary(1, 0.0)],
                    &[],
                    None,
                    &[Some(source.clone())],
                    duration,
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                report.head_evaluations,
                original.queries.load(Ordering::Relaxed)
            );
            let accepted = report.flow_heads[0].as_ref().unwrap();
            assert_eq!(accepted.source, source);
            let q = accepted.midpoint_flow.value();
            let receipt = accepted.response.power.unwrap();
            let shaft = 5.0 + 2000.0 * q;
            let hydraulic = original.original.pressure(q) * q;
            near(receipt.state.shaft_input_w, shaft, 1e-13, 1e-13);
            near(receipt.state.electrical_input_w, shaft + 2.0, 1e-13, 1e-13);
            near(
                receipt.state.reservoir_heat_w,
                0.25 * (shaft - hydraulic) + 1.0,
                1e-13,
                1e-13,
            );
            near(
                receipt.state.surroundings_heat_w,
                0.75 * (shaft - hydraulic) + 1.0,
                1e-13,
                1e-13,
            );
            near(
                receipt.state.electrical_input_w,
                hydraulic + receipt.state.reservoir_heat_w + receipt.state.surroundings_heat_w,
                1e-13,
                1e-13,
            );
            if base == -800.0 {
                assert_eq!(q, 0.0);
                assert_eq!(report.viscous_heat_j, 0.0);
                assert_eq!(accepted.work_j, 0.0);
                assert_eq!(receipt.state.electrical_input_w, 7.0);
                assert_eq!(receipt.state.reservoir_heat_w, 2.25);
            } else {
                let expected = independent([0.0; 2], duration, base, &original.original);
                near(q, expected.1, 1e-9, 1e-15);
                near(report.viscous_heat_j, expected.2, 1e-9, 1e-12);
                near(report.kinetic_energy_change_j, expected.3, 1e-9, 1e-12);
                near(accepted.work_j, hydraulic * duration, 1e-9, 1e-12);
            }
        }
    });
}

#[test]
fn g4_inconsistent_original_power_partition_refuses_before_any_native_publication() {
    with_cx(|cx| {
        let model = single();
        let mut state = model.rest_state(&[None; 2], cx).unwrap();
        let before = state.branches()[0].velocities_m_s().to_vec();
        let mut workspace = model.workspace().unwrap();
        let bad = powered(true);
        let result = model.step_with_head_laws(
            &mut state,
            &mut workspace,
            &[boundary(0, 0.0), boundary(1, 0.0)],
            &[],
            None,
            &[Some(FlowHeadSource::new(bad.clone()))],
            0.02,
            controls(),
            cx,
        );
        assert!(matches!(
            result,
            Err(TransientCarreauError::Head {
                source: FlowHeadError::Invalid(_),
                ..
            })
        ));
        assert_eq!(bad.queries.load(Ordering::Relaxed), 1);
        assert_eq!(state.accepted_time_s(), 0.0);
        assert_eq!(state.branches()[0].velocities_m_s(), before);
        model
            .step_with_head_laws(
                &mut state,
                &mut workspace,
                &[boundary(0, 0.0), boundary(1, 0.0)],
                &[],
                None,
                &[Some(FlowHeadSource::new(powered(false)))],
                0.02,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(state.accepted_time_s(), 0.02);
    });
}
