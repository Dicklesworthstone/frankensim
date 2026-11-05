//! Coupled trajectory, local irreversible heat and retained viscosity history.
use super::*;
use crate::graph::inertial::laminar::LaminarPipeSpec;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::linear_exact::PreparedCommonFlowZoh;
use fs_phs::{DescriptorPortHamiltonian, ZkDiffusiveBasis};
use fs_qty::Length;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap()
}
fn pipe(cx: &Cx<'_>, length: f64, mu: f64) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        fluid(mu),
        LaminarPipeSpec {
            length: Length::new(length),
            radius: Length::new(0.001),
            maximum_reynolds: 1000.0,
            band_rad_s: [0.01, 900.0],
            maximum_relative_wall_error: 1e-4,
        },
        &ZkDiffusiveBasis::new(16).unwrap(),
    )
    .unwrap()
}
fn near(a: f64, b: f64, scale: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * scale,
        "{a:e} != {b:e}, scale {scale:e}"
    );
}
fn balance(step: &LaminarSeriesStep, heat: &[f64]) {
    let scale = step
        .native
        .supplied
        .abs()
        .max(step.native.dissipated)
        .max(step.initial_mean_kinetic_energy_j + step.initial_profile_memory_energy_j);
    near(
        step.native.delta_h + step.native.dissipated,
        step.native.supplied,
        scale,
        3e-10,
    );
    near(heat.iter().sum(), step.native.dissipated, scale, 3e-10);
    assert!(heat.iter().all(|h| h.is_finite() && *h >= 0.0));
}

#[test]
fn splitting_a_uniform_pipe_preserves_startup_reversal_storage_and_local_heat() {
    with_cx(|cx| {
        let whole = pipe(cx, 0.4, 0.001);
        let split =
            LaminarSeries::try_new(cx, vec![pipe(cx, 0.15, 0.001), pipe(cx, 0.25, 0.001)]).unwrap();
        let mut single = whole
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut path = split
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut one_step = whole.prepare(cx, 0.025).unwrap();
        let mut path_step = split.prepare(cx, 0.025).unwrap();
        let mut heat = [0.0; 2];
        for pressure in [120.0, 120.0, -80.0, -80.0, 0.0, 0.0] {
            let ports = [
                Pressure::new(100_000.0 + pressure),
                Pressure::new(100_000.0),
            ];
            let (history, mobility) = path_step.average_flow_relation(&path).unwrap();
            let original = one_step.advance(cx, &mut single, ports).unwrap();
            let connected = path_step.advance(cx, &mut path, ports, &mut heat).unwrap();
            near(
                path.volume_flow_m3_s(),
                single.volume_flow_m3_s(),
                1e-7,
                2e-9,
            );
            near(
                path.kinetic_energy_j(),
                single.kinetic_energy_j(),
                single.kinetic_energy_j(),
                5e-9,
            );
            near(
                connected.average_volume_flow_m3_s,
                original.average_volume_flow_m3_s,
                1e-7,
                2e-9,
            );
            near(
                connected.average_volume_flow_m3_s,
                history + mobility * pressure,
                1e-7,
                2e-10,
            );
            near(
                heat[0],
                original.native.dissipated * 0.375,
                original.native.dissipated,
                5e-9,
            );
            near(
                heat[1],
                original.native.dissipated * 0.625,
                original.native.dissipated,
                5e-9,
            );
            let section_energy = (0..2)
                .map(|i| split.section_kinetic_energy_j(&path, i).unwrap())
                .sum::<f64>();
            near(
                section_energy,
                path.kinetic_energy_j(),
                path.kinetic_energy_j(),
                1e-14,
            );
            balance(&connected, &heat);
        }
    });
}

// Independent composition retains both original pipe mean coordinates and one
// massless pressure multiplier; the native exact descriptor projection enforces
// Q0=Q1. This does not use the series mean-flow elimination or its R assembly.
fn descriptor_reference(
    sections: &[LaminarPipe<FluidTransportState>],
    dt: f64,
) -> PreparedCommonFlowZoh {
    assert_eq!(sections.len(), 2);
    let dims = sections
        .iter()
        .map(|p| p.system().state_dim())
        .collect::<Vec<_>>();
    let diff = dims.iter().sum::<usize>();
    let n = diff + 1;
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; 2 * n];
    let mut q = vec![0.0; diff * diff];
    let mut start = 0;
    for (i, section) in sections.iter().enumerate() {
        let (_, local, _) = section.system().structure();
        for row in 0..dims[i] {
            q[(start + row) * diff + start + row] = 1.0;
            for column in 0..dims[i] {
                r[(start + row) * n + start + column] = local[row * dims[i] + column];
            }
        }
        let boundary = if i == 0 { 1.0 } else { -1.0 } / section.flow_scales[0];
        g[start * 2 + i] = boundary;
        j[start * n + diff] = -boundary;
        j[diff * n + start] = boundary;
        start += dims[i];
    }
    let descriptor = DescriptorPortHamiltonian::new(
        n,
        diff,
        2,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, diff).unwrap()),
    )
    .unwrap();
    PreparedCommonFlowZoh::new_owned(
        descriptor,
        dt,
        LinearSplitBudget {
            max_modes: n,
            max_omega_dt: 1.0,
        },
    )
    .unwrap()
}
fn descriptor_state(
    series: &LaminarSeries<FluidTransportState>,
    state: &LaminarSeriesState,
) -> Vec<f64> {
    let mut x = Vec::new();
    for (i, pipe) in series.sections.iter().enumerate() {
        x.push(pipe.flow_scales[0] * state.volume_flow_m3_s());
        x.extend_from_slice(&state.coordinates[series.offsets[i]..series.offsets[i + 1]]);
    }
    x.push(0.0);
    x
}

#[test]
fn unlike_viscosities_match_full_descriptor_and_continue_nonsteady_material_history() {
    with_cx(|cx| {
        let series =
            LaminarSeries::try_new(cx, vec![pipe(cx, 0.13, 0.001), pipe(cx, 0.27, 0.004)]).unwrap();
        let mut state = series
            .steady_state(VolumetricFlowRate::new(3e-8), 2.0)
            .unwrap();
        let mut heat = [0.0; 2];
        let before = state.kinetic_energy_j();
        let warm = series
            .prepare(cx, 0.017)
            .unwrap()
            .advance(
                cx,
                &mut state,
                [Pressure::new(-80.0), Pressure::new(0.0)],
                &mut heat,
            )
            .unwrap();
        assert_ne!(state.kinetic_energy_j(), before);
        assert!(
            (state.profile_memory_energy_j() - state.mean_kinetic_energy_j() / 3.0).abs()
                > 1e-5 * state.kinetic_energy_j()
        );
        balance(&warm, &heat);
        let retained = state.clone();
        let updated = series
            .with_sources(cx, vec![fluid(0.003), fluid(0.0007)])
            .unwrap();
        assert_eq!(state.coordinates(), retained.coordinates());
        assert_eq!(state.accepted_time_s(), retained.accepted_time_s());
        near(
            (0..2)
                .map(|i| updated.section_kinetic_energy_j(&state, i).unwrap())
                .sum(),
            retained.kinetic_energy_j(),
            retained.kinetic_energy_j(),
            1e-14,
        );
        let reference = descriptor_reference(updated.sections(), 0.021);
        let mut original = descriptor_state(&updated, &state);
        let mut step = updated.prepare(cx, 0.021).unwrap();
        for pressure in [70.0, -40.0, 0.0, 15.0] {
            let expected = reference.step(&original, &[pressure, 0.0]).unwrap();
            let record = step
                .advance(
                    cx,
                    &mut state,
                    [Pressure::new(pressure), Pressure::new(0.0)],
                    &mut heat,
                )
                .unwrap();
            let dim = updated.sections()[0].system().state_dim();
            let flow0 = expected.x[0] / updated.sections()[0].flow_scales[0];
            let flow1 = expected.x[dim] / updated.sections()[1].flow_scales[0];
            near(flow0, flow1, 1e-7, 2e-10);
            near(state.volume_flow_m3_s(), flow0, 1e-7, 3e-9);
            near(record.average_volume_flow_m3_s, expected.y[0], 1e-7, 3e-9);
            near(
                state.kinetic_energy_j(),
                reference.hamiltonian(&expected.x),
                state.kinetic_energy_j(),
                5e-9,
            );
            near(
                record.native.dissipated,
                expected.dissipated,
                expected.dissipated,
                5e-9,
            );
            balance(&record, &heat);
            original = expected.x;
        }
    });
}

#[test]
fn steady_section_heat_uses_each_actual_viscosity_and_harmonic_impedances_add() {
    with_cx(|cx| {
        let series =
            LaminarSeries::try_new(cx, vec![pipe(cx, 0.2, 0.001), pipe(cx, 0.2, 0.004)]).unwrap();
        let flow = 4e-8;
        let mut state = series
            .steady_state(VolumetricFlowRate::new(flow), 0.0)
            .unwrap();
        let initial = state.kinetic_energy_j();
        let mut heat = [0.0; 2];
        let dt = 0.031;
        let report = series
            .prepare(cx, dt)
            .unwrap()
            .advance(
                cx,
                &mut state,
                [
                    Pressure::new(series.steady_resistance_pa_s_m3() * flow),
                    Pressure::new(0.0),
                ],
                &mut heat,
            )
            .unwrap();
        near(state.volume_flow_m3_s(), flow, flow, 2e-10);
        near(state.kinetic_energy_j(), initial, initial, 3e-10);
        for (i, section) in series.sections().iter().enumerate() {
            let exact = section.steady_resistance_pa_s_m3() * flow * flow * dt;
            near(heat[i], exact, exact, 2e-10);
        }
        near(heat[1], 4.0 * heat[0], heat[1], 2e-10);
        balance(&report, &heat);
        for omega in [0.01, 0.7, 25.0, 500.0] {
            let z = series.impedance_at(omega).unwrap();
            let expected = series.sections()[0].impedance_at(omega).unwrap()
                + series.sections()[1].impedance_at(omega).unwrap();
            assert_eq!(z, expected);
            let model = fs_phs::harmonic::HarmonicModel::from_ode(
                series.system(),
                fs_phs::harmonic::HarmonicBudget {
                    max_dimension: 128,
                    max_dense_bytes: 4_000_000,
                    max_dense_work: 4_000_000,
                    relative_residual_tolerance: 1e-11,
                },
                || false,
            )
            .unwrap();
            let response = model
                .factor(omega, || false)
                .unwrap()
                .response(&[C64::ONE, C64::ZERO], || false)
                .unwrap();
            // Native harmonic uses exp(-i*w*t); the wall oracle uses exp(+i*w*t).
            let actual = (C64::ONE / response.output[0]).conj();
            near(actual.re, z.re, z.abs(), 2e-9);
            near(actual.im, z.im, z.abs(), 2e-9);
        }
    });
}

#[test]
fn series_refusals_leave_the_physical_history_and_local_heat_unchanged() {
    with_cx(|cx| {
        let series =
            LaminarSeries::try_new(cx, vec![pipe(cx, 0.2, 0.001), pipe(cx, 0.2, 0.004)]).unwrap();
        let mut state = series
            .steady_state(VolumetricFlowRate::new(2e-8), 4.0)
            .unwrap();
        let before = state.clone();
        let mut heat = [13.0, 17.0];
        let mut step = series.prepare(cx, 0.03).unwrap();
        assert!(
            step.advance(
                cx,
                &mut state,
                [Pressure::new(1e12), Pressure::new(0.0)],
                &mut heat
            )
            .is_err()
        );
        assert_eq!(state.coordinates(), before.coordinates());
        assert_eq!(state.accepted_time_s(), before.accepted_time_s());
        assert_eq!(heat, [13.0, 17.0]);
        assert!(
            step.advance(
                cx,
                &mut state,
                [Pressure::new(10.0), Pressure::new(0.0)],
                &mut heat[..1]
            )
            .is_err()
        );
        assert_eq!(heat, [13.0, 17.0]);
        let foreign = pipe(cx, 0.4, 0.001);
        let mut other = foreign
            .steady_state(VolumetricFlowRate::new(2e-8), 4.0)
            .unwrap();
        assert!(
            step.advance(
                cx,
                &mut other,
                [Pressure::new(10.0), Pressure::new(0.0)],
                &mut heat
            )
            .is_err()
        );
        assert!(series.with_sources(cx, vec![fluid(0.002)]).is_err());
        let different_density = FluidTransportState::try_new(1000.01, 4200.0, 0.6, 0.002).unwrap();
        assert!(
            series
                .with_sources(cx, vec![fluid(0.002), different_density])
                .is_err()
        );
        assert_eq!(state.coordinates(), before.coordinates());
        assert_eq!(state.accepted_time_s(), before.accepted_time_s());
        assert_eq!(heat, [13.0, 17.0]);
    });
}
