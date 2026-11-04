use super::super::super::frequency::tests::{loaded, oscillator};
use super::*;
use fs_la::eigen_complex::lu_complex_with_cancel;

const WORK: usize = 100_000_000;
const COLUMN: [f64; 1] = [0.003];
const FORCE: [C64; 1] = [C64::new(1e-5, -4e-6)];

pub(super) fn preparation(reference: f64) -> DgPolePreparation {
    DgPolePreparation {
        reference_rate_per_s: reference,
        maximum_dimension: 256,
        maximum_payload_bytes: 16 << 20,
        maximum_work: 100_000_000,
    }
}
pub(super) fn pole_options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 1800.0,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-11,
        max_iterations: 32,
        max_dimension: 256,
        max_dense_bytes: 16 << 20,
        max_work: 1_000_000_000,
    }
}
pub(super) fn options() -> DgPoleResponseOptions {
    DgPoleResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: 1800.0,
            pole_residual_tolerance: 1e-11,
            pole_step_tolerance: 1e-10,
            linear_residual_tolerance: 1e-11,
            max_dimension: 256,
            max_vectors: 8,
            max_dense_bytes: 16 << 20,
            max_work: 1_000_000_000,
        },
        max_observations: 4,
        max_physical_bytes: 16 << 20,
        max_physical_work: 100_000_000,
    }
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-8 * actual.abs().max(expected.abs()) + 1e-28,
        "{actual:?} != {expected:?}"
    );
}
pub(super) fn vector_error(actual: &[C64], expected: &[C64]) -> f64 {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 { error } else { error / scale }
}
fn resolve(prepared: &PreparedModalDgPoles<'_>, s: C64, force: &[C64]) -> Vec<C64> {
    let evaluation = prepared.pencil.evaluate(s, &mut || false).unwrap();
    let mut state = vec![C64::ZERO; prepared.dimension()];
    state[prepared.coordinates.modal_displacement.clone()].copy_from_slice(force);
    prepared
        .pencil
        .scale_rhs(&mut state, &mut || false)
        .unwrap();
    lu_complex_with_cancel(&evaluation.matrix, prepared.dimension(), || false)
        .unwrap()
        .solve_with_cancel(&mut state, || false)
        .unwrap();
    state
}
pub(super) fn physical_pole(
    prepared: &PreparedModalDgPoles<'_>,
    gate: &CancelGate,
) -> DgLoadedModalPole {
    let response = prepared
        .source
        .harmonic_response_with_load(1810.0, &FORCE, prepared.load, &COLUMN, WORK, gate)
        .unwrap();
    let seed = prepared.seed_from_harmonic_response(&response).unwrap();
    prepared
        .resonance_pole(C64::new(-10.0, -1810.0), &seed, &pole_options(), gate)
        .unwrap()
}

#[test]
fn full_pencil_matches_native_ac_and_residues_match_its_complex_resolvent_limit() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    for circular in [false, true] {
        let load = loaded(circular, true);
        let saved = load.state().clone();
        let branch_saved = load.branches()[0].state().clone();
        let prepared = system
            .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
            .unwrap();
        let layout = prepared.coordinates();
        assert!(layout.branch_energy[0].1.len() > 4);
        assert_eq!(layout.terminal_flow.len(), 2);
        assert_eq!(layout.channels[0].velocity_cells.len(), 8);
        assert_eq!(layout.material_memory.len(), 1);
        // The polynomial composer must agree with the independent source-owned
        // harmonic route, including every condensed-and-restored memory state.
        let ac = system
            .harmonic_response_with_load(1300.0, &FORCE, &load, &COLUMN, WORK, &gate)
            .unwrap();
        let packed = prepared.seed_from_harmonic_response(&ac).unwrap();
        let resolved = resolve(&prepared, C64::new(0.0, -1300.0), &FORCE);
        let mut groups = vec![
            layout.modal_displacement.clone(),
            layout.radiation_displacement.clone(),
            layout.duct_energy.clone(),
            layout.terminal_state.clone(),
            layout.material_memory.clone(),
        ];
        groups.extend(layout.branch_energy.iter().map(|(_, span)| span.clone()));
        groups.extend(
            layout
                .channels
                .iter()
                .map(|channel| channel.velocity_cells.clone()),
        );
        for group in groups {
            assert!(vector_error(&resolved[group.clone()], &packed[group]) < 3e-8);
        }
        let indices: Vec<_> = layout
            .terminal_flow
            .iter()
            .chain(&layout.junction_flow)
            .map(|(_, index)| *index)
            .collect();
        assert!(
            vector_error(
                &indices.iter().map(|&i| resolved[i]).collect::<Vec<_>>(),
                &indices.iter().map(|&i| packed[i]).collect::<Vec<_>>()
            ) < 3e-8
        );
        for index in layout
            .channels
            .iter()
            .filter_map(|channel| channel.pressure_drop)
        {
            close(resolved[index], packed[index]);
        }
        let pole = physical_pole(&prepared, &gate);
        let response = prepared
            .prepare_pole_response(&pole, options(), &gate)
            .unwrap();
        let mut pressure = vec![C64::ZERO; prepared.dimension()];
        pressure[layout.duct_energy.clone()]
            .iter_mut()
            .zip(&prepared.pressure_column)
            .for_each(|(value, &weight)| *value = C64::from_re(weight));
        let mut imaginary = vec![C64::ZERO; prepared.dimension()];
        imaginary[layout.material_memory.start] = C64::new(0.3, 0.7);
        let result = response
            .forcing_residues(&[&FORCE], &[&pressure, &imaginary], &gate)
            .unwrap()
            .remove(0);
        close(result.observations[0], result.inlet_pressure_pa);
        close(
            result.observations[1],
            C64::new(0.3, 0.7) * result.state[layout.material_memory.start],
        );
        assert!(
            (result.observations[1]
                - C64::new(0.3, -0.7) * result.state[layout.material_memory.start])
                .abs()
                > 0.5 * result.observations[1].abs()
        );
        close(
            result.modal_velocity_m_sqrt_kg_per_s[0],
            pole.laplace_pole_per_s * result.state[0],
        );
        close(
            result.inlet_flow_m3_s,
            result.modal_velocity_m_sqrt_kg_per_s[0].scale(COLUMN[0]),
        );
        assert_eq!(result.state.len(), prepared.dimension());
        assert!(result.linear_residual < 1e-11 && result.original_relative_residual < 1e-11);
        let mut errors = Vec::new();
        for delta in [0.1, 0.03, 0.01] {
            let near = resolve(
                &prepared,
                pole.laplace_pole_per_s + C64::from_re(delta),
                &FORCE,
            );
            let laurent: Vec<_> = near.iter().map(|value| value.scale(delta)).collect();
            errors.push(vector_error(&laurent, &result.state));
        }
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-4, "{errors:?}");
        assert_eq!(load.state(), &saved);
        assert_eq!(load.branches()[0].state(), &branch_saved);
    }
    assert_eq!(system.material_memory(), &[0.004]);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn physical_residues_preserve_row_scaling_gauge_phase_and_extreme_force_homogeneity() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let load = loaded(false, false);
    let prepared = system
        .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
        .unwrap();
    let pole = physical_pole(&prepared, &gate);
    let response = prepared
        .prepare_pole_response(&pole, options(), &gate)
        .unwrap();
    let rotated = [C64::new(0.0, 1.0) * FORCE[0]];
    let tiny = [FORCE[0].scale(1e-180)];
    let large = [FORCE[0].scale(1e180)];
    let zero = [C64::ZERO];
    let results = response
        .forcing_residues(&[&FORCE, &rotated, &tiny, &large, &zero], &[], &gate)
        .unwrap();
    for (result, scale) in results.iter().skip(1).zip([
        C64::new(0.0, 1.0),
        C64::from_re(1e-180),
        C64::from_re(1e180),
        C64::ZERO,
    ]) {
        // Remove the declared input scale first so a fixed absolute tolerance
        // cannot turn the small-force test into a vacuous zero comparison.
        if scale != C64::ZERO {
            let normalized: Vec<_> = result.state.iter().map(|&value| value / scale).collect();
            assert!(vector_error(&normalized, &results[0].state) < 3e-8);
            close(
                result.inlet_pressure_pa / scale,
                results[0].inlet_pressure_pa,
            );
            close(result.participation / scale, results[0].participation);
        } else {
            assert!(result.state.iter().all(|&value| value == C64::ZERO));
            assert_eq!(result.original_relative_residual, 0.0);
        }
    }
    let mut regauged = pole.clone();
    for value in &mut regauged.state {
        *value = *value * C64::new(-3.0, 2.0);
    }
    regauged.inlet_pressure_pa = C64::new(f64::NAN, f64::NAN);
    regauged.relative_residual = 0.0;
    let other = system
        .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(36_000.0), &gate)
        .unwrap();
    let cache = other
        .prepare_pole_response(&regauged, options(), &gate)
        .unwrap();
    let alternate = cache
        .forcing_residues(&[&FORCE], &[], &gate)
        .unwrap()
        .remove(0);
    assert!(vector_error(&alternate.state, &results[0].state) < 3e-8);
    close(alternate.inlet_pressure_pa, results[0].inlet_pressure_pa);
    let replay = response
        .forcing_residues(&[&FORCE], &[], &gate)
        .unwrap()
        .remove(0);
    assert_eq!(replay.state, results[0].state);
    assert_eq!(replay.participation, results[0].participation);
}

fn dark_oscillator() -> CoupledModalSystem {
    oscillator()
        .with_additional_material_relaxation(
            vec![fs_phs::RelaxationBranch {
                projection: vec![0.0, 0.0],
                stiffness: 3.0,
                relaxation_time_s: 0.002,
            }],
            Some(vec![0.003]),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

#[test]
fn a_dark_material_mode_remains_present_with_zero_forced_residue_and_native_refusals() {
    let gate = CancelGate::new_clock_free();
    let system = dark_oscillator();
    let load = loaded(false, false);
    let prepared = system
        .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
        .unwrap();
    let memory = prepared.coordinates().material_memory.end - 1;
    let mut seed = vec![C64::ZERO; prepared.dimension()];
    seed[memory] = C64::ONE;
    let pole = prepared
        .resonance_pole(C64::from_re(-500.0), &seed, &pole_options(), &gate)
        .unwrap();
    let response = prepared
        .prepare_pole_response(&pole, options(), &gate)
        .unwrap();
    assert_eq!(response.right_mode().state[memory], C64::ONE);
    let force = [C64::new(1.0, 2.0)];
    let residue = response
        .forcing_residues(&[&force], &[&seed], &gate)
        .unwrap()
        .remove(0);
    assert_eq!(residue.participation, C64::ZERO);
    assert!(residue.state.iter().all(|&value| value == C64::ZERO));
    assert_eq!(residue.observations, &[C64::ZERO]);
    assert_eq!(residue.original_relative_residual, 0.0);
    assert_eq!(response.coordinates().material_memory.len(), 2);
    let mut invalid_pole = pole.clone();
    invalid_pole.state[0] = C64::ONE;
    invalid_pole.relative_residual = 0.0;
    assert!(
        prepared
            .prepare_pole_response(&invalid_pole, options(), &gate)
            .is_err()
    );
    invalid_pole = pole.clone();
    invalid_pole.laplace_pole_per_s = C64::new(-500.0, -1e6);
    assert!(matches!(
        prepared.prepare_pole_response(&invalid_pole, options(), &gate),
        Err(DgPoleResponseError::Physical(_))
    ));
    let mut limits = options();
    limits.max_physical_bytes = 1;
    assert!(
        prepared
            .prepare_pole_response(&pole, limits, &gate)
            .is_err()
    );
    limits = options();
    limits.numerical.max_work = 1;
    assert!(
        prepared
            .prepare_pole_response(&pole, limits, &gate)
            .is_err()
    );
    assert!(response.forcing_residues(&[], &[], &gate).is_err());
    assert!(
        response
            .forcing_residues(&[&[C64::ONE, C64::ONE]], &[], &gate)
            .is_err()
    );
    assert!(
        response
            .forcing_residues(&[&force], &[&[C64::ONE]], &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        response.forcing_residues(&[&force], &[], &cancelled),
        Err(DgPoleResponseError::Cancelled)
    ));
    assert_eq!(
        response
            .forcing_residues(&[&force], &[&seed], &gate)
            .unwrap()[0]
            .state,
        residue.state
    );
    assert_eq!(system.material_memory(), &[0.004, 0.003]);
    assert_eq!(system.samples_rendered(), 0);
    assert_eq!(load.state().step_count(), 0);
}
