//! Independent surface moments, physical work, causal frequency/time parity.

use super::super::super::HarmonicLimits;
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const FACE: [u32; 3] = [0, 1, 2];
const RELEASE: [u32; 3] = [0, 1, 3];
const X: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
const GX: [[f64; 3]; 3] = [
    [1.0 / 10.0, 7.0 / 120.0, 1.0 / 20.0],
    [7.0 / 120.0, 2.0 / 15.0, 7.0 / 120.0],
    [1.0 / 20.0, 7.0 / 120.0, 1.0 / 10.0],
];

fn model<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn field<'m, 'f, 'c>(
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    law: impl Fn([f64; 3]) -> f64,
) -> SpatialResistiveBoundary<'m, 'f, 'c> {
    model
        .spatial_resistive_boundary(
            &[FACE],
            4,
            &mut |_, point| Ok::<_, AcousticFeError>(law(point.point.position_m)),
            &mut || false,
        )
        .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn bilinear(matrix: [[f64; 3]; 3], a: &[f64], b: &[f64]) -> f64 {
    (0..3)
        .flat_map(|i| (0..3).map(move |j| a[i] * matrix[i][j] * b[j]))
        .sum()
}

fn constant_matrix(y: f64) -> [[f64; 3]; 3] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| y * if i == j { 1.0 / 12.0 } else { 1.0 / 24.0 })
    })
}

fn gy() -> [[f64; 3]; 3] {
    let swap = [0, 2, 1];
    core::array::from_fn(|i| core::array::from_fn(|j| 0.4 * GX[swap[i]][swap[j]]))
}

fn entry(matrix: &Csr, i: usize, j: usize) -> f64 {
    let (cols, values) = matrix.row(i);
    cols.binary_search(&j).ok().map_or(0.0, |at| values[at])
}

#[test]
fn original_spatial_real_poles_have_independent_passive_harmonic_forms() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let calls = Cell::new(0);
    let c0 = field(&model, |_| 0.2);
    let cx = field(&model, |x| {
        calls.set(calls.get() + 1);
        1.0 + x[0]
    });
    let cy = field(&model, |x| 0.4 * (1.0 + x[1]));
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&cx, 0.03), (&cy, 0.17)], &mut || false)
        .unwrap();
    assert_eq!(law.faces(), &[FACE]);
    assert_eq!(law.quadrature_order(), 4);
    assert_eq!(law.relaxation_terms(), 2);
    assert_eq!(law.relaxation_time_s(1), Some(0.17));
    for frequency in [0.03, 2.0, 41.0] {
        let harmonic = law.at_frequency(frequency, &mut || false).unwrap();
        let omega = core::f64::consts::TAU * frequency;
        for i in 0..4 {
            for j in 0..4 {
                let mut expected = C64::ZERO;
                if i < 3 && j < 3 {
                    expected = C64::from_re(constant_matrix(0.2)[i][j])
                        + C64::from_re(GX[i][j]) / C64::new(1.0, omega * 0.03)
                        + C64::from_re(gy()[i][j]) / C64::new(1.0, omega * 0.17);
                }
                close(entry(&harmonic.conductance, i, j), expected.re, 4e-15);
                close(entry(&harmonic.susceptance, i, j), expected.im, 4e-15);
            }
        }
    }
    assert_eq!(
        calls.get(),
        16,
        "causal preparation and frequency reuse must not resample"
    );
    assert!(law.at_frequency(0.0, &mut || false).is_err());
    assert!(matches!(
        law.at_frequency(1.0, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
}

#[test]
fn spatial_memories_close_independent_surface_storage_heat_and_exchange() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let c0 = field(&model, |_| 0.2);
    let cx = field(&model, |x| 1.0 + x[0]);
    let cy = field(&model, |x| 0.4 * (1.0 + x[1]));
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&cx, 0.03), (&cy, 0.17)], &mut || false)
        .unwrap();
    let h = 0.001;
    let owner = model
        .prepare_time_with_relaxing_boundaries(
            h,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.memory_dofs(0), Some([0, 1, 2].as_slice()));
    assert_eq!(owner.memory_patch_and_term(1), Some((0, 1)));
    let mut state = owner.state(&[0.0; 4], &[1.0, -0.4, 2.0, 0.3]).unwrap();
    let initial_energy = owner.energy(&state).unwrap().total_j();
    let mut heat = 0.0;
    for _ in 0..80 {
        let result = owner.step(&state, &[], &mut || false).unwrap();
        let p = state
            .pressure_pa()
            .iter()
            .zip(result.state.pressure_pa())
            .map(|(a, b)| 0.5 * (a + b))
            .collect::<Vec<_>>();
        let mut memory_heat = 0.0;
        let mut memory_work = 0.0;
        let mut stored = 0.0;
        for (arm, (matrix, tau)) in [(GX, 0.03), (gy(), 0.17)].into_iter().enumerate() {
            let q0 = state.memory_pressure_pa(arm).unwrap();
            let q1 = result.state.memory_pressure_pa(arm).unwrap();
            let qm = q0
                .iter()
                .zip(q1)
                .map(|(a, b)| 0.5 * (a + b))
                .collect::<Vec<_>>();
            memory_heat += h * bilinear(matrix, &qm, &qm);
            memory_work += h * bilinear(matrix, &p, &qm);
            stored += 0.5 * tau * bilinear(matrix, q1, q1);
        }
        let instantaneous_heat = h * bilinear(constant_matrix(0.2), &p, &p);
        close(
            result.work.instantaneous_dissipation_j,
            instantaneous_heat,
            2e-13 * initial_energy,
        );
        close(
            result.work.relaxation_dissipation_j,
            memory_heat,
            2e-13 * initial_energy,
        );
        close(
            result.work.memory_work_j,
            memory_work,
            2e-13 * initial_energy,
        );
        close(
            owner.energy(&result.state).unwrap().boundary_memory_j,
            stored,
            2e-13 * initial_energy,
        );
        close(
            result.work.acoustic_balance_residual_j(),
            0.0,
            2e-11 * initial_energy,
        );
        close(
            result.work.memory_balance_residual_j(),
            0.0,
            2e-11 * initial_energy,
        );
        heat += memory_heat + instantaneous_heat;
        state = result.state;
    }
    assert!(heat > 0.1 * initial_energy);
    close(
        owner.energy(&state).unwrap().total_j() + heat,
        initial_energy,
        2e-10 * initial_energy,
    );
}

#[test]
fn causal_spatial_wall_matches_the_independent_harmonic_midpoint_relation() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let c0 = field(&model, |_| 0.05);
    let residue = field(&model, |x| 1.0 + x[0]);
    let tau = 0.07;
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, tau)], &mut || false)
        .unwrap();
    let omega = 18.0;
    let h = 0.003;
    // Exact trapezoidal identity: difference/h = i*omega*endpoint-average.
    let theta = 2.0 * (0.5_f64 * omega * h).atan();
    let harmonic_wall = law
        .at_frequency(omega / core::f64::consts::TAU, &mut || false)
        .unwrap();
    let harmonic = model
        .prepare_with_spatial_admittances(
            omega / core::f64::consts::TAU,
            &[],
            &[&harmonic_wall],
            &[],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap()
        .solve(&[(FACE, [C64::ONE; 3])], &mut || false)
        .unwrap();
    let p = &harmonic.pressure_pa;
    let psi = p.iter().map(|z| z.im / omega).collect::<Vec<_>>();
    let q = p[..3]
        .iter()
        .map(|&z| z / C64::new(1.0, omega * tau))
        .collect::<Vec<_>>();
    let owner = model
        .prepare_time_with_relaxing_boundaries(
            h,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut state = owner
        .state_with_memories(
            &psi,
            &p.iter().map(|z| z.re).collect::<Vec<_>>(),
            &[&q.iter().map(|z| z.re).collect::<Vec<_>>()],
        )
        .unwrap();
    for n in 1..=24 {
        let drive = ((n as f64 - 0.5) * theta).cos() * (0.5 * theta).cos();
        let next = owner
            .step(&state, &[(FACE, [drive; 3])], &mut || false)
            .unwrap();
        let phase = C64::new((n as f64 * theta).cos(), (n as f64 * theta).sin());
        for i in 0..4 {
            close(
                next.state.pressure_pa()[i],
                (p[i] * phase).re,
                3e-10 * (1.0 + p[i].abs()),
            );
            close(
                next.state.pressure_potential_pa_s()[i],
                (p[i] * phase).im / omega,
                3e-10 * (1.0 + p[i].abs()) / omega,
            );
        }
        for (i, value) in next.state.memory_pressure_pa(0).unwrap().iter().enumerate() {
            close(*value, (q[i] * phase).re, 3e-10 * (1.0 + q[i].abs()));
        }
        state = next.state;
    }
}

#[test]
fn midpoint_refines_to_independent_three_state_causal_dynamics() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let c0 = field(&model, |_| 0.1);
    let residue = field(&model, |x| 1.0 + x[0]);
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.04)], &mut || false)
        .unwrap();
    // Exact unit-tetra and unit-triangle moments, with only pressure vertex2
    // unconstrained. This ODE oracle does not call FE assembly or its solver.
    let derivative = |s: [f64; 3]| {
        [
            s[1],
            -1080.0 * (s[0] / 12.0 + s[1] / 120.0 + s[2] / 10.0),
            (s[1] - s[2]) / 0.04,
        ]
    };
    let duration = 0.025;
    let mut reference = [0.0, 1.0, 0.0];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(reference);
        let b = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for count in [10, 20, 40] {
        let owner = model
            .prepare_time_with_relaxing_boundaries(
                duration / count as f64,
                &[],
                &[],
                &[&law],
                &[RELEASE],
                RelaxingTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = owner.state(&[0.0; 4], &[0.0, 0.0, 1.0, 0.0]).unwrap();
        for _ in 0..count {
            state = owner.step(&state, &[], &mut || false).unwrap().state;
        }
        let actual = [
            state.pressure_potential_pa_s()[2],
            state.pressure_pa()[2],
            state.memory_pressure_pa(0).unwrap()[2],
        ];
        errors.push(
            actual
                .into_iter()
                .zip(reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(errors[0] > 1e-5, "nonvacuous refinement: {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.27 * pair[0],
            "second-order refinement: {errors:?}"
        );
    }
}

#[test]
fn release_shared_memories_remain_live_and_can_return_energy_to_fluid() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let residue = field(&model, |x| 1.0 + x[0]);
    let law = model
        .spatial_relaxing_boundary(None, &[(&residue, 0.05)], &mut || false)
        .unwrap();
    let h = 0.002;
    let owner = model
        .prepare_time_with_relaxing_boundaries(
            h,
            &[],
            &[],
            &[&law],
            &[RELEASE],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.memory_dofs(0), Some([0, 1, 2].as_slice()));
    let initial = owner
        .state_with_memories(&[0.0; 4], &[0.0; 4], &[&[1.0, 0.0, 0.0]])
        .unwrap();
    let result = owner.step(&initial, &[], &mut || false).unwrap();
    for i in [0, 1, 3] {
        assert_eq!(result.state.pressure_pa()[i], 0.0);
        assert_eq!(result.state.pressure_potential_pa_s()[i], 0.0);
    }
    close(
        result.state.memory_pressure_pa(0).unwrap()[0],
        (0.05 - 0.5 * h) / (0.05 + 0.5 * h),
        2e-15,
    );
    assert!(
        result.state.pressure_pa()[2] < -0.01,
        "shared-coordinate memory must drive free pressure through the original surface form"
    );
    assert!(result.work.memory_work_j < 0.0);
    assert!(result.work.relaxation_dissipation_j > 0.0);
    close(result.work.balance_residual_j(), 0.0, 1e-13);
}

#[test]
fn zero_residue_keeps_geometric_coordinates_without_inverting_the_form() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let zero = field(&model, |_| 0.0);
    let law = model
        .spatial_relaxing_boundary(None, &[(&zero, 0.1)], &mut || false)
        .unwrap();
    let owner = model
        .prepare_time_with_relaxing_boundaries(
            0.01,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.memory_dofs(0), Some([0, 1, 2].as_slice()));
    let state = owner
        .state_with_memories(&[0.0; 4], &[0.0; 4], &[&[1.0, -2.0, 3.0]])
        .unwrap();
    let result = owner.step(&state, &[], &mut || false).unwrap();
    assert_eq!(result.state.pressure_pa(), &[0.0; 4]);
    assert_eq!(owner.energy(&result.state).unwrap().total_j(), 0.0);
    assert_eq!(result.work.relaxation_dissipation_j, 0.0);
    close(
        result.state.memory_pressure_pa(0).unwrap()[0],
        0.095 / 0.105,
        3e-15,
    );
}

#[test]
fn causal_wall_admission_source_work_ownership_and_cancellation_are_preserved() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let model = model(&geometry);
    let coefficient = field(&model, |x| 1.0 + x[0]);
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .spatial_relaxing_boundary(None, &[(&coefficient, tau)], &mut || false)
                .is_err()
        );
    }
    assert!(
        model
            .spatial_relaxing_boundary(None, &[], &mut || false)
            .is_err()
    );
    assert!(
        model
            .spatial_relaxing_boundary(None, &vec![(&coefficient, 0.1); 65], &mut || false)
            .is_err()
    );
    let law = model
        .spatial_relaxing_boundary(None, &[(&coefficient, 0.1)], &mut || false)
        .unwrap();
    assert!(matches!(
        model.prepare_time_with_relaxing_boundaries(
            0.002,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits {
                max_memory_dofs: 2,
                ..RelaxingTimeLimits::default()
            },
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    for (resistive, release) in [
        (
            vec![ResistiveBoundary {
                faces: vec![FACE],
                admittance_m_pa_s: 0.1,
            }],
            vec![],
        ),
        (vec![], vec![FACE]),
    ] {
        assert!(
            model
                .prepare_time_with_relaxing_boundaries(
                    0.002,
                    &resistive,
                    &[],
                    &[&law],
                    &release,
                    RelaxingTimeLimits::default(),
                    &mut || false
                )
                .is_err()
        );
    }
    let owner = model
        .prepare_time_with_relaxing_boundaries(
            0.002,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let other = model
        .prepare_time_with_relaxing_boundaries(
            0.002,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner.state(&[0.0; 4], &[0.0; 4]).unwrap();
    assert!(other.step(&initial, &[], &mut || false).is_err());
    assert!(matches!(
        owner.step(&initial, &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(initial.pressure_pa(), &[0.0; 4]);
    let boundary = model
        .quadratic_velocity_source(&[(FACE, [1.0; 6])], &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let result = owner
        .step_with_boundary_sources(
            &initial,
            &[],
            &[(&boundary, 0.3)],
            &[(&volume, 0.02)],
            &mut || false,
        )
        .unwrap();
    assert!(result.work.boundary_source_work_j > 0.0);
    assert!(result.work.volume_source_work_j > 0.0);
    close(
        result.work.input_work_j,
        result.work.boundary_source_work_j + result.work.volume_source_work_j,
        2e-13,
    );
    close(result.work.balance_residual_j(), 0.0, 2e-13);
}
