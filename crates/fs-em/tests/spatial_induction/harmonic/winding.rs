//! Closed 3D current profiles, actual field drives and original copper/end work.
use super::*;
use fs_em::{
    spatial::vector_potential::{
        induction::{
            harmonic::winding::{
                SpatialHarmonicWindingControls, SpatialHarmonicWindingError,
                SpatialHarmonicWindings,
            },
            nonlinear::winding::SpatialWindingConductor,
        },
        winding::{SpatialWindingBasis, SpatialWindingConfig},
    },
    winding::end::{EndWindingInductance, EndWindingMode},
};

pub(super) struct Fixture {
    pub(super) problem: SpatialVectorMagnetostatics,
    pub(super) profiles: Vec<Vec<[f64; 3]>>,
    pub(super) eddy: Vec<Option<ResolvedConductorState>>,
    pub(super) sources: Vec<SpatialWindingConductor>,
    pub(super) velocity: Vec<[f64; 3]>,
}

pub(super) fn winding_controls() -> SpatialHarmonicWindingControls {
    SpatialHarmonicWindingControls {
        field: harmonic_controls(),
        max_windings: 4,
        maximum_response_scalars: 2_000_000,
        relative_terminal_tolerance: 1e-11,
        voltage_tolerance_v: 2e-8,
        power_tolerance_w: 5e-8,
    }
}

pub(super) fn fixture(n: usize, with_eddy: bool) -> Fixture {
    let (problem, all) = cube(n, None);
    let center = problem
        .positions_m()
        .iter()
        .position(|&p| p == [1.0 / n as f64; 3])
        .unwrap() as u32;
    let mut profiles = vec![Vec::new(), Vec::new()];
    let mut occupied = vec![false; problem.mesh().vertex_count];
    let mut eddy = all.clone();
    let mut sources = Vec::new();
    // The same closed-curl construction as the static spatial-winding oracle:
    // J/I = n curl(N_center e_z) carries one ampere around a complete return.
    // Winding two adds an x-axis loop in series, with a nonzero mutual response.
    for (cell, tet) in problem.mesh().tets.iter().enumerate() {
        let gradient = if tet.contains(&center) {
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != center).collect();
            let p = problem.positions_m()[face[0] as usize];
            let sub = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|axis| a[axis] - b[axis]);
            let normal = cross(
                sub(problem.positions_m()[face[1] as usize], p),
                sub(problem.positions_m()[face[2] as usize], p),
            );
            let direction = sub(problem.positions_m()[center as usize], p);
            let height: f64 = normal.into_iter().zip(direction).map(|(a, b)| a * b).sum();
            normal.map(|v| v / height)
        } else {
            [0.0; 3]
        };
        let values = [
            cross(gradient, [0.0, 0.0, n as f64]),
            cross(gradient, [n as f64, 0.0, n as f64]),
        ];
        let mut stranded = false;
        for (winding, value) in values.into_iter().enumerate() {
            profiles[winding].push(value);
            if value != [0.0; 3] {
                stranded = true;
                sources.push(SpatialWindingConductor {
                    winding,
                    cell,
                    fill_fraction: 0.3,
                    conductor: all[cell].as_ref().unwrap().clone(),
                });
            }
        }
        if stranded {
            for &vertex in tet {
                occupied[vertex as usize] = true;
            }
        }
        if stranded || !with_eddy {
            eddy[cell] = None;
        }
    }
    let velocity = problem
        .positions_m()
        .iter()
        .enumerate()
        .map(|(v, p)| {
            if occupied[v] || !with_eddy {
                [0.0; 3]
            } else {
                [0.2 * (p[1] + 0.1), -0.1 * (p[0] - 0.5), 0.05 * p[2]]
            }
        })
        .collect();
    Fixture {
        problem,
        profiles,
        eddy,
        sources,
        velocity,
    }
}

pub(super) fn basis(f: &Fixture, end: EndWindingInductance) -> SpatialWindingBasis<'_> {
    SpatialWindingBasis::new_with_end_modes(
        &f.problem,
        f.profiles.clone(),
        trace(&f.problem, [0.0; 3], false),
        end,
        SpatialWindingConfig {
            max_windings: 4,
            max_basis_values: 1_000_000,
            pairing_tolerance_h: 2e-8,
        },
        static_controls(),
        &mut || false,
    )
    .unwrap()
}

pub(super) fn end_model() -> EndWindingInductance {
    EndWindingInductance::new(
        vec![0.02, 0.03],
        vec![EndWindingMode {
            inductance_h: 0.07,
            turns: vec![1.0, -0.4],
        }],
    )
    .unwrap()
}

pub(super) fn impressed(profiles: &[Vec<[f64; 3]>], currents: &[C64]) -> Vec<[C64; 3]> {
    (0..profiles[0].len())
        .map(|cell| {
            std::array::from_fn(|axis| {
                profiles
                    .iter()
                    .zip(currents)
                    .fold(C64::ZERO, |sum, (p, i)| sum + i.scale(p[cell][axis]))
            })
        })
        .collect()
}

#[test]
fn no_eddy_impedance_matches_original_3d_mutual_inductance_copper_and_voltage_drive() {
    let f = fixture(2, false);
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let mut expected_resistance = [0.0; 2];
    for source in &f.sources {
        let profile = f.profiles[source.winding][source.cell];
        // Every tetrahedron occupies 1/(6*2^3) of the unit cube.
        expected_resistance[source.winding] += source.conductor.resistivity_ohm_m()
            * profile.iter().map(|v| v * v).sum::<f64>()
            / (48.0 * source.fill_fraction);
    }
    for (actual, expected) in prepared.resistance_ohm().iter().zip(expected_resistance) {
        near(*actual, expected, 2e-12);
    }
    let currents = [C64::new(0.7, -0.2), C64::new(-0.3, 0.4)];
    for omega in [0.2, 3.0, 40.0] {
        let frequency = prepared
            .at_frequency(omega, &f.velocity, &mut || false)
            .unwrap();
        assert!(
            frequency
                .unit_reports()
                .iter()
                .all(|r| r.converged_euclidean())
        );
        for row in 0..2 {
            for column in 0..2 {
                complex_near(
                    frequency.impedance_ohm()[2 * row + column],
                    C64::new(
                        if row == column {
                            expected_resistance[row]
                        } else {
                            0.0
                        },
                        -omega * winding_basis.inductance_h()[2 * row + column],
                    ),
                    2e-8,
                );
            }
        }
        let driven = frequency
            .current_response(&currents, &mut || false)
            .unwrap();
        // Cached units are only a seed: native recovery/gates still run, but
        // this linear reference fixture needs no new Krylov recurrence.
        assert_eq!(driven.field.report().iters, 0);
        let expected_voltage: Vec<_> = (0..2)
            .map(|row| {
                (0..2).fold(
                    currents[row].scale(expected_resistance[row]),
                    |sum, column| {
                        sum + C64::new(0.0, -omega * winding_basis.inductance_h()[2 * row + column])
                            * currents[column]
                    },
                )
            })
            .collect();
        for (&actual, &expected) in driven.terminal_voltages_v.iter().zip(&expected_voltage) {
            complex_near(actual, expected, 2e-8);
        }
        for cell in 0..f.problem.mesh().tets.len() {
            for axis in 0..3 {
                let expected = currents[0]
                    .scale(winding_basis.unit_field(0).unwrap().flux_density_t()[cell][axis])
                    + currents[1]
                        .scale(winding_basis.unit_field(1).unwrap().flux_density_t()[cell][axis]);
                complex_near(driven.field.flux_density_t()[cell][axis], expected, 2e-8);
            }
        }
        let mut heat = 0.0;
        for cell in &driven.winding_cells {
            let intrinsic = cell
                .current_density_a_m2
                .iter()
                .map(|j| j.norm_sq())
                .sum::<f64>();
            let power = 0.5
                * cell.conductor_volume_m3
                * cell.source.conductor.resistivity_ohm_m()
                * intrinsic;
            near(cell.mean_power_w, power, 2e-10);
            for axis in 0..3 {
                complex_near(
                    cell.electric_field_v_m[axis],
                    cell.current_density_a_m2[axis]
                        .scale(cell.source.conductor.resistivity_ohm_m()),
                    1e-14,
                );
            }
            heat += power;
        }
        near(driven.mean_winding_power_w, heat, 2e-10);
        near(driven.field.balance().mean_ohmic_power_w, 0.0, 1e-15);
        near(driven.mean_dissipated_power_w(), heat, 2e-10);
        near(driven.supplied_power_va.re, heat, 2e-9);
        near(
            driven.supplied_power_va.im,
            -2.0 * omega * driven.mean_magnetic_energy_j,
            2e-8,
        );
        near(driven.power_residual_w, 0.0, 2e-8);
        near(driven.reactive_power_residual_var, 0.0, 2e-8);
        let voltage = frequency
            .voltage_response(&expected_voltage, &mut || false)
            .unwrap();
        for (&actual, expected) in voltage.currents_a.iter().zip(currents) {
            complex_near(actual, expected, 2e-9);
        }
        assert!(voltage.maximum_terminal_voltage_defect_v < winding_controls().voltage_tolerance_v);
        assert!(
            voltage.terminal_solver_relative_residual.unwrap()
                < winding_controls().relative_terminal_tolerance
        );
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn moving_eddy_windings_keep_native_field_phase_and_separate_mechanical_power() {
    let f = fixture(3, true);
    assert!(f.velocity.iter().any(|v| *v != [0.0; 3]));
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let frequency = prepared
        .at_frequency(3.0, &f.velocity, &mut || false)
        .unwrap();
    let currents = [C64::new(0.7, -0.2), C64::new(-0.3, 0.4)];
    let driven = frequency
        .current_response(&currents, &mut || false)
        .unwrap();
    let direct = state
        .solve_harmonic(
            &impressed(&f.profiles, &currents),
            &f.velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in driven
        .field
        .flux_density_t()
        .iter()
        .flatten()
        .zip(direct.flux_density_t().iter().flatten())
    {
        complex_near(*actual, *expected, 2e-8);
    }
    near(
        driven.field.balance().mean_ohmic_power_w,
        direct.balance().mean_ohmic_power_w,
        2e-8,
    );
    near(
        driven.field.balance().mean_mechanical_power_w,
        direct.balance().mean_mechanical_power_w,
        2e-8,
    );
    assert!(direct.balance().mean_ohmic_power_w > 1e-7);
    assert!(direct.balance().mean_mechanical_power_w.abs() > 1e-12);
    // Stationary winding profiles do not make a moving conducting shield reciprocal.
    assert!((frequency.impedance_ohm()[1] - frequency.impedance_ohm()[2]).abs() > 1e-10);
    near(
        driven.supplied_power_va.re,
        driven.mean_winding_power_w
            + direct.balance().mean_ohmic_power_w
            + direct.balance().mean_mechanical_power_w,
        5e-8,
    );
    near(
        driven.supplied_power_va.im,
        direct.balance().reactive_mechanical_power_var - 6.0 * driven.mean_magnetic_energy_j,
        5e-8,
    );
    let voltage = frequency
        .voltage_response(&driven.terminal_voltages_v, &mut || false)
        .unwrap();
    for (&actual, expected) in voltage.currents_a.iter().zip(currents) {
        complex_near(actual, expected, 2e-8);
    }
    let mut bad_velocity = f.velocity.clone();
    bad_velocity[f.problem.mesh().tets[f.sources[0].cell][0] as usize] = [0.1, 0.0, 0.0];
    assert!(
        prepared
            .at_frequency(3.0, &bad_velocity, &mut || false)
            .is_err()
    );
}

#[test]
fn cancelling_original_end_modes_cannot_erase_active_winding_field_or_copper() {
    let f = fixture(2, false);
    let huge_end = EndWindingInductance::new(
        vec![0.0; 2],
        vec![EndWindingMode {
            inductance_h: 1e22,
            turns: vec![1.0, 1.0],
        }],
    )
    .unwrap();
    let winding_basis = basis(&f, huge_end);
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    let frequency = prepared
        .at_frequency(3.0, &f.velocity, &mut || false)
        .unwrap();
    let currents = [C64::new(1.0, 0.4), C64::new(-1.0, -0.4)];
    let driven = frequency
        .current_response(&currents, &mut || false)
        .unwrap();
    assert_eq!(driven.mean_end_energy_j, 0.0);
    assert!(driven.field.balance().mean_magnetic_energy_j > 1e-5);
    assert!(driven.mean_winding_power_w > 1e-5);
    for row in 0..2 {
        let expected = (0..2).fold(
            currents[row].scale(prepared.resistance_ohm()[row]),
            |sum, column| {
                sum + C64::new(
                    0.0,
                    -3.0 * winding_basis.active_inductance_h()[2 * row + column],
                ) * currents[column]
            },
        );
        complex_near(driven.terminal_voltages_v[row], expected, 2e-8);
    }
    near(
        driven.supplied_power_va.re,
        driven.mean_winding_power_w,
        2e-8,
    );
    near(
        driven.supplied_power_va.im,
        -6.0 * driven.field.balance().mean_magnetic_energy_j,
        2e-8,
    );
}

#[test]
fn harmonic_winding_source_identity_budget_and_cancel_refusals_preserve_original_history() {
    let f = fixture(2, false);
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    let mut overflow = winding_controls();
    overflow.field.maximum_cycles = usize::MAX / 8;
    for control in [
        overflow,
        SpatialHarmonicWindingControls {
            maximum_response_scalars: 0,
            ..winding_controls()
        },
        SpatialHarmonicWindingControls {
            max_windings: 1,
            ..winding_controls()
        },
    ] {
        assert!(
            SpatialHarmonicWindings::new(
                &state,
                &winding_basis,
                f.sources.clone(),
                control,
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        SpatialHarmonicWindings::new(
            &state,
            &winding_basis,
            f.sources[1..].to_vec(),
            winding_controls(),
            &mut || false
        )
        .is_err()
    );
    let mut wrong = f.sources.clone();
    wrong[0].conductor = source(2.0, Some(7.0), None).1.unwrap();
    assert!(
        SpatialHarmonicWindings::new(
            &state,
            &winding_basis,
            wrong,
            winding_controls(),
            &mut || false
        )
        .is_err()
    );
    let mut overlap = f.eddy.clone();
    overlap[f.sources[0].cell] = Some(f.sources[0].conductor.clone());
    let overlapping =
        SpatialInduction::new(&initial, overlap, preparation(), 2e-9, &mut || false).unwrap();
    assert!(
        SpatialHarmonicWindings::new(
            &overlapping,
            &winding_basis,
            f.sources.clone(),
            winding_controls(),
            &mut || false
        )
        .is_err()
    );
    let other = fixture(2, false);
    let other_basis = basis(&other, end_model());
    assert!(
        SpatialHarmonicWindings::new(
            &state,
            &other_basis,
            f.sources.clone(),
            winding_controls(),
            &mut || false
        )
        .is_err()
    );
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        winding_controls(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.at_frequency(3.0, &f.velocity, &mut || true),
        Err(SpatialHarmonicWindingError::Cancelled)
    ));
    let frequency = prepared
        .at_frequency(3.0, &f.velocity, &mut || false)
        .unwrap();
    assert!(
        frequency
            .current_response(&[C64::ONE], &mut || false)
            .is_err()
    );
    assert!(
        frequency
            .voltage_response(&[C64::new(f64::NAN, 0.0); 2], &mut || false)
            .is_err()
    );
    assert!(matches!(
        frequency.current_response(&[C64::ONE; 2], &mut || true),
        Err(SpatialHarmonicWindingError::Cancelled)
    ));
    let zero = frequency
        .current_response(&[C64::ZERO; 2], &mut || false)
        .unwrap();
    near(zero.mean_dissipated_power_w(), 0.0, 1e-15);
    assert_eq!(zero.terminal_voltages_v, vec![C64::ZERO; 2]);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn cached_unit_acceptance_does_not_expand_actual_drive_ampere_allowance() {
    let f = fixture(2, false);
    let winding_basis = basis(&f, end_model());
    let initial = f
        .problem
        .solve(
            &vec![[0.0; 3]; f.profiles[0].len()],
            &trace(&f.problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, f.eddy.clone(), preparation(), 2e-9, &mut || false)
        .unwrap();
    // Deliberately allow the zero unit seed through the algebraic stopping rule.
    // The independent full Ampere gate still has an absolute 100 A allowance.
    // Unit acceptance therefore cannot authorize a million-ampere drive.
    let mut controls = winding_controls();
    controls.field.absolute_residual_a = 1e12;
    controls.field.current_tolerance_a = 100.0;
    let prepared = SpatialHarmonicWindings::new(
        &state,
        &winding_basis,
        f.sources.clone(),
        controls,
        &mut || false,
    )
    .unwrap();
    let frequency = prepared
        .at_frequency(3.0, &f.velocity, &mut || false)
        .unwrap();
    assert!(
        frequency
            .unit_reports()
            .iter()
            .all(|report| report.iters == 0)
    );
    let small = frequency
        .current_response(&[C64::ONE, C64::ZERO], &mut || false)
        .unwrap();
    assert!(small.field.maximum_ampere_residual_a() < controls.field.current_tolerance_a);
    assert!(matches!(
        frequency.current_response(&[C64::new(1e6, 4e5), C64::ZERO], &mut || false),
        Err(SpatialHarmonicWindingError::Field(
            SpatialMagneticError::Invalid(
                "harmonic full Ampere or floating charge equation failed"
            )
        ))
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}
