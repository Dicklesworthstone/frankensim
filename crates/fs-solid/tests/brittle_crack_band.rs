//! Independent geometry, force and fracture-work checks through native solid
//! consumers. Synthetic coefficients verify the stated fixed-plane model; they
//! do not qualify a real material or establish converged crack topology.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::brittle::{BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings};
use fs_solid::brittle::{TetCrackBandGeometry, TetCrackBandGeometryError};
use fs_solid::linear::Formulation;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const YOUNG: f64 = 30.0e9;
const POISSON: f64 = 0.22;
const STRENGTH: f64 = 3.0e6;
const FRACTURE_WORK: f64 = 100.0;
const TET: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 5, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn close(actual: f64, expected: f64, relative: f64) {
    let scale = expected.abs().max(actual.abs()).max(1.0e-14);
    assert!((actual - expected).abs() <= relative * scale,
        "actual={actual:e}, expected={expected:e}, relative error={:e}",
        (actual-expected).abs()/scale);
}

fn points(length: f64) -> [[f64; 3]; 4] {
    [[0.0, 0.0, 0.0], [length, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.01]]
}

fn band(points: &[[f64; 3]; 4], normal: [f64; 3]) -> TetCrackBandGeometry {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(points, &TET, TetAssemblyBudget::standard(), 4, cx).unwrap();
        TetCrackBandGeometry::from_element(&geometry, 0, normal).unwrap()
    })
}

fn parameters() -> BrittleCrackBandParameters {
    BrittleCrackBandParameters {
        tensile_strength_pa: STRENGTH,
        shear_strength_pa: 4.0e6,
        mode_i_fracture_energy_j_m2: FRACTURE_WORK,
        mode_ii_fracture_energy_j_m2: 200.0,
        bk_exponent: 1.5,
    }
}

fn model(settings: BrittleCrackBandSettings) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(YOUNG, POISSON, 0.03).unwrap(),
        parameters(),
        settings,
    ).unwrap()
}

/// Axial extension with two genuinely free transverse strain DOFs. All shear
/// and rigid motions are constrained; y(node 2) and z(node 3) find equilibrium.
fn axial_prescribed(length: f64) -> Vec<(usize, f64)> {
    (0..12).filter(|&dof| dof != 7 && dof != 11)
        .map(|dof| (dof, if dof == 3 { length } else { 0.0 })).collect()
}

#[test]
fn directional_width_uses_projected_area_and_preserves_frame_and_scale() {
    let p = [[0.0; 3], [0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.125]];
    for axis in 0..3 {
        let mut n = [0.0; 3];
        n[axis] = 1.0;
        let geometry = band(&p, n);
        close(geometry.reference_volume_m3(), 0.5 * 0.25 * 0.125 / 6.0, 2e-14);
        close(geometry.band_width_m(), p[axis + 1][axis] / 3.0, 2e-14);
        close(geometry.reference_volume_m3() / geometry.band_width_m(), geometry.projected_area_m2(), 2e-14);
        assert_eq!(geometry.element(), 0);
        assert_eq!(geometry.nodes(), TET[0]);
    }
    let diagonal = band(&p, [1.0, 1.0, 1.0]);
    let expected_area = (0.5 * 0.25 + 0.25 * 0.125 + 0.5 * 0.125) / (2.0 * 3.0_f64.sqrt());
    close(diagonal.projected_area_m2(), expected_area, 2e-14);
    let transformed = p.map(|x| [-2.0 * x[1] + 1024.0, 2.0 * x[0] - 512.0, 2.0 * x[2] + 256.0]);
    let rotated = band(&transformed, [-1.0, 1.0, 1.0]);
    close(rotated.projected_area_m2(), 4.0 * diagonal.projected_area_m2(), 2e-14);
    close(rotated.reference_volume_m3(), 8.0 * diagonal.reference_volume_m3(), 2e-14);
    close(rotated.band_width_m(), 2.0 * diagonal.band_width_m(), 2e-14);
    close(band(&p, [1.0e308, 0.0, 0.0]).band_width_m(), 0.5 / 3.0, 2e-14);
    close(band(&p, [0.0, 1.0e-320, 0.0]).band_width_m(), 0.25 / 3.0, 2e-14);
    let reverse_normal = band(&p, [-1.0, -1.0, -1.0]);
    close(reverse_normal.projected_area_m2(), diagonal.projected_area_m2(), 2e-14);
}

#[test]
fn geometry_refusals_and_connectivity_winding_preserve_original_admission() {
    with_cx(|cx| {
        let p = points(0.03);
        let geometry = TetLoadGeometry3::new(&p, &TET, TetAssemblyBudget::standard(), 4, cx).unwrap();
        for normal in [[0.0; 3], [f64::NAN, 0.0, 1.0], [0.0, f64::INFINITY, 0.0]] {
            assert_eq!(TetCrackBandGeometry::from_element(&geometry, 0, normal), Err(TetCrackBandGeometryError::InvalidNormal));
        }
        assert_eq!(TetCrackBandGeometry::from_element(&geometry, 1, [1.0, 0.0, 0.0]),
            Err(TetCrackBandGeometryError::ElementOutOfRange { element: 1 }));
        let reverse = [[0, 1, 3, 2]];
        let reversed = TetLoadGeometry3::new(&p, &reverse, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let original_band = TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 2.0, 3.0]).unwrap();
        let reversed_band = TetCrackBandGeometry::from_element(&reversed, 0, [1.0, 2.0, 3.0]).unwrap();
        close(original_band.band_width_m(), reversed_band.band_width_m(), 2e-14);
        assert_eq!(reversed_band.nodes(), reverse[0]);
        let coplanar = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert!(TetLoadGeometry3::new(&coplanar, &TET, TetAssemblyBudget::standard(), 4, cx).is_err());
    });
}

#[test]
fn native_tetrahedron_softens_unloads_fails_and_recovers_compression_without_healing() {
    with_cx(|cx| {
        let p = points(0.03);
        let geometry = band(&p, [1.0, 0.0, 0.0]);
        let material = model(geometry.material_settings(0.03, 1e-7));
        let law = material.small_strain_law();
        let prescribed = axial_prescribed(0.03);
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * geometry.band_width_m());
        let peak = 0.002;
        let peak_stress = STRENGTH * (failure - peak) / (failure - onset);
        let path = [0.5 * onset, peak, 0.25 * peak, peak, 1.2 * failure, -0.001, 0.001];
        let expected = [0.5 * STRENGTH, peak_stress, 0.25 * peak_stress, peak_stress, 0.0, -0.001 * YOUNG, 0.0];
        let mut previous = None;
        let mut peak_damage = 0.0;
        let mut peak_dissipation = 0.0;
        for (index, (&target, &expected_stress)) in path.iter().zip(&expected).enumerate() {
            let load_path = [target];
            let problem = SmallStrainTetProblem {
                nodes_m: &p, tetrahedra: &TET, laws: SmallStrainLawField::Uniform(&law),
                prescribed_m: &prescribed, nodal_forces_n: &[], load_path: &load_path,
                budget: TetAssemblyBudget::standard(), settings: IncrementSettings::default(),
            };
            let solution = if let Some(old) = &previous {
                let old: &fs_solid::SmallStrainTetSolution<_> = old;
                let displacement: Vec<_> = old.displacement_m.iter().flatten().copied().collect();
                problem.solve_from(&displacement, &old.element_states, cx).unwrap()
            } else { problem.solve(cx).unwrap() };
            let sigma = solution.element_stress[0];
            assert!((sigma[0] - expected_stress).abs() <= 1e-7 * STRENGTH,
                "increment {index}: sigma={sigma:?}, expected {expected_stress:e}");
            for value in &sigma[1..] { assert!(value.abs() <= 1e-7 * STRENGTH); }
            close(solution.element_strain[0][0], target, 2e-12);
            assert!((solution.element_strain[0][1] + POISSON * expected_stress / YOUNG).abs() < 2e-10);
            assert!((solution.element_strain[0][2] + POISSON * expected_stress / YOUNG).abs() < 2e-10);
            let loaded_dof = prescribed.iter().position(|(dof, _)| *dof == 3).unwrap();
            let reaction = solution.increments[0].reactions_n[loaded_dof];
            assert!((reaction - expected_stress * geometry.reference_volume_m3() / 0.03).abs() < 1e-5);
            let state = &solution.element_states[0];
            assert!(state.is_valid());
            assert_eq!(*state.strain(), solution.element_strain[0]);
            if index == 0 { assert_eq!(state.damage(), 0.0); }
            if index == 1 {
                peak_damage = state.damage();
                peak_dissipation = state.dissipated_j_m3();
                assert!(peak_damage > 0.0 && peak_damage < 1.0);
                assert!(peak_dissipation > 0.0);
            }
            if index == 2 || index == 3 {
                close(state.damage(), peak_damage, 2e-12);
                close(state.dissipated_j_m3(), peak_dissipation, 2e-12);
            }
            if index >= 4 {
                assert!(state.failed());
                close(state.damage(), 1.0, 2e-12);
                close(state.dissipated_j_m3() * geometry.reference_volume_m3(),
                    FRACTURE_WORK * geometry.projected_area_m2(), 2e-10);
            }
            let displacement: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
            let force = problem.internal_force(&displacement, &solution.element_states, cx).unwrap();
            assert!(force[7].abs() < 1e-5 && force[11].abs() < 1e-5);
            previous = Some(solution);
        }
    });
}

#[test]
fn native_refused_step_keeps_accepted_brittle_history_reusable() {
    with_cx(|cx| {
        let p = points(0.03);
        let material = model(band(&p, [1.0, 0.0, 0.0]).material_settings(0.03, 1e-7));
        let law = material.small_strain_law();
        let prescribed = axial_prescribed(0.03);
        let mut problem = SmallStrainTetProblem {
            nodes_m: &p, tetrahedra: &TET, laws: SmallStrainLawField::Uniform(&law),
            prescribed_m: &prescribed, nodal_forces_n: &[], load_path: &[0.002],
            budget: TetAssemblyBudget::standard(), settings: IncrementSettings::default(),
        };
        let accepted = problem.solve(cx).unwrap();
        let state_before = accepted.element_states.clone();
        let displacement: Vec<_> = accepted.displacement_m.iter().flatten().copied().collect();
        problem.load_path = &[0.08];
        assert!(problem.solve_from(&displacement, &accepted.element_states, cx).is_err());
        assert_eq!(accepted.element_states, state_before);
        problem.load_path = &[0.001];
        let resumed = problem.solve_from(&displacement, &accepted.element_states, cx).unwrap();
        close(resumed.element_states[0].damage(), accepted.element_states[0].damage(), 2e-12);
        close(resumed.element_stress[0][0], 0.5 * accepted.element_stress[0][0], 2e-9);
        close(resumed.element_states[0].dissipated_j_m3(), accepted.element_states[0].dissipated_j_m3(), 2e-12);
    });
}

#[test]
fn different_native_element_widths_preserve_fracture_work_per_projected_area() {
    with_cx(|cx| {
        let first = points(0.03);
        let second = points(0.015).map(|p| [p[0] + 0.1, p[1], p[2]]);
        let nodes: Vec<_> = first.into_iter().chain(second).collect();
        let tetrahedra = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let original = TetLoadGeometry3::new(&nodes, &tetrahedra, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let geometry = [0, 1].map(|e| TetCrackBandGeometry::from_element(&original, e, [1.0, 0.0, 0.0]).unwrap());
        let materials = geometry.map(|g| model(g.material_settings(0.03, 1e-7)));
        let laws = materials.each_ref().map(|m| m.small_strain_law());
        let prescribed: Vec<_> = axial_prescribed(0.03).into_iter()
            .chain(axial_prescribed(0.015).into_iter().map(|(dof, value)| (dof + 12, value))).collect();
        let mut problem = SmallStrainTetProblem {
            nodes_m: &nodes, tetrahedra: &tetrahedra, laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m: &prescribed, nodal_forces_n: &[], load_path: &[0.004],
            budget: TetAssemblyBudget::standard(), settings: IncrementSettings::default(),
        };
        let partial = problem.solve(cx).unwrap();
        for e in 0..2 {
            let failure = 2.0 * FRACTURE_WORK / (STRENGTH * geometry[e].band_width_m());
            let expected = STRENGTH * (failure - 0.004) / (failure - STRENGTH / YOUNG);
            close(partial.element_stress[e][0], expected, 2e-9);
        }
        assert!(partial.element_stress[1][0] > partial.element_stress[0][0]);
        let displacement: Vec<_> = partial.displacement_m.iter().flatten().copied().collect();
        problem.load_path = &[0.016];
        let failed = problem.solve_from(&displacement, &partial.element_states, cx).unwrap();
        for e in 0..2 {
            assert!(failed.element_states[e].failed());
            close(failed.element_states[e].dissipated_j_m3() * geometry[e].reference_volume_m3(),
                FRACTURE_WORK * geometry[e].projected_area_m2(), 2e-10);
        }
        close(geometry[0].projected_area_m2(), geometry[1].projected_area_m2(), 2e-14);
        close(failed.element_states[1].dissipated_j_m3(), 2.0 * failed.element_states[0].dissipated_j_m3(), 2e-10);
    });
}

#[test]
fn native_hex_retains_eight_independent_brittle_histories_and_the_same_surface_work() {
    with_cx(|cx| {
        let length = 0.01;
        let p = [[0.0, 0.0, 0.0], [length, 0.0, 0.0], [length, length, 0.0], [0.0, length, 0.0],
            [0.0, 0.0, length], [length, 0.0, length], [length, length, length], [0.0, length, length]];
        // For this explicitly axis-aligned brick V/A = length. No tetrahedral
        // width is silently reused for a different element or quadrature rule.
        let material = model(BrittleCrackBandSettings {
            band_width_m: length, normal: [1.0, 0.0, 0.0], maximum_strain: 0.03, energy_tolerance_j_m3: 1e-7,
        });
        let law = material.small_strain_law();
        let mut prescribed = Vec::new();
        for (i, p) in p.iter().enumerate() {
            prescribed.push((3 * i, p[0]));
            if p[1] == 0.0 { prescribed.push((3 * i + 1, 0.0)); }
            if p[2] == 0.0 { prescribed.push((3 * i + 2, 0.0)); }
        }
        let problem = SmallStrainHexProblem {
            nodes_m: &p, hexahedra: &[[0, 1, 2, 3, 4, 5, 6, 7]],
            laws: SmallStrainLawField::Uniform(&law), formulation: Formulation::Standard,
            strain_limit: 0.03, rotation_limit_rad: 0.03,
            prescribed_m: &prescribed, nodal_forces_n: &[], load_path: &[0.00005, 0.002, 0.001, 0.008],
            budget: HexAssemblyBudget::default(), settings: IncrementSettings::default(),
        };
        let solution = problem.solve(cx).unwrap();
        assert_eq!(solution.element_offsets, [0, 8]);
        assert_eq!(solution.quadrature_states.len(), 8);
        for state in &solution.quadrature_states { assert!(state.failed()); }
        for stress in &solution.quadrature_stress_pa {
            assert!(stress.iter().all(|value| value.abs() < 1e-7 * STRENGTH));
        }
        let dissipated: f64 = solution.quadrature_states.iter().zip(&solution.quadrature_volume_m3)
            .map(|(state, volume)| state.dissipated_j_m3() * volume).sum();
        close(dissipated, FRACTURE_WORK * length * length, 2e-10);
    });
}
