//! Resolved solid currents against an independent complex slab-diffusion solution.

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use fs_em::{
    planar::{
        eddy::{
            EddyConductor,
            harmonic::{
                HarmonicWindingControls, HarmonicWindingError, PlanarHarmonicWindings,
                SolidWindingRegion,
            },
        },
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        WindingConfig,
        incremental::{IncrementalWindingConfig, PlanarIncrementalWindingBasis},
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_math::c64::C64;
use fs_rep_mesh::{Metric2, TriComplex2};

const HEIGHT: f64 = 0.25;
const DEPTH: f64 = 0.4;
const SIGMA: f64 = 2.0;
const STARTS: [f64; 4] = [-3.0, -1.0, 1.0, 3.0];

fn controls() -> HarmonicWindingControls {
    HarmonicWindingControls {
        relative_tolerance: 1e-10,
        absolute_residual_a_m: 1e-11,
        restart: 80,
        max_cycles: 400,
        max_krylov_scalars: 200_000,
        max_field_scalars: 200_000,
        impedance_tolerance_ohm: 1e-7,
        power_tolerance_w: 1e-8,
        current_balance_tolerance_a: 1e-9,
        max_cell_skin_ratio: 4.0,
    }
}
fn bar(x: f64) -> Option<usize> {
    STARTS
        .iter()
        .position(|&start| x > start && x < start + 1.0)
}
fn mesh(h: f64) -> TriComplex2 {
    let nx = (9.0 / h) as usize;
    let ny = (HEIGHT / h) as usize;
    let vertices: Vec<_> = (0..=nx)
        .flat_map(|x| (0..=ny).map(move |y| [-4.0 + x as f64 * h, y as f64 * h]))
        .collect();
    let mut faces = Vec::new();
    for x in 0..nx {
        for y in 0..ny {
            let a = (x * (ny + 1) + y) as u32;
            let b = a + 1;
            let c = a + (ny + 1) as u32;
            let d = c + 1;
            faces.extend([[a, c, d], [a, d, b]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "four-current-constrained-slabs",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}
fn profiles(problem: &PlanarNonlinearMagnetostatics) -> Vec<Vec<f64>> {
    let mut profiles = vec![vec![0.0; problem.mesh().faces().len()]; 2];
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        let x = face
            .iter()
            .map(|&v| problem.mesh().vertices()[v as usize][0] / 3.0)
            .sum();
        if let Some(section) = bar(x) {
            let winding = if section == 0 || section == 3 { 0 } else { 1 };
            profiles[winding][cell] = if section < 2 {
                1.0 / HEIGHT
            } else {
                -1.0 / HEIGHT
            };
        }
    }
    profiles
}
fn basis<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    profiles: Vec<Vec<f64>>,
) -> PlanarIncrementalWindingBasis<'m> {
    let boundary: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, point)| point[0] == -4.0)
        .map(|(v, _)| (v, 0.0))
        .collect();
    let mut nonlinear = common::material::controls().nonlinear;
    nonlinear.linear.max_iterations = 2000;
    PlanarIncrementalWindingBasis::new(
        problem,
        profiles,
        &boundary,
        vec![0.0; 2],
        vec![0.01, 0.02],
        &vec![0.0; problem.mesh().vertices().len()],
        IncrementalWindingConfig {
            winding: WindingConfig {
                max_windings: 4,
                max_basis_values: 30_000,
                pairing_tolerance_h: 1e-8,
                current_balance_tolerance_turns: 1e-10,
                inverse_residual_tolerance: 1e-9,
            },
            minimum_smooth_field_margin_t: 1e-8,
        },
        nonlinear,
        &mut || false,
    )
    .unwrap()
}
fn sources(
    problem: &PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
) -> (Vec<EddyConductor>, Vec<usize>) {
    let state = common::material::conductor(card, 300.0, ElectricalPropertyBasis::Resistivity);
    let mut conductors = Vec::new();
    let mut insulators = Vec::new();
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        let x = face
            .iter()
            .map(|&v| problem.mesh().vertices()[v as usize][0] / 3.0)
            .sum();
        if let Some(section) = bar(x) {
            conductors.push(EddyConductor {
                triangle: cell,
                region: 10 + section as u32,
                conductor: state.clone(),
            });
        } else {
            insulators.push(cell);
        }
    }
    (conductors, insulators)
}
fn bindings() -> Vec<SolidWindingRegion> {
    vec![
        SolidWindingRegion {
            region: 10,
            winding: 0,
            signed_turns: 1.0,
        },
        SolidWindingRegion {
            region: 11,
            winding: 1,
            signed_turns: 1.0,
        },
        SolidWindingRegion {
            region: 12,
            winding: 1,
            signed_turns: -1.0,
        },
        SolidWindingRegion {
            region: 13,
            winding: 0,
            signed_turns: -1.0,
        },
    ]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn cosh(z: C64) -> C64 {
    C64::new(z.re.cosh() * z.im.cos(), z.re.sinh() * z.im.sin())
}
fn sinh(z: C64) -> C64 {
    C64::new(z.re.sinh() * z.im.cos(), z.re.cosh() * z.im.sin())
}

// In each slab J''=-i omega mu sigma J, mu=1. The magnetic field at an
// interface is the sum of all currents to its left divided by HEIGHT.
// Thus J'(0)=gamma^2 I_before/HEIGHT and J'(1)=gamma^2 I_after/HEIGHT.
// No field matrix, impedance or recovered native field enters this oracle.
fn analytic_current(section: usize, x: f64, omega: f64, currents: [C64; 2]) -> C64 {
    let regional = [currents[0], currents[1], -currents[1], -currents[0]];
    let before = regional[..section]
        .iter()
        .fold(C64::ZERO, |sum, &i| sum + i);
    let after = before + regional[section];
    let gamma = C64::new(0.0, -omega * SIGMA).sqrt();
    let s = gamma.scale(x - STARTS[section]);
    gamma.scale(1.0 / HEIGHT)
        * ((after - before * cosh(gamma)) * cosh(s) / sinh(gamma) + before * sinh(s))
}
fn analytic_power(omega: f64, currents: [C64; 2]) -> f64 {
    let steps = 256;
    let mut integral = 0.0;
    for section in 0..4 {
        for sample in 0..=steps {
            let weight = if sample == 0 || sample == steps {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let j = analytic_current(
                section,
                STARTS[section] + sample as f64 / steps as f64,
                omega,
                currents,
            );
            integral += weight * j.norm_sq() / (3.0 * steps as f64);
        }
    }
    0.5 * HEIGHT * DEPTH / SIGMA * integral
}

#[test]
fn resolved_skin_and_proximity_currents_converge_to_independent_slab_solution() {
    let card = common::source(false, SIGMA);
    let currents = [C64::new(0.018, -0.011), C64::new(0.013, 0.009)];
    let omega = 16.0;
    let expected_power = analytic_power(omega, currents);
    let mut errors = Vec::new();
    for h in [0.25, 0.125] {
        let problem = common::material::problem(&card, mesh(h));
        let basis = basis(&problem, profiles(&problem));
        let (sources, air) = sources(&problem, &card);
        let model = PlanarHarmonicWindings::new_with_solid_regions(
            &basis,
            sources,
            &air,
            vec![],
            bindings(),
            controls(),
            &mut || false,
        )
        .unwrap();
        for &resistance in model.resistance_ohm() {
            near(resistance, 2.0 * DEPTH / (SIGMA * HEIGHT), 1e-12);
        }
        let frequency = model.at_frequency(omega, &mut || false).unwrap();
        let response = frequency.response(&currents, &mut || false).unwrap();
        assert!(response.eddy_cells.is_empty() && response.winding_cells.is_empty());
        assert_eq!(response.solid_region_currents.len(), 4);
        assert_eq!(response.mean_eddy_power_w, 0.0);
        assert_eq!(response.mean_winding_power_w, 0.0);
        near(
            response.mean_dissipated_power_w(),
            response.mean_solid_power_w,
            1e-12,
        );
        near(
            response.mean_supplied_power_w,
            response.mean_solid_power_w,
            1e-8,
        );
        near(response.power_residual_w, 0.0, 1e-8);
        near(response.reactive_power_residual_var, 0.0, 1e-8);
        let mut square_error = 0.0;
        let mut square_reference = 0.0;
        let mut largest_redistribution = 0.0_f64;
        let mut largest_heat_cross_term = 0.0_f64;
        for cell in &response.solid_cells {
            let section = (cell.binding.region - 10) as usize;
            let face = problem.mesh().faces()[cell.source.triangle];
            let dc = currents[cell.binding.winding].scale(cell.binding.signed_turns / HEIGHT);
            for (k, &v) in face.iter().enumerate() {
                let x = problem.mesh().vertices()[v as usize][0];
                let expected = analytic_current(section, x, omega, currents);
                square_error +=
                    cell.volume_m3 * (cell.current_density_a_m2[k] - expected).norm_sq() / 3.0;
                square_reference += cell.volume_m3 * expected.norm_sq() / 3.0;
                largest_redistribution =
                    largest_redistribution.max((cell.current_density_a_m2[k] - dc).abs());
                let e = cell.electric_field_v_m[k];
                assert!((cell.current_density_a_m2[k] - e.scale(SIGMA)).abs() < 1e-11);
            }
            let eddy = cell.current_density_a_m2.map(|j| j - dc);
            let sum = eddy.iter().fold(C64::ZERO, |sum, &j| sum + j);
            let separately_added = cell.volume_m3 / (24.0 * SIGMA)
                * (sum.norm_sq() + eddy.iter().map(|j| j.norm_sq()).sum::<f64>())
                + 0.5 * cell.volume_m3 / SIGMA * dc.norm_sq();
            largest_heat_cross_term =
                largest_heat_cross_term.max((separately_added - cell.mean_power_w).abs());
        }
        let relative_error = (square_error / square_reference).sqrt();
        errors.push(relative_error);
        assert!(largest_redistribution > 0.02); // Actual nonuniform J, not a scalar loss factor.
        assert!(largest_heat_cross_term > 1e-8); // Local cross terms must survive heat recovery.
        for region in &response.solid_region_currents {
            let expected = currents[region.binding.winding].scale(region.binding.signed_turns);
            assert!((region.current_a - expected).abs() < 1e-9);
            assert!(region.residual_a.abs() < 1e-9);
        }
        assert!((response.mean_supplied_power_w - expected_power).abs() / expected_power < 0.15);
    }
    assert!(errors[1] < 0.8 * errors[0], "refinement errors: {errors:?}");
    assert!(errors[1] < 0.08, "fine-grid current error: {}", errors[1]);
}

#[test]
fn neighbor_phase_changes_inner_winding_current_distribution_at_fixed_terminal_current() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh(0.125));
    let basis = basis(&problem, profiles(&problem));
    let (sources, air) = sources(&problem, &card);
    let model = PlanarHarmonicWindings::new_with_solid_regions(
        &basis,
        sources,
        &air,
        vec![],
        bindings(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let inner = C64::new(0.02, 0.005);
    let frequency = model.at_frequency(16.0, &mut || false).unwrap();
    let plus = frequency
        .response(&[C64::new(0.0, 0.015), inner], &mut || false)
        .unwrap();
    let minus = frequency
        .response(&[C64::new(0.0, -0.015), inner], &mut || false)
        .unwrap();
    let mut difference = 0.0_f64;
    for (a, b) in plus.solid_cells.iter().zip(&minus.solid_cells) {
        if a.binding.region == 11 {
            for (&ja, &jb) in a.current_density_a_m2.iter().zip(&b.current_density_a_m2) {
                difference = difference.max((ja - jb).abs());
            }
        }
    }
    assert!(difference > 0.02);
    for response in [&plus, &minus] {
        let region = response
            .solid_region_currents
            .iter()
            .find(|r| r.binding.region == 11)
            .unwrap();
        assert!((region.current_a - inner).abs() < 1e-9);
        near(
            response.mean_dissipated_power_w(),
            response.mean_supplied_power_w,
            1e-8,
        );
    }
    let low = model
        .at_frequency(1e-5, &mut || false)
        .unwrap()
        .response(&[C64::new(0.0, 0.015), inner], &mut || false)
        .unwrap();
    for cell in &low.solid_cells {
        let current = if cell.binding.winding == 0 {
            C64::new(0.0, 0.015)
        } else {
            inner
        };
        let uniform = current.scale(cell.binding.signed_turns / HEIGHT);
        for &j in &cell.current_density_a_m2 {
            assert!((j - uniform).abs() < 1e-5);
        }
    }
}

#[test]
fn mixed_solid_stranded_and_floating_cells_partition_actual_heat_once() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh(0.125));
    let basis = basis(&problem, profiles(&problem));
    let (all, _) = sources(&problem, &card);
    let state = common::material::conductor(&card, 300.0, ElectricalPropertyBasis::Resistivity);
    let mut resolved = Vec::new();
    let mut stranded = Vec::new();
    for source in all {
        if source.region == 11 || source.region == 13 {
            stranded.push(WindingConductor {
                triangle: source.triangle,
                conductor: source.conductor,
            });
        } else {
            resolved.push(source);
        }
    }
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        let x: f64 = face
            .iter()
            .map(|&v| problem.mesh().vertices()[v as usize][0] / 3.0)
            .sum();
        if x > 0.25 && x < 0.75 {
            resolved.push(EddyConductor {
                triangle: cell,
                region: 20,
                conductor: state.clone(),
            });
        }
    }
    let air: Vec<_> = (0..problem.mesh().faces().len())
        .filter(|cell| !resolved.iter().any(|s| s.triangle == *cell))
        .collect();
    let solid = bindings()
        .into_iter()
        .filter(|s| s.region == 10 || s.region == 12)
        .collect();
    let model = PlanarHarmonicWindings::new_with_solid_regions(
        &basis,
        resolved,
        &air,
        stranded,
        solid,
        controls(),
        &mut || false,
    )
    .unwrap();
    let response = model
        .at_frequency(12.0, &mut || false)
        .unwrap()
        .response(&[C64::new(0.02, 0.01), C64::new(-0.01, 0.015)], &mut || {
            false
        })
        .unwrap();
    assert!(
        response.mean_eddy_power_w > 0.0
            && response.mean_winding_power_w > 0.0
            && response.mean_solid_power_w > 0.0
    );
    assert!(response.eddy_cells.iter().all(|s| s.source.region == 20));
    for winding in 0..2 {
        // Each winding deliberately mixes one solid and one stranded section.
        assert!(response.winding_cells.iter().any(|s| s.winding == winding));
        assert!(
            response
                .solid_cells
                .iter()
                .any(|s| s.binding.winding == winding)
        );
    }
    let mut cells = std::collections::BTreeSet::new();
    for cell in response
        .eddy_cells
        .iter()
        .map(|s| s.source.triangle)
        .chain(response.winding_cells.iter().map(|s| s.source.triangle))
        .chain(response.solid_cells.iter().map(|s| s.source.triangle))
    {
        assert!(cells.insert(cell));
    }
    assert!(
        response
            .region_currents_a
            .iter()
            .all(|(_, i)| i.abs() < 1e-9)
    );
    near(
        response.mean_supplied_power_w,
        response.mean_dissipated_power_w(),
        1e-8,
    );
}

#[test]
fn solid_source_admission_refuses_unmatched_bias_profiles_topology_and_double_counting() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh(0.25));
    let basis_ok = basis(&problem, profiles(&problem));
    let (conductors, air) = sources(&problem, &card);
    for bad in [
        vec![],
        vec![bindings()[0], bindings()[0]],
        vec![SolidWindingRegion {
            region: 99,
            winding: 0,
            signed_turns: 1.0,
        }],
        vec![SolidWindingRegion {
            region: 10,
            winding: 2,
            signed_turns: 1.0,
        }],
    ] {
        assert!(
            PlanarHarmonicWindings::new_with_solid_regions(
                &basis_ok,
                conductors.clone(),
                &air,
                vec![],
                bad,
                controls(),
                &mut || false
            )
            .is_err()
        );
    }
    let duplicate = WindingConductor {
        triangle: conductors[0].triangle,
        conductor: conductors[0].conductor.clone(),
    };
    assert!(
        PlanarHarmonicWindings::new_with_solid_regions(
            &basis_ok,
            conductors.clone(),
            &air,
            vec![duplicate],
            bindings(),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut changed = profiles(&problem);
    let selected: Vec<_> = changed[0]
        .iter()
        .enumerate()
        .filter(|(_, v)| **v > 0.0)
        .map(|(cell, _)| cell)
        .take(2)
        .collect();
    changed[0][selected[0]] += 0.1;
    changed[0][selected[1]] -= 0.1; // Still current-balanced; wrong DC distribution.
    let wrong_basis = basis(&problem, changed);
    assert!(
        PlanarHarmonicWindings::new_with_solid_regions(
            &wrong_basis,
            conductors.clone(),
            &air,
            vec![],
            bindings(),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut reverse = bindings();
    reverse[0].signed_turns = -1.0;
    assert!(
        PlanarHarmonicWindings::new_with_solid_regions(
            &basis_ok,
            conductors.clone(),
            &air,
            vec![],
            reverse,
            controls(),
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PlanarHarmonicWindings::new_with_solid_regions(
            &basis_ok,
            conductors,
            &air,
            vec![],
            bindings(),
            controls(),
            &mut || true
        ),
        Err(HarmonicWindingError::Cancelled)
    ));
    assert!(
        basis_ok
            .operating_field()
            .potential_wb_m()
            .iter()
            .all(|&a| a == 0.0)
    );
}
