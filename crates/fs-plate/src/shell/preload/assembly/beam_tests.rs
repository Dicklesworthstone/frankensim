use super::tests::{Part, budget, options, seam};
use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::preload::{equilibrate_shell_materials, tests::reference_newton};
use crate::shell::stiffened::BeamSection;

fn reinforcement(first: usize, segments: usize) -> Vec<ShellBeam> {
    (0..segments)
        .map(|i| ShellBeam {
            nodes: [first + i, first + i + 1],
            offsets_m: [[0., 0.002, 0.006]; 2],
            section_z: [0., 0., 1.],
            section: BeamSection {
                young_pa: 11e9,
                shear_pa: 0.7e9,
                density_kg_m3: 450.,
                area_m2: 8e-5,
                iy_m4: 4.266666666666667e-10,
                iz_m4: 6.666666666666667e-10,
                torsion_m4: 1e-9,
            },
        })
        .collect()
}

#[test]
fn tied_eccentric_beam_stretch_matches_the_monolithic_loaded_shell() {
    let mut whole = Part::rectangle(0, 4, 2, 0., true);
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 2, 0., false);
    let whole_beams = reinforcement(5, 4);
    let left_beams = reinforcement(3, 2);
    let right_beams = reinforcement(3, 2);
    // The slave-seam force must reach both shell and beam equilibrium through
    // the original reciprocal map. Axial and transverse loads activate stretch.
    right.loads[6 * 3 + 2] = -0.25;
    right.loads[6 * 5 + 2] = -0.5;
    right.loads[6 * 5] = 10.;
    whole.loads[6 * 7 + 2] = -0.25;
    whole.loads[6 * 9 + 2] = -0.5;
    whole.loads[6 * 9] = 10.;
    let direct = equilibrate_shell_materials(
        &whole.mesh,
        &whole.sections,
        &whole.fixed,
        ShellSupport::Clamped,
        &whole_beams,
        &whole.loads,
        &[],
        options(),
        &mut reference_newton,
    )
    .unwrap();
    let mut components = [left.input(), right.input()];
    components[0].beams = &left_beams;
    components[1].beams = &right_beams;
    let mut nonlinear_force = 0.0_f64;
    let tied = equilibrate_shell_assembly(
        &components,
        &[seam(2, TieOffset::Midsurface)],
        &[],
        budget(),
        options(),
        &mut |problem, x, limit| {
            assert_eq!(problem.beams.len(), 4);
            let (x, iterations) = reference_newton(problem, x, limit)?;
            let u = problem.physical(&x);
            let mut force = vec![0.; problem.dimension()];
            for beam in problem.beams {
                beam.add_force_energy(&u, &mut force, problem.options.maximum_gradient)
                    .expect("accepted equilibrium has admitted beam stretch");
            }
            nonlinear_force =
                nonlinear_force.max(force.iter().copied().map(f64::abs).fold(0., f64::max));
            Ok((x, iterations))
        },
        || false,
    )
    .unwrap();
    assert!(
        nonlinear_force > 1e-5,
        "beam correction is inactive: {nonlinear_force}"
    );
    assert!(tied.displacements[1][6 * 5 + 2] < -1e-5);
    assert!(tied.residual_force_n < 2e-8);
    let n = direct.model.free;
    assert_eq!(tied.assembly.model.free, n);
    let mut permutation = vec![usize::MAX; n];
    for component in 0..2 {
        for local in 0..9 {
            let node = 5 * (local / 3) + local % 3 + 2 * component;
            for dof in 0..6 {
                let actual = tied.displacements[component][6 * local + dof];
                let expected = direct.displacement[6 * node + dof];
                assert!(
                    (actual - expected).abs() < 2e-10,
                    "component={component} node={local} dof={dof}: {actual} vs {expected}"
                );
                if let Some(retained) = tied.assembly.map[component][local][dof] {
                    permutation[retained] = direct.model.dof_map[6 * node + dof].unwrap();
                }
            }
        }
    }
    assert!(permutation.iter().all(|&i| i < n));
    for i in 0..n {
        for j in 0..n {
            let actual = tied.assembly.model.k.get(i, j);
            let expected = direct.model.k.get(permutation[i], permutation[j]);
            assert!(
                (actual - expected).abs() < 2e-7 + 1e-9 * expected.abs(),
                "loaded tangent ({i},{j}): {actual} vs {expected}"
            );
            let actual = tied.assembly.model.m.get(i, j);
            let expected = direct.model.m.get(permutation[i], permutation[j]);
            assert!(
                (actual - expected).abs() < 1e-13,
                "reference inertia ({i},{j}): {actual} vs {expected}"
            );
        }
    }
    assert!((tied.stored_energy_j - direct.stored_energy_j).abs() < 1e-11);
    let unloaded = assemble_stiffened_shell_materials(
        &whole.mesh,
        &whole.sections,
        &whole.fixed,
        ShellSupport::Clamped,
        &whole_beams,
    )
    .unwrap();
    let z = direct.model.dof_map[6 * 9 + 2].unwrap();
    assert!((direct.model.k.get(z, z) - unloaded.k.get(z, z)).abs() > 1e-3);
}

#[test]
fn one_beam_stiffened_part_preserves_original_equilibrium_and_pencil_bits() {
    let mut part = Part::rectangle(0, 4, 2, 0., true);
    let beams = reinforcement(5, 4);
    part.loads[6 * 9] = 10.;
    part.loads[6 * 9 + 2] = -0.5;
    let direct = equilibrate_shell_materials(
        &part.mesh,
        &part.sections,
        &part.fixed,
        ShellSupport::Clamped,
        &beams,
        &part.loads,
        &[],
        options(),
        &mut reference_newton,
    )
    .unwrap();
    let mut input = part.input();
    input.beams = &beams;
    let assembled = equilibrate_shell_assembly(
        &[input],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(assembled.displacements[0], direct.displacement);
    assert_eq!(assembled.residual_force_n, direct.residual_force_n);
    assert_eq!(assembled.stored_energy_j, direct.stored_energy_j);
    assert_eq!(assembled.iterations, direct.iterations);
    for row in 0..direct.model.free {
        assert_eq!(assembled.assembly.model.k.row(row), direct.model.k.row(row));
        assert_eq!(assembled.assembly.model.m.row(row), direct.model.m.row(row));
    }
}
