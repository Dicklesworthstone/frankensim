use super::*;
use crate::PlateSection;
use crate::shell::ShellSupport;
use crate::shell::preload::{PreloadOptions, equilibrate_stiffened_shell, tests::reference_newton};
use crate::shell::stiffened::{BeamSection, assemble_stiffened_shell};

fn section() -> BeamSection {
    BeamSection {
        young_pa: 4e7,
        shear_pa: 1.5e7,
        density_kg_m3: 600.,
        area_m2: 4e-4,
        iy_m4: 1e-8,
        iz_m4: 2e-8,
        torsion_m4: 1.5e-8,
    }
}

fn straight() -> (ShellMesh, ShellModel, ShellBeam, BeamStretch) {
    let mesh = ShellMesh::new(
        vec![[0., 0., 0.], [1.5, 0., 0.], [0., 1., 0.]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let beam = ShellBeam {
        nodes: [0, 1],
        offsets_m: [[0.; 3]; 2],
        section_z: [0., 0., 1.],
        section: section(),
    };
    let plate = PlateSection::isotropic(2e6, 0.25, 0.004, 400.).unwrap();
    let model =
        assemble_stiffened_shell(&mesh, &[plate], &[], ShellSupport::Free, &[beam]).unwrap();
    let stretch = BeamStretch::new(&mesh, &model, &beam).unwrap();
    (mesh, model, beam, stretch)
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1. + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn signed_extension_has_exact_axial_energy_force_and_geometric_work() {
    let (_, model, beam, stretch) = straight();
    let l = 1.5;
    let ea = beam.section.young_pa * beam.section.area_m2;
    for extension in [-0.02, 0., 0.02] {
        let mut u = vec![0.; model.free];
        u[model.dof_map[6].unwrap()] = extension * l;
        let strain = extension + 0.5 * extension * extension;
        let mut force = vec![0.; model.free];
        let energy = stretch.add_force_energy(&u, &mut force, 0.1).unwrap();
        near(
            energy + 0.5 * ea * l * extension * extension,
            0.5 * ea * l * strain * strain,
            1e-13,
        );
        near(
            force[6] + ea * extension,
            ea * strain * (1. + extension),
            1e-13,
        );
        near(force[0], -force[6], 1e-13);
        // An independent closed-form integral of the cubic Hermite slopes.
        let v = [0.03, -0.02, -0.01, 0.04];
        let mut probe = vec![0.; model.free];
        for (slot, value) in [1, 5, 7, 11].into_iter().zip(v) {
            probe[slot] = value;
        }
        let g = [
            [36., 3. * l, -36., 3. * l],
            [3. * l, 4. * l * l, -3. * l, -l * l],
            [-36., -3. * l, 36., -3. * l],
            [3. * l, -l * l, -3. * l, 4. * l * l],
        ];
        let work = ea * strain / (30. * l)
            * (0..4)
                .map(|i| (0..4).map(|j| v[i] * g[i][j] * v[j]).sum::<f64>())
                .sum::<f64>();
        let mut action = vec![0.; model.free];
        stretch.add_action(&u, &probe, &mut action, 0.1).unwrap();
        near(
            probe.iter().zip(action).map(|(v, a)| v * a).sum(),
            work,
            1e-13,
        );
    }
}

#[test]
fn eccentric_spatial_energy_force_tangent_and_matrix_free_action_agree() {
    let (mut mesh, _, mut beam, _) = straight();
    mesh.nodes[0] = [0.1, -0.2, 0.3];
    mesh.nodes[1] = [1.3, 0.2, 0.8];
    beam.offsets_m = [[0.01, -0.03, 0.04], [-0.02, 0.01, 0.05]];
    let plate = PlateSection::isotropic(2e6, 0.25, 0.004, 400.).unwrap();
    let model =
        assemble_stiffened_shell(&mesh, &[plate], &[], ShellSupport::Free, &[beam]).unwrap();
    let stretch = BeamStretch::new(&mesh, &model, &beam).unwrap();
    let u: Vec<_> = (0..model.free)
        .map(|i| ((i as f64 + 0.3) * 1.7).sin() * 0.002)
        .collect();
    let mut force = vec![0.; model.free];
    stretch.add_force_energy(&u, &mut force, 0.1).unwrap();
    let mut coo = Coo::new(model.free, model.free);
    stretch.add_tangent(&u, &mut coo, 0.1).unwrap();
    let k = coo.assemble();
    let h = 1e-7;
    for i in 0..12 {
        let mut plus = u.clone();
        let mut minus = u.clone();
        plus[i] += h;
        minus[i] -= h;
        let mut fp = vec![0.; model.free];
        let mut fm = fp.clone();
        let ep = stretch.add_force_energy(&plus, &mut fp, 0.1).unwrap();
        let em = stretch.add_force_energy(&minus, &mut fm, 0.1).unwrap();
        near((ep - em) / (2. * h), force[i], 2e-8);
        let mut v = vec![0.; model.free];
        v[i] = 1.;
        let mut action = v.clone();
        action.fill(0.);
        stretch.add_action(&u, &v, &mut action, 0.1).unwrap();
        for j in 0..12 {
            near((fp[j] - fm[j]) / (2. * h), action[j], 2e-8);
            near(k.get(j, i), action[j], 1e-12);
            assert_eq!(k.get(j, i).to_bits(), k.get(i, j).to_bits());
        }
    }
}

#[test]
fn overlapping_sparse_tie_rows_preserve_nonlinear_virtual_work() {
    let (_, model, _, mut stretch) = straight();
    let n = 7;
    let rows: Vec<Vec<(usize, f64)>> = (0..model.free)
        .map(|i| vec![(i % n, 0.7), ((i + 3) % n, -0.2)])
        .collect();
    let u: Vec<_> = (0..n).map(|i| 0.0004 * (i as f64 - 2.)).collect();
    let v: Vec<_> = (0..n).map(|i| 0.3 * (i as f64 + 1.)).collect();
    let expand = |x: &[f64]| {
        rows.iter()
            .map(|r| r.iter().map(|&(d, w)| w * x[d]).sum())
            .collect::<Vec<_>>()
    };
    let local = expand(&u);
    let direction = expand(&v);
    let mut full_force = vec![0.; model.free];
    let energy = stretch
        .add_force_energy(&local, &mut full_force, 0.1)
        .unwrap();
    let mut full_action = vec![0.; model.free];
    stretch
        .add_action(&local, &direction, &mut full_action, 0.1)
        .unwrap();
    stretch.transform = Some(std::sync::Arc::new(rows.clone()));
    let mut force = vec![0.; n];
    near(
        stretch.add_force_energy(&u, &mut force, 0.1).unwrap(),
        energy,
        1e-14,
    );
    let mut coo = Coo::new(n, n);
    stretch.add_tangent(&u, &mut coo, 0.1).unwrap();
    let k = coo.assemble();
    let mut action = vec![0.; n];
    stretch.add_action(&u, &v, &mut action, 0.1).unwrap();
    let mut expected_force = vec![0.; n];
    let mut expected_action = vec![0.; n];
    for (i, row) in rows.iter().enumerate() {
        for &(d, w) in row {
            expected_force[d] += w * full_force[i];
            expected_action[d] += w * full_action[i];
        }
    }
    let mut applied = vec![0.; n];
    k.spmv(&v, &mut applied);
    for i in 0..n {
        near(force[i], expected_force[i], 1e-13);
        near(action[i], expected_action[i], 1e-13);
        near(applied[i], action[i], 1e-12);
    }
}

#[test]
fn gradient_admission_covers_segment_endpoints_outside_gauss_samples() {
    let (_, model, _, stretch) = straight();
    let mut u = vec![0.; model.free];
    u[5] = 0.105;
    let values = stretch.values(&u);
    assert!(stretch.samples.iter().all(|s| {
        let g = apply(&s.rows, &values);
        g[0].hypot(g[1]).hypot(g[2]) < 0.1
    }));
    assert!(
        stretch
            .add_force_energy(&u, &mut vec![0.; model.free], 0.1)
            .is_none()
    );
    assert!(
        stretch
            .add_action(&u, &u, &mut vec![0.; model.free], 0.1)
            .is_none()
    );
    assert!(
        stretch
            .add_tangent(&u, &mut Coo::new(model.free, model.free), 0.1)
            .is_none()
    );
}

#[test]
fn spatial_refinement_converges_to_independent_sinusoidal_stretch_energy() {
    let length = 1.5;
    let extension = 0.001_f64;
    let amplitude = 0.002;
    let wave = std::f64::consts::PI / length;
    let e0 = extension + 0.5 * extension * extension;
    let slope_squared = (amplitude * wave).powi(2);
    let ea = section().young_pa * section().area_m2;
    let expected = 0.5
        * ea
        * length
        * (extension.powi(3)
            + 0.25 * extension.powi(4)
            + 0.5 * e0 * slope_squared
            + (3. / 32.) * slope_squared.powi(2));
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut nodes: Vec<_> = (0..=count)
            .map(|i| [length * i as f64 / count as f64, 0., 0.])
            .collect();
        nodes.push([0., 1., 0.]);
        let mesh =
            ShellMesh::new(nodes, (0..count).map(|i| [i, i + 1, count + 1]).collect()).unwrap();
        let plate = PlateSection::isotropic(2e6, 0.25, 0.004, 400.).unwrap();
        let beam_chain: Vec<_> = (0..count)
            .map(|i| ShellBeam {
                nodes: [i, i + 1],
                offsets_m: [[0.; 3]; 2],
                section_z: [0., 0., 1.],
                section: section(),
            })
            .collect();
        let model = assemble_stiffened_shell(
            &mesh,
            &vec![plate; count],
            &[],
            ShellSupport::Free,
            &beam_chain,
        )
        .unwrap();
        let mut u = vec![0.; model.free];
        for i in 0..=count {
            let x = mesh.nodes[i][0];
            u[6 * i] = extension * x;
            u[6 * i + 2] = amplitude * (wave * x).sin();
            u[6 * i + 4] = -amplitude * wave * (wave * x).cos();
        }
        let mut force = vec![0.; model.free];
        let energy: f64 = beam_chain
            .iter()
            .map(|beam| {
                BeamStretch::new(&mesh, &model, beam)
                    .unwrap()
                    .add_force_energy(&u, &mut force, 0.1)
                    .unwrap()
            })
            .sum();
        errors.push((energy / expected - 1.).abs());
    }
    assert!(errors[2] < 1e-5, "energy errors: {errors:?}");
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (3.8..4.2).contains(&order),
            "energy errors: {errors:?}, order {order}"
        );
    }
}

#[test]
fn reinforced_shell_equilibrates_the_independent_quartic_beam_energy() {
    let mut nodes = Vec::new();
    let mut tris = Vec::new();
    let mut fixed = Vec::new();
    for j in 0..3 {
        for i in 0..3 {
            nodes.push([i as f64 * 0.25, j as f64 * 0.25, 0.]);
            if i != 1 || j != 1 {
                fixed.push(3 * j + i);
            }
        }
    }
    for j in 0..2 {
        for i in 0..2 {
            let a = 3 * j + i;
            tris.extend([[a, a + 1, a + 4], [a, a + 4, a + 3]]);
        }
    }
    let mesh = ShellMesh::new(nodes, tris).unwrap();
    let plate = PlateSection::isotropic(2e6, 0.25, 0.004, 400.).unwrap();
    let sections = vec![plate; mesh.tris.len()];
    let beams = [[3, 4], [4, 5]].map(|nodes| ShellBeam {
        nodes,
        offsets_m: [[0.; 3]; 2],
        section_z: [0., 0., 1.],
        section: section(),
    });
    let base =
        assemble_stiffened_shell(&mesh, &sections, &fixed, ShellSupport::Clamped, &beams).unwrap();
    let w = base.dof_map[26].unwrap();
    let l: f64 = 0.25;
    let ea = section().young_pa * section().area_m2;
    // Two cubic beam fields with fixed endpoint slopes have total axial
    // energy (18/35)*EA*delta^4/L^3, by exact integration of (w')^4.
    let beam_coefficient = (18. / 35.) * ea / l.powi(3);
    let mut shell_coefficient = 0.;
    let a = 2e6 * 0.004 / (1. - 0.25_f64.powi(2));
    for (e, triangle) in mesh.tris.iter().enumerate() {
        if let Some(local) = triangle.iter().position(|&n| n == 4) {
            let geom = mesh.facet(e).unwrap();
            let [gx, gy] = geom.gradient[local];
            let q = [0.5 * gx * gx, 0.5 * gy * gy, gx * gy];
            shell_coefficient += 0.5
                * geom.area_m2
                * a
                * (q[0] * q[0] + q[1] * q[1] + 0.5 * q[0] * q[1] + 0.375 * q[2] * q[2]);
        }
    }
    let coefficient = shell_coefficient + beam_coefficient;
    for delta in [-0.004_f64, 0., 0.004] {
        let mut loads = vec![0.; 6 * mesh.nodes.len()];
        loads[26] = base.k.get(w, w) * delta + 4. * coefficient * delta.powi(3);
        let loaded = equilibrate_stiffened_shell(
            &mesh,
            &sections,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &loads,
            PreloadOptions {
                relative_force_tolerance: 1e-11,
                ..PreloadOptions::default()
            },
            &mut reference_newton,
        )
        .unwrap();
        near(loaded.displacement[26], delta, 1e-11);
        near(
            loaded.stored_energy_j,
            0.5 * base.k.get(w, w) * delta * delta + coefficient * delta.powi(4),
            1e-10,
        );
        near(
            loaded.model.k.get(w, w),
            base.k.get(w, w) + 12. * coefficient * delta * delta,
            1e-9,
        );
        if delta != 0. {
            assert!(12. * beam_coefficient * delta * delta > 0.1 * base.k.get(w, w));
        }
        for row in 0..base.free {
            assert_eq!(loaded.model.m.row(row), base.m.row(row));
            if delta == 0. {
                assert_eq!(loaded.model.k.row(row), base.k.row(row));
            }
        }
    }
    let excessive = vec![beams[0]; MAX_BONDED_BEAMS + 1];
    assert!(
        equilibrate_stiffened_shell(
            &mesh,
            &sections,
            &fixed,
            ShellSupport::Clamped,
            &excessive,
            &vec![0.; 54],
            PreloadOptions::default(),
            &mut reference_newton
        )
        .is_err()
    );
}
