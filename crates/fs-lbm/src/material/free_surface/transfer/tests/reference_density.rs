//! Original density inventories survive deformation before native liquid ingress.

use super::*;

pub(super) fn compatible_reference(
    current: &[[f64; 3]],
    tetrahedra: &[[u32; 4]],
    rho: f64,
) -> (Vec<[f64; 3]>, Vec<[f64; 4]>) {
    let mut reference = current.to_vec();
    // A conforming, nonuniform affine deformation: the six tetrahedra have
    // different J even though the resulting liquid has one uniform density.
    reference[0][0] -= 0.001;
    let volume = |points: &[[f64; 3]], indices: [u32; 4]| {
        let p = indices.map(|i| points[i as usize]);
        dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0
    };
    let densities = tetrahedra
        .iter()
        .map(|&indices| [rho * volume(current, indices) / volume(&reference, indices); 4])
        .collect();
    (reference, densities)
}

pub(super) fn density_controls() -> TetReferenceDensityTransferControls {
    TetReferenceDensityTransferControls {
        maximum_density_defect_kg_m3: 1e-8,
        maximum_volume_defect_m3: 1e-19,
    }
}

#[test]
fn original_graded_density_pushes_forward_exact_mass_moments_and_p1_storage() {
    let reference = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let density = [3.0, 7.0, 14.0, 24.0];
    let map = [[2.0, 0.5, 0.0], [0.0, 3.0, 0.25], [0.0, 0.0, 0.5]];
    let offset = [1.0, 2.0, 3.0];
    let current = reference.map(|p| add(map.map(|row| dot(row, p)), offset));
    let velocity = [
        [0.0; 3],
        [1.0, 2.0, -1.0],
        [-2.0, 0.5, 0.0],
        [0.0, -1.0, 3.0],
    ];
    let controls = TetReferenceDensityTransferControls {
        // Inspect original inventory mathematics without publishing a native
        // flow. Actual ingress below uses the strict original-liquid bound.
        maximum_density_defect_kg_m3: 100.0,
        maximum_volume_defect_m3: 0.0,
    };
    let binding = density::Binding {
        field: TetrahedralLiquidReferenceDensity {
            positions_m: &reference,
            corner_density_kg_m3: &[density],
        },
        controls,
    };
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let tet = binding
            .tetrahedron(0, [0, 1, 2, 3], current, 0.5, 4.0, cx)
            .unwrap();
        let mut report = TetReferenceDensityTransferReport::default();
        tet.add_report(&mut report, [0.0; 3]).unwrap();
        close(report.reference.mass, 2.0, 2e-14);
        close(report.deformed.mass, report.reference.mass, 0.0);
        close(report.reference.volume, 1.0 / 6.0, 2e-16);
        close(report.deformed.volume, 0.5, 0.0);
        for i in 0..4 {
            close(tet.current_density[i], density[i] / 3.0, 2e-14);
        }
        // Independent Dirichlet integration: V prod(count!)/(3+degree)! * 3!.
        // Enumerate original barycentric monomials rather than reuse F moments.
        let integral = |indices: &[usize]| {
            let mut counts = [0usize; 4];
            for &i in indices {
                counts[i] += 1;
            }
            let factorial = |n: usize| (1..=n).product::<usize>() as f64;
            counts.into_iter().map(factorial).product::<f64>() / factorial(3 + indices.len())
        };
        let mut expected_first = [0.0; 3];
        let mut expected_second = [[0.0; 3]; 3];
        let mut expected_p = [0.0; 3];
        let mut expected_l = [0.0; 3];
        let mut expected_k = 0.0;
        let h = [10.0, 20.0, 40.0, 80.0];
        let mut expected_h = 0.0;
        for a in 0..4 {
            for b in 0..4 {
                let weight = density[a] * integral(&[a, b]);
                expected_first = add(expected_first, scale(current[b], weight));
                expected_p = add(expected_p, scale(velocity[b], weight));
                expected_h += weight * h[b];
                for c in 0..4 {
                    let cubic = density[a] * integral(&[a, b, c]);
                    expected_l = add(expected_l, scale(cross(current[b], velocity[c]), cubic));
                    expected_k += 0.5 * cubic * dot(velocity[b], velocity[c]);
                    for i in 0..3 {
                        for j in 0..3 {
                            expected_second[i][j] += cubic * current[b][i] * current[c][j];
                        }
                    }
                }
            }
        }
        let mut physical = LiquidTransferMoments::default();
        tet.add_source(
            &mut physical,
            velocity[0],
            field_gradient(current, velocity).unwrap(),
            [0.0; 3],
        );
        for i in 0..3 {
            close(report.deformed.first[i], expected_first[i], 5e-14);
            close(physical.momentum_kg_m_s[i], expected_p[i], 2e-14);
            close(physical.angular_momentum_kg_m2_s[i], expected_l[i], 5e-14);
            for j in 0..3 {
                close(report.deformed.second[i][j], expected_second[i][j], 2e-13);
            }
        }
        close(physical.kinetic_energy_j, expected_k, 2e-14);
        close(tet.specific_field_integral(h), expected_h, 2e-13);
        assert!((expected_h - physical.mass_kg * 37.5).abs() > 1.0);

        // Retained kinetic work uses the original cubic density form even
        // beside a large, unchanged transport velocity.
        let transport = [1e9, -2e9, 3e9];
        let mut inventory = kinetic::Inventory::new(transport);
        inventory.add_density_tetrahedron(
            tet.original.mass,
            density,
            velocity.map(|v| add(v, transport)),
        );
        let frame = inventory.in_frame(transport).unwrap();
        close(
            frame
                .relative_kinetic_energy_j
                .into_iter()
                .fold(Dd::ZERO, |a, b| a + b)
                .to_f64(),
            expected_k,
            2e-14,
        );
        for i in 0..3 {
            close(
                frame.relative_momentum_kg_m_s[i].to_f64(),
                expected_p[i],
                2e-14,
            );
        }
    });
}

#[test]
fn different_reference_densities_and_jacobians_enter_the_same_native_liquid() {
    let material = liquid(1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let (reference, density) = compatible_reference(&p, &t, 1000.0);
    assert!(density.iter().any(|rho| rho[0] != density[0][0]));
    let snapshot = TetrahedralLiquidState {
        positions_m: &p,
        tetrahedra: &t,
        velocities_m_s: &v,
    };
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&material, 0.001, 0.001, 0.1).unwrap();
        let (mut flow, report) = scaling
            .free_surface3_from_reference_tetrahedra(
                domain(),
                [0.0; 3],
                snapshot,
                TetrahedralLiquidReferenceDensity {
                    positions_m: &reference,
                    corner_density_kg_m3: &density,
                },
                density_controls(),
                limits(),
                cx,
            )
            .unwrap();
        let inventory = report.reference_density.unwrap();
        let mass = 1000.0 * 0.004_f64.powi(3);
        close(inventory.reference.mass, mass, 2e-18);
        close(inventory.deformed.mass, mass, 2e-18);
        close(report.source.mass_kg, inventory.reference.mass, 2e-18);
        close(report.target.mass_kg, mass, 2e-18);
        assert!(inventory.reference.volume > inventory.deformed.volume);
        assert!(inventory.reference.first[0] < inventory.deformed.first[0]);
        for i in 0..3 {
            close(inventory.deformed.first[i], mass * 0.004, 2e-20);
            for j in 0..3 {
                let expected = mass
                    * (0.004 * 0.004
                        + if i == j {
                            0.004_f64.powi(2) / 12.0
                        } else {
                            0.0
                        });
                close(inventory.deformed.second[i][j], expected, 2e-22);
            }
            close(
                report.target.momentum_kg_m_s[i],
                report.source.momentum_kg_m_s[i],
                2e-20,
            );
            close(
                report.target.angular_momentum_kg_m2_s[i],
                report.source.angular_momentum_kg_m2_s[i],
                2e-22,
            );
        }
        assert!(
            inventory.maximum_current_density_defect_kg_m3
                <= density_controls().maximum_density_defect_kg_m3
        );
        assert!(
            inventory.sum_absolute_cell_mass_defect_kg + inventory.overlap_mass_partition_defect_kg
                <= limits().mass_tolerance_kg
        );
        assert!(
            inventory.sum_absolute_cell_volume_defect_m3
                + inventory.overlap_volume_partition_defect_m3
                <= density_controls().maximum_volume_defect_m3
        );
        let actual_volume: f64 = report
            .components
            .iter()
            .flat_map(|piece| &piece.target_cells)
            .map(|&[x, y, z]| {
                flow.sample_at_si(x, y, z).unwrap().unwrap().fill_fraction * 0.001_f64.powi(3)
            })
            .sum();
        close(inventory.native_occupied_volume_m3, actual_volume, 2e-22);
        close(inventory.native_minus_deformed_volume_m3, 0.0, 2e-22);
        let original_mass = flow.ledger_mass_kg().unwrap();
        flow.step().unwrap();
        close(flow.ledger_mass_kg().unwrap(), original_mass, 1e-14);
    });
}

#[test]
fn density_changes_inversions_and_reference_work_exhaustion_refuse_before_publication() {
    let material = liquid(1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let (reference, density) = compatible_reference(&p, &t, 1000.0);
    let unchanged = (p.clone(), reference.clone(), density.clone(), v.clone());
    let snapshot = TetrahedralLiquidState {
        positions_m: &p,
        tetrahedra: &t,
        velocities_m_s: &v,
    };
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&material, 0.001, 0.001, 0.1).unwrap();
        let attempt = |points: &[[f64; 3]], values: &[[f64; 4]], controls, bounds| {
            scaling.free_surface3_from_reference_tetrahedra(
                domain(),
                [0.0; 3],
                snapshot,
                TetrahedralLiquidReferenceDensity {
                    positions_m: points,
                    corner_density_kg_m3: values,
                },
                controls,
                bounds,
                cx,
            )
        };
        let mut graded = density.clone();
        graded[0][0] *= 1.01;
        assert!(attempt(&reference, &graded, density_controls(), limits()).is_err());
        let wrong: Vec<_> = density
            .iter()
            .map(|rho| rho.map(|value| value * 1.01))
            .collect();
        // Even an intentionally loose local density bound cannot manufacture
        // the missing original mass by pressure or a changed geometric fill.
        let mut loose = density_controls();
        loose.maximum_density_defect_kg_m3 = 20.0;
        assert!(attempt(&reference, &wrong, loose, limits()).is_err());
        let inverted: Vec<_> = reference.iter().map(|p| [-p[0], p[1], p[2]]).collect();
        assert!(attempt(&inverted, &density, density_controls(), limits()).is_err());
        let mut exhausted = limits();
        exhausted.max_overlap_tests = 1;
        assert!(attempt(&reference, &density, density_controls(), exhausted).is_err());
        assert!(
            attempt(
                &reference,
                &density[..density.len() - 1],
                density_controls(),
                limits()
            )
            .is_err()
        );
    });
    assert_eq!((p, reference, density, v), unchanged);
}
