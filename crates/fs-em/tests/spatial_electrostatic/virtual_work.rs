//! G1/G3: original 3D field energy, physical shape work and source boundaries.
//! Synthetic material fixtures establish numerical behavior, not qualification.

use super::*;
use fs_em::spatial::electrostatic::{ElectrostaticChargeSheet, SpatialElectrostaticSources};

fn displaced(
    original: &SpatialElectrostatics,
    motion: &[[f64; 3]],
    step: f64,
) -> SpatialElectrostatics {
    SpatialElectrostatics::new(
        original
            .positions_m()
            .iter()
            .zip(motion)
            .map(|(x, u)| std::array::from_fn(|axis| step.mul_add(u[axis], x[axis])))
            .collect(),
        original.mesh().tets.clone(),
        original.materials().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn g1_spatial_plate_work_preserves_attraction_scaling_rigid_modes_and_orientation() {
    let epsilon = 3e-11;
    for reversed in [false, true] {
        let (positions, mut tets) = cube(&[0.0, 0.25, 1.0], 1, 1);
        if reversed {
            for tet in &mut tets {
                tet.swap(1, 2);
            }
        }
        let problem = SpatialElectrostatics::new(
            positions,
            tets,
            vec![dielectric(epsilon); 12],
            limits(),
            &mut || false,
        )
        .unwrap();
        let terminals = terminal_faces(&problem);
        let capacitance = problem
            .capacitance(&terminals[0], &terminals[1], controls(), &mut || false)
            .unwrap();
        for voltage in [4.0, -4.0, 8.0] {
            let field = capacitance
                .field_at_voltage(voltage, controls(), &mut || false)
                .unwrap();
            let energy = 0.5 * epsilon * voltage * voltage;
            near(field.energy_j(), energy);
            for axis in 0..3 {
                let motion: Vec<_> = problem
                    .positions_m()
                    .iter()
                    .map(|x| {
                        let mut u = [0.0; 3];
                        u[axis] = x[axis];
                        u
                    })
                    .collect();
                near(
                    field.virtual_work(&motion, &mut || false).unwrap(),
                    if axis == 0 { -energy } else { energy },
                );
            }
            near(
                field
                    .virtual_work(problem.positions_m(), &mut || false)
                    .unwrap(),
                energy,
            );
            // Three independent infinitesimal rotations; the original cell
            // tensors are isotropic, so no frame rotation supplies work.
            for axis in 0..3 {
                let rotation: Vec<_> = problem
                    .positions_m()
                    .iter()
                    .map(|&[x, y, z]| match axis {
                        0 => [0.0, -z, y],
                        1 => [z, 0.0, -x],
                        _ => [-y, x, 0.0],
                    })
                    .collect();
                assert!(
                    field.virtual_work(&rotation, &mut || false).unwrap().abs() < 1e-12 * energy
                );
            }
            assert_eq!(
                field
                    .virtual_work(&[[2_f64.powi(50); 3]; 12], &mut || false)
                    .unwrap(),
                0.0
            );
        }
    }
}

#[test]
fn g1_spatial_interface_work_matches_true_rebuilt_fields_at_fixed_voltage_and_charge() {
    let (positions, tets) = cube(&[0.0, 0.25, 1.0], 1, 1);
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        (0..12)
            .map(|cell| dielectric([2e-11, 5e-11][cell / 6]))
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals = terminal_faces(&problem);
    let capacitance = problem
        .capacitance(&terminals[0], &terminals[1], controls(), &mut || false)
        .unwrap();
    let voltage = 7.0;
    let charge = capacitance.capacitance_f() * voltage;
    let field = capacitance
        .field_at_voltage(voltage, controls(), &mut || false)
        .unwrap();
    let motion: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|&[x, y, z]| {
            [
                0.13 * x + 0.05 * y * z,
                0.09 * y + 0.07 * x * z,
                -0.04 * z + 0.03 * x * y,
            ]
        })
        .collect();
    let work = field.virtual_work(&motion, &mut || false).unwrap();
    let h = 1e-5;
    let mut energies = [[0.0; 2]; 2];
    for (i, step) in [-h, h].into_iter().enumerate() {
        let moved = displaced(&problem, &motion, step);
        let c = moved
            .capacitance(&terminals[0], &terminals[1], controls(), &mut || false)
            .unwrap();
        energies[0][i] = c
            .field_at_voltage(voltage, controls(), &mut || false)
            .unwrap()
            .energy_j();
        energies[1][i] = c
            .field_at_charge(charge, controls(), &mut || false)
            .unwrap()
            .energy_j();
    }
    let held_voltage_derivative = (energies[0][1] - energies[0][0]) / (2.0 * h);
    let held_charge_derivative = (energies[1][1] - energies[1][0]) / (2.0 * h);
    let tolerance = 5e-7 * work.abs().max(field.energy_j() * 1e-3);
    assert!(
        (work - held_voltage_derivative).abs() < tolerance,
        "{work} != dU/dt|V = {held_voltage_derivative}"
    );
    assert!(
        (work + held_charge_derivative).abs() < tolerance,
        "{work} != -dU/dt|Q = {}",
        -held_charge_derivative
    );
}

#[test]
fn g3_spatial_unowned_sources_tensor_motion_bad_displacement_and_cancellation_refuse() {
    let (positions, tets) = cube(&[0.0, 0.25, 1.0], 1, 1);
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(3e-11); 12],
        limits(),
        &mut || false,
    )
    .unwrap();
    let motion = problem.positions_m();
    let zero = [0.0; 12];
    let rho = [1e-11; 12];
    let polarization = [[2e-11, -1e-11, 3e-11]; 12];
    let face = problem
        .mesh()
        .faces
        .iter()
        .position(|f| {
            f.iter()
                .all(|&v| problem.positions_m()[v as usize][0] == 0.25)
        })
        .unwrap();
    let sheet = [ElectrostaticChargeSheet {
        face,
        free_charge_density_c_m2: 1e-11,
    }];
    for source in 0..3 {
        let field = problem
            .solve_with_sources(
                SpatialElectrostaticSources {
                    free_charge_density_c_m3: if source == 0 { &rho } else { &zero },
                    impressed_polarization_c_m2: (source == 1).then_some(polarization.as_slice()),
                    sheets: if source == 2 { &sheet } else { &[] },
                },
                &end_voltages(&problem, 0.0, 1.0),
                controls(),
                &mut || false,
            )
            .unwrap();
        assert!(matches!(
            field.virtual_work(motion, &mut || false),
            Err(SpatialElectrostaticError::Invalid(_))
        ));
    }
    let material = source(&PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.map(|name| (name, 3e-11)));
    let tensor = SpatialElectrostatics::new_with_laws(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        vec![material; 12],
        vec![
            SpatialDielectricLaw::orthotropic([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
                .unwrap();
            12
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = tensor
        .solve(
            &zero,
            &end_voltages(&tensor, 0.0, 1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        field.virtual_work(motion, &mut || false),
        Err(SpatialElectrostaticError::Invalid(_))
    ));
    let field = problem
        .solve(
            &zero,
            &end_voltages(&problem, 0.0, 1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        field.virtual_work(&motion[..11], &mut || false),
        Err(SpatialElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        field.virtual_work(&[[f64::INFINITY; 3]; 12], &mut || false),
        Err(SpatialElectrostaticError::Invalid(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        field.virtual_work(motion, &mut || {
            polls += 1;
            polls == 5
        }),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    assert!(field.virtual_work(motion, &mut || false).is_ok());
}

#[test]
fn g3_cancelling_nodal_weak_loads_do_not_hide_imposed_displacement_work() {
    let (positions, tets) = cube(&[0.0, 1.0], 1, 1);
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(3e-11); 6],
        limits(),
        &mut || false,
    )
    .unwrap();
    // Equal-area face sources cancel at EVERY vertex. Their original nonzero
    // normal traces still have unowned geometric source work and must refuse.
    let displacement = [
        ([0, 4, 6], 1.0),
        ([4, 6, 7], -1.0),
        ([0, 2, 6], -1.0),
        ([2, 6, 7], 1.0),
    ]
    .into_iter()
    .map(|(vertices, sign)| ElectrostaticDisplacementFace {
        face: problem
            .mesh()
            .faces
            .iter()
            .position(|f| *f == vertices)
            .unwrap(),
        outward_displacement_c_m2: sign * 3e-11,
    })
    .collect();
    let boundary = SpatialElectrostaticBoundary {
        voltage: vec![],
        displacement,
    };
    let field = problem
        .solve(&[0.0; 6], &boundary, controls(), &mut || false)
        .unwrap();
    assert!(
        field
            .nodal_applied_surface_charge_c()
            .iter()
            .all(|&q| q == 0.0)
    );
    assert_eq!(field.energy_j(), 0.0);
    assert!(matches!(
        field.virtual_work(problem.positions_m(), &mut || false),
        Err(SpatialElectrostaticError::Invalid(_))
    ));
}
