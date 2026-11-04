//! G1/G3: physical energy derivatives and source-domain refusals on native P1
//! fields. Material values are synthetic analytic fixtures, not measured data.

use super::*;
use fs_em::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, PlanarElectrostaticSources,
};

fn displaced(
    original: &PlanarElectrostatics,
    motion: &[[f64; 2]],
    step: f64,
) -> PlanarElectrostatics {
    let mesh = original.mesh();
    let positions = mesh
        .vertices()
        .iter()
        .zip(motion)
        .map(|(x, u)| std::array::from_fn(|axis| step.mul_add(u[axis], x[axis])))
        .collect();
    PlanarElectrostatics::new(
        TriComplex2::from_indexed_triangles(
            "displaced original dielectric",
            positions,
            mesh.faces().to_vec(),
            mesh.metric(),
        )
        .unwrap(),
        original.materials().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn g1_plate_work_has_the_physical_sign_voltage_scaling_and_original_depth() {
    let epsilon = 3e-11;
    for (voltage, depth) in [(4.0, 0.2), (-4.0, 0.2), (8.0, 0.2), (4.0, 0.6)] {
        let problem = PlanarElectrostatics::new(
            layered(0.5, 1.0, depth),
            vec![dielectric(epsilon); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = problem
            .solve(&[0.0; 4], &common::ends(voltage), controls(), &mut || false)
            .unwrap();
        let energy = 0.5 * epsilon * depth * voltage * voltage;
        near(field.energy_j(), energy);
        let gap: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|x| [x[0], 0.0])
            .collect();
        let height: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|x| [0.0, x[1]])
            .collect();
        // Positive gap displacement separates the electrodes: force attracts.
        near(field.virtual_work(&gap, &mut || false).unwrap(), -energy);
        near(field.virtual_work(&height, &mut || false).unwrap(), energy);
        let rotation: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|x| [-x[1], x[0]])
            .collect();
        assert!(field.virtual_work(&rotation, &mut || false).unwrap().abs() < 1e-12 * energy);
        assert!(
            field
                .virtual_work(problem.mesh().vertices(), &mut || false)
                .unwrap()
                .abs()
                < 1e-12 * energy,
            "in-plane uniform scaling at fixed depth leaves C unchanged"
        );
        assert_eq!(
            field
                .virtual_work(&[[2_f64.powi(50), -2_f64.powi(50)]; 6], &mut || false)
                .unwrap(),
            0.0
        );
    }
}

#[test]
fn g1_virtual_work_matches_rebuilt_interface_fields_at_fixed_voltage_and_fixed_charge() {
    let problem = slab(dielectric(2e-11), dielectric(5e-11));
    let voltage = 7.0;
    let field = problem
        .solve(&[0.0; 4], &common::ends(voltage), controls(), &mut || false)
        .unwrap();
    let capacitance = problem
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let charge = capacitance.capacitance_f() * voltage;
    // A nonuniform deformation also moves the dielectric interface. The
    // perturbed fields must equilibrate; merely transporting E is insufficient.
    let motion: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|&[x, y]| [0.12 * x + 0.17 * x * y, -0.09 * y + 0.05 * x * (1.0 - x)])
        .collect();
    let work = field.virtual_work(&motion, &mut || false).unwrap();
    let h = 1e-5;
    let mut energies = [[0.0; 2]; 2];
    for (i, step) in [-h, h].into_iter().enumerate() {
        let moved = displaced(&problem, &motion, step);
        let c = moved
            .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        for (control, voltage) in [voltage, charge / c.capacitance_f()]
            .into_iter()
            .enumerate()
        {
            let actual = moved
                .solve(&[0.0; 4], &common::ends(voltage), controls(), &mut || false)
                .unwrap();
            energies[control][i] = actual.energy_j();
        }
    }
    let held_voltage_derivative = (energies[0][1] - energies[0][0]) / (2.0 * h);
    let held_charge_derivative = (energies[1][1] - energies[1][0]) / (2.0 * h);
    let tolerance = 3e-7 * work.abs().max(field.energy_j() * 1e-3);
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
fn g3_unowned_source_work_tensor_transport_and_invalid_displacements_refuse() {
    let problem = slab(dielectric(3e-11), dielectric(3e-11));
    let motion = problem.mesh().vertices();
    let zero = [0.0; 4];
    let rho = [1e-11; 4];
    let polarization = [[2e-11, -1e-11]; 4];
    let sheet = [ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: 1e-11,
    }];
    let trace = [ElectrostaticDisplacementTrace {
        vertices: [0, 1],
        outward_displacement_c_m2: 1e-11,
    }];
    for source in 0..4 {
        let field = problem
            .solve_with_sources(
                PlanarElectrostaticSources {
                    free_charge_density_c_m3: if source == 0 { &rho } else { &zero },
                    impressed_polarization_c_m2: (source == 1).then_some(polarization.as_slice()),
                    charge_sheets: if source == 2 { &sheet } else { &[] },
                    displacement_traces: if source == 3 { &trace } else { &[] },
                },
                &common::ends(1.0),
                controls(),
                &mut || false,
            )
            .unwrap();
        assert!(matches!(
            field.virtual_work(motion, &mut || false),
            Err(ElectrostaticError::Invalid(_))
        ));
    }
    let tensor = PlanarElectrostatics::new_with_laws(
        layered(0.5, 1.0, 0.2),
        vec![anisotropic::material([3e-11; 2]); 4],
        vec![anisotropic::law(0.3); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = tensor
        .solve(&zero, &common::ends(1.0), controls(), &mut || false)
        .unwrap();
    assert!(
        matches!(
            field.virtual_work(motion, &mut || false),
            Err(ElectrostaticError::Invalid(_))
        ),
        "equal tensor eigenvalues do not declare a material-axis transport law"
    );
    let field = problem
        .solve(&zero, &common::ends(1.0), controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.virtual_work(&motion[..5], &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        field.virtual_work(&[[f64::NAN, 0.0]; 6], &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        field.virtual_work(motion, &mut || {
            polls += 1;
            polls == 4
        }),
        Err(ElectrostaticError::Cancelled)
    ));
    // Explicitly zero source records remain the same source-free problem.
    let field = problem
        .solve_with_sources(
            PlanarElectrostaticSources {
                free_charge_density_c_m3: &zero,
                impressed_polarization_c_m2: Some(&[[0.0; 2]; 4]),
                charge_sheets: &[ElectrostaticChargeSheet {
                    vertices: [1, 4],
                    free_charge_density_c_m2: 0.0,
                }],
                displacement_traces: &[ElectrostaticDisplacementTrace {
                    vertices: [0, 1],
                    outward_displacement_c_m2: 0.0,
                }],
            },
            &common::ends(1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(field.virtual_work(motion, &mut || false).is_ok());
}
