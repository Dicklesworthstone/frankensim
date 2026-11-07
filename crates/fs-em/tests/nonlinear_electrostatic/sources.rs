//! Independent Gauss, material-storage and radial-moment source regressions.
//! Impressed polarization is prescribed separately from the original D(E) law.
use super::*;
use fs_em::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, ElectrostaticError,
    PlanarElectrostaticSources,
    axisymmetric::{AxisymmetricElectrostaticSources, AxisymmetricElectrostatics},
    nonlinear::{
        NonlinearElectricField, NonlinearElectrostaticError,
        axisymmetric::AxisymmetricNonlinearElectrostatics,
    },
};

fn slab_problem() -> PlanarNonlinearElectrostatics {
    PlanarNonlinearElectrostatics::new(
        slab_mesh(),
        vec![DielectricMaterial::Reversible(nonlinear()); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn assert_invalid<T>(result: Result<T, NonlinearElectrostaticError>) {
    assert!(matches!(
        result,
        Err(NonlinearElectrostaticError::Electric(
            ElectrostaticError::Invalid(_)
        ))
    ));
}

fn assert_source_transport_is_not_invented(field: &NonlinearElectricField<'_>) {
    let motion = [
        [0.0, 0.0],
        [1e-6, 0.0],
        [2e-6, 0.0],
        [0.0, 0.0],
        [1e-6, 0.0],
        [2e-6, 0.0],
    ];
    assert_invalid(field.shape_forces(&mut || false));
    assert_invalid(field.virtual_work(&motion, 0.0, &mut || false));
}

#[test]
fn internal_sheet_drives_nonlinear_jump_storage_and_fixed_source_tangent() {
    let problem = slab_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let sheets = [ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: 5.0 * D0,
    }];
    let sources = PlanarElectrostaticSources {
        charge_sheets: &sheets,
        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
    };
    let field = electrodes
        .solve_with_sources(
            sources,
            &[ElectrodeControl::Voltage(0.0); 2],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    // Each half has |E|=2 E0 and |D|=2.5 D0; their jump is the sheet.
    for (cell, (e, d)) in field
        .electric_field_v_m()
        .iter()
        .zip(field.displacement_c_m2())
        .enumerate()
    {
        let sign = if cell < 2 { -1.0 } else { 1.0 };
        near(e[0], sign * 2.0 * E0, E0);
        near(e[1], 0.0, E0);
        near(d[0], sign * 2.5 * D0, D0);
        near(d[1], 0.0, D0);
    }
    for (actual, expected) in field
        .relative_voltage_v()
        .iter()
        .zip([0.0, V0, 0.0, 0.0, V0, 0.0])
    {
        near(*actual, expected, V0);
    }
    for &charge in field.terminal_charge_c() {
        near(charge, -2.5 * Q0, Q0);
    }
    assert_eq!(field.free_charge_sheets(), &sheets);
    assert!(field.displacement_traces().is_empty());
    assert_eq!(
        field.displacement_c_m2(),
        field.constitutive_displacement_c_m2()
    );
    near(field.sheet_charge_c(), 5.0 * Q0, Q0);
    near(field.volume_charge_c(), 0.0, Q0);
    near(field.outward_displacement_c(), 0.0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    near(
        field.energy_j() + field.coenergy_j(),
        field.sheet_charge_c() * V0,
        U0,
    );
    near(field.energy_identity_residual_j(), 0.0, U0);
    near(field.charge_balance_residual_c(), 0.0, Q0);
    assert_source_transport_is_not_invented(&field);

    let tangent = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    let capacitance = 0.5 * EPSILON0 * AREA / GAP;
    for (&actual, expected) in
        tangent
            .matrix_f()
            .iter()
            .zip([capacitance, -capacitance, -capacitance, capacitance])
    {
        near(actual, expected, capacitance);
    }
    let charge_controlled = electrodes
        .solve_with_sources(
            sources,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(-2.5 * Q0),
            ],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&actual, &expected) in charge_controlled
        .relative_voltage_v()
        .iter()
        .zip(field.relative_voltage_v())
    {
        near(actual, expected, V0);
    }
    near(charge_controlled.terminal_charge_c()[0], -2.5 * Q0, Q0);
    near(charge_controlled.energy_j(), field.energy_j(), U0);
    near(charge_controlled.coenergy_j(), field.coenergy_j(), U0);
}

#[test]
fn impressed_polarization_open_circuit_uses_total_displacement_and_original_storage() {
    let problem = slab_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3]], &mut || false)
        .unwrap();
    let polarization = [[2.5 * D0, 0.0]; 4];
    let field = electrodes
        .solve_with_sources(
            PlanarElectrostaticSources {
                impressed_polarization_c_m2: Some(&polarization),
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        field.impressed_polarization_c_m2(),
        Some(polarization.as_slice())
    );
    for cell in 0..4 {
        near(field.electric_field_v_m()[cell][0], -2.0 * E0, E0);
        near(field.electric_field_v_m()[cell][1], 0.0, E0);
        near(
            field.constitutive_displacement_c_m2()[cell][0],
            -2.5 * D0,
            D0,
        );
        for &component in &field.displacement_c_m2()[cell] {
            near(component, 0.0, D0);
        }
    }
    for (actual, expected) in
        field
            .relative_voltage_v()
            .iter()
            .zip([0.0, V0, 2.0 * V0, 0.0, V0, 2.0 * V0])
    {
        near(*actual, expected, V0);
    }
    near(field.terminal_charge_c()[0], 0.0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    near(field.polarization_pairing_j(), -5.0 * U0, U0);
    near(
        field.energy_j() + field.coenergy_j() + field.polarization_pairing_j(),
        0.0,
        U0,
    );
    near(field.energy_identity_residual_j(), 0.0, U0);
    assert!(field.energy_j() > 0.0);
    assert_source_transport_is_not_invented(&field);
}

#[test]
fn outward_displacement_drives_nonlinear_field_with_physical_flux_sign() {
    let problem = slab_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3]], &mut || false)
        .unwrap();
    let traces = [ElectrostaticDisplacementTrace {
        vertices: [2, 5],
        outward_displacement_c_m2: 2.5 * D0,
    }];
    let field = electrodes
        .solve_with_sources(
            PlanarElectrostaticSources {
                displacement_traces: &traces,
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (e, d) in field
        .electric_field_v_m()
        .iter()
        .zip(field.displacement_c_m2())
    {
        near(e[0], 2.0 * E0, E0);
        near(e[1], 0.0, E0);
        near(d[0], 2.5 * D0, D0);
        near(d[1], 0.0, D0);
    }
    near(field.relative_voltage_v()[2], -2.0 * V0, V0);
    near(field.relative_voltage_v()[5], -2.0 * V0, V0);
    assert_eq!(field.displacement_traces(), &traces);
    near(field.outward_displacement_c(), 2.5 * Q0, Q0);
    near(
        field.terminal_charge_c()[0],
        field.outward_displacement_c(),
        Q0,
    );
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    near(
        field.energy_j() + field.coenergy_j(),
        -field.outward_displacement_c() * field.relative_voltage_v()[2],
        U0,
    );
    near(field.charge_balance_residual_c(), 0.0, Q0);
    assert_source_transport_is_not_invented(&field);
}

const ANGLE: f64 = 0.7;

fn radial_mesh() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "nonlinear radial source moments",
        vec![
            [0.0, 0.0],
            [WIDTH, 0.0],
            [0.0, GAP / 2.0],
            [WIDTH, GAP / 2.0],
            [0.0, GAP],
            [WIDTH, GAP],
        ],
        vec![[0, 1, 3], [0, 3, 2], [2, 3, 5], [2, 5, 4]],
        Metric2::axisymmetric(ANGLE).unwrap(),
    )
    .unwrap()
}

#[test]
fn radial_owner_matches_linear_sources_and_independent_exact_radial_moments() {
    let material = linear(EPSILON0);
    let original = AxisymmetricElectrostatics::new(
        radial_mesh(),
        vec![material.clone(); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let problem = AxisymmetricNonlinearElectrostatics::new(
        radial_mesh(),
        vec![DielectricMaterial::isotropic(material); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 1]], &mut || false)
        .unwrap();
    let rho = [0.7, -0.2, 0.9, 0.3].map(|value| value * D0 / GAP);
    let polarization = [[0.0, 0.15 * D0]; 4];
    let sheets = [ElectrostaticChargeSheet {
        vertices: [2, 3],
        free_charge_density_c_m2: 0.35 * D0,
    }];
    let traces = [ElectrostaticDisplacementTrace {
        vertices: [4, 5],
        outward_displacement_c_m2: 0.4 * D0,
    }];
    let sources = AxisymmetricElectrostaticSources {
        free_charge_density_c_m3: &rho,
        impressed_polarization_c_m2: Some(&polarization),
        charge_sheets: &sheets,
        displacement_traces: &traces,
    };
    let reference = original
        .solve_with_sources(
            sources,
            &[(0, 0.0), (1, 0.0)],
            ElectrostaticControls {
                linear: linear_controls(),
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    let field = electrodes
        .solve_with_sources(
            sources,
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&actual, &expected) in field
        .relative_voltage_v()
        .iter()
        .zip(reference.relative_voltage_v())
    {
        near(actual, expected, V0);
    }
    for (actual, expected) in field
        .electric_field_v_m()
        .iter()
        .zip(reference.electric_field_v_m())
    {
        for axis in 0..2 {
            near(actual[axis], expected[axis], E0);
        }
    }
    for (actual, expected) in field
        .displacement_c_m2()
        .iter()
        .zip(reference.displacement_c_m2())
    {
        for axis in 0..2 {
            near(actual[axis], expected[axis], D0);
        }
    }
    near(field.energy_j(), reference.energy_j(), U0);
    near(field.coenergy_j(), reference.energy_j(), U0);
    near(
        field.terminal_charge_c()[0],
        reference.boundary_charge_c(),
        Q0,
    );
    assert_eq!(
        field.impressed_polarization_c_m2(),
        Some(polarization.as_slice())
    );
    assert_eq!(field.free_charge_sheets(), &sheets);
    assert_eq!(field.displacement_traces(), &traces);

    // Exact revolved volumes; nodal volume shares differ from equal thirds.
    // A radial edge from 0 to R has surface shares (1/3, 2/3), not (1/2, 1/2).
    let area = 0.5 * ANGLE * WIDTH * WIDTH;
    let volumes = [
        area * GAP / 3.0,
        area * GAP / 6.0,
        area * GAP / 3.0,
        area * GAP / 6.0,
    ];
    let shares = [
        [0.25, 0.375, 0.375],
        [0.25, 0.5, 0.25],
        [0.25, 0.375, 0.375],
        [0.25, 0.5, 0.25],
    ];
    let v = field.relative_voltage_v();
    let mut volume_work = 0.0;
    let mut expected_polarization_pairing = 0.0;
    let mut expected_storage = 0.0;
    for cell in 0..4 {
        let vertices = problem.mesh().faces()[cell];
        let mean_v: f64 = vertices
            .iter()
            .zip(shares[cell])
            .map(|(&vertex, weight)| v[vertex as usize] * weight)
            .sum();
        volume_work += rho[cell] * volumes[cell] * mean_v;
        let e = field.electric_field_v_m()[cell];
        expected_polarization_pairing += volumes[cell] * polarization[cell][1] * e[1];
        expected_storage += 0.5 * EPSILON0 * volumes[cell] * (e[0] * e[0] + e[1] * e[1]);
        for axis in 0..2 {
            near(
                field.constitutive_displacement_c_m2()[cell][axis],
                EPSILON0 * e[axis],
                D0,
            );
        }
    }
    let volume_charge: f64 = rho
        .iter()
        .zip(volumes)
        .map(|(rho, volume)| rho * volume)
        .sum();
    let sheet_charge = 0.35 * D0 * area;
    let outward_charge = 0.4 * D0 * area;
    let sheet_work = sheet_charge * (v[2] / 3.0 + 2.0 * v[3] / 3.0);
    let outward_work = outward_charge * (v[4] / 3.0 + 2.0 * v[5] / 3.0);
    near(field.volume_charge_c(), volume_charge, Q0);
    near(field.sheet_charge_c(), sheet_charge, Q0);
    near(field.outward_displacement_c(), outward_charge, Q0);
    near(
        field.terminal_charge_c()[0],
        outward_charge - volume_charge - sheet_charge,
        Q0,
    );
    near(field.energy_j(), expected_storage, U0);
    near(
        field.polarization_pairing_j(),
        expected_polarization_pairing,
        U0,
    );
    near(
        field.energy_j() + field.coenergy_j(),
        volume_work + sheet_work - outward_work - expected_polarization_pairing,
        U0,
    );
    near(field.energy_identity_residual_j(), 0.0, U0);
    near(field.charge_balance_residual_c(), 0.0, Q0);
}

#[test]
fn radial_virtual_work_refuses_each_nonzero_extra_source_without_transport_law() {
    let problem = AxisymmetricNonlinearElectrostatics::new(
        radial_mesh(),
        vec![DielectricMaterial::isotropic(linear(EPSILON0)); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 1]], &mut || false)
        .unwrap();
    let polarization = [[0.0, 0.15 * D0]; 4];
    let sheets = [ElectrostaticChargeSheet {
        vertices: [2, 3],
        free_charge_density_c_m2: 0.35 * D0,
    }];
    let traces = [ElectrostaticDisplacementTrace {
        vertices: [4, 5],
        outward_displacement_c_m2: 0.4 * D0,
    }];
    let empty = AxisymmetricElectrostaticSources::volume_charge(&[0.0; 4]);
    let motion: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|x| [0.0, x[1] * 1e-3])
        .collect();
    for sources in [
        AxisymmetricElectrostaticSources {
            impressed_polarization_c_m2: Some(&polarization),
            ..empty
        },
        AxisymmetricElectrostaticSources {
            charge_sheets: &sheets,
            ..empty
        },
        AxisymmetricElectrostaticSources {
            displacement_traces: &traces,
            ..empty
        },
    ] {
        let field = electrodes
            .solve_with_sources(
                sources,
                &[ElectrodeControl::Voltage(0.0)],
                &[0.0; 6],
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_invalid(field.virtual_work(&motion, &mut || false));
    }
}

#[test]
fn source_surfaces_cannot_replace_voltage_or_floating_electrode_controls() {
    let problem = slab_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let traces = [ElectrostaticDisplacementTrace {
        vertices: [2, 5],
        outward_displacement_c_m2: D0,
    }];
    for right in [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Charge(0.0),
    ] {
        assert_invalid(electrodes.solve_with_sources(
            PlanarElectrostaticSources {
                displacement_traces: &traces,
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[ElectrodeControl::Voltage(0.0), right],
            &[0.0; 6],
            controls(),
            &mut || false,
        ));
    }
    for edges in [vec![[0, 3]], vec![[0, 5]], vec![[1, 4], [4, 1]]] {
        let sheets: Vec<_> = edges
            .into_iter()
            .map(|vertices| ElectrostaticChargeSheet {
                vertices,
                free_charge_density_c_m2: D0,
            })
            .collect();
        assert_invalid(electrodes.solve_with_sources(
            PlanarElectrostaticSources {
                charge_sheets: &sheets,
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[ElectrodeControl::Voltage(0.0); 2],
            &[0.0; 6],
            controls(),
            &mut || false,
        ));
    }
}

#[test]
fn source_solve_cancellation_never_returns_a_partial_field_or_changes_reuse() {
    let problem = slab_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let sheets = [ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: 5.0 * D0,
    }];
    let sources = PlanarElectrostaticSources {
        charge_sheets: &sheets,
        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
    };
    let mut polls = 0;
    let baseline = electrodes
        .solve_with_sources(
            sources,
            &[ElectrodeControl::Voltage(0.0); 2],
            &[0.0; 6],
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        let result = electrodes.solve_with_sources(
            sources,
            &[ElectrodeControl::Voltage(0.0); 2],
            &[0.0; 6],
            controls(),
            &mut || {
                seen += 1;
                seen >= stop
            },
        );
        assert!(matches!(
            result,
            Err(NonlinearElectrostaticError::Electric(
                ElectrostaticError::Cancelled
            ))
        ));
        assert_eq!(seen, stop);
    }
    let repeated = electrodes
        .solve_with_sources(
            sources,
            &[ElectrodeControl::Voltage(0.0); 2],
            &[0.0; 6],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(repeated.relative_voltage_v(), baseline.relative_voltage_v());
    assert_eq!(repeated.terminal_charge_c(), baseline.terminal_charge_c());
    assert_eq!(repeated.energy_j(), baseline.energy_j());
}
