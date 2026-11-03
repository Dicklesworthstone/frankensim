//! G1/G3/G4: tetrahedral dielectric fields against analytic charge, energy,
//! capacitance and refinement results. Sources are synthetic numerical fixtures,
//! not experimentally qualified dielectric data.

#[path = "spatial_electrostatic/force.rs"]
mod force;

use fs_em::{
    electrostatic::ElectrostaticControls,
    spatial::electrostatic::{
        ElectrostaticDisplacementFace, ElectrostaticVoltageFace, SpatialElectrostaticBoundary,
        SpatialElectrostaticError, SpatialElectrostatics,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        spatial::{PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, SpatialDielectricLaw},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_solver::CheckedCgConfig;

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4096,
        max_dofs: 16_384,
        max_triplets: 147_456,
    }
}

fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 2000,
            max_corrections: 2,
        },
        charge_tolerance_c: 1e-18,
        energy_tolerance_j: 1e-18,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1e-25),
        "{actual} != {expected}"
    );
}

fn near_vector(actual: [f64; 3], expected: [f64; 3]) {
    let scale = expected.into_iter().map(f64::abs).fold(1e-25, f64::max);
    for i in 0..3 {
        assert!(
            (actual[i] - expected[i]).abs() <= 2e-8 * scale,
            "{actual:?} != {expected:?}"
        );
    }
}

fn source(properties: &[(&str, f64)]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, ELECTRIC_PERMITTIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: ELECTRIC_PERMITTIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured three-dimensional dielectric".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                name,
                ELECTRIC_PERMITTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn dielectric(epsilon: f64) -> ResolvedMaterialStatePoint {
    source(&[(ELECTRIC_PERMITTIVITY_PROPERTY, epsilon)])
}

fn cube(xs: &[f64], ny: usize, nz: usize) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut positions = Vec::new();
    for &x in xs {
        for y in 0..=ny {
            for z in 0..=nz {
                positions.push([x, y as f64 / ny as f64, z as f64 / nz as f64]);
            }
        }
    }
    let index = |x: usize, y: usize, z: usize| ((x * (ny + 1) + y) * (nz + 1) + z) as u32;
    let mut tets = Vec::new();
    for x in 0..xs.len() - 1 {
        for y in 0..ny {
            for z in 0..nz {
                let a = index(x, y, z);
                let b = index(x + 1, y, z);
                let c = index(x, y + 1, z);
                let d = index(x + 1, y + 1, z);
                let e = index(x, y, z + 1);
                let f = index(x + 1, y, z + 1);
                let g = index(x, y + 1, z + 1);
                let h = index(x + 1, y + 1, z + 1);
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    (positions, tets)
}

fn terminal_faces(problem: &SpatialElectrostatics) -> [Vec<usize>; 2] {
    let mut terminals = [Vec::new(), Vec::new()];
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [1.0; 3] {
            terminals[0].push(face);
        }
        if xs == [0.0; 3] {
            terminals[1].push(face);
        }
    }
    terminals
}

fn end_voltages(
    problem: &SpatialElectrostatics,
    left: f64,
    right: f64,
) -> SpatialElectrostaticBoundary {
    let terminals = terminal_faces(problem);
    SpatialElectrostaticBoundary {
        voltage: terminals
            .iter()
            .zip([right, left])
            .flat_map(|(faces, value)| {
                faces.iter().map(move |&face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: [value; 3],
                })
            })
            .collect(),
        displacement: vec![],
    }
}

#[test]
fn layered_dielectrics_recover_series_capacitance_and_interface_displacement() {
    let (positions, tets) = cube(&[0.0, 0.25, 1.0], 1, 1);
    let sources = [dielectric(2e-11), dielectric(5e-11)];
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        (0..12).map(|cell| sources[cell / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals = terminal_faces(&problem);
    let capacitor = problem
        .capacitance(&terminals[0], &terminals[1], controls(), &mut || false)
        .unwrap();
    let expected = 1.0 / (0.25 / 2e-11 + 0.75 / 5e-11);
    near(capacitor.capacitance_f(), expected);
    let field = capacitor.unit_voltage_field();
    near(field.energy_j(), 0.5 * expected);
    for cell in 0..12 {
        let epsilon = [2e-11, 5e-11][cell / 6];
        near_vector(
            field.electric_field_v_m()[cell],
            [-expected / epsilon, 0.0, 0.0],
        );
        near_vector(field.displacement_c_m2()[cell], [-expected, 0.0, 0.0]);
        assert_eq!(&problem.materials()[cell], &sources[cell / 6]);
    }
    assert!(field.report().unwrap().converged_euclidean());
    assert!(field.charge_balance_residual_c().abs() < 1e-22);
    assert!(field.energy_balance_residual_j().abs() < 1e-22);
    assert!(field.hodge_energy_residual_j().abs() < 1e-22);

    let direct = problem
        .solve(
            &[0.0; 12],
            &end_voltages(&problem, 2.0, 8.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    near(direct.energy_j(), 18.0 * expected);
    for (v, position) in problem.positions_m().iter().enumerate() {
        let expected_voltage = 2.0
            + 6.0
                * expected
                * (position[0].min(0.25) / 2e-11 + (position[0] - 0.25).max(0.0) / 5e-11);
        near(direct.voltage_v()[v], expected_voltage);
    }
    let reused = capacitor
        .field_at_voltages([8.0, 2.0], controls(), &mut || false)
        .unwrap();
    assert!(reused.report().is_none());
    near(reused.energy_j(), direct.energy_j());
    for cell in 0..12 {
        near_vector(
            reused.electric_field_v_m()[cell],
            direct.electric_field_v_m()[cell],
        );
    }
    let charged = capacitor
        .field_at_charge(-2.5 * expected, controls(), &mut || false)
        .unwrap();
    near(charged.energy_j(), 0.5 * expected * 2.5_f64.powi(2));
    near_vector(charged.displacement_c_m2()[0], [2.5 * expected, 0.0, 0.0]);
}

#[test]
fn retained_relative_field_survives_a_large_common_voltage_reference() {
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
    let original = problem
        .solve(
            &[0.0; 12],
            &end_voltages(&problem, 0.0, 1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    let bias = 2_f64.powi(45);
    let shifted = problem
        .solve(
            &[0.0; 12],
            &end_voltages(&problem, bias, bias + 1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    // Interior absolute voltages lose fractional bits at this bias. The original
    // relative solve must still own E, D, electrode charge and stored energy.
    assert_eq!(original.relative_voltage_v(), shifted.relative_voltage_v());
    assert_eq!(original.electric_field_v_m(), shifted.electric_field_v_m());
    assert_eq!(original.displacement_c_m2(), shifted.displacement_c_m2());
    assert_eq!(
        original.nodal_boundary_charge_c(),
        shifted.nodal_boundary_charge_c()
    );
    assert_eq!(original.energy_j(), shifted.energy_j());
}

#[test]
fn volume_charge_manufactured_solution_converges_in_three_dimensions() {
    let epsilon = 2e-11;
    let potential = |p: [f64; 3]| p.into_iter().map(|x| x * (1.0 - x)).sum::<f64>();
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let xs: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let (positions, tets) = cube(&xs, n, n);
        let problem = SpatialElectrostatics::new(
            positions,
            tets.clone(),
            vec![dielectric(epsilon); tets.len()],
            limits(),
            &mut || false,
        )
        .unwrap();
        let boundary = SpatialElectrostaticBoundary {
            voltage: problem
                .exterior_faces()
                .map(|face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: problem.mesh().faces[face]
                        .map(|v| potential(problem.positions_m()[v as usize])),
                })
                .collect(),
            displacement: vec![],
        };
        let field = problem
            .solve(
                &vec![6.0 * epsilon; tets.len()],
                &boundary,
                controls(),
                &mut || false,
            )
            .unwrap();
        near(field.volume_charge_c(), 6.0 * epsilon);
        near(field.boundary_charge_c(), -6.0 * epsilon);
        let (mut potential_error, mut field_error) = (0.0, 0.0);
        for (cell, tet) in problem.mesh().tets.iter().enumerate() {
            let centre = std::array::from_fn(|j| {
                tet.iter()
                    .map(|&v| problem.positions_m()[v as usize][j] / 4.0)
                    .sum::<f64>()
            });
            let approximation: f64 = tet
                .iter()
                .map(|&v| field.voltage_v()[v as usize] / 4.0)
                .sum();
            let volume = problem.cell_volume_m3(cell).unwrap();
            potential_error += volume * (approximation - potential(centre)).powi(2);
            let exact_field = centre.map(|x| 2.0 * x - 1.0);
            field_error += volume
                * (0..3)
                    .map(|j| (field.electric_field_v_m()[cell][j] - exact_field[j]).powi(2))
                    .sum::<f64>();
        }
        errors.push((potential_error.sqrt(), field_error.sqrt()));
    }
    for pair in errors.windows(2) {
        let potential_order = (pair[0].0 / pair[1].0).log2();
        let field_order = (pair[0].1 / pair[1].1).log2();
        assert!(
            (1.8..2.2).contains(&potential_order),
            "potential order {potential_order}; {errors:?}"
        );
        assert!(
            (0.8..1.2).contains(&field_order),
            "field order {field_order}; {errors:?}"
        );
    }
}

#[test]
fn pure_displacement_boundary_closes_volume_charge_without_a_gauge_electrode() {
    let epsilon = 2e-11;
    let rho = 6.0 * epsilon;
    let (positions, tets) = cube(&[0.0, 0.5, 1.0], 2, 2);
    let cells = tets.len();
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(epsilon); cells],
        limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = SpatialElectrostaticBoundary {
        voltage: vec![],
        displacement: problem
            .exterior_faces()
            .map(|face| {
                let normal = problem.exterior_normal(face).unwrap();
                let point = problem.positions_m()[problem.mesh().faces[face][0] as usize];
                ElectrostaticDisplacementFace {
                    face,
                    outward_displacement_c_m2: rho / 3.0
                        * (0..3).map(|j| normal[j] * point[j]).sum::<f64>(),
                }
            })
            .collect(),
    };
    let field = problem
        .solve(&vec![rho; cells], &boundary, controls(), &mut || false)
        .unwrap();
    assert_eq!(field.gauge_vertices().len(), 1);
    assert_eq!(field.voltage_v()[field.gauge_vertices()[0]], 0.0);
    near(field.volume_charge_c(), rho);
    near(field.boundary_charge_c(), -rho);
    assert!(field.maximum_free_charge_residual_c() < 1e-20);
    assert!(field.charge_balance_residual_c().abs() < 1e-22);
    assert!(field.energy_j() > 0.0);
    // A mathematical reference must not absorb the excess physical charge.
    assert!(
        problem
            .solve(&vec![1.01 * rho; cells], &boundary, controls(), &mut || {
                false
            })
            .is_err()
    );
}

#[test]
fn disconnected_neumann_components_cannot_cancel_each_others_charge() {
    let (mut positions, mut tets) = cube(&[0.0, 1.0], 1, 1);
    let (second_positions, second_tets) = cube(&[0.0, 1.0], 1, 1);
    let offset = positions.len() as u32;
    positions.extend(
        second_positions
            .into_iter()
            .map(|p| [p[0], p[1] + 2.0, p[2]]),
    );
    tets.extend(second_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(2e-11); 12],
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut rho = vec![2e-11; 6];
    rho.extend([-2e-11; 6]);
    assert!(
        problem
            .solve(
                &rho,
                &SpatialElectrostaticBoundary::default(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn oriented_dielectric_and_geometry_rotate_the_full_field_together() {
    let s = 0.5_f64.sqrt();
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    let principal = [2e-11, 5e-11, 7e-11];
    let properties = std::array::from_fn::<_, 3, _>(|i| {
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[i], principal[i])
    });
    let material = source(&properties);
    let electric = [2.0, -1.0, 0.5];
    let displacement = [8.5e-11, -6.5e-11, 3.5e-11];
    let rotate = |p: [f64; 3]| [0.6 * p[0] - 0.8 * p[2], p[1], 0.8 * p[0] + 0.6 * p[2]];
    let mut energies = Vec::new();
    for rotated in [false, true] {
        let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
        let voltages: Vec<_> = positions
            .iter()
            .map(|p| -(0..3).map(|j| electric[j] * p[j]).sum::<f64>())
            .collect();
        let count = tets.len();
        let law = SpatialDielectricLaw::orthotropic(if rotated { axes.map(rotate) } else { axes })
            .unwrap();
        let problem = SpatialElectrostatics::new_with_laws(
            if rotated {
                positions.into_iter().map(rotate).collect()
            } else {
                positions
            },
            tets,
            vec![material.clone(); count],
            vec![law; count],
            limits(),
            &mut || false,
        )
        .unwrap();
        let boundary = SpatialElectrostaticBoundary {
            voltage: problem
                .exterior_faces()
                .map(|face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: problem.mesh().faces[face].map(|v| voltages[v as usize]),
                })
                .collect(),
            displacement: vec![],
        };
        let field = problem
            .solve(&vec![0.0; count], &boundary, controls(), &mut || false)
            .unwrap();
        for cell in 0..count {
            near_vector(
                field.electric_field_v_m()[cell],
                if rotated { rotate(electric) } else { electric },
            );
            near_vector(
                field.displacement_c_m2()[cell],
                if rotated {
                    rotate(displacement)
                } else {
                    displacement
                },
            );
        }
        near(
            field.energy_j(),
            0.5 * (0..3).map(|j| electric[j] * displacement[j]).sum::<f64>(),
        );
        energies.push(field.energy_j());
    }
    near(energies[0], energies[1]);
}

#[test]
fn geometry_derived_capacitance_drives_native_rc_charge_and_heat() {
    let mut next_voltage = Vec::new();
    for epsilon in [4e-6, 8e-6] {
        let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
        let cells = tets.len();
        let problem = SpatialElectrostatics::new(
            positions,
            tets,
            vec![dielectric(epsilon); cells],
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = terminal_faces(&problem);
        let capacitor = problem
            .capacitance(&faces[0], &faces[1], controls(), &mut || false)
            .unwrap();
        let c = capacitor.capacitance_f();
        near(c, epsilon); // Independent parallel-plate C=epsilon A/d.
        let circuit = assemble_circuit(&CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, capacitor.circuit_branch()),
                (1, 0, Branch::Resistor { ohms: 1000.0 }),
            ],
            transformers: vec![],
        })
        .unwrap();
        let mut state = vec![0.0; circuit.system.state_dim()];
        state[circuit.charge_index[0]] = c;
        let state = circuit.consistent_initial_state(&state, &[]).unwrap();
        let (step, defect) = circuit.step_audited(&state, &[], 1e-4).unwrap();
        let ratio = (1.0 - 1e-4 / (2.0 * 1000.0 * c)) / (1.0 + 1e-4 / (2.0 * 1000.0 * c));
        near(step.x[circuit.charge_index[0]], c * ratio);
        near(step.dissipated, 0.5 * c * (1.0 - ratio * ratio));
        let recovered = capacitor
            .field_at_charge(step.x[circuit.charge_index[0]], controls(), &mut || false)
            .unwrap();
        near(recovered.energy_j(), 0.5 * c * ratio * ratio);
        near_vector(recovered.electric_field_v_m()[0], [-ratio, 0.0, 0.0]);
        near(
            capacitor.unit_voltage_field().energy_j() - recovered.energy_j(),
            step.dissipated,
        );
        assert!(defect < 1e-14);
        next_voltage.push(step.x[circuit.charge_index[0]] / c);
    }
    assert!(next_voltage[1] > next_voltage[0]);
}

#[test]
fn invalid_electrodes_loads_and_cancelled_work_return_no_field() {
    let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
    let sources = vec![dielectric(2e-11); tets.len()];
    assert!(matches!(
        SpatialElectrostatics::new(
            positions.clone(),
            tets.clone(),
            sources.clone(),
            limits(),
            &mut || true
        ),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    let mut small = limits();
    small.max_cells = tets.len() - 1;
    assert!(
        SpatialElectrostatics::new(
            positions.clone(),
            tets.clone(),
            sources.clone(),
            small,
            &mut || false
        )
        .is_err()
    );
    let problem =
        SpatialElectrostatics::new(positions, tets, sources, limits(), &mut || false).unwrap();
    let boundary = end_voltages(&problem, 0.0, 1.0);
    assert!(matches!(
        problem.solve(&[0.0; 12], &boundary, controls(), &mut || true),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    assert!(
        problem
            .solve(&[0.0; 11], &boundary, controls(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .solve(&[f64::NAN; 12], &boundary, controls(), &mut || false)
            .is_err()
    );
    let faces = terminal_faces(&problem);
    assert!(
        problem
            .capacitance(&faces[0], &faces[0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .capacitance(&[], &faces[1], controls(), &mut || false)
            .is_err()
    );
    let mut conflicting = boundary;
    conflicting
        .displacement
        .push(ElectrostaticDisplacementFace {
            face: faces[0][0],
            outward_displacement_c_m2: 0.0,
        });
    assert!(
        problem
            .solve(&[0.0; 12], &conflicting, controls(), &mut || false)
            .is_err()
    );
}
