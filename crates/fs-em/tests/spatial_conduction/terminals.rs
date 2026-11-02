//! G1/G3/G4: total-current electrodes against independently derived series,
//! parallel and oriented-material solutions, including actual local heating.
use super::*;

pub(super) fn terminal_faces<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    right: f64,
) -> [Vec<usize>; 2] {
    let mut faces = [Vec::new(), Vec::new()];
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [right; 3] {
            faces[0].push(face);
        }
        if xs == [0.0; 3] {
            faces[1].push(face);
        }
    }
    faces
}

#[test]
fn current_driven_heterogeneous_resistance_reuses_native_heat_at_each_amplitude() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let rho = [2.0, 5.0];
        let sources = rho.map(|r: f64| {
            ResolvedConductorState::try_from_material_with_basis(
                material(
                    basis,
                    &[if basis == ElectricalPropertyBasis::Resistivity {
                        r
                    } else {
                        r.recip()
                    }],
                ),
                basis,
            )
            .unwrap()
        });
        let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
        let problem = SpatialConduction::new(
            positions,
            tets,
            (0..12).map(|i| sources[i / 6].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = terminal_faces(&problem, 2.0);
        for reference in [1.0, 0.3] {
            let resistance = problem
                .resistance_at_reference_voltage(
                    &faces[0],
                    &faces[1],
                    reference,
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(resistance.resistance_ohm(), 7.0); // sum(rho_i L_i / A)
            assert_eq!(resistance.reference_voltage_v(), reference);
            assert!(resistance.reference_field().linear.is_some());
            near(
                resistance.reference_field().power_w,
                reference * reference / 7.0,
            );
            for current in [-1.25, 0.0, 2.5] {
                let field = resistance
                    .field_at_current(
                        current,
                        controls().current_tolerance_a,
                        controls().power_tolerance_w,
                        &mut || false,
                    )
                    .unwrap();
                assert!(field.linear.is_none());
                for (v, p) in problem.positions_m().iter().enumerate() {
                    near(
                        field.voltage_v[v],
                        current * (2.0 * p[0].min(1.0) + 5.0 * (p[0] - 1.0).max(0.0)),
                    );
                }
                for cell in 0..12 {
                    near(field.current_density_a_m2[cell][0], -current);
                    near(
                        field.cell_power_w[cell],
                        current * current * rho[cell / 6] / 6.0,
                    );
                }
                near(field.power_w, 7.0 * current * current);
                near(field.supplied_power_w, field.power_w);
                let direct = problem
                    .solve_current_with_reference_voltage(
                        &faces[0],
                        &faces[1],
                        current,
                        reference,
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                assert_eq!(field.voltage_v, direct.voltage_v);
                assert_eq!(field.cell_power_w, direct.cell_power_w);
            }
            let reversed = resistance
                .field_at_voltage(-3.5, 1e-10, 1e-9, &mut || false)
                .unwrap();
            near(reversed.power_w, 1.75);
            near(reversed.current_density_a_m2[0][0], 0.5);
        }
    }
}

#[test]
fn equipotential_parallel_paths_determine_nonuniform_current_injection() {
    let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
    let (mut second, second_tets) = slab(&[0.0, 0.5, 1.0]);
    for p in &mut second {
        p[1] += 2.0;
    }
    positions.extend(second);
    tets.extend(second_tets.into_iter().map(|t| t.map(|v| v + 12)));
    let sources = [1.0, 4.0].map(|sigma| {
        ResolvedConductorState::try_from_material(material(
            ElectricalPropertyBasis::Conductivity,
            &[sigma],
        ))
        .unwrap()
    });
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..24).map(|i| sources[i / 12].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut faces = terminal_faces(&problem, 1.0);
    let resistance = problem
        .resistance(&faces[0], &faces[1], controls(), &mut || false)
        .unwrap();
    near(resistance.resistance_ohm(), 1.0 / 5.0);
    let field = problem
        .solve_current(&faces[0], &faces[1], 5.0, controls(), &mut || false)
        .unwrap();
    // Equal geometry and equipotential terminals: currents split as sigma,
    // NOT equally by area. Each branch has V=1 and I=(1,4) A.
    for cell in 0..24 {
        near(field.current_density_a_m2[cell][0], -[1.0, 4.0][cell / 12]);
    }
    near(field.cell_power_w[..12].iter().sum(), 1.0);
    near(field.cell_power_w[12..].iter().sum(), 4.0);
    for (v, p) in problem.positions_m().iter().enumerate() {
        near(field.voltage_v[v], p[0]);
    }
    near(field.power_w, 5.0);
    for terminal in &mut faces {
        terminal.reverse();
    }
    let repeated = problem
        .solve_current(&faces[0], &faces[1], 5.0, controls(), &mut || false)
        .unwrap();
    assert_eq!(field.voltage_v, repeated.voltage_v);
    assert_eq!(field.cell_power_w, repeated.cell_power_w);
    // Pick only the positive face of one branch and the return of the other:
    // both are valid boundary faces, but there is no conducting path.
    let a: Vec<_> = faces[0]
        .iter()
        .copied()
        .filter(|&f| problem.mesh().faces[f][0] < 12)
        .collect();
    let b: Vec<_> = faces[1]
        .iter()
        .copied()
        .filter(|&f| problem.mesh().faces[f][0] >= 12)
        .collect();
    for current in [0.0, 1.0] {
        assert!(
            problem
                .solve_current(&a, &b, current, controls(), &mut || false)
                .is_err()
        );
    }
}

#[test]
fn oriented_conductor_terminal_resistance_scales_with_actual_geometry() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    let source = ResolvedOrthotropicConductor3::try_from_material(
        material(ElectricalPropertyBasis::Conductivity, &[2.0, 8.0, 3.0]),
        ElectricalPropertyBasis::Conductivity,
        std::array::from_fn(|i| {
            PropertyKey::new(
                format!("principal_{i}"),
                ElectricalPropertyBasis::Conductivity.dimensions(),
            )
        }),
        axes,
    )
    .unwrap();
    for scale in [1.0, 0.25] {
        let (mut positions, tets) = slab(&[0.0, 0.5, 1.0]);
        for p in &mut positions {
            *p = [
                s * (p[0] - p[1]) * scale,
                s * (p[0] + p[1]) * scale,
                p[2] * scale,
            ];
        }
        let problem = SpatialConduction::new(
            positions,
            tets,
            vec![source.clone(); 12],
            limits(),
            &mut || false,
        )
        .unwrap();
        let positive: Vec<_> = problem
            .exterior_faces()
            .filter(|&f| {
                problem.mesh().faces[f]
                    .iter()
                    .all(|&v| (8..12).contains(&v))
            })
            .collect();
        let negative: Vec<_> = problem
            .exterior_faces()
            .filter(|&f| problem.mesh().faces[f].iter().all(|&v| v < 4))
            .collect();
        let resistance = problem
            .resistance(&positive, &negative, controls(), &mut || false)
            .unwrap();
        near(resistance.resistance_ohm(), 1.0 / (2.0 * scale));
        let field = resistance
            .field_at_current(0.75, 1e-10, 1e-9, &mut || false)
            .unwrap();
        for j in &field.current_density_a_m2 {
            near(j[0], -0.75 * s / scale.powi(2));
            near(j[1], -0.75 * s / scale.powi(2));
            near(j[2], 0.0);
        }
        near(field.power_w, 0.75f64.powi(2) / (2.0 * scale));
    }
}

#[test]
fn invalid_electrodes_and_cancelled_recovery_preserve_prepared_response() {
    let source = ResolvedConductorState::try_from_material(material(
        ElectricalPropertyBasis::Conductivity,
        &[2.0],
    ))
    .unwrap();
    let (positions, tets) = slab(&[0.0, 0.5, 1.0]);
    let problem =
        SpatialConduction::new(positions, tets, vec![source; 12], limits(), &mut || false).unwrap();
    let faces = terminal_faces(&problem, 1.0);
    let interior = (0..problem.mesh().faces.len())
        .find(|f| !problem.exterior_faces().any(|x| x == *f))
        .unwrap();
    for (positive, negative) in [
        (vec![], faces[1].clone()),
        (vec![faces[0][0]; 2], faces[1].clone()),
        (faces[0].clone(), faces[0].clone()),
        (vec![interior], faces[1].clone()),
        (vec![usize::MAX], faces[1].clone()),
    ] {
        assert!(
            problem
                .resistance(&positive, &negative, controls(), &mut || false)
                .is_err()
        );
    }
    for bad in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        assert!(
            problem
                .resistance_at_reference_voltage(&faces[0], &faces[1], bad, controls(), &mut || {
                    false
                },)
                .is_err()
        );
    }
    let response = problem
        .resistance(&faces[0], &faces[1], controls(), &mut || false)
        .unwrap();
    let before = response
        .field_at_current(1.0, 1e-10, 1e-9, &mut || false)
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        response.field_at_current(1.0, 1e-10, 1e-9, &mut || {
            calls += 1;
            calls > 30
        }),
        Err(SpatialConductionError::Cancelled)
    ));
    for current in [f64::NAN, f64::INFINITY] {
        assert!(
            response
                .field_at_current(current, 1e-10, 1e-9, &mut || false)
                .is_err()
        );
    }
    let after = response
        .field_at_current(1.0, 1e-10, 1e-9, &mut || false)
        .unwrap();
    assert_eq!(before.voltage_v, after.voltage_v);
    assert_eq!(before.cell_power_w, after.cell_power_w);
}
