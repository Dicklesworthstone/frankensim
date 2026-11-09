//! G1/G3/G4: original material interfaces on actual surfaces of revolution.
//! Manufactured SI laws check voltage jumps, radial surface integration,
//! material heat ownership and source/geometry refusal.

use super::*;
use fs_em::conduction::{
    ConductionContact,
    axisymmetric::{AxisymmetricConductanceControls, AxisymmetricConduction},
};

fn axial_mesh(inner: f64, outer: f64, lengths: [f64; 2], angle: f64) -> TriComplex2 {
    let [a, b] = lengths;
    TriComplex2::from_indexed_triangles(
        "axisymmetric axial material interface",
        vec![
            [inner, 0.0],
            [outer, 0.0],
            [outer, a],
            [inner, a],
            [inner, a],
            [outer, a],
            [outer, a + b],
            [inner, a + b],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn specimen(
    inner: f64,
    outer: f64,
    lengths: [f64; 2],
    angle: f64,
    pressure: f64,
) -> AxisymmetricConduction {
    let lower = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let upper = conductor(4.0, ElectricalPropertyBasis::Resistivity);
    let source = contacts::interface(&lower, &upper, pressure, 1.0);
    AxisymmetricConduction::new_with_contacts(
        axial_mesh(inner, outer, lengths, angle),
        vec![lower.clone(), lower, upper.clone(), upper],
        vec![ConductionContact {
            first_vertices: [3, 2],
            second_vertices: [4, 5],
            source,
        }],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn axial_contact_resistance_heat_and_current_use_actual_disk_and_annular_area() {
    for (inner, outer, angle, pressure) in [
        (0.0, 0.6, core::f64::consts::TAU, 1e5),
        (0.0, 0.6, 0.7, 2e5),
        (0.2, 0.8, 1.3, 1.5e5),
    ] {
        let lengths = [0.4, 0.9];
        let p = specimen(inner, outer, lengths, angle, pressure);
        let area = angle * (outer * outer - inner * inner) / 2.0;
        let specific = 0.5 - 0.25 * (pressure - 1e5) / 1e5;
        let exact_r = (lengths[0] / 2.0 + specific + lengths[1] / 4.0) / area;
        let resistance = p
            .resistance(&[6, 7], &[0, 1], controls(), &mut || false)
            .unwrap();
        near(resistance.resistance_ohm(), exact_r);
        assert_eq!(
            p.contacts()[0].source.resolved().query_point(),
            &[("normal_pressure".into(), pressure)],
        );
        assert_eq!(p.contact_liquid_mass_fraction_axes(), &[[None, None]]);
        for voltage in [-2.0, 0.0, 3.0] {
            let field = resistance
                .field_at_voltage(voltage, 1e-10, 1e-10, &mut || false)
                .unwrap();
            let current = voltage / exact_r;
            let density = current / area;
            let contact = &field.contact_response()[0];
            near(contact.current_a, -current);
            for jump in contact.voltage_jump_v {
                near(jump, -density * specific);
            }
            for j in contact.current_density_a_m2 {
                near(j, -density);
            }
            near(contact.power_w, current * current * specific / area);
            for cell in 0..4 {
                let sigma = if cell < 2 { 2.0 } else { 4.0 };
                near(field.current_density_a_m2()[cell][0], 0.0);
                near(field.current_density_a_m2()[cell][1], -density);
                near(
                    field.cell_power_w()[cell],
                    density * density / sigma * p.mesh().face_measure(cell).unwrap(),
                );
            }
            near(field.power_w(), current * voltage);
            near(
                field.power_w(),
                field.cell_power_w().iter().sum::<f64>() + contact.power_w,
            );
            near(field.boundary_current_a(), 0.0);
            near(field.power_balance_residual_w(), 0.0);
        }
        // The same physical contacts participate in the reusable terminal map.
        let matrix = p
            .conductance_matrix(
                &[&[0, 1], &[6, 7]],
                AxisymmetricConductanceControls {
                    maximum_terminals: 2,
                    maximum_field_solves: 1,
                    maximum_response_values: 1000,
                    maximum_dense_work: 1000,
                    reference_voltage_v: 0.5,
                    coefficient_tolerance_s: 1e-9,
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        let field = matrix
            .field_at_currents(&[-0.2, 0.2], 1e-9, 1e-9, &mut || false)
            .unwrap();
        near(field.contact_response()[0].current_a, -0.2);
        near(field.contact_response()[0].power_w, 0.04 * specific / area);
        near(field.power_w(), 0.04 * exact_r);

        // Source currents and contact jumps retain drops below a common bias.
        for bias in [0.0, 2.0_f64.powi(56), -2.0_f64.powi(56)] {
            let field = p
                .solve_with_current_traces(
                    &[(0, bias), (1, bias)],
                    &[ConductionCurrentTrace {
                        vertices: [7, 6],
                        inward_current_density_a_m2: 0.3,
                    }],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(field.contact_response()[0].current_a, -0.3 * area);
            near(field.contact_response()[0].power_w, 0.09 * specific * area);
            near(field.power_w(), 0.09 * area * area * exact_r);
            near(field.power_balance_residual_w(), 0.0);
        }
    }
}

#[test]
fn varying_contact_jump_uses_exact_radial_gauss_current_and_heat() {
    let (inner, outer, angle) = (0.0, 0.7, 1.4);
    let p = specimen(inner, outer, [0.4, 0.9], angle, 1e5);
    let field = p
        .solve(
            &[(0, 0.0), (1, 0.0), (6, 2.0), (7, 1.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let response = &field.contact_response()[0];
    let [d0, d1] = response.voltage_jump_v;
    assert!((d0 - d1).abs() > 1e-6, "nonuniform jump must be exercised");
    let mut current = 0.0;
    let mut heat = 0.0;
    let offset = 1.0 / 3.0_f64.sqrt();
    for s in [(1.0 - offset) / 2.0, (1.0 + offset) / 2.0] {
        let radius = (1.0 - s) * inner + s * outer;
        let jump = (1.0 - s) * d0 + s * d1;
        let surface_weight = 0.5 * angle * (outer - inner) * radius;
        current += surface_weight * jump / 0.5;
        heat += surface_weight * jump * jump / 0.5;
    }
    near(response.current_a, current);
    near(response.power_w, heat);
    near(
        field.power_w(),
        field.cell_power_w().iter().sum::<f64>() + heat,
    );
    near(field.power_balance_residual_w(), 0.0);

    // Reverse both paired traces: the physical field is unchanged and the
    // original endpoint response reverses its reporting order.
    let mut reversed = p.contacts()[0].clone();
    reversed.first_vertices.reverse();
    reversed.second_vertices.reverse();
    let reversed = AxisymmetricConduction::new_with_contacts(
        p.mesh().clone(),
        p.conductors().to_vec(),
        vec![reversed],
        limits(),
        &mut || false,
    )
    .unwrap();
    let other = reversed
        .solve(
            &[(0, 0.0), (1, 0.0), (6, 2.0), (7, 1.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in field.voltage_v().iter().zip(other.voltage_v()) {
        near(a, b);
    }
    near(other.contact_response()[0].voltage_jump_v[0], d1);
    near(other.contact_response()[0].voltage_jump_v[1], d0);
    near(other.contact_response()[0].current_a, current);
    near(other.contact_response()[0].power_w, heat);
}

#[test]
fn coaxial_contact_resistance_and_separate_material_heat_converge() {
    let (angle, length, interface_radius) = (1.3, 0.8, 1.0);
    let sigma = [0.25, 0.75];
    let exact_parts = [
        2.0_f64.ln() / (sigma[0] * angle * length),
        0.5 / (angle * interface_radius * length),
        1.5_f64.ln() / (sigma[1] * angle * length),
    ];
    let exact_r = exact_parts.iter().sum::<f64>();
    let controls = ConductionControls {
        linear: fs_solver::CheckedCgConfig {
            tolerance: 1e-13,
            ..config()
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    };
    let mut previous = f64::INFINITY;
    for nr in [8, 16, 32] {
        let nz = 4;
        let first = axisymmetric::rectangle(0.5, interface_radius, length, nr, nz, angle);
        let second = axisymmetric::rectangle(interface_radius, 1.5, length, nr, nz, angle);
        let offset = first.vertices().len();
        let cells = first.faces().len();
        let mut vertices = first.vertices().to_vec();
        vertices.extend_from_slice(second.vertices());
        let mut faces = first.faces().to_vec();
        faces.extend(
            second
                .faces()
                .iter()
                .map(|face| face.map(|v| v + offset as u32)),
        );
        let mesh = TriComplex2::from_indexed_triangles(
            "coaxial finite-resistance material interface",
            vertices,
            faces,
            Metric2::axisymmetric(angle).unwrap(),
        )
        .unwrap();
        let materials = sigma.map(|s| conductor(s, ElectricalPropertyBasis::Conductivity));
        let source = contacts::interface(&materials[0], &materials[1], 1e5, 1.0);
        let interfaces = (0..nz)
            .map(|j| ConductionContact {
                first_vertices: [j * (nr + 1) + nr, (j + 1) * (nr + 1) + nr],
                second_vertices: [offset + j * (nr + 1), offset + (j + 1) * (nr + 1)],
                source: source.clone(),
            })
            .collect();
        let p = AxisymmetricConduction::new_with_contacts(
            mesh,
            (0..2 * cells)
                .map(|cell| materials[usize::from(cell >= cells)].clone())
                .collect(),
            interfaces,
            limits(),
            &mut || false,
        )
        .unwrap();
        let inner: Vec<_> = (0..=nz).map(|j| j * (nr + 1)).collect();
        let outer: Vec<_> = (0..=nz).map(|j| offset + j * (nr + 1) + nr).collect();
        let response = p
            .resistance(&inner, &outer, controls, &mut || false)
            .unwrap();
        let error = (response.resistance_ohm() / exact_r - 1.0).abs();
        assert!(error < 0.6 * previous, "nr={nr}: {error} vs {previous}");
        previous = error;
        let f = response
            .field_at_voltage(1.0, 1e-9, 1e-9, &mut || false)
            .unwrap();
        let current = 1.0 / response.resistance_ohm();
        near(
            f.contact_response().iter().map(|c| c.current_a).sum(),
            current,
        );
        let heat = [
            f.cell_power_w()[..cells].iter().sum::<f64>(),
            f.contact_response().iter().map(|c| c.power_w).sum(),
            f.cell_power_w()[cells..].iter().sum(),
        ];
        for (actual, r) in heat.into_iter().zip(exact_parts) {
            assert!((actual / (current * current * r) - 1.0).abs() < 4.0 / (nr * nr) as f64);
        }
        near(f.power_w(), current);
        near(f.power_balance_residual_w(), 0.0);
    }
    assert!(previous < 5e-4);
}

#[test]
fn contact_geometry_sources_budgets_and_cancellation_refuse_without_mutation() {
    let p = specimen(0.0, 0.6, [0.4, 0.9], 1.7, 1e5);
    let contact = p.contacts()[0].clone();
    let build = |contacts, bounds| {
        AxisymmetricConduction::new_with_contacts(
            p.mesh().clone(),
            p.conductors().to_vec(),
            contacts,
            bounds,
            &mut || false,
        )
    };
    assert!(build(vec![contact.clone(), contact.clone()], limits()).is_err());
    let mut crossed = contact.clone();
    crossed.second_vertices.reverse();
    assert!(build(vec![crossed], limits()).is_err());
    let mut misplaced = contact.clone();
    misplaced.second_vertices = [7, 6];
    assert!(build(vec![misplaced], limits()).is_err());
    let mut perfect = contact.clone();
    perfect.source = contacts::interface(&p.conductors()[0], &p.conductors()[2], 1e5, 0.0);
    assert!(build(vec![perfect], limits()).is_err());
    let mut bounds = limits();
    bounds.max_triplets = p.mesh().faces().len() * 9 + 15;
    assert!(build(vec![contact.clone()], bounds).is_err());
    assert!(
        PlanarConduction::new_with_contacts(
            p.mesh().clone(),
            p.conductors().to_vec(),
            vec![contact.clone()],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let planar = TriComplex2::from_indexed_triangles(
        "wrong contact metric",
        p.mesh().vertices().to_vec(),
        p.mesh().faces().to_vec(),
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    assert!(
        AxisymmetricConduction::new_with_contacts(
            planar,
            p.conductors().to_vec(),
            vec![contact.clone()],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        AxisymmetricConduction::new_with_phase_contacts(
            p.mesh().clone(),
            p.conductors().to_vec(),
            vec![contact.clone()],
            vec![[Some("unowned_phase".into()), None]],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        AxisymmetricConduction::new_with_contacts(
            p.mesh().clone(),
            p.conductors().to_vec(),
            vec![contact],
            limits(),
            &mut || true,
        ),
        Err(ConductionError::Cancelled)
    ));
    // An internal material interface cannot also be an external electrode.
    assert!(
        p.solve(
            &[(0, 0.0), (1, 0.0), (3, 1.0), (2, 1.0)],
            controls(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        p.solve_with_current_traces(
            &[(0, 0.0), (1, 0.0)],
            &[ConductionCurrentTrace {
                vertices: [4, 5],
                inward_current_density_a_m2: 0.1
            }],
            controls(),
            &mut || false,
        )
        .is_err()
    );
    let boundary = [(0, 0.0), (1, 0.0), (6, 2.0), (7, 1.0)];
    let mut polls = 0;
    let expected = p
        .solve(&boundary, controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut attempts = 0;
    assert!(matches!(
        p.solve(&boundary, controls(), &mut || {
            attempts += 1;
            attempts == polls - 1
        }),
        Err(ConductionError::Cancelled)
    ));
    let retry = p.solve(&boundary, controls(), &mut || false).unwrap();
    assert_eq!(retry.voltage_v(), expected.voltage_v());
    assert_eq!(retry.cell_power_w(), expected.cell_power_w());
    assert_eq!(retry.contact_response(), expected.contact_response());
}
