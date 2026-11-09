//! G1/G3/G4: SI cylinder/coaxial resistance and exact radial current moments.
//! Source values are manufactured; these tests do not qualify real materials.

use super::*;
use fs_em::conduction::axisymmetric::AxisymmetricConduction;
use fs_material::conductor::planar::{
    PRINCIPAL_CONDUCTIVITY_PROPERTIES, PlanarElectricalFrame, resolve_orthotropic_conductor2,
};

fn rectangle(a: f64, b: f64, h: f64, nr: usize, nz: usize, angle: f64) -> TriComplex2 {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    for j in 0..=nz {
        for i in 0..=nr {
            vertices.push([
                a + (b - a) * i as f64 / nr as f64,
                h * j as f64 / nz as f64,
            ]);
        }
    }
    for j in 0..nz {
        for i in 0..nr {
            let v = (j * (nr + 1) + i) as u32;
            let top = v + (nr + 1) as u32;
            faces.extend([[v, v + 1, top + 1], [v, top + 1, top]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "axisymmetric conductor fixture",
        vertices,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

#[test]
fn axial_cylinder_resistance_current_and_heat_use_actual_swept_geometry() {
    for basis in [
        ElectricalPropertyBasis::Conductivity,
        ElectricalPropertyBasis::Resistivity,
    ] {
        for (scale, angle) in [(1.0, core::f64::consts::TAU), (0.2, 0.7)] {
            let (radius, length, sigma) = (0.3 * scale, 0.8 * scale, 2.5);
            let mesh = rectangle(0.0, radius, length, 2, 3, angle);
            let cells = mesh.faces().len();
            let material = conductor(sigma, basis);
            let p = AxisymmetricConduction::new(
                mesh,
                vec![material.clone(); cells],
                limits(),
                &mut || false,
            )
            .unwrap();
            assert_eq!(p.conductors()[0], material);
            let bottom = [0, 1, 2];
            let top = [9, 10, 11];
            let response = p.resistance(&top, &bottom, controls(), &mut || false).unwrap();
            let area = angle * radius * radius / 2.0;
            let expected_r = length / (sigma * area);
            near(response.resistance_ohm(), expected_r);
            assert_eq!(
                response.circuit_branch(),
                Branch::Resistor {
                    ohms: response.resistance_ohm(),
                }
            );
            for voltage in [-2.0, 0.0, 1.5] {
                let f = response
                    .field_at_voltage(voltage, 1e-10, 1e-10, &mut || false)
                    .unwrap();
                assert!(f.report().is_none());
                for (point, &v) in p.mesh().vertices().iter().zip(f.voltage_v()) {
                    near(v, voltage * point[1] / length);
                }
                for (cell, (&e, &j)) in f
                    .electric_field_v_m()
                    .iter()
                    .zip(f.current_density_a_m2())
                    .enumerate()
                {
                    near(e[0], 0.0);
                    near(e[1], -voltage / length);
                    near(j[0], 0.0);
                    near(j[1], -sigma * voltage / length);
                    near(
                        f.cell_power_w()[cell],
                        sigma * (voltage / length).powi(2)
                            * p.mesh().face_measure(cell).unwrap(),
                    );
                }
                near(f.power_w(), voltage * voltage / expected_r);
                near(f.boundary_current_a(), 0.0);
                near(f.power_balance_residual_w(), 0.0);
            }
            for current in [-0.2, 0.0, 0.1] {
                let f = response
                    .field_at_current(current, 1e-10, 1e-10, &mut || false)
                    .unwrap();
                near(f.power_w(), current * current * expected_r);
                near(
                    top.iter()
                        .map(|&v| f.nodal_boundary_current_a()[v].unwrap())
                        .sum(),
                    current,
                );
            }
            let f = response
                .field_at_power(
                    ConductionPowerSource {
                        power_w: 0.02,
                        maximum_voltage_v: 100.0,
                        maximum_current_a: 100.0,
                    },
                    1e-10,
                    1e-10,
                    &mut || false,
                )
                .unwrap();
            near(f.power_w(), 0.02);
        }
    }
}

#[test]
fn density_on_a_disk_uses_exact_radial_p1_moments() {
    let angle = 1.7;
    let radius = 0.6;
    let length = 0.9;
    let sigma = 2.0;
    let density = 0.4;
    let mesh = rectangle(0.0, radius, length, 1, 2, angle);
    let cells = mesh.faces().len();
    let p = AxisymmetricConduction::new(
        mesh,
        vec![conductor(sigma, ElectricalPropertyBasis::Conductivity); cells],
        limits(),
        &mut || false,
    )
    .unwrap();
    let total = density * angle * radius * radius / 2.0;
    for vertices in [[4, 5], [5, 4]] {
        let f = p
            .solve_with_current_traces(
                &[(0, 0.0), (1, 0.0)],
                &[ConductionCurrentTrace {
                    vertices,
                    inward_current_density_a_m2: density,
                }],
                controls(),
                &mut || false,
            )
            .unwrap();
        near(f.nodal_applied_current_a()[4].unwrap(), total / 3.0);
        near(f.nodal_applied_current_a()[5].unwrap(), 2.0 * total / 3.0);
        for (point, &v) in p.mesh().vertices().iter().zip(f.voltage_v()) {
            near(v, density * point[1] / sigma);
        }
        for &j in f.current_density_a_m2() {
            near(j[0], 0.0);
            near(j[1], -density);
        }
        near(f.power_w(), total * density * length / sigma);
        near(f.boundary_current_a(), 0.0);
        near(f.power_balance_residual_w(), 0.0);
    }
}

#[test]
fn coaxial_series_resistance_and_material_heat_converge_to_logarithmic_oracle() {
    let (angle, length) = (1.3, 0.8);
    let (sigma_a, sigma_b) = (0.25, 0.75);
    let part_r = [2.0_f64.ln() / sigma_a, 1.5_f64.ln() / sigma_b]
        .map(|value| value / (angle * length));
    let exact_r = part_r[0] + part_r[1];
    let mut previous_error = f64::INFINITY;
    let cfg = ConductionControls {
        linear: fs_solver::CheckedCgConfig {
            tolerance: 1e-13,
            ..config()
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    };
    for n in [8, 16, 32] {
        let mesh = rectangle(0.5, 1.5, length, n, n, angle);
        let side: Vec<_> = mesh
            .faces()
            .iter()
            .map(|face| {
                let r = face
                    .iter()
                    .map(|&v| mesh.vertices()[v as usize][0])
                    .sum::<f64>()
                    / 3.0;
                usize::from(r >= 1.0)
            })
            .collect();
        let materials = [sigma_a, sigma_b]
            .map(|sigma| conductor(sigma, ElectricalPropertyBasis::Conductivity));
        let p = AxisymmetricConduction::new(
            mesh,
            side.iter().map(|&i| materials[i].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let inner: Vec<_> = (0..=n).map(|j| j * (n + 1)).collect();
        let outer: Vec<_> = inner.iter().map(|v| v + n).collect();
        let response = p.resistance(&inner, &outer, cfg, &mut || false).unwrap();
        let error = (response.resistance_ohm() / exact_r - 1.0).abs();
        assert!(
            error < 0.45 * previous_error,
            "n={n}: {error}, previous={previous_error}"
        );
        previous_error = error;
        let f = response
            .field_at_voltage(1.0, 1e-9, 1e-9, &mut || false)
            .unwrap();
        let current = 1.0 / response.resistance_ohm();
        for (&[r, _], &voltage) in p.mesh().vertices().iter().zip(f.voltage_v()) {
            let remaining = if r < 1.0 {
                (1.0 / r).ln() / sigma_a + 1.5_f64.ln() / sigma_b
            } else {
                (1.5 / r).ln() / sigma_b
            } / (angle * length);
            assert!((voltage - remaining / exact_r).abs() < 2.0 / (n * n) as f64);
        }
        let mut heat = [0.0; 2];
        for (&which, &power) in side.iter().zip(f.cell_power_w()) {
            heat[which] += power;
        }
        for i in 0..2 {
            assert!(
                (heat[i] / (current * current * part_r[i]) - 1.0).abs()
                    < 2.0 / (n * n) as f64
            );
        }
        near(f.power_w(), current);
        near(f.power_balance_residual_w(), 0.0);
    }
    assert!(previous_error < 5e-4);
}

#[test]
fn rotated_meridional_conductivity_keeps_full_tensor_and_insulating_radial_walls() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut claims = ClaimSet::new();
    for (name, value) in PRINCIPAL_CONDUCTIVITY_PROPERTIES.into_iter().zip([2.0, 8.0]) {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, basis.dimensions()),
                value: PropertyValue::Scalar {
                    value,
                    dims: basis.dimensions(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic meridional tensor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured orthotropic".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let frame = PlanarElectricalFrame::new(
        card.content_hash(),
        conductor(1.0, basis).material().card_identity(),
        core::f64::consts::FRAC_PI_4,
    )
    .unwrap();
    let material = resolve_orthotropic_conductor2(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        frame,
    )
    .unwrap();
    let mesh = rectangle(0.5, 1.5, 1.0, 2, 2, 0.7);
    let cells = mesh.faces().len();
    let p = AxisymmetricConduction::new(mesh, vec![material; cells], limits(), &mut || false)
        .unwrap();
    // Sigma=[[5,-3],[-3,5]], V=0.6r+z gives J=(0,-3.2).
    // Both radial walls insulate, including the tensor cross term; div(J)=0.
    let boundary: Vec<_> = p
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, point)| point[1] == 0.0 || point[1] == 1.0)
        .map(|(i, point)| (i, 0.6 * point[0] + point[1]))
        .collect();
    let f = p.solve(&boundary, controls(), &mut || false).unwrap();
    for (&[r, z], &v) in p.mesh().vertices().iter().zip(f.voltage_v()) {
        near(v, 0.6 * r + z);
    }
    for &j in f.current_density_a_m2() {
        near(j[0], 0.0);
        near(j[1], -3.2);
    }
    near(f.power_w(), 3.2 * 0.7);
    near(f.power_balance_residual_w(), 0.0);
}

#[test]
fn axis_electrodes_wrong_metrics_budgets_and_cancellation_refuse() {
    let source = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let mesh = rectangle(0.0, 1.0, 1.0, 2, 2, 1.0);
    let materials = vec![source; mesh.faces().len()];
    assert!(PlanarConduction::new(
        mesh.clone(), materials.clone(), limits(), &mut || false,
    ).is_err());
    let planar = TriComplex2::from_indexed_triangles(
        "wrong conductor metric",
        mesh.vertices().to_vec(),
        mesh.faces().to_vec(),
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    assert!(AxisymmetricConduction::new(
        planar, materials.clone(), limits(), &mut || false,
    ).is_err());
    let mut too_small = limits();
    too_small.max_triplets = 1;
    assert!(AxisymmetricConduction::new(
        mesh.clone(), materials.clone(), too_small, &mut || false,
    ).is_err());
    assert!(matches!(
        AxisymmetricConduction::new(
            mesh.clone(), materials.clone(), limits(), &mut || true,
        ),
        Err(ConductionError::Cancelled)
    ));
    let p = AxisymmetricConduction::new(mesh, materials, limits(), &mut || false).unwrap();
    assert!(p.resistance(&[0, 3, 6], &[2, 5, 8], controls(), &mut || false).is_err());
    assert!(p.solve(&[(0, 0.0), (3, 0.0), (6, 0.0)], controls(), &mut || false).is_err());
    for density in [0.0, 1.0] {
        assert!(p.solve_with_current_traces(
            &[(2, 0.0), (5, 0.0), (8, 0.0)],
            &[ConductionCurrentTrace {
                vertices: [0, 3],
                inward_current_density_a_m2: density,
            }],
            controls(),
            &mut || false,
        ).is_err());
    }
    let boundary = [
        (0, 0.0), (1, 0.0), (2, 0.0),
        (6, 1.0), (7, 1.0), (8, 1.0),
    ];
    let expected = p.solve(&boundary, controls(), &mut || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        p.solve(&boundary, controls(), &mut || {
            polls += 1;
            polls == 20
        }),
        Err(ConductionError::Cancelled)
    ));
    let retry = p.solve(&boundary, controls(), &mut || false).unwrap();
    assert_eq!(retry.voltage_v(), expected.voltage_v());
    assert_eq!(retry.cell_power_w(), expected.cell_power_w());
}
