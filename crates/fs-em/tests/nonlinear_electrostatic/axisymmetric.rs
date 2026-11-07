//! G1/G3/G4: nonlinear radial storage, charge, tangent and virtual work.
use super::*;
use fs_em::electrostatic::{
    axisymmetric::AxisymmetricElectrostatics,
    nonlinear::{NonlinearElectrostaticError, axisymmetric::AxisymmetricNonlinearElectrostatics},
};

const RADIUS: f64 = 2e-2;
const ANGLE: f64 = core::f64::consts::TAU;

fn cylinder(n: usize, angle: f64, radial_scale: f64, axial_scale: f64) -> TriComplex2 {
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            vertices.push([
                RADIUS * radial_scale * i as f64 / n as f64,
                GAP * axial_scale * j as f64 / n as f64,
            ]);
        }
    }
    for j in 0..n {
        for i in 0..n {
            let a = (j * (n + 1) + i) as u32;
            let b = a + (n + 1) as u32;
            triangles.extend([[a, a + 1, b + 1], [a, b + 1, b]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "nonlinear revolved dielectric",
        vertices,
        triangles,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn problem(mesh: TriComplex2, material: DielectricMaterial) -> AxisymmetricNonlinearElectrostatics {
    let cells = mesh.faces().len();
    AxisymmetricNonlinearElectrostatics::new(mesh, vec![material; cells], limits(), &mut || false)
        .unwrap()
}

#[test]
fn nonlinear_revolved_slab_preserves_material_storage_total_charge_and_local_tangent() {
    let law = nonlinear();
    let response = law.evaluate(2.0 * E0).unwrap();
    for angle in [ANGLE, ANGLE / 7.0] {
        let p = problem(
            cylinder(2, angle, 1.0, 1.0),
            DielectricMaterial::Reversible(law.clone()),
        );
        let electrodes = p
            .prepare_electrodes(&[&[0, 1, 2], &[6, 7, 8]], &mut || false)
            .unwrap();
        let area = 0.5 * angle * RADIUS * RADIUS;
        let expected_q = area * response.displacement_c_m2;
        let voltage = electrodes
            .solve(
                &[0.0; 8],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(2.0 * V0),
                ],
                &[0.0; 9],
                controls(),
                &mut || false,
            )
            .unwrap();
        near(voltage.terminal_charge_c()[1], expected_q, expected_q);
        near(voltage.terminal_charge_c()[0], -expected_q, expected_q);
        near(
            voltage.energy_j(),
            area * GAP * response.energy_density_j_m3,
            U0,
        );
        near(
            voltage.coenergy_j(),
            area * GAP * response.coenergy_density_j_m3,
            U0,
        );
        assert!(voltage.coenergy_j() > voltage.energy_j());
        for (&[er, ez], &[dr, dz]) in voltage
            .electric_field_v_m()
            .iter()
            .zip(voltage.displacement_c_m2())
        {
            near(er, 0.0, E0);
            near(ez, -2.0 * E0, E0);
            near(dr, 0.0, D0);
            near(dz, -response.displacement_c_m2, D0);
        }
        let differential = voltage
            .differential_capacitance(differential_controls(), &mut || false)
            .unwrap();
        let capacitance = area * response.differential_permittivity_f_m / GAP;
        for (&actual, expected) in differential.matrix_f().iter().zip([
            capacitance,
            -capacitance,
            -capacitance,
            capacitance,
        ]) {
            near(actual, expected, capacitance);
        }
        assert!(core::ptr::eq(differential.field(), &voltage));
        // A zero guess exercises nonlinear charge control rather than Q/C
        // reconstruction from a fitted secant or a unit-voltage response.
        let charge = electrodes
            .solve(
                &[0.0; 8],
                &[
                    ElectrodeControl::Voltage(1e12),
                    ElectrodeControl::Charge(expected_q),
                ],
                &[0.0; 9],
                controls(),
                &mut || false,
            )
            .unwrap();
        near(charge.terminal_voltage_v()[1], 2.0 * V0, V0);
        near(charge.energy_j(), voltage.energy_j(), U0);
        assert!(charge.maximum_floating_charge_residual_c() < controls().charge_tolerance_c);
    }
}

#[test]
fn nonlinear_coaxial_field_and_differential_capacitance_converge_to_gauss_law() {
    // Entire annulus lies in the original affine segment D=a+s*E.
    // Gauss: D_r=k/r; V(r)=[k*ln(b/r)-a*(b-r)]/s.
    // This retains the nonlinear offset a and differs from a secant epsilon.
    let (inner, outer, height) = (RADIUS, 1.2 * RADIUS, 0.5 * RADIUS);
    let offset = 1.5 * D0;
    let slope = 0.5 * EPSILON0;
    let k = 2.7 * D0 * inner;
    let exact_v = |r: f64| (k * (outer / r).ln() - offset * (outer - r)) / slope;
    let charge = ANGLE * height * k;
    let differential = ANGLE * height * slope / (outer / inner).ln();
    let mut previous = f64::INFINITY;
    for n in [4, 8] {
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();
        for j in 0..=2 {
            for i in 0..=n {
                vertices.push([
                    inner + (outer - inner) * i as f64 / n as f64,
                    height * j as f64 / 2.0,
                ]);
            }
        }
        for j in 0..2 {
            for i in 0..n {
                let a = (j * (n + 1) + i) as u32;
                let b = a + (n + 1) as u32;
                triangles.extend([[a, a + 1, b + 1], [a, b + 1, b]]);
            }
        }
        let initial: Vec<_> = vertices
            .iter()
            .map(|x| exact_v(x[0]) - exact_v(inner))
            .collect();
        let mesh = TriComplex2::from_indexed_triangles(
            "nonlinear coaxial dielectric",
            vertices,
            triangles,
            Metric2::axisymmetric(ANGLE).unwrap(),
        )
        .unwrap();
        let p = problem(mesh, DielectricMaterial::Reversible(nonlinear()));
        let inner_trace: Vec<_> = (0..=2).map(|j| j * (n + 1)).collect();
        let outer_trace: Vec<_> = inner_trace.iter().map(|v| v + n).collect();
        let electrodes = p
            .prepare_electrodes(&[&inner_trace, &outer_trace], &mut || false)
            .unwrap();
        let field = electrodes
            .solve(
                &vec![0.0; 4 * n],
                &[
                    ElectrodeControl::Voltage(exact_v(inner)),
                    ElectrodeControl::Voltage(0.0),
                ],
                &initial,
                controls(),
                &mut || false,
            )
            .unwrap();
        let error = (field.terminal_charge_c()[0] / charge - 1.0).abs();
        assert!(
            error < 0.4 * previous,
            "n={n}, error={error:e}, previous={previous:e}"
        );
        previous = error;
        let tangent = field
            .differential_capacitance(differential_controls(), &mut || false)
            .unwrap();
        assert!((tangent.matrix_f()[0] / differential - 1.0).abs() < 0.003);
        for (x, &v) in p.mesh().vertices().iter().zip(field.relative_voltage_v()) {
            assert!((v - (exact_v(x[0]) - exact_v(inner))).abs() < 0.002 * exact_v(inner));
        }
    }
    assert!(previous < 3e-4);
}

#[test]
fn nonlinear_owner_retains_exact_radial_space_charge_moments() {
    // Lower triangle radii (0,R,R): shares (1/4,3/8,3/8).
    // Upper triangle radii (0,R,0): shares (1/4,1/2,1/4).
    let rho = 3e-6;
    let epsilon = 2e-11;
    let p = problem(
        cylinder(1, ANGLE, 1.0, 1.0),
        DielectricMaterial::isotropic(linear(epsilon)),
    );
    let electrodes = p
        .prepare_electrodes(&[&[0, 1], &[2, 3]], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[rho; 2],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(V0),
            ],
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let area = 0.5 * ANGLE * RADIUS * RADIUS;
    let charge0 = rho * p.mesh().face_measure(0).unwrap();
    let charge1 = rho * p.mesh().face_measure(1).unwrap();
    let flux = epsilon * V0 / GAP;
    let expected = [
        -area * flux / 3.0 - charge0 / 4.0 - charge1 / 4.0,
        -2.0 * area * flux / 3.0 - 3.0 * charge0 / 8.0,
        area * flux / 3.0 - charge1 / 4.0,
        2.0 * area * flux / 3.0 - 3.0 * charge0 / 8.0 - charge1 / 2.0,
    ];
    for (actual, expected) in field.nodal_boundary_charge_c().iter().zip(expected) {
        near(actual.unwrap(), expected, area * flux);
    }
    near(
        field.volume_charge_c(),
        charge0 + charge1,
        charge0 + charge1,
    );
    assert!(field.virtual_work(&[[0.0; 2]; 4], &mut || false).is_err());

    // Compare the complete solve, including a free interior row, against the
    // independently assembled Whitney-Hodge linear radial owner.
    let mesh = cylinder(2, ANGLE, 1.0, 1.0);
    let p = problem(mesh.clone(), DielectricMaterial::isotropic(linear(epsilon)));
    let native =
        AxisymmetricElectrostatics::new(mesh, vec![linear(epsilon); 8], limits(), &mut || false)
            .unwrap();
    let a = native
        .solve(
            &[rho; 8],
            &[(0, 0.0), (1, 0.0), (2, 0.0), (6, V0), (7, V0), (8, V0)],
            ElectrostaticControls {
                linear: linear_controls(),
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    let electrodes = p
        .prepare_electrodes(&[&[0, 1, 2], &[6, 7, 8]], &mut || false)
        .unwrap();
    let b = electrodes
        .solve(
            &[rho; 8],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(V0),
            ],
            &[0.0; 9],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in a.relative_voltage_v().iter().zip(b.relative_voltage_v()) {
        near(a, b, V0);
    }
    near(a.energy_j(), b.energy_j(), U0);
}

#[test]
fn radial_virtual_work_matches_resolved_voltage_and_charge_shape_derivatives() {
    let law = nonlinear();
    let p = problem(
        cylinder(2, ANGLE, 1.0, 1.0),
        DielectricMaterial::Reversible(law.clone()),
    );
    let electrodes = p
        .prepare_electrodes(&[&[0, 1, 2], &[6, 7, 8]], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[0.0; 8],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(2.0 * V0),
            ],
            &[0.0; 9],
            controls(),
            &mut || false,
        )
        .unwrap();
    let q = field.terminal_charge_c()[1];
    let radial: Vec<_> = p.mesh().vertices().iter().map(|x| [x[0], 0.0]).collect();
    let axial: Vec<_> = p.mesh().vertices().iter().map(|x| [0.0, x[1]]).collect();
    let dilation = p.mesh().vertices().to_vec();
    near(
        field.virtual_work(&radial, &mut || false).unwrap(),
        2.0 * field.coenergy_j(),
        U0,
    );
    near(
        field.virtual_work(&axial, &mut || false).unwrap(),
        -field.energy_j(),
        U0,
    );
    near(
        field.virtual_work(&dilation, &mut || false).unwrap(),
        2.0 * field.coenergy_j() - field.energy_j(),
        U0,
    );
    assert_eq!(
        field
            .virtual_work(&[[0.0, 1e10]; 9], &mut || false)
            .unwrap(),
        0.0
    );
    let h = 2e-5;
    for (radial_change, axial_change, displacement) in [
        (1.0, 0.0, &radial),
        (0.0, 1.0, &axial),
        (1.0, 1.0, &dilation),
    ] {
        for fixed_charge in [false, true] {
            let functional = |t: f64| {
                let perturbed = problem(
                    cylinder(2, ANGLE, 1.0 + t * radial_change, 1.0 + t * axial_change),
                    DielectricMaterial::Reversible(law.clone()),
                );
                let electrodes = perturbed
                    .prepare_electrodes(&[&[0, 1, 2], &[6, 7, 8]], &mut || false)
                    .unwrap();
                let boundary = [
                    ElectrodeControl::Voltage(0.0),
                    if fixed_charge {
                        ElectrodeControl::Charge(q)
                    } else {
                        ElectrodeControl::Voltage(2.0 * V0)
                    },
                ];
                let solved = electrodes
                    .solve(
                        &[0.0; 8],
                        &boundary,
                        field.relative_voltage_v(),
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                if fixed_charge {
                    -solved.energy_j()
                } else {
                    solved.coenergy_j()
                }
            };
            let derivative = (functional(h) - functional(-h)) / (2.0 * h);
            let work = field.virtual_work(displacement, &mut || false).unwrap();
            assert!(
                (derivative - work).abs() < 2e-8 * field.coenergy_j(),
                "{derivative:e} != {work:e}"
            );
        }
    }
}

#[test]
fn radial_geometry_axis_electrodes_and_work_bounds_refuse_without_partial_fields() {
    let material = DielectricMaterial::Reversible(nonlinear());
    assert!(
        AxisymmetricNonlinearElectrostatics::new(
            slab_mesh(),
            vec![material.clone(); 4],
            limits(),
            &mut || false
        )
        .is_err()
    );
    let mesh = cylinder(2, ANGLE, 1.0, 1.0);
    assert!(
        PlanarNonlinearElectrostatics::new(
            mesh.clone(),
            vec![material.clone(); 8],
            limits(),
            &mut || false
        )
        .is_err()
    );
    let p = problem(mesh, material);
    assert!(
        p.prepare_electrodes(&[&[0, 3, 6], &[2, 5, 8]], &mut || false)
            .is_err()
    );
    let electrodes = p
        .prepare_electrodes(&[&[0, 1, 2], &[6, 7, 8]], &mut || false)
        .unwrap();
    let boundary = [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Voltage(2.0 * V0),
    ];
    let mut limited = controls();
    limited.max_operator_applications = 0;
    assert!(matches!(
        electrodes.solve(&[0.0; 8], &boundary, &[0.0; 9], limited, &mut || false),
        Err(NonlinearElectrostaticError::OperatorBudget)
    ));
    assert!(
        electrodes
            .solve(&[0.0; 8], &boundary, &[0.0; 9], controls(), &mut || true)
            .is_err()
    );
    let field = electrodes
        .solve(&[0.0; 8], &boundary, &[0.0; 9], controls(), &mut || false)
        .unwrap();
    assert!(field.virtual_work(&[[1.0, 0.0]; 9], &mut || false).is_err());
    assert!(field.virtual_work(&[[0.0; 2]; 9], &mut || true).is_err());
    assert!(
        field
            .differential_capacitance(differential_controls(), &mut || true)
            .is_err()
    );
}
