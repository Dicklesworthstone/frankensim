//! G0/G1/G3/G4: exact radial source moments and analytic rotational fields.
//! Source cards are synthetic fixtures, not experimentally qualified materials.
mod common;

use common::{config, limits};
use fs_em::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticControls, ElectrostaticDisplacementTrace,
    ElectrostaticError, PlanarElectrostatics,
    axisymmetric::{AxisymmetricElectrostaticSources, AxisymmetricElectrostatics},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::planar::PlanarElectricalFrame,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        planar::{PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, PlanarDielectricLaw},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_phs::circuit::Branch;
use fs_rep_mesh::{Metric2, TriComplex2};

fn material(values: &[(&str, f64)]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, value) in values {
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
                    source: "synthetic meridional dielectric".into(),
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
            process: "numerical fixture".into(),
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
    material(&[(ELECTRIC_PERMITTIVITY_PROPERTY, epsilon)])
}

fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: config(),
        charge_tolerance_c: 1e-15,
        energy_tolerance_j: 1e-15,
    }
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() < 2e-9 * scale, "{a} != {b}, scale {scale}");
}

fn rectangle(a: f64, b: f64, h: f64, nr: usize, nz: usize, angle: f64) -> TriComplex2 {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    for j in 0..=nz {
        for i in 0..=nr {
            vertices.push([a + (b - a) * i as f64 / nr as f64, h * j as f64 / nz as f64]);
        }
    }
    for j in 0..nz {
        for i in 0..nr {
            let v = (j * (nr + 1) + i) as u32;
            let top = v + (nr + 1) as u32;
            faces.push([v, v + 1, top + 1]);
            faces.push([v, top + 1, top]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "axisymmetric dielectric fixture",
        vertices,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn problem(mesh: TriComplex2, epsilon: f64) -> AxisymmetricElectrostatics {
    let cells = mesh.faces().len();
    AxisymmetricElectrostatics::new(
        mesh,
        vec![dielectric(epsilon); cells],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn coaxial_capacitance_and_field_converge_to_logarithmic_solution() {
    let epsilon = 4e-11;
    let angle = core::f64::consts::TAU;
    let expected = angle * 2.0 * epsilon / 3.0_f64.ln();
    let mut previous_error = f64::INFINITY;
    for n in [8, 16, 32] {
        let p = problem(rectangle(1.0, 3.0, 2.0, n, n, angle), epsilon);
        let inner: Vec<_> = (0..=n).map(|j| j * (n + 1)).collect();
        let outer: Vec<_> = inner.iter().map(|v| v + n).collect();
        let cap = p
            .capacitance(&inner, &outer, controls(), &mut || false)
            .unwrap();
        let error = (cap.capacitance_f() - expected).abs() / expected;
        assert!(
            error < 0.4 * previous_error,
            "n={n}: {error}, {previous_error}"
        );
        previous_error = error;
        let field = cap.unit_voltage_field();
        for (&[r, _], &v) in p.mesh().vertices().iter().zip(field.voltage_v()) {
            let exact = (3.0 / r).ln() / 3.0_f64.ln();
            assert!((v - exact).abs() < 2.0 / (n * n) as f64);
        }
        near(field.energy_j(), 0.5 * cap.capacitance_f(), expected);
        assert_eq!(
            cap.circuit_branch(),
            Branch::Capacitor {
                farads: cap.capacitance_f()
            }
        );
    }
    assert!(previous_error < 6e-4);
}

#[test]
fn layered_coaxial_dielectrics_match_series_capacitance() {
    let (nr, nz) = (40, 10);
    let mesh = rectangle(1.0, 3.0, 1.3, nr, nz, core::f64::consts::TAU);
    let low = dielectric(2e-11);
    let high = dielectric(9e-11);
    let materials = mesh
        .faces()
        .iter()
        .map(|face| {
            let r = face
                .iter()
                .map(|&v| mesh.vertices()[v as usize][0])
                .sum::<f64>()
                / 3.0;
            if r < 2.0 { low.clone() } else { high.clone() }
        })
        .collect();
    let p = AxisymmetricElectrostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let inner: Vec<_> = (0..=nz).map(|j| j * (nr + 1)).collect();
    let outer: Vec<_> = inner.iter().map(|v| v + nr).collect();
    let cap = p
        .capacitance(&inner, &outer, controls(), &mut || false)
        .unwrap();
    let expected = core::f64::consts::TAU * 1.3 / (2.0_f64.ln() / 2e-11 + 1.5_f64.ln() / 9e-11);
    assert!((cap.capacitance_f() / expected - 1.0).abs() < 6e-4);
    near(cap.unit_voltage_field().boundary_charge_c(), 0.0, expected);
}

#[test]
fn axial_plate_field_is_exact_at_axis_and_scales_with_angular_span() {
    let epsilon = 3e-11;
    for angle in [core::f64::consts::TAU, 0.7] {
        let p = problem(rectangle(0.0, 2.0, 0.4, 2, 2, angle), epsilon);
        let cap = p
            .capacitance(&[6, 7, 8], &[0, 1, 2], controls(), &mut || false)
            .unwrap();
        // One-dimensional axial field, but the physical electrode is a disk.
        let expected = epsilon * angle * 2.0_f64.powi(2) / (2.0 * 0.4);
        near(cap.capacitance_f(), expected, expected);
        for e in cap.unit_voltage_field().electric_field_v_m() {
            near(e[0], 0.0, 2.5);
            near(e[1], -2.5, 2.5);
        }
    }
}

#[test]
fn volume_charge_uses_exact_radial_p1_moments() {
    let mesh = TriComplex2::from_indexed_triangles(
        "radial source moments",
        vec![[1.0, 0.0], [3.0, 0.0], [2.0, 2.0]],
        vec![[0, 1, 2]],
        Metric2::axisymmetric(0.7).unwrap(),
    )
    .unwrap();
    let p = problem(mesh, 2e-11);
    let rho = 1e-10;
    let field = p
        .solve(
            &[rho],
            &[(0, 0.0), (1, 0.0), (2, 0.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let total = rho * 2.0 * 0.7 * 2.0;
    // Integral lambda_i*r = area*(sum(r_j)+r_i)/12. Equal thirds is wrong.
    for (q, numerator) in field.nodal_boundary_charge_c().iter().zip([7.0, 9.0, 8.0]) {
        near(q.unwrap(), -total * numerator / 24.0, total);
    }
    near(field.volume_charge_c(), total, total);
    near(field.charge_balance_residual_c(), 0.0, total);
}

#[test]
fn charged_solid_cylinder_converges_with_natural_axis_regularity() {
    let epsilon = 2e-11;
    let rho = 4.0 * epsilon;
    let angle = core::f64::consts::TAU;
    let exact_energy = angle * rho * rho / (32.0 * epsilon);
    let mut previous_error = f64::INFINITY;
    for n in [8, 16, 32] {
        let p = problem(rectangle(0.0, 1.0, 1.0, n, n, angle), epsilon);
        let boundary: Vec<_> = (0..=n).map(|j| (j * (n + 1) + n, 0.0)).collect();
        let field = p
            .solve(&vec![rho; 2 * n * n], &boundary, controls(), &mut || false)
            .unwrap();
        let mut maximum_error = 0.0_f64;
        for (&[r, _], &v) in p.mesh().vertices().iter().zip(field.voltage_v()) {
            maximum_error = maximum_error.max((v - (1.0 - r * r)).abs());
        }
        assert!(
            maximum_error < 0.45 * previous_error,
            "n={n}: {maximum_error}, {previous_error}"
        );
        previous_error = maximum_error;
        near(field.volume_charge_c(), rho * angle / 2.0, rho * angle);
        near(field.boundary_charge_c(), -rho * angle / 2.0, rho * angle);
        assert!((field.energy_j() / exact_energy - 1.0).abs() < 2.0 / (n * n) as f64);
    }
    assert!(previous_error < 0.002);
}

#[test]
fn radial_surface_moments_reproduce_constant_axial_displacement() {
    let epsilon = 4e-11;
    let p = problem(rectangle(0.0, 2.0, 0.5, 1, 1, 1.2), epsilon);
    let rho = [0.0; 2];
    let trace = [ElectrostaticDisplacementTrace {
        vertices: [0, 1],
        outward_displacement_c_m2: -3.0 * epsilon,
    }];
    let field = p
        .solve_with_sources(
            AxisymmetricElectrostaticSources {
                displacement_traces: &trace,
                ..AxisymmetricElectrostaticSources::volume_charge(&rho)
            },
            &[(2, 0.0), (3, 0.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    for e in field.electric_field_v_m() {
        near(e[0], 0.0, 3.0);
        near(e[1], 3.0, 3.0);
    }
    near(field.voltage_v()[0], 1.5, 1.5);
    near(field.voltage_v()[1], 1.5, 1.5);
    let outward = -3.0 * epsilon * 1.2 * 2.0;
    near(field.outward_displacement_c(), outward, outward.abs());
    near(field.boundary_charge_c(), outward, outward.abs());
}

#[test]
fn charged_disk_sheet_preserves_jump_and_electrode_reactions() {
    let epsilon = 3e-11;
    let sigma = 2e-10;
    let p = problem(rectangle(0.0, 2.0, 1.0, 1, 2, 0.9), epsilon);
    let rho = [0.0; 4];
    let sheets = [ElectrostaticChargeSheet {
        vertices: [2, 3],
        free_charge_density_c_m2: sigma,
    }];
    let field = p
        .solve_with_sources(
            AxisymmetricElectrostaticSources {
                charge_sheets: &sheets,
                ..AxisymmetricElectrostaticSources::volume_charge(&rho)
            },
            &[(0, 0.0), (1, 0.0), (4, 0.0), (5, 0.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (cell, d) in field.displacement_c_m2().iter().enumerate() {
        near(d[0], 0.0, sigma);
        near(
            d[1],
            if cell < 2 { -sigma / 2.0 } else { sigma / 2.0 },
            sigma,
        );
    }
    near(
        field.voltage_v()[2],
        sigma / (4.0 * epsilon),
        sigma / epsilon,
    );
    near(
        field.voltage_v()[3],
        sigma / (4.0 * epsilon),
        sigma / epsilon,
    );
    near(field.sheet_charge_c(), sigma * 0.9 * 2.0, sigma * 2.0);
    near(
        field.boundary_charge_c(),
        -field.sheet_charge_c(),
        sigma * 2.0,
    );
}

#[test]
fn impressed_axial_polarization_has_open_circuit_field_without_free_charge() {
    let epsilon = 2e-11;
    let pz = 7e-11;
    let p = problem(rectangle(0.0, 1.4, 0.3, 1, 1, 0.8), epsilon);
    let polarization = [[0.0, pz]; 2];
    let sources = AxisymmetricElectrostaticSources {
        impressed_polarization_c_m2: Some(&polarization),
        ..AxisymmetricElectrostaticSources::volume_charge(&[0.0; 2])
    };
    let field = p
        .solve_with_sources(sources, &[(2, 0.0), (3, 0.0)], controls(), &mut || false)
        .unwrap();
    let volume = 0.8 * 1.4_f64.powi(2) * 0.3 / 2.0;
    for (e, d) in field
        .electric_field_v_m()
        .iter()
        .zip(field.displacement_c_m2())
    {
        near(e[0], 0.0, pz / epsilon);
        near(e[1], -pz / epsilon, pz / epsilon);
        near(d[0], 0.0, pz);
        near(d[1], 0.0, pz);
    }
    near(
        field.energy_j(),
        volume * pz * pz / (2.0 * epsilon),
        volume * pz * pz / epsilon,
    );
    near(
        field.polarization_pairing_j(),
        -2.0 * field.energy_j(),
        field.energy_j(),
    );
    near(field.boundary_charge_c(), 0.0, pz * volume / 0.3);
    let biased = p
        .solve_with_sources(sources, &[(2, 1e20), (3, 1e20)], controls(), &mut || false)
        .unwrap();
    assert_eq!(field.relative_voltage_v(), biased.relative_voltage_v());
    assert_eq!(field.electric_field_v_m(), biased.electric_field_v_m());
    assert_eq!(field.energy_j(), biased.energy_j());
}

#[test]
fn oblique_meridional_tensor_retains_source_frame_and_radial_divergence() {
    let m = material(&[
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[0], 2e-11),
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1], 6e-11),
    ]);
    let frame = PlanarElectricalFrame::new(
        dielectric(1e-11).card_identity(),
        dielectric(3e-11).card_identity(),
        core::f64::consts::FRAC_PI_4,
    )
    .unwrap();
    let law = PlanarDielectricLaw::Orthotropic(frame);
    let mesh = rectangle(1.0, 2.0, 1.0, 2, 2, 0.7);
    // epsilon=[[4,-2],[-2,4]]e-11; E=(0.5,1) yields D=(0,3)e-11.
    // Dr=0 makes div_cyl(D)=0, despite the nonzero radial electric field.
    let boundary: Vec<_> = mesh
        .boundary_trace()
        .unwrap()
        .vertices
        .into_iter()
        .map(|v| {
            let [r, z] = mesh.vertices()[v as usize];
            (v as usize, -0.5 * r - z)
        })
        .collect();
    let p = AxisymmetricElectrostatics::new_with_laws(
        mesh,
        vec![m.clone(); 8],
        vec![law; 8],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = p
        .solve(&[0.0; 8], &boundary, controls(), &mut || false)
        .unwrap();
    for cell in 0..8 {
        near(field.electric_field_v_m()[cell][0], 0.5, 1.0);
        near(field.electric_field_v_m()[cell][1], 1.0, 1.0);
        near(field.displacement_c_m2()[cell][0], 0.0, 3e-11);
        near(field.displacement_c_m2()[cell][1], 3e-11, 3e-11);
        assert_eq!(&p.materials()[cell], &m);
    }
    near(field.energy_j(), 0.5 * 3e-11 * 0.7 * 1.5, 3e-11);
    assert_eq!(p.permittivity_laws(), &[law; 8]);
}

#[test]
fn metric_axis_sources_budgets_and_cancellation_refuse_without_publication() {
    let mesh = rectangle(0.0, 1.0, 1.0, 1, 1, 1.0);
    assert!(matches!(
        PlanarElectrostatics::new(mesh, vec![dielectric(2e-11); 2], limits(), &mut || false,),
        Err(ElectrostaticError::Invalid(_))
    ));
    let planar = common::layered(0.5, 1.0, 0.2);
    assert!(matches!(
        AxisymmetricElectrostatics::new(planar, vec![dielectric(2e-11); 4], limits(), &mut || {
            false
        },),
        Err(ElectrostaticError::Invalid(_))
    ));
    let p = problem(rectangle(0.0, 1.0, 1.0, 1, 1, 1.0), 2e-11);
    assert!(matches!(
        p.solve(&[0.0; 2], &[(0, 0.0), (2, 0.0)], controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    let axis_flux = [ElectrostaticDisplacementTrace {
        vertices: [0, 2],
        outward_displacement_c_m2: 1e-10,
    }];
    assert!(matches!(
        p.solve_with_sources(
            AxisymmetricElectrostaticSources {
                displacement_traces: &axis_flux,
                ..AxisymmetricElectrostaticSources::volume_charge(&[0.0; 2])
            },
            &[(1, 0.0), (3, 0.0)],
            controls(),
            &mut || false,
        ),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        p.capacitance(&[2, 3], &[0, 1], controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
    let mut budget = limits();
    budget.max_cells = 1;
    assert!(matches!(
        AxisymmetricElectrostatics::new(
            rectangle(0.0, 1.0, 1.0, 1, 1, 1.0),
            vec![dielectric(2e-11); 2],
            budget,
            &mut || false,
        ),
        Err(ElectrostaticError::Invalid(_))
    ));
    let first = p
        .capacitance(&[2, 3], &[0, 1], controls(), &mut || false)
        .unwrap();
    let repeat = p
        .capacitance(&[2, 3], &[0, 1], controls(), &mut || false)
        .unwrap();
    assert_eq!(first.capacitance_f(), repeat.capacitance_f());
    assert_eq!(
        first.unit_voltage_field().electric_field_v_m(),
        repeat.unit_voltage_field().electric_field_v_m()
    );
}
