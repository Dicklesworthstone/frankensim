//! Native spatial solutions against independently reduced nonlinear specimens.
use super::*;

#[path = "response/current.rs"]
mod current;
use fs_em::conduction::response::{
    ConductionResponseControls, ConductionResponseError, PlanarConductionResponse,
};
use fs_material::conductor::response::{ConductorMagnitude, ScalarConductorResponse};

fn response(
    basis: ElectricalPropertyBasis,
    coordinate: ConductorMagnitude,
    intercept: f64,
    slope: f64,
    maximum: f64,
) -> ScalarConductorResponse {
    let card = source(
        PropertyValue::Curve {
            abscissa: "operating".into(),
            abscissa_dims: coordinate.quantity().dims(),
            knots: vec![(0.0, intercept), (maximum, intercept + slope * maximum)],
            dims: basis.dimensions(),
        },
        ValidityDomain::unconstrained().with_quantity(
            "operating",
            coordinate.quantity(),
            0.0,
            maximum,
        ),
        basis,
    );
    let point = QueryPoint::new()
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap();
    ScalarConductorResponse::new(
        &card,
        resolve(&card, &point, basis),
        &[("operating", coordinate)],
    )
    .unwrap()
}
fn policy() -> ConductionResponseControls {
    ConductionResponseControls {
        maximum_iterations: 128,
        absolute_conductivity_tolerance_s_m: 1e-12,
        relative_conductivity_tolerance: 1e-12,
        current_density_tolerance_a_m2: 1e-10,
    }
}
fn specimen(
    a: ScalarConductorResponse,
    b: ScalarConductorResponse,
    depth: f64,
) -> PlanarConductionResponse {
    PlanarConductionResponse::new(
        layered(0.5, 1.0, depth),
        vec![a.clone(), a, b.clone(), b],
        limits(),
        &mut || false,
    )
    .unwrap()
}

// For sigma=a+b|E|, independent series continuity gives J=aE+bE².
// Invert each region analytically, then enforce the imposed total voltage.
fn series_current(a: [f64; 2], b: [f64; 2], voltage: f64) -> (f64, [f64; 2]) {
    let fields = |j: f64| [0, 1].map(|i| 2.0 * j / (a[i] + (a[i] * a[i] + 4.0 * b[i] * j).sqrt()));
    let (mut low, mut high) = (0.0, 100.0);
    for _ in 0..180 {
        let j = 0.5 * (low + high);
        if fields(j).iter().sum::<f64>() * 0.5 < voltage {
            low = j;
        } else {
            high = j;
        }
    }
    let j = 0.5 * (low + high);
    (j, fields(j))
}

#[test]
fn nonlinear_material_data_redistribute_native_current_field_and_local_joule_heat() {
    let mut currents = Vec::new();
    for first in [1.0, 2.0] {
        let a = [first, 3.0];
        let b = [0.2, 0.4];
        let depth = 0.3;
        let model = specimen(
            response(
                ElectricalPropertyBasis::Conductivity,
                ConductorMagnitude::ElectricField,
                a[0],
                b[0],
                4.0,
            ),
            response(
                ElectricalPropertyBasis::Conductivity,
                ConductorMagnitude::ElectricField,
                a[1],
                b[1],
                4.0,
            ),
            depth,
        );
        // Ordinary fixed-source admission must not freeze this operating axis.
        assert!(
            PlanarConduction::new(
                model.mesh().clone(),
                model.conductors().to_vec(),
                limits(),
                &mut || false
            )
            .is_err()
        );
        let original = model.conductors().to_vec();
        let (expected_j, expected_e) = series_current(a, b, 1.5);
        for voltage in [1.5, -1.5] {
            let solved = model
                .solve_candidate(
                    &ends(voltage),
                    &[],
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let report = solved.report();
            assert!(report.field_solves > 1);
            near(report.voltage_v[1], voltage.signum() * 0.5 * expected_e[0]);
            near(report.voltage_v[4], report.voltage_v[1]);
            for triangle in 0..4 {
                let region = triangle / 2;
                near(
                    report.electric_field_v_m[triangle][0],
                    -voltage.signum() * expected_e[region],
                );
                near(
                    report.current_density_a_m2[triangle][0],
                    -voltage.signum() * expected_j,
                );
                near(report.current_density_a_m2[triangle][1], 0.0);
                near(
                    report.cell_power_w[triangle],
                    0.25 * depth * expected_e[region] * expected_j,
                );
                let evaluated = &report.evaluated[triangle].conductor;
                near(
                    evaluated.material().query().axes()["operating"],
                    expected_e[region],
                );
                assert_eq!(
                    evaluated.source_property().answer().receipt.selected,
                    original[triangle]
                        .source_property()
                        .answer()
                        .receipt
                        .selected
                );
            }
            near(report.power_w, depth * expected_j * voltage.abs());
            near(report.boundary_current_a, 0.0);
            near(report.maximum_free_current_residual_a, 0.0);
            near(report.power_balance_residual_w, 0.0);
            assert_eq!(report.source_evaluations, 4 * report.field_solves);
            let (continued, _) = solved.into_parts();
            let next = continued
                .solve_candidate(
                    &ends(voltage),
                    &[],
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(next.report().power_w, depth * expected_j * voltage.abs());
        }
        assert_eq!(model.conductors(), original);
        currents.push(expected_j);
    }
    assert!(currents[1] > currents[0]);
}

#[test]
fn current_density_drives_original_resistivity_laws_without_unit_voltage_scaling() {
    let model = specimen(
        response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            1.0,
            0.1,
            8.0,
        ),
        response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            2.0,
            0.3,
            8.0,
        ),
        0.4,
    );
    for density in [3.5_f64, -3.5] {
        let offset = 17.0;
        let rho = [1.0 + 0.1 * density.abs(), 2.0 + 0.3 * density.abs()];
        let traces = [ConductionCurrentTrace {
            vertices: [2, 5],
            inward_current_density_a_m2: density,
        }];
        let next = model
            .solve_candidate(
                &[(0, offset), (3, offset)],
                &traces,
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap();
        let field = next.report();
        near(field.voltage_v[1], offset + 0.5 * density * rho[0]);
        near(
            field.voltage_v[2],
            offset + 0.5 * density * (rho[0] + rho[1]),
        );
        for triangle in 0..4 {
            near(field.current_density_a_m2[triangle][0], -density);
            near(
                field.electric_field_v_m[triangle][0],
                -density * rho[triangle / 2],
            );
            near(
                field.cell_power_w[triangle],
                0.1 * rho[triangle / 2] * density * density,
            );
        }
        near(field.power_w, 0.2 * (rho[0] + rho[1]) * density * density);
        near(field.nodal_applied_current_a[2].unwrap(), 0.2 * density);
        near(field.power_balance_residual_w, 0.0);
    }
}

#[test]
fn joint_e_j_source_coordinates_and_mesh_rotation_retain_isotropic_response() {
    let e = ConductorMagnitude::ElectricField;
    let j = ConductorMagnitude::CurrentDensity;
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = source(
        PropertyValue::Surface2D {
            axes: ["e".into(), "j".into()],
            axis_dims: [e.quantity().dims(), j.quantity().dims()],
            samples: [[0.0, 0.0], [4.0, 0.0], [0.0, 20.0], [4.0, 20.0]]
                .into_iter()
                .map(|p| (p, 2.0 + 0.2 * p[0] + 0.05 * p[1]))
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: basis.dimensions(),
        },
        ValidityDomain::unconstrained()
            .with_quantity("e", e.quantity(), 0.0, 4.0)
            .with_quantity("j", j.quantity(), 0.0, 20.0),
        basis,
    );
    let point = QueryPoint::new()
        .with_quantity("e", e.quantity(), 0.0)
        .unwrap()
        .with_quantity("j", j.quantity(), 0.0)
        .unwrap();
    let initial = resolve(&card, &point, basis);
    assert!(ScalarConductorResponse::new(&card, initial.clone(), &[("e", e)]).is_err());
    let law = ScalarConductorResponse::new(&card, initial, &[("e", e), ("j", j)]).unwrap();
    let expected_j = (2.0 + 0.2) / (1.0 - 0.05);
    for angle in [0.0_f64, 0.7] {
        let (s, c) = angle.sin_cos();
        let original = layered(0.5, 1.0, 0.2);
        let mesh = TriComplex2::from_indexed_triangles(
            "rotated nonlinear conductor",
            original
                .vertices()
                .iter()
                .map(|p| [c * p[0] - s * p[1], s * p[0] + c * p[1]])
                .collect(),
            original.faces().to_vec(),
            Metric2::planar(0.2).unwrap(),
        )
        .unwrap();
        let model =
            PlanarConductionResponse::new(mesh, vec![law.clone(); 4], limits(), &mut || false)
                .unwrap();
        let solved = model
            .solve_candidate(&ends(1.0), &[], limits(), controls(), policy(), &mut || {
                false
            })
            .unwrap();
        let report = solved.report();
        for triangle in 0..4 {
            near(report.electric_field_v_m[triangle][0], -c);
            near(report.electric_field_v_m[triangle][1], -s);
            near(report.current_density_a_m2[triangle][0], -c * expected_j);
            near(report.current_density_a_m2[triangle][1], -s * expected_j);
            near(
                report.evaluated[triangle]
                    .conductor
                    .material()
                    .query()
                    .axes()["j"],
                expected_j,
            );
        }
        near(report.power_w, 0.2 * expected_j);
    }
}

#[test]
fn zero_response_limit_and_independent_power_admission_use_original_native_fields() {
    let a = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        2.0,
        0.0,
        4.0,
    );
    let b = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        4.0,
        0.0,
        4.0,
    );
    let model = specimen(a, b, 0.2);
    let linear = slab(
        conductor(2.0, ElectricalPropertyBasis::Conductivity),
        conductor(4.0, ElectricalPropertyBasis::Conductivity),
        0.2,
    );
    let field = linear.solve(&ends(1.0), controls(), &mut || false).unwrap();
    let next = model
        .solve_candidate(&ends(1.0), &[], limits(), controls(), policy(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(next.report().field_solves, 1);
    assert_eq!(next.report().voltage_v, field.voltage_v());
    assert_eq!(next.report().cell_power_w, field.cell_power_w());
    let nonlinear = specimen(
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            1.0,
            0.2,
            4.0,
        ),
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            3.0,
            0.4,
            4.0,
        ),
        0.2,
    );
    let loose = ConductionResponseControls {
        absolute_conductivity_tolerance_s_m: 100.0,
        current_density_tolerance_a_m2: 100.0,
        ..policy()
    };
    let next = nonlinear
        .solve_candidate(&ends(1.5), &[], limits(), controls(), loose, &mut || false)
        .unwrap();
    assert!(next.report().field_solves > 1);
    assert!(
        next.report()
            .evaluated
            .iter()
            .map(|s| s.power_defect_w.abs())
            .sum::<f64>()
            <= controls().power_tolerance_w
    );
}

#[test]
fn source_budget_and_late_cancellation_refuse_without_replacing_accepted_materials() {
    let model = specimen(
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            1.0,
            0.2,
            4.0,
        ),
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            3.0,
            0.4,
            4.0,
        ),
        0.3,
    );
    let original = model.conductors().to_vec();
    let mut calls = 0;
    let baseline = model
        .solve_candidate(&ends(1.5), &[], limits(), controls(), policy(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut cancelled_calls = 0;
    assert!(matches!(
        model.solve_candidate(&ends(1.5), &[], limits(), controls(), policy(), &mut || {
            cancelled_calls += 1;
            cancelled_calls >= calls
        }),
        Err(ConductionResponseError::Cancelled)
    ));
    assert!(matches!(
        model.solve_candidate(
            &ends(1.5),
            &[],
            limits(),
            controls(),
            ConductionResponseControls {
                maximum_iterations: 1,
                ..policy()
            },
            &mut || false
        ),
        Err(ConductionResponseError::NotConverged { .. })
    ));
    assert!(
        model
            .solve_candidate(
                &ends(10.0),
                &[],
                limits(),
                controls(),
                policy(),
                &mut || false
            )
            .is_err()
    );
    let mut foreign = original.clone();
    foreign[0] = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        2.0,
        0.2,
        4.0,
    )
    .original()
    .clone();
    assert!(
        model
            .solve_with_materials_candidate(
                foreign,
                &ends(1.5),
                &[],
                limits(),
                controls(),
                policy(),
                &mut || false
            )
            .is_err()
    );
    assert_eq!(model.conductors(), original);
    let replay = model
        .solve_candidate(&ends(1.5), &[], limits(), controls(), policy(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(replay.report().voltage_v, baseline.report().voltage_v);
    assert_eq!(replay.report().evaluated, baseline.report().evaluated);
}
