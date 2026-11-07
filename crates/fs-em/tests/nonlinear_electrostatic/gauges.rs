//! Independent electrical gauges must not become physical dielectric fields.
use super::*;

#[test]
fn opposite_component_charge_defects_cannot_cancel_the_physical_balance_gate() {
    // Each fully electroded triangle encloses one unit volume. Nodal third-
    // shares round in opposite directions for equal and opposite source loads:
    // their global charge defect vanishes while each component retains one.
    let mesh = TriComplex2::from_indexed_triangles(
        "independent-enclosed-volume-charges",
        vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [0.0, 2.0],
            [2.0, 2.0],
            [0.0, 3.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    let problem = PlanarNonlinearElectrostatics::new(
        mesh,
        vec![DielectricMaterial::isotropic(linear(EPSILON0)); 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 1, 2], &[3, 4, 5]], &mut || false)
        .unwrap();
    let boundary = [ElectrodeControl::Voltage(0.0); 2];
    let rho = [1.5e-9, -1.5e-9];
    let admitted = electrodes
        .solve(&rho, &boundary, &[0.0; 6], controls(), &mut || false)
        .unwrap();
    assert_eq!(admitted.charge_balance_residual_c(), 0.0);
    assert!(admitted.maximum_component_charge_residual_c() > 1e-26);
    let mut strict = controls();
    strict.charge_tolerance_c = 1e-26;
    match electrodes.solve(&rho, &boundary, &[0.0; 6], strict, &mut || false) {
        Err(fs_em::electrostatic::nonlinear::NonlinearElectrostaticError::Electric(
            fs_em::electrostatic::ElectrostaticError::Balance {
                charge_residual_c, ..
            },
        )) => assert!(charge_residual_c.abs() > strict.charge_tolerance_c),
        Err(other) => panic!("expected independent-component charge refusal, got {other:?}"),
        Ok(_) => panic!("opposite component defects must not cancel the charge gate"),
    }
}

#[test]
fn independently_biased_components_preserve_charge_driven_fields_and_storage() {
    let problem = series_problem();
    // These four separate traces do not electrically join the two mesh pieces.
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2], &[4, 7], &[5, 6]], &mut || false)
        .unwrap();
    let solve = |first, second| {
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(first),
                    ElectrodeControl::Charge(2.5 * Q0),
                    ElectrodeControl::Voltage(second),
                    ElectrodeControl::Charge(0.5 * Q0),
                ],
                &[0.0; 8],
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    let ordinary = solve(0.0, 0.0);
    let biased = solve(1e250, -1e250);
    assert_eq!(biased.relative_voltage_v(), ordinary.relative_voltage_v());
    assert_eq!(biased.electric_field_v_m(), ordinary.electric_field_v_m());
    assert_eq!(biased.displacement_c_m2(), ordinary.displacement_c_m2());
    assert_eq!(biased.terminal_charge_c(), ordinary.terminal_charge_c());
    assert_eq!(biased.energy_j(), ordinary.energy_j());
    assert_eq!(biased.coenergy_j(), ordinary.coenergy_j());
    assert_eq!(biased.reference_terminal(), 0);
    assert_eq!(biased.reference_voltage_v(), 1e250);
    assert_eq!(
        biased.reference_terminal_by_vertex(),
        &[0, 0, 0, 0, 2, 2, 2, 2]
    );
    assert_eq!(&biased.reference_voltage_by_vertex_v()[..4], &[1e250; 4]);
    assert_eq!(&biased.reference_voltage_by_vertex_v()[4..], &[-1e250; 4]);
    for terminal in 0..4 {
        assert_eq!(
            biased.reference_terminal_for_terminal(terminal),
            Some(terminal / 2 * 2)
        );
        assert_eq!(
            biased.reference_voltage_for_terminal_v(terminal),
            Some(if terminal < 2 { 1e250 } else { -1e250 }),
        );
    }
    near(biased.terminal_voltage_v()[1], 2.0 * V0, V0);
    near(biased.terminal_voltage_v()[3], 0.5 * V0, V0);
    assert!(biased.maximum_component_charge_residual_c() <= controls().charge_tolerance_c);

    // A reference on the first body cannot anchor the independent second body.
    assert!(
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(2.5 * Q0),
                    ElectrodeControl::Charge(-0.5 * Q0),
                    ElectrodeControl::Charge(0.5 * Q0),
                ],
                &[0.0; 8],
                controls(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn an_unrelated_huge_bias_cannot_erase_a_resolved_voltage_or_differential_response() {
    let problem = series_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2], &[4, 7], &[5, 6]], &mut || false)
        .unwrap();
    let density = D0 / GAP;
    let rho = [0.2 * density, 0.2 * density, -0.3 * density, -0.3 * density];
    let solve = |first, second| {
        electrodes
            .solve(
                &rho,
                &[
                    ElectrodeControl::Voltage(first),
                    ElectrodeControl::Voltage(first),
                    ElectrodeControl::Voltage(second),
                    ElectrodeControl::Voltage(second + V0),
                ],
                &[0.0; 8],
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    let ordinary = solve(0.0, 0.0);
    let biased = solve(1e250, -1e12);
    assert_eq!(biased.relative_voltage_v(), ordinary.relative_voltage_v());
    assert_eq!(biased.electric_field_v_m(), ordinary.electric_field_v_m());
    assert_eq!(biased.terminal_charge_c(), ordinary.terminal_charge_c());
    assert_eq!(biased.energy_j(), ordinary.energy_j());
    near(biased.terminal_voltage_v()[3], V0, V0);
    near(
        biased.terminal_charge_c()[0] + biased.terminal_charge_c()[1],
        -0.2 * Q0,
        Q0,
    );
    near(
        biased.terminal_charge_c()[2] + biased.terminal_charge_c()[3],
        0.3 * Q0,
        Q0,
    );
    assert!(biased.maximum_component_charge_residual_c() <= controls().charge_tolerance_c);

    let actual = biased
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    let baseline = ordinary
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(actual.matrix_f(), baseline.matrix_f());
    let c = EPSILON0 * AREA / GAP;
    let expected = [
        2.0 * c,
        -2.0 * c,
        0.0,
        0.0,
        -2.0 * c,
        2.0 * c,
        0.0,
        0.0,
        0.0,
        0.0,
        c,
        -c,
        0.0,
        0.0,
        -c,
        c,
    ];
    for (&observed, &expected) in actual.matrix_f().iter().zip(&expected) {
        near(observed, expected, c);
    }
}
