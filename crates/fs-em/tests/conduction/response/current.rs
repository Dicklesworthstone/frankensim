//! G1/G3/G4: total current, source response and original spatial Joule work.
use super::*;

const POSITIVE: &[usize] = &[2, 5];
const NEGATIVE: &[usize] = &[0, 3];

fn current_model(first: f64) -> PlanarConductionResponse {
    specimen(
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            first,
            0.2,
            20.0,
        ),
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            3.0,
            0.4,
            20.0,
        ),
        0.3,
    )
}
fn solve(
    model: &PlanarConductionResponse,
    current: f64,
) -> fs_em::conduction::response::ConductionResponseCandidate {
    model
        .solve_current_candidate(
            POSITIVE,
            NEGATIVE,
            current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap()
}

#[test]
fn authored_current_closes_series_fields_and_local_work_for_both_signs_and_material_swaps() {
    let mut voltages = Vec::new();
    for first in [1.0, 2.0] {
        let model = current_model(first);
        let original = model.conductors().to_vec();
        for current in [0.6_f64, -0.6, 0.0] {
            let next = solve(&model, current);
            let r = next.report();
            let j = current.abs() / 0.3;
            let e = [(first, 0.2), (3.0, 0.4)]
                .map(|(a, b)| 2.0 * j / (a + (a * a + 4.0 * b * j).sqrt()));
            let voltage = current.signum() * 0.5 * (e[0] + e[1]);
            near(r.voltage_v[2], voltage);
            near(r.voltage_v[5], voltage);
            near(r.voltage_v[1], current.signum() * 0.5 * e[0]);
            near(r.power_w, current * voltage);
            near(
                POSITIVE
                    .iter()
                    .map(|&v| r.nodal_boundary_current_a[v].unwrap())
                    .sum(),
                current,
            );
            near(
                NEGATIVE
                    .iter()
                    .map(|&v| r.nodal_boundary_current_a[v].unwrap())
                    .sum(),
                -current,
            );
            for triangle in 0..4 {
                let i = triangle / 2;
                near(r.electric_field_v_m[triangle][0], -current.signum() * e[i]);
                near(r.current_density_a_m2[triangle][0], -current.signum() * j);
                near(r.current_density_a_m2[triangle][1], 0.0);
                near(r.cell_power_w[triangle], 0.075 * e[i] * j);
                let source = &r.evaluated[triangle].conductor;
                near(source.material().query().axes()["operating"], e[i]);
                assert_eq!(
                    source.source_property().answer().receipt.selected,
                    original[triangle]
                        .source_property()
                        .answer()
                        .receipt
                        .selected
                );
            }
            if current == 0.6 {
                voltages.push(voltage);
            }
            if current != 0.0 {
                assert!(r.field_solves > 1);
                assert!(r.linear.is_none());
            }
        }
        // Re-solving at a new authored current must not scale a nonlinear field.
        let low = solve(&model, 0.3);
        let high = solve(low.problem(), 0.6);
        let fresh = solve(&model, 0.6);
        assert!((high.report().voltage_v[2] - 2.0 * low.report().voltage_v[2]).abs() > 0.01);
        near(high.report().voltage_v[2], fresh.report().voltage_v[2]);
        assert_eq!(model.conductors(), original);
    }
    assert!(voltages[1] < voltages[0]);
}

#[test]
fn current_dependent_resistivity_uses_actual_current_density_and_authored_source_work() {
    let model = specimen(
        response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            1.0,
            0.2,
            20.0,
        ),
        response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            2.0,
            0.4,
            20.0,
        ),
        0.3,
    );
    for current in [0.6_f64, -0.6] {
        let next = solve(&model, current);
        let j = current.abs() / 0.3;
        let e = [j * (1.0 + 0.2 * j), j * (2.0 + 0.4 * j)];
        let voltage = current.signum() * 0.5 * (e[0] + e[1]);
        near(next.report().voltage_v[2], voltage);
        near(next.report().power_w, current * voltage);
        for (i, evaluated) in next.report().evaluated.iter().enumerate() {
            near(
                evaluated.conductor.material().query().axes()["operating"],
                j,
            );
            near(next.report().cell_power_w[i], 0.075 * j * e[i / 2]);
        }
    }
}

#[test]
fn equipotential_electrodes_solve_nonuniform_reaction_in_parallel_materials() {
    let mesh = TriComplex2::from_indexed_triangles(
        "parallel-current-controlled-conductors",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
        ],
        vec![[0, 1, 3], [0, 3, 2], [2, 3, 5], [2, 5, 4]],
        Metric2::planar(0.3).unwrap(),
    )
    .unwrap();
    let a = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        20.0,
    );
    let b = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        3.0,
        0.4,
        20.0,
    );
    let model = PlanarConductionResponse::new(
        mesh,
        vec![a.clone(), a, b.clone(), b],
        limits(),
        &mut || false,
    )
    .unwrap();
    let next = model
        .solve_current_candidate(
            &[1, 3, 5],
            &[0, 2, 4],
            0.6,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    // I/d=(a0+a1)E+(b0+b1)E², independently of electrode nodal reactions.
    let e = 4.0 / (4.0 + (16.0_f64 + 4.0 * 0.6 * 2.0).sqrt());
    let j = [(1.0 + 0.2 * e) * e, (3.0 + 0.4 * e) * e];
    for v in [1, 3, 5] {
        near(next.report().voltage_v[v], e);
    }
    near(
        next.report().nodal_boundary_current_a[1].unwrap(),
        0.15 * j[0],
    );
    near(
        next.report().nodal_boundary_current_a[5].unwrap(),
        0.15 * j[1],
    );
    assert!(j[1] > 2.0 * j[0]);
    for t in 0..4 {
        near(next.report().current_density_a_m2[t][0], -j[t / 2]);
    }
    near(next.report().power_w, 0.6 * e);
}

#[test]
fn constant_source_limit_uses_the_original_native_current_owner() {
    let native = slab(
        conductor(2.0, ElectricalPropertyBasis::Conductivity),
        conductor(4.0, ElectricalPropertyBasis::Conductivity),
        0.3,
    );
    let model = specimen(
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            2.0,
            0.0,
            20.0,
        ),
        response(
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            4.0,
            0.0,
            20.0,
        ),
        0.3,
    );
    for current in [-0.6, 0.0, 0.6] {
        let direct = native
            .solve_current(POSITIVE, NEGATIVE, current, controls(), &mut || false)
            .unwrap();
        let next = solve(&model, current);
        assert_eq!(next.report().field_solves, 1);
        assert_eq!(next.report().voltage_v, direct.voltage_v());
        assert_eq!(next.report().cell_power_w, direct.cell_power_w());
        near(direct.voltage_v()[2], current * 1.25);
        near(direct.power_w(), current * current * 1.25);
    }
}

#[test]
fn total_current_rejection_and_late_cancellation_keep_the_original_source_owner() {
    let model = current_model(1.0);
    let original = model.conductors().to_vec();
    for (positive, negative, current) in [
        (POSITIVE, NEGATIVE, f64::NAN),
        (&[2_usize][..], NEGATIVE, 0.6),
        (POSITIVE, POSITIVE, 0.6),
        (POSITIVE, NEGATIVE, 1e6),
    ] {
        assert!(
            model
                .solve_current_candidate(
                    positive,
                    negative,
                    current,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(matches!(
        model.solve_current_candidate(
            POSITIVE,
            NEGATIVE,
            0.6,
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
    let mut polls = 0;
    let clean = model
        .solve_current_candidate(
            POSITIVE,
            NEGATIVE,
            0.6,
            limits(),
            controls(),
            policy(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut late = 0;
    assert!(matches!(
        model.solve_current_candidate(
            POSITIVE,
            NEGATIVE,
            0.6,
            limits(),
            controls(),
            policy(),
            &mut || {
                late += 1;
                late == polls
            }
        ),
        Err(ConductionResponseError::Cancelled)
    ));
    assert_eq!(model.conductors(), original);
    assert_eq!(
        solve(&model, 0.6).report().voltage_v,
        clean.report().voltage_v
    );
}

#[test]
fn independent_power_gate_remains_active_under_loose_constitutive_tolerances() {
    let model = current_model(1.0);
    let next = model
        .solve_current_candidate(
            POSITIVE,
            NEGATIVE,
            0.6,
            limits(),
            controls(),
            ConductionResponseControls {
                absolute_conductivity_tolerance_s_m: 100.0,
                current_density_tolerance_a_m2: 100.0,
                ..policy()
            },
            &mut || false,
        )
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
    near(next.report().power_w, 0.6 * next.report().voltage_v[2]);
}
