//! G1/G3/G4: actual electrode currents and source-bound dissipation on the mesh.
use super::*;
use fs_couple::piezoelectric::spatial::transient::{
    TetPiezoelectricTransient as Body,
    load::{ElectrodeDrive as Drive, ResistiveElectrodeLoad as Load},
};
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, resolve_uniform_conductor,
};
use std::sync::Arc;

fn resistor(resistance: f64, electrodes: [usize; 2]) -> Load {
    let length = 0.01;
    let area = 1e-6;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: resistance * area / length,
                dims: ELECTRICAL_RESISTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored uniform electrode resistor".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "uniform geometry".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    Load {
        electrodes,
        conductor: Arc::new(
            resolve_uniform_conductor(
                &card,
                &QueryPoint::new().with("T", 300.0).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                length,
                area,
            )
            .unwrap(),
        ),
    }
}

#[test]
fn g1_resistive_electrode_relaxation_matches_rc_decay_and_retains_its_heat_source() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let capacitance = EPSILON[2] * LENGTH;
        let v0 = 5.0;
        let resistance = 1e4;
        let tau = resistance * capacitance;
        let load = resistor(resistance, [1, 0]);
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let mut body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, 0.0, Condition::Charge(capacitance * v0)),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let initial = body.state().total_energy_j();
            let mut heat = 0.0;
            for _ in 0..count {
                let (candidate, report) = body
                    .step_with_loads_candidate(
                        0.7 * tau / count as f64,
                        &[Drive::Voltage(0.0), Drive::Current(0.0)],
                        std::slice::from_ref(&load),
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap();
                assert!(report.joule_heat_j > 0.0);
                near(report.electrical_source_work_j, 0.0, 1e-22, 0.0);
                near(report.energy_residual_j, 0.0, 1e-20, 0.0);
                assert!(Arc::ptr_eq(
                    &report.loads[0].load.conductor,
                    &load.conductor
                ));
                heat += report.joule_heat_j;
                body = candidate.into_parts().0;
            }
            let voltage = body.state().electrodes()[1].voltage_v;
            errors.push((voltage / v0 - (-0.7_f64).exp()).abs());
            near(body.state().total_energy_j() + heat, initial, 1e-19, 1e-11);
            assert_eq!(body.state().kinetic_energy_j(), 0.0);
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{errors:?}"
            );
        }
        assert!(errors[2] < 2e-5, "{errors:?}");
    });
}

#[test]
fn g1_external_current_and_voltage_supply_actual_body_charge_and_resistor_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let capacitance = EPSILON[2] * LENGTH;
        let load = resistor(1e4, [1, 0]);
        let dt = 0.01 * capacitance * load.conductor.resistance_ohm();
        for current in [-3e-4, 3e-4] {
            let body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let (candidate, report) = body
                .step_with_loads_candidate(
                    dt,
                    &[Drive::Voltage(0.0), Drive::Current(current)],
                    std::slice::from_ref(&load),
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let endpoint = candidate.body().state();
            let v = endpoint.electrodes()[1].voltage_v;
            let exact = current * load.conductor.resistance_ohm() * (1.0 - (-0.01_f64).exp());
            near(v, exact, 3e-7, 0.0);
            near(
                endpoint.electrodes()[1].charge_c + report.loads[0].transferred_charge_c,
                current * dt,
                1e-22,
                1e-11,
            );
            near(
                report.electrical_source_work_j,
                0.5 * v * current * dt,
                1e-22,
                1e-11,
            );
            assert!(report.joule_heat_j > 0.0);
        }
        let voltage = 4.0;
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Voltage(voltage)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let (candidate, report) = body
            .step_with_loads_candidate(
                dt,
                &[Drive::Voltage(0.0), Drive::Voltage(voltage)],
                std::slice::from_ref(&load),
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let expected = dt * voltage * voltage / load.conductor.resistance_ohm();
        near(report.joule_heat_j, expected, 1e-22, 1e-11);
        near(report.electrical_source_work_j, expected, 1e-22, 1e-11);
        near(
            candidate.body().state().total_energy_j(),
            body.state().total_energy_j(),
            1e-22,
            1e-11,
        );
        let mut reversed = load.clone();
        reversed.electrodes.reverse();
        let (other, opposite) = body
            .step_with_loads_candidate(
                dt,
                &[Drive::Voltage(0.0), Drive::Voltage(voltage)],
                &[reversed],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        near(
            opposite.loads[0].midpoint_current_a,
            -report.loads[0].midpoint_current_a,
            0.0,
            1e-14,
        );
        near(opposite.joule_heat_j, report.joule_heat_j, 0.0, 1e-14);
        assert_eq!(
            other.body().state().voltage_v(),
            candidate.body().state().voltage_v()
        );
    });
}

#[test]
fn g1_mesh_motion_powers_the_actual_resistive_load_without_artificial_damping() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, 0.0))
            .collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let mut displacement = zero.clone();
        displacement[node][2] = 1e-8;
        let charge = vec![0.0; n];
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &displacement,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let initial = body.state().total_energy_j();
        let mut open = body.clone();
        let load = resistor(1e4, [1, 0]);
        let mut heat = 0.0;
        for _ in 0..30 {
            let (candidate, report) = body
                .step_with_loads_candidate(
                    1e-9,
                    &[Drive::Voltage(0.0), Drive::Current(0.0)],
                    std::slice::from_ref(&load),
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            near(report.electrical_source_work_j, 0.0, 1e-22, 0.0);
            heat += report.joule_heat_j;
            body = candidate.into_parts().0;
            open = open
                .step_candidate(
                    1e-9,
                    &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        assert!(heat > 1e-6 * initial);
        assert!(body.state().kinetic_energy_j() > 0.0);
        near(body.state().total_energy_j() + heat, initial, 1e-15, 1e-9);
        near(open.state().total_energy_j(), initial, 1e-15, 1e-9);
        assert!(body.state().total_energy_j() < open.state().total_energy_j());
    });
}

#[test]
fn g4_loads_refuse_foreign_gauges_switching_and_cancellation_without_state_change() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let k = cube.nodes.len();
        let mut nodes = cube.nodes.clone();
        nodes.extend(cube.nodes.iter().map(|p| [p[0] + 3.0 * LENGTH, p[1], p[2]]));
        let mut tets = cube.tets.clone();
        tets.extend(cube.tets.iter().map(|t| t.map(|i| i + k)));
        let bottom2 = cube.bottom.iter().map(|i| i + k).collect::<Vec<_>>();
        let top2 = cube.top.iter().map(|i| i + k).collect::<Vec<_>>();
        let electrodes = [
            Electrode {
                nodes: &cube.bottom,
                condition: Condition::Voltage(0.0),
            },
            Electrode {
                nodes: &cube.top,
                condition: Condition::Charge(0.0),
            },
            Electrode {
                nodes: &bottom2,
                condition: Condition::Voltage(0.0),
            },
            Electrode {
                nodes: &top2,
                condition: Condition::Charge(0.0),
            },
        ];
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model = Model::assemble(&nodes, &tets, &materials, FRAME, budget(), cx).unwrap();
        let zero = vec![[0.0; 3]; nodes.len()];
        let charge = vec![0.0; nodes.len()];
        let fixed = (0..3 * nodes.len()).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let body = Body::new(
            &model,
            &fixed,
            &electrodes,
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let drives = [
            Drive::Voltage(0.0),
            Drive::Current(0.0),
            Drive::Voltage(0.0),
            Drive::Current(0.0),
        ];
        for endpoints in [[1, 3], [1, 1], [1, 4]] {
            assert!(
                body.step_with_loads_candidate(
                    1e-8,
                    &drives,
                    &[resistor(1e4, endpoints)],
                    &zero,
                    &charge,
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        let mut switched = drives;
        switched[1] = Drive::Voltage(0.0);
        assert!(
            body.step_with_loads_candidate(1e-8, &switched, &[], &zero, &charge, controls(), cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            body.step_with_loads_candidate(
                1e-8,
                &drives,
                &[resistor(1e4, [1, 0])],
                &zero,
                &charge,
                controls(),
                cx
            ),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().total_energy_j(), 0.0);
    });
}
