//! G1/G3/G4: original spatial motion and material C/L/R histories share one solve.
#[path = "reactive/acoustic.rs"]
mod acoustic;

use super::*;
use fs_couple::piezoelectric::spatial::transient::load::reactive::{
    ReactiveElectrodeLoad as Reactive, ReactiveElectrodeState as State,
    ReactivePiezoelectricTransient as Shunted,
};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
};

fn storage_card(name: &str, dims: fs_qty::Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored reactive piezoelectric numerical fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "uniform fixed field geometry".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(crate) fn capacitor(c: f64) -> Arc<ResolvedCapacitor> {
    let (area, gap) = (1e-4, 1e-3);
    Arc::new(
        ResolvedCapacitor::resolve(
            &storage_card(
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                c * gap / area,
            ),
            &QueryPoint::new().with("T", 300.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            area,
            gap,
        )
        .unwrap(),
    )
}

pub(crate) fn inductor(l: f64) -> Arc<ResolvedInductor> {
    let (area, length, turns) = (1e-4, 1e-2, 10);
    Arc::new(
        ResolvedInductor::resolve(
            &storage_card(
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                l * length / (area * f64::from(turns * turns)),
            ),
            &QueryPoint::new().with("T", 300.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            area,
            length,
            turns,
        )
        .unwrap(),
    )
}

fn coordinate(state: State) -> f64 {
    match state {
        State::FluxLinkage(x) | State::Charge(x) => x,
    }
}

#[test]
fn g1_source_bound_lc_resonance_has_second_order_phase_and_conserves_all_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let ceramic_c = EPSILON[2] * LENGTH;
        let c = capacitor(2.0 * ceramic_c);
        let l = inductor(2e-3);
        let total_c = ceramic_c + c.capacitance_f();
        let omega = 1.0 / (l.inductance_h() * total_c).sqrt();
        let v0 = 3.0;
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, 0.0, Condition::Charge(ceramic_c * v0)),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let initial =
                body.state().total_energy_j() + c.stored_energy_j(c.capacitance_f() * v0).unwrap();
            let mut coupled = Shunted::new(
                body,
                vec![
                    Reactive::SeriesRl {
                        electrodes: [1, 0],
                        inductor: Arc::clone(&l),
                        conductor: None,
                    },
                    Reactive::Capacitor {
                        electrodes: [1, 0],
                        capacitor: Arc::clone(&c),
                    },
                ],
                vec![
                    State::FluxLinkage(0.0),
                    State::Charge(c.capacitance_f() * v0),
                ],
                vec![],
                controls(),
                cx,
            )
            .unwrap();
            for _ in 0..count {
                let candidate = coupled
                    .step_candidate(
                        0.7 / (omega * f64::from(count)),
                        &[Drive::Voltage(0.0), Drive::Current(0.0)],
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap();
                let step = candidate.step();
                assert_eq!(step.joule_heat_j, 0.0);
                near(step.electrical_source_work_j, 0.0, 1e-24, 0.0);
                near(step.energy_residual_j, 0.0, 1e-20, 0.0);
                near(
                    candidate.model().body().state().total_energy_j() + step.reactive_energy_j[1],
                    initial,
                    1e-20,
                    1e-10,
                );
                coupled = candidate.into_parts().0;
            }
            let v = coupled.body().state().electrodes()[1].voltage_v;
            let current = coordinate(coupled.states()[0]) / l.inductance_h();
            let scale_i = v0 * (total_c / l.inductance_h()).sqrt();
            errors.push((v / v0 - 0.7_f64.cos()).hypot(current / scale_i - 0.7_f64.sin()));
            near(
                coordinate(coupled.states()[1]),
                c.capacitance_f() * v,
                1e-22,
                1e-10,
            );
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{errors:?}"
            );
        }
        assert!(errors[2] < 4e-5, "{errors:?}");
    });
}

#[test]
fn g1_rl_shunt_matches_damped_solution_and_reverses_flux_without_changing_body_or_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let c = EPSILON[2] * LENGTH;
        let l = inductor(2e-3);
        let r = resistor(2e3, [1, 0]).conductor;
        let (v0, i0) = (3.0, 1e-4);
        let alpha = r.resistance_ohm() / (2.0 * l.inductance_h());
        let omega = (1.0 / (l.inductance_h() * c) - alpha * alpha).sqrt();
        let time = 0.7 / omega;
        let expected_v = (-alpha * time).exp()
            * (v0 * (omega * time).cos() + (alpha * v0 - i0 / c) / omega * (omega * time).sin());
        let expected_i = (-alpha * time).exp()
            * (i0 * (omega * time).cos()
                + (v0 / l.inductance_h() - alpha * i0) / omega * (omega * time).sin());
        let mut outcomes = Vec::new();
        for (orientation, reference) in [(1.0, 0.0), (-1.0, 128.0)] {
            let pair = if orientation > 0.0 { [1, 0] } else { [0, 1] };
            let body = Body::new(
                &model,
                &fixed,
                &electrodes(&cube, reference, Condition::Charge(c * v0)),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let initial =
                body.state().total_energy_j() + l.stored_energy_j(l.inductance_h() * i0).unwrap();
            let mut coupled = Shunted::new(
                body,
                vec![Reactive::SeriesRl {
                    electrodes: pair,
                    inductor: Arc::clone(&l),
                    conductor: Some(Arc::clone(&r)),
                }],
                vec![State::FluxLinkage(orientation * l.inductance_h() * i0)],
                vec![],
                controls(),
                cx,
            )
            .unwrap();
            let mut heat = 0.0;
            for _ in 0..64 {
                let candidate = coupled
                    .step_candidate(
                        time / 64.0,
                        &[Drive::Voltage(reference), Drive::Current(0.0)],
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap();
                let report = &candidate.step().reactive[0];
                match &report.load {
                    Reactive::SeriesRl {
                        inductor,
                        conductor: Some(conductor),
                        ..
                    } => {
                        assert!(Arc::ptr_eq(inductor, &l));
                        assert!(Arc::ptr_eq(conductor, &r));
                    }
                    _ => panic!("original material branch lost"),
                }
                near(
                    report.joule_heat_j,
                    time / 64.0 * r.resistance_ohm() * report.midpoint_current_a.powi(2),
                    1e-25,
                    1e-12,
                );
                near(report.energy_residual_j, 0.0, 1e-20, 0.0);
                heat += candidate.step().joule_heat_j;
                near(
                    candidate.model().body().state().total_energy_j()
                        + candidate.step().reactive_energy_j[1]
                        + heat,
                    initial,
                    1e-20,
                    1e-10,
                );
                coupled = candidate.into_parts().0;
            }
            let v = coupled.body().state().electrodes()[1].voltage_v - reference;
            let i = orientation * coordinate(coupled.states()[0]) / l.inductance_h();
            near(v, expected_v, 3e-5, 0.0);
            near(i, expected_i, 3e-9, 0.0);
            assert!(heat > 0.0);
            outcomes.push((v, i, heat));
        }
        near(outcomes[0].0, outcomes[1].0, 1e-11, 1e-11);
        near(outcomes[0].1, outcomes[1].1, 1e-15, 1e-11);
        near(outcomes[0].2, outcomes[1].2, 1e-21, 1e-10);
    });
}

#[test]
fn g1_reactive_current_changes_actual_mesh_motion_and_matches_scalar_elimination() {
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
        let charge = vec![0.0; n];
        let mut u = zero.clone();
        u[node][2] = 1e-8;
        let diagonal = |matrix: &fs_sparse::Csr| {
            let (cols, values) = matrix.row(dof);
            values[cols.iter().position(|&i| i == dof).unwrap()]
        };
        let (mass, stiffness) = (
            diagonal(&model.solid().mass),
            diagonal(&model.solid().stiffness),
        );
        let (cols, values) = model.coupling().row(dof);
        let coupling = cols
            .iter()
            .zip(values)
            .filter(|(i, _)| cube.top.contains(i))
            .map(|(_, v)| *v)
            .sum::<f64>();
        let ceramic_c = EPSILON[2] * LENGTH;
        let cap = capacitor(ceramic_c);
        let coil = inductor(1e-5);
        let resistor = resistor(100.0, [1, 0]).conductor;
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &u,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let open = body.clone();
        let v0 = body.state().electrodes()[1].voltage_v;
        let c = ceramic_c + cap.capacitance_f();
        let dt = 2e-8;
        let h = 0.5 * dt;
        let a = mass + h * h * stiffness;
        let conductance = h / (coil.inductance_h() + h * resistor.resistance_ohm());
        let vm = (c * v0 - h * h * coupling * stiffness * u[node][2] / a)
            / (c + h * conductance + h * h * coupling * coupling / a);
        let speed = (-h * stiffness * u[node][2] - h * coupling * vm) / a;
        let current = conductance * vm;
        let coupled = Shunted::new(
            body,
            vec![
                Reactive::SeriesRl {
                    electrodes: [1, 0],
                    inductor: Arc::clone(&coil),
                    conductor: Some(resistor),
                },
                Reactive::Capacitor {
                    electrodes: [1, 0],
                    capacitor: Arc::clone(&cap),
                },
            ],
            vec![
                State::FluxLinkage(0.0),
                State::Charge(cap.capacitance_f() * v0),
            ],
            vec![],
            controls(),
            cx,
        )
        .unwrap();
        let candidate = coupled
            .step_candidate(
                dt,
                &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let next = candidate.model().body().state();
        near(
            next.displacement_m()[node][2],
            u[node][2] + dt * speed,
            1e-19,
            1e-10,
        );
        near(next.velocity_m_s()[node][2], 2.0 * speed, 1e-11, 1e-10);
        near(next.electrodes()[1].voltage_v, 2.0 * vm - v0, 1e-10, 1e-10);
        near(
            coordinate(candidate.model().states()[0]),
            2.0 * coil.inductance_h() * current,
            1e-20,
            1e-10,
        );
        near(
            candidate.step().reactive[0].midpoint_current_a,
            current,
            1e-12,
            1e-10,
        );
        let open_next = open
            .step_candidate(
                dt,
                &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        assert!(
            (next.displacement_m()[node][2] - open_next.body().state().displacement_m()[node][2])
                .abs()
                > 1e-14
        );
        assert!(candidate.step().joule_heat_j > 0.0);
        near(candidate.step().energy_residual_j, 0.0, 1e-18, 0.0);
    });
}

#[test]
fn g1_voltage_sources_pay_for_both_external_storage_and_all_conductor_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let cap = capacitor(2e-11);
        let coil = inductor(0.01);
        let series = resistor(100.0, [1, 0]).conductor;
        let shunt = resistor(1e4, [1, 0]);
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Voltage(2.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let coupled = Shunted::new(
            body,
            vec![
                Reactive::SeriesRl {
                    electrodes: [1, 0],
                    inductor: Arc::clone(&coil),
                    conductor: Some(Arc::clone(&series)),
                },
                Reactive::Capacitor {
                    electrodes: [1, 0],
                    capacitor: Arc::clone(&cap),
                },
            ],
            vec![
                State::FluxLinkage(0.0),
                State::Charge(2.0 * cap.capacitance_f()),
            ],
            vec![shunt],
            controls(),
            cx,
        )
        .unwrap();
        let dt = 1e-6;
        let candidate = coupled
            .step_candidate(
                dt,
                &[Drive::Voltage(0.0), Drive::Voltage(4.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let i = 0.5 * dt * 3.0 / (coil.inductance_h() + 0.5 * dt * series.resistance_ohm());
        let expected_work =
            3.0 * (2.0 * (EPSILON[2] * LENGTH + cap.capacitance_f()) + dt * (i + 3.0 / 1e4));
        near(
            candidate.step().electrical_source_work_j,
            expected_work,
            1e-22,
            1e-11,
        );
        near(
            candidate.step().joule_heat_j,
            dt * (series.resistance_ohm() * i * i + 9.0 / 1e4),
            1e-22,
            1e-11,
        );
        near(candidate.step().energy_residual_j, 0.0, 1e-20, 0.0);
    });
}

#[test]
fn g4_reactive_refusal_and_cancellation_preserve_every_accepted_coordinate() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let cap = capacitor(1e-11);
        let coil = inductor(1e-3);
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let load = Reactive::Capacitor {
            electrodes: [1, 0],
            capacitor: Arc::clone(&cap),
        };
        for states in [
            vec![],
            vec![State::FluxLinkage(0.0)],
            vec![State::Charge(1e-9)],
            vec![State::Charge(f64::NAN)],
        ] {
            assert!(
                Shunted::new(
                    body.clone(),
                    vec![load.clone()],
                    states,
                    vec![],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        for pair in [[1, 1], [1, 2]] {
            assert!(
                Shunted::new(
                    body.clone(),
                    vec![Reactive::SeriesRl {
                        electrodes: pair,
                        inductor: Arc::clone(&coil),
                        conductor: None
                    }],
                    vec![State::FluxLinkage(0.0)],
                    vec![],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        let coupled = Shunted::new(
            body,
            vec![
                load,
                Reactive::SeriesRl {
                    electrodes: [1, 0],
                    inductor: coil,
                    conductor: None,
                },
            ],
            vec![State::Charge(0.0), State::FluxLinkage(0.0)],
            vec![],
            controls(),
            cx,
        )
        .unwrap();
        for (dt, drive) in [
            (0.0, Drive::Voltage(1.0)),
            (1e-8, Drive::Current(0.0)),
            (1e-8, Drive::Voltage(1e8)),
        ] {
            assert!(
                coupled
                    .step_candidate(
                        dt,
                        &[Drive::Voltage(0.0), drive],
                        &zero,
                        &charge,
                        controls(),
                        cx
                    )
                    .is_err()
            );
            assert_eq!(
                coupled.states(),
                &[State::Charge(0.0), State::FluxLinkage(0.0)]
            );
            assert_eq!(coupled.body().state().time_s(), 0.0);
        }
        let retry = coupled
            .step_candidate(
                1e-8,
                &[Drive::Voltage(0.0), Drive::Voltage(1.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        assert!(coordinate(retry.model().states()[0]) > 0.0);
        gate.request();
        assert!(matches!(
            coupled.step_candidate(
                1e-8,
                &[Drive::Voltage(0.0), Drive::Voltage(1.0)],
                &zero,
                &charge,
                controls(),
                cx
            ),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(coupled.body().state().total_energy_j(), 0.0);
        assert_eq!(
            coupled.states(),
            &[State::Charge(0.0), State::FluxLinkage(0.0)]
        );
    });
}
