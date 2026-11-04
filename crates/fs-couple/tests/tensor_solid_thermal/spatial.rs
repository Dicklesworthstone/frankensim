//! Original sparse motion supplies heat once, independently of thermal cuts.
use super::*;
use fs_couple::tensor_material::solid::thermal::spatial::{
    SpatialViscoLoading, step_thermal_candidate,
};
use fs_solid::linear3::visco::ViscoDynamicConfig;
use fs_solid::linear3::{TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState};

fn dynamic_controls() -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: 100,
        maximum_free_dofs: 128,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-13,
            maximum_iterations: 100,
        },
        absolute_energy_tolerance_j: 1e-13,
        relative_energy_tolerance: 1e-11,
    }
}
fn load(force: &[f64]) -> SpatialViscoLoading<'_> {
    SpatialViscoLoading {
        mean_force_n: force,
        support_velocity_mps: None,
        eigenstrain: None,
    }
}

#[test]
fn sparse_kelvin_midpoint_motion_deposits_analytic_heat_and_continues() {
    scope(|cx, _| {
        let mat = material(0.01, 0.);
        let solid = source(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &TETS,
            },
            TetViscoField::Uniform(&mat),
            &FIXED,
            cx,
        );
        let dt = 0.125;
        let op = solid.dynamic_operator(dt, dynamic_controls(), cx).unwrap();
        let mut velocity = [0.; 12];
        velocity[11] = 0.1;
        let initial = op.start(solid.state(cx).unwrap(), &velocity, cx).unwrap();
        let chart = curve(1.);
        let initial_heat = thermal(&chart, &[30], vec![], 16);
        let bindings = cells(&chart, &[30]);
        let mut motion = initial.clone();
        let mut heat = initial_heat.clone();
        let (mass, stiffness, viscosity) = (1. / 30., 20., 0.2);
        let (mut q, mut v, mut deposited) = (0., 0.1, 0.);
        for iteration in 1..=2 {
            let mean_v = (mass * v - 0.5 * dt * stiffness * q)
                / (mass + 0.5 * dt * viscosity + 0.25 * dt * dt * stiffness);
            q += dt * mean_v;
            v = 2. * mean_v - v;
            let loss = dt * viscosity * mean_v * mean_v;
            deposited += loss;
            let next = step_thermal_candidate(
                &motion,
                &heat,
                &(),
                load(&[0.; 12]),
                &bindings,
                controls(),
                cx,
            )
            .unwrap();
            near(
                next.solid().material_history().displacement_m()[11],
                q,
                1e-12,
            );
            near(next.solid().velocity_mps()[11], v, 1e-12);
            near(next.report().energy.material_heat_j, loss, 1e-12);
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300. + deposited / MASS,
                1e-9,
            );
            near(next.report().energy.external_thermal_heat_j, 0., 1e-12);
            near(next.report().energy.combined_residual_j, 0., 1e-10);
            assert_eq!(next.report().energy.thermal.substeps, 4);
            assert_eq!(next.solid().steps(), iteration);
            let (m, h, _, _) = next.into_parts();
            motion = m;
            heat = h;
        }
        assert_eq!(initial.steps(), 0);
        assert_eq!(initial_heat.time_s(), 0.);
    });
}

#[test]
fn original_regional_maxwell_heat_survives_thermal_events_and_all_transport_hooks() {
    scope(|cx, _| {
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|mut p| {
            p[0] += 2.;
            p
        }));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let fixed: Vec<_> = (0..24).filter(|&d| d != 11 && d != 23).collect();
        let materials = [material(0.01, 0.), material(0., 0.5)];
        let solid = source(
            SolidMeshRef::Tet4 {
                nodes_m: &nodes,
                elements: &tets,
            },
            TetViscoField::PerElement(&materials),
            &fixed,
            cx,
        );
        let duration = 0.0625;
        let op = solid
            .dynamic_operator(duration, dynamic_controls(), cx)
            .unwrap();
        let mut velocity = [0.; 24];
        velocity[11] = 0.1;
        velocity[23] = 0.1;
        let initial = op.start(solid.state(cx).unwrap(), &velocity, cx).unwrap();
        let chart = curve(1.);
        let switch = 0.0234375;
        let ambient = Boundary {
            ambient_k: 310.,
            htc_w_m2_k: 0.7,
            ..boundary(0.)
        };
        let reservoir = Boundary {
            reservoir: Some(Reservoir {
                temperature_k: 320.,
                area_m2: 0.2,
                conductance_w_k: 9.,
            }),
            ..boundary(0.)
        };
        let heat = thermal(
            &chart,
            &[30, 10],
            vec![
                TimedEvent {
                    time_s: 0.,
                    event: Event::Contact(ThermalLink {
                        id: 55,
                        a: CellId(30),
                        b: CellId(10),
                        area_m2: 0.1,
                        conductance_w_k: 0.,
                    }),
                },
                TimedEvent {
                    time_s: 0.,
                    event: Event::Boundary(CellId(30), ambient),
                },
                TimedEvent {
                    time_s: 0.,
                    event: Event::Boundary(CellId(10), reservoir),
                },
                TimedEvent {
                    time_s: switch,
                    event: Event::Boundary(
                        CellId(30),
                        Boundary {
                            source_w: 0.004,
                            ..ambient
                        },
                    ),
                },
            ],
            16,
        );
        let mut force = [0.; 24];
        force[11] = 0.01;
        force[23] = -0.02;
        let observer = OtherPhysics::default();
        let next = step_thermal_candidate(
            &initial,
            &heat,
            &observer,
            load(&force),
            &cells(&chart, &[30, 10]),
            controls(),
            cx,
        )
        .unwrap();
        let mut direct = initial.clone();
        let direct_step = direct.advance(&force, cx).unwrap();
        let heats: Vec<_> = direct_step
            .constitutive
            .element_steps
            .iter()
            .map(|e| e.step_dissipated_j)
            .collect();
        assert!((heats[0] - heats[1]).abs() > 1e-8);
        let obs = next.observer();
        let count = next.report().energy.thermal.substeps;
        assert_eq!(count, 5);
        assert_eq!(obs.events, 1);
        assert_eq!(
            (
                obs.sources,
                obs.conductances,
                obs.transfers,
                obs.limits,
                obs.steps
            ),
            (count, count, count, count, count)
        );
        assert_eq!(obs.intervals.len(), count);
        let mut temperatures = [300., 300.];
        let mut external = 0.;
        for &(start, dt) in &obs.intervals {
            let transfer = 0.5 * (temperatures[1] - temperatures[0]) + 0.0002;
            let a =
                0.2 * (310. - temperatures[0]) + 0.001 + if start >= switch { 0.004 } else { 0. };
            let b = 0.3 * (320. - temperatures[1]) + 0.002;
            external += dt * (a + b);
            temperatures[0] += dt * (heats[0] / duration + a + transfer) / MASS;
            temperatures[1] += dt * (heats[1] / duration + b - transfer) / MASS;
        }
        for (i, id) in [30, 10].into_iter().enumerate() {
            near(
                next.thermal().temperature_k(CellId(id)).unwrap(),
                temperatures[i],
                1e-9,
            );
            assert_eq!(next.report().energy.element_heat_j[i], heats[i]);
        }
        near(
            next.report().energy.external_thermal_heat_j,
            external,
            1e-11,
        );
        near(next.report().energy.combined_residual_j, 0., 1e-10);
        assert_eq!(
            next.solid().material_history().displacement_m(),
            direct.material_history().displacement_m()
        );
        assert_eq!(next.solid().velocity_mps(), direct.velocity_mps());
        assert_eq!(next.solid().steps(), 1); // Five thermal cuts never rebuild mechanics.
        assert_eq!(observer.steps, 0);
        assert_eq!(next.thermal().links()[&55].conductance_w_k, 0.);
        assert_eq!(next.thermal().cells()[&CellId(30)].boundary.htc_w_m2_k, 0.7);
    });
}

#[derive(Clone, Default)]
struct LateRefusal {
    steps: usize,
}
impl DepositionObserver for LateRefusal {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if thermal.time_s() > 0.04 {
            Err(DepositionError::InvalidInput("late synthetic refusal"))
        } else {
            Ok(())
        }
    }
}

#[test]
fn late_refusals_and_wrong_receiver_inventory_preserve_a_retriable_whole_state() {
    scope(|cx, gate| {
        let mat = material(0.01, 0.5);
        let solid = source(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &TETS,
            },
            TetViscoField::Uniform(&mat),
            &FIXED,
            cx,
        );
        let op = solid
            .dynamic_operator(0.125, dynamic_controls(), cx)
            .unwrap();
        let mut velocity = [0.; 12];
        velocity[11] = 0.1;
        let initial = op.start(solid.state(cx).unwrap(), &velocity, cx).unwrap();
        let chart = curve(1.);
        let heat = thermal(&chart, &[30], vec![], 16);
        let binding = cells(&chart, &[30]);
        let failed = LateRefusal::default();
        assert!(matches!(
            step_thermal_candidate(
                &initial,
                &heat,
                &failed,
                load(&[0.; 12]),
                &binding,
                controls(),
                cx
            ),
            Err(SolidThermalError::Observer(_))
        ));
        let mut limited = binding.clone();
        limited[0].temperature_interval_k = [300., 300.000000001];
        assert!(
            step_thermal_candidate(
                &initial,
                &heat,
                &(),
                load(&[0.; 12]),
                &limited,
                controls(),
                cx
            )
            .is_err()
        );
        let budget_heat = thermal(&chart, &[30], vec![], 1);
        assert!(matches!(
            step_thermal_candidate(
                &initial,
                &budget_heat,
                &(),
                load(&[0.; 12]),
                &binding,
                controls(),
                cx
            ),
            Err(SolidThermalError::Thermal(DepositionError::Budget(_)))
        ));
        let copy = chart.clone();
        assert!(
            step_thermal_candidate(
                &initial,
                &heat,
                &(),
                load(&[0.; 12]),
                &cells(&copy, &[30]),
                controls(),
                cx
            )
            .is_err()
        );
        let mut wrong = mat.clone();
        wrong.density_kg_m3 = 3.;
        let wrong_source = source(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &TETS,
            },
            TetViscoField::Uniform(&wrong),
            &FIXED,
            cx,
        );
        let wrong_op = wrong_source
            .dynamic_operator(0.125, dynamic_controls(), cx)
            .unwrap();
        let wrong_state = wrong_op
            .start(wrong_source.state(cx).unwrap(), &velocity, cx)
            .unwrap();
        assert!(
            step_thermal_candidate(
                &wrong_state,
                &heat,
                &(),
                load(&[0.; 12]),
                &binding,
                controls(),
                cx
            )
            .is_err()
        );
        let retry = step_thermal_candidate(
            &initial,
            &heat,
            &(),
            load(&[0.; 12]),
            &binding,
            controls(),
            cx,
        )
        .unwrap();
        let mut direct = initial.clone();
        let step = direct.advance(&[0.; 12], cx).unwrap();
        assert_eq!(
            retry.report().energy.material_heat_j,
            step.constitutive.step_dissipated_j
        );
        assert_eq!(
            retry.solid().material_history().displacement_m(),
            direct.material_history().displacement_m()
        );
        assert_eq!(initial.time_s(), 0.);
        assert_eq!(initial.steps(), 0);
        assert_eq!(initial.velocity_mps(), velocity);
        assert_eq!(heat.cells()[&CellId(30)].specific_h_j_kg, 0.);
        assert_eq!(failed.steps, 0);
        gate.request();
        assert!(matches!(
            step_thermal_candidate(
                &initial,
                &heat,
                &(),
                load(&[0.; 12]),
                &binding,
                controls(),
                cx
            ),
            Err(SolidThermalError::Cancelled)
        ));
    });
}

#[test]
fn moving_support_and_eigenstrain_work_remain_work_while_original_loss_heats() {
    scope(|cx, _| {
        let mat = material(0.01, 0.5);
        let solid = ViscoElements::new_with_full_reference_mass(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &TETS,
            },
            TetViscoField::Uniform(&mat),
            &FIXED,
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 10000,
            },
            cx,
        )
        .unwrap();
        let op = solid
            .dynamic_operator(0.0625, dynamic_controls(), cx)
            .unwrap();
        let initial = op.start(solid.state(cx).unwrap(), &[0.; 12], cx).unwrap();
        let chart = curve(1.);
        let heat = thermal(&chart, &[30], vec![], 16);
        let strain = TetThermalStrainState::try_new(
            320.,
            300.,
            [0.001, 0., 0., 0., 0., 0.],
            mat.material_state_identity,
            ContentHash([89; 32]),
            0.05,
        )
        .unwrap();
        let eigenstrain = TetThermalStrainField::Uniform(&strain).into();
        let mut support = [0.; 12];
        support[3] = 0.01;
        for (moving, expanding) in [(true, false), (false, true), (true, true)] {
            let loading = SpatialViscoLoading {
                mean_force_n: &[0.; 12],
                support_velocity_mps: moving.then_some(support.as_slice()),
                eigenstrain: expanding.then_some(eigenstrain),
            };
            let next = step_thermal_candidate(
                &initial,
                &heat,
                &(),
                loading,
                &cells(&chart, &[30]),
                controls(),
                cx,
            )
            .unwrap();
            let mut direct = initial.clone();
            let step = match (moving, expanding) {
                (true, true) => direct.advance_with_support_velocity_and_thermal_strain(
                    &[0.; 12],
                    &support,
                    eigenstrain,
                    cx,
                ),
                (true, false) => direct.advance_with_support_velocity(&[0.; 12], &support, cx),
                _ => direct.advance_with_thermal_strain(&[0.; 12], eigenstrain, cx),
            }
            .unwrap();
            let energy = &next.report().energy;
            let work = step.external_work_j
                + step.support_work_j
                + step.constitutive.step_eigenstrain_work_j;
            near(energy.mechanical_work_j, work, 1e-12);
            near(
                energy.material_heat_j,
                step.constitutive.step_dissipated_j,
                1e-12,
            );
            near(
                energy.thermal.stored_energy_change_j,
                energy.material_heat_j,
                1e-10,
            );
            near(energy.combined_residual_j, 0., 1e-10);
            assert!(energy.material_heat_j > 0.);
            if moving {
                assert_ne!(step.support_work_j, 0.);
            }
            if expanding {
                assert_ne!(step.constitutive.step_eigenstrain_work_j, 0.);
            }
        }
    });
}
