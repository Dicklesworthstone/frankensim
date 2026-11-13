//! G1/G3/G4: true load-release inertia, reciprocal heating and accepted history.
use super::*;
use fs_couple::cut_thermoelastic::dynamic::CutThermoelasticDynamicsSettings3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_solid::cutfront::dynamics3::CutThermalMidpointSettings3;
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::ImplicitSolveConfig;
use fs_time::galpha::initialization::InitialSolveConfig;

fn controls(n: usize) -> CutThermoelasticDynamicsSettings3 {
    CutThermoelasticDynamicsSettings3 {
        initial: InitialSolveConfig { restart: n, max_cycles: 8, tolerance: 1e-12 },
        midpoint: CutThermalMidpointSettings3 {
            solve: ImplicitSolveConfig {
                newton: NewtonKrylovConfig {
                    absolute_tolerance: 1e-11, relative_tolerance: 1e-11,
                    linear_restart: n, max_linear_cycles: 8, forcing_maximum: 1e-3,
                    ..Default::default()
                },
                max_newton_iterations: 8,
            },
            maximum_steps: 16, maximum_workspace_components: 2_000_000,
            maximum_force_residual_n: 2e-8,
        },
    }
}
fn free_dimension(model: &Model) -> usize {
    3 * native!(model, fixed).iter().filter(|&&fixed| !fixed).count()
}
fn thermal_step() -> StepConfig {
    StepConfig { max_dt_s: 0.005, energy_tolerance_j: step().energy_tolerance_j }
}

#[test]
fn released_loaded_raw_and_hanging_solids_keep_real_inertia_during_reciprocal_heating() {
    for adaptive in [false, true] {
        with_cx(|cx, _| {
            let card = card(); let law = law(&card); let curve = caloric(&card);
            let context = point(300.0);
            let mut model = build(adaptive, &source(&law, 300.0));
            let mut thermal = timeline(&model, curve.curve(), false);
            let loading = forces(&model);
            let zero = vec![0.0; model.operator().n()];
            let statics = CutThermoelasticSolid3::try_new(model.source(), &thermal,
                cells(&model, &card, &context, curve.curve(), &law), &loading,
                settings(model.operator().n()), config(), cx).unwrap();
            let equilibrium_u = statics.equilibrium().solution.coefficients().to_vec();
            // Release at the SAME initial instant: q0 is retained, F0 is zero,
            // and native initial acceleration must balance the remaining stress.
            let mut solid = statics.begin_transient(model.source(), &zero, &zero,
                controls(free_dimension(&model)), cx).unwrap();
            assert_eq!(solid.state().displacement, equilibrium_u);
            assert!(solid.native_acceleration_m_s2().iter().any(|a| a.abs() > 1e-5));
            assert_eq!(solid.kinetic_energy_j(), 0.0);
            assert_eq!(solid.motion().steps, 0);
            let occupancy = vec![1.0; native!(&model, cells)];
            let original_mass = solid.state().mass.diagonal(&occupancy, proceed).unwrap();
            let original_inventory = thermal.mass_kg();
            let mut observer = Observe::default();

            for endpoint in [0.025, 0.05] {
                let dt = endpoint - solid.time_s();
                let old_q = solid.state().displacement.clone();
                let old_v = solid.native_velocity_m_s().to_vec();
                let old_a = solid.native_acceleration_m_s2().to_vec();
                let old_k = solid.kinetic_energy_j();
                let old_heat = thermal.energy_j();
                let one_way = thermal.advance_to(endpoint, thermal_step(), || Ok(())).unwrap().0;
                let previous_steps = observer.steps;
                let authored: f64 = thermal.cells().values().map(|c| c.boundary.source_w * dt).sum();
                let report = solid.advance(model.mutable(), &mut thermal, &mut observer,
                    endpoint, &zero, thermal_step(), config(), cx).unwrap();
                assert!(report.coupled.iterations > 1, "the physical reversible source must be iterated");
                assert_eq!(solid.time_s(), thermal.time_s());
                assert_eq!(solid.time_s(), endpoint);
                assert_eq!(solid.motion().history.len(), solid.step());
                assert_eq!(observer.steps, previous_steps + report.coupled.thermal_audit.substeps);
                assert!(solid.native_velocity_m_s().iter().any(|v| v.abs() > 1e-6));
                assert!(thermal.cells().iter().any(|(id, c)|
                    (c.specific_h_j_kg - one_way.cells()[id].specific_h_j_kg).abs() > 1e-8));
                assert_eq!(report.coupled.external_work_j, 0.0);
                near(report.coupled.external_heat_j, authored, 2e-9);
                near(thermal.energy_j() - old_heat + solid.kinetic_energy_j() - old_k
                    + report.coupled.internal_energy_change_j + report.coupled.numerical_work_j,
                    authored, 3e-8);
                near(report.kinetic_energy_change_j + report.coupled.native_work_j,
                    report.coupled.external_work_j, 2e-8);
                near(report.coupled.physical_bulk_work_j - report.coupled.internal_energy_change_j,
                    report.coupled.reversible_heat_j, 1e-12);
                assert!(report.coupled.combined_residual_j.abs() <= config().energy_tolerance_j);

                // Independently read the original consistent mass and actual
                // endpoint resistance. This is dynamic force balance, not a
                // residual which was incorrectly forced to static equilibrium.
                let mut inertia = vec![0.0; zero.len()];
                solid.state().mass.apply_into(&occupancy,
                    solid.native_acceleration_m_s2(), &mut inertia, proceed).unwrap();
                for i in 0..zero.len() {
                    if !native!(&model, fixed)[i / 3] {
                        near(inertia[i] + solid.state().resistance[i], 0.0, 5e-8);
                    }
                    near(solid.state().displacement[i] - old_q[i],
                        0.5 * dt * (old_v[i] + solid.native_velocity_m_s()[i]), 2e-11);
                    near(solid.native_velocity_m_s()[i] - old_v[i],
                        0.5 * dt * (old_a[i] + solid.native_acceleration_m_s2()[i]), 2e-11);
                }
                assert_eq!(solid.state().mass.diagonal(&occupancy, proceed).unwrap(), original_mass);
                assert_eq!(thermal.mass_kg(), original_inventory);
                native!(&model, validate_material_mass, &solid.state().mass).unwrap();
                native!(&model, validate_eigenstrain, &solid.state().eigenstrain).unwrap();
                let actual_g = solid.state().prestress.as_ref().unwrap();
                let recovered = native!(&model, prepare_eigenstrain_prestress,
                    &solid.state().eigenstrain, &solid.state().displacement,
                    0.08, POINTS, POINTS, proceed).unwrap();
                let direction: Vec<_> = (0..zero.len()).map(|i| (i as f64).sin()).collect();
                let mut a = zero.clone(); let mut b = zero.clone();
                actual_g.apply_into(&direction, &mut a, proceed).unwrap();
                recovered.apply_into(&direction, &mut b, proceed).unwrap();
                assert_eq!(a, b);

                // A fresh static solve on the SAME accepted hot coefficients
                // gives a different displacement; the transient did not reset
                // its real motion to instantaneous thermal equilibrium.
                let equilibrium = native!(&model, solve_thermoelastic_equilibrium,
                    ReferenceLoad3::default(), Some(zero.as_slice()), None,
                    settings(model.operator().n()), |_| proceed()).unwrap();
                assert!(solid.state().displacement.iter().zip(equilibrium.solution.coefficients())
                    .any(|(a,b)| (a-b).abs() > 1e-7));
            }
            assert_eq!(solid.step(), 2);
            near(observer.duration_s, 0.05, 1e-14);
        });
    }
}

#[test]
fn rejected_heat_iteration_late_observer_and_cancel_keep_the_original_motion_prefix() {
    with_cx(|cx, gate| {
        let card = card(); let law = law(&card); let curve = caloric(&card);
        let context = point(300.0);
        let mut model = build(true, &source(&law, 300.0));
        let mut thermal = timeline(&model, curve.curve(), false);
        let loads = forces(&model); let zero = vec![0.0; model.operator().n()];
        let statics = CutThermoelasticSolid3::try_new(model.source(), &thermal,
            cells(&model, &card, &context, curve.curve(), &law), &loads,
            settings(model.operator().n()), config(), cx).unwrap();
        let mut dyn_controls = controls(free_dimension(&model));
        dyn_controls.midpoint.maximum_steps = 1;
        let mut solid = statics.begin_transient(model.source(), &zero, &zero,
            dyn_controls, cx).unwrap();
        let before_motion = solid.motion().clone();
        let before_native = action(&model);
        let before_heat = thermal.cells().clone();
        let mut observer = Observe::default();
        let one_trial = CutThermoelasticConfig3 {maximum_iterations: 1, heat_tolerance_j: 1e-25, ..config()};
        assert!(matches!(solid.advance(model.mutable(), &mut thermal, &mut observer,
            0.025, &zero, thermal_step(), one_trial, cx),
            Err(CutThermoelasticError3::NotConverged {..})));
        assert_eq!(solid.motion(), &before_motion);
        assert_eq!(thermal.cells(), &before_heat);
        assert_eq!(action(&model), before_native);
        assert_eq!(observer, Observe::default());
        native!(&model, validate_material_mass, &solid.state().mass).unwrap();
        native!(&model, validate_eigenstrain, &solid.state().eigenstrain).unwrap();

        // A converged physical heat/motion trial can still fail its separate
        // numerical-work allowance. This late gate cannot publish motion or
        // reinterpret stabilization work as a caloric source.
        let exact_numerical = CutThermoelasticConfig3 {maximum_numerical_work_j: 0.0, ..config()};
        assert!(matches!(solid.advance(model.mutable(), &mut thermal, &mut observer,
            0.025, &zero, thermal_step(), exact_numerical, cx),
            Err(CutThermoelasticError3::NumericalWork {..})));
        assert_eq!(solid.motion(), &before_motion);
        assert_eq!(thermal.cells(), &before_heat);
        assert_eq!(action(&model), before_native);
        assert_eq!(observer, Observe::default());

        let mut late = Observe {fail_after_s: Some(0.024), ..Default::default()};
        let old_late = late.clone();
        assert!(matches!(solid.advance(model.mutable(), &mut thermal, &mut late,
            0.025, &zero, thermal_step(), config(), cx), Err(CutThermoelasticError3::Observer(_))));
        assert_eq!(late, old_late); assert_eq!(solid.motion(), &before_motion);
        assert_eq!(thermal.cells(), &before_heat); assert_eq!(action(&model), before_native);

        solid.advance(model.mutable(), &mut thermal, &mut observer,
            0.025, &zero, thermal_step(), config(), cx).unwrap();
        let prefix_motion = solid.motion().clone();
        let prefix_heat = thermal.cells().clone(); let prefix_observer = observer.clone();
        let prefix_native = action(&model);
        assert_eq!(solid.step(), 1);
        assert!(solid.advance(model.mutable(), &mut thermal, &mut observer,
            0.05, &zero, thermal_step(), config(), cx).is_err());
        assert_eq!(solid.motion(), &prefix_motion);
        assert_eq!(thermal.cells(), &prefix_heat);
        assert_eq!(observer, prefix_observer); assert_eq!(action(&model), prefix_native);
        gate.request();
        assert!(matches!(solid.advance(model.mutable(), &mut thermal, &mut observer,
            0.05, &zero, thermal_step(), config(), cx), Err(CutThermoelasticError3::Cancelled)));
        assert_eq!(solid.motion(), &prefix_motion);
        assert_eq!(thermal.cells(), &prefix_heat);
        assert_eq!(observer, prefix_observer); assert_eq!(action(&model), prefix_native);
    });
}

#[test]
fn prepared_source_conduction_changes_the_actual_inertial_and_caloric_trajectory_together() {
    use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
    use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
    use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
    with_cx(|cx, _| {
        let card = card(); let law = law(&card); let curve = caloric(&card); let context = point(300.0);
        let cold = source(&law, 300.0);
        let mut model = build(false, &cold); let mut bare_model = build(false, &cold);
        let mut thermal = timeline(&model, curve.curve(), true); let mut bare_heat = thermal.clone();
        let loads = forces(&model); let zero = vec![0.0; model.operator().n()];
        let statics = CutThermoelasticSolid3::try_new(model.source(), &thermal,
            cells(&model, &card, &context, curve.curve(), &law), &loads,
            settings(model.operator().n()), config(), cx).unwrap();
        let bare_statics = CutThermoelasticSolid3::try_new(bare_model.source(), &bare_heat,
            cells(&bare_model, &card, &context, curve.curve(), &law), &loads,
            settings(bare_model.operator().n()), config(), cx).unwrap();
        let mut solid = statics.begin_transient(model.source(), &zero, &zero,
            controls(free_dimension(&model)), cx).unwrap();
        let mut bare = bare_statics.begin_transient(bare_model.source(), &zero, &zero,
            controls(free_dimension(&bare_model)), cx).unwrap();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let material_cells: Vec<_> = (0..2).map(|i| MaterialCell {
            cell: CellId(i+10), material_id: MaterialId(0), curve: curve.curve(),
            point: binding(&card, &context),
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        }).collect();
        let links = [ConductionLink {link: 7, distance_a_m: 0.25, distance_b_m: 0.25}];
        let bulk = MaterialNetwork {cells: &material_cells, resistors: &[], links: &links};
        let electrical = ElectricalContactNetwork {bulk: &bulk, contacts: &[], electrical_bases: &[]};
        let network = MaterialHeatNetwork {electrical: &electrical, thermal_contacts: &[], radiation: None};
        let prepared = network.prepare(cx, &thermal).unwrap();
        near(prepared.resolve(cx, &thermal).unwrap().bulk.conductances_w_k[&7],
            2.0*0.73/0.5, 1e-13);
        let mut observer = Observe::default(); let mut bare_observer = Observe::default();
        for endpoint in [0.025, 0.05, 0.075, 0.1] {
            bare.advance(bare_model.mutable(), &mut bare_heat, &mut bare_observer, endpoint,
                &zero, thermal_step(), config(), cx).unwrap();
            let report = solid.advance_heat_network(model.mutable(), &mut thermal, &prepared,
                &mut observer, endpoint, &zero, thermal_step(), config(), cx).unwrap();
            assert!(report.coupled.combined_residual_j.abs() <= config().energy_tolerance_j);
            assert_eq!(solid.time_s(), thermal.time_s());
            assert_eq!(report.coupled.closed_radiation_net_j, 0.0);
        }
        assert!(thermal.temperature_k(CellId(10)).unwrap() >
            bare_heat.temperature_k(CellId(10)).unwrap() + 1e-6);
        assert!(thermal.temperature_k(CellId(11)).unwrap() <
            bare_heat.temperature_k(CellId(11)).unwrap() - 1e-6);
        assert!(solid.state().displacement.iter().zip(&bare.state().displacement)
            .any(|(a,b)| (a-b).abs() > 1e-9));
        assert_eq!(solid.step(), 4); assert_eq!(bare.step(), 4);
        assert_eq!(thermal.links()[&7].conductance_w_k, 0.0);
        near(observer.duration_s, 0.1, 1e-14);
        native!(&model, validate_material_mass, &solid.state().mass).unwrap();
    });
}
