//! Caloric and elastic source data govern the same damped inertial trajectory.
use super::*;
use fs_solid::loads3::transient::TransientConfig3;

#[test]
fn source_data_drives_damped_heating_with_reference_mass_and_endpoint_acceleration() {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let mut initial_q = [[0.0; 3]; 4]; initial_q[1][0] = 0.001;
        let mut initial_v = [[0.0; 3]; 4]; initial_v[1][0] = 0.002;
        let controls = TransientConfig3 {
            mass_damping_per_s: 0.4, stiffness_damping_s: 0.06,
            ..TransientConfig3::default()
        };
        let mut observations = Vec::new();
        for (young_end, cp) in [(36.0, 10.0), (90.0, 10.0), (36.0, 20.0)] {
            // Only source data changes. Geometry, heat input, reference mass,
            // excitation and the declared Rayleigh coefficients are identical.
            let (card, caloric) = material_with_modulus_endpoint(
                60.0, young_end, 1e-5, cp, 0.0, false);
            let thermal = initial(caloric.curve(), &[(30, 20_000.0)], 1.0);
            let context = point(300.0);
            let model = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
                vec![source(30, &card, caloric.curve(), &context)],
                TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx).unwrap();
            let loads = model.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let heated = thermal.advance_to(0.01, thermal_config(), || Ok(())).unwrap().0;
            let heated_again = heated.advance_to(0.02, thermal_config(), || Ok(())).unwrap().0;
            let states = [model.resolve(&thermal, cx).unwrap(),
                model.resolve(&heated, cx).unwrap(), model.resolve(&heated_again, cx).unwrap()];
            let problems: Vec<_> = states.iter().map(|s| s.problem(&fixed)).collect();
            let mut run = loads.begin_transient(&problems[0], Some(states[0].thermal_strains()),
                0.0, &initial_q, &initial_v, controls, cx).unwrap();
            let (mut q, mut v, mut old_young, mut old_strain) = (0.001, 0.002, 60.0, 0.0);
            let mut dissipated = 0.0;
            for index in 1..states.len() {
                let rise = 200.0 * index as f64 / cp;
                let young = 60.0 + (young_end - 60.0) * rise / 200.0;
                let strain = 1e-5 * (rise + 0.005 * rise * rise);
                let (h, m) = (0.01, 0.1);
                let (k0, k1) = (old_young / 6.0, young / 6.0);
                let cbar = 0.4 * m + 0.06 * (k0 + k1) / 2.0;
                let force_mid = (k0 * old_strain + k1 * strain) / 2.0;
                // Independent scalar momentum/kinematic equations, not a
                // second invocation of the native time integrator.
                let v1 = ((m / h - cbar / 2.0) * v - k0 * q / 2.0 + force_mid
                    - k1 * (q + h * v / 2.0) / 2.0)
                    / (m / h + cbar / 2.0 + h * k1 / 4.0);
                let q1 = q + h * (v + v1) / 2.0;
                let a1 = (k1 * (strain - q1) - (0.4 * m + 0.06 * k1) * v1) / m;
                let saved = run.state().clone();
                let (candidate, report) = states[index].advance_transient_candidate(
                    &states[index - 1], &run, &problems[index], &loads, cx).unwrap();
                near(candidate.state().displacement[0], q1, 1e-14);
                near(candidate.state().velocity[0], v1, 1e-12);
                near(candidate.state().acceleration[0], a1, 1e-10);
                let loss = h * cbar * ((v + v1) / 2.0).powi(2);
                near(report.native.damping_loss_j, loss, 1e-14);
                assert!(report.native.damping_loss_j >= 0.0);
                dissipated += report.native.damping_loss_j;
                // Isotropic expansion also stores energy in both restrained
                // transverse directions. An axial-only storage shortcut fails.
                near(report.kinetic_plus_elastic_after_j,
                    0.5 * m * v1 * v1 + young / 12.0
                        * ((q1 - strain).powi(2) + 2.0 * strain * strain), 1e-13);
                near(report.stiffness_material_work_j, (young - old_young) / 12.0
                    * ((q - old_strain) * (q1 - strain) + 2.0 * old_strain * strain), 1e-13);
                near(report.kinetic_plus_elastic_after_j - report.kinetic_plus_elastic_before_j
                    + report.native.damping_loss_j,
                    report.thermal_material_work_j + report.stiffness_material_work_j, 1e-13);
                let response = candidate.recover_response(config().mesh_frame, 0.1, cx).unwrap();
                near(response.stress.elements()[0].stress_mandel_pa()[0], young * (q1 - strain), 1e-11);
                near(states[index].materials()[0].density_kg_m3, 6.0, 1e-12);
                near(states[index].samples()[0].reference_mass_kg, 1.0, 1e-12);
                assert_eq!(states[index].samples()[0].elastic.property("young_modulus").unwrap()
                    .answer().receipt.selected, states[0].samples()[0].elastic.property("young_modulus")
                    .unwrap().answer().receipt.selected);
                assert_eq!(run.state(), &saved);
                assert_eq!(candidate.state().time_s.to_bits(), states[index].time_s().to_bits());
                (q, v, old_young, old_strain) = (q1, v1, young, strain);
                run = candidate;
            }
            assert_eq!(run.state().step, 2);
            assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg, 0.0);
            observations.push([q, v, dissipated]);
        }
        // Both a different elastic law and a different caloric law propagate
        // into motion and dissipation, without modifying a physics coefficient.
        for other in &observations[1..] {
            for axis in 0..3 { assert_ne!(observations[0][axis], other[axis]); }
        }
    });
}
