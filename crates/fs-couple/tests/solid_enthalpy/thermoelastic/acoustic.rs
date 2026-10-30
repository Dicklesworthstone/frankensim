//! Native heat -> source potential -> isentropic stiffness -> tet/wave output.
use super::*;
use fs_solid::loads3::transient::TransientConfig3;

#[test]
fn g1_accepted_heating_and_capacity_swaps_reach_native_wave_and_mass_stiffness_operators() {
    with_cx(|cx| {
        let mut cold_speeds = Vec::new();
        for cp in [0.2, 0.4] {
            let alpha = 0.0003;
            let (card, curve) = fixture(1000.0, 2.0, alpha, cp, false);
            let heat = thermal(&[&curve], &[cp], &[2.0], 0.0);
            let context = context(300.0);
            let source = SolidEnthalpy::try_new(
                &heat,
                &NODES,
                &TETS,
                vec![source(30, &card, &curve, &context)],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let body = ThermoelasticSolid::try_new(
                &source,
                &heat,
                &loads,
                &FREE,
                &[[280.0, 360.0]],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let cold = body.adiabatic_state(cx).unwrap();
            let next = body
                .advance(&heat, &(), 0.1, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let warm = next.solid.adiabatic_state(cx).unwrap();
            let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
            for (snapshot, temperature, response_identity) in [
                (&cold, 300.0, body.equilibrium().solution.stress.identity()),
                (
                    &warm,
                    300.0 + 0.2 / cp,
                    next.solid.equilibrium().solution.stress.identity(),
                ),
            ] {
                assert_eq!(snapshot.cells(), &[IDS[0]]);
                assert_eq!(snapshot.mechanical_response_identity(), response_identity);
                assert_eq!(snapshot.mesh_frame(), mechanical().mesh_frame);
                let response = &snapshot.responses()[0];
                near(response.source().temperature_k(), temperature, 2e-8);
                near(
                    response.input().caloric_heat_capacity_j_m3_k,
                    6.0 * cp,
                    1e-14,
                );
                near(response.input().reference_density_kg_m3, 6.0, 1e-14);
                let young = 1000.0 + 2.0 * (temperature - 300.0);
                // In this free, nu=0 reference state, C_epsilon =
                // rho0*Cp - 3*T*E*alpha^2 and the longitudinal modulus follows
                // from the independent isotropic K_S/K_T heat-capacity identity.
                let c_epsilon = 6.0 * cp - 3.0 * temperature * young * alpha * alpha;
                let bulk_adiabatic = young / 3.0 * (6.0 * cp) / c_epsilon;
                let longitudinal = bulk_adiabatic + 2.0 * young / 3.0;
                let waves = response
                    .acoustic()
                    .directional_waves([1.0, 0.0, 0.0])
                    .unwrap();
                near(
                    waves.modes[2].phase_speed_m_s,
                    (longitudinal / 6.0).sqrt(),
                    2e-8,
                );
                near(waves.modes[0].phase_speed_m_s, (young / 12.0).sqrt(), 2e-8);
                let native = snapshot.problem(&fixed).assemble(cx).unwrap();
                assert_eq!(native.free_dofs, vec![3]);
                near(native.total_mass_kg, 1.0, 1e-14);
                near(native.mass.get(0, 0), 0.1, 1e-14);
                near(native.stiffness.get(0, 0), VOLUME * longitudinal, 2e-7);
                // The one-DOF native natural frequency uses the same tangent,
                // not an independently tuned acoustic coefficient.
                near(
                    (native.stiffness.get(0, 0) / native.mass.get(0, 0)).sqrt(),
                    (longitudinal / 0.6).sqrt(),
                    2e-8,
                );
                let material = &snapshot.materials()[0];
                assert_eq!(
                    material.material_state_identity,
                    response.acoustic().material_state_identity().0
                );
                assert_eq!(
                    material.stiffness_mandel_pa(),
                    response.acoustic().stiffness_mandel_pa()
                );
            }
            assert_eq!(cold.time_s(), 0.0);
            assert_eq!(warm.time_s(), 0.1);
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(heat.time_s(), 0.0);
            assert_ne!(
                cold.materials()[0].material_state_identity,
                warm.materials()[0].material_state_identity
            );
            cold_speeds.push(
                cold.responses()[0]
                    .acoustic()
                    .directional_waves([1.0, 0.0, 0.0])
                    .unwrap()
                    .modes[2]
                    .phase_speed_m_s,
            );
        }
        assert!(cold_speeds[0] > cold_speeds[1]);
    });
}

#[test]
fn g3_accepted_inertial_strain_and_heat_feed_the_same_adiabatic_snapshot() {
    with_cx(|cx| {
        let (card, curve) = fixture(1000.0, 2.0, 0.0003, 0.2, false);
        let heat = thermal(&[&curve], &[0.2], &[2.0], 0.0);
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(
            &heat,
            &NODES,
            &TETS,
            vec![source(30, &card, &curve, &context)],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
            .unwrap();
        let body = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &loads,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            solve(),
            cx,
        )
        .unwrap();
        let velocity = [[0.0; 3], [0.05, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let controls = TransientConfig3 {
            maximum_steps: 100,
            maximum_displacement_m: 0.05,
            absolute_force_tolerance_n: 1e-11,
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-8,
            ..TransientConfig3::default()
        };
        let moving = body
            .begin_transient(&loads, &velocity, controls, cx)
            .unwrap();
        let next = moving
            .advance(&heat, &(), 0.02, &loads, thermal_controls(), solve(), cx)
            .unwrap();
        let acoustic = next.solid.adiabatic_state(cx).unwrap();
        let motion = next.solid.nodal_motion(cx).unwrap();
        let sample = &next.solid.state().samples()[0];
        let response = &acoustic.responses()[0];
        near(
            response.input().total_strain_mandel[0],
            motion[0][1][0],
            1e-14,
        );
        for value in &response.input().total_strain_mandel[1..] {
            near(*value, 0.0, 1e-14);
        }
        assert_eq!(
            response.input().caloric_state_identity,
            sample.phase.identity()
        );
        assert_eq!(response.source().elastic, sample.elastic);
        assert_eq!(
            acoustic.mechanical_response_identity(),
            next.solid.response().stress.identity()
        );
        assert_eq!(acoustic.time_s(), next.solid.time_s());
        assert_eq!(acoustic.materials()[0].density_kg_m3, 6.0);
        let repeated = next.solid.adiabatic_state(cx).unwrap();
        assert_eq!(acoustic.materials(), repeated.materials());
        assert_eq!(moving.step(), 0);
        assert_eq!(next.solid.step(), 1);
    });
}

#[test]
fn g4_caloric_kinks_and_one_sided_source_domains_cannot_masquerade_as_wave_tangents() {
    with_cx(|cx| {
        let (card, smooth) = fixture(1000.0, 0.0, 0.0003, 0.2, false);
        let kink = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    temperature_k: 280.0,
                    specific_enthalpy_j_kg: 0.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 6.0,
                },
                EnthalpyPhaseKnot {
                    temperature_k: 300.0,
                    specific_enthalpy_j_kg: 4.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 6.0,
                },
                EnthalpyPhaseKnot {
                    temperature_k: 500.0,
                    specific_enthalpy_j_kg: 64.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 6.0,
                },
            ],
        )
        .unwrap();
        let context = context(300.0);
        for (curve, range) in [(&kink, [280.0, 360.0]), (&smooth, [300.0, 360.0])] {
            let heat = thermal(&[curve], &[0.2], &[0.0], 0.0);
            let source = SolidEnthalpy::try_new(
                &heat,
                &NODES,
                &TETS,
                vec![source(30, &card, curve, &context)],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let body = ThermoelasticSolid::try_new(
                &source,
                &heat,
                &loads,
                &FREE,
                &[range],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let before = body.equilibrium().solution.displacement.clone();
            assert!(body.adiabatic_state(cx).is_err());
            assert_eq!(body.equilibrium().solution.displacement, before);
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(heat.time_s(), 0.0);
        }
    });
}
