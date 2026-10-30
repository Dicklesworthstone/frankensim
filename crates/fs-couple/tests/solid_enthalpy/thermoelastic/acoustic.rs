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
        let combined = next.solid.adiabatic_prestressed_state(cx).unwrap();
        assert_eq!(combined.tangent().stress(), &next.solid.response().stress);
        let only_material = acoustic.problem(&fixed).assemble(cx).unwrap();
        let sigma_xx = next.solid.response().stress.elements()[0].stress_mandel_pa()[0];
        near(
            combined.tangent().geometric_stiffness().get(0, 0),
            sigma_xx / 6.0,
            1e-10,
        );
        near(
            combined.tangent().pencil().stiffness.get(0, 0),
            only_material.stiffness.get(0, 0) + sigma_xx / 6.0,
            1e-9,
        );
        assert_eq!(
            combined.tangent().pencil().mass.to_dense(),
            only_material.mass.to_dense()
        );
        // Inertia means the accepted stress is NOT the instantaneous static
        // force/area value. The new pencil must retain that distinction.
        assert!((sigma_xx - 6.0 * 0.8).abs() > 0.01);
        assert_eq!(next.solid.nodal_motion(cx).unwrap(), motion);
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

// Independent local internal energy on an isentrope, using an authored nu=0,
// constant-E/alpha potential and the reference caloric entropy. The virtual
// Green-strain increment contributes the native initial-stress approximation.
fn independent_adiabatic_virtual_energy(
    scale: f64,
    gradient: [[f64; 3]; 3],
    base_temperature: f64,
    young: f64,
    alpha: f64,
    caloric_capacity: f64,
) -> f64 {
    let free = alpha * (base_temperature - 300.0);
    let strain: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let base = if i == j && i > 0 { free } else { 0.0 };
            base + 0.5 * scale * (gradient[i][j] + gradient[j][i])
                + 0.5 * scale * scale * (0..3).map(|k| gradient[k][i] * gradient[k][j]).sum::<f64>()
        })
    });
    let entropy = |temperature: f64| {
        let trace = (0..3).map(|i| strain[i][i]).sum::<f64>() - 3.0 * alpha * (temperature - 300.0);
        caloric_capacity * (temperature / 300.0).ln() + young * alpha * trace
    };
    let target = caloric_capacity * (base_temperature / 300.0).ln() - young * alpha * free;
    let mut low = base_temperature - 1.0;
    let mut high = base_temperature + 1.0;
    assert!(entropy(low) < target && entropy(high) > target);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if entropy(mid) < target {
            low = mid;
        } else {
            high = mid;
        }
    }
    let temperature = 0.5 * (low + high);
    let thermal = alpha * (temperature - 300.0);
    let elastic: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| strain[i][j] - if i == j { thermal } else { 0.0 })
    });
    let energy = 0.5 * young * elastic.iter().flatten().map(|v| v * v).sum::<f64>();
    let elastic_entropy = young * alpha * (0..3).map(|i| elastic[i][i]).sum::<f64>();
    VOLUME * (caloric_capacity * (temperature - 300.0) + energy + temperature * elastic_entropy)
}

#[test]
fn g1_restraint_heating_combines_adiabatic_stiffness_and_signed_prestress_in_native_modes() {
    with_cx(|cx| {
        let fixed = [0, 1, 2, 3, 6, 8, 9];
        let young = 1000.0;
        let alpha = 0.0003;
        let cp = 0.2;
        let mut eigenvalues = Vec::new();
        for power in [-2.0, 0.0, 2.0] {
            let (card, curve) = fixture(young, 0.0, alpha, cp, false);
            let heat = thermal(&[&curve], &[cp], &[power], 0.0);
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
                &fixed,
                &[[280.0, 360.0]],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let next = body
                .advance(&heat, &(), 0.1, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let response = next.solid.adiabatic_prestressed_state(cx).unwrap();
            let material = response.material_state();
            let tangent = response.tangent();
            let pencil = tangent.pencil();
            let temperature = next.thermal.temperature_k(IDS[0]).unwrap();
            let stress = -young * alpha * (temperature - 300.0);
            assert_eq!(tangent.stress(), &next.solid.equilibrium().solution.stress);
            near(
                tangent.stress().elements()[0].stress_mandel_pa()[0],
                stress,
                2e-8,
            );
            for value in &tangent.stress().elements()[0].stress_mandel_pa()[1..] {
                near(*value, 0.0, 2e-8);
            }
            let only_material = material.problem(&fixed).assemble(cx).unwrap();
            assert_eq!(
                tangent.elastic_stiffness().to_dense(),
                only_material.stiffness.to_dense()
            );
            assert_eq!(pencil.mass.to_dense(), only_material.mass.to_dense());
            near(pencil.total_mass_kg, 1.0, 1e-14);
            for (row, &dof) in pencil.free_dofs.iter().enumerate() {
                for column in 0..pencil.free_dofs.len() {
                    let expected = if row == column && (dof == 4 || dof == 5) {
                        stress / 6.0
                    } else {
                        0.0
                    };
                    near(
                        tangent.geometric_stiffness().get(row, column),
                        expected,
                        2e-8,
                    );
                    near(
                        pencil.stiffness.get(row, column),
                        only_material.stiffness.get(row, column) + expected,
                        2e-8,
                    );
                }
            }
            let q = [0.02, -0.01, 0.03, 0.01, -0.02];
            assert_eq!(q.len(), pencil.free_dofs.len());
            let mut full = [[0.0; 3]; 4];
            for (&dof, value) in pencil.free_dofs.iter().zip(q) {
                full[dof / 3][dof % 3] = value;
            }
            let gradient =
                core::array::from_fn(|i| core::array::from_fn(|j| full[j + 1][i] - full[0][i]));
            let potential = |s| {
                independent_adiabatic_virtual_energy(
                    s,
                    gradient,
                    temperature,
                    young,
                    alpha,
                    6.0 * cp,
                )
            };
            let h = 0.01;
            let expected = (-potential(-2.0 * h) + 16.0 * potential(-h) - 30.0 * potential(0.0)
                + 16.0 * potential(h)
                - potential(2.0 * h))
                / (12.0 * h * h);
            let mut applied = vec![0.0; q.len()];
            pencil.stiffness.spmv(&q, &mut applied);
            let actual: f64 = q.iter().zip(applied).map(|(a, b)| a * b).sum();
            near(actual, expected, 2e-6);
            let modes = fs_modal::eigh_gen_dense(
                &pencil.stiffness.to_dense(),
                &pencil.mass.to_dense(),
                pencil.free_dofs.len(),
            )
            .unwrap();
            eigenvalues.push(
                modes
                    .into_iter()
                    .map(|mode| mode.lambda)
                    .collect::<Vec<_>>(),
            );
            assert_eq!(material.time_s(), next.solid.time_s());
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(heat.time_s(), 0.0);
            let repeat = next.solid.adiabatic_prestressed_state(cx).unwrap();
            assert_eq!(repeat.tangent().identity(), tangent.identity());
        }
        // Thermal restraint stiffens transverse motion in tension and softens
        // it in compression. Other modes also see the distinct thermal tangent.
        assert!(
            eigenvalues[2]
                .iter()
                .zip(&eigenvalues[1])
                .any(|(hot, base)| *hot < *base - 0.01)
        );
        assert!(
            eigenvalues[0]
                .iter()
                .zip(&eigenvalues[1])
                .any(|(cold, base)| *cold > *base + 0.01)
        );
    });
}
