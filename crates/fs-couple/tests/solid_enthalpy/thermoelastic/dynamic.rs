//! G1/G3/G4: finite source heat coupled to actual inertia, no equilibrium reset.
use super::*;
use fs_solid::linear3::TetMaterialField;
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};

fn dynamic_controls() -> TransientConfig3 {
    TransientConfig3 {
        maximum_steps: 1000,
        maximum_displacement_m: 0.05,
        absolute_force_tolerance_n: 1e-11,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
        ..TransientConfig3::default()
    }
}
fn fixed() -> Vec<usize> {
    (0..12).filter(|d| *d != 3).collect()
}
fn velocity(value: f64) -> [[f64; 3]; 4] {
    [[0.0; 3], [value, 0.0, 0.0], [0.0; 3], [0.0; 3]]
}

// Independent one-DOF consistent-mass equation. The tetrahedron has rho*V=1 kg,
// so its sole free diagonal mass is 0.1 kg; k=V*E with nu=0. The scalar thermal
// root uses explicit Psi and -dPsi/dT from element(), not runtime material calls.
#[allow(clippy::too_many_arguments)]
fn endpoint(
    old: [f64; 3],
    forces: [f64; 2],
    young: f64,
    slope: f64,
    alpha: f64,
    capacity: f64,
    power: f64,
    damping: f64,
    dt: f64,
) -> [f64; 5] {
    let [t0, u0, v0] = old;
    let evaluate = |t: f64| {
        let mass = 0.1;
        let k0 = VOLUME * (young + slope * (t0 - 300.0));
        let k1 = VOLUME * (young + slope * (t - 300.0));
        let c = damping * mass;
        let force =
            0.5 * (forces[0] + forces[1]) + 0.5 * alpha * (k0 * (t0 - 300.0) + k1 * (t - 300.0));
        let u = ((2.0 * mass / dt.powi(2) + c / dt - 0.5 * k0) * u0 + 2.0 * mass / dt * v0 + force)
            / (2.0 * mass / dt.powi(2) + c / dt + 0.5 * k1);
        let v = 2.0 * (u - u0) / dt - v0;
        let before = element(t0, u0, young, slope, alpha);
        let after = element(t, u, young, slope, alpha);
        let work = 0.5 * VOLUME * (before.0 + after.0) * (u - u0);
        let heat = work - (after.1 - before.1);
        let loss = c * ((u - u0) / dt).powi(2) * dt;
        (
            capacity * (t - t0) - power * dt - heat,
            [t, u, v, heat, loss],
        )
    };
    let t = root(280.0, 359.0, |t| evaluate(t).0);
    evaluate(t).1
}

#[test]
fn g1_heated_inertia_matches_independent_coupled_scalar_equations() {
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
        let initial = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
            .unwrap();
        let applied = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 1.2)], cx)
            .unwrap();
        let fixed = fixed();
        let equilibrium = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &initial,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            solve(),
            cx,
        )
        .unwrap();
        let u0 = equilibrium
            .equilibrium()
            .solution
            .displacement
            .displacement_m()[1][0];
        let body = equilibrium
            .begin_transient(&initial, &velocity(0.08), dynamic_controls(), cx)
            .unwrap();
        let expected = endpoint(
            [300.0, u0, 0.08],
            [0.8, 1.2],
            1000.0,
            2.0,
            0.0003,
            0.2,
            2.0,
            0.0,
            0.02,
        );
        let next = body
            .advance(&heat, &(), 0.02, &applied, thermal_controls(), solve(), cx)
            .unwrap();
        let motion = next.solid.nodal_motion(cx).unwrap();
        near(
            next.thermal.temperature_k(IDS[0]).unwrap(),
            expected[0],
            3e-8,
        );
        near(motion[0][1][0], expected[1], 1e-10);
        near(motion[1][1][0], expected[2], 1e-9);
        near(next.report.reversible_heat_j, expected[3], 2e-9);
        near(
            next.report.mechanical_work_j,
            1.0 * (expected[1] - u0),
            1e-10,
        );
        near(
            next.report.kinetic_energy_change_j,
            0.05 * (expected[2].powi(2) - 0.08_f64.powi(2)),
            1e-10,
        );
        near(next.report.external_heat_j, 0.04, 1e-9);
        near(next.report.combined_residual_j, 0.0, 1e-9);
        assert_eq!(next.solid.step(), 1);
        assert_eq!(body.step(), 0);
        assert_eq!(heat.time_s(), 0.0);
    });
}

#[test]
fn g1_adiabatic_oscillation_retains_temperature_feedback_and_native_history() {
    with_cx(|cx| {
        let (card, curve) = fixture(1000.0, 0.0, 0.0008, 0.2, false);
        let mut heat = thermal(&[&curve], &[0.2], &[0.0], 0.0);
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
        let fixed = fixed();
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
            .unwrap();
        let equilibrium = ThermoelasticSolid::try_new(
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
        let u0 = equilibrium
            .equilibrium()
            .solution
            .displacement
            .displacement_m()[1][0];
        let mut body = equilibrium
            .begin_transient(&loads, &velocity(0.05), dynamic_controls(), cx)
            .unwrap();
        let initial_energy = heat.energy_j()
            + body.kinetic_energy_j()
            + body
                .potentials()
                .iter()
                .map(|p| p.internal_energy_j)
                .sum::<f64>();
        let mut expected = [300.0, u0, 0.05];
        let mut largest_temperature_change = 0.0_f64;
        for step in 1..=24 {
            let target = step as f64 * 0.005;
            let dt = target - heat.time_s();
            let reference = endpoint(expected, [0.8, 0.8], 1000.0, 0.0, 0.0008, 0.2, 0.0, 0.0, dt);
            let next = body
                .advance(&heat, &(), target, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let motion = next.solid.nodal_motion(cx).unwrap();
            near(
                next.thermal.temperature_k(IDS[0]).unwrap(),
                reference[0],
                3e-8,
            );
            near(motion[0][1][0], reference[1], 2e-10);
            near(motion[1][1][0], reference[2], 2e-9);
            assert_eq!(next.solid.step(), step);
            near(
                next.thermal.energy_j()
                    + next.solid.kinetic_energy_j()
                    + next
                        .solid
                        .potentials()
                        .iter()
                        .map(|p| p.internal_energy_j)
                        .sum::<f64>(),
                initial_energy + 0.8 * (motion[0][1][0] - u0),
                3e-8,
            );
            largest_temperature_change =
                largest_temperature_change.max((reference[0] - 300.0).abs());
            expected = [reference[0], reference[1], reference[2]];
            let (thermal_next, solid_next, _, _) = next.into_parts();
            heat = thermal_next;
            body = solid_next;
        }
        assert!(largest_temperature_change > 0.1);
        assert_eq!(body.state().samples()[0].reference_mass_kg, 1.0);
    });
}

#[test]
fn g3_uncoupled_limit_keeps_damping_export_distinct_from_cell_heat() {
    with_cx(|cx| {
        let (card, curve) = fixture(1000.0, 0.0, 0.0, 0.2, false);
        let heat = thermal(&[&curve], &[0.2], &[0.0], 0.0);
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
        let fixed = fixed();
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let equilibrium = ThermoelasticSolid::try_new(
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
        let controls = TransientConfig3 {
            mass_damping_per_s: 0.7,
            ..dynamic_controls()
        };
        let body = equilibrium
            .begin_transient(&loads, &velocity(0.1), controls, cx)
            .unwrap();
        let next = body
            .advance(&heat, &(), 0.01, &loads, thermal_controls(), solve(), cx)
            .unwrap();
        let reference = endpoint(
            [300.0, 0.0, 0.1],
            [0.0, 0.0],
            1000.0,
            0.0,
            0.0,
            0.2,
            0.0,
            0.7,
            0.01,
        );
        let motion = next.solid.nodal_motion(cx).unwrap();
        near(motion[0][1][0], reference[1], 1e-11);
        near(motion[1][1][0], reference[2], 1e-10);
        near(next.report.damping_export_j, reference[4], 1e-12);
        assert!(next.report.damping_export_j > 0.0);
        near(next.report.reversible_heat_j, 0.0, 1e-12);
        assert_eq!(next.thermal.cells(), heat.cells());
        near(
            next.solid.kinetic_energy_j()
                + next.solid.potentials()[0].internal_energy_j
                + next.report.damping_export_j,
            body.kinetic_energy_j(),
            1e-10,
        );
    });
}

#[test]
fn g4_failed_feedback_and_late_observer_leave_original_motion_and_heat_paired() {
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
        let fixed = fixed();
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
            .unwrap();
        let equilibrium = ThermoelasticSolid::try_new(
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
        let body = equilibrium
            .begin_transient(&loads, &velocity(0.01), dynamic_controls(), cx)
            .unwrap();
        let before = body.nodal_motion(cx).unwrap();
        let failed = ThermoelasticSolve {
            maximum_iterations: 1,
            ..solve()
        };
        assert!(matches!(
            body.advance(&heat, &(), 0.1, &loads, thermal_controls(), failed, cx),
            Err(ThermoelasticFeedbackError::NotConverged { .. })
        ));
        let observer = Steps {
            count: 0,
            cancel_at: 2,
        };
        assert!(
            body.advance(
                &heat,
                &observer,
                0.1,
                &loads,
                thermal_controls(),
                solve(),
                cx
            )
            .is_err()
        );
        assert_eq!(observer.count, 0);
        assert_eq!(body.nodal_motion(cx).unwrap(), before);
        assert_eq!(body.step(), 0);
        assert_eq!(heat.time_s(), 0.0);
        let a = body
            .advance(&heat, &(), 0.02, &loads, thermal_controls(), solve(), cx)
            .unwrap();
        let b = body
            .advance(&heat, &(), 0.02, &loads, thermal_controls(), solve(), cx)
            .unwrap();
        assert_eq!(
            a.solid.nodal_motion(cx).unwrap(),
            b.solid.nodal_motion(cx).unwrap()
        );
        assert_eq!(a.thermal.cells(), b.thermal.cells());
        assert!(
            a.solid
                .advance(&heat, &(), 0.04, &loads, thermal_controls(), solve(), cx)
                .is_err()
        );
    });
}

#[test]
fn g3_native_checkpoint_keeps_step_budget_and_refuses_material_changes() {
    with_cx(|cx| {
        let (card, curve) = fixture(1000.0, 0.0, 0.0, 0.2, false);
        let heat = thermal(&[&curve], &[0.2], &[0.0], 0.0);
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
        let state = source.resolve(&heat, cx).unwrap();
        let fixed = fixed();
        let problem = state.problem(&fixed);
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let controls = TransientConfig3 {
            maximum_steps: 2,
            ..dynamic_controls()
        };
        let mut direct = loads
            .begin_transient(
                &problem,
                Some(state.thermal_strains()),
                0.0,
                &[[0.0; 3]; 4],
                &velocity(0.1),
                controls,
                cx,
            )
            .unwrap();
        direct.advance(&loads, 0.01, cx).unwrap();
        let checkpoint = direct.clone().into_checkpoint(cx).unwrap();
        let mut rebound = checkpoint.bind(&problem, cx).unwrap();
        let a = direct.advance(&loads, 0.01, cx).unwrap();
        let b = rebound.advance(&loads, 0.01, cx).unwrap();
        assert_eq!(a, b);
        assert_eq!(
            direct.nodal_motion(cx).unwrap(),
            rebound.nodal_motion(cx).unwrap()
        );
        assert_eq!(rebound.state().step, 2);
        assert!(matches!(
            rebound.advance(&loads, 0.01, cx),
            Err(TransientError3::Budget)
        ));
        let mut changed = state.materials().to_vec();
        changed[0].density_kg_m3 *= 2.0;
        let bad = fs_solid::linear3::TetLinearElasticProblem {
            materials: TetMaterialField::PerElement(&changed),
            ..state.problem(&fixed)
        };
        assert!(checkpoint.bind(&bad, cx).is_err());
        let alternate: Vec<_> = (0..12).filter(|d| *d != 4).collect();
        assert!(checkpoint.bind(&state.problem(&alternate), cx).is_err());
    });
}

#[path = "dynamic/cavity.rs"]
mod cavity;
