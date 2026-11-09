//! Tiny physical increments under substantial prestress must not subtract
//! rounded absolute elastic energies. These tests use the existing public
//! static and inertial owners and independent constant-E, nu=0 invariants.
use super::*;
use fs_solid::loads3::transient::TransientConfig3;

const YOUNG: f64 = 1_073_741_824.0;
const ALPHA: f64 = 1.0 / 32768.0;
const CP: f64 = 1024.0;
const PRESTRAIN: f64 = 1.0 / 256.0;

fn precision_solve() -> ThermoelasticSolve {
    ThermoelasticSolve {
        maximum_iterations: 64,
        heat_relaxation: 1.0,
        heat_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-8,
        stress_tolerance_pa: 1e-4,
    }
}

fn precision_thermal() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-8,
    }
}

fn component_close(actual: f64, expected: f64) {
    let tolerance = 2e-10 * expected.abs() + 1e-25;
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

// C=E I in Mandel coordinates for nu=0. For constant C and alpha,
// delta S=V E alpha*(trace(delta eps)-3 alpha delta T), and Q=-Tmean delta S.
// This oracle uses input increments, never an absolute stored-energy difference
// or the material exchange implementation.
fn increment(
    before: [f64; 6],
    after: [f64; 6],
    temperature: [f64; 2],
) -> [f64; 3] {
    let mean_t = 0.5 * (temperature[0] + temperature[1]);
    let delta_t = temperature[1] - temperature[0];
    let delta: [f64; 6] = std::array::from_fn(|i| after[i] - before[i]);
    let mean_elastic: [f64; 6] = std::array::from_fn(|i| {
        0.5 * (before[i] + after[i])
            - if i < 3 { ALPHA * (mean_t - 300.0) } else { 0.0 }
    });
    let work = VOLUME * YOUNG
        * (0..6).map(|i| mean_elastic[i] * delta[i]).sum::<f64>();
    let entropy_change = VOLUME * YOUNG * ALPHA
        * (delta[0] + delta[1] + delta[2] - 3.0 * ALPHA * delta_t);
    let heat = -mean_t * entropy_change;
    [work, work - heat, heat]
}

#[test]
fn loaded_static_tiny_increment_retains_reversible_heat_and_internal_storage() {
    with_cx(|cx| {
        let (card, curve) = fixture(YOUNG, 0.0, ALPHA, CP, false);
        let heat = thermal(&[&curve], &[CP], &[0.0], 0.0);
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
        let force = VOLUME * YOUNG * PRESTRAIN;
        let initial = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, force)], cx)
            .unwrap();
        let applied = source
            .geometry()
            .assemble(
                &[],
                VolumeLoad3::None,
                &[(3, f64::from_bits(force.to_bits() + 32))],
                cx,
            )
            .unwrap();
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
        let body = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &initial,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            precision_solve(),
            cx,
        )
        .unwrap();
        let before = *body.equilibrium().solution.stress.elements()[0].total_strain_mandel();
        let old_u = body.equilibrium().solution.displacement.clone();
        let next = body
            .advance(
                &heat,
                &(),
                1.0,
                &applied,
                precision_thermal(),
                precision_solve(),
                cx,
            )
            .unwrap();
        let after = *next.solid.equilibrium().solution.stress.elements()[0].total_strain_mandel();
        let expected = increment(
            before,
            after,
            [
                heat.temperature_k(IDS[0]).unwrap(),
                next.thermal.temperature_k(IDS[0]).unwrap(),
            ],
        );
        assert!(after[0] > before[0]);
        assert!(expected[2].abs() > 1e-15 && expected[2].abs() < 1e-8);
        // The absolute interval tolerances allow this tiny heat to remain
        // unresolved in caloric storage. Constitutive component accuracy is
        // still required; these assertions do not claim a resolved tiny delta T.
        component_close(next.report.cells[0].mechanical_work_j, expected[0]);
        component_close(next.report.internal_energy_change_j, expected[1]);
        component_close(next.report.cells[0].reversible_heat_j, expected[2]);
        component_close(next.report.reversible_heat_j, expected[2]);
        assert!(next.report.cells[0].free_energy_residual_j.is_finite());
        assert_eq!(body.equilibrium().solution.displacement, old_u);
        assert_eq!(body.time_s(), 0.0);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(heat.temperature_k(IDS[0]).unwrap(), 300.0);
    });
}

#[test]
fn tiny_native_inertial_motion_retains_the_same_source_exchange_under_prestress() {
    with_cx(|cx| {
        let (card, curve) = fixture(YOUNG, 0.0, ALPHA, CP, false);
        let heat = thermal(&[&curve], &[CP], &[0.0], 0.0);
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
        let force = VOLUME * YOUNG * PRESTRAIN;
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, force)], cx)
            .unwrap();
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
        let equilibrium = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &loads,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            precision_solve(),
            cx,
        )
        .unwrap();
        let speed = 1.0 / 17_179_869_184.0;
        let velocity = [[0.0; 3], [speed, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let body = equilibrium
            .begin_transient(
                &loads,
                &velocity,
                TransientConfig3 {
                    maximum_steps: 8,
                    maximum_displacement_m: 0.05,
                    absolute_force_tolerance_n: 1e-8,
                    absolute_energy_tolerance_j: 1e-8,
                    relative_energy_tolerance: 1e-8,
                    ..TransientConfig3::default()
                },
                cx,
            )
            .unwrap();
        let before = *body.response().stress.elements()[0].total_strain_mandel();
        let old_motion = body.nodal_motion(cx).unwrap();
        let next = body
            .advance(
                &heat,
                &(),
                1.0 / 65536.0,
                &loads,
                precision_thermal(),
                precision_solve(),
                cx,
            )
            .unwrap();
        let after = *next.solid.response().stress.elements()[0].total_strain_mandel();
        let expected = increment(
            before,
            after,
            [
                heat.temperature_k(IDS[0]).unwrap(),
                next.thermal.temperature_k(IDS[0]).unwrap(),
            ],
        );
        assert!(after[0] > before[0]);
        assert!(expected[2].abs() > 1e-15 && expected[2].abs() < 1e-8);
        // Use actual accepted native strain, independently of the trajectory's
        // force tolerance. This checks finite material accounting, not a claim
        // of relative accuracy below the native physical solve tolerances.
        component_close(next.report.cells[0].mechanical_work_j, expected[0]);
        component_close(next.report.internal_energy_change_j, expected[1]);
        component_close(next.report.cells[0].reversible_heat_j, expected[2]);
        component_close(next.report.reversible_heat_j, expected[2]);
        assert_eq!(next.report.damping_export_j, 0.0);
        assert_eq!(next.solid.step(), 1);
        assert_eq!(body.step(), 0);
        assert_eq!(body.nodal_motion(cx).unwrap(), old_motion);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(heat.temperature_k(IDS[0]).unwrap(), 300.0);
    });
}
