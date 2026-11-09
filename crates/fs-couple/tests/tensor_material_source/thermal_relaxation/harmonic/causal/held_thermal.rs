//! Held source eigenstrain survives causal memory relaxation and heat handoff.

use super::*;
use fs_blake3::ContentHash;
use fs_solid::linear3::{TetThermalStrainField, TetThermalStrainState};

const FREE_ZZ: f64 = 0.001;

fn preloaded<'s>(
    solid: &'s ViscoElements,
    cx: &Cx<'_>,
) -> (ViscoHarmonicMaterialHistory<'s>, TetThermalStrainState) {
    let snapshot = TetThermalStrainState::try_new(
        310.0,
        300.0,
        [0.0, 0.0, FREE_ZZ, 0.0, 0.0, 0.0],
        solid.reference_material(0).unwrap().material_state_identity,
        ContentHash([77; 32]),
        0.02,
    )
    .unwrap();
    let mut state = solid.state(cx).unwrap();
    // Original affine total/free-strain ramps leave exactly the same mechanical
    // strain and raw Maxwell spring memory as the existing independent oracle.
    state
        .step_internal_force_with_thermal_strain(
            &nodal(PRELOAD + FREE_ZZ),
            PRE_DT,
            TetThermalStrainField::Uniform(&snapshot),
            cx,
        )
        .unwrap();
    assert!(state.eigenstrain_work_j() < 0.0);
    (state.bind_harmonic_clock(Dd::ZERO).unwrap(), snapshot)
}

#[test]
fn held_eigenstrain_memory_relaxes_into_heat_without_rebinding_its_source() {
    with_case(
        false,
        false,
        |cx, _, provider, mut thermal, solid, curve| {
            let observer = solid
                .prepare_harmonic_material_causal_observer(&[], causal_budget(), cx)
                .unwrap();
            let (mut history, snapshot) = preloaded(solid, cx);
            let initial = history.history().energy();
            let eigenstrain_work = history.history().eigenstrain_work_j();
            let mut total_heat = 0.0;
            for step in 1..=2 {
                let target = step as f64 * 0.007;
                let start = thermal.time_s();
                let h = target - start;
                let heat = 0.5
                    * VOLUME
                    * 280.0
                    * spring0().powi(2)
                    * (-2.0 * start / TAU).exp()
                    * (-(-2.0 * h / TAU).exp_m1());
                let before = history.history().energy();
                let before_cells = thermal.cells().clone();
                let candidate = causal_heat_assembly_candidate(
                    provider,
                    cx,
                    &thermal,
                    target,
                    &observer,
                    &history,
                    &[],
                    &nodal(PRELOAD + FREE_ZZ),
                    &receivers(curve),
                    solve_controls(),
                    &[0.0; 2],
                )
                .unwrap();
                let report = candidate.report();
                near(report.material.material_heat_j, heat, 1e-12);
                near(report.material.kelvin_heat_j, 0.0, 1e-14);
                near(report.material.stress_work_j, 0.0, 1e-14);
                near(report.material.stored_energy_change_j, -heat, 1e-12);
                near(report.deposited_heat_j, heat, 1e-12);
                near(report.constitutive_thermal_energy_residual_j, 0.0, 2e-10);
                near(report.caloric_energy_residual_j, 0.0, 2e-10);
                assert_eq!(
                    candidate.history().history().thermal_strains(),
                    Some(core::slice::from_ref(&snapshot))
                );
                assert_eq!(
                    candidate.history().history().eigenstrain_work_j(),
                    eigenstrain_work
                );
                assert_eq!(
                    candidate.history().history().displacement_m(),
                    &nodal(PRELOAD + FREE_ZZ)
                );
                assert!(candidate.thermal().temperature_k(CellId(0)).unwrap() > 300.0);
                assert_eq!(history.history().energy(), before);
                assert_eq!(thermal.cells(), &before_cells);
                total_heat += heat;
                let (next_history, next_thermal, _) = candidate.into_parts();
                history = next_history;
                thermal = next_thermal;
            }
            near(
                history.history().energy().stored_j,
                initial.stored_j - total_heat,
                1e-12,
            );
            near(history.history().energy().work_j, initial.work_j, 1e-12);
            near(
                history.history().energy().dissipated_j,
                initial.dissipated_j + total_heat,
                1e-12,
            );
            assert_eq!(
                history.history().thermal_strains(),
                Some(core::slice::from_ref(&snapshot))
            );
            assert_eq!(history.time_s_dd(), Dd::from_f64(thermal.time_s()));
        },
    );
}

#[test]
fn held_eigenstrain_coherent_heat_is_atomic_across_late_refusal_and_retry() {
    with_case(false, false, |cx, gate, provider, thermal, solid, curve| {
        let observer = solid
            .prepare_harmonic_material_causal_observer(&FREQUENCIES, causal_budget(), cx)
            .unwrap();
        let (history, snapshot) = preloaded(solid, cx);
        let initial = history.history().energy();
        let eigenstrain_work = history.history().eigenstrain_work_j();
        let before_cells = thermal.cells().clone();
        let fields = full_bank();
        let refs = fields.each_ref().map(|q| q.as_slice());
        let mut bias = static_bias();
        bias[11] += FREE_ZZ;
        let late = Late {
            inner: provider,
            gate,
            calls: Cell::new(0),
            fail_at: Cell::new(Some(3)),
            cancel_at: Cell::new(None),
        };
        assert!(
            causal_heat_assembly_candidate(
                &late,
                cx,
                &thermal,
                DT,
                &observer,
                &history,
                &refs,
                &bias,
                &receivers(curve),
                solve_controls(),
                &[0.0; 2],
            )
            .is_err()
        );
        assert_eq!(late.calls.get(), 3);
        assert_eq!(history.history().energy(), initial);
        assert_eq!(
            history.history().thermal_strains(),
            Some(core::slice::from_ref(&snapshot))
        );
        assert_eq!(history.history().eigenstrain_work_j(), eigenstrain_work);
        assert_eq!(
            history.history().displacement_m(),
            &nodal(PRELOAD + FREE_ZZ)
        );
        assert_eq!(thermal.cells(), &before_cells);
        late.calls.set(0);
        late.fail_at.set(None);
        let candidate = causal_heat_assembly_candidate(
            &late,
            cx,
            &thermal,
            DT,
            &observer,
            &history,
            &refs,
            &bias,
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        let [heat, work] = integrals(0.0, DT);
        let report = candidate.report();
        near(report.material.material_heat_j, heat, 1e-11);
        near(report.material.stress_work_j, work, 1e-11);
        near(
            report.material.start_stored_energy_j,
            physical_causal(0.0)[2],
            1e-11,
        );
        near(
            report.material.end_stored_energy_j,
            physical_causal(DT)[2],
            1e-11,
        );
        near(report.deposited_heat_j, heat, 1e-11);
        near(report.constitutive_thermal_energy_residual_j, 0.0, 2e-10);
        near(
            candidate.history().history().energy().work_j,
            initial.work_j + work,
            1e-11,
        );
        near(
            candidate.history().history().energy().dissipated_j,
            initial.dissipated_j + heat,
            1e-11,
        );
        assert_eq!(
            candidate.history().history().thermal_strains(),
            Some(core::slice::from_ref(&snapshot))
        );
        assert_eq!(
            candidate.history().history().eigenstrain_work_j(),
            eigenstrain_work
        );
        late.calls.set(0);
        late.cancel_at.set(Some(3));
        assert!(matches!(
            causal_heat_assembly_candidate(
                &late,
                cx,
                &thermal,
                DT,
                &observer,
                &history,
                &refs,
                &bias,
                &receivers(curve),
                solve_controls(),
                &[0.0; 2],
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(history.history().energy(), initial);
        assert_eq!(
            history.history().thermal_strains(),
            Some(core::slice::from_ref(&snapshot))
        );
        assert_eq!(history.history().eigenstrain_work_j(), eigenstrain_work);
        assert_eq!(thermal.cells(), &before_cells);
        assert_eq!(history.time_s_dd(), Dd::ZERO);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
