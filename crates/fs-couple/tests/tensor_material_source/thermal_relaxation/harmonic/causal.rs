//! Retained Maxwell memory heats the original thermal cells without a reset.
use super::*;
use fs_couple::tensor_material::solid::thermal::harmonic::causal_heat_assembly_candidate;
use fs_math::dd::Dd;
use fs_solid::linear3::visco::ViscoHarmonicMaterialHistory;

const PRELOAD: f64 = 0.004;
const PRE_DT: f64 = 0.01;
const TAU: f64 = 0.02;
fn nodal(value: f64) -> [f64; 12] {
    let mut out = [0.0; 12];
    out[11] = value;
    out
}
fn causal_budget() -> ViscoHarmonicMaterialBudget {
    ViscoHarmonicMaterialBudget {
        maximum_points: 1,
        maximum_payload_bytes: 1 << 24,
        maximum_work: 100_000_000,
    }
}
fn accepted<'s>(solid: &'s ViscoElements, cx: &Cx<'_>) -> ViscoHarmonicMaterialHistory<'s> {
    let mut history = solid.state(cx).unwrap();
    history.step_internal_force(&nodal(PRELOAD), PRE_DT, cx).unwrap();
    history.bind_harmonic_clock(Dd::ZERO).unwrap()
}
fn spring0() -> f64 {
    PRELOAD * TAU / PRE_DT * (-(-PRE_DT / TAU).exp_m1())
}
fn full_bank() -> [[C64; 12]; 2] {
    DISPLACEMENTS.map(|[a]| {
        let mut row = [C64::ZERO; 12];
        row[11] = a;
        row
    })
}
fn static_bias() -> [f64; 12] {
    nodal(PRELOAD - DISPLACEMENTS.iter().map(|[z]| z.re).sum::<f64>())
}

// Independent scalar exact initial-value solution, then real-time quadrature.
// This oracle does not use the causal FFT/material/solid implementation.
fn physical_causal(t: f64) -> [f64; 3] {
    let mut strain = static_bias()[11];
    let (mut rate, mut settled, mut start_settled) = (0.0, 0.0, 0.0);
    for (w, [z]) in FREQUENCIES.into_iter().zip(DISPLACEMENTS) {
        let (s, c) = (w * t).sin_cos();
        let x = w * TAU;
        let a = x * x / (1.0 + x * x);
        let b = x / (1.0 + x * x);
        strain += z.re * c + z.im * s;
        rate += w * (z.im * c - z.re * s);
        settled += (a * z.re + b * z.im) * c + (a * z.im - b * z.re) * s;
        start_settled += a * z.re + b * z.im;
    }
    let spring = settled + (spring0() - start_settled) * (-t / TAU).exp();
    let stress = 1000.0 * strain + 280.0 * spring + VISCOSITY * rate;
    [
        VOLUME * (VISCOSITY * rate * rate + 280.0 * spring * spring / TAU),
        VOLUME * stress * rate,
        0.5 * VOLUME * (1000.0 * strain * strain + 280.0 * spring * spring),
    ]
}
fn integrals(start: f64, duration: f64) -> [f64; 2] {
    let n = 4096;
    let mut sums = [0.0; 2];
    for i in 0..=n {
        let weight = if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
        let values = physical_causal(start + duration * i as f64 / n as f64);
        for j in 0..2 { sums[j] += weight * values[j]; }
    }
    sums.map(|v| v * duration / (3.0 * n as f64))
}

#[test]
fn retained_spring_energy_relaxes_into_heat_without_an_ac_drive() {
    with_case(false, false, |cx, _, provider, mut thermal, solid, curve| {
        let observer = solid.prepare_harmonic_material_causal_observer(&[], causal_budget(), cx).unwrap();
        let mut history = accepted(solid, cx);
        let initial = history.history().energy();
        let dt = 0.007;
        let mut temperatures = [300.0; 2];
        let mut deposited = 0.0;
        for i in 1..=3 {
            let target = i as f64 * dt;
            let start = thermal.time_s();
            let h = target - start;
            let expected = 0.5 * VOLUME * 280.0 * spring0().powi(2)
                * (-2.0 * start / TAU).exp() * (-(-2.0 * h / TAU).exp_m1());
            let before = history.history().energy();
            let before_cells = thermal.cells().clone();
            let candidate = causal_heat_assembly_candidate(
                provider, cx, &thermal, target, &observer, &history,
                &[], &nodal(PRELOAD), &receivers(curve), solve_controls(), &[0.0; 2],
            ).unwrap();
            let report = candidate.report();
            near(report.material.material_heat_j, expected, 1e-12);
            near(report.material.kelvin_heat_j, 0.0, 1e-14);
            near(report.material.stress_work_j, 0.0, 1e-14);
            near(report.material.stored_energy_change_j, -expected, 1e-12);
            near(report.deposited_heat_j, expected, 1e-12);
            near(report.constitutive_thermal_energy_residual_j, 0.0, 2e-10);
            near(report.caloric_energy_residual_j, 0.0, 2e-10);
            assert_eq!(candidate.history().time_s_dd(), Dd::from_f64(target));
            assert_eq!(history.history().energy(), before);
            assert_eq!(thermal.cells(), &before_cells);
            let sum = temperatures[0] + temperatures[1] + expected / CAPACITY;
            let difference = (temperatures[0] - temperatures[1] + expected / CAPACITY)
                / (1.0 + 2.0 * 0.01 * h / CAPACITY);
            temperatures = [0.5 * (sum + difference), 0.5 * (sum - difference)];
            for id in [0_u64, 1] {
                near(candidate.thermal().temperature_k(CellId(id)).unwrap(),
                    temperatures[id as usize], 3e-7);
            }
            deposited += expected;
            let (next_history, next_thermal, _) = candidate.into_parts();
            history = next_history;
            thermal = next_thermal;
        }
        near(history.history().energy().stored_j, initial.stored_j - deposited, 1e-12);
        near(history.history().energy().work_j, initial.work_j, 1e-12);
        near(history.history().energy().dissipated_j, initial.dissipated_j + deposited, 1e-12);
    });
}

#[test]
fn two_coherent_thermal_intervals_preserve_original_memory_and_independent_work() {
    with_case(false, false, |cx, _, provider, mut thermal, solid, curve| {
        let observer = solid.prepare_harmonic_material_causal_observer(
            &FREQUENCIES, causal_budget(), cx,
        ).unwrap();
        let fields = full_bank();
        let refs = fields.each_ref().map(|q| q.as_slice());
        let mut history = accepted(solid, cx);
        let initial = history.history().energy();
        let mut total_heat = 0.0;
        let mut total_work = 0.0;
        for i in 1..=2 {
            let target = i as f64 * DT;
            let start = thermal.time_s();
            let [heat, work] = integrals(start, target - start);
            let candidate = causal_heat_assembly_candidate(
                provider, cx, &thermal, target, &observer, &history,
                &refs, &static_bias(), &receivers(curve), solve_controls(), &[0.0; 2],
            ).unwrap();
            let report = candidate.report();
            near(report.material.material_heat_j, heat, 1e-11);
            near(report.material.stress_work_j, work, 1e-11);
            near(report.material.start_stored_energy_j, physical_causal(start)[2], 1e-11);
            near(report.material.end_stored_energy_j, physical_causal(target)[2], 1e-11);
            near(report.material.stored_energy_change_j,
                physical_causal(target)[2] - physical_causal(start)[2], 1e-11);
            near(report.deposited_heat_j, heat, 1e-11);
            near(report.heat_handoff_residual_j, 0.0, 1e-12);
            near(report.constitutive_thermal_energy_residual_j, 0.0, 2e-10);
            total_heat += heat;
            total_work += work;
            let (next_history, next_thermal, _) = candidate.into_parts();
            history = next_history;
            thermal = next_thermal;
        }
        near(total_heat, integrals(0.0, 2.0 * DT)[0], 1e-11);
        near(history.history().energy().work_j, initial.work_j + total_work, 1e-11);
        near(history.history().energy().dissipated_j, initial.dissipated_j + total_heat, 1e-11);
        assert_eq!(history.time_s_dd(), Dd::from_f64(thermal.time_s()));
    });
}

#[test]
fn clock_mismatch_late_refusal_and_cancellation_leave_both_accepted_owners_intact() {
    with_case(false, false, |cx, gate, provider, thermal, solid, curve| {
        let observer = solid.prepare_harmonic_material_causal_observer(&[], causal_budget(), cx).unwrap();
        let history = accepted(solid, cx);
        let before = history.history().energy();
        let before_cells = thermal.cells().clone();
        let shifted = history.clone().into_parts().0.bind_harmonic_clock(Dd::from_f64(1e-20)).unwrap();
        assert!(causal_heat_assembly_candidate(
            provider, cx, &thermal, DT, &observer, &shifted,
            &[], &nodal(PRELOAD), &receivers(curve), solve_controls(), &[0.0; 2],
        ).is_err());
        let late = Late {
            inner: provider, gate, calls: Cell::new(0),
            fail_at: Cell::new(Some(3)), cancel_at: Cell::new(None),
        };
        assert!(causal_heat_assembly_candidate(
            &late, cx, &thermal, DT, &observer, &history,
            &[], &nodal(PRELOAD), &receivers(curve), solve_controls(), &[0.0; 2],
        ).is_err());
        assert_eq!(history.history().energy(), before);
        assert_eq!(thermal.cells(), &before_cells);
        late.calls.set(0);
        late.fail_at.set(None);
        let retry = causal_heat_assembly_candidate(
            &late, cx, &thermal, DT, &observer, &history,
            &[], &nodal(PRELOAD), &receivers(curve), solve_controls(), &[0.0; 2],
        ).unwrap();
        assert!(retry.report().deposited_heat_j > 0.0);
        assert_eq!(history.time_s_dd(), Dd::ZERO);
        late.calls.set(0);
        late.cancel_at.set(Some(3));
        assert!(matches!(
            causal_heat_assembly_candidate(
                &late, cx, &thermal, DT, &observer, &history,
                &[], &nodal(PRELOAD), &receivers(curve), solve_controls(), &[0.0; 2],
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(history.history().energy(), before);
        assert_eq!(thermal.cells(), &before_cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
