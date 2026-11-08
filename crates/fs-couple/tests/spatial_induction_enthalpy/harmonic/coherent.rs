//! G1/G3/G4: coherent beat heat reaches the original nonlinear caloric receiver.
use super::*;
use fs_couple::circuit_enthalpy::spatial_induction::harmonic::coherent::{
    CoherentInductionHeatError, heat_assembly_candidate as coherent_candidate,
};
use fs_em::spatial::vector_potential::induction::harmonic::coherent_exposure::SpatialHarmonicInductionExposureLimits;

fn limits() -> SpatialHarmonicInductionExposureLimits {
    SpatialHarmonicInductionExposureLimits {
        maximum_lines: 4,
        maximum_pair_products: 100000,
        maximum_workspace_bytes: 1 << 20,
    }
}
const TONES: [(f64, C64); 2] = [
    (OMEGA, C64::new(0.8 * PEAK, -0.6 * PEAK)),
    (OMEGA + 1.3, C64::new(0.35 * PEAK, 0.5 * PEAK)),
];

// Independent original one-edge A-phi solution: K=4, M=3*sigma/40,
// A=J/(2*(K-i*w*M)). Integrate M*(-dA/dt)^2 in physical time;
// no harmonic-window or field-recovery calculation enters this oracle.
fn physical_power(t: f64) -> f64 {
    let mass = 3.0 / (40.0 * RHO);
    let electric = TONES
        .into_iter()
        .map(|(w, j)| {
            let denominator = 16.0 + (w * mass).powi(2);
            let a = j * C64::new(2.0 / denominator, 0.5 * w * mass / denominator);
            let e = C64::new(-w * a.im, w * a.re);
            e.re * (w * t).cos() + e.im * (w * t).sin()
        })
        .sum::<f64>();
    mass * electric * electric
}
fn integrated_heat(start: f64, dt: f64) -> f64 {
    let n = 32768;
    let h = dt / f64::from(n);
    let mut sum = physical_power(start) + physical_power(start + dt);
    for i in 1..n {
        sum += f64::from(if i % 2 == 0 { 2 } else { 4 }) * physical_power(start + f64::from(i) * h);
    }
    h * sum / 3.0
}

#[test]
fn coherent_beat_heat_drives_sourced_contacts_and_radiation_across_intervals() {
    with_case_with_alpha(false, 0.5, 0.0, |cx, _, heat, mut thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let field = initial(&problem, BASIS, 0.7, 1.0);
        let original_a = field.potential_wb().to_vec();
        let responses = TONES
            .into_iter()
            .map(|(omega, j)| {
                field
                    .solve_harmonic(
                        &[[j; 3]; 6],
                        &[[0.0; 3]; 8],
                        omega,
                        ac_controls(),
                        &mut || false,
                    )
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let bank = field
            .prepare_harmonic_exposure(&responses.iter().collect::<Vec<_>>(), limits(), &mut || {
                false
            })
            .unwrap();
        let cycle_mean: f64 = responses
            .iter()
            .map(|r| r.balance().mean_ohmic_power_w)
            .sum();
        let mut old_temperature = [300.0, 340.0];
        let mut actual_heat = Vec::new();
        for interval in 1..=2 {
            let start = thermal.time_s();
            let expected_heat = integrated_heat(start, DT);
            let expected_temperature = reference(old_temperature, 0.5, Some(expected_heat / DT));
            let old_cells = thermal.cells().clone();
            let candidate = coherent_candidate(
                &prepared,
                cx,
                &thermal,
                f64::from(interval) * DT,
                &bank,
                &regions(),
                1e-10,
                thermal_controls(),
                &[0.0; 2],
            )
            .unwrap();
            let report = candidate.report();
            near(report.deposited_heat_j, expected_heat, 3e-7);
            near(report.exposure.ohmic_heat_j, expected_heat, 3e-7);
            near(report.caloric_energy_residual_j, 0.0, 2e-8);
            near(report.other_external_heat_j, 0.0, 2e-8);
            near(
                candidate.thermal().energy_j() - thermal.energy_j(),
                expected_heat,
                3e-7,
            );
            assert!((expected_heat - DT * cycle_mean).abs() > 0.1);
            for i in 0..2 {
                near(
                    candidate.thermal().temperature_k([BODY, SINK][i]).unwrap(),
                    expected_temperature[i],
                    3e-7,
                );
            }
            assert_eq!(report.exposure.window.start_s(), start);
            assert_eq!(report.exposure.cells().len(), 6);
            assert_eq!(report.source.network.thermal_contacts.len(), 1);
            assert!(report.thermal_trials > 1);
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(field.potential_wb(), original_a);
            assert_eq!(field.time_s(), 0.0);
            actual_heat.push(report.deposited_heat_j);
            thermal = candidate.into_parts().0;
            old_temperature = expected_temperature;
        }
        assert!((actual_heat[0] - actual_heat[1]).abs() > 0.1);
        near(
            actual_heat.into_iter().sum(),
            integrated_heat(0.0, 2.0 * DT),
            6e-7,
        );
    });
}

#[test]
fn frozen_sources_refuse_temperature_dependence_even_without_heating() {
    with_case(false, 0.5, |cx, _, heat, thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let field = initial(&problem, BASIS, 0.0, 1.0);
        let response = field
            .solve_harmonic(
                &[[C64::new(0.0, 0.0); 3]; 6],
                &[[0.0; 3]; 8],
                OMEGA,
                ac_controls(),
                &mut || false,
            )
            .unwrap();
        let bank = field
            .prepare_harmonic_exposure(&[&response], limits(), &mut || false)
            .unwrap();
        let old_cells = thermal.cells().clone();
        assert!(matches!(
            coherent_candidate(
                &prepared,
                cx,
                &thermal,
                DT,
                &bank,
                &regions(),
                1e-10,
                thermal_controls(),
                &[0.0; 2],
            ),
            Err(CircuitEnthalpyError::Source(
                CoherentInductionHeatError::Material(MagneticThermalError::Invalid(
                    "coherent electrical sources must be independent of evolving temperature"
                ))
            ))
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(field.time_s(), 0.0);
    });
}

#[test]
fn late_heat_refusal_and_cancellation_leave_the_bank_and_caloric_state_reusable() {
    with_case_with_alpha(false, 0.5, 0.0, |cx, gate, heat, thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let field = initial(&problem, BASIS, 0.0, 1.0);
        let response = field
            .solve_harmonic(&drive(), &[[0.0; 3]; 8], OMEGA, ac_controls(), &mut || {
                false
            })
            .unwrap();
        let bank = field
            .prepare_harmonic_exposure(&[&response], limits(), &mut || false)
            .unwrap();
        let tagged = Tagged {
            inner: &prepared,
            gate,
            calls: Cell::new(0),
            fail_at: Cell::new(Some(3)),
            cancel_at: Cell::new(None),
        };
        let old_cells = thermal.cells().clone();
        assert!(matches!(
            coherent_candidate(
                &tagged,
                cx,
                &thermal,
                DT,
                &bank,
                &regions(),
                1e-10,
                thermal_controls(),
                &[0.0; 2],
            ),
            Err(CircuitEnthalpyError::Source(
                CoherentInductionHeatError::Heat(TaggedError::Refused)
            ))
        ));
        tagged.calls.set(0);
        tagged.fail_at.set(None);
        let first = coherent_candidate(
            &tagged,
            cx,
            &thermal,
            DT,
            &bank,
            &regions(),
            1e-10,
            thermal_controls(),
            &[0.0; 2],
        )
        .unwrap();
        let repeated = coherent_candidate(
            &prepared,
            cx,
            &thermal,
            DT,
            &bank,
            &regions(),
            1e-10,
            thermal_controls(),
            &[0.0; 2],
        )
        .unwrap();
        assert_eq!(first.thermal().cells(), repeated.thermal().cells());
        assert_eq!(
            first.report().deposited_heat_j.to_bits(),
            repeated.report().deposited_heat_j.to_bits()
        );
        tagged.calls.set(0);
        tagged.cancel_at.set(Some(4));
        assert!(matches!(
            coherent_candidate(
                &tagged,
                cx,
                &thermal,
                DT,
                &bank,
                &regions(),
                1e-10,
                thermal_controls(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(field.time_s(), 0.0);
    });
}
