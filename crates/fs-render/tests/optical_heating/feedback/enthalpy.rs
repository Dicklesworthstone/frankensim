//! Direct incident absorption follows the original latent enthalpy inventory.

use super::*;
use fs_render::volumes::heating::feedback::enthalpy::{
    EnthalpyAbsorptionTable, EnthalpyOpticalCell,
};
use fs_render::volumes::heating::feedback::nonlinear::{
    TemperatureTwoPhotonTable, ThermalTwoPhotonCell, TwoPhotonFeedback,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::{DepositionObserver, ThermalStorageBasis};

/// Manufactured sharp latent interval with an explicitly authored energy datum.
struct Plateau {
    datum: f64,
    identity: u8,
    internal: bool,
}
impl ThermalLaw for Plateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        (self.identity != 0).then_some([self.identity; 32])
    }
    fn storage_basis(&self) -> ThermalStorageBasis {
        if self.internal {
            ThermalStorageBasis::InternalEnergy
        } else {
            ThermalStorageBasis::Enthalpy
        }
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let x = h - self.datum;
        if !x.is_finite() || !(-100.0..=200.0).contains(&x) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if x < 0.0 {
            500.0 + x
        } else if x <= 100.0 {
            500.0
        } else {
            500.0 + x - 100.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static PLATEAU: Plateau = Plateau {
    datum: 0.0,
    identity: 113,
    internal: false,
};
static SHIFTED: Plateau = Plateau {
    datum: -1000.0,
    identity: 114,
    internal: false,
};
static INTERNAL: Plateau = Plateau {
    datum: 0.0,
    identity: 113,
    internal: true,
};
static UNIDENTIFIED: Plateau = Plateau {
    datum: 0.0,
    identity: 0,
    internal: false,
};

fn specimen(
    law: &'static Plateau,
    states: &[(f64, f64)],
    events: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let cells = states
        .iter()
        .enumerate()
        .map(|(index, &(latent_h, source_w))| {
            let mut c = cell((index as u64 + 1) * 11, 1.0, 1.0, source_w);
            c.specific_h_j_kg = law.datum + latent_h;
            c
        })
        .collect();
    thermal(cells, law, events)
}

fn htable(law: &Plateau, span: (f64, f64), coefficients: (f64, f64)) -> EnthalpyAbsorptionTable {
    EnthalpyAbsorptionTable::try_new(
        ContentHash([115; 32]),
        ContentHash([law.identity; 32]),
        (500.0, 500.0),
        &[
            (
                law.datum + span.0,
                MassAbsorptionSpectrum::constant(coefficients.0).unwrap(),
            ),
            (
                law.datum + span.1,
                MassAbsorptionSpectrum::constant(coefficients.1).unwrap(),
            ),
        ],
        2,
    )
    .unwrap()
}

fn hbinding(index: usize, table: usize) -> EnthalpyOpticalCell {
    EnthalpyOpticalCell {
        grid_cell: index,
        thermal_cell: CellId((index as u64 + 1) * 11),
        material_id: MaterialId(7),
        absorption_table: table,
    }
}

fn input(
    cx: &Cx<'_>,
    original: &DepositionSimulation<'_>,
    power: f64,
    tables: &[EnthalpyAbsorptionTable],
    bindings: &[EnthalpyOpticalCell],
    work: OpticalFeedbackBudget,
) -> Result<TemperatureOpticalHeating<1>, OpticalFeedbackError> {
    let density = vec![1.0; bindings.len()];
    let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
    TemperatureOpticalHeating::try_new_enthalpy(
        cx,
        &medium,
        beam(power, density.len()),
        tables,
        bindings,
        original,
        work,
        1.0e-11,
        || Ok(()),
    )
}

// rho=1, length=1, mass=1: h'=P*(1-exp[-a-b*h]).
// Transforming y=exp(b*h) gives a linear constant-coefficient ODE.
fn exact_h(initial: f64, power: f64, time: f64) -> f64 {
    initial + (-(-(0.2 + 0.01 * initial)).exp_m1() * (0.01 * power * time).exp_m1()).ln_1p() / 0.01
}

#[test]
fn g1_beam_absorption_evolves_on_the_native_latent_plateau_and_refines() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, &[(10.0, 0.0)], vec![]);
        let optics = input(
            cx,
            &original,
            20.0,
            &[htable(&PLATEAU, (0.0, 100.0), (0.2, 1.2))],
            &[hbinding(0, 0)],
            budget(),
        )
        .unwrap();
        let exact = exact_h(10.0, 20.0, 4.0);
        let mut errors = Vec::new();
        for dt in [0.25, 0.125, 0.0625] {
            let (next, accepted, audit) = original
                .advance_observed(&optics, 4.0, config(dt), || Ok(()))
                .unwrap();
            let h = next.cells()[&CellId(11)].specific_h_j_kg;
            errors.push(exact - h);
            assert!(h > 10.0 && h < exact && h < 100.0);
            assert_eq!(next.temperature_k(CellId(11)).unwrap(), 500.0);
            assert!(accepted.power().absorbed_w[0] > optics.power().absorbed_w[0]);
            assert_ne!(
                accepted.power().medium_identity,
                optics.power().medium_identity
            );
            close(h - 10.0, accepted.energy().thermal_input_j);
            close(audit.external_heat_j, accepted.energy().absorbed_j[0]);
            close(
                accepted.energy().absorbed_j[0] + accepted.energy().transmitted_j[0],
                80.0,
            );
        }
        for pair in errors.windows(2) {
            assert!((1.9..2.1).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
        let (next, accepted, report) = optics
            .advance_refined(
                &original,
                4.0,
                FeedbackRefinement {
                    enthalpy_absolute_tolerance_j_kg: 0.05,
                    temperature_absolute_tolerance_k: 1.0e-9,
                    optical_energy_absolute_tolerance_j: 0.05,
                    ..refinement()
                },
                || Ok(()),
            )
            .unwrap();
        assert!((next.cells()[&CellId(11)].specific_h_j_kg - exact).abs() < 0.055);
        assert_eq!(next.temperature_k(CellId(11)).unwrap(), 500.0);
        assert!((0.9..1.1).contains(&report.observed_order.unwrap()));
        assert!(report.normalized_difference <= 1.0);
        close(
            accepted.energy().thermal_input_j,
            report.thermal.external_heat_j,
        );
        assert_eq!(optics.energy().thermal_input_j, 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 10.0);
    });
}

#[test]
fn g3_shifted_caloric_datum_and_constant_law_preserve_physical_heat() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, &[(10.0, 0.0)], vec![]);
        let shifted = specimen(&SHIFTED, &[(10.0, 0.0)], vec![]);
        let a = input(
            cx,
            &original,
            20.0,
            &[htable(&PLATEAU, (0.0, 100.0), (0.2, 1.2))],
            &[hbinding(0, 0)],
            budget(),
        )
        .unwrap();
        let b = input(
            cx,
            &shifted,
            20.0,
            &[htable(&SHIFTED, (0.0, 100.0), (0.2, 1.2))],
            &[hbinding(0, 0)],
            budget(),
        )
        .unwrap();
        assert_ne!(a.power().medium_identity, b.power().medium_identity);
        let (a_end, a_heat, _) = original
            .advance_observed(&a, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (b_end, b_heat, _) = shifted
            .advance_observed(&b, 2.0, config(0.125), || Ok(()))
            .unwrap();
        close(
            a_end.cells()[&CellId(11)].specific_h_j_kg,
            b_end.cells()[&CellId(11)].specific_h_j_kg + 1000.0,
        );
        close(
            a_heat.energy().thermal_input_j,
            b_heat.energy().thermal_input_j,
        );
        let constant_h = input(
            cx,
            &original,
            20.0,
            &[htable(&PLATEAU, (0.0, 100.0), (0.4, 0.4))],
            &[hbinding(0, 0)],
            budget(),
        )
        .unwrap();
        let constant_t = feedback(
            cx,
            &original,
            20.0,
            &[table(450.0, 550.0, 0.4, 0.4)],
            &[binding(0, 11, 0)],
        );
        let (h_end, h_heat, _) = original
            .advance_observed(&constant_h, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (t_end, t_heat, _) = original
            .advance_observed(&constant_t, 2.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(h_end.cells(), t_end.cells());
        assert_eq!(h_heat.energy(), t_heat.energy());
    });
}

#[test]
fn g3_local_enthalpy_changes_downstream_shielding_without_changing_temperature() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, &[(10.0, 0.0), (10.0, 0.0)], vec![]);
        let optics = input(
            cx,
            &original,
            20.0,
            &[
                htable(&PLATEAU, (0.0, 100.0), (0.2, 1.2)),
                htable(&PLATEAU, (0.0, 100.0), (0.5, 0.5)),
            ],
            &[hbinding(0, 0), hbinding(1, 1)],
            budget(),
        )
        .unwrap();
        let (next, accepted, audit) = original
            .advance_observed(&optics, 4.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(next.temperature_k(CellId(11)).unwrap(), 500.0);
        assert_eq!(next.temperature_k(CellId(22)).unwrap(), 500.0);
        assert!(accepted.sources_w()[&CellId(11)] > optics.sources_w()[&CellId(11)]);
        assert!(accepted.sources_w()[&CellId(22)] < optics.sources_w()[&CellId(22)]);
        let front_heat = next.cells()[&CellId(11)].specific_h_j_kg - 10.0;
        let rear_heat = next.cells()[&CellId(22)].specific_h_j_kg - 10.0;
        close(rear_heat, (80.0 - front_heat) * -(-0.5_f64).exp_m1());
        close(front_heat + rear_heat, accepted.energy().thermal_input_j);
        close(audit.stored_energy_change_j, front_heat + rear_heat);
    });
}

fn nonlinear_input(
    cx: &Cx<'_>,
    original: &DepositionSimulation<'_>,
    work: OpticalFeedbackBudget,
) -> Result<TemperatureOpticalHeating<1>, OpticalFeedbackError> {
    let medium = grid(cx, &[1.0], MassAbsorptionSpectrum::constant(0.0).unwrap());
    let beta = TemperatureTwoPhotonTable::try_new(
        ContentHash([116; 32]),
        &[(450.0, 0.05), (550.0, 0.15)],
        (400.0, 800.0),
        10.0,
        2,
    )
    .unwrap();
    TemperatureOpticalHeating::try_new_enthalpy_two_photon(
        cx,
        &medium,
        beam(6.0, 1),
        &[htable(&PLATEAU, (0.0, 100.0), (0.2, 1.2))],
        &[hbinding(0, 0)],
        TwoPhotonFeedback {
            beam_area_m2: 1.0,
            tables: &[beta],
            cells: &[ThermalTwoPhotonCell {
                grid_cell: 0,
                two_photon_table: 0,
            }],
        },
        original,
        work,
        1.0e-11,
        || Ok(()),
    )
}

// Independent RK4 using the closed spatial solution of dP/dz=-alpha*P-beta*P².
fn ramp_reference(mut h: f64, duration: f64, start: f64, end: f64, external: f64) -> f64 {
    let rhs = |time: f64, enthalpy: f64| {
        let p = start + (end - start) * time / duration;
        let alpha = 0.2 + 0.01 * enthalpy;
        let transmitted = p * (-alpha).exp() / (1.0 + 0.1 * p * -(-alpha).exp_m1() / alpha);
        external + p - transmitted
    };
    let dt = duration / 4096.0;
    for index in 0..4096 {
        let t = f64::from(index) * dt;
        let k1 = rhs(t, h);
        let k2 = rhs(t + dt * 0.5, h + dt * 0.5 * k1);
        let k3 = rhs(t + dt * 0.5, h + dt * 0.5 * k2);
        let k4 = rhs(t + dt, h + dt * k3);
        h += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    h
}

#[test]
fn g1_enthalpy_and_two_photon_pulses_refine_and_restart_on_native_events() {
    with_cx(|cx| {
        let original = specimen(
            &PLATEAU,
            &[(10.0, 0.0)],
            vec![TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(11), boundary(1.0)),
            }],
        );
        let schedule = SpectralPowerSchedule::try_new(
            &[
                SpectralPowerSegment {
                    start_s: 0.0,
                    end_s: 1.0,
                    start_power_w: [2.0],
                    end_power_w: [6.0],
                },
                SpectralPowerSegment {
                    start_s: 1.5,
                    end_s: 2.0,
                    start_power_w: [3.0],
                    end_power_w: [3.0],
                },
            ],
            2,
        )
        .unwrap();
        let optics = nonlinear_input(cx, &original, budget())
            .unwrap()
            .with_power_schedule(&original, schedule, 2, || Ok(()))
            .unwrap();
        let h1 = ramp_reference(10.0, 1.0, 2.0, 6.0, 0.0);
        let exact = ramp_reference(h1 + 0.5, 0.5, 3.0, 3.0, 1.0);
        let mut errors = Vec::new();
        for dt in [0.125, 0.0625, 0.03125] {
            let (next, accepted, audit) = original
                .advance_observed(&optics, 2.0, config(dt), || Ok(()))
                .unwrap();
            let h = next.cells()[&CellId(11)].specific_h_j_kg;
            errors.push(exact - h);
            assert_eq!(next.temperature_k(CellId(11)).unwrap(), 500.0);
            close(accepted.energy().incident_j[0], 5.5);
            close(h - 10.0, accepted.energy().thermal_input_j + 1.0);
            close(audit.external_heat_j, h - 10.0);
            assert_eq!(audit.events, 1);
        }
        for pair in errors.windows(2) {
            assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
        let (refined, _, report) = optics
            .advance_refined(&original, 2.0, refinement(), || Ok(()))
            .unwrap();
        assert!((refined.cells()[&CellId(11)].specific_h_j_kg - exact).abs() < 0.0041);
        assert!(report.normalized_difference <= 1.0);
        let (whole, whole_heat, _) = original
            .advance_observed(&optics, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, mid_heat, _) = original
            .advance_observed(&optics, 1.5, config(0.125), || Ok(()))
            .unwrap();
        let (continued, continued_heat, _) = middle
            .advance_observed(&mid_heat, 2.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(whole.cells(), continued.cells());
        assert_eq!(whole_heat.energy(), continued_heat.energy());
    });
}

#[test]
fn g4_original_caloric_binding_dark_endpoints_work_and_cancellation_are_binding() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, &[(1.0, 10.0)], vec![]);
        let tables = [htable(&PLATEAU, (0.0, 5.0), (0.0, 0.0))];
        let bindings = [hbinding(0, 0)];
        let dark = input(cx, &original, 0.0, &tables, &bindings, budget()).unwrap();
        assert!(matches!(
            original.advance_observed(&dark, 0.5, config(0.5), || Ok(())),
            Err(OpticalFeedbackError::EnthalpyDomain {
                thermal_cell: CellId(11),
                ..
            })
        ));
        let scheduled = dark
            .clone()
            .with_power_schedule(
                &original,
                SpectralPowerSchedule::try_new(&[], 0).unwrap(),
                0,
                || Ok(()),
            )
            .unwrap();
        assert!(matches!(
            original.advance_observed(&scheduled, 0.5, config(0.5), || Ok(())),
            Err(OpticalFeedbackError::EnthalpyDomain {
                thermal_cell: CellId(11),
                ..
            })
        ));
        for incompatible in [&SHIFTED, &INTERNAL, &UNIDENTIFIED] {
            let bad = specimen(incompatible, &[(1.0, 0.0)], vec![]);
            let mut observer = dark.clone();
            assert!(matches!(
                observer.sources(&bad, &mut || Ok(())),
                Err(OpticalFeedbackError::Invalid(_))
            ));
        }
        let mut wrong = hbinding(0, 0);
        wrong.material_id = MaterialId(8);
        assert!(matches!(
            input(cx, &original, 0.0, &tables, &[wrong], budget()),
            Err(OpticalFeedbackError::Invalid(_))
        ));
        assert!(matches!(
            input(
                cx,
                &original,
                0.0,
                &tables,
                &bindings,
                OpticalFeedbackBudget {
                    max_total_temperature_knots: 1,
                    ..budget()
                },
            ),
            Err(OpticalFeedbackError::Budget)
        ));
        for work in [
            OpticalFeedbackBudget {
                max_tables: 1,
                ..budget()
            },
            OpticalFeedbackBudget {
                max_total_temperature_knots: 3,
                ..budget()
            },
        ] {
            assert!(matches!(
                nonlinear_input(cx, &original, work),
                Err(OpticalFeedbackError::Budget)
            ));
        }
        assert!(matches!(
            original.advance_observed(&dark, 0.1, config(0.1), || Err(DepositionError::Cancelled)),
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(dark.energy().thermal_input_j, 0.0);
        assert_eq!(scheduled.energy().thermal_input_j, 0.0);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 1.0);
    });
}
