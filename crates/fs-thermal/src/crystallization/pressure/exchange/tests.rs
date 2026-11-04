use super::super::super::{CrystallizationCaloric, CrystallizationConfig, CrystallizationLaw};
use super::*;

struct Phase {
    cp: f64,
    datum: f64,
    rho: f64,
}
impl CrystallizationCaloric for Phase {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        Ok(self.datum + self.cp * (t - 300.0))
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(self.cp)
    }
    fn density_kg_m3(&self, _: f64) -> Result<Option<f64>, CrystallizationError> {
        Ok(Some(self.rho))
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut identity = [0; 32];
        for (i, value) in [self.cp, self.datum, self.rho].into_iter().enumerate() {
            identity[8 * i..8 * (i + 1)].copy_from_slice(&value.to_le_bytes());
        }
        Some(identity)
    }
}
fn phases(cp: f64) -> (Phase, Phase) {
    (
        Phase {
            cp,
            datum: 0.0,
            rho: 1000.0,
        },
        Phase {
            cp,
            datum: -100_000.0,
            rho: 1250.0,
        },
    )
}
fn config(id: u8) -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [id; 32],
        avrami_exponent: 1.0,
        maximum_crystalline_mass_fraction: 0.8,
        temperature_range_k: [250.0, 650.0],
        melting_temperature_k: 500.0,
    }
}
fn inactive(_: f64) -> Result<f64, CrystallizationError> {
    Ok(0.0)
}
fn growing(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.02 * (500.0 - t).max(0.0))
}
fn parcel<'law, 'source>(
    law: &'law CrystallizationLaw<'source>,
    t: f64,
    age: f64,
    mass: f64,
    pressure: f64,
) -> IsobaricCrystallizationParcel<'law, 'source> {
    IsobaricCrystallizationParcel::new(
        law,
        law.state_at_temperature(t, age).unwrap(),
        mass,
        mass / 1000.0,
        pressure,
        IsobaricCrystallizationControls::default(),
        || false,
    )
    .unwrap()
}
fn controls(bracket: [f64; 2]) -> IsobaricExchangeControls {
    IsobaricExchangeControls {
        heat_bracket_j: bracket,
        maximum_evaluations: 96,
        residual_tolerance_j: 1e-5,
        energy_tolerance_j: 1e-6,
    }
}
fn side() -> IsobaricExchangeSideControls {
    IsobaricExchangeSideControls {
        caloric: IsobaricCrystallizationControls {
            energy_tolerance_j: 1e-6,
            ..IsobaricCrystallizationControls::default()
        },
        adaptivity: CrystallizationAdaptivity {
            temperature_tolerance_k: 1e-5,
            crystalline_mass_fraction_tolerance: 1e-8,
            exposure_tolerance: 1e-8,
            ..CrystallizationAdaptivity::default()
        },
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}
fn same_parcel(
    actual: &IsobaricCrystallizationParcel<'_, '_>,
    initial: &IsobaricCrystallizationParcel<'_, '_>,
) {
    assert_eq!(actual.state(), initial.state());
    assert_eq!(actual.volume_m3(), initial.volume_m3());
    assert_eq!(actual.time_s(), initial.time_s());
    assert_eq!(actual.mass_kg(), initial.mass_kg());
    assert_eq!(actual.reference_volume_m3(), initial.reference_volume_m3());
}

#[test]
fn g1_stiff_unlike_capacities_match_backward_euler_and_replay_bitwise() {
    let (lp, lc) = phases(1000.0);
    let (rp, rc) = phases(2000.0);
    let ll = CrystallizationLaw::new(&lp, &lc, &inactive, config(0x91)).unwrap();
    let rl = CrystallizationLaw::new(&rp, &rc, &inactive, config(0x92)).unwrap();
    let left = parcel(&ll, 450.0, 0.0, 2.0, 100_000.0);
    let right = parcel(&rl, 300.0, 0.0, 3.0, 200_000.0);
    let accepted_left = left.clone();
    let accepted_right = right.clone();
    let control = IsobaricExchangeControls {
        residual_tolerance_j: 1e-4,
        ..controls([0.0, 250_000.0])
    };
    let run = || exchange_heat(&left, &right, 1e6, 2.0, side(), side(), control, || false).unwrap();
    let (next_left, next_right, report) = run();
    let (replay_left, replay_right, replay) = run();
    let exact_q = 2e6 * 150.0 / (1.0 + 2e6 * (1.0 / 2000.0 + 1.0 / 6000.0));
    near(report.heat_left_to_right_j, exact_q, 1e-7);
    near(
        next_left.state().temperature_k(),
        450.0 - exact_q / 2000.0,
        1e-9,
    );
    near(
        next_right.state().temperature_k(),
        300.0 + exact_q / 6000.0,
        1e-9,
    );
    assert!(next_left.state().temperature_k() > next_right.state().temperature_k());
    assert!(report.residual_j.abs() <= control.residual_tolerance_j);
    assert!(report.evaluated_candidates <= control.maximum_evaluations);
    assert_eq!(report.left.supplied_heat_j, -report.right.supplied_heat_j);
    assert!(report.absolute_energy_residual_j <= control.energy_tolerance_j);
    assert_eq!(report.duration_s, 2.0);
    assert_eq!(next_left.time_s(), next_right.time_s());
    assert_eq!(report, replay);
    same_parcel(&next_left, &replay_left);
    same_parcel(&next_right, &replay_right);
    same_parcel(&left, &accepted_left);
    same_parcel(&right, &accepted_right);
}

#[test]
fn g1_recalescence_reverses_initial_temperature_order_and_zero_conductance_still_grows() {
    let (parent, crystal) = phases(1000.0);
    let ll = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0x93)).unwrap();
    let rl = CrystallizationLaw::new(&parent, &crystal, &growing, config(0x94)).unwrap();
    let left = parcel(&ll, 390.0, 0.0, 1.0, 1e5);
    let right = parcel(&rl, 380.0, 0.0, 1.0, 2e5);
    assert!(left.state().temperature_k() > right.state().temperature_k());
    let control = controls([-3000.0, 1000.0]);
    let (next_left, next_right, report) =
        exchange_heat(&left, &right, 20.0, 1.0, side(), side(), control, || false).unwrap();
    assert!(report.heat_left_to_right_j < 0.0);
    assert!(next_right.state().temperature_k() > next_left.state().temperature_k());
    assert!(next_right.state().crystalline_mass_fraction() > 0.0);
    assert!(report.right.boundary_work_j < 0.0);
    assert!(report.residual_j.abs() <= control.residual_tolerance_j);
    let replay = right
        .advance_heat(
            report.heat_left_to_right_j,
            1.0,
            side().caloric,
            side().adaptivity,
            || false,
        )
        .unwrap();
    same_parcel(&next_right, &replay.0);
    assert_eq!(report.right, replay.1);
    let (_, adiabatic, disconnected) = exchange_heat(
        &left,
        &right,
        0.0,
        1.0,
        side(),
        side(),
        IsobaricExchangeControls {
            maximum_evaluations: 1,
            ..control
        },
        || false,
    )
    .unwrap();
    assert_eq!(disconnected.heat_left_to_right_j, 0.0);
    assert_eq!(disconnected.residual_j, 0.0);
    assert_eq!(disconnected.evaluated_candidates, 1);
    let x = 0.8 * (1.0 - (3.0 * 0.8_f64.exp() - 2.0).recip());
    near(adiabatic.state().crystalline_mass_fraction(), x, 3e-6);
    near(adiabatic.state().temperature_k(), 380.0 + 100.0 * x, 3e-4);
}

#[test]
fn g1_contact_resolves_a_latent_plateau_and_external_pressure_work_once() {
    let (parent, crystal) = phases(1000.0);
    let law = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0x95)).unwrap();
    let left = parcel(&law, 550.0, 0.0, 1.0, 1e5);
    let right = parcel(&law, 500.0, f64::MAX, 1.0, 2e5)
        .with_equilibrium_remelting(law.equivalent_age_remelting().unwrap())
        .unwrap();
    let (next_left, next_right, report) = exchange_heat(
        &left,
        &right,
        100.0,
        1.0,
        side(),
        side(),
        controls([0.0, 10_000.0]),
        || false,
    )
    .unwrap();
    let exact_q = 5000.0 / 1.1;
    near(report.heat_left_to_right_j, exact_q, 1e-6);
    near(
        next_left.state().temperature_k(),
        550.0 - exact_q / 1000.0,
        1e-9,
    );
    assert_eq!(next_right.state().temperature_k(), 500.0);
    near(
        next_right.state().crystalline_mass_fraction(),
        0.8 - exact_q / 100_000.0,
        1e-10,
    );
    assert!(report.right.boundary_work_j > 0.0);
    near(
        report.right.boundary_work_j,
        2e5 * (0.0002 * exact_q / 100_000.0),
        1e-9,
    );
    near(report.combined_enthalpy_change_j, 0.0, 1e-7);
    near(
        report.combined_internal_energy_change_j,
        -report.combined_boundary_work_j,
        1e-7,
    );
    near(
        report.right.internal_energy_change_j + report.right.boundary_work_j,
        exact_q,
        1e-6,
    );
}

#[test]
fn g0_bracket_budget_source_endpoint_and_clock_failures_publish_neither_parcel() {
    let (parent, crystal) = phases(1000.0);
    let law = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0x96)).unwrap();
    let left = parcel(&law, 450.0, 0.0, 1.0, 1e5);
    let right = parcel(&law, 300.0, 0.0, 1.0, 2e5);
    let accepted_left = left.clone();
    let accepted_right = right.clone();
    let run = |control| exchange_heat(&left, &right, 100.0, 1.0, side(), side(), control, || false);
    assert!(matches!(
        run(controls([0.0, 1.0])),
        Err(CrystallizationError::Invalid(
            "isobaric exchange heat bracket does not enclose a root"
        ))
    ));
    assert!(matches!(
        run(IsobaricExchangeControls {
            maximum_evaluations: 1,
            ..controls([0.0, 20_000.0])
        }),
        Err(CrystallizationError::NonConverged)
    ));
    assert!(matches!(
        run(controls([0.0, 1e9])),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    assert!(controls([1.0, -1.0]).validate().is_err());
    assert!(
        IsobaricExchangeControls {
            maximum_evaluations: 0,
            ..controls([0.0, 1.0])
        }
        .validate()
        .is_err()
    );
    let advanced = left
        .advance_heat(0.0, 0.1, side().caloric, side().adaptivity, || false)
        .unwrap()
        .0;
    assert!(matches!(
        exchange_heat(
            &advanced,
            &right,
            100.0,
            1.0,
            side(),
            side(),
            controls([0.0, 20_000.0]),
            || false
        ),
        Err(CrystallizationError::Invalid(
            "isobaric exchange requires positive duration and equal parcel clocks"
        ))
    ));
    same_parcel(&left, &accepted_left);
    same_parcel(&right, &accepted_right);
}

#[test]
fn g4_original_source_errors_final_cancellation_and_pair_energy_budget_remain_distinct() {
    let (parent, crystal) = phases(1000.0);
    let failure = std::cell::Cell::new(false);
    let source = |_| {
        if failure.get() {
            Err(invalid("fixture kinetic source unavailable"))
        } else {
            Ok(0.0)
        }
    };
    let ll = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0x97)).unwrap();
    let rl = CrystallizationLaw::new(&parent, &crystal, &source, config(0x98)).unwrap();
    let left = parcel(&ll, 450.0, 0.0, 1.0, 1e5);
    let right = parcel(&rl, 300.0, 0.0, 1.0, 2e5);
    let accepted_left = left.clone();
    let accepted_right = right.clone();
    let control = controls([0.0, 1000.0]);
    failure.set(true);
    assert!(matches!(
        exchange_heat(&left, &right, 7.0, 0.3, side(), side(), control, || false),
        Err(CrystallizationError::Invalid(
            "fixture kinetic source unavailable"
        ))
    ));
    failure.set(false);
    let mut complete_queries = 0;
    let (_, _, report) = exchange_heat(&left, &right, 7.0, 0.3, side(), side(), control, || {
        complete_queries += 1;
        false
    })
    .unwrap();
    assert!(report.absolute_energy_residual_j > 0.0);
    assert!(report.residual_j.abs() <= control.residual_tolerance_j);
    assert!(matches!(
        exchange_heat(
            &left,
            &right,
            7.0,
            0.3,
            side(),
            side(),
            IsobaricExchangeControls {
                energy_tolerance_j: 0.0,
                ..control
            },
            || false
        ),
        Err(CrystallizationError::Invalid(
            "combined isobaric exchange energy residual exceeds its tolerance"
        ))
    ));
    assert!(complete_queries > 10);
    for stop in [1, complete_queries] {
        let mut queries = 0;
        assert!(matches!(
            exchange_heat(&left, &right, 7.0, 0.3, side(), side(), control, || {
                queries += 1;
                queries == stop
            }),
            Err(CrystallizationError::Cancelled)
        ));
    }
    same_parcel(&left, &accepted_left);
    same_parcel(&right, &accepted_right);
}

#[test]
fn g0_contact_uses_the_shared_represented_clock_increment_and_refuses_lost_time() {
    let (parent, crystal) = phases(1000.0);
    let law = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0x99)).unwrap();
    let long = CrystallizationAdaptivity {
        maximum_dt_s: 1e16,
        ..side().adaptivity
    };
    let left = parcel(&law, 450.0, 0.0, 1.0, 1e5)
        .advance_heat(0.0, 1e16, side().caloric, long, || false)
        .unwrap()
        .0;
    let right = parcel(&law, 300.0, 0.0, 1.0, 2e5)
        .advance_heat(0.0, 1e16, side().caloric, long, || false)
        .unwrap()
        .0;
    assert!(matches!(
        exchange_heat(
            &left,
            &right,
            100.0,
            1.0,
            side(),
            side(),
            controls([0.0, 50_000.0]),
            || false
        ),
        Err(CrystallizationError::Invalid(
            "isobaric exchange elapsed-time increment is unrepresentable"
        ))
    ));
    let (next_left, next_right, report) = exchange_heat(
        &left,
        &right,
        100.0,
        3.0,
        side(),
        side(),
        controls([0.0, 50_000.0]),
        || false,
    )
    .unwrap();
    assert_eq!(report.duration_s, 4.0);
    assert_eq!(report.left.duration_s, 4.0);
    assert_eq!(report.right.duration_s, 4.0);
    assert_eq!(next_left.time_s(), next_right.time_s());
    assert_eq!(next_left.time_s() - left.time_s(), report.duration_s);
    let exact_q = 400.0 * 150.0 / (1.0 + 400.0 * (1.0 / 1000.0 + 1.0 / 1000.0));
    near(report.heat_left_to_right_j, exact_q, 1e-6);
}

#[test]
fn g1_original_endpoint_conductance_matches_an_independent_radiative_cubic_oracle() {
    let (parent, crystal) = phases(1000.0);
    let ll = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0xa1)).unwrap();
    let rl = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0xa2)).unwrap();
    let left = parcel(&ll, 450.0, 0.0, 1.0, 1e5);
    let right = parcel(&rl, 300.0, 0.0, 1.0, 2e5);
    // Authored radiation-shaped fixture: qdot=alpha*(Tl^4-Tr^4).
    // Its endpoint conductance is alpha*(Tl+Tr)*(Tl^2+Tr^2).
    let alpha = 5e-8;
    let dt = 10.0;
    let mut source_queries = 0;
    let mut lowest_left: f64 = 450.0;
    let mut highest_left: f64 = 250.0;
    let control = IsobaricExchangeControls {
        residual_tolerance_j: 1e-7,
        ..controls([0.0, 20_000.0])
    };
    let (next_left, next_right, report) = exchange_heat_with_conductance(
        &left,
        &right,
        |l, r| {
            source_queries += 1;
            assert_eq!(l.time_s(), dt);
            assert_eq!(r.time_s(), dt);
            assert_eq!(l.state().crystalline_mass_fraction(), 0.0);
            assert_eq!(r.state().crystalline_mass_fraction(), 0.0);
            assert_eq!(l.volume_m3(), 0.001);
            assert_eq!(r.volume_m3(), 0.001);
            let tl = l.state().temperature_k();
            let tr = r.state().temperature_k();
            lowest_left = lowest_left.min(tl);
            highest_left = highest_left.max(tl);
            Ok(alpha * (tl + tr) * (tl * tl + tr * tr))
        },
        dt,
        side(),
        side(),
        control,
        || false,
    )
    .unwrap();
    // Equal capacities preserve M=(Tl+Tr)/2=375. With d=(Tl-Tr)/2,
    // Q=C*(75-d) gives a*d^3+(C+a*M^2)*d=C*75, a=8*alpha*dt*M.
    // This analytic hyperbolic form solves the depressed cubic independently
    // of the production bracket and its factored endpoint conductance.
    let capacity = 1000.0;
    let mean = 375.0;
    let a = 8.0 * alpha * dt * mean;
    let p = (capacity + a * mean * mean) / a;
    let r = capacity * 75.0 / a;
    let root_scale = (p / 3.0).sqrt();
    let d = 2.0 * root_scale * ((r / (2.0 * root_scale.powi(3))).asinh() / 3.0).sinh();
    let exact_q = capacity * (75.0 - d);
    near(report.heat_left_to_right_j, exact_q, 1e-6);
    near(next_left.state().temperature_k(), mean + d, 1e-9);
    near(next_right.state().temperature_k(), mean - d, 1e-9);
    assert_eq!(source_queries, report.evaluated_candidates);
    assert!(lowest_left < highest_left);
    let tl = next_left.state().temperature_k();
    let tr = next_right.state().temperature_k();
    assert_eq!(
        report.conductance_w_k,
        alpha * (tl + tr) * (tl * tl + tr * tr)
    );
    near(
        report.heat_left_to_right_j,
        alpha * dt * (tl.powi(4) - tr.powi(4)),
        1e-7,
    );
    let frozen_g = alpha * (450.0 + 300.0) * (450.0 * 450.0 + 300.0 * 300.0);
    let frozen = exchange_heat(&left, &right, frozen_g, dt, side(), side(), control, || {
        false
    })
    .unwrap();
    assert!((frozen.2.heat_left_to_right_j - report.heat_left_to_right_j).abs() > 1.0);
    assert!(report.conductance_w_k < frozen_g);
    assert!(report.absolute_energy_residual_j <= control.energy_tolerance_j);
}

#[test]
fn g4_endpoint_conductance_source_failures_and_final_cancellation_publish_neither_body() {
    let (parent, crystal) = phases(1000.0);
    let law = CrystallizationLaw::new(&parent, &crystal, &inactive, config(0xa3)).unwrap();
    let left = parcel(&law, 450.0, 0.0, 1.0, 1e5);
    let right = parcel(&law, 300.0, 0.0, 1.0, 2e5);
    let accepted_left = left.clone();
    let accepted_right = right.clone();
    let control = controls([0.0, 1000.0]);
    let source_error = CrystallizationError::TemperatureOutOfRange {
        temperature_k: 450.0,
        minimum_k: 250.0,
        maximum_k: 440.0,
    };
    let mut queries = 0;
    let failure = exchange_heat_with_conductance(
        &left,
        &right,
        |_, _| {
            queries += 1;
            if queries == 2 {
                Err(source_error)
            } else {
                Ok(1.0)
            }
        },
        1.0,
        side(),
        side(),
        control,
        || false,
    );
    assert_eq!(failure.err(), Some(source_error));
    assert_eq!(queries, 2);
    // Validation is per proposed endpoint, not only at the first source query.
    for invalid_g in [f64::NAN, f64::INFINITY, -1.0] {
        let mut queries = 0;
        let failure = exchange_heat_with_conductance(
            &left,
            &right,
            |_, _| {
                queries += 1;
                Ok(if queries == 2 { invalid_g } else { 1.0 })
            },
            1.0,
            side(),
            side(),
            control,
            || false,
        );
        assert_eq!(
            failure.err(),
            Some(invalid(
                "isobaric exchange conductance must be finite and nonnegative"
            ))
        );
        assert_eq!(queries, 2);
    }
    let cancelled_by_source = std::cell::Cell::new(false);
    let failure = exchange_heat_with_conductance(
        &left,
        &right,
        |_, _| {
            cancelled_by_source.set(true);
            Ok(1.0)
        },
        1.0,
        side(),
        side(),
        control,
        || cancelled_by_source.get(),
    );
    assert_eq!(failure.err(), Some(CrystallizationError::Cancelled));
    let mut complete_checks = 0;
    let with_source = exchange_heat_with_conductance(
        &left,
        &right,
        |_, _| Ok(1.0),
        1.0,
        side(),
        side(),
        control,
        || {
            complete_checks += 1;
            false
        },
    )
    .unwrap();
    let with_fixed =
        exchange_heat(&left, &right, 1.0, 1.0, side(), side(), control, || false).unwrap();
    assert_eq!(with_source.2, with_fixed.2);
    same_parcel(&with_source.0, &with_fixed.0);
    same_parcel(&with_source.1, &with_fixed.1);
    let mut checks = 0;
    let failure = exchange_heat_with_conductance(
        &left,
        &right,
        |_, _| Ok(1.0),
        1.0,
        side(),
        side(),
        control,
        || {
            checks += 1;
            checks == complete_checks
        },
    );
    assert_eq!(failure.err(), Some(CrystallizationError::Cancelled));
    same_parcel(&left, &accepted_left);
    same_parcel(&right, &accepted_right);
}
