use super::super::{CrystallizationCaloric, CrystallizationConfig};
use super::*;

#[derive(Clone, Copy)]
struct Phase {
    cp: f64,
    datum: f64,
    rho: f64,
    slope: f64,
}

impl CrystallizationCaloric for Phase {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        Ok(self.datum + self.cp * (t - 300.0))
    }

    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(self.cp)
    }

    fn density_kg_m3(&self, t: f64) -> Result<Option<f64>, CrystallizationError> {
        Ok(Some(self.rho + self.slope * (t - 300.0)))
    }

    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut identity = [0; 32];
        for (index, value) in [self.cp, self.datum, self.rho, self.slope]
            .into_iter()
            .enumerate()
        {
            identity[8 * index..8 * (index + 1)].copy_from_slice(&value.to_le_bytes());
        }
        Some(identity)
    }
}

fn phases() -> (Phase, Phase) {
    (
        Phase {
            cp: 1000.0,
            datum: 0.0,
            rho: 1000.0,
            slope: 0.0,
        },
        Phase {
            cp: 1000.0,
            datum: -100_000.0,
            rho: 1250.0,
            slope: 0.0,
        },
    )
}

fn config() -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [0x85; 32],
        avrami_exponent: 1.0,
        maximum_crystalline_mass_fraction: 0.8,
        temperature_range_k: [250.0, 650.0],
        melting_temperature_k: 500.0,
    }
}

fn rate(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.002 * (500.0 - t).max(0.0))
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}

#[test]
fn g1_adiabatic_crystallization_contracts_using_the_exact_original_adaptive_path() {
    let (parent, crystal) = phases();
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let isobaric = IsobaricCrystallization::new(&law, 100_000.0).unwrap();
    let initial = law.state_at_temperature(380.0, 0.0).unwrap();
    let controls = IsobaricCrystallizationControls::default();
    let adaptivity = CrystallizationAdaptivity {
        maximum_dt_s: 2.0,
        temperature_tolerance_k: 1e-5,
        crystalline_mass_fraction_tolerance: 1e-8,
        exposure_tolerance: 1e-8,
        ..CrystallizationAdaptivity::default()
    };
    let step = isobaric
        .advance_heat(initial, 2.0, 0.0, 2.0, controls, adaptivity, None, || false)
        .unwrap();
    let original = law
        .advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            2.0,
            controls.caloric,
            adaptivity,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(step.adaptive, Some(original));
    assert_eq!(step.state, original.state);
    // n=1 and fixed h: T=380+100X, exp(z)=3 exp(.08 t)-2.
    let expected_x = 0.8 * (1.0 - (3.0 * (0.08_f64 * 2.0).exp() - 2.0).recip());
    near(step.state.crystalline_mass_fraction(), expected_x, 3e-6);
    near(step.state.temperature_k(), 380.0 + 100.0 * expected_x, 3e-4);
    near(
        step.source_volume_change_m3,
        2.0 * step.state.crystalline_mass_fraction() * (1.0 / 1250.0 - 1.0 / 1000.0),
        1e-18,
    );
    assert_eq!(step.enthalpy_change_j, 0.0);
    assert!(step.volume_change_m3 < 0.0 && step.boundary_work_j < 0.0);
    assert!(step.internal_energy_change_j > 0.0);
    near(
        step.internal_energy_change_j,
        -isobaric.pressure_pa() * step.source_volume_change_m3,
        1e-12,
    );
    assert!(step.absolute_energy_residual_j <= controls.energy_tolerance_j);
    assert_eq!(step.duration_s, 2.0);
}

#[test]
fn g1_instantaneous_partial_and_full_melting_return_pressure_work_without_extra_heat() {
    let (parent, crystal) = phases();
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let model = law.equivalent_age_remelting().unwrap();
    let isobaric = IsobaricCrystallization::new(&law, 200_000.0).unwrap();
    let initial = law.state_at_temperature(500.0, f64::MAX).unwrap();
    let controls = IsobaricCrystallizationControls::default();
    assert!(matches!(
        isobaric.apply_heat(initial, 2.0, 120_000.0, controls, None, || false),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    let partial = isobaric
        .apply_heat(initial, 2.0, 120_000.0, controls, Some(model), || false)
        .unwrap();
    assert_eq!(partial.state.temperature_k(), 500.0);
    near(partial.state.crystalline_mass_fraction(), 0.2, 1e-13);
    near(partial.source_volume_change_m3, 0.00024, 1e-16);
    near(partial.boundary_work_j, 48.0, 1e-10);
    near(partial.internal_energy_change_j, 119_952.0, 1e-9);
    let full = isobaric
        .apply_heat(partial.state, 2.0, 140_000.0, controls, Some(model), || {
            false
        })
        .unwrap();
    assert_eq!(full.state.crystalline_mass_fraction(), 0.0);
    assert_eq!(full.state.exposure(), 0.0);
    near(full.state.temperature_k(), 550.0, 1e-9);
    near(full.source_volume_change_m3, 0.00008, 1e-16);
    near(full.boundary_work_j, 16.0, 1e-10);
    near(full.internal_energy_change_j, 139_984.0, 1e-9);
    for report in [partial, full] {
        assert_eq!(report.duration_s, 0.0);
        assert_eq!(report.adaptive, None);
        assert_eq!(report.enthalpy_change_j, report.supplied_heat_j);
    }
    let cooled = isobaric
        .apply_heat(
            partial.state,
            2.0,
            -200_000.0,
            controls,
            Some(model),
            || false,
        )
        .unwrap();
    near(cooled.state.temperature_k(), 400.0, 1e-9);
    assert_eq!(cooled.state.exposure(), partial.state.exposure());
    assert!(law.rate_s_inv(cooled.state.temperature_k()).unwrap() > 0.0);
    let held = isobaric
        .apply_heat(cooled.state, 2.0, 0.0, controls, Some(model), || false)
        .unwrap();
    assert_eq!(held.state, cooled.state);
    assert_eq!(held.volume_change_m3, 0.0);
    assert_eq!(held.enthalpy_change_j, 0.0);
}

#[test]
fn g1_original_temperature_dependent_densities_use_mass_fractions_and_reverse() {
    let (mut parent, mut crystal) = phases();
    parent.slope = -0.1;
    crystal.slope = -0.2;
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let isobaric = IsobaricCrystallization::new(&law, 300_000.0).unwrap();
    let initial = law
        .state_at_temperature(400.0, -(1.0_f64 - 0.3 / 0.8).ln())
        .unwrap();
    let controls = IsobaricCrystallizationControls::default();
    let heated = isobaric
        .apply_heat(initial, 3.0, 150_000.0, controls, None, || false)
        .unwrap();
    let x = initial.crystalline_mass_fraction();
    let expected_initial = 3.0 * ((1.0 - x) / 990.0 + x / 1230.0);
    let expected_final = 3.0 * ((1.0 - x) / 985.0 + x / 1220.0);
    near(heated.initial_volume_m3, expected_initial, 1e-17);
    near(heated.final_volume_m3, expected_final, 1e-16);
    near(
        heated.source_volume_change_m3,
        expected_final - expected_initial,
        1e-16,
    );
    assert!(heated.source_volume_change_m3 > 0.0);
    assert_eq!(heated.state.exposure(), initial.exposure());
    let cooled = isobaric
        .apply_heat(heated.state, 3.0, -150_000.0, controls, None, || false)
        .unwrap();
    near(cooled.state.temperature_k(), initial.temperature_k(), 1e-9);
    near(
        cooled.source_volume_change_m3,
        -heated.source_volume_change_m3,
        1e-16,
    );
    near(
        cooled.internal_energy_change_j,
        -heated.internal_energy_change_j,
        1e-9,
    );
    assert_eq!(cooled.state.exposure(), initial.exposure());
}

#[test]
fn g3_published_volume_rounding_is_bounded_even_at_zero_pressure_and_charged_to_work() {
    let parent = Phase {
        cp: 1.0,
        datum: 0.0,
        rho: 1.0,
        slope: f64::EPSILON,
    };
    let crystal = Phase {
        cp: 1.0,
        datum: -0.1,
        rho: 1.2,
        slope: 0.0,
    };
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(300.0, 0.0).unwrap();
    let controls = IsobaricCrystallizationControls {
        caloric: CrystallizationSolve {
            specific_enthalpy_tolerance_j_kg: 1e-14,
            temperature_tolerance_k: 1e-13,
            ..CrystallizationSolve::default()
        },
        energy_tolerance_j: 1e10,
        volume_tolerance_m3: 1e5,
    };
    let vacuum = IsobaricCrystallization::new(&law, 0.0).unwrap();
    let accepted = vacuum
        .apply_heat(initial, 1e20, 1e20, controls, None, || false)
        .unwrap();
    assert_eq!(accepted.state.temperature_k(), 301.0);
    let volume_defect = accepted.volume_change_m3 - accepted.source_volume_change_m3;
    assert!(volume_defect.abs() > 1000.0);
    assert_eq!(
        accepted.volume_change_m3,
        accepted.final_volume_m3 - accepted.initial_volume_m3
    );
    assert_eq!(
        vacuum.apply_heat(
            initial,
            1e20,
            1e20,
            IsobaricCrystallizationControls {
                volume_tolerance_m3: 1.0,
                ..controls
            },
            None,
            || false,
        ),
        Err(invalid(
            "isobaric volume representation exceeds its tolerance"
        ))
    );
    let loaded = IsobaricCrystallization::new(&law, 1e5).unwrap();
    assert_eq!(
        loaded.apply_heat(
            initial,
            1e20,
            1e20,
            IsobaricCrystallizationControls {
                energy_tolerance_j: 1.0,
                ..controls
            },
            None,
            || false,
        ),
        Err(invalid("isobaric energy residual exceeds its tolerance"))
    );
    let reported = loaded
        .apply_heat(initial, 1e20, 1e20, controls, None, || false)
        .unwrap();
    assert!(reported.absolute_energy_residual_j >= 1e5 * volume_defect.abs());
    assert!(reported.absolute_energy_residual_j >= reported.energy_residual_j.abs());
    assert_eq!(reported.boundary_work_j, 1e5 * reported.volume_change_m3);
}

struct WithoutDensity(Phase);
impl CrystallizationCaloric for WithoutDensity {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        self.0.specific_enthalpy_j_kg(t)
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.0.caloric_identity()
    }
}

#[test]
fn g0_original_source_admission_underflow_and_late_cancellation_refuse_atomically() {
    let (parent, crystal) = phases();
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let isobaric = IsobaricCrystallization::new(&law, 1e5).unwrap();
    let initial = law.state_at_temperature(400.0, 0.4).unwrap();
    let controls = IsobaricCrystallizationControls::default();
    for pressure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            IsobaricCrystallization::new(&law, pressure),
            Err(CrystallizationError::Invalid(_))
        ));
    }
    let missing = WithoutDensity(parent);
    let missing_law = CrystallizationLaw::new(&missing, &crystal, &rate, config()).unwrap();
    assert!(matches!(
        IsobaricCrystallization::new(&missing_law, 1e5),
        Err(CrystallizationError::Invalid(
            "original parent-phase density unavailable"
        ))
    ));
    let unidentified = |t| parent.specific_enthalpy_j_kg(t);
    let unidentified_law =
        CrystallizationLaw::new(&unidentified, &crystal, &rate, config()).unwrap();
    assert!(matches!(
        IsobaricCrystallization::new(&unidentified_law, 1e5),
        Err(CrystallizationError::Invalid(_))
    ));
    assert!(
        IsobaricCrystallizationControls {
            energy_tolerance_j: f64::NAN,
            ..controls
        }
        .validate()
        .is_err()
    );
    assert!(
        IsobaricCrystallizationControls {
            volume_tolerance_m3: -1.0,
            ..controls
        }
        .validate()
        .is_err()
    );
    assert_eq!(
        isobaric.apply_heat(initial, f64::MAX, f64::from_bits(1), controls, None, || {
            false
        }),
        Err(invalid("isobaric specific heat increment underflows"))
    );
    assert!(
        isobaric
            .volume_m3(initial, 0.0, controls.caloric, || false)
            .is_err()
    );
    let saturated = law.state_at_temperature(500.0, f64::MAX).unwrap();
    let subnormal_pressure = IsobaricCrystallization::new(&law, f64::from_bits(1)).unwrap();
    assert_eq!(
        subnormal_pressure.apply_heat(
            saturated,
            1.0,
            60_000.0,
            controls,
            Some(law.equivalent_age_remelting().unwrap()),
            || false
        ),
        Err(invalid("required isobaric product underflows"))
    );
    let mut complete_queries = 0;
    isobaric
        .apply_heat(initial, 1.0, 1000.0, controls, None, || {
            complete_queries += 1;
            false
        })
        .unwrap();
    assert!(complete_queries > 5);
    for cancellation_query in [1, complete_queries] {
        let mut queries = 0;
        assert_eq!(
            isobaric.apply_heat(initial, 1.0, 1000.0, controls, None, || {
                queries += 1;
                queries == cancellation_query
            }),
            Err(CrystallizationError::Cancelled)
        );
        assert_eq!(initial, law.state_at_temperature(400.0, 0.4).unwrap());
    }
}

#[test]
fn g3_pressure_work_and_internal_energy_changes_survive_a_large_common_enthalpy_datum() {
    let run = |offset| {
        let (mut parent, mut crystal) = phases();
        parent.datum += offset;
        crystal.datum += offset;
        let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
        let isobaric = IsobaricCrystallization::new(&law, 200_000.0).unwrap();
        let initial = law.state_at_temperature(500.0, f64::MAX).unwrap();
        isobaric
            .apply_heat(
                initial,
                2.0,
                120_000.0,
                IsobaricCrystallizationControls {
                    caloric: CrystallizationSolve {
                        specific_enthalpy_tolerance_j_kg: 0.001,
                        temperature_tolerance_k: 1e-6,
                        ..CrystallizationSolve::default()
                    },
                    energy_tolerance_j: 0.01,
                    volume_tolerance_m3: 1e-12,
                },
                Some(law.equivalent_age_remelting().unwrap()),
                || false,
            )
            .unwrap()
    };
    let original = run(0.0);
    let shifted = run(1e12);
    assert_eq!(original.enthalpy_change_j, 120_000.0);
    assert_eq!(shifted.enthalpy_change_j, original.enthalpy_change_j);
    near(
        shifted.state.crystalline_mass_fraction(),
        original.state.crystalline_mass_fraction(),
        4e-9,
    );
    near(
        shifted.source_volume_change_m3,
        original.source_volume_change_m3,
        2e-12,
    );
    near(shifted.boundary_work_j, original.boundary_work_j, 4e-7);
    near(
        shifted.internal_energy_change_j,
        original.internal_energy_change_j,
        4e-7,
    );
}
