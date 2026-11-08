use super::*;
use crate::thermal::transport::expansion::{FreeExpansionSourceQuery, FreeExpansionStrainSources};
use fs_thermal::crystallization::{
    CrystallizationAdaptivity,
    pressure::{IsobaricCrystallization, IsobaricCrystallizationStep},
};

fn density_key(phase: usize) -> PropertyKey {
    PropertyKey::new(
        ["parent_density", "crystal_density"][phase],
        Dims([-3, 1, 0, 0, 0, 0]),
    )
}

fn strain_key(phase: usize) -> PropertyKey {
    PropertyKey::new(
        ["parent_free_length", "crystal_free_length"][phase],
        Dims::NONE,
    )
}

fn fixture(pressure_qualified: bool) -> Fixture {
    let mut fixture = Fixture::standard();
    let mut claims = fixture.card.claims().clone();
    for phase in 0..2 {
        for (key, value, interpolation, range) in [
            (
                density_key(phase),
                PropertyValue::Scalar {
                    value: [1000.0, 1250.0][phase],
                    dims: Dims([-3, 1, 0, 0, 0, 0]),
                },
                InterpolationPolicy::ConstantWithinValidity,
                [400.0, 400.0],
            ),
            (
                strain_key(phase),
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, [-0.01, -0.005][phase]),
                        (500.0, [0.01, 0.005][phase]),
                    ],
                    dims: Dims::NONE,
                },
                InterpolationPolicy::LinearInside,
                [300.0, 500.0],
            ),
        ] {
            let mut validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                range[0],
                range[1],
            );
            if pressure_qualified {
                validity = validity.with_quantity(
                    "p",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    100000.0,
                    200000.0,
                );
            }
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value,
                    validity,
                    interpolation,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic matched phase density and free-length oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    fixture.card = MaterialCard::assemble(fixture.card.id().clone(), claims, vec![]).unwrap();
    fixture
}

fn density(fixture: &Fixture, phase: usize, pressure: f64) -> PreparedFreeExpansionDensity {
    PreparedFreeExpansionDensity::try_new(
        &fixture.card,
        FreeExpansionSourceQuery {
            density: &density_key(phase),
            density_selection: SELECTION,
            reference: &point(400.0, pressure),
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &strain_key(phase),
                selection: SELECTION,
            },
            temperature_axis: "T",
            lower: &point(300.0, pressure),
            upper: &point(500.0, pressure),
        },
    )
    .unwrap()
}

fn expanding_pair<'a>(
    fixture: &'a Fixture,
    parent: &'a IsobaricHeatCapacityLaw,
    crystal: &'a IsobaricHeatCapacityLaw,
    parent_density: &'a PreparedFreeExpansionDensity,
    crystal_density: &'a PreparedFreeExpansionDensity,
) -> Result<IsobaricCrystallizationPair<'a>, ThermalTransportError> {
    IsobaricCrystallizationPair::try_new_with_free_expansion(
        &fixture.card,
        IsobaricCrystallizationQuery {
            parent,
            crystal,
            enthalpy_gap: &fixture.gap_key,
            reference: &fixture.reference,
            selection: SELECTION,
            temperature_range_k: [300.0, 500.0],
        },
        parent_density,
        crystal_density,
    )
}

fn controls() -> IsobaricCrystallizationControls {
    IsobaricCrystallizationControls {
        caloric: CrystallizationSolve {
            maximum_iterations: 128,
            specific_enthalpy_tolerance_j_kg: 1e-7,
            temperature_tolerance_k: 1e-10,
        },
        energy_tolerance_j: 1e-6,
        volume_tolerance_m3: 1e-15,
    }
}

fn rate(_: f64) -> Result<f64, CrystallizationError> {
    Ok(0.7)
}

// Independent analytic specific volumes from the two authored engineering
// length laws, with their common 400 K density reference. Fractions are masses.
fn volume(mass: f64, temperature: f64, fraction: f64) -> f64 {
    let parent = (1.0 + 1e-4 * (temperature - 400.0)).powi(3) / 1000.0;
    let crystal = (1.0 + 5e-5 * (temperature - 400.0)).powi(3) / 1250.0;
    mass * ((1.0 - fraction) * parent + fraction * crystal)
}

fn assert_account(step: IsobaricCrystallizationStep, initial_volume: f64, final_volume: f64) {
    let change = final_volume - initial_volume;
    close(step.initial_volume_m3, initial_volume, 3e-15);
    close(step.final_volume_m3, final_volume, 3e-15);
    close(step.volume_change_m3, change, 4e-15);
    close(step.source_volume_change_m3, change, 4e-15);
    close(step.boundary_work_j, 100000.0 * change, 4e-10);
    close(step.enthalpy_change_j, step.supplied_heat_j, 1e-6);
    close(
        step.internal_energy_change_j + step.boundary_work_j,
        step.supplied_heat_j,
        1e-6,
    );
    assert!(step.absolute_energy_residual_j <= 1e-6);
}

#[test]
fn source_density_uses_rebased_phase_temperature_without_changing_exact_cp() {
    let fixture = fixture(true);
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 450.0);
    let parent_density = density(&fixture, 0, 100000.0);
    let crystal_density = density(&fixture, 1, 100000.0);
    let plain = fixture.pair(&parent, &crystal).unwrap();
    let pair = expanding_pair(
        &fixture,
        &parent,
        &crystal,
        &parent_density,
        &crystal_density,
    )
    .unwrap();
    assert_ne!(pair.identity(), plain.identity());
    assert_eq!(pair.pressure_pa(), 100000.0);
    for (phase, (branch, original)) in [
        (pair.parent(), plain.parent()),
        (pair.crystal(), plain.crystal()),
    ]
    .into_iter()
    .enumerate()
    {
        assert_ne!(branch.identity(), original.identity());
        assert_eq!(
            ThermalLaw::caloric_identity(branch),
            CrystallizationCaloric::caloric_identity(branch)
        );
        assert!(branch.free_expansion().is_some());
        for temperature in [300.0, 350.0, 400.0, 450.0, 500.0] {
            let h = branch.enthalpy_at_temperature(temperature).unwrap();
            assert_eq!(h, original.enthalpy_at_temperature(temperature).unwrap());
            let expected = volume(1.0, temperature, phase as f64).recip();
            close(
                branch.density_at_temperature(temperature).unwrap().unwrap(),
                expected,
                2e-12,
            );
            close(
                CrystallizationCaloric::density_kg_m3(branch, temperature)
                    .unwrap()
                    .unwrap(),
                expected,
                2e-12,
            );
            close(
                branch.bulk_density_kg_m3(h).unwrap().unwrap(),
                expected,
                2e-12,
            );
        }
        assert!(branch.density_at_temperature(299.0).is_err());
        assert!(branch.bulk_density_kg_m3(f64::NAN).is_err());
    }
    let law = plain
        .law(&rate, ContentHash([3; 32]), 1.0, 1.0, 550.0)
        .unwrap();
    assert!(IsobaricCrystallization::new(&law, plain.pressure_pa()).is_err());
}

#[test]
fn native_phase_growth_and_heat_conserve_mass_volume_pressure_work_and_history() {
    let fixture = fixture(true);
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 450.0);
    let parent_density = density(&fixture, 0, 100000.0);
    let crystal_density = density(&fixture, 1, 100000.0);
    let pair = expanding_pair(
        &fixture,
        &parent,
        &crystal,
        &parent_density,
        &crystal_density,
    )
    .unwrap();
    let law = pair
        .law(&rate, ContentHash([3; 32]), 1.0, 1.0, 550.0)
        .unwrap();
    let state = law.state_at_temperature(400.0, 0.0).unwrap();
    let mass = 2.5;
    // A distinct mechanical reference volume must not replace m*v(T,X).
    let initial = pair
        .isobaric_parcel(&law, state, mass, 0.004, controls(), || false)
        .unwrap();
    close(initial.volume_m3(), volume(mass, 400.0, 0.0), 0.0);
    let (grown, growth) = initial
        .advance_heat(
            0.0,
            1.0,
            controls(),
            CrystallizationAdaptivity::default(),
            || false,
        )
        .unwrap();
    let fraction = -(-0.7_f64).exp_m1();
    // The original affine Cp laws give h(T,X)=(1.5-0.5X)*T²-240000.
    let temperature = (240000.0 / (1.5 - 0.5 * fraction)).sqrt();
    close(grown.state().crystalline_mass_fraction(), fraction, 3e-13);
    close(grown.state().temperature_k(), temperature, 2e-9);
    assert_account(
        growth,
        initial.volume_m3(),
        volume(mass, temperature, fraction),
    );
    assert!(growth.boundary_work_j < 0.0);
    assert!(growth.internal_energy_change_j > 0.0);

    let heat = 20000.0;
    let (heated, heating) = grown.apply_heat(heat, controls(), || false).unwrap();
    let heated_temperature = ((240000.0 + heat / mass) / (1.5 - 0.5 * fraction)).sqrt();
    close(heated.state().temperature_k(), heated_temperature, 2e-9);
    assert_eq!(
        heated.state().crystalline_mass_fraction(),
        grown.state().crystalline_mass_fraction()
    );
    assert_eq!(heated.time_s(), grown.time_s());
    assert_account(
        heating,
        grown.volume_m3(),
        volume(mass, heated_temperature, fraction),
    );
    assert!(heating.boundary_work_j > 0.0);

    let (continued, continuation) = heated
        .advance_heat(
            0.0,
            0.5,
            controls(),
            CrystallizationAdaptivity::default(),
            || false,
        )
        .unwrap();
    let final_fraction = -(-1.05_f64).exp_m1();
    let final_temperature = ((240000.0 + heat / mass) / (1.5 - 0.5 * final_fraction)).sqrt();
    close(
        continued.state().crystalline_mass_fraction(),
        final_fraction,
        3e-13,
    );
    close(continued.state().temperature_k(), final_temperature, 2e-9);
    assert_account(
        continuation,
        heated.volume_m3(),
        volume(mass, final_temperature, final_fraction),
    );
    assert_eq!(continued.time_s(), 1.5);
    for parcel in [&initial, &grown, &heated, &continued] {
        assert_eq!(parcel.mass_kg(), mass);
        assert_eq!(parcel.reference_volume_m3(), 0.004);
        assert_eq!(parcel.pressure_pa(), pair.pressure_pa());
    }
    close(
        continued.internal_energy_j().unwrap() - initial.internal_energy_j().unwrap()
            + growth.boundary_work_j
            + heating.boundary_work_j
            + continuation.boundary_work_j,
        heat,
        2e-6,
    );
    assert_eq!(initial.state(), state);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn density_binding_refuses_another_pressure_or_phase_history_and_late_cancellation_is_atomic() {
    let fixture = fixture(true);
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 450.0);
    let parent_density = density(&fixture, 0, 100000.0);
    let crystal_density = density(&fixture, 1, 100000.0);
    let wrong_pressure = density(&fixture, 0, 200000.0);
    assert!(
        expanding_pair(
            &fixture,
            &parent,
            &crystal,
            &wrong_pressure,
            &crystal_density
        )
        .is_err()
    );
    let pair = expanding_pair(
        &fixture,
        &parent,
        &crystal,
        &parent_density,
        &crystal_density,
    )
    .unwrap();
    let plain = fixture.pair(&parent, &crystal).unwrap();
    let plain_law = plain
        .law(&rate, ContentHash([3; 32]), 1.0, 1.0, 550.0)
        .unwrap();
    let plain_state = plain_law.state_at_temperature(400.0, 0.0).unwrap();
    assert!(
        pair.isobaric_parcel(&plain_law, plain_state, 2.5, 0.004, controls(), || false)
            .is_err()
    );
    let swapped = expanding_pair(
        &fixture,
        &parent,
        &crystal,
        &crystal_density,
        &parent_density,
    )
    .unwrap();
    let swapped_law = swapped
        .law(&rate, ContentHash([3; 32]), 1.0, 1.0, 550.0)
        .unwrap();
    assert_ne!(pair.identity(), swapped.identity());
    let swapped_state = swapped_law.state_at_temperature(400.0, 0.0).unwrap();
    assert!(
        pair.isobaric_parcel(&swapped_law, swapped_state, 2.5, 0.004, controls(), || {
            false
        })
        .is_err()
    );

    let law = pair
        .law(&rate, ContentHash([3; 32]), 1.0, 1.0, 550.0)
        .unwrap();
    let state = law.state_at_temperature(400.0, 0.3).unwrap();
    let parcel = pair
        .isobaric_parcel(&law, state, 2.5, 0.004, controls(), || false)
        .unwrap();
    let original_volume = parcel.volume_m3();
    let original_energy = parcel.internal_energy_j().unwrap();
    let mut visits = 0;
    parcel
        .apply_heat(1000.0, controls(), || {
            visits += 1;
            false
        })
        .unwrap();
    assert!(visits > 3);
    let mut cancelled_visits = 0;
    let cancelled = parcel.apply_heat(1000.0, controls(), || {
        cancelled_visits += 1;
        cancelled_visits == visits
    });
    assert!(matches!(cancelled, Err(CrystallizationError::Cancelled)));
    assert_eq!(parcel.state(), state);
    assert_eq!(parcel.volume_m3(), original_volume);
    assert_eq!(parcel.internal_energy_j().unwrap(), original_energy);
    assert_eq!(parcel.mass_kg(), 2.5);
    assert_eq!(parcel.reference_volume_m3(), 0.004);
    assert_eq!(parcel.time_s(), 0.0);
}

#[test]
fn phase_density_requires_explicit_source_pressure_support() {
    let fixture = fixture(false);
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 450.0);
    let parent_density = density(&fixture, 0, 100000.0);
    let crystal_density = density(&fixture, 1, 100000.0);
    assert!(
        expanding_pair(
            &fixture,
            &parent,
            &crystal,
            &parent_density,
            &crystal_density
        )
        .is_err()
    );
}
