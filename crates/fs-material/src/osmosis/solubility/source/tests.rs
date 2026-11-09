//! Source admission exercised through the original finite material owners.

use super::*;
use crate::osmosis::{
    AqueousSolution,
    solubility::{FiniteSolubilitySolution, kinetics::KineticSolubilitySolution},
};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, LawId, LawParameter, MaterialStateId, Provenance};
use fs_qty::{Time, Volume};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::kinetics::SolubilityKineticStateV1,
};

fn species(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn domain(solvent: &str, solute: &str, low: f64, high: f64) -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            TEMPERATURE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            low,
            high,
        )
        .with_quantity(
            format!("source-solvent:{solvent}"),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        )
        .with_quantity(
            format!("source-solute:{solute}"),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        )
}

fn model(law: &str, parameters: &[(&str, f64, Dims)], bounds: [f64; 2]) -> ConstitutiveModelCard {
    let equilibrium = law == NEUTRAL_SOLUBILITY_LAW;
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        parameters: parameters
            .iter()
            .map(|&(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        state_schema_version: u32::from(!equilibrium),
        initial_state: if equilibrium {
            InitialStatePolicy::ZeroInternalState
        } else {
            InitialStatePolicy::RequiresDeclaredState
        },
        validity: domain("water", "A", bounds[0], bounds[1]),
        sources: vec![hash_bytes(
            b"synthetic neutral-solubility integration fixture",
        )],
        provenance: Provenance {
            source: "synthetic law fixture, not measured material calibration".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn equilibrium_model(enthalpy: f64) -> ConstitutiveModelCard {
    model(
        NEUTRAL_SOLUBILITY_LAW,
        &[
            (REFERENCE_TEMPERATURE_PARAMETER, 300.0, Temperature::DIMS),
            (
                REFERENCE_SOLUBILITY_PARAMETER,
                100.0,
                AmountConcentration::DIMS,
            ),
            (
                DISSOLUTION_ENTHALPY_PARAMETER,
                enthalpy,
                MolarEnergyQuantityV1::DIMS,
            ),
        ],
        [280.0, 350.0],
    )
}

fn particle_model() -> ConstitutiveModelCard {
    model(
        PARTICLE_FILM_KINETICS_LAW,
        &[
            (
                TRANSFER_CONDUCTANCE_PARAMETER,
                1e-3,
                VolumetricFlowRate::DIMS,
            ),
            (REFERENCE_SOLID_PARAMETER, 0.1, Amount::DIMS),
        ],
        [290.0, 330.0],
    )
}

fn fixed_model(policy: f64) -> ConstitutiveModelCard {
    model(
        FIXED_FILM_KINETICS_LAW,
        &[
            (
                TRANSFER_CONDUCTANCE_PARAMETER,
                1e-3,
                VolumetricFlowRate::DIMS,
            ),
            (PRECIPITATION_POLICY_PARAMETER, policy, Dims::NONE),
        ],
        [290.0, 330.0],
    )
}

fn card(models: Vec<ConstitutiveModelCard>) -> Arc<MaterialCard> {
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "explicit water/A pair".into(),
                phase: "neutral aqueous and pure solid".into(),
                process: "prescribed synthetic bath and film".into(),
                revision: 0,
            },
            ClaimSet::new(),
            models,
        )
        .unwrap(),
    )
}

fn admit(
    source: ConstitutiveModelCard,
) -> Result<NeutralSolubilitySource, NeutralSolubilitySourceError> {
    let identity = source.content_hash();
    NeutralSolubilitySource::from_material(
        card(vec![source]),
        identity,
        species("water"),
        species("A"),
    )
}

fn bind(
    equilibrium: &NeutralSolubilitySource,
    source: ConstitutiveModelCard,
) -> Result<NeutralSolubilityKineticsSource, NeutralSolubilitySourceError> {
    let identity = source.content_hash();
    equilibrium.with_kinetics(card(vec![source]), identity)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_source_drives_temperature_drying_and_redissolution() {
    let original = equilibrium_model(12_000.0);
    let material = card(vec![original.clone()]);
    let source = NeutralSolubilitySource::from_material(
        material.clone(),
        original.content_hash(),
        species("water"),
        species("A"),
    )
    .unwrap();
    assert!(Arc::ptr_eq(source.material(), &material));
    assert_eq!(source.model_card(), &original);
    assert_eq!(source.law().phase_heat_capacities(), None);
    let solution = FiniteSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(source.law().clone(), Amount::new(0.15))],
    )
    .unwrap();
    near(solution.phase_partitions()[0].dissolved.value(), 0.1);
    near(solution.phase_partitions()[0].solid.value(), 0.05);
    let dried = solution.with_solvent_volume_m3(5e-4).unwrap();
    near(dried.phase_partitions()[0].dissolved.value(), 0.05);
    near(dried.phase_partitions()[0].solid.value(), 0.1);
    let diluted = dried.with_solvent_volume_m3(2e-3).unwrap();
    near(diluted.phase_partitions()[0].dissolved.value(), 0.15);
    assert_eq!(diluted.phase_partitions()[0].solid.value(), 0.0);
    near(
        diluted.solution().osmotic_pressure_pa().unwrap(),
        R * 300.0 * 75.0,
    );
    let warmed = solution.with_temperature(Temperature::new(310.0)).unwrap();
    let expected = 0.1 * (12_000.0 / R * (1.0 / 300.0 - 1.0 / 310.0)).exp();
    near(warmed.phase_partitions()[0].dissolved.value(), expected);
    near(warmed.phase_partitions()[0].total_solute.value(), 0.15);
    assert!(solution.with_temperature(Temperature::new(351.0)).is_err());
}

#[test]
fn original_particle_reference_survives_native_dissolution_and_volume_change() {
    let source = admit(equilibrium_model(0.0)).unwrap();
    let original = particle_model();
    let kinetic = bind(&source, original.clone()).unwrap();
    assert_eq!(kinetic.model_card(), &original);
    assert_eq!(kinetic.equilibrium_source(), &source);
    let mut solution = KineticSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(
            kinetic.kinetics().clone(),
            SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.1)).unwrap(),
        )],
    )
    .unwrap();
    solution.advance(Time::new(1.0), || false).unwrap();
    // For total = V*c_sat, exact integration of ds/dt = -G_ref/V *
    // (s/s_ref)^(2/3) * s gives this independent finite-bath solution.
    let solid = 0.1 / (1.0_f64 + 2.0 / 3.0).powf(1.5);
    near(solution.last_steps()[0].state().solid().value(), solid);
    near(solution.solution().solute_moles()[0], 0.1 - solid);
    near(
        solution.models()[0]
            .conductance_at(Amount::new(solid))
            .unwrap()
            .value(),
        1e-3 * (solid / 0.1).powf(2.0 / 3.0),
    );
    let mut dried = solution.with_solvent_volume_m3(2.5e-4).unwrap();
    dried.advance(Time::new(0.2), || false).unwrap();
    let grown = dried.last_steps()[0].state();
    assert!(grown.solid().value() > solid);
    near(grown.dissolved().value() + grown.solid().value(), 0.1);
    assert_eq!(
        dried.models()[0].particle_reference_solid(),
        Some(Amount::new(0.1))
    );
    near(
        dried.solution().solute_moles()[0],
        grown.dissolved().value(),
    );
    assert!(dried.with_temperature(Temperature::new(340.0)).is_err());
    assert_eq!(kinetic.model_card(), &original);
}

#[test]
fn source_policy_controls_precipitation_without_inventing_a_surface() {
    let equilibrium = admit(equilibrium_model(0.0)).unwrap();
    let existing = bind(&equilibrium, fixed_model(0.0)).unwrap();
    let seeded = bind(&equilibrium, fixed_model(1.0)).unwrap();
    let particle = bind(&equilibrium, particle_model()).unwrap();
    let initial =
        SolubilityKineticStateV1::from_phases(Amount::new(0.2), Amount::new(0.0)).unwrap();
    for source in [&existing, &particle] {
        let step = source
            .kinetics()
            .advance(
                Temperature::new(300.0),
                Volume::new(1e-3),
                &initial,
                Time::new(1.0),
            )
            .unwrap();
        assert_eq!(step.state(), &initial);
        assert_eq!(step.transferred_to_solution().value(), 0.0);
    }
    let step = seeded
        .kinetics()
        .advance(
            Temperature::new(300.0),
            Volume::new(1e-3),
            &initial,
            Time::new(1.0),
        )
        .unwrap();
    near(
        step.state().dissolved().value(),
        0.1 + 0.1 * (-1.0_f64).exp(),
    );
    assert!(step.state().solid().value() > 0.0);
    near(
        step.state().dissolved().value() + step.state().solid().value(),
        0.2,
    );
}

#[test]
fn independent_support_intersection_rebases_the_same_physical_curve() {
    let equilibrium = admit(equilibrium_model(20_000.0)).unwrap();
    let mut original = particle_model();
    original.validity = domain("water", "A", 310.0, 360.0);
    let kinetic = bind(&equilibrium, original.clone()).unwrap();
    assert_eq!(kinetic.model_card(), &original);
    assert_eq!(
        kinetic.kinetics().law().temperature_domain(),
        [Temperature::new(310.0), Temperature::new(350.0)]
    );
    assert_eq!(
        kinetic.kinetics().law().reference_temperature(),
        Temperature::new(310.0)
    );
    assert_eq!(
        kinetic.equilibrium_source().law().reference_temperature(),
        Temperature::new(300.0)
    );
    for temperature in [310.0, 320.0, 350.0] {
        let t = Temperature::new(temperature);
        near(
            kinetic
                .kinetics()
                .law()
                .saturation_concentration(t)
                .unwrap()
                .value(),
            equilibrium
                .law()
                .saturation_concentration(t)
                .unwrap()
                .value(),
        );
    }
    for temperature in [300.0, 351.0] {
        assert!(
            kinetic
                .kinetics()
                .law()
                .saturation_concentration(Temperature::new(temperature))
                .is_err()
        );
    }
    original.validity = domain("water", "A", 360.0, 370.0);
    assert!(bind(&equilibrium, original).is_err());
}

#[test]
fn full_model_identity_exact_pair_and_declared_dimensions_are_required() {
    let original = equilibrium_model(0.0);
    let material = card(vec![original.clone()]);
    for (solvent, solute) in [("ethanol", "A"), ("water", "B"), ("A", "water")] {
        assert!(
            NeutralSolubilitySource::from_material(
                material.clone(),
                original.content_hash(),
                species(solvent),
                species(solute)
            )
            .is_err()
        );
    }
    assert!(
        NeutralSolubilitySource::from_material(
            material.clone(),
            original.canonical_parameters_hash().unwrap(),
            species("water"),
            species("A")
        )
        .is_err()
    );
    assert!(
        NeutralSolubilitySource::from_material(
            card(vec![original.clone(), original.clone()]),
            original.content_hash(),
            species("water"),
            species("A")
        )
        .is_err()
    );
    let mut wrong_units = original.clone();
    wrong_units
        .parameters
        .get_mut(REFERENCE_SOLUBILITY_PARAMETER)
        .unwrap()
        .dims = Dims::NONE;
    assert!(admit(wrong_units).is_err());
    let mut unsupported = original.clone();
    unsupported.parameters.insert(
        "mass_fraction".into(),
        LawParameter {
            value: 0.1,
            dims: Dims::NONE,
        },
    );
    assert!(admit(unsupported).is_err());
    let mut missing = original.clone();
    missing.parameters.remove(DISSOLUTION_ENTHALPY_PARAMETER);
    assert!(admit(missing).is_err());
    let mut missing_sources = original.clone();
    missing_sources.sources.clear();
    assert!(admit(missing_sources).is_err());
    let mut wrong_schema = original;
    wrong_schema.state_schema_version = 1;
    assert!(admit(wrong_schema).is_err());
}

#[test]
fn unsupported_environment_and_kinetic_reinterpretations_refuse() {
    let mut wrong_temperature = equilibrium_model(0.0);
    wrong_temperature.validity = wrong_temperature.validity.with_quantity(
        TEMPERATURE_AXIS,
        QuantitySpec::dimensional(Temperature::DIMS),
        280.0,
        350.0,
    );
    assert!(admit(wrong_temperature).is_err());
    let mut unowned_condition = equilibrium_model(0.0);
    unowned_condition.validity = unowned_condition.validity.with("pressure", 1e5, 2e5);
    assert!(admit(unowned_condition).is_err());
    let equilibrium = admit(equilibrium_model(0.0)).unwrap();
    assert!(bind(&equilibrium, fixed_model(0.5)).is_err());
    let mut foreign = particle_model();
    foreign.validity = domain("water", "B", 290.0, 330.0);
    assert!(bind(&equilibrium, foreign).is_err());
    let mut diffusivity = particle_model();
    diffusivity
        .parameters
        .get_mut(TRANSFER_CONDUCTANCE_PARAMETER)
        .unwrap()
        .dims = Dims([2, 0, -1, 0, 0, 0]);
    assert!(bind(&equilibrium, diffusivity).is_err());
    let mut no_geometry = particle_model();
    no_geometry
        .parameters
        .get_mut(REFERENCE_SOLID_PARAMETER)
        .unwrap()
        .value = 0.0;
    assert!(bind(&equilibrium, no_geometry).is_err());
    let mut no_initial_state = particle_model();
    no_initial_state.initial_state = InitialStatePolicy::ZeroInternalState;
    assert!(bind(&equilibrium, no_initial_state).is_err());
    let mut inactive = particle_model();
    inactive
        .parameters
        .get_mut(TRANSFER_CONDUCTANCE_PARAMETER)
        .unwrap()
        .value = 0.0;
    assert_eq!(
        bind(&equilibrium, inactive)
            .unwrap()
            .kinetics()
            .conductance()
            .value(),
        0.0
    );
}
