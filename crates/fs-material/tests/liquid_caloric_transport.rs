//! Manufactured nonideal liquids exercise the original caloric and material
//! owners together. The sources are analytic fixtures, not empirical liquids.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::solution::caloric::{BinaryLiquidCaloricError, BinaryLiquidCaloricLaw};
use fs_material::liquid::solution::{BinaryLiquidSolutionLaw, ResolvedBinaryLiquidEquilibrium};
use fs_material::liquid::{
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dimensionless, MolarMass, QuantitySpec, Temperature};
use fs_thermal::enthalpy::ThermalLaw;
use fs_thermochem::solution::caloric::BinaryCaloricSolveV1;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const PRESSURE: f64 = 100_000.0;
const MASSES: [f64; 2] = [0.018, 0.10];

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), p)
        .unwrap()
}
fn domain(p: f64) -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            250.0,
            650.0,
        )
        .with_quantity("p", quantity(QuantityKind::Pressure), p, p)
}
fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn capacity_curve(bias: f64) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: [250.0, 400.0, 650.0]
            .into_iter()
            .map(|t| (t, 1000.0 + t + bias))
            .collect(),
        dims: SPECIFIC_HEAT_CAPACITY_DIMS,
    }
}
fn pure_source(
    species: &str,
    pressure: f64,
    datum: f64,
    cp_bias: f64,
) -> (SpeciesId, IsobaricHeatCapacityLaw) {
    let key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: capacity_curve(cp_bias),
            validity: domain(pressure),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(&format!("pure liquid {species} Cp")),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: species.into(),
            phase: "liquid".into(),
            process: "caloric fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = point(250.0, pressure);
    let upper = point(650.0, pressure);
    let heat_capacity = IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &lower,
            upper: &upper,
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        },
        datum,
    )
    .unwrap();
    (SpeciesId::new(species).unwrap(), heat_capacity)
}
fn sources() -> [(SpeciesId, IsobaricHeatCapacityLaw); 2] {
    [
        pure_source("A", PRESSURE, 300.0, 0.0),
        pure_source("B", PRESSURE, 300.0, 0.0),
    ]
}
fn mixture_card(upper_x: f64, cp_bias: f64, marker: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let composition = quantity(QuantityKind::Composition(CompositionBasis::MoleFraction));
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, 1100.0, 800.0],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [0.0; 3],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            [0.6, 0.8, 0.2],
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            [1e-6, 6e-6, 3e-6],
        ),
    ] {
        let value = if name == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            capacity_curve(cp_bias)
        } else {
            PropertyValue::Curve {
                abscissa: "x_a".into(),
                abscissa_dims: Dimensionless::DIMS,
                knots: vec![(0.0, values[0]), (0.5, values[1]), (1.0, values[2])],
                dims: spec.dims(),
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value,
                validity: domain(PRESSURE).with_quantity("x_a", composition, 0.0, upper_x),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(marker),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "A+B".into(),
            phase: "liquid".into(),
            process: marker.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn transport(card: &MaterialCard, b: f64) -> BinaryLiquidSolutionLaw<'_> {
    let reference = point(500.0, PRESSURE)
        .with_quantity(
            "x_a",
            quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
            0.4,
        )
        .unwrap();
    let state = resolve_liquid_state_from_kinematic(
        card,
        &reference,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        b,
        Temperature::new(250.0),
        Temperature::new(650.0),
    )
    .unwrap();
    BinaryLiquidSolutionLaw::try_new(pair, MASSES.map(MolarMass::new), card, state, "T", "x_a")
        .unwrap()
}
fn law<'a, 'card>(liquid: &'a BinaryLiquidSolutionLaw<'card>) -> BinaryLiquidCaloricLaw<'a, 'card> {
    BinaryLiquidCaloricLaw::try_new(liquid, sources(), "p").unwrap()
}
fn controls() -> BinaryCaloricSolveV1 {
    BinaryCaloricSolveV1 {
        temperature_bracket: [Temperature::new(250.0), Temperature::new(650.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
        coexistence: BinaryCoexistenceOptionsV1::default(),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} +/- {tolerance:e}"
    );
}
fn sensible(t: f64) -> f64 {
    1000.0 * t + 0.5 * t * t
}
fn inverse_sensible(h: f64) -> f64 {
    2.0 * h / (1000.0 + (1_000_000.0 + 2.0 * h).sqrt())
}
fn verify_inventory_and_sources(
    state: &ResolvedBinaryLiquidEquilibrium,
    liquid: &BinaryLiquidSolutionLaw<'_>,
    expected: [f64; 2],
) {
    assert_eq!(state.total_component_moles(), expected);
    let expected_mass = expected[0] * MASSES[0] + expected[1] * MASSES[1];
    near(state.total_mass_kg(), expected_mass, 1e-15);
    near(
        state.phases().iter().map(|p| p.mass_kg()).sum(),
        expected_mass,
        1e-15,
    );
    near(
        state.phases().iter().map(|p| p.volume_m3()).sum(),
        state.total_volume_m3(),
        1e-18,
    );
    for (i, &amount) in expected.iter().enumerate() {
        near(
            state.phases().iter().map(|p| p.component_moles()[i]).sum(),
            amount,
            1e-14,
        );
    }
    for phase in state.phases() {
        let fluid = phase.transport_state().unwrap();
        let thermal = phase.thermal_transport_state().unwrap();
        assert_eq!(thermal.transport(), fluid);
        assert_eq!(thermal.temperature_k(), state.temperature().value());
        near(
            fluid.specific_heat_j_kg_k(),
            1000.0 + state.temperature().value(),
            1e-10,
        );
        near(
            phase.mass_kg() / phase.volume_m3(),
            fluid.density_kg_m3(),
            1e-9,
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("p"),
            Some(PRESSURE)
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("T"),
            Some(state.temperature().value())
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("x_a"),
            Some(phase.mole_fractions()[0])
        );
        for property in phase.liquid().material().properties() {
            let original = liquid
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
    }
}

#[test]
fn unequal_temperature_mixing_retains_variable_cp_excess_enthalpy_and_actual_transport() {
    let card = mixture_card(1.0, 0.0, "shared mixture source");
    for b in [0.0, 200.0, -200.0] {
        let liquid = transport(&card, b);
        let law = law(&liquid);
        let left = law
            .equilibrate(
                Temperature::new(300.0),
                [1.0, 0.0],
                controls().coexistence,
                || false,
            )
            .unwrap();
        let right = law
            .equilibrate(
                Temperature::new(500.0),
                [0.0, 3.0],
                controls().coexistence,
                || false,
            )
            .unwrap();
        let initial = [left.clone(), right.clone()];
        let result = law.mix(&left, &right, 0.0, controls(), || false).unwrap();
        let final_state = result.equilibrium();
        // Independent integral of Cp=1000+T, weighted by conserved MASS.
        // Four moles at x_A=1/4 contribute excess H=4 R B x_A x_B once.
        let excess = 4.0 * R * b * 0.25 * 0.75;
        let expected = inverse_sensible(
            (MASSES[0] * sensible(300.0) + 3.0 * MASSES[1] * sensible(500.0) - excess)
                / (MASSES[0] + 3.0 * MASSES[1]),
        );
        near(final_state.temperature().value(), expected, 2e-9);
        // Cp=1000+T integrates to 1000 ln(Tf/Ti)+(Tf-Ti) in
        // entropy. The initial parcels are pure; all ideal mixing entropy
        // belongs to the accepted four-mole final composition.
        let entropy = MASSES[0] * (1000.0 * (expected / 300.0).ln() + expected - 300.0)
            + 3.0 * MASSES[1] * (1000.0 * (expected / 500.0).ln() + expected - 500.0)
            - R * (0.25_f64.ln() + 3.0 * 0.75_f64.ln());
        near(result.entropy_change_j_per_k(), entropy, 5e-9);
        assert!(result.entropy_change_j_per_k() > 0.0);
        assert_eq!(
            result.entropy_change_j_per_k(),
            result.balance().entropy_change_j_per_k(),
        );
        assert!((final_state.temperature().value() - 450.0).abs() > 20.0);
        assert_eq!(final_state.phases().len(), 1);
        near(
            result.balance().energy_residual_j(),
            0.0,
            controls().energy_tolerance_j,
        );
        verify_inventory_and_sources(final_state, &liquid, [1.0, 3.0]);
        let fluid = final_state.phases()[0].transport_state().unwrap();
        near(fluid.density_kg_m3(), 1050.0, 1e-10);
        near(fluid.kinematic_viscosity_m2_s(), 3.5e-6, 1e-18);
        assert!(
            (final_state.total_volume_m3()
                - left.equilibrium().total_volume_m3()
                - right.equilibrium().total_volume_m3())
            .abs()
                > 1e-6
        );
        assert_eq!(
            result.final_state().caloric_state(),
            result.balance().final_state()
        );
        assert_eq!([left, right], initial);
        assert_eq!(law.sources()[0].source_line().knots.len(), 3);
    }
}

#[test]
fn adiabatic_demixing_heats_then_supplied_heat_remixes_the_same_finite_inventory() {
    let card = mixture_card(1.0, 0.0, "shared mixture source");
    let liquid = transport(&card, 900.0);
    let law = law(&liquid);
    let feed = [2.0, 1.0];
    let initial = law
        .homogeneous(Temperature::new(300.0), feed, || false)
        .unwrap();
    assert!(initial.caloric_state().equilibrium().is_none());
    assert!(initial.homogeneous_state().thermodynamic_factor() < 0.0);
    assert_eq!(initial.component_moles(), feed);
    near(initial.mass_kg(), 2.0 * MASSES[0] + MASSES[1], 1e-15);
    near(
        initial.volume_m3(),
        initial.mass_kg() / initial.material().density_kg_m3(),
        1e-18,
    );
    assert_eq!(
        initial.thermal_transport_state().unwrap().temperature_k(),
        300.0
    );
    let separated = law
        .relax_homogeneous(&initial, 0.0, controls(), || false)
        .unwrap();
    let cold = separated.equilibrium();
    assert_eq!(cold.phases().len(), 2);
    assert!(cold.temperature().value() > 300.0);
    // Independent 60-digit solution of the symmetric binodal equation and
    // the integrated pure-Cp energy balance, with the source constants above.
    near(cold.temperature().value(), 317.827_033_032_885_9, 2e-9);
    near(
        cold.phases()[0].mole_fractions()[0],
        0.088_734_605_892_491_94,
        2e-12,
    );
    near(
        separated.final_state().caloric_state().enthalpy_j(),
        initial.caloric_state().enthalpy_j(),
        controls().energy_tolerance_j,
    );
    near(
        separated.balance().enthalpy_change_j(),
        0.0,
        controls().energy_tolerance_j,
    );
    verify_inventory_and_sources(cold, &liquid, feed);
    let x = cold.phases()[0].mole_fractions();
    let mass = 2.0 * MASSES[0] + MASSES[1];
    let ideal_entropy = |x: f64| -R * (x * x.ln() + (1.0 - x) * (1.0 - x).ln());
    let expected_t = 317.827_033_032_885_9;
    let expected_x = 0.088_734_605_892_491_94;
    let demixing_entropy = mass * (1000.0 * (expected_t / 300.0_f64).ln() + expected_t - 300.0)
        + 3.0 * (ideal_entropy(expected_x) - ideal_entropy(2.0 / 3.0));
    near(separated.entropy_change_j_per_k(), demixing_entropy, 3e-9);
    assert!(separated.entropy_change_j_per_k() > 0.0);
    // Demixing loses ideal configurational entropy, but its released excess
    // enthalpy heats the same closed inventory enough to increase total S.
    assert!(ideal_entropy(expected_x) < ideal_entropy(2.0 / 3.0));
    // Direct source equation closes the heat needed to reach a homogeneous
    // 520 K liquid; changing excess enthalpy cannot be omitted or counted twice.
    let heat = mass * (sensible(520.0) - sensible(cold.temperature().value()))
        + 3.0 * R * 900.0 * ((2.0 / 3.0) * (1.0 / 3.0) - x[0] * x[1]);
    let warmed = law
        .add_heat(separated.final_state(), heat, controls(), || false)
        .unwrap();
    let warm = warmed.equilibrium();
    near(warm.temperature().value(), 520.0, 2e-9);
    assert_eq!(warm.phases().len(), 1);
    near(
        warmed.balance().enthalpy_change_j(),
        heat,
        controls().energy_tolerance_j,
    );
    verify_inventory_and_sources(warm, &liquid, feed);
    let warming_entropy = mass * (1000.0 * (520.0 / expected_t).ln() + 520.0 - expected_t)
        + 3.0 * (ideal_entropy(2.0 / 3.0) - ideal_entropy(expected_x));
    near(warmed.entropy_change_j_per_k(), warming_entropy, 3e-9);
    let cooled = law
        .add_heat(warmed.final_state(), -heat, controls(), || false)
        .unwrap();
    near(
        cooled.equilibrium().temperature().value(),
        cold.temperature().value(),
        3e-9,
    );
    assert_eq!(cooled.equilibrium().phases().len(), 2);
    verify_inventory_and_sources(cooled.equilibrium(), &liquid, feed);
    near(
        cooled.entropy_change_j_per_k(),
        -warmed.entropy_change_j_per_k(),
        5e-9,
    );
    assert!(cooled.entropy_change_j_per_k() < 0.0);
    // The new endpoint retains the existing equilibrium and thermal-storage APIs.
    let isothermal = liquid
        .mix_isothermal(cold, cold, controls().coexistence, || false)
        .unwrap();
    near(isothermal.heat_into_liquid_j(), 0.0, 1e-10);
    let storage = law.thermal_storage(feed, controls(), || false).unwrap();
    let h = storage
        .specific_enthalpy_j_kg(cold.temperature(), || false)
        .unwrap();
    near(
        storage.temperature_k(h).unwrap(),
        cold.temperature().value(),
        3e-9,
    );
    let warm_h = storage
        .specific_enthalpy_j_kg(warm.temperature(), || false)
        .unwrap();
    let storage_entropy = storage
        .specific_entropy_change_j_kg_k(h, warm_h)
        .unwrap()
        .unwrap();
    near(
        storage_entropy * mass,
        warmed.entropy_change_j_per_k(),
        1e-8,
    );
    assert!(initial.caloric_state().equilibrium().is_none());
    assert_eq!(initial.temperature().value(), 300.0);
}

#[test]
fn metastable_feed_preserves_initial_mixing_energy_before_heat_or_relaxation() {
    let card = mixture_card(1.0, 0.0, "shared mixture source");
    let liquid = transport(&card, 900.0);
    let law = law(&liquid);
    let feed = [1.0, 9.0];
    let initial = law
        .homogeneous(Temperature::new(300.0), feed, || false)
        .unwrap();
    // Locally stable but globally inside the miscibility gap: no implicit
    // phase split may replace the homogeneous energy at initialization.
    assert!(initial.homogeneous_state().thermodynamic_factor() > 0.0);
    assert!(initial.caloric_state().equilibrium().is_none());
    let isothermal = law
        .equilibrate(
            Temperature::new(300.0),
            feed,
            controls().coexistence,
            || false,
        )
        .unwrap();
    assert_eq!(isothermal.equilibrium().phases().len(), 2);
    assert!(initial.caloric_state().enthalpy_j() > isothermal.caloric_state().enthalpy_j());
    let relaxed = law
        .relax_homogeneous(&initial, 0.0, controls(), || false)
        .unwrap();
    assert_eq!(relaxed.equilibrium().phases().len(), 2);
    assert!(relaxed.equilibrium().temperature().value() > 300.0);
    near(
        relaxed.balance().enthalpy_change_j(),
        0.0,
        controls().energy_tolerance_j,
    );
    verify_inventory_and_sources(relaxed.equilibrium(), &liquid, feed);
    // Initial and final homogeneous compositions match, so their original
    // excess enthalpies cancel. All supplied heat is integrated sensible Cp.
    let mass = MASSES[0] + 9.0 * MASSES[1];
    let heat = mass * (sensible(520.0) - sensible(300.0));
    let heated = law
        .relax_homogeneous(&initial, heat, controls(), || false)
        .unwrap();
    near(heated.equilibrium().temperature().value(), 520.0, 2e-9);
    assert_eq!(heated.equilibrium().phases().len(), 1);
    near(
        heated.balance().enthalpy_change_j(),
        heat,
        controls().energy_tolerance_j,
    );
    verify_inventory_and_sources(heated.equilibrium(), &liquid, feed);
    assert_eq!(initial.temperature().value(), 300.0);
    assert!(initial.caloric_state().equilibrium().is_none());
}

#[test]
fn different_species_pressure_datums_or_material_sources_cannot_enter_the_same_heat_account() {
    let card = mixture_card(1.0, 0.0, "shared mixture source");
    let liquid = transport(&card, 900.0);
    let mut reversed = sources();
    reversed.swap(0, 1);
    assert!(BinaryLiquidCaloricLaw::try_new(&liquid, reversed, "p").is_err());
    for (p, datum) in [(200_000.0, 300.0), (PRESSURE, 310.0)] {
        let components = [
            pure_source("A", PRESSURE, 300.0, 0.0),
            pure_source("B", p, datum, 0.0),
        ];
        assert!(BinaryLiquidCaloricLaw::try_new(&liquid, components, "p").is_err());
    }
    let first = law(&liquid);
    let initial = first
        .homogeneous(Temperature::new(300.0), [1.0, 1.0], || false)
        .unwrap();
    let foreign_card = mixture_card(1.0, 0.0, "different original source");
    let foreign_liquid = transport(&foreign_card, 900.0);
    let foreign = law(&foreign_liquid);
    assert!(matches!(
        foreign.relax_homogeneous(&initial, 1.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Invalid(_))
    ));
    let different_liquid = transport(&card, 800.0);
    let changed_pair = law(&different_liquid);
    assert!(
        changed_pair
            .relax_homogeneous(&initial, 1.0, controls(), || false)
            .is_err()
    );
    let changed_calorics = BinaryLiquidCaloricLaw::try_new(
        &liquid,
        [
            pure_source("A", PRESSURE, 300.0, 0.0),
            pure_source("B", PRESSURE, 300.0, 25.0),
        ],
        "p",
    )
    .unwrap();
    assert!(matches!(
        changed_calorics.relax_homogeneous(&initial, 0.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Invalid(_))
    ));
    // Identically reconstructed immutable owners retain existing admission.
    let reconstructed_liquid = transport(&card, 900.0);
    let reconstructed = law(&reconstructed_liquid);
    assert_eq!(
        first
            .relax_homogeneous(&initial, 0.0, controls(), || false)
            .unwrap(),
        reconstructed
            .relax_homogeneous(&initial, 0.0, controls(), || false)
            .unwrap()
    );
    let original = initial.clone();
    let invalid_domain = BinaryCaloricSolveV1 {
        temperature_bracket: [Temperature::new(249.0), Temperature::new(650.0)],
        ..controls()
    };
    assert!(
        first
            .relax_homogeneous(&initial, 1.0, invalid_domain, || false)
            .is_err()
    );
    assert!(
        first
            .relax_homogeneous(&initial, 1e20, controls(), || false)
            .is_err()
    );
    assert_eq!(initial, original);
}

#[test]
fn unsupported_initial_or_second_phase_source_refuses_without_replacing_the_feed() {
    let limited_card = mixture_card(0.8, 0.0, "limited composition support");
    let liquid = transport(&limited_card, 900.0);
    let limited = law(&liquid);
    assert!(
        limited
            .homogeneous(Temperature::new(300.0), [9.0, 1.0], || false)
            .is_err()
    );
    assert!(
        limited
            .homogeneous(Temperature::new(249.0), [1.0, 1.0], || false)
            .is_err()
    );
    assert!(
        limited
            .homogeneous(Temperature::new(300.0), [-1.0, 2.0], || false)
            .is_err()
    );
    let initial = limited
        .homogeneous(Temperature::new(300.0), [1.0, 1.0], || false)
        .unwrap();
    let original = initial.clone();
    assert!(matches!(
        limited.relax_homogeneous(&initial, 0.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Material(_))
    ));
    assert_eq!(initial, original);
}

#[test]
fn cancellation_after_final_material_refresh_leaves_the_original_feed_intact() {
    let card = mixture_card(1.0, 0.0, "shared mixture source");
    let liquid = transport(&card, 900.0);
    let law = law(&liquid);
    assert!(matches!(
        law.homogeneous(Temperature::new(300.0), [2.0, 1.0], || true),
        Err(BinaryLiquidCaloricError::Cancelled)
    ));
    let initial = law
        .homogeneous(Temperature::new(300.0), [2.0, 1.0], || false)
        .unwrap();
    let original = initial.clone();
    let mut polls = 0;
    let accepted = law
        .relax_homogeneous(&initial, 0.0, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    let mut calls = 0;
    assert!(matches!(
        law.relax_homogeneous(&initial, 0.0, controls(), || {
            calls += 1;
            calls == polls
        }),
        Err(BinaryLiquidCaloricError::Cancelled)
    ));
    assert_eq!(calls, polls);
    assert_eq!(initial, original);
    assert_eq!(
        law.relax_homogeneous(&initial, 0.0, controls(), || false)
            .unwrap(),
        accepted
    );
}
