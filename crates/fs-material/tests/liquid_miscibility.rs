//! G0/G3 source-backed liquid miscibility. Synthetic constitutive data are
//! deliberately nonideal; these fixtures are not measurements of real liquids.
#[path = "liquid_miscibility/caloric.rs"]
mod caloric;
#[path = "liquid_miscibility/gradient.rs"]
mod gradient;
#[path = "liquid_miscibility/mixing.rs"]
mod mixing;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::solution::{BinaryLiquidError, BinaryLiquidSolutionLaw};
use fs_material::liquid::{
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, ResolvedLiquidState,
    resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dimensionless, MolarMass, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn composition(basis: CompositionBasis) -> QuantitySpec {
    quantity(QuantityKind::Composition(basis))
}
fn card(basis: CompositionBasis, upper_x: f64) -> MaterialCard {
    card_with_density(basis, upper_x, [1000.0, 1100.0, 800.0])
}
fn card_with_density(basis: CompositionBasis, upper_x: f64, density: [f64; 3]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            density,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [4000.0, 3500.0, 2000.0],
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
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value: PropertyValue::Curve {
                    abscissa: "x_a".into(),
                    abscissa_dims: Dimensionless::DIMS,
                    knots: vec![(0.0, values[0]), (0.5, values[1]), (1.0, values[2])],
                    dims: spec.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("x_a", composition(basis), 0.0, upper_x)
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        290.0,
                        600.0,
                    )
                    .with("pressure", 100_000.0, 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic binary solution transport".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A+B solution".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reference(card: &MaterialCard, basis: CompositionBasis) -> ResolvedLiquidState {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 500.0)
        .unwrap()
        .with_quantity("x_a", composition(basis), 0.4)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap();
    resolve_liquid_state_from_kinematic(card, &point, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}
fn model() -> BinaryRegularSolutionV1 {
    BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        900.0,
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap()
}
fn law(card: &MaterialCard) -> BinaryLiquidSolutionLaw<'_> {
    BinaryLiquidSolutionLaw::try_new(
        model(),
        [MolarMass::new(0.018), MolarMass::new(0.10)],
        card,
        reference(card, CompositionBasis::MoleFraction),
        "T",
        "x_a",
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5e-13 * b.abs().max(1e-20),
        "{a:e} != {b:e}"
    );
}

#[test]
fn cooling_splits_and_heating_remixes_the_same_conserved_components() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let feed = [2.0, 3.0];
    let cold = law
        .equilibrate(
            Temperature::new(300.0),
            feed,
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    let warm = law
        .equilibrate(
            Temperature::new(500.0),
            feed,
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    assert_eq!(cold.phases().len(), 2);
    assert_eq!(warm.phases().len(), 1);
    for result in [&cold, &warm] {
        assert_eq!(result.total_component_moles(), feed);
        near(result.total_mass_kg(), 0.336);
        near(result.phases().iter().map(|p| p.mass_kg()).sum(), 0.336);
        for (i, &amount) in feed.iter().enumerate() {
            near(
                result.phases().iter().map(|p| p.component_moles()[i]).sum(),
                amount,
            );
        }
        near(
            result.phases().iter().map(|p| p.volume_m3()).sum(),
            result.total_volume_m3(),
        );
        for phase in result.phases() {
            near(
                phase.mass_kg(),
                phase.component_moles()[0] * 0.018 + phase.component_moles()[1] * 0.10,
            );
            near(
                phase.volume_m3(),
                phase.mass_kg() / phase.liquid().density_kg_m3(),
            );
            near(
                phase.component_moles()[0] / phase.component_moles().iter().sum::<f64>(),
                phase.mole_fractions()[0],
            );
            near(phase.mass_fractions().iter().sum(), 1.0);
        }
    }
    // Component mass fractions and phase volume fractions must not inherit the
    // thermodynamic lever rule's molar basis.
    let low = &cold.phases()[0];
    assert!((low.mass_fractions()[0] - low.mole_fractions()[0]).abs() > 0.02);
    let mole_share = low.component_moles().iter().sum::<f64>() / 5.0;
    assert!((low.mass_kg() / cold.total_mass_kg() - mole_share).abs() > 0.1);
    assert!((low.volume_m3() / cold.total_volume_m3() - mole_share).abs() > 0.05);
}

#[test]
fn dilute_species_fractions_agree_with_their_conserved_inventories() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    for (feed, minority) in [([1.0, 1.5e-16], 1), ([1.5e-16, 1.0], 0)] {
        let result = law
            .equilibrate(
                Temperature::new(500.0),
                feed,
                BinaryCoexistenceOptionsV1::default(),
                || false,
            )
            .unwrap();
        assert_eq!(result.phases().len(), 1);
        assert_eq!(result.species(), law.thermodynamics().species());
        let phase = &result.phases()[0];
        assert_eq!(phase.species(), result.species());
        assert_eq!(phase.component_moles(), feed);
        near(
            phase.mole_fractions()[minority],
            feed[minority] / (feed[0] + feed[1]),
        );
        assert_eq!(
            phase.mole_fractions(),
            result.thermodynamic_equilibrium().feed().mole_fractions()
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("x_a"),
            Some(phase.mole_fractions()[0])
        );
    }
    // The material source has one x_A axis. It cannot represent a nonzero B
    // once x_A rounds to exactly one, even though the paired L1 law can.
    assert!(matches!(
        law.equilibrate(
            Temperature::new(500.0),
            [1.0, 1e-20],
            BinaryCoexistenceOptionsV1::default(),
            || false
        ),
        Err(BinaryLiquidError::Arithmetic(_))
    ));
}

#[test]
fn source_composition_reaches_flow_heat_and_concentration_without_property_averaging() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let result = law
        .equilibrate(
            Temperature::new(500.0),
            [1.0, 1.0],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    let phase = &result.phases()[0];
    let fluid = phase.transport_state().unwrap();
    // Density and viscosity at the actual middle knot, deliberately not the
    // average of the pure-liquid endpoints (900 kg/m3 and 2e-6 m2/s).
    assert_eq!(fluid.density_kg_m3(), 1100.0);
    near(fluid.kinematic_viscosity_m2_s(), 6e-6);
    near(fluid.dynamic_viscosity_pa_s(), 1100.0 * 6e-6);
    assert_eq!(fluid.specific_heat_j_kg_k(), 3500.0);
    assert_eq!(fluid.thermal_conductivity_w_m_k(), 0.8);
    near(phase.molar_concentration_mol_m3(), 1100.0 / 0.059);
    let thermal = phase.thermal_transport_state().unwrap();
    assert_eq!(thermal.transport(), fluid);
    assert_eq!(thermal.temperature_k(), 500.0);
    assert_eq!(thermal.volumetric_expansion_per_k(), None);
    for p in phase.liquid().material().properties() {
        let original = law
            .reference_state()
            .material()
            .property(p.requirement().name())
            .unwrap();
        assert_eq!(
            p.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
    assert_eq!(
        phase.liquid().material().state_coordinate("pressure"),
        Some(100_000.0)
    );
    assert_eq!(phase.liquid().material().state_coordinate("x_a"), Some(0.5));
}

#[test]
fn source_domain_loss_cannot_publish_one_of_two_required_phases() {
    let card = card(CompositionBasis::MoleFraction, 0.8);
    let law = law(&card);
    let before = law.reference_state().clone();
    assert!(matches!(
        law.equilibrate(
            Temperature::new(300.0),
            [2.0, 3.0],
            BinaryCoexistenceOptionsV1::default(),
            || false
        ),
        Err(BinaryLiquidError::Source(_))
    ));
    assert_eq!(law.reference_state(), &before);
    assert!(
        law.equilibrate(
            Temperature::new(500.0),
            [2.0, 3.0],
            BinaryCoexistenceOptionsV1::default(),
            || false
        )
        .is_ok()
    );
}

#[test]
fn pure_components_are_exact_and_invalid_inventories_refuse() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    for (feed, fractions) in [([2.0, 0.0], [1.0, 0.0]), ([0.0, 3.0], [0.0, 1.0])] {
        let result = law
            .equilibrate(
                Temperature::new(300.0),
                feed,
                BinaryCoexistenceOptionsV1::default(),
                || false,
            )
            .unwrap();
        assert_eq!(result.phases().len(), 1);
        assert_eq!(result.phases()[0].component_moles(), feed);
        assert_eq!(result.phases()[0].mole_fractions(), fractions);
        assert_eq!(result.phases()[0].mass_fractions(), fractions);
    }
    for feed in [
        [0.0, 0.0],
        [-1.0, 1.0],
        [f64::NAN, 1.0],
        [1.0, f64::INFINITY],
        [f64::MAX, f64::MAX],
        [1.0, 1e-100],
    ] {
        assert!(
            law.equilibrate(
                Temperature::new(300.0),
                feed,
                BinaryCoexistenceOptionsV1::default(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        law.equilibrate(
            Temperature::new(200.0),
            [1.0, 1.0],
            BinaryCoexistenceOptionsV1::default(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        law.equilibrate(
            Temperature::new(300.0),
            [1.0, 1.0],
            BinaryCoexistenceOptionsV1::default(),
            || true
        ),
        Err(BinaryLiquidError::Cancelled)
    ));
}

#[test]
fn composition_basis_and_original_card_are_enforced() {
    let wrong_basis = card(CompositionBasis::MassFraction, 1.0);
    assert!(matches!(
        BinaryLiquidSolutionLaw::try_new(
            model(),
            [MolarMass::new(0.018), MolarMass::new(0.10)],
            &wrong_basis,
            reference(&wrong_basis, CompositionBasis::MassFraction),
            "T",
            "x_a",
        ),
        Err(BinaryLiquidError::Invalid(_))
    ));
    let original = card(CompositionBasis::MoleFraction, 1.0);
    let different = card(CompositionBasis::MoleFraction, 0.8);
    assert!(
        BinaryLiquidSolutionLaw::try_new(
            model(),
            [MolarMass::new(0.018), MolarMass::new(0.10)],
            &different,
            reference(&original, CompositionBasis::MoleFraction),
            "T",
            "x_a",
        )
        .is_err()
    );
}
