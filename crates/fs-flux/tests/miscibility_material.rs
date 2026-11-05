//! G0/G1/G3/G4: conserved fields retain original homogeneous liquid sources.
//! All material curves below are synthetic numerical fixtures, not measurements.

use fs_evidence::ValidityDomain;
use fs_flux::miscibility::advection::{
    BinaryMiscibilityAdvectionControls, BinaryMiscibilityAdvectionFace,
};
use fs_flux::miscibility::material::BinaryMiscibilityMaterialError;
use fs_flux::miscibility::{
    BinaryMiscibilityCell, BinaryMiscibilityFace, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::solution::BinaryLiquidSolutionLaw;
use fs_material::liquid::solution::gradient::BinaryLiquidSquareGradientLaw;
use fs_material::liquid::{
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{AmountConcentration, Density, Dimensionless, MolarMass, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn composition() -> QuantitySpec {
    quantity(QuantityKind::Composition(CompositionBasis::MoleFraction))
}
fn card(upper_x: f64, middle_density: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, middle_density, 180.0],
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
                    .with_quantity("x_a", composition(), 0.0, upper_x)
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
                    source: "synthetic mixture field".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A+B".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn pair(interaction_k: f64) -> BinaryRegularSolutionV1 {
    BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        interaction_k,
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap()
}
fn liquid(card: &MaterialCard) -> BinaryLiquidSolutionLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("x_a", composition(), 0.4)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap();
    let reference = resolve_liquid_state_from_kinematic(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    BinaryLiquidSolutionLaw::try_new(
        pair(900.0),
        [MolarMass::new(0.018), MolarMass::new(0.10)],
        card,
        reference,
        "T",
        "x_a",
    )
    .unwrap()
}
fn gradient(
    temperature: f64,
    concentration: f64,
    kappa: f64,
    interaction_k: f64,
) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        pair(interaction_k),
        Temperature::new(temperature),
        AmountConcentration::new(concentration),
        kappa,
    )
    .unwrap()
}
fn model() -> BinarySquareGradientV1 {
    gradient(300.0, 10_000.0, 1e-6, 900.0)
}
fn network(model: BinarySquareGradientV1, diffusivity: f64) -> BinaryMiscibilityNetwork {
    BinaryMiscibilityNetwork::new(
        model,
        vec![BinaryMiscibilityCell { volume_m3: 1e-4 }; 2],
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: 1e-2,
            distance_m: 1e-2,
            diffusivity_m2_s: diffusivity,
        }],
        1 << 20,
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_stirring_refreshes_nonlinear_properties_and_keeps_spinodal_cells_homogeneous() {
    // rho=c*(x*M_A+(1-x)*M_B) supplies equal partial molar volumes while
    // unequal molar masses make the local densities change under transport.
    let card = card(1.0, 590.0);
    let liquid = liquid(&card);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let net = network(model(), 0.0);
    let initial = net
        .state_from_mole_fractions(vec![[1.0, 0.0], [0.0, 1.0]], || false)
        .unwrap();
    let field = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 5e-5,
                reverse_m3_s: 5e-5,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let candidate = field
        .transport_candidate(
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    let final_field = candidate.field();
    assert!(std::ptr::eq(final_field.source(), &source));
    assert_eq!(final_field.state().elapsed_s(), 1.0);
    assert_eq!(final_field.maximum_resolved_cells(), 2);
    near(
        final_field
            .cells()
            .iter()
            .map(|cell| cell.phase().mass_kg())
            .sum(),
        0.118,
    );
    for cell in final_field.cells() {
        assert_eq!(cell.phase().component_moles(), [0.5, 0.5]);
        assert_eq!(cell.phase().mole_fractions(), [0.5, 0.5]);
        assert!(cell.thermodynamics().thermodynamic_factor() < 0.0);
        near(cell.phase().volume_m3(), 1e-4);
        let actual = cell.phase().transport_state().unwrap();
        near(actual.density_kg_m3(), 590.0);
        near(actual.kinematic_viscosity_m2_s(), 6e-6);
        near(actual.dynamic_viscosity_pa_s(), 590.0 * 6e-6);
        near(actual.specific_heat_j_kg_k(), 3500.0);
        near(actual.thermal_conductivity_w_m_k(), 0.8);
        assert_eq!(
            cell.phase()
                .thermal_transport_state()
                .unwrap()
                .temperature_k(),
            300.0
        );
        for property in cell.phase().liquid().material().properties() {
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
    // Bulk equilibrium would split this same inventory; the spatial state must not.
    assert_eq!(
        liquid
            .equilibrate(
                Temperature::new(300.0),
                [0.5, 0.5],
                BinaryCoexistenceOptionsV1::default(),
                || false
            )
            .unwrap()
            .phases()
            .len(),
        2
    );
    assert_eq!(field.state(), &initial);
    assert_eq!(field.cells()[0].phase().component_moles(), [1.0, 0.0]);
    let (accepted, step) = candidate.into_parts();
    assert!(std::ptr::eq(accepted.source(), field.source()));
    assert_eq!(step.duration_s, 1.0);
}

#[test]
fn g0_full_gradient_owner_budget_and_source_volume_are_enforced() {
    let card = card(1.0, 590.0);
    let liquid = liquid(&card);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    for altered in [
        gradient(301.0, 10_000.0, 1e-6, 900.0),
        gradient(300.0, 20_000.0, 1e-6, 900.0),
        gradient(300.0, 10_000.0, 2e-6, 900.0),
        gradient(300.0, 10_000.0, 1e-6, 950.0),
    ] {
        let net = network(altered, 0.0);
        let state = net
            .state_from_mole_fractions(vec![[0.4, 0.6]; 2], || false)
            .unwrap();
        assert!(matches!(
            net.resolve_liquid_field(&state, &source, 2, || false),
            Err(BinaryMiscibilityMaterialError::Invalid(_))
        ));
    }
    let net = network(model(), 0.0);
    let state = net
        .state_from_mole_fractions(vec![[0.5, 0.5]; 2], || false)
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        net.resolve_liquid_field(&state, &source, 1, || {
            calls += 1;
            false
        }),
        Err(BinaryMiscibilityMaterialError::CellBudget {
            required: 2,
            maximum: 1
        })
    ));
    assert_eq!(calls, 1, "the budget must refuse before source queries");
    let other = network(model(), 1e-5);
    assert!(matches!(
        other.resolve_liquid_field(&state, &source, 2, || false),
        Err(BinaryMiscibilityMaterialError::Solver(_))
    ));
    let incompatible = self::card(1.0, 700.0);
    let bad_liquid = self::liquid(&incompatible);
    let bad_source = BinaryLiquidSquareGradientLaw::try_new(&bad_liquid, model(), 1e-12).unwrap();
    assert!(matches!(
        net.resolve_liquid_field(&state, &bad_source, 2, || false),
        Err(BinaryMiscibilityMaterialError::Source { cell: 0, .. })
    ));
}

#[test]
fn g4_final_diffusion_source_domain_failure_discards_the_complete_candidate() {
    let card = card(0.6, 590.0);
    let liquid = liquid(&card);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let net = network(model(), 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
        .unwrap();
    let field = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let before_cells = field.cells().to_vec();
    let flow = net
        .prepare_advection(vec![BinaryMiscibilityAdvectionFace::default()], 0.0, || {
            false
        })
        .unwrap();
    let numerical = net
        .transport_candidate(
            &initial,
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    assert!(
        numerical.state().mole_fractions()[1][0] > 0.6,
        "spinodal diffusion must leave the source domain"
    );
    assert!(matches!(
        field.transport_candidate(
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false
        ),
        Err(BinaryMiscibilityMaterialError::Source { cell: 1, .. })
    ));
    assert_eq!(field.state(), &initial);
    assert_eq!(field.cells(), before_cells);
}

#[test]
fn g4_cancellation_during_binding_and_final_source_refresh_preserves_the_original() {
    let card = card(1.0, 590.0);
    let liquid = liquid(&card);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let net = network(model(), 0.0);
    let state = net
        .state_from_mole_fractions(vec![[0.2, 0.8], [0.8, 0.2]], || false)
        .unwrap();
    let mut polls = 0;
    let field = net
        .resolve_liquid_field(&state, &source, 2, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, 3, polls] {
        let mut calls = 0;
        assert!(
            net.resolve_liquid_field(&state, &source, 2, || {
                calls += 1;
                calls == stop
            })
            .is_err()
        );
    }
    let before_cells = field.cells().to_vec();
    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 1e-5,
                reverse_m3_s: 1e-5,
            }],
            0.0,
            || false,
        )
        .unwrap();
    polls = 0;
    field
        .transport_candidate(
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls - 3, polls] {
        let mut calls = 0;
        assert!(
            field
                .transport_candidate(
                    &flow,
                    0.1,
                    BinaryMiscibilityAdvectionControls::default(),
                    BinaryMiscibilitySolve::default(),
                    || {
                        calls += 1;
                        calls == stop
                    }
                )
                .is_err()
        );
        assert_eq!(field.state(), &state);
        assert_eq!(field.cells(), before_cells);
    }
}
