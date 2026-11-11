//! Source-resolved ternary cross diffusion and atomic original-source refusals.
//! All molecular coefficients and material surfaces are synthetic test data.

use super::*;
use fs_flux::miscibility::multicomponent::adaptive::MulticomponentMiscibilityTimeControl;
use fs_flux::miscibility::multicomponent::material::MulticomponentMiscibilityMaterialError as Error;
use fs_flux::miscibility::multicomponent::{
    MulticomponentMiscibilityCell, MulticomponentMiscibilityFace, MulticomponentMiscibilityField,
    MulticomponentMiscibilityNetwork, MulticomponentMiscibilitySolve,
};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidComponent, MulticomponentLiquidSolutionLaw,
};
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;

const MASSES: [f64; 3] = [0.002, 0.028, 0.044];

#[path = "multicomponent/nonisothermal.rs"]
mod nonisothermal;
#[path = "multicomponent/phase_field.rs"]
mod phase_field;
#[path = "multicomponent/transport.rs"]
mod transport;

fn species(index: usize) -> SpeciesId {
    SpeciesId::new(["A", "B", "C"][index]).unwrap()
}

fn model(order: [usize; 3]) -> MulticomponentRegularSolutionV1 {
    let pairs = [0.0, 0.2, -0.3, 0.2, 0.0, 0.1, -0.3, 0.1, 0.0];
    MulticomponentRegularSolutionV1::new(
        order.map(species).into(),
        order
            .into_iter()
            .flat_map(|i| order.map(|j| pairs[3 * i + j]))
            .collect(),
        vec![0.0; 9],
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap()
}

fn source_card(maximum_x_a: f64, density_perturbation: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    // At c=1000 mol/m³, rho=44-42*x_A-16*x_B gives the original molecular
    // masses and fixed volume exactly. The perturbation is zero on the line
    // through BOTH initial cells, but cross diffusion leaves that line.
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            [
                44.0 - 6.3 * density_perturbation,
                2.0 + 13.7 * density_perturbation,
                28.0 - 5.3 * density_perturbation,
            ],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [1500.0, 2500.0, 4000.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            [0.1, 0.4, 0.8],
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            [1e-6, 2e-6, 4e-6],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value: PropertyValue::Surface2D {
                    axes: ["x_a".into(), "x_b".into()],
                    axis_dims: [Dimensionless::DIMS; 2],
                    samples: [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
                        .into_iter()
                        .zip(values)
                        .collect(),
                    triangles: vec![[0, 1, 2]],
                    dims: spec.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("x_a", composition(), 0.0, maximum_x_a)
                    .with_quantity("x_b", composition(), 0.0, 1.0)
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        290.0,
                        600.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        100_000.0,
                        100_000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic ternary liquid; not measured mixture data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A+B+C".into(),
            phase: "liquid".into(),
            process: "ternary numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source(
    card: &MaterialCard,
    thermodynamics: MulticomponentRegularSolutionV1,
    masses: [f64; 3],
) -> MulticomponentLiquidSolutionLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("x_a", composition(), 0.3)
        .unwrap()
        .with_quantity("x_b", composition(), 0.3)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    let reference = resolve_liquid_state_from_kinematic(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MulticomponentLiquidSolutionLaw::try_new(
        thermodynamics,
        [2, 0, 1]
            .map(|i| MulticomponentLiquidComponent {
                species: species(i),
                molar_mass: MolarMass::new(masses[i]),
                mole_fraction_axis: [Some("x_a"), Some("x_b"), None][i].map(str::to_owned),
            })
            .into(),
        card,
        reference,
        "T",
    )
    .unwrap()
}

fn numerical() -> MulticomponentMiscibilityField {
    let transport = MaxwellStefanV1::new(
        [1, 2, 0]
            .map(|i| (species(i), MolarMass::new(MASSES[i])))
            .into(),
    )
    .unwrap();
    MulticomponentMiscibilityNetwork::new(
        model([2, 1, 0]),
        transport,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        vec![MulticomponentMiscibilityCell { volume_m3: 0.001 }; 2],
        vec![MulticomponentMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_matrix: 0,
        }],
        vec![vec![0.0, 1e-6, 1e-4, 1e-6, 0.0, 1e-4, 1e-4, 1e-4, 0.0]],
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
    .bind(vec![vec![0.3, 0.3, 0.4], vec![0.31, 0.1, 0.59]], || false)
    .unwrap()
}

fn time() -> MulticomponentMiscibilityTimeControl {
    MulticomponentMiscibilityTimeControl {
        initial_step_s: 0.1,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.1,
        absolute_mole_fraction_tolerance: 1e-5,
        relative_mole_fraction_tolerance: 0.0,
        maximum_attempts: 1000,
    }
}

#[test]
fn ternary_cross_diffusion_refreshes_original_properties_and_replays_after_cancellation() {
    let card = source_card(1.0, 0.0);
    // Independent construction/permutation of the SAME complete law is valid.
    let source = source(&card, model([1, 0, 2]), MASSES);
    let numerical = numerical();
    let field = numerical
        .resolve_liquid_field(&source, 2, 1e-12, || false)
        .unwrap();
    let before = field.field().clone();
    let before_cells = field.cells().to_vec();
    let controls = MulticomponentMiscibilitySolve::default();
    let mut polls = 0;
    let candidate = field
        .step_candidate(0.1, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    assert!(std::ptr::eq(end.source(), &source));
    assert_eq!(end.field().elapsed_s(), 0.1);
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert_eq!(end.maximum_relative_volume_mismatch(), 1e-12);
    // Geometry, all authored pair matrices and original numerical owners persist.
    assert_eq!(end.network(), numerical.network());
    // A moves uphill due to the full cross-component friction/force system.
    assert!(end.field().cell_component_moles(0).unwrap()[0] < 0.3);
    assert!(end.field().cell_component_moles(1).unwrap()[0] > 0.31);
    for (cell, resolved) in end.cells().iter().enumerate() {
        let amounts = end.field().cell_component_moles(cell).unwrap();
        let fractions = end.field().cell_mole_fractions(cell).unwrap();
        assert_eq!(resolved.component_moles(), amounts);
        assert_eq!(resolved.species(), source.species());
        assert_eq!(resolved.thermodynamics().model(), source.thermodynamics());
        assert_eq!(
            resolved.thermodynamics().temperature(),
            Temperature::new(300.0)
        );
        near(resolved.volume_m3(), 0.001);
        near(resolved.molar_concentration_mol_m3(), 1000.0);
        let transport = resolved.transport_state().unwrap();
        near(
            transport.density_kg_m3(),
            44.0 - 42.0 * fractions[0] - 16.0 * fractions[1],
        );
        near(
            transport.specific_heat_j_kg_k(),
            1500.0 + 1000.0 * fractions[0] + 2500.0 * fractions[1],
        );
        let material = resolved.liquid().material();
        assert_eq!(material.card_identity(), card.content_hash());
        assert_eq!(material.state_coordinate("T"), Some(300.0));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        near(material.state_coordinate("x_a").unwrap(), fractions[0]);
        near(material.state_coordinate("x_b").unwrap(), fractions[1]);
    }
    for species in 0..3 {
        near(
            end.field().cell_component_moles(0).unwrap()[species]
                + end.field().cell_component_moles(1).unwrap()[species],
            before.cell_component_moles(0).unwrap()[species]
                + before.cell_component_moles(1).unwrap()[species],
        );
    }
    let advanced = end
        .advance_candidate(0.1, time(), controls, || false)
        .unwrap();
    assert!(advanced.advance().accepted_intervals > 0);
    assert!(advanced.advance().maximum_accepted_error_ratio <= 1.0);
    let (accepted, report) = advanced.into_parts();
    assert_eq!(accepted.field().elapsed_s(), 0.2);
    assert_eq!(report.duration_s, 0.1);
    assert!(std::ptr::eq(accepted.source(), &source));
    for (cell, resolved) in accepted.cells().iter().enumerate() {
        assert_eq!(
            resolved.component_moles(),
            accepted.field().cell_component_moles(cell).unwrap()
        );
        near(resolved.volume_m3(), 0.001);
    }
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(
            field
                .step_candidate(0.1, controls, || {
                    calls += 1;
                    calls == stop
                })
                .is_err()
        );
        assert_eq!(field.field(), &before);
        assert_eq!(field.cells(), before_cells.as_slice());
    }
    let replay = field.step_candidate(0.1, controls, || false).unwrap();
    assert_eq!(replay.field().field(), candidate.field().field());
    assert_eq!(replay.field().cells(), candidate.field().cells());
    assert_eq!(replay.step(), candidate.step());
}

#[test]
fn ternary_material_binding_and_late_source_refusals_preserve_original_state() {
    let card = source_card(1.0, 0.0);
    let source = source(&card, model([0, 1, 2]), MASSES);
    let numerical = numerical();
    assert!(matches!(
        numerical.resolve_liquid_field(&source, 1, 1e-12, || false),
        Err(Error::CellBudget {
            required: 2,
            maximum: 1
        })
    ));
    let other_model = MulticomponentRegularSolutionV1::new(
        [0, 1, 2].map(species).into(),
        vec![0.0; 9],
        vec![0.0; 9],
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let other_law = self::source(&card, other_model, MASSES);
    assert!(matches!(
        numerical.resolve_liquid_field(&other_law, 2, 1e-12, || false),
        Err(Error::Invalid(_))
    ));
    let wrong_masses = self::source(&card, model([0, 1, 2]), [0.002, 0.028, 0.045]);
    assert!(matches!(
        numerical.resolve_liquid_field(&wrong_masses, 2, 1e-12, || false),
        Err(Error::Invalid(_))
    ));
    let controls = MulticomponentMiscibilitySolve::default();
    for (maximum_x_a, density_perturbation) in [(0.31, 0.0), (1.0, 1.0)] {
        let card = source_card(maximum_x_a, density_perturbation);
        let source = self::source(&card, model([2, 0, 1]), MASSES);
        let field = numerical
            .resolve_liquid_field(&source, 2, 1e-12, || false)
            .unwrap();
        let original_cells = field.cells().to_vec();
        let original = field.field().clone();
        // The complete numerical step succeeds. Only original property support
        // or density compatibility refuses the material endpoint afterwards.
        field
            .field()
            .step_candidate(0.1, controls, || false)
            .unwrap();
        for error in [
            field.step_candidate(0.1, controls, || false).unwrap_err(),
            field
                .advance_candidate(0.1, time(), controls, || false)
                .unwrap_err(),
        ] {
            if density_perturbation == 0.0 {
                assert!(matches!(error, Error::Source { cell: 1, .. }));
            } else {
                assert!(matches!(error, Error::VolumeMismatch { .. }));
            }
        }
        assert_eq!(field.field(), &original);
        assert_eq!(field.cells(), original_cells.as_slice());
    }
}
