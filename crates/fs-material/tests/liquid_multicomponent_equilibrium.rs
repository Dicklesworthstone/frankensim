//! Finite phase inventories resolve their own original mixture-property sources.
//! Numerical fixtures are authored regular solutions, not measured solvent data.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidComponent, MulticomponentLiquidError, MulticomponentLiquidSolutionLaw,
};
use fs_material::liquid::{
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dimensionless, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::solution::multicomponent::equilibrium::{
    MulticomponentEquilibriumCandidateV1, MulticomponentEquilibriumControlV1,
    MulticomponentEquilibriumSeedV1,
};

const MASSES: [f64; 3] = [0.018, 0.046, 0.100];
const DENSITIES: [f64; 4] = [800.0, 1000.0, 1400.0, 1200.0];
const HEAT_CAPACITIES: [f64; 4] = [1000.0, 2000.0, 4000.0, 3000.0];
const CONDUCTIVITIES: [f64; 4] = [0.1, 0.2, 0.7, 0.9];
const VISCOSITIES: [f64; 4] = [1e-6, 2e-6, 3e-6, 9e-6];

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn mole_fraction() -> QuantitySpec {
    quantity(QuantityKind::Composition(CompositionBasis::MoleFraction))
}
fn species(index: usize) -> SpeciesId {
    SpeciesId::new(["A", "B", "C"][index]).unwrap()
}
fn model(chi: f64) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        (0..3).map(species).collect(),
        vec![0.0, chi, chi, chi, 0.0, chi, chi, chi, 0.0],
        vec![0.0; 9],
        Temperature::new(250.0),
        Temperature::new(650.0),
    )
    .unwrap()
}
fn card(upper_a: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            DENSITIES,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            HEAT_CAPACITIES,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            CONDUCTIVITIES,
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            VISCOSITIES,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value: PropertyValue::Surface2D {
                    axes: ["x_a".into(), "x_b".into()],
                    axis_dims: [Dimensionless::DIMS; 2],
                    samples: [
                        [0.0, 0.0],
                        [upper_a, 0.0],
                        [0.0, 1.0],
                        [upper_a / 3.0, 1.0 / 3.0],
                    ]
                    .into_iter()
                    .zip(values)
                    .collect(),
                    triangles: vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
                    dims: spec.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("x_a", mole_fraction(), 0.0, upper_a)
                    .with_quantity("x_b", mole_fraction(), 0.0, 1.0)
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        600.0,
                    )
                    .with_quantity(
                        "p",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        100_000.0,
                        100_000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored ternary property surface; not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "A+B+C".into(),
            phase: "liquid".into(),
            process: "authored phase-source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn law(
    card: &MaterialCard,
    thermodynamics: MulticomponentRegularSolutionV1,
) -> MulticomponentLiquidSolutionLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 350.0)
        .unwrap()
        .with_quantity("x_a", mole_fraction(), 0.2)
        .unwrap()
        .with_quantity("x_b", mole_fraction(), 0.3)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
        .unwrap();
    let reference = resolve_liquid_state_from_kinematic(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let components = [2, 0, 1]
        .into_iter()
        .map(|i| MulticomponentLiquidComponent {
            species: species(i),
            molar_mass: MolarMass::new(MASSES[i]),
            mole_fraction_axis: [Some("x_a"), Some("x_b"), None][i].map(str::to_owned),
        })
        .collect();
    MulticomponentLiquidSolutionLaw::try_new(thermodynamics, components, card, reference, "T")
        .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-12 * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn source_value(fractions: &[f64], values: [f64; 4]) -> f64 {
    // Independent barycentric evaluation of the three authored triangles.
    let minimum = fractions.iter().copied().fold(f64::INFINITY, f64::min);
    (fractions[2] - minimum) * values[0]
        + (fractions[0] - minimum) * values[1]
        + (fractions[1] - minimum) * values[2]
        + 3.0 * minimum * values[3]
}

fn candidate(
    model: &MulticomponentRegularSolutionV1,
    temperature: f64,
    amounts: [f64; 3],
    split: bool,
) -> MulticomponentEquilibriumCandidateV1 {
    let total: f64 = amounts.iter().sum();
    let feed = amounts.map(|n| n / total);
    let seed = if split {
        // Independent symmetric ternary common-potential solution for chi=4:
        // rich r=1-2q, lean q, ln(r/q)=4(r-q). Start near the actual selected
        // three-phase candidate so these tests isolate material integration.
        let (mut lo, mut hi) = (1e-12, 0.1);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            let residual = ((1.0_f64 - 2.0 * q) / q).ln() - 4.0 * (1.0 - 3.0 * q);
            if residual > 0.0 {
                lo = q;
            } else {
                hi = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let r = 1.0 - 2.0 * q;
        MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![vec![r, q, q], vec![q, r, q], vec![q, q, r]],
            phase_weights: feed.map(|x| (x - q) / (r - q)).into(),
        }
    } else {
        MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![feed.into()],
            phase_weights: vec![1.0],
        }
    };
    model
        .equilibrium_candidate(
            Temperature::new(temperature),
            &amounts,
            &seed,
            MulticomponentEquilibriumControlV1::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn g1_three_phase_inventory_resolves_each_actual_property_surface() {
    let card = card(1.0);
    let law = law(&card, model(4.0));
    let amounts = [2.0, 3.0, 5.0];
    let selected = candidate(law.thermodynamics(), 400.0, amounts, true);
    assert_eq!(selected.phases().len(), 3);
    let result = law.resolve_equilibrium(selected.clone(), || false).unwrap();
    assert_eq!(result.equilibrium_candidate(), &selected);
    assert_eq!(result.species(), law.species());
    assert_eq!(result.components(), law.components());
    assert_eq!(result.total_component_moles(), amounts);
    assert_eq!(result.temperature(), Temperature::new(400.0));
    near(result.total_moles(), 10.0);

    let mut volume = 0.0;
    let mut mass = 0.0;
    for (phase, accepted) in result.phases().iter().zip(selected.phases()) {
        assert_eq!(phase.component_moles(), accepted.component_moles());
        let x = phase.mole_fractions();
        for (&resolved, &selected) in x.iter().zip(accepted.thermodynamics().mole_fractions()) {
            near(resolved, selected);
        }
        let density = source_value(x, DENSITIES);
        let viscosity = source_value(x, VISCOSITIES);
        near(phase.liquid().density_kg_m3(), density);
        near(
            phase.liquid().specific_heat_j_kg_k(),
            source_value(x, HEAT_CAPACITIES),
        );
        near(
            phase.liquid().thermal_conductivity_w_m_k(),
            source_value(x, CONDUCTIVITIES),
        );
        near(phase.liquid().kinematic_viscosity_m2_s(), viscosity);
        near(phase.liquid().dynamic_viscosity_pa_s(), density * viscosity);
        let phase_mass: f64 = phase
            .component_moles()
            .iter()
            .zip(MASSES)
            .map(|(n, m)| n * m)
            .sum();
        near(phase.mass_kg(), phase_mass);
        near(phase.volume_m3(), phase_mass / density);
        mass += phase_mass;
        volume += phase_mass / density;
        assert_eq!(
            phase.liquid().material().state_coordinate("x_a"),
            Some(x[0])
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("x_b"),
            Some(x[1])
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("p"),
            Some(100_000.0)
        );
        assert_eq!(
            phase.liquid().material().card_identity(),
            card.content_hash()
        );
        for property in phase.liquid().material().properties() {
            let original = law
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
        let thermal = phase.thermal_transport_state().unwrap();
        assert_eq!(thermal.transport(), phase.transport_state().unwrap());
        assert_eq!(thermal.temperature_k(), 400.0);
    }
    near(
        result.total_mass_kg(),
        amounts.into_iter().zip(MASSES).map(|(n, m)| n * m).sum(),
    );
    near(result.total_mass_kg(), mass);
    near(result.total_volume_m3(), volume);
    for (i, amount) in amounts.into_iter().enumerate() {
        near(
            result
                .phases()
                .iter()
                .map(|phase| phase.component_moles()[i])
                .sum(),
            amount,
        );
    }
    let homogeneous = law
        .homogeneous(Temperature::new(400.0), &amounts, || false)
        .unwrap();
    assert!(
        (result.total_volume_m3() - homogeneous.volume_m3()).abs() > 0.01 * homogeneous.volume_m3()
    );
    assert_ne!(
        result.phases()[0].transport_state().unwrap(),
        result.phases()[1].transport_state().unwrap()
    );
}

#[test]
fn g1_one_phase_keeps_exact_inventory_absences_and_explicit_traces() {
    let card = card(1.0);
    let law = law(&card, model(0.0));
    for amounts in [[2.0, 3.0, 5.0], [0.0, 3.0, 7.0], [1e-20, 0.0, 10.0]] {
        let selected = candidate(law.thermodynamics(), 400.0, amounts, false);
        let result = law.resolve_equilibrium(selected, || false).unwrap();
        let direct = law
            .homogeneous(Temperature::new(400.0), &amounts, || false)
            .unwrap();
        assert_eq!(result.phases(), [direct.clone()]);
        assert_eq!(result.total_component_moles(), amounts);
        assert_eq!(result.total_mass_kg(), direct.mass_kg());
        assert_eq!(result.total_volume_m3(), direct.volume_m3());
    }
}

#[test]
fn g3_scaling_finite_feed_preserves_phase_compositions_and_source_properties() {
    let card = card(1.0);
    let law = law(&card, model(4.0));
    let base = law
        .resolve_equilibrium(
            candidate(law.thermodynamics(), 400.0, [2.0, 3.0, 5.0], true),
            || false,
        )
        .unwrap();
    for scale in [1e-90, 1e90] {
        let amounts = [2.0 * scale, 3.0 * scale, 5.0 * scale];
        let result = law
            .resolve_equilibrium(
                candidate(law.thermodynamics(), 400.0, amounts, true),
                || false,
            )
            .unwrap();
        assert_eq!(result.total_component_moles(), amounts);
        assert_eq!(result.phases().len(), base.phases().len());
        near(result.total_mass_kg(), base.total_mass_kg() * scale);
        near(result.total_volume_m3(), base.total_volume_m3() * scale);
        for (phase, original) in result.phases().iter().zip(base.phases()) {
            for (&x, &y) in phase.mole_fractions().iter().zip(original.mole_fractions()) {
                near(x, y);
            }
            near(
                phase.liquid().density_kg_m3(),
                original.liquid().density_kg_m3(),
            );
            near(
                phase.liquid().kinematic_viscosity_m2_s(),
                original.liquid().kinematic_viscosity_m2_s(),
            );
            near(phase.mass_kg(), original.mass_kg() * scale);
            near(phase.volume_m3(), original.volume_m3() * scale);
            near(
                phase.molar_concentration_mol_m3(),
                original.molar_concentration_mol_m3(),
            );
        }
    }
}

#[test]
fn g0_foreign_owners_and_any_unsupported_phase_refuse_atomically() {
    let full_card = card(1.0);
    let full_law = law(&full_card, model(4.0));
    let selected = candidate(full_law.thermodynamics(), 400.0, [2.0, 3.0, 5.0], true);
    let foreign = law(&full_card, model(3.5));
    assert!(matches!(
        foreign.resolve_equilibrium(selected.clone(), || false),
        Err(MulticomponentLiquidError::Invalid(_))
    ));

    let narrow_card = card(0.8);
    let narrow = law(&narrow_card, model(4.0));
    let reference = narrow.reference_state().clone();
    let supported = selected
        .phases()
        .iter()
        .filter(|phase| {
            narrow
                .homogeneous(Temperature::new(400.0), phase.component_moles(), || false)
                .is_ok()
        })
        .count();
    assert_eq!(
        supported, 2,
        "only the A-rich phase is outside the source table"
    );
    assert!(matches!(
        narrow.resolve_equilibrium(selected.clone(), || false),
        Err(MulticomponentLiquidError::Source(_))
    ));
    assert_eq!(narrow.reference_state(), &reference);
    assert_eq!(selected.component_moles(), [2.0, 3.0, 5.0]);

    let hot = candidate(full_law.thermodynamics(), 625.0, [2.0, 3.0, 5.0], true);
    assert!(matches!(
        full_law.resolve_equilibrium(hot, || false),
        Err(MulticomponentLiquidError::Source(_))
    ));
    let ideal = law(&full_card, model(0.0));
    let trace = candidate(ideal.thermodynamics(), 400.0, [0.5, 0.5, 1e-20], false);
    assert!(matches!(
        ideal.resolve_equilibrium(trace, || false),
        Err(MulticomponentLiquidError::Arithmetic(_))
    ));
}

#[test]
fn g4_every_phase_and_final_accounting_cancellation_publish_no_partial_result() {
    let card = card(1.0);
    let law = law(&card, model(4.0));
    let selected = candidate(law.thermodynamics(), 400.0, [2.0, 3.0, 5.0], true);
    let reference = law.reference_state().clone();
    let mut polls = 0;
    let complete = law
        .resolve_equilibrium(selected.clone(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut current = 0;
        assert!(
            law.resolve_equilibrium(selected.clone(), || {
                current += 1;
                current == stop
            })
            .is_err(),
            "poll {stop}"
        );
        assert_eq!(law.reference_state(), &reference);
    }
    assert_eq!(
        law.resolve_equilibrium(selected, || false).unwrap(),
        complete
    );
}

#[test]
fn g1_all_explicit_composition_axes_resolve_a_trace_sensitive_original_source() {
    let thermodynamics = model(0.0);
    let amounts = [0.5, 0.5, 1e-20];
    let selected = candidate(&thermodynamics, 400.0, amounts, false);
    let dependent_card = card(1.0);
    let dependent = law(&dependent_card, thermodynamics.clone());
    assert!(matches!(
        dependent.resolve_equilibrium(selected.clone(), || false),
        Err(MulticomponentLiquidError::Arithmetic(_))
    ));

    // This intentionally trace-sensitive authored curve makes dropping x_c
    // observable in the actual transport response, not only in stored metadata.
    let mut claims = ClaimSet::new();
    for (name, spec, base) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            800.0,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            1000.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            0.1,
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            1e-6,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value: PropertyValue::Curve {
                    abscissa: "x_c".into(),
                    abscissa_dims: Dimensionless::DIMS,
                    knots: vec![(0.0, base), (1e-20, 2.0 * base), (2e-20, 3.0 * base)],
                    dims: spec.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("x_a", mole_fraction(), 0.0, 1.0)
                    .with_quantity("x_b", mole_fraction(), 0.0, 1.0)
                    .with_quantity("x_c", mole_fraction(), 0.0, 2e-20)
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        600.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored trace-source resolution regression; not measured data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let explicit_card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "A+B+C".into(),
            phase: "liquid".into(),
            process: "trace-sensitive authored source".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 350.0)
        .unwrap()
        .with_quantity("x_a", mole_fraction(), 0.5)
        .unwrap()
        .with_quantity("x_b", mole_fraction(), 0.5)
        .unwrap()
        .with_quantity("x_c", mole_fraction(), 0.0)
        .unwrap();
    let reference = resolve_liquid_state_from_kinematic(
        &explicit_card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let bindings = (0..3)
        .map(|i| MulticomponentLiquidComponent {
            species: species(i),
            molar_mass: MolarMass::new(MASSES[i]),
            mole_fraction_axis: Some(["x_a", "x_b", "x_c"][i].into()),
        })
        .collect();
    let explicit = MulticomponentLiquidSolutionLaw::try_new(
        thermodynamics,
        bindings,
        &explicit_card,
        reference,
        "T",
    )
    .unwrap();
    let resolved = explicit.resolve_equilibrium(selected, || false).unwrap();
    let phase = &resolved.phases()[0];
    assert_eq!(phase.component_moles(), amounts);
    assert_eq!(phase.mole_fractions()[2], 1e-20);
    assert_eq!(
        phase.liquid().material().state_coordinate("x_c"),
        Some(1e-20)
    );
    near(explicit.reference_state().density_kg_m3(), 800.0);
    near(phase.liquid().density_kg_m3(), 1600.0);
    near(phase.liquid().specific_heat_j_kg_k(), 2000.0);
    near(phase.liquid().thermal_conductivity_w_m_k(), 0.2);
    near(phase.liquid().kinematic_viscosity_m2_s(), 2e-6);
    assert_eq!(
        phase.liquid().material().card_identity(),
        explicit_card.content_hash()
    );
}
