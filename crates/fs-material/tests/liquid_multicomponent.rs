//! G0/G1/G3/G4 liquid-source integration, using authored nonideal surface data.
//! These numerical fixtures are not measured solvent-pair fits.

#[path = "liquid_multicomponent/caloric.rs"]
mod caloric;

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
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, ResolvedLiquidState,
    resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{
    AmountConcentration, Density, Dimensionless, MolarMass, Pressure, QuantitySpec, Temperature,
};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::transport::maxwell_stefan::{MaxwellStefanMolarSolveV1, MaxwellStefanV1};

const MASSES: [f64; 3] = [0.018, 0.046, 0.100];
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn species(index: usize) -> SpeciesId {
    SpeciesId::new(["A", "B", "C"][index]).unwrap()
}
fn composition(basis: CompositionBasis) -> QuantitySpec {
    quantity(QuantityKind::Composition(basis))
}
fn card(basis: CompositionBasis, upper_a: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            [800.0, 1000.0, 1400.0, 1200.0],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [1000.0, 2000.0, 4000.0, 3000.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            [0.1, 0.2, 0.7, 0.9],
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            [1e-6, 2e-6, 3e-6, 9e-6],
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
                    .with_quantity("x_a", composition(basis), 0.0, upper_a)
                    .with_quantity("x_b", composition(basis), 0.0, 1.0)
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
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic ternary composition response; not measurements".into(),
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
            process: "authored surface".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reference(card: &MaterialCard, basis: CompositionBasis) -> ResolvedLiquidState {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 350.0)
        .unwrap()
        .with_quantity("x_a", composition(basis), 0.2)
        .unwrap()
        .with_quantity("x_b", composition(basis), 0.3)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
        .unwrap();
    resolve_liquid_state_from_kinematic(card, &point, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}
fn model(order: [usize; 3]) -> MulticomponentRegularSolutionV1 {
    let canonical = [0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    MulticomponentRegularSolutionV1::new(
        order.map(species).into(),
        order
            .into_iter()
            .flat_map(|i| order.map(|j| canonical[3 * i + j]))
            .collect(),
        vec![0.0; 9],
        Temperature::new(250.0),
        Temperature::new(650.0),
    )
    .unwrap()
}
fn components(order: [usize; 3]) -> Vec<MulticomponentLiquidComponent> {
    order
        .into_iter()
        .map(|i| MulticomponentLiquidComponent {
            species: species(i),
            molar_mass: MolarMass::new(MASSES[i]),
            mole_fraction_axis: [Some("x_a"), Some("x_b"), None][i].map(str::to_owned),
        })
        .collect()
}
fn law(card: &MaterialCard) -> MulticomponentLiquidSolutionLaw<'_> {
    MulticomponentLiquidSolutionLaw::try_new(
        model([2, 0, 1]),
        components([1, 2, 0]),
        card,
        reference(card, CompositionBasis::MoleFraction),
        "T",
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_ternary_inventory_uses_the_actual_composition_surface_and_sources() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let amounts = [2.0, 3.0, 5.0];
    let cell = law
        .homogeneous(Temperature::new(450.0), &amounts, || false)
        .unwrap();
    assert_eq!(cell.species(), [species(0), species(1), species(2)]);
    assert_eq!(cell.component_moles(), amounts);
    assert_eq!(cell.mole_fractions(), [0.2, 0.3, 0.5]);
    // At (0.2,0.3), barycentric weights are C:0.3, B:0.1, interior:0.6.
    // This is deliberately different from a pure-component mixing rule.
    near(cell.liquid().density_kg_m3(), 1100.0);
    near(cell.liquid().specific_heat_j_kg_k(), 2500.0);
    near(cell.liquid().thermal_conductivity_w_m_k(), 0.64);
    near(cell.liquid().kinematic_viscosity_m2_s(), 6e-6);
    near(cell.liquid().dynamic_viscosity_pa_s(), 1100.0 * 6e-6);
    let mass = 2.0 * 0.018 + 3.0 * 0.046 + 5.0 * 0.100;
    near(cell.mass_kg(), mass);
    near(cell.volume_m3(), mass / 1100.0);
    near(cell.molar_concentration_mol_m3(), 10.0 * 1100.0 / mass);
    for i in 0..3 {
        near(cell.mass_fractions()[i], amounts[i] * MASSES[i] / mass);
    }
    let transport = cell.transport_state().unwrap();
    let thermal = cell.thermal_transport_state().unwrap();
    assert_eq!(transport, thermal.transport());
    assert_eq!(thermal.temperature_k(), 450.0);
    assert_eq!(
        cell.liquid().material().state_coordinate("p"),
        Some(100_000.0)
    );
    for property in cell.liquid().material().properties() {
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
    assert_eq!(
        cell.liquid().material().card_identity(),
        card.content_hash()
    );
    let without_third = law
        .homogeneous(Temperature::new(450.0), &[2.0, 3.0, 0.0], || false)
        .unwrap();
    near(without_third.liquid().density_kg_m3(), 1240.0);
    assert_ne!(without_third.transport_state().unwrap(), transport);
}

#[test]
fn g1_unstable_ternary_keeps_homogeneous_chemical_forces_for_maxwell_stefan() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let cell = law
        .homogeneous(Temperature::new(400.0), &[1.0; 3], || false)
        .unwrap();
    near(
        cell.thermodynamics()
            .reduced_curvature(&[1.0, -1.0, 0.0])
            .unwrap(),
        -2.0,
    );
    assert_eq!(cell.mole_fractions(), [1.0 / 3.0; 3]);
    near(cell.liquid().density_kg_m3(), 1200.0);
    let gradient = [0.01, 0.0, -0.01];
    let driving = cell
        .thermodynamics()
        .driving_forces(&gradient, || false)
        .unwrap();
    // Cross interactions drive B even though its own concentration gradient is zero.
    assert!(driving.values()[1] > 0.0);
    let transport = MaxwellStefanV1::new(
        cell.components()
            .iter()
            .map(|c| (c.species.clone(), c.molar_mass))
            .collect(),
    )
    .unwrap();
    let mut workspace = transport.workspace().unwrap();
    let diffusivity = 1e-9;
    let pairs = [
        0.0,
        diffusivity,
        diffusivity,
        diffusivity,
        0.0,
        diffusivity,
        diffusivity,
        diffusivity,
        0.0,
    ];
    let diffusion = cell
        .thermodynamics()
        .molar_diffusion(
            &transport,
            &gradient,
            AmountConcentration::new(cell.molar_concentration_mol_m3()),
            &pairs,
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    for i in 0..3 {
        // Equal binary diffusivities have the independent exact molar-frame solution.
        near(
            diffusion.molar_fluxes()[i],
            -cell.molar_concentration_mol_m3() * diffusivity * driving.values()[i],
        );
    }
    assert!(diffusion.drag_entropy_production_w_m3_k().unwrap() >= 0.0);
}

#[test]
fn g3_species_bound_mass_and_source_axes_survive_every_input_permutation() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let reference = reference(&card, CompositionBasis::MoleFraction);
    let baseline = law(&card)
        .homogeneous(Temperature::new(400.0), &[2.0, 3.0, 5.0], || false)
        .unwrap();
    for order in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let permuted = MulticomponentLiquidSolutionLaw::try_new(
            model(order),
            components(order),
            &card,
            reference.clone(),
            "T",
        )
        .unwrap();
        assert_eq!(
            permuted
                .homogeneous(Temperature::new(400.0), &[2.0, 3.0, 5.0], || false)
                .unwrap(),
            baseline
        );
    }
    for scale in [1e-100, 1e100] {
        let scaled = law(&card)
            .homogeneous(
                Temperature::new(400.0),
                &[2.0 * scale, 3.0 * scale, 5.0 * scale],
                || false,
            )
            .unwrap();
        near(scaled.mass_kg(), baseline.mass_kg() * scale);
        near(scaled.volume_m3(), baseline.volume_m3() * scale);
        near(
            scaled.molar_concentration_mol_m3(),
            baseline.molar_concentration_mol_m3(),
        );
    }
}

#[test]
fn g0_misbound_species_and_composition_conventions_refuse() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let reference = reference(&card, CompositionBasis::MoleFraction);
    let check = |bindings| {
        MulticomponentLiquidSolutionLaw::try_new(
            model([0, 1, 2]),
            bindings,
            &card,
            reference.clone(),
            "T",
        )
        .is_err()
    };
    let valid = components([0, 1, 2]);
    for variant in 0..7 {
        let mut bad = valid.clone();
        match variant {
            0 => bad[0].species = SpeciesId::new("D").unwrap(),
            1 => bad[1].species = bad[0].species.clone(),
            2 => bad[1].mole_fraction_axis = Some("x_a".into()),
            3 => bad[0].mole_fraction_axis = None,
            4 => bad[2].mole_fraction_axis = Some("x_c".into()),
            5 => bad[0].mole_fraction_axis = Some("T".into()),
            _ => bad[0].molar_mass = MolarMass::new(f64::NAN),
        }
        assert!(check(bad), "variant {variant}");
    }
    assert!(check(valid[..2].to_vec()));
    let redundant_point = reference
        .material()
        .query()
        .with_quantity("x_c", composition(CompositionBasis::MoleFraction), 0.5)
        .unwrap();
    let redundant = resolve_liquid_state_from_kinematic(
        &card,
        &redundant_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        MulticomponentLiquidSolutionLaw::try_new(
            model([0, 1, 2]),
            valid.clone(),
            &card,
            redundant,
            "T",
        )
        .is_err()
    );
    let wrong_card = self::card(CompositionBasis::MassFraction, 1.0);
    assert!(
        MulticomponentLiquidSolutionLaw::try_new(
            model([0, 1, 2]),
            valid.clone(),
            &wrong_card,
            reference.clone(),
            "T"
        )
        .is_err()
    );
    assert!(
        MulticomponentLiquidSolutionLaw::try_new(
            model([0, 1, 2]),
            valid,
            &wrong_card,
            self::reference(&wrong_card, CompositionBasis::MassFraction),
            "T"
        )
        .is_err()
    );
}
#[test]
fn g0_pure_limits_and_explicit_traces_do_not_become_invented_components() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    for amounts in [
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [1e-20, 0.0, 1.0],
    ] {
        let cell = law
            .homogeneous(Temperature::new(400.0), &amounts, || false)
            .unwrap();
        assert_eq!(cell.component_moles(), amounts);
        for i in 0..3 {
            assert_eq!(cell.mole_fractions()[i] > 0.0, amounts[i] > 0.0);
            assert_eq!(cell.mass_fractions()[i] > 0.0, amounts[i] > 0.0);
            assert_eq!(
                cell.thermodynamics().reduced_chemical_potentials()[i].is_some(),
                amounts[i] > 0.0
            );
        }
    }
    // C has no explicit source coordinate: this positive trace cannot survive
    // the independent x_a+x_b chart and must not be silently published as absent.
    assert!(matches!(
        law.homogeneous(Temperature::new(400.0), &[0.5, 0.5, 1e-20], || false),
        Err(MulticomponentLiquidError::Arithmetic(_))
    ));
    // Rounded independent coordinates can also imply a positive C that is
    // absent or orders of magnitude larger than the actual trace. Absolute
    // epsilon tolerances would hide both source/thermodynamics mismatches.
    for dependent in [0.0, 1e-20] {
        assert!(matches!(
            law.homogeneous(
                Temperature::new(400.0),
                &[0.9334141514912707, 0.13043915609467183, dependent],
                || false,
            ),
            Err(MulticomponentLiquidError::Arithmetic(_))
        ));
    }
    for amounts in [
        [0.0; 3],
        [-1.0, 1.0, 1.0],
        [f64::NAN, 1.0, 1.0],
        [f64::MAX; 3],
        [f64::MIN_POSITIVE, 1e308, 0.0],
    ] {
        assert!(
            law.homogeneous(Temperature::new(400.0), &amounts, || false)
                .is_err()
        );
    }
    assert!(
        law.homogeneous(Temperature::new(400.0), &[1.0; 2], || false)
            .is_err()
    );
}

#[test]
fn g4_unsupported_sources_and_every_cancellation_leave_inputs_unchanged() {
    let card = card(CompositionBasis::MoleFraction, 0.8);
    let law = law(&card);
    let original = law.reference_state().clone();
    let amounts = [2.0, 3.0, 5.0];
    let mut polls = 0;
    let completed = law
        .homogeneous(Temperature::new(400.0), &amounts, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut current = 0;
        assert!(
            law.homogeneous(Temperature::new(400.0), &amounts, || {
                current += 1;
                current == stop
            })
            .is_err()
        );
    }
    assert!(matches!(
        law.homogeneous(Temperature::new(400.0), &[9.0, 0.0, 1.0], || false),
        Err(MulticomponentLiquidError::Source(_))
    ));
    assert!(matches!(
        law.homogeneous(Temperature::new(625.0), &amounts, || false),
        Err(MulticomponentLiquidError::Source(_))
    ));
    assert!(matches!(
        law.homogeneous(Temperature::new(700.0), &amounts, || false),
        Err(MulticomponentLiquidError::Thermodynamics(_))
    ));
    assert_eq!(law.reference_state(), &original);
    assert_eq!(amounts, [2.0, 3.0, 5.0]);
    assert_eq!(
        law.homogeneous(Temperature::new(400.0), &amounts, || false)
            .unwrap(),
        completed
    );
}
