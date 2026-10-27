//! Synthetic balanced isomers exercise the real kinetics, caloric and momentum
//! owners. No empirical combustion-rate or engine validation claim is made.
mod thermal;
mod network;
mod implicit;
mod arrhenius;

use super::*;
use crate::cavity_phs::piston::material::MaterialGasPiston;
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::SpeciesTransportMaterial;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, DynViscosity, MolarMass, Pressure, QuantitySpec};
use fs_thermochem::equilibrium::{
    IdealGasReactionEquilibriumV1, MassActionKineticOrderConventionV1,
    ReactionProgressRateQuantityV1, ReactionProgressRateScaleV1, ThermodynamicReverseRateClosureV1,
};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReactionId, ReferenceEquationOfStateV1, StandardStateConventionV1, StandardStatePhaseV1,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
    verify_conservation,
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic reactive-gas source fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn species(
    name: &str,
    molar_mass: f64,
    cp_r: f64,
    formation_j_mol: f64,
) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[7] = formation_j_mol / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: provenance(),
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(molar_mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic-isomers").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn transport(name: &str, factor: f64) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (name, dims, base) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, factor * base), (2000.0, 2.0 * factor * base)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn catalog(release: f64, product_cp_r: f64, product_mass: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, 0.0),
                species("B", product_mass, product_cp_r, -release),
                species("I", 0.004, 2.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.0),
            transport("B", 1.5),
            transport("I", 0.7),
        ],
        None,
    )
    .unwrap()
}
fn kinetics(catalog: &ThermochemicalGasModel, scale: f64) -> DeclaredStoichiometricMassActionV1 {
    kinetics_with_atoms(catalog, scale, 1)
}
fn kinetics_with_atoms(
    catalog: &ThermochemicalGasModel,
    scale: f64,
    product_atoms: i128,
) -> DeclaredStoichiometricMassActionV1 {
    let species = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let element = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, product_atoms]],
    )
    .unwrap();
    let reaction = ReactionId::new("synthetic-conserved-reaction").unwrap();
    let matrix = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-product_atoms], vec![1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&element, &matrix, &charge).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        matrix,
        certificate,
        reaction,
        catalog.thermodynamics().components()[..2].to_vec(),
    )
    .unwrap();
    DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(scale)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    )
}
fn reaction(release: f64, cp_r: f64, rate: f64) -> GasReaction {
    let catalog = catalog(release, cp_r, 0.028);
    let kinetics = kinetics(&catalog, rate);
    GasReaction::new(catalog, kinetics).unwrap()
}
fn gas(reaction: &GasReaction, fractions: [f64; 3]) -> ClosedGasParcel {
    let model = reaction
        .catalog()
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        )
        .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            600.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, 0.001).unwrap()
}
fn chemical_controls() -> GasReactionSolve {
    GasReactionSolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        maximum_consumed_fraction: 0.2,
        species_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
    }
}
fn mechanical_controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.0007, 0.0015],
        gas: chemical_controls().gas,
        momentum_tolerance_n_s: 1e-11,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn drive() -> PistonInput {
    PistonInput {
        duration_s: 0.00005,
        heat_j: 0.0,
        exterior_pressure_pa: 100000.0,
        external_force_n: 0.0,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} within {tolerance}"
    );
}

#[test]
fn g1_formation_energy_creates_temperature_pressure_and_products_without_extra_heat() {
    let reaction = reaction(4000.0, 4.0, 1000.0);
    let before = gas(&reaction, [0.75, 0.0, 0.25]);
    let candidate = reaction
        .parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    let after = candidate.gas();
    let report = candidate.step();
    near(
        report.extent_mol,
        1000.0 * 0.75 * 0.001 * drive().duration_s,
        1e-17,
    );
    assert_eq!(report.species[1].initial_mol, 0.0);
    assert!(report.species[1].final_mol > 0.0);
    let n: Vec<_> = report.species.iter().map(|s| s.final_mol).collect();
    let before_u = before.mass_kg() * before.state().specific_internal_energy_j_kg();
    let expected_t = (before_u + n[1] * 4000.0) / (R * (2.5 * n[0] + 3.0 * n[1] + 1.5 * n[2]));
    near(after.state().state().temperature, expected_t, 1e-7);
    near(
        after.state().state().pressure,
        n.iter().sum::<f64>() * R * expected_t / after.volume_m3(),
        1e-6,
    );
    assert!(after.state().state().temperature > before.state().state().temperature);
    assert!(report.entropy_change_j_k > 0.0);
    assert_eq!(after.mass_kg(), before.mass_kg());
    assert_eq!(after.volume_m3(), before.volume_m3());
    assert_eq!(after.point().axes().get("specimen-lot"), Some(&42.0));
    assert!(
        report.stored_energy_change_j.abs()
            <= before.mass_kg() * chemical_controls().gas.energy_tolerance_j_kg
    );
    assert_eq!(
        after.model().thermodynamics().components(),
        reaction.catalog().thermodynamics().components()
    );
    assert_ne!(
        after.state().state().sound_speed,
        before.state().state().sound_speed
    );
    assert_ne!(
        after.state().state().dynamic_viscosity,
        before.state().state().dynamic_viscosity
    );
    assert_ne!(
        after.state().state().thermal_conductivity,
        before.state().state().thermal_conductivity
    );
}

#[test]
fn g1_molecule_count_changes_pressure_without_inventing_heat_or_mass() {
    // 2 A -> B has equal total molar cv on each side, and zero formation
    // offset. Thus T stays fixed while pressure follows the changing N/V.
    let catalog = catalog(0.0, 6.0, 0.056);
    let kinetics = kinetics_with_atoms(&catalog, 1000.0, 2);
    let reaction = GasReaction::new(catalog, kinetics).unwrap();
    let initial = gas(&reaction, [1.0, 0.0, 0.0]);
    let total = initial.mass_kg() / 0.028;
    let candidate = reaction
        .parcel_candidate(&initial, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    let step = candidate.step();
    near(step.extent_mol, 1000.0 * 0.001 * drive().duration_s, 1e-17);
    near(
        step.species[0].final_mol,
        total - 2.0 * step.extent_mol,
        1e-14,
    );
    near(step.species[1].final_mol, step.extent_mol, 1e-14);
    assert_eq!(candidate.gas().mass_kg(), initial.mass_kg());
    near(candidate.gas().state().state().temperature, 600.0, 1e-7);
    near(
        candidate.gas().state().state().pressure,
        100000.0 * (total - step.extent_mol) / total,
        1e-6,
    );
}

#[test]
fn g3_reverse_reaction_cools_and_absent_reagents_do_not_need_trace_seeds() {
    let reaction = reaction(4000.0, 3.5, 1000.0);
    let products = gas(&reaction, [0.0, 1.0, 0.0]);
    let candidate = reaction
        .parcel_candidate(&products, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    assert!(candidate.step().extent_mol < 0.0);
    assert!(candidate.gas().state().state().temperature < products.state().state().temperature);
    assert!(candidate.step().species[0].final_mol > 0.0);
    let inert = gas(&reaction, [0.0, 0.0, 1.0]);
    let candidate = reaction
        .parcel_candidate(&inert, drive().duration_s, chemical_controls(), || false)
        .unwrap();
    assert_eq!(candidate.gas(), &inert);
    assert_eq!(candidate.step().extent_mol, 0.0);
    assert_eq!(candidate.step().stored_energy_change_j, 0.0);
    assert_eq!(candidate.step().entropy_change_j_k, 0.0);
}

#[test]
fn g1_constant_temperature_isomer_extent_converges_at_first_order() {
    let reaction = reaction(0.0, 3.5, 1000.0);
    let initial = gas(&reaction, [1.0, 0.0, 0.0]);
    let amount = initial.mass_kg() / 0.028;
    let duration = 0.002;
    let exact_b =
        0.5 * amount * (1.0 - (-2.0 * 1000.0 * initial.volume_m3() * duration / amount).exp());
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let mut gas = initial.clone();
        let mut b = 0.0;
        for _ in 0..steps {
            let candidate = reaction
                .parcel_candidate(&gas, duration / steps as f64, chemical_controls(), || false)
                .unwrap();
            b = candidate.step().species[1].final_mol;
            gas = candidate.into_parts().0;
        }
        near(gas.state().state().temperature, 600.0, 1e-7);
        errors.push((b - exact_b).abs());
    }
    assert!(
        errors
            .windows(2)
            .all(|e| e[0] / e[1] > 1.95 && e[0] / e[1] < 2.05),
        "{errors:?}"
    );
}

fn body(density: f64) -> fs_mbd_material::ResolvedBody {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, density), (400.0, density)],
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "specimen".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    resolve_body(
        &[RegionBinding {
            key: "wall",
            material: &card,
            point: &point,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}

#[test]
fn g3_data_only_chemistry_dilution_and_solid_inertia_drive_conservative_motion() {
    let run = |release, rate, dilution, density| {
        let law = reaction(release, 3.5, rate);
        let mut piston = MaterialGasPiston::new(
            gas(&law, [1.0 - dilution, 0.0, dilution]),
            body(density),
            0.01,
            0.0,
        )
        .unwrap();
        let original = piston.clone();
        let mut work = 0.0;
        for _ in 0..64 {
            let candidate = piston
                .reactive_candidate(
                    &law,
                    drive(),
                    chemical_controls(),
                    mechanical_controls(),
                    || false,
                )
                .unwrap();
            let report = candidate.step();
            assert!(
                report.total_energy_residual_j.abs() <= mechanical_controls().energy_tolerance_j
            );
            assert!(
                report.reaction.maximum_species_residual_mol
                    <= chemical_controls().species_tolerance_mol
            );
            work += report.piston.gas.work_by_gas_j;
            piston = candidate.into_parts().0;
        }
        assert_eq!(piston.body(), original.body());
        assert_eq!(piston.gas().mass_kg(), original.gas().mass_kg());
        let gas_change = piston.gas().mass_kg()
            * (piston.gas().state().specific_internal_energy_j_kg()
                - original.gas().state().specific_internal_energy_j_kg());
        let kinetic = 0.5 * piston.body().inertia().mass() * piston.piston().velocity_m_s().powi(2);
        let exterior =
            drive().exterior_pressure_pa * (piston.gas().volume_m3() - original.gas().volume_m3());
        near(
            gas_change + kinetic + exterior,
            0.0,
            64.0 * mechanical_controls().energy_tolerance_j,
        );
        (piston, work)
    };
    let (base, work) = run(4000.0, 1000.0, 0.0, 8000.0);
    let (energetic, _) = run(8000.0, 1000.0, 0.0, 8000.0);
    let (faster, _) = run(4000.0, 2000.0, 0.0, 8000.0);
    let (dilute, _) = run(4000.0, 1000.0, 0.5, 8000.0);
    let (light, light_work) = run(4000.0, 1000.0, 0.0, 1000.0);
    assert!(work > 0.0);
    assert!(base.piston().velocity_m_s() > 0.0);
    assert!(energetic.piston().velocity_m_s() > base.piston().velocity_m_s());
    assert!(faster.piston().velocity_m_s() > base.piston().velocity_m_s());
    assert!(dilute.piston().velocity_m_s() < base.piston().velocity_m_s());
    assert!(light.piston().velocity_m_s() > base.piston().velocity_m_s());
    assert!(light_work > work);
}

#[test]
fn g3_source_conflict_mass_mismatch_and_unresolved_consumption_refuse() {
    let law = reaction(4000.0, 3.5, 1000.0);
    let before = gas(&law, [1.0, 0.0, 0.0]);
    let other = reaction(5000.0, 3.5, 1000.0);
    // Absent products still have one retained source catalog: source mismatch
    // cannot be hidden until the two catalogs become active in the same gas.
    assert!(GasReaction::new(other.catalog().clone(), law.kinetics().clone()).is_err());
    let changed_a = catalog(4000.0, 3.5, 0.03);
    assert!(GasReaction::new(changed_a.clone(), kinetics(&changed_a, 1000.0)).is_err());
    let mixture = gas(&law, [0.5, 0.5, 0.0]);
    assert!(
        other
            .parcel_candidate(&mixture, drive().duration_s, chemical_controls(), || false)
            .is_err()
    );
    assert!(matches!(
        law.parcel_candidate(&before, 1.0, chemical_controls(), || false),
        Err(GasReactionError::ConsumptionLimit { .. })
    ));
    let unsupported = GasReactionSolve {
        gas: GasEnergySolve {
            maximum_temperature_k: 600.01,
            ..chemical_controls().gas
        },
        ..chemical_controls()
    };
    assert!(
        law.parcel_candidate(&before, drive().duration_s, unsupported, || false)
            .is_err()
    );
    assert_eq!(before, gas(&law, [1.0, 0.0, 0.0]));
}

#[test]
fn g4_early_and_late_cancellation_preserve_all_accepted_owners_and_ports() {
    let law = reaction(4000.0, 3.5, 1000.0);
    let piston =
        MaterialGasPiston::new(gas(&law, [1.0, 0.0, 0.0]), body(8000.0), 0.01, 0.0).unwrap();
    let before = piston.clone();
    let mut count = 0;
    piston
        .reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    assert!(count > 4);
    for cancel_at in [1, count] {
        let mut calls = 0;
        let result = piston.reactive_candidate(
            &law,
            drive(),
            chemical_controls(),
            mechanical_controls(),
            || {
                calls += 1;
                calls == cancel_at
            },
        );
        assert!(matches!(result, Err(GasReactionError::Cancelled)));
        assert_eq!(piston, before);
        assert_eq!(
            piston.piston().ports.next_tick(),
            before.piston().ports.next_tick()
        );
    }
}

#[test]
fn material_bound_fixed_reactor_changes_sound_and_transport_without_wall_work() {
    use crate::cavity_phs::piston::PistonSupport;
    let run = |release, density| {
        let law = reaction(release, 3.5, 1000.0);
        let mut wall = MaterialGasPiston::new(
            gas(&law, [1.0, 0.0, 0.0]), body(density), 0.01, 0.0,
        ).unwrap().with_fixed_support().unwrap();
        let original = wall.clone();
        let controls = PistonSolve {
            volume_bracket_m3: [wall.gas().volume_m3(); 2], ..mechanical_controls()
        };
        for _ in 0..16 {
            let candidate = wall.reactive_candidate(
                &law, drive(), chemical_controls(), controls, || false,
            ).unwrap();
            let report = candidate.step();
            assert!(report.reaction.extent_mol > 0.0);
            assert!(report.total_energy_residual_j.abs() <= controls.energy_tolerance_j);
            assert!(report.piston.support_impulse_n_s < 0.0);
            assert_eq!(report.piston.support_work_j, 0.0);
            assert_eq!(report.piston.gas.supplied_heat_j, 0.0);
            assert_eq!(report.piston.gas.work_by_gas_j, 0.0);
            wall = candidate.into_parts().0;
            assert_eq!(wall.piston().support(), PistonSupport::Fixed);
            assert_eq!(wall.piston().velocity_m_s(), 0.0);
            assert_eq!(wall.gas().volume_m3().to_bits(), original.gas().volume_m3().to_bits());
            assert_eq!(wall.body(), original.body());
        }
        assert_eq!(wall.gas().mass_kg(), original.gas().mass_kg());
        near(wall.gas().mass_kg() * (wall.gas().state().specific_internal_energy_j_kg()
            - original.gas().state().specific_internal_energy_j_kg()), 0.0,
            16.0 * controls.energy_tolerance_j);
        assert!(wall.gas().state().state().temperature > original.gas().state().state().temperature);
        assert!(wall.gas().state().state().sound_speed > original.gas().state().state().sound_speed);
        assert_ne!(wall.gas().state().state().dynamic_viscosity, original.gas().state().state().dynamic_viscosity);
        wall
    };
    let base = run(4000.0, 8000.0);
    let light = run(4000.0, 1000.0);
    assert_ne!(base.body().inertia().mass(), light.body().inertia().mass());
    assert_eq!(base.gas(), light.gas());
    let energetic = run(8000.0, 8000.0);
    assert!(energetic.gas().state().state().temperature > base.gas().state().state().temperature);
    assert!(energetic.gas().state().state().sound_speed > base.gas().state().state().sound_speed);
}

