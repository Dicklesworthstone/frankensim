//! Manufactured exact-chemistry salt/acid equilibria, not measured solubility.

mod stoichiometry;

use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::acid_base::{
            AcidDissociation, AcidFamily, SpeciationOptions, WaterAutoprotolysis,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};

const REFERENCE: &str = "source-elemental-reference-authored-salt-test";
const T: f64 = 298.15;
const KA: f64 = 1e-7;
const KW: f64 = 1e-14;
const H: &str = "H+(aq)";
const OH: &str = "OH[-](aq)";
const WATER: &str = "H2O(l)";
const METAL: &str = "Fe[3+](aq)";
const FORMS: [&str; 4] = [
    "H3AsO4(aq)", "H2AsO4[-](aq)", "HAsO4[2-](aq)", "AsO4[3-](aq)",
];
const SOLID: &str = "FeAsO4.2H2O(s)";
const DRY_SOLID: &str = "authored-FeAsO4(s)";
const FAMILY: &str = "arsenate";

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn never_cancelled() -> bool {
    false
}

fn point(phase: &str, charge: i128) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static,
            )),
            T,
        ).unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        ).unwrap();
    let mut flags = vec![REFERENCE, "source-standard-state", phase];
    if phase == "source-phase-aqueous" {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
        if charge != 0 {
            flags.extend(["source-individual-ion", "source-ionic-reference-hydrogen-zero"]);
        }
    } else {
        flags.push("source-pure-species");
    }
    for flag in flags {
        at = at.with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0).unwrap();
    }
    at
}

fn pack(name: &str, gibbs: f64, point: &QueryPoint) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in point.axes() {
        validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
    }
    let dims = MolarEnergyQuantityV1::DIMS;
    let mut claim = PropertyClaim {
        key: PropertyKey::new("standard-molar-formation-gibbs-energy", dims),
        value: PropertyValue::Scalar { value: gibbs, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("authored salt equilibrium fixture: {name}"),
            license: "project fixture".into(),
            artifact: None,
        },
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    claim.observations.push(claims.register_observation(ObservationDataset {
        specimen: name.into(),
        method: "independent analytical equilibrium construction".into(),
        artifact: hash,
        caveats: "not measured iron/arsenic thermochemistry".into(),
        provenance: claim.provenance.clone(),
    }).unwrap());
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        name, "fixture-v1", hash, "project fixture", claims, vec![], vec![],
    ).unwrap()
}

#[derive(Clone)]
struct Record {
    name: &'static str,
    // H, O, As, Fe, C, N. Exact atom/charge identities remain independently checked.
    atoms: [i128; 6],
    charge: i128,
    gibbs: f64,
    point: QueryPoint,
    pack: NormalizedPack,
}
impl Record {
    fn new(
        name: &'static str, atoms: [i128; 6], charge: i128, gibbs: f64, phase: &str,
    ) -> Self {
        let point = point(phase, charge);
        let pack = pack(name, gibbs, &point);
        Self { name, atoms, charge, gibbs, point, pack }
    }
}

#[derive(Clone)]
struct Fixture {
    records: Vec<Record>,
}
impl Fixture {
    fn new(kdiss: f64) -> Self {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T;
        let step = -rt * KA.ln();
        let solid = 3.0 * step + rt * kdiss.ln();
        let aqueous = "source-phase-aqueous";
        let mut records = vec![
            Record::new(H, [1, 0, 0, 0, 0, 0], 1, 0.0, aqueous),
            Record::new(OH, [1, 1, 0, 0, 0, 0], -1, -rt * KW.ln(), aqueous),
            Record::new(WATER, [2, 1, 0, 0, 0, 0], 0, 0.0, "source-phase-liquid"),
            Record::new(METAL, [0, 0, 0, 1, 0, 0], 3, 0.0, aqueous),
            Record::new(SOLID, [4, 6, 1, 1, 0, 0], 0, solid, "source-phase-solid"),
            Record::new(DRY_SOLID, [0, 4, 1, 1, 0, 0], 0, solid, "source-phase-solid"),
            Record::new("HCN(aq)", [1, 0, 0, 0, 1, 1], 0, 0.0, aqueous),
            Record::new("CN[-](aq)", [0, 0, 0, 0, 1, 1], -1, step, aqueous),
        ];
        for (j, name) in FORMS.into_iter().enumerate() {
            records.push(Record::new(
                name, [3 - j as i128, 4, 1, 0, 0, 0],
                -(j as i128), j as f64 * step, aqueous,
            ));
        }
        Self { records }
    }

    fn index(&self, name: &str) -> usize {
        self.records.iter().position(|record| record.name == name).unwrap()
    }

    fn sources<'a>(&'a self, reaction: &'a ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
        reaction.terms().iter().map(|term| {
            let record = &self.records[self.index(term.species.as_str())];
            FormationGibbsSource {
                species: &term.species, pack: &record.pack,
                at: &record.point, reference_axis: REFERENCE,
            }
        }).collect()
    }

    fn reaction(&self, name: &str, terms: &[(&str, i128)]) -> ConservedReaction {
        let records: Vec<_> = terms.iter()
            .map(|(name, _)| &self.records[self.index(name)]).collect();
        let species: Vec<_> = records.iter().map(|record| id(record.name)).collect();
        let reaction = ReactionId::new(name).unwrap();
        let elemental = ElementalMatrix::new(
            ["H", "O", "As", "Fe", "C", "N"].into_iter()
                .map(|element| ElementId::new(element).unwrap()).collect(),
            species.clone(),
            (0..6).map(|row| records.iter().map(|record| record.atoms[row]).collect()).collect(),
        ).unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            species.clone(), vec![reaction.clone()],
            terms.iter().map(|(_, coefficient)| vec![*coefficient]).collect(),
        ).unwrap();
        let charges = ChargeVector::new(
            species, records.iter().map(|record| record.charge).collect(),
        ).unwrap();
        ConservedReaction::new(&elemental, &stoichiometric, &charges, reaction).unwrap()
    }

    fn aqueous(&self, buffered: bool) -> IdealAcidBaseSystem {
        let reaction = MineralReaction::WaterDissociation.build().unwrap();
        let water = WaterAutoprotolysis::from_sources(
            reaction.clone(), id(WATER), id(H), id(OH), &self.sources(&reaction),
        ).unwrap();
        let mut steps = Vec::new();
        for (j, kind) in [
            MineralReaction::ArsenicAcidDissociation,
            MineralReaction::DihydrogenArsenateDissociation,
            MineralReaction::HydrogenArsenateDissociation,
        ].into_iter().enumerate() {
            let reaction = kind.build().unwrap();
            steps.push(AcidDissociation::from_sources(
                reaction.clone(), id(FORMS[j]), id(FORMS[j + 1]), id(H),
                &self.sources(&reaction),
            ).unwrap());
        }
        let mut families = vec![AcidFamily::new(FAMILY, steps).unwrap()];
        if buffered {
            let reaction = MineralReaction::HydrogenCyanideDissociation.build().unwrap();
            families.push(AcidFamily::new("buffer", vec![
                AcidDissociation::from_sources(
                    reaction.clone(), id("HCN(aq)"), id("CN[-](aq)"), id(H),
                    &self.sources(&reaction),
                ).unwrap(),
            ]).unwrap());
        }
        IdealAcidBaseSystem::new(water, families).unwrap()
    }

    fn model(&self, buffered: bool) -> SaltPrecipitation {
        let reaction = MineralReaction::ScoroditePrecipitation.build().unwrap();
        SaltPrecipitation::from_sources(
            self.aqueous(buffered), reaction.clone(), id(METAL), id(SOLID), FAMILY,
            &self.sources(&reaction),
        ).unwrap()
    }
}

fn options() -> PrecipitationOptions {
    PrecipitationOptions {
        aqueous: SpeciationOptions {
            ph_tolerance: Dimensionless::new(1e-13),
            absolute_charge_tolerance: Amount::new(1e-26),
            relative_charge_tolerance: Dimensionless::new(1e-12),
            max_iterations: 256,
            ..SpeciationOptions::default()
        },
        ..PrecipitationOptions::default()
    }
}
fn ion(name: &str, charge: i32, amount: f64) -> SpectatorIon {
    SpectatorIon { species: id(name), charge, amount: Amount::new(amount) }
}
fn total(name: &str, amount: f64) -> (String, Amount) {
    (name.into(), Amount::new(amount))
}
fn counterions(net: f64) -> Vec<SpectatorIon> {
    if net > 0.0 {
        vec![ion("Cl[-](aq)", -1, net)]
    } else if net < 0.0 {
        vec![ion("Na[+](aq)", 1, -net)]
    } else {
        vec![]
    }
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "actual {actual:.17e}, expected {expected:.17e}",
    );
}
fn solve(
    model: &SaltPrecipitation, metal: f64, family: f64, spectators: &[SpectatorIon],
) -> SaltSpeciation {
    model.speciate(
        Amount::new(metal), &[total(FAMILY, family)], spectators,
        Amount::new(1.0), options(), never_cancelled,
    ).unwrap()
}
fn closed(result: &SaltSpeciation) {
    assert!(
        result.cation_balance_residual.value().abs()
            <= 64.0 * f64::EPSILON * result.cation_total.value()
    );
    assert!(
        result.family_balance_residual.value().abs()
            <= 64.0 * f64::EPSILON * result.family_total.value()
    );
    let family = result.aqueous.families.iter().find(|f| f.id == FAMILY).unwrap();
    near(
        family.forms.iter().map(|form| form.amount.value()).sum(),
        result.dissolved_family.value(), 3e-14, 0.0,
    );
    assert!(result.aqueous.max_log_mass_action_residual.value() < 1e-11);
    if result.status == PrecipitationStatus::SolidPresent {
        assert!(
            result.log_saturation_ratio.unwrap().value().abs()
                + result.log_saturation_numerical_spread.unwrap().value()
                <= options().log_saturation_tolerance.value()
        );
    }
}

#[test]
fn analytical_four_form_salt_root_closes_each_component_and_shared_charge() {
    // h=Ka=1e-7 makes all four fractions 1/4 and mean family charge -3/2.
    // Free Fe=1e-5 and dissolved As=2e-5 are then neutral without spectators.
    // Kdiss=(1e-5)*(2e-5/4)=5e-11, with a separately chosen 1e-4 solid.
    let model = Fixture::new(5e-11).model(false);
    let result = solve(&model, 1.1e-4, 1.2e-4, &[]);
    assert_eq!(result.status, PrecipitationStatus::SolidPresent);
    near(result.dissolved_cation.value(), 1e-5, 3e-9, 0.0);
    near(result.dissolved_family.value(), 2e-5, 3e-9, 0.0);
    near(result.solid_amount.value(), 1e-4, 3e-9, 0.0);
    near(result.aqueous.ph.value(), 7.0, 0.0, 2e-9);
    for form in &result.aqueous.families[0].forms {
        near(form.fraction.value(), 0.25, 3e-9, 0.0);
        near(form.amount.value(), 5e-6, 4e-9, 0.0);
    }
    assert!(result.aqueous.charge_balance_residual.value().abs() < 1e-16);
    assert!(
        result.dissolved_limiting_bracket[0].value() <= result.dissolved_cation.value()
            && result.dissolved_cation.value() <= result.dissolved_limiting_bracket[1].value()
    );
    assert!(result.aqueous_solves >= result.iterations + 2);
    assert!(result.aqueous_iterations > result.aqueous_solves);
    assert_eq!(result.model().family_id(), FAMILY);
    assert_eq!(result.model().anion(), &id(FORMS[3]));
    assert_eq!(result.model().cation(), &id(METAL));
    assert_eq!(result.model().solid(), &id(SOLID));
    assert_eq!(result.model().reaction_network().reactions().len(), 5);
    assert_eq!(result.model().standard().terms().len(), 4);
    closed(&result);
}

#[test]
fn either_component_can_limit_and_equal_totals_keep_the_same_formula_basis() {
    for (metal, family) in [(1e-5, 2e-5), (2e-5, 1e-5), (1e-5, 1e-5)] {
        let model = Fixture::new(metal * family / 4.0).model(false);
        let result = solve(
            &model, metal + 1e-4, family + 1e-4,
            &counterions(3.0 * metal - 1.5 * family),
        );
        near(result.dissolved_cation.value(), metal, 4e-9, 0.0);
        near(result.dissolved_family.value(), family, 4e-9, 0.0);
        near(result.solid_amount.value(), 1e-4, 4e-9, 0.0);
        near(result.aqueous.ph.value(), 7.0, 0.0, 3e-9);
        closed(&result);
    }
}

#[test]
fn another_buffer_shares_the_recomputed_ph_and_input_order_is_irrelevant() {
    let model = Fixture::new(5e-11).model(true);
    let mut families = vec![total("buffer", 3e-5), total(FAMILY, 1.2e-4)];
    let mut spectators = vec![
        ion("Na[+](aq)", 1, 1.5e-5),
        ion("K[+](aq)", 1, 4e-6),
        ion("Cl[-](aq)", -1, 4e-6),
    ];
    let mut before_ions = spectators.clone();
    before_ions.push(ion(METAL, 3, 1.1e-4));
    let before = model.aqueous_model().speciate(
        &families, &before_ions, Amount::new(1.0), options().aqueous, never_cancelled,
    ).unwrap();
    let first = model.speciate(
        Amount::new(1.1e-4), &families, &spectators,
        Amount::new(1.0), options(), never_cancelled,
    ).unwrap();
    assert!(before.ph.value() > first.aqueous.ph.value() + 0.1);
    near(first.aqueous.ph.value(), 7.0, 0.0, 3e-9);
    let buffer = first.aqueous.families.iter().find(|f| f.id == "buffer").unwrap();
    for form in &buffer.forms {
        near(form.amount.value(), 1.5e-5, 4e-9, 0.0);
    }
    families.reverse();
    spectators.reverse();
    let second = model.speciate(
        Amount::new(1.1e-4), &families, &spectators,
        Amount::new(1.0), options(), never_cancelled,
    ).unwrap();
    assert_eq!(first.dissolved_cation, second.dissolved_cation);
    assert_eq!(first.dissolved_family, second.dissolved_family);
    assert_eq!(first.aqueous.ph, second.aqueous.ph);
    assert_eq!(first.log_saturation_ratio, second.log_saturation_ratio);
    closed(&first);
}

#[test]
fn scaling_all_extensive_amounts_and_basis_preserves_activities_and_partition() {
    let model = Fixture::new(5e-11).model(false);
    for scale in [1e-6, 1.0, 1e6] {
        let result = model.speciate(
            Amount::new(scale * 1.1e-4), &[total(FAMILY, scale * 1.2e-4)], &[],
            Amount::new(scale), options(), never_cancelled,
        ).unwrap();
        near(result.dissolved_cation.value() / scale, 1e-5, 4e-9, 0.0);
        near(result.dissolved_family.value() / scale, 2e-5, 4e-9, 0.0);
        near(result.solid_amount.value() / scale, 1e-4, 4e-9, 0.0);
        near(result.aqueous.ph.value(), 7.0, 0.0, 3e-9);
        closed(&result);
    }
}

#[test]
fn empty_and_absent_phases_preserve_the_entire_other_component() {
    let model = Fixture::new(5e-11).model(false);
    for (metal, family) in [(0.0, 0.0), (0.0, 0.1), (0.1, 0.0)] {
        let result = solve(
            &model, metal, family, &counterions(3.0 * metal - 1.5 * family),
        );
        assert_eq!(result.status, PrecipitationStatus::EmptyComponent);
        assert_eq!(result.dissolved_cation.value(), metal);
        assert_eq!(result.dissolved_family.value(), family);
        assert_eq!(result.solid_amount.value(), 0.0);
        assert_eq!(result.log_saturation_ratio, None);
        assert_eq!(result.log_saturation_numerical_spread, None);
        assert_eq!(result.aqueous_solves, 1);
        closed(&result);
    }
    let absent = solve(&model, 1e-9, 1e-9, &[]);
    assert_eq!(absent.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(absent.dissolved_cation.value(), 1e-9);
    assert_eq!(absent.dissolved_family.value(), 1e-9);
    assert_eq!(absent.solid_amount.value(), 0.0);
    assert!(
        absent.log_saturation_ratio.unwrap().value()
            + absent.log_saturation_numerical_spread.unwrap().value() < 0.0
    );
    assert_eq!(absent.iterations, 0);
    closed(&absent);
}

#[test]
fn logarithmic_partition_preserves_trace_in_either_limiting_component() {
    let model = Fixture::new(2.5e-41).model(false);
    for (metal, family) in [(1e-40, 1.0), (1.0, 1e-40)] {
        let result = solve(
            &model, metal + 1.0, family + 1.0,
            &counterions(3.0 * metal - 1.5 * family),
        );
        assert_eq!(result.status, PrecipitationStatus::SolidPresent);
        near(result.dissolved_cation.value(), metal, 4e-9, 0.0);
        near(result.dissolved_family.value(), family, 4e-9, 0.0);
        assert_eq!(result.solid_amount.value(), 1.0);
        let residual = if metal < family {
            result.cation_balance_residual
        } else {
            result.family_balance_residual
        };
        near(residual.value(), 1e-40, 4e-9, 0.0);
        assert!(result.dissolved_limiting_bracket[0].value() > 0.0);
        closed(&result);
    }
}

#[test]
fn optional_hydrate_water_may_be_absent_but_never_changes_formula_unit_counts() {
    let fixture = Fixture::new(5e-11);
    let reaction = fixture.reaction(
        "authored-anhydrous-salt",
        &[(METAL, -1), (FORMS[3], -1), (DRY_SOLID, 1)],
    );
    let model = SaltPrecipitation::from_sources(
        fixture.aqueous(false), reaction.clone(), id(METAL), id(DRY_SOLID), FAMILY,
        &fixture.sources(&reaction),
    ).unwrap();
    let result = solve(&model, 1.1e-4, 1.2e-4, &[]);
    near(result.dissolved_cation.value(), 1e-5, 4e-9, 0.0);
    near(result.dissolved_family.value(), 2e-5, 4e-9, 0.0);
    assert_eq!(result.model().standard().terms().len(), 3);
    closed(&result);
}

#[test]
fn invalid_topology_and_nonidentical_sources_cannot_enter_the_shared_model() {
    let fixture = Fixture::new(5e-11);
    let phase = MineralReaction::ScoroditePrecipitation.build().unwrap();
    for family in ["missing", "buffer"] {
        assert!(SaltPrecipitation::from_sources(
            fixture.aqueous(true), phase.clone(), id(METAL), id(SOLID), family,
            &fixture.sources(&phase),
        ).is_err());
    }
    for factor in [-1, 2] {
        let reaction = fixture.reaction(
            "invalid-salt-direction-or-scale",
            &[(METAL, -factor), (FORMS[3], -factor), (WATER, -2 * factor), (SOLID, factor)],
        );
        assert!(SaltPrecipitation::from_sources(
            fixture.aqueous(false), reaction.clone(), id(METAL), id(SOLID), FAMILY,
            &fixture.sources(&reaction),
        ).is_err());
    }
    let mut foreign = fixture.clone();
    let water = foreign.index(WATER);
    foreign.records[water].pack = pack(
        "foreign-water-source", 0.0, &foreign.records[water].point,
    );
    assert!(SaltPrecipitation::from_sources(
        fixture.aqueous(false), phase.clone(), id(METAL), id(SOLID), FAMILY,
        &foreign.sources(&phase),
    ).is_err());

    let solid = foreign.index(SOLID);
    foreign = fixture.clone();
    foreign.records[solid].point = point("source-phase-aqueous", 0);
    foreign.records[solid].pack = pack(
        SOLID, foreign.records[solid].gibbs, &foreign.records[solid].point,
    );
    assert!(SaltPrecipitation::from_sources(
        fixture.aqueous(false), phase.clone(), id(METAL), id(SOLID), FAMILY,
        &foreign.sources(&phase),
    ).is_err());

    foreign = fixture.clone();
    let anion = foreign.index(FORMS[3]);
    foreign.records[anion].atoms[1] += 1;
    foreign.records[solid].atoms[1] += 1;
    let redefined = foreign.reaction(
        "redefined-salt-anion",
        &[(METAL, -1), (FORMS[3], -1), (WATER, -2), (SOLID, 1)],
    );
    assert!(matches!(
        SaltPrecipitation::from_sources(
            fixture.aqueous(false), redefined.clone(), id(METAL), id(SOLID), FAMILY,
            &fixture.sources(&redefined),
        ),
        Err(AcidBaseError::Network(_))
    ));
}

#[test]
fn invalid_totals_duplicates_controls_and_search_floor_refuse() {
    let model = Fixture::new(5e-11).model(false);
    let families = [total(FAMILY, 1.2e-4)];
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(model.speciate(
            Amount::new(invalid), &families, &[], Amount::new(1.0),
            options(), never_cancelled,
        ).is_err());
        assert!(model.speciate(
            Amount::new(1.1e-4), &[total(FAMILY, invalid)], &[], Amount::new(1.0),
            options(), never_cancelled,
        ).is_err());
    }
    for families in [
        vec![],
        vec![total(FAMILY, 1.2e-4), total(FAMILY, 1.2e-4)],
        vec![total(FAMILY, 1.2e-4), total("unknown", 0.0)],
    ] {
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            options(), never_cancelled,
        ).is_err());
    }
    for name in [METAL, SOLID, FORMS[3]] {
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[ion(name, 1, 0.0)],
            Amount::new(1.0), options(), never_cancelled,
        ).is_err());
    }
    for invalid in [0.0, f64::NAN, f64::INFINITY] {
        let mut settings = options();
        settings.log_saturation_tolerance = Dimensionless::new(invalid);
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            settings, never_cancelled,
        ).is_err());
        let mut settings = options();
        settings.log_activity_tolerance = Dimensionless::new(invalid);
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            settings, never_cancelled,
        ).is_err());
        let mut settings = options();
        settings.minimum_dissolved_amount = Amount::new(invalid);
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            settings, never_cancelled,
        ).is_err());
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(invalid),
            options(), never_cancelled,
        ).is_err());
    }
    let mut settings = options();
    settings.minimum_dissolved_amount = Amount::new(5e-5);
    assert!(matches!(
        model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            settings, never_cancelled,
        ),
        Err(AcidBaseError::OutsideBracket)
    ));
    for budget in [0, 1, 257] {
        let mut settings = options();
        settings.max_iterations = budget;
        assert!(model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0),
            settings, never_cancelled,
        ).is_err());
    }
}

#[test]
fn cancellation_and_unresolved_inner_ph_never_publish_a_partition() {
    let model = Fixture::new(5e-11).model(false);
    let families = [total(FAMILY, 1.2e-4)];
    assert!(matches!(
        model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0), options(), || true,
        ),
        Err(AcidBaseError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        model.speciate(
            Amount::new(1.1e-4), &families, &[], Amount::new(1.0), options(),
            || { calls += 1; calls == 20 },
        ),
        Err(AcidBaseError::Cancelled)
    ));
    assert_eq!(calls, 20);

    // Choose an independent non-grid pH, then construct K and the spectator
    // charge directly. A deliberately loose inner solve cannot certify it at
    // a much tighter outer saturation tolerance.
    let h = 10_f64.powf(-7.37);
    let weight = KA / h;
    let weights = [1.0, weight, weight * weight, weight * weight * weight];
    let denominator: f64 = weights.iter().sum();
    let mean_charge = -(weights[1] + 2.0 * weights[2] + 3.0 * weights[3]) / denominator;
    let metal = 1e-5;
    let family = 2e-5;
    let kdiss = metal * family * weights[3] / denominator;
    let ions = counterions(h - KW / h + 3.0 * metal + mean_charge * family);
    let model = Fixture::new(kdiss).model(false);
    let mut settings = options();
    settings.aqueous.ph_tolerance = Dimensionless::new(1e-2);
    settings.aqueous.absolute_charge_tolerance = Amount::new(1e-6);
    settings.aqueous.relative_charge_tolerance = Dimensionless::new(0.0);
    settings.log_saturation_tolerance = Dimensionless::new(1e-14);
    assert!(model.speciate(
        Amount::new(metal + 1e-4), &[total(FAMILY, family + 1e-4)], &ions,
        Amount::new(1.0), settings, never_cancelled,
    ).is_err());
}
