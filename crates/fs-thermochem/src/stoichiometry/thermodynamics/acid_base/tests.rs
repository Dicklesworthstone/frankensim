//! Independent ideal-solution fixtures; none of these authored Gibbs values
//! represents measured acid constants, a calibrated solution, or plant recovery.

use super::*;
use crate::{ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dimensionless, Dims, QuantitySpec};

const REFERENCE: &str = "source-elemental-reference-authored-acid-base-test";
const PROPERTY: &str = "standard-molar-formation-gibbs-energy";
const TEMPERATURE: f64 = 298.15;
const H: &str = "H+(aq)";
const OH: &str = "OH[-](aq)";
const WATER: &str = "H2O(l)";
const HA: &str = "authored-HA(aq)";
const A: &str = "authored-A[-](aq)";
const H2L: &str = "authored-H2L(aq)";
const HL: &str = "authored-HL[-](aq)";
const L: &str = "authored-L[2-](aq)";

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn options() -> SpeciationOptions {
    SpeciationOptions {
        ph_bracket: [Dimensionless::new(-2.0), Dimensionless::new(16.0)],
        ph_tolerance: Dimensionless::new(1e-12),
        absolute_charge_tolerance: Amount::new(1e-18),
        relative_charge_tolerance: Dimensionless::new(1e-12),
        max_iterations: 256,
    }
}

fn point(phase: &str, ionic: bool, reference: &str) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![reference, "source-standard-state", phase];
    if phase == "source-phase-aqueous" {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
        if ionic {
            flags.extend([
                "source-ionic-reference-hydrogen-zero",
                "source-individual-ion",
            ]);
        }
    } else {
        flags.push("source-pure-species");
    }
    for flag in flags {
        at = at
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    at
}

fn pack(name: &str, gibbs: f64, at: &QueryPoint, omit: Option<&str>) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in at.axes() {
        if omit != Some(axis.as_str()) {
            validity = validity.with_quantity(axis, at.axis_quantities()[axis], *value, *value);
        }
    }
    let dims = crate::MolarEnergyQuantityV1::DIMS;
    let mut claim = PropertyClaim {
        key: PropertyKey::new(PROPERTY, dims),
        value: PropertyValue::Scalar { value: gibbs, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("authored acid-base arithmetic fixture: {name}"),
            license: "project fixture".into(),
            artifact: None,
        },
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("synthetic exact-composition state {name}"),
            method: "independent equilibrium known-answer construction".into(),
            artifact: hash,
            caveats: "not experimental thermochemistry or process calibration".into(),
            provenance: claim.provenance.clone(),
        })
        .unwrap();
    claim.observations.push(observation);
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        name,
        "test",
        hash,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

struct Definition {
    name: &'static str,
    // H, C, N, O. The family labels are deliberately synthetic, but their
    // integer compositions and charges are checked by the real chemistry owner.
    atoms: [i128; 4],
    charge: i128,
    gibbs: f64,
}

struct Fixture {
    definitions: Vec<Definition>,
    points: Vec<QueryPoint>,
    packs: Vec<NormalizedPack>,
    reference: &'static str,
}

impl Fixture {
    fn new(kw: f64, ka: f64, poly: [f64; 2]) -> Self {
        Self::with_reference(kw, ka, poly, REFERENCE)
    }

    fn with_reference(kw: f64, ka: f64, poly: [f64; 2], reference: &'static str) -> Self {
        let rt = crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE;
        let definitions = vec![
            Definition {
                name: H,
                atoms: [1, 0, 0, 0],
                charge: 1,
                gibbs: 0.0,
            },
            Definition {
                name: OH,
                atoms: [1, 0, 0, 1],
                charge: -1,
                gibbs: -rt * kw.ln(),
            },
            Definition {
                name: WATER,
                atoms: [2, 0, 0, 1],
                charge: 0,
                gibbs: 0.0,
            },
            Definition {
                name: HA,
                atoms: [1, 1, 1, 0],
                charge: 0,
                gibbs: 0.0,
            },
            Definition {
                name: A,
                atoms: [0, 1, 1, 0],
                charge: -1,
                gibbs: -rt * ka.ln(),
            },
            Definition {
                name: H2L,
                atoms: [2, 1, 0, 3],
                charge: 0,
                gibbs: 0.0,
            },
            Definition {
                name: HL,
                atoms: [1, 1, 0, 3],
                charge: -1,
                gibbs: -rt * poly[0].ln(),
            },
            Definition {
                name: L,
                atoms: [0, 1, 0, 3],
                charge: -2,
                gibbs: -rt * (poly[0].ln() + poly[1].ln()),
            },
        ];
        let points: Vec<_> = definitions
            .iter()
            .map(|d| {
                point(
                    if d.name == WATER {
                        "source-phase-liquid"
                    } else {
                        "source-phase-aqueous"
                    },
                    d.charge != 0,
                    reference,
                )
            })
            .collect();
        let packs = definitions
            .iter()
            .zip(&points)
            .map(|(d, at)| pack(d.name, d.gibbs, at, None))
            .collect();
        Self {
            definitions,
            points,
            packs,
            reference,
        }
    }

    fn index(&self, name: &str) -> usize {
        self.definitions
            .iter()
            .position(|d| d.name == name)
            .unwrap()
    }

    fn reaction(&self, name: &str, terms: &[(&str, i128)]) -> ConservedReaction {
        let definitions: Vec<_> = terms
            .iter()
            .map(|(name, _)| &self.definitions[self.index(name)])
            .collect();
        let species: Vec<_> = definitions.iter().map(|d| id(d.name)).collect();
        let reaction = ReactionId::new(name).unwrap();
        let elemental = ElementalMatrix::new(
            ["H", "C", "N", "O"]
                .into_iter()
                .map(|s| ElementId::new(s).unwrap())
                .collect(),
            species.clone(),
            (0..4)
                .map(|row| definitions.iter().map(|d| d.atoms[row]).collect())
                .collect(),
        )
        .unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            species.clone(),
            vec![reaction.clone()],
            terms.iter().map(|(_, n)| vec![*n]).collect(),
        )
        .unwrap();
        let charge =
            ChargeVector::new(species, definitions.iter().map(|d| d.charge).collect()).unwrap();
        ConservedReaction::new(&elemental, &stoichiometric, &charge, reaction).unwrap()
    }

    fn sources<'a>(&'a self, reaction: &'a ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let index = self.index(term.species.as_str());
                FormationGibbsSource {
                    species: &term.species,
                    pack: &self.packs[index],
                    at: &self.points[index],
                    reference_axis: self.reference,
                }
            })
            .collect()
    }

    fn step(&self, name: &str, acid: &str, base: &str) -> Result<AcidDissociation> {
        let reaction = self.reaction(name, &[(acid, -1), (H, 1), (base, 1)]);
        AcidDissociation::from_sources(
            reaction.clone(),
            id(acid),
            id(base),
            id(H),
            &self.sources(&reaction),
        )
    }

    fn water(&self) -> Result<WaterAutoprotolysis> {
        let reaction = self.reaction("authored-water-ionization", &[(WATER, -1), (H, 1), (OH, 1)]);
        WaterAutoprotolysis::from_sources(
            reaction.clone(),
            id(WATER),
            id(H),
            id(OH),
            &self.sources(&reaction),
        )
    }

    fn mono(&self) -> AcidFamily {
        AcidFamily::new("monoprotic", vec![self.step("authored-ka", HA, A).unwrap()]).unwrap()
    }

    fn poly(&self) -> AcidFamily {
        AcidFamily::new(
            "diprotic",
            vec![
                self.step("authored-k1", H2L, HL).unwrap(),
                self.step("authored-k2", HL, L).unwrap(),
            ],
        )
        .unwrap()
    }

    fn system(&self, families: Vec<AcidFamily>) -> IdealAcidBaseSystem {
        IdealAcidBaseSystem::new(self.water().unwrap(), families).unwrap()
    }
}

fn ion(name: &str, charge: i32, amount: f64) -> SpectatorIon {
    SpectatorIon {
        species: id(name),
        charge,
        amount: Amount::new(amount),
    }
}

fn total(name: &str, amount: f64) -> (String, Amount) {
    (name.into(), Amount::new(amount))
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn family<'a>(report: &'a AcidBaseSpeciation, name: &str) -> &'a FamilySpeciation {
    report.families.iter().find(|f| f.id == name).unwrap()
}

fn form<'a>(
    report: &'a AcidBaseSpeciation,
    family_name: &str,
    species: &str,
) -> &'a AcidBaseSpecies {
    family(report, family_name)
        .forms
        .iter()
        .find(|f| f.species == id(species))
        .unwrap()
}

fn solve(
    system: &IdealAcidBaseSystem,
    totals: &[(String, Amount)],
    ions: &[SpectatorIon],
    scale: f64,
) -> AcidBaseSpeciation {
    system
        .speciate(totals, ions, Amount::new(scale), options(), || false)
        .unwrap()
}

#[test]
fn pure_water_and_strong_spectator_ions_match_independent_quadratics() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![]);
    let scale = 2.5;
    let neutral = solve(&system, &[], &[], scale);
    near(neutral.ph.value(), 7.0, 0.0, 2e-11);
    near(neutral.proton_amount.value(), scale * 1e-7, 2e-11, 0.0);
    near(neutral.hydroxide_amount.value(), scale * 1e-7, 2e-11, 0.0);
    assert!(neutral.families.is_empty());

    let q = 0.002;
    let product = scale * scale * 1e-14;
    let larger = (q + (q * q + 4.0 * product).sqrt()) / 2.0;
    let smaller = product / larger;
    for (charge, expected_h, expected_oh) in [(-1, larger, smaller), (1, smaller, larger)] {
        let report = solve(&system, &[], &[ion("declared-spectator", charge, q)], scale);
        near(report.proton_amount.value(), expected_h, 2e-10, 0.0);
        near(report.hydroxide_amount.value(), expected_oh, 2e-10, 0.0);
        near(
            report.proton_activity.value(),
            expected_h / scale,
            2e-10,
            0.0,
        );
        near(report.ph.value(), -(expected_h / scale).log10(), 0.0, 2e-10);
        near(report.charge_balance_residual.value(), 0.0, 0.0, 3e-14);
        assert_eq!(
            report.model().water().reaction().id(),
            system.water().reaction().id()
        );
    }
}

#[test]
fn mixed_families_share_one_proton_activity_and_retain_each_analytical_total() {
    // At h=1e-7, HA:A = 1:1 and H2L:HL:L = 1:2:1. The
    // spectator charge .02/2 + .03 balances these families exactly; water
    // contributes equal H and OH. These ratios precede any root algorithm.
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![fixture.mono(), fixture.poly()]);
    let totals = [total("monoprotic", 0.02), total("diprotic", 0.03)];
    let ions = [ion("authored-Na+", 1, 0.01), ion("authored-Ca2+", 2, 0.015)];
    let report = solve(&system, &totals, &ions, 1.0);
    near(report.ph.value(), 7.0, 0.0, 2e-10);
    for (family_name, species, fraction, amount) in [
        ("monoprotic", HA, 0.5, 0.01),
        ("monoprotic", A, 0.5, 0.01),
        ("diprotic", H2L, 0.25, 0.0075),
        ("diprotic", HL, 0.5, 0.015),
        ("diprotic", L, 0.25, 0.0075),
    ] {
        let actual = form(&report, family_name, species);
        near(actual.fraction.value(), fraction, 2e-10, 0.0);
        near(actual.amount.value(), amount, 2e-10, 0.0);
        near(actual.activity.value(), amount, 2e-10, 0.0);
    }
    for f in &report.families {
        near(
            f.forms.iter().map(|s| s.amount.value()).sum(),
            f.total.value(),
            5e-15,
            1e-18,
        );
    }
    near(report.charge_balance_residual.value(), 0.0, 0.0, 1e-13);
    for step in report.model().families().iter().flat_map(|f| f.steps()) {
        assert_eq!(step.standard().terms().len(), 3);
        assert!(
            step.standard()
                .terms()
                .iter()
                .all(|t| !t.answer.receipt.source_hashes.is_empty())
        );
    }
}

#[test]
fn input_permutation_and_extensive_rescaling_preserve_the_same_solution() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![fixture.mono(), fixture.poly()]);
    let reversed = fixture.system(vec![fixture.poly(), fixture.mono()]);
    let totals = [total("monoprotic", 0.02), total("diprotic", 0.03)];
    let ions = [ion("authored-Na+", 1, 0.01), ion("authored-Ca2+", 2, 0.015)];
    let baseline = solve(&system, &totals, &ions, 1.0);
    let permuted = solve(
        &reversed,
        &[totals[1].clone(), totals[0].clone()],
        &[ions[1].clone(), ions[0].clone()],
        1.0,
    );
    assert_eq!(baseline.ph.value().to_bits(), permuted.ph.value().to_bits());
    let factor = 7.5;
    let scaled = solve(
        &system,
        &[
            total("monoprotic", 0.02 * factor),
            total("diprotic", 0.03 * factor),
        ],
        &[
            ion("authored-Na+", 1, 0.01 * factor),
            ion("authored-Ca2+", 2, 0.015 * factor),
        ],
        factor,
    );
    near(scaled.ph.value(), baseline.ph.value(), 0.0, 2e-10);
    for f in &baseline.families {
        for species in &f.forms {
            let other = form(&scaled, &f.id, species.species.as_str());
            let permuted_species = form(&permuted, &f.id, species.species.as_str());
            assert_eq!(
                species.amount.value().to_bits(),
                permuted_species.amount.value().to_bits()
            );
            near(other.activity.value(), species.activity.value(), 3e-10, 0.0);
            near(
                other.amount.value(),
                factor * species.amount.value(),
                3e-10,
                0.0,
            );
        }
    }

    // Changing solvent basis alone IS dilution, unlike scaling the complete
    // inventory. A weak-base salt approaches neutral water on dilution.
    let mono = fixture.system(vec![fixture.mono()]);
    let concentrated = solve(
        &mono,
        &[total("monoprotic", 0.01)],
        &[ion("authored-Na+", 1, 0.01)],
        1.0,
    );
    let diluted = solve(
        &mono,
        &[total("monoprotic", 0.01)],
        &[ion("authored-Na+", 1, 0.01)],
        100.0,
    );
    assert!(concentrated.ph.value() - diluted.ph.value() > 0.8);
    assert!(
        form(&diluted, "monoprotic", A).fraction.value()
            < form(&concentrated, "monoprotic", A).fraction.value()
    );
    near(family(&diluted, "monoprotic").total.value(), 0.01, 0.0, 0.0);
}

#[test]
fn counterion_is_explicit_and_tiny_free_base_is_not_replaced_by_total_or_zero() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![fixture.mono()]);
    let acid = solve(&system, &[total("monoprotic", 0.01)], &[], 1.0);
    let salt = solve(
        &system,
        &[total("monoprotic", 0.01)],
        &[ion("authored-Na+", 1, 0.01)],
        1.0,
    );
    assert!(salt.ph.value() - acid.ph.value() > 4.0);
    assert!(form(&acid, "monoprotic", A).fraction.value() < 0.01);
    assert!(form(&salt, "monoprotic", A).fraction.value() > 0.99);
    let absent = solve(&system, &[total("monoprotic", 0.0)], &[], 1.0);
    near(absent.ph.value(), 7.0, 0.0, 2e-10);
    assert!(
        family(&absent, "monoprotic")
            .forms
            .iter()
            .all(|s| s.amount.value() == 0.0)
    );

    let trace_fixture = Fixture::new(1e-14, 1e-100, [2e-7, 5e-8]);
    let trace_system = trace_fixture.system(vec![trace_fixture.mono()]);
    let trace = solve(
        &trace_system,
        &[total("monoprotic", 0.01)],
        &[ion("authored-Cl-", -1, 0.001)],
        1.0,
    );
    let expected = 0.01 * 1e-100 / (trace.proton_activity.value() + 1e-100);
    let actual = form(&trace, "monoprotic", A).amount.value();
    assert!(actual > 0.0);
    near(actual, expected, 3e-10, 0.0);
}

#[test]
fn original_phase_reference_and_shared_species_evidence_cannot_be_substituted() {
    let good = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let other = Fixture::with_reference(
        1e-14,
        1e-7,
        [2e-7, 5e-8],
        "source-elemental-reference-other-authored-test",
    );
    assert!(IdealAcidBaseSystem::new(good.water().unwrap(), vec![other.mono()]).is_err());

    let mut solid_acid = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let i = solid_acid.index(HA);
    solid_acid.points[i] = point("source-phase-solid", false, REFERENCE);
    solid_acid.packs[i] = pack(HA, 0.0, &solid_acid.points[i], None);
    assert!(solid_acid.step("authored-solid-acid", HA, A).is_err());

    let mut missing_ionic = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let i = missing_ionic.index(H);
    missing_ionic.packs[i] = pack(
        H,
        0.0,
        &missing_ionic.points[i],
        Some("source-ionic-reference-hydrogen-zero"),
    );
    assert!(missing_ionic.water().is_err());

    let mut altered_proton = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let i = altered_proton.index(H);
    altered_proton.packs[i] = pack(
        "different-authored-proton-source",
        0.0,
        &altered_proton.points[i],
        None,
    );
    // Both individual equations still have legitimate complete sources, but
    // one shared species cannot silently choose different original evidence.
    assert!(IdealAcidBaseSystem::new(good.water().unwrap(), vec![altered_proton.mono()]).is_err());

    let mut inconsistent_atoms = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let hl = inconsistent_atoms.index(HL);
    let l = inconsistent_atoms.index(L);
    inconsistent_atoms.definitions[hl].atoms[1] = 2;
    inconsistent_atoms.definitions[l].atoms[1] = 2;
    // Each step independently conserves atoms. The intermediate's conflicting
    // C count is exposed only when those original definitions are joined.
    assert!(
        AcidFamily::new(
            "inconsistent-intermediate",
            vec![
                good.step("authored-k1", H2L, HL).unwrap(),
                inconsistent_atoms.step("authored-k2", HL, L).unwrap(),
            ]
        )
        .is_err()
    );
}

#[test]
fn only_nonoverlapping_conserved_single_proton_chains_are_admitted() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    assert!(AcidFamily::new("empty", vec![]).is_err());
    assert!(
        AcidFamily::new(
            "disconnected",
            vec![
                fixture.step("authored-ka", HA, A).unwrap(),
                fixture.step("authored-k2", HL, L).unwrap()
            ]
        )
        .is_err()
    );
    assert!(
        AcidFamily::new(
            "reversed",
            vec![
                fixture.step("authored-k2", HL, L).unwrap(),
                fixture.step("authored-k1", H2L, HL).unwrap()
            ]
        )
        .is_err()
    );
    let skipped = fixture.reaction("authored-two-proton-step", &[(H2L, -1), (H, 2), (L, 1)]);
    assert!(
        AcidDissociation::from_sources(
            skipped.clone(),
            id(H2L),
            id(L),
            id(H),
            &fixture.sources(&skipped)
        )
        .is_err()
    );
    let scaled = fixture.reaction("authored-scaled-step", &[(HA, -2), (H, 2), (A, 2)]);
    assert!(
        AcidDissociation::from_sources(
            scaled.clone(),
            id(HA),
            id(A),
            id(H),
            &fixture.sources(&scaled)
        )
        .is_err()
    );
    let reversed_water =
        fixture.reaction("authored-water-reversed", &[(WATER, 1), (H, -1), (OH, -1)]);
    assert!(
        WaterAutoprotolysis::from_sources(
            reversed_water.clone(),
            id(WATER),
            id(H),
            id(OH),
            &fixture.sources(&reversed_water)
        )
        .is_err()
    );
    let alias = AcidFamily::new(
        "overlapping-alias",
        vec![fixture.step("authored-alias-ka", HA, A).unwrap()],
    )
    .unwrap();
    assert!(
        IdealAcidBaseSystem::new(fixture.water().unwrap(), vec![fixture.mono(), alias]).is_err()
    );
    assert!(
        IdealAcidBaseSystem::new(
            fixture.water().unwrap(),
            vec![fixture.mono(), fixture.mono()]
        )
        .is_err()
    );
}

#[test]
fn complete_totals_and_spectators_refuse_duplicate_or_double_counted_inventory() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![fixture.mono()]);
    for totals in [
        vec![],
        vec![total("unknown", 0.01)],
        vec![total("monoprotic", 0.01), total("monoprotic", 0.02)],
        vec![total("monoprotic", -1.0)],
        vec![total("monoprotic", f64::NAN)],
    ] {
        assert!(
            system
                .speciate(&totals, &[], Amount::new(1.0), options(), || false)
                .is_err()
        );
    }
    let totals = [total("monoprotic", 0.01)];
    for spectators in [
        vec![ion(H, 1, 0.01)],
        vec![ion(OH, -1, 0.01)],
        vec![ion(A, -1, 0.01)],
        vec![ion("authored-Na+", 1, 0.01), ion("authored-Na+", 1, 0.01)],
        vec![ion("authored-Na+", 1, -0.01)],
        vec![ion("authored-Na+", 1, f64::INFINITY)],
    ] {
        assert!(
            system
                .speciate(&totals, &spectators, Amount::new(1.0), options(), || false)
                .is_err()
        );
    }
}

#[test]
fn finite_search_bounds_work_and_cancellation_publish_no_partial_speciation() {
    let fixture = Fixture::new(1e-14, 1e-7, [2e-7, 5e-8]);
    let system = fixture.system(vec![fixture.mono()]);
    let totals = [total("monoprotic", 0.01)];
    let ions = [ion("authored-Na+", 1, 0.01)];
    let before = solve(&system, &totals, &ions, 1.0);
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            system
                .speciate(&totals, &ions, Amount::new(scale), options(), || false)
                .is_err()
        );
    }
    let mut invalid = options();
    invalid.ph_bracket = [Dimensionless::new(10.0), Dimensionless::new(3.0)];
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), invalid, || false)
            .is_err()
    );
    invalid = options();
    invalid.ph_bracket = [Dimensionless::new(3.0), Dimensionless::new(4.0)];
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), invalid, || false)
            .is_err()
    );
    invalid = options();
    invalid.ph_bracket[0] = Dimensionless::new(f64::NAN);
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), invalid, || false)
            .is_err()
    );
    invalid = options();
    invalid.ph_tolerance = Dimensionless::new(0.0);
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), invalid, || false)
            .is_err()
    );
    invalid = options();
    invalid.max_iterations = 1;
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), invalid, || false)
            .is_err()
    );
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), options(), || true)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        system
            .speciate(&totals, &ions, Amount::new(1.0), options(), || {
                polls += 1;
                polls > 2
            })
            .is_err()
    );
    assert!(polls > 2);
    let after = solve(&system, &totals, &ions, 1.0);
    assert_eq!(before.ph.value().to_bits(), after.ph.value().to_bits());
    assert_eq!(totals[0].1.value(), 0.01);
    assert_eq!(ions[0].amount.value(), 0.01);
}
