//! G0/G3/G4: independent manufactured phase states and finite-inventory refusals.
use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::minerals::MineralReaction,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const REFERENCE: &str = "source-elemental-reference-authored-copper-sulfide-test";
const SOLID: &str = "Cu2S(s-nbs-alpha)";
const CU: &str = "Cu[+](aq)";
const COMPLEX: &str = "Cu(CN)2[-](aq)";
const CN: &str = "CN[-](aq)";
const H: &str = "H+(aq)";
const HS: &str = "HS[-](aq)";
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

#[derive(Clone)]
struct Record {
    id: SpeciesId,
    pack: NormalizedPack,
    point: QueryPoint,
}
impl Record {
    fn new(name: &str, charge: i32, gibbs: f64) -> Self {
        let phase = if name == SOLID {
            "source-phase-solid"
        } else if name == "H2O(l)" {
            "source-phase-liquid"
        } else {
            "source-phase-aqueous"
        };
        let mut point = QueryPoint::new()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                298.15,
            )
            .unwrap()
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                100_000.0,
            )
            .unwrap();
        let mut flags = vec![
            REFERENCE,
            "source-standard-state",
            phase,
            if phase == "source-phase-aqueous" {
                "source-aqueous-unit-activity-infinite-dilution"
            } else {
                "source-pure-species"
            },
        ];
        if charge != 0 {
            flags.extend([
                "source-individual-ion",
                "source-ionic-reference-hydrogen-zero",
            ]);
        }
        for flag in flags {
            point = point
                .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
                .unwrap();
        }
        let mut validity = ValidityDomain::unconstrained();
        for (axis, value) in point.axes() {
            validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
        }
        let mut claim = PropertyClaim {
            key: PropertyKey::new(
                "standard-molar-formation-gibbs-energy",
                MolarEnergyQuantityV1::DIMS,
            ),
            value: PropertyValue::Scalar {
                value: if gibbs == 0.0 { 0.0 } else { gibbs },
                dims: MolarEnergyQuantityV1::DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: format!("Authored Cu2S equilibrium fixture {name}"),
                license: "project fixture".into(),
                artifact: None,
            },
        };
        let hash = claim.content_hash();
        let mut claims = ClaimSet::new();
        claim.observations.push(
            claims
                .register_observation(ObservationDataset {
                    specimen: name.into(),
                    method: "manufactured standard Gibbs energy".into(),
                    artifact: hash,
                    caveats: "not measured thermochemical data".into(),
                    provenance: claim.provenance.clone(),
                })
                .unwrap(),
        );
        claims.insert_claim(claim).unwrap();
        Self {
            id: id(name),
            point,
            pack: NormalizedPack::new(
                name,
                "fixture-v1",
                hash,
                "project fixture",
                claims,
                vec![],
                vec![],
            )
            .unwrap(),
        }
    }
}
fn refs<'a>(records: &'a [Record], reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .map(|t| {
            let r = records.iter().find(|r| r.id == t.species).unwrap();
            FormationGibbsSource {
                species: &r.id,
                pack: &r.pack,
                at: &r.point,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn association() -> ConservedReaction {
    let species = vec![id(CU), id(CN), id(COMPLEX)];
    let elements = ElementalMatrix::new(
        vec![
            ElementId::new("Cu").unwrap(),
            ElementId::new("C").unwrap(),
            ElementId::new("N").unwrap(),
        ],
        species.clone(),
        vec![vec![1, 0, 1], vec![0, 1, 2], vec![0, 1, 2]],
    )
    .unwrap();
    let reaction = ReactionId::new("authored-copper-dicyanide").unwrap();
    let matrix = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-1], vec![-2], vec![1]],
    )
    .unwrap();
    ConservedReaction::new(
        &elements,
        &matrix,
        &ChargeVector::new(species, vec![1, -1, -1]).unwrap(),
        reaction,
    )
    .unwrap()
}
fn model(records: &[Record]) -> CopperSulfidePrecipitation {
    let water = MineralReaction::WaterDissociation.build().unwrap();
    let water = WaterAutoprotolysis::from_sources(
        water.clone(),
        id("H2O(l)"),
        id(H),
        id("OH[-](aq)"),
        &refs(records, &water),
    )
    .unwrap();
    let families = [
        (
            "cyanide",
            MineralReaction::HydrogenCyanideDissociation,
            "HCN(aq)",
            CN,
        ),
        (
            "sulfide",
            MineralReaction::HydrogenSulfideDissociation,
            "H2S(aq)",
            HS,
        ),
    ]
    .into_iter()
    .map(|(name, recipe, acid, base)| {
        let r = recipe.build().unwrap();
        AcidFamily::new(
            name,
            vec![
                AcidDissociation::from_sources(
                    r.clone(),
                    id(acid),
                    id(base),
                    id(H),
                    &refs(records, &r),
                )
                .unwrap(),
            ],
        )
        .unwrap()
    })
    .collect();
    let aqueous = IdealAcidBaseSystem::new(water, families).unwrap();
    let r = association();
    let copper = MetalLigandFamily::new(
        "copper",
        vec![
            LigandAssociation::from_sources(
                r.clone(),
                id(CU),
                id(CN),
                id(COMPLEX),
                &refs(records, &r),
            )
            .unwrap(),
        ],
    )
    .unwrap();
    let ligand = IdealLigandSystem::new(aqueous, "cyanide", vec![copper]).unwrap();
    let phase = MineralReaction::CopperSulfideHydrosulfidePrecipitation
        .build()
        .unwrap();
    CopperSulfidePrecipitation::from_sources(
        ligand,
        phase.clone(),
        "copper",
        "sulfide",
        id(SOLID),
        &refs(records, &phase),
    )
    .unwrap()
}

struct Known {
    records: Vec<Record>,
    copper: f64,
    sulfide: f64,
    cyanide: f64,
    sodium: f64,
    free_copper: f64,
    complex: f64,
    hs: f64,
    solid: f64,
}
impl Known {
    fn new(free_copper: f64, hs: f64, solid: f64, phase_k_factor: f64) -> Self {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
        let h = 1e-6;
        let cn = 1e-4;
        let complex = 1e6 * cn * cn * free_copper;
        let hcn = cn * h / 1e-6;
        let h2s = hs * h / 1e-7;
        let k = phase_k_factor * h / (free_copper * free_copper * hs);
        let records = vec![
            Record::new(H, 1, 0.0),
            Record::new("OH[-](aq)", -1, -rt * 1e-14_f64.ln()),
            Record::new("H2O(l)", 0, 0.0),
            Record::new(CN, -1, 0.0),
            Record::new("HCN(aq)", 0, rt * 1e-6_f64.ln()),
            Record::new(HS, -1, 0.0),
            Record::new("H2S(aq)", 0, rt * 1e-7_f64.ln()),
            Record::new(CU, 1, 0.0),
            Record::new(COMPLEX, -1, -rt * 1e6_f64.ln()),
            Record::new(SOLID, 0, -rt * k.ln()),
        ];
        Self {
            records,
            copper: free_copper + complex + 2.0 * solid,
            sulfide: hs + h2s + solid,
            cyanide: cn + hcn + 2.0 * complex,
            sodium: cn + hs + complex + 1e-14 / h - h - free_copper,
            free_copper,
            complex,
            hs,
            solid,
        }
    }
    fn run(&self, options: CopperSulfideOptions) -> Result<CopperSulfideSpeciation> {
        model(&self.records).speciate(
            Amount::new(self.cyanide),
            &[("copper".into(), Amount::new(self.copper))],
            &[("sulfide".into(), Amount::new(self.sulfide))],
            &ions(self.sodium),
            Amount::new(1.0),
            options,
            || false,
        )
    }
}
fn ions(net: f64) -> Vec<SpectatorIon> {
    vec![SpectatorIon {
        species: id(if net >= 0.0 { "Na[+](aq)" } else { "Cl[-](aq)" }),
        charge: if net >= 0.0 { 1 } else { -1 },
        amount: Amount::new(net.abs()),
    }]
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "{actual:e} != {expected:e}"
    );
}
fn check_known(f: &Known, a: &CopperSulfideSpeciation) {
    near(a.ligand.aqueous.ph.value(), 6.0, 2e-8);
    near(a.ligand.free_ligand.value(), 1e-4, 2e-7);
    near(
        a.ligand.metals[0].members[0].amount.value(),
        f.free_copper,
        2e-7,
    );
    near(
        a.ligand.metals[0].members[1].amount.value(),
        f.complex,
        3e-7,
    );
    near(a.dissolved_sulfide.value(), 11.0 * f.hs, 3e-7);
    if f.solid > 0.0 {
        near(a.solid_amount.value(), f.solid, 2e-7);
    }
    assert!(a.copper_balance_residual.value().abs() <= 1e-14 * f.copper);
    assert!(a.sulfide_balance_residual.value().abs() <= 1e-14 * f.sulfide);
    let q = a.ligand.aqueous.proton_amount.value() - a.ligand.aqueous.hydroxide_amount.value()
        + a.ligand
            .aqueous
            .spectators
            .iter()
            .map(|s| f64::from(s.charge) * s.amount.value())
            .sum::<f64>()
        - a.ligand
            .aqueous
            .families
            .iter()
            .map(|f| f.forms[1].amount.value())
            .sum::<f64>();
    assert!(q.abs() < 2e-14, "independent charge {q:e}");
}

#[test]
fn independent_known_states_with_either_limiting_component() {
    for hs in [1e-4, 1e-7] {
        let f = Known::new(1e-4, hs, 0.001, 1.0);
        let a = f.run(CopperSulfideOptions::default()).unwrap();
        assert_eq!(a.status, PrecipitationStatus::SolidPresent);
        assert!(
            a.log_saturation_ratio.unwrap().value().abs()
                + a.log_saturation_numerical_spread.unwrap().value()
                <= 1e-7
        );
        check_known(&f, &a);
    }
}

#[test]
fn absent_solid_keeps_all_components_and_resolves_undersaturation() {
    let f = Known::new(1e-4, 1e-4, 0.0, 0.1);
    let a = f.run(CopperSulfideOptions::default()).unwrap();
    assert_eq!(a.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(a.solid_amount.value(), 0.0);
    assert!(a.log_saturation_ratio.unwrap().value() < -2.0);
    check_known(&f, &a);
}

#[test]
fn near_complete_precipitation_keeps_positive_independent_copper_trace() {
    let f = Known::new(1e-25, 1e-4, 1.0, 1.0);
    let a = f.run(CopperSulfideOptions::default()).unwrap();
    assert_eq!(a.copper_total.value() - 2.0 * a.solid_amount.value(), 0.0);
    assert!(a.dissolved_copper.value() > 1e-26);
    check_known(&f, &a);
}

#[test]
fn zero_components_have_no_solid_or_defined_saturation() {
    let f = Known::new(1e-4, 1e-4, 0.001, 1.0);
    let m = model(&f.records);
    for (cu, s, net) in [(0.0, 0.0, 0.0), (1e-4, 0.0, -1e-4), (0.0, 1e-4, 1e-4)] {
        let a = m
            .speciate(
                Amount::new(0.0),
                &[("copper".into(), Amount::new(cu))],
                &[("sulfide".into(), Amount::new(s))],
                &ions(net),
                Amount::new(1.0),
                CopperSulfideOptions::default(),
                || false,
            )
            .unwrap();
        assert_eq!(a.status, PrecipitationStatus::EmptyComponent);
        assert_eq!(a.solid_amount.value(), 0.0);
        assert!(a.log_saturation_ratio.is_none());
    }
}

#[test]
fn solvent_basis_scaling_preserves_activities_and_phase_fractions() {
    let f = Known::new(1e-4, 1e-4, 0.001, 1.0);
    let a = f.run(CopperSulfideOptions::default()).unwrap();
    let scale = 100.0;
    let b = model(&f.records)
        .speciate(
            Amount::new(scale * f.cyanide),
            &[("copper".into(), Amount::new(scale * f.copper))],
            &[("sulfide".into(), Amount::new(scale * f.sulfide))],
            &ions(scale * f.sodium),
            Amount::new(scale),
            CopperSulfideOptions::default(),
            || false,
        )
        .unwrap();
    near(
        b.ligand.aqueous.ph.value(),
        a.ligand.aqueous.ph.value(),
        2e-8,
    );
    near(b.solid_amount.value(), scale * a.solid_amount.value(), 2e-7);
    near(
        b.dissolved_copper.value(),
        scale * a.dissolved_copper.value(),
        2e-7,
    );
    near(
        b.ligand.free_ligand.value(),
        scale * a.ligand.free_ligand.value(),
        2e-7,
    );
}

#[test]
fn shared_source_changes_and_wrong_phase_topology_refuse() {
    let f = Known::new(1e-4, 1e-4, 0.001, 1.0);
    let m = model(&f.records);
    let mut changed = f.records.clone();
    *changed.iter_mut().find(|r| r.id == id(CU)).unwrap() = Record::new(CU, 1, 1000.0);
    assert!(
        CopperSulfidePrecipitation::from_sources(
            m.ligand_model().clone(),
            m.reaction().clone(),
            "copper",
            "sulfide",
            id(SOLID),
            &refs(&changed, m.reaction())
        )
        .is_err()
    );
    let wrong = MineralReaction::CopperCyanideHydrosulfideRecovery
        .build()
        .unwrap();
    assert!(
        CopperSulfidePrecipitation::from_sources(
            m.ligand_model().clone(),
            wrong,
            "copper",
            "sulfide",
            id(SOLID),
            &[]
        )
        .is_err()
    );
}

#[test]
fn missing_duplicate_solid_and_outside_ph_inputs_refuse() {
    let f = Known::new(1e-4, 1e-4, 0.001, 1.0);
    let m = model(&f.records);
    let metals = vec![("copper".into(), Amount::new(f.copper))];
    let acids = vec![("sulfide".into(), Amount::new(f.sulfide))];
    assert!(
        m.speciate(
            Amount::new(f.cyanide),
            &metals,
            &[],
            &ions(f.sodium),
            Amount::new(1.0),
            CopperSulfideOptions::default(),
            || false
        )
        .is_err()
    );
    assert!(
        m.speciate(
            Amount::new(f.cyanide),
            &metals,
            &[acids[0].clone(), acids[0].clone()],
            &ions(f.sodium),
            Amount::new(1.0),
            CopperSulfideOptions::default(),
            || false
        )
        .is_err()
    );
    let duplicate = [SpectatorIon {
        species: id(SOLID),
        charge: 0,
        amount: Amount::new(1.0),
    }];
    assert!(
        m.speciate(
            Amount::new(f.cyanide),
            &metals,
            &acids,
            &duplicate,
            Amount::new(1.0),
            CopperSulfideOptions::default(),
            || false
        )
        .is_err()
    );
    let zero_cu = [("copper".into(), Amount::new(0.0))];
    let zero_s = [("sulfide".into(), Amount::new(0.0))];
    for net in [0.01, -0.1] {
        assert!(
            m.speciate(
                Amount::new(0.0),
                &zero_cu,
                &zero_s,
                &ions(net),
                Amount::new(1.0),
                CopperSulfideOptions::default(),
                || false
            )
            .is_err()
        );
    }
}

#[test]
fn explicit_trace_floor_work_limits_and_cancellation_refuse() {
    let f = Known::new(1e-4, 1e-4, 0.001, 1.0);
    assert_eq!(
        f.run(CopperSulfideOptions {
            max_iterations: 1,
            ..CopperSulfideOptions::default()
        })
        .unwrap_err(),
        LigandSpeciationError::Budget
    );
    assert_eq!(
        f.run(CopperSulfideOptions {
            max_ligand_solves: 1,
            ..CopperSulfideOptions::default()
        })
        .unwrap_err(),
        LigandSpeciationError::Budget
    );
    assert_eq!(
        f.run(CopperSulfideOptions {
            minimum_dissolved_amount: Amount::new(1e-3),
            ..CopperSulfideOptions::default()
        })
        .unwrap_err(),
        LigandSpeciationError::OutsideBracket
    );
    let m = model(&f.records);
    let mut calls = 0;
    let error = m
        .speciate(
            Amount::new(f.cyanide),
            &[("copper".into(), Amount::new(f.copper))],
            &[("sulfide".into(), Amount::new(f.sulfide))],
            &ions(f.sodium),
            Amount::new(1.0),
            CopperSulfideOptions::default(),
            || {
                calls += 1;
                calls > 100
            },
        )
        .unwrap_err();
    assert_eq!(error, LigandSpeciationError::Cancelled);
    assert_eq!(calls, 101);
}
