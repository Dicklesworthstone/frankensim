//! Test-only authored Gibbs values for exact atom/charge and equilibrium checks.
//! No value here represents measured cyanide, calcium or hydroxide solubility.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId, SpeciesId,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    stoichiometry::{
        ConservedReaction,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, WaterAutoprotolysis,
                precipitation::{HydrolysisPrecipitation, salt::SaltPrecipitation},
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-authored-electrolyte-material-test";
pub const H: &str = "H+(aq)";
pub const OH: &str = "OH[-](aq)";
pub const WATER: &str = "H2O(l)";
pub const HC: &str = "HCN(aq)";
pub const CN: &str = "CN[-](aq)";
pub const CA: &str = "Ca[2+](aq)";
pub const SALT: &str = "authored-Ca(CN)2(s)";
pub const HYDROXIDE: &str = "Ca(OH)2(s)";
pub const FAMILY: &str = "uncomplexed-cyanide";

pub fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

struct Record {
    name: &'static str,
    atoms: [i128; 5], // H, C, N, O, Ca
    charge: i128,
    point: QueryPoint,
    pack: NormalizedPack,
}
struct Fixture {
    records: Vec<Record>,
}
impl Fixture {
    fn new(temperature: f64, salt_kdiss: f64, hydroxide_kprec: f64) -> Self {
        let rt = R * temperature;
        let cn_gibbs = -rt * 1e-7_f64.ln();
        let records = [
            (H, [1, 0, 0, 0, 0], 1, 0.0, "source-phase-aqueous"),
            (
                OH,
                [1, 0, 0, 1, 0],
                -1,
                -rt * 1e-14_f64.ln(),
                "source-phase-aqueous",
            ),
            (WATER, [2, 0, 0, 1, 0], 0, 0.0, "source-phase-liquid"),
            (HC, [1, 1, 1, 0, 0], 0, 0.0, "source-phase-aqueous"),
            (CN, [0, 1, 1, 0, 0], -1, cn_gibbs, "source-phase-aqueous"),
            (CA, [0, 0, 0, 0, 1], 2, 0.0, "source-phase-aqueous"),
            (
                SALT,
                [0, 2, 2, 0, 1],
                0,
                2.0 * cn_gibbs + rt * salt_kdiss.ln(),
                "source-phase-solid",
            ),
            (
                HYDROXIDE,
                [2, 0, 0, 2, 1],
                0,
                -rt * hydroxide_kprec.ln(),
                "source-phase-solid",
            ),
        ]
        .into_iter()
        .map(|(name, atoms, charge, gibbs, phase)| {
            let point = point(temperature, phase, charge != 0);
            let pack = pack(name, gibbs, &point);
            Record {
                name,
                atoms,
                charge,
                point,
                pack,
            }
        })
        .collect();
        Self { records }
    }
    fn record(&self, name: &str) -> &Record {
        self.records
            .iter()
            .find(|record| record.name == name)
            .unwrap()
    }
    fn sources<'a>(&'a self, reaction: &'a ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let record = self.record(term.species.as_str());
                FormationGibbsSource {
                    species: &term.species,
                    pack: &record.pack,
                    at: &record.point,
                    reference_axis: REFERENCE,
                }
            })
            .collect()
    }
    fn reaction(&self, name: &str, terms: &[(&str, i128)]) -> ConservedReaction {
        let records: Vec<_> = terms.iter().map(|(name, _)| self.record(name)).collect();
        let species: Vec<_> = records.iter().map(|record| id(record.name)).collect();
        let reaction = ReactionId::new(name).unwrap();
        let elemental = ElementalMatrix::new(
            ["H", "C", "N", "O", "Ca"]
                .into_iter()
                .map(|element| ElementId::new(element).unwrap())
                .collect(),
            species.clone(),
            (0..5)
                .map(|row| records.iter().map(|record| record.atoms[row]).collect())
                .collect(),
        )
        .unwrap();
        let stoichiometry = StoichiometricMatrix::new(
            species.clone(),
            vec![reaction.clone()],
            terms
                .iter()
                .map(|(_, coefficient)| vec![*coefficient])
                .collect(),
        )
        .unwrap();
        let charges = ChargeVector::new(
            species,
            records.iter().map(|record| record.charge).collect(),
        )
        .unwrap();
        ConservedReaction::new(&elemental, &stoichiometry, &charges, reaction).unwrap()
    }
    fn aqueous(&self, include_cyanide: bool) -> IdealAcidBaseSystem {
        let reaction = MineralReaction::WaterDissociation.build().unwrap();
        let water = WaterAutoprotolysis::from_sources(
            reaction.clone(),
            id(WATER),
            id(H),
            id(OH),
            &self.sources(&reaction),
        )
        .unwrap();
        let mut families = Vec::new();
        if include_cyanide {
            let reaction = MineralReaction::HydrogenCyanideDissociation
                .build()
                .unwrap();
            let step = AcidDissociation::from_sources(
                reaction.clone(),
                id(HC),
                id(CN),
                id(H),
                &self.sources(&reaction),
            )
            .unwrap();
            families.push(AcidFamily::new(FAMILY, vec![step]).unwrap());
        }
        IdealAcidBaseSystem::new(water, families).unwrap()
    }
}

pub fn salt(temperature: f64, kdiss: f64) -> SaltPrecipitation {
    let fixture = Fixture::new(temperature, kdiss, 1e-9);
    let reaction = fixture.reaction(
        "authored-calcium-cyanide-precipitation",
        &[(CA, -1), (CN, -2), (SALT, 1)],
    );
    SaltPrecipitation::from_sources(
        fixture.aqueous(true),
        reaction.clone(),
        id(CA),
        id(SALT),
        FAMILY,
        &fixture.sources(&reaction),
    )
    .unwrap()
}

pub fn hydrolysis(temperature: f64, kprec: f64) -> HydrolysisPrecipitation {
    let fixture = Fixture::new(temperature, 4e-6, kprec);
    let reaction = MineralReaction::CalciumHydroxidePrecipitation
        .build()
        .unwrap();
    HydrolysisPrecipitation::from_sources(
        fixture.aqueous(false),
        reaction.clone(),
        id(CA),
        id(HYDROXIDE),
        &fixture.sources(&reaction),
    )
    .unwrap()
}

fn point(temperature: f64, phase: &str, ionic: bool) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![REFERENCE, "source-standard-state", phase];
    if phase == "source-phase-aqueous" {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
        if ionic {
            flags.extend([
                "source-individual-ion",
                "source-ionic-reference-hydrogen-zero",
            ]);
        }
    } else {
        flags.push("source-pure-species");
    }
    for flag in flags {
        point = point
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    point
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
            source: format!("authored electrolyte material fixture: {name}"),
            license: "project test fixture".into(),
            artifact: None,
        },
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    claim.observations.push(
        claims
            .register_observation(ObservationDataset {
                specimen: name.into(),
                method: "independently manufactured mineral/acid equilibrium".into(),
                artifact: hash,
                caveats: "not measured solubility, pH or thermal data".into(),
                provenance: claim.provenance.clone(),
            })
            .unwrap(),
    );
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        name,
        "fixture-v1",
        hash,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
