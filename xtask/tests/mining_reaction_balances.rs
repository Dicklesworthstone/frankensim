//! G0/G1/G3: real mineral/reagent packs feed the shared reaction-heat closure.
//! A caller's labels cannot supply missing source reference-state authority.
#![deny(unsafe_code)]

use std::{path::Path, process::Command};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dims, Mass, MolarMass, QuantitySpec};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{FormationEnthalpySource, ReactionBalanceError, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-janaf";
const SLUGS: [&str; 3] = [
    "calcium-oxide-solid-nist-mining",
    "water-liquid-janaf-mining",
    "calcium-hydroxide-solid-nist-mining",
];

fn compile(slug: &str, omit_formation_axis: Option<&str>) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "frankensim-reaction-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&out).unwrap();
    let source = root.join(format!("data/matdb/seed-v1/{slug}"));
    let manifest = if let Some(axis) = omit_formation_axis {
        // Intentional mutation: the offered datum keeps its numeric value but
        // loses one source constraint. Extra query labels cannot repair that.
        std::fs::copy(source.join("manifest.tsv"), out.join("manifest.tsv")).unwrap();
        let text = std::fs::read_to_string(source.join("properties.tsv")).unwrap();
        let text = text
            .lines()
            .filter(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                !(fields.first() == Some(&"validity")
                    && fields.get(1) == Some(&"formation")
                    && fields.get(2) == Some(&axis))
            })
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        std::fs::write(out.join("properties.tsv"), text).unwrap();
        out.join("manifest.tsv")
    } else {
        source.join("manifest.tsv")
    };
    let output = out.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(manifest)
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap()
}

fn point(liquid: bool) -> QueryPoint {
    phase_point(if liquid {
        "source-phase-liquid"
    } else {
        "source-phase-solid"
    })
}

fn phase_point(phase_axis: &str) -> QueryPoint {
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
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (phase_axis, 1.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn ids() -> [SpeciesId; 3] {
    ["CaO(s)", "H2O(l)", "Ca(OH)2(s)"].map(|name| SpeciesId::new(name).unwrap())
}

fn sources<'a, const N: usize>(
    packs: &'a [NormalizedPack; N],
    ids: &'a [SpeciesId; N],
    points: &'a [QueryPoint; N],
) -> [FormationEnthalpySource<'a>; N] {
    std::array::from_fn(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    })
}

#[test]
fn g1_g3_sourced_slaking_heat_reaches_conserved_reaction_and_mass_basis() {
    let packs = SLUGS.map(|slug| compile(slug, None));
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let heat = reaction
        .standard_enthalpy_change(&sources(&packs, &ids, &points))
        .unwrap();
    // Independently read source values: -986090 - (-635090) - (-285830).
    assert_eq!(heat.molar().value(), -65_170.0);
    assert_eq!(heat.terms().len(), 3);
    for term in heat.terms() {
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(!term.answer.receipt.source_hashes.is_empty());
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
    }
    let mass = packs[0]
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            &points[0],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &ids[0],
            Mass::new(56.077),
            MolarMass::new(mass.evidence.value.value),
        )
        .unwrap();
    assert!((extent.value() - 1000.0).abs() < 1e-10);
    assert!((heat.at_extent(extent).unwrap().value() + 65_170_000.0).abs() < 1e-7);
    let transfers = reaction.amounts(extent).unwrap();
    let water = transfers
        .iter()
        .find(|term| term.term.species == ids[1])
        .unwrap();
    assert_eq!(water.term.coefficient, -1);
    assert!((water.amount.value() - 1000.0).abs() < 1e-10);

    // Caller ordering never changes the canonical summation/receipt order.
    let mut reversed = sources(&packs, &ids, &points);
    reversed.reverse();
    let other = reaction.standard_enthalpy_change(&reversed).unwrap();
    assert_eq!(
        other.molar().value().to_bits(),
        heat.molar().value().to_bits()
    );
    assert_eq!(
        other.terms().iter().map(|t| &t.term).collect::<Vec<_>>(),
        heat.terms().iter().map(|t| &t.term).collect::<Vec<_>>()
    );
    assert_eq!(heat.at_extent(Amount::new(0.0)).unwrap().value(), 0.0);
}

#[test]
fn g0_source_must_declare_reference_axes_not_only_accept_caller_labels() {
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    for omitted in [
        "temperature",
        "pressure",
        "reference_convention",
        "standard_state",
        "phase",
    ] {
        let packs = [
            compile(SLUGS[0], Some(omitted)),
            compile(SLUGS[1], None),
            compile(SLUGS[2], None),
        ];
        // The generic material owner may answer an underconstrained datum.
        // Reaction thermochemistry must additionally prove its common reference.
        packs[0]
            .claims()
            .query_typed(
                &PropertyKey::with_quantity(
                    "standard-molar-formation-enthalpy",
                    QuantitySpec::dimensional(Dims([2, 1, -2, 0, 0, -1])),
                ),
                &points[0],
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!(
            matches!(
                reaction.standard_enthalpy_change(&sources(&packs, &ids, &points)),
                Err(ReactionBalanceError::Invalid(_))
            ),
            "accepted source lacking {omitted}"
        );
    }
}

#[test]
fn g3_reaction_refuses_incomplete_duplicate_or_incompatible_source_contexts() {
    let packs = SLUGS.map(|slug| compile(slug, None));
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let mut offered = sources(&packs, &ids, &points);
    assert!(reaction.standard_enthalpy_change(&offered[..2]).is_err());
    offered[2].species = &ids[0];
    assert!(matches!(
        reaction.standard_enthalpy_change(&offered),
        Err(ReactionBalanceError::Species { .. })
    ));
    let mut offered = sources(&packs, &ids, &points);
    offered[0].reference_axis = "source-phase-solid";
    assert!(reaction.standard_enthalpy_change(&offered).is_err());
    for changed in [
        point(true),
        point(false).with("temperature", 300.0).unwrap(),
        point(false).with("pressure", 101_325.0).unwrap(),
        point(false).with(REFERENCE, 0.0).unwrap(),
        point(false).with("source-pure-species", 0.0).unwrap(),
    ] {
        let changed_points = [changed, point(true), point(false)];
        assert!(
            reaction
                .standard_enthalpy_change(&sources(&packs, &ids, &changed_points))
                .is_err()
        );
    }
}

#[test]
fn g1_janaf_magnesite_calcination_closes_carbon_dioxide_and_standard_heat() {
    let slugs = [
        "magnesite-solid-janaf-mining",
        "magnesium-oxide-solid-janaf-mining",
        "carbon-dioxide-gas-janaf-mining",
    ];
    let packs = slugs.map(|slug| compile(slug, None));
    let ids = ["MgCO3(s)", "MgO(s)", "CO2(g)"].map(|name| SpeciesId::new(name).unwrap());
    let points = [point(false), point(false), phase_point("source-phase-gas")];
    let reaction = MineralReaction::MagnesiteCalcination.build().unwrap();
    let heat = reaction
        .standard_enthalpy_change(&sources(&packs, &ids, &points))
        .unwrap();
    // JANAF C-079, Mg-008, C-095: products - reactant, literal 298.15 K rows.
    assert_eq!(heat.molar().value(), 116_926.0);
    let molar_mass = packs[0]
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            &points[0],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    packs[0]
        .claims()
        .verify_receipt(&molar_mass.receipt)
        .unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &ids[0],
            Mass::new(84.3139),
            MolarMass::new(molar_mass.evidence.value.value),
        )
        .unwrap();
    assert!((extent.value() - 1000.0).abs() < 1e-10);
    assert!((heat.at_extent(extent).unwrap().value() - 116_926_000.0).abs() < 1e-7);
    let amounts = reaction.amounts(extent).unwrap();
    let carbon_dioxide = amounts
        .iter()
        .find(|amount| amount.term.species == ids[2])
        .unwrap();
    assert_eq!(carbon_dioxide.term.coefficient, 1);
    assert!((carbon_dioxide.amount.value() - 1000.0).abs() < 1e-10);
    let co2_mass = packs[2]
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            &points[2],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!(
        (carbon_dioxide.amount.value() * co2_mass.evidence.value.value - 44.0095).abs() < 1e-12
    );
    for term in heat.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(!term.answer.receipt.source_hashes.is_empty());
    }
    // Reference-point chemical energy cannot be silently promoted to kiln duty.
    let hot_points = points
        .clone()
        .map(|point| point.with("temperature", 900.0).unwrap());
    assert!(
        reaction
            .standard_enthalpy_change(&sources(&packs, &ids, &hot_points))
            .is_err()
    );
}

#[test]
fn g1_janaf_magnetite_oxidation_uses_four_mineral_moles_per_extent_and_real_oxygen() {
    let slugs = [
        "magnetite-solid-janaf-mining",
        "oxygen-gas-janaf-mining",
        "hematite-solid-janaf-mining",
    ];
    let packs = slugs.map(|slug| compile(slug, None));
    let ids = ["Fe3O4(s)", "O2(g)", "Fe2O3(s)"].map(|name| SpeciesId::new(name).unwrap());
    let points = [point(false), phase_point("source-phase-gas"), point(false)];
    let reaction = MineralReaction::MagnetiteOxidationToHematite
        .build()
        .unwrap();
    let heat = reaction
        .standard_enthalpy_change(&sources(&packs, &ids, &points))
        .unwrap();
    // JANAF Fe-032, Fe-030, O-029: 6*(-825503) - 4*(-1120894) - 0.
    assert_eq!(heat.molar().value(), -469_442.0);
    assert_eq!(
        heat.at_extent(Amount::new(0.25)).unwrap().value(),
        -117_360.5
    );
    let extent = reaction
        .limiting_extent(&[
            (ids[0].clone(), Amount::new(1000.0)),
            (ids[1].clone(), Amount::new(200.0)),
        ])
        .unwrap();
    assert_eq!(extent.value(), 200.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (species, amount) in [(0, 800.0), (1, 200.0), (2, 1200.0)] {
        assert_eq!(
            amounts
                .iter()
                .find(|term| term.term.species == ids[species])
                .unwrap()
                .amount
                .value(),
            amount
        );
    }
    assert_eq!(heat.at_extent(extent).unwrap().value(), -93_888_400.0);
    // Zero elemental formation enthalpy still has an actual source dependency.
    let oxygen = heat
        .terms()
        .iter()
        .find(|term| term.term.species == ids[1])
        .unwrap();
    assert_eq!(oxygen.answer.evidence.value.value, 0.0);
    packs[1]
        .claims()
        .verify_receipt(&oxygen.answer.receipt)
        .unwrap();
    assert!(!oxygen.answer.receipt.source_hashes.is_empty());
    let wrong_phase = [point(false), point(true), point(false)];
    assert!(
        reaction
            .standard_enthalpy_change(&sources(&packs, &ids, &wrong_phase))
            .is_err()
    );
}

#[test]
fn g1_janaf_copper_oxide_reductions_close_reductant_metal_and_heat() {
    for (reaction, oxide, slug, expected_heat, oxide_molar_mass, copper_coefficient) in [
        (
            MineralReaction::CupricOxideCarbonMonoxideReduction,
            "CuO(s)",
            "copper-ii-oxide-solid-janaf-mining",
            -126_932.0,
            0.079545,
            1.0,
        ),
        (
            MineralReaction::CuprousOxideCarbonMonoxideReduction,
            "Cu2O(s)",
            "copper-i-oxide-solid-janaf-mining",
            -112_288.0,
            0.143091,
            2.0,
        ),
    ] {
        let packs = [
            slug,
            "carbon-monoxide-gas-janaf-mining",
            "copper-solid-janaf-mining",
            "carbon-dioxide-gas-janaf-mining",
        ]
        .map(|name| compile(name, None));
        let ids = [oxide, "CO(g)", "Cu(s)", "CO2(g)"].map(|id| SpeciesId::new(id).unwrap());
        let points = [
            point(false),
            phase_point("source-phase-gas"),
            point(false),
            phase_point("source-phase-gas"),
        ];
        let reaction = reaction.build().unwrap();
        let heat = reaction
            .standard_enthalpy_change(&sources(&packs, &ids, &points))
            .unwrap();
        // Independent JANAF rows: CO2=-393522, CO=-110527, Cu=0,
        // CuO=-156063 or Cu2O=-170707, all J/mol at 298.15 K and 1 bar.
        assert_eq!(heat.molar().value(), expected_heat);
        assert_eq!(heat.terms().len(), 4);
        for term in heat.terms() {
            let pack = packs
                .iter()
                .find(|pack| pack.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert_eq!(
                term.answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(!term.answer.receipt.source_hashes.is_empty());
        }
        let masses = std::array::from_fn::<_, 4, _>(|i| {
            let answer = packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::with_quantity(
                        "molar-mass",
                        QuantitySpec::dimensional(MolarMass::DIMS),
                    ),
                    &points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            packs[i].claims().verify_receipt(&answer.receipt).unwrap();
            answer.evidence.value.value
        });
        assert!((masses[0] - oxide_molar_mass).abs() < 1e-15 * oxide_molar_mass);
        let extent = reaction
            .extent_from_reacted_mass(
                &ids[0],
                Mass::new(1000.0 * oxide_molar_mass),
                MolarMass::new(masses[0]),
            )
            .unwrap();
        assert!((extent.value() - 1000.0).abs() < 2e-12);
        assert!((heat.at_extent(extent).unwrap().value() - 1000.0 * expected_heat).abs() < 2e-7);
        let amounts = reaction.amounts(extent).unwrap();
        let amount_of = |index: usize| {
            amounts
                .iter()
                .find(|amount| amount.term.species == ids[index])
                .unwrap()
                .amount
                .value()
        };
        assert!((amount_of(1) * masses[1] - 28.0101).abs() < 1e-11);
        assert!((amount_of(2) * masses[2] - 63.546 * copper_coefficient).abs() < 1e-11);
        assert!((amount_of(3) * masses[3] - 44.0095).abs() < 1e-11);

        // Rounded published masses differ by 0.4 g per kmol extent. Preserve
        // that visible mismatch rather than changing a source to force closure.
        let residual = amount_of(2) * masses[2] + amount_of(3) * masses[3]
            - amount_of(0) * masses[0]
            - amount_of(1) * masses[1];
        assert!((residual - 0.0004).abs() < 1e-11);

        let inventory = [
            (ids[0].clone(), Amount::new(1000.0)),
            (ids[1].clone(), Amount::new(250.0)),
            (ids[2].clone(), Amount::new(0.0)),
            (ids[3].clone(), Amount::new(0.0)),
        ];
        let limit = reaction.limiting_extent(&inventory[..2]).unwrap();
        assert_eq!(limit.value(), 250.0);
        let after = reaction.advance_inventory(&inventory, limit).unwrap();
        let remaining = |index: usize| {
            after
                .iter()
                .find(|(id, _)| *id == ids[index])
                .unwrap()
                .1
                .value()
        };
        assert_eq!(remaining(0), 750.0);
        assert_eq!(remaining(1), 0.0);
        assert_eq!(remaining(2), 250.0 * copper_coefficient);
        assert_eq!(remaining(3), 250.0);
        assert!(
            reaction
                .advance_inventory(&inventory, Amount::new(251.0))
                .is_err()
        );

        // A zero elemental formation heat is still a required real source;
        // neither omitting copper nor using a hot reference is admissible.
        assert!(
            reaction
                .standard_enthalpy_change(&sources(&packs, &ids, &points)[..3])
                .is_err()
        );
        let hot = points
            .clone()
            .map(|point| point.with("temperature", 1000.0).unwrap());
        assert!(
            reaction
                .standard_enthalpy_change(&sources(&packs, &ids, &hot))
                .is_err()
        );
    }
}

#[test]
fn g1_janaf_cuprous_oxidation_retains_oxygen_extent_and_standard_heat() {
    let packs = [
        "copper-i-oxide-solid-janaf-mining",
        "oxygen-gas-janaf-mining",
        "copper-ii-oxide-solid-janaf-mining",
    ]
    .map(|slug| compile(slug, None));
    let ids = ["Cu2O(s)", "O2(g)", "CuO(s)"].map(|id| SpeciesId::new(id).unwrap());
    let points = [point(false), phase_point("source-phase-gas"), point(false)];
    let reaction = MineralReaction::CuprousOxideOxidation.build().unwrap();
    let heat = reaction
        .standard_enthalpy_change(&sources(&packs, &ids, &points))
        .unwrap();
    // 4*(-156063) - 2*(-170707) - 0, per written integer equation.
    assert_eq!(heat.molar().value(), -282_838.0);
    let extent = reaction
        .limiting_extent(&[
            (ids[0].clone(), Amount::new(1000.0)),
            (ids[1].clone(), Amount::new(600.0)),
        ])
        .unwrap();
    assert_eq!(extent.value(), 500.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (index, expected) in [(0, 1000.0), (1, 500.0), (2, 2000.0)] {
        assert_eq!(
            amounts
                .iter()
                .find(|a| a.term.species == ids[index])
                .unwrap()
                .amount
                .value(),
            expected
        );
    }
    assert_eq!(heat.at_extent(extent).unwrap().value(), -141_419_000.0);
    for term in heat.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
    }
    let wrong_phase = [point(false), point(false), point(false)];
    assert!(
        reaction
            .standard_enthalpy_change(&sources(&packs, &ids, &wrong_phase))
            .is_err()
    );
}
