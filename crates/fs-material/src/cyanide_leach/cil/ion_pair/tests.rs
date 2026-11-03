use super::*;

fn declared_pair(
    ligand_count: i128,
    counterion_charge: i128,
    counterion_gold: i128,
    scale: i128,
) -> (ConservedReaction, SpeciesId, SpeciesId) {
    let counterion = SpeciesId::new("declared-counterion").unwrap();
    let pair = SpeciesId::new("declared-adsorbate").unwrap();
    let species = vec![
        SpeciesId::new(COMPLEX).unwrap(),
        counterion.clone(),
        pair.clone(),
    ];
    let id = ReactionId::new("declared-gold-pair").unwrap();
    let elemental = ElementalMatrix::new(
        ["Au", "C", "N", "Na"]
            .into_iter()
            .map(|e| ElementId::new(e).unwrap())
            .collect(),
        species.clone(),
        vec![
            vec![1, counterion_gold, 1 + counterion_gold],
            vec![ligand_count, 0, ligand_count],
            vec![ligand_count, 0, ligand_count],
            vec![0, 1, 1],
        ],
    )
    .unwrap();
    let stoichiometric = StoichiometricMatrix::new(
        species.clone(),
        vec![id.clone()],
        vec![vec![-scale], vec![-scale], vec![scale]],
    )
    .unwrap();
    let charge =
        ChargeVector::new(species, vec![-1, counterion_charge, counterion_charge - 1]).unwrap();
    (
        ConservedReaction::new(&elemental, &stoichiometric, &charge, id).unwrap(),
        counterion,
        pair,
    )
}

#[test]
fn sodium_pair_retains_exact_admitted_atoms_charge_and_extent() {
    let pair = IonPairTransfer::sodium().unwrap();
    let reaction = pair.reaction();
    assert_eq!(reaction.species_charge(pair.counterion()), Some(1));
    assert_eq!(reaction.species_charge(pair.adsorbed_pair()), Some(0));
    let complex = SpeciesId::new(COMPLEX).unwrap();
    assert_eq!(reaction.species_charge(&complex), Some(-1));
    let atoms: Vec<_> = reaction
        .species_elements(pair.adsorbed_pair())
        .unwrap()
        .iter()
        .map(|(element, count)| (element.as_str(), *count))
        .collect();
    assert_eq!(atoms, [("Au", 1), ("C", 2), ("N", 2), ("Na", 1)]);
    let unknown = SpeciesId::new("not-in-reaction").unwrap();
    assert_eq!(reaction.species_charge(&unknown), None);
    assert_eq!(reaction.species_elements(&unknown), None);
    for amount in reaction.amounts(fs_qty::Amount::new(0.25)).unwrap() {
        assert_eq!(amount.amount.value(), 0.25);
    }
}

#[test]
fn chemical_admission_uses_artifacts_instead_of_species_name_suffixes() {
    let (reaction, cation, adsorbed) = declared_pair(2, 1, 0, 1);
    let admitted =
        IonPairTransfer::from_reaction(reaction, cation.clone(), adsorbed.clone()).unwrap();
    assert_eq!(admitted.counterion(), &cation);
    assert_eq!(admitted.adsorbed_pair(), &adsorbed);
    // Every following equation conserves atoms and charge, but is not the
    // specific neutral 1:1 Au(CN)2-/monovalent-cation model admitted by CIL.
    for (ligands, charge, cation_gold, scale) in [
        (2, 2, 0, 1), // conserved, but cation is divalent and pair remains charged
        (1, 1, 0, 1), // conserved counterfeit formula for original Au(CN)2-
        (2, 1, 1, 1), // a second gold atom hidden in the countercation
        (2, 1, 0, 2), // different written extent convention
    ] {
        let (reaction, cation, adsorbed) = declared_pair(ligands, charge, cation_gold, scale);
        assert!(IonPairTransfer::from_reaction(reaction, cation, adsorbed).is_err());
    }
}

#[test]
fn pair_roles_require_distinct_present_species() {
    let pair = IonPairTransfer::sodium().unwrap();
    assert!(
        IonPairTransfer::from_reaction(
            pair.reaction().clone(),
            pair.counterion().clone(),
            pair.counterion().clone(),
        )
        .is_err()
    );
    assert!(
        IonPairTransfer::from_reaction(
            pair.reaction().clone(),
            SpeciesId::new("absent-counterion").unwrap(),
            pair.adsorbed_pair().clone(),
        )
        .is_err()
    );
}
