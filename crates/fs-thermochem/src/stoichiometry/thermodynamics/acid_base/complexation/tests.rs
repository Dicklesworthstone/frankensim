//! G0/G3 simultaneous known states from authored exact-chemistry/Gibbs fixtures.
use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const REFERENCE: &str = "source-elemental-reference-authored-metal-ligand-test";
const H: &str = "H+(aq)";
const OH: &str = "OH[-](aq)";
const WATER: &str = "H2O(l)";
const HL: &str = "HL(aq)";
const L: &str = "L[-](aq)";
const T: f64 = 298.15;
const KA: f64 = 1e-7;
const KW: f64 = 1e-14;
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
#[derive(Clone)]
struct Record {
    name: String,
    atoms: Vec<(&'static str, i128)>,
    charge: i128,
    gibbs: f64,
    point: QueryPoint,
    pack: NormalizedPack,
}
fn point(liquid: bool, charge: i128) -> QueryPoint {
    let mut p = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            T,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![
        REFERENCE,
        "source-standard-state",
        if liquid {
            "source-phase-liquid"
        } else {
            "source-phase-aqueous"
        },
        if liquid {
            "source-pure-species"
        } else {
            "source-aqueous-unit-activity-infinite-dilution"
        },
    ];
    if charge != 0 {
        flags.extend([
            "source-individual-ion",
            "source-ionic-reference-hydrogen-zero",
        ]);
    }
    for flag in flags {
        p = p
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    p
}
fn pack(name: &str, g: f64, p: &QueryPoint) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (a, v) in p.axes() {
        validity = validity.with_quantity(a, p.axis_quantities()[a], *v, *v);
    }
    let mut claim = PropertyClaim {
        key: PropertyKey::new(
            "standard-molar-formation-gibbs-energy",
            MolarEnergyQuantityV1::DIMS,
        ),
        value: PropertyValue::Scalar {
            value: if g == 0.0 { 0.0 } else { g },
            dims: MolarEnergyQuantityV1::DIMS,
        },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("Authored analytical metal/ligand fixture {name}"),
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
                method: "manufactured thermodynamic constant".into(),
                artifact: hash,
                caveats: "not measured copper/cyanide data".into(),
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
impl Record {
    fn new(name: &str, atoms: Vec<(&'static str, i128)>, charge: i128, gibbs: f64) -> Self {
        let point = point(name == WATER, charge);
        let pack = pack(name, gibbs, &point);
        Self {
            name: name.into(),
            atoms,
            charge,
            gibbs,
            point,
            pack,
        }
    }
}
#[derive(Clone)]
struct Fixture {
    records: Vec<Record>,
    qs: Vec<usize>,
}
impl Fixture {
    fn new(qs: &[usize], ks: &[f64]) -> Self {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T;
        let ligand_g = -rt * KA.ln();
        let mut records = vec![
            Record::new(H, vec![("H", 1)], 1, 0.0),
            Record::new(OH, vec![("H", 1), ("O", 1)], -1, -rt * KW.ln()),
            Record::new(WATER, vec![("H", 2), ("O", 1)], 0, 0.0),
            Record::new(HL, vec![("H", 1), ("C", 1), ("N", 1)], 0, 0.0),
            Record::new(L, vec![("C", 1), ("N", 1)], -1, ligand_g),
            Record::new("HB(aq)", vec![("H", 1), ("S", 1)], 0, 0.0),
            Record::new("B[-](aq)", vec![("S", 1)], -1, ligand_g),
            Record::new("M0", vec![("Cu", 1)], 1, 0.0),
            Record::new("Z0", vec![("Zn", 1)], 2, 0.0),
            Record::new(
                "Z1",
                vec![("Zn", 1), ("C", 1), ("N", 1)],
                1,
                ligand_g - rt * 20.0_f64.ln(),
            ),
        ];
        let (mut n, mut g) = (0usize, 0.0);
        for (i, (&q, &k)) in qs.iter().zip(ks).enumerate() {
            n += q;
            g += (q as f64) * ligand_g - rt * k.ln();
            records.push(Record::new(
                &format!("M{}", i + 1),
                vec![("Cu", 1), ("C", n as i128), ("N", n as i128)],
                1 - n as i128,
                g,
            ));
        }
        Self {
            records,
            qs: qs.to_vec(),
        }
    }
    fn reaction(&self, name: &str, terms: &[(&str, i128)]) -> ConservedReaction {
        let members: Vec<_> = terms
            .iter()
            .map(|(n, _)| self.records.iter().find(|r| r.name == *n).unwrap())
            .collect();
        let species: Vec<_> = members.iter().map(|r| id(&r.name)).collect();
        let elements = ["H", "C", "N", "O", "Cu", "Zn", "S"];
        let a = ElementalMatrix::new(
            elements
                .iter()
                .map(|e| ElementId::new(*e).unwrap())
                .collect(),
            species.clone(),
            elements
                .iter()
                .map(|e| {
                    members
                        .iter()
                        .map(|r| {
                            r.atoms
                                .iter()
                                .find(|(name, _)| name == e)
                                .map_or(0, |(_, n)| *n)
                        })
                        .collect()
                })
                .collect(),
        )
        .unwrap();
        let rid = ReactionId::new(name).unwrap();
        let matrix = StoichiometricMatrix::new(
            species.clone(),
            vec![rid.clone()],
            terms.iter().map(|(_, n)| vec![*n]).collect(),
        )
        .unwrap();
        let charges =
            ChargeVector::new(species, members.iter().map(|r| r.charge).collect()).unwrap();
        ConservedReaction::new(&a, &matrix, &charges, rid).unwrap()
    }
    fn sources<'a>(&'a self, reaction: &'a ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
        reaction
            .terms()
            .iter()
            .map(|t| {
                let r = self
                    .records
                    .iter()
                    .find(|r| r.name == t.species.as_str())
                    .unwrap();
                FormationGibbsSource {
                    species: &t.species,
                    pack: &r.pack,
                    at: &r.point,
                    reference_axis: REFERENCE,
                }
            })
            .collect()
    }
    fn aqueous(&self, buffer: bool) -> IdealAcidBaseSystem {
        let w = self.reaction("water", &[(WATER, -1), (H, 1), (OH, 1)]);
        let water = WaterAutoprotolysis::from_sources(
            w.clone(),
            id(WATER),
            id(H),
            id(OH),
            &self.sources(&w),
        )
        .unwrap();
        let acid = |name: &str, ha: &str, a: &str| {
            let r = self.reaction(name, &[(ha, -1), (H, 1), (a, 1)]);
            AcidFamily::new(
                name,
                vec![
                    AcidDissociation::from_sources(
                        r.clone(),
                        id(ha),
                        id(a),
                        id(H),
                        &self.sources(&r),
                    )
                    .unwrap(),
                ],
            )
            .unwrap()
        };
        let mut acids = vec![acid("ligand", HL, L)];
        if buffer {
            acids.push(acid("buffer", "HB(aq)", "B[-](aq)"));
        }
        IdealAcidBaseSystem::new(water, acids).unwrap()
    }
    fn family(&self, zinc: bool) -> MetalLigandFamily {
        let (prefix, qs) = if zinc {
            ("Z", vec![1])
        } else {
            ("M", self.qs.clone())
        };
        let mut steps = vec![];
        for (i, q) in qs.iter().enumerate() {
            let parent = format!("{prefix}{i}");
            let complex = format!("{prefix}{}", i + 1);
            let r = self.reaction(
                &format!("{prefix}-association-{i}"),
                &[(&parent, -1), (L, -(*q as i128)), (&complex, 1)],
            );
            steps.push(
                LigandAssociation::from_sources(
                    r.clone(),
                    id(&parent),
                    id(L),
                    id(&complex),
                    &self.sources(&r),
                )
                .unwrap(),
            );
        }
        MetalLigandFamily::new(if zinc { "zinc" } else { "copper" }, steps).unwrap()
    }
    fn model(&self, buffer: bool, zinc: bool) -> IdealLigandSystem {
        let mut families = vec![self.family(false)];
        if zinc {
            families.push(self.family(true));
        }
        IdealLigandSystem::new(self.aqueous(buffer), "ligand", families).unwrap()
    }
}
fn ion(name: &str, z: i32, n: f64) -> SpectatorIon {
    SpectatorIon {
        species: id(name),
        charge: z,
        amount: Amount::new(n),
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-8 * expected.abs().max(1e-100),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn solve(
    model: &IdealLigandSystem,
    ligand: f64,
    metals: &[(&str, f64)],
    other: &[(&str, f64)],
    ions: &[SpectatorIon],
    basis: f64,
) -> LigandSpeciation {
    model
        .speciate(
            Amount::new(ligand),
            &metals
                .iter()
                .map(|(id, n)| (id.to_string(), Amount::new(*n)))
                .collect::<Vec<_>>(),
            &other
                .iter()
                .map(|(id, n)| (id.to_string(), Amount::new(*n)))
                .collect::<Vec<_>>(),
            ions,
            Amount::new(basis),
            LigandSpeciationOptions::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn manufactured_one_ligand_state_closes_charge_and_finite_component_balances() {
    let fx = Fixture::new(&[1], &[4.0]);
    let model = fx.model(false, false);
    let (l, h, mt) = (0.2, 1e-6, 0.3);
    let root = mt / (1.0 + 4.0 * l);
    let complex = mt - root;
    let total = l + l * h / KA + complex;
    let sodium = l + KW / h - h - root;
    let r = solve(
        &model,
        total,
        &[("copper", mt)],
        &[],
        &[ion("Na+", 1, sodium)],
        1.0,
    );
    close(r.free_ligand.value(), l);
    close(r.aqueous.ph.value(), 6.0);
    close(r.metals[0].members[0].amount.value(), root);
    close(r.metals[0].members[1].amount.value(), complex);
    close(r.bound_ligand.value(), complex);
    close(r.uncomplexed_ligand.value(), l + l * h / KA);
    assert!(
        r.ligand_balance_residual.value().abs() + r.ligand_balance_numerical_spread.value()
            <= 1e-15 + 1e-10 * total
    );
    assert!(r.aqueous.charge_balance_residual.value().abs() < 1e-10);
    assert!(r.max_log_mass_action_residual.value() < 1e-8);
    assert!(r.metals[0].members[1].charge == 0);
    assert!(
        !r.aqueous
            .spectators
            .iter()
            .any(|s| s.species.as_str() == "M1")
    );
}

#[test]
fn strongly_cancelled_charge_uses_component_elimination_before_final_aqueous_solve() {
    let fx = Fixture::new(&[1, 1], &[10.0, 5.0]);
    let model = fx.model(false, false);
    let options = LigandSpeciationOptions {
        log_ligand_tolerance: Dimensionless::new(1e-10),
        absolute_ligand_tolerance: Amount::new(1e-11),
        relative_ligand_tolerance: Dimensionless::new(1e-10),
        log_mass_action_tolerance: Dimensionless::new(1e-9),
        ..LigandSpeciationOptions::default()
    };
    let r = model
        .speciate(
            Amount::new(2.0),
            &[("copper".into(), Amount::new(0.2))],
            &[],
            &[ion("Na+", 1, 1.0), ion("Cl-", -1, 0.2)],
            Amount::new(1.0),
            options,
            || false,
        )
        .unwrap();
    // Independent high-precision root of ligand balance and charge, with
    // weights [1, 10*l, 50*l^2], Ka=1e-7 and Kw=1e-14. The water charge
    // is ~8.6e-8 mol after cancellation of mol-scale nonwater ions.
    assert!((r.free_ligand.value() - 0.658_655_599_643_603).abs() < 2e-10);
    assert!((r.aqueous.ph.value() - 6.818_658_425_716_493).abs() < 2e-9);
    assert!((r.bound_ligand.value() - 0.341_344_486_315_231_66).abs() < 2e-10);
    assert!(r.ligand_balance_residual.value().abs() < 2e-14);
    assert_eq!(r.ligand_balance_numerical_spread.value(), 0.0);
    assert!(r.aqueous.charge_balance_residual.value().abs() < 1e-11);
    assert!(r.max_log_mass_action_residual.value() < 1e-9);
    assert!(
        det::ln(r.free_ligand_bracket[1].value()) - det::ln(r.free_ligand_bracket[0].value())
            <= options.log_ligand_tolerance.value()
    );
    // Cancellation must not turn the same analytical component difference
    // into a spurious zero or inaccurate trace pool.
    assert!(resolved_component_difference(1.0, 1.0, 1e-9).is_none());
    assert!(resolved_component_difference(1.0, 1.0 - 1e-12, 1e-9).is_none());
}

#[test]
fn two_then_one_then_one_ligands_and_buffer_share_one_free_pool_and_ph() {
    let fx = Fixture::new(&[2, 1, 1], &[1e7, 1e3, 100.0]);
    let model = fx.model(true, false);
    let (l, h, mt, buffer) = (1e-3, 1e-7, 0.002, 0.001);
    let weights = [1.0, 10.0, 10.0, 1.0];
    let counts = [0, 2, 3, 4];
    let bound = mt * 54.0 / 22.0;
    let total = l + l * h / KA + bound;
    let sodium = l + buffer * 0.5 - (mt - bound) + KW / h - h;
    let r = solve(
        &model,
        total,
        &[("copper", mt)],
        &[("buffer", buffer)],
        &[ion("Na+", 1, sodium)],
        1.0,
    );
    close(r.free_ligand.value(), l);
    close(r.aqueous.ph.value(), 7.0);
    close(r.bound_ligand.value(), bound);
    for (j, m) in r.metals[0].members.iter().enumerate() {
        close(m.amount.value(), mt * weights[j] / 22.0);
        assert_eq!(m.ligand_count, counts[j]);
        assert_eq!(m.charge, 1 - counts[j] as i32);
    }
    let b = r
        .aqueous
        .families
        .iter()
        .find(|f| f.id == "buffer")
        .unwrap();
    close(b.forms[0].amount.value(), buffer * 0.5);
    assert_eq!(r.model().families()[0].steps()[0].ligand_count(), 2);
}

#[test]
fn different_metals_compete_independently_of_input_order_and_extensive_scale() {
    let fx = Fixture::new(&[1], &[10.0]);
    let model = fx.model(false, true);
    let base = solve(
        &model,
        0.5,
        &[("copper", 0.2), ("zinc", 0.3)],
        &[],
        &[ion("Cl-", -1, 0.4)],
        1.0,
    );
    close(base.free_ligand.value(), 0.1);
    close(base.aqueous.ph.value(), 7.0);
    close(base.metals[0].bound_ligand.value(), 0.1);
    close(base.metals[1].bound_ligand.value(), 0.2);
    let reversed = IdealLigandSystem::new(
        fx.aqueous(false),
        "ligand",
        vec![fx.family(true), fx.family(false)],
    )
    .unwrap();
    for scale in [1e-8, 2.0, 1e8] {
        let r = solve(
            &reversed,
            0.5 * scale,
            &[("zinc", 0.3 * scale), ("copper", 0.2 * scale)],
            &[],
            &[ion("Cl-", -1, 0.4 * scale)],
            scale,
        );
        close(r.free_ligand.value() / scale, base.free_ligand.value());
        close(r.bound_ligand.value() / scale, base.bound_ligand.value());
        close(r.aqueous.ph.value(), base.aqueous.ph.value());
    }
    let r = solve(
        &reversed,
        0.5,
        &[("zinc", 0.3), ("copper", 0.2)],
        &[],
        &[ion("Cl-", -1, 0.4)],
        1.0,
    );
    assert_eq!(
        r.free_ligand.value().to_bits(),
        base.free_ligand.value().to_bits()
    );
}

#[test]
fn zero_and_trace_pools_keep_original_aqueous_behavior_and_metal_stocks() {
    let fx = Fixture::new(&[1], &[10.0]);
    let model = fx.model(false, false);
    let empty = solve(
        &model,
        0.0,
        &[("copper", 0.2)],
        &[],
        &[ion("Cl-", -1, 0.2)],
        1.0,
    );
    assert_eq!(empty.free_ligand.value(), 0.0);
    assert_eq!(empty.bound_ligand.value(), 0.0);
    assert_eq!(empty.metals[0].members[0].amount.value(), 0.2);
    assert_eq!(empty.metals[0].members[1].amount.value(), 0.0);
    close(empty.aqueous.ph.value(), 7.0);
    let ions = [ion("Na+", 1, 0.2)];
    let no_metal = solve(&model, 0.2, &[("copper", 0.0)], &[], &ions, 1.0);
    let original = fx
        .aqueous(false)
        .speciate(
            &[("ligand".into(), Amount::new(0.2))],
            &ions,
            Amount::new(1.0),
            LigandSpeciationOptions::default().aqueous,
            || false,
        )
        .unwrap();
    assert_eq!(
        no_metal.aqueous.ph.value().to_bits(),
        original.ph.value().to_bits()
    );
    assert_eq!(no_metal.iterations, 0);
    let trace = solve(&model, 0.2, &[("copper", 1e-30)], &[], &ions, 1.0);
    assert!(
        trace.metals[0]
            .members
            .iter()
            .all(|m| m.amount.value() > 0.0)
    );
    close(
        trace.metals[0]
            .members
            .iter()
            .map(|m| m.amount.value())
            .sum(),
        1e-30,
    );
}

#[test]
fn insufficient_inner_accuracy_brackets_work_and_cancellation_refuse_atomically() {
    let model = Fixture::new(&[1], &[4.0]).model(false, false);
    let metals = vec![("copper".into(), Amount::new(0.3))];
    let ions = vec![ion("Na+", 1, 0.03333234333333335)];
    let before = metals.clone();
    let run = |options| {
        model.speciate(
            Amount::new(2.3333333333333335),
            &metals,
            &[],
            &ions,
            Amount::new(1.0),
            options,
            || false,
        )
    };
    for o in [
        LigandSpeciationOptions {
            max_iterations: 1,
            ..Default::default()
        },
        LigandSpeciationOptions {
            max_aqueous_solves: 1,
            ..Default::default()
        },
        LigandSpeciationOptions {
            minimum_free_ligand: Amount::new(2.0),
            ..Default::default()
        },
    ] {
        assert!(run(o).is_err());
    }
    let loose = LigandSpeciationOptions {
        aqueous: SpeciationOptions {
            ph_tolerance: Dimensionless::new(0.5),
            absolute_charge_tolerance: Amount::new(100.0),
            relative_charge_tolerance: Dimensionless::new(0.01),
            ..Default::default()
        },
        log_ligand_tolerance: Dimensionless::new(1e-12),
        absolute_ligand_tolerance: Amount::new(1e-14),
        relative_ligand_tolerance: Dimensionless::new(0.0),
        log_mass_action_tolerance: Dimensionless::new(1e-12),
        ..Default::default()
    };
    assert!(run(loose).is_err());
    let mut checks = 0;
    assert!(matches!(
        model.speciate(
            Amount::new(2.3333333333333335),
            &metals,
            &[],
            &ions,
            Amount::new(1.0),
            Default::default(),
            || {
                checks += 1;
                checks > 15
            }
        ),
        Err(LigandSpeciationError::Cancelled)
    ));
    assert_eq!(metals, before);
    assert!(
        model
            .speciate(
                Amount::new(1.0),
                &[],
                &[],
                &ions,
                Amount::new(1.0),
                Default::default(),
                || false
            )
            .is_err()
    );
    assert!(
        model
            .speciate(
                Amount::new(1.0),
                &metals,
                &[],
                &[ion("M0", 1, 0.1)],
                Amount::new(1.0),
                Default::default(),
                || false
            )
            .is_err()
    );
}

#[test]
fn reused_sources_exact_roles_and_nonoverlapping_components_are_required() {
    let fx = Fixture::new(&[1, 1], &[10.0, 5.0]);
    let original = fx.aqueous(false);
    let mut changed = fx.clone();
    let ligand = changed.records.iter_mut().find(|r| r.name == L).unwrap();
    ligand.gibbs += 1.0;
    ligand.pack = pack(L, ligand.gibbs, &ligand.point);
    assert!(
        IdealLigandSystem::new(original.clone(), "ligand", vec![changed.family(false)]).is_err()
    );
    let family = fx.family(false);
    assert!(
        MetalLigandFamily::new("broken", family.steps().iter().rev().cloned().collect()).is_err()
    );
    assert!(
        IdealLigandSystem::new(original.clone(), "ligand", vec![family.clone(), family]).is_err()
    );
    let r = fx.reaction("scaled", &[("M0", -2), (L, -2), ("M1", 2)]);
    assert!(
        LigandAssociation::from_sources(r.clone(), id("M0"), id(L), id("M1"), &fx.sources(&r))
            .is_err()
    );
    assert!(IdealLigandSystem::new(original, "absent", vec![fx.family(false)]).is_err());
    let mut wrong = fx.clone();
    let metal = wrong.records.iter_mut().find(|r| r.name == "M0").unwrap();
    metal.point = point(true, 0);
    metal.pack = pack("M0", metal.gibbs, &metal.point);
    let r = wrong.reaction("wrong-phase", &[("M0", -1), (L, -1), ("M1", 1)]);
    assert!(
        LigandAssociation::from_sources(r.clone(), id("M0"), id(L), id("M1"), &wrong.sources(&r))
            .is_err()
    );
}
