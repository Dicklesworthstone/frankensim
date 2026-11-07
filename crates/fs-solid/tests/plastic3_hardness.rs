//! Native spatial J2 history drives independently sourced local hardness/wear.
//! Synthetic source cards and prescribed surface quadrature isolate coupling;
//! these checks neither calibrate a material nor solve a contact-pressure field.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::{
    ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA, ResolvedJ2Hardness,
};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY, ResolvedJ2Plasticity,
    resolve_j2_plasticity,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};
use fs_solid::plastic3::hardness::{J2HardnessField, J2HardnessObservationError};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};
use fs_tribo::{
    ArchardLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, LocalHardnessArchardSample,
    WearState,
};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const YIELD_TOLERANCE_PA: f64 = 1e-3;
const HARDNESS_SLOPE_PA: f64 = 2e10;
const VIRGIN_HARDNESS_PA: [f64; 2] = [1.6e9, 8e8];

fn material(region: usize) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: format!("synthetic spatial hardness region {region}; no measured claim"),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 2e11),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, [5e8, 2e8][region]),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 1e9),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("synthetic specimen {region}"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", &[region as u8]),
            caveats: "numerical regression only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic spherical indenter, radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic projected-area protocol revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    let base = VIRGIN_HARDNESS_PA[region];
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: ALPHA.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, base), (0.08, base + 0.08 * HARDNESS_SLOPE_PA)],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 300.0)
                .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.08),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic region {region}"),
                phase: "solid".into(),
                process: "synthetic indentation and J2 tests".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        key,
    )
}

struct Fixture {
    cards: [MaterialCard; 2],
    sources: [ResolvedJ2Hardness; 2],
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    laws: Vec<ResolvedJ2Plasticity>,
    prescribed: Vec<(usize, f64)>,
}

impl Fixture {
    fn new() -> Self {
        let [(first, key0), (second, key1)] = [material(0), material(1)];
        let cards = [first, second];
        let sources = [(key0, &cards[0]), (key1, &cards[1])].map(|(key, card)| {
            let law = resolve_j2_plasticity(
                card,
                &QueryPoint::new().with("T", 300.0).unwrap(),
                SELECT,
                HardeningPropertyBasis::PlasticStrainModulus,
                0.05,
            )
            .unwrap();
            ResolvedJ2Hardness::new(law, card, key, SELECT, 0.06).unwrap()
        });
        let mut nodes = Vec::new();
        for z in 0..=1 {
            for y in 0..=1 {
                for x in 0..=2 {
                    nodes.push([f64::from(x), f64::from(y), f64::from(z)]);
                }
            }
        }
        let mut tets = Vec::new();
        for x in 0..2 {
            let c = [x, x + 1, x + 3, x + 4, x + 6, x + 7, x + 9, x + 10];
            for tet in [
                [0, 1, 3, 7],
                [0, 1, 7, 5],
                [0, 2, 7, 3],
                [0, 2, 6, 7],
                [0, 4, 5, 7],
                [0, 4, 7, 6],
            ] {
                tets.push(tet.map(|i| c[i]));
            }
        }
        let laws = (0..12)
            .map(|element| sources[element / 6].plastic().clone())
            .collect();
        let mut prescribed = Vec::new();
        for (node, x) in nodes.iter().enumerate() {
            if x[0] == 0.0 {
                prescribed.push((3 * node, 0.0));
            }
            if x[0] == 2.0 {
                prescribed.push((3 * node, 2.0));
            }
            if x[1] == 0.0 {
                prescribed.push((3 * node + 1, 0.0));
            }
            if x[2] == 0.0 {
                prescribed.push((3 * node + 2, 0.0));
            }
        }
        Self {
            cards,
            sources,
            nodes,
            tets,
            laws,
            prescribed,
        }
    }

    fn problem<'a>(&'a self, path: &'a [f64]) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
        SmallStrainTetProblem {
            nodes_m: &self.nodes,
            tetrahedra: &self.tets,
            laws: SmallStrainLawField::PerElement(&self.laws),
            prescribed_m: &self.prescribed,
            nodal_forces_n: &[],
            load_path: path,
            budget: TetAssemblyBudget::standard(),
            settings: IncrementSettings::default(),
        }
    }

    fn source_field(&self) -> Vec<&ResolvedJ2Hardness> {
        (0..12).map(|element| &self.sources[element / 6]).collect()
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 5,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn native_bonded_solid_history_supplies_local_hardness_and_archard_wear() {
    let fixture = Fixture::new();
    let sources = fixture.source_field();
    let path = [0.0005, 0.0035, 0.0034, -0.0035, 0.0035];
    let mut fields = Vec::new();
    for end in 1..=path.len() {
        let problem = fixture.problem(&path[..end]);
        let mechanical = with_cx(|cx| problem.solve(cx)).unwrap();
        let before = mechanical.clone();
        let field = with_cx(|cx| {
            problem.observe_hardness(
                &mechanical,
                J2HardnessField::PerElement(&sources),
                YIELD_TOLERANCE_PA,
                cx,
            )
        })
        .unwrap();
        for element in 0..12 {
            near(
                field[element],
                VIRGIN_HARDNESS_PA[element / 6]
                    + HARDNESS_SLOPE_PA * mechanical.element_states[element].alpha,
            );
            let observed = sources[element]
                .observe_accepted(
                    &fixture.laws[element],
                    &mechanical.element_strain[element],
                    &mechanical.element_states[element],
                    YIELD_TOLERANCE_PA,
                )
                .unwrap();
            let typed = sources[element]
                .resolve_hardness_state(&fixture.cards[element / 6], &observed)
                .unwrap();
            near(typed.hardness_pa(), field[element]);
            assert_eq!(typed.key(), sources[element].hardness_source().key());
            assert_eq!(
                typed.resolved().state_coordinate(ALPHA),
                Some(mechanical.element_states[element].alpha)
            );
        }
        assert_eq!(mechanical, before);
        fields.push(field);
    }
    assert!(
        fields[1]
            .iter()
            .zip(&fields[0])
            .any(|(after, before)| after > before)
    );
    for (unloaded, loaded) in fields[2].iter().zip(&fields[1]) {
        near(*unloaded, *loaded);
    }
    assert!(
        fields[3]
            .iter()
            .zip(&fields[2])
            .any(|(after, before)| after > before)
    );
    assert!(
        fields[4]
            .iter()
            .zip(&fields[3])
            .any(|(after, before)| after > before)
    );

    let interface = InterfaceSystemRef::new(
        "fixture/spatial-solid->counterface",
        "fixture/prescribed-local-sliding",
        "fixture/independent-indentation-tests",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let work: Vec<_> = (0..12)
        .map(|element| {
            // Explicit local p * area * unsigned slip; no pressure is inferred from
            // a bulk equivalent stress and no slip is inferred from plastic strain.
            let pressure_pa = 2e6;
            let area_m2 = 1e-3;
            let slip_m = 1e-3 * f64::from(element + 1);
            pressure_pa * area_m2 * slip_m
        })
        .collect();
    let wear = |hardness: &[f64]| {
        ArchardLaw::advance_local_hardness_pressure_weighted_sliding(
            &interface,
            &mut WearState::default(),
            hardness
                .iter()
                .zip(&work)
                .map(
                    |(&hardness_pa, &pressure_weighted_sliding_j)| LocalHardnessArchardSample {
                        wear_coefficient: 1e-5,
                        pressure_weighted_sliding_j,
                        hardness_pa,
                    },
                ),
        )
        .unwrap()
        .increment_m3
    };
    let worn = wear(&fields[4]);
    near(
        worn,
        fields[4]
            .iter()
            .zip(&work)
            .map(|(h, work)| 1e-5 * work / h)
            .sum(),
    );
    assert!(
        worn < wear(&fields[0]),
        "native work hardening must reduce otherwise equal local wear"
    );
    let averaged = 1e-5 * work.iter().sum::<f64>() / (fields[4].iter().sum::<f64>() / 12.0);
    assert!(
        (worn - averaged).abs() > 1e-4 * worn,
        "averaging hardness loses the local wear response"
    );
}

#[test]
fn late_source_refusal_keeps_native_solution_and_complete_previous_hardness_field() {
    let fixture = Fixture::new();
    let problem = fixture.problem(&[0.0005, 0.0035]);
    let mechanical = with_cx(|cx| problem.solve(cx)).unwrap();
    let before = mechanical.clone();
    let mut sources = fixture.source_field();
    let previous = with_cx(|cx| {
        problem.observe_hardness(
            &mechanical,
            J2HardnessField::PerElement(&sources),
            YIELD_TOLERANCE_PA,
            cx,
        )
    })
    .unwrap();
    let alpha = mechanical.element_states[11].alpha;
    assert!(alpha > 0.0);
    let limited = ResolvedJ2Hardness::new(
        fixture.sources[1].plastic().clone(),
        &fixture.cards[1],
        fixture.sources[1].hardness_source().key().clone(),
        SELECT,
        0.5 * alpha,
    )
    .unwrap();
    sources[11] = &limited;
    assert!(matches!(
        with_cx(|cx| problem.observe_hardness(
            &mechanical,
            J2HardnessField::PerElement(&sources),
            YIELD_TOLERANCE_PA,
            cx,
        )),
        Err(J2HardnessObservationError::Element { element: 11, .. })
    ));
    sources[11] = &fixture.sources[0];
    assert!(matches!(
        with_cx(|cx| problem.observe_hardness(
            &mechanical,
            J2HardnessField::PerElement(&sources),
            YIELD_TOLERANCE_PA,
            cx,
        )),
        Err(J2HardnessObservationError::Element { element: 11, .. })
    ));
    assert!(matches!(
        with_cx(|cx| problem.observe_hardness(
            &mechanical,
            J2HardnessField::PerElement(&sources[..11]),
            YIELD_TOLERANCE_PA,
            cx,
        )),
        Err(J2HardnessObservationError::Invalid(_))
    ));
    sources[11] = &fixture.sources[1];
    assert_eq!(
        with_cx(|cx| problem.observe_hardness(
            &mechanical,
            J2HardnessField::PerElement(&sources),
            YIELD_TOLERANCE_PA,
            cx,
        ))
        .unwrap(),
        previous
    );
    assert_eq!(mechanical, before);
}
