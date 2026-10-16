//! G0/G1 analytic domain boundaries. These are synthetic module fixtures,
//! never a measured specimen or an instrument-fidelity claim.
use fs_couple::description::*;
use fs_couple::material_assembly::MaterialSource;
use fs_couple::selector::RegistrySnapshot;
use fs_duct::LossModel;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::*;
use fs_scenario::{EntityCatalog, EntityDeclaration, Environment, Scenario, WORLD};

const M: Dims = Dims([1, 0, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TAU: f64 = core::f64::consts::TAU;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn q(name: &str, value: f64, dims: Dims) -> (ParameterKey, QtyAny) {
    (ParameterKey::new(name, dims), QtyAny::new(value, dims))
}
fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", PA, 8e9),
        ("density", Dims([-3, 1, 0, 0, 0, 0]), 700.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "test-fixture:authored scalars, no specimen".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test-fixture:synthetic".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
struct Fixture {
    scenario: Scenario,
    input: AssemblyInput,
    registry: RegistrySnapshot,
    card: MaterialCard,
    point: QueryPoint,
}
impl Fixture {
    fn new(specs: Vec<(PartKind, &str, Vec<(ParameterKey, QtyAny)>)>) -> Self {
        let scenario = Scenario::new("test-fixture:validity", 9, Environment::earth_lab());
        let mut catalog = EntityCatalog::new();
        let root = catalog
            .declare(EntityDeclaration::assembly("test-fixture:root"))
            .unwrap();
        let parts = specs
            .into_iter()
            .enumerate()
            .map(|(i, (kind, generator, params))| {
                let entity = catalog
                    .declare(EntityDeclaration::part(
                        root,
                        &format!("test-fixture:part-{i}"),
                    ))
                    .unwrap();
                PartDecl {
                    entity,
                    kind,
                    geometry: GeometrySource::Parametric {
                        generator: format!("test-fixture:{generator}"),
                        params,
                    },
                    materials: vec![],
                    uncertainty: vec![],
                }
            })
            .collect();
        Self {
            scenario,
            input: AssemblyInput {
                root,
                catalog,
                parts,
                interfaces: vec![],
                configurations: vec![],
                observers: vec![],
                reference_state: Some(ReferenceState::Assembled),
            },
            registry: RegistrySnapshot::parse(include_str!("../../../instrument-claims.json"))
                .unwrap(),
            card: card(),
            point: QueryPoint::new(),
        }
    }
    fn description(&self) -> AssemblyDescription {
        AssemblyDescription::new(&self.scenario, self.input.clone(), DEFAULT_ASSEMBLY_BUDGET)
            .unwrap()
    }
    fn add_material(&mut self, index: usize) {
        self.input.parts[index].materials.push(MaterialRef {
            card: self.card.content_hash(),
            query_point: vec![],
            frame: WORLD,
            layer: None,
            unknown: vec![],
        });
    }
    fn claim(&self, band_hz: (f64, f64)) -> ClaimDomain<'_> {
        let modules = self
            .input
            .parts
            .iter()
            .map(|p| {
                let GeometrySource::Parametric { generator, .. } = &p.geometry else {
                    panic!("fixture");
                };
                let (model, image, qoi) = match generator.as_str() {
                    "test-fixture:tone-hole.v1" => (
                        DescriptionModel::CompactToneHole,
                        "tmm",
                        "fingering-peak-cents",
                    ),
                    "test-fixture:circular-opening.v1" => (
                        DescriptionModel::CompactRadiation {
                            flange: fs_phs::MouthFlange::Unflanged,
                        },
                        "phs-chain",
                        "td-oracle-agreement",
                    ),
                    "test-fixture:lumped-reed.v1" => (
                        DescriptionModel::LumpedReed,
                        "char-line",
                        "quarter-wave-lock",
                    ),
                    "test-fixture:circular-string.v1" => (
                        DescriptionModel::OscillatingCylinderString,
                        "phs-modal-bank",
                        "truncation-authority",
                    ),
                    _ => (
                        DescriptionModel::CylindricalDuct {
                            loss: LossModel::Bessel,
                        },
                        "tmm",
                        "fingering-peak-cents",
                    ),
                };
                ModuleDomain {
                    part: p.entity,
                    model,
                    claim: RegistryClaim {
                        filling: if p.kind == PartKind::String {
                            "string".into()
                        } else {
                            "wind-reed".into()
                        },
                        image: image.into(),
                        qoi: qoi.into(),
                    },
                    domain: ValidityDomain::unconstrained(),
                    materials: p
                        .materials
                        .iter()
                        .map(|_| MaterialSource {
                            card: &self.card,
                            point: &self.point,
                            selection: MaterialPropertySelection::SingleClaimOnly,
                        })
                        .collect(),
                }
            })
            .collect();
        let mut c = ClaimDomain {
            band_hz,
            dynamics: DynamicsEnvelope {
                max_mouth_pressure_pa: 20.0,
                max_bow_force_n: 0.0,
                max_pluck_force_n: 0.0,
                lowest_excitation_hz: 100.0,
                max_string_displacement_m: 0.0,
            },
            gas: gas(),
            registry: &self.registry,
            modules,
        };
        // Only declare relevant axes: unknown extra axes are deliberately
        // required by the validator, so a generic all-model box cannot pass.
        for m in &mut c.modules {
            let mut b = ValidityDomain::unconstrained().with("frequency_hz", 0.0, 50_000.0);
            match m.model {
                DescriptionModel::CylindricalDuct { .. } => {
                    b = b.with("evanescent_relative_error", 0.0, 0.01).with(
                        "length_over_shock_distance",
                        0.0,
                        0.1,
                    );
                }
                DescriptionModel::CompactToneHole => {
                    b = b.with("edge_displacement_ratio", 0.0, 0.1);
                }
                DescriptionModel::CompactRadiation { .. } => {
                    b = b.with("ka", 0.0, 1.0);
                }
                DescriptionModel::LumpedReed => {
                    b = b.with("mouth_pressure_pa", 0.0, 6000.0);
                }
                DescriptionModel::OscillatingCylinderString => {
                    b = b
                        .with("bow_force_n", 0.0, 10.0)
                        .with("pluck_force_n", 0.0, 10.0);
                }
            }
            m.domain = b;
        }
        c
    }
}
fn cylinder() -> Fixture {
    Fixture::new(vec![(
        PartKind::Duct,
        "cylinder.v1",
        vec![
            q("radius", 0.0073, M),
            q("length", 0.6, M),
            q("evanescent_relative_error", 0.001, Dims::NONE),
        ],
    )])
}
fn hole(radius: f64) -> (PartKind, &'static str, Vec<(ParameterKey, QtyAny)>) {
    (
        PartKind::Duct,
        "tone-hole.v1",
        vec![
            q("hole_radius", radius, M),
            q("bore_radius", 0.0073, M),
            q("chimney_height", 0.004, M),
            q("edge_radius", 0.001, M),
        ],
    )
}
fn reed() -> (PartKind, &'static str, Vec<(ParameterKey, QtyAny)>) {
    (
        PartKind::ValvePattern,
        "lumped-reed.v1",
        vec![
            q("rest_opening", 0.0004, M),
            q("width", 0.013, M),
            q("closing_pressure", 6000.0, PA),
            q("attack", 0.008, Dims([0, 0, 1, 0, 0, 0])),
            q("mass", 0.0, Dims([0, 1, 0, 0, 0, 0])),
            q("stiffness", 0.0, Dims([0, 1, -2, 0, 0, 0])),
            q("damping_ratio", 0.35, Dims::NONE),
        ],
    )
}
fn refused(f: &Fixture, band: (f64, f64)) -> ValidityError {
    let d = f.description();
    let c = f.claim(band);
    validate(&d, &c).unwrap_err()
}
fn set(f: &mut Fixture, index: usize, name: &str, value: f64) {
    let GeometrySource::Parametric { params, .. } = &mut f.input.parts[index].geometry else {
        panic!("fixture");
    };
    params
        .iter_mut()
        .find(|(k, _)| k.name == name)
        .unwrap()
        .1
        .value = value;
}

#[test]
fn dv_001_duct_band_above_asymmetric_cutoff_refuses() {
    let f = cylinder();
    let expected = 1.841 * gas().sound_speed / (2.0 * core::f64::consts::PI * 0.0073);
    let ValidityError::OutOfDomain { part, axis } = refused(&f, (0.0, 20_000.0)) else {
        panic!("typed cutoff refusal");
    };
    assert_eq!(part, f.input.parts[0].entity);
    assert_eq!(axis.module, "fs-duct/cylindrical");
    assert_eq!(axis.axis, "first_asymmetric_cutoff_hz");
    assert_eq!(axis.value, 20_000.0);
    assert!((axis.domain.1 / expected - 1.0).abs() <= 1e-12);
    assert!((13_700.0..13_900.0).contains(&expected));
    let d = f.description();
    let c = f.claim((0.0, 13_000.0));
    assert_eq!(validate(&d, &c).unwrap().claim_band_hz(), (0.0, 13_000.0));
}
#[test]
fn dv_002_cutoff_boundary_cap_and_cap_plus_one() {
    let f = cylinder();
    let d = f.description();
    let cutoff = 1.841 * gas().sound_speed / (TAU * 0.0073);
    let c = f.claim((0.0, cutoff));
    validate(&d, &c).unwrap();
    assert!(
        matches!(refused(&f, (0.0, cutoff * (1.0 + 1e-9))), ValidityError::OutOfDomain { axis, .. } if axis.axis == "first_asymmetric_cutoff_hz")
    );
    assert!(matches!(
        refused(&f, (0.0, f64::from_bits(cutoff.to_bits() + 1))),
        ValidityError::OutOfDomain { .. }
    ));
}
#[test]
fn dv_003_compact_radiation_ka_ceiling() {
    // Choose a = c/tau, so frequency itself is ka (no rounded inverse).
    let f = Fixture::new(vec![(
        PartKind::Duct,
        "circular-opening.v1",
        vec![q("radius", gas().sound_speed / TAU, M)],
    )]);
    let d = f.description();
    let mut c = f.claim((0.0, 1.0));
    c.dynamics.lowest_excitation_hz = 0.5;
    let a = validate(&d, &c).unwrap();
    assert_eq!(
        a.limitations()[0]
            .axes
            .iter()
            .find(|a| a.axis == "ka")
            .unwrap()
            .value,
        1.0
    );
    let mut next = c.clone();
    next.band_hz.1 = f64::from_bits(1.0f64.to_bits() + 1);
    assert!(
        matches!(validate(&d, &next), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "ka" && axis.value > 1.0 && axis.domain.1 == 1.0)
    );
    c.modules[0].domain = c.modules[0].domain.clone().with("ka", 0.0, 0.5);
    c.band_hz.1 = 0.5;
    validate(&d, &c).unwrap();
    c.band_hz.1 = 0.5 * (1.0 + 1e-9);
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "ka" && axis.domain.1 == 0.5)
    );
}
#[test]
fn dv_004_tone_hole_not_narrower_than_bore_refuses() {
    let mut f = Fixture::new(vec![hole(f64::from_bits(0.0073f64.to_bits() - 1))]);
    let d = f.description();
    let c = f.claim((0.0, 13_000.0));
    validate(&d, &c).unwrap();
    for r in [0.0073, f64::from_bits(0.0073f64.to_bits() + 1)] {
        set(&mut f, 0, "hole_radius", r);
        assert!(
            matches!(refused(&f, (0.0, 13_000.0)), ValidityError::OutOfDomain { axis, .. } if axis.axis == "hole_radius" && axis.value == r)
        );
    }
    set(&mut f, 0, "hole_radius", -1.0);
    assert!(matches!(
        refused(&f, (0.0, 13_000.0)),
        ValidityError::OwnerRefusal { .. }
    ));
}
#[test]
fn dv_005_unknown_required_property_refuses() {
    let mut f = Fixture::new(vec![reed()]);
    f.add_material(0);
    f.input.parts[0].materials[0]
        .unknown
        .push(UnknownProperty::Unknown {
            property: ParameterKey::new("young_modulus", PA),
            bound: None,
        });
    assert!(
        matches!(refused(&f, (0.0, 13_000.0)), ValidityError::UnknownRequiredProperty { part, property, .. } if part == f.input.parts[0].entity && property == "young_modulus")
    );
    // Missing actual source also refuses; Unknown cannot be replaced by a
    // scalar in geometry or a citation-only hash.
    f.input.parts[0].materials[0].unknown.clear();
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    c.modules[0].materials.clear();
    assert!(matches!(
        validate(&d, &c),
        Err(ValidityError::InvalidBinding { .. })
    ));
}
#[test]
fn dv_006_slit_channel_not_admitted_as_cylindrical_duct() {
    let mut f = cylinder();
    let GeometrySource::Parametric { generator, .. } = &mut f.input.parts[0].geometry else {
        panic!("fixture");
    };
    *generator = "test-fixture:slit.v1".into();
    assert!(
        matches!(refused(&f, (0.0, 13_000.0)), ValidityError::OutOfDomain { axis, .. } if axis.axis == "circular_cross_section" && axis.refinement == Some("M2.4-W1c / M3.2"))
    );
}
#[test]
fn dv_007_missing_axis_refuses_not_unconstrained() {
    let f = cylinder();
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    c.modules[0].domain = c.modules[0].domain.clone().with("Re", 0.0, 100.0);
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::MissingAxis { part, axis, .. }) if part == f.input.parts[0].entity && axis == "Re")
    );
    c.modules[0].domain = ValidityDomain::unconstrained();
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::MissingAxis { axis, .. }) if axis == "frequency_hz")
    );
    let mut missing = cylinder();
    let GeometrySource::Parametric { params, .. } = &mut missing.input.parts[0].geometry else {
        panic!("fixture");
    };
    params.retain(|(k, _)| k.name != "radius");
    assert!(
        matches!(refused(&missing, (0.0, 13_000.0)), ValidityError::MissingAxis { axis, .. } if axis == "radius")
    );
}
#[test]
fn dv_008_valid_fixture_admits_with_limitations() {
    let mut f = cylinder();
    // A clarinet-shaped synthetic topology: cylinder, compact hole and valve.
    // Rebuild together so EntityIds and catalog ownership remain exact.
    let GeometrySource::Parametric { params, .. } = f.input.parts[0].geometry.clone() else {
        panic!("fixture");
    };
    f = Fixture::new(vec![
        (PartKind::Duct, "cylinder.v1", params),
        hole(0.003),
        reed(),
    ]);
    f.add_material(2);
    let d = f.description();
    let c = f.claim((0.0, 13_000.0));
    let a = validate(&d, &c).unwrap();
    let b = validate(&d, &c).unwrap();
    assert!(std::ptr::eq(a.description(), &d));
    assert!(std::ptr::eq(a.claim(), &c));
    assert_eq!(a.claim_band_hz(), (0.0, 13_000.0));
    assert_eq!(a.limitations(), b.limitations());
    assert_eq!(a.limitations().len(), 3);
    assert!(a.limitations().iter().any(|l| l.statement.contains("1-D")));
    assert!(
        a.limitations()
            .iter()
            .any(|l| l.statement.contains("Compact tone-hole"))
    );
    assert!(
        a.limitations()
            .iter()
            .any(|l| l.statement.contains("Lumped reed"))
    );
    assert_eq!(a.limitations()[2].material_states.len(), 1);
    assert_eq!(a.limitations()[0].claim.image, "tmm");
    // An absent claim row cannot turn an otherwise valid model into Admission.
    let mut bad = c.clone();
    bad.modules[0].claim.qoi = "not-registered".into();
    assert!(matches!(
        validate(&d, &bad),
        Err(ValidityError::RegistryRefusal { .. })
    ));
}

#[test]
fn dv_009_dynamics_drive_edge_shock_and_string_drag_domains() {
    let f = Fixture::new(vec![hole(0.003)]);
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    let a = validate(&d, &c).unwrap();
    let ratio = a.limitations()[0]
        .axes
        .iter()
        .find(|a| a.axis == "edge_displacement_ratio")
        .unwrap()
        .value;
    let oracle = 20.0 / (gas().density * gas().sound_speed * TAU * 100.0 * 0.001);
    assert!((ratio / oracle - 1.0).abs() < 1e-12);
    c.modules[0].domain = c.modules[0]
        .domain
        .clone()
        .with("edge_displacement_ratio", 0.0, ratio);
    validate(&d, &c).unwrap();
    c.dynamics.max_mouth_pressure_pa *= 1.0 + 1e-9;
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "edge_displacement_ratio")
    );

    let f = cylinder();
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    let a = validate(&d, &c).unwrap();
    let ratio = a.limitations()[0]
        .axes
        .iter()
        .find(|a| a.axis == "length_over_shock_distance")
        .unwrap()
        .value;
    let g = gas();
    let oracle = 0.6 * ((g.gamma + 1.0) / 2.0) * (TAU * 13_000.0) * 20.0
        / (g.density * g.sound_speed.powi(3));
    assert!((ratio / oracle - 1.0).abs() < 1e-12);
    c.modules[0].domain =
        c.modules[0]
            .domain
            .clone()
            .with("length_over_shock_distance", 0.0, ratio);
    validate(&d, &c).unwrap();
    c.dynamics.max_mouth_pressure_pa *= 1.0 + 1e-9;
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "length_over_shock_distance")
    );

    let mut f = Fixture::new(vec![(
        PartKind::String,
        "circular-string.v1",
        vec![q("radius", 0.0005, M)],
    )]);
    f.add_material(0);
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    c.dynamics.max_string_displacement_m = 0.001 / TAU;
    validate(&d, &c).unwrap();
    c.dynamics.max_string_displacement_m *= 1.0 + 1e-9;
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "keulegan_carpenter")
    );
    // A stricter caller policy intersects, and cannot be skipped merely
    // because the owner has already evaluated this axis.
    c.dynamics.max_string_displacement_m = 0.00075 / TAU;
    c.modules[0].domain = c.modules[0]
        .domain
        .clone()
        .with("keulegan_carpenter", 0.0, 0.5);
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::OutOfDomain { axis, .. }) if axis.axis == "keulegan_carpenter" && axis.domain.1 == 0.5)
    );
}

#[test]
fn dv_010_configurations_and_binding_fail_closed() {
    let mut f = cylinder();
    let part = f.input.parts[0].entity;
    f.input.configurations.push(ConfigurationState {
        name: "wider".into(),
        assignments: vec![(
            ParameterPath {
                entity: part,
                parameter: ParameterKey::new("radius", M),
            },
            QtyAny::new(0.03, M),
        )],
        provenance: StateProvenance::Declared {
            rationale: "test-fixture:configuration changes transverse cutoff".into(),
        },
    });
    assert!(
        matches!(refused(&f, (0.0, 13_000.0)), ValidityError::OutOfDomain { axis, .. } if axis.axis == "first_asymmetric_cutoff_hz")
    );
    f.input.configurations[0].assignments[0].1.value = 0.0073;
    let d = f.description();
    let mut c = f.claim((0.0, 13_000.0));
    let a = validate(&d, &c).unwrap();
    assert_eq!(a.limitations().len(), 2);
    assert_eq!(a.limitations()[1].configuration.as_deref(), Some("wider"));
    c.modules.push(c.modules[0].clone());
    assert!(matches!(
        validate(&d, &c),
        Err(ValidityError::InvalidBinding { .. })
    ));
    c.modules.clear();
    assert!(matches!(
        validate(&d, &c),
        Err(ValidityError::UnsupportedPart { .. })
    ));
    let mut c = f.claim((0.0, 13_000.0));
    c.modules[0].domain = c.modules[0].domain.clone().with_quantity(
        "frequency_hz",
        fs_qty::QuantitySpec::dimensional(M),
        0.0,
        50_000.0,
    );
    assert!(
        matches!(validate(&d, &c), Err(ValidityError::DimensionMismatch { axis, .. }) if axis == "frequency_hz")
    );
    let mut c = f.claim((0.0, 13_000.0));
    c.gas.characteristic_impedance *= 2.0;
    assert!(matches!(
        validate(&d, &c),
        Err(ValidityError::InvalidBinding {
            module: "fs-material/gas",
            ..
        })
    ));
}
