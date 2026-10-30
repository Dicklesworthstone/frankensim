//! Sourced expansion drives native Kirchhoff storage, Joule heat and field motion.
//! Synthetic material tables establish numerical causality, not material fidelity.

#[path = "expansion/boundary.rs"]
mod boundary;

use super::*;
use fs_couple::circuit_enthalpy::CircuitEnthalpyError;
use fs_couple::material_circuit::MaterialCircuit;
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue,
    Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::field_storage::{
    MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

const CELL: CellId = CellId(17);
const R0: f64 = 2.0;
const CP: f64 = 2.0;

fn temperature(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}

fn claim(key: PropertyKey, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            dims: key.quantity().dims(),
            knots,
        },
        key,
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            290.0,
            380.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic native expansion fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn material_expansion(
    basis: ElectricalPropertyBasis,
    alphas: [f64; 3],
) -> (MaterialCard, ResolvedHeatCapacityCurve, [PropertyKey; 3]) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new(basis.property_name(), basis.dimensions());
    let keys = std::array::from_fn(|i| PropertyKey::new(format!("length_change_{i}"), Dims::NONE));
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (heat.clone(), CP),
        (density.clone(), 1.0),
        (
            electrical,
            if basis == ElectricalPropertyBasis::Resistivity {
                R0
            } else {
                1.0 / R0
            },
        ),
    ] {
        claims
            .insert_claim(claim(key, vec![(290.0, value), (380.0, value)]))
            .unwrap();
    }
    for i in 0..3 {
        // A nonzero fitted datum must cancel through the common source owner.
        claims
            .insert_claim(claim(
                keys[i].clone(),
                vec![
                    (290.0, 0.2 - 12.0 * alphas[i]),
                    (380.0, 0.2 + 96.0 * alphas[i]),
                ],
            ))
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "native source coupling fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &temperature(290.0),
            upper: &temperature(380.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric, keys)
}

fn path(card: &MaterialCard, keys: &[PropertyKey; 3], high: f64) -> PreparedFreeExpansion {
    PreparedFreeExpansion::try_new(
        card,
        FreeExpansionGeometryQuery {
            reference: &temperature(300.0),
            temperature_axis: "T",
            lower: &temperature(290.0),
            upper: &temperature(high),
            strains: FreeExpansionStrainSources::FixedPrincipal {
                frame: fs_blake3::hash_bytes(b"native conductor principal axes"),
                strains: [&keys[0], &keys[1], &keys[2]],
                selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
            },
        },
    )
    .unwrap()
}

fn thermal(curve: &ResolvedHeatCapacityCurve, htc: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: 10.0 * CP,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    htc_w_m2_k: htc,
                    ..boundary()
                },
            }),
        }],
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn cell<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: CELL,
        material_id: MaterialId(0),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }
}
fn resistor(branch: usize) -> CellResistor {
    CellResistor {
        destination: fs_couple::circuit_enthalpy::ResistorCell { branch, cell: CELL },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}
fn resistance(alphas: [f64; 3], axis: usize, t: f64) -> f64 {
    let s = alphas.map(|a| 1.0 + a * (t - 300.0));
    R0 * s[axis] / (s[(axis + 1) % 3] * s[(axis + 2) % 3])
}
fn root(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lo) <= 0.0 && f(hi) >= 0.0);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_native_voltage_and_current_heating_use_actual_principal_expansion() {
    let alphas = [0.015, 0.005, 0.02];
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for axis in 0..3 {
            for voltage_drive in [false, true] {
                let (card, caloric, keys) = material_expansion(basis, alphas);
                let context = temperature(300.0);
                let cells = [cell(&card, &caloric, &context)];
                let resistors = [resistor(0)];
                let network = MaterialNetwork {
                    cells: &cells,
                    resistors: &resistors,
                    links: &[],
                };
                let mut heat = thermal(&caloric, 0.0);
                let reference = heat.cells()[&CELL].geometry;
                let expanded = network
                    .prepare_with_bases(&heat, &[(0, basis)], || false)
                    .unwrap()
                    .with_free_expansion(
                        &heat,
                        &[ConductorFreeExpansion {
                            branch: 0,
                            source: path(&card, &keys, 380.0),
                            longitudinal_axis: axis,
                        }],
                        || false,
                    )
                    .unwrap();
                let drive = if voltage_drive {
                    Branch::VoltageSource { port: 0 }
                } else {
                    Branch::CurrentSource { port: 0 }
                };
                let amplitude: f64 = if voltage_drive { 4.0 } else { 2.0 };
                let mut circuit = MaterialCircuit::new(
                    CircuitGraph {
                        node_count: 2,
                        branches: vec![(1, 0, Branch::Resistor { ohms: f64::NAN }), (1, 0, drive)],
                        transformers: vec![],
                    },
                    expanded.resolve(&heat, || false).unwrap().resistors,
                    vec![],
                    vec![],
                    &[amplitude],
                )
                .unwrap();
                for step in 1..=2 {
                    let before = heat.temperature_k(CELL).unwrap();
                    let power = |t| {
                        if voltage_drive {
                            amplitude.powi(2) / resistance(alphas, axis, t)
                        } else {
                            amplitude.powi(2) * resistance(alphas, axis, t)
                        }
                    };
                    let expected = root(before, 379.0, |t| CP * (t - before) - 0.25 * power(t));
                    let next = expanded
                        .step_candidate(
                            &heat,
                            f64::from(step) * 0.25,
                            &circuit,
                            &[amplitude],
                            config(),
                            &[0.0],
                            || false,
                        )
                        .unwrap();
                    let actual = next.thermal().temperature_k(CELL).unwrap();
                    assert!((actual - expected).abs() < 1e-7, "{actual} vs {expected}");
                    let reported = &next.report().electrical.resistors[0].material.conductor;
                    near(reported.resistance_ohm(), resistance(alphas, axis, actual));
                    near(reported.length_m(), 1.0 + alphas[axis] * (actual - 300.0));
                    near(
                        reported.area_m2(),
                        (1.0 + alphas[(axis + 1) % 3] * (actual - 300.0))
                            * (1.0 + alphas[(axis + 2) % 3] * (actual - 300.0)),
                    );
                    assert_eq!(reported.source_basis(), basis);
                    assert_eq!(reported.material().card_identity(), card.content_hash());
                    near(
                        next.report().electrical.native.supplied,
                        CP * (actual - before),
                    );
                    near(
                        next.report().deposited_joule_j,
                        next.report().electrical.heat_j,
                    );
                    assert_eq!(next.thermal().cells()[&CELL].mass_kg, 1.0);
                    assert_eq!(next.thermal().cells()[&CELL].geometry, reference);
                    let (new_heat, new_circuit, _) = next.into_parts();
                    heat = new_heat;
                    circuit = new_circuit;
                }
            }
        }
    }
}

#[test]
fn g1_expansion_changes_finite_inductor_current_and_closes_native_energy() {
    let run = |alpha| {
        let basis = ElectricalPropertyBasis::Conductivity;
        let alphas = [alpha; 3];
        let (card, caloric, keys) = material_expansion(basis, alphas);
        let context = temperature(300.0);
        let cells = [cell(&card, &caloric, &context)];
        let resistors = [resistor(1)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut heat = thermal(&caloric, 0.0);
        let expanded = network
            .prepare_with_bases(&heat, &[(1, basis)], || false)
            .unwrap()
            .with_free_expansion(
                &heat,
                &[ConductorFreeExpansion {
                    branch: 1,
                    source: path(&card, &keys, 380.0),
                    longitudinal_axis: 0,
                }],
                || false,
            )
            .unwrap();
        let l = 0.8;
        let dt = 0.2;
        let mut circuit = MaterialCircuit::new(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Inductor { henries: l }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            expanded.resolve(&heat, || false).unwrap().resistors,
            vec![4.0],
            vec![],
            &[],
        )
        .unwrap();
        let energy = heat.energy_j() + 4.0_f64.powi(2) / (2.0 * l);
        for step in 1..=4 {
            let old_t = heat.temperature_k(CELL).unwrap();
            let old_phi = circuit.fluxes_wb()[0];
            let phi = |t| {
                let lambda = dt * resistance(alphas, 0, t) / (2.0 * l);
                old_phi * (1.0 - lambda) / (1.0 + lambda)
            };
            let expected_t = root(old_t, 379.0, |t| {
                CP * (t - old_t) - (old_phi * old_phi - phi(t) * phi(t)) / (2.0 * l)
            });
            let next = expanded
                .step_candidate(
                    &heat,
                    f64::from(step) * dt,
                    &circuit,
                    &[],
                    config(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let t = next.thermal().temperature_k(CELL).unwrap();
            assert!((t - expected_t).abs() < 1e-7);
            near(next.circuit().fluxes_wb()[0], phi(t));
            near(
                next.thermal().energy_j() + next.report().electrical.stored_after_j,
                energy,
            );
            near(next.report().total_energy_residual_j, 0.0);
            let (a, b, _) = next.into_parts();
            heat = a;
            circuit = b;
        }
        (circuit.fluxes_wb()[0], heat.temperature_k(CELL).unwrap())
    };
    let constant = run(0.0);
    let expanding = run(0.02);
    assert!(
        expanding.0 > constant.0 && expanding.1 < constant.1,
        "lower expansion-dependent R must retain more finite magnetic storage"
    );
}

fn scalar_card(name: &str, dims: Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut source = claim(
        PropertyKey::new(name, dims),
        vec![(290.0, value), (380.0, value)],
    );
    source.value = PropertyValue::Scalar { value, dims };
    source.validity = ValidityDomain::unconstrained();
    source.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claims.insert_claim(source).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "synthetic solid".into(),
            process: "mechanical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn body() -> ResolvedBody {
    let card = scalar_card("density", Density::DIMS, 100.0);
    resolve_body(
        &[RegionBinding {
            key: "guided magnetic pole",
            material: &card,
            point: &QueryPoint::new(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        }],
        1,
        || false,
    )
    .unwrap()
}

#[test]
fn g3_source_only_expansion_reaches_reciprocal_magnetic_motion_and_heat() {
    let run = |alpha| {
        let basis = ElectricalPropertyBasis::Resistivity;
        let (card, caloric, keys) = material_expansion(basis, [alpha; 3]);
        let context = temperature(300.0);
        let cells = [cell(&card, &caloric, &context)];
        let resistors = [resistor(1)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut heat = thermal(&caloric, 0.0);
        let expanded = network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion(
                &heat,
                &[ConductorFreeExpansion {
                    branch: 1,
                    source: path(&card, &keys, 380.0),
                    longitudinal_axis: 0,
                }],
                || false,
            )
            .unwrap();
        let magnetic = scalar_card(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        );
        let coil = MaterialStorage {
            branch: 0,
            element: MaterialStorageElement::Inductor(
                ResolvedInductor::resolve(
                    &magnetic,
                    &QueryPoint::new(),
                    MaterialPropertySelection::SingleClaimOnly,
                    1.0,
                    1.0,
                    1,
                )
                .unwrap(),
            ),
        };
        let circuit = MaterialCircuit::new_with_storage(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Inductor { henries: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            expanded.resolve(&heat, || false).unwrap().resistors,
            vec![coil],
            vec![4.0],
            vec![],
            &[],
        )
        .unwrap();
        let mut moving = MaterialFieldMotion::new(circuit, 0, body(), 0.0, [0.5, 1.5]).unwrap();
        let stored = |motion: &MaterialFieldMotion| {
            motion
                .inductor()
                .unwrap()
                .stored_energy_j(motion.circuit().fluxes_wb()[0])
                .unwrap()
                + 0.5 * motion.body().inertia().mass() * motion.velocity_m_s().powi(2)
        };
        let total = heat.energy_j() + stored(&moving);
        for step in 1..=6 {
            let next = expanded
                .field_motion_candidate(
                    &heat,
                    f64::from(step) * 0.1,
                    &moving,
                    &[],
                    0.0,
                    config(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            near(next.thermal().energy_j() + stored(next.motion()), total);
            near(
                next.report().deposited_joule_j,
                next.report().electrical.heat_j,
            );
            let (t, m, _) = next.into_parts();
            heat = t;
            moving = m;
        }
        (
            moving.gap_m(),
            moving.circuit().fluxes_wb()[0],
            heat.temperature_k(CELL).unwrap(),
        )
    };
    let fixed = run(0.0);
    let expanding = run(0.02);
    assert!(
        expanding.0 < fixed.0 && expanding.1 > fixed.1 && expanding.2 < fixed.2,
        "free resistor expansion must change actual electromagnetic attraction via native current"
    );
}

#[test]
fn g4_source_boundary_and_late_cancellation_refuse_without_advancing_histories() {
    let (card, caloric, keys) = material_expansion(ElectricalPropertyBasis::Resistivity, [0.02; 3]);
    let context = temperature(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor(0)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let heat = thermal(&caloric, 0.0);
    let binding = ConductorFreeExpansion {
        branch: 0,
        source: path(&card, &keys, 380.0),
        longitudinal_axis: 0,
    };
    let expanded = network
        .prepare(&heat, || false)
        .unwrap()
        .with_free_expansion(&heat, &[binding.clone()], || false)
        .unwrap();
    let circuit = MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        expanded.resolve(&heat, || false).unwrap().resistors,
        vec![],
        vec![],
        &[4.0],
    )
    .unwrap();
    let before = heat.cells().clone();
    let source_before = circuit.resistors().to_vec();
    let mut polls = 0;
    let expected = expanded
        .step_candidate(&heat, 0.25, &circuit, &[4.0], config(), &[0.0], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut actual = 0;
    assert!(matches!(
        expanded.step_candidate(&heat, 0.25, &circuit, &[4.0], config(), &[0.0], || {
            actual += 1;
            actual == polls
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(heat.cells(), &before);
    assert_eq!(circuit.resistors(), &source_before);
    assert_eq!(heat.time_s(), 0.0);
    assert!(circuit.charges_c().is_empty() && circuit.fluxes_wb().is_empty());
    let replay = expanded
        .step_candidate(&heat, 0.25, &circuit, &[4.0], config(), &[0.0], || false)
        .unwrap();
    assert_eq!(replay.thermal().cells(), expected.thermal().cells());
    assert_eq!(replay.circuit().resistors(), expected.circuit().resistors());
    let narrow = network
        .prepare(&heat, || false)
        .unwrap()
        .with_free_expansion(
            &heat,
            &[ConductorFreeExpansion {
                source: path(&card, &keys, 300.1),
                ..binding.clone()
            }],
            || false,
        )
        .unwrap();
    assert!(
        narrow
            .step_candidate(&heat, 0.25, &circuit, &[4.0], config(), &[0.0], || false)
            .is_err()
    );
    assert!(
        network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion(
                &heat,
                &[ConductorFreeExpansion {
                    longitudinal_axis: 3,
                    ..binding.clone()
                }],
                || false
            )
            .is_err()
    );
    let cooling = thermal(&caloric, 1.0);
    assert!(
        network
            .prepare(&cooling, || false)
            .unwrap()
            .with_free_expansion(&cooling, &[binding.clone()], || false)
            .is_err()
    );
    let (foreign, _, foreign_keys) =
        material_expansion(ElectricalPropertyBasis::Resistivity, [0.01; 3]);
    assert!(
        network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion(
                &heat,
                &[ConductorFreeExpansion {
                    source: path(&foreign, &foreign_keys, 380.0),
                    ..binding
                }],
                || false
            )
            .is_err()
    );
    assert_eq!(heat.cells(), &before);
    assert_eq!(circuit.resistors(), &source_before);
}
