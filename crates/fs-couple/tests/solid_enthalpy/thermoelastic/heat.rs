//! G1/G3/G4: original thermal material data changes loaded solid motion.
//! The two-variable thermal/Newmark reference uses explicit analytic equations,
//! not the runtime's source queries, radiation solve or heat fixed point.
use super::*;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, GrayDiffuseEnclosure, RadiationSurface, SURFACE_EMISSIVITY_PROPERTY,
    SurfaceEmissivity, ViewFactorMatrix,
    fixtures::{box_grid, on_box_face},
    mesh::ConductionMesh,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::*,
};
use fs_couple::solid_enthalpy::thermoelastic::heat::ThermoelasticHeatError;
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SelectionPolicy, SurfaceSpec, SystemContext};
use fs_material::{
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::{
        conduction::{LayerConductionConfig, network::*, source::ResolvedConductionLayer},
        interface::ThermalInterfacePoint,
        source::ThermalMaterialPoint,
    },
};
use fs_solid::loads3::transient::TransientConfig3;

const YOUNG: [f64; 2] = [1000.0, 1600.0];
const SLOPE: [f64; 2] = [2.0, -1.0];
const ALPHA: [f64; 2] = [0.0002, 0.0003];
const CAPACITY: [f64; 2] = [0.2, 0.3];
const POWER: [f64; 2] = [2.0, 0.1];
const BETA: f64 = 0.002;

fn transport_card(i: usize, emissivity: f64, conductivity: f64) -> MaterialCard {
    let original = fixture(YOUNG[i], SLOPE[i], ALPHA[i], CAPACITY[i], false).0;
    let mut claims = original.claims().clone();
    for (name, dims, value, slope) in [
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            conductivity,
            BETA * conductivity,
        ),
        (
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            emissivity,
            0.0005,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, value - 20.0 * slope),
                        (500.0, value + 200.0 * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    500.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coupled solid transport".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn with_network<R>(
    emissivity: f64,
    conductivity_scale: f64,
    resistance: f64,
    integrated: bool,
    heat_pressure: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &DepositionSimulation<'_>,
        &SolidEnthalpy<'_>,
        &PreparedHeatNetwork<'_, '_>,
    ) -> R,
) -> R {
    with_cx(|cx| {
        let cards = [
            transport_card(0, emissivity, 0.4 * conductivity_scale),
            transport_card(1, 0.6, 0.8 * conductivity_scale),
        ];
        let curves = [0, 1].map(|i| {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                cards[i].content_hash(),
                SolidLiquidPhase::Solid,
                [280.0, 500.0]
                    .map(|t| EnthalpyPhaseKnot {
                        temperature_k: t,
                        specific_enthalpy_j_kg: CAPACITY[i] * (t - 280.0),
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: 6.0,
                    })
                    .to_vec(),
            )
            .unwrap()
        });
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let context = context(300.0)
            .with_quantity("p", pressure, 100000.0)
            .unwrap();
        let heat_context = context
            .clone()
            .with_quantity("p", pressure, heat_pressure)
            .unwrap();
        let thermal = thermal(&[&curves[0], &curves[1]], &CAPACITY, &POWER, 0.01);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                source(30, &cards[0], &curves[0], &context),
                source(10, &cards[1], &curves[1], &context),
            ],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        let points = [0, 1].map(|i| ThermalMaterialPoint {
            material: &cards[i],
            context: &heat_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        });
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [0, 1].map(|i| MaterialCell {
            cell: IDS[i],
            material_id: MaterialId(i as u64),
            curve: &curves[i],
            point: points[i],
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 9,
            distance_a_m: 0.5,
            distance_b_m: 0.5,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(RP, RD),
                value: PropertyValue::Curve {
                    abscissa: "aT".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 0.98 * resistance), (500.0, 1.2 * resistance)],
                    dims: RD,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("aT", quantity(), 280.0, 500.0)
                    .with_quantity("bT", quantity(), 280.0, 500.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored thermal interface".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        let interface = InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: cards[0].id().clone(),
                texture_frame: "a".into(),
            },
            SurfaceSpec {
                material: cards[1].id().clone(),
                texture_frame: "b".into(),
            },
            SystemContext {
                medium: "fixed solid contact".into(),
                third_body: None,
                environment: "authored".into(),
                history: "fixed".into(),
            },
            claims,
            vec![],
        )
        .unwrap();
        let interface_context = QueryPoint::new()
            .with_quantity("aT", quantity(), 300.0)
            .unwrap()
            .with_quantity("bT", quantity(), 300.0)
            .unwrap();
        let contacts = if integrated {
            vec![]
        } else {
            vec![MaterialContactLink {
                link: 9,
                point: ThermalInterfacePoint {
                    interface: &interface,
                    materials: [&cards[0], &cards[1]],
                    context: &interface_context,
                    temperature_axes: ["aT", "bT"],
                    liquid_fraction_axes: [None, None],
                    mean_temperature_axis: None,
                },
                sides: [0, 1].map(|i| ContactCell {
                    cell: IDS[i],
                    material_id: MaterialId(i as u64),
                    curve: &curves[i],
                }),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            }]
        };
        let (complex, positions) = box_grid([1, 1, 1], [1.0; 3]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let surfaces = [0, 1].map(|i| {
            let phase = curves[i]
                .state_at_specific_enthalpy(20.0 * CAPACITY[i])
                .unwrap();
            let emissivity = SurfaceEmissivity::from_thermal_state(
                &format!("solid-{i}"),
                points[i],
                phase,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
            .0;
            RadiationSurface::new(
                &mesh,
                format!("surface-{i}"),
                |face| on_box_face(face.centroid[0], i as f64),
                emissivity,
            )
            .unwrap()
        });
        // Infinite-plane view-factor reference; cube faces supply area/trace only.
        let enclosure = GrayDiffuseEnclosure::new(
            surfaces.to_vec(),
            ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        )
        .unwrap();
        let radiating = [0, 1].map(|i| RadiatingCell {
            surface: i,
            cell: IDS[i],
            material_id: MaterialId(i as u64),
            curve: &curves[i],
            material: points[i],
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
        let radiation = MaterialRadiation {
            enclosure: &enclosure,
            cells: &radiating,
            power_tolerance_w: 1e-10,
            radiosity_tolerance_w_m2: 1e-10,
        };
        let network = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &contacts,
            radiation: Some(&radiation),
        };
        let layers = [0, 1].map(|i| {
            ResolvedConductionLayer::new(
                &curves[i],
                points[i],
                &key,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                SolidLiquidPhase::Solid,
                [280.0, 360.0],
                || Ok(()),
            )
            .unwrap()
        });
        let integrated_layers = PreparedLayerConduction::new(
            &thermal,
            vec![LayerLinkBinding {
                link: 9,
                layers: [&layers[0], &layers[1]],
                distances_m: [0.5; 2],
                contact_area_resistance_m2_k_w: 0.0,
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        let prepared = network.prepare(cx, &thermal).unwrap();
        let prepared = if integrated {
            prepared
                .with_integrated_conduction(
                    cx,
                    &thermal,
                    &integrated_layers,
                    LayerConductionConfig {
                        max_iterations: 80,
                        absolute_heat_tolerance_w: 1e-11,
                        relative_heat_tolerance: 1e-12,
                    },
                )
                .unwrap()
        } else {
            prepared
        };
        run(cx, &thermal, &source, &prepared)
    })
}

fn flux(t: [f64; 2], emissivity: f64, scale: f64, resistance: f64, integrated: bool) -> f64 {
    let k = [0.4 * scale, 0.8 * scale];
    let conduction = if integrated {
        let primitive = t.map(|x| (x - 300.0) + 0.5 * BETA * (x - 300.0).powi(2));
        0.5 * (primitive[0] - primitive[1]) / (0.5 / k[0] + 0.5 / k[1])
    } else {
        let actual = [0, 1].map(|i| k[i] * (1.0 + BETA * (t[i] - 300.0)));
        0.5 * (t[0] - t[1])
            / (0.5 / actual[0] + 0.5 / actual[1] + resistance * (1.0 + 0.001 * (t[0] - 300.0)))
    };
    let e = [
        emissivity + 0.0005 * (t[0] - 300.0),
        0.6 + 0.0005 * (t[1] - 300.0),
    ];
    let radiation =
        5.670374419e-8 * (t[0].powi(4) - t[1].powi(4)) / (1.0 / e[0] + 1.0 / e[1] - 1.0);
    conduction + radiation
}

#[allow(clippy::too_many_arguments)]
fn reference(
    old_t: [f64; 2],
    old_u: f64,
    old_v: f64,
    forces: [f64; 2],
    emissivity: f64,
    scale: f64,
    resistance: f64,
    integrated: bool,
    inertial: bool,
) -> ([f64; 2], f64, f64, [f64; 2]) {
    let dt: f64 = 0.1;
    let evaluate = |t: [f64; 2]| {
        let u = if inertial {
            let old_k = [0, 1].map(|i| VOLUME * (YOUNG[i] + SLOPE[i] * (old_t[i] - 300.0)));
            let k = [0, 1].map(|i| VOLUME * (YOUNG[i] + SLOPE[i] * (t[i] - 300.0)));
            let thermal_force = (0..2)
                .map(|i| ALPHA[i] * (old_k[i] * (old_t[i] - 300.0) + k[i] * (t[i] - 300.0)))
                .sum::<f64>();
            // Two unit-mass tetrahedra: retained consistent diagonal mass = 0.2 kg.
            ((0.4 / dt.powi(2) - 0.5 * old_k.iter().sum::<f64>()) * old_u
                + 0.4 / dt * old_v
                + 0.5 * (forces[0] + forces[1] + thermal_force))
                / (0.4 / dt.powi(2) + 0.5 * k.iter().sum::<f64>())
        } else {
            displacement(t, forces[1], YOUNG, SLOPE, ALPHA)
        };
        let v = if inertial {
            2.0 * (u - old_u) / dt - old_v
        } else {
            0.0
        };
        let heat = [0, 1].map(|i| {
            let a = element(old_t[i], old_u, YOUNG[i], SLOPE[i], ALPHA[i]);
            let b = element(t[i], u, YOUNG[i], SLOPE[i], ALPHA[i]);
            0.5 * VOLUME * (a.0 + b.0) * (u - old_u) - (b.1 - a.1)
        });
        let mut predicted = old_t;
        for _ in 0..4 {
            let q = flux(predicted, emissivity, scale, resistance, integrated);
            predicted = [0, 1].map(|i| {
                predicted[i]
                    + 0.025 / CAPACITY[i] * (POWER[i] + heat[i] / dt + if i == 0 { -q } else { q })
            });
        }
        (
            [0, 1].map(|i| CAPACITY[i] * (t[i] - predicted[i])),
            u,
            v,
            heat,
        )
    };
    let mut t = old_t;
    for _ in 0..40 {
        let (r, u, v, heat) = evaluate(t);
        if r[0].abs().max(r[1].abs()) < 1e-12 {
            return (t, u, v, heat);
        }
        let mut jac = [[0.0; 2]; 2];
        for col in 0..2 {
            let mut probe = t;
            probe[col] += 1e-4;
            let next = evaluate(probe).0;
            for row in 0..2 {
                jac[row][col] = (next[row] - r[row]) / 1e-4;
            }
        }
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        t[0] -= (jac[1][1] * r[0] - jac[0][1] * r[1]) / det;
        t[1] -= (jac[0][0] * r[1] - jac[1][0] * r[0]) / det;
    }
    panic!("independent coupled thermal/motion reference did not converge");
}

fn run_case(e: f64, k: f64, r: f64, integrated: bool, inertial: bool) -> ([f64; 2], f64) {
    with_network(
        e,
        k,
        r,
        integrated,
        100000.0,
        |cx, accepted, source, prepared| {
            let fixed: Vec<_> = (0..15).filter(|d| *d != 3).collect();
            let initial = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
                .unwrap();
            let body = ThermoelasticSolid::try_new(
                source,
                accepted,
                &initial,
                &fixed,
                &[[280.0, 360.0]; 2],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let mut stationary = Some(body);
            let mut dynamic = if inertial {
                let mut velocity = [[0.0; 3]; 5];
                velocity[1][0] = 0.02;
                Some(
                    stationary
                        .take()
                        .unwrap()
                        .begin_transient(
                            &initial,
                            &velocity,
                            TransientConfig3 {
                                maximum_steps: 8,
                                maximum_displacement_m: 0.05,
                                absolute_force_tolerance_n: 1e-11,
                                absolute_energy_tolerance_j: 1e-10,
                                relative_energy_tolerance: 1e-8,
                                ..TransientConfig3::default()
                            },
                            cx,
                        )
                        .unwrap(),
                )
            } else {
                None
            };
            let mut thermal = accepted.clone();
            let mut record = Steps {
                count: 0,
                cancel_at: usize::MAX,
            };
            let mut old_t = [300.0; 2];
            let mut old_u = displacement(old_t, 0.8, YOUNG, SLOPE, ALPHA);
            let mut old_v = if inertial { 0.02 } else { 0.0 };
            let mut old_force = 0.8;
            for (n, force) in [1.0, 0.6, 0.9].into_iter().enumerate() {
                let expected = reference(
                    old_t,
                    old_u,
                    old_v,
                    [old_force, force],
                    e,
                    k,
                    r,
                    integrated,
                    inertial,
                );
                let loads = source
                    .geometry()
                    .assemble(&[], VolumeLoad3::None, &[(3, force)], cx)
                    .unwrap();
                let target = (n + 1) as f64 * 0.1;
                let (next_thermal, next_record, temperatures, u, v, work, heat, audit) = if inertial
                {
                    let (next, audit) = dynamic
                        .as_ref()
                        .unwrap()
                        .advance_heat_network_candidate(
                            &thermal,
                            prepared,
                            &record,
                            target,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    let motion = next.solid.nodal_motion(cx).unwrap();
                    assert_eq!(next.solid.step(), n + 1);
                    let t = [0, 1].map(|i| next.thermal.temperature_k(IDS[i]).unwrap());
                    let work = next.report.mechanical_work_j;
                    let heat = next
                        .report
                        .cells
                        .iter()
                        .map(|x| x.reversible_heat_j)
                        .collect::<Vec<_>>();
                    near(
                        next.report.kinetic_energy_change_j,
                        0.1 * (motion[1][1][0].powi(2) - old_v.powi(2)),
                        2e-9,
                    );
                    let (next_thermal, body, observers, _) = next.into_parts();
                    dynamic = Some(body);
                    (
                        next_thermal,
                        observers.second,
                        t,
                        motion[0][1][0],
                        motion[1][1][0],
                        work,
                        heat,
                        audit,
                    )
                } else {
                    let (next, audit) = stationary
                        .as_ref()
                        .unwrap()
                        .advance_heat_network_candidate(
                            &thermal,
                            prepared,
                            &record,
                            target,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    let u = next
                        .solid
                        .equilibrium()
                        .solution
                        .displacement
                        .displacement_m()[1][0];
                    let t = [0, 1].map(|i| next.thermal.temperature_k(IDS[i]).unwrap());
                    let work = next.report.mechanical_work_j;
                    let heat = next
                        .report
                        .cells
                        .iter()
                        .map(|x| x.reversible_heat_j)
                        .collect::<Vec<_>>();
                    let (next_thermal, body, observers, _) = next.into_parts();
                    stationary = Some(body);
                    (next_thermal, observers.second, t, u, 0.0, work, heat, audit)
                };
                for i in 0..2 {
                    near(temperatures[i], expected.0[i], 4e-8);
                    near(heat[i], expected.3[i], 2e-9);
                }
                near(u, expected.1, 2e-10);
                near(v, expected.2, 2e-9);
                near(work, 0.5 * (old_force + force) * (u - old_u), 2e-10);
                near(audit.external_heat_j, 0.1 * POWER.iter().sum::<f64>(), 1e-9);
                near(audit.thermal_residual_j, 0.0, solve().energy_tolerance_j);
                near(audit.combined_residual_j, 0.0, solve().energy_tolerance_j);
                assert!(next_record.count >= (n + 1) * 4);
                thermal = next_thermal;
                record = next_record;
                old_t = temperatures;
                old_u = u;
                old_v = v;
                old_force = force;
            }
            assert_eq!(accepted.time_s(), 0.0);
            (old_t, old_u)
        },
    )
}

#[test]
fn g1_material_transport_changes_reciprocal_equilibrium_and_continuing_inertia() {
    for inertial in [false, true] {
        let baseline = run_case(0.2, 1.0, 0.2, false, inertial);
        for changed in [
            run_case(0.8, 1.0, 0.2, false, inertial),
            run_case(0.2, 2.0, 0.2, false, inertial),
            run_case(0.2, 1.0, 2.0, false, inertial),
            run_case(0.2, 1.0, 0.2, true, inertial),
        ] {
            assert!((changed.0[0] - baseline.0[0]).abs() > 1e-4);
            assert!((changed.0[1] - baseline.0[1]).abs() > 1e-4);
            assert!((changed.1 - baseline.1).abs() > 1e-8);
        }
    }
}

#[test]
fn complete_assembly_preserves_equilibrium_and_actual_newmark_motion() {
    with_network(0.2, 1.0, 0.2, false, 100000.0, |cx, thermal, source, prepared| {
        let fixed: Vec<_> = (0..15).filter(|d| *d != 3).collect();
        let initial = source.geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx).unwrap();
        let loads = source.geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 1.0)], cx).unwrap();
        let body = ThermoelasticSolid::try_new(
            source, thermal, &initial, &fixed, &[[280.0, 360.0]; 2],
            mechanical(), solve(), cx,
        ).unwrap();
        let record = Steps { count: 0, cancel_at: usize::MAX };
        let saved_cells = thermal.cells().clone();
        let original_displacement = body.equilibrium().solution.displacement.clone();
        // A borrowed trait object also exercises the public ?Sized observer alias.
        let assembly: &dyn HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError> = prepared;
        let (native, native_audit) = body.advance_heat_network_candidate(
            thermal, prepared, &record, 0.1, &loads, thermal_controls(), solve(), cx,
        ).unwrap();
        let (complete, audit) = body.advance_heat_assembly_candidate(
            thermal, assembly, &record, 0.1, &loads, thermal_controls(), solve(), cx,
        ).unwrap();
        assert_eq!(complete.thermal.cells(), native.thermal.cells());
        assert_eq!(complete.solid.equilibrium().solution.displacement,
            native.solid.equilibrium().solution.displacement);
        assert_eq!(complete.solid.time_s(), native.solid.time_s());
        assert_eq!(complete.report.mechanical_work_j, native.report.mechanical_work_j);
        assert_eq!(complete.observer.second, native.observer.second);
        assert_eq!(complete.observer.first.sample().unwrap().radiation,
            native.observer.first.sample().unwrap().radiation);
        assert_eq!(audit, native_audit);
        assert!(complete.report.cells.iter().any(|cell| cell.reversible_heat_j != 0.0));
        assert_ne!(complete.solid.equilibrium().solution.displacement, original_displacement);
        assert_eq!(body.equilibrium().solution.displacement, original_displacement);
        assert_eq!(body.time_s(), 0.0);
        let mut velocity = [[0.0; 3]; 5];
        velocity[1][0] = 0.02;
        let moving = body.begin_transient(
            &initial,
            &velocity,
            TransientConfig3 {
                maximum_steps: 8,
                maximum_displacement_m: 0.05,
                absolute_force_tolerance_n: 1e-11,
                absolute_energy_tolerance_j: 1e-10,
                relative_energy_tolerance: 1e-8,
                ..TransientConfig3::default()
            },
            cx,
        ).unwrap();
        let original_motion = moving.nodal_motion(cx).unwrap();
        let (native, native_audit) = moving.advance_heat_network_candidate(
            thermal, prepared, &record, 0.1, &loads, thermal_controls(), solve(), cx,
        ).unwrap();
        let (complete, audit) = moving.advance_heat_assembly_candidate(
            thermal, assembly, &record, 0.1, &loads, thermal_controls(), solve(), cx,
        ).unwrap();
        assert_eq!(complete.thermal.cells(), native.thermal.cells());
        assert_eq!(complete.solid.nodal_motion(cx).unwrap(), native.solid.nodal_motion(cx).unwrap());
        assert_eq!(complete.solid.step(), 1);
        assert_eq!(complete.solid.time_s(), native.solid.time_s());
        assert_eq!(complete.report.kinetic_energy_change_j, native.report.kinetic_energy_change_j);
        assert_eq!(complete.report.mechanical_work_j, native.report.mechanical_work_j);
        assert_eq!(complete.observer.second, native.observer.second);
        assert_eq!(complete.observer.first.sample().unwrap().radiation,
            native.observer.first.sample().unwrap().radiation);
        assert_eq!(audit, native_audit);
        assert_ne!(complete.solid.nodal_motion(cx).unwrap(), original_motion);
        assert_eq!(moving.nodal_motion(cx).unwrap(), original_motion);
        assert_eq!(moving.step(), 0);
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(record.count, 0);
    });
}

#[test]
fn g4_shared_source_admission_and_late_failures_preserve_every_original_history() {
    for pressure in [100000.0, 100001.0] {
        with_network(
            0.2,
            1.0,
            0.2,
            false,
            pressure,
            |cx, thermal, source, prepared| {
                let fixed: Vec<_> = (0..15).filter(|d| *d != 3).collect();
                let loads = source
                    .geometry()
                    .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
                    .unwrap();
                let body = ThermoelasticSolid::try_new(
                    source,
                    thermal,
                    &loads,
                    &fixed,
                    &[[280.0, 360.0]; 2],
                    mechanical(),
                    solve(),
                    cx,
                )
                .unwrap();
                let cells = thermal.cells().clone();
                let initial = body.equilibrium().solution.displacement.clone();
                assert!(prepared.resolve(cx, thermal).is_ok());
                let record = Steps {
                    count: 0,
                    cancel_at: usize::MAX,
                };
                if pressure != 100000.0 {
                    assert!(matches!(
                        body.advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx
                        ),
                        Err(ThermoelasticHeatError::Invalid(
                            "heat and reciprocal solid material conditions disagree"
                        ))
                    ));
                    assert!(matches!(
                        body.advance_heat_assembly_candidate(
                            thermal, prepared, &record, 0.1, &loads,
                            thermal_controls(), solve(), cx,
                        ),
                        Err(HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Invalid(
                            "heat and reciprocal solid material conditions disagree"
                        )))
                    ));
                } else {
                    let first = body
                        .advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    let cancel = Steps {
                        count: 0,
                        cancel_at: 2,
                    };
                    assert!(
                        body.advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &cancel,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx
                        )
                        .is_err()
                    );
                    assert!(matches!(
                        body.advance_heat_assembly_candidate(
                            thermal, prepared, &cancel, 0.1, &loads,
                            thermal_controls(), solve(), cx,
                        ),
                        Err(HeatAssemblyAdvanceError::Coupled(_))
                    ));
                    assert_eq!(cancel.count, 0);
                    let mut short = solve();
                    short.maximum_iterations = 1;
                    assert!(
                        body.advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            short,
                            cx
                        )
                        .is_err()
                    );
                    let narrow = ThermoelasticSolid::try_new(
                        source,
                        thermal,
                        &loads,
                        &fixed,
                        &[[280.0, 300.01]; 2],
                        mechanical(),
                        solve(),
                        cx,
                    )
                    .unwrap();
                    assert!(
                        narrow
                            .advance_heat_network_candidate(
                                thermal,
                                prepared,
                                &record,
                                0.1,
                                &loads,
                                thermal_controls(),
                                solve(),
                                cx
                            )
                            .is_err()
                    );
                    let retry = body
                        .advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(retry.0.thermal.cells(), first.0.thermal.cells());
                    assert_eq!(
                        retry.0.solid.equilibrium().solution.displacement,
                        first.0.solid.equilibrium().solution.displacement
                    );
                    assert_eq!(retry.0.observer.second, first.0.observer.second);
                }
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(body.equilibrium().solution.displacement, initial);
                assert_eq!(body.time_s(), 0.0);
                assert_eq!(record.count, 0);
                let mut velocity = [[0.0; 3]; 5];
                velocity[1][0] = 0.02;
                let moving = body
                    .begin_transient(
                        &loads,
                        &velocity,
                        TransientConfig3 {
                            maximum_steps: 8,
                            maximum_displacement_m: 0.05,
                            absolute_force_tolerance_n: 1e-11,
                            absolute_energy_tolerance_j: 1e-10,
                            relative_energy_tolerance: 1e-8,
                            ..TransientConfig3::default()
                        },
                        cx,
                    )
                    .unwrap();
                let motion = moving.nodal_motion(cx).unwrap();
                let cancel = Steps {
                    count: 0,
                    cancel_at: 2,
                };
                if pressure != 100000.0 {
                    assert!(matches!(
                        moving.advance_heat_assembly_candidate(
                            thermal, prepared, &record, 0.1, &loads,
                            thermal_controls(), solve(), cx,
                        ),
                        Err(HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Invalid(
                            "heat and reciprocal solid material conditions disagree"
                        )))
                    ));
                } else {
                    assert!(matches!(
                        moving.advance_heat_assembly_candidate(
                            thermal, prepared, &cancel, 0.1, &loads,
                            thermal_controls(), solve(), cx,
                        ),
                        Err(HeatAssemblyAdvanceError::Coupled(_))
                    ));
                }
                assert!(
                    moving
                        .advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &cancel,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx
                        )
                        .is_err()
                );
                assert_eq!(moving.nodal_motion(cx).unwrap(), motion);
                assert_eq!(moving.step(), 0);
                assert_eq!(moving.time_s(), 0.0);
                assert_eq!(cancel.count, 0);
                assert_eq!(thermal.cells(), &cells);
                if pressure == 100000.0 {
                    let first = moving
                        .advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    let retry = moving
                        .advance_heat_network_candidate(
                            thermal,
                            prepared,
                            &record,
                            0.1,
                            &loads,
                            thermal_controls(),
                            solve(),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(retry.0.thermal.cells(), first.0.thermal.cells());
                    assert_eq!(
                        retry.0.solid.nodal_motion(cx).unwrap(),
                        first.0.solid.nodal_motion(cx).unwrap()
                    );
                    assert_eq!(retry.0.observer.second, first.0.observer.second);
                }
            },
        );
    }
}
