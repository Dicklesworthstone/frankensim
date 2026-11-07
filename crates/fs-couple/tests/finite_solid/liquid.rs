//! Original accepted thermal/reference mass enters the existing moving liquid.
//! Synthetic source cards and analytic cubes; these are conservation fixtures.
use super::*;
use fs_couple::finite_solid::liquid::*;
use fs_couple::finite_solid::thermal::MeanThermalReference;
use fs_ga::Vec3;
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyControls;
use fs_lbm::material::free_surface::transfer::TetLiquidTransferLimits;
use fs_lbm::material::free_surface::transfer::TetReferenceDensityTransferControls;
use fs_lbm::material::free_surface::transfer::enthalpy::TetEnthalpyTransferLimits;
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain3, plan_free_surface_from_material};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, resolve_liquid_surface_state,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
use fs_mbd_material::density_field::{DensityVertexBinding, TetrahedralDensityBody};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::sync::Arc;

const DX: f64 = 0.001;
const DT: f64 = 0.0002;
const SIDE: f64 = 0.004;
const RHO: f64 = 1000.0;
const TIME: f64 = 0.125;
const H: f64 = 200000.0;
const CUBE_TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 3, 2, 7],
    [0, 2, 6, 7],
    [0, 6, 4, 7],
    [0, 4, 5, 7],
    [0, 5, 1, 7],
];
const SPIN: [f64; 3] = [1.7, -2.4, 3.1];

fn card(reference_density_scale: f64) -> MaterialCard {
    let base = source_card(RHO, 1.0, false, false);
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        (
            "density",
            Density::DIMS,
            vec![
                (300.0, reference_density_scale * RHO * 1.25_f64.powi(3)),
                (350.0, reference_density_scale * RHO * 1.125_f64.powi(3)),
                (390.0, RHO),
                (500.0, RHO),
            ],
        ),
        (
            "dynamic_viscosity",
            Dims([-1, 1, -1, 0, 0, 0]),
            vec![(300.0, 0.1), (500.0, 0.02)],
        ),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            vec![(300.0, 2000.0), (500.0, 2000.0)],
        ),
        (
            "thermal_conductivity",
            Dims([1, 1, -3, -1, 0, 0]),
            vec![(300.0, 0.5), (500.0, 0.5)],
        ),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    300.0,
                    500.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: base.models()[0].provenance.clone(),
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                "mean_engineering_thermal_expansion_coefficient",
                INVERSE_TEMPERATURE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: 0.0025,
                dims: INVERSE_TEMPERATURE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T_reference", temperature, 300.0, 350.0)
                .with_quantity("T_endpoint", temperature, 390.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: base.models()[0].provenance.clone(),
        })
        .unwrap();
    let mut models = base.models().to_vec();
    models[0].validity =
        ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 500.0);
    MaterialCard::assemble(base.id().clone(), claims, models).unwrap()
}
fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 180000.0,
                temperature_k: 390.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: RHO,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 400000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: RHO,
            },
        ],
    )
    .unwrap()
}
fn with_sources(
    reference_density_scale: f64,
    f: impl FnOnce(
        &MaterialCard,
        &LiquidTemperatureLaw<'_>,
        &ResolvedLiquidSurfaceState,
        &EquilibriumEnthalpyPhaseCurve,
    ),
) {
    let card = card(reference_density_scale);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card);
    f(&card, &law, &surface, &caloric);
}
fn material(card: &MaterialCard, reference_k: f64) -> FiniteSolidMaterial {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let query = QueryPoint::new()
        .with_quantity("T_reference", temperature, reference_k)
        .unwrap()
        .with_quantity("T_endpoint", temperature, 400.0)
        .unwrap();
    let reference = MeanThermalReference::resolve(
        card,
        &PropertyKey::new(
            "mean_engineering_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
        ),
        &query,
        "T_reference",
        "T_endpoint",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    FiniteSolidMaterial::resolve_mean_thermal(
        card,
        0,
        reference,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn release_controls() -> FiniteSolidLiquidReleaseControls {
    FiniteSolidLiquidReleaseControls {
        transfer: TetEnthalpyTransferLimits {
            mechanical: TetLiquidTransferLimits {
                max_tetrahedra: 16,
                max_overlap_tests: 30000,
                max_wet_cells: 4096,
                relative_volume_tolerance: 1e-10,
                mass_tolerance_kg: 1e-13,
                momentum_tolerance_kg_m_s: 1e-15,
                angular_momentum_tolerance_kg_m2_s: 1e-17,
                maximum_centroid_shift_m: 1e-10,
                maximum_projection_energy_change_j: 1e-6,
            },
            max_source_queries: 2048,
            enthalpy_tolerance_j: 1e-9,
        },
        thermal: FreeSurfaceEnthalpyControls::default(),
        density: TetReferenceDensityTransferControls {
            maximum_density_defect_kg_m3: 1e-8,
            maximum_volume_defect_m3: 1e-16,
        },
        maximum_cauchy_stress_pa: 1e-6,
        maximum_elastic_storage_j: 1e-14,
        maximum_inertial_force_n: 1e-10,
        source_kinetic_tolerance_j: 1e-14,
        source_centroid_tolerance_m: 1e-11,
        source_temperature_tolerance_k: 1e-10,
        source_reference_volume_tolerance_m3: 1e-17,
    }
}
fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [16; 3],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 4096,
    }
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn translation(piece: usize) -> [f64; 3] {
    [0.01 * (piece + 1) as f64, -0.006, 0.003]
}
fn center(piece: usize) -> [f64; 3] {
    [0.004 + 0.008 * piece as f64, 0.008, 0.008]
}
struct Geometry {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    displacement: Vec<[f64; 3]>,
    velocity: Vec<[f64; 3]>,
    reference_volumes: Vec<f64>,
}
fn geometry(stretches: &[f64]) -> Geometry {
    let mut out = Geometry {
        nodes: Vec::new(),
        tets: Vec::new(),
        displacement: Vec::new(),
        velocity: Vec::new(),
        reference_volumes: Vec::new(),
    };
    for (piece, &stretch) in stretches.iter().enumerate() {
        for corner in 0..8 {
            let local: [f64; 3] = core::array::from_fn(|axis| {
                if corner & (1 << axis) == 0 {
                    -0.5 * SIDE
                } else {
                    0.5 * SIDE
                }
            });
            out.nodes.push(core::array::from_fn(|axis| {
                center(piece)[axis] + local[axis] / stretch
            }));
            out.displacement.push(local.map(|x| x - x / stretch));
            let spin = cross(SPIN, local);
            out.velocity.push(core::array::from_fn(|axis| {
                translation(piece)[axis] + spin[axis]
            }));
        }
        out.tets
            .extend(CUBE_TETS.map(|tet| tet.map(|i| i + 8 * piece)));
        out.reference_volumes
            .extend([((SIDE / stretch).powi(3)) / 6.0; 6]);
    }
    out
}
fn source_cells(geometry: &Geometry, reference_density: &[f64]) -> Vec<DepositedCell> {
    geometry
        .tets
        .iter()
        .enumerate()
        .map(|(i, tet)| DepositedCell {
            id: CellId(10 + i as u64),
            segment_id: i as u64,
            body_id: 77,
            material_id: MaterialId(0),
            mass_kg: reference_density[i] * geometry.reference_volumes[i],
            specific_h_j_kg: H,
            geometry: CellGeometry {
                start_m: geometry.nodes[tet[0]],
                end_m: geometry.nodes[tet[1]],
                volume_m3: geometry.reference_volumes[i],
                surface_area_m2: SIDE * SIDE,
            },
            boundary: Boundary {
                ambient_k: 400.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        })
        .collect()
}
fn timeline<'a>(
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    cells: Vec<DepositedCell>,
    time: f64,
    future: bool,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells
        .iter()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if future {
        let mut cell = cells[0].clone();
        cell.id = CellId(99);
        events.push(TimedEvent {
            time_s: 1.0,
            event: Event::Birth(cell),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), caloric as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 64,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        time,
        StepConfig {
            max_dt_s: 0.125,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn body_problem<'a>(
    geometry: &'a Geometry,
    materials: &'a [FiniteSolidMaterial],
) -> HyperTetProblem<'a, FiniteSolidMaterial> {
    HyperTetProblem {
        nodes_m: &geometry.nodes,
        tetrahedra: &geometry.tets,
        material: HyperelasticLawField::PerElement(materials),
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

#[test]
fn g1_unequal_original_densities_and_expansions_retain_mass_spin_heat_and_accepted_clock() {
    with_sources(1.0, |card, law, surface, caloric| {
        with_cx(|gate, cx| {
            let geometry = geometry(&[1.25, 1.125]);
            let materials: Vec<_> = [300.0, 350.0]
                .into_iter()
                .flat_map(|t| vec![material(card, t); 6])
                .collect();
            let queries: Vec<_> = [300.0, 350.0]
                .into_iter()
                .flat_map(|t| vec![point(t); 8])
                .collect();
            let vertices: Vec<_> = geometry
                .nodes
                .iter()
                .zip(&queries)
                .map(|(p, query)| DensityVertexBinding {
                    position_m: Vec3::new(p[0], p[1], p[2]),
                    material: card,
                    point: query,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                })
                .collect();
            let density =
                TetrahedralDensityBody::resolve(&vertices, &geometry.tets, vertices.len(), || {
                    false
                })
                .unwrap();
            assert!(density.vertices()[0].density_kg_m3() > density.vertices()[8].density_kg_m3());
            let body = FiniteSolidBody::new_with_reference_density_field(
                body_problem(&geometry, &materials),
                Arc::new(density),
                &geometry.displacement,
                &geometry.velocity,
                TIME,
                Default::default(),
                cx,
            )
            .unwrap();
            let density: Vec<_> = materials
                .iter()
                .map(FiniteSolidMaterial::reference_density_kg_m3)
                .collect();
            let cells = source_cells(&geometry, &density);
            let ids: Vec<_> = cells.iter().map(|c| c.id).collect();
            let thermal = timeline(caloric, cells, TIME, false);
            let before = (
                body.native_state().clone(),
                thermal.cells().clone(),
                thermal.time_s(),
            );
            let scaling = plan_free_surface_from_material(surface, DX, DT, 0.1).unwrap();
            let (mut liquid, report) = body
                .release_relaxed_liquid_candidate(
                    &thermal,
                    &ids,
                    scaling,
                    domain(),
                    [0.0; 3],
                    law,
                    caloric,
                    release_controls(),
                    cx,
                    gate,
                )
                .unwrap();
            // Independent current cuboid integrals; rho_ref*V_ref = rho_current*V_current.
            let piece_mass = RHO * SIDE.powi(3);
            let inertia = piece_mass * SIDE * SIDE / 6.0;
            let mut expected_p = [0.0; 3];
            let mut expected_l = [0.0; 3];
            let mut expected_k = 0.0;
            for piece in 0..2 {
                let velocity = translation(piece);
                let momentum = velocity.map(|v| piece_mass * v);
                let orbital = cross(center(piece), momentum);
                for axis in 0..3 {
                    expected_p[axis] += momentum[axis];
                    expected_l[axis] += orbital[axis] + inertia * SPIN[axis];
                    expected_k += 0.5 * piece_mass * velocity[axis].powi(2)
                        + 0.5 * inertia * SPIN[axis].powi(2);
                }
            }
            for inventory in [report.solid, report.transfer.mechanical.source] {
                near(inventory.mass_kg, 2.0 * piece_mass, 2e-18);
                near(inventory.kinetic_energy_j, expected_k, 2e-20);
                for axis in 0..3 {
                    near(inventory.momentum_kg_m_s[axis], expected_p[axis], 2e-19);
                    near(
                        inventory.angular_momentum_kg_m2_s[axis],
                        expected_l[axis],
                        2e-21,
                    );
                }
            }
            let density = report
                .transfer
                .mechanical
                .reference_density
                .as_ref()
                .unwrap();
            let original_volume: f64 = geometry.reference_volumes.iter().sum();
            near(density.reference.volume, original_volume, 1e-20);
            near(density.deformed.volume, 2.0 * SIDE.powi(3), 1e-20);
            near(density.native_occupied_volume_m3, 2.0 * SIDE.powi(3), 1e-20);
            assert!(density.reference.volume < 0.8 * density.deformed.volume);
            assert!(density.maximum_current_density_defect_kg_m3 < 1e-8);
            near(report.thermal_enthalpy_j, 2.0 * piece_mass * H, 1e-11);
            near(
                liquid.total_enthalpy_j().unwrap(),
                2.0 * piece_mass * H,
                1e-9,
            );
            assert_eq!(report.transfer.mechanical.components.len(), 2);
            assert!(report.retained_kinetic_change_j().to_f64().abs() > 1e-13);
            near(
                report.retained_mechanical_change_j().to_f64(),
                report.retained_kinetic_change_j().to_f64() - report.elastic_storage_j,
                1e-24,
            );
            // Observe actual native cells, independent of either transfer receipt.
            let (mut native_mass, mut native_h, mut native_volume) = (0.0, 0.0, 0.0);
            let (mut native_p, mut native_l) = ([0.0; 3], [0.0; 3]);
            for piece in &report.transfer.mechanical.components {
                for &xyz in &piece.target_cells {
                    let sample = liquid
                        .flow()
                        .sample_at_si(xyz[0], xyz[1], xyz[2])
                        .unwrap()
                        .unwrap();
                    let volume = sample.fill_fraction * DX.powi(3);
                    let mass = sample.density_kg_m3 * volume;
                    let momentum = sample.velocity_m_s.map(|v| mass * v);
                    let angular = cross(xyz.map(|i| (i as f64 + 0.5) * DX), momentum);
                    native_volume += volume;
                    native_mass += mass;
                    native_h += mass * liquid.specific_enthalpy_at(xyz).unwrap();
                    for axis in 0..3 {
                        native_p[axis] += momentum[axis];
                        native_l[axis] += angular[axis];
                    }
                }
            }
            near(native_mass, 2.0 * piece_mass, 2e-18);
            near(native_volume, 2.0 * SIDE.powi(3), 1e-20);
            near(native_h, 2.0 * piece_mass * H, 1e-9);
            for axis in 0..3 {
                near(native_p[axis], expected_p[axis], 2e-19);
                near(native_l[axis], expected_l[axis], 2e-21);
            }
            near(liquid.physical_time_s(), TIME, 0.0);
            let step = liquid.step(gate).unwrap();
            near(step.enthalpy_before_j, 2.0 * piece_mass * H, 1e-9);
            near(step.external_heat_j, 0.0, 0.0);
            near(step.viscous_heat_j, 0.0, 0.0);
            near(liquid.physical_time_s(), TIME + DT, 0.0);
            assert_eq!(
                (
                    body.native_state().clone(),
                    thermal.cells().clone(),
                    thermal.time_s()
                ),
                before
            );
        });
    });
}

#[test]
fn g4_original_chart_clock_cell_inventory_ports_and_shared_enthalpy_cannot_be_replaced() {
    with_sources(1.0, |card, law, surface, caloric| {
        with_cx(|gate, cx| {
            let geometry = geometry(&[1.25]);
            let materials = vec![material(card, 300.0); 6];
            let body = FiniteSolidBody::new(
                body_problem(&geometry, &materials),
                &geometry.displacement,
                &geometry.velocity,
                TIME,
                Default::default(),
                cx,
            )
            .unwrap();
            let cells = source_cells(&geometry, &[materials[0].reference_density_kg_m3(); 6]);
            let ids: Vec<_> = cells.iter().map(|c| c.id).collect();
            let thermal = timeline(caloric, cells.clone(), TIME, false);
            let before = (body.native_state().clone(), thermal.cells().clone());
            let attempt = |thermal: &DepositionSimulation<'_>,
                           chart: &EquilibriumEnthalpyPhaseCurve,
                           controls: FiniteSolidLiquidReleaseControls| {
                body.release_relaxed_liquid_candidate(
                    thermal,
                    &ids,
                    plan_free_surface_from_material(surface, DX, DT, 0.1).unwrap(),
                    domain(),
                    [0.0; 3],
                    law,
                    chart,
                    controls,
                    cx,
                    gate,
                )
                .map(|_| ())
                .map_err(|error| error.to_string())
            };
            let equal_but_foreign_chart = chart(card);
            assert!(attempt(&thermal, &equal_but_foreign_chart, release_controls()).is_err());
            let late = timeline(caloric, cells.clone(), TIME + 0.125, false);
            assert!(attempt(&late, caloric, release_controls()).is_err());
            let future = timeline(caloric, cells.clone(), TIME, true);
            assert!(attempt(&future, caloric, release_controls()).is_err());
            let mut wrong_mass = cells.clone();
            wrong_mass[0].mass_kg *= 1.01;
            assert!(
                attempt(
                    &timeline(caloric, wrong_mass, TIME, false),
                    caloric,
                    release_controls()
                )
                .is_err()
            );
            let mut shared_h = cells.clone();
            shared_h[1].specific_h_j_kg += 2000.0;
            let mut controls = release_controls();
            controls.source_temperature_tolerance_k = 2.0;
            let error =
                attempt(&timeline(caloric, shared_h, TIME, false), caloric, controls).unwrap_err();
            assert!(error.contains("incident thermal cells"), "{error}");
            let mut heated = cells.clone();
            heated[0].boundary.source_w = 1e-6;
            assert!(
                attempt(
                    &timeline(caloric, heated, TIME, false),
                    caloric,
                    release_controls()
                )
                .is_err()
            );
            let mut projection = release_controls();
            projection
                .transfer
                .mechanical
                .maximum_projection_energy_change_j = 0.0;
            assert!(attempt(&thermal, caloric, projection).is_err());
            assert_eq!(
                (body.native_state().clone(), thermal.cells().clone()),
                before
            );
            // The original scalar-density body is still available for a successful retry.
            attempt(&thermal, caloric, release_controls()).unwrap();
            gate.request();
            assert!(matches!(
                body.release_relaxed_liquid_candidate(
                    &thermal,
                    &ids,
                    plan_free_surface_from_material(surface, DX, DT, 0.1).unwrap(),
                    domain(),
                    [0.0; 3],
                    law,
                    caloric,
                    release_controls(),
                    cx,
                    gate,
                ),
                Err(FiniteSolidLiquidReleaseError::Cancelled)
            ));
            assert_eq!(
                (body.native_state().clone(), thermal.cells().clone()),
                before
            );
        });
    });
}

#[test]
fn g4_retained_stress_supports_and_incompatible_physical_density_refuse_without_state_loss() {
    for reference_density_scale in [1.0, 1.1] {
        with_sources(reference_density_scale, |card, law, surface, caloric| {
            with_cx(|gate, cx| {
                let geometry = geometry(&[1.25]);
                let materials = vec![material(card, 300.0); 6];
                let cells = source_cells(&geometry, &[materials[0].reference_density_kg_m3(); 6]);
                let ids: Vec<_> = cells.iter().map(|c| c.id).collect();
                let thermal = timeline(caloric, cells, TIME, false);
                let mut displacement = geometry.displacement.clone();
                if reference_density_scale == 1.0 {
                    // A further elastic stretch has real native stress and stored energy.
                    for (u, x) in displacement.iter_mut().zip(&geometry.nodes) {
                        u[0] += 0.01 * 1.25 * (x[0] - center(0)[0]);
                    }
                }
                let body = FiniteSolidBody::new(
                    body_problem(&geometry, &materials),
                    &displacement,
                    &geometry.velocity,
                    TIME,
                    Default::default(),
                    cx,
                )
                .unwrap();
                let before = (body.native_state().clone(), thermal.cells().clone());
                let error = body
                    .release_relaxed_liquid_candidate(
                        &thermal,
                        &ids,
                        plan_free_surface_from_material(surface, DX, DT, 0.1).unwrap(),
                        domain(),
                        [0.0; 3],
                        law,
                        caloric,
                        release_controls(),
                        cx,
                        gate,
                    )
                    .err()
                    .expect("physical incompatibility must refuse");
                let reason = error.to_string();
                assert!(
                    reason.contains(if reference_density_scale == 1.0 {
                        "elastic"
                    } else {
                        "density"
                    }),
                    "{reason}"
                );
                assert_eq!(
                    (body.native_state().clone(), thermal.cells().clone()),
                    before
                );
                if reference_density_scale == 1.0 {
                    let fixed = [(0, geometry.displacement[0][0])];
                    let mut velocity = geometry.velocity.clone();
                    velocity[0][0] = 0.0;
                    let mut problem = body_problem(&geometry, &materials);
                    problem.prescribed_m = &fixed;
                    let supported = FiniteSolidBody::new(
                        problem,
                        &geometry.displacement,
                        &velocity,
                        TIME,
                        Default::default(),
                        cx,
                    )
                    .unwrap();
                    let error = supported
                        .release_relaxed_liquid_candidate(
                            &thermal,
                            &ids,
                            plan_free_surface_from_material(surface, DX, DT, 0.1).unwrap(),
                            domain(),
                            [0.0; 3],
                            law,
                            caloric,
                            release_controls(),
                            cx,
                            gate,
                        )
                        .err()
                        .expect("support owner cannot be silently removed");
                    assert!(error.to_string().contains("supports"));
                }
            });
        });
    }
}
