//! Two sourced adhesive disks share finite inertia and finite thermal storage.
//! The independent reference eliminates the three caloric balances and solves
//! only the two physical face temperatures using analytic uniform-disk torque.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::spinning::SourceAdhesiveSlidingPatch;
use fs_contact::interface_binding::adhesive_shear::*;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, JkrResponse, SourceJkrSphereContact};
use fs_couple::friction_heat::adhesive::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::multiple::RigidAdhesiveIslandControls;
use fs_mbd::tangential_friction::adhesive::{RigidAdhesiveContact, RigidAdhesiveControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::{Dims, Pressure, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::InputAuthority;
use fs_tribo::adhesive::spinning::contact::{AdhesiveWrenchConfig, AdhesiveWrenchMode};
use fs_tribo::partial_slip::TangentFrame;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 0.05;
const HALF_LENGTH: f64 = 0.0625;
const MASS: [f64; 3] = [2.0, 3.0, 4.0];
const INERTIA: [f64; 3] = [0.08, 0.12, 0.15];
const SPIN: [f64; 3] = [4.0, 0.0, 1.0];
const CAPACITY: [f64; 3] = [0.03, 0.02, 0.04];
const K: [f64; 3] = [2.0, 1.0, 3.0];
const STRENGTH: [f64; 2] = [2000.0, 1500.0];
const SLOPE: [f64; 2] = [20.0, 15.0];
const LINK: [u64; 2] = [20, 10];
const PATCH: [&str; 2] = ["lower-junction", "upper-junction"];
const INITIAL: [f64; 3] = [302.0, 292.0, 310.0];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn claim(name: &str, dims: Dims, values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (320.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured adhesive heat fixture: {name}; not measured data"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

struct Sources {
    normal: [BoundNormalInterfaceState; 2],
    shear: [SourceAdhesiveShear; 2],
    response: [JkrResponse; 2],
}
impl Sources {
    fn new() -> Self {
        let build = |i: usize| {
            let a = common::solid("adhesive-heat-middle", 1000.0, 1e6);
            let b = common::solid(&format!("adhesive-heat-outer-{i}"), 1000.0, 2e6);
            let mut claims = ClaimSet::new();
            claims
                .insert_claim(claim(
                    ADHESION_ENERGY_PROPERTY,
                    ADHESION_ENERGY_DIMS,
                    [1.0; 2],
                ))
                .unwrap();
            let values = [STRENGTH[i] + 20.0 * SLOPE[i], STRENGTH[i] - 20.0 * SLOPE[i]];
            claims
                .insert_claim(claim(
                    STATIC_JUNCTION_SHEAR_PROPERTY,
                    Pressure::DIMS,
                    values.map(|v| 1.5 * v),
                ))
                .unwrap();
            claims
                .insert_claim(claim(
                    KINETIC_JUNCTION_SHEAR_PROPERTY,
                    Pressure::DIMS,
                    values,
                ))
                .unwrap();
            let card = InterfaceSystemCard::assemble(
                SurfaceSpec {
                    material: a.id().clone(),
                    texture_frame: "middle-face".into(),
                },
                SurfaceSpec {
                    material: b.id().clone(),
                    texture_frame: "outer-face".into(),
                },
                SystemContext {
                    medium: "dry".into(),
                    third_body: None,
                    environment: "manufactured".into(),
                    history: "initial".into(),
                },
                claims,
                vec![],
            )
            .unwrap();
            let query = QueryPoint::new().with("T", 300.0).unwrap();
            let selection = MaterialPropertySelection::SingleClaimOnly;
            let sa = resolve_isotropic_solid_state_point(&a, &query, selection.clone()).unwrap();
            let sb = resolve_isotropic_solid_state_point(&b, &query, selection.clone()).unwrap();
            let dry =
                bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
            let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
            let adhesion = resolve_interface_state_point(
                &card,
                &query,
                &[ScalarPropertyRequirement::try_new(
                    ADHESION_ENERGY_PROPERTY,
                    ADHESION_ENERGY_DIMS,
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap()],
                selection.clone(),
            )
            .unwrap();
            let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
            let shear = SourceAdhesiveShear::resolve(&normal, selection).unwrap();
            let owner = SourceJkrSphereContact::new(
                normal.clone(),
                1.0,
                1e-7,
                JkrEnvelope {
                    maximum_indentation_m: 0.01,
                    half_space_depth_m: 1.0,
                    layer_thickness_m: 1.0,
                    maximum_patch_to_radius: 0.1,
                    maximum_patch_to_depth: 0.1,
                    maximum_patch_to_layer: 0.1,
                    maximum_nominal_strain: 0.1,
                    maximum_repulsive_pressure_to_yield: 0.1,
                    maximum_process_zone_to_patch: 0.01,
                    minimum_tabor_parameter: 5.0,
                },
            )
            .unwrap();
            let a0 = (4.5 * PI / normal.elastic().reduced_modulus_pa()).cbrt();
            let indentation = if i == 0 {
                a0 * a0 - (2.0 * PI * a0 / normal.elastic().reduced_modulus_pa()).sqrt()
            } else {
                0.0
            };
            let response = owner.evaluate_attached(indentation).unwrap();
            (normal, shear, response)
        };
        let [(n0, s0, r0), (n1, s1, r1)] = [build(0), build(1)];
        Self {
            normal: [n0, n1],
            shear: [s0, s1],
            response: [r0, r1],
        }
    }
    fn patches(&self) -> [SourceAdhesiveSlidingPatch<'_>; 2] {
        std::array::from_fn(|i| {
            self.shear[i]
                .sliding_patch_jkr(&self.normal[i], &self.response[i])
                .unwrap()
        })
    }
    fn areas(&self) -> [f64; 2] {
        self.response.map(|r| r.contact_area_m2)
    }
    fn radii(&self) -> [f64; 2] {
        self.response.map(|r| r.contact_radius_m)
    }
}

struct Island {
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
}
impl Island {
    fn new(moving: bool) -> Self {
        Self {
            properties: std::array::from_fn(|i| {
                MassProperties::new(
                    MASS[i],
                    Vec3::ZERO,
                    Vec3::new(INERTIA[i], INERTIA[i], INERTIA[i]),
                )
                .unwrap()
            }),
            bodies: std::array::from_fn(|i| {
                RigidBodyState::new(
                    Pose::new(
                        Vec3::new(0.0, 0.0, [0.0, -2.0 * HALF_LENGTH, 2.0 * HALF_LENGTH][i]),
                        UnitQuaternion::IDENTITY,
                    )
                    .unwrap(),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, if moving { INERTIA[i] * SPIN[i] } else { 0.0 }),
                )
                .unwrap()
            }),
        }
    }
    fn contacts<'a>(
        &self,
        sources: &'a Sources,
        routes: &'a [FrictionHeatRoute; 2],
    ) -> [RigidAdhesiveIslandHeatContact<'a>; 2] {
        let patches = sources.patches();
        std::array::from_fn(|i| {
            let arm = [-HALF_LENGTH, HALF_LENGTH][i];
            RigidAdhesiveIslandHeatContact {
                patch_id: PATCH[i],
                temperature_axis: "T",
                route: &routes[i],
                temperature_guess_k: 300.0,
                contact: RigidAdhesiveContact {
                    binding: RigidContactBinding {
                        body_a: 0,
                        body_b: i + 1,
                        arm_a_body: Vec3::new(0.0, 0.0, arm),
                        arm_b_body: Vec3::new(0.0, 0.0, -arm),
                    },
                    // A shared world-z spin convention keeps the independent
                    // angular response W_ij=1/I_middle+delta_ij/I_outer explicit.
                    frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                    source: patches[i],
                },
            }
        })
    }
}

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 {
            Ok(300.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        from: f64,
        to: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(from)?;
        self.temperature_k(to)?;
        Ok(Some(((to - from) / t).ln_1p()))
    }
}
static CALORIC: Caloric = Caloric;

fn thermal(initial: [f64; 3], areas: [f64; 2]) -> DepositionSimulation<'static> {
    thermal_with_geometry(initial, areas, HALF_LENGTH, areas[0].max(areas[1]))
}

fn thermal_with_geometry(
    initial: [f64; 3],
    areas: [f64; 2],
    half_length: f64,
    face: f64,
) -> DepositionSimulation<'static> {
    let endpoints = [
        [-half_length, half_length],
        [-3.0 * half_length, -half_length],
        [half_length, 3.0 * half_length],
    ];
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: CAPACITY[i],
                specific_h_j_kg: initial[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, endpoints[i][0]],
                    end_m: [0.0, 0.0, endpoints[i][1]],
                    volume_m3: face * 2.0 * half_length,
                    surface_area_m2: 2.0 * face + 8.0 * face.sqrt() * half_length,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINK[i],
                a: CellId(0),
                b: CellId(i as u64 + 1),
                area_m2: areas[i],
                conductance_w_k: series_conductance_w_k(
                    areas[i],
                    [half_length; 2],
                    [K[0], K[i + 1]],
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &CALORIC as &dyn ThermalLaw),
            (MaterialId(1), &CALORIC as &dyn ThermalLaw),
            (MaterialId(2), &CALORIC as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn routes(thermal: &DepositionSimulation<'_>) -> [FrictionHeatRoute; 2] {
    std::array::from_fn(|i| {
        FrictionHeatRoute::from_link(thermal, LINK[i], [K[0], K[i + 1]]).unwrap()
    })
}

fn controls() -> RigidAdhesiveIslandHeatControls {
    RigidAdhesiveIslandHeatControls {
        mechanical: RigidAdhesiveIslandControls {
            mechanical: RigidAdhesiveControls {
                maximum_bodies: 3,
                geometry_tolerance_m: 1e-10,
                velocity_tolerance_m_s: 1e-10,
                angular_velocity_tolerance_rad_s: 1e-9,
                linear_momentum_tolerance_ns: 1e-10,
                angular_momentum_tolerance_nms: 1e-10,
                energy_tolerance_j: 1e-10,
                contact: AdhesiveWrenchConfig {
                    maximum_static_iterations: 96,
                    maximum_iterations: 256,
                    maximum_evaluations: 512,
                    force_tolerance_n: 2e-9,
                    moment_tolerance_nm: 2e-10,
                    displacement_tolerance_m: 2e-13,
                    energy_tolerance_j: 2e-11,
                    relative_matrix_tolerance: 1e-13,
                    maximum_force_n: 100.0,
                    maximum_spin_moment_nm: 10.0,
                },
            },
            maximum_contacts: 2,
            maximum_matrix_terms: 36,
            maximum_sweeps: 96,
            maximum_evaluations: 8192,
        },
        maximum_cells: 3,
        maximum_unknowns: 5,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.05,
        heat: FrictionHeatTolerances {
            cell_j: 1e-9,
            thermal_j: 1e-9,
            interface_j: 1e-9,
            total_j: 1e-9,
        },
        combined_energy_tolerance_j: 2e-9,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-11,
            linear_restart: 5,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn transport(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<RigidAdhesiveIslandHeatSource<Vec<f64>>, &'static str> {
    assert_eq!(ids, &[10, 20]);
    Ok(RigidAdhesiveIslandHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| {
                let i = LINK.iter().position(|id| *id == link_id).unwrap();
                RigidAdhesiveIslandHeatConductivity {
                    link_id,
                    conductivity_w_m_k: [K[0], K[i + 1]],
                }
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperatures.to_vec(),
        },
    })
}
fn run<'a, 'm>(
    island: &Island,
    contacts: &[RigidAdhesiveIslandHeatContact<'a>],
    thermal: &DepositionSimulation<'m>,
    controls: RigidAdhesiveIslandHeatControls,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidAdhesiveIslandHeatSource<Vec<f64>>, &'static str>,
    gate: &CancelGate,
) -> Result<
    HeatedRigidAdhesiveIslandCandidate<'a, 'm, Vec<f64>>,
    RigidAdhesiveIslandHeatError<&'static str>,
> {
    propose_rigid_adhesive_island_heat(
        &island.bodies,
        &island.properties,
        contacts,
        thermal,
        thermal.time_s(),
        DT,
        controls,
        &[0.0; 3],
        resolve,
        gate,
    )
}

struct Reference {
    faces: [f64; 2],
    cells: [f64; 3],
    heat: [f64; 2],
    moment: [f64; 2],
}
fn reference(initial: [f64; 3], areas: [f64; 2], radii: [f64; 2], moving: bool) -> Reference {
    let d: [[f64; 2]; 2] = std::array::from_fn(|i| {
        [
            DT * K[0] * areas[i] / HALF_LENGTH,
            DT * K[i + 1] * areas[i] / HALF_LENGTH,
        ]
    });
    let base: [f64; 3] = std::array::from_fn(|i| CAPACITY[i] * (initial[i] - 300.0));
    let denom = [
        CAPACITY[0] + d[0][0] + d[1][0],
        CAPACITY[1] + d[0][1],
        CAPACITY[2] + d[1][1],
    ];
    let matrix: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            -d[i][0] * d[j][0] / denom[0]
                + if i == j {
                    d[i][0] + d[i][1] - d[i][1] * d[i][1] / denom[i + 1]
                } else {
                    0.0
                }
        })
    });
    let rhs: [f64; 2] = std::array::from_fn(|i| {
        d[i][0] * base[0] / denom[0] + d[i][1] * base[i + 1] / denom[i + 1]
    });
    let w: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| 1.0 / INERTIA[0] + if i == j { 1.0 / INERTIA[i + 1] } else { 0.0 })
    });
    let factor: [f64; 2] = std::array::from_fn(|i| (2.0 / 3.0) * areas[i] * radii[i]);
    let slope: [f64; 2] = std::array::from_fn(|i| if moving { factor[i] * SLOPE[i] } else { 0.0 });
    let mut theta = [0.0; 2];
    for _ in 0..20 {
        let moment: [f64; 2] = std::array::from_fn(|i| {
            if moving {
                factor[i] * (STRENGTH[i] - SLOPE[i] * theta[i])
            } else {
                0.0
            }
        });
        let travel: [f64; 2] = std::array::from_fn(|i| {
            if moving {
                DT * (SPIN[0] - SPIN[i + 1])
                    - 0.5 * DT * DT * (w[i][0] * moment[0] + w[i][1] * moment[1])
            } else {
                0.0
            }
        });
        let residual: [f64; 2] = std::array::from_fn(|i| {
            matrix[i][0] * theta[0] + matrix[i][1] * theta[1] - rhs[i] - moment[i] * travel[i]
        });
        let jac: [[f64; 2]; 2] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                matrix[i][j] + if i == j { slope[i] * travel[i] } else { 0.0 }
                    - 0.5 * DT * DT * moment[i] * w[i][j] * slope[j]
            })
        });
        let determinant = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        assert!(determinant > 0.0);
        theta[0] -= (jac[1][1] * residual[0] - jac[0][1] * residual[1]) / determinant;
        theta[1] -= (jac[0][0] * residual[1] - jac[1][0] * residual[0]) / determinant;
    }
    let moment: [f64; 2] = std::array::from_fn(|i| {
        if moving {
            factor[i] * (STRENGTH[i] - SLOPE[i] * theta[i])
        } else {
            0.0
        }
    });
    let heat = std::array::from_fn(|i| {
        moment[i]
            * (DT * (SPIN[0] - SPIN[i + 1])
                - 0.5 * DT * DT * (w[i][0] * moment[0] + w[i][1] * moment[1]))
    });
    Reference {
        faces: theta.map(|t| t + 300.0),
        moment,
        heat,
        cells: [
            300.0 + (base[0] + d[0][0] * theta[0] + d[1][0] * theta[1]) / denom[0],
            300.0 + (base[1] + d[0][1] * theta[0]) / denom[1],
            300.0 + (base[2] + d[1][1] * theta[1]) / denom[2],
        ],
    }
}
fn energy(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum()
}
fn kinetic(bodies: &[RigidBodyState]) -> f64 {
    bodies
        .iter()
        .enumerate()
        .map(|(i, body)| {
            0.5 * (body.linear_momentum_world().norm_squared() / MASS[i]
                + body.angular_momentum_body().norm_squared() / INERTIA[i])
        })
        .sum()
}
fn compare(
    candidate: &HeatedRigidAdhesiveIslandCandidate<'_, '_, Vec<f64>>,
    island: &Island,
    before: &DepositionSimulation<'_>,
    sources: &Sources,
    moving: bool,
) {
    let expected = reference(INITIAL, sources.areas(), sources.radii(), moving);
    assert_eq!(candidate.mechanical.body_indices, [0, 1, 2]);
    for i in 0..2 {
        let native = &candidate.mechanical.report.mechanical.contacts[i];
        assert_eq!(native.patch_id, PATCH[i]);
        assert_eq!(
            native.contact.mode,
            if moving {
                AdhesiveWrenchMode::Sliding
            } else {
                AdhesiveWrenchMode::Sticking
            }
        );
        close(native.contact.force_n[0], 0.0, 2e-9);
        close(native.contact.force_n[1], 0.0, 2e-9);
        close(native.contact.spin_moment_nm, -expected.moment[i], 1e-9);
        close(native.dissipated_energy_j, expected.heat[i], 1e-9);
        let heat = candidate
            .heat
            .interfaces
            .iter()
            .find(|h| h.link_id == LINK[i])
            .unwrap();
        close(heat.interface_temperature_k, expected.faces[i], 3e-7);
        close(heat.dissipated_j, expected.heat[i], 1e-9);
        for side in 0..2 {
            let body = if side == 0 { 0 } else { i + 1 };
            let g = K[body] * sources.areas()[i] / HALF_LENGTH;
            close(heat.conductances_w_k[side], g, 1e-12);
            close(
                heat.body_heat_j[side],
                DT * g * (expected.faces[i] - expected.cells[body]),
                1e-9,
            );
        }
        let source = &candidate
            .friction_sources
            .iter()
            .find(|s| s.patch_id == PATCH[i])
            .unwrap()
            .friction_source;
        assert_eq!(source.original(), sources.patches()[i]);
        close(
            source.friction_source().state_coordinate("T").unwrap(),
            heat.interface_temperature_k,
            1e-12,
        );
        close(
            source.law().parameters().kinetic_shear_strength_pa,
            STRENGTH[i] - SLOPE[i] * (heat.interface_temperature_k - 300.0),
            2e-9,
        );
        for (old, new) in sources.shear[i]
            .resolved()
            .properties()
            .iter()
            .zip(source.friction_source().properties())
        {
            assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
        }
    }
    assert!((expected.faces[0] - expected.faces[1]).abs() > 1.0);
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected.cells[i],
            3e-7,
        );
        let change = DT
            * if i == 0 {
                -expected.moment[0] - expected.moment[1]
            } else {
                expected.moment[i - 1]
            };
        close(
            candidate.mechanical.bodies_after[i]
                .angular_momentum_body()
                .z,
            island.bodies[i].angular_momentum_body().z + change,
            2e-10,
        );
        assert_eq!(
            candidate.mechanical.bodies_after[i].pose(),
            island.bodies[i].pose()
        );
    }
    close(
        energy(&candidate.thermal) - energy(before),
        expected.heat.iter().sum(),
        2e-9,
    );
    close(
        kinetic(&candidate.mechanical.bodies_after) - kinetic(&island.bodies)
            + energy(&candidate.thermal)
            - energy(before),
        0.0,
        3e-9,
    );
    assert!(candidate.combined_energy_residual_j.abs() <= controls().combined_energy_tolerance_j);
    assert!(candidate.summed_energy_residual_j.abs() <= controls().combined_energy_tolerance_j);
    assert_eq!(candidate.heat.external_heat_j, 0.0);
    assert_eq!(candidate.thermal.time_s(), DT);
}

#[test]
fn decimal_prismatic_faces_admit_only_arithmetic_roundoff_not_physical_overcoverage() {
    let sources = Sources::new();
    let areas = sources.areas();
    let face = areas[0].max(areas[1]);
    let half: f64 = 0.05;
    let volume = face * (2.0 * half);
    let reconstructed_length = (-half - (-3.0 * half)).abs();
    // Independent reproduction of the original strict-greater-than refusal.
    assert!(face > volume / reconstructed_length);
    let admitted = thermal_with_geometry(INITIAL, areas, half, face);
    for i in 0..2 {
        let route = FrictionHeatRoute::from_link(&admitted, LINK[i], [K[0], K[i + 1]])
            .expect("roundoff-consistent original physical face");
        assert_eq!(route.contact_area_m2(), areas[i]);
    }
    let oversized = thermal_with_geometry(INITIAL, areas, half, 0.99 * face);
    assert!(FrictionHeatRoute::from_link(&oversized, LINK[0], [K[0], K[1]]).is_err());
}

#[test]
fn zero_and_tensile_adhesive_spin_heat_close_shared_inertia_source_feedback_and_finite_caloric_balance()
 {
    let sources = Sources::new();
    assert!(sources.response[0].normal_force_n.abs() < 1e-12);
    assert!(sources.response[1].normal_force_n < 0.0);
    let original = (sources.normal.clone(), sources.shear.clone());
    let island = Island::new(true);
    let saved = island.bodies;
    let thermal = thermal(INITIAL, sources.areas());
    let routes = routes(&thermal);
    let contacts = island.contacts(&sources, &routes);
    let candidate = run(
        &island,
        &contacts,
        &thermal,
        controls(),
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    compare(&candidate, &island, &thermal, &sources, true);
    assert!(candidate.nonlinear.converged);
    assert!(
        candidate
            .heat
            .interfaces
            .iter()
            .flat_map(|h| h.body_heat_j)
            .any(|q| q < 0.0)
    );
    let mut reversed = contacts;
    reversed.reverse();
    let retry = run(
        &island,
        &reversed,
        &thermal,
        controls(),
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.mechanical, candidate.mechanical);
    assert_eq!(retry.source, candidate.source);
    assert_eq!(retry.friction_sources, candidate.friction_sources);
    assert_eq!(island.bodies, saved);
    assert_eq!(sources.normal, original.0);
    assert_eq!(sources.shear, original.1);
}

#[test]
fn stationary_adhesive_disks_still_conduct_with_signed_heat_and_no_mechanical_dissipation() {
    let sources = Sources::new();
    let island = Island::new(false);
    let thermal = thermal(INITIAL, sources.areas());
    let routes = routes(&thermal);
    let candidate = run(
        &island,
        &island.contacts(&sources, &routes),
        &thermal,
        controls(),
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    compare(&candidate, &island, &thermal, &sources, false);
    assert_eq!(candidate.heat.dissipated_j, 0.0);
    assert_eq!(candidate.mechanical.bodies_after, island.bodies);
    assert!(
        candidate
            .heat
            .interfaces
            .iter()
            .all(|h| h.body_heat_j[0] * h.body_heat_j[1] < 0.0)
    );
}

#[test]
fn source_domain_geometry_and_late_cancellation_refuse_without_publishing_any_owner() {
    let sources = Sources::new();
    let original_sources = (sources.normal.clone(), sources.shear.clone());
    let island = Island::new(true);
    let saved = island.bodies;
    let thermal = thermal(INITIAL, sources.areas());
    let before_energy = energy(&thermal);
    let routes = routes(&thermal);
    let contacts = island.contacts(&sources, &routes);
    let baseline = run(
        &island,
        &contacts,
        &thermal,
        controls(),
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    let mut out_of_domain = contacts;
    out_of_domain[0].temperature_guess_k = 330.0;
    assert!(matches!(
        run(
            &island,
            &out_of_domain,
            &thermal,
            controls(),
            transport,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveIslandHeatError::FrictionSource(_))
    ));
    let smaller = sources.areas().map(|a| 0.9 * a);
    let wrong_thermal = self::thermal(INITIAL, smaller);
    let wrong_routes = self::routes(&wrong_thermal);
    assert!(matches!(
        run(
            &island,
            &island.contacts(&sources, &wrong_routes),
            &wrong_thermal,
            controls(),
            transport,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveIslandHeatError::Invalid(_))
    ));
    let calls = Cell::new(0usize);
    let gate = CancelGate::new();
    let cancelled = run(
        &island,
        &contacts,
        &thermal,
        controls(),
        |state, ids, temperatures| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.trials {
                gate.request();
            }
            transport(state, ids, temperatures)
        },
        &gate,
    );
    assert!(matches!(
        cancelled,
        Err(RigidAdhesiveIslandHeatError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(calls.get(), baseline.trials);
    let retry = run(
        &island,
        &contacts,
        &thermal,
        controls(),
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.mechanical, baseline.mechanical);
    assert_eq!(retry.friction_sources, baseline.friction_sources);
    assert_eq!(island.bodies, saved);
    assert_eq!(sources.normal, original_sources.0);
    assert_eq!(sources.shear, original_sources.1);
    assert_eq!(energy(&thermal), before_energy);
    assert_eq!(thermal.time_s(), 0.0);
}
