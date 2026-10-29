//! G1/G3/G4: original full articulated response and one finite caloric network.
//! Curves and caloric data below are explicit synthetic fixtures, not measurements.

use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    SourceThermalCoulomb, SourceThermalCoulombState, ThermalFrictionAxes,
};
use fs_couple::friction_heat::articulated::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{Mat3, Quat, Se3, So3, Twist, Vec3};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::tangential::planar::ArticulatedPlanarTangentialContact;
use fs_mbd::articulated_contacts::tangential::planar::multiple::ArticulatedMultiTangentialConfig;
use fs_mbd::articulated_contacts::tangential::planar::transport::ArticulatedFramedTangentialSite;
use fs_mbd::articulated_contacts::tangential::{
    ArticulatedTangentialConfig, ArticulatedTangentialDrive,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
};
use fs_qty::{Dims, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::tangential::multiple::MultipleTangentialConfig;
use fs_tribo::tangential::transport::{
    FramedTangentialState, ObjectiveTangentialLaw, TangentialTransportConfig,
};
use fs_tribo::{
    FrictionRegime, InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 0.01;
const A: f64 = 0.01;
const AREA: f64 = 0.001;
const MASS: f64 = 2.0;
const MOMENT: f64 = 4.0 * A * A;
const CAP: [f64; 3] = [0.03, 0.02, 0.04];
const K: [f64; 3] = [60.0, 40.0, 100.0];
const MU: [f64; 2] = [0.4, 0.35];
const SLOPE: [f64; 2] = [0.01, 0.015];
const SUPPORT: [f64; 2] = [0.2, -0.1];
const LINK: [u64; 2] = [20, 10];

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-9 * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn close_vec(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn world(frame: TangentFrame, v: [f64; 2]) -> Vec3 {
    vector(frame.longitudinal()).scale(v[0]) + vector(frame.lateral()).scale(v[1])
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

fn thermal(temperatures: [f64; 3], rotation: So3) -> DepositionSimulation<'static> {
    let endpoints = [[-A, A], [-3.0 * A, -A], [A, 3.0 * A]];
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAP[i],
                specific_h_j_kg: temperatures[i] - 300.0,
                geometry: CellGeometry {
                    start_m: array(
                        rotation
                            .rotate(Vec3::new(0.0, 0.0, endpoints[i][0]))
                            .unwrap(),
                    ),
                    end_m: array(
                        rotation
                            .rotate(Vec3::new(0.0, 0.0, endpoints[i][1]))
                            .unwrap(),
                    ),
                    volume_m3: AREA * 2.0 * A,
                    surface_area_m2: 2.0 * AREA + 8.0 * AREA.sqrt() * A,
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
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(AREA, [A; 2], [K[0], K[i + 1]]).unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
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
fn bindings(routes: &[FrictionHeatRoute; 2]) -> [ArticulatedFrictionHeatRoute<'_>; 2] {
    std::array::from_fn(|i| ArticulatedFrictionHeatRoute {
        contact_index: i,
        route: &routes[i],
        temperature_guess_k: 300.0,
    })
}

fn curve(name: &str, i: usize, values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, Dims::NONE),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (320.0, values[1])],
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("synthetic articulated thermal fixture {i}/{name}; not measured data"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn law(i: usize) -> SourceThermalCoulomb {
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic fixture".into(),
            revision: 0,
        },
        texture_frame: format!("texture-{name}"),
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(curve("static", i, [0.9, 0.7])).unwrap();
    claims
        .insert_claim(curve(
            "kinetic",
            i,
            [MU[i] + 20.0 * SLOPE[i], MU[i] - 20.0 * SLOPE[i]],
        ))
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        surface("middle"),
        surface(&format!("outer-{i}")),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        QueryPoint::new().with("T", 300.0).unwrap(),
        requirement("static"),
        requirement("kinetic"),
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap()
}
fn model() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "middle",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            MASS,
            Vec3::default(),
            Mat3 {
                m: [MOMENT, 0.0, 0.0, 0.0, MOMENT, 0.0, 0.0, 0.0, 0.5 * MOMENT],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}
fn system(model: &ArticulatedModel, rotation: So3, moving: bool) -> ArticulatedContactSystem<'_> {
    ArticulatedContactSystem {
        model,
        base: FreeFloatingBaseState::new(
            Se3::from_parts(rotation, Vec3::default()).unwrap(),
            if moving {
                Twist::new(Vec3::new(0.0, 50.0, 0.0), Vec3::new(4.0, 0.0, 0.0))
            } else {
                Twist::default()
            },
        ),
        joint_position: &[],
        joint_velocity: &[],
    }
}
fn contacts(rotation: So3, moving: bool) -> [ArticulatedPlanarTangentialContact; 2] {
    std::array::from_fn(|i| {
        let point = Vec3::new(0.0, 0.0, [-A, A][i]);
        ArticulatedPlanarTangentialContact {
            point: LinkContactPoint {
                system_index: 0,
                link_index: 0,
                point_body: point,
            },
            other: ContactCounterpart::PrescribedPlane {
                point_world: rotation.rotate(point).unwrap(),
                velocity_world: rotation
                    .rotate(Vec3::new(if moving { SUPPORT[i] } else { 0.0 }, 0.0, 0.0))
                    .unwrap(),
            },
            frame_world: TangentFrame::new(
                array(
                    rotation
                        .rotate(Vec3::new(0.0, 0.0, [1.0, -1.0][i]))
                        .unwrap(),
                ),
                array(rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap()),
            )
            .unwrap(),
        }
    })
}
fn drive() -> ArticulatedTangentialDrive {
    ArticulatedTangentialDrive {
        dt_s: DT,
        normal_force_n: 10.0,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
    }
}
fn sites<'a, L: ObjectiveTangentialLaw>(
    contacts: &[ArticulatedPlanarTangentialContact; 2],
    laws: &'a [L; 2],
    history: &'a [FramedTangentialState<L::State>; 2],
) -> [ArticulatedFramedTangentialSite<'a, L>; 2] {
    std::array::from_fn(|i| ArticulatedFramedTangentialSite {
        contact: contacts[i],
        law: &laws[i],
        history: &history[i],
        drive: drive(),
        common_spin_angle_rad: 0.0,
    })
}
fn histories<L: ObjectiveTangentialLaw>(
    contacts: &[ArticulatedPlanarTangentialContact; 2],
    laws: &[L; 2],
) -> [FramedTangentialState<L::State>; 2] {
    std::array::from_fn(|i| {
        FramedTangentialState::new(&laws[i], laws[i].initial_state(), contacts[i].frame_world)
            .unwrap()
    })
}
fn controls() -> ArticulatedFrictionHeatControls {
    ArticulatedFrictionHeatControls {
        mechanical: ArticulatedMultiTangentialConfig {
            mechanical: ArticulatedTangentialConfig {
                solve: TangentialSolveConfig {
                    max_iterations: 96,
                    maximum_force_n: 100.0,
                    force_absolute_tolerance_n: 1e-12,
                    force_relative_tolerance: 1e-12,
                },
                maximum_step_s: 0.1,
                maximum_slip_speed_m_s: 10.0,
                geometry_tolerance_m: 1e-10,
                displacement_tolerance_m: 1e-11,
                velocity_tolerance_m_s: 1e-10,
                energy_tolerance_j: 1e-10,
            },
            coupling: MultipleTangentialConfig {
                maximum_contacts: 2,
                maximum_matrix_terms: 16,
                maximum_sweeps: 64,
                maximum_law_evaluations: 4096,
                local: MatrixTangentialConfig {
                    maximum_evaluations: 256,
                    displacement_tolerance_m: 1e-13,
                    energy_tolerance_j: 1e-12,
                    matrix_relative_tolerance: 1e-12,
                },
                total_energy_tolerance_j: 1e-11,
            },
        },
        transport: TangentialTransportConfig {
            absolute_energy_tolerance_j: 1e-13,
            relative_energy_tolerance: 1e-12,
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
fn source(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<ArticulatedFrictionHeatSource<Vec<f64>>, &'static str> {
    assert_eq!(ids, &[10, 20]);
    Ok(ArticulatedFrictionHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| ArticulatedInterfaceHeatSource {
                link_id,
                conductivity_w_m_k: [K[0], K[if link_id == 20 { 1 } else { 2 }]],
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperatures.to_vec(),
        },
    })
}

// Independent two-temperature Newton solution: eliminate the three finite cells,
// use the rigid Gram response W=[[.75,.25],[.25,.75]], then analytic dry heat.
fn reference() -> ([f64; 2], [f64; 3], [f64; 2], [f64; 2]) {
    let d = K.map(|k| DT * k * AREA / A);
    let cross = -d[0] * d[0] / (CAP[0] + 2.0 * d[0]);
    let matrix = [
        [d[0] + d[1] * CAP[1] / (CAP[1] + d[1]) + cross, cross],
        [cross, d[0] + d[2] * CAP[2] / (CAP[2] + d[2]) + cross],
    ];
    let w = [[0.75, 0.25], [0.25, 0.75]];
    let speed = [3.5 - SUPPORT[0], 4.5 - SUPPORT[1]];
    let slopes = SLOPE.map(|x| 10.0 * x);
    let mut theta = [0.0; 2];
    let mut heat = [0.0; 2];
    for _ in 0..12 {
        let c: [f64; 2] = std::array::from_fn(|i| 10.0 * MU[i] - slopes[i] * theta[i]);
        let mut residual = [0.0; 2];
        let mut jac = matrix;
        for i in 0..2 {
            let j = 1 - i;
            heat[i] = c[i] * (DT * speed[i] - 0.5 * DT * DT * (w[i][i] * c[i] + w[i][j] * c[j]));
            residual[i] = matrix[i][0] * theta[0] + matrix[i][1] * theta[1] - heat[i];
            jac[i][i] += slopes[i]
                * (DT * speed[i] - DT * DT * w[i][i] * c[i] - 0.5 * DT * DT * w[i][j] * c[j]);
            jac[i][j] -= 0.5 * DT * DT * w[i][j] * c[i] * slopes[j];
        }
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        theta[0] -= (jac[1][1] * residual[0] - jac[0][1] * residual[1]) / det;
        theta[1] -= (jac[0][0] * residual[1] - jac[1][0] * residual[0]) / det;
    }
    let cell = [
        d[0] * (theta[0] + theta[1]) / (CAP[0] + 2.0 * d[0]),
        d[1] * theta[0] / (CAP[1] + d[1]),
        d[2] * theta[1] / (CAP[2] + d[2]),
    ];
    (
        theta.map(|x| 300.0 + x),
        cell.map(|x| 300.0 + x),
        heat,
        std::array::from_fn(|i| 10.0 * MU[i] - slopes[i] * theta[i]),
    )
}
fn energy(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
        .sum()
}
fn audit<S, P>(
    candidate: &HeatedArticulatedFrictionCandidate<'_, S, P>,
    before: &DepositionSimulation<'_>,
) {
    let m = &candidate.mechanical.mechanics;
    close(candidate.combined_energy_residual_j, 0.0);
    close(candidate.heat.total_energy_residual_j, 0.0);
    close(candidate.mechanical.accepted_energy_balance_residual_j, 0.0);
    let shear: f64 = m
        .contact_energy_after_j
        .iter()
        .zip(&m.contact_energy_before_j)
        .map(|(a, b)| a - b)
        .sum();
    let transport: f64 = candidate
        .mechanical
        .transport_energy_residuals_j
        .iter()
        .sum();
    close(
        m.kinetic_energy_change_j + shear + transport + energy(&candidate.thermal)
            - energy(before)
            - m.support_work_j
            - candidate.heat.external_heat_j,
        0.0,
    );
    for (i, heat) in candidate.heat.interfaces.iter().enumerate() {
        close(heat.body_heat_j.iter().sum(), heat.dissipated_j);
        close(
            heat.dissipated_j,
            m.contacts[candidate.contact_indices[i]].dissipated_j,
        );
    }
}

#[test]
fn original_temperature_claims_and_shared_body_cell_match_independent_joint_reference() {
    let model = model();
    let body = system(&model, So3::identity(), true);
    let contact = contacts(So3::identity(), true);
    let laws = [law(0), law(1)];
    let history = histories(&contact, &laws);
    let saved = history.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let candidate = propose_articulated_friction_heat(
        &[body],
        &sites(&contact, &laws, &history),
        &thermal,
        0.0,
        &bindings(&routes),
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::default(),
    )
    .unwrap();
    let (surface, cells, heat, force) = reference();
    assert_eq!(candidate.contact_indices, [1, 0]);
    assert_eq!(
        candidate.source,
        candidate
            .heat
            .interfaces
            .iter()
            .map(|h| h.interface_temperature_k)
            .collect::<Vec<_>>()
    );
    assert!(candidate.trials > 2);
    assert!(candidate.nonlinear.converged);
    for i in 0..2 {
        let h = &candidate.heat.interfaces[1 - i];
        let m = &candidate.mechanical.mechanics.contacts[i];
        close(h.interface_temperature_k, surface[i]);
        close(h.dissipated_j, heat[i]);
        close(m.traction_n[0], -force[i]);
        close(m.traction_n[1], 0.0);
        assert_eq!(m.regime, FrictionRegime::Sliding);
        assert!(force[i] < 10.0 * MU[i] - 0.1);
        close(
            m.state
                .state()
                .constitutive_state()
                .contact_temperature_k()
                .unwrap(),
            surface[i],
        );
        for property in m.state.state().friction_source().properties() {
            let original = laws[i]
                .initial_source()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            laws[i]
                .interface()
                .card()
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            cells[i],
        );
    }
    let response = &candidate.mechanical.mechanics.responses[0];
    close(
        response.base_after.twist_body.linear.x,
        4.0 - DT * (force[0] + force[1]) / MASS,
    );
    close(
        response.base_after.twist_body.angular.y,
        50.0 + DT * A * (force[0] - force[1]) / MOMENT,
    );
    close(candidate.mechanical.mechanics.mobility_kg_inverse[2], 0.25);
    audit(&candidate, &thermal);
    let entropy = candidate.heat.entropy.unwrap();
    assert!(entropy.friction_generation_j_k > 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(history, saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(body.base.twist_body.linear.x, 4.0);
    assert_eq!(candidate.target_time_s, DT);
    assert_eq!(candidate.thermal.time_s(), DT);
}

#[test]
fn world_rotation_tangent_gauge_and_binding_order_preserve_actual_heat_and_recoil() {
    let model = model();
    let laws = [law(0), law(1)];
    let rotations = [
        So3::identity(),
        So3::try_from_quat(Quat {
            w: 0.5,
            x: 0.5,
            y: 0.5,
            z: 0.5,
        })
        .unwrap(),
    ];
    let mut baseline: Option<
        HeatedArticulatedFrictionCandidate<'_, SourceThermalCoulombState, Vec<f64>>,
    > = None;
    for (pass, rotation) in rotations.into_iter().enumerate() {
        let body = system(&model, rotation, true);
        let accepted_contact = contacts(rotation, true);
        let mut history = histories(&accepted_contact, &laws);
        for i in 0..2 {
            let old = laws[i]
                .advance(
                    history[i].state(),
                    &TangentialStepInput {
                        normal_n: 10.0,
                        slip_increment_m: [0.003, -0.004],
                        dt_s: DT,
                        contact_temperature_k: Some(300.0),
                        nominal_pressure_pa: None,
                        compliance_m_per_n: 0.0,
                        solve: controls().mechanical.mechanical.solve,
                    },
                )
                .unwrap()
                .state;
            history[i] =
                FramedTangentialState::new(&laws[i], old, accepted_contact[i].frame_world).unwrap();
        }
        let saved = history.clone();
        let mut contact = accepted_contact;
        if pass == 1 {
            for i in 0..2 {
                contact[i].frame_world = contact[i].frame_world.rotated([0.43, -0.27][i]).unwrap();
            }
        }
        let thermal = thermal([300.0; 3], rotation);
        let routes = routes(&thermal);
        let mut bound = bindings(&routes);
        if pass == 1 {
            bound.reverse();
        }
        let candidate = propose_articulated_friction_heat(
            &[body],
            &sites(&contact, &laws, &history),
            &thermal,
            0.0,
            &bound,
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::default(),
        )
        .unwrap();
        if let Some(reference) = &baseline {
            let reference: &HeatedArticulatedFrictionCandidate<'_, _, Vec<f64>> = reference;
            close_vec(
                candidate.mechanical.mechanics.responses[0]
                    .base_after
                    .twist_body
                    .linear,
                reference.mechanical.mechanics.responses[0]
                    .base_after
                    .twist_body
                    .linear,
            );
            close_vec(
                candidate.mechanical.mechanics.responses[0]
                    .base_after
                    .twist_body
                    .angular,
                reference.mechanical.mechanics.responses[0]
                    .base_after
                    .twist_body
                    .angular,
            );
            for i in 0..2 {
                close_vec(
                    world(
                        contact[i].frame_world,
                        candidate.mechanical.mechanics.contacts[i].traction_n,
                    ),
                    rotation
                        .rotate(world(
                            contacts(So3::identity(), true)[i].frame_world,
                            reference.mechanical.mechanics.contacts[i].traction_n,
                        ))
                        .unwrap(),
                );
                close_vec(
                    world(
                        contact[i].frame_world,
                        candidate.mechanical.mechanics.contacts[i]
                            .state
                            .state()
                            .constitutive_state()
                            .slip_m(),
                    ),
                    rotation
                        .rotate(world(
                            contacts(So3::identity(), true)[i].frame_world,
                            reference.mechanical.mechanics.contacts[i]
                                .state
                                .state()
                                .constitutive_state()
                                .slip_m(),
                        ))
                        .unwrap(),
                );
                close(
                    candidate.heat.interfaces[i].interface_temperature_k,
                    reference.heat.interfaces[i].interface_temperature_k,
                );
                close(
                    candidate.heat.interfaces[i].dissipated_j,
                    reference.heat.interfaces[i].dissipated_j,
                );
            }
        }
        audit(&candidate, &thermal);
        assert_eq!(history, saved);
        baseline = Some(candidate);
    }
}

fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, p);
        b.swap(k, p);
        for i in k + 1..n {
            let r = a[i][k] / a[k][k];
            for j in k..n {
                a[i][j] -= r * a[k][j];
            }
            b[i] -= r * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = (b[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}

#[test]
fn stationary_loaded_contacts_retain_pure_shared_cell_conduction_and_entropy() {
    let model = model();
    let body = system(&model, So3::identity(), false);
    let contact = contacts(So3::identity(), false);
    let laws = [law(0), law(1)];
    let history = histories(&contact, &laws);
    let old = [310.0, 290.0, 300.0];
    let thermal = thermal(old, So3::identity());
    let routes = routes(&thermal);
    let candidate = propose_articulated_friction_heat(
        &[body],
        &sites(&contact, &laws, &history),
        &thermal,
        0.0,
        &bindings(&routes),
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::default(),
    )
    .unwrap();
    let g = [1, 2].map(|i| DT * AREA / (A / K[0] + A / K[i]));
    let reference = solve_dense(
        vec![
            vec![CAP[0] + g[0] + g[1], -g[0], -g[1]],
            vec![-g[0], CAP[1] + g[0], 0.0],
            vec![-g[1], 0.0, CAP[2] + g[1]],
        ],
        (0..3).map(|i| CAP[i] * old[i]).collect(),
    );
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            reference[i],
        );
    }
    for i in 0..2 {
        let h = &candidate.heat.interfaces[1 - i];
        close(
            h.interface_temperature_k,
            (K[0] * reference[0] + K[i + 1] * reference[i + 1]) / (K[0] + K[i + 1]),
        );
        assert!(h.body_heat_j[0] < 0.0);
        assert!(h.body_heat_j[1] > 0.0);
        assert_eq!(h.dissipated_j, 0.0);
    }
    assert_eq!(
        candidate.mechanical.mechanics.responses[0].base_after,
        body.base
    );
    let entropy = candidate.heat.entropy.unwrap();
    assert_eq!(entropy.friction_generation_j_k, 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert_eq!(entropy.network_conduction_generation_j_k, 0.0);
    close(energy(&candidate.thermal), energy(&thermal));
    audit(&candidate, &thermal);
}

#[test]
fn topology_geometry_source_and_temperature_ownership_refuse_before_candidate_publication() {
    let model = model();
    let body = system(&model, So3::identity(), true);
    let contact = contacts(So3::identity(), true);
    let laws = [law(0), law(1)];
    let history = histories(&contact, &laws);
    let saved = history.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let bound = bindings(&routes);
    let call = |site: &[ArticulatedFramedTangentialSite<'_, SourceThermalCoulomb>],
                bindings: &[ArticulatedFrictionHeatRoute<'_>],
                clock| {
        propose_articulated_friction_heat(
            &[body],
            site,
            &thermal,
            clock,
            bindings,
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::default(),
        )
    };
    let original = sites(&contact, &laws, &history);
    assert!(matches!(
        call(&original, &bound[..1], 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut bad = bound;
    bad[1].contact_index = 0;
    assert!(matches!(
        call(&original, &bad, 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut bad = bound;
    bad[1].route = bad[0].route;
    assert!(matches!(
        call(&original, &bad, 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    assert!(matches!(
        call(&original, &bound, DT),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut bad = sites(&contact, &laws, &history);
    bad[1].drive.contact_temperature_k = Some(300.0);
    assert!(matches!(
        call(&bad, &bound, 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut bad_contact = contact;
    bad_contact[0].point.point_body.x = 0.001;
    if let ContactCounterpart::PrescribedPlane { point_world, .. } = &mut bad_contact[0].other {
        point_world.x = 0.001;
    }
    assert!(matches!(
        call(&sites(&bad_contact, &laws, &history), &bound, 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut bad_contact = contact;
    bad_contact[0].frame_world = TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    assert!(matches!(
        call(&sites(&bad_contact, &laws, &history), &bound, 0.0),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let mut wrong_source = sites(&contact, &laws, &history);
    wrong_source[0].law = &laws[1];
    assert!(matches!(
        call(&wrong_source, &bound, 0.0),
        Err(ArticulatedFrictionHeatError::Mechanical(_))
    ));
    let mut outside = bound;
    outside[0].temperature_guess_k = 330.0;
    assert!(matches!(
        call(&original, &outside, 0.0),
        Err(ArticulatedFrictionHeatError::Mechanical(_))
    ));
    let wrong_ids = propose_articulated_friction_heat(
        &[body],
        &original,
        &thermal,
        0.0,
        &bound,
        controls(),
        &[0.0; 3],
        |t, ids, temps| {
            let mut s = source(t, ids, temps)?;
            s.interfaces.swap(0, 1);
            Ok(s)
        },
        &CancelGate::default(),
    );
    assert!(matches!(
        wrong_ids,
        Err(ArticulatedFrictionHeatError::<&'static str>::Invalid(_))
    ));
    assert_eq!(history, saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(energy(&thermal), 0.0);
}

#[test]
fn final_source_failure_cancellation_and_trial_budget_leave_exact_loaded_retry() {
    let model = model();
    let body = system(&model, So3::identity(), true);
    let contact = contacts(So3::identity(), true);
    let laws = [law(0), law(1)];
    let history = histories(&contact, &laws);
    let saved = history.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let bound = bindings(&routes);
    let site = sites(&contact, &laws, &history);
    let calls = Cell::new(0usize);
    let baseline = propose_articulated_friction_heat(
        &[body],
        &site,
        &thermal,
        0.0,
        &bound,
        controls(),
        &[0.0; 3],
        |t, ids, temps| {
            calls.set(calls.get() + 1);
            source(t, ids, temps)
        },
        &CancelGate::default(),
    )
    .unwrap();
    assert_eq!(calls.get(), baseline.trials);
    assert!(baseline.trials > 2);
    let calls = Cell::new(0usize);
    let failed = propose_articulated_friction_heat(
        &[body],
        &site,
        &thermal,
        0.0,
        &bound,
        controls(),
        &[0.0; 3],
        |t, ids, temps| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.trials {
                Err("late original source")
            } else {
                source(t, ids, temps)
            }
        },
        &CancelGate::default(),
    );
    assert!(matches!(
        failed,
        Err(ArticulatedFrictionHeatError::Source("late original source"))
    ));
    assert_eq!(calls.get(), baseline.trials);
    let calls = Cell::new(0usize);
    let gate = CancelGate::default();
    let cancelled = propose_articulated_friction_heat(
        &[body],
        &site,
        &thermal,
        0.0,
        &bound,
        controls(),
        &[0.0; 3],
        |t, ids, temps| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.trials {
                gate.request();
            }
            source(t, ids, temps)
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), baseline.trials);
    let limited = ArticulatedFrictionHeatControls {
        maximum_trials: baseline.trials - 1,
        ..controls()
    };
    assert!(matches!(
        propose_articulated_friction_heat(
            &[body],
            &site,
            &thermal,
            0.0,
            &bound,
            limited,
            &[0.0; 3],
            source,
            &CancelGate::default()
        ),
        Err(ArticulatedFrictionHeatError::TrialBudget)
    ));
    let bounded = ArticulatedFrictionHeatControls {
        maximum_unknowns: 4,
        ..controls()
    };
    assert!(matches!(
        propose_articulated_friction_heat(
            &[body],
            &site,
            &thermal,
            0.0,
            &bound,
            bounded,
            &[0.0; 3],
            source,
            &CancelGate::default()
        ),
        Err(ArticulatedFrictionHeatError::Invalid(_))
    ));
    let retry = propose_articulated_friction_heat(
        &[body],
        &site,
        &thermal,
        0.0,
        &bound,
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::default(),
    )
    .unwrap();
    assert_eq!(
        retry.mechanical.mechanics.contacts,
        baseline.mechanical.mechanics.contacts
    );
    assert_eq!(
        retry.mechanical.mechanics.responses[0].base_after,
        baseline.mechanical.mechanics.responses[0].base_after
    );
    assert_eq!(retry.source, baseline.source);
    for id in thermal.cells().keys() {
        assert_eq!(
            retry.thermal.temperature_k(*id).unwrap(),
            baseline.thermal.temperature_k(*id).unwrap()
        );
    }
    assert_eq!(history, saved);
    assert_eq!(thermal.time_s(), 0.0);
    audit(&retry, &thermal);
}

#[test]
fn retained_elastic_release_stays_mechanical_and_frozen_plastic_temperature_refuses() {
    let model = model();
    let body = system(&model, So3::identity(), false);
    let contact = contacts(So3::identity(), false);
    let laws = [
        ThirdBodyLayerContact::new(400.0, 100.0, 0.4, 300.0, 1.0).unwrap(),
        ThirdBodyLayerContact::new(700.0, 100.0, 0.7, 300.0, 1.0).unwrap(),
    ];
    let history = std::array::from_fn(|i| {
        let pre = laws[i]
            .advance(
                &laws[i].initial_state(),
                &TangentialStepInput {
                    normal_n: 10.0,
                    slip_increment_m: [[0.004, -0.002], [-0.003, 0.005]][i],
                    dt_s: DT,
                    contact_temperature_k: Some(300.0),
                    nominal_pressure_pa: None,
                    compliance_m_per_n: 0.0,
                    solve: controls().mechanical.mechanical.solve,
                },
            )
            .unwrap()
            .state;
        FramedTangentialState::new(&laws[i], pre, contact[i].frame_world).unwrap()
    });
    let saved = history.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let mut site = sites(&contact, &laws, &history);
    site[0].common_spin_angle_rad = 0.3;
    site[1].common_spin_angle_rad = -0.2;
    let candidate = propose_articulated_friction_heat(
        &[body],
        &site,
        &thermal,
        0.0,
        &bindings(&routes),
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::default(),
    )
    .unwrap();
    assert!(candidate.mechanical.mechanics.kinetic_energy_change_j > 0.0);
    assert_eq!(candidate.heat.dissipated_j, 0.0);
    assert_eq!(energy(&candidate.thermal), 0.0);
    assert!(
        candidate
            .mechanical
            .mechanics
            .contacts
            .iter()
            .all(|step| step.stored_change_j < 0.0)
    );
    audit(&candidate, &thermal);
    assert_eq!(history, saved);
    // A native frozen-temperature plastic law must not be silently continued
    // through the hotter interface found by a thermal root iteration.
    let plastic = [
        ThirdBodyLayerContact::new(400.0, 0.5, 0.4, 300.0, 1.0).unwrap(),
        ThirdBodyLayerContact::new(700.0, 0.7, 0.7, 300.0, 1.0).unwrap(),
    ];
    let contact = contacts(So3::identity(), true);
    let state = histories(&contact, &plastic);
    let saved = state.clone();
    let result = propose_articulated_friction_heat(
        &[system(&model, So3::identity(), true)],
        &sites(&contact, &plastic, &state),
        &thermal,
        0.0,
        &bindings(&routes),
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::default(),
    );
    assert!(matches!(
        result,
        Err(ArticulatedFrictionHeatError::Mechanical(_))
    ));
    assert_eq!(state, saved);
    assert_eq!(thermal.time_s(), 0.0);
}
