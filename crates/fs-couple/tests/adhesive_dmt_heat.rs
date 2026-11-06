//! Actual driven DMT body recoil and the original finite-cell transport owner.
//! The material, caloric and conductivity inputs are manufactured test data.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::{SourceDmtSphereContact, coupled::DmtCoupledState};
use fs_couple::friction_heat::adhesive::dmt::*;
use fs_couple::friction_heat::cohesive::CohesiveHeatError;
use fs_couple::friction_heat::network::FrictionHeatNetworkStep;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::normal_motion::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 5e-10;
const LENGTH: f64 = 1e-6;
const MASS: [f64; 2] = [2e-15, 3e-15];
const CP: f64 = 1000.0;
const K: [f64; 2] = [2.0, 1.0];
const KINETIC: f64 = 2e6;
const LINK: u64 = 37;
const BODY_IDS: [u64; 2] = [90, 40];
const MATERIAL_IDS: [MaterialId; 2] = [MaterialId(10), MaterialId(20)];

fn close(a: f64, b: f64, absolute: f64) {
    common::close(a, b, absolute, 5e-9);
}

fn controls() -> RigidDmtHeatControls {
    let mut contact = common::config();
    contact.coordinate_absolute_tolerance_m = 1e-23;
    contact.coordinate_relative_tolerance = 2e-12;
    contact.force_absolute_tolerance_n = 1e-20;
    contact.force_relative_tolerance = 1e-10;
    RigidDmtHeatControls {
        mechanical: RigidDmtControls {
            contact,
            maximum_bodies: 3,
            geometry_tolerance_m: 1e-20,
            displacement_tolerance_m: 2e-21,
            velocity_tolerance_m_s: 1e-12,
            angular_velocity_tolerance_rad_s: 1e-6,
            linear_momentum_tolerance_ns: 1e-29,
            angular_momentum_tolerance_nms: 1e-34,
            energy_tolerance_j: 2e-26,
        },
        thermal: RigidDmtThermalControls {
            maximum_cells: 2,
            maximum_interfaces: 1,
            maximum_unknowns: 3,
            maximum_trials: 128,
            maximum_newton_iterations: 12,
            difference_step_j: 1e-21,
            temperature_scale_j_per_k: 1e-18,
            heat: FrictionHeatTolerances {
                cell_j: 2e-26,
                thermal_j: 4e-26,
                interface_j: 2e-26,
                total_j: 4e-26,
            },
            combined_energy_tolerance_j: 8e-26,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-26,
                relative_tolerance: 1e-10,
                linear_restart: 3,
                forcing_minimum: 1e-6,
                ..Default::default()
            },
            audit_entropy: true,
        },
        maximum_fractional_area_change: 0.1,
    }
}

struct Fixture {
    normal: SourceDmtSphereContact,
    shear: SourceAdhesiveShear,
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
    approach: f64,
}
impl Fixture {
    fn new(approach: f64, velocity: Vec3) -> Self {
        let (normal, shear) = common::sources(0.05, 1.5 * KINETIC, KINETIC);
        let a =
            MassProperties::new(MASS[0], Vec3::ZERO, Vec3::new(0.8e-27, 1.1e-27, 1.4e-27)).unwrap();
        let b =
            MassProperties::new(MASS[1], Vec3::ZERO, Vec3::new(1.5e-27, 1.7e-27, 2e-27)).unwrap();
        let body = |position, momentum| {
            RigidBodyState::new(
                Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
                momentum,
                Vec3::ZERO,
            )
            .unwrap()
        };
        Self {
            normal,
            shear,
            // The contact's A/B order is 2/0; body 1 must remain untouched.
            bodies: [
                body(Vec3::new(0.0, 0.0, LENGTH), Vec3::ZERO),
                body(Vec3::new(5.0 * LENGTH, 0.0, 0.0), Vec3::ZERO),
                body(Vec3::new(0.0, 0.0, -LENGTH), velocity.scale(MASS[0])),
            ],
            properties: [b, a, a],
            approach,
        }
    }
    fn contact(&self) -> RigidDmtContact<'_> {
        RigidDmtContact {
            binding: RigidContactBinding {
                body_a: 2,
                body_b: 0,
                arm_a_body: Vec3::new(0.0, 0.0, LENGTH + 0.5 * self.approach),
                arm_b_body: Vec3::new(0.0, 0.0, -LENGTH - 0.5 * self.approach),
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            reference_approach_m: 0.0,
            normal: &self.normal,
            shear: &self.shear,
        }
    }
    fn prepared(&self) -> PreparedRigidDmtPair<'_> {
        PreparedRigidDmtPair::new(
            &self.bodies,
            &self.properties,
            self.contact(),
            DT,
            controls().mechanical,
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn history(&self) -> DmtCoupledState {
        self.prepared().initialize([0.0; 2]).unwrap()
    }
    fn area(&self) -> f64 {
        self.normal
            .evaluate(self.history().position_m()[0])
            .unwrap()
            .contact_area_m2
    }
    fn balanced_normal(&self) -> [Wrench; 2] {
        let force = self
            .normal
            .evaluate(self.history().position_m()[0])
            .unwrap()
            .normal_force_n;
        [force, -force].map(|z| Wrench {
            force_world: Vec3::new(0.0, 0.0, z),
            torque_body: Vec3::ZERO,
        })
    }
    fn heat_contact<'a>(&'a self, route: &'a FrictionHeatRoute) -> RigidDmtHeatContact<'a> {
        let materials = self
            .normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        RigidDmtHeatContact {
            contact: self.contact(),
            bodies: std::array::from_fn(|i| RigidDmtHeatBody {
                body_id: BODY_IDS[i],
                thermal_material: MATERIAL_IDS[i],
                mechanical_material: &materials[i],
            }),
            route,
            temperature_guess_k: 300.0,
        }
    }
}

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / CP;
        if t.is_finite() && (250.0..=350.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        CP.recip()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        from: f64,
        to: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(from)?;
        self.temperature_k(to)?;
        Ok(Some(CP * ((to - from) / (CP * t)).ln_1p()))
    }
}
static CALORIC: Caloric = Caloric;

fn thermal(area: f64, initial: [f64; 2]) -> DepositionSimulation<'static> {
    thermal_with_ends(area, initial, [[-2.0 * LENGTH, 0.0], [0.0, 2.0 * LENGTH]])
}

fn thermal_with_ends(
    area: f64,
    initial: [f64; 2],
    ends: [[f64; 2]; 2],
) -> DepositionSimulation<'static> {
    from_events(thermal_events(area, initial, ends))
}

fn thermal_events(area: f64, initial: [f64; 2], ends: [[f64; 2]; 2]) -> Vec<TimedEvent> {
    let mut events: Vec<_> = (0..2)
        .map(|i| {
            let volume = MASS[i] / 1000.0;
            let face = volume / (2.0 * LENGTH);
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: i as u64,
                    body_id: BODY_IDS[i],
                    material_id: MATERIAL_IDS[i],
                    mass_kg: MASS[i],
                    specific_h_j_kg: CP * (initial[i] - 300.0),
                    geometry: CellGeometry {
                        start_m: [0.0, 0.0, ends[i][0]],
                        end_m: [0.0, 0.0, ends[i][1]],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * face + 8.0 * face.sqrt() * LENGTH,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: series_conductance_w_k(area, [LENGTH; 2], K).unwrap(),
        }),
    });
    events
}

fn from_events(events: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[
            (MATERIAL_IDS[0], &CALORIC as &dyn ThermalLaw),
            (MATERIAL_IDS[1], &CALORIC as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 6,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-26,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn transport(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<RigidDmtHeatSource<f64>, &'static str> {
    assert_eq!(ids, &[LINK]);
    Ok(RigidDmtHeatSource {
        interfaces: vec![RigidDmtHeatConductivity {
            link_id: LINK,
            conductivity_w_m_k: K,
        }],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperatures[0],
        },
    })
}

fn run(
    fixture: &Fixture,
    thermal: &DepositionSimulation<'static>,
    route: &FrictionHeatRoute,
    wrenches: [Wrench; 2],
    controls: RigidDmtHeatControls,
) -> Result<HeatedRigidDmtCandidate<'static, f64>, RigidDmtHeatError<&'static str>> {
    propose_rigid_dmt_heat(
        &fixture.bodies,
        &fixture.properties,
        fixture.heat_contact(route),
        fixture.history(),
        wrenches,
        thermal,
        0.0,
        DT,
        controls,
        &[0.0; 2],
        transport,
        &CancelGate::new(),
    )
}

#[test]
fn tensile_real_contact_deposits_actual_sliding_loss_and_matches_analytic_finite_bodies() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.0));
    let initial = [302.0, 292.0];
    let thermal = thermal(fixture.area(), initial);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let mut limits = controls();
    limits.maximum_fractional_area_change = 1e-8;
    let candidate = run(
        &fixture,
        &thermal,
        &route,
        fixture.balanced_normal(),
        limits,
    )
    .unwrap();
    let report = &candidate.report;
    let mechanical = &report.mechanical;
    assert!(mechanical.contact.path.contact.normal_force_n < 0.0);
    assert!(mechanical.dissipated_energy_j > 0.0);
    let mobility =
        MASS[0].recip() + MASS[1].recip() + LENGTH * LENGTH / 1.1e-27 + LENGTH * LENGTH / 1.7e-27;
    let resistance = KINETIC * fixture.area();
    let displacement = DT * (0.12 - 0.5 * DT * mobility * resistance);
    close(
        mechanical.dissipated_energy_j,
        resistance * displacement,
        2e-27,
    );
    close(mechanical.normal_storage_change_j, 0.0, 2e-27);
    assert_eq!(
        report.heat.dissipated_j,
        mechanical.contact.friction_dissipation_j
    );

    // Eliminate each original finite-body backward-Euler balance analytically.
    // No heat fractions or mechanical/thermal adapter calls enter this reference.
    let g = K.map(|k| k * fixture.area() / LENGTH);
    let c = MASS.map(|mass| CP * mass);
    let r = std::array::from_fn::<_, 2, _>(|i| DT * g[i] * c[i] / (c[i] + DT * g[i]));
    let face = 300.0
        + (resistance * displacement + r[0] * (initial[0] - 300.0) + r[1] * (initial[1] - 300.0))
            / (r[0] + r[1]);
    let interface = &report.heat.interfaces[0];
    close(interface.interface_temperature_k, face, 1e-7);
    assert_eq!(report.source, interface.interface_temperature_k);
    assert!(interface.body_heat_j[0] < 0.0);
    assert!(interface.body_heat_j[1] > 0.0);
    for i in 0..2 {
        let q = r[i] * (face - initial[i]);
        close(interface.conductances_w_k[i], g[i], 1e-24);
        close(interface.body_heat_j[i], q, 2e-26);
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            initial[i] + q / c[i],
            1e-9,
        );
    }
    close(
        report.combined_energy_residual_j,
        0.0,
        limits.thermal.combined_energy_tolerance_j,
    );
    close(
        report.summed_energy_residual_j,
        0.0,
        limits.thermal.combined_energy_tolerance_j,
    );
    assert!(report.heat.entropy.is_some());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(candidate.thermal.time_s(), DT);
    assert_eq!(candidate.bodies_after[0].pose(), fixture.bodies[2].pose());
    assert_eq!(candidate.bodies_after[1].pose(), fixture.bodies[0].pose());
}

#[test]
fn changing_attached_area_is_bounded_and_reversible_normal_work_is_not_heat() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.0, 0.0, 0.02));
    let thermal = thermal(fixture.area(), [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let candidate = run(
        &fixture,
        &thermal,
        &route,
        fixture.balanced_normal(),
        controls(),
    )
    .unwrap();
    let report = &candidate.report;
    let area = report.areas;
    assert!(area.endpoint_m2 > area.incoming_m2);
    assert!(area.maximum_fractional_change > 0.001);
    close(
        area.path_mean_m2,
        0.5 * (area.incoming_m2 + area.endpoint_m2),
        1e-29,
    );
    assert_eq!(area.path_extrema_m2, [area.incoming_m2, area.endpoint_m2]);
    assert_ne!(report.mechanical.normal_storage_change_j, 0.0);
    assert_ne!(report.mechanical.external_work_j, 0.0);
    assert_eq!(report.mechanical.dissipated_energy_j, 0.0);
    assert_eq!(report.heat.dissipated_j, 0.0);
    assert_eq!(report.heat.interfaces[0].body_heat_j, [0.0; 2]);
    assert_eq!(report.heat.thermal.stored_energy_change_j, 0.0);
    // Exact native mean area drives transport; endpoint area drives continuation.
    close(
        report.heat.interfaces[0].conductances_w_k[0],
        K[0] * area.path_mean_m2 / LENGTH,
        1e-24,
    );
    assert_eq!(report.route_after.contact_area_m2(), area.endpoint_m2);
    let endpoint = candidate.thermal.contact_patch(LINK).unwrap();
    assert_eq!(endpoint.link.area_m2, area.endpoint_m2);
    assert_eq!(endpoint.covered_area_m2, [area.endpoint_m2; 2]);
    FrictionHeatNetworkStep::new(
        &candidate.thermal,
        2.0 * DT,
        &[&report.route_after],
        1,
        || Ok(()),
    )
    .unwrap();
    assert!(
        FrictionHeatNetworkStep::new(&candidate.thermal, 2.0 * DT, &[&route], 1, || Ok(()))
            .is_err()
    );
    let mut strict = controls();
    strict.maximum_fractional_area_change = 0.5 * area.maximum_fractional_change;
    assert!(matches!(
        run(
            &fixture,
            &thermal,
            &route,
            fixture.balanced_normal(),
            strict
        ),
        Err(RigidDmtHeatError::AreaChange { .. })
    ));
}

#[test]
fn evolving_real_area_integrates_contact_and_ambient_heat_without_changing_other_links() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.02));
    let incoming = fixture.area();
    let initial = [302.0, 292.0];
    let ambient = [280.0, 330.0];
    let htc = [1e6, 2e6];
    let mut events = thermal_events(
        incoming,
        initial,
        [[-2.0 * LENGTH, 0.0], [0.0, 2.0 * LENGTH]],
    );
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            let i = cell.id.0 as usize;
            cell.boundary.ambient_k = ambient[i];
            cell.boundary.htc_w_m2_k = htc[i];
        }
    }
    let other = ThermalLink {
        id: LINK + 1,
        a: CellId(0),
        b: CellId(1),
        area_m2: 0.25 * incoming,
        conductance_w_k: 2e-10,
    };
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(other),
    });
    let thermal = from_events(events);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let mut limits = controls();
    limits.thermal.newton.relative_tolerance = 1e-14;
    let candidate = run(
        &fixture,
        &thermal,
        &route,
        fixture.balanced_normal(),
        limits,
    )
    .unwrap();
    let report = &candidate.report;
    let endpoint = report.mechanical.contact.path.contact.contact_area_m2;
    let mean = 0.5 * (incoming + endpoint);
    let loss = report.mechanical.contact.friction_dissipation_j;
    assert!(endpoint > incoming);
    assert!(loss > 0.0);

    // Independent finite-body Schur elimination with exactly integrated area.
    // The unchanged secondary link exchanges opposite heat; ambient exposure
    // subtracts the SAME mean contact area used by the two half conductances.
    let g = K.map(|k| k * mean / LENGTH);
    let capacity = MASS.map(|m| CP * m);
    let coupling = DT * (g[0] * g[1] / (g[0] + g[1]) + other.conductance_w_k);
    let boundary = std::array::from_fn::<_, 2, _>(|i| {
        DT * htc[i]
            * (thermal.cells()[&CellId(i as u64)].geometry.surface_area_m2 - mean - other.area_m2)
    });
    let diagonal = std::array::from_fn::<_, 2, _>(|i| capacity[i] + coupling + boundary[i]);
    let rhs = std::array::from_fn::<_, 2, _>(|i| {
        capacity[i] * (initial[i] - 300.0)
            + boundary[i] * (ambient[i] - 300.0)
            + g[i] / (g[0] + g[1]) * loss
    });
    let theta0 = (rhs[0] + coupling * rhs[1] / diagonal[1])
        / (diagonal[0] - coupling * coupling / diagonal[1]);
    let theta = [theta0, (rhs[1] + coupling * theta0) / diagonal[1]];
    let face = (g[0] * theta[0] + g[1] * theta[1] + loss / DT) / (g[0] + g[1]);
    let heat = &report.heat.interfaces[0];
    common::close(heat.interface_temperature_k, 300.0 + face, 2e-10, 1e-13);
    for i in 0..2 {
        common::close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            300.0 + theta[i],
            2e-10,
            1e-13,
        );
        common::close(heat.conductances_w_k[i], g[i], 1e-24, 1e-13);
        common::close(
            heat.body_heat_j[i],
            DT * g[i] * (face - theta[i]),
            2e-26,
            1e-12,
        );
        let cell = &thermal.cells()[&CellId(i as u64)];
        close(
            candidate.thermal.exposed_area_m2(cell.id).unwrap(),
            cell.geometry.surface_area_m2 - endpoint - other.area_m2,
            1e-28,
        );
    }
    let ambient_heat: f64 = (0..2)
        .map(|i| boundary[i] * (ambient[i] - 300.0 - theta[i]))
        .sum();
    common::close(report.heat.external_heat_j, ambient_heat, 4e-26, 1e-13);
    assert_eq!(candidate.thermal.links()[&(LINK + 1)], other);
    assert_eq!(
        candidate.thermal.contact_patch(LINK + 1),
        thermal.contact_patch(LINK + 1)
    );
    assert_eq!(thermal.contact_patch(LINK).unwrap().link.area_m2, incoming);
    assert_eq!(
        candidate.thermal.contact_patch(LINK).unwrap().link.area_m2,
        endpoint
    );
    // The returned route also passes the original material/geometry binding.
    let rebound = FrictionHeatRoute::from_link(&candidate.thermal, LINK, K).unwrap();
    assert_eq!(
        report.route_after.contact_area_m2(),
        rebound.contact_area_m2()
    );
}

#[test]
fn mean_and_endpoint_face_caps_and_pending_topology_preserve_every_owner() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.02));
    let native = fixture
        .prepared()
        .advance_with_wrenches(
            fixture.history(),
            fixture.balanced_normal(),
            &CancelGate::new(),
        )
        .unwrap();
    let incoming = fixture.area();
    let endpoint = native.report.contact.path.contact.contact_area_m2;
    let mean = 0.5 * (incoming + endpoint);
    for case in 0..3 {
        let mut events = thermal_events(
            incoming,
            [300.0; 2],
            [[-2.0 * LENGTH, 0.0], [0.0, 2.0 * LENGTH]],
        );
        if case < 2 {
            // Incoming area fits. First mean, then only endpoint, exceeds the
            // original physical face; total surface alone cannot admit it.
            let face = if case == 0 {
                0.5 * (incoming + mean)
            } else {
                0.5 * (mean + endpoint)
            };
            for event in &mut events {
                if let Event::Birth(cell) = &mut event.event {
                    cell.geometry.volume_m3 = 2.0 * LENGTH * face;
                }
            }
        } else {
            events.push(TimedEvent {
                time_s: DT,
                event: Event::Separation(LINK),
            });
        }
        let mut thermal = from_events(events);
        let original = thermal.clone();
        let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
        let mut bodies = fixture.bodies;
        let mut history = fixture.history();
        let calls = Cell::new(0usize);
        let result = apply_rigid_dmt_heat(
            &mut bodies,
            &fixture.properties,
            fixture.heat_contact(&route),
            &mut history,
            fixture.balanced_normal(),
            &mut thermal,
            0.0,
            DT,
            controls(),
            &[0.0; 2],
            |state, ids, temperatures| {
                calls.set(calls.get() + 1);
                transport(state, ids, temperatures)
            },
            &CancelGate::new(),
        );
        assert!(result.is_err());
        assert_eq!(bodies, fixture.bodies);
        assert_eq!(history, fixture.history());
        assert_eq!(thermal.cells(), original.cells());
        assert_eq!(thermal.links(), original.links());
        assert_eq!(thermal.contact_patch(LINK), original.contact_patch(LINK));
        assert_eq!(thermal.time_s(), 0.0);
        if case < 2 {
            assert_eq!(calls.get(), 0);
        } else {
            assert!(calls.get() > 0);
        }
    }
}

#[test]
fn driven_candidate_matches_original_native_pair_and_apply_publishes_original_indices() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.01));
    let mut thermal = thermal(fixture.area(), [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let mut wrenches = fixture.balanced_normal();
    wrenches[0].force_world.x = 3e-9;
    wrenches[0].torque_body.y = 2e-15;
    wrenches[1].force_world.x = -1e-9;
    let native = fixture
        .prepared()
        .advance_with_wrenches(fixture.history(), wrenches, &CancelGate::new())
        .unwrap();
    let mut bodies = fixture.bodies;
    let mut history = fixture.history();
    let report = apply_rigid_dmt_heat(
        &mut bodies,
        &fixture.properties,
        fixture.heat_contact(&route),
        &mut history,
        wrenches,
        &mut thermal,
        0.0,
        DT,
        controls(),
        &[0.0; 2],
        transport,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.mechanical, native.report);
    assert_eq!(history, native.state_after);
    assert_eq!(bodies[2], native.bodies_after[0]);
    assert_eq!(bodies[0], native.bodies_after[1]);
    assert_eq!(bodies[1], fixture.bodies[1]);
    assert_eq!(thermal.time_s(), DT);
    assert_eq!(report.heat.dissipated_j, native.report.dissipated_energy_j);
    assert_ne!(report.mechanical.external_work_j, 0.0);
}

#[test]
fn source_thermal_and_cancellation_refusals_preserve_bodies_history_and_all_cells() {
    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.0));
    let initial = thermal(fixture.area(), [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&initial, LINK, K).unwrap();
    for failure in 0..3 {
        let mut bodies = fixture.bodies;
        let mut history = fixture.history();
        let mut thermal = initial.clone();
        let gate = CancelGate::new();
        let calls = Cell::new(0usize);
        let result = apply_rigid_dmt_heat(
            &mut bodies,
            &fixture.properties,
            fixture.heat_contact(&route),
            &mut history,
            fixture.balanced_normal(),
            &mut thermal,
            0.0,
            DT,
            controls(),
            &[0.0; 2],
            |state, ids, temperatures| {
                calls.set(calls.get() + 1);
                if failure == 0 {
                    return Err("source fixture refusal");
                }
                let mut source = transport(state, ids, temperatures)?;
                if failure == 1 {
                    // Both original receiver conductivities vanish: physical refusal.
                    source.interfaces[0].conductivity_w_m_k = [0.0; 2];
                } else {
                    gate.request();
                }
                Ok(source)
            },
            &gate,
        );
        if failure == 0 {
            assert!(matches!(
                result,
                Err(RigidDmtHeatError::Heat(CohesiveHeatError::Source(
                    "source fixture refusal"
                )))
            ));
        } else if failure == 1 {
            assert!(matches!(
                result,
                Err(RigidDmtHeatError::Heat(CohesiveHeatError::Thermal(_)))
            ));
        } else {
            assert!(matches!(
                result,
                Err(RigidDmtHeatError::Heat(CohesiveHeatError::Thermal(
                    DepositionError::Cancelled
                )))
            ));
        }
        assert!(calls.get() > 0);
        assert_eq!(bodies, fixture.bodies);
        assert_eq!(history, fixture.history());
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), initial.cells());
        assert_eq!(thermal.contact_patch(LINK), initial.contact_patch(LINK));
    }
}

#[test]
fn gap_touch_opening_area_and_source_receiver_mismatches_are_atomic() {
    for approach in [-0.2 * common::Z, 0.0] {
        let fixture = Fixture::new(approach, Vec3::ZERO);
        // A supplied finite face cannot turn attraction or geometric touch into contact.
        let thermal = thermal(1e-16, [300.0; 2]);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
        assert!(matches!(
            run(&fixture, &thermal, &route, [Wrench::ZERO; 2], controls()),
            Err(RigidDmtHeatError::Invalid(_))
        ));
    }
    let opening = Fixture::new(0.02 * common::Z, Vec3::new(0.0, 0.0, -0.3));
    let open_thermal = thermal(opening.area(), [300.0; 2]);
    let open_route = FrictionHeatRoute::from_link(&open_thermal, LINK, K).unwrap();
    assert!(matches!(
        run(
            &opening,
            &open_thermal,
            &open_route,
            opening.balanced_normal(),
            controls()
        ),
        Err(RigidDmtHeatError::Invalid(_))
    ));

    let fixture = Fixture::new(0.2 * common::Z, Vec3::new(0.12, 0.0, 0.0));
    let thermal = thermal(fixture.area(), [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    for mismatch in 0..3 {
        let mut bodies = fixture.bodies;
        let mut history = fixture.history();
        let mut thermal = thermal.clone();
        let cells = thermal.cells().clone();
        let mut contact = fixture.heat_contact(&route);
        match mismatch {
            0 => contact.bodies[0].body_id = 123,
            1 => contact.bodies.swap(0, 1),
            _ => contact.bodies[0].thermal_material = MATERIAL_IDS[1],
        }
        let calls = Cell::new(0usize);
        let result = apply_rigid_dmt_heat(
            &mut bodies,
            &fixture.properties,
            contact,
            &mut history,
            fixture.balanced_normal(),
            &mut thermal,
            0.0,
            DT,
            controls(),
            &[0.0; 2],
            |state, ids, temperatures| {
                calls.set(calls.get() + 1);
                transport(state, ids, temperatures)
            },
            &CancelGate::new(),
        );
        assert!(matches!(result, Err(RigidDmtHeatError::Invalid(_))));
        assert_eq!(calls.get(), 0);
        assert_eq!(bodies, fixture.bodies);
        assert_eq!(history, fixture.history());
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
    // This is a valid original thermal route: opposing faces, the same centre,
    // unchanged A/B IDs/materials and area. Its cell interiors occupy the wrong
    // physical sides of the native A-to-B normal and must not exchange roles.
    let mut reversed = thermal_with_ends(
        fixture.area(),
        [302.0, 292.0],
        [[0.0, 2.0 * LENGTH], [-2.0 * LENGTH, 0.0]],
    );
    let reversed_route = FrictionHeatRoute::from_link(&reversed, LINK, K).unwrap();
    let original_cells = reversed.cells().clone();
    let mut bodies = fixture.bodies;
    let mut history = fixture.history();
    let calls = Cell::new(0usize);
    let result = apply_rigid_dmt_heat(
        &mut bodies,
        &fixture.properties,
        fixture.heat_contact(&reversed_route),
        &mut history,
        fixture.balanced_normal(),
        &mut reversed,
        0.0,
        DT,
        controls(),
        &[0.0; 2],
        |state, ids, temperatures| {
            calls.set(calls.get() + 1);
            transport(state, ids, temperatures)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        result,
        Err(RigidDmtHeatError::Invalid(
            "thermal A/B interiors must retain the native A-to-B normal ordering"
        ))
    ));
    assert_eq!(calls.get(), 0);
    assert_eq!(bodies, fixture.bodies);
    assert_eq!(history, fixture.history());
    assert_eq!(reversed.cells(), &original_cells);
    assert_eq!(reversed.time_s(), 0.0);
    let wrong_thermal = thermal_with_wrong_area(&fixture);
    let wrong_route = FrictionHeatRoute::from_link(&wrong_thermal, LINK, K).unwrap();
    assert!(matches!(
        run(
            &fixture,
            &wrong_thermal,
            &wrong_route,
            fixture.balanced_normal(),
            controls()
        ),
        Err(RigidDmtHeatError::Invalid(_))
    ));
}

fn thermal_with_wrong_area(fixture: &Fixture) -> DepositionSimulation<'static> {
    thermal(1.01 * fixture.area(), [300.0; 2])
}
