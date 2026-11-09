//! Analytic damping, original nonlinear calorics, shared receivers and clocks.

use fs_couple::rigid_damping_heat::{
    RigidDampingThermalConfig, RigidDampingThermalError, RigidDampingThermalInput,
    advance_candidate,
};
use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
        network::{ViscousNetworkLimits, ViscousNetworkLink},
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads, viscous::network::ViscousNetworkMotionControls,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, HeatTransfer, MaterialId, Reservoir,
    StepConfig, ThermalLaw, ThermalLink, ThermalStorageBasis, TimedEvent,
};
use std::collections::BTreeMap;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, allowance {tolerance:e}"
    );
}

// Synthetic calorics: e(T) = a*(T-300) + b*(T-300)^2/2.
// H and fixed-volume U keep distinct declared storage bases.
struct Caloric {
    a: f64,
    b: f64,
    maximum: f64,
    basis: ThermalStorageBasis,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
        if !e.is_finite() || e < 0.0 || e > self.maximum {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + 2.0 * e / (self.a + (self.a * self.a + 2.0 * self.b * e).sqrt()))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.a
    }
}
fn laws() -> [Caloric; 2] {
    [
        Caloric {
            a: 0.5,
            b: 0.2,
            maximum: 1e6,
            basis: ThermalStorageBasis::Enthalpy,
        },
        Caloric {
            a: 0.8,
            b: 0.3,
            maximum: 1e6,
            basis: ThermalStorageBasis::InternalEnergy,
        },
    ]
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w,
        reservoir: Some(Reservoir {
            area_m2: 0.1,
            conductance_w_k: 0.0,
            temperature_k: 300.0,
        }),
    }
}
fn birth(id: u64, time_s: f64, source_w: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id + 10,
            body_id: id + 20,
            material_id: MaterialId(id),
            mass_kg: 2.0 + id as f64,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0,
                surface_area_m2: 6.0,
            },
            boundary: boundary(source_w),
        }),
    }
}
fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 16,
        max_links: 16,
        max_events: 16,
        max_substeps: 128,
    }
}
fn thermal<'m>(laws: &'m [Caloric; 2], events: bool, activate: bool) -> DepositionSimulation<'m> {
    let mut schedule = vec![
        birth(0, 0.0, if events { 2.0 } else { 0.0 }),
        birth(1, 0.0, 0.0),
    ];
    if events {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 17,
                a: CellId(0),
                b: CellId(1),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        schedule.push(TimedEvent {
            time_s: 0.35,
            event: Event::Boundary(CellId(0), boundary(4.0)),
        });
    }
    let initial = DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        budget(),
    )
    .unwrap();
    if activate {
        initial
            .advance_to(0.0, config().thermal, || Ok(()))
            .unwrap()
            .0
    } else {
        initial
    }
}
fn config() -> RigidDampingThermalConfig {
    RigidDampingThermalConfig {
        mechanical: ViscousNetworkMotionControls {
            source: ViscousConnectorControls {
                impulse: InternalWrenchImpulseControls {
                    angular_origin_world_m: Vec3::ZERO,
                    linear_tolerance_ns: 1e-10,
                    angular_tolerance_nms: 1e-10,
                    work_tolerance_j: 1e-10,
                },
                linear_speed_tolerance_m_s: 1e-10,
                angular_speed_tolerance_rad_s: 1e-10,
            },
            source_substeps: 1,
            limits: ViscousNetworkLimits {
                maximum_bodies: 16,
                maximum_connectors: 16,
                maximum_pair_relaxations: 10_000,
                maximum_live_bytes: 1_000_000,
            },
        },
        thermal: StepConfig {
            max_dt_s: 0.125,
            energy_tolerance_j: 1e-10,
        },
        maximum_heat_bindings: 8,
        energy_tolerance_j: 1e-9,
    }
}
fn axial(bodies: [usize; 2], c: f64) -> ViscousNetworkLink {
    let mut coefficients = [[0.0; 6]; 6];
    coefficients[0][0] = c;
    ViscousNetworkLink {
        bodies,
        arms_body_m: [Vec3::ZERO; 2],
        damping: PassiveBushingMatrix::try_new(
            coefficients,
            BushingCoordinateScales {
                translation_m: 1.0,
                rotation_rad: 1.0,
            },
        )
        .unwrap(),
    }
}
struct Bodies {
    states: [RigidBodyState; 4],
    properties: [MassProperties; 4],
    links: [ViscousNetworkLink; 2],
    loads: [Wrench; 4],
    gravity: Gravity,
}
impl Bodies {
    fn new(c: [f64; 2]) -> Self {
        let properties =
            [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(); 4];
        let states = core::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(Vec3::new(i as f64, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(2.0 * [1.0, -1.0, 0.5, -0.5][i], 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap()
        });
        Self {
            states,
            properties,
            links: [axial([0, 1], c[0]), axial([2, 3], c[1])],
            loads: [Wrench::ZERO; 4],
            gravity: Gravity::ZERO,
        }
    }
    fn input<'a>(&'a self, routes: &'a [BTreeMap<CellId, f64>]) -> RigidDampingThermalInput<'a> {
        RigidDampingThermalInput {
            time_s: 0.0,
            states: &self.states,
            properties: &self.properties,
            links: &self.links,
            loads: HeldExternalLoads {
                wrenches: &self.loads,
                gravity: self.gravity,
            },
            heat_fractions: routes,
        }
    }
}
fn routes() -> [BTreeMap<CellId, f64>; 2] {
    [
        BTreeMap::from([(CellId(0), 0.25), (CellId(1), 0.5)]),
        BTreeMap::from([(CellId(0), 0.6), (CellId(1), 0.2)]),
    ]
}
fn exact_heat(t: f64) -> [f64; 2] {
    // Two 2-kg masses: relative mobility 1/kg. Initial relative kinetic
    // energies are 2 J and 0.5 J, each decaying at twice its velocity rate.
    [2.0 * -(-1.4 * t).exp_m1(), 0.5 * -(-2.2 * t).exp_m1()]
}

#[derive(Clone, Debug, Default, PartialEq)]
struct Observe {
    supplied: bool,
    interval_watts: Option<f64>,
    fail_at: Option<usize>,
    hooks: [usize; 9],
    intervals: Vec<(f64, f64, f64)>,
    events: Vec<f64>,
}
impl DepositionObserver for Observe {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.hooks[0] += 1;
        Ok(())
    }
    fn requires_step_interval(&self) -> bool {
        self.supplied || self.interval_watts.is_some()
    }
    fn next_transition_time_s(
        &mut self,
        t: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.hooks[1] += 1;
        Ok((self.supplied && t.time_s() < 0.7).then_some(0.7))
    }
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.hooks[2] += 1;
        Ok(if self.supplied {
            BTreeMap::from([(CellId(0), 1.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.hooks[3] += 1;
        Ok(if self.supplied {
            BTreeMap::from([(17, 0.2)])
        } else {
            BTreeMap::new()
        })
    }
    fn ambient_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.hooks[4] += 1;
        Ok(if self.supplied {
            BTreeMap::from([(CellId(0), 0.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.hooks[5] += 1;
        Ok(if self.supplied {
            BTreeMap::from([(CellId(1), 0.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn heat_transfers(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.hooks[6] += 1;
        Ok(if self.supplied {
            vec![HeatTransfer {
                a: CellId(0),
                b: CellId(1),
                power_w: 0.1,
            }]
        } else {
            Vec::new()
        })
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.hooks[7] += 1;
        Ok(if self.supplied {
            proposed.min(0.5)
        } else {
            proposed
        })
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        _: f64,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.hooks[8] += 1;
        Ok(if let Some(watts) = self.interval_watts {
            BTreeMap::from([(CellId(0), watts)])
        } else if self.supplied {
            BTreeMap::from([(CellId(1), 0.25)])
        } else {
            BTreeMap::new()
        })
    }
    fn event(&mut self, _: &Event, t: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events.push(t.time_s());
        Ok(())
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        t: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.intervals.push((start, dt, t.time_s()));
        if self.fail_at == Some(self.intervals.len()) {
            return Err(DepositionError::InvalidInput("synthetic observer refusal"));
        }
        Ok(())
    }
}

#[test]
fn exact_connector_losses_heat_shared_nonlinear_h_and_u_receivers_once() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, false, true);
    let routes = routes();
    let out = advance_candidate(bodies.input(&routes), &initial, &(), 0.5, config(), || {
        false
    })
    .unwrap();
    let heat = exact_heat(0.5);
    for i in 0..2 {
        near(out.report().connector_heat_j[i], heat[i], 5e-14);
    }
    let deposited = [
        0.25 * heat[0] + 0.6 * heat[1],
        0.5 * heat[0] + 0.2 * heat[1],
    ];
    for i in 0..2 {
        let id = CellId(i as u64);
        let h = deposited[i] / initial.cells()[&id].mass_kg;
        near(out.thermal().cells()[&id].specific_h_j_kg, h, 3e-14);
        near(
            out.thermal().temperature_k(id).unwrap(),
            laws[i].temperature_k(h).unwrap(),
            6e-14,
        );
        let mut expected = initial.cells()[&id].clone();
        expected.specific_h_j_kg = out.thermal().cells()[&id].specific_h_j_kg;
        assert_eq!(out.thermal().cells()[&id], expected);
    }
    near(
        out.report().unallocated_energy_j,
        0.25 * heat[0] + 0.2 * heat[1],
        3e-14,
    );
    near(out.report().deposited_heat_j, deposited.iter().sum(), 5e-14);
    near(
        out.report().thermal.external_heat_j,
        out.report().deposited_heat_j,
        3e-14,
    );
    near(out.report().external_thermal_heat_j, 0.0, 3e-14);
    near(out.report().combined_energy_residual_j, 0.0, 3e-13);
    assert_eq!(out.report().thermal.substeps, 4);
    assert_eq!(out.report().pair_relaxations, 24);
    assert_eq!(out.thermal().time_s(), 0.5);
    for i in 0..4 {
        let rate = if i < 2 { 0.7 } else { 1.1 };
        near(
            out.states()[i].linear_momentum_world().x,
            bodies.states[i].linear_momentum_world().x * (-rate * 0.5_f64).exp(),
            4e-14,
        );
    }
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn event_clipped_intervals_forward_transport_and_external_sources_on_one_clock() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, true, true);
    let routes = routes();
    let observer = Observe {
        supplied: true,
        ..Default::default()
    };
    let mut controls = config();
    controls.thermal.max_dt_s = 0.6;
    let out = advance_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        1.0,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(out.observer().events, [0.35]);
    assert_eq!(
        out.observer().intervals,
        [(0.0, 0.35, 0.35), (0.35, 0.35, 0.7), (0.7, 1.0 - 0.7, 1.0)]
    );
    assert_eq!(out.observer().hooks, [1, 3, 3, 3, 3, 3, 3, 3, 3]);
    for (actual, exact) in out.report().connector_heat_j.iter().zip(exact_heat(1.0)) {
        near(*actual, exact, 5e-14);
    }
    near(
        out.report().external_thermal_heat_j,
        2.0 * 0.35 + 4.0 * 0.65 + 1.0 + 0.25,
        2e-13,
    );
    near(out.report().combined_energy_residual_j, 0.0, 5e-13);
    assert_eq!(out.thermal().cells()[&CellId(0)].boundary.source_w, 4.0);
    assert_eq!(out.thermal().links()[&17].conductance_w_k, 0.0);
    assert_eq!(observer.intervals.len(), 0);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn zero_damping_reproduces_native_loaded_motion_and_unmodified_thermal_owner() {
    let mut bodies = Bodies::new([0.0, 0.0]);
    let laws = laws();
    let initial = thermal(&laws, false, true);
    let routes = routes();
    bodies.loads = [Wrench {
        force_world: Vec3::new(0.2, -0.1, 0.3),
        torque_body: Vec3::new(0.01, 0.02, -0.03),
    }; 4];
    bodies.gravity = Gravity::new(Vec3::new(0.0, -9.81, 0.0)).unwrap();
    let observer = Observe::default();
    let (heat, observed, audit) = initial
        .advance_observed(&observer, 0.5, config().thermal, || Ok(()))
        .unwrap();
    let out = advance_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.5,
        config(),
        || false,
    )
    .unwrap();
    let mut expected = bodies.states;
    for (_, dt, _) in &observed.intervals {
        for i in 0..4 {
            expected[i] = RigidBodyIntegrator::new(bodies.gravity)
                .step_cayley(expected[i], bodies.properties[i], bodies.loads[i], *dt)
                .unwrap()
                .state_after;
        }
    }
    assert_eq!(out.states(), expected);
    assert_eq!(out.thermal().cells(), heat.cells());
    assert_eq!(out.observer(), &observed);
    assert_eq!(out.report().thermal, audit);
    assert_eq!(out.report().connector_heat_j, [0.0, 0.0]);
    assert_eq!(out.report().dissipated_energy_j, 0.0);
    assert_eq!(out.report().deposited_heat_j, 0.0);
    assert_eq!(out.report().pair_relaxations, 0);
    near(out.report().combined_energy_residual_j, 0.0, 2e-12);
}

#[test]
fn same_time_births_finish_before_receiver_admission() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, false, false);
    let routes = routes();
    let observer = Observe::default();
    let out = advance_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.25,
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(out.observer().events, [0.0, 0.0]);
    assert_eq!(out.report().thermal.events, 2);
    assert_eq!(out.report().thermal.substeps, 2);
    assert_eq!(out.thermal().cells().len(), 2);
    near(out.report().connector_heat_j[0], exact_heat(0.25)[0], 4e-14);
    assert!(initial.cells().is_empty());
}

#[test]
fn empty_receiver_network_keeps_motion_clock_and_unallocated_loss() {
    let bodies = Bodies::new([0.7, 1.1]);
    let initial = DepositionSimulation::new(vec![], &[], budget()).unwrap();
    let routes = [BTreeMap::new(), BTreeMap::new()];
    let out = advance_candidate(bodies.input(&routes), &initial, &(), 0.5, config(), || {
        false
    })
    .unwrap();
    assert!(out.thermal().cells().is_empty());
    assert_eq!(out.report().thermal.substeps, 4);
    near(
        out.report().unallocated_energy_j,
        exact_heat(0.5).iter().sum(),
        6e-14,
    );
    assert_eq!(out.report().deposited_heat_j, 0.0);
    near(out.report().combined_energy_residual_j, 0.0, 3e-13);
    assert_eq!(out.thermal().time_s(), 0.5);
}

#[test]
fn partition_clock_budget_and_not_yet_born_receivers_refuse() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, false, true);
    let mut invalid_routes = routes();
    invalid_routes[0] = BTreeMap::from([(CellId(0), 1.0), (CellId(1), 2.0_f64.powi(-60))]);
    assert!(matches!(
        advance_candidate(
            bodies.input(&invalid_routes),
            &initial,
            &(),
            0.5,
            config(),
            || false
        ),
        Err(RigidDampingThermalError::Invalid(_))
    ));
    invalid_routes[0] = BTreeMap::from([(CellId(7), 0.5)]);
    assert!(
        advance_candidate(
            bodies.input(&invalid_routes),
            &initial,
            &(),
            0.5,
            config(),
            || false
        )
        .is_err()
    );
    let routes = routes();
    let mut input = bodies.input(&routes);
    input.time_s = 0.1;
    assert!(advance_candidate(input, &initial, &(), 0.5, config(), || false).is_err());
    let mut controls = config();
    controls.maximum_heat_bindings = 3;
    assert!(
        advance_candidate(bodies.input(&routes), &initial, &(), 0.5, controls, || {
            false
        })
        .is_err()
    );
    let future = DepositionSimulation::new(
        vec![birth(0, 0.0, 0.0), birth(1, 0.2, 0.0)],
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        budget(),
    )
    .unwrap();
    assert!(
        advance_candidate(bodies.input(&routes), &future, &(), 0.5, config(), || false).is_err()
    );
    assert!(future.cells().is_empty());
    assert_eq!(future.time_s(), 0.0);
}

#[test]
fn late_observer_caloric_and_final_cancellation_refusals_are_transactional_and_retryable() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, false, true);
    let routes = routes();
    let observer = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    let cells = initial.cells().clone();
    assert!(matches!(
        advance_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            config(),
            || false
        ),
        Err(RigidDampingThermalError::Observer(_))
    ));
    assert!(observer.intervals.is_empty());
    assert_eq!(initial.cells(), &cells);
    assert_eq!(initial.time_s(), 0.0);
    let clean = Observe::default();
    let mut polls = 0;
    let expected = advance_candidate(
        bodies.input(&routes),
        &initial,
        &clean,
        0.5,
        config(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut again = 0;
    assert!(matches!(
        advance_candidate(
            bodies.input(&routes),
            &initial,
            &clean,
            0.5,
            config(),
            || {
                again += 1;
                again == polls
            }
        ),
        Err(RigidDampingThermalError::Cancelled)
    ));
    let retry = advance_candidate(
        bodies.input(&routes),
        &initial,
        &clean,
        0.5,
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(retry.observer(), expected.observer());
    let mut narrow = self::laws();
    narrow[0].maximum = 1e-6;
    let cold = thermal(&narrow, false, true);
    let saved = cold.cells().clone();
    assert!(matches!(
        advance_candidate(bodies.input(&routes), &cold, &clean, 0.5, config(), || {
            false
        }),
        Err(RigidDampingThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(cold.cells(), &saved);
    assert_eq!(cold.time_s(), 0.0);
}

#[test]
fn a_large_outside_interval_source_cannot_erase_physical_damper_heat() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = thermal(&laws, false, true);
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(0), 1.0)]),
    ];
    let observer = Observe {
        interval_watts: Some(1e20),
        ..Default::default()
    };
    // The original cell would also refuse this large input on its caloric
    // domain. The independent handoff gate must catch erased damper joules
    // BEFORE thermal integration, while all caller histories remain unchanged.
    assert!(matches!(
        advance_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            config(),
            || false
        ),
        Err(RigidDampingThermalError::Energy { .. })
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty());
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
}
