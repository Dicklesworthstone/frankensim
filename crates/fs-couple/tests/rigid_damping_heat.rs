//! Analytic damping, original nonlinear calorics, shared receivers and clocks.

#[path = "rigid_damping_heat/shift.rs"]
mod calibrated_shift;

#[path = "rigid_damping_heat/driven.rs"]
mod driven;

#[path = "rigid_damping_heat/driven_network.rs"]
mod driven_network;

#[path = "rigid_damping_heat/relaxation.rs"]
mod relaxation;

use fs_couple::rigid_damping_heat::temperature::{
    TemperatureDampingLaw, advance_candidate as advance_temperature,
};
use fs_couple::rigid_damping_heat::{
    ImplicitThermalScheme, RigidDampingImplicitConfig, RigidDampingThermalConfig,
    RigidDampingThermalError, RigidDampingThermalInput, advance_candidate,
    advance_implicit_candidate,
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
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
    ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionState,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
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
    fn specific_entropy_change_j_kg_k(
        &self,
        initial: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial_t = self.temperature_k(initial)?;
        self.temperature_k(final_h)?;
        let delta_t = 2.0 * (final_h - initial)
            / ((self.a * self.a + 2.0 * self.b * final_h).sqrt()
                + (self.a * self.a + 2.0 * self.b * initial).sqrt());
        Ok(Some(
            self.b * delta_t + (self.a - 300.0 * self.b) * (delta_t / initial_t).ln_1p(),
        ))
    }
}
impl ConductivityLaw for Caloric {
    fn conductivity_w_m_k(&self, e: f64) -> Result<f64, DepositionError> {
        self.temperature_k(e)?;
        Ok(1000.0)
    }
}
impl ImplicitConductionLaw for Caloric {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [0.0, self.maximum]
    }
    fn state_at_specific_enthalpy(
        &self,
        e: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let temperature_k = self.temperature_k(e)?;
        Ok(ImplicitConductionState {
            temperature_k,
            temperature_slope_kg_k_j: 1.0 / (self.a + self.b * (temperature_k - 300.0)),
            conductivity_w_m_k: self.conductivity_w_m_k(e)?,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
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
fn opposing_large_sources_preserve_physical_damper_heat() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let large_watts = 2.0_f64.powi(80);
    let initial = DepositionSimulation::new(
        vec![birth(0, 0.0, -large_watts), birth(1, 0.0, 0.0)],
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        budget(),
    )
    .unwrap()
    .advance_to(0.0, config().thermal, || Ok(()))
    .unwrap()
    .0;
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(0), 1.0)]),
    ];
    let observer = Observe {
        interval_watts: Some(large_watts),
        ..Default::default()
    };
    // Authored cooling and outside heating cancel, leaving independently
    // integrated physical damping. The interval wrapper must preserve that
    // small contribution before the native owner reduces all source terms.
    let candidate = advance_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.5,
        config(),
        || false,
    )
    .unwrap();
    let heat: f64 = exact_heat(0.5).into_iter().sum();
    near(candidate.report().deposited_heat_j, heat, 5e-14);
    near(candidate.report().external_thermal_heat_j, 0.0, 5e-14);
    near(
        candidate.thermal().cells()[&CellId(0)].specific_h_j_kg,
        heat / 2.0,
        5e-14,
    );
    near(candidate.report().combined_energy_residual_j, 0.0, 3e-13);
    assert_eq!(
        candidate.thermal().cells()[&CellId(0)].boundary.source_w,
        -large_watts
    );
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty());
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
}

// Synthetic memoryless source, with an independent closed-form continuum
// solution below. This is not a fitted or measured material assertion.
struct LinearTemperatureDamping {
    beta_per_k: f64,
    maximum_k: f64,
}
impl TemperatureDampingLaw for LinearTemperatureDamping {
    type Error = &'static str;
    fn at_temperature(
        &self,
        _connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, Self::Error> {
        if !(300.0..=self.maximum_k).contains(&temperature_k) {
            return Err("synthetic source temperature domain");
        }
        let c = reference.native_coefficients()[0][0]
            * (1.0 + self.beta_per_k * (temperature_k - 300.0));
        let mut matrix = [[0.0; 6]; 6];
        matrix[0][0] = c;
        PassiveBushingMatrix::try_new(
            matrix,
            BushingCoordinateScales {
                translation_m: 1.0,
                rotation_rad: 1.0,
            },
        )
        .map_err(|_| "synthetic source passive form")
    }
}
fn linear_calorics() -> [Caloric; 2] {
    core::array::from_fn(|i| Caloric {
        a: 0.5,
        b: 0.0,
        maximum: 1e6,
        basis: if i == 0 {
            ThermalStorageBasis::Enthalpy
        } else {
            ThermalStorageBasis::InternalEnergy
        },
    })
}

fn stiff_thermal<'m>(laws: &'m [Caloric; 2], maximum_steps: usize) -> DepositionSimulation<'m> {
    let mut limits = budget();
    limits.max_substeps = maximum_steps;
    let mut schedule = vec![
        birth(0, 0.0, 0.0),
        birth(1, 0.0, 0.0),
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 17,
                a: CellId(0),
                b: CellId(1),
                area_m2: 1.0,
                conductance_w_k: 1000.0,
            }),
        },
    ];
    for event in &mut schedule {
        if let Event::Birth(cell) = &mut event.event {
            // The independent reference is an isolated two-cell system.
            // Unlike the general fixture, it has no external reservoir owner.
            cell.boundary.reservoir = None;
        }
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        limits,
    )
    .unwrap()
    .advance_to(0.0, config().thermal, || Ok(()))
    .unwrap()
    .0
}

fn stiff_graph<'m>(
    laws: &'m [Caloric; 2],
    initial: &DepositionSimulation<'_>,
) -> ImplicitConductionGraph<'m> {
    let conduction = ConductionObserver::new(
        initial,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        &[ConductiveLink {
            link_id: 17,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    ImplicitConductionGraph::new(
        initial,
        &conduction,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
}

fn implicit_config(scheme: ImplicitThermalScheme) -> RigidDampingImplicitConfig {
    let original = config();
    RigidDampingImplicitConfig {
        mechanical: original.mechanical,
        thermal: ImplicitConductionControls {
            maximum_dt_s: 0.25,
            cell_tolerance_j: 1e-10,
            total_tolerance_j: 1e-9,
            ..Default::default()
        },
        scheme,
        maximum_heat_bindings: original.maximum_heat_bindings,
        energy_tolerance_j: 1e-8,
    }
}

#[test]
fn stiff_conduction_and_rigid_damping_share_be_and_sdirk_intervals() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = linear_calorics();
    let initial = stiff_thermal(&laws, 8);
    let graph = stiff_graph(&laws, &initial);
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(0), 1.0)]),
    ];
    let observer = Observe::default();
    // G=1000 W/K gives a sub-millisecond explicit stability restriction;
    // the same physical system fits in two admitted implicit intervals.
    assert!(matches!(
        advance_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            config(),
            || false
        ),
        Err(RigidDampingThermalError::Thermal(DepositionError::Budget(
            _
        )))
    ));
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (out, native) = advance_implicit_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            &graph,
            implicit_config(scheme),
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(native.energy, out.report().thermal);
        assert_eq!(native.energy.substeps, 2);
        assert_eq!(native.rejected_intervals, 0);
        assert_eq!(
            out.observer().intervals,
            [(0.0, 0.25, 0.25), (0.25, 0.25, 0.5)]
        );
        assert_eq!(out.observer().hooks[0], 1);
        assert_eq!(out.observer().hooks[8], 2);
        assert_eq!(out.report().pair_relaxations, 12);
        assert!(native.entropy.has_unclassified_sources);
        assert!(native.entropy.thermal_storage_change_j_k > 0.0);
        assert!(native.entropy.conduction_generation_j_k >= 0.0);

        // Independent two-capacity discrete solution. C0=1 and C1=1.5 J/K;
        // d=T0-T1 obeys d'=-G(1/C0+1/C1)d + mean_power/C0. Each
        // mechanical interval supplies its exact native loss only once.
        let mut total_energy = 0.0;
        let mut difference = 0.0;
        let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
        for &(start, dt, end) in &out.observer().intervals {
            let heat: f64 = exact_heat(end)
                .into_iter()
                .zip(exact_heat(start))
                .map(|(after, before)| after - before)
                .sum();
            total_energy += heat;
            let stiffness = dt * 1000.0 * (1.0 + 1.0 / 1.5);
            difference = match scheme {
                ImplicitThermalScheme::BackwardEuler => (difference + heat) / (1.0 + stiffness),
                ImplicitThermalScheme::Sdirk2 => {
                    let first = (difference + gamma * heat) / (1.0 + gamma * stiffness);
                    (difference + heat - (1.0 - gamma) * stiffness * first)
                        / (1.0 + gamma * stiffness)
                }
            };
        }
        let cell_zero_energy = total_energy / 2.5 + (1.5 / 2.5) * difference;
        near(
            out.thermal().cells()[&CellId(0)].specific_h_j_kg,
            cell_zero_energy / 2.0,
            2e-9,
        );
        near(
            out.thermal().cells()[&CellId(1)].specific_h_j_kg,
            (total_energy - cell_zero_energy) / 3.0,
            2e-9,
        );
        near(out.report().deposited_heat_j, total_energy, 1e-13);
        near(out.report().combined_energy_residual_j, 0.0, 1e-9);
        for (i, expected_heat) in exact_heat(0.5).into_iter().enumerate() {
            near(out.report().connector_heat_j[i], expected_heat, 1e-13);
        }
        for i in 0..4 {
            let rate = if i < 2 { 0.7 } else { 1.1 };
            near(
                out.states()[i].linear_momentum_world().x,
                bodies.states[i].linear_momentum_world().x * (-rate * 0.5_f64).exp(),
                1e-13,
            );
        }
    }
    assert_eq!(initial.time_s(), 0.0);
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    assert_eq!(observer, Observe::default());
}

#[test]
fn implicit_retries_replace_mechanical_receipts_and_retain_attempted_work() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = laws();
    let initial = stiff_thermal(&laws, 128);
    let graph = stiff_graph(&laws, &initial);
    let routes = routes();
    let observer = Observe::default();
    let mut controls = implicit_config(ImplicitThermalScheme::Sdirk2);
    controls.thermal.maximum_newton_iterations = 1;
    controls.thermal.cell_tolerance_j = 1e-5;
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (out, native) = advance_implicit_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.25,
        &graph,
        controls,
        &mut workspace,
        || false,
    )
    .unwrap();
    assert!(native.rejected_intervals > 0);
    let attempts = native.energy.substeps + native.rejected_intervals;
    assert_eq!(out.observer().hooks[0], 1);
    assert_eq!(out.observer().hooks[8], attempts);
    assert_eq!(out.observer().intervals.len(), native.energy.substeps);
    assert_eq!(out.report().pair_relaxations, 6 * attempts);
    assert!(out.report().pair_relaxations > 6 * native.energy.substeps);
    let heat = exact_heat(0.25);
    for (i, expected) in heat.into_iter().enumerate() {
        near(out.report().connector_heat_j[i], expected, 2e-12);
    }
    near(
        out.report().deposited_heat_j,
        0.75 * heat[0] + 0.8 * heat[1],
        2e-12,
    );
    near(out.report().combined_energy_residual_j, 0.0, 1e-8);
    for i in 0..4 {
        let rate = if i < 2 { 0.7 } else { 1.1 };
        near(
            out.states()[i].linear_momentum_world().x,
            bodies.states[i].linear_momentum_world().x * (-rate * 0.25_f64).exp(),
            2e-12,
        );
    }
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(observer, Observe::default());
}

#[test]
fn implicit_rigid_refusal_and_late_cancellation_publish_no_partial_state() {
    let bodies = Bodies::new([0.7, 1.1]);
    let laws = linear_calorics();
    let initial = stiff_thermal(&laws, 8);
    let saved_cells = initial.cells().clone();
    let graph = stiff_graph(&laws, &initial);
    let routes = routes();
    let controls = implicit_config(ImplicitThermalScheme::Sdirk2);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let refusing = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        advance_implicit_candidate(
            bodies.input(&routes),
            &initial,
            &refusing,
            0.5,
            &graph,
            controls,
            &mut workspace,
            || false
        ),
        Err(RigidDampingThermalError::Observer(
            DepositionError::InvalidInput("synthetic observer refusal")
        ))
    ));
    assert!(refusing.intervals.is_empty());
    let observer = Observe::default();
    let mut polls = 0;
    let (expected, expected_report) = advance_implicit_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.5,
        &graph,
        controls,
        &mut workspace,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut again = 0;
    assert!(matches!(
        advance_implicit_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            &graph,
            controls,
            &mut workspace,
            || {
                again += 1;
                again == polls
            }
        ),
        Err(RigidDampingThermalError::Cancelled)
    ));
    let mut insufficient = controls;
    insufficient.thermal.maximum_steps = 1;
    assert!(matches!(
        advance_implicit_candidate(
            bodies.input(&routes),
            &initial,
            &observer,
            0.5,
            &graph,
            insufficient,
            &mut workspace,
            || false
        ),
        Err(RigidDampingThermalError::ImplicitThermal(
            ImplicitConductionError::Thermal(DepositionError::Budget(_))
        ))
    ));
    let (retry, retry_report) = advance_implicit_candidate(
        bodies.input(&routes),
        &initial,
        &observer,
        0.5,
        &graph,
        controls,
        &mut workspace,
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(retry.observer(), expected.observer());
    assert_eq!(retry_report, expected_report);
    assert_eq!(initial.cells(), &saved_cells);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(observer, Observe::default());
}

#[test]
fn two_way_feedback_matches_independent_nonlinear_decay_with_first_order_refinement() {
    let bodies = Bodies::new([0.7, 0.0]);
    let laws = linear_calorics();
    let initial = thermal(&laws, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let mut input = bodies.input(&routes);
    input.links = &bodies.links[..1];
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    // m1=m2=2 kg gives relative mobility 1/kg and E0=2 J. Cell 0 has
    // capacity 2*0.5=1 J/K, so T=300+(E0-E). Then E'=-a E+b E².
    // This continuum reference is independent of the split source algorithm.
    let e0 = 2.0;
    let a = 2.0 * 0.7 * (1.0 + law.beta_per_k * e0);
    let b = 2.0 * 0.7 * law.beta_per_k;
    let target = 0.5;
    let exponential = (-a * target).exp();
    let exact = a * e0 * exponential / (a - b * e0 * (1.0 - exponential));
    let mut errors = Vec::new();
    for steps in [16, 32, 64, 128] {
        let mut controls = config();
        controls.thermal.max_dt_s = target / steps as f64;
        let candidate = advance_temperature(
            input,
            &initial,
            &(),
            &law,
            &[CellId(0)],
            target,
            controls,
            || false,
        )
        .unwrap();
        let first = candidate.states()[0]
            .center_of_mass_velocity_world(bodies.properties[0])
            .unwrap();
        let second = candidate.states()[1]
            .center_of_mass_velocity_world(bodies.properties[1])
            .unwrap();
        let energy = 0.5 * (first.x - second.x).powi(2);
        let heat = candidate.report().dissipated_energy_j;
        let temperature = candidate.thermal().temperature_k(CellId(0)).unwrap();
        near(energy + heat, e0, 1e-12);
        near(temperature - 300.0, heat, 7e-14);
        near(candidate.report().combined_energy_residual_j, 0.0, 1e-11);
        assert_eq!(candidate.temperatures_k(), &[temperature]);
        near(
            candidate.endpoint_links()[0].damping.native_coefficients()[0][0],
            0.7 * (1.0 + law.beta_per_k * (temperature - 300.0)),
            0.0,
        );
        assert_eq!(candidate.law_evaluations(), steps + 1);
        assert!(
            energy > exact,
            "lagged hardening must underdamp this reference"
        );
        errors.push(energy - exact);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.8..2.2).contains(&ratio),
            "first-order ladder: {errors:?}"
        );
    }
    assert!(
        errors[3] < 0.005,
        "fine-grid physical energy errors: {errors:?}"
    );
    assert_eq!(bodies.links[0].damping.native_coefficients()[0][0], 0.7);
    assert_eq!(initial.temperature_k(CellId(0)).unwrap(), 300.0);
}

#[test]
fn temperature_feedback_continuation_uses_original_reference_without_compounding() {
    let bodies = Bodies::new([0.7, 0.0]);
    let laws = linear_calorics();
    let initial = thermal(&laws, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let mut input = bodies.input(&routes);
    input.links = &bodies.links[..1];
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    let cells = [CellId(0)];
    let mut controls = config();
    controls.thermal.max_dt_s = 0.03125;
    let whole =
        advance_temperature(input, &initial, &(), &law, &cells, 0.5, controls, || false).unwrap();
    let first =
        advance_temperature(input, &initial, &(), &law, &cells, 0.25, controls, || false).unwrap();
    let continued = advance_temperature(
        RigidDampingThermalInput {
            time_s: first.report().time_s,
            states: first.states(),
            ..input
        },
        first.thermal(),
        first.observer(),
        &law,
        &cells,
        0.5,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(continued.states(), whole.states());
    assert_eq!(continued.thermal().cells(), whole.thermal().cells());
    assert_eq!(continued.endpoint_links(), whole.endpoint_links());
    assert_eq!(continued.temperatures_k(), whole.temperatures_k());
    near(
        first.report().dissipated_energy_j + continued.report().dissipated_energy_j,
        whole.report().dissipated_energy_j,
        5e-15,
    );
}

#[test]
fn explicit_temperature_coordinate_is_independent_of_heat_routes_and_event_clipped() {
    let bodies = Bodies::new([0.7, 0.0]);
    let laws = linear_calorics();
    let initial = DepositionSimulation::new(
        vec![
            birth(0, 0.0, 0.0),
            birth(1, 0.0, 0.0),
            TimedEvent {
                time_s: 0.35,
                event: Event::Boundary(CellId(1), boundary(3.0)),
            },
        ],
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        budget(),
    )
    .unwrap();
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let mut input = bodies.input(&routes);
    input.links = &bodies.links[..1];
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    let mut controls = config();
    controls.thermal.max_dt_s = 0.35;
    let observer = Observe::default();
    let candidate = advance_temperature(
        input,
        &initial,
        &observer,
        &law,
        &[CellId(1)],
        1.0,
        controls,
        || false,
    )
    .unwrap();
    // Heat goes only to cell 0; the explicitly distinct source coordinate is
    // cell 1 (capacity 1.5 J/K). Its boundary supplies 3 W only after t=.35.
    // Source temperatures are 300, 300 and 300.7 K on the three intervals.
    let integrated_c: f64 = 0.7 * (0.7 + (1.0 + 2.0 * 0.7) * 0.3);
    let expected_velocity = (-integrated_c).exp();
    near(
        candidate.states()[0]
            .center_of_mass_velocity_world(bodies.properties[0])
            .unwrap()
            .x,
        expected_velocity,
        2e-14,
    );
    near(candidate.temperatures_k()[0], 301.3, 6e-14);
    near(
        candidate.report().external_thermal_heat_j,
        3.0 * 0.65,
        2e-14,
    );
    near(candidate.report().combined_energy_residual_j, 0.0, 2e-12);
    assert_eq!(candidate.observer().intervals.len(), 3);
    assert_eq!(candidate.observer().intervals[0], (0.0, 0.35, 0.35));
    assert_eq!(candidate.observer().intervals[1], (0.35, 0.35, 0.7));
    assert_eq!(candidate.observer().events, [0.0, 0.0, 0.35]);
    assert_eq!(candidate.law_evaluations(), 3); // unchanged 300 K is reused
    assert_eq!(observer, Observe::default());
    assert!(initial.cells().is_empty());
}

#[test]
fn constant_temperature_and_zero_resistance_recover_original_native_results() {
    let laws = linear_calorics();
    let initial = thermal(&laws, false, true);
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    for zero in [false, true] {
        let mut bodies = Bodies::new(if zero { [0.0; 2] } else { [0.7, 1.1] });
        if zero {
            bodies.loads[0] = Wrench {
                force_world: Vec3::new(0.2, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            };
            bodies.gravity = Gravity::new(Vec3::new(0.0, -0.3, 0.0)).unwrap();
        }
        let routes = [BTreeMap::new(), BTreeMap::new()];
        let input = bodies.input(&routes);
        let observer = Observe {
            interval_watts: zero.then_some(0.2),
            ..Default::default()
        };
        let frozen =
            advance_candidate(input, &initial, &observer, 0.5, config(), || false).unwrap();
        let feedback = advance_temperature(
            input,
            &initial,
            &observer,
            &law,
            &[CellId(0), CellId(0)],
            0.5,
            config(),
            || false,
        )
        .unwrap();
        assert_eq!(feedback.states(), frozen.states());
        assert_eq!(feedback.thermal().cells(), frozen.thermal().cells());
        assert_eq!(feedback.observer(), frozen.observer());
        assert_eq!(
            feedback.report().dissipated_energy_j,
            frozen.report().dissipated_energy_j
        );
        assert_eq!(
            feedback.report().combined_energy_residual_j,
            frozen.report().combined_energy_residual_j
        );
        if zero {
            assert_eq!(feedback.report().pair_relaxations, 0);
            assert_eq!(feedback.report().dissipated_energy_j, 0.0);
        } else {
            assert_eq!(feedback.law_evaluations(), 2);
            assert_eq!(feedback.endpoint_links(), &bodies.links);
        }
    }
}

#[test]
fn temperature_source_endpoint_refusal_and_cancellation_leave_originals_retryable() {
    let bodies = Bodies::new([0.7, 0.0]);
    let laws = linear_calorics();
    let initial = thermal(&laws, false, true);
    let saved = initial.cells().clone();
    let observer = Observe::default();
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let mut input = bodies.input(&routes);
    input.links = &bodies.links[..1];
    let mut controls = config();
    controls.thermal.max_dt_s = 0.25;
    let narrow = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 300.01,
    };
    // The one whole interval is mechanically and calorically admissible. Its
    // final stored temperature lies outside the SOURCE domain, so the final
    // zero-work endpoint synchronization must reject all candidate histories.
    assert!(matches!(
        advance_temperature(
            input,
            &initial,
            &observer,
            &narrow,
            &[CellId(0)],
            0.25,
            controls,
            || false,
        ),
        Err(RigidDampingThermalError::Source {
            connector: 0,
            error: "synthetic source temperature domain"
        })
    ));
    assert_eq!(initial.cells(), &saved);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(observer, Observe::default());
    let broad = LinearTemperatureDamping {
        maximum_k: 303.0,
        ..narrow
    };
    let mut polls = 0;
    let expected = advance_temperature(
        input,
        &initial,
        &observer,
        &broad,
        &[CellId(0)],
        0.25,
        controls,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut repeated = 0;
    assert!(matches!(
        advance_temperature(
            input,
            &initial,
            &observer,
            &broad,
            &[CellId(0)],
            0.25,
            controls,
            || {
                repeated += 1;
                repeated == polls
            },
        ),
        Err(RigidDampingThermalError::Cancelled)
    ));
    for cells in [&[][..], &[CellId(99)][..]] {
        assert!(
            advance_temperature(
                input,
                &initial,
                &observer,
                &broad,
                cells,
                0.25,
                controls,
                || false,
            )
            .is_err()
        );
    }
    let retry = advance_temperature(
        input,
        &initial,
        &observer,
        &broad,
        &[CellId(0)],
        0.25,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(retry.observer(), expected.observer());
    assert_eq!(retry.endpoint_links(), expected.endpoint_links());
    assert_eq!(initial.cells(), &saved);
    assert_eq!(observer, Observe::default());
}
