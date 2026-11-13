//! Start-state material feedback uses the original momentum and thermal owners.

use super::*;
use std::cell::Cell;

use fs_flux::MaterialFluxError;
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};

/// Explicit authored numerical source, not a calibrated real-fluid claim.
/// Keeping this owner in the receipt preserves its query and source identifier.
struct TemperatureSource {
    source_key: &'static str,
    query_k: f64,
    momentum: FluidTransportState,
    thermal: ThermalFluidTransportState,
    momentum_reads: Cell<usize>,
    thermal_reads: Cell<usize>,
}
impl FluidTransport for TemperatureSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        assert_eq!(self.momentum_reads.get(), 0, "momentum source queried twice");
        self.momentum_reads.set(1);
        Ok(self.momentum)
    }
}
impl ThermalFluidTransport for TemperatureSource {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        assert_eq!(self.thermal_reads.get(), 0, "thermal source queried twice");
        self.thermal_reads.set(1);
        Ok(self.thermal)
    }
}
fn source(temperature_k: f64, viscosity: f64) -> TemperatureSource {
    let momentum = FluidTransportState::try_new(RHO, CP, 0.5, viscosity).unwrap();
    TemperatureSource {
        source_key: "manufactured-viscosity-curve",
        query_k: temperature_k,
        momentum,
        thermal: ThermalFluidTransportState::try_new(momentum, temperature_k, None).unwrap(),
        momentum_reads: Cell::new(0),
        thermal_reads: Cell::new(0),
    }
}
fn viscosity(temperature_k: f64) -> f64 {
    MU / (1.0 + (temperature_k - 300.0) / 20.0)
}
fn boundaries<'a>(
    mesh: &TriMesh,
    body: &[usize],
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
) -> Vec<MaterialStokesBoundary<'a>> {
    (0..mesh.edges.len())
        .filter(|e| mesh.is_boundary(*e) && !body.contains(e))
        .map(|edge| {
            if mesh.edges[edge].normal[1] < -0.5 {
                MaterialStokesBoundary::Velocity { edge, value: &zero }
            } else {
                MaterialStokesBoundary::TotalTraction { edge, value: traction }
            }
        })
        .collect()
}

#[test]
fn g1_temperature_source_changes_native_drag_and_heat_then_resamples_accepted_temperature() {
    let law = Caloric { maximum_h: 1000.0 };
    for n in [1, 2] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let baseline = fluid();
        let edges = top_edges(&mesh);
        let mut first_loads = Vec::new();
        let mut first_heats = Vec::new();
        for start_k in [300.0, 320.0] {
            let mut owner = motion(&system, &baseline, &edges, config(), shear, SPEED);
            let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, cell| {
                cell.specific_h_j_kg = CP * (start_k - 300.0);
            });
            let queries = Cell::new(0);
            let mut expected_k = start_k;
            let mut time_s = 0.0;
            for (interval, dt) in [0.02, 0.03].into_iter().enumerate() {
                let mu = viscosity(expected_k);
                let stress = mu * SPEED;
                let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
                let exterior = boundaries(&mesh, &edges, &traction);
                let step = owner.step_with_temperature_feedback(
                    &mut thermal, binding(&ids), &[], 1e-8,
                    |temperature| {
                        queries.set(queries.get() + 1);
                        Ok(source(temperature, viscosity(temperature)))
                    },
                    dt, &zero, &exterior, &[None], DEPTH * stress,
                    GuidedStokesControls::default(), thermal_config(), || false,
                ).unwrap();
                assert_eq!(step.material.source_key, "manufactured-viscosity-curve");
                assert_eq!(step.material.query_k, step.sample.temperature_k);
                assert_eq!(step.material.momentum_reads.get(), 1);
                assert_eq!(step.material.thermal_reads.get(), 1);
                assert_eq!(step.transport.transport(), *step.heat.motion.flow.transport());
                near(step.sample.temperature_k, expected_k, 2e-9);
                near(step.sample.minimum_temperature_k, expected_k, 2e-9);
                near(step.sample.maximum_temperature_k, expected_k, 2e-9);
                near(step.heat.motion.speed_m_s, SPEED, 2e-9);
                near(step.heat.motion.fluid_force_n, -DEPTH * stress, 2e-9);
                near(step.heat.deposited_viscous_heat_j, dt * DEPTH * mu * SPEED * SPEED, 2e-10);
                near(step.heat.deposited_slip_heat_j, 0.0, 1e-14);
                near(step.heat.thermal.external_heat_j, step.heat.deposited_viscous_heat_j, 1e-12);
                if interval == 0 {
                    first_loads.push(step.heat.motion.fluid_force_n);
                    first_heats.push(step.heat.deposited_viscous_heat_j);
                }
                expected_k += dt * mu * SPEED * SPEED / (RHO * CP);
                time_s += dt;
                for &id in &ids {
                    near(thermal.temperature_k(id).unwrap(), expected_k, 2e-9);
                }
                near(owner.state().displacement_m(), time_s * SPEED, 2e-9);
                assert_eq!(owner.state().time_s(), time_s);
                assert_eq!(thermal.time_s(), time_s);
                // Feedback does not silently replace the source used by ordinary steps.
                assert_eq!(owner.material(), &baseline);
            }
            assert_eq!(queries.get(), 2);
        }
        // The 320 K source has half the 300 K viscosity. Exact Couette strain
        // therefore has half the actual native wall load and physical heat.
        near(first_loads[1] / first_loads[0], 0.5, 2e-8);
        near(first_heats[1] / first_heats[0], 0.5, 2e-8);
    }
}

#[test]
fn g3_constant_feedback_matches_frozen_flow_slip_and_every_observer_history() {
    let law = Caloric { maximum_h: 100.0 };
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let baseline = fluid();
    let edges = top_edges(&mesh);
    for with_slip in [false, true] {
        let gamma = if with_slip { SPEED / 2.0 } else { SPEED };
        let slip = SPEED - gamma;
        let field = |x: [f64; 2]| [slip + gamma * x[1], 0.0];
        let stress = MU * gamma;
        let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
        let exterior = if with_slip {
            slip_exterior(&mesh, &edges, MU, &zero, &traction)
        } else {
            boundaries(&mesh, &edges, &traction)
        };
        let routes = slip_routes(&exterior, 0.25);
        let (original, ids) = initial_thermal_with_events(
            &mesh, &law, true,
            |i, cell| { if i == 0 { cell.specific_h_j_kg = 1.0; } },
            vec![TimedEvent { time_s: 0.0, event: Event::Birth(wall_cell(DEPTH)) }],
        );
        let mut frozen = motion(&system, &baseline, &edges, config(), field, SPEED);
        let mut frozen_thermal = original.clone();
        let mut frozen_observer = ClockedHeat::new(ids[0]);
        let expected = frozen.step_with_slip_heating_observed(
            &mut frozen_thermal, &mut frozen_observer, binding(&ids), &routes,
            0.05, &zero, &exterior, &[None], DEPTH * stress,
            GuidedStokesControls::default(), thermal_config(), || false,
        ).unwrap();
        let mut refreshed = motion(&system, &baseline, &edges, config(), field, SPEED);
        let mut thermal = original.clone();
        let mut observer = ClockedHeat::new(ids[0]);
        let step = refreshed.step_with_temperature_feedback_observed(
            &mut thermal, &mut observer, binding(&ids), &routes, 1.0,
            |temperature| Ok(source(temperature, MU)),
            0.05, &zero, &exterior, &[None], DEPTH * stress,
            GuidedStokesControls::default(), thermal_config(), || false,
        ).unwrap();
        assert_eq!(refreshed.state(), frozen.state());
        assert_eq!(thermal.cells(), frozen_thermal.cells());
        assert_eq!(thermal.links(), frozen_thermal.links());
        assert_eq!(thermal.time_s(), frozen_thermal.time_s());
        assert_eq!(observer, frozen_observer);
        assert_eq!(step.heat.thermal, expected.thermal);
        assert_eq!(step.heat.deposited_viscous_heat_j, expected.deposited_viscous_heat_j);
        assert_eq!(step.heat.deposited_slip_heat_j, expected.deposited_slip_heat_j);
        assert_eq!(step.sample.minimum_temperature_k, 300.0);
        assert_eq!(step.sample.maximum_temperature_k, 300.25);
        assert_eq!(step.sample.temperature_k, 300.125);
        assert_eq!(step.material.momentum_reads.get(), 1);
        assert_eq!(step.material.thermal_reads.get(), 1);
    }
}

#[derive(Clone, Copy, Debug)]
enum BadSource {
    Spread,
    NonfiniteSpread,
    Density,
    Temperature,
    Views,
    Refusal,
}

#[test]
fn g4_source_and_temperature_admission_preserve_all_owners_before_flow() {
    let law = Caloric { maximum_h: 100.0 };
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let baseline = fluid();
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let (original, ids) = initial_thermal(&mesh, &law, true, |i, cell| {
        if i == 0 { cell.specific_h_j_kg = 1.0; }
    });
    for case in [
        BadSource::Spread, BadSource::NonfiniteSpread, BadSource::Density,
        BadSource::Temperature, BadSource::Views, BadSource::Refusal,
    ] {
        let mut owner = motion(&system, &baseline, &edges, config(), shear, SPEED);
        let old_motion = owner.state().clone();
        let mut thermal = original.clone();
        let mut observer = ClockedHeat::new(ids[0]);
        let old_observer = observer.clone();
        let queries = Cell::new(0);
        let force_calls = Cell::new(0);
        let force = |_: [f64; 2]| {
            force_calls.set(force_calls.get() + 1);
            [0.0; 2]
        };
        let spread = match case {
            BadSource::Spread => 0.1,
            BadSource::NonfiniteSpread => f64::NAN,
            _ => 1.0,
        };
        let result = owner.step_with_temperature_feedback_observed(
            &mut thermal, &mut observer, binding(&ids), &[], spread,
            |temperature| {
                queries.set(queries.get() + 1);
                let mut material = source(temperature, MU);
                match case {
                    BadSource::Density => {
                        // Even one representable density increment must not
                        // silently alter the stored incompressible momentum.
                        let rho = f64::from_bits(RHO.to_bits() + 1);
                        material.momentum = FluidTransportState::try_new(rho, CP, 0.5, MU).unwrap();
                        material.thermal = ThermalFluidTransportState::try_new(
                            material.momentum, temperature, None,
                        ).unwrap();
                    }
                    BadSource::Temperature => material = source(temperature + 1.0, MU),
                    BadSource::Views => {
                        material.thermal = ThermalFluidTransportState::try_new(
                            FluidTransportState::try_new(RHO, CP, 0.5, 2.0 * MU).unwrap(),
                            temperature, None,
                        ).unwrap();
                    }
                    BadSource::Refusal => return Err(MaterialError::State {
                        what: "manufactured source domain".into(),
                    }),
                    _ => {}
                }
                Ok(material)
            },
            0.05, &force, &exterior, &[None], DEPTH * MU * SPEED,
            GuidedStokesControls::default(), thermal_config(), || false,
        );
        if matches!(case, BadSource::Refusal) {
            assert!(matches!(result, Err(GuidedStokesObservedHeatError::Heat(
                GuidedStokesHeatError::Motion(GuidedStokesError::Flow(
                    MaterialFluxError::Material(MaterialError::State { .. })
                ))
            ))));
        } else {
            assert!(matches!(result, Err(GuidedStokesObservedHeatError::Heat(
                GuidedStokesHeatError::Binding(_)
            ))), "{case:?}");
        }
        let expected_queries = if matches!(case, BadSource::Spread | BadSource::NonfiniteSpread) {
            0
        } else {
            1
        };
        assert_eq!(queries.get(), expected_queries, "{case:?}");
        assert_eq!(force_calls.get(), 0, "{case:?}");
        assert_eq!(owner.state(), &old_motion, "{case:?}");
        assert_eq!(thermal.time_s(), original.time_s(), "{case:?}");
        assert_eq!(thermal.cells(), original.cells(), "{case:?}");
        assert_eq!(thermal.links(), original.links(), "{case:?}");
        assert_eq!(observer, old_observer, "{case:?}");
    }
}

#[test]
fn g4_feedback_discards_late_observer_refusal_and_final_cancel_then_retries_once() {
    let law = Caloric { maximum_h: 100.0 };
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let baseline = fluid();
    let edges = top_edges(&mesh);
    let (original, ids) = initial_thermal(&mesh, &law, true, |i, cell| {
        if i == 0 { cell.specific_h_j_kg = 1.0; }
    });
    let stress = viscosity(300.125) * SPEED;
    let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
    let exterior = boundaries(&mesh, &edges, &traction);
    let advance = |
        owner: &mut GuidedStokesMotion<'_, '_, '_, FluidTransportState>,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut ClockedHeat,
        cancelled: &mut dyn FnMut() -> bool,
    | {
        owner.step_with_temperature_feedback_observed(
            thermal, observer, binding(&ids), &[], 1.0,
            |temperature| Ok(source(temperature, viscosity(temperature))),
            0.05, &zero, &exterior, &[None], DEPTH * stress,
            GuidedStokesControls::default(), thermal_config(), cancelled,
        )
    };
    let mut control = motion(&system, &baseline, &edges, config(), shear, SPEED);
    let mut control_thermal = original.clone();
    let mut control_observer = ClockedHeat::new(ids[0]);
    let mut polls = 0;
    let accepted = advance(
        &mut control, &mut control_thermal, &mut control_observer, &mut || {
            polls += 1;
            false
        },
    ).unwrap();
    assert!(control_observer.starts_and_durations.len() > 3);
    for refuse in [true, false] {
        let mut retry = motion(&system, &baseline, &edges, config(), shear, SPEED);
        let old_motion = retry.state().clone();
        let mut thermal = original.clone();
        let mut observer = ClockedHeat::new(ids[0]);
        observer.refuse_at_s = refuse.then_some(0.03);
        let old_observer = observer.clone();
        let mut attempts = 0;
        let result = advance(&mut retry, &mut thermal, &mut observer, &mut || {
            attempts += 1;
            !refuse && attempts == polls
        });
        if refuse {
            assert!(matches!(result, Err(GuidedStokesObservedHeatError::Observer(
                CoupledHeatRefusal::Endpoint
            ))));
        } else {
            assert!(matches!(result, Err(GuidedStokesObservedHeatError::Heat(
                GuidedStokesHeatError::Motion(GuidedStokesError::Cancelled)
            ))));
        }
        assert_eq!(retry.state(), &old_motion);
        assert_eq!(thermal.time_s(), original.time_s());
        assert_eq!(thermal.cells(), original.cells());
        assert_eq!(thermal.links(), original.links());
        assert_eq!(observer, old_observer);
        observer.refuse_at_s = None;
        let step = advance(&mut retry, &mut thermal, &mut observer, &mut || false).unwrap();
        assert_eq!(retry.state(), control.state());
        assert_eq!(thermal.cells(), control_thermal.cells());
        assert_eq!(thermal.time_s(), control_thermal.time_s());
        assert_eq!(observer, control_observer);
        assert_eq!(step.heat.thermal, accepted.heat.thermal);
        assert_eq!(step.sample, accepted.sample);
    }
}

#[test]
fn g4_feedback_caloric_refusal_keeps_initial_material_and_short_retry_deposits_once() {
    let law = Caloric { maximum_h: 0.0002 };
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let baseline = fluid();
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let mut owner = motion(&system, &baseline, &edges, config(), shear, SPEED);
    let old_motion = owner.state().clone();
    let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
    let original = thermal.clone();
    let queries = Cell::new(0);
    let resolve = |temperature| {
        queries.set(queries.get() + 1);
        assert_eq!(temperature, 300.0);
        Ok(source(temperature, viscosity(temperature)))
    };
    let refused = owner.step_with_temperature_feedback(
        &mut thermal, binding(&ids), &[], 1e-8, resolve,
        0.05, &zero, &exterior, &[None], DEPTH * MU * SPEED,
        GuidedStokesControls::default(), thermal_config(), || false,
    );
    assert!(matches!(refused, Err(GuidedStokesHeatError::Thermal(
        DepositionError::MaterialDomain
    ))));
    assert_eq!(owner.state(), &old_motion);
    assert_eq!(thermal.cells(), original.cells());
    assert_eq!(thermal.time_s(), 0.0);
    let step = owner.step_with_temperature_feedback(
        &mut thermal, binding(&ids), &[], 1e-8, resolve,
        0.01, &zero, &exterior, &[None], DEPTH * MU * SPEED,
        GuidedStokesControls::default(), thermal_config(), || false,
    ).unwrap();
    assert_eq!(queries.get(), 2);
    assert_eq!(step.material.momentum_reads.get(), 1);
    assert_eq!(step.material.thermal_reads.get(), 1);
    assert_eq!(owner.state().time_s(), 0.01);
    assert_eq!(thermal.time_s(), 0.01);
    near(step.heat.deposited_viscous_heat_j, 0.01 * DEPTH * MU * SPEED * SPEED, 2e-10);
    for id in ids {
        near(thermal.cells()[&id].specific_h_j_kg, 0.01 * MU * SPEED * SPEED / RHO, 2e-9);
    }
}
