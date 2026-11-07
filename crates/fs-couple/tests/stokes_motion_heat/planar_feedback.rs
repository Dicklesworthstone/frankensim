//! G1/G3/G4: evolved fluid temperature reaches planar force, torque and heat.

use super::{
    CP, DEPTH, MU, RHO, SPEED, WALL, Caloric, ClockedHeat, CoupledHeatRefusal,
    binding, fluid, initial_thermal, initial_thermal_with_events, near, slip_exterior,
    slip_routes, thermal_config, top_edges, wall_cell, zero,
};
use super::planar::planar_motion;
use fs_couple::stokes_motion::heating::{
    PlanarStokesHeatError, PlanarStokesObservedHeatError,
};
use fs_couple::stokes_motion::planar::{
    PlanarStokesControls, PlanarStokesError, PlanarStokesMotion,
};
use fs_flux::{FluxSystem, MaterialStokesBoundary, TriMesh};
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{DepositionSimulation, Event, TimedEvent};

/// Authored numerical law with retained query metadata; no measured-fluid claim.
struct TemperatureSource {
    query_k: f64,
    label: &'static str,
    state: ThermalFluidTransportState,
}
impl FluidTransport for TemperatureSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermalFluidTransport for TemperatureSource {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        Ok(self.state)
    }
}
fn source(temperature_k: f64) -> Result<TemperatureSource, MaterialError> {
    let viscosity = MU * (1.0 + 0.1 * (temperature_k - 300.0));
    let transport = FluidTransportState::try_new(RHO, CP, 0.5, viscosity)?;
    Ok(TemperatureSource {
        query_k: temperature_k,
        label: "authored temperature-dependent viscosity",
        state: ThermalFluidTransportState::try_new(transport, temperature_k, Some(0.0))?,
    })
}

#[test]
fn g1_planar_temperature_changes_resolved_shear_and_preserves_zero_heat_rigid_rotation() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let reference = fluid();
    let edges = top_edges(&mesh);
    let law = Caloric { maximum_h: 100.0 };
    for temperature in [300.0, 308.0] {
        for omega in [-SPEED, 0.12] {
            let field = |x: [f64; 2]| [SPEED * x[1], omega * (x[0] - 0.5)];
            let mu = source(temperature).unwrap().state.transport().dynamic_viscosity_pa_s();
            let stress = mu * (SPEED + omega);
            let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|edge| mesh.is_boundary(*edge) && !edges.contains(edge))
                .map(|edge| {
                    if mesh.edges[edge].normal[1] < -0.5 {
                        MaterialStokesBoundary::Velocity { edge, value: &field }
                    } else {
                        MaterialStokesBoundary::TotalTraction { edge, value: &traction }
                    }
                })
                .collect();
            let mut owner = planar_motion(&system, &reference, &edges, field, [SPEED, 0.0, omega]);
            let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, cell| {
                cell.specific_h_j_kg = CP * (temperature - 300.0);
            });
            let step = owner.step_with_temperature_feedback(
                &mut thermal, binding(&ids), &[], 0.0, source, 0.05, &zero,
                &boundaries, &[None], [DEPTH * mu * (SPEED + omega), 0.0], 0.0,
                PlanarStokesControls::default(), thermal_config(), || false,
            ).unwrap();
            assert_eq!(step.material.query_k, temperature);
            assert_eq!(step.material.label, "authored temperature-dependent viscosity");
            assert_eq!(step.sample.temperature_k, temperature);
            assert_eq!(step.transport.transport().dynamic_viscosity_pa_s(), mu);
            assert_eq!(owner.material(), &reference, "ordinary binding remains explicit");
            near(step.heat.motion.velocity_m_s[0], SPEED, 2e-8);
            near(step.heat.motion.angular_velocity_rad_s, omega, 2e-8);
            near(step.heat.motion.fluid_force_n[0], -DEPTH * mu * (SPEED + omega), 2e-9);
            let heat_density = 0.05 * mu * (SPEED + omega).powi(2);
            near(step.heat.deposited_viscous_heat_j, DEPTH * heat_density, 2e-10);
            near(step.heat.thermal.stored_energy_change_j, DEPTH * heat_density, 2e-10);
            for id in ids {
                near(thermal.temperature_k(id).unwrap(),
                    temperature + heat_density / (RHO * CP), 2e-9);
            }
            assert_eq!(owner.state().time_s(), thermal.time_s());
        }
    }
}

#[test]
fn g4_planar_temperature_feedback_reuses_native_slip_and_observer_transaction() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = fluid();
    let edges = top_edges(&mesh);
    let law = Caloric { maximum_h: 100.0 };
    let temperature = 305.0;
    let warm = source(temperature).unwrap().state.transport();
    let mu = warm.dynamic_viscosity_pa_s();
    // Fixed exterior slip resistance beta=MU, while the fluid viscosity changes.
    // tau=mu*gamma=beta*(SPEED-gamma) independently fixes both losses.
    let gamma = MU * SPEED / (mu + MU);
    let slip = SPEED - gamma;
    let field = |x: [f64; 2]| [slip + gamma * x[1], 0.0];
    let stress = mu * gamma;
    let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
    let boundaries = slip_exterior(&mesh, &edges, MU, &zero, &traction);
    let routes = slip_routes(&boundaries, 0.25);
    let (original, ids) = initial_thermal_with_events(
        &mesh, &law, true,
        |_, cell| { cell.specific_h_j_kg = CP * (temperature - 300.0); },
        vec![TimedEvent { time_s: 0.0, event: Event::Birth(wall_cell(DEPTH)) }],
    );
    let advance = |
        owner: &mut PlanarStokesMotion<'_, '_, '_, FluidTransportState>,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut ClockedHeat,
        cancelled: &mut dyn FnMut() -> bool,
    | owner.step_with_temperature_feedback_observed(
        thermal, observer, binding(&ids), &routes, 0.0, source, 0.05, &zero,
        &boundaries, &[None], [DEPTH * stress, 0.0], 0.0,
        PlanarStokesControls::default(), thermal_config(), cancelled,
    );

    let mut owner = planar_motion(&system, &reference, &edges, field, [SPEED, 0.0, 0.0]);
    let mut thermal = original.clone();
    let mut observer = ClockedHeat::new(ids[0]);
    let mut polls = 0;
    let accepted = advance(&mut owner, &mut thermal, &mut observer, &mut || {
        polls += 1;
        false
    }).unwrap();
    let bulk_heat = 0.05 * DEPTH * mu * gamma * gamma;
    let slip_heat = 0.05 * DEPTH * MU * slip * slip;
    near(accepted.heat.deposited_viscous_heat_j, bulk_heat, 2e-10);
    near(accepted.heat.deposited_slip_heat_j, slip_heat, 2e-10);
    near(thermal.cells()[&WALL].specific_h_j_kg,
        0.75 * slip_heat / thermal.cells()[&WALL].mass_kg, 2e-9);

    // Exact parity with the ordinary native path bound to that same warm source.
    let mut direct = planar_motion(&system, &warm, &edges, field, [SPEED, 0.0, 0.0]);
    let mut direct_thermal = original.clone();
    let mut direct_observer = ClockedHeat::new(ids[0]);
    direct.step_with_slip_heating_observed(
        &mut direct_thermal, &mut direct_observer, binding(&ids), &routes, 0.05,
        &zero, &boundaries, &[None], [DEPTH * stress, 0.0], 0.0,
        PlanarStokesControls::default(), thermal_config(), || false,
    ).unwrap();
    assert_eq!(owner.state(), direct.state());
    assert_eq!(thermal.cells(), direct_thermal.cells());
    assert_eq!(observer, direct_observer);

    for late_observer_refusal in [true, false] {
        let mut retry = planar_motion(&system, &reference, &edges, field, [SPEED, 0.0, 0.0]);
        let before = retry.state().clone();
        let mut retry_thermal = original.clone();
        let mut retry_observer = ClockedHeat::new(ids[0]);
        retry_observer.refuse_at_s = late_observer_refusal.then_some(0.03);
        let original_observer = retry_observer.clone();
        let mut calls = 0;
        let refused = advance(&mut retry, &mut retry_thermal, &mut retry_observer, &mut || {
            calls += 1;
            !late_observer_refusal && calls == polls
        });
        if late_observer_refusal {
            assert!(matches!(refused, Err(PlanarStokesObservedHeatError::Observer(
                CoupledHeatRefusal::Endpoint
            ))));
        } else {
            assert!(matches!(refused, Err(PlanarStokesObservedHeatError::Heat(
                PlanarStokesHeatError::Motion(PlanarStokesError::Cancelled)
            ))));
        }
        assert_eq!(retry.state(), &before);
        assert_eq!(retry_thermal.cells(), original.cells());
        assert_eq!(retry_thermal.links(), original.links());
        assert_eq!(retry_thermal.time_s(), original.time_s());
        assert_eq!(retry_observer, original_observer);
        retry_observer.refuse_at_s = None;
        advance(&mut retry, &mut retry_thermal, &mut retry_observer, &mut || false).unwrap();
        assert_eq!(retry.state(), owner.state());
        assert_eq!(retry_thermal.cells(), thermal.cells());
        assert_eq!(retry_thermal.time_s(), thermal.time_s());
        assert_eq!(retry_observer, observer);
    }
}
