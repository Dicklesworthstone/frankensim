//! Independent cooling oracle, material substitution and atomic physical refusals.
#[path = "heated_creeping_heat/common.rs"]
mod common;
use common::*;
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_material::gas::{GasSpec, GasState};
use fs_mbd::{Gravity, Vec3};
use fs_mbd_material::creeping_flow::CreepingSphere;
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;

#[test]
fn heating_and_cooling_converge_to_the_independent_spherical_conduction_solution() {
    let material = card(CP, 0.0, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    for (initial, reservoir) in [(300.0, 320.0), (320.0, 300.0)] {
        let reference = reference(&material, initial, RADIUS);
        let medium = fluid(reservoir, 10.0, None);
        let sphere =
            CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
        let conductance = 4.0 * std::f64::consts::PI * 10.0 * RADIUS;
        let capacity = reference.inertia().mass() * CP;
        let exact = reservoir + (initial - reservoir) * (-conductance / capacity).exp();
        let mut errors = Vec::new();
        for dt in [0.04, 0.02, 0.01] {
            let original = body(
                &reference,
                &material,
                &chart,
                CP,
                initial,
                RADIUS,
                dt,
                Gravity::ZERO,
                0.0,
            );
            let (next, report) = original
                .advance_creeping_sphere_in_thermal_fluid(
                    &sphere,
                    exchange(&solid),
                    flow_controls(),
                    1.0,
                    config(dt),
                    || Ok(()),
                )
                .unwrap();
            let actual = next.thermal().temperature_k(CellId(0)).unwrap();
            assert!(actual > initial.min(reservoir) && actual < initial.max(reservoir));
            errors.push((actual - exact).abs());
            close(report.before.conductance_w_k, conductance, 1e-12);
            close(
                report.reservoir_heat_into_body_j,
                capacity * (actual - initial),
                2e-8,
            );
            close(
                report.coupled.thermal.external_heat_j,
                report.reservoir_heat_into_body_j,
                1e-10,
            );
            close(report.coupled.combined_residual_j, 0.0, 1e-8);
            assert_eq!(report.coupled.viscous_heat_j, 0.0);
        }
        assert!(errors[1] < 0.52 * errors[0]);
        assert!(errors[2] < 0.52 * errors[1]);
        assert!(errors[2] < 1e-4);
    }
}

#[test]
fn material_heat_changes_expansion_inertia_rotation_and_later_drag_without_new_physics() {
    let mut outcomes = Vec::new();
    for (cp, expansion, fluid_k) in [
        (CP, 0.002, 10.0),
        (CP, 0.002, 20.0),
        (2.0 * CP, 0.002, 10.0),
        (CP, 0.0, 10.0),
    ] {
        let material = card(cp, expansion, 1000.0);
        let chart = caloric(&material);
        let solid = transport(&material, 1e5, 360.0);
        let reference = reference(&material, 300.0, RADIUS);
        let medium = fluid(320.0, fluid_k, None);
        let sphere =
            CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
        let original = body(
            &reference,
            &material,
            &chart,
            cp,
            300.0,
            RADIUS,
            0.01,
            Gravity::ZERO,
            0.2,
        );
        let (next, report) = original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                1.0,
                config(0.01),
                || Ok(()),
            )
            .unwrap();
        let temperature = next.thermal().temperature_k(CellId(0)).unwrap();
        let radius = report.after.radius_m;
        close(
            radius,
            RADIUS * (1.0 + expansion * (temperature - 300.0)),
            1e-12,
        );
        close(
            report.after.conductance_w_k / report.before.conductance_w_k,
            radius / RADIUS,
            1e-12,
        );
        close(
            next.configuration().inertia().mass(),
            reference.inertia().mass(),
            1e-12,
        );
        close(
            report.coupled.thermal.external_heat_j,
            report.reservoir_heat_into_body_j,
            1e-10,
        );
        assert!(report.coupled.viscous_heat_j > 0.0);
        close(report.coupled.combined_residual_j, 0.0, 1e-8);
        let omega = next
            .principal_frame()
            .properties()
            .angular_velocity_body_checked(next.state().angular_momentum_body())
            .unwrap()
            .z;
        outcomes.push((
            temperature,
            radius,
            omega,
            report.coupled.load_after.rotary_torque_body_nm.z,
        ));
    }
    assert!(
        outcomes[1].0 > outcomes[0].0
            && outcomes[1].1 > outcomes[0].1
            && outcomes[1].2 < outcomes[0].2
    );
    assert!(
        outcomes[2].0 < outcomes[0].0
            && outcomes[2].1 < outcomes[0].1
            && outcomes[2].2 > outcomes[0].2
    );
    assert_eq!(outcomes[3].1, RADIUS);
    assert!(outcomes[3].2 > outcomes[0].2);
    assert_ne!(outcomes[0].3, outcomes[3].3);
}

#[derive(Clone)]
struct Heater {
    watts: f64,
    steps: usize,
    refuse_at: usize,
}
impl DepositionObserver for Heater {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(CellId(0), self.watts)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.steps == self.refuse_at {
            return Err(DepositionError::InvalidInput(
                "late independent owner refusal",
            ));
        }
        Ok(())
    }
}

#[test]
fn source_heaters_share_substeps_and_late_failures_publish_no_partial_state() {
    let material = card(CP, 0.002, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let medium = fluid(320.0, 10.0, None);
    let sphere =
        CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        0.01,
        Gravity::ZERO,
        0.2,
    );
    let heater = Heater {
        watts: 100.0,
        steps: 0,
        refuse_at: usize::MAX,
    };
    let (next, advanced, report) = original
        .advance_creeping_sphere_in_thermal_fluid_with(
            &heater,
            &sphere,
            exchange(&solid),
            flow_controls(),
            0.2,
            config(0.01),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(advanced.steps, report.coupled.thermal.substeps);
    close(
        report.coupled.thermal.external_heat_j,
        report.reservoir_heat_into_body_j + 20.0,
        1e-9,
    );
    close(next.thermal().time_s(), 0.2, 0.0);
    let failing = Heater {
        refuse_at: 3,
        ..heater.clone()
    };
    let original_state = original.state();
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid_with(
                &failing,
                &sphere,
                exchange(&solid),
                flow_controls(),
                0.2,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(failing.steps, 0);
    assert_eq!(original.state(), original_state);
    assert_eq!(original.thermal().time_s(), 0.0);
    assert_eq!(original.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
    let limited_source = transport(&material, 1e5, 300.01);
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&limited_source),
                flow_controls(),
                1.0,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(original.state(), original_state);
    assert_eq!(original.thermal().time_s(), 0.0);
}

#[test]
fn ideal_gas_pressure_and_temperature_control_the_same_heat_and_motion_regimes() {
    let material = card(CP, 0.002, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let radius = 0.0002;
    let reference = reference(&material, 300.0, radius);
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap();
    let original = body(
        &reference, &material, &chart, CP, 300.0, radius, 1e-6, gravity, 0.0,
    );
    let mut samples = Vec::new();
    for pressure in [100_000.0, 200_000.0] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 301.0, pressure).unwrap();
        let sphere = CreepingSphere::new(&gas, radius, gravity, Vec3::ZERO, flow_limits()).unwrap();
        let (_, audit) = original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                1e-5,
                config(1e-6),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(audit.fluid.temperature_k(), gas.temperature);
        assert_eq!(audit.fluid.transport(), sphere.transport());
        assert_eq!(
            audit.fluid.volumetric_expansion_per_k(),
            Some(1.0 / gas.temperature)
        );
        assert!(audit.reservoir_heat_into_body_j > 0.0);
        samples.push(audit.before);
    }
    close(
        samples[1].absolute_rayleigh,
        4.0 * samples[0].absolute_rayleigh,
        1e-12,
    );
    close(
        samples[1].thermal_diffusion_ratio,
        2.0 * samples[0].thermal_diffusion_ratio,
        1e-12,
    );
    close(
        samples[1].conductance_w_k,
        samples[0].conductance_w_k,
        1e-15,
    );
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 301.0, 500_000.0).unwrap();
    let sphere = CreepingSphere::new(&gas, radius, gravity, Vec3::ZERO, flow_limits()).unwrap();
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                1e-5,
                config(1e-6),
                || Ok(()),
            )
            .is_err()
    );
}

#[test]
fn unrelated_sources_unknown_liquid_buoyancy_and_invalid_lumped_regimes_refuse() {
    let material = card(CP, 0.002, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let wrong_context = transport(&material, 2e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let medium = fluid(320.0, 10.0, None);
    let sphere =
        CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        0.01,
        Gravity::ZERO,
        0.0,
    );
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&wrong_context),
                flow_controls(),
                0.1,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    for (peclet, biot, diffusion, contrast) in [
        (0.1, 0.0001, 0.1, 0.1),
        (0.1, 0.1, 1e-12, 0.1),
        (0.1, 0.1, 0.1, 0.001),
    ] {
        let mut restrictive = exchange(&solid);
        restrictive.limits.maximum_peclet = peclet;
        restrictive.limits.maximum_biot = biot;
        restrictive.limits.maximum_diffusion_ratio = diffusion;
        restrictive.limits.maximum_relative_temperature_contrast = contrast;
        assert!(
            original
                .advance_creeping_sphere_in_thermal_fluid(
                    &sphere,
                    restrictive,
                    flow_controls(),
                    0.1,
                    config(0.01),
                    || Ok(()),
                )
                .is_err()
        );
    }
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.01)).unwrap();
    let falling = body(
        &reference, &material, &chart, CP, 300.0, RADIUS, 0.01, gravity, 0.0,
    );
    let unknown = CreepingSphere::new(&medium, RADIUS, gravity, Vec3::ZERO, flow_limits()).unwrap();
    assert!(
        falling
            .advance_creeping_sphere_in_thermal_fluid(
                &unknown,
                exchange(&solid),
                flow_controls(),
                0.1,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    let equilibrium = fluid(300.0, 10.0, None);
    let balanced =
        CreepingSphere::new(&equilibrium, RADIUS, gravity, Vec3::ZERO, flow_limits()).unwrap();
    let (_, report) = falling
        .advance_creeping_sphere_in_thermal_fluid(
            &balanced,
            exchange(&solid),
            flow_controls(),
            0.01,
            config(0.01),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.reservoir_heat_into_body_j, 0.0);
    let heater = Heater {
        watts: 100.0,
        steps: 0,
        refuse_at: usize::MAX,
    };
    assert!(
        falling
            .advance_creeping_sphere_in_thermal_fluid_with(
                &heater,
                &balanced,
                exchange(&solid),
                flow_controls(),
                0.01,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(falling.thermal().time_s(), 0.0);
}

#[test]
fn small_reynolds_does_not_hide_large_thermal_advection_or_duplicate_boundary_ownership() {
    let material = card(CP, 0.0, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let high_pr = fluid(320.0, 1e-6, None);
    let sphere =
        CreepingSphere::new(&high_pr, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        0.01,
        Gravity::ZERO,
        0.2,
    );
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                0.1,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    let medium = fluid(320.0, 10.0, None);
    let sphere =
        CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let duplicate = body_with_events(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        0.01,
        Gravity::ZERO,
        0.0,
        vec![TimedEvent {
            time_s: 0.03,
            event: Event::Boundary(
                CellId(0),
                Boundary {
                    ambient_k: 310.0,
                    htc_w_m2_k: 1.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            ),
        }],
    );
    assert!(
        duplicate
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                0.03,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(duplicate.thermal().time_s(), 0.0);
    let mut polls = 0;
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                1.0,
                config(0.01),
                || {
                    polls += 1;
                    if polls >= 150 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            )
            .is_err()
    );
    assert_eq!(original.thermal().time_s(), 0.0);
    assert_eq!(original.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
}

struct DisagreeingFluid(ThermalFluidTransportState);
impl FluidTransport for DisagreeingFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.0.transport_state()
    }
}
impl ThermalFluidTransport for DisagreeingFluid {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        Ok(fluid(self.0.temperature_k(), 20.0, None))
    }
}
#[test]
fn source_cannot_silently_change_transport_between_mechanical_and_thermal_bindings() {
    let material = card(CP, 0.0, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let source = DisagreeingFluid(fluid(320.0, 10.0, None));
    let sphere =
        CreepingSphere::new(&source, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        0.01,
        Gravity::ZERO,
        0.0,
    );
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                1.0,
                config(0.01),
                || Ok(()),
            )
            .is_err()
    );
}
