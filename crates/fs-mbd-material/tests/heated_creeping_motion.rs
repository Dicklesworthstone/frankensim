//! G1/G3/G4: source expansion reaches the original viscous trajectory.
//! Synthetic coefficients exercise the admitted model, not measured materials.
#[path = "heating/common.rs"]
mod common;
use common::{budget, caloric, card, config, point, region};
use fs_ga::{Se3, So3, Twist, Vec3};
use fs_matdb::MaterialCard;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::creeping_flow::{
    CreepingFlowError, CreepingFlowLimits, CreepingSphere, CreepingStepControls,
};
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::HeatedBody;
use fs_mbd_material::thermal::heating::motion::creeping_flow::HeatedCreepingControls;
use fs_mbd_material::thermal::heating::motion::{
    HeatedMotionControls, HeatedMotionError, HeatedRigidBody,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;

const RADIUS: f64 = 0.1;
const DT: f64 = 0.03125;
const CP: f64 = 10.0;
const RATE: f64 = 16.0;
const DENSITY: f64 = 2.0;

fn close(actual: f64, expected: f64) {
    let tolerance = 5e-15 + 2e-12 * actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn gravity() -> Gravity {
    Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -0.1)).unwrap()
}
fn controls() -> HeatedMotionControls {
    HeatedMotionControls {
        maximum_step_s: DT,
        relative_tensor_tolerance: 1e-12,
        transfer: ConfigurationTransferControls {
            relative_momentum_tolerance: 1e-12,
            maximum_parameter_work_j: 1.0,
        },
        maximum_numerical_energy_change_j: 1e-12,
    }
}
fn flow_controls() -> HeatedCreepingControls {
    HeatedCreepingControls {
        maximum_relaxation_fraction: 0.5,
        maximum_step_energy_defect_j: 1e-13,
    }
}
fn limits() -> CreepingFlowLimits {
    CreepingFlowLimits {
        maximum_reynolds: 0.05,
        maximum_diffusion_ratio: 0.05,
    }
}
fn flow() -> fs_mbd::Vec3 {
    fs_mbd::Vec3::new(0.02, 0.0, 0.0)
}
fn twist() -> Twist {
    Twist::new(Vec3::new(0.0, 0.0, 0.2), Vec3::new(0.01, -0.02, 0.03))
}
fn fluid(rho: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 1000.0, 0.1, mu).unwrap()
}
fn reference_body(material: &MaterialCard, geometry: VolumeMoments) -> ResolvedBody {
    resolve_body(
        &[RegionBinding {
            key: "a",
            material,
            point: &point(300.0, 1e5),
            geometry,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn body<'a>(
    reference: &'a ResolvedBody,
    material: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    rate: f64,
    control: HeatedMotionControls,
    gravity: Gravity,
    twist: Twist,
) -> HeatedRigidBody<'a> {
    let mut thermal = region(
        material,
        chart,
        "a",
        0,
        reference.inertia().mass() * CP * rate,
    );
    thermal.geometry.volume_m3 = reference.regions()[0].geometry().volume_m3();
    thermal.geometry.surface_area_m2 = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    thermal.geometry.start_m = [0.0, 0.0, -RADIUS];
    thermal.geometry.end_m = [0.0, 0.0, RADIUS];
    // This reference-coordinate anchor is deliberately off centre. The
    // original transfer must still preserve the world centre and momenta.
    thermal.fixed_point_m = Vec3::new(0.02, -0.01, 0.03);
    let thermal = HeatedBody::new(reference, vec![thermal], vec![], budget(), config(), || {
        Ok(())
    })
    .unwrap();
    HeatedRigidBody::new(
        thermal,
        Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, -0.2)).unwrap(),
        twist,
        RigidBodyIntegrator::new(gravity),
        control,
        || Ok(()),
    )
    .unwrap()
}

/// Independent midpoint solution of m v'=m g-rho_f V g-c(v-U),
/// L'=-c_rot L/I, followed by the declared fixed-momentum radius change.
struct Reference {
    velocity: [f64; 3],
    position: [f64; 3],
    spin: f64,
    heat: f64,
    flow_work: f64,
    inertia_work: f64,
    displacement_work: f64,
}
fn reference(alpha: f64, rho: f64, mu: f64, steps: usize) -> Reference {
    let volume = |r: f64| 4.0 * std::f64::consts::PI * r * r * r / 3.0;
    let mass = DENSITY * volume(RADIUS);
    let inertia = |r: f64| 0.4 * mass * r * r;
    let mut result = Reference {
        velocity: [0.01, -0.02, 0.03],
        position: [0.0, 0.0, -0.2],
        spin: inertia(RADIUS) * 0.2,
        heat: 0.0,
        flow_work: 0.0,
        inertia_work: 0.0,
        displacement_work: 0.0,
    };
    for i in 0..steps {
        let radius = RADIUS * (1.0 + alpha * RATE * (i as f64 * DT));
        let next_radius = RADIUS * (1.0 + alpha * RATE * ((i + 1) as f64 * DT));
        let c = 6.0 * std::f64::consts::PI * mu * radius;
        let rotary = 8.0 * std::f64::consts::PI * mu * radius.powi(3);
        let fluid_mass = rho * volume(radius);
        let angular_mid =
            (result.spin / inertia(radius)) / (1.0 + 0.5 * DT * rotary / inertia(radius));
        for axis in 0..3 {
            let ambient = [0.02, 0.0, 0.0][axis];
            let acceleration = [0.0, 0.0, -0.1][axis] * (1.0 - fluid_mass / mass);
            let relative_mid = (result.velocity[axis] - ambient + 0.5 * DT * acceleration)
                / (1.0 + 0.5 * DT * c / mass);
            let midpoint = relative_mid + ambient;
            result.position[axis] += DT * midpoint;
            result.velocity[axis] = 2.0 * midpoint - result.velocity[axis];
            result.heat += DT * c * relative_mid * relative_mid;
            result.flow_work += -DT * c * relative_mid * ambient;
        }
        result.spin -= DT * rotary * angular_mid;
        result.heat += DT * rotary * angular_mid * angular_mid;
        result.inertia_work +=
            0.5 * result.spin * result.spin * (1.0 / inertia(next_radius) - 1.0 / inertia(radius));
        result.displacement_work +=
            rho * (volume(next_radius) - volume(radius)) * (-0.1 * result.position[2]);
    }
    result
}

#[test]
fn source_expansion_viscosity_and_density_reach_actual_motion_and_physical_accounts() {
    let mut outcomes = Vec::new();
    for (alpha, rho, mu) in [
        (0.0, 0.002, 0.001),
        (0.02, 0.002, 0.001),
        (0.02, 0.002, 0.002),
        (0.02, 0.004, 0.001),
    ] {
        let material = card(CP, alpha, 500.0);
        let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
        let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
        let body = body(
            &source,
            &material,
            &chart,
            RATE,
            controls(),
            gravity(),
            twist(),
        );
        let reservoir = fluid(rho, mu);
        let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
        let (next, audit) = body
            .advance_creeping_sphere(&sphere, flow_controls(), 4.0 * DT, config(), || Ok(()))
            .unwrap();
        let expected = reference(alpha, rho, mu, 4);
        let v = next
            .state()
            .center_of_mass_velocity_world(next.principal_frame().properties())
            .unwrap();
        let x = next.state().pose().position_world();
        let spin = next
            .state()
            .pose()
            .orientation()
            .rotate_body_to_world(next.state().angular_momentum_body());
        for (actual, value) in [v.x, v.y, v.z, x.x, x.y, x.z, spin.z].into_iter().zip(
            expected
                .velocity
                .into_iter()
                .chain(expected.position)
                .chain([expected.spin]),
        ) {
            close(actual, value);
        }
        let radius = RADIUS * (1.0 + alpha * RATE * 4.0 * DT);
        close(audit.radius_after_m, radius);
        close(audit.viscous_heat_j, expected.heat);
        close(audit.background_work_j, expected.flow_work);
        close(
            audit.displacement_parameter_work_j,
            expected.displacement_work,
        );
        close(
            audit.parameter_work_j,
            expected.inertia_work + expected.displacement_work,
        );
        assert!(audit.absolute_parameter_work_j + 1e-16 >= audit.parameter_work_j.abs());
        close(audit.combined_residual_j, 0.0);
        close(
            audit.thermal.external_heat_j,
            source.inertia().mass() * CP * RATE * 4.0 * DT,
        );
        close(next.thermal().energy_j(), audit.thermal.external_heat_j);
        close(next.thermal().temperature_k(CellId(0)).unwrap(), 302.0);
        assert_eq!(
            next.configuration().inertia().mass(),
            source.inertia().mass()
        );
        assert!(audit.viscous_heat_j > 0.0);
        let omega = next
            .principal_frame()
            .properties()
            .angular_velocity_body_checked(next.state().angular_momentum_body())
            .unwrap();
        close(
            audit.load_after.drag_force_world_n.x,
            -6.0 * std::f64::consts::PI * mu * radius * (v.x - flow().x),
        );
        close(
            audit.load_after.rotary_torque_body_nm.z,
            -8.0 * std::f64::consts::PI * mu * radius.powi(3) * omega.z,
        );
        close(
            audit.load_after.buoyancy_force_world_n.z,
            rho * 4.0 * std::f64::consts::PI * radius.powi(3) / 3.0 * 0.1,
        );
        outcomes.push((v, spin, audit.viscous_heat_j));
    }
    for i in 1..outcomes.len() {
        assert_ne!(outcomes[i], outcomes[i - 1]);
    }
}

#[test]
fn no_expansion_is_the_original_native_step_and_reference_rebuild_restarts_exactly() {
    let reservoir = fluid(0.002, 0.001);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    for alpha in [0.0, 0.02] {
        let material = card(CP, alpha, 500.0);
        let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
        let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
        let body = body(
            &source,
            &material,
            &chart,
            RATE,
            controls(),
            gravity(),
            twist(),
        );
        let (whole, audit) = body
            .advance_creeping_sphere(&sphere, flow_controls(), 8.0 * DT, config(), || Ok(()))
            .unwrap();
        let (half, first) = body
            .advance_creeping_sphere(&sphere, flow_controls(), 4.0 * DT, config(), || Ok(()))
            .unwrap();
        let (split, second) = half
            .advance_creeping_sphere(&sphere, flow_controls(), 8.0 * DT, config(), || Ok(()))
            .unwrap();
        assert_eq!(whole.state(), split.state());
        assert_eq!(whole.configuration(), split.configuration());
        assert_eq!(whole.thermal().cells(), split.thermal().cells());
        close(
            audit.viscous_heat_j,
            first.viscous_heat_j + second.viscous_heat_j,
        );
        close(
            audit.parameter_work_j,
            first.parameter_work_j + second.parameter_work_j,
        );
        close(
            audit.radius_after_m,
            RADIUS * (1.0 + alpha * RATE * 8.0 * DT),
        );
        if alpha == 0.0 {
            let original = source.principal_body(1e-12, || false).unwrap();
            let mut state = body.state();
            for _ in 0..8 {
                state = sphere
                    .step(
                        &original,
                        state,
                        CreepingStepControls {
                            duration_s: DT,
                            maximum_relaxation_fraction: flow_controls()
                                .maximum_relaxation_fraction,
                            maximum_energy_defect_j: flow_controls().maximum_step_energy_defect_j,
                        },
                        || false,
                    )
                    .unwrap()
                    .dynamics
                    .state_after;
            }
            assert_eq!(whole.state(), state);
            assert_eq!(audit.parameter_work_j, 0.0);
        }
    }
}

#[derive(Clone)]
struct PrescribedHeat {
    watts: f64,
    time_s: f64,
    supplied_j: f64,
    refuse_after_s: Option<f64>,
}
impl DepositionObserver for PrescribedHeat {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        Ok(BTreeMap::from([(CellId(0), self.watts)]))
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        assert_eq!(self.time_s, start);
        self.time_s = thermal.time_s();
        self.supplied_j += self.watts * dt;
        if self.refuse_after_s.is_some_and(|time| self.time_s > time) {
            return Err(DepositionError::InvalidInput(
                "test heat source refused endpoint",
            ));
        }
        Ok(())
    }
}

#[test]
fn generated_heat_and_fluid_loaded_motion_share_substeps_and_atomic_publication() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let direct = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let coupled = body(
        &source,
        &material,
        &chart,
        0.0,
        controls(),
        gravity(),
        twist(),
    );
    let reservoir = fluid(0.002, 0.001);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let input = PrescribedHeat {
        watts: source.inertia().mass() * CP * RATE,
        time_s: 0.0,
        supplied_j: 0.0,
        refuse_after_s: None,
    };
    let (expected, _) = direct
        .advance_creeping_sphere(&sphere, flow_controls(), 4.0 * DT, config(), || Ok(()))
        .unwrap();
    let (actual, generated, audit) = coupled
        .advance_creeping_sphere_with(&input, &sphere, flow_controls(), 4.0 * DT, config(), || {
            Ok(())
        })
        .unwrap();
    assert_eq!(actual.state(), expected.state());
    close(actual.thermal().energy_j(), expected.thermal().energy_j());
    close(generated.time_s, actual.thermal().time_s());
    close(generated.supplied_j, audit.thermal.external_heat_j);
    close(audit.combined_residual_j, 0.0);
    let failure = PrescribedHeat {
        refuse_after_s: Some(DT),
        ..input.clone()
    };
    assert!(
        coupled
            .advance_creeping_sphere_with(
                &failure,
                &sphere,
                flow_controls(),
                4.0 * DT,
                config(),
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(input.time_s, 0.0);
    assert_eq!(failure.supplied_j, 0.0);
    assert_eq!(coupled.thermal().time_s(), 0.0);
    assert_eq!(coupled.state(), direct.state());
}

#[test]
fn updated_radius_rechecks_regime_and_late_failures_preserve_both_histories() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let initial = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let reservoir = fluid(0.002, 0.001);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let (expected, audit) = initial
        .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || Ok(()))
        .unwrap();
    let mut count = 0;
    initial
        .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || {
            count += 1;
            Ok(())
        })
        .unwrap();
    let mut visited = 0;
    let failure = initial.advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || {
        visited += 1;
        if visited == count {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        failure,
        Err(HeatedMotionError::Heating(
            fs_mbd_material::thermal::heating::HeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    let mut tight = controls();
    tight.transfer.maximum_parameter_work_j = 0.5 * audit.absolute_parameter_work_j;
    let limited = body(&source, &material, &chart, RATE, tight, gravity(), twist());
    assert!(
        limited
            .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(limited.thermal().time_s(), 0.0);
    assert_eq!(limited.state(), initial.state());
    let (again, _) = initial
        .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || Ok(()))
        .unwrap();
    assert_eq!(expected.state(), again.state());
    assert_eq!(initial.thermal().time_s(), 0.0);

    // Start inside the creeping envelope, then let the SOURCE expansion make
    // the same body leave it at its very first newly accepted thermal state.
    let fast_heat = body(
        &source,
        &material,
        &chart,
        100.0,
        controls(),
        Gravity::ZERO,
        Twist::new(Vec3::default(), Vec3::new(0.03, 0.0, 0.0)),
    );
    let bounded = CreepingSphere::new(
        &reservoir,
        RADIUS,
        Gravity::ZERO,
        fs_mbd::Vec3::ZERO,
        CreepingFlowLimits {
            maximum_reynolds: 0.0121,
            ..limits()
        },
    )
    .unwrap();
    let original = source.principal_body(1e-12, || false).unwrap();
    assert!(bounded.load(&original, fast_heat.state(), || false).is_ok());
    assert!(matches!(
        fast_heat.advance_creeping_sphere(&bounded, flow_controls(), DT, config(), || Ok(())),
        Err(HeatedMotionError::CreepingFlow(
            CreepingFlowError::Reynolds { .. }
        ))
    ));
    assert_eq!(fast_heat.thermal().energy_j(), 0.0);

    let dense = fluid(0.006, 0.001);
    let bounded = CreepingSphere::new(
        &dense,
        RADIUS,
        Gravity::ZERO,
        fs_mbd::Vec3::ZERO,
        CreepingFlowLimits {
            maximum_diffusion_ratio: 0.046,
            ..limits()
        },
    )
    .unwrap();
    let quiet_heat = body(
        &source,
        &material,
        &chart,
        100.0,
        controls(),
        Gravity::ZERO,
        Twist::zero(),
    );
    assert!(matches!(
        quiet_heat.advance_creeping_sphere(&bounded, flow_controls(), DT, config(), || Ok(())),
        Err(HeatedMotionError::CreepingFlow(
            CreepingFlowError::Unsteady { .. }
        ))
    ));

    let wrong_gravity =
        CreepingSphere::new(&reservoir, RADIUS, Gravity::ZERO, flow(), limits()).unwrap();
    assert!(
        initial
            .advance_creeping_sphere(&wrong_gravity, flow_controls(), DT, config(), || Ok(()))
            .is_err()
    );
    let wrong_radius =
        CreepingSphere::new(&reservoir, 1.1 * RADIUS, gravity(), flow(), limits()).unwrap();
    assert!(
        initial
            .advance_creeping_sphere(&wrong_radius, flow_controls(), DT, config(), || Ok(()))
            .is_err()
    );
    let wrong = reference_body(
        &material,
        VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.1)).unwrap(),
    );
    let wrong = body(
        &wrong,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    assert!(
        wrong
            .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || Ok(()))
            .is_err()
    );
}

#[test]
fn radius_rebuilds_keep_the_original_fluid_snapshot_and_refuse_missing_expansion_data() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            let count = self.0.get();
            self.0.set(count + 1);
            // A second query would change the mechanics; an immutable reservoir
            // must not sample this later state when only its sphere expands.
            FluidTransportState::try_new(0.002, 1000.0, 0.1, 0.001 * (1 + count) as f64)
        }
    }
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let initial = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let reservoir = Counted(std::cell::Cell::new(0));
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let (half, _) = initial
        .advance_creeping_sphere(&sphere, flow_controls(), 2.0 * DT, config(), || Ok(()))
        .unwrap();
    half.advance_creeping_sphere(&sphere, flow_controls(), 4.0 * DT, config(), || Ok(()))
        .unwrap();
    assert_eq!(reservoir.0.get(), 1);
    let short = card(CP, 0.02, 300.75);
    let source = reference_body(&short, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&short, 1e5, SolidLiquidPhase::Solid);
    let initial = body(
        &source,
        &short,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    assert!(
        initial
            .advance_creeping_sphere(&sphere, flow_controls(), 4.0 * DT, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.thermal().time_s(), 0.0);
    assert_eq!(initial.thermal().energy_j(), 0.0);
}
