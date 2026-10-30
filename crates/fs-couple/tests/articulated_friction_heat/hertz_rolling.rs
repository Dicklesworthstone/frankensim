//! Independent planar mechanics and finite-cell reference for the joint owner.
//! Original material/source owners use explicitly synthetic data throughout.

use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_couple::friction_heat::hertz::RigidHertzHeatSource;
use fs_couple::friction_heat::hertz_rolling::*;
use fs_couple::friction_heat::rolling::RigidRollingHeatBinding;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::tangential_friction::hertz_rolling::RigidHertzRollingControls;
use fs_mbd::tangential_friction::rolling::RigidRollingHistory;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::sliding_patch::contact::{HertzWrenchConfig, HertzWrenchMode};
use fs_tribo::{ApplicabilityRange, HertzSpherePlane, InputAuthority};
use std::cell::Cell;
use std::collections::BTreeMap;

const R: f64 = 0.1;
const D: f64 = 0.05;
const INDENT: f64 = 0.001;
const TAU: f64 = 0.001;
const DT: f64 = 0.01;
const SLIP: f64 = 0.2;
const ROLL: f64 = 0.05;
const CP: [f64; 2] = [0.002, 0.001];
const SIDE_G: [f64; 2] = [0.2, 0.1];
const IDS: [CellId; 2] = [CellId(11), CellId(37)];
const MATERIALS: [MaterialId; 2] = [MaterialId(17), MaterialId(29)];
const BODY_IDS: [u64; 2] = [101, 203];
const LINK: u64 = 73;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}

struct Sensible(f64);
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + h / self.0;
        if temperature.is_finite() && temperature > 0.0 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.0.recip()
    }
}
static SPHERE_CALORIC: Sensible = Sensible(CP[0]);
static SUPPORT_CALORIC: Sensible = Sensible(CP[1]);

#[derive(Clone)]
struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    source: SourceThermalCoulomb,
    rolling: RollingLossLaw,
    patch: RollingPatchReceipt,
    work: RollingWorkOwnership,
    tau: f64,
    slip: f64,
    roll: f64,
}
impl Fixture {
    fn new(tau: f64, slip: f64, roll: f64) -> Self {
        let source = super::law(0);
        let normal = HertzSpherePlane {
            effective_radius: R,
            reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
        }
        .response(source.interface().interface(), INDENT)
        .unwrap();
        let patch = RollingPatchReceipt::new(
            "joint thermal Hertz sphere patch",
            "original Hertz sphere/plane",
            "synthetic original normal geometry",
            InputAuthority::SyntheticFixture,
            normal.normal_force_n,
            core::f64::consts::PI * normal.contact_radius_m.powi(2),
            PatchCurvature::Principal {
                first_per_m: R.recip(),
                second_per_m: R.recip(),
            },
        )
        .unwrap();
        let mass = 4.0 * core::f64::consts::PI * 1000.0 * R.powi(3) / 3.0;
        let inertia = 0.4 * mass * R * R;
        let properties = [
            MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
            MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
        ];
        let bodies = [
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, R), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(mass * (slip + R * roll), 0.0, 0.0),
                Vec3::new(0.0, inertia * roll, 0.0),
            )
            .unwrap(),
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, -D), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap(),
        ];
        let rolling = RollingLossLaw::ViscoelasticSphere(
            ViscoelasticSphereRollingCard::new(
                VISCOELASTIC_SPHERE_LAW_ID,
                "synthetic common relaxation time; not calibrated",
                InputAuthority::SyntheticFixture,
                ViscoelasticSphereParameters {
                    sphere_surface: DeformingSphereSurface::SurfaceA,
                    radius_m: R,
                    young_modulus_pa: 1e6,
                    poisson_ratio: 0.25,
                    density_kg_m3: 1000.0,
                    shear_viscosity_pa_s: tau * 1e6 / 2.5,
                    bulk_viscosity_pa_s: tau * 1e6 / 1.5,
                },
                QuasistaticSphereLimits {
                    maximum_indentation_radius_ratio: 0.1,
                    maximum_speed_wave_ratio: 0.1,
                    maximum_relaxation_ratio: 0.1,
                },
                RollingLossApplicability::new(
                    ApplicabilityRange::new(280.0, 320.0).unwrap(),
                    ApplicabilityRange::new(0.0, 100.0).unwrap(),
                )
                .unwrap(),
            )
            .unwrap(),
        );
        Self {
            bodies,
            properties,
            source,
            rolling,
            patch,
            work: RollingWorkOwnership::new(
                "joint thermal Hertz sphere patch",
                "joint thermal interval",
                "original tangent rolling vector",
                RollingLossChannel::RollingHysteresis,
            )
            .unwrap(),
            tau,
            slip,
            roll,
        }
    }
    fn contact(&self) -> RigidHertzRollingHeatContact<'_> {
        RigidHertzRollingHeatContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: Vec3::new(0.0, 0.0, -R),
                arm_b_body: Vec3::new(0.0, 0.0, D),
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            source: &self.source,
            patch: &self.patch,
            interface: self.source.interface().interface(),
            rolling: &self.rolling,
            excitation_frequency_hz: 10.0,
            work: &self.work,
            nominal_pressure_pa: None,
        }
    }
    fn conductivity(&self) -> [f64; 2] {
        [
            SIDE_G[0] * R / self.patch.contact_area_m2(),
            SIDE_G[1] * D / self.patch.contact_area_m2(),
        ]
    }
}

fn bulk_binding() -> RigidRollingHeatBinding {
    RigidRollingHeatBinding {
        body_index: 0,
        cell_id: IDS[0],
        thermal_body_id: BODY_IDS[0],
        material_id: MATERIALS[0],
    }
}
fn thermal(
    f: &Fixture,
    temperatures: [f64; 2],
    duplicate_sphere: bool,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: BODY_IDS[i],
                material_id: MATERIALS[i],
                mass_kg: f.properties[i].mass(),
                specific_h_j_kg: CP[i] * (temperatures[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, if i == 0 { 0.0 } else { -2.0 * D }],
                    end_m: [0.0, 0.0, if i == 0 { 2.0 * R } else { 0.0 }],
                    volume_m3: f.properties[i].mass() / 1000.0,
                    surface_area_m2: 4.0 * core::f64::consts::PI * R * R,
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
    if duplicate_sphere {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(99),
                segment_id: 2,
                body_id: BODY_IDS[0],
                material_id: MATERIALS[0],
                mass_kg: 0.001,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [1.0, 0.0, 0.0],
                    end_m: [2.0, 0.0, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: 0.001,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: IDS[0],
            b: IDS[1],
            area_m2: f.patch.contact_area_m2(),
            conductance_w_k: series_conductance_w_k(
                f.patch.contact_area_m2(),
                [R, D],
                f.conductivity(),
            )
            .unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MATERIALS[0], &SPHERE_CALORIC as &dyn ThermalLaw),
            (MATERIALS[1], &SUPPORT_CALORIC as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 1,
            max_events: 4,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn route(f: &Fixture, t: &DepositionSimulation<'_>) -> FrictionHeatRoute {
    FrictionHeatRoute::from_link(t, LINK, f.conductivity()).unwrap()
}
fn controls() -> RigidHertzRollingHeatControls {
    RigidHertzRollingHeatControls {
        mechanical: RigidHertzRollingControls {
            mechanical: RigidHertzControls {
                maximum_bodies: 2,
                geometry_tolerance_m: 1e-10,
                velocity_tolerance_m_s: 1e-10,
                angular_velocity_tolerance_rad_s: 1e-9,
                linear_momentum_tolerance_ns: 1e-10,
                angular_momentum_tolerance_nms: 1e-10,
                energy_tolerance_j: 2e-10,
                contact: HertzWrenchConfig {
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
            maximum_rolling_moment_nm: 10.0,
            rolling_moment_tolerance_nm: 2e-10,
        },
        maximum_cells: 3,
        maximum_trials: 256,
        maximum_newton_iterations: 25,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.01,
        heat: FrictionHeatTolerances {
            cell_j: 1e-10,
            thermal_j: 1e-10,
            interface_j: 1e-10,
            total_j: 2e-10,
        },
        combined_energy_tolerance_j: 3e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 3,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
    }
}

fn source(
    f: &Fixture,
    t: &DepositionSimulation<'_>,
    face: f64,
    external: [f64; 2],
) -> Result<RigidHertzHeatSource<[f64; 3]>, &'static str> {
    Ok(RigidHertzHeatSource {
        conductivity_w_m_k: f.conductivity(),
        transport: StepCoefficients {
            sources_w: BTreeMap::from([(IDS[0], external[0]), (IDS[1], external[1])]),
            conductances_w_k: BTreeMap::new(),
            payload: [
                face,
                t.temperature_k(IDS[0]).unwrap(),
                t.temperature_k(IDS[1]).unwrap(),
            ],
        },
    })
}

// Independent one-plane force/torque equations. No production Schur matrix or
// native mechanical report is used in this reference. The original symmetric
// sphere has I=2mR²/5; the support's authored y moment is 0.13 kg m².
#[derive(Clone, Copy)]
struct Reference {
    face: f64,
    temperatures: [f64; 2],
    force: f64,
    moment: f64,
    surface_heat: f64,
    bulk_heat: f64,
}
fn reference(f: &Fixture, old: [f64; 2], external: [f64; 2], dt: f64) -> Reference {
    let mass = f.properties[0].mass();
    let ia = 0.4 * mass * R * R;
    let ib = 0.13;
    let a = mass.recip() + 1.0 / 7.0 + R * R / ia + D * D / ib;
    let b = -R / ia + D / ib;
    let c = ia.recip() + ib.recip();
    let damping = f.tau * f.patch.normal_load_n() * R;
    let capacities = [mass * CP[0], 7.0 * CP[1]];
    let at = |face: f64| {
        let force = if f.slip == 0.0 && f.roll == 0.0 {
            0.0
        } else {
            -(super::MU[0] - super::SLOPE[0] * (face - 300.0)) * f.patch.normal_load_n()
        };
        let omega_mid = (f.roll + 0.5 * dt * b * force) / (1.0 + 0.5 * dt * damping * c);
        let moment = -damping * omega_mid;
        let slip_mid = f.slip + 0.5 * dt * (a * force + b * moment);
        let surface_heat = -dt * force * slip_mid;
        let bulk_heat = dt * damping * omega_mid * omega_mid;
        let temperatures = std::array::from_fn(|i| {
            (capacities[i] * old[i]
                + dt * SIDE_G[i] * face
                + dt * external[i]
                + if i == 0 { bulk_heat } else { 0.0 })
                / (capacities[i] + dt * SIDE_G[i])
        });
        let residual = dt
            * (0..2)
                .map(|i| SIDE_G[i] * (face - temperatures[i]))
                .sum::<f64>()
            - surface_heat;
        (
            residual,
            Reference {
                face,
                temperatures,
                force,
                moment,
                surface_heat,
                bulk_heat,
            },
        )
    };
    let (mut lo, mut hi) = (280.0, 320.0);
    assert!(at(lo).0 < 0.0 && at(hi).0 > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (lo + hi);
        if at(middle).0 > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    at(0.5 * (lo + hi)).1
}

fn energy_change(before: &DepositionSimulation<'_>, after: &DepositionSimulation<'_>) -> f64 {
    before
        .cells()
        .iter()
        .map(|(id, cell)| cell.mass_kg * (after.cells()[id].specific_h_j_kg - cell.specific_h_j_kg))
        .sum()
}
fn pair_energy(f: &Fixture, states: &[RigidBodyState; 2]) -> f64 {
    (0..2)
        .map(|i| states[i].kinetic_energy(f.properties[i]).unwrap())
        .sum()
}

#[test]
fn joint_surface_and_bulk_heat_match_independent_planar_thermal_root_and_recoil() {
    let f = Fixture::new(TAU, SLIP, ROLL);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let history = RigidRollingHistory::new();
    let expected = reference(&f, [300.0; 2], [0.0; 2], DT);
    let result = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(result.mechanical.hertz.mode, HertzWrenchMode::Sliding);
    assert!(expected.surface_heat > 1e-3 && expected.bulk_heat > 1e-8);
    near(result.mechanical.hertz.force_n[0], expected.force, 2e-9);
    near(
        result.mechanical.rolling_angular_impulse_world_nms.y,
        DT * expected.moment,
        2e-12,
    );
    near(
        result.mechanical.surface_heat_j,
        expected.surface_heat,
        5e-11,
    );
    near(result.mechanical.bulk_heat_j, expected.bulk_heat, 5e-12);
    near(result.source[0], expected.face, 1e-7);
    for i in 0..2 {
        near(
            result.thermal.temperature_k(IDS[i]).unwrap(),
            expected.temperatures[i],
            1e-7,
        );
    }
    let mass = f.properties[0].mass();
    let ia = 0.4 * mass * R * R;
    near(
        result.bodies_after[0].linear_momentum_world().x,
        mass * (SLIP + R * ROLL) + DT * expected.force,
        2e-11,
    );
    near(
        result.bodies_after[1].linear_momentum_world().x,
        -DT * expected.force,
        2e-11,
    );
    near(
        result.bodies_after[0].angular_momentum_body().y,
        ia * ROLL + DT * (-R * expected.force + expected.moment),
        2e-12,
    );
    near(
        result.bodies_after[1].angular_momentum_body().y,
        DT * (-D * expected.force - expected.moment),
        2e-12,
    );
    near(
        pair_energy(&f, &result.bodies_after) - pair_energy(&f, &f.bodies)
            + energy_change(&t, &result.thermal),
        0.0,
        3e-10,
    );
    near(
        result.rolling_history.state().cumulative_heat_j(),
        expected.bulk_heat,
        5e-12,
    );
    assert_eq!(result.rolling_history.state().accepted_steps(), 1);
    near(result.external_heat_j, 0.0, 1e-15);
    assert_eq!(history, RigidRollingHistory::new());
    assert_eq!(t.time_s(), 0.0);
}

#[test]
fn actual_interface_and_sphere_temperatures_keep_distinct_sources_and_external_heat() {
    // Initially zero rolling must remain coupled: shear creates angular motion
    // and bulk loss. The original source still belongs to the hot interface.
    let f = Fixture::new(TAU, SLIP, 0.0);
    let initial = [299.0, 301.0];
    let external = [0.15, -0.05];
    let t = thermal(&f, initial, false);
    let route = route(&f, &t);
    let expected = reference(&f, initial, external, DT);
    let result = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &RigidRollingHistory::new(),
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, external),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(result.mechanical.bulk_heat_j > 0.0);
    assert!((expected.face - expected.temperatures[0]).abs() > 1.0);
    near(result.sphere_temperature_before_k, initial[0], 1e-12);
    near(
        result.sphere_temperature_after_k,
        expected.temperatures[0],
        1e-7,
    );
    near(
        result.mechanical.rolling_temperature_k,
        expected.temperatures[0],
        1e-7,
    );
    for (actual, expected) in result.source.iter().zip([
        expected.face,
        expected.temperatures[0],
        expected.temperatures[1],
    ]) {
        near(*actual, expected, 1e-7);
    }
    let a = (f.patch.contact_area_m2() / core::f64::consts::PI).sqrt();
    let exact_source = f
        .source
        .hertz_at(result.source[0], None, f.patch.normal_load_n(), a)
        .unwrap();
    assert_eq!(result.friction_source, exact_source);
    assert_eq!(result.friction_source.owner_identity(), f.source.identity());
    near(result.mechanical.hertz.force_n[0], expected.force, 2e-9);
    near(
        result.mechanical.surface_heat_j,
        expected.surface_heat,
        5e-11,
    );
    near(result.mechanical.bulk_heat_j, expected.bulk_heat, 5e-12);
    near(
        result.external_heat_j,
        DT * external.iter().sum::<f64>(),
        1e-14,
    );
    near(
        pair_energy(&f, &result.bodies_after) - pair_energy(&f, &f.bodies)
            + energy_change(&t, &result.thermal)
            - DT * external.iter().sum::<f64>(),
        0.0,
        3e-10,
    );
    let restored = f
        .rolling
        .restore_checkpoint(
            &f.patch,
            f.source.interface().interface(),
            &f.work,
            result.rolling_history.checkpoint().unwrap(),
        )
        .unwrap();
    assert_eq!(&restored, result.rolling_history.state());
}

#[test]
fn stationary_and_zero_viscosity_limits_keep_original_conduction_and_heat_channels() {
    for f in [Fixture::new(TAU, 0.0, 0.0), Fixture::new(0.0, SLIP, ROLL)] {
        let initial = [295.0, 305.0];
        let t = thermal(&f, initial, false);
        let route = route(&f, &t);
        let expected = reference(&f, initial, [0.0; 2], DT);
        let result = propose_rigid_hertz_rolling_heat(
            &f.bodies,
            &f.properties,
            &RigidRollingHistory::new(),
            f.contact(),
            &t,
            &route,
            bulk_binding(),
            t.time_s(),
            DT,
            controls(),
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(result.mechanical.bulk_heat_j, 0.0);
        near(
            result.mechanical.surface_heat_j,
            expected.surface_heat,
            5e-11,
        );
        for i in 0..2 {
            near(
                result.thermal.temperature_k(IDS[i]).unwrap(),
                expected.temperatures[i],
                1e-7,
            );
        }
        assert!(result.thermal.temperature_k(IDS[0]).unwrap() > initial[0]);
        assert!(result.thermal.temperature_k(IDS[1]).unwrap() < initial[1]);
        if f.slip == 0.0 {
            assert_eq!(result.mechanical.surface_heat_j, 0.0);
            assert_eq!(result.bodies_after, f.bodies);
        } else {
            assert!(result.mechanical.surface_heat_j > 0.0);
        }
        near(
            pair_energy(&f, &result.bodies_after) - pair_energy(&f, &f.bodies)
                + energy_change(&t, &result.thermal),
            0.0,
            3e-10,
        );
    }
}

#[test]
fn source_geometry_whole_sphere_inventory_clock_and_effort_limits_refuse_atomically() {
    let f = Fixture::new(TAU, SLIP, ROLL);
    let history = RigidRollingHistory::new();
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let saved = t.cells().clone();
    for bad in [
        RigidRollingHeatBinding {
            body_index: 1,
            ..bulk_binding()
        },
        RigidRollingHeatBinding {
            thermal_body_id: 999,
            ..bulk_binding()
        },
        RigidRollingHeatBinding {
            material_id: MATERIALS[1],
            ..bulk_binding()
        },
        RigidRollingHeatBinding {
            cell_id: IDS[1],
            ..bulk_binding()
        },
    ] {
        assert!(
            propose_rigid_hertz_rolling_heat(
                &f.bodies,
                &f.properties,
                &history,
                f.contact(),
                &t,
                &route,
                bad,
                0.0,
                DT,
                controls(),
                &[0.0; 2],
                300.0,
                |state, face| source(&f, state, face, [0.0; 2]),
                &CancelGate::new(),
            )
            .is_err()
        );
    }
    let mut wrong_geometry = f.contact();
    wrong_geometry.binding.arm_a_body.x = 0.01;
    let mut wrong_pressure = f.contact();
    wrong_pressure.nominal_pressure_pa = Some(1.0);
    for contact in [wrong_geometry, wrong_pressure] {
        assert!(
            propose_rigid_hertz_rolling_heat(
                &f.bodies,
                &f.properties,
                &history,
                contact,
                &t,
                &route,
                bulk_binding(),
                0.0,
                DT,
                controls(),
                &[0.0; 2],
                300.0,
                |state, face| source(&f, state, face, [0.0; 2]),
                &CancelGate::new(),
            )
            .is_err()
        );
    }
    for (accepted, guess) in [(1.0, 300.0), (0.0, 321.0)] {
        assert!(
            propose_rigid_hertz_rolling_heat(
                &f.bodies,
                &f.properties,
                &history,
                f.contact(),
                &t,
                &route,
                bulk_binding(),
                accepted,
                DT,
                controls(),
                &[0.0; 2],
                guess,
                |state, face| source(&f, state, face, [0.0; 2]),
                &CancelGate::new(),
            )
            .is_err()
        );
    }
    let mut small_effort = controls();
    small_effort.mechanical.maximum_rolling_moment_nm = 1e-8;
    assert!(
        propose_rigid_hertz_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            &t,
            &route,
            bulk_binding(),
            0.0,
            DT,
            small_effort,
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .is_err()
    );
    // A finite midpoint sticking reaction would reverse this slow initial
    // slip. The joint owner must refuse the unresolved event, not invent an
    // endpoint stop or turn the resulting numerical loss into either heat.
    let slow = Fixture::new(TAU, 0.001, 0.0);
    assert!(
        propose_rigid_hertz_rolling_heat(
            &slow.bodies,
            &slow.properties,
            &history,
            slow.contact(),
            &t,
            &route,
            bulk_binding(),
            0.0,
            DT,
            controls(),
            &[0.0; 2],
            300.0,
            |state, face| source(&slow, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .is_err()
    );
    let outside_bulk = thermal(&f, [321.0, 300.0], false);
    let outside_route = self::route(&f, &outside_bulk);
    assert!(
        propose_rigid_hertz_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            &outside_bulk,
            &outside_route,
            bulk_binding(),
            0.0,
            DT,
            controls(),
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .is_err()
    );
    let duplicate = thermal(&f, [300.0; 2], true);
    let duplicate_route = self::route(&f, &duplicate);
    assert!(
        propose_rigid_hertz_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            &duplicate,
            &duplicate_route,
            bulk_binding(),
            0.0,
            DT,
            controls(),
            &[0.0; 3],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(history, RigidRollingHistory::new());
    assert_eq!(t.cells(), &saved);
    assert_eq!(t.time_s(), 0.0);
}

#[test]
fn loaded_joint_history_late_source_failure_cancellation_and_budget_have_exact_retry() {
    let mut f = Fixture::new(TAU, SLIP, ROLL);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let first = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &RigidRollingHistory::new(),
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    f.bodies = first.bodies_after;
    let history = first.rolling_history;
    let t = first.thermal;
    let saved_history = history.clone();
    let saved_bodies = f.bodies;
    let saved_cells = t.cells().clone();
    let saved_source = f.source.clone();
    let calls = Cell::new(0usize);
    let reference = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| {
            calls.set(calls.get() + 1);
            source(&f, state, face, [0.0; 2])
        },
        &CancelGate::new(),
    )
    .unwrap();
    let final_call = calls.get();
    assert!(final_call > 2);
    assert!(history.state().cumulative_heat_j() > 0.0);
    assert_eq!(reference.rolling_history.state().accepted_steps(), 2);
    let calls = Cell::new(0usize);
    let failed = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                Err("late joint thermal source")
            } else {
                source(&f, state, face, [0.0; 2])
            }
        },
        &CancelGate::new(),
    );
    assert!(failed.is_err());
    assert_eq!(calls.get(), final_call);
    let calls = Cell::new(0usize);
    let gate = CancelGate::new();
    let cancelled = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                gate.request();
            }
            source(&f, state, face, [0.0; 2])
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), final_call);
    assert!(
        propose_rigid_hertz_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            &t,
            &route,
            bulk_binding(),
            t.time_s(),
            DT,
            RigidHertzRollingHeatControls {
                maximum_trials: reference.trials - 1,
                ..controls()
            },
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .is_err()
    );
    let retry = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.mechanical, reference.mechanical);
    assert_eq!(retry.rolling_history, reference.rolling_history);
    assert_eq!(retry.friction_source, reference.friction_source);
    assert_eq!(retry.source, reference.source);
    assert_eq!(retry.thermal.cells(), reference.thermal.cells());
    assert_eq!(retry.thermal.time_s(), 2.0 * DT);
    assert_eq!(history, saved_history);
    assert_eq!(f.bodies, saved_bodies);
    assert_eq!(t.cells(), &saved_cells);
    assert_eq!(f.source, saved_source);
}
