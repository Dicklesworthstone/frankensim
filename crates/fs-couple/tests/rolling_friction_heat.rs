//! Independent finite-body and finite-caloric references for bulk rolling heat.
//! Synthetic material parameters exercise the original owners, not calibration.

use fs_couple::friction_heat::rolling::*;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::rolling::{
    RigidRollingContact, RigidRollingControls, RigidRollingHistory,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::{
    ApplicabilityRange, FrictionLaw, HertzSpherePlane, InputAuthority, InterfaceMedium,
    InterfaceSystemRef,
};
use std::cell::Cell;
use std::collections::BTreeMap;

const R: f64 = 0.1;
const D: f64 = 0.05;
const TAU: f64 = 0.01;
const RATE: [f64; 2] = [0.03, -0.02];
const DT: f64 = 0.2;
const G: f64 = 0.05;
const CP: [f64; 2] = [0.01, 0.02];
const IDS: [CellId; 2] = [CellId(11), CellId(37)];
const MATERIALS: [MaterialId; 2] = [MaterialId(17), MaterialId(29)];
const THERMAL_BODIES: [u64; 2] = [101, 203];
const LINK: u64 = 73;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn z_cross(v: Vec3) -> Vec3 {
    Vec3::new(-v.y, v.x, 0.0)
}

struct Sensible(f64);
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / self.0;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.0
    }
}
static SPHERE_CALORIC: Sensible = Sensible(CP[0]);
static SUPPORT_CALORIC: Sensible = Sensible(CP[1]);

struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let sensible = if h < 0.0 { h } else { (h - 0.001).max(0.0) };
        let t = 300.0 + sensible / CP[0];
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP[0]
    }
}
static LATENT: Latent = Latent;

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: RollingLossLaw,
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    friction: FrictionLaw,
    work: RollingWorkOwnership,
}
impl Fixture {
    fn contact(&self, temperature_k: f64) -> RigidRollingContact<'_> {
        RigidRollingContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: Vec3::new(0.0, 0.0, -R),
                arm_b_body: Vec3::new(0.0, 0.0, D),
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            patch: &self.patch,
            interface: &self.interface,
            law: &self.law,
            static_friction: &self.friction,
            temperature_k,
            excitation_frequency_hz: 10.0,
            work: &self.work,
        }
    }
}
fn law(tau: f64, surface: DeformingSphereSurface, source: &str) -> RollingLossLaw {
    let young = 1e6;
    let nu = 0.25;
    RollingLossLaw::ViscoelasticSphere(
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            source,
            InputAuthority::SyntheticFixture,
            ViscoelasticSphereParameters {
                sphere_surface: surface,
                radius_m: R,
                young_modulus_pa: young,
                poisson_ratio: nu,
                density_kg_m3: 1000.0,
                shear_viscosity_pa_s: tau * young / (2.0 * (1.0 + nu)),
                bulk_viscosity_pa_s: tau * young / (3.0 * (1.0 - 2.0 * nu)),
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
    )
}
fn fixture(tau: f64, rate: [f64; 2], indentation: f64) -> Fixture {
    let interface = InterfaceSystemRef::new(
        "synthetic whole sphere/finite support",
        "original bulk rolling history",
        "independent synthetic interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let normal = HertzSpherePlane {
        effective_radius: R,
        reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
    }
    .response(&interface, indentation)
    .unwrap();
    let patch = RollingPatchReceipt::new(
        "original bulk rolling patch",
        "original Hertz",
        "synthetic normal source",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        core::f64::consts::PI * normal.contact_radius_m.powi(2),
        PatchCurvature::Principal {
            first_per_m: 1.0 / R,
            second_per_m: 1.0 / R,
        },
    )
    .unwrap();
    let mass = 4.0 * core::f64::consts::PI * 1000.0 * R.powi(3) / 3.0;
    let inertia = 0.4 * mass * R * R;
    let properties = [
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
        MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
    ];
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.47).unwrap();
    let omega = rotation.rotate_body_to_world(Vec3::new(rate[0], rate[1], 0.0));
    Fixture {
        bodies: [
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, R), UnitQuaternion::IDENTITY).unwrap(),
                z_cross(omega).scale(-mass * R),
                omega.scale(inertia),
            )
            .unwrap(),
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, -D), rotation).unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap(),
        ],
        properties,
        law: law(
            tau,
            DeformingSphereSurface::SurfaceA,
            "synthetic original viscosity",
        ),
        patch,
        interface,
        friction: FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: 0.2,
        },
        work: RollingWorkOwnership::new(
            "original bulk rolling patch",
            "shared interval",
            "native tangent rolling vector",
            RollingLossChannel::RollingHysteresis,
        )
        .unwrap(),
    }
}
fn thermal(
    f: &Fixture,
    temperatures: [f64; 2],
    conductance: f64,
    sphere_law: &'static dyn ThermalLaw,
) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: THERMAL_BODIES[i],
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
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: IDS[0],
            b: IDS[1],
            area_m2: 0.001,
            conductance_w_k: conductance,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MATERIALS[0], sphere_law), (MATERIALS[1], &SUPPORT_CALORIC)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
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
fn binding() -> RigidRollingHeatBinding {
    RigidRollingHeatBinding {
        body_index: 0,
        cell_id: IDS[0],
        thermal_body_id: THERMAL_BODIES[0],
        material_id: MATERIALS[0],
    }
}
fn controls() -> RigidRollingHeatControls {
    RigidRollingHeatControls {
        mechanical: RigidRollingControls {
            maximum_bodies: 3,
            maximum_substeps: 4096,
            maximum_step_s: DT / 4.0,
            geometry_tolerance_m: 1e-12,
            velocity_tolerance_m_s: 1e-11,
            angular_velocity_tolerance_rad_s: 1e-10,
            linear_momentum_tolerance_ns: 1e-11,
            angular_momentum_tolerance_nms: 1e-11,
            energy_tolerance_j: 1e-12,
        },
        maximum_cells: 2,
        maximum_trials: 256,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-6,
        cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10,
        combined_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 2,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
    }
}
fn source(t: &DepositionSimulation<'_>) -> Result<StepCoefficients<[f64; 2]>, &'static str> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: IDS.map(|id| t.temperature_k(id).unwrap()),
    })
}
fn propose<'m>(
    f: &Fixture,
    history: &RigidRollingHistory,
    thermal: &DepositionSimulation<'m>,
) -> RigidRollingHeatCandidate<'m, [f64; 2]> {
    propose_rigid_rolling_heat(
        &f.bodies,
        &f.properties,
        history,
        f.contact(thermal.temperature_k(IDS[0]).unwrap()),
        thermal,
        thermal.time_s(),
        binding(),
        DT,
        controls(),
        &[0.0; 2],
        source,
        &CancelGate::new(),
    )
    .unwrap()
}

// Independently eliminate workless point force in each support principal axis.
// aF+b(z cross M)=0, so k=1/Is+1/Ip-b²/a. Proportional viscosity gives A=tau.
fn mechanical_reference(f: &Fixture, dt: f64, steps: usize) -> ([RigidBodyState; 2], f64) {
    let is = f.properties[0].principal_inertia_body().x;
    let rotation = f.bodies[1].pose().orientation();
    let mut force = Vec3::ZERO;
    let mut couple = Vec3::ZERO;
    let mut heat = 0.0;
    for (i, ip) in [0.08, 0.13].into_iter().enumerate() {
        let a =
            1.0 / f.properties[0].mass() + 1.0 / f.properties[1].mass() + R * R / is + D * D / ip;
        let b = R / is - D / ip;
        let k = 1.0 / is + 1.0 / ip - b * b / a;
        let alpha = 0.5 * dt / steps as f64 * TAU * R * f.patch.normal_load_n() * k;
        let after = RATE[i] * ((1.0 - alpha) / (1.0 + alpha)).powi(steps as i32);
        let j = (after - RATE[i]) / k;
        let direction = rotation.rotate_body_to_world(if i == 0 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        });
        couple = couple.add(direction.scale(j));
        force = force.add(z_cross(direction).scale(-b * j / a));
        heat += 0.5 * (RATE[i] * RATE[i] - after * after) / k;
    }
    (
        [
            RigidBodyState::new(
                f.bodies[0].pose(),
                f.bodies[0].linear_momentum_world().add(force),
                f.bodies[0]
                    .angular_momentum_body()
                    .add(z_cross(force).scale(-R))
                    .add(couple),
            )
            .unwrap(),
            RigidBodyState::new(
                f.bodies[1].pose(),
                f.bodies[1].linear_momentum_world().sub(force),
                rotation.rotate_world_to_body(z_cross(force).scale(-D).sub(couple)),
            )
            .unwrap(),
        ],
        heat,
    )
}
fn caloric_reference(f: &Fixture, old: [f64; 2], heat: f64, g: f64) -> [f64; 2] {
    let c = [
        f.properties[0].mass() * CP[0],
        f.properties[1].mass() * CP[1],
    ];
    let gdt = DT * g;
    let determinant = c[0] * c[1] + gdt * (c[0] + c[1]);
    let difference = old[0] - old[1];
    [
        old[0] + ((c[1] + gdt) * heat - gdt * c[1] * difference) / determinant,
        old[1] + (gdt * heat + gdt * c[0] * difference) / determinant,
    ]
}
fn energy_change(a: &DepositionSimulation<'_>, b: &DepositionSimulation<'_>) -> f64 {
    a.cells()
        .iter()
        .map(|(id, c)| c.mass_kg * (b.cells()[id].specific_h_j_kg - c.specific_h_j_kg))
        .sum()
}

#[test]
fn native_both_body_recoil_and_single_bulk_heat_match_independent_caloric_conduction() {
    let f = fixture(TAU, RATE, 0.001);
    let thermal = thermal(&f, [300.0; 2], G, &SPHERE_CALORIC);
    let history = RigidRollingHistory::new();
    let candidate = propose(&f, &history, &thermal);
    let report = candidate.report();
    let (expected, heat) = mechanical_reference(&f, DT, 4);
    assert_eq!(report.mechanical.substeps, 4);
    assert!(report.mechanical.rolling_mobility[0][1].abs() > 1.0);
    for i in 0..2 {
        close_vec(
            candidate.bodies()[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            3e-16,
        );
        close_vec(
            candidate.bodies()[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            3e-17,
        );
    }
    close(report.mechanical.dissipated_heat_j, heat, 3e-18);
    let expected_t = caloric_reference(&f, [300.0; 2], heat, G);
    for i in 0..2 {
        close(
            candidate.thermal().temperature_k(IDS[i]).unwrap(),
            expected_t[i],
            2e-10,
        );
        assert_eq!(
            report.source[i],
            candidate.thermal().temperature_k(IDS[i]).unwrap()
        );
    }
    let plane_heat = f.properties[1].mass() * candidate.thermal().cells()[&IDS[1]].specific_h_j_kg;
    close(plane_heat, DT * G * (expected_t[0] - expected_t[1]), 2e-12);
    assert!(plane_heat > 0.0);
    assert!(plane_heat < heat);
    close(energy_change(&thermal, candidate.thermal()), heat, 2e-12);
    close(
        report.mechanical.kinetic_after_j - report.mechanical.kinetic_before_j
            + energy_change(&thermal, candidate.thermal()),
        0.0,
        2e-12,
    );
    assert_eq!(report.external_heat_j, 0.0);
    assert_eq!(report.binding, binding());
    assert_eq!(candidate.body_indices(), [0, 1]);
    assert_eq!(candidate.history().state().accepted_steps(), 4);
    close(candidate.history().state().cumulative_heat_j(), heat, 3e-18);
    assert_eq!(
        report.mechanical.final_step.generalized_work.ownership,
        f.work
    );
    assert_eq!(
        candidate.thermal().contact_patch(LINK),
        thermal.contact_patch(LINK)
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(history, RigidRollingHistory::new());
}

#[test]
fn stationary_zero_viscosity_and_released_patch_keep_ordinary_conduction() {
    for (tau, rate, indentation) in [(TAU, [0.0; 2], 0.001), (0.0, RATE, 0.001), (TAU, RATE, 0.0)] {
        let f = fixture(tau, rate, indentation);
        let initial = [310.0, 290.0];
        let thermal = thermal(&f, initial, G, &SPHERE_CALORIC);
        let candidate = propose(&f, &RigidRollingHistory::new(), &thermal);
        let expected = caloric_reference(&f, initial, 0.0, G);
        assert_eq!(*candidate.bodies(), f.bodies);
        assert_eq!(candidate.report().mechanical.dissipated_heat_j, 0.0);
        for i in 0..2 {
            close(
                candidate.thermal().temperature_k(IDS[i]).unwrap(),
                expected[i],
                3e-10,
            );
        }
        assert!(expected[0] < initial[0]);
        assert!(expected[1] > initial[1]);
        close(energy_change(&thermal, candidate.thermal()), 0.0, 3e-11);
        assert_eq!(
            candidate.thermal().contact_patch(LINK),
            thermal.contact_patch(LINK)
        );
    }
}

#[test]
fn original_latent_enthalpy_stores_bulk_work_without_a_fabricated_heat_capacity() {
    let f = fixture(TAU, RATE, 0.001);
    let thermal = thermal(&f, [300.0; 2], 0.0, &LATENT);
    let candidate = propose(&f, &RigidRollingHistory::new(), &thermal);
    let (_, heat) = mechanical_reference(&f, DT, 4);
    assert!(heat > 0.0);
    assert_eq!(candidate.report().sphere_temperature_after_k, 300.0);
    assert_eq!(candidate.report().sphere_temperature_before_k, 300.0);
    close(energy_change(&thermal, candidate.thermal()), heat, 2e-12);
    close(
        candidate.thermal().cells()[&IDS[0]].specific_h_j_kg,
        heat / f.properties[0].mass(),
        2e-13,
    );
    assert_eq!(candidate.thermal().cells()[&IDS[1]].specific_h_j_kg, 0.0);
    assert!(candidate.thermal().uses_material_law(MATERIALS[0], &LATENT));
    assert_eq!(candidate.history().state().accepted_steps(), 4);
}

#[test]
fn other_heat_and_endpoint_conductance_keep_separate_physical_energy_accounts() {
    let f = fixture(TAU, RATE, 0.001);
    let thermal = thermal(&f, [300.0; 2], G, &SPHERE_CALORIC);
    let history = RigidRollingHistory::new();
    let powers = [0.012, -0.003];
    let conductance = 2.0 * G;
    let candidate = propose_rigid_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(300.0),
        &thermal,
        0.0,
        binding(),
        DT,
        controls(),
        &[0.0; 2],
        |t| -> Result<_, &'static str> {
            Ok(StepCoefficients {
                sources_w: BTreeMap::from([(IDS[0], powers[0]), (IDS[1], powers[1])]),
                conductances_w_k: BTreeMap::from([(LINK, conductance)]),
                payload: IDS.map(|id| t.temperature_k(id).unwrap()),
            })
        },
        &CancelGate::new(),
    )
    .unwrap();
    let report = candidate.report();
    let external = DT * (powers[0] + powers[1]);
    let c = [
        f.properties[0].mass() * CP[0],
        f.properties[1].mass() * CP[1],
    ];
    let q = [
        report.mechanical.dissipated_heat_j + DT * powers[0],
        DT * powers[1],
    ];
    let gd = DT * conductance;
    let determinant = c[0] * c[1] + gd * (c[0] + c[1]);
    let expected = [
        300.0 + ((c[1] + gd) * q[0] + gd * q[1]) / determinant,
        300.0 + (gd * q[0] + (c[0] + gd) * q[1]) / determinant,
    ];
    for i in 0..2 {
        close(
            candidate.thermal().temperature_k(IDS[i]).unwrap(),
            expected[i],
            2e-10,
        );
    }
    close(report.external_heat_j, external, 1e-17);
    close(
        energy_change(&thermal, candidate.thermal()),
        external + report.mechanical.dissipated_heat_j,
        2e-12,
    );
    close(
        report.thermal.external_heat_j,
        external + report.mechanical.dissipated_heat_j,
        1e-17,
    );
    close(report.combined_energy_residual_j, 0.0, 2e-12);
    // Endpoint transport updates are not published as authored network inputs.
    assert_eq!(
        candidate.thermal().contact_patch(LINK),
        thermal.contact_patch(LINK)
    );
    for id in IDS {
        assert_eq!(
            candidate.thermal().cells()[&id].boundary,
            thermal.cells()[&id].boundary
        );
    }
}

#[test]
fn surface_b_and_original_nonparticipant_indices_keep_heat_in_the_deforming_sphere() {
    let f = fixture(TAU, RATE, 0.001);
    let thermal = thermal(&f, [300.0; 2], G, &SPHERE_CALORIC);
    let baseline = propose(&f, &RigidRollingHistory::new(), &thermal);
    let spectator = RigidBodyState::new(
        Pose::new(Vec3::new(7.0, -11.0, 3.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1e18, -2e18, 3e18),
        Vec3::new(1e17, 2e17, -3e17),
    )
    .unwrap();
    let scene = [spectator, f.bodies[1], f.bodies[0]];
    let properties = [
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(),
        f.properties[1],
        f.properties[0],
    ];
    let swapped_law = law(
        TAU,
        DeformingSphereSurface::SurfaceB,
        "synthetic original B viscosity",
    );
    let candidate = propose_rigid_rolling_heat(
        &scene,
        &properties,
        &RigidRollingHistory::new(),
        RigidRollingContact {
            binding: RigidContactBinding {
                body_a: 1,
                body_b: 2,
                arm_a_body: Vec3::new(0.0, 0.0, D),
                arm_b_body: Vec3::new(0.0, 0.0, -R),
            },
            frame: TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap(),
            law: &swapped_law,
            ..f.contact(300.0)
        },
        &thermal,
        0.0,
        RigidRollingHeatBinding {
            body_index: 2,
            ..binding()
        },
        DT,
        controls(),
        &[0.0; 2],
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(candidate.body_indices(), [1, 2]);
    assert_eq!(candidate.report().mechanical.deforming_body, 2);
    assert_eq!(candidate.report().binding.cell_id, IDS[0]);
    close(
        candidate.report().mechanical.dissipated_heat_j,
        baseline.report().mechanical.dissipated_heat_j,
        3e-18,
    );
    for i in 0..2 {
        close_vec(
            candidate.bodies()[1 - i].linear_momentum_world(),
            baseline.bodies()[i].linear_momentum_world(),
            3e-16,
        );
        close_vec(
            candidate.bodies()[1 - i].angular_momentum_body(),
            baseline.bodies()[i].angular_momentum_body(),
            3e-17,
        );
        close(
            candidate.thermal().temperature_k(IDS[i]).unwrap(),
            baseline.thermal().temperature_k(IDS[i]).unwrap(),
            1e-12,
        );
    }
    assert_eq!(scene[0], spectator);
    assert_eq!([scene[2], scene[1]], f.bodies);
}

#[test]
fn source_temperature_and_original_bulk_cell_admission_refuse_without_mutation() {
    let f = fixture(TAU, RATE, 0.001);
    let thermal = thermal(&f, [300.0; 2], G, &SPHERE_CALORIC);
    let history = RigidRollingHistory::new();
    let saved_cells = thermal.cells().clone();
    for bad in [
        RigidRollingHeatBinding {
            body_index: 1,
            ..binding()
        },
        RigidRollingHeatBinding {
            cell_id: CellId(999),
            ..binding()
        },
        RigidRollingHeatBinding {
            thermal_body_id: 999,
            ..binding()
        },
        RigidRollingHeatBinding {
            material_id: MATERIALS[1],
            ..binding()
        },
        RigidRollingHeatBinding {
            cell_id: IDS[1],
            thermal_body_id: THERMAL_BODIES[1],
            material_id: MATERIALS[1],
            ..binding()
        },
    ] {
        let result = propose_rigid_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(300.0),
            &thermal,
            0.0,
            bad,
            DT,
            controls(),
            &[0.0; 2],
            source,
            &CancelGate::new(),
        );
        assert!(matches!(result, Err(RigidRollingHeatError::Invalid(_))));
    }
    // The bound cell already stores the full native sphere mass. A second
    // active cell with its body ID would duplicate that caloric inventory.
    let mut extra = saved_cells[&IDS[1]].clone();
    extra.body_id = THERMAL_BODIES[0];
    let duplicate_inventory = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(saved_cells[&IDS[0]].clone()),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(extra),
            },
        ],
        &[
            (MATERIALS[0], &SPHERE_CALORIC),
            (MATERIALS[1], &SUPPORT_CALORIC),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 100,
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
    .0;
    assert!(matches!(
        propose_rigid_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(300.0),
            &duplicate_inventory,
            0.0,
            binding(),
            DT,
            controls(),
            &[0.0; 2],
            source,
            &CancelGate::new(),
        ),
        Err(RigidRollingHeatError::Invalid(
            "whole-sphere thermal body must have exactly one active caloric cell"
        ))
    ));
    let stale_temperature = propose_rigid_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(301.0),
        &thermal,
        0.0,
        binding(),
        DT,
        controls(),
        &[0.0; 2],
        source,
        &CancelGate::new(),
    );
    assert!(matches!(
        stale_temperature,
        Err(RigidRollingHeatError::Invalid(_))
    ));
    // Initial T and mechanical state are admitted, but the source tries to heat
    // the actual finite sphere beyond the original 320 K material envelope.
    let too_hot = propose_rigid_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(300.0),
        &thermal,
        0.0,
        binding(),
        DT,
        controls(),
        &[0.0; 2],
        |t| -> Result<_, &'static str> {
            Ok(StepCoefficients {
                sources_w: BTreeMap::from([(IDS[0], 100.0)]),
                ..source(t)?
            })
        },
        &CancelGate::new(),
    );
    assert!(matches!(too_hot, Err(RigidRollingHeatError::Mechanical(_))));
    let bounded = RigidRollingHeatControls {
        maximum_cells: 1,
        ..controls()
    };
    assert!(matches!(
        propose_rigid_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(300.0),
            &thermal,
            0.0,
            binding(),
            DT,
            bounded,
            &[0.0; 2],
            source,
            &CancelGate::new(),
        ),
        Err(RigidRollingHeatError::Invalid(_))
    ));
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(history, RigidRollingHistory::new());
}

#[test]
fn loaded_history_final_failure_cancellation_and_budget_have_exact_retry() {
    let mut f = fixture(TAU, RATE, 0.001);
    let initial = thermal(&f, [300.0; 2], G, &SPHERE_CALORIC);
    let first = propose(&f, &RigidRollingHistory::new(), &initial);
    let (_, accepted_bodies, history, thermal, first_report) = first.into_parts();
    f.bodies = accepted_bodies;
    let saved_bodies = f.bodies;
    let saved_history = history.clone();
    let saved_cells = thermal.cells().clone();
    let calls = Cell::new(0usize);
    let run = |resolve: &dyn Fn(
        &DepositionSimulation<'static>,
    ) -> Result<StepCoefficients<[f64; 2]>, &'static str>,
               c,
               gate: &CancelGate| {
        propose_rigid_rolling_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(thermal.temperature_k(IDS[0]).unwrap()),
            &thermal,
            thermal.time_s(),
            binding(),
            DT,
            c,
            &[0.0; 2],
            resolve,
            gate,
        )
    };
    let baseline = run(
        &|t| {
            calls.set(calls.get() + 1);
            source(t)
        },
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(calls.get(), baseline.report().trials);
    let final_call = calls.get();
    assert!(final_call > 2);
    calls.set(0);
    let failed = run(
        &|t| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                Err("late original source")
            } else {
                source(t)
            }
        },
        controls(),
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidRollingHeatError::Source("late original source"))
    ));
    assert_eq!(calls.get(), final_call);
    let gate = CancelGate::new();
    calls.set(0);
    let cancelled = run(
        &|t| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                gate.request();
            }
            source(t)
        },
        controls(),
        &gate,
    );
    assert!(matches!(
        cancelled,
        Err(RigidRollingHeatError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(calls.get(), final_call);
    assert!(matches!(
        run(
            &source,
            RigidRollingHeatControls {
                maximum_trials: final_call - 1,
                ..controls()
            },
            &CancelGate::new()
        ),
        Err(RigidRollingHeatError::TrialBudget)
    ));
    let replacement = law(
        TAU,
        DeformingSphereSurface::SurfaceA,
        "different source must not reset history",
    );
    let changed = propose_rigid_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        RigidRollingContact {
            law: &replacement,
            ..f.contact(thermal.temperature_k(IDS[0]).unwrap())
        },
        &thermal,
        thermal.time_s(),
        binding(),
        DT,
        controls(),
        &[0.0; 2],
        source,
        &CancelGate::new(),
    );
    assert!(matches!(changed, Err(RigidRollingHeatError::Mechanical(_))));
    let retry = run(&source, controls(), &CancelGate::new()).unwrap();
    assert_eq!(retry.bodies(), baseline.bodies());
    assert_eq!(retry.history(), baseline.history());
    assert_eq!(retry.thermal().cells(), baseline.thermal().cells());
    assert_eq!(retry.report().source, baseline.report().source);
    assert_eq!(retry.report().mechanical, baseline.report().mechanical);
    close(
        retry.history().state().cumulative_heat_j(),
        first_report.mechanical.dissipated_heat_j + retry.report().mechanical.dissipated_heat_j,
        3e-18,
    );
    assert_eq!(retry.history().state().accepted_steps(), 8);
    assert_eq!(f.bodies, saved_bodies);
    assert_eq!(history, saved_history);
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), DT);
    assert_eq!(retry.thermal().time_s(), 2.0 * DT);
}
