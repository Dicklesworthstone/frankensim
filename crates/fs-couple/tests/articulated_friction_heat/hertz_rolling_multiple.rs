//! G1/G3/G4: one finite sphere between two finite supports, two original faces.
//! All material/source data are explicit synthetic fixtures. The reference
//! solves independent force/torque and three-capacity equations, not production
//! Schur matrices or mechanical heat reports.

use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_couple::friction_heat::hertz::multiple::{
    RigidHertzIslandHeatConductivity, RigidHertzIslandHeatSource,
};
use fs_couple::friction_heat::hertz_rolling::RigidHertzRollingHeatContact;
use fs_couple::friction_heat::hertz_rolling::multiple::*;
use fs_couple::friction_heat::rolling::RigidRollingHeatBinding;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::tangential_friction::hertz_rolling::RigidHertzRollingControls;
use fs_mbd::tangential_friction::hertz_rolling::multiple::RigidHertzRollingIslandControls;
use fs_mbd::tangential_friction::rolling::RigidRollingHistory;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::sliding_patch::contact::HertzWrenchConfig;
use fs_tribo::{ApplicabilityRange, HertzSpherePlane, InputAuthority};
use std::cell::Cell;
use std::collections::BTreeMap;

const R: f64 = 0.1;
const D: f64 = 0.05;
const DT: f64 = 0.01;
const TAU: f64 = 0.001;
const CP: [f64; 3] = [0.002, 0.001, 0.0015];
const G: [[f64; 2]; 2] = [[0.2, 0.1], [0.25, 0.15]];
const IDS: [CellId; 3] = [CellId(11), CellId(37), CellId(53)];
const MATS: [MaterialId; 3] = [MaterialId(17), MaterialId(29), MaterialId(43)];
const BODY: [u64; 3] = [101, 203, 307];
const LINKS: [u64; 2] = [73, 31];
const PATCHES: [&str; 2] = ["bottom sphere patch", "top sphere patch"];
const SIGNS: [f64; 2] = [-1.0, 1.0];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
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
        self.0.recip()
    }
}
static CALORIC: [Sensible; 3] = [Sensible(CP[0]), Sensible(CP[1]), Sensible(CP[2])];

#[derive(Clone)]
struct Fixture {
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
    sources: [SourceThermalCoulomb; 2],
    patches: [RollingPatchReceipt; 2],
    work: [RollingWorkOwnership; 2],
    rolling: RollingLossLaw,
    tau: f64,
    speed: f64,
    omega: f64,
}
impl Fixture {
    fn new(tau: f64, speed: f64, omega: f64) -> Self {
        let sources = [super::law(0), super::law(1)];
        let normal = HertzSpherePlane {
            effective_radius: R,
            reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
        }
        .response(sources[0].interface().interface(), 0.001)
        .unwrap();
        let patches = std::array::from_fn(|i| {
            RollingPatchReceipt::new(
                PATCHES[i],
                "original Hertz sphere/plane",
                "independent opposite sphere pole",
                InputAuthority::SyntheticFixture,
                normal.normal_force_n,
                core::f64::consts::PI * normal.contact_radius_m.powi(2),
                PatchCurvature::Principal {
                    first_per_m: R.recip(),
                    second_per_m: R.recip(),
                },
            )
            .unwrap()
        });
        let work = std::array::from_fn(|i| {
            RollingWorkOwnership::new(
                PATCHES[i],
                "same original thermal interval",
                "tangent rolling vector",
                RollingLossChannel::RollingHysteresis,
            )
            .unwrap()
        });
        let mass = 4.0 * core::f64::consts::PI * 1000.0 * R.powi(3) / 3.0;
        let inertia = 0.4 * mass * R * R;
        let properties = [
            MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
            MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
            MassProperties::new(5.0, Vec3::ZERO, Vec3::new(0.09, 0.17, 0.21)).unwrap(),
        ];
        let bodies = std::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(
                    Vec3::new(0.0, 0.0, if i == 0 { 0.0 } else { SIGNS[i - 1] * (R + D) }),
                    UnitQuaternion::IDENTITY,
                )
                .unwrap(),
                if i == 0 {
                    Vec3::new(mass * speed, 0.0, 0.0)
                } else {
                    Vec3::ZERO
                },
                if i == 0 {
                    Vec3::new(0.0, inertia * omega, 0.0)
                } else {
                    Vec3::ZERO
                },
            )
            .unwrap()
        });
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
            sources,
            patches,
            work,
            rolling,
            tau,
            speed,
            omega,
        }
    }
    fn contact(&self, i: usize) -> RigidHertzRollingHeatContact<'_> {
        RigidHertzRollingHeatContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: i + 1,
                arm_a_body: Vec3::new(0.0, 0.0, SIGNS[i] * R),
                arm_b_body: Vec3::new(0.0, 0.0, -SIGNS[i] * D),
            },
            frame: TangentFrame::new([0.0, 0.0, -SIGNS[i]], [1.0, 0.0, 0.0]).unwrap(),
            source: &self.sources[i],
            patch: &self.patches[i],
            interface: self.sources[i].interface().interface(),
            rolling: &self.rolling,
            excitation_frequency_hz: 10.0,
            work: &self.work[i],
            nominal_pressure_pa: None,
        }
    }
    fn conductivity(&self, i: usize) -> [f64; 2] {
        [
            G[i][0] * R / self.patches[i].contact_area_m2(),
            G[i][1] * D / self.patches[i].contact_area_m2(),
        ]
    }
    fn heat_contacts<'a>(
        &'a self,
        routes: &'a [FrictionHeatRoute; 2],
        histories: &'a [RigidRollingHistory; 2],
    ) -> [RigidHertzRollingIslandHeatContact<'a>; 2] {
        std::array::from_fn(|i| RigidHertzRollingIslandHeatContact {
            patch_id: PATCHES[i],
            contact: self.contact(i),
            history: &histories[i],
            route: &routes[i],
            bulk_binding: bulk_binding(),
            temperature_guess_k: 300.0,
        })
    }
}
fn bulk_binding() -> RigidRollingHeatBinding {
    RigidRollingHeatBinding {
        body_index: 0,
        cell_id: IDS[0],
        thermal_body_id: BODY[0],
        material_id: MATS[0],
    }
}
fn thermal(
    f: &Fixture,
    temperatures: [f64; 3],
    duplicate_sphere: bool,
) -> DepositionSimulation<'static> {
    let ends = [[-R, R], [-R - 2.0 * D, -R], [R, R + 2.0 * D]];
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: BODY[i],
                material_id: MATS[i],
                mass_kg: f.properties[i].mass(),
                specific_h_j_kg: CP[i] * (temperatures[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, ends[i][0]],
                    end_m: [0.0, 0.0, ends[i][1]],
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
        let mut other = match &events[0].event {
            Event::Birth(c) => c.clone(),
            _ => unreachable!(),
        };
        other.id = CellId(99);
        other.segment_id = 99;
        other.mass_kg = 0.001;
        other.geometry.start_m = [1.0, 0.0, 0.0];
        other.geometry.end_m = [2.0, 0.0, 0.0];
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(other),
        });
    }
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: IDS[0],
                b: IDS[i + 1],
                area_m2: f.patches[i].contact_area_m2(),
                conductance_w_k: series_conductance_w_k(
                    f.patches[i].contact_area_m2(),
                    [R, D],
                    f.conductivity(i),
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MATS[0], &CALORIC[0] as &dyn ThermalLaw),
            (MATS[1], &CALORIC[1]),
            (MATS[2], &CALORIC[2]),
        ],
        DepositionBudget {
            max_cells: 4,
            max_links: 2,
            max_events: 6,
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
fn routes(f: &Fixture, t: &DepositionSimulation<'_>) -> [FrictionHeatRoute; 2] {
    std::array::from_fn(|i| FrictionHeatRoute::from_link(t, LINKS[i], f.conductivity(i)).unwrap())
}
fn controls() -> RigidHertzRollingIslandHeatControls {
    RigidHertzRollingIslandHeatControls {
        mechanical: RigidHertzRollingIslandControls {
            mechanical: RigidHertzRollingControls {
                mechanical: RigidHertzControls {
                    maximum_bodies: 3,
                    geometry_tolerance_m: 1e-10,
                    velocity_tolerance_m_s: 1e-9,
                    angular_velocity_tolerance_rad_s: 1e-8,
                    linear_momentum_tolerance_ns: 1e-10,
                    angular_momentum_tolerance_nms: 1e-10,
                    energy_tolerance_j: 5e-10,
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
            maximum_contacts: 2,
            maximum_matrix_terms: 100,
            maximum_sweeps: 256,
            maximum_evaluations: 100_000,
        },
        maximum_cells: 4,
        maximum_unknowns: 6,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.01,
        heat: FrictionHeatTolerances {
            cell_j: 2e-10,
            thermal_j: 2e-10,
            interface_j: 2e-10,
            total_j: 5e-10,
        },
        combined_energy_tolerance_j: 7e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 5,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
    }
}
fn source(
    f: &Fixture,
    endpoint: &DepositionSimulation<'_>,
    links: &[u64],
    faces: &[f64],
    outside: [f64; 3],
) -> Result<RigidHertzIslandHeatSource<Vec<f64>>, &'static str> {
    Ok(RigidHertzIslandHeatSource {
        interfaces: links
            .iter()
            .map(|&id| {
                let i = LINKS.iter().position(|&x| x == id).unwrap();
                RigidHertzIslandHeatConductivity {
                    link_id: id,
                    conductivity_w_m_k: f.conductivity(i),
                }
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: IDS.into_iter().zip(outside).collect(),
            conductances_w_k: BTreeMap::new(),
            payload: faces
                .iter()
                .copied()
                .chain(IDS.map(|id| endpoint.temperature_k(id).unwrap()))
                .collect(),
        },
    })
}

#[derive(Clone, Debug)]
struct Reference {
    face: [f64; 2],
    temperature: [f64; 3],
    force: [f64; 2],
    moment: [f64; 2],
    surface: [f64; 2],
    bulk: [f64; 2],
    momentum: [f64; 3],
    angular: [f64; 3],
}
fn reference(f: &Fixture, old: [f64; 3], outside: [f64; 3], dt: f64) -> Reference {
    let masses = f.properties.map(|p| p.mass());
    let inertia = [0.4 * masses[0] * R * R, 0.13, 0.17];
    let capacity: [f64; 3] = std::array::from_fn(|i| masses[i] * CP[i]);
    let damping: [f64; 2] = std::array::from_fn(|i| f.tau * f.patches[i].normal_load_n() * R);
    let at = |face: [f64; 2]| {
        let force: [f64; 2] = std::array::from_fn(|i| {
            if f.speed == 0.0 && f.omega == 0.0 {
                0.0
            } else {
                -(super::MU[i] - super::SLOPE[i] * (face[i] - 300.0)) * f.patches[i].normal_load_n()
            }
        });
        let sphere_torque = (0..2).map(|i| SIGNS[i] * R * force[i]).sum::<f64>();
        let rhs: [f64; 2] = std::array::from_fn(|i| {
            f.omega
                + 0.5 * dt * (sphere_torque / inertia[0] - SIGNS[i] * D * force[i] / inertia[i + 1])
        });
        let a00 = 1.0 + 0.5 * dt * damping[0] * (inertia[0].recip() + inertia[1].recip());
        let a01 = 0.5 * dt * damping[1] / inertia[0];
        let a10 = 0.5 * dt * damping[0] / inertia[0];
        let a11 = 1.0 + 0.5 * dt * damping[1] * (inertia[0].recip() + inertia[2].recip());
        let determinant = a00 * a11 - a01 * a10;
        let relative = [
            (a11 * rhs[0] - a01 * rhs[1]) / determinant,
            (a00 * rhs[1] - a10 * rhs[0]) / determinant,
        ];
        let moment: [f64; 2] = std::array::from_fn(|i| -damping[i] * relative[i]);
        let momentum = [
            masses[0] * f.speed + dt * force.iter().sum::<f64>(),
            -dt * force[0],
            -dt * force[1],
        ];
        let angular = [
            inertia[0] * f.omega + dt * (sphere_torque + moment.iter().sum::<f64>()),
            dt * (SIGNS[0] * D * force[0] - moment[0]),
            dt * (SIGNS[1] * D * force[1] - moment[1]),
        ];
        let sphere_v = 0.5 * (f.speed + momentum[0] / masses[0]);
        let sphere_w = 0.5 * (f.omega + angular[0] / inertia[0]);
        let surface: [f64; 2] = std::array::from_fn(|i| {
            let slip = sphere_v - 0.5 * momentum[i + 1] / masses[i + 1]
                + SIGNS[i] * (R * sphere_w + 0.5 * D * angular[i + 1] / inertia[i + 1]);
            -dt * force[i] * slip
        });
        let bulk: [f64; 2] = std::array::from_fn(|i| dt * damping[i] * relative[i].powi(2));
        let temperature = [
            (capacity[0] * old[0]
                + dt * (G[0][0] * face[0] + G[1][0] * face[1] + outside[0])
                + bulk.iter().sum::<f64>())
                / (capacity[0] + dt * (G[0][0] + G[1][0])),
            (capacity[1] * old[1] + dt * (G[0][1] * face[0] + outside[1]))
                / (capacity[1] + dt * G[0][1]),
            (capacity[2] * old[2] + dt * (G[1][1] * face[1] + outside[2]))
                / (capacity[2] + dt * G[1][1]),
        ];
        let residual: [f64; 2] = std::array::from_fn(|i| {
            dt * (G[i][0] * (face[i] - temperature[0]) + G[i][1] * (face[i] - temperature[i + 1]))
                - surface[i]
        });
        (
            residual,
            Reference {
                face,
                temperature,
                force,
                moment,
                surface,
                bulk,
                momentum,
                angular,
            },
        )
    };
    let inner = |first| {
        let (mut lo, mut hi) = (280.0, 320.0);
        assert!(at([first, lo]).0[1] < 0.0 && at([first, hi]).0[1] > 0.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if at([first, mid]).0[1] > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        at([first, 0.5 * (lo + hi)])
    };
    let (mut lo, mut hi) = (280.0, 320.0);
    assert!(inner(lo).0[0] < 0.0 && inner(hi).0[0] > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if inner(mid).0[0] > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    inner(0.5 * (lo + hi)).1
}

fn verify(
    f: &Fixture,
    before: &DepositionSimulation<'_>,
    result: &HeatedRigidHertzRollingIslandCandidate<'_, Vec<f64>>,
    r: &Reference,
    outside: [f64; 3],
    dt: f64,
) {
    assert_eq!(result.mechanical.body_indices, vec![0, 1, 2]);
    assert_eq!(result.bulk_heat.len(), 1);
    assert_eq!(result.bulk_heat[0].binding, bulk_binding());
    near(
        result.bulk_heat[0].dissipated_energy_j,
        r.bulk.iter().sum(),
        8e-11,
    );
    for i in 0..2 {
        let native = &result.mechanical.report.contacts[i];
        assert_eq!(native.patch_id, PATCHES[i]);
        near(native.hertz.force_n[0], r.force[i], 4e-8);
        near(
            native.rolling_angular_impulse_world_nms.y,
            dt * r.moment[i],
            2e-11,
        );
        near(native.surface_heat_j, r.surface[i], 2e-10);
        near(native.bulk_heat_j, r.bulk[i], 8e-11);
        near(native.rolling_temperature_k, r.temperature[0], 2e-7);
        let history = &result.mechanical.rolling_histories[i];
        assert_eq!(history.patch_id, PATCHES[i]);
        near(
            history.history.state().cumulative_heat_j(),
            r.bulk[i],
            8e-11,
        );
        assert_eq!(history.history.state().accepted_steps(), 1);
        let restored = f
            .rolling
            .restore_checkpoint(
                &f.patches[i],
                f.sources[i].interface().interface(),
                &f.work[i],
                history.history.checkpoint().unwrap(),
            )
            .unwrap();
        assert_eq!(&restored, history.history.state());
    }
    for i in 0..3 {
        near(
            result.thermal.temperature_k(IDS[i]).unwrap(),
            r.temperature[i],
            2e-7,
        );
        near(
            result.mechanical.bodies_after[i].linear_momentum_world().x,
            r.momentum[i],
            5e-10,
        );
        near(
            result.mechanical.bodies_after[i].angular_momentum_body().y,
            r.angular[i],
            5e-11,
        );
    }
    for (j, i) in [1, 0].into_iter().enumerate() {
        near(result.source[j], r.face[i], 2e-7);
        assert_eq!(result.friction_sources[j].patch_id, PATCHES[i]);
        let radius = (f.patches[i].contact_area_m2() / core::f64::consts::PI).sqrt();
        assert_eq!(
            result.friction_sources[j].friction_source,
            f.sources[i]
                .hertz_at(result.source[j], None, f.patches[i].normal_load_n(), radius,)
                .unwrap()
        );
    }
    let dk: f64 = (0..3)
        .map(|i| {
            result.mechanical.bodies_after[i]
                .kinetic_energy(f.properties[i])
                .unwrap()
                - f.bodies[i].kinetic_energy(f.properties[i]).unwrap()
        })
        .sum();
    let dh: f64 = before
        .cells()
        .iter()
        .map(|(id, c)| c.mass_kg * (result.thermal.cells()[id].specific_h_j_kg - c.specific_h_j_kg))
        .sum();
    near(
        result.external_heat_j,
        dt * outside.iter().sum::<f64>(),
        1e-14,
    );
    near(dk + dh - dt * outside.iter().sum::<f64>(), 0.0, 7e-10);
    near(
        dh - dt * outside.iter().sum::<f64>()
            - r.surface.iter().sum::<f64>()
            - r.bulk.iter().sum::<f64>(),
        0.0,
        7e-10,
    );
    assert_eq!(before.time_s(), 0.0);
}

#[test]
fn shared_sphere_patches_solve_joint_recoil_and_one_original_bulk_inventory() {
    let f = Fixture::new(TAU, 0.2, 0.05);
    let initial = [299.0, 301.0, 303.0];
    let t = thermal(&f, initial, false);
    let routes = routes(&f, &t);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let result = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &f.heat_contacts(&routes, &histories),
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| source(&f, s, ids, faces, [0.0; 3]),
        &CancelGate::new(),
    )
    .unwrap();
    let r = reference(&f, initial, [0.0; 3], DT);
    assert!(r.surface.iter().all(|&d| d > 1e-3));
    assert!(r.bulk.iter().all(|&d| d > 1e-9));
    verify(&f, &t, &result, &r, [0.0; 3], DT);
    assert_eq!(
        histories,
        [RigidRollingHistory::new(), RigidRollingHistory::new()]
    );
}

#[test]
fn canonical_patch_order_preserves_each_face_source_and_signed_outside_heat() {
    let f = Fixture::new(TAU, 0.2, 0.0);
    let initial = [300.0, 298.0, 302.0];
    let outside = [0.15, -0.04, -0.03];
    let t = thermal(&f, initial, false);
    let routes = routes(&f, &t);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let mut contacts = f.heat_contacts(&routes, &histories);
    let first = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &contacts,
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| source(&f, s, ids, faces, outside),
        &CancelGate::new(),
    )
    .unwrap();
    contacts.reverse();
    let second = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &contacts,
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| source(&f, s, ids, faces, outside),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(first.mechanical, second.mechanical);
    assert_eq!(first.bulk_heat, second.bulk_heat);
    assert_eq!(first.friction_sources, second.friction_sources);
    assert_eq!(first.source, second.source);
    for id in IDS {
        assert_eq!(first.thermal.cells()[&id], second.thermal.cells()[&id]);
    }
    verify(
        &f,
        &t,
        &first,
        &reference(&f, initial, outside, DT),
        outside,
        DT,
    );
    assert!(
        first.mechanical.report.bulk_heat_j > 0.0,
        "shear induces rolling from rest"
    );
}

#[test]
fn stationary_and_lossless_material_limits_keep_original_conduction_and_work_channels() {
    for f in [Fixture::new(TAU, 0.0, 0.0), Fixture::new(0.0, 0.2, 0.05)] {
        let initial = [295.0, 305.0, 300.0];
        let t = thermal(&f, initial, false);
        let routes = routes(&f, &t);
        let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
        let result = propose_rigid_hertz_rolling_island_heat(
            &f.bodies,
            &f.properties,
            &f.heat_contacts(&routes, &histories),
            &t,
            0.0,
            DT,
            controls(),
            &[0.0; 3],
            |s, ids, faces| source(&f, s, ids, faces, [0.0; 3]),
            &CancelGate::new(),
        )
        .unwrap();
        verify(
            &f,
            &t,
            &result,
            &reference(&f, initial, [0.0; 3], DT),
            [0.0; 3],
            DT,
        );
        assert_eq!(result.mechanical.report.bulk_heat_j, 0.0);
        if f.speed == 0.0 {
            assert_eq!(result.mechanical.report.surface_heat_j, 0.0);
            assert_eq!(result.mechanical.bodies_after, f.bodies);
            assert!(result.thermal.temperature_k(IDS[0]).unwrap() > initial[0]);
            assert!(result.thermal.temperature_k(IDS[1]).unwrap() < initial[1]);
        }
    }
}

#[test]
fn conflicting_bulk_inventory_links_and_wrong_original_interface_refuse_before_publication() {
    let f = Fixture::new(TAU, 0.2, 0.05);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    for malformed in 0..6 {
        let t = thermal(&f, [300.0; 3], malformed == 0);
        let routes = routes(&f, &t);
        let mut contacts = f.heat_contacts(&routes, &histories);
        match malformed {
            0 => {}
            1 => contacts[1].bulk_binding.cell_id = IDS[2],
            2 => contacts[1].bulk_binding.thermal_body_id = BODY[2],
            3 => contacts[1].route = &routes[0],
            4 => contacts[1].contact.interface = f.sources[0].interface().interface(),
            5 => contacts[1].bulk_binding.body_index = 2,
            _ => unreachable!(),
        }
        assert!(
            propose_rigid_hertz_rolling_island_heat(
                &f.bodies,
                &f.properties,
                &contacts,
                &t,
                0.0,
                DT,
                controls(),
                &vec![0.0; t.cells().len()],
                |s, ids, faces| source(&f, s, ids, faces, [0.0; 3]),
                &CancelGate::new(),
            )
            .is_err(),
            "malformation {malformed}"
        );
        assert_eq!(t.time_s(), 0.0);
    }
    assert_eq!(
        histories,
        [RigidRollingHistory::new(), RigidRollingHistory::new()]
    );
}

#[test]
fn source_failure_trial_budget_cancellation_and_conductivity_identity_leave_exact_retry() {
    let f = Fixture::new(TAU, 0.2, 0.05);
    let t = thermal(&f, [300.0; 3], false);
    let routes = routes(&f, &t);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let contacts = f.heat_contacts(&routes, &histories);
    let run = |control, gate: &CancelGate| {
        propose_rigid_hertz_rolling_island_heat(
            &f.bodies,
            &f.properties,
            &contacts,
            &t,
            0.0,
            DT,
            control,
            &[0.0; 3],
            |s, ids, faces| source(&f, s, ids, faces, [0.0; 3]),
            gate,
        )
    };
    let baseline = run(controls(), &CancelGate::new()).unwrap();
    let calls = Cell::new(0usize);
    let failed = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &contacts,
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| {
            calls.set(calls.get() + 1);
            if calls.get() == 3 {
                Err("late source failure")
            } else {
                source(&f, s, ids, faces, [0.0; 3])
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidHertzRollingIslandHeatError::Source(
            "late source failure"
        ))
    ));
    let wrong_order = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &contacts,
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| {
            let mut x = source(&f, s, ids, faces, [0.0; 3])?;
            x.interfaces.reverse();
            Ok::<_, &'static str>(x)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        wrong_order,
        Err(RigidHertzRollingIslandHeatError::Invalid(_))
    ));
    let mut bounded = controls();
    bounded.maximum_trials = 1;
    assert!(matches!(
        run(bounded, &CancelGate::new()),
        Err(RigidHertzRollingIslandHeatError::TrialBudget)
    ));
    let mut bounded = controls();
    bounded.maximum_unknowns = 4;
    assert!(matches!(
        run(bounded, &CancelGate::new()),
        Err(RigidHertzRollingIslandHeatError::Invalid(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(run(controls(), &cancelled).is_err());
    let retry = run(controls(), &CancelGate::new()).unwrap();
    assert_eq!(baseline.mechanical, retry.mechanical);
    assert_eq!(baseline.source, retry.source);
    assert_eq!(baseline.bulk_heat, retry.bulk_heat);
    assert_eq!(
        histories,
        [RigidRollingHistory::new(), RigidRollingHistory::new()]
    );
    assert_eq!(t.time_s(), 0.0);
}

#[test]
fn two_distinct_spheres_keep_each_bulk_temperature_inventory_and_loss_destination() {
    // Patch order is lower/upper, native sphere order is upper/lower, and
    // thermal link order is upper/lower. This exercises a real nonconstant
    // patch-to-sphere mapping while the shared finite support couples recoil.
    let f = Fixture::new(TAU, 0.2, 0.0);
    let properties = [f.properties[0], f.properties[1], f.properties[0]];
    let positions = [R + D, 0.0, -R - D];
    let initial_speed = [0.22, 0.0, 0.17];
    let bodies: [RigidBodyState; 3] = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, positions[i]), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(properties[i].mass() * initial_speed[i], 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    let initial_temperature = [298.0, 301.0, 304.0];
    let ends = [[D, D + 2.0 * R], [-D, D], [-D - 2.0 * R, -D]];
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: BODY[i],
                material_id: if i == 1 { MATS[1] } else { MATS[0] },
                mass_kg: properties[i].mass(),
                specific_h_j_kg: CP[usize::from(i == 1)] * (initial_temperature[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, ends[i][0]],
                    end_m: [0.0, 0.0, ends[i][1]],
                    volume_m3: properties[i].mass() / 1000.0,
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
    for (i, sphere) in [2, 0].into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: IDS[sphere],
                b: IDS[1],
                area_m2: f.patches[i].contact_area_m2(),
                conductance_w_k: series_conductance_w_k(
                    f.patches[i].contact_area_m2(),
                    [R, D],
                    f.conductivity(i),
                )
                .unwrap(),
            }),
        });
    }
    let t = DepositionSimulation::new(
        events,
        &[
            (MATS[0], &CALORIC[0] as &dyn ThermalLaw),
            (MATS[1], &CALORIC[1]),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
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
    .0;
    let routes = routes(&f, &t);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let contacts = std::array::from_fn::<_, 2, _>(|i| {
        let sphere = [2, 0][i];
        let mut contact = f.contact(i);
        contact.binding = RigidContactBinding {
            body_a: sphere,
            body_b: 1,
            arm_a_body: Vec3::new(0.0, 0.0, -SIGNS[i] * R),
            arm_b_body: Vec3::new(0.0, 0.0, SIGNS[i] * D),
        };
        contact.frame = TangentFrame::new([0.0, 0.0, SIGNS[i]], [1.0, 0.0, 0.0]).unwrap();
        RigidHertzRollingIslandHeatContact {
            patch_id: PATCHES[i],
            contact,
            history: &histories[i],
            route: &routes[i],
            bulk_binding: RigidRollingHeatBinding {
                body_index: sphere,
                cell_id: IDS[sphere],
                thermal_body_id: BODY[sphere],
                material_id: MATS[0],
            },
            temperature_guess_k: 300.0,
        }
    });
    let result = propose_rigid_hertz_rolling_island_heat(
        &bodies,
        &properties,
        &contacts,
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |s, ids, faces| source(&f, s, ids, faces, [0.0; 3]),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        result
            .bulk_heat
            .iter()
            .map(|r| r.binding.body_index)
            .collect::<Vec<_>>(),
        [0, 2]
    );
    let endpoint = &result.mechanical.bodies_after;
    let inertia = [
        0.4 * properties[0].mass() * R * R,
        0.13,
        0.4 * properties[2].mass() * R * R,
    ];
    let mut sphere_losses = [0.0; 3];
    let mut contact_heat = [0.0; 3];
    for (i, sphere) in [2, 0].into_iter().enumerate() {
        let native = &result.mechanical.report.contacts[i];
        let face = result.source[1 - i];
        let sphere_t = result.thermal.temperature_k(IDS[sphere]).unwrap();
        let support_t = result.thermal.temperature_k(IDS[1]).unwrap();
        near(native.rolling_temperature_k, sphere_t, 1e-12);
        assert!((sphere_t - support_t).abs() > 0.5);
        let relative = 0.5
            * (endpoint[sphere].angular_momentum_body().y / inertia[sphere]
                - endpoint[1].angular_momentum_body().y / inertia[1]);
        let bulk = TAU * f.patches[i].normal_load_n() * R * DT * relative * relative;
        near(native.bulk_heat_j, bulk, 1e-12);
        assert!(bulk > 0.0);
        sphere_losses[sphere] += bulk;
        contact_heat[sphere] += DT * G[i][0] * (face - sphere_t);
        contact_heat[1] += DT * G[i][1] * (face - support_t);
        let force =
            -(super::MU[i] - super::SLOPE[i] * (face - 300.0)) * f.patches[i].normal_load_n();
        near(native.hertz.force_n[0], force, 3e-9);
        let midpoint_slip = 0.5
            * (initial_speed[sphere]
                + endpoint[sphere].linear_momentum_world().x / properties[sphere].mass()
                - endpoint[1].linear_momentum_world().x / properties[1].mass()
                - SIGNS[i] * R * endpoint[sphere].angular_momentum_body().y / inertia[sphere]
                - SIGNS[i] * D * endpoint[1].angular_momentum_body().y / inertia[1]);
        near(native.surface_heat_j, -DT * force * midpoint_slip, 2e-10);
    }
    for i in 0..3 {
        let caloric = properties[i].mass()
            * (result.thermal.cells()[&IDS[i]].specific_h_j_kg
                - t.cells()[&IDS[i]].specific_h_j_kg);
        near(caloric, contact_heat[i] + sphere_losses[i], 3e-10);
    }
    for bulk in &result.bulk_heat {
        near(
            bulk.dissipated_energy_j,
            sphere_losses[bulk.binding.body_index],
            1e-12,
        );
        near(
            bulk.temperature_before_k,
            initial_temperature[bulk.binding.body_index],
            1e-12,
        );
    }
    let change: f64 = (0..3)
        .map(|i| {
            endpoint[i].kinetic_energy(properties[i]).unwrap()
                - bodies[i].kinetic_energy(properties[i]).unwrap()
                + properties[i].mass()
                    * (result.thermal.cells()[&IDS[i]].specific_h_j_kg
                        - t.cells()[&IDS[i]].specific_h_j_kg)
        })
        .sum();
    near(change, 0.0, 7e-10);
    assert_eq!(
        histories,
        [RigidRollingHistory::new(), RigidRollingHistory::new()]
    );
}
