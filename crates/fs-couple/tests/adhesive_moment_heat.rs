//! Original adhesive moment memory, finite calorics and explicit loss conversion.
//! Constitutive thresholds and thermal coefficients are manufactured test data.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::rolling::{JkrRollingControls, SourceJkrRollingContact};
use fs_contact::normal_patch::adhesive::torsion::{JkrTorsionControls, SourceJkrTorsionContact};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_couple::friction_heat::adhesive::moment::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::{rolling::*, torsion::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1e-4;
const HALF_LENGTH: f64 = 1e-4;
const CELL_MASS: [f64; 2] = [2e-7, 3e-7];
const BODY_IDS: [u64; 2] = [70, 20];
const LINK: u64 = 39;

enum Source {
    Rolling(SourceJkrRollingContact),
    Torsion(SourceJkrTorsionContact),
}
impl Source {
    fn new(torsion: bool, tensile: bool) -> Self {
        let a = common::solid("moment-A", 1000.0, 1e7);
        let b = common::solid("moment-B", 2000.0, 2e7);
        let normal = SourceJkrSphereContact::new(
            common::bind(&a, &b, None, 0.05).normal_state().clone(),
            0.01,
            1e-9,
            JkrEnvelope {
                maximum_indentation_m: 1e-5,
                half_space_depth_m: 0.1,
                layer_thickness_m: 0.1,
                maximum_patch_to_radius: 0.1,
                maximum_patch_to_depth: 0.1,
                maximum_patch_to_layer: 0.1,
                maximum_nominal_strain: 0.1,
                maximum_repulsive_pressure_to_yield: 0.1,
                maximum_process_zone_to_patch: 0.01,
                minimum_tabor_parameter: 5.0,
            },
        )
        .unwrap();
        let load = if tensile {
            -0.75 * PI * 0.05 * 0.01
        } else {
            0.0
        };
        let approach = normal
            .evaluate_force_controlled(load)
            .unwrap()
            .indentation_m;
        if torsion {
            Self::Torsion(
                SourceJkrTorsionContact::new(
                    normal,
                    approach,
                    JkrTorsionControls {
                        critical_twist_angle_rad: 1e-5,
                        maximum_critical_twist_angle_rad: 0.01,
                        maximum_mode_iii_energy_release_to_adhesion: 0.01,
                    },
                )
                .unwrap(),
            )
        } else {
            Self::Rolling(
                SourceJkrRollingContact::new(
                    normal,
                    approach,
                    JkrRollingControls {
                        critical_rolling_displacement_m: 1e-8,
                        maximum_critical_displacement_to_contact_radius: 0.01,
                    },
                )
                .unwrap(),
            )
        }
    }
    fn normal(&self) -> &SourceJkrSphereContact {
        match self {
            Self::Rolling(s) => s.normal_source(),
            Self::Torsion(s) => s.normal_source(),
        }
    }
    fn area(&self) -> f64 {
        match self {
            Self::Rolling(s) => s.normal_response().contact_area_m2,
            Self::Torsion(s) => s.normal_response().contact_area_m2,
        }
    }
    fn angle(&self) -> f64 {
        match self {
            Self::Rolling(s) => s.law().critical_rolling_angle_rad(),
            Self::Torsion(s) => s.law().parameters().critical_twist_angle_rad,
        }
    }
    fn moment(&self) -> f64 {
        match self {
            Self::Rolling(s) => s.law().critical_moment_n_m(),
            Self::Torsion(s) => s.law().critical_moment_n_m(),
        }
    }
    fn scale(&self) -> f64 {
        self.angle() * self.moment()
    }
    fn axis(&self) -> Vec3 {
        match self {
            Self::Rolling(_) => Vec3::new(1.0, 0.0, 0.0),
            Self::Torsion(_) => Vec3::new(0.0, 0.0, 1.0),
        }
    }
    fn capacity(&self) -> [f64; 2] {
        [2.0 * self.scale(), 3.0 * self.scale()]
    }
    fn conductivity(&self) -> [f64; 2] {
        [4.0, 1.0].map(|v| v * self.scale() * HALF_LENGTH / (DT * self.area()))
    }
    fn controls(&self) -> AdhesiveMomentHeatControls {
        let e = self.scale();
        AdhesiveMomentHeatControls {
            mechanical: AdhesiveMomentMechanicalControls {
                maximum_bodies: 3,
                maximum_step_s: 10.0 * DT,
                maximum_angle_increment_rad: 100.0 * self.angle(),
                maximum_iterations: 96,
                angle_tolerance_rad: 1e-18,
                geometry_tolerance_m: 1e-15,
                angular_velocity_tolerance_rad_s: 1e-10,
                linear_momentum_tolerance_ns: 1e-22,
                angular_momentum_tolerance_nms: 1e-22,
                energy_tolerance_j: 1e-10 * e,
            },
            thermal: AdhesiveMomentThermalControls {
                maximum_cells: 2,
                maximum_interfaces: 1,
                maximum_unknowns: 3,
                maximum_trials: 192,
                maximum_newton_iterations: 12,
                difference_step_j: 1e-4 * e,
                temperature_scale_j_per_k: e,
                heat: FrictionHeatTolerances {
                    cell_j: 2e-7 * e,
                    thermal_j: 4e-7 * e,
                    interface_j: 2e-7 * e,
                    total_j: 4e-7 * e,
                },
                combined_energy_tolerance_j: 8e-7 * e,
                newton: NewtonKrylovConfig {
                    absolute_tolerance: 1e-9 * e,
                    relative_tolerance: 1e-10,
                    linear_restart: 3,
                    forcing_minimum: 1e-6,
                    ..Default::default()
                },
                audit_entropy: true,
            },
        }
    }
    fn curves(&self) -> [EquilibriumEnthalpyPhaseCurve; 2] {
        let elastic = self.normal().normal_state().elastic();
        let ids = [
            elastic.surface_a_state().resolved().card_identity(),
            elastic.surface_b_state().resolved().card_identity(),
        ];
        std::array::from_fn(|i| {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ids[i],
                SolidLiquidPhase::Solid,
                [280.0, 380.0]
                    .map(|t| EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: (t - 280.0) * self.capacity()[i] / CELL_MASS[i],
                        temperature_k: t,
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: [1000.0, 2000.0][i],
                    })
                    .to_vec(),
            )
            .unwrap()
        })
    }
    fn thermal<'t>(
        &self,
        curves: &'t [EquilibriumEnthalpyPhaseCurve; 2],
        initial: [f64; 2],
        area: f64,
    ) -> DepositionSimulation<'t> {
        let mut events: Vec<_> = (0..2)
            .map(|i| {
                let volume = CELL_MASS[i] / [1000.0, 2000.0][i];
                let face = volume / (2.0 * HALF_LENGTH);
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id: CellId(i as u64),
                        segment_id: i as u64,
                        body_id: BODY_IDS[i],
                        material_id: MaterialId(i as u64),
                        mass_kg: CELL_MASS[i],
                        specific_h_j_kg: (initial[i] - 280.0) * self.capacity()[i] / CELL_MASS[i],
                        geometry: CellGeometry {
                            start_m: [0.0, 0.0, if i == 0 { -2.0 * HALF_LENGTH } else { 0.0 }],
                            end_m: [0.0, 0.0, if i == 0 { 0.0 } else { 2.0 * HALF_LENGTH }],
                            volume_m3: volume,
                            surface_area_m2: 2.0 * face + 8.0 * face.sqrt() * HALF_LENGTH,
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
                conductance_w_k: series_conductance_w_k(
                    area,
                    [HALF_LENGTH; 2],
                    self.conductivity(),
                )
                .unwrap(),
            }),
        });
        DepositionSimulation::new(
            events,
            &[(MaterialId(0), &curves[0]), (MaterialId(1), &curves[1])],
            DepositionBudget {
                max_cells: 2,
                max_links: 1,
                max_events: 3,
                max_substeps: 4,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: DT,
                energy_tolerance_j: self.scale() * 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn resolve(
        &self,
        ids: &[u64],
        temperatures: &[f64],
        other_j: [f64; 2],
    ) -> Result<AdhesiveMomentHeatSource<f64>, &'static str> {
        assert_eq!(ids, &[LINK]);
        Ok(AdhesiveMomentHeatSource {
            interfaces: vec![AdhesiveMomentHeatConductivity {
                link_id: LINK,
                conductivity_w_m_k: self.conductivity(),
            }],
            transport: StepCoefficients {
                sources_w: other_j
                    .into_iter()
                    .enumerate()
                    .filter(|(_, q)| *q != 0.0)
                    .map(|(i, q)| (CellId(i as u64), q / DT))
                    .collect(),
                conductances_w_k: BTreeMap::new(),
                payload: temperatures[0],
            },
        })
    }
}

struct Pair {
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
}
impl Pair {
    fn new(source: &Source, rate: f64) -> Self {
        let a = MassProperties::new(2e-5, Vec3::ZERO, Vec3::new(0.8e-9, 1e-9, 1.2e-9)).unwrap();
        let b = MassProperties::new(3e-5, Vec3::ZERO, Vec3::new(1.5e-9, 2e-9, 2.5e-9)).unwrap();
        let properties = [b, a, a];
        let bodies = std::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(
                    match i {
                        0 => Vec3::new(0.0, 0.0, 0.01),
                        2 => Vec3::new(0.0, 0.0, -0.01),
                        _ => Vec3::new(1.0, 0.0, 0.0),
                    },
                    UnitQuaternion::IDENTITY,
                )
                .unwrap(),
                Vec3::ZERO,
                properties[i]
                    .angular_momentum_body_checked(source.axis().scale(if i == 2 {
                        rate
                    } else {
                        0.0
                    }))
                    .unwrap(),
            )
            .unwrap()
        });
        Self { bodies, properties }
    }
    fn native<'s>(&self, source: &'s Source) -> NativeAdhesiveMoment<'s> {
        let binding = RigidContactBinding {
            body_a: 2,
            body_b: 0,
            arm_a_body: Vec3::new(0.0, 0.0, 0.01),
            arm_b_body: Vec3::new(0.0, 0.0, -0.01),
        };
        let frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
        let gate = CancelGate::new_clock_free();
        match source {
            Source::Rolling(s) => {
                let contact = RigidAdhesiveRollingContact {
                    binding,
                    frame,
                    source: s,
                };
                let state = contact
                    .initial_state(
                        &self.bodies,
                        &self.properties,
                        source.controls().mechanical,
                        &gate,
                    )
                    .unwrap();
                NativeAdhesiveMoment::Rolling { contact, state }
            }
            Source::Torsion(s) => {
                let contact = RigidAdhesiveTorsionContact {
                    binding,
                    frame,
                    source: s,
                };
                let state = contact
                    .initial_state(
                        &self.bodies,
                        &self.properties,
                        source.controls().mechanical,
                        &gate,
                    )
                    .unwrap();
                NativeAdhesiveMoment::Torsion { contact, state }
            }
        }
    }
}

fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-8 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn history_angle(native: NativeAdhesiveMoment<'_>) -> f64 {
    match native {
        NativeAdhesiveMoment::Rolling { state, .. } => state.accumulated_angle_rad(),
        NativeAdhesiveMoment::Torsion { state, .. } => state.accumulated_angle_rad(),
    }
}
fn accounting(report: &AdhesiveMomentMechanicalReport) -> [f64; 5] {
    match report {
        AdhesiveMomentMechanicalReport::Rolling(r) => [
            r.kinetic_energy_change_j,
            r.contact.rolling().storage_change_j,
            r.contact.rolling().irreversible_loss_j,
            r.external_work_j,
            r.contact.rolling().storage_before_j,
        ],
        AdhesiveMomentMechanicalReport::Torsion(r) => [
            r.kinetic_energy_change_j,
            r.contact.torsion().storage_change_j,
            r.contact.torsion().irreversible_loss_j,
            r.external_work_j,
            r.contact.torsion().storage_before_j,
        ],
    }
}
fn heat_owner<'s>(
    source: &'s Source,
    pair: &Pair,
    curves: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
    fraction: f64,
) -> AdhesiveMomentHeat<'s> {
    AdhesiveMomentHeat::new(
        pair.native(source),
        BODY_IDS,
        [&curves[0], &curves[1]],
        FrictionHeatRoute::from_link(thermal, LINK, source.conductivity()).unwrap(),
        thermal,
        AdhesiveMomentHeatPartition::new(
            fraction,
            Some("irreversible junction restructuring".into()),
        )
        .unwrap(),
    )
    .unwrap()
}
fn exact_heat(
    source: &Source,
    before: [f64; 2],
    generated: f64,
    other: [f64; 2],
) -> (f64, [f64; 2], [f64; 2]) {
    let capacity = source.capacity();
    let conductance = source
        .conductivity()
        .map(|k| k * source.area() / HALF_LENGTH);
    let beta: [f64; 2] =
        std::array::from_fn(|i| DT * conductance[i] / (capacity[i] + DT * conductance[i]));
    let alpha: [f64; 2] = std::array::from_fn(|i| beta[i] * capacity[i]);
    let face = 300.0
        + (generated
            + (0..2)
                .map(|i| alpha[i] * (before[i] - 300.0) + beta[i] * other[i])
                .sum::<f64>())
            / alpha.iter().sum::<f64>();
    let after = std::array::from_fn(|i| {
        300.0
            + (capacity[i] * (before[i] - 300.0) + DT * conductance[i] * (face - 300.0) + other[i])
                / (capacity[i] + DT * conductance[i])
    });
    let heat = std::array::from_fn(|i| DT * conductance[i] * (face - after[i]));
    (face, after, heat)
}

#[test]
fn driven_tensile_rolling_and_torsion_heat_actual_cells_and_retain_noncaloric_loss() {
    for torsion in [false, true] {
        let source = Source::new(torsion, true);
        let curves = source.curves();
        for fraction in [0.0, 0.35, 1.0] {
            let mut pair = Pair::new(&source, 8.0 * source.angle() / DT);
            let initial = [320.0, 300.0];
            let mut thermal = source.thermal(&curves, initial, source.area());
            let mut history = heat_owner(&source, &pair, &curves, &thermal, fraction);
            let untouched = pair.bodies[1];
            let other = [0.0, 0.5 * source.scale()];
            let wrenches = [
                Wrench {
                    force_world: Vec3::ZERO,
                    torque_body: source.axis().scale(0.2 * source.moment()),
                },
                Wrench::ZERO,
            ];
            let report = history
                .apply(
                    &mut pair.bodies,
                    &pair.properties,
                    wrenches,
                    &mut thermal,
                    DT,
                    source.controls(),
                    &[0.0; 2],
                    310.0,
                    |_, ids, temperatures| source.resolve(ids, temperatures, other),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            let [kinetic, storage, loss, work, _] = accounting(&report.mechanical);
            assert!(loss > 0.0 && storage > 0.0 && work > 0.0);
            close(kinetic + storage + loss, work, 2e-9 * source.scale());
            close(
                report.losses.generated_heat_j(),
                fraction * loss,
                1e-13 * source.scale(),
            );
            close(
                report.losses.noncaloric_loss_j(),
                (1.0 - fraction) * loss,
                1e-13 * source.scale(),
            );
            assert_eq!(
                report.noncaloric_destination.as_deref(),
                Some("irreversible junction restructuring")
            );
            let (face, expected, expected_heat) =
                exact_heat(&source, initial, fraction * loss, other);
            let interface = &report.heat.interfaces[0];
            close(interface.interface_temperature_k, face, 1e-6);
            for i in 0..2 {
                close(
                    thermal.temperature_k(CellId(i as u64)).unwrap(),
                    expected[i],
                    1e-6,
                );
                close(
                    interface.body_heat_j[i],
                    expected_heat[i],
                    4e-7 * source.scale(),
                );
            }
            assert!(interface.body_heat_j[0] < 0.0 && interface.body_heat_j[1] > 0.0);
            close(
                report.heat.external_heat_j,
                other[1],
                1e-10 * source.scale(),
            );
            close(
                report.heat.thermal.stored_energy_change_j,
                fraction * loss + other[1],
                4e-7 * source.scale(),
            );
            close(
                report.heat.entropy.unwrap().friction_generation_j_k,
                fraction * loss / face,
                1e-8 * source.scale(),
            );
            assert!(report.combined_energy_residual_j.abs() < 8e-7 * source.scale());
            assert_eq!(pair.bodies[1], untouched);
            assert_eq!(history.time_s(), DT);
            assert_eq!(thermal.time_s(), DT);
            assert_eq!(history.losses(), report.losses);
            let previous = history.losses();
            let angle = history_angle(history.native());
            let next = history
                .apply(
                    &mut pair.bodies,
                    &pair.properties,
                    wrenches,
                    &mut thermal,
                    DT,
                    source.controls(),
                    &[0.0; 2],
                    face,
                    |_, ids, temperatures| source.resolve(ids, temperatures, [0.0; 2]),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            assert!(history_angle(history.native()) > angle);
            close(
                history.losses().irreversible_loss_j(),
                previous.irreversible_loss_j() + next.losses.irreversible_loss_j(),
                1e-13 * source.scale(),
            );
            close(
                history.losses().noncaloric_loss_j(),
                previous.noncaloric_loss_j() + next.losses.noncaloric_loss_j(),
                1e-13 * source.scale(),
            );
            close(
                next.heat.thermal.stored_energy_change_j,
                next.losses.generated_heat_j(),
                4e-7 * source.scale(),
            );
        }
    }
}

#[test]
fn elastic_storage_and_actual_external_work_are_not_heat() {
    for torsion in [false, true] {
        let source = Source::new(torsion, false);
        let curves = source.curves();
        let mut pair = Pair::new(&source, 0.2 * source.angle() / DT);
        let mut thermal = source.thermal(&curves, [300.0; 2], source.area());
        let mut history = heat_owner(&source, &pair, &curves, &thermal, 1.0);
        let force = 0.2 * (source.scale() * pair.properties[2].mass()).sqrt() / DT;
        let wrenches = [
            Wrench {
                force_world: Vec3::new(0.0, 0.0, force),
                torque_body: source.axis().scale(0.1 * source.moment()),
            },
            Wrench::ZERO,
        ];
        let report = history
            .apply(
                &mut pair.bodies,
                &pair.properties,
                wrenches,
                &mut thermal,
                DT,
                source.controls(),
                &[0.0; 2],
                300.0,
                |_, ids, temperatures| source.resolve(ids, temperatures, [0.0; 2]),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let [kinetic, storage, loss, work, _] = accounting(&report.mechanical);
        assert!(storage > 0.0 && work > 0.0);
        close(kinetic + storage, work, 2e-10 * source.scale());
        assert_eq!(loss, 0.0);
        assert_eq!(report.losses, AdhesiveMomentLoss::default());
        assert_eq!(report.heat.dissipated_j, 0.0);
        close(
            report.heat.thermal.stored_energy_change_j,
            0.0,
            1e-12 * source.scale(),
        );
        for i in 0..2 {
            close(thermal.temperature_k(CellId(i)).unwrap(), 300.0, 1e-10);
        }

        // Bind the actual loaded native state, then elastically unload it. No
        // state is rebuilt from angle or reports, and prior storage is not heat.
        let mut unloading = Pair::new(&source, -0.05 * source.angle() / DT);
        let loaded = history.native();
        let loaded_angle = history_angle(loaded);
        let mut continued = AdhesiveMomentHeat::new(
            loaded,
            BODY_IDS,
            [&curves[0], &curves[1]],
            history.route().clone(),
            &thermal,
            history.partition().clone(),
        )
        .unwrap();
        assert_eq!(history_angle(continued.native()), loaded_angle);
        let released = continued
            .apply(
                &mut unloading.bodies,
                &unloading.properties,
                [Wrench::ZERO; 2],
                &mut thermal,
                DT,
                source.controls(),
                &[0.0; 2],
                300.0,
                |_, ids, temperatures| source.resolve(ids, temperatures, [0.0; 2]),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let [kinetic, storage, loss, work, before] = accounting(&released.mechanical);
        assert!(before > 0.0 && storage < 0.0 && kinetic > 0.0);
        assert_eq!(loss, 0.0);
        assert_eq!(work, 0.0);
        assert_eq!(released.heat.dissipated_j, 0.0);
        close(kinetic, -storage, 2e-10 * source.scale());
    }
}

#[test]
fn thermal_refusal_and_late_cancellation_leave_all_owners_retryable() {
    for torsion in [false, true] {
        let source = Source::new(torsion, true);
        let curves = source.curves();
        let mut pair = Pair::new(&source, 8.0 * source.angle() / DT);
        let mut thermal = source.thermal(&curves, [300.0; 2], source.area());
        let mut history = heat_owner(&source, &pair, &curves, &thermal, 0.4);
        let bodies = pair.bodies;
        let angle = history_angle(history.native());
        let enthalpy: Vec<_> = thermal
            .cells()
            .values()
            .map(|c| c.specific_h_j_kg)
            .collect();
        let calls = Cell::new(0);
        let pilot = history
            .propose(
                &pair.bodies,
                &pair.properties,
                [Wrench::ZERO; 2],
                &thermal,
                DT,
                source.controls(),
                &[0.0; 2],
                300.0,
                |_, ids, temperatures| {
                    calls.set(calls.get() + 1);
                    source.resolve(ids, temperatures, [0.0; 2])
                },
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let final_call = calls.get();
        assert!(final_call > 2);
        assert_eq!(final_call, pilot.report().thermal_trials);
        for cancel in [false, true] {
            let gate = CancelGate::new_clock_free();
            calls.set(0);
            let refused = history.apply(
                &mut pair.bodies,
                &pair.properties,
                [Wrench::ZERO; 2],
                &mut thermal,
                DT,
                source.controls(),
                &[0.0; 2],
                300.0,
                |_, ids, temperatures| {
                    calls.set(calls.get() + 1);
                    if cancel && calls.get() == final_call {
                        gate.request();
                    }
                    let other = if cancel {
                        [0.0; 2]
                    } else {
                        [1e6 * source.scale(), 0.0]
                    };
                    source.resolve(ids, temperatures, other)
                },
                &gate,
            );
            assert!(refused.is_err());
            assert_eq!(pair.bodies, bodies);
            assert_eq!(history_angle(history.native()), angle);
            assert_eq!(history.losses(), AdhesiveMomentLoss::default());
            assert_eq!(history.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(
                thermal
                    .cells()
                    .values()
                    .map(|c| c.specific_h_j_kg)
                    .collect::<Vec<_>>(),
                enthalpy
            );
            if cancel {
                assert_eq!(calls.get(), final_call);
            }
        }
        let retry = history
            .apply(
                &mut pair.bodies,
                &pair.properties,
                [Wrench::ZERO; 2],
                &mut thermal,
                DT,
                source.controls(),
                &[0.0; 2],
                300.0,
                |_, ids, temperatures| source.resolve(ids, temperatures, [0.0; 2]),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(retry.losses, pilot.report().losses);
        assert_eq!(retry.mechanical, pilot.report().mechanical);
        assert_eq!(pair.bodies[2], pilot.bodies_after()[0]);
        assert_eq!(pair.bodies[0], pilot.bodies_after()[1]);
    }
}

#[test]
fn original_source_area_and_caloric_assignment_cannot_be_replaced() {
    assert!(AdhesiveMomentHeatPartition::new(0.9, None).is_err());
    assert!(AdhesiveMomentHeatPartition::new(0.5, Some(" ".into())).is_err());
    assert!(AdhesiveMomentHeatPartition::new(f64::NAN, None).is_err());
    for torsion in [false, true] {
        let source = Source::new(torsion, true);
        let other = Source::new(torsion, false);
        let curves = source.curves();
        let pair = Pair::new(&source, source.angle() / DT);
        let thermal = source.thermal(&curves, [300.0; 2], source.area());
        let route = FrictionHeatRoute::from_link(&thermal, LINK, source.conductivity()).unwrap();
        let partition = AdhesiveMomentHeatPartition::new(1.0, None).unwrap();
        let native = pair.native(&source);
        let foreign = match (native, pair.native(&other)) {
            (
                NativeAdhesiveMoment::Rolling { state, .. },
                NativeAdhesiveMoment::Rolling { contact, .. },
            ) => NativeAdhesiveMoment::Rolling { contact, state },
            (
                NativeAdhesiveMoment::Torsion { state, .. },
                NativeAdhesiveMoment::Torsion { contact, .. },
            ) => NativeAdhesiveMoment::Torsion { contact, state },
            _ => unreachable!(),
        };
        assert!(
            AdhesiveMomentHeat::new(
                foreign,
                BODY_IDS,
                [&curves[0], &curves[1]],
                route.clone(),
                &thermal,
                partition.clone()
            )
            .is_err()
        );
        assert!(
            AdhesiveMomentHeat::new(
                native,
                BODY_IDS,
                [&curves[1], &curves[0]],
                route.clone(),
                &thermal,
                partition.clone()
            )
            .is_err()
        );
        assert!(
            AdhesiveMomentHeat::new(
                native,
                [BODY_IDS[1], BODY_IDS[0]],
                [&curves[0], &curves[1]],
                route,
                &thermal,
                partition.clone()
            )
            .is_err()
        );
        let narrower = source.thermal(&curves, [300.0; 2], 0.5 * source.area());
        let narrowed_route =
            FrictionHeatRoute::from_link(&narrower, LINK, source.conductivity()).unwrap();
        assert!(
            AdhesiveMomentHeat::new(
                native,
                BODY_IDS,
                [&curves[0], &curves[1]],
                narrowed_route,
                &narrower,
                partition
            )
            .is_err()
        );
    }
}
