//! G1/G3/G4: original Hertz spin events and two finite caloric receivers.
//! The material curves below are explicit synthetic analytical fixtures.
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{SourceThermalCoulomb, ThermalFrictionAxes};
use fs_couple::friction_heat::sphere_plane::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::sphere_plane::hertz_support::{SpherePlaneHertzState, SpherePlaneHertzSupportControls};
use fs_mbd::sphere_plane::{SpherePlane, SpherePlaneConfig, SpherePlaneState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::InputAuthority;
use fs_tribo::sliding_patch::contact::HertzWrenchConfig;
use std::cell::Cell;
use std::collections::BTreeMap;

const R: f64 = 0.5;
const A: f64 = 0.1;
const I: f64 = 0.2;
const N: f64 = 20.0;
const LENGTH: f64 = 0.02;
const IDS: [CellId; 2] = [CellId(20), CellId(10)];
const BODIES: [u64; 2] = [11, 22];

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let ta = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / ta).ln_1p()))
    }
}
static CALORIC: Sensible = Sensible;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn friction(static_pair: [f64; 2], kinetic_pair: [f64; 2]) -> SourceThermalCoulomb {
    let mut claims = ClaimSet::new();
    for (name, values) in [("static", static_pair), ("kinetic", kinetic_pair)] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(290.0, values[0]), (330.0, values[1])],
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    290.0,
                    330.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic finite spin brake; not calibration".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        texture_frame: "isotropic".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("sphere"),
        surface("plane"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        QueryPoint::new()
            .with_quantity("T", temperature_quantity(), 300.0)
            .unwrap(),
        requirement("static"),
        requirement("kinetic"),
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap()
}
fn scene(normal: Vec3) -> SpherePlane {
    SpherePlane {
        properties: MassProperties::new(2.0, Vec3::ZERO, Vec3::new(I, I, I)).unwrap(),
        radius: R,
        point_world: Vec3::ZERO,
        normal_world: normal,
        velocity_world: Vec3::ZERO,
        force_world: normal.scale(-N),
        restitution: 0.0,
        friction_coefficient: 0.0,
    }
}
fn state(normal: Vec3, spin: f64) -> SpherePlaneHertzState {
    SpherePlaneHertzState::new(
        SpherePlaneState::new(
            RigidBodyState::new(
                Pose::new(normal.scale(R), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                normal.scale(I * spin),
            )
            .unwrap(),
            0.0,
        )
        .unwrap(),
    )
}
fn thermal(
    normal: Vec3,
    capacity: [f64; 2],
    conductance: [f64; 2],
    temperature: [f64; 2],
) -> DepositionSimulation<'static> {
    let area = core::f64::consts::PI * A * A;
    let mut events = Vec::new();
    for i in 0..2 {
        let signed = if i == 0 { 1.0 } else { -1.0 };
        let end = normal.scale(signed * LENGTH);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: BODIES[i],
                material_id: MaterialId(i as u64),
                mass_kg: capacity[i],
                specific_h_j_kg: temperature[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [end.x, end.y, end.z],
                    volume_m3: area * LENGTH,
                    surface_area_m2: 2.0 * area + 2.0 * core::f64::consts::PI * A * LENGTH,
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
            id: 7,
            a: IDS[0],
            b: IDS[1],
            area_m2: area,
            conductance_w_k: 1.0 / (1.0 / conductance[0] + 1.0 / conductance[1]),
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &CALORIC as &dyn ThermalLaw),
            (MaterialId(1), &CALORIC as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn conductivity(g: [f64; 2]) -> [f64; 2] {
    g.map(|x| x * 0.5 * LENGTH / (core::f64::consts::PI * A * A))
}
fn controls() -> SpherePlaneHertzHeatControls {
    SpherePlaneHertzHeatControls {
        mechanical: SpherePlaneHertzSupportControls {
            motion: SpherePlaneConfig::default(),
            wrench: HertzWrenchConfig {
                maximum_static_iterations: 96,
                maximum_iterations: 256,
                maximum_evaluations: 4096,
                force_tolerance_n: 2e-9,
                moment_tolerance_nm: 2e-10,
                displacement_tolerance_m: 2e-13,
                energy_tolerance_j: 2e-10,
                relative_matrix_tolerance: 1e-13,
                maximum_force_n: 30.0,
                maximum_spin_moment_nm: 5.0,
            },
            maximum_support_phase_s: 10.0,
            maximum_support_phases: 1024,
            velocity_tolerance_m_s: 2e-10,
            spin_tolerance_rad_s: 2e-10,
            work_tolerance_j: 2e-10,
            maximum_absolute_static_work_residual_j: 2e-9,
        },
        maximum_cells: 2,
        maximum_trials: 256,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-5,
        temperature_scale_j_per_k: 1.0,
        geometry_tolerance_m: 1e-12,
        heat: FrictionHeatTolerances {
            cell_j: 2e-9,
            thermal_j: 2e-9,
            interface_j: 2e-9,
            total_j: 2e-9,
        },
        combined_energy_tolerance_j: 5e-9,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-10,
            relative_tolerance: 1e-11,
            linear_restart: 3,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn contact(law: &SourceThermalCoulomb) -> SpherePlaneHertzHeatContact<'_> {
    SpherePlaneHertzHeatContact {
        source: law,
        contact_radius_m: A,
        thermal_body_ids: BODIES,
        nominal_pressure_pa: None,
    }
}
fn source(
    thermal: &DepositionSimulation<'_>,
    g: [f64; 2],
    external: [f64; 2],
) -> Result<SpherePlaneHertzHeatSource<[f64; 2]>, &'static str> {
    Ok(SpherePlaneHertzHeatSource {
        conductivity_w_m_k: conductivity(g),
        transport: StepCoefficients {
            sources_w: IDS.into_iter().zip(external).collect(),
            conductances_w_k: BTreeMap::new(),
            payload: IDS.map(|id| thermal.temperature_k(id).unwrap()),
        },
    })
}
#[allow(clippy::too_many_arguments)]
fn run(
    scene: &SpherePlane,
    initial: SpherePlaneHertzState,
    law: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'static>,
    g: [f64; 2],
    dt: f64,
    torque: f64,
    external: [f64; 2],
) -> HeatedSpherePlaneHertzCandidate<'static, [f64; 2]> {
    let route = FrictionHeatRoute::from_link(thermal, 7, conductivity(g)).unwrap();
    propose_sphere_plane_hertz_heat(
        scene,
        initial,
        scene.normal_world.scale(torque),
        contact(law),
        thermal,
        &route,
        dt,
        controls(),
        &[0.0; 2],
        300.0,
        |t, _| source(t, g, external),
        &CancelGate::new(),
    )
    .unwrap()
}
fn spin(result: &HeatedSpherePlaneHertzCandidate<'_, [f64; 2]>, normal: Vec3) -> f64 {
    let b = result.mechanical.motion.state_after.body();
    b.pose()
        .orientation()
        .rotate_body_to_world(b.angular_momentum_body())
        .dot(normal)
        / I
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; tolerance {tolerance:e}"
    );
}
fn audit(
    scene: &SpherePlane,
    before: SpherePlaneHertzState,
    heat_before: f64,
    r: &HeatedSpherePlaneHertzCandidate<'_, [f64; 2]>,
) {
    let m = &r.mechanical.motion;
    let delta_k = before
        .motion()
        .body()
        .kinetic_energy_change_to(scene.properties, m.state_after.body())
        .unwrap();
    close(
        delta_k + r.thermal.energy_j() - heat_before,
        r.heat.external_heat_j + m.external_torque_work + m.static_contact_work_residual,
        5e-9,
    );
    close(
        r.heat.dissipated_j,
        r.mechanical.sustained_friction_loss_j,
        0.0,
    );
    close(r.heat.dissipated_j, m.dissipated_energy, 0.0);
    close(r.combined_energy_residual_j, 0.0, 5e-9);
    assert_eq!(r.source, r.body_temperatures_after_k);
    assert_eq!(m.state_after.time_s(), r.thermal.time_s());
    assert_eq!(m.state_after.body().linear_momentum_world(), Vec3::ZERO);
    assert_eq!(m.external_force_work, 0.0);
    assert_eq!(m.plane_work, 0.0);
}

#[test]
fn actual_spin_stop_and_unequal_finite_body_heat_match_independent_closed_forms() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let s = scene(n);
    let initial = state(n, 4.0);
    let law = friction([0.6; 2], [0.2; 2]);
    let capacity = [0.5, 1.5];
    let g = [1.0, 3.0];
    let t0 = [310.0, 290.0];
    let t = thermal(n, capacity, g, t0);
    let r = run(&s, initial, &law, &t, g, 5.0, 0.0, [0.0; 2]);
    let moment = 3.0 * core::f64::consts::PI / 16.0 * N * A * 0.2;
    let stop = I * 4.0 / moment;
    let loss = 0.5 * I * 16.0;
    close(r.mechanical.transitions[0].time_s, stop, 2e-12);
    close(r.mechanical.support_spin_angle_rad, 2.0 * stop, 2e-12);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let expected_axis = UnitQuaternion::from_axis_angle(n, 2.0 * stop)
        .unwrap()
        .rotate_body_to_world(axis);
    let actual_axis = r
        .mechanical
        .motion
        .state_after
        .body()
        .pose()
        .orientation()
        .rotate_body_to_world(axis);
    close(actual_axis.x, expected_axis.x, 2e-12);
    close(actual_axis.y, expected_axis.y, 2e-12);
    close(actual_axis.z, expected_axis.z, 2e-12);
    close(spin(&r, n), 0.0, 2e-13);
    assert!(r.mechanical.state_after.patch_sticking());
    close(r.heat.dissipated_j, loss, 2e-12);
    let a = [0, 1].map(|i| capacity[i] / (capacity[i] + 5.0 * g[i]));
    let face =
        (loss / 5.0 + g[0] * a[0] * t0[0] + g[1] * a[1] * t0[1]) / (g[0] * a[0] + g[1] * a[1]);
    for i in 0..2 {
        close(
            r.body_temperatures_after_k[i],
            a[i] * t0[i] + (1.0 - a[i]) * face,
            5e-9,
        );
    }
    close(r.heat.interfaces[0].interface_temperature_k, face, 5e-9);
    assert!(r.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(r.heat.interfaces[0].body_heat_j[1] > loss);
    assert!(r.heat.entropy.is_some());
    assert_eq!(r.friction_source.owner_identity(), law.identity());
    assert_eq!(
        r.friction_source.friction_source().state_coordinate("T"),
        Some(r.heat.interfaces[0].interface_temperature_k)
    );
    for property in r.friction_source.friction_source().properties() {
        let original = law
            .initial_source()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
        law.interface()
            .card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    audit(&s, initial, t.energy_j(), &r);

    // The retained tiny native residual at arrest is not converted into heat
    // by the holding motor. Ordinary hot-to-cold conduction still continues.
    let held = run(
        &s,
        r.mechanical.state_after,
        &law,
        &r.thermal,
        g,
        0.7,
        0.3,
        [0.0; 2],
    );
    assert!(held.mechanical.state_after.patch_sticking());
    assert_eq!(held.heat.dissipated_j, 0.0);
    assert!(held.mechanical.absolute_static_work_residual_j > 0.0);
    close(held.thermal.energy_j(), r.thermal.energy_j(), 5e-9);
    audit(&s, r.mechanical.state_after, r.thermal.energy_j(), &held);
}

// Independent continuously coupled spin/bulk-temperature flow. A scalar
// algebraic elimination supplies face T, not the production thermal/Newton or
// native event solver. Here mu_k(T)=.45-.015(T-300), before the first stop.
fn reference() -> [f64; 5] {
    let capacity = [0.05, 0.15];
    let g = [0.2, 0.6];
    let factor = 3.0 * core::f64::consts::PI / 16.0 * N * A;
    let rate = |x: [f64; 5]| {
        let face = (g[0] * x[1] + g[1] * x[2] + factor * x[0] * (0.45 + 300.0 * 0.015))
            / (g[0] + g[1] + factor * x[0] * 0.015);
        let moment = factor * (0.45 - 0.015 * (face - 300.0));
        [
            -moment / I,
            g[0] * (face - x[1]) / capacity[0],
            g[1] * (face - x[2]) / capacity[1],
            x[0],
            moment * x[0],
        ]
    };
    let mut x = [4.0, 300.0, 300.0, 0.0, 0.0];
    let h = 1.0 / 8192.0;
    for _ in 0..8192 {
        let a = rate(x);
        let b = rate(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rate(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rate(std::array::from_fn(|i| x[i] + h * c[i]));
        x = std::array::from_fn(|i| x[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    x
}
#[test]
fn sourced_thermal_fade_converges_to_independent_continuous_spin_and_heat() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let s = scene(n);
    let law = friction([0.8, 0.2], [0.6, 0.0]);
    let expected = reference();
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut m = state(n, 4.0);
        let mut t = thermal(n, [0.05, 0.15], [0.2, 0.6], [300.0; 2]);
        let mut angle = 0.0;
        let mut heat = 0.0;
        for _ in 0..steps {
            let r = run(
                &s,
                m,
                &law,
                &t,
                [0.2, 0.6],
                1.0 / steps as f64,
                0.0,
                [0.0; 2],
            );
            audit(&s, m, t.energy_j(), &r);
            angle += r.mechanical.support_spin_angle_rad;
            heat += r.heat.dissipated_j;
            m = r.mechanical.state_after;
            t = r.thermal;
        }
        let omega = m.motion().body().angular_momentum_body().z / I;
        let actual = [
            omega,
            t.temperature_k(IDS[0]).unwrap(),
            t.temperature_k(IDS[1]).unwrap(),
            angle,
            heat,
        ];
        let error = actual
            .into_iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        eprintln!("thermal windows={steps}; physical state maximum error={error:e}");
        errors.push(error);
        // Heating weakens the same source enough to leave more spin than its
        // initial isothermal coefficient would predict over this duration.
        assert!(omega > 4.0 - (3.0 * core::f64::consts::PI / 16.0 * N * A * 0.45) / I);
    }
    assert!(errors[1] < 0.58 * errors[0]);
    assert!(errors[2] < 0.58 * errors[1]);
    assert!(errors[2] < 0.07);
}

#[test]
fn external_heat_weakens_static_retention_and_native_motor_work_drives_breakaway() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let s = scene(n);
    let law = friction([0.8, 0.2], [0.6, 0.0]);
    let m = state(n, 0.0);
    let t = thermal(n, [0.05, 0.15], [0.2, 0.6], [300.0; 2]);
    let held = run(&s, m, &law, &t, [0.2, 0.6], 0.1, 0.65, [0.0; 2]);
    assert!(held.mechanical.state_after.patch_sticking());
    assert_eq!(held.heat.dissipated_j, 0.0);
    let driven = run(
        &s,
        held.mechanical.state_after,
        &law,
        &held.thermal,
        [0.2, 0.6],
        0.25,
        0.65,
        [2.0, 6.0],
    );
    assert!(!driven.mechanical.state_after.patch_sticking());
    assert!(spin(&driven, n) > 0.0);
    assert!(driven.heat.dissipated_j > 0.0);
    assert!(driven.mechanical.motion.external_torque_work > driven.heat.dissipated_j);
    assert!(
        driven.friction_source.law().parameters().static_mu
            < held.friction_source.law().parameters().static_mu
    );
    close(driven.heat.external_heat_j, 2.0, 2e-12);
    audit(
        &s,
        held.mechanical.state_after,
        held.thermal.energy_j(),
        &driven,
    );
}

#[test]
fn opposed_motor_reversal_and_signed_axis_permutation_preserve_both_thermal_receivers() {
    let law = friction([0.6; 2], [0.2; 2]);
    let mut reference = None;
    for n in [Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, -1.0, 0.0)] {
        let s = scene(n);
        let initial = state(n, 1.0);
        let t = thermal(n, [0.5, 1.5], [1.0, 3.0], [300.0; 2]);
        let r = run(&s, initial, &law, &t, [1.0, 3.0], 0.4, -1.0, [0.0; 2]);
        let moment = 3.0 * core::f64::consts::PI / 16.0 * N * A * 0.2;
        let zero = I / (1.0 + moment);
        let after = (-1.0 + moment) * (0.4 - zero) / I;
        close(r.mechanical.transitions[0].time_s, zero, 1e-12);
        close(spin(&r, n), after, 1e-12);
        assert!(!r.mechanical.state_after.patch_sticking());
        audit(&s, initial, t.energy_j(), &r);
        let observed = [
            spin(&r, n),
            r.mechanical.support_spin_angle_rad,
            r.heat.dissipated_j,
            r.body_temperatures_after_k[0],
            r.body_temperatures_after_k[1],
        ];
        if let Some(expected) = reference {
            for (a, b) in observed.into_iter().zip(expected) {
                close(a, b, 2e-12);
            }
        }
        reference = Some(observed);
    }
}

#[test]
fn late_source_cancellation_and_geometry_refusals_preserve_the_complete_warmed_checkpoint() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let s = scene(n);
    let law = friction([0.8, 0.2], [0.6, 0.0]);
    let start = state(n, 4.0);
    let t = thermal(n, [0.05, 0.15], [0.2, 0.6], [300.0; 2]);
    let warm = run(&s, start, &law, &t, [0.2, 0.6], 0.1, 0.0, [0.0; 2]);
    let m = warm.mechanical.state_after;
    let t = warm.thermal;
    let saved = t.cells().clone();
    let saved_source = law.clone();
    let route = FrictionHeatRoute::from_link(&t, 7, conductivity([0.2, 0.6])).unwrap();
    let accepted = run(&s, m, &law, &t, [0.2, 0.6], 0.2, 0.0, [0.0; 2]);
    assert!(accepted.trials > 2);
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let calls = Cell::new(0usize);
        let failed = propose_sphere_plane_hertz_heat(
            &s,
            m,
            Vec3::ZERO,
            contact(&law),
            &t,
            &route,
            0.2,
            controls(),
            &[0.0; 2],
            300.0,
            |endpoint, _| {
                calls.set(calls.get() + 1);
                if calls.get() == accepted.trials {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("late material refusal");
                    }
                }
                source(endpoint, [0.2, 0.6], [0.0; 2])
            },
            &gate,
        );
        if cancel {
            assert!(matches!(
                failed,
                Err(SpherePlaneHertzHeatError::Thermal(
                    DepositionError::Cancelled
                ))
            ));
        } else {
            assert!(matches!(
                failed,
                Err(SpherePlaneHertzHeatError::Source("late material refusal"))
            ));
        }
        assert_eq!(calls.get(), accepted.trials);
        assert_eq!(t.cells(), &saved);
        assert_eq!(law, saved_source);
    }
    for case in 0..5 {
        let mut c = contact(&law);
        let mut sc = s;
        let mut control = controls();
        let torque = if case == 3 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::ZERO
        };
        match case {
            0 => c.contact_radius_m *= 0.9,
            1 => c.thermal_body_ids.swap(0, 1),
            2 => sc.velocity_world = Vec3::new(0.1, 0.0, 0.0),
            4 => control.maximum_trials = 1,
            _ => {}
        }
        assert!(
            propose_sphere_plane_hertz_heat(
                &sc,
                m,
                torque,
                c,
                &t,
                &route,
                0.2,
                control,
                &[0.0; 2],
                300.0,
                |e, _| source(e, [0.2, 0.6], [0.0; 2]),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(t.cells(), &saved);
    }
    let retry = run(&s, m, &law, &t, [0.2, 0.6], 0.2, 0.0, [0.0; 2]);
    assert_eq!(
        retry.mechanical.state_after,
        accepted.mechanical.state_after
    );
    assert_eq!(retry.thermal.cells(), accepted.thermal.cells());
    assert_eq!(retry.friction_source, accepted.friction_source);
}
