//! Actual native rigid shear/spin, retained storage and finite caloric cells.
//! Coefficients are analytical fixtures, not material-source claims.
use fs_contact::tangential::nonsmooth::{
    NonsmoothContact, NonsmoothControls, NonsmoothTangentialState,
};
use fs_couple::friction_heat::rigid::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{
    RigidContactBinding, RigidFrictionControls, apply_rigid_friction,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::*;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 0.01;
const ARM: f64 = 0.01;
const ELL: f64 = 0.01;
const AREA: f64 = core::f64::consts::PI * 0.02 * 0.02;
const NORMAL: f64 = 10.0;
const STIFFNESS: f64 = 10000.0;
const CAPACITIES: [f64; 2] = [0.02, 0.06];
const CONDUCTIVITY: [f64; 2] = [50.0, 100.0];
const MASSES: [f64; 2] = [2.0, 3.0];

struct Caloric;
impl ThermalLaw for Caloric {
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
}
static CALORIC: Caloric = Caloric;
fn law(t: f64) -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "synthetic/rigid-finite-heat",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4 - 0.02 * (t - 300.0),
            tangential_stiffness_n_per_m: STIFFNESS,
            torsional_stiffness_nm_per_rad: STIFFNESS * ELL * ELL,
            torsional_capacity_factor: ELL / 0.0002_f64.sqrt(),
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .unwrap()
}
fn patch() -> NormalPatchView {
    NormalPatchView::new(
        "patch",
        "synthetic/normal",
        "synthetic/normal-source",
        NormalPatchAuthority::SyntheticFixture,
        NORMAL,
        0.02,
        0.02,
        0.0002,
    )
    .unwrap()
}
fn interface() -> PartialSlipInterface {
    PartialSlipInterface::new(
        "body-a->body-b",
        "retained-shear-spin",
        "synthetic/interface",
        NormalPatchAuthority::SyntheticFixture,
    )
    .unwrap()
}
fn contact<'a>(
    law: &'a PartialSlipLaw,
    patch: &'a NormalPatchView,
    interface: &'a PartialSlipInterface,
    interval: &str,
) -> NonsmoothContact<'a> {
    NonsmoothContact {
        law,
        patch,
        interface,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        common_spin_rad_per_s: 0.0,
        torsional_length_m: ELL,
        creepage_reference_speed_mps: 1.0,
        work: GeneralizedWorkOwnership::new("patch", interval, "x", "y", "spin").unwrap(),
    }
}
fn binding() -> RigidContactBinding {
    RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: Vec3::new(0.0, 0.0, ARM),
        arm_b_body: Vec3::new(0.0, 0.0, -ARM),
    }
}
fn mechanics() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let properties = MASSES.map(|mass| {
        MassProperties::new(
            mass,
            Vec3::ZERO,
            Vec3::new(mass * ARM * ARM, mass * ARM * ARM, 0.5 * mass * ELL * ELL),
        )
        .unwrap()
    });
    let bodies = std::array::from_fn(|i| {
        let z = if i == 0 { -ARM } else { ARM };
        let omega = if i == 0 { 250.0 } else { -250.0 };
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, z), UnitQuaternion::IDENTITY).unwrap(),
            if i == 0 {
                Vec3::new(3.0 * MASSES[i], 4.0 * MASSES[i], 0.0)
            } else {
                Vec3::ZERO
            },
            Vec3::new(0.0, 0.0, properties[i].principal_inertia_body().z * omega),
        )
        .unwrap()
    });
    (bodies, properties)
}
fn thermal() -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..2 {
        let z = if i == 0 { -2.0 * ARM } else { 0.0 };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAPACITIES[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, z],
                    end_m: [0.0, 0.0, z + 2.0 * ARM],
                    volume_m3: AREA * 2.0 * ARM,
                    surface_area_m2: 2.0 * AREA + 2.0 * core::f64::consts::PI * 0.02 * 2.0 * ARM,
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
            a: CellId(0),
            b: CellId(1),
            area_m2: AREA,
            conductance_w_k: series_conductance_w_k(AREA, [ARM; 2], CONDUCTIVITY).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn controls() -> RigidFrictionHeatControls {
    let newton = NewtonKrylovConfig {
        absolute_tolerance: 1e-11,
        relative_tolerance: 1e-11,
        linear_restart: 3,
        forcing_minimum: 1e-6,
        ..Default::default()
    };
    RigidFrictionHeatControls {
        mechanical: RigidFrictionControls {
            maximum_bodies: 2,
            maximum_contacts: 1,
            geometry_tolerance_m: 1e-12,
            velocity_tolerance_mps: 1e-9,
            linear_momentum_tolerance_ns: 1e-10,
            angular_momentum_tolerance_nms: 1e-10,
            energy_tolerance_j: 1e-9,
            contact: NonsmoothControls {
                maximum_evaluations: 2000,
                maximum_newton_iterations: 30,
                maximum_mechanical_dofs: 12,
                difference_step_mps: 1e-5,
                velocity_tolerance_mps: 1e-10,
                cone_tolerance: 1e-10,
                energy_tolerance_j: 1e-9,
                newton,
            },
        },
        maximum_cells: 2,
        maximum_trials: 256,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.05,
        heat: FrictionHeatTolerances {
            cell_j: 1e-9,
            thermal_j: 1e-9,
            interface_j: 1e-9,
            total_j: 1e-9,
        },
        combined_energy_tolerance_j: 2e-9,
        newton,
    }
}
fn source(
    _: &DepositionSimulation<'_>,
    t: f64,
) -> Result<RigidFrictionHeatSource<f64>, DepositionError> {
    Ok(RigidFrictionHeatSource {
        law: law(t),
        conductivity_w_m_k: CONDUCTIVITY,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: t,
        },
    })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e} (tol {tolerance:e})"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn energy(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
        .sum()
}

#[test]
fn simultaneous_shear_spin_and_finite_heat_match_independent_closed_form() {
    let (mut bodies, properties) = mechanics();
    let before = bodies;
    let mut history = NonsmoothTangentialState::new(1, 16).unwrap();
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
    let original = law(300.0);
    let patch = patch();
    let interface = interface();
    let c = contact(&original, &patch, &interface, "first");
    let result = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &c,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();

    // ORIGINAL arms and inertias give isotropic mobility in [vx,vy,ell*w]:
    // W = 1/mA+1/mB+rA^2/I_Ay+rB^2/I_By = ell^2(1/I_Az+1/I_Bz).
    // In gross slip each radial Jenkins region has force C*f_i, stiffness k*w_i.
    // Plastic heat H=C*dt*s-C^2(dt^2*W+sum(f_i^2/w_i)/k).
    // This and the TWO finite caloric/conduction equations are an independent
    // scalar reduction; no production mechanical/thermal solver is used below.
    let mobility = 2.0 * (1.0 / MASSES[0] + 1.0 / MASSES[1]);
    let speed = 50.0_f64.sqrt();
    let regional = 0.3_f64.powi(2) / 0.6 + 0.7_f64.powi(2) / 0.4;
    let admittance = std::array::from_fn::<_, 2, _>(|i| {
        let dg = DT * CONDUCTIVITY[i] * AREA / ARM;
        CAPACITIES[i] * dg / (CAPACITIES[i] + dg)
    });
    let total_admittance: f64 = admittance.iter().sum();
    let heat_at = |t: f64| {
        let capacity = NORMAL * (0.4 - 0.02 * (t - 300.0));
        capacity * DT * speed - capacity * capacity * (DT * DT * mobility + regional / STIFFNESS)
    };
    let (mut low, mut high) = (300.0, 315.0);
    for _ in 0..60 {
        let mid = 0.5 * (low + high);
        if heat_at(mid) > total_admittance * (mid - 300.0) {
            low = mid;
        } else {
            high = mid;
        }
    }
    let expected_t = 0.5 * (low + high);
    let capacity = NORMAL * (0.4 - 0.02 * (expected_t - 300.0));
    let direction = [3.0 / speed, 4.0 / speed, 5.0 / speed];
    let endpoint_speed = speed - DT * mobility * capacity;
    let impulse = Vec3::new(
        -DT * capacity * direction[0],
        -DT * capacity * direction[1],
        0.0,
    );
    let couple = Vec3::new(0.0, 0.0, -DT * capacity * ELL * direction[2]);
    close(result.heat.interface_temperature_k, expected_t, 2e-7);
    close(result.source, expected_t, 2e-7);
    assert!(expected_t > 303.0 && expected_t < 305.0);
    assert!(result.law.parameters().kinetic_mu < 0.34);
    assert_eq!(
        result.mechanical.contacts[0].law.state,
        PartialSlipStateKind::GrossSlide
    );
    for i in 0..3 {
        close(
            result.mechanical.contacts[0].endpoint_velocity_mps[i],
            endpoint_speed * direction[i],
            2e-9,
        );
    }
    close_vec(
        bodies[0].linear_momentum_world(),
        before[0].linear_momentum_world().add(impulse),
        1e-10,
    );
    close_vec(
        bodies[1].linear_momentum_world(),
        before[1].linear_momentum_world().sub(impulse),
        1e-10,
    );
    close_vec(
        bodies[0].angular_momentum_body(),
        before[0]
            .angular_momentum_body()
            .add(binding().arm_a_body.cross(impulse))
            .add(couple),
        1e-11,
    );
    close_vec(
        bodies[1].angular_momentum_body(),
        before[1]
            .angular_momentum_body()
            .sub(binding().arm_b_body.cross(impulse))
            .sub(couple),
        1e-11,
    );
    close_vec(
        result.mechanical.linear_momentum_change_ns,
        Vec3::ZERO,
        1e-10,
    );
    close_vec(
        result.mechanical.angular_momentum_change_nms,
        Vec3::ZERO,
        1e-10,
    );
    close(result.heat.dissipated_j, heat_at(expected_t), 1e-9);
    close(energy(&thermal), heat_at(expected_t), 1e-9);
    close(
        result.heat.body_heat_j[0],
        admittance[0] * (expected_t - 300.0),
        1e-9,
    );
    close(
        result.heat.body_heat_j[1],
        admittance[1] * (expected_t - 300.0),
        1e-9,
    );
    let stored = 0.5 * capacity * capacity * regional / STIFFNESS;
    close(result.mechanical.friction.storage_change_j, stored, 1e-9);
    close(
        result.mechanical.friction.elastic_numerical_loss_j,
        stored,
        1e-9,
    );
    close(
        result.mechanical.friction.mechanical_numerical_loss_j,
        0.5 * DT * DT * mobility * capacity * capacity,
        1e-9,
    );
    assert!(
        result.mechanical.kinetic_before_j - result.mechanical.kinetic_after_j > energy(&thermal)
    );
    close(result.mechanical.work_energy_residual_j, 0.0, 1e-9);
    close(result.combined_energy_residual_j, 0.0, 2e-9);
    assert_eq!(result.mechanical.friction.coefficient_work_j, 0.0);
    assert_eq!(
        history.contact_checkpoint("patch").unwrap().law(),
        &result.law
    );
    assert_eq!(history.time_s(), thermal.time_s());
    assert_eq!(history.version(), 1);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), before[i].pose());
    }
}

#[test]
fn retained_temperature_continuation_is_atomic_on_source_stiffness_and_budget_refusals() {
    let (mut bodies, properties) = mechanics();
    let mut history = NonsmoothTangentialState::new(1, 32).unwrap();
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
    let original = law(300.0);
    let patch = patch();
    let interface = interface();
    let first = contact(&original, &patch, &interface, "first");
    let accepted = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &first,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let before = bodies;
    let history_before = history.clone();
    let thermal_before = thermal.clone();
    let second = contact(&accepted.law, &patch, &interface, "second");
    let calls = Cell::new(0usize);
    let rejected = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &second,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        |state, t| {
            calls.set(calls.get() + 1);
            if calls.get() == 3 {
                Err(DepositionError::MaterialDomain)
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        rejected,
        Err(RigidFrictionHeatError::Source(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(calls.get(), 3);
    assert_eq!(bodies, before);
    assert_eq!(history, history_before);
    assert_eq!(thermal.cells(), thermal_before.cells());
    assert_eq!(thermal.time_s(), thermal_before.time_s());

    let stiffness_refusal = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &second,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        |state, t| {
            let mut value = source(state, t)?;
            let mut p = value.law.parameters();
            p.tangential_stiffness_n_per_m *= 0.9;
            value.law =
                PartialSlipLaw::new(value.law.model_id(), value.law.source_id(), p).unwrap();
            Ok::<_, DepositionError>(value)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        stiffness_refusal,
        Err(RigidFrictionHeatError::UnbalancedCoefficientWork)
    ));
    assert_eq!(bodies, before);
    assert_eq!(history, history_before);
    assert_eq!(thermal.cells(), thermal_before.cells());
    let mut bounded = controls();
    bounded.maximum_trials = 1;
    let budget = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &second,
        binding(),
        &mut thermal,
        &route,
        DT,
        bounded,
        &[0.0; 2],
        accepted.source,
        source,
        &CancelGate::new(),
    );
    assert!(matches!(budget, Err(RigidFrictionHeatError::TrialBudget)));
    assert_eq!(bodies, before);
    assert_eq!(history, history_before);
    assert_eq!(thermal.cells(), thermal_before.cells());
    let gate = CancelGate::new();
    let cancelled = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &second,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        |state, t| {
            gate.request();
            source(state, t)
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(bodies, before);
    assert_eq!(history, history_before);
    assert_eq!(thermal.cells(), thermal_before.cells());

    let retry = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &second,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let mut fresh_bodies = before;
    let mut fresh_history = history_before.clone();
    let mut fresh_thermal = thermal_before.clone();
    let fresh = apply_rigid_friction_heat(
        &mut fresh_bodies,
        &properties,
        &mut fresh_history,
        &second,
        binding(),
        &mut fresh_thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, fresh_bodies);
    assert_eq!(history, fresh_history);
    assert_eq!(thermal.cells(), fresh_thermal.cells());
    assert_eq!(retry.source, fresh.source);
    assert!(retry.source > accepted.source);
    assert!(retry.law.parameters().kinetic_mu < accepted.law.parameters().kinetic_mu);
    assert_eq!(retry.mechanical.friction.coefficient_work_j, 0.0);
    assert_eq!(history.version(), 2);
    // Rerun the ORIGINAL mechanical owner from the exact accepted memory at the
    // solved temperature, checking retained history and its physical work.
    let mut reference_bodies = before;
    let mut reference_history = history_before;
    let resolved = contact(&retry.law, &patch, &interface, "second");
    let mechanical = apply_rigid_friction(
        &mut reference_bodies,
        &properties,
        &mut reference_history,
        &[resolved],
        &[binding()],
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, reference_bodies);
    assert_eq!(history, reference_history);
    assert_eq!(retry.heat.dissipated_j, mechanical.friction.physical_heat_j);
    close(
        energy(&thermal) - energy(&thermal_before),
        mechanical.friction.physical_heat_j,
        1e-9,
    );
    close(retry.combined_energy_residual_j, 0.0, 2e-9);
}

#[test]
fn geometry_mismatch_and_duplicate_work_preserve_every_accepted_owner() {
    let (mut bodies, properties) = mechanics();
    let before = bodies;
    let mut history = NonsmoothTangentialState::new(1, 8).unwrap();
    let original_history = history.clone();
    let mut thermal = thermal();
    let original_thermal = thermal.clone();
    let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
    let law = law(300.0);
    let patch = patch();
    let interface = interface();
    let c = contact(&law, &patch, &interface, "first");
    let mut wrong = binding();
    wrong.arm_a_body.x += 1e-4;
    let refusal = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &c,
        wrong,
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        source,
        &CancelGate::new(),
    );
    assert!(matches!(refusal, Err(RigidFrictionHeatError::Invalid(_))));
    assert_eq!(bodies, before);
    assert_eq!(history, original_history);
    assert_eq!(thermal.cells(), original_thermal.cells());
    let accepted = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &c,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let before = bodies;
    let old_history = history.clone();
    let old_thermal = thermal.clone();
    let refusal = apply_rigid_friction_heat(
        &mut bodies,
        &properties,
        &mut history,
        &c,
        binding(),
        &mut thermal,
        &route,
        DT,
        controls(),
        &[0.0; 2],
        accepted.source,
        source,
        &CancelGate::new(),
    );
    assert!(matches!(
        refusal,
        Err(RigidFrictionHeatError::Mechanical(_))
    ));
    assert_eq!(bodies, before);
    assert_eq!(history, old_history);
    assert_eq!(thermal.cells(), old_thermal.cells());
    assert_eq!(thermal.time_s(), old_thermal.time_s());
}

#[test]
fn sticking_and_partial_slip_keep_elastic_energy_out_of_physical_heat() {
    for (speed, kind) in [
        (0.01, PartialSlipStateKind::Sticking),
        (0.1, PartialSlipStateKind::PartialSlip),
    ] {
        let (initial, properties) = mechanics();
        let scale = speed / 50.0_f64.sqrt();
        let mut bodies = initial.map(|body| {
            RigidBodyState::new(
                body.pose(),
                body.linear_momentum_world().scale(scale),
                body.angular_momentum_body().scale(scale),
            )
            .unwrap()
        });
        let mut history = NonsmoothTangentialState::new(1, 8).unwrap();
        let mut thermal = thermal();
        let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
        let law = law(300.0);
        let patch = patch();
        let interface = interface();
        let c = contact(&law, &patch, &interface, "low-speed");
        let result = apply_rigid_friction_heat(
            &mut bodies,
            &properties,
            &mut history,
            &c,
            binding(),
            &mut thermal,
            &route,
            DT,
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(result.mechanical.contacts[0].law.state, kind);
        assert!(result.mechanical.friction.storage_change_j > 0.0);
        assert!(result.mechanical.friction.elastic_numerical_loss_j > 0.0);
        assert!(result.mechanical.friction.mechanical_numerical_loss_j > 0.0);
        close(result.combined_energy_residual_j, 0.0, 2e-9);
        if kind == PartialSlipStateKind::Sticking {
            assert_eq!(result.heat.dissipated_j, 0.0);
            assert_eq!(energy(&thermal), 0.0);
            assert_eq!(result.heat.interface_temperature_k, 300.0);
        } else {
            // Only the weaker region yields; the harder region's endpoint
            // spring force remains k_h*dt*v. The same native inertia gives W.
            let mobility = 2.0 * (1.0 / MASSES[0] + 1.0 / MASSES[1]);
            let weak_capacity = 0.3 * NORMAL * result.law.parameters().kinetic_mu;
            let velocity = (speed - DT * mobility * weak_capacity)
                / (1.0 + DT * DT * mobility * 0.4 * STIFFNESS);
            let heat = weak_capacity * (DT * velocity - weak_capacity / (0.6 * STIFFNESS));
            assert!(heat > 0.0);
            close(result.heat.dissipated_j, heat, 1e-9);
            close(energy(&thermal), heat, 1e-9);
        }
    }
}
