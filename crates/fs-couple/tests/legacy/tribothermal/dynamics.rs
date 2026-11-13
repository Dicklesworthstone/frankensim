//! Force-driven material contact, using the same source and finite heat fixtures.
#[path = "dynamics/collinear.rs"]
mod collinear;
#[path = "dynamics/planar.rs"]
mod planar;
#[path = "dynamics/directional.rs"]
mod directional;

use super::*;
use fs_ga::Vec3;
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::Density;
use fs_couple::tribothermal::dynamics::*;
use fs_thermal::enthalpy::{DepositionBudget, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent};

fn motion_controls() -> MotionSolve {
    let mut sliding = controls();
    // Resolve the native compliant return map more tightly than the outer
    // source/impulse gate, so its bisection error cannot mask source feedback.
    sliding.tangential.max_iterations = 64;
    sliding.tangential.force_absolute_tolerance_n = 1e-12;
    sliding.tangential.force_relative_tolerance = 1e-12;
    MotionSolve {
        sliding, maximum_motion_trials: 80, maximum_newton_iterations: 20,
        difference_step_mps: 1e-5, impulse_tolerance_n_s: 1e-8,
        total_energy_tolerance_j: 1e-6,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-10, relative_tolerance: 1e-12,
            linear_restart: 2, ..NewtonKrylovConfig::default() },
    }
}

fn inventory<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2], mass: [f64; 2]) -> DepositionSimulation<'a> {
    inventory_with_geometry(curves, mass, mass) // Existing fixtures have density 1 kg/m³.
}

fn inventory_with_geometry<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2], mass: [f64; 2], volume: [f64; 2],
) -> DepositionSimulation<'a> {
    let base = common::initial(curves, 0.0, None);
    let mut events: Vec<_> = CELLS.into_iter().zip(mass).zip(volume).map(|((id, mass_kg), volume_m3)| {
        let mut cell = base.cells()[&id].clone();
        cell.mass_kg = mass_kg;
        cell.geometry.volume_m3 = volume_m3;
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(base.links()[&7].clone()) });
    DepositionSimulation::new(events,
        &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
          (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 20 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
}

#[test]
fn original_body_masses_and_force_work_close_momentum_heat_and_kinetic_storage() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    // This reference is Coulomb with equal static/kinetic capacity; distinct
    // stiction thresholds are exercised independently below.
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = inventory(&curves, [2.0, 3.0]);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.4, 0.6, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    assert_eq!(pair.masses_kg(), [2.0, 3.0]);
    let next = pair.step_candidate(&thermal, 0.1,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[3.0, 0.0], [-1.0, 0.0]] },
        motion_controls(), [0.8, 0.0], &[0.0; 2], other, || false).unwrap();
    let expected = [1.0 + 0.1 * (3.0 - 2.0) / 2.0, 0.1 * (-1.0 + 2.0) / 3.0];
    for (v, exact) in next.pair().velocities_mps().into_iter().zip(expected) {
        close(v[0], exact); assert_eq!(v[1], 0.0);
    }
    let r = next.report();
    let heat = 2.0 * 0.1 * (expected[0] - expected[1]);
    let kinetic = 0.5 * (2.0 * (expected[0] * expected[0] - 1.0) + 3.0 * expected[1] * expected[1]);
    let numerical = 0.5 * (2.0 * (expected[0] - 1.0).powi(2) + 3.0 * expected[1].powi(2));
    close(r.kinetic_energy_change_j, kinetic);
    close(r.numerical_dissipation_j, numerical);
    close(r.external_work_j, 0.1 * (3.0 * expected[0] - expected[1]));
    close(r.sliding.deposited_heat_j, heat);
    close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + 0.4 * heat / 4.0);
    close(next.thermal().temperature_k(CELLS[1]).unwrap(), 300.0 + 0.6 * heat / 9.0);
    close(r.total_energy_residual_j, 0.0);
    assert!(r.body_momentum_residual_n_s.iter().flatten().all(|r| r.abs() < 1e-8));
    assert!(r.motion_trials > 1 && r.friction_trials > r.sliding.friction_trials);
    // The scheme's kinetic defect is independently visible, and not body heat.
    assert!(r.numerical_dissipation_j > 1e-3);
    assert!((r.sliding.thermal.stored_energy_change_j - heat - numerical).abs() > 1e-3);
    close(next.pair().kinetic_energy_j().unwrap() - pair.kinetic_energy_j().unwrap(), kinetic);
    close(next.pair().contact().work().dissipated_work_j(), heat);
    assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
}

#[test]
fn heat_capacity_data_changes_actual_slowing_through_endpoint_material_friction() {
    let mut speeds = Vec::new();
    for cp_a in [0.1, 1.0] {
        let (a, ca) = common::material(1.0, cp_a, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b]; let curves = [ca, cb]; let point = context();
        // mu(Ta) = 0.2 + 0.1*(Ta-300), explicitly within its 299..310 K range.
        // Static capacity 4 N is below the 5 N reaction needed to arrest this
        // interval, and above the reached temperature-dependent kinetic force.
        let card = interface_with_static(&cards, 0.4, "Ta", Temperature::DIMS, [(299.0, 0.1), (310.0, 1.2)], 2000.0);
        let thermal = common::initial(&curves, 0.0, None);
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
            motion_controls(), [0.6, 0.0], &[0.0; 2], other, || false).unwrap();
        // Independent elimination: mu=.2/(1-k*s), s=1-2*mu, k=.05/Ca.
        // Select the positive branch continuous from k=0 without cancellation.
        let k = 0.05 / cp_a;
        let exact = 1.2 / (1.0 + k + ((1.0 + k).powi(2) - 2.4 * k).sqrt());
        let mu = (1.0 - exact) / 2.0;
        let velocity = next.pair().velocities_mps();
        let speed = velocity[0][0] - velocity[1][0];
        close(speed, exact);
        close(velocity[0][0], 0.5 * (1.0 + exact));
        close(velocity[1][0], 0.5 * (1.0 - exact));
        assert_eq!([velocity[0][1], velocity[1][1]], [0.0; 2]);
        close(next.report().sliding.material.speed_mps, speed);
        close(next.report().sliding.friction.traction_n[0], 10.0 * mu);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + 0.5 * mu * exact / cp_a);
        close(next.pair().contact().accumulated_distance_m(), 0.1 * exact);
        close(next.report().sliding.deposited_heat_j, mu * exact);
        close(next.report().numerical_dissipation_j, mu * mu);
        close(next.report().total_energy_residual_j, 0.0);
        assert_eq!(pair.contact().binding().selection, next.pair().contact().binding().selection);
        speeds.push(speed);
    }
    assert!(speeds[1] - speeds[0] > 0.1);
}

#[test]
fn two_coordinate_friction_surface_drives_actual_motion_and_finite_heat() {
    let (a, ca) = common::material(1.0, 0.5, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = common::initial(&curves, 0.0, None);
    let template = interface_with_static(
        &cards,
        0.4,
        "speed",
        SPEED,
        [(0.0, 0.2), (2.0, 0.2)],
        2000.0,
    );
    let mut endpoints = Vec::new();
    for beta in [0.0, 0.15] {
        let mut claims = ClaimSet::new();
        for (_, original) in template.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == KINETIC_FRICTION_PROPERTY {
                claim.value = PropertyValue::Surface2D {
                    axes: ["Ta".into(), "speed".into()],
                    axis_dims: [Temperature::DIMS, SPEED],
                    samples: [[299.0, 0.0], [310.0, 0.0], [310.0, 2.0], [299.0, 2.0]]
                        .into_iter()
                        .map(|[t, s]| ([t, s], 0.2 + 0.1 * (t - 300.0) + beta * s))
                        .collect(),
                    triangles: vec![[0, 1, 2], [0, 2, 3]],
                    dims: Dims::NONE,
                };
                claim.validity = claim.validity.with_quantity(
                    "Ta",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    299.0,
                    310.0,
                );
            }
            claims.insert_claim(claim).unwrap();
        }
        let card = InterfaceSystemCard::assemble(
            template.surface_a().clone(),
            template.surface_b().clone(),
            template.context().clone(),
            claims,
            vec![],
        )
        .unwrap();
        let contact = MaterialSlidingContact::new(
            &thermal,
            binding(
                &cards,
                &curves,
                &card,
                &point,
                HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            ),
            10.0,
            1.0,
            InputAuthority::SyntheticFixture,
        )
        .unwrap();
        let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
        let next = pair
            .step_candidate(
                &thermal,
                0.1,
                MotionDrive {
                    normal_load_n: 10.0,
                    body_forces_n: [[0.0; 2]; 2],
                },
                motion_controls(),
                [0.5, 0.0],
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        // mA=mB=1, N=10, dt=.1 imply s=1-2*mu and heat=mu*s.
        // Ca=.5 and half the heat to A give Ta-300=mu*s. Thus
        // .1*s^2-(1.1+2*beta)*s+.6=0; take its positive small root.
        let b = 1.1_f64 + 2.0 * beta;
        let speed = 1.2 / (b + (b * b - 0.24).sqrt());
        let mu = 0.5 * (1.0 - speed);
        let heat = mu * speed;
        let ta = 300.0 + heat;
        let tb = 300.0 + heat / 6.0;
        let velocities = next.pair().velocities_mps();
        close(velocities[0][0], 1.0 - mu);
        close(velocities[1][0], mu);
        assert_eq!([velocities[0][1], velocities[1][1]], [0.0; 2]);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), ta);
        close(next.thermal().temperature_k(CELLS[1]).unwrap(), tb);
        let report = next.report();
        let material = &report.sliding.material.resolved;
        let property = material.property(KINETIC_FRICTION_PROPERTY).unwrap();
        close(property.value_si(), 0.2 + 0.1 * (ta - 300.0) + beta * speed);
        close(property.value_si(), mu);
        close(material.state_coordinate("Ta").unwrap(), ta);
        close(
            material.state_coordinate("speed").unwrap(),
            velocities[0][0] - velocities[1][0],
        );
        assert_eq!(
            property.source_axes(),
            &[("Ta".into(), Temperature::DIMS), ("speed".into(), SPEED)]
        );
        assert_eq!(
            property.answer().receipt.selected,
            card.claims_for(KINETIC_FRICTION_PROPERTY)[0].0
        );
        assert_eq!(
            pair.contact().binding().selection,
            next.pair().contact().binding().selection
        );
        close(report.sliding.friction.traction_n[0], 10.0 * mu);
        close(report.sliding.deposited_heat_j, heat);
        close(report.sliding.thermal.stored_energy_change_j, heat);
        close(report.numerical_dissipation_j, mu * mu);
        close(report.total_energy_residual_j, 0.0);
        assert!(
            report
                .body_momentum_residual_n_s
                .iter()
                .flatten()
                .all(|r| r.abs() < 1e-8)
        );
        close(next.pair().contact().accumulated_distance_m(), 0.1 * speed);
        assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
        assert_eq!(pair.contact().work().dissipated_work_j(), 0.0);
        endpoints.push((velocities[0][0] - velocities[1][0], ta));
    }
    assert!(endpoints[0].0 - endpoints[1].0 > 0.1);
    assert!((endpoints[0].1 - endpoints[1].1).abs() > 1e-4);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn late_failures_and_budget_refuse_before_any_motion_or_heat_is_published() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let original_cells = thermal.cells().clone();
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    let drive = MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] };
    let polls = Cell::new(0);
    let good = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), [0.8, 0.0], &[0.0; 2], other,
        || { polls.set(polls.get() + 1); false }).unwrap();
    let last = polls.get(); polls.set(0);
    assert!(matches!(pair.step_candidate(&thermal, 0.1, drive, motion_controls(), [0.8, 0.0], &[0.0; 2], other,
        || { polls.set(polls.get() + 1); polls.get() == last }), Err(SlidingError::Cancelled)));
    let mut limited = motion_controls(); limited.maximum_motion_trials = good.report().motion_trials - 1;
    assert!(matches!(pair.step_candidate(&thermal, 0.1, drive, limited, [0.8, 0.0], &[0.0; 2], other, || false),
        Err(SlidingError::MotionBudget)));
    let refused = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), [0.8, 0.0], &[0.0; 2], |state| {
        if state.temperature_k(CELLS[0]).unwrap() > 300.01 { return Err("thermal transport refused"); }
        Ok(StepCoefficients { sources_w: BTreeMap::new(), conductances_w_k: BTreeMap::new(), payload: () })
    }, || false);
    assert!(matches!(refused, Err(SlidingError::Source("thermal transport refused"))));
    let changed_inventory = inventory(&curves, [2.0, 1.0]);
    assert!(matches!(pair.step_candidate(&changed_inventory, 0.1, drive, motion_controls(), [0.8, 0.0], &[0.0; 2], other, || false),
        Err(SlidingError::Invalid(_))));
    assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
    assert_eq!(pair.contact().work().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    let retry = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), [0.8, 0.0], &[0.0; 2], other, || false).unwrap();
    assert_eq!(retry.thermal().cells(), good.thermal().cells());
    assert_eq!(retry.pair().velocities_mps(), good.pair().velocities_mps());
    assert_eq!(retry.pair().contact().state(), good.pair().contact().state());
    let (thermal, pair, _) = retry.into_parts();
    let continued = pair.step_candidate(&thermal, 0.2, drive, motion_controls(), [0.4, 0.0],
        &[0.0; 2], other, || false).unwrap();
    close(continued.pair().velocities_mps()[0][0], 0.6);
    close(continued.pair().velocities_mps()[1][0], 0.4);
    close(continued.pair().contact().accumulated_distance_m(), 0.08);
    close(continued.pair().contact().work().dissipated_work_j(), 0.16);
    close(continued.report().total_energy_residual_j, 0.0);
}

#[test]
fn native_static_capacity_drives_rest_below_at_and_above_breakaway() {
    use fs_tribo::FrictionRegime;
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.4, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 0.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.0; 2]; 2]).unwrap();
    for force in [6.0, 8.0, 10.0] {
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [[force, 0.0], [0.0; 2]] },
            motion_controls(), [0.3, 0.0], &[0.0; 2], other, || false).unwrap();
        let r = next.report();
        let v = next.pair().velocities_mps();
        if force <= 8.0 {
            // Native static traction F/2 gives the common acceleration F/(mA+mB).
            assert_eq!(r.sliding.friction.regime, FrictionRegime::Sticking);
            assert_eq!(v[0], v[1]);
            close(v[0][0], 0.05 * force);
            assert_eq!([v[0][1], v[1][1]], [0.0; 2]);
            close(r.sliding.friction.traction_n[0], force / 2.0);
            assert_eq!(r.sliding.friction.slip_increment_m, [0.0; 2]);
            assert_eq!(r.sliding.material.speed_mps, 0.0);
            assert_eq!(r.sliding.deposited_heat_j, 0.0);
            assert_eq!(next.pair().contact().accumulated_distance_m(), 0.0);
            assert_eq!(next.thermal().cells(), thermal.cells());
        } else {
            assert_eq!(r.sliding.friction.regime, FrictionRegime::Sliding);
            close(v[0][0], 0.8); close(v[1][0], 0.2);
            close(r.sliding.friction.traction_n[0], 2.0);
            close(r.sliding.friction.slip_increment_m[0], 0.06);
            close(r.sliding.deposited_heat_j, 0.12);
        }
        assert!(r.body_momentum_residual_n_s.iter().flatten().all(|x| x.abs() < 1e-8));
        assert!(r.source_motion_residual_n_s.iter().all(|x| x.abs() < 1e-8));
        close(r.total_energy_residual_j, 0.0);
    }
}

#[test]
fn signed_source_speed_closes_against_the_native_solved_slip() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.4, "speed", SPEED, [(0.0, 0.1), (2.0, 0.5)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    for sign in [1.0, -1.0] {
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::new(&thermal, contact, [[sign, 0.0], [0.0; 2]]).unwrap();
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
            motion_controls(), [sign * 0.5, 0.0], &[0.0; 2], other, || false).unwrap();
        // s=1-2*(.1+.2*s), so |s|=4/7; querying released speed would differ.
        let r = next.report();
        close(next.pair().velocities_mps()[0][0] - next.pair().velocities_mps()[1][0], sign * 4.0 / 7.0);
        close(r.sliding.material.speed_mps, 4.0 / 7.0);
        close(r.sliding.friction.traction_n[0], sign * 15.0 / 7.0);
        close(r.sliding.material.accumulated_distance_m, 0.1 * 4.0 / 7.0);
        close(r.sliding.deposited_heat_j, 6.0 / 49.0);
        close(r.total_energy_residual_j, 0.0);
        assert!(r.source_motion_residual_n_s.iter().all(|x| x.abs() < 1e-8));
    }
}

#[test]
fn exact_native_stop_then_opposite_restart_retains_heat_and_signed_history() {
    use fs_tribo::FrictionRegime;
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    let stopped = pair.step_candidate(&thermal, 0.25,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
        motion_controls(), [0.2, 0.0], &[0.0; 2], other, || false).unwrap();
    assert_eq!(stopped.report().sliding.friction.regime, FrictionRegime::Sticking);
    assert_eq!(stopped.pair().velocities_mps(), [[0.5, 0.0]; 2]);
    assert_eq!(stopped.report().sliding.material.speed_mps, 0.0);
    assert_eq!(stopped.report().sliding.deposited_heat_j, 0.0);
    close(stopped.report().sliding.friction.traction_n[0], 2.0);
    close(stopped.report().kinetic_energy_change_j, -0.25);
    close(stopped.report().numerical_dissipation_j, 0.25);
    close(stopped.report().total_energy_residual_j, 0.0);
    // Endpoint stopping is BE numerical loss, not invented frictional heating.
    assert_eq!(stopped.thermal().cells(), thermal.cells());
    let (thermal, pair, _) = stopped.into_parts();
    let reversed = pair.step_candidate(&thermal, 0.35,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[-10.0, 0.0], [0.0; 2]] },
        motion_controls(), [-0.5, 0.0], &[0.0; 2], other, || false).unwrap();
    let r = reversed.report();
    assert_eq!(r.sliding.friction.regime, FrictionRegime::Sliding);
    close(reversed.pair().velocities_mps()[0][0], -0.3);
    close(reversed.pair().velocities_mps()[1][0], 0.3);
    close(r.sliding.friction.traction_n[0], -2.0);
    close(r.sliding.friction.slip_increment_m[0], -0.06);
    close(r.sliding.material.speed_mps, 0.6);
    close(r.sliding.material.accumulated_distance_m, 0.06);
    close(reversed.pair().contact().state().slip_m()[0], -0.06);
    close(reversed.pair().contact().work().dissipated_work_j(), 0.12);
    close(r.sliding.deposited_heat_j, 0.12);
    close(r.total_energy_residual_j, 0.0);
    // A zero-speed source limit is actually queried; no nonzero speed substitute.
    let unsupported = interface_with_static(&cards, 0.2, "speed", SPEED,
        [(0.1, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &unsupported, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 0.5, InputAuthority::SyntheticFixture).unwrap();
    let bounded = SlidingPair::new(&thermal, contact, [[0.5, 0.0]; 2]).unwrap();
    assert!(matches!(bounded.step_candidate(&thermal, 0.35,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
        motion_controls(), [0.0; 2], &[0.0; 2], other, || false), Err(SlidingError::Material(_))));
    assert_eq!(bounded.velocities_mps(), [[0.5, 0.0]; 2]);
    assert_eq!(bounded.contact().work().dissipated_work_j(), 0.0);
}

#[test]
fn oblique_force_driven_sliding_closes_both_impulses_and_uses_speed_norm() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED,
        [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let masses = [2.0, 3.0];
    let thermal = inventory(&curves, masses);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.4, 0.6, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let initial = [[1.0, 0.0], [0.0; 2]];
    let forces = [[0.0, 10.0], [-3.0, -3.0]];
    let pair = SlidingPair::new(&thermal, contact, initial).unwrap();
    let dt = 0.1;
    let next = pair.step_candidate(&thermal, dt,
        MotionDrive { normal_load_n: 10.0, body_forces_n: forces },
        motion_controls(), [0.9, 0.5], &[0.0; 2], other, || false).unwrap();
    // Isotropic Coulomb force is radial to the released relative velocity.
    // The applied y force is perpendicular to the initial relative motion.
    let released: [f64; 2] = core::array::from_fn(|axis|
        initial[0][axis] - initial[1][axis]
            + dt * (forces[0][axis] / masses[0] - forces[1][axis] / masses[1]));
    let free_speed = released[0].hypot(released[1]);
    let traction = released.map(|component| 2.0 * component / free_speed);
    let expected: [[f64; 2]; 2] = core::array::from_fn(|body| core::array::from_fn(|axis|
        initial[body][axis] + dt * (forces[body][axis]
            + if body == 0 { -traction[axis] } else { traction[axis] }) / masses[body]));
    let speed = free_speed - dt * (1.0 / masses[0] + 1.0 / masses[1]) * 2.0;
    let heat = 2.0 * dt * speed;
    let mut kinetic = 0.0;
    let mut numerical = 0.0;
    let mut work = 0.0;
    let report = next.report();
    for body in 0..2 {
        for axis in 0..2 {
            let velocity = next.pair().velocities_mps()[body][axis];
            close(velocity, expected[body][axis]);
            close(report.body_displacements_m[body][axis], dt * velocity);
            let signed_traction = if body == 0 { -traction[axis] } else { traction[axis] };
            close(masses[body] * (velocity - initial[body][axis]),
                dt * (forces[body][axis] + signed_traction));
            kinetic += 0.5 * masses[body] * (velocity.powi(2) - initial[body][axis].powi(2));
            numerical += 0.5 * masses[body] * (velocity - initial[body][axis]).powi(2);
            work += dt * forces[body][axis] * velocity;
        }
    }
    for axis in 0..2 { close(report.sliding.friction.traction_n[axis], traction[axis]); }
    close(report.sliding.material.speed_mps, speed);
    assert!((speed - (expected[0][0] - expected[1][0]).abs()).abs() > 0.1);
    close(next.pair().contact().accumulated_distance_m(), dt * speed);
    close(report.sliding.deposited_heat_j, heat);
    close(report.kinetic_energy_change_j, kinetic);
    close(report.numerical_dissipation_j, numerical);
    close(report.external_work_j, work);
    close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + 0.4 * heat / 4.0);
    close(next.thermal().temperature_k(CELLS[1]).unwrap(), 300.0 + 0.6 * heat / 9.0);
    assert!(report.body_momentum_residual_n_s.iter().flatten().all(|x| x.abs() < 1e-8));
    assert!(report.source_motion_residual_n_s.iter().all(|x| x.abs() < 1e-8));
    close(kinetic + heat + numerical - work, 0.0);
    close(report.total_energy_residual_j, 0.0);
    assert_eq!(pair.velocities_mps(), initial);
}

#[test]
fn rotating_tangent_coordinates_preserves_speed_heat_and_exact_native_stopping() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.4, "speed", SPEED,
        [(0.0, 0.1), (2.0, 0.5)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    let velocities: [[f64; 2]; 2] = [[1.0, 0.2], [-0.2, 0.4]];
    let forces = [[0.0, 4.0], [-1.0, -2.0]];
    let run = |rotated: bool| {
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.2_f64.hypot(-0.2), InputAuthority::SyntheticFixture).unwrap();
        let initial = if rotated { velocities.map(rotate) } else { velocities };
        let drive = if rotated { forces.map(rotate) } else { forces };
        let guess = if rotated { rotate([0.8, 0.2]) } else { [0.8, 0.2] };
        SlidingPair::new(&thermal, contact, initial).unwrap().step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: drive },
            motion_controls(), guess, &[0.0; 2], other, || false).unwrap()
    };
    let ordinary = run(false);
    let rotated = run(true);
    let a = ordinary.report(); let b = rotated.report();
    for body in 0..2 {
        let expected_velocity = rotate(ordinary.pair().velocities_mps()[body]);
        let expected_displacement = rotate(a.body_displacements_m[body]);
        for axis in 0..2 {
            close(rotated.pair().velocities_mps()[body][axis], expected_velocity[axis]);
            close(b.body_displacements_m[body][axis], expected_displacement[axis]);
        }
    }
    let rotated_traction = rotate(a.sliding.friction.traction_n);
    let rotated_slip = rotate(a.sliding.friction.slip_increment_m);
    for axis in 0..2 {
        close(b.sliding.friction.traction_n[axis], rotated_traction[axis]);
        close(b.sliding.friction.slip_increment_m[axis], rotated_slip[axis]);
    }
    // s = |v_free| - 2*(.1 + .2*s) for two unit masses and dt=.1.
    let speed = (1.3_f64.hypot(0.4) - 0.2) / 1.4;
    close(a.sliding.material.speed_mps, speed);
    close(b.sliding.material.speed_mps, speed);
    close(a.sliding.friction.traction_n[0].hypot(a.sliding.friction.traction_n[1]),
        10.0 * (0.1 + 0.2 * speed));
    close(a.sliding.deposited_heat_j, (0.1 + 0.2 * speed) * speed);
    close(a.sliding.deposited_heat_j, b.sliding.deposited_heat_j);
    close(ordinary.pair().contact().accumulated_distance_m(), rotated.pair().contact().accumulated_distance_m());
    close(a.kinetic_energy_change_j, b.kinetic_energy_change_j);
    close(a.numerical_dissipation_j, b.numerical_dissipation_j);
    close(a.external_work_j, b.external_work_j);
    for report in [a, b] {
        assert!(report.body_momentum_residual_n_s.iter().flatten().all(|x| x.abs() < 1e-8));
        assert!(report.source_motion_residual_n_s.iter().all(|x| x.abs() < 1e-8));
        close(report.total_energy_residual_j, 0.0);
    }
    for cell in CELLS {
        close(ordinary.thermal().temperature_k(cell).unwrap(), rotated.thermal().temperature_k(cell).unwrap());
    }
    // The native static active set is exactly zero in both components, even
    // when neither the initial motion nor the force is tangent-axis aligned.
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.6, 0.8], [0.0; 2]]).unwrap();
    let stopped = pair.step_candidate(&thermal, 0.25,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
        motion_controls(), [0.1, 0.1], &[0.0; 2], other, || false).unwrap();
    let report = stopped.report();
    assert_eq!(report.sliding.friction.regime, fs_tribo::FrictionRegime::Sticking);
    assert_eq!(stopped.pair().velocities_mps(), [[0.3, 0.4]; 2]);
    assert_eq!(report.sliding.friction.slip_increment_m, [0.0; 2]);
    assert_eq!(report.sliding.material.speed_mps, 0.0);
    assert_eq!(report.sliding.deposited_heat_j, 0.0);
    assert_eq!(stopped.thermal().cells(), thermal.cells());
    close(report.numerical_dissipation_j, 0.25);
    close(report.total_energy_residual_j, 0.0);
}

// Synthetic density/Cp data; fixed unit-cube geometry uses the existing compiler.
fn sourced_body(density: f64, cp: f64) -> (MaterialCard, ResolvedHeatCapacityCurve, ResolvedBody) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [(heat.clone(), cp), (rho.clone(), density)] {
        let dims = key.quantity().dims();
        claims.insert_claim(PropertyClaim { key, value: PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, value), (2000.0, value)], dims,
        }, validity: ValidityDomain::unconstrained().with_quantity("T",
            QuantitySpec::dimensional(Temperature::DIMS), 200.0, 2000.0).with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic sourced sliding body".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "authored test".into(), revision: 0 }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &heat, density: &rho, temperature_axis: "T", lower: &common::point(200.0),
        upper: &common::point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    let body = resolve_body(&[RegionBinding { key: "complete body",
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        material: &card, point: &common::point(300.0),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }], 1, || false).unwrap();
    (card, caloric, body)
}

#[test]
fn density_source_swap_changes_inertia_and_heating_at_unchanged_geometry() {
    let mut outcomes = Vec::new();
    for density in [1.0, 2.0] {
        let (a, ca, ba) = sourced_body(density, 2.0);
        let (b, cb, bb) = sourced_body(1.0, 3.0);
        let cards = [a, b]; let curves = [ca, cb]; let bodies = [ba, bb]; let point = context();
        let mass = bodies.each_ref().map(|body| body.inertia().mass());
        let thermal = inventory_with_geometry(&curves, mass, [1.0; 2]);
        let original_cells = thermal.cells().clone();
        let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::from_material_bodies(&thermal, contact, bodies, [[1.0, 0.0], [0.0; 2]]).unwrap();
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
            motion_controls(), [0.8, 0.0], &[0.0; 2], other, || false).unwrap();
        let slip = 0.1 * (0.8 - 0.2 / density);
        close(next.pair().masses_kg()[0], density);
        close(next.pair().velocities_mps()[0][0], 1.0 - 0.2 / density);
        close(next.pair().velocities_mps()[1][0], 0.2);
        close(next.pair().contact().accumulated_distance_m(), slip);
        close(next.report().sliding.deposited_heat_j, 2.0 * slip);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + slip / (2.0 * density));
        close(next.report().kinetic_energy_change_j, -0.18 + 0.02 / density);
        close(next.report().numerical_dissipation_j, 0.02 + 0.02 / density);
        close(next.pair().kinetic_energy_j().unwrap() - pair.kinetic_energy_j().unwrap(),
            next.report().kinetic_energy_change_j);
        assert!(next.report().total_energy_residual_j.abs() < 1e-6);
        assert!(next.report().body_momentum_residual_n_s.iter().flatten().all(|r| r.abs() < 1e-8));
        let retained = next.pair().material_bodies().unwrap();
        assert_eq!(retained[0].regions()[0].material().card_identity(), cards[0].content_hash());
        assert!(std::ptr::eq(pair.material_bodies().unwrap(), retained));
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
        outcomes.push((next.pair().velocities_mps()[0][0], next.thermal().temperature_k(CELLS[0]).unwrap()));
    }
    assert!(outcomes[1].0 > outcomes[0].0);
    assert!(outcomes[1].1 < outcomes[0].1);
}

#[test]
fn sourced_motion_refuses_foreign_material_mass_and_volume_inventories() {
    let (a, ca, ba) = sourced_body(2.0, 2.0);
    let (b, cb, bb) = sourced_body(1.0, 3.0);
    let cards = [a, b]; let curves = [ca, cb]; let bodies = [ba, bb]; let point = context();
    let thermal = inventory_with_geometry(&curves, [2.0, 1.0], [1.0; 2]);
    let original_cells = thermal.cells().clone();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    // Equal source mass/volume cannot substitute a different original material.
    let (_, _, foreign) = sourced_body(2.0, 4.0);
    assert!(matches!(SlidingPair::from_material_bodies(&thermal, contact.clone(),
        [foreign, bodies[1].clone()], [[1.0, 0.0], [0.0; 2]]), Err(SlidingError::Invalid(_))));
    let wrong_mass = inventory_with_geometry(&curves, [3.0, 1.0], [1.0; 2]);
    assert!(matches!(SlidingPair::from_material_bodies(&wrong_mass, contact.clone(),
        bodies.clone(), [[1.0, 0.0], [0.0; 2]]), Err(SlidingError::Invalid(_))));
    let wrong_volume = inventory_with_geometry(&curves, [2.0, 1.0], [2.0, 1.0]);
    assert!(matches!(SlidingPair::from_material_bodies(&wrong_volume, contact.clone(),
        bodies, [[1.0, 0.0], [0.0; 2]]), Err(SlidingError::Invalid(_))));
    let authored = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    assert!(authored.material_bodies().is_none());
    assert_eq!(authored.masses_kg(), [2.0, 1.0]);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(authored.contact().time_s(), 0.0);
}
