//! Constitutive closure on the actual native pressure/wall/body heat timeline.
//! The synthetic affine source gives independent scalar roots, not a measured
//! non-Newtonian calibration or a continuum convergence claim.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::rheology::QuasisteadyRheologyControls;

fn iteration() -> QuasisteadyRheologyControls {
    QuasisteadyRheologyControls {
        max_iterations: 128,
        absolute_rate_tolerance_per_s: 1e-11,
        relative_rate_tolerance: 1e-10,
        relative_transport_tolerance: 1e-10,
        relaxation: 0.5,
    }
}

fn pressure_channel(sign: f64) -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Pressure {
        pressure_pa: 100000.0 + sign * 30.0,
    };
    faces[1] = FvBoundary::Pressure {
        pressure_pa: 100000.0,
    };
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-10;
    config
}

// Two transverse rows have identical u and no internal transverse shear.
// Their native wall half-cell quadrature is gamma_rms²=2u², and the exact
// discrete pressure balance is 2(a+b*sqrt(2)*|u|)u=sign*10.
fn channel_speed(a: f64, b: f64) -> f64 {
    10.0 / (a + (a * a + 20.0 * b * 2.0_f64.sqrt()).sqrt())
}

#[derive(Clone, Debug, Default, PartialEq)]
struct History {
    starts: usize,
    sources: usize,
    conductances: usize,
    intervals: Vec<f64>,
}

impl DepositionObserver for History {
    type Error = DepositionError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        *self = Self {
            starts: 1,
            ..Self::default()
        };
        Ok(())
    }

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.sources += 1;
        Ok(BTreeMap::new())
    }

    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.conductances += 1;
        Ok(BTreeMap::new())
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.intervals.push(dt);
        Ok(())
    }
}

#[test]
fn original_rate_data_closes_pressure_flow_and_heat_independently_of_seed_and_direction() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    let mut material_speeds = Vec::new();
    for b in [0.0, 300.0, 1000.0] {
        let card = rheological_card(0.0, b, 2.0);
        let liquid = law(&card);
        let caloric = caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
        let u = channel_speed(100.0, b);
        let gamma = 2.0_f64.sqrt() * u;
        for sign in [-1.0, 1.0] {
            let config = pressure_channel(sign);
            for seed in [0.0, 0.2] {
                let observer = source
                    .quasisteady_feedback(
                        &old,
                        &domain,
                        &config,
                        History::default(),
                        &gate,
                        256,
                        vec![seed; 6],
                        iteration(),
                    )
                    .unwrap()
                    .with_pressure_work_heating(controls())
                    .unwrap();
                let dt = 0.125;
                let (next, accepted, audit) = old
                    .advance_observed(&observer, dt, thermal_step(), || Ok(()))
                    .unwrap();
                let applied = accepted.applied().unwrap();
                let power = applied.viscous_power.as_ref().unwrap();
                for i in 0..=3 {
                    for j in 0..2 {
                        near(applied.flow.field.flux_x(i, j, 0), sign * u);
                    }
                }
                for (index, sample) in applied.samples.iter().enumerate() {
                    near(applied.shear_rate_per_s.as_ref().unwrap()[index], gamma);
                    near(accepted.endpoint_shear_rate_per_s().unwrap()[index], gamma);
                    near(sample.material.dynamic_viscosity_pa_s(), 100.0 + b * gamma);
                    near(
                        sample.material.thermal_conductivity_w_m_k(),
                        0.01 + 0.05 * gamma,
                    );
                    assert_eq!(
                        sample.material.material().card_identity(),
                        card.content_hash()
                    );
                    near(power.cell_heat_w[index], 10.0 * u);
                    near(next.cells()[&sample.cell].specific_h_j_kg, dt * 10.0 * u);
                }
                near(power.supplied_pressure_power_w, 60.0 * u);
                near(power.molecular_heat_w, 60.0 * u);
                near(power.heat_w, 60.0 * u);
                near(power.residual_w, 0.0);
                near(audit.external_heat_j, dt * 60.0 * u);
                near(next.energy_j() - old.energy_j(), audit.external_heat_j);
                near(
                    accepted.viscous_heating_audit().unwrap().generated_heat_j,
                    audit.external_heat_j,
                );
                let convergence = accepted.quasisteady_rheology_report().unwrap();
                assert!(
                    convergence.iterations > 1
                        && convergence.iterations <= iteration().max_iterations
                );
                assert_eq!(accepted.momentum_solves(), convergence.iterations);
                assert!(
                    convergence.maximum_relative_transport_mismatch
                        <= iteration().relative_transport_tolerance
                );
                assert_eq!(
                    accepted.inner(),
                    &History {
                        starts: 1,
                        sources: 1,
                        conductances: 1,
                        intervals: vec![dt],
                    }
                );
                assert!(observer.applied().is_none());
                assert_eq!(observer.inner(), &History::default());
                assert_eq!(observer.momentum_solves(), 0);
            }
        }
        material_speeds.push(u);
    }
    assert!(material_speeds[0] > material_speeds[1] && material_speeds[1] > material_speeds[2]);
}

#[test]
fn thermal_feedback_recloses_the_source_and_accounts_stream_enthalpy_once() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = pressure_channel(1.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    let b = 300.0;
    let card = rheological_card(-5.0, b, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let observer = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            History::default(),
            &gate,
            256,
            vec![0.0; 6],
            iteration(),
        )
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap();
    let dt = 0.5;
    let u1 = channel_speed(100.0, b);
    let h1 = dt * 10.0 * u1;
    let t1 = 300.0 + h1 / 100.0;
    let (warm, first, _) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    for cell in warm.cells().values() {
        near(cell.specific_h_j_kg, h1);
    }
    let a2 = 100.0 - 5.0 * (t1 - 300.0);
    let u2 = channel_speed(a2, b);
    let gamma2 = 2.0_f64.sqrt() * u2;
    let (next, second, audit) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let applied = second.applied().unwrap();
    near(applied.flow.field.flux_x(1, 0, 0), u2);
    assert!(u2 > u1);
    for sample in &applied.samples {
        near(sample.temperature_k, t1);
        near(sample.material.dynamic_viscosity_pa_s(), a2 + b * gamma2);
        near(
            sample.material.thermal_conductivity_w_m_k(),
            0.01 + 0.05 * gamma2,
        );
    }
    let pressure = 60.0 * u2;
    let stream = -2.0 * u2 * h1;
    near(second.audit().external_enthalpy_j, dt * stream);
    near(
        second
            .viscous_heating_audit()
            .unwrap()
            .supplied_pressure_work_j,
        dt * pressure,
    );
    near(
        second.viscous_heating_audit().unwrap().generated_heat_j,
        dt * pressure,
    );
    near(audit.external_heat_j, dt * (pressure + stream));
    near(next.energy_j() - warm.energy_j(), dt * (pressure + stream));
    assert_eq!(second.inner().sources, 1);
    assert_eq!(second.inner().intervals, vec![dt]);
    let (replayed, replay_owner, replay_audit) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(replayed.cells(), next.cells());
    assert_eq!(replay_owner.applied(), second.applied());
    assert_eq!(replay_audit, audit);
    assert_eq!(first.inner().intervals, vec![dt]);
}

#[test]
fn converged_rate_conductivity_enters_the_original_face_heat_exchange() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = pressure_channel(1.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    let b = 300.0;
    let card = rheological_card(0.0, b, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let seed = initial(&domain, 1.0, &caloric);
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = if cell.id.0 < 3 { 100.0 } else { 0.0 };
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 100,
            a: CellId(1),
            b: CellId(4),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let old = DepositionSimulation::new(
        events,
        &[(MAT, caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 64,
            max_events: 7,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_step(), || Ok(()))
    .unwrap()
    .0;
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let observer = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            (),
            &gate,
            256,
            vec![0.0; 6],
            iteration(),
        )
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap()
        .with_fluid_conduction(&old, &[100])
        .unwrap();
    let dt = 0.125;
    let u = channel_speed(100.0, b);
    let k = 0.01 + 0.05 * 2.0_f64.sqrt() * u;
    let (next, accepted, audit) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    // Half-cell resistances .5/k + .5/k give G=k. Interior x neighbours
    // have identical initial h, so these two cells have zero advective change.
    near(
        next.cells()[&CellId(1)].specific_h_j_kg,
        100.0 + dt * (10.0 * u - k),
    );
    near(
        next.cells()[&CellId(4)].specific_h_j_kg,
        dt * (10.0 * u + k),
    );
    near(
        next.energy_j() - old.energy_j(),
        dt * (60.0 * u - 100.0 * u),
    );
    near(audit.external_heat_j, dt * (60.0 * u - 100.0 * u));
    near(accepted.audit().external_enthalpy_j, -dt * 100.0 * u);
    assert_eq!(next.links()[&100].conductance_w_k, 999.0);
}

#[test]
fn body_forced_closed_material_rheology_uses_the_same_heat_owner() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let b = 100.0;
    let card = rheological_card(0.0, b, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let observer = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            (),
            &gate,
            256,
            vec![0.1; 4],
            iteration(),
        )
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap();
    let dt = 0.01;
    let (next, accepted, audit) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    // Four-face circulation: D=24mu*q², W=400q, gamma_rms=sqrt(6)*q.
    let q = (100.0 / 3.0) / (100.0 + (10000.0 + (200.0 / 3.0) * b * 6.0_f64.sqrt()).sqrt());
    let applied = accepted.applied().unwrap();
    near(applied.flow.field.flux_x(1, 1, 0), q);
    for rate in applied.shear_rate_per_s.as_ref().unwrap() {
        near(*rate, 6.0_f64.sqrt() * q);
    }
    let power = applied.viscous_power.as_ref().unwrap();
    near(power.supplied_body_power_w, 400.0 * q);
    near(power.supplied_pressure_power_w, 0.0);
    near(power.heat_w, 400.0 * q);
    near(next.energy_j() - old.energy_j(), dt * 400.0 * q);
    near(audit.external_heat_j, dt * 400.0 * q);
}

#[test]
fn constant_source_preserves_native_pressure_motion_and_thermal_results() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = pressure_channel(1.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let input = input(&ids, &liquid, &caloric);
    let source = LiquidFlowRheology::try_new(input, "rate").unwrap();
    let rheological = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            (),
            &gate,
            256,
            vec![0.0; 6],
            iteration(),
        )
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap();
    let standard = input
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap();
    let (a, ah, aa) = old
        .advance_observed(&rheological, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    let (b, bh, ba) = old
        .advance_observed(&standard, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(ah.applied().unwrap().flow, bh.applied().unwrap().flow);
    assert_eq!(ah.viscous_heating_audit(), bh.viscous_heating_audit());
    assert_eq!(aa, ba);
}

#[test]
fn source_exit_iteration_budget_and_late_cancel_preserve_all_accepted_histories() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = pressure_channel(1.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    for (max_rate, max_iterations, max_solves) in [(0.02, 128, 256), (2.0, 1, 256), (2.0, 128, 1)] {
        let card = rheological_card(0.0, 300.0, max_rate);
        let liquid = law(&card);
        let caloric = caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
        let observer = source
            .quasisteady_feedback(
                &old,
                &domain,
                &config,
                History::default(),
                &gate,
                max_solves,
                vec![0.0; 6],
                QuasisteadyRheologyControls {
                    max_iterations,
                    ..iteration()
                },
            )
            .unwrap()
            .with_pressure_work_heating(controls())
            .unwrap();
        assert!(
            old.advance_observed(&observer, 0.125, thermal_step(), || Ok(()))
                .is_err()
        );
        assert_eq!(observer.momentum_solves(), 0);
        assert!(observer.applied().is_none());
        assert_eq!(observer.inner(), &History::default());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.energy_j(), 0.0);
    }

    let card = rheological_card(-5.0, 300.0, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let saved = old.cells().clone();
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let observer = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            History::default(),
            &gate,
            256,
            vec![0.0; 6],
            iteration(),
        )
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap();
    let mut calls = 0;
    let (expected, expected_owner, _) = old
        .advance_observed(&observer, 0.125, thermal_step(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut visited = 0;
    let failed = old.advance_observed(&observer, 0.125, thermal_step(), || {
        visited += 1;
        if visited == calls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(failed.is_err());
    assert_eq!(old.cells(), &saved);
    assert_eq!(observer.inner(), &History::default());
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    let (retried, retried_owner, _) = old
        .advance_observed(&observer, 0.125, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(retried.cells(), expected.cells());
    assert_eq!(retried_owner.applied(), expected_owner.applied());
    assert_eq!(retried_owner.inner(), expected_owner.inner());
}
