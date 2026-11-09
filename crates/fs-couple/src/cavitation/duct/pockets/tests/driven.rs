//! Original DG response, independent column identities and driven contact events.
use super::super::driven::DuctReservoirResistance;
use super::*;
use fs_duct::dg::PortInput;

fn bind(
    duct: DgDuct<DgFluidMaterial<CavitationState>>,
    initial: DgState,
    cavity: DuctPocket,
    drive: PortId,
) -> CavitatingDuctPockets<CavitationState> {
    CavitatingDuctPockets::new_with_driven_port(
        duct,
        initial,
        0.0,
        &[cavity],
        &[],
        &[],
        drive,
        limits(),
    )
    .unwrap()
}

#[test]
fn either_endpoint_drive_recovers_native_response_and_uniform_compression() {
    let dt = 1e-4;
    let inflow = area() * 0.001;
    for inlet_drive in [true, false] {
        let duct = native(
            1,
            DgBoundary::Port(PortId(9)),
            DgBoundary::Port(PortId(2)),
            &[],
        );
        // q_x=-Q/L and uniform p(t)=rho*c^2*Q*t/(A L) solve the
        // continuum equations exactly and lie in this native P2 space.
        let initial = duct.project_initial_state(
            |_| 0.0,
            |x| {
                if inlet_drive {
                    inflow * (1.0 - x)
                } else {
                    -inflow * x
                }
            },
        );
        let before = initial.stored_energy_j();
        let drive = PortId(if inlet_drive { 9 } else { 2 });
        let cavity = pocket(if inlet_drive { 2 } else { 9 }, 0.0);
        let propagator = duct.propagator(dt, 1e-12).unwrap();
        let mut reference = initial.clone();
        let reference_step = propagator
            .step(
                &mut reference,
                &[PortInput {
                    id: drive,
                    effort: inflow,
                }],
            )
            .unwrap();
        propagator.sync_nodal(&mut reference).unwrap();
        let mut model = bind(duct, initial, cavity, drive);
        let step = model
            .advance_with_port_flow(dt, &numerics(), inflow, &mut || false)
            .unwrap();
        let pressure = RHO * 1500.0_f64.powi(2) * inflow * dt / area();
        for (&actual, &native) in model
            .state()
            .duct()
            .pressure()
            .iter()
            .zip(reference.pressure())
        {
            near(actual, pressure, 1e-7);
            near(actual, native, 1e-7);
        }
        for (&actual, &native) in model
            .state()
            .duct()
            .volume_flow()
            .iter()
            .zip(reference.volume_flow())
        {
            near(actual, native, 1e-13);
        }
        let source = step.drive.unwrap();
        let work = 0.5 * pressure * inflow * dt;
        assert_eq!(source.port, drive);
        assert_eq!(step.pockets[0].mode, DuctCavityMode::LiquidFilled);
        assert_eq!(step.pockets[0].volume_m3, 0.0);
        near(source.end_absolute_pressure_pa, REFERENCE + pressure, 1e-7);
        near(
            source.acoustic_pressure_impulse_pa_s,
            0.5 * pressure * dt,
            1e-11,
        );
        near(source.acoustic_work_j, work, 1e-13);
        near(
            source.acoustic_work_j,
            reference_step.ledger.port_work_step_j,
            1e-13,
        );
        near(
            source.reference_pressure_transport_work_j,
            REFERENCE * inflow * dt,
            1e-15,
        );
        near(step.total_stored_energy_j - before, work, 1e-13);
        near(step.energy_balance_residual_j, 0.0, 1e-13);
        assert_eq!(
            step.duct
                .port_outputs
                .iter()
                .map(|p| p.id)
                .collect::<Vec<_>>(),
            [PortId(9), PortId(2)]
        );
    }
}

#[test]
fn driven_open_pocket_retains_exact_column_velocity_and_opposite_port_work() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    let inflow = area() * 0.01;
    let volume = 1e-7;
    let dt = 1e-3;
    let initial = duct.project_initial_state(|_| VAPOR - REFERENCE, |_| -inflow);
    let energy = initial.stored_energy_j() + (REFERENCE - VAPOR) * volume;
    let mut model = bind(duct, initial, pocket(9, volume), PortId(2));
    assert!(model.drive.as_ref().unwrap().pressure_rates[0].abs() > 0.01 * model.gram[0]);
    let step = model
        .advance_with_port_flow(dt, &numerics(), inflow, &mut || false)
        .unwrap();
    near(step.pockets[0].end_inflow_m3_s, -inflow, 2e-12);
    near(step.pockets[0].volume_m3, volume - inflow * dt, 1e-14);
    near(step.pockets[0].end_absolute_pressure_pa, VAPOR, 1e-6);
    for &q in model.state().duct().volume_flow() {
        near(q, -inflow, 2e-12);
    }
    let source = step.drive.unwrap();
    near(
        source.acoustic_work_j,
        (VAPOR - REFERENCE) * inflow * dt,
        1e-11,
    );
    near(step.duct.ledger.port_work_step_j, 0.0, 1e-11);
    near(
        step.total_stored_energy_j - energy,
        source.acoustic_work_j,
        1e-11,
    );
    near(step.energy_balance_residual_j, 0.0, 1e-11);
}

#[test]
fn reversing_the_drive_resolves_initial_rupture_and_subsequent_collapse() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    // Only the declared pocket is at saturation. Other liquid stations retain
    // a positive pressure margin during the withdrawal/return interval.
    let initial = duct.project_initial_state(|x| VAPOR - REFERENCE + 1e5 * x * (2.0 - x), |_| 0.0);
    let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
    assert_eq!(model.state().mode(), DuctCavityMode::LiquidFilled);
    let inflow = area() * 0.1;
    let opening = model
        .advance_with_port_flow(1e-6, &numerics(), -inflow, &mut || false)
        .unwrap();
    assert_eq!(opening.pockets[0].rupture_events, 1);
    assert_eq!(opening.pockets[0].last_rupture_time_s, Some(0.0));
    assert!(opening.pockets[0].volume_m3 > 100.0 * numerics().absolute_volume_tolerance_m3);
    assert_eq!(opening.pockets[0].mode, DuctCavityMode::VaporCavity);
    let closing = model
        .advance_with_port_flow(3e-6, &numerics(), inflow, &mut || false)
        .unwrap();
    assert_eq!(closing.pockets[0].collapse_events, 1);
    assert_eq!(closing.pockets[0].mode, DuctCavityMode::LiquidFilled);
    assert!(closing.pockets[0].volume_m3.abs() <= numerics().absolute_volume_tolerance_m3);
    assert!(closing.pockets[0].last_collapse_time_s.unwrap() > 1e-6);
    near(opening.energy_balance_residual_j, 0.0, 1e-10);
    near(closing.energy_balance_residual_j, 0.0, 1e-10);
}

#[test]
fn inadmissible_drive_and_interruption_preserve_the_entire_accepted_state() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    let initial = duct.zero_state();
    let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
    let before = model.state().clone();
    for invalid in [
        f64::NAN,
        area() * 1500.0 * limits().maximum_interface_mach * 1.01,
    ] {
        assert!(
            model
                .advance_with_port_flow(1e-4, &numerics(), invalid, &mut || false)
                .is_err()
        );
        assert_eq!(model.state(), &before);
    }
    let mut polls = 0;
    let error = model
        .advance_with_port_flow(1e-4, &numerics(), area() * 0.001, &mut || {
            polls += 1;
            polls >= 5
        })
        .unwrap_err();
    assert_eq!(error, DuctCavityError::Cancelled);
    assert_eq!(model.state(), &before);
}

#[test]
fn time_and_pressure_dependent_flow_recovers_an_exact_acoustic_polynomial() {
    let duration = 1e-3;
    let initial_flow = area() * 1e-4;
    let acceleration = area() * 0.1;
    let stiffness = RHO * 1500.0_f64.powi(2) / area();
    let spatial_pressure = RHO * acceleration / (2.0 * area());
    let volume_at = |t: f64| initial_flow * t + 0.5 * acceleration * t * t;
    let volume = volume_at(duration);
    // These p(x,t), q(x,t) satisfy BOTH continuum acoustic equations and the
    // sealed inlet exactly, and lie inside the original P2 discretization.
    // Feedback vanishes on that exact trajectory but must use CURRENT p.
    for feedback in [0.0, 1e-10] {
        let duct = native(
            1,
            DgBoundary::Port(PortId(9)),
            DgBoundary::Port(PortId(2)),
            &[],
        );
        let positions = duct.node_positions();
        let initial =
            duct.project_initial_state(|x| spatial_pressure * x * x, |x| -initial_flow * x);
        let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
        let law = |time: f64, pressure: f64| {
            let expected_pressure = REFERENCE + spatial_pressure + stiffness * volume_at(time);
            Ok(initial_flow + acceleration * time + feedback * (expected_pressure - pressure))
        };
        let step = model
            .advance_with_port_flow_law(duration, &numerics(), &law, &mut || false)
            .unwrap();
        let source = step.drive.unwrap();
        near(source.initial_inflow_m3_s, initial_flow, 1e-16);
        near(
            source.inflow_m3_s,
            initial_flow + acceleration * duration,
            2e-13,
        );
        near(source.transported_volume_m3, volume, 2e-16);
        near(source.average_inflow_m3_s, volume / duration, 2e-13);
        assert!((source.transported_volume_m3 - source.inflow_m3_s * duration).abs() > 1e-11);
        near(
            source.acoustic_work_j,
            spatial_pressure * volume + 0.5 * stiffness * volume * volume,
            1e-12,
        );
        near(
            source.acoustic_pressure_impulse_pa_s,
            spatial_pressure * duration
                + stiffness
                    * (0.5 * initial_flow * duration.powi(2)
                        + acceleration * duration.powi(3) / 6.0),
            1e-8,
        );
        for (i, &x) in positions.iter().enumerate() {
            near(
                model.state().duct().pressure()[i],
                spatial_pressure * x * x + stiffness * volume,
                2e-3,
            );
            near(
                model.state().duct().volume_flow()[i],
                -(initial_flow + acceleration * duration) * x,
                2e-13,
            );
        }
        assert_eq!(step.pockets[0].mode, DuctCavityMode::LiquidFilled);
        assert_eq!(source.reservoir, None);
        near(step.energy_balance_residual_j, 0.0, 1e-11);
    }
}

#[test]
fn passive_reservoir_recovers_the_distributed_decay_and_lumped_rc_limit() {
    // The exact continuum mode has p=P cosh(lambda*x/c) exp(-lambda*t),
    // q=A P/(rho*c) sinh(lambda*x/c) exp(-lambda*t), and
    // tanh(lambda*L/c)=Zc/R. R>>Zc tends to tau=R*A*L/(rho*c^2).
    let sound = 1500.0;
    let impedance = RHO * sound / area();
    let resistance = 10.0 * impedance;
    let z = (impedance / resistance).atanh();
    let decay = sound * z;
    let duration = 1.0 / decay;
    let attenuation = (-decay * duration).exp();
    let compliance = area() / (RHO * sound * sound);
    let reservoir = DuctReservoirResistance::try_new(REFERENCE, resistance).unwrap();
    for amplitude in [-1000.0, 1000.0] {
        let duct = native(
            4,
            DgBoundary::Port(PortId(9)),
            DgBoundary::Port(PortId(2)),
            &[],
        );
        let positions = duct.node_positions();
        let initial = duct.project_initial_state(
            |x| amplitude * (z * x).cosh(),
            |x| amplitude / impedance * (z * x).sinh(),
        );
        let before = initial.stored_energy_j();
        let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
        let mut controls = numerics();
        controls.absolute_energy_tolerance_j = 1e-13;
        controls.absolute_volume_tolerance_m3 = 1e-16;
        let step = model
            .advance_with_reservoir(duration, &controls, reservoir, &mut || false)
            .unwrap();
        let source = step.drive.unwrap();
        let balance = source.reservoir.unwrap();
        let initial_flow = -amplitude * z.cosh() / resistance;
        let volume = initial_flow * (1.0 - attenuation) / decay;
        let loss = resistance * initial_flow.powi(2) * (1.0 - attenuation.powi(2)) / (2.0 * decay);
        near(source.transported_volume_m3, volume, volume.abs() * 2e-4);
        near(
            source.inflow_m3_s,
            initial_flow * attenuation,
            initial_flow.abs() * 2e-4,
        );
        near(balance.dissipated_energy_j, loss, loss * 2e-4);
        near(
            balance.supplied_absolute_work_j,
            REFERENCE * source.transported_volume_m3,
            1e-16,
        );
        assert_eq!(balance.supplied_reduced_work_j, 0.0);
        assert!(balance.dissipated_energy_j > 0.0);
        near(source.acoustic_work_j, -balance.dissipated_energy_j, 1e-13);
        near(balance.energy_balance_residual_j, 0.0, 1e-13);
        near(
            step.total_stored_energy_j - before + balance.dissipated_energy_j,
            0.0,
            2e-12,
        );
        for (i, &x) in positions.iter().enumerate() {
            near(
                model.state().duct().pressure()[i],
                amplitude * (z * x).cosh() * attenuation,
                0.02,
            );
        }
        let rc_pressure = amplitude * (-duration / (resistance * compliance)).exp();
        near(
            model.state().duct().pressure()[0],
            rc_pressure,
            0.005 * rc_pressure.abs(),
        );
        assert_eq!(step.pockets[0].mode, DuctCavityMode::LiquidFilled);
    }
}

#[test]
fn reservoir_withdrawal_and_return_resolve_rupture_and_collapse_on_the_same_clock() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    let initial = duct.project_initial_state(|x| VAPOR - REFERENCE + 1e5 * x * (2.0 - x), |_| 0.0);
    let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
    let resistance = 72_000.0 / (area() * 0.1);
    let withdraw = DuctReservoirResistance::try_new(30_000.0, resistance).unwrap();
    let supply = DuctReservoirResistance::try_new(174_000.0, resistance).unwrap();
    let opening = model
        .advance_with_reservoir(1e-6, &numerics(), withdraw, &mut || false)
        .unwrap();
    assert_eq!(opening.pockets[0].rupture_events, 1);
    assert_eq!(opening.pockets[0].last_rupture_time_s, Some(0.0));
    assert_eq!(opening.pockets[0].mode, DuctCavityMode::VaporCavity);
    assert!(opening.drive.unwrap().transported_volume_m3 < 0.0);
    let closing = model
        .advance_with_reservoir(3e-6, &numerics(), supply, &mut || false)
        .unwrap();
    assert_eq!(closing.pockets[0].collapse_events, 1);
    assert_eq!(closing.pockets[0].mode, DuctCavityMode::LiquidFilled);
    assert!(closing.pockets[0].last_collapse_time_s.unwrap() > 1e-6);
    assert!(closing.drive.unwrap().transported_volume_m3 > 0.0);
    for step in [&opening, &closing] {
        let source = step.drive.unwrap();
        let reservoir = source.reservoir.unwrap();
        assert!(reservoir.dissipated_energy_j > 0.0);
        near(reservoir.energy_balance_residual_j, 0.0, 1e-11);
        near(step.energy_balance_residual_j, 0.0, 1e-10);
    }
}

#[test]
fn hydraulic_callback_failure_invalid_reservoir_and_cancellation_are_atomic() {
    let duct = native(
        1,
        DgBoundary::Port(PortId(9)),
        DgBoundary::Port(PortId(2)),
        &[],
    );
    let initial = duct.zero_state();
    let mut model = bind(duct, initial, pocket(9, 0.0), PortId(2));
    let before = model.state().clone();
    let expected = DuctCavityError::OutsideDomain("caller hydraulic domain");
    let law = |time: f64, _: f64| {
        if time > 0.0 {
            Err(expected.clone())
        } else {
            Ok(0.0)
        }
    };
    assert_eq!(
        model
            .advance_with_port_flow_law(1e-4, &numerics(), &law, &mut || false)
            .unwrap_err(),
        expected
    );
    assert_eq!(model.state(), &before);
    let invalid = DuctReservoirResistance::try_new(0.5 * VAPOR, 1e10).unwrap();
    assert!(
        model
            .advance_with_reservoir(1e-4, &numerics(), invalid, &mut || false)
            .is_err()
    );
    assert_eq!(model.state(), &before);
    let excessive = DuctReservoirResistance::try_new(
        REFERENCE + 1.1 * limits().maximum_liquid_strain * RHO * 1500.0_f64.powi(2),
        1e10,
    )
    .unwrap();
    assert!(
        model
            .advance_with_reservoir(1e-4, &numerics(), excessive, &mut || false)
            .is_err()
    );
    assert_eq!(model.state(), &before);
    let reservoir = DuctReservoirResistance::try_new(REFERENCE + 1000.0, 1e10).unwrap();
    let mut polls = 0;
    let error = model
        .advance_with_reservoir(1e-4, &numerics(), reservoir, &mut || {
            polls += 1;
            polls >= 5
        })
        .unwrap_err();
    assert_eq!(error, DuctCavityError::Cancelled);
    assert_eq!(model.state(), &before);
    near(
        reservoir.inflow_m3_s(-500.0).unwrap(),
        (REFERENCE + 1500.0) / 1e10,
        1e-20,
    );
    assert!(DuctReservoirResistance::try_new(REFERENCE, 0.0).is_err());
    assert!(DuctReservoirResistance::try_new(f64::NAN, 1e10).is_err());
}
