//! Original DG response, independent column identities and driven contact events.
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
