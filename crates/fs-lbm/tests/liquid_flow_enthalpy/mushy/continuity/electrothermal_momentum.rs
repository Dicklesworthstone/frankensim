//! G3/G4: current, phase, native inertia and heat share the accepted timeline.
use super::*;

#[test]
fn current_driven_melting_preserves_native_momentum_and_its_separate_energy_accounts() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let electric = setup(1.0);
    let links: Vec<_> = (0..12).collect();
    let mut speeds = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let old = initial(&domain, &curve);
        let flow = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            FlowEnthalpyInput {
                cell_ids: &flow_ids,
                material_id: MAT,
                caloric: &curve,
                inlet_specific_h_j_kg: [None; 6],
                max_cells: 9,
                max_connections: 64,
            },
            &config,
            CycleHeat,
            &gate,
            1000,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_fluid_conduction(&old, &links)
        .unwrap()
        .with_unsteady_viscous_heating(controls())
        .unwrap();
        let (solid, solid_flow, solid_current) = advance_joule_enthalpy_with_material_flow(
            &old,
            0.1,
            step(),
            current_input(&domain, &electric, &ids),
            &flow,
            resistivity,
        )
        .unwrap();
        assert_eq!(solid_flow.unsteady_time_s(), Some(solid.time_s()));
        assert_eq!(solid_flow.momentum_solves(), 0);
        assert!(solid_current.electrode_energy_j > 0.0);
        for sample in solid_flow.endpoint() {
            assert!(sample.transport.mobile_transport().is_none());
            assert!(
                sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
            );
            assert!(
                sample
                    .material
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .is_none()
            );
        }
        near(
            solid.energy_j() - old.energy_j(),
            90.0 + solid_current.electrode_energy_j,
        );

        let (liquid, liquid_flow, melting_current) = advance_joule_enthalpy_with_material_flow(
            &solid,
            1.5,
            step(),
            current_input(&domain, &electric, &ids),
            &solid_flow,
            resistivity,
        )
        .unwrap();
        assert_eq!(liquid_flow.unsteady_time_s(), Some(liquid.time_s()));
        assert!(
            liquid_flow
                .endpoint()
                .iter()
                .all(|s| s.material.phase_state().liquid_mass_fraction() == 1.0)
        );
        assert!(
            liquid_flow
                .unsteady_flow()
                .unwrap()
                .velocity_m_s
                .iter()
                .flatten()
                .all(|&v| v == 0.0)
        );
        near(
            liquid.energy_j() - old.energy_j(),
            1350.0 + solid_current.electrode_energy_j + melting_current.electrode_energy_j,
        );

        // The actuator starts only after melting, on the original native owner.
        let force = vec![
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0; 3],
            [0.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        let driven = liquid_flow
            .with_unsteady_body_force_heating(force, controls())
            .unwrap();
        // The actual binary interval from 1.525 to 1.55 is slightly longer
        // than 0.025. Admit each intended interval whole instead of appending
        // a roundoff-sized native inertia step after a nominal 25 ms step.
        let motion_step = StepConfig {
            max_dt_s: 0.03,
            ..step()
        };
        let (moving, moving_flow, moving_current) = advance_joule_enthalpy_with_material_flow(
            &liquid,
            1.525,
            motion_step,
            current_input(&domain, &electric, &ids),
            &driven,
            resistivity,
        )
        .unwrap();
        let work = moving_flow.unsteady_heating_audit().unwrap();
        assert!(work.final_kinetic_energy_j > 0.0);
        assert!(work.generated_heat_j > 0.0);
        assert!(work.backward_euler_loss_j > 0.0);
        near(work.initial_kinetic_energy_j, 0.0);
        near(
            moving.energy_j() - liquid.energy_j(),
            moving_current.electrode_energy_j + work.generated_heat_j,
        );
        near(
            moving.energy_j() - liquid.energy_j()
                + work.final_kinetic_energy_j
                + work.backward_euler_loss_j,
            moving_current.electrode_energy_j
                + work.supplied_body_work_j
                + work.supplied_wall_work_j,
        );
        let speed = moving_flow
            .unsteady_flow()
            .unwrap()
            .velocity_m_s
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        assert!(speed > 0.0);
        speeds.push(speed);

        let saved_motion = moving_flow.unsteady_flow().unwrap().clone();
        let coast = moving_flow
            .with_unsteady_body_force_heating(vec![[0.0; 3]; 9], controls())
            .unwrap();
        assert_eq!(coast.unsteady_flow(), Some(&saved_motion));
        let (cooled, coast_flow, coast_current) = advance_joule_enthalpy_with_material_flow(
            &moving,
            1.55,
            motion_step,
            current_input(&domain, &electric, &ids),
            &coast,
            resistivity,
        )
        .unwrap();
        let decay = coast_flow.unsteady_heating_audit().unwrap();
        near(decay.initial_kinetic_energy_j, work.final_kinetic_energy_j);
        assert!(decay.final_kinetic_energy_j > 0.0);
        assert!(decay.final_kinetic_energy_j < decay.initial_kinetic_energy_j);
        near(decay.supplied_body_work_j, 0.0);
        near(decay.supplied_wall_work_j, 0.0);
        near(
            decay.initial_kinetic_energy_j - decay.final_kinetic_energy_j,
            decay.generated_heat_j + decay.backward_euler_loss_j,
        );
        near(
            cooled.energy_j() - moving.energy_j(),
            coast_current.electrode_energy_j + decay.generated_heat_j,
        );
        near(cooled.mass_kg(), old.mass_kg());
        assert_eq!(coast_flow.unsteady_time_s(), Some(cooled.time_s()));
        assert_eq!(coast_flow.endpoint().len(), 9);
        for sample in coast_flow.endpoint() {
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
            assert_eq!(
                sample.specific_enthalpy_j_kg.to_bits(),
                cooled.cells()[&sample.cell].specific_h_j_kg.to_bits()
            );
        }

        // Current's endpoint refusal comes after the candidate native step.
        // Its rejection must preserve motion, heat and source history together.
        let saved_cells = cooled.cells().clone();
        let saved_flow = coast_flow.unsteady_flow().unwrap().clone();
        let saved_heat = coast_flow.unsteady_heating_audit();
        let mut starved = current_input(&domain, &electric, &ids);
        starved.max_current_solves = 1;
        assert!(matches!(
            advance_joule_enthalpy_with_material_flow(
                &cooled,
                1.575,
                motion_step,
                starved,
                &coast_flow,
                resistivity,
            ),
            Err(MaterialFlowElectricalError::Electrical(
                JouleEnthalpyError::CurrentBudget
            ))
        ));
        assert_eq!(cooled.cells(), &saved_cells);
        assert_eq!(coast_flow.unsteady_flow(), Some(&saved_flow));
        assert_eq!(coast_flow.unsteady_heating_audit(), saved_heat);
        assert_eq!(coast_flow.unsteady_time_s(), Some(1.55));
        assert_eq!(old.time_s(), 0.0);
        assert!(flow.applied().is_none());
    }
    // Native inertia stays in the equation: doubling mu reduces the speed,
    // rather than imposing a quasisteady inverse-mu rescaling of old motion.
    assert!(speeds[1] < speeds[0]);
}
