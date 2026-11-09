//! Both physical heat owners feed actual spatial conduction and phase density.

use super::*;
use crate::conjugate::{FaceContactResistance, FvBoundary, MaterialViscosityControls, Voxel};

fn molecular_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 16,
        max_sweeps: 400,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-9,
        maximum_angular_impulse_defect_kg_m2_s: 1e-9,
        maximum_energy_defect_j: 1e-7,
    }
}

fn limits(cells: usize, faces: usize) -> ConductiveResistiveThermalVolumeMomentumControls {
    let mut limits = controls();
    limits
        .resistance
        .thermal
        .volume
        .thermal
        .volume
        .transport
        .max_cells = cells;
    limits.resistance.thermal.mechanics.reconstruction.max_cells = cells;
    limits.resistance.thermal.mechanics.forcing.max_cells = cells;
    limits.max_conductive_faces = faces;
    limits
}

fn prepare<'m>(
    card: &'m MaterialCard,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
    context: &QueryPoint,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'m> {
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        &PropertyKey::with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
        ),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

/// Fourier/contact heat and both dissipations change v(h); the finite wall
/// starts from its original inventory in EVERY rejected Darcy heat trial.
#[test]
fn finite_wall_and_both_physical_heats_drive_unequal_density_volume() {
    let (card, curve) = fixture_with_densities([4.0, 2.0]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let context = context();
    let conductivity = prepare(&card, &curve, &context, [40.0, 90.0]);
    let wall_conductivity = prepare(&card, &curve, &context, [110.0, 130.0]);
    let mut domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    domain.set(0, 1, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let initial = curve.state_at_specific_enthalpy(50.0).unwrap();
    let phases = [Some(initial), None];
    let mass = [initial.bulk_density_kg_m3(), 0.0];
    let momentum = [[0.1 * mass[0], 0.0, 0.0], [0.0; 3]];
    let heat = [mass[0] * 10.0, 0.0];
    let force = [[0.0; 3]; 2];
    let mut input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    let mut faces = [FlowFace::Wall; 6];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    for side in 0..2 {
        faces[side] = FlowFace::Free;
        boundaries[side] = FvBoundary::Outlet;
        boundaries[2 + side] = FvBoundary::wall();
        input.volume.boundary_inflow[side] = Some(initial);
        input.boundary_velocity_m_s[side] = Some([0.0; 3]);
    }
    let px = [Some(0.0), Some(0.0), None, None];
    let py = [None; 3];
    let pz = [None; 4];
    let reservoir = Some((
        faces,
        StaggeredBoundaryPressure {
            x_pa: &px,
            y_pa: &py,
            z_pa: &pz,
        },
    ));
    let solids = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall_conductivity,
        mass_kg: 2.0,
        phase: curve.state_at_specific_enthalpy(116.5).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contact_resistance = 16.0 / 15.0;
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: contact_resistance,
    })];
    let limits = limits(2, 1);
    let solve = ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let run = |workspace: &mut ImplicitThermalVolumeWorkspace,
               solve,
               gate: &CancelGate,
               cancel_endpoint: bool| {
        let mut queries = 0;
        field.advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths(
            &domain, reservoir, &conductivity, ThermalVolumeBoundaries::default(),
            &solids, &contacts, &law, input, boundaries, Some(&[0.0; 2]), limits,
            molecular_controls(), solve, workspace, gate,
            |_, phase| {
                queries += 1;
                if cancel_endpoint && queries == 2 { gate.request(); }
                source(&card, &law, phase)
            },
        )
    };
    let result = run(&mut workspace, solve, &CancelGate::new_clock_free(), false).unwrap();
    let combined = &result.evolution;
    let native = &combined.evolution.evolution;
    let volume = &native.thermal.volume;
    assert!(combined.viscosity.molecular_heat_j > 0.0);
    assert!(native.resistance.resistance_heat_j > 0.0);
    let dissipated = combined.viscosity.molecular_heat_j + native.resistance.resistance_heat_j;
    let wall_heat = result.conductive_enthalpy_increment_j[0];
    assert!(wall_heat > 0.0);
    let deposited = heat[0] + dissipated + wall_heat;
    close(volume.report.requested_enthalpy_increment_j, deposited);
    // The latent chart has v(h)=1/4+h/400. All three PHYSICAL heat sources
    // drive the native volume equation, with no local material mass source.
    close(volume.volume_source_m3_s[0], deposited / (400.0 * 0.25));
    let endpoint = volume.phase_states[0].unwrap();
    assert!(endpoint.bulk_density_kg_m3() < initial.bulk_density_kg_m3());
    close(
        volume.transport.mass_kg[0] / endpoint.bulk_density_kg_m3(),
        1.0,
    );
    assert_eq!(volume.phase_states[1], None);
    assert_eq!(volume.transport.mass_kg[1], 0.0);
    assert_eq!(native.momentum.momentum_kg_m_s[1], [0.0; 3]);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);

    // Independently reconstruct the endpoint Fourier law from the two real
    // half cells plus contact R, using PRE-advection fluid caloric state.
    let h_star = 50.0 + deposited / mass[0];
    let k_fluid = 0.125 + 0.25 * h_star / 100.0;
    let wall = &result.solid_endpoints[0];
    let conductance = 1.0 / (0.5 / k_fluid + contact_resistance + 0.5 / 0.375);
    close(
        wall_heat,
        0.25 * conductance * (wall.phase.temperature_k() - 500.0),
    );
    close(wall.mass_kg, 2.0);
    close(
        wall.mass_kg * (wall.phase.specific_enthalpy_j_kg() - 116.5),
        -wall_heat,
    );
    close(result.report.solid_enthalpy_change_j, -wall_heat);
    close(result.report.solid_to_fluid_heat_j, wall_heat);
    close(
        volume.transport.mass_kg[0] * endpoint.specific_enthalpy_j_kg() - mass[0] * 50.0
            + result.report.solid_enthalpy_change_j,
        heat[0] + dissipated + volume.transport.report.carried_enthalpy_j,
    );
    close(
        volume.transport.mass_kg[0] - mass[0],
        volume.transport.report.incoming_mass_kg - volume.transport.report.outgoing_mass_kg,
    );
    assert!(native.coupling.coupling_iterations > 2);
    assert_eq!(result.work.heat_trials, native.coupling.coupling_iterations);
    assert!(result.work.residual_evaluations > result.report.conduction.solve.residual_evaluations);
    assert_eq!(
        result.report.conduction.absolute_enthalpy_defect_j,
        combined.absolute_enthalpy_defect_j
    );
    assert!(combined.absolute_energy_defect_j <= molecular_controls().maximum_energy_defect_j);
    close(combined.viscosity.energy_defect_j, 0.0);
    close(native.resistance.energy_defect_j, 0.0);

    // Fail after a real heat solve, then exhaust aggregate work. Neither may
    // publish a viscous kick or warm-start physical wall history on retry.
    assert!(matches!(
        run(&mut workspace, solve, &CancelGate::new_clock_free(), true),
        Err(ChtError::Cancelled)
    ));
    let exhausted = ImplicitThermalVolumeSolveControls {
        maximum_residual_evaluations: result.work.residual_evaluations - 1,
        ..solve
    };
    assert!(
        run(
            &mut workspace,
            exhausted,
            &CancelGate::new_clock_free(),
            false
        )
        .is_err()
    );
    assert_eq!(
        run(&mut workspace, solve, &CancelGate::new_clock_free(), false).unwrap(),
        result
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(phases, [Some(initial), None]);
    assert_eq!(mass, [initial.bulk_density_kg_m3(), 0.0]);
    assert_eq!(momentum[0][0], 0.1 * mass[0]);
    close(solids[0].phase.specific_enthalpy_j_kg(), 116.5);
}

/// Zero Fourier exchange must preserve the existing composed physical step.
#[test]
fn isothermal_conduction_preserves_the_original_combined_step() {
    let (card, curve) = fixture();
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
    let context = context();
    let conductivity = prepare(&card, &curve, &context, [40.0, 90.0]);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 4];
    let phases = [Some(curve.state_at_specific_enthalpy(50.0).unwrap()); 4];
    let momentum = [
        [0.04, -0.04, 0.0],
        [0.04, 0.04, 0.0],
        [-0.04, -0.04, 0.0],
        [-0.04, 0.04, 0.0],
    ];
    let input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &[0.0; 4],
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 4],
        boundary_velocity_m_s: [None; 6],
    };
    let limits = limits(4, 4);
    let baseline = field
        .advance_viscous_phase_resistance_thermal_volume_momentum(
            &domain,
            None,
            &law,
            input,
            [FvBoundary::wall(); 6],
            None,
            limits.resistance,
            molecular_controls(),
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(4).unwrap();
    let result = field.advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths(
        &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[], &law,
        input, [FvBoundary::wall(); 6], None, limits, molecular_controls(),
        ImplicitThermalVolumeSolveControls::default(), &mut workspace,
        &CancelGate::new_clock_free(), |_, phase| source(&card, &law, phase),
    ).unwrap();
    assert_eq!(result.conductive_enthalpy_increment_j, vec![0.0; 4]);
    assert_eq!(result.evolution.evolution, baseline.evolution);
    assert_eq!(result.evolution.viscosity, baseline.viscosity);
    assert_eq!(
        result.evolution.viscosity_sources,
        baseline.viscosity_sources
    );
    assert!(result.solid_endpoints.is_empty());
    assert!(result.evolution.absolute_enthalpy_defect_j >= baseline.absolute_enthalpy_defect_j);
}
