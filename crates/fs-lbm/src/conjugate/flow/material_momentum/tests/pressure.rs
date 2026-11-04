//! Reservoir traction, actual density and one accepted material transport.

use super::*;

#[test]
fn pressure_drives_unequal_density_mass_enthalpy_and_momentum_once() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let mass = [2.0, 4.0];
    let enthalpy = [10.0, 50.0];
    let zero = [[0.0; 3]; 2];
    let sources = [MaterialVolumeSource::None; 2];
    let source_velocity = [None; 2];
    let mut inlet = [None; 6];
    inlet[0] = Some(AdvectedMaterial {
        density_kg_m3: 8.0,
        specific_enthalpy_j_kg: 100.0,
    });
    let mut inlet_velocity = [None; 6];
    inlet_velocity[0] = Some([6.0, 0.0, 0.0]);
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &sources,
            boundary_inflow: inlet,
        },
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &zero,
        source_velocity_m_s: &source_velocity,
        boundary_velocity_m_s: inlet_velocity,
    };
    for reference in [0.0, 1e12] {
        let x = [Some(reference + 12.0), None, Some(reference)];
        let transverse = [None; 4];
        let profile = StaggeredBoundaryPressure {
            x_pa: &x,
            y_pa: &transverse,
            z_pa: &transverse,
        };
        let result = field
            .advance_material_momentum_with_pressure_profile(
                &domain,
                open_x(),
                profile,
                input,
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Boundary/interior dual masses are [1,3,2]. Summing their momentum
        // equations gives u = dt*(p_left-p_right)/sum(M) = 1, independently
        // of the projection. The physical cell impulse is [2,4], exactly once.
        for face in 0..3 {
            close(result.transport.projection.field.flux_x(face, 0, 0), 1.0);
        }
        close(result.force_projection.pressure_impulse_kg_m_s[0], 6.0);
        close(result.feedback.cell_momentum_change_kg_m_s[0], 6.0);
        close(result.forcing.applied_work_j, 0.0);
        close(
            result
                .transport
                .projection
                .report
                .open_boundary_pressure_work_j,
            -6.0,
        );
        close(result.force_projection.kinetic_energy_change_j, 3.0);
        close(result.force_projection.energy_defect_j, 0.0);
        // Backward-Euler density: 1.5*M0=2+4; 1.5*M1=4+0.5*M0.
        // Thus accepted route masses are inlet=4, interior=2, outlet=2.
        assert_eq!(result.transport.report.mass_routes, 3);
        for &value in &result.transport.mass_kg {
            close(value, 4.0);
        }
        close(result.transport.report.incoming_mass_kg, 4.0);
        close(result.transport.report.outgoing_mass_kg, 2.0);
        // (4+2)h0=2*10+4*100; (4+2)h1=4*50+2*h0.
        close(result.transport.specific_enthalpy_j_kg[0], 70.0);
        close(result.transport.specific_enthalpy_j_kg[1], 170.0 / 3.0);
        close(result.transport.report.enthalpy_change_j, 860.0 / 3.0);
        close(result.transport.report.carried_enthalpy_j, 860.0 / 3.0);
        // The same mass routes carry kicked P=[2,4] and inlet velocity=6:
        // 6*u0=2+4*6; 6*u1=4+2*u0. No extra pressure impulse is inserted.
        close(result.momentum.momentum_kg_m_s[0][0], 52.0 / 3.0);
        close(result.momentum.momentum_kg_m_s[1][0], 76.0 / 9.0);
        close(
            result.momentum.report.carried_momentum_kg_m_s[0],
            178.0 / 9.0,
        );
        close(result.momentum.report.momentum_after_kg_m_s[0], 232.0 / 9.0);
        assert_eq!(result.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
    }
    assert_eq!(field, original);
    assert_eq!(mass, [2.0, 4.0]);
    assert_eq!(enthalpy, [10.0, 50.0]);
}

#[test]
fn hydrostatic_reservoirs_balance_actual_unequal_density_without_transport() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 4.0];
    let enthalpy = [10.0, 50.0];
    let zero = [[0.0; 3]; 2];
    let sources = [MaterialVolumeSource::None; 2];
    let source_velocity = [None; 2];
    let mut inlet = [None; 6];
    for f in 0..2 {
        inlet[f] = Some(AdvectedMaterial {
            density_kg_m3: mass[f],
            specific_enthalpy_j_kg: enthalpy[f],
        });
    }
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &sources,
            boundary_inflow: inlet,
        },
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [2.0, 0.0, 0.0],
        force_density_n_m3: &zero,
        source_velocity_m_s: &source_velocity,
        boundary_velocity_m_s: [Some([0.0; 3]); 6],
    };
    // dp/dx=rho*g: total pressure rise is (2+4)*2=12 Pa. The two
    // cell-centre pressures are 2 and 8 Pa; no reference density is used.
    let x = [Some(1e12), None, Some(1e12 + 12.0)];
    let transverse = [None; 4];
    let profile = StaggeredBoundaryPressure {
        x_pa: &x,
        y_pa: &transverse,
        z_pa: &transverse,
    };
    let result = field
        .advance_material_momentum_with_pressure_profile(
            &domain,
            open_x(),
            profile,
            input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for c in 0..2 {
        close(result.transport.mass_kg[c], mass[c]);
        close(result.transport.specific_enthalpy_j_kg[c], enthalpy[c]);
        for &velocity in &result.momentum.velocity_m_s[c] {
            close(velocity, 0.0);
        }
    }
    close(result.forcing.applied_impulse_kg_m_s[0], 6.0);
    close(result.force_projection.pressure_impulse_kg_m_s[0], -6.0);
    close(result.feedback.cell_momentum_change_kg_m_s[0], 0.0);
    close(result.transport.report.enthalpy_change_j, 0.0);
    close(result.force_projection.energy_defect_j, 0.0);

    let empty_x = [None; 3];
    let empty = StaggeredBoundaryPressure {
        x_pa: &empty_x,
        ..profile
    };
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        field
            .advance_material_momentum(&domain, open_x(), input, controls(), &gate)
            .unwrap(),
        field
            .advance_material_momentum_with_pressure_profile(
                &domain,
                open_x(),
                empty,
                input,
                controls(),
                &gate,
            )
            .unwrap(),
    );

    let invalid_x = [Some(0.0), Some(0.0), Some(12.0)];
    assert!(
        field
            .advance_material_momentum_with_pressure_profile(
                &domain,
                open_x(),
                StaggeredBoundaryPressure {
                    x_pa: &invalid_x,
                    ..profile
                },
                input,
                controls(),
                &gate,
            )
            .is_err()
    );
    let mut late = controls();
    late.momentum.max_iterations = 0;
    assert!(
        field
            .advance_material_momentum_with_pressure_profile(
                &domain,
                open_x(),
                profile,
                input,
                late,
                &gate,
            )
            .is_err()
    );
    assert_eq!(
        field
            .advance_material_momentum_with_pressure_profile(
                &domain,
                open_x(),
                profile,
                input,
                controls(),
                &gate,
            )
            .unwrap(),
        result
    );
}
