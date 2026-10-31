//! G1/G3/G4: an independently reduced circulation senses signed cross-axis drag.
//! The source curves are manufactured fixtures, not measured porous materials.
use super::*;
use fs_lbm::conjugate::FlowField;

fn axes(plane: [usize; 3], sign: f64) -> [[f64; 3]; 3] {
    let [a, b, c] = plane;
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let mut axes = [[0.0; 3]; 3];
    axes[0][a] = h;
    axes[0][b] = sign * h;
    axes[1][a] = -sign * h;
    axes[1][b] = h;
    axes[2][c] = 1.0;
    axes
}

fn layers<'a>(
    source: &'a MaterialCard,
    domain: &VoxelDomain,
    plane: [usize; 3],
    orientation: f64,
) -> Vec<Option<PermeabilityTemperatureLaw<'a>>> {
    (0..domain.cell_count())
        .map(|cell| {
            let p = domain.coords(cell);
            let sign = if p[plane[0]] == p[plane[1]] {
                -orientation
            } else {
                orientation
            };
            Some(permeability(source, axes(plane, sign)))
        })
        .collect()
}

fn flux(field: &FlowField, axis: usize, p: [usize; 3]) -> f64 {
    match axis {
        0 => field.flux_x(p[0], p[1], p[2]),
        1 => field.flux_y(p[0], p[1], p[2]),
        _ => field.flux_z(p[0], p[1], p[2]),
    }
}

#[test]
fn rotated_material_changes_cross_axis_motion_and_heat_in_every_coordinate_plane() {
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let skeleton = skeleton([0.5, 1.0, 2.0], false, 810.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    for plane in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
        let [a, b, c] = plane;
        let mut n = [2; 3];
        n[c] = 1;
        let domain = VoxelDomain::new(n[0], n[1], n[2], 1.0).unwrap();
        let old = initial(&domain, 1.0, &caloric);
        let mut config = driven_box();
        config.faces = [FvBoundary::wall(); 6];
        config.faces[2 * c] = FvBoundary::Symmetry;
        config.faces[2 * c + 1] = FvBoundary::Symmetry;
        let mut wall = [0.0; 3];
        wall[a] = 1.0;
        config.faces[2 * b + 1] = FvBoundary::Wall { velocity: wall };
        let mut speeds = Vec::new();
        for orientation in [-1.0, 1.0] {
            let observer = input(&ids, &liquid, &caloric)
                .feedback(&old, &domain, &config, (), &gate, 10)
                .unwrap()
                .with_porous_materials(layers(&skeleton, &domain, plane, orientation))
                .unwrap()
                .with_viscous_heating(controls())
                .unwrap();
            let (next, accepted, thermal) = old
                .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
                .unwrap();
            let epoch = accepted.applied().unwrap();
            let power = epoch.viscous_power.as_ref().unwrap();
            // Four incompressible internal faces leave one amplitude u. The
            // 45-degree layers have B_aa=B_bb=150, |B_ab|=50. Their alternating
            // signs give D_Darcy=2*(300+orientation*50)*u²; a uniform B_ab
            // would cancel and could not detect a missing off-diagonal action.
            let resistance = 300.0 + orientation * 50.0;
            let u = 100.0 / (1200.0 + resistance);
            let mut p = [0; 3];
            p[a] = 1;
            p[b] = 1;
            near(flux(&epoch.flow.field, a, p), u);
            p[b] = 0;
            near(flux(&epoch.flow.field, a, p), -u);
            p[a] = 0;
            p[b] = 1;
            near(flux(&epoch.flow.field, b, p), u);
            p[a] = 1;
            near(flux(&epoch.flow.field, b, p), -u);
            near(power.darcy_heat_w, 2.0 * resistance * u * u);
            near(
                power.molecular_heat_w,
                100.0 * (24.0 * u * u - 4.0 * u + 2.0),
            );
            near(power.supplied_wall_power_w, 200.0 * (1.0 - u));
            near(power.heat_w, power.supplied_wall_power_w);
            near(thermal.external_heat_j, 0.01 * power.heat_w);
            near(next.energy_j() - old.energy_j(), thermal.external_heat_j);
            near(power.residual_w, 0.0);
            for cell in 0..4 {
                let p = domain.coords(cell);
                let sign = if p[a] == p[b] {
                    -orientation
                } else {
                    orientation
                };
                let matrix = epoch.darcy_resistance_pa_s_m2.as_ref().unwrap()[cell];
                near(matrix[a][a], 150.0);
                near(matrix[b][b], 150.0);
                near(matrix[c][c], 50.0);
                near(matrix[a][b], sign * 50.0);
                near(matrix[b][a], sign * 50.0);
                near(matrix[a][c], 0.0);
                near(matrix[b][c], 0.0);
                near(
                    power.cell_darcy_heat_w.as_ref().unwrap()[cell],
                    0.5 * resistance * u * u,
                );
                let sample = &epoch.porous_samples[cell].material;
                assert_eq!(sample.material().card_identity(), skeleton.content_hash());
                assert_eq!(
                    sample.permeability().material_axes_in_world(),
                    &axes(plane, sign)
                );
                near(
                    next.cells()[&CellId(cell as u64)].specific_h_j_kg,
                    0.01 * power.cell_heat_w[cell],
                );
            }
            speeds.push(u);
        }
        assert!(
            speeds[0] > speeds[1],
            "the same material's orientation changes flow"
        );
    }
}

#[test]
fn rotated_thermal_feedback_preserves_the_native_transient_work_account() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, true);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let skeleton = skeleton([0.5, 1.0, 2.0], true, 810.0);
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(layers(&skeleton, &domain, [0, 1, 2], 1.0))
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let dt = 0.5;
    let a0 = amplitude(0.0, dt, 100.0, 350.0, 100.0);
    let heat0 = dt * 3100.0 * a0 * a0;
    let t1 = 300.0 + heat0 / 400.0;
    let (warm, first, _) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    near(warm.energy_j() - old.energy_j(), heat0);
    near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a0);
    let mu1 = 100.0 * (1.0 - 0.09 * (t1 - 300.0));
    let factor = 1.0 + 0.1 * (t1 - 300.0);
    let resistance = 3.5 * mu1 / factor;
    let a1 = amplitude(a0, dt, mu1, resistance, 100.0);
    let (next, second, thermal) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let step = second.unsteady_mechanical_step().unwrap();
    let audit = second.unsteady_heating_audit().unwrap();
    near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
    near(step.power.darcy_heat_w, 2.0 * resistance * a1 * a1);
    near(step.power.molecular_heat_w, 24.0 * mu1 * a1 * a1);
    near(audit.generated_heat_j, dt * step.power.heat_w);
    near(audit.initial_kinetic_energy_j, 2.0 * a0 * a0);
    near(audit.final_kinetic_energy_j, 2.0 * a1 * a1);
    near(audit.backward_euler_loss_j, 2.0 * (a1 - a0).powi(2));
    near(
        audit.supplied_body_work_j,
        audit.generated_heat_j + audit.final_kinetic_energy_j - audit.initial_kinetic_energy_j
            + audit.backward_euler_loss_j,
    );
    near(thermal.external_heat_j, audit.generated_heat_j);
    near(next.energy_j() - warm.energy_j(), audit.generated_heat_j);
    for cell in 0..4 {
        near(warm.temperature_k(CellId(cell)).unwrap(), t1);
        let epoch = second.applied().unwrap();
        let sample = &epoch.porous_samples[cell as usize].material;
        near(sample.material().query().axes()["T"], t1);
        near(sample.permeability().principal_m2()[0], 0.5 * factor);
        near(
            epoch.samples[cell as usize]
                .transport
                .dynamic_viscosity_pa_s(),
            mu1,
        );
        let matrix = epoch.darcy_resistance_pa_s_m2.as_ref().unwrap()[cell as usize];
        near(matrix[0][0], 1.5 * mu1 / factor);
        near(matrix[0][1].abs(), 0.5 * mu1 / factor);
        let endpoint = &second.endpoint_porous_samples()[cell as usize];
        near(
            endpoint.material.material().query().axes()["T"],
            next.temperature_k(endpoint.cell).unwrap(),
        );
        let initial = &observer.endpoint_porous_samples()[cell as usize].material;
        for (a, b) in initial
            .material()
            .properties()
            .iter()
            .zip(sample.material().properties())
        {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
    }
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
}

#[test]
fn rotated_source_exit_and_cancel_preserve_motion_materials_and_heat() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let skeleton = skeleton([0.5, 1.0, 2.0], false, 300.01);
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(layers(&skeleton, &domain, [0, 1, 2], 1.0))
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let sources = observer.endpoint_porous_samples().to_vec();
    let flow = observer.unsteady_flow().unwrap().clone();
    assert!(matches!(
        old.advance_observed(&observer, 0.5, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Material { .. })
    ));
    assert_eq!(observer.unsteady_flow(), Some(&flow));
    assert_eq!(observer.endpoint_porous_samples(), sources);
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
    near(old.energy_j(), 0.0);
    let (retry, accepted, _) = old
        .advance_observed(&observer, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    assert!(retry.energy_j() > 0.0);
    assert_eq!(accepted.unsteady_heating_audit().unwrap().substeps, 1);
    gate.request();
    assert!(
        old.advance_observed(&observer, 0.001, thermal_step(), || Ok(()))
            .is_err()
    );
    assert_eq!(observer.unsteady_flow(), Some(&flow));
    assert_eq!(observer.endpoint_porous_samples(), sources);
    near(old.energy_j(), 0.0);
}

#[test]
fn cross_axis_drag_at_a_pressure_opening_refuses_before_publication() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    config.faces[1] = FvBoundary::Outlet;
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let skeleton = skeleton([0.5, 1.0, 2.0], false, 810.0);
    let old = initial(&domain, 1.0, &caloric);
    let original = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap();
    // No heating attachment: this refusal belongs to the tensor/opening
    // momentum reduction, not the independent closed-box heating admission.
    let rotated = original
        .clone()
        .with_porous_materials(layers(&skeleton, &domain, [0, 1, 2], 1.0))
        .unwrap();
    let sources = rotated.endpoint_porous_samples().to_vec();
    assert!(matches!(
        old.advance_observed(&rotated, 0.001, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Momentum(ChtError::InvalidInput {
            field: "simple.darcy_resistance_pa_s_m2",
            ..
        }))
    ));
    assert!(rotated.applied().is_none());
    assert_eq!(rotated.endpoint_porous_samples(), sources);
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);

    let diagonal = original
        .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
        .unwrap();
    let (next, accepted, report) = old
        .advance_observed(&diagonal, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), old.cells());
    near(report.external_heat_j, 0.0);
    near(accepted.applied().unwrap().flow.field.flux_x(2, 0, 0), 0.0);
    assert!(accepted.applied().unwrap().viscous_power.is_none());
}
