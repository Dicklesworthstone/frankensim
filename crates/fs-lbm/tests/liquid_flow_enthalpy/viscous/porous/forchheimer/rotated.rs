//! Independent eight-tuple work reduction for rotated, calibrated porous inertia.
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
    speed_limit: f64,
) -> Vec<Option<PermeabilityTemperatureLaw<'a>>> {
    (0..domain.cell_count())
        .map(|cell| {
            let p = domain.coords(cell);
            let sign = if p[plane[0]] == p[plane[1]] {
                -orientation
            } else {
                orientation
            };
            Some(nonlinear_law(source, axes(plane, sign), speed_limit))
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

fn cubic_work_coefficient(beta: [f64; 3], orientation: f64) -> f64 {
    // Each cell has one moving and one fixed face per planar component.
    // In its material axes the four distinct pairs project to 0, |u|/sqrt(2)
    // or sqrt(2)*|u|. Summing their cubes, with the 1/4 pair weights and
    // four cells, gives this coefficient of total dissipated |u|^3.
    let active = if orientation > 0.0 { beta[0] } else { beta[1] };
    (beta[0] + beta[1] + 4.0 * active) * std::f64::consts::FRAC_1_SQRT_2
}

#[test]
fn rotated_material_data_controls_signed_flow_and_cubic_heat_in_all_planes() {
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    for plane in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
        let [a, b, c] = plane;
        let mut n = [2; 3];
        n[c] = 1;
        let domain = VoxelDomain::new(n[0], n[1], n[2], 1.0).unwrap();
        let old = initial(&domain, 1.0, &caloric);
        let mut speeds = Vec::new();
        for scale in [1.0, 2.0] {
            let beta = [4000.0, 1000.0, 2000.0].map(|v| scale * v);
            let skeleton = nonlinear_skeleton(beta, false, 810.0);
            for orientation in [-1.0, 1.0] {
                for wall in [-1.0, 1.0] {
                    let mut config = driven_box();
                    config.faces = [FvBoundary::wall(); 6];
                    config.faces[2 * c] = FvBoundary::Symmetry;
                    config.faces[2 * c + 1] = FvBoundary::Symmetry;
                    let mut velocity = [0.0; 3];
                    velocity[a] = wall;
                    config.faces[2 * b + 1] = FvBoundary::Wall { velocity };
                    let observer = input(&ids, &liquid, &caloric)
                        .feedback(&old, &domain, &config, (), &gate, 10)
                        .unwrap()
                        .with_porous_materials(layers(&skeleton, &domain, plane, orientation, 1.0))
                        .unwrap()
                        .with_viscous_heating(controls())
                        .unwrap();
                    let (next, accepted, thermal) = old
                        .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
                        .unwrap();
                    let epoch = accepted.applied().unwrap();
                    let power = epoch.viscous_power.as_ref().unwrap();
                    let resistance = 300.0 + orientation * 50.0;
                    let cubic = cubic_work_coefficient(beta, orientation);
                    let u = nonlinear_amplitude(1200.0 + resistance, cubic / 2.0, 100.0 * wall);
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
                    near(power.forchheimer_heat_w, cubic * u.abs().powi(3));
                    near(
                        power.molecular_heat_w,
                        100.0 * (24.0 * u * u - 4.0 * wall * u + 2.0 * wall * wall),
                    );
                    near(power.supplied_wall_power_w, 200.0 * wall * (wall - u));
                    near(power.heat_w, power.supplied_wall_power_w);
                    near(
                        power.heat_w,
                        power.molecular_heat_w + power.darcy_heat_w + power.forchheimer_heat_w,
                    );
                    near(power.residual_w, 0.0);
                    near(thermal.external_heat_j, 0.01 * power.heat_w);
                    near(next.energy_j() - old.energy_j(), thermal.external_heat_j);
                    for cell in 0..4 {
                        near(
                            power.cell_forchheimer_heat_w.as_ref().unwrap()[cell],
                            cubic * u.abs().powi(3) / 4.0,
                        );
                        near(
                            power.cell_darcy_heat_w.as_ref().unwrap()[cell],
                            resistance * u * u / 2.0,
                        );
                        let p = domain.coords(cell);
                        let sign = if p[a] == p[b] {
                            -orientation
                        } else {
                            orientation
                        };
                        let native = epoch.forchheimer_resistance.as_ref().unwrap()[cell]
                            .as_ref()
                            .unwrap();
                        assert_eq!(native.principal_kg_m4(), &beta);
                        assert_eq!(native.material_axes_in_world(), &axes(plane, sign));
                        let sample = &epoch.porous_samples[cell].material;
                        assert_eq!(sample.material().card_identity(), skeleton.content_hash());
                        assert_eq!(sample.forchheimer_beta_per_m(), Some(&beta));
                        near(
                            next.cells()[&CellId(cell as u64)].specific_h_j_kg,
                            0.01 * power.cell_heat_w[cell],
                        );
                    }
                    speeds.push(u.abs());
                }
            }
        }
        for i in 0..4 {
            assert!(speeds[i + 4] < speeds[i]);
        }
        assert!(
            speeds[0] > speeds[2],
            "material orientation changes nonlinear resistance"
        );
    }
}

#[test]
fn rotating_nonlinear_sources_evolve_on_the_same_transient_energy_timeline() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, true);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let beta = [4000.0, 1000.0, 2000.0];
    let skeleton = nonlinear_skeleton(beta, true, 810.0);
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(layers(&skeleton, &domain, [0, 1, 2], 1.0, 1.0))
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let dt = 0.5;
    let cubic0 = cubic_work_coefficient(beta, 1.0);
    let u0 = nonlinear_amplitude(4.0 / dt + 3100.0, cubic0, 400.0);
    let heat0 = dt * (3100.0 * u0 * u0 + cubic0 * u0.powi(3));
    let t1 = 300.0 + heat0 / 400.0;
    let (warm, first, _) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), u0);
    near(warm.energy_j() - old.energy_j(), heat0);
    let mu1 = 100.0 * (1.0 - 0.09 * (t1 - 300.0));
    let k_factor = 1.0 + 0.1 * (t1 - 300.0);
    let beta_factor = 1.0 + 0.2 * (t1 - 300.0);
    let resistance = 3.5 * mu1 / k_factor;
    let cubic1 = cubic0 * beta_factor;
    let u1 = nonlinear_amplitude(
        4.0 / dt + 24.0 * mu1 + 2.0 * resistance,
        cubic1,
        4.0 * u0 / dt + 400.0,
    );
    let (next, second, thermal) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let step = second.unsteady_mechanical_step().unwrap();
    let audit = second.unsteady_heating_audit().unwrap();
    near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), u1);
    near(step.power.darcy_heat_w, 2.0 * resistance * u1 * u1);
    near(step.power.forchheimer_heat_w, cubic1 * u1.powi(3));
    near(step.power.molecular_heat_w, 24.0 * mu1 * u1 * u1);
    near(audit.generated_heat_j, dt * step.power.heat_w);
    near(audit.forchheimer_heat_j, dt * cubic1 * u1.powi(3));
    near(audit.initial_kinetic_energy_j, 2.0 * u0 * u0);
    near(audit.final_kinetic_energy_j, 2.0 * u1 * u1);
    near(audit.backward_euler_loss_j, 2.0 * (u1 - u0).powi(2));
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
        near(sample.permeability().principal_m2()[0], 0.5 * k_factor);
        near(
            sample.forchheimer_beta_per_m().unwrap()[0],
            beta[0] * beta_factor,
        );
        near(
            epoch.samples[cell as usize]
                .transport
                .dynamic_viscosity_pa_s(),
            mu1,
        );
        let native = epoch.forchheimer_resistance.as_ref().unwrap()[cell as usize]
            .as_ref()
            .unwrap();
        near(native.principal_kg_m4()[0], beta[0] * beta_factor);
        assert_eq!(
            native.material_axes_in_world(),
            sample.permeability().material_axes_in_world()
        );
        let endpoint = &second.endpoint_porous_samples()[cell as usize].material;
        near(
            endpoint.material().query().axes()["T"],
            next.temperature_k(CellId(cell)).unwrap(),
        );
        let original = &observer.endpoint_porous_samples()[cell as usize].material;
        for (new, old) in endpoint
            .material()
            .properties()
            .iter()
            .zip(original.material().properties())
        {
            assert_eq!(new.answer().receipt.selected, old.answer().receipt.selected);
        }
    }
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
}

#[test]
fn principal_tuple_speed_bound_refuses_even_when_every_grid_face_is_below_it() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let beta = [4000.0, 1000.0, 2000.0];
    let skeleton = nonlinear_skeleton(beta, false, 810.0);
    let old = initial(&domain, 1.0, &caloric);
    let u = nonlinear_amplitude(3108.0, cubic_work_coefficient(beta, 1.0), 400.0);
    let maximum = 1.2 * u;
    assert!(u < maximum && std::f64::consts::SQRT_2 * u > maximum);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(layers(&skeleton, &domain, [0, 1, 2], 1.0, maximum))
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let flow = observer.unsteady_flow().unwrap().clone();
    let sources = observer.endpoint_porous_samples().to_vec();
    assert!(
        old.advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .is_err()
    );
    assert_eq!(observer.unsteady_flow(), Some(&flow));
    assert_eq!(observer.endpoint_porous_samples(), sources);
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
    assert!(observer.applied().is_none());
    near(old.energy_j(), 0.0);
    let (retry, accepted, _) = old
        .advance_observed(&observer, 0.00001, thermal_step(), || Ok(()))
        .unwrap();
    assert!(retry.energy_j() > 0.0);
    assert_eq!(accepted.unsteady_heating_audit().unwrap().substeps, 1);
}

#[test]
fn nonlinear_orientation_alone_requires_a_clear_fluid_layer_at_openings() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let mut config = throughflow();
    config.faces[0] = FvBoundary::Inlet {
        velocity: [0.01, 0.0, 0.0],
    };
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let original = nonlinear_skeleton([4000.0, 1000.0, 2000.0], false, 810.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        if K_NAMES.contains(&claim.key.name()) {
            claim.value = PropertyValue::Scalar {
                value: 1.0,
                dims: INTRINSIC_PERMEABILITY_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let skeleton = card(claims);
    let old = initial(&domain, 1.0, &caloric);
    let base = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap();
    let frame = axes([0, 1, 2], 1.0);
    let direct = base
        .clone()
        .with_porous_materials(vec![Some(nonlinear_law(&skeleton, frame, 1.0)); 6])
        .unwrap();
    let sources = direct.endpoint_porous_samples().to_vec();
    // K is isotropic, so this cannot be the earlier Darcy cross-axis refusal.
    assert!(matches!(
        old.advance_observed(&direct, 0.001, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Momentum(ChtError::InvalidInput {
            field: "simple.forchheimer_resistance",
            ..
        }))
    ));
    assert!(direct.applied().is_none());
    assert_eq!(direct.endpoint_porous_samples(), sources);
    near(old.energy_j(), 0.0);
    let insert = (0..6)
        .map(|cell| (domain.coords(cell)[0] == 1).then(|| nonlinear_law(&skeleton, frame, 1.0)))
        .collect();
    let separated = base.with_porous_materials(insert).unwrap();
    let (next, accepted, thermal) = old
        .advance_observed(&separated, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    let epoch = accepted.applied().unwrap();
    near(
        epoch.flow.field.flux_x(3, 0, 0) + epoch.flow.field.flux_x(3, 1, 0),
        0.02,
    );
    assert!(epoch.viscous_power.is_none());
    assert_eq!(next.cells(), old.cells());
    near(thermal.external_heat_j, 0.0);
    assert_eq!(
        epoch
            .forchheimer_resistance
            .as_ref()
            .unwrap()
            .iter()
            .filter(|c| c.is_some())
            .count(),
        2
    );
}
