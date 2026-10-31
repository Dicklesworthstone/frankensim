//! G1/G3/G4: calibrated nonlinear porous drag, heat and source feedback.
//! Synthetic coefficients test the discrete physical model, not a real filter.
use super::*;

#[path = "forchheimer/rotated.rs"]
mod rotated;

const BETA_DIMS: Dims = Dims([-1, 0, 0, 0, 0, 0]);
const BETA_NAMES: [&str; 3] = [
    "forchheimer_beta_1",
    "forchheimer_beta_2",
    "forchheimer_beta_3",
];

fn nonlinear_skeleton(beta: [f64; 3], changing: bool, upper_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (names, dims, values, slope) in [
        (K_NAMES, INTRINSIC_PERMEABILITY_DIMS, [0.5, 1.0, 2.0], 0.1),
        (BETA_NAMES, BETA_DIMS, beta, 0.2),
    ] {
        for (name, value) in names.into_iter().zip(values) {
            let mut claim = scalar(name, dims, value);
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", absolute(), 300.0, upper_k)
                .with_quantity("pressure", pressure_quantity(), 100000.0, 200000.0);
            if changing {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, value),
                        (upper_k, value * (1.0 + slope * (upper_k - 300.0))),
                    ],
                    dims,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(claim).unwrap();
        }
    }
    card(claims)
}

fn nonlinear_law(
    card: &MaterialCard,
    axes: [[f64; 3]; 3],
    speed_limit: f64,
) -> PermeabilityTemperatureLaw<'_> {
    let keys = K_NAMES.map(|name| PropertyKey::new(name, INTRINSIC_PERMEABILITY_DIMS));
    let beta_keys = BETA_NAMES.map(|name| PropertyKey::new(name, BETA_DIMS));
    let requirements = keys
        .iter()
        .chain(&beta_keys)
        .map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        })
        .collect::<Vec<_>>();
    let point = QueryPoint::new()
        .with_quantity("T", absolute(), 300.0)
        .unwrap()
        .with_quantity("pressure", pressure_quantity(), 100000.0)
        .unwrap();
    let resolved = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PermeabilityTemperatureLaw::try_new(card, resolved, keys, axes, "T")
        .unwrap()
        .with_forchheimer(beta_keys, speed_limit)
        .unwrap()
}

// Independent one-degree-of-freedom reduction of the four free MAC faces:
// L*a + Q*|a|a = F. Stable signed root; no solver coefficients are inspected.
fn nonlinear_amplitude(linear: f64, quadratic: f64, force: f64) -> f64 {
    2.0 * force / (linear + (linear * linear + 4.0 * quadratic * force.abs()).sqrt())
}

#[test]
fn material_beta_and_signed_axis_permutation_drive_motion_and_cubic_heat() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut speeds = Vec::new();
    for scale in [1.0, 2.0] {
        let beta = [4000.0, 1000.0, 2000.0].map(|b| scale * b);
        let skeleton = nonlinear_skeleton(beta, false, 810.0);
        for (axes, grid_k, grid_beta) in [
            (GRID_ALIGNED_AXES, [0.5, 1.0, 2.0], beta),
            (
                [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
                [2.0, 1.0, 0.5],
                [beta[2], beta[1], beta[0]],
            ),
        ] {
            for wall in [-1.0, 1.0] {
                let mut config = driven_box();
                config.faces[3] = FvBoundary::Wall {
                    velocity: [wall, 0.0, 0.0],
                };
                let observer = input(&ids, &liquid, &caloric)
                    .feedback(&old, &domain, &config, (), &gate, 10)
                    .unwrap()
                    .with_porous_materials(vec![Some(nonlinear_law(&skeleton, axes, 1.0)); 4])
                    .unwrap()
                    .with_viscous_heating(controls())
                    .unwrap();
                let (next, accepted, thermal) = old
                    .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
                    .unwrap();
                let epoch = accepted.applied().unwrap();
                let power = epoch.viscous_power.as_ref().unwrap();
                let linear = 100.0 / grid_k[0] + 100.0 / grid_k[1];
                // rho=1; sourced beta has NO implicit factor of one half.
                let quadratic = grid_beta[0] + grid_beta[1];
                let a = nonlinear_amplitude(1200.0 + linear, quadratic, 100.0 * wall);
                near(epoch.flow.field.flux_x(1, 1, 0), a);
                near(epoch.flow.field.flux_x(1, 0, 0), -a);
                near(epoch.flow.field.flux_y(0, 1, 0), a);
                near(epoch.flow.field.flux_y(1, 1, 0), -a);
                near(power.darcy_heat_w, 2.0 * linear * a * a);
                near(power.forchheimer_heat_w, 2.0 * quadratic * a.abs().powi(3));
                near(
                    power.molecular_heat_w,
                    100.0 * (24.0 * a * a - 4.0 * wall * a + 2.0),
                );
                near(power.supplied_wall_power_w, 200.0 * wall * (wall - a));
                near(
                    power.heat_w,
                    power.molecular_heat_w + power.darcy_heat_w + power.forchheimer_heat_w,
                );
                near(power.heat_w, power.supplied_wall_power_w);
                near(thermal.external_heat_j, 0.01 * power.heat_w);
                near(next.energy_j() - old.energy_j(), thermal.external_heat_j);
                near(
                    accepted.viscous_heating_audit().unwrap().forchheimer_heat_j,
                    0.01 * power.forchheimer_heat_w,
                );
                near(power.residual_w, 0.0);
                for cell in 0..4 {
                    near(
                        power.cell_forchheimer_heat_w.as_ref().unwrap()[cell],
                        0.5 * quadratic * a.abs().powi(3),
                    );
                    let resistance = epoch.forchheimer_resistance.as_ref().unwrap()[cell]
                        .as_ref()
                        .unwrap();
                    assert_eq!(resistance.principal_kg_m4(), &beta);
                    assert_eq!(resistance.material_axes_in_world(), &axes);
                    assert_eq!(
                        epoch.porous_samples[cell].material.forchheimer_beta_per_m(),
                        Some(&beta)
                    );
                    assert_eq!(
                        epoch.porous_samples[cell]
                            .material
                            .material()
                            .card_identity(),
                        skeleton.content_hash()
                    );
                }
                speeds.push(a.abs());
            }
        }
    }
    // Only original beta data changed. Every orientation slows consistently.
    for i in 0..4 {
        assert!(speeds[i + 4] < speeds[i]);
    }
}

#[test]
fn nonlinear_porous_heat_updates_original_sources_and_closes_transient_energy() {
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
        .with_porous_materials(vec![
            Some(nonlinear_law(&skeleton, GRID_ALIGNED_AXES, 1.0));
            4
        ])
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let dt = 0.5;
    let a0 = nonlinear_amplitude(4.0 / dt + 3000.0, 10000.0, 400.0);
    let heat0 = dt * (3000.0 * a0 * a0 + 10000.0 * a0.powi(3));
    let t1 = 300.0 + heat0 / 400.0;
    let (warm, first, _) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a0);
    near(warm.energy_j() - old.energy_j(), heat0);
    let mu1 = 100.0 * (1.0 - 0.09 * (t1 - 300.0));
    let k_factor = 1.0 + 0.1 * (t1 - 300.0);
    let beta_factor = 1.0 + 0.2 * (t1 - 300.0);
    let linear = 3.0 * mu1 / k_factor;
    let quadratic = 5000.0 * beta_factor;
    let a1 = nonlinear_amplitude(
        4.0 / dt + 24.0 * mu1 + 2.0 * linear,
        2.0 * quadratic,
        4.0 * a0 / dt + 400.0,
    );
    let (next, second, thermal) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let step = second.unsteady_mechanical_step().unwrap();
    let audit = second.unsteady_heating_audit().unwrap();
    near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
    near(step.power.darcy_heat_w, 2.0 * linear * a1 * a1);
    near(step.power.forchheimer_heat_w, 2.0 * quadratic * a1.powi(3));
    near(step.power.molecular_heat_w, 24.0 * mu1 * a1 * a1);
    near(audit.generated_heat_j, dt * step.power.heat_w);
    near(audit.forchheimer_heat_j, dt * step.power.forchheimer_heat_w);
    near(
        audit.generated_heat_j,
        audit.molecular_heat_j + audit.darcy_heat_j + audit.forchheimer_heat_j,
    );
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
        let source = &second.applied().unwrap().porous_samples[cell as usize].material;
        near(source.material().query().axes()["T"], t1);
        near(source.permeability().principal_m2()[0], 0.5 * k_factor);
        near(
            source.forchheimer_beta_per_m().unwrap()[0],
            beta[0] * beta_factor,
        );
        let endpoint = &second.endpoint_porous_samples()[cell as usize].material;
        near(
            endpoint.material().query().axes()["T"],
            next.temperature_k(CellId(cell)).unwrap(),
        );
        assert_eq!(endpoint.material().card_identity(), skeleton.content_hash());
        near(endpoint.material().query().axes()["pressure"], 100000.0);
    }
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
}

#[test]
fn zero_beta_recovers_rotated_darcy_exactly_and_positive_rotated_beta_is_admitted() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let old = initial(&domain, 1.0, &caloric);
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let axes = [[h, h, 0.0], [-h, h, 0.0], [0.0, 0.0, 1.0]];
    let skeleton = nonlinear_skeleton([0.0; 3], false, 810.0);
    let original = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap();
    let darcy = original
        .clone()
        .with_porous_materials(vec![Some(permeability(&skeleton, axes)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let inertial = original
        .clone()
        .with_porous_materials(vec![Some(nonlinear_law(&skeleton, axes, 1.0)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let (linear_next, linear, _) = old
        .advance_observed(&darcy, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    let (zero_next, zero, _) = old
        .advance_observed(&inertial, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(linear_next.cells(), zero_next.cells());
    assert_eq!(
        linear.applied().unwrap().flow.field,
        zero.applied().unwrap().flow.field
    );
    assert_eq!(
        zero.applied()
            .unwrap()
            .viscous_power
            .as_ref()
            .unwrap()
            .forchheimer_heat_w,
        0.0
    );
    assert!(zero.applied().unwrap().forchheimer_resistance.is_none());
    let positive = nonlinear_skeleton([1.0; 3], false, 810.0);
    let rotated = original
        .with_porous_materials(vec![Some(nonlinear_law(&positive, axes, 1.0)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let (_, accepted, _) = old
        .advance_observed(&rotated, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    assert!(
        accepted
            .applied()
            .unwrap()
            .viscous_power
            .as_ref()
            .unwrap()
            .forchheimer_heat_w
            > 0.0
    );
}

#[test]
fn calibration_speed_and_source_domains_fail_without_advancing_any_owner() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let old = initial(&domain, 1.0, &caloric);
    for (upper_k, speed_limit) in [(810.0, 0.01), (300.01, 1.0)] {
        let skeleton = nonlinear_skeleton([4000.0, 1000.0, 2000.0], false, upper_k);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_porous_materials(vec![
                Some(nonlinear_law(
                    &skeleton,
                    GRID_ALIGNED_AXES,
                    speed_limit
                ));
                4
            ])
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
        let mut calls = 0usize;
        let (retry, accepted, _) = old
            .advance_observed(&observer, 0.00001, thermal_step(), || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        assert!(retry.energy_j() > 0.0);
        assert_eq!(accepted.unsteady_heating_audit().unwrap().substeps, 1);
        let mut replay_calls = 0usize;
        assert!(
            old.advance_observed(&observer, 0.00001, thermal_step(), || {
                replay_calls += 1;
                if replay_calls == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
        );
        assert_eq!(observer.unsteady_flow(), Some(&flow));
        assert_eq!(observer.endpoint_porous_samples(), sources);
        near(old.energy_j(), 0.0);
    }
}
