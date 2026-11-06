//! G1/G3/G4: source-bound vector-speed inertia on the original MAC heat owner.
//! Synthetic coefficients; these tests are discrete oracles, not filter measurements.
use super::*;

const VECTOR_BETA: &str = "isotropic_forchheimer_beta";

fn vector_skeleton(beta: f64, changing: bool, upper_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (K_NAMES[0], INTRINSIC_PERMEABILITY_DIMS, 1.0, 0.1),
        (K_NAMES[1], INTRINSIC_PERMEABILITY_DIMS, 1.0, 0.1),
        (K_NAMES[2], INTRINSIC_PERMEABILITY_DIMS, 1.0, 0.1),
        (VECTOR_BETA, BETA_DIMS, beta, 0.2),
    ] {
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
    card(claims)
}

fn vector_law(
    card: &MaterialCard,
    axes: [[f64; 3]; 3],
    speed_limit: f64,
) -> PermeabilityTemperatureLaw<'_> {
    let keys = K_NAMES.map(|name| PropertyKey::new(name, INTRINSIC_PERMEABILITY_DIMS));
    let beta = PropertyKey::new(VECTOR_BETA, BETA_DIMS);
    let requirements = keys
        .iter()
        .chain([&beta])
        .map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        })
        .collect::<Vec<_>>();
    let point = QueryPoint::new()
        .with_quantity("T", absolute(), 300.0)
        .unwrap()
        .with_quantity("pressure", pressure_quantity(), 100000.0)
        .unwrap();
    let material = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PermeabilityTemperatureLaw::try_new(card, material, keys, axes, "T")
        .unwrap()
        .with_isotropic_forchheimer(beta, speed_limit)
        .unwrap()
}

fn cubic_coefficient(beta: f64) -> f64 {
    // Each planar cell has tuple speeds 0, |a|, |a|, sqrt(2)*|a|, with
    // probability 1/4 each. Four cells give this total coefficient of |a|^3.
    2.0 * (1.0 + 2.0_f64.sqrt()) * beta
}

#[test]
fn vector_speed_circulation_rotates_and_source_beta_changes_the_actual_flow() {
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
        let mut previous_speed = f64::INFINITY;
        for beta in [1000.0, 4000.0] {
            let skeleton = vector_skeleton(beta, false, 810.0);
            let cubic = cubic_coefficient(beta);
            let mut reference = None;
            for axes in [
                GRID_ALIGNED_AXES,
                [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]],
            ] {
                for sign in [-1.0, 1.0] {
                    let mut config = driven_box();
                    config.faces = [FvBoundary::wall(); 6];
                    config.faces[2 * c] = FvBoundary::Symmetry;
                    config.faces[2 * c + 1] = FvBoundary::Symmetry;
                    let mut velocity = [0.0; 3];
                    velocity[a] = sign;
                    config.faces[2 * b + 1] = FvBoundary::Wall { velocity };
                    let observer = input(&ids, &liquid, &caloric)
                        .feedback(&old, &domain, &config, (), &gate, 10)
                        .unwrap()
                        .with_porous_materials(vec![Some(vector_law(&skeleton, axes, 1.0)); 4])
                        .unwrap()
                        .with_viscous_heating(controls())
                        .unwrap();
                    let (next, accepted, thermal) = old
                        .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
                        .unwrap();
                    let epoch = accepted.applied().unwrap();
                    let power = epoch.viscous_power.as_ref().unwrap();
                    let u = nonlinear_amplitude(1400.0, cubic / 2.0, 100.0 * sign);
                    let mut at = [0; 3];
                    at[a] = 1;
                    at[b] = 1;
                    let actual = match a {
                        0 => epoch.flow.field.flux_x(at[0], at[1], at[2]),
                        1 => epoch.flow.field.flux_y(at[0], at[1], at[2]),
                        _ => epoch.flow.field.flux_z(at[0], at[1], at[2]),
                    };
                    near(actual, u);
                    near(power.forchheimer_heat_w, cubic * u.abs().powi(3));
                    near(power.darcy_heat_w, 400.0 * u * u);
                    near(
                        power.molecular_heat_w,
                        100.0 * (24.0 * u * u - 4.0 * sign * u + 2.0),
                    );
                    near(power.supplied_wall_power_w, 200.0 * sign * (sign - u));
                    near(power.heat_w, power.supplied_wall_power_w);
                    near(
                        power.heat_w,
                        power.molecular_heat_w + power.darcy_heat_w + power.forchheimer_heat_w,
                    );
                    near(power.residual_w, 0.0);
                    near(next.energy_j() - old.energy_j(), thermal.external_heat_j);
                    near(thermal.external_heat_j, 0.01 * power.heat_w);
                    for cell in 0..4 {
                        near(
                            power.cell_forchheimer_heat_w.as_ref().unwrap()[cell],
                            cubic * u.abs().powi(3) / 4.0,
                        );
                        let native = epoch.forchheimer_resistance.as_ref().unwrap()[cell]
                            .as_ref()
                            .unwrap();
                        assert_eq!(native.isotropic_kg_m4(), Some(beta));
                        assert_eq!(native.principal_kg_m4(), None);
                        assert_eq!(native.material_axes_in_world(), None);
                        let sample = &epoch.porous_samples[cell].material;
                        assert_eq!(sample.isotropic_forchheimer_beta_per_m(), Some(beta));
                        assert_eq!(sample.material().card_identity(), skeleton.content_hash());
                    }
                    if sign > 0.0 {
                        if let Some(previous) = &reference {
                            assert_eq!(previous, &epoch.flow.field);
                        } else {
                            reference = Some(epoch.flow.field.clone());
                        }
                    }
                }
            }
            let speed = nonlinear_amplitude(1400.0, cubic / 2.0, 100.0);
            assert!(speed < previous_speed);
            previous_speed = speed;
        }
    }
}

#[test]
fn vector_drag_heat_refreshes_original_sources_and_closes_native_transient_energy() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, true);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let skeleton = vector_skeleton(2000.0, true, 810.0);
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(vec![Some(vector_law(&skeleton, GRID_ALIGNED_AXES, 1.0)); 4])
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let dt = 0.5;
    let cubic0 = cubic_coefficient(2000.0);
    let u0 = nonlinear_amplitude(4.0 / dt + 2800.0, cubic0, 400.0);
    let heat0 = dt * (2800.0 * u0 * u0 + cubic0 * u0.powi(3));
    let t1 = 300.0 + heat0 / 400.0;
    let (warm, first, _) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), u0);
    near(warm.energy_j() - old.energy_j(), heat0);
    let mu1 = 100.0 * (1.0 - 0.09 * (t1 - 300.0));
    let k1 = 1.0 + 0.1 * (t1 - 300.0);
    let beta1 = 2000.0 * (1.0 + 0.2 * (t1 - 300.0));
    let linear1 = 24.0 * mu1 + 4.0 * mu1 / k1;
    let cubic1 = cubic_coefficient(beta1);
    let u1 = nonlinear_amplitude(4.0 / dt + linear1, cubic1, 4.0 * u0 / dt + 400.0);
    let (next, second, thermal) = warm
        .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let step = second.unsteady_mechanical_step().unwrap();
    let audit = second.unsteady_heating_audit().unwrap();
    near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), u1);
    near(step.power.forchheimer_heat_w, cubic1 * u1.powi(3));
    near(step.power.darcy_heat_w, 4.0 * mu1 / k1 * u1 * u1);
    near(step.power.molecular_heat_w, 24.0 * mu1 * u1 * u1);
    near(
        audit.generated_heat_j,
        dt * (linear1 * u1 * u1 + cubic1 * u1.powi(3)),
    );
    near(audit.forchheimer_heat_j, dt * cubic1 * u1.powi(3));
    near(audit.initial_kinetic_energy_j, 2.0 * u0 * u0);
    near(audit.final_kinetic_energy_j, 2.0 * u1 * u1);
    near(audit.backward_euler_loss_j, 2.0 * (u1 - u0).powi(2));
    near(
        audit.supplied_body_work_j,
        audit.generated_heat_j + audit.final_kinetic_energy_j - audit.initial_kinetic_energy_j
            + audit.backward_euler_loss_j,
    );
    near(next.energy_j() - warm.energy_j(), audit.generated_heat_j);
    near(thermal.external_heat_j, audit.generated_heat_j);
    for cell in 0..4 {
        let source = &second.applied().unwrap().porous_samples[cell].material;
        near(source.material().query().axes()["T"], t1);
        near(source.permeability().principal_m2()[0], k1);
        near(source.isotropic_forchheimer_beta_per_m().unwrap(), beta1);
        let endpoint = &second.endpoint_porous_samples()[cell].material;
        near(
            endpoint.material().query().axes()["T"],
            next.temperature_k(CellId(cell as u64)).unwrap(),
        );
        assert_eq!(endpoint.material().card_identity(), skeleton.content_hash());
        near(endpoint.material().query().axes()["pressure"], 100000.0);
    }
}

#[test]
fn vector_zero_is_exact_darcy_and_full_speed_refusal_is_atomic() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let liquid_card = source(100.0, false);
    let liquid = liquid(&liquid_card);
    let caloric = liquid_caloric(&liquid_card);
    let old = initial(&domain, 1.0, &caloric);
    let skeleton = vector_skeleton(0.0, false, 810.0);
    let original = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap();
    let darcy = original
        .clone()
        .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let vector = original
        .clone()
        .with_porous_materials(vec![Some(vector_law(&skeleton, GRID_ALIGNED_AXES, 1.0)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let (linear_next, linear, _) = old
        .advance_observed(&darcy, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    let (zero_next, zero, _) = old
        .advance_observed(&vector, 0.01, thermal_step(), || Ok(()))
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

    let skeleton = vector_skeleton(2000.0, false, 810.0);
    let u = nonlinear_amplitude(1400.0, cubic_coefficient(2000.0) / 2.0, 100.0);
    // Every individual face is within 1.1*u, but the same native tuples
    // have sqrt(2)*u speed. A component test or a cell average would admit it.
    let limited = original
        .with_porous_materials(vec![
            Some(vector_law(&skeleton, GRID_ALIGNED_AXES, 1.1 * u));
            4
        ])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let sources = limited.endpoint_porous_samples().to_vec();
    assert!(
        old.advance_observed(&limited, 0.01, thermal_step(), || Ok(()))
            .is_err()
    );
    assert!(limited.applied().is_none());
    assert_eq!(limited.endpoint_porous_samples(), sources);
    assert_eq!(
        limited.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
    assert_eq!(old.energy_j(), 0.0);
}
