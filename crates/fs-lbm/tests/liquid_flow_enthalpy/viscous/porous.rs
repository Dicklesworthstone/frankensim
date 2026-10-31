//! Original skeleton sources share native momentum, temperature and Darcy heat.
use super::*;
use fs_lbm::conjugate::{
    UnsteadyHeatingAudit, creeping_flow_with_darcy_heat, simple_flow_with_material_resistance,
};
use fs_material::porous::permeability::INTRINSIC_PERMEABILITY_DIMS;
use fs_material::porous::permeability::temperature::PermeabilityTemperatureLaw;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

#[path = "porous/rotated.rs"]
mod rotated;
#[path = "porous/forchheimer.rs"]
mod forchheimer;

const K_NAMES: [&str; 3] = ["intrinsic_k_1", "intrinsic_k_2", "intrinsic_k_3"];
const GRID_ALIGNED_AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn skeleton(principal: [f64; 3], changing: bool, upper_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, k) in K_NAMES.into_iter().zip(principal) {
        let mut claim = scalar(name, INTRINSIC_PERMEABILITY_DIMS, k);
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", absolute(), 300.0, upper_k)
            .with_quantity("pressure", pressure_quantity(), 100000.0, 200000.0);
        if changing {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, k), (upper_k, k * (1.0 + 0.1 * (upper_k - 300.0)))],
                dims: INTRINSIC_PERMEABILITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

fn permeability(card: &MaterialCard, axes: [[f64; 3]; 3]) -> PermeabilityTemperatureLaw<'_> {
    let keys = K_NAMES.map(|name| PropertyKey::new(name, INTRINSIC_PERMEABILITY_DIMS));
    let requirements = keys
        .iter()
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
    PermeabilityTemperatureLaw::try_new(card, material, keys, axes, "T").unwrap()
}

fn forces(f: f64) -> Vec<[f64; 3]> {
    vec![[-f, f, 0.0], [-f, -f, 0.0], [f, f, 0.0], [f, -f, 0.0]]
}

fn amplitude(old: f64, dt: f64, mu: f64, beta_xy: f64, force: f64) -> f64 {
    // Independent four-face constrained MAC reduction. Only x/y faces move;
    // each directional coefficient occurs twice in the mechanical quadratic.
    (4.0 * old / dt + 4.0 * force) / (4.0 / dt + 24.0 * mu + 2.0 * beta_xy)
}

#[test]
fn original_axis_permutation_changes_native_drag_and_darcy_heat_once() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let fluid_card = source(100.0, false);
    let liquid = liquid(&fluid_card);
    let caloric = liquid_caloric(&fluid_card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let principal = [0.5, 1.0, 2.0];
    let skeleton = skeleton(principal, false, 810.0);
    let mut speeds = Vec::new();
    for (axes, grid_k) in [
        (GRID_ALIGNED_AXES, principal),
        (
            [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
            [2.0, 1.0, 0.5],
        ),
    ] {
        let law = permeability(&skeleton, axes);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_porous_materials(vec![Some(law); 4])
            .unwrap()
            .with_viscous_heating(controls())
            .unwrap();
        let dt = 0.01;
        let (next, accepted, thermal) = old
            .advance_observed(&observer, dt, thermal_step(), || Ok(()))
            .unwrap();
        let epoch = accepted.applied().unwrap();
        let beta = grid_k.map(|k| 100.0 / k);
        let beta_xy = beta[0] + beta[1];
        let a = 100.0 / (1200.0 + beta_xy);
        let power = epoch.viscous_power.as_ref().unwrap();
        near(epoch.flow.field.flux_x(1, 1, 0), a);
        near(epoch.flow.field.flux_y(0, 1, 0), a);
        near(power.darcy_heat_w, 2.0 * beta_xy * a * a);
        near(
            power.molecular_heat_w,
            100.0 * (24.0 * a * a - 4.0 * a + 2.0),
        );
        near(power.heat_w, power.molecular_heat_w + power.darcy_heat_w);
        near(power.supplied_wall_power_w, power.heat_w);
        near(thermal.external_heat_j, dt * power.heat_w);
        near(next.energy_j() - old.energy_j(), thermal.external_heat_j);
        near(
            accepted.viscous_heating_audit().unwrap().darcy_heat_j,
            dt * power.darcy_heat_w,
        );
        assert!(epoch.mushy_drag_pa_s_m2.is_none());
        assert_eq!(
            epoch.darcy_resistance_pa_s_m2.as_ref().unwrap(),
            &vec![[[beta[0], 0.0, 0.0], [0.0, beta[1], 0.0], [0.0, 0.0, beta[2]]]; 4]
        );
        for c in 0..4 {
            near(
                power.cell_darcy_heat_w.as_ref().unwrap()[c],
                0.5 * beta_xy * a * a,
            );
            near(
                next.cells()[&CellId(c as u64)].specific_h_j_kg,
                dt * power.cell_heat_w[c],
            );
            let sampled = &epoch.porous_samples[c].material;
            assert_eq!(sampled.material().card_identity(), skeleton.content_hash());
            assert_eq!(sampled.permeability().principal_m2(), &principal);
            assert_eq!(sampled.permeability().material_axes_in_world(), &axes);
            let endpoint = &accepted.endpoint_porous_samples()[c];
            near(
                endpoint.material.material().query().axes()["T"],
                next.temperature_k(endpoint.cell).unwrap(),
            );
        }
        speeds.push(a);
    }
    assert!(speeds[1] > speeds[0]);
}

#[test]
fn isotropic_source_preserves_the_original_scalar_native_solution() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let skeleton = skeleton([0.5; 3], false, 810.0);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let (_, accepted, _) = old
        .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 100.0,
    };
    let (flow, power) = creeping_flow_with_darcy_heat(
        &domain,
        &fluid,
        &[100.0; 4],
        &[200.0; 4],
        &config,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(&accepted.applied().unwrap().flow, &flow);
    assert_eq!(
        accepted.applied().unwrap().viscous_power.as_ref(),
        Some(&power)
    );

    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
        .unwrap();
    let (next, accepted, _) = old
        .advance_observed(&observer, 0.01, thermal_step(), || Ok(()))
        .unwrap();
    let original = simple_flow_with_material_resistance(
        &domain,
        &fluid,
        &[100.0; 4],
        &[200.0; 4],
        None,
        &config,
        &gate,
    )
    .unwrap();
    assert_eq!(accepted.applied().unwrap().flow, original);
    assert!(accepted.applied().unwrap().viscous_power.is_none());
    near(next.energy_j(), old.energy_j());
}

#[test]
fn actual_temperature_changes_both_fluid_viscosity_and_pinned_skeleton_permeability() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut speeds = Vec::new();
    for (thinning, changing_k) in [(false, false), (true, false), (false, true), (true, true)] {
        let card = source(100.0, thinning);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let skeleton = skeleton([0.5, 1.0, 2.0], changing_k, 810.0);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
            .unwrap()
            .with_unsteady_body_force_heating(forces(100.0), controls())
            .unwrap();
        let dt = 0.5;
        let (warm, first, _) = old
            .advance_observed(&observer, dt, thermal_step(), || Ok(()))
            .unwrap();
        let a0 = amplitude(0.0, dt, 100.0, 300.0, 100.0);
        let heat0 = dt * (2400.0 + 600.0) * a0 * a0;
        let t1 = 300.0 + heat0 / 400.0;
        near(warm.energy_j() - old.energy_j(), heat0);
        let mu1 = if thinning {
            100.0 * (1.0 - 0.09 * (t1 - 300.0))
        } else {
            100.0
        };
        let k_factor = if changing_k {
            1.0 + 0.1 * (t1 - 300.0)
        } else {
            1.0
        };
        let beta_xy = 3.0 * mu1 / k_factor;
        let a1 = amplitude(a0, dt, mu1, beta_xy, 100.0);
        let (next, second, thermal) = warm
            .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
            .unwrap();
        let step = second.unsteady_mechanical_step().unwrap();
        let audit = second.unsteady_heating_audit().unwrap();
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
        near(step.power.darcy_heat_w, 2.0 * beta_xy * a1 * a1);
        near(step.power.molecular_heat_w, 24.0 * mu1 * a1 * a1);
        near(audit.generated_heat_j, dt * step.power.heat_w);
        near(audit.backward_euler_loss_j, 2.0 * (a1 - a0).powi(2));
        near(audit.initial_kinetic_energy_j, 2.0 * a0 * a0);
        near(audit.final_kinetic_energy_j, 2.0 * a1 * a1);
        near(
            audit.supplied_body_work_j,
            audit.generated_heat_j + audit.final_kinetic_energy_j - audit.initial_kinetic_energy_j
                + audit.backward_euler_loss_j,
        );
        near(thermal.external_heat_j, audit.generated_heat_j);
        near(next.energy_j() - warm.energy_j(), audit.generated_heat_j);
        for sample in &second.applied().unwrap().samples {
            near(sample.temperature_k, t1);
            near(sample.transport.dynamic_viscosity_pa_s(), mu1);
        }
        for sample in &second.applied().unwrap().porous_samples {
            near(sample.material.material().query().axes()["T"], t1);
            near(
                sample.material.permeability().principal_m2()[0],
                0.5 * k_factor,
            );
        }
        for row in second
            .applied()
            .unwrap()
            .darcy_resistance_pa_s_m2
            .as_ref()
            .unwrap()
        {
            near(row[0][0], 2.0 * mu1 / k_factor);
            near(row[1][1], mu1 / k_factor);
            near(row[2][2], 0.5 * mu1 / k_factor);
        }
        assert_eq!(
            observer.unsteady_heating_audit(),
            Some(UnsteadyHeatingAudit::default())
        );
        speeds.push(a1);
    }
    assert!(speeds[1] > speeds[0]);
    assert!(speeds[2] > speeds[0]);
    assert!(speeds[3] > speeds[1] && speeds[3] > speeds[2]);
}

#[test]
fn sealed_axes_and_endpoint_source_exits_publish_no_candidate() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    {
        let skeleton = skeleton([0.0, 1.0, 2.0], false, 810.0);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap();
        assert!(
            observer
                .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
                .is_err()
        );
    }
    let skeleton = skeleton([0.5, 1.0, 2.0], false, 300.01);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_porous_materials(vec![Some(permeability(&skeleton, GRID_ALIGNED_AXES)); 4])
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let source_before = observer.endpoint_porous_samples().to_vec();
    let flow_before = observer.unsteady_flow().unwrap().clone();
    // The start and native solve are admitted; generated heat exits only the
    // original skeleton's endpoint support, while the liquid remains valid.
    assert!(matches!(
        old.advance_observed(&observer, 0.5, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Material { .. })
    ));
    assert_eq!(observer.endpoint_porous_samples(), source_before);
    assert_eq!(observer.unsteady_flow(), Some(&flow_before));
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
    assert!(observer.applied().is_none());
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
    assert_eq!(observer.unsteady_flow(), Some(&flow_before));
    assert_eq!(observer.endpoint_porous_samples(), source_before);
}
