//! Independent native strain, temperature/rate feedback and complete heat flow.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::environment::{
    FiniteGasContact, advance_finite_gas_with_material_flow,
};
use fs_lbm::conjugate::enthalpy::advection::liquid::rheology::LiquidFlowRheology;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_qty::Frequency;

#[path = "../../nsf3_thermochemical/source.rs"]
mod gas_source;

fn point_at(t: f64) -> QueryPoint {
    point(t, 100000.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap()
}

fn rheological_card(thermal_slope: f64, rate_slope: f64, max_rate: f64) -> MaterialCard {
    let original = source(100.0, false);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let viscosity = claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY;
        if viscosity || claim.key.name() == THERMAL_CONDUCTIVITY_PROPERTY {
            let dims = claim.key.dims();
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [300.0, 0.0],
                    [310.0, 0.0],
                    [300.0, max_rate],
                    [310.0, max_rate],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        if viscosity {
                            100.0 + thermal_slope * (p[0] - 300.0) + rate_slope * p[1]
                        } else {
                            0.01 + 0.05 * p[1]
                        },
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claim.validity = claim.validity.with_quantity(
                "rate",
                QuantitySpec::dimensional(Frequency::DIMS),
                0.0,
                max_rate,
            );
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &point_at(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &point_at(300.0),
            upper: &point_at(310.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-8,
    )
    .unwrap()
}

fn force() -> Vec<[f64; 3]> {
    vec![
        [-100.0, 100.0, 0.0],
        [-100.0, -100.0, 0.0],
        [100.0, 100.0, 0.0],
        [100.0, -100.0, 0.0],
    ]
}

fn closed() -> SimpleConfig {
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    config
}

#[test]
fn source_temperature_and_native_symmetric_strain_drive_the_same_momentum_and_heat() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut outcomes = Vec::new();
    for (thermal_slope, rate_slope) in [(0.0, 0.0), (-5.0, 20.0), (-5.0, -20.0)] {
        let card = rheological_card(thermal_slope, rate_slope, 2.0);
        let liquid = law(&card);
        let caloric = caloric(&card);
        let input = input(&ids, &liquid, &caloric);
        let old = initial(&domain, 1.0, &caloric);
        // The frozen adapter must not silently admit the same evolving table.
        assert!(
            input
                .feedback(&old, &domain, &config, (), &gate, 10)
                .is_err()
        );
        let source = LiquidFlowRheology::try_new(input, "rate").unwrap();
        let flow = source
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_unsteady_body_force_heating(force(), controls())
            .unwrap();
        assert_eq!(flow.endpoint_shear_rate_per_s().unwrap(), &[0.0; 4]);
        let dt = 0.5;
        let (warm, first, thermal) = old
            .advance_observed(&flow, dt, thermal_step(), || Ok(()))
            .unwrap();
        // Four native face masses, full symmetric stress, and fixed half-cell
        // walls: q1=400/(4/dt+24mu), gamma_rms^2=6q1^2 in each unit cell.
        let q1 = 400.0 / (4.0 / dt + 2400.0);
        let gamma1 = 6.0_f64.sqrt() * q1;
        let heat1 = dt * 2400.0 * q1 * q1;
        let t1 = 300.0 + heat1 / 400.0;
        let mu1 = 100.0 + thermal_slope * (t1 - 300.0) + rate_slope * gamma1;
        near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q1);
        near(thermal.external_heat_j, heat1);
        for (index, sample) in first.endpoint().iter().enumerate() {
            near(first.endpoint_shear_rate_per_s().unwrap()[index], gamma1);
            near(sample.temperature_k, t1);
            near(sample.material.dynamic_viscosity_pa_s(), mu1);
            near(
                sample.material.thermal_conductivity_w_m_k(),
                0.01 + 0.05 * gamma1,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        let (next, second, audit) = warm
            .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
            .unwrap();
        let q2 = (4.0 * q1 / dt + 400.0) / (4.0 / dt + 24.0 * mu1);
        let heat2 = dt * 24.0 * mu1 * q2 * q2;
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q2);
        near(audit.external_heat_j, heat2);
        near(next.energy_j() - warm.energy_j(), heat2);
        for rate in second.applied().unwrap().shear_rate_per_s.as_ref().unwrap() {
            near(*rate, gamma1);
        }
        let mechanical = second.unsteady_mechanical_step().unwrap();
        near(mechanical.power.molecular_heat_w, 24.0 * mu1 * q2 * q2);
        near(mechanical.final_kinetic_energy_j, 2.0 * q2 * q2);
        near(mechanical.backward_euler_loss_j, 2.0 * (q2 - q1).powi(2));
        near(
            dt * mechanical.power.supplied_body_power_w,
            heat2 + mechanical.final_kinetic_energy_j - mechanical.initial_kinetic_energy_j
                + mechanical.backward_euler_loss_j,
        );
        outcomes.push((q2, heat2));
        assert_eq!(flow.unsteady_time_s(), Some(0.0));
    }
    assert!(outcomes[1].0 < outcomes[0].0 && outcomes[2].0 > outcomes[0].0);
    assert!((outcomes[1].1 - outcomes[2].1).abs() > 0.1);
}

#[test]
fn graded_native_rates_include_normal_reciprocal_shear_and_half_cell_walls() {
    let x = [1.0, 2.0];
    let y = [0.5, 1.0];
    let depth = 0.75;
    let domain = VoxelDomain::graded([x.to_vec(), y.to_vec(), vec![depth]]).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let card = rheological_card(-5.0, 20.0, 4.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let old = initial(&domain, 1.0, &caloric);
    let flow = source
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_viscous_heating(controls())
        .unwrap();
    let (_, accepted, _) = old
        .advance_observed(&flow, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    let q = accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0);
    let lx = 0.5 * (x[0] + x[1]);
    let ly = 0.5 * (y[0] + y[1]);
    let shear = q * ((1.0 / y[0] + 1.0 / y[1]) / ly - (1.0 / x[0] + 1.0 / x[1]) / lx) / depth;
    assert!(shear.abs() > 1e-8);
    let edge_integral = lx * ly * depth * shear * shear;
    let area = (x[0] + x[1]) * (y[0] + y[1]);
    for c in 0..4 {
        let (i, j) = (c % 2, c / 2);
        let volume = x[i] * y[j] * depth;
        let normal = 4.0 * q * q / volume;
        let edge = edge_integral * x[i] * y[j] / area;
        let u = if j == 0 { -q } else { q } / (y[j] * depth);
        let wall = if j == 0 { 0.0 } else { 1.0 };
        let wx = depth * x[i] / y[j] * (wall - u).powi(2);
        let wy = depth * y[j] / x[i] * (q / (x[i] * depth)).powi(2);
        near(
            accepted.endpoint_shear_rate_per_s().unwrap()[c],
            ((normal + edge + wx + wy) / volume).sqrt(),
        );
    }
}

#[test]
fn constant_newtonian_data_preserves_native_flow_and_heat_exactly() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let card = source(100.0, false);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let input = input(&ids, &liquid, &caloric);
    let old = initial(&domain, 1.0, &caloric);
    let source = LiquidFlowRheology::try_new(input, "rate").unwrap();
    let rheological = source
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    let standard = input
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    let (a, a_history, a_audit) = old
        .advance_observed(&rheological, 1.0, thermal_step(), || Ok(()))
        .unwrap();
    let (b, b_history, b_audit) = old
        .advance_observed(&standard, 1.0, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(a_history.unsteady_flow(), b_history.unsteady_flow());
    assert_eq!(
        a_history.unsteady_heating_audit(),
        b_history.unsteady_heating_audit()
    );
    assert_eq!(a_audit, b_audit);
}

#[test]
fn actual_endpoint_shear_outside_source_support_rolls_back_the_complete_candidate() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let card = rheological_card(-5.0, 20.0, 0.2);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let old = initial(&domain, 1.0, &caloric);
    let flow = source
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    let refused = old.advance_observed(&flow, 0.5, thermal_step(), || Ok(()));
    assert!(matches!(refused, Err(MaterialFlowError::Material { .. })));
    assert_eq!(old.energy_j(), 0.0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(flow.endpoint_shear_rate_per_s().unwrap(), &[0.0; 4]);
    assert_eq!(flow.unsteady_time_s(), Some(0.0));
    assert!(flow.unsteady_step().is_none());
    gate.request();
    assert!(matches!(
        old.advance_observed(&flow, 0.5, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
    ));
}

#[test]
fn finite_atmosphere_reuses_live_rheology_for_heat_storage_and_motion() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let card = rheological_card(-5.0, 20.0, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let old = initial(&domain, 1.0, &caloric);
    let flow = source
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    let gas = IsochoricGasParcel::new(
        gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(305.0, 100000.0),
        1.0,
    )
    .unwrap();
    let contacts: Vec<_> = (0..4)
        .map(|i| FiniteGasContact {
            cell: CellId(i),
            area_m2: 1.0,
            distances_m: [0.1; 2],
        })
        .collect();
    let gas_controls = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 690.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let (warm, first, gas, first_report) = advance_finite_gas_with_material_flow(
        &old,
        0.5,
        thermal_step(),
        &flow,
        &gas,
        &contacts,
        gas_controls,
    )
    .unwrap();
    assert!(first_report.heat_into_condensed_j > 0.0);
    let (next, second, _, report) = advance_finite_gas_with_material_flow(
        &warm,
        1.0,
        thermal_step(),
        &first,
        &gas,
        &contacts,
        gas_controls,
    )
    .unwrap();
    let epoch = second.applied().unwrap();
    for (i, sample) in epoch.samples.iter().enumerate() {
        let rate = epoch.shear_rate_per_s.as_ref().unwrap()[i];
        assert!(rate > 0.0);
        near(
            sample.material.dynamic_viscosity_pa_s(),
            100.0 - 5.0 * (sample.temperature_k - 300.0) + 20.0 * rate,
        );
        near(
            report.last_applied[i].condensed_conductivity_w_m_k,
            0.01 + 0.05 * rate,
        );
    }
    near(
        report.other_thermal_external_heat_j,
        second.unsteady_heating_audit().unwrap().generated_heat_j,
    );
    assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
    assert!(next.energy_j() > warm.energy_j());
    assert_eq!(second.unsteady_time_s(), Some(1.0));
}
