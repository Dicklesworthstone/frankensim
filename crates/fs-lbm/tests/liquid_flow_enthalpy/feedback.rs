//! Data-only liquid temperature feedback through the existing physical owners.
use super::*;
use fs_lbm::conjugate::FlowResistance;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;

fn feedback_card(scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.01,
        ),
    ] {
        claims.insert_claim(scalar(name, dims, value)).unwrap();
    }
    let mut viscosity = scalar(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, scale);
    viscosity.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, scale), (340.0, 0.2 * scale), (810.0, 0.1 * scale)],
        dims: DYNAMIC_VISCOSITY_DIMS,
    };
    viscosity.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(viscosity).unwrap();
    card(claims)
}

fn porous_config() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.resistances.push(FlowResistance::Volume {
        lo: [0, 0, 0],
        hi: [4, 1, 1],
        permeability_m2: [1.0; 3],
        inertial_per_m: [0.0; 3],
    });
    config.internal_fans.push(InternalFan {
        patch: FacePatch {
            axis: 0,
            index: 2,
            lo: [0, 0, 0],
            hi: [4, 1, 1],
        },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 4.0), (100.0, 4.0)]).unwrap(),
    });
    config
}

fn viscosity(scale: f64, h: f64) -> f64 {
    let t = 300.0 + h / 2.0;
    if t <= 340.0 {
        scale * (1.0 - 0.8 * (t - 300.0) / 40.0)
    } else {
        scale * (0.2 - 0.1 * (t - 340.0) / 470.0)
    }
}

#[test]
fn source_temperature_changes_actual_flow_each_substep_and_closes_enthalpy() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut outputs = Vec::new();
    for scale in [1.0, 2.0] {
        let card = feedback_card(scale);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
            .unwrap();
        let (next, accepted, audit) = old
            .advance_observed(&observer, 1.0, config(), || Ok(()))
            .unwrap();
        let mut h = [0.0; 4];
        let mut exported = 0.0;
        let mut mass_through = 0.0;
        let mut last_q = 0.0;
        for &dt in &accepted.inner().durations {
            // Independent Darcy series-resistance + donor-cell energy balance.
            // The native MAC solver owns the actual flow used by the test.
            let q = 4.0 / h.iter().map(|&h| viscosity(scale, h)).sum::<f64>();
            exported += dt * q * h[3];
            mass_through += dt * q;
            let old_h = h;
            for c in 0..4 {
                let upstream = if c == 0 { 0.0 } else { old_h[c - 1] };
                let heater = if c == 1 { 100.0 } else { 0.0 };
                h[c] += dt * (heater + q * (upstream - old_h[c]));
            }
            last_q = q;
        }
        assert!(accepted.momentum_solves() > 2);
        assert_eq!(accepted.momentum_solves(), accepted.inner().durations.len());
        assert_eq!(accepted.audit().substeps, accepted.inner().durations.len());
        let applied = accepted.applied().unwrap();
        near(applied.flow.field.flux_x(2, 0, 0), last_q);
        assert!(last_q > 1.0 / scale);
        assert!(applied.time_s < next.time_s());
        for (c, expected) in h.iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, *expected);
        }
        for sample in accepted.endpoint() {
            near(sample.specific_enthalpy_j_kg, h[sample.voxel]);
            near(sample.temperature_k, 300.0 + h[sample.voxel] / 2.0);
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                viscosity(scale, h[sample.voxel]),
            );
            let receipt = &sample
                .material
                .material()
                .property(DYNAMIC_VISCOSITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt;
            assert_eq!(
                receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
        near(accepted.audit().inlet_mass_kg, mass_through);
        near(accepted.audit().outlet_mass_kg, mass_through);
        near(accepted.audit().external_enthalpy_j, -exported);
        near(next.energy_j() - old.energy_j(), 100.0 - exported);
        near(
            audit.external_heat_j,
            accepted.inner().supplied_j - exported,
        );
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(observer.momentum_solves(), 0);
        assert!(observer.applied().is_none());
        assert_eq!(observer.inner(), &Heater::default());
        outputs.push((last_q, h[1]));
    }
    assert!(outputs[1].0 < outputs[0].0);
    assert!(outputs[1].1 > outputs[0].1);
}

#[test]
fn fixed_viscosity_feedback_matches_existing_frozen_owner() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let live = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
        .unwrap();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let flow = fluid.solve_flow(&domain, &config_flow, &gate).unwrap();
    let frozen = flow
        .advect_liquid_enthalpy(
            &old,
            input(&ids, &liquid, &caloric),
            Heater::default(),
            &gate,
        )
        .unwrap();
    let (a, a_observer, _) = old
        .advance_observed(&live, 1.0, config(), || Ok(()))
        .unwrap();
    let (b, b_observer, _) = old
        .advance_observed(&frozen, 1.0, config(), || Ok(()))
        .unwrap();
    for &id in ids.iter().flatten() {
        near(
            a.cells()[&id].specific_h_j_kg,
            b.cells()[&id].specific_h_j_kg,
        );
    }
    near(
        a_observer.audit().external_enthalpy_j,
        b_observer.advection().audit().external_enthalpy_j,
    );
}

#[test]
fn source_support_budget_inner_history_and_late_cancellation_are_atomic() {
    let card = card(claims(1.0, 1.0, true));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&observer, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Material { .. })
    ));
    let limited = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 1)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&limited, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::Budget(
            "live material momentum solves"
        )))
    ));
    let history = Heater {
        reject_step: Some(1),
        ..Heater::default()
    };
    let rejecting = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, history.clone(), &gate, 100)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&rejecting, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Inner(DepositionError::InvalidInput(_)))
    ));
    let mut checks = 0;
    let result = old.advance_observed(&observer, 1.0, config(), || {
        checks += 1;
        if checks == 30 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
            | Err(MaterialFlowError::Inner(DepositionError::Cancelled))
    ));
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    assert_eq!(observer.audit().substeps, 0);
    assert_eq!(observer.inner(), &Heater::default());
    assert_eq!(rejecting.inner(), &history);
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
    near(old.mass_kg(), 4.0);
}

#[test]
fn live_binding_keeps_original_caloric_context_and_borrowed_law() {
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let pressure = source_caloric(
        &card,
        200000.0,
        SPECIFIC_HEAT_CAPACITY_PROPERTY,
        DENSITY_PROPERTY,
        SolidLiquidPhase::Liquid,
    );
    let equal_curve = caloric.clone();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    for other in [&pressure, &equal_curve] {
        assert!(
            input(&ids, &liquid, other)
                .feedback(&old, &domain, &config_flow, (), &gate, 100,)
                .is_err()
        );
    }
}
