//! Original-owner pressure budgets, calibration limits and conserved gas flow.
use super::*;

fn calibration(beta: f64) -> ForchheimerCalibration {
    ForchheimerCalibration::try_new(
        beta,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
        "synthetic calibrated beta [1/m], superficial-velocity convention",
    )
    .unwrap()
}

fn inertial_barrier(k: f64, thickness: f64, beta: f64) -> DarcyGasBarrier {
    barrier(k, thickness)
        .with_forchheimer_calibrations(&[calibration(beta)])
        .unwrap()
}

#[test]
fn zero_beta_retains_exact_original_darcy_and_slip_results() {
    let original = barrier(1e-12, 0.01);
    let zero = original
        .with_forchheimer_calibrations(&[calibration(0.0)])
        .unwrap();
    for (left, right) in [
        (200000.0, 100000.0),
        (100000.0, 200000.0),
        (150000.0, 150000.0),
    ] {
        assert_eq!(
            zero.mass_flow_profile(&gas(), left, right, 0.02).unwrap(),
            original
                .mass_flow_profile(&gas(), left, right, 0.02)
                .unwrap()
        );
    }
    assert_eq!(
        zero.linear_volume_resistance(&gas(), 0.02).unwrap(),
        original.linear_volume_resistance(&gas(), 0.02).unwrap()
    );
    let m = model(vec![0.75, 0.25], 0.0);
    let slip = manufactured_laminate(&m);
    let zero = slip
        .with_forchheimer_calibrations(&vec![calibration(0.0); slip.layers().len()])
        .unwrap();
    let state = m.evaluate(&point(300.0, 150000.0)).unwrap();
    assert_eq!(
        slip.mass_flow_profile_for_model(state.state(), &m, 240000.0, 80000.0, 0.01)
            .unwrap(),
        zero.mass_flow_profile_for_model(state.state(), &m, 240000.0, 80000.0, 0.01)
            .unwrap()
    );
}

#[test]
fn manufactured_layered_inertial_gas_closes_every_pressure_and_mass_budget() {
    let g = gas();
    let rt = g.specific_gas_constant * g.temperature;
    let properties = [(1e-10, 0.02, 2e6), (4e-10, 0.01, 0.0), (2e-10, 0.03, 5e5)];
    let mut layers = Vec::new();
    let mut pressures = [100000.0_f64; 4];
    let flux = 1.1;
    for (i, &(k, thickness, beta)) in properties.iter().enumerate().rev() {
        let loss = thickness * (g.dynamic_viscosity * flux / k + beta * flux * flux);
        pressures[i] = (pressures[i + 1] * pressures[i + 1] + 2.0 * rt * loss).sqrt();
    }
    for &(k, thickness, beta) in &properties {
        layers.push(
            layer(k, thickness)
                .with_forchheimer_calibration(calibration(beta))
                .unwrap(),
        );
    }
    let b = DarcyGasBarrier::try_new(layers, [1.0, 0.0, 0.0], support()).unwrap();
    let profile = b
        .mass_flow_profile(&g, pressures[0], pressures[3], 0.015)
        .unwrap();
    let flow = profile.flow;
    near(flow.mass_flux_left_to_right_kg_m2_s, flux);
    assert_eq!(profile.iterations, 0);
    for (&actual, &expected) in profile
        .interface_pressures_pa
        .as_ref()
        .unwrap()
        .iter()
        .zip(&pressures)
    {
        near(actual, expected);
    }
    near(
        flow.left_volume_rate_m3_s * pressures[0] / rt,
        flow.mass_rate_left_to_right_kg_s,
    );
    near(
        flow.right_volume_rate_m3_s * pressures[3] / rt,
        flow.mass_rate_left_to_right_kg_s,
    );
    near(
        flow.drag_dissipation_w,
        flow.mass_rate_left_to_right_kg_s * rt * (pressures[0] / pressures[3]).ln(),
    );
    let expected_conductance = 0.015
        / (g.dynamic_viscosity * b.geometric_resistance_per_m().unwrap()
            + b.inertial_geometric_resistance() * flux);
    near(flow.mean_volume_conductance_m3_pa_s, expected_conductance);
    let reverse = b
        .mass_flow_profile(&g, pressures[3], pressures[0], 0.015)
        .unwrap();
    assert_eq!(
        reverse.flow.mass_rate_left_to_right_kg_s,
        -flow.mass_rate_left_to_right_kg_s
    );
    assert_eq!(reverse.flow.drag_dissipation_w, flow.drag_dissipation_w);
    let faces = reverse.interface_pressures_pa.unwrap();
    for (i, &(k, thickness, beta)) in properties.iter().enumerate() {
        let pressure_budget =
            (faces[i + 1] - faces[i]) * (0.5 * faces[i + 1] + 0.5 * faces[i]) / rt;
        near(
            pressure_budget,
            thickness * (g.dynamic_viscosity * flux / k + beta * flux * flux),
        );
    }
}

#[test]
fn actual_inertial_speed_is_admitted_without_using_the_excessive_darcy_speed() {
    let domain = PorousGasSupport {
        maximum_superficial_speed_m_s: 1.0,
        ..support()
    };
    let original =
        DarcyGasBarrier::try_new(vec![layer(1e-8, 0.01)], [1.0, 0.0, 0.0], domain).unwrap();
    let law = ForchheimerCalibration::try_new(
        1e7,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        domain,
        "synthetic inertial support",
    )
    .unwrap();
    let b = original.with_forchheimer_calibrations(&[law]).unwrap();
    assert!(
        original
            .mass_flow(&gas(), 150100.0, 149900.0, 0.01)
            .is_err()
    );
    let flow = b.mass_flow(&gas(), 150100.0, 149900.0, 0.01).unwrap();
    assert!(flow.right_volume_rate_m3_s / 0.01 < 1.0);
    let doubled = b.mass_flow(&gas(), 150100.0, 149900.0, 0.02).unwrap();
    assert_eq!(
        doubled.mass_flux_left_to_right_kg_m2_s,
        flow.mass_flux_left_to_right_kg_m2_s
    );
    assert_eq!(
        doubled.mass_rate_left_to_right_kg_s,
        2.0 * flow.mass_rate_left_to_right_kg_s
    );
    assert!(b.linear_volume_resistance(&gas(), 0.01).is_err());

    // Synthetic numerical-range case: actual inertia regularizes the flow
    // even though the unopposed viscous conductance cannot be represented.
    let extreme = inertial_barrier(1e306, 0.01, 1e6);
    let flow = extreme.mass_flow(&gas(), 150100.0, 149900.0, 0.01).unwrap();
    assert!(flow.mean_volume_conductance_m3_pa_s.is_finite());
    assert!(flow.mean_volume_conductance_m3_pa_s > 0.0);
    near(
        flow.mean_volume_conductance_m3_pa_s,
        0.01 / (gas().dynamic_viscosity * extreme.geometric_resistance_per_m().unwrap()
            + extreme.inertial_geometric_resistance() * flow.mass_flux_left_to_right_kg_m2_s),
    );
}

#[test]
fn source_and_calibration_support_stay_attached_and_intersect() {
    let source = material(INTRINSIC_PERMEABILITY_DIMS);
    let original = PorousGasLayer::from_material_isotropic(
        source.clone(),
        "intrinsic_permeability",
        0.01,
        None,
    )
    .unwrap();
    let domain = PorousGasSupport {
        temperature_k: (290.0, 320.0),
        pressure_pa: (120000.0, 300000.0),
        maximum_superficial_speed_m_s: 0.1,
    };
    let law = ForchheimerCalibration::try_new(
        2e6,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        domain,
        "explicit synthetic calibration source",
    )
    .unwrap();
    let calibrated = original.with_forchheimer_calibration(law.clone()).unwrap();
    assert_eq!(calibrated.source(), Some(&source));
    assert_eq!(calibrated.permeability(), original.permeability());
    assert_eq!(calibrated.forchheimer_calibration(), Some(&law));
    let b = DarcyGasBarrier::try_new(vec![calibrated], [1.0, 0.0, 0.0], support()).unwrap();
    assert_eq!(b.support().temperature_k, (290.0, 320.0));
    assert_eq!(b.support().pressure_pa, (120000.0, 200000.0));
    assert_eq!(b.support().maximum_superficial_speed_m_s, 0.1);
    assert!(b.mass_flow(&gas(), 201000.0, 150000.0, 0.01).is_err());
    assert!(b.mass_flow(&gas(), 150000.0, 119000.0, 0.01).is_err());
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 325.0, 150000.0).unwrap();
    assert!(b.mass_flow(&hot, 151000.0, 149000.0, 0.01).is_err());
    assert!(b.with_forchheimer_calibrations(&[]).is_err());
}

#[test]
fn confined_anisotropy_is_exact_and_unqualified_directions_refuse() {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let tensor = IntrinsicPermeability::orthotropic(
        [1e-10, 9e-10, 1e-10],
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let original_layer =
        PorousGasLayer::try_new(tensor, 0.01, 0.0, "synthetic rotated pore tensor").unwrap();
    let original =
        DarcyGasBarrier::try_new(vec![original_layer], [1.0, 0.0, 0.0], support()).unwrap();
    assert!(
        original
            .with_forchheimer_calibrations(&[calibration(1e6)])
            .is_err()
    );
    let confined = original
        .with_flow_constraint(DarcyFlowConstraint::ConfinedNormalFlow)
        .unwrap();
    let law = ForchheimerCalibration::try_new(
        1e6,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::ConfinedNormalFlow,
        support(),
        "synthetic confined calibration",
    )
    .unwrap();
    let confined = confined.with_forchheimer_calibrations(&[law]).unwrap();
    near(
        confined.geometric_resistance_per_m().unwrap(),
        0.01 / (1.8e-10),
    );
    assert!(confined.mass_flow(&gas(), 160000.0, 140000.0, 0.01).is_ok());
    assert!(
        confined
            .with_flow_constraint(DarcyFlowConstraint::NormalPressureGradient)
            .is_err()
    );
    let wrong = ForchheimerCalibration::try_new(
        1e6,
        [0.0, 1.0, 0.0],
        DarcyFlowConstraint::ConfinedNormalFlow,
        support(),
        "different specimen direction",
    )
    .unwrap();
    assert!(confined.with_forchheimer_calibrations(&[wrong]).is_err());
}

#[test]
fn frozen_coefficients_share_continuum_geometry_and_preserve_seals() {
    let b = inertial_barrier(1e-10, 0.02, 2e6);
    let g = gas();
    let geometry = b.continuum_geometry().unwrap();
    assert_eq!(geometry.inertial_geometric_resistance(), 40000.0);
    let (linear, quadratic) = b.volume_resistance_coefficients(&g, 0.01).unwrap().unwrap();
    near(
        linear,
        g.dynamic_viscosity * geometry.geometric_resistance_per_m().unwrap() / 0.01,
    );
    near(quadratic, g.density * 40000.0 / 0.01_f64.powi(2));
    let equilibrium = b.mass_flow_profile(&g, 150000.0, 150000.0, 0.01).unwrap();
    near(
        equilibrium.flow.mean_volume_conductance_m3_pa_s,
        1.0 / linear,
    );
    assert_eq!(equilibrium.interface_pressures_pa, Some(vec![150000.0; 2]));
    let sealed = inertial_barrier(0.0, 0.02, 2e6);
    let flow = sealed
        .mass_flow_profile(&g, 200000.0, 100000.0, 0.01)
        .unwrap();
    assert_eq!(flow.flow.mass_rate_left_to_right_kg_s, 0.0);
    assert_eq!(flow.interface_pressures_pa, None);
    assert_eq!(
        sealed.volume_resistance_coefficients(&g, 0.01).unwrap(),
        None
    );
    assert_eq!(b.volume_resistance_coefficients(&g, 0.0).unwrap(), None);
    assert_eq!(
        b.mass_flow_profile(&g, 200000.0, 100000.0, 0.0)
            .unwrap()
            .interface_pressures_pa,
        None
    );
    let mut incoherent = g;
    incoherent.density *= 2.0;
    assert!(b.volume_resistance_coefficients(&incoherent, 0.01).is_err());
}

#[test]
fn invalid_coefficients_and_unqualified_combined_laws_refuse() {
    for beta in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ForchheimerCalibration::try_new(
                beta,
                [1.0, 0.0, 0.0],
                DarcyFlowConstraint::NormalPressureGradient,
                support(),
                "invalid beta"
            )
            .is_err()
        );
    }
    assert!(
        ForchheimerCalibration::try_new(
            1.0,
            [1.0, 0.0, 0.0],
            DarcyFlowConstraint::NormalPressureGradient,
            support(),
            " "
        )
        .is_err()
    );
    assert!(
        slip_layer(1e-12, 0.01, 1000.0)
            .with_forchheimer_calibration(calibration(1.0))
            .is_err()
    );
    use crate::porous::wet::{GasDiffusivityLaw, WetPoreGasLaw, WetPoreSupport};
    let wet = WetPoreGasLaw::try_new(
        0.0,
        0.0,
        0.5,
        GasDiffusivityLaw::ConstantTransportTortuosity { factor: 1.0 },
        WetPoreSupport {
            liquid_saturation: (0.0, 1.0),
            total_porosity: (0.5, 0.5),
            mean_gas_pore_radius_m: (1e-5, 1e-5),
        },
        "synthetic wet law",
    )
    .unwrap();
    let state = wet.sample(0.5, 0.5, 1e-5).unwrap();
    assert!(
        layer(1e-12, 0.01)
            .with_wet_state(state.clone())
            .unwrap()
            .with_forchheimer_calibration(calibration(1.0))
            .is_err()
    );
    assert!(
        inertial_barrier(1e-12, 0.01, 1.0)
            .with_wet_states(&[state])
            .is_err()
    );
    use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::dusty_gas::{
        DustyGasMedium, DustyGasPores,
    };
    let pores = DustyGasPores::try_new(1e-5, 0.5, 1.0, "synthetic dusty pores").unwrap();
    assert!(DustyGasMedium::try_new(inertial_barrier(1e-12, 0.01, 1.0), pores).is_err());
}
