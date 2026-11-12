//! G1/G3: original wall-property laws with prescribed bulk mass flux.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::{
    film::GasMassFlowState as FilmMassFlowState,
    sphere::ForcedSphereGasSurface,
    tube::{ForcedTubeGasSurface, PrescribedMassFlowTubeGasSurface},
};
use fs_couple::finite_gas_thermal::{FiniteGasThermalReport, WallConvectionSolve};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::UniformEnthalpyStep;
use fs_material::thermochemical_gas::ThermochemicalGasModel;

const TUBE_D: f64 = 0.01;
const TUBE_L: f64 = 1.0;
const SPHERE_D: f64 = 0.1;
const SPHERE_STREAM_AREA: f64 = 0.005;

fn tube_area() -> f64 {
    core::f64::consts::PI * TUBE_D * TUBE_D / 4.0
}

fn tube(model: &ThermochemicalGasModel) -> ForcedTubeGasSurface<'_> {
    ForcedTubeGasSurface {
        model,
        area_m2: core::f64::consts::PI * TUBE_D * TUBE_L,
        diameter: Length::new(TUBE_D),
        heated_length: Length::new(TUBE_L),
        speed_m_s: 0.0, // Explicit mass flow replaces the authored speed.
    }
}

fn tube_flow(
    surface: ForcedTubeGasSurface<'_>,
    rate: f64,
    laminar: bool,
) -> PrescribedMassFlowTubeGasSurface<'_> {
    if laminar {
        surface
            .laminar()
            .with_prescribed_mass_flow(rate, tube_area())
    } else {
        surface.with_prescribed_mass_flow(rate, tube_area())
    }
    .unwrap()
}

fn sphere(model: &ThermochemicalGasModel) -> ForcedSphereGasSurface<'_> {
    ForcedSphereGasSurface {
        model,
        area_m2: 0.02,
        diameter: Length::new(SPHERE_D),
        speed_m_s: 0.0,
    }
}

// Independent formulas use mass flux directly; the adapter's velocity is not
// an input to the reference Reynolds number.
fn tube_reference(
    bulk: &impl FluidTransport,
    wall: &impl FluidTransport,
    rate: f64,
    laminar: bool,
) -> (f64, f64) {
    let p = bulk.transport_state().unwrap();
    let w = wall.transport_state().unwrap();
    let re = rate * TUBE_D / (tube_area() * p.dynamic_viscosity_pa_s());
    let ratio = p.dynamic_viscosity_pa_s() / w.dynamic_viscosity_pa_s();
    let nu = if laminar {
        1.86 * (re * p.prandtl() * TUBE_D / TUBE_L).cbrt() * ratio.powf(0.14)
    } else {
        0.027 * re.powf(0.8) * p.prandtl().cbrt() * ratio.powf(0.14)
    };
    (re, nu * p.thermal_conductivity_w_m_k() / TUBE_D)
}

fn sphere_reference(
    bulk: &impl FluidTransport,
    wall: &impl FluidTransport,
    rate: f64,
) -> (f64, f64) {
    let p = bulk.transport_state().unwrap();
    let w = wall.transport_state().unwrap();
    let re = rate * SPHERE_D / (SPHERE_STREAM_AREA * p.dynamic_viscosity_pa_s());
    let ratio = p.dynamic_viscosity_pa_s() / w.dynamic_viscosity_pa_s();
    let nu = 2.0
        + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0)) * p.prandtl().powf(0.4) * ratio.powf(0.25);
    (re, nu * p.thermal_conductivity_w_m_k() / SPHERE_D)
}

fn moving_gas(temperature: f64) -> GasPiston {
    GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(1.0, true),
            source::point(temperature, 100000.0),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap()
}

fn wall_solve() -> WallConvectionSolve {
    WallConvectionSolve {
        temperature_increment_bracket_k: [-5.0, 5.0],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

fn isochoric_balance<S>(
    old: &IsochoricGasParcel,
    thermal: UniformEnthalpyStepInput<'_>,
    next: &UniformEnthalpyStep<FiniteGasThermalReport<S>>,
) {
    let r = &next.report;
    let du = old.mass_kg()
        * (r.next_gas.state().specific_internal_energy_j_kg()
            - old.state().specific_internal_energy_j_kg());
    let dh = thermal.mass_kg
        * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
    near(du, r.solid_to_gas_heat_j + r.external_gas_heat_j, 1e-7);
    near(
        dh,
        thermal.internal_heat_j + r.radiation_heat_j - r.solid_to_gas_heat_j,
        1e-7,
    );
    near(r.combined_energy_residual_j, 0.0, 1e-7);
    assert!(r.entropy.convection_production_j_k >= 0.0);
    near(r.entropy.balance_residual_j_k, 0.0, 1e-12);
    assert_eq!(r.next_gas.mass_kg(), old.mass_kg());
}

#[test]
fn tube_mass_flux_keeps_both_original_cards_wall_sources_and_piston_energy() {
    let old = moving_gas(400.0);
    let saved = old.clone();
    let curve = source::body_curve();
    // The bank's initial mixture must not replace the actual pure-A state.
    let bank = source::model_transport(0.25, true);
    let authored = tube(&bank);
    let thermal = UniformEnthalpyStepInput {
        surface_area_m2: authored.area_m2,
        // The free piston's gas spring, sqrt(rho c^2 A^2 / (m V)) = 118.3 rad/s
        // (1 kg, 0.01 m2, 1 L of ideal gas at 1 bar, gamma 1.4), is admitted by fs-phs only for
        // rate*dt <= 0.04 (MIDPOINT_ONE_CENT_ANGULAR_STEP): dt <= 3.38e-4 s. 0.3 ms gives 0.0355.
        duration_s: 0.0003,
        ..body(&curve, 120000.0)
    };
    let piston_input = PistonInput {
        duration_s: thermal.duration_s,
        ..input()
    };
    let heat = HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        ..exchange()
    };
    for laminar in [true, false] {
        let initial_re = if laminar { 1600.0 } else { 20000.0 };
        let rate = initial_re
            * old
                .gas()
                .state()
                .transport_state()
                .unwrap()
                .dynamic_viscosity_pa_s()
            * tube_area()
            / TUBE_D;
        let flow = tube_flow(authored, rate, laminar);
        assert_eq!(flow.mass_flow_kg_s(), rate);
        assert_eq!(flow.flow_area_m2(), tube_area());
        let mut values = Vec::new();
        for (temperature, pressure) in [(400.0, 100000.0), (400.0, 200000.0), (450.0, 100000.0)] {
            let bulk = old
                .gas()
                .model()
                .evaluate(&source::point(temperature, pressure))
                .unwrap();
            let resolved = flow.resolve(&bulk, thermal.initial).unwrap();
            let sample = &resolved.source;
            // The pre-existing public film path names this same diagnostic type.
            let velocity: FilmMassFlowState = sample.prescribed_mass_flow.unwrap();
            let wall = sample.surface_gas.as_ref().unwrap();
            let evaluated = sample.interface.evaluation.as_ref().unwrap();
            let (re, h) = tube_reference(&bulk, wall, rate, laminar);
            assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
            near(
                velocity.speed_m_s * bulk.state().density * tube_area(),
                rate,
                2e-14 * rate,
            );
            assert!((wall.state().density / bulk.state().density - 1.0).abs() > 0.02);
            assert_eq!(wall.state().temperature, thermal.initial.temperature_k());
            assert_eq!(wall.state().pressure, pressure);
            assert_eq!(wall.source_packs(), bulk.source_packs());
            assert_eq!(wall.components().len(), 1);
            assert_eq!(
                wall.components()[0].card_identity(),
                bulk.components()[0].card_identity()
            );
            assert_eq!(sample.interface.gas, bulk);
            assert_eq!(evaluated.transport, bulk.transport_state().unwrap());
            assert_eq!(
                evaluated.nusselt.card().id,
                if laminar {
                    CorrelationId::SiederTateLaminar
                } else {
                    CorrelationId::SiederTateTurbulent
                }
            );
            near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
            near(evaluated.coefficient.value.value(), h, 2e-12 * h);
            near(resolved.conductance_w_k, authored.area_m2 * h, 2e-12 * h);
            values.push((resolved.conductance_w_k, velocity.speed_m_s));
        }
        near(values[0].0, values[1].0, 2e-12 * values[0].0);
        near(values[1].1, 0.5 * values[0].1, 2e-12 * values[0].1);
        assert!((values[2].0 - values[0].0).abs() > 1e-4);
        let fixed = ForcedTubeGasSurface {
            speed_m_s: values[0].1,
            ..authored
        };
        let fixed = if laminar {
            fixed.laminar().resolve(old.gas().state(), thermal.initial)
        } else {
            fixed.resolve(old.gas().state(), thermal.initial)
        }
        .unwrap();
        assert!(fixed.source.prescribed_mass_flow.is_none());
        near(fixed.conductance_w_k, values[0].0, 2e-12 * values[0].0);
        assert!(
            authored
                .with_prescribed_mass_flow(rate, 0.5 * tube_area())
                .is_err()
        );
        let invalid_rate = if laminar { rate * 2.0 } else { rate / 10.0 };
        assert!(
            tube_flow(authored, invalid_rate, laminar)
                .resolve(old.gas().state(), thermal.initial)
                .is_err()
        );

        let next = flow
            .finite_body_candidate(&old, piston_input, thermal, controls(), heat, || false)
            .unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        let bulk = &sample.interface.gas;
        let wall = sample.surface_gas.as_ref().unwrap();
        let velocity = sample.prescribed_mass_flow.unwrap();
        assert_eq!(bulk, &report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
        assert!(
            (velocity.bulk_density_kg_m3 / old.gas().state().state().density - 1.0).abs() > 1e-10
        );
        near(
            velocity.speed_m_s * bulk.state().density * tube_area(),
            rate,
            2e-14 * rate,
        );
        let (_, h) = tube_reference(bulk, wall, rate, laminar);
        near(report.conductance_w_k, authored.area_m2 * h, 2e-12 * h);
        near(
            report.exchanged_heat_j,
            thermal.duration_s
                * report.conductance_w_k
                * (report.thermal_midpoint.temperature_k() - bulk.state().temperature),
            heat.source_tolerance_j,
        );
        let du = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside = piston_input.exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + kinetic + outside, 0.0, 1e-7);
        near(dh + report.exchanged_heat_j, 0.0, 1e-7);
        assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());
        let wrong_bank = source::model_transport(1.0, false);
        let wrong = ForcedTubeGasSurface {
            model: &wrong_bank,
            ..authored
        };
        assert!(
            tube_flow(wrong, rate, laminar)
                .resolve(old.gas().state(), thermal.initial)
                .is_err()
        );
        let insulated = tube_flow(
            ForcedTubeGasSurface {
                area_m2: 0.0,
                ..wrong
            },
            rate,
            laminar,
        )
        .resolve(old.gas().state(), thermal.initial)
        .unwrap();
        assert_eq!(insulated.conductance_w_k, 0.0);
        assert!(insulated.source.surface_gas.is_none());
        assert!(insulated.source.interface.evaluation.is_none());
        assert_eq!(old, saved);
    }
}

#[test]
fn sphere_stream_area_and_bulk_density_preserve_whitaker_and_latent_piston_balance() {
    let old = moving_gas(600.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let authored = sphere(old.gas().model());
    let rate = 1000.0
        * old
            .gas()
            .state()
            .transport_state()
            .unwrap()
            .dynamic_viscosity_pa_s()
        * SPHERE_STREAM_AREA
        / SPHERE_D;
    let flow = authored
        .with_prescribed_mass_flow(rate, SPHERE_STREAM_AREA)
        .unwrap();
    assert_eq!(flow.mass_flow_kg_s(), rate);
    assert_eq!(flow.flow_area_m2(), SPHERE_STREAM_AREA);
    assert_ne!(
        SPHERE_STREAM_AREA,
        core::f64::consts::PI * SPHERE_D * SPHERE_D / 4.0
    );
    let thermal = UniformEnthalpyStepInput {
        // The free piston's gas spring, sqrt(rho c^2 A^2 / (m V)) = 118.3 rad/s
        // (1 kg, 0.01 m2, 1 L of ideal gas at 1 bar, gamma 1.4), is admitted by fs-phs only for
        // rate*dt <= 0.04 (MIDPOINT_ONE_CENT_ANGULAR_STEP): dt <= 3.38e-4 s. 0.3 ms gives 0.0355.
        duration_s: 0.0003,
        ..body(&curve, 180000.0)
    };
    let resolved = flow.resolve(old.gas().state(), thermal.initial).unwrap();
    let wall = resolved.source.surface_gas.as_ref().unwrap();
    let (re, h) = sphere_reference(old.gas().state(), wall, rate);
    let velocity = resolved.source.prescribed_mass_flow.unwrap();
    assert_eq!(
        velocity.bulk_density_kg_m3,
        old.gas().state().state().density
    );
    assert!((wall.state().density / velocity.bulk_density_kg_m3 - 1.0).abs() > 0.3);
    let evaluated = resolved.source.interface.evaluation.as_ref().unwrap();
    assert_eq!(
        evaluated.transport,
        old.gas().state().transport_state().unwrap()
    );
    near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
    near(evaluated.coefficient.value.value(), h, 2e-12 * h);
    let same_flux = authored
        .with_prescribed_mass_flow(2.0 * rate, 2.0 * SPHERE_STREAM_AREA)
        .unwrap()
        .resolve(old.gas().state(), thermal.initial)
        .unwrap();
    assert_eq!(same_flux.conductance_w_k, resolved.conductance_w_k);
    let fixed = ForcedSphereGasSurface {
        speed_m_s: velocity.speed_m_s,
        ..authored
    }
    .resolve(old.gas().state(), thermal.initial)
    .unwrap();
    assert_eq!(fixed.conductance_w_k, resolved.conductance_w_k);
    assert!(fixed.source.prescribed_mass_flow.is_none());
    assert!(
        authored
            .with_prescribed_mass_flow(0.0, SPHERE_STREAM_AREA)
            .is_err()
    );
    assert!(authored.with_prescribed_mass_flow(rate, f64::NAN).is_err());
    let cold = old
        .gas()
        .model()
        .evaluate(&source::point(400.0, 100000.0))
        .unwrap();
    assert!(flow.resolve(&cold, thermal.initial).is_err()); // No viscosity-ratio clamp.
    let wrong_bank = source::model_transport(1.0, false);
    let wrong = ForcedSphereGasSurface {
        model: &wrong_bank,
        ..authored
    };
    assert!(
        wrong
            .with_prescribed_mass_flow(rate, SPHERE_STREAM_AREA)
            .unwrap()
            .resolve(old.gas().state(), thermal.initial)
            .is_err()
    );

    let piston_input = PistonInput {
        duration_s: thermal.duration_s,
        ..input()
    };
    let heat = HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        ..exchange()
    };
    let next = flow
        .finite_body_candidate(&old, piston_input, thermal, controls(), heat, || false)
        .unwrap();
    let report = next.step();
    let sample = &next.thermal().report;
    let bulk = &sample.interface.gas;
    let wall = sample.surface_gas.as_ref().unwrap();
    let velocity = sample.prescribed_mass_flow.unwrap();
    assert_eq!(bulk, &report.piston.midpoint);
    assert_eq!(sample.interface.body, report.thermal_midpoint);
    assert_eq!(
        wall.state().temperature,
        report.thermal_midpoint.temperature_k()
    );
    assert_eq!(wall.source_packs(), bulk.source_packs());
    assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
    assert!((velocity.bulk_density_kg_m3 / old.gas().state().state().density - 1.0).abs() > 1e-10);
    near(
        velocity.speed_m_s * bulk.state().density * SPHERE_STREAM_AREA,
        rate,
        2e-14 * rate,
    );
    let (_, h) = sphere_reference(bulk, wall, rate);
    near(report.conductance_w_k, authored.area_m2 * h, 2e-12 * h);
    near(
        report.exchanged_heat_j,
        thermal.duration_s
            * report.conductance_w_k
            * (report.thermal_midpoint.temperature_k() - bulk.state().temperature),
        heat.source_tolerance_j,
    );
    assert!(report.exchanged_heat_j < 0.0);
    assert_eq!(next.thermal().state.temperature_k(), 450.0);
    assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
    let du = old.gas().mass_kg()
        * (next.piston().gas().state().specific_internal_energy_j_kg()
            - old.gas().state().specific_internal_energy_j_kg());
    let dh = thermal.mass_kg
        * (next.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
    let outside = piston_input.exterior_pressure_pa
        * (next.piston().gas().volume_m3() - old.gas().volume_m3());
    near(du + dh + kinetic + outside, 0.0, 1e-7);
    near(dh + report.exchanged_heat_j, 0.0, 1e-7);
    assert_eq!(next.piston().gas().mass_kg(), old.gas().mass_kg());
    let insulated = ForcedSphereGasSurface {
        area_m2: 0.0,
        ..wrong
    }
    .with_prescribed_mass_flow(rate, SPHERE_STREAM_AREA)
    .unwrap()
    .finite_body_candidate(&old, piston_input, thermal, controls(), heat, || false)
    .unwrap();
    assert_eq!(insulated.step().exchanged_heat_j, 0.0);
    assert_eq!(insulated.thermal().state, thermal.initial);
    assert!(insulated.thermal().report.surface_gas.is_none());
    assert!(insulated.thermal().report.interface.evaluation.is_none());
    assert_eq!(old, saved);
}

#[test]
fn isochoric_tube_mass_flux_retains_original_density_and_actual_wall_feedback() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(400.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let authored = tube(old.model());
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.02,
            surface_area_m2: authored.area_m2,
            ..body(&curve, 120000.0)
        };
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-4.0, 4.0];
        for laminar in [true, false] {
            let initial_re = if laminar { 1600.0 } else { 20000.0 };
            let rate = initial_re
                * old
                    .state()
                    .transport_state()
                    .unwrap()
                    .dynamic_viscosity_pa_s()
                * tube_area()
                / TUBE_D;
            let next = tube_flow(authored, rate, laminar)
                .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
                .unwrap();
            let report = &next.report;
            let feedback = &report.convection.source;
            let sample = &feedback.source;
            let bulk = report.next_gas.state();
            let wall = sample.surface_gas.as_ref().unwrap();
            let velocity = sample.prescribed_mass_flow.unwrap();
            assert_eq!(&sample.interface.gas, bulk);
            assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
            near(
                velocity.bulk_density_kg_m3,
                old.state().state().density,
                2e-12,
            );
            near(
                velocity.speed_m_s,
                rate / (old.state().state().density * tube_area()),
                2e-12 * velocity.speed_m_s,
            );
            assert!((wall.state().density / velocity.bulk_density_kg_m3 - 1.0).abs() > 0.02);
            assert!(bulk.state().temperature > 400.0);
            assert!(next.state.temperature_k() < 420.0);
            assert_eq!(wall.state().temperature, feedback.temperature_k);
            assert_eq!(wall.state().pressure, bulk.state().pressure);
            assert_eq!(wall.source_packs(), bulk.source_packs());
            near(
                feedback
                    .endpoint_source
                    .surface_gas
                    .as_ref()
                    .unwrap()
                    .state()
                    .temperature,
                next.state.temperature_k(),
                1e-12,
            );
            let (re, h) = tube_reference(bulk, wall, rate, laminar);
            let evaluated = sample.interface.evaluation.as_ref().unwrap();
            assert_eq!(evaluated.transport, bulk.transport_state().unwrap());
            near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
            near(report.convection.coefficient_w_per_m2_k, h, 2e-12 * h);
            near(
                report.solid_to_gas_heat_j,
                thermal.duration_s
                    * thermal.surface_area_m2
                    * h
                    * (next.state.temperature_k() - bulk.state().temperature),
                environment.exchange.source_tolerance_j,
            );
            assert!(feedback.temperature_residual_k.abs() <= wall_solve().temperature_tolerance_k);
            isochoric_balance(&old, thermal, &next);
            assert_eq!(old, saved);
        }
    });
}

#[test]
fn isochoric_sphere_mass_flux_preserves_latent_heat_and_entropy() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(600.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.02,
            ..body(&curve, 180000.0)
        };
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-4.0, 4.0];
        let rate = 1000.0
            * old
                .state()
                .transport_state()
                .unwrap()
                .dynamic_viscosity_pa_s()
            * SPHERE_STREAM_AREA
            / SPHERE_D;
        let next = sphere(old.model())
            .with_prescribed_mass_flow(rate, SPHERE_STREAM_AREA)
            .unwrap()
            .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
            .unwrap();
        let report = &next.report;
        let feedback = &report.convection.source;
        let sample = &feedback.source;
        let bulk = report.next_gas.state();
        let wall = sample.surface_gas.as_ref().unwrap();
        let velocity = sample.prescribed_mass_flow.unwrap();
        assert_eq!(&sample.interface.gas, bulk);
        assert_eq!(velocity.bulk_density_kg_m3, bulk.state().density);
        near(
            velocity.bulk_density_kg_m3,
            old.state().state().density,
            2e-12,
        );
        near(
            velocity.speed_m_s,
            rate / (old.state().state().density * SPHERE_STREAM_AREA),
            2e-12 * velocity.speed_m_s,
        );
        assert!((wall.state().density / velocity.bulk_density_kg_m3 - 1.0).abs() > 0.3);
        assert!(bulk.state().temperature < 600.0);
        assert_eq!(wall.state().temperature, feedback.temperature_k);
        assert_eq!(wall.state().pressure, bulk.state().pressure);
        assert_eq!(wall.source_packs(), bulk.source_packs());
        near(
            feedback
                .endpoint_source
                .surface_gas
                .as_ref()
                .unwrap()
                .state()
                .temperature,
            next.state.temperature_k(),
            1e-12,
        );
        let (re, h) = sphere_reference(bulk, wall, rate);
        let evaluated = sample.interface.evaluation.as_ref().unwrap();
        assert_eq!(evaluated.transport, bulk.transport_state().unwrap());
        near(evaluated.nusselt.groups()["Re"], re, 2e-12 * re);
        near(report.convection.coefficient_w_per_m2_k, h, 2e-12 * h);
        near(
            report.solid_to_gas_heat_j,
            thermal.duration_s
                * thermal.surface_area_m2
                * h
                * (next.state.temperature_k() - bulk.state().temperature),
            environment.exchange.source_tolerance_j,
        );
        assert!(report.solid_to_gas_heat_j < 0.0);
        assert_eq!(next.state.temperature_k(), 450.0);
        assert!(next.state.liquid_mass_fraction() > 0.5);
        near(
            thermal.mass_kg * 60000.0 * (next.state.liquid_mass_fraction() - 0.5),
            -report.solid_to_gas_heat_j,
            1e-7,
        );
        isochoric_balance(&old, thermal, &next);
        assert_eq!(old, saved);
    });
}
