//! G0/G1/G3/G4: finite wall storage, molecular recoil, frame and atomic refusal.
//! Manufactured caloric/transport inputs do not validate a physical heat law.

mod exchange;
mod trajectory;

use super::*;
use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
use crate::finite_gas_thermal::{ConvectionCoefficient, ConvectionError};
use crate::rocket_blowdown::effusion::thermal::{
    EffusionWallCandidate, EffusionWallError, EffusionWallSolve,
};
use fs_blake3::ContentHash;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};

const WALL_MASS: f64 = 0.005;
const WALL_CP: f64 = 1000.0;
const WALL_AREA: f64 = 0.01;

fn curve(densities: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([41; 32]),
        SolidLiquidPhase::Solid,
        [200.0, 1800.0]
            .into_iter()
            .zip(densities)
            .map(|(temperature_k, density)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: WALL_CP * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: density,
            })
            .collect(),
    )
    .unwrap()
}

fn wall(curve: &EquilibriumEnthalpyPhaseCurve, temperature: f64) -> UniformEnthalpyStepInput<'_> {
    let initial = curve
        .state_at_specific_enthalpy(WALL_CP * (temperature - 200.0))
        .unwrap();
    UniformEnthalpyStepInput {
        curve,
        initial,
        mass_kg: WALL_MASS,
        volume_m3: WALL_MASS / initial.bulk_density_kg_m3(),
        surface_area_m2: WALL_AREA,
        internal_heat_j: 0.0,
        duration_s: DT,
    }
}

fn body(h: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: T0,
        radiation_temperature_k: T0,
        convection_w_per_m2_k: h,
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-10,
        maximum_thermal_residual_j: 1e-9,
    }
}

fn solve() -> EffusionWallSolve {
    EffusionWallSolve {
        gas: endpoint(),
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-50.0, 50.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
        maximum_gas_temperature_increment_k: 200.0,
        recoil: controls(),
        combined_energy_tolerance_j: 1e-6,
    }
}

fn free() -> EffusionRocketEnvironment {
    EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    }
}

fn scope(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x455754,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn check_total<S>(
    before: &GasRocketState,
    wall: UniformEnthalpyStepInput<'_>,
    candidate: &EffusionWallCandidate<S>,
) {
    let after = candidate.endpoint();
    let step = candidate.step();
    let wall_energy = wall.mass_kg
        * (candidate.wall().specific_enthalpy_j_kg() - wall.initial.specific_enthalpy_j_kg());
    near(
        internal(after) - internal(before) + wall_energy + kinetic(after) - kinetic(before)
            + step.recoil.exhaust_energy_world_j
            - step.recoil.exhaust.external_work_j
            - step.external_thermal_heat_j,
        0.0,
        3e-7,
    );
    vector_near(
        after
            .vehicle()
            .velocity_m_s
            .scale(mass(after))
            .sub(before.vehicle().velocity_m_s.scale(mass(before)))
            .add(step.recoil.exhaust_momentum_world_kg_m_s)
            .sub(step.recoil.exhaust.external_impulse_n_s),
        Vec3::ZERO,
        3e-12,
    );
    near(
        step.combined_energy_residual_j,
        0.0,
        solve().combined_energy_tolerance_j,
    );
    assert!(
        step.absolute_energy_residual_sum_j + step.energy_uncertainty_bound_j
            <= solve().combined_energy_tolerance_j
    );
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(
        after.reservoir().volume_m3(),
        before.reservoir().volume_m3()
    );
}

#[test]
fn finite_wall_heating_changes_actual_escape_and_recoil_with_one_closed_energy_account() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let law = interface(before.reservoir().model());
        let normal = Vec3::new(1.0, 0.0, 0.0);
        let outlets = [outlet(&law, normal)];
        let candidate = before
            .advance_effusion_with_wall(cx, free(), &outlets, input, &body(100.0), 0.0, solve())
            .unwrap();
        let heat_capacity_gas = before.reservoir().mass_kg() * 1.5 * R / 0.028;
        let heat_capacity_wall = WALL_MASS * WALL_CP;
        let conductance_time = DT * 100.0 * WALL_AREA;
        let q = conductance_time * (600.0 - T0)
            / (1.0 + conductance_time * (1.0 / heat_capacity_wall + 1.0 / heat_capacity_gas));
        let heated_t = T0 + q / heat_capacity_gas;
        near(candidate.step().thermal.report.solid_to_gas_heat_j, q, 3e-8);
        near(
            candidate
                .step()
                .thermal
                .report
                .next_gas
                .state()
                .state()
                .temperature,
            heated_t,
            2e-7,
        );
        near(
            candidate.wall().temperature_k(),
            600.0 - q / heat_capacity_wall,
            2e-8,
        );
        near(
            candidate.step().recoil.expelled_mass_kg,
            escaped_moles(0.028, 1.0) * 0.028 * (heated_t / T0).sqrt(),
            1e-14,
        );
        vector_near(
            candidate.step().recoil.relative_momentum_kg_m_s,
            normal.scale(AREA * P0 * DT * heated_t / (2.0 * T0)),
            3e-12,
        );
        assert!(candidate.endpoint().reservoir().state().state().temperature < heated_t);
        assert!(candidate.endpoint().vehicle().velocity_m_s.x < 0.0);
        assert_eq!(candidate.step().duration_s, DT);
        assert_eq!(candidate.step().recoil.duration_s, DT);
        assert_eq!(candidate.step().external_thermal_heat_j, 0.0);
        check_total(&before, input, &candidate);
        assert_eq!(before, saved);
    });
}

#[test]
fn insulated_wall_is_exact_original_recoil_and_closed_outlets_still_exchange_heat() {
    scope(|cx, _| {
        let before = initial(MIXTURE, Vec3::new(10.0, -5.0, 2.0));
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(0.6, 0.8, 0.0))];
        let direct = before
            .advance_effusion_with_environment(free(), &outlets, DT, endpoint(), controls(), || {
                false
            })
            .unwrap();
        let insulated = before
            .advance_effusion_with_wall(cx, free(), &outlets, input, &body(0.0), 0.0, solve())
            .unwrap();
        assert_eq!(insulated.endpoint(), direct.endpoint());
        assert_eq!(&insulated.step().recoil, direct.step());
        assert_eq!(*insulated.wall(), input.initial);
        let closed = before
            .advance_effusion_with_wall(cx, free(), &[], input, &body(100.0), 0.0, solve())
            .unwrap();
        assert_eq!(closed.step().recoil.expelled_mass_kg, 0.0);
        assert_eq!(
            closed.endpoint().reservoir().mass_kg(),
            before.reservoir().mass_kg()
        );
        assert!(closed.endpoint().reservoir().state().state().temperature > T0);
        assert!(closed.wall().temperature_k() < input.initial.temperature_k());
        check_total(&before, input, &closed);
    });
}

#[test]
fn gas_wall_radiation_and_mechanical_sources_are_counted_once() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::new(10.0, -5.0, 3.0));
        let chart = curve([1000.0; 2]);
        let mut input = wall(&chart, 600.0);
        input.internal_heat_j = 2.0;
        let mut boundary = body(100.0);
        boundary.transport = LumpedThermalTransport::try_declared(100.0, 0.2).unwrap();
        boundary.radiation_temperature_k = 650.0;
        let environment = EffusionRocketEnvironment {
            gravity_world_m_s2: Vec3::new(0.0, -9.81, 0.0),
            external_force_world_n: Vec3::new(0.02, 0.0, -0.01),
            ..free()
        };
        let law = interface(before.reservoir().model());
        let candidate = before
            .advance_effusion_with_wall(
                cx,
                environment,
                &[outlet(&law, Vec3::new(1.0, 0.0, 0.0))],
                input,
                &boundary,
                3.0,
                solve(),
            )
            .unwrap();
        assert!(candidate.step().thermal.report.radiation_heat_j > 0.0);
        near(
            candidate.step().external_thermal_heat_j,
            5.0 + candidate.step().thermal.report.radiation_heat_j,
            1e-14,
        );
        assert_ne!(candidate.step().recoil.exhaust.external_work_j, 0.0);
        check_total(&before, input, &candidate);
    });
}

#[test]
fn common_velocity_shift_preserves_wall_heat_and_relative_molecular_recoil() {
    scope(|cx, _| {
        let at_rest = initial(MIXTURE, Vec3::ZERO);
        let boost = Vec3::new(37.0, -19.0, 11.0);
        let moving = initial(MIXTURE, boost);
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let law = interface(at_rest.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(0.0, 1.0, 0.0))];
        let run = |state: &GasRocketState| {
            state
                .advance_effusion_with_wall(cx, free(), &outlets, input, &body(100.0), 0.0, solve())
                .unwrap()
        };
        let a = run(&at_rest);
        let b = run(&moving);
        assert_eq!(a.wall(), b.wall());
        assert_eq!(a.step().thermal, b.step().thermal);
        assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
        assert_eq!(
            a.step().recoil.relative_momentum_kg_m_s,
            b.step().recoil.relative_momentum_kg_m_s
        );
        vector_near(
            b.endpoint()
                .vehicle()
                .velocity_m_s
                .sub(a.endpoint().vehicle().velocity_m_s),
            boost,
            1e-12,
        );
        check_total(&moving, input, &b);
    });
}

#[test]
fn endpoint_convection_is_retained_and_late_escape_or_cancellation_is_atomic() {
    scope(|cx, gate| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let law_h = |gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState| {
            Ok::<_, &'static str>(ConvectionCoefficient {
                coefficient_w_per_m2_k: gas.state().temperature / 4.0,
                source: gas.state().temperature,
            })
        };
        let candidate = before
            .advance_effusion_with_wall_convection(
                cx,
                free(),
                &outlets,
                input,
                &body(0.0),
                0.0,
                solve(),
                law_h,
            )
            .unwrap();
        let thermal = &candidate.step().thermal.report;
        assert_eq!(
            thermal.convection.source,
            thermal.next_gas.state().state().temperature
        );
        assert_eq!(
            thermal.convection.coefficient_w_per_m2_k,
            thermal.convection.source / 4.0
        );
        assert!(thermal.solid_to_gas_heat_j > 0.0);
        check_total(&before, input, &candidate);

        let mut too_much = solve();
        too_much.recoil.transport.exchange.maximum_species_fraction = 1e-12;
        let mut calls = 0;
        let failed = before.advance_effusion_with_wall_convection(
            cx,
            free(),
            &outlets,
            input,
            &body(100.0),
            0.0,
            too_much,
            |gas| {
                calls += 1;
                law_h(gas)
            },
        );
        assert!(calls > 1);
        assert!(matches!(failed, Err(EffusionWallError::Effusion(_))));
        let rejected = before.advance_effusion_with_wall_convection(
            cx,
            free(),
            &outlets,
            input,
            &body(100.0),
            0.0,
            solve(),
            |_| Err::<ConvectionCoefficient<()>, _>("unsupported source"),
        );
        assert!(matches!(
            rejected,
            Err(EffusionWallError::Thermal(ConvectionError::Convection(
                "unsupported source"
            )))
        ));

        let cancelled = before.advance_effusion_with_wall_convection(
            cx,
            free(),
            &outlets,
            input,
            &body(100.0),
            0.0,
            solve(),
            |gas| {
                gate.request();
                law_h(gas)
            },
        );
        assert!(cancelled.is_err());
        assert_eq!(before, saved);
        assert_eq!(input.initial.temperature_k(), 600.0);
    });
}

#[test]
fn unmodeled_expansion_extra_wall_mass_and_inconsistent_geometry_are_refused() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        for invalid in [
            UniformEnthalpyStepInput {
                mass_kg: 0.02,
                ..input
            },
            UniformEnthalpyStepInput {
                volume_m3: input.volume_m3 * 2.0,
                ..input
            },
            UniformEnthalpyStepInput {
                duration_s: 0.0,
                ..input
            },
        ] {
            assert!(matches!(
                before.advance_effusion_with_wall(
                    cx,
                    free(),
                    &[],
                    invalid,
                    &body(100.0),
                    0.0,
                    solve(),
                ),
                Err(EffusionWallError::Invalid(_))
            ));
        }
        let expanding = curve([1000.0, 950.0]);
        assert!(matches!(
            before.advance_effusion_with_wall(
                cx,
                free(),
                &[],
                wall(&expanding, 600.0),
                &body(100.0),
                0.0,
                solve(),
            ),
            Err(EffusionWallError::Invalid(_))
        ));
        let mut invalid_budget = solve();
        invalid_budget.combined_energy_tolerance_j = f64::NAN;
        assert!(matches!(
            before.advance_effusion_with_wall(
                cx,
                free(),
                &[],
                input,
                &body(100.0),
                0.0,
                invalid_budget,
            ),
            Err(EffusionWallError::Invalid(_))
        ));
    });
}
