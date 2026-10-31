//! G0/G3/G4: real finite-water storage through the original body heat owner.
//! Gas packs and the body are manufactured; saturation uses the IAPWS source.

use fs_blake3::ContentHash;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, ThermalPistonError};
use fs_couple::finite_gas_thermal::{ConvectionCoefficient, condensation::*};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::QueryPoint;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
    condensation::{WaterCondensationModel, WetGasParcel},
};
use fs_material::thermochemical_gas::{
    energy::{GasEnergySolve, IsochoricGasParcel},
    packs::from_packs,
};
use fs_qty::{
    QuantitySpec,
    chemistry::SpeciesId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod source;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn wet_parcel(temperature_k: f64, water_fraction: f64, liquid_kg: f64) -> WetGasParcel {
    let carrier = source::packs("carrier", 0.028, 3.5, 200.0);
    let water = source::packs("water", 0.01801528, 4.0, 200.0);
    let water_id = SpeciesId::new("water").unwrap();
    let model = from_packs(
        &[carrier.binding(), water.binding()],
        Composition::new(
            CompositionBasis::MoleFraction,
            vec![1.0 - water_fraction, water_fraction],
        )
        .unwrap(),
        Some(&water_id),
    )
    .unwrap();
    let condensation = WaterCondensationModel::new(model.clone(), water_id, 300.0).unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            temperature_k,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
        .unwrap();
    WetGasParcel::new(
        condensation,
        IsochoricGasParcel::new(model, point, 0.001).unwrap(),
        liquid_kg,
    )
    .unwrap()
}
fn controls() -> WetGasHeatSolve {
    WetGasHeatSolve {
        temperature_bracket_k: [274.0, 323.0],
        energy_tolerance_j: 1e-9,
        maximum_iterations: 100,
    }
}
fn water_mass(parcel: &WetGasParcel) -> f64 {
    let gas = parcel.gas().state().state();
    parcel.liquid_mass_kg()
        + gas.water_mole_fraction * gas.pressure * parcel.gas().volume_m3() * 0.01801528
            / (R * gas.temperature)
}
fn carrier_moles(parcel: &WetGasParcel) -> f64 {
    let gas = parcel.gas().state().state();
    (1.0 - gas.water_mole_fraction) * gas.pressure * parcel.gas().volume_m3()
        / (R * gas.temperature)
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
                kernel_id: 0x434f4e44,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([41; 32]),
        SolidLiquidPhase::Solid,
        [270.0, 340.0]
            .map(|t| EnthalpyPhaseKnot {
                temperature_k: t,
                specific_enthalpy_j_kg: 1000.0 * (t - 270.0),
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap()
}
fn boundary() -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 999.0,
        radiation_temperature_k: 285.0,
        convection_w_per_m2_k: 100.0,
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-8,
        maximum_thermal_residual_j: 1e-8,
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-20.0, 0.0],
        source_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 100,
    }
}

#[test]
fn cooling_crosses_dew_point_and_reheating_evaporates_the_same_finite_water() {
    let warm = wet_parcel(310.0, 0.05, 0.0);
    let initial = warm.clone();
    let cold = warm.at_temperature(290.0, || false).unwrap();
    assert!(cold.condensed_mass_kg > 0.0);
    let heat = cold.total_internal_energy_change_j;
    assert!(heat < 0.0);
    let solved = heat_candidate(&warm, heat, controls(), || false).unwrap();
    near(solved.phase.parcel.temperature_k(), 290.0, 2e-8);
    near(water_mass(&solved.phase.parcel), water_mass(&warm), 1e-16);
    near(
        carrier_moles(&solved.phase.parcel),
        carrier_moles(&warm),
        1e-15,
    );
    assert!(solved.energy_residual_j.abs() <= controls().energy_tolerance_j);
    let mut frozen = warm.gas().clone();
    frozen
        .add_heat(
            heat,
            GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 100,
            },
            || false,
        )
        .unwrap();
    assert!(
        solved.phase.parcel.temperature_k() > frozen.state().state().temperature + 10.0,
        "latent storage must materially resist cooling"
    );

    let reheated = heat_candidate(&solved.phase.parcel, -heat, controls(), || false).unwrap();
    near(reheated.phase.parcel.temperature_k(), 310.0, 2e-8);
    assert_eq!(reheated.phase.parcel.liquid_mass_kg(), 0.0);
    near(water_mass(&reheated.phase.parcel), water_mass(&warm), 1e-16);
    near(
        reheated
            .phase
            .parcel
            .gas()
            .state()
            .state()
            .water_mole_fraction,
        0.05,
        1e-14,
    );
    assert_eq!(warm, initial);
    assert_eq!(
        solved,
        heat_candidate(&warm, heat, controls(), || false).unwrap()
    );
}

#[test]
fn cold_finite_body_collects_latent_heat_and_continuation_closes_all_three_stores() {
    scope(|cx, _| {
        let curve = curve();
        let boundary = boundary();
        let mut parcel = wet_parcel(310.0, 0.05, 0.0);
        let original_water = water_mass(&parcel);
        let original_carrier = carrier_moles(&parcel);
        let initial_energy = parcel.total_internal_energy_j();
        let mut state = curve.state_at_specific_enthalpy(15000.0).unwrap();
        let initial_h = state.specific_enthalpy_j_kg();
        let mut total_heat = 0.0;
        for _ in 0..8 {
            let input = UniformEnthalpyStepInput {
                curve: &curve,
                initial: state,
                mass_kg: 0.01,
                volume_m3: 1e-5,
                surface_area_m2: 0.01,
                duration_s: 1.0,
                internal_heat_j: 0.0,
            };
            let saved = parcel.clone();
            let environment = CondensingGasThermalEnvironment {
                parcel: &parcel,
                body: &boundary,
                gas_heat_j: 0.0,
                gas_solve: controls(),
                exchange: exchange(),
                maximum_gas_temperature_increment_k: 30.0,
            };
            let step = environment.advance(cx, input).unwrap();
            let report = &step.report;
            assert!(report.solid_to_gas_heat_j < 0.0);
            assert!(report.gas.phase.condensed_mass_kg > 0.0);
            assert!(step.state.temperature_k() > state.temperature_k());
            assert!(report.gas.phase.parcel.temperature_k() < parcel.temperature_k());
            near(
                report.solid_to_gas_heat_j,
                input.duration_s
                    * input.surface_area_m2
                    * boundary.convection_w_per_m2_k
                    * (step.state.temperature_k() - report.gas.phase.parcel.temperature_k()),
                1e-8,
            );
            assert!(report.combined_energy_residual_j.abs() < 1e-7);
            assert!(report.solid_energy_residual_j.abs() < 1e-7);
            near(water_mass(&report.gas.phase.parcel), original_water, 1e-16);
            near(
                carrier_moles(&report.gas.phase.parcel),
                original_carrier,
                1e-15,
            );
            assert_eq!(parcel, saved);
            total_heat += report.solid_to_gas_heat_j;
            state = step.state;
            parcel = step.report.gas.phase.parcel;
        }
        near(
            parcel.total_internal_energy_j() - initial_energy,
            total_heat,
            1e-7,
        );
        near(
            0.01 * (state.specific_enthalpy_j_kg() - initial_h) + total_heat,
            0.0,
            1e-7,
        );
        near(
            0.01 * (state.specific_enthalpy_j_kg() - initial_h) + parcel.total_internal_energy_j()
                - initial_energy,
            0.0,
            1e-7,
        );
        assert!(parcel.liquid_mass_kg() > 1e-5);
    });
}

#[test]
fn convection_uses_the_actual_depleted_vapor_and_preserves_pack_identity() {
    scope(|cx, _| {
        let curve = curve();
        let boundary = boundary();
        let parcel = wet_parcel(310.0, 0.05, 0.0);
        let environment = CondensingGasThermalEnvironment {
            parcel: &parcel,
            body: &boundary,
            gas_heat_j: 0.0,
            gas_solve: controls(),
            exchange: exchange(),
            maximum_gas_temperature_increment_k: 30.0,
        };
        let input = UniformEnthalpyStepInput {
            curve: &curve,
            initial: curve.state_at_specific_enthalpy(15000.0).unwrap(),
            mass_kg: 0.01,
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
            duration_s: 1.0,
            internal_heat_j: 0.0,
        };
        let step = environment
            .advance_with_convection(cx, input, |gas| {
                Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 50.0 + 1000.0 * gas.state().water_mole_fraction,
                    source: gas.clone(),
                })
            })
            .unwrap();
        let next = step.report.gas.phase.parcel.gas().state();
        assert_eq!(&step.report.convection.source, next);
        assert_eq!(next.source_packs(), parcel.gas().state().source_packs());
        assert!(
            next.state().water_mole_fraction < parcel.gas().state().state().water_mole_fraction
        );
        assert_eq!(
            step.report.convection.coefficient_w_per_m2_k,
            50.0 + 1000.0 * next.state().water_mole_fraction
        );
        assert!(step.report.combined_energy_residual_j.abs() < 1e-7);
    });
}

#[test]
fn inadequate_energy_brackets_cancellation_and_refused_convection_publish_nothing() {
    let parcel = wet_parcel(310.0, 0.05, 0.0);
    let saved = parcel.clone();
    let too_narrow = WetGasHeatSolve {
        temperature_bracket_k: [309.0, 311.0],
        ..controls()
    };
    assert!(matches!(
        heat_candidate(&parcel, -20.0, too_narrow, || false),
        Err(Error::Thermal(ThermalPistonError::OutsideBracket { .. }))
    ));
    assert!(matches!(
        heat_candidate(&parcel, 0.0, controls(), || true),
        Err(Error::Thermal(ThermalPistonError::Cancelled))
    ));
    for invalid in [f64::NAN, f64::INFINITY] {
        assert!(heat_candidate(&parcel, invalid, controls(), || false).is_err());
    }
    scope(|cx, gate| {
        let curve = curve();
        let boundary = boundary();
        let input = UniformEnthalpyStepInput {
            curve: &curve,
            initial: curve.state_at_specific_enthalpy(15000.0).unwrap(),
            mass_kg: 0.01,
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
            duration_s: 1.0,
            internal_heat_j: 0.0,
        };
        let environment = CondensingGasThermalEnvironment {
            parcel: &parcel,
            body: &boundary,
            gas_heat_j: 0.0,
            gas_solve: controls(),
            exchange: exchange(),
            maximum_gas_temperature_increment_k: 30.0,
        };
        assert!(matches!(
            environment.advance_with_convection(
                cx,
                input,
                |_| Err::<ConvectionCoefficient<()>, _>("physical correlation refused")
            ),
            Err(Error::Convection("physical correlation refused"))
        ));
        gate.request();
        assert!(matches!(
            environment.advance(cx, input),
            Err(Error::Thermal(ThermalPistonError::Cancelled))
        ));
    });
    assert_eq!(parcel, saved);
}
