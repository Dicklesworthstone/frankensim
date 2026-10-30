//! Authored powers bind the existing caloric inventories and joint DG transaction.
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::thermal_wall::FiniteWallThermalSolve;
use crate::bernoulli_aperture::dg::thermal_wall::gas::FiniteGasWallThermalSolve;
use crate::bernoulli_aperture::performance::gas_wall::{
    FiniteGasWallHeatProgram, GasWallConduction, GasWallHeatControl,
};
use crate::bernoulli_aperture::performance::thermal_wall::{
    FiniteWallHeatProgram, WallHeatControl,
};
use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
use crate::dg_port::finite_gas::{FiniteGasDgLoss, FiniteGasDgSolve, admitted_volume};
use crate::render::RenderError;
use fs_matdb::{QueryPoint, SelectionPolicy};
use fs_material::gas_species::input::{
    ThermochemicalPhase, ThermochemicalSources, read_frozen_thermochemical_model,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve,
    SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_phs::linear_exact::DenseExponentialBudget;
use std::sync::Arc;

#[path = "../../../../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const DT: f64 = 1.0 / 192000.0;
const WALL_HEAT: &str = concat!(
    "finite_wall solid\n",
    "wall_caloric specific_heat_capacity density T 200 800 0 1e-5 single-claim\n",
    "wall_inventory 0.00004 0.00000001 0.0004 50000\n",
    "wall_feedback -0.0002 0.0002 1e-13 1e-11 80 0.0001\n",
    "wall_dissipation none\n",
    "wall_heat 3\nwall_power 0 0\nwall_power 7 600\nwall_power 19 -300\n",
);
const INTERFACE: &str = concat!(
    "gas_wall_conduction thermal_conductivity 0.0004 0.00002 0.00002\n",
    "gas_wall_feedback -0.001 0.001 1e-13 1e-11 80 20000\n",
);
const NATIVE: &str = "wall_dissipation native 512 67108864 40000000000";

fn thermal_input(wall: &str, native: bool) -> String {
    input(wall)
        .replace("audio 192000 64 20", "audio 192000 24 20")
        .replace(
            "observation inlet",
            &format!("{WALL_HEAT}observation inlet"),
        )
        .replace(
            "wall_dissipation none",
            if native {
                NATIVE
            } else {
                "wall_dissipation none"
            },
        )
}

fn joint_cold(temperature: f64) -> String {
    let gas = format!("thermochemical {temperature} 101325 mole 1 custom-alpha 1");
    input(WALL)
        .replace("audio 192000 64 20", "audio 192000 24 20")
        .replace("mixture 293.15 101325 mole 2 N2 0.7 O2 0.3", &gas)
        .replace(
            DG,
            &format!(
                "dg_bore 1 open 2 2 8 50 2000 0.01 1e-11 256 4 gas_field 2 mole 1 custom-alpha\n\
             dg_segment cylinder 0.001 0.05\n\
             gas_station 0 {temperature} 1\n\
             gas_station 0.05 {temperature} 1"
            ),
        )
}

fn joint_input(temperature: f64) -> String {
    let heat = format!(
        "finite_gas thermochemical {temperature} 101325 mole 1 custom-alpha 1\n\
         gas_energy 200 900 1e-8 128\n\
         gas_feedback -1e-5 1e-5 1e-15 1e-11 80 80 spatial-viscothermal 0.01\n\
         gas_heat 3\nheat 0 0\nheat 3 2\nheat 11 -1\n"
    );
    joint_cold(temperature).replace(
        "observation inlet",
        &format!("{heat}{WALL_HEAT}{INTERFACE}observation inlet"),
    )
}

fn joint_loaded(
    text: &str,
    material: &NormalizedMaterialCardPack,
    block: usize,
) -> Result<PlateValvePerformance, PlateValveInputError> {
    let gas = gas_support::packs("custom-alpha", 0.028, 3.5, 200.0);
    let bindings = [gas.binding()];
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    PlateValvePerformance::from_bytes_with_all_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        Some(&registry),
        |_| Ok(material.clone()),
    )
}

fn source_point(cold: &PlateValvePerformance) -> QueryPoint {
    cold.dg_wall_layout()
        .unwrap()
        .uniform_record()
        .unwrap()
        .wall()
        .material()
        .query()
        .clone()
}

fn caloric(material: &NormalizedMaterialCardPack, point: &QueryPoint) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        material.card(),
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point.clone().with("T", 200.0).unwrap(),
            upper: &point.clone().with("T", 800.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap()
}

fn body(curve: &EquilibriumEnthalpyPhaseCurve) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(50000.0).unwrap(),
        mass_kg: 0.00004,
        volume_m3: 0.00000001,
        surface_area_m2: 0.0004,
        internal_heat_j: 0.0,
        duration_s: DT,
    }
}

fn wall_solve(native: bool) -> FiniteWallThermalSolve {
    FiniteWallThermalSolve {
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-0.0002, 0.0002],
            source_tolerance_j: 1e-13,
            energy_tolerance_j: 1e-11,
            maximum_iterations: 80,
        },
        maximum_absolute_parameter_work_j: 0.0001,
        dissipation_heat: native.then_some(DenseExponentialBudget {
            max_dimension: 512,
            max_dense_bytes: 67108864,
            max_multiply_terms: 40000000000,
        }),
    }
}

fn wall_program(native: bool) -> FiniteWallHeatProgram {
    FiniteWallHeatProgram {
        controls: [(0, 0.0), (7, 600.0), (19, -300.0)]
            .into_iter()
            .map(|(sample, power_w)| WallHeatControl { sample, power_w })
            .collect(),
        coupling: wall_solve(native),
    }
}

fn joint_program() -> FiniteGasWallHeatProgram {
    FiniteGasWallHeatProgram {
        controls: [
            (0, 0.0, 0.0),
            (3, 2.0, 0.0),
            (7, 2.0, 600.0),
            (11, -1.0, 600.0),
            (19, -1.0, -300.0),
        ]
        .into_iter()
        .map(|(sample, gas_power_w, wall_power_w)| GasWallHeatControl {
            sample,
            gas_power_w,
            wall_power_w,
        })
        .collect(),
        coupling: FiniteGasWallThermalSolve {
            gas: FiniteGasDgSolve {
                work_bracket_j: [-1e-5, 1e-5],
                work_tolerance_j: 1e-15,
                energy_tolerance_j: 1e-11,
                maximum_iterations: 80,
                maximum_temperature_increment_k: 80.0,
                loss: FiniteGasDgLoss::SpatialViscothermal {
                    maximum_relative_error: 0.01,
                },
            },
            caloric: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            wall: wall_solve(false),
            exchange: HeatExchangeSolve {
                heat_bracket_j: [-0.001, 0.001],
                source_tolerance_j: 1e-13,
                energy_tolerance_j: 1e-11,
                maximum_iterations: 80,
            },
            maximum_material_trials: 20000,
        },
    }
}

fn same_physics(a: &AperturePerformance, b: &AperturePerformance) {
    assert_eq!(a.samples_rendered(), b.samples_rendered());
    assert_eq!(a.system().aperture().state(), b.system().aperture().state());
    assert_eq!(a.system().aperture().spec(), b.system().aperture().spec());
    assert_eq!(dg(a).duct_state(), dg(b).duct_state());
    assert_eq!(
        dg(a).load().material_wall_source(),
        dg(b).load().material_wall_source()
    );
    assert_eq!(
        dg(a).load().duct().gas_field(),
        dg(b).load().duct().gas_field()
    );
}

#[test]
fn authored_wall_heat_matches_manual_attachment_with_original_caloric_receipts() {
    let material = wall_pack(1e8, None);
    for native in [false, true] {
        let cold_text = input(WALL).replace("audio 192000 64 20", "audio 192000 24 20");
        let cold = loaded(&cold_text, &material, 1).unwrap();
        let point = source_point(&cold);
        let caloric = caloric(&material, &point);
        let program = wall_program(native);
        let mut direct = cold
            .into_renderer()
            .with_finite_wall_heat(
                body(caloric.curve()),
                Arc::new(material.card().clone()),
                "T".into(),
                program.clone(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap();
        let mut actual = loaded(&thermal_input(WALL, native), &material, 1)
            .unwrap()
            .into_renderer();
        let initial = actual.finite_wall_body().unwrap();
        assert_eq!(initial.curve, caloric.curve());
        assert_eq!(initial.curve.source_receipts(), caloric.receipts());
        assert!(!caloric.receipts().is_empty());
        assert_eq!(initial.initial.temperature_k(), 300.0);
        assert_eq!(initial.mass_kg, 0.00004);
        assert_eq!(actual.finite_wall_material(), Some(material.card()));
        assert_eq!(actual.finite_wall_heat_program(), Some(&program));
        let mut native_heat = 0.0;
        for sample in 0..24 {
            let mut a = [f64::NAN];
            let mut b = [f64::NAN];
            actual.block(&mut a).unwrap();
            direct.block(&mut b).unwrap();
            assert_eq!(a, b);
            same_physics(&actual, &direct);
            assert_eq!(actual.finite_wall_state(), direct.finite_wall_state());
            assert_eq!(
                actual.last_finite_wall_step(),
                direct.last_finite_wall_step()
            );
            let power = program
                .controls
                .iter()
                .rev()
                .find(|event| event.sample <= sample)
                .unwrap();
            let report = actual.last_finite_wall_step().unwrap();
            assert_eq!(report.external_heat_j, power.power_w * DT);
            native_heat += report.wall_heat_j;
            assert!(report.body_energy_residual_j.abs() <= 1e-11);
            let source = dg(&actual)
                .load()
                .material_wall_source()
                .unwrap()
                .uniform_wall()
                .unwrap();
            assert_eq!(
                source.material().query().axes()["T"],
                report.body_state.temperature_k()
            );
            assert_eq!(source.material().query().axes()["omega"], 0.0);
            assert_eq!(
                source.material().card_identity(),
                material.card().content_hash()
            );
        }
        assert_eq!(native_heat > 0.0, native);
        assert!(actual.finite_wall_state().unwrap().temperature_k() > 300.0);
    }
    let regional = loaded(&thermal_input(REGIONAL, false), &material, 1).unwrap();
    assert_eq!(
        regional
            .renderer()
            .finite_wall_body()
            .unwrap()
            .curve
            .material_card_identity(),
        material.card().content_hash()
    );
}

#[test]
fn joint_authored_powers_merge_on_the_source_clock_and_use_both_actual_temperatures() {
    let material = wall_pack(1e8, None);
    for temperature in [280.0, 360.0] {
        let cold = joint_loaded(&joint_cold(temperature), &material, 1).unwrap();
        let point = source_point(&cold);
        let caloric = caloric(&material, &point);
        let packs = gas_support::packs("custom-alpha", 0.028, 3.5, 200.0);
        let bindings = [packs.binding()];
        let registry = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let declaration = format!("thermochemical {temperature} 101325 mole 1 custom-alpha 1");
        let (model, gas_point) =
            read_frozen_thermochemical_model(&mut declaration.split_ascii_whitespace(), &registry)
                .unwrap();
        let probe = ClosedGasParcel::new(model.clone(), gas_point.clone(), 1.0).unwrap();
        let volume =
            admitted_volume::<GasEnergyError>(dg(cold.renderer()).load(), &probe, &mut || false)
                .unwrap();
        let gas = ClosedGasParcel::new(model, gas_point, volume).unwrap();
        let conduction = GasWallConduction {
            context: point,
            conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            selection: MaterialPropertySelection::SingleClaimOnly,
            liquid_fraction_axis: None,
            area_m2: 0.0004,
            gas_distance_m: 0.00002,
            body_distance_m: 0.00002,
        };
        let program = joint_program();
        let mut direct = cold
            .into_renderer()
            .with_finite_gas_wall_heat(
                gas,
                body(caloric.curve()),
                Arc::new(material.card().clone()),
                "T".into(),
                conduction.clone(),
                program.clone(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap();
        let mut actual = joint_loaded(&joint_input(temperature), &material, 1)
            .unwrap()
            .into_renderer();
        assert_eq!(actual.finite_gas_wall_heat_program(), Some(&program));
        assert_eq!(actual.finite_gas_wall_conduction(), Some(&conduction));
        assert_eq!(
            actual.finite_gas_wall_body().unwrap().curve,
            caloric.curve()
        );
        assert!(actual.finite_wall_state().is_none());
        assert!(actual.finite_gas().is_none());
        let mut pressure = Vec::new();
        for sample in 0..24 {
            let mut a = [f64::NAN];
            let mut b = [f64::NAN];
            actual.block(&mut a).unwrap();
            direct.block(&mut b).unwrap();
            assert_eq!(a, b);
            same_physics(&actual, &direct);
            assert_eq!(actual.finite_gas_wall_gas(), direct.finite_gas_wall_gas());
            assert_eq!(
                actual.finite_gas_wall_state(),
                direct.finite_gas_wall_state()
            );
            let report = actual.last_finite_gas_wall_step().unwrap();
            let reference = direct.last_finite_gas_wall_step().unwrap();
            pressure.push(a[0]);
            assert_eq!(report.coupled, reference.coupled);
            assert_eq!(report.exchanged_heat_j, reference.exchanged_heat_j);
            assert_eq!(report.interface, reference.interface);
            let held = program
                .controls
                .iter()
                .rev()
                .find(|event| event.sample <= sample)
                .unwrap();
            assert_eq!(report.external_gas_heat_j, held.gas_power_w * DT);
            assert_eq!(report.external_body_heat_j, held.wall_power_w * DT);
            if sample == 0 {
                assert_eq!(
                    report.exchanged_heat_j.is_sign_positive(),
                    temperature < 300.0
                );
                assert!(report.exchanged_heat_j.abs() > 1e-6);
            }
            assert_eq!(
                report.interface.source.card_identity(),
                material.card().content_hash()
            );
            assert_eq!(
                report.interface.source.query().axes()["T"],
                report.body_midpoint.temperature_k()
            );
            assert_eq!(report.interface.source.query().axes()["omega"], 0.0);
            let conductivity = report
                .interface
                .source
                .property("thermal_conductivity")
                .unwrap()
                .value_si();
            assert!(
                (conductivity - (3.0 + 0.01 * report.body_midpoint.temperature_k())).abs() < 1e-13
            );
            assert!(report.combined_energy_residual_j.abs() <= 1e-11);
        }
        assert_eq!(
            actual.applied_gas_wall_heat_controls(),
            program.controls.as_slice()
        );
        let (chunked, final_state) = rendered(
            joint_loaded(&joint_input(temperature), &material, 7).unwrap(),
            7,
        );
        assert_eq!(chunked, pressure);
        same_physics(&actual, &final_state);
        assert_eq!(
            actual.finite_gas_wall_gas(),
            final_state.finite_gas_wall_gas()
        );
        assert_eq!(
            actual.finite_gas_wall_state(),
            final_state.finite_gas_wall_state()
        );
    }
}

#[test]
fn finite_thermal_sources_bindings_controls_and_shared_budgets_refuse_without_fallback() {
    let material = wall_pack(1e8, None);
    let wall = thermal_input(WALL, false);
    for text in [
        wall.replace("finite_wall solid", "finite_wall liquid"),
        wall.replace("density T 200", "density omega 200"),
        wall.replace("wall_inventory 0.00004", "wall_inventory -0.00004"),
        wall.replace("wall_power 7 600", "wall_power 0 600"),
        wall.replace("wall_power 19 -300", "wall_power 24 -300"),
        wall.replace("wall_power 0 0", "wall_power 1 0"),
        wall.replace("wall_power 7 600", "wall_power 7 5e-324"),
        wall.replace("1e-11 80 0.0001", "1e-11 0 0.0001"),
        wall.replace(WALL, "dg_wall rigid\n"),
        wall.replace(WALL, ""),
        wall.replace(
            "observation inlet",
            &format!("{INTERFACE}observation inlet"),
        ),
        thermal_input(REGIONAL, false)
            .replace("wall.matpack 300 0 0.00004", "wall.matpack 300 1 0.00004"),
    ] {
        assert!(
            loaded(&text, &material, 1).is_err(),
            "accepted invalid thermal declaration"
        );
    }
    let error = loaded(&wall, &wall_pack(1e8, Some("specific_heat_capacity")), 1)
        .err()
        .unwrap();
    assert!(std::error::Error::source(&error).is_some());
    let joint = joint_input(280.0);
    for text in [
        joint.replace(INTERFACE, ""),
        joint.replace("heat 3 2", "heat 0 2"),
        joint.replace("heat 11 -1", "heat 24 -1"),
        joint.replace("heat 0 0", "heat 1 0"),
        joint.replace("compile_limits 10000 128", "compile_limits 10000 10"),
        joint.replace("1e-11 80 20000", "1e-11 80 0"),
    ] {
        assert!(joint_loaded(&text, &material, 1).is_err());
    }
    assert!(
        joint_loaded(
            &joint.replace("compile_limits 10000 128", "compile_limits 10000 11"),
            &material,
            1,
        )
        .is_ok()
    );
    assert!(joint_loaded(&joint, &wall_pack(1e8, Some("thermal_conductivity")), 1).is_err());
    let different = thermal_input(REGIONAL, false).replace(
        "0.025 0.05 isotropic wall.matpack",
        "0.025 0.05 isotropic other.matpack",
    );
    assert!(
        PlateValvePerformance::from_bytes_with_all_sources(
            different.as_bytes(),
            1,
            &CancelGate::new_clock_free(),
            None,
            |path| Ok(if path == "wall.matpack" {
                material.clone()
            } else {
                wall_pack(4e8, None)
            }),
        )
        .is_err()
    );
}

#[test]
fn cancelling_authored_wall_or_joint_heat_preserves_all_accepted_histories_and_cursors() {
    let material = wall_pack(1e8, None);
    for joint in [false, true] {
        let mut player = if joint {
            joint_loaded(&joint_input(280.0), &material, 7).unwrap()
        } else {
            loaded(&thermal_input(WALL, false), &material, 7).unwrap()
        }
        .into_renderer();
        player.block(&mut [0.0; 7]).unwrap();
        let before = dg(&player).clone();
        let wall = player.finite_wall_state().copied();
        let joint_wall = player.finite_gas_wall_state().copied();
        let gas = player.finite_gas_wall_gas().cloned();
        let powers = player.applied_wall_heat_controls().to_vec();
        let paired = player.applied_gas_wall_heat_controls().to_vec();
        let pressure = player.applied_controls().to_vec();
        if joint {
            player.finite_gas_wall_cancel_gate().unwrap().request();
        } else {
            player.finite_wall_cancel_gate().unwrap().request();
        }
        let mut output = [123.0];
        assert!(matches!(
            player.block(&mut output),
            Err(RenderError::Voice(AcousticRealizeError::Cancelled))
        ));
        assert_eq!(output, [123.0]);
        assert_eq!(player.samples_rendered(), 7);
        assert_eq!(dg(&player).duct_state(), before.duct_state());
        assert_eq!(dg(&player).aperture().state(), before.aperture().state());
        assert_eq!(
            dg(&player).load().material_wall_source(),
            before.load().material_wall_source()
        );
        assert_eq!(
            dg(&player).load().duct().gas_field(),
            before.load().duct().gas_field()
        );
        assert_eq!(player.finite_wall_state().copied(), wall);
        assert_eq!(player.finite_gas_wall_state().copied(), joint_wall);
        assert_eq!(player.finite_gas_wall_gas(), gas.as_ref());
        assert_eq!(player.applied_wall_heat_controls(), powers.as_slice());
        assert_eq!(player.applied_gas_wall_heat_controls(), paired.as_slice());
        assert_eq!(player.applied_controls(), pressure.as_slice());
        assert!(matches!(
            player.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    }
}
