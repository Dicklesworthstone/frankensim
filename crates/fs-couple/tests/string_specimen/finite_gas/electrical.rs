//! G1/G3/G4: ideal electrical work reaches the same vibrating solid and finite gas.
use super::*;
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, layers::MaterialConductionLayers};
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;

const RATE: u32 = 48000;
pub(super) const CP: f64 = 0.001; // Deliberately tiny manufactured capacity; no measured material claim.

pub(super) fn conducting_string(
    cx: &fs_exec::Cx<'_>,
    gas: &fs_material::gas::GasState,
    temperature: f64,
    pluck: bool,
    final_resistivity: f64,
) -> (
    ThermalMaterialStringRuntime,
    MaterialCard,
    QueryPoint,
    PropertyKey,
) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, initial, final_value, constant) in [
        ("density", Density::DIMS, 1000.0, 1000.0, true),
        ("young_modulus", Pressure::DIMS, 2e9, 2e9, true),
        (
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            1e7,
            9e7,
            false,
        ),
        (
            ELECTRICAL_RESISTIVITY_PROPERTY,
            ELECTRICAL_RESISTIVITY_DIMS,
            1e-6,
            final_resistivity,
            false,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            20.0,
            20.0,
            true,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if constant {
                    PropertyValue::Scalar {
                        value: initial,
                        dims,
                    }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(300.0, initial), (340.0, final_value)],
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 340.0)
                    .with(
                        "omega",
                        if name == "thermal_conductivity" {
                            0.0
                        } else {
                            1.0
                        },
                        20000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if constant {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite-gas electrically heated wire".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 1.0)
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let specimen = with_uniform_circular_material_state(loss_template(), 0.0015, &material)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap();
    let mechanical = LinearMaterialStringRuntime::try_new_with_gas(
        cx,
        specimen,
        pluck.then_some(Pluck {
            station_frac: 0.4,
            height_m: 1e-4,
        }),
        gas,
        1.0,
        RATE,
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: CP * 40.0,
                temperature_k: 340.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            // The shared solid/liquid chart must cover its terminal phase.
            // These tests operate at or below the fully solid endpoint above.
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: CP * 40.0 + 0.1,
                temperature_k: 340.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    let live = ThermalMaterialStringRuntime::try_new(
        cx,
        mechanical,
        card.clone(),
        curve,
        CP * (temperature - 300.0),
        100.0,
    )
    .unwrap();
    (
        live,
        card,
        point,
        PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    )
}

pub(super) fn controls() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-2e-5, 2e-5],
        source_tolerance_j: 1e-13,
        energy_tolerance_j: 2e-10,
        maximum_iterations: 64,
    }
}

#[test]
fn g1_ideal_current_and_voltage_heat_one_source_conductor_and_finite_gas_once() {
    with_string_cx(|cx, _| {
        let gas = parcel(0.028, 100000.0);
        let (initial, _, _, _) = conducting_string(cx, gas.state().state(), 300.0, false, 9e-6);
        let specimen = initial.mechanical().specimen();
        let r0 = 1e-6 * specimen.string().length_m / specimen.area_m2();
        let mass = specimen.mass_kg();
        let current: f64 = 0.05;
        let heat = FiniteGasStringHeat {
            string_heat_j: 3e-9,
            gas_heat_j: 1e-6,
        };
        let joule = current.powi(2) * r0 / f64::from(RATE);
        let expected_temperature = 300.0 + (joule + heat.string_heat_j) / (mass * CP);
        let gas_cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
        for drive in [
            OhmicDrive::Current(current),
            OhmicDrive::Current(-current),
            OhmicDrive::Voltage(current * r0),
            OhmicDrive::Voltage(-current * r0),
        ] {
            let mut live = initial.clone();
            let mut next_gas = gas.clone();
            let frame = live
                .step_with_drive_and_isochoric_gas(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut next_gas,
                    heat,
                    solve(),
                )
                .unwrap();
            close(frame.joule_heat_j, joule);
            close(live.thermal().temperature_k(), expected_temperature);
            close(frame.conductor_before.resistance_ohm(), r0);
            close(
                frame.conductor_after.resistance_ohm(),
                r0 * (1.0 + 0.2 * (expected_temperature - 300.0)),
            );
            assert_eq!(frame.external_string_heat_j, heat.string_heat_j);
            assert_eq!(
                frame.coupled.string.external_heat_j,
                heat.string_heat_j + frame.joule_heat_j
            );
            assert_eq!(frame.coupled.drag_heat_j, 0.0);
            assert_eq!(
                frame.coupled.string.vibration.dissipation.material_heat_j,
                0.0
            );
            assert!(
                (next_gas.state().state().temperature
                    - 300.0
                    - heat.gas_heat_j / (gas.mass_kg() * gas_cv))
                    .abs()
                    < 2e-10
            );
            let stored = mass * live.thermal().specific_enthalpy_j_kg()
                + gas.mass_kg()
                    * (next_gas.state().specific_internal_energy_j_kg()
                        - gas.state().specific_internal_energy_j_kg());
            assert!(
                (stored - joule - heat.string_heat_j - heat.gas_heat_j).abs()
                    <= frame.coupled.energy_tolerance_j
            );
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);
            assert_eq!(live.mechanical().ambient(), next_gas.state().state());
            assert_eq!(
                frame.conductor_after.material().query_point(),
                live.mechanical().specimen().material().query_point()
            );
            match drive {
                OhmicDrive::Current(value) => assert_eq!(frame.final_current_a, value),
                OhmicDrive::Voltage(value) => close(
                    frame.final_current_a,
                    value / frame.conductor_after.resistance_ohm(),
                ),
            }
        }
    });
}

#[test]
fn g3_powered_conduction_changes_live_pressure_and_closes_all_three_storages() {
    with_string_cx(|cx, _| {
        let mut gas = parcel(0.028, 100000.0);
        let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
        gas.add_heat(gas.mass_kg() * cv * 20.0, solve().gas, || false)
            .unwrap();
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let specimen = initial.mechanical().specimen();
        let area = 2.0
            * core::f64::consts::PI
            * specimen.radius_m()
            * (specimen.string().length_m + specimen.radius_m());
        let r0 = 3e-6 * specimen.string().length_m / specimen.area_m2();
        let layers = MaterialConductionLayers {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
            conductivity: &key,
            selection: MaterialPropertySelection::SingleClaimOnly,
            area_m2: area,
            gas_distance_m: 0.01,
            body_distance_m: specimen.radius_m() / 2.0,
        };
        let mut powered = initial.clone();
        let mut unpowered = initial.clone();
        let mut reversed = initial.clone();
        let mut powered_gas = gas.clone();
        let mut control_gas = gas.clone();
        let mut reversed_gas = gas.clone();
        let initial_energy = energy(&initial);
        let mut input = 0.0;
        let mut numerical = 0.0;
        let mut allowance = 0.0;
        let mut pressure_difference = 0.0;
        for epoch in 0..128 {
            let heat = FiniteGasStringHeat {
                string_heat_j: 1e-9,
                gas_heat_j: 1e-7,
            };
            let frame = powered
                .step_with_drive_and_isochoric_gas_conduction(
                    cx,
                    epoch,
                    &[0.0],
                    OhmicDrive::Voltage(0.5 * r0),
                    &mut powered_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                )
                .unwrap();
            let opposite = reversed
                .step_with_drive_and_isochoric_gas_conduction(
                    cx,
                    epoch,
                    &[0.0],
                    OhmicDrive::Voltage(-0.5 * r0),
                    &mut reversed_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                )
                .unwrap();
            let reference = unpowered
                .step_with_drive_and_isochoric_gas_conduction(
                    cx,
                    epoch,
                    &[0.0],
                    OhmicDrive::Voltage(0.0),
                    &mut control_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                )
                .unwrap();
            let exchange = &frame.coupled;
            let coupled = &exchange.coupled;
            assert_eq!(coupled.string, opposite.coupled.coupled.string);
            assert_eq!(powered_gas, reversed_gas);
            assert_eq!(frame.current_a, -opposite.current_a);
            assert_eq!(frame.final_current_a, -opposite.final_current_a);
            close(
                frame.joule_heat_j,
                (0.5 * r0) * frame.current_a / f64::from(RATE),
            );
            close(
                frame.conductor_after.resistance_ohm(),
                1e-6 * specimen.string().length_m / specimen.area_m2()
                    * (1.0 + 0.2 * (powered.thermal().temperature_k() - 300.0)),
            );
            close(
                exchange.conductance_w_k,
                area / (0.01 / powered_gas.state().state().thermal_conductivity
                    + layers.body_distance_m / 20.0),
            );
            assert_eq!(
                exchange.external_string_heat_j,
                heat.string_heat_j + frame.joule_heat_j
            );
            assert_eq!(frame.external_string_heat_j, heat.string_heat_j);
            assert_eq!(
                exchange.source.card_identity(),
                frame.conductor_after.material().card_identity()
            );
            assert_eq!(
                exchange.source.query_point(),
                frame.conductor_after.material().query_point()
            );
            assert!(coupled.energy_balance_residual_j.abs() <= coupled.energy_tolerance_j);
            assert!(exchange.source_residual_j.abs() <= controls().source_tolerance_j);
            input += frame.joule_heat_j
                + heat.string_heat_j
                + heat.gas_heat_j
                + coupled.string.vibration.acoustic.input_work_j;
            numerical += coupled.string.vibration.dissipation.roundoff_residual_j
                + coupled.string.vibration.dissipation.authored_loss_j;
            allowance += coupled.energy_tolerance_j;
            let stored = energy(&powered) - initial_energy
                + specimen.mass_kg()
                    * (powered.thermal().specific_enthalpy_j_kg()
                        - initial.thermal().specific_enthalpy_j_kg())
                + gas.mass_kg()
                    * (powered_gas.state().specific_internal_energy_j_kg()
                        - gas.state().specific_internal_energy_j_kg());
            assert!((stored + numerical - input).abs() <= allowance);
            pressure_difference += (coupled.string.vibration.acoustic.observer_pressure_pa
                - reference
                    .coupled
                    .coupled
                    .string
                    .vibration
                    .acoustic
                    .observer_pressure_pa)
                .abs();
            assert_eq!(powered.mechanical().epoch(), epoch + 1);
            assert_eq!(powered.mechanical().accepted_samples(), epoch + 1);
            assert_eq!(powered.mechanical().ambient(), powered_gas.state().state());
        }
        assert!(powered.thermal().temperature_k() > unpowered.thermal().temperature_k());
        assert!(
            powered.mechanical().modes()[0].damping_ratio
                > unpowered.mechanical().modes()[0].damping_ratio
        );
        assert!(pressure_difference > 1e-10);
        assert!(powered_gas.state().state().temperature > control_gas.state().state().temperature);
    });
}

#[test]
fn g4_late_electrical_refusal_preserves_both_states_and_sound_history() {
    with_string_cx(|cx, gate| {
        let gas = parcel(0.028, 100000.0);
        // The thermal/material path admits zero as a retained auxiliary value;
        // the actual conductor must refuse it at the endpoint, after gas heat.
        let (mut live, _, _, _) = conducting_string(cx, gas.state().state(), 300.0, false, 0.0);
        let original = live.clone();
        let mut next_gas = gas.clone();
        let heat = FiniteGasStringHeat {
            string_heat_j: live.mechanical().specimen().mass_kg() * CP * 40.0,
            gas_heat_j: 1e-6,
        };
        assert!(matches!(
            live.step_with_drive_and_isochoric_gas(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.0),
                &mut next_gas,
                heat,
                solve()
            ),
            Err(
                fs_couple::acoustic_realize::finite_gas::FiniteGasStringError::String(
                    ThermalStringError::Electrical(_)
                )
            )
        ));
        assert_eq!(live.thermal(), original.thermal());
        assert_eq!(live.mechanical().states(), original.mechanical().states());
        assert_eq!(
            live.mechanical().specimen(),
            original.mechanical().specimen()
        );
        assert_eq!(live.mechanical().epoch(), 0);
        assert_eq!(live.mechanical().accepted_samples(), 0);
        assert_eq!(next_gas, gas);
        let mut retry = original.clone();
        let mut retry_gas = gas.clone();
        let valid = FiniteGasStringHeat {
            string_heat_j: heat.string_heat_j / 40.0,
            gas_heat_j: heat.gas_heat_j,
        };
        let accepted = live
            .step_with_drive_and_isochoric_gas(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.0),
                &mut next_gas,
                valid,
                solve(),
            )
            .unwrap();
        let replay = retry
            .step_with_drive_and_isochoric_gas(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.0),
                &mut retry_gas,
                valid,
                solve(),
            )
            .unwrap();
        assert_eq!(accepted.coupled.string, replay.coupled.string);
        assert_eq!(next_gas, retry_gas);
        for (epoch, drive) in [
            (0, OhmicDrive::Current(0.0)),
            (1, OhmicDrive::Voltage(f64::NAN)),
        ] {
            assert!(
                live.step_with_drive_and_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0],
                    drive,
                    &mut next_gas,
                    valid,
                    solve()
                )
                .is_err()
            );
        }
        gate.request();
        assert!(
            live.step_with_drive_and_isochoric_gas(
                cx,
                1,
                &[0.0],
                OhmicDrive::Current(0.0),
                &mut next_gas,
                valid,
                solve()
            )
            .is_err()
        );
        assert_eq!(live.thermal(), retry.thermal());
        assert_eq!(live.mechanical().states(), retry.mechanical().states());
        assert_eq!(live.mechanical().epoch(), 1);
        assert_eq!(next_gas, retry_gas);
    });
}
