//! Heat-transfer physics over the original binary phase-equilibrium inventory.
#![cfg(all(feature = "material", feature = "experimental-euler"))]

// Share only the authored source fixture; all reference outcomes below come
// from independent cubic/fugacity calorics, not a copied production EOS.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, model, source};
use fs_gas::material::phase_thermal::{
    IsobaricMixtureThermalControl, MixtureThermalControl, exchange_isobaric_mixture_parcel,
    exchange_mixture_parcel,
};
use fs_gas::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallDrag, WallExchangeError, WallThermal,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    equilibrium::{PhaseEnthalpyControls, volume::PhaseEnergyVolumeControls},
    parcel::{IsochoricMixtureEquilibriumParcel, isobaric::IsobaricMixtureEquilibriumParcel},
};
use fs_qty::{Pressure, Temperature, Time};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
fn equilibrium() -> PhaseEnergyVolumeControls {
    PhaseEnergyVolumeControls::new(
        [Pressure::new(2e6), Pressure::new(8e6)],
        PhaseEnthalpyControls::new([Temperature::new(175.0), Temperature::new(325.0)]),
    )
}
fn parcel(model: ThermochemicalPengRobinsonModel, mass: f64) -> IsochoricMixtureEquilibriumParcel {
    IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
        model,
        mass,
        Temperature::new(220.0),
        Pressure::new(3e6),
        equilibrium().phase.phase_options(),
        || false,
    )
    .unwrap()
}
fn controls(bracket: [f64; 2]) -> MixtureThermalControl {
    MixtureThermalControl::new(bracket, equilibrium())
}

fn isobaric_controls(bracket: [f64; 2]) -> IsobaricMixtureThermalControl {
    IsobaricMixtureThermalControl::new(
        bracket,
        PhaseEnthalpyControls::new([Temperature::new(175.0), Temperature::new(285.0)]),
    )
}

fn isobaric_parcel(
    model: ThermochemicalPengRobinsonModel,
    mass: f64,
    temperature: f64,
) -> IsobaricMixtureEquilibriumParcel {
    IsobaricMixtureEquilibriumParcel::from_pressure_temperature(
        model,
        mass,
        Temperature::new(temperature),
        Pressure::new(3e6),
        isobaric_controls([-1.0, 1.0]).equilibrium.phase_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_isobaric_bath_heat_crosses_bubble_point_with_actual_expansion_work() {
    let mut gas = isobaric_parcel(source(), 2.0 * MIXTURE_MASS, 180.0);
    let before = gas.clone();
    // Independent 50-digit cubic/fugacity references retained in the original
    // PH battery: at 3 MPa, 180 -> 220 K crosses the 205.7464 K bubble point.
    let heat = gas.mass_kg() * (-16_172.250_165_027_297 + 203_472.419_617_408_38);
    let stored = gas.mass_kg() * (-39_099.945_545_128_86 + 209_511.548_518_424_48);
    let mut bath = WallThermal::Reservoir {
        temperature_k: 300.0,
        conductance_w_per_k: 1.0,
    };
    let report = exchange_isobaric_mixture_parcel(
        &mut gas,
        &mut bath,
        None,
        Time::new(heat / (300.0 - 220.0)),
        isobaric_controls([0.0, 1.05 * heat]),
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 220.0, 1e-6);
    assert!(gas.state().phase().two_phase().is_some());
    close(
        gas.volume_m3(),
        gas.mass_kg() / 130.846_120_827_461_9,
        3e-12,
    );
    close(report.heat.heat_to_gas_j, heat, 2e-5);
    close(report.equilibrium.enthalpy_change_j, heat, 2e-5);
    close(report.equilibrium.stored_energy_change_j, stored, 2e-5);
    let work = 3e6 * (gas.volume_m3() - before.volume_m3());
    assert!(work > 0.0);
    assert_eq!(report.equilibrium.boundary_work_j, work);
    assert_eq!(report.heat.mechanical_work_to_gas_j, -work);
    close(stored + work, heat, 2e-5);
    close(report.heat.energy_residual_j, 0.0, 2e-5);
    close(report.heat.thermal_residual_j, 0.0, 2e-5);
    assert_eq!(
        report.heat.external_heat_to_gas_j,
        report.heat.heat_to_gas_j
    );
    assert!(report.equilibrium.entropy_change_j_k > heat / 300.0);
    assert_eq!(gas.mass_kg(), before.mass_kg());
    assert_eq!(gas.model(), before.model());
    assert_eq!(gas.pressure_pa(), 3e6);
    assert_eq!(report.pressure_evaluations, 0);
    assert!(report.pt_evaluations > report.equilibrium.pt_evaluations);
}

#[test]
fn g1_isobaric_finite_wall_and_latent_wall_close_heat_and_pressure_work() {
    let mut gas = isobaric_parcel(source(), 0.25 * MIXTURE_MASS, 220.0);
    let initial = gas.clone();
    // Independent 220/280 K PH references cross the 250.6892 K dew point.
    let heat = gas.mass_kg() * (346_713.375_192_810_4 + 16_172.250_165_027_297);
    let stored = gas.mass_kg() * (264_458.439_200_782_34 + 39_099.945_545_128_86);
    let wall_capacity = 100.0;
    let expected_wall = 400.0 - heat / wall_capacity;
    let mut wall = WallThermal::Lumped {
        temperature_k: 400.0,
        heat_capacity_j_per_k: wall_capacity,
        conductance_w_per_k: 1.0,
    };
    let heating = exchange_isobaric_mixture_parcel(
        &mut gas,
        &mut wall,
        None,
        Time::new(heat / (expected_wall - 280.0)),
        isobaric_controls([0.0, 1.001 * heat]),
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 280.0, 2e-6);
    assert!(gas.state().phase().two_phase().is_none());
    close(
        gas.volume_m3(),
        gas.mass_kg() / 36.471_975_375_322_785,
        3e-12,
    );
    close(heating.endpoint_temperature_k[1], expected_wall, 2e-7);
    close(heating.heat.wall_energy_change_j, -heat, 2e-5);
    close(heating.equilibrium.stored_energy_change_j, stored, 2e-5);
    close(
        stored + heating.equilibrium.boundary_work_j - heat,
        0.0,
        2e-5,
    );
    assert_eq!(heating.heat.external_heat_to_gas_j, 0.0);
    close(heating.heat.energy_residual_j, 0.0, 2e-5);

    // Reverse the mixture phase change into an original enthalpy wall. Its
    // temperature stays at 215 K while its latent inventory increases by Q.
    let source = latent_wall([0.0, 20_000.0]);
    let mut wall = source.initial_contact(1.0, WallDrag::None).unwrap().thermal;
    let cooling = exchange_isobaric_mixture_parcel(
        &mut gas,
        &mut wall,
        Some(&source),
        Time::new(heat / (220.0 - 215.0)),
        isobaric_controls([-1.05 * heat, 0.0]),
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 220.0, 2e-6);
    close(gas.volume_m3(), initial.volume_m3(), 3e-12);
    assert!(gas.state().phase().two_phase().is_some());
    assert_eq!(cooling.endpoint_temperature_k[1], 215.0);
    assert!(cooling.heat.mechanical_work_to_gas_j > 0.0);
    let WallThermal::Caloric {
        specific_energy_j_kg,
        ..
    } = wall
    else {
        panic!("latent wall lost its original energy coordinate")
    };
    close(specific_energy_j_kg, 5000.0 + heat, 2e-5);
    close(cooling.heat.heat_to_gas_j, -heat, 2e-5);
    close(cooling.heat.wall_energy_change_j, heat, 2e-5);
    close(cooling.heat.energy_residual_j, 0.0, 2e-5);
    close(
        cooling.equilibrium.boundary_work_j + heating.equilibrium.boundary_work_j,
        0.0,
        2e-5,
    );
}

#[test]
fn g3_isobaric_heat_transfer_preserves_large_species_reference_invariance() {
    let mass = 2.0 * MIXTURE_MASS;
    let mut plain = isobaric_parcel(source(), mass, 220.0);
    let mut shifted = isobaric_parcel(
        model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]),
        mass,
        220.0,
    );
    let wall = WallThermal::Lumped {
        temperature_k: 230.0,
        heat_capacity_j_per_k: 100.0,
        conductance_w_per_k: 1.0,
    };
    let mut first_wall = wall;
    let mut shifted_wall = wall;
    let a = exchange_isobaric_mixture_parcel(
        &mut plain,
        &mut first_wall,
        None,
        Time::new(1e-3),
        isobaric_controls([0.0, 0.02]),
        || false,
    )
    .unwrap();
    let b = exchange_isobaric_mixture_parcel(
        &mut shifted,
        &mut shifted_wall,
        None,
        Time::new(1e-3),
        isobaric_controls([0.0, 0.02]),
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(first_wall, shifted_wall);
    assert_eq!(plain.volume_m3(), shifted.volume_m3());
    assert_eq!(
        plain.state().temperature_k(),
        shifted.state().temperature_k()
    );
    assert!(a.heat.heat_to_gas_j > 0.0);
    assert!(a.equilibrium.boundary_work_j > 0.0);
}

#[test]
fn g4_isobaric_heat_refusal_and_final_cancellation_publish_neither_inventory() {
    let initial = isobaric_parcel(source(), 2.0 * MIXTURE_MASS, 220.0);
    let original_wall = WallThermal::Lumped {
        temperature_k: 300.0,
        heat_capacity_j_per_k: 100.0,
        conductance_w_per_k: 1.0,
    };
    let normal = isobaric_controls([0.0, 1000.0]);
    let mut gas = initial.clone();
    let mut wall = original_wall;
    let mut limited = normal;
    limited.maximum_pt_evaluations = 1;
    let rest = exchange_isobaric_mixture_parcel(
        &mut gas,
        &mut wall,
        None,
        Time::new(0.0),
        limited,
        || false,
    )
    .unwrap();
    assert_eq!(rest.pt_evaluations, 0);
    assert_eq!(gas, initial);
    assert_eq!(wall, original_wall);
    assert!(matches!(
        exchange_isobaric_mixture_parcel(
            &mut gas,
            &mut wall,
            None,
            Time::new(10.0),
            limited,
            || false
        ),
        Err(WallExchangeError::MixtureThermalWork { .. })
    ));
    assert_eq!(gas, initial);
    assert_eq!(wall, original_wall);
    let mut outside = normal;
    outside.equilibrium.enthalpy.temperature_bracket =
        [Temperature::new(175.0), Temperature::new(220.01)];
    assert!(
        exchange_isobaric_mixture_parcel(
            &mut gas,
            &mut wall,
            None,
            Time::new(10.0),
            outside,
            || false
        )
        .is_err()
    );
    assert_eq!(gas, initial);
    assert_eq!(wall, original_wall);
    let mut polls = 0;
    let expected = exchange_isobaric_mixture_parcel(
        &mut gas,
        &mut wall,
        None,
        Time::new(10.0),
        normal,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut retried = initial.clone();
    let mut retry_wall = original_wall;
    let mut seen = 0;
    assert_eq!(
        exchange_isobaric_mixture_parcel(
            &mut retried,
            &mut retry_wall,
            None,
            Time::new(10.0),
            normal,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(retried, initial);
    assert_eq!(retry_wall, original_wall);
    assert_eq!(
        exchange_isobaric_mixture_parcel(
            &mut retried,
            &mut retry_wall,
            None,
            Time::new(10.0),
            normal,
            || false
        )
        .unwrap(),
        expected
    );
    assert_eq!(retried, gas);
    assert_eq!(retry_wall, wall);
}

#[test]
fn g1_held_bath_matches_independent_isochoric_binary_endpoint() {
    // Independent simultaneous PR cubic/fugacity roots and species energy
    // integrals at fixed Vbar = 0.00017620507091992501 m3/mol:
    // 220 K, 3 MPa -> 230 K, 3512821.2981070513 Pa; du=37282.2427837216 J/kg.
    let mut gas = parcel(source(), 2.0 * MIXTURE_MASS);
    let before = gas.clone();
    let expected_heat = gas.mass_kg() * 37_282.242_783_721_6;
    let dt = expected_heat / (300.0 - 230.0);
    let mut bath = WallThermal::Reservoir {
        temperature_k: 300.0,
        conductance_w_per_k: 1.0,
    };
    let old_bath = bath;
    let report = exchange_mixture_parcel(
        &mut gas,
        &mut bath,
        None,
        Time::new(dt),
        controls([0.0, 2000.0]),
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 230.0, 2e-6);
    close(gas.state().pressure_pa(), 3_512_821.298_107_051_3, 0.2);
    close(report.heat.heat_to_gas_j, expected_heat, 2e-5);
    assert_eq!(
        report.heat.external_heat_to_gas_j,
        report.heat.heat_to_gas_j
    );
    assert_eq!(report.heat.wall_energy_change_j, 0.0);
    assert_eq!(report.heat.energy_residual_j, 0.0);
    assert!(report.heat.thermal_residual_j.abs() < 1.3e-5);
    assert!(report.equilibrium.entropy_change_j_k > report.heat.heat_to_gas_j / 300.0);
    assert_eq!(gas.mass_kg(), before.mass_kg());
    assert_eq!(gas.volume_m3(), before.volume_m3());
    assert_eq!(gas.model(), before.model());
    assert_eq!(bath, old_bath);
    assert!(report.pt_evaluations > report.equilibrium.pt_evaluations);
    assert!(report.pt_evaluations <= controls([0.0, 2000.0]).maximum_pt_evaluations);
    assert!(
        report.heat_evaluations < 24,
        "bounded secant should avoid an unnecessary bisection ladder"
    );
}

#[test]
fn g1_finite_wall_supplies_phase_crossing_heat_without_extra_latent_energy() {
    // Same independent isochore at 280 K: p=6938572.708280131 Pa and
    // du=218043.6284587497 J/kg; this endpoint is homogeneous.
    let mut gas = parcel(source(), 2.0 * MIXTURE_MASS);
    let before = gas.clone();
    let q = gas.mass_kg() * 218_043.628_458_749_7;
    let wall_capacity = 1000.0;
    let expected_wall = 400.0 - q / wall_capacity;
    let mut wall = WallThermal::Lumped {
        temperature_k: 400.0,
        heat_capacity_j_per_k: wall_capacity,
        conductance_w_per_k: 1.0,
    };
    let dt = q / (expected_wall - 280.0);
    assert!(gas.state().phase().two_phase().is_some());
    // This longer phase-crossing heat bracket needs more original incipient
    // phase work; stationarity/fugacity and all physical tolerances stay fixed.
    let mut control = controls([0.0, 10_400.0]);
    control
        .equilibrium
        .phase
        .enthalpy
        .flash_options
        .maximum_iterations = 256;
    let report = exchange_mixture_parcel(
        &mut gas,
        &mut wall,
        None,
        Time::new(dt),
        // Keep the upper heat endpoint inside the declared 8 MPa U,V
        // pressure envelope; 11 kJ would require more than 8 MPa.
        control,
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 280.0, 2e-6);
    close(gas.state().pressure_pa(), 6_938_572.708_280_131, 0.2);
    assert!(gas.state().phase().two_phase().is_none());
    close(report.endpoint_temperature_k[1], expected_wall, 1e-7);
    close(report.heat.heat_to_gas_j, q, 3e-5);
    close(report.heat.wall_energy_change_j, -q, 3e-5);
    assert!(report.heat.energy_residual_j.abs() < 2e-5);
    assert_eq!(report.heat.external_heat_to_gas_j, 0.0);
    close(
        gas.internal_energy_j() - before.internal_energy_j() + report.heat.wall_energy_change_j,
        0.0,
        3e-5,
    );
    assert_eq!(gas.mass_kg(), before.mass_kg());
    assert_eq!(gas.volume_m3(), before.volume_m3());
}

#[derive(Clone)]
struct LatentWall;
impl ThermalLaw for LatentWall {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=20_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 1000.0 {
            205.0 + h / 100.0
        } else if h <= 10_000.0 {
            215.0
        } else {
            215.0 + (h - 10_000.0) / 100.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.01
    }
}
fn latent_wall(range: [f64; 2]) -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(
        LatentWall,
        DepositedCell {
            id: CellId(37),
            segment_id: 2,
            body_id: 4,
            material_id: MaterialId(9),
            mass_kg: 1.0,
            specific_h_j_kg: 5000.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.001,
                surface_area_m2: 0.01,
            },
            boundary: ThermalBoundary {
                ambient_k: 215.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        range,
    )
    .unwrap()
}

#[test]
fn g1_cooling_reaches_original_latent_wall_inventory_at_constant_temperature() {
    let mut gas = parcel(source(), 2.0 * MIXTURE_MASS);
    let q = gas.mass_kg() * 37_282.242_783_721_6;
    gas.add_heat(q, equilibrium(), || false).unwrap();
    let wall_source = latent_wall([0.0, 20_000.0]);
    let mut wall = wall_source
        .initial_contact(1.0, WallDrag::None)
        .unwrap()
        .thermal;
    let report = exchange_mixture_parcel(
        &mut gas,
        &mut wall,
        Some(&wall_source),
        Time::new(q / 5.0),
        controls([-1900.0, 0.0]),
        || false,
    )
    .unwrap();
    close(gas.state().temperature_k(), 220.0, 2e-6);
    close(gas.state().pressure_pa(), 3e6, 0.2);
    assert_eq!(report.endpoint_temperature_k[1], 215.0);
    let WallThermal::Caloric {
        specific_energy_j_kg,
        mass_kg,
        ..
    } = wall
    else {
        panic!("original energy coordinate lost")
    };
    assert_eq!(mass_kg, 1.0);
    close(specific_energy_j_kg, 5000.0 + q, 3e-5);
    close(report.heat.heat_to_gas_j, -q, 3e-5);
    close(report.heat.wall_energy_change_j, q, 3e-5);
    assert!(report.heat.energy_residual_j.abs() < 2e-5);
}

#[test]
fn g3_species_energy_datums_and_extensive_rescaling_preserve_thermal_dynamics() {
    let mass = 2.0 * MIXTURE_MASS;
    let mut plain = parcel(source(), mass);
    let mut shifted = parcel(model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]), mass);
    let mut scaled = parcel(source(), 3.0 * mass);
    let original = WallThermal::Lumped {
        temperature_k: 300.0,
        heat_capacity_j_per_k: 100.0,
        conductance_w_per_k: 1.0,
    };
    let mut a = original;
    let mut b = original;
    let mut c = WallThermal::Lumped {
        temperature_k: 300.0,
        heat_capacity_j_per_k: 300.0,
        conductance_w_per_k: 3.0,
    };
    let first = exchange_mixture_parcel(
        &mut plain,
        &mut a,
        None,
        Time::new(10.0),
        controls([0.0, 1000.0]),
        || false,
    )
    .unwrap();
    let second = exchange_mixture_parcel(
        &mut shifted,
        &mut b,
        None,
        Time::new(10.0),
        controls([0.0, 1000.0]),
        || false,
    )
    .unwrap();
    assert_eq!(first, second);
    assert_eq!(a, b);
    assert_eq!(
        plain.state().temperature_k(),
        shifted.state().temperature_k()
    );
    let mut scale_controls = controls([0.0, 3000.0]);
    scale_controls.thermal.absolute_energy_tolerance_j *= 3.0;
    let third = exchange_mixture_parcel(
        &mut scaled,
        &mut c,
        None,
        Time::new(10.0),
        scale_controls,
        || false,
    )
    .unwrap();
    close(
        scaled.state().temperature_k(),
        plain.state().temperature_k(),
        2e-6,
    );
    close(
        third.endpoint_temperature_k[1],
        first.endpoint_temperature_k[1],
        2e-6,
    );
    close(
        third.heat.heat_to_gas_j,
        3.0 * first.heat.heat_to_gas_j,
        3e-5,
    );
}

#[test]
fn g0_exact_rest_and_invalid_sources_preserve_both_inventories() {
    let mut gas = parcel(source(), 2.0 * MIXTURE_MASS);
    let before = gas.clone();
    let mut capped = controls([-100.0, 100.0]);
    capped.maximum_heat_evaluations = 1;
    capped.maximum_pt_evaluations = 1;
    for (dt, tw, g) in [(0.0, 300.0, 1.0), (10.0, 220.0, 1.0), (10.0, 300.0, 0.0)] {
        let mut wall = WallThermal::Lumped {
            temperature_k: tw,
            heat_capacity_j_per_k: 100.0,
            conductance_w_per_k: g,
        };
        let old_wall = wall;
        let report =
            exchange_mixture_parcel(&mut gas, &mut wall, None, Time::new(dt), capped, || false)
                .unwrap();
        assert_eq!(gas, before);
        assert_eq!(wall, old_wall);
        assert_eq!(report.pt_evaluations, 0);
        assert_eq!(report.heat.heat_to_gas_j, 0.0);
        assert_eq!(report.heat.thermal_residual_j, 0.0);
    }
    let source = latent_wall([4999.0, 5001.0]);
    let mut wall = source.initial_contact(1.0, WallDrag::None).unwrap().thermal;
    let old_wall = wall;
    assert!(
        exchange_mixture_parcel(
            &mut gas,
            &mut wall,
            Some(&source),
            Time::new(100.0),
            controls([-100.0, 100.0]),
            || false
        )
        .is_err()
    );
    assert_eq!(gas, before);
    assert_eq!(wall, old_wall);
    assert!(
        exchange_mixture_parcel(&mut gas, &mut wall, None, Time::new(1.0), capped, || false)
            .is_err()
    );
    assert_eq!(gas, before);
    assert_eq!(wall, old_wall);
}

#[test]
fn g4_total_work_and_final_cancellation_are_atomic_and_replayable() {
    let initial = parcel(source(), 2.0 * MIXTURE_MASS);
    let original_wall = WallThermal::Reservoir {
        temperature_k: 300.0,
        conductance_w_per_k: 1.0,
    };
    let normal = controls([0.0, 1000.0]);
    for (max_heat, max_pt) in [(1, normal.maximum_pt_evaluations), (132, 1)] {
        let mut gas = initial.clone();
        let mut wall = original_wall;
        let limited = MixtureThermalControl {
            maximum_heat_evaluations: max_heat,
            maximum_pt_evaluations: max_pt,
            ..normal
        };
        let error =
            exchange_mixture_parcel(&mut gas, &mut wall, None, Time::new(10.0), limited, || {
                false
            })
            .unwrap_err();
        assert!(matches!(
            error,
            WallExchangeError::MixtureThermalWork { .. }
        ));
        assert_eq!(gas, initial);
        assert_eq!(wall, original_wall);
    }
    let mut expected = initial.clone();
    let mut expected_wall = original_wall;
    let mut polls = 0;
    let accepted = exchange_mixture_parcel(
        &mut expected,
        &mut expected_wall,
        None,
        Time::new(10.0),
        normal,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut gas = initial.clone();
    let mut wall = original_wall;
    let mut repeated = 0;
    assert_eq!(
        exchange_mixture_parcel(&mut gas, &mut wall, None, Time::new(10.0), normal, || {
            repeated += 1;
            repeated == polls
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(repeated, polls);
    assert_eq!(gas, initial);
    assert_eq!(wall, original_wall);
    let retried =
        exchange_mixture_parcel(&mut gas, &mut wall, None, Time::new(10.0), normal, || false)
            .unwrap();
    assert_eq!(gas, expected);
    assert_eq!(wall, expected_wall);
    assert_eq!(accepted, retried);
}
