//! Actual finite-gas accounts around changing spherical geometry.
use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use fs_thermal::evaporation::{EvaporationProperties, MOLAR_GAS_CONSTANT};
use fs_thermochem::SpeciesId;

const D: f64 = 2e-5;
fn initial() -> SphericalDropletState {
    SphericalDropletState {
        mass_kg: 1e-6,
        specific_enthalpy_j_kg: liquid(300.0).specific_enthalpy_j_kg,
    }
}
fn heat_controls() -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [
            liquid(200.0).specific_enthalpy_j_kg,
            liquid(900.0).specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn evaporation_and_condensation_update_size_and_conserve_original_species_and_energy() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let vapor = SpeciesId::new("B").unwrap();
    let carrier = SpeciesId::new("A").unwrap();
    for (partial_pressure, gas_t) in [(0.0, 300.0), (5000.0, 320.0)] {
        let original = gas(&source, gas_t, partial_pressure, 1e-4);
        let saved = original.clone();
        let state = initial();
        let dt = 0.1;
        let candidate = advance_droplet_with_gas(
            &Phase,
            &binding,
            state,
            &original,
            &source,
            temperature,
            D,
            dt,
            controls(),
            || false,
        )
        .unwrap();
        let [before, after] = *candidate.geometry();
        let step = candidate.step();
        let p = original.state().state().pressure;
        let pv = original.species_partial_pressure_pa(&vapor).unwrap();
        // Independent integrated stagnant-sphere diffusion, not a second call
        // to the production surface or coefficient helper.
        let request = dt
            * 4.0
            * core::f64::consts::PI
            * before.radius_m()
            * D
            * (0.018015 * p / (MOLAR_GAS_CONSTANT * gas_t))
            * ((p - pv) / (p - 1500.0)).ln();
        near(step.requested_transferred_mass_kg, request, 2e-22);
        assert_eq!(step.transferred_mass_kg > 0.0, partial_pressure == 0.0);
        assert_eq!(
            after.radius_m() < before.radius_m(),
            step.transferred_mass_kg > 0.0
        );
        assert_eq!(
            after.area_m2() < before.area_m2(),
            step.transferred_mass_kg > 0.0
        );
        near(
            after.radius_m(),
            (3.0 * candidate.droplet().mass_kg / (4.0 * core::f64::consts::PI * 1000.0)).cbrt(),
            2e-17,
        );
        near(
            candidate.gas().species_mass_kg(&vapor).unwrap() + candidate.droplet().mass_kg,
            original.species_mass_kg(&vapor).unwrap() + state.mass_kg,
            controls().gas.mass_tolerance_kg,
        );
        near(
            candidate.gas().species_mass_kg(&carrier).unwrap(),
            original.species_mass_kg(&carrier).unwrap(),
            controls().gas.mass_tolerance_kg,
        );
        near(
            step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
            0.0,
            controls().gas.energy_tolerance_j,
        );
        near(
            step.energy_residual_j,
            0.0,
            controls().gas.energy_tolerance_j,
        );
        assert_eq!(candidate.gas().volume_m3(), original.volume_m3());
        assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&11.0));
        assert!(std::ptr::eq(candidate.catalog(), &source));
        assert_eq!(candidate.phase_step(), candidate.step());
        assert_eq!(candidate.heat(), None);
        let next = advance_droplet_with_gas(
            &Phase,
            &binding,
            candidate.droplet(),
            candidate.gas(),
            &source,
            temperature,
            D,
            dt,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(next.geometry()[0], after);
        assert_ne!(
            next.report().mass_transfer_coefficient_m_s,
            candidate.report().mass_transfer_coefficient_m_s
        );
        assert_eq!(original, saved);
    }
}

struct ExpandingPhase;
impl EvaporationLaw for ExpandingPhase {
    fn temperature_range_k(&self) -> [f64; 2] {
        Phase.temperature_range_k()
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let mut value = Phase.properties(t)?;
        value.liquid_density_kg_m3 = 1000.0 / (1.0 + 0.001 * (t - 300.0));
        Ok(value)
    }
}

#[test]
fn sealed_diffusion_still_exchanges_reciprocal_heat_and_uses_final_material_density() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 340.0, 1000.0, 1e-4);
    let state = initial();
    let candidate = advance_droplet_with_gas_and_heat(
        &ExpandingPhase,
        &binding,
        state,
        &original,
        &source,
        temperature,
        0.0,
        80.0,
        0.1,
        controls(),
        heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate.droplet().mass_kg, state.mass_kg);
    assert_eq!(candidate.gas().mass_kg(), original.mass_kg());
    assert_eq!(candidate.gas().model(), original.model());
    assert_eq!(candidate.step().transferred_mass_kg, 0.0);
    assert_eq!(candidate.step().energy_into_gas_j, 0.0);
    assert_eq!(candidate.step().recession_m, 0.0);
    assert!(candidate.heat().unwrap().heat_into_liquid_j > 0.0);
    let final_t = candidate.step().surface_temperature_k[1];
    assert!(final_t > 300.0);
    let density = 1000.0 / (1.0 + 0.001 * (final_t - 300.0));
    assert_eq!(candidate.geometry()[1].liquid_density_kg_m3(), density);
    assert!(candidate.geometry()[1].radius_m() > candidate.geometry()[0].radius_m());
    near(
        candidate.geometry()[1].radius_m(),
        (3.0 * state.mass_kg / (4.0 * core::f64::consts::PI * density)).cbrt(),
        2e-17,
    );
    near(
        candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    near(
        candidate.step().energy_residual_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    assert_eq!(candidate.phase_step().surface_temperature_k[0], final_t);

    let phase = advance_droplet_with_gas(
        &ExpandingPhase,
        &binding,
        state,
        &original,
        &source,
        temperature,
        D,
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    let zero_heat = advance_droplet_with_gas_and_heat(
        &ExpandingPhase,
        &binding,
        state,
        &original,
        &source,
        temperature,
        D,
        0.0,
        0.1,
        controls(),
        heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(phase.droplet(), zero_heat.droplet());
    assert_eq!(phase.gas(), zero_heat.gas());
    assert_eq!(phase.geometry(), zero_heat.geometry());
    assert_eq!(phase.step(), zero_heat.step());
    assert_eq!(zero_heat.heat(), Some(&GasFilmHeatStep::default()));
}

#[test]
fn excessive_withdrawal_invalid_diffusion_and_cancellation_leave_inputs_available_to_retry() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 300.0, 0.0, 1e-4);
    let saved = original.clone();
    let state = initial();
    assert!(matches!(
        advance_droplet_with_gas(
            &Phase,
            &binding,
            state,
            &original,
            &source,
            temperature,
            D,
            1e6,
            controls(),
            || false,
        ),
        Err(GasEvaporationError::StepSizeLimit { .. })
    ));
    assert!(
        advance_droplet_with_gas(
            &Phase,
            &binding,
            state,
            &original,
            &source,
            temperature,
            -D,
            0.1,
            controls(),
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        advance_droplet_with_gas(
            &Phase,
            &binding,
            state,
            &original,
            &source,
            temperature,
            D,
            0.1,
            controls(),
            || true,
        ),
        Err(GasEvaporationError::Cancelled)
    ));
    assert_eq!(original, saved);
    assert_eq!(state, initial());
    assert!(
        advance_droplet_with_gas(
            &Phase,
            &binding,
            state,
            &original,
            &source,
            temperature,
            D,
            0.1,
            controls(),
            || false,
        )
        .is_ok()
    );
}

#[test]
fn explicit_sherwood_uses_original_surface_owner_and_retains_true_diffusivity() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let state = initial();
    for partial_pressure in [0.0, 5000.0] {
        let original = gas(&source, 300.0, partial_pressure, 1e-4);
        let saved = original.clone();
        let geometry = state.geometry(&Phase, 300.0).unwrap();
        let sherwood = 8.0;
        let coefficient = sherwood * D / geometry.diameter_m();
        let direct = advance_surface_with_gas(
            &Phase,
            &binding,
            LiquidSurfaceState {
                mass_kg: state.mass_kg,
                specific_enthalpy_j_kg: state.specific_enthalpy_j_kg,
                area_m2: geometry.area_m2(),
            },
            &original,
            &source,
            temperature,
            coefficient,
            0.1,
            controls(),
            || false,
        ).unwrap();
        let candidate = advance_droplet_with_gas_at_sherwood(
            &Phase, &binding, state, &original, &source, temperature,
            D, sherwood, 0.1, controls(), || false,
        ).unwrap();
        assert_eq!(candidate.report().binary_diffusivity_m2_s, D);
        assert_eq!(candidate.report().sherwood_number, sherwood);
        near(candidate.report().mass_transfer_coefficient_m_s, coefficient, 1e-16);
        assert_eq!(candidate.gas(), direct.gas());
        assert_eq!(candidate.droplet().mass_kg, direct.liquid().mass_kg);
        assert_eq!(
            candidate.droplet().specific_enthalpy_j_kg,
            direct.liquid().specific_enthalpy_j_kg,
        );
        assert_eq!(candidate.step(), direct.step());
        assert_eq!(candidate.step().transferred_mass_kg > 0.0, partial_pressure == 0.0);
        near(candidate.step().energy_residual_j, 0.0, controls().gas.energy_tolerance_j);
        assert_eq!(original, saved);

        let stagnant = advance_droplet_with_gas(
            &Phase, &binding, state, &original, &source, temperature,
            D, 0.1, controls(), || false,
        ).unwrap();
        let explicit_two = advance_droplet_with_gas_at_sherwood(
            &Phase, &binding, state, &original, &source, temperature,
            D, 2.0, 0.1, controls(), || false,
        ).unwrap();
        assert_eq!(stagnant, explicit_two);
        assert_eq!(
            stagnant.report().mass_transfer_coefficient_m_s.to_bits(),
            geometry.mass_transfer_coefficient_m_s(D).unwrap().to_bits(),
        );
        assert!(candidate.step().requested_transferred_mass_kg.abs()
            > stagnant.step().requested_transferred_mass_kg.abs());
    }
}

#[test]
fn sherwood_heat_path_keeps_independent_sensible_heat_and_refuses_invalid_transport() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let state = initial();
    let original = gas(&source, 320.0, 0.0, 1e-4);
    let saved = original.clone();
    let baseline = advance_droplet_with_gas_and_heat(
        &Phase, &binding, state, &original, &source, temperature,
        0.0, 80.0, 0.1, controls(), heat_controls(), || false,
    ).unwrap();
    let forced = advance_droplet_with_gas_and_heat_at_sherwood(
        &Phase, &binding, state, &original, &source, temperature,
        0.0, 8.0, 80.0, 0.1, controls(), heat_controls(), || false,
    ).unwrap();
    assert_eq!(forced.droplet(), baseline.droplet());
    assert_eq!(forced.gas(), baseline.gas());
    assert_eq!(forced.step(), baseline.step());
    assert_eq!(forced.heat(), baseline.heat());
    assert_eq!(forced.droplet().mass_kg, state.mass_kg);
    assert_eq!(forced.report().binary_diffusivity_m2_s, 0.0);
    assert_eq!(forced.report().sherwood_number, 8.0);
    assert_eq!(forced.report().mass_transfer_coefficient_m_s, 0.0);

    for sherwood in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(advance_droplet_with_gas_at_sherwood(
            &Phase, &binding, state, &original, &source, temperature,
            D, sherwood, 0.1, controls(), || false,
        ).is_err());
    }
    assert_eq!(original, saved);
}
