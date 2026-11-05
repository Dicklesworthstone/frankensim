//! G0/G1/G3: finite shared inventory, analytical pressure, storage and refusal.

use super::*;
use crate::sublimation::{MOLAR_GAS_CONSTANT, SublimationProperties};

const TEMPERATURE: f64 = 230.0;
const MOLAR_MASS: f64 = 0.018;
const LATENT: f64 = 2.8e6;
const SUBSTRATE_CP: f64 = 500.0;
const DEPOSIT_CP: f64 = 2000.0;
const VOLUME: f64 = 0.001;

// The emission coefficient times saturation pressure is constant. This
// separates an independent linear shared-pressure oracle from latent heating.
struct ManufacturedLaw;
impl SublimationLaw for ManufacturedLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let pressure = 10.0 * (t / TEMPERATURE).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure / (2.0 * t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}
static LAW: ManufacturedLaw = ManufacturedLaw;

struct ClausiusLaw;
impl SublimationLaw for ClausiusLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let exponent = 500_000.0 * MOLAR_MASS / MOLAR_GAS_CONSTANT;
        let pressure = 100.0 * (exponent * (1.0 / 250.0 - 1.0 / t)).exp();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * exponent / t / t,
            latent_heat_j_kg: 500_000.0,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}
static CLAUSIUS: ClausiusLaw = ClausiusLaw;

fn deposit_enthalpy(t: f64) -> Result<f64, SublimationError> {
    Ok(DEPOSIT_CP * t)
}
fn surface(area: f64, mass: f64, t: f64) -> SupportedChamberSurface<'static> {
    SupportedChamberSurface {
        state: SupportedSurfaceState {
            substrate_mass_kg: 0.02,
            substrate_specific_enthalpy_j_kg: SUBSTRATE_CP * t,
            deposit_mass_kg: mass,
            area_m2: area,
        },
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 1.0,
            mass_transfer_coefficient_m_s: None,
        },
        law: &LAW,
        deposit_enthalpy: &deposit_enthalpy,
    }
}
fn reservoir(pressure: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: pressure / beta(),
        volume_m3: VOLUME,
        temperature_k: TEMPERATURE,
    }
}
fn beta() -> f64 {
    MOLAR_GAS_CONSTANT * TEMPERATURE / MOLAR_MASS / VOLUME
}
fn mass(reservoir: VaporReservoir) -> f64 {
    let VaporReservoir::Finite { vapor_mass_kg, .. } = reservoir else {
        panic!("expected finite reservoir")
    };
    vapor_mass_kg
}
fn config() -> ChamberSolverConfig {
    ChamberSolverConfig {
        surface_solver: crate::sublimation::SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 256,
        },
        relative_mass_tolerance: 1e-12,
        absolute_mass_tolerance_kg: 1e-18,
        max_iterations: 128,
        max_surfaces: 8,
    }
}
fn run(
    surfaces: &[SupportedChamberSurface<'_>],
    bath: VaporReservoir,
    dt: f64,
) -> SupportedChamberStep {
    advance_supported_chamber(
        surfaces,
        bath,
        dt,
        config(),
        |_, h| Ok(h / SUBSTRATE_CP),
        || Ok(()),
    )
    .unwrap()
}
fn close(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn conservation(
    old: &[SupportedChamberSurface<'_>],
    bath: VaporReservoir,
    step: &SupportedChamberStep,
) {
    let old_mass = mass(bath) + old.iter().map(|s| s.state.deposit_mass_kg).sum::<f64>();
    let new_mass = mass(step.reservoir)
        + step
            .surfaces
            .iter()
            .map(|s| s.state.deposit_mass_kg)
            .sum::<f64>();
    close(new_mass, old_mass, 1e-18, 1e-13);
    close(
        step.substrate_energy_change_j + step.deposit_energy_change_j + step.vapor_energy_j,
        0.0,
        5e-9,
        0.0,
    );
    close(step.energy_residual_j, 0.0, 5e-9, 0.0);
    close(step.mass_residual_kg, 0.0, 1e-18, 0.0);
    for (initial, final_state) in old.iter().zip(&step.surfaces) {
        assert_eq!(
            final_state.state.substrate_mass_kg,
            initial.state.substrate_mass_kg
        );
        assert_eq!(final_state.phase.reservoir, step.reservoir);
        close(
            final_state.state.deposit_mass_kg - initial.state.deposit_mass_kg
                + final_state.phase.transferred_mass_kg,
            0.0,
            1e-18,
            0.0,
        );
    }
}

#[test]
fn g1_bare_substrates_share_one_inventory_and_match_independent_pressure_algebra() {
    let mut surfaces = [
        surface(0.001, 0.0, TEMPERATURE),
        surface(0.003, 0.0, TEMPERATURE),
    ];
    surfaces[1].kinetics.accommodation_coefficient = 0.4;
    surfaces[1].kinetics.mass_transfer_coefficient_m_s = Some(0.02);
    let dt = 0.02;
    let bath = reservoir(25.0);
    // Independent Hertz--Knudsen plus dilute-film series resistance.
    let g = surfaces.map(|s| {
        let kinetic = s.kinetics.accommodation_coefficient
            * (MOLAR_MASS / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * TEMPERATURE))
                .sqrt();
        let coefficient = match s.kinetics.mass_transfer_coefficient_m_s {
            None => kinetic,
            Some(h) => 1.0 / (1.0 / kinetic + MOLAR_GAS_CONSTANT * TEMPERATURE / (h * MOLAR_MASS)),
        };
        s.state.area_m2 * coefficient
    });
    let rate = beta() * (g[0] + g[1]);
    let expected_pressure = (25.0 + dt * rate * 10.0) / (1.0 + dt * rate);
    let step = run(&surfaces, bath, dt);
    close(
        step.reservoir.partial_pressure_pa(MOLAR_MASS).unwrap(),
        expected_pressure,
        1e-9,
        1e-11,
    );
    for (index, result) in step.surfaces.iter().enumerate() {
        let deposit = dt * g[index] * (expected_pressure - 10.0);
        close(result.state.deposit_mass_kg, deposit, 1e-18, 1e-10);
        close(
            result.phase.surface_temperature_k,
            TEMPERATURE + deposit * LATENT / (0.02 * SUBSTRATE_CP),
            2e-10,
            0.0,
        );
        assert!(result.state.deposit_mass_kg > 0.0);
        assert!(result.phase.recession_m < 0.0);
    }
    conservation(&surfaces, bath, &step);
}

#[test]
fn g1_warm_coating_supplies_bare_cold_substrate_in_initial_vacuum() {
    let mut surfaces = [surface(0.001, 1e-3, 260.0), surface(0.001, 0.0, 225.0)];
    for surface in &mut surfaces {
        surface.law = &CLAUSIUS;
    }
    let bath = reservoir(0.0);
    let step = run(&surfaces, bath, 0.1);
    assert!(step.surfaces[0].phase.transferred_mass_kg > 0.0);
    assert!(step.surfaces[1].phase.transferred_mass_kg < 0.0);
    assert!(step.surfaces[1].state.deposit_mass_kg > 0.0);
    assert!(step.surfaces[0].phase.surface_temperature_k < 260.0);
    assert!(step.surfaces[1].phase.surface_temperature_k > 225.0);
    assert!(mass(step.reservoir) > 0.0);
    conservation(&surfaces, bath, &step);
    let reversed = run(&[surfaces[1], surfaces[0]], bath, 0.1);
    assert_eq!(step.reservoir, reversed.reservoir);
    assert_eq!(step.surfaces[0], reversed.surfaces[1]);
    assert_eq!(step.surfaces[1], reversed.surfaces[0]);
}

#[test]
fn g0_bare_vacuum_is_stationary_and_exhaustion_never_spends_substrate_mass() {
    let bare = [
        surface(0.001, 0.0, TEMPERATURE),
        surface(0.002, 0.0, TEMPERATURE),
    ];
    let step = run(&bare, reservoir(0.0), 1e4);
    assert_eq!(step.reservoir, reservoir(0.0));
    for (initial, result) in bare.iter().zip(&step.surfaces) {
        assert_eq!(initial.state, result.state);
        assert_eq!(result.phase.transferred_mass_kg, 0.0);
    }
    let thin = [surface(0.001, 1e-10, TEMPERATURE)];
    let exhausted = run(&thin, reservoir(0.0), 1e4);
    assert_eq!(exhausted.surfaces[0].state.deposit_mass_kg, 0.0);
    assert_eq!(exhausted.surfaces[0].state.substrate_mass_kg, 0.02);
    close(mass(exhausted.reservoir), 1e-10, 1e-25, 0.0);
    conservation(&thin, reservoir(0.0), &exhausted);
}

#[test]
fn g3_independent_substrate_references_and_common_deposit_reference_cancel() {
    let original = [surface(0.001, 0.0, 230.0), surface(0.002, 2e-5, 240.0)];
    let bath = reservoir(30.0);
    let baseline = run(&original, bath, 0.05);
    let offsets = [1e6, -3e6];
    let deposit_offset = -4e6;
    let shifted_deposit = |t| Ok(DEPOSIT_CP * t + deposit_offset);
    let mut shifted = original;
    for (index, surface) in shifted.iter_mut().enumerate() {
        surface.state.substrate_specific_enthalpy_j_kg += offsets[index];
        surface.deposit_enthalpy = &shifted_deposit;
    }
    let result = advance_supported_chamber(
        &shifted,
        bath,
        0.05,
        config(),
        |index, h| Ok((h - offsets[index]) / SUBSTRATE_CP),
        || Ok(()),
    )
    .unwrap();
    close(
        mass(result.reservoir),
        mass(baseline.reservoir),
        1e-18,
        1e-10,
    );
    for (a, b) in result.surfaces.iter().zip(&baseline.surfaces) {
        close(
            a.phase.surface_temperature_k,
            b.phase.surface_temperature_k,
            1e-8,
            0.0,
        );
        close(
            a.state.deposit_mass_kg,
            b.state.deposit_mass_kg,
            1e-18,
            1e-10,
        );
    }
    close(
        result.vapor_energy_j - baseline.vapor_energy_j,
        result.transferred_mass_kg * deposit_offset,
        1e-8,
        0.0,
    );
    conservation(&shifted, bath, &result);
}

#[test]
fn g1_pressure_history_converges_at_first_order_without_spending_vapor_twice() {
    let initial = [
        surface(0.001, 0.0, TEMPERATURE),
        surface(0.002, 0.0, TEMPERATURE),
    ];
    let duration = 0.002;
    let conductance = 0.003
        * (MOLAR_MASS / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * TEMPERATURE)).sqrt();
    let exact = 10.0 + 15.0 * (-beta() * conductance * duration).exp();
    let mut errors = Vec::new();
    for count in [4, 8, 16, 32] {
        let mut surfaces = initial;
        let mut bath = reservoir(25.0);
        for _ in 0..count {
            let step = run(&surfaces, bath, duration / count as f64);
            conservation(&surfaces, bath, &step);
            for (surface, result) in surfaces.iter_mut().zip(&step.surfaces) {
                surface.state = result.state;
            }
            bath = step.reservoir;
        }
        errors.push((bath.partial_pressure_pa(MOLAR_MASS).unwrap() - exact).abs());
    }
    assert!(errors[0] > 1e-3);
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.8..2.1).contains(&ratio),
            "errors={errors:?}, ratio={ratio}"
        );
    }
}

#[test]
fn g0_caps_bad_species_and_cancellation_refuse_atomically() {
    let surfaces = [
        surface(0.001, 0.0, TEMPERATURE),
        surface(0.002, 0.0, TEMPERATURE),
    ];
    let bath = reservoir(20.0);
    for control in [
        ChamberSolverConfig {
            max_surfaces: 1,
            ..config()
        },
        ChamberSolverConfig {
            max_iterations: 1,
            ..config()
        },
    ] {
        assert!(
            advance_supported_chamber(
                &surfaces,
                bath,
                0.1,
                control,
                |_, h| Ok(h / SUBSTRATE_CP),
                || Ok(()),
            )
            .is_err()
        );
    }
    let mut checks = 0;
    let cancelled = advance_supported_chamber(
        &surfaces,
        bath,
        0.1,
        config(),
        |_, h| Ok(h / SUBSTRATE_CP),
        || {
            checks += 1;
            if checks == 12 {
                Err(SublimationError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(cancelled, Err(SublimationError::Cancelled));
    assert_eq!(surfaces[0].state.deposit_mass_kg, 0.0);
    assert_eq!(bath, reservoir(20.0));

    struct OtherSpecies;
    impl SublimationLaw for OtherSpecies {
        fn temperature_range_k(&self) -> [f64; 2] {
            [100.0, 400.0]
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            Ok(SublimationProperties {
                molar_mass_kg_mol: 0.044,
                ..LAW.properties(t)?
            })
        }
    }
    let mut incompatible = surfaces;
    incompatible[1].law = &OtherSpecies;
    assert!(matches!(
        advance_supported_chamber(
            &incompatible,
            bath,
            0.1,
            config(),
            |_, h| Ok(h / SUBSTRATE_CP),
            || Ok(()),
        ),
        Err(SublimationError::Invalid(_))
    ));
    let empty = run(&[], bath, 1.0);
    assert_eq!(empty.reservoir, bath);
    assert!(empty.surfaces.is_empty());
    assert_eq!(empty.pressure_iterations, 0);
}
