//! Manufactured original PR/NASA sources: physical identities and refinement,
//! not empirical validation of a particular gas or vessel.

use super::*;
use crate::thermochemical_gas::energy::GasEnergyError;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::nonideal::energy::tests::bound_with_reference;
use crate::thermochemical_gas::nonideal::energy::parcel::transfer::RealGasMassTransferSolve;
use crate::thermochemical_gas::tests::{model_in_domains, point};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

fn parcel(reference: f64) -> IsochoricRealGasParcel {
    IsochoricRealGasParcel::new(
        bound_with_reference(reference),
        point(300.0, 5.0e6).with("material-lot", 12.0).unwrap(),
        0.01, PengRobinsonPhaseV1::Vapor,
    ).unwrap()
}
fn controls() -> RealGasNozzleBlowdownSolve {
    let entropy = GasEntropySolve {
        minimum_temperature_k: 220.0, maximum_temperature_k: 450.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128,
    };
    let energy = GasEnergySolve {
        minimum_temperature_k: 220.0, maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
    };
    RealGasNozzleBlowdownSolve {
        nozzle: RealGasNozzleStepSolve {
            entropy,
            transfer: RealGasMassTransferSolve {
                donor: energy, receiver: energy,
                mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-7,
            },
            minimum_density_fraction: 0.5,
            pressure_tolerance_pa: 1e-5, sonic_tolerance_m2_s2: 1e-5,
            maximum_iterations: 128, maximum_donor_fraction: 0.1,
        },
        reservoir: RealGasBlowdownSolve { entropy, mass_tolerance_kg: 1e-12 },
        maximum_iterations: 128,
    }
}
fn nozzle() -> GasNozzle {
    GasNozzle { area_m2: 1e-6, discharge_coefficient: 0.8 }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance={tolerance:.3e}");
}

#[test]
fn g1_dense_finite_exhaust_matches_an_independent_isentrope_and_energy_oracle() {
    let original = parcel(0.0);
    let saved = original.clone();
    let mass = original.mass_kg() * (1.0 - 1.0 / 1.1);
    let candidate = original.exhaust_mass_adiabatic(
        mass, controls().reservoir, || false,
    ).unwrap();
    let endpoint = candidate.endpoint();
    let step = candidate.step();
    // Independent 75-digit PR pressure/NASA-Cv entropy integration already
    // used by the volume-work oracle: rho1=rho0/1.1, -m0*Delta u=4658.361642714088 J.
    close(endpoint.state().temperature_k(), 287.505_188_450_419_14, 2e-8);
    close(endpoint.state().pressure_pa(), 4_343_285.645_717_975, 1e-3);
    close(
        step.exhausted_enthalpy_j - mass * original.state().specific_internal_energy_j_kg(),
        4_658.361_642_714_088 / 1.1,
        3e-5,
    );
    close(endpoint.mass_kg() + mass, original.mass_kg(), 1e-14);
    assert_eq!(endpoint.volume_m3(), original.volume_m3());
    assert_eq!(endpoint.density_kg_m3(), endpoint.mass_kg() / endpoint.volume_m3());
    assert_eq!(endpoint.model(), original.model());
    assert_eq!(endpoint.point().axes().get("material-lot"), Some(&12.0));
    assert_eq!(original, saved);
    assert_eq!(step.species.len(), 2);
    // Original mole fractions .75/.25 and molar masses .028/.004 give these
    // mass fractions, rather than naively using mole fractions as mass fluxes.
    close(step.species[0].mass_kg, mass * 0.75 * 0.028 / 0.022, 1e-14);
    close(step.species[1].mass_kg, mass * 0.25 * 0.004 / 0.022, 1e-14);
    assert!(step.mass_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg);
    assert!(step.maximum_species_mass_residual_kg <= controls().reservoir.mass_tolerance_kg);
    assert!(step.specific_entropy_residual_j_kg_k.abs()
        <= controls().reservoir.entropy.entropy_tolerance_j_kg_k);
    let mean = step.mean_stagnation_enthalpy_j_kg.unwrap();
    assert!(mean < original.state().specific_enthalpy_j_kg());
    assert!(mean > endpoint.state().specific_enthalpy_j_kg());
    // Finite withdrawal must not export an unchanged initial enthalpy.
    assert!((step.exhausted_enthalpy_j
        - mass * original.state().specific_enthalpy_j_kg()).abs() > 100.0);
}

fn dilute_constant_cp_parcel() -> IsochoricRealGasParcel {
    let source = model_in_domains(
        vec![0.75, 0.25], 0.0, [200.0, 2000.0], [200.0, 2000.0], [0.01, 1e7],
    );
    let model = ThermochemicalPengRobinsonModel::new(
        source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6), acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6), acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    ).unwrap();
    IsochoricRealGasParcel::new(
        model, point(300.0, 0.1), 0.01, PengRobinsonPhaseV1::Vapor,
    ).unwrap()
}

#[test]
fn g1_g3_nozzle_blowdown_converges_quadratically_to_the_analytic_dilute_limit() {
    let original = dilute_constant_cp_parcel();
    let gamma: f64 = 3.25 / 2.25;
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.022;
    let geometry = nozzle();
    let rate_over_initial_mass = geometry.discharge_coefficient * geometry.area_m2
        / original.volume_m3() * (gamma * r * 300.0).sqrt()
        * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let power = 0.5 * (gamma - 1.0);
    let duration = 4.0;
    // dm/dt=-k*m0*(m/m0)^((gamma+1)/2), solved independently in closed form.
    let expected_ratio = (1.0 + power * rate_over_initial_mass * duration).powf(-1.0 / power);
    let mut solve = controls();
    solve.reservoir.mass_tolerance_kg = 1e-20;
    let mut errors = Vec::new();
    for intervals in [1_u32, 2, 4] {
        let mut state = original.clone();
        let mut expelled = 0.0;
        let mut energy = 0.0;
        for _ in 0..intervals {
            let candidate = state.nozzle_blowdown(
                0.0, geometry, duration / f64::from(intervals), solve, || false,
            ).unwrap();
            assert!(candidate.midpoint_flow().choked);
            assert!(candidate.mass_quadrature_residual_kg.abs()
                <= solve.reservoir.mass_tolerance_kg);
            let (next, step, _) = candidate.into_parts();
            expelled += step.expelled_mass_kg;
            energy += step.exhausted_enthalpy_j;
            state = next;
        }
        let ratio = state.mass_kg() / original.mass_kg();
        errors.push((ratio - expected_ratio).abs());
        close(state.mass_kg() + expelled, original.mass_kg(), 2e-20);
        close(
            energy,
            original.mass_kg() * original.state().specific_internal_energy_j_kg()
                - state.mass_kg() * state.state().specific_internal_energy_j_kg(),
            1e-12,
        );
        assert!((state.state().temperature_k() / (300.0 * ratio.powf(gamma - 1.0)) - 1.0).abs() < 3e-8);
        assert!((state.state().pressure_pa() / (0.1 * ratio.powf(gamma)) - 1.0).abs() < 3e-8);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.8..4.2).contains(&ratio), "time-refinement errors {errors:?}");
    }
    assert!(errors[2] < 3e-6);
}

#[test]
fn g1_dense_midpoint_rate_tracks_pressure_cooling_and_choking() {
    let gas = parcel(0.0);
    let solve = controls();
    // Independent scalar PR/NASA entropy/enthalpy calculation, solving the
    // midpoint mass equation separately for each physical back pressure.
    for (back_pressure, expected_mass, expected_t, expected_p) in [
        (0.0, 0.001_715_920_040_081_015_2, 299.519_004_941_611_34, 4_973_468.564_275_586),
        (1e6, 0.001_715_920_040_081_015_2, 299.519_004_941_611_34, 4_973_468.564_275_586),
        (4.8e6, 0.000_677_000_806_694_197_6, 299.810_334_751_982_05, 4_989_525.934_778_947),
    ] {
        let initial = gas.nozzle_flow_to(back_pressure, nozzle(), solve.nozzle, || false).unwrap();
        let candidate = gas.nozzle_blowdown(back_pressure, nozzle(), 0.2, solve, || false).unwrap();
        let endpoint = candidate.reservoir().endpoint();
        let step = candidate.reservoir().step();
        close(step.expelled_mass_kg, expected_mass, 2e-12);
        close(endpoint.state().temperature_k(), expected_t, 3e-8);
        close(endpoint.state().pressure_pa(), expected_p, 2e-3);
        assert!(endpoint.mass_kg() < gas.mass_kg());
        assert!(endpoint.state().temperature_k() < gas.state().temperature_k());
        assert!(endpoint.state().pressure_pa() < gas.state().pressure_pa());
        assert!(endpoint.state().pressure_pa() > back_pressure);
        assert!(step.expelled_mass_kg < 0.2 * initial.mass_flow_kg_s);
        assert_eq!(candidate.midpoint_flow().choked, back_pressure < 4e6);
        close(step.expelled_mass_kg, 0.2 * candidate.midpoint_flow().mass_flow_kg_s,
            solve.reservoir.mass_tolerance_kg);
        assert!(candidate.mass_quadrature_residual_kg.abs()
            <= solve.reservoir.mass_tolerance_kg);
    }
}

#[test]
fn g3_large_formation_references_do_not_change_depletion_or_pressure_driven_time() {
    let original = parcel(0.0);
    for back_pressure in [0.0, 1e6, 4.8e6] {
        let plain = original.nozzle_blowdown(back_pressure, nozzle(), 0.2, controls(), || false).unwrap();
        for reference in [-1e24, 1e24] {
            let shifted = parcel(reference).nozzle_blowdown(
                back_pressure, nozzle(), 0.2, controls(), || false,
            ).unwrap();
            let a = plain.reservoir().endpoint();
            let b = shifted.reservoir().endpoint();
            assert_eq!(a.mass_kg(), b.mass_kg());
            assert_eq!(a.density_kg_m3(), b.density_kg_m3());
            assert_eq!(a.state().temperature_k(), b.state().temperature_k());
            assert_eq!(a.state().pressure_pa(), b.state().pressure_pa());
            assert_eq!(plain.midpoint_flow().mass_flow_kg_s, shifted.midpoint_flow().mass_flow_kg_s);
            assert_eq!(plain.midpoint_flow().velocity_m_s, shifted.midpoint_flow().velocity_m_s);
            assert_eq!(plain.mass_quadrature_residual_kg, shifted.mass_quadrature_residual_kg);
            assert_eq!(plain.iterations, shifted.iterations);
            assert_eq!(plain.reservoir().step().species, shifted.reservoir().step().species);
            assert_ne!(plain.reservoir().step().exhausted_enthalpy_j,
                shifted.reservoir().step().exhausted_enthalpy_j);
        }
    }
}

#[test]
fn g4_zero_flow_refusals_and_late_cancellation_publish_no_partial_depletion() {
    let gas = parcel(0.0);
    let saved = gas.clone();
    let solve = controls();
    let identity = gas.exhaust_mass_adiabatic(0.0, solve.reservoir, || false).unwrap();
    assert_eq!(identity.endpoint(), &gas);
    assert_eq!(identity.step().exhausted_enthalpy_j, 0.0);
    for (back, area) in [(gas.state().pressure_pa(), 1e-6), (0.0, 0.0)] {
        let candidate = gas.nozzle_blowdown(
            back, GasNozzle { area_m2: area, ..nozzle() }, 1.0, solve, || false,
        ).unwrap();
        assert_eq!(candidate.reservoir().endpoint(), &gas);
        assert_eq!(candidate.reservoir().step().expelled_mass_kg, 0.0);
        assert_eq!(candidate.iterations, 0);
    }
    for mass in [-1.0, gas.mass_kg(), f64::NAN, 0.5 * solve.reservoir.mass_tolerance_kg] {
        assert!(gas.exhaust_mass_adiabatic(mass, solve.reservoir, || false).is_err());
    }
    let mut narrow = solve.reservoir;
    narrow.entropy.minimum_temperature_k = 295.0;
    assert!(gas.exhaust_mass_adiabatic(0.09 * gas.mass_kg(), narrow, || false).is_err());
    assert!(gas.nozzle_blowdown(1e6, nozzle(), 100.0, solve, || false).is_err());
    assert!(gas.nozzle_blowdown(5.1e6, nozzle(), 0.2, solve, || false).is_err());
    let crossing = gas.nozzle_blowdown(4.8e6, nozzle(), 6.0, solve, || false);
    assert!(matches!(crossing, Err(RealGasBlowdownError::Invalid(
        "real-gas blowdown crosses back-pressure equilibrium; reduce the interval"
    ))), "{crossing:?}");
    assert!(matches!(
        gas.nozzle_blowdown(1e6, nozzle(), 0.2,
            RealGasNozzleBlowdownSolve { maximum_iterations: 1, ..solve }, || false),
        Err(RealGasBlowdownError::NotConverged { iterations: 1, .. })
    ));
    let mut polls = 0;
    let expected = gas.nozzle_blowdown(1e6, nozzle(), 0.2, solve, || {
        polls += 1;
        false
    }).unwrap();
    let mut seen = 0;
    assert!(matches!(
        gas.nozzle_blowdown(1e6, nozzle(), 0.2, solve, || {
            seen += 1;
            seen == polls
        }),
        Err(RealGasBlowdownError::Energy(GasEnergyError::Cancelled))
    ));
    assert_eq!(seen, polls);
    assert_eq!(gas, saved);
    assert_eq!(gas.nozzle_blowdown(1e6, nozzle(), 0.2, solve, || false).unwrap(), expected);
}
