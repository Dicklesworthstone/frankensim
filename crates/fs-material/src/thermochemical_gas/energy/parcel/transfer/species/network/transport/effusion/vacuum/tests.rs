//! G1/G3/G4: molecular escape, caloric cooling, refinement and atomic admission.
//! Synthetic gas coefficients test the equations, not empirical gas accuracy.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::{
        effusion::{EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
        network::{
            SpeciesNetworkError, SpeciesNetworkSolve, transport::effusion::EffusionNetworkLink,
        },
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.028;
const VOLUME: f64 = 0.001;
const AREA: f64 = 1e-7;

fn model(rows: &[(&str, f64, f64)], fractions: Vec<f64>, datum: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .map(|&(id, m, cp)| species_with_energy_reference(id, m, cp, 0.0, datum))
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _, _)| transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1.0, 1e7]))
            .collect(),
        None,
    )
    .unwrap()
}

fn pure(cp_r: f64, datum: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model(&[("A", M, cp_r)], vec![1.0], datum),
        point(t, p),
        VOLUME,
    )
    .unwrap()
}

fn support() -> EffusiveSupport {
    EffusiveSupport {
        temperature_k: (200.0, 1800.0),
        pressure_pa: (1.0, 1e7),
    }
}

fn interface(
    gas: &ThermochemicalGasModel,
    domain: EffusiveSupport,
    probability: f64,
) -> EffusiveInterface {
    let rows = gas
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                probability,
                "synthetic reciprocal energy-independent aperture",
            )
            .unwrap()
        })
        .collect();
    EffusiveInterface::try_new(
        gas.clone(),
        rows,
        domain,
        EffusiveRegime::try_new(
            1e-8,
            1e-11,
            1e-4,
            "synthetic thin-array collisionless support; no measured regime claim",
        )
        .unwrap(),
    )
    .unwrap()
}

fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases
        .keys()
        .map(|&id| {
            (
                id,
                GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
            )
        })
        .collect()
}

fn controls() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}

fn outlet(interface: &EffusiveInterface, area_m2: f64) -> VacuumEffusionLink<'_> {
    VacuumEffusionLink {
        id: 7,
        parcel: 0,
        interface,
        area_m2,
    }
}

fn rate(area: f64, p: f64, mass: f64, t: f64) -> f64 {
    area * p / (2.0 * core::f64::consts::PI * mass * R * t).sqrt()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn vacuum_monatomic_and_diatomic_finite_cooling_is_reference_invariant() {
    for cp_r in [2.5, 3.5] {
        let mut temperatures = Vec::new();
        for datum in [0.0, 1e12] {
            let gases = BTreeMap::from([(0, pure(cp_r, datum, 400.0, 200.0))]);
            let saved = gases.clone();
            let surface = interface(gases[&0].model(), support(), 1.0);
            let result = step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[outlet(&surface, AREA)],
                1.0,
                &endpoints(&gases),
                controls(),
                || false,
            )
            .unwrap();
            let escaped_moles = rate(AREA, 200.0, M, 400.0);
            let escaped_mass = M * escaped_moles;
            let n0 = 200.0 * VOLUME / (R * 400.0);
            let f = escaped_moles / n0;
            // Molecular energy is (Cv+R/2)T, hence 2RT or 3RT per mole.
            let factor = (cp_r - 0.5) / (cp_r - 1.0);
            let expected_t = 400.0 * (1.0 - factor * f) / (1.0 - f);
            let gas = &result.parcels()[&0];
            near(gas.mass_kg(), (n0 - escaped_moles) * M, 1e-18);
            near(gas.state().state().temperature, expected_t, 2e-8);
            near(
                gas.state().state().pressure,
                200.0 * (1.0 - factor * f),
                2e-8,
            );
            assert_eq!(gas.volume_m3(), VOLUME);
            let flux = &result.step().outlets[0].flux;
            near(flux.molar_rate_mol_s, escaped_moles, 1e-18);
            near(flux.mass_rate_kg_s, escaped_mass, 1e-20);
            near(
                flux.additional_energy_rate_w,
                -0.5 * escaped_moles * R * 400.0,
                1e-13,
            );
            let boundary = &result.step().exchange.boundaries[0];
            assert_eq!(boundary.transfer.incoming_mass_kg, 0.0);
            near(boundary.transfer.outgoing_mass_kg, escaped_mass, 1e-20);
            near(boundary.transfer.net_mass_into_gas_kg, -escaped_mass, 1e-20);
            near(
                boundary.transfer.relative_energy_into_gas_j,
                0.5 * escaped_moles * R * 400.0,
                1e-13,
            );
            if datum == 0.0 {
                let energy = escaped_moles * (cp_r - 0.5) * R * 400.0;
                near(flux.energy_rate_w, energy, 1e-12);
                near(boundary.energy_into_gas_j, -energy, 1e-12);
            }
            near(result.step().exchange.total_mass_residual_kg, 0.0, 1e-13);
            near(result.step().exchange.total_energy_residual_j, 0.0, 1e-8);
            assert_eq!(gases, saved);
            temperatures.push(gas.state().state().temperature);
        }
        near(temperatures[0], temperatures[1], 2e-8);
    }
}

#[test]
fn vacuum_uses_original_temperature_dependent_calorics() {
    let (cp_r, slope, t0, p0) = (3.5, 0.01, 1200.0, 200.0);
    let n0 = p0 * VOLUME / (R * t0);
    let area = 0.08 * n0 / rate(1.0, p0, M, t0);
    let loss = rate(area, p0, M, t0);
    let f = loss / n0;
    // NASA h/R=aT+bT^2/2 and u/R=(a-1)T+bT^2/2. Solve the
    // independent quadratic caloric equation after removing molecular energy.
    let u0_over_r = (cp_r - 1.0) * t0 + 0.5 * slope * t0 * t0;
    let crossing_over_r = (cp_r - 0.5) * t0 + 0.5 * slope * t0 * t0;
    let target = (u0_over_r - f * crossing_over_r) / (1.0 - f);
    let expected_t =
        2.0 * target / ((cp_r - 1.0) + ((cp_r - 1.0).powi(2) + 2.0 * slope * target).sqrt());
    let frozen_cv_t = t0 - 0.5 * f * t0 / ((1.0 - f) * (cp_r - 1.0 + slope * t0));
    assert!((expected_t - frozen_cv_t).abs() > 1e-3);
    for datum in [0.0, 1e12] {
        let model = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species_with_energy_reference("A", M, cp_r, slope, datum)],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![transport_in_domain(
                "A",
                1.7e-5,
                0.025,
                [200.0, 2000.0],
                [1.0, 1e7],
            )],
            None,
        )
        .unwrap();
        let gases = BTreeMap::from([(
            0,
            ClosedGasParcel::new(model, point(t0, p0), VOLUME).unwrap(),
        )]);
        let surface = interface(gases[&0].model(), support(), 1.0);
        let result = step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[outlet(&surface, area)],
            1.0,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        near(
            result.parcels()[&0].state().state().temperature,
            expected_t,
            2e-8,
        );
        near(result.parcels()[&0].mass_kg(), (n0 - loss) * M, 1e-18);
        near(result.step().exchange.total_energy_residual_j, 0.0, 1e-8);
        if datum == 0.0 {
            near(
                result.step().outlets[0].flux.energy_rate_w,
                loss * R * crossing_over_r,
                1e-12,
            );
        }
    }
}

#[test]
fn vacuum_fractionates_mixtures_and_closes_each_species_and_energy() {
    let mixture = model(&[("A", M, 3.5), ("B", 0.004, 2.5)], vec![0.5, 0.5], 0.0);
    let gases = BTreeMap::from([(
        0,
        ClosedGasParcel::new(mixture, point(400.0, 200.0), VOLUME).unwrap(),
    )]);
    let surface = interface(gases[&0].model(), support(), 1.0);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[outlet(&surface, AREA)],
        1.0,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let original_moles = 100.0 * VOLUME / (R * 400.0);
    let loss = [rate(AREA, 100.0, M, 400.0), rate(AREA, 100.0, 0.004, 400.0)];
    let final_moles = [original_moles - loss[0], original_moles - loss[1]];
    let escaped_energy = R * 400.0 * (3.0 * loss[0] + 2.0 * loss[1]);
    let stored_energy = R * 400.0 * original_moles * (2.5 + 1.5);
    let expected_t =
        (stored_energy - escaped_energy) / (R * (2.5 * final_moles[0] + 1.5 * final_moles[1]));
    let gas = &result.parcels()[&0];
    for (i, (id, mass)) in [("A", M), ("B", 0.004)].into_iter().enumerate() {
        near(
            gas.species_mass_kg(&SpeciesId::new(id).unwrap()).unwrap(),
            final_moles[i] * mass,
            1e-18,
        );
    }
    near(gas.state().state().temperature, expected_t, 2e-8);
    near(
        result.step().outlets[0].flux.energy_rate_w,
        escaped_energy,
        1e-12,
    );
    near(
        result.step().exchange.parcels[0].energy_change_j,
        -escaped_energy,
        1e-9,
    );
    near(loss[1] / loss[0], (M / 0.004).sqrt(), 1e-12);
    assert!(gas.model().thermodynamics().mole_composition().fractions()[0] > 0.5);
    assert_eq!(result.step().exchange.boundaries.len(), 2);
}

#[test]
fn vacuum_time_refinement_converges_to_adiabatic_blowdown() {
    let initial = pure(2.5, 0.0, 600.0, 200.0);
    let m0 = initial.mass_kg();
    let surface = interface(initial.model(), support(), 1.0);
    let k0 = AREA / VOLUME * (R * 600.0 / (2.0 * core::f64::consts::PI * M)).sqrt();
    let duration = 0.2 / k0;
    // Cv=3R/2: T/T0=(1+k0*t/6)^-2, m/m0=(1+k0*t/6)^-6.
    let denominator = 1.0 + k0 * duration / 6.0;
    let expected_t = 600.0 / denominator.powi(2);
    let expected_mass = m0 / denominator.powi(6);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut gases = BTreeMap::from([(0, initial.clone())]);
        let solves = endpoints(&gases);
        for _ in 0..steps {
            gases = step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[outlet(&surface, AREA)],
                duration / f64::from(steps),
                &solves,
                controls(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0;
        }
        errors.push(
            (gases[&0].state().state().temperature / expected_t - 1.0)
                .abs()
                .max((gases[&0].mass_kg() / expected_mass - 1.0).abs()),
        );
    }
    assert!(errors[2] < 5e-4, "{errors:?}");
    for pair in errors.windows(2) {
        assert!(
            pair[1] > 0.45 * pair[0] && pair[1] < 0.55 * pair[0],
            "{errors:?}"
        );
    }
}

#[test]
fn vacuum_and_internal_beams_share_original_gross_species_allowance() {
    let gas = pure(2.5, 0.0, 400.0, 200.0);
    let gases = BTreeMap::from([(0, gas.clone()), (1, gas)]);
    let saved = gases.clone();
    let surface = interface(gases[&0].model(), support(), 1.0);
    let area = 0.06 * gases[&0].mass_kg() / (M * rate(1.0, 200.0, M, 400.0));
    let internal = [EffusionNetworkLink {
        id: 4,
        parcels: [0, 1],
        interface: &surface,
        area_m2: area,
    }];
    let paths = GasTransportPaths {
        effusion: &internal,
        ..Default::default()
    };
    let solves = endpoints(&gases);
    let vacuum = [outlet(&surface, area)];
    step_candidate(&gases, paths, &[], 1.0, &solves, controls(), || false).unwrap();
    step_candidate(
        &gases,
        GasTransportPaths::default(),
        &vacuum,
        1.0,
        &solves,
        controls(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        step_candidate(&gases, paths, &vacuum, 1.0, &solves, controls(), || false),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { id: 0, .. }
        ))
    ));
    assert_eq!(gases, saved);
}

#[test]
fn vacuum_rejects_final_source_support_escape_and_changed_caloric_identity() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 200.0))]);
    let saved = gases.clone();
    let solves = endpoints(&gases);
    for domain in [
        EffusiveSupport {
            temperature_k: (399.0, 1800.0),
            ..support()
        },
        EffusiveSupport {
            pressure_pa: (199.0, 1e7),
            ..support()
        },
    ] {
        let surface = interface(gases[&0].model(), domain, 1.0);
        // The original state/rate is valid; this finite cooling/depletion is not.
        gases[&0]
            .effusive_flux_to_vacuum(&surface, AREA, || false)
            .unwrap();
        assert!(
            step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[outlet(&surface, AREA)],
                1.0,
                &solves,
                controls(),
                || false
            )
            .is_err()
        );
    }
    let changed = pure(2.5, 1e12, 400.0, 200.0);
    let substituted = interface(changed.model(), support(), 1.0);
    assert!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[outlet(&substituted, AREA)],
            1.0,
            &solves,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(gases, saved);
}

#[test]
fn zero_area_bypasses_inactive_sources_and_zero_transmission_spends_nothing() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 200.0))]);
    let foreign = pure(3.5, 1e12, 800.0, 200.0);
    let inactive = interface(
        foreign.model(),
        EffusiveSupport {
            temperature_k: (700.0, 1000.0),
            ..support()
        },
        1.0,
    );
    let sealed = interface(gases[&0].model(), support(), 0.0);
    for (surface, area) in [(&inactive, 0.0), (&sealed, AREA)] {
        let result = step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[outlet(surface, area)],
            1.0,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.parcels(), &gases);
        assert!(result.step().exchange.boundaries.is_empty());
        assert_eq!(result.step().outlets[0].flux.mass_rate_kg_s, 0.0);
        assert_eq!(result.step().outlets[0].flux.energy_rate_w, 0.0);
    }
}

#[test]
fn vacuum_outlet_order_and_invalid_or_cancelled_proposals_are_atomic() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 200.0))]);
    let saved = gases.clone();
    let surface = interface(gases[&0].model(), support(), 1.0);
    let first = outlet(&surface, 0.3 * AREA);
    let second = VacuumEffusionLink {
        id: 2,
        area_m2: 0.2 * AREA,
        ..first
    };
    let solves = endpoints(&gases);
    let run = |links: &[VacuumEffusionLink<'_>]| {
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            links,
            1.0,
            &solves,
            controls(),
            || false,
        )
    };
    assert_eq!(
        run(&[first, second]).unwrap(),
        run(&[second, first]).unwrap()
    );
    assert!(run(&[first, first]).is_err());
    assert!(
        run(&[VacuumEffusionLink {
            parcel: 99,
            ..first
        }])
        .is_err()
    );
    assert!(
        run(&[VacuumEffusionLink {
            area_m2: -1.0,
            ..first
        }])
        .is_err()
    );
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[first],
                duration,
                &solves,
                controls(),
                || false
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[first, second],
            1.0,
            &solves,
            controls(),
            || {
                polls += 1;
                polls > 8
            }
        ),
        Err(GasTransportError::Cancelled)
    ));
    assert_eq!(gases, saved);
}
