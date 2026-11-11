//! Actual-mixture free paths from explicit, reciprocal binary collision data.
//!
//! For a dilute mixture of equilibrium Maxwellians at a common temperature,
//! constant total pair cross sections give
//! `nu_i = sum_j n_j sigma_ij <|v_i-v_j|>` and
//! `lambda_i = <|v_i|>/nu_i`.
//! Their Maxwellian speed ratio is `sqrt(1 + M_i/M_j)`, hence
//! `lambda_i = k_B T / (p sum_j x_j sigma_ij sqrt(1 + M_i/M_j))`.
//! A pure gas recovers `1/(sqrt(2) n sigma)`. All collision partners count,
//! including species that cannot pass the aperture. No cross section is inferred
//! from a species name, viscosity, diffusivity or Lennard-Jones diameter.
//!
//! This is an explicitly authored hard-sphere collision-frequency approximation,
//! not a rigorous free-flight probability of flux-selected molecules. Energy-
//! dependent scattering, non-Maxwellian counterstreaming and in-channel density
//! fields need other owners. The original global aperture and drift declarations
//! remain required; the optional actual-state check only narrows admission.

use super::{
    ClosedGasParcel, EffusiveSupport, GasEnergyError, MAX_EFFUSIVE_SPECIES, SpeciesId, Sum,
    ThermochemicalGasModel, admit_calibration, invalid, poll, positive, positive_product_ratio,
    valid_provenance,
};
use fs_math::det;
use std::sync::Arc;

const BOLTZMANN_J_K: f64 = 1.380_649e-23;

/// One explicitly declared symmetric, speed-independent TOTAL collision cross
/// section [m²]. It is constant throughout the model's authored T/p support.
/// Self pairs are required. Momentum-transfer and optical cross sections cannot
/// substitute for this total collision cross section without an explicit model.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveCollisionPair {
    species: [SpeciesId; 2],
    cross_section_m2: f64,
    provenance: String,
}
impl EffusiveCollisionPair {
    pub fn try_new(
        mut species: [SpeciesId; 2],
        cross_section_m2: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = provenance.into();
        if !positive(cross_section_m2) || !valid_provenance(&provenance) {
            return Err(invalid(
                "effusive collision pair needs a positive total cross section and bounded provenance",
            ));
        }
        species.sort();
        Ok(Self {
            species,
            cross_section_m2,
            provenance,
        })
    }
    pub const fn species(&self) -> &[SpeciesId; 2] {
        &self.species
    }
    pub const fn cross_section_m2(&self) -> f64 {
        self.cross_section_m2
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
}

/// Immutable original sources with a complete symmetric pair table.
///
/// Catalog fractions are not a privileged calibration mixture. Actual fractions,
/// pressure and temperature come exclusively from each accepted gas parcel.
/// Construction bounds quadratic work to 32 species and precomputes the pair
/// speed factors once. Evaluation uses fixed stack storage without allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveCollisionModel {
    gas_catalog: Arc<ThermochemicalGasModel>,
    pairs: Arc<[EffusiveCollisionPair]>,
    weighted_cross_sections_m2: Arc<[f64]>,
    support: EffusiveSupport,
    provenance: String,
}
impl EffusiveCollisionModel {
    pub fn try_new(
        gas_catalog: ThermochemicalGasModel,
        mut pairs: Vec<EffusiveCollisionPair>,
        support: EffusiveSupport,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        support.validate()?;
        let provenance = provenance.into();
        let components = gas_catalog.thermodynamics().components();
        let n = components.len();
        if n == 0 || n > MAX_EFFUSIVE_SPECIES || !valid_provenance(&provenance) {
            return Err(invalid(
                "effusive collision model needs a bounded original gas catalog and provenance",
            ));
        }
        if pairs.len() != n * (n + 1) / 2 {
            return Err(invalid(
                "effusive collision data must cover every unordered original pair including self pairs",
            ));
        }
        pairs.sort_by(|a, b| a.species.cmp(&b.species));
        let mut coefficients = vec![0.0; n * n];
        let mut row = 0;
        for (i, left) in components.iter().enumerate() {
            for (j, right) in components.iter().enumerate().skip(i) {
                let pair = &pairs[row];
                if pair.species[0] != *left.species() || pair.species[1] != *right.species() {
                    return Err(invalid(
                        "effusive collision table has a missing, repeated or foreign species pair",
                    ));
                }
                let left_mass = left.molar_mass().value();
                let right_mass = right.molar_mass().value();
                for (a, b, mass, other_mass) in
                    [(i, j, left_mass, right_mass), (j, i, right_mass, left_mass)]
                {
                    let speed_ratio = det::sqrt(1.0 + mass / other_mass);
                    coefficients[a * n + b] =
                        positive_product_ratio(&[pair.cross_section_m2, speed_ratio], &[])?;
                }
                row += 1;
            }
        }
        Ok(Self {
            gas_catalog: Arc::new(gas_catalog),
            pairs: pairs.into(),
            weighted_cross_sections_m2: coefficients.into(),
            support,
            provenance,
        })
    }

    pub fn gas_catalog(&self) -> &ThermochemicalGasModel {
        &self.gas_catalog
    }
    pub fn pairs(&self) -> &[EffusiveCollisionPair] {
        &self.pairs
    }
    pub const fn support(&self) -> EffusiveSupport {
        self.support
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Re-evaluate the actual original parcel, admitting source identities and
    /// the collision source's own support. Absent species have no reported path;
    /// every present species collides with every present partner regardless of
    /// its aperture transmission. No temperature or composition is reconstructed.
    pub fn evaluate_parcel(
        &self,
        parcel: &ClosedGasParcel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusiveCollisionEstimate, GasEnergyError> {
        poll(&mut cancelled)?;
        let state = parcel.state().state();
        self.support.admit(state.temperature, state.pressure)?;
        let mixture = parcel.model().thermodynamics();
        let components = mixture.components();
        if components.len() > MAX_EFFUSIVE_SPECIES {
            return Err(invalid(
                "effusive collision query exceeds the species work bound",
            ));
        }
        admit_calibration(parcel.model(), self.gas_catalog())?;
        let catalog = self.gas_catalog.thermodynamics().components();
        let n = catalog.len();
        let mut indices = [0; MAX_EFFUSIVE_SPECIES];
        for (i, component) in components.iter().enumerate() {
            poll(&mut cancelled)?;
            indices[i] = catalog
                .binary_search_by(|candidate| candidate.species().cmp(component.species()))
                .map_err(|_| invalid("effusive collision query lost an original species"))?;
        }
        let mut result = EffusiveCollisionEstimate {
            temperature_k: state.temperature,
            pressure_pa: state.pressure,
            species_count: n,
            mean_free_paths_m: [None; MAX_EFFUSIVE_SPECIES],
            minimum_mean_free_path_m: f64::INFINITY,
        };
        let fractions = mixture.mole_composition().fractions();
        for (i, &fraction) in fractions.iter().enumerate() {
            poll(&mut cancelled)?;
            if fraction == 0.0 {
                continue;
            }
            if !positive(fraction) {
                return Err(invalid(
                    "effusive collision query needs nonnegative finite mole fractions",
                ));
            }
            let row = indices[i] * n;
            let mut cross_section = Sum::default();
            for (j, &partner_fraction) in fractions.iter().enumerate() {
                if partner_fraction != 0.0 {
                    cross_section.add(positive_product_ratio(
                        &[
                            partner_fraction,
                            self.weighted_cross_sections_m2[row + indices[j]],
                        ],
                        &[],
                    )?)?;
                }
            }
            let path = positive_product_ratio(
                &[BOLTZMANN_J_K, state.temperature],
                &[state.pressure, cross_section.value],
            )?;
            result.mean_free_paths_m[indices[i]] = Some(path);
            result.minimum_mean_free_path_m = result.minimum_mean_free_path_m.min(path);
        }
        if !positive(result.minimum_mean_free_path_m) {
            return Err(invalid(
                "effusive collision query has no representable active species path",
            ));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

/// Mean-path estimates in the immutable collision catalog's canonical order.
/// An absent species has `None`; it is never assigned a pressure floor or a
/// fictitious trace population. These are model estimates, not error bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveCollisionEstimate {
    temperature_k: f64,
    pressure_pa: f64,
    species_count: usize,
    mean_free_paths_m: [Option<f64>; MAX_EFFUSIVE_SPECIES],
    minimum_mean_free_path_m: f64,
}
impl EffusiveCollisionEstimate {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }
    pub fn mean_free_paths_m(&self) -> &[Option<f64>] {
        &self.mean_free_paths_m[..self.species_count]
    }
    pub const fn minimum_mean_free_path_m(&self) -> f64 {
        self.minimum_mean_free_path_m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::{
        energy::{
            GasEnergySolve,
            parcel::transfer::species::{
                effusion::{
                    EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission,
                    drift::EffusiveDriftInterface, momentum::EffusiveMomentumInterface,
                },
                network::{
                    SpeciesNetworkSolve,
                    transport::{
                        GasTransportError, GasTransportPaths, GasTransportSolve,
                        effusion::{EffusionNetworkError, EffusionNetworkLink},
                        sample_transfers, step_candidate_with_paths,
                    },
                },
            },
        },
        tests::{model, point, species, transport},
    };
    use fs_qty::Amount;
    use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};
    use std::collections::BTreeMap;

    fn support() -> EffusiveSupport {
        EffusiveSupport {
            temperature_k: (200.0, 2000.0),
            pressure_pa: (1000.0, 1e7),
        }
    }
    fn pair(a: &str, b: &str, cross_section: f64) -> EffusiveCollisionPair {
        EffusiveCollisionPair::try_new(
            [SpeciesId::new(a).unwrap(), SpeciesId::new(b).unwrap()],
            cross_section,
            "synthetic constant total hard-sphere collision cross section",
        )
        .unwrap()
    }
    fn pairs() -> Vec<EffusiveCollisionPair> {
        vec![
            pair("A", "A", 4e-19),
            pair("A", "B", 3e-19),
            pair("B", "B", 2e-19),
        ]
    }
    fn collision(bank: ThermochemicalGasModel) -> EffusiveCollisionModel {
        EffusiveCollisionModel::try_new(
            bank,
            pairs(),
            support(),
            "synthetic Maxwellian hard-sphere collision model",
        )
        .unwrap()
    }
    fn parcel(fractions: [f64; 2], temperature: f64, pressure: f64) -> ClosedGasParcel {
        let bank = model(vec![0.5, 0.5], 0.0)
            .with_mole_amounts(&fractions.map(Amount::new))
            .unwrap();
        ClosedGasParcel::new(bank, point(temperature, pressure), 0.001).unwrap()
    }
    fn interface(bank: ThermochemicalGasModel, transmissions: [f64; 2]) -> EffusiveInterface {
        let rows = bank
            .thermodynamics()
            .components()
            .iter()
            .zip(transmissions)
            .map(|(s, probability)| {
                EffusiveSpeciesTransmission::try_new(
                    s.species().clone(),
                    probability,
                    "synthetic reciprocal velocity-independent aperture transmission",
                )
                .unwrap()
            })
            .collect();
        EffusiveInterface::try_new(
            bank,
            rows,
            support(),
            EffusiveRegime::try_new(
                1e-9,
                0.0,
                1e-7,
                "original authored thin-aperture qualification; not measured accuracy",
            )
            .unwrap(),
        )
        .unwrap()
    }
    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-12 * expected.abs(),
            "{actual:.17e} != {expected:.17e}"
        );
    }

    #[test]
    fn g1_pure_limit_and_temperature_pressure_scaling_use_actual_gas() {
        let law = collision(model(vec![0.5, 0.5], 0.0));
        for (t, p) in [(300.0, 1e5), (600.0, 1e5), (300.0, 2e5)] {
            let gas = parcel([1.0, 0.0], t, p);
            let estimate = law.evaluate_parcel(&gas, || false).unwrap();
            let expected = BOLTZMANN_J_K * t / (p * 4e-19 * 2.0_f64.sqrt());
            assert_eq!(estimate.temperature_k(), t);
            assert_eq!(estimate.pressure_pa(), p);
            near(estimate.mean_free_paths_m()[0].unwrap(), expected);
            near(estimate.minimum_mean_free_path_m(), expected);
            assert_eq!(estimate.mean_free_paths_m()[1], None);
        }
    }

    #[test]
    fn g1_unequal_mass_mixture_uses_all_original_binary_partners() {
        let law = collision(model(vec![0.5, 0.5], 0.0));
        let gas = parcel([0.25, 0.75], 400.0, 1e5);
        let estimate = law.evaluate_parcel(&gas, || false).unwrap();
        let density = 1e5 / (BOLTZMANN_J_K * 400.0);
        let expected_a = 1.0
            / (density
                * (0.25 * 4e-19 * 2.0_f64.sqrt()
                    + 0.75 * 3e-19 * (1.0_f64 + 0.028 / 0.004).sqrt()));
        let expected_b = 1.0
            / (density
                * (0.25 * 3e-19 * (1.0_f64 + 0.004 / 0.028).sqrt()
                    + 0.75 * 2e-19 * 2.0_f64.sqrt()));
        near(estimate.mean_free_paths_m()[0].unwrap(), expected_a);
        near(estimate.mean_free_paths_m()[1].unwrap(), expected_b);
        near(
            estimate.minimum_mean_free_path_m(),
            expected_a.min(expected_b),
        );

        let reordered = EffusiveCollisionModel::try_new(
            model(vec![0.5, 0.5], 0.0),
            vec![
                pair("B", "B", 2e-19),
                pair("B", "A", 3e-19),
                pair("A", "A", 4e-19),
            ],
            support(),
            "synthetic Maxwellian hard-sphere collision model",
        )
        .unwrap();
        assert_eq!(law, reordered);
        assert_eq!(estimate, reordered.evaluate_parcel(&gas, || false).unwrap());
    }

    #[test]
    fn g0_complete_pair_identity_source_and_domain_are_required() {
        let bank = model(vec![0.5, 0.5], 0.0);
        for rows in [
            vec![pair("A", "A", 4e-19), pair("B", "B", 2e-19)],
            vec![
                pair("A", "A", 4e-19),
                pair("B", "B", 2e-19),
                pair("B", "B", 3e-19),
            ],
            vec![
                pair("A", "A", 4e-19),
                pair("A", "C", 3e-19),
                pair("B", "B", 2e-19),
            ],
        ] {
            assert!(
                EffusiveCollisionModel::try_new(bank.clone(), rows, support(), "test").is_err()
            );
        }
        for cross_section in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                EffusiveCollisionPair::try_new(
                    [SpeciesId::new("A").unwrap(), SpeciesId::new("A").unwrap()],
                    cross_section,
                    "test",
                )
                .is_err()
            );
        }
        let changed_source = collision(model(vec![0.5, 0.5], 0.001));
        assert!(
            interface(bank.clone(), [1.0, 1.0])
                .with_collision_model(changed_source)
                .is_err()
        );

        let narrow = EffusiveCollisionModel::try_new(
            bank.clone(),
            pairs(),
            EffusiveSupport {
                temperature_k: (300.0, 500.0),
                pressure_pa: (1e5, 2e5),
            },
            "explicit narrow synthetic collision support",
        )
        .unwrap();
        let law = interface(bank, [1.0, 1.0])
            .with_collision_model(narrow)
            .unwrap();
        assert!(
            law.admit_parcel(&parcel([0.5, 0.5], 400.0, 1e5), || false)
                .is_ok()
        );
        // Both are native gas states with the same p/T ratio and thus mean
        // path. Collision-source support is still independently necessary.
        assert!(
            law.admit_parcel(&parcel([0.5, 0.5], 800.0, 2e5), || false)
                .is_err()
        );
    }

    #[test]
    fn g3_blocked_species_still_collide_and_active_path_length_is_not_ignored() {
        let bank = model(vec![0.5, 0.5], 0.0);
        let collisions = EffusiveCollisionModel::try_new(
            bank.clone(),
            vec![
                pair("A", "A", 1e-20),
                pair("A", "B", 9e-18),
                pair("B", "B", 1e-20),
            ],
            support(),
            "synthetic large unlike-species collision area",
        )
        .unwrap();
        let original = interface(bank, [1.0, 0.0]);
        let law = original.clone().with_collision_model(collisions).unwrap();
        let pure_a = parcel([1.0, 0.0], 400.0, 1e5);
        let mixed = parcel([0.5, 0.5], 400.0, 1e5);
        law.admit_parcel(&pure_a, || false).unwrap();
        assert!(
            pure_a
                .effusive_flux_to_vacuum(&law, 1e-12, || false)
                .is_ok()
        );
        assert!(
            mixed
                .effusive_flux_to_vacuum(&original, 1e-12, || false)
                .is_ok()
        );
        assert!(
            mixed
                .effusive_flux_to_vacuum(&law, 1e-12, || false)
                .is_err()
        );
        assert!(
            pure_a
                .effusive_flux_to(&mixed, &law, 1e-12, || false)
                .is_err()
        );
        assert!(
            mixed
                .effusive_flux_to(&pure_a, &law, 1e-12, || false)
                .is_err()
        );
        // The same immutable model checks an independently authored long path.
        law.admit_collision_length(&pure_a, 1e-8, || false).unwrap();
        assert!(law.admit_collision_length(&pure_a, 1e-4, || false).is_err());
        assert!(
            law.admit_collision_length(&pure_a, f64::NAN, || false)
                .is_err()
        );
        assert_eq!(law.regime(), original.regime());
        assert_eq!(law.support(), original.support());
    }

    #[test]
    fn g3_zero_area_is_inactive_but_zero_transmission_and_drift_keep_the_guard() {
        let bank = model(vec![0.5, 0.5], 0.0);
        let law = interface(bank.clone(), [0.0, 0.0])
            .with_collision_model(collision(bank))
            .unwrap();
        let dense = parcel([0.5, 0.5], 400.0, 1e7);
        assert!(dense.effusive_flux_to_vacuum(&law, 0.0, || false).is_ok());
        assert!(dense.effusive_flux_to(&dense, &law, 0.0, || false).is_ok());
        assert!(
            dense
                .effusive_flux_to_vacuum(&law, 1e-12, || false)
                .is_err()
        );
        assert!(
            dense
                .effusive_flux_to(&dense, &law, 1e-12, || false)
                .is_err()
        );
        let momentum = EffusiveMomentumInterface::try_new(
            law,
            "synthetic nondeflecting angle-independent transmission",
        )
        .unwrap();
        let drift = EffusiveDriftInterface::try_new(
            momentum,
            100.0,
            "synthetic collisionless drift support; collision guard still applies",
        )
        .unwrap();
        assert!(
            dense
                .drifting_effusive_flux(&drift, 0.0, [10.0, 0.0, 0.0], [1.0, 0.0, 0.0], || false)
                .is_ok()
        );
        assert!(
            dense
                .drifting_effusive_flux(&drift, 1e-12, [10.0, 0.0, 0.0], [1.0, 0.0, 0.0], || false)
                .is_err()
        );
    }

    #[test]
    fn g0_momentum_qualification_cannot_borrow_a_rebound_collision_source() {
        let bank = model(vec![0.5, 0.5], 0.0);
        let original = interface(bank.clone(), [1.0, 1.0]);
        let bound = original
            .clone()
            .with_collision_model(collision(bank.clone()))
            .unwrap();
        let replacement = original
            .clone()
            .with_collision_model(collision(bank))
            .unwrap();
        let momentum =
            EffusiveMomentumInterface::try_new(bound.clone(), "original angular source").unwrap();
        let gas = parcel([0.5, 0.5], 400.0, 1e5);
        let admitted = gas
            .effusive_flux_to_vacuum(&bound, 1e-12, || false)
            .unwrap();
        momentum.vacuum_flux(&admitted, || false).unwrap();
        // Even equal model values are a new source declaration, just as a new
        // interface cannot reuse another interface's angular qualification.
        for other in [original, replacement] {
            let flux = gas
                .effusive_flux_to_vacuum(&other, 1e-12, || false)
                .unwrap();
            assert!(momentum.vacuum_flux(&flux, || false).is_err());
        }
    }

    #[test]
    fn g4_final_native_network_states_retain_the_collision_guard_atomically() {
        let bank = model(vec![0.5, 0.5], 0.0);
        let law = interface(bank.clone(), [1.0, 1.0])
            .with_collision_model(collision(bank))
            .unwrap();
        let gas = parcel([0.5, 0.5], 400.0, 1e5);
        let original = BTreeMap::from([(0, gas.clone()), (1, gas)]);
        let links = [EffusionNetworkLink {
            id: 7,
            parcels: [0, 1],
            interface: &law,
            area_m2: 1e-12,
        }];
        let paths = GasTransportPaths {
            effusion: &links,
            ..Default::default()
        };
        let energy = GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        };
        let controls = GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        };
        let endpoints = BTreeMap::from([(0, energy), (1, energy)]);
        let candidate =
            step_candidate_with_paths(&original, paths, 1e-6, &endpoints, controls, || false)
                .unwrap();
        let sampled =
            sample_transfers(&original, paths, 1e-6, &endpoints, controls, || false).unwrap();
        let mut compressed = original.clone();
        compressed.insert(1, parcel([0.5, 0.5], 400.0, 1e7));
        for result in [
            candidate
                .step()
                .admit_endpoint_driving_forces(&compressed, || false),
            sampled
                .rates
                .admit_endpoint_driving_forces(&compressed, || false),
        ] {
            assert!(matches!(
                result,
                Err(GasTransportError::Effusion(EffusionNetworkError::Face {
                    id: 7,
                    ..
                }))
            ));
        }
        assert_eq!(candidate.parcels(), &original);
        let before = original.clone();
        let mut calls = 0;
        assert!(matches!(
            law.admit_parcel(&original[&0], || {
                calls += 1;
                calls >= 6
            }),
            Err(GasEnergyError::Cancelled)
        ));
        assert_eq!(original, before);
    }

    #[test]
    fn g0_equal_mass_mixtures_recover_the_pure_limit_at_the_work_bound() {
        let make_bank = |count: usize| {
            let names: Vec<_> = (0..count).map(|i| format!("S{i:02}")).collect();
            // fs-thermochem mixture.rs:349 sorts the fractions into canonical SpeciesId
            // order and :366-367 refuse unless their `Iterator::sum` is bit-exact 1.0
            // (d571a4691d); 33 copies of fl(1/33) sum to 1 - 6 * 2^-53. S00, S01, ...
            // already sort in that order, so the last fraction closes the sum as 1 - s,
            // with s the sum of the others taken the same way. With the partial sum s in
            // [0.5, 2], Sterbenz's lemma makes 1 - s exact, so s + (1 - s) == 1 exactly.
            // Here s = 31/32 (32 names, every fraction stays 1/32) or 0.9696969696969691
            // (33 names).
            let mut fractions = vec![1.0 / count as f64; count - 1];
            let others: f64 = fractions.iter().sum();
            fractions.push(1.0 - others);
            assert_eq!(fractions.iter().sum::<f64>().to_bits(), 1.0_f64.to_bits());
            ThermochemicalGasModel::new(
                FrozenIdealGasMixtureModelV1::new(
                    names
                        .iter()
                        .map(|name| species(name, 0.028, 3.5, 0.0))
                        .collect(),
                    Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
                )
                .unwrap(),
                names
                    .iter()
                    .map(|name| transport(name, 1.7e-5, 0.025))
                    .collect(),
                None,
            )
            .unwrap()
        };
        let bank = make_bank(MAX_EFFUSIVE_SPECIES);
        let mut rows = Vec::new();
        for (i, left) in bank.thermodynamics().components().iter().enumerate() {
            for right in bank.thermodynamics().components().iter().skip(i) {
                rows.push(pair(
                    left.species().as_str(),
                    right.species().as_str(),
                    4e-19,
                ));
            }
        }
        let law = EffusiveCollisionModel::try_new(
            bank.clone(),
            rows,
            support(),
            "equal-mass synthetic gas",
        )
        .unwrap();
        let gas = ClosedGasParcel::new(bank, point(400.0, 1e5), 0.001).unwrap();
        let estimate = law.evaluate_parcel(&gas, || false).unwrap();
        assert_eq!(estimate.mean_free_paths_m().len(), MAX_EFFUSIVE_SPECIES);
        for &path in estimate.mean_free_paths_m() {
            near(
                path.unwrap(),
                BOLTZMANN_J_K * 400.0 / (1e5 * 4e-19 * 2.0_f64.sqrt()),
            );
        }
        assert!(
            EffusiveCollisionModel::try_new(
                make_bank(MAX_EFFUSIVE_SPECIES + 1),
                Vec::new(),
                support(),
                "over-bound synthetic gas",
            )
            .is_err()
        );
    }
}
