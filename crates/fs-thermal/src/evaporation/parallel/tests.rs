//! Independent parallel fluxes and conserved single-inventory thermal feedback.
use super::*;
use crate::evaporation::{
    KineticEvaporation, LiquidSurfaceState, MOLAR_GAS_CONSTANT, SolverConfig,
    advance_surface_with_heat, surface_flux,
};
use std::cell::Cell;

struct Liquid {
    property_queries: Cell<usize>,
    film_queries: Cell<usize>,
    species: [String; 1],
    fractions: [f64; 1],
}
impl Liquid {
    fn new() -> Self {
        Self {
            property_queries: Cell::new(0),
            film_queries: Cell::new(0),
            species: ["original-liquid".into()],
            fractions: [1.0],
        }
    }
}
fn phase(t: f64) -> EvaporationProperties {
    EvaporationProperties {
        saturation_pressure_pa: 4000.0 + 100.0 * (t - 300.0),
        saturation_pressure_derivative_pa_k: 100.0,
        latent_heat_j_kg: 2.4e6 - 1200.0 * (t - 300.0),
        liquid_density_kg_m3: 1000.0 - 0.2 * (t - 300.0),
        molar_mass_kg_mol: 0.018,
    }
}
impl EvaporationLaw for Liquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.property_queries.set(self.property_queries.get() + 1);
        Ok(phase(t))
    }
    fn evaluate_film(
        &self,
        p: EvaporationProperties,
        t: f64,
        film: GasFilm,
        kinetics: InterfaceKinetics,
    ) -> Result<SurfaceFlux, EvaporationError> {
        self.film_queries.set(self.film_queries.get() + 1);
        crate::evaporation::flux_with_kinetics(p, t, film, kinetics)
    }
}
impl ThermalLaw for Liquid {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([41; 32])
    }
    fn validate_cell(&self, _: &DepositedCell) -> Result<(), DepositionError> {
        Err(DepositionError::MaterialDomain)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        Ok(300.0 + h / 4000.0)
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        Ok(Some(phase(self.temperature_k(h)?).liquid_density_kg_m3))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 4000.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some(
            4000.0 * ((300.0 + b / 4000.0) / (300.0 + a / 4000.0)).ln(),
        ))
    }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        Some(CaloricComposition {
            species: &self.species,
            mass_fractions: &self.fractions,
            enthalpy_reference: "original-liquid-enthalpy",
        })
    }
    fn advected_entropy_anchor(&self) -> Option<CaloricEntropyAnchor<'_>> {
        Some(CaloricEntropyAnchor {
            specific_h_j_kg: 123.0,
            reduced_specific_entropy_j_kg_k: 45.0,
            entropy_reference: "original-liquid-entropy",
        })
    }
}

fn patches() -> [EvaporationPatch; 3] {
    [
        EvaporationPatch {
            area_m2: 0.002,
            mass_transfer_coefficient_m_s: 1e-4,
        },
        EvaporationPatch {
            area_m2: 0.008,
            mass_transfer_coefficient_m_s: 0.2,
        },
        // This coefficient would overflow the film evaluation. No exposed area
        // means no transport or source evaluation, but the coefficient is valid.
        EvaporationPatch {
            area_m2: 0.0,
            mass_transfer_coefficient_m_s: f64::MAX,
        },
    ]
}
fn film(vapor: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa: vapor,
        temperature_k: 310.0,
        // The declared parallel patches replace this otherwise valid value.
        mass_transfer_coefficient_m_s: 0.0,
    }
}
fn hk(alpha: f64) -> InterfaceKinetics {
    InterfaceKinetics::HertzKnudsen {
        accommodation_coefficient: alpha,
    }
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

// Independent normal-pressure bisection of the two physical flux equations.
// This does not call the production flux, scaled kinetic solver or wrapper.
fn reference(t: f64, bath: GasFilm, coefficient: f64, kinetics: InterfaceKinetics) -> f64 {
    let p = phase(t);
    let c = coefficient * p.molar_mass_kg_mol * bath.total_pressure_pa
        / (MOLAR_GAS_CONSTANT * bath.temperature_k);
    let factor = match kinetics {
        InterfaceKinetics::Equilibrium => {
            return c
                * ((bath.total_pressure_pa - bath.vapor_partial_pressure_pa)
                    / (bath.total_pressure_pa - p.saturation_pressure_pa))
                    .ln();
        }
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient,
        } => accommodation_coefficient,
        InterfaceKinetics::Schrage {
            accommodation_coefficient: alpha,
        } => 2.0 * alpha / (2.0 - alpha),
    };
    if factor == 0.0 || coefficient == 0.0 {
        return 0.0;
    }
    let k = factor
        * (p.molar_mass_kg_mol / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * t)).sqrt();
    let mut low = p.saturation_pressure_pa.min(bath.vapor_partial_pressure_pa);
    let mut high = p.saturation_pressure_pa.max(bath.vapor_partial_pressure_pa);
    for _ in 0..100 {
        let interface = 0.5 * (low + high);
        let residual = c
            * ((bath.total_pressure_pa - bath.vapor_partial_pressure_pa)
                / (bath.total_pressure_pa - interface))
                .ln()
            - k * (p.saturation_pressure_pa - interface);
        if residual < 0.0 {
            low = interface;
        } else {
            high = interface;
        }
    }
    k * (p.saturation_pressure_pa - 0.5 * (low + high))
}
fn reference_rate(t: f64, bath: GasFilm, kinetics: InterfaceKinetics) -> f64 {
    patches()
        .iter()
        .filter(|patch| patch.area_m2 > 0.0)
        .map(|patch| {
            patch.area_m2 * reference(t, bath, patch.mass_transfer_coefficient_m_s, kinetics)
        })
        .sum()
}

#[test]
fn heterogeneous_films_sum_exact_physics_in_both_directions_at_the_actual_temperature() {
    let source = Liquid::new();
    let patches = patches();
    for kinetics in [
        InterfaceKinetics::Equilibrium,
        hk(1e-5),
        InterfaceKinetics::Schrage {
            accommodation_coefficient: 1e-5,
        },
    ] {
        let selected = KineticEvaporation::try_new(&source, kinetics).unwrap();
        let parallel = ParallelFilmEvaporation::try_new(&selected, &patches).unwrap();
        for vapor in [500.0, 10_000.0] {
            for t in [290.0, 300.0, 315.0] {
                source.property_queries.set(0);
                source.film_queries.set(0);
                let actual = surface_flux(&parallel, t, film(vapor)).unwrap();
                close(
                    actual.mass_flux_kg_m2_s * parallel.area_m2(),
                    reference_rate(t, film(vapor), kinetics),
                    1e-17,
                    2e-11,
                );
                assert_eq!(
                    source.property_queries.get(),
                    1,
                    "one original property query per state"
                );
                assert_eq!(
                    source.film_queries.get(),
                    2,
                    "zero-area patch is not evaluated"
                );
                assert_eq!(
                    actual.latent_heat_flux_w_m2,
                    actual.mass_flux_kg_m2_s * phase(t).latent_heat_j_kg
                );
                assert_eq!(
                    actual.outward_normal_velocity_m_s,
                    -actual.mass_flux_kg_m2_s / phase(t).liquid_density_kg_m3
                );
            }
        }
        if kinetics != InterfaceKinetics::Equilibrium {
            let averaged = reference(300.0, film(500.0), 0.16002, kinetics);
            let actual = surface_flux(&parallel, 300.0, film(500.0))
                .unwrap()
                .mass_flux_kg_m2_s;
            assert!(
                (averaged - actual).abs() > 0.05 * actual.abs(),
                "averaging coefficients must be distinguishable from parallel nonlinear films"
            );
        }
    }
}

#[test]
fn kinetic_and_parallel_wrappers_compose_in_either_order_and_keep_the_tangent() {
    let source = Liquid::new();
    let patches = patches();
    let selected = KineticEvaporation::try_new(&source, hk(1e-5)).unwrap();
    let parallel_after = ParallelFilmEvaporation::try_new(&selected, &patches).unwrap();
    let parallel_before = ParallelFilmEvaporation::try_new(&source, &patches).unwrap();
    let kinetic_after = KineticEvaporation::try_new(&parallel_before, hk(1e-5)).unwrap();
    let law: &dyn EvaporationLaw = &kinetic_after;
    for vapor in [500.0, 10_000.0] {
        let bath = film(vapor);
        let actual = surface_flux(law, 300.0, bath).unwrap();
        assert_eq!(actual, surface_flux(&parallel_after, 300.0, bath).unwrap());
        let delta = 0.002;
        let derivative = (reference_rate(300.0 + delta, bath, hk(1e-5))
            - reference_rate(300.0 - delta, bath, hk(1e-5)))
            / (2.0 * delta * parallel_after.area_m2());
        close(
            actual.mass_flux_derivative_kg_m2_s_k,
            derivative,
            1e-14,
            2e-5,
        );
    }
    assert!(std::ptr::eq(parallel_before.source(), &source));
    assert!(std::ptr::eq(parallel_before.patches(), &patches[..]));
}

#[test]
fn uniform_or_single_exposed_patch_preserves_the_original_flux_bitwise() {
    let source = Liquid::new();
    let selected = KineticEvaporation::try_new(&source, hk(1e-5)).unwrap();
    for patches in [
        [
            EvaporationPatch {
                area_m2: 2.0,
                mass_transfer_coefficient_m_s: 0.003,
            },
            EvaporationPatch {
                area_m2: 0.0,
                mass_transfer_coefficient_m_s: 1.0,
            },
        ],
        [
            EvaporationPatch {
                area_m2: 2.0,
                mass_transfer_coefficient_m_s: 0.003,
            },
            EvaporationPatch {
                area_m2: 3.0,
                mass_transfer_coefficient_m_s: 0.003,
            },
        ],
    ] {
        let parallel = ParallelFilmEvaporation::try_new(&selected, &patches).unwrap();
        source.film_queries.set(0);
        let actual = surface_flux(&parallel, 300.0, film(500.0)).unwrap();
        assert_eq!(source.film_queries.get(), 1);
        let expected = surface_flux(
            &selected,
            300.0,
            GasFilm {
                mass_transfer_coefficient_m_s: 0.003,
                ..film(500.0)
            },
        )
        .unwrap();
        assert_eq!(actual, expected);
    }
}

#[test]
fn patch_validation_seals_and_source_failures_are_not_bypassed() {
    let source = Liquid::new();
    assert!(ParallelFilmEvaporation::try_new(&source, &[]).is_err());
    assert!(
        ParallelFilmEvaporation::try_new(
            &source,
            &[EvaporationPatch {
                area_m2: 0.0,
                mass_transfer_coefficient_m_s: 0.01,
            }]
        )
        .is_err()
    );
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ParallelFilmEvaporation::try_new(
                &source,
                &[EvaporationPatch {
                    area_m2: invalid,
                    mass_transfer_coefficient_m_s: 0.01,
                }]
            )
            .is_err()
        );
        assert!(
            ParallelFilmEvaporation::try_new(
                &source,
                &[EvaporationPatch {
                    area_m2: 1.0,
                    mass_transfer_coefficient_m_s: invalid,
                }]
            )
            .is_err()
        );
    }
    assert!(
        ParallelFilmEvaporation::try_new(
            &source,
            &[EvaporationPatch {
                area_m2: f64::MAX,
                mass_transfer_coefficient_m_s: 0.01,
            }; 2]
        )
        .is_err()
    );
    let sealed = KineticEvaporation::try_new(&source, hk(0.0)).unwrap();
    let patches = patches();
    let parallel = ParallelFilmEvaporation::try_new(&sealed, &patches).unwrap();
    let flux = surface_flux(&parallel, 300.0, film(500.0)).unwrap();
    assert_eq!(flux.mass_flux_kg_m2_s, 0.0);
    assert_eq!(flux.mass_flux_derivative_kg_m2_s_k, 0.0);
    assert!(matches!(
        surface_flux(&parallel, 279.0, film(500.0)),
        Err(EvaporationError::TemperatureOutOfRange { .. })
    ));
    assert!(matches!(
        surface_flux(
            &parallel,
            300.0,
            GasFilm {
                total_pressure_pa: 4000.0,
                ..film(500.0)
            }
        ),
        Err(EvaporationError::Boiling { .. })
    ));
    assert!(surface_flux(&parallel, 300.0, film(100_000.0)).is_err());
}

#[test]
fn area_and_flux_accumulation_retain_small_contributions_without_intermediate_overflow() {
    let source = Liquid::new();
    let areas = [
        EvaporationPatch {
            area_m2: 1e16,
            mass_transfer_coefficient_m_s: 1.0,
        },
        EvaporationPatch {
            area_m2: 1.0,
            mass_transfer_coefficient_m_s: 1.0,
        },
        EvaporationPatch {
            area_m2: 1.0,
            mass_transfer_coefficient_m_s: 1.0,
        },
    ];
    let parallel = ParallelFilmEvaporation::try_new(&source, &areas).unwrap();
    assert_eq!(parallel.area_m2(), 1e16 + 2.0);
    assert_eq!(weighted(1e300, 1e200, 2e200).unwrap(), 5e299);
    close(weighted(1e300, 1e-200, 1e300).unwrap(), 1e-200, 0.0, 2e-15);
    assert!(weighted(f64::from_bits(1), 1.0, 1e300).is_err());
}

#[test]
fn scalar_heat_and_phase_close_one_inventory_using_the_parallel_endpoint_flux() {
    let source = Liquid::new();
    let patches = patches();
    let selected = KineticEvaporation::try_new(&source, hk(1e-5)).unwrap();
    let parallel = ParallelFilmEvaporation::try_new(&selected, &patches).unwrap();
    for vapor in [500.0, 10_000.0] {
        let initial = LiquidSurfaceState {
            mass_kg: 0.1,
            specific_enthalpy_j_kg: 0.0,
            area_m2: parallel.area_m2(),
        };
        let duration = 100.0;
        let target = 305.0;
        let expected_dm = reference_rate(target, film(vapor), hk(1e-5)) * duration;
        let supplied = initial.mass_kg * 4000.0 * (target - 300.0)
            + expected_dm * phase(target).latent_heat_j_kg;
        let next = advance_surface_with_heat(
            &parallel,
            film(vapor),
            initial,
            duration,
            supplied,
            SolverConfig {
                specific_enthalpy_bracket_j_kg: [-80_000.0, 120_000.0],
                relative_tolerance: 1e-11,
                absolute_energy_tolerance_j: 1e-9,
                max_iterations: 128,
            },
            |h| Ok(300.0 + h / 4000.0),
            || Ok(()),
        )
        .unwrap();
        close(next.surface_temperature_k, target, 2e-8, 0.0);
        close(next.transferred_mass_kg, expected_dm, 1e-13, 1e-9);
        assert_eq!(
            next.transferred_mass_kg,
            initial.mass_kg - next.liquid.mass_kg
        );
        assert!(next.energy_residual_j.abs() < 1e-7);
        close(
            next.liquid.mass_kg * next.liquid.specific_enthalpy_j_kg + next.vapor_energy_j,
            supplied,
            1e-7,
            0.0,
        );
        assert_eq!(initial.mass_kg, 0.1);
    }
}

#[test]
fn original_caloric_identity_density_entropy_and_species_capabilities_survive() {
    use crate::enthalpy::{Boundary, CellGeometry, CellId, MaterialId};
    let source = Liquid::new();
    let patches = patches();
    let parallel = ParallelFilmEvaporation::try_new(&source, &patches).unwrap();
    assert_eq!(parallel.storage_basis(), source.storage_basis());
    assert_eq!(parallel.caloric_identity(), source.caloric_identity());
    assert_eq!(parallel.temperature_k(4000.0), source.temperature_k(4000.0));
    assert_eq!(
        parallel.maximum_temperature_slope(),
        source.maximum_temperature_slope()
    );
    assert_eq!(
        parallel.bulk_density_kg_m3(4000.0),
        source.bulk_density_kg_m3(4000.0)
    );
    assert_eq!(
        parallel.specific_entropy_change_j_kg_k(0.0, 4000.0),
        source.specific_entropy_change_j_kg_k(0.0, 4000.0)
    );
    let composition = parallel.advected_composition().unwrap();
    assert!(std::ptr::eq(composition.species, &source.species[..]));
    assert!(std::ptr::eq(
        composition.mass_fractions,
        &source.fractions[..]
    ));
    assert_eq!(composition.enthalpy_reference, "original-liquid-enthalpy");
    let anchor = parallel.advected_entropy_anchor().unwrap();
    assert_eq!(anchor.specific_h_j_kg, 123.0);
    assert_eq!(anchor.reduced_specific_entropy_j_kg_k, 45.0);
    assert_eq!(anchor.entropy_reference, "original-liquid-entropy");
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 0.1,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.1, 0.0, 0.0],
            volume_m3: 1e-4,
            surface_area_m2: 0.01,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    assert_eq!(parallel.validate_cell(&cell), source.validate_cell(&cell));
}
