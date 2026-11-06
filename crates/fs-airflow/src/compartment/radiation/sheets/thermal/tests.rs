use super::*;
use crate::compartment::radiation::spectral::{BandSurfaceEmissivity, SpectralBand};
use crate::compartment::radiation::{AbsorptionCoefficient, OpticalPath};
use fs_conduction::radiation::{
    BlackbodyBand, EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY,
    ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance, sheet::ReciprocalSheetOptics,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Length};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent,
};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / 1000.0;
        if (250.0..=2200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}
static SENSIBLE: Sensible = Sensible;

struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=1e6).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 10_000.0 {
            300.0 + h / 1000.0
        } else if h <= 110_000.0 {
            310.0
        } else {
            310.0 + (h - 110_000.0) / 1000.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}
static LATENT: Latent = Latent;

fn state(
    law: &'static dyn ThermalLaw,
    inventory: &[(f64, f64)],
    conductance: Option<f64>,
    last_source: f64,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = inventory
        .iter()
        .enumerate()
        .map(|(i, &(mass, h))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: i as u64 + 1,
                body_id: i as u64 + 1,
                material_id: MaterialId(1),
                mass_kg: mass,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: mass / 1000.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: if i + 1 == inventory.len() {
                        last_source
                    } else {
                        0.0
                    },
                },
            }),
        })
        .collect();
    if let Some(g) = conductance {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(1),
                b: CellId(2),
                area_m2: 1.0,
                conductance_w_k: g,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn enclosure() -> ParticipatingEnclosure {
    let geometry = ViewFactorMatrix::admit(
        vec![1.0; 4],
        vec![
            vec![0.0, 1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
            vec![0.0, 0.0, 1.0, 0.0],
        ],
        ViewFactorEvidence::Analytic {
            geometry: "two black infinite plates separated by a two-sided sheet".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    ParticipatingEnclosure::new(
        geometry,
        [[0, 1], [2, 3]]
            .map(|surfaces| OpticalPath {
                surfaces,
                exchange_area: Area::new(1.0),
                segments: vec![],
            })
            .to_vec(),
        "analytic reciprocal vacuum gaps",
        1e-12,
        || false,
    )
    .unwrap()
}

fn sheet(a: [f64; 2], transmission: f64) -> SheetPair {
    SheetPair {
        surfaces: [1, 2],
        optics: ReciprocalSheetOptics::new(
            a,
            transmission,
            [250.0, 2200.0],
            "declared analytic sheet coefficients",
        )
        .unwrap(),
    }
}
fn bindings(baths: [f64; 2]) -> [SheetThermalSurface; 4] {
    [
        SheetThermalSurface::Prescribed(Temperature::new(baths[0])),
        SheetThermalSurface::Stored(CellId(1)),
        SheetThermalSurface::Stored(CellId(1)),
        SheetThermalSurface::Prescribed(Temperature::new(baths[1])),
    ]
}
fn absorption() -> [LayerAbsorption; 2] {
    [0.0; 2].map(|v| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(v),
            "declared vacuum optical fixture",
        )
        .unwrap()
    })
}
fn black(surface: usize, temperature: Temperature) -> Result<SurfaceEmissivity, ConductionError> {
    assert!(
        surface == 0 || surface == 3,
        "sheet faces must not query opaque gray cards"
    );
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 1.0,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 2200.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored black boundary".into(),
                license: "internal-test".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic black bath".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card(
        "bath",
        &card,
        temperature.value(),
        SelectionPolicy::SingleClaimOnly,
    )
}
fn gray(
    thermal: &DepositionSimulation<'static>,
    sheets: &[SheetPair],
    baths: [f64; 2],
    dt: f64,
    controls: SheetThermalSolve,
) -> SheetThermalCandidate<'static, 'static> {
    let enclosure = enclosure();
    enclosure
        .candidate_sheet_thermal(
            thermal,
            Time::new(dt),
            &bindings(baths),
            [Temperature::new(300.0); 2],
            SheetThermalOptics::Gray {
                absorption: &absorption(),
                sheets,
                solve: RadiationSolve::default(),
            },
            black,
            controls,
            &mut enclosure.sheet_thermal_workspace(thermal).unwrap(),
            || false,
        )
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, allowed {tolerance}"
    );
}
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
fn gray_heat(t: f64, baths: [f64; 2], a: [f64; 2]) -> f64 {
    STEFAN_BOLTZMANN_W_M2_K4
        * (a[0] * (baths[0].powi(4) - t.powi(4)) + a[1] * (baths[1].powi(4) - t.powi(4)))
}

#[test]
fn actual_sheet_inventory_matches_independent_stiff_backward_euler_balance() {
    let initial = state(&SENSIBLE, &[(0.05, 0.0)], None, 0.0);
    let baths = [900.0, 300.0];
    let a = [0.3, 0.1];
    let dt = 10.0;
    let result = gray(
        &initial,
        &[sheet(a, 0.5)],
        baths,
        dt,
        SheetThermalSolve::default(),
    );
    let expected = root(300.0, 900.0, |t| {
        50.0 * (t - 300.0) - dt * gray_heat(t, baths, a)
    });
    let t = result.thermal.temperature_k(CellId(1)).unwrap();
    near(t, expected, 5e-8);
    near(
        result.report.audit.stored_energy_change_j,
        50.0 * (t - 300.0),
        1e-8,
    );
    near(
        result.report.cell_radiant_power[&CellId(1)].value(),
        gray_heat(t, baths, a),
        1e-7,
    );
    near(
        result.report.audit.stored_energy_change_j,
        dt * gray_heat(t, baths, a),
        1e-6,
    );
    assert_eq!(
        result.report.surface_temperatures[1],
        result.report.surface_temperatures[2]
    );
    assert_eq!(result.report.cell_radiant_power.len(), 1);
    assert!(result.report.combined_energy_residual_j.abs() < 1e-5);
    assert!(300.0 + dt * gray_heat(300.0, baths, a) / 50.0 > 2200.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(result.thermal.time_s(), dt);
}

#[test]
fn lossless_transmission_carries_heat_between_baths_without_heating_sheet_and_equilibrium_is_exact()
{
    let initial = state(&SENSIBLE, &[(0.05, 0.0)], None, 0.0);
    let result = gray(
        &initial,
        &[sheet([0.0; 2], 1.0)],
        [900.0, 300.0],
        10.0,
        SheetThermalSolve::default(),
    );
    assert_eq!(result.thermal.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(result.report.cell_radiant_power[&CellId(1)].value(), 0.0);
    assert_eq!(result.report.audit.stored_energy_change_j, 0.0);
    let transfer = STEFAN_BOLTZMANN_W_M2_K4 * (900.0_f64.powi(4) - 300.0_f64.powi(4));
    near(
        result.report.prescribed_surface_heat_power[0].value(),
        -transfer,
        1e-7,
    );
    near(
        result.report.prescribed_surface_heat_power[3].value(),
        transfer,
        1e-7,
    );
    let equilibrium = state(&SENSIBLE, &[(0.05, 300_000.0)], None, 0.0);
    let result = gray(
        &equilibrium,
        &[sheet([0.3, 0.1], 0.5)],
        [600.0; 2],
        100.0,
        SheetThermalSolve::default(),
    );
    assert_eq!(
        result.thermal.cells()[&CellId(1)].specific_h_j_kg,
        300_000.0
    );
    assert_eq!(result.report.audit.stored_energy_change_j, 0.0);
}

#[test]
fn original_unbound_cell_conductive_link_and_authored_source_remain_in_full_residual() {
    let initial = state(&SENSIBLE, &[(0.05, 0.0), (0.1, 0.0)], Some(20.0), 30.0);
    let dt = 5.0;
    let baths = [900.0, 300.0];
    let a = [0.3, 0.1];
    let second = |first| (100.0 * 300.0 + dt * (20.0 * first + 30.0)) / (100.0 + dt * 20.0);
    let expected = root(300.0, 900.0, |t| {
        50.0 * (t - 300.0) - dt * (gray_heat(t, baths, a) - 20.0 * (t - second(t)))
    });
    let result = gray(
        &initial,
        &[sheet(a, 0.5)],
        baths,
        dt,
        SheetThermalSolve::default(),
    );
    near(
        result.thermal.temperature_k(CellId(1)).unwrap(),
        expected,
        5e-8,
    );
    near(
        result.thermal.temperature_k(CellId(2)).unwrap(),
        second(expected),
        5e-8,
    );
    assert_eq!(result.report.cell_radiant_power.len(), 1);
    near(
        result.report.audit.stored_energy_change_j,
        dt * (result.report.cell_radiant_power[&CellId(1)].value() + 30.0),
        1e-6,
    );
    assert_eq!(result.thermal.links()[&1], initial.links()[&1]);
    assert_eq!(initial.temperature_k(CellId(2)).unwrap(), 300.0);
}

#[test]
fn native_latent_heat_is_stored_at_constant_sheet_temperature() {
    let initial = state(&LATENT, &[(1.0, 20_000.0)], None, 0.0);
    let result = gray(
        &initial,
        &[sheet([0.2; 2], 0.5)],
        [400.0, 310.0],
        10.0,
        SheetThermalSolve::default(),
    );
    assert_eq!(result.thermal.temperature_k(CellId(1)).unwrap(), 310.0);
    let gain = result.thermal.cells()[&CellId(1)].specific_h_j_kg - 20_000.0;
    near(
        gain,
        10.0 * gray_heat(310.0, [400.0, 310.0], [0.2; 2]),
        1e-6,
    );
    assert!(gain > 0.0 && gain < 90_000.0);
    near(result.report.audit.stored_energy_change_j, gain, 1e-9);
}

fn spectral() -> SpectralOptics<'static> {
    let wavelength_bounds = [[1e-7, 4e-6], [4e-6, 1e-3]];
    let values = [0.1, 0.8];
    let surface_bands = values
        .map(|a| {
            [1.0, a, a, 1.0]
                .map(|value| {
                    BandSurfaceEmissivity::new(
                        value,
                        [250.0, 2200.0].map(Temperature::new),
                        "declared spectral surface fixture",
                    )
                    .unwrap()
                })
                .to_vec()
        })
        .to_vec();
    SpectralOptics::new_with_sheet_bands(
        wavelength_bounds
            .map(|bounds| SpectralBand {
                wavelengths: bounds.map(Length::new),
                absorption: absorption(),
            })
            .to_vec(),
        [None; 2],
        surface_bands,
        vec![vec![sheet([0.1; 2], 0.8)], vec![sheet([0.8; 2], 0.1)]],
    )
    .unwrap()
}

#[test]
fn selective_spectral_opacity_evolves_original_temperature_with_exact_band_sources() {
    let initial = state(&SENSIBLE, &[(0.05, 100_000.0)], None, 0.0);
    let enclosure = enclosure();
    let optics = spectral();
    let dt = 5.0;
    let result = enclosure
        .candidate_sheet_thermal(
            &initial,
            Time::new(dt),
            &bindings([900.0, 300.0]),
            [Temperature::new(300.0); 2],
            SheetThermalOptics::Spectral {
                optics: &optics,
                solve: SpectralSolve::default(),
            },
            |_, _| panic!("spectral bands own the sources; no gray cards"),
            SheetThermalSolve::default(),
            &mut enclosure.sheet_thermal_workspace(&initial).unwrap(),
            || false,
        )
        .unwrap();
    let bands = [
        (BlackbodyBand::new(1e-7, 4e-6).unwrap(), 0.1),
        (BlackbodyBand::new(4e-6, 1e-3).unwrap(), 0.8),
    ];
    let heat = |t| {
        bands
            .iter()
            .map(|(band, a)| {
                a * (band.evaluate(900.0).unwrap().power_w_m2
                    + band.evaluate(300.0).unwrap().power_w_m2
                    - 2.0 * band.evaluate(t).unwrap().power_w_m2)
            })
            .sum::<f64>()
    };
    let expected = root(400.0, 900.0, |t| 50.0 * (t - 400.0) - dt * heat(t));
    let t = result.thermal.temperature_k(CellId(1)).unwrap();
    near(t, expected, 5e-8);
    near(
        result.report.audit.stored_energy_change_j,
        dt * heat(t),
        1e-6,
    );
    let WallRadiation::Spectral(report) = &result.report.radiation else {
        panic!("spectral report");
    };
    assert!(report.surface_emissivities.is_empty());
    assert_eq!(report.optics.sheet_bands().unwrap().len(), 2);
    assert_eq!(report.factorizations, 2);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.coverage[1].temperature, Temperature::new(t));
    assert_eq!(report.coverage[2].temperature, Temperature::new(t));
    assert!(report.coverage[1].longwave_fraction > 0.0);
    let restricted = SpectralSolve {
        maximum_omitted_blackbody_fraction: 0.0,
        ..SpectralSolve::default()
    };
    assert!(matches!(
        enclosure.candidate_sheet_thermal(
            &initial,
            Time::new(dt),
            &bindings([900.0, 300.0]),
            [Temperature::new(300.0); 2],
            SheetThermalOptics::Spectral {
                optics: &optics,
                solve: restricted
            },
            |_, _| panic!("no gray source"),
            SheetThermalSolve::default(),
            &mut enclosure.sheet_thermal_workspace(&initial).unwrap(),
            || false,
        ),
        Err(SheetThermalError::Radiation(
            RadiationError::SpectralTruncation { .. }
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn small_inventory_shrinks_oversized_probe_without_extrapolating_optical_source() {
    let initial = state(&SENSIBLE, &[(1e-12, 200_000.0)], None, 0.0);
    let pair = SheetPair {
        surfaces: [1, 2],
        optics: ReciprocalSheetOptics::new(
            [0.3, 0.1],
            0.5,
            [499.0, 501.0],
            "narrow-temperature-domain fixture",
        )
        .unwrap(),
    };
    let controls = SheetThermalSolve {
        cell_energy_tolerance_j: 1e-17,
        total_energy_tolerance_j: 1e-17,
        combined_energy_tolerance_j: 1e-17,
        ..SheetThermalSolve::default()
    };
    let result = gray(&initial, &[pair], [900.0, 300.0], 1e-14, controls);
    let expected = root(500.0, 501.0, |t| {
        1e-9 * (t - 500.0) - 1e-14 * gray_heat(t, [900.0, 300.0], [0.3, 0.1])
    });
    near(
        result.thermal.temperature_k(CellId(1)).unwrap(),
        expected,
        2e-8,
    );
    assert!(result.report.evaluations > 20);
    assert!(result.report.maximum_cell_residual_j < 1e-17);
}

#[test]
fn cancellation_budgets_and_split_sheet_ownership_preserve_input_and_reusable_work() {
    let initial = state(&SENSIBLE, &[(0.05, 0.0), (0.1, 0.0)], Some(20.0), 0.0);
    let enclosure = enclosure();
    let mut work = enclosure.sheet_thermal_workspace(&initial).unwrap();
    let pairs = [sheet([0.3, 0.1], 0.5)];
    let absorption = absorption();
    let optics = SheetThermalOptics::Gray {
        absorption: &absorption,
        sheets: &pairs,
        solve: RadiationSolve::default(),
    };
    let surfaces = bindings([900.0, 300.0]);
    let mut polls = 0;
    assert!(matches!(
        enclosure.candidate_sheet_thermal(
            &initial,
            Time::new(5.0),
            &surfaces,
            [Temperature::new(300.0); 2],
            optics,
            black,
            SheetThermalSolve::default(),
            &mut work,
            || {
                polls += 1;
                polls == 20
            },
        ),
        Err(SheetThermalError::Cancelled)
    ));
    assert!(matches!(
        enclosure.candidate_sheet_thermal(
            &initial,
            Time::new(5.0),
            &surfaces,
            [Temperature::new(300.0); 2],
            optics,
            black,
            SheetThermalSolve {
                maximum_evaluations: 1,
                ..SheetThermalSolve::default()
            },
            &mut work,
            || false,
        ),
        Err(SheetThermalError::Budget {
            evaluations: 1,
            maximum_evaluations: 1
        })
    ));
    let mut split = surfaces;
    split[2] = SheetThermalSurface::Stored(CellId(2));
    assert!(matches!(
        enclosure.candidate_sheet_thermal(
            &initial,
            Time::new(5.0),
            &split,
            [Temperature::new(300.0); 2],
            optics,
            black,
            SheetThermalSolve::default(),
            &mut work,
            || false,
        ),
        Err(SheetThermalError::Invalid(_))
    ));
    let accepted = enclosure
        .candidate_sheet_thermal(
            &initial,
            Time::new(5.0),
            &surfaces,
            [Temperature::new(300.0); 2],
            optics,
            black,
            SheetThermalSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(accepted.thermal.temperature_k(CellId(1)).unwrap() > 300.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(initial.cells()[&CellId(2)].specific_h_j_kg, 0.0);
}
