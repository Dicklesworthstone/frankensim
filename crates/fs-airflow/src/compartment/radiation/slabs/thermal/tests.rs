use super::*;
use crate::compartment::radiation::{AbsorptionCoefficient, OpticalPath, RadiationError};
use crate::compartment::radiation::spectral::{BandSurfaceEmissivity, SpectralBand};
use fs_conduction::radiation::{
    BlackbodyBand, EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY,
    ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance,
    slab::{AbsorbingLayer, AbsorbingSlab, SlabQuadrature},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Length};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / 1000.0;
        if (250.0..=2200.0).contains(&t) { Ok(t) } else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.001 }
}
static SENSIBLE: Sensible = Sensible;
struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=1e6).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(if h < 10_000.0 { 300.0 + h / 1000.0 }
            else if h <= 110_000.0 { 310.0 }
            else { 310.0 + (h - 110_000.0) / 1000.0 })
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.001 }
}
static LATENT: Latent = Latent;

fn state(
    inventory: &[(f64, f64, u64)],
    links: &[(u64, u64, f64)],
    last_source: f64,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = inventory.iter().enumerate().map(|(i, &(mass, h, material))| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i as u64 + 1), segment_id: i as u64 + 1, body_id: 1,
            material_id: MaterialId(material), mass_kg: mass, specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [0.0, 0.0, i as f64 * 0.01],
                end_m: [0.0, 0.0, (i + 1) as f64 * 0.01],
                volume_m3: 0.01, surface_area_m2: 2.04,
            },
            boundary: Boundary {
                ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None,
                source_w: if i + 1 == inventory.len() { last_source } else { 0.0 },
            },
        }),
    }).collect();
    for (i, &(a, b, conductance)) in links.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64 + 1, a: CellId(a), b: CellId(b),
                area_m2: 1.0, conductance_w_k: conductance,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &SENSIBLE as &dyn ThermalLaw), (MaterialId(2), &LATENT)],
        DepositionBudget { max_cells: 8, max_links: 8, max_events: 16, max_substeps: 100 },
    ).unwrap().advance_to(
        0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-9 }, || Ok(()),
    ).unwrap().0
}
fn enclosure() -> ParticipatingEnclosure {
    let geometry = ViewFactorMatrix::admit(
        vec![1.0; 4],
        vec![
            vec![0.0, 1.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0], vec![0.0, 0.0, 1.0, 0.0],
        ],
        ViewFactorEvidence::Analytic {
            geometry: "black infinite plates surrounding matched-index absorbing material".into(),
        },
        ViewFactorTolerance::default(),
    ).unwrap();
    ParticipatingEnclosure::new(
        geometry,
        [[0, 1], [2, 3]].map(|surfaces| OpticalPath {
            surfaces, exchange_area: Area::new(1.0), segments: vec![],
        }).to_vec(),
        "declared reciprocal vacuum gaps", 1e-12, || false,
    ).unwrap()
}
fn slab(depths: &[f64]) -> SlabPair {
    SlabPair {
        surfaces: [1, 2],
        optics: AbsorbingSlab::new(
            depths.iter().map(|&tau| AbsorbingLayer::new(
                0.01, tau / 0.01, [250.0, 2200.0], "declared analytic slab layer",
            ).unwrap()).collect(),
            SlabQuadrature {
                absolute_coefficient_error: 1e-12,
                relative_coefficient_error: 1e-11,
                ..SlabQuadrature::default()
            }, || false,
        ).unwrap(),
    }
}
fn bindings(baths: [f64; 2]) -> [SlabThermalSurface; 4] {
    [
        SlabThermalSurface::Prescribed(Temperature::new(baths[0])),
        SlabThermalSurface::Stored(CellId(1)),
        SlabThermalSurface::Stored(CellId(2)),
        SlabThermalSurface::Prescribed(Temperature::new(baths[1])),
    ]
}
fn absorption() -> [LayerAbsorption; 2] {
    [0.0; 2].map(|value| LayerAbsorption::new(
        AbsorptionCoefficient::new(value), "declared transparent gas fixture",
    ).unwrap())
}
fn black(surface: usize, temperature: Temperature) -> Result<SurfaceEmissivity, ConductionError> {
    assert!(surface == 0 || surface == 3, "material ports must not query gray surface cards");
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: 1.0, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 2200.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored black bath".into(), license: "internal-test".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic black bath".into(), phase: "solid".into(),
            process: "test".into(), revision: 0,
        }, claims, vec![],
    ).unwrap();
    SurfaceEmissivity::from_card("bath", &card, temperature.value(), SelectionPolicy::SingleClaimOnly)
}
fn gray(
    initial: &DepositionSimulation<'static>,
    depths: &[f64],
    baths: [f64; 2],
    dt: f64,
) -> SlabThermalCandidate<'static, 'static> {
    let enclosure = enclosure();
    enclosure.candidate_slab_thermal(
        initial, Time::new(dt), &bindings(baths), &[vec![CellId(1), CellId(2)]],
        [Temperature::new(300.0); 2],
        SlabThermalOptics::Gray {
            absorption: &absorption(), slabs: &[slab(depths)], solve: RadiationSolve::default(),
        },
        black, SlabThermalSolve::default(), &mut enclosure.slab_thermal_workspace(initial).unwrap(),
        || false,
    ).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}; allowed {tolerance}");
}
fn blackbody(t: f64) -> f64 { STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4) }

// Independent composite Simpson integration of the two directional formal
// ray solutions. Each direction retains its mu through EVERY material layer.
fn ray_heat(depths: &[f64], emission: &[f64], baths: [f64; 2]) -> Vec<f64> {
    const PANELS: usize = 4096;
    let mut heat = vec![0.0; depths.len()];
    for sample in 1..=PANELS {
        let mu = sample as f64 / PANELS as f64;
        let simpson = if sample == PANELS { 1.0 } else if sample % 2 == 0 { 2.0 } else { 4.0 };
        let weight = simpson * 2.0 * mu / (3.0 * PANELS as f64);
        let mut incoming = baths[0];
        for i in 0..depths.len() {
            let absorption = -(-depths[i] / mu).exp_m1();
            let transmission = (-depths[i] / mu).exp();
            heat[i] += weight * absorption * (incoming - emission[i]);
            incoming = incoming * transmission + emission[i] * absorption;
        }
        incoming = baths[1];
        for i in (0..depths.len()).rev() {
            let absorption = -(-depths[i] / mu).exp_m1();
            let transmission = (-depths[i] / mu).exp();
            heat[i] += weight * absorption * (incoming - emission[i]);
            incoming = incoming * transmission + emission[i] * absorption;
        }
    }
    heat
}

#[test]
fn unequal_layer_temperatures_heat_native_cells_and_all_original_conductive_links() {
    let initial = state(
        &[(0.05, 100_000.0, 1), (0.08, 400_000.0, 1), (0.1, 0.0, 1)],
        &[(1, 2, 5.0), (2, 3, 3.0)], 30.0,
    );
    let depths = [0.4, 0.8];
    let dt = 2.0;
    let result = gray(&initial, &depths, [900.0, 300.0], dt);
    let t = [CellId(1), CellId(2), CellId(3)]
        .map(|id| result.thermal.temperature_k(id).unwrap());
    assert!((t[0] - t[1]).abs() > 1.0);
    assert!(t[2] > 300.0);
    let reference = ray_heat(&depths, &[blackbody(t[0]), blackbody(t[1])],
        [blackbody(900.0), blackbody(300.0)]);
    let WallRadiation::Gray(optical) = &result.report.radiation else { panic!("gray report"); };
    for i in 0..2 {
        near(optical.layer_heat_power[0][i].value(), reference[i], 3e-6);
        near(result.report.cell_radiant_power[&CellId(i as u64 + 1)].value(), reference[i], 3e-6);
    }
    near(50.0 * (t[0] - 400.0), dt * (reference[0] + 5.0 * (t[1] - t[0])), 1e-5);
    near(80.0 * (t[1] - 700.0),
        dt * (reference[1] + 5.0 * (t[0] - t[1]) + 3.0 * (t[2] - t[1])), 1e-5);
    near(100.0 * (t[2] - 300.0), dt * (3.0 * (t[1] - t[2]) + 30.0), 1e-6);
    assert_eq!(optical.surface_heat_power[1].value(), 0.0);
    assert_eq!(optical.surface_heat_power[2].value(), 0.0);
    assert_eq!(result.report.cell_radiant_power.len(), 2);
    assert_eq!(optical.layer_temperatures[0], [t[0], t[1]].map(Temperature::new));
    near(result.report.audit.stored_energy_change_j, dt * (reference[0] + reference[1] + 30.0), 1e-5);
    assert!(result.report.combined_energy_residual_j.abs() < 1e-4);
    assert_eq!(result.thermal.links(), initial.links());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 400.0);
}

#[test]
fn transparent_material_has_no_heat_at_unequal_temperatures_and_equilibrium_stays_stationary() {
    let initial = state(&[(0.05, 50_000.0, 1), (0.08, 400_000.0, 1)], &[], 0.0);
    let result = gray(&initial, &[0.0; 2], [900.0, 300.0], 10.0);
    for id in [CellId(1), CellId(2)] {
        assert_eq!(result.thermal.cells()[&id].specific_h_j_kg, initial.cells()[&id].specific_h_j_kg);
        assert_eq!(result.report.cell_radiant_power[&id].value(), 0.0);
    }
    assert_eq!(result.report.audit.stored_energy_change_j, 0.0);
    near(result.report.prescribed_surface_heat_power[0].value(),
        blackbody(300.0) - blackbody(900.0), 1e-7);
    let equilibrium = state(&[(0.05, 300_000.0, 1), (0.08, 300_000.0, 1)], &[(1, 2, 5.0)], 0.0);
    let equal = gray(&equilibrium, &[0.4, 0.8], [600.0; 2], 10.0);
    assert_eq!(equal.report.audit.stored_energy_change_j, 0.0);
    assert_eq!(equal.thermal.temperature_k(CellId(1)).unwrap(), 600.0);
    assert_eq!(equal.thermal.temperature_k(CellId(2)).unwrap(), 600.0);
}

#[test]
fn one_native_latent_layer_absorbs_heat_while_neighbor_keeps_its_own_temperature() {
    let initial = state(&[(0.2, 30_000.0, 2), (0.08, 100_000.0, 1)], &[], 0.0);
    let result = gray(&initial, &[0.5, 0.7], [500.0, 300.0], 0.1);
    assert_eq!(result.thermal.temperature_k(CellId(1)).unwrap(), 310.0);
    assert!(result.thermal.temperature_k(CellId(2)).unwrap() > 310.0);
    let gain = 0.2 * (result.thermal.cells()[&CellId(1)].specific_h_j_kg - 30_000.0);
    assert!(gain > 0.0);
    near(gain, 0.1 * result.report.cell_radiant_power[&CellId(1)].value(), 1e-6);
}

fn spectral() -> SpectralOptics<'static> {
    let bands = [[1e-7, 4e-6], [4e-6, 1e-3]].map(|limits| SpectralBand {
        wavelengths: limits.map(Length::new), absorption: absorption(),
    }).to_vec();
    let surfaces = (0..2).map(|_| (0..2).map(|_| BandSurfaceEmissivity::new(
        1.0, [250.0, 2200.0].map(Temperature::new), "declared spectral black bath",
    ).unwrap()).collect()).collect();
    SpectralOptics::new_with_slab_bands(
        bands, [None; 2], surfaces,
        vec![vec![slab(&[0.2, 1.4])], vec![slab(&[2.0, 0.4])]],
    ).unwrap()
}

#[test]
fn spectral_layer_opacity_heats_original_storage_and_retains_missing_spectrum_admission() {
    let initial = state(&[(0.05, 100_000.0, 1), (0.08, 300_000.0, 1)], &[(1, 2, 3.0)], 0.0);
    let enclosure = enclosure();
    let optics = spectral();
    let dt = 2.0;
    let result = enclosure.candidate_slab_thermal(
        &initial, Time::new(dt), &bindings([900.0, 300.0]), &[vec![CellId(1), CellId(2)]],
        [Temperature::new(300.0); 2],
        SlabThermalOptics::Spectral { optics: &optics, solve: SpectralSolve::default() },
        |_, _| panic!("spectral slab has genuine band sources"),
        SlabThermalSolve::default(), &mut enclosure.slab_thermal_workspace(&initial).unwrap(),
        || false,
    ).unwrap();
    let t = [CellId(1), CellId(2)].map(|id| result.thermal.temperature_k(id).unwrap());
    let mut expected = [0.0; 2];
    for (limits, depths) in [([1e-7, 4e-6], [0.2, 1.4]), ([4e-6, 1e-3], [2.0, 0.4])] {
        let band = BlackbodyBand::new(limits[0], limits[1]).unwrap();
        let emission = t.map(|t| band.evaluate(t).unwrap().power_w_m2);
        let baths = [900.0, 300.0].map(|t| band.evaluate(t).unwrap().power_w_m2);
        let heat = ray_heat(&depths, &emission, baths);
        for i in 0..2 { expected[i] += heat[i]; }
    }
    let WallRadiation::Spectral(report) = &result.report.radiation else { panic!("spectral report"); };
    for i in 0..2 { near(report.layer_heat_power[0][i].value(), expected[i], 3e-6); }
    near(50.0 * (t[0] - 400.0), dt * (expected[0] + 3.0 * (t[1] - t[0])), 1e-5);
    near(80.0 * (t[1] - 600.0), dt * (expected[1] + 3.0 * (t[0] - t[1])), 1e-5);
    assert_eq!(report.factorizations, 2);
    assert_eq!(report.wavelength_evaluations, 0);
    assert!(report.surface_emissivities.is_empty());
    assert_eq!(report.layer_temperatures[0], t.map(Temperature::new));
    assert_eq!(report.optics.slab_bands().unwrap().len(), 2);
    assert!(matches!(
        enclosure.candidate_slab_thermal(
            &initial, Time::new(dt), &bindings([900.0, 300.0]), &[vec![CellId(1), CellId(2)]],
            [Temperature::new(300.0); 2],
            SlabThermalOptics::Spectral { optics: &optics, solve: SpectralSolve {
                maximum_omitted_blackbody_fraction: 0.0, ..SpectralSolve::default()
            } },
            |_, _| panic!("no gray cards"), SlabThermalSolve::default(),
            &mut enclosure.slab_thermal_workspace(&initial).unwrap(), || false,
        ),
        Err(SlabThermalError::Radiation(RadiationError::SpectralTruncation { .. }))
    ));
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn layer_mapping_refuses_wrong_ports_duplicate_or_missing_cells_and_keeps_work_reusable() {
    let initial = state(&[(0.05, 100_000.0, 1), (0.08, 300_000.0, 1)], &[], 0.0);
    let enclosure = enclosure();
    let slabs = [slab(&[0.4, 0.8])];
    let absorption = absorption();
    let optics = SlabThermalOptics::Gray {
        absorption: &absorption, slabs: &slabs, solve: RadiationSolve::default(),
    };
    let surfaces = bindings([900.0, 300.0]);
    let mut work = enclosure.slab_thermal_workspace(&initial).unwrap();
    for ids in [vec![CellId(1)], vec![CellId(1), CellId(1)], vec![CellId(1), CellId(99)], vec![CellId(2), CellId(1)]] {
        assert!(matches!(
            enclosure.candidate_slab_thermal(
                &initial, Time::new(2.0), &surfaces, &[ids], [Temperature::new(300.0); 2],
                optics, black, SlabThermalSolve::default(), &mut work, || false,
            ),
            Err(SlabThermalError::Invalid(_))
        ));
    }
    let ids = [vec![CellId(1), CellId(2)]];
    assert!(matches!(
        enclosure.candidate_slab_thermal(
            &initial, Time::new(2.0), &surfaces, &ids, [Temperature::new(300.0); 2],
            optics, black, SlabThermalSolve::default(), &mut work, || true,
        ),
        Err(SlabThermalError::Cancelled)
    ));
    assert!(matches!(
        enclosure.candidate_slab_thermal(
            &initial, Time::new(2.0), &surfaces, &ids, [Temperature::new(300.0); 2],
            optics, black, SlabThermalSolve { maximum_evaluations: 1, ..SlabThermalSolve::default() },
            &mut work, || false,
        ),
        Err(SlabThermalError::Budget { evaluations: 1, maximum_evaluations: 1 })
    ));
    let result = enclosure.candidate_slab_thermal(
        &initial, Time::new(2.0), &surfaces, &ids, [Temperature::new(300.0); 2],
        optics, black, SlabThermalSolve::default(), &mut work, || false,
    ).unwrap();
    assert_eq!(result.thermal.time_s(), 2.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 400.0);
}
