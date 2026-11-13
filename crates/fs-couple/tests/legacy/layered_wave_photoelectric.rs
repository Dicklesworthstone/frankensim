#![cfg(feature = "layered-wave-heating")]
//! G0/G1/G3 native cell/escape-bin coupling on explicitly synthetic materials.

use fs_couple::layered_wave_heating::{
    LayeredWaveDrive, LayeredWaveHeating, LayeredWaveHeatingError,
    photoelectric::{
        LayeredPhotoelectricError, LocalThermalization, PreparedLayeredPhotoelectricHeating,
    },
};
use fs_em::plane_wave::{
    DielectricResponse, Polarization,
    layers::{IlluminationBudget, photoelectric::EscapeFace},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
    },
    photoelectric::{
        ELEMENTARY_CHARGE_C, LIGHT_M_S, PLANCK_J_S,
        source::{
            PhotoelectricKnot,
            depth::{DepthPhotoelectricMaterial, PhotoelectricDepthBin},
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
        SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_NAME, ScalarAdmissibility, ScalarPropertyRequirement,
        THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY as K_NAME,
        resolve_material_state_point,
    },
    thermal::{
        conduction::LayerConductionConfig,
        transport::{
            stack::{
                LayerStackSourceQuery,
                transient::{
                    BoundLayeredThermalBody, LayeredThermalIds, PreparedLayeredThermalBody,
                    TransientLayerSourceQuery,
                    implicit::{LayeredThermalSolveControls, LayeredThermalWorkspace},
                },
            },
            storage::{HeatStorageAmount, HeatStorageComponentQuery},
        },
    },
};
use fs_math::c64::C64;
use fs_qty::{
    Dims, Length, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{DepositionBudget, DepositionSimulation, Event, StepConfig};

const WIDTH: f64 = 12e-9;
const AREA: f64 = 0.001;
const WORK: f64 = 3e-19;
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, 0]);
const DT: f64 = 0.001;
const INPUT: [C64; 2] = [C64::new(100.0, -30.0), C64::new(-25.0, 15.0)];

fn at(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn properties(vacuum: bool, first_yield: f64, deep_yield: f64) -> Vec<(&'static str, Dims, f64)> {
    let basis = ElectricalPropertyBasis::Conductivity;
    vec![
        (
            EPS_NAME,
            EPS,
            8.854_187_812_8e-12 * if vacuum { 1.0 } else { 2.0 },
        ),
        (MU_NAME, MU, 1.256_637_062_12e-6),
        (
            basis.property_name(),
            basis.dimensions(),
            if vacuum { 0.0 } else { 1e5 },
        ),
        (K_NAME, THERMAL_CONDUCTIVITY_DIMS, 1e-8),
        (CP_NAME, SPECIFIC_HEAT_CAPACITY_DIMS, 200.0),
        ("work", ENERGY, WORK),
        ("film_width", Length::DIMS, WIDTH),
        ("first_yield", Dims::NONE, first_yield),
        ("deep_yield", Dims::NONE, deep_yield),
        ("first_mean", Dims::NONE, 0.3),
        ("deep_mean", Dims::NONE, 0.6),
    ]
}
fn card(name: &str, vacuum: bool, first_yield: f64, deep_yield: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in properties(vacuum, first_yield, deep_yield) {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    280.0,
                    360.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic depth-photoemission/thermal fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "authored numerical film".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn state(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    let requirements: Vec<_> = properties(false, 0.0, 0.0)
        .into_iter()
        .map(|(name, dims, _)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        &at(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn film(
    card: &MaterialCard,
    boundary: f64,
) -> (ResolvedConductorState, DepthPhotoelectricMaterial) {
    let resolved = state(card);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        resolved.clone(),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    let bins: Vec<_> = [
        (boundary, "first_yield", "first_mean"),
        (1.0, "deep_yield", "deep_mean"),
    ]
    .map(
        |(end_fraction, yield_key, mean_key)| PhotoelectricDepthBin {
            end_fraction,
            knots: [200.0, 2000.0]
                .map(|wavelength_nm| PhotoelectricKnot {
                    wavelength_nm,
                    absorbed_quantum_yield_property: PropertyKey::new(yield_key, Dims::NONE),
                    mean_excess_fraction_property: PropertyKey::new(mean_key, Dims::NONE),
                })
                .to_vec(),
        },
    )
    .to_vec();
    let source = DepthPhotoelectricMaterial::try_from_resolved(
        resolved,
        PropertyKey::new("work", ENERGY),
        PropertyKey::new("film_width", Length::DIMS),
        &bins,
        1800.0,
    )
    .unwrap();
    (conductor, source)
}
fn body(cards: &[&MaterialCard], counts: &[usize]) -> PreparedLayeredThermalBody {
    let lower = at(280.0);
    let upper = at(360.0);
    let conductivity = PropertyKey::new(K_NAME, THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new(CP_NAME, SPECIFIC_HEAT_CAPACITY_DIMS);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let mass = 1000.0 * AREA * WIDTH;
    let queries: Vec<_> = cards
        .iter()
        .zip(counts)
        .enumerate()
        .map(|(i, (&card, &cell_count))| TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                id: 100 + i as u64,
                card,
                conductivity: &conductivity,
                thickness_m: WIDTH,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            capacity: HeatStorageComponentQuery {
                id: 200 + i as u64,
                card,
                capacity: &capacity,
                amount: HeatStorageAmount::MassKg(mass),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            mass_kg: mass,
            cell_count,
        })
        .collect();
    PreparedLayeredThermalBody::try_new(&queries, &vec![0.0; cards.len() - 1], AREA, 300.0).unwrap()
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 10,
        first_material: 30,
        first_link: 50,
        body_id: 70,
    }
}
fn conduction() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 1e-12,
    }
}
fn controls() -> LayeredThermalSolveControls {
    LayeredThermalSolveControls {
        cell_tolerance_j: 1e-13,
        total_tolerance_j: 1e-13,
        ..Default::default()
    }
}
fn initial<'a>(bound: &'a BoundLayeredThermalBody<'_>, source_w: f64) -> DepositionSimulation<'a> {
    let n = bound.prepared_body().cell_count();
    let mut events = bound
        .birth_events(0.0, &vec![300.0; n], conduction(), || Ok(()))
        .unwrap();
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            cell.boundary.source_w = source_w;
        }
    }
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-13,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn omega(wavelength_m: f64) -> f64 {
    core::f64::consts::TAU * LIGHT_M_S / wavelength_m
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tol={tolerance:e}"
    );
}
fn relative(a: f64, b: f64) {
    near(a, b, 3e-10 * b.abs().max(1e-300));
}

#[test]
fn g1_cell_bin_intersections_preserve_full_optics_and_deposit_only_residual_heat() {
    let exterior = state(&card("vacuum", true, 0.0, 0.0));
    let cards = [
        card("emitter", false, 0.8, 0.02),
        card("ordinary absorber", false, 0.8, 0.02),
    ];
    let (metal, source) = film(&cards[0], 0.41);
    let (backing, _) = film(&cards[1], 0.41);
    let native = body(&[&cards[0], &cards[1]], &[3, 2]);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let start = initial(&bound, 2e-5);
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[
            DielectricResponse::Ohmic(&metal),
            DielectricResponse::Ohmic(&backing),
        ],
        || false,
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let wave = owner.prepare(0.73, polarization, || false).unwrap();
        let photo = PreparedLayeredPhotoelectricHeating::try_new(
            &wave,
            0,
            &source,
            EscapeFace::EntranceVacuum,
            LocalThermalization::AtAbsorptionPosition,
            || false,
        )
        .unwrap();
        assert_eq!(photo.prepared_interval_count(), 4);
        let mut workspace = LayeredThermalWorkspace::new(5).unwrap();
        let result = photo
            .advance_to(
                &start,
                DT,
                LayeredWaveDrive::Coherent(INPUT),
                conduction(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        let optical = owner
            .optical_layers()
            .prepare_illumination(0.73, polarization)
            .unwrap();
        let mut original = [0.0; 5];
        optical.evaluate(INPUT, &mut original).unwrap();
        let mut expected: [f64; 5] = original.map(|p| p * AREA * DT);
        expected[..3].fill(0.0);
        let (mut emitted, mut removal, mut kinetic) = (0.0, 0.0, 0.0);
        let photon = PLANCK_J_S * LIGHT_M_S / 250e-9;
        let width = WIDTH / 3.0;
        // The authored 41% boundary cuts only cell 1, at 23% of its width.
        for (cell, a, b, yield_, mean) in [
            (0, 0.0, width, 0.8, 0.3),
            (1, 0.0, 0.23 * width, 0.8, 0.3),
            (1, 0.23 * width, width, 0.02, 0.6),
            (2, 0.0, width, 0.02, 0.6),
        ] {
            let energy = optical
                .prepare_depth_interval(cell, a, b)
                .unwrap()
                .evaluate(INPUT)
                .unwrap()
                * AREA
                * DT;
            let count = energy * yield_ / photon;
            emitted += count;
            removal += count * WORK;
            kinetic += count * mean * (photon - WORK);
            expected[cell] += energy * (1.0 - yield_) + count * (1.0 - mean) * (photon - WORK);
        }
        for ((id, actual), (&expected, original_id)) in result
            .cell_energy_j
            .iter()
            .zip(expected.iter().zip(owner.cell_ids()))
        {
            assert_eq!(id, original_id);
            relative(*actual, expected);
        }
        relative(
            result.optical.absorbed_j,
            original.iter().sum::<f64>() * AREA * DT,
        );
        relative(
            result.photoelectric.absorbed_emitter_j,
            original[..3].iter().sum::<f64>() * AREA * DT,
        );
        relative(result.photoelectric.removal_j, removal);
        relative(result.photoelectric.kinetic_j, kinetic);
        relative(result.photoelectric.emitted_electrons, emitted);
        relative(
            result.photoelectric.surface_charge_c,
            emitted * ELEMENTARY_CHARGE_C,
        );
        assert_eq!(
            result.photoelectric.surface_charge_c,
            -result.photoelectric.emitted_charge_c
        );
        relative(result.photoelectric.residual_heat_j, expected.iter().sum());
        near(result.photoelectric.partition_residual_j, 0.0, 1e-17);
        near(result.optical.residual_j, 0.0, 1e-17);
        near(
            result.thermal.energy_j() - start.energy_j(),
            expected.iter().sum::<f64>() + 5.0 * 2e-5 * DT,
            3e-13,
        );
        assert!(result.optical.absorbed_j > result.photoelectric.residual_heat_j);
        assert!(std::ptr::eq(result.source, &source));
        assert_eq!(start.time_s(), 0.0);
    }
}

#[test]
fn g0_zero_yield_and_subthreshold_light_recover_native_optical_heating() {
    let exterior = state(&card("vacuum", true, 0.0, 0.0));
    for (yield_, wavelength) in [(0.0, 250e-9), (0.8, 1000e-9)] {
        let card = card("emitter", false, yield_, yield_);
        let (metal, source) = film(&card, 0.41);
        let native = body(&[&card], &[2]);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let start = initial(&bound, 0.0);
        let owner = LayeredWaveHeating::try_new(
            &bound,
            &exterior,
            DielectricResponse::Lossless(&exterior),
            omega(wavelength),
            &[DielectricResponse::Ohmic(&metal)],
            || false,
        )
        .unwrap();
        let wave = owner.prepare(1.0, Polarization::Te, || false).unwrap();
        let photo = PreparedLayeredPhotoelectricHeating::try_new(
            &wave,
            0,
            &source,
            EscapeFace::EntranceVacuum,
            LocalThermalization::AtAbsorptionPosition,
            || false,
        )
        .unwrap();
        let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
        let result = photo
            .advance_to(
                &start,
                DT,
                LayeredWaveDrive::Coherent(INPUT),
                conduction(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        let original = wave
            .advance_to(
                &start,
                DT,
                LayeredWaveDrive::Coherent(INPUT),
                conduction(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(result.photoelectric.emitted_electrons, 0.0);
        assert_eq!(result.photoelectric.removal_j, 0.0);
        assert_eq!(result.photoelectric.kinetic_j, 0.0);
        relative(
            result.photoelectric.residual_heat_j,
            original.optical.absorbed_j,
        );
        for ((id, heat), (_, expected)) in result.cell_energy_j.iter().zip(&original.cell_energy_j)
        {
            relative(*heat, *expected);
            near(
                result.thermal.temperature_k(*id).unwrap(),
                original.thermal.temperature_k(*id).unwrap(),
                1e-9,
            );
        }
    }
}

#[test]
fn g1_exit_face_reverses_profile_on_original_cells_and_buried_source_refuses() {
    let exterior = state(&card("vacuum", true, 0.0, 0.0));
    let cards = [
        card("first", false, 0.8, 0.02),
        card("last", false, 0.8, 0.02),
    ];
    let (front, wrong_source) = film(&cards[0], 0.25);
    let (metal, source) = film(&cards[1], 0.25);
    let native = body(&[&cards[0], &cards[1]], &[1, 4]);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let start = initial(&bound, 0.0);
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[
            DielectricResponse::Ohmic(&front),
            DielectricResponse::Ohmic(&metal),
        ],
        || false,
    )
    .unwrap();
    let wave = owner.prepare(1.0, Polarization::Te, || false).unwrap();
    assert!(
        PreparedLayeredPhotoelectricHeating::try_new(
            &wave,
            0,
            &wrong_source,
            EscapeFace::ExitVacuum,
            LocalThermalization::AtAbsorptionPosition,
            || false
        )
        .is_err()
    );
    assert!(
        PreparedLayeredPhotoelectricHeating::try_new(
            &wave,
            1,
            &wrong_source,
            EscapeFace::ExitVacuum,
            LocalThermalization::AtAbsorptionPosition,
            || false
        )
        .is_err()
    );
    // A positive near-face bin must not vanish when reversing its coordinate.
    let (_, collapsed) = film(&cards[1], 1e-20);
    assert!(matches!(
        PreparedLayeredPhotoelectricHeating::try_new(
            &wave,
            1,
            &collapsed,
            EscapeFace::ExitVacuum,
            LocalThermalization::AtAbsorptionPosition,
            || false,
        ),
        Err(LayeredPhotoelectricError::Invalid(_))
    ));
    let photo = PreparedLayeredPhotoelectricHeating::try_new(
        &wave,
        1,
        &source,
        EscapeFace::ExitVacuum,
        LocalThermalization::AtAbsorptionPosition,
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(5).unwrap();
    let result = photo
        .advance_to(
            &start,
            DT,
            LayeredWaveDrive::Coherent(INPUT),
            conduction(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let mut original = [0.0; 5];
    wave.cell_power_into(LayeredWaveDrive::Coherent(INPUT), &mut original, || false)
        .unwrap();
    let count = (0.02 * original[1..4].iter().sum::<f64>() + 0.8 * original[4]) * DT
        / (PLANCK_J_S * LIGHT_M_S / 250e-9);
    relative(result.photoelectric.emitted_electrons, count);
    relative(result.cell_energy_j[0].1, original[0] * DT);
    assert_eq!(result.escape_face, EscapeFace::ExitVacuum);
    assert_eq!(result.emitting_layer_index, 1);
}

#[test]
fn g3_native_covariance_global_work_budget_and_cancelled_candidate_remain_atomic() {
    let exterior = state(&card("vacuum", true, 0.0, 0.0));
    let card = card("emitter", false, 0.8, 0.02);
    let (metal, source) = film(&card, 0.41);
    let native = body(&[&card], &[2]);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let start = initial(&bound, 0.0);
    let energy = start.energy_j();
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricResponse::Ohmic(&metal)],
        || false,
    )
    .unwrap();
    let wave = owner.prepare(1.0, Polarization::Tm, || false).unwrap();
    let photo = PreparedLayeredPhotoelectricHeating::try_new(
        &wave,
        0,
        &source,
        EscapeFace::EntranceVacuum,
        LocalThermalization::AtAbsorptionPosition,
        || false,
    )
    .unwrap();
    let factors = [[INPUT[0], C64::ZERO], [C64::ZERO, INPUT[1]]];
    let pairs = 2 * (2 + photo.prepared_interval_count());
    let budget = IlluminationBudget {
        maximum_modes: 2,
        maximum_layer_mode_pairs: pairs,
    };
    let drive = LayeredWaveDrive::Coherence {
        factors: &factors,
        budget,
    };
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let result = photo
        .advance_to(
            &start,
            DT,
            drive,
            conduction(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let independent = photo
        .advance_to(
            &start,
            DT,
            LayeredWaveDrive::Independent(INPUT),
            conduction(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    relative(
        result.photoelectric.emitted_electrons,
        independent.photoelectric.emitted_electrons,
    );
    relative(
        result.photoelectric.residual_heat_j,
        independent.photoelectric.residual_heat_j,
    );
    assert!(matches!(
        photo.advance_to(
            &start,
            DT,
            LayeredWaveDrive::Coherence {
                factors: &factors,
                budget: IlluminationBudget {
                    maximum_layer_mode_pairs: pairs - 1,
                    ..budget
                }
            },
            conduction(),
            controls(),
            &mut workspace,
            || false
        ),
        Err(LayeredPhotoelectricError::Invalid(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        photo.advance_to(
            &start,
            DT,
            drive,
            conduction(),
            controls(),
            &mut workspace,
            || {
                polls += 1;
                polls >= 5
            }
        ),
        Err(LayeredPhotoelectricError::Heating(
            LayeredWaveHeatingError::Cancelled
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.energy_j(), energy);
    let retry = photo
        .advance_to(
            &start,
            DT,
            drive,
            conduction(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.photoelectric, retry.photoelectric);
    assert_eq!(result.thermal.energy_j(), retry.thermal.energy_j());
}
