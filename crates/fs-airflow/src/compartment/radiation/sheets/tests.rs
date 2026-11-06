use super::super::{AbsorptionCoefficient, OpticalPath, OpticalSegment, STEFAN_BOLTZMANN_W_M2_K4};
use super::*;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Length};

fn emissivity(value: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 1500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "analytic opaque bath".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "opaque bath".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("bath", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

fn sheet(a: [f64; 2], t: f64) -> SheetPair {
    SheetPair {
        surfaces: [1, 2],
        optics: ReciprocalSheetOptics::new(a, t, [200.0, 1500.0], "analytic diffuse sheet")
            .unwrap(),
    }
}

fn absorption(k: [f64; 2]) -> [LayerAbsorption; 2] {
    k.map(|value| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(value),
            "explicit uniform gas absorption",
        )
        .unwrap()
    })
}

fn geometry(areas: [f64; 2], gas: bool) -> ParticipatingEnclosure {
    let factors = vec![
        vec![0.0, 1.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0],
        vec![0.0, 0.0, 1.0, 0.0],
    ];
    let matrix = ViewFactorMatrix::admit(
        vec![areas[0], areas[0], areas[1], areas[1]],
        factors,
        ViewFactorEvidence::Analytic {
            geometry: "two opposed black-bath/sheet pairs".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let paths = [[0, 1], [2, 3]]
        .into_iter()
        .enumerate()
        .map(|(side, surfaces)| OpticalPath {
            surfaces,
            exchange_area: Area::new(areas[side]),
            segments: if gas {
                vec![OpticalSegment {
                    layer: if side == 0 {
                        Layer::Lower
                    } else {
                        Layer::Upper
                    },
                    length: Length::new(if side == 0 { 1.0 } else { 2.0 }),
                }]
            } else {
                Vec::new()
            },
        })
        .collect();
    ParticipatingEnclosure::new(
        matrix,
        paths,
        "declared paired-chamber paths",
        1e-12,
        || false,
    )
    .unwrap()
}

fn blackbody(t: f64) -> f64 {
    STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn run(model: &ParticipatingEnclosure, pair: SheetPair, temperatures: [f64; 4]) -> RadiationReport {
    model
        .evaluate_with_sheets(
            &temperatures.map(Temperature::new),
            &[emissivity(1.0), emissivity(1.0)],
            &[pair],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap()
}

fn conserved(report: &RadiationReport) {
    let total: f64 = report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .map(|q| q.value())
        .sum();
    assert!(total.abs() < 1e-7, "net power {total:e}");
    assert_eq!(total, report.energy_residual_w);
}

#[test]
fn transparent_sheet_transfers_between_baths_without_heating() {
    let model = geometry([2.0; 2], false);
    let source = sheet([0.0; 2], 1.0);
    let report = run(&model, source.clone(), [900.0, 500.0, 500.0, 300.0]);
    let heat = 2.0 * (blackbody(900.0) - blackbody(300.0));
    close(report.surface_heat_power[0].value(), -heat);
    close(report.surface_heat_power[3].value(), heat);
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    close(report.surface_radiosity_w_m2[1], blackbody(300.0));
    close(report.surface_radiosity_w_m2[2], blackbody(900.0));
    assert_eq!(report.sheets, vec![source]);
    assert_eq!(report.surface_emissivities.len(), 2);
    conserved(&report);
}

#[test]
fn asymmetric_sheet_matches_face_absorption_and_outgoing_energy() {
    let model = geometry([1.7; 2], false);
    let a = [0.2, 0.4];
    let t = 0.3;
    let source = sheet(a, t);
    let r = source.optics.reflectance();
    let report = run(&model, source, [900.0, 500.0, 500.0, 300.0]);
    let [hot, wall, cold] = [900.0, 500.0, 300.0].map(blackbody);
    let j = [
        a[0] * wall + r[0] * hot + t * cold,
        a[1] * wall + r[1] * cold + t * hot,
    ];
    close(report.surface_radiosity_w_m2[1], j[0]);
    close(report.surface_radiosity_w_m2[2], j[1]);
    close(
        report.surface_heat_power[1].value(),
        1.7 * a[0] * (hot - wall),
    );
    close(
        report.surface_heat_power[2].value(),
        1.7 * a[1] * (cold - wall),
    );
    close(report.surface_heat_power[0].value(), 1.7 * (j[0] - hot));
    close(report.surface_heat_power[3].value(), 1.7 * (j[1] - cold));
    conserved(&report);
}

#[test]
fn asymmetric_sheet_at_equilibrium_has_no_net_heat() {
    let model = geometry([1.0; 2], false);
    let report = run(&model, sheet([0.1, 0.65], 0.25), [700.0; 4]);
    for power in report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
    {
        assert_eq!(power.value(), 0.0);
    }
    for &j in &report.surface_radiosity_w_m2 {
        close(j, blackbody(700.0));
    }
    conserved(&report);
}

#[test]
fn zero_transmission_recovers_original_opaque_solver() {
    let model = geometry([1.0; 2], false);
    let temperatures = [900.0, 500.0, 500.0, 300.0].map(Temperature::new);
    let cards = [0.65, 0.2, 0.6, 0.8].map(emissivity);
    let old = model
        .evaluate(
            &temperatures,
            &cards,
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    let new = model
        .evaluate_with_sheets(
            &temperatures,
            &[cards[0].clone(), cards[3].clone()],
            &[sheet([0.2, 0.6], 0.0)],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    for (a, b) in new.surface_heat_power.iter().zip(&old.surface_heat_power) {
        close(a.value(), b.value());
    }
    for (&a, &b) in new
        .surface_radiosity_w_m2
        .iter()
        .zip(&old.surface_radiosity_w_m2)
    {
        close(a, b);
    }
    assert!(old.sheets.is_empty());
    conserved(&new);
}

#[test]
fn transparent_sheet_carries_radiation_into_the_opposite_gas() {
    let model = geometry([1.3; 2], true);
    let report = model
        .evaluate_with_sheets(
            &[900.0, 500.0, 500.0, 300.0].map(Temperature::new),
            &[emissivity(1.0), emissivity(1.0)],
            &[sheet([0.0; 2], 1.0)],
            [600.0, 350.0].map(Temperature::new),
            &absorption([0.7, 0.9]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    let [hot, cold, lower, upper] = [900.0, 300.0, 600.0, 350.0].map(blackbody);
    let [u, v] = [(-0.7_f64).exp(), (-1.8_f64).exp()];
    let j_left = v * cold + (1.0 - v) * upper;
    let j_right = u * hot + (1.0 - u) * lower;
    close(report.surface_radiosity_w_m2[1], j_left);
    close(report.surface_radiosity_w_m2[2], j_right);
    close(
        report.gas_heat_power[0].value(),
        1.3 * (1.0 - u) * (hot + j_left - 2.0 * lower),
    );
    close(
        report.gas_heat_power[1].value(),
        1.3 * (1.0 - v) * (j_right + cold - 2.0 * upper),
    );
    close(
        report.surface_heat_power[0].value(),
        1.3 * (u * j_left + (1.0 - u) * lower - hot),
    );
    close(
        report.surface_heat_power[3].value(),
        1.3 * (v * j_right + (1.0 - v) * upper - cold),
    );
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    conserved(&report);
}

#[test]
fn incompatible_pair_geometry_state_and_source_axes_refuse() {
    for (areas, temperatures, pairs, cards) in [
        ([1.0, 2.0], [500.0; 4], vec![sheet([0.2; 2], 0.5)], 2),
        (
            [1.0; 2],
            [500.0, 500.0, 501.0, 500.0],
            vec![sheet([0.2; 2], 0.5)],
            2,
        ),
        (
            [1.0; 2],
            [500.0, 190.0, 190.0, 500.0],
            vec![sheet([0.2; 2], 0.5)],
            2,
        ),
        ([1.0; 2], [500.0; 4], vec![sheet([0.2; 2], 0.5); 2], 0),
        ([1.0; 2], [500.0; 4], vec![sheet([0.2; 2], 0.5)], 4),
    ] {
        let model = geometry(areas, false);
        assert!(
            model
                .evaluate_with_sheets(
                    &temperatures.map(Temperature::new),
                    &vec![emissivity(1.0); cards],
                    &pairs,
                    [Temperature::new(300.0); 2],
                    &absorption([0.0; 2]),
                    RadiationSolve::default(),
                    &mut model.workspace().unwrap(),
                    || false,
                )
                .is_err()
        );
    }
}

#[test]
fn lossless_transmitting_or_reflecting_cycle_refuses() {
    let model = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: Vec::new(),
        }],
        "closed sheet-only cycle",
        1e-12,
        || false,
    )
    .unwrap();
    for t in [0.0, 0.6, 1.0] {
        let mut pair = sheet([0.0; 2], t);
        pair.surfaces = [0, 1];
        let result = model.evaluate_with_sheets(
            &[Temperature::new(500.0); 2],
            &[],
            &[pair],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        );
        assert!(matches!(
            result,
            Err(RadiationError::Invalid(
                "closed lossless sheet/surface component has indeterminate radiosity"
            ))
        ));
    }
}
