use super::*;
use super::super::{AbsorptionCoefficient, OpticalPath, OpticalSegment, STEFAN_BOLTZMANN_W_M2_K4};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance, slab::{AbsorbingLayer, SlabQuadrature},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Length};

fn emissivity(value: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "analytic gray reservoir".into(), license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(), phase: "solid".into(),
            process: "opaque reservoir".into(), revision: 0,
        }, claims, Vec::new(),
    ).unwrap();
    SurfaceEmissivity::from_card("bath", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

fn slab(depths: &[f64]) -> SlabPair {
    SlabPair {
        surfaces: [1, 2],
        optics: AbsorbingSlab::new(
            depths.iter().map(|&depth| AbsorbingLayer::new(
                1.0, depth, [200.0, 1500.0], "explicit constant absorption coefficient",
            ).unwrap()).collect(),
            SlabQuadrature::default(), || false,
        ).unwrap(),
    }
}

fn absorption(k: [f64; 2]) -> [LayerAbsorption; 2] {
    k.map(|k| LayerAbsorption::new(
        AbsorptionCoefficient::new(k), "explicit uniform gas coefficient",
    ).unwrap())
}

fn geometry(areas: [f64; 2], gas: bool) -> ParticipatingEnclosure {
    let matrix = ViewFactorMatrix::admit(
        vec![areas[0], areas[0], areas[1], areas[1]],
        vec![vec![0.0, 1.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0],
             vec![0.0, 0.0, 0.0, 1.0], vec![0.0, 0.0, 1.0, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "two diffuse slab ports facing reservoirs".into() },
        ViewFactorTolerance::default(),
    ).unwrap();
    let paths = [[0, 1], [2, 3]].into_iter().enumerate().map(|(side, surfaces)| OpticalPath {
        surfaces, exchange_area: Area::new(areas[side]),
        segments: if gas { vec![OpticalSegment {
            layer: if side == 0 { Layer::Lower } else { Layer::Upper },
            length: Length::new(1.0),
        }] } else { Vec::new() },
    }).collect();
    ParticipatingEnclosure::new(matrix, paths, "authored diffuse-port paths", 1e-12, || false).unwrap()
}

fn blackbody(t: f64) -> f64 { STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4) }
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}");
}
fn conserved(report: &RadiationReport) {
    let total: f64 = report.surface_heat_power.iter().chain(report.gas_heat_power.iter())
        .chain(report.layer_heat_power.iter().flatten()).map(|q| q.value()).sum();
    assert!(total.abs() < 1e-7, "net heat {total:e}");
    assert_eq!(report.energy_residual_w, total);
}
fn run(model: &ParticipatingEnclosure, pair: SlabPair, temperatures: &[f64]) -> RadiationReport {
    model.evaluate_with_slabs(
        &[900.0, temperatures[0], temperatures[temperatures.len()-1], 250.0].map(Temperature::new),
        &[1.0; 2].map(emissivity), &[pair],
        &[temperatures.iter().copied().map(Temperature::new).collect()],
        [Temperature::new(300.0); 2], &absorption([0.0; 2]), RadiationSolve::default(),
        &mut model.workspace().unwrap(), || false,
    ).unwrap()
}

#[test]
fn nonisothermal_two_layer_heat_matches_independent_exponential_integral_moments() {
    // Exact diffuse transmission: 2E3(0.5), 2E3(1), independently tabulated.
    let half = 0.443_208_728_550_356_9;
    let full = 0.219_383_934_395_520_27;
    let near = 1.0 - half;
    let far = half - full;
    let internal = 1.0 - 2.0 * half + full;
    for with_gas in [false, true] {
        let model = geometry([1.3; 2], with_gas);
        let pair = slab(&[0.5, 0.5]);
        let report = model.evaluate_with_slabs(
            &[900.0, 700.0, 400.0, 250.0].map(Temperature::new),
            &[1.0; 2].map(emissivity), &[pair.clone()],
            &[vec![Temperature::new(700.0), Temperature::new(400.0)]],
            [600.0, 350.0].map(Temperature::new),
            &absorption(if with_gas { [0.7, 0.9] } else { [0.0; 2] }),
            RadiationSolve::default(), &mut model.workspace().unwrap(), || false,
        ).unwrap();
        let [hot, first, second, cold, lower, upper] =
            [900.0, 700.0, 400.0, 250.0, 600.0, 350.0].map(blackbody);
        let u = if with_gas { (-0.7_f64).exp() } else { 1.0 };
        let v = if with_gas { (-0.9_f64).exp() } else { 1.0 };
        let h_left = u * hot + (1.0 - u) * lower;
        let h_right = v * cold + (1.0 - v) * upper;
        let j_left = full * h_right + near * first + far * second;
        let j_right = full * h_left + far * first + near * second;
        close(report.surface_radiosity_w_m2[1], j_left);
        close(report.surface_radiosity_w_m2[2], j_right);
        close(report.layer_heat_power[0][0].value(), 1.3 * (
            near * (h_left - first) + far * (h_right - first) + internal * (second - first)));
        close(report.layer_heat_power[0][1].value(), 1.3 * (
            far * (h_left - second) + near * (h_right - second) + internal * (first - second)));
        close(report.gas_heat_power[0].value(), 1.3 * (1.0 - u) * (hot + j_left - 2.0 * lower));
        close(report.gas_heat_power[1].value(), 1.3 * (1.0 - v) * (cold + j_right - 2.0 * upper));
        close(report.surface_heat_power[0].value(), 1.3 * (u * j_left + (1.0 - u) * lower - hot));
        close(report.surface_heat_power[3].value(), 1.3 * (v * j_right + (1.0 - v) * upper - cold));
        assert_eq!(report.surface_heat_power[1].value(), 0.0);
        assert_eq!(report.surface_heat_power[2].value(), 0.0);
        assert_eq!(report.slabs, vec![pair]);
        conserved(&report);
    }
}

#[test]
fn transparent_layers_pass_radiation_and_receive_exactly_zero_heat() {
    let model = geometry([2.0; 2], false);
    let report = run(&model, slab(&[0.0; 3]), &[300.0, 800.0, 500.0]);
    for q in &report.layer_heat_power[0] { assert_eq!(q.value(), 0.0); }
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    close(report.surface_heat_power[0].value(), 2.0 * (blackbody(250.0)-blackbody(900.0)));
    close(report.surface_heat_power[3].value(), -report.surface_heat_power[0].value());
    conserved(&report);
}

#[test]
fn opaque_layers_emit_locally_and_exchange_heat_internally() {
    let model = geometry([1.0; 2], false);
    let report = run(&model, slab(&[1000.0; 2]), &[700.0, 400.0]);
    let [hot, first, second, cold] = [900.0, 700.0, 400.0, 250.0].map(blackbody);
    close(report.surface_radiosity_w_m2[1], first);
    close(report.surface_radiosity_w_m2[2], second);
    close(report.layer_heat_power[0][0].value(), hot + second - 2.0*first);
    close(report.layer_heat_power[0][1].value(), first + cold - 2.0*second);
    conserved(&report);
}

#[test]
fn isothermal_subdivision_preserves_full_stack_angular_memory() {
    let model = geometry([1.0; 2], false);
    let whole = run(&model, slab(&[1.0]), &[600.0]);
    let split = run(&model, slab(&[0.5, 0.5]), &[600.0, 600.0]);
    for (a, b) in split.surface_heat_power.iter().zip(&whole.surface_heat_power) { close(a.value(), b.value()); }
    for (&a, &b) in split.surface_radiosity_w_m2.iter().zip(&whole.surface_radiosity_w_m2) { close(a,b); }
    close(split.layer_heat_power[0].iter().map(|q|q.value()).sum(), whole.layer_heat_power[0][0].value());
    conserved(&split);
}

#[test]
fn full_thermal_equilibrium_has_zero_layer_gas_and_surface_heat() {
    let model = geometry([1.0; 2], true);
    let report = model.evaluate_with_slabs(
        &[Temperature::new(700.0); 4], &[0.6; 2].map(emissivity), &[slab(&[0.2, 3.0])],
        &[vec![Temperature::new(700.0); 2]], [Temperature::new(700.0); 2],
        &absorption([0.3, 0.8]), RadiationSolve::default(),
        &mut model.workspace().unwrap(), || false,
    ).unwrap();
    for q in report.surface_heat_power.iter().chain(report.gas_heat_power.iter())
        .chain(report.layer_heat_power.iter().flatten()) { assert_eq!(q.value(),0.0); }
    conserved(&report);
}

#[test]
fn slab_source_geometry_temperature_and_lossless_axes_are_enforced() {
    for (areas, face, temperatures) in [
        ([1.0, 2.0], [900.0,700.0,400.0,250.0], vec![700.0,400.0]),
        ([1.0;2], [900.0,701.0,400.0,250.0], vec![700.0,400.0]),
        ([1.0;2], [900.0,700.0,400.0,250.0], vec![700.0]),
        ([1.0;2], [900.0,700.0,190.0,250.0], vec![700.0,190.0]),
    ] {
        let model = geometry(areas, false);
        assert!(model.evaluate_with_slabs(
            &face.map(Temperature::new), &[1.0; 2].map(emissivity), &[slab(&[0.5;2])],
            &[temperatures.into_iter().map(Temperature::new).collect()],
            [Temperature::new(300.0);2], &absorption([0.0;2]), RadiationSolve::default(),
            &mut model.workspace().unwrap(), ||false,
        ).is_err());
    }
    let closed = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath { surfaces:[0,1], exchange_area:Area::new(1.0), segments:Vec::new() }],
        "lossless closed ports", 1e-12, ||false,
    ).unwrap();
    let mut pair = slab(&[0.0]);
    pair.surfaces = [0,1];
    assert!(closed.evaluate_with_slabs(
        &[Temperature::new(500.0);2], &[], &[pair], &[vec![Temperature::new(500.0)]],
        [Temperature::new(300.0);2], &absorption([0.0;2]), RadiationSolve::default(),
        &mut closed.workspace().unwrap(), ||false,
    ).is_err());
}
