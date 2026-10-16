//! G1/G3 contact-aware Robin response checks against series resistance and
//! independent perturbed production FEM solves. Contact resistance stays fixed.
use super::*;
use crate::{ConductivityModel, ConductionMesh, InterfaceFacePair, InterfaceResistance,
    InterfaceSurface, ThermalBoundaryBuilder, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use crate::fixtures::{box_grid, on_box_face};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 71, kernel_id: 718, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}

fn model<T>(references: [f64; 2], h: [f64; 2],
    f: impl FnOnce(ConductionProblem<'_>, &ThermalInterfaces) -> T) -> T {
    let mut positions = Vec::new();
    let mut tets = Vec::new();
    for side in 0..2 {
        let (complex, points) = box_grid([2, 1, 1], [1.0, 1.0, 1.0]);
        let offset = positions.len() as u32;
        tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
        positions.extend(points.into_iter().map(|[x, y, z]| [x + side as f64, y, z]));
    }
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |face| on_box_face(face.centroid[0], 0.0), ThermalBc::robin(h[0], references[0]).unwrap()).unwrap()
        .region("cold", |face| on_box_face(face.centroid[0], 2.0), ThermalBc::robin(h[1], references[1]).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP, RD),
        value: PropertyValue::Scalar { value: 0.1, dims: RD }, validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(), provenance: Provenance { source: "declared contact derivative fixture".into(),
            license: "internal-test-use".into(), artifact: None } }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId { chemistry: name.into(),
        phase: "solid".into(), process: "as-fixtured".into(), revision: 0 }, texture_frame: "declared".into() };
    let card = InterfaceSystemCard::assemble(side("a"), side("b"), SystemContext {
        medium: "declared".into(), third_body: None, environment: "declared".into(), history: "declared".into(),
    }, claims, Vec::new()).unwrap();
    let resistance = InterfaceResistance::from_card("bond", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).unwrap();
    let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap().into_iter().map(|pair| {
        if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 { pair }
        else { InterfaceFacePair { side_a: pair.side_b, side_b: pair.side_a } }
    }).collect();
    let interfaces = ThermalInterfaces::new(&mesh, &boundary,
        vec![InterfaceSurface::new("bond", pairs, resistance).unwrap()]).unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    f(ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source }, &interfaces)
}

fn config() -> SolveConfig {
    let mut c = SolveConfig::default();
    c.linear.tolerance = 1e-12;
    c.stop.residual_rtol = 1e-12;
    c.stop.step_atol = 0.0;
    c.initial = crate::InitialGuess::Uniform(315.0);
    c
}
fn linear(cx: &Cx<'_>, references: [f64; 2], h: [f64; 2]) -> RobinLinearization {
    model(references, h, |p, i| RobinLinearization::new_with_interfaces(cx, p, i, config(), &["hot", "cold"]).unwrap())
}
fn close(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a:.16e} != {b:.16e}"); }
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

#[test]
fn contact_primal_and_tangent_match_series_resistance_and_jump() {
    with_cx(|cx| {
        let l = linear(cx, [330.0, 300.0], [20.0, 40.0]);
        let resistance = 0.05 + 0.1 + 0.1 + 0.1 + 0.025;
        let heat = 30.0 / resistance;
        let means = l.wall_means(cx, &l.primal().temperature).unwrap();
        close(means[0], 330.0 - heat/20.0, 1e-7);
        close(means[1], 300.0 + heat/40.0, 1e-7);
        model([330.0, 300.0], [20.0, 40.0], |_, interfaces| {
            let flux = interfaces.fluxes(&l.primal().temperature).unwrap();
            close(flux[0].heat_rate_a_to_b_w, heat, 1e-7);
            close(flux[0].mean_jump_k, heat * 0.1, 1e-7);
        });
        let mut d = l.zero_direction(); d.references_k[0] = 1.0;
        let tangent = l.apply(cx, &d).unwrap();
        close(tangent.mean_wall_temperatures_k[0], 1.0 - 1.0/resistance/20.0, 1e-9);
        close(tangent.mean_wall_temperatures_k[1], 1.0/resistance/40.0, 1e-9);
        close(tangent.heat_rates_w[0], -1.0/resistance, 1e-8);
        close(tangent.heat_rates_w[1], 1.0/resistance, 1e-8);
    });
}

#[test]
fn contact_log_h_derivative_matches_perturbed_fem_and_full_transpose() {
    with_cx(|cx| {
        let l = linear(cx, [330.0, 300.0], [20.0, 40.0]);
        let mut d = l.zero_direction(); d.log_htc[1] = 1.0;
        let tangent = l.apply(cx, &d).unwrap();
        let delta = 1e-4_f64;
        let plus = linear(cx, [330.0, 300.0], [20.0, 40.0 * delta.exp()]);
        let minus = linear(cx, [330.0, 300.0], [20.0, 40.0 * (-delta).exp()]);
        for ((&got, &a), &b) in tangent.temperature_k.iter().zip(&plus.primal().temperature).zip(&minus.primal().temperature) {
            close(got, (a-b)/(2.0*delta), 2e-6);
        }
        d.references_k = vec![0.3, -0.7]; d.log_htc = vec![0.2, -0.4];
        d.nodal_load_w[5] = 0.13;
        let tangent = l.apply(cx, &d).unwrap();
        let nodal: Vec<_> = (0..d.nodal_load_w.len()).map(|i| (i as f64 - 4.0)*0.01).collect();
        let walls = [0.7, -0.2]; let heats = [0.03, 0.04];
        let g = l.pullback(cx, &nodal, &walls, &heats).unwrap();
        let left = dot(&nodal, &tangent.temperature_k) + dot(&walls, &tangent.mean_wall_temperatures_k)
            + dot(&heats, &tangent.heat_rates_w);
        let right = dot(&g.references, &d.references_k) + dot(&g.log_htc, &d.log_htc) + dot(&g.nodal_load, &d.nodal_load_w);
        close(left, right, 1e-8);
    });
}

#[test]
fn omitting_contact_still_refuses_and_cancellation_is_preserved() {
    with_cx(|cx| model([330.0, 300.0], [20.0, 40.0], |p, _| {
        assert!(RobinLinearization::new(cx, p, config(), &["hot", "cold"]).is_err());
    }));
    let gate = CancelGate::new_clock_free(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 71, kernel_id: 718, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        model([330.0, 300.0], [20.0, 40.0], |p, i| {
            assert!(matches!(RobinLinearization::new_with_interfaces(&cx, p, i, config(), &["hot", "cold"]), Err(ConductionError::Cancelled { .. })));
        });
    });
}
