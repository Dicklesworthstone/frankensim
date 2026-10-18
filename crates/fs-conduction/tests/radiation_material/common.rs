#![allow(dead_code)]
//! Synthetic source data and trace carriers, not measured radiative geometry.
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity,
    ViewFactorMatrix, EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::mesh::ConductionMesh;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::*;
use fs_evidence::ValidityDomain;
use fs_qty::{Dims, Pressure};

pub fn card(e0: f64, e1: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0,e0),(1.0,e1)], dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T",300.0,500.0)
            .with_quantity("pressure",fs_qty::QuantitySpec::dimensional(Pressure::DIMS),90000.0,110000.0)
            .with("lot",7.0,7.0).with("melt",0.0,1.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic phase emissivity".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "chart-defined".into(),
        process: "radiation fixture".into(), revision: 0 },claims,vec![]).unwrap()
}
pub fn curve(card: &MaterialCard, latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),
        [(-100.0,300.0,0.0),(0.0,400.0,0.0),(latent,400.0,1.0),(latent+100.0,500.0,1.0)]
        .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
            liquid_mass_fraction:f, bulk_density_kg_m3:1.0 }).to_vec()).unwrap()
}
pub fn context() -> QueryPoint {
    QueryPoint::new().with("T",400.0).unwrap().with("melt",0.0).unwrap().with("lot",7.0).unwrap()
        .with_quantity("pressure",fs_qty::QuantitySpec::dimensional(Pressure::DIMS),100000.0).unwrap()
}
pub fn binding<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ThermalMaterialPoint<'a> {
    ThermalMaterialPoint { material: card, context, temperature_axis:"T", liquid_fraction_axis:Some("melt") }
}
pub fn emissivity(card: &MaterialCard, curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> SurfaceEmissivity {
    SurfaceEmissivity::from_thermal_state("fixture",binding(card,&context()),
        curve.state_at_specific_enthalpy(h).unwrap(),MaterialPropertySelection::SingleClaimOnly).unwrap().0
}
pub fn enclosure(emissivities: [SurfaceEmissivity;2]) -> GrayDiffuseEnclosure {
    let (complex, positions) = box_grid([1,1,1],[1.0,1.0,1.0]);
    let mesh = ConductionMesh::new(complex,positions).unwrap();
    let surfaces = emissivities.into_iter().enumerate().map(|(i,e)| {
        RadiationSurface::new(&mesh,format!("surface-{i}"),
            |face| on_box_face(face.centroid[0],i as f64),e).unwrap()
    }).collect();
    // Algebraic infinite-plate limiting test. Cube faces carry trace/area only;
    // these factors are not claimed as finite-cube ray geometry.
    GrayDiffuseEnclosure::new(surfaces,ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap()).unwrap()
}
pub fn near(a:f64,b:f64) { assert!((a-b).abs()<1e-9* b.abs().max(1.0),"{a} != {b}"); }
