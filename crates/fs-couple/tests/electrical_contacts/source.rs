//! Synthetic interface fixtures, never material-name defaults.
use fs_evidence::ValidityDomain;
use fs_material::electrical_interface::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::*;
use fs_qty::{Dims,Temperature};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;

pub fn point() -> QueryPoint {
    QueryPoint::new().with("Ta",300.0).unwrap().with("Tb",300.0).unwrap()
        .with("pressure",100000.0).unwrap().with("melt_a",0.0).unwrap().with("melt_b",0.0).unwrap()
}
pub fn card(a:&MaterialCard,b:&MaterialCard,r:f64,slope:f64) -> InterfaceSystemCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve {abscissa:"Ta".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,r*(1.0-100.0*slope)),(2000.0,r*(1.0+1700.0*slope))],dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained().with("Ta",200.0,2000.0).with("Tb",200.0,2000.0)
            .with("pressure",100000.0,100000.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"synthetic electrical contact".into(),license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    InterfaceSystemCard::assemble(SurfaceSpec {material:a.id().clone(),texture_frame:"first finish".into()},
        SurfaceSpec {material:b.id().clone(),texture_frame:"second finish".into()},
        SystemContext {medium:"authored contact".into(),third_body:None,environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap()
}
pub fn contact(id:u64,branch:usize,card:&InterfaceSystemCard,area:f64) -> SeriesElectricalContact {
    SeriesElectricalContact::new(id,branch,resolve_electrical_interface_state(card,&point(),
        MaterialPropertySelection::SingleClaimOnly).unwrap(),area).unwrap()
}

/// Manufactured R'' = r [1 + gamma (Ta - 300) + beta * drive].
/// A true two-coordinate surface exercises coupled temperature/electrical queries.
#[allow(clippy::too_many_arguments)]
pub fn response_card(a: &MaterialCard, b: &MaterialCard, axis: &str, axis_dims: Dims,
    r: f64, gamma: f64, beta: f64, limit: f64) -> InterfaceSystemCard {
    let template = card(a, b, r, 0.0);
    let dims = ELECTRICAL_CONTACT_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Surface2D {
            axes: ["Ta".into(), axis.into()], axis_dims: [Temperature::DIMS, axis_dims],
            samples: [[200.0, -limit], [500.0, -limit], [500.0, limit], [200.0, limit]]
                .into_iter().map(|[t, drive]| ([t, drive], r * (1.0 + gamma * (t - 300.0) + beta * drive))).collect(),
            triangles: vec![[0, 1, 2], [0, 2, 3]], dims,
        },
        validity: ValidityDomain::unconstrained().with("Ta", 200.0, 500.0)
            .with("Tb", 200.0, 500.0).with("pressure", 100000.0, 100000.0)
            .with_quantity(axis, fs_qty::QuantitySpec::dimensional(axis_dims), -limit, limit),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "synthetic joint electrical/thermal contact surface".into(),
            license: "CC0-1.0".into(), artifact: None,
        },
    }).unwrap();
    InterfaceSystemCard::assemble(template.surface_a().clone(), template.surface_b().clone(),
        template.context().clone(), claims, vec![]).unwrap()
}

pub fn response_contact(id: u64, branch: usize, card: &InterfaceSystemCard, area: f64,
    axis: &str, dims: Dims) -> SeriesElectricalContact {
    let point = point().with_quantity(axis, fs_qty::QuantitySpec::dimensional(dims), 0.0).unwrap();
    SeriesElectricalContact::new(id, branch, resolve_electrical_interface_state(card, &point,
        MaterialPropertySelection::SingleClaimOnly).unwrap(), area).unwrap()
}
pub fn scalar_card(property:&str,dims:Dims,value:f64) -> MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(property,dims),value:PropertyValue::Scalar {value,dims},
        validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance {source:"synthetic field storage".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"test".into(),phase:"solid".into(),process:"authored".into(),revision:0},claims,vec![]).unwrap()
}
