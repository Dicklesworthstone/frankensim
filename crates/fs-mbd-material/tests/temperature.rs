//! Material-bound temperatures drive the original expansion and dynamics owners.
mod common;
use common::{close, density_card};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection};
use fs_mbd_material::{MaterialInertiaError, RegionBinding, VolumeMoments, resolve_body};
use fs_mbd_material::thermal::temperature::{RegionTemperature, expand_body_at_temperatures};
use fs_qty::{Pressure, QuantitySpec, Temperature};

fn card(alpha: f64) -> MaterialCard {
    let density = density_card("synthetic", 6.0, 2.0);
    let mut claims = density.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0,alpha),(300.0,2.0*alpha),(400.0,4.0*alpha)], dims: INVERSE_TEMPERATURE_DIMS },
        validity: ValidityDomain::unconstrained().with("T",200.0,400.0).with("pressure",100000.0,100000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source:"synthetic alpha, not measurements".into(),
            license:"CC0-1.0".into(), artifact:None },
    }).unwrap();
    MaterialCard::assemble(density.id().clone(),claims,vec![]).unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap()
        .with_quantity("pressure",QuantitySpec::dimensional(Pressure::DIMS),100000.0).unwrap()
}
fn sample<'a>(card: &'a MaterialCard, key: &'static str, t: f64) -> RegionTemperature<'a> {
    RegionTemperature { key, card, temperature_k:t, selection:MaterialPropertySelection::SingleClaimOnly,
        fixed_point_m:Vec3::new(0.0,0.0,0.0) }
}
#[test]
fn temperature_samples_preserve_full_source_context_and_mass() {
    let card = card(1e-4);
    let point = point();
    let reference = resolve_body(&[RegionBinding { key:"body", material:&card, point:&point,
        geometry:VolumeMoments::cuboid(Vec3::new(1.0,2.0,3.0)).unwrap(),
        selection:MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap();
    let hot = expand_body_at_temperatures(&reference,&[sample(&card,"body",400.0)],||false).unwrap();
    let receipt = hot.regions()[0].expansion();
    close(receipt.free_linear_strain(),0.03);
    assert_eq!(receipt.reference().resolved().query_point(),reference.regions()[0].material().query_point());
    assert_eq!(receipt.current().resolved().axis_quantities(),point.axis_quantities());
    assert_eq!(receipt.current().resolved().query_point().iter().find(|(k,_)|k=="pressure").unwrap().1,100000.0);
    assert_eq!(hot.inertia().mass().to_bits(),reference.inertia().mass().to_bits());
    close(hot.regions()[0].geometry().volume_m3(),6.0*1.03_f64.powi(3));
    let reset = expand_body_at_temperatures(&reference,&[sample(&card,"body",300.0)],||false).unwrap();
    assert_eq!(reset.inertia(),reference.inertia());
    assert_eq!(expand_body_at_temperatures(&reference,&[sample(&card,"body",400.0)],||false).unwrap(),hot);
}
#[test]
fn heterogeneous_temperatures_reorder_without_leaving_stale_regions() {
    let card = card(1e-4);
    let point = point();
    let geometry = VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap();
    let region = |key| RegionBinding { key, material:&card, point:&point, geometry,
        selection:MaterialPropertySelection::SingleClaimOnly };
    let reference = resolve_body(&[region("a"),region("b")],2,||false).unwrap();
    let a = sample(&card,"a",400.0); let b = sample(&card,"b",300.0);
    let hot = expand_body_at_temperatures(&reference,&[a.clone(),b.clone()],||false).unwrap();
    assert_eq!(hot,expand_body_at_temperatures(&reference,&[b,a.clone()],||false).unwrap());
    assert_eq!(hot.regions()[1].inertia(),reference.regions()[1].inertia());
    assert!(expand_body_at_temperatures(&reference,&[a.clone(),a],||false).is_err());
}
#[test]
fn source_refusal_and_final_cancellation_publish_no_configuration() {
    let a = card(1e-4); let foreign = card(2e-4); let p = point();
    let reference = resolve_body(&[RegionBinding { key:"body",material:&a,point:&p,
        geometry:VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap(),
        selection:MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap();
    let unchanged = reference.clone();
    for t in [f64::NAN,0.0,-1.0,401.0] {
        assert!(expand_body_at_temperatures(&reference,&[sample(&a,"body",t)],||false).is_err());
    }
    assert!(expand_body_at_temperatures(&reference,&[sample(&foreign,"body",350.0)],||false).is_err());
    assert!(expand_body_at_temperatures(&reference,&[],||false).is_err());
    assert!(expand_body_at_temperatures(&reference,&[sample(&a,"unknown",350.0)],||false).is_err());
    let mut total = 0;
    expand_body_at_temperatures(&reference,&[sample(&a,"body",350.0)],||{total+=1;false}).unwrap();
    let mut n = 0;
    assert!(matches!(expand_body_at_temperatures(&reference,&[sample(&a,"body",350.0)],||{n+=1;n==total}),
        Err(MaterialInertiaError::Cancelled)));
    assert_eq!(reference,unchanged);
}
