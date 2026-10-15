//! Existing source integrator -> conserved reference mass -> existing dynamics.
mod common;
use common::{close, density_card, point};
use fs_evidence::ValidityDomain;
use fs_ga::{Mat3, Se3, Vec3, Wrench};
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    integrate_isotropic_thermal_expansion};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_dynamics};
use fs_mbd_material::{MaterialInertiaError, RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_mbd_material::thermal::{ThermalExpansionBinding, expand_body};
use fs_qty::Temperature;

fn thermal_card(name: &str, alpha: f64) -> MaterialCard {
    let density = density_card(name, 6.0, 2.0);
    let mut claims = density.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, alpha), (300.0, 2.0*alpha), (400.0, 4.0*alpha)],
            dims: INVERSE_TEMPERATURE_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T",200.0,400.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance { source: "synthetic expansion fixture".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(density.id().clone(), claims, vec![]).unwrap()
}
fn body(card: &MaterialCard) -> ResolvedBody {
    resolve_body(&[RegionBinding {
        key: "body", geometry: VolumeMoments::cuboid(Vec3::new(1.0,2.0,3.0)).unwrap(),
        material: card, point: &point(300.0), selection: MaterialPropertySelection::SingleClaimOnly,
    }],1,||false).unwrap()
}
fn path(card: &MaterialCard, from: f64, to: f64) -> IntegratedIsotropicThermalExpansion {
    integrate_isotropic_thermal_expansion(card,&point(from),&point(to),
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn binding(path: &IntegratedIsotropicThermalExpansion) -> ThermalExpansionBinding<'_> {
    ThermalExpansionBinding { key: "body", expansion: path, fixed_point_m: Vec3::new(0.0,0.0,0.0) }
}
fn matrix_close(a: Mat3, b: Mat3, scale: f64) {
    for (x,y) in a.m.into_iter().zip(b.m) { close(x,y*scale); }
}

#[test]
fn source_curve_drives_geometry_without_requerying_density_or_changing_mass() {
    let card = thermal_card("A",1e-4);
    let reference = body(&card);
    let receipt = path(&card,300.0,400.0);
    // Integral of linearly varying alpha, NOT alpha(end)*delta_T.
    close(receipt.free_linear_strain(),0.03);
    let s = 1.03_f64;
    let expanded = expand_body(&reference,&[binding(&receipt)],||false).unwrap();
    assert_eq!(expanded.inertia().mass().to_bits(),reference.inertia().mass().to_bits());
    let region = &expanded.regions()[0];
    close(region.geometry().volume_m3(),6.0*s.powi(3));
    close(region.density_kg_m3(),4.0/s.powi(3));
    matrix_close(region.inertia().inertia_com(),reference.inertia().inertia_com(),s*s);
    matrix_close(region.geometry().inertia_m5(),reference.regions()[0].geometry().inertia_m5(),s.powi(5));
    assert_eq!(region.reference_index(),0);
    assert_eq!(region.expansion(),&receipt);
    assert!(std::ptr::eq(expanded.reference(),&reference));
    // The deliberately different density(T) table is not queried at 400 K to
    // turn a fixed 24 kg inventory into 12 kg. Derived rho is labelled as such.
    close(reference.inertia().mass(),24.0);
}

fn rotor_acceleration(inertia: SpatialInertia) -> f64 {
    let model = ArticulatedModel::new(vec![
        Link::new("root",None,Se3::identity(),JointModel::FIXED,inertia),
        Link::new("rotor",Some(0),Se3::identity(),
            JointModel::revolute(Vec3::new(1.0,0.0,0.0),None).unwrap(),inertia),
    ]).unwrap();
    forward_dynamics(&model,BaseState::stationary(Se3::identity()),&[0.0],&[0.0],&[1.0],
        Vec3::new(0.0,0.0,0.0),&[Wrench::zero();2]).unwrap().generalized_acceleration[0]
}

#[test]
fn material_only_alpha_swap_changes_actual_rotor_response_with_identical_reference_mass() {
    let a = thermal_card("A",1e-4);
    let b = thermal_card("B",2e-4);
    let reference_a = body(&a);
    let reference_b = body(&b);
    assert_eq!(reference_a.inertia(),reference_b.inertia());
    let path_a = path(&a,300.0,400.0);
    let path_b = path(&b,300.0,400.0);
    let hot_a = expand_body(&reference_a,&[binding(&path_a)],||false).unwrap();
    let hot_b = expand_body(&reference_b,&[binding(&path_b)],||false).unwrap();
    let s_a = 1.0+path_a.free_linear_strain();
    let s_b = 1.0+path_b.free_linear_strain();
    close(rotor_acceleration(hot_a.inertia())/rotor_acceleration(hot_b.inertia()),(s_b/s_a).powi(2));
    assert_eq!(hot_a.inertia().mass(),hot_b.inertia().mass());
}

#[test]
fn repeated_configuration_uses_the_reference_not_incremental_strain() {
    let card = thermal_card("repeat",1e-4);
    let reference = body(&card);
    let receipt = path(&card,300.0,400.0);
    let a = expand_body(&reference,&[binding(&receipt)],||false).unwrap();
    let b = expand_body(&reference,&[binding(&receipt)],||false).unwrap();
    assert_eq!(a,b);
    let unchanged = path(&card,300.0,300.0);
    let reset = expand_body(&reference,&[binding(&unchanged)],||false).unwrap();
    assert_eq!(reset.inertia(),reference.inertia());
    assert_eq!(reset.regions()[0].geometry(),reference.regions()[0].geometry());
    let cooling = path(&card,300.0,200.0);
    let cold = expand_body(&reference,&[binding(&cooling)],||false).unwrap();
    assert!(cold.regions()[0].geometry().volume_m3() < reference.regions()[0].geometry().volume_m3());
    assert_eq!(cold.inertia().mass(),reference.inertia().mass());
}

#[test]
fn explicit_fixed_point_changes_centroid_not_centroidal_inertia() {
    let card = thermal_card("pivot",1e-4);
    let reference = body(&card);
    let receipt = path(&card,300.0,400.0);
    let mut b = binding(&receipt);
    b.fixed_point_m = Vec3::new(10.0,0.0,0.0);
    let moved = expand_body(&reference,&[b],||false).unwrap();
    close(moved.inertia().center_of_mass().x,-0.3);
    close(moved.inertia().center_of_mass().y,0.0);
    let centered = expand_body(&reference,&[binding(&receipt)],||false).unwrap();
    assert_eq!(moved.inertia().inertia_com(),centered.inertia().inertia_com());
    assert_eq!(moved.regions()[0].fixed_point_m(),b.fixed_point_m);
}

#[test]
fn foreign_card_and_wrong_reference_state_are_refused() {
    let card = thermal_card("source",1e-4);
    let reference = body(&card);
    // Same material name but changed content must not rebind the mass inventory.
    let foreign = thermal_card("source",2e-4);
    for wrong in [path(&foreign,300.0,400.0),path(&card,290.0,400.0)] {
        assert!(matches!(expand_body(&reference,&[binding(&wrong)],||false),
            Err(MaterialInertiaError::Invalid(_))));
    }
}

#[test]
fn missing_duplicate_nonfinite_and_nonpositive_configurations_fail_closed() {
    let card = thermal_card("bad",1e-4);
    let reference = body(&card);
    let receipt = path(&card,300.0,400.0);
    assert!(expand_body(&reference,&[],||false).is_err());
    assert!(expand_body(&reference,&[binding(&receipt),binding(&receipt)],||false).is_err());
    let mut unknown = binding(&receipt); unknown.key = "missing";
    assert!(expand_body(&reference,&[unknown],||false).is_err());
    let mut bad_point = binding(&receipt); bad_point.fixed_point_m.x = f64::NAN;
    assert!(expand_body(&reference,&[bad_point],||false).is_err());
    let contracting = thermal_card("inverted",-0.01);
    let baseline = body(&contracting);
    let inverted = path(&contracting,300.0,400.0);
    assert!(expand_body(&baseline,&[binding(&inverted)],||false).is_err());
    assert!(matches!(expand_body(&reference,&[binding(&receipt)],||true),Err(MaterialInertiaError::Cancelled)));
    let mut calls = 0;
    assert!(matches!(expand_body(&reference,&[binding(&receipt)],||{ calls+=1; calls==3 }),
        Err(MaterialInertiaError::Cancelled)));
}

#[test]
fn composite_expansion_reuses_parallel_axis_owner_and_is_binding_order_independent() {
    let card = thermal_card("composite",1e-4);
    let cube = VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap();
    let left = VolumeMoments::new(cube.volume_m3(),Vec3::new(-2.0,0.0,0.0),cube.inertia_m5()).unwrap();
    let right = VolumeMoments::new(cube.volume_m3(),Vec3::new(2.0,0.0,0.0),cube.inertia_m5()).unwrap();
    let p = point(300.0);
    let reference = resolve_body(&[
        RegionBinding { key: "left", geometry: left, material: &card, point: &p,
            selection: MaterialPropertySelection::SingleClaimOnly },
        RegionBinding { key: "right", geometry: right, material: &card, point: &p,
            selection: MaterialPropertySelection::SingleClaimOnly },
    ],2,||false).unwrap();
    let receipt = path(&card,300.0,400.0);
    let mut a = binding(&receipt); a.key = "left";
    let mut b = binding(&receipt); b.key = "right";
    let expanded = expand_body(&reference,&[a,b],||false).unwrap();
    let reordered = expand_body(&reference,&[b,a],||false).unwrap();
    assert_eq!(expanded,reordered);
    assert_eq!(expanded.inertia().mass(),reference.inertia().mass());
    matrix_close(expanded.inertia().inertia_com(),reference.inertia().inertia_com(),1.03_f64.powi(2));
    assert!(expand_body(&reference,&[a,a],||false).is_err());
}

#[test]
fn complete_mechanism_thermal_update_preserves_topology_and_drives_existing_solver() {
    use fs_mbd_material::articulated::{MaterialLink,resolve_articulated};
    use fs_mbd_material::thermal::{LinkThermalExpansion,expand_articulated};
    let card = thermal_card("mechanism",1e-4);
    let p = point(300.0);
    let regions = [RegionBinding { key: "body", material: &card, point: &p,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0,2.0,3.0)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly }];
    let reference = resolve_articulated(&[
        MaterialLink { name: "mount", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &regions },
        MaterialLink { name: "rotor", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: JointModel::revolute(Vec3::new(1.0,0.0,0.0),None).unwrap(),regions: &regions },
    ],2,||false).unwrap();
    let original = reference.clone();
    let cold = path(&card,300.0,300.0);
    let hot = path(&card,300.0,400.0);
    let fixed_regions = [binding(&cold)];
    let rotor_regions = [binding(&hot)];
    let fixed = LinkThermalExpansion { link: "mount", regions: &fixed_regions };
    let rotor = LinkThermalExpansion { link: "rotor", regions: &rotor_regions };
    let configured = expand_articulated(&reference,&[rotor,fixed],||false).unwrap();
    let reordered = expand_articulated(&reference,&[fixed,rotor],||false).unwrap();
    assert_eq!(configured,reordered);
    assert_eq!(reference,original);
    assert!(std::ptr::eq(configured.reference(),&reference));
    for (a,b) in configured.model().links().iter().zip(reference.model().links()) {
        assert_eq!(a.name(),b.name()); assert_eq!(a.parent(),b.parent());
        assert_eq!(a.joint(),b.joint()); assert_eq!(a.parent_from_child_zero(),b.parent_from_child_zero());
        assert_eq!(a.inertia().mass(),b.inertia().mass());
    }
    assert_eq!(configured.bodies()[0].inertia(),reference.bodies()[0].inertia());
    let acceleration = |model: &ArticulatedModel| forward_dynamics(model,
        BaseState::stationary(Se3::identity()),&[0.0],&[0.0],&[1.0],
        Vec3::new(0.0,0.0,0.0),&[Wrench::zero();2]).unwrap().generalized_acceleration[0];
    close(acceleration(reference.model())/acceleration(configured.model()),1.03_f64.powi(2));
    assert!(expand_articulated(&reference,&[fixed],||false).is_err());
    assert!(expand_articulated(&reference,&[fixed,fixed],||false).is_err());
    let foreign = LinkThermalExpansion { link: "foreign", regions: &rotor_regions };
    assert!(expand_articulated(&reference,&[fixed,foreign],||false).is_err());
    // Discover and then cancel at the final checkpoint, after all body work.
    let mut checkpoints = 0;
    expand_articulated(&reference,&[fixed,rotor],||{ checkpoints+=1; false }).unwrap();
    let mut calls = 0;
    assert!(matches!(expand_articulated(&reference,&[fixed,rotor],||{
        calls+=1; calls==checkpoints
    }),Err(MaterialInertiaError::Cancelled)));
    assert_eq!(reference,original);
}
