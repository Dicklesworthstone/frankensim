mod common;
use common::*;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{MaterialInertiaError, RegionBinding, VolumeMoments, resolve_body};
use fs_qty::{Density, Pressure};

fn region<'a>(key: &'a str, material: &'a MaterialCard, point: &'a QueryPoint,
    geometry: VolumeMoments) -> RegionBinding<'a> {
    RegionBinding { key, material, point, geometry,
        selection: MaterialPropertySelection::SingleClaimOnly }
}
fn cube_at(center: Vec3) -> VolumeMoments {
    VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap()
        .transformed(Se3::from_parts(So3::identity(), center).unwrap()).unwrap()
}

#[test]
fn density_only_material_changes_mass_inertia_momentum_and_energy_together() {
    let point = point(300.0);
    let light = density_card("material-a", 2.0, 2.0);
    let heavy = density_card("material-b", 6.0, 6.0);
    let geometry = cube_at(Vec3::new(0.4, -0.3, 0.2));
    let a = resolve_body(&[region("body", &light, &point, geometry)], 1, || false).unwrap();
    let b = resolve_body(&[region("body", &heavy, &point, geometry)], 1, || false).unwrap();
    close(a.inertia().mass(), 2.0);
    close(b.inertia().mass(), 6.0);
    assert_eq!(a.inertia().center_of_mass(), b.inertia().center_of_mass());
    for (x, y) in a.inertia().matrix().unwrap().m.iter().zip(b.inertia().matrix().unwrap().m) {
        close(3.0*x, y);
    }
    let velocity = Twist::new(Vec3::new(0.2, -0.4, 0.8), Vec3::new(1.0, 2.0, 3.0));
    close(3.0*a.inertia().kinetic_energy(velocity).unwrap(), b.inertia().kinetic_energy(velocity).unwrap());
    assert_ne!(a.regions()[0].material(), b.regions()[0].material());
    close(b.regions()[0].material().property("density").unwrap().value_si(), 6.0);
    // Neither card contains elastic, yield or inertial overrides.
}

#[test]
fn asymmetric_material_assignment_moves_balance_and_full_tensor() {
    let point = point(300.0);
    let light = density_card("light", 1.0, 1.0);
    let heavy = density_card("heavy", 3.0, 3.0);
    let geometry_a = cube_at(Vec3::new(-2.0, -1.0, 0.0));
    let geometry_b = cube_at(Vec3::new(2.0, 1.0, 0.0));
    let bindings = [region("left", &light, &point, geometry_a), region("right", &heavy, &point, geometry_b)];
    let body = resolve_body(&bindings, 2, || false).unwrap();
    let inertia = body.inertia();
    close(inertia.mass(), 4.0);
    close(inertia.center_of_mass().x, 1.0);
    close(inertia.center_of_mass().y, 0.5);
    let tensor = inertia.inertia_com().m;
    for (i, expected) in [(0, 2.0/3.0+3.0), (4, 2.0/3.0+12.0), (8, 2.0/3.0+15.0), (1, -6.0), (3, -6.0)] {
        close(tensor[i], expected);
    }
    let mut reversed = bindings.clone(); reversed.reverse();
    assert_eq!(body, resolve_body(&reversed, 2, || false).unwrap());
    let exchanged = resolve_body(&[region("left", &heavy, &point, geometry_a),
        region("right", &light, &point, geometry_b)], 2, || false).unwrap();
    close(exchanged.inertia().center_of_mass().x, -1.0);
    close(exchanged.inertia().center_of_mass().y, -0.5);
}

#[test]
fn temperature_query_re_resolves_inertia_instead_of_using_a_cached_mass() {
    let material = density_card("temperature-law", 2.0, 4.0);
    let geometry = cube_at(Vec3::new(0.0, 0.0, 0.0));
    let cold = resolve_body(&[region("body", &material, &point(250.0), geometry)], 1, || false).unwrap();
    let hot = resolve_body(&[region("body", &material, &point(350.0), geometry)], 1, || false).unwrap();
    close(cold.inertia().mass(), 2.5); close(hot.inertia().mass(), 3.5);
    close(hot.inertia().inertia_com().m[0]/cold.inertia().inertia_com().m[0], 1.4);
    assert_ne!(cold.regions()[0].material(), hot.regions()[0].material());
}

#[test]
fn tetrahedral_geometry_matches_independent_simplex_integrals() {
    let vertices = [Vec3::new(0.0,0.0,0.0), Vec3::new(1.0,0.0,0.0),
        Vec3::new(0.0,1.0,0.0), Vec3::new(0.0,0.0,1.0)];
    let geometry = VolumeMoments::tetrahedron(vertices).unwrap();
    close(geometry.volume_m3(), 1.0/6.0);
    assert_eq!(geometry.centroid_m(), Vec3::new(0.25,0.25,0.25));
    for i in [0,4,8] { close(geometry.inertia_m5().m[i], 1.0/80.0); }
    for i in [1,2,3,5,6,7] { close(geometry.inertia_m5().m[i], 1.0/480.0); }
    let mut reversed = vertices; reversed.swap(1,2);
    let reordered = VolumeMoments::tetrahedron(reversed).unwrap();
    close(geometry.volume_m3(), reordered.volume_m3());
    for (a,b) in geometry.inertia_m5().m.iter().zip(reordered.inertia_m5().m) { close(*a,b); }
    assert!(VolumeMoments::tetrahedron([vertices[0];4]).is_err());
}

#[test]
fn placement_rotates_inertia_without_an_extra_centroidal_parallel_axis_term() {
    let geometry = VolumeMoments::cuboid(Vec3::new(1.0, 2.0, 3.0)).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.0,0.0,core::f64::consts::FRAC_PI_2))).unwrap();
    let transformed = geometry.transformed(Se3::from_parts(rotation, Vec3::new(20.0,30.0,40.0)).unwrap()).unwrap();
    close(transformed.volume_m3(), geometry.volume_m3());
    close(transformed.centroid_m().x, 20.0);
    close(transformed.centroid_m().y, 30.0);
    close(transformed.inertia_m5().m[0], geometry.inertia_m5().m[4]);
    close(transformed.inertia_m5().m[4], geometry.inertia_m5().m[0]);
    close(transformed.inertia_m5().m[8], geometry.inertia_m5().m[8]);
}

#[test]
fn central_merge_is_not_corrupted_by_a_large_common_coordinate_origin() {
    let material = density_card("uniform", 3.0, 3.0); let p = point(300.0);
    let make = |shift: f64| resolve_body(&[
        region("a", &material, &p, cube_at(Vec3::new(shift-2.0,0.0,0.0))),
        region("b", &material, &p, cube_at(Vec3::new(shift+2.0,0.0,0.0)))], 2, || false).unwrap();
    let a = make(0.0); let b = make(1e8);
    close(b.inertia().center_of_mass().x, 1e8);
    for (x,y) in a.inertia().inertia_com().m.iter().zip(b.inertia().inertia_com().m) { close(*x,y); }
}

#[test]
fn source_failures_do_not_silently_choose_density_or_extrapolate() {
    let geometry = cube_at(Vec3::new(0.0,0.0,0.0));
    let missing = card("missing", "not-density", Density::DIMS, &[(2.0,2.0)]);
    let wrong = card("wrong-dimension", "density", Pressure::DIMS, &[(2.0,2.0)]);
    let conflict = card("conflicting", "density", Density::DIMS, &[(2.0,2.0),(3.0,3.0)]);
    for material in [&missing,&wrong,&conflict] {
        assert!(matches!(resolve_body(&[region("body", material, &point(300.0), geometry)],
            1, || false), Err(MaterialInertiaError::Material(_))));
    }
    let good = density_card("good", 2.0,2.0);
    assert!(matches!(resolve_body(&[region("body", &good, &point(450.0), geometry)], 1, || false),
        Err(MaterialInertiaError::Material(_))));
}

#[test]
fn count_identity_geometry_and_late_cancellation_fail_before_publication() {
    let material = density_card("material", 1.0,1.0); let p = point(300.0);
    let geometry = cube_at(Vec3::new(0.0,0.0,0.0));
    let binding = region("body", &material, &p, geometry);
    assert!(resolve_body(&[], 1, || false).is_err());
    assert!(resolve_body(&[binding.clone()], 0, || false).is_err());
    assert!(resolve_body(&[binding.clone()], 4097, || false).is_err());
    assert!(resolve_body(&[binding.clone(),binding.clone()], 2, || false).is_err());
    let mut polls = 0;
    assert_eq!(resolve_body(&[binding], 1, || { polls += 1; polls == 3 }), Err(MaterialInertiaError::Cancelled));
    assert!(VolumeMoments::cuboid(Vec3::new(0.0,1.0,1.0)).is_err());
    assert!(VolumeMoments::new(1.0, Vec3::new(0.0,0.0,0.0), Mat3 {m:[0.0;9]}).is_err());
}
