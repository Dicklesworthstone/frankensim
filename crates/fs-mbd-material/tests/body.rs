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

/// Volume, Ixx and Izz of the prism of straight tetrahedra on the n-gon
/// inscribed in a cylinder of radius `big_r` and length `l`, each as its
/// relative shortfall from the cylinder's closed form.
fn inscribed_prism_shortfall(n: u32, big_r: f64, l: f64) -> [f64; 3] {
    let at = |k: u32, z: f64| {
        let angle = 2.0 * core::f64::consts::PI * f64::from(k) / f64::from(n);
        Vec3::new(big_r * angle.cos(), big_r * angle.sin(), z)
    };
    let (h, mut volume, mut inertia) = (l / 2.0, 0.0, [0.0; 9]);
    for k in 0..n {
        let [a, b, c] = [Vec3::new(0.0, 0.0, -h), at(k, -h), at(k + 1, -h)];
        let [d, e, f] = [Vec3::new(0.0, 0.0, h), at(k, h), at(k + 1, h)];
        for tet in [[a, b, c, d], [b, c, d, e], [c, d, e, f]] {
            let part = VolumeMoments::tetrahedron(tet).unwrap();
            let (v, g) = (part.volume_m3(), part.centroid_m());
            // Parallel axis to the cylinder's centre, which is the prism's.
            let (xy, xz, yz) = (-g.x * g.y, -g.x * g.z, -g.y * g.z);
            let [xx, yy, zz] = [g.x * g.x, g.y * g.y, g.z * g.z];
            let shift = [yy + zz, xy, xz, xy, xx + zz, yz, xz, yz, xx + yy];
            volume += v;
            for ((sum, own), moved) in inertia.iter_mut().zip(part.inertia_m5().m).zip(shift) {
                *sum += own + v * moved;
            }
        }
    }
    let exact = VolumeMoments::cylinder(big_r, l).unwrap();
    let moments = exact.inertia_m5().m;
    let exact = [exact.volume_m3(), moments[0], moments[8]];
    let mesh = [volume, inertia[0], inertia[8]];
    core::array::from_fn(|i| (exact[i] - mesh[i]) / exact[i])
}

#[test]
fn cylinder_moments_match_the_closed_form_the_tube_difference_and_a_converging_mesh() {
    let (r, big_r, l) = (0.3, 0.5, 0.7);
    let solid = VolumeMoments::cylinder(big_r, l).unwrap();
    let v = core::f64::consts::PI * big_r * big_r * l;
    close(solid.volume_m3(), v);
    assert_eq!(solid.centroid_m(), Vec3::new(0.0, 0.0, 0.0));
    let moments = solid.inertia_m5().m;
    close(moments[8], v * big_r * big_r / 2.0);
    for i in [0, 4] {
        close(moments[i], v * (3.0 * big_r * big_r + l * l) / 12.0);
    }
    for i in [1, 2, 3, 5, 6, 7] {
        assert!(moments[i].abs() <= 0.0);
    }
    // A tube is the outer cylinder less the inner one; both are centred
    // there, so their volumes and centroidal moments subtract.
    let tube = VolumeMoments::hollow_cylinder(r, big_r, l).unwrap();
    let core = VolumeMoments::cylinder(r, l).unwrap();
    close(tube.volume_m3(), solid.volume_m3() - core.volume_m3());
    let pairs = moments.into_iter().zip(core.inertia_m5().m);
    for (got, (outer, inner)) in tube.inertia_m5().m.into_iter().zip(pairs) {
        close(got, outer - inner);
    }
    let same = VolumeMoments::hollow_cylinder(0.0, big_r, l).unwrap();
    assert_eq!(same, solid);
    // The inscribed prism of straight tetrahedra falls short of the
    // cylinder at second order in 1/n, in each of the three quantities.
    let coarse = inscribed_prism_shortfall(32, big_r, l);
    let fine = inscribed_prism_shortfall(64, big_r, l);
    for (c, f) in coarse.into_iter().zip(fine) {
        let order = (c / f).log2();
        assert!(c > 0.0 && f > 0.0 && (order - 2.0).abs() < 0.02, "{c} {f}");
    }
    for (inner, outer, length) in [
        (0.5, 0.5, 1.0),
        (-0.1, 0.5, 1.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.5, 0.0),
        (f64::NAN, 0.5, 1.0),
        (0.0, f64::INFINITY, 1.0),
    ] {
        assert!(VolumeMoments::hollow_cylinder(inner, outer, length).is_err());
    }
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
