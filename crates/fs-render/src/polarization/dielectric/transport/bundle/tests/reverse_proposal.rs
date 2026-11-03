//! Reverse proposals remain separate from physical polarized transmission.

use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoherentDielectricCoating, CoherentPeriodicCoating,
    ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{DielectricCoating, DielectricSurface};

fn provenance(name: &[u8]) -> ConductorSource {
    ConductorSource::try_new(
        fs_blake3::hash_domain("native-bundle-reverse-proposal", name),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

#[test]
fn g1_reverse_ball_proposal_is_scalar_product_not_polarized_power() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.35, 0.0, -3.0),
            dir: unit(Vec3::new(0.05, 0.0, 1.0)),
        };
        let lens = glass(1.5, 0.17);
        let reference = reference_ball(ray, 1.5, 2.5);
        let path = through_ball(cx, ray, &lens, 550.0, 2.5);
        let mean = 0.5 * (reference.transmission_s + reference.transmission_p);
        let probability = mean * mean;
        let physical = 0.5 *
            (reference.transmission_s.powi(2) + reference.transmission_p.powi(2));
        close(path.reverse_branch_log_probability().exp(), probability, 3e-12);
        close(path.interface_log_power().exp(), physical, 3e-12);
        assert!((physical - probability).abs() > 1e-6);
        assert!(path.propagation().propagation().optical_depth > 0.0);
        // Ordinary travel in the original material preserves the proposal.
        let air = glass(1.0, 0.0);
        let continued = path.continue_in(&HomogeneousIndexField::new(&air),
            0.25, BundleOptions::default()).unwrap();
        assert_eq!(continued.reverse_branch_log_probability(),
            path.reverse_branch_log_probability());
    });
}

#[test]
fn g3_absorbing_asymmetric_stack_uses_native_reverse_survival_proposal() {
    let layers = [(1.25, 0.0, 70e-9), (2.15, 0.35, 45e-9), (1.4, 0.0, 115e-9)]
        .map(|(index, extinction, thickness)| {
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(index).unwrap(),
                extinction, thickness).unwrap()
        });
    let coated = DielectricSurface::SMOOTH.with_coating(
        DielectricCoating::try_new(CoherentPeriodicCoating::try_new(
            &layers, 1, provenance(b"asymmetric-absorbing-stack"),
        ).unwrap()).unwrap(),
    );
    let incident = glass(1.0, 0.0);
    let transmitted = glass(1.7, 0.0);
    let cosine = 0.8_f64;
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let ray = Ray { origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -cosine) };
    let first = propagate_bundle(&HomogeneousIndexField::new(&incident), ray,
        Vec3::new(0.0, -1.0, 0.0), 550.0, 0.0, BundleOptions::default()).unwrap();
    let original = PolarizedRayBundle::unpolarized(first).unwrap();
    let target = HomogeneousIndexField::new(&transmitted).sample(ray.origin, 550.0).unwrap();
    let forward = coated.interface(550.0, false).unwrap()
        .fresnel_polarized(cosine, 1.0, 1.7).unwrap();
    let forward_q = forward.scalar.event_probabilities();
    assert!(forward.scalar.absorption_residual() > 0.01);
    // Independent complex-admittance Airy recursion, with normal-power
    // normalization at the two exterior media, gives these scalar proposals.
    close(forward_q[1], 0.733_728_280_909_517_9, 3e-12);
    for event in [DielectricEvent::Transmission, DielectricEvent::Reflection] {
        let result = original.smooth_surface_interface(plane_normal(normal), target,
            &coated, event, InterfaceOptions::default()).unwrap();
        let reverse_direction = result.propagation().propagation().ray.dir.scale(-1.0);
        let (from, to) = if event == DielectricEvent::Transmission {
            (transmitted, incident)
        } else { (incident, transmitted) };
        let reverse_beam = PolarizedBeam::try_new(ray.origin,
            PolarizationFrame::try_new(reverse_direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
            550.0, from, Stokes::unpolarized(1.0)).unwrap();
        let native = [0.0, 0.5, 1.0 - f64::EPSILON].into_iter()
            .map(|u| reverse_beam.sample_interface(normal, to, coated, [0.2, 0.7, u])
                .unwrap().unwrap())
            .find(|sample| sample.event == event).unwrap();
        close(result.reverse_branch_log_probability().exp(), native.pdf, 3e-12);
        let index = usize::from(event == DielectricEvent::Transmission);
        if event == DielectricEvent::Transmission {
            close(native.pdf, 0.693_698_153_925_260_4, 3e-12);
            assert!((native.pdf - forward_q[index]).abs() > 1e-4,
                "absorbing asymmetric stack must retain different reverse proposal");
            close(result.interface_log_power().exp(), forward.transmission.0[0][0], 3e-12);
            let wrong_order = coated.interface(550.0, false).unwrap()
                .fresnel(reverse_direction.z, 1.7, 1.0).unwrap()
                .event_probabilities()[1];
            assert!((native.pdf - wrong_order).abs() > 1e-4);
        } else {
            close(native.pdf, forward_q[index], 3e-12);
            close(result.interface_log_power().exp(), forward.reflection.0[0][0], 3e-12);
        }
    }
    assert_eq!(original.reverse_branch_log_probability(), 0.0);
}

#[test]
fn g1_transparent_film_and_deterministic_mirror_keep_distinct_proposals() {
    let medium = glass(1.0, 0.0);
    let field = HomogeneousIndexField::new(&medium);
    let ray = Ray { origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0) };
    let first = propagate_bundle(&field, ray, Vec3::new(1.0, 0.0, 0.0),
        550.0, 1.0, BundleOptions::default()).unwrap();
    let original = PolarizedRayBundle::unpolarized(first).unwrap();
    let endpoint = original.propagation().propagation().ray.origin;
    let target = field.sample(endpoint, 550.0).unwrap();
    let normal = plane_normal(Vec3::new(0.0, 0.0, -1.0));
    let film = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.5).unwrap(), 550e-9 / 6.0,
        provenance(b"quarter-wave-film"),
    ).unwrap();
    let reflected = original.smooth_coated_interface(normal, target, &film,
        DielectricEvent::Reflection, InterfaceOptions::default()).unwrap();
    let reflectance = ((1.0 - 2.25_f64) / (1.0 + 2.25)).powi(2);
    close(reflected.reverse_branch_log_probability().exp(), reflectance, 3e-12);
    let mirror = reflected.smooth_conductor_interface(
        plane_normal(Vec3::new(0.0, 0.0, 1.0)),
        &ConductorOptics::representative_tungsten(), InterfaceOptions::default(),
    ).unwrap();
    assert_eq!(mirror.reverse_branch_log_probability(),
        reflected.reverse_branch_log_probability());
    assert!(mirror.interface_log_power() < reflected.interface_log_power());
}
