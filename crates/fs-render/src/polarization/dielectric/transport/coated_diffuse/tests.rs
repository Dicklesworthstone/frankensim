use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use fs_geom::Point3;

fn beam(direction: Vec3, state: Stokes) -> PolarizedBeam {
    beam_in_medium(direction, state, 1.0)
}

fn beam_in_medium(direction: Vec3, state: Stokes, index: f64) -> PolarizedBeam {
    let glass = DielectricGlass::new(CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR, GlassProvenance::Custom);
    let axis = if direction.y == 0.0 { Vec3::new(0.0, 1.0, 0.0) } else { Vec3::new(1.0, 0.0, 0.0) };
    PolarizedBeam::try_new(Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, axis).unwrap(), 550.0, glass, state).unwrap()
}

#[test]
fn g1_absorbing_coats_share_scalar_power_pdfs_and_connections_in_either_medium() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for (outside, inside) in [(1.0, 1.0), (1.0, 1.5), (1.5, 1.0)] {
        let coating = SmoothDiffuseCoating::try_new(outside, inside).unwrap()
            .with_absorption(0.001, 700.0).unwrap();
        let input = beam_in_medium(Vec3::new(0.6, 0.0, -0.8), Stokes::unpolarized(2.0), outside);
        let masses = coating.branch_probabilities(0.8).unwrap();
        let draw = if masses[0] <= masses[1] { masses[0]+0.5*masses[1] } else { 0.5*masses[1] };
        for unit in [0.0, 0.3, 0.8] {
            let event = input.sample_coated_diffuse(normal, coating, 0.7, [unit, 0.2, draw]).unwrap();
            assert!(!event.delta);
            let direction = event.beam.frame().direction();
            close(event.beam.power().i*event.pdf/(2.0*direction.z),
                coating.diffuse_brdf(0.8, direction.z, 0.7).unwrap());
            close(event.pdf, masses[1]*coating.diffuse_direction_pdf(direction.z).unwrap());
            let source = PolarizationFrame::try_new(direction.scale(-1.0), event.beam.frame().horizontal()).unwrap();
            let response = input.evaluate_coated_diffuse_connection(normal, coating, 0.7,
                source, Stokes::unpolarized(1.0)).unwrap();
            close(response.factor, coating.diffuse_brdf(0.8, direction.z, 0.7).unwrap());
            close(response.pdf, event.pdf);
        }
        let a = beam_in_medium(Vec3::new(0.6, 0.0, -0.8),
            Stokes { i: 1.0, q: 0.2, u: 0.3, v: 0.4 }, outside);
        let b = beam_in_medium(Vec3::new(0.0, 0.6, -0.8),
            Stokes { i: 1.0, q: -0.3, u: 0.1, v: 0.2 }, outside);
        let forward = a.evaluate_coated_diffuse_connection(normal, coating, 0.7, b.frame(), b.power()).unwrap();
        let reverse = b.evaluate_coated_diffuse_connection(normal, coating, 0.7, a.frame(), a.power()).unwrap();
        close(forward.factor, reverse.factor);
    }
}

#[test]
fn g1_absorbing_coat_furnace_retains_polarized_entry_and_repeated_layer_loss() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap()
        .with_absorption(0.001, 400.0).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let fresnel = fresnel_dielectric_polarized(0.8, 1.0, 1.5).unwrap();
    let masses = coating.branch_probabilities(0.8).unwrap();
    let rho = 0.7;
    let denominator = (1.0-rho)+rho*(coating.escape_probability()+coating.layer_absorption_probability());
    let downward = (-coating.normal_optical_depth()/fresnel.scalar.transmitted_cosine).exp();
    for q in [-1.0, 0.0, 1.0] {
        let input = beam(Vec3::new(0.6, 0.0, -0.8), Stokes { i: 1.0, q, u: 0.0, v: 0.0 });
        let mirror = input.sample_coated_diffuse(normal, coating, rho, [0.3, 0.2, 0.0]).unwrap();
        let mut diffuse = 0.0;
        const N: usize = 8192;
        for i in 0..N {
            let event = input.sample_coated_diffuse(normal, coating, rho,
                [(i as f64+0.5)/N as f64, 0.2, 0.9]).unwrap();
            diffuse += event.beam.power().i/N as f64;
        }
        let entered = fresnel
            .transmission
            .apply_diagonal_jones(input.power())
            .unwrap()
            .i;
        let expected = fresnel
            .reflection
            .apply_diagonal_jones(input.power())
            .unwrap()
            .i
            + entered*downward*rho*coating.escape_probability()/denominator;
        assert!((masses[0]*mirror.beam.power().i+masses[1]*diffuse-expected).abs() < 4e-6);
    }
}

#[test]
fn g1_immersed_coat_scalar_agreement_and_external_tir_support() {
    let coating = SmoothDiffuseCoating::try_new(1.5, 1.0).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let input = beam_in_medium(Vec3::new(0.6, 0.0, -0.8), Stokes::unpolarized(2.0), 1.5);
    let masses = coating.branch_probabilities(0.8).unwrap();
    let diffuse_draw = if masses[0] <= masses[1] { masses[0]+0.5*masses[1] } else { 0.5*masses[1] };
    for unit in [0.0, 0.2, 0.8, 0.999] {
        let event = input.sample_coated_diffuse(normal, coating, 0.7, [unit, 0.3, diffuse_draw]).unwrap();
        assert!(!event.delta);
        let direction = event.beam.frame().direction();
        assert!(direction.z > (1.0_f64 - (1.0_f64/1.5).powi(2)).sqrt());
        let density = masses[1]*coating.diffuse_direction_pdf(direction.z).unwrap();
        close(event.pdf, density);
        close(event.beam.power().i*event.pdf/(2.0*direction.z),
            coating.diffuse_brdf(0.8, direction.z, 0.7).unwrap());
        let source = PolarizationFrame::try_new(direction.scale(-1.0), event.beam.frame().horizontal()).unwrap();
        let response = input.evaluate_coated_diffuse_connection(normal, coating, 0.7,
            source, Stokes::unpolarized(1.0)).unwrap();
        close(response.factor, coating.diffuse_brdf(0.8, direction.z, 0.7).unwrap());
        close(response.pdf, event.pdf);
    }
    let trapped = beam_in_medium(Vec3::new(3.0_f64.sqrt()/2.0, 0.0, -0.5),
        Stokes { i: 1.0, q: 0.0, u: 0.6, v: 0.8 }, 1.5);
    for draw in [0.0, 0.5, 0.999] {
        let mirror = trapped.sample_coated_diffuse(normal, coating, 0.0, [0.2, 0.3, draw]).unwrap();
        assert!(mirror.delta);
        assert_eq!(mirror.pdf, 1.0);
        close(mirror.beam.power().i, 1.0);
        close(mirror.beam.power().degree(), 1.0);
    }
    let blocked_source = input.evaluate_coated_diffuse_connection(normal, coating, 1.0,
        trapped.frame(), Stokes::unpolarized(1.0)).unwrap();
    assert_eq!(blocked_source.factor, 0.0);
    assert_eq!(blocked_source.pdf, 0.0);
    let blocked_camera = trapped.evaluate_coated_diffuse_connection(normal, coating, 1.0,
        input.frame(), Stokes::unpolarized(1.0)).unwrap();
    assert_eq!(blocked_camera.factor, 0.0);
    assert_eq!(blocked_camera.pdf, 0.0);
}

#[test]
fn g1_immersed_coat_white_furnace_retains_polarized_incident_power() {
    let coating = SmoothDiffuseCoating::try_new(1.5, 1.0).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for q in [-1.0, 0.0, 1.0] {
        let input = beam_in_medium(Vec3::new(0.6, 0.0, -0.8),
            Stokes { i: 1.0, q, u: 0.0, v: 0.0 }, 1.5);
        let masses = coating.branch_probabilities(0.8).unwrap();
        let mirror_draw = if masses[0] <= masses[1] { 0.5*masses[0] } else { masses[1]+0.5*masses[0] };
        let diffuse_draw = if masses[0] <= masses[1] { masses[0]+0.5*masses[1] } else { 0.5*masses[1] };
        let mirror = input.sample_coated_diffuse(normal, coating, 1.0, [0.3, 0.2, mirror_draw]).unwrap();
        let mut diffuse = 0.0;
        const N: usize = 16384;
        for i in 0..N {
            let event = input.sample_coated_diffuse(normal, coating, 1.0,
                [(i as f64+0.5)/N as f64, 0.2, diffuse_draw]).unwrap();
            diffuse += event.beam.power().i/N as f64;
        }
        assert!((masses[0]*mirror.beam.power().i + masses[1]*diffuse - 1.0).abs() < 4e-6);
    }
}

fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-11*a.abs().max(b.abs()).max(1.0), "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_coated_diffuse_matches_scalar_response_and_mixed_proposal() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let input = beam(Vec3::new(0.6, 0.0, -0.8), Stokes::unpolarized(2.0));
    let sample = input.sample_coated_diffuse(normal, coating, 0.6, [0.36, 0.27, 0.9]).unwrap();
    assert!(!sample.delta);
    let outgoing = sample.beam.frame().direction();
    let cosine = outgoing.z;
    close(sample.pdf, coating.branch_probabilities(0.8).unwrap()[1]*cosine/core::f64::consts::PI);
    close(sample.beam.power().i*sample.pdf/(2.0*cosine), coating.diffuse_brdf(0.8, cosine, 0.6).unwrap());
    assert_eq!(sample.beam.medium(), input.medium());
    assert_eq!(sample.beam.ray().origin, input.ray().origin);
    let source_frame = PolarizationFrame::try_new(outgoing.scale(-1.0), sample.beam.frame().horizontal()).unwrap();
    let response = input.evaluate_coated_diffuse_connection(normal, coating, 0.6,
        source_frame, Stokes::unpolarized(1.0)).unwrap();
    close(response.factor, coating.diffuse_brdf(0.8, cosine, 0.6).unwrap());
    close(response.pdf, sample.pdf);
    let mirror = input.sample_coated_diffuse(normal, coating, 0.6, [0.36, 0.27, 0.0]).unwrap();
    assert!(mirror.delta);
    close(mirror.pdf, coating.branch_probabilities(0.8).unwrap()[0]);
    close(mirror.beam.power().i*mirror.pdf/2.0,
        crate::dielectric::fresnel_dielectric(0.8, 1.0, 1.5).unwrap().reflectance);
}

#[test]
fn g1_base_depolarizes_but_exit_retains_fresnel_polarization() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let cosine = 1.0 / 1.5_f64.hypot(1.0);
    for q in [-1.0, 0.0, 1.0] {
        let input = beam(Vec3::new(0.0, 0.0, -1.0), Stokes { i: 1.0, q, u: 0.0, v: 0.0 });
        let sample = input.sample_coated_diffuse(normal, coating, 0.7,
            [1.0-cosine*cosine, 0.25, 0.9]).unwrap();
        let state = sample.beam.power();
        let native = fresnel_dielectric_polarized(sample.beam.frame().direction().z, 1.0, 1.5).unwrap();
        close(state.q/state.i, native.transmission.0[0][1]/native.transmission.0[0][0]);
        assert!(state.q < 0.0);
        assert_eq!(state.u, 0.0);
        assert_eq!(state.v, 0.0);
        let normal_exit = input.sample_coated_diffuse(normal, coating, 0.7, [0.0, 0.25, 0.9]).unwrap();
        close(normal_exit.beam.power().degree(), 0.0);
    }
}

#[test]
fn g3_reciprocal_source_and_analyser_retain_full_stokes_frames() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let a = beam(Vec3::new(0.6, 0.0, -0.8), Stokes { i: 1.0, q: 0.2, u: 0.3, v: 0.4 });
    let b = beam(Vec3::new(0.0, 0.8, -0.6), Stokes { i: 1.0, q: -0.3, u: 0.1, v: 0.2 });
    let forward = a.evaluate_coated_diffuse_connection(normal, coating, 0.6, b.frame(), b.power()).unwrap();
    let reverse = b.evaluate_coated_diffuse_connection(normal, coating, 0.6, a.frame(), a.power()).unwrap();
    close(forward.factor, reverse.factor);
    let rotated = a.with_reference_axis(Vec3::new(1.0, 0.0, 0.0)).unwrap();
    close(forward.factor, rotated.evaluate_coated_diffuse_connection(normal.scale(-3.0), coating,
        0.6, b.frame(), b.power()).unwrap().factor);
}

#[test]
fn g1_white_furnace_conserves_polarized_incident_power() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for q in [-1.0, 0.0, 1.0] {
        let input = beam(Vec3::new(0.6, 0.0, -0.8), Stokes { i: 1.0, q, u: 0.0, v: 0.0 });
        let masses = coating.branch_probabilities(0.8).unwrap();
        let mirror = input.sample_coated_diffuse(normal, coating, 1.0, [0.3, 0.2, 0.0]).unwrap();
        let mut diffuse = 0.0;
        const N: usize = 16384;
        for i in 0..N {
            let sample = input.sample_coated_diffuse(normal, coating, 1.0,
                [(i as f64+0.5)/N as f64, 0.2, 0.9]).unwrap();
            diffuse += sample.beam.power().i/N as f64;
        }
        assert!((masses[0]*mirror.beam.power().i + masses[1]*diffuse - 1.0).abs() < 4e-6);
    }
}

#[test]
fn g0_matched_black_medium_and_sample_limits() {
    let input = beam(Vec3::new(0.6, 0.0, -0.8), Stokes { i: 1.0, q: 0.2, u: 0.3, v: 0.4 });
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let matched = SmoothDiffuseCoating::try_new(1.0, 1.0).unwrap();
    let actual = input.sample_coated_diffuse(normal, matched, 0.4, [0.2, 0.6, 0.3]).unwrap();
    let original = input.sample_depolarizing_diffuse(normal, 0.4, [0.2, 0.6]).unwrap();
    assert_eq!(actual.beam, original.beam);
    assert_eq!(actual.pdf, original.pdf_solid_angle);
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    assert_eq!(input.sample_coated_diffuse(normal, coating, 0.0, [0.2, 0.6, 0.9]).unwrap().beam.power().i, 0.0);
    assert_eq!(input.sample_coated_diffuse(normal, SmoothDiffuseCoating::try_new(1.33, 1.5).unwrap(),
        0.4, [0.2, 0.6, 0.3]), Err(BeamError::MediumMismatch));
    for samples in [[-0.1, 0.2, 0.3], [0.1, 1.0, 0.3], [0.1, 0.2, f64::NAN]] {
        assert_eq!(input.sample_coated_diffuse(normal, coating, 0.4, samples), Err(BeamError::InvalidDiffuseSample));
    }
}

#[test]
fn g0_coated_diffuse_normalizes_grazing_transmission_before_application() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let cosine = 1e-310;
    let input = beam(Vec3::new(1.0, 0.0, -cosine), Stokes::unpolarized(1.0));
    let mass = coating.branch_probabilities(cosine).unwrap();
    assert!(mass[1] > 0.0 && mass[1] < f64::MIN_POSITIVE);
    assert!((coating.diffuse_weight(1.0).unwrap() / mass[1]).is_infinite());
    let sample = input.sample_coated_diffuse(normal, coating, 1.0,
        [0.3, 0.2, 0.5 * mass[1]]).unwrap();
    assert!(!sample.delta);
    assert!(sample.pdf > 0.0);
    let exit = crate::dielectric::fresnel_dielectric(
        sample.beam.frame().direction().z, 1.0, 1.5).unwrap();
    let expected = coating.diffuse_weight(1.0).unwrap() * exit.transmittance;
    assert!((sample.beam.power().i / expected - 1.0).abs() < 3e-12);
}

#[test]
fn g0_coated_diffuse_keeps_weak_conditional_power() {
    let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for (cosine, power) in [(1e-20, 1e-310), (0.01, f64::from_bits(8))] {
        let input = beam(Vec3::new((1.0_f64-cosine*cosine).sqrt(), 0.0, -cosine),
            Stokes::unpolarized(power));
        let mass = coating.branch_probabilities(cosine).unwrap();
        assert!(mass[0] > mass[1]);
        let sample = input.sample_coated_diffuse(normal, coating, 1.0,
            [0.3, 0.2, 0.5 * mass[1]]).unwrap();
        assert!(!sample.delta);
        let exit = crate::dielectric::fresnel_dielectric(
            sample.beam.frame().direction().z, 1.0, 1.5).unwrap();
        let expected = power * (coating.diffuse_weight(1.0).unwrap() * exit.transmittance);
        let actual = sample.beam.power().i;
        assert!(actual > 0.0, "recoverable conditional power was erased: {cosine}, {power}");
        assert!((actual / expected - 1.0).abs() < 3e-12
            || (actual - expected).abs() <= f64::from_bits(2), "{actual} != {expected}");
    }
}
