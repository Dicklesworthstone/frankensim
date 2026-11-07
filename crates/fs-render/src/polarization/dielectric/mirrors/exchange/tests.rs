use super::*;
use crate::conductor::{
    CoatingIorSample, CoherentAbsorbingCoating, ConductorDataStatus, ConductorIorSample,
    ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_geom::Point3;

fn glass(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn metal(eta: f64) -> ConductorOptics {
    let source = ConductorSource::try_new(
        fs_blake3::hash_domain("org.frankensim.test.mirror-exchange", b"synthetic"),
        ConductorDataStatus::Representative,
    ).unwrap();
    ConductorOptics::try_new(core::array::from_fn(|i| {
        ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, 3.0).unwrap()
    }), source).unwrap()
}

fn beam(direction: Vec3, power: Stokes, medium: DielectricGlass) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(1.0, 2.0, 3.0),
        PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        550.0,
        medium,
        power,
    ).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-12 * expected.abs().max(1e-300),
        "actual={actual:.17e}, expected={expected:.17e}");
}

#[test]
fn g1_mirror_balance_matches_normal_pressure_and_independent_oblique_fresnel() {
    let optics = metal(0.2);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let input = beam(Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0), glass(1.0));
    let output = input.reflect_conductor_with_balance(normal, &optics).unwrap();
    let r = (0.8 * 0.8 + 9.0) / (1.2 * 1.2 + 9.0);
    close(output.reflection().beam.power().i, r);
    close(output.absorbed_power_w(), 1.0 - r);
    let force = output.vacuum_force_n().unwrap();
    assert_eq!((force.x, force.y), (0.0, 0.0));
    close(force.z, -(1.0 + r) / LIGHT_M_S);

    // Independent complex Fresnel evaluation at cos(theta)=0.6, n=0.2+3i:
    // q=sqrt(n²-sin²(theta)), rs=(cos-q)/(cos+q),
    // rp=(n²*cos-q)/(n²*cos+q). These are not native R/A self-agreement.
    for (q, r, fx, fz) in [
        (2.0, 0.954_819_051_576_397_3, 0.036_144_758_738_882_12, -1.172_891_430_945_838_5),
        (-2.0, 0.884_672_979_286_979_7, 0.092_261_616_570_416_2, -1.130_803_787_572_187_7),
    ] {
        let input = beam(Vec3::new(0.8, 0.0, -0.6), Stokes::try_new(2.0, q, 0.0, 0.0).unwrap(), glass(1.0));
        let output = input.reflect_conductor_with_balance(normal, &optics).unwrap();
        assert_eq!(*output.reflection(), input.reflect_conductor(normal, &optics).unwrap());
        close(output.reflection().beam.power().i, 2.0 * r);
        close(output.absorbed_power_w(), 2.0 * (1.0 - r));
        let force = output.vacuum_force_n().unwrap();
        close(force.x, 2.0 * fx / LIGHT_M_S);
        assert_eq!(force.y, 0.0);
        close(force.z, 2.0 * fz / LIGHT_M_S);
    }
}

#[test]
fn g1_mirror_balance_retains_weak_loss_and_finite_force_at_large_power() {
    let optics = metal(1e-20);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for power in [1.0, f64::MAX] {
        let input = beam(Vec3::new(0.0, 0.0, -1.0),
            Stokes::try_new(power, power, 0.0, 0.0).unwrap(), glass(1.0));
        let output = input.reflect_conductor_with_balance(normal, &optics).unwrap();
        assert_eq!(output.reflection().beam.power().i, power);
        // Positive Poynting flux at normal incidence: A=4 eta/((1+eta)²+k²).
        close(output.absorbed_power_w(), power * 4e-21);
        close(output.vacuum_force_n().unwrap().z, -2.0 * (power / LIGHT_M_S));
    }
    let input = beam(Vec3::new(0.8, 0.0, -0.6), Stokes::unpolarized(1.0), glass(1.0));
    let output = input.reflect_conductor_with_balance(normal, &optics).unwrap();
    assert!(output.absorbed_power_w() > 0.0);
    // The tangential force must retain the tiny physical loss, irrespective
    // of whether the reflected intensity rounds to 1 or an adjacent float.
    close(output.vacuum_force_n().unwrap().x, 0.8 * output.absorbed_power_w() / LIGHT_M_S);
    assert!(output.vacuum_force_n().unwrap().x > 0.0);
}

#[test]
fn g1_mirror_balance_includes_coherent_absorbing_film_without_double_loss() {
    let bare = metal(0.2);
    let film = CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.45, 0.07).unwrap()
        }),
        100e-9,
        bare.source(),
    ).unwrap();
    let optics = bare.with_absorbing_coating(film);
    let input = beam(Vec3::new(0.0, 0.0, -1.0), Stokes::try_new(3.0, 0.0, 1.8, 2.4).unwrap(), glass(1.0));
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let output = input.reflect_conductor_with_balance(normal, &optics).unwrap();
    // Independent complex Airy result for this absorbing film/substrate.
    let r = 0.515_088_324_939_222_6;
    assert_eq!(*output.reflection(), input.reflect_conductor(normal, &optics).unwrap());
    close(output.absorbed_power_w(), 3.0 * (1.0 - r));
    close(output.vacuum_force_n().unwrap().z, -3.0 * (1.0 + r) / LIGHT_M_S);
}

#[test]
fn g3_mirror_exchange_is_reference_covariant_and_preserves_host_authority() {
    let optics = metal(0.2);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let input = beam(Vec3::new(0.8, 0.0, -0.6), Stokes::try_new(2.0, 0.4, -0.7, 0.2).unwrap(), glass(1.0));
    let original = input.reflect_conductor_with_balance(normal, &optics).unwrap();
    let rotated = input.with_reference_axis(Vec3::new(1.0, 2.0, 3.0)).unwrap();
    let changed = rotated.reflect_conductor_with_balance(normal.scale(-f64::MAX), &optics).unwrap();
    close(changed.absorbed_power_w(), original.absorbed_power_w());
    let a = original.vacuum_force_n().unwrap();
    let b = changed.vacuum_force_n().unwrap();
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);

    let immersed = beam(input.frame().direction(), input.power(), glass(1.33));
    let exchange = immersed.reflect_conductor_with_balance(normal, &optics).unwrap();
    assert_eq!(exchange.vacuum_force_n(), None);
    assert!(exchange.absorbed_power_w() > 0.0);
    assert_eq!(*exchange.reflection(), immersed.reflect_conductor(normal, &optics).unwrap());
    close(exchange.absorbed_power_w() + exchange.reflection().beam.power().i, input.power().i);
}

#[test]
fn g0_mirror_balance_handles_dark_beams_and_refuses_invalid_geometry() {
    let optics = metal(0.2);
    let input = beam(Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(0.0), glass(1.0));
    let output = input.reflect_conductor_with_balance(Vec3::new(0.0, 0.0, 1.0), &optics).unwrap();
    assert_eq!(output.absorbed_power_w(), 0.0);
    assert_eq!(output.vacuum_force_n(), Some(Vec3::new(0.0, 0.0, 0.0)));
    for normal in [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)] {
        assert_eq!(input.reflect_conductor_with_balance(normal, &optics).unwrap_err(),
            input.reflect_conductor(normal, &optics).unwrap_err());
    }
}
