use super::*;
use crate::compose::compose_chain;
use crate::{RadiusProfile, RiccatiSteps, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;
use fs_math::det;
use std::sync::Arc;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn material() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-20),
        "{a:?} != {b:?}"
    );
}

#[test]
fn g0_lining_preserves_original_source_fluid_and_phasor() {
    let material = material();
    let air = air();
    let omega = core::f64::consts::TAU * 700.0;
    let lining = PorousLining::try_new(&material, &air, omega, 0.025).unwrap();
    let source = material
        .surface_impedance_rigid_fluid(&air, omega, 0.025)
        .unwrap();
    assert!(std::ptr::eq(lining.material(), &material));
    assert_eq!(
        lining.response().fluid(),
        &air.thermoviscous_state().unwrap()
    );
    assert_eq!(lining.thickness_m(), 0.025);
    assert_eq!(lining.surface_impedance_pa_s_m(), source.conj());
    // Planar local impedance is radius-independent; the actual perimeter is
    // applied exactly once by the existing duct geometry owner.
    for radius in [0.01, 0.02, 0.04] {
        assert_eq!(
            lining.specific_impedance(radius, omega).unwrap(),
            source.conj()
        );
    }
}

#[test]
fn g1_lined_cylinder_matches_independent_line_equation_and_wall_power() {
    let material = material();
    let air = air();
    let omega = core::f64::consts::TAU * 1_000.0;
    let (radius, length) = (0.02, 0.4);
    let lining = PorousLining::try_new(&material, &air, omega, 0.025).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let perimeter = core::f64::consts::TAU * radius;
    let ys = lining.surface_impedance_pa_s_m().recip();
    // Independent pressure/volume-flow balance, exp(-i omega t):
    // dp/dx = -z U; dU/dx = -y p; k^2 = -z y.
    let z = C64::new(0.0, -omega * air.density / area);
    let y = C64::new(0.0, -omega * area / (air.density * air.sound_speed.powi(2)))
        + ys.scale(perimeter);
    let k = (-z * y).sqrt();
    let zc = C64::from_re(omega * air.density / area) / k;
    assert!(k.re > 0.0 && k.im > 0.0 && zc.re > 0.0);
    let e =
        C64::new(det::cos(k.re * length), det::sin(k.re * length)).scale(det::exp(-k.im * length));
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let options = ComposeOptions::default();
    let prepared = lining
        .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
        .unwrap();
    let solved = compose_chain(&prepared.blocks, Some(zc), options).unwrap();
    near(solved.z_in, zc, 3e-12);
    near(solved.u_out_over_u_in, e, 3e-12);
    near(solved.p_out_over_u_in, zc * e, 3e-12);
    // With U_in=1, missing axial power equals integrated wall heat. This
    // fails if specific impedance is not conjugated or its perimeter is lost.
    let input_power = 0.5 * solved.z_in.re;
    let output_power = 0.5 * (solved.p_out_over_u_in * solved.u_out_over_u_in.conj()).re;
    let wall_power =
        0.5 * perimeter * ys.re * zc.norm_sq() * (-(-2.0 * k.im * length).exp_m1()) / (2.0 * k.im);
    assert!(wall_power > 0.0 && output_power > 0.0 && output_power < input_power);
    assert!((input_power - output_power - wall_power).abs() < 2e-11 * input_power);
}

#[test]
fn g3_lining_is_unchanged_by_segmenting_or_using_a_constant_profile() {
    let material = material();
    let omega = core::f64::consts::TAU * 700.0;
    let lining = PorousLining::try_new(&material, &air(), omega, 0.025).unwrap();
    let (radius, length) = (0.02, 0.4);
    let ducts = [
        Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        Duct {
            segments: vec![
                Segment::Cylinder {
                    radius,
                    length: length / 4.0
                };
                4
            ],
        },
        Duct {
            segments: vec![Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(radius, 0.0, length).unwrap()),
                steps: RiccatiSteps(8),
            }],
        },
    ];
    let options = ComposeOptions::default();
    let load = C64::new(2e5, -4e4);
    let mut reference = None;
    for duct in ducts {
        let chain = lining
            .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
            .unwrap();
        let got = compose_chain(&chain.blocks, Some(load), options).unwrap();
        if let Some((z, transfer)) = reference {
            near(got.z_in, z, 2e-11);
            near(got.p_out_over_u_in, transfer, 2e-11);
        } else {
            reference = Some((got.z_in, got.p_out_over_u_in));
        }
    }
}

#[test]
fn g0_source_support_geometry_and_cancellation_remain_explicit() {
    let material = material();
    let omega = core::f64::consts::TAU * 700.0;
    let lining = PorousLining::try_new(&material, &air(), omega, 0.025).unwrap();
    assert!(PorousLining::try_new(&material, &air(), 1.0, 0.025).is_err());
    for thickness in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(PorousLining::try_new(&material, &air(), omega, thickness).is_err());
    }
    let mut stale = air();
    stale.pressure *= 2.0;
    assert!(PorousLining::try_new(&material, &stale, omega, 0.025).is_err());
    for radius in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(lining.specific_impedance(radius, omega).is_err());
    }
    for frequency in [0.0, omega * 1.01, f64::NAN, f64::INFINITY] {
        assert!(lining.specific_impedance(0.02, frequency).is_err());
    }
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.02,
            length: 0.4,
        }],
    };
    assert!(matches!(
        lining.duct_chain(
            &duct,
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        }))
    ));
}
