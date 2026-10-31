use super::*;
use crate::compose::compose_chain;
use crate::{RadiusProfile, RiccatiSteps, Segment};
use core::cell::Cell;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;
use fs_material::porous::layers::Layer;
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

#[test]
fn g1_layered_lining_reproduces_the_original_single_layer_chain() {
    let law = material();
    let air = air();
    let omega = core::f64::consts::TAU * 700.0;
    let original = PorousLining::try_new(&law, &air, omega, 0.025).unwrap();
    let layers = RigidBackedLayers::new(
        &air,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.025,
        }],
    )
    .unwrap();
    let lining = PorousLayersLining::try_new(&layers).unwrap();
    assert!(std::ptr::eq(lining.layers(), &layers));
    near(
        lining.specific_impedance(0.02, omega).unwrap(),
        original.specific_impedance(0.02, omega).unwrap(),
        3e-14,
    );
    near(
        lining.surface_admittance_m_pa_s(),
        original.specific_admittance(0.02, omega).unwrap(),
        3e-14,
    );
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.02,
            length: 0.4,
        }],
    };
    let options = ComposeOptions::default();
    let before = original
        .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
        .unwrap();
    let after = lining
        .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
        .unwrap();
    let load = Some(C64::new(2e5, -4e4));
    let a = compose_chain(&before.blocks, load, options).unwrap();
    let b = compose_chain(&after.blocks, load, options).unwrap();
    near(a.z_in, b.z_in, 3e-12);
    near(a.p_out_over_u_in, b.p_out_over_u_in, 3e-12);
}

#[test]
fn g1_cavity_backed_layers_drive_the_original_pressure_flow_equations() {
    let law = material();
    let air = air();
    let omega = core::f64::consts::TAU * 700.0;
    let stack = RigidBackedLayers::new(
        &air,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
            Layer::FluidGap { thickness_m: 0.04 },
        ],
    )
    .unwrap();
    let lining = PorousLayersLining::try_new(&stack).unwrap();
    let bare = PorousLining::try_new(&law, &air, omega, 0.02).unwrap();
    assert!(
        (lining.surface_admittance_m_pa_s() - bare.specific_admittance(0.02, omega).unwrap()).abs()
            > 1e-5
    );
    let (radius, length) = (0.02, 0.4);
    let area = core::f64::consts::PI * radius * radius;
    let ys = stack
        .reflection(1.0)
        .unwrap()
        .normalized_surface_admittance()
        .unwrap()
        .conj()
        / C64::from_re(air.density * air.sound_speed);
    let series = C64::new(0.0, -omega * air.density / area);
    let shunt = C64::new(0.0, -omega * area / (air.density * air.sound_speed.powi(2)))
        + ys.scale(core::f64::consts::TAU * radius);
    let k = (-series * shunt).sqrt();
    let zc = C64::from_re(omega * air.density / area) / k;
    let expected =
        C64::new(det::cos(k.re * length), det::sin(k.re * length)).scale(det::exp(-k.im * length));
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let options = ComposeOptions::default();
    let prepared = lining
        .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
        .unwrap();
    let result = compose_chain(&prepared.blocks, Some(zc), options).unwrap();
    assert!(k.im > 0.0 && ys.re > 0.0);
    near(result.z_in, zc, 3e-12);
    near(result.u_out_over_u_in, expected, 3e-12);
    let output = (result.p_out_over_u_in * result.u_out_over_u_in.conj()).re;
    assert!(output > 0.0 && output < result.z_in.re);
}

#[test]
fn g0_empty_and_zero_depth_layers_are_exact_rigid_admittance() {
    let law = material();
    let air = air();
    let omega = core::f64::consts::TAU * 700.0;
    let empty = RigidBackedLayers::new(&air, omega, &[]).unwrap();
    let zero = RigidBackedLayers::new(
        &air,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.0,
            },
            Layer::FluidGap { thickness_m: 0.0 },
        ],
    )
    .unwrap();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.02,
                length: 0.2,
            },
            Segment::Cone {
                inlet_radius: 0.02,
                outlet_radius: 0.03,
                length: 0.2,
            },
        ],
    };
    let options = ComposeOptions::default();
    for stack in [&empty, &zero] {
        let lining = PorousLayersLining::try_new(stack).unwrap();
        assert_eq!(lining.specific_admittance(0.02, omega).unwrap(), C64::ZERO);
        assert!(lining.specific_impedance(0.02, omega).is_err());
        let rigid = fluid::duct_chain(
            &duct,
            stack.fluid(),
            omega,
            WallModel::Rigid,
            ConeSlices::Fixed(8),
            options,
            || false,
        )
        .unwrap();
        let actual = lining
            .duct_chain(&duct, ConeSlices::Fixed(8), options, || false)
            .unwrap();
        let load = Some(C64::new(2e5, -4e4));
        let rigid = compose_chain(&rigid.blocks, load, options).unwrap();
        let actual = compose_chain(&actual.blocks, load, options).unwrap();
        assert_eq!(actual.z_in, rigid.z_in);
        assert_eq!(actual.p_out_over_u_in, rigid.p_out_over_u_in);
    }
}

struct CountedLiningFluid {
    state: ThermoviscousFluidState,
    calls: Cell<usize>,
}
impl FluidAcoustics for CountedLiningFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedLiningFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for CountedLiningFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn g0_layered_lining_keeps_one_source_sample_frequency_and_cancellation() {
    let fluid = CountedLiningFluid {
        state: air().thermoviscous_state().unwrap(),
        calls: Cell::new(0),
    };
    let law = material();
    let omega = core::f64::consts::TAU * 700.0;
    let stack = RigidBackedLayers::new(
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
            Layer::FluidGap { thickness_m: 0.04 },
        ],
    )
    .unwrap();
    let lining = PorousLayersLining::try_new(&stack).unwrap();
    for radius in [0.01, 0.02, 0.04] {
        lining.specific_admittance(radius, omega).unwrap();
    }
    for radius in [0.0, -0.1, f64::NAN] {
        assert!(lining.specific_admittance(radius, omega).is_err());
    }
    for frequency in [0.0, omega * 1.01, f64::NAN, f64::INFINITY] {
        assert!(lining.specific_admittance(0.02, frequency).is_err());
    }
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.02,
            length: 0.4,
        }],
    };
    lining
        .duct_chain(
            &duct,
            ConeSlices::Fixed(8),
            ComposeOptions::default(),
            || false,
        )
        .unwrap();
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
    assert_eq!(fluid.calls.get(), 1);
    assert_eq!(lining.layers().fluid(), &fluid.state);
}

#[test]
fn g0_default_admittance_keeps_existing_impedance_admission() {
    struct Legacy(C64);
    impl WallLaw for Legacy {
        fn specific_impedance(&self, _: f64, _: f64) -> Result<C64, DuctError> {
            Ok(self.0)
        }
    }
    let z = C64::new(20.0, -3.0);
    assert_eq!(Legacy(z).specific_admittance(0.02, 1.0).unwrap(), z.recip());
    for z in [
        C64::ZERO,
        C64::new(-1.0, 2.0),
        C64::new(f64::NAN, 1.0),
        C64::new(1.0, f64::INFINITY),
    ] {
        assert!(Legacy(z).specific_admittance(0.02, 1.0).is_err());
    }
}
