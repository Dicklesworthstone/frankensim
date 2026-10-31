//! Physics controls for the source-bound multilayer reflection owner.
use core::cell::Cell;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_material::porous::layers::{Layer, MAX_LAYERS, RigidBackedLayers};
use fs_math::c64::C64;

fn fluid() -> ThermoviscousFluidState {
    // Explicit analytic fixture, not measured ambient data.
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    ).unwrap()
}
fn foam() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual-expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}

#[test]
fn single_layer_reproduces_original_owner_in_actual_gas_states() {
    let law = foam();
    for temperature in [293.15, 323.15] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap();
        for frequency in [150.0, 500.0, 1000.0, 1600.0] {
            let omega = core::f64::consts::TAU*frequency;
            let stack = RigidBackedLayers::new(&gas, omega,
                &[Layer::Porous { law: &law, thickness_m: 0.06 }]).unwrap();
            let sample = law.effective_fluid(&gas, omega).unwrap();
            assert_eq!(stack.material_response(0), Some(&sample));
            let zs = sample.surface_impedance_rigid(0.06).unwrap();
            let z0 = C64::from_re(stack.fluid().acoustic().characteristic_impedance_pa_s_m());
            let response = stack.reflection(1.0).unwrap();
            near(response.pressure_reflection, (zs-z0)/(zs+z0), 3e-14);
            near(response.normalized_surface_admittance().unwrap(), z0/zs, 3e-14);
        }
    }
}

// Independent conventional complex transfer matrices. These intentionally use
// growing cos/sin, only in a well-conditioned oracle, not the production E form.
fn cos_sin(z: C64) -> (C64, C64) {
    (
        C64::new(z.re.cos()*z.im.cosh(), -z.re.sin()*z.im.sinh()),
        C64::new(z.re.sin()*z.im.cosh(), z.re.cos()*z.im.sinh()),
    )
}
fn matrix_reflection(stack: &RigidBackedLayers<'_>, cosine: f64) -> C64 {
    let z0 = stack.fluid().acoustic().characteristic_impedance_pa_s_m();
    let kt = stack.exterior_wavenumber_rad_m()*(1.0-cosine*cosine).sqrt();
    let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    for i in 0..stack.layer_count() {
        let (kn, ratio, d) = match stack.layer(i).unwrap() {
            Layer::FluidGap { thickness_m } =>
                (C64::from_re(stack.exterior_wavenumber_rad_m()*cosine), C64::ONE, thickness_m),
            Layer::Porous { thickness_m, .. } => {
                let (zc, k) = stack.material_response(i).unwrap().characteristic();
                let mut kn = (k*k-C64::from_re(kt*kt)).sqrt();
                if kn.im > 0.0 { kn = -kn; }
                (kn, (zc*k/kn).scale(cosine/z0), thickness_m)
            }
        };
        let (cos, sin) = cos_sin(kn.scale(d));
        let j_sin = C64::new(0.0, 1.0)*sin;
        let layer = [[cos, ratio*j_sin], [j_sin/ratio, cos]];
        let mut product = [[C64::ZERO; 2]; 2];
        for row in 0..2 {
            for col in 0..2 {
                product[row][col] = m[row][0]*layer[0][col] + m[row][1]*layer[1][col];
            }
        }
        m = product;
    }
    (m[0][0]-m[1][0])/(m[0][0]+m[1][0])
}

#[test]
fn ordered_porous_gap_stack_matches_independent_transfer_matrices() {
    let fluid = fluid();
    let first = foam();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous { law: &first, thickness_m: 0.03 },
        Layer::FluidGap { thickness_m: 0.05 },
        Layer::Porous { law: &second, thickness_m: 0.02 },
    ];
    for frequency in [150.0, 500.0, 1000.0, 1600.0] {
        let omega = core::f64::consts::TAU*frequency;
        let stack = RigidBackedLayers::new(&fluid, omega, &layers).unwrap();
        let reverse = RigidBackedLayers::new(&fluid, omega,
            &[layers[2], layers[1], layers[0]]).unwrap();
        for cosine in [1.0, 0.85, 0.5, 0.05] {
            let actual = stack.reflection(cosine).unwrap();
            let expected = matrix_reflection(&stack, cosine);
            near(actual.pressure_reflection, expected, 4e-14);
            assert!((actual.absorbed_fraction-(1.0-expected.norm_sq())).abs() < 4e-14);
            assert!(actual.absorbed_fraction >= 0.0);
        }
        assert!((stack.reflection(1.0).unwrap().pressure_reflection
            -reverse.reflection(1.0).unwrap().pressure_reflection).abs() > 1e-4);
    }
}

#[test]
fn subdivision_and_zero_thickness_preserve_coherent_reflection() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU*800.0;
    let whole = RigidBackedLayers::new(&fluid, omega,
        &[Layer::Porous { law: &law, thickness_m: 0.06 }]).unwrap();
    let split = RigidBackedLayers::new(&fluid, omega, &[
        Layer::Porous { law: &law, thickness_m: 0.01 },
        Layer::FluidGap { thickness_m: 0.0 },
        Layer::Porous { law: &law, thickness_m: 0.02 },
        Layer::Porous { law: &law, thickness_m: 0.03 },
    ]).unwrap();
    for cosine in [1.0, 0.5, 0.05] {
        near(split.reflection(cosine).unwrap().pressure_reflection,
            whole.reflection(cosine).unwrap().pressure_reflection, 4e-14);
    }
}

#[test]
fn lossless_gaps_preserve_round_trip_phase_poles_and_zero_work() {
    let fluid = fluid();
    // k0=1 rad/m makes quarter/half-wave controls explicit.
    let omega = fluid.acoustic().sound_speed_m_s();
    for d in [0.0, 0.1, core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
        let stack = RigidBackedLayers::new(&fluid, omega,
            &[Layer::FluidGap { thickness_m: d }]).unwrap();
        for cosine in [1.0, 0.5, 1e-200] {
            let response = stack.reflection(cosine).unwrap();
            let phase = 2.0*cosine*d;
            near(response.pressure_reflection, C64::new(phase.cos(), -phase.sin()), 3e-15);
            assert_eq!(response.absorbed_fraction, 0.0);
            if let Ok(beta) = response.normalized_surface_admittance() {
                assert_eq!(beta.re, 0.0);
            }
        }
    }
    let empty = RigidBackedLayers::new(&fluid, omega, &[]).unwrap();
    assert_eq!(empty.reflection(0.25).unwrap().pressure_reflection, C64::ONE);
    assert_eq!(empty.reflection(0.25).unwrap().normalized_surface_admittance().unwrap(), C64::ZERO);
}

#[test]
fn attenuating_thick_limit_and_extreme_grazing_keep_finite_positive_work() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU*500.0;
    // Conventional cosh matrices overflow; the attenuating half-space limit is finite.
    let thick = RigidBackedLayers::new(&fluid, omega,
        &[Layer::Porous { law: &law, thickness_m: 1e308 }]).unwrap();
    let (zc, k) = thick.material_response(0).unwrap().characteristic();
    let z0 = fluid.acoustic().characteristic_impedance_pa_s_m();
    for cosine in [1.0, 0.5, 0.01] {
        let kt = thick.exterior_wavenumber_rad_m()*(1.0-cosine*cosine).sqrt();
        let mut kn = (k*k-C64::from_re(kt*kt)).sqrt();
        if kn.im > 0.0 { kn = -kn; }
        let ratio = (zc*k/kn).scale(cosine/z0);
        let actual = thick.reflection(cosine).unwrap();
        near(actual.pressure_reflection, (ratio-C64::ONE)/(ratio+C64::ONE), 3e-14);
        assert!(actual.absorbed_fraction > 0.0);
    }
    let stack = RigidBackedLayers::new(&fluid, omega, &[
        Layer::Porous { law: &law, thickness_m: 0.02 },
        Layer::Porous { law: &law, thickness_m: 0.04 },
        Layer::FluidGap { thickness_m: 0.1 },
    ]).unwrap();
    let baseline = stack.reflection(1e-100).unwrap();
    for cosine in [1e-200, 1e-300] {
        let response = stack.reflection(cosine).unwrap();
        assert!(response.absorbed_fraction > 0.0);
        assert!((response.absorbed_fraction/cosine
            -baseline.absorbed_fraction/1e-100).abs() < 2e-13);
        assert!((response.pressure_reflection.im/cosine
            -baseline.pressure_reflection.im/1e-100).abs() < 2e-13);
    }
}

#[test]
fn weak_absorption_survives_unit_reflected_magnitude_roundoff() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU*500.0;
    let d = 1e-100;
    let stack = RigidBackedLayers::new(&fluid, omega,
        &[Layer::Porous { law: &law, thickness_m: d }]).unwrap();
    let response = stack.reflection(1.0).unwrap();
    let bulk = stack.material_response(0).unwrap().dynamic_bulk_modulus_pa();
    let beta = C64::new(0.0, omega*d*fluid.acoustic().characteristic_impedance_pa_s_m())/bulk;
    assert_eq!(response.pressure_reflection.norm_sq(), 1.0);
    assert!(response.absorbed_fraction > 0.0);
    assert!((response.absorbed_fraction/(4.0*beta.re)-1.0).abs() < 2e-13);
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    samples: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.samples.set(self.samples.get()+1);
        Ok(self.state)
    }
}

#[test]
fn all_layers_share_one_actual_snapshot_and_keep_source_domain() {
    let fluid = CountedFluid { state: fluid(), samples: Cell::new(0) };
    let law = foam();
    let layers = [
        Layer::Porous { law: &law, thickness_m: 0.03 },
        Layer::FluidGap { thickness_m: 0.05 },
        Layer::Porous { law: &law, thickness_m: 0.01 },
    ];
    let stack = RigidBackedLayers::new(&fluid, core::f64::consts::TAU*500.0, &layers).unwrap();
    for cosine in [1.0, 0.5, 0.01] { stack.reflection(cosine).unwrap(); }
    assert_eq!(fluid.samples.get(), 1);
    assert_eq!(stack.material_response(0).unwrap().fluid(), &fluid.state);
    assert_eq!(stack.material_response(2).unwrap().fluid(), &fluid.state);
    assert!(stack.material_response(1).is_none());
    assert!(matches!(stack.layer(0), Some(Layer::Porous { law: source, .. }) if core::ptr::eq(source, &law)));
    assert!(stack.layer(3).is_none());
    for cosine in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(stack.reflection(cosine).is_err());
    }
    assert!(RigidBackedLayers::new(&fluid, core::f64::consts::TAU*149.0,
        &[Layer::Porous { law: &law, thickness_m: 0.0 }]).is_err());
    assert!(RigidBackedLayers::new(&fluid, 1.0,
        &[Layer::FluidGap { thickness_m: -0.1 }]).is_err());
    assert!(RigidBackedLayers::new(&fluid, 1.0,
        &[Layer::FluidGap { thickness_m: 0.0 }; MAX_LAYERS+1]).is_err());
    for omega in [0.0, f64::NAN, f64::INFINITY] {
        assert!(RigidBackedLayers::new(&fluid, omega, &[]).is_err());
    }
}
