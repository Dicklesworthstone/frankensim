use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::interface::fluid::FluidDriveLimits;
use fs_material::elastic::acoustic::viscoelastic::interface::FluidViscoelasticInterface;
use fs_material::elastic::acoustic::viscoelastic::interface::layer::{
    FluidViscoelasticLayer, ViscoelasticLayerWorkBudget,
};
use fs_material::elastic::acoustic::viscoelastic::{ViscoelasticSlabError, ViscoelasticSlabSource};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::{c64::C64, det};
use std::cell::Cell;

fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
) -> ViscoelasticSlabSource<'a> {
    ViscoelasticSlabSource {
        bulk,
        shear,
        density_kg_m3: 4.0,
        frequency_band_rad_s: (10.0, 100.0),
        maximum_absolute_strain: 0.1,
        source_identity: ContentHash([9; 32]),
    }
}

fn limits() -> [FluidDriveLimits; 2] {
    [FluidDriveLimits::try_new(0.05, 0.05).unwrap(); 2]
}

fn budget() -> ViscoelasticLayerWorkBudget {
    ViscoelasticLayerWorkBudget {
        maximum_segments: 4096,
        relative_tolerance: 1e-10,
        absolute_tolerance_w_m2: 0.0,
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn exp_minus_i(z: C64) -> C64 {
    let (s, c) = det::sin_cos(z.re);
    C64::new(c, -s).scale(det::exp(z.im))
}

struct CountedFluid {
    state: FluidAcousticState,
    calls: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn g1_normal_complex_fabry_perot_and_original_sources() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let left = CountedFluid {
        state: FluidAcousticState::try_new(1.0, 2.0).unwrap(),
        calls: Cell::new(0),
    };
    let right = CountedFluid {
        state: FluidAcousticState::try_new(2.0, 3.0).unwrap(),
        calls: Cell::new(0),
    };
    for (index, h) in [1e-9, 0.27, 3.0].into_iter().enumerate() {
        let layer = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&left, &right],
            20.0,
            1.0,
            h,
            2,
            &mut || false,
        )
        .unwrap();
        assert!(core::ptr::eq(layer.source().bulk, &bulk));
        assert!(core::ptr::eq(layer.source().shear, &shear));
        assert_eq!(layer.source().source_identity, ContentHash([9; 32]));
        let [k, g] = layer.moduli_pa();
        let modulus = k + g.scale(4.0 / 3.0);
        let impedance = modulus.sqrt().scale(2.0);
        let q = C64::from_re(2.0) / modulus.sqrt();
        let phase = exp_minus_i(q.scale(20.0 * h));
        let r01 = (impedance - C64::from_re(2.0)) / (impedance + C64::from_re(2.0));
        let r12 = (C64::from_re(6.0) - impedance) / (C64::from_re(6.0) + impedance);
        let denominator = C64::ONE + r01 * r12 * phase * phase;
        let expected = [
            (r01 + r12 * phase * phase) / denominator,
            (C64::ONE - r12) * (C64::ONE - r01) * phase / denominator,
            (C64::ONE + r01) * (C64::ONE + r12) * phase / denominator,
            (-r12 - r01 * phase * phase) / denominator,
        ];
        for (actual, expected) in layer.pressure_scattering().into_iter().zip(expected) {
            near(actual, expected, 7e-14);
        }
        let power_wave = layer.power_normalized_scattering().unwrap();
        near(power_wave[1], power_wave[2], 4e-14);
        let driven = layer
            .drive(
                [C64::new(0.02, 0.01), C64::new(0.005, -0.003)],
                limits(),
                &mut || false,
            )
            .unwrap();
        let faces = driven.solid_face_fields().unwrap();
        for (side, field) in faces.into_iter().enumerate() {
            near(
                field.traction_pa[2],
                -driven.fluid_face_pressure_pa()[side],
                3e-15,
            );
            near(
                field.velocity_m_s[2],
                driven.fluid_face_velocity_m_s()[side][2],
                3e-15,
            );
            near(field.traction_pa[0], C64::ZERO, 3e-15);
            near(field.velocity_m_s[0], C64::ZERO, 3e-15);
        }
        let work = driven
            .absorption_between(0.0, h, budget(), &mut || false)
            .unwrap();
        let heat = work.bulk_absorbed_power_w_m2 + work.shear_absorbed_power_w_m2;
        let incident: f64 = driven.power().incident_power_w_m2.iter().sum();
        assert!((heat - driven.power().net_absorbed_power_w_m2).abs() < 2e-12 * incident);
        assert_eq!(left.calls.get(), index + 1);
        assert_eq!(right.calls.get(), index + 1);
    }
}

#[test]
fn g1_lossless_half_wave_resonance_and_critical_roots() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let h = core::f64::consts::PI * 5.0 / 20.0;
    let layer = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&fluid, &fluid],
        20.0,
        1.0,
        h,
        0,
        &mut || false,
    )
    .unwrap();
    near(layer.pressure_scattering()[0], C64::ZERO, 2e-14);
    near(layer.pressure_scattering()[2], -C64::ONE, 2e-14);
    // Choose binary-exact p=1/2. P critical: rho/M=1/4; S critical: rho/G=1/4.
    let cosine = 0.8;
    let speed = 2.0 * det::sqrt((1.0 - cosine) * (1.0 + cosine));
    let critical_fluid = FluidAcousticState::try_new(1.0, speed).unwrap();
    for (kind, g) in [(0, 6.0), (1, 16.0)] {
        let k = GeneralizedMaxwell::new(8.0, vec![]).unwrap();
        let g = GeneralizedMaxwell::new(g, vec![]).unwrap();
        let layer = FluidViscoelasticLayer::try_new(
            source(&k, &g),
            [&critical_fluid, &critical_fluid],
            20.0,
            cosine,
            0.13,
            0,
            &mut || false,
        )
        .unwrap();
        assert_eq!(layer.tangential_slowness_s_m(), 0.5);
        assert_eq!(layer.normal_slownesses_s_m().unwrap()[kind], C64::ZERO);
        let s = layer.pressure_scattering();
        assert!((s[0].norm_sq() + s[2].norm_sq() - 1.0).abs() < 2e-13);
        let driven = layer
            .drive([C64::from_re(0.001), C64::ZERO], limits(), &mut || false)
            .unwrap();
        let field = driven.field(0.065, &mut || false).unwrap();
        assert_eq!(
            field.bulk_heat_density_w_m3 + field.shear_heat_density_w_m3,
            0.0
        );
        for nearby in [cosine - 1e-8, cosine + 1e-8] {
            let near_layer = FluidViscoelasticLayer::try_new(
                source(&k, &g),
                [&critical_fluid, &critical_fluid],
                20.0,
                nearby,
                0.13,
                0,
                &mut || false,
            )
            .unwrap();
            near(near_layer.pressure_scattering()[0], s[0], 2e-6);
            near(near_layer.pressure_scattering()[2], s[2], 2e-6);
        }
    }
}

#[test]
fn g1_normal_scattering_remains_defined_at_unforced_shear_resonances() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(16.0, vec![]).unwrap();
    let left = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let right = FluidAcousticState::try_new(2.0, 3.0).unwrap();
    // rho=4 and G=16 give c_S=2 exactly. These are unforced free shear
    // thickness resonances; the driven longitudinal solution stays unique.
    for order in [1.0, 2.0, 7.0] {
        let thickness = order * core::f64::consts::PI / 10.0;
        let layer = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&left, &right],
            20.0,
            1.0,
            thickness,
            0,
            &mut || false,
        )
        .unwrap();
        let [k, g] = layer.moduli_pa();
        let root = (k + g.scale(4.0 / 3.0)).sqrt();
        let impedance = root.scale(2.0);
        let phase = exp_minus_i((C64::from_re(2.0) / root).scale(20.0 * thickness));
        let r01 = (impedance - C64::from_re(2.0)) / (impedance + C64::from_re(2.0));
        let r12 = (C64::from_re(6.0) - impedance) / (C64::from_re(6.0) + impedance);
        let denominator = C64::ONE + r01 * r12 * phase * phase;
        near(
            layer.pressure_scattering()[0],
            (r01 + r12 * phase * phase) / denominator,
            5e-14,
        );
        near(
            layer.pressure_scattering()[2],
            (C64::ONE + r01) * (C64::ONE + r12) * phase / denominator,
            5e-14,
        );
        let driven = layer
            .drive(
                [C64::new(0.02, 0.01), C64::new(-0.01, 0.005)],
                limits(),
                &mut || false,
            )
            .unwrap();
        for fraction in [0.0, 0.2, 0.5, 0.9, 1.0] {
            let field = driven.field(fraction * thickness, &mut || false).unwrap();
            assert_eq!(field.velocity_m_s[0], C64::ZERO);
            assert_eq!(field.traction_pa[0], C64::ZERO);
            assert_eq!(field.strain_mandel[5], C64::ZERO);
        }
        let power = driven.power();
        assert!(
            power.net_absorbed_power_w_m2.abs()
                < 5e-14 * power.incident_power_w_m2.iter().sum::<f64>()
        );
    }
}

#[test]
fn g1_evanescent_ps_waves_tunnel_and_carry_coherent_flux() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let thin = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&fluid, &fluid],
        20.0,
        0.4,
        0.03,
        0,
        &mut || false,
    )
    .unwrap();
    assert!(
        thin.normal_slownesses_s_m()
            .unwrap()
            .iter()
            .all(|q| q.re == 0.0 && q.im < 0.0)
    );
    let s = thin.pressure_scattering();
    assert!(s[2].norm_sq() > 1e-3);
    assert!((s[0].norm_sq() + s[2].norm_sq() - 1.0).abs() < 3e-13);
    let driven = thin
        .drive([C64::from_re(0.01), C64::ZERO], limits(), &mut || false)
        .unwrap();
    let transmitted = driven.power().outgoing_power_w_m2[1];
    assert!(transmitted > 0.0);
    for z in [0.0, 0.007, 0.015, 0.03] {
        let field = driven.field(z, &mut || false).unwrap();
        assert!((field.normal_power_w_m2 / transmitted - 1.0).abs() < 3e-12);
        assert_eq!(
            field.bulk_heat_density_w_m3 + field.shear_heat_density_w_m3,
            0.0
        );
    }
    let thick = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&fluid, &fluid],
        20.0,
        0.4,
        40.0,
        0,
        &mut || false,
    )
    .unwrap();
    assert!(thick.pressure_scattering()[2].norm_sq() < 1e-100);
    assert!((thick.pressure_scattering()[0].norm_sq() - 1.0).abs() < 2e-14);
}

#[test]
fn g1_coherent_oblique_heat_closes_depth_flux_and_subdivision() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(20.0, 0.02)]).unwrap();
    let left = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let right = FluidAcousticState::try_new(1.4, 2.5).unwrap();
    for cosine in [0.98, 0.9, 0.8] {
        let layer = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&left, &right],
            20.0,
            cosine,
            0.7,
            2,
            &mut || false,
        )
        .unwrap();
        let s = layer.power_normalized_scattering().unwrap();
        near(s[1], s[2], 3e-14);
        let driven = layer
            .drive(
                [C64::new(0.02, 0.006), C64::new(-0.01, 0.004)],
                limits(),
                &mut || false,
            )
            .unwrap();
        let work = driven
            .absorption_between(0.0, 0.7, budget(), &mut || false)
            .unwrap();
        let first = driven
            .absorption_between(0.0, 0.23, budget(), &mut || false)
            .unwrap();
        let second = driven
            .absorption_between(0.23, 0.7, budget(), &mut || false)
            .unwrap();
        let total = work.bulk_absorbed_power_w_m2 + work.shear_absorbed_power_w_m2;
        assert!(work.bulk_absorbed_power_w_m2 > 0.0 && work.shear_absorbed_power_w_m2 > 0.0);
        assert!((total / driven.power().net_absorbed_power_w_m2 - 1.0).abs() < 2e-11);
        assert!(
            ((first.bulk_absorbed_power_w_m2 + second.bulk_absorbed_power_w_m2)
                / work.bulk_absorbed_power_w_m2
                - 1.0)
                .abs()
                < 2e-12
        );
        assert!(
            ((first.shear_absorbed_power_w_m2 + second.shear_absorbed_power_w_m2)
                / work.shear_absorbed_power_w_m2
                - 1.0)
                .abs()
                < 2e-12
        );
        for z in [0.1, 0.35, 0.6] {
            let field = driven.field(z, &mut || false).unwrap();
            let before = driven.field(z - 1e-5, &mut || false).unwrap();
            let after = driven.field(z + 1e-5, &mut || false).unwrap();
            let divergence = (before.normal_power_w_m2 - after.normal_power_w_m2) / 2e-5;
            let heat = field.bulk_heat_density_w_m3 + field.shear_heat_density_w_m3;
            assert!((divergence / heat - 1.0).abs() < 2e-8);
        }
        assert_eq!(
            work,
            driven
                .absorption_between(0.0, 0.7, budget(), &mut || false)
                .unwrap()
        );
    }
}

#[test]
fn g1_weak_thin_odd_standing_wave_heat_survives_power_roundoff() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(1e-12, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(1e-12, 0.05)]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let mut heating = [0.0; 2];
    for (index, h) in [1e-5, 2e-5].into_iter().enumerate() {
        let layer = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&fluid, &fluid],
            20.0,
            1.0,
            h,
            2,
            &mut || false,
        )
        .unwrap();
        let driven = layer
            .drive(
                [C64::from_re(0.02), C64::from_re(-0.02)],
                limits(),
                &mut || false,
            )
            .unwrap();
        let work = driven
            .absorption_between(0.0, h, budget(), &mut || false)
            .unwrap();
        heating[index] = work.bulk_absorbed_power_w_m2 + work.shear_absorbed_power_w_m2;
        let [k, g] = layer.moduli_pa();
        let modulus = k + g.scale(4.0 / 3.0);
        // Leading independent inertia balance: sigma=i*rho*omega*v*(z-h/2),
        // v -> 2*Pin/Zfluid. Integral |z-h/2|² dz = h³/12.
        let expected = 0.5 * 20.0 * modulus.im * (4.0_f64 * 20.0 * 0.02).powi(2)
            / modulus.norm_sq()
            * h
            * h
            * h
            / 12.0;
        assert!(heating[index] > 0.0);
        assert!((heating[index] / expected - 1.0).abs() < 2e-6);
    }
    assert!((heating[1] / heating[0] - 8.0).abs() < 2e-6);
}

#[test]
fn g3_opaque_absorbing_layer_recovers_original_half_space() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(20.0, 0.02)]).unwrap();
    let left = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let right = FluidAcousticState::try_new(1.4, 2.5).unwrap();
    let half = FluidViscoelasticInterface::try_new(
        source(&bulk, &shear),
        &left,
        20.0,
        0.9,
        2,
        &mut || false,
    )
    .unwrap();
    let layer = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&left, &right],
        20.0,
        0.9,
        1000.0,
        2,
        &mut || false,
    )
    .unwrap();
    near(
        layer.pressure_scattering()[0],
        half.pressure_reflection(),
        2e-14,
    );
    assert!(layer.pressure_scattering()[2].abs() < 1e-100);
    let a = half
        .drive(C64::from_re(0.01), limits()[0], &mut || false)
        .unwrap();
    let b = layer
        .drive([C64::from_re(0.01), C64::ZERO], limits(), &mut || false)
        .unwrap();
    for z in [0.0, 0.3, 1.0] {
        let fa = a.field(z, &mut || false).unwrap();
        let fb = b.field(z, &mut || false).unwrap();
        for i in 0..3 {
            near(fa.velocity_m_s[i], fb.velocity_m_s[i], 2e-15);
        }
        assert!((fa.bulk_heat_density_w_m3 / fb.bulk_heat_density_w_m3 - 1.0).abs() < 2e-12);
        assert!((fa.shear_heat_density_w_m3 / fb.shear_heat_density_w_m3 - 1.0).abs() < 2e-12);
    }
}

#[test]
fn g3_critical_regular_basis_crossover_preserves_scattering_fields_and_work() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(20.0, 0.02)]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let probe = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&fluid, &fluid],
        20.0,
        0.8,
        1.0,
        2,
        &mut || false,
    )
    .unwrap();
    for q in probe.normal_slownesses_s_m().unwrap() {
        let crossing = 0.5 / (20.0 * q.abs());
        let before = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&fluid, &fluid],
            20.0,
            0.8,
            crossing * (1.0 - 1e-10),
            2,
            &mut || false,
        )
        .unwrap();
        let after = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&fluid, &fluid],
            20.0,
            0.8,
            crossing * (1.0 + 1e-10),
            2,
            &mut || false,
        )
        .unwrap();
        for (a, b) in before
            .pressure_scattering()
            .into_iter()
            .zip(after.pressure_scattering())
        {
            near(a, b, 2e-9);
        }
        let input = [C64::new(0.02, 0.006), C64::new(-0.01, 0.004)];
        let a = before.drive(input, limits(), &mut || false).unwrap();
        let b = after.drive(input, limits(), &mut || false).unwrap();
        let fa = a.field(0.4 * before.thickness_m(), &mut || false).unwrap();
        let fb = b.field(0.4 * after.thickness_m(), &mut || false).unwrap();
        for i in 0..3 {
            near(fa.velocity_m_s[i], fb.velocity_m_s[i], 2e-11);
        }
        assert!((fa.bulk_heat_density_w_m3 / fb.bulk_heat_density_w_m3 - 1.0).abs() < 2e-8);
        assert!((fa.shear_heat_density_w_m3 / fb.shear_heat_density_w_m3 - 1.0).abs() < 2e-8);
        let wa = a
            .absorption_between(0.0, before.thickness_m(), budget(), &mut || false)
            .unwrap();
        let wb = b
            .absorption_between(0.0, after.thickness_m(), budget(), &mut || false)
            .unwrap();
        assert!((wa.bulk_absorbed_power_w_m2 / wb.bulk_absorbed_power_w_m2 - 1.0).abs() < 2e-8);
        assert!((wa.shear_absorbed_power_w_m2 / wb.shear_absorbed_power_w_m2 - 1.0).abs() < 2e-8);
    }
}

#[test]
fn g0_evanescent_exit_zero_thickness_and_bounded_refusals() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![]).unwrap();
    let left = CountedFluid {
        state: FluidAcousticState::try_new(1.0, 2.0).unwrap(),
        calls: Cell::new(0),
    };
    let right = FluidAcousticState::try_new(1.4, 5.0).unwrap();
    for h in [0.0, 0.2] {
        let layer = FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&left, &right],
            20.0,
            0.8,
            h,
            0,
            &mut || false,
        )
        .unwrap();
        assert!(layer.exterior_normal_cosines()[1].im < 0.0);
        assert!(layer.power_normalized_scattering().is_none());
        assert!(
            layer
                .drive([C64::ZERO, C64::from_re(0.01)], limits(), &mut || false)
                .is_err()
        );
        let driven = layer
            .drive([C64::from_re(0.01), C64::ZERO], limits(), &mut || false)
            .unwrap();
        assert_eq!(driven.power().outgoing_power_w_m2[1], 0.0);
        assert!(driven.fluid_face_pressure_pa()[1].abs() > 0.0);
        assert!((layer.pressure_scattering()[0].norm_sq() - 1.0).abs() < 3e-14);
        if h == 0.0 {
            assert!(driven.solid_face_fields().is_none());
            assert!(driven.field(0.0, &mut || false).is_err());
            assert_eq!(
                driven
                    .absorption_between(0.0, 0.0, budget(), &mut || false)
                    .unwrap()
                    .bulk_absorbed_power_w_m2,
                0.0
            );
            let y0 = C64::from_re(0.8 / 2.0);
            let y1 = layer.exterior_normal_cosines()[1].scale(1.0 / 7.0);
            near(layer.pressure_scattering()[0], (y0 - y1) / (y0 + y1), 2e-15);
        }
        assert!(matches!(
            layer.drive([C64::from_re(100.0), C64::ZERO], limits(), &mut || false),
            Err(ViscoelasticSlabError::OutsideSourceSupport)
        ));
        assert!(matches!(
            driven.field(0.1, &mut || true),
            Err(ViscoelasticSlabError::Cancelled)
        ));
        assert!(matches!(
            driven.absorption_between(0.0, h, budget(), &mut || true),
            Err(ViscoelasticSlabError::Cancelled)
        ));
    }
    let before = left.calls.get();
    let sine = det::sqrt((1.0 - 0.8) * (1.0 + 0.8));
    let cutoff_fluid = FluidAcousticState::try_new(1.4, 2.0 / sine).unwrap();
    let cutoff = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&left, &cutoff_fluid],
        20.0,
        0.8,
        0.2,
        0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(cutoff.exterior_normal_cosines()[1], C64::ZERO);
    assert!(cutoff.power_normalized_scattering().is_none());
    assert!(
        cutoff
            .drive([C64::ZERO, C64::from_re(0.01)], limits(), &mut || false)
            .is_err()
    );
    let cutoff_drive = cutoff
        .drive([C64::from_re(0.01), C64::ZERO], limits(), &mut || false)
        .unwrap();
    assert_eq!(cutoff_drive.power().outgoing_power_w_m2[1], 0.0);
    assert!(cutoff_drive.power().net_absorbed_power_w_m2.abs() < 1e-17);
    let _same = FluidViscoelasticLayer::try_new(
        source(&bulk, &shear),
        [&left, &left],
        20.0,
        0.8,
        0.2,
        0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(left.calls.get(), before + 2);
    assert!(matches!(
        FluidViscoelasticLayer::try_new(
            source(&bulk, &shear),
            [&left, &left],
            20.0,
            0.8,
            0.2,
            0,
            &mut || true
        ),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    let lossy = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let layer = FluidViscoelasticLayer::try_new(
        source(&lossy, &shear),
        [&left, &left],
        20.0,
        0.8,
        100.0,
        1,
        &mut || false,
    )
    .unwrap();
    let driven = layer
        .drive([C64::from_re(0.01), C64::ZERO], limits(), &mut || false)
        .unwrap();
    let too_small = ViscoelasticLayerWorkBudget {
        maximum_segments: 1,
        ..budget()
    };
    assert!(matches!(
        driven.absorption_between(0.0, 100.0, too_small, &mut || false),
        Err(ViscoelasticSlabError::Budget)
    ));
}
