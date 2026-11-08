use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::fluid::{
    FluidDriveLimits, FluidElasticChannel, FluidElasticInterface,
};
use fs_material::elastic::acoustic::viscoelastic::interface::FluidViscoelasticInterface;
use fs_material::elastic::acoustic::viscoelastic::{ViscoelasticSlabError, ViscoelasticSlabSource};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
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
fn limits() -> FluidDriveLimits {
    FluidDriveLimits::try_new(0.05, 0.05).unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}"
    );
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
fn g1_normal_incidence_matches_original_complex_impedance_without_requery() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    for omega in [10.0, 20.0, 80.0] {
        let fluid = CountedFluid {
            state: FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            calls: Cell::new(0),
        };
        let interface = FluidViscoelasticInterface::try_new(
            source(&bulk, &shear),
            &fluid,
            omega,
            1.0,
            2,
            &mut || false,
        )
        .unwrap();
        assert_eq!(fluid.calls.get(), 1);
        assert!(core::ptr::eq(interface.source().bulk, &bulk));
        assert!(core::ptr::eq(interface.source().shear, &shear));
        assert_eq!(interface.source().source_identity, ContentHash([9; 32]));
        let moduli = [&bulk, &shear].map(|law| {
            let (storage, loss) = law.modulus(omega);
            C64::new(storage, loss)
        });
        assert_eq!(interface.moduli_pa(), moduli);
        let longitudinal = moduli[0] + moduli[1].scale(4.0 / 3.0);
        let solid_impedance = longitudinal.scale(4.0).sqrt();
        let fluid_impedance = C64::from_re(20.0);
        let reflection = (solid_impedance - fluid_impedance) / (solid_impedance + fluid_impedance);
        near(interface.pressure_reflection(), reflection, 3e-15);
        let incident = C64::new(1.0, 0.2);
        let driven = interface.drive(incident, limits(), &mut || false).unwrap();
        let field = driven.solid_boundary_field();
        near(driven.reflected_pressure_pa(), incident * reflection, 3e-15);
        near(
            field.velocity_m_s[2],
            incident.scale(2.0) / (solid_impedance + fluid_impedance),
            3e-15,
        );
        near(field.traction_pa[2], -driven.fluid_pressure_pa(), 3e-15);
        near(
            field.volumetric_strain,
            -driven.fluid_pressure_pa() / longitudinal,
            3e-15,
        );
        assert_eq!(field.velocity_m_s[0], C64::ZERO);
        assert_eq!(field.traction_pa[0], C64::ZERO);
        let expected_bulk = 0.5 * omega * moduli[0].im * field.volumetric_strain.norm_sq();
        let expected_shear = omega * moduli[1].im * (2.0 / 3.0) * field.volumetric_strain.norm_sq();
        assert!((field.bulk_heat_density_w_m3 / expected_bulk - 1.0).abs() < 2e-14);
        assert!((field.shear_heat_density_w_m3 / expected_shear - 1.0).abs() < 2e-14);
        assert!(driven.power().balance_residual_w_m2.abs() < 2e-16);
        assert_eq!(
            driven.field(0.7, &mut || false).unwrap(),
            driven.field(0.7, &mut || false).unwrap()
        );
        assert_eq!(fluid.calls.get(), 1);
    }
}

#[test]
fn g0_outgoing_root_retains_weak_attenuation_when_squared_slowness_underflows() {
    let bulk = GeneralizedMaxwell::new(1e200, vec![(2.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(1e100, vec![]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let mut material = source(&bulk, &shear);
    material.density_kg_m3 = 1.0;
    let interface =
        FluidViscoelasticInterface::try_new(material, &fluid, 20.0, 1.0, 1, &mut || false).unwrap();
    // M*=1e200+i Pa to working precision: Im(rho/M) underflows, whereas
    // Im(q)=-5e-301 s/m remains representable and must not become lossless.
    let q = interface.normal_slownesses_s_m()[0];
    assert!((q.re / 1e-100 - 1.0).abs() < 5e-15);
    assert!((q.im / -5e-301 - 1.0).abs() < 5e-15);
}

#[test]
fn g3_lossless_oblique_limit_matches_existing_elastic_contact_through_critical_angles() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    let solid =
        ElasticAcousticState::from_bulk_and_shear(4.0, 60.0, 30.0, ContentHash([9; 32]), 0.1)
            .unwrap();
    let elastic = FluidElasticInterface::new(&fluid, &solid).unwrap();
    let incident_pressure = C64::from_re(0.01);
    for cosine in [1.0, 0.98, 0.9, 0.7, 0.4, 0.001] {
        let interface = FluidViscoelasticInterface::try_new(
            source(&bulk, &shear),
            &fluid,
            20.0,
            cosine,
            0,
            &mut || false,
        )
        .unwrap();
        let driven = interface
            .drive(incident_pressure, limits(), &mut || false)
            .unwrap();
        let old = elastic
            .prepare_incidence(FluidElasticChannel::Fluid, cosine)
            .unwrap()
            .drive(
                [C64::from_re(0.005), C64::ZERO, C64::ZERO, C64::ZERO],
                limits(),
            )
            .unwrap();
        // The elastic owner uses exp(-i omega t); conjugate its complete field.
        near(
            driven.fluid_pressure_pa(),
            old.fluid_pressure_pa.conj(),
            3e-14,
        );
        let field = driven.solid_boundary_field();
        for i in 0..3 {
            near(
                field.velocity_m_s[i],
                old.boundary_velocity_m_s[1][i].conj(),
                3e-14,
            );
            near(field.traction_pa[i], old.solid_traction_pa[i].conj(), 3e-14);
        }
        assert_eq!(field.bulk_heat_density_w_m3, 0.0);
        assert_eq!(field.shear_heat_density_w_m3, 0.0);
        let [p, s] = interface.normal_slownesses_s_m();
        assert!(p.im <= 0.0 && s.im <= 0.0);
        if cosine < 0.6 {
            assert!(p.im < 0.0 && s.im < 0.0);
            assert!((interface.pressure_reflection().norm_sq() - 1.0).abs() < 4e-14);
        }
        assert!(
            (driven.power().solid_boundary_power_w_m2 - old.boundary_normal_flux_w_m2[1]).abs()
                < 3e-16
        );
    }
}

#[test]
fn g1_oblique_reflection_matches_independent_solid_surface_impedance() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(20.0, 0.02)]).unwrap();
    let fluid = FluidAcousticState::try_new(1.0, 2.0).unwrap();
    for cosine in [1.0, 0.9, 0.65, 0.4] {
        let interface = FluidViscoelasticInterface::try_new(
            source(&bulk, &shear),
            &fluid,
            20.0,
            cosine,
            2,
            &mut || false,
        )
        .unwrap();
        // Eliminate the transmitted SV amplitude with sigma_xz=0, then
        // Z_surface=-sigma_zz/v_z = G²/(rho qP) [(qS²-p²)²+4p²qPqS].
        // This scalar impedance route is independent of the production LU.
        let [k, g] = interface.moduli_pa();
        let p2 = (1.0 - cosine * cosine) / 4.0;
        let qp = (C64::from_re(4.0) / (k + g.scale(4.0 / 3.0)) - C64::from_re(p2)).sqrt();
        let qs2 = C64::from_re(4.0) / g - C64::from_re(p2);
        let qs = qs2.sqrt();
        let a = qs2 - C64::from_re(p2);
        let impedance = (g * g).scale(0.25) * (a * a + (qp * qs).scale(4.0 * p2)) / qp;
        let scaled = impedance.scale(cosine);
        near(
            interface.pressure_reflection(),
            (scaled - C64::from_re(2.0)) / (scaled + C64::from_re(2.0)),
            2e-14,
        );
    }
}

#[test]
fn g1_total_ps_flux_closes_independent_bulk_shear_work_in_the_absorbing_half_space() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(40.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(20.0, 0.02)]).unwrap();
    for (speed, cosine) in [(4.0, 0.65), (2.0, 0.5)] {
        let fluid = FluidAcousticState::try_new(1.0, speed).unwrap();
        let interface = FluidViscoelasticInterface::try_new(
            source(&bulk, &shear),
            &fluid,
            20.0,
            cosine,
            2,
            &mut || false,
        )
        .unwrap();
        let driven = interface
            .drive(C64::new(0.05, 0.015), limits(), &mut || false)
            .unwrap();
        let boundary = driven.solid_boundary_field();
        near(boundary.traction_pa[0], C64::ZERO, 2e-15);
        near(boundary.traction_pa[2], -driven.fluid_pressure_pa(), 2e-15);
        near(
            boundary.velocity_m_s[2],
            driven.fluid_velocity_m_s()[2],
            2e-15,
        );
        assert!(boundary.velocity_m_s[0].abs() > 1e-5);
        let waves = interface.normal_slownesses_s_m();
        assert!(waves.iter().all(|q| q.re > 0.0 && q.im < 0.0));
        let decay = waves
            .iter()
            .map(|q| -20.0 * q.im)
            .fold(f64::INFINITY, f64::min);
        let step = 1e-4 / (20.0 * (waves[0].abs() + waves[1].abs()));
        for depth in [0.1 / decay, 0.7 / decay, 2.3 / decay] {
            let field = driven.field(depth, &mut || false).unwrap();
            let left = driven.field(depth - step, &mut || false).unwrap();
            let right = driven.field(depth + step, &mut || false).unwrap();
            let trace = field.strain_mandel[0] + field.strain_mandel[1] + field.strain_mandel[2];
            let deviator_norm2: f64 = field
                .strain_mandel
                .iter()
                .enumerate()
                .map(|(i, &e)| (if i < 3 { e - trace.scale(1.0 / 3.0) } else { e }).norm_sq())
                .sum();
            let [k, g] = interface.moduli_pa();
            let work = 10.0 * k.im * trace.norm_sq() + 20.0 * g.im * deviator_norm2;
            let reported = field.bulk_heat_density_w_m3 + field.shear_heat_density_w_m3;
            assert!((reported / work - 1.0).abs() < 5e-13);
            let flux_loss = (left.normal_power_w_m2 - right.normal_power_w_m2) / (2.0 * step);
            assert!(
                (flux_loss / work - 1.0).abs() < 2e-7,
                "flux loss={flux_loss}, work={work}"
            );
        }
        // Independent composite Simpson integral of the original constitutive
        // work. Both outgoing modes decay here; no power reaches infinity.
        let end = 14.0 / decay;
        let n = 4096;
        let h = end / n as f64;
        let mut integral = 0.0;
        for i in 0..=n {
            let field = driven.field(i as f64 * h, &mut || false).unwrap();
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * (field.bulk_heat_density_w_m3 + field.shear_heat_density_w_m3);
        }
        integral *= h / 3.0;
        assert!((integral / driven.power().solid_boundary_power_w_m2 - 1.0).abs() < 2e-7);
        assert!(driven.power().balance_residual_w_m2.abs() < 3e-16);
    }
}

#[test]
fn g0_g4_material_fluid_strain_frequency_and_cancellation_refusals_are_atomic() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = source(&bulk, &shear);
    for cosine in [0.0, -1.0, 1.1, f64::NAN] {
        assert!(
            FluidViscoelasticInterface::try_new(original, &fluid, 20.0, cosine, 2, &mut || false)
                .is_err()
        );
    }
    for omega in [0.0, 9.0, 101.0, f64::NAN] {
        assert!(matches!(
            FluidViscoelasticInterface::try_new(original, &fluid, omega, 0.8, 2, &mut || false),
            Err(ViscoelasticSlabError::OutsideSourceSupport)
        ));
    }
    assert!(matches!(
        FluidViscoelasticInterface::try_new(original, &fluid, 20.0, 0.8, 1, &mut || false),
        Err(ViscoelasticSlabError::Budget)
    ));
    assert!(matches!(
        FluidViscoelasticInterface::try_new(original, &fluid, 20.0, 0.8, 2, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    for density in [0.0, -1.0, f64::NAN] {
        let mut bad = original;
        bad.density_kg_m3 = density;
        assert!(
            FluidViscoelasticInterface::try_new(bad, &fluid, 20.0, 0.8, 2, &mut || false).is_err()
        );
    }
    let interface =
        FluidViscoelasticInterface::try_new(original, &fluid, 20.0, 0.8, 2, &mut || false).unwrap();
    assert!(matches!(
        interface.drive(C64::from_re(1000.0), limits(), &mut || false),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    assert!(matches!(
        interface.drive(
            C64::ONE,
            FluidDriveLimits::try_new(1e-10, 0.05).unwrap(),
            &mut || false
        ),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    let mut calls = 0;
    assert!(matches!(
        interface.drive(C64::ONE, limits(), &mut || {
            calls += 1;
            calls == 2
        }),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    let good = interface.drive(C64::ONE, limits(), &mut || false).unwrap();
    assert!(matches!(
        good.field(0.1, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(good.field(-0.1, &mut || false).is_err());
    assert!(good.field(f64::NAN, &mut || false).is_err());
    assert!(good.field(0.1, &mut || false).is_ok());
    let zero = interface.drive(C64::ZERO, limits(), &mut || false).unwrap();
    assert_eq!(zero.power().incident_power_w_m2, 0.0);
    assert_eq!(
        zero.field(0.1, &mut || false)
            .unwrap()
            .bulk_heat_density_w_m3,
        0.0
    );
}
