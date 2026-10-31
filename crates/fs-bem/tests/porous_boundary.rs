//! Material-to-BEM reflection, phase and passive-work regressions.
//! These tests exercise the existing JCA(L) and Robin owners; no new absorber
//! fit or BEM accuracy claim is introduced.
use fs_bem::helmholtz::bc::{self, Admittance, PlaneWave, PorousRigidBoundary};
use fs_bem::helmholtz::Formulation;
use fs_bem::panel3d::SpherePanels;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(density: f64, sound_speed: f64, viscosity: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(density, sound_speed).unwrap(),
        FluidTransportState::try_new(density, 1005.0, 0.0257, viscosity).unwrap(),
        1.4,
    )
    .unwrap()
}

fn law() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap()
}

#[test]
fn material_reflection_preserves_phase_and_independent_known_answers() {
    let law = law();
    let fluid = fluid(1.2, 340.0, 1.81e-5);
    // Independent 75-digit Decimal evaluation of the published JCAL equations.
    // Coth was evaluated from complex exponentials, rather than the material
    // owner's tanh/tan quotient. Synthetic fluid tuple above; 60 mm backing.
    // Values use the BEM exp(-i omega t) convention.
    for (frequency, real, imaginary, absorption) in [
        (150.0, 0.8438535316431461, 0.41061635992745694, 0.119305422093315),
        (500.0, 0.1499064065842737, 0.5847212235088352, 0.6356291600433212),
        (1000.0, -0.18212577259495563, 0.09485375838418858, 0.9578329674770845),
        (1600.0, 0.14566508996865872, -0.18087374506259907, 0.9460663699114524),
    ] {
        let boundary = PorousRigidBoundary::new(&law, &fluid, TAU * frequency, 0.06).unwrap();
        let reflection = boundary.plane_wave_reflection(1.0).unwrap();
        assert!((reflection - C64::new(real, imaginary)).abs() < 1e-10);
        assert!((1.0 - reflection.norm_sq() - absorption).abs() < 1e-10);
        // A dropped conjugation preserves absorption but reverses this phase.
        assert!((reflection - C64::new(real, -imaginary)).abs() > 0.1);
        assert!(std::ptr::eq(boundary.law(), &law));
        assert_eq!(*boundary.response().fluid(), fluid);
        assert_eq!(boundary.thickness_m(), 0.06);
        assert_eq!(boundary.medium().density, 1.2);
        assert_eq!(boundary.medium().sound_speed, 340.0);
        assert_eq!(boundary.wavenumber_rad_m(), TAU * frequency / 340.0);
    }
}

#[test]
fn oblique_reflection_closes_normal_power_without_losing_phase() {
    let law = law();
    let fluid = fluid(1.2, 340.0, 1.81e-5);
    for frequency in [150.0, 500.0, 1000.0, 1600.0] {
        let boundary = PorousRigidBoundary::new(&law, &fluid, TAU * frequency, 0.06).unwrap();
        let zs = boundary.response().surface_impedance_rigid(0.06).unwrap().conj();
        let z0 = 1.2 * 340.0;
        for cosine in [1.0, 0.5, 0.1, 1e-6] {
            let r = boundary.plane_wave_reflection(cosine).unwrap();
            // Independent surface-work calculation: |p_total|^2 Re(Y)/2
            // divided by incident normal intensity cosine/(2 rho c).
            let work_fraction = (C64::ONE + r).norm_sq()
                * zs.re / zs.norm_sq() * z0 / cosine;
            let reflected_fraction = r.norm_sq();
            assert!(work_fraction >= 0.0 && work_fraction <= 1.0 + 1e-12);
            assert!((work_fraction + reflected_fraction - 1.0).abs() < 2e-10);
        }
        assert!((boundary.plane_wave_reflection(1e-6).unwrap() + C64::ONE).abs() < 2e-5);
    }
}

#[test]
fn material_changes_move_both_boundary_response_and_exterior_medium() {
    let law = law();
    let air = fluid(1.2, 340.0, 1.81e-5);
    // Manufactured liquid-like state, not a measured qualification of this foam.
    let liquid = fluid(1000.0, 1400.0, 1e-3);
    let first = PorousRigidBoundary::new(&law, &air, TAU * 500.0, 0.06).unwrap();
    let other = PorousRigidBoundary::new(&law, &liquid, TAU * 500.0, 0.06).unwrap();
    assert_eq!(other.medium().density, 1000.0);
    assert_eq!(other.medium().sound_speed, 1400.0);
    assert_eq!(*other.response().fluid(), liquid);
    assert_eq!(other.wavenumber_rad_m(), TAU * 500.0 / 1400.0);
    assert!((other.normalized_admittance() - first.normalized_admittance()).abs() > 1e-3);
    let thin = PorousRigidBoundary::new(&law, &air, TAU * 500.0, 0.002).unwrap();
    assert!((thin.plane_wave_reflection(1.0).unwrap()
        - first.plane_wave_reflection(1.0).unwrap()).abs() > 0.1);
}

#[test]
fn porous_scattering_uses_original_robin_rows_and_total_wall_work() {
    let law = law();
    let fluid = fluid(1.2, 340.0, 1.81e-5);
    let boundary = PorousRigidBoundary::new(&law, &fluid, TAU * 150.0, 0.03).unwrap();
    let surface = SpherePanels::icosphere(0.5, 1).unwrap();
    let incident = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.8, 0.2),
    };
    let solved = boundary.solve_scattering(&surface, incident, Formulation::BurtonMiller).unwrap();
    // Construct beta independently from real/imaginary specific impedance.
    let z = law.surface_impedance_rigid_fluid(&fluid, TAU * 150.0, 0.03).unwrap();
    let denominator = z.re * z.re + z.im * z.im;
    let beta = C64::new(408.0 * z.re / denominator, 408.0 * z.im / denominator);
    let coefficients = vec![beta; surface.areas().len()];
    let reference = bc::solve_scattering_bc(
        &surface, boundary.wavenumber_rad_m(), boundary.medium(),
        Admittance::Local { beta: &coefficients, declared_active: false },
        incident, Formulation::BurtonMiller,
    ).unwrap();
    for (&actual, &expected) in solved.scattered.pressure.iter().zip(&reference.scattered.pressure) {
        assert!((actual - expected).abs() < 1e-10 * expected.abs().max(1.0));
    }
    let surface_work: f64 = solved.scattered.pressure.iter()
        .zip(&solved.incident_pressure_pa).zip(surface.areas())
        .map(|((&scattered, &incident), &area)| {
            0.5 * (scattered + incident).norm_sq() * area * z.re / denominator
        }).sum();
    assert!(surface_work > 0.0);
    assert!((solved.absorbed_power_w - surface_work).abs() < 1e-12 * surface_work);
    assert!(!solved.declared_active);

    let drive = vec![C64::from_re(0.01); surface.areas().len()];
    let radiated = boundary.solve_radiation(&surface, &drive, Formulation::BurtonMiller).unwrap();
    for ((&pressure, &velocity), &driving) in radiated.radiation.pressure.iter()
        .zip(&radiated.radiation.velocity).zip(&drive) {
        let absorbed_velocity = pressure / z.conj();
        assert!((velocity + absorbed_velocity - driving).abs() < 1e-12);
    }
    assert!(radiated.absorbed_power_w > 0.0);
    let balance = radiated.driving_power_w
        - radiated.absorbed_power_w - radiated.radiation.radiated_power;
    assert!(balance.abs() < 1e-10 * radiated.driving_power_w.abs().max(1e-30));
}

#[test]
fn unsupported_sources_geometry_and_incidence_refuse() {
    let law = law();
    let fluid = fluid(1.2, 340.0, 1.81e-5);
    for frequency in [0.0, 149.0, 1601.0, f64::NAN, f64::INFINITY] {
        assert!(PorousRigidBoundary::new(&law, &fluid, TAU * frequency, 0.06).is_err());
    }
    for thickness in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        assert!(PorousRigidBoundary::new(&law, &fluid, TAU * 500.0, thickness).is_err());
    }
    let boundary = PorousRigidBoundary::new(&law, &fluid, TAU * 500.0, 0.06).unwrap();
    for cosine in [0.0, -0.1, 1.001, f64::NAN, f64::INFINITY] {
        assert!(boundary.plane_wave_reflection(cosine).is_err());
    }
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    assert!(boundary.solve_radiation(&surface, &[], Formulation::BurtonMiller).is_err());
}
