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

#[test]
fn layer_boundary_reduces_to_single_porous_coating_in_its_original_fluid() {
    use fs_bem::helmholtz::bc::PorousLayersBoundary;
    use fs_material::porous::layers::{Layer, RigidBackedLayers};

    let law = law();
    for fluid in [
        fluid(1.2, 340.0, 1.81e-5),
        fluid(1000.0, 1400.0, 1e-3),
    ] {
        let layers = [Layer::Porous { law: &law, thickness_m: 0.06 }];
        let stack = RigidBackedLayers::new(&fluid, TAU * 500.0, &layers).unwrap();
        let layered = PorousLayersBoundary::new(&stack).unwrap();
        let single = PorousRigidBoundary::new(&law, &fluid, TAU * 500.0, 0.06).unwrap();
        let scale = single.normalized_admittance().abs().max(1.0);
        assert!((layered.normalized_admittance() - single.normalized_admittance()).abs()
            < 1e-12 * scale);
        assert!(std::ptr::eq(layered.layers(), &stack));
        assert_eq!(*layered.layers().fluid(), fluid);
        assert_eq!(layered.medium(), single.medium());
        assert_eq!(layered.wavenumber_rad_m(), single.wavenumber_rad_m());
        let beta = layered.normalized_admittance();
        let reflection = (C64::ONE - beta) / (C64::ONE + beta);
        assert!((reflection - stack.reflection(1.0).unwrap().pressure_reflection.conj()).abs()
            < 1e-12);
    }
}

#[test]
fn porous_gap_stack_drives_dense_and_iterative_robin_owners_with_physical_work() {
    use fs_bem::helmholtz::bc::PorousLayersBoundary;
    use fs_bem::helmholtz::{IterativeOptions, Preconditioner};
    use fs_exec::CancelGate;
    use fs_material::porous::layers::{Layer, RigidBackedLayers};

    let law = law();
    let fluid = fluid(1.2, 340.0, 1.81e-5);
    let layers = [
        Layer::Porous { law: &law, thickness_m: 0.03 },
        Layer::FluidGap { thickness_m: 0.05 },
    ];
    let stack = RigidBackedLayers::new(&fluid, TAU * 150.0, &layers).unwrap();
    let boundary = PorousLayersBoundary::new(&stack).unwrap();
    let normal = stack.reflection(1.0).unwrap();
    let beta = boundary.normalized_admittance();
    // Independent conventional cos/sin transfer matrices for the porous
    // slab and fluid gap, evaluated in the BEM exp(-i omega t) convention.
    let expected_beta = C64::new(0.027714988312696837, -0.25759369018286377);
    let expected_reflection = C64::new(0.831032207187226, 0.45894275013679325);
    assert!((beta - expected_beta).abs() < 1e-12);
    assert!((beta - normal.normalized_surface_admittance().unwrap().conj()).abs() < 1e-13);
    assert!(beta.re > 0.0);
    assert!((beta - normal.normalized_surface_admittance().unwrap()).abs() > 0.01);
    let reflected = (C64::ONE - beta) / (C64::ONE + beta);
    assert!((reflected - expected_reflection).abs() < 1e-12);
    assert!((reflected - normal.pressure_reflection.conj()).abs() < 1e-12);
    assert!((1.0 - reflected.norm_sq() - normal.absorbed_fraction).abs() < 1e-12);

    let surface = SpherePanels::icosphere(0.5, 1).unwrap();
    let incident = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.8, 0.2),
    };
    let solved = boundary.solve_scattering(&surface, incident, Formulation::BurtonMiller).unwrap();
    let surface_work: f64 = solved.scattered.pressure.iter()
        .zip(&solved.incident_pressure_pa).zip(surface.areas())
        .map(|((&scattered, &incident), &area)| {
            0.5 * (scattered + incident).norm_sq() * area * beta.re / 408.0
        }).sum();
    assert!(surface_work > 0.0);
    assert!((solved.absorbed_power_w - surface_work).abs() < 1e-12 * surface_work);
    assert!(!solved.declared_active);

    let drive = vec![C64::from_re(0.01); surface.areas().len()];
    let radiation = boundary.solve_radiation(&surface, &drive, Formulation::BurtonMiller).unwrap();
    for ((&pressure, &velocity), &driving) in radiation.radiation.pressure.iter()
        .zip(&radiation.radiation.velocity).zip(&drive) {
        assert!((velocity + (beta * pressure).scale(1.0 / 408.0) - driving).abs() < 1e-12);
    }
    let balance = radiation.driving_power_w
        - radiation.absorbed_power_w - radiation.radiation.radiated_power;
    assert!(radiation.absorbed_power_w > 0.0);
    assert!(balance.abs() < 1e-10 * radiation.driving_power_w.abs().max(1e-30));

    // One FMM leaf evaluates every pair directly; this isolates the adapter's
    // iterative Robin integration from interpolation-order error.
    let options = IterativeOptions {
        tol: 1e-10,
        restart: 40,
        max_iters: 100,
        leaf_capacity: 128,
        precond: Preconditioner::NearFieldBlockJacobi { max_block_panels: 80 },
        ..IterativeOptions::default()
    };
    let (iterative, report, _) = boundary.solve_scattering_iterative(
        &surface, incident, Formulation::BurtonMiller, &options,
        &CancelGate::new_clock_free(),
    ).unwrap();
    assert!(report.rel_residual < 1e-9);
    let error: f64 = iterative.scattered.pressure.iter().zip(&solved.scattered.pressure)
        .map(|(&actual, &expected)| (actual - expected).norm_sq()).sum();
    let norm: f64 = solved.scattered.pressure.iter().map(|p| p.norm_sq()).sum();
    assert!(error < 1e-14 * norm);
    assert!((iterative.absorbed_power_w - surface_work).abs() < 1e-7 * surface_work);
}

#[test]
fn empty_layer_boundary_is_rigid_and_keeps_iterative_cancellation() {
    use fs_bem::helmholtz::bc::PorousLayersBoundary;
    use fs_bem::helmholtz::{IterativeError, IterativeOptions};
    use fs_exec::CancelGate;
    use fs_material::porous::layers::RigidBackedLayers;

    let fluid = fluid(1.7, 290.0, 2.1e-5);
    let stack = RigidBackedLayers::new(&fluid, TAU * 150.0, &[]).unwrap();
    let boundary = PorousLayersBoundary::new(&stack).unwrap();
    assert_eq!(boundary.normalized_admittance(), C64::ZERO);
    assert_eq!(boundary.medium().density, 1.7);
    assert_eq!(boundary.medium().sound_speed, 290.0);
    assert_eq!(boundary.wavenumber_rad_m(), TAU * 150.0 / 290.0);

    let surface = SpherePanels::icosphere(0.5, 1).unwrap();
    let incident = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    };
    let solved = boundary.solve_scattering(&surface, incident, Formulation::BurtonMiller).unwrap();
    let rigid = bc::solve_scattering_bc(
        &surface, boundary.wavenumber_rad_m(), boundary.medium(),
        Admittance::Rigid, incident, Formulation::BurtonMiller,
    ).unwrap();
    assert_eq!(solved.scattered.pressure, rigid.scattered.pressure);
    assert_eq!(solved.absorbed_power_w, 0.0);

    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        boundary.solve_scattering_iterative(
            &surface, incident, Formulation::BurtonMiller, &IterativeOptions::default(), &gate,
        ),
        Err(IterativeError::Cancelled { matvecs: 0 })
    ));
}
