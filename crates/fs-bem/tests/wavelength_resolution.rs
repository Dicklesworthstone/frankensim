//! G0/G3: spatial phase must be resolved along every retained panel direction.
use fs_bem::helmholtz::iterative::{IterativeError, IterativeOptions, IterativeRadiation};
use fs_bem::helmholtz::laplace::complex::{self, ComplexError, ComplexOptions};
use fs_bem::helmholtz::molecular::{MolecularError, MolecularExterior};
use fs_bem::helmholtz::transmission::{FluidInterface, Options, TransmissionError};
use fs_bem::helmholtz::{Formulation, HelmholtzError, Medium, solve_radiation};
use fs_bem::panel3d::SpherePanels;
use fs_bem::reciprocal_radiation::open_screen::{OpenOptions, OpenPanels, PreparedScreen};
use fs_bem::reciprocal_radiation::{Budget, Prepared};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn slender_tetrahedron(scale: f64) -> SpherePanels {
    let vertices = [[0., 0., 0.], [2., 0., 0.], [0., 0.001, 0.], [0., 0., 0.001]];
    let faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
    SpherePanels::from_triangles(
        faces
            .map(|face| face.map(|i| vertices[i].map(|x| scale * x)))
            .to_vec(),
    )
    .unwrap()
}

fn complex_options(n: usize) -> ComplexOptions {
    ComplexOptions {
        max_panels: n,
        max_bytes: 64 * n * n,
        max_work: 512 * (n * n) as u64,
    }
}

#[test]
fn wavelength_resolution_rejects_long_thin_panels_in_dense_and_iterative_paths() {
    let gate = CancelGate::new_clock_free();
    for scale in [0.01, 1.0, 100.0] {
        let surface = slender_tetrahedron(scale);
        let k = PI / scale;
        let old_size = surface
            .areas()
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
            .sqrt();
        // The former area rule sees > 50 panels per wavelength. In fact each
        // long edge spans a whole wavelength, so no constant trace resolves it.
        assert!(TAU / (k * old_size) > 50.0);
        assert!((surface.resolution_length() / scale - 4.000001_f64.sqrt()).abs() < 1e-14);
        let velocity = vec![C64::ONE; surface.areas().len()];
        let expected = TAU / (k * surface.resolution_length());
        for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
            let error =
                solve_radiation(&surface, k, Medium::air(), &velocity, formulation).unwrap_err();
            assert!(
                matches!(error, HelmholtzError::TooCoarse { panels_per_wavelength }
                if (panels_per_wavelength - expected).abs() < 1e-14)
            );
            assert!(matches!(
                IterativeRadiation::prepare(&surface, k, Medium::air(), formulation,
                    &IterativeOptions::default(), &gate),
                Err(IterativeError::Helmholtz(HelmholtzError::TooCoarse { panels_per_wavelength }))
                    if (panels_per_wavelength - expected).abs() < 1e-14
            ));
        }
    }
}

#[test]
fn wavelength_resolution_limits_attenuation_and_molecular_propagation_too() {
    let surface = slender_tetrahedron(1.0);
    let options = complex_options(surface.areas().len());
    // Phase alone is resolved; the short attenuation length is not. The
    // existing |k| criterion must use the true geometric span as well.
    assert!(TAU / (0.2 * surface.resolution_length()) > 6.0);
    assert!(matches!(
        complex::assemble_with_cancel(
            &surface,
            C64::new(0.2, PI),
            Formulation::PlainCbie,
            &options,
            || false
        ),
        Err(ComplexError::Boundary(HelmholtzError::TooCoarse { .. }))
    ));
    let law = read_relaxation(&mut "declared 1 synthetic 0".split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let omega = PI * gas.sound_speed;
    let source = law
        .evaluate(&RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap())
        .unwrap();
    let exterior = MolecularExterior::new(source, [0.5 * omega, 2.0 * omega]).unwrap();
    let velocity = vec![C64::ONE; surface.areas().len()];
    assert!(matches!(
        exterior.solve_batch(
            &surface,
            omega,
            &[&velocity],
            Formulation::PlainCbie,
            &options,
            || false
        ),
        Err(MolecularError::Operator(ComplexError::Boundary(
            HelmholtzError::TooCoarse { .. }
        )))
    ));
}

#[test]
fn wavelength_resolution_checks_the_shorter_material_wavelength_at_an_interface() {
    let surface = slender_tetrahedron(1.0);
    let fast = FluidAcousticState::try_new(1.0, 10.0).unwrap();
    let slow = FluidAcousticState::try_new(1.0, 1.0).unwrap();
    assert!(TAU / ((PI / 10.0) * surface.resolution_length()) > 6.0);
    let gate = CancelGate::new_clock_free();
    for (outside, inside) in [(&fast, &slow), (&slow, &fast)] {
        assert!(matches!(
            FluidInterface::new(&surface, outside, inside, PI, Options::default(), &gate),
            Err(TransmissionError::Helmholtz(
                HelmholtzError::TooCoarse { .. }
            ))
        ));
    }
}

#[test]
fn wavelength_resolution_is_shared_by_reciprocal_closed_and_open_surfaces() {
    let closed = slender_tetrahedron(1.0);
    assert!(matches!(
        Prepared::new(
            &closed,
            PI,
            Budget {
                maximum_unknowns: 32,
                maximum_cache_bytes: 1_000_000,
                maximum_work: 1_000_000,
                angular_refinement: 1
            },
            || false
        ),
        Err(HelmholtzError::TooCoarse { .. })
    ));
    let open = OpenPanels::new(
        vec![
            [0., 0., 0.],
            [-1., -0.001, 0.],
            [1., -0.001, 0.],
            [1., 0.001, 0.],
            [-1., 0.001, 0.],
        ],
        vec![[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
    )
    .unwrap();
    assert_eq!(open.resolution_length(), 2.0);
    assert!(matches!(
        PreparedScreen::new(
            &open,
            PI,
            OpenOptions {
                edge_level: 0,
                maximum_edge_level: 0,
                triangle_order: 4,
                angular_refinement: 1,
                maximum_unknowns: 32,
                maximum_cache_bytes: 1_000_000,
                maximum_work: 1_000_000,
                workers: 1
            },
            || false
        ),
        Err(HelmholtzError::TooCoarse { .. })
    ));
}

#[test]
fn wavelength_resolution_preserves_resolved_solutions_and_the_exact_boundary() {
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let legacy = SpherePanels::new(
        surface.centroids().to_vec(),
        surface.normals().to_vec(),
        surface.areas().to_vec(),
    )
    .unwrap();
    let expected_legacy = surface
        .areas()
        .iter()
        .copied()
        .fold(0.0_f64, f64::max)
        .sqrt();
    assert_eq!(
        legacy.resolution_length().to_bits(),
        expected_legacy.to_bits()
    );
    let velocity = vec![C64::ONE; surface.areas().len()];
    let solve = |panels: &SpherePanels, k| {
        solve_radiation(
            panels,
            k,
            Medium::air(),
            &velocity,
            Formulation::BurtonMiller,
        )
    };
    // The default Burton-Miller arm uses identical centroid/equivalent-disc
    // kernels with or without retained triangles. Plain CBIE deliberately
    // uses a different triangle-integrated operator and is not compared here.
    let triangle = solve(&surface, 0.5).unwrap();
    let centroid = solve(&legacy, 0.5).unwrap();
    assert_eq!(triangle.pressure, centroid.pressure);
    assert_eq!(
        triangle.radiated_power.to_bits(),
        centroid.radiated_power.to_bits()
    );
    let boundary = TAU / (6.0 * surface.resolution_length());
    assert!(solve(&surface, boundary).is_ok());
    assert!(matches!(
        solve(&surface, boundary * (1.0 + 1e-12)),
        Err(HelmholtzError::TooCoarse { .. })
    ));
}
