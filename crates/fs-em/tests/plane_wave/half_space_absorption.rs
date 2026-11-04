use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumError, AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide,
        PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_math::det;

fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn relative(a: f64, b: f64, tolerance: f64) {
    near(a, b, tolerance * a.abs().max(b.abs()).max(1e-300));
}
fn conductor(epsilon: f64, mu: f64, sigma: f64) -> ResolvedConductorState {
    let basis = ElectricalPropertyBasis::Conductivity;
    ResolvedConductorState::try_from_material_with_basis(
        source(
            &[
                (EPS_NAME, EPS, epsilon),
                (MU_NAME, MU, mu),
                (basis.property_name(), basis.dimensions(), sigma),
            ],
            ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            QueryPoint::new().with("T", 300.0).unwrap(),
        ),
        basis,
    )
    .unwrap()
}

#[test]
fn coherent_substrate_heat_closes_volume_poynting_and_film_power_balances() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, 1.2 * E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Debye(&lossy),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: 0.23,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [2.0, 3.0],
    };
    let mut input = vec![[C64::ZERO; 2]; grid.width * grid.height];
    input[2] = [C64::new(0.9, 0.3), C64::new(0.2, -0.7)];
    input[21] = [C64::new(-0.4, 0.5), C64::new(0.8, 0.2)];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let begin = 0.1;
    let end = 0.9;
    let intervals = 128;
    let mut integrated = 0.0;
    for j in 0..=intervals {
        let depth = begin + (end - begin) * j as f64 / intervals as f64;
        let plane = beam
            .transmitted_absorption_field(depth, limits(), || false)
            .unwrap();
        for (e, &heat) in plane
            .field()
            .electric_v_per_m()
            .iter()
            .zip(plane.absorbed_power_density_w_m3())
        {
            let expected = -0.5
                * 1e9
                * stack.exit_permittivity_f_m().im
                * e.iter().map(|v| v.norm_sq()).sum::<f64>();
            relative(heat, expected, 5e-15);
        }
        let weight = if j == 0 || j == intervals {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integrated += weight * plane.absorbed_power_per_depth_w_m();
    }
    integrated *= (end - begin) / (3.0 * intervals as f64);
    let exact = beam
        .transmitted_absorbed_power_w(begin, end, || false)
        .unwrap();
    relative(integrated, exact, 3e-10);
    let entering = beam
        .field(AngularSpectrumSide::Transmitted, begin, || false)
        .unwrap();
    let leaving = beam
        .field(AngularSpectrumSide::Transmitted, end, || false)
        .unwrap();
    relative(
        entering.normal_power_w() - leaving.normal_power_w(),
        exact,
        1e-12,
    );
    let all_substrate = beam
        .transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
        .unwrap();
    relative(all_substrate, beam.power().transmitted_w, 2e-15);
    relative(
        beam.power().reflected_w + beam.power().film_absorbed_w + all_substrate,
        beam.power().incident_w,
        1e-12,
    );
    let surface = beam
        .transmitted_absorption_field(0.0, limits(), || false)
        .unwrap();
    let maximum = surface
        .absorbed_power_density_w_m3()
        .iter()
        .copied()
        .fold(0.0, f64::max);
    let minimum = surface
        .absorbed_power_density_w_m3()
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    assert!(
        maximum > 100.0 * minimum,
        "coherent spatial absorption must retain its fringes"
    );
}

#[test]
fn oblique_ohmic_substrate_deposits_joule_heat_with_normal_skin_decay() {
    let incident = material(E, M);
    let sigma = 0.02;
    let metal = conductor(3.0 * E, 1.4 * M, sigma);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Ohmic(&metal), 1e9, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [2.0, 3.0],
    };
    for polarization in [0, 1] {
        let mut input = [[C64::ZERO; 2]; 2];
        input[0][polarization] = C64::new(0.8, -0.3);
        input[1][polarization] = -input[0][polarization];
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
        let kx = core::f64::consts::PI / grid.spacing_m[0];
        let epsilon = C64::new(3.0 * E, -sigma / 1e9);
        let kz = (epsilon.scale(1e18 * 1.4 * M) - C64::from_re(kx * kx)).sqrt();
        let depth = 0.6;
        let plane = beam
            .transmitted_absorption_field(depth, limits(), || false)
            .unwrap();
        for (e, &heat) in plane
            .field()
            .electric_v_per_m()
            .iter()
            .zip(plane.absorbed_power_density_w_m3())
        {
            relative(
                heat,
                0.5 * sigma * e.iter().map(|v| v.norm_sq()).sum::<f64>(),
                4e-15,
            );
        }
        relative(
            plane.absorbed_power_per_depth_w_m(),
            -2.0 * kz.im * beam.power().transmitted_w * det::exp(2.0 * kz.im * depth),
            4e-13,
        );
        let tail = beam
            .transmitted_absorbed_power_w(depth, f64::INFINITY, || false)
            .unwrap();
        relative(tail, plane.field().normal_power_w(), 4e-13);
    }
}

#[test]
fn weak_absorption_is_positive_below_endpoint_flux_roundoff_and_opaque_tail_is_zero() {
    let incident = material(1.0, 1.0);
    let sigma = 1e-14;
    let weak = conductor(1.0, 1.0, sigma);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Ohmic(&weak), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0, 1.0],
    };
    let beam =
        PlanarAngularSpectrum::prepare(&stack, grid, &[[C64::ONE, C64::ZERO]], limits(), || false)
            .unwrap();
    let width = 1e-6;
    let absorbed = beam
        .transmitted_absorbed_power_w(0.0, width, || false)
        .unwrap();
    assert!(absorbed > 0.0);
    relative(absorbed, 0.5 * sigma * width, 1e-13);
    let before = beam
        .field(AngularSpectrumSide::Transmitted, 0.0, || false)
        .unwrap();
    let after = beam
        .field(AngularSpectrumSide::Transmitted, width, || false)
        .unwrap();
    relative(before.normal_power_w(), after.normal_power_w(), 5e-15);
    assert_eq!(
        beam.transmitted_absorbed_power_w(0.0, 0.0, || false)
            .unwrap(),
        0.0
    );
    assert_eq!(
        beam.transmitted_absorbed_power_w(1e300, f64::INFINITY, || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn evanescent_lossless_substrate_has_zero_heat_despite_field_decay() {
    let incident = material(4.0, 1.0);
    let exit = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [2.0, 1.0],
    };
    let input = [[C64::ONE, C64::ONE], [-C64::ONE, -C64::ONE]];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let surface = beam
        .transmitted_absorption_field(0.0, limits(), || false)
        .unwrap();
    let deeper = beam
        .transmitted_absorption_field(1.0, limits(), || false)
        .unwrap();
    assert!(
        deeper.field().electric_v_per_m()[0][0].abs()
            < surface.field().electric_v_per_m()[0][0].abs()
    );
    assert_eq!(surface.absorbed_power_density_w_m3(), &[0.0, 0.0]);
    assert_eq!(deeper.absorbed_power_per_depth_w_m(), 0.0);
    assert_eq!(
        beam.transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn substrate_heat_admits_budgets_and_depths_before_work_and_obeys_cancellation() {
    let incident = material(1.0, 1.0);
    let exit = conductor(1.0, 1.0, 0.2);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Ohmic(&exit), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0, 1.0],
    };
    let beam =
        PlanarAngularSpectrum::prepare(&stack, grid, &[[C64::ONE, C64::ZERO]], limits(), || false)
            .unwrap();
    let mut admitted = limits();
    admitted.maximum_workspace_bytes = beam.transmitted_absorption_workspace_bytes().unwrap();
    assert!(
        beam.transmitted_absorption_field(0.0, admitted, || false)
            .is_ok()
    );
    admitted.maximum_workspace_bytes -= 1;
    assert!(
        beam.transmitted_absorption_field(0.0, admitted, || false)
            .is_err()
    );
    for (start, end) in [
        (f64::NAN, 1.0),
        (0.0, f64::NAN),
        (-1.0, 0.0),
        (1.0, 0.0),
        (f64::INFINITY, f64::INFINITY),
        (0.0, f64::NEG_INFINITY),
    ] {
        assert!(
            beam.transmitted_absorbed_power_w(start, end, || false)
                .is_err()
        );
    }
    assert!(matches!(
        beam.transmitted_absorbed_power_w(0.0, 1.0, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        beam.transmitted_absorption_field(0.0, limits(), || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(
        beam.transmitted_absorption_field(0.0, limits(), || false)
            .is_ok()
    );
}
