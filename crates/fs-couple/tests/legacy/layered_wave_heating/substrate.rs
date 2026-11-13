use super::*;
use fs_couple::layered_wave_heating::substrate::{
    PreparedSubstrateWaveHeating, SubstrateWaveHeatingLimits,
};
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumGrid,
    AngularSpectrumLimits, PlanarAngularSpectrum,
};

fn field_limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 64,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn source_limits() -> AngularSpectrumCoherenceLimits {
    AngularSpectrumCoherenceLimits {
        maximum_modes: 8,
        maximum_point_mode_pairs: 512,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn heat_limits() -> SubstrateWaveHeatingLimits {
    SubstrateWaveHeatingLimits {
        maximum_cells: 128,
        maximum_mode_cell_pairs: 4096,
        maximum_workspace_bytes: 1 << 20,
    }
}
fn grid(spacing_x: f64) -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [spacing_x, AREA / (2.0 * spacing_x)],
    }
}
fn beam<'s, 'm>(
    stack: &'s PlanarDielectricLayers<'m>,
    grid: AngularSpectrumGrid,
    amplitude: f64,
    component: usize,
) -> PlanarAngularSpectrum<'s, 'm> {
    let mut input = [[C64::ZERO; 2]; 2];
    input[0][component] = C64::from_re(amplitude);
    input[1][component] = C64::from_re(-amplitude);
    PlanarAngularSpectrum::prepare(stack, grid, &input, field_limits(), || false).unwrap()
}

/// G1: independent oblique Fresnel amplitudes and exponential volume integrals
/// feed the original constant-capacity/conductivity finite-volume equations.
#[test]
fn substrate_cells_drive_native_heat_with_oblique_maxwell_power_once() {
    let exterior = state(&card("substrate incidence", 1.0, 0.0), false, 100_000.0);
    let card = card("substrate", 3.0, 0.7);
    let model = debye(&card, 100_000.0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior, DielectricResponse::Debye(&model), OMEGA, &[],
    ).unwrap();
    let width = 0.4 / k();
    let start = 0.2 / k();
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let old_cells = initial.cells().clone();
    let amplitude = 2e4;
    let dt = 0.02;
    let grid = grid(2.0);
    let kx = core::f64::consts::PI / grid.spacing_m[0];
    let epsilon = C64::new(3.0 * E, -0.7 * E);
    let kz = (epsilon.scale(OMEGA * OMEGA * M) - C64::from_re(kx * kx)).sqrt();
    let incident_kz = (OMEGA * OMEGA * E * M - kx * kx).sqrt();
    for component in [0, 1] {
        let beam = beam(&stack, grid, amplitude, component);
        let heat = PreparedSubstrateWaveHeating::try_new(
            &bound, &beam, start, heat_limits(), || false,
        ).unwrap();
        assert!(std::ptr::eq(heat.source(), &beam));
        assert!(std::ptr::eq(heat.thermal_body(), &bound));
        let (incident_y, exit_y) = if component == 0 {
            (OMEGA * E / incident_kz, epsilon.scale(OMEGA) / kz)
        } else {
            (incident_kz / (OMEGA * M), kz.scale(1.0 / (OMEGA * M)))
        };
        let transmission = C64::from_re(2.0 * incident_y)
            / (C64::from_re(incident_y) + exit_y);
        let interface_w = 0.5 * AREA * amplitude * amplitude
            * transmission.norm_sq() * exit_y.re;
        near(beam.power().transmitted_w, interface_w, 2e-12 * interface_w);
        let dx = width / 4.0;
        let expected: Vec<_> = (0..4).map(|i| {
            interface_w * (2.0 * kz.im * (start + i as f64 * dx)).exp()
                * -(2.0 * kz.im * dx).exp_m1()
        }).collect();
        for (i, cell) in heat.cells().iter().enumerate() {
            assert_eq!(cell.cell_id, bound.cell_ids()[i]);
            near(cell.depth_m[0], start + i as f64 * dx, 2e-16);
            near(cell.depth_m[1], start + (i + 1) as f64 * dx, 2e-16);
            near(cell.absorbed_w, expected[i], 3e-12 * expected[i]);
        }
        near(heat.power().entering_w, interface_w * (2.0 * kz.im * start).exp(),
            3e-12 * interface_w);
        near(heat.power().leaving_w, interface_w * (2.0 * kz.im * (start + width)).exp(),
            3e-12 * interface_w);
        assert!(heat.power().residual_w.abs() < 4e-12 * interface_w);
        let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
        let step = heat.advance_to(
            &initial, dt, conduction(), LayeredThermalSolveControls::default(),
            &mut workspace, || false,
        ).unwrap();
        let powers: Vec<_> = expected.iter().zip(authored).map(|(&p, q)| p + q).collect();
        let temperatures = thermal_reference(&powers, width, AREA, dt);
        for (i, cell) in heat.cells().iter().enumerate() {
            near(step.thermal.temperature_k(cell.cell_id).unwrap(), temperatures[i], 3e-10);
            assert_eq!(step.thermal.cells()[&cell.cell_id].boundary,
                initial.cells()[&cell.cell_id].boundary);
            near(step.cell_energy_j[i].1, expected[i] * dt, 4e-12 * expected[i] * dt);
        }
        near(step.heat_report.energy.external_heat_j,
            step.optical.absorbed_j + authored.iter().sum::<f64>() * dt, 1e-8);
        near(step.heat_report.energy.stored_energy_change_j,
            step.heat_report.energy.external_heat_j, 2e-8);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.cells(), &old_cells);
    }
}

/// G3: native depth refinement preserves integrated heat without introducing
/// artificial optical films. The original film solver's 64-layer cap does not
/// restrict an analytic homogeneous-substrate thermal mesh.
#[test]
fn substrate_refinement_and_partial_coherence_preserve_positive_heat() {
    let exterior = state(&card("refinement incidence", 1.0, 0.0), false, 100_000.0);
    let card = card("refined substrate", 2.0, 0.5);
    let model = debye(&card, 100_000.0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior, DielectricResponse::Debye(&model), OMEGA, &[],
    ).unwrap();
    let factors = [beam(&stack, grid(2.0), 1.0, 0), beam(&stack, grid(2.0), -1.0, 0)];
    let ensemble = AngularSpectrumCoherence::try_new(
        &factors[0], &factors, source_limits(), || false,
    ).unwrap();
    let width = 0.6 / k();
    let mut totals = Vec::new();
    for count in [3, 96] {
        let native = body(&[&card], &[width], &[count], AREA);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let coherent = PreparedSubstrateWaveHeating::try_new(
            &bound, &factors[0], 0.0, heat_limits(), || false,
        ).unwrap();
        let mixed = PreparedSubstrateWaveHeating::try_new(
            &bound, &ensemble, 0.0, heat_limits(), || false,
        ).unwrap();
        assert_eq!(mixed.cells().len(), count);
        assert_eq!(mixed.preparation_mode_cell_pairs(), (count + 2) * 4);
        for (cell, original) in mixed.cells().iter().zip(coherent.cells()) {
            assert!(cell.absorbed_w > 0.0);
            near(cell.absorbed_w, 2.0 * original.absorbed_w, 4e-15 * cell.absorbed_w);
        }
        // The two factors have opposite fields but are explicitly uncorrelated.
        // Adding their fields first would erase this strictly positive heat.
        totals.push(mixed.power().absorbed_w);
        let empty = [];
        let darkness = AngularSpectrumCoherence::try_new(
            &factors[0], &empty, source_limits(), || false,
        ).unwrap();
        let dark = PreparedSubstrateWaveHeating::try_new(
            &bound, &darkness, 0.0, heat_limits(), || false,
        ).unwrap();
        assert_eq!(dark.power().absorbed_w, 0.0);
        assert!(dark.cells().iter().all(|cell| cell.absorbed_w == 0.0));
    }
    near(totals[0], totals[1], 4e-13 * totals[0]);
}

/// G1: transmitted flux through a lossless column is not deposited as heat;
/// total internal reflection also creates no heat from evanescent decay.
#[test]
fn lossless_propagation_and_evanescent_transmission_do_not_heat_substrate_cells() {
    let card = card("lossless substrate", 1.0, 0.0);
    let exit = state(&card, false, 100_000.0);
    let native = body(&[&card], &[0.2], &[3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    for (ratio, spacing) in [(1.0, 2.0), (4.0, 0.7)] {
        let incident = state(&super::card("lossless incidence", ratio, 0.0), false, 100_000.0);
        let stack = PlanarDielectricLayers::try_new(
            &incident, DielectricResponse::Lossless(&exit), OMEGA, &[],
        ).unwrap();
        let beam = beam(&stack, grid(spacing), 1.0, 1);
        let heat = PreparedSubstrateWaveHeating::try_new(
            &bound, &beam, 0.4, heat_limits(), || false,
        ).unwrap();
        assert_eq!(heat.power().absorbed_w, 0.0);
        assert_eq!(heat.power().entering_w, beam.power().transmitted_w);
        assert_eq!(heat.power().leaving_w, beam.power().transmitted_w);
        assert!(heat.cells().iter().all(|cell| cell.absorbed_w == 0.0));
        if ratio == 1.0 {
            assert!(heat.power().entering_w > 0.0);
        }
    }
}

/// G0/G4: no geometry/source substitution, unbounded preparation or accepted
/// thermal mutation is permitted on admission or cancellation failure.
#[test]
fn substrate_source_geometry_budgets_and_cancellation_are_atomic() {
    let exterior = state(&card("admission incidence", 1.0, 0.0), false, 100_000.0);
    let card = card("admitted substrate", 2.0, 0.5);
    let model = debye(&card, 100_000.0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior, DielectricResponse::Debye(&model), OMEGA, &[],
    ).unwrap();
    let beam = beam(&stack, grid(2.0), 1.0, 0);
    let native = body(&[&card], &[0.2], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let heat = PreparedSubstrateWaveHeating::try_new(
        &bound, &beam, 0.0, heat_limits(), || false,
    ).unwrap();
    let mut exact = heat_limits();
    exact.maximum_cells = 4;
    exact.maximum_mode_cell_pairs = heat.preparation_mode_cell_pairs();
    exact.maximum_workspace_bytes = heat.workspace_bytes();
    assert!(PreparedSubstrateWaveHeating::try_new(&bound, &beam, 0.0, exact, || false).is_ok());
    for which in 0..3 {
        let mut short = exact;
        match which {
            0 => short.maximum_cells -= 1,
            1 => short.maximum_mode_cell_pairs -= 1,
            _ => short.maximum_workspace_bytes -= 1,
        }
        assert!(PreparedSubstrateWaveHeating::try_new(&bound, &beam, 0.0, short, || false).is_err());
    }
    for start in [-1.0, f64::NAN, f64::INFINITY, 1e300] {
        assert!(PreparedSubstrateWaveHeating::try_new(
            &bound, &beam, start, heat_limits(), || false,
        ).is_err());
    }
    let foreign = super::card("different substrate card", 2.0, 0.5);
    for (thermal_card, area) in [(&foreign, AREA), (&card, AREA * 2.0)] {
        let other = body(&[thermal_card], &[0.2], &[4], area);
        let other_laws = other.thermal_laws(false).unwrap();
        let other_bound = other.bind(&other_laws, ids()).unwrap();
        assert!(PreparedSubstrateWaveHeating::try_new(
            &other_bound, &beam, 0.0, heat_limits(), || false,
        ).is_err());
    }
    let wrong_pressure = debye(&card, 120_000.0);
    let wrong_stack = PlanarDielectricLayers::try_new(
        &exterior, DielectricResponse::Debye(&wrong_pressure), OMEGA, &[],
    ).unwrap();
    let wrong_beam = self::beam(&wrong_stack, grid(2.0), 1.0, 0);
    assert!(PreparedSubstrateWaveHeating::try_new(
        &bound, &wrong_beam, 0.0, heat_limits(), || false,
    ).is_err());
    assert!(matches!(PreparedSubstrateWaveHeating::try_new(
        &bound, &beam, 0.0, heat_limits(), || true,
    ), Err(LayeredWaveHeatingError::Cancelled)));
    let initial = initial(&bound, &[0.0; 4]);
    let old = initial.cells().clone();
    let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
    for stop in [true, false] {
        let mut calls = 0;
        let result = heat.advance_to(
            &initial, 0.1, conduction(), LayeredThermalSolveControls::default(),
            &mut workspace, || { calls += 1; stop || calls > 12 },
        );
        assert!(matches!(result, Err(LayeredWaveHeatingError::Cancelled)));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.cells(), &old);
    }
}

/// G1: positive local absorption survives when endpoint fluxes round equal.
/// Large fields are an authored numerical fixture, not nonlinear material data.
#[test]
fn weak_substrate_loss_reaches_thermal_storage_below_flux_subtraction_roundoff() {
    let exterior = state(&card("weak incidence", 1.0, 0.0), false, 100_000.0);
    let loss = 1e-20;
    let amplitude = 1e12;
    let card = card("weak substrate", 1.0, loss);
    let model = debye(&card, 100_000.0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior, DielectricResponse::Debye(&model), OMEGA, &[],
    ).unwrap();
    let beam = PlanarAngularSpectrum::prepare(
        &stack,
        AngularSpectrumGrid { width: 1, height: 1, spacing_m: [AREA, 1.0] },
        &[[C64::from_re(amplitude), C64::ZERO]],
        field_limits(),
        || false,
    ).unwrap();
    let width = 0.2;
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let heat = PreparedSubstrateWaveHeating::try_new(
        &bound, &beam, 0.0, heat_limits(), || false,
    ).unwrap();
    let expected = 0.5 * OMEGA * E * loss * amplitude * amplitude * AREA * width;
    assert_eq!(heat.power().entering_w, heat.power().leaving_w);
    assert!(heat.cells().iter().all(|cell| cell.absorbed_w > 0.0));
    near(heat.power().absorbed_w, expected, 2e-12 * expected);
    let initial = initial(&bound, &[0.0; 4]);
    let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
    let step = heat.advance_to(
        &initial, 1.0, conduction(), LayeredThermalSolveControls::default(),
        &mut workspace, || false,
    ).unwrap();
    near(step.optical.absorbed_j, expected, 2e-12 * expected);
    near(step.heat_report.energy.external_heat_j, expected, 1e-8);
    near(step.heat_report.energy.stored_energy_change_j, expected, 3e-8);
    for cell in heat.cells() {
        assert!(step.thermal.temperature_k(cell.cell_id).unwrap() > 300.0);
    }
}
