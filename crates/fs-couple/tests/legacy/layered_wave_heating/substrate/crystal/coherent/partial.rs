//! Partial spatial coherence and spectral independence reach native cells.
use super::*;
use fs_couple::layered_wave_heating::substrate::crystal::
    PreparedPartiallyCoherentCrystalHeating;
use fs_em::plane_wave::uniaxial::oriented::absorbing::
    AbsorbingCrystalPartiallyCoherentIllumination;

#[test]
fn partial_coherence_visibility_drives_native_cell_temperatures() {
    let outside = state(&card("partial footprint incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("partial footprint crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let period = 0.8;
    let p = WAVELENGTH_M / (2.0 * period);
    let reference_beams = beams(&interface, [p, -p], [C64::from_re(1e5); 2]);
    let extent = [0.2, AREA / 0.2];
    let width = 2.0 * WAVELENGTH_M;
    let start = 0.1 * WAVELENGTH_M;
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let dt = 2e-9;
    let kx = core::f64::consts::PI / period;
    for gamma in [-1.0_f64, -0.6, 0.0, 0.65, 1.0] {
        // Factorize [[1, gamma], [gamma, 1]] at fixed per-beam power.
        // Positive/negative gamma moves the bright/dark footprints; zero
        // gamma removes their cross-beam fringe without removing mode loss.
        let a = 1e5 * ((1.0 + gamma) * 0.5).sqrt();
        let b = 1e5 * ((1.0 - gamma) * 0.5).sqrt();
        // The native field owner requires a nonzero incident amplitude.
        // A zero statistical factor is represented by its empty field bank.
        let plus = if a > 0.0 {
            beams(&interface, [p, -p], [C64::from_re(a); 2])
                .into_iter().collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let minus = if b > 0.0 {
            beams(&interface, [p, -p], [C64::from_re(b), C64::from_re(-b)])
                .into_iter().collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let factors = [
            AbsorbingCrystalCoherentIllumination::try_new(&interface, &plus, 2, || false)
                .unwrap(),
            AbsorbingCrystalCoherentIllumination::try_new(&interface, &minus, 2, || false)
                .unwrap(),
        ];
        let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
            &interface, &factors, 2, || false,
        ).unwrap();
        let mut totals = Vec::new();
        for center in [0.0, 0.5 * period] {
            let heat = PreparedPartiallyCoherentCrystalHeating::try_new(
                &bound, &illumination,
                [center - 0.5 * extent[0], -0.5 * extent[1]], extent,
                start, limits(), || false,
            ).unwrap();
            assert!(core::ptr::eq(heat.source(), &illumination));
            assert!(core::ptr::eq(heat.thermal_body(), &bound));
            assert_eq!(heat.preparation_mode_cell_pairs(), illumination.mode_pair_count() * 4);
            // Independent closed-form transverse integral of the stated
            // covariance, multiplied by the original one-beam depth loss.
            let contrast = 2.0 + 2.0 * gamma * (2.0 * kx * center).cos()
                * (kx * extent[0]).sin() / (kx * extent[0]);
            let mut powers = Vec::new();
            for (i, cell) in heat.cells().iter().enumerate() {
                let single = reference_beams[0]
                    .slab_absorption(AREA, cell.depth_m[0], cell.depth_m[1])
                    .unwrap().absorbed_power_w();
                near(cell.absorbed_w, contrast * single, 3e-10 * single);
                powers.push(contrast * single + authored[i]);
            }
            let step = heat.advance_to(
                &initial, dt, conduction(), LayeredThermalSolveControls::default(),
                &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
            ).unwrap();
            let reference = thermal_reference(&powers, width, AREA, dt);
            for (i, cell) in heat.cells().iter().enumerate() {
                near(step.thermal.temperature_k(cell.cell_id).unwrap(), reference[i], 5e-8);
                near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-13);
            }
            near(step.heat_report.energy.external_heat_j,
                step.optical.absorbed_j + authored.iter().sum::<f64>() * dt, 1e-9);
            near(step.heat_report.energy.stored_energy_change_j,
                step.heat_report.energy.external_heat_j, 2e-9);
            for face in 0..6 {
                near(step.optical.outward_face_j[face],
                    heat.power().outward_face_w[face] * dt, 1e-13);
            }
            assert!(step.optical.residual_j.abs() < 3e-9 * step.optical.absorbed_j);
            totals.push(heat.power().absorbed_w);
        }
        if gamma == 0.0 {
            near(totals[0], totals[1], 3e-11 * totals[0]);
        } else {
            assert_eq!(totals[0] > totals[1], gamma > 0.0);
        }
    }
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn independent_spectral_angular_factors_keep_cell_heat_and_six_face_work() {
    let outside = state(&card("spectral footprint incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("spectral footprint crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let fields: [[PreparedAbsorbingCrystalField<'_, '_>; 2]; 2] =
        core::array::from_fn(|factor| {
            let wavelength = [500e-9, 650e-9][factor];
            let p = wavelength / (2.0 * [0.8, 0.5][factor]);
            let amplitudes = [
                [C64::new(0.8e5, 0.2e5), C64::new(-0.3e5, 0.6e5)],
                [C64::new(0.4e5, -0.7e5), C64::new(0.5e5, 0.1e5)],
            ][factor];
            core::array::from_fn(|beam| {
                interface.prepare_absorbing_field_at_tangential_index_in_plane(
                    core::f64::consts::TAU * C0 / wavelength,
                    if beam == 0 { p } else { -p }, [1.0, 0.0],
                    [C64::ZERO, amplitudes[beam]],
                ).unwrap()
            })
        });
    let factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[0], 2, || false)
            .unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[1], 2, || false)
            .unwrap(),
    ];
    let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors, 2, || false,
    ).unwrap();
    let width = 2.0 * WAVELENGTH_M;
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let origin = [-0.1, -1.5];
    let extent = [0.2, 3.0];
    let first = PreparedCoherentCrystalHeating::try_new(
        &bound, &factors[0], origin, extent, 0.0, limits(), || false,
    ).unwrap();
    let second = PreparedCoherentCrystalHeating::try_new(
        &bound, &factors[1], origin, extent, 0.0, limits(), || false,
    ).unwrap();
    let heat = PreparedPartiallyCoherentCrystalHeating::try_new(
        &bound, &illumination, origin, extent, 0.0, limits(), || false,
    ).unwrap();
    for face in 0..6 {
        let expected = first.power().outward_face_w[face] + second.power().outward_face_w[face];
        near(heat.power().outward_face_w[face], expected,
            3e-13 * heat.power().absorbed_w);
    }
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let mut powers = Vec::new();
    for (i, cell) in heat.cells().iter().enumerate() {
        let expected = first.cells()[i].absorbed_w + second.cells()[i].absorbed_w;
        near(cell.absorbed_w, expected, 3e-13 * expected);
        powers.push(expected + authored[i]);
    }
    let dt = 2e-9;
    let step = heat.advance_to(
        &initial, dt, conduction(), LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
    ).unwrap();
    for (i, expected) in thermal_reference(&powers, width, AREA, dt).iter().enumerate() {
        near(step.thermal.temperature_k(bound.cell_ids()[i]).unwrap(), *expected, 5e-8);
    }
    near(step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * dt, 1e-9);
    assert!(step.optical.residual_j.abs() < 3e-9 * step.optical.absorbed_j);
}

#[test]
fn partial_crystal_heat_bounds_work_preserves_sources_and_cancels_atomically() {
    let outside = state(&card("bounded partial incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("bounded partial crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let fields = beams(&interface, [0.25, -0.25], [C64::ONE, C64::new(0.0, 1.0)]);
    let factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields, 2, || false)
            .unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false)
            .unwrap(),
    ];
    let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors, 2, || false,
    ).unwrap();
    let native = body(&[&card], &[WAVELENGTH_M], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let mut preparation_polls = 0;
    let heat = PreparedPartiallyCoherentCrystalHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || {
            preparation_polls += 1;
            false
        },
    ).unwrap();
    assert_eq!(heat.preparation_mode_cell_pairs(), 16 * 4);
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 4,
        maximum_mode_cell_pairs: heat.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: heat.workspace_bytes(),
    };
    assert!(PreparedPartiallyCoherentCrystalHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, exact, || false,
    ).is_ok());
    for which in 0..3 {
        let mut short = exact;
        match which {
            0 => short.maximum_cells -= 1,
            1 => short.maximum_mode_cell_pairs -= 1,
            _ => short.maximum_workspace_bytes -= 1,
        }
        assert!(PreparedPartiallyCoherentCrystalHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, short, || false,
        ).is_err());
    }
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut polls = 0;
        assert!(matches!(PreparedPartiallyCoherentCrystalHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || {
                polls += 1;
                polls == stop
            },
        ), Err(LayeredWaveHeatingError::Cancelled)));
    }
    let other_card = crystal_card("other partial crystal", [0.12, 0.42]);
    let other_native = body(&[&other_card], &[WAVELENGTH_M], &[4], AREA);
    let other_laws = other_native.thermal_laws(false).unwrap();
    let other = other_native.bind(&other_laws, ids()).unwrap();
    assert!(PreparedPartiallyCoherentCrystalHeating::try_new(
        &other, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    ).is_err());
    let initial = initial(&bound, &[1.0, 2.0, 3.0, 4.0]);
    let original = initial.cells().clone();
    assert!(matches!(heat.advance_to(
        &initial, 1e-9, conduction(), LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(), || true,
    ), Err(LayeredWaveHeatingError::Cancelled)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);

    let empty = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &[], 0, || false,
    ).unwrap();
    let dark = PreparedPartiallyCoherentCrystalHeating::try_new(
        &bound, &empty, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    ).unwrap();
    assert_eq!(dark.preparation_mode_cell_pairs(), 0);
    assert_eq!(dark.power().absorbed_w, 0.0);
    assert_eq!(dark.power().outward_face_w, [0.0; 6]);
    let dt = 1e-6;
    let step = dark.advance_to(
        &initial, dt, conduction(), LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
    ).unwrap();
    assert_eq!(step.optical.absorbed_j, 0.0);
    for (i, expected) in thermal_reference(
        &[1.0, 2.0, 3.0, 4.0], WAVELENGTH_M, AREA, dt,
    ).iter().enumerate() {
        near(step.thermal.temperature_k(bound.cell_ids()[i]).unwrap(), *expected, 5e-8);
    }
    near(step.heat_report.energy.external_heat_j, 10.0 * dt, 1e-12);
}
