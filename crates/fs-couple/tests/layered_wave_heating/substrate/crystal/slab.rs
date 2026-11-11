//! Finite-crystal standing waves feed the same physical thermal cells.
use super::*;
use fs_couple::layered_wave_heating::substrate::crystal::PreparedCrystalSlabHeating;
use fs_em::plane_wave::uniaxial::oriented::{
    absorbing::slab::PreparedAbsorbingSlabField,
    slab::{OrientedUniaxialSlab, SlabRegion},
};

#[path = "slab/illumination.rs"]
mod illumination;

fn film_density_integral(field: &PreparedAbsorbingSlabField<'_, '_, '_>, a: f64, b: f64) -> f64 {
    let n = 1024;
    let dz = (b - a) / n as f64;
    let mut sum = 0.0;
    for j in 0..=n {
        let z = if j == n { b } else { a + dz * j as f64 };
        let weight = if j == 0 || j == n {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight
            * field
                .field_at(SlabRegion::Crystal, [0.0, 0.0, z])
                .unwrap()
                .absorbed_power_density_w_m3;
    }
    AREA * dz * sum / 3.0
}

#[test]
fn finite_crystal_cell_integrals_drive_native_temperature_once() {
    let incident = state(&card("film incidence", 1.0, 0.0), false, 100_000.0);
    let exit = state(&card("film exit", 4.0, 0.0), false, 100_000.0);
    let material = crystal_card("heated crystal film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&exit),
        2.0 * WAVELENGTH_M,
    )
    .unwrap();
    let field = film
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::new(1.4e5, 0.4e5), C64::new(-0.6e5, 1.0e5)],
        )
        .unwrap();
    assert!(
        field
            .response()
            .backward
            .iter()
            .any(|wave| wave.electric_v_m.iter().any(|v| v.abs() > 1.0))
    );
    let width = 1.5 * WAVELENGTH_M;
    let start = 0.25 * WAVELENGTH_M;
    let native = body(&[&material], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let heat =
        PreparedCrystalSlabHeating::try_new(&bound, &field, start, limits(), || false).unwrap();
    assert!(std::ptr::eq(heat.source(), &field));
    assert!(std::ptr::eq(heat.thermal_body(), &bound));
    assert_eq!(heat.preparation_mode_cell_pairs(), 16 * (4 + 2));
    let mut powers = Vec::new();
    for (i, cell) in heat.cells().iter().enumerate() {
        assert_eq!(cell.cell_id, bound.cell_ids()[i]);
        let expected = film_density_integral(&field, cell.depth_m[0], cell.depth_m[1]);
        assert!(expected > 0.0);
        near(cell.absorbed_w, expected, 2e-9 * expected);
        powers.push(expected + authored[i]);
    }
    assert!(heat.power().residual_w.abs() < 3e-11 * heat.power().entering_w);
    let dt = 2e-9;
    let step = heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(4).unwrap(),
            || false,
        )
        .unwrap();
    let temperatures = thermal_reference(&powers, width, AREA, dt);
    for (i, cell) in heat.cells().iter().enumerate() {
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            temperatures[i],
            3e-8,
        );
        near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-14);
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * dt,
        1e-9,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-9,
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn coherent_two_face_crystal_heat_and_reverse_entry_drive_native_energy_once() {
    let exterior = state(&card("two-face exterior", 1.0, 0.0), false, 100_000.0);
    let material = crystal_card("two-face heated film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&exterior),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 0.0, 0.0],
    )
    .unwrap();
    let width = 0.125 * WAVELENGTH_M;
    let film =
        OrientedUniaxialSlab::try_new(&entrance, DielectricResponse::Lossless(&exterior), width)
            .unwrap();
    let native = body(&[&material], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[0.0; 4]);
    let original = initial.cells().clone();
    let drive = [C64::ZERO, C64::new(5e4, 0.0)];
    let opposite = [C64::ZERO, C64::new(-5e4, 0.0)];
    let dark = [C64::ZERO; 2];
    let mut total_heat = Vec::new();
    for (left, right, both_faces) in [
        (drive, drive, true),
        (drive, opposite, true),
        (dark, drive, false),
    ] {
        let field = film
            .prepare_absorbing_field_from_both_sides(
                core::f64::consts::TAU * C0 / WAVELENGTH_M,
                1.0,
                left,
                right,
            )
            .unwrap();
        let heat =
            PreparedCrystalSlabHeating::try_new(&bound, &field, 0.0, limits(), || false).unwrap();
        assert_eq!(heat.preparation_mode_cell_pairs(), 16 * (4 + 2));
        let top_w = AREA
            * field
                .field_at(SlabRegion::Crystal, [0.0, 0.0, 0.0])
                .unwrap()
                .field
                .poynting_w_m2[2];
        let bottom_w = AREA
            * field
                .field_at(SlabRegion::Crystal, [0.0, 0.0, width])
                .unwrap()
                .field
                .poynting_w_m2[2];
        assert!(bottom_w < 0.0);
        if both_faces {
            assert!(top_w > 0.0);
            near(heat.power().entering_w, top_w - bottom_w, 1e-9);
            assert_eq!(heat.power().leaving_w, 0.0);
        } else {
            assert!(top_w < 0.0);
            assert_eq!(heat.power().entering_w, -bottom_w);
            assert_eq!(heat.power().leaving_w, -top_w);
        }
        assert!(heat.power().absorbed_w > 0.0);
        assert!(heat.power().residual_w.abs() < 3e-11 * heat.power().entering_w);
        total_heat.push(heat.power().absorbed_w);
        for cell in heat.cells() {
            let expected = film_density_integral(&field, cell.depth_m[0], cell.depth_m[1]);
            near(cell.absorbed_w, expected, 2e-9 * expected);
        }
        let dt = 2e-9;
        let step = heat
            .advance_to(
                &initial,
                dt,
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut LayeredThermalWorkspace::new(4).unwrap(),
                || false,
            )
            .unwrap();
        for (i, cell) in heat.cells().iter().enumerate() {
            near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-14);
        }
        near(step.optical.entering_j, heat.power().entering_w * dt, 1e-14);
        near(step.optical.leaving_j, heat.power().leaving_w * dt, 1e-14);
        near(
            step.heat_report.energy.external_heat_j,
            step.optical.absorbed_j,
            1e-9,
        );
        near(
            step.heat_report.energy.stored_energy_change_j,
            step.optical.absorbed_j,
            2e-9,
        );
    }
    // G1: equal input powers with different relative phase produce different
    // coherent dielectric work, which survives the native thermal coupling.
    assert!((total_heat[0] - total_heat[1]).abs() > 0.1 * (total_heat[0] + total_heat[1]));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn finite_crystal_refinement_preserves_heat_and_transparent_films_do_not_heat() {
    let incident = state(
        &card("film refinement incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let exit = state(&card("film refinement exit", 4.0, 0.0), false, 100_000.0);
    for loss in [[0.12, 0.42], [0.0, 0.0]] {
        let material = crystal_card("refined crystal film", loss);
        let [o, e] = axes(&material, 100_000.0);
        let entrance = OrientedUniaxialInterface::try_new(
            DielectricResponse::Lossless(&incident),
            UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
            [1.0, 2.0, -0.6],
        )
        .unwrap();
        let film = OrientedUniaxialSlab::try_new(
            &entrance,
            DielectricResponse::Lossless(&exit),
            2.0 * WAVELENGTH_M,
        )
        .unwrap();
        let field = film
            .prepare_absorbing_field(
                core::f64::consts::TAU * C0 / WAVELENGTH_M,
                0.6,
                [C64::ONE, C64::new(0.3, 0.7)],
            )
            .unwrap();
        let mut totals = Vec::new();
        for count in [1, 4, 16] {
            let native = body(&[&material], &[WAVELENGTH_M], &[count], AREA);
            let laws = native.thermal_laws(false).unwrap();
            let bound = native.bind(&laws, ids()).unwrap();
            let heat = PreparedCrystalSlabHeating::try_new(
                &bound,
                &field,
                0.25 * WAVELENGTH_M,
                limits(),
                || false,
            )
            .unwrap();
            assert_eq!(heat.preparation_mode_cell_pairs(), 16 * (count + 2));
            assert!(heat.cells().iter().all(|cell| cell.absorbed_w >= 0.0));
            totals.push(heat.power().absorbed_w);
            if loss == [0.0, 0.0] {
                assert_eq!(heat.power().absorbed_w, 0.0);
                assert!(heat.power().entering_w > 0.0);
                near(heat.power().leaving_w, heat.power().entering_w, 2e-15);
            }
        }
        for &total in &totals[1..] {
            near(total, totals[0], 2e-12 * totals[0]);
        }
    }
}

#[test]
fn finite_crystal_heat_enforces_optical_extent_material_and_atomic_cancellation() {
    let incident = state(
        &card("film admission incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let material = crystal_card("admitted crystal film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&incident),
        2.0 * WAVELENGTH_M,
    )
    .unwrap();
    let field = film
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::ONE, C64::ZERO],
        )
        .unwrap();
    let native = body(&[&material], &[WAVELENGTH_M], &[3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let heat =
        PreparedCrystalSlabHeating::try_new(&bound, &field, 0.0, limits(), || false).unwrap();
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 3,
        maximum_mode_cell_pairs: heat.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: heat.workspace_bytes(),
    };
    assert!(PreparedCrystalSlabHeating::try_new(&bound, &field, 0.0, exact, || false).is_ok());
    assert!(
        PreparedCrystalSlabHeating::try_new(
            &bound,
            &field,
            0.0,
            SubstrateWaveHeatingLimits {
                maximum_mode_cell_pairs: exact.maximum_mode_cell_pairs - 1,
                ..exact
            },
            || false
        )
        .is_err()
    );
    // The thermal body's right face would extend past the actual film.
    assert!(
        PreparedCrystalSlabHeating::try_new(&bound, &field, 1.25 * WAVELENGTH_M, limits(), || {
            false
        },)
        .is_err()
    );
    let other = crystal_card("different film material", [0.12, 0.42]);
    let wrong_native = body(&[&other], &[WAVELENGTH_M], &[3], AREA);
    let wrong_laws = wrong_native.thermal_laws(false).unwrap();
    let wrong = wrong_native.bind(&wrong_laws, ids()).unwrap();
    assert!(PreparedCrystalSlabHeating::try_new(&wrong, &field, 0.0, limits(), || false).is_err());
    assert!(matches!(
        PreparedCrystalSlabHeating::try_new(&bound, &field, 0.0, limits(), || true,),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    let initial = initial(&bound, &[0.0; 3]);
    let original = initial.cells().clone();
    assert!(matches!(
        heat.advance_to(
            &initial,
            1e-9,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(3).unwrap(),
            || true,
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}
