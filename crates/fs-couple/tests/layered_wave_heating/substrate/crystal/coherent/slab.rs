//! Original finite-film standing waves and angular coherence heat native cells.
use super::*;
use fs_couple::layered_wave_heating::substrate::{
    SubstrateCellPower,
    crystal::{
        CoherentCrystalWavePower, CoherentCrystalWaveStep,
        PreparedCoherentCrystalSlabHeating, PreparedPartiallyCoherentCrystalSlabHeating,
    },
};
use fs_em::plane_wave::uniaxial::oriented::{
    absorbing::slab::{
        AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination,
        PreparedAbsorbingSlabField,
    },
    slab::{OrientedUniaxialSlab, SlabRegion},
};

const FILM_DEPTH: f64 = 2.0 * WAVELENGTH_M;
const BODY_DEPTH: f64 = 1.5 * WAVELENGTH_M;
const BODY_START: f64 = 0.25 * WAVELENGTH_M;
const PERIOD: f64 = 0.8;
const DT: f64 = 2e-9;

fn film_beams<'f, 's, 'm>(
    film: &'f OrientedUniaxialSlab<'s, 'm>,
    amplitudes: [C64; 2],
) -> [PreparedAbsorbingSlabField<'f, 's, 'm>; 2] {
    let p = WAVELENGTH_M / (2.0 * PERIOD);
    core::array::from_fn(|i| {
        film.prepare_absorbing_field_at_tangential_index_in_plane(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            [p, -p][i],
            [1.0, 0.0],
            [C64::ZERO, amplitudes[i]],
        ).unwrap()
    })
}

// Independent Simpson integration of the actual film point field. It retains
// the backward eigenwaves but does not call the production volume-box kernel.
fn film_cell_watts(field: &PreparedAbsorbingSlabField<'_, '_, '_>, a: f64, b: f64) -> f64 {
    let n = 1024;
    let dz = (b - a) / n as f64;
    let mut sum = 0.0;
    for j in 0..=n {
        let z = if j == n { b } else { a + dz * j as f64 };
        let weight = if j == 0 || j == n { 1.0 } else if j % 2 == 0 { 2.0 } else { 4.0 };
        sum += weight * field.field_at(SlabRegion::Crystal, [0.0, 0.0, z])
            .unwrap().absorbed_power_density_w_m3;
    }
    AREA * dz * sum / 3.0
}

fn check_native_step(
    step: &CoherentCrystalWaveStep<'_>,
    cells: &[SubstrateCellPower],
    power: CoherentCrystalWavePower,
    total_cell_watts: &[f64],
    authored: &[f64],
) {
    let reference = thermal_reference(total_cell_watts, BODY_DEPTH, AREA, DT);
    for (i, cell) in cells.iter().enumerate() {
        near(step.thermal.temperature_k(cell.cell_id).unwrap(), reference[i], 5e-8);
        near(step.cell_energy_j[i].1, cell.absorbed_w * DT, 1e-13);
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * DT,
        1e-9,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-9,
    );
    for face in 0..6 {
        near(step.optical.outward_face_j[face], power.outward_face_w[face] * DT, 1e-13);
    }
    assert!(step.optical.residual_j.abs() < 3e-9 * step.optical.absorbed_j);
}

#[test]
fn finite_film_angular_fringes_and_internal_reflection_drive_temperature_once() {
    let outside = state(&card("angular film incidence", 1.0, 0.0), false, 100_000.0);
    let exit = state(&card("angular film exit", 4.0, 0.0), false, 100_000.0);
    let material = crystal_card("angular absorbing film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &interface, DielectricResponse::Lossless(&exit), FILM_DEPTH,
    ).unwrap();
    let components = film_beams(&film, [C64::from_re(1e5); 2]);
    assert!(components[0].response().backward.iter()
        .any(|wave| wave.electric_v_m.iter().any(|v| v.abs() > 1.0)));
    let illumination = AbsorbingSlabCoherentIllumination::try_new(
        &film, &components, 2, || false,
    ).unwrap();
    let native = body(&[&material], &[BODY_DEPTH], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let extent = [0.2, AREA / 0.2];
    let kx = core::f64::consts::PI / PERIOD;
    let mut totals = Vec::new();
    for center in [0.0, 0.5 * PERIOD] {
        let heat = PreparedCoherentCrystalSlabHeating::try_new(
            &bound, &illumination,
            [center - 0.5 * extent[0], -0.5 * extent[1]], extent,
            BODY_START, limits(), || false,
        ).unwrap();
        assert!(core::ptr::eq(heat.source(), &illumination));
        assert!(core::ptr::eq(heat.thermal_body(), &bound));
        assert_eq!(heat.preparation_mode_cell_pairs(), 64 * 4);
        // Equal-|kz| TE waves have the SAME finite-film standing wave.
        // Only the independently integrated transverse fringe factor changes.
        let contrast = 2.0 + 2.0 * (2.0 * kx * center).cos()
            * (kx * extent[0]).sin() / (kx * extent[0]);
        let mut powers = Vec::new();
        for (i, cell) in heat.cells().iter().enumerate() {
            assert_eq!(cell.cell_id, bound.cell_ids()[i]);
            let single = film_cell_watts(&components[0], cell.depth_m[0], cell.depth_m[1]);
            near(cell.absorbed_w, contrast * single, 2e-9 * single);
            powers.push(contrast * single + authored[i]);
        }
        let step = heat.advance_to(
            &initial, DT, conduction(), LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
        ).unwrap();
        check_native_step(&step, heat.cells(), heat.power(), &powers, &authored);
        totals.push(heat.power().absorbed_w);
    }
    assert!(totals[0] > 10.0 * totals[1]);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn independent_film_groups_preserve_standing_waves_and_covariance_visibility() {
    let outside = state(&card("partial film incidence", 1.0, 0.0), false, 100_000.0);
    let exit = state(&card("partial film exit", 4.0, 0.0), false, 100_000.0);
    let material = crystal_card("partially coherent film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &interface, DielectricResponse::Lossless(&exit), FILM_DEPTH,
    ).unwrap();
    let reference = film_beams(&film, [C64::from_re(1e5); 2]);
    let native = body(&[&material], &[BODY_DEPTH], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let extent = [0.2, AREA / 0.2];
    let kx = core::f64::consts::PI / PERIOD;
    for gamma in [-0.65_f64, 0.0, 0.65] {
        // Factorize [[1,gamma],[gamma,1]] without changing per-beam power.
        // Independent groups retain their OWN angular and standing-wave terms.
        let a = 1e5 * ((1.0 + gamma) * 0.5).sqrt();
        let b = 1e5 * ((1.0 - gamma) * 0.5).sqrt();
        let plus = film_beams(&film, [C64::from_re(a); 2]);
        let minus = film_beams(&film, [C64::from_re(b), C64::from_re(-b)]);
        let groups = [
            AbsorbingSlabCoherentIllumination::try_new(&film, &plus, 2, || false).unwrap(),
            AbsorbingSlabCoherentIllumination::try_new(&film, &minus, 2, || false).unwrap(),
        ];
        let illumination = AbsorbingSlabIllumination::try_new(
            &film, &groups, 2, 4, || false,
        ).unwrap();
        let mut totals = Vec::new();
        for center in [0.0, 0.5 * PERIOD] {
            let heat = PreparedPartiallyCoherentCrystalSlabHeating::try_new(
                &bound, &illumination,
                [center - 0.5 * extent[0], -0.5 * extent[1]], extent,
                BODY_START, limits(), || false,
            ).unwrap();
            assert!(core::ptr::eq(heat.source(), &illumination));
            assert!(core::ptr::eq(heat.thermal_body(), &bound));
            assert_eq!(heat.preparation_mode_cell_pairs(), 128 * 4);
            let contrast = 2.0 + 2.0 * gamma * (2.0 * kx * center).cos()
                * (kx * extent[0]).sin() / (kx * extent[0]);
            let mut powers = Vec::new();
            for (i, cell) in heat.cells().iter().enumerate() {
                let single = film_cell_watts(&reference[0], cell.depth_m[0], cell.depth_m[1]);
                near(cell.absorbed_w, contrast * single, 2e-9 * single);
                powers.push(contrast * single + authored[i]);
            }
            let step = heat.advance_to(
                &initial, DT, conduction(), LayeredThermalSolveControls::default(),
                &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
            ).unwrap();
            check_native_step(&step, heat.cells(), heat.power(), &powers, &authored);
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
fn finite_film_boxes_preserve_extent_material_budgets_and_atomic_refusal() {
    let outside = state(&card("bounded angular film incidence", 1.0, 0.0), false, 100_000.0);
    let material = crystal_card("bounded angular film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    ).unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &interface, DielectricResponse::Lossless(&outside), FILM_DEPTH,
    ).unwrap();
    let fields = film_beams(&film, [C64::ONE, C64::new(0.0, 1.0)]);
    let groups = [
        AbsorbingSlabCoherentIllumination::try_new(&film, &fields, 2, || false).unwrap(),
        AbsorbingSlabCoherentIllumination::try_new(&film, &[], 0, || false).unwrap(),
    ];
    let illumination = AbsorbingSlabIllumination::try_new(
        &film, &groups, 2, 2, || false,
    ).unwrap();
    let native = body(&[&material], &[WAVELENGTH_M], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let mut preparation_polls = 0;
    // The last physical cell ends at the original optical exit face.
    let heat = PreparedPartiallyCoherentCrystalSlabHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], WAVELENGTH_M, limits(), || {
            preparation_polls += 1;
            false
        },
    ).unwrap();
    assert_eq!(heat.cells().last().unwrap().depth_m[1], FILM_DEPTH);
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 4,
        maximum_mode_cell_pairs: heat.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: heat.workspace_bytes(),
    };
    assert!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], WAVELENGTH_M, exact, || false,
    ).is_ok());
    for which in 0..3 {
        let mut short = exact;
        match which {
            0 => short.maximum_cells -= 1,
            1 => short.maximum_mode_cell_pairs -= 1,
            _ => short.maximum_workspace_bytes -= 1,
        }
        assert!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], WAVELENGTH_M, short, || false,
        ).is_err());
    }
    for depth in [-WAVELENGTH_M, 1.25 * WAVELENGTH_M, f64::NAN] {
        assert!(PreparedCoherentCrystalSlabHeating::try_new(
            &bound, &groups[0], [0.0; 2], [1.0, AREA], depth, limits(), || false,
        ).is_err());
        assert!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], depth, limits(), || false,
        ).is_err());
    }
    assert!(PreparedCoherentCrystalSlabHeating::try_new(
        &bound, &groups[0], [0.0; 2], [1.0, AREA * 0.5], 0.0, limits(), || false,
    ).is_err());
    let wrong_card = crystal_card("other angular film", [0.12, 0.42]);
    let wrong_native = body(&[&wrong_card], &[WAVELENGTH_M], &[4], AREA);
    let wrong_laws = wrong_native.thermal_laws(false).unwrap();
    let wrong = wrong_native.bind(&wrong_laws, ids()).unwrap();
    assert!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
        &wrong, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    ).is_err());
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut polls = 0;
        assert!(matches!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], WAVELENGTH_M, limits(), || {
                polls += 1;
                polls == stop
            },
        ), Err(LayeredWaveHeatingError::Cancelled)));
    }
    let initial = initial(&bound, &[1.0, 2.0, 3.0, 4.0]);
    let original = initial.cells().clone();
    assert!(matches!(heat.advance_to(
        &initial, DT, conduction(), LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(), || true,
    ), Err(LayeredWaveHeatingError::Cancelled)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);

    // Darkness does not create permission to put thermal cells outside the film.
    let empty = AbsorbingSlabIllumination::try_new(&film, &[], 0, 0, || false).unwrap();
    assert!(PreparedPartiallyCoherentCrystalSlabHeating::try_new(
        &bound, &empty, [0.0; 2], [1.0, AREA], 1.25 * WAVELENGTH_M, limits(), || false,
    ).is_err());
    assert!(PreparedCoherentCrystalSlabHeating::try_new(
        &bound, &groups[1], [0.0; 2], [1.0, AREA], 1.25 * WAVELENGTH_M, limits(), || false,
    ).is_err());
    let dark = PreparedPartiallyCoherentCrystalSlabHeating::try_new(
        &bound, &empty, [0.0; 2], [1.0, AREA], WAVELENGTH_M, limits(), || false,
    ).unwrap();
    assert_eq!(dark.preparation_mode_cell_pairs(), 0);
    assert_eq!(dark.power().absorbed_w, 0.0);
    assert_eq!(dark.power().outward_face_w, [0.0; 6]);
    let step = dark.advance_to(
        &initial, DT, conduction(), LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(), || false,
    ).unwrap();
    assert_eq!(step.optical.absorbed_j, 0.0);
    near(step.heat_report.energy.external_heat_j, 10.0 * DT, 1e-14);
}
