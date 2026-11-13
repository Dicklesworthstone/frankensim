//! G1/G3/G4: coherent finite-volume heat reaches the original thermal owner.

use super::*;
use fs_couple::layered_wave_heating::substrate::crystal::PreparedCoherentCrystalHeating;
use fs_em::plane_wave::uniaxial::oriented::absorbing::{
    AbsorbingCrystalCoherentIllumination, AbsorbingCrystalBoxLimits,
};

fn beams<'s, 'm>(
    interface: &'s OrientedUniaxialInterface<'m>,
    tangential: [f64; 2],
    amplitudes: [C64; 2],
) -> [PreparedAbsorbingCrystalField<'s, 'm>; 2] {
    core::array::from_fn(|i| {
        interface
            .prepare_absorbing_field_at_tangential_index_in_plane(
                core::f64::consts::TAU * C0 / WAVELENGTH_M,
                tangential[i],
                [1.0, 0.0],
                [C64::ZERO, amplitudes[i]],
            )
            .unwrap()
    })
}

#[test]
fn coherent_crystal_footprint_fringes_drive_native_temperature_once() {
    let outside = state(&card("coherent footprint incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("coherent footprint crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    )
    .unwrap();
    // Nearly parallel beams make macroscopic fringes at an optical carrier.
    // The thermal area is unchanged between a bright and a dark footprint.
    let period = 0.8;
    let p = WAVELENGTH_M / (2.0 * period);
    let components = beams(&interface, [p, -p], [C64::from_re(1e5); 2]);
    let illumination =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || false)
            .unwrap();
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
    let mut totals = Vec::new();
    for center in [0.0, 0.5 * period] {
        let origin = [center - 0.5 * extent[0], -0.5 * extent[1]];
        let heat = PreparedCoherentCrystalHeating::try_new(
            &bound,
            &illumination,
            origin,
            extent,
            start,
            limits(),
            || false,
        )
        .unwrap();
        assert!(core::ptr::eq(heat.source(), &illumination));
        assert!(core::ptr::eq(heat.thermal_body(), &bound));
        assert_eq!(heat.preparation_mode_cell_pairs(), 16 * 4);
        // Integrating |exp(-ikx x)+exp(+ikx x)|² over the footprint
        // gives this analytic contrast, independent of the production box sum.
        let contrast = 2.0
            + 2.0 * (2.0 * kx * center).cos() * (kx * extent[0]).sin()
                / (kx * extent[0]);
        let mut powers = Vec::new();
        for (i, cell) in heat.cells().iter().enumerate() {
            let single = components[0]
                .slab_absorption(AREA, cell.depth_m[0], cell.depth_m[1])
                .unwrap()
                .absorbed_power_w();
            near(cell.absorbed_w, contrast * single, 3e-10 * single);
            powers.push(contrast * single + authored[i]);
        }
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
        let reference = thermal_reference(&powers, width, AREA, dt);
        for (i, cell) in heat.cells().iter().enumerate() {
            near(step.thermal.temperature_k(cell.cell_id).unwrap(), reference[i], 5e-8);
            near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-13);
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
        for face in 0..6 {
            near(
                step.optical.outward_face_j[face],
                heat.power().outward_face_w[face] * dt,
                1e-13,
            );
        }
        assert!(step.optical.residual_j.abs() < 3e-9 * step.optical.absorbed_j);
        totals.push(heat.power().absorbed_w);
    }
    assert!(totals[0] > 10.0 * totals[1]);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

// Independent tensor-product Simpson integration of local physical fields.
// The production adapter uses the native analytic box integrals and never
// deposits this sampled density or substitutes flux difference for heat.
fn point_integrals(
    illumination: &AbsorbingCrystalCoherentIllumination<'_, '_, '_>,
    origin: [f64; 3],
    extent: [f64; 3],
) -> (f64, [f64; 6]) {
    let n = 24;
    let weight = |i| if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
    let dx = extent.map(|v| v / n as f64);
    let side = OrientedCrystalFieldSide::Transmitted;
    let mut density = 0.0;
    for i in 0..=n {
        for j in 0..=n {
            for k in 0..=n {
                let point = [
                    origin[0] + i as f64 * dx[0],
                    origin[1] + j as f64 * dx[1],
                    origin[2] + k as f64 * dx[2],
                ];
                density += weight(i) * weight(j) * weight(k)
                    * illumination.field(side, point, || false).unwrap()
                        .absorbed_power_density_w_m3;
            }
        }
    }
    let heat = density * dx[0] * dx[1] * dx[2] / 27.0;
    let mut faces = [0.0; 6];
    for axis in 0..3 {
        let transverse = [(axis + 1) % 3, (axis + 2) % 3];
        for boundary in 0..2 {
            let mut sum = 0.0;
            for i in 0..=n {
                for j in 0..=n {
                    let mut point = origin;
                    point[axis] += boundary as f64 * extent[axis];
                    point[transverse[0]] += i as f64 * dx[transverse[0]];
                    point[transverse[1]] += j as f64 * dx[transverse[1]];
                    sum += weight(i) * weight(j)
                        * illumination.poynting_w_m2(side, point, || false).unwrap()[axis];
                }
            }
            let sign = if boundary == 0 { -1.0 } else { 1.0 };
            faces[2 * axis + boundary] =
                sign * sum * dx[transverse[0]] * dx[transverse[1]] / 9.0;
        }
    }
    (heat, faces)
}

#[test]
fn coherent_crystal_column_keeps_lateral_flux_and_native_depth_refinement() {
    let outside = state(&card("finite box incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("finite box crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    )
    .unwrap();
    let components = beams(&interface, [0.25, 0.75], [C64::ONE, C64::new(0.4, 0.9)]);
    let illumination =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || false)
            .unwrap();
    let origin = [0.13 * WAVELENGTH_M, -0.17 * WAVELENGTH_M, 0.1 * WAVELENGTH_M];
    let extent = [0.31 * WAVELENGTH_M, 0.4 * WAVELENGTH_M, 0.9 * WAVELENGTH_M];
    let area = extent[0] * extent[1];
    let (heat_reference, flux_reference) = point_integrals(&illumination, origin, extent);
    let scale = flux_reference.iter().map(|v| v.abs()).sum::<f64>();
    let mut totals = Vec::new();
    for count in [1, 4, 16] {
        let native = body(&[&card], &[extent[2]], &[count], area);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let heat = PreparedCoherentCrystalHeating::try_new(
            &bound,
            &illumination,
            [origin[0], origin[1]],
            [extent[0], extent[1]],
            origin[2],
            limits(),
            || false,
        )
        .unwrap();
        assert_eq!(heat.cells().len(), count);
        near(heat.power().absorbed_w, heat_reference, 2e-6 * heat_reference);
        for face in 0..6 {
            near(heat.power().outward_face_w[face], flux_reference[face], 2e-6 * scale);
        }
        assert!(heat.power().residual_w.abs() < 3e-10 * scale);
        let lateral = heat.power().outward_face_w[..4].iter().sum::<f64>();
        assert!(lateral.abs() > 1e-3 * heat.power().absorbed_w);
        // Normal faces alone are an incomplete optical energy account.
        let normal_only = -heat.power().outward_face_w[4]
            - heat.power().outward_face_w[5] - heat.power().absorbed_w;
        near(normal_only, lateral, 3e-10 * scale);
        totals.push(heat.power().absorbed_w);
    }
    for &total in &totals[1..] {
        near(total, totals[0], 3e-11 * totals[0]);
    }
}

#[test]
fn coherent_crystal_cancellation_deposits_only_authored_heat() {
    let outside = state(&card("dark crystal incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("dark coherent crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    )
    .unwrap();
    let components = beams(
        &interface,
        [0.3, 0.3],
        [C64::from_re(1e5), C64::from_re(-1e5)],
    );
    let illumination =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || false)
            .unwrap();
    let width = WAVELENGTH_M;
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let heat = PreparedCoherentCrystalHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    )
    .unwrap();
    assert_eq!(heat.power().absorbed_w, 0.0);
    assert_eq!(heat.power().outward_face_w, [0.0; 6]);
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let dt = 1e-6;
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
    assert_eq!(step.optical.absorbed_j, 0.0);
    for (i, temperature) in thermal_reference(&authored, width, AREA, dt).into_iter().enumerate() {
        near(step.thermal.temperature_k(bound.cell_ids()[i]).unwrap(), temperature, 5e-8);
    }
    near(step.heat_report.energy.external_heat_j, authored.iter().sum::<f64>() * dt, 1e-12);
    let empty = AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false)
        .unwrap();
    let dark = PreparedCoherentCrystalHeating::try_new(
        &bound, &empty, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    )
    .unwrap();
    assert_eq!(dark.preparation_mode_cell_pairs(), 0);
    assert_eq!(dark.power(), heat.power());
}

#[test]
fn coherent_crystal_heat_preserves_geometry_source_budgets_and_atomic_cancellation() {
    let outside = state(&card("bounded crystal incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("bounded coherent crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [0.0, 0.0, 1.0],
    )
    .unwrap();
    let components = beams(&interface, [0.25, -0.25], [C64::ONE; 2]);
    let illumination =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || false)
            .unwrap();
    let native = body(&[&card], &[WAVELENGTH_M], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let mut preparation_polls = 0;
    let prepared = PreparedCoherentCrystalHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || {
            preparation_polls += 1;
            false
        },
    )
    .unwrap();
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 4,
        maximum_mode_cell_pairs: prepared.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: prepared.workspace_bytes(),
    };
    assert!(PreparedCoherentCrystalHeating::try_new(
        &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, exact, || false,
    ).is_ok());
    for which in 0..3 {
        let mut short = exact;
        match which {
            0 => short.maximum_cells -= 1,
            1 => short.maximum_mode_cell_pairs -= 1,
            _ => short.maximum_workspace_bytes -= 1,
        }
        assert!(PreparedCoherentCrystalHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, short, || false,
        ).is_err());
    }
    for extent in [[1.0, AREA * 0.5], [1.0, 0.0], [f64::NAN, AREA]] {
        assert!(PreparedCoherentCrystalHeating::try_new(
            &bound, &illumination, [0.0; 2], extent, 0.0, limits(), || false,
        ).is_err());
    }
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut polls = 0;
        assert!(matches!(PreparedCoherentCrystalHeating::try_new(
            &bound, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || {
                polls += 1;
                polls == stop
            },
        ), Err(LayeredWaveHeatingError::Cancelled)));
    }
    let other_card = crystal_card("different coherent crystal", [0.12, 0.42]);
    let other_native = body(&[&other_card], &[WAVELENGTH_M], &[4], AREA);
    let other_laws = other_native.thermal_laws(false).unwrap();
    let other = other_native.bind(&other_laws, ids()).unwrap();
    assert!(PreparedCoherentCrystalHeating::try_new(
        &other, &illumination, [0.0; 2], [1.0, AREA], 0.0, limits(), || false,
    ).is_err());
    let initial = initial(&bound, &[0.0; 4]);
    let original = initial.cells().clone();
    assert!(matches!(prepared.advance_to(
        &initial,
        1e-9,
        conduction(),
        LayeredThermalSolveControls::default(),
        &mut LayeredThermalWorkspace::new(4).unwrap(),
        || true,
    ), Err(LayeredWaveHeatingError::Cancelled)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
    // The integrated box remains bound to the same illumination object.
    let check = illumination.box_absorption(
        fs_em::plane_wave::uniaxial::oriented::absorbing::AbsorbingCrystalBox::try_new(
            [0.0; 3], [1.0, AREA, WAVELENGTH_M],
        ).unwrap(),
        AbsorbingCrystalBoxLimits {
            maximum_mode_pairs: 16,
            maximum_workspace_bytes: 1 << 20,
        },
        || false,
    ).unwrap();
    assert!(core::ptr::eq(check.illumination(), &illumination));
}

#[path = "coherent/partial.rs"]
mod partial;

#[path = "coherent/slab.rs"]
mod slab;
