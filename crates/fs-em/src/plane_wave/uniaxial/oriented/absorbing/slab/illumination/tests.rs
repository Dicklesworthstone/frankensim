//! Direct numerical regressions through the original solved Maxwell fields.
use super::super::super::super::Incidence;
use super::super::{
    AbsorbingCrystalField, AbsorbingSlabResponse, C0, FieldSum, SlabRegion, crystal, kernel,
};
use fs_math::{c64::C64, det};

fn solved(electric: [C64; 2], index: f64, plane: [f64; 2]) -> AbsorbingSlabResponse {
    let mu = 1e-6;
    let omega = core::f64::consts::TAU * C0 / 1e-6;
    let axis = [1.0, 2.0, 3.0].map(|value| value / det::sqrt(14.0));
    let indices = [C64::new(1.6, -0.08), C64::new(2.0, -0.2)];
    let loss_f_m = indices.map(|n| -2.0 * n.re * n.im / (mu * C0 * C0));
    let mut response = kernel::solve(
        crystal::Input {
            omega,
            ni: 1.0,
            indices,
            mui: mu,
            mut_: mu,
            axis: [
                axis[0] * plane[0] + axis[1] * plane[1],
                axis[1] * plane[0] - axis[0] * plane[1],
                axis[2],
            ],
            incidence: Incidence::TangentialIndex(index),
            electric: [
                electric[0].scale(plane[0]) + electric[1].scale(plane[1]),
                electric[1].scale(plane[0]) - electric[0].scale(plane[1]),
            ],
            loss_f_m,
        },
        1.2,
        0.8 * mu,
        1.3e-6,
    )
    .unwrap();
    response.rotate(plane, axis).unwrap();
    response
}

fn coherent(
    responses: &[AbsorbingSlabResponse],
    region: SlabRegion,
    point: [f64; 3],
) -> AbsorbingCrystalField {
    let mut sum = FieldSum::at_point(point).unwrap();
    for response in responses {
        response.add_waves(region, &mut sum).unwrap();
    }
    responses[0]
        .with_heating(region, sum.finish().unwrap())
        .unwrap()
}

fn near(a: f64, b: f64, floor: f64) {
    assert!(
        (a - b).abs() <= 4e-11 * b.abs().max(floor),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_coherent_film_groups_preserve_tensor_interference() {
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    let a = solved(electric, 0.6, [1.0, 0.0]);
    let opposite = solved(electric.map(|value| -value), 0.6, [1.0, 0.0]);
    for (region, z) in [
        (SlabRegion::Incident, -0.2e-6),
        (SlabRegion::Reflected, -0.2e-6),
        (SlabRegion::TotalIncidentSide, -0.2e-6),
        (SlabRegion::Crystal, 0.37 * a.thickness_m),
        (SlabRegion::Transmitted, 1.4 * a.thickness_m),
    ] {
        let point = [0.3e-6, -0.2e-6, z];
        let single = a.field_at(region, point).unwrap();
        let bright = coherent(&[a, a], region, point);
        let dark = coherent(&[a, opposite], region, point);
        for i in 0..3 {
            assert!(
                (bright.field.electric_v_m[i] - single.field.electric_v_m[i].scale(2.0)).abs()
                    < 2e-12
            );
            assert!(
                (bright.field.magnetic_a_m[i] - single.field.magnetic_a_m[i].scale(2.0)).abs()
                    < 2e-14
            );
            assert!(dark.field.electric_v_m[i].abs() < 2e-12);
            assert!(dark.field.magnetic_a_m[i].abs() < 2e-14);
            near(
                bright.field.poynting_w_m2[i],
                4.0 * single.field.poynting_w_m2[i],
                1e-6,
            );
        }
        let heat = single.absorbed_power_density_w_m3;
        near(bright.absorbed_power_density_w_m3, 4.0 * heat, 1.0);
        assert!(dark.absorbed_power_density_w_m3 < 1e-20 * heat.max(1.0));
        let mut independent = crystal::StatisticsSum::new();
        independent.add(single);
        independent.add(opposite.field_at(region, point).unwrap());
        near(
            independent.finish().unwrap().absorbed_power_density_w_m3,
            2.0 * heat,
            1.0,
        );
    }
}

#[test]
fn g3_independent_film_factors_are_unitary_basis_invariant() {
    // Equal TRANSVERSE p/s fields at c=0.8 have tangential amplitudes 0.8,1.
    let linear =
        [[C64::from_re(0.8), C64::ZERO], [C64::ZERO, C64::ONE]].map(|e| solved(e, 0.6, [1.0, 0.0]));
    let circular = [
        [C64::from_re(0.8), C64::new(0.0, 1.0)],
        [C64::new(0.0, 0.8), C64::ONE],
    ]
    .map(|e| {
        solved(
            e.map(|z| z.scale(core::f64::consts::FRAC_1_SQRT_2)),
            0.6,
            [1.0, 0.0],
        )
    });
    for z in [0.0, 0.3 * linear[0].thickness_m, linear[0].thickness_m] {
        let statistics = |responses: &[AbsorbingSlabResponse]| {
            let mut sum = crystal::StatisticsSum::new();
            for response in responses {
                sum.add(
                    response
                        .field_at(SlabRegion::Crystal, [0.0, 0.0, z])
                        .unwrap(),
                );
            }
            sum.finish().unwrap()
        };
        let a = statistics(&linear);
        let b = statistics(&circular);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (a.electric_coherency_v2_m2[i][j] - b.electric_coherency_v2_m2[i][j]).abs()
                        < 3e-12
                );
            }
            near(
                a.magnetic_rms_squared_a2_m2[i],
                b.magnetic_rms_squared_a2_m2[i],
                1e-8,
            );
            near(a.poynting_w_m2[i], b.poynting_w_m2[i], 1e-6);
        }
        assert!(a.absorbed_power_density_w_m3 > 0.0);
        near(
            a.absorbed_power_density_w_m3,
            b.absorbed_power_density_w_m3,
            1.0,
        );
    }
}

#[test]
fn g1_crossed_film_beams_use_total_global_tensor_work() {
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    let a = solved(electric, 0.4, [0.6, 0.8]);
    let b = solved(electric.map(|z| z * C64::new(0.3, 0.7)), -0.3, [0.8, -0.6]);
    let norm = a.axis.iter().map(|v| v * v).sum::<f64>();
    let mut interferes = false;
    for x in [-0.8e-6, 0.0, 0.7e-6] {
        let point = [x, 0.1e-6, 0.6 * a.thickness_m];
        let fa = a.field_at(SlabRegion::Crystal, point).unwrap();
        let fb = b.field_at(SlabRegion::Crystal, point).unwrap();
        let actual = coherent(&[a, b], SlabRegion::Crystal, point);
        let e: [C64; 3] =
            core::array::from_fn(|i| fa.field.electric_v_m[i] + fb.field.electric_v_m[i]);
        let h: [C64; 3] =
            core::array::from_fn(|i| fa.field.magnetic_a_m[i] + fb.field.magnetic_a_m[i]);
        let dot = (0..3).fold(C64::ZERO, |sum, i| sum + e[i].scale(a.axis[i]));
        let parallel = dot.norm_sq() / norm;
        let perpendicular = (0..3)
            .map(|i| (e[i] - dot.scale(a.axis[i] / norm)).norm_sq())
            .sum::<f64>();
        let heat = 0.5 * a.omega * (a.loss_f_m[0] * perpendicular + a.loss_f_m[1] * parallel);
        near(actual.absorbed_power_density_w_m3, heat, 1.0);
        for i in 0..3 {
            assert!((actual.field.electric_v_m[i] - e[i]).abs() < 3e-12);
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            near(
                actual.field.poynting_w_m2[i],
                0.5 * (e[j] * h[k].conj() - e[k] * h[j].conj()).re,
                1e-6,
            );
        }
        interferes |= (heat - fa.absorbed_power_density_w_m3 - fb.absorbed_power_density_w_m3)
            .abs()
            > 1e-3 * heat;
    }
    assert!(interferes);
}
