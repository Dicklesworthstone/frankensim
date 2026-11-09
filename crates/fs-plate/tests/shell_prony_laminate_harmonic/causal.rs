//! Causal CST/DKT layer histories, original reference volumes and real heat.
use super::*;
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl, TensorPlaneStressHarmonicBudget, TensorStressControlBudget,
};
use fs_math::dd::Dd;
use fs_plate::shell::viscosity::harmonic::{
    TensorShellMaterialBudget, TensorShellMaterialCausalBudget,
};

const FREQUENCIES: [f64; 2] = [2.4, 3.17];
fn causal_budget() -> TensorShellMaterialCausalBudget {
    TensorShellMaterialCausalBudget {
        field: TensorShellMaterialBudget {
            maximum_points: 24,
            maximum_payload_bytes: 1 << 24,
            maximum_work: 2_000_000_000,
        },
        maximum_frequencies: 4,
        control: TensorStressControlBudget::default(),
        temporal: TensorPlaneStressHarmonicBudget {
            maximum_payload_bytes: 4_000_000,
            maximum_scalar_work: 5_000_000,
        },
    }
}
fn motion() -> [[C64; 6]; 2] {
    let first = [2e-4, -3e-4, 1.3e-4, -2.1e-4, 0.8e-4, 0.0].map(|v| C64::new(0.0, v));
    [first, first.map(|v| v.scale(0.4))]
}

// Independent physical free-node CST gradients and public DKT curvature, then
// the layer's signed z and engineering-strain rotation. No prepared observer
// rows, point strains or internal history fields enter this calculation.
fn principal(q: &[C64; 6], bary: [f64; 3], z: f64, angle: f64) -> [C64; 3] {
    let (rows, _) = dkt_curvature_quadrature(&[0.0, 0.4, 0.0], &[0.0, 0.0, 0.3], 0).unwrap();
    let mid = rows.map(|r| {
        core::array::from_fn::<_, 3, _>(|i| {
            q[2].scale(r[i][6]) + q[3].scale(r[i][8]) - q[4].scale(r[i][7])
        })
    });
    let curvature: [C64; 3] = core::array::from_fn(|i| {
        (0..3).fold(C64::ZERO, |sum, v| {
            sum + (mid[v][i] + mid[(v + 2) % 3][i] - mid[(v + 1) % 3][i]).scale(bary[v])
        })
    });
    let e = [
        curvature[0].scale(z),
        q[1].scale(1.0 / 0.3) + curvature[1].scale(z),
        q[0].scale(1.0 / 0.3) + curvature[2].scale(z),
    ];
    let (s, c) = angle.sin_cos();
    [
        e[0].scale(c * c) + e[1].scale(s * s) + e[2].scale(c * s),
        e[0].scale(s * s) + e[1].scale(c * c) - e[2].scale(c * s),
        e[0].scale(-2.0 * c * s) + e[1].scale(2.0 * c * s) + e[2].scale(c * c - s * s),
    ]
}

// DIFFERENT positive spatial rule (triangle vertices+centroid, layer Simpson)
// and an independent temporal owner: exact original mixed-control linear ramps
// approximate the prescribed sine, never the new harmonic temporal kernels.
fn ramp_reference(section: &LaminateSection, q: &[[C64; 6]; 2], duration: f64) -> [f64; 3] {
    let triangle = [
        ([1.0, 0.0, 0.0], 1.0 / 12.0),
        ([0.0, 1.0, 0.0], 1.0 / 12.0),
        ([0.0, 0.0, 1.0], 1.0 / 12.0),
        ([1.0 / 3.0; 3], 0.75),
    ];
    let mut total = [0.0; 3];
    for (layer, geometry) in section.layers().iter().zip(section.layer_integrals()) {
        let law = match &layer.law {
            LayerLaw::Prony(law) => law.clone(),
            LayerLaw::Elastic(law) => {
                TensorProny::new_mandel(law.stiffness(), vec![], 0.03).unwrap()
            }
            LayerLaw::PlaneStress { .. } => unreachable!("fixture retains full tensor sources"),
        };
        for (bary, area_weight) in triangle {
            for (fraction, height_weight) in [(0.0, 1.0 / 6.0), (0.5, 2.0 / 3.0), (1.0, 1.0 / 6.0)]
            {
                let z = geometry.z_bottom_m + fraction * layer.thickness_m;
                let strains = q.map(|q| principal(&q, bary, z, layer.orientation_rad));
                let mut history = law.history();
                let control = PreparedTensorMixedControl::new(
                    &history,
                    [false, false, true, false, true, true],
                    TensorStressControlBudget::default(),
                    || false,
                )
                .unwrap();
                let steps = 1024;
                for i in 1..=steps {
                    let t = duration * i as f64 / steps as f64;
                    let e: [f64; 3] = core::array::from_fn(|j| {
                        (0..2)
                            .map(|k| strains[k][j].im * (FREQUENCIES[k] * t).sin())
                            .sum()
                    });
                    control
                        .advance(
                            &mut history,
                            [e[0], e[1], 0.0, 0.5 * e[2], 0.0, 0.0],
                            duration / steps as f64,
                            || false,
                        )
                        .unwrap();
                }
                let volume = 0.06 * area_weight * height_weight * layer.thickness_m;
                for (sum, value) in total.iter_mut().zip([
                    history.stored_energy(),
                    history.dissipated(),
                    history.work(),
                ]) {
                    *sum += volume * value;
                }
            }
        }
    }
    total
}

#[test]
fn startup_matches_original_layer_ramps_and_independent_positive_spatial_rule() {
    with_cx(|cx| {
        let section = laminate();
        let shell = source(section.clone(), CstDktMass::Consistent, cx).unwrap();
        let observer = shell
            .prepare_material_causal_observer(causal_budget(), cx)
            .unwrap();
        let initial = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let q = motion();
        let refs = q.each_ref().map(|q| q.as_slice());
        let step = observer
            .advance_window(&initial, &FREQUENCIES, &refs, &[0.0; 6], 0.29, cx)
            .unwrap();
        let expected = ramp_reference(&section, &q, 0.29);
        for (actual, reference) in [
            step.field.end_stored_energy_j,
            step.field.material_heat_j,
            step.field.stress_work_j,
        ]
        .into_iter()
        .zip(expected)
        {
            real_near(actual, reference, 3e-6);
        }
        assert!(step.field.material_heat_j > 0.0);
        assert_eq!(step.field.kelvin_heat_j, 0.0);
        assert!(
            step.field.energy_balance_residual_j.abs() < 3e-10 * step.field.stress_work_j.abs()
        );
        assert_eq!(step.field.points.len(), 12);
        let mut volumes = [0.0; 2];
        for point in &step.field.points {
            let layer = &section.layers()[point.layer];
            volumes[point.layer] += point.volume_m3;
            for k in 0..2 {
                let strain = principal(
                    &q[k],
                    point.triangle_barycentric,
                    point.thickness_position_m,
                    layer.orientation_rad,
                );
                for j in 0..3 {
                    near(point.in_plane_strain[k][j], strain[j], 3e-12);
                }
            }
            for j in [2, 4, 5] {
                assert_eq!(point.end_stress[j], 0.0);
            }
            if point.layer == 0 {
                assert_eq!(point.material.material_heat_j_m3, 0.0);
            }
        }
        real_near(volumes[0], 0.06 * 0.08, 3e-13);
        real_near(volumes[1], 0.06 * 0.02, 3e-13);
        assert_eq!(initial.full_displacement(), &[0.0; 18]);
        assert_eq!(initial.stored_energy_j(), 0.0);
        assert_eq!(initial.time_s_dd(), Dd::ZERO);
        assert!(core::ptr::eq(step.history.source(), &shell));
    });
}

#[test]
fn partitions_and_empty_carrier_holds_preserve_actual_layer_memory() {
    with_cx(|cx| {
        let shell = source(laminate(), CstDktMass::Consistent, cx).unwrap();
        let observer = shell
            .prepare_material_causal_observer(causal_budget(), cx)
            .unwrap();
        let history = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let q = motion();
        let refs = q.each_ref().map(|q| q.as_slice());
        let whole = observer
            .advance_window(&history, &FREQUENCIES, &refs, &[0.0; 6], 0.375, cx)
            .unwrap();
        let left = observer
            .advance_window(&history, &FREQUENCIES, &refs, &[0.0; 6], 0.125, cx)
            .unwrap();
        let right = observer
            .advance_window(&left.history, &FREQUENCIES, &refs, &[0.0; 6], 0.25, cx)
            .unwrap();
        assert_eq!(whole.history.time_s_dd(), right.history.time_s_dd());
        assert_eq!(&right.history.full_displacement()[..12], &[0.0; 12]);
        for (a, b) in whole
            .history
            .full_displacement()
            .iter()
            .zip(right.history.full_displacement())
        {
            real_near(*a, *b, 3e-12);
        }
        for (a, b) in [
            (
                whole.field.material_heat_j,
                left.field.material_heat_j + right.field.material_heat_j,
            ),
            (
                whole.field.stress_work_j,
                left.field.stress_work_j + right.field.stress_work_j,
            ),
            (
                whole.field.stored_energy_change_j,
                left.field.stored_energy_change_j + right.field.stored_energy_change_j,
            ),
        ] {
            real_near(a, b, 3e-10);
        }
        let held_bias = &right.history.full_displacement()[12..];
        let held = observer
            .advance_window(&right.history, &[], &[], held_bias, 0.2, cx)
            .unwrap();
        assert!(
            held.field.material_heat_j > 0.0,
            "empty AC must retain prior Maxwell heat"
        );
        assert_eq!(held.field.stress_work_j, 0.0);
        real_near(
            held.field.stored_energy_change_j,
            -held.field.material_heat_j,
            3e-10,
        );
        real_near(
            held.history.dissipated_energy_j(),
            right.history.dissipated_energy_j() + held.field.material_heat_j,
            3e-12,
        );
        assert_eq!(
            held.history.full_displacement(),
            right.history.full_displacement()
        );
        assert_eq!(
            held.history.time_s_dd(),
            right.history.time_s_dd() + Dd::from_f64(0.2)
        );
        assert!(
            observer
                .advance_window(&right.history, &[], &[], &[0.0; 6], 0.2, cx)
                .is_err()
        );
        assert_eq!(history.dissipated_energy_j(), 0.0);
    });
}

fn elastic_norm(strain: [C64; 3], law: &OrthotropicElastic) -> f64 {
    let c = law.stiffness();
    let e = strain.map(|v| v.im);
    let zz = -(c[2][0] * e[0] + c[2][1] * e[1]) / c[2][2];
    (e[0] * e[0] + e[1] * e[1] + zz * zz + 0.5 * e[2] * e[2]).sqrt()
}

#[test]
fn layer_face_refusals_and_cancellation_publish_no_partial_history() {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let key = fs_exec::StreamKey {
            seed: 59,
            kernel_id: 69,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(
            &gate,
            arena,
            key,
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let elastic = elastic(0.03);
        let section = LaminateSection::try_new(
            vec![layer(LayerLaw::Elastic(elastic.clone()), 0.1, 4.0, 0.31)],
            ReferenceSurface::Offset(0.03),
        )
        .unwrap();
        let shell = source(section.clone(), CstDktMass::Consistent, &cx).unwrap();
        let observer = shell
            .prepare_material_causal_observer(causal_budget(), &cx)
            .unwrap();
        let initial = observer.start_relaxed(Dd::ZERO, &cx).unwrap();
        let mut q = [C64::ZERO; 6];
        q[3] = C64::new(0.0, 1.0);
        let geometry = &section.layer_integrals()[0];
        let (mut boundary, mut interior) = (0.0_f64, 0.0_f64);
        for bary in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            for z in [geometry.z_bottom_m, geometry.z_top_m] {
                boundary = boundary.max(elastic_norm(principal(&q, bary, z, 0.31), &elastic));
            }
        }
        let center = 0.5 * (geometry.z_bottom_m + geometry.z_top_m);
        for bary in [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] {
            for sign in [-1.0, 1.0] {
                let z = center + sign * 0.1 / (2.0 * 3.0_f64.sqrt());
                interior = interior.max(elastic_norm(principal(&q, bary, z, 0.31), &elastic));
            }
        }
        assert!(boundary > 1.01 * interior);
        let scale = 0.03 / (boundary * interior).sqrt();
        assert!(scale * interior < 0.03 && scale * boundary > 0.03);
        assert!(
            observer
                .advance_window(
                    &initial,
                    &[2.4],
                    &[&q.map(|v| v.scale(scale))],
                    &[0.0; 6],
                    0.01,
                    &cx
                )
                .is_err()
        );
        let safe = motion();
        let refs = safe.each_ref().map(|q| q.as_slice());
        let accepted = observer
            .advance_window(&initial, &FREQUENCIES, &refs, &[0.0; 6], 0.2, &cx)
            .unwrap();
        let other = shell
            .prepare_material_causal_observer(causal_budget(), &cx)
            .unwrap();
        assert!(
            other
                .advance_window(&initial, &FREQUENCIES, &refs, &[0.0; 6], 0.2, &cx)
                .is_err()
        );
        assert!(
            observer
                .advance_window(
                    &initial,
                    &FREQUENCIES,
                    &[&safe[0][..5], &safe[1]],
                    &[0.0; 6],
                    0.2,
                    &cx
                )
                .is_err()
        );
        // Drilling produces no material strain. Exact duplicate cancellation
        // must also remove huge phasors from the full-nodal continuity scale;
        // otherwise a spurious static drilling jump could hide below it.
        let mut huge = [C64::ZERO; 6];
        huge[5] = C64::new(1e100, -2e100);
        let opposite = huge.map(|v| -v);
        let mut jump = [0.0; 6];
        jump[5] = 1e-4;
        assert!(
            observer
                .advance_window(&initial, &[2.4; 2], &[&huge, &opposite], &jump, 0.2, &cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            observer.advance_window(&initial, &FREQUENCIES, &refs, &[0.0; 6], 0.2, &cx),
            Err(TensorShellHarmonicError::Plate(
                fs_plate::PlateError::Cancelled
            ))
        ));
        assert_eq!(initial.time_s_dd(), Dd::ZERO);
        assert_eq!(initial.full_displacement(), &[0.0; 18]);
        assert_eq!(initial.stress_work_j(), 0.0);
        let fresh = CancelGate::new_clock_free();
        let fresh_cx = Cx::new(
            &fresh,
            arena,
            key,
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let retry = observer
            .advance_window(&initial, &FREQUENCIES, &refs, &[0.0; 6], 0.2, &fresh_cx)
            .unwrap();
        assert_eq!(retry.field, accepted.field);
    });
}
