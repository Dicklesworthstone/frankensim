//! Native layered-shell fatigue, using synthetic data and independent references.
//! These fixtures check implementation, not measured component lifetime.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::miner::CriticalPlaneMinerLimits;
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::fatigue::strain_life::source::{
    CmbPowerBasis, CmbProperty, CmbSourceRequest, resolve_cmb,
};
use fs_material::fatigue::strain_life::{
    CmbParameters, PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainStressSample,
};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};
use fs_material::plastic::source::{
    HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_plate::shell::{DRILLING_ALPHA, ShellMesh};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::IncrementSettings;
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;
use fs_solid::plastic_shell::fatigue::{
    ShellFatigueConfig, ShellFatigueError, ShellFatigueGauge, ShellFatiguePath,
};
use fs_solid::plastic_shell::{PlasticShellBudget, PlasticShellSections, SmallStrainShellProblem};

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const THICKNESS: f64 = 0.1;
const SIGMA_F: f64 = 10.0;
const EPSILON_F: f64 = 0.04;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 101,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} ({tolerance:e})"
    );
}
fn temperature() -> QuantitySpec {
    QuantitySpec::dimensional(Temperature::DIMS)
}
fn point(quantity: QuantitySpec) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity, 300.0)
        .unwrap()
}
fn request() -> CmbSourceRequest {
    CmbSourceRequest {
        power_basis: CmbPowerBasis::Reversals,
        temperature_axis: "T".into(),
        stress_ratio_axis: "calibration-stress-ratio".into(),
        surface_and_geometry: "synthetic smooth shell coupon".into(),
        environment: "synthetic fixed environment".into(),
        load_basis: "total mechanical normal strain amplitude".into(),
    }
}
fn card(revision: u32, quantity: QuantitySpec) -> MaterialCard {
    card_with_properties(revision, quantity, &[])
}
fn card_with_properties(
    revision: u32,
    quantity: QuantitySpec,
    additional: &[(&str, Dims, f64)],
) -> MaterialCard {
    let thermal = ValidityDomain::unconstrained().with_quantity("T", quantity, 300.0, 300.0);
    let mut claims = ClaimSet::new();
    let request = request();
    for property in CmbProperty::ALL {
        let key = property.property_key();
        let value = [E, SIGMA_F, -0.5, EPSILON_F, -0.5][property.index()];
        let validity = if property == CmbProperty::YoungsModulus {
            thermal.clone()
        } else {
            thermal
                .clone()
                .with_quantity(
                    &request.stress_ratio_axis,
                    QuantitySpec::dimensional(Dims::NONE),
                    -1.0,
                    -1.0,
                )
                .with_quantity(
                    request.power_basis.life_axis(),
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                    2e12,
                )
        };
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: if property == CmbProperty::YoungsModulus {
                        "synthetic elastic measurement"
                    } else {
                        "one synthetic CMB fit"
                    }
                    .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    for (name, dims, value) in [
        ("poisson_ratio", Dims::NONE, NU),
        ("yield_stress", Pressure::DIMS, 1.0),
        ("isotropic_hardening_modulus", Pressure::DIMS, 20.0),
        ("kinematic_hardening_modulus", Pressure::DIMS, 80.0),
    ]
    .into_iter()
    .chain(additional.iter().copied())
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: thermal.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic cyclic shell material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shell alloy".into(),
            phase: "solid".into(),
            process: "synthetic cyclic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    card
}
fn material(card: &MaterialCard, quantity: QuantitySpec) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &point(quantity),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
    )
    .unwrap()
}
fn sourced_curve(card: &MaterialCard) -> Arc<StrainLifeCurve> {
    let request = request();
    let query = point(temperature())
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
        )
        .unwrap()
        .with_quantity(
            request.power_basis.life_axis(),
            QuantitySpec::dimensional(Dims::NONE),
            100.0,
        )
        .unwrap();
    Arc::new(
        resolve_cmb(
            card,
            &query,
            &request,
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 5],
        )
        .unwrap(),
    )
}
fn curve(minimum_cycles: f64) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic equal-exponent CMB reference, not measured".into(),
                material: "synthetic shell alloy".into(),
                process: "synthetic cyclic coupon".into(),
                surface_and_geometry: "synthetic smooth shell coupon".into(),
                environment: "synthetic fixed environment".into(),
                load_basis: "total mechanical normal strain amplitude".into(),
                temperature_k: 300.0,
            },
            CmbParameters {
                youngs_modulus_pa: E,
                fatigue_strength_coefficient_pa: SIGMA_F,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: EPSILON_F,
                fatigue_ductility_exponent: -0.5,
            },
            minimum_cycles,
            1e12,
        )
        .unwrap(),
    )
}
fn section<L: SourceBoundFatigueLaw>(law: L) -> LayeredPlaneStressSection<L> {
    LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.5 * THICKNESS,
            top_z_m: 0.5 * THICKNESS,
            subdivisions: 1,
            law,
        }],
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 48,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.05,
        },
        SectionBudget {
            maximum_layers: 2,
            maximum_quadrature_points: 12,
        },
    )
    .unwrap()
}
fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
fn problem<'a, L: SourceBoundFatigueLaw>(
    mesh: &'a ShellMesh,
    section: &'a LayeredPlaneStressSection<L>,
    prescribed: &'a [(usize, f64)],
    loads: &'a [(usize, f64)],
    factors: &'a [f64],
) -> SmallStrainShellProblem<'a, L> {
    SmallStrainShellProblem {
        mesh,
        sections: PlasticShellSections::Uniform(section),
        coordinate_length_m: 0.4,
        rotation_limit_rad: 0.05,
        drilling_alpha: DRILLING_ALPHA,
        prescribed,
        nodal_loads: loads,
        load_path: factors,
        budget: PlasticShellBudget::default(),
        settings: IncrementSettings {
            max_newton_iterations: 48,
            ..IncrementSettings::default()
        },
    }
}
fn limits() -> CriticalPlaneMinerLimits {
    CriticalPlaneMinerLimits {
        rainflow: RainflowLimits {
            maximum_samples: 128,
            maximum_residue: 24,
        },
        maximum_planes: 8,
        maximum_projection_work: 1024,
    }
}
fn config() -> ShellFatigueConfig {
    ShellFatigueConfig {
        maximum_gauges: 8,
        maximum_total_residue: 1024,
        maximum_total_projection_work: 4096,
    }
}
fn gauge<L: SourceBoundFatigueLaw>(
    law: &L,
    point: usize,
    policy: StrainLifeMeanStress,
) -> ShellFatigueGauge {
    ShellFatigueGauge {
        element: 0,
        surface_point: 0,
        thickness_point: point,
        material_state_identity: law.fatigue_source().identity(),
        curve: curve(0.5),
        mean_stress: policy,
        normals_in_facet_frame: vec![[1.0, 0.0, 0.0]],
        limits: limits(),
    }
}

#[path = "plastic_shell_fatigue/cyclic.rs"]
mod cyclic;
#[path = "plastic_shell_fatigue/loading.rs"]
mod loading;
#[path = "plastic_shell_fatigue/shear.rs"]
mod shear;
fn membrane_prescribed() -> Vec<(usize, f64)> {
    // epsilon_x is prescribed, but the free y translation supplies the actual
    // uniaxial Poisson/plastic contraction through native shell equilibrium.
    (0..18)
        .filter(|&d| d != 13)
        .map(|d| (d, if d == 6 { 1.0 } else { 0.0 }))
        .collect()
}
fn elastic_field(mesh: &ShellMesh, epsilon: f64, curvature: f64) -> Vec<(usize, f64)> {
    let frame = mesh.facet(0).unwrap().frame;
    let origin = mesh.nodes[0];
    let mut prescribed = Vec::new();
    for (node, point) in mesh.nodes.iter().enumerate() {
        let x: f64 = (0..3).map(|c| frame[0][c] * (point[c] - origin[c])).sum();
        let y: f64 = (0..3).map(|c| frame[1][c] * (point[c] - origin[c])).sum();
        let ky = -NU * curvature;
        let displacement = [
            epsilon * x,
            -NU * epsilon * y,
            -0.5 * (curvature * x * x + ky * y * y),
        ];
        let rotation = [-ky * y, curvature * x, 0.0];
        for c in 0..3 {
            prescribed.push((
                6 * node + c,
                (0..3).map(|a| frame[a][c] * displacement[a]).sum(),
            ));
            prescribed.push((
                6 * node + 3 + c,
                (0..3).map(|a| frame[a][c] * rotation[a]).sum(),
            ));
        }
    }
    prescribed
}

#[test]
fn native_elastic_cycle_matches_hand_cmb_and_prefix_preview_preserves_reversals() {
    with_cx(&CancelGate::new(), |cx| {
        let card = card(0, temperature());
        let law = material(&card, temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let amplitude = 0.0008;
        let factors = [amplitude, 0.0];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut g = gauge(&law, 1, StrainLifeMeanStress::Uncorrected);
        g.curve = sourced_curve(&card);
        let mut run = ShellFatiguePath::new(&p, vec![g], config(), cx).unwrap();
        assert_eq!(run.sample_count(), 1);
        assert_eq!(run.next_step(), 0);
        assert!(run.curve(0).unwrap().source().is_some());
        assert_eq!(
            run.material_source(0).unwrap().identity(),
            law.material().identity()
        );
        assert_eq!(run.finish(cx), Err(ShellFatigueError::PathIncomplete));
        assert!(run.advance_next(cx).unwrap());
        close(
            run.samples()[0].total_mechanical_strain[0],
            amplitude,
            1e-14,
        );
        close(run.samples()[0].stress_pa[0], E * amplitude, 1e-10);
        let before = run.plane_summaries(0).unwrap();
        let prefix = run.preview_finish(cx).unwrap();
        assert!(prefix[0].summary.finished);
        assert_eq!(run.plane_summaries(0).unwrap(), before);
        assert_eq!(run.results(), vec![None]);
        assert!(!run.is_finished());
        let coefficient = SIGMA_F / E + EPSILON_F;
        let expected = 2.0 * (0.5 * amplitude / coefficient).powi(2);
        close(prefix[0].summary.miner_damage, 0.5 * expected, 1e-14);
        let mut resumed = run.clone();
        assert!(run.advance_next(cx).unwrap());
        assert!(resumed.advance_next(cx).unwrap());
        assert!(!run.advance_next(cx).unwrap());
        let preview = run.preview_finish(cx).unwrap();
        let result = run.finish(cx).unwrap();
        assert_eq!(result, preview);
        assert_eq!(result, resumed.finish(cx).unwrap());
        assert_eq!(run.solution(), resumed.solution());
        assert_eq!(result[0].samples, 3);
        close(result[0].summary.miner_damage, expected, 1e-14);
        let mut direct = StrainLifeMiner::new(
            curve(0.5),
            StrainLifeMeanStress::Uncorrected,
            limits().rainflow,
        )
        .unwrap();
        for strain in [0.0, amplitude, 0.0] {
            direct.push(strain).unwrap();
        }
        close(
            result[0].summary.miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-14,
        );
        assert_eq!(run.finish(cx).unwrap(), result);
        assert_eq!(run.advance_next(cx), Err(ShellFatigueError::Finished));
    });
}

#[test]
fn cyclic_shell_j2_stress_and_morrow_match_independent_prager_history() {
    with_cx(&CancelGate::new(), |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let factors = [0.006, -0.006, 0.006, -0.006, 0.006];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut run = ShellFatiguePath::new(
            &p,
            vec![gauge(&law, 0, StrainLifeMeanStress::Morrow)],
            config(),
            cx,
        )
        .unwrap();
        let mut direct =
            PairedStrainLifeMiner::new(curve(0.5), StrainLifeMeanStress::Morrow, limits().rainflow)
                .unwrap();
        direct
            .push(StrainStressSample {
                mechanical_strain: 0.0,
                normal_stress_pa: 0.0,
            })
            .unwrap();
        let (mut plastic, mut accumulated) = (0.0_f64, 0.0_f64);
        for strain in factors {
            let relative_trial = E * (strain - plastic) - 80.0 * plastic;
            let excess = relative_trial.abs() - (1.0 + 20.0 * accumulated);
            if excess > 0.0 {
                let increment = excess / (E + 20.0 + 80.0);
                plastic += relative_trial.signum() * increment;
                accumulated += increment;
            }
            let stress = E * (strain - plastic);
            run.advance_next(cx).unwrap();
            let sample = run.samples()[0];
            close(sample.total_mechanical_strain[0], strain, 1e-12);
            close(
                sample.total_mechanical_strain[1],
                -NU * stress / E - 0.5 * plastic,
                2e-10,
            );
            close(
                sample.total_mechanical_strain[2],
                sample.total_mechanical_strain[1],
                2e-10,
            );
            close(sample.stress_pa[0], stress, 2e-8);
            close(sample.stress_pa[1], 0.0, 2e-8);
            close(sample.stress_pa[2], 0.0, 2e-8);
            assert!((sample.stress_pa[0] - E * strain).abs() > 0.5);
            for state in &run.solution().section_states {
                for point in &state.points {
                    close(point.history.plastic_strain[0], plastic, 2e-10);
                    close(point.history.alpha, accumulated, 2e-10);
                }
            }
            direct
                .push(StrainStressSample {
                    mechanical_strain: strain,
                    normal_stress_pa: stress,
                })
                .unwrap();
        }
        assert!(accumulated > 0.02);
        close(
            run.finish(cx).unwrap()[0].summary.miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-10,
        );
        let native = p.solve(cx).unwrap();
        assert_eq!(run.solution().coordinates, native.coordinates);
        assert_eq!(run.solution().section_states, native.section_states);
    });
}

#[test]
fn finite_plane_search_uses_native_facet_axes_and_selects_largest_damage() {
    with_cx(&CancelGate::new(), |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let flat = mesh();
        let rotated = ShellMesh::new(
            vec![[0.2, -0.4, 0.7], [0.2, 0.6, 0.7], [0.2, -0.4, 1.7]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let amplitude = 0.0006;
        for mesh in [&flat, &rotated] {
            let fixed = elastic_field(mesh, amplitude, 0.0);
            let factors = [1.0, -1.0, 1.0, 0.0];
            let p = problem(mesh, &section, &fixed, &[], &factors);
            let mut g = gauge(&law, 2, StrainLifeMeanStress::Uncorrected);
            g.normals_in_facet_frame = vec![
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let mut run = ShellFatiguePath::new(&p, vec![g], config(), cx).unwrap();
            let mut direct = StrainLifeMiner::new(
                curve(0.5),
                StrainLifeMeanStress::Uncorrected,
                limits().rainflow,
            )
            .unwrap();
            direct.push(0.0).unwrap();
            for factor in factors {
                run.advance_next(cx).unwrap();
                direct.push(factor * amplitude).unwrap();
                close(
                    run.samples()[0].total_mechanical_strain[0],
                    factor * amplitude,
                    1e-14,
                );
            }
            let result = &run.finish(cx).unwrap()[0];
            assert_eq!(result.selected_plane, 1);
            assert_eq!(result.normal_in_frame, [1.0, 0.0, 0.0]);
            assert_eq!(result.evaluated_planes, 4);
            assert_eq!(result.projection_work, 20);
            close(
                result.summary.miner_damage,
                direct.finish().unwrap().miner_damage,
                1e-13,
            );
            let summaries = run.plane_summaries(0).unwrap();
            close(
                summaries[0].miner_damage,
                NU * NU * summaries[1].miner_damage,
                1e-13,
            );
        }
    });
}

#[test]
fn bending_gauges_keep_original_depth_and_signed_mean_stress() {
    with_cx(&CancelGate::new(), |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let curvature = 0.012;
        let fixed = elastic_field(&mesh, 0.0, curvature);
        let p = problem(&mesh, &section, &fixed, &[], &[1.0, 0.0]);
        let bottom = gauge(&law, 0, StrainLifeMeanStress::Morrow);
        let mut top = gauge(&law, 2, StrainLifeMeanStress::Morrow);
        top.surface_point = 2;
        let mut run = ShellFatiguePath::new(&p, vec![bottom, top], config(), cx).unwrap();
        run.advance_next(cx).unwrap();
        for (index, thickness_point) in [0, 2].into_iter().enumerate() {
            let strain = curvature * section.points()[thickness_point].z_m;
            close(
                run.samples()[index].total_mechanical_strain[0],
                strain,
                1e-14,
            );
            close(run.samples()[index].stress_pa[0], E * strain, 1e-10);
        }
        assert!(run.samples()[0].stress_pa[0] < 0.0);
        assert!(run.samples()[1].stress_pa[0] > 0.0);
        run.advance_next(cx).unwrap();
        let results = run.finish(cx).unwrap();
        for (index, thickness_point) in [0, 2].into_iter().enumerate() {
            let peak_strain = curvature * section.points()[thickness_point].z_m;
            let coefficient = (SIGMA_F - 0.5 * E * peak_strain) / E + EPSILON_F;
            let expected = 2.0 * (0.5 * peak_strain.abs() / coefficient).powi(2);
            close(results[index].summary.miner_damage, expected, 1e-13);
        }
        assert!(results[1].summary.miner_damage > results[0].summary.miner_damage);
    });
}

#[test]
fn late_gauge_refusal_and_cancellation_preserve_native_plastic_state_and_counters() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let factors = [0.006, -0.006, 0.006];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let wide = gauge(&law, 0, StrainLifeMeanStress::Uncorrected);
        let mut narrow = gauge(&law, 2, StrainLifeMeanStress::Uncorrected);
        narrow.curve = curve(1000.0);
        let mut run = ShellFatiguePath::new(&p, vec![wide, narrow], config(), cx).unwrap();
        run.advance_next(cx).unwrap();
        let saved = run.clone();
        for _ in 0..2 {
            assert!(matches!(
                run.advance_next(cx),
                Err(ShellFatigueError::Gauge { gauge: 1, .. })
            ));
            assert_eq!(run.solution(), saved.solution());
            assert_eq!(run.samples(), saved.samples());
            assert_eq!(run.velocity(), saved.velocity());
            assert_eq!(run.sample_count(), 2);
            assert_eq!(run.next_step(), 1);
            for gauge in 0..2 {
                assert_eq!(run.plane_summaries(gauge), saved.plane_summaries(gauge));
            }
        }
        gate.request();
        assert_eq!(run.advance_next(cx), Err(ShellFatigueError::Cancelled));
        assert_eq!(run.preview_finish(cx), Err(ShellFatigueError::Cancelled));
        assert_eq!(run.solution(), saved.solution());
        assert_eq!(run.results(), saved.results());
    });
}

#[test]
fn source_receipts_absolute_temperature_and_actual_layer_point_are_required() {
    with_cx(&CancelGate::new(), |cx| {
        let card0 = card(0, temperature());
        let law0 = material(&card0, temperature());
        let card1 = card(1, temperature());
        let law1 = material(&card1, temperature());
        let section = LayeredPlaneStressSection::new(
            vec![
                PlaneStressLayer {
                    bottom_z_m: -0.05,
                    top_z_m: 0.0,
                    subdivisions: 1,
                    law: law0.clone(),
                },
                PlaneStressLayer {
                    bottom_z_m: 0.0,
                    top_z_m: 0.05,
                    subdivisions: 1,
                    law: law1.clone(),
                },
            ],
            ThicknessQuadrature::Gauss3,
            ReductionConfig {
                max_iterations: 48,
                stress_tolerance_pa: 1e-12,
                strain_limit: 0.05,
            },
            SectionBudget {
                maximum_layers: 2,
                maximum_quadrature_points: 12,
            },
        )
        .unwrap();
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let p = problem(&mesh, &section, &fixed, &[], &[0.0005, 0.0]);
        let mut correct = gauge(&law0, 0, StrainLifeMeanStress::Uncorrected);
        correct.curve = sourced_curve(&card0);
        let mut run = ShellFatiguePath::new(&p, vec![correct.clone()], config(), cx).unwrap();
        while run.advance_next(cx).unwrap() {}
        let expected = run.finish(cx).unwrap();
        assert!(expected[0].summary.miner_damage > 0.0);
        let borrowed_sections = [&section];
        let mut borrowed = problem(&mesh, &section, &fixed, &[], &[0.0005, 0.0]);
        borrowed.sections = PlasticShellSections::PerElementBorrowed(&borrowed_sections);
        let mut borrowed_run =
            ShellFatiguePath::new(&borrowed, vec![correct.clone()], config(), cx).unwrap();
        while borrowed_run.advance_next(cx).unwrap() {}
        assert_eq!(borrowed_run.solution(), run.solution());
        assert_eq!(borrowed_run.finish(cx).unwrap(), expected);
        let mut wrong = correct.clone();
        wrong.material_state_identity = ContentHash([91; 32]);
        assert!(matches!(
            ShellFatiguePath::new(&p, vec![wrong], config(), cx),
            Err(ShellFatigueError::SourceMismatch { .. })
        ));
        let mut wrong_layer = correct.clone();
        wrong_layer.thickness_point = 3;
        wrong_layer.material_state_identity = law1.material().identity();
        // Equal values and labels from another card do not bind the CMB fit.
        assert!(matches!(
            ShellFatiguePath::new(&p, vec![wrong_layer], config(), cx),
            Err(ShellFatigueError::SourceMismatch { .. })
        ));
        for (surface, depth) in [(3, 0), (0, section.point_count())] {
            let mut invalid = correct.clone();
            invalid.surface_point = surface;
            invalid.thickness_point = depth;
            assert!(matches!(
                ShellFatiguePath::new(&p, vec![invalid], config(), cx),
                Err(ShellFatigueError::Invalid(_))
            ));
        }
        // A "T" axis typed as a pressure is not an absolute temperature. The J2
        // source admits only fixed external conditions, so it refuses this card
        // before any shell section or gauge can bind it.
        let pressure = QuantitySpec::dimensional(Pressure::DIMS);
        assert!(matches!(
            resolve_j2_plasticity(
                &card(0, pressure),
                &point(pressure),
                MaterialPropertySelection::SingleClaimOnly,
                HardeningPropertyBasis::MixedPlasticStrainModuli,
                0.05,
            ),
            Err(J2SourceError::Invalid(reason)) if reason.contains("fixed external conditions")
        ));
    });
}

#[test]
fn combined_residue_and_complete_path_projection_bounds_include_initial_sample() {
    with_cx(&CancelGate::new(), |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let p = problem(&mesh, &section, &fixed, &[], &[0.0005, 0.0]);
        let mut g = gauge(&law, 0, StrainLifeMeanStress::Uncorrected);
        g.normals_in_facet_frame.push([0.0, 1.0, 0.0]);
        for controls in [
            ShellFatigueConfig {
                maximum_total_residue: 47,
                ..config()
            },
            ShellFatigueConfig {
                maximum_total_projection_work: 5,
                ..config()
            },
            ShellFatigueConfig {
                maximum_gauges: 0,
                ..config()
            },
        ] {
            assert!(matches!(
                ShellFatiguePath::new(&p, vec![g.clone()], controls, cx),
                Err(ShellFatigueError::Budget)
            ));
        }
        let exact = ShellFatigueConfig {
            maximum_total_residue: 48,
            maximum_total_projection_work: 6,
            ..config()
        };
        let mut run = ShellFatiguePath::new(&p, vec![g.clone()], exact, cx).unwrap();
        while run.advance_next(cx).unwrap() {}
        assert_eq!(run.finish(cx).unwrap()[0].projection_work, 6);
        g.limits.rainflow.maximum_residue = usize::MAX;
        assert!(matches!(
            ShellFatiguePath::new(&p, vec![g], config(), cx),
            Err(ShellFatigueError::Budget)
        ));
    });
}

#[test]
fn dynamic_motion_energy_and_prefix_resume_match_original_native_shell_steps() {
    with_cx(&CancelGate::new(), |cx| {
        let law = material(&card(0, temperature()), temperature());
        let section = section(law.clone());
        let mesh = mesh();
        let fixed: Vec<_> = (0..18)
            .filter(|&d| d != 6 && d != 13)
            .map(|d| (d, 0.0))
            .collect();
        let loads = [(6, 0.08)];
        let factors = [1.0, 0.0, 0.0, -0.5, 0.0, 0.5, 0.0, 0.0];
        let p = problem(&mesh, &section, &fixed, &loads, &factors);
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1.0), cx)
            .unwrap();
        let replacement = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(2.0), cx)
            .unwrap();
        let mut run = ShellFatiguePath::new(
            &p,
            vec![gauge(&law, 0, StrainLifeMeanStress::Uncorrected)],
            config(),
            cx,
        )
        .unwrap();
        // A real static plastic preload may start the same native dynamic path.
        run.advance_next(cx).unwrap();
        let mut native = problem(&mesh, &section, &fixed, &loads, &[1.0])
            .solve(cx)
            .unwrap();
        let mut velocity = vec![[0.0; 6]; 3];
        assert_eq!(run.solution(), &native);
        assert!(native.section_states[0].points[0].history.alpha > 0.0);
        let mut direct = StrainLifeMiner::new(
            curve(0.5),
            StrainLifeMeanStress::Uncorrected,
            limits().rainflow,
        )
        .unwrap();
        direct.push(0.0).unwrap();
        direct.push(native.coordinates[1][0]).unwrap();
        let dt = 0.01;
        let mut resumed = None;
        for (index, factor) in factors.iter().enumerate().skip(1) {
            let scaled_loads: Vec<_> = loads.iter().map(|&(d, v)| (d, factor * v)).collect();
            let endpoint = problem(&mesh, &section, &fixed, &scaled_loads, &[1.0]);
            let mut step = endpoint
                .advance_backward_euler(Some(&native), &velocity, &reference, dt, cx)
                .unwrap();
            let energy = run
                .advance_next_backward_euler(&reference, dt, cx)
                .unwrap()
                .unwrap();
            assert_eq!(energy, step.inertia);
            step.solution.increments[0].load_factor = *factor;
            assert_eq!(run.solution(), &step.solution);
            assert_eq!(run.velocity(), step.velocity);
            direct.push(step.solution.coordinates[1][0]).unwrap();
            native = step.solution;
            velocity = step.velocity;
            if index == 1 {
                let saved = run.clone();
                assert!(matches!(
                    run.advance_next(cx),
                    Err(ShellFatigueError::Invalid(_))
                ));
                assert!(matches!(
                    run.advance_next_backward_euler(&replacement, dt, cx),
                    Err(ShellFatigueError::Invalid(_))
                ));
                assert_eq!(run.solution(), saved.solution());
                assert_eq!(run.next_step(), saved.next_step());
                run.preview_finish(cx).unwrap();
                resumed = Some(saved);
            } else if let Some(restarted) = &mut resumed {
                assert_eq!(
                    restarted
                        .advance_next_backward_euler(&reference, dt, cx)
                        .unwrap(),
                    Some(energy)
                );
                assert_eq!(restarted.solution(), run.solution());
                assert_eq!(restarted.samples(), run.samples());
            }
        }
        close(run.elapsed_time_s(), 7.0 * dt, 1e-15);
        assert_eq!(run.sample_count(), 9);
        assert!(
            run.advance_next_backward_euler(&reference, dt, cx)
                .unwrap()
                .is_none()
        );
        let result = run.finish(cx).unwrap();
        assert_eq!(result, resumed.unwrap().finish(cx).unwrap());
        close(
            result[0].summary.miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-12,
        );
    });
}
