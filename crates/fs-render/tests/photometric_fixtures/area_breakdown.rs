//! Attribution of the original finite-area reference, including its numerical scope.
use super::*;
use fs_evidence::NumericalKind;
use fs_render::photometry::IlluminanceEstimate;
use fs_render::photometry::area::reference::occlusion::{
    AreaOcclusionBudget, AreaOcclusionError, AreaOcclusionTermination, OccludedAreaReference,
};
use fs_render::photometry::area::reference::{
    AreaReferenceBudget, AreaReferenceError, UnoccludedAreaReference,
};
use fs_render::photometry::source::{LUMINOUS_EFFICACY_LM_W, photopic_weight_nm};
use fs_render::photometry::transport::PhotometricMedium;
use fs_render::volumes::transport::enclosure::GridEnclosureBudget;

fn check_component(value: &IlluminanceEstimate) {
    assert_eq!(value.numerical.kind, NumericalKind::Enclosure);
    assert!(value.numerical.lo.is_finite() && value.numerical.hi.is_finite());
    assert!(0.0 <= value.numerical.lo && value.numerical.lo <= value.lux);
    assert!(value.lux <= value.numerical.hi);
    assert_eq!(value.total_error_bound_lux, None);
}

fn check_zero(value: &IlluminanceEstimate) {
    check_component(value);
    assert_eq!(
        (value.lux, value.numerical.lo, value.numerical.hi),
        (0.0, 0.0, 0.0)
    );
}

// Reconstruct only the public native source binding, without another integrator.
fn rebound(source: &PhotometricAreaEmitter, primitive: usize) -> PhotometricAreaEmitter {
    assert!(source.projector().is_none());
    PhotometricAreaEmitter::try_new(
        primitive,
        source.flux_lm(),
        1.0,
        source.spectrum().clone(),
        source.sidedness(),
    )
    .unwrap()
}

#[test]
fn transformed_reordered_fixtures_keep_emitter_and_fixture_finite_attribution() {
    with_cx(|cx, _| {
        let sideways =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [2.0, -1.0, 0.0])
                .unwrap();
        let turn =
            RigidTransform::try_new([0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2], [-2.0, 0.0, 0.5])
                .unwrap();
        let mut dark = finite_fixture(2, turn);
        dark.surfaces[1].source = dark.surfaces[1].source.dimmed(0.0).unwrap();
        let rig = assemble(
            &[
                fixture(90, placement([0.0; 3]), vec![isotropic([0.0; 3], 1e6, 1.0)]),
                finite_fixture(7, sideways),
                fixture(13, placement([0.0; 3]), Vec::new()),
                dark,
            ],
            cx,
        );
        let mut primitives = vec![Primitive {
            shape: Shape::Mesh(square(0.3)),
            material: gray(),
            emission: None,
        }];
        rig.append_geometry_to(&mut primitives, cx).unwrap();
        primitives.swap(1, 4);
        let extra = PhotometricAreaEmitter::try_new(
            0,
            12.0 * PI,
            0.5,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap();
        let reference = UnoccludedAreaReference::new_with_rig(
            &primitives,
            &rig,
            core::slice::from_ref(&extra),
            AreaReferenceBudget::default(),
            cx,
        )
        .unwrap();
        let receiver = probe([4.0, 0.0, 4.0], [-1.0, 0.0, -1.0]);
        let result = reference.illuminance_breakdown(cx, &receiver).unwrap();
        assert_eq!(result.total, reference.illuminance(cx, &receiver).unwrap());
        assert_eq!(result.total.total_error_bound_lux, None);
        assert_eq!(
            reference
                .source_bindings()
                .iter()
                .map(|b| b.primitive_index)
                .collect::<Vec<_>>(),
            [4, 2, 3, 1, 0]
        );
        assert_eq!(result.emitters.len(), 5);
        assert_eq!(
            result
                .fixtures
                .iter()
                .map(|f| f.fixture)
                .collect::<Vec<_>>(),
            [FixtureId(90), FixtureId(7), FixtureId(13), FixtureId(2)]
        );
        for (index, component) in result.emitters.iter().enumerate() {
            assert_eq!(component.source, index);
            assert_eq!(
                component.identity,
                rig.surface_emitter_identities().get(index).copied()
            );
            check_component(&component.estimate);
            let source = rig.surface_sources().get(index).unwrap_or(&extra);
            let single = [rebound(
                source,
                reference.source_bindings()[index].primitive_index,
            )];
            let control_reference = UnoccludedAreaReference::new(
                &primitives,
                &single,
                AreaReferenceBudget::default(),
                cx,
            )
            .unwrap();
            let control = control_reference.illuminance(cx, &receiver).unwrap();
            assert_eq!(component.estimate.numerical, control.numerical);
            assert_eq!(component.estimate.lux.to_bits(), control.lux.to_bits());
            let flat = control_reference
                .illuminance_breakdown(cx, &receiver)
                .unwrap();
            assert_eq!(flat.total, control);
            assert!(flat.fixtures.is_empty());
            assert_eq!(flat.emitters.len(), 1);
            assert_eq!(flat.emitters[0].source, 0);
            assert_eq!(flat.emitters[0].identity, None);
            assert_eq!(flat.emitters[0].estimate, component.estimate);
        }
        assert!(result.emitters[0].estimate.lux > 0.0 && result.emitters[1].estimate.lux > 0.0);
        assert!(result.emitters[2].estimate.lux > 0.0 && result.emitters[4].estimate.lux > 0.0);
        check_zero(&result.emitters[3].estimate);
        check_zero(&result.fixtures[0].estimate); // Bright point light is outside this component.
        check_zero(&result.fixtures[2].estimate); // Empty fixtures still have stable rows.
        for (fixture_index, range) in [(1, 0..2), (3, 2..4)] {
            let estimate = &result.fixtures[fixture_index].estimate;
            check_component(estimate);
            let sum: f64 = result.emitters[range].iter().map(|e| e.estimate.lux).sum();
            assert!(estimate.numerical.lo <= sum && sum <= estimate.numerical.hi);
            close(estimate.lux, sum);
        }
        let assigned: f64 = result.fixtures.iter().map(|f| f.estimate.lux).sum();
        close(assigned + result.emitters[4].estimate.lux, result.total.lux);
        assert!(
            assigned < result.total.lux,
            "supplemental light must remain unassigned"
        );
    });
}

#[test]
fn one_adaptive_shadow_and_spectral_extinction_solve_retains_each_source_interval() {
    with_cx(|cx, _| {
        let rig = assemble(&[finite_fixture(31, placement([0.0; 3]))], cx);
        // At z=1 this strip shadows source-plane x in [-.8,-.7] on the left
        // aperture. The right aperture remains wholly visible.
        let mut primitives = vec![Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [-0.4, -0.5, 1.0],
                    [-0.35, -0.5, 1.0],
                    [-0.35, 0.5, 1.0],
                    [-0.4, 0.5, 1.0],
                ],
                vec![[0, 1, 2], [0, 2, 3]],
            )),
            material: gray(),
            emission: None,
        }];
        rig.append_geometry_to(&mut primitives, cx).unwrap();
        let density = [1.0, 2.0];
        let grid = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 1.0]),
                MassAbsorptionSpectrum::table(&[(380.0, 0.02), (780.0, 0.2)]).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 2,
                    max_steps_per_segment: 32,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.04).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let receiver = probe([0.0, 0.0, 2.0], [0.0, 0.0, -1.0]);
        let fine_budget = AreaOcclusionBudget {
            absolute_error_lux: 0.01,
            max_patch_evaluations: 4096,
            ..AreaOcclusionBudget::default()
        };
        let coarse_budget = AreaOcclusionBudget {
            max_patch_evaluations: 4,
            ..fine_budget
        };
        for spectral in [false, true] {
            let medium = if spectral {
                PhotometricMedium::Scattering(&grid)
            } else {
                PhotometricMedium::Homogeneous(HomogeneousMedium::vacuum())
            };
            let reference = OccludedAreaReference::new_with_rig_and_medium(
                &primitives,
                &rig,
                &[],
                medium.clone(),
                AreaReferenceBudget::default(),
                GridEnclosureBudget::default(),
                cx,
            )
            .unwrap();
            let coarse = reference
                .illuminance_breakdown(cx, &receiver, coarse_budget)
                .unwrap();
            let fine = reference
                .illuminance_breakdown(cx, &receiver, fine_budget)
                .unwrap();
            assert_eq!(
                coarse.total,
                reference.illuminance(cx, &receiver, coarse_budget).unwrap()
            );
            assert_eq!(
                fine.total,
                reference.illuminance(cx, &receiver, fine_budget).unwrap()
            );
            assert_eq!(
                coarse.total.termination,
                AreaOcclusionTermination::PatchBudget
            );
            assert_eq!(coarse.total.evaluated_patches, 4);
            assert!(coarse.total.uncertain_visibility_patches > 0);
            assert!(fine.total.evaluated_patches > 4);
            assert!(fine.total.evaluated_patches <= fine_budget.max_patch_evaluations);
            assert_eq!(fine.total.total_error_bound_lux, None);
            assert_eq!(fine.emitters.len(), 2);
            assert_eq!(fine.fixtures.len(), 1);
            assert_eq!(fine.fixtures[0].fixture, FixtureId(31));
            check_component(&fine.fixtures[0].estimate);
            // Resolved clear leaves leave the adaptive frontier, but their
            // contributions must survive both the initial and refined report.
            assert!(coarse.emitters[1].estimate.numerical.lo > 0.0);
            assert!(fine.emitters[1].estimate.numerical.lo > 0.0);
            let before = coarse.emitters[0].estimate.numerical;
            let after = fine.emitters[0].estimate.numerical;
            assert!(
                after.hi - after.lo < before.hi - before.lo,
                "shared refinement must tighten the shadowed emitter's contribution"
            );
            for (index, component) in fine.emitters.iter().enumerate() {
                check_component(&component.estimate);
                assert_eq!(component.source, index);
                assert_eq!(
                    component.identity,
                    Some(rig.surface_emitter_identities()[index])
                );
                assert!(
                    component.estimate.numerical.lo >= coarse.emitters[index].estimate.numerical.lo
                );
                assert!(
                    component.estimate.numerical.hi <= coarse.emitters[index].estimate.numerical.hi
                );
                let single = [rebound(
                    &rig.surface_sources()[index],
                    reference.source_bindings()[index].primitive_index,
                )];
                // Keep all original emitter bodies and blockers. An independent
                // source solve may spend its budget differently, so require overlap.
                let control = OccludedAreaReference::new_with_medium(
                    &primitives,
                    &single,
                    medium.clone(),
                    AreaReferenceBudget::default(),
                    GridEnclosureBudget::default(),
                    cx,
                )
                .unwrap()
                .illuminance(cx, &receiver, fine_budget)
                .unwrap();
                assert!(component.estimate.numerical.lo <= control.numerical.hi);
                assert!(control.numerical.lo <= component.estimate.numerical.hi);
            }
            if spectral {
                assert_eq!(
                    fine.total.physical_grid_identity,
                    Some(grid.content_identity())
                );
                assert_eq!(fine.total.absorption_per_m, None);
                assert_eq!(fine.total.scattering_per_m, None);
                assert!(fine.total.varying_transmittance_patches > 0);
            } else {
                assert_eq!(fine.total.physical_grid_identity, None);
                assert_eq!(fine.total.varying_transmittance_patches, 0);
            }
            let sum: f64 = fine.emitters.iter().map(|e| e.estimate.lux).sum();
            assert!(
                fine.fixtures[0].estimate.numerical.lo <= sum
                    && sum <= fine.fixtures[0].estimate.numerical.hi
            );
            // A component has no separate tolerance guarantee or in-scattering term.
            assert_eq!(fine.fixtures[0].estimate.total_error_bound_lux, None);
        }
    });
}

#[test]
fn breakdown_caps_late_refusals_and_cancellation_leave_references_reusable() {
    with_cx(|cx, gate| {
        let mut authored = finite_fixture(71, placement([0.0; 3]));
        let elevated = &authored.geometry[1].instance;
        authored.geometry[1].instance = GeometryInstance::try_new(
            elevated.object_id(),
            elevated.geometry_identity(),
            elevated.geometry().clone(),
            placement([0.75, 0.0, 1.0]),
        )
        .unwrap();
        let rig = assemble(
            &[authored, fixture(72, placement([0.0; 3]), Vec::new())],
            cx,
        );
        let mut primitives = Vec::new();
        rig.append_geometry_to(&mut primitives, cx).unwrap();
        let receiver = probe([0.0, 0.0, 3.0], [0.0, 0.0, -1.0]);
        let coplanar = probe([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]);
        let short = AreaReferenceBudget {
            max_outputs: 4,
            ..AreaReferenceBudget::default()
        };
        let enough = AreaReferenceBudget {
            max_outputs: 5,
            ..short
        };
        let budget = AreaOcclusionBudget::default();
        let capped =
            UnoccludedAreaReference::new_with_rig(&primitives, &rig, &[], short, cx).unwrap();
        let capped_shadow =
            OccludedAreaReference::new_with_rig(&primitives, &rig, &[], short, cx).unwrap();
        let original = capped.illuminance(cx, &receiver).unwrap();
        let shadow_original = capped_shadow.illuminance(cx, &receiver, budget).unwrap();
        assert!(matches!(
            capped.illuminance_breakdown(cx, &receiver),
            Err(AreaReferenceError::Budget)
        ));
        assert!(matches!(
            capped_shadow.illuminance_breakdown(cx, &receiver, budget),
            Err(AreaOcclusionError::Budget)
        ));
        // Admission counts the total, both finite emitters, and both fixtures.
        let reference =
            UnoccludedAreaReference::new_with_rig(&primitives, &rig, &[], enough, cx).unwrap();
        let shadow =
            OccludedAreaReference::new_with_rig(&primitives, &rig, &[], enough, cx).unwrap();
        let expected = reference.illuminance_breakdown(cx, &receiver).unwrap();
        let shadow_expected = shadow.illuminance_breakdown(cx, &receiver, budget).unwrap();
        assert_eq!(expected.total, original);
        assert_eq!(shadow_expected.total, shadow_original);
        assert!(matches!(
            reference.illuminance_breakdown(cx, &coplanar),
            Err(AreaReferenceError::UnresolvedSourcePlane { triangle: 2 })
        ));
        assert!(matches!(
            shadow.illuminance_breakdown(cx, &coplanar, budget),
            Err(AreaOcclusionError::Reference(
                AreaReferenceError::UnresolvedSourcePlane { triangle: 2 }
            ))
        ));
        // A later-source refusal publishes no partial report and cannot poison reuse.
        assert_eq!(
            reference.illuminance_breakdown(cx, &receiver).unwrap(),
            expected
        );
        assert_eq!(
            shadow.illuminance_breakdown(cx, &receiver, budget).unwrap(),
            shadow_expected
        );
        gate.request();
        assert!(matches!(
            reference.illuminance_breakdown(cx, &receiver),
            Err(AreaReferenceError::Cancelled)
        ));
        assert!(matches!(
            shadow.illuminance_breakdown(cx, &receiver, budget),
            Err(AreaOcclusionError::Cancelled)
        ));
        with_cx(|retry, _| {
            assert_eq!(
                reference.illuminance_breakdown(retry, &receiver).unwrap(),
                expected
            );
            assert_eq!(
                shadow
                    .illuminance_breakdown(retry, &receiver, budget)
                    .unwrap(),
                shadow_expected
            );
        });
    });
}

// Closed-form integral of h²/r⁴ over a horizontal rectangle, independent of
// the native projected-triangle enclosure. Inputs and output use metres and sr.
fn rectangle_factor(x: f64, y: f64, half: f64, height: f64) -> f64 {
    let primitive = |x: f64, y: f64| {
        let ax = (x * x + height * height).sqrt();
        let by = (y * y + height * height).sqrt();
        0.5 * (x / ax * (y / ax).atan() + y / by * (x / by).atan())
    };
    primitive(x + half, y + half) - primitive(x - half, y + half) - primitive(x + half, y - half)
        + primitive(x - half, y - half)
}

#[test]
fn three_hundred_fixtures_keep_six_hundred_calibrated_finite_components() {
    with_cx(|cx, _| {
        // This is functional coverage of realistic counts, not a timing or
        // measured physical-error assertion. Unoccluded reference scope is explicit.
        let authored: Vec<_> = (0..300)
            .map(|index| {
                finite_fixture(
                    index as u64 + 1,
                    placement([(index % 20) as f64 - 9.5, (index / 20) as f64 - 7.0, 0.0]),
                )
            })
            .collect();
        let rig = assemble(&authored, cx);
        let mut primitives = Vec::new();
        rig.append_geometry_to(&mut primitives, cx).unwrap();
        let reference = UnoccludedAreaReference::new_with_rig(
            &primitives,
            &rig,
            &[],
            AreaReferenceBudget {
                max_outputs: 901,
                ..AreaReferenceBudget::default()
            },
            cx,
        )
        .unwrap();
        let receiver = probe([0.0, 0.0, 10.0], [0.0, 0.0, -1.0]);
        let result = reference.illuminance_breakdown(cx, &receiver).unwrap();
        assert_eq!(result.total, reference.illuminance(cx, &receiver).unwrap());
        assert_eq!(result.emitters.len(), 600);
        assert_eq!(result.fixtures.len(), 300);
        let mut independent_total = 0.0;
        for (index, component) in result.emitters.iter().enumerate() {
            check_component(&component.estimate);
            assert_eq!(component.source, index);
            assert_eq!(
                component.identity,
                Some(FixtureEmitterIdentity {
                    fixture: FixtureId(index as u64 / 2 + 1),
                    emitter: EmitterId(index as u64 % 2 + 11),
                })
            );
            let fixture_index = index / 2;
            let x = (fixture_index % 20) as f64 - 9.5 + (index % 2) as f64 * 1.5 - 0.75;
            let y = (fixture_index / 20) as f64 - 7.0;
            let source = &rig.surface_sources()[index];
            let photopic_mass: f64 = source
                .spectrum()
                .bins_w_per_lumen_nm()
                .iter()
                .enumerate()
                .map(|(bin, value)| value * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W)
                .sum();
            let expected =
                source.flux_lm() / (0.25 * PI) * photopic_mass * rectangle_factor(x, y, 0.25, 10.0);
            close(component.estimate.lux, expected);
            independent_total += expected;
        }
        for (index, fixture) in result.fixtures.iter().enumerate() {
            check_component(&fixture.estimate);
            assert_eq!(fixture.fixture, FixtureId(index as u64 + 1));
            close(
                fixture.estimate.lux,
                result.emitters[2 * index].estimate.lux
                    + result.emitters[2 * index + 1].estimate.lux,
            );
        }
        close(result.total.lux, independent_total);
        assert_eq!(result.total.total_error_bound_lux, None);
    });
}
