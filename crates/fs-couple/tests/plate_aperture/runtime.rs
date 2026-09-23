use super::*;
use fs_couple::bernoulli_aperture::dynamic::{ApertureDrive, ApertureState, ApertureTerminal, DynamicAperture};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, TubeFrame, UniformTubeSpec};
use fs_dcontact::Obstacle;

const DT: f64 = 1e-5;
fn tube_spec() -> UniformTubeSpec {
    UniformTubeSpec {
        length_m: 32.0 * (343.0 * DT), radius_m: 0.007, sound_speed_m_s: 343.0,
        terminal_reflection: -0.8, max_length_error_m: 1e-12, max_wave_memory_bytes: 1 << 20,
    }
}
fn lay() -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic generalized plate-slit contact; not measured cane".into())
        .unwrap().with_internal_loss(5.0).unwrap()
}
fn valve(r: PlateApertureReduction, steps: u64, state: ApertureState) -> DynamicAperture {
    DynamicAperture::from_plate(r, 1.2, tube_spec().characteristic_impedance(1.2).unwrap(),
        DT, steps, state, lay()).unwrap()
}
fn rest(r: &PlateApertureReduction) -> ApertureState {
    ApertureState { opening_m: r.options().rest_opening_m, opening_velocity_m_s: 0.0 }
}
fn tube(r: PlateApertureReduction, steps: u64) -> ApertureTube {
    let state = rest(&r);
    ApertureTube::new(valve(r, steps, state), tube_spec()).unwrap()
}
fn drive(sample: usize, pressure: f64) -> TubeDrive {
    TubeDrive { upstream_pressure_pa: if sample < 256 { pressure } else { 0.0 }, body_flow_m3_s: 0.0 }
}
fn signature(f: TubeFrame) -> [u64; 9] {
    [f.aperture.state.opening_m, f.aperture.state.opening_velocity_m_s,
        f.aperture.bore_pressure_pa, f.aperture.swept_flow_m3_s, f.aperture.jet_flow_m3_s,
        f.waveguide.incoming_pressure_pa, f.waveguide.terminal_pressure_pa,
        f.stored_energy_j, f.dissipated_energy_j].map(f64::to_bits)
}

#[test]
fn bound_plate_uses_the_existing_junction_and_conserves_coupled_pressure_work() {
    let r = reduction(4e9, 900.0);
    let area = r.pressure_area_m2();
    let pressure = r.closing_pressure_pa() * 0.25;
    let spec = r.dynamic_spec(1.2, tube_spec().characteristic_impedance(1.2).unwrap(), DT, 1024);
    let independent = DynamicAperture::new(spec, rest(&r), lay()).unwrap();
    assert!(independent.plate_reduction().is_none());
    let mut reference = ApertureTube::new(independent, tube_spec()).unwrap();
    let mut actual = tube(r, 1024);
    for sample in 0..1024 {
        let old = actual.aperture().state();
        let before = actual.stored_energy_j();
        let f = actual.step(drive(sample, pressure)).unwrap();
        let expected = reference.step(drive(sample, pressure)).unwrap();
        assert_eq!(signature(f), signature(expected), "no second solver or observer gain");
        let vm = f64::midpoint(old.opening_velocity_m_s, f.aperture.state.opening_velocity_m_s);
        near(f.aperture.swept_flow_m3_s, -area * vm, 1e-12);
        let scale = before + f.stored_energy_j + f.dissipated_energy_j + f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs() <= 3e-10 * scale.max(f64::MIN_POSITIVE));
        assert!(f.aperture.flow_residual_m3_s.abs() <= 1e-12);
        let plate = actual.aperture().plate_reduction().unwrap();
        let tip = plate.options().slit_edges[0][0];
        let motion = plate.nodal_motion(tip, f.aperture.state.opening_m,
            f.aperture.state.opening_velocity_m_s).unwrap();
        near(motion.displacement_rotation[0], plate.shape_per_opening()[tip][0]
            * (f.aperture.state.opening_m - plate.options().rest_opening_m), 1e-12);
        // The support is physically stationary, not driven by an amplitude display.
        let root = plate.chart().mesh.nodes.iter().position(|p| p.0 == 0.0).unwrap();
        let fixed = plate.nodal_motion(root, f.aperture.state.opening_m,
            f.aperture.state.opening_velocity_m_s).unwrap();
        assert!(fixed.displacement_rotation.iter().chain(&fixed.velocity_rotation_rate).all(|x| *x == 0.0));
    }
    assert!(actual.stored_energy_j() > 0.0, "release preserves ringdown instead of resetting");
}

#[test]
fn equal_static_closure_but_different_specimen_mass_changes_the_coupled_trajectory() {
    let a = reduction(4e9, 900.0);
    let b = reduction(4e9, 1800.0);
    near(a.closing_pressure_pa(), b.closing_pressure_pa(), 2e-5);
    let pressure = a.closing_pressure_pa() * 0.25;
    let mut light = tube(a, 512);
    let mut heavy = tube(b, 512);
    let (mut opening_change, mut pressure_change, mut swept_change) = (0.0_f64, 0.0_f64, 0.0_f64);
    for sample in 0..512 {
        let (a, b) = (light.step(drive(sample, pressure)).unwrap(), heavy.step(drive(sample, pressure)).unwrap());
        opening_change = opening_change.max((a.aperture.state.opening_m-b.aperture.state.opening_m).abs());
        pressure_change = pressure_change.max((a.aperture.bore_pressure_pa-b.aperture.bore_pressure_pa).abs());
        swept_change = swept_change.max((a.aperture.swept_flow_m3_s-b.aperture.swept_flow_m3_s).abs());
    }
    assert!(opening_change > 1e-9 && pressure_change > 1e-6 && swept_change > 1e-12);
}

#[test]
fn slope_refusal_preserves_both_the_valve_and_tube_history_for_exact_retry() {
    let (chart, mut options) = fixture(4e9, 900.0);
    options.max_slope = 1e-6; // A deliberately tight declared linear domain.
    let plate = PlateApertureReduction::from_chart(chart, options, &CancelGate::new()).unwrap();
    let pressure = plate.closing_pressure_pa();
    let mut reference = tube(plate.clone(), 128);
    let mut candidate = tube(plate, 128);
    let gentle = TubeDrive { upstream_pressure_pa: pressure*1e-7, body_flow_m3_s: 0.0 };
    // Build nonzero traveling-wave history before testing transactional refusal.
    for _ in 0..70 { reference.step(gentle).unwrap(); candidate.step(gentle).unwrap(); }
    let before = (candidate.aperture().state(), candidate.stored_energy_j().to_bits());
    let err = candidate.step(TubeDrive { upstream_pressure_pa: pressure*100.0, body_flow_m3_s: 0.0 }).unwrap_err();
    assert!(matches!(err, fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription {
        what: "plate aperture exceeds its declared linear-slope domain"
    }), "expected the specimen gate, not a solver failure: {err:?}");
    assert_eq!((candidate.aperture().state(), candidate.stored_energy_j().to_bits()), before);
    assert_eq!(candidate.aperture().accepted_steps(), 70);
    for _ in 0..32 { assert_eq!(signature(candidate.step(gentle).unwrap()), signature(reference.step(gentle).unwrap())); }
}

#[test]
fn cancellation_and_budget_extension_keep_specimen_state_and_all_propagation_memory() {
    let r = reduction(4e9, 900.0);
    let inputs: Vec<_> = (0..513).map(|n| drive(n, 0.25*r.closing_pressure_pa())).collect();
    let mut baseline = tube(r.clone(), 513);
    let expected: Vec<_> = inputs.iter().map(|&d| signature(baseline.step(d).unwrap())).collect();
    let mut resumed = tube(r, 37);
    let sentinel = TubeFrame { stored_energy_j: -1.0, ..TubeFrame::default() };
    let mut frames = vec![sentinel; 513];
    let active = CancelGate::new_clock_free();
    let progress = resumed.advance_block(&inputs, &mut frames, &active).unwrap();
    assert_eq!(progress.completed, 37);
    assert_eq!(progress.terminal, ApertureTerminal::BudgetExhausted);
    let state = resumed.aperture().state();
    resumed.extend_step_budget(513).unwrap();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    let progress = resumed.advance_block(&inputs[37..], &mut frames[37..], &cancelled).unwrap();
    assert_eq!(progress.completed, 0);
    assert_eq!(progress.terminal, ApertureTerminal::Cancelled);
    assert_eq!(resumed.aperture().state(), state);
    assert!(frames[37..].iter().all(|f| *f == sentinel));
    for (data, out) in inputs[37..].chunks(17).zip(frames[37..].chunks_mut(17)) {
        resumed.advance_block(data, out, &active).unwrap();
    }
    assert_eq!(frames.into_iter().map(signature).collect::<Vec<_>>(), expected);
    assert!(resumed.aperture().plate_reduction().is_some());
}

#[test]
fn specimen_initial_state_and_motion_queries_cannot_bypass_linear_geometry_limits() {
    let plate = reduction(4e9, 900.0);
    assert!(DynamicAperture::from_plate(plate.clone(), 1.2, 1e6, DT, 128,
        ApertureState { opening_m: 1.0, opening_velocity_m_s: 0.0 }, lay()).is_err());
    assert!(plate.nodal_motion(usize::MAX, plate.options().rest_opening_m, 0.0).is_err());
    assert!(plate.nodal_motion(0, plate.options().rest_opening_m, f64::NAN).is_err());
    let mut valve = valve(plate.clone(), 4, rest(&plate));
    // Zero pressure does not prescribe displacement after release.
    let start = valve.state();
    assert!(valve.step(ApertureDrive { upstream_pressure_pa: f64::NAN, ..ApertureDrive::default() }).is_err());
    assert_eq!(valve.state(), start);
    assert_eq!(valve.accepted_steps(), 0);
}

#[test]
fn resolved_regional_material_receipts_survive_nonlinear_tube_composition() {
    use fs_couple::thin_plate::{compile_plate_material_chart, PlateMaterialModel,
        PlateRegionMaterial, PlateThicknessConstraint};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::state_point::{resolve_material_state_point, MaterialPropertySelection,
        ScalarAdmissibility, ScalarPropertyRequirement};
    use fs_qty::{Density, Dims, Pressure, QuantitySpec};
    let properties = [("density", Density::DIMS, 900.0),
        ("young_modulus", Pressure::DIMS, 4e9), ("poisson_ratio", Dims::NONE, 0.3)];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 293.15, 293.15),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            provenance: Provenance { source: "synthetic plate-slit material binding test".into(),
                license: "CC0-1.0".into(), artifact: None }, observations: vec![],
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-plate".into(), phase: "solid".into(),
        process: "authored homogeneous specimen".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let requirements: Vec<_> = properties.iter().map(|(name, dims, _)|
        ScalarPropertyRequirement::try_with_quantity(*name, QuantitySpec::dimensional(*dims),
            ScalarAdmissibility::Finite).unwrap()).collect();
    let material = resolve_material_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap(), &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let (chart, options) = fixture(4e9, 900.0);
    let count = chart.mesh.tris.len();
    let bound = compile_plate_material_chart(chart.mesh, chart.boundary_nodes, vec![PlateRegionMaterial {
        region: fs_plate::PlateRegion { name: "vamp".into(), triangle_indices: (0..count).collect() },
        material: material.clone(), model: PlateMaterialModel::Isotropic,
        thickness_constraint: PlateThicknessConstraint::FixedThickness(0.0003), material_angle_rad: 0.0,
    }]).unwrap();
    let identity = bound.specimen_identity();
    let r = PlateApertureReduction::from_material_chart(bound, options, &CancelGate::new()).unwrap();
    let pressure = 0.25*r.closing_pressure_pa();
    let mut model = tube(r, 128);
    for n in 0..128 { model.step(drive(n, pressure)).unwrap(); }
    let retained = model.aperture().plate_reduction().unwrap().material_chart().unwrap();
    assert_eq!(retained.specimen_identity(), identity);
    assert_eq!(retained.regions()[0].input.material, material);
    assert_eq!(retained.regions()[0].input.region.name, "vamp");
}
