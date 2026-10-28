use super::*;
use fs_couple::bernoulli_aperture::modal::material::ModalPlateOptions;
use fs_couple::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA;
use fs_couple::thin_plate::{
    PlateChartRadiation, PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint,
    ResolvedPlateChart, compile_plate_material_chart,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

fn source(rho: f64, eta: f64) -> ResolvedMaterialStatePoint {
    source_in_band(rho, eta, 100000.)
}
fn source_in_band(rho: f64, eta: f64, upper: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", rho, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        (ETA, eta, DynViscosity::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290., 300.)
                    .with("omega", 0., upper),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic multimode material specimen".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 293.15)
        .unwrap()
        .with("omega", 100.)
        .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn plate(rho: f64, eta: f64) -> (ResolvedPlateChart, ModalPlateOptions) {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let supports = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p.0 == 0.)
        .map(|(i, _)| i)
        .collect();
    let edges: Vec<_> = mesh
        .boundary_edges()
        .into_iter()
        .filter(|&(a, b)| {
            (mesh.nodes[a].0 - 0.025).abs() < 1e-12 && (mesh.nodes[b].0 - 0.025).abs() < 1e-12
        })
        .map(|(a, b)| [a, b])
        .collect();
    let mut force = vec![0.; mesh.nodes.len()];
    let node = edges[0][0];
    force[node] = 1.;
    let elements = (0..mesh.tris.len()).collect();
    let chart = compile_plate_material_chart(
        mesh,
        supports,
        vec![PlateRegionMaterial {
            region: PlateRegion {
                name: "blade".into(),
                triangle_indices: elements,
            },
            material: source(rho, eta),
            model: PlateMaterialModel::Isotropic,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.0003),
            material_angle_rad: 0.,
        }],
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    (
        chart,
        ModalPlateOptions {
            reduction: PlateChartRadiation {
                assembly: AssemblyOptions {
                    pretension: 0.,
                    support: EdgeSupport::Clamped,
                },
                eigenvalue_window: (1., (std::f64::consts::TAU * 3000.).powi(2)),
                n_modes: 3,
                damping_ratio: 0.,
                unit_force_weights: force,
            },
            slit_edges: edges,
            maximum_trace_slope: 0.1,
            maximum_slit_deviation_m: 2e-5,
            maximum_angular_step: 0.25,
            max_nodes: 100,
            max_triangles: 100,
            max_projection_terms: 100000,
        },
    )
}
fn bound(plates: Vec<(ResolvedPlateChart, ModalPlateOptions)>, z: f64) -> ModalAperture {
    let n = plates.iter().map(|(_, o)| o.reduction.n_modes).sum();
    ModalAperture::from_material_plates(
        plates,
        4e-4,
        super::network::gas(),
        z,
        1. / f64::from(MECHANICS_RATE),
        4096 * CLOCK_REFINEMENT as u64,
        limits(),
        vec![BladeState::default(); n],
        contact(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-24),
        "{a:e} != {b:e}"
    );
}

#[test]
fn material_basis_preserves_original_mass_loss_and_signed_work_ports() {
    let (chart, options) = plate(900., 2e4);
    let original = chart.specimen_identity();
    let (bodies, _) = chart
        .certified_radiators_with_force_ports(&options.reduction, &[], 0)
        .unwrap();
    let mut a = bound(vec![(chart, options)], 1e6);
    let retained = &a.material_plates()[0];
    assert_eq!(retained.source().specimen_identity(), original);
    assert!(retained.source().bending_viscosity().is_some());
    assert_eq!(retained.mode_range(), 0..3);
    let mesh = &retained.source().chart().mesh;
    for (i, (m, b)) in a.spec().modes.iter().zip(&bodies).enumerate() {
        near(m.stiffness_n_m / m.mass_kg, b.omega * b.omega, 2e-5);
        near(m.damping_n_s_m / m.mass_kg, 2. * b.zeta * b.omega, 2e-5);
        let mut area = 0.;
        for tri in &mesh.tris {
            let [x, y, z] = tri.map(|j| mesh.nodes[j]);
            let da = 0.5 * ((y.0 - x.0) * (z.1 - x.1) - (y.1 - x.1) * (z.0 - x.0));
            area += da
                * tri
                    .iter()
                    .map(|&j| retained.displacement_shapes()[j][i])
                    .sum::<f64>()
                / 3.;
        }
        assert!((m.pressure_area_m2 - area).abs() < 2.5e-16);
        let force: f64 = retained
            .options()
            .reduction
            .unit_force_weights
            .iter()
            .zip(retained.displacement_shapes())
            .map(|(f, p)| f * p[i])
            .sum();
        assert!((m.force_projection[0] - force).abs() < 1e-10);
    }
    for _ in 0..32 {
        a.step(drive(5.), [0.; 2]).unwrap();
    }
    let root = a.material_plates()[0].source().chart().boundary_nodes[0];
    assert_eq!(a.nodal_motion(0, root).unwrap(), BladeState::default());
    assert!(a.nodal_motion(1, root).is_err());
    assert!(a.nodal_motion(0, usize::MAX).is_err());
}

#[test]
fn density_and_loss_data_change_independent_material_modes_without_retuning() {
    let a = bound(vec![plate(900., 2e4)], 1e6);
    let heavy = bound(vec![plate(1800., 2e4)], 1e6);
    let lossy = bound(vec![plate(900., 4e4)], 1e6);
    for ((a, b), c) in a
        .spec()
        .modes
        .iter()
        .zip(&heavy.spec().modes)
        .zip(&lossy.spec().modes)
    {
        near(
            (b.stiffness_n_m / b.mass_kg) / (a.stiffness_n_m / a.mass_kg),
            0.5,
            2e-4,
        );
        near(
            c.stiffness_n_m / c.mass_kg,
            a.stiffness_n_m / a.mass_kg,
            2e-4,
        );
        near(
            c.damping_n_s_m / c.mass_kg,
            2. * a.damping_n_s_m / a.mass_kg,
            2e-4,
        );
    }
    assert_ne!(
        a.material_plates()[0].source().specimen_identity(),
        heavy.material_plates()[0].source().specimen_identity()
    );
}

#[test]
fn material_pair_changes_actual_pressure_and_retains_sources_across_callbacks() {
    use super::network::{gas, sides};
    use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
    use fs_couple::bernoulli_aperture::two_sided::{
        LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, combined_inlet_state,
    };
    use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
    use fs_couple::pcm_wav::baffled::openings::duct::{
        DuctOpening, DuctRadiationSpec, DuctReconstruction,
    };
    use fs_couple::pcm_wav::observation::PressureRenderer;
    let render = |rho| {
        let ducts = sides();
        let (_, z) = combined_inlet_state(&ducts).unwrap();
        let valve = bound(vec![plate(900., 2e4), plate(rho, 2e4)], z);
        let system = TwoSidedAperture::new(
            valve.into(),
            ducts,
            TwoSidedLimits {
                maximum_energy_j: 1.,
                energy_absolute_tolerance_j: 1e-16,
                energy_relative_tolerance: 1e-7,
                flow_absolute_tolerance_m3_s: 1e-16,
                flow_relative_tolerance: 1e-8,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        TwoSidedPressure::with_radiation(
            system,
            DuctRadiationSpec {
                openings: vec![DuctOpening {
                    network: 1,
                    node: 1,
                    center_m: [0.; 2],
                }],
                exterior: None,
                receiver: CircularOutletReceiver {
                    position_m: [0., 0., 0.25],
                    radial_rings: 4,
                    angular_points: 16,
                    maximum_frequency_hz: 1000.,
                },
                reconstruction: DuctReconstruction::Differentiating {
                    mechanics_rate_hz: MECHANICS_RATE,
                    output_rate_hz: 96000,
                },
            },
            64,
        )
        .unwrap()
    };
    let mut reference = vec![0.; 512];
    for (rho, size) in [(900., 1), (900., 37), (1800., 37)] {
        let mut p = render(rho);
        p.set_drive(TwoSidedDrive {
            source_pressure_pa: 5.,
            blade_forces_n: [0.; 2],
        })
        .unwrap();
        let mut out = vec![0.; 512];
        for chunk in out.chunks_mut(size) {
            p.block(chunk).unwrap();
        }
        if size == 1 {
            reference = out;
        } else if rho == 900. {
            assert_eq!(out, reference);
        } else {
            assert!(
                out.iter()
                    .zip(&reference)
                    .any(|(a, b)| (a - b).abs() > 1e-12)
            );
        }
        let LoadedAperture::Modal(a) = p.system().valve() else {
            panic!("material modal valve")
        };
        assert_eq!(a.material_plates().len(), 2);
        assert_eq!(a.modal_states().len(), 6);
        assert_eq!(a.gas(), Some(&gas()));
    }
}

#[test]
fn material_plate_admission_and_summed_motion_fail_without_publishing_state() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(
        ModalAperture::from_material_plates(
            vec![plate(900., 2e4)],
            4e-4,
            super::network::gas(),
            1e6,
            1. / f64::from(MECHANICS_RATE),
            4096,
            limits(),
            vec![BladeState::default(); 3],
            contact(),
            &gate
        )
        .is_err()
    );
    let mut input = plate(900., 2e4);
    input.1.max_projection_terms = 1;
    assert!(
        ModalAperture::from_material_plates(
            vec![input],
            4e-4,
            super::network::gas(),
            1e6,
            1. / f64::from(MECHANICS_RATE),
            4096,
            limits(),
            vec![BladeState::default(); 3],
            contact(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut input = plate(900., 2e4);
    input.1.maximum_trace_slope = 1e-30;
    let mut a = bound(vec![input], 1e6);
    let before = a.modal_states().to_vec();
    assert!(a.step(drive(5.), [0.; 2]).is_err());
    assert_eq!(a.modal_states(), before);
    assert_eq!(a.accepted_steps(), 0);
}

fn one_mode_options(o: &ModalPlateOptions)
    -> fs_couple::bernoulli_aperture::plate::PlateApertureOptions
{
    fs_couple::bernoulli_aperture::plate::PlateApertureOptions {
        assembly: o.reduction.assembly.clone(),
        eigenvalue_window: o.reduction.eigenvalue_window,
        mode_index: 0, eigensolver: fs_plate::SliceOptions::default(),
        max_nodes: o.max_nodes, max_triangles: o.max_triangles,
        slit_edges: o.slit_edges.clone(), rest_opening_m: 4e-4,
        damping_ratio: o.reduction.damping_ratio,
        max_slit_mode_variation: 0.25, max_slope: o.maximum_trace_slope,
    }
}

#[test]
fn material_loss_reaches_single_symmetric_and_independent_reed_motion() {
    use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
    use fs_couple::bernoulli_aperture::plate::PlateApertureReduction;
    let gate = CancelGate::new_clock_free();
    let (chart, options) = plate(900., 2e4);
    let changed = chart.with_region_material("blade", &source(900., 8e4)).unwrap();
    let a = PlateApertureReduction::from_material_chart(chart, one_mode_options(&options), &gate).unwrap();
    let b = PlateApertureReduction::from_material_chart(changed, one_mode_options(&options), &gate).unwrap();
    assert_eq!(a.mass_kg(), b.mass_kg());
    assert_eq!(a.stiffness_n_m(), b.stiffness_n_m());
    assert_eq!(a.pressure_area_m2(), b.pressure_area_m2());
    assert_eq!(a.options().damping_ratio, 0.0);
    let omega = (a.stiffness_n_m() / a.mass_kg()).sqrt();
    near(a.damping_ratio(), 0.5 * 2e4 / 4e9 * omega, 1e-7);
    near(b.damping_ratio(), 4.0 * a.damping_ratio(), 1e-12);
    let symmetric = a.clone().into_symmetric_opposed_pair().unwrap();
    assert_eq!(symmetric.damping_ratio(), a.damping_ratio());
    let gas = super::network::gas();
    let dt = 1.0 / f64::from(MECHANICS_RATE);
    let initial = ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 };
    let scalar = |p| DynamicAperture::from_plate(p, gas.density, 1e6, dt, 2048,
        initial, contact()).unwrap();
    let mut single = scalar(a.clone());
    let mut single_changed = scalar(b.clone());
    let mut symmetric = scalar(symmetric);
    let pair = |blades| PairedAperture::from_plates(blades, gas, 1e6, dt, 2048,
        limits(), [BladeState::default(); 2], contact()).unwrap();
    let mut identical = pair([a.clone(), a.clone()]);
    let mut unequal = pair([a, b]);
    near(unequal.spec().blades[1].damping_n_s_m,
        4.0 * unequal.spec().blades[0].damping_n_s_m, 1e-12);
    let (mut single_pressure_change, mut pair_pressure_change) = (0.0_f64, 0.0_f64);
    for n in 0..2048 {
        let input = drive(if n < 1024 { 5.0 } else { 0.0 });
        let low = single.step(input).unwrap();
        let high = single_changed.step(input).unwrap();
        let equivalent = symmetric.step(input).unwrap();
        let two = identical.step(input).unwrap();
        let different = unequal.step(input).unwrap();
        single_pressure_change = single_pressure_change.max((low.bore_pressure_pa - high.bore_pressure_pa).abs());
        pair_pressure_change = pair_pressure_change.max((two.aperture.bore_pressure_pa - different.aperture.bore_pressure_pa).abs());
        assert!((two.aperture.state.opening_m - equivalent.state.opening_m).abs() < 1e-11);
        assert!((two.aperture.bore_pressure_pa - equivalent.bore_pressure_pa).abs() < 1e-6);
        assert_eq!(two.blades[0], two.blades[1]);
        for frame in [low, high, equivalent, two.aperture, different.aperture] {
            assert!(frame.dissipated_energy_j >= 0.0);
            assert!(frame.balance_residual_j().abs() < 1e-13);
        }
    }
    assert!(single_pressure_change > 1e-6 && pair_pressure_change > 1e-6);
    assert_ne!(unequal.blade_states()[0], unequal.blade_states()[1]);
}

#[test]
fn selected_reed_material_loss_refuses_authored_damping_and_out_of_band_modes() {
    use fs_couple::bernoulli_aperture::plate::PlateApertureReduction;
    let gate = CancelGate::new_clock_free();
    // Explicit zero viscosity still owns its source band and excludes a second
    // authored material-loss law, exactly as in the other plate reductions.
    let (chart, options) = plate(900., 0.0);
    let mut duplicate = one_mode_options(&options);
    duplicate.damping_ratio = 0.02;
    assert!(PlateApertureReduction::from_material_chart(chart.clone(), duplicate, &gate).is_err());
    let narrow = chart.with_region_material("blade", &source_in_band(900., 0.0, 100.)).unwrap();
    assert!(PlateApertureReduction::from_material_chart(narrow, one_mode_options(&options), &gate).is_err());
}

#[path = "material/spatial.rs"]
mod spatial;
#[path = "material/pressure.rs"]
mod pressure;
#[path = "material/wetted.rs"]
mod wetted;
#[path = "material/tensor.rs"]
mod tensor;
