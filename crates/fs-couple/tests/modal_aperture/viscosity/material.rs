//! Actual material tensor -> DKT modes -> moving aperture -> physical pressure.
//! Synthetic data checks coupling; it does not identify an instrument material.
use super::super::network::{gas, sides};
use super::super::*;
use fs_couple::bernoulli_aperture::modal::material::ModalPlateOptions;
use fs_couple::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
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
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as PROPS;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

type Input = (
    ResolvedPlateChart,
    ModalPlateOptions,
    PlateViscosity,
    PlateViscosityLimits,
);

mod keyed;

#[test]
fn mutual_exterior_retains_original_material_modes_viscosity_and_source_receipts() {
    use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeSection};
    use fs_couple::bernoulli_aperture::paired::BladePropagator;
    use fs_couple::bernoulli_aperture::radiation::environment::multiport::{
        MimoRadiationLimits, MimoRadiationModel,
    };
    use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
    use fs_couple::bernoulli_aperture::two_sided::{
        LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, ValveSide,
        combined_inlet_state,
    };
    use fs_couple::pcm_wav::baffled::openings::{CircularOpening, duct::{
        DuctOpening, DuctRadiationSpec, DuctReconstruction,
    }};
    use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
    use fs_math::c64::C64;
    use fs_vfit::{mimo::MimoModel, model::PoleTerm};
    let dt = 1.0 / 48_000.0;
    let gate = CancelGate::new_clock_free();
    let mut ducts = sides();
    ducts[1].graph.nodes[1] = NetworkNode::Junction;
    for (node, delay, radius) in [(2, 12.0, 0.006), (3, 20.0, 0.005)] {
        ducts[1].graph.nodes.push(NetworkNode::Termination { reflection: -0.8 });
        ducts[1].graph.sections.push(TubeSection { nodes: [1, node],
            length_m: delay * gas().sound_speed * dt, radius_m: radius, max_length_error_m: 1e-12 });
        ducts[1].gases.push(gas());
    }
    let (g, z) = combined_inlet_state(&ducts).unwrap();
    let mut inputs = vec![input(1.0, 900.0), input(0.0, 1200.0)];
    for plate in &mut inputs { plate.1.maximum_angular_step = 2.0; }
    let identities: Vec<_> = inputs.iter().map(|p| p.0.specimen_identity()).collect();
    let valve = ModalAperture::from_viscous_material_plates_with_propagator(
        inputs, 4e-4, g, z, dt, 256, limits(), vec![BladeState::default(); 6],
        contact(), BladePropagator::ExactZoh, &gate,
    ).unwrap();
    let original_spec = valve.spec().clone();
    let original_viscosity = valve.viscosity_matrix_n_s_m().unwrap().to_vec();
    let mut radiation = DuctRadiationSpec {
        openings: vec![DuctOpening { network: 0, node: 2, center_m: [0.0, 0.0] },
            DuctOpening { network: 0, node: 3, center_m: [0.05, 0.0] }],
        exterior: Some(g),
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 0.25], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        reconstruction: DuctReconstruction::FiniteDifference { sample_rate_hz: 48_000 },
    };
    let d = vec![1e5, 2e4, 2e4, 1.5e5];
    let source = MimoRadiationModel {
        declaration_identity: [3; 32], realization_identity: [4; 32], exterior: g,
        apertures: radiation.openings.iter().zip([0.006, 0.005]).map(|(o, radius_m)|
            CircularOpening { center_m: o.center_m, radius_m }).collect(),
        band_hz: [1.0, 24_000.0], impedance: MimoModel { ports: 2,
            poles: vec![PoleTerm::Real { pole: -3000.0, residue: 0.0 }],
            residues: vec![d.iter().map(|v| C64::from_re(3000.0 * v)).collect()],
            d, e: vec![0.0; 4], reciprocal: true },
    };
    let system = TwoSidedAperture::new(valve.into(), ducts, TwoSidedLimits {
        maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16,
        energy_relative_tolerance: 1e-7, flow_absolute_tolerance_m3_s: 1e-16,
        flow_relative_tolerance: 1e-8,
    }, &gate).unwrap().with_mimo_environment(ValveSide::Downstream, source,
        radiation.clone(), MimoRadiationLimits { max_states: 16, max_memory_bytes: 1 << 20,
            sampled_passivity_points: 257 }, &gate).unwrap();
    for opening in &mut radiation.openings { opening.network = 1; }
    let mut pressure = TwoSidedPressure::with_radiation(system, radiation, 64).unwrap();
    pressure.set_drive(TwoSidedDrive { source_pressure_pa: 5.0,
        blade_forces_n: [0.0; 2] }).unwrap();
    let mut out = [0.0; 256];
    let mut exchanged = 0.0;
    for block in out.chunks_mut(37) {
        pressure.block(block).unwrap();
        let frame = pressure.last_frame().unwrap();
        assert!(frame.balance_residual_j().abs() < 1e-12);
        exchanged += frame.coupled_terminal_work_j.abs();
    }
    assert!(exchanged > 0.0 && out.iter().any(|p| p.abs() > 1e-12));
    let LoadedAperture::Modal(valve) = pressure.system().valve() else {
        panic!("original material-modal owner");
    };
    assert_eq!(valve.spec(), &original_spec);
    assert_eq!(valve.viscosity_matrix_n_s_m().unwrap(), original_viscosity);
    assert_eq!(valve.accepted_steps(), 256);
    assert_eq!(valve.blade_propagator(), BladePropagator::ExactZoh);
    for (plate, identity) in valve.material_plates().iter().zip(identities) {
        assert_eq!(plate.source().specimen_identity(), identity);
        assert_eq!(plate.viscosity().unwrap().chart_identity(), identity);
    }
}

#[test]
fn exact_material_tensor_reaches_both_ducts_and_microphone_at_the_declared_audio_clock() {
    use fs_couple::bernoulli_aperture::paired::BladePropagator;
    use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
    use fs_couple::bernoulli_aperture::two_sided::{
        LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, combined_inlet_state,
    };
    use fs_couple::pcm_wav::baffled::openings::duct::{
        DuctOpening, DuctRadiationSpec, DuctReconstruction,
    };
    use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
    let render = |loss: f64, block_size| {
        let ducts = sides();
        let (g, z) = combined_inlet_state(&ducts).unwrap();
        let mut inputs = vec![input(loss, 900.0), input(0.0, 1200.0)];
        for p in &mut inputs {
            p.1.maximum_angular_step = 2.0;
        }
        let identities: Vec<_> = inputs.iter().map(|p| p.0.specimen_identity()).collect();
        let valve = ModalAperture::from_viscous_material_plates_with_propagator(
            inputs,
            4e-4,
            g,
            z,
            1.0 / 48_000.0,
            512,
            limits(),
            vec![BladeState::default(); 6],
            contact(),
            BladePropagator::ExactZoh,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        assert_eq!(valve.blade_propagator(), BladePropagator::ExactZoh);
        for (plate, identity) in valve.material_plates().iter().zip(identities) {
            assert_eq!(plate.source().specimen_identity(), identity);
            assert_eq!(plate.viscosity().unwrap().chart_identity(), identity);
        }
        assert!(valve.spec().modes.iter().all(|m| m.damping_n_s_m == 0.0));
        let c = valve.viscosity_matrix_n_s_m().unwrap();
        assert!((0..3).any(|i| (0..3).any(|j| i != j && c[6 * i + j] != 0.0)));
        assert!((3..6).all(|i| (0..6).all(|j| c[6 * i + j] == 0.0 && c[6 * j + i] == 0.0)));
        // Exact selection admits these unchanged material modes at 48 kHz;
        // reconstructing the same mechanics through midpoint still refuses.
        assert!(
            ModalAperture::new(
                valve.spec().clone(),
                vec![BladeState::default(); 6],
                contact()
            )
            .is_err()
        );
        let system = TwoSidedAperture::new(
            valve.into(),
            ducts,
            TwoSidedLimits {
                maximum_energy_j: 1.0,
                energy_absolute_tolerance_j: 1e-16,
                energy_relative_tolerance: 1e-7,
                flow_absolute_tolerance_m3_s: 1e-16,
                flow_relative_tolerance: 1e-8,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let mut pressure = TwoSidedPressure::with_radiation(
            system,
            DuctRadiationSpec {
                openings: vec![DuctOpening {
                    network: 1,
                    node: 1,
                    center_m: [0.0; 2],
                }],
                exterior: None,
                receiver: CircularOutletReceiver {
                    position_m: [0.0, 0.0, 0.25],
                    radial_rings: 4,
                    angular_points: 16,
                    maximum_frequency_hz: 1000.0,
                },
                reconstruction: DuctReconstruction::FiniteDifference {
                    sample_rate_hz: 48_000,
                },
            },
            64,
        )
        .unwrap();
        pressure
            .set_drive(TwoSidedDrive {
                source_pressure_pa: 5.0,
                blade_forces_n: [0.0; 2],
            })
            .unwrap();
        let mut out = vec![0.0; 512];
        for block in out[..256].chunks_mut(block_size) {
            pressure.block(block).unwrap();
        }
        let first = pressure.last_frame().unwrap();
        assert!(first.upstream.incoming_pressure_pa != 0.0);
        assert!(first.downstream.incoming_pressure_pa != 0.0);
        pressure.set_drive(TwoSidedDrive::default()).unwrap();
        for block in out[256..].chunks_mut(block_size) {
            pressure.block(block).unwrap();
        }
        assert!(out.iter().all(|p| p.is_finite()));
        assert!(out.iter().any(|p| p.abs() > 1e-12));
        assert!(pressure.last_frame().unwrap().balance_residual_j().abs() < 1e-12);
        let LoadedAperture::Modal(valve) = pressure.system().valve() else {
            panic!("original modal source");
        };
        assert_eq!(valve.accepted_steps(), 512);
        assert_eq!(valve.material_plates().len(), 2);
        assert!(valve.stored_energy_j() > 0.0);
        out
    };
    let baseline = render(1.0, 1);
    assert_eq!(baseline, render(1.0, 37));
    let different_material = render(4.0, 37);
    assert!(
        baseline
            .iter()
            .zip(different_material)
            .any(|(a, b)| (a - b).abs() > 1e-12)
    );
}

fn source(rho: f64, scale: f64) -> ResolvedMaterialStatePoint {
    let mut values = vec![
        ("density", rho, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
    ];
    for (name, v) in PROPS.into_iter().zip([2e4, 4e4, 1.5e4, 3e3, 2e3, -1e3]) {
        values.push((name, scale * v, DynViscosity::DIMS));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290., 300.)
                    .with("omega", 0., 100000.),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic full-viscosity blade".into(),
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
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 100.)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn input(scale: f64, rho: f64) -> Input {
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
    force[edges[0][0]] = 1.;
    let mut regions = vec![];
    for part in 0..2 {
        let triangles = mesh
            .tris
            .iter()
            .enumerate()
            .filter(|(_, t)| {
                let x = t.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.;
                (x < 0.0125) == (part == 0)
            })
            .map(|(i, _)| i)
            .collect();
        regions.push(PlateRegionMaterial {
            region: PlateRegion {
                name: format!("part-{part}"),
                triangle_indices: triangles,
            },
            material: source(rho, scale * if part == 0 { 1.0 } else { 4.0 }),
            model: PlateMaterialModel::Isotropic,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.0003),
            material_angle_rad: 0.35 * part as f64,
        });
    }
    let chart = compile_plate_material_chart(mesh, supports, regions).unwrap();
    let options = ModalPlateOptions {
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
    };
    let law = PlateViscosity::from_chart(&chart, PROPS, &CancelGate::new_clock_free()).unwrap();
    (
        chart,
        options,
        law,
        PlateViscosityLimits {
            max_modes: 3,
            max_projection_terms: 100000,
            reference_resistance_n_s_m: 1.0,
        },
    )
}
fn make(
    plates: Vec<Input>,
    z: f64,
) -> Result<ModalAperture, fs_couple::acoustic_realize::AcousticRealizeError> {
    let n = plates.iter().map(|p| p.1.reduction.n_modes).sum();
    ModalAperture::from_viscous_material_plates(
        plates,
        4e-4,
        gas(),
        z,
        1. / f64::from(MECHANICS_RATE),
        4096 * CLOCK_REFINEMENT as u64,
        limits(),
        vec![BladeState::default(); n],
        contact(),
        &CancelGate::new_clock_free(),
    )
}
#[test]
fn material_full_matrix_uses_same_eigenbasis_and_work_congruence_as_original_owner() {
    let original = input(1., 900.);
    let identity = original.0.specimen_identity();
    let (image, _) = original
        .0
        .certified_viscous_radiators_with_force_ports(
            &original.1.reduction,
            &[],
            0,
            &original.2,
            original.3,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let a = make(vec![original], 1e6).unwrap();
    let c = a.viscosity_matrix_n_s_m().unwrap();
    let plate = &a.material_plates()[0];
    assert_eq!(plate.source().specimen_identity(), identity);
    assert_eq!(plate.viscosity().unwrap().chart_identity(), identity);
    assert_eq!(plate.viscosity_matrix_n_s_m().unwrap(), c);
    assert!(a.spec().modes.iter().all(|m| m.damping_n_s_m == 0.));
    let mut cross = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let expected = image.damping().matrix_per_s()[3 * i + j]
                * a.spec().modes[i].mass_kg.sqrt()
                * a.spec().modes[j].mass_kg.sqrt();
            let scale = (c[3 * i + i] * c[3 * j + j]).sqrt().max(1e-30);
            assert!((c[3 * i + j] - expected).abs() < 1e-8 * scale);
            if i != j {
                cross = cross.max(c[3 * i + j].abs() / scale);
            }
        }
    }
    assert!(
        cross > 1e-4,
        "heterogeneous anisotropic loss must not become diagonal"
    );
    let twice = make(vec![input(2., 900.)], 1e6).unwrap();
    for (x, y) in c.iter().zip(twice.viscosity_matrix_n_s_m().unwrap()) {
        assert!((2. * x - y).abs() < 1e-10 * c.iter().map(|v| v.abs()).fold(0.0, f64::max));
    }
    for (x, y) in a.spec().modes.iter().zip(&twice.spec().modes) {
        assert_eq!(x.mass_kg, y.mass_kg);
        assert_eq!(x.stiffness_n_m, y.stiffness_n_m);
    }
}
#[test]
fn tensor_material_loss_changes_actual_multirate_pressure_without_resetting_modes() {
    use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
    use fs_couple::bernoulli_aperture::two_sided::{
        LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, combined_inlet_state,
    };
    use fs_couple::pcm_wav::baffled::openings::duct::{
        DuctOpening, DuctRadiationSpec, DuctReconstruction,
    };
    use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
    let render = |loss, size| {
        let ducts = sides();
        let (_, z) = combined_inlet_state(&ducts).unwrap();
        let valve = make(vec![input(loss, 900.), input(0., 1200.)], z).unwrap();
        let c = valve.viscosity_matrix_n_s_m().unwrap();
        assert!((3..6).all(|i| (0..6).all(|j| c[6 * i + j] == 0. && c[6 * j + i] == 0.)));
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
        let mut p = TwoSidedPressure::with_radiation(
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
        .unwrap();
        p.set_drive(TwoSidedDrive {
            source_pressure_pa: 5.,
            blade_forces_n: [0.; 2],
        })
        .unwrap();
        let mut out = vec![0.; 1024];
        for block in out[..512].chunks_mut(size) {
            p.block(block).unwrap();
        }
        p.set_drive(TwoSidedDrive::default()).unwrap();
        for block in out[512..].chunks_mut(size) {
            p.block(block).unwrap();
        }
        assert!(out.iter().all(|p| p.is_finite()));
        assert!(out.iter().any(|p| p.abs() > 1e-12));
        let LoadedAperture::Modal(a) = p.system().valve() else {
            panic!("modal material source")
        };
        assert!(a.material_plates()[0].viscosity().is_some());
        assert!(a.stored_energy_j() > 0.);
        assert_eq!(a.accepted_steps(), 1024 * CLOCK_REFINEMENT as u64);
        out
    };
    let a = render(1., 1);
    assert_eq!(a, render(1., 37));
    assert!(
        a.iter()
            .zip(render(4., 37))
            .any(|(x, y)| (x - y).abs() > 1e-12)
    );
}
#[test]
fn explicit_zero_tensor_keeps_original_source_image_and_invalid_bindings_refuse() {
    let zero = input(0., 900.);
    let mut original = ModalAperture::from_material_plates(
        vec![(zero.0.clone(), zero.1.clone())],
        4e-4,
        gas(),
        1e6,
        1. / f64::from(MECHANICS_RATE),
        4096 * CLOCK_REFINEMENT as u64,
        limits(),
        vec![BladeState::default(); 3],
        contact(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let mut selected = make(vec![zero], 1e6).unwrap();
    for _ in 0..64 {
        assert_eq!(
            original.step(drive(5.), [0.; 2]).unwrap(),
            selected.step(drive(5.), [0.; 2]).unwrap()
        );
    }
    let mut stale = input(1., 1800.);
    stale.2 = input(1., 900.).2;
    assert!(make(vec![stale], 1e6).is_err());
    let mut duplicate = input(1., 900.);
    duplicate.1.reduction.damping_ratio = 0.1;
    assert!(make(vec![duplicate], 1e6).is_err());
    let mut exhausted = input(1., 900.);
    exhausted.3.max_projection_terms = 1;
    assert!(make(vec![exhausted], 1e6).is_err());
}

fn single_options(o: &ModalPlateOptions) -> fs_couple::bernoulli_aperture::plate::PlateApertureOptions {
    fs_couple::bernoulli_aperture::plate::PlateApertureOptions {
        assembly: o.reduction.assembly, eigenvalue_window: o.reduction.eigenvalue_window,
        mode_index: 0, eigensolver: fs_plate::SliceOptions::default(),
        max_nodes: o.max_nodes, max_triangles: o.max_triangles,
        slit_edges: o.slit_edges.clone(), rest_opening_m: 4e-4,
        damping_ratio: o.reduction.damping_ratio, max_slit_mode_variation: 0.25,
        max_slope: o.maximum_trace_slope,
    }
}
fn single(input: Input) -> Result<fs_couple::bernoulli_aperture::plate::PlateApertureReduction,
    fs_couple::acoustic_realize::AcousticRealizeError>
{
    fs_couple::bernoulli_aperture::plate::PlateApertureReduction::from_viscous_material_chart(
        input.0, single_options(&input.1), input.2, input.3, &CancelGate::new_clock_free())
}

#[test]
fn anisotropic_material_swap_reaches_single_and_paired_reeds_without_retuning() {
    use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
    use fs_couple::bernoulli_aperture::paired::BladePropagator;
    let gate = CancelGate::new_clock_free();
    let original = input(1.0, 900.0);
    let mut reference_options = original.1.reduction.clone();
    reference_options.n_modes = 1;
    let reference = original.0.certified_viscous_radiators(
        &reference_options, &original.2, original.3, &gate).unwrap();
    let mut replacement = original.clone();
    // Change only one region's anisotropic loss tensor, keeping its elasticity,
    // density, thickness, angle, mesh and every playing control unchanged.
    replacement.0 = original.0.with_region_material("part-1", &source(900.0, 16.0)).unwrap();
    replacement.2 = original.2.rebind(&replacement.0, &gate).unwrap();
    let projection_limits = original.3;
    let a = single(original).unwrap();
    let b = single(replacement).unwrap();
    assert_eq!(a.mass_kg(), b.mass_kg());
    assert_eq!(a.stiffness_n_m(), b.stiffness_n_m());
    assert_eq!(a.pressure_area_m2(), b.pressure_area_m2());
    let expected = reference.damping().matrix_per_s()[0] / (2.0 * reference.radiators()[0].omega);
    assert!((a.damping_ratio() - expected).abs() < 1e-8 * expected);
    assert!(b.damping_ratio() > a.damping_ratio());
    for p in [&a, &b] {
        assert_eq!(p.options().damping_ratio, 0.0);
        assert_eq!(p.viscosity_limits(), Some(projection_limits));
        assert_eq!(p.viscosity().unwrap().chart_identity(), p.material_chart().unwrap().specimen_identity());
    }
    let symmetric = a.clone().into_symmetric_opposed_pair().unwrap();
    assert_eq!(symmetric.damping_ratio(), a.damping_ratio());
    assert_eq!(symmetric.viscosity().unwrap().chart_identity(), a.viscosity().unwrap().chart_identity());
    let mut pair = PairedAperture::from_plates([a.clone(), b.clone()], gas(), 1e6,
        1.0 / f64::from(MECHANICS_RATE), 64, limits(), [BladeState::default(); 2], contact()).unwrap();
    assert!(pair.spec().blades[1].damping_n_s_m > pair.spec().blades[0].damping_n_s_m);
    for p in pair.plates().unwrap() { assert!(p.viscosity().is_some()); }
    for _ in 0..64 {
        let frame = pair.step(drive(5.0)).unwrap().aperture;
        assert!(frame.dissipated_energy_j >= 0.0 && frame.balance_residual_j().abs() < 1e-13);
    }
    // Both single and symmetric two-blade sources use the selected tensor in
    // the actual pressure-coupled exact propagator at the ordinary audio clock.
    for opposed in [false, true] {
        let prepare = |p: fs_couple::bernoulli_aperture::plate::PlateApertureReduction| {
            let p = if opposed { p.into_symmetric_opposed_pair().unwrap() } else { p };
            DynamicAperture::from_plate_with_propagator(p, gas().density, 1e6, 1.0 / 48000.0,
                512, ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 },
                contact(), BladePropagator::ExactZoh).unwrap()
        };
        let (mut low, mut high) = (prepare(a.clone()), prepare(b.clone()));
        let mut pressure_change = 0.0_f64;
        for n in 0..512 {
            let input = drive(if n < 256 { 5.0 } else { 0.0 });
            let x = low.step(input).unwrap();
            let y = high.step(input).unwrap();
            pressure_change = pressure_change.max((x.bore_pressure_pa - y.bore_pressure_pa).abs());
            for frame in [x, y] {
                assert!(frame.dissipated_energy_j >= 0.0 && frame.balance_residual_j().abs() < 1e-13);
            }
        }
        assert!(pressure_change > 1e-6);
    }
}

#[test]
fn single_mode_tensor_preserves_zero_loss_and_refuses_stale_or_duplicate_sources() {
    use fs_couple::bernoulli_aperture::plate::PlateApertureReduction;
    let zero = input(0.0, 900.0);
    let plain = PlateApertureReduction::from_material_chart(zero.0.clone(), single_options(&zero.1),
        &CancelGate::new_clock_free()).unwrap();
    let selected = single(zero).unwrap();
    assert!(selected.viscosity().unwrap().is_zero());
    assert_eq!(plain.damping_ratio(), selected.damping_ratio());
    assert_eq!(plain.mass_kg(), selected.mass_kg());
    assert_eq!(plain.stiffness_n_m(), selected.stiffness_n_m());
    let mut stale = input(1.0, 1800.0);
    stale.2 = input(1.0, 900.0).2;
    assert!(single(stale).is_err());
    let mut duplicate = input(1.0, 900.0);
    duplicate.1.reduction.damping_ratio = 0.01;
    assert!(single(duplicate).is_err());
    let mut exhausted = input(1.0, 900.0);
    exhausted.3.max_projection_terms = 1;
    assert!(single(exhausted).is_err());
}
