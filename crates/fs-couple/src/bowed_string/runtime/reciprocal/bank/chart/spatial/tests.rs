//! Manufactured source/mesh regressions, not measured instrument validation.
use super::*;
use fs_evidence::ValidityDomain;
use fs_scenario::gesture::GestureTarget;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::GasSpec;
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as V;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_scenario::{HelmholtzCavity, gesture::{GestureTrack, GestureValue}};
use crate::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use crate::string_specimen::{with_uniform_circular_material_state,
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use crate::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint,
    compile_plate_material_chart};

fn scoped(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena, fs_exec::StreamKey {
        seed: 7, kernel_id: 83, tile: 0, iteration: 0,
    }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate));
}
fn material(viscosity: f64) -> ResolvedMaterialStatePoint {
    let mut values = vec![("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, 200e9), ("poisson_ratio", Dims::NONE, 0.3),
        ("young_modulus_1", Pressure::DIMS, 8e9), ("young_modulus_2", Pressure::DIMS, 0.8e9),
        ("poisson_ratio_12", Dims::NONE, 0.3), ("shear_modulus_12", Pressure::DIMS, 0.6e9),
        (ETA, DynViscosity::DIMS, 10.0)];
    values.extend(V.into_iter().zip([1e6, 7e6, 2e6, 0.4e6, 0.1e6, -0.2e6])
        .map(|(key, value)| (key, DynViscosity::DIMS, value * viscosity)));
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic spatial body".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn sources() -> Vec<MaterialStringSource> {
    [3, 4].map(|n_modes| MaterialStringSource::Homogeneous(
        with_uniform_circular_material_state(PrestressedString {
            length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
            width_m: 0.01, n_modes, damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None, relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        }, 0.0003, &material(0.0)).unwrap().with_kelvin_voigt_bending_loss().unwrap())).to_vec()
}
fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap() }
fn cavity() -> MaterialBodyCavity {
    MaterialBodyCavity { cavity: HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3: 0.004, neck_radius_m: 0.01, neck_length_m: 0.02,
        inner_end: fs_scenario::NeckEnd::Flanged, outer_end: fs_scenario::NeckEnd::Flanged,
    }, gas(), CavityNeckLoss::Resistance(100.0)).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.5 }
}
fn point(node: usize) -> Vec<f64> { let mut w = vec![0.0; 30]; w[node] = 1.0; w }
fn footprints() -> PlateBridgeFootprints {
    PlateBridgeFootprints { weights: vec![point(14), point(15)], max_projection_terms: 100_000 }
}
fn body() -> SharedMaterialPlateChartBody {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 5, 4);
    let material = material(1.0);
    let regions = ["left", "right"].into_iter().enumerate().map(|(r, name)| PlateRegionMaterial {
        region: PlateRegion { name: name.into(), triangle_indices: mesh.tris.iter().enumerate()
            .filter_map(|(e, tri)| {
                let right = tri.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0 > 0.2;
                (right == (r == 1)).then_some(e)
            }).collect() },
        material: material.clone(), model: PlateMaterialModel::Orthotropic12,
        thickness_constraint: PlateThicknessConstraint::FixedThickness(0.002),
        material_angle_rad: if r == 0 { 0.22 } else { -0.3 },
    }).collect();
    let specimen = compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(5, 4), regions).unwrap();
    SharedMaterialPlateChartBody { specimen, reduction: PlateChartRadiation {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
        eigenvalue_window: (1.0, 1e8), n_modes: 4, damping_ratio: 0.0, unit_force_weights: point(14),
    }, bridges: vec![SharedPlateBridge { station_fraction: 0.8, stiffness_n_m: 3000.0,
        damping_n_s_m: 0.1, rest_extension_m: 0.0 }; 2],
        max_plate_dofs: 90, max_plate_elements: 40,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0 }
}
fn options() -> StringBankOptions {
    StringBankOptions { sample_rate_hz: 48_000, bow: StringBowLocation { string: 0, station_fraction: 0.2 },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig { max_modes: 32, max_connections: 16, max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 },
        friction: ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 } }
}
fn limits() -> PlateViscosityLimits { PlateViscosityLimits {
    max_modes: 8, max_projection_terms: 1_000_000, reference_resistance_n_s_m: 1.0,
} }
fn law() -> FrictionLaw { FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } }
fn drive() -> FrictionDrive { FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.0001 } }
fn build(cx: &Cx<'_>, gate: &CancelGate, b: SharedMaterialPlateChartBody,
    ports: Option<PlateBridgeFootprints>, tensor: bool) -> ReciprocalStringBank {
    let viscosity = tensor.then(|| (PlateViscosity::from_chart(&b.specimen, V, gate).unwrap(), limits()));
    ReciprocalStringBank::build_material_plate_chart(cx, sources(), gas(), b, Some(cavity()), law(),
        options(), gate, viscosity, ports).unwrap()
}

#[test]
fn repeated_reference_footprints_preserve_both_original_material_paths_exactly() {
    scoped(|cx, gate| {
        for tensor in [false, true] {
            let b = body(); let weights = vec![b.reduction.unit_force_weights.clone(); 2];
            let mut a = build(cx, gate, b.clone(), None, tensor);
            let mut b = build(cx, gate, b, Some(PlateBridgeFootprints { weights, max_projection_terms: 100_000 }), tensor);
            assert!(a.source().material_chart_bridges().is_none());
            for source in [a.source(), b.source()] {
                assert_eq!(source.bridge_surface_materials(0).unwrap(),
                    [material(0.0).material().clone(), material(1.0).material().clone()]);
            }
            let zero = vec![0.0; a.system().mode_count()];
            for _ in 0..128 { assert_eq!(a.step(&zero, drive()).unwrap(), b.step(&zero, drive()).unwrap()); }
        }
    });
}

#[test]
fn spatial_connections_full_damping_and_radiation_share_the_original_reduced_basis() {
    scoped(|cx, gate| {
        let b = body(); let input = footprints();
        let viscosity = PlateViscosity::from_chart(&b.specimen, V, gate).unwrap();
        let (image, ports) = b.specimen.certified_viscous_radiators_with_force_ports(
            &b.reduction, &input.weights, input.max_projection_terms, &viscosity, limits(), gate).unwrap();
        let bank = ReciprocalStringBank::from_viscous_material_plate_chart_with_force_ports(cx, sources(), gas(),
            b, input.clone(), viscosity, limits(), Some(cavity()), law(), options(), gate).unwrap();
        assert_eq!(bank.system().components().len(), 4); // two strings, ONE body, ONE neck
        assert_eq!(bank.source().material_chart_bridges(), Some(&input));
        assert_eq!(bank.source().bridge_surface_materials(1).unwrap(),
            [material(0.0).material().clone(), material(1.0).material().clone()]);
        for (i, expected) in ports.iter().enumerate() {
            assert_eq!(&bank.source().strings()[i].connections[0].body.shapes, expected);
        }
        assert_ne!(ports[0], ports[1]);
        assert_eq!(bank.source().material_chart_viscosity().unwrap().damping().matrix_per_s(), image.damping().matrix_per_s());
        assert_eq!(bank.source().body_connections().len(), 2 + image.damping().channel_count());
        assert_eq!(bank.source().gas(), bank.source().material_cavity().unwrap().cavity.gas());
        let weights = bank.system().acceleration_pressure_weights().unwrap();
        for (i, r) in image.radiators().iter().enumerate() {
            assert_eq!(weights[7 + i], r.radiate_mass_normalized(1.0, gas().density, 1.0));
        }
    });
}

#[test]
fn opposite_string_loads_do_not_cancel_at_distinct_physical_body_contacts() {
    scoped(|cx, gate| {
        let mut s = sources(); s[1] = s[0].clone();
        let make = |ports| ReciprocalStringBank::build_material_plate_chart(cx, s.clone(), gas(), body(),
            Some(cavity()), law(), options(), gate, None, ports).unwrap();
        let mut common = make(None); let mut spatial = make(Some(footprints()));
        let mut forces = vec![0.0; spatial.system().mode_count()];
        forces[0] = 1e-4; forces[3] = -1e-4;
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        for _ in 0..128 {
            for bank in [&mut common, &mut spatial] {
                let frame = bank.step(&forces, released).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
        let motion = |bank: &ReciprocalStringBank| bank.system().components()[2].states().iter()
            .map(|s| s.displacement_m_sqrt_kg.abs()).sum::<f64>();
        assert!(motion(&spatial) > 1e-15);
        assert!(motion(&common) < 1e-6 * motion(&spatial));
    });
}

#[test]
fn moving_only_one_cold_body_contact_changes_sympathetic_motion_and_physical_pressure() {
    scoped(|cx, gate| {
        let mut changed = footprints(); changed.weights[1] = point(8);
        let mut a = build(cx, gate, body(), Some(footprints()), true);
        let mut b = build(cx, gate, body(), Some(changed), true);
        assert_eq!(a.system().components()[2].modes(), b.system().components()[2].modes());
        assert_eq!(a.source().material_chart_viscosity().unwrap().damping().matrix_per_s(),
            b.source().material_chart_viscosity().unwrap().damping().matrix_per_s());
        let zero = vec![0.0; a.system().mode_count()];
        let (mut movement, mut pressure) = (0.0, 0.0);
        for _ in 0..256 {
            let pa = a.step(&zero, drive()).unwrap().observer_pressure_pa;
            let pb = b.step(&zero, drive()).unwrap().observer_pressure_pa;
            pressure += (pa - pb).abs();
            movement += (a.system().components()[1].states()[0].displacement_m_sqrt_kg
                - b.system().components()[1].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(movement > 1e-18 && pressure > 1e-14);
    });
}

#[test]
fn spatial_sources_survive_schedule_partitioning_and_failed_samples() {
    scoped(|cx, gate| {
        let mut expected = None;
        for block in [1, 37, 128] {
            let bank = build(cx, gate, body(), Some(footprints()), true);
            let schedule = GestureSchedule::try_new(700, vec![GestureTrack { id: "bow".into(),
                target: GestureTarget::BowStroke { string: 0 }, initial: GestureValue::Bow {
                    velocity_m_per_s: drive().speed_m_s, normal_force_n: drive().normal_force_n, station: 0.2 }, events: vec![],
            }]).unwrap();
            let mut p = bank.into_scheduled(&schedule, "bow", FrictionGestureConfig {
                sample_rate_hz: 48_000, samples: 129, max_block: 128, max_work: 1000,
                max_events: 32, station_fraction: 0.2,
            }).unwrap();
            assert_eq!(p.source().material_chart_bridges(), Some(&footprints()));
            assert!(p.source().material_chart_viscosity().is_some());
            let mut out = vec![0.0; 129];
            for chunk in out.chunks_mut(block) { p.renderer_mut().block(chunk).unwrap(); }
            if let Some(reference) = &expected { assert_eq!(&out, reference); } else { expected = Some(out); }
        }
        let mut bank = build(cx, gate, body(), Some(footprints()), true);
        let zero = vec![0.0; bank.system().mode_count()]; bank.step(&zero, drive()).unwrap();
        let states: Vec<_> = bank.system().components().iter().map(|m| m.states().to_vec()).collect();
        let count = bank.system().samples_rendered();
        assert!(bank.step(&zero, FrictionDrive { speed_m_s: f64::NAN, ..drive() }).is_err());
        let stop = CancelGate::new_clock_free(); stop.request();
        assert!(bank.step_under_gate(&zero, drive(), &stop).is_err());
        assert_eq!(bank.system().samples_rendered(), count);
        assert_eq!(bank.system().components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>(), states);
        assert_eq!(bank.source().material_chart_bridges(), Some(&footprints()));
    });
}

#[test]
fn missing_invalid_detached_and_unfunded_footprints_refuse_without_fallback() {
    scoped(|cx, gate| {
        for tensor in [false, true] { for defect in 0..6 {
            let b = body(); let mut ports = footprints();
            match defect {
                0 => { ports.weights.pop(); },
                1 => { ports.weights[1].pop(); },
                2 => ports.weights[1][0] = f64::NAN,
                3 => ports.weights[1][0] = 1.0, // sum two, never normalized
                4 => ports.weights[1] = point(0), // pinned body point has no motion
                _ => ports.max_projection_terms = 0,
            }
            let viscous = tensor.then(|| (PlateViscosity::from_chart(&b.specimen, V, gate).unwrap(), limits()));
            assert!(ReciprocalStringBank::build_material_plate_chart(cx, sources(), gas(), b,
                Some(cavity()), law(), options(), gate, viscous, Some(ports)).is_err());
        } }
    });
}
