//! Synthetic material/mesh regressions, not measured instrument validation.
use fs_couple::bowed_string::runtime::reciprocal::{InteriorStringConnection, plate::MaterialBodyCavity};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember, MaterialStringSource, ReciprocalStringBank, StringBankOptions, StringBowLocation,
    chart::SharedMaterialPlateChartBody, plate::SharedPlateBridge,
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig};
use fs_couple::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use fs_couple::string_specimen::{with_uniform_circular_material_state, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_couple::thin_plate::{PlateChartRadiation, PlateMaterialModel, PlateRegionMaterial,
    PlateThicknessConstraint, certified_chart_radiators, compile_plate_material_chart};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_math::{c64::C64, det::sqrt};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{HelmholtzCavity, PrestressedString,
    gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena, fs_exec::StreamKey {
        seed: 7, kernel_id: 83, tile: 0, iteration: 0,
    }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate))
}
fn material(rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, 0.3), (ETA, DynViscosity::DIMS, 10.0)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic regional body".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn sources() -> Vec<MaterialStringSource> {
    [4, 3].into_iter().map(|n_modes| {
        MaterialStringSource::Homogeneous(with_uniform_circular_material_state(PrestressedString {
            length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
            width_m: 0.01, n_modes, damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None, relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        }, 0.0003, &material(7800.0, 200e9)).unwrap().with_kelvin_voigt_bending_loss().unwrap())
    }).collect()
}
fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn cavity(gas: GasState) -> MaterialBodyCavity {
    MaterialBodyCavity { cavity: HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3: 0.004, neck_radius_m: 0.01, neck_length_m: 0.02,
        inner_end: fs_scenario::NeckEnd::Flanged, outer_end: fs_scenario::NeckEnd::Flanged,
    }, gas, CavityNeckLoss::RadiationAtResonance).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.5 }
}
fn body() -> SharedMaterialPlateChartBody {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 5, 4);
    let regions = ["left", "right"].into_iter().enumerate().map(|(r, name)| PlateRegionMaterial {
        region: PlateRegion { name: name.into(), triangle_indices: mesh.tris.iter().enumerate()
            .filter_map(|(e, tri)| {
                let right = tri.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0 > 0.2;
                (right == (r == 1)).then_some(e)
            }).collect() },
        material: material(if r == 0 { 7800.0 } else { 3900.0 }, 200e9),
        model: PlateMaterialModel::Isotropic,
        thickness_constraint: PlateThicknessConstraint::FixedThickness(0.001), material_angle_rad: 0.0,
    }).collect();
    let mut force = vec![0.0; mesh.nodes.len()];
    force[2 * 6 + 2] = 0.75; // two physical interior nodes, not modal coefficients
    force[2 * 6 + 3] = 0.25;
    let specimen = compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(5, 4), regions).unwrap();
    SharedMaterialPlateChartBody { specimen, reduction: PlateChartRadiation {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
        eigenvalue_window: (1.0, 1e8), n_modes: 3, damping_ratio: 0.002, unit_force_weights: force,
    }, bridges: vec![SharedPlateBridge { station_fraction: 0.8, stiffness_n_m: 3000.0,
        damping_n_s_m: 0.1, rest_extension_m: 0.0 }; 2],
        max_plate_dofs: 90, max_plate_elements: 40,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0 }
}
fn options() -> StringBankOptions {
    StringBankOptions { sample_rate_hz: 48_000, bow: StringBowLocation { string: 0, station_fraction: 0.2 },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig { max_modes: 32, max_connections: 8, max_setup_terms: 16_384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 },
        friction: ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 } }
}
fn law() -> FrictionLaw { FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } }
fn drive() -> FrictionDrive { FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 } }
fn build(cx: &Cx<'_>, gate: &CancelGate, body: SharedMaterialPlateChartBody, gas: GasState) -> ReciprocalStringBank {
    ReciprocalStringBank::from_material_plate_chart(cx, sources(), gas, body, Some(cavity(gas)), law(), options(), gate).unwrap()
}

#[test]
fn regional_lowering_matches_the_existing_network_and_physical_pressure_exactly() {
    with_cx(|cx, gate| {
        let body = body(); let gas = gas(101_325.0); let cavity = cavity(gas);
        let radiators = certified_chart_radiators(body.specimen.chart(), &body.reduction).unwrap();
        let bridge: Vec<_> = radiators.iter().map(|r| r.drive_participation / sqrt(r.mass_kg)).collect();
        let areas: Vec<_> = radiators.iter().map(|r| r.area_m2 / sqrt(r.mass_kg)).collect();
        let receiver = ModalAcousticTimeModel::try_new(48_000, radiators.iter().map(|r| ModalAcousticMode {
            angular_frequency_rad_s: r.omega, damping_ratio: r.zeta, pressure_per_modal_velocity: C64::ZERO,
        }).collect(), body.modal_budget).unwrap();
        let members = sources().into_iter().zip(&body.bridges).map(|(source, b)| MaterialStringMember {
            source, connections: vec![InteriorStringConnection { station_fraction: b.station_fraction,
                body: ModalAttachment { component: 0, shapes: bridge.clone() }, stiffness_n_m: b.stiffness_n_m,
                damping_n_s_m: b.damping_n_s_m, rest_extension_m: b.rest_extension_m }],
        }).collect();
        let raw = ReciprocalStringBank::with_body_connections(cx, members, gas,
            vec![receiver, cavity.cavity.neck_component(48_000, cavity.modal_budget).unwrap()],
            cavity.cavity.connections(0, 1, &areas).unwrap(), law(), options(), gate).unwrap();
        let mut weights = vec![0.0; 7]; // deliberately unequal string mode counts
        weights.extend(radiators.iter().map(|r| r.radiate_mass_normalized(1.0, gas.density, body.listener_m)));
        weights.push(gas.density * cavity.cavity.neck_volume_shape_m2_per_sqrt_kg()
            / (core::f64::consts::TAU * cavity.listener_m));
        let mut reference = raw.into_parts().1.with_acceleration_pressure(weights.clone()).unwrap();
        let mut actual = build(cx, gate, body, gas);
        assert_eq!(actual.system().acceleration_pressure_weights().unwrap(), weights);
        assert_eq!(actual.system().components().len(), 4);
        assert_eq!(actual.system().components()[0].modes().len(), 4);
        assert_eq!(actual.system().components()[1].modes().len(), 3);
        let zero = vec![0.0; actual.system().mode_count()];
        for _ in 0..128 { assert_eq!(actual.step(&zero, drive()).unwrap(), reference.step(&zero, drive()).unwrap()); }
        assert!(actual.system().components()[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-10));
        let before = actual.system().total_energy_j().unwrap();
        for _ in 0..128 {
            let frame = actual.step(&zero, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }).unwrap();
            assert_eq!(frame.drive_work_j, 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(actual.system().total_energy_j().unwrap() <= before + 1e-10);
    });
}

#[test]
fn replacing_one_material_region_changes_actual_string_motion_and_keeps_other_sources() {
    with_cx(|cx, gate| {
        let a_body = body(); let mut b_body = a_body.clone();
        let left = a_body.specimen.regions()[0].input.material.identity();
        b_body.specimen = b_body.specimen.with_region_material("right", &material(15_600.0, 100e9)).unwrap();
        assert_ne!(a_body.specimen.specimen_identity(), b_body.specimen.specimen_identity());
        assert_eq!(left, b_body.specimen.regions()[0].input.material.identity());
        let mut a = build(cx, gate, a_body, gas(101_325.0));
        let mut b = build(cx, gate, b_body, gas(101_325.0));
        assert_eq!(a.cards()[0].fundamental_hz(), b.cards()[0].fundamental_hz());
        assert_ne!(a.system().components()[2].modes(), b.system().components()[2].modes());
        let zero = vec![0.0; a.system().mode_count()]; let mut difference = 0.0;
        for _ in 0..256 {
            for model in [&mut a, &mut b] {
                let frame = model.step(&zero, drive()).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(difference > 1e-14, "regional material did not react back: {difference}");
        assert_eq!(b.source().material_chart().unwrap().specimen.regions()[0].input.material.identity(), left);
    });
}

#[test]
fn one_gas_reaches_strings_body_observation_and_shared_cavity() {
    with_cx(|cx, gate| {
        let low = gas(50_662.5); let high = gas(101_325.0);
        let mut a = build(cx, gate, body(), low); let mut b = build(cx, gate, body(), high);
        assert_ne!(a.cards()[0].zetas, b.cards()[0].zetas);
        assert_eq!(a.source().gas(), a.source().material_cavity().unwrap().cavity.gas());
        let wa = a.system().acceleration_pressure_weights().unwrap();
        let wb = b.system().acceleration_pressure_weights().unwrap();
        for k in 7..wa.len()-1 { assert!((wa[k] - 0.5 * wb[k]).abs() < 1e-12); }
        let zero = vec![0.0; a.system().mode_count()]; let mut difference = 0.0;
        for _ in 0..128 {
            difference += (a.step(&zero, drive()).unwrap().observer_pressure_pa
                - b.step(&zero, drive()).unwrap().observer_pressure_pa).abs();
        }
        assert!(difference > 1e-12);
    });
}

#[test]
fn scheduler_keeps_regional_sources_and_exact_sample_partitioning() {
    with_cx(|cx, gate| {
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 }, events: vec![],
        }]).unwrap();
        let config = FrictionGestureConfig { sample_rate_hz: 48_000, samples: 129,
            max_block: 128, max_work: 1000, max_events: 16, station_fraction: 0.2 };
        let mut direct = build(cx, gate, body(), gas(101_325.0));
        let identity = direct.source().material_chart().unwrap().specimen.specimen_identity();
        let zero = vec![0.0; direct.system().mode_count()];
        let expected: Vec<_> = (0..129).map(|_| direct.step(&zero, drive()).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1, 37, 128] {
            let mut performance = build(cx, gate, body(), gas(101_325.0)).into_scheduled(&schedule, "bow", config).unwrap();
            assert_eq!(performance.source().material_chart().unwrap().specimen.specimen_identity(), identity);
            assert!(performance.source().material_plate().is_none());
            let mut output = vec![0.0; 129];
            for block in output.chunks_mut(partition) { performance.renderer_mut().block(block).unwrap(); }
            assert_eq!(output.iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
        }
    });
}

#[test]
fn invalid_mesh_modal_bridge_and_medium_inputs_refuse() {
    with_cx(|cx, gate| {
        for defect in 0..12 {
            let mut body = body(); let mut o = options(); let medium = gas(101_325.0);
            let mut c = cavity(medium);
            match defect {
                0 => body.max_plate_dofs = 89,
                1 => body.max_plate_elements = 39,
                2 => body.reduction.n_modes = 0,
                3 => body.reduction.damping_ratio = f64::NAN,
                4 => body.reduction.eigenvalue_window = (2.0, 1.0),
                5 => { body.reduction.unit_force_weights.fill(0.0); body.reduction.unit_force_weights[0] = 1.0; }
                6 => body.reduction.unit_force_weights[0] = 2.0,
                7 => { body.bridges.pop(); }
                8 => body.bridges[0].station_fraction = 1.0,
                9 => body.bridges[0].damping_n_s_m = -0.1,
                10 => o.network.max_modes = 7,
                _ => c = cavity(gas(50_662.5)),
            }
            assert!(ReciprocalStringBank::from_material_plate_chart(cx, sources(), medium,
                body, Some(c), law(), o, gate).is_err(), "defect {defect} was accepted");
        }
    });
}

#[test]
fn cancelled_step_leaves_source_and_every_mechanical_component_unchanged() {
    with_cx(|cx, gate| {
        let mut bank = build(cx, gate, body(), gas(101_325.0));
        let zero = vec![0.0; bank.system().mode_count()];
        bank.step(&zero, drive()).unwrap();
        let states: Vec<_> = bank.system().components().iter().map(|m| m.states().to_vec()).collect();
        let identity = bank.source().material_chart().unwrap().specimen.specimen_identity();
        let count = bank.system().samples_rendered();
        let stop = CancelGate::new_clock_free(); stop.request();
        assert!(bank.step_under_gate(&zero, drive(), &stop).is_err());
        assert_eq!(bank.system().samples_rendered(), count);
        assert_eq!(bank.source().material_chart().unwrap().specimen.specimen_identity(), identity);
        assert_eq!(bank.system().components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>(), states);
        assert!(ReciprocalStringBank::from_material_plate_chart(cx, sources(), gas(101_325.0),
            body(), None, law(), options(), &stop).is_err());
    });
}

#[path = "material_plate_chart_bank/multiple.rs"]
mod multiple;

#[path = "material_plate_chart_bank/relaxation.rs"]
mod relaxation;

#[path = "material_plate_chart_bank/radiation.rs"]
mod radiation;
