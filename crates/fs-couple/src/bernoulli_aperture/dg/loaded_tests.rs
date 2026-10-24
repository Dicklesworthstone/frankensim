//! Physical coupling regressions, not instrument-fidelity certificates.
use super::*;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicApertureSpec};
use crate::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use crate::bernoulli_aperture::modal::ModalAperture;
use crate::bernoulli_aperture::modal::pressure::MaterialPressureShape;
use crate::bernoulli_aperture::modal::relaxation::InitialModalMemory;
use crate::bernoulli_aperture::paired::{BladeState, PairedAperture, PairedApertureLimits};
use crate::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance;
use crate::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use fs_blake3::hash_domain;
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use fs_qty::{Density, Dims, Pressure};
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn load(gas: &GasState, dt: f64) -> ApertureDgLoad {
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        gas, DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    let terminal = PreparedMimoLoad::new(&MimoModel {
        ports: 1, poles: vec![PoleTerm::Real { pole: -2000.0, residue: 0.0 }],
        residues: vec![vec![C64::from_re(2e11)]], d: vec![1e8], e: vec![0.0],
        reciprocal: true,
    }, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt], sampled_passivity_points: 257,
        max_states: 4, max_memory_bytes: 1 << 20,
    }).unwrap();
    ApertureDgLoad::new_with_terminal(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal).unwrap()
}

fn no_contact() -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 0.0, 2.0,
        "manufactured loaded DG coupling: explicit zero contact".into()).unwrap()
}

fn single(load: &ApertureDgLoad, propagator: BladePropagator) -> DynamicAperture {
    // The low stiffness also admits Midpoint, so refusal below tests the DG
    // coupling requirement rather than the source's phase-resolution guard.
    let stiffness = 5.0;
    DynamicAperture::new_with_propagator(DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 0.0002, width_m: 0.001,
            closing_pressure_pa: stiffness * 0.0002 / 1e-5,
        },
        mass_kg: 1e-5, stiffness_n_m: stiffness, damping_ratio: 0.02,
        density_kg_m3: load.duct().inlet_gas().density,
        impedance_pa_s_m3: load.impedance_pa_s_m3(), time_step_s: load.propagator().dt(),
        max_steps: 256,
    }, ApertureState { opening_m: 0.0002, opening_velocity_m_s: 0.0 },
        no_contact(), propagator).unwrap()
}

fn assert_balance(frame: &DgApertureFrame) {
    let scale = frame.stored_energy_j.abs() + frame.storage_change_j.abs()
        + frame.inferred_total_dissipation_j.abs() + frame.upstream_work_j.abs()
        + frame.body_work_j.abs() + frame.aperture.mechanical_work_j.abs()
        + frame.terminal.as_ref().map_or(0.0, |t| t.port_work_j.abs());
    assert!(frame.balance_residual_j().abs() <= 1e-8 * scale.max(1e-24));
    assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-24));
}

#[test]
fn loaded_single_preserves_original_dg_arithmetic_and_strict_cold_binding() {
    let load = load(&gas(), 1.0 / 96000.0);
    let mut original = ApertureDg::new(single(&load, BladePropagator::ExactZoh), load.clone()).unwrap();
    let mut loaded = LoadedApertureDg::new_loaded(single(&load, BladePropagator::ExactZoh), load.clone()).unwrap();
    assert!(LoadedApertureDg::new_loaded(single(&load, BladePropagator::Midpoint), load.clone()).is_err());
    for step in 0..128 {
        let drive = TubeDrive {
            upstream_pressure_pa: if step < 64 { 0.5 } else { 0.0 },
            body_flow_m3_s: if step < 16 { 2e-10 } else { 0.0 },
        };
        let force = if step < 32 { 1e-6 } else { 0.0 };
        let expected = original.step_with_force(drive, force).unwrap();
        let actual = loaded.step_with_forces(drive, [force, 0.0]).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(loaded.duct_state(), original.duct_state());
        assert_eq!(loaded.load().terminal().unwrap().pole_states(), original.load().terminal().unwrap().pole_states());
        let LoadedAperture::Single(aperture) = loaded.aperture() else { unreachable!() };
        assert_eq!(aperture.state(), original.aperture().state());
        assert_balance(&actual);
    }
    let state = loaded.duct_state().clone();
    let poles = loaded.load().terminal().unwrap().pole_states().to_vec();
    assert!(loaded.step_with_forces(TubeDrive::default(), [0.0, 1.0]).is_err());
    assert_eq!(loaded.duct_state(), &state);
    assert_eq!(loaded.load().terminal().unwrap().pole_states(), poles);
    assert!(LoadedApertureDg::new_loaded(single(&load, BladePropagator::ExactZoh), loaded.load().clone()).is_err());
    let foreign = self::load(&gas(), 1.0 / 192000.0);
    assert!(LoadedApertureDg::new_loaded(single(&foreign, BladePropagator::ExactZoh), load).is_err());
    loaded.extend_step_budget(512).unwrap();
    assert_eq!(loaded.aperture().max_steps(), 512);
}

fn plate(modulus: f64, density: f64) -> PlateApertureReduction {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let root = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0).map(|(i, _)| i).collect();
    let edges = mesh.boundary_edges().into_iter().filter(|&(a, b)|
        (mesh.nodes[a].0 - 0.025).abs() < 1e-12 && (mesh.nodes[b].0 - 0.025).abs() < 1e-12)
        .map(|(a, b)| [a, b]).collect();
    let section = PlateSection::isotropic(modulus, 0.3, 0.0003, density).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, section, root, vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: edges, rest_opening_m: 0.0002, damping_ratio: 0.0,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new_clock_free()).unwrap()
}

fn paired(plates: [PlateApertureReduction; 2], load: &ApertureDgLoad) -> PairedAperture {
    let laws = std::array::from_fn::<_, 2, _>(|i| PlateRelaxationSpec {
        regions: vec![PlateRelaxationRegion {
            triangles: (0..plates[i].chart().mesh.tris.len()).collect(),
            material: GeneralizedMaxwell::new([4e9, 8e9][i],
                [vec![(2e9, 0.001), (1e9, 0.01)], vec![(3e9, 0.003)]][i].clone()).unwrap(),
            poisson_ratio: 0.3, band_hz: (0.0, 1000.0),
            provenance: "manufactured independent DG blade material".into(),
        }], max_branches: 64, max_dt_over_tau: 0.1, max_angular_step: 0.1,
    });
    let [first, second] = laws;
    PairedAperture::from_plates_with_propagator(plates, *load.duct().inlet_gas(),
        load.impedance_pa_s_m3(), load.propagator().dt(), 256,
        PairedApertureLimits { maximum_abs_displacement_m: 0.001, maximum_abs_velocity_m_s: 10.0,
            maximum_energy_j: 1.0, max_evaluations: 512, velocity_tolerance_m_s: 1e-15,
            relative_tolerance: 1e-13 },
        [BladeState { displacement_m: -3e-6, velocity_m_s: 0.0 },
         BladeState { displacement_m: 2e-6, velocity_m_s: 0.0 }], no_contact(),
        BladePropagator::ExactZoh).unwrap()
        .with_blade_relaxation(0, first, InitialApertureMemory::Unrelaxed).unwrap()
        .with_blade_relaxation(1, second, InitialApertureMemory::Relaxed).unwrap()
}

fn paired_memory(source: &PairedAperture) -> [Vec<f64>; 2] {
    std::array::from_fn(|i| source.blade_relaxation(i).unwrap().memory_sqrt_j().to_vec())
}

#[test]
fn independent_blades_and_maxwell_histories_share_the_original_dg_step() {
    let load = load(&gas(), 1.0 / 96000.0);
    let plates = [plate(4e9, 900.0), plate(8e9, 1300.0)];
    let mut reference = paired(plates.clone(), &load);
    let initial_memory = paired_memory(&reference);
    let mut acoustic = load.clone();
    let mut coupled = LoadedApertureDg::new_loaded(paired(plates, &load), load).unwrap();
    for step in 0..128 {
        let drive = TubeDrive { upstream_pressure_pa: if step < 64 { 0.5 } else { 0.0 },
            body_flow_m3_s: if step < 16 { 1e-10 } else { 0.0 } };
        let forces = if step < 32 { [1e-6, -2e-6] } else { [0.0; 2] };
        let relation = acoustic.port_response().unwrap();
        let expected = reference.step_with_forces(ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: relation.history_pressure_pa / 2.0,
            body_flow_m3_s: drive.body_flow_m3_s,
        }, forces).unwrap();
        let duct = acoustic.preview_step(expected.aperture.bore_flow_m3_s).unwrap().accept();
        let actual = coupled.step_with_forces(drive, forces).unwrap();
        assert_eq!(actual.aperture, expected.aperture);
        assert_eq!(actual.duct, duct);
        assert_eq!(coupled.duct_state(), acoustic.state());
        assert_eq!(actual.terminal.as_ref(), acoustic.last_terminal_frame());
        let LoadedAperture::Paired(retained) = coupled.aperture() else { unreachable!() };
        assert_eq!(retained.blade_states(), reference.blade_states());
        assert_eq!(paired_memory(retained), paired_memory(&reference));
        assert_balance(&actual);
    }
    let LoadedAperture::Paired(retained) = coupled.aperture() else { unreachable!() };
    assert_ne!(retained.blade_states()[0], retained.blade_states()[1]);
    let memory = paired_memory(retained);
    assert_eq!((memory[0].len(), memory[1].len()), (2, 1));
    assert_ne!(memory, initial_memory);
    let states = *retained.blade_states();
    let duct = coupled.duct_state().clone();
    let poles = coupled.load().terminal().unwrap().pole_states().to_vec();
    for (pressure, forces, spatial) in [
        (f64::NAN, [0.0; 2], [0.0; 2]),
        (0.0, [0.0, f64::INFINITY], [0.0; 2]),
        (0.0, [0.0; 2], [1.0, 0.0]),
    ] {
        assert!(coupled.step_with_forces_and_pressures(TubeDrive {
            upstream_pressure_pa: pressure, body_flow_m3_s: 0.0,
        }, forces, spatial).is_err());
        let LoadedAperture::Paired(retained) = coupled.aperture() else { unreachable!() };
        assert_eq!(*retained.blade_states(), states);
        assert_eq!(paired_memory(retained), memory);
        assert_eq!(retained.accepted_steps(), 128);
        assert_eq!(coupled.duct_state(), &duct);
        assert_eq!(coupled.load().terminal().unwrap().pole_states(), poles);
    }
}

fn material_pack() -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let artifact = hash_domain("org.frankensim.test.loaded-dg-modal", b"manufactured Maxwell plate");
    let provenance = Provenance { source: "manufactured DG material regression".into(),
        license: "CC0-1.0".into(), artifact: None };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured six-mode reed".into(), method: "physical regression".into(),
        artifact, caveats: "not a measured material".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [
        ("density", 900.0, Density::DIMS), ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE), ("boundary_modulus", 1e8, Pressure::DIMS),
        ("boundary_time", 0.001, Dims([0, 0, 1, 0, 0, 0])),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: provenance.clone(), observations: vec![observation.clone()],
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "manufactured".into(), revision: 0 }, claims, vec![]).unwrap();
    let normalized = NormalizedPack::new("synthetic-reed", "synthetic-v1", artifact,
        "CC0-1.0", card.claims().clone(), vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(card.id().clone(), normalized).unwrap()
}

fn material_source() -> MaterialModalValvePerformance {
    // Reuse the existing source geometry and its material-card ingest boundary.
    // Only its original plate images are retained by the DG fixture below.
    let text = include_str!("../../../examples/material_modal_valve/paired.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
        .replace("loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000",
            "loss proportional-maxwell relaxed 1 boundary_modulus boundary_time");
    MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material_pack())).unwrap()
}

fn modal(source: &ModalAperture, load: &ApertureDgLoad) -> ModalAperture {
    let spec = source.spec();
    let mut modal = ModalAperture::from_material_plates_with_propagator(
        source.material_plates().iter().map(|p| (p.source().clone(), p.options().clone())).collect(),
        spec.rest_opening_m, *load.duct().inlet_gas(), load.impedance_pa_s_m3(),
        load.propagator().dt(), 256, spec.limits,
        source.modal_states().iter().enumerate().map(|(i, _)| BladeState {
            displacement_m: 1e-7 / (i + 1) as f64, velocity_m_s: 0.0,
        }).collect(), no_contact(), BladePropagator::ExactZoh, &CancelGate::new_clock_free(),
    ).unwrap().with_material_relaxation(
        &[("boundary_modulus".into(), "boundary_time".into())],
        InitialModalMemory::Unrelaxed, &CancelGate::new_clock_free(),
    ).unwrap();
    for blade in 0..2 {
        let plate = &modal.material_plates()[blade];
        let shape = MaterialPressureShape::from_plate(plate,
            (0..plate.displacement_shapes().len()).map(|i| if i % 2 == 0 { 1.0 } else { -0.5 }).collect(),
            "manufactured signed pressure trace in original source-node order".into());
        modal = modal.with_material_pressure_shape(blade, shape, 1_000_000,
            &CancelGate::new_clock_free()).unwrap();
    }
    modal
}

#[test]
fn every_material_mode_memory_and_pressure_port_survives_dg_coupling_and_rollback() {
    let source = material_source();
    let LoadedAperture::Modal(seed) = source.renderer().source().system().valve() else { unreachable!() };
    let load = load(seed.gas().unwrap(), seed.spec().time_step_s);
    let mut reference = modal(seed, &load);
    let mut acoustic = load.clone();
    let mut coupled = LoadedApertureDg::new_loaded(modal(seed, &load), load.clone()).unwrap();
    let mut pressure_work = 0.0;
    for step in 0..128 {
        let drive = TubeDrive { upstream_pressure_pa: if (16..64).contains(&step) { 0.5 } else { 0.0 },
            body_flow_m3_s: 0.0 };
        let forces = if (16..32).contains(&step) { [1e-6, -2e-6] } else { [0.0; 2] };
        let pressures = if step < 64 { [3.0, -2.0] } else { [0.0; 2] };
        let relation = acoustic.port_response().unwrap();
        let expected = reference.step_with_pressures(ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: relation.history_pressure_pa / 2.0,
            body_flow_m3_s: drive.body_flow_m3_s,
        }, forces, pressures).unwrap();
        let duct = acoustic.preview_step(expected.aperture.bore_flow_m3_s).unwrap().accept();
        let actual = coupled.step_with_forces_and_pressures(drive, forces, pressures).unwrap();
        assert_eq!(actual.aperture, expected.aperture);
        assert_eq!(actual.spatial_pressure, expected.spatial_pressure);
        assert_eq!(actual.duct, duct);
        assert_eq!(coupled.duct_state(), acoustic.state());
        assert_eq!(actual.terminal.as_ref(), acoustic.last_terminal_frame());
        let LoadedAperture::Modal(retained) = coupled.aperture() else { unreachable!() };
        assert_eq!(retained.modal_states(), reference.modal_states());
        assert_eq!(retained.material_relaxation().unwrap().memory_sqrt_j(),
            reference.material_relaxation().unwrap().memory_sqrt_j());
        let work: f64 = actual.spatial_pressure.iter().map(|p| p.work_j).sum();
        pressure_work += work.abs();
        if forces == [0.0; 2] {
            assert!((work - actual.aperture.mechanical_work_j).abs() < 1e-22);
        }
        assert_balance(&actual);
    }
    assert!(pressure_work > 0.0);
    let LoadedAperture::Modal(retained) = coupled.aperture() else { unreachable!() };
    let states = retained.modal_states().to_vec();
    let memory = retained.material_relaxation().unwrap().memory_sqrt_j().to_vec();
    assert_eq!(states.len(), 6);
    assert_eq!(memory.len(), states.len());
    assert!(memory.iter().flatten().any(|z| *z != 0.0));
    assert!(retained.material_pressure_ports().iter().all(Option::is_some));
    let duct = coupled.duct_state().clone();
    let poles = coupled.load().terminal().unwrap().pole_states().to_vec();
    for (pressure, forces, spatial) in [
        (f64::NAN, [0.0; 2], [0.0; 2]),
        (0.0, [f64::INFINITY, 0.0], [0.0; 2]),
        (0.0, [0.0; 2], [0.0, f64::NAN]),
    ] {
        assert!(coupled.step_with_forces_and_pressures(TubeDrive {
            upstream_pressure_pa: pressure, body_flow_m3_s: 0.0,
        }, forces, spatial).is_err());
        let LoadedAperture::Modal(retained) = coupled.aperture() else { unreachable!() };
        assert_eq!(retained.modal_states(), states);
        assert_eq!(retained.material_relaxation().unwrap().memory_sqrt_j(), memory);
        assert_eq!(retained.accepted_steps(), 128);
        assert_eq!(coupled.duct_state(), &duct);
        assert_eq!(coupled.load().terminal().unwrap().pole_states(), poles);
    }
    scheduled_modal_parity(seed, &load);
}

fn scheduled_modal_parity(seed: &ModalAperture, load: &ApertureDgLoad) {
    use super::loaded::performance::LoadedDgPerformance;
    use crate::bernoulli_aperture::performance::{ApertureObservation, AperturePerformanceConfig};
    use crate::bernoulli_aperture::performance::force::ApertureForceEvent;
    use crate::bernoulli_aperture::two_sided::observer::schedule::spatial_pressures::MaterialPressureEvent;
    use crate::pcm_wav::observation::PressureRenderer;
    use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

    let build = || {
        let source = LoadedApertureDg::new_loaded(modal(seed, load), load.clone()).unwrap();
        let schedule = GestureSchedule::try_new(768000, vec![GestureTrack {
            id: "blowing".into(), target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.5), events: vec![],
        }]).unwrap();
        let mut performance = LoadedDgPerformance::new(source, ApertureObservation::Inlet,
            schedule, AperturePerformanceConfig {
                sample_rate_hz: 768000, samples: 128, max_block: 37,
                max_compile_work: 1 << 20, max_controls: 1024,
            }).unwrap();
        for blade in 0..2 {
            performance = performance.with_material_plate_force_controls(blade, vec![
                ApertureForceEvent { sample: 0, port: 0, force_n: [1e-6, -2e-6][blade] },
                ApertureForceEvent { sample: [32, 16][blade], port: 0, force_n: 0.0 },
            ]).unwrap().with_material_pressure_controls(blade, vec![
                MaterialPressureEvent { sample: 0, pressure_pa: [3.0, -2.0][blade] },
                MaterialPressureEvent { sample: 64, pressure_pa: 0.0 },
            ]).unwrap();
        }
        performance
    };
    let mut blocks = build();
    let mut singles = build();
    let mut direct = LoadedApertureDg::new_loaded(modal(seed, load), load.clone()).unwrap();
    let pending = blocks.pending_controls().to_vec();
    let pressure_pending = blocks.pending_material_pressure_controls(0).to_vec();
    let force_pending = blocks.pending_plate_force_controls(1).to_vec();
    let cold_duct = blocks.system().duct_state().clone();
    let mut rejected = [123.0; 38];
    assert!(blocks.block(&mut rejected).is_err());
    assert_eq!(rejected, [123.0; 38]);
    assert_eq!(blocks.system().duct_state(), &cold_duct);
    assert_eq!(blocks.system().aperture().accepted_steps(), 0);
    assert_eq!(blocks.pending_controls(), pending);
    assert_eq!(blocks.pending_material_pressure_controls(0), pressure_pending);
    assert_eq!(blocks.pending_plate_force_controls(1), force_pending);
    assert_eq!(blocks.held_generalized_forces_n(), [0.0; 2]);
    assert_eq!(blocks.held_material_pressures_pa(), [0.0; 2]);
    assert!(blocks.last_frame().is_none());

    let mut expected = Vec::with_capacity(128);
    for sample in 0..128 {
        let frame = direct.step_with_forces_and_pressures(
            TubeDrive { upstream_pressure_pa: 0.5, body_flow_m3_s: 0.0 },
            [if sample < 32 { 1e-6 } else { 0.0 }, if sample < 16 { -2e-6 } else { 0.0 }],
            if sample < 64 { [3.0, -2.0] } else { [0.0; 2] },
        ).unwrap();
        let mut observed = [0.0];
        singles.block(&mut observed).unwrap();
        assert_eq!(observed[0].to_bits(), frame.aperture.bore_pressure_pa.to_bits());
        assert_eq!(singles.last_frame(), Some(&frame));
        expected.push(observed[0]);
    }
    let mut actual = vec![0.0; 128];
    for block in actual.chunks_mut(37) { blocks.block(block).unwrap(); }
    assert_eq!(actual, expected);
    assert_eq!(blocks.last_frame(), singles.last_frame());
    assert_eq!(blocks.system().duct_state(), direct.duct_state());
    assert_eq!(blocks.applied_controls(), singles.applied_controls());
    for blade in 0..2 {
        assert!(blocks.pending_plate_force_controls(blade).is_empty());
        assert!(blocks.pending_material_pressure_controls(blade).is_empty());
        assert_eq!(blocks.applied_plate_force_controls(blade), singles.applied_plate_force_controls(blade));
        assert_eq!(blocks.applied_material_pressure_controls(blade), singles.applied_material_pressure_controls(blade));
    }
    let LoadedAperture::Modal(a) = blocks.system().aperture() else { unreachable!() };
    let LoadedAperture::Modal(b) = direct.aperture() else { unreachable!() };
    assert_eq!(a.modal_states(), b.modal_states());
    assert_eq!(a.material_relaxation().unwrap().memory_sqrt_j(),
        b.material_relaxation().unwrap().memory_sqrt_j());
}
