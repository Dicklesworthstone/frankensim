//! One material key on each acoustic side of the original three valve owners.
use super::*;
use fs_couple::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::network::boundaries::{BoundaryAperture, BoundaryDrive};
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_couple::bernoulli_aperture::viscous::ViscousSlit;
use fs_material::visco::GeneralizedMaxwell;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

const SIDES: [ValveSide; 2] = [ValveSide::Upstream, ValveSide::Downstream];

#[derive(Clone, Copy, Debug)]
enum Family { Single, Paired, Modal }
const FAMILIES: [Family; 3] = [Family::Single, Family::Paired, Family::Modal];

fn geometry(temperature: f64) -> ApertureDuctSpec {
    let regional = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap();
    let gases = vec![gas(), gas(), regional];
    let sections = [([0, 1], 8, 0.007), ([1, 2], 12, 0.006), ([1, 3], 6, 0.003)]
        .into_iter().zip(&gases).map(|((nodes, delay, radius_m), gas)| TubeSection {
            nodes, length_m: delay as f64 * gas.sound_speed * DT,
            radius_m, max_length_error_m: 1e-12,
        }).collect();
    ApertureDuctSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
                NetworkNode::Termination { reflection: -0.8 }, NetworkNode::Driven],
            sections, sound_speed_m_s: gases[0].sound_speed, max_wave_memory_bytes: 1 << 20,
        }, gases, viscothermal: vec![],
    }
}

fn unloaded(family: Family, steps: u64) -> TwoSidedAperture {
    let ducts = [geometry(303.15), geometry(313.15)];
    let (g, z) = combined_inlet_state(&ducts).unwrap();
    let valve: LoadedAperture = match family {
        Family::Single => DynamicAperture::new_with_propagator(DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013,
                closing_pressure_pa: 4000.0 },
            mass_kg: 1e-4, stiffness_n_m: 50.0, damping_ratio: 0.1,
            density_kg_m3: g.density, impedance_pa_s_m3: z,
            time_step_s: DT, max_steps: steps,
        }, ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 },
            contact(), BladePropagator::ExactZoh).unwrap().into(),
        Family::Paired => paired(z, g, steps, 1.7).into(),
        Family::Modal => {
            // These retained modes and loaded damping lie inside the original
            // midpoint owner's 0.04-radian band at this test's 48 kHz clock.
            let modes = vec![
                ApertureMode { mass_kg: 1e-4, stiffness_n_m: 50.0,
                    damping_n_s_m: 0.005, pressure_area_m2: 5e-6,
                    gap_projection: 1.0, force_projection: [1.0, 0.0] },
                ApertureMode { mass_kg: 8e-5, stiffness_n_m: 90.0,
                    damping_n_s_m: 0.004, pressure_area_m2: -3e-6,
                    gap_projection: -0.2, force_projection: [-0.3, 0.7] },
            ];
            ModalAperture::new(ModalApertureSpec {
                rest_opening_m: 4e-4, width_m: 0.013, modes,
                density_kg_m3: g.density, impedance_pa_s_m3: z,
                time_step_s: DT, max_steps: steps,
                limits: PairedApertureLimits {
                    maximum_abs_displacement_m: 0.01, maximum_abs_velocity_m_s: 100.0,
                    maximum_energy_j: 1.0, max_evaluations: 256,
                    velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11,
                },
            }, vec![BladeState::default(); 2], contact()).unwrap().into()
        }
    };
    TwoSidedAperture::new(valve, ducts, limits(), &CancelGate::new_clock_free()).unwrap()
}

fn plate() -> PlateApertureReduction {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let chart = PlateChart::with_boundary_and_regions(mesh,
        PlateSection::isotropic(4e9, 0.3, 0.0003, 900.0).unwrap(),
        vec![0, 7, 14], vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)),
        mode_index: 0, eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: vec![[6, 13], [13, 20]], rest_opening_m: 0.0002,
        damping_ratio: 0.0, max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new_clock_free()).unwrap()
}

fn boundary(plate: &PlateApertureReduction, gas: GasState, z: f64, steps: u64) -> BoundaryAperture {
    let relaxation = PlateRelaxationSpec { regions: vec![PlateRelaxationRegion {
        triangles: (0..plate.chart().mesh.tris.len()).collect(),
        material: GeneralizedMaxwell::new(4e9, vec![(2e9, 0.001), (1e9, 0.01)]).unwrap(),
        poisson_ratio: 0.3, band_hz: (0.0, 1000.0),
        provenance: "manufactured two-sided key material; not measured instrument data".into(),
    }], max_branches: 2, max_dt_over_tau: 0.1, max_angular_step: 0.1 };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e10, 2.0,
        "manufactured two-sided key contact".into()).unwrap().with_internal_loss(0.5).unwrap();
    let aperture = DynamicAperture::from_plate_with_propagator(plate.clone(), gas.density,
        z, DT, steps, ApertureState { opening_m: plate.options().rest_opening_m,
            opening_velocity_m_s: 0.0 }, lay, BladePropagator::ExactZoh).unwrap()
        .with_plate_relaxation(relaxation, InitialApertureMemory::Relaxed).unwrap();
    BoundaryAperture { node: 3, aperture }
}

fn model(plate: &PlateApertureReduction, family: Family, steps: u64) -> TwoSidedAperture {
    let mut model = unloaded(family, steps);
    for side in SIDES {
        let key = boundary(plate, model.section_gases(side)[2],
            model.represented_sections(side)[2].impedance_pa_s_m3, steps);
        model = model.with_boundary_apertures(side, vec![key], 1,
            &CancelGate::new_clock_free()).unwrap();
    }
    model
}

fn key_drives(sample: u64) -> [[BoundaryDrive; 1]; 2] {
    std::array::from_fn(|side| [BoundaryDrive {
        node: 3,
        exterior_pressure_pa: if sample < 512 { 1.0 + side as f64 } else { 0.0 },
        opening_force_n: if (128..768).contains(&sample) { -0.02 - 0.01 * side as f64 } else { 0.0 },
    }])
}

#[test]
fn single_paired_and_modal_sources_keep_both_material_keys_in_the_complete_balance() {
    let plate = plate();
    let gate = CancelGate::new_clock_free();
    for family in FAMILIES {
        let mut model = model(&plate, family, 2048);
        let mut closed = [false; 2];
        let mut reopened = [false; 2];
        let mut work = [0.0; 2];
        for side in SIDES {
            assert_eq!(model.boundary_apertures(side).unwrap().gas(3), Some(&model.section_gases(side)[2]));
            assert_ne!(model.section_gases(side)[0], model.section_gases(side)[2]);
        }
        for sample in 0..2048 {
            let before = model.stored_energy_j();
            let keys = key_drives(sample);
            let frame = model.step_with_boundaries_under_gate(
                drive(if sample < 512 { 40.0 } else { 0.0 }), [&keys[0], &keys[1]], &gate).unwrap();
            let accounts = frame.boundaries.map(Option::unwrap);
            assert_eq!(frame.upstream, accounts[0].network);
            assert_eq!(frame.downstream, accounts[1].network);
            let scale = before + frame.stored_energy_j + frame.dissipated_energy_j.abs()
                + frame.source_work_j.abs() + frame.aperture.mechanical_work_j.abs()
                + accounts.iter().map(|f| f.exterior_work_j.abs() + f.mechanical_work_j.abs()).sum::<f64>();
            assert!(frame.balance_residual_j().abs() <= 1e-16 + 1e-7 * scale, "{family:?}: {frame:?}");
            assert!((frame.storage_change_j - (frame.stored_energy_j - before)).abs() <= 1e-16 + 1e-7 * scale);
            assert!((model.stored_energy_j() - frame.stored_energy_j).abs()
                <= 16.0 * f64::EPSILON * frame.stored_energy_j.max(1e-20));
            for (index, side) in SIDES.into_iter().enumerate() {
                let key = model.boundary_apertures(side).unwrap().aperture_frame(3).unwrap();
                assert_eq!(key.step, frame.aperture.step);
                assert_eq!(model.node_frame(side, 3).unwrap().pressure_pa, key.bore_pressure_pa);
                work[index] += accounts[index].mechanical_work_j.abs();
                closed[index] |= (128..768).contains(&sample) && key.midpoint_opening_m == 0.0;
                reopened[index] |= sample >= 768 && key.midpoint_opening_m > 1e-4;
            }
            if sample >= 768 {
                assert!(frame.stored_energy_j <= before + 1e-16 + 1e-7 * before);
            }
        }
        assert!(closed.into_iter().all(|x| x), "both material keys close for {family:?}");
        assert!(reopened.into_iter().all(|x| x), "both material keys reopen for {family:?}");
        assert!(work.into_iter().all(|x| x > 0.0));
        for side in SIDES {
            assert!(model.boundary_apertures(side).unwrap().aperture(3).unwrap()
                .relaxation().unwrap().dissipated_energy_j() > 0.0);
        }
    }
}

#[test]
fn downstream_refusal_cancellation_and_horizon_resume_preserve_all_source_and_key_states() {
    let plate = plate();
    let gate = CancelGate::new_clock_free();
    for family in FAMILIES {
        let mut actual = model(&plate, family, 64);
        let mut reference = model(&plate, family, 96);
        let keys = key_drives(256);
        for _ in 0..32 {
            assert_eq!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap(),
                reference.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap());
        }
        let energy = actual.stored_energy_j();
        let frames: Vec<_> = SIDES.into_iter().flat_map(|side|
            (0..4).map(move |node| (side, node))).map(|(side, node)|
            (side, node, *actual.node_frame(side, node).unwrap())).collect();
        let keys_before: Vec<_> = SIDES.into_iter().map(|side| {
            let key = actual.boundary_apertures(side).unwrap().aperture(3).unwrap();
            (key.state(), key.relaxation().unwrap().memory_sqrt_j().to_vec())
        }).collect();
        // The upstream valve/key/network preview succeeds before these downstream
        // failures. None of its candidate states may become an accepted state.
        for bad in [vec![], vec![BoundaryDrive { node: 2, ..keys[1][0] }],
            vec![BoundaryDrive { opening_force_n: f64::NAN, ..keys[1][0] }]] {
            assert!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &bad], &gate).is_err());
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &cancelled).is_err());
        assert!(actual.step_under_gate(drive(40.0), &gate).is_err());
        assert_eq!(actual.stored_energy_j(), energy);
        assert_eq!(actual.valve().accepted_steps(), 32);
        for (side, node, frame) in frames { assert_eq!(actual.node_frame(side, node), Some(&frame)); }
        for (index, side) in SIDES.into_iter().enumerate() {
            let key = actual.boundary_apertures(side).unwrap().aperture(3).unwrap();
            assert_eq!(key.state(), keys_before[index].0);
            assert_eq!(key.relaxation().unwrap().memory_sqrt_j(), keys_before[index].1);
            assert_eq!(key.accepted_steps(), 32);
        }
        for _ in 32..64 {
            assert_eq!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap(),
                reference.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap());
        }
        assert!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).is_err());
        assert!(actual.extend_step_budget(64).is_err());
        actual.extend_step_budget(96).unwrap();
        for side in SIDES {
            assert_eq!(actual.boundary_apertures(side).unwrap().aperture(3).unwrap().spec().max_steps, 96);
        }
        for _ in 64..96 {
            assert_eq!(actual.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap(),
                reference.step_with_boundaries_under_gate(drive(40.0), [&keys[0], &keys[1]], &gate).unwrap());
        }
    }
}

#[test]
fn material_key_binding_rejects_stale_regional_transport_and_lifetime() {
    let plate = plate();
    for mismatch in [false, true] {
        let cold = unloaded(Family::Modal, 64);
        let actual = cold.section_gases(ValveSide::Downstream)[2];
        let z = cold.represented_sections(ValveSide::Downstream)[2].impedance_pa_s_m3;
        let mut key = boundary(&plate, actual, z, if mismatch { 64 } else { 63 });
        if mismatch {
            let mut stale = actual;
            stale.dynamic_viscosity *= 2.0;
            stale.prandtl *= 2.0;
            stale.validate_coherence().unwrap();
            key.aperture.set_viscous_slit(ViscousSlit::from_gas(&stale, 0.0003).unwrap()).unwrap();
        }
        assert!(cold.with_boundary_apertures(ValveSide::Downstream, vec![key], 1,
            &CancelGate::new_clock_free()).is_err());
    }
}
