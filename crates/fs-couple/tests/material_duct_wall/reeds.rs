//! The same material wall owner in both actual reed feedback networks.
use super::*;
use fs_couple::bernoulli_aperture::{BernoulliAperture,
    dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec},
    paired::{ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits,
        PairedApertureSpec},
    fractional::FractionalTubePolicy,
    two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive,
        TwoSidedLimits, ValveSide, combined_inlet_state},
    two_sided::observer::TwoSidedPressure,
    viscothermal::ViscothermalSectionSpec,
    wall::material::{MaterialWallInput, WallRelaxation}};
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_dcontact::Obstacle;
use fs_material::visco::GeneralizedMaxwell;

const DT: f64 = 1.0 / 48_000.0;
const SIDES: [ValveSide; 2] = [ValveSide::Upstream, ValveSide::Downstream];
fn sources() -> [ApertureDuctSpec; 2] {
    std::array::from_fn(|_| ApertureDuctSpec { graph: graph(), gases: vec![gas(); 2], viscothermal: vec![] })
}
fn input(young: f64) -> MaterialWallInput {
    MaterialWallInput { assignments: vec![assignment(1000.0, young, 2.0)], band_hz: BAND, max_walls: 1 }
}
fn memory(scale: f64) -> MaterialWallInput {
    let mut w = input(1e6);
    w.assignments[0].loss = WallMaterialLoss::Maxwell(WallRelaxation::new(
        GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003 * scale), (4e5, 0.0005 * scale)]).unwrap(),
        "synthetic fixed-temperature duct wall; not measured wood or cane".into(), BAND, 1.0, 0.8).unwrap());
    w
}
fn policy() -> [FractionalTubePolicy; 2] {
    [FractionalTubePolicy { maximum_frequency_hz: 1000.0, maximum_summed_phase_error_rad: 0.1 }; 2]
}
fn limits() -> TwoSidedLimits {
    TwoSidedLimits { maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16,
        energy_relative_tolerance: 1e-7, flow_absolute_tolerance_m3_s: 1e-16, flow_relative_tolerance: 1e-8 }
}
fn valve(sources: &[ApertureDuctSpec; 2], paired: bool) -> LoadedAperture {
    let (g, z) = combined_inlet_state(sources).unwrap();
    let b = ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 500.0,
        damping_n_s_m: 0.03, pressure_area_m2: 1.25e-4 };
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic material-duct reed contact".into()).unwrap();
    if paired {
        let mut blades = [b; 2]; blades[1].mass_kg *= 1.7;
        // Both images sit outside the midpoint one-cent band (DEF-04): exact blocks.
        PairedAperture::new_with_propagator(PairedApertureSpec {
            rest_opening_m: 4e-4, width_m: 0.013, blades, density_kg_m3: g.density,
            impedance_pa_s_m3: z, time_step_s: DT, max_steps: 2048,
            limits: PairedApertureLimits { maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0, max_evaluations: 256,
                velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
        }, [BladeState::default(); 2], contact, BladePropagator::ExactZoh).unwrap().into()
    } else {
        DynamicAperture::new_with_propagator(DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013,
                closing_pressure_pa: 0.5 * b.stiffness_n_m * 4e-4 / b.pressure_area_m2 },
            mass_kg: 0.5 * b.mass_kg, stiffness_n_m: 0.5 * b.stiffness_n_m,
            damping_ratio: b.damping_n_s_m / (2.0 * (b.mass_kg * b.stiffness_n_m).sqrt()),
            density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: 2048,
        }, ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 }, contact,
            BladePropagator::ExactZoh).unwrap().into()
    }
}
fn model(s: [ApertureDuctSpec; 2], w: [Option<MaterialWallInput>; 2], paired: bool,
    fractional: bool) -> TwoSidedAperture {
    TwoSidedAperture::with_material_walls(valve(&s, paired), s, w, limits(),
        fractional.then(policy), &CancelGate::new_clock_free()).unwrap()
}
fn drive(p: f64) -> TwoSidedDrive { TwoSidedDrive { source_pressure_pa: p, blade_forces_n: [0.0; 2] } }
fn snapshot(m: &TwoSidedAperture) -> Vec<f64> {
    let mut result = vec![m.stored_energy_j(), m.valve().accepted_steps() as f64];
    for side in SIDES {
        let w = m.wall_state(side, 1).unwrap();
        let n = m.node_frame(side, 1).unwrap();
        result.extend([w.displacement_m, w.velocity_m_s, w.stored_energy_j,
            m.wall_relaxation_energy_j(side, 1).unwrap(), n.pressure_pa, n.load_flow_m3_s,
            n.stored_energy_j, n.absorbed_energy_j]);
    }
    result
}

#[test]
fn absent_and_empty_wall_selections_retain_both_legacy_images_exactly() {
    let gate = CancelGate::new_clock_free();
    for paired in [false, true] {
        for fractional in [false, true] {
            let s = sources();
            let mut original = if fractional {
                TwoSidedAperture::with_fractional_transits(valve(&s, paired), s, limits(), policy(), &gate).unwrap()
            } else { TwoSidedAperture::new(valve(&s, paired), s, limits(), &gate).unwrap() };
            let mut absent = model(sources(), [None, None], paired, fractional);
            let empty = MaterialWallInput { assignments: vec![], band_hz: BAND, max_walls: 0 };
            let mut selected = model(sources(), [Some(empty.clone()), None], paired, fractional);
            assert_eq!(selected.wall_input(ValveSide::Upstream), Some(&empty));
            assert!(absent.wall_input(ValveSide::Upstream).is_none());
            assert!(selected.material_walls(ValveSide::Upstream).is_empty());
            assert!(selected.wall_state(ValveSide::Upstream, 1).is_err());
            for k in 0..128 {
                let d = drive(if k < 64 { 50.0 } else { 0.0 });
                let expected = original.step_under_gate(d, &gate).unwrap();
                assert_eq!(absent.step_under_gate(d, &gate).unwrap(), expected);
                assert_eq!(selected.step_under_gate(d, &gate).unwrap(), expected);
            }
        }
    }
}

#[test]
fn either_duct_material_changes_single_and_independent_blade_motion_and_actual_outlet_flow() {
    let gate = CancelGate::new_clock_free();
    for paired in [false, true] {
        for changed_side in 0..2 {
            let base = [Some(input(1e6)), Some(input(1e6))];
            let mut changed = base.clone(); changed[changed_side] = Some(input(4e6));
            let mut a = model(sources(), base.clone(), paired, false);
            let mut b = model(sources(), changed.clone(), paired, false);
            let delay = a.represented_sections(SIDES[changed_side])[0].one_way_samples;
            let mut changed_motion = false; let mut changed_flow = false;
            for side in SIDES {
                assert_eq!(a.source(side).graph, graph());
                assert_eq!(a.wall_input(side), base[side as usize].as_ref());
                assert_eq!(a.material_walls(side)[0].length_m(), 0.15);
                assert_eq!(a.graph(side).nodes[1], a.material_walls(side)[0].shunt());
            }
            for k in 0..512 {
                let d = drive(if k < 128 { 50.0 } else { 0.0 });
                let x = a.step_under_gate(d, &gate).unwrap();
                let y = b.step_under_gate(d, &gate).unwrap();
                assert!(x.balance_residual_j().abs() < 1e-12);
                assert!(y.balance_residual_j().abs() < 1e-12);
                if k < 2 * delay { assert_eq!(x.aperture.state, y.aperture.state); }
                changed_motion |= x.aperture.state != y.aperture.state;
                changed_flow |= a.node_frame(ValveSide::Downstream, 2).unwrap().net_flow_into_node_m3_s
                    != b.node_frame(ValveSide::Downstream, 2).unwrap().net_flow_into_node_m3_s;
            }
            assert!(changed_motion && changed_flow, "wall data must react back, not just label output");
            for side in SIDES {
                assert!(a.wall_state(side, 1).unwrap().stored_energy_j > 0.0);
                close(a.wall_state(side, 1).unwrap().stored_energy_j, a.node_frame(side, 1).unwrap().stored_energy_j);
            }
        }
    }
}

#[test]
fn independent_wall_memories_share_the_energy_window_and_refusal_rollback() {
    let gate = CancelGate::new_clock_free();
    let mut a = model(sources(), [Some(memory(1.0)), Some(memory(3.0))], true, true);
    let mut work = 0.0; let mut loss = 0.0;
    for k in 0..512 {
        let f = a.step_under_gate(drive(if k < 128 { 50.0 } else { 0.0 }), &gate).unwrap();
        work += f.source_work_j; loss += f.dissipated_energy_j;
        assert!(f.balance_residual_j().abs() < 1e-12);
        for side in SIDES {
            close(a.wall_state(side, 1).unwrap().stored_energy_j + a.wall_relaxation_energy_j(side, 1).unwrap(),
                a.node_frame(side, 1).unwrap().stored_energy_j);
        }
    }
    assert!(loss > 0.0);
    assert!((a.stored_energy_j() + loss - work).abs() < 1e-11);
    for side in SIDES { assert!(a.wall_relaxation_energy_j(side, 1).unwrap() > 0.0); }
    let old = snapshot(&a);
    assert!(a.step_under_gate(drive(f64::NAN), &gate).is_err());
    assert!(a.set_terminal_reflection(ValveSide::Upstream, 1, 0.0).is_err());
    assert_eq!(snapshot(&a), old);
    gate.request();
    assert!(a.step_under_gate(drive(50.0), &gate).is_err());
    assert_eq!(snapshot(&a), old);
    assert!(a.wall_state(ValveSide::Upstream, 0).is_err());
}

#[test]
fn regional_gas_and_fractional_loss_keep_the_original_wall_geometry_and_source_records() {
    let mut s = sources();
    for side in &mut s {
        side.graph.sections[0].length_m = 0.1501;
        side.graph.sections[1].length_m = 0.1502;
        side.gases[1] = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, gas().pressure).unwrap();
        side.viscothermal.push(ViscothermalSectionSpec { section: 0, cells: 4,
            minimum_frequency_hz: 300.0, maximum_frequency_hz: 500.0 });
    }
    let original = s.clone();
    let w = [Some(memory(1.0)), Some(input(4e6))];
    let mut a = model(s, w.clone(), true, true);
    for (i, side) in SIDES.into_iter().enumerate() {
        assert_eq!(a.source(side).graph, original[i].graph);
        assert_eq!(a.source(side).gases, original[i].gases);
        assert_eq!(a.source(side).viscothermal, original[i].viscothermal);
        assert_eq!(a.wall_input(side), w[i].as_ref());
        close(a.material_walls(side)[0].length_m(), 0.15015);
        assert_eq!(a.graph(side).nodes[1], a.material_walls(side)[0].shunt());
        let r = &a.viscothermal_sections(side)[0];
        close(r.represented_length_m, 0.1501);
        assert!(a.fractional_sections(side).is_some());
        for gas in &a.section_gases(side)[r.section_range[0]..r.section_range[1]] { assert_eq!(*gas, original[i].gases[0]); }
        assert_eq!(*a.section_gases(side).last().unwrap(), original[i].gases[1]);
    }
    let gate = CancelGate::new_clock_free();
    for _ in 0..256 { a.step_under_gate(drive(50.0), &gate).unwrap(); }
    for side in SIDES { assert!(a.wall_state(side, 1).unwrap().stored_energy_j > 0.0); }
}

#[test]
fn original_pressure_renderer_retains_both_materials_and_matches_direct_flow_across_callbacks() {
    let receiver = CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 4,
        angular_points: 16, maximum_frequency_hz: 1000.0 };
    let mut signals = Vec::new();
    for scale in [1.0, 3.0] {
        let w = [Some(memory(scale)), Some(input(1e6))];
        let mut direct = model(sources(), w.clone(), true, true);
        let mut observer = BaffledPressure::circular_outlet(0.006, RATE, receiver,
            RayleighMedium { density: gas().density, sound_speed: gas().sound_speed }).unwrap();
        let gate = CancelGate::new_clock_free();
        let mut expected = Vec::new();
        for _ in 0..512 {
            direct.step_under_gate(drive(50.0), &gate).unwrap();
            expected.push(observer.step(&[direct.node_frame(ValveSide::Downstream, 2).unwrap().net_flow_into_node_m3_s]).unwrap());
        }
        for partition in [1, 37, 512] {
            let mut p = TwoSidedPressure::new(model(sources(), w.clone(), true, true),
                ValveSide::Downstream, 2, receiver, RATE, 512).unwrap();
            p.set_drive(drive(50.0)).unwrap();
            let mut actual = [0.0; 512];
            for block in actual.chunks_mut(partition) { p.block(block).unwrap(); }
            assert_eq!(actual.as_slice(), expected.as_slice());
            assert_eq!(snapshot(p.system()), snapshot(&direct));
            assert_eq!(p.system().wall_input(ValveSide::Upstream), w[0].as_ref());
        }
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        signals.push(expected);
    }
    assert_ne!(signals[0], signals[1]);
}

#[test]
fn a_bad_second_material_or_unresolved_wall_never_becomes_a_partly_coupled_source() {
    let gate = CancelGate::new_clock_free();
    let mut invalids = Vec::new();
    let mut w = input(1e6); w.assignments[0].node = 2; invalids.push(w);
    let mut w = input(1e6); w.max_walls = 0; invalids.push(w);
    let mut w = input(1e6); w.band_hz[1] = f64::from(RATE) / 2.0; invalids.push(w);
    let mut w = input(1e6); w.assignments[0].material = state(&properties(1000.0, 1e6, 2.0)[..3], None, None); invalids.push(w);
    invalids.push(memory(1e-6)); // Original memory exists but its time scale cannot be resolved.
    for w in invalids {
        let s = sources();
        assert!(TwoSidedAperture::with_material_walls(valve(&s, true), s,
            [Some(input(1e6)), Some(w)], limits(), Some(policy()), &gate).is_err());
    }
    gate.request(); let s = sources();
    assert!(TwoSidedAperture::with_material_walls(valve(&s, true), s,
        [Some(input(1e6)), Some(input(1e6))], limits(), None, &gate).is_err());
}
