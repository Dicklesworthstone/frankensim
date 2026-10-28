//! Physical-owner composition checks, not a resolved-instrument qualification.
use super::performance::NsfDgPerformance;
use super::*;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::observation::PressureRenderer;
use fs_duct::dg::PortId;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGasMedium, NsfSpec, ThermochemicalNsfGas,
    WallThermal,
};
use fs_material::gas::GasSpec;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use std::sync::Arc;
#[path = "../../../fs-lbm/tests/nsf3_thermochemical/source.rs"]
mod source;

const RATE: u32 = 96000;
const CELLS: usize = 32;

fn limits(substeps: usize) -> NsfDgLimits {
    NsfDgLimits {
        substeps,
        maximum_cells: CELLS,
        maximum_cell_steps: CELLS * substeps * 64,
        maximum_evaluations: 64,
        maximum_steps: 128,
        incoming_bracket_m_s: [-0.05, 0.05],
        pressure_absolute_tolerance_pa: 1e-7,
        pressure_relative_tolerance: 1e-9,
        // These tests audit the original accounts and transaction. The separate
        // zero-allowance refusal below exercises projection admission. A loose
        // work budget here is deliberately NOT a reflection/accuracy claim.
        projection_absolute_tolerance_j: 0.0,
        projection_relative_tolerance: 1.0,
    }
}
fn zone<G: Into<NsfGasMedium>>(
    gas: G,
    reference: GasState,
    high: bool,
    substeps: usize,
) -> NsfGrid3<D3q27> {
    let velocity = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
    let c = (reference.gamma / 3.0).sqrt();
    let sign = if high { 1.0 } else { -1.0 };
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: reference.pressure,
        reference_temperature_k: reference.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims: [CELLS, 1, 1],
            periodic: [false, true, true],
            gas,
            temperature_k: reference.temperature,
            pressure_pa: reference.pressure,
            dx_m: velocity / f64::from(RATE) / substeps as f64,
            temperature_window_k: [250.0, 600.0],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored numerical coupling fixture; not gas validation",
            maximum_cells: CELLS,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |x| {
            let distance = if high { CELLS - 1 - x[0] } else { x[0] };
            let z = (distance as f64 - 5.0) / 8.0;
            let pulse = 1e-5 * fs_math::det::exp(-z * z);
            (
                1.0 + pulse / reference.gamma,
                [sign * c * pulse / reference.gamma, 0.0, 0.0],
                reference.temperature * (1.0 + (reference.gamma - 1.0) / reference.gamma * pulse),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(NsfCharacteristicAxis {
        axis: 0,
        faces: [face; 2],
        maximum_relative_perturbation: 0.05,
    })
    .unwrap()
}
fn load(
    zone: &NsfGrid3<D3q27>,
    gas: GasState,
    area_factor: f64,
    clock_factor: f64,
    outlet: DgBoundary,
) -> DgPortLoad {
    let radius = (zone.law().dx_m.powi(2) * area_factor / core::f64::consts::PI).sqrt();
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.03,
            }],
        },
        &gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        outlet,
    )
    .unwrap();
    DgPortLoad::new(
        duct,
        clock_factor / f64::from(RATE),
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
    )
    .unwrap()
}
fn system(pressure: f64, high: bool, substeps: usize, outlet: DgBoundary) -> NsfDg<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, pressure).unwrap();
    let zone = zone(gas, reference, high, substeps);
    let load = load(&zone, reference, 1.0, 1.0, outlet);
    NsfDg::new(
        zone,
        load,
        high,
        limits(substeps),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn both_oriented_faces_drive_original_dg_and_keep_extensive_flux_separate() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for substeps in [1, 3] {
            let mut system = system(100000.0, high, substeps, DgBoundary::Matched);
            let mass = system.zone().total_mass_kg();
            let energy = system.zone().total_energy_j();
            let frame = system.step(0.0, &gate).unwrap();
            let face = frame.faces[usize::from(high)];
            assert!(frame.inlet_pressure_pa > 0.0 && frame.inlet_flow_m3_s > 0.0);
            assert_eq!(system.zone().steps(), substeps);
            assert_eq!(system.load().state().step_count(), 1);
            assert!(frame.pressure_residual_pa.abs() < 1e-6);
            let dm = frame.faces.iter().map(|f| f.mass_out_kg).sum::<f64>();
            let de = frame
                .faces
                .iter()
                .map(|f| f.total_energy_out_j)
                .sum::<f64>();
            assert!((frame.flow_mass_change_kg + dm).abs() < 3e-12 * mass);
            assert!((frame.flow_energy_change_j + de).abs() < 3e-12 * energy);
            assert!((frame.flow.open_boundary_energy_j + de).abs() < 3e-12 * energy);
            let dt = system.load().propagator().dt();
            let identity = frame.projection_covariance_j
                + frame.pressure_residual_pa * (frame.inlet_flow_m3_s * dt)
                - frame.port_work_residual_j;
            let scale = face.acoustic_pressure_work_out_j.abs()
                + (frame.inlet_pressure_pa * frame.inlet_flow_m3_s * dt).abs();
            assert!((frame.acoustic_work_residual_j - identity).abs() < 1e-12 * scale);
            assert_eq!(
                frame.mass_projection_difference_kg,
                face.mass_out_kg - system.reference_gas().density * (frame.inlet_flow_m3_s * dt)
            );
        }
    }
}

#[test]
fn live_bore_reflection_and_material_pressure_change_the_same_resolved_source() {
    let gate = CancelGate::new_clock_free();
    let mut open = system(100000.0, true, 1, DgBoundary::Matched);
    let mut closed = system(100000.0, true, 1, DgBoundary::Closed);
    let mut dense = system(200000.0, true, 1, DgBoundary::Matched);
    let a = open.step(0.0, &gate).unwrap();
    let b = dense.step(0.0, &gate).unwrap();
    assert_eq!(
        dense.reference_gas().density,
        2.0 * open.reference_gas().density
    );
    assert!(b.inlet_pressure_pa > 1.5 * a.inlet_pressure_pa);
    closed.step(0.0, &gate).unwrap();
    for _ in 1..32 {
        open.step(0.0, &gate).unwrap();
        closed.step(0.0, &gate).unwrap();
    }
    assert_ne!(
        open.zone().state_hash(),
        closed.zone().state_hash(),
        "downstream radiation/reflection must feed back into original resolved populations"
    );
}

#[test]
fn finite_walled_geometry_changes_resolved_flow_and_the_same_live_bore() {
    let gate = CancelGate::new_clock_free();
    let build = |constricted: bool| {
        let gas = GasSpec::dry_air_ussa1976();
        let reference = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
        let velocity_unit = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: reference.pressure,
            reference_temperature_k: reference.temperature,
            reference_velocity_m_s: [0.0; 3],
            incoming_characteristic_m_s: 0.0,
        };
        let zone = NsfGrid3::<D3q27>::new(
            NsfSpec {
                dims: [12, 3, 3],
                periodic: [false; 3],
                gas,
                temperature_k: reference.temperature,
                pressure_pa: reference.pressure,
                dx_m: velocity_unit / f64::from(RATE),
                temperature_window_k: [250.0, 600.0],
                maximum_mach: 0.2,
                maximum_temperature_ratio: 1.0,
                bulk_to_shear_ratio: 0.73,
                bulk_ratio_source: "authored finite-passage composition fixture",
                maximum_cells: 108,
                wall: WallThermal::Adiabatic,
                workers: 1,
                transport_points: 17,
            },
            |x| constricted && x[0] == 5 && x[1] != 1,
            |_| (1.0, [0.0; 3], 300.0),
        )
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
        let load = load(&zone, reference, 9.0, 1.0, DgBoundary::Matched);
        let mut budget = limits(1);
        budget.maximum_cells = 108;
        budget.maximum_cell_steps = 108 * budget.maximum_evaluations;
        budget.incoming_bracket_m_s = [-0.5, 0.5];
        NsfDg::new(zone, load, true, budget, &gate).unwrap()
    };
    let mut empty = build(false);
    let mut aperture = build(true);
    let mut largest_pressure_change: f64 = 0.0;
    for _ in 0..24 {
        let a = empty.step(-0.1, &gate).unwrap();
        let b = aperture.step(-0.1, &gate).unwrap();
        largest_pressure_change =
            largest_pressure_change.max((a.inlet_pressure_pa - b.inlet_pressure_pa).abs());
        for frame in [&a, &b] {
            assert!(
                (frame.flow_mass_change_kg
                    + frame.faces.iter().map(|f| f.mass_out_kg).sum::<f64>())
                .abs()
                    < 1e-12 * empty.zone().total_mass_kg()
            );
            assert!(
                (frame.flow_energy_change_j
                    + frame
                        .faces
                        .iter()
                        .map(|f| f.total_energy_out_j)
                        .sum::<f64>())
                .abs()
                    < 1e-12 * empty.zone().total_energy_j()
            );
        }
    }
    assert!(
        largest_pressure_change > 1e-6,
        "solid geometry must reach the physical DG pressure"
    );
    assert_ne!(empty.load().state(), aperture.load().state());
    assert_eq!(aperture.zone().steps(), 24);
    assert_eq!(aperture.load().state().step_count(), 24);
}

#[test]
fn complete_gas_area_clock_and_trial_work_are_cold_admission_constraints() {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let changed = GasState::try_new(&gas, 300.0, 120000.0).unwrap();
    for (state, area, clock, maximum_cells, maximum_work) in [
        (changed, 1.0, 1.0, CELLS, CELLS * 64),
        (reference, 1.01, 1.0, CELLS, CELLS * 64),
        (reference, 1.0, 2.0, CELLS, CELLS * 64),
        (reference, 1.0, 1.0, CELLS - 1, CELLS * 64),
        (reference, 1.0, 1.0, CELLS, CELLS * 64 - 1),
    ] {
        let zone = zone(gas.clone(), reference, true, 1);
        let load = load(&zone, state, area, clock, DgBoundary::Matched);
        let mut limits = limits(1);
        limits.maximum_cells = maximum_cells;
        limits.maximum_cell_steps = maximum_work;
        assert!(NsfDg::new(zone, load, true, limits, &CancelGate::new_clock_free()).is_err());
    }
}

#[test]
fn failed_projection_bracket_and_cancellation_publish_neither_owner() {
    let gate = CancelGate::new_clock_free();
    let mut system = system(100000.0, true, 3, DgBoundary::Matched);
    let before = (
        system.zone().state_hash(),
        system.zone().cumulative_ledger(),
        *system.zone().characteristic_axis().unwrap(),
        system.load().state().clone(),
    );
    system.limits.projection_relative_tolerance = 0.0;
    assert!(matches!(
        system.step(0.0, &gate),
        Err(NsfDgError::Projection { .. })
    ));
    system.limits.projection_relative_tolerance = 1.0;
    system.limits.incoming_bracket_m_s = [1.0, 2.0];
    assert!(matches!(
        system.step(0.0, &gate),
        Err(NsfDgError::NoBracket { .. })
    ));
    system.limits.incoming_bracket_m_s = [-0.05, 0.05];
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.step(0.0, &cancelled),
        Err(NsfDgError::Cancelled)
    ));
    assert_eq!(system.accepted_steps(), 0);
    assert_eq!(
        before,
        (
            system.zone().state_hash(),
            system.zone().cumulative_ledger(),
            *system.zone().characteristic_axis().unwrap(),
            system.load().state().clone()
        )
    );
    system.step(0.0, &gate).unwrap();
    assert_eq!(system.accepted_steps(), 1);
}

#[test]
fn thermochemical_source_and_its_pressure_domain_survive_the_live_coupling() {
    let gas = ThermochemicalNsfGas::new(
        source::model(-4000.0, 100500.0, 1.0),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    let reference = *gas.model().evaluate(gas.context()).unwrap().state();
    let zone = zone(gas, reference, true, 1);
    let load = load(&zone, reference, 1.0, 1.0, DgBoundary::Matched);
    let gate = CancelGate::new_clock_free();
    let mut system = NsfDg::new(zone, load, true, limits(1), &gate).unwrap();
    system.step(0.0, &gate).unwrap();
    assert!(matches!(
        system.zone().medium(),
        NsfGasMedium::Thermochemical(_)
    ));
    assert!(
        system
            .zone()
            .thermochemical_state([CELLS - 1, 0, 0])
            .unwrap()
            .is_some()
    );
    let before = (
        system.zone().state_hash(),
        system.zone().cumulative_ledger(),
        *system.zone().characteristic_axis().unwrap(),
        system.load().state().clone(),
    );
    assert!(matches!(
        system.step(-20.0, &gate),
        Err(NsfDgError::Flow(_))
    ));
    assert_eq!(
        before,
        (
            system.zone().state_hash(),
            system.zone().cumulative_ledger(),
            *system.zone().characteristic_axis().unwrap(),
            system.load().state().clone()
        )
    );
    system.step(0.0, &gate).unwrap();
}

#[test]
fn actual_pressure_and_pcm_are_callback_partition_invariant() {
    let make = || {
        NsfDgPerformance::new(
            system(100000.0, true, 1, DgBoundary::Matched),
            ApertureObservation::Inlet,
            RATE,
            vec![0.0; 64],
            32,
            Arc::new(CancelGate::new_clock_free()),
        )
        .unwrap()
    };
    let mut direct = system(100000.0, true, 1, DgBoundary::Matched);
    let gate = CancelGate::new_clock_free();
    let expected: Vec<_> = (0..64)
        .map(|_| direct.step(0.0, &gate).unwrap().inlet_pressure_pa)
        .collect();
    let mut single = make();
    let mut partitioned = make();
    let mut a = vec![0.0; 64];
    let mut b = vec![0.0; 64];
    for sample in &mut a {
        single.block(std::slice::from_mut(sample)).unwrap();
    }
    let mut start = 0;
    for count in [16, 17, 31] {
        partitioned.block(&mut b[start..start + count]).unwrap();
        start += count;
    }
    assert_eq!(a, expected);
    assert_eq!(b, expected);
    assert!(a.iter().any(|p| p.abs() > 1e-4));
    assert_eq!(
        single.system().zone().state_hash(),
        direct.zone().state_hash()
    );
    assert_eq!(partitioned.system().load().state(), direct.load().state());
    assert_eq!(
        crate::pcm_wav::encode_pcm16_wav(&a, RATE, 100.0).unwrap(),
        crate::pcm_wav::encode_pcm16_wav(&b, RATE, 100.0).unwrap()
    );
    assert!(partitioned.block(&mut [0.0]).is_err());
}

#[test]
fn cartesian_face_area_admission_does_not_accumulate_resolution_roundoff() {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let dx: f64 = 5e-5;
    let area = (dx * dx) * 1024.0;
    let naive = (0..1024).map(|_| dx * dx).sum::<f64>();
    assert!((naive - area).abs() > 32.0 * f64::EPSILON * area);
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: reference.pressure,
        reference_temperature_k: reference.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    let zone = NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims: [4, 32, 32],
            periodic: [false, true, true],
            gas,
            temperature_k: 300.0,
            pressure_pa: 100000.0,
            dx_m: dx,
            temperature_window_k: [250.0, 600.0],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored Cartesian area admission fixture",
            maximum_cells: 4096,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap()
    .with_characteristic_axis(NsfCharacteristicAxis {
        axis: 0,
        faces: [face; 2],
        maximum_relative_perturbation: 0.05,
    })
    .unwrap();
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: (area / core::f64::consts::PI).sqrt(),
                length: 0.03,
            }],
        },
        &reference,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Matched,
    )
    .unwrap();
    let load = DgPortLoad::new(duct, zone.law().dt_s, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap();
    let mut limits = limits(1);
    limits.maximum_cells = 4096;
    limits.maximum_cell_steps = 4096 * limits.maximum_evaluations;
    let coupled = NsfDg::new(zone, load, true, limits, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(coupled.area_m2(), area);
}
