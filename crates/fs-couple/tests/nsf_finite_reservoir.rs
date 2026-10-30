//! Finite source ownership through the original NASA9, NSF and DG solvers.
//! Synthetic material and coarse grids exercise composition, not flute fidelity.
use fs_couple::dg_port::DgPortLoad;
use fs_couple::nsf_dg::{NsfDg, NsfDgDrive, NsfDgLimits, NsfFiniteReservoirLimits};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfBoundaryDrive, NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGrid3,
    NsfOpenState, NsfReservoirAperture, NsfSpec, ThermochemicalNsfGas, WallThermal,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[path = "../../fs-lbm/tests/nsf3_thermochemical/source.rs"]
mod source;

const RATE: f64 = 96_000.0;
const CELLS: usize = 128;
const VOLUME: f64 = 1e-5;

fn source_limits() -> NsfFiniteReservoirLimits {
    NsfFiniteReservoirLimits {
        maximum_relative_mass_exchange: 0.1,
        maximum_relative_energy_exchange: 0.1,
        mass_absolute_tolerance_kg: 1e-18,
        energy_absolute_tolerance_j: 1e-11,
        relative_tolerance: 1e-9,
    }
}

fn bare(pressure: f64, upstream: bool, substeps: usize) -> NsfDg<D3q27> {
    let gate = CancelGate::new_clock_free();
    let model = source::model(-4000.0, 200_000.0, 1.0);
    let medium = ThermochemicalNsfGas::new(
        model.clone(), source::point(300.0, 100_000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    ).unwrap();
    let gas = *model.evaluate(medium.context()).unwrap().state();
    let initial = *model.evaluate(&source::point(300.0, pressure)).unwrap().state();
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt()
        / RATE / substeps as f64;
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: gas.pressure,
        reference_temperature_k: gas.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    let zone = NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims: [8, 4, 4], periodic: [false; 3], gas: medium,
            temperature_k: gas.temperature, pressure_pa: gas.pressure, dx_m: dx,
            temperature_window_k: [250.0, 600.0], maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored finite source fixture; not a gas validation",
            maximum_cells: CELLS, wall: WallThermal::Adiabatic,
            workers: 1, transport_points: 17,
        },
        |_| false, |_| (1.0, [0.0; 3], gas.temperature),
    ).unwrap().with_characteristic_axis(NsfCharacteristicAxis {
        axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
    }).unwrap().with_reservoir_aperture(NsfReservoirAperture {
        axis: 1, high: false,
        cells: vec![[3, 0, 1], [3, 0, 2], [4, 0, 1], [4, 0, 2]],
        state: NsfOpenState {
            density_kg_m3: initial.density, velocity_m_s: [0.0; 3],
            temperature_k: initial.temperature,
        },
    }).unwrap();
    let load = |outlet| {
        let duct = DgDuct::new(
            &Duct { segments: vec![Segment::Cylinder {
                radius: (16.0 * dx * dx / core::f64::consts::PI).sqrt(), length: 0.03,
            }] },
            &gas, DgMesh { elements: 2, order: 2 },
            DgBoundary::Port(PortId(0)), outlet,
        ).unwrap();
        DgPortLoad::new(duct, 1.0 / RATE, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
    };
    let owner = NsfDg::new(zone, load(DgBoundary::Matched), true, NsfDgLimits {
        substeps, maximum_cells: CELLS, maximum_cell_steps: CELLS * substeps * 128,
        maximum_evaluations: 128, maximum_steps: 128, incoming_bracket_m_s: [-0.5, 0.5],
        pressure_absolute_tolerance_pa: 1e-7, pressure_relative_tolerance: 1e-9,
        // Loose projection admission isolates source conservation and ownership.
        projection_absolute_tolerance_j: 0.0, projection_relative_tolerance: 1.0,
    }, &gate).unwrap();
    if upstream { owner.with_upstream_load(load(DgBoundary::Closed), &gate).unwrap() }
    else { owner }
}

macro_rules! snapshot {
    ($owner:expr) => {{
        let owner = &$owner;
        (owner.zone().state_hash(), owner.zone().cumulative_ledger(),
            owner.zone().reservoir_aperture_at(0).unwrap().clone(),
            owner.zone().last_reservoir_exchanges().to_vec(),
            owner.load().state().clone(), owner.upstream_load().map(|load| load.state().clone()),
            owner.finite_reservoir().cloned(), owner.accepted_steps())
    }};
}

#[test]
fn finite_source_reuses_original_fluxes_and_feeds_accepted_state_into_next_interval() {
    let gate = CancelGate::new_clock_free();
    for upstream in [false, true] {
        for substeps in [1, 3] {
            for pressure in [99_998.0, 100_002.0] {
                let mut finite = bare(pressure, upstream, substeps)
                    .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
                let initial = *finite.finite_reservoir().unwrap().state().state();
                let mut prescribed = bare(pressure, upstream, substeps);
                let original = prescribed.step_with_drive(NsfDgDrive::ReservoirAperture {
                    aperture: 0, pressure_pa: initial.pressure, temperature_k: initial.temperature,
                }, &gate).unwrap();
                let first = finite.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
                assert_eq!(finite.zone().state_hash(), prescribed.zone().state_hash());
                assert_eq!(finite.load().state(), prescribed.load().state());
                assert_eq!(finite.upstream_load().map(|load| load.state()),
                    prescribed.upstream_load().map(|load| load.state()));
                assert_eq!(first.reservoir_exchanges, original.reservoir_exchanges);
                let account = first.finite_reservoir.as_ref().unwrap();
                let exchange = first.reservoir_exchanges[0];
                assert_eq!(account.source_exchange, exchange);
                assert_eq!(account.aperture, 0);
                assert_eq!(account.volume_m3, VOLUME);
                assert!(account.mass_change_kg * (100_000.0 - pressure) > 0.0);
                assert!((account.mass_change_kg - exchange.mass_out_kg).abs() <= 1e-18);
                assert!((account.sensible_energy_change_j - exchange.total_energy_out_j).abs() <= 1e-9);
                // The finite source receives the conservative FV energy, not the
                // much smaller acoustic pressure-work projection.
                assert!((exchange.total_energy_out_j - exchange.acoustic_pressure_work_out_j).abs() > 1e-11);
                assert!(account.mass_residual_kg.abs() <= 1e-18);
                assert!(account.energy_residual_j.abs() <= 1e-9);
                assert!(exchange.mass_exchange_variation_kg >= exchange.mass_out_kg.abs());
                assert!(exchange.energy_exchange_variation_j >= exchange.total_energy_out_j.abs());
                let before_second = *finite.finite_reservoir().unwrap().state().state();
                assert_eq!(before_second.pressure, account.final_pressure_pa);
                assert_eq!(before_second.temperature, account.final_temperature_k);
                assert_ne!(before_second.pressure, initial.pressure);
                let second = finite.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
                let next = second.finite_reservoir.as_ref().unwrap();
                assert_eq!(next.initial_pressure_pa, account.final_pressure_pa);
                assert_eq!(next.initial_temperature_k, account.final_temperature_k);
                let held = finite.zone().reservoir_aperture_at(0).unwrap().state;
                assert_eq!(held.temperature_k, before_second.temperature);
                assert!((held.density_kg_m3 - before_second.density).abs()
                    <= 32.0 * f64::EPSILON * before_second.density);
                assert_eq!(finite.accepted_steps(), 2);
            }
        }
    }
}

#[test]
fn quiescent_source_does_not_create_heat_or_consume_inventory() {
    let gate = CancelGate::new_clock_free();
    let mut owner = bare(100_000.0, true, 3)
        .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
    let original = owner.finite_reservoir().unwrap().clone();
    for _ in 0..3 {
        let frame = owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
        let account = frame.finite_reservoir.unwrap();
        assert!(account.mass_change_kg.abs() < 1e-16);
        assert!(account.sensible_energy_change_j.abs() < 1e-10);
    }
    let final_state = owner.finite_reservoir().unwrap();
    assert!((final_state.mass_kg() - original.mass_kg()).abs() < 1e-16);
    assert!((final_state.state().state().temperature - 300.0).abs() < 1e-8);
}

#[test]
fn refusal_and_cancellation_preserve_finite_gas_fluid_and_both_ducts() {
    let gate = CancelGate::new_clock_free();
    let mut owner = bare(100_002.0, true, 1)
        .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
    owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    let accepted = snapshot!(owner);
    for drive in [
        NsfDgDrive::Boundary(NsfBoundaryDrive::IncomingCharacteristic(0.0)),
        NsfDgDrive::UpstreamPressureRise(0.0),
        NsfDgDrive::Reservoir { pressure_pa: 100_000.0, temperature_k: 300.0 },
        NsfDgDrive::ReservoirAperture { aperture: 0, pressure_pa: 100_000.0, temperature_k: 300.0 },
        NsfDgDrive::FiniteReservoir { aperture: 1 },
    ] {
        assert!(owner.step_with_drive(drive, &gate).is_err());
        assert_eq!(snapshot!(owner), accepted);
    }
    // Public convenience entry points must enforce the same source ownership.
    assert!(owner.step_with_reservoir_aperture(0, 100_000.0, 300.0, &gate).is_err());
    assert_eq!(snapshot!(owner), accepted);
    assert!(owner.step_with_upstream_pressure_rise(0.0, &gate).is_err());
    assert_eq!(snapshot!(owner), accepted);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &cancelled).is_err());
    assert_eq!(snapshot!(owner), accepted);
    owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    assert_eq!(owner.accepted_steps(), 2);
}

#[test]
fn finite_inventory_limits_refuse_without_publishing_a_partial_interval() {
    let gate = CancelGate::new_clock_free();
    for upstream in [false, true] {
        let mut owner = bare(100_002.0, upstream, 3)
            .with_finite_reservoir(0, 1e-18, source_limits(), &gate).unwrap();
        let before = snapshot!(owner);
        assert!(owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).is_err());
        assert_eq!(snapshot!(owner), before);
    }
    for volume in [0.0, -VOLUME, f64::NAN, f64::INFINITY] {
        assert!(bare(100_002.0, false, 1)
            .with_finite_reservoir(0, volume, source_limits(), &gate).is_err());
    }
}

#[path = "nsf_finite_reservoir/actuation.rs"]
mod actuation;
