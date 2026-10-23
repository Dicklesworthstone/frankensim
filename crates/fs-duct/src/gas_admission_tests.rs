//! Every gas consumer admits the same material-owned physical state.
use crate::dg::{DgBoundary, DgDuct, DgMesh};
use crate::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary,
    network_port_impedance,
};
use crate::*;
use fs_material::gas::{GasSpec, GasState};

fn air(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn duct() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.01,
            length: 0.5,
        }],
    }
}

fn network() -> DuctNetwork {
    DuctNetwork {
        nodes: vec![
            NetworkNode {
                id: "inlet".into(),
                boundary: NodeBoundary::Port {
                    label: "inlet".into(),
                },
            },
            NetworkNode {
                id: "outlet".into(),
                boundary: NodeBoundary::Termination(NetworkLoad::PressureRelease),
            },
        ],
        edges: vec![NetworkEdge {
            id: "tube".into(),
            from: 0,
            to: 1,
            duct: duct(),
            corrections: JunctionCorrections::Disabled,
        }],
    }
}

fn dg(gas: &GasState) -> Result<DgDuct, DuctError> {
    DgDuct::new(
        &duct(),
        gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
}

#[test]
fn stale_thermodynamic_snapshots_refuse_across_frequency_time_and_radiation() {
    let mutations: [fn(&mut GasState); 6] = [
        |gas| gas.pressure *= 2.0,
        |gas| gas.temperature *= 1.1,
        |gas| gas.gamma = 1.3,
        |gas| gas.density *= 2.0,
        |gas| {
            gas.sound_speed *= 1.1;
            gas.characteristic_impedance = gas.density * gas.sound_speed;
        },
        |gas| {
            gas.specific_heat_cp *= 1.1;
            gas.prandtl = gas.dynamic_viscosity * gas.specific_heat_cp / gas.thermal_conductivity;
        },
    ];
    let omega = std::f64::consts::TAU * 233.0;
    for mutate in mutations {
        let mut gas = air(101_325.0);
        mutate(&mut gas);
        assert!(segment_wave(&gas, 0.01, omega, LossModel::Lossless).is_err());
        assert!(
            input_impedance(
                &duct(),
                &gas,
                omega,
                LossModel::Lossless,
                Termination::IdealOpen
            )
            .is_err()
        );
        assert!(
            input_impedance_load(&duct(), &gas, omega, LossModel::Lossless, C64::ZERO, None)
                .is_err()
        );
        assert!(network_port_impedance(&network(), &gas, omega, LossModel::Lossless).is_err());
        assert!(dg(&gas).is_err());
        assert!(
            compact_radiation_impedance(Termination::UnflangedOpen, &gas, 0.01, omega).is_err()
        );
        assert!(flanged_pipe_impedance(&gas, 0.01, omega, 1).is_err());
        let (_, witness) =
            fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.0).unwrap();
        assert!(absorbed_spherical_pressure(1.0, &gas, &witness, omega, 1.0, 0.0).is_err());
    }
}

#[test]
fn reevaluated_pressure_changes_impedance_and_viscous_loss_without_spurious_detuning() {
    let low = air(101_325.0);
    let high = air(202_650.0);
    let omega = std::f64::consts::TAU * 233.0;
    // At fixed temperature and composition, ideal-gas sound speed does not
    // depend on pressure. Inertia, compressibility and viscous loss still do.
    assert_eq!(low.sound_speed, high.sound_speed);
    assert_eq!(2.0 * low.density, high.density);
    let a = segment_wave(&low, 0.01, omega, LossModel::Lossless).unwrap();
    let b = segment_wave(&high, 0.01, omega, LossModel::Lossless).unwrap();
    assert_eq!(a.wavenumber, b.wavenumber);
    assert_eq!(
        a.characteristic_impedance.scale(2.0),
        b.characteristic_impedance
    );
    let za = input_impedance(
        &duct(),
        &low,
        omega,
        LossModel::Lossless,
        Termination::IdealOpen,
    )
    .unwrap()
    .impedance;
    let zb = network_port_impedance(&network(), &high, omega, LossModel::Lossless)
        .unwrap()
        .z[0];
    assert!((zb - za.scale(2.0)).abs() <= 1.0e-10 * zb.abs());
    let a = segment_wave(&low, 0.01, omega, LossModel::WideTube).unwrap();
    let b = segment_wave(&high, 0.01, omega, LossModel::WideTube).unwrap();
    assert!(b.wavenumber.im > 0.0 && b.wavenumber.im < a.wavenumber.im);
    assert_eq!(dg(&low).unwrap().gas(), &low);
    assert_eq!(dg(&high).unwrap().gas(), &high);
}
