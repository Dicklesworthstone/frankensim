use super::*;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_material::porous::{PorousAmbient, PorousCard, gas::PorousGasLaw};
use fs_qty::semantic::{Composition, CompositionBasis};

const RADIUS: f64 = 0.02;
const LENGTH: f64 = 0.24;
const DEPTH: f64 = 0.05;

fn material() -> PorousGasLaw {
    // The existing Nguyen card explicitly gives k0=31.47e-10 m2 and
    // sigma=1.81e-5/k0. Preserve that characterization, not current gas drag.
    PorousGasLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}

fn area() -> f64 { core::f64::consts::PI * RADIUS * RADIUS }

fn porous(material: &PorousGasLaw) -> NetworkLoad {
    NetworkLoad::PorousRigid {
        material: material.clone(), thickness_m: DEPTH, area_m2: area(),
    }
}

fn two_ports() -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("input"), port("outlet")],
        edges: vec![edge("pipe", 0, 1, RADIUS, LENGTH)],
    }
}

fn reference_load(material: &PorousGasLaw, gas: &GasState, omega: f64) -> C64 {
    // Independent binding to the ORIGINAL material owner, followed by the
    // analytic pressure/volume-flow and phasor conversion.
    let mut card = material.reference_card().clone();
    card.flow_resistivity = gas.dynamic_viscosity / material.viscous_permeability_m2();
    let ambient = PorousAmbient {
        rho0: gas.density, c0: gas.sound_speed, eta: gas.dynamic_viscosity,
        gamma: gas.gamma, p0: gas.pressure, prandtl: gas.prandtl,
    };
    let surface = card.surface_impedance_rigid(ambient, omega, DEPTH).unwrap();
    C64::new(surface.re / area(), -surface.im / area())
}

#[test]
fn porous_both_network_load_paths_match_analytic_line_and_absorbed_power() {
    let gas = gas();
    let material = material();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for (rung, frequency) in [100.0, 550.0, 950.0].into_iter().enumerate() {
            let start = Instant::now();
            let omega = core::f64::consts::TAU * frequency;
            let load = reference_load(&material, &gas, omega);
            assert!(load.re > 0.0);
            if frequency == 100.0 {
                // A rigid-backed thin layer is compliance-like: positive
                // imaginary impedance under the NETWORK's exp(-i wt).
                assert!(load.im > 0.0);
            }
            let expected = cylinder_input(&gas, RADIUS, LENGTH, omega, loss, load);
            let net = two_ports();
            let response = network_port_impedance(&net, &gas, omega, loss).unwrap();
            let drive = response.drive("input", &[PortTermination {
                label: "outlet".into(), load: porous(&material),
            }], &gas).unwrap();
            let error = relative(drive.impedance, expected);
            assert!(error < 1e-11, "port load: {error}");
            let outlet = response.ports.iter().position(|p| p == "outlet").unwrap();
            let outgoing = -drive.port_inflows[outlet];
            assert!(relative(drive.port_pressures[outlet], load * outgoing) < 1e-11);
            let absorbed = 0.5 * load.re * outgoing.norm_sq();
            let input = 0.5 * drive.impedance.re;
            assert!(absorbed > 0.0);
            let balance = (input - absorbed) / input;
            if matches!(loss, LossModel::Lossless) {
                assert!(balance.abs() < 1e-11, "power defect: {balance}");
            } else {
                assert!(balance > 0.0, "duct wall loss must remain separate");
            }
            let mut terminated = net.clone();
            terminated.nodes[1].boundary = NodeBoundary::Termination(porous(&material));
            let assembled = network_port_impedance(&terminated, &gas, omega, loss).unwrap();
            assert!(relative(assembled.z[0], expected) < 1e-11);
            row("porous_gas_boundary_line_power", &net, &response, error, 1e-11, rung, start);
        }
    }
}

#[test]
fn changing_only_gas_or_pore_data_changes_solved_impedance_and_absorbed_power() {
    let base = gas();
    let composition = Composition::new(
        CompositionBasis::MoleFraction, vec![0.6, 0.2, 0.2],
    ).unwrap();
    let states = [
        base,
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 323.15, base.pressure).unwrap(),
        GasState::try_new(&GasSpec::dry_air_ussa1976(), base.temperature, 200_000.0).unwrap(),
        resolve_gas_mixture(
            &[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::CarbonDioxide],
            &composition, base.temperature, base.pressure,
        ).unwrap(),
    ];
    let original = material();
    let mut changed = original.reference_card().clone();
    changed.flow_resistivity *= 2.0;
    changed.provenance = "analytic sensitivity fixture with half the intrinsic permeability";
    let changed = PorousGasLaw::try_new(changed, 1.81e-5).unwrap();
    let solve = |gas: &GasState, material: &PorousGasLaw| {
        let omega = core::f64::consts::TAU * 550.0;
        let response = network_port_impedance(&two_ports(), gas, omega, LossModel::Lossless).unwrap();
        let drive = response.drive("input", &[PortTermination {
            label: "outlet".into(), load: porous(material),
        }], gas).unwrap();
        let outlet = response.ports.iter().position(|p| p == "outlet").unwrap();
        let absorbed = 0.5 * reference_load(material, gas, omega).re
            * drive.port_inflows[outlet].norm_sq();
        (drive.impedance, absorbed)
    };
    let baseline = solve(&base, &original);
    for (gas, material) in states[1..].iter().map(|gas| (gas, &original))
        .chain(core::iter::once((&base, &changed)))
    {
        let actual = solve(gas, material);
        assert!(relative(actual.0, baseline.0) > 1e-3);
        assert!((actual.1 / baseline.1 - 1.0).abs() > 1e-3);
        assert!((actual.1 / (0.5 * actual.0.re) - 1.0).abs() < 1e-11);
    }
}

#[test]
fn porous_source_geometry_and_stale_gas_refuse_in_both_solver_paths() {
    let gas = gas();
    let material = material();
    let omega = core::f64::consts::TAU * 550.0;
    let net = two_ports();
    let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    let bad = [
        NetworkLoad::PorousRigid { material: material.clone(), thickness_m: DEPTH, area_m2: 0.0 },
        NetworkLoad::PorousRigid { material: material.clone(), thickness_m: DEPTH, area_m2: f64::NAN },
        NetworkLoad::PorousRigid { material: material.clone(), thickness_m: -DEPTH, area_m2: area() },
    ];
    for load in bad {
        assert!(response.drive("input", &[PortTermination {
            label: "outlet".into(), load: load.clone(),
        }], &gas).is_err());
        let mut terminated = net.clone();
        terminated.nodes[1].boundary = NodeBoundary::Termination(load);
        assert!(network_port_impedance(&terminated, &gas, omega, LossModel::Lossless).is_err());
    }
    for frequency in [99.0, 4001.0] {
        let omega = core::f64::consts::TAU * frequency;
        let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
        assert!(response.drive("input", &[PortTermination {
            label: "outlet".into(), load: porous(&material),
        }], &gas).is_err());
        let mut terminated = net.clone();
        terminated.nodes[1].boundary = NodeBoundary::Termination(porous(&material));
        assert!(network_port_impedance(&terminated, &gas, omega, LossModel::Lossless).is_err());
    }
    let changed_gas = GasState::try_new(
        &GasSpec::dry_air_ussa1976(), 323.15, gas.pressure,
    ).unwrap();
    assert!(response.drive("input", &[PortTermination {
        label: "outlet".into(), load: porous(&material),
    }], &changed_gas).is_err());
}
