//! Independent finite-volume molecular/thermal port balances and decay laws.
use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use crate::multicomponent::tests::{gas as source_gas, source};
use crate::multicomponent::{Boundaries, Boundary, Cell, Conserved, MixtureEuler1d, StepControl};
use fs_qty::{DynViscosity, Length};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

pub(in crate::multicomponent::network) fn gas(
    count: usize,
    variable: bool,
    formation: f64,
) -> MixtureGas {
    let mut models = vec![source(
        "A",
        R,
        3.5,
        if variable { 0.1 } else { 0.0 },
        formation,
    )];
    if count == 2 {
        models.push(source(
            "B",
            if variable { 2.0 * R } else { R },
            if variable { 4.3 } else { 3.5 },
            if variable { 0.05 } else { 0.0 },
            -2.0 * formation,
        ));
    }
    source_gas(models)
}

pub(in crate::multicomponent::network) fn law(
    gas: &MixtureGas,
    mu: f64,
    bulk: f64,
    conductivity: f64,
    diffusivity: f64,
) -> ConstantTransport {
    let n = gas.species_count();
    let mut binary = vec![diffusivity; n * n];
    for k in 0..n {
        binary[k * n + k] = 0.0;
    }
    ConstantTransport::new(
        gas.source().clone(),
        "manufactured-network-molecular-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        binary,
    )
    .unwrap()
}

const AREA: f64 = 0.2;
const WIDTH: f64 = 0.4;
const PLENUM_DISTANCE: f64 = 0.15;

// Symmetry makes the two equal branch cells one homogeneous physical pool,
// coupled through two equal ports to the other pool. Interior molecular flux
// remains zero, so an independent two-volume exponential is available.
fn two_pools(
    gas: &MixtureGas,
    branch: Conserved,
    plenum: Conserved,
    volume: f64,
) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(WIDTH),
        AREA,
        vec![branch.clone(), branch],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(
        vec![branch],
        vec![super::super::Plenum::new(volume, plenum).unwrap()],
        [
            super::super::BranchEnd::Left,
            super::super::BranchEnd::Right,
        ]
        .map(|end| super::super::PlenumConnection {
            branch: 0,
            end,
            plenum: 0,
        })
        .to_vec(),
        || false,
    )
    .unwrap()
}

pub(in crate::multicomponent::network) fn one_port(
    gas: &MixtureGas,
    reversed: bool,
) -> MixtureNetwork1d {
    let densities = [0.8, 0.2];
    let plenum_densities = [0.2, 0.8];
    let pressure = |d: &[f64], t| {
        d.iter()
            .zip(gas.source().models())
            .map(|(rho, m)| rho * R / m.molar_mass().value())
            .sum::<f64>()
            * t
    };
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(WIDTH),
        AREA,
        vec![
            gas.conserved(
                &densities,
                if reversed { -0.3 } else { 0.3 },
                pressure(&densities, 2.3)
            )
            .unwrap();
            4
        ],
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(
        vec![branch],
        vec![
            super::super::Plenum::new(
                0.07,
                gas.conserved(&plenum_densities, 0.0, pressure(&plenum_densities, 1.1))
                    .unwrap(),
            )
            .unwrap(),
        ],
        vec![super::super::PlenumConnection {
            branch: 0,
            end: if reversed {
                super::super::BranchEnd::Left
            } else {
                super::super::BranchEnd::Right
            },
            plenum: 0,
        }],
        || false,
    )
    .unwrap()
}

pub(in crate::multicomponent::network) fn owned(q: Cell<'_>) -> Conserved {
    Conserved {
        partial_densities: q.partial_densities.to_vec(),
        momentum: q.momentum,
        energy: q.energy,
    }
}
pub(in crate::multicomponent::network) fn fields(
    network: &MixtureNetwork1d,
) -> (Vec<Vec<Conserved>>, Vec<Conserved>) {
    (
        network
            .branches()
            .iter()
            .map(|b| {
                (0..b.cell_count())
                    .map(|i| owned(b.cell(i).unwrap()))
                    .collect()
            })
            .collect(),
        network.plenums().iter().map(|p| owned(p.state())).collect(),
    )
}
pub(in crate::multicomponent::network) fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:.16e}, expected={b:.16e}, tolerance={tolerance:.3e}"
    );
}

fn binary_error(step: f64) -> f64 {
    let g = gas(2, false, 0.0);
    let volume = 0.07;
    let mut network = two_pools(
        &g,
        g.conserved(&[0.8, 0.2], 0.0, 2.0).unwrap(),
        g.conserved(&[0.2, 0.8], 0.0, 2.0).unwrap(),
        volume,
    );
    let mut op = NetworkDissipationOperator::new(
        &network,
        law(&g, 0.0, 0.0, 0.0, 0.03),
        vec![Length::new(PLENUM_DISTANCE); 2],
        || false,
    )
    .unwrap();
    let duration = 0.4;
    let before = network.totals();
    let report = op
        .advance_for(
            &mut network,
            Time::new(duration),
            DissipationControl {
                maximum_step_s: step,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap();
    let branch_volume = 2.0 * AREA * WIDTH;
    let conductance = 2.0 * AREA * 0.03 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let rate = conductance * (1.0 / branch_volume + 1.0 / volume);
    let mean = (branch_volume * 0.8 + volume * 0.2) / (branch_volume + volume);
    let decay = (-rate * duration).exp();
    let expected_branch = mean + (0.8 - mean) * decay;
    let expected_plenum = mean + (0.2 - mean) * decay;
    let mut error = 0.0_f64;
    for i in 0..2 {
        let q = network.branches()[0].cell(i).unwrap();
        error = error.max((q.partial_densities[0] - expected_branch).abs());
        close(q.partial_densities.iter().sum(), 1.0, 3e-13);
        close(g.primitive(q, || false).unwrap().temperature, 2.0, 3e-11);
    }
    let q = network.plenums()[0].state();
    error = error.max((q.partial_densities[0] - expected_plenum).abs());
    close(q.partial_densities.iter().sum(), 1.0, 3e-13);
    close(g.primitive(q, || false).unwrap().temperature, 2.0, 3e-11);
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(report.elapsed.value(), duration);
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
    for k in 0..2 {
        close(
            network.totals().species_masses[k],
            before.species_masses[k],
            2e-13,
        );
        close(report.conservation_residual.species_masses[k], 0.0, 2e-13);
    }
    close(report.absolute_energy_residual_j, 0.0, 2e-12);
    let ports = op.accepted_port_exchanges().unwrap();
    assert!(
        ports
            .iter()
            .all(|p| p.species_masses_to_plenum[0] > 0.0 && p.species_masses_to_plenum[1] < 0.0)
    );
    close(
        ports[0].species_masses_to_plenum[0],
        ports[1].species_masses_to_plenum[0],
        2e-14,
    );
    error
}

#[test]
fn binary_diffusion_crosses_both_oriented_ports_and_refines_to_two_volume_decay() {
    let coarse = binary_error(0.04);
    let fine = binary_error(0.02);
    eprintln!("network binary two-volume errors: {coarse:e} / {fine:e}");
    assert!(fine < 4e-6, "fine binary error={fine:e}");
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

fn heat_error(step: f64) -> f64 {
    let g = gas(1, false, 0.0);
    let volume = 0.07;
    let mut network = two_pools(
        &g,
        g.conserved(&[1.0], 0.0, 2.0).unwrap(),
        g.conserved(&[2.0], 0.0, 8.0).unwrap(),
        volume,
    );
    let initial = network.totals();
    let mut op = NetworkDissipationOperator::new(
        &network,
        law(&g, 0.0, 0.0, 0.05, 0.0),
        vec![Length::new(PLENUM_DISTANCE); 2],
        || false,
    )
    .unwrap();
    let duration = 0.6;
    let report = op
        .advance_for(
            &mut network,
            Time::new(duration),
            DissipationControl {
                maximum_step_s: step,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap();
    let capacity_branch = (2.0 * AREA * WIDTH) * 2.5;
    let capacity_plenum = volume * 2.0 * 2.5;
    let conductance = 2.0 * AREA * 0.05 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let rate = conductance * (1.0 / capacity_branch + 1.0 / capacity_plenum);
    let equilibrium =
        (2.0 * capacity_branch + 4.0 * capacity_plenum) / (capacity_branch + capacity_plenum);
    let decay = (-rate * duration).exp();
    let expected_branch = equilibrium + (2.0 - equilibrium) * decay;
    let expected_plenum = equilibrium + (4.0 - equilibrium) * decay;
    let mut error = 0.0_f64;
    for i in 0..2 {
        let q = network.branches()[0].cell(i).unwrap();
        assert_eq!(q.partial_densities, [1.0]);
        error = error.max((g.primitive(q, || false).unwrap().temperature - expected_branch).abs());
    }
    let q = network.plenums()[0].state();
    assert_eq!(q.partial_densities, [2.0]);
    error = error.max((g.primitive(q, || false).unwrap().temperature - expected_plenum).abs());
    assert_eq!(network.time().value(), 0.0);
    close(network.totals().energy, initial.energy, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    for port in op.accepted_port_exchanges().unwrap() {
        assert_eq!(port.mass_to_plenum(), 0.0);
        assert!(port.energy_to_plenum < 0.0);
    }
    error
}

#[test]
fn fourier_heat_crosses_ports_with_actual_lengths_volumes_and_heat_capacities() {
    let coarse = heat_error(0.06);
    let fine = heat_error(0.03);
    eprintln!("network thermal two-volume errors: {coarse:e} / {fine:e}");
    assert!(fine < 5e-6, "fine thermal error={fine:e}");
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn molecular_ports_close_species_enthalpy_and_support_impulse_with_reference_and_orientation_invariance()
 {
    let mut reference_fields = None;
    let mut reference_ports: Option<Vec<super::super::PortExchange>> = None;
    for (formation, reversed) in [(0.0, false), (1e12, false), (0.0, true)] {
        let g = gas(2, true, formation);
        let mut network = one_port(&g, reversed);
        let before = fields(&network);
        let mut op = NetworkDissipationOperator::new(
            &network,
            law(&g, 0.04, 0.02, 0.05, 0.03),
            vec![Length::new(0.1)],
            || false,
        )
        .unwrap();
        let report = op
            .advance_for(
                &mut network,
                Time::new(0.04),
                DissipationControl {
                    maximum_step_s: 0.005,
                    ..DissipationControl::default()
                },
                || false,
            )
            .unwrap();
        let port = op.accepted_port_exchanges().unwrap()[0];
        let branch_report = op.accepted_branch_reports().unwrap()[0];
        for k in 0..2 {
            assert_eq!(
                port.species_masses_to_plenum[k],
                -branch_report.boundary_change.species_masses[k]
            );
            close(
                0.07 * (network.plenums()[0].state().partial_densities[k]
                    - before.1[0].partial_densities[k]),
                port.species_masses_to_plenum[k],
                3e-14,
            );
            close(report.conservation_residual.species_masses[k], 0.0, 5e-14);
        }
        assert_eq!(port.energy_to_plenum, -branch_report.boundary_change.energy);
        assert_eq!(
            port.momentum_to_plenum,
            -branch_report.boundary_change.momentum
        );
        assert_eq!(network.plenums()[0].state().momentum, 0.0);
        close(port.mass_to_plenum(), 0.0, 5e-14);
        assert!(port.momentum_to_plenum.abs() > 1e-5);
        assert!(port.species_masses_to_plenum[0].abs() > 1e-5);
        assert!(port.energy_to_plenum.abs() > 1e-5);
        close(report.conservation_residual.energy, 0.0, 3e-12);
        if formation == 0.0 {
            // Reconstruct each owner's original physical energy, separately
            // from reported residuals or the reference-energy identity.
            let mut branch_energy_change = 0.0;
            let mut impulse = 0.0;
            for i in 0..4 {
                let q = network.branches()[0].cell(i).unwrap();
                let old = before.0[0][i].as_cell();
                branch_energy_change += AREA
                    * WIDTH
                    * (g.absolute_energy_density(q).unwrap()
                        - g.absolute_energy_density(old).unwrap());
                impulse += AREA * WIDTH * (q.momentum - old.momentum);
            }
            let plenum_energy_change = 0.07
                * (g.absolute_energy_density(network.plenums()[0].state())
                    .unwrap()
                    - g.absolute_energy_density(before.1[0].as_cell()).unwrap());
            let exchanged = port.absolute_energy_to_plenum(&g).unwrap();
            close(branch_energy_change, -exchanged, 3e-12);
            close(plenum_energy_change, exchanged, 3e-12);
            close(impulse, -port.momentum_to_plenum, 2e-13);
            close(report.absolute_energy_residual_j, 0.0, 3e-12);
        }
        let current = fields(&network);
        if !reversed {
            if let Some(reference) = &reference_fields {
                assert_eq!(
                    &current, reference,
                    "formation offsets changed molecular dynamics"
                );
                assert_eq!(
                    op.accepted_port_exchanges().unwrap(),
                    reference_ports.as_ref().unwrap().as_slice()
                );
            } else {
                reference_fields = Some(current);
                reference_ports = Some(op.accepted_port_exchanges().unwrap().to_vec());
            }
        } else {
            let reference = reference_fields.as_ref().unwrap();
            for i in 0..4 {
                let q = &current.0[0][i];
                let original = &reference.0[0][3 - i];
                for k in 0..2 {
                    close(q.partial_densities[k], original.partial_densities[k], 3e-12);
                }
                close(q.momentum, -original.momentum, 3e-12);
                close(q.energy, original.energy, 3e-12);
            }
            for k in 0..2 {
                close(
                    current.1[0].partial_densities[k],
                    reference.1[0].partial_densities[k],
                    3e-12,
                );
            }
            close(current.1[0].energy, reference.1[0].energy, 3e-12);
            let original = reference_ports.as_ref().unwrap()[0];
            close(port.momentum_to_plenum, -original.momentum_to_plenum, 3e-12);
        }
    }
}

#[test]
fn compact_plenum_limits_molecular_substeps_and_invalid_geometry_or_source_refuses() {
    let g = gas(2, false, 0.0);
    let make = |volume| {
        two_pools(
            &g,
            g.conserved(&[0.5, 0.5], 0.0, 2.0).unwrap(),
            g.conserved(&[0.5, 0.5], 0.0, 2.0).unwrap(),
            volume,
        )
    };
    let mut counts = Vec::new();
    for volume in [0.07, 0.0004] {
        let mut network = make(volume);
        let before = fields(&network);
        let mut op = NetworkDissipationOperator::new(
            &network,
            law(&g, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(PLENUM_DISTANCE); 2],
            || false,
        )
        .unwrap();
        let report = op
            .advance_for(
                &mut network,
                Time::new(0.04),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        counts.push(report.accepted_steps);
        assert_eq!(fields(&network), before);
        assert!(report.maximum_parabolic_number <= 0.40000000000001);
    }
    assert_eq!(counts[0], 1);
    assert!(counts[1] >= 10, "compact-volume substeps={counts:?}");
    let network = make(0.07);
    assert!(matches!(
        NetworkDissipationOperator::new(
            &network,
            law(&g, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(0.0); 2],
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert!(matches!(
        NetworkDissipationOperator::new(
            &network,
            law(&g, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(PLENUM_DISTANCE)],
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    let other = gas(2, true, 0.0);
    assert!(matches!(
        NetworkDissipationOperator::new(
            &network,
            law(&other, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(PLENUM_DISTANCE); 2],
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
}

#[test]
fn late_molecular_cancellation_preserves_every_owner_clock_and_accepted_report_then_reuses_work() {
    let g = gas(2, true, 5.0);
    let make = || {
        let mut n = one_port(&g, false);
        n.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let mut op = NetworkDissipationOperator::new(
            &n,
            law(&g, 0.04, 0.02, 0.05, 0.03),
            vec![Length::new(0.1)],
            || false,
        )
        .unwrap();
        op.advance_for(
            &mut n,
            Time::new(0.01),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
        // A later successful A step leaves real reports for the refused D
        // interval to preserve, rather than comparing only empty receipts.
        n.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        (n, op)
    };
    let control = DissipationControl {
        maximum_step_s: 0.005,
        ..DissipationControl::default()
    };
    let (mut probe, mut probe_op) = make();
    let mut polls = 0;
    probe_op
        .advance_for(&mut probe, Time::new(0.04), control, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 100);
    let (mut network, mut op) = make();
    let before = fields(&network);
    let clock = network.time();
    let advection_report = network.last_report();
    let advection_ports = network.accepted_port_exchanges().unwrap().to_vec();
    let advection_branches = network.accepted_branch_reports().unwrap().to_vec();
    let report = op.last_report();
    let ports = op.accepted_port_exchanges().unwrap().to_vec();
    let branches = op.accepted_branch_reports().unwrap().to_vec();
    let mut visits = 0;
    let result = op.advance_for(&mut network, Time::new(0.04), control, || {
        visits += 1;
        visits >= polls
    });
    assert_eq!(result, Err(DissipationError::Cancelled));
    assert_eq!(fields(&network), before);
    assert_eq!(network.time(), clock);
    assert_eq!(network.last_report(), advection_report);
    assert_eq!(network.accepted_port_exchanges().unwrap(), advection_ports);
    assert_eq!(
        network.accepted_branch_reports().unwrap(),
        advection_branches
    );
    assert_eq!(op.last_report(), report);
    assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(op.accepted_branch_reports().unwrap(), branches);
    op.advance_for(&mut network, Time::new(0.04), control, || false)
        .unwrap();
    assert_eq!(network.time(), probe.time());
    let a = fields(&network);
    let b = fields(&probe);
    for (x, y) in
        a.0.iter()
            .flatten()
            .chain(a.1.iter())
            .zip(b.0.iter().flatten().chain(b.1.iter()))
    {
        for k in 0..2 {
            close(x.partial_densities[k], y.partial_densities[k], 3e-12);
        }
        close(x.momentum, y.momentum, 3e-12);
        close(x.energy, y.energy, 3e-11);
    }
    assert!(network.last_report().is_none());
}
