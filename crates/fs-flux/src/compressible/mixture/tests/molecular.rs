use super::*;
use crate::compressible::mixture::molecular::{
    DissipationControl, DissipationError, MolecularReport, MolecularTransport2d, PlanarTransport,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, DissipationOperator, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

pub(super) fn law(
    gas: &MixtureGas,
    diffusion: f64,
    shear: f64,
    bulk: f64,
    conductivity: f64,
) -> ConstantTransport {
    let n = gas.species_count();
    let mut binary = vec![diffusion; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        gas.source().clone(),
        "manufactured-planar-molecular",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(shear),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        binary,
    )
    .unwrap()
}
pub(super) fn prepare(
    flow: &MixtureEuler2d,
    law: ConstantTransport,
) -> MolecularTransport2d<ConstantTransport> {
    let source = PlanarTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport2d::required_workspace_bytes(flow, &source).unwrap();
    let operator = MolecularTransport2d::new(flow, source, bytes, || false).unwrap();
    assert_eq!(operator.workspace_bytes() + flow.workspace_bytes(), bytes);
    operator
}
fn check(report: MolecularReport) {
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
    assert!(report.maximum_mass_frame_residual_kg_m2_s < 1e-11);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-11);
    close(report.absolute_energy_residual_j, 0.0, 3e-9);
}
fn equal_mass_gas() -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 4.5, 0.0, -50.0),
        source("absent-C", R, 3.5, 0.0, 2.0),
    ])
}

fn diffusion_mode(n: usize, diffusion: f64) -> (f64, f64) {
    let dx = 1.0 / n as f64;
    let z = std::f64::consts::PI * dx;
    let average = (z.sin() / z).powi(2);
    let wave = |x: [f64; 2]| (2.0 * std::f64::consts::PI * (x[0] + x[1])).sin();
    let grid = CartesianGrid::new(n, n, dx, dx, 0.6).unwrap();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.2 * average * wave(x);
        ([a, 1.0 - a, 0.0], [0.0; 2], 1.0)
    });
    let before = flow.states().to_vec();
    let mut operator = prepare(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let duration = 0.05;
    check(
        operator
            .advance(
                &mut flow,
                Time::new(duration),
                DissipationControl::default(),
                || false,
            )
            .unwrap(),
    );
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert_ne!(flow.states(), before);
    let decay = (-diffusion * 8.0 * std::f64::consts::PI.powi(2) * duration).exp();
    let mut error = 0.0;
    let mut projection = 0.0;
    let mut norm = 0.0;
    for index in 0..grid.len() {
        let q = flow.cell(index).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        let mode = wave(grid.cell_center(index).unwrap());
        error += (q.partial_densities[0] - (0.5 + 0.2 * average * decay * mode)).abs();
        projection += (q.partial_densities[0] - 0.5) * mode;
        norm += mode * mode;
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 2e-12);
        // Unequal original Cv and formation offsets: correct partial enthalpy
        // transport must change E_ref while retaining this isothermal solution.
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
    }
    (error / grid.len() as f64, projection / norm)
}

#[test]
fn g1_oblique_species_diffusion_carries_original_partial_enthalpy_and_refines() {
    let coarse = diffusion_mode(10, 0.02);
    let fine = diffusion_mode(20, 0.02);
    eprintln!(
        "2D Maxwell-Stefan oblique L1 10/20: {} / {}",
        coarse.0, fine.0
    );
    assert!(fine.0 < 0.3 * coarse.0);
    assert!(fine.0 < 1e-4);
    let changed_data = diffusion_mode(20, 0.04);
    assert!(
        changed_data.1 < fine.1,
        "only binary diffusion data changed"
    );
}

#[test]
fn g1_complete_planar_stress_matches_independent_oblique_vector_symbols() {
    let n = 12;
    let h = 1.0 / n as f64;
    let mu = 0.015;
    let bulk = 0.005;
    let amplitude = 0.02;
    let t: f64 = 0.02;
    let dt: f64 = 0.001;
    let wave = |x: [f64; 2]| (2.0 * std::f64::consts::PI * (x[0] + x[1])).sin();
    let theta = 2.0 * std::f64::consts::PI * h;
    let diagonal = (2.0 * mu + bulk + mu / 3.0) * 4.0 * (0.5 * theta).sin().powi(2) / (h * h);
    let cross = (bulk + mu / 3.0) * (theta.sin() / h).powi(2);
    for sign in [-1.0, 1.0] {
        let grid = CartesianGrid::new(n, n, h, h, 0.6).unwrap();
        let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
            let u = amplitude * wave(x);
            ([1.0, 0.0, 0.0], [u, sign * u], 1.0)
        });
        let kinetic = |flow: &MixtureEuler2d| {
            (0..flow.cell_count())
                .map(|i| {
                    let q = flow.cell(i).unwrap();
                    0.5 * (q.momentum[0].powi(2) + q.momentum[1].powi(2))
                })
                .sum::<f64>()
                * flow.grid.cell_volume()
        };
        let before_kinetic = kinetic(&flow);
        let before_energy = flow.totals().energy;
        let mut operator = prepare(&flow, law(flow.gas(), 0.001, mu, bulk, 0.0));
        let report = operator
            .advance(
                &mut flow,
                Time::new(t),
                DissipationControl {
                    maximum_step_s: dt,
                    ..DissipationControl::default()
                },
                || false,
            )
            .unwrap();
        check(report);
        let eigenvalue = diagonal + sign * cross;
        let exact =
            (1.0 - eigenvalue * dt + 0.5 * (eigenvalue * dt).powi(2)).powi((t / dt).round() as i32);
        for index in 0..grid.len() {
            let p = flow
                .cell(index)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            let u = amplitude * wave(grid.cell_center(index).unwrap()) * exact;
            close(p.velocity[0], u, 3e-12);
            close(p.velocity[1], sign * u, 3e-12);
        }
        close(flow.totals().energy, before_energy, 2e-12);
        let after_kinetic = kinetic(&flow);
        assert!(after_kinetic < before_kinetic);
        let thermal_gain = (0..grid.len())
            .map(|i| {
                let q = flow.cell(i).unwrap();
                let p = q.primitive(flow.gas(), || false).unwrap();
                2.5 * (p.temperature - 1.0) * grid.cell_volume()
            })
            .sum::<f64>();
        close(thermal_gain, before_kinetic - after_kinetic, 3e-12);
    }
}

#[test]
fn g3_evolving_planar_molecular_flow_matches_original_axial_owner_and_rotation() {
    let g = standard(0.0);
    let n = 12;
    let h = 1.0 / n as f64;
    let initial = |x: f64| {
        let angle = 2.0 * std::f64::consts::PI * x;
        let a = 0.4 + 0.1 * angle.sin();
        let b = 1.0 - a;
        let t = 1.0 + 0.1 * angle.cos();
        ([a, b, 0.0], 0.1 * angle.sin(), (a + 0.5 * b) * t)
    };
    let mut horizontal = field(
        g.clone(),
        CartesianGrid::new(n, 2, h, 0.2, 0.7).unwrap(),
        PERIODIC,
        |p| {
            let (rho, u, pressure) = initial(p[0]);
            (rho, [u, 0.3], pressure)
        },
    );
    let mut vertical = field(
        g.clone(),
        CartesianGrid::new(2, n, 0.2, h, 0.7).unwrap(),
        PERIODIC,
        |p| {
            let (rho, u, pressure) = initial(p[1]);
            (rho, [0.3, u], pressure)
        },
    );
    let states = (0..n)
        .map(|i| {
            let (rho, u, p) = initial((i as f64 + 0.5) * h);
            g.conserved(&rho, u, p).unwrap()
        })
        .collect();
    let mut axial = MixtureEuler1d::new(
        g.clone(),
        Length::new(h),
        0.14,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut ah =
        DissipationOperator::new(&axial, law(&g, 0.03, 0.02, 0.01, 0.03), || false).unwrap();
    let mut hh = prepare(&horizontal, law(&g, 0.03, 0.02, 0.01, 0.03));
    let mut vh = prepare(&vertical, law(&g, 0.03, 0.02, 0.01, 0.03));
    let controls = DissipationControl {
        maximum_step_s: 0.0005,
        ..DissipationControl::default()
    };
    ah.advance_for(&mut axial, Time::new(0.002), controls, || false)
        .unwrap();
    check(
        hh.advance(&mut horizontal, Time::new(0.002), controls, || false)
            .unwrap(),
    );
    check(
        vh.advance(&mut vertical, Time::new(0.002), controls, || false)
            .unwrap(),
    );
    for index in 0..n {
        let a = axial.cell(index).unwrap();
        for (planar, index, normal) in [(&horizontal, index, 0), (&vertical, 2 * index, 1)] {
            let q = planar.cell(index).unwrap();
            for k in 0..3 {
                close(q.partial_densities[k], a.partial_densities[k], 3e-10);
            }
            close(q.momentum[normal], a.momentum, 3e-10);
            close(q.momentum[1 - normal], 0.3, 3e-10);
            close(q.energy - 0.5 * 0.3_f64.powi(2), a.energy, 3e-10);
        }
    }
}

fn wall_field() -> MixtureEuler2d {
    field(
        standard(0.0),
        CartesianGrid::new(4, 3, 0.25, 1.0 / 3.0, 0.6).unwrap(),
        Boundaries2d { x: WALLS, y: WALLS },
        |p| {
            (
                [0.6, 0.4, 0.0],
                [0.15 + 0.05 * p[0], -0.2 + 0.1 * p[1]],
                0.8 + 0.1 * p[0],
            )
        },
    )
}
pub(super) fn impulses(
    operator: &MolecularTransport2d<ConstantTransport>,
) -> Vec<([f64; 3], [f64; 2], f64)> {
    SIDES
        .into_iter()
        .flat_map(|side| operator.boundary_impulses(side))
        .map(|q| {
            (
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                ],
                q.momentum,
                q.energy,
            )
        })
        .collect()
}

#[test]
fn g4_planar_molecular_walls_close_and_requested_duration_is_atomic() {
    let create = || {
        let flow = wall_field();
        let op = prepare(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.02));
        (flow, op)
    };
    let (mut reference, mut expected_op) = create();
    expected_op
        .advance(
            &mut reference,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = expected_op
        .advance(
            &mut reference,
            Time::new(0.002),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    check(expected);
    let mut outward = [0.0; 2];
    for (species, momentum, energy) in impulses(&expected_op) {
        assert_eq!(species, [0.0; 3]);
        assert_eq!(energy, 0.0);
        for axis in 0..2 {
            outward[axis] += momentum[axis];
        }
    }
    for axis in 0..2 {
        close(
            expected.after.momentum[axis] - expected.before.momentum[axis] + outward[axis],
            0.0,
            2e-12,
        );
    }
    let (mut flow, mut operator) = create();
    operator
        .advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let before = flow.states().to_vec();
    let before_impulses = impulses(&operator);
    let before_report = *operator.last_step().unwrap();
    let mut attempted = 0;
    assert_eq!(
        operator
            .advance(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || {
                    attempted += 1;
                    attempted >= polls
                }
            )
            .unwrap_err(),
        DissipationError::Cancelled
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert_eq!(operator.last_step(), Some(&before_report));
    assert_eq!(impulses(&operator), before_impulses);
    let failed = operator.advance(
        &mut flow,
        Time::new(0.002),
        DissipationControl {
            maximum_step_s: 0.0001,
            maximum_steps: 1,
            ..DissipationControl::default()
        },
        || false,
    );
    assert!(matches!(failed, Err(DissipationError::WorkLimit { .. })));
    assert_eq!(flow.states(), before);
    assert_eq!(impulses(&operator), before_impulses);
    assert_eq!(
        operator
            .advance(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || false
            )
            .unwrap(),
        expected
    );
    assert_eq!(flow.states(), reference.states());
}

#[test]
fn g0_planar_molecular_budget_capacity_and_boundary_ownership_are_explicit() {
    let mut flow = wall_field();
    // Carry actual caller excess capacity through successive O(1) publications.
    let mut states = Vec::with_capacity(flow.states.len() + 37);
    states.extend_from_slice(&flow.states);
    let extra = (states.capacity() - flow.states.capacity()) * size_of::<f64>();
    flow.states = states;
    flow.workspace_bytes += extra;
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.01, 0.02, 0.0, 0.03),
        || false,
    )
    .unwrap();
    let required = MolecularTransport2d::required_workspace_bytes(&flow, &source).unwrap();
    assert!(matches!(
        MolecularTransport2d::new(&flow, source, required - 1, || false),
        Err(DissipationError::Invalid(_))
    ));
    let mut op = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.0, 0.03));
    let combined = flow.workspace_bytes() + op.workspace_bytes();
    for _ in 0..3 {
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
        assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), combined);
    }
    flow.boundaries.x = OPEN;
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.01, 0.02, 0.0, 0.03),
        || false,
    )
    .unwrap();
    assert!(matches!(
        MolecularTransport2d::new(&flow, source, MAX_MIXTURE_WORKSPACE_BYTES, || false),
        Err(DissipationError::Invalid(_))
    ));
}
