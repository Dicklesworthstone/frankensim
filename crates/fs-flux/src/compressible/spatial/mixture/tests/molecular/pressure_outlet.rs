//! G0/G1/G4: declared diffusive outlet closure and original advective transfer.

use super::*;
use crate::compressible::spatial::mixture::molecular::{
    AdaptiveDissipationControl, MolecularEuler3d, MolecularStepControl, MolecularStepReport3d,
    PressureOutletPolicy, StabilizedDissipationControl,
};

fn transport(
    flow: &MixtureEuler3d,
    native: ConstantTransport,
    policy: PressureOutletPolicy,
) -> Result<MolecularTransport3d<ConstantTransport>, DissipationError> {
    let source = SpatialTransport::new(flow.gas().clone(), native, || false)?;
    let bytes = MolecularTransport3d::required_adaptive_workspace_bytes(flow, &source)?;
    MolecularTransport3d::new_with_pressure_outlet_policy(flow, source, policy, bytes, || false)
}

fn outlets(pressure_pa: f64) -> Boundaries3d {
    Boundaries3d {
        x: AxisBoundaries::Ends {
            lower: BoundaryCondition::PressureOutlet { pressure_pa },
            upper: BoundaryCondition::PressureOutlet { pressure_pa },
        },
        ..PERIODIC
    }
}

fn zero_outlet_transfer(operator: &MolecularTransport3d<ConstantTransport>) {
    for side in [BoundarySide::XMin, BoundarySide::XMax] {
        for load in operator.boundary_impulses(side) {
            assert!(load.species_masses.iter().all(|&value| value == 0.0));
            assert_eq!(load.momentum, [0.0; 3]);
            assert_eq!(load.energy, 0.0);
        }
    }
}

fn coupled_balance(report: &MolecularStepReport3d) {
    assert_eq!(report.elapsed, report.advection.elapsed);
    assert_eq!(report.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(report.second.elapsed, report.first.elapsed);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-11);
    close(report.absolute_energy_residual_j, 0.0, 3e-9);
}

#[test]
fn g0_pressure_only_molecular_exchange_requires_explicit_opt_in() {
    let grid = CartesianGrid::new(4, 1, 1, 0.25, 1.0, 1.0).unwrap();
    let flow = field(equal_mass_gas(), grid, outlets(1.0), |_| {
        ([0.4, 0.6, 0.0], [0.0; 3], 1.0)
    });
    let native = || law(flow.gas(), 0.03, 0.02, 0.01, 0.02);
    let source = SpatialTransport::new(flow.gas().clone(), native(), || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap();
    assert!(MolecularTransport3d::new(&flow, source, bytes, || false).is_err());
    assert_eq!(
        PressureOutletPolicy::default(),
        PressureOutletPolicy::Refuse
    );
    assert!(transport(&flow, native(), PressureOutletPolicy::Refuse).is_err());
    let admitted = transport(&flow, native(), PressureOutletPolicy::ZeroDiffusiveFlux).unwrap();
    assert_eq!(
        admitted.pressure_outlet_policy(),
        PressureOutletPolicy::ZeroDiffusiveFlux
    );

    // Opting into pressure outlets must not invent a bath on a plain open end.
    let open = field(
        equal_mass_gas(),
        grid,
        Boundaries3d {
            x: OPEN,
            ..PERIODIC
        },
        |_| ([0.4, 0.6, 0.0], [0.0; 3], 1.0),
    );
    assert!(
        transport(
            &open,
            law(open.gas(), 0.03, 0.02, 0.01, 0.02),
            PressureOutletPolicy::ZeroDiffusiveFlux
        )
        .is_err()
    );

    // The same selection changes no periodic or reflecting-wall arithmetic.
    for boundaries in [
        PERIODIC,
        Boundaries3d {
            x: WALLS,
            ..PERIODIC
        },
    ] {
        let make = || {
            field(equal_mass_gas(), grid, boundaries, |x| {
                let mode = (std::f64::consts::TAU * x[0]).cos();
                (
                    [0.5 + 0.1 * mode, 0.5 - 0.1 * mode, 0.0],
                    [0.0, 0.1 * mode, 0.0],
                    1.0,
                )
            })
        };
        let mut reference = make();
        let mut opted = make();
        let mut strict = transport(
            &reference,
            law(reference.gas(), 0.03, 0.02, 0.01, 0.02),
            PressureOutletPolicy::Refuse,
        )
        .unwrap();
        let mut enabled = transport(
            &opted,
            law(opted.gas(), 0.03, 0.02, 0.01, 0.02),
            PressureOutletPolicy::ZeroDiffusiveFlux,
        )
        .unwrap();
        let before = strict
            .advance(
                &mut reference,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let after = enabled
            .advance(
                &mut opted,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(after, before);
        assert_eq!(opted.states(), reference.states());
    }
}

#[test]
fn g1_outlet_species_and_vorticity_decay_match_independent_neumann_symbols() {
    let n = 8;
    let h = 1.0 / n as f64;
    let grid = CartesianGrid::new(n, 1, 1, h, 1.0, 1.0).unwrap();
    let mut flow = field(equal_mass_gas(), grid, outlets(1.0), |x| {
        let mode = (std::f64::consts::PI * x[0]).cos();
        (
            [0.5 + 0.15 * mode, 0.5 - 0.15 * mode, 0.0],
            [0.0, 0.2 * mode, 0.0],
            1.0,
        )
    });
    let mut curl_work = flow
        .vorticity_workspace(flow.required_vorticity_workspace_bytes().unwrap())
        .unwrap();
    let initial_curl = flow.vorticity(&mut curl_work, || false).unwrap().to_vec();
    let mut operator = transport(
        &flow,
        law(flow.gas(), 0.03, 0.02, 0.005, 0.0),
        PressureOutletPolicy::ZeroDiffusiveFlux,
    )
    .unwrap();
    let dt = 0.002;
    let report = operator
        .advance(
            &mut flow,
            Time::new(dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    check(report);
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(flow.time().value(), 0.0);
    zero_outlet_transfer(&operator);
    // Cell-centered cosine modes have zero physical face derivative at both
    // exits. The finite-volume eigenvalue is 4*kappa*sin(pi/2n)^2/h^2.
    let decay = |diffusivity: f64| {
        let z = 4.0 * diffusivity * (0.5 * std::f64::consts::PI / n as f64).sin().powi(2) * dt
            / (h * h);
        1.0 - z + 0.5 * z * z
    };
    for i in 0..n {
        let mode = (std::f64::consts::PI * grid.cell_center(i).unwrap()[0]).cos();
        let q = flow.cell(i).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(
            q.partial_densities[0],
            0.5 + 0.15 * decay(0.03) * mode,
            3e-12,
        );
        close(
            q.partial_densities[1],
            0.5 - 0.15 * decay(0.03) * mode,
            3e-12,
        );
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 3e-12);
        close(p.velocity[1], 0.2 * decay(0.02) * mode, 3e-12);
    }
    for (initial, final_state) in initial_curl
        .iter()
        .zip(flow.vorticity(&mut curl_work, || false).unwrap())
    {
        close(
            final_state.vorticity[2],
            initial.vorticity[2] * decay(0.02),
            3e-12,
        );
    }
}

#[test]
fn g1_stabilized_and_adaptive_outlet_stages_retain_uniform_mixture() {
    for adaptive in [false, true] {
        let grid = CartesianGrid::new(4, 1, 1, 0.25, 1.0, 1.0).unwrap();
        let mut flow = field(equal_mass_gas(), grid, outlets(1.0), |_| {
            ([0.4, 0.6, 0.0], [0.0, 0.2, -0.1], 1.0)
        });
        let original = flow.states().to_vec();
        let mut operator = transport(
            &flow,
            law(flow.gas(), 0.03, 0.02, 0.01, 0.02),
            PressureOutletPolicy::ZeroDiffusiveFlux,
        )
        .unwrap();
        let report = if adaptive {
            operator
                .advance_adaptive(
                    &mut flow,
                    Time::new(0.3),
                    AdaptiveDissipationControl::default(),
                    || false,
                )
                .unwrap()
                .stabilized
        } else {
            operator
                .advance_stabilized(
                    &mut flow,
                    Time::new(0.3),
                    StabilizedDissipationControl::default(),
                    || false,
                )
                .unwrap()
        };
        assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
        for value in report.transport.conservation_residual.species_masses {
            close(value, 0.0, 3e-11);
        }
        for value in report.transport.conservation_residual.momentum {
            close(value, 0.0, 3e-11);
        }
        close(report.transport.conservation_residual.energy, 0.0, 3e-11);
        close(report.transport.absolute_energy_residual_j, 0.0, 3e-9);
        zero_outlet_transfer(&operator);
        for (&actual, expected) in flow.states().iter().zip(original) {
            close(actual, expected, 3e-12);
        }
    }
}

#[test]
fn g1_coupled_outlet_retains_all_original_advective_components_and_bath() {
    let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap();
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let opposite = BoundarySide::ALL[2 * axis + usize::from(side.sign() < 0.0)];
        let outlet = BoundaryCondition::PressureOutlet { pressure_pa: 1.0 };
        let mut axes = [AxisBoundaries::Periodic; 3];
        axes[axis] = AxisBoundaries::Ends {
            lower: if side.sign() < 0.0 {
                outlet
            } else {
                BoundaryCondition::Transmissive
            },
            upper: if side.sign() > 0.0 {
                outlet
            } else {
                BoundaryCondition::Transmissive
            },
        };
        let boundaries = Boundaries3d {
            x: axes[0],
            y: axes[1],
            z: axes[2],
        };
        let densities = [0.35, 0.65, 0.0];
        let mut velocity = [0.1, -0.2, 0.3];
        velocity[axis] = 0.4 * side.sign();
        let plain = field(standard(0.4), grid, boundaries, |_| {
            (densities, velocity, 1.0)
        });
        let initial = plain.states().to_vec();
        let energy = plain.cell(0).unwrap().energy;
        let bath = ReservoirBoundary {
            side: opposite,
            partial_densities: &densities,
            velocity,
            pressure_pa: 1.0,
        };
        let bytes = MixtureEuler3d::required_workspace_bytes_with_reservoirs(plain.gas(), &grid, 1)
            .unwrap();
        let flow = MixtureEuler3d::new_with_reservoirs(
            plain.gas().clone(),
            grid,
            initial.clone(),
            boundaries,
            &[bath],
            bytes,
            || false,
        )
        .unwrap();
        let operator = transport(
            &flow,
            law(flow.gas(), 0.03, 0.02, 0.01, 0.02),
            PressureOutletPolicy::ZeroDiffusiveFlux,
        )
        .unwrap();
        let bytes = MolecularEuler3d::required_workspace_bytes(&flow, &operator).unwrap();
        let mut owner = MolecularEuler3d::new(flow, operator, bytes, || false).unwrap();
        let report = owner
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        coupled_balance(&report);
        assert_eq!(
            owner.molecular().pressure_outlet_policy(),
            PressureOutletPolicy::ZeroDiffusiveFlux
        );
        for (&actual, expected) in owner.flow().states().iter().zip(initial) {
            close(actual, expected, 3e-11);
        }
        let transfer = report.elapsed.value() * grid.cell_volume() / grid.spacing()[axis];
        for load in owner.boundary_impulses(side) {
            for (species, density) in densities.into_iter().enumerate() {
                close(
                    load.species_masses[species],
                    0.4 * density * transfer,
                    3e-12,
                );
            }
            for component in 0..3 {
                let pressure = if component == axis { side.sign() } else { 0.0 };
                close(
                    load.momentum[component],
                    (0.4 * velocity[component] + pressure) * transfer,
                    3e-12,
                );
            }
            close(load.energy, 0.4 * (energy + 1.0) * transfer, 3e-11);
        }
        for load in owner.molecular().boundary_impulses(side) {
            assert!(load.species_masses.iter().all(|&value| value == 0.0));
            assert_eq!(load.momentum, [0.0; 3]);
            assert_eq!(load.energy, 0.0);
        }
    }
}

fn draining_shear() -> MolecularEuler3d<ConstantTransport> {
    let grid = CartesianGrid::new(8, 1, 1, 0.125, 1.0, 1.0).unwrap();
    let flow = field(equal_mass_gas(), grid, outlets(0.8), |x| {
        let mode = (std::f64::consts::PI * x[0]).cos();
        (
            [0.5 + 0.15 * mode, 0.5 - 0.15 * mode, 0.0],
            [0.0, 0.2 * mode, 0.1 * mode],
            1.0,
        )
    });
    let operator = transport(
        &flow,
        law(flow.gas(), 0.03, 0.02, 0.005, 0.02),
        PressureOutletPolicy::ZeroDiffusiveFlux,
    )
    .unwrap();
    let bytes = MolecularEuler3d::required_workspace_bytes(&flow, &operator).unwrap();
    MolecularEuler3d::new(flow, operator, bytes, || false).unwrap()
}

#[test]
fn g4_sheared_species_drain_and_late_cancellation_preserve_coupled_accounting() {
    let mut reference = draining_shear();
    let mut polls = 0;
    let expected = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    coupled_balance(&expected);
    assert!(expected.boundary_change.mass() < 0.0);
    assert!(expected.boundary_change.species_masses[0] < 0.0);
    assert!(expected.boundary_change.species_masses[1] < 0.0);
    zero_outlet_transfer(reference.molecular());

    let mut owner = draining_shear();
    let original = owner.flow().states().to_vec();
    let time = owner.flow().time();
    let mut seen = 0;
    assert_eq!(
        owner
            .step(Time::new(0.001), MolecularStepControl::default(), || {
                seen += 1;
                seen == polls
            })
            .unwrap_err(),
        DissipationError::Cancelled
    );
    assert_eq!(owner.flow().states(), original);
    assert_eq!(owner.flow().time(), time);
    assert!(owner.last_step().is_none());
    assert!(owner.molecular().last_step().is_none());
    zero_outlet_transfer(owner.molecular());
    let recovered = owner
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(owner.flow().states(), reference.flow().states());
}
