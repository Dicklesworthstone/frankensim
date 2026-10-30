//! Combined accepted pressure and viscous traction, including tensile normals.
use super::*;
use fs_flux::compressible::spatial::{
    Boundaries3d, BoundarySide as Side3d, CartesianGrid as Grid3d, Euler3d, State as State3d,
    molecular::MolecularEuler3d,
};
use fs_flux::compressible::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide, CartesianGrid, Euler2d, State,
    molecular::{MolecularEuler2d, MolecularStepControl},
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

const AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const REFLECTING: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};

fn law(gas: &FrozenGas, viscosity: f64) -> ConstantTransport {
    let models = gas.source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0]).unwrap(),
        models,
    )
    .unwrap();
    ConstantTransport::new(
        source,
        "manufactured-normal-traction",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 100.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(viscosity),
            bulk_viscosity: DynViscosity::new(0.15 * viscosity),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.05 * viscosity),
        },
        vec![0.0],
    )
    .unwrap()
}

#[test]
fn complete_planar_molecular_traction_preserves_tension_and_the_original_p1_load() {
    let gas = gas();
    let transport = law(&gas, 0.2);
    let state = State::from_primitive(&gas, 1.0, [1.0, 0.0], 1.0).unwrap();
    let flow = Euler2d::new(
        gas,
        CartesianGrid::new(4, 1, 0.5, 1.0, 1.0).unwrap(),
        vec![state; 4],
        Boundaries2d {
            x: REFLECTING,
            y: REFLECTING,
        },
        || false,
    )
    .unwrap();
    let mut flow = MolecularEuler2d::new(flow, transport, || false).unwrap();
    let axes = [AXES[0], AXES[1]];
    assert!(
        WallTractionImpulse::from_molecular_euler2d(&flow, BoundarySide::XMin, 0, axes, 0.0)
            .is_err()
    );
    let report = flow
        .step(Time::new(1e-4), MolecularStepControl::default(), || false)
        .unwrap();
    let traction =
        WallTractionImpulse::from_molecular_euler2d(&flow, BoundarySide::XMin, 0, axes, 0.25)
            .unwrap();
    // The admitted synthetic viscosity gives tau_nn > p at this separating
    // rigid slip wall. The total traction is tensile despite positive gas p.
    assert!(traction.gas_traction_impulse_pa_s() < 0.0);
    assert!(traction.mean_gas_traction_pa() < 0.0);
    assert_eq!(traction.duration_s(), report.elapsed.value());
    assert_eq!(traction.end_time_s(), flow.flow().time().value());
    let complete = flow.boundary_impulses(BoundarySide::XMin)[0].momentum[0];
    let advective = flow.flow().boundary_impulses(BoundarySide::XMin)[0].momentum[0];
    assert!((complete - advective).abs() > 1e-5);
    close(traction.gas_traction_impulse_pa_s(), -complete);

    let mut left_mesh = mesh();
    for point in &mut left_mesh.nodes {
        point[0] = 0.0;
    }
    let left = ShockLoadSurface::new(
        &left_mesh,
        vec![selection(ShellFace::Negative)],
        budget(),
        || false,
    )
    .unwrap();
    let load = left
        .project_tractions(&[traction], [0.0, 0.5, 0.5], || false)
        .unwrap();
    let expected = complete + 0.25 * report.elapsed.value();
    assert!(expected > 0.0);
    close(load.resultant_impulse_n_s[0], expected);
    close(load.exterior_impulse_n_s[0], 0.25 * report.elapsed.value());
    let displacement: Vec<_> = left_mesh
        .nodes
        .iter()
        .map(|p| [1.0 + 2.0 * p[1] + 3.0 * p[2], 0.0, 0.0])
        .collect();
    close(
        load.virtual_work_j(&displacement).unwrap(),
        expected / report.elapsed.value() * 3.5,
    );
    let mut support = [0.0; 2];
    for side in [
        BoundarySide::XMin,
        BoundarySide::XMax,
        BoundarySide::YMin,
        BoundarySide::YMax,
    ] {
        for face in 0..flow.flow().grid().face_count(side) {
            let impulse =
                WallTractionImpulse::from_molecular_euler2d(&flow, side, face, axes, 0.0).unwrap();
            for (axis, sum) in support.iter_mut().enumerate() {
                *sum += impulse.outward_normal()[axis]
                    * impulse.gas_traction_impulse_pa_s()
                    * impulse.area_m2();
            }
        }
    }
    for axis in 0..2 {
        close(
            report.after.momentum[axis] - report.before.momentum[axis] + support[axis],
            0.0,
        );
    }
}

#[test]
fn volumetric_molecular_support_uses_the_complete_split_and_retains_it_on_cancellation() {
    let gas = gas();
    let transport = law(&gas, 0.2);
    let state = State3d::from_primitive(&gas, 1.0, [0.7, -0.2, 0.3], 1.0).unwrap();
    let flow = Euler3d::new(
        gas,
        Grid3d::new(4, 1, 1, 0.5, 1.0, 1.0).unwrap(),
        vec![state; 4],
        Boundaries3d {
            x: REFLECTING,
            y: REFLECTING,
            z: REFLECTING,
        },
        || false,
    )
    .unwrap();
    let mut flow = MolecularEuler3d::new(flow, transport, || false).unwrap();
    let mut calls = 0;
    let report = flow
        .step(Time::new(1e-4), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    let mut support = [0.0; 3];
    for side in Side3d::ALL {
        for face in 0..flow.flow().grid().face_count(side) {
            let impulse =
                WallTractionImpulse::from_molecular_euler3d(&flow, side, face, AXES, 0.0).unwrap();
            for (axis, sum) in support.iter_mut().enumerate() {
                *sum += impulse.outward_normal()[axis]
                    * impulse.gas_traction_impulse_pa_s()
                    * impulse.area_m2();
            }
        }
    }
    for axis in 0..3 {
        close(
            report.after.momentum[axis] - report.before.momentum[axis] + support[axis],
            0.0,
        );
    }
    let side = Side3d::XMax;
    let impulse = WallTractionImpulse::from_molecular_euler3d(&flow, side, 0, AXES, 0.25).unwrap();
    let raw = flow.boundary_impulses(side)[0].momentum[0];
    let middle = flow.flow().boundary_impulses(side)[0].momentum[0];
    assert!((raw - middle).abs() > 1e-5);
    let load = surface(ShellFace::Positive)
        .project_tractions(&[impulse], [2.0, 0.5, 0.5], || false)
        .unwrap();
    close(
        load.resultant_impulse_n_s[0],
        raw - 0.25 * report.elapsed.value(),
    );
    assert_eq!(
        surface(ShellFace::Positive)
            .project_tractions(&[impulse], [0.0; 3], || true)
            .unwrap_err(),
        Error::Cancelled
    );
    assert!(WallTractionImpulse::from_molecular_euler3d(&flow, side, 1, AXES, 0.0).is_err());
    assert!(
        WallTractionImpulse::from_molecular_euler3d(&flow, side, 0, [AXES[0]; 3], 0.0).is_err()
    );
    let states = flow.flow().states().to_vec();
    let accepted = *flow.last_step().unwrap();
    let mut polls = 0;
    assert!(
        flow.step(Time::new(1e-4), MolecularStepControl::default(), || {
            polls += 1;
            polls == calls / 2
        })
        .is_err()
    );
    assert!(polls > 1);
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.last_step(), Some(&accepted));
    assert_eq!(
        WallTractionImpulse::from_molecular_euler3d(&flow, side, 0, AXES, 0.25).unwrap(),
        impulse
    );
}
