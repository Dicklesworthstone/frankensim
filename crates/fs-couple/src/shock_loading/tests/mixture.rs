//! G0/G3/G4: actual composition/reaction shock transfers into the original shell port.
use super::*;
use fs_flux::compressible::mixture::{
    MixtureEuler2d,
    molecular::{DissipationControl, MolecularTransport2d, PlanarTransport},
    reactive::{ReactingEuler2d, ReactiveStepControl as PlanarReactiveControl},
    state_into as planar_state,
};
use fs_flux::compressible::spatial::mixture::{
    MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler3d,
    molecular::{MolecularTransport3d, SpatialTransport},
    reactive::{ReactingEuler3d, ReactiveStepControl},
    state_into,
};
use fs_flux::compressible::spatial::{Boundaries3d, BoundarySide, CartesianGrid};
use fs_flux::compressible::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide as PlanarSide,
    CartesianGrid as PlanarGrid,
};
use fs_gas::chemistry::IsochoricKineticsSolve;
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_gas::multicomponent::reactive::CellChemistry;
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId};

const AXES: [[f64; 3]; 3] = [
    [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    [-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
];
const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const CLOSED: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};

fn source(cp_b: f64, heat: f64) -> IdealGasEquilibriumV1 {
    let frozen = gas();
    let original = &frozen.source().components()[0];
    let models = [("A", 3.5, 0.0), ("B", cp_b, -heat)]
        .map(|(name, cp, formation)| {
            let mut card = original.regions()[0].card().clone();
            card.parameters.get_mut("a2").unwrap().value = cp;
            card.parameters.get_mut("a7").unwrap().value = formation;
            Nasa9StandardStateModelV1::new(
                SpeciesId::new(name).unwrap(),
                MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
                original.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .to_vec();
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(),
        models,
    )
    .unwrap()
}

fn mixture(source: IdealGasEquilibriumV1) -> MixtureGas {
    MixtureGas::new(
        source,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(100.0)),
        || false,
    )
    .unwrap()
}

fn field(
    source: IdealGasEquilibriumV1,
    boundaries: Boundaries3d,
    nonuniform: bool,
) -> MixtureEuler3d {
    let gas = mixture(source);
    let grid = CartesianGrid::new(4, 2, 2, 0.5, 0.5, 0.5).unwrap();
    let mut states = vec![0.0; 16 * 6];
    for (i, q) in states.chunks_exact_mut(6).enumerate() {
        let y = (i / 4) % 2;
        let z = i / 8;
        let (densities, velocity, pressure) = if nonuniform {
            let a = 0.15 + 0.1 * (y + z) as f64;
            (
                [a, 1.0 - a],
                [1.0, 0.0, 0.0],
                1.0 + 2.0 * y as f64 + 3.0 * z as f64,
            )
        } else {
            ([1.0, 0.0], [0.0; 3], 1.0)
        };
        state_into(&gas, &densities, velocity, pressure, q).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| AXES[j][i] * v[j]).sum())
}
fn translated(v: [f64; 3]) -> [f64; 3] {
    let v = rotate(v);
    core::array::from_fn(|i| v[i] + [0.25, -0.5, 1.0][i])
}
fn footprint() -> ShockLoadSurface {
    let nodes = (0..9)
        .map(|i| translated([2.0, 0.5 * (i % 3) as f64, 0.5 * (i / 3) as f64]))
        .collect();
    let mut triangles = Vec::new();
    let patches = (0..4)
        .map(|face| {
            let node = face % 2 + 3 * (face / 2);
            triangles.push([node, node + 1, node + 4]);
            triangles.push([node, node + 4, node + 3]);
            (2 * face..2 * face + 2)
                .map(|triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Positive,
                })
                .collect()
        })
        .collect();
    let mesh = ShellMesh::new(nodes, triangles).unwrap();
    ShockLoadSurface::new(
        &mesh,
        patches,
        SurfaceBudget {
            geometry: ProfileBudget {
                max_nodes: 9,
                max_triangles: 8,
                max_feature_evaluations: 0,
            },
            maximum_patches: 4,
            maximum_preparation_terms: 400,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn evolving_species_footprint_retains_resultant_both_moments_and_affine_work() {
    let mut flow = field(source(4.5, 0.0), CLOSED, true);
    assert!(
        WallPressureImpulse::from_mixture_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0).is_err()
    );
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let exterior = 0.2;
    let impulses = (0..4)
        .map(|face| {
            WallPressureImpulse::from_mixture_euler3d(
                &flow,
                BoundarySide::XMax,
                face,
                AXES,
                exterior,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let load = footprint()
        .project(&impulses, translated([2.0, 0.5, 0.5]), || false)
        .unwrap();
    let dt = report.elapsed.value();
    let (mut total, mut moment, mut work) = (0.0, [0.0; 3], 0.0);
    let mut magnitudes = Vec::new();
    for (face, raw) in flow.boundary_impulses(BoundarySide::XMax).enumerate() {
        assert!(raw.species_masses.iter().all(|mass| *mass == 0.0));
        assert_eq!(raw.energy, 0.0);
        let y = 0.25 + 0.5 * (face % 2) as f64;
        let z = 0.25 + 0.5 * (face / 2) as f64;
        let impulse = raw.momentum[0] - exterior * 0.25 * dt;
        magnitudes.push(impulse);
        total += impulse;
        moment[1] += (z - 0.5) * impulse;
        moment[2] -= (y - 0.5) * impulse;
        work += impulse / dt * (1.0 + 2.0 * y + 3.0 * z);
    }
    assert!(magnitudes[1] > magnitudes[0] && magnitudes[2] > magnitudes[0]);
    assert!(moment[1] > 0.0 && moment[2] < 0.0);
    for (a, b) in load
        .resultant_impulse_n_s
        .into_iter()
        .zip(rotate([total, 0.0, 0.0]))
    {
        close(a, b);
    }
    for (a, b) in load.angular_impulse_n_m_s.into_iter().zip(rotate(moment)) {
        close(a, b);
    }
    let displacement = (0..9)
        .map(|i| rotate([1.0 + (i % 3) as f64 + 1.5 * (i / 3) as f64, 0.0, 0.0]))
        .collect::<Vec<_>>();
    close(load.virtual_work_j(&displacement).unwrap(), work);
    let mut support = [0.0; 3];
    for side in BoundarySide::ALL {
        for face in 0..flow.grid().face_count(side) {
            let p =
                WallPressureImpulse::from_mixture_euler3d(&flow, side, face, AXES, 0.0).unwrap();
            for (sum, n) in support.iter_mut().zip(p.outward_normal()) {
                *sum += n * p.area_m2() * p.gas_pressure_impulse_pa_s();
            }
        }
    }
    let change = rotate(core::array::from_fn(|i| {
        report.after.momentum[i] - report.before.momentum[i]
    }));
    for i in 0..3 {
        close(change[i] + support[i], 0.0);
    }

    let before = flow.states().to_vec();
    let accepted = impulses[0];
    let mut polls = 0;
    assert!(
        flow.step(Time::new(0.001), controls(), || {
            polls += 1;
            polls == 80
        })
        .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(
        WallPressureImpulse::from_mixture_euler3d(&flow, BoundarySide::XMax, 0, AXES, exterior)
            .unwrap(),
        accepted
    );
}

#[test]
fn planar_species_wall_loads_close_momentum_and_reject_transport_and_bad_geometry() {
    let gas = mixture(source(4.5, 0.0));
    let grid = PlanarGrid::new(2, 3, 0.2, 0.3, 0.4).unwrap();
    let mut states = vec![0.0; 6 * 5];
    for q in states.chunks_exact_mut(5) {
        planar_state(&gas, &[0.4, 0.6], [0.5, -0.3], 1.0, q).unwrap();
    }
    let bytes = MixtureEuler2d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler2d::new(
        gas,
        grid,
        states,
        Boundaries2d { x: WALL, y: WALL },
        bytes,
        || false,
    )
    .unwrap();
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let axes = [AXES[0], AXES[1]];
    let mut support = [0.0; 3];
    for side in [
        PlanarSide::XMin,
        PlanarSide::XMax,
        PlanarSide::YMin,
        PlanarSide::YMax,
    ] {
        for face in 0..flow.grid().face_count(side) {
            let p =
                WallPressureImpulse::from_mixture_euler2d(&flow, side, face, axes, 0.0).unwrap();
            for (sum, n) in support.iter_mut().zip(p.outward_normal()) {
                *sum += n * p.area_m2() * p.gas_pressure_impulse_pa_s();
            }
        }
    }
    let change = rotate([
        report.after.momentum[0] - report.before.momentum[0],
        report.after.momentum[1] - report.before.momentum[1],
        0.0,
    ]);
    for i in 0..3 {
        close(change[i] + support[i], 0.0);
    }
    assert!(
        WallPressureImpulse::from_mixture_euler2d(&flow, PlanarSide::XMax, 3, axes, 0.0).is_err()
    );
    assert!(
        WallPressureImpulse::from_mixture_euler2d(&flow, PlanarSide::XMax, 0, [AXES[0]; 2], 0.0)
            .is_err()
    );
    assert!(
        WallPressureImpulse::from_mixture_euler2d(&flow, PlanarSide::XMax, 0, axes, -1.0).is_err()
    );

    let chemistry = CellChemistry::new(
        flow.gas().clone(),
        kinetics(flow.gas().source().clone()),
        chemistry(),
        || false,
    )
    .unwrap();
    let mut reacting =
        ReactingEuler2d::new(flow, chemistry, MAX_MIXTURE_WORKSPACE_BYTES, || false).unwrap();
    assert!(
        WallPressureImpulse::from_reacting_euler2d(&reacting, PlanarSide::XMax, 0, axes, 0.0)
            .is_err()
    );
    let report = reacting
        .step(
            Time::new(0.001),
            PlanarReactiveControl {
                transport: controls(),
                ..PlanarReactiveControl::default()
            },
            || false,
        )
        .unwrap();
    assert!(
        report
            .chemistry
            .species_change_kg
            .iter()
            .any(|mass| mass.abs() > 1e-10)
    );
    let mut support = [0.0; 3];
    for side in [
        PlanarSide::XMin,
        PlanarSide::XMax,
        PlanarSide::YMin,
        PlanarSide::YMax,
    ] {
        for face in 0..reacting.flow().grid().face_count(side) {
            let p = WallPressureImpulse::from_reacting_euler2d(&reacting, side, face, axes, 0.0)
                .unwrap();
            let traction =
                WallTractionImpulse::from_reacting_euler2d(&reacting, side, face, axes, 0.0)
                    .unwrap();
            assert_eq!(
                traction.gas_traction_impulse_pa_s(),
                p.gas_pressure_impulse_pa_s(),
            );
            assert_eq!(p.duration_s(), report.elapsed.value());
            assert_eq!(p.end_time_s(), reacting.flow().time().value());
            for (sum, n) in support.iter_mut().zip(p.outward_normal()) {
                *sum += n * p.area_m2() * p.gas_pressure_impulse_pa_s();
            }
        }
    }
    let change = rotate([
        report.after.momentum[0] - report.before.momentum[0],
        report.after.momentum[1] - report.before.momentum[1],
        0.0,
    ]);
    for i in 0..3 {
        close(change[i] + support[i], 0.0);
    }

    for pair in [
        AxisBoundaries::Periodic,
        AxisBoundaries::Ends {
            lower: BoundaryCondition::Transmissive,
            upper: BoundaryCondition::Transmissive,
        },
    ] {
        let mut flow = field(source(4.5, 0.0), Boundaries3d { x: pair, ..CLOSED }, false);
        flow.step(Time::new(0.001), controls(), || false).unwrap();
        for side in [BoundarySide::XMin, BoundarySide::XMax] {
            assert!(WallPressureImpulse::from_mixture_euler3d(&flow, side, 0, AXES, 0.0).is_err());
        }
        assert!(
            WallPressureImpulse::from_mixture_euler3d(&flow, BoundarySide::YMax, 0, AXES, 0.0)
                .is_ok()
        );
    }
}

fn kinetics(source: IdealGasEquilibriumV1) -> IdealGasKineticsV1 {
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("shock-support-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(50.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    IdealGasKineticsV1::new(
        source,
        "manufactured-reacting-shock-support",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap()
}
fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 100.0,
        energy_tolerance_j: 1e-10,
        relative_tolerance: 2e-6,
        absolute_amount_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 0.01,
        maximum_steps: 30000,
        ..IsochoricKineticsSolve::default()
    }
}
fn reacting() -> ReactingEuler3d {
    let source = source(3.5, 15.0);
    ReactingEuler3d::new(
        field(source.clone(), CLOSED, false),
        kinetics(source),
        chemistry(),
        MAX_MIXTURE_WORKSPACE_BYTES,
        || false,
    )
    .unwrap()
}

#[test]
fn reacting_wall_loading_uses_the_complete_accepted_interval_and_not_endpoint_pressure() {
    let mut flow = reacting();
    assert!(
        WallPressureImpulse::from_reacting_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0)
            .is_err()
    );
    let control = ReactiveStepControl {
        transport: controls(),
        ..ReactiveStepControl::default()
    };
    let report = flow.step(Time::new(0.005), control, || false).unwrap();
    assert!(report.chemistry.species_change_kg[1] > 0.0);
    let mut support = [0.0; 3];
    for side in BoundarySide::ALL {
        for face in 0..flow.flow().grid().face_count(side) {
            let p =
                WallPressureImpulse::from_reacting_euler3d(&flow, side, face, AXES, 0.0).unwrap();
            assert_eq!(p.duration_s(), report.elapsed.value());
            assert_eq!(p.end_time_s(), flow.flow().time().value());
            let raw = flow.boundary_impulses(side).nth(face).unwrap();
            close(
                p.gas_pressure_impulse_pa_s() * p.area_m2(),
                side.sign() * raw.momentum[side.axis()],
            );
            for (sum, n) in support.iter_mut().zip(p.outward_normal()) {
                *sum += n * p.gas_pressure_impulse_pa_s() * p.area_m2();
            }
        }
    }
    let change = rotate(core::array::from_fn(|i| {
        report.after.momentum[i] - report.before.momentum[i]
    }));
    for i in 0..3 {
        close(change[i] + support[i], 0.0);
    }
    let impulses = (0..4)
        .map(|face| {
            WallPressureImpulse::from_reacting_euler3d(&flow, BoundarySide::XMax, face, AXES, 0.2)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let final_pressure = flow
        .flow()
        .cell(3)
        .unwrap()
        .primitive(flow.flow().gas(), || false)
        .unwrap()
        .pressure;
    assert!(
        final_pressure > 1.1 * impulses[0].mean_gas_pressure_pa(),
        "final p={final_pressure}, interval p={}",
        impulses[0].mean_gas_pressure_pa()
    );
    let load = footprint()
        .project(&impulses, translated([2.0, 0.5, 0.5]), || false)
        .unwrap();
    let direct = flow
        .boundary_impulses(BoundarySide::XMax)
        .map(|face| face.momentum[0])
        .sum::<f64>()
        - 0.2 * report.elapsed.value();
    for (a, b) in load
        .resultant_impulse_n_s
        .into_iter()
        .zip(rotate([direct, 0.0, 0.0]))
    {
        close(a, b);
    }
    let displacement = vec![AXES[0]; 9];
    close(
        load.virtual_work_j(&displacement).unwrap(),
        direct / report.elapsed.value(),
    );

    let mut twin = reacting();
    twin.step(Time::new(0.005), control, || false).unwrap();
    let mut completed_polls = 0;
    twin.step(Time::new(0.005), control, || {
        completed_polls += 1;
        false
    })
    .unwrap();
    assert!(completed_polls > 10);
    let states = flow.flow().states().to_vec();
    let accepted = impulses[0];
    let mut polls = 0;
    assert!(
        flow.step(Time::new(0.005), control, || {
            polls += 1;
            polls == completed_polls
        })
        .is_err()
    );
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(
        WallPressureImpulse::from_reacting_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.2)
            .unwrap(),
        accepted
    );
}

#[test]
fn spatial_reacting_molecular_traction_retains_signed_support_and_three_dimensional_footprint() {
    let input = field(source(3.5, 15.0), CLOSED, true);
    let law = ConstantTransport::new(
        input.gas().source().clone(),
        "manufactured-spatial-reacting-wall-traction",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 100.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(4.0),
            bulk_viscosity: DynViscosity::new(0.03),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.01),
        },
        vec![0.0, 0.02, 0.02, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(input.gas().clone(), law, || false).unwrap();
    let mut molecular =
        MolecularTransport3d::new(&input, source, MAX_MIXTURE_WORKSPACE_BYTES, || false).unwrap();
    let mechanism = kinetics(input.gas().source().clone());
    let mut flow = ReactingEuler3d::new(
        input,
        mechanism,
        chemistry(),
        MAX_MIXTURE_WORKSPACE_BYTES,
        || false,
    )
    .unwrap();
    let bytes = flow.workspace_bytes() + molecular.workspace_bytes();
    let report = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(1e-4),
            ReactiveStepControl {
                transport: controls(),
                ..ReactiveStepControl::default()
            },
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    assert!(report.chemistry.species_change_kg[1] > 0.0);
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    let left = WallTractionImpulse::from_reacting_euler3d(&flow, BoundarySide::XMin, 0, AXES, 0.0)
        .unwrap();
    let right = WallTractionImpulse::from_reacting_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0)
        .unwrap();
    assert!(left.gas_traction_impulse_pa_s() < 0.0);
    assert!(right.gas_traction_impulse_pa_s() > 0.0);
    // Even the compressive side is total stress; a positive value does not
    // make a molecular traction interchangeable with thermodynamic pressure.
    for side in [BoundarySide::XMin, BoundarySide::XMax] {
        assert!(WallPressureImpulse::from_reacting_euler3d(&flow, side, 0, AXES, 0.0).is_err());
    }
    let mut support = [0.0; 3];
    for side in BoundarySide::ALL {
        for face in 0..flow.flow().grid().face_count(side) {
            let traction =
                WallTractionImpulse::from_reacting_euler3d(&flow, side, face, AXES, 0.0).unwrap();
            let raw = flow.boundary_impulses(side).nth(face).unwrap();
            assert_eq!(traction.duration_s(), report.elapsed.value());
            assert_eq!(traction.end_time_s(), flow.flow().time().value());
            close(
                traction.gas_traction_impulse_pa_s() * traction.area_m2(),
                side.sign() * raw.momentum[side.axis()],
            );
            for axis in 0..3 {
                support[axis] += traction.outward_normal()[axis]
                    * traction.area_m2()
                    * traction.gas_traction_impulse_pa_s();
            }
        }
    }
    let change = rotate(core::array::from_fn(|axis| {
        report.after.momentum[axis] - report.before.momentum[axis]
    }));
    for axis in 0..3 {
        close(change[axis] + support[axis], 0.0);
        close(
            report.boundary_change.momentum[axis],
            report.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
        );
    }
    let exterior = 0.2;
    let impulses = (0..4)
        .map(|face| {
            WallTractionImpulse::from_reacting_euler3d(
                &flow,
                BoundarySide::XMax,
                face,
                AXES,
                exterior,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let load = footprint()
        .project_tractions(&impulses, translated([2.0, 0.5, 0.5]), || false)
        .unwrap();
    let dt = report.elapsed.value();
    let (mut total, mut moment, mut work) = (0.0, [0.0; 3], 0.0);
    for (face, raw) in flow.boundary_impulses(BoundarySide::XMax).enumerate() {
        let middle = flow
            .flow()
            .boundary_impulses(BoundarySide::XMax)
            .nth(face)
            .unwrap();
        assert!((raw.momentum[0] - middle.momentum[0]).abs() > 1e-5);
        let y = 0.25 + 0.5 * (face % 2) as f64;
        let z = 0.25 + 0.5 * (face / 2) as f64;
        let net = raw.momentum[0] - exterior * 0.25 * dt;
        total += net;
        moment[1] += (z - 0.5) * net;
        moment[2] -= (y - 0.5) * net;
        work += net / dt * (1.0 + 2.0 * y + 3.0 * z);
    }
    assert!(moment[1] > 0.0 && moment[2] < 0.0);
    let resultant = rotate([total, 0.0, 0.0]);
    let angular = rotate(moment);
    for axis in 0..3 {
        close(load.resultant_impulse_n_s[axis], resultant[axis]);
        close(load.angular_impulse_n_m_s[axis], angular[axis]);
        close(
            load.exterior_impulse_n_s[axis],
            -AXES[0][axis] * exterior * dt,
        );
    }
    let displacement = (0..9)
        .map(|i| rotate([1.0 + (i % 3) as f64 + 1.5 * (i / 3) as f64, 0.0, 0.0]))
        .collect::<Vec<_>>();
    close(load.virtual_work_j(&displacement).unwrap(), work);
}

#[test]
fn reacting_molecular_support_retains_both_signed_stress_halves_and_refuses_pressure_alias() {
    let gas = mixture(source(3.5, 15.0));
    let grid = PlanarGrid::new(4, 1, 0.5, 1.0, 1.0).unwrap();
    let mut states = vec![0.0; 4 * 5];
    for q in states.chunks_exact_mut(5) {
        planar_state(&gas, &[0.7, 0.3], [1.0, 0.0], 1.0, q).unwrap();
    }
    let bytes = MixtureEuler2d::required_workspace_bytes(&gas, &grid).unwrap();
    let field = MixtureEuler2d::new(
        gas,
        grid,
        states,
        Boundaries2d { x: WALL, y: WALL },
        bytes,
        || false,
    )
    .unwrap();
    let law = ConstantTransport::new(
        field.gas().source().clone(),
        "manufactured-reacting-wall-traction",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 100.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.2),
            bulk_viscosity: DynViscosity::new(0.03),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.01),
        },
        vec![0.0, 0.02, 0.02, 0.0],
    )
    .unwrap();
    let source = PlanarTransport::new(field.gas().clone(), law, || false).unwrap();
    let mut operator =
        MolecularTransport2d::new(&field, source, MAX_MIXTURE_WORKSPACE_BYTES, || false).unwrap();
    let chemistry = CellChemistry::new(
        field.gas().clone(),
        kinetics(field.gas().source().clone()),
        chemistry(),
        || false,
    )
    .unwrap();
    let mut flow =
        ReactingEuler2d::new(field, chemistry, MAX_MIXTURE_WORKSPACE_BYTES, || false).unwrap();
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let controls = PlanarReactiveControl {
        transport: controls(),
        ..PlanarReactiveControl::default()
    };
    let bytes = flow.workspace_bytes() + operator.workspace_bytes();
    let report = flow
        .step_with_dissipation(
            &mut operator,
            Time::new(1e-4),
            controls,
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    assert!(report.chemistry.species_change_kg[1] > 0.0);
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    let left =
        WallTractionImpulse::from_reacting_euler2d(&flow, PlanarSide::XMin, 0, axes, 0.25).unwrap();
    let right =
        WallTractionImpulse::from_reacting_euler2d(&flow, PlanarSide::XMax, 0, axes, 0.0).unwrap();
    // The separating wall has tau_nn > p, while the opposite wall has a
    // compressive total traction. Neither molecular interval is pure pressure.
    assert!(left.gas_traction_impulse_pa_s() < 0.0);
    assert!(right.gas_traction_impulse_pa_s() > 0.0);
    for side in [PlanarSide::XMin, PlanarSide::XMax] {
        assert!(WallPressureImpulse::from_reacting_euler2d(&flow, side, 0, axes, 0.0).is_err());
    }
    let complete = flow
        .boundary_impulses(PlanarSide::XMin)
        .next()
        .unwrap()
        .momentum[0];
    let middle = flow
        .flow()
        .boundary_impulses(PlanarSide::XMin)
        .next()
        .unwrap()
        .momentum[0];
    assert!((complete - middle).abs() > 1e-5);
    close(left.gas_traction_impulse_pa_s(), -complete);
    assert_eq!(left.duration_s(), report.elapsed.value());
    assert_eq!(left.end_time_s(), flow.flow().time().value());

    let mut left_mesh = mesh();
    for point in &mut left_mesh.nodes {
        point[0] = 0.0;
    }
    let surface = ShockLoadSurface::new(
        &left_mesh,
        vec![selection(ShellFace::Negative)],
        budget(),
        || false,
    )
    .unwrap();
    let load = surface
        .project_tractions(&[left], [0.0, 0.5, 0.5], || false)
        .unwrap();
    let expected = complete + 0.25 * report.elapsed.value();
    close(load.resultant_impulse_n_s[0], expected);
    close(load.exterior_impulse_n_s[0], 0.25 * report.elapsed.value());
    let displacement = left_mesh
        .nodes
        .iter()
        .map(|p| [1.0 + 2.0 * p[1] + 3.0 * p[2], 0.0, 0.0])
        .collect::<Vec<_>>();
    close(
        load.virtual_work_j(&displacement).unwrap(),
        expected / report.elapsed.value() * 3.5,
    );

    let mut support = [0.0; 2];
    for side in [
        PlanarSide::XMin,
        PlanarSide::XMax,
        PlanarSide::YMin,
        PlanarSide::YMax,
    ] {
        for face in 0..flow.flow().grid().face_count(side) {
            let traction =
                WallTractionImpulse::from_reacting_euler2d(&flow, side, face, axes, 0.0).unwrap();
            let raw = flow.boundary_impulses(side).nth(face).unwrap();
            assert!(raw.species_masses.iter().all(|mass| *mass == 0.0));
            assert_eq!(raw.energy, 0.0);
            for axis in 0..2 {
                support[axis] += traction.outward_normal()[axis]
                    * traction.area_m2()
                    * traction.gas_traction_impulse_pa_s();
            }
        }
    }
    for axis in 0..2 {
        close(
            report.after.momentum[axis] - report.before.momentum[axis] + support[axis],
            0.0,
        );
        close(
            report.boundary_change.momentum[axis],
            report.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
        );
    }
    let states = flow.flow().states().to_vec();
    let molecular_report = *operator.last_step().unwrap();
    let mut polls = 0;
    assert!(
        flow.step_with_dissipation(
            &mut operator,
            Time::new(1e-4),
            controls,
            DissipationControl::default(),
            bytes,
            || {
                polls += 1;
                polls == 100
            },
        )
        .is_err()
    );
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(operator.last_step(), Some(&molecular_report));
    assert_eq!(
        WallTractionImpulse::from_reacting_euler2d(&flow, PlanarSide::XMin, 0, axes, 0.25).unwrap(),
        left
    );
}
