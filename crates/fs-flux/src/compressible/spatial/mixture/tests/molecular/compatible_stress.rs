use super::*;
use crate::compressible::spatial::mixture::molecular::{
    AdaptiveDissipationControl, PressureOutletPolicy, StabilizedDissipationControl,
};
use fs_gas::multicomponent::dissipation::{SoretDufourTransport, TransportLaw};
use fs_qty::Pressure;
use fs_thermochem::{
    MassSpecificEnergyQuantityV1, equilibrium::composition::IdealGasEquilibriumV1,
};

#[derive(Clone, Copy)]
enum Material {
    Constant,
    Composition,
    WallComposition,
    CellPeaks,
    FacePeak,
}

struct VariableTransport {
    native: ConstantTransport,
    material: Material,
}
impl TransportLaw for VariableTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.native.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value =
            self.native
                .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        let midpoint = (1.25..1.75).contains(&temperature.value());
        value.dynamic_viscosity = DynViscosity::new(match self.material {
            Material::Constant => 1.0,
            Material::Composition | Material::WallComposition => 1.0 + 299.0 * fractions[0],
            Material::CellPeaks => {
                if midpoint {
                    1.0
                } else {
                    300.0
                }
            }
            Material::FacePeak => {
                if midpoint {
                    1e30
                } else {
                    1.0
                }
            }
        });
        if matches!(self.material, Material::WallComposition) {
            value.bulk_viscosity = DynViscosity::new(0.1 + 0.2 * fractions[0]);
        }
        Ok(value)
    }
}
type OriginalTransport = SoretDufourTransport<VariableTransport>;

fn source(flow: &MixtureEuler3d, material: Material) -> SpatialTransport<OriginalTransport> {
    let heats = flow
        .gas()
        .source()
        .models()
        .iter()
        .enumerate()
        .map(|(i, model)| {
            (
                model.species().clone(),
                MassSpecificEnergyQuantityV1::new([0.2, -0.3, 0.0][i]),
            )
        })
        .collect();
    let native = VariableTransport {
        native: law(flow.gas(), 0.03, 1.0, 0.0, 0.2),
        material,
    };
    let cross = SoretDufourTransport::new(
        native,
        "compatible-mixture-source",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        heats,
    )
    .unwrap();
    SpatialTransport::new(flow.gas().clone(), cross, || false).unwrap()
}

fn operator(flow: &MixtureEuler3d, material: Material) -> MolecularTransport3d<OriginalTransport> {
    let native = source(flow, material);
    let base = MolecularTransport3d::required_workspace_bytes(flow, &native).unwrap();
    let bytes = MolecularTransport3d::required_adaptive_workspace_bytes(flow, &native).unwrap();
    let op = MolecularTransport3d::new(flow, native, bytes, || false).unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), base);
    op
}

fn vortex() -> MixtureEuler3d {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    field(
        equal_mass_gas(),
        CartesianGrid::new(4, 4, 1, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            let i = x[0] as usize;
            let j = x[1] as usize;
            let species_a = if i == 0 && j == 0 { 1.0 } else { 0.0 };
            (
                [species_a, 1.0 - species_a, 0.0],
                [-a[i] * b[j], a[j] * b[i], 0.0],
                1.0,
            )
        },
    )
}

fn kinetic(flow: &MixtureEuler3d) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            0.5 * q.momentum.iter().map(|p| p * p).sum::<f64>()
                / q.partial_densities.iter().sum::<f64>()
                * flow.grid().cell_volume()
        })
        .sum()
}

#[test]
fn g1_original_composition_viscosity_changes_curl_without_creating_kinetic_energy() {
    for integrator in 0..3 {
        let mut flow = vortex();
        let mut op = operator(&flow, Material::Composition);
        assert!(op.uses_compatible_molecular_stress());
        let before = flow.states().to_vec();
        let initial_kinetic = kinetic(&flow);
        let mut curl = flow
            .vorticity_workspace(flow.required_vorticity_workspace_bytes().unwrap())
            .unwrap();
        let initial_curl = flow.vorticity(&mut curl, || false).unwrap().to_vec();
        let dt = Time::new(1e-6);
        let report = match integrator {
            0 => op
                .advance(&mut flow, dt, DissipationControl::default(), || false)
                .unwrap(),
            1 => {
                op.advance_stabilized(
                    &mut flow,
                    dt,
                    StabilizedDissipationControl::default(),
                    || false,
                )
                .unwrap()
                .transport
            }
            _ => {
                op.advance_adaptive(&mut flow, dt, AdaptiveDissipationControl::default(), || {
                    false
                })
                .unwrap()
                .stabilized
                .transport
            }
        };
        check(report);
        assert!(kinetic(&flow) < initial_kinetic);
        assert_ne!(flow.cell(0).unwrap().partial_densities, &before[..3]);
        for i in 0..flow.cell_count() {
            assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0);
        }
        let changed_curl = flow.vorticity(&mut curl, || false).unwrap();
        assert!(
            initial_curl
                .iter()
                .zip(changed_curl)
                .any(|(a, b)| a.vorticity != b.vorticity)
        );
        assert_eq!(flow.time().value(), 0.0);
    }
}

fn alternating() -> MixtureEuler3d {
    field(
        equal_mass_gas(),
        CartesianGrid::new(2, 2, 2, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            let high = (x[0] as usize + x[1] as usize + x[2] as usize) % 2 == 0;
            let sign = if high { 1.0 } else { -1.0 };
            (
                [0.5 + 0.1 * sign, 0.5 - 0.1 * sign, 0.0],
                [1.0 + 0.1 * sign, -0.7 - 0.2 * sign, 0.3 + 0.05 * sign],
                if high { 2.0 } else { 1.0 },
            )
        },
    )
}

#[test]
fn g0_original_scalar_transport_survives_extreme_unused_face_stress_and_actual_cell_rate_bounds() {
    let mut flow = alternating();
    let mut cells = operator(&flow, Material::CellPeaks);
    let report = cells
        .advance(
            &mut flow,
            Time::new(1e-6),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    check(report);
    // Midpoint mu is1, while every cell mu is300: the actual mechanical bound
    // is at least 4*(20/3)*300*(1/dx²+1/dy²+1/dz²), before multiplying bydt.
    assert!(report.maximum_parabolic_number >= 0.024 * (1.0 - 2e-12));

    let mut reference = alternating();
    let mut peaked = alternating();
    let mut base = operator(&reference, Material::Constant);
    let mut extreme = operator(&peaked, Material::FacePeak);
    let expected = base
        .advance(
            &mut reference,
            Time::new(1e-5),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let actual = extreme
        .advance(
            &mut peaked,
            Time::new(1e-5),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    check(actual);
    close(
        actual.maximum_parabolic_number,
        expected.maximum_parabolic_number,
        2e-11,
    );
    for (&a, &b) in peaked.states().iter().zip(reference.states()) {
        close(a, b, 3e-11);
    }
    // This equality includes original Maxwell--Stefan/Soret species flux,
    // referenced partial enthalpies, Fourier heat and reciprocal Dufour heat.
    assert_eq!(actual.accepted_steps, expected.accepted_steps);
}

#[test]
fn g4_cell_sampling_budget_boundary_scope_and_late_cancellation_remain_transactional() {
    let mut reference = vortex();
    let native = source(&reference, Material::Composition);
    let bytes = MolecularTransport3d::required_workspace_bytes(&reference, &native).unwrap();
    assert!(MolecularTransport3d::new(&reference, native, bytes - 1, || false).is_err());
    let mut base = operator(&reference, Material::Composition);
    let mut calls = 0;
    let expected = base
        .advance(
            &mut reference,
            Time::new(1e-6),
            DissipationControl::default(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = vortex();
    let before = flow.states().to_vec();
    let time = flow.time();
    let mut op = operator(&flow, Material::Composition);
    let mut seen = 0;
    assert_eq!(
        op.advance(
            &mut flow,
            Time::new(1e-6),
            DissipationControl::default(),
            || {
                seen += 1;
                seen == calls
            }
        )
        .unwrap_err(),
        DissipationError::Cancelled
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), time);
    assert!(op.last_step().is_none());
    let actual = op
        .advance(
            &mut flow,
            Time::new(1e-6),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), reference.states());

    let grid = CartesianGrid::new(2, 2, 2, 1.0, 1.0, 1.0).unwrap();
    let wall = field(
        equal_mass_gas(),
        grid,
        Boundaries3d {
            x: WALLS,
            ..PERIODIC
        },
        |_| ([0.5, 0.5, 0.0], [0.0; 3], 1.0),
    );
    let free = operator(&wall, Material::Constant);
    assert!(free.uses_compatible_molecular_stress());
    let prescribed = free
        .with_wall(
            BoundarySide::XMin,
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
            },
        )
        .unwrap();
    assert!(prescribed.uses_compatible_molecular_stress());
    let finite = prescribed
        .with_wall(
            BoundarySide::XMin,
            MolecularWall3d::NavierSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
                slip_length_m: 0.1,
                fluid_heat_fraction: 0.5,
            },
        )
        .unwrap();
    assert!(!finite.uses_compatible_molecular_stress());
    let outlet = field(
        equal_mass_gas(),
        grid,
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::PressureOutlet { pressure_pa: 1.0 },
                upper: BoundaryCondition::PressureOutlet { pressure_pa: 1.0 },
            },
            ..PERIODIC
        },
        |_| ([0.5, 0.5, 0.0], [0.0; 3], 1.0),
    );
    let native = source(&outlet, Material::Constant);
    let bytes = MolecularTransport3d::required_workspace_bytes(&outlet, &native).unwrap();
    let open = MolecularTransport3d::new_with_pressure_outlet_policy(
        &outlet,
        native,
        PressureOutletPolicy::ZeroDiffusiveFlux,
        bytes,
        || false,
    )
    .unwrap();
    assert!(open.uses_compatible_molecular_stress());
}

#[test]
fn original_composition_shear_and_bulk_obey_the_prescribed_wall_work_budget_on_every_axis() {
    let profile = [0.2, 0.3, 0.1, 0.3];
    let wave = [0.0, 1.0, 0.0, -1.0];
    let duration = 1e-8;
    for inactive in 0..3 {
        let a = (inactive + 1) % 3;
        let b = (inactive + 2) % 3;
        let shape: [usize; 3] = std::array::from_fn(|axis| if axis == inactive { 1 } else { 4 });
        let axes: [AxisBoundaries; 3] = std::array::from_fn(|axis| {
            if axis == inactive {
                AxisBoundaries::Periodic
            } else {
                WALLS
            }
        });
        for (lid_speed, dissipation) in [(0.0, 10.267166666666667), (0.4, 10.107166666666667)] {
            let mut flow = field(
                equal_mass_gas(),
                CartesianGrid::new(shape[0], shape[1], shape[2], 1.0, 1.0, 1.0).unwrap(),
                Boundaries3d {
                    x: axes[0],
                    y: axes[1],
                    z: axes[2],
                },
                |x| {
                    let (i, j) = (x[a] as usize, x[b] as usize);
                    let species_a = if i == 0 && j == 0 { 1.0 } else { 0.0 };
                    let mut velocity = [0.0; 3];
                    velocity[a] = -profile[i] * wave[j];
                    velocity[b] = profile[j] * wave[i];
                    ([species_a, 1.0 - species_a, 0.0], velocity, 1.0)
                },
            );
            let mut op = operator(&flow, Material::WallComposition);
            let lid = BoundarySide::ALL[2 * b + 1];
            for side in BoundarySide::ALL {
                if side.axis() == inactive {
                    continue;
                }
                let mut velocity = [0.0; 3];
                if side == lid {
                    velocity[a] = lid_speed;
                }
                op = op
                    .with_wall(
                        side,
                        MolecularWall3d::NoSlip {
                            velocity,
                            thermal: WallThermalCondition::Adiabatic,
                        },
                    )
                    .unwrap();
            }
            assert!(op.uses_compatible_molecular_stress());
            let before = kinetic(&flow);
            let cells = flow.states().to_vec();
            let report = op
                .advance(
                    &mut flow,
                    Time::new(duration),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap();
            check(report);
            let mut outward_work = 0.0;
            for side in BoundarySide::ALL {
                if side.axis() == inactive {
                    continue;
                }
                let speed = if side == lid { lid_speed } else { 0.0 };
                for load in op.boundary_impulses(side) {
                    assert_eq!(load.species_masses, [0.0; 3]);
                    close(load.energy, speed * load.momentum[a], 2e-17);
                    outward_work += load.energy;
                }
            }
            // Independently summed cell form: mu=1+299*x_A, kappa=.1+.2*x_A.
            // Rotating this extruded fixture checks derivative-first source
            // caches and normal/tangent momentum permutations on all six sides.
            close(
                (kinetic(&flow) - before + outward_work) / duration,
                -dissipation,
                0.003,
            );
            assert_ne!(flow.cell(0).unwrap().partial_densities, &cells[..3]);
            for i in 0..flow.cell_count() {
                assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0);
            }
            assert_eq!(flow.time().value(), 0.0);
        }
    }
}

#[test]
fn heated_no_slip_keeps_source_soret_dufour_and_enthalpy_with_large_unused_wall_viscosity() {
    let initial = || {
        field(
            equal_mass_gas(),
            CartesianGrid::new(2, 2, 2, 1.0, 1.0, 1.0).unwrap(),
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let y = 0.4 + 0.03 * (x[0] + x[1]);
                ([y, 1.0 - y, 0.0], [0.0; 3], 1.0)
            },
        )
    };
    let mut actual = initial();
    let mut expected = initial();
    let before = actual.states().to_vec();
    let mut spike = operator(&actual, Material::FacePeak);
    let mut constant = operator(&expected, Material::Constant);
    for side in BoundarySide::ALL {
        let wall = MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
        };
        spike = spike.with_wall(side, wall).unwrap();
        constant = constant.with_wall(side, wall).unwrap();
    }
    let report = spike
        .advance_stabilized(
            &mut actual,
            Time::new(1e-5),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap()
        .transport;
    check(report);
    check(
        constant
            .advance_stabilized(
                &mut expected,
                Time::new(1e-5),
                StabilizedDissipationControl::default(),
                || false,
            )
            .unwrap()
            .transport,
    );
    assert_eq!(actual.states(), expected.states());
    assert_ne!(actual.states(), before);
    assert!(report.maximum_parabolic_number < 0.002);
    for side in BoundarySide::ALL {
        for (left, right) in spike
            .boundary_impulses(side)
            .zip(constant.boundary_impulses(side))
        {
            assert_eq!(left.species_masses, right.species_masses);
            assert_eq!(left.momentum, right.momentum);
            assert_eq!(left.energy, right.energy);
        }
    }
}
