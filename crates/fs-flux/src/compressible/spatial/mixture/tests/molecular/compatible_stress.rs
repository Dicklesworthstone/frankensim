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
            Material::Composition => 1.0 + 299.0 * fractions[0],
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
    assert!(!prescribed.uses_compatible_molecular_stress());
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
