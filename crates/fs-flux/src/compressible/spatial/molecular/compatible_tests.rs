use super::*;
use fs_gas::multicomponent::dissipation::DissipationError;

#[derive(Clone, Copy)]
enum Profile {
    Constant,
    Linear,
    CellPeaks,
    FacePeak,
}

struct SampledViscosity {
    native: ConstantTransport,
    profile: Profile,
}

impl TransportLaw for SampledViscosity {
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
        let t = temperature.value();
        let midpoint = (1.25..1.75).contains(&t);
        value.dynamic_viscosity = DynViscosity::new(match self.profile {
            Profile::Constant => 1.0,
            Profile::Linear => 1.0 + 299.0 * (t - 1.0).max(0.0),
            Profile::CellPeaks => {
                if midpoint {
                    1.0
                } else {
                    300.0
                }
            }
            Profile::FacePeak => {
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

fn source(gas: &FrozenGas, profile: Profile, conductivity: f64) -> SampledViscosity {
    SampledViscosity {
        native: law(gas, 1.0, 0.0, conductivity),
        profile,
    }
}

fn operator(
    flow: &Euler3d,
    profile: Profile,
    conductivity: f64,
    modeled: bool,
) -> MolecularTransport3d<SampledViscosity> {
    let native = source(flow.gas(), profile, conductivity);
    #[cfg(feature = "experimental-sgs")]
    if modeled {
        use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};
        return MolecularTransport3d::with_sgs(
            flow,
            native,
            SubgridConfig {
                model: SgsModel::Smagorinsky { cs: 0.18 },
                filter_width_factor: 1.0,
                turbulent_prandtl: 0.8,
                maximum_eddy_to_molecular: 20.0,
                usage: SgsUse::Exploration,
            },
            || false,
        )
        .unwrap();
    }
    assert!(!modeled);
    MolecularTransport3d::new(flow, native, || false).unwrap()
}

fn prepare<L: TransportLaw>(op: &mut MolecularTransport3d<L>, flow: &Euler3d) -> f64 {
    spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        None,
        flow.states(),
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.viscosity,
        &mut op.work.flux,
        &mut op.work.binary,
        #[cfg(feature = "experimental-sgs")]
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap()
}

fn vortex() -> Euler3d {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    field(
        gas(0.0),
        CartesianGrid::new(4, 4, 1, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            let i = x[0] as usize;
            let j = x[1] as usize;
            [
                1.0,
                -a[i] * b[j],
                a[j] * b[i],
                0.0,
                if i == 0 && j == 0 { 2.0 } else { 1.0 },
            ]
        },
    )
}

#[test]
fn g1_variable_molecular_viscosity_dissipates_the_original_planar_and_spatial_vortex() {
    // The old midpoint-coefficient cross stress gives dK/dt=+3.613333... .
    // Original cell sampling and compatible mechanics give -2.366666... .
    let mut space = vortex();
    let source_gas = space.gas().clone();
    let cells = space
        .states()
        .iter()
        .map(|q| {
            let p = q.primitive(&source_gas, || false).unwrap();
            planar::State::from_primitive(
                &source_gas,
                p.density,
                [p.velocity[0], p.velocity[1]],
                p.pressure,
            )
            .unwrap()
        })
        .collect();
    let mut plane = planar::Euler2d::new(
        source_gas,
        planar::CartesianGrid::new(4, 4, 1.0, 1.0, 1.0).unwrap(),
        cells,
        planar::Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    let mut planar_op = planar::molecular::MolecularTransport2d::new(
        &plane,
        source(plane.gas(), Profile::Linear, 0.0),
        || false,
    )
    .unwrap();
    let mut spatial_op = operator(&space, Profile::Linear, 0.0, false);
    let before = partition(&space);
    let duration = Time::new(1e-6);
    let planar_report = planar_op
        .advance(&mut plane, duration, MolecularControl::default(), || false)
        .unwrap();
    let spatial_report = spatial_op
        .advance(&mut space, duration, MolecularControl::default(), || false)
        .unwrap();
    assert_eq!(planar_report.accepted_steps, 1);
    assert_eq!(spatial_report.accepted_steps, 1);
    let after = partition(&space);
    assert!(after.0 < before.0);
    assert!(after.1 > before.1);
    close(after.0 + after.1, before.0 + before.1, 3e-12);
    close(spatial_report.conservation_residual.energy, 0.0, 3e-12);
    close(planar_report.conservation_residual.energy, 0.0, 3e-12);
    for (two, three) in plane.states().iter().zip(space.states()) {
        close(two.momentum[0], three.momentum[0], 3e-12);
        close(two.momentum[1], three.momentum[1], 3e-12);
        close(two.energy, three.energy, 3e-12);
        assert_eq!(three.momentum[2], 0.0);
    }
}

fn alternating() -> Euler3d {
    field(
        gas(0.0),
        CartesianGrid::new(2, 2, 2, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            let high = (x[0] as usize + x[1] as usize + x[2] as usize) % 2 == 0;
            let sign = if high { 1.0 } else { -1.0 };
            [
                1.0,
                0.1 * sign,
                -0.2 * sign,
                0.05 * sign,
                if high { 2.0 } else { 1.0 },
            ]
        },
    )
}

#[test]
fn g0_nonlinear_material_bounds_use_cells_and_extreme_face_peaks_keep_stress_heat_and_sgs() {
    let flow = alternating();
    let mut cells = operator(&flow, Profile::CellPeaks, 0.2, false);
    let rate = prepare(&mut cells, &flow);
    // All old midpoint samples are one, but every actual cell sample is300.
    close(rate, 4.0 * (20.0 / 3.0) * 300.0 * 3.0, 2e-10);
    for modeled in [false, true] {
        #[cfg(not(feature = "experimental-sgs"))]
        if modeled {
            continue;
        }
        let mut reference = operator(&flow, Profile::Constant, 0.2, modeled);
        let mut peaked = operator(&flow, Profile::FacePeak, 0.2, modeled);
        let reference_rate = prepare(&mut reference, &flow);
        let actual_rate = prepare(&mut peaked, &flow);
        close(actual_rate, reference_rate, 2e-12);
        for axis in 0..3 {
            for (actual, expected) in peaked.work.flux[axis]
                .iter()
                .zip(&reference.work.flux[axis])
            {
                close(actual.normal_momentum, expected.normal_momentum, 2e-13);
                for tangent in 0..2 {
                    close(
                        actual.transverse_momentum[tangent],
                        expected.transverse_momentum[tangent],
                        2e-13,
                    );
                }
                close(actual.energy, expected.energy, 2e-13);
            }
        }
    }
}

#[test]
fn g0_compatible_scope_preserves_prescribed_boundaries_and_constant_law_bits() {
    let flow = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap(),
        WALLS,
        |_| [1.0, 0.0, 0.0, 0.0, 1.0],
    );
    let free = operator(&flow, Profile::Linear, 0.0, false);
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
    assert!(
        !operator(&flow, Profile::Linear, 0.0, false)
            .with_reservoir_exchange()
            .uses_compatible_molecular_stress()
    );

    let mut compatible_flow = vortex();
    let mut native_flow = vortex();
    let mut compatible = operator(&compatible_flow, Profile::Constant, 0.02, false);
    // Empty bath exchange on a periodic grid selects the retained native
    // mechanical path without changing a physical face or its source data.
    let mut native =
        operator(&native_flow, Profile::Constant, 0.02, false).with_reservoir_exchange();
    let expected = native
        .advance(
            &mut native_flow,
            Time::new(1e-5),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let actual = compatible
        .advance(
            &mut compatible_flow,
            Time::new(1e-5),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(compatible_flow.states(), native_flow.states());
}

#[test]
fn g4_variable_material_late_cancellation_preserves_flow_and_replays_native_sources() {
    let mut reference = vortex();
    let mut reference_op = operator(&reference, Profile::Linear, 0.01, false);
    let mut calls = 0;
    let expected = reference_op
        .advance(
            &mut reference,
            Time::new(1e-6),
            MolecularControl::default(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = vortex();
    let original = flow.states().to_vec();
    let time = flow.time();
    let mut op = operator(&flow, Profile::Linear, 0.01, false);
    let mut seen = 0;
    assert_eq!(
        op.advance(
            &mut flow,
            Time::new(1e-6),
            MolecularControl::default(),
            || {
                seen += 1;
                seen == calls
            }
        )
        .unwrap_err(),
        MolecularError::Cancelled
    );
    assert_eq!(flow.states(), original);
    assert_eq!(flow.time(), time);
    for side in BoundarySide::ALL {
        assert!(
            op.boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
    let actual = op
        .advance(
            &mut flow,
            Time::new(1e-6),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), reference.states());
}
