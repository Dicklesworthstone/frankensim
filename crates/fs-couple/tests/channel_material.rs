//! G1/G3/G4: original fluid-source integration, shared geometry arithmetic,
//! retained liquid momentum and atomic material continuation. Raw coefficient
//! fixtures test the mechanism; the IF97 case exercises the original water owner.

use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::bernoulli_aperture::slit_channel::{
    SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
    annular::{AnnularChannel, AnnularChannelSpec},
    circular::{CircularChannel, CircularChannelSpec},
    material::{MaterialChannel, MaterialChannelError, MaterialChannelGeometry},
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::transport::If97TransportState;
use std::cell::Cell;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap()
}

fn geometries(exit: SlitExit) -> [MaterialChannelGeometry; 3] {
    [
        MaterialChannelGeometry::Slit(SlitChannelSpec {
            length_m: 0.1,
            width_m: 0.006,
            height_m: 0.001,
            transverse_cells: 24,
            sample_period_s: 0.001,
            exit,
        }),
        MaterialChannelGeometry::Circular(CircularChannelSpec {
            length_m: 0.1,
            radius_m: 0.0005,
            radial_cells: 24,
            sample_period_s: 0.001,
            exit,
        }),
        MaterialChannelGeometry::Annular(AnnularChannelSpec {
            length_m: 0.1,
            inner_radius_m: 0.0002,
            outer_radius_m: 0.0005,
            radial_cells: 24,
            sample_period_s: 0.001,
            exit,
        }),
    ]
}

fn balance(step: SlitChannelStep) {
    let scale = step.pressure_work_j.abs()
        + step.wall_work_j.abs()
        + step.stored_delta_j.abs()
        + step.viscous_dissipation_j
        + step.numerical_dissipation_j
        + step.exit_work_j;
    assert!(step.energy_residual_j.abs() <= 3e-12 * scale, "{step:?}");
    assert!(step.viscous_dissipation_j >= 0.0 && step.numerical_dissipation_j >= 0.0);
}

// Deliberately not Clone: material ownership is independent of warm stepping.
#[derive(Debug)]
struct CountedSource {
    name: &'static str,
    coefficients: FluidTransportState,
    calls: Cell<usize>,
}
impl FluidTransport for CountedSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        Ok(self.coefficients)
    }
}

#[test]
fn g3_original_source_is_sampled_once_and_borrowed_without_clone() {
    for geometry in geometries(SlitExit::StaticPressure) {
        let source = CountedSource {
            name: "manufactured Newtonian liquid",
            coefficients: fluid(0.002),
            calls: Cell::new(0),
        };
        let mut channel = MaterialChannel::try_new(&source, geometry).unwrap();
        assert!(std::ptr::eq(*channel.source(), &source));
        assert_eq!(channel.transport(), source.coefficients);
        assert_eq!(channel.sample_period_s(), 0.001);
        assert_eq!(source.calls.get(), 1);
        let expected_volume = match geometry {
            MaterialChannelGeometry::Slit(spec) => spec.length_m * spec.width_m * spec.height_m,
            MaterialChannelGeometry::Circular(spec) => {
                std::f64::consts::PI * spec.radius_m.powi(2) * spec.length_m
            }
            MaterialChannelGeometry::Annular(spec) => {
                std::f64::consts::PI
                    * (spec.outer_radius_m - spec.inner_radius_m)
                    * (spec.outer_radius_m + spec.inner_radius_m)
                    * spec.length_m
            }
        };
        let mass: f64 = channel
            .continuous_operator()
            .unwrap()
            .cell_masses_kg()
            .iter()
            .sum();
        assert!((mass / (1000.0 * expected_volume) - 1.0).abs() < 5e-15);
        for pressure in [7.0, 7.0, -3.0, 0.0] {
            balance(
                channel
                    .step_with_wall_velocities(pressure, [0.0, 0.03])
                    .unwrap(),
            );
        }
        channel.pressure_load().unwrap();
        channel
            .preview_at_flow_with_wall_velocities(-1e-9, [0.0, 0.03])
            .unwrap();
        assert_eq!(
            source.calls.get(),
            1,
            "time steps and trials must not resample source"
        );
        assert_eq!(channel.source().name, "manufactured Newtonian liquid");
        let owned = CountedSource {
            name: "owned non-Clone source",
            coefficients: fluid(0.002),
            calls: Cell::new(0),
        };
        let mut owned_channel = MaterialChannel::try_new(owned, geometry).unwrap();
        balance(owned_channel.step(10.0).unwrap());
        assert_eq!(owned_channel.source().calls.get(), 1);
    }
}

#[test]
fn g3_viscosity_change_preserves_native_history_and_changes_future_diffusion() {
    for geometry in geometries(SlitExit::StaticPressure) {
        let mut old = MaterialChannel::try_new(
            CountedSource {
                name: "original liquid",
                coefficients: fluid(0.002),
                calls: Cell::new(0),
            },
            geometry,
        )
        .unwrap();
        for _ in 0..40 {
            old.step(50.0).unwrap();
        }
        let before = old.velocities_m_s().to_vec();
        let energy = old.stored_energy_j();
        let last = old.last_step();
        let mut next = old
            .with_source(CountedSource {
                name: "higher-viscosity liquid",
                coefficients: fluid(0.008),
                calls: Cell::new(0),
            })
            .unwrap();
        assert_eq!(next.source().calls.get(), 1);
        assert_eq!(old.source().calls.get(), 1);
        assert_eq!(next.source().name, "higher-viscosity liquid");
        assert_eq!(next.velocities_m_s(), before);
        assert_eq!(next.stored_energy_j(), energy);
        assert_eq!(
            next.last_step(),
            last,
            "a zero-duration source proposal fabricates no exchange"
        );
        assert_eq!(next.geometry(), geometry);
        assert_ne!(next.pressure_load().unwrap(), old.pressure_load().unwrap());
        // The generic source type may change while the original history survives.
        let mut same = old.with_source(fluid(0.002)).unwrap();
        for pressure in [4.0, -2.0, 0.0] {
            assert_eq!(same.step(pressure).unwrap(), old.step(pressure).unwrap());
        }
        // Compare both viscosities from exactly the same original field.
        let mut fast = same.with_source(fluid(0.008)).unwrap();
        for _ in 0..8 {
            balance(same.step(0.0).unwrap());
            balance(fast.step(0.0).unwrap());
        }
        assert!(fast.stored_energy_j() < 0.99 * same.stored_energy_j());
        assert_eq!(
            next.velocities_m_s(),
            before,
            "another continuation must not mutate a proposal"
        );
        balance(next.step(0.0).unwrap());
    }
}

#[derive(Debug)]
struct RefusedSource;
impl FluidTransport for RefusedSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Err(MaterialError::State {
            what: "original source rejected this state point".into(),
        })
    }
}

#[test]
fn g4_source_density_and_cancellation_refusals_preserve_retry() {
    for geometry in geometries(SlitExit::StaticPressure) {
        let mut channel = MaterialChannel::try_new(fluid(0.002), geometry).unwrap();
        for pressure in [12.0, -3.0, 8.0] {
            channel.step(pressure).unwrap();
        }
        let mut clean = channel.clone();
        let velocities = channel.velocities_m_s().to_vec();
        let last = channel.last_step();
        let energy = channel.stored_energy_j();
        assert_eq!(
            channel.with_source(RefusedSource).unwrap_err(),
            MaterialChannelError::Material(MaterialError::State {
                what: "original source rejected this state point".into(),
            }),
        );
        let different_density = FluidTransportState::try_new(1001.0, 4200.0, 0.6, 0.008).unwrap();
        assert_eq!(
            channel.with_source(different_density).unwrap_err(),
            MaterialChannelError::DensityChange {
                previous_kg_m3: 1000.0,
                proposed_kg_m3: 1001.0,
            },
        );
        let mut successful_polls = 0;
        channel
            .with_source_controlled(fluid(0.008), &mut || {
                successful_polls += 1;
                false
            })
            .unwrap();
        assert!(successful_polls >= 4);
        for cancel_at in [1, successful_polls] {
            let mut polls = 0;
            let refused = channel.with_source_controlled(fluid(0.008), &mut || {
                polls += 1;
                polls == cancel_at
            });
            assert!(matches!(
                refused,
                Err(MaterialChannelError::Channel(
                    AcousticRealizeError::Cancelled
                ))
            ));
            assert_eq!(channel.velocities_m_s(), velocities);
            assert_eq!(channel.last_step(), last);
            assert_eq!(channel.stored_energy_j(), energy);
            assert_eq!(channel.transport(), fluid(0.002));
        }
        for pressure in [0.0, -8.0, 1.0] {
            assert_eq!(
                channel.step(pressure).unwrap(),
                clean.step(pressure).unwrap()
            );
        }
    }
}

#[test]
fn g3_material_flow_preview_uses_original_pressure_and_wall_work() {
    for geometry in geometries(SlitExit::StaticPressure) {
        let mut channel = MaterialChannel::try_new(fluid(0.002), geometry).unwrap();
        for _ in 0..12 {
            channel.step(10.0).unwrap();
        }
        let before = channel.velocities_m_s().to_vec();
        let last = channel.last_step();
        let walls = [0.0, 0.05];
        let flow = -1e-9;
        let load = channel.pressure_load_with_wall_velocities(walls).unwrap();
        let trial = channel
            .preview_at_flow_with_wall_velocities(flow, walls)
            .unwrap();
        let step = trial.last_step();
        balance(step);
        assert!((step.flow_m3_s - flow).abs() < 2e-21);
        assert_eq!(
            step.pressure_work_j,
            load.pressure_drop_pa(flow).unwrap() * step.flow_m3_s * channel.sample_period_s()
        );
        assert!(step.wall_work_j > 0.0);
        assert_eq!(channel.velocities_m_s(), before);
        assert_eq!(channel.last_step(), last);
        assert!(
            channel
                .step_with_wall_velocities(0.0, [0.0, f64::NAN])
                .is_err()
        );
        assert_eq!(channel.velocities_m_s(), before);
        assert_eq!(channel.last_step(), last);
    }
    for geometry in geometries(SlitExit::UniformKineticHead) {
        let mut channel = MaterialChannel::try_new(fluid(0.002), geometry).unwrap();
        balance(channel.step_with_wall_velocities(0.0, [0.0, 0.05]).unwrap());
        assert!(channel.pressure_load().is_err());
        assert!(channel.continuous_operator().is_err());
        assert!(channel.preview_at_flow(1e-9).is_err());
    }
}

enum OriginalGasChannel {
    Slit(SlitChannel),
    Circular(CircularChannel),
    Annular(AnnularChannel),
}
impl OriginalGasChannel {
    fn new(gas: &GasState, geometry: MaterialChannelGeometry) -> Self {
        match geometry {
            MaterialChannelGeometry::Slit(spec) => {
                Self::Slit(SlitChannel::from_gas(gas, spec).unwrap())
            }
            MaterialChannelGeometry::Circular(spec) => {
                Self::Circular(CircularChannel::from_gas(gas, spec).unwrap())
            }
            MaterialChannelGeometry::Annular(spec) => {
                Self::Annular(AnnularChannel::from_gas(gas, spec).unwrap())
            }
        }
    }
    fn step(&mut self, pressure: f64, walls: [f64; 2]) -> SlitChannelStep {
        match self {
            Self::Slit(owner) => owner.step_with_wall_velocities(pressure, walls).unwrap(),
            Self::Circular(owner) => owner.step_with_wall_velocities(pressure, walls).unwrap(),
            Self::Annular(owner) => owner.step_with_wall_velocities(pressure, walls).unwrap(),
        }
    }
    fn velocities(&self) -> &[f64] {
        match self {
            Self::Slit(owner) => owner.velocities_m_s(),
            Self::Circular(owner) => owner.velocities_m_s(),
            Self::Annular(owner) => owner.velocities_m_s(),
        }
    }
}

#[test]
fn g3_gas_adapter_and_material_adapter_remain_bit_identical() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    for exit in [SlitExit::StaticPressure, SlitExit::UniformKineticHead] {
        for geometry in geometries(exit) {
            let mut original = OriginalGasChannel::new(&gas, geometry);
            let mut material = MaterialChannel::try_new(gas, geometry).unwrap();
            for (pressure, upper_wall) in [(8.0, 0.0), (8.0, 0.04), (-3.0, -0.02), (0.0, 0.01)] {
                let walls = [0.0, upper_wall];
                assert_eq!(
                    material.step_with_wall_velocities(pressure, walls).unwrap(),
                    original.step(pressure, walls)
                );
                assert_eq!(material.velocities_m_s(), original.velocities());
            }
        }
    }
}

#[test]
fn g1_real_if97_liquid_uses_original_source_and_converges_to_poiseuille() {
    // This validates channel/source integration and the independent Poiseuille
    // limit; the industrial water constitutive law has its own source checks.
    let water = If97TransportState::at_pressure_temperature(
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
        || false,
    )
    .unwrap();
    let transport = water.transport_state().unwrap();
    let (length, radius, pressure) = (0.1, 0.0005_f64, 10.0);
    let exact = pressure * radius.powi(2) / (8.0 * transport.dynamic_viscosity_pa_s() * length);
    let mut errors = Vec::new();
    for radial_cells in [16, 32, 64] {
        let spec = CircularChannelSpec {
            length_m: length,
            radius_m: radius,
            radial_cells,
            sample_period_s: 0.005,
            exit: SlitExit::StaticPressure,
        };
        let mut channel =
            MaterialChannel::try_new(water, MaterialChannelGeometry::Circular(spec)).unwrap();
        assert_eq!(channel.source(), &water);
        assert_eq!(channel.transport(), transport);
        assert_eq!(
            channel.source().thermodynamics().temperature().value(),
            300.0
        );
        let steady = channel.steady_bulk_speed_m_s(pressure).unwrap();
        errors.push((steady - exact).abs());
        if radial_cells == 64 {
            for _ in 0..500 {
                balance(channel.step(pressure).unwrap());
            }
            assert!((channel.last_step().bulk_speed_m_s / steady - 1.0).abs() < 2e-11);
        }
    }
    assert!(
        errors[0] > exact * 1e-4,
        "spatial refinement must be observable"
    );
    for pair in errors.windows(2) {
        assert!(
            (1.95..2.05).contains(&(pair[0] / pair[1]).log2()),
            "{errors:?}"
        );
    }
}
