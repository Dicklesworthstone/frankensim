//! Actual gas-source epochs on the original material and native thermal path.
use fs_exec::CancelGate;
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::boussinesq::{BoussinesqConfig, BoussinesqEndpointStep};
use fs_lbm::material::{
    MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3,
    plan_scaling_from_material,
};
use fs_lbm::thermal::sources::transport::RelaxationField;
use fs_lbm::thermal::sources::{HeatSourceLimits, ThermalSourceLimits};
use fs_material::gas::{GasSpec, GasState};
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};
use std::cell::Cell;

fn setup() -> (GasSpec, FluidTransportState, MaterialSlab, MaterialSlab3) {
    let spec = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&spec, 305.0, 100_000.0).unwrap();
    let retained = reference.transport_state().unwrap();
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.001, 0.02).unwrap();
    let mut two = scaling
        .thermal_slab(MaterialThermalSlab {
            nx: 2,
            ny: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap();
    let mut three = scaling
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 1,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap();
    // A physical non-equilibrium stress makes the momentum response depend
    // on the originally sampled molecular viscosity during this first step.
    two.lattice.grid.f[6][0] -= 0.002;
    two.lattice.grid.f[6][1] += 0.001;
    two.lattice.grid.f[6][3] += 0.001;
    three.lattice.grid.f[6][0] -= 0.002;
    three.lattice.grid.f[6][1] += 0.001;
    three.lattice.grid.f[6][2] += 0.001;
    (spec, retained, two, three)
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: 16,
        temperature_range_k: [290.0, 330.0],
        balance_tolerance_j: 1e-8,
    }
}
fn config<const D: usize>(gravity: [f64; D]) -> BoussinesqConfig<D> {
    BoussinesqConfig {
        gravity_m_s2: gravity,
        max_relative_density_departure: 0.05,
        max_relative_specific_heat_departure: 0.01,
    }
}
fn power() -> [f64; 16] {
    let mut values = [0.0; 16];
    values[2..14].fill(1e6);
    values
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[derive(Debug)]
struct GasTicket {
    gas: GasState,
    temperature_k: f64,
    query: usize,
    reads: Cell<usize>,
}
impl FluidTransport for GasTicket {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        assert_eq!(
            self.reads.get(),
            1,
            "an original source response is resolved once"
        );
        self.gas.transport_state()
    }
}
fn ticket(spec: &GasSpec, temperature_k: f64, query: usize) -> Result<GasTicket, MaterialError> {
    Ok(GasTicket {
        gas: GasState::try_new(spec, temperature_k, 100_000.0)?,
        temperature_k,
        query,
        reads: Cell::new(0),
    })
}

fn check_report<const D: usize>(
    report: &BoussinesqEndpointStep<GasTicket>,
    queries: &[(usize, u64)],
    accepted: &[f64],
    force: &[[f64; D]],
    external: [f64; D],
    gravity: [f64; D],
    retained: FluidTransportState,
) {
    assert_eq!(report.transport.applied.len(), 12);
    assert_eq!(report.force_endpoint.len(), 12);
    assert_eq!(report.transport.endpoint.len(), 12);
    assert_eq!(queries.len(), 36);
    for (query, sample) in report
        .transport
        .applied
        .iter()
        .chain(&report.force_endpoint)
        .chain(&report.transport.endpoint)
        .enumerate()
    {
        assert_eq!(
            queries[query],
            (sample.cell, sample.temperature_k.to_bits())
        );
        let original = sample.material.material();
        assert_eq!(original.query, query);
        assert_eq!(
            original.temperature_k.to_bits(),
            sample.temperature_k.to_bits()
        );
        assert_eq!(original.reads.get(), 1);
        assert_eq!(sample.transport.density_kg_m3(), retained.density_kg_m3());
        assert_eq!(
            sample.transport.specific_heat_j_kg_k(),
            retained.specific_heat_j_kg_k()
        );
        assert_eq!(
            sample.transport,
            *sample.material.reference_storage_transport()
        );
        let raw = sample.material.resolved_transport();
        close(
            raw.density_kg_m3(),
            retained.density_kg_m3() * 305.0 / sample.temperature_k,
            1e-12,
        );
        close(
            sample.momentum_relaxation,
            0.5 + 3.0 * raw.dynamic_viscosity_pa_s() / retained.density_kg_m3() * 10.0,
            1e-14,
        );
    }
    let mut nonlinear_difference = 0.0_f64;
    for ((initial, predicted), endpoint) in report
        .transport
        .applied
        .iter()
        .zip(&report.force_endpoint)
        .zip(&report.transport.endpoint)
    {
        assert_eq!(
            endpoint.temperature_k.to_bits(),
            accepted[endpoint.cell].to_bits()
        );
        assert_eq!(endpoint.material.material().query, 24 + endpoint.cell - 2);
        // Independent constant-pressure gas law: average the TWO source
        // densities, not rho at the mean temperature or a fitted beta.
        let rho0 = retained.density_kg_m3();
        let departure = 0.5 * rho0 * (305.0 / initial.temperature_k - 1.0)
            + 0.5 * rho0 * (305.0 / predicted.temperature_k - 1.0);
        let midpoint_density =
            rho0 * (305.0 / (0.5 * initial.temperature_k + 0.5 * predicted.temperature_k) - 1.0);
        nonlinear_difference = nonlinear_difference.max((departure - midpoint_density).abs());
        for axis in 0..D {
            let expected =
                (external[axis] + departure * gravity[axis]) / rho0 * 0.001_f64.powi(2) / 0.01;
            close(force[initial.cell][axis], expected, 2e-16);
        }
    }
    assert!(nonlinear_difference > 1e-7);
    close(
        report.transport.heat.source_j,
        12.0 * 1e6 * 1e-6 * 0.001,
        1e-14,
    );
    assert!(report.transport.heat.defect_j.abs() < 1e-10);
}

#[test]
fn nonlinear_gas_density_has_distinct_force_and_exact_accepted_source_epochs() {
    let (spec, retained, mut two, mut three) = setup();
    let (_, _, mut control2, mut control3) = setup();
    let gate = CancelGate::new();
    let external2 = [1.0, -0.5];
    let external3 = [1.0, -0.5, 0.25];
    let gravity2 = [-9.81, 0.0];
    let gravity3 = [3.0, -4.0, -8.0];
    two.lattice.gbeta = 0.25;
    three.lattice.gbeta = 0.25;
    let mut queries2 = Vec::new();
    let report2 = two
        .step_with_boussinesq_material_transport_endpoint_averaged(
            &power(),
            &[external2; 16],
            0.01,
            limits(),
            config(gravity2),
            &gate,
            |cell, temperature, reference| {
                assert_eq!(*reference, retained);
                let query = queries2.len();
                queries2.push((cell, temperature.to_bits()));
                ticket(&spec, temperature, query)
            },
        )
        .unwrap();
    let mut queries3 = Vec::new();
    let report3 = three
        .step_with_boussinesq_material_transport_endpoint_averaged(
            &power(),
            &[external3; 16],
            limits(),
            config(gravity3),
            &gate,
            |cell, temperature, reference| {
                assert_eq!(*reference, retained);
                let query = queries3.len();
                queries3.push((cell, temperature.to_bits()));
                ticket(&spec, temperature, query)
            },
        )
        .unwrap();
    let accepted2: Vec<_> = two
        .lattice
        .gpop
        .iter()
        .map(|p| two.temperature_si(p.iter().sum()).unwrap())
        .collect();
    let accepted3: Vec<_> = three
        .lattice
        .gpop
        .iter()
        .map(|p| three.temperature_si(p.iter().sum()).unwrap())
        .collect();
    check_report(
        &report2,
        &queries2,
        &accepted2,
        &two.lattice.grid.fext,
        external2,
        gravity2,
        retained,
    );
    check_report(
        &report3,
        &queries3,
        &accepted3,
        &three.lattice.grid.fext,
        external3,
        gravity3,
        retained,
    );
    assert_eq!(two.lattice.gbeta, 0.25);
    assert_eq!(three.lattice.gbeta, 0.25);
    // Replay only in this independent test oracle: the ordinary native step
    // receives the ORIGINAL initial rates and the independently checked held
    // average force. Production must produce exactly this single heat/flow pass.
    let heat = power().map(|p| {
        (((p / retained.density_kg_m3()) / retained.specific_heat_j_kg_k()) * 0.001) / 10.0
    });
    let native_limits = ThermalSourceLimits {
        max_cells: 16,
        temperature_range: [-1.0, 3.0],
        balance_tolerance: 1e-10,
    };
    let mut momentum2 = [1.0; 16];
    let mut momentum3 = [1.0; 16];
    let mut thermal2 = [1.0; 16];
    let mut thermal3 = [1.0; 16];
    for sample in &report2.transport.applied {
        momentum2[sample.cell] = sample.momentum_relaxation;
        thermal2[sample.cell] = sample.thermal_relaxation;
    }
    for sample in &report3.transport.applied {
        momentum3[sample.cell] = sample.momentum_relaxation;
        thermal3[sample.cell] = sample.thermal_relaxation;
    }
    control2
        .lattice
        .step_with_sources_and_transport(
            &heat,
            &two.lattice.grid.fext,
            RelaxationField {
                momentum: &momentum2,
                thermal: &thermal2,
            },
            native_limits,
        )
        .unwrap();
    control3
        .lattice
        .step_with_sources_and_transport(
            &heat,
            &three.lattice.grid.fext,
            RelaxationField {
                momentum: &momentum3,
                thermal: &thermal3,
            },
            native_limits,
        )
        .unwrap();
    for (actual, expected) in two
        .lattice
        .grid
        .f
        .iter()
        .flatten()
        .zip(control2.lattice.grid.f.iter().flatten())
        .chain(
            two.lattice
                .gpop
                .iter()
                .flatten()
                .zip(control2.lattice.gpop.iter().flatten()),
        )
    {
        close(*actual, *expected, 3e-14);
    }
    for (actual, expected) in three
        .lattice
        .grid
        .f
        .iter()
        .flatten()
        .zip(control3.lattice.grid.f.iter().flatten())
        .chain(
            three
                .lattice
                .gpop
                .iter()
                .flatten()
                .zip(control3.lattice.gpop.iter().flatten()),
        )
    {
        close(*actual, *expected, 3e-14);
    }
    assert_eq!(two.lattice.grid.tau, control2.lattice.grid.tau);
    assert_eq!(three.lattice.grid.tau, control3.lattice.grid.tau);
}

#[test]
fn force_query_and_exact_final_source_refusals_keep_the_whole_material_step_atomic() {
    let (spec, _, _, mut fluid) = setup();
    let before = (
        fluid.lattice.grid.f.clone(),
        fluid.lattice.gpop.clone(),
        fluid.lattice.grid.fext.clone(),
        fluid.lattice.grid.tau.clone(),
        fluid.lattice.last_wall_exchange(),
    );
    for (refuse_at, cancel) in [(12, false), (24, false), (24, true)] {
        let gate = CancelGate::new();
        let mut calls = 0;
        let refused = fluid.step_with_boussinesq_material_transport_endpoint_averaged(
            &power(),
            &[[0.0; 3]; 16],
            limits(),
            config([0.0, 0.0, -9.81]),
            &gate,
            |_, temperature, _| {
                let query = calls;
                calls += 1;
                if query == refuse_at {
                    if cancel {
                        gate.request();
                    } else {
                        return Err(MaterialError::State {
                            what: "original endpoint source refusal".into(),
                        });
                    }
                }
                ticket(&spec, temperature, query)
            },
        );
        assert_eq!(calls, refuse_at + 1);
        if cancel {
            assert!(matches!(refused, Err(MaterialTransportError::Cancelled)));
        } else {
            let error = refused.unwrap_err();
            assert!(matches!(
                error,
                MaterialTransportError::State { cell: 2, .. }
            ));
            assert!(
                error
                    .to_string()
                    .contains("original endpoint source refusal")
            );
        }
        assert_eq!(
            (
                fluid.lattice.grid.f.clone(),
                fluid.lattice.gpop.clone(),
                fluid.lattice.grid.fext.clone(),
                fluid.lattice.grid.tau.clone(),
                fluid.lattice.last_wall_exchange()
            ),
            before
        );
    }
    let (_, _, _, mut clean) = setup();
    for state in [&mut fluid, &mut clean] {
        state
            .step_with_boussinesq_material_transport_endpoint_averaged(
                &power(),
                &[[0.0; 3]; 16],
                limits(),
                config([0.0, 0.0, -9.81]),
                &CancelGate::new(),
                |_, temperature, _| ticket(&spec, temperature, 0),
            )
            .unwrap();
    }
    assert_eq!(fluid.lattice.grid.f, clean.lattice.grid.f);
    assert_eq!(fluid.lattice.gpop, clean.lattice.gpop);
    assert_eq!(fluid.lattice.grid.fext, clean.lattice.grid.fext);
}
