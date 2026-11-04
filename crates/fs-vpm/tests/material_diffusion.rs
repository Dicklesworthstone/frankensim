//! G0/G2/G3 material substitutions in the actual conservative PSE owner.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_vpm::diffusion::{
    MaterialPseDiffusion, MaterialPseError, PseDiffusion, PseKernel, diffuse_with_cx,
    simulate_viscous_with_cx,
};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, VpmRun, total_circulation};
use std::cell::Cell;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let clock = VirtualClock::new();
    let stream = StreamKey {
        seed: 0x4D41_545F_5053_45,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            stream,
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn diffuse(particles: &[VortexParticle], areas: &[f64], dt: f64, law: PseDiffusion) -> VpmRun {
    with_cx(|cx| {
        diffuse_with_cx(
            cx,
            particles,
            areas,
            dt,
            1,
            law,
            VpmBudget::new(particles.len(), 1, u64::MAX, usize::MAX),
        )
    })
    .unwrap()
}

#[test]
fn g2_air_state_changes_real_diffusion_with_the_exact_pair_relaxation() {
    let air = GasSpec::dry_air_ussa1976();
    let states = [
        GasState::try_new(&air, 293.15, 101_325.0).unwrap(),
        GasState::try_new(&air, 373.15, 101_325.0).unwrap(),
        GasState::try_new(&air, 293.15, 202_650.0).unwrap(),
    ];
    let particles = [
        VortexParticle::new([0.0, 0.0], 2e-6),
        VortexParticle::new([0.001, 0.0], 0.0),
    ];
    let areas = [1e-6; 2];
    let width: f64 = 0.003;
    let dt = 0.2;
    let mut received = [0.0; 3];
    for (i, state) in states.iter().enumerate() {
        let binding = MaterialPseDiffusion::new(state, width).unwrap();
        assert!(std::ptr::eq(binding.source(), state));
        assert_eq!(binding.law().kernel(), PseKernel::Gaussian);
        assert_eq!(binding.law().support_radius_m(), None);
        let transport = binding.transport();
        let nu = transport.dynamic_viscosity_pa_s() / transport.density_kg_m3();
        assert_eq!(binding.law().kinematic_viscosity_m2_s(), nu);
        let result = diffuse(&particles, &areas, dt, binding.law());
        // The two-site ODE has decay rate k*(A0+A1), independently of the
        // native symmetric pair-sweep implementation and material conversion.
        let k = 4.0 * nu / (std::f64::consts::PI * width.powi(4))
            * (-0.001_f64.powi(2) / width.powi(2)).exp();
        let expected = 1e-6 * -(-k * (areas[0] + areas[1]) * dt).exp_m1();
        received[i] = result.particles[1].circulation;
        assert!((received[i] / expected - 1.0).abs() < 2e-13);
        assert!((total_circulation(&result.particles) - 2e-6).abs() < 1e-20);
    }
    // Heating raises mu and lowers rho; doubling pressure leaves mu fixed
    // while doubling rho. Both changes reach circulation through actual nu.
    assert!(received[1] > received[0]);
    assert!(received[2] < received[0]);
    let base = states[0].transport_state().unwrap();
    let pressurized = states[2].transport_state().unwrap();
    assert_eq!(
        base.dynamic_viscosity_pa_s(),
        pressurized.dynamic_viscosity_pa_s()
    );
    assert_eq!(
        pressurized.kinematic_viscosity_m2_s(),
        0.5 * base.kinematic_viscosity_m2_s()
    );
}

#[derive(Debug)]
struct SampledFluid {
    state: Cell<FluidTransportState>,
    queries: Cell<usize>,
}
impl FluidTransport for SampledFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state.get())
    }
}

#[test]
fn g0_frozen_snapshot_and_fresh_material_preparation_share_native_consumers() {
    let fluid = SampledFluid {
        state: Cell::new(FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.1).unwrap()),
        queries: Cell::new(0),
    };
    let source: &dyn FluidTransport = &fluid;
    let first = MaterialPseDiffusion::compact_wendland(source, 1.0).unwrap();
    assert!(std::ptr::eq(first.source(), source));
    assert_eq!(fluid.queries.get(), 1);
    let before = *first.transport();
    fluid
        .state
        .set(FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.2).unwrap());
    let fresh = MaterialPseDiffusion::compact_wendland(source, 1.0).unwrap();
    assert_eq!(fluid.queries.get(), 2);
    assert_eq!(*first.transport(), before);
    assert_eq!(
        fresh.law().kinematic_viscosity_m2_s(),
        2.0 * first.law().kinematic_viscosity_m2_s()
    );
    assert_eq!(fresh.law().kernel(), PseKernel::CompactWendland);
    assert_eq!(fresh.law().support_radius_m(), Some(1.0));

    let particles = [
        VortexParticle::new([0.0, 0.0], 0.3),
        VortexParticle::new([0.4, 0.0], -0.1),
    ];
    let areas = [0.2, 0.4];
    let old = diffuse(&particles, &areas, 0.1, first.law());
    let changed = diffuse(&particles, &areas, 0.1, fresh.law());
    assert!(changed.particles[0].circulation < old.particles[0].circulation);

    // This is the original RK4/D(dt/2) composition, not a new material solver.
    // Source conversion cannot change its work/byte charges or numerical path.
    let run = |law| {
        with_cx(|cx| {
            simulate_viscous_with_cx(
                cx,
                &particles,
                &areas,
                0.01,
                4,
                0.05,
                law,
                VpmBudget::new(2, 4, u64::MAX, usize::MAX),
            )
        })
        .unwrap()
    };
    let raw =
        PseDiffusion::compact_wendland(fresh.transport().kinematic_viscosity_m2_s(), 1.0).unwrap();
    let result = run(fresh.law());
    assert_eq!(result, run(raw));
    assert!((total_circulation(&result.particles) - total_circulation(&particles)).abs() < 1e-14);
    assert_eq!(
        fluid.queries.get(),
        2,
        "native stepping must not query the source again"
    );
}

#[test]
fn g0_compact_material_binding_preserves_the_authored_support_cutoff() {
    let fluid = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.1).unwrap();
    let binding = MaterialPseDiffusion::compact_wendland(&fluid, 0.5).unwrap();
    let particles = [
        VortexParticle::new([0.0, 0.0], 1.0),
        VortexParticle::new([0.5, 0.0], 0.0),
    ];
    let result = diffuse(&particles, &[1.0; 2], 10.0, binding.law());
    assert_eq!(result.particles, particles);
    let gaussian = MaterialPseDiffusion::new(&fluid, 0.5).unwrap();
    assert!(diffuse(&particles, &[1.0; 2], 10.0, gaussian.law()).particles[1].circulation > 0.0);
}

#[test]
fn g0_original_material_and_native_input_refusals_remain_typed() {
    #[derive(Debug)]
    struct Unavailable(MaterialError);
    impl FluidTransport for Unavailable {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            Err(self.0.clone())
        }
    }
    let source = Unavailable(MaterialError::State {
        what: "outside the source temperature range".into(),
    });
    let error = MaterialPseDiffusion::new(&source, 1.0).unwrap_err();
    assert_eq!(error, MaterialPseError::Material(source.0.clone()));
    assert_eq!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<MaterialError>(),
        Some(&source.0)
    );
    assert_eq!(
        MaterialPseDiffusion::compact_wendland(&source, 1.0).unwrap_err(),
        error
    );

    let state = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.1).unwrap();
    for width in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let expected = PseDiffusion::new(state.kinematic_viscosity_m2_s(), width).unwrap_err();
        let error = MaterialPseDiffusion::new(&state, width).unwrap_err();
        assert_eq!(
            std::error::Error::source(&error)
                .unwrap()
                .downcast_ref::<VpmError>(),
            Some(&expected)
        );
        assert_eq!(error, MaterialPseError::Diffusion(expected));
    }
}
