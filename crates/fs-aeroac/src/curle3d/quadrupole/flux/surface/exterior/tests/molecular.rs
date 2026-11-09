use super::*;
use crate::curle3d::quadrupole::flux::tests::molecular::{frozen_law, law};
use fs_flux::compressible::spatial::{mixture::molecular as multi, molecular as single};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, DissipationError, TransportCoefficients, TransportLaw,
};
use fs_qty::Pressure;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

struct CountedLaw {
    inner: ConstantTransport,
    evaluations: usize,
}
impl TransportLaw for CountedLaw {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.inner.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.evaluations += 1;
        self.inner.evaluate(temperature, pressure, fractions, binary, cancelled)
    }
}

fn molecular_flow(mu: f64, bulk: f64) -> single::MolecularEuler3d<ConstantTransport> {
    let native = frozen_flow(gas(6), grid());
    let transport = frozen_law(native.gas(), mu, bulk);
    single::MolecularEuler3d::new(native, transport, || false).unwrap()
}

fn record_molecular<L: TransportLaw>(
    flow: &mut single::MolecularEuler3d<L>,
    samples: usize,
) -> ExteriorSurfaceWelchStream {
    let source_config = settings(flow.flow().gas());
    let clock = SourceSampleClock::new(flow.flow().time().value(), DT, 1e-17).unwrap();
    ExteriorSurfaceWelchStream::from_molecular_frozen(
        flow, BOX, &EXTERIOR, source_config, clock, samples, config(), limits(0, samples),
        || false,
    ).unwrap()
}

#[test]
fn g1_joint_molecular_source_evaluates_native_stress_once_and_keeps_coherent_factors() {
    let native = frozen_flow(gas(6), grid());
    let transport = CountedLaw {
        inner: frozen_law(native.gas(), 0.12, 0.07),
        evaluations: 0,
    };
    let mut flow = single::MolecularEuler3d::new(native, transport, || false).unwrap();
    let before = flow.law().evaluations;
    drop(flow.stress_field(|| false).unwrap());
    let one_evaluation = flow.law().evaluations - before;
    assert!(one_evaluation > 64);
    let before = flow.law().evaluations;
    let mut stream = record_molecular(&mut flow, 4);
    assert_eq!(flow.law().evaluations - before, one_evaluation);
    assert!(stream.surface().includes_viscous_stress());
    assert!(stream.volume().includes_viscous_stress());
    let mut history = vec![frame(stream.surface(), stream.volume())];
    for _ in 1..4 {
        flow.step(Time::new(DT), single::MolecularStepControl::default(), || false).unwrap();
        let accepted = flow.flow().states().to_vec();
        let before = flow.law().evaluations;
        stream.sample_molecular_frozen(&mut flow, || false).unwrap();
        assert_eq!(flow.law().evaluations - before, one_evaluation);
        assert_eq!(flow.flow().states(), accepted);
        assert_eq!(stream.surface().sample_time_s(), flow.flow().time().value());
        assert_eq!(stream.volume().sample_time_s(), flow.flow().time().value());
        history.push(frame(stream.surface(), stream.volume()));
    }
    let direct = batch(&history, 0);
    assert_factors(stream.finish().unwrap().spectrum(), &direct);
}

#[test]
fn g4_molecular_joint_excludes_enclosed_sources_and_rolls_back_partial_and_late_samples() {
    let mut flow = molecular_flow(0.12, 0.07);
    let source_config = settings(flow.flow().gas());
    assert!(ExteriorSurfaceWelchStream::from_molecular_frozen(
        &mut flow, BOX, &[21], source_config,
        SourceSampleClock::new(0.0, DT, 1e-17).unwrap(), 4, config(), limits(0, 4), || false,
    ).is_err());
    let mut stream = record_molecular(&mut flow, 4);
    let mut reference = record_molecular(&mut flow, 4);
    let mut inviscid = recorder(flow.flow(), 4);
    for n in 1..3 {
        flow.step(Time::new(DT), single::MolecularStepControl::default(), || false).unwrap();
        if n == 1 {
            assert!(matches!(
                stream.sample_frozen(flow.flow(), || false),
                Err(ExteriorSurfaceRecordingError::Source(CompressibleVolumeError3::SourceMismatch)),
            ));
            assert!(matches!(
                inviscid.sample_molecular_frozen(&mut flow, || false),
                Err(ExteriorSurfaceRecordingError::Source(CompressibleVolumeError3::SourceMismatch)),
            ));
        }
        stream.sample_molecular_frozen(&mut flow, || false).unwrap();
        reference.sample_molecular_frozen(&mut flow, || false).unwrap();
    }
    let previous = frame(stream.surface(), stream.volume());
    let time = stream.surface().sample_time_s();
    let budget = limits(0, 4);
    let mut surface_probe = PreparedPermeableSurface3::from_molecular_frozen(
        &mut flow, BOX, source_config, budget.surface, || false,
    ).unwrap();
    let mut volume_probe = PreparedCompressibleVolume3::from_molecular_frozen(
        &mut flow, &EXTERIOR, source_config, budget.maximum_volume_workspace_bytes, || false,
    ).unwrap();
    flow.step(Time::new(DT), single::MolecularStepControl::default(), || false).unwrap();
    let states = flow.flow().states().to_vec();
    let report = flow.last_step().copied();
    let mut native_polls = 0;
    let field = flow.stress_field(|| {
        native_polls += 1;
        false
    }).unwrap();
    let mut surface_polls = 0;
    surface_probe.sample_stress(&field, &mut || {
        surface_polls += 1;
        false
    }).unwrap();
    let mut volume_polls = 0;
    volume_probe.sample_stress(&field, &mut || {
        volume_polls += 1;
        false
    }).unwrap();
    drop(field);
    // One joint cadence poll precedes the single native evaluation. Cancel at
    // the volume's final poll: surface has published, volume has only staged.
    let cutoff = 1 + native_polls + surface_polls + volume_polls;
    let mut calls = 0;
    assert!(matches!(
        stream.sample_molecular_frozen(&mut flow, || {
            calls += 1;
            calls == cutoff
        }),
        Err(ExteriorSurfaceRecordingError::Source(CompressibleVolumeError3::Cancelled)),
    ));
    assert_eq!(calls, cutoff);
    assert_eq!(stream.accepted_sample_count(), 3);
    assert_eq!(frame(stream.surface(), stream.volume()), previous);
    assert_eq!(stream.surface().sample_time_s(), time);
    assert_eq!(stream.volume().sample_time_s(), time);

    let mut total_polls = 0;
    reference.sample_molecular_frozen(&mut flow, || {
        total_polls += 1;
        false
    }).unwrap();
    assert!(total_polls > cutoff);
    let mut calls = 0;
    assert!(stream.sample_molecular_frozen(&mut flow, || {
        calls += 1;
        calls == total_polls
    }).is_err());
    assert_eq!(calls, total_polls);
    assert_eq!(stream.accepted_sample_count(), 3);
    assert_eq!(stream.completed_segment_count(), 0);
    assert_eq!(frame(stream.surface(), stream.volume()), previous);
    assert_eq!(stream.surface().sample_time_s(), time);
    assert_eq!(stream.volume().sample_time_s(), time);
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.last_step().copied(), report);
    stream.sample_molecular_frozen(&mut flow, || false).unwrap();
    assert_factors(stream.finish().unwrap().spectrum(), reference.finish().unwrap().spectrum());
}

#[test]
fn g0_zero_transport_joint_source_retains_original_euler_record_bit_for_bit() {
    let mut flow = molecular_flow(0.0, 0.0);
    let mut full = record_molecular(&mut flow, 4);
    let mut inviscid = recorder(flow.flow(), 4);
    for n in 0..4 {
        if n != 0 {
            flow.step(Time::new(DT), single::MolecularStepControl::default(), || false).unwrap();
            full.sample_molecular_frozen(&mut flow, || false).unwrap();
            inviscid.sample_frozen(flow.flow(), || false).unwrap();
        }
        assert_eq!(frame(full.surface(), full.volume()), frame(inviscid.surface(), inviscid.volume()));
    }
    assert_factors(full.finish().unwrap().spectrum(), inviscid.finish().unwrap().spectrum());
}

#[test]
fn g1_molecular_mixture_surface_volume_record_uses_one_original_species_state() {
    let source_config = settings(&gas(6));
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    ).unwrap();
    let species = gas.species_count();
    let mut ambient = vec![0.0; species];
    ambient[6] = source_config.acoustic.density_kg_m3();
    let transport = law(gas.source().clone(), 0.12, 0.07);
    let mut states = vec![0.0; 64 * (species + 4)];
    for (cell, row) in states.chunks_exact_mut(species + 4).enumerate() {
        let x = (cell % 4) as f64;
        let z = (cell / 16) as f64;
        let mut partial = vec![0.0; species];
        partial[0] = 0.1 + 0.01 * x;
        partial[6] = 1.6 - 0.01 * x;
        state_into(
            &gas, &partial, [0.2 + 0.1 * x, -0.3, 0.5 - 0.05 * z],
            source_config.reference_pressure_pa + 100.0 * (x - z), row,
        ).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid()).unwrap();
    let native = MixtureEuler3d::new(gas, grid(), states, PERIODIC, bytes, || false).unwrap();
    let source = multi::SpatialTransport::new(native.gas().clone(), transport, || false).unwrap();
    let bytes = multi::MolecularTransport3d::required_workspace_bytes(&native, &source).unwrap();
    let operator = multi::MolecularTransport3d::new(&native, source, bytes, || false).unwrap();
    let bytes = multi::MolecularEuler3d::required_workspace_bytes(&native, &operator).unwrap();
    let mut flow = multi::MolecularEuler3d::new(native, operator, bytes, || false).unwrap();
    let budget = limits(species, 4);
    let mut stream = ExteriorSurfaceWelchStream::from_molecular_mixture(
        &mut flow, &ambient, BOX, &EXTERIOR, source_config,
        SourceSampleClock::new(0.0, DT, 1e-17).unwrap(), 4, config(), budget,
        multi::DissipationControl::default(), || false,
    ).unwrap();
    let mut history = Vec::new();
    for n in 0..4 {
        if n != 0 {
            flow.step(Time::new(DT), multi::MolecularStepControl::default(), || false).unwrap();
            stream.sample_molecular_mixture(
                &mut flow, multi::DissipationControl::default(), || false,
            ).unwrap();
        }
        let before = flow.flow().states().to_vec();
        let surface = PreparedPermeableSurface3::from_molecular_mixture(
            &mut flow, &ambient, BOX, source_config, budget.surface,
            multi::DissipationControl::default(), || false,
        ).unwrap();
        let volume = PreparedCompressibleVolume3::from_molecular_mixture(
            &mut flow, &EXTERIOR, source_config, budget.maximum_volume_workspace_bytes,
            multi::DissipationControl::default(), || false,
        ).unwrap();
        assert_eq!(stream.surface().reference_partial_densities_kg_m3(), ambient);
        assert_eq!(stream.surface().mixture_source(), Some(flow.flow().gas()));
        assert_eq!(stream.volume().mixture_source(), Some(flow.flow().gas()));
        assert_eq!(frame(stream.surface(), stream.volume()), frame(&surface, &volume));
        assert_eq!(flow.flow().states(), before);
        history.push(frame(&surface, &volume));
    }
    let expected = batch(&history, species);
    assert_factors(stream.finish().unwrap().spectrum(), &expected);
}
