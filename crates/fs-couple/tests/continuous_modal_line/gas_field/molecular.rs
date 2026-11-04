//! All local gas reservoirs and original load states on the exact shared clock.
use super::*;
use fs_duct::modal::time_domain::ModalMolecularBudget;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn molecular_budget() -> ModalMolecularBudget {
    ModalMolecularBudget {
        maximum_states: 124,
        maximum_payload_bytes: 32 << 20,
        maximum_preparation_work: 1_000_000_000,
    }
}

fn law(active: bool) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["A", "B"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if active {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 1000.0 + 300.0 * i as f64,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 101325.0 * (2e-5 + i as f64 * 8e-5),
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 1000.0,
                        },
                    }]
                } else {
                    vec![]
                },
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn empty_spatial_reservoir_law_preserves_every_existing_runtime_sample_bit() {
    let profile = RadiusProfile::exponential(0.025, 3.0, 0.04).unwrap();
    let field = field(profile.length(), false);
    let mut old = ContinuousModalLine::new_field(
        &profile,
        &field,
        ModeCount(2),
        2,
        ModalEnd::Ports,
        1.0 / 96000.0,
        budget(),
        || false,
    )
    .unwrap();
    let mut new = ContinuousModalLine::new_field_molecular(
        &profile,
        &field,
        law(false),
        ModeCount(2),
        2,
        ModalEnd::Ports,
        1.0 / 96000.0,
        molecular_budget(),
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(old.state().len(), new.state().len());
    for step in 0..24 {
        let input = [70.0 * det::sin(step as f64 * 0.3), 20.0, -30.0, 4.0];
        let a = old.push(&input).unwrap();
        let b = new.push(&input).unwrap();
        assert_eq!(a.physical_heat.to_bits(), b.physical_heat.to_bits());
        assert_eq!(
            snapshot(&old)
                .into_iter()
                .map(f64::to_bits)
                .collect::<Vec<_>>(),
            snapshot(&new)
                .into_iter()
                .map(f64::to_bits)
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn local_molecular_and_mouth_memory_preserve_heat_work_and_clock_subdivision() {
    let profile = RadiusProfile::exponential(0.025, 2.0, 0.04).unwrap();
    let field = field(profile.length(), false);
    let outlet = *field.state_at(profile.length()).unwrap().state();
    let duration = 1.0 / 12000.0;
    let mut results = Vec::new();
    for subdivision in [1, 2, 4] {
        let load = coupled_load(profile.outlet_radius(), outlet, &[1.0, 0.6]);
        let mut line = ContinuousModalLine::new_field_molecular_loaded(
            &profile,
            &field,
            law(true),
            ModeCount(2),
            2,
            load,
            duration / subdivision as f64,
            molecular_budget(),
            budget(),
            || false,
        )
        .unwrap();
        assert_eq!(line.source().state_dimension(), 36);
        assert_eq!(line.state().len(), 38);
        assert_eq!(line.load_state().len(), 2);
        assert_eq!(line.source().gas_field(), Some(&field));
        let (mut incident, mut outgoing, mut heat) = (0.0, 0.0, 0.0);
        for _ in 0..subdivision {
            let record = line.push(&[75.0, 15.0]).unwrap();
            incident += record.incident_energy;
            outgoing += record.outgoing_energy;
            heat += record.physical_heat;
            assert_eq!(record.external_work, 0.0);
        }
        let local_heat = line
            .source()
            .molecular_heat_per_process_w(line.duct_state())
            .unwrap();
        assert!(local_heat.iter().sum::<f64>() > 0.0);
        assert!(heat > 0.0);
        close(
            line.stored_energy_j() + heat + outgoing,
            incident,
            4e-8,
            1e-14,
        );
        results.push((line.state().to_vec(), heat, outgoing));
    }
    for result in &results[1..] {
        for (a, b) in result.0.iter().zip(&results[0].0) {
            close(*a, *b, 3e-7, 1e-12);
        }
        close(result.1, results[0].1, 3e-7, 1e-14);
        close(result.2, results[0].2, 3e-7, 1e-14);
    }
}
