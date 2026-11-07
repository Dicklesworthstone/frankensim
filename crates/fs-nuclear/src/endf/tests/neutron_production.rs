//! Exclusive neutron-production data on the original finite-inventory owner.
//! All coefficients and decay laws are manufactured; no isotope-data validation.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::irradiation::{GroupFlux, IrradiationModel, TabulatedIrradiationData};
use crate::{EnergyChannels, FissionModel};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn curve(kind: NeutronReaction, temperature: f64, sigma_barns: f64) -> PointwiseCrossSection {
    let selected = selection(temperature, kind);
    parse(
        &fixture(
            selected,
            2,
            &[(0.0, 0.0), (2e6, sigma_barns), (2e7, sigma_barns)],
            &[(3, 1)],
        ),
        selected,
    )
    .unwrap()
}

fn table(kind: NeutronReaction) -> Arc<CrossSectionTable> {
    Arc::new(
        collapse_temperature_table(
            "manufactured-exclusive-neutron-production".into(),
            &[curve(kind, 300.0, 2.0), curve(kind, 600.0, 6.0)],
            groups(&[2e6, 2e7]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn data(table: Arc<CrossSectionTable>, emitted: u16) -> BurnupData {
    let parent = table.nuclide();
    let daughter = NuclideId::new(
        parent.atomic_number(),
        parent.mass_number() + 1 - emitted,
        0,
    )
    .unwrap();
    BurnupData {
        data_id: "manufactured-residual-inventory".into(),
        groups: Arc::clone(table.groups()),
        temperature_min_k: 300.0,
        temperature_max_k: 600.0,
        nuclides: [parent, daughter]
            .map(|id| BurnupNuclide {
                id,
                law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
            })
            .to_vec(),
        reactions: vec![BurnupReaction {
            data_id: "manufactured-exclusive-neutron-removal".into(),
            parent: 0,
            rate: BurnupRate::TabulatedNeutronProduction {
                cross_sections: table,
            },
            products: vec![BurnupProduct {
                nuclide: 1,
                yield_per_event: 1.0,
            }],
            // MF3 Q remains metadata. The transport energy owner must explicitly
            // declare incident/outgoing kinetic energies and signed prompt Q.
            energy_per_event: EnergyChannels::default(),
        }],
    }
}

fn budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-18,
    }
}

#[test]
fn g1_mt16_mt17_mt37_keep_exclusive_identity_and_finite_target_depletion() {
    for (kind, mt, emitted) in [
        (NeutronReaction::TwoNeutronProduction, 16, 2),
        (NeutronReaction::ThreeNeutronProduction, 17, 3),
        (NeutronReaction::FourNeutronProduction, 37, 4),
    ] {
        assert_eq!(reaction_mt(kind), mt);
        assert_eq!(kind.emitted_neutrons_per_event(), Some(emitted));
        assert!(!kind.is_nonfission_absorption());
        let source = curve(kind, 300.0, 2.0);
        assert_eq!(source.reaction(), kind);
        close(source.reaction_q_j(), -1e6 * EV_J);
        let coefficients = table(kind);
        close(
            coefficients.cross_section_at_group(450.0, 0).unwrap(),
            4e-28,
        );
        let model = Arc::new(
            BurnupModel::new(data(coefficients, emitted as u16), BurnupLimits::default()).unwrap(),
        );
        let flux = GroupFlux::new(Arc::clone(model.groups()), vec![1e28]).unwrap();
        let initial = BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0]).unwrap();
        let (next, report) = initial
            .trial(0.25, &flux, 450.0, budget(), || Ok(()))
            .unwrap();
        let events = 100.0 * -(-1.0_f64).exp_m1();
        close(next.populations()[0], 100.0 * (-1.0_f64).exp());
        close(next.populations()[1], events);
        close(report.neutron_producing_reactions, events);
        close(report.reaction_counts[0], events);
        assert_eq!([report.fissions, report.captures, report.decays], [0.0; 3]);
        assert_eq!(report.generated_energy_j, 0.0);
        assert_eq!(report.neutron_production_energy_j, 0.0);
        let masses = model
            .nuclides()
            .iter()
            .map(|n| f64::from(n.id.mass_number()));
        let retained: f64 = masses.zip(next.populations()).map(|(a, n)| a * n).sum();
        close(
            retained + (f64::from(emitted) - 1.0) * events,
            235.0 * 100.0,
        );
        assert_eq!(initial.populations(), &[100.0, 0.0]);
    }
}

#[test]
fn g3_neutron_production_uses_original_continuous_supply_and_restart_accounts() {
    let mut authored = data(table(NeutronReaction::TwoNeutronProduction), 2);
    authored.reactions[0].rate = BurnupRate::NeutronProduction {
        cross_sections_m2: vec![4e-28],
    };
    // An independent explicit event-energy source verifies the disjoint report;
    // it is not copied from the negative-Q ENDF metadata or a stopping model.
    authored.reactions[0].energy_per_event = EnergyChannels {
        deposited_j: 2e-13,
        radiation_j: 3e-13,
        escaped_j: 1e-13,
    };
    let model = Arc::new(BurnupModel::new(authored, BurnupLimits::default()).unwrap());
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![1e28]).unwrap();
    let initial = BurnupInventory::new(model, vec![100.0, 0.0]).unwrap();
    let (full, report) = initial
        .trial_with_production(0.25, &flux, 450.0, &[20.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let (half, first) = initial
        .trial_with_production(0.125, &flux, 450.0, &[20.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let (resumed, second) = half
        .trial_with_production(0.125, &flux, 450.0, &[20.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let remaining = 100.0 * (-1.0_f64).exp() + 20.0 / 4.0 * -(-1.0_f64).exp_m1();
    let events = 100.0 + 20.0 * 0.25 - remaining;
    close(full.populations()[0], remaining);
    close(full.populations()[1], events);
    for (actual, expected) in resumed.populations().iter().zip(full.populations()) {
        close(*actual, *expected);
    }
    close(report.burnup.neutron_producing_reactions, events);
    close(
        first.burnup.neutron_producing_reactions + second.burnup.neutron_producing_reactions,
        events,
    );
    close(report.burnup.neutron_production_energy_j, events * 6e-13);
    close(report.burnup.released.deposited_j, events * 2e-13);
    assert_eq!(report.burnup.capture_energy_j, 0.0);
    assert_eq!(report.burnup.fission_energy_j, 0.0);
}

#[test]
fn g0_neutron_production_tables_cannot_be_silently_used_as_absorption_or_fission() {
    let production = table(NeutronReaction::TwoNeutronProduction);
    for rate in [
        BurnupRate::TabulatedCapture {
            cross_sections: Arc::clone(&production),
        },
        BurnupRate::TabulatedFission {
            cross_sections: Arc::clone(&production),
        },
    ] {
        let mut wrong = data(Arc::clone(&production), 2);
        wrong.reactions[0].rate = rate;
        assert!(BurnupModel::new(wrong, BurnupLimits::default()).is_err());
    }
    for kind in [
        NeutronReaction::Capture,
        NeutronReaction::AlphaProduction,
        NeutronReaction::Fission,
    ] {
        assert_eq!(kind.emitted_neutrons_per_event(), None);
        assert!(BurnupModel::new(data(table(kind), 2), BurnupLimits::default()).is_err());
    }
    let fission = table(NeutronReaction::Fission);
    assert!(
        IrradiationModel::new_tabulated(TabulatedIrradiationData {
            data_id: "no-hidden-production-in-absorption".into(),
            fission,
            capture: production,
            temperature_min_k: 300.0,
            temperature_max_k: 600.0,
            fission_model: Arc::new(
                FissionModel::new(
                    "inactive-fission".into(),
                    EnergyChannels {
                        deposited_j: 1.0e-11,
                        ..EnergyChannels::default()
                    },
                    vec![],
                )
                .unwrap()
            ),
            capture_energy_per_event: EnergyChannels::default(),
        })
        .is_err()
    );
    let selected = selection(300.0, NeutronReaction::TwoNeutronProduction);
    let original = fixture(selected, 2, &[(0.0, 0.0), (2e7, 2.0)], &[(2, 2)]);
    let section = (i32::from(selected.material_number), 3, 16);
    let breakup = original.replacen(
        &control(2e8, -1e6, [0, 0, 1, 2], section),
        &control(2e8, -1e6, [0, 22, 1, 2], section),
        1,
    );
    assert_ne!(original, breakup);
    assert_eq!(
        parse(&breakup, selected),
        Err(EndfError::Unsupported("MF3 breakup flag LR must be zero"))
    );
}
