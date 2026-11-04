//! Manufactured normal-dependent heat transfer on the original moving PLIC drop.
//! This verifies geometry and coupling, not a correlation for curved surfaces.

use super::*;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::{
    FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceGasHeatTransfer,
    FreeSurfaceHeatTransferCoefficient,
};
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

fn coefficient(normal: [f64; 3]) -> f64 {
    5.0 * (1.0 + normal[0] * normal[0] + 0.5 * normal[1] * normal[1])
}

struct NormalLaw {
    initial: BTreeMap<[usize; 3], [u64; 3]>,
    uniform: bool,
    refuse_evolved: bool,
    calls: Cell<usize>,
    evolved_calls: Cell<usize>,
}

impl FreeSurfaceGasConvection for NormalLaw {
    fn coefficient(
        &self,
        query: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        self.calls.set(self.calls.get() + 1);
        assert_eq!(query.surface.face.liquid, query.surface.face.gas);
        assert!(query.surface.face.area_m2 > 0.0);
        near(
            query.outward_normal.into_iter().map(|x| x * x).sum::<f64>(),
            1.0,
            2e-14,
        );
        assert_eq!(
            query.film.state().temperature,
            0.5 * query.surface_temperature_k + 0.5 * query.bulk.state().temperature
        );
        let bits = query.outward_normal.map(f64::to_bits);
        if self.initial.get(&query.surface.face.liquid) != Some(&bits) {
            self.evolved_calls.set(self.evolved_calls.get() + 1);
            if self.refuse_evolved {
                return Err(FreeSurfaceEnthalpyError::Convection(
                    "manufactured refusal at the evolved native interface".into(),
                ));
            }
        }
        let h = coefficient(query.outward_normal);
        Ok(Some(Evidence {
            value: FreeSurfaceHeatTransferCoefficient::new(h),
            qoi: h,
            numerical: NumericalCertificate::estimate(h, h),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::of_bytes(
                b"manufactured PLIC normal-dependent heat transfer",
            ),
            adjoint_ref: None,
        }))
    }

    fn admits_uniform_boundary_evolution(&self) -> bool {
        self.uniform
    }
}

#[test]
fn reconstructed_convection_follows_plic_owners_and_rolls_back_a_late_refusal() {
    let gas = gas_source::model(0.0, 250000.0, 1.0);
    let point = gas_source::point(450.0, 100000.0);
    let controls = FreeSurfaceConductionControls {
        maximum_dt_s: DT,
        ..heat_controls()
    };
    with_moving_drop(0.0, [0.01, 0.0, 0.0], [0.0; 3], |thermal, gate| {
        let mut changed_owners = false;
        for interval in 0..2 {
            let patches = thermal.local_interface_patches(gate).unwrap();
            assert!(
                patches
                    .iter()
                    .any(|p| p.outward_normal.iter().all(|n| n.abs() < 0.99))
            );
            let initial: BTreeMap<_, _> = patches
                .iter()
                .map(|p| (p.liquid, p.outward_normal.map(f64::to_bits)))
                .collect();
            let owners: BTreeSet<_> = initial.keys().copied().collect();
            let boundaries: Vec<_> = patches
                .iter()
                .map(|patch| FreeSurfaceInterfaceContact {
                    surface: patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0),
                    model: &gas,
                    context: &point,
                })
                .collect();
            let law = NormalLaw {
                initial,
                uniform: true,
                refuse_evolved: false,
                calls: Cell::new(0),
                evolved_calls: Cell::new(0),
            };
            let energy = thermal.total_enthalpy_j().unwrap();
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            if interval == 0 {
                let populations = thermal.flow().lattice().grid().state_bits();
                let mut refusing = NormalLaw {
                    uniform: false,
                    ..law
                };
                // New-owner evolution needs an explicit constitutive capability.
                assert!(
                    thermal
                        .step_with_surface_convection(
                            &[],
                            &[],
                            controls,
                            &mut [],
                            &[],
                            &[],
                            &boundaries,
                            uniform_solve(),
                            &refusing,
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(refusing.calls.get(), 0);
                refusing.uniform = true;
                refusing.refuse_evolved = true;
                assert!(matches!(
                    thermal.step_with_surface_convection(
                        &[],
                        &[],
                        controls,
                        &mut [],
                        &[],
                        &[],
                        &boundaries,
                        uniform_solve(),
                        &refusing,
                        false,
                        gate,
                    ),
                    Err(FreeSurfaceEnthalpyError::Convection(_))
                ));
                assert!(
                    refusing.evolved_calls.get() > 0,
                    "refusal must occur after geometry evolves"
                );
                // Do not print the full D3Q19 population array on a failure.
                assert!(populations == thermal.flow().lattice().grid().state_bits());
                assert_eq!(
                    energy.to_bits(),
                    thermal.total_enthalpy_j().unwrap().to_bits()
                );
                assert_eq!(
                    mass.to_bits(),
                    thermal.flow().ledger_mass_kg().unwrap().to_bits()
                );
                refusing.refuse_evolved = false;
                refusing.evolved_calls.set(0);
                let report = thermal
                    .step_with_surface_convection(
                        &[],
                        &[],
                        controls,
                        &mut [],
                        &[],
                        &[],
                        &boundaries,
                        uniform_solve(),
                        &refusing,
                        false,
                        gate,
                    )
                    .unwrap()
                    .conduction;
                check_accepted(
                    thermal,
                    gate,
                    &refusing,
                    &owners,
                    energy,
                    mass,
                    &report,
                    &mut changed_owners,
                );
            } else {
                let report = thermal
                    .step_with_surface_convection(
                        &[],
                        &[],
                        controls,
                        &mut [],
                        &[],
                        &[],
                        &boundaries,
                        uniform_solve(),
                        &law,
                        false,
                        gate,
                    )
                    .unwrap()
                    .conduction;
                check_accepted(
                    thermal,
                    gate,
                    &law,
                    &owners,
                    energy,
                    mass,
                    &report,
                    &mut changed_owners,
                );
            }
        }
        assert!(
            changed_owners,
            "native motion must expose new PLIC heat owners"
        );
    });
}

fn check_accepted(
    thermal: &fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gate: &CancelGate,
    law: &NormalLaw,
    owners: &BTreeSet<[usize; 3]>,
    energy: f64,
    mass: f64,
    report: &fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceConductionReport,
    changed_owners: &mut bool,
) {
    let exchange = report.surface.as_ref().unwrap();
    let live = thermal.local_interface_patches(gate).unwrap();
    let current: BTreeMap<_, _> = live.iter().map(|p| (p.liquid, p)).collect();
    *changed_owners |= *owners != current.keys().copied().collect();
    assert_eq!(exchange.native_endpoint.len(), live.len());
    assert!(law.evolved_calls.get() > 0);
    for sample in exchange
        .last_applied
        .iter()
        .chain(&exchange.thermal_endpoint)
    {
        assert_eq!(
            sample.convection_outward_normal.unwrap().map(f64::to_bits),
            law.initial[&sample.surface.face.liquid]
        );
    }
    for sample in &exchange.native_endpoint {
        let patch = current[&sample.surface.face.liquid];
        assert_eq!(
            sample.convection_outward_normal.unwrap().map(f64::to_bits),
            patch.outward_normal.map(f64::to_bits)
        );
        assert!(
            (sample.surface.face.area_m2 - patch.area_m2).abs() <= 64.0 * f64::EPSILON * DX * DX
        );
    }
    for sample in exchange
        .last_applied
        .iter()
        .chain(&exchange.thermal_endpoint)
        .chain(&exchange.native_endpoint)
    {
        let h = coefficient(sample.convection_outward_normal.unwrap());
        let FreeSurfaceGasHeatTransfer::Convection(Some(evidence)) = &sample.gas_heat_transfer
        else {
            panic!("missing convection evidence")
        };
        near(evidence.value.value(), h, 2e-14);
        let gl = 0.5 / (0.25 * DX);
        let ts = (gl * sample.liquid_temperature_k + h * 450.0) / (gl + h);
        near(sample.surface_temperature_k, ts, 2e-11);
        near(
            sample.gas_conductance_w_k,
            h * sample.surface.face.area_m2,
            2e-14,
        );
    }
    assert!(exchange.heat_from_gas_j > 0.0);
    assert_eq!(exchange.heat_from_radiation_j, 0.0);
    near(report.flow.external_heat_j, exchange.heat_from_gas_j, 2e-9);
    assert!((thermal.total_enthalpy_j().unwrap() - energy - exchange.heat_from_gas_j).abs() < 1e-7);
    near(thermal.flow().ledger_mass_kg().unwrap(), mass, 3e-12);
}
