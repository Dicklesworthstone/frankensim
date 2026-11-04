//! Manufactured curved-surface h checks geometry plumbing, not a plate model.
use super::*;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::{
    FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceGasHeatTransfer,
    FreeSurfaceHeatTransferCoefficient,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use std::cell::RefCell;

struct Seen {
    liquid: [usize; 3],
    normal: [f64; 3],
    area_m2: f64,
    temperatures: [f64; 3],
}
#[derive(Default)]
struct Manufactured(RefCell<Vec<Seen>>);
fn coefficient(normal: [f64; 3]) -> f64 {
    // Strictly positive for every unit normal; deliberately not an empirical law.
    200.0 + 40.0 * normal[0] + 30.0 * normal[1] + 20.0 * normal[2]
}
impl FreeSurfaceGasConvection for Manufactured {
    fn coefficient(
        &self,
        q: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        self.0.borrow_mut().push(Seen {
            liquid: q.surface.face.liquid,
            normal: q.outward_normal,
            area_m2: q.surface.face.area_m2,
            temperatures: [
                q.surface_temperature_k,
                q.bulk.state().temperature,
                q.film.state().temperature,
            ],
        });
        let h = coefficient(q.outward_normal);
        Ok(Some(Evidence {
            value: FreeSurfaceHeatTransferCoefficient::new(h),
            qoi: h,
            numerical: NumericalCertificate::estimate(h, h),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::of_bytes(
                b"manufactured curved normal-dependent h; not calibrated",
            ),
            adjoint_ref: None,
        }))
    }
}

#[test]
fn g1_height_graph_convection_uses_held_and_live_normals_with_one_finite_gas_heat() {
    let card = thermal_card(-1e-5, 500.0);
    with_graph(&card, 1, 1, |thermal, gate| {
        let held = thermal.height_graph_interface_patches(gate).unwrap();
        assert_eq!(held.len(), 64);
        assert!(
            held.iter()
                .all(|p| p.outward_normal[0] != 0.0 && p.outward_normal[2] != 0.0)
        );
        let surfaces: Vec<_> = held
            .iter()
            .map(|p| p.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0))
            .collect();
        let mut gas = IsochoricGasParcel::new(
            gas_source::model(0.0, 250000.0, 1000.0),
            gas_source::point(450.0, 100000.0),
            1e-4,
        )
        .unwrap();
        let liquid_before = thermal.total_enthalpy_j().unwrap();
        let gas_before = gas.state().specific_internal_energy_j_kg();
        let mass_before = thermal.flow().ledger_mass_kg().unwrap();
        let law = Manufactured::default();
        let report = thermal
            .step_with_finite_gas_surface_convection(
                &[],
                &[],
                surface_controls(),
                &mut [],
                &[],
                &[],
                &mut gas,
                &surfaces,
                GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 690.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                FreeSurfaceInterfaceSolve::default(),
                &law,
                None,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        let live = thermal.height_graph_interface_patches(gate).unwrap();
        assert!(live.iter().any(|p| {
            held.iter()
                .find(|h| h.liquid == p.liquid)
                .is_none_or(|h| h.outward_normal != p.outward_normal)
        }));
        let calls = law.0.borrow();
        assert!(
            calls
                .iter()
                .all(|c| c.temperatures[2] == 0.5 * c.temperatures[0] + 0.5 * c.temperatures[1])
        );
        for (samples, geometry) in [
            (&exchange.last_applied, &held),
            (&exchange.thermal_endpoint, &held),
            (&exchange.native_endpoint, &live),
        ] {
            assert_eq!(samples.len(), geometry.len());
            for sample in samples {
                let patch = geometry
                    .iter()
                    .find(|p| p.liquid == sample.surface.face.liquid)
                    .unwrap();
                assert_eq!(sample.surface.face.gas, patch.gas);
                assert_eq!(sample.convection_outward_normal, Some(patch.outward_normal));
                // Public reconstruction uses (area*dx)*dx; the native solve
                // uses area*(dx*dx). Only their conversion roundoff may differ.
                assert!(
                    (sample.surface.face.area_m2 - patch.area_m2).abs()
                        <= 64.0 * f64::EPSILON * DX * DX
                );
                assert!(calls.iter().any(|c| c.liquid == patch.liquid
                    && c.normal == patch.outward_normal
                    && c.area_m2 == sample.surface.face.area_m2
                    && (c.temperatures[0] - sample.surface_temperature_k).abs() < 1e-10));
                let FreeSurfaceGasHeatTransfer::Convection(Some(evidence)) =
                    &sample.gas_heat_transfer
                else {
                    panic!("lost manufactured coefficient evidence");
                };
                let h = coefficient(patch.outward_normal);
                near(evidence.value.value(), h, 1e-12);
                assert!(evidence.model.cards.is_empty());
                let gl = sample.liquid.thermal_conductivity_w_m_k() * patch.area_m2 / (0.25 * DX);
                let gg = h * patch.area_m2;
                let ts =
                    (gl * sample.liquid_temperature_k + gg * sample.gas_temperature_k) / (gl + gg);
                near(sample.gas_conductance_w_k, gg, 2e-12);
                near(sample.surface_temperature_k, ts, 2e-10);
                near(
                    sample.heat_from_gas_w,
                    gg * (sample.gas_temperature_k - ts),
                    2e-10,
                );
                near(
                    sample.liquid_decay_conductance_w_k,
                    gl * gg / (gl + gg),
                    2e-12,
                );
            }
        }
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - gas_before);
        assert!(exchange.heat_from_gas_j > 0.0 && gas_change < 0.0);
        let finite = report.finite_gas.as_ref().unwrap();
        assert!((gas_change + exchange.heat_from_gas_j).abs() <= finite.gas_energy_tolerance_j);
        let liquid_change = thermal.total_enthalpy_j().unwrap() - liquid_before;
        assert!(
            (liquid_change - exchange.heat_from_gas_j).abs()
                <= surface_controls().energy_tolerance_j
        );
        near(report.flow.external_heat_j, exchange.heat_from_gas_j, 2e-9);
        assert!(
            finite.combined_thermal_balance_defect_j.abs()
                <= surface_controls().energy_tolerance_j + finite.gas_energy_tolerance_j
        );
        assert_eq!(exchange.heat_from_radiation_j, 0.0);
        near(thermal.flow().ledger_mass_kg().unwrap(), mass_before, 2e-12);
    });
}
