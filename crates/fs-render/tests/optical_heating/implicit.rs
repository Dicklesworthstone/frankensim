//! Real surface and volume absorption through the original implicit heat owner.
use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{
    StatefulSurfaceIlluminationSource, SurfaceHeatingBudget, SurfaceHeatingError,
    SurfaceHeatingTolerance, SurfaceIlluminationSource, SurfaceIncidence, SurfaceOpticalHeating,
    SurfacePowerResponse,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::ThermalLink;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionError,
    ImplicitConductionGraph, ImplicitConductionLaw, ImplicitConductionLinearControls,
    ImplicitConductionObservedError, ImplicitConductionState, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use std::sync::Arc;

// Reuse the original test caloric chart; supply its exact entropy and tangent
// plus an explicitly constant spatial conductivity. No second heat integrator.
struct NativeLaw;
static LAW: NativeLaw = NativeLaw;
impl ThermalLaw for NativeLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        SENSIBLE.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        SENSIBLE.maximum_temperature_slope()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some(
            4.0 * (self.temperature_k(b)? / self.temperature_k(a)?).ln(),
        ))
    }
}
impl ConductivityLaw for NativeLaw {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(40.0)
    }
}
impl ImplicitConductionLaw for NativeLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [0.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: SENSIBLE.maximum_temperature_slope(),
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

fn native(
    cells: Vec<DepositedCell>,
) -> (DepositionSimulation<'static>, ConductionObserver<'static>) {
    let linked = cells.len() == 2;
    let events = if linked {
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(11),
                b: CellId(12),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        }]
    } else {
        vec![]
    };
    let original = thermal(cells, &LAW, events);
    let geometry = if linked {
        vec![ConductiveLink {
            link_id: 1,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }]
    } else {
        vec![]
    };
    let conduction =
        ConductionObserver::new(&original, &[(MaterialId(7), &LAW)], &geometry, &[], || {
            Ok(())
        })
        .unwrap();
    (original, conduction)
}
fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..ImplicitConductionControls::default()
    }
}
fn provenance(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}
fn metal(k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        [
            380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
        ]
        .map(|w| ConductorIorSample::try_new(w, 1.0, k).unwrap()),
        provenance(231),
    )
    .unwrap()
}
fn segment(start_s: f64, end_s: f64, start: f64, end: f64) -> SpectralPowerSegment<1> {
    SpectralPowerSegment {
        start_s,
        end_s,
        start_power_w: [start],
        end_power_w: [end],
    }
}
fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 1,
        max_total_power_segments: 4,
        max_source_evaluations: 64,
    }
}
fn tolerance() -> SurfaceHeatingTolerance {
    SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    }
}

#[test]
fn surface_pulses_heat_stiff_chain_and_graph_with_exact_fluence_and_robin_loss() {
    let mut second = cell(12, 1.0, 1.0, 0.0);
    second.geometry.start_m = [1.0, 0.0, 0.0];
    second.geometry.end_m = [2.0, 0.0, 0.0];
    second.boundary.htc_w_m2_k = 1.0; // Five uncovered m²: G_bath=5 W/K.
    let (original, conduction) = native(vec![cell(11, 1.0, 1.0, 0.0), second]);
    let response =
        SurfacePowerResponse::from_conductor([580.0], &metal(2.0), [1.0], 1.0, [0.5]).unwrap();
    close(response.absorbed_fraction()[0], 0.5);
    let input = SurfaceOpticalHeating::try_new(
        &[SurfaceIlluminationSource {
            source_id: 1,
            thermal_cell: CellId(11),
            material_id: MaterialId(7),
            temperature_bounds_k: (300.0, 600.0),
            response,
            power: SpectralPowerSchedule::try_new(
                &[
                    segment(0.25, 1.25, 8.0, 24.0),
                    segment(1.5, 2.0, 12.0, 12.0),
                ],
                2,
            )
            .unwrap(),
        }],
        &original,
        budget(),
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&original, &conduction, &[(MaterialId(7), &LAW)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &original,
        &conduction,
        &[(MaterialId(7), &LAW)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (hot, optics, audit) = chain
        .advance_observed_to(
            &original,
            &input,
            2.5,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    // Independent closed-form 2x2 BE resolvent: C=4, G_link=40, G_bath=5.
    // The pulse interval of 1 s is over fourteen explicit monotone limits.
    let (mut x, mut y, mut bath_loss) = (0.0, 0.0, 0.0);
    for (dt, absorbed_w) in [(0.25, 0.0), (1.0, 8.0), (0.25, 0.0), (0.5, 6.0), (0.5, 0.0)] {
        let a = 4.0 + dt * 40.0;
        let d = 4.0 + dt * 45.0;
        let g = dt * 40.0;
        let determinant = a * d - g * g;
        let rhs = [4.0 * x + dt * absorbed_w, 4.0 * y];
        x = (d * rhs[0] + g * rhs[1]) / determinant;
        y = (g * rhs[0] + a * rhs[1]) / determinant;
        bath_loss += dt * 5.0 * y;
    }
    close(hot.temperature_k(CellId(11)).unwrap(), 300.0 + x);
    close(hot.temperature_k(CellId(12)).unwrap(), 300.0 + y);
    close(optics.energy().incident_j[0], 22.0);
    close(optics.energy().reflected_j[0], 11.0);
    close(optics.energy().transmitted_j[0], 0.0);
    close(optics.energy().absorbed_j[0], 11.0);
    close(optics.energy().thermal_input_j, 11.0);
    close(optics.energy().residual_j[0], 0.0);
    close(audit.energy.external_heat_j, 11.0 - bath_loss);
    close(audit.energy.stored_energy_change_j, 4.0 * (x + y));
    assert_eq!(audit.energy.substeps, 5);
    assert_eq!(optics.source_evaluations(), 5);
    let (graph_hot, graph_optics, graph_audit) = graph
        .advance_observed_to(
            &original,
            &input,
            2.5,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    for id in [CellId(11), CellId(12)] {
        close(
            graph_hot.temperature_k(id).unwrap(),
            hot.temperature_k(id).unwrap(),
        );
    }
    assert_eq!(graph_optics.energy(), optics.energy());
    close(
        graph_audit.energy.external_heat_j,
        audit.energy.external_heat_j,
    );
    let (split, split_optics, _) = chain
        .advance_observed_to(
            &original,
            &input,
            1.25,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let (restarted, restarted_optics, _) = chain
        .advance_observed_to(
            &split,
            &split_optics,
            2.5,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(restarted.cells(), hot.cells());
    assert_eq!(restarted_optics.energy(), optics.energy());
    assert_eq!(restarted_optics.model_identity(), input.model_identity());
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
}

#[test]
fn stateful_surface_imex_uses_accepted_start_optics_and_rolls_back_late_refusals() {
    let (original, conduction) = native(vec![cell(11, 1.0, 1.0, 0.0)]);
    let chain =
        ImplicitConductionChain::new(&original, &conduction, &[(MaterialId(7), &LAW)], || Ok(()))
            .unwrap();
    let make_input = |upper_k: f64| {
        let law = SurfaceConductorLaw::try_temperature(
            provenance(232),
            &[(300.0, metal(1.0)), (upper_k, metal(1.0 + upper_k - 300.0))],
            2,
        )
        .unwrap();
        SurfaceOpticalHeating::try_new_with_stateful(
            &[],
            &[StatefulSurfaceIlluminationSource {
                source_id: 1,
                thermal_cell: CellId(11),
                material_id: MaterialId(7),
                temperature_bounds_k: (300.0, 600.0),
                incidence: SurfaceIncidence {
                    wavelengths_nm: [580.0],
                    incident_eta: [1.0],
                    incidence_cosine: 1.0,
                    s_power_fraction: [0.5],
                },
                optics: Arc::new(law),
                power: SpectralPowerSchedule::try_new(&[segment(0.0, 1.0, 8.0, 8.0)], 1).unwrap(),
            }],
            &original,
            budget(),
            2,
            tolerance(),
            || Ok(()),
        )
        .unwrap()
    };
    let input = make_input(302.0);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut complete_checkpoints = 0;
    let (hot, optics, audit) = chain
        .advance_observed_to(
            &original,
            &input,
            1.0,
            controls(0.5),
            &mut workspace,
            || {
                complete_checkpoints += 1;
                Ok(())
            },
        )
        .unwrap();
    // n=1, k=1+(T-300): A=4/(4+k²), evaluated at each accepted start.
    let first_rise = 0.8;
    let final_rise = first_rise + 4.0 / (4.0 + 1.8 * 1.8);
    close(hot.temperature_k(CellId(11)).unwrap(), 300.0 + final_rise);
    close(optics.energy().incident_j[0], 8.0);
    close(optics.energy().absorbed_j[0], 4.0 * final_rise);
    close(optics.energy().reflected_j[0], 8.0 - 4.0 * final_rise);
    close(
        audit.energy.external_heat_j,
        optics.energy().thermal_input_j,
    );
    assert_eq!(optics.source_evaluations(), 2);
    assert!(optics.response_evaluations() > 0);

    // First half-step is supported; the second accepted endpoint leaves the
    // optical chart. No first-half thermal or optical prefix may escape.
    let bounded = make_input(300.9);
    assert!(matches!(
        chain.advance_observed_to(
            &original,
            &bounded,
            1.0,
            controls(0.5),
            &mut workspace,
            || Ok(())
        ),
        Err(ImplicitConductionObservedError::Observer(
            SurfaceHeatingError::TemperatureDomain { .. }
        ))
    ));
    let mut checkpoints = 0;
    let cancelled = chain.advance_observed_to(
        &original,
        &input,
        1.0,
        controls(0.5),
        &mut workspace,
        || {
            checkpoints += 1;
            if checkpoints == complete_checkpoints - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(ImplicitConductionObservedError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        )) | Err(ImplicitConductionObservedError::Observer(
            SurfaceHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    for unchanged in [&input, &bounded] {
        assert_eq!(unchanged.time_s(), 0.0);
        assert_eq!(unchanged.energy().thermal_input_j, 0.0);
        assert_eq!(unchanged.energy().incident_j, [0.0]);
        assert_eq!(unchanged.source_evaluations(), 0);
        assert_eq!(unchanged.response_evaluations(), 0);
    }
}

#[test]
fn deposited_beer_lambert_absorption_enters_implicit_bath_balance_once() {
    with_cx(|cx| {
        let mut receiver = cell(11, 1.0, 1.0, 2.0);
        receiver.boundary.htc_w_m2_k = 0.5; // Six m²: G_bath=3 W/K.
        let (original, conduction) = native(vec![receiver]);
        let deposition = report(
            cx,
            &[1.0],
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            [580.0],
            [12.0],
        );
        let input =
            OpticalHeating::try_new(deposition, &[(0, CellId(11))], &original, 1e-11, || Ok(()))
                .unwrap();
        let chain =
            ImplicitConductionChain::new(&original, &conduction, &[(MaterialId(7), &LAW)], || {
                Ok(())
            })
            .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (hot, optics, audit) = chain
            .advance_observed_to(
                &original,
                &input,
                8.0,
                controls(8.0),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        let absorbed_w = 12.0 * -(-0.7_f64).exp_m1();
        let rise = 8.0 * (2.0 + absorbed_w) / (4.0 + 8.0 * 3.0);
        close(hot.temperature_k(CellId(11)).unwrap(), 300.0 + rise);
        close(audit.energy.external_heat_j, 4.0 * rise);
        close(audit.energy.stored_energy_change_j, 4.0 * rise);
        close(optics.energy().incident_j[0], 96.0);
        close(optics.energy().transmitted_j[0], 96.0 * (-0.7_f64).exp());
        close(optics.energy().absorbed_j[0], 8.0 * absorbed_w);
        close(optics.energy().scattered_j[0], 0.0);
        close(optics.energy().thermal_input_j, 8.0 * absorbed_w);
        close(
            audit.energy.external_heat_j,
            8.0 * (2.0 + absorbed_w - 3.0 * rise),
        );
        assert_eq!(audit.energy.substeps, 1);
        assert_eq!(hot.cells()[&CellId(11)].boundary.source_w, 2.0);
        assert_eq!(
            optics.power().medium_identity,
            input.power().medium_identity
        );
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(input.energy().thermal_input_j, 0.0);
    });
}
