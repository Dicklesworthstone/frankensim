//! Natural-strain events at fixed in-plane configuration and material memory.
//!
//! Each original point first changes its prescribed reference at fixed full F,
//! then restores transverse plane-stress equilibrium through the native
//! zero-duration update. This is an explicit ordered parameter/mechanics split.
//! Both endpoints retain the same viscous metrics and physical clock.
//!
//! Exact elastic volume needs pressure actuation when det(G) changes. With the
//! native compression-positive convention P = P_remainder - p cof(F), its
//! endpoint reference work is +p Delta det(G) per original volume. It is a
//! physical parameter input, kept separate from the assessed endpoint error.
//! The private off-constraint remainder evaluation is never an accepted solid.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::visco_hyper::ViscoHyperelastic;

use super::{
    Endpoint, FiniteViscoPlaneError, FiniteViscoPlaneSolid, FiniteViscoStepControls, PointHistory,
    PointUpdate, ProposedHistory, ViscoelasticLawField, checkpoint, gate, invalid, material_error,
};

/// Original-volume energy accounts for one zero-duration reference event [J].
/// These are increments; cumulative native point ledgers remain with each
/// retained state. No term represents relaxation heat or a caloric model.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteViscoPlaneReferenceEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed native reference change at fixed full deformation and memory.
    pub reference_work_j: f64,
    /// Exact-volume parameter actuation, +p_endpoint Delta det(G), integrated
    /// over original volume. Zero for finite-bulk plane stress. Physical total
    /// reference input is this plus `reference_work_j`, counted once.
    pub constraint_reference_work_j: f64,
    /// Native zero-duration transverse loading increment. The exact remainder
    /// ledger includes constraint actuation here; it is not independently an
    /// additional external mechanical input or a physical heat account.
    pub native_transverse_loading_work_j: f64,
    /// Actual full endpoint Piola contracted with the transverse F increment.
    /// Free transverse traction makes this zero up to the local solve accuracy.
    pub transverse_endpoint_work_j: f64,
    /// Fixed-F reference work + native loading work - physical storage change.
    pub material_balance_error_j: f64,
    /// Transverse endpoint work + constraint reference work - native loading.
    /// Signed endpoint/splitting error, separately gated and never made heat.
    pub endpoint_work_defect_j: f64,
}

/// A complete admissible local reference transition. The in-plane state can
/// acquire a global force imbalance; the next `advance` resolves the caller's
/// actual supports and loads with these same retained material histories.
#[derive(Debug)]
pub struct FiniteViscoPlaneReferenceChange<'a> {
    pub solid: FiniteViscoPlaneSolid<'a>,
    pub energy: FiniteViscoPlaneReferenceEnergy,
}

struct ReferencePoint {
    material: PointUpdate,
    fixed_reference_work_pa: f64,
    constraint_reference_work_pa: f64,
    transverse_endpoint_work_pa: f64,
}

fn change_point(
    target: &ViscoHyperelastic,
    previous: &ViscoHyperelastic,
    accepted: PointHistory<'_>,
    local: PlaneStressSettings,
) -> Result<ReferencePoint, MaterialError> {
    let f = accepted.gradient();
    match accepted {
        PointHistory::Full(state) => {
            let response = target.change_plane_stress_reference_from(
                previous,
                state,
                Some([f[2], f[5], f[8]]),
                local,
            )?;
            let material = response.material;
            Ok(ReferencePoint {
                material: PointUpdate {
                    stored_pa: material.state.stored,
                    history: ProposedHistory::Full(material.state),
                    first_piola_pa: material.first_piola_pa,
                    pressure_pa: None,
                    loading_work_pa: material.loading_work_pa,
                    dissipated_pa: material.dissipated_pa,
                    incremental_potential_pa: material.incremental_potential_pa,
                },
                fixed_reference_work_pa: response.reference_work_pa,
                constraint_reference_work_pa: response.constraint_reference_work_pa,
                transverse_endpoint_work_pa: response.transverse_endpoint_work_pa,
            })
        }
        PointHistory::Remainder(state) => {
            let response = target.change_incompressible_plane_stress_reference_from(
                previous,
                state,
                Some([f[2], f[5]]),
                local,
            )?;
            let material = response.material;
            Ok(ReferencePoint {
                material: PointUpdate {
                    stored_pa: material.state.stored,
                    history: ProposedHistory::Remainder(material.state),
                    first_piola_pa: response.plane_stress.full_piola,
                    pressure_pa: Some(response.plane_stress.pressure_pa),
                    loading_work_pa: material.loading_work_pa,
                    dissipated_pa: material.dissipated_pa,
                    incremental_potential_pa: material.incremental_potential_pa,
                },
                fixed_reference_work_pa: response.reference_work_pa,
                constraint_reference_work_pa: response.constraint_reference_work_pa,
                transverse_endpoint_work_pa: response.transverse_endpoint_work_pa,
            })
        }
    }
}

impl<'a> FiniteViscoPlaneSolid<'a> {
    /// Change prescribed natural strain while retaining in-plane displacement,
    /// original quadrature, thickness, material memory and elapsed time.
    ///
    /// Only the native reference deformation, relaxation times and local
    /// constitutive controls may change. Equilibrium/ordered branch energies
    /// and elastic domains must remain the same; the original native owner
    /// enforces this for every point. Transverse shear/thickness are solved on
    /// the target law. Every point must admit, and both global energy gates
    /// must pass, before a successor is returned.
    ///
    /// The fixed-F reference update and zero-time transverse equilibration are
    /// an ordered split. Exact pressure-volume actuation is physical reference
    /// input; the remaining signed endpoint error uses the caller's explicit
    /// allowance. Refining the prescribed reference path controls that error.
    /// This event does not resolve an unseen heating path, thermal transport,
    /// entropy or branch relaxation. Positive physical duration belongs to
    /// successive calls to `advance` with the actual endpoint loads/supports.
    pub fn change_reference<'b>(
        &self,
        target_material: ViscoelasticLawField<'b>,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoPlaneReferenceChange<'b>, FiniteViscoPlaneError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        if ![
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "finite nonnegative reference-event work allowances required",
            ));
        }
        if let ViscoelasticLawField::PerElement(laws) = target_material
            && laws.len() != self.mesh.elems.len()
        {
            return Err(invalid(
                "one target viscoelastic law per original element required",
            ));
        }
        let mut endpoint = Endpoint::new(self.mode, self.points.len(), self.stored_energy_j);
        let mut energy = FiniteViscoPlaneReferenceEnergy {
            stored_before_j: self.stored_energy_j,
            ..Default::default()
        };
        for (index, point) in self.points.iter().enumerate() {
            checkpoint(cx)?;
            let update = change_point(
                target_material.at(point.element),
                self.material.at(point.element),
                self.histories.at(index),
                self.local,
            )
            .map_err(|source| material_error(point, source))?;
            let volume_m3 = point.reference_area_m2 * self.reference_thickness_m;
            energy.reference_work_j += volume_m3 * update.fixed_reference_work_pa;
            energy.constraint_reference_work_j += volume_m3 * update.constraint_reference_work_pa;
            energy.transverse_endpoint_work_j += volume_m3 * update.transverse_endpoint_work_pa;
            if ![
                energy.reference_work_j,
                energy.constraint_reference_work_j,
                energy.transverse_endpoint_work_j,
            ]
            .iter()
            .all(|value| value.is_finite())
            {
                return Err(invalid("reference-event work is unrepresentable"));
            }
            endpoint.push(volume_m3, update.material)?;
        }
        if endpoint.energy.dissipated_j != 0.0 {
            return Err(invalid(
                "a zero-duration reference event cannot generate relaxation loss",
            ));
        }
        energy.stored_after_j = endpoint.energy.stored_after_j;
        energy.native_transverse_loading_work_j = endpoint.energy.native_loading_work_j;
        energy.material_balance_error_j = energy.reference_work_j
            + energy.native_transverse_loading_work_j
            - (energy.stored_after_j - energy.stored_before_j);
        energy.endpoint_work_defect_j = energy.transverse_endpoint_work_j
            + energy.constraint_reference_work_j
            - energy.native_transverse_loading_work_j;
        gate(
            "reference-event material balance",
            energy.material_balance_error_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        gate(
            "reference-event endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteViscoPlaneReferenceChange {
            solid: FiniteViscoPlaneSolid {
                mesh: self.mesh,
                material: target_material,
                mode: self.mode,
                reference_thickness_m: self.reference_thickness_m,
                local: self.local,
                displacement_m: Arc::clone(&self.displacement_m),
                points: Arc::clone(&self.points),
                offsets: Arc::clone(&self.offsets),
                histories: endpoint.histories,
                first_piola_pa: endpoint.first_piola_pa,
                pressure_pa: endpoint.pressure_pa,
                stored_energy_j: energy.stored_after_j,
                time_s: self.time_s,
            },
            energy,
        })
    }
}
