//! Original one-link geometry and storage, independent of scalar contact work.
use super::super::adhesive_events::native::projected_position;
use super::{
    compatible, finite, invalid, limit, CoupledModalSystem, JkrNormalEventModalConfig,
    JkrNormalEventModalError, ModalAcousticTimeModel,
};

#[derive(Clone, Copy)]
pub(super) struct ElasticConnection {
    pub stiffness_n_per_m: f64,
    pub zero_force_indentation_m: f64,
    rest_extension_m: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct SpringObservation {
    pub extension_m: f64,
    pub energy_j: f64,
    pub reaction_n: f64,
}

pub(super) fn admit(
    network: &CoupledModalSystem,
    initial_indentation_m: f64,
) -> Result<(Vec<f64>, ElasticConnection), JkrNormalEventModalError> {
    if network.connections.len() != 1 || network.columns.len() != 1 {
        return Err(
            invalid("normal elastic adhesion needs exactly one original connection").into(),
        );
    }
    let link = &network.connections[0];
    if !link.stiffness_n_m.is_finite()
        || link.stiffness_n_m <= 0.0
        || link.damping_n_s_m != 0.0
        || link.constant_loss_reference_omega_rad_s.is_some()
    {
        return Err(
            invalid("normal elastic adhesion needs a positive undamped original spring").into(),
        );
    }
    let column = &network.columns[0];
    if column.len() != network.mode_count() || column.iter().any(|x| !x.is_finite()) {
        return Err(invalid("original spring column must cover every finite native mode").into());
    }
    let extension =
        finite(projected_position(network.components(), column)? - link.rest_extension_m)?;
    limit(
        "initial original normal spring force",
        finite(link.stiffness_n_m * extension)?.abs(),
        network.config.maximum_abs_connection_force_n,
    )?;
    Ok((
        column.clone(),
        ElasticConnection {
            stiffness_n_per_m: link.stiffness_n_m,
            zero_force_indentation_m: finite(initial_indentation_m - extension)?,
            rest_extension_m: link.rest_extension_m,
        },
    ))
}

impl ElasticConnection {
    pub(super) fn observe(
        self,
        models: &[ModalAcousticTimeModel],
        column: &[f64],
        indentation_m: f64,
        config: JkrNormalEventModalConfig,
    ) -> Result<SpringObservation, JkrNormalEventModalError> {
        let extension = finite(projected_position(models, column)? - self.rest_extension_m)?;
        compatible(
            "original normal spring extension",
            extension,
            finite(indentation_m - self.zero_force_indentation_m)?,
            config.coordinate_absolute_tolerance_m,
            config.coordinate_relative_tolerance,
        )?;
        Ok(SpringObservation {
            extension_m: extension,
            energy_j: finite(0.5 * self.stiffness_n_per_m * extension * extension)?,
            reaction_n: finite(-self.stiffness_n_per_m * extension)?,
        })
    }

    pub(super) fn storage_change(
        self,
        before: SpringObservation,
        after: SpringObservation,
    ) -> Result<f64, JkrNormalEventModalError> {
        Ok(finite(
            0.5 * self.stiffness_n_per_m
                * (after.extension_m - before.extension_m)
                * (after.extension_m + before.extension_m),
        )?)
    }
}
