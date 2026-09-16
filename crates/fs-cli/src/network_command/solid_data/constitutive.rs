//! Constant anisotropic material declarations for the actual cooling producer.
//! All tensor entries and principal-axis rows use the mesh's Cartesian frame.
//! Constitutive validation and tensor construction belong to fs-conduction.

use super::*;
use fs_conduction::material::ConductivityTable;

#[derive(Debug, Clone)]
pub(super) enum Conductivity {
    Isotropic(f64),
    Tensor([[f64; 3]; 3]),
    Orthotropic { axes: [[f64; 3]; 3], principal: [f64; 3] },
}

impl Conductivity {
    pub(super) fn parse(row: &J) -> Result<Self> {
        let model = match (
            row.get("conductivity_w_m_k"),
            row.get("conductivity_tensor_w_m_k"),
            row.get("orthotropic"),
        ) {
            (Some(value), None, None) => Self::Isotropic(positive(value, "material.conductivity_w_m_k")?),
            (None, Some(value), None) => Self::Tensor(matrix(value, "conductivity_tensor_w_m_k")?),
            (None, None, Some(value)) => {
                object(value, &["principal_axes", "conductivity_w_m_k"], "material.orthotropic")?;
                let axes = matrix(get(value, "principal_axes")?, "orthotropic.principal_axes")?;
                let principal = triple(get(value, "conductivity_w_m_k")?, "orthotropic.conductivity_w_m_k")?;
                Self::Orthotropic { axes, principal }
            }
            _ => return Err(bad("each material requires exactly one of conductivity_w_m_k, conductivity_tensor_w_m_k, or orthotropic")),
        };
        // Refuse invalid constitutive data before mesh assignment or any solve.
        model.model()?;
        Ok(model)
    }

    pub(super) fn model(&self) -> Result<ConductivityModel> {
        match *self {
            Self::Isotropic(k) => ConductivityModel::isotropic_declared(k).map_err(producer),
            Self::Tensor(k) => ConductivityModel::constant_tensor(k).map_err(producer),
            Self::Orthotropic { axes, principal: k } => {
                let model = ConductivityModel::orthotropic(axes, [
                    ConductivityTable::declared(k[0]).map_err(producer)?,
                    ConductivityTable::declared(k[1]).map_err(producer)?,
                    ConductivityTable::declared(k[2]).map_err(producer)?,
                ]).map_err(producer)?;
                // Constant tables have no temperature dependence. This extra
                // admission catches an unrepresentable rotated tensor as well
                // as malformed axes, without implementing a rival SPD check.
                let tensor = model.tensor_at(0.0).map_err(producer)?;
                ConductivityModel::constant_tensor(tensor).map_err(producer)?;
                Ok(model)
            }
        }
    }

    /// Compatibility slot in Request: ignored whenever ElementMaterials exists.
    /// This is NOT an effective isotropic coefficient used by any solve.
    pub(super) fn inactive_scalar(&self) -> f64 {
        match *self {
            Self::Isotropic(k) => k,
            Self::Tensor(k) => k[0][0],
            Self::Orthotropic { principal, .. } => principal[0],
        }
    }

    pub(super) fn render_fields(&self) -> Result<String> {
        match *self {
            Self::Isotropic(k) => Ok(format!("\"conductivity_w_m_k\":{}", num(k)?)),
            Self::Tensor(k) => Ok(format!(
                "\"conductivity_tensor_w_m_k\":{},\"coordinate_frame\":\"mesh-cartesian\"",
                matrix_json(k)?,
            )),
            Self::Orthotropic { axes, principal } => Ok(format!(
                "\"orthotropic\":{{\"principal_axes\":{},\"conductivity_w_m_k\":{}}},\"coordinate_frame\":\"mesh-cartesian\",\"resolved_conductivity_tensor_w_m_k\":{}",
                matrix_json(axes)?, numbers(&principal)?,
                matrix_json(self.model()?.tensor_at(0.0).map_err(producer)?)?,
            )),
        }
    }
}

fn triple(value: &J, name: &str) -> Result<[f64; 3]> {
    let values = array(value, name, 3)?;
    if values.len() != 3 { return Err(bad(format!("{name} requires exactly three entries"))); }
    Ok([number(&values[0], name)?, number(&values[1], name)?, number(&values[2], name)?])
}

fn matrix(value: &J, name: &str) -> Result<[[f64; 3]; 3]> {
    let rows = array(value, name, 3)?;
    if rows.len() != 3 { return Err(bad(format!("{name} requires exactly three rows"))); }
    Ok([triple(&rows[0], name)?, triple(&rows[1], name)?, triple(&rows[2], name)?])
}

fn matrix_json(value: [[f64; 3]; 3]) -> Result<String> {
    Ok(format!("[{},{},{}]", numbers(&value[0])?, numbers(&value[1])?, numbers(&value[2])?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(text: &str) -> Result<Conductivity> { Conductivity::parse(&J::parse(text).unwrap()) }

    #[test]
    fn principal_axis_rows_generate_the_full_rotated_tensor() {
        let model = parsed(r#"{"orthotropic":{"principal_axes":[[0.6,0.8,0],[-0.8,0.6,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}}"#).unwrap();
        let expected = [[8.48, 8.64, 0.0], [8.64, 13.52, 0.0], [0.0, 0.0, 1.0]];
        let actual = model.model().unwrap().tensor_at(300.0).unwrap();
        for (row, expected) in actual.iter().zip(expected) {
            for (&value, expected) in row.iter().zip(expected) {
                assert!((value - expected).abs() < 1.0e-12);
            }
        }
        assert!(!model.model().unwrap().is_temperature_dependent());
        let receipt = J::parse(&format!("{{{}}}", model.render_fields().unwrap())).unwrap();
        assert_eq!(receipt.str_field("coordinate_frame"), Some("mesh-cartesian"));
        assert!(receipt.get("resolved_conductivity_tensor_w_m_k").is_some());
    }

    #[test]
    fn constant_tensor_and_orthotropic_forms_agree() {
        let tensor = parsed(r#"{"conductivity_tensor_w_m_k":[[2,0,0],[0,20,0],[0,0,1]]}"#).unwrap();
        let axes = parsed(r#"{"orthotropic":{"principal_axes":[[0,1,0],[-1,0,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}}"#).unwrap();
        assert_eq!(tensor.model().unwrap().tensor_at(300.0).unwrap(), axes.model().unwrap().tensor_at(300.0).unwrap());
    }

    #[test]
    fn malformed_ambiguous_and_nonphysical_tensors_refuse() {
        for text in [
            r#"{}"#,
            r#"{"conductivity_w_m_k":10,"conductivity_tensor_w_m_k":[[1,0,0],[0,1,0],[0,0,1]]}"#,
            r#"{"conductivity_tensor_w_m_k":[[1,0],[0,1]]}"#,
            r#"{"conductivity_tensor_w_m_k":[[1,2,0],[0,1,0],[0,0,1]]}"#,
            r#"{"conductivity_tensor_w_m_k":[[1,2,0],[2,1,0],[0,0,1]]}"#,
            r#"{"conductivity_tensor_w_m_k":[[1,0,0],[0,0,0],[0,0,1]]}"#,
            r#"{"orthotropic":{"principal_axes":[[1,0,0],[1,0,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}}"#,
            r#"{"orthotropic":{"principal_axes":[[2,0,0],[0,1,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}}"#,
            r#"{"orthotropic":{"principal_axes":[[1,0,0],[0,1,0],[0,0,1]],"conductivity_w_m_k":[20,-2,1]}}"#,
            r#"{"orthotropic":{"principal_axes":[[1,0,0],[0,1,0],[0,0,1]],"conductivity_w_m_k":[20,2,1],"angle":90}}"#,
        ] { assert!(parsed(text).is_err(), "accepted {text}"); }
    }
}
