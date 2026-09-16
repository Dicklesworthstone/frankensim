//! Conductivity uncertainty changes the actual material declaration consumed
//! by cooling-network, not a surrogate scalar resistance or post-hoc interval.
//! Principal-value sampling keeps the declared orthonormal axes fixed.

use super::*;

#[derive(Debug, Clone)]
pub(crate) enum Target {
    Uniform,
    Isotropic(String),
    Principal { material: String, axis: usize },
}

impl Target {
    /// None means the ordinary target parser owns this kind.
    pub(super) fn parse(value: &J) -> Result<Option<Self>> {
        let target = match field(value, "kind")?.as_str() {
            Some("solid-conductivity") => {
                object(value, &["kind"], "target")?;
                Self::Uniform
            }
            Some("material-conductivity") => {
                object(value, &["kind", "material"], "target")?;
                Self::Isotropic(string(field(value, "material")?, "material")?)
            }
            Some("material-principal-conductivity") => {
                object(value, &["kind", "material", "axis"], "target")?;
                Self::Principal {
                    material: string(field(value, "material")?, "material")?,
                    axis: integer(field(value, "axis")?, "principal axis", 2)?,
                }
            }
            _ => return Ok(None),
        };
        Ok(Some(target))
    }

    pub(super) fn name(&self) -> String {
        match self {
            Self::Uniform => "solid.conductivity_w_m_k".into(),
            Self::Isotropic(material) => format!("material[{material}].conductivity_w_m_k"),
            Self::Principal { material, axis } => format!("material[{material}].orthotropic.conductivity_w_m_k[{axis}]"),
        }
    }

    pub(super) fn render(&self) -> String {
        match self {
            Self::Uniform => "{\"kind\":\"solid-conductivity\"}".into(),
            Self::Isotropic(material) => format!("{{\"kind\":\"material-conductivity\",\"material\":{}}}", quote(material)),
            Self::Principal { material, axis } => format!("{{\"kind\":\"material-principal-conductivity\",\"material\":{},\"axis\":{axis}}}", quote(material)),
        }
    }

    pub(super) fn validate(&self, base: &J) -> Result<()> {
        match self {
            Self::Uniform => {
                let solid = field(base, "solid")?;
                if solid.get("materials").is_some() || solid.get("element_materials").is_some() {
                    return Err(bad("solid-conductivity requires the uniform scalar material mode"));
                }
                positive(field(solid, "conductivity_w_m_k")?, "solid conductivity")?;
            }
            Self::Isotropic(name) => {
                let material = material(base, name)?;
                if material.get("orthotropic").is_some() || material.get("conductivity_tensor_w_m_k").is_some()
                    || material.get("conductivity_curve").is_some()
                {
                    return Err(bad("material-conductivity requires an isotropic scalar material; do not scalarize a tensor or temperature law"));
                }
                positive(field(material, "conductivity_w_m_k")?, "material conductivity")?;
            }
            Self::Principal { material: name, axis } => {
                let material = material(base, name)?;
                if material.get("conductivity_w_m_k").is_some() || material.get("conductivity_tensor_w_m_k").is_some()
                    || material.get("conductivity_curve").is_some()
                {
                    return Err(bad("material-principal-conductivity requires an explicit orthotropic declaration"));
                }
                let values = array(field(field(material, "orthotropic")?, "conductivity_w_m_k")?, "principal conductivities", 3)?;
                if values.len() != 3 || *axis >= 3 {
                    return Err(bad("orthotropic material requires three principal conductivities and axis 0, 1 or 2"));
                }
                for value in values { positive(value, "principal conductivity")?; }
            }
        }
        Ok(())
    }

    pub(super) fn apply(&self, base: &mut J, value: f64) -> Result<()> {
        if !(value.is_finite() && value > 0.0) {
            return Err(model_failure("sampled conductivity must be finite and positive; samples are never clipped or skipped"));
        }
        match self {
            Self::Uniform => set_path_number(base, &["solid", "conductivity_w_m_k"], value),
            Self::Isotropic(name) => set_member_number(material_mut(base, name)?, "conductivity_w_m_k", value),
            Self::Principal { material, axis } => {
                let model = material_mut(base, material)?;
                let values = array_mut_path(model, &["orthotropic", "conductivity_w_m_k"])?;
                if values.len() != 3 || *axis >= 3 { return Err(bad("sample material principal axis disappeared")); }
                values[*axis] = J::Number { value, raw: value.to_string() };
                Ok(())
            }
        }
    }
}

fn material<'a>(base: &'a J, name: &str) -> Result<&'a J> {
    let solid = field(base, "solid")?;
    if solid.get("conductivity_w_m_k").is_some() {
        return Err(bad("named material uncertainty requires the element-material mode"));
    }
    field(solid, "element_materials")?;
    let rows = array(field(solid, "materials")?, "materials", 4096)?;
    let mut matches = rows.iter().filter(|row| row.str_field("name") == Some(name));
    let found = matches.next().ok_or_else(|| bad(format!("unknown uncertain material {name}")))?;
    if matches.next().is_some() { return Err(bad(format!("duplicate uncertain material {name}"))); }
    Ok(found)
}

fn material_mut<'a>(base: &'a mut J, name: &str) -> Result<&'a mut J> {
    // Re-admit before mutation; malformed identity must not select an arbitrary row.
    material(base, name)?;
    array_mut_path(base, &["solid", "materials"])?
        .iter_mut().find(|row| row.str_field("name") == Some(name))
        .ok_or_else(|| bad(format!("sample material {name} disappeared")))
}

#[cfg(test)]
mod tests {
    use super::*;

    const BASE: &str = r#"{"solid":{"materials":[{"name":"metal","conductivity_w_m_k":20},{"name":"board","orthotropic":{"principal_axes":[[0.6,0.8,0],[-0.8,0.6,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}}],"element_materials":["metal","board"]}}"#;

    #[test]
    fn principal_sampling_changes_one_eigenvalue_not_the_axes_or_assignment() {
        let mut base = J::parse(BASE).unwrap();
        let old = base.clone();
        let target = Target::parse(&J::parse(r#"{"kind":"material-principal-conductivity","material":"board","axis":1}"#).unwrap()).unwrap().unwrap();
        target.validate(&base).unwrap();
        target.apply(&mut base, 3.0).unwrap();
        let before = material(&old, "board").unwrap();
        let after = material(&base, "board").unwrap();
        assert_eq!(before.path(&["orthotropic", "principal_axes"]), after.path(&["orthotropic", "principal_axes"]));
        let values = after.path(&["orthotropic", "conductivity_w_m_k"]).unwrap().as_array().unwrap();
        assert_eq!(values.iter().map(J::as_f64).collect::<Vec<_>>(), vec![Some(20.0), Some(3.0), Some(1.0)]);
        assert_eq!(old.path(&["solid", "element_materials"]), base.path(&["solid", "element_materials"]));
        assert_eq!(material(&old, "metal").unwrap(), material(&base, "metal").unwrap());
        assert!(J::parse(&target.render()).is_ok());
    }

    #[test]
    fn material_modes_and_invalid_axes_refuse_without_mutating_the_model() {
        let mut base = J::parse(BASE).unwrap();
        assert!(Target::Uniform.validate(&base).is_err());
        assert!(Target::Isotropic("board".into()).validate(&base).is_err());
        assert!(Target::Principal { material: "metal".into(), axis: 0 }.validate(&base).is_err());
        assert!(Target::Isotropic("missing".into()).validate(&base).is_err());
        assert!(Target::parse(&J::parse(r#"{"kind":"material-principal-conductivity","material":"board","axis":3}"#).unwrap()).is_err());
        assert!(Target::parse(&J::parse(r#"{"kind":"material-conductivity","material":"metal","axis":0}"#).unwrap()).is_err());
        let curve = J::parse(r#"{"solid":{"materials":[{"name":"metal","conductivity_curve":{"temperature_k":[250,400],"conductivity_w_m_k":[17,2]}}],"element_materials":["metal"]}}"#).unwrap();
        assert!(Target::Isotropic("metal".into()).validate(&curve).is_err());
        let original = base.clone();
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Target::Isotropic("metal".into()).apply(&mut base, value).is_err());
            assert_eq!(base, original);
        }
    }

    #[test]
    fn uniform_and_named_scalar_targets_update_the_actual_operator_input() {
        let mut uniform = J::parse(r#"{"solid":{"conductivity_w_m_k":10}}"#).unwrap();
        Target::Uniform.validate(&uniform).unwrap();
        Target::Uniform.apply(&mut uniform, 12.0).unwrap();
        assert_eq!(uniform.path(&["solid", "conductivity_w_m_k"]).and_then(J::as_f64), Some(12.0));
        let mut heterogeneous = J::parse(BASE).unwrap();
        let target = Target::Isotropic("metal".into());
        target.validate(&heterogeneous).unwrap();
        target.apply(&mut heterogeneous, 25.0).unwrap();
        assert_eq!(material(&heterogeneous, "metal").unwrap().f64_field("conductivity_w_m_k"), Some(25.0));
    }
}
