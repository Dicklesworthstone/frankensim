//! Radiation controls modify the declared physical input, not the solved heat.
//! Emissivity and surroundings temperature are distinct from convective h and
//! air inlet temperature. Neither changes source footprints or surface faces.
use super::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum Target {
    Emissivity(String),
    AmbientTemperature(String),
}

impl Target {
    pub(super) fn parse(value: &J) -> Result<Option<Self>> {
        let kind = value.str_field("kind");
        if !matches!(kind, Some("radiation-emissivity" | "radiation-ambient-temperature")) {
            return Ok(None);
        }
        object(value, &["kind", "surface"], "radiation target")?;
        let surface = string(field(value, "surface")?, "radiation target surface")?;
        Ok(Some(if kind == Some("radiation-emissivity") {
            Self::Emissivity(surface)
        } else { Self::AmbientTemperature(surface) }))
    }

    fn surface(&self) -> &str {
        match self { Self::Emissivity(s) | Self::AmbientTemperature(s) => s }
    }
    fn key(&self) -> &'static str {
        match self { Self::Emissivity(_) => "emissivity", Self::AmbientTemperature(_) => "ambient_temperature_k" }
    }
    pub(super) fn name(&self) -> String {
        format!("radiation.surface[{}].{}", self.surface(), self.key())
    }
    pub(super) fn unit(&self) -> &'static str {
        match self { Self::Emissivity(_) => "1", Self::AmbientTemperature(_) => "K" }
    }
    pub(super) fn render(&self) -> String {
        let kind = match self { Self::Emissivity(_) => "radiation-emissivity",
            Self::AmbientTemperature(_) => "radiation-ambient-temperature" };
        format!("{{\"kind\":{},\"surface\":{}}}", quote(kind), quote(self.surface()))
    }

    fn row(&self, base: &J) -> Result<usize> {
        let rows = array_path(base, &["radiation", "surfaces"])?;
        let indices = rows.iter().enumerate().filter(|(_, row)| row.str_field("surface") == Some(self.surface()))
            .map(|(index, _)| index).collect::<Vec<_>>();
        if indices.len() != 1 {
            return Err(bad(format!("radiation target {} must identify exactly one declared patch", self.surface())));
        }
        Ok(indices[0])
    }
    fn admitted(&self, value: f64) -> bool {
        value.is_finite() && value > 0.0 && (!matches!(self, Self::Emissivity(_)) || value <= 1.0)
    }
    pub(super) fn validate(&self, base: &J) -> Result<()> {
        let index = self.row(base)?;
        let rows = array_path(base, &["radiation", "surfaces"])?;
        let value = number(field(&rows[index], self.key())?, self.key())?;
        if !self.admitted(value) { return Err(bad("base radiation target is outside its physical domain")); }
        Ok(())
    }
    pub(super) fn validate_support(&self, distribution: &Distribution) -> Result<()> {
        if let Distribution::Uniform { lo, hi } = *distribution {
            if !self.admitted(lo) || !self.admitted(hi) {
                return Err(bad("uniform radiation support must have emissivity in (0,1] and temperature above zero"));
            }
        }
        Ok(())
    }
    pub(super) fn apply(&self, base: &mut J, value: f64) -> Result<()> {
        if !self.admitted(value) {
            return Err(model_failure(format!("sampled {}={} leaves its physical domain; never clip or redraw", self.name(), value)));
        }
        let index = self.row(base)?;
        let rows = array_mut_path(base, &["radiation", "surfaces"])?;
        set_member_number(&mut rows[index], self.key(), value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/radiative-contact-hotspot.json"));

    #[test]
    fn samples_change_only_the_named_radiation_control() {
        let original = J::parse(BASE).unwrap();
        for target in [Target::Emissivity("first-face".into()), Target::AmbientTemperature("last-face".into())] {
            let mut base = original.clone();
            target.validate(&base).unwrap();
            let value = if target.unit() == "1" { 0.7 } else { 310.0 };
            target.apply(&mut base, value).unwrap();
            assert_eq!(base.get("solid"), original.get("solid"));
            assert_eq!(base.get("hydraulics"), original.get("hydraulics"));
            let index = target.row(&base).unwrap();
            assert_eq!(array_path(&base,&["radiation","surfaces"]).unwrap()[index].f64_field(target.key()), Some(value));
            let restored = Target::parse(&J::parse(&target.render()).unwrap()).unwrap().unwrap();
            assert_eq!(restored, target);
        }
    }

    #[test]
    fn invalid_support_and_samples_refuse_without_mutation() {
        let target = Target::Emissivity("first-face".into());
        assert!(target.validate_support(&Distribution::Uniform { lo: 0.8, hi: 1.1 }).is_err());
        assert!(target.validate_support(&Distribution::Uniform { lo: 0.0, hi: 0.8 }).is_err());
        assert!(target.validate_support(&Distribution::Uniform { lo: 0.1, hi: 1.0 }).is_ok());
        let original = J::parse(BASE).unwrap();
        for value in [0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
            let mut base = original.clone();
            assert!(target.apply(&mut base,value).is_err());
            assert_eq!(base,original);
        }
    }

    #[test]
    fn unknown_or_ambiguous_patch_never_selects_an_arbitrary_owner() {
        let mut base = J::parse(BASE).unwrap();
        assert!(Target::Emissivity("missing".into()).validate(&base).is_err());
        let rows = array_mut_path(&mut base,&["radiation","surfaces"]).unwrap();
        rows.push(rows[0].clone());
        assert!(Target::Emissivity("first-face".into()).validate(&base).is_err());
    }
}
