#!/usr/bin/env python3
"""Offline FrankenTextiles -> existing FrankenSim matdb-pack source adapter.

No simulation dependency, network fetch, new database, or invented material law.
Input: ftx-textile-catalog/1 or one authored Material record. Output: the existing
material-tsv-v1 source plus manifest, compiled by `cargo run -p xtask -- matdb-pack`.
Original JSON bytes are retained verbatim and also bound into the source envelope.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

MAX_BYTES = 8 * 1024 * 1024
# name -> (native property, SI unit). Deliberately no bulk->sheet, coupon->sheet,
# permeability-test->Darcy, skin-friction->arbitrary-interface, or mixture rules.
PROPERTIES = {
    "areal_density": ("areal_density", "kg/m2"),
    "thickness": ("thickness", "m"), "width": ("width", "m"),
    "bulk_density": ("textile_bulk_density", "kg/m3"),
    "warp_stiffness": ("warp_stiffness", "N/m"),
    "weft_stiffness": ("weft_stiffness", "N/m"),
    "poisson_warp_weft": ("poisson_warp_weft", "1"),
    "shear_stiffness": ("shear_stiffness", "N/m"),
    "bending_rigidity": ("bending_rigidity", "N*m"),
    "bending_warp": ("bending_warp", "N*m"),
    "bending_weft": ("bending_weft", "N*m"),
    "bending_coupling": ("bending_coupling", "N*m"),
    "bending_twist": ("bending_twist", "N*m"),
    "specimen_flexural_rigidity": ("specimen_flexural_rigidity", "N*m2"),
    "friction_skin": ("friction_skin", "1"), "friction_self": ("friction_self", "1"),
    "thermal_conductivity": ("thermal_conductivity", "W/(m*K)"),
    "thermal_resistance": ("thermal_resistance", "m2*K/W"),
    "specific_heat": ("specific_heat_capacity", "J/(kg*K)"),
    "emissivity": ("emissivity", "1"),
    "air_permeability": ("air_permeability", "m/s"),
    "vapor_resistance": ("vapor_resistance", "m2*Pa/W"),
    "fiber_diameter": ("fiber_diameter", "m"), "staple_length": ("staple_length", "m"),
    "yarn_linear_density": ("yarn_linear_density", "kg/m"),
    "thread_density_warp": ("thread_density_warp", "1/m"),
    "thread_density_weft": ("thread_density_weft", "1/m"),
}
CONVERSIONS = {
    "g/m2": ("kg/m2", .001), "oz/yd2": ("kg/m2", .028349523125 / .9144**2),
    "mm": ("m", .001), "cm": ("m", .01), "um": ("m", 1e-6), "in": ("m", .0254),
    "g/cm3": ("kg/m3", 1000.0), "kN/m": ("N/m", 1000.0),
    "mN*m": ("N*m", .001), "lbf*in2": ("N*m2", 4.4482216152605 * .0254**2),
    "%": ("1", .01), "clo": ("m2*K/W", .155), "tog": ("m2*K/W", .1),
    "tex": ("kg/m", 1e-6), "denier": ("kg/m", 1e-3 / 9000.0),
    "1/in": ("1/m", 1.0 / .0254), "1/cm": ("1/m", 100.0), "1/mm": ("1/m", 1000.0),
}


def wire_unit(unit: str) -> str:
    # fs-qty deliberately uses a strict left-to-right unit grammar.
    return {"W/(m*K)":"W/m/K", "J/(kg*K)":"J/kg/K", "1/m":"m^-1"}.get(unit,unit)


def fail(message: str) -> None:
    raise ValueError(message)


def text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 8192:
        fail(f"{label}: nonempty bounded text required")
    return value


def number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        fail(f"{label}: finite number required")
    result = float(value)
    if not math.isfinite(result):
        fail(f"{label}: nonfinite number")
    return result


def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def metadata(value: Any) -> str:
    # Escaping is intentional: JSON strings cannot inject TSV records.
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), allow_nan=False)


def conditions(value: Any) -> dict[str, str]:
    if not isinstance(value, dict) or len(value) > 64:
        fail("conditions: bounded text dictionary required")
    for k, v in value.items():
        text(k, "condition key"); text(v, "condition value")
    return value


def source_for(quantity: dict[str, Any], sources: dict[str, Any], allow_estimates: bool) -> dict[str, Any]:
    evidence = quantity.get("evidence")
    if evidence not in ("measured", "published", "estimate"):
        fail("unknown or missing evidence class")
    if evidence == "estimate" and not allow_estimates:
        fail("estimated properties require --allow-estimates")
    source_id = text(quantity.get("source"), "source id")
    source = sources.get(source_id)
    if not isinstance(source, dict):
        fail(f"missing source {source_id!r}")
    for key in ("reference", "accessed_on", "note"):
        text(source.get(key), f"source {key}")
    text(quantity.get("method"), "method")
    conditions(quantity.get("conditions", {}))
    return source


def normalized(value: Any, source_unit: str, target_unit: str) -> float:
    v = number(value, "property value")
    unit, scale = CONVERSIONS.get(source_unit, (source_unit, 1.0))
    if unit != target_unit:
        fail(f"unsupported conversion {source_unit!r} -> {target_unit!r}; dimensions are not interchangeable")
    result = v * scale
    if not math.isfinite(result) or (v != 0.0 and result == 0.0):
        fail("unit conversion is unrepresentable")
    return result


def render(raw: bytes, material_id: str, *, citation: str, license_text: str,
           redistribution_terms: str, allow_estimates: bool = False,
           interpolate_temperature: bool = False, selected: set[str] | None = None,
           expected_sha256: str | None = None,
           nonlinear_membrane: bool = False) -> tuple[bytes, bytes, dict[str, Any]]:
    """Validate and render completely before creating any destination files."""
    if len(raw) > MAX_BYTES:
        fail("catalog exceeds 8 MiB")
    digest = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and expected_sha256.lower() != digest:
        fail("source SHA-256 mismatch")
    for value in (citation, license_text, redistribution_terms):
        text(value, "publication terms")
        if any(c in value for c in "\t\r\n\0"):
            fail("manifest citation/license/terms must be single-line text")
    try:
        doc = json.loads(raw, object_pairs_hook=object_pairs, parse_constant=lambda x: fail(f"invalid JSON number {x}"))
    except (UnicodeError, RecursionError, json.JSONDecodeError) as exc:
        fail(f"invalid bounded JSON: {exc}")
    if not isinstance(doc, dict):
        fail("expected an authored Material or ftx-textile-catalog/1")
    if "schema" in doc:
        if doc["schema"] != "ftx-textile-catalog/1" or set(doc) != {"schema", "materials"}:
            fail("only ftx-textile-catalog/1 is admitted; resolve legacy/snapshot formats in FrankenTextiles first")
        records = doc["materials"]
    else:
        records = [doc]
    if not isinstance(records, list) or not 1 <= len(records) <= 4096:
        fail("expected 1..4096 material records")
    seen: set[str] = set(); matches = []
    for record in records:
        if not isinstance(record, dict): fail("invalid material record")
        key = text(record.get("id"), "material id")
        if key in seen: fail("duplicate material id")
        seen.add(key)
        if key == material_id: matches.append(record)
    if len(matches) != 1:
        fail(f"exact material id {material_id!r} not present")
    material = matches[0]
    text(material.get("name"), "material name")
    if material.get("form") not in ("sheet", "fiber", "yarn", "woven", "knit", "nonwoven", "hide", "fill", "laminate"):
        fail("unknown material form")
    props, sources = material.get("properties", {}), material.get("sources", {})
    if not isinstance(props, dict) or len(props) > 256 or not isinstance(sources, dict) or len(sources) > 128:
        fail("bounded properties and sources dictionaries required")
    if selected is not None and not selected <= PROPERTIES.keys():
        fail("--properties names use the supported ORIGINAL FrankenTextiles property vocabulary")
    chosen = {k: v for k, v in props.items() if k in PROPERTIES and (selected is None or k in selected)}
    membrane_keys = {"warp_stiffness", "weft_stiffness", "poisson_warp_weft", "shear_stiffness"}
    if nonlinear_membrane:
        if selected is not None and selected & membrane_keys:
            fail("nonlinear membrane replaces scalar membrane properties; do not explicitly request both")
        if material.get("membrane_relaxation") is not None:
            fail("a history-dependent membrane cannot become an equilibrium stress curve")
        if material.get("membrane_response") is None:
            fail("--nonlinear-membrane requires an explicit membrane_response")
        chosen = {k:v for k,v in chosen.items() if k not in membrane_keys}
    mechanical = set(chosen) & membrane_keys
    if mechanical and (material.get("membrane_response") is not None or material.get("membrane_relaxation") is not None):
        fail("nonlinear/history-dependent membrane must not be flattened to equilibrium scalars; select nonmechanical properties explicitly")
    if material["form"] in ("fiber", "yarn", "fill") and (mechanical or nonlinear_membrane):
        fail("fiber/yarn/fill cannot be imported as an effective sheet law")
    curvature = [v for k, v in props.items() if k.startswith("preferred_curvature_")]
    if any(number(v.get("value"), "preferred curvature") != 0.0 for v in curvature):
        if any(k.startswith("bending_") for k in chosen):
            fail("nonzero preferred curvature needs an explicit rest-shape consumer, not a flat-plate substitution")
    lines = ["frankensim.matdb-source.v1", f"# FrankenTextiles input SHA-256 {digest}",
             f"# exact material id {metadata(material_id)}", "# Scalar constants are frozen at authored specimen/conditioning, not a T/rate/composition law."]
    # Native source-file hashes now transitively bind EVERY original JSON byte,
    # including composition, methods, bounds, sample SD and unconsumed fields.
    encoded = base64.b64encode(raw).decode("ascii")
    lines += ["# original_json_base64 " + encoded[i:i+8192] for i in range(0, len(encoded), 8192)]
    imported: list[str] = []
    def claim(key: str, quantity: dict[str, Any], value_fields: list[str], validity: list[str],
              native_property: str | None = None) -> None:
        source = source_for(quantity, sources, allow_estimates)
        index = len(imported); obs, cid = f"o{index}", f"c{index}"
        qualifier = {"record_id":material_id,"form":material["form"],"quantity":quantity,"source":source,
                     "source_sha256":digest,"projection":"frozen specimen; no implicit construction homogenization"}
        lines.append("\t".join(["observation",obs,metadata(material_id),metadata(quantity["method"]),metadata(qualifier)]))
        # Published labels/estimates are not promoted to measured observation backing.
        refs = obs if quantity["evidence"] == "measured" else "-"
        lines.append("\t".join([value_fields[0],cid,refs,native_property or PROPERTIES[key][0],*value_fields[1:]]))
        lines.append("\t".join(["uncertainty",cid,"unstated","-","-","-","-"]))
        lines.append("\t".join(["validity",cid,"ftx_characterization","0","0","1"]))
        if validity: lines.append("\t".join(["validity",cid,*validity]))
        imported.append(key)
    for key, q in sorted(chosen.items()):
        if not isinstance(q, dict): fail(f"{key}: quantity object required")
        value = normalized(q.get("value"), text(q.get("unit"), "unit"), PROPERTIES[key][1])
        if q.get("bounds") is not None:
            b = q["bounds"]
            if not isinstance(b, list) or len(b) != 2 or not number(b[0], "lower bound") <= number(q["value"], "value") <= number(b[1], "upper bound"):
                fail(f"{key}: invalid source bounds")
        claim(key,q,["scalar",repr(value),wire_unit(PROPERTIES[key][1]),"constant"],[])
    if nonlinear_membrane:
        response = material["membrane_response"]
        if not isinstance(response, dict) or response.get("model") != "tabulated_green_strain/1":
            fail("unsupported nonlinear membrane model; require second-Piola resultants versus Green strain")
        if response.get("interpolation") != "piecewise_linear":
            fail("only explicitly authored piecewise_linear response is supported; cubic/default interpolation is not silently changed")
        if set(response) != {"model", "interpolation", "warp", "weft", "shear", "coupling"}:
            fail("nonlinear membrane has missing or unsupported fields")
        for direction in ("warp", "weft", "shear"):
            curve = response[direction]
            if not isinstance(curve, dict) or set(curve) - {"points", "source", "evidence", "method", "conditions"}:
                fail("invalid or unsupported stress-curve fields")
            source_for(curve, sources, allow_estimates)
            points = curve.get("points")
            if not isinstance(points, list) or not 3 <= len(points) <= 128:
                fail("stress curve requires 3..128 signed strain/resultant pairs")
            knots = []
            for pair in points:
                if not isinstance(pair, list) or len(pair) != 2:
                    fail("stress curve needs [Green strain, second-Piola resultant N/m] pairs")
                x, y = (number(v, "stress-curve coordinate") for v in pair)
                if knots and (x <= knots[-1][0] or y <= knots[-1][1]):
                    fail("stress and strain knots must both strictly increase")
                if knots:
                    dx, dy = x-knots[-1][0], y-knots[-1][1]
                    if not math.isfinite(dx) or not math.isfinite(dy) or not math.isfinite(dy/dx) or dy/dx <= 0:
                        fail("unrepresentable stress-curve segment")
                knots.append((x,y))
            if (0.0,0.0) not in knots or knots[0][0] >= 0 or knots[-1][0] <= 0:
                fail("stress curve must include a stress-free origin and both strain signs")
            q = dict(curve, model=response["model"], interpolation="piecewise_linear")
            wire = ",".join(f"{x!r}:{y!r}" for x,y in knots)
            claim(f"membrane_response.{direction}",q,
                  ["curve","sheet_green_strain","1","N/m",wire,"linear"],
                  ["sheet_green_strain",repr(knots[0][0]),repr(knots[-1][0]),"1"],
                  f"sheet_second_piola_{direction}")
        coupling = response["coupling"]
        if not isinstance(coupling, dict): fail("explicit source-backed cross coupling required")
        value = normalized(coupling.get("value"), text(coupling.get("unit"),"coupling unit"),"N/m")
        bounds = coupling.get("bounds")
        if bounds is not None and (not isinstance(bounds,list) or len(bounds)!=2
                or not number(bounds[0],"coupling lower") <= number(coupling["value"],"coupling") <= number(bounds[1],"coupling upper")):
            fail("invalid coupling source bounds")
        claim("membrane_response.coupling",coupling,["scalar",repr(value),"N/m","constant"],[],
              "sheet_green_strain_coupling")
    thermal = material.get("thermal_response")
    if thermal is not None:
        if not isinstance(thermal, dict): fail("thermal_response must be an object")
        specimen = text(thermal.get("specimen"), "thermal specimen")
        conditioning = text(thermal.get("conditioning"), "thermal conditioning")
        for field, key in (("conductivity","thermal_conductivity"),("specific_heat","specific_heat")):
            series = thermal.get(field)
            if series is None or (selected is not None and key not in selected): continue
            if key in props or (key == "thermal_conductivity" and "thermal_resistance" in props):
                fail("ambiguous scalar and temperature-dependent thermal law")
            if not isinstance(series, dict) or series.get("polynomial") is not None:
                fail("polynomial thermal fits need their original analytic owner; this adapter never silently samples them")
            source_for(series, sources, allow_estimates)
            sc = conditions(series.get("conditions", {}))
            for name, expected in (("specimen",specimen),("conditioning",conditioning)):
                if name in sc and sc[name] != expected: fail(f"conflicting thermal {name}")
            # A nominal supplier thickness cannot become the HFM specimen's geometry.
            for geom in ("areal_density", "thickness"):
                if geom in chosen:
                    gc = conditions(chosen[geom].get("conditions", {}))
                    if gc.get("specimen") != specimen or gc.get("conditioning") != conditioning:
                        fail("thermal companion geometry/mass must match specimen and conditioning; use --properties for transport-only records")
            points = series.get("points")
            if not isinstance(points, list) or not 2 <= len(points) <= 4096:
                fail("thermal table requires 2..4096 ordered measured/reported points")
            knots = []
            for p in points:
                if not isinstance(p, dict): fail("invalid thermal point")
                t = number(p.get("temperature_c"), "temperature_c") + 273.15
                v = normalized(p.get("mean"), text(series.get("unit"), "unit"), PROPERTIES[key][1])
                if t <= 0 or v <= 0 or (knots and t <= knots[-1][0]): fail("nonpositive or unordered thermal observations")
                if p.get("standard_deviation") is not None and number(p["standard_deviation"],"sample standard deviation") < 0:
                    fail("negative sample standard deviation")
                knots.append((t,v))
            if interpolate_temperature and not allow_estimates:
                fail("temperature interpolation is an estimate; --allow-estimates is required")
            q = dict(series); q["conditions"] = dict(sc,specimen=specimen,conditioning=conditioning)
            q["interpolation_admission"] = "explicit linear estimate" if interpolate_temperature else "exact knots only"
            encoded_points = ",".join(f"{t!r}:{v!r}" for t,v in knots)
            claim(key,q,["curve","T","K",wire_unit(PROPERTIES[key][1]),encoded_points,
                         "linear" if interpolate_temperature else "tabulated"],
                  ["T",repr(knots[0][0]),repr(knots[-1][0]),"K"])
    if not imported: fail("no selected supported observations; missing properties are never defaulted")
    if selected is not None and not selected <= set(imported):
        fail(f"requested properties are missing: {sorted(selected-set(imported))}")
    slug = re.sub(r"[^A-Za-z0-9_-]", "-", material_id)[:64]
    pack_id = f"ftx-{slug}-{digest[:16]}"
    manifest = "\n".join(["frankensim.matdb-manifest.v1", f"pack_id\t{pack_id}",
        f"redistribution\tpermitted\t{redistribution_terms}", f"citation\t{citation}; FTX record {metadata(material_id)}; input SHA-256 {digest}",
        f"license\t{license_text}","source\tfranken-textiles\tproperties.tsv\tmaterial-tsv-v1", ""])
    report = {"pack_id":pack_id,"source_sha256":digest,"material_id":material_id,"imported":imported,
              "source_only_properties":sorted(set(props)-set(imported)),
              "query": {"ftx_characterization":0},"temperature_interpolation":interpolate_temperature}
    if nonlinear_membrane:
        report["nonlinear_membrane"] = "tabulated_green_strain/1:piecewise_linear"
        report["query"]["sheet_green_strain"] = 0
    return manifest.encode(), ("\n".join(lines)+"\n").encode(), report


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path); p.add_argument("material_id"); p.add_argument("output", type=Path)
    p.add_argument("--citation", required=True); p.add_argument("--license", required=True)
    p.add_argument("--redistribution-terms", required=True)
    p.add_argument("--source-sha256"); p.add_argument("--allow-estimates", action="store_true")
    p.add_argument("--interpolate-temperature", action="store_true")
    p.add_argument("--nonlinear-membrane", action="store_true", help="retain an explicitly piecewise-linear Green-strain membrane law")
    p.add_argument("--properties", help="comma-separated ORIGINAL FTX property keys; omit for all supported observations")
    args = p.parse_args()
    try:
        with args.input.open("rb") as f: raw = f.read(MAX_BYTES+1)
        chosen = None if args.properties is None else set(args.properties.split(","))
        manifest, source, report = render(raw,args.material_id,citation=args.citation,license_text=args.license,
            redistribution_terms=args.redistribution_terms,allow_estimates=args.allow_estimates,
            interpolate_temperature=args.interpolate_temperature,selected=chosen,expected_sha256=args.source_sha256,nonlinear_membrane=args.nonlinear_membrane)
        args.output.mkdir(parents=False, exist_ok=False)
        # Manifest is last: an interrupted directory is not a compilable source pack.
        for name, contents in (("original.json",raw),("properties.tsv",source),
                               ("selection.json",(json.dumps(report,indent=2)+"\n").encode()),("manifest.tsv",manifest)):
            with (args.output/name).open("xb") as f: f.write(contents)
        print(json.dumps(report,indent=2))
        return 0
    except (OSError, ValueError, TypeError, OverflowError) as exc:
        print(f"textile import refused: {exc}",file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
