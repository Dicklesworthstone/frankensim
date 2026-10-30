import { el } from "./util";
import type { NativeStringProfiles, NativeCoatingStrainGuard } from "./native-wind-client.mjs";

/** Optional physical profiles remain original attached data, never browser-generated coefficients. */
export function createNativeStringProfiles(changed: () => void) {
  const fields = [
    ["bodyPath", "Body geometry/material profile"],
    ["stopPath", "Finger, fret and contact profile"],
    ["contactThermalPath", "Contact thermal profile"],
  ] as const;
  const inputs = fields.map(([key, label]) => ({key, input: el("select", {"aria-label": label}), label}));
  const coatingEnabled = el("input", {type: "checkbox"});
  const maximumBytes = el("input", {type: "text", inputmode: "numeric", "aria-label": "Coating strain maximum bytes"});
  const maximumWork = el("input", {type: "text", inputmode: "numeric", "aria-label": "Coating strain maximum work"});
  const coatingLimits = el("div", {class: "nw-gas-component"},
    el("label", {}, "Maximum bytes", maximumBytes), el("label", {}, "Maximum work", maximumWork));
  const refreshCoating = (): void => {
    coatingLimits.hidden = !coatingEnabled.checked;
    maximumBytes.disabled = maximumWork.disabled = !coatingEnabled.checked;
  };
  const element = el("fieldset", {class: "nw-gas-fields"},
    el("legend", {}, "String specimen profiles"),
    el("p", {class: "small muted"}, "Select any optional body, contact or thermal profiles from the attached files. All materials referenced inside these profiles must also be attached. Leave a choice empty to use the performance's existing description."),
    el("div", {class: "nw-gas-component"}, ...inputs.map(({input, label}) => el("label", {}, label, input))),
    el("label", {}, coatingEnabled, " Check coating strain during playback"),
    el("p", {class: "small muted"}, "For a coated curved or assembled shell body, check the sourced coating limits at every mechanical endpoint. Supply explicit memory and work budgets; the native renderer admits the body and coating material."),
    coatingLimits);
  let paths: readonly string[] = [];
  const values = (): NativeStringProfiles => Object.fromEntries(inputs
    .filter(({input}) => input.value).map(({key, input}) => [key, input.value]));
  const refresh = (): void => {
    for (const {input} of inputs) {
      const selected = input.value;
      input.replaceChildren(el("option", {value: ""}, "No additional profile"),
        ...paths.map(path => el("option", {value: path}, path)));
      if (selected && !paths.includes(selected)) input.append(el("option", {value: selected}, `Missing: ${selected}`));
      input.value = selected;
    }
  };
  inputs.forEach(({input}) => input.addEventListener("change", changed));
  coatingEnabled.addEventListener("change", () => { refreshCoating(); changed(); });
  [maximumBytes, maximumWork].forEach(input => input.addEventListener("input", changed));
  refresh();
  refreshCoating();
  return {
    element,
    read(): NativeStringProfiles | undefined {
      const selected = values();
      if (Object.values(selected).some(path => !paths.includes(path))) throw new Error("Select an attached file for each string specimen profile.");
      return Object.keys(selected).length ? selected : undefined;
    },
    readCoatingStrainGuard(): NativeCoatingStrainGuard | undefined {
      if (!coatingEnabled.checked) return undefined;
      if (!values().bodyPath) throw new Error("Select an attached body profile for coating strain checks.");
      const limits = [maximumBytes, maximumWork].map(input => {
        const text = input.value.trim();
        const value = Number(text);
        if (!/^[0-9]+$/.test(text) || !Number.isSafeInteger(value) || value <= 0)
          throw new Error("Coating strain memory and work budgets must be positive safe integers.");
        return value;
      });
      return {maximumBytes: limits[0]!, maximumWork: limits[1]!};
    },
    identity: (): string => JSON.stringify({profiles: values(), coatingEnabled: coatingEnabled.checked,
      maximumBytes: maximumBytes.value, maximumWork: maximumWork.value}),
    load(profiles?: NativeStringProfiles, guard?: NativeCoatingStrainGuard): void {
      inputs.forEach(({key, input}) => { input.value = profiles?.[key] ?? ""; });
      coatingEnabled.checked = guard !== undefined;
      maximumBytes.value = guard ? String(guard.maximumBytes) : "";
      maximumWork.value = guard ? String(guard.maximumWork) : "";
      refreshCoating();
    },
    setAttachments(next: readonly string[]): void { paths = next; refresh(); },
    setDisabled(disabled: boolean): void { element.disabled = disabled; },
    reset(): void {
      inputs.forEach(({input}) => { input.value = ""; });
      coatingEnabled.checked = false;
      maximumBytes.value = maximumWork.value = "";
      refreshCoating();
    },
  };
}
