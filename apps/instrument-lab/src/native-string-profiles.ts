import { el } from "./util";
import type { NativeStringProfiles } from "./native-wind-client.mjs";

/** Optional physical profiles remain original attached data, never browser-generated coefficients. */
export function createNativeStringProfiles(changed: () => void) {
  const fields = [
    ["bodyPath", "Body geometry/material profile"],
    ["stopPath", "Finger, fret and contact profile"],
    ["contactThermalPath", "Contact thermal profile"],
  ] as const;
  const inputs = fields.map(([key, label]) => ({key, input: el("select", {"aria-label": label}), label}));
  const element = el("fieldset", {class: "nw-gas-fields"},
    el("legend", {}, "String specimen profiles"),
    el("p", {class: "small muted"}, "Select any optional body, contact or thermal profiles from the attached files. All materials referenced inside these profiles must also be attached. Leave a choice empty to use the performance's existing description."),
    el("div", {class: "nw-gas-component"}, ...inputs.map(({input, label}) => el("label", {}, label, input))));
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
  refresh();
  return {
    element,
    read(): NativeStringProfiles | undefined {
      const selected = values();
      if (Object.values(selected).some(path => !paths.includes(path))) throw new Error("Select an attached file for each string specimen profile.");
      return Object.keys(selected).length ? selected : undefined;
    },
    identity: (): string => JSON.stringify(values()),
    setAttachments(next: readonly string[]): void { paths = next; refresh(); },
    setDisabled(disabled: boolean): void { element.disabled = disabled; },
    reset(): void { inputs.forEach(({input}) => { input.value = ""; }); },
  };
}
