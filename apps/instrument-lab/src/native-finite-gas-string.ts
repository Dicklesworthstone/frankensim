import { el } from "./util";
import type { NativeFiniteGasStringSources } from "./native-wind-client.mjs";

/** Bind original attached files; the native runner owns materials and gas state. */
export function createNativeFiniteGasStringControls(changed: () => void) {
  const fields = [["solidPath", "String material pack"], ["gasRecordPath", "Initial gas record"]] as const;
  const inputs = fields.map(([key, label]) => ({key, label, input:el("select", {"aria-label":label})}));
  const element = el("fieldset", {class:"nw-gas-fields"},
    el("legend", {}, "Heated string and finite atmosphere"),
    el("p", {class:"small muted"}, "Choose the attached solid material pack and initial gas record. Attach every sleeve and nozzle reservoir named in the scenario, then associate its species, NASA-9 and transport packs below with the frozen phase policy. The scenario must declare a 48 kHz clock."),
    el("div", {class:"nw-gas-component"}, ...inputs.map(({input, label}) => el("label", {}, label, input))));
  let paths: readonly string[] = [];
  const values = (): NativeFiniteGasStringSources => ({solidPath:inputs[0]!.input.value, gasRecordPath:inputs[1]!.input.value});
  const refresh = (): void => {
    for (const {input} of inputs) {
      const selected = input.value;
      input.replaceChildren(el("option", {value:""}, "Select an attached file…"),
        ...paths.map(path => el("option", {value:path}, path)));
      if (selected && !paths.includes(selected)) input.append(el("option", {value:selected}, `Missing: ${selected}`));
      input.value = selected;
    }
  };
  inputs.forEach(({input}) => input.addEventListener("change", changed));
  refresh();
  return {
    element,
    read(): NativeFiniteGasStringSources {
      const selected = values();
      if (Object.values(selected).some(path => !path || !paths.includes(path)))
        throw new Error("Select attached string-material and initial-gas files for the heated string.");
      return selected;
    },
    identity: (): string => JSON.stringify(values()),
    load(sources?: NativeFiniteGasStringSources): void {
      inputs.forEach(({key, input}) => { input.value = sources?.[key] ?? ""; });
    },
    setAttachments(next: readonly string[]): void { paths = next; refresh(); },
    setDisabled(disabled: boolean): void { element.disabled = disabled; },
    reset(): void { inputs.forEach(({input}) => { input.value = ""; }); },
  };
}
