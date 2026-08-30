import SwiftUI

@main
struct FrankenSimApp: App {
    var body: some Scene {
        WindowGroup {
            SimulationStudioView()
                .background(CatalystWindowFreedom())
#if targetEnvironment(macCatalyst)
                .frame(minWidth: 760, minHeight: 560)
#endif
        }
#if targetEnvironment(macCatalyst)
        .defaultSize(width: 1420, height: 920)
        .windowResizability(.contentMinSize)
#endif
        .commands {
            SidebarCommands()
            CommandMenu("Simulation") {
                Button("Run Current Kernel") {
                    NotificationCenter.default.post(name: .runSimulation, object: nil)
                }
                .keyboardShortcut("r", modifiers: .command)
            }
        }
    }
}

extension Notification.Name {
    static let runSimulation = Notification.Name("FrankenSim.runSimulation")
}
