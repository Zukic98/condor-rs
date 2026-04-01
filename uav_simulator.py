# uav_simulator.py
import subprocess
import json
import matplotlib.pyplot as plt
import os

# Configuration for the simulation
NUM_DRONE_CLUSTERS = 5  # Number of times we run the aggregation simulation
RUST_BINARY_PATH = "./target/release/examples/simulator_bridge"

def run_rust_aggregator():
    """
    Executes the compiled Rust binary and parses the JSON output.
    This simulates a Cluster Head (drone) aggregating signatures and 
    the Base Station verifying them.
    """
    if not os.path.exists(RUST_BINARY_PATH):
        print(f"Error: Rust binary not found at {RUST_BINARY_PATH}")
        print("Please run: cargo build --release --example simulator_bridge")
        return None

    try:
        # Run the binary and capture the output (stdout)
        result = subprocess.run(
            [RUST_BINARY_PATH], 
            capture_output=True, 
            text=True, 
            check=True
        )
        
        # Parse the JSON string printed by Rust
        data = json.loads(result.stdout.strip())
        return data

    except subprocess.CalledProcessError as e:
        print(f"Rust execution failed: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON from Rust: {e}")
        print(f"Raw output was: {result.stdout}")
        return None

def main():
    print(f"--- Starting UAV Network Simulation ({NUM_DRONE_CLUSTERS} iterations) ---")
    
    prove_times = []
    verify_times = []
    
    for i in range(1, NUM_DRONE_CLUSTERS + 1):
        print(f"[*] Simulating Cluster {i} aggregation...")
        metrics = run_rust_aggregator()
        
        if metrics and metrics.get("is_valid"):
            prove_times.append(metrics["prove_time_ms"])
            verify_times.append(metrics["verify_time_ms"])
            print(f"    -> Success: Prove = {metrics['prove_time_ms']:.2f} ms | "
                  f"Verify = {metrics['verify_time_ms']:.2f} ms")
        else:
            print("    -> Simulation failed or proof invalid.")

    # Data Visualization (Metrics for the Master's Thesis)
    if prove_times and verify_times:
        print("\n[*] Generating performance graph...")
        
        clusters = list(range(1, len(prove_times) + 1))
        
        plt.figure(figsize=(10, 6))
        plt.plot(clusters, prove_times, marker='o', linestyle='-', color='b', label='Aggregation Time (Prover)')
        plt.plot(clusters, verify_times, marker='s', linestyle='--', color='g', label='Verification Time')
        
        plt.title('LaBRADOR Signature Aggregation Performance in UAV Networks')
        plt.xlabel('Drone Cluster (Iteration)')
        plt.ylabel('Time (Milliseconds)')
        plt.xticks(clusters)
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.7)
        
        # Save the plot as an image and show it
        plt.savefig("uav_performance_graph.png")
        print("[+] Graph saved as 'uav_performance_graph.png'.")
        # plt.show() # Uncomment this if you are running in an environment with a GUI

if __name__ == "__main__":
    main()