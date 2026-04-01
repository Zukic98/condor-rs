# uav_simulator.py
import subprocess
import json
import matplotlib.pyplot as plt
import os
import time  # NEW: For simulating network delays
import random  # NEW: For adding realistic variance to network delays

NUM_DRONE_CLUSTERS = 5  # Number of simulated clusters/iterations
RUST_BINARY_PATH = "./target/release/examples/simulator_bridge"

# New Network Configuration Parameters (times in milliseconds)
# Assuming a local wireless link (e.g., Wi-Fi / LoRa)
AVG_LATENCY_CM_TO_CH_MS = 10.0   
LATENCY_CM_CH_VARIANCE = 3.0    # Random +/- variance

# Assuming a long-range wireless link (e.g., LTE / Satellite)
AVG_LATENCY_CH_TO_BS_MS = 70.0   
LATENCY_CH_BS_VARIANCE = 15.0   # Random +/- variance

# --- HELPER FUNCTIONS ---

def simulate_network_delay_ms(avg_ms, variance_ms):
    """Simulates realistic wireless network delay by blocking execution."""
    # Add random variance to make it realistic
    actual_delay = avg_ms + random.uniform(-variance_ms, variance_ms)
    # Convert ms to seconds for time.sleep
    time.sleep(actual_delay / 1000.0) 
    return actual_delay

def run_rust_crypto_logic():
    """Executes the Rust binary and parses the mathematical metrics."""
    if not os.path.exists(RUST_BINARY_PATH):
        print(f"Error: Rust binary not found at {RUST_BINARY_PATH}")
        print("Please run: cargo build --release --example simulator_bridge")
        return None

    try:
        # Run binary (Computation only, as implemented in simulator_bridge.rs)
        result = subprocess.run([RUST_BINARY_PATH], capture_output=True, text=True, check=True)
        # Parse the JSON string printed by Rust
        return json.loads(result.stdout.strip())
    except Exception as e:
        print(f"Execution or parsing failed: {e}")
        return None

# --- MAIN SIMULATION ---

def main():
    print(f"--- Starting REALISTIC UAV Network Simulation ({NUM_DRONE_CLUSTERS} iterations) ---")
    
    total_system_times = []
    pure_prove_times = []
    pure_verify_times = []
    
    for i in range(1, NUM_DRONE_CLUSTERS + 1):
        print(f"\n[*] Simulating Drone Cluster {i} workflow:")
        
        # STAGE 1: Network - Cluster Members (CM) send signatures to Cluster Head (CH)
        # We assume they send somewhat in parallel, so we just add base local latency
        print("    [Net] CMs sending signatures to CH...")
        cm_to_ch_latency = simulate_network_delay_ms(AVG_LATENCY_CM_TO_CH_MS, LATENCY_CM_CH_VARIANCE)
        print(f"    -> Local Network Delay: {cm_to_ch_latency:.2f} ms")

        # STAGE 2: Computation - Rust Cluster Head aggregates signatures
        # (This is what your old graph showed)
        metrics = run_rust_crypto_logic()
        
        if not metrics or not metrics.get("is_valid"):
            print("    -> CRITICAL ERROR: Proof generation failed or invalid!")
            break

        pure_prove = metrics["prove_time_ms"]
        pure_verify = metrics["verify_time_ms"]
        
        print(f"    [CPU] CH Aggregation (Prove) Time: {pure_prove:.2f} ms")

        # STAGE 3: Network - Cluster Head sends Proof to Base Station (BS)
        print("    [Net] CH sending aggregate Proof to Base Station...")
        ch_to_bs_latency = simulate_network_delay_ms(AVG_LATENCY_CH_TO_BS_MS, LATENCY_CH_BS_VARIANCE)
        print(f"    -> Long-Range Network Delay: {ch_to_bs_latency:.2f} ms")

        # STAGE 4: Computation - Base Station verifies proof
        print(f"    [CPU] BS Verification Time: {pure_verify:.2f} ms")

        # Calculate Total System Time for this iteration
        total_time = cm_to_ch_latency + pure_prove + ch_to_bs_latency + pure_verify
        print(f"    [+] TOTAL Workflow Time: {total_time:.2f} ms")

        total_system_times.append(total_time)
        pure_prove_times.append(pure_prove)
        pure_verify_times.append(pure_verify)

    # Data Visualization (Updated to show Total vs. Pure Math)
    if total_system_times:
        print("\n[*] Generating updated performance graph...")
        clusters = list(range(1, len(total_system_times) + 1))
        
        plt.figure(figsize=(10, 6))
        
        # Plot Pure Computation Component (Old graph data)
        plt.plot(clusters, pure_prove_times, marker='o', linestyle='-', color='#ADD8E6', label='Aggregation CPU Time')
        
        # Plot TOTAL System Time (including network)
        plt.plot(clusters, total_system_times, marker='^', linestyle='-', color='r', label='Total Workflow Time (Network + CPU)')
        
        # Fill area to emphasize the impact of networking
        plt.fill_between(clusters, pure_prove_times, total_system_times, color='r', alpha=0.1)

        plt.title('Impact of Network Latency on LaBRADOR Performance in UAV Networks')
        plt.xlabel('Drone Cluster (Iteration)')
        plt.ylabel('Time (Milliseconds)')
        plt.xticks(clusters)
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.7)
        
        plt.savefig("realistic_uav_performance.png")
        print("[+] Graph saved as 'realistic_uav_performance.png'.")
        # plt.show() # Uncomment for GUI environments

if __name__ == "__main__":
    main()